use bio::io::fasta;
use bstr::io::BufReadExt;
use clap::Parser;
use gfa::{
    gfa::{name_conversion::NameMap, Orientation, GFA},
    optfields::OptField,
    parser::GFAParser,
    writer::write_gfa,
};
use petgraph::{graph::NodeIndex, Undirected};
use rayon::prelude::*;
use speedytree::{DistanceMatrix, NeighborJoiningSolver, RapidBtrees};
use std::{
    collections::{HashMap, HashSet, VecDeque},
    error::Error,
    io::{stdout, BufReader},
    ops::{Deref, DerefMut},
};

mod minigraph;
mod post_process;

// #[global_allocator]
// static ALLOCATOR: ConstLimit<System, 5_000_000_000> = ConstLimit::new(System);

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    #[arg(long, short)]
    fasta: String,
    #[arg(short)]
    k: usize,
    #[arg(long, short)]
    graph: Option<String>,
    #[arg(long)]
    minigraph_only: bool,
}

type Graph = GFA<Vec<u8>, Vec<OptField>>;

trait GraphExt {
    fn longest_path_oriented(&self) -> Vec<(&[u8], bool)>;
    fn reconstruct_path_seq(&self, path: &[(&[u8], bool)]) -> String;
    fn from_bytes(data: &[u8]) -> Graph;
    fn write(&self, sink: &mut impl std::io::Write) -> std::io::Result<()>;
    fn mass_rename(&mut self);
}

impl GraphExt for Graph {
    fn longest_path_oriented(&self) -> Vec<(&[u8], bool)> {
        let n = self.segments.len();
        let size = n * 2;

        let mut adj = vec![Vec::new(); size];
        let name_map = NameMap::build_from_gfa(self);

        for link in &self.links {
            let from = name_map
                .map_name(&link.from_segment)
                .expect("Segment name not found");
            let to = name_map
                .map_name(&link.to_segment)
                .expect("Segment name not found");

            if !(from < n) {
                eprintln!(
                    "from index {} out of range (n = {}) for segment {:?}",
                    from,
                    n,
                    std::str::from_utf8(&link.from_segment),
                );
                continue;
            }
            if !(to < n) {
                eprintln!(
                    "to index {} out of range (n = {}) for segment {:?}",
                    to,
                    n,
                    std::str::from_utf8(&link.from_segment),
                );
                continue;
            }

            // skip backward links for now
            if (link.from_orient == Orientation::Backward)
                || (link.to_orient == Orientation::Backward)
            {
                continue;
            }

            let from = from + n * (link.from_orient == Orientation::Forward) as usize;
            let to = to + n * (link.to_orient == Orientation::Forward) as usize;

            if from >= size || to >= size {
                eprintln!(
                    "Skipping invalid link from {} to {} (size = {})",
                    from, to, size
                );
                continue;
            }

            adj[from].push(to);
        }

        let mut indeg = vec![0usize; size];
        for u in 0..size {
            for &v in &adj[u] {
                indeg[v] += 1;
            }
        }

        let mut queue = VecDeque::new();
        let mut indeg_kahn = indeg.clone();

        for u in 0..size {
            if indeg_kahn[u] == 0 {
                queue.push_back(u);
            }
        }

        let mut topo = Vec::with_capacity(size);
        while let Some(u) = queue.pop_front() {
            topo.push(u);
            for &v in &adj[u] {
                indeg_kahn[v] -= 1;
                if indeg_kahn[v] == 0 {
                    queue.push_back(v);
                }
            }
        }

        let path = {
            let mut dist = vec![1usize; size];
            let mut parent = vec![None; size];

            for &u in &topo {
                for &v in &adj[u] {
                    if dist[u] + 1 > dist[v] {
                        dist[v] = dist[u] + 1;
                        parent[v] = Some(u);
                    }
                }
            }

            let end = (0..size).max_by_key(|&i| dist[i]).unwrap();

            let mut path = Vec::new();
            let mut cur = Some(end);
            while let Some(u) = cur {
                path.push(u);
                cur = parent[u];
            }
            path.reverse();
            path
        };

        path.into_iter()
            .map(|idx| {
                let segment_idx = idx % n;
                let forward = idx >= n;
                (&self.segments[segment_idx].name[..], forward)
            })
            .collect()
    }
    /*
         // safe and slow version
        fn longest_path_oriented(&self) -> Vec<(&[u8], bool)> {
            let n = self.segments.len();
            let size = n * 2;

            let mut adj = vec![Vec::new(); size];
            let name_map = NameMap::build_from_gfa(self);

            for link in &self.links {
                let from = name_map
                    .map_name(&link.from_segment)
                    .expect("Segment name not found");
                let to = name_map
                    .map_name(&link.to_segment)
                    .expect("Segment name not found");

                assert!(
                    from < n,
                    "from index {} out of range (n = {}) for segment {:?}",
                    from,
                    n,
                    std::str::from_utf8(&link.from_segment)
                );
                assert!(
                    to < n,
                    "to index {} out of range (n = {}) for segment {:?}",
                    to,
                    n,
                    std::str::from_utf8(&link.to_segment)
                );

                let from = from + n * (link.from_orient == Orientation::Forward) as usize;
                let to = to + n * (link.to_orient == Orientation::Forward) as usize;

                adj[from].push(to);
            }

            let mut indeg = vec![0usize; size];
            for u in 0..size {
                for &v in &adj[u] {
                    indeg[v] += 1;
                }
            }

            let mut queue = VecDeque::new();
            let mut indeg_kahn = indeg.clone();

            for u in 0..size {
                if indeg_kahn[u] == 0 {
                    queue.push_back(u);
                }
            }

            let mut topo = Vec::with_capacity(size);
            while let Some(u) = queue.pop_front() {
                topo.push(u);
                for &v in &adj[u] {
                    indeg_kahn[v] -= 1;
                    if indeg_kahn[v] == 0 {
                        queue.push_back(v);
                    }
                }
            }

            let path = if topo.len() == size {
                let mut dist = vec![1usize; size];
                let mut parent = vec![None; size];

                for &u in &topo {
                    for &v in &adj[u] {
                        if dist[u] + 1 > dist[v] {
                            dist[v] = dist[u] + 1;
                            parent[v] = Some(u);
                        }
                    }
                }

                let end = (0..size).max_by_key(|&i| dist[i]).unwrap();

                let mut path = Vec::new();
                let mut cur = Some(end);
                while let Some(u) = cur {
                    path.push(u);
                    cur = parent[u];
                }
                path.reverse();
                path
            } else {
                eprintln!("Graph has cycles, using DFS for longest path");

                fn dfs(
                    u: usize,
                    adj: &Vec<Vec<usize>>,
                    visited: &mut Vec<bool>,
                    path: &mut Vec<usize>,
                    best: &mut Vec<usize>,
                ) {
                    visited[u] = true;
                    path.push(u);

                    if path.len() > best.len() {
                        *best = path.clone();
                    }

                    for &v in &adj[u] {
                        if !visited[v] {
                            dfs(v, adj, visited, path, best);
                        }
                    }

                    path.pop();
                    visited[u] = false;
                }

                let mut visited = vec![false; size];
                let mut path = Vec::new();
                let mut best = Vec::new();

                for u in 0..size {
                    dfs(u, &adj, &mut visited, &mut path, &mut best);
                }
                best
            };

            path.into_iter()
                .map(|idx| {
                    let segment_idx = idx % n;
                    let forward = idx >= n;
                    (&self.segments[segment_idx].name[..], forward)
                })
                .collect()
        }
    */

    /// Convert a path of (index, orientation) to a full sequence
    fn reconstruct_path_seq(&self, path: &[(&[u8], bool)]) -> String {
        let name_map = NameMap::build_from_gfa(self);
        let mut seq = String::new();
        for &(name, forward) in path {
            let idx = name_map
                .map_name(name)
                .expect("Segment name not found in graph");
            let s = std::str::from_utf8(&self.segments[idx].sequence).unwrap();
            if forward {
                seq.push_str(s);
            } else {
                seq.push_str(&String::from_utf8(revcomp(s)).unwrap());
            }
        }
        seq
    }

    fn from_bytes(data: &[u8]) -> Self {
        let parser = GFAParser::new();
        let tolerance = Default::default();
        let mut gfa: GFA<Vec<u8>, Vec<OptField>> = GFA::new();
        let lines = BufReader::new(data).byte_lines();
        for line in lines {
            match parser.parse_gfa_line(line.unwrap().as_ref()) {
                Ok(parsed) => gfa.insert_line(parsed),
                Err(err) if err.can_safely_continue(&tolerance) => (),
                Err(err) => panic!("Error parsing GFA line: {}", err),
            };
        }
        gfa
    }

    fn write(&self, sink: &mut impl std::io::Write) -> std::io::Result<()> {
        let mut b = String::new();
        write_gfa(self, &mut b);
        sink.write_all(b.as_bytes())
    }

    fn mass_rename(&mut self) {
        let mut id = 0;
        let mut name_map = HashMap::new();
        for segment in self.segments.iter_mut() {
            let old_name = segment.name.clone();
            let new_name = format!("s{}", id);
            name_map.insert(old_name, new_name.clone());
            segment.name = new_name.into_bytes();
            id += 1;
        }
        let mut i = 0;
        while i < self.links.len() {
            let link = &mut self.links[i];
            match (
                name_map.get(&link.from_segment),
                name_map.get(&link.to_segment),
            ) {
                (Some(from_segment), Some(to_segment)) => {
                    link.from_segment = from_segment.clone().into_bytes();
                    link.to_segment = to_segment.clone().into_bytes();
                }
                _ => {
                    // link points to a non existing segment
                    self.links.swap_remove(i);
                    continue;
                }
            }
            i += 1;
        }
    }
}

#[derive(Debug, Clone, Default)]
enum TreeNode {
    Leaf(String),
    Internal(Graph),
    #[default]
    Empty,
}

type Tree = petgraph::Graph<TreeNode, f64, Undirected>;

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();
    let Args {
        fasta, k, graph, minigraph_only, ..
    } = args;

    let data = fasta::Reader::from_file(&fasta)?;
    let records: Vec<fasta::Record> = data.records().collect::<Result<_, _>>()?;

    if let Some(graph_path) = graph {
        let data = std::fs::read(&graph_path)?;
        let graph = Graph::from_bytes(&data);
        let new_graph = post_process::post_process_graph(graph, &records);
        new_graph.write(&mut stdout())?;
        return Ok(());
    }

    let distance_matrix = {
        let kmers_set: Vec<_> = (&records)
            .into_par_iter()
            .enumerate()
            .map(|(i, record)| {
                eprintln!("Processing record {}", i);
                extract_kmers(record.seq(), k)
            })
            .collect();

        let labels = records
            .iter()
            .map(|r| r.id().to_string())
            .collect::<Vec<_>>();

        let matrix = compute_distance_matrix(&kmers_set);
        DistanceMatrix::build(matrix, labels.clone())?
    };

    eprintln!("Building tree");
    let nj_tree = NeighborJoiningSolver::<RapidBtrees>::default(distance_matrix).solve()?;

    let mut tree = Tree::new_undirected();
    for node in nj_tree.raw_nodes() {
        let new_label = if node.weight.is_empty() {
            TreeNode::Empty
        } else {
            let index = records.iter().position(|r| r.id() == node.weight).unwrap();
            TreeNode::Leaf(String::from_utf8(records[index].seq().to_vec()).unwrap())
        };
        tree.add_node(new_label);
    }
    for edge in nj_tree.raw_edges() {
        tree.add_edge(edge.source(), edge.target(), edge.weight);
    }

    let root = tree
        .node_indices()
        .find(|node| tree.neighbors(*node).count() == 3)
        .unwrap();

    eprintln!("Aligning");
    dfs(&mut tree, root, None, minigraph_only);

    let Some(root_node) = tree.remove_node(root) else {
        panic!("Root node not found");
    };
    match root_node {
        TreeNode::Leaf(_) => panic!("Final node is a leaf"),
        TreeNode::Empty => panic!("Final node is empty"),
        TreeNode::Internal(graph) => {
            let g = post_process::post_process_graph(graph, &records);
            g.write(&mut stdout())?;
        }
    }

    Ok(())
}

fn dfs(tree: &mut Tree, node: NodeIndex, parent: Option<NodeIndex>, minigraph_only: bool) {
    let children: Vec<_> = tree
        .neighbors(node)
        .filter(|&n| Some(n) != parent)
        .collect();

    for &child in &children {
        dfs(tree, child, Some(node), minigraph_only);
    }

    let new_graph = {
        let label = tree.node_weight(node).unwrap();
        match label {
            TreeNode::Leaf(_) | TreeNode::Internal(_) => None,
            TreeNode::Empty => {
                if children.len() == 2 {
                    let child1 = tree.node_weight(children[0]).unwrap().clone();
                    let child2 = tree.node_weight(children[1]).unwrap().clone();
                    Some(align(&child1, &child2, minigraph_only))
                } else if children.len() == 3 {
                    let child1 = tree.node_weight(children[0]).unwrap().clone();
                    let child2 = tree.node_weight(children[1]).unwrap().clone();
                    let child3 = tree.node_weight(children[2]).unwrap().clone();

                    let new_graph = align(&child1, &child2, minigraph_only);
                    Some(align(
                        &TreeNode::Internal(new_graph),
                        &child3,
                        minigraph_only,
                    ))
                } else {
                    unreachable!()
                }
            }
        }
    };
    if let Some(new_graph) = new_graph {
        *tree.node_weight_mut(node).unwrap() = TreeNode::Internal(new_graph);
    }
}

fn align(a: &TreeNode, b: &TreeNode, minigraph_only: bool) -> Graph {
    match (a, b) {
        (TreeNode::Empty, _) | (_, TreeNode::Empty) => unreachable!(),
        (TreeNode::Leaf(a), TreeNode::Internal(b)) | (TreeNode::Internal(b), TreeNode::Leaf(a)) => {
            eprintln!("Sequence Graph alignment");
            minigraph::align_seq_graph(a, b)
        }
        (TreeNode::Internal(a), TreeNode::Internal(b)) => {
            eprintln!("Graph Graph alignment");
            minigraph::align_graph(a, b, minigraph_only)
        }
        (TreeNode::Leaf(a), TreeNode::Leaf(b)) => {
            eprintln!("Sequence Sequence alignment");
            minigraph::align_seq_seq(a, b)
        }
    }
}

fn compute_distance_matrix(kmers_set: &[HashSet<Vec<u8>>]) -> Vec<Vec<f64>> {
    let mut distance_matrix = vec![vec![1.0; kmers_set.len()]; kmers_set.len()];
    let indexes = (0..kmers_set.len())
        .flat_map(|i| {
            (0..kmers_set.len()).filter_map(move |j| if i >= j { None } else { Some((i, j)) })
        })
        .collect::<Vec<_>>();

    let distances = indexes
        .into_par_iter()
        .map(|(i, j)| {
            let a = &kmers_set[i];
            let b = &kmers_set[j];
            let intersection = a.intersection(b).count();
            let union = a.len() + b.len() - intersection;
            let jaccard = intersection as f64 / union as f64;
            (i, j, jaccard)
        })
        .collect::<Vec<_>>();

    for (i, j, jaccard) in distances {
        distance_matrix[i][j] = 1. - jaccard;
        distance_matrix[j][i] = 1. - jaccard;
    }
    distance_matrix
}

fn extract_kmers(seq: &[u8], k: usize) -> HashSet<Vec<u8>> {
    let mut kmers = HashSet::new();
    for i in 0..seq.len() - k {
        // make kmer canonical
        let kmer = &seq[i..i + k];
        let rc = revcomp(kmer);
        // pick lexicographically smaller
        let canonical = if kmer <= rc.as_slice() {
            kmer.to_vec()
        } else {
            rc
        };

        kmers.insert(canonical);
    }
    kmers
}
fn revcomp(seq: impl AsRef<[u8]>) -> Vec<u8> {
    seq.as_ref()
        .iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            other => other,
        })
        .collect()
}

#[derive(Debug, Clone)]
pub struct AdjList(Vec<Vec<usize>>);

impl DerefMut for AdjList {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}
impl Deref for AdjList {
    type Target = Vec<Vec<usize>>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl AdjList {
    pub fn new(from: usize) -> Self {
        Self(vec![Vec::new(); from])
    }

    #[inline]
    pub fn is_edge(&self, i: usize, j: usize) -> bool {
        self[i].contains(&j)
    }

    #[inline]
    pub fn add_edge(&mut self, i: usize, j: usize) {
        if self[i].contains(&j) {
            return;
        }
        self[i].push(j);
    }

    pub fn remove_edge(&mut self, i: usize, j: usize) {
        self[i].retain(|&x| x != j);
    }
}
