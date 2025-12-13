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
    collections::HashSet,
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
}

type Graph = GFA<Vec<u8>, Vec<OptField>>;

trait GraphExt {
    fn longest_path_oriented(&self) -> Vec<(&[u8], bool)>;
    fn reconstruct_path_seq(&self, path: &[(&[u8], bool)]) -> String;
    fn contract_empty_vertices(&mut self);
    fn from_bytes(data: &[u8]) -> Graph;
    fn write(&self, sink: &mut impl std::io::Write) -> std::io::Result<()>;
}

impl GraphExt for Graph {
    fn longest_path_oriented(&self) -> Vec<(&[u8], bool)> {
        let n = self.segments.len();

        // Build direct edges with linear name matching:
        // (from_idx, to_idx, to_orient)
        let mut edges: Vec<(usize, usize, bool)> = Vec::new();

        let name_map = NameMap::build_from_gfa(self);
        for link in &self.links {
            let from_idx = name_map
                .map_name(&link.from_segment)
                .expect("Segment name not found in graph");
            let to_idx = name_map
                .map_name(&link.to_segment)
                .expect("Segment name not found in graph");
            edges.push((from_idx, to_idx, link.to_orient == Orientation::Forward));
        }

        let mut best_path: Vec<(usize, bool)> = Vec::new();

        // DFS-from-every-node using explicit stack (no recursion)
        for start in 0..n {
            // stack frame:
            // (node, orient, next-edge-index, path, visited)
            let mut stack: Vec<(usize, bool, usize, Vec<(usize, bool)>, Vec<bool>)> = Vec::new();

            let mut visited = vec![false; n];
            visited[start] = true;

            stack.push((start, true, 0, vec![(start, true)], visited));

            while let Some((node, orient, mut e_i, path, visited)) = stack.pop() {
                if path.len() > best_path.len() {
                    best_path = path.clone();
                }

                // scan entire edge list to find outgoing edges
                while e_i < edges.len() {
                    let (f, t, t_orient) = edges[e_i];
                    e_i += 1;

                    if f == node && !visited[t] {
                        let mut next_path = path.clone();
                        next_path.push((t, t_orient));

                        let mut next_visit = visited.clone();
                        next_visit[t] = true;

                        // Push continuation frame to resume after this edge
                        stack.push((node, orient, e_i, path, visited));

                        // Push next node
                        stack.push((t, t_orient, 0, next_path, next_visit));
                        break;
                    }
                }
            }
        }

        best_path
            .iter()
            .map(|&(idx, o)| (self.segments[idx].name.as_ref(), o))
            .collect()
    }

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
                seq.push_str(&revcomp_str(s));
            }
        }
        seq
    }

    fn contract_empty_vertices(&mut self) {
        for i in 0..self.segments.len() {
            if !self.segments[i].sequence.is_empty() {
                continue;
            }

            let name = self.segments[i].name.clone();
            // Find incoming and outgoing links
            let incoming: Vec<_> = self
                .links
                .iter()
                .filter(|link| link.to_segment == name)
                .cloned()
                .collect();
            let outgoing: Vec<_> = self
                .links
                .iter()
                .filter(|link| link.from_segment == name)
                .cloned()
                .collect();

            // Create new links bypassing the empty vertex
            for in_link in &incoming {
                for out_link in &outgoing {
                    let new_link = gfa::gfa::Link {
                        from_segment: in_link.from_segment.clone(),
                        from_orient: in_link.from_orient,
                        to_segment: out_link.to_segment.clone(),
                        to_orient: out_link.to_orient,
                        overlap: b"0M".to_vec(),
                        optional: Vec::new(),
                    };
                    self.links.push(new_link);
                }
            }
        }
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
        fasta, k, graph, ..
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
    dfs(&mut tree, root, None);

    match tree.node_weight(root).unwrap() {
        TreeNode::Leaf(_) => panic!("Final node is a leaf"),
        TreeNode::Empty => panic!("Final node is empty"),
        TreeNode::Internal(graph) => {
            // post_process::post_process_graph(graph.clone(), &records);
            graph.write(&mut stdout())?;
        }
    }

    Ok(())
}

fn dfs(tree: &mut Tree, node: NodeIndex, parent: Option<NodeIndex>) {
    let children: Vec<_> = tree
        .neighbors(node)
        .filter(|&n| Some(n) != parent)
        .collect();

    for &child in &children {
        dfs(tree, child, Some(node));
    }

    let new_graph = {
        let label = tree.node_weight(node).unwrap();
        match label {
            TreeNode::Leaf(_) | TreeNode::Internal(_) => None,
            TreeNode::Empty => {
                if children.len() == 2 {
                    let child1 = tree.node_weight(children[0]).unwrap().clone();
                    let child2 = tree.node_weight(children[1]).unwrap().clone();
                    Some(align(&child1, &child2))
                } else if children.len() == 3 {
                    let child1 = tree.node_weight(children[0]).unwrap().clone();
                    let child2 = tree.node_weight(children[1]).unwrap().clone();
                    let child3 = tree.node_weight(children[2]).unwrap().clone();

                    let new_graph = align(&child1, &child2);
                    Some(align(&TreeNode::Internal(new_graph), &child3))
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

fn align(a: &TreeNode, b: &TreeNode) -> Graph {
    match (a, b) {
        (TreeNode::Empty, _) | (_, TreeNode::Empty) => unreachable!(),
        (TreeNode::Leaf(a), TreeNode::Internal(b)) | (TreeNode::Internal(b), TreeNode::Leaf(a)) => {
            eprintln!("Sequence Graph alignment");
            minigraph::align_seq_graph(a, b)
        }
        (TreeNode::Internal(a), TreeNode::Internal(b)) => {
            eprintln!("Graph Graph alignment");
            minigraph::align_graph(a, b)
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
        let rc = revcomp_bytes(kmer);
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
fn revcomp_bytes(seq: &[u8]) -> Vec<u8> {
    seq.iter()
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

fn revcomp_str(seq: &str) -> String {
    String::from_utf8(revcomp_bytes(seq.as_bytes())).unwrap()
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
