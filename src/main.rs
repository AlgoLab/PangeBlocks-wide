use bio::io::fasta;
use bstr::io::BufReadExt;
use clap::Parser;
use gfa::{gfa::GFA, optfields::OptField, parser::GFAParser};
use limit_alloc::ConstLimit;
use petgraph::{graph::NodeIndex, Undirected};
use rayon::prelude::*;
use speedytree::{DistanceMatrix, NeighborJoiningSolver, RapidBtrees};
use std::{
    alloc::System,
    collections::{HashMap, HashSet},
    error::Error,
    io::{self, stdout, BufReader},
    ops::{Deref, DerefMut},
};

mod minigraph;

#[global_allocator]
static ALLOCATOR: ConstLimit<System, 5_000_000_000> = ConstLimit::new(System);

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    #[arg(long, short)]
    fasta: String,
    #[arg(short)]
    k: usize,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Graph {
    vertices: Vec<String>,
    adj_list: Vec<Vec<usize>>,
}

impl Graph {
    /// Finds the longest simple (non-cyclic) path
    fn longest_path(&self) -> Vec<usize> {
        let n = self.vertices.len();
        let mut stack = vec![false; n];
        let mut memo: Vec<Option<Vec<usize>>> = vec![None; n];

        let mut best_path = Vec::new();

        for v in 0..n {
            let path = Self::dfs(v, &self.adj_list, &mut stack, &mut memo);
            if path.len() > best_path.len() {
                best_path = path;
            }
        }

        best_path
    }

    fn dfs(
        node: usize,
        adj: &Vec<Vec<usize>>,
        stack: &mut Vec<bool>,
        memo: &mut Vec<Option<Vec<usize>>>,
    ) -> Vec<usize> {
        if let Some(cached) = &memo[node] {
            return cached.clone();
        }

        if stack[node] {
            return vec![node];
        }

        stack[node] = true;

        let mut best = vec![node];
        for &next in &adj[node] {
            if !stack[next] {
                let mut candidate = Self::dfs(next, adj, stack, memo);
                if candidate.len() + 1 > best.len() {
                    let mut new_path = vec![node];
                    new_path.append(&mut candidate);
                    best = new_path;
                }
            }
        }

        stack[node] = false;
        memo[node] = Some(best.clone());
        best
    }

    pub fn contract_empty_vertices(&mut self) {
        // compute reverse adj_list
        let mut reverse_adj_list = vec![Vec::new(); self.vertices.len()];
        for (u, neighbors) in self.adj_list.iter().enumerate() {
            for &v in neighbors {
                reverse_adj_list[v].push(u);
            }
        }

        // collect empty vertices
        let empty: HashSet<usize> = self
            .vertices
            .iter()
            .enumerate()
            .filter(|(_, v)| v.is_empty())
            .map(|(i, _)| i)
            .collect();

        // create mapping old->new index for non-empty vertices
        let mut old_to_new = HashMap::new();
        let mut new_vertices = Vec::new();
        let mut next_idx = 0;

        for (i, v) in self.vertices.iter().enumerate() {
            if !empty.contains(&i) {
                old_to_new.insert(i, next_idx);
                new_vertices.push(v.clone());
                next_idx += 1;
            }
        }

        // build new adjacency structures
        let mut new_adj: Vec<HashSet<usize>> = vec![HashSet::new(); new_vertices.len()];
        let mut new_rev: Vec<HashSet<usize>> = vec![HashSet::new(); new_vertices.len()];

        for i in 0..self.vertices.len() {
            if empty.contains(&i) {
                // Stitch predecessors to successors
                for &pred in &reverse_adj_list[i] {
                    if empty.contains(&pred) {
                        continue;
                    }
                    for &succ in &self.adj_list[i] {
                        if empty.contains(&succ) {
                            continue;
                        }
                        let np = old_to_new[&pred];
                        let ns = old_to_new[&succ];
                        if np != ns {
                            new_adj[np].insert(ns);
                            new_rev[ns].insert(np);
                        }
                    }
                }
            } else {
                // Copy normal edges that don't go through empty nodes
                let ni = old_to_new[&i];
                for &succ in &self.adj_list[i] {
                    if !empty.contains(&succ) {
                        let ns = old_to_new[&succ];
                        new_adj[ni].insert(ns);
                        new_rev[ns].insert(ni);
                    }
                }
            }
        }

        self.vertices = new_vertices;
        self.adj_list = new_adj
            .into_iter()
            .map(|s| s.into_iter().collect())
            .collect();
    }

    pub fn from_gfa(data: &[u8]) -> Graph {
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

        let mut graph = Graph::default();
        let mut id_to_idx: HashMap<Vec<u8>, usize> = HashMap::new();

        // Segments → vertices
        for (i, seg) in gfa.segments.iter().enumerate() {
            let sequence = String::from_utf8(seg.sequence.clone()).unwrap();
            graph.vertices.push(sequence);
            id_to_idx.insert(seg.name.clone(), i);
            graph.adj_list.push(Vec::new());
        }

        // Links → adjacency list
        for link in gfa.links.iter() {
            let from = *id_to_idx.get(&link.from_segment).unwrap();
            let to = *id_to_idx.get(&link.to_segment).unwrap();
            graph.adj_list[from].push(to);
        }
        graph
    }

    fn write_to_gfa(&self, sink: &mut impl io::Write) -> io::Result<()> {
        writeln!(sink, "H\tVN:Z:1.0")?;

        for (i, v) in self.vertices.iter().enumerate() {
            writeln!(sink, "S\t{}\t{}", i + 1, v)?;
        }

        for (i, neighbors) in self.adj_list.iter().enumerate() {
            for &j in neighbors {
                writeln!(sink, "L\t{}\t+\t{}\t+\t0M", i + 1, j + 1)?;
            }
        }

        Ok(())
    }
}

#[derive(Debug, Clone, Default)]
enum TreeNode {
    Leaf(String),
    Internal(Graph),
    #[default]
    Empty,
}

impl std::fmt::Display for TreeNode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            TreeNode::Leaf(s) => write!(f, "Leaf({})", s.len()),
            TreeNode::Internal(g) => {
                write!(f, "Graph({})", g.vertices.len())
            }
            TreeNode::Empty => write!(f, "Empty"),
        }
    }
}

type Tree = petgraph::Graph<TreeNode, f64, Undirected>;

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();
    let Args { fasta, k, .. } = args;

    let data = fasta::Reader::from_file(&fasta)?;
    let records: Vec<fasta::Record> = data.records().take(50).collect::<Result<_, _>>()?;

    let distance_matrix = {
        let kmers_set: Vec<_> = (&records)
            .into_par_iter()
            .enumerate()
            .map(|(i, record)| {
                println!("Processing record {}", i);
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

    println!("Building tree");
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

    println!("Aligning");
    dfs(&mut tree, root, None);

    match tree.node_weight(root).unwrap() {
        TreeNode::Leaf(_) => panic!("Final node is a leaf"),
        TreeNode::Empty => panic!("Final node is empty"),
        TreeNode::Internal(graph) => {
            graph.write_to_gfa(&mut stdout()).unwrap();
        }
    }

    Ok(())
}

fn dfs(tree: &mut Tree, node: NodeIndex, parent: Option<NodeIndex>) {
    let children: Vec<_> = tree
        .neighbors(node)
        .filter(|&n| Some(n) != parent)
        .collect();

    // recursively process children in parallel
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
            println!("Sequence Graph alignment");
            minigraph::align_seq_graph(a, b)
        }
        (TreeNode::Internal(a), TreeNode::Internal(b)) => {
            println!("Graph Graph alignment");
            minigraph::align_graph(a, b)
        }
        (TreeNode::Leaf(a), TreeNode::Leaf(b)) => {
            println!("Sequence Sequence alignment");
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
fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b'N',
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
