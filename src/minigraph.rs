use gax::gaf::GafRecord;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::{io::Write, process::Command};
use tempfile::NamedTempFile;

static QUERY_COUNTER: AtomicUsize = AtomicUsize::new(0);

pub fn align_graph(a: &crate::Graph, b: &crate::Graph) -> crate::Graph {
    // choose a representative sequence from graph a
    let repr_path = a.longest_path();
    let repr = repr_path
        .iter()
        .map(|&i| &a.vertices[i])
        .cloned()
        .collect::<Vec<_>>()
        .join("");

    // call minigraph and parse the output
    let mut graph: crate::Graph = align_seq_graph(&repr, b);

    // align repr on new graph
    let mut repr_file = NamedTempFile::new().unwrap();

    let id = QUERY_COUNTER.fetch_add(1, Ordering::SeqCst);
    write!(repr_file, ">query{}\n{}", id, repr).unwrap();

    let mut graph_file = NamedTempFile::new().unwrap();
    graph.write_to_gfa(&mut graph_file).unwrap();

    repr_file.flush().unwrap();
    graph_file.flush().unwrap();

    let out_file = NamedTempFile::new().unwrap();
    Command::new("GraphAligner")
        .args(["--seeds-minimizer-ignore-frequent", "0.0"])
        .args(["--seeds-minimizer-density", "-1"])
        .args(["--seeds-extend-density", "-1"])
        .args(["--seeds-mxm-length", "15"])
        .args(["--bandwidth", "200"])
        .args(["--tangle-effort", "-1"])
        .args(["--precise-clipping", "0.9"])
        .arg("-g")
        .arg(graph_file.path())
        .arg("-f")
        .arg(repr_file.path())
        .arg("-a")
        .arg(out_file.path())
        .args(["-x", "vg"])
        .output()
        .expect("failed to execute GraphAligner");

    let gaf_file = File::open(out_file.path()).unwrap();
    let gaf = gax::gaf::parse(gaf_file).unwrap();

    // compute reverse adj_list
    let n = a.vertices.len();
    let mut a_reverse_adj_list = vec![Vec::new(); n];
    for (u, neighbors) in a.adj_list.iter().enumerate() {
        for &v in neighbors {
            a_reverse_adj_list[v].push(u);
        }
    }

    let mut mapping = build_mapping(a, &graph, &repr_path, &gaf);
    let mut old_to_new: HashMap<usize, usize> = HashMap::new();
    let mut added_back: HashSet<usize> = HashSet::new();
    let mut added_forw: HashSet<usize> = HashSet::new();

    for (&node, mapping_i) in repr_path.iter().zip(0..mapping.len()) {
        let Some(MapEntry {
            start_node,
            start_offset,
            end_node,
            end_offset,
        }) = mapping[mapping_i]
        else {
            continue;
        };

        if !(start_node < graph.vertices.len()
            && end_node < graph.vertices.len()
            && start_offset <= graph.vertices[start_node].len()
            && end_offset <= graph.vertices[end_node].len())
        {
            println!("Skipping node {} due to invalid mapping", node);
            println!("{:#?}", mapping[mapping_i]);
            println!("vertices len {}", graph.vertices.len());
            println!("start node len {}", graph.vertices[start_node].len());
            println!("end node len {}", graph.vertices[end_node].len());
            panic!();
        }

        let new_root = if start_offset > 0 {
            // split the start_node in the new graph at start_offset
            split_node(&mut graph, start_node, start_offset, &mut mapping)
        } else {
            start_node
        };
        old_to_new.insert(node, new_root);

        // explore backwards in the old graph and add edges
        let mut stack = vec![node];
        while let Some(current) = stack.pop() {
            for &prev in &a_reverse_adj_list[current] {
                if repr_path.contains(&prev) || added_back.contains(&prev) {
                    continue;
                }
                added_back.insert(prev);
                let new_node = if let Some(&new_node) = old_to_new.get(&prev) {
                    new_node
                } else {
                    // node needs to be added
                    graph.vertices.push(a.vertices[prev].clone());
                    graph.adj_list.push(vec![]);
                    old_to_new.insert(prev, graph.vertices.len() - 1);
                    graph.vertices.len() - 1
                };
                graph.adj_list[new_node].push(old_to_new[&current]);
                stack.push(prev);
            }
        }

        let new_root = if end_offset < graph.vertices[end_node].len() {
            // split the end_node in the new graph at end_offset
            split_node(&mut graph, end_node, end_offset, &mut mapping)
        } else {
            end_node
        };
        old_to_new.insert(node, new_root);

        // explore forwards
        let mut stack = vec![node];
        while let Some(current) = stack.pop() {
            for &next in &a.adj_list[current] {
                if repr_path.contains(&next) || added_forw.contains(&next) {
                    continue;
                }
                added_forw.insert(next);
                let new_node = if let Some(&new_node) = old_to_new.get(&next) {
                    new_node
                } else {
                    // node needs to be added
                    graph.vertices.push(a.vertices[next].clone());
                    graph.adj_list.push(vec![]);
                    old_to_new.insert(next, graph.vertices.len() - 1);
                    graph.vertices.len() - 1
                };
                graph.adj_list[old_to_new[&current]].push(new_node);
                stack.push(next);
            }
        }
    }

    // remove empty nodes
    graph.contract_empty_vertices();

    graph
}

#[derive(Clone, Copy, Debug)]
struct MapEntry {
    start_node: usize,
    start_offset: usize,
    end_node: usize,
    end_offset: usize,
}

fn build_mapping(
    old_graph: &crate::Graph,
    new_graph: &crate::Graph,
    repr_path: &[usize],
    gaf: &[GafRecord],
) -> Vec<Option<MapEntry>> {
    let mut mapping: Vec<Option<MapEntry>> = vec![None; repr_path.len()];

    // compute cumulative start positions of old nodes
    let mut i = 0;
    let repr_starts: Vec<usize> = repr_path
        .iter()
        .map(|&idx| {
            let start = i;
            i += old_graph.vertices[idx].len();
            start
        })
        .collect();

    for ((&node, &query_pos), entry) in repr_path
        .iter()
        .zip(repr_starts.iter())
        .zip(mapping.iter_mut())
    {
        let mut start = None;
        for record in gaf {
            start = map_query_pos_to_node(record, query_pos, new_graph);
            if start.is_some() {
                break;
            }
        }

        let mut end = None;
        for record in gaf {
            end = map_query_pos_to_node(
                record,
                query_pos + old_graph.vertices[node].len(),
                new_graph,
            );
            if end.is_some() {
                break;
            }
        }

        if let (Some((start_node, start_offset)), Some((end_node, end_offset))) = (start, end) {
            *entry = Some(MapEntry {
                start_node,
                start_offset,
                end_node,
                end_offset,
            });
        };
    }

    mapping
}

fn map_query_pos_to_node(
    record: &GafRecord,
    target_query_pos: usize,
    graph: &crate::Graph,
) -> Option<(usize, usize)> {
    if target_query_pos < record.query_start as usize
        || target_query_pos > record.query_end as usize
    {
        return None;
    }

    let cigar_ops = parse_cigar(&record.opt_fields["cg"].1);

    let mut current_query_pos = record.query_start as usize;
    let mut path_idx = 0;
    let mut node_offset = 0;

    for op in &cigar_ops {
        let mut node_idx = record.path[path_idx].name.parse::<usize>().unwrap() - 1;
        let mut node_len = graph.vertices[node_idx].len();

        match *op {
            CigarOp::Match(len) | CigarOp::Mismatch(len) => {
                let query_start = current_query_pos;
                let query_end = current_query_pos + len;

                if target_query_pos < query_end {
                    let offset_into_op = target_query_pos - query_start;
                    let mut ref_offset = node_offset + offset_into_op;

                    // walk forward across nodes if needed
                    while ref_offset >= node_len {
                        ref_offset -= node_len;
                        path_idx += 1;
                        if path_idx >= record.path.len() {
                            return None;
                        }
                        node_idx = record.path[path_idx].name.parse::<usize>().unwrap() - 1;
                        node_len = graph.vertices[node_idx].len();
                    }

                    return Some((node_idx, ref_offset));
                }

                current_query_pos = query_end;
                node_offset += len;

                while node_offset >= node_len {
                    node_offset -= node_len;
                    path_idx += 1;
                    if path_idx >= record.path.len() {
                        return None;
                    }
                    node_idx = record.path[path_idx].name.parse::<usize>().unwrap() - 1;
                    node_len = graph.vertices[node_idx].len();
                }
            }
            CigarOp::Insertion(len) => {
                current_query_pos += len;
                if current_query_pos > target_query_pos {
                    return None;
                }
            }
            CigarOp::Deletion(len) => {
                node_offset += len;
                while node_offset >= node_len {
                    node_offset -= node_len;
                    path_idx += 1;
                    if path_idx >= record.path.len() {
                        return None;
                    }
                    node_idx = record.path[path_idx].name.parse::<usize>().unwrap() - 1;
                    node_len = graph.vertices[node_idx].len();
                }
            }
        }
    }

    None
}

pub fn align_seq_graph(a: &str, b: &crate::Graph) -> crate::Graph {
    let mut query = NamedTempFile::new().unwrap();

    let id = QUERY_COUNTER.fetch_add(1, Ordering::SeqCst);
    write!(query, ">query{}\n{}", id, a).unwrap();

    let mut graph = NamedTempFile::new().unwrap();
    b.write_to_gfa(&mut graph).unwrap();

    query.flush().unwrap();
    graph.flush().unwrap();

    let out = Command::new("minigraph")
        .args(["-cxggs"])
        .arg(graph.path())
        .arg(query.path())
        .output()
        .expect("failed to execute minigraph");
    crate::Graph::from_gfa(&out.stdout[..])
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum CigarOp {
    Match(usize),
    Mismatch(usize),
    Insertion(usize),
    Deletion(usize),
}

fn parse_cigar(cigar: &str) -> Vec<CigarOp> {
    let mut ops = Vec::new();
    let mut num = 0;

    for c in cigar.chars() {
        if c.is_ascii_digit() {
            num = num * 10 + c.to_digit(10).unwrap() as usize;
        } else {
            let op = match c {
                '=' => CigarOp::Match(num),
                'X' => CigarOp::Mismatch(num),
                'I' => CigarOp::Insertion(num),
                'D' => CigarOp::Deletion(num),
                _ => panic!("Unexpected CIGAR op: {}", c),
            };
            ops.push(op);
            num = 0;
        }
    }

    ops
}

pub fn align_seq_seq(a: &str, b: &str) -> crate::Graph {
    let gfa = crate::Graph {
        vertices: vec![a.to_string()],
        adj_list: vec![vec![]],
    };
    align_seq_graph(b, &gfa)
}

/// Split a node in [0, split_point) and [split_point, node.len()).
/// Returns indices new node index if split happened.
fn split_node(
    graph: &mut crate::Graph,
    node_index: usize,
    split_point: usize,
    mapping: &mut [Option<MapEntry>],
) -> usize {
    let original = graph.vertices[node_index].clone();
    graph.vertices[node_index] = original[..split_point].to_string();
    graph.vertices.push(original[split_point..].to_string());
    let new_idx = graph.vertices.len() - 1;

    let old_adj = graph.adj_list[node_index].clone();
    graph.adj_list.push(old_adj);
    graph.adj_list[node_index] = vec![new_idx];

    for entry in mapping.iter_mut().filter_map(|e| e.as_mut()) {
        if entry.start_node == node_index && entry.start_offset >= split_point {
            entry.start_node = new_idx;
            entry.start_offset -= split_point;
        }
        if entry.end_node == node_index && entry.end_offset >= split_point {
            entry.end_node = new_idx;
            entry.end_offset -= split_point;
        }
    }

    new_idx
}
