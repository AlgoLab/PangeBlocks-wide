use gax::gaf::GafRecord;
use gfa::cigar::{CIGAROp, CIGAR};
use gfa::gfa::name_conversion::NameMap;
use gfa::gfa::{Header, Link, Segment};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::{io::Write, process::Command};
use tempfile::NamedTempFile;

use crate::{Graph, GraphExt};

pub fn align_graph(a: &Graph, b: &Graph) -> Graph {
    // choose a representative sequence from graph a
    let repr_path = a.longest_path_oriented();
    let repr = a.reconstruct_path_seq(&repr_path);

    eprintln!(
        "Starting minigraph alignment with seq of length {}",
        repr.len()
    );

    // call minigraph and parse the output
    let mut graph: Graph = align_seq_graph(&repr, b);
    eprintln!("Minigraph is done building the graph");

    let temp_dir = tempfile::tempdir().unwrap();

    // write representative sequence to file
    let repr_filepath = temp_dir.path().join("repr.fa");
    let mut repr_file = File::create(&repr_filepath).unwrap();
    write!(repr_file, ">query\n{}\n", repr).unwrap();
    repr_file.flush().unwrap();

    // write graph to file
    let graph_filepath = temp_dir.path().join("graph.gfa");
    let mut graph_file = File::create(&graph_filepath).unwrap();
    graph.write(&mut graph_file).unwrap();
    graph_file.flush().unwrap();

    let out_path = temp_dir.path().join("out.gaf");

    let output = Command::new("GraphAligner")
        .arg("-g")
        .arg(&graph_filepath)
        .arg("-f")
        .arg(&repr_filepath)
        .arg("-a")
        .arg(&out_path)
        .args(["-x", "vg"])
        .output()
        .expect("failed to execute GraphAligner");

    eprintln!(
        "GraphAligner stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    eprintln!("GraphAligner is done aligning");

    let gaf_file = File::open(out_path).unwrap();
    let gaf = gax::gaf::parse(gaf_file).unwrap();
    assert!(!gaf.is_empty());

    if gaf.iter().any(|r| r.strand == '-') {
        panic!("GraphAligner produced reverse strand alignments, which is not supported yet.");
    }

    // where is the representitive sequence mapped in the old graph?
    let mapping_repr_to_old: Vec<(&[u8], usize)> = build_mapping_seq_to_old(a, repr_path, &repr);

    // println!("old mapping:");
    // for (s, i) in mapping_repr_to_old.iter() {
    //     let s_str = std::str::from_utf8(s).unwrap();
    //     println!("  ({},{})", s_str, i);
    // }

    let mapping_repr_to_new: Vec<Option<(Vec<u8>, usize)>> =
        build_mapping_seq_to_new(repr, &graph, gaf)
            .into_iter()
            .map(|x| x.map(|(s, o)| (s.to_vec(), o)))
            .collect();

    // println!("new mapping:");
    // for x in mapping_repr_to_new.iter() {
    //     match x {
    //         Some((s, i)) => {
    //             let s_str = std::str::from_utf8(s).unwrap();
    //             println!("  ({},{})", s_str, i);
    //         }
    //         None => {
    //             println!("  None");
    //         }
    //     }
    // }

    merge_graph_a_into_b(a, &mut graph, &mapping_repr_to_old, &mapping_repr_to_new);

    graph.mass_rename();
    graph
}

fn merge_graph_a_into_b(
    old_graph: &Graph,
    new_graph: &mut Graph,
    mapping_repr_to_old: &[(&[u8], usize)],
    mapping_repr_to_new: &[Option<(Vec<u8>, usize)>],
) {
    assert_eq!(mapping_repr_to_old.len(), mapping_repr_to_new.len());

    // A vertex → B vertex
    let mut a_to_b: HashMap<&[u8], &[u8]> = HashMap::new();

    // get vertex mapping from aligned positions
    for ((a_vid, _), b_vid) in mapping_repr_to_old.iter().zip(mapping_repr_to_new.iter()) {
        let Some((b_vid, _)) = b_vid else {
            continue;
        };
        a_to_b.entry(*a_vid).or_insert(b_vid);
    }

    // create missing vertices in B
    for a_segment in old_graph.segments.iter() {
        if !a_to_b.contains_key(&a_segment.name.as_ref()) {
            let mut b_segment = a_segment.clone();
            b_segment.name.extend_from_slice(b"_a");
            assert!(new_graph.segments.iter().all(|s| s.name != b_segment.name));
            new_graph.segments.push(b_segment);
            a_to_b.insert(a_segment.name.as_ref(), a_segment.name.as_ref());
        }
    }

    // copy edges
    let mut seen_edges = HashSet::new();
    for edge in &old_graph.links {
        let Some(from_b) = a_to_b.get(&edge.from_segment.as_ref()) else {
            continue;
        };
        let Some(to_b) = a_to_b.get(&edge.to_segment.as_ref()) else {
            continue;
        };

        // Avoid duplicate edges
        if seen_edges.insert((from_b, to_b, edge.overlap.clone())) {
            let link = Link {
                from_segment: from_b.to_vec(),
                from_orient: edge.from_orient,
                to_segment: to_b.to_vec(),
                to_orient: edge.to_orient,
                overlap: edge.overlap.clone(),
                optional: Vec::new(),
            };
            new_graph.links.push(link);
        }
    }
}

fn build_mapping_seq_to_new(
    repr: String,
    graph: &Graph,
    gaf: Vec<GafRecord>,
) -> Vec<Option<(&[u8], usize)>> {
    // build `mapping_from_repr_to_new` using the GAF
    let mut mapping_from_repr_to_new: Vec<Option<(&[u8], usize)>> = vec![None; repr.len()];
    let name_map_new = NameMap::build_from_gfa(graph);
    for record in gaf {
        let mut cigar: CIGAR = parse_cigar(&record.opt_fields["cg"].1);
        for step in record.path {
            // build mapping for this step using current_cigar
            let seq_idx = name_map_new
                .map_name(&step.name)
                .expect("Name not found in new graph");
            let seq = graph.segments[seq_idx].sequence.as_slice();

            let mut seq_pos = 0;
            let mut repr_pos = record.query_start as usize;

            let step_length: usize = match (step.start, step.end) {
                (Some(start), Some(end)) => (end - start) as usize, // segId
                _ => seq.len(),                                     // stableIntv or stableId
            };

            let (current_cigar, rest) = cigar.split_at(step_length);
            cigar = rest;
            for (len, op) in current_cigar.iter() {
                let len = len as usize;
                match op {
                    CIGAROp::M | CIGAROp::X => {
                        // Match | Mismatch
                        for _ in 0..len {
                            mapping_from_repr_to_new[repr_pos] = Some((seq, seq_pos));
                            repr_pos += 1;
                            seq_pos += 1;
                        }
                    }
                    CIGAROp::I => {
                        // Insertion
                        repr_pos += len;
                    }
                    CIGAROp::D => {
                        // Deletion
                        seq_pos += len;
                    }
                    CIGAROp::N => {
                        // Reference skip
                        seq_pos += len;
                    }
                    CIGAROp::S => {
                        // Soft clip
                        repr_pos += len;
                    }
                    CIGAROp::H | CIGAROp::P | CIGAROp::E => {
                        // Hard clip | Padding | Extension
                    }
                }
            }
        }
        // assert!(cigar.is_empty(), "CIGAR not fully consumed");
        if !cigar.is_empty() {
            eprintln!(
                "Warning: CIGAR not fully consumed for one GAF record. Cigar left={:?}",
                cigar
            );
        }
    }
    mapping_from_repr_to_new
}

fn build_mapping_seq_to_old<'graph>(
    a: &'graph Graph,
    repr_path: Vec<(&'graph [u8], bool)>,
    repr: &String,
) -> Vec<(&'graph [u8], usize)> {
    let name_map_old = NameMap::build_from_gfa(a);
    let mut mapping_from_repr_to_old: Vec<(&[u8], usize)> = Vec::new();
    for &(name, direction) in &repr_path {
        let seq_indx = name_map_old.map_name(name).expect("Name not found");
        let seq = a.segments[seq_indx].sequence.as_slice();
        let range = 0..seq.len();
        if direction {
            for i in range {
                mapping_from_repr_to_old.push((seq, i));
            }
        } else {
            for i in range.rev() {
                mapping_from_repr_to_old.push((seq, i));
            }
        }
    }
    assert_eq!(repr.len(), mapping_from_repr_to_old.len());
    mapping_from_repr_to_old
}

pub fn align_seq_graph(a: &str, b: &Graph) -> Graph {
    let mut query = NamedTempFile::new().unwrap();

    write!(query, ">query\n{}", a).unwrap();

    let mut graph = NamedTempFile::new().unwrap();
    b.write(&mut graph).unwrap();

    query.flush().unwrap();
    graph.flush().unwrap();

    let out = Command::new("minigraph")
        .arg("-cxggs")
        .arg(graph.path())
        .arg(query.path())
        .output()
        .expect("failed to execute minigraph");

    eprintln!("Minigraph stderr: {}", String::from_utf8_lossy(&out.stderr));

    Graph::from_bytes(&out.stdout[..])
}

pub fn align_seq_seq(a: &str, b: &str) -> Graph {
    let graph = Graph {
        segments: vec![Segment {
            name: "seq".into(),
            sequence: a.as_bytes().to_vec(),
            optional: Vec::new(),
        }],
        links: Vec::new(),
        paths: Vec::new(),
        containments: Vec::new(),
        header: Header::default(),
    };

    // Align `b` onto this one-segment graph
    align_seq_graph(b, &graph)
}

fn parse_cigar(cigar: &str) -> CIGAR {
    CIGAR::from_bytestring(cigar.as_bytes()).expect("Failed to parse CIGAR")
}
