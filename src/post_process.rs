use std::{fs::File, io::Write, process::Command};

use bio::io::fasta::Record;
use gfa::gfa::{Link, Orientation, Segment};

use crate::{Graph, GraphExt};

pub fn post_process_graph(mut graph: Graph, records: &[Record]) -> Graph {
    let temp_dir = tempfile::tempdir().unwrap();

    let graph_path = temp_dir.path().join("graph.gfa");
    let mut graph_file = File::create(&graph_path).unwrap();
    graph.write(&mut graph_file).unwrap();
    graph_file.flush().unwrap();

    for record in records {
        let out_path = temp_dir.path().join("out.gaf");
        let repr_path = temp_dir.path().join("repr.fa");

        let mut repr_file = File::create(&repr_path).unwrap();
        write!(
            repr_file,
            ">{}\n{}\n",
            record.id(),
            std::str::from_utf8(record.seq()).unwrap()
        )
        .unwrap();
        repr_file.flush().unwrap();

        Command::new("/home/matteo/miniconda3/bin/GraphAligner")
            .args(["--seeds-minimizer-ignore-frequent", "0.0"])
            .args(["--seeds-minimizer-density", "-1"])
            .args(["--seeds-extend-density", "-1"])
            .args(["--seeds-mxm-length", "15"])
            .args(["--bandwidth", "200"])
            .args(["--tangle-effort", "-1"])
            .args(["--precise-clipping", "0.9"])
            .arg("-g")
            .arg(&graph_path)
            .arg("-f")
            .arg(&repr_path)
            .arg("-a")
            .arg(&out_path)
            .args(["-x", "vg"])
            .output()
            .expect("failed to execute GraphAligner");

        let gaf_file = File::open(out_path).unwrap();
        let mut gaf = gax::gaf::parse(gaf_file).unwrap();
        assert!(!gaf.is_empty());

        // find all holes between alignments
        gaf.sort_by_key(|aln| aln.query_start);
        for window in gaf.windows(2) {
            let first = &window[0];
            let second = &window[1];
            if first.query_end < second.query_start {
                // hole found
                // check if there is already a link from first to second
                let from = first.path.last().unwrap();
                let to = second.path.first().unwrap();
                let from_orient = if from.is_reverse {
                    Orientation::Backward
                } else {
                    Orientation::Forward
                };
                let to_orient = if to.is_reverse {
                    Orientation::Backward
                } else {
                    Orientation::Forward
                };
                let from_segment = from.name.as_bytes();
                let to_segment = from.name.as_bytes();

                if graph.links.iter().any(|link| {
                    link.from_segment == from_segment
                        && link.from_orient == from_orient
                        && link.to_segment == to_segment
                        && link.to_orient == to_orient
                }) {
                    continue;
                }

                // get missing sequence
                let missing_seq =
                    &record.seq()[first.query_end as usize..second.query_start as usize];

                let missing_segment_name =
                    format!("missing_{}_{}", first.query_end, second.query_start);
                let missing_segment = Segment::new(missing_segment_name.as_bytes(), missing_seq);
                graph.segments.push(missing_segment.clone());

                // add links: from -> missing -> to
                graph.links.push(Link::new(
                    from_segment,
                    from_orient,
                    missing_segment_name.as_bytes(),
                    Orientation::Forward,
                    b"0M",
                ));

                graph.links.push(Link::new(
                    missing_segment_name.as_bytes(),
                    Orientation::Forward,
                    to_segment,
                    to_orient,
                    b"0M",
                ));
            }
        }
    }

    graph
}
