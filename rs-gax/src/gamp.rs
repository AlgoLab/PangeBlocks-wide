use crate::framing::{self, vg, FramingError};
use std::{
    fs::File,
    io::{Read, Write},
};

pub fn parse(data: impl Read) -> Result<Vec<vg::MultipathAlignment>, FramingError> {
    framing::parse::<vg::MultipathAlignment>(data)
}

pub fn parse_from_file(
    path: impl AsRef<std::path::Path>,
) -> Result<Vec<vg::MultipathAlignment>, FramingError> {
    let f = File::open(path)?;
    parse(f)
}

pub fn write(alignments: &[vg::MultipathAlignment], mut out_file: impl Write) -> Result<(), FramingError> {
    framing::write::<vg::MultipathAlignment>(alignments, &mut out_file)
}

pub fn write_to_file(
    alignments: &[vg::MultipathAlignment],
    path: impl AsRef<std::path::Path>,
) -> Result<(), FramingError> {
    let f = File::create(path)?;
    write(alignments, f)
}

impl From<vg::Alignment> for vg::MultipathAlignment {
    fn from(value: vg::Alignment) -> Self {
        Self {
            sequence: value.sequence,
            quality: value.quality,
            name: value.name,
            sample_name: value.sample_name,
            read_group: value.read_group,
            subpath: vec![vg::Subpath {
                path: value.path,
                score: value.score,
                ..Default::default()
            }],
            mapping_quality: value.mapping_quality,
            annotation: value.annotation,
            ..Default::default()
        }
    }
}

#[cfg(test)]
mod tests {
    use prost_types::{value::Kind, Value};

    use super::*;
    use std::fs::File;

    #[test]
    fn gamp_read() {
        let in_file = "data/example.gamp";
        let f = File::open(in_file).unwrap();
        let alignments: Vec<vg::MultipathAlignment> = parse(f).unwrap();
        let first = alignments[0].clone();

        let name = "FBtr0342963_e_1536_X_294766";
        let mapping_quality = 60;
        let sequence = "TAGATAAAAAATAAACGGAAAATTTGTTATTTCTTTCGTACATGGTAAAGAATCTTTTTTTACTTGTGTTTCTGTGATTTGAGTGTTTGAAAAATTTAAC";
        let start = vec![0];

        assert_eq!(first.name, name);
        assert_eq!(first.mapping_quality, mapping_quality);
        assert_eq!(first.sequence, sequence);
        assert_eq!(first.start, start);
    }

    #[test]
    fn gamp_write() {
        let out_file = "data/example.out.gamp";
        let of = File::create(out_file).unwrap();
        let alignment = vg::MultipathAlignment {
            name: "test".into(),
            mapping_quality: 1000,
            sequence: "AAAAATAAACGG".into(),
            start: vec![0],
            ..Default::default()
        };
        let alignments: Vec<vg::MultipathAlignment> = vec![alignment.clone()];
        write(&alignments, of).unwrap();

        let in_file = "data/example.out.gamp";
        let f = File::open(in_file).unwrap();
        let alignments: Vec<vg::MultipathAlignment> = parse(f).unwrap();
        let first = alignments[0].clone();

        assert_eq!(first, alignment);
    }

    #[test]
    fn gamp_edit() {
        let in_file = "data/example.gamp";
        let f = File::open(in_file).unwrap();

        let alignments: Vec<vg::MultipathAlignment> = parse(f).unwrap();
        let mut alignment = alignments[0].clone();

        alignment.name = "new_name".into();
        alignment.annotation.as_mut().unwrap().fields.insert(
            "new_value".into(),
            Value {
                kind: Some(Kind::NumberValue(10.0)),
            },
        );

        let out_file = "data/example.out.gamp";
        let of = File::create(out_file).unwrap();
        write(&(vec![alignment.clone()]), of).unwrap();

        let in_file = "data/example.out.gamp";
        let f = File::open(in_file).unwrap();
        let alignments: Vec<vg::MultipathAlignment> = parse(f).unwrap();
        let first = alignments[0].clone();

        assert_eq!(first, alignment);
    }
}
