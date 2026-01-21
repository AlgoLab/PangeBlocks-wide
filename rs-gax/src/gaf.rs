use crate::{graph::GFAExt, vg, ConversionError};
use gfa::gfa::GFA;
use prost_types::value::Kind;
use pyo3::FromPyObject;
use std::{
    collections::HashMap,
    fs::File,
    io::{Read, Write},
};

#[derive(thiserror::Error, Debug)]
#[error(transparent)]
pub enum GafError {
    Io(#[from] std::io::Error),
    ParseInt(#[from] std::num::ParseIntError),
    #[error("Missing start in interval step")]
    MissingStart,
    #[error("Missing end in interval step")]
    MissingEnd,
    #[error("Not enough tokens in line")]
    MissingToken,
}

pub fn parse(data: impl Read) -> Result<Vec<GafRecord>, GafError> {
    let mut string = String::new();
    let mut reader = std::io::BufReader::new(data);
    reader.read_to_string(&mut string)?;
    string.lines().map(|line| GafRecord::parse(line)).collect()
}

pub fn parse_from_file(path: impl AsRef<std::path::Path>) -> Result<Vec<GafRecord>, GafError> {
    let f = File::open(path)?;
    parse(f)
}

pub fn write(records: &Vec<GafRecord>, mut out_file: impl Write) -> Result<(), GafError> {
    for record in records {
        record.write(&mut out_file)?;
    }
    Ok(())
}

pub fn write_to_file(
    records: &Vec<GafRecord>,
    path: impl AsRef<std::path::Path>,
) -> Result<(), GafError> {
    let f = File::create(path)?;
    write(records, f)
}

/**
 * We allow pretty much any field to be set as "*" in the GAF, which gets mapped to a -1
 * in the numeric fields (as there are no valid negative values)
 */
const MISSING_INT: i64 = -1;
const MISSING_STRING: &str = "*";

/**
 * One step of a GAF path
 */
#[derive(Debug, PartialEq, Eq, Clone, FromPyObject)]
pub struct GafStep {
    pub name: String,       // Either a path name or segment/node name (see above)
    pub is_reverse: bool,   // In reverse orientation ('<' in GAF)
    pub is_stable: bool,    // True if it's a stable path name (as opposed to segment/node name)
    pub is_interval: bool, // True if it's an interval of a stable path (false if it's the whole path)
    pub start: Option<i64>, // 0-based start (inclusive). only defined if is_stable and is_interval are true
    pub end: Option<i64>, // 0-based end (inclusive). only defined if is_stable and is_interval are true
}

impl GafStep {
    /*
     * Write a GAF Step to a stream
     */
    fn write(&self, mut f: impl Write) -> Result<(), GafError> {
        if !self.is_stable || self.is_interval {
            if self.is_reverse {
                write!(f, "<")?;
            } else {
                write!(f, ">")?;
            }
        }
        write!(f, "{}", self.name)?;
        if self.is_interval {
            write!(
                f,
                ":{}-{}",
                self.start.ok_or(GafError::MissingStart)?,
                self.end.ok_or(GafError::MissingEnd)?
            )?;
        }
        Ok(())
    }
}

/**
 * One line of GAF as described here: https://github.com/lh3/gfatools/blob/master/doc/rGFA.md
 */

#[derive(Debug, Clone, PartialEq, Eq, FromPyObject, Default)]
pub struct GafRecord {
    pub query_name: String, // Query sequence name
    pub query_length: i64,  // Query sequence length
    pub query_start: i64,   // 0-based, closed
    pub query_end: i64,     // 0-based, open
    pub path_length: i64,
    pub path_start: i64,    // Start position on the path (0-based)
    pub path_end: i64,      // End position on the path (0-based)
    pub matches: i64,       // Number of residue matches
    pub block_length: i64,  // Alignment block length
    pub mapq: i32,          // Mapping quality (0-255; 255 for missing)
    pub strand: char,       // strand relative to the path + or -
    pub path: Vec<GafStep>, // the path

    // Map a tag name to its type and value
    // ex: "de:f:0.2183" in the GAF would appear as opt_fields["de"] = ("f", "0.2183")
    pub opt_fields: HashMap<String, (String, String)>,
}

fn number_or_missing(token: &str) -> Result<i64, GafError> {
    Ok(if token == "*" {
        MISSING_INT
    } else {
        token.parse::<i64>()?
    })
}

impl GafRecord {
    /**
     * Parse a single GAF record
     */
    pub fn parse(line: &str) -> Result<Self, GafError> {
        let mut split = line.split('\t');

        let mut token = split.next().ok_or(GafError::MissingToken)?;
        let query_name = token.to_string();
        token = split.next().ok_or(GafError::MissingToken)?;
        let query_length = number_or_missing(token)?;
        token = split.next().ok_or(GafError::MissingToken)?;
        let query_start = number_or_missing(token)?;
        token = split.next().ok_or(GafError::MissingToken)?;
        let query_end = number_or_missing(token)?;
        token = split.next().ok_or(GafError::MissingToken)?;
        let strand = token.chars().next().unwrap();

        token = split.next().ok_or(GafError::MissingToken)?;
        let mut path = Vec::new();
        if token.to_string().starts_with(['<', '>']) {
            // orientIntv
            let mut splits: Vec<_> = token.match_indices(['<', '>']).map(|(i, _)| i).collect();
            splits.push(token.len());

            for step_token in splits
                .windows(2)
                .map(|indexes| &token[indexes[0]..indexes[1]])
            {
                let is_reverse = &step_token[0..1] == "<";
                let s = match step_token.find(':') {
                    Some(colon) => {
                        let Some(dash) = step_token[colon..].find('-') else {
                            panic!("Error parsing GAF range of {}", step_token)
                        };
                        let start = step_token[colon + 1..colon + dash].parse::<i64>()?;
                        let end = step_token[colon + 1 + dash..].parse::<i64>()?;
                        // stableIntv
                        GafStep {
                            name: (step_token[1..colon - 1]).to_string(),
                            is_reverse,
                            is_stable: true,
                            is_interval: true,
                            start: Some(start),
                            end: Some(end),
                        }
                    }
                    None => {
                        // segId
                        GafStep {
                            name: (step_token[1..]).to_string(),
                            is_reverse,
                            is_stable: false,
                            is_interval: false,
                            start: None,
                            end: None,
                        }
                    }
                };
                path.push(s);
            }
        } else {
            // stableId
            path.push(GafStep {
                name: token.to_string(),
                is_reverse: false,
                is_stable: true,
                is_interval: false,
                start: None,
                end: None,
            });
        }

        token = split.next().ok_or(GafError::MissingToken)?;
        let path_length = number_or_missing(token)?;
        token = split.next().ok_or(GafError::MissingToken)?;
        let path_start = number_or_missing(token)?;
        token = split.next().ok_or(GafError::MissingToken)?;
        let path_end = number_or_missing(token)?;
        token = split.next().ok_or(GafError::MissingToken)?;
        let matches = number_or_missing(token)?;
        token = split.next().ok_or(GafError::MissingToken)?;
        let block_length = number_or_missing(token)?;

        token = split.next().ok_or(GafError::MissingToken)?;
        let mapq = if token == MISSING_STRING {
            MISSING_INT as _
        } else {
            token.parse::<i32>()?
        };

        let mut opt_fields = HashMap::new();
        for opt_token in split.collect::<Vec<&str>>() {
            let key: &str = &opt_token[..2];
            let typ: &str = &opt_token[3..4];
            let value: &str = &opt_token[5..];
            if opt_fields
                .insert(key.to_string(), (typ.to_string(), value.to_string()))
                .is_some()
            {
                panic!("Duplicate optional field found: {}", key)
            }
        }

        Ok(Self {
            query_name,
            query_length,
            query_start,
            query_end,
            path_length,
            path_start,
            path_end,
            matches,
            block_length,
            mapq,
            strand,
            path,
            opt_fields,
        })
    }

    /**
     * Write a GAF record to a stream
     */
    pub fn write(&self, mut f: impl Write) -> Result<(), GafError> {
        if self.query_name.is_empty() {
            write!(f, "{}\t", MISSING_STRING)?;
        } else {
            write!(f, "{}\t", self.query_name)?;
        }
        write!(
            f,
            "{}\t{}\t{}\t{}\t",
            self.query_length, self.query_start, self.query_end, self.strand
        )?;
        if self.path.is_empty() {
            write!(
                f,
                "{}\t{}\t{}\t{}\t{}\t{}",
                MISSING_STRING,
                MISSING_STRING,
                MISSING_STRING,
                MISSING_STRING,
                MISSING_STRING,
                MISSING_STRING
            )?;
        } else {
            for step in self.path.iter() {
                step.write(&mut f)?;
            }
            write!(
                f,
                "\t{}\t{}\t{}\t{}\t{}\t",
                self.path_length, self.path_start, self.path_end, self.matches, self.block_length
            )?
        }

        if self.mapq as i64 == MISSING_INT {
            write!(f, "{}", 255)?;
        } else {
            write!(f, "{}", self.mapq)?;
        }
        for (key, value) in &self.opt_fields {
            write!(f, "\t{}:{}:{}", key, value.0, value.1)?;
        }
        writeln!(f)?;
        Ok(())
    }

    pub fn iter_cigar(&self) -> Vec<Cigar> {
        if self.opt_fields.contains_key("cs") {
            self.iter_cs()
                .into_iter()
                .map(|cs| {
                    let first_char = cs.chars().next().unwrap();
                    match first_char {
                        ':' => Cigar {
                            cat: first_char,
                            length: cs[1..].parse::<usize>().unwrap(),
                            query: "".into(),
                            target: "".into(),
                        },
                        '+' => {
                            let query = &cs[1..];
                            Cigar {
                                cat: first_char.into(),
                                length: query.len(),
                                query: query.into(),
                                target: "".into(),
                            }
                        }
                        '-' => {
                            let target = &cs[1..];
                            Cigar {
                                cat: first_char.into(),
                                length: target.len(),
                                query: "".into(),
                                target: target.into(),
                            }
                        }
                        '*' => Cigar {
                            cat: first_char.into(),
                            length: 1,
                            query: cs[2..3].into(),
                            target: cs[1..2].into(),
                        },
                        _ => unreachable!(),
                    }
                })
                .collect()
        } else {
            self.iter_cg()
        }
    }

    pub fn iter_cs(&self) -> Vec<&str> {
        let Some(cigar_pair) = self.opt_fields.get("cs") else { return vec![]; };
        let cs_cigar = &cigar_pair.1;
        let mut splits = cs_cigar
            .match_indices([':', '*', '-', '+'])
            .map(|(i, _)| i)
            .collect::<Vec<_>>();
        splits.push(cs_cigar.len());
        splits
            .windows(2)
            .map(|indexes| &cs_cigar[indexes[0]..indexes[1]])
            .collect()
    }

    pub fn iter_cg(&self) -> Vec<Cigar> {
        let Some(cigar_pair) = self.opt_fields.get("cg") else { return vec![]; };
        let cg_cigar = &cigar_pair.1;
        let mut splits = cg_cigar
            .match_indices(['M', 'I', 'D', 'N', 'S', 'H', 'P', 'X', '='])
            .map(|(i, _)| i)
            .collect::<Vec<_>>();
        splits.push(cg_cigar.len());
        splits
            .windows(2)
            .map(|indexes| {
                let cat = cg_cigar[indexes[0]..indexes[1]].chars().next().unwrap();
                Cigar {
                    length: 1,
                    cat,
                    query: "".into(),
                    target: "".into(),
                }
            })
            .collect()
    }

    pub fn convert_from_gam(
        value: &vg::Alignment,
        graph: &GFA<usize, ()>,
    ) -> Result<Self, ConversionError> {
        let mut query_name = value.name.clone();
        if query_name.is_empty() {
            query_name = MISSING_STRING.into();
        }

        let mut gaf = GafRecord {
            query_name,
            query_length: value.sequence.len() as _,
            mapq: value.mapping_quality,
            ..Default::default()
        };

        if let Some(path) = value.path.clone() {
            if !path.mapping.is_empty() {
                gaf.query_start = 0;
                gaf.query_end = value.sequence.len() as _;
                gaf.strand = '+';
                gaf.path_length = 0;
                gaf.path_start = 0; // missing
                gaf.matches = 0;

                let mut cs_cigar_str = "".to_string();
                let mut running_match_length = 0;
                let mut running_deletion = false;
                let mut total_to_len = 0;
                let mut prev_offset = 0;
                for (i, mapping) in path.mapping.iter().enumerate() {
                    let position = mapping
                        .position
                        .as_ref()
                        .ok_or(ConversionError::MissingPosition)?;
                    let start_offset_on_node = position.offset;
                    let mut offset = start_offset_on_node;
                    let node_to_segment_offset = 0;
                    let node_length = graph.node_to_length(position.node_id);
                    let mut node_seq = "".to_string();
                    let mut skip_step = false;
                    let mut _prev_range = (0, false, 0, 0);

                    if i > 0 && start_offset_on_node > 0 {
                        let prev_position = &path.mapping[i - 1]
                            .position
                            .as_ref()
                            .ok_or(ConversionError::MissingPosition)?;
                        if start_offset_on_node == prev_offset
                            && position.node_id == prev_position.node_id
                            && position.is_reverse == prev_position.is_reverse
                        {
                            skip_step = true;
                        } else {
                            if node_seq.is_empty() {
                                node_seq =
                                    graph.node_to_sequence(position.node_id, position.is_reverse);
                            }

                            let mut del_start_offset = 0;
                            if position.node_id == prev_position.node_id {
                                del_start_offset = prev_offset;
                            }
                            if start_offset_on_node > del_start_offset {
                                if running_match_length > 0 {
                                    cs_cigar_str += ":";
                                    cs_cigar_str += &running_match_length.to_string();
                                }
                                if !running_deletion {
                                    cs_cigar_str += "-";
                                }
                                cs_cigar_str += &node_seq[del_start_offset as usize
                                    ..(start_offset_on_node - del_start_offset) as usize];
                                running_deletion = true;
                            }
                        }
                    }

                    for edit in &mapping.edit {
                        if edit.is_match() {
                            gaf.matches += edit.from_length as i64;
                            running_match_length += edit.from_length as i64;
                            running_deletion = false;
                        } else {
                            if running_match_length > 0 {
                                cs_cigar_str += ":";
                                cs_cigar_str += &running_match_length.to_string();
                                running_match_length = 0;
                            }
                            if edit.is_sub() {
                                if node_seq.is_empty() {
                                    node_seq = graph
                                        .node_to_sequence(position.node_id, position.is_reverse);
                                }

                                for i in 0..edit.from_length as i64 {
                                    cs_cigar_str += "*";
                                    cs_cigar_str +=
                                        &node_seq[(offset + i) as usize..(offset + i + 1) as usize];
                                    cs_cigar_str += &edit.sequence[i as usize..i as usize + 1];
                                }
                                running_deletion = false;
                            } else if edit.is_deletion() {
                                if node_seq.is_empty() {
                                    node_seq = graph
                                        .node_to_sequence(position.node_id, position.is_reverse);
                                }

                                if !running_deletion {
                                    cs_cigar_str += "-";
                                }

                                cs_cigar_str += &node_seq
                                    [offset as usize..(offset + edit.from_length as i64) as usize];
                                running_deletion = true;
                            } else if edit.is_insertion() {
                                cs_cigar_str += "+";
                                cs_cigar_str += &edit.sequence;
                                running_deletion = false;
                            }
                        }
                        offset += edit.from_length as i64;
                        total_to_len += edit.to_length as i64;
                    }

                    // range
                    let range = (
                        position.node_id,
                        position.is_reverse,
                        start_offset_on_node,
                        offset - start_offset_on_node,
                    );

                    if i == 0 {
                        gaf.path_start = range.2;
                    } else if i + 1 == path.mapping.len()
                        && i > 0
                        && path.mapping[i].edit.len() == 1
                        && path.mapping[i].edit[0].is_insertion()
                    {
                        skip_step = true;
                    }

                    if i < path.mapping.len() - 1 && offset != node_length as i64 {
                        let next_position = path.mapping[i + 1]
                            .position
                            .as_ref()
                            .ok_or(ConversionError::MissingPosition)?;
                        if position.node_id != next_position.node_id
                            || position.is_reverse != next_position.is_reverse
                        {
                            if node_seq.is_empty() {
                                node_seq =
                                    graph.node_to_sequence(position.node_id, position.is_reverse);
                            }
                            if running_match_length > 0 {
                                cs_cigar_str += ":";
                                cs_cigar_str += &running_match_length.to_string();
                                running_match_length = 0;
                            }
                            if !running_deletion {
                                cs_cigar_str += "-";
                            }
                            cs_cigar_str += &node_seq[offset as usize..];
                            running_deletion = true;
                        } else {
                            skip_step = true;
                        }
                    }

                    if !skip_step {
                        gaf.path_length += node_length as i64;

                        gaf.path.push(GafStep {
                            name: range.0.to_string(),
                            is_stable: false,
                            is_reverse: range.1,
                            is_interval: false,
                            start: None,
                            end: None,
                        });
                    }

                    if i == path.mapping.len() - 1 {
                        gaf.path_end = gaf.path_start;

                        let offset_on_path_visit = offset + node_to_segment_offset;
                        if gaf.path_length > offset_on_path_visit {
                            gaf.path_end =
                                gaf.path_length - 1 - (node_length as i64 - offset_on_path_visit);
                        }
                    }
                    _prev_range = range;
                    prev_offset = offset;
                }

                if running_match_length > 0 {
                    cs_cigar_str += ":";
                    cs_cigar_str += &running_match_length.to_string();
                }

                if gaf.query_length == 0 && total_to_len > 0 {
                    gaf.query_length = total_to_len;
                    gaf.query_end = total_to_len;
                }

                gaf.block_length = gaf.query_length.max(gaf.path_end - gaf.path_start);

                gaf.opt_fields
                    .insert("cs".to_string(), ("Z".to_string(), cs_cigar_str));

                if value.identity > 0.0 {
                    let identity = ((1. - value.identity) * 10_000. + 0.5).floor() / 10000.;
                    gaf.opt_fields
                        .insert("dv".to_string(), ("f".to_string(), identity.to_string()));
                }

                if value.score > 0 {
                    gaf.opt_fields
                        .insert("AS".to_string(), ("i".to_string(), value.score.to_string()));
                }

                if !value.quality.is_empty() {
                    gaf.opt_fields.insert(
                        "bq".to_string(),
                        (
                            "Z".to_string(),
                            string_quality_short_to_char(&value.quality),
                        ),
                    );
                }

                if let Some(annotation) = value.annotation.clone() {
                    if annotation.fields.contains_key("proper_pair") {
                        if let Some(Kind::BoolValue(is_properly_paired)) =
                            annotation.fields["proper_pair"].kind.as_ref()
                        {
                            gaf.opt_fields.insert(
                                "pd".to_string(),
                                ("b".to_string(), is_properly_paired.to_string()),
                            );
                        }
                        if let Some(Kind::StringValue(support)) =
                            annotation.fields["support"].kind.as_ref()
                        {
                            gaf.opt_fields
                                .insert("AD".to_string(), ("i".to_string(), support.clone()));
                        }
                    }
                }
            }
        }
        if let Some(fragment_next) = value.fragment_next.clone() {
            gaf.opt_fields
                .insert("fn".to_string(), ("Z".to_string(), fragment_next.name));
        }
        if let Some(fragment_prev) = value.fragment_prev.clone() {
            gaf.opt_fields
                .insert("fp".to_string(), ("Z".to_string(), fragment_prev.name));
        }
        Ok(gaf)
    }
}

fn string_quality_short_to_char(quality: &[u8]) -> String {
    quality.iter().map(|byte| (byte + 33) as char).collect()
}

#[derive(Debug, Clone, Default)]
pub struct Cigar {
    pub cat: char,
    pub length: usize,
    pub query: String,
    pub target: String,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{convert_gam_to_gaf, gam};
    use gfa::parser::GFAParser;

    #[test]
    fn gaf_read() -> Result<(), Box<dyn std::error::Error>> {
        let line: &str = "read2\t7\t0\t7\t-\t>chr1:5-8>foo:8-16\t11\t1\t8\t7\t7\t60\tcg:Z:7M";
        let rec: GafRecord = GafRecord::parse(line)?;
        assert_eq!(rec.query_name, "read2");
        assert_eq!(rec.query_length, 7);
        assert_eq!(rec.strand, '-');
        Ok(())
    }

    #[test]
    fn convert_from_gam() -> Result<(), Box<dyn std::error::Error>> {
        use pretty_assertions::assert_eq;
        let graph: GFA<usize, ()> = GFAParser::new().parse_file("data/convert.gfa")?;
        let gam = gam::parse_from_file("data/convert.gam")?;
        let gaf = parse_from_file("data/convert.gaf")?;

        let generated_gaf = convert_gam_to_gaf(&gam, &graph)?;

        assert_eq!(gaf.len(), generated_gaf.len());

        let mut match_count = 0;
        let mut to_check = generated_gaf.clone();
        for item in &gaf {
            if let Some(index) = to_check.iter().position(|e| e == item) {
                to_check.remove(index);
                match_count += 1;
            }
        }
        assert_eq!(gaf.len(), match_count);
        Ok(())
    }
}
