use bytes::Buf;
use flate2::{read::MultiGzDecoder, Compression, GzBuilder};
use prost::encoding::{decode_varint, encode_varint};
use std::io::prelude::*;
use std::io::Cursor;

// Include the `vg` module, which is generated from vg.proto.
#[allow(clippy::all)]
pub mod vg {
    include!(concat!(env!("OUT_DIR"), "/vg.rs"));
}

const MAX_GROUP_SIZE: usize = 1000;

pub(crate) trait SupportedFormat: prost::Message + Default + Clone {
    fn type_tag() -> String;
}

impl SupportedFormat for vg::Alignment {
    fn type_tag() -> String {
        "GAM".to_string()
    }
}

impl SupportedFormat for vg::MultipathAlignment {
    fn type_tag() -> String {
        "MGAM".to_string()
    }
}

pub(crate) fn parse<Message: SupportedFormat>(data: impl Read) -> Result<Vec<Message>, FramingError> {
    // Decompress data
    let mut decoder = MultiGzDecoder::new(data);
    let mut data = vec![];
    decoder.read_to_end(&mut data)?;
    let mut cursor = Cursor::new(data);

    let mut alignments = vec![];
    while cursor.remaining() != 0 {
        // Read number of messages in the group
        let number_messages = decode_varint(&mut cursor)? - 1;

        // Read type tag
        let type_tag_len = decode_varint(&mut cursor)?;
        let mut type_tag = vec![0; type_tag_len as usize];
        cursor.read_exact(&mut type_tag)?;
        let type_tag = String::from_utf8(type_tag)?;
        // Should always be equal to GAM/MGAM
        if type_tag != Message::type_tag() {
            return Err(FramingError::InvalidTypeTag(type_tag, Message::type_tag()));
        }

        // Read all messages in the group
        for _ in 0..number_messages {
            let message_len = decode_varint(&mut cursor)?;
            let mut buffer = vec![0; message_len as _];
            cursor.read_exact(&mut buffer)?;
            let mut tmp = &buffer[..];
            let alignment = Message::decode(&mut tmp)?;
            alignments.push(alignment);
        }
    }
    Ok(alignments)
}

pub(crate) fn write<Message: SupportedFormat>(
    alignments: &[Message],
    mut out_file: impl Write,
) -> Result<(), FramingError> {
    let alignments = alignments.to_owned();
    let mut buf = vec![];

    for group in alignments.chunks(MAX_GROUP_SIZE) {
        // This step may be optional
        write_group(group, &mut buf)?;
    }

    // Compress data
    // FIXME For big files vg uses multi streams
    // This is not currently supported by flate2, see this PR:
    // https://github.com/rust-lang/flate2-rs/pull/325
    let mut encoder = GzBuilder::new().write(&mut out_file, Compression::new(9));
    encoder.write_all(&buf)?;
    Ok(())
}

fn write_group<Message: SupportedFormat>(
    alignments: &[Message],
    mut file: impl Write,
) -> Result<(), FramingError> {
    let mut buf = vec![];
    // Write number of messages in the group
    encode_varint(alignments.len() as u64 + 1, &mut buf);
    // Write type tag
    let type_tag = Message::type_tag();
    encode_varint(type_tag.len() as _, &mut buf);
    file.write_all(&buf)?;
    file.write_all(type_tag.as_bytes())?;

    // Write all messages
    for alignment in alignments {
        let mut buf = vec![];
        Message::encode(alignment, &mut buf)?;
        let mut buf_len = vec![];
        encode_varint(buf.len() as _, &mut buf_len);
        // Write message length
        file.write_all(&buf_len)?;
        // Write message
        file.write_all(&buf)?;
    }
    Ok(())
}

#[derive(thiserror::Error, Debug)]
#[error(transparent)]
pub enum FramingError {
    Io(#[from] std::io::Error),
    Utf8(#[from] std::string::FromUtf8Error),
    ProstDecode(#[from] prost::DecodeError),
    ProstEncode(#[from] prost::EncodeError),
    #[error("Type tag is {0}, expected \"{1}\"")]
    InvalidTypeTag(String, String),
}
