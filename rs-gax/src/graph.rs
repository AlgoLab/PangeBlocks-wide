use gfa::gfa::GFA;

pub trait GFAExt {
    fn node_to_length(&self, node_id: i64) -> usize;
    fn node_to_sequence(&self, node_id: i64, is_reverse: bool) -> String;
}

impl GFAExt for GFA<usize, ()> {
    fn node_to_length(&self, node_id: i64) -> usize {
        let node = self
            .segments
            .iter()
            .find(|n| n.name == node_id as usize)
            .unwrap();
        node.sequence.len()
    }

    fn node_to_sequence(&self, node_id: i64, is_reverse: bool) -> String {
        let node = self
            .segments
            .iter()
            .find(|n| n.name == node_id as usize)
            .unwrap()
            .sequence
            .clone();
        String::from_utf8(if is_reverse {
            node.into_iter().rev().collect()
        } else {
            node
        })
        .unwrap()
    }
}
