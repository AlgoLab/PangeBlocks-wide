use gfa::{gfa::GFA, parser::GFAParser};
use pyo3::{pyclass, pyfunction, PyResult};

#[pyclass(name = "Graph")]
pub struct GFAWrapper {
    pub graph: GFA<usize, ()>,
}

#[pyfunction]
pub(crate) fn load_graph(path: &str) -> PyResult<GFAWrapper> {
    Ok(GFAWrapper {
        graph: GFAParser::new().parse_file(path).unwrap(),
    })
}
