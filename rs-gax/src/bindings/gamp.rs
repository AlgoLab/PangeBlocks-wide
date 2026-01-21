use super::dict::{pydict_to_struct, struct_to_pydict};
use crate::{gamp, vg};
use pyo3::prelude::*;

// IntoPyObject
impl IntoPy<PyObject> for vg::MultipathAlignment {
    fn into_py(self, py: Python<'_>) -> PyObject {
        let def = py.import("gax").unwrap();
        let o = def.getattr("MultipathAlignment").unwrap().call0().unwrap();
        o.setattr("sequence", self.sequence).unwrap();
        o.setattr("quality", self.quality).unwrap();
        o.setattr("name", self.name).unwrap();
        o.setattr("sample_name", self.sample_name).unwrap();
        o.setattr("read_group", self.read_group).unwrap();
        o.setattr("subpath", self.subpath.into_py(py)).unwrap();
        o.setattr("mapping_quality", self.mapping_quality).unwrap();
        o.setattr("start", self.start).unwrap();
        o.setattr("paired_read_name", self.paired_read_name)
            .unwrap();
        o.setattr(
            "annotation",
            self.annotation.map_or_else(
                || py.None(),
                |inner| struct_to_pydict(py, &inner).unwrap().into_py(py),
            ),
        )
        .unwrap();
        o.into_py(py)
    }
}

impl IntoPy<PyObject> for vg::Subpath {
    fn into_py(self, py: Python<'_>) -> PyObject {
        let def = py.import("gax").unwrap();
        let o = def.getattr("Subpath").unwrap().call0().unwrap();
        o.setattr("path", self.path.into_py(py)).unwrap();
        o.setattr("next", self.next.into_py(py)).unwrap();
        o.setattr("score", self.score).unwrap();
        o.setattr("connection", self.connection.into_py(py))
            .unwrap();
        o.into_py(py)
    }
}

impl IntoPy<PyObject> for vg::Connection {
    fn into_py(self, py: Python<'_>) -> PyObject {
        let def = py.import("gax").unwrap();
        let o = def.getattr("Connection").unwrap().call0().unwrap();
        o.setattr("next", self.next).unwrap();
        o.setattr("score", self.score).unwrap();
        o.into_py(py)
    }
}

// FromPyObject
impl FromPyObject<'_> for vg::MultipathAlignment {
    fn extract(ob: &PyAny) -> PyResult<Self> {
        Ok(Self {
            sequence: ob.getattr("sequence")?.extract()?,
            quality: ob.getattr("quality")?.extract()?,
            name: ob.getattr("name")?.extract()?,
            sample_name: ob.getattr("sample_name")?.extract()?,
            read_group: ob.getattr("read_group")?.extract()?,
            subpath: ob.getattr("subpath")?.extract()?,
            mapping_quality: ob.getattr("mapping_quality")?.extract()?,
            start: ob.getattr("start")?.extract()?,
            paired_read_name: ob.getattr("paired_read_name")?.extract()?,
            annotation: Some(pydict_to_struct(ob.getattr("annotation")?.extract()?)?),
        })
    }
}

impl FromPyObject<'_> for vg::Subpath {
    fn extract(ob: &PyAny) -> PyResult<Self> {
        Ok(Self {
            path: ob.getattr("path")?.extract()?,
            next: ob.getattr("next")?.extract()?,
            score: ob.getattr("score")?.extract()?,
            connection: ob.getattr("connection")?.extract()?,
        })
    }
}

impl FromPyObject<'_> for vg::Connection {
    fn extract(ob: &PyAny) -> PyResult<Self> {
        Ok(Self {
            next: ob.getattr("next")?.extract()?,
            score: ob.getattr("score")?.extract()?,
        })
    }
}

#[pyfunction]
fn parse(file_name: &str) -> PyResult<Vec<PyObject>> {
    let gamp = gamp::parse_from_file(file_name)?;
    Python::with_gil(|py| {
        let gamp: Vec<_> = gamp.iter().map(|o| o.clone().into_py(py)).collect();
        Ok(gamp)
    })
}

#[pyfunction]
fn write(gamp: Vec<PyObject>, file_name: &str) -> PyResult<()> {
    Python::with_gil(|py| -> PyResult<_> {
        let records = gamp
            .iter()
            .map(|o| -> PyResult<_> { o.extract::<vg::MultipathAlignment>(py) })
            .collect::<PyResult<Vec<_>>>()?;
        gamp::write_to_file(&records, file_name)?;
        Ok(())
    })
}

pub(crate) fn submodule(py: Python<'_>) -> PyResult<&PyModule> {
    let module = PyModule::new(py, "gamp")?;
    module.add_function(wrap_pyfunction!(parse, module)?)?;
    module.add_function(wrap_pyfunction!(write, module)?)?;
    Ok(module)
}
