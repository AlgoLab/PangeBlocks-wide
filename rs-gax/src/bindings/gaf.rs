use crate::gaf;
use pyo3::prelude::*;

impl IntoPy<PyObject> for gaf::GafRecord {
    fn into_py(self, py: Python) -> PyObject {
        let def = py.import("gax").unwrap();
        let o = def.getattr("GafRecord").unwrap().call0().unwrap();
        o.setattr("query_name", self.query_name).unwrap();
        o.setattr("query_length", self.query_length).unwrap();
        o.setattr("query_start", self.query_start).unwrap();
        o.setattr("query_end", self.query_end).unwrap();
        o.setattr("path_length", self.path_length).unwrap();
        o.setattr("path_start", self.path_start).unwrap();
        o.setattr("path_end", self.path_end).unwrap();
        o.setattr("matches", self.matches).unwrap();
        o.setattr("block_length", self.block_length).unwrap();
        o.setattr("mapq", self.mapq).unwrap();
        o.setattr("strand", self.strand).unwrap();
        o.setattr("path", self.path.into_py(py)).unwrap();
        o.setattr("opt_fields", self.opt_fields).unwrap();
        o.into_py(py)
    }
}

impl IntoPy<PyObject> for gaf::GafStep {
    fn into_py(self, py: Python) -> PyObject {
        let def = py.import("gax").unwrap();
        let o = def.getattr("GafStep").unwrap().call0().unwrap();
        o.setattr("name", self.name).unwrap();
        o.setattr("is_reverse", self.is_reverse).unwrap();
        o.setattr("is_stable", self.is_stable).unwrap();
        o.setattr("is_interval", self.is_interval).unwrap();
        o.setattr("start", self.start).unwrap();
        o.setattr("end", self.end).unwrap();
        o.into_py(py)
    }
}

#[pyfunction]
fn parse(file_name: &str) -> PyResult<Vec<PyObject>> {
    let gaf = gaf::parse_from_file(file_name);
    Python::with_gil(|py| {
        let gaf: Vec<_> = gaf.iter().map(|o| o.clone().into_py(py)).collect();
        Ok(gaf)
    })
}

#[pyfunction]
fn write(gafs: Vec<PyObject>, file_name: &str) -> PyResult<()> {
    Python::with_gil(|py| -> PyResult<_> {
        let records = gafs
            .iter()
            .map(|o| -> PyResult<_> { o.extract::<gaf::GafRecord>(py) })
            .collect::<PyResult<_>>()?;
        gaf::write_to_file(&records, file_name)?;
        Ok(())
    })
}

pub(crate) fn submodule(py: Python<'_>) -> PyResult<&PyModule> {
    let module = PyModule::new(py, "gaf")?;
    module.add_function(wrap_pyfunction!(parse, module)?)?;
    module.add_function(wrap_pyfunction!(write, module)?)?;
    Ok(module)
}
