use super::dict::{pydict_to_struct, struct_to_pydict};
use crate::{gam, vg};
use pyo3::prelude::*;

// IntoPyObject
impl IntoPy<PyObject> for vg::Alignment {
    fn into_py(self, py: Python) -> PyObject {
        let def = py.import("gax").unwrap();
        let o = def.getattr("Alignment").unwrap().call0().unwrap();
        o.setattr("sequence", self.sequence).unwrap();
        o.setattr("path", self.path.into_py(py)).unwrap();
        o.setattr("name", self.name).unwrap();
        o.setattr("quality", self.quality).unwrap();
        o.setattr("mapping_quality", self.mapping_quality).unwrap();
        o.setattr("score", self.score).unwrap();
        o.setattr("query_position", self.query_position).unwrap();
        o.setattr("sample_name", self.sample_name).unwrap();
        o.setattr("read_group", self.read_group).unwrap();
        o.setattr(
            "fragment_prev",
            self.fragment_prev
                .map_or_else(|| py.None(), |inner| inner.into_py(py)),
        )
        .unwrap();
        o.setattr(
            "fragment_next",
            self.fragment_next
                .map_or_else(|| py.None(), |inner| inner.into_py(py)),
        )
        .unwrap();
        o.setattr("is_secondary", self.is_secondary).unwrap();
        o.setattr("identity", self.identity).unwrap();
        o.setattr("fragment", self.fragment.into_py(py)).unwrap();
        o.setattr("locus", self.locus.into_py(py)).unwrap();
        o.setattr("refpos", self.refpos.into_py(py)).unwrap();
        o.setattr("read_paired", self.read_paired).unwrap();
        o.setattr("read_mapped", self.read_mapped).unwrap();
        o.setattr("mate_unmapped", self.mate_unmapped).unwrap();
        o.setattr("read_on_reverse_strand", self.read_on_reverse_strand)
            .unwrap();
        o.setattr("mate_on_reverse_strand", self.mate_on_reverse_strand)
            .unwrap();
        o.setattr("soft_clipped", self.soft_clipped).unwrap();
        o.setattr("discordant_insert_size", self.discordant_insert_size)
            .unwrap();
        o.setattr("uniqueness", self.uniqueness).unwrap();
        o.setattr("correct", self.correct).unwrap();
        o.setattr("secondary_score", self.secondary_score).unwrap();
        o.setattr("fragment_score", self.fragment_score).unwrap();
        o.setattr(
            "mate_mapped_to_disjoint_subgraph",
            self.mate_mapped_to_disjoint_subgraph,
        )
        .unwrap();
        o.setattr(
            "fragment_length_distribution",
            self.fragment_length_distribution,
        )
        .unwrap();
        o.setattr("time_used", self.time_used).unwrap();
        o.setattr("to_correct", self.to_correct.into_py(py))
            .unwrap();
        o.setattr("correctly_mapped", self.correctly_mapped)
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

impl IntoPy<PyObject> for vg::Path {
    fn into_py(self, py: Python) -> PyObject {
        let def = py.import("gax").unwrap();
        let o = def.getattr("Path").unwrap().call0().unwrap();
        o.setattr("name", self.name).unwrap();
        o.setattr("mapping", self.mapping.into_py(py)).unwrap();
        o.setattr("is_circular", self.is_circular).unwrap();
        o.setattr("length", self.length).unwrap();
        o.into_py(py)
    }
}

impl IntoPy<PyObject> for vg::Locus {
    fn into_py(self, py: Python) -> PyObject {
        let def = py.import("gax").unwrap();
        let o = def.getattr("Locus").unwrap().call0().unwrap();
        o.setattr("name", self.name).unwrap();
        o.setattr("allele", self.allele.into_py(py)).unwrap();
        o.setattr("support", self.support.into_py(py)).unwrap();
        o.setattr("genotype", self.genotype.into_py(py)).unwrap();
        o.setattr("overall_support", self.overall_support.into_py(py))
            .unwrap();
        o.setattr("allele_log_likelihood", self.allele_log_likelihood)
            .unwrap();
        o.into_py(py)
    }
}

impl IntoPy<PyObject> for vg::Position {
    fn into_py(self, py: Python) -> PyObject {
        let def = py.import("gax").unwrap();
        let o = def.getattr("Position").unwrap().call0().unwrap();
        o.setattr("node_id", self.node_id).unwrap();
        o.setattr("offset", self.offset).unwrap();
        o.setattr("is_reverse", self.is_reverse).unwrap();
        o.setattr("name", self.name).unwrap();
        o.into_py(py)
    }
}

impl IntoPy<PyObject> for vg::Mapping {
    fn into_py(self, py: Python) -> PyObject {
        let def = py.import("gax").unwrap();
        let o = def.getattr("Mapping").unwrap().call0().unwrap();
        o.setattr("position", self.position.into_py(py)).unwrap();
        o.setattr("edit", self.edit.into_py(py)).unwrap();
        o.setattr("rank", self.rank).unwrap();
        o.into_py(py)
    }
}

impl IntoPy<PyObject> for vg::Support {
    fn into_py(self, py: Python) -> PyObject {
        let def = py.import("gax").unwrap();
        let o = def.getattr("Support").unwrap().call0().unwrap();
        o.setattr("quality", self.quality).unwrap();
        o.setattr("forward", self.forward).unwrap();
        o.setattr("reverse", self.reverse).unwrap();
        o.setattr("left", self.left).unwrap();
        o.setattr("right", self.right).unwrap();
        o.into_py(py)
    }
}

impl IntoPy<PyObject> for vg::Genotype {
    fn into_py(self, py: Python) -> PyObject {
        let def = py.import("gax").unwrap();
        let o = def.getattr("Genotype").unwrap().call0().unwrap();
        o.setattr("allele", self.allele).unwrap();
        o.setattr("is_phased", self.is_phased).unwrap();
        o.setattr("likelihood", self.likelihood).unwrap();
        o.setattr("log_likelihood", self.log_likelihood).unwrap();
        o.setattr("log_prior", self.log_prior).unwrap();
        o.setattr("log_posterior", self.log_posterior).unwrap();
        o.into_py(py)
    }
}

impl IntoPy<PyObject> for vg::Edit {
    fn into_py(self, py: Python) -> PyObject {
        let def = py.import("gax").unwrap();
        let o = def.getattr("Edit").unwrap().call0().unwrap();
        o.setattr("from_length", self.from_length).unwrap();
        o.setattr("to_length", self.to_length).unwrap();
        o.setattr("sequence", self.sequence).unwrap();
        o.into_py(py)
    }
}

// FromPyObject
impl FromPyObject<'_> for vg::Alignment {
    fn extract(ob: &PyAny) -> PyResult<Self> {
        let fragment_prev = ob.getattr("fragment_prev")?;
        let fragment_prev = if !fragment_prev.is_none() {
            Some(Box::new(fragment_prev.extract::<vg::Alignment>()?))
        } else {
            None
        };
        let fragment_next = ob.getattr("fragment_next")?;
        let fragment_next = if !fragment_next.is_none() {
            Some(Box::new(fragment_next.extract::<vg::Alignment>()?))
        } else {
            None
        };

        let annotation = ob.getattr("annotation")?;
        let annotation = if !annotation.is_none() {
            Some(pydict_to_struct(annotation.extract()?)?)
        } else {
            None
        };

        Ok(Self {
            sequence: ob.getattr("sequence")?.extract()?,
            path: ob.getattr("path")?.extract()?,
            name: ob.getattr("name")?.extract()?,
            quality: ob.getattr("quality")?.extract()?,
            mapping_quality: ob.getattr("mapping_quality")?.extract()?,
            score: ob.getattr("score")?.extract()?,
            query_position: ob.getattr("query_position")?.extract()?,
            sample_name: ob.getattr("sample_name")?.extract()?,
            read_group: ob.getattr("read_group")?.extract()?,
            fragment_prev,
            fragment_next,
            is_secondary: ob.getattr("is_secondary")?.extract()?,
            identity: ob.getattr("identity")?.extract()?,
            fragment: ob.getattr("fragment")?.extract()?,
            locus: ob.getattr("locus")?.extract()?,
            refpos: ob.getattr("refpos")?.extract()?,
            read_paired: ob.getattr("read_paired")?.extract()?,
            read_mapped: ob.getattr("read_mapped")?.extract()?,
            mate_unmapped: ob.getattr("mate_unmapped")?.extract()?,
            read_on_reverse_strand: ob.getattr("read_on_reverse_strand")?.extract()?,
            mate_on_reverse_strand: ob.getattr("mate_on_reverse_strand")?.extract()?,
            soft_clipped: ob.getattr("soft_clipped")?.extract()?,
            discordant_insert_size: ob.getattr("discordant_insert_size")?.extract()?,
            uniqueness: ob.getattr("uniqueness")?.extract()?,
            correct: ob.getattr("correct")?.extract()?,
            secondary_score: ob.getattr("secondary_score")?.extract()?,
            fragment_score: ob.getattr("fragment_score")?.extract()?,
            mate_mapped_to_disjoint_subgraph: ob
                .getattr("mate_mapped_to_disjoint_subgraph")?
                .extract()?,
            fragment_length_distribution: ob.getattr("fragment_length_distribution")?.extract()?,
            time_used: ob.getattr("time_used")?.extract()?,
            to_correct: ob.getattr("to_correct")?.extract()?,
            correctly_mapped: ob.getattr("correctly_mapped")?.extract()?,
            annotation,
        })
    }
}

impl FromPyObject<'_> for vg::Path {
    fn extract(ob: &PyAny) -> PyResult<Self> {
        Ok(Self {
            name: ob.getattr("name")?.extract()?,
            mapping: ob.getattr("mapping")?.extract()?,
            is_circular: ob.getattr("is_circular")?.extract()?,
            length: ob.getattr("length")?.extract()?,
        })
    }
}

impl FromPyObject<'_> for vg::Locus {
    fn extract(ob: &PyAny) -> PyResult<Self> {
        Ok(Self {
            name: ob.getattr("name")?.extract()?,
            allele: ob.getattr("allele")?.extract()?,
            support: ob.getattr("support")?.extract()?,
            genotype: ob.getattr("genotype")?.extract()?,
            overall_support: ob.getattr("overall_support")?.extract()?,
            allele_log_likelihood: ob.getattr("allele_log_likelihood")?.extract()?,
        })
    }
}

impl FromPyObject<'_> for vg::Position {
    fn extract(ob: &PyAny) -> PyResult<Self> {
        Ok(Self {
            node_id: ob.getattr("node_id")?.extract()?,
            offset: ob.getattr("offset")?.extract()?,
            is_reverse: ob.getattr("is_reverse")?.extract()?,
            name: ob.getattr("name")?.extract()?,
        })
    }
}

impl FromPyObject<'_> for vg::Mapping {
    fn extract(ob: &PyAny) -> PyResult<Self> {
        Ok(Self {
            position: ob.getattr("position")?.extract()?,
            edit: ob.getattr("edit")?.extract()?,
            rank: ob.getattr("rank")?.extract()?,
        })
    }
}

impl FromPyObject<'_> for vg::Support {
    fn extract(ob: &PyAny) -> PyResult<Self> {
        Ok(Self {
            quality: ob.getattr("quality")?.extract()?,
            forward: ob.getattr("forward")?.extract()?,
            reverse: ob.getattr("reverse")?.extract()?,
            left: ob.getattr("left")?.extract()?,
            right: ob.getattr("right")?.extract()?,
        })
    }
}

impl FromPyObject<'_> for vg::Genotype {
    fn extract(ob: &PyAny) -> PyResult<Self> {
        Ok(Self {
            allele: ob.getattr("allele")?.extract()?,
            is_phased: ob.getattr("is_phased")?.extract()?,
            likelihood: ob.getattr("likelihood")?.extract()?,
            log_likelihood: ob.getattr("log_likelihood")?.extract()?,
            log_prior: ob.getattr("log_prior")?.extract()?,
            log_posterior: ob.getattr("log_posterior")?.extract()?,
        })
    }
}

impl FromPyObject<'_> for vg::Edit {
    fn extract(ob: &PyAny) -> PyResult<Self> {
        Ok(Self {
            from_length: ob.getattr("from_length")?.extract()?,
            to_length: ob.getattr("to_length")?.extract()?,
            sequence: ob.getattr("sequence")?.extract()?,
        })
    }
}

#[pyfunction(name = "parse")]
fn parse(file_name: &str) -> PyResult<Vec<PyObject>> {
    let gam = gam::parse_from_file(file_name)?;
    Python::with_gil(|py| {
        let gam: Vec<_> = gam.iter().map(|o| o.clone().into_py(py)).collect();
        Ok(gam)
    })
}

#[pyfunction(name = "write")]
fn write(gams: Vec<PyObject>, file_name: &str) -> PyResult<()> {
    Python::with_gil(|py| -> PyResult<_> {
        let records = gams
            .iter()
            .map(|o| -> PyResult<_> { o.extract::<vg::Alignment>(py) })
            .collect::<PyResult<Vec<_>>>()?;
        gam::write_to_file(&records, file_name)?;
        Ok(())
    })
}

pub(crate) fn submodule(py: Python<'_>) -> PyResult<&PyModule> {
    let module = PyModule::new(py, "gam")?;
    module.add_function(wrap_pyfunction!(parse, module)?)?;
    module.add_function(wrap_pyfunction!(write, module)?)?;
    Ok(module)
}
