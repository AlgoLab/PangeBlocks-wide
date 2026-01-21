use prost_types::{value::Kind, ListValue, Struct, Value};
use pyo3::{
    prelude::*,
    types::{PyDict, PyList},
};

pub(crate) fn listvalue_to_pylist<'a>(py: Python<'a>, v: &ListValue) -> PyResult<&'a PyList> {
    let list = PyList::empty(py);
    for item in v.values.iter() {
        let value = value_to_pyany(py, item)?;
        list.append(value)?;
    }
    Ok(list)
}

pub(crate) fn struct_to_pydict<'a>(py: Python<'a>, s: &Struct) -> PyResult<&'a PyDict> {
    let dict = PyDict::new(py);
    for (key, value) in s.fields.iter() {
        let value = value_to_pyany(py, value)?;
        dict.set_item(key, value)?;
    }
    Ok(dict)
}

pub(crate) fn pydict_to_struct(dict: &PyDict) -> PyResult<Struct> {
    let mut s = Struct::default();
    for (key, value) in dict.iter() {
        let key = key.extract::<String>()?;
        let value = pyany_to_value(value)?;
        s.fields.insert(key, value);
    }
    Ok(s)
}

pub(crate) fn pylist_to_listvalue(list: &PyList) -> PyResult<ListValue> {
    let mut l = ListValue::default();
    for item in list.iter() {
        let value = pyany_to_value(item);
        l.values.push(value?);
    }
    Ok(l)
}

pub(crate) fn pyany_to_value(object: &PyAny) -> PyResult<Value> {
    let kind = if let Ok(object) = object.extract::<bool>() {
        Kind::BoolValue(object)
    } else if let Ok(object) = object.extract::<f64>() {
        Kind::NumberValue(object)
    } else if let Ok(object) = object.extract::<String>() {
        Kind::StringValue(object)
    } else if let Ok(object) = object.downcast::<PyDict>() {
        Kind::StructValue(pydict_to_struct(object)?)
    } else if let Ok(object) = object.downcast::<PyList>() {
        Kind::ListValue(pylist_to_listvalue(object)?)
    } else {
        Kind::NullValue(0)
    };

    Ok(Value { kind: Some(kind) })
}

pub(crate) fn value_to_pyany(py: Python, value: &Value) -> PyResult<Py<PyAny>> {
    let value = match &value.kind {
        Some(Kind::NullValue(_)) => py.None(),
        Some(Kind::NumberValue(v)) => v.into_py(py),
        Some(Kind::StringValue(v)) => v.into_py(py),
        Some(Kind::BoolValue(v)) => v.into_py(py),
        Some(Kind::StructValue(ref v)) => struct_to_pydict(py, v)?.into_py(py),
        Some(Kind::ListValue(ref v)) => listvalue_to_pylist(py, v)?.into_py(py),
        None => py.None(),
    };
    Ok(value)
}
