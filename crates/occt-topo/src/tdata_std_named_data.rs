//! `TDataStd_NamedData` (TKLCAF/TDataStd) — the named attribute container the
//! glTF reader fills from a node's or a mesh's `extras` object
//! (`RWGltf_ExtrasParser`, `RWGltf_GltfJsonParser.cxx:183-319`) and attaches to
//! the document label with `setShapeNamedData` (`RWGltf_CafReader.cxx:658`).
//!
//! Only the members the extras parser writes are modelled: named integers,
//! reals and strings, plus arrays of integers and reals. OCCT keeps them in
//! `TDataStd_IntegerArray` / `TDataStd_RealArray` child attributes; the maps
//! below stand in for those attributes and are keyed by the attribute name the
//! parser builds (nested `extras` objects are flattened with `.` separators,
//! `cxx:194-197`).

use std::collections::BTreeMap;

/// `TDataStd_NamedData`, restricted to what `RWGltf_ExtrasParser` fills.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct TDataStdNamedData {
    integers: BTreeMap<String, i32>,
    reals: BTreeMap<String, f64>,
    strings: BTreeMap<String, String>,
    arrays_of_integers: BTreeMap<String, Vec<i32>>,
    arrays_of_reals: BTreeMap<String, Vec<f64>>,
}

impl TDataStdNamedData {
    pub fn new() -> Self {
        Self::default()
    }

    /// `TDataStd_NamedData::SetInteger`.
    pub fn set_integer(&mut self, name: &str, value: i32) {
        self.integers.insert(name.to_string(), value);
    }

    /// `TDataStd_NamedData::GetInteger`.
    pub fn integer(&self, name: &str) -> Option<i32> {
        self.integers.get(name).copied()
    }

    /// `TDataStd_NamedData::SetReal`.
    pub fn set_real(&mut self, name: &str, value: f64) {
        self.reals.insert(name.to_string(), value);
    }

    /// `TDataStd_NamedData::GetReal`.
    pub fn real(&self, name: &str) -> Option<f64> {
        self.reals.get(name).copied()
    }

    /// `TDataStd_NamedData::SetString`.
    pub fn set_string(&mut self, name: &str, value: &str) {
        self.strings.insert(name.to_string(), value.to_string());
    }

    /// `TDataStd_NamedData::GetString`.
    pub fn string(&self, name: &str) -> Option<&str> {
        self.strings.get(name).map(String::as_str)
    }

    /// `TDataStd_NamedData::SetArrayOfIntegers`.
    pub fn set_array_of_integers(&mut self, name: &str, values: Vec<i32>) {
        self.arrays_of_integers.insert(name.to_string(), values);
    }

    /// `TDataStd_NamedData::GetArrayOfIntegers`.
    pub fn array_of_integers(&self, name: &str) -> Option<&[i32]> {
        self.arrays_of_integers.get(name).map(Vec::as_slice)
    }

    /// `TDataStd_NamedData::SetArrayOfReals`.
    pub fn set_array_of_reals(&mut self, name: &str, values: Vec<f64>) {
        self.arrays_of_reals.insert(name.to_string(), values);
    }

    /// `TDataStd_NamedData::GetArrayOfReals`.
    pub fn array_of_reals(&self, name: &str) -> Option<&[f64]> {
        self.arrays_of_reals.get(name).map(Vec::as_slice)
    }

    pub fn integers(&self) -> impl Iterator<Item = (&str, i32)> {
        self.integers.iter().map(|(k, v)| (k.as_str(), *v))
    }

    pub fn reals(&self) -> impl Iterator<Item = (&str, f64)> {
        self.reals.iter().map(|(k, v)| (k.as_str(), *v))
    }

    pub fn strings(&self) -> impl Iterator<Item = (&str, &str)> {
        self.strings.iter().map(|(k, v)| (k.as_str(), v.as_str()))
    }

    pub fn arrays_of_integers(&self) -> impl Iterator<Item = (&str, &[i32])> {
        self.arrays_of_integers.iter().map(|(k, v)| (k.as_str(), v.as_slice()))
    }

    pub fn arrays_of_reals(&self) -> impl Iterator<Item = (&str, &[f64])> {
        self.arrays_of_reals.iter().map(|(k, v)| (k.as_str(), v.as_slice()))
    }

    pub fn is_empty(&self) -> bool {
        self.integers.is_empty()
            && self.reals.is_empty()
            && self.strings.is_empty()
            && self.arrays_of_integers.is_empty()
            && self.arrays_of_reals.is_empty()
    }
}
