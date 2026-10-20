/// Sorted set of real parameters used to track boundary U/V values. Mirrors
/// OCCT's `IMeshData::IMapOfReal`; `f64` does not implement `Ord` (NaN), so the
/// values are kept in a sorted `Vec` ordered by `partial_cmp`.
#[derive(Debug, Clone, Default)]
pub struct ParamSet {
    pub(in crate::meshing::range_splitter) values: Vec<f64>,
}

impl ParamSet {
    /// Creates an empty parameter set.
    pub fn new() -> Self {
        Self { values: Vec::new() }
    }

    /// Inserts `v` keeping the set sorted; duplicates are ignored.
    pub fn insert(&mut self, v: f64) {
        if let Err(pos) = self
            .values
            .binary_search_by(|x| x.partial_cmp(&v).unwrap_or(std::cmp::Ordering::Equal))
        {
            self.values.insert(pos, v);
        }
    }

    /// True when `v` is present.
    pub fn contains(&self, v: &f64) -> bool {
        self.values
            .binary_search_by(|x| x.partial_cmp(v).unwrap_or(std::cmp::Ordering::Equal))
            .is_ok()
    }

    /// Number of stored values.
    pub fn len(&self) -> usize {
        self.values.len()
    }

    /// True when the set is empty.
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    /// Clears the set.
    pub fn clear(&mut self) {
        self.values.clear();
    }

    /// Iterator over the sorted values.
    pub fn iter(&self) -> std::slice::Iter<'_, f64> {
        self.values.iter()
    }

    /// The sorted values as a slice.
    pub fn values(&self) -> &[f64] {
        &self.values
    }
}
