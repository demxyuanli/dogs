use std::collections::HashMap;
use std::fmt;
use std::ops::Index;

/// Keyed slots with the traversal order of an OCCT list.
///
/// `BRep_TEdge::myCurves` is an `NCollection_List` (`BRep_TEdge.hxx:81`) and
/// `BRep_Builder::UpdateCurves` (`BRep_Builder.cxx:104-167`) removes the
/// representation it is about to replace and appends the new one at the end.
#[derive(Clone)]
pub struct OrderedSlots<V> {
    order: Vec<usize>,
    map: HashMap<usize, V>,
}

impl<V> Default for OrderedSlots<V> {
    fn default() -> Self {
        Self {
            order: Vec::new(),
            map: HashMap::new(),
        }
    }
}

impl<V> OrderedSlots<V> {
    /// Replace a slot and append it to the traversal order.
    pub fn set(&mut self, key: usize, value: V) -> Option<V> {
        self.order.retain(|k| *k != key);
        self.order.push(key);
        self.map.insert(key, value)
    }

    /// Insert a slot while preserving its existing traversal position.
    pub fn insert(&mut self, key: usize, value: V) -> Option<V> {
        if !self.map.contains_key(&key) {
            self.order.push(key);
        }
        self.map.insert(key, value)
    }

    pub fn get(&self, key: &usize) -> Option<&V> {
        self.map.get(key)
    }

    pub fn get_mut(&mut self, key: &usize) -> Option<&mut V> {
        self.map.get_mut(key)
    }

    pub fn contains_key(&self, key: &usize) -> bool {
        self.map.contains_key(key)
    }

    pub fn remove(&mut self, key: &usize) -> Option<V> {
        let out = self.map.remove(key);
        if out.is_some() {
            self.order.retain(|k| k != key);
        }
        out
    }

    /// Re-key in place while preserving the traversal position.
    pub fn rename_key(&mut self, from: usize, to: usize) {
        if from == to {
            return;
        }
        let Some(value) = self.map.remove(&from) else {
            return;
        };
        match self.order.iter().position(|k| *k == from) {
            Some(pos) => self.order[pos] = to,
            None => self.order.push(to),
        }
        self.map.insert(to, value);
        let Self { order, map } = self;
        let mut seen = std::collections::HashSet::new();
        order.retain(|k| map.contains_key(k) && seen.insert(*k));
    }

    pub fn iter(&self) -> impl Iterator<Item = (&usize, &V)> {
        self.order
            .iter()
            .filter_map(move |k| self.map.get(k).map(|v| (k, v)))
    }

    pub fn keys(&self) -> impl Iterator<Item = &usize> {
        self.order.iter().filter(|k| self.map.contains_key(k))
    }

    /// Return the same order with remapped keys.
    pub fn map_keys(&self, f: impl Fn(usize) -> usize) -> Self
    where
        V: Clone,
    {
        let mut out = Self::default();
        for (k, v) in self.iter() {
            out.insert(f(*k), v.clone());
        }
        out
    }

    pub fn for_each_mut(&mut self, mut f: impl FnMut(&usize, &mut V)) {
        let Self { order, map } = self;
        for k in order.iter() {
            if let Some(v) = map.get_mut(k) {
                f(k, v);
            }
        }
    }

    pub fn len(&self) -> usize {
        self.map.len()
    }

    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }
}

impl<V> Index<&usize> for OrderedSlots<V> {
    type Output = V;

    fn index(&self, key: &usize) -> &V {
        self.map.get(key).expect("ordered slot key")
    }
}

impl<'a, V> IntoIterator for &'a OrderedSlots<V> {
    type Item = (&'a usize, &'a V);
    type IntoIter = Box<dyn Iterator<Item = (&'a usize, &'a V)> + 'a>;

    fn into_iter(self) -> Self::IntoIter {
        Box::new(self.iter())
    }
}

impl<V: fmt::Debug> fmt::Debug for OrderedSlots<V> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_map().entries(self.iter()).finish()
    }
}
