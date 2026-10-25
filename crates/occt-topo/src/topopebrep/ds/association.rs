//! `TopOpeBRepDS_Association` — a symmetric "these interferences belong
//! together" relation, keyed on interference identity.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_Association.hxx:30-70`,
//! `TopOpeBRepDS_Association.cxx:22-108`.
//!
//! OCCT keys `myMap` on the `handle<TopOpeBRepDS_Interference>`, i.e. on the
//! pointer, and OCCT's hasher for handles compares that pointer. The port
//! therefore stores `(handle, list)` pairs and compares with
//! [`std::sync::Arc::ptr_eq`], which is the same identity relation without
//! requiring a `Hash` impl that would silently key on the interference value.
//!
//! The membership test inside a list is `HasSameGeometry` (`cxx:26-36`), not
//! identity.

use super::interference::InterferenceHandle;
use super::interference_iterator::InterferenceList;

/// `static bool Contains(LI, I)` (`cxx:26-36`).
fn contains(li: &InterferenceList, i: &InterferenceHandle) -> bool {
    let ih = i.read().expect("interference lock");
    li.iter().any(|other| {
        let oh = other.read().expect("interference lock");
        ih.has_same_geometry(&oh)
    })
}

/// `TopOpeBRepDS_Association` (`hxx:30-70`).
#[derive(Clone, Debug, Default)]
pub struct Association {
    /// `myMap`: `DataMap<handle, List>` with pointer identity.
    map: Vec<(InterferenceHandle, InterferenceList)>,
}

impl Association {
    /// `TopOpeBRepDS_Association()` (`cxx:42`): `= default`.
    pub fn new() -> Self {
        Self::default()
    }

    /// Position of `i` in `myMap`, by handle identity.
    fn position(&self, i: &InterferenceHandle) -> Option<usize> {
        self.map
            .iter()
            .position(|(k, _)| std::sync::Arc::ptr_eq(k, i))
    }

    /// `Associate(I, K)` (`cxx:46-67`): both directions, each appended only
    /// when the other side is not already bound to it by geometry.
    pub fn associate(&mut self, i: &InterferenceHandle, k: &InterferenceHandle) {
        match self.position(i) {
            None => self.map.push((i.clone(), vec![k.clone()])),
            Some(p) => {
                if !contains(&self.map[p].1, k) {
                    self.map[p].1.push(k.clone());
                }
            }
        }
        match self.position(k) {
            None => self.map.push((k.clone(), vec![i.clone()])),
            Some(p) => {
                if !contains(&self.map[p].1, i) {
                    self.map[p].1.push(i.clone());
                }
            }
        }
    }

    /// `Associate(I, LI)` (`cxx:72-81`).
    pub fn associate_list(&mut self, i: &InterferenceHandle, li: &InterferenceList) {
        for k in li {
            self.associate(i, k);
        }
    }

    /// `HasAssociation(I)` (`cxx:85-88`).
    pub fn has_association(&self, i: &InterferenceHandle) -> bool {
        self.position(i).is_some()
    }

    /// `Associated(I)` (`cxx:92-101`).
    ///
    /// OCCT returns a `static` empty list for an unbound handle; `None` is that
    /// same empty result, and a write through it is discarded there too.
    pub fn associated(&mut self, i: &InterferenceHandle) -> Option<&mut InterferenceList> {
        match self.position(i) {
            Some(p) => Some(&mut self.map[p].1),
            None => None,
        }
    }

    /// `AreAssociated(I, K)` (`cxx:105-110`).
    pub fn are_associated(&self, i: &InterferenceHandle, k: &InterferenceHandle) -> bool {
        match self.position(i) {
            Some(p) => contains(&self.map[p].1, k),
            None => false,
        }
    }
}
