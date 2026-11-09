//! `BRepAlgo_AsDes` (`ModelingAlgorithms/TKBool/BRepAlgo/BRepAlgo_AsDes.cxx`,
//! OCCT 8.0.0): the two-way shape/descendant table the offset and boolean
//! algorithms use to walk "what came from what".
//!
//! Two `NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>,
//! TopTools_ShapeMapHasher>` maps are stored here as `HashMap<ShapeMapKey,
//! Vec<TopoShape>>`: the hasher's equality is `TopoDS_Shape::IsSame` (TShape +
//! location, orientation ignored), which is exactly [`ShapeMapKey`], and every
//! list operation used below (`Append`, `InsertBefore`, ordered iteration) is
//! order-preserving, so a `Vec` reproduces `NCollection_List`.
//!
//! Port channels:
//! * [`AsDes::change_descendant`] answers `Option<&mut Vec<TopoShape>>`. OCCT
//!   returns a function-local `static NCollection_List` for an unbound key, so
//!   a caller's `ChangeDescendant(S).Append(X)` on an unbound `S` silently
//!   drops `X`; `None` is that dead end.
//! * [`AsDes::remove`] answers `Result<(), String>` for the two
//!   `Standard_ConstructionError` raises of `cxx:286-292`.

use std::collections::{HashMap, HashSet};

use crate::shape::TopoShape;
use crate::topo_tools::{is_same_occurrence, shape_map_key, ShapeMapKey};

/// `BRepAlgo_AsDes`.
#[derive(Debug, Clone, Default)]
pub struct AsDes {
    /// `up`: `SS -> list of S` (the ascendants of `SS`).
    up: HashMap<ShapeMapKey, Vec<TopoShape>>,
    /// `down`: `S -> list of SS` (the descendants of `S`).
    down: HashMap<ShapeMapKey, Vec<TopoShape>>,
}

impl AsDes {
    /// `BRepAlgo_AsDes::BRepAlgo_AsDes()` (`cxx:28`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `BRepAlgo_AsDes::Clear()` (`cxx:62-66`).
    pub fn clear(&mut self) {
        self.up.clear();
        self.down.clear();
    }

    /// `BRepAlgo_AsDes::Add(S, SS)` (`cxx:32-47`): `SS` becomes a descendant of
    /// `S`. Both bind their key with an empty list when unbound, then append, so
    /// the lists keep call order and may hold duplicates (`IsSame` repeats are
    /// not filtered here, only by [`AsDes::replace`]).
    pub fn add(&mut self, s: &TopoShape, ss: &TopoShape) {
        self.down
            .entry(shape_map_key(s))
            .or_default()
            .push(ss.clone());
        self.up
            .entry(shape_map_key(ss))
            .or_default()
            .push(s.clone());
    }

    /// `BRepAlgo_AsDes::Add(S, SS)` with an `NCollection_List` (`cxx:51-58`):
    /// the per-element overload, in list order.
    pub fn add_list(&mut self, s: &TopoShape, ss: &[TopoShape]) {
        for shape in ss {
            self.add(s, shape);
        }
    }

    /// `BRepAlgo_AsDes::HasAscendant(S)` (`cxx:70-73`): the key is *bound* in
    /// `up`, even when its list is empty.
    pub fn has_ascendant(&self, s: &TopoShape) -> bool {
        self.up.contains_key(&shape_map_key(s))
    }

    /// `BRepAlgo_AsDes::HasDescendant(S)` (`cxx:77-80`).
    pub fn has_descendant(&self, s: &TopoShape) -> bool {
        self.down.contains_key(&shape_map_key(s))
    }

    /// `BRepAlgo_AsDes::Ascendant(S)` (`cxx:85-91`): an unbound key answers the
    /// static empty list, i.e. the empty slice.
    pub fn ascendant(&self, s: &TopoShape) -> &[TopoShape] {
        self.up
            .get(&shape_map_key(s))
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// `BRepAlgo_AsDes::Descendant(S)` (`cxx:96-102`).
    pub fn descendant(&self, s: &TopoShape) -> &[TopoShape] {
        self.down
            .get(&shape_map_key(s))
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// `BRepAlgo_AsDes::ChangeDescendant(S)` (`cxx:107-114`). `None` is OCCT's
    /// static empty list for an unbound key, where a caller's append is lost.
    pub fn change_descendant(&mut self, s: &TopoShape) -> Option<&mut Vec<TopoShape>> {
        self.down.get_mut(&shape_map_key(s))
    }

    /// `BRepAlgo_AsDes::HasCommonDescendant(S1, S2, LC)` (`cxx:166-188`): clears
    /// `lc`, then for every `DS1` of `S1` appends it once per `S2` occurrence
    /// among `DS1`'s ascendants. Returns `!lc.is_empty()`.
    ///
    /// `LC` is an out-parameter only, so no `RefCell` is needed.
    pub fn has_common_descendant(
        &self,
        s1: &TopoShape,
        s2: &TopoShape,
        lc: &mut Vec<TopoShape>,
    ) -> bool {
        lc.clear();
        if self.has_descendant(s1) && self.has_descendant(s2) {
            for ds1 in self.descendant(s1) {
                for ads1 in self.ascendant(ds1) {
                    // `ADS1.IsSame(S2)` (`cxx:181`).
                    if is_same_occurrence(ads1, s2) {
                        lc.push(ds1.clone());
                    }
                }
            }
        }
        !lc.is_empty()
    }

    /// `BRepAlgo_AsDes::Replace(OldS, NewS)` (`cxx:222-278`): the two passes of
    /// `cxx:225` (`up` first, then `down`).
    ///
    /// Per pass: back-replace inside every list that mentions `OldS`
    /// (`cxx:232`), then either move `OldS`'s filtered list onto `NewS`
    /// (`cxx:234-253`) or append its still-missing elements to `NewS`'s list
    /// (`cxx:255-274`), and finally unbind `OldS` (`cxx:276`).
    pub fn replace(&mut self, old_s: &TopoShape, new_s: &TopoShape) {
        for pass in 0..2 {
            // `aMap = !i ? up : down` (`cxx:226-228`).
            let old_key = shape_map_key(old_s);
            let new_key = shape_map_key(new_s);
            // `ChangeSeek(OldS)`; `nullptr` continues (`cxx:229-232`).
            let Some(old_list) = (if pass == 0 {
                self.up.get(&old_key)
            } else {
                self.down.get(&old_key)
            })
            .cloned() else {
                continue;
            };
            // `BackReplace(OldS, NewS, *pLSOld, InUp)` with `InUp = i != 0`
            // (`cxx:231`).
            let in_up = pass != 0;
            for s in &old_list {
                let key = shape_map_key(s);
                let target = if in_up {
                    self.up.get_mut(&key)
                } else {
                    self.down.get_mut(&key)
                };
                if let Some(list) = target {
                    replace_in_list(old_s, new_s, list);
                }
            }

            let map = if pass == 0 {
                &mut self.up
            } else {
                &mut self.down
            };
            if !map.contains_key(&new_key) {
                // `cxx:238-253`: filter the list by `IsSame`, then bind it to
                // `NewS`.
                let mut list = old_list;
                dedupe_in_place(&mut list);
                map.insert(new_key, list);
            } else if let Some(new_list) = map.get_mut(&new_key) {
                // `cxx:255-274`: append only the elements `NewS`'s list does
                // not hold yet.
                let mut seen: HashSet<ShapeMapKey> = new_list.iter().map(shape_map_key).collect();
                for a_s in &old_list {
                    if seen.insert(shape_map_key(a_s)) {
                        new_list.push(a_s.clone());
                    }
                }
            }
            // `aMap.UnBind(OldS)` (`cxx:276`).
            map.remove(&old_key);
        }
    }

    /// `BRepAlgo_AsDes::Remove(SS)` (`cxx:282-299`): only a shape with no
    /// descendants and at least one ascendant may be removed; `SS` is dropped
    /// from every ascendant's descendant list and unbound from `up`.
    ///
    /// `Err` is the two `Standard_ConstructionError` raises (`cxx:286-293`).
    pub fn remove(&mut self, ss: &TopoShape) -> Result<(), String> {
        let key = shape_map_key(ss);
        if self.down.contains_key(&key) {
            return Err("Standard_ConstructionError: BRepAlgo_AsDes::Remove".to_string());
        }
        if !self.up.contains_key(&key) {
            return Err("Standard_ConstructionError: BRepAlgo_AsDes::Remove".to_string());
        }
        for asc in self.up.get(&key).cloned().unwrap_or_default() {
            if let Some(list) = self.down.get_mut(&shape_map_key(&asc)) {
                remove_in_list(ss, list);
            }
        }
        self.up.remove(&key);
        Ok(())
    }
}

/// `ReplaceInList(OldS, NewS, L)` (`cxx:119-146`): `aMS` is a
/// `NCollection_Map<TopoDS_Shape>`, i.e. `IsSame` keys ignoring orientation;
/// every `IsSame(OldS)` entry is replaced in place by `NewS.Oriented(O)` of the
/// same orientation when that occurrence is not already present, and dropped
/// otherwise. `InsertBefore` + `Remove` at the same iterator is an in-place
/// substitution, so the port writes the new element at the index and drops the
/// old slot.
fn replace_in_list(old_s: &TopoShape, new_s: &TopoShape, list: &mut Vec<TopoShape>) {
    let mut seen: HashSet<ShapeMapKey> = list.iter().map(shape_map_key).collect();
    let mut i = 0;
    while i < list.len() {
        if is_same_occurrence(&list[i], old_s) {
            let oriented = new_s.oriented(list[i].orientation());
            if seen.insert(shape_map_key(&oriented)) {
                // `L.InsertBefore(NewS.Oriented(O), it)` then `L.Remove(it)`.
                list[i] = oriented;
                i += 1;
            } else {
                // `L.Remove(it)` without an insert.
                list.remove(i);
            }
        } else {
            i += 1;
        }
    }
}

/// `RemoveInList(S, L)` (`cxx:150-163`): drop the first `IsSame(S)` entry only.
fn remove_in_list(s: &TopoShape, list: &mut Vec<TopoShape>) {
    if let Some(i) = list.iter().position(|v| is_same_occurrence(v, s)) {
        list.remove(i);
    }
}

/// The `!aMS.Add(aIt.Value())` filter of `cxx:240-252`: keep the first
/// occurrence of every `IsSame` key.
fn dedupe_in_place(list: &mut Vec<TopoShape>) {
    let mut seen: HashSet<ShapeMapKey> = HashSet::new();
    list.retain(|v| seen.insert(shape_map_key(v)));
}
