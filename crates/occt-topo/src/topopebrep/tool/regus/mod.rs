//! `TopOpeBRepTool_REGUS` — regularize the shells of a solid.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepTool/TopOpeBRepTool_REGUS.cxx`
//! (OCCT 8.0.0, 1039 lines) and `TopOpeBRepTool_REGUS.hxx`.
//!
//! Layout of the port:
//! * `mod.rs` — the state (`hxx`), the constructor and `Init` (`cxx:88-118`),
//!   the accessors (`cxx:122-160`) and the face/edge map helpers.
//! * [`map_s`] — `MapS` (`cxx:164-227`).
//! * [`split_faces`] — `WireToFace` (`cxx:230-258`), `SplitF` (`cxx:262-397`)
//!   and `SplitFaces` (`cxx:400-471`).
//! * [`regu`] — `FUN_update` (`cxx:473-517`) and `REGU` (`cxx:520-666`).
//! * [`block`] — `InitBlock` (`cxx:670-706`) and `NextinBlock` (`cxx:709-851`).
//! * [`nearest_f`] — `FUN_vectors` (`cxx:853-887`) and `NearestF`
//!   (`cxx:891-1039`).
//!
//! OCCT's `NCollection_DataMap` / `NCollection_Map` walk in bucket order; the
//! port keeps insertion order (`*_order` mirrors) so a run is reproducible.

use std::collections::{HashMap, HashSet};

use crate::shape::TopoShape;
use crate::topo_tools::{shape_map_key, ShapeMapKey};

mod block;
mod map_s;
pub(crate) mod nearest_f;
mod regu;
mod split_faces;

/// `class TopOpeBRepTool_REGUS` (`hxx:35-88`).
pub struct Regus {
    /// `hasnewsplits` (`hxx:76`).
    has_new_splits: bool,
    /// `myFsplits` (`hxx:70-71`): ancestor face -> its splits.
    my_f_splits: HashMap<ShapeMapKey, Vec<TopoShape>>,
    /// `myOshNsh` (`hxx:72-73`): old shell -> its new shells.
    my_osh_nsh: HashMap<ShapeMapKey, Vec<TopoShape>>,
    /// `mynF` (`hxx:78`).
    my_n_f: i32,
    /// `myoldnF` (`hxx:79`).
    my_old_n_f: i32,
    /// `myS` (`hxx:80`); `None` is a null `TopoDS_Shape`.
    my_s: Option<TopoShape>,
    /// `mymapeFs` (`hxx:82-83`): edge -> its untouched ancestor faces.
    my_mape_fs: HashMap<ShapeMapKey, Vec<TopoShape>>,
    /// Insertion order of `my_mape_fs`.
    my_mape_fs_order: Vec<ShapeMapKey>,
    /// `mymapeFsstatic` (`hxx:84-85`): edge -> all its ancestor faces in `myS`.
    my_mape_fs_static: HashMap<ShapeMapKey, Vec<TopoShape>>,
    /// `mymapemult` (`hxx:86`): edges with more than two ancestor faces.
    my_mape_mult: HashSet<ShapeMapKey>,
    /// `myedstoconnect` (`hxx:87`): the block boundary still to connect.
    my_eds_to_connect: EdgeSet,
    /// `mylFinBlock` (`hxx:88`): the faces of the block being built.
    my_l_fin_block: Vec<TopoShape>,
    /// `myf` (`hxx:81`): the face being examined.
    my_f: Option<TopoShape>,
}

impl Regus {
    /// `TopOpeBRepTool_REGUS()` (`cxx:88-102`).
    pub fn new() -> Self {
        Self {
            has_new_splits: false,
            my_f_splits: HashMap::new(),
            my_osh_nsh: HashMap::new(),
            my_n_f: 0,
            my_old_n_f: 0,
            my_s: None,
            my_mape_fs: HashMap::new(),
            my_mape_fs_order: Vec::new(),
            my_mape_fs_static: HashMap::new(),
            my_mape_mult: HashSet::new(),
            my_eds_to_connect: EdgeSet::new(),
            my_l_fin_block: Vec::new(),
            my_f: None,
        }
    }

    /// `Init(S)` (`cxx:106-118`): `myFsplits` / `myOshNsh` are *not* cleared,
    /// they accumulate over the shells of one solid (`RegularizeS.cxx:55-72`).
    pub fn init(&mut self, s: &TopoShape) {
        self.has_new_splits = false;
        self.my_n_f = 0;
        self.my_old_n_f = 0;
        self.my_s = Some(s.clone());
        self.my_mape_fs.clear();
        self.my_mape_fs_order.clear();
        self.my_mape_fs_static.clear();
        self.my_mape_mult.clear();
        self.my_eds_to_connect.clear();
    }

    /// `S()` (`cxx:122-125`).
    pub fn s(&self) -> Option<&TopoShape> {
        self.my_s.as_ref()
    }

    /// `SetFsplits(Fsplits)` (`cxx:128-135`).
    pub fn set_f_splits(&mut self, f_splits: &HashMap<ShapeMapKey, Vec<TopoShape>>) {
        self.my_f_splits = f_splits.clone();
    }

    /// `GetFsplits(Fsplits)` (`cxx:137-143`).
    pub fn f_splits(&self) -> &HashMap<ShapeMapKey, Vec<TopoShape>> {
        &self.my_f_splits
    }

    /// `SetOshNsh(OshNsh)` (`cxx:146-152`).
    pub fn set_osh_nsh(&mut self, osh_nsh: &HashMap<ShapeMapKey, Vec<TopoShape>>) {
        self.my_osh_nsh = osh_nsh.clone();
    }

    /// `GetOshNsh(OshNsh)` (`cxx:154-160`).
    pub fn osh_nsh(&self) -> &HashMap<ShapeMapKey, Vec<TopoShape>> {
        &self.my_osh_nsh
    }

    /// `Hasnewsplits()` (`hxx:52`).
    pub fn has_new_splits(&self) -> bool {
        self.has_new_splits
    }

    /// `mymapeFs.Bind` / `mymapeFs.ChangeFind(e).Append(f)` (`cxx:186-197`).
    fn bind_face_to_edge(&mut self, e: &TopoShape, f: &TopoShape) {
        let key = shape_map_key(e);
        match self.my_mape_fs.get_mut(&key) {
            Some(lof) => lof.push(f.clone()),
            None => {
                self.my_mape_fs.insert(key.clone(), vec![f.clone()]);
                self.my_mape_fs_order.push(key);
            }
        }
    }

    /// `mymapeFs.IsBound(e)`.
    fn is_bound(&self, e: &TopoShape) -> bool {
        self.my_mape_fs.contains_key(&shape_map_key(e))
    }

    /// `mymapeFs.Find(e)`.
    fn find(&self, e: &TopoShape) -> Option<&Vec<TopoShape>> {
        self.my_mape_fs.get(&shape_map_key(e))
    }

    /// `mymapeFs.ChangeFind(e)`.
    fn change_find(&mut self, e: &TopoShape) -> Option<&mut Vec<TopoShape>> {
        self.my_mape_fs.get_mut(&shape_map_key(e))
    }

    /// `mymapeFs.UnBind(e)`.
    fn unbind(&mut self, e: &TopoShape) {
        self.my_mape_fs.remove(&shape_map_key(e));
    }

    /// `mymapeFs` keys, in insertion order (`NCollection_DataMap::Iterator`).
    fn edge_keys(&self) -> Vec<ShapeMapKey> {
        self.my_mape_fs_order
            .iter()
            .filter(|k| self.my_mape_fs.contains_key(*k))
            .cloned()
            .collect()
    }

    /// `mymapeFsstatic.Find(e)` (`cxx:769`).
    fn find_static(&self, e: &TopoShape) -> Option<&Vec<TopoShape>> {
        self.my_mape_fs_static.get(&shape_map_key(e))
    }
}

impl Default for Regus {
    fn default() -> Self {
        Self::new()
    }
}

/// `NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>` that keeps the
/// insertion order of its keys, so `NCollection_Map::Iterator` walks in a
/// reproducible order.
pub(crate) struct EdgeSet {
    order: Vec<ShapeMapKey>,
    items: HashMap<ShapeMapKey, TopoShape>,
}

impl EdgeSet {
    fn new() -> Self {
        Self {
            order: Vec::new(),
            items: HashMap::new(),
        }
    }

    fn clear(&mut self) {
        self.order.clear();
        self.items.clear();
    }

    fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    fn len(&self) -> usize {
        self.items.len()
    }

    fn contains(&self, e: &TopoShape) -> bool {
        self.items.contains_key(&shape_map_key(e))
    }

    fn add(&mut self, e: &TopoShape) {
        let key = shape_map_key(e);
        if self.items.insert(key.clone(), e.clone()).is_none() {
            self.order.push(key);
        }
    }

    fn remove(&mut self, e: &TopoShape) {
        self.items.remove(&shape_map_key(e));
    }

    /// Keys in insertion order, skipping the ones removed meanwhile.
    fn items(&self) -> Vec<TopoShape> {
        self.order
            .iter()
            .filter_map(|k| self.items.get(k).cloned())
            .collect()
    }
}
