//! `BRepTools_Modifier` (`BRepTools_Modifier.hxx:53-176`,
//! `BRepTools_Modifier.cxx:1-852`): replays a [`Modification`] over a whole
//! shape tree.
//!
//! `perform` is the reference's four preparation passes followed by the single
//! recursive rebuild:
//! 1. `CreateNewVertices` (`cxx:651-678`) — bind every vertex the modification
//!    gives a new point to.
//! 2. `FillNewCurveInfo` (`cxx:680-703`) — collect the `NewCurve` answers per
//!    edge.
//! 3. `FillNewSurfaceInfo` (`cxx:705-762`) — collect the `NewSurface` answers
//!    per face and record the faces whose sub-shapes are *not* going to change
//!    (`myNonUpdFace`).
//! 4. `CreateOtherVertices` (`cxx:764-832`) — bind the remaining vertices,
//!    sharing the original vertex when nothing below it changes.
//! 5. `Rebuild` (`cxx:226-649`) — the recursive walk.
//!
//! UNPORTED: `Message_ProgressScope` / `Message_ProgressRange` — the port has
//! no progress reporting, so the two "processing was broken" early returns
//! (`cxx:132-136`, `cxx:364-368`) have no counterpart.

use std::collections::{HashMap, HashSet};

use occt_core::kernel::geomabs::Shape as GeomAbsShape;

use crate::abs::{Orientation, ShapeType};
use crate::shape::TopoShape;
use crate::topo_tools::{explorer_shapes, is_same_occurrence, shape_map_key, ShapeMapKey};

use super::modification::{Modification, NewCurveInfo, NewSurfaceInfo};

/// `BRepTools_Modifier` (`BRepTools_Modifier.hxx:53-176`).
///
/// The maps are the reference's `myMap` / `myNCInfo` / `myNSInfo` /
/// `myNonUpdFace` / `myHasNewGeom`, keyed by `TopTools_ShapeMapHasher`
/// (`TShape` + location, [`ShapeMapKey`]).
pub struct Modifier {
    /// `myMap` — `NCollection_DataMap<TopoDS_Shape, TopoDS_Shape>`. `None` is
    /// a bound-but-null entry, which is how `Put` marks "seen, not rebuilt
    /// yet" (`BRepTools_Modifier.cxx:211-222`).
    pub(super) my_map: HashMap<ShapeMapKey, Option<TopoShape>>,
    /// `myShape`.
    pub(super) my_shape: Option<TopoShape>,
    /// `myDone`.
    my_done: bool,
    /// `myNCInfo`.
    pub(super) my_nc_info: HashMap<ShapeMapKey, NewCurveInfo>,
    /// `myNSInfo`.
    pub(super) my_ns_info: HashMap<ShapeMapKey, NewSurfaceInfo>,
    /// `myNonUpdFace`.
    pub(super) my_non_upd_face: HashSet<ShapeMapKey>,
    /// `myHasNewGeom`.
    pub(super) my_has_new_geom: HashSet<ShapeMapKey>,
    /// `myMutableInput`.
    pub(super) my_mutable_input: bool,
}
impl Default for Modifier {
    fn default() -> Self {
        Self::new(false)
    }
}

impl Modifier {
    /// `BRepTools_Modifier(bool)` (`cxx:53-59`).
    pub fn new(mutable_input: bool) -> Self {
        Self {
            my_map: HashMap::new(),
            my_shape: None,
            my_done: false,
            my_nc_info: HashMap::new(),
            my_ns_info: HashMap::new(),
            my_non_upd_face: HashSet::new(),
            my_has_new_geom: HashSet::new(),
            my_mutable_input: mutable_input,
        }
    }

    /// `BRepTools_Modifier(const TopoDS_Shape&)` (`cxx:61-68`).
    pub fn on_shape(shape: &TopoShape) -> Self {
        let mut m = Self::new(false);
        m.init(shape);
        m
    }

    /// `BRepTools_Modifier(const TopoDS_Shape&, const handle<BRepTools_Modification>&)`
    /// (`cxx:71-80`).
    pub fn on_shape_with(shape: &TopoShape, modification: &dyn Modification) -> Result<Self, String> {
        let mut m = Self::new(false);
        m.init(shape);
        m.perform(modification)?;
        Ok(m)
    }

    /// `Init(S)` (`cxx:83-88`). `myMap` is not cleared, exactly as in the
    /// reference.
    pub fn init(&mut self, shape: &TopoShape) {
        self.my_shape = Some(shape.clone());
        self.my_done = false;
        self.put(shape);
    }

    /// `IsDone()` (`BRepTools_Modifier.lxx:34-37`).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// `IsMutableInput()` (`cxx:844-847`).
    pub fn is_mutable_input(&self) -> bool {
        self.my_mutable_input
    }

    /// `SetMutableInput(bool)` (`cxx:849-852`).
    pub fn set_mutable_input(&mut self, v: bool) {
        self.my_mutable_input = v;
    }

    /// `ModifiedShape(S)` (`BRepTools_Modifier.lxx:21-28`).
    ///
    /// `Err` covers OCCT's two failure modes: an unbound shape throws
    /// `Standard_NoSuchObject` (`lxx:23-26`), and a bound-but-null entry is a
    /// shape that has not been rebuilt yet (only reachable before `perform`
    /// completes).
    pub fn modified_shape(&self, s: &TopoShape) -> Result<TopoShape, String> {
        match self.my_map.get(&shape_map_key(s)) {
            None => Err("BRepTools_Modifier::ModifiedShape: shape is not bound".into()),
            Some(None) => Err("BRepTools_Modifier::ModifiedShape: bound null shape".into()),
            Some(Some(shape)) => Ok(shape.clone()),
        }
    }

    /// `Put(S)` (`cxx:211-222`): bind `S` and, recursively, every shape of its
    /// tree with a null entry. The walk is `TopoDS_Iterator(S, false)`.
    pub fn put(&mut self, s: &TopoShape) {
        let key = shape_map_key(s);
        if self.my_map.contains_key(&key) {
            return;
        }
        self.my_map.insert(key, None);
        for child in children_keep_orientation(s) {
            self.put(&child);
        }
    }

    /// `Perform(M)` (`cxx:96-207`).
    pub fn perform(&mut self, m: &dyn Modification) -> Result<(), String> {
        let shape = self
            .my_shape
            .clone()
            .ok_or_else(|| "BRepTools_Modifier::Perform: null shape".to_string())?;

        // `cxx:114-116`.
        let mve = keyed_ancestors(&shape, ShapeType::Vertex, ShapeType::Edge);
        let mef = keyed_ancestors(&shape, ShapeType::Edge, ShapeType::Face);

        // `cxx:118-127`.
        self.create_new_vertices(&mve, m);
        self.fill_new_curve_info(&mef, m);
        self.fill_new_surface_info(m);
        if !self.my_mutable_input {
            self.create_other_vertices(&mve, &mef, m);
        }

        // `cxx:129-130`.
        let mut new_geom = false;
        self.rebuild(&shape, m, &mut new_geom)?;

        // `cxx:138-152`: the top shape's orientation. A reversed top face is
        // reversed in the map; every other case copies the orientation over.
        let key = shape_map_key(&shape);
        if let Some(Some(top)) = self.my_map.get_mut(&key) {
            if shape.shape_type() == ShapeType::Face && shape.orientation() == Orientation::Reversed {
                top.reverse();
            } else {
                let ori = shape.orientation();
                top.set_orientation(ori);
            }
        }

        // `cxx:166-199`: the continuity of every rebuilt edge between the two
        // faces that carry it. `aMEF.FindKey(ii)` walks the edge keys and the
        // ancestor list gives `F1` (first) and `F2` (second); the reference
        // stops looking after the second one (`cxx:175-186`) and skips edges
        // with `CurE.IsSame(NewE)` (`cxx:170`).
        //
        // UNPORTED: the reference's `TopoDS::Edge(myMap(CurE))` /
        // `TopoDS::Face(myMap(F1))` dereference the map entry with no null
        // test, so an unbound shape is undefined behaviour there
        // (`TopoDS.hxx:110-113`, `TopoDS_Shape.hxx:140`). Every edge and face
        // of `myShape` is bound by `put` + `rebuild`, so the port's `continue`
        // is unreachable in practice.
        for (cur_e, faces) in mef.iter() {
            let Some(new_e) = self.mapped(cur_e) else {
                continue;
            };
            if is_same_occurrence(cur_e, &new_e) {
                continue;
            }
            let mut f1: Option<&TopoShape> = None;
            let mut f2: Option<&TopoShape> = None;
            for f in faces {
                if f1.is_none() {
                    f1 = Some(f);
                } else {
                    f2 = Some(f);
                    break;
                }
            }
            let (Some(f1), Some(f2)) = (f1, f2) else {
                continue;
            };
            let new_f1 = match self.mapped_or_self(f1) {
                Some(s) => s,
                None => continue,
            };
            let new_f2 = match self.mapped_or_self(f2) {
                Some(s) => s,
                None => continue,
            };
            let new_cont = m.continuity(cur_e, f1, f2, &new_e, &new_f1, &new_f2);
            if new_cont > GeomAbsShape::C0 {
                // `aBB.Continuity(NewE, newf1, newf2, Newcont)`
                // (`BRep_Builder.cxx:1012-1021`).
                if let (Some(e), Some(nf1), Some(nf2)) = (
                    crate::shape::Edge::wrap(new_e.clone()),
                    crate::shape::Face::wrap(new_f1),
                    crate::shape::Face::wrap(new_f2),
                ) {
                    crate::tshape_continuity::set_stored_continuity(&e, &nf1, &nf2, new_cont as u8);
                }
            }
        }

        self.my_done = true;
        Ok(())
    }

    /// `myMap(S)` read back: `None` when unbound or bound null.
    pub(super) fn mapped(&self, s: &TopoShape) -> Option<TopoShape> {
        self.my_map.get(&shape_map_key(s)).and_then(|e| e.clone())
    }

    /// `(myMap.IsBound(S) ? myMap(S) : S)`, with a null result replaced by `S`
    /// (`BRepTools_Modifier.cxx:440-444`, `cxx:455-459`).
    pub(super) fn mapped_or_self(&self, s: &TopoShape) -> Option<TopoShape> {
        Some(self.mapped(s).unwrap_or_else(|| s.clone()))
    }
}

/// `TopoDS_Iterator(S, false)` (`TopoDS_Iterator.cxx:38-55`): the direct
/// children with the parent's location moved into them and their own
/// orientations kept (`theCumOri = false`, `theCumLoc` left at its `true`
/// default). Lives in [`crate::topo_tools`] because the `BRep_Builder` port
/// shares it.
pub(super) use crate::topo_tools::children_keep_orientation;

/// `TopExp::MapShapesAndAncestors(S, TS, TA, M)` (`TopExp.cxx:80-120`) as an
/// insertion-ordered list of `(sub-shape, ancestors)` pairs.
///
/// `topo_tools::map_shapes_and_ancestors` returns the same content in a
/// `HashMap`, which loses `NCollection_IndexedDataMap`'s insertion order; the
/// reference walks these maps by index (`FindKey(i)`), so the order is kept
/// here. The passes are the reference's: every `TA` ancestor contributes one
/// list entry per `TS` child (`cxx:96-105`), then every `TS` sub-shape without
/// a `TA` ancestor is added with an empty list (`cxx:110-119`).
pub(super) fn keyed_ancestors(
    shape: &TopoShape,
    ts: ShapeType,
    ta: ShapeType,
) -> Vec<(TopoShape, Vec<TopoShape>)> {
    let mut order: Vec<TopoShape> = Vec::new();
    let mut index: HashMap<ShapeMapKey, usize> = HashMap::new();
    let mut ancestors: Vec<Vec<TopoShape>> = Vec::new();
    for anc in explorer_shapes(shape, ta) {
        for sub in explorer_shapes(&anc, ts) {
            let key = shape_map_key(&sub);
            let i = match index.get(&key) {
                Some(&i) => i,
                None => {
                    let i = order.len();
                    index.insert(key, i);
                    order.push(sub);
                    ancestors.push(Vec::new());
                    i
                }
            };
            ancestors[i].push(anc.clone());
        }
    }
    for sub in crate::topo_tools::explorer_shapes_avoiding(shape, ts, ta) {
        let key = shape_map_key(&sub);
        if !index.contains_key(&key) {
            index.insert(key, order.len());
            order.push(sub);
            ancestors.push(Vec::new());
        }
    }
    order.into_iter().zip(ancestors).collect()
}
