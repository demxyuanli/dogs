//! `LocOpe_SplitShape` (`ModelingAlgorithms/TKFeat/LocOpe/LocOpe_SplitShape`).
//!
//! The class keeps one *descendant list* per shape of the original model
//! (`myMap`). `Put` fills the map with the original tree, every splitting
//! operation replaces a list by the pieces it produced, and `Rebuild` assembles
//! the answer bottom-up. Only `Rebuild`, `DescendantShapes` and `LeftOf` read
//! the map; the `Add` overloads underneath only rewrite it.
//!
//! Module layout (each mirrors a block of the reference):
//! * [`add_vertex`] - `Add(V, P, E)` (`cxx:143-249`).
//! * [`add_wires`] - `Add(Lwires, F)` (`cxx:256-598`) and `Add(W, F)` (`cxx:602-654`).
//! * [`add_closed_wire`] - `AddClosedWire` (`cxx:658-732`).
//! * [`add_open_wire`] - `AddOpenWire` (`cxx:781-1271`).
//! * [`is_inside`] - the two `IsInside` overloads (`cxx:1459-1593`).
//! * [`direction`] - `GetDirection` / `ChoixUV` / `ChooseDirection` /
//!   `SameUV` / `checkOverlapping` (`cxx:75-96`, `735-777`, `1595-1776`).

mod add_closed_wire;
mod add_open_wire;
mod add_vertex;
mod add_wires;
mod direction;
mod is_inside;

use crate::abs::{Orientation, ShapeType};
use crate::builder::TopoBuilder;
use crate::iterator::cumulated_children;
use crate::locope::shape_maps::{ShapeSet, ShapeValueMap};
use crate::shape::{Edge, Face, TopoShape, Wire};
use crate::topo_tools::is_same;

pub use direction::{
    check_overlapping, choose_direction, choix_uv, get_direction, same_uv,
};

/// `LocOpe_SplitShape` (`hxx:36-93`).
#[derive(Debug, Clone)]
pub struct SplitShape {
    /// `myDone` (`hxx:87`).
    my_done: bool,
    /// `myShape` (`hxx:88`). `None` is the reference's null `TopoDS_Shape`.
    my_shape: Option<TopoShape>,
    /// `myMap` (`hxx:89`):
    /// `NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>>`.
    my_map: ShapeValueMap<Vec<TopoShape>>,
    /// `myDblE` (`hxx:90`): the edges doubled to close a split face.
    my_dbl_e: ShapeSet,
    /// `myLeft` (`hxx:91`): `LeftOf`'s answer.
    my_left: Vec<TopoShape>,
}

impl Default for SplitShape {
    fn default() -> Self {
        Self::new()
    }
}

impl SplitShape {
    /// `LocOpe_SplitShape()` (`lxx:22-25`): `myDone = false`, no shape.
    pub fn new() -> Self {
        Self {
            my_done: false,
            my_shape: None,
            my_map: ShapeValueMap::new(),
            my_dbl_e: ShapeSet::new(),
            my_left: Vec::new(),
        }
    }

    /// `LocOpe_SplitShape(const TopoDS_Shape& S)` (`lxx:31-34`).
    pub fn with_shape(s: &TopoShape) -> Self {
        let mut this = Self::new();
        this.init(s);
        this
    }

    /// `Init(const TopoDS_Shape& S)` (`cxx:94-101`).
    pub fn init(&mut self, s: &TopoShape) {
        self.my_done = false;
        self.my_shape = Some(s.clone());
        self.my_dbl_e.clear();
        self.my_map.clear();
        let shape = s.clone();
        self.put(&shape);
    }

    /// `Shape()` (`lxx:37-40`).
    pub fn shape(&self) -> Option<&TopoShape> {
        self.my_shape.as_ref()
    }

    /// `CanSplit(const TopoDS_Edge& E)` (`cxx:105-139`).
    pub fn can_split(&self, e: &Edge) -> bool {
        // `cxx:107-111`.
        if self.my_done {
            return false;
        }
        // `cxx:112-116`.
        if self.my_map.is_empty() {
            return false;
        }
        // `cxx:117-120`.
        if !self.my_map.is_bound(&e.0) {
            return false;
        }

        // `cxx:122-137`: the edge must not belong to a wire that was already
        // rebuilt.
        for (key, value) in self.my_map.iter() {
            if key.shape_type() == ShapeType::Wire && !value.is_empty() {
                for sub in crate::topo_tools::explorer_shapes(key, ShapeType::Edge) {
                    if is_same(&sub, &e.0) {
                        return false;
                    }
                }
            }
        }
        true
    }

    /// `DescendantShapes(const TopoDS_Shape& S)` (`cxx:1337-1352`).
    pub fn descendant_shapes(&mut self, s: &TopoShape) -> Vec<TopoShape> {
        // `cxx:1339-1343`.
        if !self.my_done {
            let Some(shape) = self.my_shape.clone() else {
                // A null `myShape` was never `Put`, so `myMap(myShape)` would
                // raise `Standard_NoSuchObject`.
                self.my_done = true;
                return Vec::new();
            };
            self.rebuild(&shape);
            self.my_done = true;
        }
        // `cxx:1350`: the reference hands out its own list.
        self.my_map.find(s).cloned().unwrap_or_default()
    }

    /// `Put(const TopoDS_Shape& S)` (`cxx:1355-1373`): bind `S` and, unless it
    /// is a vertex, recurse into its children. A vertex seeds its list with
    /// itself (`cxx:1369-1372`).
    fn put(&mut self, s: &TopoShape) {
        if self.my_map.is_bound(s) {
            return;
        }
        self.my_map.bind(s, Vec::new());
        if s.shape_type() != ShapeType::Vertex {
            for child in cumulated_children(s) {
                self.put(&child);
            }
        } else if let Some(list) = self.my_map.change(s) {
            list.push(s.clone());
        }
    }

    /// `Rebuild(const TopoDS_Shape& S)` (`cxx:1411-1455`).
    ///
    /// Returns whether anything below `S` was rebuilt. `myMap(S)` raises
    /// `Standard_NoSuchObject` in OCCT when `S` was never `Put`; the port
    /// answers `false` for the top-level call from [`Self::descendant_shapes`]
    /// on a null shape, and otherwise follows the reference.
    fn rebuild(&mut self, s: &TopoShape) -> bool {
        // `cxx:1412-1415`.
        if s.shape_type() == ShapeType::Face {
            if let Some(f) = Face::wrap(s.clone()) {
                update_tolerances(&f, &self.my_map);
            }
        }

        // `cxx:1416-1422`.
        if let Some(list) = self.my_map.find(s) {
            if !list.is_empty() {
                return !is_same(&list[0], s);
            }
        }

        // `cxx:1424-1428`: `TopoDS_Iterator it(S)` walks the children with the
        // parent orientation composed in.
        let children = cumulated_children(s);
        let mut rebuild = false;
        for child in &children {
            rebuild = self.rebuild(child) || rebuild;
        }

        if rebuild {
            // `cxx:1431-1448`.
            let b = TopoBuilder::new();
            let mut result = s.empty_copied();
            for child in &children {
                let orient = child.orientation();
                let descendants = self.my_map.find(child).cloned().unwrap_or_default();
                for d in &descendants {
                    let mut shaped = d.clone();
                    shaped.set_orientation(Orientation::compose(orient, d.orientation()));
                    b.add(&mut result, &shaped);
                }
            }
            // `cxx:1442-1445`: `Closed` is assigned for wires and shells only,
            // from `BRep_Tool::IsClosed(result)`, i.e. the flag the empty copy
            // left at `false`.
            if result.shape_type() == ShapeType::Wire || result.shape_type() == ShapeType::Shell {
                result.set_closed(brep_tool_is_closed(&result));
            }
            if let Some(list) = self.my_map.change(s) {
                list.push(result);
            }
        } else {
            // `cxx:1450-1453`.
            if let Some(list) = self.my_map.change(s) {
                list.push(s.clone());
            }
        }
        rebuild
    }

    /// `LeftOf(const TopoDS_Wire& W, const TopoDS_Face& F)` (`cxx:1273-1335`).
    /// `Err` is the reference's `Standard_NoSuchObject` (`cxx:1276-1278`,
    /// `1290-1293`).
    pub fn left_of(&mut self, w: &Wire, f: &Face) -> Result<Vec<TopoShape>, String> {
        // `cxx:1276-1278`.
        let Some(my_shape) = self.my_shape.clone() else {
            return Err("Standard_NoSuchObject: LocOpe_SplitShape::LeftOf".into());
        };

        // `cxx:1281-1293`: the occurrence of `F` inside `myShape`.
        let Some(the_face) = crate::topo_tools::explorer_shapes(&my_shape, ShapeType::Face)
            .into_iter()
            .find(|s| is_same(s, &f.0))
        else {
            return Err("Standard_NoSuchObject: LocOpe_SplitShape::LeftOf".into());
        };
        self.my_left.clear();

        let or_face = the_face.orientation();
        let descendants = self.my_map.find(&the_face).cloned().unwrap_or_default();

        // `cxx:1300-1331`.
        for edg in crate::topo_tools::explorer_shapes(&w.0, ShapeType::Edge) {
            let Some(edg) = Edge::wrap(edg) else {
                continue;
            };
            for fac in &descendants {
                let Some(mut fac) = Face::wrap(fac.clone()) else {
                    continue;
                };
                fac.0.set_orientation(or_face);
                let mut found = false;
                for edgbis in crate::topo_tools::explorer_shapes(&fac.0, ShapeType::Edge) {
                    if is_same(&edgbis, &edg.0) && edgbis.orientation() == edg.0.orientation() {
                        // `cxx:1313-1323`: append the face once.
                        if !self.my_left.iter().any(|s| is_same(s, &fac.0)) {
                            self.my_left.push(fac.0.clone());
                        }
                        found = true;
                        break;
                    }
                }
                if found {
                    break;
                }
            }
        }
        Ok(self.my_left.clone())
    }

    /// `myMap(K)` read side: the bound list, or an empty one.
    pub(super) fn map_of(&self, key: &TopoShape) -> Vec<TopoShape> {
        self.my_map.find(key).cloned().unwrap_or_default()
    }

    /// `myMap.IsBound(K)`.
    pub(super) fn map_bound(&self, key: &TopoShape) -> bool {
        self.my_map.is_bound(key)
    }

    /// `myMap.Bind(K, V)`.
    pub(super) fn map_bind(&mut self, key: &TopoShape, value: Vec<TopoShape>) {
        self.my_map.bind(key, value);
    }

    /// `myMap(K) = V`, i.e. `ChangeFind` on a bound key.
    pub(super) fn set_map_of(&mut self, key: &TopoShape, value: Vec<TopoShape>) {
        if let Some(list) = self.my_map.change(key) {
            *list = value;
        }
    }

    /// `myMap.UnBind(K)`.
    pub(super) fn map_remove(&mut self, key: &TopoShape) {
        self.my_map.unbind(key);
    }

    /// `myDone`.
    pub(super) fn done(&self) -> bool {
        self.my_done
    }

    /// `myDone = ...`.
    #[allow(dead_code)]
    pub(super) fn set_done(&mut self, done: bool) {
        self.my_done = done;
    }

    /// `myDblE.Contains(E)`.
    #[allow(dead_code)]
    pub(super) fn is_double_edge(&self, e: &TopoShape) -> bool {
        self.my_dbl_e.contains(e)
    }

    /// `myDblE.Add(E)`.
    pub(super) fn add_double_edge(&mut self, e: &TopoShape) {
        self.my_dbl_e.add(e);
    }
}

/// `TopoDS_Shape::EmptyCopied()` on an edge, with the *geometry* `BRep_TEdge::EmptyCopy`
/// (`BRep_TEdge.cxx:105-133`) carries over: the curve, the range, the
/// tolerance, the degenerated / same-parameter / same-range flags and every
/// `CurveOnSurface` representation (p-curves included) - only the children are
/// dropped. This port keeps all of that in [`crate::tgeometry::EdgeGeom`], so
/// cloning the slot is exactly that empty copy.
pub(super) fn empty_copied_edge(e: &Edge) -> Edge {
    let geom = crate::tgeometry::GeometryRegistry::global().edge_geom(&e.0);
    let copy = Edge::new();
    if let Some(geom) = geom {
        crate::tgeometry::GeometryRegistry::global().set_edge(&copy.0, geom);
    }
    copy
}

/// `TopoDS_Shape::EmptyCopied()` on a face, with what `BRep_TFace::EmptyCopy`
/// (`BRep_TFace.cxx:36-43`) carries over: the surface, the location and the
/// tolerance - not the wires, not the triangulation, and not
/// `myNaturalRestriction` (a fresh `BRep_TFace` defaults it to `false`,
/// `BRep_TFace.cxx:28-32`).
pub(super) fn empty_copied_face(f: &Face) -> Face {
    let registry = crate::tgeometry::GeometryRegistry::global();
    let mut copy = Face::wrap(f.0.empty_copied()).expect("an empty copy of a face is a face");
    if let Some(geom) = registry.face_geom(&f.0) {
        let mut new_geom = crate::tgeometry::FaceGeom::new(geom.surface);
        new_geom.tolerance = geom.tolerance;
        new_geom.natural_restriction = false;
        registry.set_face(&copy.0, new_geom);
    }
    copy.0.set_orientation(f.0.orientation());
    copy
}

/// `updateToleraces(theFace, theMap)` (`cxx:1375-1407`): for every edge of the
/// face that was split into more than one descendant, raise each descendant's
/// tolerance to the deviation `ShapeAnalysis_Edge::CheckSameParameter` measures
/// on the face.
fn update_tolerances(face: &Face, the_map: &ShapeValueMap<Vec<TopoShape>>) {
    for edge in crate::topo_tools::explorer_shapes(&face.0, ShapeType::Edge) {
        let Some(e) = Edge::wrap(edge) else {
            continue;
        };
        // `cxx:1383-1387`.
        let Some(l_edges) = the_map.find(&e.0) else {
            continue;
        };
        // `cxx:1388-1392`.
        if l_edges.len() <= 1 {
            continue;
        }
        // `cxx:1394-1405`.
        for cur in l_edges {
            let Some(cur) = Edge::wrap(cur.clone()) else {
                continue;
            };
            if let Some(check) = crate::shape_analysis_edge::check_same_parameter(&cur, Some(face)) {
                if check.done {
                    // `BRep_Builder::UpdateEdge(E, Tol)` only raises
                    // (`BRep_TEdge::UpdateTolerance`).
                    crate::tgeometry::GeometryRegistry::global()
                        .raise_edge_tolerance(&cur.0, check.maxdev);
                }
            }
        }
    }
}

/// `BRep_Tool::IsClosed(const TopoDS_Shape& S)` (`BRep_Tool.cxx:1726-1765`):
/// the `Closed` flag, which `TopoDS_Shape::EmptyCopied` does not carry over.
fn brep_tool_is_closed(s: &TopoShape) -> bool {
    match s.shape_type() {
        ShapeType::Wire | ShapeType::Shell => s.is_closed(),
        _ => false,
    }
}
