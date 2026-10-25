//! `TopOpeBRepDS_DataStructure` — the central store of the classic boolean
//! engine: points, curves, surfaces, shapes, section edges, same-domain
//! relations and the per-shape state maps.
//!
//! Sources:
//! - `TKBool/TopOpeBRepDS/TopOpeBRepDS_DataStructure.hxx:52-349`
//! - `TKBool/TopOpeBRepDS/TopOpeBRepDS_DataStructure.cxx:33-1355`
//!
//! OCCT keeps `NCollection_DataMap<int, *Data>` for the three geometry tables
//! (keyed 1..NbX by incrementing counters, see `AddSurface`/`AddCurve`/
//! `AddPoint`) and `NCollection_IndexedDataMap<TopoDS_Shape, ...>` for the
//! shapes. The Rust translation keeps `BTreeMap<i32, _>` for the geometry
//! tables (a `DataMap` has no iteration order to preserve) and
//! [`IndexedShapeData`] for the shape tables (which do: OCCT indexes them
//! 1..Extent).
//!
//! Three OCCT members are "returned by reference even when absent": the
//! `myEmpty*` scratch values (`myEmptyShape`, `myEmptyPoint`,
//! `myEmptySurface`, `myEmptyCurve`, `myEmptyListOfInterference`,
//! `myEmptyListOfShape`, `myEmptyGSurface`) and the `static` dummies of
//! `ChangeMapOfShapeWithState` / `GetShapeWithState`. Where OCCT hands back
//! such a value the Rust side returns `None`, except for the interference and
//! same-domain lists, which stay borrowed slices (`&[]`) and a mutable
//! swallow-list, because `AddShapeInterference` relies on appending to a
//! throw-away list when the shape is absent (`cxx:804-809`).

mod index_map;
mod shape_with_state_map;
mod shapes;

use std::collections::{BTreeMap, HashMap};
use std::sync::Arc;

use occt_geom::Surface as GeomSurface;

use crate::shape::TopoShape;
use crate::topo_tools::ShapeMapKey;

use super::curve::Curve;
use super::curve_data::CurveData;
use super::interference::InterferenceHandle;
use super::point::Point;
use super::point_data::PointData;
use super::shape_data::ShapeData;
use super::shape_with_state::ShapeWithState;
use super::surface::Surface;
use super::surface_data::SurfaceData;

pub use index_map::{IndexedShapeData, IndexedShapeSet};
pub use shape_with_state_map::ShapeWithStateMap;

/// Source: `TopOpeBRepDS_DataStructure.hxx:320-347`.
#[derive(Clone)]
pub struct DataStructure {
    /// `myNbSurfaces`
    nb_surfaces: i32,
    /// `mySurfaces`
    surfaces: BTreeMap<i32, SurfaceData>,
    /// `myNbCurves`
    nb_curves: i32,
    /// `myCurves`
    curves: BTreeMap<i32, CurveData>,
    /// `myNbPoints`
    nb_points: i32,
    /// `myPoints`
    points: BTreeMap<i32, PointData>,
    /// `myShapes`
    shapes: IndexedShapeData<ShapeData>,
    /// `mySectionEdges`
    section_edges: IndexedShapeSet,
    /// `myEmptyListOfInterference` — what
    /// `ChangeXInterferences` appends to when the index is unbound
    /// (`cxx:384-393`, `:412-424`, `:442-453`).
    empty_interferences: Vec<InterferenceHandle>,
    /// `myEmptyShape` — no longer a value: the accessors return `None`.
    /// `myEmptyPoint` / `myEmptySurface` / `myEmptyCurve` are the scratch
    /// values returned by `Point` / `ChangeSurface` / `ChangeCurve` when the
    /// index is unbound (`cxx:1052-1106`, `:1128-1144`).
    ///
    /// `myEmptyPoint` is only reachable through an index that is within
    /// `1..=myNbPoints` but unbound, which the three `Bind` sites cannot
    /// produce, so the field stays unread; it is kept because OCCT keeps it.
    #[allow(dead_code)]
    empty_point: Point,
    /// `myEmptySurface`
    empty_surface: Surface,
    /// `myEmptyCurve`
    empty_curve: Curve,
    /// `myNewSurface` — `myEmptyGSurface` (the null handle) is `None`.
    new_surface: HashMap<ShapeMapKey, Arc<dyn GeomSurface>>,
    /// `myIsfafa`
    isfafa: bool,
    /// `myI` — an unread scratch member in OCCT too (`hxx:340`; nothing in
    /// `DataStructure.cxx` reads it after the ctor sets it to 0).
    #[allow(dead_code)]
    i: i32,
    /// `myMapOfShapeWithStateObj`
    map_of_shape_with_state_obj: IndexedShapeData<ShapeWithState>,
    /// `myMapOfShapeWithStateTool`
    map_of_shape_with_state_tool: IndexedShapeData<ShapeWithState>,
    /// `myMapOfRejectedShapesObj`
    map_of_rejected_shapes_obj: IndexedShapeSet,
    /// `myMapOfRejectedShapesTool`
    map_of_rejected_shapes_tool: IndexedShapeSet,
}

impl Default for DataStructure {
    /// `TopOpeBRepDS_DataStructure()` (`cxx:33-40`).
    fn default() -> Self {
        Self::new()
    }
}

impl DataStructure {
    /// `TopOpeBRepDS_DataStructure()` (`cxx:33-40`): the counters and
    /// `myIsfafa`/`myI` are in the init list, everything else default.
    pub fn new() -> Self {
        Self {
            nb_surfaces: 0,
            surfaces: BTreeMap::new(),
            nb_curves: 0,
            curves: BTreeMap::new(),
            nb_points: 0,
            points: BTreeMap::new(),
            shapes: IndexedShapeData::new(),
            section_edges: IndexedShapeSet::new(),
            empty_interferences: Vec::new(),
            empty_point: Point::new(),
            empty_surface: Surface::new(),
            empty_curve: Curve::new(),
            new_surface: HashMap::new(),
            isfafa: false,
            i: 0,
            map_of_shape_with_state_obj: IndexedShapeData::new(),
            map_of_shape_with_state_tool: IndexedShapeData::new(),
            map_of_rejected_shapes_obj: IndexedShapeSet::new(),
            map_of_rejected_shapes_tool: IndexedShapeSet::new(),
        }
    }

    /// `Init()` (`cxx:44-61`).
    pub fn init(&mut self) {
        self.nb_surfaces = 0;
        self.nb_curves = 0;
        self.nb_points = 0;
        self.surfaces.clear();
        self.curves.clear();
        self.points.clear();
        self.shapes.clear();
        // cxx:54-57 (modified by NIZHNY-MZV  Tue Apr 18 16:33:26 2000)
        self.map_of_shape_with_state_obj.clear();
        self.map_of_shape_with_state_tool.clear();
        self.map_of_rejected_shapes_obj.clear();
        self.map_of_rejected_shapes_tool.clear();
        self.init_section_edges();
    }

    /// `NbSurfaces()` (`cxx:1005-1010`).
    pub fn nb_surfaces(&self) -> i32 {
        self.nb_surfaces
    }

    /// `NbCurves()` (`cxx:1013-1018`).
    pub fn nb_curves(&self) -> i32 {
        self.nb_curves
    }

    /// `ChangeNbCurves(N)` (`cxx:1021-1026`).
    pub fn change_nb_curves(&mut self, n: i32) {
        self.nb_curves = n;
    }

    /// `NbPoints()` (`cxx:1028-1033`).
    pub fn nb_points(&self) -> i32 {
        self.nb_points
    }

    /// `Isfafa()` (`cxx:1278-1282`).
    pub fn isfafa(&self) -> bool {
        self.isfafa
    }

    /// `Isfafa(bool)` (`cxx:1271-1275`).
    pub fn set_isfafa(&mut self, isfafa: bool) {
        self.isfafa = isfafa;
    }

    // --- surfaces ---------------------------------------------------------

    /// `AddSurface(S)` (`cxx:64-73`).
    pub fn add_surface(&mut self, s: &Surface) -> i32 {
        self.nb_surfaces += 1;
        self.surfaces
            .insert(self.nb_surfaces, SurfaceData::with_surface(s.clone()));
        self.nb_surfaces
    }

    /// `RemoveSurface(I)` (`cxx:77-82`): `mySurfaces.UnBind(I)`.
    pub fn remove_surface(&mut self, i: i32) {
        self.surfaces.remove(&i);
    }

    /// `KeepSurface(I)` (`cxx:84-90`).
    pub fn keep_surface(&self, i: i32) -> bool {
        self.surface(i).keep()
    }

    /// `KeepSurface(TopOpeBRepDS_Surface&)` (`cxx:93-98`).
    pub fn keep_surface_of(&self, s: &Surface) -> bool {
        s.keep()
    }

    /// `ChangeKeepSurface(I, FindKeep)` (`cxx:101-106`).
    pub fn change_keep_surface(&mut self, i: i32, find_keep: bool) {
        self.change_surface(i).change_keep(find_keep);
    }

    /// `ChangeKeepSurface(TopOpeBRepDS_Surface&, FindKeep)` (`cxx:109-113`).
    pub fn change_keep_surface_of(&self, s: &mut Surface, find_keep: bool) {
        s.change_keep(find_keep);
    }

    /// `Surface(I)` (`cxx:1052-1063`): `myEmptySurface` when unbound.
    pub fn surface(&self, i: i32) -> &Surface {
        match self.surfaces.get(&i) {
            Some(sd) => sd.surface(),
            None => &self.empty_surface,
        }
    }

    /// `ChangeSurface(I)` (`cxx:1066-1077`): `myEmptySurface` when unbound,
    /// so a write through it is discarded exactly as in OCCT.
    pub fn change_surface(&mut self, i: i32) -> &mut Surface {
        match self.surfaces.get_mut(&i) {
            Some(sd) => sd.surface_mut(),
            None => &mut self.empty_surface,
        }
    }

    /// `mySurfaces.IsBound(I)` — `TopOpeBRepDS_SurfaceExplorer` is a `friend`
    /// class (`TopOpeBRepDS_DataStructure.hxx:311`).
    pub(crate) fn is_surface_bound(&self, i: i32) -> bool {
        self.surfaces.contains_key(&i)
    }

    // --- curves -----------------------------------------------------------

    /// `AddCurve(C)` (`cxx:116-130`).
    ///
    /// OCCT casts the `const&` away to store the new index in the caller's
    /// curve (`cxx:120-121`), hence `&mut` here; the stored copy is then
    /// taken (`TopOpeBRepDS_CurveData CD(C)`).
    pub fn add_curve(&mut self, c: &mut Curve) -> i32 {
        self.nb_curves += 1;
        c.change_ds_index(self.nb_curves);
        self.curves
            .insert(self.nb_curves, CurveData::with_curve(c.clone()));
        self.nb_curves
    }

    /// `RemoveCurve(I)` (`cxx:132-152`).
    pub fn remove_curve(&mut self, i: i32) {
        let Some(cd) = self.curves.get_mut(&i) else {
            // `ChangeCurve(I)` yields `myEmptyCurve`: no shapes, no
            // interferences, so every step below is a no-op in OCCT too.
            return;
        };
        let (s1, s2) = match cd.curve().get_shapes() {
            (Some(a), Some(b)) => (Some(a.clone()), Some(b.clone())),
            (Some(a), None) => (Some(a.clone()), None),
            (None, Some(b)) => (None, Some(b.clone())),
            (None, None) => (None, None),
        };
        let i1 = cd.curve().get_sci1().cloned();
        let i2 = cd.curve().get_sci2().cloned();
        cd.curve_mut().change_keep(false);
        if let (Some(s1), Some(i1)) = (s1, i1) {
            self.remove_shape_interference(&s1, &i1);
        }
        if let (Some(s2), Some(i2)) = (s2, i2) {
            self.remove_shape_interference(&s2, &i2);
        }
        // cxx:150: "NYI : nullify interferences I1,I2"
    }

    /// `KeepCurve(I)` (`cxx:155-160`).
    pub fn keep_curve(&self, i: i32) -> bool {
        self.curve(i).keep()
    }

    /// `KeepCurve(const TopOpeBRepDS_Curve&)` (`cxx:164-168`).
    pub fn keep_curve_of(&self, c: &Curve) -> bool {
        c.keep()
    }

    /// `ChangeKeepCurve(I, FindKeep)` (`cxx:172-176`).
    pub fn change_keep_curve(&mut self, i: i32, find_keep: bool) {
        self.change_curve(i).change_keep(find_keep);
    }

    /// `ChangeKeepCurve(TopOpeBRepDS_Curve&, FindKeep)` (`cxx:180-183`).
    pub fn change_keep_curve_of(&self, c: &mut Curve, find_keep: bool) {
        c.change_keep(find_keep);
    }

    /// `Curve(I)` (`cxx:1080-1092`): `myEmptyCurve` when unbound.
    pub fn curve(&self, i: i32) -> &Curve {
        match self.curves.get(&i) {
            Some(cd) => cd.curve(),
            None => &self.empty_curve,
        }
    }

    /// `ChangeCurve(I)` (`cxx:1096-1106`): `myEmptyCurve` when unbound.
    pub fn change_curve(&mut self, i: i32) -> &mut Curve {
        match self.curves.get_mut(&i) {
            Some(cd) => cd.curve_mut(),
            None => &mut self.empty_curve,
        }
    }

    /// `myCurves.IsBound(I)` — `TopOpeBRepDS_CurveExplorer` is a `friend`
    /// class (`TopOpeBRepDS_DataStructure.hxx:312`).
    pub(crate) fn is_curve_bound(&self, i: i32) -> bool {
        self.curves.contains_key(&i)
    }

    // --- points -----------------------------------------------------------

    /// `AddPoint(PDS)` (`cxx:187-192`).
    pub fn add_point(&mut self, pds: &Point) -> i32 {
        self.nb_points += 1;
        self.points
            .insert(self.nb_points, PointData::with_point(pds.clone()));
        self.nb_points
    }

    /// `AddPointSS(PDS, S1, S2)` (`cxx:196-202`): the two shapes are ignored.
    pub fn add_point_ss(&mut self, pds: &Point, _s1: &TopoShape, _s2: &TopoShape) -> i32 {
        self.add_point(pds)
    }

    /// `RemovePoint(I)` (`cxx:206-211`).
    pub fn remove_point(&mut self, i: i32) {
        // `ChangePoint(I)` is a `friend` view of `PointData::myPoint`.
        if let Some(pd) = self.points.get_mut(&i) {
            pd.point.change_keep(false);
        }
    }

    /// `KeepPoint(I)` (`cxx:214-219`).
    pub fn keep_point(&self, i: i32) -> bool {
        match self.point(i) {
            Some(p) => p.keep(),
            None => false,
        }
    }

    /// `KeepPoint(const TopOpeBRepDS_Point&)` (`cxx:223-227`).
    pub fn keep_point_of(&self, p: &Point) -> bool {
        p.keep()
    }

    /// `ChangeKeepPoint(I, FindKeep)` (`cxx:231-235`).
    pub fn change_keep_point(&mut self, i: i32, find_keep: bool) {
        if let Some(p) = self.change_point(i) {
            p.change_keep(find_keep);
        }
    }

    /// `ChangeKeepPoint(TopOpeBRepDS_Point&, FindKeep)` (`cxx:239-242`).
    pub fn change_keep_point_of(&self, p: &mut Point, find_keep: bool) {
        p.change_keep(find_keep);
    }

    /// `Point(I)` (`cxx:1109-1124`).
    ///
    /// OCCT throws `Standard_ProgramError` for `I` outside `1..=myNbPoints`
    /// and returns `myEmptyPoint` for an in-range but unbound index; both
    /// cases are `None` here.
    pub fn point(&self, i: i32) -> Option<&Point> {
        if i < 1 || i > self.nb_points {
            return None;
        }
        self.points.get(&i).map(|pd| &pd.point)
    }

    /// `ChangePoint(I)` (`cxx:1128-1144`); same range rule as [`Self::point`].
    pub fn change_point(&mut self, i: i32) -> Option<&mut Point> {
        if i < 1 || i > self.nb_points {
            return None;
        }
        self.points.get_mut(&i).map(|pd| &mut pd.point)
    }

    /// `myPoints.IsBound(I)` — `TopOpeBRepDS_PointExplorer` is a `friend`
    /// class (`TopOpeBRepDS_DataStructure.hxx:313`).
    pub(crate) fn is_point_bound(&self, i: i32) -> bool {
        self.points.contains_key(&i)
    }

    // --- new surfaces -----------------------------------------------------

    /// `SetNewSurface(F, S)` (`cxx:1228-1233`).
    pub fn set_new_surface(&mut self, f: &TopoShape, s: Arc<dyn GeomSurface>) {
        self.new_surface.insert(crate::topo_tools::shape_map_key(f), s);
    }

    /// `HasNewSurface(F)` (`cxx:1236-1241`).
    pub fn has_new_surface(&self, f: &TopoShape) -> bool {
        self.new_surface
            .contains_key(&crate::topo_tools::shape_map_key(f))
    }

    /// `NewSurface(F)` (`cxx:1244-1252`): `myEmptyGSurface` (the null handle)
    /// is `None`.
    pub fn new_surface(&self, f: &TopoShape) -> Option<&Arc<dyn GeomSurface>> {
        self.new_surface
            .get(&crate::topo_tools::shape_map_key(f))
    }
}
