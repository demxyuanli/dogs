//! `TopOpeBRep_FaceEdgeIntersector` — the intersection of a face and an edge.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRep/TopOpeBRep_FaceEdgeIntersector.hxx:44-137`,
//! `.cxx:26-352`.
//!
//! Translation notes:
//! - `myCurve` is a `GeomAdaptor_Curve` loaded with
//!   `BRep_Tool::Curve(myEdge, loc, f, l)->Transformed(loc)` (`cxx:97-103`);
//!   [`BRepTool::edge_curve_world`] is the same curve with the edge location
//!   applied, and `f` / `l` stay the edge's own range.
//! - `BRepIntCurveSurface_Inter` (`cxx:112-127`) is
//!   [`crate::brep_int_curve_surface_inter::BRepIntCurveSurfaceInter`]; the
//!   loop appends `FEINT.Point()` and maps `FEINT.State()` to `0` for `IN` and
//!   `1` for `ON`.
//! - `mySequenceOfState` stores that integer, so `State()` (`cxx:213-217`) is a
//!   read of the stored byte, not a re-classification.
//! - `Index()` returns `myPointIndex` only under `OCCT_DEBUG` (`cxx:333-339`);
//!   the port has no debug build, so it returns `0` like the release branch.
//! - `myNullShape` / `myNullVertex` (`hxx:135-136`) are `None` here.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d};
use occt_core::precision::INTERSECTION;

use crate::abs::{Orientation, ShapeType, State};
use crate::brep_int_curve_surface_inter::BRepIntCurveSurfaceInter;
use crate::brep_tool::BRepTool;
use crate::intcurvesurface::IntersectionPoint;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topo_tools::explorer_shapes_avoiding;
use crate::topopebrep::ds::transition::Transition;
use crate::topopebrep::tool::shape_tool::tolerance_of;

/// `class TopOpeBRep_FaceEdgeIntersector` (`hxx:44-137`).
pub struct FaceEdgeIntersector {
    /// `myFace` (`hxx:113`).
    face: Face,
    /// `myEdge` (`hxx:114`).
    edge: Edge,
    /// `myTol` (`hxx:115`).
    tol: f64,
    /// `myForceTolerance` (`hxx:116`).
    force_tolerance: bool,
    /// `myIntersectionDone` (`hxx:119`).
    intersection_done: bool,
    /// `mySequenceOfPnt` (`hxx:120`).
    sequence_of_pnt: Vec<IntersectionPoint>,
    /// `mySequenceOfState` (`hxx:121`): `0` = `IN`, `1` = `ON`.
    sequence_of_state: Vec<i32>,
    /// `myPointIndex` (`hxx:122`), 1-based as in OCCT.
    point_index: usize,
    /// `myNbPoints` (`hxx:123`).
    nb_points: usize,
}

impl Default for FaceEdgeIntersector {
    fn default() -> Self {
        Self {
            face: Face::new(),
            edge: Edge::new(),
            tol: 0.0,
            force_tolerance: false,
            intersection_done: false,
            sequence_of_pnt: Vec::new(),
            sequence_of_state: Vec::new(),
            point_index: 0,
            nb_points: 0,
        }
    }
}

impl FaceEdgeIntersector {
    /// `TopOpeBRep_FaceEdgeIntersector()` (`cxx:55-58`).
    pub fn new() -> Self {
        let mut i = Self::default();
        i.reset_intersection();
        i
    }

    /// `ResetIntersection()` (`cxx:62-68`).
    fn reset_intersection(&mut self) {
        self.sequence_of_pnt.clear();
        self.sequence_of_state.clear();
        self.nb_points = 0;
        self.intersection_done = false;
    }

    /// `Perform(SF, SE)` (`cxx:72-131`).
    pub fn perform(&mut self, sf: &TopoShape, se: &TopoShape) {
        self.reset_intersection();
        if !self.force_tolerance {
            self.shape_tolerances(sf, se);
        }
        // myTol = BRep_Tool::Tolerance(TopoDS::Edge(SE)) (cxx:76).
        self.tol = BRepTool::edge_tolerance(&Edge::wrap(se.clone()).unwrap_or_else(Edge::new));

        // myFace.Orientation(FORWARD); myEdge.Orientation(FORWARD) (cxx:79-82).
        self.face = Face::wrap(sf.clone()).unwrap_or_else(Face::new);
        self.face.0.set_orientation(Orientation::Forward);
        self.edge = Edge::wrap(se.clone()).unwrap_or_else(Edge::new);
        self.edge.0.set_orientation(Orientation::Forward);

        // myCurve.Load(C->Transformed(loc), f, l) (cxx:97-103).
        let curve = BRepTool::edge_curve_world(&self.edge);
        let (f, l) = BRepTool::edge_parameters(&self.edge);

        let mut feint = BRepIntCurveSurfaceInter::new();
        let Some(curve) = curve else {
            self.intersection_done = true;
            return;
        };
        let c: Arc<dyn occt_geom::Curve> = curve;
        feint.init_curve_range(c, f, l);

        while feint.more() {
            if let Ok(p) = feint.point() {
                self.sequence_of_pnt.push(p);
                let i = if feint.state() == Ok(State::In) { 0 } else { 1 };
                self.sequence_of_state.push(i);
            }
            feint.next();
        }

        self.nb_points = self.sequence_of_pnt.len();
        self.intersection_done = true;
    }

    /// `IsEmpty()` (`cxx:135-139`).
    pub fn is_empty(&self) -> bool {
        self.nb_points == 0
    }

    /// `Shape(I)` (`cxx:143-158`): `I` = 1 for the face, 2 for the edge.
    pub fn shape(&self, index: i32) -> Result<&TopoShape, String> {
        if index == 1 {
            Ok(&self.face.0)
        } else if index == 2 {
            Ok(&self.edge.0)
        } else {
            Err("TopOpeBRep_FaceEdgeIntersector::Shape".to_string())
        }
    }

    /// `ForceTolerance(Tol)` (`cxx:162-169`).
    pub fn force_tolerance(&mut self, tol: f64) {
        self.tol = tol;
        self.force_tolerance = true;
    }

    /// `Tolerance()` (`cxx:173-176`).
    pub fn tolerance(&self) -> f64 {
        self.tol
    }

    /// `NbPoints()` (`cxx:180-184`).
    pub fn nb_points(&self) -> i32 {
        self.nb_points as i32
    }

    /// `InitPoint()` (`cxx:188-191`).
    pub fn init_point(&mut self) {
        self.point_index = 1;
    }

    /// `MorePoint()` (`cxx:195-199`).
    pub fn more_point(&self) -> bool {
        self.point_index <= self.nb_points
    }

    /// `NextPoint()` (`cxx:203-206`).
    pub fn next_point(&mut self) {
        self.point_index += 1;
    }

    /// The point `mySequenceOfPnt(myPointIndex)`.
    fn current(&self) -> Option<&IntersectionPoint> {
        if self.point_index == 0 || self.point_index > self.sequence_of_pnt.len() {
            return None;
        }
        self.sequence_of_pnt.get(self.point_index - 1)
    }

    /// `Value()` (`cxx:210-215`).
    pub fn value(&self) -> GpPnt {
        match self.current() {
            Some(p) => p.pnt(),
            None => GpPnt::new(0.0, 0.0, 0.0),
        }
    }

    /// `Parameter()` (`cxx:219-224`): the parameter on the edge.
    pub fn parameter(&self) -> f64 {
        match self.current() {
            Some(p) => p.param(),
            None => 0.0,
        }
    }

    /// `UVPoint(P)` (`cxx:228-234`).
    pub fn uv_point(&self) -> Result<GpPnt2d, String> {
        match self.current() {
            Some(p) => Ok(GpPnt2d::new(p.u(), p.v())),
            None => Err("TopOpeBRep_FaceEdgeIntersector::UVPoint".to_string()),
        }
    }

    /// `State()` (`cxx:238-243`): the stored `0` / `1` means `IN` / `ON`.
    pub fn state(&self) -> State {
        match self.sequence_of_state.get(self.point_index.wrapping_sub(1)) {
            Some(0) => State::In,
            Some(_) => State::On,
            None => State::Unknown,
        }
    }

    /// `Transition(Index, FaceOrientation)` (`cxx:247-320`).
    pub fn transition(&self, index: i32, face_orientation: Orientation) -> Result<Transition, String> {
        let ip = self
            .current()
            .ok_or_else(|| "TopOpeBRep_FaceEdgeIntersector::Transition".to_string())?;

        if index == 2 {
            // Edge IN <=> entering the material of the face.
            let (st_b, st_a) = match ip.state() {
                crate::intcurvesurface::State::In => (State::Out, State::In),
                crate::intcurvesurface::State::Out => (State::In, State::Out),
                _ => (State::In, State::In),
            };
            let mut tr = Transition::default();
            let (on_b, on_a) = (ShapeType::Face, ShapeType::Face);
            if face_orientation == Orientation::Forward {
                tr.set(st_b, st_a, on_b, on_a);
            } else if face_orientation == Orientation::Reversed {
                tr.set(st_a, st_b, on_a, on_b);
            } else if face_orientation == Orientation::External {
                tr.set(State::Out, State::Out, on_a, on_b);
            } else if face_orientation == Orientation::Internal {
                tr.set(State::In, State::In, on_a, on_b);
            }
            return Ok(tr);
        }

        if index == 1 {
            // The face is always inside the face.
            let (on_b, on_a) = (ShapeType::Face, ShapeType::Face);
            let mut tr = Transition::default();
            tr.set(State::In, State::In, on_b, on_a);
            return Ok(tr);
        }

        Err("FEINT Transition Index".to_string())
    }

    /// `IsVertex(S, P, Tol, VR)` (`cxx:324-345`).
    pub fn is_vertex_on_shape(
        &self,
        s: &TopoShape,
        p: &GpPnt,
        tol: f64,
    ) -> Option<Vertex> {
        let tol2 = tol * tol;
        let mut found: Option<Vertex> = None;
        for ss in explorer_shapes_avoiding(s, ShapeType::Vertex, ShapeType::Shape) {
            let Some(vv) = Vertex::wrap(ss) else {
                continue;
            };
            let pv = BRepTool::vertex_point(&vv);
            if p.square_distance(&pv) < tol2 {
                found = Some(vv);
            }
        }
        found
    }

    /// `IsVertex(I, VR)` (`cxx:349-363`): the current value tested against the
    /// face (`I == 1`) or the edge (`I == 2`).
    pub fn is_vertex(&self, i: i32) -> Option<Vertex> {
        let p = self.value();
        if i == 1 {
            self.is_vertex_on_shape(&self.face.0, &p, self.tol)
        } else if i == 2 {
            self.is_vertex_on_shape(&self.edge.0, &p, self.tol)
        } else {
            None
        }
    }

    /// `Index()` (`cxx:367-374`): `OCCT_DEBUG` only, so `0`.
    pub fn index(&self) -> i32 {
        0
    }

    /// `ShapeTolerances(S1, S2)` (`cxx:378-390`).
    fn shape_tolerances(&mut self, s1: &TopoShape, s2: &TopoShape) {
        self.tol = self
            .tolerance_max(s1, ShapeType::Edge)
            .max(self.tolerance_max(s2, ShapeType::Edge));
        self.force_tolerance = false;
    }

    /// `ToleranceMax(S, T)` (`cxx:394-410`): `Precision::Intersection()` when
    /// `S` has no sub-shape of type `T`.
    fn tolerance_max(&self, s: &TopoShape, t: ShapeType) -> f64 {
        let subs = explorer_shapes_avoiding(s, t, ShapeType::Shape);
        if subs.is_empty() {
            return INTERSECTION;
        }
        let mut tol = f64::NEG_INFINITY;
        for sub in subs {
            if let Ok(v) = tolerance_of(&sub) {
                tol = tol.max(v);
            }
        }
        tol
    }
}

impl std::fmt::Debug for FaceEdgeIntersector {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("FaceEdgeIntersector")
            .field("tol", &self.tol)
            .field("nb_points", &self.nb_points)
            .field("intersection_done", &self.intersection_done)
            .finish()
    }
}
