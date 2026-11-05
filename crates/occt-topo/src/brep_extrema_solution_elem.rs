//! `BRepExtrema_SolutionElem` and `BRepExtrema_SupportType`.
//!
//! Source: `ModelingAlgorithms/TKTopAlgo/BRepExtrema/BRepExtrema_SolutionElem.hxx:23-143`
//! and `BRepExtrema_SupportType.hxx` (OCCT 8.0.0). Both are header-only.
//!
//! One element of the solution sequence of `BRepExtrema_DistanceSS`: the
//! distance, the 3D support point, what kind of support it is (vertex / on edge
//! / in face), the support sub-shape and the location parameters.
//!
//! # Translation notes
//!
//! * The four OCCT constructors (`hxx:30-37`, `:44-56`, `:64-75`, `:85-101`)
//!   differ only in which support sub-shape and which parameters they write;
//!   the port keeps them as four constructors so every call site keeps its OCCT
//!   spelling.
//! * OCCT leaves `myVertex` / `myEdge` / `myFace` at the *default* `TopoDS`
//!   instance (a null shape) for the constructors that do not take them, and
//!   `myPar1` / `myPar2` at `0`. The port stores the same: `None` for the null
//!   shape, `0.0` for the parameters.
//! * `Vertex()` / `Edge()` / `Face()` are unchecked `TopoDS::Vertex`-style
//!   reads in OCCT (they return the stored handle without a type test); the
//!   port returns `Option` because the stored shape can be null.

use occt_core::gp::GpPnt;

use crate::shape::{Edge, Face, TopoShape, Vertex};

/// `BRepExtrema_SupportType` (`BRepExtrema_SupportType.hxx:21-27`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SupportType {
    /// `BRepExtrema_IsVertex`: the solution is a vertex.
    IsVertex,
    /// `BRepExtrema_IsOnEdge`: the solution belongs to an edge.
    IsOnEdge,
    /// `BRepExtrema_IsInFace`: the solution is inside a face.
    IsInFace,
}

/// `BRepExtrema_SolutionElem` (`BRepExtrema_SolutionElem.hxx:23-143`).
#[derive(Debug, Clone)]
pub struct SolutionElem {
    /// `myDist` (`hxx:134`).
    my_dist: f64,
    /// `myPoint` (`hxx:135`).
    my_point: GpPnt,
    /// `mySupType` (`hxx:136`).
    my_sup_type: SupportType,
    /// `myVertex` (`hxx:137`).
    my_vertex: Option<TopoShape>,
    /// `myEdge` (`hxx:138`).
    my_edge: Option<TopoShape>,
    /// `myFace` (`hxx:139`).
    my_face: Option<TopoShape>,
    /// `myPar1` (`hxx:140`).
    my_par1: f64,
    /// `myPar2` (`hxx:141`).
    my_par2: f64,
}

impl Default for SolutionElem {
    /// `BRepExtrema_SolutionElem()` (`hxx:30-37`): distance `0`, point at the
    /// origin, support kind `IsVertex`, null sub-shapes, parameters `0`.
    fn default() -> Self {
        Self {
            my_dist: 0.0,
            my_point: GpPnt::zero(),
            my_sup_type: SupportType::IsVertex,
            my_vertex: None,
            my_edge: None,
            my_face: None,
            my_par1: 0.0,
            my_par2: 0.0,
        }
    }
}

impl SolutionElem {
    /// `BRepExtrema_SolutionElem(theDist, thePoint, theSolType, theVertex)`
    /// (`hxx:44-56`).
    pub fn of_vertex(
        the_dist: f64,
        the_point: &GpPnt,
        the_sol_type: SupportType,
        the_vertex: &Vertex,
    ) -> Self {
        Self {
            my_dist: the_dist,
            my_point: *the_point,
            my_sup_type: the_sol_type,
            my_vertex: Some(the_vertex.0.clone()),
            ..Self::default()
        }
    }

    /// `BRepExtrema_SolutionElem(theDist, thePoint, theSolType, theEdge,
    /// theParam)` (`hxx:64-75`).
    pub fn of_edge(
        the_dist: f64,
        the_point: &GpPnt,
        the_sol_type: SupportType,
        the_edge: &Edge,
        the_param: f64,
    ) -> Self {
        Self {
            my_dist: the_dist,
            my_point: *the_point,
            my_sup_type: the_sol_type,
            my_edge: Some(the_edge.0.clone()),
            my_par1: the_param,
            ..Self::default()
        }
    }

    /// `BRepExtrema_SolutionElem(theDist, thePoint, theSolType, theFace, theU,
    /// theV)` (`hxx:85-101`).
    pub fn of_face(
        the_dist: f64,
        the_point: &GpPnt,
        the_sol_type: SupportType,
        the_face: &Face,
        the_u: f64,
        the_v: f64,
    ) -> Self {
        Self {
            my_dist: the_dist,
            my_point: *the_point,
            my_sup_type: the_sol_type,
            my_face: Some(the_face.0.clone()),
            my_par1: the_u,
            my_par2: the_v,
            ..Self::default()
        }
    }

    /// `Dist()` (`hxx:103`).
    pub fn dist(&self) -> f64 {
        self.my_dist
    }

    /// `Point()` (`hxx:106`).
    pub fn point(&self) -> &GpPnt {
        &self.my_point
    }

    /// `SupportKind()` (`hxx:112`).
    pub fn support_kind(&self) -> SupportType {
        self.my_sup_type
    }

    /// `Vertex()` (`hxx:116`): the stored vertex, `None` when the constructor
    /// did not take one (OCCT returns the null handle).
    pub fn vertex(&self) -> Option<Vertex> {
        self.my_vertex.as_ref().and_then(|s| Vertex::wrap(s.clone()))
    }

    /// `Edge()` (`hxx:119`).
    pub fn edge(&self) -> Option<Edge> {
        self.my_edge.as_ref().and_then(|s| Edge::wrap(s.clone()))
    }

    /// `Face()` (`hxx:122`).
    pub fn face(&self) -> Option<Face> {
        self.my_face.as_ref().and_then(|s| Face::wrap(s.clone()))
    }

    /// `EdgeParameter(theParam)` (`hxx:125-127`).
    pub fn edge_parameter(&self) -> f64 {
        self.my_par1
    }

    /// `FaceParameter(theU, theV)` (`hxx:129-133`).
    pub fn face_parameter(&self) -> (f64, f64) {
        (self.my_par1, self.my_par2)
    }
}
