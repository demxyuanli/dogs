//! `BRepAdaptor_Curve2d` — an edge's p-curve on a face seen as a 2D curve.
//!
//! Source: `ModelingData/TKBRep/BRepAdaptor/BRepAdaptor_Curve2d.hxx:41-79`,
//! `.cxx:19-77`.
//!
//! Translation notes:
//! - `Geom2dAdaptor_Curve` has no port class; the port's [`Curve2d`] trait
//!   already carries `D0` / `D1` / `FirstParameter` / `LastParameter` /
//!   `Resolution`, so [`Curve2dAdaptor`] stores the `Arc<dyn Curve2d>` that
//!   `Initialize` loads and forwards the adaptor queries to it.
//! - `Initialize(E, F)` (`cxx:51-62`) reads `BRep_Tool::CurveOnSurface(E, F,
//!   first, last)` and only loads the adaptor when that curve is not null
//!   (`cxx:55-59`); [`crate::boptools_2d::curve_on_surface_range`] is the port
//!   of that call. A null p-curve therefore leaves the adaptor unloaded, which
//!   is what a default-constructed `BRepAdaptor_Curve2d` is.
//! - `Edge()` / `Face()` (`cxx:69-77`) hand back the two stored shapes; the
//!   default constructor leaves both null, mirrored by `Edge::new()` /
//!   `Face::new()` (a null `TopoDS_Edge` / `TopoDS_Face`).

use std::sync::Arc;

use occt_core::gp::{GpPnt2d, GpVec2d};
use occt_geom2d::curve::Curve2d;

use crate::shape::{Edge, Face};

/// `class BRepAdaptor_Curve2d : public Geom2dAdaptor_Curve` (`hxx:41-79`).
#[derive(Clone)]
pub struct Curve2dAdaptor {
    /// `Geom2dAdaptor_Curve::myCurve` — `None` when no p-curve was loaded.
    curve: Option<Arc<dyn Curve2d>>,
    /// `Geom2dAdaptor_Curve::myFirst`.
    first: f64,
    /// `Geom2dAdaptor_Curve::myLast`.
    last: f64,
    /// `myEdge` (`hxx:76`).
    edge: Edge,
    /// `myFace` (`hxx:77`).
    face: Face,
}

impl Default for Curve2dAdaptor {
    /// `BRepAdaptor_Curve2d::BRepAdaptor_Curve2d() = default` (`cxx:27`).
    fn default() -> Self {
        Self::new()
    }
}

impl Curve2dAdaptor {
    /// `BRepAdaptor_Curve2d()` (`cxx:27`): no curve, no edge, no face.
    pub fn new() -> Self {
        Self {
            curve: None,
            first: 0.0,
            last: 0.0,
            edge: Edge::new(),
            face: Face::new(),
        }
    }

    /// `BRepAdaptor_Curve2d(E, F)` (`cxx:31-34`).
    pub fn of_edge_face(e: &Edge, f: &Face) -> Self {
        let mut a = Self::new();
        a.initialize(e, f);
        a
    }

    /// An adaptor over an already-extracted p-curve.
    ///
    /// `Initialize` (`cxx:51-62`) reads the p-curve out of
    /// `BRep_Tool::CurveOnSurface`; `GeomInt` keeps the p-curve and its range
    /// instead of the (edge, face) pair, so the same adaptor is assembled from
    /// those stored pieces. `edge` is the arc's edge (`cxx:69-72` still answers
    /// it), the face stays null because no p-curve lookup happens.
    pub fn from_stored(e: &Edge, curve: Arc<dyn Curve2d>, first: f64, last: f64) -> Self {
        Self {
            curve: Some(curve),
            first,
            last,
            edge: e.clone(),
            face: Face::new(),
        }
    }

    /// `Initialize(E, F)` (`cxx:51-62`).
    pub fn initialize(&mut self, e: &Edge, f: &Face) {
        self.edge = e.clone();
        self.face = f.clone();
        if let Some((pc, first, last)) = crate::boptools_2d::curve_on_surface_range(e, f) {
            // `Geom2dAdaptor_Curve::Load(aPCurve, aFirs, aLast)`.
            self.curve = Some(pc);
            self.first = first;
            self.last = last;
        }
    }

    /// `Edge()` (`cxx:69-72`).
    pub fn edge(&self) -> &Edge {
        &self.edge
    }

    /// `Face()` (`cxx:74-77`).
    pub fn face(&self) -> &Face {
        &self.face
    }

    /// `Geom2dAdaptor_Curve::Curve()`: the loaded p-curve, if any.
    pub fn curve(&self) -> Option<&Arc<dyn Curve2d>> {
        self.curve.as_ref()
    }

    /// `Adaptor2d_Curve2d::FirstParameter()`.
    pub fn first_parameter(&self) -> f64 {
        self.first
    }

    /// `Adaptor2d_Curve2d::LastParameter()`.
    pub fn last_parameter(&self) -> f64 {
        self.last
    }

    /// `Adaptor2d_Curve2d::D0(U, P)`.
    pub fn d0(&self, u: f64) -> Option<GpPnt2d> {
        self.curve.as_ref().map(|c| c.d0(u))
    }

    /// `Adaptor2d_Curve2d::D1(U, P, D)`.
    pub fn d1(&self, u: f64) -> Option<(GpPnt2d, GpVec2d)> {
        self.curve.as_ref().map(|c| c.d1(u))
    }

    /// `Adaptor2d_Curve2d::Resolution(Ruv)` (`Geom2dAdaptor_Curve.cxx:1186-1219`
    /// for the loaded curve's own type). An unloaded adaptor has no resolution.
    pub fn resolution(&self, ruv: f64) -> Option<f64> {
        self.curve.as_ref().map(|c| c.resolution(ruv))
    }
}
