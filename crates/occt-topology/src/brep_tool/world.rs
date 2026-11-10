use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_geom::{Curve, Surface};

use crate::shape::{Edge, Face, Vertex};
use crate::tgeometry::GeometryRegistry;

use super::BRepTool;

impl BRepTool {
    /// World-coordinate vertex point: the registered local point transformed
    /// by the shape's accumulated `location` (`TopoDS_Shape::Location`).
    pub fn vertex_point_world(v: &Vertex) -> GpPnt {
        let t = v.0.location().transformation();
        Self::vertex_point(v).transformed(&t)
    }

    /// World-coordinate edge curve: the registered local curve transformed by
    /// the edge's location.
    pub fn edge_curve_world(e: &Edge) -> Option<Arc<dyn Curve>> {
        let t = e.0.location().transformation();
        Self::edge_curve(e).map(|c| Arc::from(c.transformed(&t)))
    }

    /// World-coordinate edge curve as `BRepAdaptor_Curve(E)` evaluates it
    /// (`BRepAdaptor_Curve.cxx`: `myCurve` materialized from
    /// `BRep_Tool::Curve(E, L)` when the edge has a `BRep_Curve3D`, otherwise
    /// the `Adaptor3d_CurveOnSurface` built from the pcurve, and `myTrsf = L`
    /// applied to every returned point). `edge_curve_world` cannot serve here:
    /// it answers `None` on a pcurve-only edge (`curve3d_null`), which the
    /// adaptor still evaluates.
    pub fn edge_curve_adaptor_world(e: &Edge) -> Option<Arc<dyn Curve>> {
        let geom = GeometryRegistry::global().edge_geom(&e.0)?;
        let loc = e.0.location();
        if loc.is_identity() {
            return Some(geom.curve);
        }
        Some(Arc::from(geom.curve.transformed(&loc.transformation())))
    }

    /// World-coordinate face surface: the registered local surface transformed
    /// by the face's location.
    pub fn face_surface_world(f: &Face) -> Option<Arc<dyn Surface>> {
        let t = f.0.location().transformation();
        Self::face_surface(f).map(|s| Arc::from(s.transformed(&t)))
    }

    /// World-coordinate face UV bounds (from the transformed surface).
    pub fn uv_bounds_world(f: &Face) -> (f64, f64, f64, f64) {
        match Self::face_surface_world(f) {
            Some(s) => {
                let (u1, u2) = s.u_range();
                let (v1, v2) = s.v_range();
                (u1, u2, v1, v2)
            }
            None => (f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY),
        }
    }
}
