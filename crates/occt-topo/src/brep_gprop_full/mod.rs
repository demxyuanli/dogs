//! Exact BRepGProp mass properties — full port of `TKTopAlgo/BRepGProp`.
//!
//! Computes linear (edge-length), surface (area) and volume global properties
//! of a `TopoShape` by integrating over the exact boundary geometry, mirroring
//! the OCCT `BRepGProp` entry points:
//!
//! * [`linear_properties`] — Gauss quadrature along every edge's curve
//!   (`BRepGProp_Cinert` + `BRepGProp_EdgeTool`). The edge curve is a
//!   `dyn Curve`; d0/d1 are used directly, d2 numerically when needed.
//! * [`surface_properties`] — Gauss quadrature over each face's UV domain
//!   (`BRepGProp_Sinert` + `BRepGProp_Gauss`). Returns the total area.
//! * [`volume_properties`] — the divergence-theorem surface integral over each
//!   face (`BRepGProp_Vinert` + `BRepGProp_Gauss`), giving the signed volume.
//! * [`volume_properties_gk`] — the same volume via adaptive Gauss–Kronrod
//!   (`BRepGProp_VinertGK` + `UFunction` / `TFunction`).
//!
//! Trimmed faces (bounded by wires) are integrated with the boundary line
//! integral (Green's theorem) over the pcurves of the wire edges, exactly like
//! OCCT. The face's UV bounds come from the pcurves (`BRepTools::UVBounds`),
//! which makes the integration well-defined even for surfaces (planes,
//! cylinders) whose natural parameter ranges are unbounded. The wire's UV
//! orientation sign is recovered from the signed area of the UV boundary
//! polygon, so the direct (positively-oriented) surface integral is obtained
//! regardless of how the wires were wound. Faces whose wire repeats an edge
//! (the cylinder lateral face's double seam) are integrated directly over the
//! UV bounding rectangle.
//!
//! The surface `d1` derivative is computed numerically for surfaces whose port
//! `d1` returns zero vectors (cylinder / sphere / torus / cone), following the
//! Phase 13 convention that `Surface` carries no `d2`.
mod prelude {

pub(crate) use std::f64::consts::PI;
pub(crate) use std::sync::Arc;

pub(crate) use occt_core::gp::{GpMat, GpPnt, GpPnt2d, GpVec, GpVec2d, GpXyz};
pub(crate) use occt_geom::{Curve, Surface};
pub(crate) use occt_math::gauss::gauss_legendre;

pub(crate) use crate::brep_surface::classify_surface;
pub(crate) use crate::brep_tool::BRepTool;
pub(crate) use crate::shape::{Edge, Face, TopoShape};
pub(crate) use crate::topo_tools_full::{edges_of, edges_of_wire, faces_of, vertices_of, wires_of_face};

}

use prelude::*;


// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::{BRepPrimBox, BRepPrimCylinder, BRepPrimSphere};

    const TOL: f64 = 1e-5;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn box_surface_volume() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let (props, area) = surface_properties(&b.solid.0).unwrap();
        assert!(approx(area, 6.0, 1e-9), "area {area}");
        assert!(approx(props.dim, 6.0, 1e-9));
        let c = props.center();
        assert!(c.distance(&GpPnt::new(0.5, 0.5, 0.5)) < 1e-9, "center {:?}", c);

        let v = volume_properties(&b.solid.0).unwrap();
        assert!(approx(v.dim, 1.0, 1e-9), "volume {}", v.dim);
        let cv = v.center();
        assert!(cv.distance(&GpPnt::new(0.5, 0.5, 0.5)) < 1e-9, "volume center {:?}", cv);
    }

    #[test]
    fn box_surface_volume_2x3x4() {
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        let (props, area) = surface_properties(&b.solid.0).unwrap();
        assert!(approx(area, 52.0, 1e-9), "area {area}");
        let _ = props;
        let v = volume_properties(&b.solid.0).unwrap();
        assert!(approx(v.dim, 24.0, 1e-9), "volume {}", v.dim);
        assert!(v.center().distance(&GpPnt::new(1.0, 1.5, 2.0)) < 1e-9);
    }

    #[test]
    fn box_matches_analytic() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let (_, area) = surface_properties(&b.solid.0).unwrap();
        let v = volume_properties(&b.solid.0).unwrap();
        // gprop_analytic reports 6 and 1 for the unit box.
        assert!(approx(area, crate::gprop_analytic::analytic_surface_area(&b.solid.0).unwrap(), 1e-9));
        assert!(approx(v.dim, crate::gprop_analytic::analytic_volume(&b.solid.0).unwrap(), 1e-9));
    }

    #[test]
    fn sphere_surface_volume() {
        let s = BRepPrimSphere::make_sphere(1.0);
        let (props, area) = surface_properties(&s.solid.0).unwrap();
        assert!(approx(area, 4.0 * PI, 1e-6), "area {area}");
        let _ = props;
        let v = volume_properties(&s.solid.0).unwrap();
        assert!(approx(v.dim, 4.0 / 3.0 * PI, 1e-6), "volume {}", v.dim);
    }

    #[test]
    fn cylinder_surface_volume() {
        let c = BRepPrimCylinder::make_cylinder(1.0, 2.0);
        let (props, area) = surface_properties(&c.solid.0).unwrap();
        assert!(approx(area, 6.0 * PI, 1e-6), "area {area} (want {})", 6.0 * PI);
        let _ = props;
        let v = volume_properties(&c.solid.0).unwrap();
        assert!(approx(v.dim, 2.0 * PI, 1e-6), "volume {} (want {})", v.dim, 2.0 * PI);
    }

    #[test]
    fn linear_properties_box() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let l = linear_properties(&b.solid.0).unwrap();
        // 12 distinct edges, each of length 1.
        assert!(approx(l.dim, 12.0, 1e-9), "edge length {}", l.dim);
        assert!(l.center().distance(&GpPnt::new(0.5, 0.5, 0.5)) < 1e-9);
    }

    #[test]
    fn adaptive_box_sphere() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let (props, err) = surface_properties_adaptive(&b.solid.0, 1e-4).unwrap();
        assert!(approx(props.dim, 6.0, 1e-4), "adaptive box area {}", props.dim);
        assert!(err >= 0.0);
        let (vprops, _) = volume_properties_adaptive(&b.solid.0, 1e-4).unwrap();
        assert!(approx(vprops.dim, 1.0, 1e-4), "adaptive box volume {}", vprops.dim);

        let s = BRepPrimSphere::make_sphere(1.0);
        let (sprops, _) = surface_properties_adaptive(&s.solid.0, 1e-3).unwrap();
        assert!(approx(sprops.dim, 4.0 * PI, 1e-2), "adaptive sphere area {}", sprops.dim);
        let (svprops, _) = volume_properties_adaptive(&s.solid.0, 1e-3).unwrap();
        assert!(approx(svprops.dim, 4.0 / 3.0 * PI, 1e-2), "adaptive sphere volume {}", svprops.dim);
    }

    #[test]
    fn volume_properties_gk_box() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let v = volume_properties_gk(&b.solid.0).unwrap();
        assert!(approx(v.dim, 1.0, 1e-3), "GK volume {}", v.dim);
    }

    #[test]
    fn mesh_props_box_matches() {
        use crate::brep_extrema::mesh_faces;
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let mut nodes = Vec::new();
        let mut tris = Vec::new();
        for (fv, ft) in mesh_faces(&b.solid.0, 8, 8) {
            let base = nodes.len();
            nodes.extend(fv);
            for (a, b2, c) in ft {
                tris.push(occt_core::poly::triangulation::Triangle::new(base + a, base + b2, base + c));
            }
        }
        let mesh = occt_core::poly::triangulation::Triangulation::new(nodes, tris);
        let props = mesh_props(&mesh, &GpPnt::zero(), false, true);
        // The mesh orientation from `mesh_faces` may be inward; the signed
        // volume is taken in absolute value (as `brep_gprop` does).
        assert!(approx(props.dim.abs(), 1.0, 1e-2), "mesh volume {}", props.dim);
    }

    #[test]
    fn curve_tool_orders() {
        let b = crate::builder::TopoBuilder::new();
        let e = b.make_edge_segment(&GpPnt::zero(), &GpPnt::new(3.0, 4.0, 0.0));
        let c = BRepTool::edge_curve(&e).unwrap();
        let (a, bb) = BRepTool::edge_parameters(&e);
        assert_eq!(curve_integration_order(c.as_ref(), a, bb), 2);
    }
}

mod p01;
mod p02;
mod p03;
mod p04;
pub use p01::*;
pub use p02::*;
pub use p03::*;
pub use p04::*;
