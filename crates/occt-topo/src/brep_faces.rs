//! Phase 4 module: brep_faces — face utilities.
//!
//! Analytic queries and reconstructions on a `Face`: outer wire, surface
//! area, normal, world centroid, planarity, perimeter, a degree-1 B-spline
//! surface approximation, and UV sub-rectangle splitting.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpTrsf, GpVec, GpXyz};
use occt_geom::Surface;

use crate::brep_measure;
use crate::brep_surface;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Face, Wire};
use crate::topo_tools_full::{edge_vertices, edges_of_wire, vertex_position, wires_of_face};

/// Degree-1 tensor-product B-spline surface (bilinear patches), used by
/// [`make_surface_from_face`]. Poles are stored row-major (`n_u × n_v`).
#[derive(Clone)]
pub struct GeomBSplineSurface {
    pub poles: Vec<GpPnt>,
    pub n_u: usize,
    pub n_v: usize,
    pub knots_u: Vec<f64>,
    pub knots_v: Vec<f64>,
    pub degree_u: usize,
    pub degree_v: usize,
}

impl GeomBSplineSurface {
    pub fn new(
        poles: Vec<GpPnt>,
        n_u: usize,
        n_v: usize,
        knots_u: Vec<f64>,
        knots_v: Vec<f64>,
        degree_u: usize,
        degree_v: usize,
    ) -> Result<Self, String> {
        if poles.len() != n_u * n_v {
            return Err("GeomBSplineSurface: poles.len() must equal n_u * n_v".into());
        }
        Ok(Self { poles, n_u, n_v, knots_u, knots_v, degree_u, degree_v })
    }
}

impl Surface for GeomBSplineSurface {
    fn d0(&self, u: f64, v: f64) -> GpPnt {
        occt_core::bspl::surface::eval_surface(
            &self.poles,
            None,
            self.n_u,
            self.n_v,
            &self.knots_u,
            &self.knots_v,
            self.degree_u,
            self.degree_v,
            u,
            v,
        )
    }

    fn d1(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
        // ponytail: central finite differences — robust for degree-1 surfaces
        // (the BSplSLib D1 evaluator is still calibrated for degree >= 2).
        let h = 1e-6;
        let p0 = self.d0(u, v);
        let du = GpVec::from_pnts(&p0, &self.d0(u + h, v)).divided(h);
        let dv = GpVec::from_pnts(&p0, &self.d0(u, v + h)).divided(h);
        (p0, du, dv)
    }

    fn u_range(&self) -> (f64, f64) { (0.0, 1.0) }
    fn v_range(&self) -> (f64, f64) { (0.0, 1.0) }
    fn continuity(&self) -> u8 { 1 }
    fn transform(&mut self, t: &GpTrsf) {
        for p in &mut self.poles {
            *p = p.transformed(t);
        }
    }
    fn clone_dyn(&self) -> Box<dyn Surface> { Box::new(self.clone()) }
}

/// The first (outer) boundary wire of a face, if any.
pub fn face_outer_wire(face: &Face) -> Option<Wire> {
    wires_of_face(face).into_iter().next()
}

/// UV rectangle spanned by the face's outer-wire vertices projected onto its
/// world surface. Returns `None` for unbounded faces (no wire) or degenerate
/// projections.
fn face_uv_rect(face: &Face) -> Option<(f64, f64, f64, f64)> {
    let surf = BRepTool::face_surface_world(face)?;
    let wire = face_outer_wire(face)?;
    let mut umin = f64::INFINITY;
    let mut umax = f64::NEG_INFINITY;
    let mut vmin = f64::INFINITY;
    let mut vmax = f64::NEG_INFINITY;
    for e in edges_of_wire(&wire) {
        let (a, z) = edge_vertices(&e);
        for v in [a, z].into_iter().flatten() {
            let p = vertex_position(&v);
            let (u, vv) = brep_surface::surface_closest_params(surf.as_ref(), &p, 64, 64);
            umin = umin.min(u);
            umax = umax.max(u);
            vmin = vmin.min(vv);
            vmax = vmax.max(vv);
        }
    }
    if !(umin.is_finite() && umax.is_finite() && vmin.is_finite() && vmax.is_finite()) {
        return None;
    }
    if umax <= umin || vmax <= vmin {
        return None;
    }
    Some((umin, umax, vmin, vmax))
}

/// Grid-triangulated surface area of the face. Bounded surfaces delegate to
/// [`brep_measure::face_area`]; unbounded planar surfaces are triangulated over
/// the UV rectangle spanned by the outer wire instead of the infinite natural
/// surface.
pub fn face_area_exact(face: &Face, nu: usize, nv: usize) -> f64 {
    let surf = BRepTool::face_surface(face);
    let bounded = surf
        .as_ref()
        .map(|s| {
            let (u0, u1) = s.u_range();
            let (v0, v1) = s.v_range();
            u0.is_finite() && u1.is_finite() && v0.is_finite() && v1.is_finite()
        })
        .unwrap_or(false);
    if bounded {
        return brep_measure::face_area(face, nu, nv);
    }
    let (umin, umax, vmin, vmax) = match face_uv_rect(face) {
        Some(r) => r,
        None => return f64::INFINITY,
    };
    let Some(s) = &surf else { return f64::INFINITY };
    let nu = nu.max(2);
    let nv = nv.max(2);
    let mut area = 0.0;
    for i in 0..nu {
        for j in 0..nv {
            let p00 = s.d0(umin + (umax - umin) * i as f64 / nu as f64, vmin + (vmax - vmin) * j as f64 / nv as f64);
            let p10 = s.d0(umin + (umax - umin) * (i + 1) as f64 / nu as f64, vmin + (vmax - vmin) * j as f64 / nv as f64);
            let p01 = s.d0(umin + (umax - umin) * i as f64 / nu as f64, vmin + (vmax - vmin) * (j + 1) as f64 / nv as f64);
            let p11 = s.d0(umin + (umax - umin) * (i + 1) as f64 / nu as f64, vmin + (vmax - vmin) * (j + 1) as f64 / nv as f64);
            let tri = |a: &GpPnt, b: &GpPnt, c: &GpPnt| {
                0.5 * GpVec::from_pnts(a, b).xyz().crossed(&GpVec::from_pnts(a, c).xyz()).modulus()
            };
            area += tri(&p00, &p10, &p11) + tri(&p00, &p11, &p01);
        }
    }
    area
}

/// Unit surface normal of the face at its natural `(u, v)` parameters.
pub fn face_normal_at(face: &Face, u: f64, v: f64) -> Option<GpVec> {
    brep_surface::face_normal(face, u, v)
}

/// Centroid of the face region in world coordinates: the mean of a `16 × 16`
/// grid sampled over the outer-wire UV rectangle of the face's world surface.
pub fn face_centroid_world(face: &Face) -> Option<GpPnt> {
    let surf = BRepTool::face_surface_world(face)?;
    let (umin, umax, vmin, vmax) = face_uv_rect(face)?;
    let (nu, nv) = (16usize, 16usize);
    let mut acc = GpXyz::zero();
    let mut count = 0.0;
    for i in 0..nu {
        for j in 0..nv {
            let u = umin + (umax - umin) * (i as f64 + 0.5) / nu as f64;
            let v = vmin + (vmax - vmin) * (j as f64 + 0.5) / nv as f64;
            acc = acc.added(&surf.d0(u, v).coord);
            count += 1.0;
        }
    }
    if count == 0.0 {
        None
    } else {
        Some(GpPnt::from_xyz(&acc.divided(count)))
    }
}

/// Whether the face's surface is planar.
pub fn face_is_planar_face(face: &Face) -> bool {
    brep_surface::face_is_planar(face)
}

/// Total length of the face's outer-wire edges.
pub fn face_wire_length(face: &Face) -> f64 {
    face_outer_wire(face)
        .map(|w| brep_measure::wire_length(&w, 32))
        .unwrap_or(0.0)
}

/// Approximate the face surface with a degree-1 bilinear B-spline surface:
/// the sampled `nu × nv` grid of face-surface points becomes the pole grid.
pub fn make_surface_from_face(face: &Face, nu: usize, nv: usize) -> Result<GeomBSplineSurface, String> {
    let nu = nu.max(2);
    let nv = nv.max(2);
    let surf = BRepTool::face_surface_world(face).ok_or("make_surface_from_face: no face surface")?;
    let (umin, umax, vmin, vmax) =
        face_uv_rect(face).ok_or("make_surface_from_face: face has no bounded wire")?;
    let mut poles = Vec::with_capacity(nu * nv);
    for i in 0..nu {
        for j in 0..nv {
            let u = umin + (umax - umin) * i as f64 / (nu - 1) as f64;
            let v = vmin + (vmax - vmin) * j as f64 / (nv - 1) as f64;
            poles.push(surf.d0(u, v));
        }
    }
    let knots_u = occt_core::bspl::knots::build_uniform_knots(nu, 1);
    let knots_v = occt_core::bspl::knots::build_uniform_knots(nv, 1);
    GeomBSplineSurface::new(poles, nu, nv, knots_u, knots_v, 1, 1)
}

/// New face on the same surface, trimmed to the `(u0..u1, v0..v1)` sub-rectangle.
/// The boundary wire is four line segments between the four surface corners.
pub fn split_face_by_uv(face: &Face, u0: f64, u1: f64, v0: f64, v1: f64) -> Result<Face, String> {
    let surf = BRepTool::face_surface(face).ok_or("split_face_by_uv: no face surface")?;
    let p00 = surf.d0(u0, v0);
    let p10 = surf.d0(u1, v0);
    let p11 = surf.d0(u1, v1);
    let p01 = surf.d0(u0, v1);
    let b = TopoBuilder::new();
    let e0 = b.make_edge_segment(&p00, &p10);
    let e1 = b.make_edge_segment(&p10, &p11);
    let e2 = b.make_edge_segment(&p11, &p01);
    let e3 = b.make_edge_segment(&p01, &p00);
    let wire = b.make_wire(&[e0, e1, e2, e3]);
    wire.set_closed(true);
    Ok(b.make_face(surf, &[wire]))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep_builder_api::make_face_from_polygon;
    use crate::primitives::BRepPrimBox;
    use crate::topo_tools_full::{edges_of_wire, faces_of};
    use occt_core::gp::GpPnt;

    fn approx(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-6 * b.abs().max(1.0)
    }

    fn unit_box_top_face() -> Face {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        faces_of(&b.solid.0)
            .into_iter()
            .find(|f| BRepTool::face_surface(f).map(|s| (s.d0(0.0, 0.0).z() - 1.0).abs() < 1e-9).unwrap_or(false))
            .expect("top face")
    }

    #[test]
    fn outer_wire_and_length() {
        let top = unit_box_top_face();
        let w = face_outer_wire(&top).expect("outer wire");
        assert_eq!(edges_of_wire(&w).len(), 4);
        assert!(approx(face_wire_length(&top), 4.0));
    }

    #[test]
    fn area_and_normal() {
        let top = unit_box_top_face();
        assert!(approx(face_area_exact(&top, 16, 16), 1.0), "area {}", face_area_exact(&top, 16, 16));
        let n = face_normal_at(&top, 0.5, 0.5).expect("normal");
        assert!((n.xyz().modulus() - 1.0).abs() < 1e-9, "unit normal");
        assert!(n.z() > 0.9, "top face normal points +Z");
    }

    #[test]
    fn centroid_world_and_planarity() {
        let top = unit_box_top_face();
        let c = face_centroid_world(&top).expect("centroid");
        assert!((c.x() - 0.5).abs() < 1e-6, "cx {}", c.x());
        assert!((c.y() - 0.5).abs() < 1e-6, "cy {}", c.y());
        assert!((c.z() - 1.0).abs() < 1e-6, "cz {}", c.z());
        assert!(face_is_planar_face(&top));
    }

    #[test]
    fn surface_from_plane_face() {
        let sq = [
            GpPnt::new(0., 0., 0.),
            GpPnt::new(1., 0., 0.),
            GpPnt::new(1., 1., 0.),
            GpPnt::new(0., 1., 0.),
        ];
        let face = make_face_from_polygon(&sq).expect("square face");
        let s = make_surface_from_face(&face, 4, 4).expect("bspline surface");
        assert_eq!(s.n_u, 4);
        assert_eq!(s.n_v, 4);
        let p = s.d0(0.5, 0.5);
        assert!(p.distance(&GpPnt::new(0.5, 0.5, 0.0)) < 1e-6, "d0 {p:?}");
        let corner = s.d0(0.0, 0.0);
        assert!(corner.distance(&GpPnt::zero()) < 1e-6);
    }

    #[test]
    fn split_face_by_uv_creates_bounded_wire() {
        let sq = [
            GpPnt::new(0., 0., 0.),
            GpPnt::new(2., 0., 0.),
            GpPnt::new(2., 2., 0.),
            GpPnt::new(0., 2., 0.),
        ];
        let face = make_face_from_polygon(&sq).expect("square face");
        let split = split_face_by_uv(&face, 0.25, 0.75, 0.25, 0.75).expect("split");
        let w = face_outer_wire(&split).expect("outer wire");
        assert_eq!(edges_of_wire(&w).len(), 4);
        assert!(face_is_planar_face(&split));
    }
}
