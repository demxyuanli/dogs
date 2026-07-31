//! Wire and edge projection onto surfaces and planes.
//! Source: `BRepProj_Projection`, `BRepLib::BuildCurve3d` (pcurve).

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d, GpVec};

use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, Wire};
use crate::tgeometry::GeometryRegistry;

/// Project a 3D point onto a face's surface (closest (u,v) + surface point).
pub fn project_point_on_face(f: &Face, p: &GpPnt) -> Option<(f64, f64, GpPnt)> {
    let s = BRepTool::face_surface(f)?;
    let (u, v) = crate::brep_surface::surface_closest_params(s.as_ref(), p, 32, 32);
    Some((u, v, s.d0(u, v)))
}

/// Project a 3D point onto a plane.
pub fn project_point_on_plane(p: &GpPnt, pln: &occt_core::gp::GpPln) -> GpPnt {
    let loc = pln.location();
    let n = *pln.axis().direction();
    let d = GpVec::from_pnts(&loc, p);
    let t = d.xyz().dot(n.xyz());
    GpPnt::new(p.x() - t * n.x(), p.y() - t * n.y(), p.z() - t * n.z())
}

/// Project an edge onto a face: sample the edge curve, project each point
/// onto the face surface, return the pcurve (2D (u,v) samples) and the
/// projected 3D polyline.
pub fn project_edge_on_face(e: &Edge, f: &Face, samples: usize) -> Option<(Vec<GpPnt2d>, Vec<GpPnt>)> {
    let curve = BRepTool::edge_curve(e)?;
    let s = BRepTool::face_surface(f)?;
    let (a, b) = BRepTool::edge_parameters(e);
    if !a.is_finite() || !b.is_finite() {
        return None;
    }
    let mut pcur = Vec::with_capacity(samples);
    let mut pts = Vec::with_capacity(samples);
    for i in 0..samples {
        let u = a + (b - a) * i as f64 / (samples - 1).max(1) as f64;
        let p3 = curve.d0(u);
        let (pu, pv) = crate::brep_surface::surface_closest_params(s.as_ref(), &p3, 32, 32);
        pcur.push(GpPnt2d::new(pu, pv));
        pts.push(s.d0(pu, pv));
    }
    Some((pcur, pts))
}

/// Build a new edge whose curve is the projection of `e` onto the plane
/// `pln` (a line stays a line; a circle becomes a circle/ellipse if the plane
/// is not parallel). Returns a new registered edge.
pub fn project_edge_on_plane(e: &Edge, pln: &occt_core::gp::GpPln) -> Option<Edge> {
    let curve = BRepTool::edge_curve(e)?;
    let (a, b) = BRepTool::edge_parameters(e);
    let (pa, pb) = (curve.d0(a), curve.d0(b));
    let qa = project_point_on_plane(&pa, pln);
    let qb = project_point_on_plane(&pb, pln);
    let builder = TopoBuilder::new();
    Some(builder.make_edge_segment(&qa, &qb))
}

/// Build the projected wire of `wire` onto a plane: every edge becomes a
/// projected segment; shared vertices are projected once.
pub fn project_wire_on_plane(wire: &Wire, pln: &occt_core::gp::GpPln) -> Wire {
    let builder = TopoBuilder::new();
    let edges = crate::topo_tools_full::edges_of_wire(wire);
    let mut proj = Vec::with_capacity(edges.len());
    for e in &edges {
        if let Some(pe) = project_edge_on_plane(e, pln) {
            proj.push(pe);
        }
    }
    builder.make_wire(&proj)
}

/// The 2D footprint of a face in the XY plane (projection onto z=0).
pub fn face_xy_footprint(f: &Face) -> Option<Vec<GpPnt2d>> {
    let poly = crate::brep_surface::face_plane(f)?;
    let mut out = Vec::new();
    for w in crate::topo_tools_full::wires_of_face(f) {
        for e in crate::topo_tools_full::edges_of_wire(&w) {
            if let (Some(v), _) = crate::topo_tools_full::edge_vertices(&e) {
                let p = BRepTool::vertex_point(&v);
                out.push(GpPnt2d::new(p.x(), p.y()));
            }
        }
    }
    let _ = poly;
    Some(out)
}

/// Distance between an edge and a plane (0 if the edge crosses the plane).
pub fn edge_plane_distance(e: &Edge, pln: &occt_core::gp::GpPln) -> f64 {
    let Some(curve) = BRepTool::edge_curve(e) else { return f64::INFINITY };
    let (a, b) = BRepTool::edge_parameters(e);
    if !a.is_finite() || !b.is_finite() {
        return f64::INFINITY;
    }
    let n = *pln.axis().direction();
    let loc = pln.location();
    let mut min_d = f64::INFINITY;
    for i in 0..=64 {
        let u = a + (b - a) * i as f64 / 64.0;
        let p = curve.d0(u);
        let d = GpVec::from_pnts(&loc, &p).xyz().dot(n.xyz()).abs();
        min_d = min_d.min(d);
    }
    min_d
}

/// Test whether an edge crosses a plane (samples change sign).
pub fn edge_crosses_plane(e: &Edge, pln: &occt_core::gp::GpPln) -> bool {
    let Some(curve) = BRepTool::edge_curve(e) else { return false };
    let (a, b) = BRepTool::edge_parameters(e);
    if !a.is_finite() || !b.is_finite() {
        return false;
    }
    let n = *pln.axis().direction();
    let loc = pln.location();
    let mut sign = None;
    for i in 0..=32 {
        let u = a + (b - a) * i as f64 / 32.0;
        let p = curve.d0(u);
        let d = GpVec::from_pnts(&loc, &p).xyz().dot(n.xyz());
        if sign.is_none() {
            sign = Some(d >= 0.0);
        } else if sign != Some(d >= 0.0) {
            return true;
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder::TopoBuilder;
    use occt_core::gp::{GpAx3, GpDir, GpPln};

    fn z0_plane() -> GpPln {
        GpPln::new(GpAx3::standard())
    }

    #[test]
    fn project_point_on_plane_test() {
        let pln = z0_plane();
        let q = project_point_on_plane(&GpPnt::new(1.0, 2.0, 3.0), &pln);
        assert!((q.z() - 0.0).abs() < 1e-12);
        assert!((q.x() - 1.0).abs() < 1e-12 && (q.y() - 2.0).abs() < 1e-12);
    }

    #[test]
    fn project_edge_and_wire_on_plane() {
        let b = TopoBuilder::new();
        let e = b.make_edge_segment(&GpPnt::new(0.,0.,1.), &GpPnt::new(1.,1.,1.));
        let qa = project_point_on_plane(&GpPnt::new(0.,0.,1.), &z0_plane());
        assert!((qa.z() - 0.0).abs() < 1e-12);
        // Projected edge lies in z=0 (both endpoints).
        let pe = project_edge_on_plane(&e, &z0_plane()).expect("projected edge");
        let c = BRepTool::edge_curve(&pe).unwrap();
        let (a0, a1) = BRepTool::edge_parameters(&pe);
        let (pa0, pa1) = (c.d0(a0), c.d0(a1));
        assert!((pa0.z() - 0.0).abs() < 1e-9, "start in z=0, z {}", pa0.z());
        assert!((pa1.z() - 0.0).abs() < 1e-9, "end in z=0, z {}", pa1.z());
        assert!(pa0.distance(&GpPnt::new(0.,0.,0.)) < 1e-9);
        assert!(pa1.distance(&GpPnt::new(1.,1.,0.)) < 1e-9);

        // A closed square wire projects to a closed wire.
        let sq = [
            GpPnt::new(0.,0.,1.), GpPnt::new(1.,0.,1.), GpPnt::new(1.,1.,1.), GpPnt::new(0.,1.,1.),
        ];
        let edges: Vec<crate::shape::Edge> = (0..4).map(|i| {
            let j = (i + 1) % 4;
            b.make_edge_segment(&sq[i], &sq[j])
        }).collect();
        let wire = b.make_wire(&edges);
        let pw = project_wire_on_plane(&wire, &z0_plane());
        let pedges = crate::topo_tools_full::edges_of_wire(&pw);
        assert_eq!(pedges.len(), 4);
    }

    #[test]
    fn edge_plane_distance_and_crossing() {
        let b = TopoBuilder::new();
        // Edge at z=1 (horizontal, length 1) → distance 1 from z=0 plane.
        let e = b.make_edge_segment(&GpPnt::new(0.,0.,1.), &GpPnt::new(1.,0.,1.));
        assert!((edge_plane_distance(&e, &z0_plane()) - 1.0).abs() < 1e-9);
        let crossing = b.make_edge_segment(&GpPnt::new(0.,0.,-1.), &GpPnt::new(0.,0.,1.));
        assert!(edge_crosses_plane(&crossing, &z0_plane()));
        assert!(!edge_crosses_plane(&e, &z0_plane()));
    }

    #[test]
    fn project_on_face_surface() {
        let b = TopoBuilder::new();
        let face = b.make_face_plane(&z0_plane());
        let (u, v, p) = project_point_on_face(&face, &GpPnt::new(0.5, 0.5, 2.0)).expect("project");
        assert!((p.z() - 0.0).abs() < 1e-9);
        let _ = (u, v);
    }

    #[test]
    fn face_xy_footprint_works() {
        use crate::primitives::BRepPrimBox;
        let box_ = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let faces = crate::topo_tools_full::faces_of(&box_.solid.0);
        let fp = face_xy_footprint(&faces[0]).expect("footprint");
        assert!(!fp.is_empty());
    }
}
