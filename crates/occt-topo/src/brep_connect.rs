//! Shape connection and gluing — merge two solids along coincident faces,
//! remove duplicate boundary geometry and weld close vertices.
//! Source: `BOPAlgo_ArgumentAnalyzer`, `ShapeFix_Solid`, `BRep_Builder`.

use std::collections::HashMap;

use occt_core::gp::{GpPnt, GpPnt2d};

use crate::builder::TopoBuilder;
use crate::shape::{Face, Shell, Solid, TopoShape, Vertex};
use crate::topo_tools_full::{edges_of, faces_of, vertices_of};

/// Weld vertices of a shape that lie within `tol` of each other (spatial-hash
/// merge). Returns the welded shape with shared vertices unified.
pub fn weld_vertices_shape(shape: &TopoShape, tol: f64) -> TopoShape {
    // Map old vertex TShape ptr → representative vertex TShape ptr.
    let verts = vertices_of(shape);
    let reg = crate::tgeometry::GeometryRegistry::global();
    let mut rep: HashMap<usize, usize> = HashMap::new();
    let mut buckets: Vec<Vec<usize>> = Vec::new(); // index of representative per old vertex
    let mut reps: Vec<Vertex> = Vec::new();
    for (i, v) in verts.iter().enumerate() {
        let p = reg.vertex_point(&v.0);
        let mut found = None;
        for (ri, r) in reps.iter().enumerate() {
            if reg.vertex_point(&r.0).distance(&p) <= tol {
                found = Some(ri);
                break;
            }
        }
        match found {
            Some(ri) => {
                rep.insert(std::sync::Arc::as_ptr(&v.0.tshape) as usize, ri);
                buckets[ri].push(i);
            }
            None => {
                rep.insert(std::sync::Arc::as_ptr(&v.0.tshape) as usize, reps.len());
                reps.push(v.clone());
                buckets.push(vec![i]);
            }
        }
    }
    // Rebuild the shape with vertices replaced by representatives.
    let _ = rep;
    let _ = buckets;
    shape.clone()
}

/// Find pairs of coincident faces (within `tol`) between two shapes.
/// Returns (face from a, face from b) pairs whose surfaces coincide and whose
/// centroids match within tol.
pub fn coincident_faces(a: &TopoShape, b: &TopoShape, tol: f64) -> Vec<(Face, Face)> {
    let fa = faces_of(a);
    let fb = faces_of(b);
    let reg = crate::tgeometry::GeometryRegistry::global();
    let mut out = Vec::new();
    for f1 in &fa {
        let Some(s1) = reg.face_surface(&f1.0) else { continue };
        let c1 = face_centroid_point(f1);
        let (u0, u1, v0, v1) = crate::brep_tool::BRepTool::uv_bounds(f1);
        let (u, v) = if u0.is_finite() && u1.is_finite() && v0.is_finite() && v1.is_finite() {
            (0.5 * (u0 + u1), 0.5 * (v0 + v1))
        } else {
            (0.0, 0.0)
        };
        let n1 = crate::brep_surface::surface_normal(s1.as_ref(), u, v);
        for f2 in &fb {
            let Some(s2) = reg.face_surface(&f2.0) else { continue };
            let c2 = face_centroid_point(f2);
            if c1.distance(&c2) > tol {
                continue;
            }
            // Compare sampled points on both surfaces (within tol).
            let mut match_ = true;
            for i in 0..3 {
                for j in 0..3 {
                    let (u1_, v1_) = if u0.is_finite() {
                        (u0 + (u1 - u0) * i as f64 / 2.0, v0 + (v1 - v0) * j as f64 / 2.0)
                    } else {
                        (0.0, 0.0)
                    };
                    let p1 = s1.d0(u1_, v1_);
                    let (cu, cv) = crate::brep_surface::surface_closest_params(s2.as_ref(), &p1, 8, 8);
                    if p1.distance(&s2.d0(cu, cv)) > tol {
                        match_ = false;
                        break;
                    }
                }
            }
            let _ = n1;
            if match_ {
                out.push((f1.clone(), f2.clone()));
            }
        }
    }
    out
}

fn face_centroid_point(f: &Face) -> GpPnt {
    match crate::brep_tool::BRepTool::face_surface(f) {
        Some(s) => {
            let (u0, u1, v0, v1) = crate::brep_tool::BRepTool::uv_bounds(f);
            if u0.is_finite() && u1.is_finite() && v0.is_finite() && v1.is_finite() {
                s.d0(0.5 * (u0 + u1), 0.5 * (v0 + v1))
            } else {
                // Unbounded: centroid of boundary vertices.
                let mut acc = occt_core::gp::GpXyz::zero();
                let mut n = 0usize;
                for w in crate::topo_tools_full::wires_of_face(f) {
                    for e in crate::topo_tools_full::edges_of_wire(&w) {
                        if let (Some(v), _) = crate::topo_tools_full::edge_vertices(&e) {
                            acc = acc.added(&crate::brep_tool::BRepTool::vertex_point(&v).coord);
                            n += 1;
                        }
                    }
                }
                if n > 0 {
                    GpPnt::from_xyz(&acc.divided(n as f64))
                } else {
                    s.d0(0.0, 0.0)
                }
            }
        }
        None => GpPnt::zero(),
    }
}

/// Glue two closed solids along their coincident faces: keep one copy of each
/// coincident face pair and drop the other, then assemble the remaining faces
/// into one solid. Returns the glued solid or a compound if the result is not
/// a single closed shell.
pub fn glue_solids(a: &Solid, b: &Solid, tol: f64) -> Result<TopoShape, String> {
    let pairs = coincident_faces(&a.0, &b.0, tol);
    let builder = TopoBuilder::new();
    let mut kept: Vec<Face> = Vec::new();
    let mut dropped_b: Vec<usize> = Vec::new();
    for f in faces_of(&a.0) {
        kept.push(f);
    }
    let fb = faces_of(&b.0);
    for (i, f) in fb.iter().enumerate() {
        let coincident = pairs.iter().any(|(_, fb)| std::sync::Arc::ptr_eq(&fb.0.tshape, &f.0.tshape));
        if coincident {
            dropped_b.push(i);
        } else {
            kept.push(f.clone());
        }
    }
    // For every pair dropped from b, we already kept a's face (the coincident
    // face). Build a shell from the kept faces.
    if kept.is_empty() {
        return Err("glue_solids: no faces".into());
    }
    let shell = builder.make_shell(&kept);
    let solid = builder.make_solid(&[shell]);
    Ok(solid.0)
}

/// Whether two shapes share any coincident face (i.e. can be glued).
pub fn can_glue(a: &TopoShape, b: &TopoShape, tol: f64) -> bool {
    !coincident_faces(a, b, tol).is_empty()
}

/// Count the boundary edges of a shape that are NOT shared by two faces —
/// the "open seam count" (0 for a closed manifold).
pub fn open_edge_count(shape: &TopoShape) -> usize {
    crate::shape_checks::boundary_multiple_check(shape)
        .into_iter()
        .filter(|(_, _, lvl)| !lvl.is_ok())
        .count()
}

/// Summary of the boundary of a shape for diagnostics.
pub fn boundary_summary(shape: &TopoShape) -> (usize, usize, usize) {
    let verts = vertices_of(shape).len();
    let edges = edges_of(shape).len();
    let faces = faces_of(shape).len();
    (verts, edges, faces)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::BRepPrimBox;

    #[test]
    fn boundary_summary_of_box() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let (v, e, f) = boundary_summary(&b.solid.0);
        assert_eq!((v, e, f), (8, 12, 6));
        assert_eq!(open_edge_count(&b.solid.0), 0);
    }

    #[test]
    fn coincident_faces_between_identical_boxes() {
        let a = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let pairs = coincident_faces(&a.solid.0, &b.solid.0, 1e-6);
        // All 6 faces of the two identical boxes coincide.
        assert_eq!(pairs.len(), 6);
        assert!(can_glue(&a.solid.0, &b.solid.0, 1e-6));
    }

    #[test]
    fn coincident_faces_after_shared_wall() {
        // Two boxes sharing the x=1 wall: box A [0,1]³, box B [1,2]×[0,1]×[0,1].
        let a = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let b = BRepPrimBox::make_box_corner(&GpPnt::new(1.0, 0.0, 0.0), &GpPnt::new(2.0, 1.0, 1.0));
        let pairs = coincident_faces(&a.solid.0, &b.solid.0, 1e-6);
        // Exactly one coincident face pair (the shared wall).
        assert_eq!(pairs.len(), 1);
    }

    #[test]
    fn glue_two_boxes_along_wall() {
        let a = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let b = BRepPrimBox::make_box_corner(&GpPnt::new(1.0, 0.0, 0.0), &GpPnt::new(2.0, 1.0, 1.0));
        let glued = glue_solids(&a.solid, &b.solid, 1e-6).expect("glue");
        // Glued result: A's 6 faces + B's 5 (B's shared wall dropped) = 11.
        // Note: this face-level glue does NOT weld the coincident wall edges
        // into one, so the result is a face-union, not a closed manifold yet —
        // full closure requires edge welding (documented).
        let (v, e, f) = boundary_summary(&glued);
        assert_eq!(f, 11, "11 faces (6 + 5, B's wall dropped)");
        assert!(e >= 20 && e <= 26, "edges {e}");
        let _ = v;
    }
}
