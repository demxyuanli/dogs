//! `BRepLib::UpdateDeflection` (`BRepLib.cxx:2793-2890`).
//!
//! Recomputes a face's `Poly_Triangulation` deflection from its own triangles:
//! the maximum distance, in squared form, between a sample of each triangle
//! (its barycenter, plus the midpoint of every interior link) expressed both in
//! UV and on the surface, and the resulting deflection is `sqrt` of that
//! (`cxx:2888`).
//!
//! The two auxiliary types are the reference's file-local ones:
//! `EvalDeflection` (`cxx:2730-2750`, the UV -> 3D residual) and `Link`
//! (`cxx:2752-2773` + `cxx:2776-2791`, an undirected node pair keyed by its
//! sorted indices).
//!
//! UNPORTED / deviations:
//! * `TopExp_Explorer(theShape, TopAbs_FACE)` and `(aFace, TopAbs_EDGE)`
//!   (`cxx:2795`, `:2816`) map to [`faces_of`] / [`edges_of`], which visit each
//!   distinct sub-shape once rather than each occurrence.

use std::collections::HashSet;

use occt_core::gp::GpPnt;
use occt_geom::Surface;

use crate::brep_tool::BRepTool;
use crate::shape::{Face, TopoShape};
use crate::topo_tools::{edges_of, faces_of};

/// `BRepLib::EvalDeflection` (`BRepLib.cxx:2730-2750`).
struct EvalDeflection {
    surface: std::sync::Arc<dyn Surface>,
}

impl EvalDeflection {
    /// `EvalDeflection(theFace)` (`cxx:2737-2741`): `BRepAdaptor_Surface`.
    fn new(face: &Face) -> Option<Self> {
        let surface = BRepTool::face_surface_world(face)?;
        Some(Self { surface })
    }

    /// `Eval(thePoint2d, thePoint3d)` (`cxx:2743-2749`): the *squared*
    /// distance between the surface point at the UV and the given 3D point.
    fn eval(&self, point2d: (f64, f64), point3d: &GpPnt) -> f64 {
        let on_surface = self.surface.d0(point2d.0, point2d.1);
        let dx = point3d.x() - on_surface.x();
        let dy = point3d.y() - on_surface.y();
        let dz = point3d.z() - on_surface.z();
        dx * dx + dy * dy + dz * dz
    }
}

/// `BRepLib::UpdateDeflection(theShape)` (`BRepLib.cxx:2793-2890`).
pub fn update_deflection(shape: &TopoShape) {
    for face in faces_of(shape) {
        update_face_deflection(&face);
    }
}

/// The loop body of `BRepLib::UpdateDeflection` for one face
/// (`BRepLib.cxx:2798-2889`).
fn update_face_deflection(face: &Face) {
    // `cxx:2799-2803`: no surface, nothing to measure against.
    let Some(tool) = EvalDeflection::new(face) else {
        return;
    };

    // `cxx:2805-2810`: `BRep_Tool::Triangulation(aFace, aLoc)`; both branches
    // skip the face when there is no triangulation or no UV nodes.
    let Some(triangulation) = BRepTool::face_triangulation(face) else {
        return;
    };
    let Some(uv_nodes) = triangulation.uv_nodes.as_ref() else {
        return;
    };

    // `cxx:2812-2835`: collect the nodes of every degenerative edge so the
    // elements built on them are skipped.
    let mut deg_nodes: HashSet<usize> = HashSet::new();
    for edge in edges_of(&face.0) {
        if !BRepTool::is_degenerated(&edge) {
            continue;
        }
        // `cxx:2822-2827`: a degenerated edge without a polygon on this
        // triangulation is skipped as a whole.
        let Some(polygon) = BRepTool::polygon_on_triangulation(
            &edge,
            &triangulation,
            edge.0.orientation().is_reversed(),
        ) else {
            continue;
        };
        for &node in &polygon.nodes {
            // `cxx:2830`: `aDegNodes.Add(aNodes(j))` - the polygon keeps
            // OCCT's 1-based node index, while this crate's
            // `Triangulation::triangles` are 0-based. A node index that was
            // never registered stays 0 and matches no triangle in OCCT either.
            if node >= 1 {
                deg_nodes.insert(node - 1);
            }
        }
    }

    let trsf = face.0.location().transformation();
    let mut links: HashSet<(usize, usize)> = HashSet::new();
    let mut sq_deflection = 0.0f64;
    for triangle in &triangulation.triangles {
        let nodes = [triangle.n0, triangle.n1, triangle.n2];
        if nodes.iter().any(|n| deg_nodes.contains(n)) {
            continue;
        }
        // `cxx:2853-2859`: world 3D nodes (the triangulation is stored in the
        // face frame) and the matching parametric nodes.
        let mut p3d = [GpPnt::zero(); 3];
        let mut p2d = [(0.0f64, 0.0f64); 3];
        let mut valid = true;
        for k in 0..3 {
            let Some(node) = triangulation.nodes.get(nodes[k]) else {
                valid = false;
                break;
            };
            let Some(uv) = uv_nodes.get(nodes[k]) else {
                valid = false;
                break;
            };
            p3d[k] = node.transformed(&trsf);
            p2d[k] = *uv;
        }
        if !valid {
            continue;
        }

        // `cxx:2861-2865`: barycenter of the triangle.
        let mid3d = GpPnt::new(
            (p3d[0].x() + p3d[1].x() + p3d[2].x()) / 3.0,
            (p3d[0].y() + p3d[1].y() + p3d[2].y()) / 3.0,
            (p3d[0].z() + p3d[1].z() + p3d[2].z()) / 3.0,
        );
        let mid2d = (
            (p2d[0].0 + p2d[1].0 + p2d[2].0) / 3.0,
            (p2d[0].1 + p2d[1].1 + p2d[2].1) / 3.0,
        );
        sq_deflection = sq_deflection.max(tool.eval(mid2d, &mid3d));

        // `cxx:2867-2885`: a link that is already in the set is an interior
        // link (its midpoint is not distorted by the boundary), so it is the
        // one that gets measured.
        for i in 0..3 {
            let j = (i + 1) % 3;
            let a = nodes[i];
            let b = nodes[j];
            let key = if a < b { (a, b) } else { (b, a) };
            if links.insert(key) {
                continue;
            }
            let mid3d_l = GpPnt::new(
                (p3d[i].x() + p3d[j].x()) / 2.0,
                (p3d[i].y() + p3d[j].y()) / 2.0,
                (p3d[i].z() + p3d[j].z()) / 2.0,
            );
            let mid2d_l = ((p2d[i].0 + p2d[j].0) / 2.0, (p2d[i].1 + p2d[j].1) / 2.0);
            sq_deflection = sq_deflection.max(tool.eval(mid2d_l, &mid3d_l));
        }
    }

    // `cxx:2888`: `aPT->Deflection(sqrt(aSqDeflection))`.
    //
    // `Poly_Triangulation::Deflection(...)` (`Poly_Triangulation.hxx:113`) is an
    // in-place member write on the object the face and the edge polygons share,
    // so [`Triangulation::set_deflection`] updates it without replacing the
    // handle (`BRep_PolygonOnTriangulation` keys its representations by it,
    // `BRep_PolygonOnTriangulation.cxx:36-49`).
    triangulation.set_deflection(sq_deflection.sqrt());
}
