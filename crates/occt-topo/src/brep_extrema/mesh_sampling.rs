use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_geom::Surface;

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::shape::{Face, TopoShape};
use crate::topexp::Explorer;

use super::*;

// ---------------------------------------------------------------------------
// Fallback shape → face-mesh sampling.
//
// `crate::shape_mesh::mesh_shape` is being written concurrently; until it
// lands, brep_gprop and is_inside use this local grid triangulator. Each face
// is sampled on a regular `(nu+1)×(nv+1)` grid in parametric space. Faces
// with an unbounded parametric range (e.g. planes) get a window inferred from
// their boundary vertices projected onto the surface.
// ---------------------------------------------------------------------------

/// Triangulate every face of `shape` into `(vertices, triangles)` pairs.
///
/// Vertices are `GpPnt`s; triangles are index triples into the corresponding
/// vertex list. Faces without a registered surface (or without any boundary
/// from which to infer an unbounded UV window) are skipped.
pub fn mesh_faces(
    shape: &TopoShape,
    nu: usize,
    nv: usize,
) -> Vec<(Vec<GpPnt>, Vec<(usize, usize, usize)>)> {
    let mut out = Vec::new();
    let mut ex = Explorer::new(shape, ShapeType::Face);
    while ex.more() {
        let fs = ex.current().clone();
        ex.next();
        if let Some(face) = Face::wrap(fs) {
            if let Some(surface) = BRepTool::face_surface(&face) {
                if let Some((u1, u2, v1, v2)) = face_uv_window(&face) {
                    let (verts, mut tris) = sample_grid(&surface, u1, u2, v1, v2, nu, nv);
                    // The grid follows the surface parameterization; a REVERSED
                    // face has its material on the opposite side, so its
                    // triangles must wind the other way (OCCT `BRep_Tool::
                    // Triangulation` + the face orientation the caller applies).
                    if face.orientation().is_reversed() {
                        for t in &mut tris {
                            std::mem::swap(&mut t.1, &mut t.2);
                        }
                    }
                    out.push((verts, tris));
                }
            }
        }
    }
    out
}

/// Sample the surface on a regular parametric grid.
fn sample_grid(
    surface: &Arc<dyn Surface>,
    u1: f64,
    u2: f64,
    v1: f64,
    v2: f64,
    nu: usize,
    nv: usize,
) -> (Vec<GpPnt>, Vec<(usize, usize, usize)>) {
    let (nu, nv) = (nu.max(1), nv.max(1));
    let mut verts = Vec::with_capacity((nu + 1) * (nv + 1));
    for i in 0..=nu {
        let u = u1 + (u2 - u1) * (i as f64 / nu as f64);
        for j in 0..=nv {
            let v = v1 + (v2 - v1) * (j as f64 / nv as f64);
            verts.push(surface.d0(u, v));
        }
    }
    let stride = nv + 1;
    let mut tris = Vec::with_capacity(2 * nu * nv);
    for i in 0..nu {
        for j in 0..nv {
            let a = i * stride + j;
            let b = a + 1;
            let c = a + stride;
            let d = c + 1;
            tris.push((a, b, c));
            tris.push((b, d, c));
        }
    }
    (verts, tris)
}

