//! Incremental deflection mesh — a port of `BRepMesh_IncrementalMesh`.
//!
//! Meshes every face of a shape with a deflection-bounded adaptive
//! tessellation: planar faces are triangulated exactly from their boundary
//! polygon; curved faces use recursive UV subdivision that refines each cell
//! until the surface's midpoint deviation from the cell's bilinear patch falls
//! below the requested deflection (or a fixed depth cap is reached). Edge
//! polylines are available through `wireframe::edge_to_polyline` /
//! `curve_approx::curve_to_polyline`; this module does not re-emit them into the
//! triangle soup (the face tessellation is self-contained).

use std::collections::HashMap;

use occt_core::gp::{GpPnt, GpXyz};
use occt_core::poly::triangulation::Triangle;
use occt_geom::Curve;

use crate::brep_surface::{face_is_planar, face_plane, surface_closest_params, surface_normal};
use crate::brep_tool::BRepTool;
use crate::mesh::ShapeMesh;
use crate::shape::{Face, TopoShape};
use crate::topo_tools_full::{edges_of_wire, faces_of, wires_of_face};
use crate::wireframe::{edge_to_polyline, face_to_triangles};

/// Result of an incremental meshing pass.
#[derive(Debug, Clone)]
pub struct IncrementalMesh {
    pub mesh: ShapeMesh,
    /// Vertex count contributed by each face (in face traversal order).
    pub nodes_per_face: Vec<usize>,
    /// Deflection actually applied (clamped to a positive minimum).
    pub deflection_used: f64,
    /// Number of cell subdivisions performed while meshing curved faces.
    pub iterations: usize,
}

/// Max recursion depth for curved-face cell subdivision.
const MAX_FACE_DEPTH: usize = 9;
/// Max recursion depth for `adaptive_curve_polyline`.
const MAX_CURVE_DEPTH: usize = 24;

/// Mesh a whole shape with a deflection-based adaptive tessellation.
///
/// Planar faces are triangulated exactly from their boundary polygon (a box
/// face becomes two triangles); curved faces use the adaptive UV-refinement
/// path. Faces without registrable geometry yield an empty contribution.
pub fn incremental_mesh(shape: &TopoShape, deflection: f64) -> Result<IncrementalMesh, String> {
    let def = deflection.max(1e-9);
    let mut vertices: Vec<GpPnt> = Vec::new();
    let mut triangles: Vec<Triangle> = Vec::new();
    let mut nodes_per_face: Vec<usize> = Vec::new();
    let mut iterations = 0usize;

    for face in faces_of(shape) {
        let (vs, ts) = mesh_face(&face, def, &mut iterations);
        nodes_per_face.push(vs.len());
        let offset = vertices.len();
        vertices.extend(vs);
        for t in ts {
            triangles.push(Triangle::new(offset + t.n0, offset + t.n1, offset + t.n2));
        }
    }

    Ok(IncrementalMesh {
        mesh: ShapeMesh {
            vertices,
            triangles,
            source_shape: shape.shape_type(),
        },
        nodes_per_face,
        deflection_used: def,
        iterations,
    })
}

/// Mesh one face: planar → exact boundary-polygon triangulation, otherwise the
/// adaptive UV subdivision.
fn mesh_face(face: &Face, deflection: f64, iterations: &mut usize) -> (Vec<GpPnt>, Vec<Triangle>) {
    if face_is_planar(face) {
        return planar_face_mesh(face, deflection);
    }
    adaptive_face_mesh(face, deflection, iterations)
}

/// Exact triangulation of a planar face from its boundary polygon.
///
/// The boundary is collected by discretizing each wire edge with
/// [`edge_to_polyline`] (edge curves may be stored in either direction, so the
/// collected points are first deduplicated and then ordered by angle around the
/// face centroid, which recovers the boundary order of a convex polygon). If
/// the resulting polygon is convex it is fan-triangulated with outward
/// orientation; otherwise we fall back to the uniform grid of
/// [`face_to_triangles`].
fn planar_face_mesh(face: &Face, deflection: f64) -> (Vec<GpPnt>, Vec<Triangle>) {
    // Distinct boundary points. The boundary is sampled with the given
    // deflection (not a hard-coded 1e-6 — that would over-plate a circular disk
    // cap into thousands of segments instead of the ~tens OCCT uses).
    let bound_def = deflection.max(0.01);
    let mut pts: Vec<GpPnt> = Vec::new();
    for wire in wires_of_face(face) {
        for e in edges_of_wire(&wire) {
            for p in edge_to_polyline(&e, bound_def) {
                if !pts.iter().any(|q| q.distance(&p) < 1e-9) {
                    pts.push(p);
                }
            }
        }
    }
    if pts.len() < 3 {
        return face_to_triangles(face, 1e-6);
    }

    let pln = face_plane(face).unwrap_or_else(occt_core::gp::GpPln::default);
    let n = *pln.axis().direction().xyz();
    let xdir = *pln.x_axis().direction().xyz();
    let ydir = xdir.crossed(&n);
    let centroid = pts
        .iter()
        .fold(GpPnt::new(0.0, 0.0, 0.0), |acc, p| {
            GpPnt::new(acc.x() + p.x(), acc.y() + p.y(), acc.z() + p.z())
        });
    let inv = 1.0 / pts.len() as f64;
    let c = GpPnt::new(centroid.x() * inv, centroid.y() * inv, centroid.z() * inv);

    // Order the boundary points around the centroid (radial sweep).
    pts.sort_by(|p, q| {
        let ap = p.coord.subtracted(&c.coord);
        let aq = q.coord.subtracted(&c.coord);
        let ang_p = ap.dot(&ydir).atan2(ap.dot(&xdir));
        let ang_q = aq.dot(&ydir).atan2(aq.dot(&xdir));
        ang_p.partial_cmp(&ang_q).unwrap_or(std::cmp::Ordering::Equal)
    });

    if !is_convex(&pts, &n) {
        return face_to_triangles(face, 1e-6);
    }

    let mut tris = Vec::with_capacity(pts.len() - 2);
    for i in 1..pts.len() - 1 {
        tris.push(orient3(&pts, 0, i, i + 1, &n));
    }
    (pts, tris)
}

/// Adaptive curved-face meshing: recursive UV refinement.
///
/// A cell is split into four when the surface point at its UV midpoint deviates
/// from the bilinear interpolant of the four corners by more than `deflection`
/// and the depth cap has not been reached. Final cells emit two triangles
/// wound so their normal agrees with the surface normal at the cell midpoint.
fn adaptive_face_mesh(
    face: &Face,
    deflection: f64,
    iterations: &mut usize,
) -> (Vec<GpPnt>, Vec<Triangle>) {
    let Some(surface) = BRepTool::face_surface(face) else {
        return face_to_triangles(face, deflection);
    };
    let (u0, u1, v0, v1) = crate::wireframe::face_uv_bounds(face, surface.as_ref());
    if !(u0.is_finite() && u1.is_finite() && v0.is_finite() && v1.is_finite())
        || u1 <= u0
        || v1 <= v0
    {
        return face_to_triangles(face, deflection);
    }

    let def = deflection.max(1e-9);
    // OCCT route: the periodic-u direction (cylinder/cone/sphere/torus side)
    // subdivides by the `GCPnts_TangentialDeflection::ArcAngularStep` chord
    // step `2·acos(1 − deflection/radius)`, bounded by a cap so a tiny
    // deflection doesn't explode; the v direction subdivides by the actual
    // chord deviation (a cylinder's axial v-direction is straight, deviation
    // ~0, so it gets few segments — matching OCCT's ~2-segment side walls).
    let u_periodic = surface.is_u_periodic();
    let u_step: Option<f64> = if u_periodic {
        // Use the maximum u-ring radius sampled across the v range — for a
        // cone it grows 0→base, for a sphere it peaks at the equator and
        // collapses to 0 at the poles. Sampling avoids both the mid-range
        // underestimate (cone) and the pole-collapse (sphere, where the
        // v-range endpoints are the poles with radius 0).
        let mut r = 0.0f64;
        for k in 0..=8 {
            let vk = v0 + (v1 - v0) * k as f64 / 8.0;
            let rk = surface.d0(0.0, vk).distance(&surface.d0(std::f64::consts::PI, vk)) * 0.5;
            r = r.max(rk);
        }
        let r = r.max(1e-9);
        let ang_step = crate::meshing::range_splitter::arc_angular_step(r, def, 1.0, 0.0);
        Some(ang_step.max((u1 - u0) / 96.0))
    } else {
        None
    };
    // Cylinder/cone side walls have a straight v-direction generatrix; their v
    // band is a single segment (OCCT Delaunay does not refine a straight line).
    let straight_v = matches!(
        crate::pcurve_full::classify_surface_kind(surface.as_ref()),
        crate::brep_surface::SurfaceKind::Cylinder | crate::brep_surface::SurfaceKind::Cone
    );
    let mut verts: Vec<GpPnt> = Vec::new();
    let mut tris: Vec<Triangle> = Vec::new();
    let mut stack = vec![(u0, u1, v0, v1, 0usize)];
    // Shared-node table: a cell boundary point at (u, v) maps to one index, so
    // neighbouring cells share their common-edge vertices (OCCT Delaunay node
    // sharing). Without it every cell emitted its own 4 vertices — ~70% of the
    // mesh was duplicate points and the mesh was not watertight.
    let scale = 1.0 / def.max(1e-9);
    let mut nodes: HashMap<(i64, i64), usize> = HashMap::new();

    // Shared-node lookup: return the vertex index for (u, v), inserting it the
    // first time. `verts`/`nodes` are passed explicitly so the closure doesn't
    // hold conflicting borrows with the triangle assembly below.
    fn get_node(
        verts: &mut Vec<GpPnt>,
        nodes: &mut HashMap<(i64, i64), usize>,
        surface: &dyn occt_geom::Surface,
        scale: f64,
        u: f64,
        v: f64,
    ) -> usize {
        let key = ((u * scale).round() as i64, (v * scale).round() as i64);
        if let Some(&i) = nodes.get(&key) {
            i
        } else {
            let i = verts.len();
            verts.push(surface.d0(u, v));
            nodes.insert(key, i);
            i
        }
    }

    while let Some((cu0, cu1, cv0, cv1, depth)) = stack.pop() {
        let (mu, mv) = (0.5 * (cu0 + cu1), 0.5 * (cv0 + cv1));
        let p00 = surface.d0(cu0, cv0);
        let p10 = surface.d0(cu1, cv0);
        let p01 = surface.d0(cu0, cv1);
        let p11 = surface.d0(cu1, cv1);
        // u direction: fixed angular step for periodic surfaces, else chord
        // deviation. v direction: chord deviation (straight directions stay
        // coarse). Subdivide while either is not converged.
        let u_done = match u_step {
            Some(step) => (cu1 - cu0) <= step,
            None => {
                let dev = point_segment_dist(&surface.d0(mu, cv0), &p00, &p10)
                    .max(point_segment_dist(&surface.d0(mu, cv1), &p01, &p11));
                dev <= def
            }
        };
        let v_done = if straight_v {
            // Cylinder/cone side walls: the v direction is a straight
            // generatrix, so it subdivides to a small fixed number of bands
            // (OCCT's Delaunay gives the straight direction a minimal mesh
            // quality density — a cylinder's ~2, a cone's ~4), not the
            // hundreds the chord-deviation chase produced.
            (cv1 - cv0) <= (v1 - v0) / 4.0
        } else {
            let dev = point_segment_dist(&surface.d0(cu0, mv), &p00, &p01)
                .max(point_segment_dist(&surface.d0(cu1, mv), &p10, &p11))
                // Cell-centre deviation from the diagonals — catches high-order
                // curvature (a sphere's latitude arcs) that edge-midpoints miss.
                .max(point_segment_dist(&surface.d0(mu, mv), &p00, &p11))
                .max(point_segment_dist(&surface.d0(mu, mv), &p10, &p01));
            dev <= def
        };

        if !(u_done && v_done) && depth < MAX_FACE_DEPTH {
            // Subdivide only the unconverged direction(s) — OCCT splits U and V
            // independently (a straight cylinder v-direction stays coarse while
            // the angular u-direction refines).
            *iterations += 1;
            match (u_done, v_done) {
                (true, false) => {
                    stack.push((cu0, cu1, mv, cv1, depth + 1));
                    stack.push((cu0, cu1, cv0, mv, depth + 1));
                }
                (false, true) => {
                    stack.push((mu, cu1, cv0, cv1, depth + 1));
                    stack.push((cu0, mu, cv0, cv1, depth + 1));
                }
                (false, false) => {
                    stack.push((mu, cu1, mv, cv1, depth + 1));
                    stack.push((cu0, mu, mv, cv1, depth + 1));
                    stack.push((mu, cu1, cv0, mv, depth + 1));
                    stack.push((cu0, mu, cv0, mv, depth + 1));
                }
                (true, true) => {}
            }
        } else {
            let nvec = surface_normal(surface.as_ref(), mu, mv);
            let n = nvec.xyz();
            let i0 = get_node(&mut verts, &mut nodes, surface.as_ref(), scale, cu0, cv0);
            let i1 = get_node(&mut verts, &mut nodes, surface.as_ref(), scale, cu1, cv0);
            let i2 = get_node(&mut verts, &mut nodes, surface.as_ref(), scale, cu0, cv1);
            let i3 = get_node(&mut verts, &mut nodes, surface.as_ref(), scale, cu1, cv1);
            let t1 = orient3(&verts, i0, i1, i2, &n);
            let t2 = orient3(&verts, i1, i3, i2, &n);
            if tri_area2(&verts, &t1) > 1e-24 {
                tris.push(t1);
            }
            if tri_area2(&verts, &t2) > 1e-24 {
                tris.push(t2);
            }
        }
    }
    (verts, tris)
}

/// Maximum deviation of a face's triangle mesh from its surface.
///
/// Each triangle's vertices are projected back to surface parameters; the
/// surface is then sampled on a barycentric `samples`-point grid per triangle
/// and the furthest distance to the triangle plane is returned.
pub fn mesh_deflection_error(mesh: &ShapeMesh, face: &Face, samples: usize) -> f64 {
    let Some(surface) = BRepTool::face_surface(face) else {
        return 0.0;
    };
    let (u0, u1) = surface.u_range();
    let period = if surface.is_u_periodic() && u1.is_finite() && u0.is_finite() {
        u1 - u0
    } else {
        0.0
    };
    let samples = samples.max(2);
    let mut max_dev = 0.0f64;
    for t in &mesh.triangles {
        let a = &mesh.vertices[t.n0];
        let b = &mesh.vertices[t.n1];
        let c = &mesh.vertices[t.n2];
        let n = b.coord.subtracted(&a.coord).crossed(&c.coord.subtracted(&a.coord));
        if n.square_modulus() < 1e-30 {
            continue;
        }
        let n = n.divided(n.modulus());
        let (ua, va) = surface_closest_params(surface.as_ref(), a, 24, 24);
        let (mut ub, vb) = surface_closest_params(surface.as_ref(), b, 24, 24);
        let (mut uc, vc) = surface_closest_params(surface.as_ref(), c, 24, 24);
        // Periodic-u surfaces report seam points at either end of the range
        // (e.g. u=0 or u=2π on a sphere). Unwrap so barycentric interpolation
        // stays on the same sheet and does not wrap to the far side.
        if period > 0.0 {
            let half = 0.5 * period;
            if (ub - ua).abs() > half {
                ub += if ub > ua { -period } else { period };
            }
            if (uc - ua).abs() > half {
                uc += if uc > ua { -period } else { period };
            }
        }
        for i in 0..=samples {
            for j in 0..=(samples - i) {
                let w = i as f64 / samples as f64;
                let v = j as f64 / samples as f64;
                let u = 1.0 - w - v;
                let uu = u * ua + w * ub + v * uc;
                let vv = u * va + w * vb + v * vc;
                let p = surface.d0(uu, vv);
                let d = p.coord.subtracted(&a.coord).dot(&n).abs();
                if d > max_dev {
                    max_dev = d;
                }
            }
        }
    }
    max_dev
}

/// Iterative midpoint subdivision of a curve into a deflection-bounded
/// polyline. Bounded to [`MAX_CURVE_DEPTH`] recursion levels. Curves with an
/// unbounded parameter range fall back to 64 uniform samples over the natural
/// `[0, 1]` window.
pub fn adaptive_curve_polyline(c: &dyn Curve, deflection: f64) -> Vec<GpPnt> {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !(a.is_finite() && b.is_finite() && b >= a) {
        return occt_geom::curve_approx::curve_to_polyline_uniform(c, 64);
    }
    if b - a < 1e-12 {
        return vec![c.d0(a)];
    }
    let def = deflection.max(1e-12);
    let mut pts: Vec<GpPnt> = Vec::new();
    let mut last_param: Option<f64> = None;
    let mut stack = vec![(a, b, 0usize)];
    while let Some((lo, hi, depth)) = stack.pop() {
        let mid = 0.5 * (lo + hi);
        let (pa, pm, pb) = (c.d0(lo), c.d0(mid), c.d0(hi));
        if point_segment_dist(&pm, &pa, &pb) > def && depth < MAX_CURVE_DEPTH {
            stack.push((mid, hi, depth + 1));
            stack.push((lo, mid, depth + 1));
        } else {
            if last_param != Some(lo) {
                pts.push(pa);
                last_param = Some(lo);
            }
            if last_param != Some(hi) {
                pts.push(pb);
                last_param = Some(hi);
            }
        }
    }
    if pts.is_empty() {
        pts.push(c.d0(a));
        pts.push(c.d0(b));
    }
    pts
}

/// Mesh quality summary: `(triangle_count, degenerate_count, min_edge, max_edge)`.
pub fn mesh_quality_report(mesh: &ShapeMesh) -> (usize, usize, f64, f64) {
    let total = mesh.triangles.len();
    let mut deg = 0usize;
    let mut min_len = f64::INFINITY;
    let mut max_len = 0.0f64;
    for t in &mesh.triangles {
        let a = &mesh.vertices[t.n0];
        let b = &mesh.vertices[t.n1];
        let c = &mesh.vertices[t.n2];
        let l1 = a.distance(b);
        let l2 = b.distance(c);
        let l3 = c.distance(a);
        min_len = min_len.min(l1).min(l2).min(l3);
        max_len = max_len.max(l1).max(l2).max(l3);
        if b.coord.subtracted(&a.coord).crossed(&c.coord.subtracted(&a.coord)).square_modulus() < 1e-24 {
            deg += 1;
        }
    }
    if total == 0 {
        min_len = 0.0;
    }
    (total, deg, min_len, max_len)
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Distance from `p` to the segment `a..b`.
fn point_segment_dist(p: &GpPnt, a: &GpPnt, b: &GpPnt) -> f64 {
    let ab = b.coord.subtracted(&a.coord);
    let len2 = ab.square_modulus();
    if len2 <= f64::EPSILON {
        return p.distance(a);
    }
    let ap = p.coord.subtracted(&a.coord);
    let t = (ap.dot(&ab) / len2).clamp(0.0, 1.0);
    let proj = a.coord.added(&ab.multiplied(t));
    p.coord.subtracted(&proj).modulus()
}

/// Squared double area of a triangle.
fn tri_area2(pts: &[GpPnt], t: &Triangle) -> f64 {
    let ab = pts[t.n1].coord.subtracted(&pts[t.n0].coord);
    let ac = pts[t.n2].coord.subtracted(&pts[t.n0].coord);
    ab.crossed(&ac).square_modulus()
}

/// Wind a triangle so its normal agrees with `n`.
fn orient3(pts: &[GpPnt], i0: usize, i1: usize, i2: usize, n: &GpXyz) -> Triangle {
    let pa = pts[i0].coord;
    let pb = pts[i1].coord;
    let pc = pts[i2].coord;
    let nn = pb.subtracted(&pa).crossed(&pc.subtracted(&pa));
    if n.dot(&nn) < 0.0 {
        Triangle::new(i0, i2, i1)
    } else {
        Triangle::new(i0, i1, i2)
    }
}

/// Whether a (planar, non-self-intersecting) polygon is convex when seen along
/// `n` (all consecutive-triple signed areas share a sign).
fn is_convex(pts: &[GpPnt], n: &GpXyz) -> bool {
    let mut sign: Option<f64> = None;
    let m = pts.len();
    for i in 0..m {
        let a = pts[i].coord;
        let b = pts[(i + 1) % m].coord;
        let c = pts[(i + 2) % m].coord;
        let d = b.subtracted(&a).crossed(&c.subtracted(&a)).dot(n);
        if d.abs() > 1e-12 {
            match sign {
                None => sign = Some(d.signum()),
                Some(s) if s != d.signum() => return false,
                _ => {}
            }
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    use occt_core::gp::{GpAx3, GpSphere as GpSphereCore};
    use occt_geom::{GeomCircle, GeomSphere, GeomTrimmedCurve};
    use crate::builder::TopoBuilder;
    use crate::mesh::mesh_surface_area;
    use crate::primitives::{BRepPrimBox, BRepPrimSphere};

    fn unit_sphere_face() -> Face {
        let b = TopoBuilder::new();
        b.make_face(
            Arc::new(GeomSphere::new(
                GpSphereCore::new(GpAx3::standard(), 2.0).unwrap(),
            )),
            &[],
        )
    }

    #[test]
    fn sphere_mesh_area_within_five_percent() {
        let face = unit_sphere_face();
        let (pts, tris) = adaptive_face_mesh(&face, 0.05, &mut 0);
        let area: f64 = tris
            .iter()
            .map(|t| 0.5 * tri_area2(&pts, t).sqrt())
            .sum();
        let expect = 4.0 * std::f64::consts::PI * 4.0;
        assert!(
            (area - expect).abs() / expect < 0.05,
            "sphere area {area} vs {expect}"
        );
        let mesh = ShapeMesh {
            vertices: pts,
            triangles: tris,
            source_shape: crate::abs::ShapeType::Face,
        };
        let err = mesh_deflection_error(&mesh, &face, 6);
        assert!(err < 0.05, "deflection error {err}");
    }

    #[test]
    fn box_meshes_exactly() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let im = incremental_mesh(&b.solid.0, 0.1).expect("mesh box");
        assert_eq!(im.mesh.triangles.len(), 12, "expected 12 triangles");
        assert!(
            (mesh_surface_area(&im.mesh) - 6.0).abs() < 1e-9,
            "area {}",
            mesh_surface_area(&im.mesh)
        );
        let (total, deg, _, _) = mesh_quality_report(&im.mesh);
        assert_eq!(total, 12);
        assert_eq!(deg, 0);
    }

    #[test]
    fn adaptive_curve_semicircle() {
        let circ = GeomCircle::new(occt_core::gp::GpCirc::new(
            occt_core::gp::GpAx2::standard(),
            2.0,
        ));
        let pi = std::f64::consts::PI;
        let trimmed = GeomTrimmedCurve::new(Arc::new(circ), 0.0, pi);
        let deflection = 0.05;
        let poly = adaptive_curve_polyline(&trimmed, deflection);
        assert!(poly.len() >= 2, "poly len {}", poly.len());
        assert!(
            poly.first().unwrap().distance(&GpPnt::new(2.0, 0.0, 0.0)) < 1e-6,
            "first endpoint {:?}",
            poly.first()
        );
        assert!(
            poly.last().unwrap().distance(&GpPnt::new(-2.0, 0.0, 0.0)) < 1e-6,
            "last endpoint {:?}",
            poly.last()
        );
        // Chord deviation of the polyline from the true curve. Each polyline
        // vertex sits on the circle at angle u = atan2(y, x); the trimmed curve
        // is reparameterized to [0,1], so its parameter is angle/π. Sample the
        // arc between consecutive vertices and measure distance to the chord.
        let mut max_dev = 0.0f64;
        for w in poly.windows(2) {
            let ua = w[0].y().atan2(w[0].x());
            let ub = w[1].y().atan2(w[1].x());
            if ub < ua {
                continue;
            }
            for k in 1..16 {
                let u = (ua + (ub - ua) * k as f64 / 16.0) / pi;
                let dev = point_segment_dist(&trimmed.d0(u), &w[0], &w[1]);
                if dev > max_dev {
                    max_dev = dev;
                }
            }
        }
        assert!(max_dev < deflection, "chord deviation {max_dev}");
    }

    #[test]
    fn deflection_error_decreases_with_finer_deflection() {
        let face = unit_sphere_face();
        let (pts1, tris1) = adaptive_face_mesh(&face, 0.2, &mut 0);
        let (pts2, tris2) = adaptive_face_mesh(&face, 0.02, &mut 0);
        let m1 = ShapeMesh {
            vertices: pts1,
            triangles: tris1,
            source_shape: crate::abs::ShapeType::Face,
        };
        let m2 = ShapeMesh {
            vertices: pts2,
            triangles: tris2,
            source_shape: crate::abs::ShapeType::Face,
        };
        let e1 = mesh_deflection_error(&m1, &face, 4);
        let e2 = mesh_deflection_error(&m2, &face, 4);
        assert!(
            e2 < e1,
            "expected finer mesh to have smaller deflection error: {e2} vs {e1}"
        );
    }

    #[test]
    fn full_sphere_incremental() {
        let s = BRepPrimSphere::make_sphere(2.0);
        let im = incremental_mesh(&s.solid.0, 0.05).expect("mesh sphere");
        assert!(im.mesh.triangles.len() > 0);
        let (_, deg, _, _) = mesh_quality_report(&im.mesh);
        // Poles produce a few degenerate triangles; keep the count small.
        assert!(deg <= 2, "degenerate count {deg}");
    }
}
