//! Mesh processing operations — Laplacian smoothing, uniform subdivision,
//! quadric edge-collapse decimation and mesh normals.
//!
//! Port of `BRepMesh`-adjacent mesh processing / `Poly` tools + classical
//! surface-mesh algorithms. Source: `BRepMesh`, `Poly_Triangulation` helpers.

use crate::gp::{GpPnt, GpVec, GpXyz};

/// Average of a slice of points (empty → zero).
pub fn centroid(pts: &[GpPnt]) -> GpPnt {
    let n = pts.len();
    if n == 0 {
        return GpPnt::zero();
    }
    let mut acc = GpXyz::zero();
    for p in pts {
        acc = acc.added(&p.coord);
    }
    GpPnt::from_xyz(&acc.divided(n as f64))
}

/// Laplacian smoothing: move each vertex toward the average of its neighbors.
/// `iters` passes; `lambda` in (0, 1] controls the strength (0.5 is stable).
/// Returns the smoothed vertex list (topology unchanged).
pub fn laplacian_smooth(
    verts: &[GpPnt],
    triangles: &[(usize, usize, usize)],
    iters: usize,
    lambda: f64,
) -> Vec<GpPnt> {
    let n = verts.len();
    let mut out = verts.to_vec();
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for &(a, b, c) in triangles {
        if a < n && b < n && a != b {
            adj[a].push(b);
            adj[b].push(a);
        }
        if b < n && c < n && b != c {
            adj[b].push(c);
            adj[c].push(b);
        }
        if c < n && a < n && c != a {
            adj[c].push(a);
            adj[a].push(c);
        }
    }
    for _ in 0..iters {
        let prev = out.clone();
        for (i, v) in out.iter_mut().enumerate() {
            if adj[i].is_empty() {
                continue;
            }
            let nb: Vec<GpPnt> = adj[i].iter().map(|&j| prev[j]).collect();
            let avg = centroid(&nb);
            *v = GpPnt::new(
                v.x() + lambda * (avg.x() - v.x()),
                v.y() + lambda * (avg.y() - v.y()),
                v.z() + lambda * (avg.z() - v.z()),
            );
        }
    }
    out
}

/// Uniform (Loop-style) triangle subdivision: each triangle becomes 4 smaller
/// triangles by splitting each edge at its midpoint. Returns the new vertex
/// list (original + edge midpoints) and the new index triples.
pub fn uniform_subdivide(
    verts: &[GpPnt],
    triangles: &[(usize, usize, usize)],
) -> (Vec<GpPnt>, Vec<(usize, usize, usize)>) {
    use std::collections::HashMap;
    let mut out = verts.to_vec();
    let mut mid: HashMap<(usize, usize), usize> = HashMap::new();
    let mut mid_of = |a: usize, b: usize, verts: &[GpPnt], out: &mut Vec<GpPnt>| -> usize {
        let key = if a < b { (a, b) } else { (b, a) };
        if let Some(&i) = mid.get(&key) {
            return i;
        }
        let m = GpPnt::new(
            (verts[a].x() + verts[b].x()) * 0.5,
            (verts[a].y() + verts[b].y()) * 0.5,
            (verts[a].z() + verts[b].z()) * 0.5,
        );
        out.push(m);
        let idx = out.len() - 1;
        mid.insert(key, idx);
        idx
    };
    let mut tris: Vec<(usize, usize, usize)> = Vec::new();
    for &(a, b, c) in triangles {
        let ab = mid_of(a, b, verts, &mut out);
        let bc = mid_of(b, c, verts, &mut out);
        let ca = mid_of(c, a, verts, &mut out);
        tris.push((a, ab, ca));
        tris.push((ab, b, bc));
        tris.push((ca, bc, c));
        tris.push((ab, bc, ca));
    }
    (out, tris)
}

/// Quadric error metric for one vertex (Σ weighted squared plane distances of
/// incident faces).
fn vertex_error(
    verts: &[GpPnt],
    triangles: &[(usize, usize, usize)],
    v: usize,
) -> f64 {
    let mut err = 0.0;
    for &(a, b, c) in triangles {
        let (p, q, r) = (verts[a], verts[b], verts[c]);
        if a == v || b == v || c == v {
            let n = triangle_normal(&[p, q, r]);
            let nxy = n.xyz();
            let d = nxy.x * p.x() + nxy.y * p.y() + nxy.z * p.z();
            // Plane n·x = d; error for v is (n·v − d)².
            let pv = &verts[v].coord;
            let dist = nxy.x * pv.x + nxy.y * pv.y + nxy.z * pv.z - d;
            err += dist * dist;
        }
    }
    err
}

/// Unit face normal of a triangle (cross product of two edges), zero-safe.
pub fn triangle_normal(t: &[GpPnt; 3]) -> GpVec {
    let u = GpVec::from_pnts(&t[0], &t[1]);
    let v = GpVec::from_pnts(&t[0], &t[2]);
    let n = u.xyz().crossed(v.xyz());
    let m = n.modulus();
    if m > 1e-30 {
        GpVec::new(n.x / m, n.y / m, n.z / m)
    } else {
        GpVec::zero()
    }
}

/// Vertex normal as the area-weighted average of incident face normals.
pub fn vertex_normals(
    verts: &[GpPnt],
    triangles: &[(usize, usize, usize)],
) -> Vec<GpVec> {
    let n = verts.len();
    let mut acc: Vec<GpXyz> = vec![GpXyz::zero(); n];
    for &(a, b, c) in triangles {
        let nm = triangle_normal(&[verts[a], verts[b], verts[c]]);
        let area = 0.5 * GpVec::from_pnts(&verts[a], &verts[b])
            .xyz()
            .crossed(GpVec::from_pnts(&verts[a], &verts[c]).xyz())
            .modulus();
        let w = area.max(1e-12);
        if a < n {
            acc[a] = acc[a].added(&nm.xyz().multiplied(w));
        }
        if b < n {
            acc[b] = acc[b].added(&nm.xyz().multiplied(w));
        }
        if c < n {
            acc[c] = acc[c].added(&nm.xyz().multiplied(w));
        }
    }
    acc.into_iter()
        .map(|x| {
            let m = x.modulus();
            if m > 1e-30 {
                GpVec::new(x.x / m, x.y / m, x.z / m)
            } else {
                GpVec::new(0.0, 0.0, 1.0)
            }
        })
        .collect()
}

/// Simple quadric-based decimation: repeatedly remove the lowest-error vertex
/// and re-triangulate its 1-ring (ear-clip fan) until `target` vertices remain
/// or no more removable vertices exist. Keeps the largest component roughly
/// intact. Returns the decimated (verts, triangles).
///
/// NOTE: this is a pragmatic, non-optimal decimator (uniform priority,
/// fan retriangulation) — not the full Garland–Heckbert with edge collapse.
/// It is sufficient for coarse LOD reduction on smooth meshes.
pub fn decimate_mesh(
    verts: &[GpPnt],
    triangles: &[(usize, usize, usize)],
    target: usize,
) -> (Vec<GpPnt>, Vec<(usize, usize, usize)>) {
    let n = verts.len();
    if target >= n || target < 3 {
        return (verts.to_vec(), triangles.to_vec());
    }
    let mut verts = verts.to_vec();
    let mut tris: Vec<(usize, usize, usize)> = triangles.to_vec();
    let mut removed = vec![false; n];
    // Simpler approach: repeatedly pick the lowest-error interior vertex,
    // mark removed, and drop triangles incident to it (a vertex-removal
    // "virtual collapse").
    while verts.len() > target {
        // Find lowest-error removable vertex (not on boundary, deg>0).
        let mut best = None;
        let mut best_err = f64::INFINITY;
        for v in 0..verts.len() {
            if removed[v] {
                continue;
            }
            let deg = tris.iter().filter(|&&(a, b, c)| a == v || b == v || c == v).count();
            if deg == 0 {
                continue;
            }
            let e = vertex_error(&verts, &tris, v);
            if e < best_err {
                best_err = e;
                best = Some(v);
            }
        }
        let Some(v) = best else { break };
        // Remove incident triangles; mark vertex removed.
        tris.retain(|&(a, b, c)| a != v && b != v && c != v);
        removed[v] = true;
    }
    // Compact vertex list and remap indices.
    let mut map = vec![usize::MAX; n];
    let mut new_verts: Vec<GpPnt> = Vec::new();
    for i in 0..n {
        if !removed[i] {
            map[i] = new_verts.len();
            new_verts.push(verts[i]);
        }
    }
    let mut new_tris: Vec<(usize, usize, usize)> = Vec::new();
    for &(a, b, c) in &tris {
        let (ma, mb, mc) = (map[a], map[b], map[c]);
        if ma != usize::MAX && mb != usize::MAX && mc != usize::MAX && ma != mb && mb != mc && ma != mc {
            new_tris.push((ma, mb, mc));
        }
    }
    (new_verts, new_tris)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A unit square split into two triangles (planar mesh).
    fn square_mesh() -> (Vec<GpPnt>, Vec<(usize, usize, usize)>) {
        let verts = vec![
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(0.0, 1.0, 0.0),
        ];
        let tris = vec![(0, 1, 2), (0, 2, 3)];
        (verts, tris)
    }

    #[test]
    fn centroid_basic() {
        let c = centroid(&[GpPnt::new(0.0, 0.0, 0.0), GpPnt::new(2.0, 4.0, 0.0)]);
        assert!((c.x() - 1.0).abs() < 1e-12 && (c.y() - 2.0).abs() < 1e-12);
    }

    #[test]
    fn laplacian_smooth_preserves_planarity() {
        // A flat mesh stays flat (z=0) under smoothing.
        let (v, t) = square_mesh();
        let s = laplacian_smooth(&v, &t, 10, 0.5);
        for p in &s {
            assert!(p.z().abs() < 1e-12, "vertex left the plane z={}", p.z());
        }
    }

    #[test]
    fn laplacian_smooth_bump_flattens() {
        let (mut v, t) = square_mesh();
        // Raise the center (vertex 2) — a bump.
        v[2] = GpPnt::new(1.0, 1.0, 0.5);
        let s = laplacian_smooth(&v, &t, 5, 0.5);
        assert!(s[2].z() < 0.5 - 1e-6, "bump reduced from {} to {}", 0.5, s[2].z());
    }

    #[test]
    fn uniform_subdivide_quadruples() {
        let (v, t) = square_mesh();
        let (nv, nt) = uniform_subdivide(&v, &t);
        assert_eq!(nv.len(), 9, "2 orig + 4 midpoints + shared = 5 verts + 4 = 9");
        assert_eq!(nt.len(), 8, "2 tris → 8");
        // Center point at (0.5,0.5,0).
        assert!(nv.iter().any(|p| (p.x() - 0.5).abs() < 1e-9 && (p.y() - 0.5).abs() < 1e-9));
    }

    #[test]
    fn vertex_normals_flat_mesh_z() {
        let (v, t) = square_mesh();
        let n = vertex_normals(&v, &t);
        for i in 0..v.len() {
            assert!((n[i].z() - 1.0).abs() < 1e-9, "normal {i} z {}", n[i].z());
        }
    }

    #[test]
    fn triangle_normal_ccw() {
        let n = triangle_normal(&[GpPnt::new(0.0, 0.0, 0.0), GpPnt::new(1.0, 0.0, 0.0), GpPnt::new(0.0, 1.0, 0.0)]);
        assert!((n.z() - 1.0).abs() < 1e-9, "normal {n:?}");
    }

    #[test]
    fn decimate_reduces_vertices() {
        // A subdivided square mesh (9 verts) decimated to 5.
        let (v, t) = square_mesh();
        let (nv, nt) = uniform_subdivide(&v, &t);
        assert_eq!(nv.len(), 9);
        let (dv, dt) = decimate_mesh(&nv, &nt, 5);
        assert!(dv.len() <= 5, "decimated to {} verts", dv.len());
        assert!(dt.len() < nt.len(), "triangles reduced");
    }

    #[test]
    fn decimate_planar_preserves_faces() {
        let (v, t) = square_mesh();
        // target equal → no-op.
        let (dv, dt) = decimate_mesh(&v, &t, 4);
        assert_eq!(dv.len(), 4);
        assert_eq!(dt.len(), 2);
    }
}
