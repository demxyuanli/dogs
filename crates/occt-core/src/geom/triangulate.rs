//! Polygon triangulation (ear clipping) + mesh utilities.
use crate::gp::GpPnt;

/// Orientation test in XY plane: >0 = CCW, <0 = CW.
fn orient2d(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> f64 {
    (b.x() - a.x()) * (c.y() - a.y()) - (b.y() - a.y()) * (c.x() - a.x())
}

/// Is point inside triangle (CCW triangle, inclusive edges)?
fn point_in_triangle(p: &GpPnt, a: &GpPnt, b: &GpPnt, c: &GpPnt) -> bool {
    let d1 = orient2d(a, b, p);
    let d2 = orient2d(b, c, p);
    let d3 = orient2d(c, a, p);
    let has_neg = d1 < 0.0 || d2 < 0.0 || d3 < 0.0;
    let has_pos = d1 > 0.0 || d2 > 0.0 || d3 > 0.0;
    !(has_neg && has_pos)
}

/// Triangulate a simple polygon (in XY plane) via ear clipping.
/// polygon: vertices in order (CW or CCW). Returns triangles as index triples.
/// Assumes polygon is simple (non-self-intersecting).
pub fn triangulate_polygon(polygon: &[GpPnt]) -> Vec<(usize, usize, usize)> {
    let n = polygon.len();
    if n < 3 { return vec![]; }
    if n == 3 { return vec![(0, 1, 2)]; }

    // Determine orientation and reverse if CW
    let mut area = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        area += polygon[i].x() * polygon[j].y() - polygon[j].x() * polygon[i].y();
    }
    let ccw = area > 0.0;

    // Work on index list
    let mut indices: Vec<usize> = (0..n).collect();
    if !ccw { indices.reverse(); }

    let mut triangles = Vec::new();
    let mut guard = 0usize;
    while indices.len() > 3 && guard < n * n {
        guard += 1;
        let m = indices.len();
        let mut clipped = false;
        for k in 0..m {
            let i0 = indices[(k + m - 1) % m];
            let i1 = indices[k];
            let i2 = indices[(k + 1) % m];
            let a = polygon[i0]; let b = polygon[i1]; let c = polygon[i2];

            // Ear if corner is convex and contains no other polygon vertex
            if orient2d(&a, &b, &c) <= 0.0 { continue; }
            let mut is_ear = true;
            for &vi in &indices {
                if vi == i0 || vi == i1 || vi == i2 { continue; }
                if point_in_triangle(&polygon[vi], &a, &b, &c) { is_ear = false; break; }
            }
            if is_ear {
                triangles.push((i0, i1, i2));
                indices.remove(k);
                clipped = true;
                break;
            }
        }
        if !clipped { break; } // degenerate polygon
    }
    if indices.len() == 3 {
        triangles.push((indices[0], indices[1], indices[2]));
    }
    triangles
}

/// Triangulate and return the total area (for validation).
pub fn triangulation_area(polygon: &[GpPnt]) -> f64 {
    let tris = triangulate_polygon(polygon);
    let mut area = 0.0;
    for (a, b, c) in tris {
        area += 0.5 * orient2d(&polygon[a], &polygon[b], &polygon[c]).abs();
    }
    area
}

/// Simple quad mesh → triangle mesh (split each quad into 2 triangles).
/// quads: 4-vertex index groups.
pub fn quads_to_triangles(quads: &[usize]) -> Vec<(usize, usize, usize)> {
    let mut tris = Vec::new();
    for q in quads.chunks(4) {
        if q.len() == 4 {
            tris.push((q[0], q[1], q[2]));
            tris.push((q[0], q[2], q[3]));
        }
    }
    tris
}

/// Count edges in a triangle mesh (with dedup) — for Euler characteristic.
pub fn count_unique_edges(triangles: &[(usize, usize, usize)]) -> usize {
    let mut edges = std::collections::HashSet::new();
    for &(a, b, c) in triangles {
        for (u, v) in [(a, b), (b, c), (c, a)] {
            edges.insert(if u < v { (u, v) } else { (v, u) });
        }
    }
    edges.len()
}

/// Compute signed volume of a closed triangle mesh (origin-based tetrahedra).
pub fn mesh_signed_volume(verts: &[GpPnt], tris: &[(usize, usize, usize)]) -> f64 {
    let mut vol = 0.0;
    for &(i, j, k) in tris {
        let a = verts[i].coord; let b = verts[j].coord; let c = verts[k].coord;
        vol += a.dot_cross(&b, &c) / 6.0;
    }
    vol
}

/// Average edge length of a mesh (scaling metric).
pub fn average_edge_length(verts: &[GpPnt], tris: &[(usize, usize, usize)]) -> f64 {
    let mut total = 0.0; let mut count = 0usize;
    for &(a, b, c) in tris {
        total += verts[a].coord.subtracted(&verts[b].coord).modulus();
        total += verts[b].coord.subtracted(&verts[c].coord).modulus();
        total += verts[c].coord.subtracted(&verts[a].coord).modulus();
        count += 3;
    }
    if count == 0 { 0.0 } else { total / count as f64 }
}

/// Axis-aligned bounding box of a vertex set.
pub fn vertices_bbox(verts: &[GpPnt]) -> crate::bnd::BndBox {
    let mut b = crate::bnd::BndBox::new();
    for v in verts { b.add_point(v); }
    b
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn triangulate_square() {
        let sq = vec![
            GpPnt::new(0.,0.,0.), GpPnt::new(1.,0.,0.), GpPnt::new(1.,1.,0.), GpPnt::new(0.,1.,0.),
        ];
        let tris = triangulate_polygon(&sq);
        assert_eq!(tris.len(), 2);
        assert!((triangulation_area(&sq) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn triangulate_concave() {
        // L-shape: concave polygon
        let l = vec![
            GpPnt::new(0.,0.,0.), GpPnt::new(2.,0.,0.), GpPnt::new(2.,1.,0.),
            GpPnt::new(1.,1.,0.), GpPnt::new(1.,2.,0.), GpPnt::new(0.,2.,0.),
        ];
        let tris = triangulate_polygon(&l);
        assert_eq!(tris.len(), 4);
        // L-shape area = 2*1 + 1*1 = 3
        assert!((triangulation_area(&l) - 3.0).abs() < 1e-12);
    }

    #[test]
    fn unique_edges() {
        let tris = vec![(0, 1, 2), (0, 2, 3)];
        assert_eq!(count_unique_edges(&tris), 5);
    }
}
