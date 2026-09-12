use super::prelude::*;
use super::*;


// ---------------------------------------------------------------------------
// 2D polygon triangulation
// ---------------------------------------------------------------------------

/// Triangulate a simple 2D polygon by ear clipping.
///
/// Wraps the shared ear-clipping implementation in
/// [`crate::geom::polygon_ops::triangulate_polygon2d`], which accepts a simple
/// polygon in either winding order and repeatedly clips convex "ears"
/// (corners that contain no other polygon vertex) until a single triangle
/// remains. The resulting triangle index triples reference positions in `pts`.
///
/// * `pts` — polygon vertices in boundary order (CCW or CW).
///
/// Returns `Some(triangles)` on success, or `None` for fewer than three
/// vertices, degenerate (zero-area) or self-intersecting input.
pub fn polygon_triangulate_ear(pts: &[GpPnt2d]) -> Option<Vec<(usize, usize, usize)>> {
    crate::geom::polygon_ops::triangulate_polygon2d(pts)
}
// ---------------------------------------------------------------------------
// Vertex-connectivity
// ---------------------------------------------------------------------------

/// Compute the per-vertex incident triangle lists of a mesh.
///
/// For each vertex index this returns the list of triangles that reference it,
/// i.e. the 1-ring of faces around the vertex. This is the data structure
/// `Poly_Connect::TriangleToNodes` / `Poly_Triangulation` adjacency tables
/// expose, and it is the starting point for vertex-normal averaging, boundary
/// extraction and local mesh editing.
///
/// * `verts` — the mesh vertex array; its length bounds the result.
/// * `triangles` — the mesh as index triples.
///
/// Returns one `Vec<usize>` per vertex; vertices that are not referenced by
/// any triangle get an empty list. Triangle indices are pushed in ascending
/// triangle order, so each list is naturally sorted.
pub fn mesh_to_polygon_indices(
    verts: &[GpPnt],
    triangles: &[(usize, usize, usize)],
) -> Vec<Vec<usize>> {
    let mut out: Vec<Vec<usize>> = vec![Vec::new(); verts.len()];
    for (ti, &(a, b, c)) in triangles.iter().enumerate() {
        if a < verts.len() {
            out[a].push(ti);
        }
        if b < verts.len() {
            out[b].push(ti);
        }
        if c < verts.len() {
            out[c].push(ti);
        }
    }
    out
}
