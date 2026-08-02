//! Port of OCCT Delabella — Wave 4 BRepMesh.
//!
//! `BRepMesh_DelabellaBaseMeshAlgo` triangulates a face's UV point set with the
//! third-party "Delabella" Delaunay library. This port replaces the C++ library
//! with an independent **incremental-insertion** Delaunay triangulator written
//! directly against the [`DelaunDataStructure`] contract, so the module is
//! self-contained and does not depend on the Wave 2 [`Delaun`] mesher.
//!
//! The triangulation uses the split-and-flip (Lawson) incremental-insertion
//! variant of Bowyer–Watson: each point is inserted into the triangle that
//! contains it (splitting it, or the two triangles across the edge when the
//! point falls on an edge), then the affected edges are flipped until every
//! circumcircle is empty. Unlike the classic cavity formulation, this stays
//! valid for cocircular point sets (e.g. regular grids) and for points that
//! land exactly on an existing edge. The empty-circle predicate mirrors the
//! OCCT Delaunay core (and `occt-core::geom::delaunay`).
//!
//! Sources (`src/ModelingAlgorithms/TKMesh/BRepMesh/`):
//! - `BRepMesh_DelabellaBaseMeshAlgo.{hxx,cxx}`
//! - `BRepMesh_DelabellaMeshAlgoFactory.{hxx,cxx}`
//! - `delabella.cpp` (the triangulation engine, replaced by
//!   [`lawson_triangulate`]).

use std::collections::HashSet;

use occt_core::gp::GpPnt2d;

use super::data_model::{MeshModel, MeshStatus};
use super::delaun_data::*;
use super::delaun_types::*;
use super::mesh_algo::{build_data_structure, BaseMeshAlgo};
use super::parameters::MeshParameters;

/// Delaunay triangulation of a face's UV point set via the Delabella algorithm.
///
/// Source: `BRepMesh_DelabellaBaseMeshAlgo` (a `CustomBaseMeshAlgo` in OCCT).
/// The algorithm collects the structure's live nodes, triangulates them with
/// [`lawson_triangulate`] and writes the resulting `Free` triangles back into
/// the structure. The finished structure is queried through
/// [`DelabellaBaseMeshAlgo::result`].
pub struct DelabellaBaseMeshAlgo {
    structure: Option<DelaunDataStructure>,
}

impl DelabellaBaseMeshAlgo {
    /// Creates an idle algorithm.
    pub fn new() -> Self {
        Self { structure: None }
    }

    /// The triangulated data structure produced by the last successful run.
    pub fn result(&self) -> Option<&DelaunDataStructure> {
        self.structure.as_ref()
    }

    /// Triangulates the UV point set of the given data structure and writes the
    /// resulting triangles back into it.
    ///
    /// Mirrors `BRepMesh_DelabellaBaseMeshAlgo::buildBaseTriangulation`. Fails
    /// when fewer than three live nodes exist or the point set is degenerate
    /// (coincident / collinear).
    pub fn triangulate(
        &mut self,
        mut structure: DelaunDataStructure,
        _parameters: &MeshParameters,
    ) -> Result<(), String> {
        // Collect live nodes, keeping their 1-based structure indices.
        let mut indices: Vec<i32> = Vec::new();
        let mut points: Vec<GpPnt2d> = Vec::new();
        for index in 1..=structure.nb_nodes() as i32 {
            let node = structure.get_node(index);
            if node.state != VertexState::Deleted {
                indices.push(index);
                points.push(node.location);
            }
        }
        let nodes_nb = points.len();
        if nodes_nb < 3 {
            return Err(format!(
                "DelabellaBaseMeshAlgo::triangulate: fewer than 3 live nodes ({nodes_nb})"
            ));
        }

        let triangles = lawson_triangulate(&points);
        if triangles.is_empty() {
            return Err("DelabellaBaseMeshAlgo::triangulate: degenerate point set".to_string());
        }

        // Replace the previous domain with the freshly computed triangulation.
        // Constraint (Frontier/Fixed) links survive; the new triangles reuse
        // them when an edge coincides with the face boundary.
        structure.clear_domain();
        for &(a, b, c) in &triangles {
            let va = indices[a];
            let vb = indices[b];
            let vc = indices[c];
            let lab = structure.add_link(va, vb, VertexState::Free);
            let lbc = structure.add_link(vb, vc, VertexState::Free);
            let lca = structure.add_link(vc, va, VertexState::Free);
            structure.add_element(DelaunTriangle::new([lab, lbc, lca], [va, vb, vc]));
        }

        self.structure = Some(structure);
        Ok(())
    }
}

impl Default for DelabellaBaseMeshAlgo {
    fn default() -> Self {
        Self::new()
    }
}

impl BaseMeshAlgo for DelabellaBaseMeshAlgo {
    fn perform(&mut self, model: &mut MeshModel, parameters: &MeshParameters) -> Result<(), String> {
        // OCCT invokes the base algorithm once per face. The model-level port
        // triangulates the first face; the produced structure is queried
        // through `result()`.
        const FIRST_FACE: usize = 0;
        match build_data_structure(model, FIRST_FACE)
            .and_then(|structure| self.triangulate(structure, parameters))
        {
            Ok(()) => {
                if let Ok(face) = model.face_mut(FIRST_FACE) {
                    // The face's triangulation now matches the requested deflection.
                    face.unset_status(MeshStatus::OUTDATED);
                }
                Ok(())
            }
            Err(error) => {
                if let Ok(face) = model.face_mut(FIRST_FACE) {
                    face.set_status(MeshStatus::FAILURE);
                }
                Err(error)
            }
        }
    }
}

/// Factory providing the Delabella-based triangulation algorithm.
///
/// Source: `BRepMesh_DelabellaMeshAlgoFactory`. OCCT selects a concrete
/// algorithm per surface type; this port returns the single Delabella algorithm
/// for every configuration (the Delabella member of the
/// `IMeshTools_MeshAlgoFactory` family).
pub struct DelabellaMeshAlgoFactory;

impl DelabellaMeshAlgoFactory {
    /// Creates a [`DelabellaBaseMeshAlgo`] for the given parameters.
    pub fn new(_parameters: &MeshParameters) -> Box<dyn BaseMeshAlgo> {
        Box::new(DelabellaBaseMeshAlgo::new())
    }
}

/// Delaunay triangulation of `points` via incremental insertion (split-and-flip).
///
/// A super-triangle circumscribing the bounding box is grown first, then each
/// point is inserted one at a time: the triangle containing it is split (or the
/// two triangles across the edge when it lies on an edge), and the affected
/// edges are flipped (Lawson legalization) until no circumcircle strictly
/// contains another point. Triangles still touching the super-triangle are
/// dropped, leaving the Delaunay triangulation of the point set in CCW order.
///
/// Returns triangle corner triples as indices into `points`. Degenerate inputs
/// (fewer than three points, all coincident or all collinear) yield an empty
/// result.
fn lawson_triangulate(points: &[GpPnt2d]) -> Vec<(usize, usize, usize)> {
    let n = points.len();
    if n < 3 {
        return Vec::new();
    }

    // Bounding box and enclosing super-triangle.
    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    for p in points {
        min_x = min_x.min(p.x());
        min_y = min_y.min(p.y());
        max_x = max_x.max(p.x());
        max_y = max_y.max(p.y());
    }
    let dmax = (max_x - min_x).max(max_y - min_y);
    if dmax <= f64::EPSILON {
        return Vec::new(); // all points coincident
    }
    let xmid = (min_x + max_x) * 0.5;
    let ymid = (min_y + max_y) * 0.5;
    let big = 100.0 * dmax;
    let super_tri: [GpPnt2d; 3] = [
        GpPnt2d::new(xmid - big, ymid - big),
        GpPnt2d::new(xmid + big, ymid - big),
        GpPnt2d::new(xmid, ymid + big),
    ];

    // Local vertices: super-triangle first (0..2), then the input points (3..).
    let mut verts: Vec<GpPnt2d> = super_tri.to_vec();
    verts.extend_from_slice(points);
    let mut tris: Vec<(usize, usize, usize)> = vec![(0, 1, 2)];

    for pi in 0..n {
        let real = 3 + pi;
        let p = verts[real];
        let Some(ti) = tris
            .iter()
            .position(|&(a, b, c)| point_in_triangle(&verts[a], &verts[b], &verts[c], &p))
        else {
            continue; // outside the super-triangle (should not happen)
        };
        let (a, b, c) = tris[ti];

        // When the point falls on an edge, split that edge across its two
        // adjacent triangles instead of splitting one triangle (which would
        // create a degenerate triangle).
        let on_edge = [(a, b), (b, c), (c, a)]
            .iter()
            .copied()
            .find(|&(x, y)| on_segment(&verts[x], &verts[y], &p));

        let mut stack: Vec<(usize, usize, usize)> = Vec::new();
        if let Some((x, y)) = on_edge {
            let z = [a, b, c].iter().copied().find(|&v| v != x && v != y).unwrap();
            if let Some((oti, o)) = find_adjacent(&tris, x, y, z) {
                let t1 = (x, z, real);
                let t2 = (y, z, real);
                let t3 = (x, real, o);
                let t4 = (y, real, o);
                let (lo, hi) = if ti < oti { (ti, oti) } else { (oti, ti) };
                tris[lo] = t1;
                tris[hi] = t2;
                tris.push(t3);
                tris.push(t4);
                stack.extend([t1, t2, t3, t4]);
            } else {
                // On a boundary edge: split the single triangle.
                tris[ti] = (x, z, real);
                tris.push((y, z, real));
                stack.extend([(x, z, real), (y, z, real)]);
            }
        } else {
            // Strictly inside: split into three triangles.
            tris[ti] = (a, b, real);
            tris.push((b, c, real));
            tris.push((c, a, real));
            stack.extend([(a, b, real), (b, c, real), (c, a, real)]);
        }
        legalize(&mut tris, &mut stack, &verts);
    }

    // Drop triangles touching the super-triangle, map to input indices.
    let mut out: Vec<(usize, usize, usize)> = tris
        .into_iter()
        .filter(|&(a, b, c)| a >= 3 && b >= 3 && c >= 3)
        .map(|(a, b, c)| (a - 3, b - 3, c - 3))
        .collect();

    // Drop degenerate (zero-area) triangles and deduplicate.
    out.retain(|&(a, b, c)| signed_area2d(&points[a], &points[b], &points[c]) != 0.0);
    let mut seen = HashSet::new();
    out.retain(|&(a, b, c)| seen.insert((a.min(b).min(c), a.max(b).max(c), a + b + c)));
    out
}

/// Lawson edge legalization: pops a triangle and checks all three of its edges
/// against the opposite vertex of the adjacent triangle, flipping the edge when
/// the circumcircle is violated and re-checking the two replacement triangles.
fn legalize(
    tris: &mut Vec<(usize, usize, usize)>,
    stack: &mut Vec<(usize, usize, usize)>,
    verts: &[GpPnt2d],
) {
    while let Some(t) = stack.pop() {
        let (x, y, z) = t;
        for (e1, e2, opp) in [(x, y, z), (y, z, x), (z, x, y)] {
            if let Some((oti, o)) = find_adjacent(tris, e1, e2, opp) {
                if in_circumcircle(&verts[e1], &verts[e2], &verts[opp], &verts[o]) {
                    // Flip edge (e1,e2): replace (e1,e2,opp) and (e1,e2,o).
                    let ej = tris
                        .iter()
                        .position(|t| tri_matches(t, e1, e2, opp))
                        .expect("flip target triangle present");
                    let n1 = (e1, opp, o);
                    let n2 = (e2, opp, o);
                    let (lo, hi) = if ej < oti { (ej, oti) } else { (oti, ej) };
                    tris[lo] = n1;
                    tris[hi] = n2;
                    stack.push(n1);
                    stack.push(n2);
                    break;
                }
            }
        }
    }
}

/// Finds the triangle sharing edge `(x, y)` whose third vertex is not `z`;
/// returns its index and the opposite vertex. `None` for a boundary edge.
fn find_adjacent(
    tris: &[(usize, usize, usize)],
    x: usize,
    y: usize,
    z: usize,
) -> Option<(usize, usize)> {
    for (ti, &(a, b, c)) in tris.iter().enumerate() {
        let verts = [a, b, c];
        if verts.contains(&x) && verts.contains(&y) && !verts.contains(&z) {
            let opp = verts.iter().copied().find(|&v| v != x && v != y).unwrap();
            return Some((ti, opp));
        }
    }
    None
}

/// Whether the three index triples denote the same triangle (order-insensitive).
fn tri_matches(t: &(usize, usize, usize), x: usize, y: usize, z: usize) -> bool {
    let mut v = [t.0, t.1, t.2];
    v.sort_unstable();
    let mut w = [x, y, z];
    w.sort_unstable();
    v == w
}

/// Twice the signed area of (a, b, c); positive = CCW.
fn signed_area2d(a: &GpPnt2d, b: &GpPnt2d, c: &GpPnt2d) -> f64 {
    (b.x() - a.x()) * (c.y() - a.y()) - (b.y() - a.y()) * (c.x() - a.x())
}

/// Is `p` inside (or on the boundary of) triangle (a, b, c)?
fn point_in_triangle(a: &GpPnt2d, b: &GpPnt2d, c: &GpPnt2d, p: &GpPnt2d) -> bool {
    let d1 = signed_area2d(a, b, p);
    let d2 = signed_area2d(b, c, p);
    let d3 = signed_area2d(c, a, p);
    let has_neg = d1 < 0.0 || d2 < 0.0 || d3 < 0.0;
    let has_pos = d1 > 0.0 || d2 > 0.0 || d3 > 0.0;
    !(has_neg && has_pos)
}

/// Is `p` on the segment from `a` to `b` (collinear and within the extent)?
fn on_segment(a: &GpPnt2d, b: &GpPnt2d, p: &GpPnt2d) -> bool {
    let (dx, dy) = (b.x() - a.x(), b.y() - a.y());
    if dx.abs() < 1e-12 && dy.abs() < 1e-12 {
        return false;
    }
    let cross = dx * (p.y() - a.y()) - dy * (p.x() - a.x());
    if cross.abs() > 1e-9 {
        return false;
    }
    let dot = (p.x() - a.x()) * dx + (p.y() - a.y()) * dy;
    let len2 = dx * dx + dy * dy;
    dot >= -1e-9 && dot <= len2 + 1e-9
}

/// Oriented in-circle determinant of `p` against the triangle (a, b, c).
///
/// The sign is positive when `p` is inside the circumcircle of a CCW triangle,
/// negative outside, zero cocircular; `None` for collinear (degenerate)
/// triangles. Mirrors the `InCircle` predicate of the OCCT Delaunay core.
fn circumcircle_det(a: &GpPnt2d, b: &GpPnt2d, c: &GpPnt2d, p: &GpPnt2d) -> Option<f64> {
    let orient = signed_area2d(a, b, c);
    if orient == 0.0 {
        return None;
    }
    let (ax, ay) = (a.x() - p.x(), a.y() - p.y());
    let (bx, by) = (b.x() - p.x(), b.y() - p.y());
    let (cx, cy) = (c.x() - p.x(), c.y() - p.y());
    let det = (ax * ax + ay * ay) * (bx * cy - by * cx)
        - (bx * bx + by * by) * (ax * cy - ay * cx)
        + (cx * cx + cy * cy) * (ax * by - ay * bx);
    Some(det * orient)
}

/// Fourth power of the max translated coordinate of the four points — the scale
/// used to size the in-circle tolerances (the determinant has units L^4).
fn translated_scale4(a: &GpPnt2d, b: &GpPnt2d, c: &GpPnt2d, p: &GpPnt2d) -> f64 {
    let (ax, ay) = (a.x() - p.x(), a.y() - p.y());
    let (bx, by) = (b.x() - p.x(), b.y() - p.y());
    let (cx, cy) = (c.x() - p.x(), c.y() - p.y());
    let s = ax
        .abs()
        .max(ay.abs())
        .max(bx.abs())
        .max(by.abs())
        .max(cx.abs())
        .max(cy.abs());
    s * s * s * s
}

/// In-circle test used during legalization: is `p` strictly inside the
/// circumcircle of (a, b, c)? A tiny positive margin absorbs float noise so
/// near-cocircular points do not trigger endless flips.
fn in_circumcircle(a: &GpPnt2d, b: &GpPnt2d, c: &GpPnt2d, p: &GpPnt2d) -> bool {
    let Some(value) = circumcircle_det(a, b, c, p) else {
        return false;
    };
    value > 1e-12 * translated_scale4(a, b, c, p)
}

/// Strict in-circle test with a small positive margin: is `p` clearly inside
/// the circumcircle of (a, b, c)? Collinear triangles return false.
fn circumcircle_contains(a: &GpPnt2d, b: &GpPnt2d, c: &GpPnt2d, p: &GpPnt2d) -> bool {
    let Some(value) = circumcircle_det(a, b, c, p) else {
        return false;
    };
    value > 1e-9 * translated_scale4(a, b, c, p)
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use super::*;

    /// Builds an `m x n` grid of Free nodes, returning the structure.
    fn grid_structure(m: usize, n: usize) -> DelaunDataStructure {
        let mut structure = DelaunDataStructure::new(16);
        for i in 0..m {
            for j in 0..n {
                structure.add_node(DelaunVertex::new_parametric(i as f64, j as f64, VertexState::Free));
            }
        }
        structure
    }

    /// Number of boundary vertices of the triangulation — vertices incident to
    /// an edge that appears in exactly one triangle.
    fn boundary_vertices(ds: &DelaunDataStructure) -> usize {
        let mut edge_count: HashMap<(i32, i32), usize> = HashMap::new();
        for &id in ds.elements_of_domain() {
            let verts = ds.get_element(id).vertex_indices;
            for k in 0..3 {
                let (u, v) = (verts[k], verts[(k + 1) % 3]);
                *edge_count.entry((u.min(v), u.max(v))).or_insert(0) += 1;
            }
        }
        let mut boundary: HashSet<i32> = HashSet::new();
        for (&(u, v), &cnt) in &edge_count {
            if cnt == 1 {
                boundary.insert(u);
                boundary.insert(v);
            }
        }
        boundary.len()
    }

    #[test]
    fn grid_triangulates_to_2n_minus_2_minus_h() {
        // 3x3 grid: the Delaunay triangulation keeps the eight perimeter points
        // on the boundary, so h = 8 and 2N - 2 - h = 2*9 - 2 - 8 = 8 triangles.
        let structure = grid_structure(3, 3);
        let nodes_nb = structure.nb_nodes();
        assert_eq!(nodes_nb, 9);

        let mut algo = DelabellaBaseMeshAlgo::new();
        algo.triangulate(structure, &MeshParameters::default()).expect("triangulate grid");

        let ds = algo.result().expect("result stored");
        let triangles = ds.elements_of_domain().len();
        let hull_nb = boundary_vertices(ds);
        assert_eq!(
            triangles,
            2 * nodes_nb - 2 - hull_nb,
            "2N-2-h triangles expected with h={hull_nb}",
            hull_nb = hull_nb
        );
        assert_eq!(hull_nb, 8, "3x3 grid triangulation keeps 8 perimeter boundary vertices");

        // Every triangle references three live, distinct nodes.
        for &id in ds.elements_of_domain() {
            let verts = ds.get_element(id).vertex_indices;
            assert!(verts[0] != verts[1] && verts[1] != verts[2] && verts[0] != verts[2]);
            for &w in &verts {
                assert_ne!(ds.get_node(w).state, VertexState::Deleted);
            }
        }
    }

    #[test]
    fn triangulation_satisfies_empty_circle_property() {
        // A non-symmetric set (no cocircular quadruples) with interior points.
        let pts = [
            (0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), // hull
            (0.6, 0.4), (1.5, 0.7), (0.5, 1.4), (1.3, 1.7), (1.0, 1.0), // interior
        ];
        let mut structure = DelaunDataStructure::new(16);
        let mut nodes: Vec<(i32, GpPnt2d)> = Vec::new();
        for &(u, v) in &pts {
            let idx = structure.add_node(DelaunVertex::new_parametric(u, v, VertexState::Free));
            nodes.push((idx, structure.get_node(idx).location));
        }

        let mut algo = DelabellaBaseMeshAlgo::new();
        algo.triangulate(structure, &MeshParameters::default()).expect("triangulate");

        let ds = algo.result().expect("result stored");
        let hull_nb = boundary_vertices(ds);
        assert_eq!(
            ds.elements_of_domain().len(),
            2 * pts.len() - 2 - hull_nb,
            "h = {hull_nb}",
            hull_nb = hull_nb
        );

        // Empty-circle property: no other point lies strictly inside the
        // circumcircle of any triangle.
        for &id in ds.elements_of_domain() {
            let tri = ds.get_element(id);
            let [va, vb, vc] = tri.vertex_indices;
            let a = ds.get_node(va).location;
            let b = ds.get_node(vb).location;
            let c = ds.get_node(vc).location;
            for &(idx, p) in &nodes {
                if idx == va || idx == vb || idx == vc {
                    continue;
                }
                assert!(
                    !circumcircle_contains(&a, &b, &c, &p),
                    "point {idx} lies strictly inside the circumcircle of triangle {id}",
                    idx = idx,
                    id = id
                );
            }
        }
    }

    #[test]
    fn factory_creates_delabella_algo() {
        let params = MeshParameters::default();
        let mut algo = DelabellaMeshAlgoFactory::new(&params);
        // The factory returns a usable Delabella algorithm: an empty model has
        // no face to triangulate, so `perform` reports the failure cleanly
        // instead of panicking.
        let mut model = MeshModel::default();
        assert!(
            algo.perform(&mut model, &params).is_err(),
            "an empty model has no face to triangulate"
        );
    }

    #[test]
    fn triangulate_errors_with_fewer_than_three_nodes() {
        let mut structure = DelaunDataStructure::new(16);
        structure.add_node(DelaunVertex::new_parametric(0.0, 0.0, VertexState::Free));
        structure.add_node(DelaunVertex::new_parametric(1.0, 0.0, VertexState::Free));

        let mut algo = DelabellaBaseMeshAlgo::new();
        let error = algo
            .triangulate(structure, &MeshParameters::default())
            .expect_err("fewer than 3 nodes must fail");
        assert!(error.contains("fewer than 3"), "unexpected error: {error}", error = error);
        assert!(algo.result().is_none());
    }

    #[test]
    fn square_triangulates_to_two_triangles() {
        // Unit square, four hull vertices => 2N - 2 - h = 8 - 2 - 4 = 2 triangles.
        let mut structure = DelaunDataStructure::new(16);
        for &(u, v) in &[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)] {
            structure.add_node(DelaunVertex::new_parametric(u, v, VertexState::Free));
        }

        let mut algo = DelabellaBaseMeshAlgo::new();
        algo.triangulate(structure, &MeshParameters::default()).expect("triangulate square");

        let ds = algo.result().expect("result stored");
        assert_eq!(ds.elements_of_domain().len(), 2);
        // The square area is covered exactly: total triangle area == 1.
        let total: f64 = ds
            .elements_of_domain()
            .iter()
            .map(|&id| {
                let tri = ds.get_element(id);
                let [va, vb, vc] = tri.vertex_indices;
                signed_area2d(
                    &ds.get_node(va).location,
                    &ds.get_node(vb).location,
                    &ds.get_node(vc).location,
                )
                .abs()
                    * 0.5
            })
            .sum();
        assert!((total - 1.0).abs() < 1e-12, "area = {total}", total = total);
    }
}
