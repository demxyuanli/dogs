//! Port of OCCT fast_discret — Wave 3 BRepMesh.
//!
//! `FastDiscret` is the legacy fast-discretization algorithm: an edge becomes a
//! deflection-adaptive polyline ([`CurveTessellator`], reusing
//! `edge_discret::CurveTessellator`), a face becomes a UV grid + Delaunay
//! triangle soup ([`FaceDiscret`] + [`Delaun`]), followed by a 3D
//! vertex-welding pass.
//!
//! Source: `BRepMesh_FastDiscret.hxx` (deprecated in modern OCCT in favour of
//! `IMeshTools`; the algorithm contract is reconstructed from the modern
//! `BRepMesh_EdgeDiscret` / `BRepMesh_FaceDiscret` split).

use std::collections::HashMap;
use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_core::poly::triangulation::Triangle;
use occt_geom::{Curve, Surface};

use super::delaun::Delaun;
use super::delaun_types::{DelaunVertex, VertexState};
use super::edge_discret::{CurveTessellator, MeshFace};
use super::face_discret::FaceDiscret;
use super::parameters::MeshParameters;

/// Legacy fast discretizer — edge → deflection-adaptive polyline, face → UV
/// grid + Delaunay triangle soup, with 3D vertex welding.
#[derive(Debug, Clone)]
pub struct FastDiscret {
    params: MeshParameters,
}

impl FastDiscret {
    /// Builds a discretizer under the given meshing parameters.
    pub fn new(params: MeshParameters) -> Self {
        Self { params }
    }

    /// The parameters in effect.
    pub fn parameters(&self) -> &MeshParameters {
        &self.params
    }

    /// Discretizes a 3D curve into a deflection-adaptive polyline.
    ///
    /// The curve is seeded uniformly and every chord refined while it deviates
    /// from the curve by more than the parameters' deflection; both endpoints
    /// are always present. Mirrors `BRepMesh_CurveTessellator` /
    /// `GCPnts_UniformDeflection`.
    pub fn discretize_edge(&self, curve: &Arc<dyn Curve>) -> Vec<GpPnt> {
        let (a, b) = (curve.first_parameter(), curve.last_parameter());
        if !(a.is_finite() && b.is_finite() && b > a) {
            return Vec::new();
        }
        let def = self.params.deflection.max(1e-12);
        CurveTessellator::from_range(curve.clone(), a, b, def, 2)
            .points()
            .to_vec()
    }

    /// Discretizes a UV-polygon face into 3D vertices + triangles.
    ///
    /// The boundary UV points (kept verbatim) are combined with a
    /// deflection-bounded interior UV grid ([`FaceDiscret::discretize_face`]),
    /// triangulated by [`Delaun`], and lifted onto `surface`; coincident 3D
    /// vertices are welded.
    pub fn discretize_face(
        &self,
        face: &MeshFace,
        surface: &dyn Surface,
    ) -> (Vec<GpPnt>, Vec<Triangle>) {
        let uv_pts = FaceDiscret::new(self.params.clone()).discretize_face(face);
        if uv_pts.len() < 3 {
            return (Vec::new(), Vec::new());
        }

        let vertices: Vec<DelaunVertex> = uv_pts
            .iter()
            .map(|p| DelaunVertex::new(*p, GpPnt::zero(), 0, VertexState::Free))
            .collect();
        let delaun = Delaun::new_vertices(&vertices);
        let ds = delaun.result();

        let mut out_v: Vec<GpPnt> = Vec::new();
        let mut node_to_v: HashMap<i32, usize> = HashMap::new();
        let mut out_t: Vec<Triangle> = Vec::new();
        for &id in ds.elements_of_domain() {
            if ds.element_movability(id) == VertexState::Deleted {
                continue;
            }
            let nodes = ds.element_nodes(&ds.get_element(id));
            let mut corner = [0usize; 3];
            for (k, &n) in nodes.iter().enumerate() {
                let vi = *node_to_v.entry(n).or_insert_with(|| {
                    let loc = ds.get_node(n).location;
                    out_v.push(surface.d0(loc.x(), loc.y()));
                    out_v.len() - 1
                });
                corner[k] = vi;
            }
            out_t.push(Triangle::new(corner[0], corner[1], corner[2]));
        }

        if out_t.is_empty() {
            return (Vec::new(), Vec::new());
        }
        let weld_tol = (self.params.deflection * 0.01).max(1e-12);
        Self::weld_vertices(out_v, &out_t, weld_tol)
    }

    /// Welds coincident 3D vertices within `tolerance` into a unique set and
    /// remaps the triangle indices to the surviving vertices.
    ///
    /// ponytail: O(n²) scan; switch to a spatial hash when vertex counts grow.
    pub fn weld_vertices(
        vertices: Vec<GpPnt>,
        triangles: &[Triangle],
        tolerance: f64,
    ) -> (Vec<GpPnt>, Vec<Triangle>) {
        let mut unique: Vec<GpPnt> = Vec::with_capacity(vertices.len());
        let mut remap: Vec<usize> = Vec::with_capacity(vertices.len());
        for p in &vertices {
            let idx = unique
                .iter()
                .position(|q| q.distance(p) <= tolerance)
                .unwrap_or_else(|| {
                    unique.push(*p);
                    unique.len() - 1
                });
            remap.push(idx);
        }
        let tris = triangles
            .iter()
            .map(|t| Triangle::new(remap[t.n0], remap[t.n1], remap[t.n2]))
            .collect();
        (unique, tris)
    }
}

/// Maximum distance of a polyline point from the chord of its segment, sampled
/// at 8 interior curve parameters per segment. Used by tests to verify the
/// deflection target.
#[cfg(test)]
fn max_chord_error(curve: &dyn Curve, points: &[GpPnt]) -> f64 {
    let n = points.len();
    if n < 2 {
        return 0.0;
    }
    let (a, b) = (curve.first_parameter(), curve.last_parameter());
    if !(a.is_finite() && b.is_finite() && b > a) {
        return 0.0;
    }
    // Chord-length parametrization maps each polyline vertex to a curve parameter.
    let mut s = Vec::with_capacity(n);
    s.push(0.0);
    for i in 1..n {
        s.push(s[i - 1] + points[i].distance(&points[i - 1]));
    }
    let total = *s.last().unwrap();
    if total <= 1e-12 {
        return 0.0;
    }

    let mut max = 0.0f64;
    for i in 0..n - 1 {
        let ua = a + (b - a) * s[i] / total;
        let ub = a + (b - a) * s[i + 1] / total;
        for k in 1..8 {
            let u = ua + (ub - ua) * k as f64 / 8.0;
            let d = point_segment_dist_3d(&curve.d0(u), &points[i], &points[i + 1]);
            if d > max {
                max = d;
            }
        }
    }
    max
}

/// Distance from `p` to the segment `a..b`.
#[cfg(test)]
fn point_segment_dist_3d(p: &GpPnt, a: &GpPnt, b: &GpPnt) -> f64 {
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

#[cfg(test)]
mod tests {
    use super::*;
    use occt_core::gp::{GpAx2, GpAx3, GpCirc, GpPln, GpPnt2d};
    use occt_geom::{GeomCircle, GeomPlane};

    fn box_uv_face() -> MeshFace {
        MeshFace::new(vec![
            GpPnt2d::new(0.0, 0.0),
            GpPnt2d::new(1.0, 0.0),
            GpPnt2d::new(1.0, 1.0),
            GpPnt2d::new(0.0, 1.0),
            GpPnt2d::new(0.0, 0.0),
        ])
    }

    #[test]
    fn circle_edge_chord_deviation_within_target() {
        let radius = 2.0;
        let circ = GeomCircle::new(GpCirc::new(GpAx2::standard(), radius));
        let curve: Arc<dyn Curve> = Arc::new(circ);
        let params = MeshParameters { deflection: 0.05, ..MeshParameters::default() };
        let pts = FastDiscret::new(params).discretize_edge(&curve);

        assert!(pts.len() >= 2, "polyline length {}", pts.len());
        // Endpoints of a full circle coincide.
        assert!(pts.first().unwrap().distance(pts.last().unwrap()) < 1e-9);
        // Every chord deviates from the true circle by at most the deflection.
        let err = max_chord_error(curve.as_ref(), &pts);
        assert!(err <= 0.05 + 1e-6, "chord deviation {err} exceeds target 0.05");
    }

    #[test]
    fn finer_deflection_gives_tighter_chord() {
        let circ = GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0));
        let curve: Arc<dyn Curve> = Arc::new(circ);
        let coarse = FastDiscret::new(MeshParameters { deflection: 0.2, ..MeshParameters::default() })
            .discretize_edge(&curve);
        let fine = FastDiscret::new(MeshParameters { deflection: 0.01, ..MeshParameters::default() })
            .discretize_edge(&curve);
        assert!(fine.len() >= coarse.len(), "finer deflection must not lose points");
        let err_fine = max_chord_error(curve.as_ref(), &fine);
        assert!(err_fine <= 0.01 + 1e-6, "fine chord deviation {err_fine}");
    }

    #[test]
    fn face_grid_produces_vertices_and_triangles() {
        let face = box_uv_face();
        let surface = GeomPlane::new(GpPln::new(GpAx3::standard()));
        let params = MeshParameters { deflection: 0.25, ..MeshParameters::default() };
        let (verts, tris) = FastDiscret::new(params).discretize_face(&face, &surface);

        assert!(verts.len() > 0, "vertex count {}", verts.len());
        assert!(tris.len() > 0, "triangle count {}", tris.len());
        // All vertices on the z=0 plane.
        assert!(verts.iter().all(|p| p.z().abs() < 1e-9));
        // Triangle indices stay in range.
        for t in &tris {
            assert!(t.n0 < verts.len() && t.n1 < verts.len() && t.n2 < verts.len());
        }
    }

    #[test]
    fn weld_vertices_collapses_coincident_points() {
        let v = vec![
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(1e-10, 0.0, 0.0), // coincides with vertex 0
        ];
        let tris = vec![Triangle::new(0, 1, 2)];
        let (uv, ut) = FastDiscret::weld_vertices(v, &tris, 1e-6);
        assert_eq!(uv.len(), 2, "two unique vertices after welding");
        assert_eq!(ut[0].n0, ut[0].n2, "coincident vertices remapped to one index");
        assert!(ut[0].n0 < 2 && ut[0].n1 < 2 && ut[0].n2 < 2);
    }
}
