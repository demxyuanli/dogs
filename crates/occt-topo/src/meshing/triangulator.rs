//! Port of OCCT triangulator — Wave 3 BRepMesh.
//!
//! `Triangulator` is the high-level face-meshing orchestrator. For each face it
//! takes the boundary UV points — either from a discrete model's edge pcurves
//! ([`MeshTool::extract_model_face`]) or from a UV polygon
//! ([`MeshTool::extract_face`]) — runs the 2D Delaunay triangulation
//! ([`Delaun`]), refines interior triangles whose deviation from the surface
//! exceeds the deflection, and evaluates the surface at every UV node to build a
//! 3D triangle mesh.
//!
//! Sources:
//! - `BRepMesh_Triangulator.hxx/.cxx` (per-face triangulation entry).
//! - `BRepMesh_DelaunayBaseMeshAlgo.hxx` / `BRepMesh_DelaunayDeflectionControlMeshAlgo.hxx`
//!   (per-face orchestration + deflection-refinement role of the modern
//!   pipeline; `mesh_algo.rs` / `deflection_control.rs` / `node_insertion.rs`
//!   are Wave-3 stubs, so this module plays that role for the mesh algorithms
//!   below).

use std::collections::{HashMap, HashSet};

use occt_core::gp::{GpPnt, GpPnt2d};
use occt_core::poly::triangulation::Triangle;
use occt_geom::Surface;

use super::data_model::MeshModel;
use super::delaun::Delaun;
use super::delaun_types::{DelaunVertex, VertexState};
use super::edge_discret::MeshFace as UvFace;
use super::mesh_tool::{FaceMeshData, MeshTool};
use super::parameters::MeshParameters;

/// Result of triangulating one face: 3D nodes, triangles and the UV of every node.
#[derive(Debug, Clone)]
pub struct FaceTriangulation {
    /// Model index of the triangulated face.
    pub face_index: usize,
    /// 3D mesh nodes (surface evaluated at the UV nodes).
    pub vertices: Vec<GpPnt>,
    /// Triangles referencing `vertices` by 0-based index.
    pub triangles: Vec<Triangle>,
    /// UV coordinates of every node (`vertices[i] == surface.d0(uv[i])`).
    pub uv: Vec<GpPnt2d>,
}

/// Maximum number of deflection-refinement passes.
const MAX_REFINE_PASSES: usize = 4;

/// High-level triangulation orchestrator (port of `BRepMesh_Triangulator`).
///
/// Holds the [`MeshParameters`] in effect and exposes the per-face Delaunay
/// triangulation pipeline: boundary UV points → `Delaun` → deflection
/// refinement → triangle mesh.
pub struct Triangulator {
    params: MeshParameters,
}

impl Triangulator {
    /// Builds a triangulator under the given meshing parameters.
    pub fn new(params: MeshParameters) -> Self {
        Self { params }
    }

    /// The parameters in effect.
    pub fn parameters(&self) -> &MeshParameters {
        &self.params
    }

    /// Triangulates a face given its UV vertices + constraint edges onto `surface`.
    ///
    /// Boundary vertices (the endpoints of `data.constraint_edges`) are tagged
    /// `Frontier`; for a convex face the Delaunay triangulation preserves the
    /// outer wire as hull edges. Interior triangles whose deviation from the
    /// surface exceeds the face deflection are refined by inserting their UV
    /// barycenters and re-triangulating.
    pub fn triangulate(
        &self,
        surface: &dyn Surface,
        data: &FaceMeshData,
    ) -> Result<FaceTriangulation, String> {
        if data.vertices.is_empty() {
            return Err("Triangulator::triangulate: no UV vertices supplied".to_string());
        }
        if data.constraint_edges.len() < 3 {
            return Err("Triangulator::triangulate: too few constraint edges".to_string());
        }

        let deflection = self.face_deflection();
        let mut interior: Vec<GpPnt2d> = Vec::new();
        let mut delaun = self.run_delaun(data, &interior)?;
        for _ in 0..MAX_REFINE_PASSES {
            let deviating = self.deviating_centers(&delaun, surface, deflection);
            if deviating.is_empty() {
                break;
            }
            interior.extend(deviating);
            delaun = self.run_delaun(data, &interior)?;
        }
        self.build_triangulation(&delaun, surface, 0)
    }

    /// Triangulates a UV-polygon face (convenience over
    /// [`MeshTool::extract_face`] + [`Self::triangulate`]).
    pub fn triangulate_face_polygon(
        &self,
        surface: &dyn Surface,
        face: &UvFace,
        tolerance: f64,
    ) -> Result<FaceTriangulation, String> {
        let data = MeshTool::new(self.params.clone())
            .extract_face(face, tolerance)
            .map_err(|e| format!("Triangulator::triangulate_face_polygon: {e}"))?;
        let mut tri = self.triangulate(surface, &data)?;
        tri.face_index = face.id;
        Ok(tri)
    }

    /// High-level entry: triangulates every face of a discrete model.
    ///
    /// The boundary UV points are taken from the edge pcurves of each face
    /// (`MeshTool::extract_model_face`), so the edge discretization step of the
    /// pipeline must have populated them first. Faces with no surface or no UV
    /// data yield an error.
    pub fn triangulate_model(&self, model: &MeshModel) -> Result<Vec<FaceTriangulation>, String> {
        let tool = MeshTool::new(self.params.clone());
        let mut out: Vec<FaceTriangulation> = Vec::with_capacity(model.faces_nb());
        for i in 0..model.faces_nb() {
            let face = model.face(i)?;
            let surface = face.surface().ok_or_else(|| {
                format!("Triangulator::triangulate_model: face {i} has no surface")
            })?;
            let data = tool
                .extract_model_face(model, i, 1e-6)
                .map_err(|e| format!("Triangulator::triangulate_model: {e}"))?;
            let tri = self.triangulate(surface.as_ref(), &data)?;
            out.push(FaceTriangulation { face_index: i, ..tri });
        }
        if out.is_empty() {
            return Err("Triangulator::triangulate_model: model has no faces".to_string());
        }
        Ok(out)
    }

    /// Deflection driving the interior refinement: the interior deflection when
    /// set, otherwise the boundary deflection.
    fn face_deflection(&self) -> f64 {
        if self.params.deflection_interior > 0.0 {
            self.params.deflection_interior
        } else {
            self.params.deflection
        }
    }

    /// Builds the Delaunay input over the face's boundary UV vertices plus the
    /// given interior refinement points and runs the triangulation.
    ///
    /// Boundary vertices (those on a constraint edge) are tagged `Frontier`;
    /// interior points are `Free`. The Delaunay triangulation of a point set
    /// keeps the convex hull edges, so the outer wire of a convex face is
    /// preserved without explicit constraint links.
    ///
    /// ponytail: constraint edges are not registered as `Frontier` links (the
    /// `Delaun::new_with_data` pre-linked path currently produces a gapped
    /// triangulation for a plain square). For non-convex faces / holes the
    /// boundary must be enforced with links — revisit when a curved/hole face
    /// needs it.
    fn run_delaun(&self, data: &FaceMeshData, interior: &[GpPnt2d]) -> Result<Delaun, String> {
        let boundary: HashSet<usize> =
            data.constraint_edges.iter().flat_map(|&(a, b)| [a, b]).collect();

        let mut vertices: Vec<DelaunVertex> = Vec::with_capacity(data.vertices.len() + interior.len());
        for (i, v) in data.vertices.iter().enumerate() {
            let state = if boundary.contains(&i) {
                VertexState::Frontier
            } else {
                VertexState::Free
            };
            vertices.push(DelaunVertex::new(
                GpPnt2d::new(v.u, v.v),
                GpPnt::zero(),
                i as i32,
                state,
            ));
        }
        for &uv in interior {
            vertices.push(DelaunVertex::new(uv, GpPnt::zero(), 0, VertexState::Free));
        }

        Ok(Delaun::new_vertices(&vertices))
    }

    /// UV barycenters of the triangles whose deviation from the surface exceeds
    /// `deflection`. Deviation is measured at the triangle UV barycenter against
    /// the linearly-interpolated corner values.
    ///
    /// ponytail: minimal local stand-in for the empty `deflection_control` stub
    /// (OCCT `DelaunayDeflectionControlMeshAlgo` iteratively splits long edges).
    /// Inserting the triangle barycenter and re-triangulating is equivalent for
    /// the planar/low-curvature cases this port targets; edge-midpoint splitting
    /// can be added when curved surfaces need it.
    fn deviating_centers(
        &self,
        delaun: &Delaun,
        surface: &dyn Surface,
        deflection: f64,
    ) -> Vec<GpPnt2d> {
        let ds = delaun.result();
        let ids: Vec<i32> = ds.elements_of_domain().iter().copied().collect();
        let mut out: Vec<GpPnt2d> = Vec::new();
        for &id in &ids {
            if ds.element_movability(id) == VertexState::Deleted {
                continue;
            }
            let nodes = ds.element_nodes(&ds.get_element(id));
            let [p0, p1, p2] = [
                ds.get_node(nodes[0]).location,
                ds.get_node(nodes[1]).location,
                ds.get_node(nodes[2]).location,
            ];
            let c = GpPnt2d::new(
                (p0.x() + p1.x() + p2.x()) / 3.0,
                (p0.y() + p1.y() + p2.y()) / 3.0,
            );
            let a = surface.d0(p0.x(), p0.y());
            let b = surface.d0(p1.x(), p1.y());
            let cc = surface.d0(p2.x(), p2.y());
            let interp = GpPnt::new(
                (a.x() + b.x() + cc.x()) / 3.0,
                (a.y() + b.y() + cc.y()) / 3.0,
                (a.z() + b.z() + cc.z()) / 3.0,
            );
            if surface.d0(c.x(), c.y()).distance(&interp) > deflection {
                out.push(c);
            }
        }
        out
    }

    /// Extracts the triangle mesh from a finished Delaunay structure: every live
    /// element becomes a [`Triangle`], each node index maps to a UV point whose
    /// 3D image on `surface` is the mesh vertex.
    fn build_triangulation(
        &self,
        delaun: &Delaun,
        surface: &dyn Surface,
        face_index: usize,
    ) -> Result<FaceTriangulation, String> {
        let ds = delaun.result();
        let mut node_to_vertex: HashMap<i32, usize> = HashMap::new();
        let mut vertices: Vec<GpPnt> = Vec::new();
        let mut uv: Vec<GpPnt2d> = Vec::new();
        let mut triangles: Vec<Triangle> = Vec::new();

        let ids: Vec<i32> = ds.elements_of_domain().iter().copied().collect();
        for &id in &ids {
            if ds.element_movability(id) == VertexState::Deleted {
                continue;
            }
            let nodes = ds.element_nodes(&ds.get_element(id));
            let mut corner = [0usize; 3];
            for (k, &n) in nodes.iter().enumerate() {
                let vi = *node_to_vertex.entry(n).or_insert_with(|| {
                    let loc = ds.get_node(n).location;
                    uv.push(loc);
                    vertices.push(surface.d0(loc.x(), loc.y()));
                    vertices.len() - 1
                });
                corner[k] = vi;
            }
            triangles.push(Triangle::new(corner[0], corner[1], corner[2]));
        }

        if triangles.is_empty() {
            return Err("Triangulator::triangulate: Delaunay produced no triangles".to_string());
        }
        Ok(FaceTriangulation { face_index, vertices, triangles, uv })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::meshing::model_builder::ModelBuilder;
    use crate::primitives::BRepPrimBox;
    use crate::tgeometry::GeometryRegistry;
    use occt_core::gp::{GpAx3, GpPln};
    use occt_geom::GeomPlane;

    fn box_uv_face() -> UvFace {
        UvFace::new(vec![
            GpPnt2d::new(0.0, 0.0),
            GpPnt2d::new(1.0, 0.0),
            GpPnt2d::new(1.0, 1.0),
            GpPnt2d::new(0.0, 1.0),
            GpPnt2d::new(0.0, 0.0),
        ])
    }

    #[test]
    fn box_face_triangulates_to_uv_triangles() {
        let params = MeshParameters::default();
        let surface = GeomPlane::new(GpPln::new(GpAx3::standard()));
        let tri = Triangulator::new(params)
            .triangulate_face_polygon(&surface, &box_uv_face(), 1e-6)
            .expect("box face triangulated");

        // A unit square triangulates into two UV triangles.
        assert!(tri.triangles.len() >= 2, "triangles {}", tri.triangles.len());
        assert_eq!(tri.uv.len(), 4, "four UV corners");
        // The 3D vertices lie on the z=0 plane at the UV coordinates.
        for (p, uv) in tri.vertices.iter().zip(&tri.uv) {
            assert!((p.z() - 0.0).abs() < 1e-9, "vertex not on plane: {p:?}");
            assert!(p.distance(&GpPnt::new(uv.x(), uv.y(), 0.0)) < 1e-9, "UV mismatch {p:?} vs {uv:?}");
        }
        // The UV triangles cover the unit square (total area 1.0).
        let area: f64 = tri
            .triangles
            .iter()
            .map(|t| {
                let a = &tri.uv[t.n0];
                let b = &tri.uv[t.n1];
                let c = &tri.uv[t.n2];
                0.5 * ((b.x() - a.x()) * (c.y() - a.y()) - (c.x() - a.x()) * (b.y() - a.y()))
                    .abs()
            })
            .sum();
        assert!((area - 1.0).abs() < 1e-9, "UV area {area}");
    }

    #[test]
    fn box_face_keeps_boundary_edges() {
        // The four frontier constraint edges of the square must appear in the
        // emitted mesh (Delaunay must not flip a boundary edge).
        let params = MeshParameters::default();
        let surface = GeomPlane::new(GpPln::new(GpAx3::standard()));
        let tri = Triangulator::new(params)
            .triangulate_face_polygon(&surface, &box_uv_face(), 1e-6)
            .expect("box face triangulated");

        let mut edges: Vec<(usize, usize)> = tri
            .triangles
            .iter()
            .flat_map(|t| [(t.n0, t.n1), (t.n1, t.n2), (t.n2, t.n0)])
            .map(|(a, b)| (a.min(b), a.max(b)))
            .collect();
        edges.sort_unstable();
        edges.dedup();
        assert!(
            edges.contains(&(0, 1)) && edges.contains(&(1, 2)) && edges.contains(&(2, 3)) && edges.contains(&(0, 3)),
            "boundary edges missing: {edges:?}"
        );
    }

    #[test]
    fn triangulate_model_reports_unpopulated_pcurves() {
        // A fresh model has empty pcurves until the edge discretizers run; the
        // model-level orchestrator must surface that instead of panicking.
        let shape = BRepPrimBox::make_box(1.0, 1.0, 1.0).solid.0.clone();
        let params = MeshParameters::default();
        let model = ModelBuilder::build_model(&shape, &params).expect("model built");
        let err = Triangulator::new(params).triangulate_model(&model);
        assert!(err.is_err(), "unpopulated pcurves must yield an error");
        GeometryRegistry::global().clear_shape(&shape);
    }
}
