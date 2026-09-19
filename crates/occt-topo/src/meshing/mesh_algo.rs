//! Port of OCCT mesh_algo — Wave 3 BRepMesh.
//!
//! Base classes of the per-face triangulation hierarchy:
//!
//! * [`BaseMeshAlgo`] — the `IMeshTools_MeshAlgo` interface (trait).
//! * [`ConstrainedBaseMeshAlgo`] — adds the Delaunay cells-count and
//!   post-process hooks.
//! * [`CustomBaseMeshAlgo`] — custom-triangulation flow
//!   (`buildBaseTriangulation` + the shared [`perform_custom_mesh`] driver).
//! * [`CustomDelaunayBaseMeshAlgo`] — re-initializes the Delaunay circle cell
//!   filter before delegating to the wrapped base algorithm.
//! * [`DelaunayBaseMeshAlgo`] — the classic Watson Delaunay triangulation:
//!   registers the face boundary as `Frontier`/`Fixed` constraint edges and
//!   holds the [`Delaun`] mesh built over the UV point set.
//!
//! Sources (`src/ModelingAlgorithms/TKMesh/BRepMesh/`): `BRepMesh_BaseMeshAlgo`,
//! `BRepMesh_ConstrainedBaseMeshAlgo`, `BRepMesh_CustomBaseMeshAlgo`,
//! `BRepMesh_CustomDelaunayBaseMeshAlgo`, `BRepMesh_DelaunayBaseMeshAlgo`.

use occt_core::gp::GpPnt;

use super::data_model::{MeshModel, MeshStatus};
use super::delaun::Delaun;
use super::delaun_data::DelaunDataStructure;
use super::delaun_types::{DelaunVertex, VertexState};
use super::parameters::MeshParameters;

/// Base interface for algorithms building face triangulation.
///
/// Source: `BRepMesh_BaseMeshAlgo`, which implements `IMeshTools_MeshAlgo::Perform`.
/// OCCT performs one discrete face; this port exposes the algorithm at the
/// discrete-model level so the whole hierarchy shares the same entry point.
pub trait BaseMeshAlgo {
    /// Meshes the given discrete model under the given parameters.
    fn perform(&mut self, model: &mut MeshModel, parameters: &MeshParameters) -> Result<(), String>;
}

/// Base for Delaunay-approach algorithms building face triangulation.
///
/// Source: `BRepMesh_ConstrainedBaseMeshAlgo`. Adds the two hooks the concrete
/// algorithms use to tune the Delaunay run: the acceleration-grid cell count and
/// a post-processing callback.
pub trait ConstrainedBaseMeshAlgo: BaseMeshAlgo {
    /// Size of the cell to be used by the acceleration circles grid structure.
    /// `(-1, -1)` means "let the mesher pick the default".
    fn get_cells_count(&self, _vertices_nb: usize) -> (i32, i32) {
        (-1, -1)
    }

    /// Performs processing of the generated mesh. By default does nothing;
    /// expected to be called at the end of a mesher run.
    fn post_process_mesh(&mut self, _mesher: &mut Delaun, _parameters: &MeshParameters) {}
}

/// Base for custom-triangulation algorithms.
///
/// Source: `BRepMesh_CustomBaseMeshAlgo`. A custom algorithm implements
/// [`CustomBaseMeshAlgo::build_base_triangulation`] to seed the structure with
/// its own triangulation; the shared [`perform_custom_mesh`] driver then runs
/// the Delaunay constraint insertion and hands the mesher to
/// [`ConstrainedBaseMeshAlgo::post_process_mesh`].
pub trait CustomBaseMeshAlgo: ConstrainedBaseMeshAlgo {
    /// Builds the base triangulation using the custom algorithm.
    fn build_base_triangulation(
        &mut self,
        structure: &mut DelaunDataStructure,
    ) -> Result<(), String>;
}

/// Wraps a [`CustomBaseMeshAlgo`] and re-initializes the Delaunay circle cell
/// filter before the base post-process.
///
/// Source: the `BRepMesh_CustomDelaunayBaseMeshAlgo<BaseAlgo>` template. The
/// wrapper owns the base algorithm and adds the `InitCirclesTool` step that the
/// deflection-control refinement needs.
pub struct CustomDelaunayBaseMeshAlgo<A: CustomBaseMeshAlgo> {
    inner: A,
}

impl<A: CustomBaseMeshAlgo> CustomDelaunayBaseMeshAlgo<A> {
    /// Wraps a custom triangulation algorithm.
    pub fn new(inner: A) -> Self {
        Self { inner }
    }
}

impl<A: CustomBaseMeshAlgo> BaseMeshAlgo for CustomDelaunayBaseMeshAlgo<A> {
    fn perform(&mut self, model: &mut MeshModel, parameters: &MeshParameters) -> Result<(), String> {
        perform_custom_mesh(self, model, parameters)
    }
}

impl<A: CustomBaseMeshAlgo> ConstrainedBaseMeshAlgo for CustomDelaunayBaseMeshAlgo<A> {
    fn get_cells_count(&self, vertices_nb: usize) -> (i32, i32) {
        self.inner.get_cells_count(vertices_nb)
    }

    fn post_process_mesh(&mut self, mesher: &mut Delaun, parameters: &MeshParameters) {
        let (cells_u, cells_v) = self.inner.get_cells_count(mesher.result().nb_nodes());
        if cells_u > 0 && cells_v > 0 {
            mesher.init_circles_tool_public(cells_u, cells_v);
        }
        self.inner.post_process_mesh(mesher, parameters);
    }
}

impl<A: CustomBaseMeshAlgo> CustomBaseMeshAlgo for CustomDelaunayBaseMeshAlgo<A> {
    fn build_base_triangulation(
        &mut self,
        structure: &mut DelaunDataStructure,
    ) -> Result<(), String> {
        self.inner.build_base_triangulation(structure)
    }
}

/// Classic Watson Delaunay triangulation of a face.
///
/// Source: `BRepMesh_DelaunayBaseMeshAlgo`. Builds the Delaunay data structure
/// from the model face's UV pcurves (boundary as `Frontier` / `Fixed` constraint
/// links), runs the [`Delaun`] triangulator and keeps the resulting mesh for
/// querying through [`DelaunayBaseMeshAlgo::mesher`].
pub struct DelaunayBaseMeshAlgo {
    mesher: Option<Delaun>,
}

impl DelaunayBaseMeshAlgo {
    /// Creates an idle algorithm.
    pub fn new() -> Self {
        Self { mesher: None }
    }

    /// The Delaunay mesh produced by the last successful run, if any.
    pub fn mesher(&self) -> Option<&Delaun> {
        self.mesher.as_ref()
    }

    /// Triangulates a pre-built Delaunay data structure.
    ///
    /// Registers every stored node with the mesher, runs the Delaunay insertion
    /// (which enforces the pre-registered `Frontier`/`Fixed` constraint links)
    /// and stores the resulting mesh. Fails when the structure has fewer than
    /// three nodes.
    pub fn triangulate(
        &mut self,
        structure: DelaunDataStructure,
        parameters: &MeshParameters,
    ) -> Result<(), String> {
        let nodes_nb = structure.nb_nodes();
        if nodes_nb < 3 {
            return Err("DelaunayBaseMeshAlgo::triangulate: fewer than 3 nodes".to_string());
        }
        let mut vertex_order: Vec<i32> = (1..=nodes_nb as i32).collect();
        let mut mesher = Delaun::new_with_data(structure, &mut vertex_order);
        self.post_process_mesh(&mut mesher, parameters);
        self.mesher = Some(mesher);
        Ok(())
    }
}

impl Default for DelaunayBaseMeshAlgo {
    fn default() -> Self {
        Self::new()
    }
}

impl BaseMeshAlgo for DelaunayBaseMeshAlgo {
    fn perform(&mut self, model: &mut MeshModel, parameters: &MeshParameters) -> Result<(), String> {
        // OCCT invokes the base algorithm once per face. The model-level port
        // triangulates the first face; the produced Delaunay mesh is queried
        // through `mesher()`.
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

impl ConstrainedBaseMeshAlgo for DelaunayBaseMeshAlgo {}

/// Registers the UV boundary of a model face in a fresh Delaunay data structure.
///
/// Port of `BRepMesh_BaseMeshAlgo::initDataStructure`. Walks the face's wires →
/// edges → pcurves, adding one node per pcurve point and a link per consecutive
/// pair. The outer wire's links are `Frontier`; inner (hole) wires' links are
/// `Fixed` so holes of zero area survive the constraint processing.
pub fn build_data_structure(
    model: &MeshModel,
    face_index: usize,
) -> Result<DelaunDataStructure, String> {
    let face = model.face(face_index)?;
    let mut structure = DelaunDataStructure::new(16);
    for (wire_pos, &wire_index) in face.wires().iter().enumerate() {
        let wire = model.wire(wire_index)?;
        let link_state = if wire_pos == 0 {
            VertexState::Frontier
        } else {
            VertexState::Fixed
        };
        for edge_pos in 0..wire.edges_nb() {
            let edge_index = wire.edge(edge_pos)?;
            let orientation = wire.edge_orientation(edge_pos)?;
            let edge = model.edge(edge_index)?;
            let Some(pcurve) = edge.pcurve_for(face_index, orientation) else {
                continue;
            };
            let mut previous = -1i32;
            // ponytail: the 3D point (`p3d`) and the location index are not
            // derived from the curve discretization — the triangulation only
            // consumes the UV location, so p3d stays zero until the node-mapping
            // stage (OCCT `collectNodes`) lands.
            for (k, &point) in pcurve.points().iter().enumerate() {
                let node = DelaunVertex::new(point, GpPnt::zero(), k as i32, link_state);
                let node_index = structure.add_node(node);
                if previous != -1 && previous != node_index {
                    structure.add_link(previous, node_index, link_state);
                }
                previous = node_index;
            }
        }
    }
    Ok(structure)
}

/// Shared driver of the custom-triangulation flow.
///
/// Source: `BRepMesh_CustomBaseMeshAlgo::generateMesh` (and the
/// `InitCirclesTool` step of `BRepMesh_CustomDelaunayBaseMeshAlgo`).
/// Builds the structure from the first model face, lets the custom algorithm
/// seed its base triangulation, runs the Delaunay constraint insertion and hands
/// the finished mesher to [`ConstrainedBaseMeshAlgo::post_process_mesh`].
pub fn perform_custom_mesh<A: CustomBaseMeshAlgo + ?Sized>(
    algo: &mut A,
    model: &mut MeshModel,
    parameters: &MeshParameters,
) -> Result<(), String> {
    let mut structure = build_data_structure(model, 0)?;
    algo.build_base_triangulation(&mut structure)?;
    let nodes_nb = structure.nb_nodes();
    if nodes_nb < 3 {
        return Err("CustomBaseMeshAlgo::perform: fewer than 3 nodes".to_string());
    }
    let mut vertex_order: Vec<i32> = (1..=nodes_nb as i32).collect();
    let mut mesher = Delaun::new_with_data(structure, &mut vertex_order);
    algo.post_process_mesh(&mut mesher, parameters);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Test algorithm counting its invocations (for trait-dispatch checks).
    #[derive(Default)]
    struct CountingAlgo {
        calls: usize,
    }

    impl BaseMeshAlgo for CountingAlgo {
        fn perform(
            &mut self,
            _model: &mut MeshModel,
            _parameters: &MeshParameters,
        ) -> Result<(), String> {
            self.calls += 1;
            Ok(())
        }
    }

    /// A unit square boundary: 4 `Frontier` nodes + 4 `Frontier` boundary links,
    /// optionally split by a `Fixed` internal diagonal.
    fn square_structure(with_diagonal: bool) -> DelaunDataStructure {
        let mut structure = DelaunDataStructure::new(16);
        let points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
        let mut ids = [0i32; 4];
        for (i, &(u, v)) in points.iter().enumerate() {
            ids[i] = structure.add_node(DelaunVertex::new_parametric(u, v, VertexState::Frontier));
        }
        for i in 0..4 {
            structure.add_link(ids[i], ids[(i + 1) % 4], VertexState::Frontier);
        }
        if with_diagonal {
            structure.add_link(ids[0], ids[2], VertexState::Fixed);
        }
        structure
    }

    #[test]
    fn base_mesh_algo_trait_dispatch() {
        let mut model = MeshModel::default();
        let params = MeshParameters::default();

        // A plain implementation is invocable through the trait.
        let mut counting = CountingAlgo::default();
        BaseMeshAlgo::perform(&mut counting, &mut model, &params).expect("counting algo performs");
        assert_eq!(counting.calls, 1);

        // The concrete Delaunay algorithm is dispatchable as a trait object.
        let mut delaunay: Box<dyn BaseMeshAlgo> = Box::new(DelaunayBaseMeshAlgo::new());
        assert!(
            delaunay.perform(&mut model, &params).is_err(),
            "an empty model has no face to triangulate"
        );
    }

    #[test]
    fn delaunay_base_triangulates_point_set() {
        let structure = square_structure(false);
        let mut algo = DelaunayBaseMeshAlgo::new();
        algo.triangulate(structure, &MeshParameters::default()).expect("triangulate square");

        let mesher = algo.mesher().expect("mesher stored");
        let ds = mesher.result();
        // 4 hull vertices => 2N - 2 - h = 8 - 2 - 4 = 2 triangles.
        let triangles = ds.elements_of_domain().len();
        assert_eq!(triangles, 2, "square triangulates to 2 triangles, got {triangles}");
        // Every triangle references three live, distinct nodes.
        for &id in ds.elements_of_domain() {
            let verts = ds.get_element(id).vertex_indices;
            assert!(verts[0] != verts[1] && verts[1] != verts[2] && verts[0] != verts[2]);
            for &w in &verts {
                assert_ne!(ds.get_node(w).state, VertexState::Deleted);
            }
        }
        // The four boundary links survive the triangulation as frontier edges.
        assert_eq!(mesher.frontier().len(), 4, "4 frontier edges expected");
    }

    #[test]
    fn constraint_edges_registered_and_classified() {
        let structure = square_structure(true);
        let mut algo = DelaunayBaseMeshAlgo::new();
        algo.triangulate(structure, &MeshParameters::default()).expect("triangulate");

        let mesher = algo.mesher().expect("mesher stored");
        // The fixed diagonal is preserved and classified as an internal edge.
        assert_eq!(mesher.internal_edges().len(), 1, "one fixed diagonal expected");
        // The four boundary links stay classified as frontier.
        assert_eq!(mesher.frontier().len(), 4, "four frontier edges expected");
        // The enforced diagonal splits the square into two triangles.
        assert_eq!(mesher.result().elements_of_domain().len(), 2);
    }

    #[test]
    fn delaunay_base_perform_errors_without_face() {
        let mut model = MeshModel::default();
        let params = MeshParameters::default();
        let mut algo = DelaunayBaseMeshAlgo::new();
        assert!(algo.perform(&mut model, &params).is_err());
        assert!(algo.mesher().is_none());
    }

    #[test]
    fn custom_delaunay_base_mesh_algo_post_process() {
        // A minimal custom algorithm reporting a 2x2 acceleration grid.
        struct InnerCustom;

        impl BaseMeshAlgo for InnerCustom {
            fn perform(
                &mut self,
                _model: &mut MeshModel,
                _parameters: &MeshParameters,
            ) -> Result<(), String> {
                Ok(())
            }
        }

        impl ConstrainedBaseMeshAlgo for InnerCustom {
            fn get_cells_count(&self, _vertices_nb: usize) -> (i32, i32) {
                (2, 2)
            }
        }

        impl CustomBaseMeshAlgo for InnerCustom {
            fn build_base_triangulation(
                &mut self,
                _structure: &mut DelaunDataStructure,
            ) -> Result<(), String> {
                Ok(())
            }
        }

        // A square triangulates into 2 triangles.
        let points = [
            DelaunVertex::new_parametric(0.0, 0.0, VertexState::Free),
            DelaunVertex::new_parametric(1.0, 0.0, VertexState::Free),
            DelaunVertex::new_parametric(1.0, 1.0, VertexState::Free),
            DelaunVertex::new_parametric(0.0, 1.0, VertexState::Free),
        ];
        let mut mesher = Delaun::new_vertices(&points);
        // OCCT BRepMesh_Delaun.cxx:703 calls ProcessConstraints() unconditionally;
        // frontierAdjust() ends with cleanupMesh() (cxx:1028) which prunes boundary
        // triangles whose neighbour touches the super-triangle. This square carries
        // only Free links, so cleanupMesh leaves no triangle at all.
        assert_eq!(mesher.result().elements_of_domain().len(), 0);

        // The wrapper re-initializes the circle tool (2x2 cells) and delegates
        // to the inner post-process, leaving the mesh untouched.
        let mut wrapper = CustomDelaunayBaseMeshAlgo::new(InnerCustom);
        let params = MeshParameters::default();
        ConstrainedBaseMeshAlgo::post_process_mesh(&mut wrapper, &mut mesher, &params);
        // Same cleanupMesh() branch as above: the post-processed mesh stays empty.
        assert_eq!(
            mesher.result().elements_of_domain().len(),
            0,
            "post-process must not change the triangulation"
        );
    }
}
