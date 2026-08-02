//! Port of OCCT factories — Wave 4 BRepMesh.
//!
//! Sources (`src/ModelingAlgorithms/TKMesh/BRepMesh/`):
//! - `BRepMesh_DiscretFactory` — top-level algorithm factory: holds the
//!   [`MeshParameters`] and creates the per-face 2D triangulation algorithm
//!   by [`MeshAlgoType`].
//! - `BRepMesh_DiscretAlgoFactory` — registry of discretization factories
//!   (the built-in `BRepMesh_IncrementalMeshFactory` registers "FastDiscret").
//! - `BRepMesh_MeshAlgoFactory` / `BRepMesh_DelabellaMeshAlgoFactory` — the
//!   per-surface-type algorithm factories, ported as the name↔type mapping and
//!   the [`MeshAlgoType`]-driven [`create_mesh_algo`] dispatcher.
//! - `BRepMesh_IncrementalMeshFactory` — creates the [`IncrementalMesh`]
//!   pipeline entry point.

use crate::shape::TopoShape;

use super::delabella::*;
use super::incremental_mesh::*;
use super::mesh_algo::*;
use super::parameters::{MeshAlgoType, MeshParameters};

/// Creates the 2D triangulation algorithm for the given algorithm type.
///
/// Source: the `GetAlgo` dispatch of `BRepMesh_MeshAlgoFactory` /
/// `BRepMesh_DelabellaMeshAlgoFactory`, collapsed to the algorithm-type level.
///
/// * [`MeshAlgoType::Delaunay`] → the classic Watson [`DelaunayBaseMeshAlgo`].
/// * [`MeshAlgoType::Delabella`] → the Delabella algorithm
///   ([`super::delabella::DelabellaBaseMeshAlgo`]).
/// * [`MeshAlgoType::DFPreOscar`] / [`MeshAlgoType::DFPresta`] (legacy
///   divide-and-fix variants) and [`MeshAlgoType::Unknown`] → the default
///   Watson algorithm.
pub fn create_mesh_algo(algo_type: MeshAlgoType) -> Box<dyn BaseMeshAlgo> {
    match algo_type {
        MeshAlgoType::Delaunay => Box::new(DelaunayBaseMeshAlgo::new()),
        MeshAlgoType::Delabella => Box::new(DelabellaBaseMeshAlgo::new()),
        MeshAlgoType::DFPreOscar | MeshAlgoType::DFPresta | MeshAlgoType::Unknown => {
            Box::new(DelaunayBaseMeshAlgo::new())
        }
    }
}

/// A built-in discretization factory, identified by its registry name.
/// Source: `BRepMesh_DiscretAlgoFactory`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DiscretAlgoFactory {
    name: &'static str,
}

impl DiscretAlgoFactory {
    /// Name the built-in incremental-mesh factory is registered under.
    pub const DEFAULT_NAME: &'static str = "FastDiscret";

    /// A factory with the given registry name.
    pub const fn new(name: &'static str) -> Self {
        Self { name }
    }

    /// The registry name of this factory.
    pub const fn name(&self) -> &'static str {
        self.name
    }

    /// The default (first registered) factory.
    pub const fn default_factory() -> Self {
        Self::new(Self::DEFAULT_NAME)
    }

    /// Finds a built-in factory by name (case-insensitive); `None` when the
    /// name is not registered.
    pub fn find_factory(name: &str) -> Option<Self> {
        if name.eq_ignore_ascii_case(Self::DEFAULT_NAME) {
            Some(Self::default_factory())
        } else {
            None
        }
    }
}

/// Per-algorithm-type factory. Source: `BRepMesh_MeshAlgoFactory` /
/// `BRepMesh_DelabellaMeshAlgoFactory`. Provides the name↔type mapping and the
/// [`MeshAlgoType`]-driven algorithm creation.
pub struct MeshAlgoFactory;

impl MeshAlgoFactory {
    /// Name of the built-in algorithm for the given type.
    pub fn name(algo_type: MeshAlgoType) -> &'static str {
        match algo_type {
            MeshAlgoType::Unknown => "Unknown",
            MeshAlgoType::Delaunay => "Delaunay",
            MeshAlgoType::DFPreOscar => "DFPreOscar",
            MeshAlgoType::DFPresta => "DFPresta",
            MeshAlgoType::Delabella => "Delabella",
        }
    }

    /// The algorithm type for the given name (case-insensitive), accepting the
    /// OCCT aliases ("Watson"/"0" and "1"). Returns [`MeshAlgoType::Unknown`]
    /// for unregistered names.
    pub fn find(name: &str) -> MeshAlgoType {
        match name.to_ascii_lowercase().as_str() {
            "delaunay" | "watson" | "0" => MeshAlgoType::Delaunay,
            "dfpreoscar" => MeshAlgoType::DFPreOscar,
            "dfpresta" => MeshAlgoType::DFPresta,
            "delabella" | "1" => MeshAlgoType::Delabella,
            _ => MeshAlgoType::Unknown,
        }
    }

    /// Creates the algorithm for the given type.
    pub fn create(algo_type: MeshAlgoType) -> Box<dyn BaseMeshAlgo> {
        create_mesh_algo(algo_type)
    }
}

/// Top-level factory for retrieving the meshing algorithm. Source:
/// `BRepMesh_DiscretFactory`.
///
/// Holds the [`MeshParameters`] — their [`MeshAlgoType`] selects the 2D
/// triangulation algorithm — and creates the algorithm on demand.
#[derive(Debug, Clone, PartialEq)]
pub struct DiscretFactory {
    parameters: MeshParameters,
}

impl DiscretFactory {
    /// Creates a factory for the given parameters.
    pub fn new(parameters: MeshParameters) -> Self {
        Self { parameters }
    }

    /// The meshing parameters.
    pub fn parameters(&self) -> &MeshParameters {
        &self.parameters
    }

    /// Mutable access to the meshing parameters.
    pub fn change_parameters(&mut self) -> &mut MeshParameters {
        &mut self.parameters
    }

    /// The 2D triangulation algorithm for the current [`MeshAlgoType`].
    pub fn discret(&self) -> Box<dyn BaseMeshAlgo> {
        create_mesh_algo(self.parameters.mesh_algo)
    }
}

/// Factory for creating [`IncrementalMesh`] instances. Source:
/// `BRepMesh_IncrementalMeshFactory` (registered under "FastDiscret").
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct IncrementalMeshFactory;

impl IncrementalMeshFactory {
    /// Registry name of this factory.
    pub const NAME: &'static str = "FastDiscret";

    /// Creates an incremental mesh for the given shape and parameters.
    ///
    /// The returned algorithm is configured but not yet run — the caller
    /// performs it (OCCT `CreateAlgorithm`, which the caller then `Perform`s).
    pub fn create_incremental_mesh(shape: &TopoShape, params: MeshParameters) -> IncrementalMesh {
        let mut mesh = IncrementalMesh::new();
        mesh.set_shape(shape);
        *mesh.change_parameters() = params;
        mesh
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::data_model::MeshModel;
    use crate::primitives::BRepPrimBox;

    fn unit_box() -> TopoShape {
        BRepPrimBox::make_box(1.0, 1.0, 1.0).solid.0.clone()
    }

    fn params_with(algo: MeshAlgoType) -> MeshParameters {
        MeshParameters { mesh_algo: algo, ..MeshParameters::default() }
    }

    #[test]
    fn discret_factory_creates_algo_by_type() {
        // Delaunay and the legacy divide-and-fix variants all select the
        // Watson DelaunayBaseMeshAlgo, which fails on an empty model (no face
        // with discretized pcurve points to triangulate).
        for algo in [
            MeshAlgoType::Delaunay,
            MeshAlgoType::DFPreOscar,
            MeshAlgoType::DFPresta,
            MeshAlgoType::Unknown,
        ] {
            let mut discret = DiscretFactory::new(params_with(algo)).discret();
            let mut empty = MeshModel::default();
            assert!(
                discret.perform(&mut empty, &MeshParameters::default()).is_err(),
                "{algo:?} must dispatch to DelaunayBaseMeshAlgo"
            );
        }

        // Delabella dispatches to the real DelabellaBaseMeshAlgo from `delabella`
        // (also a real algorithm — fails on an empty model), proving the no-op
        // placeholder of the earlier Wave-4 stub was replaced.
        let mut discret = DiscretFactory::new(params_with(MeshAlgoType::Delabella)).discret();
        let mut empty = MeshModel::default();
        assert!(
            discret.perform(&mut empty, &MeshParameters::default()).is_err(),
            "Delabella must dispatch to the real DelabellaBaseMeshAlgo"
        );
    }

    #[test]
    fn mesh_algo_factory_name_mapping_roundtrips() {
        for (ty, name) in [
            (MeshAlgoType::Delaunay, "Delaunay"),
            (MeshAlgoType::DFPreOscar, "DFPreOscar"),
            (MeshAlgoType::DFPresta, "DFPresta"),
            (MeshAlgoType::Delabella, "Delabella"),
            (MeshAlgoType::Unknown, "Unknown"),
        ] {
            assert_eq!(MeshAlgoFactory::name(ty), name);
            assert_eq!(MeshAlgoFactory::find(name), ty);
        }
        // OCCT aliases and case-insensitivity.
        assert_eq!(MeshAlgoFactory::find("watson"), MeshAlgoType::Delaunay);
        assert_eq!(MeshAlgoFactory::find("0"), MeshAlgoType::Delaunay);
        assert_eq!(MeshAlgoFactory::find("1"), MeshAlgoType::Delabella);
        assert_eq!(MeshAlgoFactory::find("DELABELLA"), MeshAlgoType::Delabella);
        assert_eq!(MeshAlgoFactory::find("nope"), MeshAlgoType::Unknown);
    }

    #[test]
    fn incremental_mesh_factory_creates_configured_mesh() {
        let shape = unit_box();
        let params = MeshParameters {
            mesh_algo: MeshAlgoType::Delaunay,
            deflection: 0.05,
            ..MeshParameters::default()
        };
        let mut mesh = IncrementalMeshFactory::create_incremental_mesh(&shape, params);
        assert!(mesh.shape().is_some(), "shape set");
        assert_eq!(mesh.parameters().mesh_algo, MeshAlgoType::Delaunay);
        assert_eq!(mesh.parameters().deflection, 0.05);
        assert!(!mesh.is_done(), "configured but not performed");
        assert!(mesh.perform().is_ok(), "performs on demand");
        assert!(mesh.is_done(), "IsDone set after perform");
    }

    #[test]
    fn discret_algo_factory_finds_fast_discret() {
        let default = DiscretAlgoFactory::default_factory();
        assert_eq!(default.name(), "FastDiscret");
        assert_eq!(DiscretAlgoFactory::find_factory("FastDiscret"), Some(default));
        assert_eq!(DiscretAlgoFactory::find_factory("fastdiscret"), Some(default));
        assert_eq!(DiscretAlgoFactory::find_factory("Delabella"), None);
        assert_eq!(IncrementalMeshFactory::NAME, "FastDiscret");
    }
}
