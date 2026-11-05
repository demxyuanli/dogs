//! `ChFi2d` — fillets and chamfers on a planar wire.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi2d/` (OCCT 8.0.0):
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `ChFi2d.cxx` / `ChFi2d.hxx` | 92 / 64 | [`free_functions`] |
//! | `ChFi2d_ConstructionError.hxx` | 41 | [`ConstructionError`] |
//! | `ChFi2d_ChamferAPI.cxx` | 160 | [`chamfer_api`] |
//! | `ChFi2d_AnaFilletAlgo.cxx` | 1101 | [`ana_fillet_algo`] |
//! | `ChFi2d_FilletAlgo.cxx` | 977 | [`fillet_algo`] |
//! | `ChFi2d_FilletAPI.cxx` | 136 | [`fillet_api`] |
//! | `ChFi2d_Builder.cxx` / `ChFi2d_Builder_0.cxx` | 1289 / 890 | [`builder`] |
//!
//! Divergences that are not heuristics:
//! - `ChFi2d_ConstructionError` is a C++ unscoped enum with no explicit values,
//!   so the port keeps the declaration order in [`ConstructionError`],
//! - OCCT raises `Standard_TypeMismatch` on a `TopoDS::Edge` cast of a non-edge
//!   and dereferences null handles without checking; those arms are the `Err`
//!   of the ported functions.

pub mod ana_fillet_algo;
pub mod builder;
pub mod chamfer_api;
pub mod fillet_algo;
pub mod fillet_api;
pub mod free_functions;

pub use ana_fillet_algo::{AnaFilletAlgo, AnaFilletResult};
pub use builder::ChFi2dBuilder;
pub use chamfer_api::ChamferApi;
pub use fillet_algo::{FilletAlgo, FilletAlgoResult};
pub use fillet_api::FilletApi;
pub use free_functions::{common_vertex, find_connected_edges};

/// `ChFi2d_ConstructionError` (`ChFi2d_ConstructionError.hxx:22-39`).
///
/// The declaration order is load-bearing: OCCT's enum has no explicit values
/// and `ChFi2d_Builder` compares the answer against `ChFi2d_Ready` /
/// `ChFi2d_IsDone`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConstructionError {
    /// The face is not planar.
    NotPlanar,
    /// The face is null.
    NoFace,
    /// The two faces used for the initialisation are uncompatible.
    InitialisationError,
    /// The parameters as distances or angle for chamfer are less or equal to
    /// zero.
    ParametersError,
    /// The initialization has been successful.
    Ready,
    IsDone,
    /// The algorithm could not find a solution.
    ComputationError,
    /// The vertex given to locate the fillet or the chamfer is not connected
    /// to 2 edges.
    ConnexionError,
    /// The two edges connected to the vertex are tangent.
    TangencyError,
    /// The first edge is degenerated.
    FirstEdgeDegenerated,
    /// The last edge is degenerated.
    LastEdgeDegenerated,
    /// The two edges are degenerated.
    BothEdgesDegenerated,
    /// One or the two edges connected to the vertex is a fillet or a chamfer,
    /// or is not a line or a circle.
    NotAuthorized,
}
