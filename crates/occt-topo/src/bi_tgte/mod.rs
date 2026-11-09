//! `BiTgte` package - rolling-ball (filler) blend of a shape.
//!
//! Source: `ModelingAlgorithms/TKOffset/BiTgte/` (OCCT 8.0.0).
//!
//! | OCCT class | lines | port |
//! |---|---|---|
//! | `BiTgte_CurveOnEdge` | `cxx:40-302`, `hxx:45-141` | `curve_on_edge.rs` |
//! | `BiTgte_CurveOnVertex` | `cxx:40-253`, `hxx:45-135` | `curve_on_vertex.rs` |
//! | `BiTgte_Blend` | `cxx:1-2664`, `hxx:1-216` | `blend/` |
//! | `BiTgte_ContactType` | `hxx:20-28` | [`ContactType`] in `mod.rs` |
//!
//! `BiTgte_CurveOnEdge` / `BiTgte_CurveOnVertex` are the two `Adaptor3d_Curve`
//! views `BiTgte_Blend` rolls on: the projection of the `EonF` curve onto the
//! rolled edge, and the constant vertex curve.
//!
//! `BiTgte_Blend` is split by OCCT method over `blend/`: the container and the
//! accessors in `blend/mod.rs`, one file per stage (`compute_centers.rs`,
//! `compute_surfaces.rs`, `compute_shape.rs`, `intersect.rs`), the file-local
//! statics in `blend/statics.rs`, `MakeCurve` in `blend/make_curve.rs`, and
//! `myMapSF` in `blend/offset_map.rs`.

mod blend;
mod curve_on_edge;
mod curve_on_vertex;

pub use blend::Blend;
pub use curve_on_edge::CurveOnEdge;
pub use curve_on_vertex::CurveOnVertex;

/// `BiTgte_ContactType` (`BiTgte_ContactType.hxx:20-28`). The declaration
/// order is the `Standard_Integer` value of each enumerator.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum ContactType {
    FaceFace = 0,
    FaceEdge = 1,
    FaceVertex = 2,
    EdgeEdge = 3,
    EdgeVertex = 4,
    VertexVertex = 5,
}
