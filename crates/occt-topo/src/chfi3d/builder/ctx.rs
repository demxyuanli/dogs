//! The slice of `ChFi3d_Builder` state that the subclass virtuals read.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder.hxx:77-848`
//! (OCCT 8.0.0).
//!
//! In C++ the `ChFi3d_ChBuilder` / `ChFi3d_FilBuilder` overrides are member
//! functions of the *same* object, so `SimulSurf` / `PerformFirstSection` /
//! `PerformSurf` read the base members directly. The port keeps the base as a
//! field of the subclass (see
//! [`ChFi3dFilBuilder`](super::fil_builder::ChFi3dFilBuilder)), and the base
//! drivers hold it mutably while they dispatch to
//! [`ChFi3dSurfOps`](super::surf_ops::ChFi3dSurfOps). Passing the whole base
//! would alias that borrow, so the drivers hand the subclasses this
//! immutable view of exactly the members the four bodies read:
//!
//! | OCCT member | used by |
//! |---|---|
//! | `tolapp3d` (`hxx:744`) | `SimulSurf` (`ChFi3d_FilBuilder.cxx:726-741`), `PerformFirstSection` (`cxx:1517`) |
//! | `myVEMap` (`hxx:752`) | `SearchFace` (`ChFi3d_Builder_2.cxx:1569`) |
//! | `myEFMap` (`hxx:748`) | `SearchFace` (`_2.cxx:1570-1619`, `:1624`) |
//! | `tolesp` (`hxx:742`) | the `SimulSurf` overloads take `TolGuide` from their own arguments; kept for the `PerformSurf` bodies |
//! | `tol2d` (`hxx:743`) | `PerformSurf` overloads (`TolGuide` of the edge / face arms) |
//! | `tolapp2d` (`hxx:745`) | `SimulData(SurfRst)` (`ChFi3d_Builder_6.cxx:1391-1445`), reached by the two edge / face `SimulSurf` overloads |
//!
//! `SearchFace` itself (`ChFi3d_Builder_2.cxx:1523-1627`) reads nothing else,
//! so it lives on this view and [`ChFi3dBuilder::search_face`] forwards to it.

use crate::chfi_ds::ChFiDSMap;

use super::ChFi3dBuilder;

/// The immutable base-class slice the subclass virtuals receive instead of
/// `this`.
#[derive(Clone, Copy)]
pub struct ChFi3dCtx<'a> {
    /// `tolapp3d` (`hxx:744`): the 3D approximation tolerance.
    pub tolapp3d: f64,
    /// `tolapp2d` (`hxx:745`): the 2D approximation tolerance.
    pub tolapp2d: f64,
    /// `tol2d` (`hxx:743`).
    pub tol2d: f64,
    /// `tolesp` (`hxx:742`): the guideline tolerance chosen by `SetParams`.
    pub tolesp: f64,
    /// `myVEMap` (`hxx:752`).
    pub ve_map: &'a ChFiDSMap,
    /// `myEFMap` (`hxx:748`).
    pub ef_map: &'a ChFiDSMap,
}

impl ChFi3dBuilder {
    /// The [`ChFi3dCtx`] view of this builder, handed to
    /// [`ChFi3dSurfOps`](super::surf_ops::ChFi3dSurfOps) for one call.
    pub fn ctx(&self) -> ChFi3dCtx<'_> {
        ChFi3dCtx {
            tolapp3d: self.tolapp3d,
            tolapp2d: self.tolapp2d,
            tol2d: self.tol2d,
            tolesp: self.tolesp,
            ve_map: &self.my_ve_map,
            ef_map: &self.my_ef_map,
        }
    }
}
