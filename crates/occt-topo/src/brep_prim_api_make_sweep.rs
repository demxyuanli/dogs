//! `BRepPrimAPI_MakeSweep` — the abstract root of the swept primitives.
//!
//! Source: `ModelingAlgorithms/TKPrim/BRepPrimAPI/BRepPrimAPI_MakeSweep.hxx:42-52`
//! (OCCT 8.0.0). The class has no `.cxx`: the header only declares two pure
//! virtuals on top of `BRepBuilderAPI_MakeShape` (`hxx:48`, `hxx:51`).
//!
//! The port has no `BRepBuilderAPI_MakeShape` base object — every wrapper keeps
//! its own `my_shape` / `my_done` / `my_generated` triple — so this trait is the
//! only shared part left to translate: the two queries every sweep answers.

use crate::shape::TopoShape;

/// `BRepPrimAPI_MakeSweep` (`BRepPrimAPI_MakeSweep.hxx:42-52`).
pub trait MakeSweep {
    /// `FirstShape()` (`hxx:48`): the shape swept from the profile start.
    fn first_shape(&self) -> Option<TopoShape>;

    /// `LastShape()` (`hxx:51`): the shape swept from the profile end.
    fn last_shape(&self) -> Option<TopoShape>;
}
