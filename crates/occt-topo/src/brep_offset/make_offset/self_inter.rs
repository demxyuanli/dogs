//! `BRepOffset_MakeOffset::SelfInter` (`BRepOffset_MakeOffset.cxx:2154-2167`).
//!
//! OCCT's body is a single `throw Standard_NotImplemented()` (`:2166`); the
//! `OCCT_DEBUG` clock block above it (`:2157-2164`) writes to `std::cout` only
//! and is not part of the algorithm. The `Modif` argument is not even named in
//! the reference signature, so nothing is consumed before the throw.
//!
//! Port note: the crate reports `Standard_NotImplemented` with
//! [`panic!`] (the same convention as `blend/chamfer/chamfer.rs:263`), so
//! callers get the identical "this entry point does not exist" behaviour
//! instead of a silent `Ok`.

use super::BRepOffsetMakeOffset;

impl BRepOffsetMakeOffset {
    /// `SelfInter(Modif)` (`cxx:2154-2167`): always `Standard_NotImplemented`.
    pub fn self_inter(&self) -> ! {
        panic!("Standard_NotImplemented: BRepOffset_MakeOffset::SelfInter");
    }
}
