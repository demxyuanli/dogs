//! `ChFi3d_Builder`'s corner-side entry points.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder.hxx:388`
//! (`SetRegul`), `:560-576` (the five corner pure virtuals), `:564` /
//! `:780-801` (the five non-virtual corner members, of which `MoreSurfdata` is
//! ported - see below) (OCCT 8.0.0).
//!
//! [`super::ChFi3dBuilder::compute`], `perform_fillet_on_vertex` and
//! `extent_analyse` call these through `this->`. In OCCT the six pure virtuals
//! are implemented by `ChFi3d_FilBuilder` / `ChFi3d_ChBuilder`, which are not
//! ported yet; the six base members are `ChFi3d_Builder_C1.cxx` bodies whose
//! `ChFiKPart` / `BRepBlend` dependencies are not ported either. Both groups
//! are therefore stated as a trait the caller implements (or fails), so the
//! orchestration layer keeps the OCCT call graph without a stub that pretends
//! to work. `MoreSurfdata` is **not** part of the trait: it is ported as
//! [`super::ChFi3dBuilder::more_surfdata`] (`find_face.rs` /
//! `more_surfdata.rs`), together with the `ChFi3d_Builder_C1.cxx` file statics
//! in [`crate::chfi3d::builder_c1`]. The trait carries the OCCT signatures with
//! `occ::handle` / `TopoDS_*` mapped the way the rest of `chfi3d` does:
//!
//! - `const occ::handle<ChFiDS_Stripe>&` → [`StripeHandle`],
//! - `const NCollection_List<occ::handle<ChFiDS_Stripe>>&` (the stripe list of
//!   one vertex, `myVDataMap(Index)`) → `&[StripeHandle]`,
//! - `const TopoDS_Vertex&` → `&TopoShape` (the port keeps vertices as shapes
//!   in `ChFiDS_StripeMap`),
//! - the OCCT `void` bodies that raise are `Result<(), String>`.
//!
//! [`ChFi3dSurfOps`]: super::ChFi3dSurfOps

use crate::chfi_ds::StripeHandle;
use crate::shape::TopoShape;

/// The corner-side entry points of `ChFi3d_Builder`.
pub trait ChFi3dCornerOps {
    /// `SetRegul()` (`hxx:388`): the regularities are coded after the cut.
    fn set_regul(&mut self) -> Result<(), String>;
    /// `PerformTwoCorner(Index)` (`hxx:560`).
    fn perform_two_corner(&mut self, index: usize) -> Result<(), String>;

    /// `PerformThreeCorner(Index)` (`hxx:562`).
    fn perform_three_corner(&mut self, index: usize) -> Result<(), String>;

    /// `ExtentOneCorner(V, S)` (`hxx:566-568`).
    fn extent_one_corner(&mut self, v: &TopoShape, s: StripeHandle) -> Result<(), String>;

    /// `ExtentTwoCorner(V, LS)` (`hxx:569-572`).
    fn extent_two_corner(&mut self, v: &TopoShape, ls: &[StripeHandle]) -> Result<(), String>;

    /// `ExtentThreeCorner(V, LS)` (`hxx:573-576`).
    fn extent_three_corner(&mut self, v: &TopoShape, ls: &[StripeHandle]) -> Result<(), String>;

    /// `PerformOneCorner(Index, PrepareOnSame = false)` (`hxx:780`, body
    /// `ChFi3d_Builder_C1.cxx:611-1647`).
    ///
    /// Not a virtual in OCCT: the base-class body is ported as the inherent
    /// [`super::ChFi3dBuilder::perform_one_corner`] (`builder/perform_one_corner`,
    /// which also reports the two `PerformIntersectionAtEnd` early returns).
    /// Callers still reach it through this trait because `PerformOneCorner` is
    /// called from `PerformFilletOnVertex` with `this`, and the subclasses that
    /// implement this trait (`ChFi3d_FilBuilder` / `ChFi3d_ChBuilder`) are not
    /// ported yet; their implementations must forward to the inherent method.
    fn perform_one_corner(&mut self, index: usize, prepare_on_same: bool) -> Result<(), String>;

    /// `PerformMoreSurfdata(Index)` (`hxx:784`, body
    /// `ChFi3d_Builder_C1.cxx:3771-4531`), `UNPORTED`.
    fn perform_more_surfdata(&mut self, index: usize) -> Result<(), String>;

    /// `PerformIntersectionAtEnd(Index)` (`hxx:786`, body
    /// `ChFi3d_Builder_C1.cxx:1828-3770`).
    ///
    /// `cxx:1828-2136` (the prologue: counts, sewing edge, regular adjacent
    /// edges, the `PerformMoreSurfdata` / `PerformMoreThreeCorner` dispatches,
    /// `ChFi3d_BoundSurf`, the end-curve orientation and `checkShrink`) is
    /// ported as the inherent
    /// [`super::ChFi3dBuilder::perform_intersection_at_end`]; `cxx:2136-3770`
    /// is UNPORTED and answers [`super::EndIntersectionExit::RestUnported`].
    /// Implementations must forward to the inherent method and dispatch its
    /// exits the way `PerformOneCorner` dispatches
    /// [`super::OneCornerExit::IntersectionAtEnd`].
    fn perform_intersection_at_end(&mut self, index: usize) -> Result<(), String>;

    /// `PerformMoreThreeCorner(Index, nbcourb)` (`hxx:564`, body
    /// `ChFi3d_Builder_C1.cxx:1237-...`), `UNPORTED`.
    ///
    /// Not a virtual in OCCT either (no `= 0` on `hxx:564`); it sits here with
    /// the same injected-implementation rationale.
    fn perform_more_three_corner(&mut self, index: usize, nbcourb: usize)
        -> Result<(), String>;
}
