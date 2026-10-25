//! The `TheMultiLineTool` template parameter of `Approx_ComputeLine` /
//! `Approx_BSplComputeLine`.
//!
//! OCCT writes these two algorithms as `.gxx` templates parameterised by
//! `TheMultiLine` (the point container) and `TheMultiLineTool` (a class of
//! static functions over it); each package picks a pair with `#define`s:
//!
//! | instantiation | line | tool |
//! |---|---|---|
//! | `AppDef_Compute` / `AppDef_BSplineCompute` | `AppDef_MultiLine` | `AppDef_MyLineTool` |
//! | `BRepApprox_TheComputeLineOfApprox` / `…TheComputeLineBezierOfApprox` | `BRepApprox_TheMultiLineOfApprox` | `BRepApprox_TheMultiLineToolOfApprox` |
//! | `GeomInt_TheComputeLineOfApprox` / `…` | `GeomInt_TheMultiLineOfApprox` | `GeomInt_TheMultiLineToolOfApprox` |
//!
//! C++ needs the tool as a separate class because it has no traits; every one
//! of its members is a pure function of its `TheMultiLine` argument
//! (`ApproxInt_MultiLineTool.lxx:28-186`), so the port attaches them to the
//! line type itself as this trait. The method set and the per-method
//! `TheMultiLine` argument are unchanged; the algorithms become generic over
//! `L: LineTool` instead of over two template arguments.
//!
//! Translation notes:
//! - The `NCollection_Array1` out-parameters become returned `Vec`s whose
//!   lengths are the `NbP3d` / `NbP2d` of the line, matching the OCCT loops
//!   that fill exactly those slots.
//! - `MakeMLOneMorePoint` returns `bool` plus a `TheMultiLine&` out-parameter
//!   in OCCT; `Option<Self>` carries the same pair (`None` = `false`, and the
//!   out-parameter is then undefined in OCCT, as documented at
//!   `AppDef_MyLineTool.cxx:164-171`).
//! - `Curvature` is not called by either algorithm (`Approx_ComputeLine.gxx`
//!   and `Approx_BSplComputeLine.gxx` only use `Value` / `Tangency` /
//!   `WhatStatus` / `MakeML*`); it is here because `ApproxInt_MultiLineTool`
//!   declares it (`lxx:112-145`).

use occt_core::approx::ApproxStatus;
use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};

/// `TheMultiLineTool` (`ApproxInt_MultiLineTool.lxx:28-186`), implemented by
/// the `TheMultiLine` it reads.
pub trait LineTool: Clone {
    /// `NbP2d` (`lxx:28-32`).
    fn nb_p2d(&self) -> i32;

    /// `NbP3d` (`lxx:35-39`).
    fn nb_p3d(&self) -> i32;

    /// `FirstPoint` (`lxx:41-46`).
    fn first_point(&self) -> i32;

    /// `LastPoint` (`lxx:49-54`).
    fn last_point(&self) -> i32;

    /// `Value(ML, Index, TabPnt)` (`lxx:56-62`).
    fn value_3d(&self, index: i32) -> Vec<GpPnt>;

    /// `Value(ML, Index, TabPnt2d)` (`lxx:65-71`).
    fn value_2d(&self, index: i32) -> Vec<GpPnt2d>;

    /// `Value(ML, Index, TabPnt, TabPnt2d)` (`lxx:74-81`).
    fn value(&self, index: i32) -> (Vec<GpPnt>, Vec<GpPnt2d>);

    /// `Tangency(ML, Index, TabVec)` (`lxx:84-91`).
    fn tangency_3d(&self, index: i32) -> Option<Vec<GpVec>>;

    /// `Tangency(ML, Index, TabVec2d)` (`lxx:94-100`).
    fn tangency_2d(&self, index: i32) -> Option<Vec<GpVec2d>>;

    /// `Tangency(ML, Index, TabVec, TabVec2d)` (`lxx:103-110`).
    fn tangency(&self, index: i32) -> Option<(Vec<GpVec>, Vec<GpVec2d>)>;

    /// `Curvature(ML, Index, TabVec)` (`lxx:112-120`): `false` in the
    /// reference implementation.
    fn curvature_3d(&self, index: i32) -> Option<Vec<GpVec>>;

    /// `Curvature(ML, Index, TabVec2d)` (`lxx:122-130`).
    fn curvature_2d(&self, index: i32) -> Option<Vec<GpVec2d>>;

    /// `Curvature(ML, Index, TabVec, TabVec2d)` (`lxx:132-143`).
    fn curvature(&self, index: i32) -> Option<(Vec<GpVec>, Vec<GpVec2d>)>;

    /// `WhatStatus(ML, I1, I2)` (`lxx:147-158`).
    fn what_status(&self, i1: i32, i2: i32) -> ApproxStatus;

    /// `MakeMLBetween(ML, I1, I2, NbPMin)` (`lxx:162-167`).
    fn make_ml_between(&self, i1: i32, i2: i32, nb_p_min: i32) -> Self
    where
        Self: Sized;

    /// `MakeMLOneMorePoint(ML, I1, I2, indbad, OtherLine)` (`lxx:171-180`).
    fn make_ml_one_more_point(&self, i1: i32, i2: i32, indbad: i32) -> Option<Self>
    where
        Self: Sized;
}

