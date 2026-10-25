//! `ApproxInt_MultiLineTool` — the `TheMultiLineTool` the `BRepApprox` /
//! `GeomInt` packages instantiate `Approx_ComputeLine` with.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/ApproxInt/ApproxInt_MultiLineTool.lxx:22-186`
//! (OCCT 8.0.0).
//!
//! Every member of the OCCT class is a one-line forward to the matching
//! `ApproxInt_MultiLine` method (`lxx:29-179`), and the algorithms only ever
//! reach it through the `TheMultiLineTool` template parameter. The Rust port
//! attaches those forwards to the line type itself through
//! [`occt_geom::approx::line_tool::LineTool`] (see that trait for why), so
//! this file is the `BRepApprox_TheMultiLineToolOfApprox` /
//! `GeomInt_TheMultiLineToolOfApprox` instantiation rather than a separate
//! `struct`.
//!
//! Translation notes:
//! - `Value(ML, Index, TabPnt)` writes only `TabPnt(1)` (`gxx:178-183`); the
//!   port returns the written slot. `NbP3d` is 1 for the `BRepApprox` and
//!   `GeomInt` lines (`gxx:56-75` sets it from the constructor argument), so
//!   the shorter vector is the same array the callers see.
//! - `Tangency` / `Curvature` return `bool` plus `NCollection_Array1` out
//!   parameters; the port returns `Option` (`false` = `None`) and the
//!   zero-filled array the OCCT bodies leave behind on failure (`gxx:234-239`,
//!   `:279-287`) as `None`'s counterpart, so a failed call cannot be told
//!   apart from a zero tangent by the caller — the algorithms test the flag.
//! - `MakeMLOneMorePoint(ML, I1, I2, indbad, OtherLine)` writes `OtherLine`
//!   and returns `true` (`gxx:631-…`); `Option<Self>` carries the pair.
//! - `Dump` (`lxx:181-184`) is not ported: it is `#ifdef OCCT_DEBUG` only and
//!   `Approx_ComputeLine` reaches it through `DUMP` macros.

use occt_core::approx::ApproxStatus;
use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_geom::approx::line_tool::LineTool;

use crate::approx_int::multi_line::{LineSource, MultiLine};

impl<L: LineSource> LineTool for MultiLine<L> {
    /// `NbP2d` (`lxx:29-32`).
    fn nb_p2d(&self) -> i32 {
        MultiLine::nb_p2d(self)
    }

    /// `NbP3d` (`lxx:36-39`).
    fn nb_p3d(&self) -> i32 {
        MultiLine::nb_p3d(self)
    }

    /// `FirstPoint` (`lxx:43-46`).
    fn first_point(&self) -> i32 {
        MultiLine::first_point(self)
    }

    /// `LastPoint` (`lxx:50-53`).
    fn last_point(&self) -> i32 {
        MultiLine::last_point(self)
    }

    /// `Value(ML, Index, TabPnt)` (`lxx:57-62`).
    fn value_3d(&self, index: i32) -> Vec<GpPnt> {
        vec![MultiLine::value_p3d(self, index)]
    }

    /// `Value(ML, Index, TabPnt2d)` (`lxx:66-71`).
    fn value_2d(&self, index: i32) -> Vec<GpPnt2d> {
        MultiLine::value_p2d(self, index)
    }

    /// `Value(ML, Index, TabPnt, TabPnt2d)` (`lxx:75-81`).
    fn value(&self, index: i32) -> (Vec<GpPnt>, Vec<GpPnt2d>) {
        let (p, p2d) = MultiLine::value(self, index);
        (vec![p], p2d)
    }

    /// `Tangency(ML, Index, TabVec)` (`lxx:85-90`).
    fn tangency_3d(&self, index: i32) -> Option<Vec<GpVec>> {
        let (ret, tg) = MultiLine::tangency(self, index);
        if ret {
            Some(vec![tg])
        } else {
            None
        }
    }

    /// `Tangency(ML, Index, TabVec2d)` (`lxx:94-99`).
    fn tangency_2d(&self, index: i32) -> Option<Vec<GpVec2d>> {
        let (ret, tg) = MultiLine::tangency_2d(self, index);
        if ret {
            Some(tg)
        } else {
            None
        }
    }

    /// `Tangency(ML, Index, TabVec, TabVec2d)` (`lxx:103-109`).
    fn tangency(&self, index: i32) -> Option<(Vec<GpVec>, Vec<GpVec2d>)> {
        let (ret, tg, tg2d) = MultiLine::tangency_both(self, index);
        if ret {
            Some((vec![tg], tg2d))
        } else {
            None
        }
    }

    /// `Curvature(ML, Index, TabVec)` (`lxx:113-120`): `return false;`.
    fn curvature_3d(&self, _index: i32) -> Option<Vec<GpVec>> {
        None
    }

    /// `Curvature(ML, Index, TabVec2d)` (`lxx:124-131`): `return false;`.
    fn curvature_2d(&self, _index: i32) -> Option<Vec<GpVec2d>> {
        None
    }

    /// `Curvature(ML, Index, TabVec, TabVec2d)` (`lxx:135-144`):
    /// `return false;`.
    fn curvature(&self, _index: i32) -> Option<(Vec<GpVec>, Vec<GpVec2d>)> {
        None
    }

    /// `WhatStatus(ML, I1, I2)` (`lxx:148-157`): `I1` / `I2` are unnamed and
    /// unused; the status comes from the line.
    fn what_status(&self, _i1: i32, _i2: i32) -> ApproxStatus {
        MultiLine::what_status(self)
    }

    /// `MakeMLBetween(ML, I1, I2, NbPMin)` (`lxx:161-168`).
    fn make_ml_between(&self, i1: i32, i2: i32, nb_p_min: i32) -> Self {
        MultiLine::make_ml_between(self, i1, i2, nb_p_min)
    }

    /// `MakeMLOneMorePoint(ML, I1, I2, indbad, OtherLine)` (`lxx:172-179`).
    fn make_ml_one_more_point(&self, i1: i32, i2: i32, indbad: i32) -> Option<Self> {
        let (made, other) = MultiLine::make_ml_one_more_point(self, i1, i2, indbad);
        if made {
            Some(other)
        } else {
            None
        }
    }
}
