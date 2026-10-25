//! `AppDef_MyLineTool` — the `LineTool` adapter that lets the generic
//! approximation algorithms read an `AppDef_MultiLine`.
//!
//! Source: `ModelingData/TKGeomBase/AppDef/AppDef_MyLineTool.hxx:30-100` and
//! `AppDef_MyLineTool.cxx:25-244` (OCCT 8.0.0).
//!
//! Translation notes:
//! - The OCCT tool has overloaded `Value` / `Tangency` / `Curvature`
//!   static members; the caller picks the overload from `NbP2d` / `NbP3d`
//!   (`Approx_BSplComputeLine.gxx:82-86`, `:119-123`). The port keeps the
//!   overloads as separate functions (`value_3d`, `value_2d`, `value`, ...)
//!   filled from the same 1-based `Index`.
//! - The `tabPt` / `tabV` out-parameters become returned vectors whose
//!   length is the pole count of the addressed column, matching the OCCT
//!   loops that write exactly `nbp3d` / `nbp2d` slots (`cxx:50-55`).
//! - `MakeMLOneMorePoint` (`cxx:164-171`) never writes its `ML` out-parameter
//!   and always returns `false`; the port drops the out-parameter and
//!   returns `false` with the same reasoning.
//! - `MakeMLBetween` (`cxx:157-162`) returns its input unchanged ("stub" in
//!   OCCT); the port clones it.

use occt_core::approx::ApproxStatus;
use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};

use crate::app_def::multi_line::MultiLine;
use crate::approx::line_tool::LineTool;

/// `AppDef_MyLineTool` (`AppDef_MyLineTool.hxx:30-100`).
pub struct MyLineTool;

impl MyLineTool {
    /// `FirstPoint` (`cxx:23-26`): always 1.
    pub fn first_point(_ml: &MultiLine) -> i32 {
        1
    }

    /// `LastPoint` (`cxx:28-31`): `ML.NbMultiPoints()`.
    pub fn last_point(ml: &MultiLine) -> i32 {
        ml.nb_multi_points() as i32
    }

    /// `NbP2d` (`cxx:33-36`): `ML.Value(1).NbPoints2d()`.
    pub fn nb_p2d(ml: &MultiLine) -> i32 {
        ml.value(1).nb_points2d() as i32
    }

    /// `NbP3d` (`cxx:38-41`): `ML.Value(1).NbPoints()`.
    pub fn nb_p3d(ml: &MultiLine) -> i32 {
        ml.value(1).nb_points() as i32
    }

    /// `Value(ML, MPointIndex, tabPt)` (`cxx:43-55`).
    pub fn value_3d(ml: &MultiLine, mpoint_index: i32) -> Vec<GpPnt> {
        let mpc = ml.value(mpoint_index);
        (1..=mpc.nb_points() as i32)
            .map(|i| *mpc.point(i))
            .collect()
    }

    /// `Value(ML, MPointIndex, tabPt2d)` (`cxx:57-68`): `Index` counts the
    /// 3d points first.
    pub fn value_2d(ml: &MultiLine, mpoint_index: i32) -> Vec<GpPnt2d> {
        let mpc = ml.value(mpoint_index);
        let nbp3d = mpc.nb_points() as i32;
        (1..=mpc.nb_points2d() as i32)
            .map(|i| *mpc.point2d(nbp3d + i))
            .collect()
    }

    /// `Value(ML, MPointIndex, tabPt, tabPt2d)` (`cxx:70-86`).
    pub fn value(ml: &MultiLine, mpoint_index: i32) -> (Vec<GpPnt>, Vec<GpPnt2d>) {
        let mpc = ml.value(mpoint_index);
        let nbp3d = mpc.nb_points() as i32;
        let pts = (1..=nbp3d).map(|i| *mpc.point(i)).collect();
        let pts2d = (1..=mpc.nb_points2d() as i32)
            .map(|i| *mpc.point2d(nbp3d + i))
            .collect();
        (pts, pts2d)
    }

    /// `Tangency(ML, MPointIndex, tabV)` (`cxx:88-106`): `None` when the
    /// column carries no tangency.
    pub fn tangency_3d(ml: &MultiLine, mpoint_index: i32) -> Option<Vec<GpVec>> {
        let mpc = ml.value(mpoint_index);
        if mpc.is_tangency_point() {
            Some((1..=mpc.nb_points() as i32).map(|i| mpc.tang(i)).collect())
        } else {
            None
        }
    }

    /// `Tangency(ML, MPointIndex, tabV2d)` (`cxx:108-127`).
    pub fn tangency_2d(ml: &MultiLine, mpoint_index: i32) -> Option<Vec<GpVec2d>> {
        let mpc = ml.value(mpoint_index);
        if mpc.is_tangency_point() {
            let nbp3d = mpc.nb_points() as i32;
            Some(
                (1..=mpc.nb_points2d() as i32)
                    .map(|i| mpc.tang2d(nbp3d + i))
                    .collect(),
            )
        } else {
            None
        }
    }

    /// `Tangency(ML, MPointIndex, tabV, tabV2d)` (`cxx:129-154`).
    pub fn tangency(
        ml: &MultiLine,
        mpoint_index: i32,
    ) -> Option<(Vec<GpVec>, Vec<GpVec2d>)> {
        let mpc = ml.value(mpoint_index);
        if mpc.is_tangency_point() {
            let nbp3d = mpc.nb_points() as i32;
            let v = (1..=nbp3d).map(|i| mpc.tang(i)).collect();
            let v2d = (1..=mpc.nb_points2d() as i32)
                .map(|i| mpc.tang2d(nbp3d + i))
                .collect();
            Some((v, v2d))
        } else {
            None
        }
    }

    /// `MakeMLBetween(ML, I1, I2, theNbPoint)` (`cxx:157-162`): the OCCT body
    /// is `return theML; // stub`.
    pub fn make_ml_between(ml: &MultiLine, _i1: i32, _i2: i32, _the_nb_point: i32) -> MultiLine {
        ml.clone()
    }

    /// `MakeMLOneMorePoint(ML, I1, I2, indbad, ML)` (`cxx:164-171`): always
    /// `false`.
    pub fn make_ml_one_more_point(
        _ml: &MultiLine,
        _i1: i32,
        _i2: i32,
        _indbad: i32,
    ) -> bool {
        false
    }

    /// `WhatStatus(ML, I1, I2)` (`cxx:173-176`): `Approx_NoPointsAdded`.
    pub fn what_status(_ml: &MultiLine, _i1: i32, _i2: i32) -> ApproxStatus {
        ApproxStatus::NoPointsAdded
    }

    /// `Curvature(ML, MPointIndex, tabV)` (`cxx:178-197`).
    pub fn curvature_3d(ml: &MultiLine, mpoint_index: i32) -> Option<Vec<GpVec>> {
        let mpc = ml.value(mpoint_index);
        if mpc.is_curvature_point() {
            Some((1..=mpc.nb_points() as i32).map(|i| mpc.curv(i)).collect())
        } else {
            None
        }
    }

    /// `Curvature(ML, MPointIndex, tabV2d)` (`cxx:199-215`).
    pub fn curvature_2d(ml: &MultiLine, mpoint_index: i32) -> Option<Vec<GpVec2d>> {
        let mpc = ml.value(mpoint_index);
        if mpc.is_curvature_point() {
            let nbp3d = mpc.nb_points() as i32;
            Some(
                (1..=mpc.nb_points2d() as i32)
                    .map(|i| mpc.curv2d(nbp3d + i))
                    .collect(),
            )
        } else {
            None
        }
    }

    /// `Curvature(ML, MPointIndex, tabV, tabV2d)` (`cxx:217-244`).
    pub fn curvature(
        ml: &MultiLine,
        mpoint_index: i32,
    ) -> Option<(Vec<GpVec>, Vec<GpVec2d>)> {
        let mpc = ml.value(mpoint_index);
        if mpc.is_curvature_point() {
            let nbp3d = mpc.nb_points() as i32;
            let v = (1..=nbp3d).map(|i| mpc.curv(i)).collect();
            let v2d = (1..=mpc.nb_points2d() as i32)
                .map(|i| mpc.curv2d(nbp3d + i))
                .collect();
            Some((v, v2d))
        } else {
            None
        }
    }
}

/// `AppDef_MyLineTool` as the `TheMultiLineTool` of `Approx_ComputeLine` /
/// `Approx_BSplComputeLine` (the `AppDef_Compute` / `AppDef_BSplineCompute`
/// instantiation). `Approx_ComputeLine.lxx` has no stand-alone tool of its own;
/// the functions are attached to `AppDef_MultiLine` itself because each of them
/// is a pure function of its line argument.
impl LineTool for MultiLine {
    fn nb_p2d(&self) -> i32 {
        MyLineTool::nb_p2d(self)
    }

    fn nb_p3d(&self) -> i32 {
        MyLineTool::nb_p3d(self)
    }

    fn first_point(&self) -> i32 {
        MyLineTool::first_point(self)
    }

    fn last_point(&self) -> i32 {
        MyLineTool::last_point(self)
    }

    fn value_3d(&self, index: i32) -> Vec<GpPnt> {
        MyLineTool::value_3d(self, index)
    }

    fn value_2d(&self, index: i32) -> Vec<GpPnt2d> {
        MyLineTool::value_2d(self, index)
    }

    fn value(&self, index: i32) -> (Vec<GpPnt>, Vec<GpPnt2d>) {
        MyLineTool::value(self, index)
    }

    fn tangency_3d(&self, index: i32) -> Option<Vec<GpVec>> {
        MyLineTool::tangency_3d(self, index)
    }

    fn tangency_2d(&self, index: i32) -> Option<Vec<GpVec2d>> {
        MyLineTool::tangency_2d(self, index)
    }

    fn tangency(&self, index: i32) -> Option<(Vec<GpVec>, Vec<GpVec2d>)> {
        MyLineTool::tangency(self, index)
    }

    fn curvature_3d(&self, index: i32) -> Option<Vec<GpVec>> {
        MyLineTool::curvature_3d(self, index)
    }

    fn curvature_2d(&self, index: i32) -> Option<Vec<GpVec2d>> {
        MyLineTool::curvature_2d(self, index)
    }

    fn curvature(&self, index: i32) -> Option<(Vec<GpVec>, Vec<GpVec2d>)> {
        MyLineTool::curvature(self, index)
    }

    fn what_status(&self, i1: i32, i2: i32) -> ApproxStatus {
        MyLineTool::what_status(self, i1, i2)
    }

    fn make_ml_between(&self, i1: i32, i2: i32, nb_p_min: i32) -> MultiLine {
        MyLineTool::make_ml_between(self, i1, i2, nb_p_min)
    }

    fn make_ml_one_more_point(&self, i1: i32, i2: i32, indbad: i32) -> Option<MultiLine> {
        // `cxx:164-171` returns `false` and never writes `ML`.
        if MyLineTool::make_ml_one_more_point(self, i1, i2, indbad) {
            Some(self.clone())
        } else {
            None
        }
    }
}
