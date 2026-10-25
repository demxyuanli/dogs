//! `ChFiDS_FaceInterference` — the face/fillet interference of one spine point.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiDS/ChFiDS_FaceInterference.hxx:30-78`,
//! `.lxx:16-105`, `.cxx:22-57` (OCCT 8.0.0).

use std::sync::Arc;

use occt_geom2d::curve::Curve2d;

use crate::abs::Orientation;

/// Source: `ChFiDS_FaceInterference.hxx:31-78`.
///
/// Constructor initial values from `cxx:22-28`: `firstParam` / `lastParam` are
/// `0.0`, `lineindex` is `0`, `LineTransition` is `FORWARD`.
#[derive(Clone)]
pub struct FaceInterference {
    /// `firstParam`
    first_param: f64,
    /// `lastParam`
    last_param: f64,
    /// `pCurveOnFace`
    p_curve_on_face: Option<Arc<dyn Curve2d>>,
    /// `pCurveOnSurf`
    p_curve_on_surf: Option<Arc<dyn Curve2d>>,
    /// `lineindex`
    lineindex: i32,
    /// `LineTransition`
    line_transition: Orientation,
}

impl Default for FaceInterference {
    /// `ChFiDS_FaceInterference()` (`cxx:22-28`).
    fn default() -> Self {
        Self {
            first_param: 0.0,
            last_param: 0.0,
            p_curve_on_face: None,
            p_curve_on_surf: None,
            lineindex: 0,
            line_transition: Orientation::Forward,
        }
    }
}

/// `dyn Curve2d` has no `Debug`, so the two pcurve handles print as set/unset
/// rather than by content. Nothing in the port prints them.
impl std::fmt::Debug for FaceInterference {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("FaceInterference")
            .field("first_param", &self.first_param)
            .field("last_param", &self.last_param)
            .field("has_p_curve_on_face", &self.p_curve_on_face.is_some())
            .field("has_p_curve_on_surf", &self.p_curve_on_surf.is_some())
            .field("lineindex", &self.lineindex)
            .field("line_transition", &self.line_transition)
            .finish()
    }
}

impl FaceInterference {
    /// `ChFiDS_FaceInterference()`.
    pub fn new() -> Self {
        Self::default()
    }

    /// `SetInterference(LineIndex, Trans, PCurv1, PCurv2)` (`lxx:19-27`).
    pub fn set_interference(
        &mut self,
        line_index: i32,
        trans: Orientation,
        p_curv1: Option<Arc<dyn Curve2d>>,
        p_curv2: Option<Arc<dyn Curve2d>>,
    ) {
        self.lineindex = line_index;
        self.line_transition = trans;
        self.p_curve_on_face = p_curv1;
        self.p_curve_on_surf = p_curv2;
    }

    /// `SetTransition(Trans)` (`cxx:42-45`).
    pub fn set_transition(&mut self, trans: Orientation) {
        self.line_transition = trans;
    }

    /// `SetFirstParameter(U1)` (`lxx:39-42`).
    pub fn set_first_parameter(&mut self, u1: f64) {
        self.first_param = u1;
    }

    /// `SetLastParameter(U1)` (`lxx:46-49`).
    pub fn set_last_parameter(&mut self, u1: f64) {
        self.last_param = u1;
    }

    /// `SetParameter(U1, IsFirst)` (`cxx:30-40`).
    pub fn set_parameter(&mut self, u1: f64, is_first: bool) {
        if is_first {
            self.set_first_parameter(u1);
        } else {
            self.set_last_parameter(u1);
        }
    }

    /// `LineIndex()` (`lxx:54-57`).
    pub fn line_index(&self) -> i32 {
        self.lineindex
    }

    /// `SetLineIndex(I)` (`lxx:32-35`).
    pub fn set_line_index(&mut self, i: i32) {
        self.lineindex = i;
    }

    /// `Transition()` (`lxx:62-65`).
    pub fn transition(&self) -> Orientation {
        self.line_transition
    }

    /// `PCurveOnFace()` (`lxx:67-70`).
    pub fn p_curve_on_face(&self) -> Option<&Arc<dyn Curve2d>> {
        self.p_curve_on_face.as_ref()
    }

    /// `PCurveOnSurf()` (`lxx:74-77`).
    pub fn p_curve_on_surf(&self) -> Option<&Arc<dyn Curve2d>> {
        self.p_curve_on_surf.as_ref()
    }

    /// `ChangePCurveOnFace()` (`lxx:82-85`).
    pub fn set_p_curve_on_face(&mut self, c: Option<Arc<dyn Curve2d>>) {
        self.p_curve_on_face = c;
    }

    /// `ChangePCurveOnSurf()` (`lxx:88-91`).
    pub fn set_p_curve_on_surf(&mut self, c: Option<Arc<dyn Curve2d>>) {
        self.p_curve_on_surf = c;
    }

    /// `FirstParameter()` (`lxx:96-99`).
    pub fn first_parameter(&self) -> f64 {
        self.first_param
    }

    /// `LastParameter()` (`lxx:104-107`).
    pub fn last_parameter(&self) -> f64 {
        self.last_param
    }

    /// `Parameter(IsFirst)` (`cxx:47-57`).
    pub fn parameter(&self, is_first: bool) -> f64 {
        if is_first {
            self.first_parameter()
        } else {
            self.last_parameter()
        }
    }
}
