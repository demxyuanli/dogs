//! `BRepExtrema_ExtFF` — the extrema between two faces.
//!
//! Source: `ModelingAlgorithms/TKTopAlgo/BRepExtrema/BRepExtrema_ExtFF.cxx:30-117`
//! and `.hxx:25-77` (OCCT 8.0.0).
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `(F1, F2)` | `cxx:32-35` | [`ExtFF::new`] |
//! | `Initialize(F2)` | `cxx:40-55` | [`ExtFF::initialize`] |
//! | `Perform(F1, F2)` | `cxx:59-116` | [`ExtFF::perform`] |
//! | `IsDone` / `IsParallel` / `NbExt` / `SquareDistance` | `hxx:44-61` | the same names |
//! | `ParameterOnFace1/2` / `PointOnFace1/2` | `cxx:70-105` | the same names |
//!
//! # Translation notes
//!
//! * `myExtSS` is `Extrema_ExtSS`, already ported as
//!   [`occt_geom::extrema_ss::ExtremaExtSS`] with the matching
//!   `initialize` / `perform` / `is_done` / `is_parallel` / `nb_ext` /
//!   `square_distance` / `points` API.
//! * `Surf.GetType() == GeomAbs_OtherSurface` (`cxx:42`, `cxx:66`) is the
//!   non-geometric guard; the port tests for a registered face surface.
//! * `Tol = min(Tol(F), Confusion)` then `min(UResolution, VResolution)`
//!   raised to `PConfusion` (`cxx:48-51` / `cxx:71-74`) uses
//!   [`occt_geom::approx_same_parameter::u_resolution`] /
//!   [`occt_geom::approx_same_parameter::v_resolution`] (the
//!   `GeomAdaptor_Surface` resolutions).
//! * `BRepClass_FaceClassifier::Perform(F, Puv, Tol)` (`cxx:100` / `cxx:107`)
//!   is [`crate::brep_class::FaceClassifier::perform`]; accepted states `ON` /
//!   `IN` (`cxx:102` / `cxx:109`).
//! * The parallel arm keeps only `SquareDistance(1)` and allocates no points
//!   (`cxx:82-85`), exactly as `Extrema_ExtSS` reports `NbExt() == 1` with an
//!   empty `myPOnS1` / `myPOnS2` for that arm — [`ExtremaExtSS::points`] then
//!   has nothing to return, so [`ExtFF::point_on_face1`] reports the origin.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d};
use occt_core::precision::Precision;
use occt_geom::extrema_ss::ExtremaExtSS;
use occt_geom::extrema_surf::ExtremaPOnSurf;
use occt_geom::Surface;

use crate::brep_class::FaceClassifier;
use crate::fclass2d::FaceState;
use crate::brep_tool::BRepTool;
use crate::shape::Face;

/// `BRepExtrema_ExtFF` (`BRepExtrema_ExtFF.hxx:25-77`).
pub struct ExtFF {
    /// `myHS` (`hxx:75`): F2's surface, `None` for the `OtherSurface` guard.
    my_hs: Option<Arc<dyn Surface>>,
    /// `myExtSS`'s `(U1, U2, V1, V2, Tol)` initialisation (`cxx:52-54`).
    my_window: (f64, f64, f64, f64, f64),
    /// `mySqDist` (`hxx:73`).
    my_sq_dist: Vec<f64>,
    /// `myPointsOnS1` (`hxx:74`).
    my_points_on_s1: Vec<ExtremaPOnSurf>,
    /// `myPointsOnS2` (`hxx:74`).
    my_points_on_s2: Vec<ExtremaPOnSurf>,
    /// `myExtSS.IsDone()` / `IsParallel()`.
    my_done: bool,
    my_parallel: bool,
}

impl ExtFF {
    /// `BRepExtrema_ExtFF(F1, F2)` (`cxx:32-35`).
    pub fn new(f1: &Face, f2: &Face) -> Self {
        let mut this = Self {
            my_hs: None,
            my_window: (0.0, 0.0, 0.0, 0.0, 0.0),
            my_sq_dist: Vec::new(),
            my_points_on_s1: Vec::new(),
            my_points_on_s2: Vec::new(),
            my_done: false,
            my_parallel: false,
        };
        this.initialize(f2);
        this.perform(f1, f2);
        this
    }

    /// `Initialize(F2)` (`cxx:40-55`).
    pub fn initialize(&mut self, f2: &Face) {
        // `cxx:41-45`.
        let surf = match BRepTool::face_surface(f2) {
            Some(s) => s,
            None => return,
        };
        // `cxx:48-51`.
        let tol = BRepTool::face_tolerance(f2).min(Precision::CONFUSION);
        let tol = occt_geom::approx_same_parameter::u_resolution(surf.as_ref(), tol)
            .min(occt_geom::approx_same_parameter::v_resolution(surf.as_ref(), tol))
            .max(Precision::PCONFUSION);
        // `cxx:52-53`.
        let (u1, u2, v1, v2) = BRepTool::uv_bounds(f2);
        self.my_hs = Some(surf);
        self.my_window = (u1, u2, v1, v2, tol);
    }

    /// `Perform(F1, F2)` (`cxx:59-116`).
    pub fn perform(&mut self, f1: &Face, f2: &Face) {
        // `cxx:60-62`.
        self.my_sq_dist.clear();
        self.my_points_on_s1.clear();
        self.my_points_on_s2.clear();
        self.my_done = false;
        self.my_parallel = false;

        // `cxx:64-68`.
        let surf1 = match BRepTool::face_surface(f1) {
            Some(s) => s,
            None => return,
        };
        let surf2 = match &self.my_hs {
            Some(s) => s.clone(),
            None => return,
        };

        // `cxx:71-74`.
        let tol1 = BRepTool::face_tolerance(f1).min(Precision::CONFUSION);
        let tol1 = occt_geom::approx_same_parameter::u_resolution(surf1.as_ref(), tol1)
            .min(occt_geom::approx_same_parameter::v_resolution(surf1.as_ref(), tol1))
            .max(Precision::PCONFUSION);
        // `cxx:75-77`.
        let (u1, u2, v1, v2) = BRepTool::uv_bounds(f1);

        // `myExtSS.Perform(*HS1, U1, U2, V1, V2, Tol1)` (`cxx:78`).
        let mut ex = ExtremaExtSS::default();
        // `Initialize(F2)` already stored S2's window by `cxx:52-54`; replay it
        // here with the same values.
        ex.initialize(
            surf2.as_ref(),
            self.my_window.0,
            self.my_window.1,
            self.my_window.2,
            self.my_window.3,
            self.my_window.4,
        );
        ex.perform(
            surf1.as_ref(),
            surf2.as_ref(),
            u1,
            u2,
            v1,
            v2,
            tol1,
        );
        self.my_done = ex.is_done();
        self.my_parallel = ex.is_parallel();
        if !ex.is_done() {
            return;
        }

        // `cxx:82-85`.
        if ex.is_parallel() {
            self.my_sq_dist.push(ex.square_distance(1));
            return;
        }

        // `cxx:88-115`.
        let mut classifier = FaceClassifier::new();
        let tol2 = BRepTool::face_tolerance(f2);
        for i in 1..=ex.nb_ext() {
            let (p1, p2) = ex.points(i as i32);
            let (pu1, pv1) = p1.parameter();
            classifier.perform(f1, GpPnt2d::new(pu1, pv1), tol1);
            let state1 = classifier.state();
            // `cxx:102-103`.
            if state1 == FaceState::On || state1 == FaceState::In {
                let (pu2, pv2) = p2.parameter();
                classifier.perform(f2, GpPnt2d::new(pu2, pv2), tol2);
                let state2 = classifier.state();
                // `cxx:109-114`.
                if state2 == FaceState::On || state2 == FaceState::In {
                    self.my_sq_dist.push(ex.square_distance(i as i32));
                    self.my_points_on_s1.push(p1);
                    self.my_points_on_s2.push(p2);
                }
            }
        }
    }

    /// `IsDone()` (`hxx:44`).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// `IsParallel()` (`hxx:47`).
    pub fn is_parallel(&self) -> bool {
        self.my_parallel
    }

    /// `NbExt()` (`hxx:50`).
    pub fn nb_ext(&self) -> usize {
        self.my_sq_dist.len()
    }

    /// `SquareDistance(N)` (`hxx:53`).
    pub fn square_distance(&self, n: usize) -> f64 {
        self.my_sq_dist.get(n.wrapping_sub(1)).copied().unwrap_or(0.0)
    }

    /// `ParameterOnFace1(N)` / `ParameterOnFace2(N)`: the `(u, v)` of the `N`th
    /// point on face 1 / face 2.
    pub fn parameter_on_face1(&self, n: usize) -> (f64, f64) {
        self.my_points_on_s1
            .get(n.wrapping_sub(1))
            .map(|p| p.parameter())
            .unwrap_or((0.0, 0.0))
    }

    /// `ParameterOnFace2(N)`.
    pub fn parameter_on_face2(&self, n: usize) -> (f64, f64) {
        self.my_points_on_s2
            .get(n.wrapping_sub(1))
            .map(|p| p.parameter())
            .unwrap_or((0.0, 0.0))
    }

    /// `PointOnFace1(N)`.
    pub fn point_on_face1(&self, n: usize) -> GpPnt {
        self.my_points_on_s1
            .get(n.wrapping_sub(1))
            .map(|p| p.value())
            .unwrap_or_else(GpPnt::zero)
    }

    /// `PointOnFace2(N)`.
    pub fn point_on_face2(&self, n: usize) -> GpPnt {
        self.my_points_on_s2
            .get(n.wrapping_sub(1))
            .map(|p| p.value())
            .unwrap_or_else(GpPnt::zero)
    }
}
