//! `BRepExtrema_ExtPF` — the extrema between a vertex and a face.
//!
//! Source: `ModelingAlgorithms/TKTopAlgo/BRepExtrema/BRepExtrema_ExtPF.cxx:30-100`
//! and `.hxx:27-73` (OCCT 8.0.0).
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `(V, F, Flag, Algo)` | `cxx:32-37` | [`ExtPF::new`] |
//! | `Initialize(F, Flag, Algo)` | `cxx:42-66` | [`ExtPF::initialize`] |
//! | `Perform(V, F)` | `cxx:70-100` | [`ExtPF::perform`] |
//!
//! # Translation notes
//!
//! * `BRepAdaptor_Surface(TheFace, false)` (`cxx:48`) + `GetType() ==
//!   GeomAbs_OtherSurface` (`cxx:50`) is the "non-geometric face" guard; the
//!   port tests for a registered surface, which is the same condition (a face
//!   whose surface the registry does not have is the `OtherSurface` case of a
//!   meshed-only face).
//! * `myExtPS.SetFlag(Flag)` / `SetAlgo(Algo)` (`cxx:61-62`) are recorded and
//!   forwarded to [`ExtPs`]. `Algo` only selects the gradient-based start
//!   points inside `Extrema_GenExtPS`, which the ported engine always uses
//!   (`Extrema_ExtAlgo_Grad` is the only algorithm OCCT's grid refinement
//!   implements in this build).
//! * `BRepClass_FaceClassifier::Perform(TheFace, Puv, Tol)` (`cxx:91`) is
//!   [`crate::brep_class::FaceClassifier::perform`]`, and the accepted states
//!   are `ON` / `IN` (`cxx:93`) — `FaceState::On` / `FaceState::In`.
//! * `mySqDist` / `myPoints` only receive the solutions the classifier accepts,
//!   exactly as in `cxx:93-97`.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d};
use occt_core::precision::Precision;
use occt_geom::extrema_surf::{ExtPs, ExtremaExtAlgo, ExtremaExtFlag};
use occt_geom::Surface;

use crate::brep_class::FaceClassifier;
use crate::fclass2d::FaceState;
use crate::brep_tool::BRepTool;
use crate::shape::{Face, Vertex};

/// `BRepExtrema_ExtPF` (`BRepExtrema_ExtPF.hxx:27-73`).
pub struct ExtPF<'s> {
    /// `mySurf` (`hxx:68`): the face's surface, `None` for the
    /// `GeomAbs_OtherSurface` guard of `cxx:50`.
    my_surf: Option<Arc<dyn Surface>>,
    /// `myExtPS`'s `(U1, U2, V1, V2, TolU, TolV)` initialisation
    /// (`cxx:59-65`).
    my_window: (f64, f64, f64, f64, f64, f64),
    /// `mySqDist` (`hxx:69`).
    my_sq_dist: Vec<f64>,
    /// `myPoints` (`hxx:70`): `Extrema_POnSurf` values.
    my_points: Vec<(f64, f64, GpPnt)>,
    /// `myFlag` / `myAlgo` (`cxx:61-62`).
    my_flag: ExtremaExtFlag,
    my_algo: ExtremaExtAlgo,
    /// `myExtPS.IsDone()` (`hxx:49`).
    my_done: bool,
    /// Phantom so the borrow of the surface lifetime is explicit at the call
    /// site, as in OCCT where `mySurf` refers to the face.
    _marker: std::marker::PhantomData<&'s ()>,
}

impl<'s> ExtPF<'s> {
    /// `BRepExtrema_ExtPF(TheVertex, TheFace, TheFlag, TheAlgo)` (`cxx:32-37`).
    pub fn new(
        the_vertex: &Vertex,
        the_face: &Face,
        the_flag: ExtremaExtFlag,
        the_algo: ExtremaExtAlgo,
    ) -> Self {
        let mut this = Self {
            my_surf: None,
            my_window: (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            my_sq_dist: Vec::new(),
            my_points: Vec::new(),
            my_flag: the_flag,
            my_algo: the_algo,
            my_done: false,
            _marker: std::marker::PhantomData,
        };
        this.initialize(the_face, the_flag, the_algo);
        this.perform(the_vertex, the_face);
        this
    }

    /// `Initialize(TheFace, TheFlag, TheAlgo)` (`cxx:42-66`).
    pub fn initialize(
        &mut self,
        the_face: &Face,
        the_flag: ExtremaExtFlag,
        the_algo: ExtremaExtAlgo,
    ) {
        // `cxx:48`.
        let surf = match BRepTool::face_surface(the_face) {
            Some(s) => s,
            // `cxx:50-53`.
            None => return,
        };
        // `cxx:55-58`.
        let tol = BRepTool::face_tolerance(the_face).min(Precision::CONFUSION);
        let tol_u = occt_geom::approx_same_parameter::u_resolution(surf.as_ref(), tol)
            .max(Precision::PCONFUSION);
        let tol_v = occt_geom::approx_same_parameter::v_resolution(surf.as_ref(), tol)
            .max(Precision::PCONFUSION);
        // `cxx:59-60`.
        let (u1, u2, v1, v2) = BRepTool::uv_bounds(the_face);
        self.my_surf = Some(surf);
        self.my_window = (u1, u2, v1, v2, tol_u, tol_v);
        self.my_flag = the_flag;
        self.my_algo = the_algo;
    }

    /// `Perform(TheVertex, TheFace)` (`cxx:70-100`).
    pub fn perform(&mut self, the_vertex: &Vertex, the_face: &Face) {
        // `cxx:71-72`.
        self.my_sq_dist.clear();
        self.my_points.clear();

        // `cxx:74`.
        let p = BRepTool::vertex_point(the_vertex);
        let surf = match &self.my_surf {
            // `cxx:75-78`.
            Some(s) => s.clone(),
            None => return,
        };

        // `myExtPS.Perform(P)` (`cxx:80`).
        let (u1, u2, v1, v2, tolu, tolv) = self.my_window;
        let mut ex = ExtPs::new();
        ex.set_flag(self.my_flag);
        ex.set_algo(self.my_algo);
        ex.initialize(surf.as_ref(), u1, u2, v1, v2, tolu, tolv);
        ex.perform(&p);
        self.my_done = ex.is_done();

        // `cxx:82-99`.
        if ex.is_done() {
            let mut classifier = FaceClassifier::new();
            let tol = BRepTool::face_tolerance(the_face);
            for i in 1..=ex.nb_ext() {
                let (pu, pv, _pt) = ex.point(i);
                classifier.perform(the_face, GpPnt2d::new(pu, pv), tol);
                let state = classifier.state();
                // `cxx:93-97`.
                if state == FaceState::On || state == FaceState::In {
                    self.my_sq_dist.push(ex.square_distance(i));
                    self.my_points.push((pu, pv, _pt));
                }
            }
        }
    }

    /// `IsDone()` (`hxx:49`): `myExtPS.IsDone()`.
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// `NbExt()` (`hxx:48`).
    pub fn nb_ext(&self) -> usize {
        self.my_sq_dist.len()
    }

    /// `SquareDistance(N)` (`hxx:51`).
    pub fn square_distance(&self, n: usize) -> f64 {
        self.my_sq_dist.get(n.wrapping_sub(1)).copied().unwrap_or(0.0)
    }

    /// `Point(N)` (`hxx:54`): the `Extrema_POnSurf` at `N`.
    pub fn point(&self, n: usize) -> (f64, f64, GpPnt) {
        self.my_points
            .get(n.wrapping_sub(1))
            .copied()
            .unwrap_or((0.0, 0.0, GpPnt::zero()))
    }

    /// `Parameter(N, U, V)` (`hxx:60-61`).
    pub fn parameter(&self, n: usize) -> (f64, f64) {
        let (u, v, _) = self.point(n);
        (u, v)
    }
}
