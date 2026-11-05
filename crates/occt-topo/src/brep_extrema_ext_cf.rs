//! `BRepExtrema_ExtCF` — the extrema between an edge and a face.
//!
//! Source: `ModelingAlgorithms/TKTopAlgo/BRepExtrema/BRepExtrema_ExtCF.cxx:28-116`
//! and `.hxx:27-76` (OCCT 8.0.0).
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `(E, F)` | `cxx:30-33` | [`ExtCF::new`] |
//! | `Initialize(E, F)` | `cxx:38-58` | [`ExtCF::initialize`] |
//! | `Perform(E, F2)` | `cxx:63-116` | [`ExtCF::perform`] |
//! | `IsDone` / `IsParallel` / `NbExt` / `SquareDistance` | `hxx:45-62` | the same names |
//! | `ParameterOnEdge` / `PointOnEdge` / `ParameterOnFace` / `PointOnFace` | `cxx:73-120` | the same names |
//!
//! `myExtCS` is `Extrema_ExtCS`
//! (`ModelingData/TKGeomBase/Extrema/Extrema_ExtCS.cxx:1-630`), ported as
//! [`occt_geom::extrema_cs::ExtCS`]; `Initialize` feeds it the face's surface,
//! its `BRepTools::UVBounds` window and the two tolerances (`cxx:50-57`), and
//! `Perform` drives it with the edge's curve over `BRep_Tool::Range(E, U1, U2)`
//! (`cxx:73-78`). [`ExtCS`] borrows its surface for the lifetime of the search,
//! so the engine is rebuilt around the window stored by [`ExtCF::initialize`]
//! on every `Perform` — `Extrema_ExtCS::Perform` resets every result sequence
//! (`cxx:108-122`), so the rebuilt engine carries exactly the state OCCT's
//! member would.
//!
//! Two behaviours are inherited, not emulated, and stay **UNPORTED** in the
//! `Extrema_ExtCS` family (see `occt_geom::extrema_cs::ext_cs`):
//!
//! * the `GeomAbs_Line` bounding-box clip of `[U1, U2]` (`Extrema_ExtCS.cxx:158-197`)
//!   needs `BndLib_AddSurface`, which lives in this crate; the generic arm then
//!   searches the caller's whole edge range, a superset of OCCT's window;
//! * the `Circ`/`Cylinder` and `Circ`/`Sphere` analytic arms of `Extrema_ExtElCS`
//!   answer `IsDone() == false` (they need the unported `Extrema_ExtElC`
//!   curve/curve arms), which is the `Extrema_ExtCS.cxx:298-415` fall-through to
//!   the generic engine.
//!
//! `P2.Parameter(U1, U2)` at `cxx:104` is the *getter*
//! (`Extrema_POnSurf.hxx:66-72`): it overwrites the `BRep_Tool::Range` values
//! with the surface parameters, and `gp_Pnt2d Puv(U1, U2)` (`cxx:105`) then
//! classifies the extrema's own `(u, v)` — which is what this port does.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d};
use occt_core::precision::Precision;
use occt_geom::extrema_cs::ExtCS;
use occt_geom::Surface;

use crate::brep_tool::BRepTool;
use crate::brep_top_adaptor_fclass2d::BrepTopAdaptorFClass2d;
use crate::fclass2d::FaceState;
use crate::shape::{Edge, Face};

/// `BRepExtrema_ExtCF` (`BRepExtrema_ExtCF.hxx:27-76`).
pub struct ExtCF {
    /// `myHS` (`hxx:75`): F's surface, `None` for the `OtherSurface` guard.
    my_hs: Option<Arc<dyn Surface>>,
    /// `myExtCS`'s `(U1, U2, V1, V2, aTolC, aTolS)` initialisation
    /// (`cxx:50-57`), replayed on every `Perform`.
    my_window: (f64, f64, f64, f64, f64, f64),
    /// `mySqDist` (`hxx:73`).
    my_sq_dist: Vec<f64>,
    /// `myPointsOnC` (`hxx:74`): `(u, point)` on the edge.
    my_points_on_c: Vec<(f64, GpPnt)>,
    /// `myPointsOnS` (`hxx:74`): `(u, v, point)` on the face.
    my_points_on_s: Vec<(f64, f64, GpPnt)>,
    /// `myExtCS.IsDone()` (`hxx:45`).
    my_done: bool,
    /// `myExtCS.IsParallel()` (`hxx:48`).
    my_parallel: bool,
}

impl ExtCF {
    /// `BRepExtrema_ExtCF(E, F)` (`cxx:30-33`).
    pub fn new(e: &Edge, f: &Face) -> Self {
        let mut this = Self {
            my_hs: None,
            my_window: (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            my_sq_dist: Vec::new(),
            my_points_on_c: Vec::new(),
            my_points_on_s: Vec::new(),
            my_done: false,
            my_parallel: false,
        };
        this.initialize(e, f);
        this.perform(e, f);
        this
    }

    /// `Initialize(E, F)` (`cxx:38-58`).
    pub fn initialize(&mut self, e: &Edge, f: &Face) {
        // `cxx:39-43`: `BRepAdaptor_Surface Surf(F)`, then the guard
        // `Surf.GetType() == GeomAbs_OtherSurface || !BRep_Tool::IsGeometric(E)`.
        let surf = match BRepTool::face_surface(f) {
            Some(s) => s,
            None => return,
        };
        if BRepTool::edge_curve(e).is_none() {
            return;
        }
        // `cxx:45-48`: `aTolS = min(Tol(F), Confusion)`, then
        // `min(UResolution, VResolution)`.
        let tol = BRepTool::face_tolerance(f).min(Precision::CONFUSION);
        let tol_s = occt_geom::approx_same_parameter::u_resolution(surf.as_ref(), tol)
            .min(occt_geom::approx_same_parameter::v_resolution(surf.as_ref(), tol))
            .max(Precision::PCONFUSION);
        // `cxx:50-53`: `aTolC = min(Tol(E), Confusion)`, then `aC.Resolution`.
        let tol_c = {
            let t = BRepTool::edge_tolerance(e).min(Precision::CONFUSION);
            match BRepTool::edge_curve(e) {
                Some(c) => c.resolution(t).max(Precision::PCONFUSION),
                None => Precision::PCONFUSION,
            }
        };
        // `cxx:55-57`: `BRepTools::UVBounds(F, U1, U2, V1, V2)` then
        // `myExtCS.Initialize(...)`.
        let (u1, u2, v1, v2) = BRepTool::uv_bounds(f);
        self.my_hs = Some(surf);
        self.my_window = (u1, u2, v1, v2, tol_c, tol_s);
    }

    /// `Perform(E, F2)` (`cxx:63-116`).
    pub fn perform(&mut self, e: &Edge, f2: &Face) {
        // `cxx:64-66`.
        self.my_sq_dist.clear();
        self.my_points_on_s.clear();
        self.my_points_on_c.clear();
        self.my_done = false;
        self.my_parallel = false;

        // `cxx:68-71`.
        let surf = match &self.my_hs {
            Some(s) => s.clone(),
            None => return,
        };
        let curve = match BRepTool::edge_curve(e) {
            Some(c) => c,
            None => return,
        };
        // `cxx:73-74`.
        let (u1, u2) = BRepTool::edge_parameters(e);

        // `myExtCS.Perform(*HC, U1, U2)` (`cxx:78`).
        let w = self.my_window;
        let mut ext = ExtCS::new();
        ext.initialize_range(surf.as_ref(), w.0, w.1, w.2, w.3, w.4, w.5);
        ext.perform(curve.as_ref(), u1, u2);

        // `cxx:80-83`.
        self.my_done = ext.is_done();
        if !self.my_done {
            return;
        }

        // `cxx:85-88`.
        self.my_parallel = ext.is_parallel();
        if self.my_parallel {
            self.my_sq_dist.push(ext.square_distance(1));
            return;
        }

        // `cxx:90-115`: classify every reported point.
        let tol = BRepTool::face_tolerance(f2);
        let classifier = BrepTopAdaptorFClass2d::new(f2, tol);
        // `cxx:100`: `isAdjustPeriodic = false`.
        let is_adjust_periodic = false;
        for i in 1..=ext.nb_ext() {
            let (p1, p2) = ext.points(i);
            // `cxx:104-105`: `P2.Parameter(U1, U2)` overwrites `U1` / `U2`,
            // so the classified 2D point is the extrema's own `(u, v)`.
            let (su, sv) = p2.parameter();
            let state = classifier.perform(GpPnt2d::new(su, sv), is_adjust_periodic);
            // `cxx:106-114`.
            if state == FaceState::On || state == FaceState::In {
                self.my_sq_dist.push(ext.square_distance(i));
                self.my_points_on_c.push((p1.parameter(), *p1.value()));
                self.my_points_on_s.push((su, sv, p2.value()));
            }
        }
    }

    /// `IsDone()` (`hxx:45`): `myExtCS.IsDone()`.
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// `IsParallel()` (`hxx:48`): `myExtCS.IsParallel()`, which raises
    /// `StdFail_NotDone` when the engine is not done (`Extrema_ExtCS.cxx:534-544`).
    pub fn is_parallel(&self) -> bool {
        assert!(self.my_done, "Extrema_ExtCS::IsParallel: StdFail_NotDone");
        self.my_parallel
    }

    /// `NbExt()` (`hxx:51`).
    pub fn nb_ext(&self) -> usize {
        self.my_sq_dist.len()
    }

    /// `SquareDistance(N)` (`hxx:54`).
    pub fn square_distance(&self, n: usize) -> f64 {
        self.my_sq_dist.get(n.wrapping_sub(1)).copied().unwrap_or(0.0)
    }

    /// `ParameterOnEdge(N)`.
    pub fn parameter_on_edge(&self, n: usize) -> f64 {
        self.my_points_on_c
            .get(n.wrapping_sub(1))
            .map(|p| p.0)
            .unwrap_or(0.0)
    }

    /// `PointOnEdge(N)`.
    pub fn point_on_edge(&self, n: usize) -> GpPnt {
        self.my_points_on_c
            .get(n.wrapping_sub(1))
            .map(|p| p.1)
            .unwrap_or_else(GpPnt::zero)
    }

    /// `ParameterOnFace(N, U, V)`.
    pub fn parameter_on_face(&self, n: usize) -> (f64, f64) {
        self.my_points_on_s
            .get(n.wrapping_sub(1))
            .map(|p| (p.0, p.1))
            .unwrap_or((0.0, 0.0))
    }

    /// `PointOnFace(N)`.
    pub fn point_on_face(&self, n: usize) -> GpPnt {
        self.my_points_on_s
            .get(n.wrapping_sub(1))
            .map(|p| p.2)
            .unwrap_or_else(GpPnt::zero)
    }
}
