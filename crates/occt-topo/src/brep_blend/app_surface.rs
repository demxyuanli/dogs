//! `BRepBlend_AppSurface` — the `AppBlend_Approx` the facing builder stores a
//! blending surface from.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_AppSurface.hxx:33-122`,
//! `.cxx:26-135` and `.lxx:28-94` (OCCT 8.0.0). The `AppBlend_Approx` base is
//! `AppBlend_Approx.hxx:31-111`; every accessor below forwards to the wrapped
//! `Approx_SweepApproximation` exactly as the base declares it.
//!
//! Translation notes:
//! - `Approx_SweepApproximation approx` (`hxx:122`) is
//!   [`ApproxSweepApproximation`]. OCCT builds it in the member-initialisation
//!   list (`cxx:34`, `approx(Func)`) and only then runs the step (1) checks;
//!   that constructor stores the handle and nothing else
//!   (`Approx_SweepApproximation.cxx:60-68`), so the port runs the checks first
//!   and then builds the approximator.
//! - `NCollection_Array1/2` out-parameters become the slices the approximator
//!   already exposes ([`AppSurface::surface`], [`AppSurface::curve_2d`]), the
//!   same mapping `AppBlend_AppSurf` uses in `occt_geom::geom_fill`.
//! - `StdFail_NotDone` of the accessors becomes a panic, matching
//!   `AppSurf::require_done`, because the wrapped approximator reports
//!   "not done" with `Option`.
//! - `Dump(Standard_OStream&)` answers the text instead of writing it, the
//!   convention of the other `Dump` ports in the tree.
//! - `NCollection_Array1<gp_Pnt> P(1, NbPolSect)` / `... P2d(1, Nb2d)` of
//!   `cxx:47-53` only size the arrays handed to `D2` / `D1`; the port's
//!   `d1` / `d2` return their values, so the arrays and the `Nb2d == 0 -> 1`
//!   adjustment (`cxx:42-46`) have no counterpart.

use std::sync::Arc;

use occt_core::approx::ApproxSweepFunction;
use occt_core::gp::{GpPnt, GpPnt2d};
use occt_core::kernel::geomabs::Shape;
use occt_geom::ApproxSweepApproximation;

/// `BRepBlend_AppSurface` (`BRepBlend_AppSurface.hxx:33-122`).
pub struct AppSurface {
    /// `approx` (`hxx:122`).
    approx: ApproxSweepApproximation,
}

impl AppSurface {
    /// `BRepBlend_AppSurface(Funct, First, Last, Tol3d, Tol2d, TolAngular,
    /// Continuity, Degmax, Segmax)` (`cxx:24-73`).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        funct: Arc<dyn ApproxSweepFunction>,
        first: f64,
        last: f64,
        tol3d: f64,
        tol2d: f64,
        tol_angular: f64,
        continuity: Shape,
        degmax: i32,
        segmax: i32,
    ) -> Self {
        let mut continuity = continuity;

        // (1) Verification de la possibilite de derivation (`cxx:39-70`).
        if continuity != Shape::C0 {
            if continuity == Shape::C2 && funct.d2(first, first, last).is_err() {
                continuity = Shape::C1;
            }
            if continuity == Shape::C1 && funct.d1(first, first, last).is_err() {
                continuity = Shape::C0;
            }
        }

        // (2) Approximation (`cxx:72`): `Tol3d` is passed twice, for the
        // boundary tolerance and the surface tolerance.
        let mut approx = ApproxSweepApproximation::new(funct);
        approx.perform(
            first,
            last,
            tol3d,
            tol3d,
            tol2d,
            tol_angular,
            continuity,
            degmax,
            segmax,
        );
        Self { approx }
    }

    /// `IsDone()` (`lxx:28-31`).
    pub fn is_done(&self) -> bool {
        self.approx.is_done()
    }

    /// `SurfShape(UDegree, VDegree, NbUPoles, NbVPoles, NbUKnots, NbVKnots)`
    /// (`cxx:76-86`).
    pub fn surf_shape(&self) -> (i32, i32, i32, i32, i32, i32) {
        self.require_done();
        let (udeg, vdeg, nb_u, nb_v, nb_uknots, nb_vknots) = self
            .approx
            .surf_shape()
            .expect("StdFail_NotDone: BRepBlend_AppSurface");
        (
            udeg,
            vdeg,
            nb_u as i32,
            nb_v as i32,
            nb_uknots as i32,
            nb_vknots as i32,
        )
    }

    /// `Surface(TPoles, TWeights, TUKnots, TVKnots, TUMults, TVMults)`
    /// (`cxx:88-96`): one accessor per out-parameter, in the same order.
    pub fn surface(
        &self,
    ) -> (
        &[Vec<GpPnt>],
        &[Vec<f64>],
        &[f64],
        &[f64],
        &[i32],
        &[i32],
    ) {
        self.require_done();
        (
            self.approx.surf_poles(),
            self.approx.surf_weights(),
            self.approx.surf_u_knots(),
            self.approx.surf_v_knots(),
            self.approx.surf_u_mults(),
            self.approx.surf_v_mults(),
        )
    }

    /// `UDegree()` (`lxx:32-35`).
    pub fn u_degree(&self) -> i32 {
        self.require_done();
        self.approx.u_degree()
    }

    /// `VDegree()` (`lxx:37-40`).
    pub fn v_degree(&self) -> i32 {
        self.require_done();
        self.approx.v_degree()
    }

    /// `SurfPoles()` (`lxx:42-45`).
    pub fn surf_poles(&self) -> &[Vec<GpPnt>] {
        self.require_done();
        self.approx.surf_poles()
    }

    /// `SurfWeights()` (`lxx:47-50`).
    pub fn surf_weights(&self) -> &[Vec<f64>] {
        self.require_done();
        self.approx.surf_weights()
    }

    /// `SurfUKnots()` (`lxx:52-55`).
    pub fn surf_u_knots(&self) -> &[f64] {
        self.approx.surf_u_knots()
    }

    /// `SurfVKnots()` (`lxx:57-60`).
    pub fn surf_v_knots(&self) -> &[f64] {
        self.approx.surf_v_knots()
    }

    /// `SurfUMults()` (`lxx:62-65`).
    pub fn surf_u_mults(&self) -> &[i32] {
        self.approx.surf_u_mults()
    }

    /// `SurfVMults()` (`lxx:67-70`).
    pub fn surf_v_mults(&self) -> &[i32] {
        self.approx.surf_v_mults()
    }

    /// `MaxErrorOnSurf()` (`cxx:98-101`): the maximum error in the surface
    /// approximation.
    pub fn max_error_on_surf(&self) -> f64 {
        self.approx
            .max_error_on_surf()
            .expect("StdFail_NotDone: BRepBlend_AppSurface")
    }

    /// `NbCurves2d()` (`lxx:71-74`).
    pub fn nb_curves_2d(&self) -> i32 {
        self.approx.nb_curves_2d() as i32
    }

    /// `Curves2dShape(Degree, NbPoles, NbKnots)` (`cxx:103-107`).
    pub fn curves_2d_shape(&self) -> (i32, i32, i32) {
        self.require_done();
        let (degree, nb_poles, nb_knots) = self
            .approx
            .curves_2d_shape()
            .expect("Standard_DomainError: BRepBlend_AppSurface::Curves2dShape");
        (degree, nb_poles as i32, nb_knots as i32)
    }

    /// `Curve2d(Index, TPoles, TKnots, TMults)` (`cxx:109-114`), `Index`
    /// `1`-based.
    pub fn curve_2d(&self, index: i32) -> (&[GpPnt2d], &[f64], &[i32]) {
        self.require_done();
        (
            self.curve_2d_poles(index),
            self.approx.curves_2d_knots(),
            self.approx.curves_2d_mults(),
        )
    }

    /// `Curves2dDegree()` (`lxx:76-79`).
    pub fn curves_2d_degree(&self) -> i32 {
        self.require_done();
        self.approx.curves_2d_degree()
    }

    /// `Curve2dPoles(Index)` (`lxx:81-85`), `Index` `1`-based.
    pub fn curve_2d_poles(&self, index: i32) -> &[GpPnt2d] {
        self.require_done();
        self.approx
            .curve_2d_poles((index - 1).max(0) as usize)
            .expect("Standard_OutOfRange: BRepBlend_AppSurface::Curve2dPoles")
    }

    /// `Curves2dKnots()` (`lxx:87-90`).
    pub fn curves_2d_knots(&self) -> &[f64] {
        self.require_done();
        self.approx.curves_2d_knots()
    }

    /// `Curves2dMults()` (`lxx:92-95`).
    pub fn curves_2d_mults(&self) -> &[i32] {
        self.require_done();
        self.approx.curves_2d_mults()
    }

    /// `TolReached(Tol3d, Tol2d)` (`cxx:121-132`): `Tol2d` is the largest
    /// `Max2dError` over the curves.
    pub fn tol_reached(&self) -> (f64, f64) {
        let tol3d = self.max_error_on_surf();
        let mut tol2d = 0.0f64;
        for ii in 1..=self.nb_curves_2d() {
            tol2d = tol2d.max(self.max_2d_error(ii));
        }
        (tol3d, tol2d)
    }

    /// `Max2dError(Index)` (`cxx:116-119`), `Index` `1`-based.
    pub fn max_2d_error(&self, index: i32) -> f64 {
        self.approx.max_2d_error((index - 1).max(0) as usize)
    }

    /// `TolCurveOnSurf(Index)` (`cxx:130-134`), `Index` `1`-based.
    pub fn tol_curve_on_surf(&self, index: i32) -> f64 {
        self.approx
            .tol_curve_on_surf((index - 1).max(0) as usize)
    }

    /// `Dump(Standard_OStream&)` (`cxx:135-137`): the text
    /// `Approx_SweepApproximation::Dump` writes.
    pub fn dump(&self) -> String {
        self.approx.dump()
    }

    /// `throw StdFail_NotDone()` of the `AppBlend_Approx` accessors
    /// (`AppBlend_Approx.hxx:38`, via the wrapped approximator's own checks).
    fn require_done(&self) {
        if !self.approx.is_done() {
            panic!("StdFail_NotDone: BRepBlend_AppSurface");
        }
    }
}
