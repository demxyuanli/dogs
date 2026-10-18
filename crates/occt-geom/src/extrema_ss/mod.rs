//! Surface-surface extrema. Port of `Extrema_ExtSS`, `Extrema_ExtElSS`,
//! `Extrema_GenExtSS`, `Extrema_FuncExtSS` (TKGeomBase).
//!
//! `Extrema_ExtSS::Perform` (`Extrema_ExtSS.cxx:100-255`) switches on
//! `S1.GetType()` (`cxx:117-119`):
//!
//! * `GeomAbs_Plane` x `GeomAbs_Plane` -> `Extrema_ExtElSS::Perform(S1.Plane(),
//!   myS2->Plane())` (`cxx:122-129`), the only analytic arm;
//! * `GeomAbs_Plane` x anything else -> `Extrema_GenExtSS(S1, *myS2, NbU, NbV,
//!   TolS1, TolS2)` with `NbU = NbV = 20` (`cxx:130-175`);
//! * anything else (the `default:` arm) -> the same `Extrema_GenExtSS`
//!   (`cxx:179-224`).
//!
//! Both `GenExtSS` arms then normalise the reported parameters into the period
//! of a periodic surface (`ElCLib::InPeriod`) and keep only the solutions that
//! lie inside the stored window within `TolS1` / `TolS2`
//! (`cxx:136-171` / `cxx:188-219`).
//!
//! OCCT has **no** analytic sphere-sphere or sphere-plane surface extrema:
//! `Extrema_ExtElSS::Perform(const gp_Pln&, const gp_Sphere&)`
//! (`Extrema_ExtElSS.cxx:62-69`) and `Perform(const gp_Sphere&, const gp_Sphere&)`
//! (`cxx:77-83`) set the flags and then `throw Standard_NotImplemented()`. An
//! earlier revision of this port carried closed-form arms for those pairs; they
//! had no `.cxx` counterpart and are gone.
//!
//! The consumer-facing entry point of a `Geom_Surface` pair is
//! `GeomAPI_ExtremaSurfaceSurface`, ported as
//! [`crate::geom_api::ExtremaSurfaceSurface`] (`GeomAPI_ExtremaSurfaceSurface.cxx:57-231`),
//! which keeps the index of the smallest `SquareDistance` in `myIndex`
//! (`cxx:75-89`) and reads it from `NearestPoints` / `LowerDistanceParameters`
//! / `LowerDistance`. `ExtremaExtSS` itself exposes only `Points(N, P1, P2)`
//! (`Extrema_ExtSS.hxx:92`, `Extrema_ExtSS.cxx:292-300`); `myPOnS1` /
//! `myPOnS2` remain private (`hxx:99-102`).

mod ext_el_ss;
mod func_ext_ss;
mod gen_ext_ss;

pub use ext_el_ss::ExtremaExtElSS;
pub use func_ext_ss::{ExtremaFuncDistSS, ExtremaFuncExtSS};
pub use gen_ext_ss::ExtremaGenExtSS;

use occt_core::elib::clib;
use occt_core::gp::GpPnt;
use occt_core::precision::{INFINITE, PCONFUSION};

use crate::extrema::ExtremaPair;
use crate::extrema_surf::{classify_plane, ExtremaPOnSurf};
use crate::surface::Surface;

/// `NbU` / `NbV` of `Extrema_ExtSS::Perform` (`Extrema_ExtSS.cxx:118`).
const NB_SAMPLE: usize = 20;

/// An adaptor parameter range as OCCT sees it.
///
/// OCCT's unbounded elementary surfaces answer `±Precision::Infinite()`
/// (`Geom_Plane::Bounds`, `Geom_Plane.cxx:178-185`; `Precision.hxx:371`,
/// `2.e100`), and `Extrema_ExtSS` feeds those numbers straight into
/// `Extrema_GenExtSS`'s sample spacing (`Extrema_ExtSS.cxx:131-180` ->
/// `Extrema_GenExtSS.cxx:190-197`). The port encodes the same span as `±inf`
/// (`plane.rs:42-43`, `cylinder.rs:49`, `cone.rs:58`,
/// `surface_of_linear_extrusion.rs:160`), so translate it back before
/// arithmetic.
pub(super) fn occt_range(r: (f64, f64)) -> (f64, f64) {
    (
        if r.0.is_finite() { r.0 } else { -INFINITE },
        if r.1.is_finite() { r.1 } else { INFINITE },
    )
}

/// `Extrema_ExtSS` (`Extrema_ExtSS.hxx:30-140`).
///
/// The three result sequences are kept separate, as in OCCT: the
/// parallel-plane arm appends to `mySqDist` only (`Extrema_ExtSS.cxx:228-235`),
/// so `NbExt()` can be 1 while `myPOnS1` / `myPOnS2` hold nothing.
#[derive(Debug, Clone)]
pub struct ExtremaExtSS {
    /// `myDone`.
    done: bool,
    /// `myIsPar`.
    is_par: bool,
    uinf1: f64,
    usup1: f64,
    vinf1: f64,
    vsup1: f64,
    uinf2: f64,
    usup2: f64,
    vinf2: f64,
    vsup2: f64,
    tol_s1: f64,
    tol_s2: f64,
    /// `mySqDist`.
    sq_dist: Vec<f64>,
    /// `myPOnS1`.
    p_on_s1: Vec<ExtremaPOnSurf>,
    /// `myPOnS2`.
    p_on_s2: Vec<ExtremaPOnSurf>,
}

impl Default for ExtremaExtSS {
    /// `Extrema_ExtSS()` (`Extrema_ExtSS.cxx:26-42`).
    fn default() -> Self {
        ExtremaExtSS {
            done: false,
            is_par: false,
            uinf1: 0.0,
            usup1: 0.0,
            vinf1: 0.0,
            vsup1: 0.0,
            uinf2: 0.0,
            usup2: 0.0,
            vinf2: 0.0,
            vsup2: 0.0,
            tol_s1: 0.0,
            tol_s2: 0.0,
            sq_dist: Vec::new(),
            p_on_s1: Vec::new(),
            p_on_s2: Vec::new(),
        }
    }
}

impl ExtremaExtSS {
    /// `Extrema_ExtSS(S1, S2, TolS1, TolS2)` (`Extrema_ExtSS.cxx:44-61`):
    /// `Initialize(S2, S2.First..LastVParameter(), TolS2)` then
    /// `Perform(S1, S1.First..LastVParameter(), TolS1)`.
    pub fn new(s1: &dyn Surface, s2: &dyn Surface, tol_s1: f64, tol_s2: f64) -> Self {
        let mut me = ExtremaExtSS::default();
        let (uinf2, usup2) = occt_range(s2.u_range());
        let (vinf2, vsup2) = occt_range(s2.v_range());
        me.initialize(s2, uinf2, usup2, vinf2, vsup2, tol_s2);
        let (uinf1, usup1) = occt_range(s1.u_range());
        let (vinf1, vsup1) = occt_range(s1.v_range());
        me.perform(s1, s2, uinf1, usup1, vinf1, vsup1, tol_s1);
        me
    }

    /// `Extrema_ExtSS::Initialize(S2, Uinf2, Usup2, Vinf2, Vsup2, TolS2)`
    /// (`Extrema_ExtSS.cxx:84-98`): stores the window and tolerance of S2 and
    /// clears `myIsPar`.
    pub fn initialize(
        &mut self,
        _s2: &dyn Surface,
        uinf2: f64,
        usup2: f64,
        vinf2: f64,
        vsup2: f64,
        tol_s2: f64,
    ) {
        self.is_par = false;
        self.uinf2 = uinf2;
        self.usup2 = usup2;
        self.vinf2 = vinf2;
        self.vsup2 = vsup2;
        self.tol_s2 = tol_s2;
    }

    /// `Extrema_ExtSS::Perform(S1, Uinf1, Usup1, Vinf1, Vsup1, TolS1)`
    /// (`Extrema_ExtSS.cxx:100-255`).
    #[allow(clippy::too_many_arguments)]
    pub fn perform(
        &mut self,
        s1: &dyn Surface,
        s2: &dyn Surface,
        uinf1: f64,
        usup1: f64,
        vinf1: f64,
        vsup1: f64,
        tol_s1: f64,
    ) {
        self.uinf1 = uinf1;
        self.usup1 = usup1;
        self.vinf1 = vinf1;
        self.vsup1 = vsup1;
        self.tol_s1 = tol_s1;
        self.sq_dist.clear();
        self.p_on_s1.clear();
        self.p_on_s2.clear();

        // switch (myS1type) — Extrema_ExtSS.cxx:117-179.
        let s1_plane = classify_plane(s1);
        let s2_plane = classify_plane(s2);
        if let (Some(p1), Some(p2)) = (&s1_plane, &s2_plane) {
            // case GeomAbs_Plane: case GeomAbs_Plane: cxx:122-129.
            let ext_el = ExtremaExtElSS::new_plane_plane(p1, p2);
            self.done = ext_el.is_done();
            if self.done {
                self.is_par = ext_el.is_parallel();
                if self.is_par {
                    // cxx:229-233: only the square distance is appended.
                    self.sq_dist.push(ext_el.square_distance(1));
                } else {
                    // cxx:235-253.
                    for i in 1..=ext_el.nb_ext() {
                        // Extrema_ExtElSS::Perform(gp_Pln, gp_Pln) never fills
                        // myPOnS1 / myPOnS2 (cxx:50-52), so NbExt() is 0 here
                        // and this loop does not run.
                        let _ = i;
                    }
                }
            }
            return;
        }

        // default: Extrema_GenExtSS(S1, *myS2, NbU, NbV, mytolS1, mytolS2)
        // — cxx:131-175 (Plane x non-Plane) and cxx:180-221 (non-Plane).
        let gen = ExtremaGenExtSS::with_surfaces(s1, s2, NB_SAMPLE, NB_SAMPLE, tol_s1, self.tol_s2);
        self.done = gen.is_done();
        if !self.done {
            return;
        }
        let nb_ext = gen.nb_ext();
        for i in 1..=nb_ext {
            let ps1 = gen.point_on_s1(i as i32);
            let ps2 = gen.point_on_s2(i as i32);
            let (mut u1, mut v1) = ps1.parameter();
            let (mut u2, mut v2) = ps2.parameter();
            if s1.is_u_periodic() {
                u1 = clib::in_period(u1, self.uinf1, self.uinf1 + s1.u_period());
            }
            if s1.is_v_periodic() {
                v1 = clib::in_period(v1, self.vinf1, self.vinf1 + s1.v_period());
            }
            if s2.is_u_periodic() {
                u2 = clib::in_period(u2, self.uinf2, self.uinf2 + s2.u_period());
            }
            if s2.is_v_periodic() {
                v2 = clib::in_period(v2, self.vinf2, self.vinf2 + s2.v_period());
            }
            if (self.uinf1 - u1) <= self.tol_s1
                && (u1 - self.usup1) <= self.tol_s1
                && (self.vinf1 - v1) <= self.tol_s1
                && (v1 - self.vsup1) <= self.tol_s1
                && (self.uinf2 - u2) <= self.tol_s2
                && (u2 - self.usup2) <= self.tol_s2
                && (self.vinf2 - v2) <= self.tol_s2
                && (v2 - self.vsup2) <= self.tol_s2
            {
                self.sq_dist.push(gen.square_distance(i as i32));
                self.p_on_s1.push(ExtremaPOnSurf::new(u1, v1, ps1.value()));
                self.p_on_s2.push(ExtremaPOnSurf::new(u2, v2, ps2.value()));
            }
        }
    }

    /// `IsDone()` (`Extrema_ExtSS.cxx:257-260`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `IsParallel()` (`Extrema_ExtSS.cxx:262-271`): throws `StdFail_NotDone`
    /// when `!IsDone()`.
    pub fn is_parallel(&self) -> bool {
        assert!(self.is_done(), "Extrema_ExtSS::IsParallel: StdFail_NotDone");
        self.is_par
    }

    /// `NbExt()` (`Extrema_ExtSS.cxx:282-287`): `mySqDist.Length()`.
    pub fn nb_ext(&self) -> usize {
        assert!(self.is_done(), "Extrema_ExtSS::NbExt: StdFail_NotDone");
        self.sq_dist.len()
    }

    /// `SquareDistance(N)` (`Extrema_ExtSS.cxx:273-280`).
    pub fn square_distance(&self, n: i32) -> f64 {
        assert!(
            n >= 1 && n as usize <= self.nb_ext(),
            "Extrema_ExtSS::SquareDistance: Standard_OutOfRange"
        );
        self.sq_dist[(n - 1) as usize]
    }

    /// `Points(N, P1, P2)` (`Extrema_ExtSS.cxx:292-300`) — the only public
    /// accessor of `myPOnS1` / `myPOnS2` (`Extrema_ExtSS.hxx:92`; the two
    /// sequences themselves are private, `hxx:99-102`).
    ///
    /// The parallel-plane arm never allocates `myPOnS1` / `myPOnS2`
    /// (`Extrema_ExtElSS.cxx:50-52`), exactly as `NbExt()` there counts only
    /// `mySqDist`; indexing past the sequences reports the same condition
    /// OCCT raises `Standard_OutOfRange` for.
    pub fn points(&self, n: i32) -> (ExtremaPOnSurf, ExtremaPOnSurf) {
        assert!(
            n >= 1 && n as usize <= self.nb_ext(),
            "Extrema_ExtSS::Points: Standard_OutOfRange"
        );
        (
            self.p_on_s1[(n - 1) as usize],
            self.p_on_s2[(n - 1) as usize],
        )
    }
}

/// Surface-surface pair: `p1` on `s1` (params `u1`,`v1`), `p2` on `s2`
/// (params `u2`,`v2`).
fn ss_pair(p1: GpPnt, u1: f64, v1: f64, p2: GpPnt, u2: f64, v2: f64, d: f64) -> ExtremaPair {
    ExtremaPair {
        p1,
        p2,
        distance: d,
        u1,
        v1: Some(v1),
        u2,
        v2: Some(v2),
    }
}

/// Every surface-surface extremum that carries both points.
///
/// The `GeomAPI_ExtremaSurfaceSurface` counterpart of this enumeration is the
/// `Points(i, P1, P2)` loop over `i = 1..=NbExtrema()`
/// (`GeomAPI_ExtremaSurfaceSurface.cxx:162-172`, `cxx:147-157`); the
/// `myIndex`-based answers (`NearestPoints`, `LowerDistance`, ...) live in
/// [`crate::geom_api::ExtremaSurfaceSurface`].
///
/// The parallel-plane arm reports `NbExt() == 1` with no allocated points
/// (`Extrema_ExtSS.cxx:228-235`), so it contributes no entry here; call
/// [`ExtremaExtSS`] directly to observe `is_parallel()` /
/// `square_distance(1)` for that case.
pub fn surface_surface_extrema_all(s1: &dyn Surface, s2: &dyn Surface) -> Vec<ExtremaPair> {
    // The OCCT entry point for a `Geom_Surface` pair is
    // `GeomAPI_ExtremaSurfaceSurface::Init`, which runs
    // `Extrema_ExtSS(TheSurface1, TheSurface2, Tol, Tol)` with
    // `Tol = Precision::PConfusion()` (`GeomAPI_ExtremaSurfaceSurface.cxx:60-69`).
    let ext = ExtremaExtSS::new(s1, s2, PCONFUSION, PCONFUSION);
    if !ext.is_done() {
        return Vec::new();
    }
    (1..=ext.nb_ext())
        .filter_map(|n| {
            let p1 = ext.p_on_s1.get(n - 1)?;
            let p2 = ext.p_on_s2.get(n - 1)?;
            let (u1, v1) = p1.parameter();
            let (u2, v2) = p2.parameter();
            Some(ss_pair(
                p1.value(),
                u1,
                v1,
                p2.value(),
                u2,
                v2,
                // `Extrema_ExtSS::SquareDistance(N)`; `ExtremaPair::distance`
                // in this port is the Euclidean distance.
                ext.square_distance(n as i32).sqrt(),
            ))
        })
        .collect()
}

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
