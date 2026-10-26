//! `GeomFill_Sweep` — the geometric sweep built either as a KPart (an
//! elementary surface identified from the location / section laws) or as a
//! global approximation of the composed location x section function.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_Sweep.hxx:39-169` and
//! `GeomFill_Sweep.cxx` (OCCT 8.0.0): constructor `cxx:104-118`, `SetDomain`
//! `cxx:122-130`, `SetTolerance` `cxx:134-143`, `Build` `cxx:179-234`,
//! `Build2d` `cxx:236-243`, `BuildAll` `cxx:248-388`, `BuildKPart`
//! `cxx:577-1167`, accessors `cxx:1169-1251`.
//!
//! Translation notes:
//! - `mySec` / `myLoc` / `mySurface` / `myCurve2d` become `Option<Arc<..>>` /
//!   `Vec<Option<Arc<..>>>`; `CError(1..2, 1..n)` becomes `Vec<[f64; 2]>`
//!   (`[u_error, v_error]`).
//! - `GeomFill_Sweep::BuildProduct` is dead code in OCCT 8.0.0: its only call
//!   site is commented out (`cxx:216-219`) and its body leaves
//!   `mySurface = BSurf` with the tensor product commented out
//!   (`cxx:395-476`). It is therefore not ported; `isProduct` is always false,
//!   exactly as in OCCT.
//! - The `myForceApproxC1` block (`cxx:288-327`) runs
//!   `GeomConvert_ApproxSurface` (`crate::convert_approx_surface`) with
//!   `Tol = 1.e-4`, `C1` in U and V, degree 14, 16 segments, precision 1, and
//!   replaces the surface and the two 2d lines when it returns a result.
//!   In OCCT this block is reached only through `BRepFill_Pipe` /
//!   `BRepFill_Sweep` / `BRepOffsetAPI_MakePipeShell`, none of which are
//!   ported, so the branch has no reachable caller here yet.
//! - `KPartState` collects the locals of `BuildKPart` that the translation and
//!   rotation parts share (`S`, `Ok`, `UFirst`...); see `k_part.rs`.

use std::sync::Arc;

use occt_core::gp::{GpAx2d, GpDir2d, GpPnt, GpPnt2d};
use occt_core::kernel::geomabs::{CurveType, Shape};
use occt_geom2d::bspline_curve::Geom2dBSplineCurve;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::line::Geom2dLine;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use crate::approx_sweep_approximation::ApproxSweepApproximation;
use crate::bspline_surface::GeomBSplineSurface;
use crate::convert_approx_surface::GeomConvertApproxSurface;
use crate::geom_fill::location_law::LocationLaw;
use crate::geom_fill::section_law::SectionLaw;
use crate::geom_fill::sweep_function::SweepFunction;
use crate::surface::Surface;

mod k_part;

/// `GeomFill_ApproxStyle` (`GeomFill_ApproxStyle.hxx`): which function the
/// approximation is applied to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ApproxStyle {
    /// `GeomFill_Section`: the composed Function Location x Section.
    Section,
    /// `GeomFill_Location`: the location law, then algeabric composition.
    Location,
}

/// `GeomFill_Sweep` (`GeomFill_Sweep.hxx:39-169`).
pub struct Sweep {
    /// `First` (`hxx:149`).
    pub(crate) first: f64,
    /// `Last` (`hxx:150`).
    pub(crate) last: f64,
    /// `SFirst` (`hxx:151`).
    pub(crate) s_first: f64,
    /// `SLast` (`hxx:152`).
    pub(crate) s_last: f64,
    /// `Tol3d` (`hxx:153`).
    pub(crate) tol3d: f64,
    /// `BoundTol` (`hxx:154`).
    pub(crate) bound_tol: f64,
    /// `Tol2d` (`hxx:155`).
    pub(crate) tol2d: f64,
    /// `TolAngular` (`hxx:156`).
    pub(crate) tol_angular: f64,
    /// `SError` (`hxx:157`).
    pub(crate) s_error: f64,
    /// `myForceApproxC1` (`hxx:158`).
    force_approx_c1: bool,
    /// `myLoc` (`hxx:159`).
    pub(crate) loc: Arc<dyn LocationLaw>,
    /// `mySec` (`hxx:160`).
    pub(crate) sec: Option<Arc<dyn SectionLaw>>,
    /// `mySurface` (`hxx:161`).
    surface: Option<Arc<dyn Surface>>,
    /// `myCurve2d` (`hxx:162`), 0-based (OCCT's `Value(1)` is `[0]`).
    curve2d: Vec<Option<Arc<dyn Curve2d>>>,
    /// `CError` (`hxx:163`): `CError(1, i)` / `CError(2, i)` in `[0]` / `[1]`.
    c_error: Vec<[f64; 2]>,
    /// `done` (`hxx:164`).
    done: bool,
    /// `myExchUV` (`hxx:165`).
    exch_uv: bool,
    /// `isUReversed` (`hxx:166`).
    u_reversed: bool,
    /// `isVReversed` (`hxx:167`).
    v_reversed: bool,
    /// `myKPart` (`hxx:168`).
    k_part: bool,
}

impl Sweep {
    /// `GeomFill_Sweep(Location, WithKpart = true)` (`cxx:104-118`).
    pub fn new(location: Arc<dyn LocationLaw>) -> Self {
        Self::with_k_part(location, true)
    }

    /// `GeomFill_Sweep(Location, WithKpart)` (`cxx:104-118`).
    pub fn with_k_part(location: Arc<dyn LocationLaw>, with_k_part: bool) -> Self {
        let (first, last) = location.get_domain();
        Self {
            first,
            last,
            // `SFirst = SLast = 30.081996` marks "not set yet" (`cxx:115`).
            s_first: 30.081996,
            s_last: 30.081996,
            // `SetTolerance(1.e-4)` (`cxx:111`).
            tol3d: 1.0e-4,
            bound_tol: 1.0,
            tol2d: 1.0e-5,
            tol_angular: 1.0,
            // `SError = RealLast()` (`cxx:117`).
            s_error: f64::MAX,
            force_approx_c1: false,
            loc: location,
            sec: None,
            surface: None,
            curve2d: Vec::new(),
            c_error: Vec::new(),
            done: false,
            exch_uv: false,
            u_reversed: false,
            v_reversed: false,
            k_part: with_k_part,
        }
    }

    /// `GeomFill_Sweep::SetDomain` (`cxx:122-130`).
    pub fn set_domain(&mut self, loc_first: f64, loc_last: f64, section_first: f64, section_last: f64) {
        self.first = loc_first;
        self.last = loc_last;
        self.s_first = section_first;
        self.s_last = section_last;
    }

    /// `GeomFill_Sweep::SetTolerance` (`cxx:134-143`).
    pub fn set_tolerance(
        &mut self,
        tolerance3d: f64,
        bound_tolerance: f64,
        tolerance2d: f64,
        tolerance_angular: f64,
    ) {
        self.tol3d = tolerance3d;
        self.bound_tol = bound_tolerance;
        self.tol2d = tolerance2d;
        self.tol_angular = tolerance_angular;
    }

    /// `GeomFill_Sweep::SetForceApproxC1` (`cxx:152-155`).
    pub fn set_force_approx_c1(&mut self, force_approx_c1: bool) {
        self.force_approx_c1 = force_approx_c1;
    }

    /// `myForceApproxC1` (`hxx:158`).
    pub fn force_approx_c1(&self) -> bool {
        self.force_approx_c1
    }

    /// `GeomFill_Sweep::ExchangeUV` (`cxx:159-162`).
    pub fn exchange_uv(&self) -> bool {
        self.exch_uv
    }

    /// `GeomFill_Sweep::UReversed` (`cxx:164-167`).
    pub fn u_reversed(&self) -> bool {
        self.u_reversed
    }

    /// `GeomFill_Sweep::VReversed` (`cxx:171-174`).
    pub fn v_reversed(&self) -> bool {
        self.v_reversed
    }

    /// `GeomFill_Sweep::Build` (`cxx:179-234`).
    pub fn build(
        &mut self,
        section: Arc<dyn SectionLaw>,
        methode: ApproxStyle,
        continuity: Shape,
        degmax: i32,
        segmax: i32,
    ) {
        self.done = false;
        self.exch_uv = false;
        self.u_reversed = false;
        self.v_reversed = false;
        self.sec = Some(section);

        if self.s_first == self.s_last && self.s_last == 30.081996 {
            if let Some(sec) = &self.sec {
                let (f, l) = sec.get_domain();
                self.s_first = f;
                self.s_last = l;
            }
        }

        let mut is_k_part = false;
        // `isProduct` stays false: the only assignment is commented out in
        // OCCT (`cxx:216-219`).
        let is_product = false;

        if self.k_part {
            is_k_part = self.build_k_part();
        }

        if !is_k_part {
            self.exch_uv = false;
            self.u_reversed = false;
            self.v_reversed = false;
        }

        if !is_k_part && methode == ApproxStyle::Location {
            // `GeomFill_Sweep.cxx:213-219`: the body only reads
            // `mySec->BSplineSurface()` and then calls BuildProduct, which is
            // commented out there; nothing observable happens.
            if let Some(sec) = &self.sec {
                let _has_bspline_surface = sec.bspline_surface().is_some();
            }
        }

        if is_k_part || is_product {
            // Approximation of the 2d curves.
            self.done = self.build_2d(continuity, degmax, segmax);
        } else {
            // Global approximation.
            self.done = self.build_all(continuity, degmax, segmax);
        }
    }

    /// `GeomFill_Sweep::Build2d` (`cxx:236-243`): succeeds only when the
    /// location law carries no 2d curve.
    fn build_2d(&mut self, _continuity: Shape, _degmax: i32, _segmax: i32) -> bool {
        self.loc.nb_2d_curves() == 0
    }

    /// `GeomFill_Sweep::BuildAll` (`cxx:248-388`).
    fn build_all(&mut self, continuity: Shape, degmax: i32, segmax: i32) -> bool {
        let sec = match &self.sec {
            Some(s) => s.clone(),
            None => return false,
        };

        let func = Arc::new(SweepFunction::new(
            sec.clone(),
            self.loc.clone(),
            self.first,
            self.s_first,
            (self.s_last - self.s_first) / (self.last - self.first),
        ));
        let mut approx = ApproxSweepApproximation::new(func);
        approx.perform(
            self.first,
            self.last,
            self.tol3d,
            self.bound_tol,
            self.tol2d,
            self.tol_angular,
            continuity,
            degmax,
            segmax,
        );
        if !approx.is_done() {
            return false;
        }

        // `cxx:266-285`: the surface. OCCT's `Surface()` hands back
        // `tabPoles(1..Num3DSS, 1..NbPoles)` (section, pole); the B-spline
        // surface wants (pole inside the section, section).
        let (u_degree, v_degree, nb_u_poles, nb_v_poles, _nb_u_knots, _nb_v_knots) =
            match approx.surf_shape() {
                Some(v) => v,
                None => return false,
            };
        let poles_tab = approx.surf_poles();
        let weights_tab = approx.surf_weights();
        let mut poles = vec![vec![GpPnt::zero(); nb_v_poles]; nb_u_poles];
        let mut weights = vec![vec![1.0f64; nb_v_poles]; nb_u_poles];
        for v in 0..nb_v_poles {
            for u in 0..nb_u_poles {
                poles[u][v] = poles_tab[v][u];
                weights[u][v] = weights_tab[v][u];
            }
        }
        let knots_u = expand_knots(approx.surf_u_knots(), approx.surf_u_mults());
        let knots_v = expand_knots(approx.surf_v_knots(), approx.surf_v_mults());
        let rational = sec.is_rational();
        let mut bs = if rational {
            match GeomBSplineSurface::rational(
                poles,
                weights,
                knots_u,
                knots_v,
                u_degree as usize,
                v_degree as usize,
            ) {
                Ok(s) => s,
                Err(_) => return false,
            }
        } else {
            match GeomBSplineSurface::new(
                poles,
                knots_u,
                knots_v,
                u_degree as usize,
                v_degree as usize,
            ) {
                Ok(s) => s,
                Err(_) => return false,
            }
        };
        // `Geom_BSplineSurface(..., IsUPeriodic)` (`cxx:286`).
        bs.set_u_periodic_flag(sec.is_u_periodic());
        self.s_error = approx.max_error_on_surf().unwrap_or(self.s_error);

        // `cxx:288-327`: force a C1 approximation with
        // `GeomConvert_ApproxSurface` when the direct result is not `C1`.
        if self.force_approx_c1 && !matches!(bs.is_cnv(1), Ok(true)) {
            let the_tol = 1.0e-4;
            let convert = GeomConvertApproxSurface::new(
                &bs,
                the_tol,
                Shape::C1,
                Shape::C1,
                14,
                14,
                16,
                1,
            );
            if convert.has_result() {
                if let Some(converted) = convert.into_surface() {
                    // `cxx:304-307`: two lines of constant U at the first and
                    // the last U knot, both trimmed to the last V knot.
                    let u_first = converted.uknots().first().copied().unwrap_or(0.0);
                    let u_last = converted.uknots().last().copied().unwrap_or(0.0);
                    let v_last = converted.vknots().last().copied().unwrap_or(0.0);
                    if let Ok(dir) = GpDir2d::new(0.0, 1.0) {
                        let first = Geom2dTrimmedCurve::new(
                            Arc::new(Geom2dLine::new(GpAx2d::new(
                                GpPnt2d::new(u_first, 0.0),
                                dir,
                            ))),
                            0.0,
                            v_last,
                        );
                        let last = Geom2dTrimmedCurve::new(
                            Arc::new(Geom2dLine::new(GpAx2d::new(
                                GpPnt2d::new(u_last, 0.0),
                                dir,
                            ))),
                            0.0,
                            v_last,
                        );
                        // `myCurve2d = HArray1(1, 2)` (`cxx:302`) replaces the
                        // wider array of `cxx:330-386`.
                        self.curve2d = vec![Some(Arc::new(first)), Some(Arc::new(last))];
                        self.c_error = vec![[0.0, 0.0], [0.0, 0.0]];
                        bs = converted;
                        self.s_error = the_tol;
                    }
                }
            }
        }
        self.surface = Some(Arc::new(bs));

        // `cxx:330-386`: the 2d curves.
        if self.curve2d.is_empty() {
            let nb_2d = 2 + self.loc.trace_number();
            let mut curves: Vec<Option<Arc<dyn Curve2d>>> = vec![None; nb_2d.max(0) as usize];
            let mut errors = vec![[0.0f64; 2]; nb_2d.max(0) as usize];

            let ideb = if self.loc.has_first_restriction() { 1 } else { 2 };
            let mut ifin = 1 + self.loc.trace_number();
            if self.loc.has_last_restriction() {
                ifin += 1;
            }

            let mut kk = 1i32;
            let mut ii = ideb;
            while ii <= ifin {
                if let Some(p2d) = approx.curve_2d(kk as usize - 1) {
                    let degree = approx.curves_2d_degree().max(1) as usize;
                    let knots = expand_knots(approx.curves_2d_knots(), approx.curves_2d_mults());
                    let xs: Vec<f64> = p2d.iter().map(|p| p.x()).collect();
                    let ys: Vec<f64> = p2d.iter().map(|p| p.y()).collect();
                    if let Ok(c) = Geom2dBSplineCurve::new(xs, ys, knots, degree) {
                        curves[(ii - 1) as usize] = Some(Arc::new(c));
                        let err = approx.max_2d_error(kk as usize - 1);
                        errors[(ii - 1) as usize] = [err, err];
                    }
                }
                ii += 1;
                kk += 1;
            }

            // If the restriction curves were not computed, take the iso edges.
            if !self.loc.has_first_restriction() {
                let p = GpPnt2d::new(approx.surf_u_knots().first().copied().unwrap_or(0.0), 0.0);
                let line = Geom2dLine::new(GpAx2d::new(p, GpDir2d::new(0.0, 1.0).unwrap()));
                curves[0] = Some(Arc::new(Geom2dTrimmedCurve::new(
                    Arc::new(line),
                    self.first,
                    self.last,
                )));
                errors[0] = [0.0, 0.0];
            }
            if !self.loc.has_last_restriction() {
                let p = GpPnt2d::new(approx.surf_u_knots().last().copied().unwrap_or(0.0), 0.0);
                let line = Geom2dLine::new(GpAx2d::new(p, GpDir2d::new(0.0, 1.0).unwrap()));
                let last = curves.len() - 1;
                curves[last] = Some(Arc::new(Geom2dTrimmedCurve::new(
                    Arc::new(line),
                    self.first,
                    self.last,
                )));
                errors[last] = [0.0, 0.0];
            }

            self.curve2d = curves;
            self.c_error = errors;
        }

        true
    }

    /// `GeomFill_Sweep::IsDone` (`cxx:1169-1172`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `GeomFill_Sweep::ErrorOnSurface` (`cxx:1176-1179`).
    pub fn error_on_surface(&self) -> f64 {
        self.s_error
    }

    /// `GeomFill_Sweep::ErrorOnRestriction` (`cxx:1184-1198`): `(UError,
    /// VError)` of the first or last restriction.
    pub fn error_on_restriction(&self, is_first: bool) -> (f64, f64) {
        let ind = if is_first { 0 } else { self.c_error.len() - 1 };
        (self.c_error[ind][0], self.c_error[ind][1])
    }

    /// `GeomFill_Sweep::ErrorOnTrace` (`cxx:1203-1212`).
    ///
    /// `Standard_OutOfRange` (`cxx:1207-1210`) becomes `None`.
    pub fn error_on_trace(&self, index_of_trace: i32) -> Option<(f64, f64)> {
        let ind = index_of_trace + 1;
        if index_of_trace > self.loc.trace_number() {
            return None;
        }
        let e = self.c_error.get(ind as usize - 1)?;
        Some((e[0], e[1]))
    }

    /// `GeomFill_Sweep::Surface` (`cxx:1217-1220`).
    pub fn surface(&self) -> Option<Arc<dyn Surface>> {
        self.surface.clone()
    }

    /// `GeomFill_Sweep::Restriction` (`cxx:1225-1231`).
    pub fn restriction(&self, is_first: bool) -> Option<Arc<dyn Curve2d>> {
        if is_first {
            return self.curve2d.first().cloned().flatten();
        }
        self.curve2d.last().cloned().flatten()
    }

    /// `GeomFill_Sweep::NumberOfTrace` (`cxx:1236-1239`).
    pub fn number_of_trace(&self) -> i32 {
        self.loc.trace_number()
    }

    /// `GeomFill_Sweep::Trace` (`cxx:1244-1251`).
    ///
    /// `Standard_OutOfRange` (`cxx:1248-1251`) becomes `None`.
    pub fn trace(&self, index_of_trace: i32) -> Option<Arc<dyn Curve2d>> {
        let ind = index_of_trace + 1;
        if index_of_trace > self.loc.trace_number() {
            return None;
        }
        self.curve2d.get(ind as usize - 1).cloned().flatten()
    }
}

/// Expand OCCT's `(Knots, Mults)` pair into the flat knot vector the ported
/// B-spline constructors take (`Geom_BSplineSurface` / `Geom2d_BSplineCurve`
/// receive them separately).
pub(crate) fn expand_knots(knots: &[f64], mults: &[i32]) -> Vec<f64> {
    let mut flat = Vec::new();
    for (k, m) in knots.iter().zip(mults.iter()) {
        for _ in 0..(*m).max(0) {
            flat.push(*k);
        }
    }
    flat
}

/// `GeomAbs_CurveType` of an adaptor curve (`GeomAdaptor_Curve::GetType`,
/// `GeomAdaptor_Curve.cxx:382-410`): the concrete types the KPart search
/// branches on.
pub(crate) fn adaptor_curve_type(curve: &Arc<dyn crate::curve::Curve>) -> CurveType {
    if curve.is_line() {
        CurveType::Line
    } else if curve.gp_circ().is_some() {
        CurveType::Circle
    } else if curve.gp_ellipse().is_some() {
        CurveType::Ellipse
    } else if curve.gp_hyperbola().is_some() {
        CurveType::Hyperbola
    } else if curve.gp_parabola().is_some() {
        CurveType::Parabola
    } else if curve.bspline_poles().is_some() {
        CurveType::BSplineCurve
    } else if curve.bezier_poles().is_some() {
        CurveType::BezierCurve
    } else if curve.offset_curve().is_some() {
        CurveType::OffsetCurve
    } else {
        CurveType::OtherCurve
    }
}
