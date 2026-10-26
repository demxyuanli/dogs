//! `Approx_CurveOnSurface` — approximation of a curve lying on a surface.
//!
//! Source: `ModelingData/TKGeomBase/Approx/Approx_CurveOnSurface.cxx`
//! (`Perform` `cxx:355-551`, `isIsoLine` `cxx:591-677`,
//! `buildC3dOnIsoLine` `cxx:682-811`), `.hxx:25-140`.
//!
//! The 2D curve is an `Adaptor2d_Curve2d` in OCCT; here it is an
//! `occt_geom2d::Curve2d`, whose type queries are read through
//! [`occt_geom2d::adaptor_curve::type2d`] (the `Geom2dAdaptor_Curve::load`
//! table). `Trim(First, Last, PConfusion)` (`cxx:387`) is
//! [`TrimmedAdaptor2d`], which — like `Geom2dAdaptor_Curve::Trim`
//! (`Geom2dAdaptor_Curve.cxx:577-584`) — keeps the curve and moves only the
//! reported range.
//!
//! UNPORTED, both inside `Perform`:
//! - `Adaptor3d_HSurfaceTool::IsSurfG1` (`Adaptor3d_HSurfaceTool.cxx:118-191`),
//!   read only when `UContinuity()` / `VContinuity()` is `GeomAbs_C0`
//!   (`cxx:436-458`). Its single `true` arm is `Geom_BSplineCurve::IsG1`
//!   (`Geom_BSplineCurve_1.cxx:64-...`) over a B-spline iso curve, and that
//!   curve predicate (with `BSplCLib::LocalD1`) is not translated. OCCT returns
//!   `false` (`cxx:190`) for every surface without such an iso curve, which is
//!   the answer taken here; the two `min(1.e-3, ...)` clamps of `cxx:438-458`
//!   are applied unconditionally instead of behind the predicate.
//! - the `fonct->Trim(StartEnd[0], StartEnd[1], PConfusion)` re-trim of the
//!   evaluators (`cxx:91-96`, `:188-192`, `:270-275`) is not needed:
//!   `Adaptor3d_CurveOnSurface::Trim` (`Adaptor3d_CurveOnSurface.cxx:1133-1141`)
//!   and `Geom2dAdaptor_Curve::Trim` only move the domain, they do not
//!   reparameterise, so the value at `par` is the same as on the full curve.

use std::sync::Arc;

use occt_core::adv_approx::{ApproxAFunction1dPair3d, ApproxAFunction3d, DichoCutting, PrefAndRec};
use occt_core::bspl::knots::reparametrize;
use occt_core::convert::conic_to_bspline::ParameterisationType;
use occt_core::gp::{GpDir2d, GpVec2d};
use occt_core::kernel::geomabs::{CurveType, Shape};
use occt_core::precision::{ANGULAR, CONFUSION, PCONFUSION};
use occt_geom2d::adaptor_curve::type2d;
use occt_geom2d::adaptor_trim::TrimmedAdaptor2d;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::Geom2dBSplineCurve;

use crate::approx_same_parameter::{u_resolution, v_resolution};
use crate::bspline_curve::GeomBSplineCurve;
use crate::convert_bspl::curve_to_bspline_curve;
use crate::curve::Curve;
use crate::curve_on_surface::CurveOnSurface;
use crate::rectangular_trimmed::GeomRectangularTrimmedSurface;
use crate::surface::Surface;
use crate::trimmed::GeomTrimmedCurveBasis;

/// `Approx_CurveOnSurface` (`hxx:25-140`).
pub struct ApproxCurveOnSurface {
    c2d: Arc<dyn Curve2d>,
    surf: Arc<dyn Surface>,
    first: f64,
    last: f64,
    tol: f64,
    curve2d: Option<Geom2dBSplineCurve>,
    curve3d: Option<GeomBSplineCurve>,
    is_done: bool,
    has_result: bool,
    error3d: f64,
    error2du: f64,
    error2dv: f64,
}

/// `(isU, param, isForward)` as `isIsoLine` reports it (`hxx:97`).
type IsoLine = (bool, f64, bool);

impl ApproxCurveOnSurface {
    /// `Approx_CurveOnSurface(theC2D, theSurf, theFirst, theLast, theTol)`
    /// (`cxx:336-352`): stores the inputs, `Perform` is the caller's.
    pub fn new(
        c2d: Arc<dyn Curve2d>,
        surf: Arc<dyn Surface>,
        first: f64,
        last: f64,
        tol: f64,
    ) -> Self {
        Self {
            c2d,
            surf,
            first,
            last,
            tol,
            curve2d: None,
            curve3d: None,
            is_done: false,
            has_result: false,
            error3d: 0.0,
            error2du: 0.0,
            error2dv: 0.0,
        }
    }

    /// `IsDone` (`cxx:553-556`).
    pub fn is_done(&self) -> bool {
        self.is_done
    }

    /// `HasResult` (`cxx:558-561`).
    pub fn has_result(&self) -> bool {
        self.has_result
    }

    /// `Curve3d` (`cxx:563-566`).
    pub fn curve3d(&self) -> Option<&GeomBSplineCurve> {
        self.curve3d.as_ref()
    }

    /// `Curve2d` (`cxx:568-571`).
    pub fn curve2d(&self) -> Option<&Geom2dBSplineCurve> {
        self.curve2d.as_ref()
    }

    /// `MaxError3d` (`cxx:573-576`).
    pub fn max_error3d(&self) -> f64 {
        self.error3d
    }

    /// `MaxError2dU` (`cxx:578-581`).
    pub fn max_error2du(&self) -> f64 {
        self.error2du
    }

    /// `MaxError2dV` (`cxx:583-586`).
    pub fn max_error2dv(&self) -> f64 {
        self.error2dv
    }

    /// `Perform(theMaxSegments, theMaxDegree, theContinuity, theOnly3d,
    /// theOnly2d)` (`cxx:355-551`).
    pub fn perform(
        &mut self,
        max_segments: i32,
        max_degree: i32,
        continuity: u8,
        only3d: bool,
        only2d: bool,
    ) -> Result<(), &'static str> {
        self.is_done = false;
        self.has_result = false;
        self.error2du = 0.0;
        self.error2dv = 0.0;
        self.error3d = 0.0;
        self.curve2d = None;
        self.curve3d = None;

        // `cxx:367-370`.
        if only3d && only2d {
            return Err("Standard_ConstructionError");
        }
        // `cxx:372-386`: G1 -> C1, G2 -> C2, above C2 -> C2 (the
        // `AdvApprox_ApproxAFunction` restriction).
        let a_cont: u8 = match continuity {
            1 => 2,
            3 => 4,
            c if c > 4 => 4,
            c => c,
        };

        // `cxx:387`: `myC2D->Trim(myFirst, myLast, Precision::PConfusion())`.
        let trimmed: Arc<TrimmedAdaptor2d> =
            Arc::new(TrimmedAdaptor2d::new(self.c2d.clone(), self.first, self.last));

        // `cxx:388-399`.
        if only3d {
            if let Some((is_u, param, is_forward)) = Self::is_iso_line(trimmed.as_ref()) {
                if self.build_c3d_on_iso_line(trimmed.as_ref(), is_u, param, is_forward)? {
                    self.is_done = true;
                    self.has_result = true;
                    return Ok(());
                }
            }
        }

        // `cxx:401`: `Adaptor3d_CurveOnSurface HCOnS(TrimmedC2D, mySurf)`.
        let hcons = CurveOnSurface::new(trimmed.clone(), self.surf.clone());

        // `cxx:427-469`: 2d tolerances, halved surface resolutions, and the
        // `GeomAbs_C0` clamps. `Adaptor3d_Surface::UContinuity` /
        // `VContinuity` share one value in this port (see
        // `crate::curve_on_surface::CurveOnSurface::continuity`).
        let mut tol_u = u_resolution(self.surf.as_ref(), self.tol) / 2.0;
        let mut tol_v = v_resolution(self.surf.as_ref(), self.tol) / 2.0;
        if self.surf.continuity() == 0 {
            // `cxx:436-447` (`IsSurfG1(mySurf, true, Angular/Confusion)`), see
            // the module note: the predicate is UNPORTED and answers false.
            tol_u = 1.0e-3_f64.min(1.0e3 * tol_u);
            tol_u = 1.0e-3_f64.min(1.0e2 * tol_u);
            // `cxx:449-459`.
            tol_v = 1.0e-3_f64.min(1.0e3 * tol_v);
            tol_v = 1.0e-3_f64.min(1.0e2 * tol_v);
        }
        // `cxx:463-470`: `ThreeDTol->Init(myTol / 2)`.
        let tol_3d = self.tol / 2.0;

        // `cxx:472-500`: the cutter depends on where the requested continuity
        // exceeds the curve / surface continuity.
        let c2d_cont = trimmed.continuity();
        let surf_cont = self.surf.continuity();
        let dicho = a_cont <= c2d_cont && a_cont <= surf_cont;
        let cut_pref: Option<PrefAndRec> = if dicho {
            None
        } else if a_cont == 2 {
            // `cxx:479-489`: PrefAndRec(CutPnts_C1, CutPnts_C2), weight 5.
            let cut_c1 = hcons.parameter_intervals(2);
            let cut_c2 = hcons.parameter_intervals(4);
            Some(PrefAndRec::new(cut_c1, cut_c2, 5.0).map_err(|_| "AdvApprox_PrefAndRec")?)
        } else {
            // `cxx:490-500`: PrefAndRec(CutPnts_C2, CutPnts_C3).
            let cut_c2 = hcons.parameter_intervals(4);
            let cut_c3 = hcons.parameter_intervals(5);
            Some(PrefAndRec::new(cut_c2, cut_c3, 5.0).map_err(|_| "AdvApprox_PrefAndRec")?)
        };
        let dicho_cut = DichoCutting;
        let cut: &dyn occt_core::adv_approx::Cutting = match &cut_pref {
            Some(p) => p,
            None => &dicho_cut,
        };
        let cont_shape = match a_cont {
            0 => Shape::C0,
            2 => Shape::C1,
            4 => Shape::C2,
            _ => return Err("Standard_ConstructionError"),
        };

        // `cxx:79-142` / `:174-227` / `:254-308`: the evaluator orders are the
        // 2D value/derivatives of the pcurve followed by the 3D ones of the
        // curve on surface (`Result[0..1]`, `Result[2..4]`).
        if only3d {
            let eval = |_a: f64, _b: f64, par: f64, order: i32, res: &mut [f64]| -> i32 {
                if res.len() != 3 {
                    return 1;
                }
                match order {
                    0 => {
                        let p = hcons.d0(par);
                        res[0] = p.x();
                        res[1] = p.y();
                        res[2] = p.z();
                    }
                    1 => {
                        let (_, d) = hcons.d1(par);
                        res[0] = d.x();
                        res[1] = d.y();
                        res[2] = d.z();
                    }
                    2 => {
                        let (_, _, d2) = hcons.d2(par);
                        res[0] = d2.x();
                        res[1] = d2.y();
                        res[2] = d2.z();
                    }
                    _ => {
                        res.fill(0.0);
                        return 3;
                    }
                }
                0
            };
            // `cxx:405-424`, `cxx:462-518`: Num1DSS = Num2DSS = 0, Num3DSS = 1.
            let a_approx = ApproxAFunction3d::approx(
                self.first,
                self.last,
                cont_shape,
                max_degree,
                max_segments,
                tol_3d,
                cut,
                &eval,
            )
            .map_err(|_| "AdvApprox_ApproxAFunction")?;
            self.is_done = a_approx.done;
            self.has_result = a_approx.has_result;
            if self.has_result {
                self.curve3d =
                    Some(GeomBSplineCurve::from_poles_knots_mults(
                        a_approx.poles.clone(),
                        a_approx.knots.clone(),
                        a_approx.mults.clone(),
                        a_approx.degree as usize,
                    )?);
                self.error3d = a_approx.max_error;
            }
            return Ok(());
        }

        let eval = |_a: f64, _b: f64, par: f64, order: i32, res: &mut [f64]| -> i32 {
            if res.len() != if only2d { 2 } else { 5 } {
                return 1;
            }
            match order {
                0 => {
                    let uv = trimmed.d0(par);
                    if only2d {
                        res[0] = uv.x();
                        res[1] = uv.y();
                    } else {
                        let p = hcons.d0(par);
                        res[0] = uv.x();
                        res[1] = uv.y();
                        res[2] = p.x();
                        res[3] = p.y();
                        res[4] = p.z();
                    }
                }
                1 => {
                    let (_, duv) = trimmed.d1(par);
                    if only2d {
                        res[0] = duv.x();
                        res[1] = duv.y();
                    } else {
                        let (_, d) = hcons.d1(par);
                        res[0] = duv.x();
                        res[1] = duv.y();
                        res[2] = d.x();
                        res[3] = d.y();
                        res[4] = d.z();
                    }
                }
                2 => {
                    let (_, _, d2uv) = trimmed.d2(par);
                    if only2d {
                        res[0] = d2uv.x();
                        res[1] = d2uv.y();
                    } else {
                        let (_, _, d2) = hcons.d2(par);
                        res[0] = d2uv.x();
                        res[1] = d2uv.y();
                        res[2] = d2.x();
                        res[3] = d2.y();
                        res[4] = d2.z();
                    }
                }
                _ => {
                    res.fill(0.0);
                    return 3;
                }
            }
            0
        };

        // `cxx:403-424` + `cxx:472-518`: Num1DSS = 2 and Num3DSS = 0 (`only2d`)
        // or 1.
        let a_approx = ApproxAFunction1dPair3d::approx(
            self.first,
            self.last,
            cont_shape,
            max_degree,
            max_segments,
            tol_u,
            tol_v,
            if only2d { None } else { Some(tol_3d) },
            cut,
            &eval,
        )
        .map_err(|_| "AdvApprox_ApproxAFunction")?;
        self.is_done = a_approx.done;
        self.has_result = a_approx.has_result;
        if self.has_result {
            // `cxx:528-551`: one knot vector, one degree, one pole count.
            let degree = a_approx.degree as usize;
            let flat = occt_core::bspl::banded_interp::knot_sequence(
                &a_approx.knots,
                &a_approx.mults,
                a_approx.degree,
            );
            self.curve2d = Some(Geom2dBSplineCurve::new(
                a_approx.poles_1d_u.clone(),
                a_approx.poles_1d_v.clone(),
                flat,
                degree,
            )?);
            self.error2du = a_approx.max_error_u;
            self.error2dv = a_approx.max_error_v;
            if let Some(poles3d) = a_approx.poles_3d.clone() {
                self.curve3d = Some(GeomBSplineCurve::from_poles_knots_mults(
                    poles3d,
                    a_approx.knots.clone(),
                    a_approx.mults.clone(),
                    degree,
                )?);
                self.error3d = a_approx.max_error_3d.unwrap_or(0.0);
            }
        }
        Ok(())
    }

    /// `isIsoLine` (`cxx:591-677`): a line with natural parameterisation,
    /// either `GeomAbs_Line`, a degree-1 two-pole B-spline, or a degree-1
    /// two-pole Bezier. `(isU, param, isForward)` where `isU == false` means
    /// horizontal (`V` constant).
    fn is_iso_line(c2d: &dyn Curve2d) -> Option<IsoLine> {
        let (loc, dir) = match type2d(c2d) {
            CurveType::Line => {
                let lin = c2d.gp_lin2d()?;
                (lin.location(), *lin.direction())
            }
            CurveType::BSplineCurve => {
                // `cxx:610-629`.
                if c2d.bspline_degree() != Some(1) {
                    return None;
                }
                let poles = c2d.poles2d()?;
                if poles.len() != 2 {
                    return None;
                }
                let v = GpVec2d::new(poles[1].x() - poles[0].x(), poles[1].y() - poles[0].y());
                if v.square_magnitude() < CONFUSION {
                    return None;
                }
                (poles[0], GpDir2d::from_vec2d(&v).ok()?)
            }
            CurveType::BezierCurve => {
                // `cxx:630-649`.
                if c2d.bezier_nb_poles() != Some(2) {
                    return None;
                }
                let poles = c2d.poles2d()?;
                if poles.len() != 2 {
                    return None;
                }
                let v = GpVec2d::new(poles[1].x() - poles[0].x(), poles[1].y() - poles[0].y());
                if v.square_magnitude() < CONFUSION {
                    return None;
                }
                (poles[0], GpDir2d::from_vec2d(&v).ok()?)
            }
            _ => return None,
        };

        // `cxx:656-673`: horizontal -> V = const, vertical -> U = const.
        let dx = GpDir2d { x: 1.0, y: 0.0 };
        let dy = GpDir2d { x: 0.0, y: 1.0 };
        if dir.is_parallel(&dx, ANGULAR) {
            Some((false, loc.y(), dir.dot(&dx) > 0.0))
        } else if dir.is_parallel(&dy, ANGULAR) {
            Some((true, loc.x(), dir.dot(&dy) > 0.0))
        } else {
            None
        }
    }

    /// `buildC3dOnIsoLine` (`cxx:682-811`).
    fn build_c3d_on_iso_line(
        &mut self,
        c2d: &dyn Curve2d,
        is_u: bool,
        param: f64,
        is_forward: bool,
    ) -> Result<bool, &'static str> {
        // `cxx:688-697`: `occ::down_cast<GeomAdaptor_Surface>` and the sphere
        // exclusion. Every surface of this port is a `Geom_Surface`, so the
        // down-cast succeeds; the type query is the sphere arm.
        if self.surf.gp_sphere().is_some() {
            return Ok(false);
        }

        // `cxx:699-707`.
        let a_surf = self.surf.clone();
        let (u1, u2) = a_surf.u_range();
        let (v1, v2) = a_surf.v_range();
        let f2d = c2d.d0(c2d.first_parameter());
        let l2d = c2d.d0(c2d.last_parameter());

        let mut is_to_trim = true;
        let mut work_surf = a_surf;
        let (p1, p2);
        if is_u {
            // `cxx:709-740`.
            let v1_param = f2d.y().min(l2d.y());
            let v2_param = f2d.y().max(l2d.y());
            if v2_param < v1 - self.tol || v1_param > v2 + self.tol {
                return Ok(false);
            } else if !v1.is_finite() || !v2.is_finite() {
                if (v2_param - v1_param).abs() < PCONFUSION {
                    return Ok(false);
                }
                work_surf = Arc::new(
                    GeomRectangularTrimmedSurface::uv(
                        work_surf.clone(),
                        u1,
                        u2,
                        v1_param,
                        v2_param,
                    )
                    .map_err(|_| "Geom_RectangularTrimmedSurface")?,
                );
                is_to_trim = false;
                p1 = v1_param;
                p2 = v2_param;
            } else {
                p1 = v1_param.max(v1);
                p2 = v2_param.min(v2);
                if (p2 - p1).abs() < PCONFUSION {
                    return Ok(false);
                }
            }
        } else {
            // `cxx:741-772`.
            let u1_param = f2d.x().min(l2d.x());
            let u2_param = f2d.x().max(l2d.x());
            if u2_param < u1 - self.tol || u1_param > u2 + self.tol {
                return Ok(false);
            } else if !u1.is_finite() || !u2.is_finite() {
                if (u2_param - u1_param).abs() < PCONFUSION {
                    return Ok(false);
                }
                work_surf = Arc::new(
                    GeomRectangularTrimmedSurface::uv(
                        work_surf.clone(),
                        u1_param,
                        u2_param,
                        v1,
                        v2,
                    )
                    .map_err(|_| "Geom_RectangularTrimmedSurface")?,
                );
                is_to_trim = false;
                p1 = u1_param;
                p2 = u2_param;
            } else {
                p1 = u1_param.max(u1);
                p2 = u2_param.min(u2);
                if (p2 - p1).abs() < PCONFUSION {
                    return Ok(false);
                }
            }
        }

        // `cxx:735` / `:768`: `aC3d = aSurf->UIso(theParam)` / `VIso`.
        let iso = if is_u {
            work_surf.u_iso_curve(param)
        } else {
            work_surf.v_iso_curve(param)
        };
        // `Geom_Surface::UIso/VIso` is pure virtual and returns a curve for
        // every translated surface type; `None` here means the type is not
        // translated, which OCCT cannot report.
        let iso = iso.ok_or("Geom_Surface::UIso/VIso UNPORTED")?;
        let a_c3d: Arc<dyn Curve> = if is_to_trim {
            Arc::new(GeomTrimmedCurveBasis::new(iso, p1, p2))
        } else {
            iso
        };

        // `cxx:775-779`: `GeomConvert::CurveToBSplineCurve(aC3d,
        // Convert_QuasiAngular)`, then `Reverse()` when the line runs backwards.
        let mut c3d = curve_to_bspline_curve(a_c3d.as_ref(), ParameterisationType::QuasiAngular)
            .map_err(|_| "GeomConvert::CurveToBSplineCurve")?;
        if !is_forward {
            c3d.reverse();
        }

        // `cxx:781-786`: `BSplCLib::Reparametrize(First, Last, aKnots)` +
        // `SetKnots`.
        let mut knots = c3d.distinct_knots().to_vec();
        reparametrize(c2d.first_parameter(), c2d.last_parameter(), &mut knots);
        c3d.set_knots(&knots)?;

        // `cxx:788-806`: 23-sample deviation between the 3d curve and the
        // surface point of the 2d curve.
        let mut error = 0.0_f64;
        const NB_PNT: i32 = 23;
        for idx in 0..=NB_PNT {
            let par = self.first + ((self.last - self.first) * idx as f64) / NB_PNT as f64;
            let p2 = c2d.d0(par);
            let p3 = c3d.d0(par);
            let ps = self.surf.d0(p2.x(), p2.y());
            error = error.max(p3.square_distance(&ps));
        }
        self.error3d = error.sqrt();
        self.curve3d = Some(c3d);
        // `cxx:808-810`: keep the general-purpose path when the tolerance is
        // not reached (isoline on a sphere and similar).
        Ok(self.error3d <= self.tol)
    }
}

/// `isIsoLine` of a trimmed 2D curve (`cxx:591-677`) is a free function as
/// well, so the caller does not need an installed approximation.
pub fn is_iso_line_of(c2d: &dyn Curve2d) -> Option<IsoLine> {
    ApproxCurveOnSurface::is_iso_line(c2d)
}
