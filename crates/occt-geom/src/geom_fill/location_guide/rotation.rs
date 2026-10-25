//! `GeomFill_LocationGuide::SetRotation` and `InitX`.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_LocationGuide.cxx:128-467`
//! (`InGoodPeriod`, `SetRotation`) and `cxx:1352-1425` (`InitX`), OCCT 8.0.0.
//!
//! Translation notes:
//! - `SetRotation` is called from the constructor path only through `Set` and
//!   `SetCurve`; it fills `myPoles2d`, `Uf` / `Ul`, `Inf` / `Sup`, `R` and
//!   `myStatus`, and returns `LastAngle`.
//! - The `Extrema_ExtCS` step is [`curve_surface_extrema_all`]; OCCT scans every
//!   extremum and keeps the first square distance strictly below the running
//!   minimum (`cxx:375-389`), which is what the loop below does.
//! - The `SOS` branch (`cxx:333-371`) leaves `Pc` default-constructed, so
//!   `w = Pc.Parameter()` below is `0.0` and `theU` / `theV` come from
//!   `IntCurveSurface_IntersectionPoint::SetValues(P, RR(2), RR(3), RR(1))`,
//!   i.e. `RR(2)` and `RR(3)`. Kept verbatim.
//! - `a1` (`cxx:392-393`) is converted with `InGoodPeriod` and then never read
//!   (`Angle` is taken from `theU` at `cxx:420`); it is kept as `_a1` so the
//!   `Extrema` branch stays recognisable.
//! - `InitX` uses `myPoles2d->RowLength()`, which is **2** for the
//!   `(1, 2, 1, myNbPts)` array built at `cxx:174-175`; the binary search of
//!   `cxx:1369-1389` therefore never runs and the seed interpolates between the
//!   first two samples. That OCCT quirk is reproduced on purpose.

use std::sync::Arc;

use occt_core::elib::clib::in_period;
use occt_core::gp::{GpAx1, GpAx3, GpDir, GpMat, GpPnt, GpVec, GpXyz};
use occt_core::precision::CONFUSION;
use occt_math::{MathFunctionSetRoot, MathVector};

use super::LocationGuide;
use crate::bspline_curve::GeomBSplineCurve;
use crate::curve::Curve;
use crate::extrema_surf::curve_surface_extrema_all;
use crate::geom_fill::function_guide::FunctionGuide;
use crate::geom_fill::pipe_error::GeomFillPipeError;
use crate::surface_of_revolution::GeomSurfaceOfRevolution;
use crate::trimmed::GeomTrimmedCurveBasis;

/// `InGoodPeriod(Prec, Period, Current)` (`cxx:128-142`).
pub fn in_good_period(prec: f64, period: f64, current: &mut f64) {
    let diff = *current - prec;
    let nb = (diff / period).trunc() as i32;
    *current -= nb as f64 * period;
    let diff = *current - prec;
    if diff > period / 2.0 {
        *current -= period;
    } else if diff < -period / 2.0 {
        *current += period;
    }
}

/// `math_FunctionSetRoot` solve of the guide / revolution-surface system
/// (`cxx:339-360` and `cxx:617-645`).
///
/// Returns `(root, function_set_error_norm)`; `None` is `!IsDone()`. OCCT reads
/// `Result.FunctionSetErrors().Norm()` only in the `SOS` arm, the caller checks
/// the norm itself.
#[allow(clippy::too_many_arguments)]
pub fn solve_guide(
    sec: &Arc<dyn crate::geom_fill::section_law::SectionLaw>,
    guide: &Arc<dyn Curve>,
    u_on_s: f64,
    p: &GpPnt,
    dir: &GpXyz,
    dx: &GpXyz,
    seed: [f64; 3],
    inf: [f64; 3],
    sup: [f64; 3],
    tol: [f64; 3],
    iter: i32,
) -> Option<([f64; 3], f64)> {
    let mut e = FunctionGuide::new(Arc::clone(sec), Arc::clone(guide), u_on_s);
    e.set_param(p, dir, dx).ok()?;

    let mut x = MathVector::new(1, 3);
    let mut lo = MathVector::new(1, 3);
    let mut hi = MathVector::new(1, 3);
    let mut tolerance = MathVector::new(1, 3);
    for i in 1..=3 {
        x.set_value(i, seed[i - 1]);
        lo.set_value(i, inf[i - 1]);
        hi.set_value(i, sup[i - 1]);
        tolerance.set_value(i, tol[i - 1]);
    }
    let mut root = MathFunctionSetRoot::with_iterations(&e, iter);
    root.set_tolerance(&tolerance);
    root.perform_with_bounds(&mut e, &x, &lo, &hi, false);
    if !root.is_done() {
        return None;
    }
    let sol = root.root();
    let norm = root.function_set_errors().norm();
    let mut out = [0.0; 3];
    for i in 1..=3 {
        out[i - 1] = sol.value(i);
    }
    Some((out, norm))
}

impl LocationGuide {
    /// `GeomFill_LocationGuide::SetRotation` (`cxx:184-467`). Returns
    /// `LastAngle`; an abandoned contact sets `myStatus` and keeps the angle.
    pub(crate) fn set_rotation(&self, prec_angle: f64) -> Result<f64, String> {
        let mut st = self.lock();
        let Some(curve) = st.curve.clone() else {
            return Err(
                "Standard_ConstructionError: GeomFill_LocationGuide::The path is not set !!".into(),
            );
        };
        let Some(sec) = st.sec.clone() else {
            return Err("Standard_NullObject: GeomFill_LocationGuide::SetRotation".into());
        };
        let Some(guide) = st.guide.clone() else {
            return Err("Standard_NullObject: GeomFill_LocationGuide::SetRotation".into());
        };
        let law = st.law.clone();

        // Fixed frame (`cxx:191`).
        let rep = GpAx3::new(
            GpPnt::zero(),
            GpDir::new(0.0, 0.0, 1.0)?,
            &GpDir::new(1.0, 0.0, 0.0)?,
        )?;

        let f = curve.first_parameter();
        let l = curve.last_parameter();
        let uperiodic = sec.is_u_periodic();
        let delta_g = (guide.last_parameter() - guide.first_parameter()) / 5.0;

        // Section shape and range (`cxx:217-241`).
        let mut tol = CONFUSION;
        let (isconst, israt, uf, ul, mults, knots, deg) = if sec.is_constant(&mut tol) {
            let my_section = sec.constant_section()?;
            let uf = my_section.first_parameter();
            let ul = my_section.last_parameter();
            (true, false, uf, ul, Vec::new(), Vec::new(), 0)
        } else {
            let (nb_poles, nb_knots, deg) = sec.section_shape();
            let _ = nb_poles;
            let knots = sec.knots();
            if knots.len() != nb_knots as usize || knots.is_empty() {
                return Err("GeomFill_LocationGuide::SetRotation: inconsistent section law".into());
            }
            let mults = sec.mults();
            let uf = knots[0];
            let ul = knots[knots.len() - 1];
            let israt = sec.is_rational();
            (false, israt, uf, ul, mults, knots, deg)
        };
        st.uf = uf;
        st.ul = ul;

        // Computation bounds (`cxx:243-257`).
        let delta = guide.last_parameter() - guide.first_parameter();
        st.inf[0] = guide.first_parameter() - delta / 10.0;
        st.sup[0] = guide.last_parameter() + delta / 10.0;
        st.inf[1] = -std::f64::consts::PI;
        st.sup[1] = 3.0 * std::f64::consts::PI;
        let delta = ul - uf;
        st.inf[2] = uf - delta / 10.0;
        st.sup[2] = ul + delta / 10.0;
        let uperiod = if uperiodic { ul - uf } else { 0.0 };

        let nb_pts = st.nb_pts;
        let mut old_angle = 0.0;
        let mut cur_angle = prec_angle;

        for ii in 1..=nb_pts {
            let t = ((nb_pts - ii) as f64 * f + (ii - 1) as f64 * l) / (nb_pts - 1) as f64;
            let p = curve.d0(t);
            let Ok((tangent, normal, binormal)) = law.d0(t) else {
                st.status = law.error_status();
                return Ok(cur_angle);
            };

            // `D = T`, transformed by `Trans` when `WithTrans` (`cxx:274-281`).
            let mut d = tangent;
            if st.with_trans {
                let m = GpMat::from_cols(normal.xyz(), binormal.xyz(), tangent.xyz())
                    .multiply(&st.trans);
                d = GpVec::from_xyz(&m.column(3));
            }
            let axe = GpAx1::new(p, GpDir::from_xyz(&d.xyz())?);

            let n3 = GpAx3::new(p, GpDir::from_xyz(&d.xyz())?, &GpDir::from_xyz(&normal.xyz())?)?;

            // `Transfo.SetTransformation(N3, Rep)` (`cxx:288`).
            let mut transfo = occt_core::gp::GpTrsf::identity();
            transfo.set_transformation_from_to(&n3, &rep);

            // Section at `t` (`cxx:290-316`).
            let section_curve: Box<dyn Curve> = if !isconst {
                let u = st.first_s + (t - curve.first_parameter()) * st.ratio;
                let (poles, weights) = sec.d0(u)?;
                let flat = occt_core::bspl::banded_interp::knot_sequence(&knots, &mults, deg);
                if israt {
                    Box::new(GeomBSplineCurve::from_parts(
                        poles,
                        Some(weights),
                        flat,
                        deg.max(0) as usize,
                        sec.is_u_periodic(),
                    ))
                } else {
                    Box::new(GeomBSplineCurve::from_parts(
                        poles,
                        None,
                        flat,
                        deg.max(0) as usize,
                        sec.is_u_periodic(),
                    ))
                }
            } else {
                let base = sec.constant_section()?;
                base.clone_dyn()
            };
            let mut s: Box<dyn Curve> = Box::new(GeomTrimmedCurveBasis::new(
                Arc::from(section_curve),
                uf,
                ul,
            ));
            s.transform(&transfo);
            let s: Arc<dyn Curve> = Arc::from(s);
            let revol = GeomSurfaceOfRevolution::new(s, axe);

            // `Extrema_ExtCS(*myGuide, GArevol, Confusion, Confusion)`
            // (`cxx:324`), then the minimum square distance scan.
            let extrema = curve_surface_extrema_all(guide.as_ref(), &revol);
            let best = extrema
                .iter()
                .enumerate()
                .fold(None::<(usize, f64)>, |acc, (idx, e)| {
                    let d2 = e.distance * e.distance;
                    match acc {
                        Some((_, best_d2)) if d2 >= best_d2 => acc,
                        _ => Some((idx, d2)),
                    }
                });

            let (the_u, the_v, w);
            if best.is_none() {
                // `SOS` branch (`cxx:333-371`): keep `Pc` default (parameter 0).
                if ii <= 1 {
                    st.status = GeomFillPipeError::ImpossibleContact;
                    return Ok(cur_angle);
                }
                let prev_w = st.poles2d[ii - 2][0].y();
                let prev_angle = st.poles2d[ii - 2][1].x();
                let prev_v = st.poles2d[ii - 2][1].y();
                let seed = [prev_w, prev_angle, prev_v];
                let u_on_s = st.first_s + (t - curve.first_parameter()) * st.ratio;
                let dir = d.xyz();
                let dx = normal.xyz();
                let solved = solve_guide(
                    &sec,
                    &guide,
                    u_on_s,
                    &p,
                    &dir,
                    &dx,
                    seed,
                    st.inf,
                    st.sup,
                    st.tol_res,
                    100,
                );
                let Some((root, norm)) = solved else {
                    st.status = GeomFillPipeError::ImpossibleContact;
                    return Ok(cur_angle);
                };
                if !(norm < st.tol_res[0] * st.tol_res[0]) {
                    st.status = GeomFillPipeError::ImpossibleContact;
                    return Ok(cur_angle);
                }
                the_u = root[1];
                the_v = root[2];
                w = 0.0; // `Pc` stays default-constructed: `Parameter() == 0.0`
            } else {
                let (idx, _) = best.unwrap();
                let e = &extrema[idx];
                the_u = e.u2;
                the_v = e.v2.unwrap_or(0.0);
                w = e.u1;
                // `a1 = theU; InGoodPeriod(CurAngle, 2*M_PI, a1);` (`cxx:392-393`)
                // is dead: `Angle` below is taken from `theU`.
                let mut a1 = the_u;
                in_good_period(cur_angle, 2.0 * std::f64::consts::PI, &mut a1);
                let _ = a1;
            };

            // Guide parameter window (`cxx:396-417`).
            let mut w = w;
            if ii > 1 {
                let prev_w = st.poles2d[ii - 2][0].y();
                let mut diff = w - prev_w;
                if diff.abs() > delta_g && guide.is_periodic() {
                    in_good_period(prev_w, guide.period(), &mut w);
                    diff = w - prev_w;
                }
                let _ = diff;
            }

            // Angle window (`cxx:420-435`).
            let mut angle = the_u;
            if ii > 1 {
                let mut diff = angle - old_angle;
                if diff.abs() > std::f64::consts::PI {
                    in_good_period(old_angle, 2.0 * std::f64::consts::PI, &mut angle);
                    diff = angle - old_angle;
                }
                let _ = diff;
            }

            // Section parameter window (`cxx:438-454`).
            let mut v = the_v;
            if ii > 1 && uperiodic {
                let prev_v = st.poles2d[ii - 2][1].y();
                in_good_period(prev_v, uperiod, &mut v);
            }

            st.poles2d[ii - 1] = [
                occt_core::gp::GpPnt2d::new(t, w),
                occt_core::gp::GpPnt2d::new(angle, v),
            ];
            old_angle = angle;
            cur_angle = angle;
        }

        st.rotation = true;
        Ok(cur_angle)
    }

    /// `GeomFill_LocationGuide::InitX` (`cxx:1352-1425`).
    pub(crate) fn init_x(&self, param: f64) {
        let mut st = self.lock();
        let mut ideb = 1usize;
        // `myPoles2d->RowLength()` is 2 (`cxx:174-175`, `hxx:193-197`).
        let mut ifin = 2usize;

        let valeur = st.poles2d[ideb - 1][0].x();
        if param == valeur {
            ifin = ideb + 1;
        }
        let valeur = st.poles2d[ifin - 1][0].x();
        if param == valeur {
            ideb = ifin - 1;
        }
        while ideb + 1 != ifin {
            let idemi = (ideb + ifin) / 2;
            let valeur = st.poles2d[idemi - 1][0].x();
            if valeur < param {
                ideb = idemi;
            } else if valeur > param {
                ifin = idemi;
            } else {
                ideb = idemi;
                ifin = ideb + 1;
            }
        }

        let t1 = st.poles2d[ideb - 1][0].x();
        let t2 = st.poles2d[ifin - 1][0].x();
        let diff = t2 - t1;
        let w1 = st.poles2d[ideb - 1][0].y();
        let w2 = st.poles2d[ifin - 1][0].y();
        let p1 = st.poles2d[ideb - 1][1];
        let p2 = st.poles2d[ifin - 1][1];

        if diff > 1.0e-7 {
            let b = (param - t1) / diff;
            let a = (t2 - param) / diff;
            st.x[0] = a * w1 + b * w2;
            st.x[1] = a * p1.x() + b * p2.x();
            st.x[2] = a * p1.y() + b * p2.y();
        } else {
            st.x[0] = (w1 + w2) / 2.0;
            st.x[1] = (p1.x() + p2.x()) / 2.0;
            st.x[2] = (p1.y() + p2.y()) / 2.0;
        }

        if let Some(guide) = &st.guide {
            if guide.is_periodic() {
                st.x[0] = in_period(st.x[0], guide.first_parameter(), guide.last_parameter());
            }
        }
        st.x[1] = in_period(st.x[1], 0.0, 2.0 * std::f64::consts::PI);
        if let Some(sec) = &st.sec {
            if sec.is_u_periodic() {
                st.x[2] = in_period(st.x[2], st.uf, st.ul);
            }
        }
    }
}
