//! `AppDef_Variational::ACR`, `SplitCurve` and the file-local `NearIndex` /
//! `GettingKnots` helpers.
//!
//! Source: `AppDef_Variational.cxx:1829-2072` (`ACR`, `NearIndex`,
//! `GettingKnots`) and `cxx:2074-2127` (`SplitCurve`) (OCCT 8.0.0).
//!
//! Translation notes:
//! - Both `Ti` and the curve knots are 1-based arrays in OCCT; the port keeps
//!   that layout. `NearIndex` therefore reports 1-based indices.
//! - `GettingKnots` fills `NewKnots(NbElm+1 .. NbMax)`; the caller passes a
//!   buffer whose upper bound is `myMaxSegment`.

use std::cell::RefCell;
use std::rc::Rc;

use occt_math::fem_tool::Curve;
use occt_math::vector::MathVector;

use super::Variational;

/// `NearIndex(T, TabPar, Eps, Flag)` (`cxx:1962-2005`), returning
/// `(index, flag)` with 1-based `index`.
pub(super) fn near_index(t: f64, tab_par: &[f64], eps: f64) -> (i32, i32) {
    let loi = 1i32;
    let upi = tab_par.len() as i32 - 1;
    let flag = 0;

    if t < tab_par[loi as usize] {
        return (loi, -1);
    }
    if t > tab_par[upi as usize] {
        return (upi, 1);
    }

    let mut ibeg = loi;
    let mut ifin = upi;

    while ibeg + 1 != ifin {
        let imidl = (ibeg + ifin) / 2;
        if t >= tab_par[ibeg as usize] && t <= tab_par[imidl as usize] {
            ifin = imidl;
        } else {
            ibeg = imidl;
        }
    }

    if (t - tab_par[ifin as usize]).abs() < eps {
        return (ifin, flag);
    }

    (ibeg, flag)
}

/// `GettingKnots(TabPar, InCurve, Deg, NbElm, NewKnots)` (`cxx:2011-2072`).
fn getting_knots(
    tab_par: &[f64],
    in_curve: &Curve,
    deg: i32,
    nb_elm: &mut i32,
    new_knots: &mut [f64],
) {
    const EPS: f64 = 1.0e-12;

    let old_knots = in_curve.knots();
    let nb_max_old = in_curve.nb_elements() as i32;
    let nb_max = new_knots.len() as i32 - 1;
    let mut el = 0i32;
    let mut i1 = 1usize;
    let mut i0 = 0usize;

    while *nb_elm < nb_max && el < nb_max_old {
        el += 1;
        i0 += 1;
        i1 += 1;

        if in_curve.degree(el as usize) == deg {
            *nb_elm += 1;

            let (mut ipt1, flag1) = near_index(old_knots[i0], tab_par, EPS);
            if flag1 != 0 {
                ipt1 = 1;
            }
            let (mut ipt2, flag2) = near_index(old_knots[i1], tab_par, EPS);
            if flag2 != 0 {
                ipt2 = tab_par.len() as i32 - 1;
            }

            if ipt2 - ipt1 >= 1 {
                let ipt = (ipt1 + ipt2) / 2;
                let tpar = if 2 * ipt == ipt1 + ipt2 {
                    2.0 * tab_par[ipt as usize]
                } else {
                    tab_par[ipt as usize] + tab_par[(ipt + 1) as usize]
                };
                new_knots[*nb_elm as usize] =
                    (old_knots[i0] + old_knots[i1] + tpar) / 4.0;
            } else {
                new_knots[*nb_elm as usize] = (old_knots[i0] + old_knots[i1]) / 2.0;
            }
        }
    }
}

/// `Curve->Length(First, Last, Len)` on a 1-based knot array.
fn curve_length(curve: &Rc<RefCell<Curve>>, knots: &[f64], first: f64, last: f64) -> f64 {
    let mut c = curve.borrow_mut();
    c.knots_mut().copy_from_slice(knots);
    c.length_of(first, last)
}

impl Variational {
    /// `ACR(Curve, Ti, Decima)` (`cxx:1829-1960`).
    pub(super) fn acr(
        &self,
        curve: &mut Rc<RefCell<Curve>>,
        ti: &mut MathVector,
        decima: i32,
    ) {
        const EPS: f64 = 1.0e-8;

        let nbr_pnt = ti.len() as i32;
        let ti_first = ti.lower() as i32;
        let ti_last = ti.upper() as i32;
        let kfirst = 1i32;
        let mut knots: Vec<f64> = curve.borrow().knots().to_vec();
        let klast = knots.len() as i32 - 1;
        let nb_cntr = self.nb_pass_points + self.nb_tang_points + self.nb_curv_points;

        // (1) Compute the curve length.
        let cb_long =
            curve_length(curve, &knots, ti.value(ti_first as usize), ti.value(ti_last as usize));

        // (2) Place the arc-length parameterization in Ti.
        if nbr_pnt >= 2 {
            // (2.0) Initialisation.
            let delta_t = (ti.value(ti_last as usize) - ti.value(ti_first as usize))
                / decima as f64;
            let mut vtest = ti.value(ti_first as usize) + delta_t;

            let mut icnt = 0i32;
            let mut pcnt;
            if nb_cntr > 0 {
                pcnt = self.typ_constraints[1] - self.first_point + ti_first;
                icnt = 1;
            } else {
                pcnt = ti_last + 1;
            }

            let mut u_old = 0.0;
            let mut t_old = ti.value(ti_first as usize);
            let mut p_old = ti_first;

            let mut ielm = kfirst;
            let mut i_old = ielm;

            ti.set_value(ti_first as usize, 0.0);

            let mut ipnt = ti_first + 1;
            while ipnt <= ti_last {
                while icnt <= nb_cntr && pcnt < ipnt {
                    icnt += 1;
                    pcnt = self.typ_constraints[(2 * icnt - 1) as usize] - self.first_point
                        + ti_first;
                }

                let mut tpara = ti.value(ipnt as usize);

                if tpara >= vtest || pcnt == ipnt {
                    if ti.value(ti_last as usize) - tpara <= 1.0e-2 * delta_t {
                        ipnt = ti_last;
                        tpara = ti.value(ipnt as usize);
                    }
                    // (2.2), (2.3) Compute the curve length.
                    let mut u_new =
                        curve_length(curve, &knots, ti.value(ti_first as usize), tpara);
                    u_new /= cb_long;

                    while knots[(ielm + 1) as usize] < tpara && ielm < klast - 1 {
                        ielm += 1;
                    }

                    // (2.4) Update the splitting parameters.
                    let dt_inv = 1.0 / (tpara - t_old);
                    let du = u_new - u_old;

                    for ii in (i_old + 1)..=ielm {
                        let ratio = (knots[ii as usize] - t_old) * dt_inv;
                        knots[ii as usize] = u_old + ratio * du;
                    }

                    // (2.5) Update the point parameters.
                    for ii in (p_old + 1)..ipnt {
                        let ratio = (ti.value(ii as usize) - t_old) * dt_inv;
                        ti.set_value(ii as usize, u_old + ratio * du);
                    }

                    ti.set_value(ipnt as usize, u_new);

                    u_old = u_new;
                    i_old = ielm;
                    t_old = tpara;
                    p_old = ipnt;
                }

                // New parametric threshold for decimation.
                if tpara >= vtest {
                    vtest += ((tpara - vtest + EPS) / delta_t).ceil() * delta_t;
                    if vtest > 1.0 - EPS {
                        vtest = 1.0;
                    }
                }

                ipnt += 1;
            }
        }

        // Adjust the extreme values.
        ti.set_value(ti_first as usize, 0.0);
        ti.set_value(ti_last as usize, 1.0);
        let mut ii = ti_last - 1;
        while ti.value(ii as usize) > knots[klast as usize] {
            ti.set_value(ii as usize, 1.0);
            ii -= 1;
        }
        knots[kfirst as usize] = 0.0;
        knots[klast as usize] = 1.0;

        curve.borrow_mut().knots_mut().copy_from_slice(&knots);
    }

    /// `SplitCurve(InCurve, Ti, CurveTol, OutCurve, iscut)`
    /// (`cxx:2074-2127`). Returns `(OutCurve, iscut)`.
    pub(super) fn split_curve(
        &self,
        in_curve: &Rc<RefCell<Curve>>,
        ti: &MathVector,
        curve_tol: f64,
    ) -> (Rc<RefCell<Curve>>, bool) {
        let nb_elm_old = in_curve.borrow().nb_elements() as i32;

        if nb_elm_old >= self.max_segment {
            return (in_curve.clone(), false);
        }

        let mut nb_elm = nb_elm_old;
        let mut new_knots = vec![0.0f64; self.max_segment.max(1) as usize + 1];
        {
            let c = in_curve.borrow();
            let work_degree = c.base().work_degree();
            getting_knots(ti.as_slice(), &c, work_degree, &mut nb_elm, &mut new_knots);
            getting_knots(ti.as_slice(), &c, work_degree - 1, &mut nb_elm, &mut new_knots);
        }

        if nb_elm > nb_elm_old {
            let in_knots: Vec<f64> = in_curve.borrow().knots().to_vec();
            let mut out = {
                let c = in_curve.borrow();
                Curve::new(c.dimension() as i32, nb_elm, c.base(), curve_tol)
            };
            {
                let out_knots = out.knots_mut();
                let out_lower = 1i32;
                for i in 1..in_knots.len() {
                    out_knots[i] = in_knots[i];
                }
                for i in (nb_elm_old + 1)..=nb_elm {
                    out_knots[(i + out_lower) as usize] = new_knots[i as usize];
                }
                out_knots[1..=nb_elm as usize + 1]
                    .sort_by(|a, b| a.partial_cmp(b).unwrap());
            }
            (Rc::new(RefCell::new(out)), true)
        } else {
            (in_curve.clone(), false)
        }
    }
}
