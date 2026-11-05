//! `Hermit::Solution` for a 3D `Geom_BSplineCurve` — the rational weight
//! function of a B-spline approximated by a Hermite cubic, returned as the 2D
//! law curve the caller multiplies into the homogeneous poles.
//!
//! Source: `ModelingData/TKGeomBase/Hermit/Hermit.cxx` (OCCT 8.0.0):
//! * `HermiteCoeff(Geom_BSplineCurve)` `:36-89`;
//! * `SignDenom` `:146-158`, `Polemax` `:160-190`;
//! * `PolyTest(Geom_BSplineCurve)` `:192-518`;
//! * `InsertKnots` `:849-865`, `MovePoles` `:867-882`;
//! * `Solution(Geom_BSplineCurve)` `:884-987`.
//!
//! `Hermit.cxx` carries two instantiations of the same code — one for
//! `Geom_BSplineCurve` and one for `Geom2d_BSplineCurve` (`:91-143`,
//! `:520-846`, `:994-1100`). This module is the 3D instantiation; the 2D one is
//! `occt-geom2d::hermit`. `PolyTest` only ever reads `LocateU`, `NbKnots` and
//! `Knot(i)` from the curve (`cxx:203-263`), which is why the body ports
//! unchanged.
//!
//! The returned law curve is 2D: `GeomConvert::MultNumandDenom`
//! (`GeomConvert.cxx:496-583`) samples its `Y` coordinate only
//! (`GeomConvert_law_evaluator` `:463-492`).
//!
//! Divergences:
//! * `throw Standard_DimensionError("Hermit Impossible Tolerance")` becomes
//!   `Err`, the failure channel of this crate.
//! * `Geom2d_BSplineCurve(Poles, Knots, Mults, 3)` with four poles at
//!   `x = 0` is built through [`Geom2dBSplineCurve::new`] with the flat knot
//!   sequence `{0, 0, 0, 0, 1, 1, 1, 1}`.

use crate::bspline_curve::GeomBSplineCurve;
use occt_core::bspl::curve_dn::curve_dn;
use occt_core::bspl::knots::reparametrize;
use occt_core::gp::GpPnt;
use occt_core::precision::CONFUSION;
use occt_geom2d::bspline_curve::Geom2dBSplineCurve;

const DEFAULT_TOL_POLES: f64 = 0.000001;
const DEFAULT_TOL_KNOTS: f64 = 0.000001;

/// `HermiteCoeff` for `Geom_BSplineCurve` (`Hermit.cxx:36-84`): the weight
/// function `w(u)` and its derivative at `u = 0` and `u = 1`, as the four
/// Hermite coefficients `{1/w0, -w'0/w0^2, -w'1/w1^2, 1/w1}`.
fn hermite_coeff(bs: &GeomBSplineCurve) -> Result<[f64; 4], String> {
    let Some(weights) = bs.weights.as_deref() else {
        return Err("Hermit::HermiteCoeff: no weights".into());
    };
    let (mut knots, mults) = bs.distinct_knots_and_mults();
    reparametrize(0.0, 1.0, &mut knots);
    let degree = bs.degree() as i32;
    let periodic = bs.periodic();
    let index0 = bs.first_u_knot_index();
    let index1 = bs.last_u_knot_index() - 1;
    // The weight function of a rational B-spline is itself a polynomial
    // B-spline of the same degree with the weights as scalar poles.
    let poles: Vec<GpPnt> = weights.iter().map(|w| GpPnt::new(*w, 0.0, 0.0)).collect();
    let denom0 = curve_dn(
        0.0, 0, index0, degree, periodic, &poles, None, &knots, Some(&mults),
    )
    .x();
    let deriv0 = curve_dn(
        0.0, 1, index0, degree, periodic, &poles, None, &knots, Some(&mults),
    )
    .x();
    let denom1 = curve_dn(
        1.0, 0, index1, degree, periodic, &poles, None, &knots, Some(&mults),
    )
    .x();
    let deriv1 = curve_dn(
        1.0, 1, index1, degree, periodic, &poles, None, &knots, Some(&mults),
    )
    .x();
    if denom0 == 0.0 || denom1 == 0.0 {
        return Err("Hermit::HermiteCoeff: zero denominator".into());
    }
    Ok([
        1.0 / denom0,
        -deriv0 / (denom0 * denom0),
        -deriv1 / (denom1 * denom1),
        1.0 / denom1,
    ])
}

/// `SignDenom` (`Hermit.cxx:146-153`): `true` when `Herm(0) >= 0`.
fn sign_denom(poles: &[f64; 4]) -> bool {
    poles[0] >= 0.0
}

/// `Polemax` (`Hermit.cxx:160-185`): indices of the smallest and largest pole.
fn polemax(poles: &[f64; 4]) -> (usize, usize) {
    let mut min = 0usize;
    let mut max = 0usize;
    let mut ymin = poles[0];
    let mut ymax = poles[0];
    for (i, &y) in poles.iter().enumerate().skip(1) {
        if y < ymin {
            ymin = y;
            min = i;
        }
        if y > ymax {
            ymax = y;
            max = i;
        }
    }
    (min, max)
}

/// `PolyTest` for `Geom_BSplineCurve` (`Hermit.cxx:192-511`): the knots `U4`
/// (first half) and `U5` (second half) to insert into the Hermite polynomial so
/// that the four control values become positive. `TolKnots` is the unnamed
/// second tolerance parameter in OCCT, so the arity of the original is kept
/// and only `TolPoles` is read.
#[allow(clippy::too_many_arguments)]
fn poly_test(
    herm: &[f64; 4],
    bs: &GeomBSplineCurve,
    boucle: &mut i32,
    tol_poles: f64,
    _tol_knots: f64,
    ux: f64,
    uy: f64,
) -> Result<(f64, f64), String> {
    let mut u4 = 0.0;
    let mut u5 = 1.0;
    let mut i1 = 0i32;
    let mut i2 = 0i32;
    let mut i3 = 0i32;
    let mut i4 = 0i32;
    let mut cas = 0i32;
    let mut mark = 0i32;
    let mut dercas = 0i32;

    if ux != 1.0 {
        let (a, b) = bs.locate_u(ux, 0.0);
        i1 = a;
        i2 = b;
        if uy != 0.0 {
            let (c, d) = bs.locate_u(uy, 0.0);
            i3 = c;
            i4 = d;
        }
    }

    // `cxx:224-267`: the distinct knots plus the knots to be inserted.
    let mut knots: Vec<f64> = bs.distinct_knots().to_vec();
    if i1 == i2 {
        if !(i3 == i4 || i3 == 0) {
            knots.push(uy);
        }
    } else if i3 == i4 || i3 == 0 {
        knots.push(ux);
    } else {
        knots.push(ux);
        knots.push(uy);
    }
    knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let nb = knots.len();

    let mut polesinit = [
        herm[0],
        herm[0] + herm[1] / 3.0,
        herm[3] - herm[2] / 3.0,
        herm[3],
    ];

    // `cxx:280-388`: tolerance loop on the poles.
    if tol_poles != 0.0 {
        let (min, max) = polemax(&polesinit);
        let polemin = polesinit[min];
        let polemax = polesinit[max];
        if polemax >= (1.0 / tol_poles) * polemin || (polemin == 0.0 && polemax >= 1.0 / tol_poles)
        {
            if polesinit[0] >= (1.0 / tol_poles) * polesinit[3]
                || polesinit[0] <= tol_poles * polesinit[3]
            {
                return Err("Hermit Impossible Tolerance".into());
            }
            if max == 0 || max == 3 {
                for p in &mut polesinit {
                    *p -= tol_poles * polemax;
                }
            } else if max == 1 || max == 2 {
                if min == 0 || min == 3 {
                    for p in &mut polesinit {
                        *p -= (1.0 / tol_poles) * polemin;
                    }
                } else {
                    if tol_poles * polemax < polesinit[0] && tol_poles * polemax < polesinit[3] {
                        for p in &mut polesinit {
                            *p -= tol_poles * polemax;
                        }
                        mark = 1;
                    }
                    if (1.0 / tol_poles) * polemin > polesinit[0]
                        && (1.0 / tol_poles) * polemin > polesinit[3]
                        && mark == 0
                    {
                        for p in &mut polesinit {
                            *p -= (1.0 / tol_poles) * polemin;
                        }
                        mark = 1;
                    }
                    if mark == 0 {
                        let pole0 = polesinit[0];
                        let pole3 = polesinit[3];
                        if pole0 < pole3 {
                            let a = (pole3 / pole0).log10();
                            if *boucle == 2 {
                                let scale = pole3 * 10f64.powf(-0.5 * tol_poles.log10() - a / 2.0);
                                for p in &mut polesinit {
                                    *p -= scale;
                                }
                            } else if *boucle == 1 {
                                let scale = pole0 * 10f64.powf(a / 2.0 + 0.5 * tol_poles.log10());
                                for p in &mut polesinit {
                                    *p -= scale;
                                }
                                dercas = 1;
                            }
                        }
                        if pole0 > pole3 {
                            let a = (pole0 / pole3).log10();
                            if *boucle == 2 {
                                let scale = pole0 * 10f64.powf(-0.5 * tol_poles.log10() - a / 2.0);
                                for p in &mut polesinit {
                                    *p -= scale;
                                }
                            } else if *boucle == 1 {
                                let scale = pole3 * 10f64.powf(a / 2.0 + 0.5 * tol_poles.log10());
                                for p in &mut polesinit {
                                    *p -= scale;
                                }
                                dercas = 1;
                            }
                        }
                    }
                }
            }
        }
    } // end of the tolerance loop

    // `cxx:390-396`: sign of the polynomial.
    if !sign_denom(&polesinit) {
        for p in &mut polesinit {
            *p = -*p;
        }
    }

    let locate = |u: f64, knots: &[f64]| -> (i32, f64) {
        occt_core::bspl::locate::locate_parameter_range(
            knots,
            u,
            false,
            1,
            knots.len() as i32,
            0.0,
            1.0,
        )
    };

    // `cxx:398-422`: positivity loop, one negative interior pole.
    if polesinit[1] < 0.0 && polesinit[2] >= 0.0 {
        let mut us1 = polesinit[0] / (polesinit[0] - polesinit[1]);
        if *boucle == 2 {
            us1 *= knots[1];
        }
        if *boucle == 1 && ux != 0.0 {
            us1 *= ux;
        }
        let (li, nu) = locate(us1, &knots);
        i1 = li;
        let us1 = nu;
        u4 = if i1 < 2 { us1 } else { knots[(i1 - 1) as usize] };
    }

    // `cxx:424-447`.
    if polesinit[1] >= 0.0 && polesinit[2] < 0.0 {
        let mut us2 = polesinit[2] / (polesinit[2] - polesinit[3]);
        if *boucle == 2 {
            // `knots(knots.Length() - 1)` is 1-based second-to-last.
            us2 = knots[nb - 2] + us2 * (1.0 - knots[nb - 2]);
        }
        if *boucle == 1 && ux != 0.0 {
            us2 = uy + us2 * (1.0 - uy);
        }
        let (li, nu) = locate(us2, &knots);
        i1 = li;
        let us2 = nu;
        u5 = if i1 >= (nb as i32 - 1) {
            us2
        } else {
            knots[i1 as usize]
        };
    }

    if dercas == 1 {
        *boucle += 1;
    }

    // `cxx:452-511`: both interior poles negative.
    if polesinit[1] < 0.0 && polesinit[2] < 0.0 {
        let mut us1 = polesinit[0] / (polesinit[0] - polesinit[1]);
        let mut us2 = polesinit[2] / (polesinit[2] - polesinit[3]);
        if *boucle != 0 && ux != 0.0 {
            us1 *= ux;
            us2 = uy + us2 * (1.0 - uy);
        }
        if us2 <= us1 {
            let (li, nu) = locate(us1, &knots);
            i1 = li;
            us1 = nu;
            if knots[(i1 - 1) as usize] >= us2 {
                u4 = knots[(i1 - 1) as usize];
            } else {
                if i1 >= 2 {
                    u4 = knots[(i1 - 1) as usize];
                    let (li3, nu2) = locate(us2, &knots);
                    i3 = li3;
                    us2 = nu2;
                    if i3 < (nb as i32 - 1) {
                        u5 = knots[i3 as usize];
                        cas = 1;
                    }
                }
                if cas == 0 {
                    u4 = (us1 + us2) / 2.0;
                }
            }
        } else {
            let (li, nu) = locate(us1, &knots);
            i1 = li;
            us1 = nu;
            u4 = if i1 >= 2 { knots[(i1 - 1) as usize] } else { us1 };
            let (li3, nu2) = locate(us2, &knots);
            i3 = li3;
            us2 = nu2;
            u5 = if i3 < (nb as i32 - 1) {
                knots[i3 as usize]
            } else {
                us2
            };
        }
    }

    Ok((u4, u5))
}

/// `InsertKnots` (`Hermit.cxx:849-860`).
fn insert_knots_hermit(bs: &mut Geom2dBSplineCurve, u4: f64, u5: f64) -> Result<(), String> {
    if u4 != 0.0 {
        bs.insert_knot(u4, 1, 0.0)
            .map_err(|e| format!("Hermit::InsertKnots: {e}"))?;
    }
    if u5 != 1.0 && u5 != u4 {
        bs.insert_knot(u5, 1, 0.0)
            .map_err(|e| format!("Hermit::InsertKnots: {e}"))?;
    }
    Ok(())
}

/// `MovePoles` (`Hermit.cxx:867-880`): collapse the interior poles onto `y(0)`.
fn move_poles(bs: &mut Geom2dBSplineCurve) {
    let y0 = bs.ys[0];
    let n = bs.ys.len();
    if n <= 4 {
        return;
    }
    for y in bs.ys.iter_mut().take(n - 2).skip(2) {
        *y = y0;
    }
}

/// The four-pole degree-3 Hermite polynomial in B-spline form
/// (`Hermit.cxx:901-911`).
fn make_hermite_bs(herm: &[f64; 4]) -> Result<Geom2dBSplineCurve, String> {
    let xs = vec![0.0; 4];
    let ys = vec![
        herm[0],
        herm[0] + herm[1] / 3.0,
        herm[3] - herm[2] / 3.0,
        herm[3],
    ];
    Geom2dBSplineCurve::new(xs, ys, hermite_flat_knots(), 3)
        .map_err(|e| format!("Hermit::Solution: {e}"))
}

/// `Knots = {0, 1}` with `Multiplicities = {4, 4}` (`Hermit.cxx:907-909`).
fn hermite_flat_knots() -> Vec<f64> {
    vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]
}

/// Result of the control flow shared by `Hermit::Solution` (`Hermit.cxx:884-987`)
/// and `Hermit::Solutionbis` (`Hermit.cxx:1102-1206`): the untouched Hermite
/// B-spline `BS1`, the knots to be inserted into it, and the `BS2` knots that
/// close both functions.
struct HermitPrepared {
    bs1: Geom2dBSplineCurve,
    knots: Vec<f64>,
    upos1: f64,
    upos2: f64,
    utol1: f64,
    utol2: f64,
    uint1: f64,
    uint2: f64,
}

/// `Hermit.cxx:902-968` / `Hermit.cxx:1112-1185`: everything up to (and
/// including) the knots-tolerance check, which the two public entries share
/// verbatim.
fn hermit_prepare(
    bs: &GeomBSplineCurve,
    tol_poles: f64,
    tol_knots: f64,
) -> Result<HermitPrepared, String> {
    let mut herm = hermite_coeff(bs)?;
    let bs1 = make_hermite_bs(&herm)?;
    let mut bs2 = make_hermite_bs(&herm)?;
    let mut zeroboucle = 0i32;
    let (upos1, upos2) = poly_test(
        &herm,
        bs,
        &mut zeroboucle,
        CONFUSION,
        CONFUSION,
        1.0,
        0.0,
    )?;
    insert_knots_hermit(&mut bs2, upos1, upos2)?;

    let (ux, uy) = if upos1 != 0.0 {
        if upos2 != 1.0 {
            (upos1.min(upos2), upos1.max(upos2))
        } else {
            (upos1, upos1)
        }
    } else {
        (upos2, upos2)
    };

    let n2 = bs2.ys.len();
    herm[0] = bs2.ys[0];
    herm[1] = 3.0 * (bs2.ys[1] - bs2.ys[0]);
    herm[2] = 3.0 * (bs2.ys[n2 - 1] - bs2.ys[n2 - 2]);
    herm[3] = bs2.ys[n2 - 1];

    let mut boucle = 1i32;
    let (mut utol1, mut utol2) = poly_test(&herm, bs, &mut boucle, tol_poles, tol_knots, ux, uy)?;
    insert_knots_hermit(&mut bs2, utol1, utol2)?;

    let mut uint1 = 0.0;
    let mut uint2 = 1.0;
    if boucle == 2 {
        let n2 = bs2.ys.len();
        herm[0] = bs2.ys[0];
        herm[1] = 3.0 * (bs2.ys[1] - bs2.ys[0]);
        herm[2] = 3.0 * (bs2.ys[n2 - 1] - bs2.ys[n2 - 2]);
        herm[3] = bs2.ys[n2 - 1];
        if utol1 == 0.0 {
            uint2 = utol2;
            let (a, b) = poly_test(&herm, bs, &mut boucle, tol_poles, tol_knots, uint2, 0.0)?;
            utol1 = a;
            utol2 = b;
        } else {
            uint1 = utol1;
            let (a, b) = poly_test(&herm, bs, &mut boucle, tol_poles, tol_knots, uint1, 0.0)?;
            utol1 = a;
            utol2 = b;
        }
        insert_knots_hermit(&mut bs2, utol1, utol2)?;
    }

    // `cxx:965-968`: `BS2->Knot(2)` / `BS2->Knot(BS2->NbKnots() - 1)` over the
    // stored (distinct) knots.
    let k2 = bs2.distinct_knots_and_mults().0;
    if k2.len() < 2 {
        return Err("Hermit Impossible Tolerance".into());
    }
    if k2[1] < tol_knots || k2[k2.len() - 2] > (1.0 - tol_knots) {
        return Err("Hermit Impossible Tolerance".into());
    }

    Ok(HermitPrepared {
        bs1,
        knots: k2,
        upos1,
        upos2,
        utol1,
        utol2,
        uint1,
        uint2,
    })
}

/// `Hermit::Solution(const occ::handle<Geom_BSplineCurve>&, TolPoles, TolKnots)`
/// (`Hermit.cxx:884-987`).
pub fn solution(
    bs: &GeomBSplineCurve,
    tol_poles: f64,
    tol_knots: f64,
) -> Result<Geom2dBSplineCurve, String> {
    let mut prep = hermit_prepare(bs, tol_poles, tol_knots)?;
    let klo = prep.knots[1];
    let khi = prep.knots[prep.knots.len() - 2];

    // `cxx:970-984`.
    if prep.upos2 == 1.0 && prep.utol2 == 1.0 && prep.uint2 == 1.0 {
        insert_knots_hermit(&mut prep.bs1, klo, 1.0)?;
    } else if prep.upos1 == 0.0 && prep.utol1 == 0.0 && prep.uint1 == 0.0 {
        insert_knots_hermit(&mut prep.bs1, khi, 1.0)?;
    } else {
        insert_knots_hermit(&mut prep.bs1, khi, klo)?;
    }
    move_poles(&mut prep.bs1);
    Ok(prep.bs1)
}

/// `Hermit::Solution` with the default tolerances (`Hermit.hxx:58-61`).
pub fn solution_default(bs: &GeomBSplineCurve) -> Result<Geom2dBSplineCurve, String> {
    solution(bs, DEFAULT_TOL_POLES, DEFAULT_TOL_KNOTS)
}

/// `Hermit::Solutionbis(const occ::handle<Geom_BSplineCurve>&, double& Knotmin,
/// double& Knotmax, TolPoles, TolKnots)` (`Hermit.cxx:1102-1206`).
///
/// `cxx:1193-1204` only assigns the bound that the test on the finally
/// inserted knots justifies, so both out-parameters keep the caller's value in
/// the other branch.
pub fn solutionbis(
    bs: &GeomBSplineCurve,
    knotmin: &mut f64,
    knotmax: &mut f64,
    tol_poles: f64,
    tol_knots: f64,
) -> Result<(), String> {
    let prep = hermit_prepare(bs, tol_poles, tol_knots)?;
    let klo = prep.knots[1];
    let khi = prep.knots[prep.knots.len() - 2];

    // `cxx:1187-1205`.
    if prep.upos2 == 1.0 && prep.utol2 == 1.0 && prep.uint2 == 1.0 {
        *knotmin = klo;
    } else if prep.upos1 == 0.0 && prep.utol1 == 0.0 && prep.uint1 == 0.0 {
        *knotmax = khi;
    } else {
        *knotmin = klo;
        *knotmax = khi;
    }
    Ok(())
}

/// `Hermit::Solutionbis` with the default tolerances (`Hermit.hxx:50-59`).
pub fn solutionbis_default(
    bs: &GeomBSplineCurve,
    knotmin: &mut f64,
    knotmax: &mut f64,
) -> Result<(), String> {
    solutionbis(bs, knotmin, knotmax, DEFAULT_TOL_POLES, DEFAULT_TOL_KNOTS)
}
