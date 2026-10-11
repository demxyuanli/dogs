//! `Hermit::Solution` for `Geom2d_BSplineCurve`.
//! Source: `Hermit.cxx` (`HermiteCoeff` `:91-139`, `PolyTest` `:520-842`,
//! `InsertKnots` `:849-860`, `MovePoles` `:867-880`, `Solution` `:994-1098`).

use crate::bspline_curve::Geom2dBSplineCurve;
use occt_core::bspl::curve_dn::dn;
use occt_core::bspl::knots::reparametrize;
use occt_core::gp::{GpPnt, GpPnt2d};
use occt_core::precision::CONFUSION;

const DEFAULT_TOL_POLES: f64 = 0.000001;
const DEFAULT_TOL_KNOTS: f64 = 0.000001;

/// `HermiteCoeff` for `Geom2d_BSplineCurve` (`Hermit.cxx:91-139`).
fn hermite_coeff(bs: &Geom2dBSplineCurve) -> Result<[f64; 4], &'static str> {
    let Some(weights) = bs.weights.as_deref() else {
        return Err("Hermit::HermiteCoeff: no weights");
    };
    let (mut knots, mults) = bs.distinct_knots_and_mults();
    reparametrize(0.0, 1.0, &mut knots);
    let degree = bs.degree() as i32;
    let periodic = bs.periodic;
    let index0 = bs.first_u_knot_index();
    let index1 = bs.last_u_knot_index() - 1;
    let poles: Vec<GpPnt> = weights.iter().map(|w| GpPnt::new(*w, 0.0, 0.0)).collect();
    let denom0 = dn(
        0.0,
        0,
        index0,
        degree,
        periodic,
        &poles,
        None,
        &knots,
        Some(&mults),
    )
    .x();
    let deriv0 = dn(
        0.0,
        1,
        index0,
        degree,
        periodic,
        &poles,
        None,
        &knots,
        Some(&mults),
    )
    .x();
    let denom1 = dn(
        1.0,
        0,
        index1,
        degree,
        periodic,
        &poles,
        None,
        &knots,
        Some(&mults),
    )
    .x();
    let deriv1 = dn(
        1.0,
        1,
        index1,
        degree,
        periodic,
        &poles,
        None,
        &knots,
        Some(&mults),
    )
    .x();
    if denom0 == 0.0 || denom1 == 0.0 {
        return Err("Hermit::HermiteCoeff: zero denominator");
    }
    Ok([
        1.0 / denom0,
        -deriv0 / (denom0 * denom0),
        -deriv1 / (denom1 * denom1),
        1.0 / denom1,
    ])
}

/// `SignDenom` (`Hermit.cxx:146-153`).
fn sign_denom(poles: &[f64; 4]) -> bool {
    poles[0] >= 0.0
}

/// `Polemax` (`Hermit.cxx:160-185`).
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

/// `PolyTest` for `Geom2d_BSplineCurve` (`Hermit.cxx:520-842`).
fn poly_test(
    herm: &[f64; 4],
    bs: &Geom2dBSplineCurve,
    boucle: &mut i32,
    tol_poles: f64,
    ux: f64,
    uy: f64,
) -> Result<(f64, f64), &'static str> {
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

    let (bs_knots, _) = bs.distinct_knots_and_mults();
    let nb = bs_knots.len();
    let mut knots: Vec<f64> = if i1 == i2 {
        if i3 == i4 || i3 == 0 {
            bs_knots.clone()
        } else {
            let mut k = bs_knots.clone();
            k.push(uy);
            k
        }
    } else if i3 == i4 || i3 == 0 {
        let mut k = bs_knots.clone();
        k.push(ux);
        k
    } else {
        let mut k = bs_knots.clone();
        k.push(ux);
        k.push(uy);
        k
    };
    knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let mut polesinit = [
        herm[0],
        herm[0] + herm[1] / 3.0,
        herm[3] - herm[2] / 3.0,
        herm[3],
    ];

    if tol_poles != 0.0 {
        let (min, max) = polemax(&polesinit);
        let polemin = polesinit[min];
        let polemax = polesinit[max];
        if polemax >= (1.0 / tol_poles) * polemin
            || (polemin == 0.0 && polemax >= 1.0 / tol_poles)
        {
            if polesinit[0] >= (1.0 / tol_poles) * polesinit[3]
                || polesinit[0] <= tol_poles * polesinit[3]
            {
                return Err("Hermit Impossible Tolerance");
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
                                let scale =
                                    pole3 * 10f64.powf(-0.5 * tol_poles.log10() - a / 2.0);
                                for p in &mut polesinit {
                                    *p -= scale;
                                }
                            } else if *boucle == 1 {
                                let scale =
                                    pole0 * 10f64.powf(a / 2.0 + 0.5 * tol_poles.log10());
                                for p in &mut polesinit {
                                    *p -= scale;
                                }
                                dercas = 1;
                            }
                        }
                        if pole0 > pole3 {
                            let a = (pole0 / pole3).log10();
                            if *boucle == 2 {
                                let scale =
                                    pole0 * 10f64.powf(-0.5 * tol_poles.log10() - a / 2.0);
                                for p in &mut polesinit {
                                    *p -= scale;
                                }
                            } else if *boucle == 1 {
                                let scale =
                                    pole3 * 10f64.powf(a / 2.0 + 0.5 * tol_poles.log10());
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
    }

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
        u4 = if i1 < 2 {
            us1
        } else {
            knots[(i1 - 1) as usize]
        };
    }

    if polesinit[1] >= 0.0 && polesinit[2] < 0.0 {
        let mut us2 = polesinit[2] / (polesinit[2] - polesinit[3]);
        if *boucle == 2 {
            // `knots(Length - 1)` is 1-based second-to-last (`Hermit.cxx:758`).
            us2 = knots[knots.len() - 2] + us2 * (1.0 - knots[knots.len() - 2]);
        }
        if *boucle == 1 && ux != 0.0 {
            us2 = uy + us2 * (1.0 - uy);
        }
        let (li, nu) = locate(us2, &knots);
        i1 = li;
        let us2 = nu;
        u5 = if i1 >= (knots.len() as i32 - 1) {
            us2
        } else {
            knots[i1 as usize]
        };
    }

    if dercas == 1 {
        *boucle += 1;
    }

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
            u4 = if i1 >= 2 {
                knots[(i1 - 1) as usize]
            } else {
                us1
            };
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
fn insert_knots_hermit(bs: &mut Geom2dBSplineCurve, u4: f64, u5: f64) -> Result<(), &'static str> {
    if u4 != 0.0 {
        bs.insert_knot(u4, 1, 0.0)?;
    }
    if u5 != 1.0 && u5 != u4 {
        bs.insert_knot(u5, 1, 0.0)?;
    }
    Ok(())
}

/// `MovePoles` (`Hermit.cxx:867-880`).
fn move_poles(bs: &mut Geom2dBSplineCurve) {
    let y0 = bs.ys[0];
    let n = bs.nb_poles();
    if n <= 4 {
        return;
    }
    for i in 2..(n - 2) {
        bs.ys[i] = y0;
    }
}

fn make_hermite_bs(herm: &[f64; 4]) -> Result<Geom2dBSplineCurve, &'static str> {
    let xs = vec![0.0, 0.0, 0.0, 0.0];
    let ys = vec![
        herm[0],
        herm[0] + herm[1] / 3.0,
        herm[3] - herm[2] / 3.0,
        herm[3],
    ];
    let flat = vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
    Geom2dBSplineCurve::new(xs, ys, flat, 3)
}

/// `Hermit::Solution` for `Geom2d_BSplineCurve` (`Hermit.cxx:994-1098`).
pub fn solution(
    bs: &Geom2dBSplineCurve,
    tol_poles: f64,
    tol_knots: f64,
) -> Result<Geom2dBSplineCurve, &'static str> {
    let mut herm = hermite_coeff(bs)?;
    let mut bs1 = make_hermite_bs(&herm)?;
    let mut bs2 = make_hermite_bs(&herm)?;
    let mut zeroboucle = 0i32;
    let (upos1, upos2) = poly_test(&herm, bs, &mut zeroboucle, CONFUSION, 1.0, 0.0)?;
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

    let n2 = bs2.nb_poles();
    herm[0] = bs2.ys[0];
    herm[1] = 3.0 * (bs2.ys[1] - bs2.ys[0]);
    herm[2] = 3.0 * (bs2.ys[n2 - 1] - bs2.ys[n2 - 2]);
    herm[3] = bs2.ys[n2 - 1];

    let mut boucle = 1i32;
    let (mut utol1, mut utol2) = poly_test(&herm, bs, &mut boucle, tol_poles, ux, uy)?;
    insert_knots_hermit(&mut bs2, utol1, utol2)?;

    let mut uint1 = 0.0;
    let mut uint2 = 1.0;
    if boucle == 2 {
        let n2 = bs2.nb_poles();
        herm[0] = bs2.ys[0];
        herm[1] = 3.0 * (bs2.ys[1] - bs2.ys[0]);
        herm[2] = 3.0 * (bs2.ys[n2 - 1] - bs2.ys[n2 - 2]);
        herm[3] = bs2.ys[n2 - 1];
        if utol1 == 0.0 {
            uint2 = utol2;
            let (a, b) = poly_test(&herm, bs, &mut boucle, tol_poles, uint2, 0.0)?;
            utol1 = a;
            utol2 = b;
        } else {
            uint1 = utol1;
            let (a, b) = poly_test(&herm, bs, &mut boucle, tol_poles, uint1, 0.0)?;
            utol1 = a;
            utol2 = b;
        }
        insert_knots_hermit(&mut bs2, utol1, utol2)?;
    }

    let (k2, _) = bs2.distinct_knots_and_mults();
    if k2.len() < 2 {
        return Err("Hermit Impossible Tolerance");
    }
    if k2[1] < tol_knots || k2[k2.len() - 2] > (1.0 - tol_knots) {
        return Err("Hermit Impossible Tolerance");
    }

    if upos2 == 1.0 && utol2 == 1.0 && uint2 == 1.0 {
        insert_knots_hermit(&mut bs1, k2[1], 1.0)?;
    } else if upos1 == 0.0 && utol1 == 0.0 && uint1 == 0.0 {
        insert_knots_hermit(&mut bs1, k2[k2.len() - 2], 1.0)?;
    } else {
        insert_knots_hermit(&mut bs1, k2[k2.len() - 2], k2[1])?;
    }
    move_poles(&mut bs1);
    Ok(bs1)
}

/// `Hermit::Solution` with default tolerances (`Hermit.hxx:58-61`).
pub fn solution_default(bs: &Geom2dBSplineCurve) -> Result<Geom2dBSplineCurve, &'static str> {
    solution(bs, DEFAULT_TOL_POLES, DEFAULT_TOL_KNOTS)
}
