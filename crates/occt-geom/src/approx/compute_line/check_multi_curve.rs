//! `CheckMultiCurve` (`Approx_ComputeLine.gxx:134-428`): the loop detector
//! that rejects an approximation whose Bezier poles double back on
//! themselves.
//!
//! Translation notes:
//! - OCCT writes `indbads` as `int indbads[4]` (`gxx:151`) and only ever scans
//!   `indbads[1..3]` (`gxx:416-421`). The port sizes the vector
//!   `max(4, NbCur + 1)` so a `NbCur > 3` input does not index out of bounds;
//!   the scan itself is unchanged.
//! - `gp_Vec::Normalize` / `gp_Vec2d::Normalize` are
//!   [`GpVec::normalized`] / `GpVec2d::normalized`.
//! - `gp_Vec2d.Dot` has no `occt-core` counterpart, so `gp_XY::Dot` is
//!   written out.

use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_core::precision::{epsilon, SQUARE_CONFUSION, REAL_SMALL};

use crate::approx::line_tool::LineTool;
use crate::app_par_curves::multi_curve::MultiCurve;

/// `CheckMultiCurve` (`gxx:134-428`).
pub(crate) fn check_multi_curve<L: LineTool>(
    the_multi_curve: &MultiCurve,
    the_line: &L,
    the_indfirst: i32,
    the_indlast: i32,
    the_indbad: &mut i32,
) -> bool {
    let nbp3d = the_line.nb_p3d();
    // `nbp2d` (`gxx:141`) only sizes the temporary `tabP`/`tabP2d` arrays,
    // which the port reads straight from `MyLineTool`, so it has no read use.
    let _nbp2d = the_line.nb_p2d();

    // 2*2
    let coeff = 4.0;
    let min_scal_prod = -0.9;
    let sq_tol3d = SQUARE_CONFUSION;

    // only simple cases are analysed
    if nbp3d > 1 {
        return true;
    }

    *the_indbad = 0;
    let nb_cur = the_multi_curve.nb_curves() as i32;
    let mut indbads = vec![0i32; (nb_cur + 1).max(4) as usize];
    let mut loop_found = false;

    if the_multi_curve.dimension(1) == 3 {
        let a_poles = the_multi_curve.curve3d(1);
        let upper = a_poles.len() as i32;

        let mut first_vec = GpVec::zero();
        let mut indp = 2i32;
        while indp <= upper {
            first_vec = GpVec::from_pnts(&a_poles[0], &a_poles[(indp - 1) as usize]);
            indp += 1;
            let a_length = first_vec.magnitude();
            if a_length > REAL_SMALL {
                first_vec = first_vec.divided(a_length);
                break;
            }
        }
        let mut mid_pnt = a_poles[(indp - 1) as usize];
        while indp <= upper {
            let second_raw = GpVec::from_pnts(&mid_pnt, &a_poles[(indp - 1) as usize]);
            let a_length = second_raw.magnitude();
            if a_length <= REAL_SMALL {
                indp += 1;
                continue;
            }
            let second_vec = second_raw.divided(a_length);
            let scal_prod = first_vec.dot(&second_vec);
            if scal_prod < min_scal_prod {
                loop_found = true;
                break;
            }
            first_vec = second_vec;
            mid_pnt = a_poles[(indp - 1) as usize];
            indp += 1;
        }
        // Check: may be it is a real loop.
        if loop_found {
            'outer: for first_ind in the_indfirst..=the_indlast - 2 {
                let first_pnt = first_point3d(the_line, first_ind);
                for k in first_ind + 1..the_indlast {
                    let pnt1 = first_point3d(the_line, k);
                    let pnt2 = first_point3d(the_line, k + 1);
                    if first_pnt.square_distance(&pnt1) <= sq_tol3d
                        || first_pnt.square_distance(&pnt2) <= sq_tol3d
                    {
                        loop_found = false;
                        break 'outer;
                    }
                    let vec1 = GpVec::from_pnts(&first_pnt, &pnt1).normalized();
                    let vec2 = GpVec::from_pnts(&first_pnt, &pnt2).normalized();
                    let scal_prod = vec1.dot(&vec2);
                    if scal_prod < min_scal_prod {
                        loop_found = false;
                        break 'outer;
                    }
                }
            }
        }
        if loop_found {
            // search <indbad>
            let mut max_sq_dist = 0.0f64;
            let mut min_sq_dist = f64::MAX;
            for k in the_indfirst + 1..=the_indlast {
                let prev_pnt = first_point3d(the_line, k - 1);
                let cur_pnt = first_point3d(the_line, k);
                let a_sq_dist = prev_pnt.square_distance(&cur_pnt);
                if a_sq_dist > max_sq_dist {
                    max_sq_dist = a_sq_dist;
                    indbads[1] = k;
                }
                if a_sq_dist > REAL_SMALL && a_sq_dist < min_sq_dist {
                    min_sq_dist = a_sq_dist;
                }
            }
            let relation = max_sq_dist / min_sq_dist;
            if relation < coeff {
                loop_found = false;
            } else {
                for indcur in 2..=nb_cur {
                    let mut max_sq_dist = 0.0f64;
                    for k in the_indfirst + 1..=the_indlast {
                        let prev = the_line.value_2d(k - 1);
                        let cur = the_line.value_2d(k);
                        let idx = (indcur - 2) as usize;
                        let a_sq_dist = prev[idx].square_distance(&cur[idx]);
                        if a_sq_dist > max_sq_dist {
                            max_sq_dist = a_sq_dist;
                            indbads[indcur as usize] = k;
                        }
                    }
                }
            }
        }
    } else {
        // 2d case
        let a_poles2d = the_multi_curve.curve2d(1);
        let a_sq_norm_toler = epsilon(1.0) * epsilon(1.0);
        let mut first_vec = vec2d_from_pnts(&a_poles2d[0], &a_poles2d[1]);
        let mut a_vec_sq_norm = square_magnitude2d(&first_vec);
        if a_vec_sq_norm < a_sq_norm_toler {
            *the_indbad = the_indfirst + 1;
            return false;
        }

        // OCCT divides by `sqrt(aSqNormToler)` here (`gxx:319`), not by
        // `sqrt(aVecSqNorm)`; kept verbatim.
        first_vec = first_vec.divided(a_sq_norm_toler.sqrt());
        let mut mid_pnt = a_poles2d[1];
        for k in 3..=a_poles2d.len() as i32 {
            // `SecondVec.SetXY(aPoles2d(k).XY() - MidPnt.XY())` (`gxx:325`).
            let second_raw = vec2d_from_pnts(&mid_pnt, &a_poles2d[(k - 1) as usize]);
            a_vec_sq_norm = square_magnitude2d(&second_raw);
            if a_vec_sq_norm < a_sq_norm_toler {
                *the_indbad = the_indfirst + k - 1;
                return false;
            }
            let second_vec = second_raw.divided(a_vec_sq_norm.sqrt());
            let scal_prod = first_vec.x() * second_vec.x() + first_vec.y() * second_vec.y();
            if scal_prod < min_scal_prod {
                loop_found = true;
                break;
            }
            first_vec = second_vec;
            mid_pnt = a_poles2d[(k - 1) as usize];
        }
        // Check: may be it is a real loop.
        if loop_found {
            'outer: for first_ind in the_indfirst..=the_indlast - 2 {
                let first_pnt = first_point2d(the_line, first_ind);
                for k in first_ind + 1..the_indlast {
                    let pnt1 = first_point2d(the_line, k);
                    let pnt2 = first_point2d(the_line, k + 1);
                    if first_pnt.square_distance(&pnt1) <= sq_tol3d
                        || first_pnt.square_distance(&pnt2) <= sq_tol3d
                    {
                        loop_found = false;
                        break 'outer;
                    }
                    let vec1 = vec2d_from_pnts(&first_pnt, &pnt1)
                        .normalized()
                        .unwrap_or(GpVec2d::zero());
                    let vec2 = vec2d_from_pnts(&first_pnt, &pnt2)
                        .normalized()
                        .unwrap_or(GpVec2d::zero());
                    let scal_prod = vec1.x() * vec2.x() + vec1.y() * vec2.y();
                    if scal_prod < min_scal_prod {
                        loop_found = false;
                        break 'outer;
                    }
                }
            }
        }
        if loop_found {
            // search <indbad>
            for indcur in 1..=nb_cur {
                let mut max_sq_dist = 0.0f64;
                let mut min_sq_dist = f64::MAX;
                for k in the_indfirst + 1..=the_indlast {
                    let prev = the_line.value_2d(k - 1);
                    let cur = the_line.value_2d(k);
                    let idx = (indcur - 1) as usize;
                    let a_sq_dist = prev[idx].square_distance(&cur[idx]);
                    if a_sq_dist > max_sq_dist {
                        max_sq_dist = a_sq_dist;
                        indbads[indcur as usize] = k;
                    }
                    if a_sq_dist > REAL_SMALL && a_sq_dist < min_sq_dist {
                        min_sq_dist = a_sq_dist;
                    }
                }
                let relation = max_sq_dist / min_sq_dist;
                if relation < coeff {
                    loop_found = false;
                }
            }
        }
    }

    // Define <indbad>
    for i in 1..=3usize {
        if indbads[i] != 0 {
            *the_indbad = indbads[i];
            break;
        }
    }

    if !loop_found {
        *the_indbad = 0;
    }

    !loop_found
}

/// `LineTool::Value(theLine, i, tabP); tabP(1)` (`gxx:231-233`).
fn first_point3d<L: LineTool>(line: &L, index: i32) -> GpPnt {
    line.value_3d(index)[0]
}

/// `LineTool::Value(theLine, i, tabP2d); tabP2d(1)` (`gxx:359-361`).
fn first_point2d<L: LineTool>(line: &L, index: i32) -> GpPnt2d {
    line.value_2d(index)[0]
}

/// `gp_Vec2d(P1, P2)` (`gp_Vec2d.hxx`).
fn vec2d_from_pnts(p1: &GpPnt2d, p2: &GpPnt2d) -> GpVec2d {
    GpVec2d::new(p2.x() - p1.x(), p2.y() - p1.y())
}

/// `gp_Vec2d::SquareMagnitude` (`gp_Vec2d.hxx`).
fn square_magnitude2d(v: &GpVec2d) -> f64 {
    v.x() * v.x() + v.y() * v.y()
}



