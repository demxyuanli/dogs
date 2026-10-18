//! `Extrema_GGenExtCC::Perform` — the general (non-elementary) curve/curve
//! extrema engine, i.e. what `Extrema_ECC` *is*
//! (`Extrema_ECC.hxx:23-28`: `using Extrema_ECC = Extrema_GGenExtCC<…>`), which
//! `Extrema_ExtCC` falls back to for every pair that is not an elementary
//! `Extrema_ExtElC` case (`Extrema_ExtCC.cxx:248-305`).
//!
//! Source: `Extrema_GGenExtCC.hxx` (920 lines) — `SetParams` (`:363-376`),
//! `Perform` (`:453-802`), the helpers `Extrema_GGenExtCC_comp` (`:135-150`),
//! `Extrema_GGenExtCC_ChangeIntervals` (`:152-202`),
//! `Extrema_GGenExtCC_PointsInspector` (`:204-246`) and
//! `Extrema_GGenExtCC_ProjPOnC` (`:248-263`), plus the accessors (`:855-918`).
//!
//! Grading of the parameter box follows OCCT exactly: `GeomAbs_C2` intervals
//! first (`hxx:462-464`), `GeomAbs_C1` when their product exceeds 100
//! (`hxx:466-471`), an optional length-ratio subdivision with `mult = 20`
//! (`hxx:473-531`), a 3-way split for closed single-interval curves
//! (`hxx:532-541`), then a `math_GlobOptMin` per interval pair
//! (`hxx:617-649`) with local parameters, the running best/`aSameTol*aValueTol`
//! bookkeeping (`hxx:650-667`) and a cell-filter dedup of the accepted points
//! (`hxx:665-686`), ending with the parallel-detection pass (`hxx:690-802`).
//!
//! Port deviations, all marked where they occur: the cell filter's *container*
//! is a plain accepted-point list (OCCT uses `NCollection_CellFilter`, a spatial
//! hash) while its acceptance semantics — no point within `aCellSize` — are
//! reproduced exactly: OCCT's query box is `[p - aCellSize, p + aCellSize]`
//! (`hxx:675-680`) with `aCellSize` also the cell size, so every point within
//! `aCellSize` of `p` lies in a queried cell and the dedup result is the same as
//! the global scan; the curve length used by the ratio test is
//! `gcpnts::curve_length` = `GCPnts_AbscissaPoint::Length(C)` (`A15/T-51`).

use super::prelude::*;
use super::GlobOptFuncCCC2;
use crate::extrema_pc::{extrema_ext_pc_range, ExtPcSolution};

use occt_core::precision::{Precision, INFINITE, PCONFUSION};
use occt_math::globoptmin::GlobOptMin;
use occt_math::vector::MathVector;

/// `GeomAbs_Shape` numeric values (OCCT enum order, as used by the port's
/// `Curve::nb_intervals`/`parameter_intervals`).
const GEOM_ABS_C1: u8 = 2;
const GEOM_ABS_C2: u8 = 4;

/// Port of `Extrema_GGenExtCC` (`Extrema_GGenExtCC.hxx:36-126`).
pub struct GGenExtCC<'a> {
    c1: &'a dyn Curve,
    c2: &'a dyn Curve,
    low_border: (f64, f64),
    upp_border: (f64, f64),
    /// `myCurveMinTol`, default `Precision::PConfusion()` (`hxx:283`).
    curve_min_tol: f64,
    /// `myIsFindSingleSolution` (`hxx:117`).
    single_solution: bool,
    done: bool,
    parallel: bool,
    points1: Vec<f64>,
    points2: Vec<f64>,
}

impl<'a> GGenExtCC<'a> {
    /// `Extrema_GGenExtCC()` + `SetParams` (`hxx:279-376`).
    pub fn new(c1: &'a dyn Curve, c2: &'a dyn Curve, low: (f64, f64), upp: (f64, f64)) -> Self {
        Self {
            c1,
            c2,
            low_border: low,
            upp_border: upp,
            curve_min_tol: PCONFUSION,
            single_solution: false,
            done: false,
            parallel: false,
            points1: Vec::new(),
            points2: Vec::new(),
        }
    }

    /// `SetTolerance` (`hxx:393-396`).
    pub fn set_tolerance(&mut self, tol: f64) {
        self.curve_min_tol = tol;
    }

    /// `SetSingleSolutionFlag` (`hxx:413-416`).
    pub fn set_single_solution_flag(&mut self, flag: bool) {
        self.single_solution = flag;
    }

    /// `GetSingleSolutionFlag`.
    pub fn single_solution_flag(&self) -> bool {
        self.single_solution
    }

    /// `IsDone()` (`hxx:813-822`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `IsParallel()` (`hxx:105`).
    pub fn is_parallel(&self) -> bool {
        self.parallel
    }

    /// `NbExt()` (`hxx:855-866`); OCCT throws when not done, the port returns 0.
    pub fn nb_ext(&self) -> usize {
        if !self.done {
            return 0;
        }
        self.points1.len()
    }

    /// `SquareDistance(theN)` (`hxx:877-892`), `theN` 1-based.
    pub fn square_distance(&self, n: usize) -> f64 {
        let u = self.points1[n - 1];
        let v = self.points2[n - 1];
        self.c1.d0(u).square_distance(&self.c2.d0(v))
    }

    /// `Points(theN, P1, P2)` (`hxx:903-918`): `(u, C1(u), v, C2(v))`.
    pub fn points(&self, n: usize) -> (f64, GpPnt, f64, GpPnt) {
        let u = self.points1[n - 1];
        let v = self.points2[n - 1];
        (u, self.c1.d0(u), v, self.c2.d0(v))
    }

    /// `Perform()` (`hxx:453-802`).
    pub fn perform(&mut self) -> Result<(), String> {
        self.done = false;
        self.parallel = false;
        let c1 = self.c1;
        let c2 = self.c2;

        // --- intervals (`hxx:461-471`) ------------------------------------
        // OCCT reads them off the *restricted* adaptors, whose range is the
        // border box (`Extrema_ExtCC.cxx:180` passes `myInf`/`mySup` as
        // `myLowBorder`/`myUppBorder`), so a bounded piece of an unbounded
        // curve yields the border interval, not `±Precision::Infinite()`.
        let mut cont = GEOM_ABS_C2;
        let mut iv1 = restricted_intervals(c1, cont, self.low_border.0, self.upp_border.0);
        let mut iv2 = restricted_intervals(c2, cont, self.low_border.1, self.upp_border.1);
        let mut nb1 = (iv1.len() - 1) as i32;
        let mut nb2 = (iv2.len() - 1) as i32;
        if nb1 * nb2 > 100 {
            cont = GEOM_ABS_C1;
            iv1 = restricted_intervals(c1, cont, self.low_border.0, self.upp_border.0);
            iv2 = restricted_intervals(c2, cont, self.low_border.1, self.upp_border.1);
            nb1 = (iv1.len() - 1) as i32;
            nb2 = (iv2.len() - 1) as i32;
        }

        // --- length-ratio subdivision (`hxx:473-511`) ---------------------
        let mut indmax: i32 = -1;
        let mut indmin: i32 = -1;
        let mult = 20.0f64;
        let mut an_l = [0.0f64; 2];
        if c1.first_parameter().is_finite()
            && c1.last_parameter().is_finite()
            && c2.first_parameter().is_finite()
            && c2.last_parameter().is_finite()
        {
            // `anL[0] = GCPnts_AbscissaPoint::Length(C1)` (`hxx:479-480`); the
            // port's `gcpnts::curve_length` is that same dispatch
            // (`GCPnts_AbscissaPoint.cxx:305-308`, A15/T-51): the type dispatch
            // at `:371-423` and `math_GaussSingleIntegration` on the speed.
            an_l[0] = crate::gcpnts::curve_length(c1);
            an_l[1] = crate::gcpnts::curve_length(c2);
            if an_l[0] / f64::from(nb1) > mult * an_l[1] / f64::from(nb2) {
                indmax = 0;
                indmin = 1;
            } else if an_l[1] / f64::from(nb2) > mult * an_l[0] / f64::from(nb1) {
                indmax = 1;
                indmin = 0;
            }
        }
        let mut nb_int_opt = 0i32;
        if indmax >= 0 {
            nb_int_opt =
                (an_l[indmax as usize] * f64::from([nb1, nb2][indmin as usize]) / an_l[indmin as usize] / (mult / 4.0)) as i32 + 1;
            if nb_int_opt > 100 || nb_int_opt < [nb1, nb2][indmax as usize] {
                indmax = -1;
            } else if nb_int_opt * [nb1, nb2][indmin as usize] > 100 {
                nb_int_opt = 100 / [nb1, nb2][indmin as usize];
                if nb_int_opt < [nb1, nb2][indmax as usize] {
                    indmax = -1;
                }
            }
        }

        // OCCT's `Adaptor3d_Curve` reports an unbounded curve as
        // `±Precision::Infinite()` (`Geom_Line::FirstParameter()`), a large but
        // *finite* number, so `math_GlobOptMin`'s midpoint `(A+B)/2` is 0; the
        // port's unbounded curves report IEEE ±inf, which would make that
        // midpoint NaN. Substitute the OCCT value (`precision::INFINITE`).
        for v in iv1.iter_mut().chain(iv2.iter_mut()) {
            if !v.is_finite() {
                *v = if *v < 0.0 { -INFINITE } else { INFINITE };
            }
        }
        if indmax >= 0 {
            if indmax == 0 {
                change_intervals(&mut iv1, nb_int_opt);
                nb1 = (iv1.len() - 1) as i32;
            } else {
                change_intervals(&mut iv2, nb_int_opt);
                nb2 = (iv2.len() - 1) as i32;
            }
        }
        if is_closed(c1) && nb1 == 1 {
            change_intervals(&mut iv1, 3);
            nb1 = (iv1.len() - 1) as i32;
        }
        if is_closed(c2) && nb2 == 1 {
            change_intervals(&mut iv2, 3);
            nb2 = (iv2.len() - 1) as i32;
        }

        // --- Lipschitz constant (`hxx:543-616`) ---------------------------
        let max_lc = 10000.0f64;
        let mut lc = 100.0f64;
        let max_der1 = 1.0 / c1.resolution(1.0);
        let max_der2 = 1.0 / c2.resolution(1.0);
        let mut max_der = max_der1.max(max_der2) * std::f64::consts::SQRT_2;
        if lc > max_der {
            lc = max_der;
        }
        let mut is_const_locked = false;
        const A_CR: f64 = 0.001;
        if max_der1 / max_der < A_CR || max_der2 / max_der < A_CR {
            is_const_locked = true;
        }
        if max_der > max_lc {
            lc = max_lc;
            is_const_locked = true;
        }
        if c1.gp_line().is_some() {
            max_der = 1.0 / c2.resolution(1.0);
            if lc > max_der {
                is_const_locked = true;
                lc = max_der;
            }
        }
        if c2.gp_line().is_some() {
            max_der = 1.0 / c1.resolution(1.0);
            if lc > max_der {
                is_const_locked = true;
                lc = max_der;
            }
        }

        let func = GlobOptFuncCCC2::new(c1, c2);
        if lc < max_lc || max_der > max_lc {
            let mut max_g = 0.0f64;
            let n1 = 21;
            let n2 = 21;
            let dt1 = (c1.last_parameter() - c1.first_parameter()) / f64::from(n1 - 1);
            let dt2 = (c2.last_parameter() - c2.first_parameter()) / f64::from(n2 - 1);
            let mut t1 = c1.first_parameter();
            for _i1 in 1..=n1 {
                let mut t2 = c2.first_parameter();
                for _i2 in 1..=n2 {
                    if let Some((_f, (g1, g2))) = func.values_uv(t1, t2) {
                        let a_mod = g1 * g1 + g2 * g2;
                        max_g = max_g.max(a_mod);
                    }
                    t2 += dt2;
                }
                t1 += dt1;
            }
            max_g = max_g.sqrt();
            if max_g > max_der {
                lc = max_g.min(max_lc);
                is_const_locked = true;
            }
            if max_g > 100.0 * max_lc {
                lc = 100.0 * max_lc;
                is_const_locked = true;
            } else if max_g < 0.1 * max_der {
                is_const_locked = true;
            }
        }

        // --- finder + interval sweep (`hxx:617-688`) ----------------------
        let mut finder = GlobOptMin::new();
        finder.set_lip_const_state(is_const_locked);
        finder.set_continuity(if cont == GEOM_ABS_C2 { 2 } else { 1 });
        let disc_tol = 1.0e-2;
        let value_tol = 1.0e-2;
        let same_tol = self.curve_min_tol / disc_tol;
        finder.set_tolerance(disc_tol, same_tol);
        finder.set_functional_minimal_value(0.0);
        finder.find_single_solution = self.single_solution;
        // `math_GlobOptMin aFinder(&aFunc, myLowBorder, myUppBorder, aLC)`
        // (`hxx:617`): the global box and `ComputeInitSol` are installed ONCE,
        // before the interval sweep; each pair then only `SetLocalParams` +
        // `Perform` (`hxx:648-649`) and the best value accumulates.
        // Same `Precision::Infinite()` substitution as the intervals above: OCCT
        // hands `math_GlobOptMin` the adaptors' `±Precision::Infinite()` bounds.
        let lo = [border(self.low_border.0), border(self.low_border.1)];
        let up = [border(self.upp_border.0), border(self.upp_border.1)];
        let lower = MathVector::from_slice(&lo);
        let upper = MathVector::from_slice(&up);
        finder.set_global_params(&func, &lower, &upper, lc)?;

        // `aCellSize` (`hxx:626-629`).
        let span = (iv1[iv1.len() - 1] - iv1[0]).max(iv2[iv2.len() - 1] - iv2[0]);
        let cell_size = (span * PCONFUSION / (2.0 * std::f64::consts::SQRT_2)).max(PCONFUSION);

        // The cell filter's container: OCCT uses `NCollection_CellFilter`, the
        // port an accepted-point list; the acceptance test (`Inspect`,
        // `hxx:231-240`, strict `< aCellSize^2`) is reproduced exactly and the
        // `hxx:675-680` query box `[p - aCellSize, p + aCellSize]` covers every
        // point the predicate can accept, so the dedup outcome is identical.
        let mut pnts: Vec<(f64, f64)> = Vec::new();
        let mut f = f64::MAX;

        for i in 0..nb1 as usize {
            for j in 0..nb2 as usize {
                let la = MathVector::from_slice(&[iv1[i], iv2[j]]);
                let lb = MathVector::from_slice(&[iv1[i + 1], iv2[j + 1]]);
                finder.set_local_params(&la, &lb);
                finder
                    .perform_local(&func)
                    .map_err(|e| format!("Extrema_GGenExtCC: {e}"))?;

                let curr_f = finder.minimal_value();
                if curr_f >= f + same_tol * value_tol {
                    continue;
                }
                if curr_f > f - same_tol * value_tol {
                    if curr_f < f {
                        f = curr_f;
                    }
                } else {
                    f = curr_f;
                    pnts.clear();
                }

                for k in 0..finder.nb_extrema() {
                    let sol = finder.point(k);
                    let pnt = (sol.value(1), sol.value(2));
                    // `PointsInspector::Inspect` (`hxx:231-240`).
                    let mut is_find = false;
                    for &(x, y) in &pnts {
                        let dx = pnt.0 - x;
                        let dy = pnt.1 - y;
                        if dx * dx + dy * dy < cell_size * cell_size {
                            is_find = true;
                            break;
                        }
                    }
                    if !is_find {
                        pnts.push(pnt);
                    }
                }
            }
        }

        // --- results (`hxx:690-802`) --------------------------------------
        let nb_sol = pnts.len();
        if nb_sol == 0 {
            self.done = false;
            return Ok(());
        }
        self.done = true;
        if nb_sol == 1 {
            self.points1.push(pnts[0].0);
            self.points2.push(pnts[0].1);
            return Ok(());
        }

        // `Extrema_GGenExtCC_comp` (`hxx:135-150`): by X, then by Y.
        pnts.sort_by(|a, b| {
            a.0.partial_cmp(&b.0)
                .unwrap_or(Ordering::Equal)
                .then(a.1.partial_cmp(&b.1).unwrap_or(Ordering::Equal))
        });

        let mut solutions: Vec<usize> = Vec::new();
        let mut b_save_solution = true;
        let mut b_dirs_coinside = true;
        let mut b_different_solutions = false;
        let mut is_parallel = true;
        for an_idx in 0..nb_sol - 1 {
            let current = pnts[an_idx];
            let next = pnts[an_idx + 1];
            let mid = ((current.0 + next.0) * 0.5, (current.1 + next.1) * 0.5);
            let val = func.value_uv(mid.0, mid.1).unwrap_or(f64::MAX);
            if (val - f).abs() < CONFUSION {
                if b_save_solution {
                    solutions.push(an_idx);
                    b_save_solution = false;
                }
            } else {
                is_parallel = false;
                solutions.push(an_idx);
                b_save_solution = true;
            }
            if !b_different_solutions {
                if next.0 > current.0 {
                    if next.1 > current.1 {
                        b_different_solutions = true;
                        b_dirs_coinside = true;
                    } else if next.1 < current.1 {
                        b_different_solutions = true;
                        b_dirs_coinside = false;
                    }
                }
            }
        }
        solutions.push(nb_sol - 1);

        if !b_different_solutions {
            is_parallel = false;
        }

        if is_parallel {
            // `hxx:765-783`: confirm with the point/curve projections.
            let a_t1 = [self.low_border.0, self.upp_border.0];
            let a_t2 = [
                if b_dirs_coinside { self.low_border.1 } else { self.upp_border.1 },
                if b_dirs_coinside { self.upp_border.1 } else { self.low_border.1 },
            ];
            for i_t in 0..2 {
                if !is_parallel {
                    break;
                }
                let dist1 = proj_p_on_c(c2, &c1.d0(a_t1[i_t]), self.low_border.1, self.upp_border.1);
                let dist2 = proj_p_on_c(c1, &c2.d0(a_t2[i_t]), self.low_border.0, self.upp_border.0);
                is_parallel = (dist1.min(dist2) - f * f).abs() < CONFUSION;
            }
        }

        if is_parallel {
            self.points1.push(pnts[0].0);
            self.points2.push(pnts[0].1);
            self.parallel = true;
        } else {
            for idx in solutions {
                self.points1.push(pnts[idx].0);
                self.points2.push(pnts[idx].1);
            }
        }
        Ok(())
    }
}

/// The interval array of `GeomAdaptor_Curve(C, lo, up)` for the continuity
/// `cont`: the curve's own intervals intersected with the border box (the
/// endpoints of the box are always kept).
fn restricted_intervals(c: &dyn Curve, cont: u8, lo: f64, up: f64) -> Vec<f64> {
    let (lo, up) = if lo <= up { (lo, up) } else { (up, lo) };
    let mut out: Vec<f64> = Vec::new();
    for v in c.parameter_intervals(cont) {
        let v = v.clamp(lo, up);
        if out.last().map_or(true, |&l| v > l) {
            out.push(v);
        }
    }
    if out.len() < 2 {
        return vec![lo, up];
    }
    if out[0] > lo {
        out.insert(0, lo);
    }
    let n = out.len();
    if out[n - 1] < up {
        out.push(up);
    }
    out
}

/// `Precision::Infinite()` substitution for a port curve bound (see `perform`).
fn border(v: f64) -> f64 {
    if v.is_finite() {
        v
    } else if v < 0.0 {
        -INFINITE
    } else {
        INFINITE
    }
}

/// `Extrema_GGenExtCC_ChangeIntervals` (`Extrema_GGenExtCC.hxx:152-202`).
fn change_intervals(ints: &mut Vec<f64>, nb_ints: i32) {
    let a_nb_ints = (ints.len() - 1) as i32;
    let mut nb_add = nb_ints - a_nb_ints;
    if a_nb_ints == 1 {
        let first = ints[0];
        let last = ints[ints.len() - 1];
        let dt = (last - first) / f64::from(nb_ints);
        let mut out = Vec::with_capacity(nb_ints as usize + 1);
        out.push(first);
        let mut t = first + dt;
        for _i in 2..=nb_ints {
            out.push(t);
            t += dt;
        }
        out.push(last);
        // OCCT writes `theNbInts + 1` entries with the last one at index
        // `theNbInts + 1`; the loop above fills indices 2..=theNbInts.
        *ints = out;
        return;
    }
    let mut out: Vec<f64> = ints.clone();
    let mut nb_last = ints.len() as i32;
    while nb_add > 0 {
        let mut an_l_int_max = -1.0f64;
        let mut max_ind = -1i32;
        for i in 1..nb_last {
            let an_l = out[i as usize] - out[(i - 1) as usize];
            if an_l > an_l_int_max {
                an_l_int_max = an_l;
                max_ind = i;
            }
        }
        let t = (out[max_ind as usize] + out[(max_ind - 1) as usize]) / 2.0;
        for i in (max_ind..nb_last).rev() {
            out[(i + 1) as usize] = out[i as usize];
        }
        nb_last += 1;
        nb_add -= 1;
        out[max_ind as usize] = t;
    }
    ints.clear();
    ints.extend_from_slice(&out[..nb_last as usize]);
}

/// `Extrema_GGenExtCC_ProjPOnC` (`Extrema_GGenExtCC.hxx:248-263`): the smallest
/// **squared** distance of `p` to any of the range-restricted `Extrema_ExtPC`
/// solutions, or `RealLast()` when there is none.
fn proj_p_on_c(c: &dyn Curve, p: &GpPnt, u1: f64, u2: f64) -> f64 {
    let sols: Vec<ExtPcSolution> = extrema_ext_pc_range(c, p, u1, u2);
    let mut dist = f64::MAX;
    for s in sols {
        if s.sq_dist < dist {
            dist = s.sq_dist;
        }
    }
    dist
}

/// `GeomAdaptor_Curve::IsClosed()` (`GeomAdaptor_Curve.cxx:576-585`): an
/// infinite range is never closed; otherwise the two ends coincide within
/// `Precision::Confusion()`. The curve here is the adaptor, not the underlying
/// `Geom_Curve`, so the predicate is endpoint coincidence -- a closed but
/// non-periodic B-spline is closed, while a periodic curve trimmed to an arc is
/// not. The 2-D twin is `extrema2d::curve2d_tool::adaptor_is_closed2d`
/// (`Geom2dAdaptor_Curve.cxx:588-600`).
fn is_closed(c: &dyn Curve) -> bool {
    let (first, last) = (c.first_parameter(), c.last_parameter());
    if Precision::is_positive_infinite(last) || Precision::is_negative_infinite(first) {
        return false;
    }
    c.d0(first).distance(&c.d0(last)) <= CONFUSION
}
