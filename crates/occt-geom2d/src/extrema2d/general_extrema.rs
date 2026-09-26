//! `Extrema_GGenExtCC::Perform` -- the general (non-elementary) curve/curve
//! extrema engine, i.e. what `Extrema_ECC2d` *is* (`Extrema_ECC2d.hxx:27-34`:
//! `using Extrema_ECC2d = Extrema_GGenExtCC<Adaptor2d_Curve2d, Extrema_Curve2dTool, ...>`),
//! which `Extrema_ExtCC2d` falls back to for every pair that is not an
//! elementary `Extrema_ExtElC2d` case (`Extrema_ExtCC2d.cxx:167-178`,
//! `:229-240`, `:295-306`, `:360-371`, `:416-427`, `:435-451`).
//!
//! Source: `Extrema_GGenExtCC.hxx` (920 lines) -- `SetParams` (`:363-376`),
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
//! hash) while its acceptance semantics -- no point within `aCellSize` -- are
//! reproduced; `Extrema_GGenExtCC_ProjPOnC` uses the port's point/curve extrema;
//! `GCPnts_AbscissaPoint::Length` in the ratio test is
//! `crate::curve_ops::curve2d_length` (the faithful port).

use super::prelude::*;

use super::curve2d_tool::{adaptor_intervals2d, adaptor_is_closed2d, adaptor_resolution2d};
use super::curve_curve::point_curve_extrema2d_all;
use super::glob_opt_func::GlobOptFuncCCC2;

use occt_math::globoptmin::GlobOptMin;
use occt_math::vector::MathVector;

/// `GeomAbs_Shape` numeric values (OCCT enum order, as used by the port's
/// `parameter_intervals`): `GeomAbs_C1` = 2, `GeomAbs_C2` = 4
/// (`occt_core::kernel::geomabs::Shape`).
const GEOM_ABS_C1: u8 = 2;
const GEOM_ABS_C2: u8 = 4;

/// Port of `Extrema_GGenExtCC` (`Extrema_GGenExtCC.hxx:36-126`), the 2D
/// instantiation `Extrema_ECC2d`.
pub(super) struct GGenExtCC<'a> {
    c1: &'a dyn Curve2d,
    c2: &'a dyn Curve2d,
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
    pub(super) fn new(
        c1: &'a dyn Curve2d,
        c2: &'a dyn Curve2d,
        low: (f64, f64),
        upp: (f64, f64),
    ) -> Self {
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
    pub(super) fn set_tolerance(&mut self, tol: f64) {
        self.curve_min_tol = tol;
    }

    /// `SetSingleSolutionFlag` (`hxx:413-416`).
    pub(super) fn set_single_solution_flag(&mut self, flag: bool) {
        self.single_solution = flag;
    }

    /// `GetSingleSolutionFlag`.
    pub(super) fn single_solution_flag(&self) -> bool {
        self.single_solution
    }

    /// `IsDone()` (`hxx:813-822`).
    pub(super) fn is_done(&self) -> bool {
        self.done
    }

    /// `IsParallel()` (`hxx:833-844`).
    pub(super) fn is_parallel(&self) -> bool {
        self.parallel
    }

    /// `NbExt()` (`hxx:855-866`); OCCT throws when not done, the port returns 0.
    pub(super) fn nb_ext(&self) -> usize {
        if !self.done {
            return 0;
        }
        self.points1.len()
    }

    /// `SquareDistance(theN)` (`hxx:877-892`), `theN` 1-based.
    pub(super) fn square_distance(&self, n: usize) -> f64 {
        let u = self.points1[n - 1];
        let v = self.points2[n - 1];
        self.c1.d0(u).square_distance(&self.c2.d0(v))
    }

    /// `Points(theN, P1, P2)` (`hxx:903-918`): `(u, C1(u), v, C2(v))`.
    pub(super) fn points(&self, n: usize) -> (f64, GpPnt2d, f64, GpPnt2d) {
        let u = self.points1[n - 1];
        let v = self.points2[n - 1];
        (u, self.c1.d0(u), v, self.c2.d0(v))
    }

    /// `Perform()` (`hxx:453-802`).
    pub(super) fn perform(&mut self) -> Result<(), String> {
        self.done = false;
        self.parallel = false;
        let c1 = self.c1;
        let c2 = self.c2;

        // --- intervals (`hxx:461-471`) ------------------------------------
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
            an_l[0] = crate::curve_ops::curve2d_length(c1); // `hxx:479`
            an_l[1] = crate::curve_ops::curve2d_length(c2); // `hxx:480`
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
            nb_int_opt = (an_l[indmax as usize] * f64::from([nb1, nb2][indmin as usize])
                / an_l[indmin as usize]
                / (mult / 4.0)) as i32
                + 1; // `hxx:495` (`RealToInt` is a saturating `as i32`)
            if nb_int_opt > 100 || nb_int_opt < [nb1, nb2][indmax as usize] {
                indmax = -1; // `hxx:496-499`
            } else if nb_int_opt * [nb1, nb2][indmin as usize] > 100 {
                nb_int_opt = 100 / [nb1, nb2][indmin as usize]; // `hxx:502-504`
                if nb_int_opt < [nb1, nb2][indmax as usize] {
                    indmax = -1;
                }
            }
        }

        // OCCT's `Adaptor2d_Curve2d` reports an unbounded curve through
        // `Precision::Infinite()` (=`2e100`, `Precision.hxx:371`) -- e.g.
        // `Geom2d_Line::FirstParameter()` -- so `math_GlobOptMin`'s midpoint
        // `(A+B)/2` is 0; the port's unbounded curves report IEEE +/-inf, which
        // would make that midpoint NaN and the ratio `0/0`. Substitute the OCCT
        // value (`occt_core::precision::INFINITE`).
        for v in iv1.iter_mut().chain(iv2.iter_mut()) {
            if !v.is_finite() {
                *v = if *v < 0.0 { -INFINITE } else { INFINITE };
            }
        }
        if indmax >= 0 {
            if indmax == 0 {
                change_intervals(&mut iv1, nb_int_opt); // `hxx:523`
                nb1 = (iv1.len() - 1) as i32;
            } else {
                change_intervals(&mut iv2, nb_int_opt); // `hxx:528`
                nb2 = (iv2.len() - 1) as i32;
            }
        }
        if adaptor_is_closed2d(c1) && nb1 == 1 {
            change_intervals(&mut iv1, 3); // `hxx:532-536`
            nb1 = (iv1.len() - 1) as i32;
        }
        if adaptor_is_closed2d(c2) && nb2 == 1 {
            change_intervals(&mut iv2, 3); // `hxx:537-541`
            nb2 = (iv2.len() - 1) as i32;
        }

        // --- Lipschitz constant (`hxx:543-616`) ---------------------------
        const A_MAX_LC: f64 = 10000.0; // `hxx:543`
        let mut a_lc = 100.0f64; // `hxx:544`
        let a_max_der1 = 1.0 / adaptor_resolution2d(c1, 1.0); // `hxx:545`
        let a_max_der2 = 1.0 / adaptor_resolution2d(c2, 1.0); // `hxx:546`
        let mut a_max_der = a_max_der1.max(a_max_der2) * std::f64::consts::SQRT_2; // `hxx:547`
        if a_lc > a_max_der {
            a_lc = a_max_der; // `hxx:548-549`
        }
        let mut is_const_locked = false;
        const A_CR: f64 = 0.001; // `hxx:552`
        if a_max_der1 / a_max_der < A_CR || a_max_der2 / a_max_der < A_CR {
            is_const_locked = true; // `hxx:553-556`
        }
        if a_max_der > A_MAX_LC {
            a_lc = A_MAX_LC;
            is_const_locked = true; // `hxx:557-561`
        }
        if c1.is_line() {
            a_max_der = 1.0 / adaptor_resolution2d(c2, 1.0); // `hxx:562-564`
            if a_lc > a_max_der {
                is_const_locked = true;
                a_lc = a_max_der; // `hxx:565-569`
            }
        }
        if c2.is_line() {
            a_max_der = 1.0 / adaptor_resolution2d(c1, 1.0); // `hxx:571-573`
            if a_lc > a_max_der {
                is_const_locked = true;
                a_lc = a_max_der; // `hxx:574-578`
            }
        }

        let func = GlobOptFuncCCC2::new(c1, c2); // `hxx:581`
        if a_lc < A_MAX_LC || a_max_der > A_MAX_LC {
            let mut a_max_g = 0.0f64;
            let n1 = 21; // `hxx:587`
            let n2 = 21;
            let dt1 = (c1.last_parameter() - c1.first_parameter()) / f64::from(n1 - 1);
            let dt2 = (c2.last_parameter() - c2.first_parameter()) / f64::from(n2 - 1);
            let mut t1 = c1.first_parameter();
            for _i1 in 1..=n1 {
                let mut t2 = c2.first_parameter();
                for _i2 in 1..=n2 {
                    if let Some((_f, (g1, g2))) = func.values(t1, t2) {
                        let a_mod = g1 * g1 + g2 * g2;
                        a_max_g = a_max_g.max(a_mod); // `hxx:597-598`
                    }
                    t2 += dt2;
                }
                t1 += dt1;
            }
            a_max_g = a_max_g.sqrt(); // `hxx:601`
            if a_max_g > a_max_der {
                a_lc = a_max_g.min(A_MAX_LC);
                is_const_locked = true; // `hxx:602-606`
            }
            if a_max_g > 100.0 * A_MAX_LC {
                a_lc = 100.0 * A_MAX_LC;
                is_const_locked = true; // `hxx:607-611`
            } else if a_max_g < 0.1 * a_max_der {
                is_const_locked = true; // `hxx:612-615`
            }
        }

        // --- finder + interval sweep (`hxx:617-688`) ----------------------
        let mut finder = GlobOptMin::new();
        finder.set_lip_const_state(is_const_locked); // `hxx:618`
        finder.set_continuity(if cont == GEOM_ABS_C2 { 2 } else { 1 }); // `hxx:619`
        let a_disc_tol = 1.0e-2; // `hxx:620`
        let a_value_tol = 1.0e-2; // `hxx:621`
        let a_same_tol = self.curve_min_tol / a_disc_tol; // `hxx:622`
        finder.set_tolerance(a_disc_tol, a_same_tol); // `hxx:623`
        finder.set_functional_minimal_value(0.0); // `hxx:624`
        finder.find_single_solution = self.single_solution;
        // `math_GlobOptMin aFinder(&aFunc, myLowBorder, myUppBorder, aLC)`
        // (`hxx:617`): the global box and `ComputeInitSol` are installed ONCE,
        // before the interval sweep; each pair then only `SetLocalParams` +
        // `Perform` (`hxx:648-649`) and the best value accumulates. Same
        // `Precision::Infinite()` substitution as the intervals above.
        let lo = [border(self.low_border.0), border(self.low_border.1)];
        let up = [border(self.upp_border.0), border(self.upp_border.1)];
        let lower = MathVector::from_slice(&lo);
        let upper = MathVector::from_slice(&up);
        let obj = |x: &MathVector| func.value(x.value(1), x.value(2)).unwrap_or(f64::MAX);
        finder.set_global_params(&obj, &lower, &upper, a_lc)?;

        // `aCellSize` (`hxx:626-629`).
        let span = (iv1[iv1.len() - 1] - iv1[0]).max(iv2[iv2.len() - 1] - iv2[0]);
        let a_cell_size = (span * PCONFUSION / (2.0 * std::f64::consts::SQRT_2)).max(PCONFUSION);

        // The cell filter's container: OCCT uses `NCollection_CellFilter`, the
        // port an accepted-point list; the acceptance test (`Inspect`,
        // `hxx:231-240`) is the same squared-distance comparison.
        let mut a_pnts: Vec<(f64, f64)> = Vec::new();
        let mut a_f = f64::MAX;

        for i in 0..nb1 as usize {
            for j in 0..nb2 as usize {
                let la = MathVector::from_slice(&[iv1[i], iv2[j]]);
                let lb = MathVector::from_slice(&[iv1[i + 1], iv2[j + 1]]);
                finder.set_local_params(&la, &lb); // `hxx:648`
                finder
                    .perform_local(&obj) // `hxx:649`
                    .map_err(|e| format!("Extrema_GGenExtCC: {e}"))?;

                let a_curr_f = finder.minimal_value(); // `hxx:651`
                if a_curr_f >= a_f + a_same_tol * a_value_tol {
                    continue; // `hxx:652-655`
                }
                if a_curr_f > a_f - a_same_tol * a_value_tol {
                    if a_curr_f < a_f {
                        a_f = a_curr_f; // `hxx:657-662`
                    }
                } else {
                    a_f = a_curr_f; // `hxx:663-667`
                    a_pnts.clear();
                }

                for k in 0..finder.nb_extrema() {
                    let sol = finder.point(k); // `hxx:672`
                    let pnt = (sol.value(1), sol.value(2));
                    // `PointsInspector::Inspect` (`hxx:231-240`).
                    let mut is_find = false;
                    for &(x, y) in &a_pnts {
                        let dx = pnt.0 - x;
                        let dy = pnt.1 - y;
                        if dx * dx + dy * dy < a_cell_size * a_cell_size {
                            is_find = true;
                            break;
                        }
                    }
                    if !is_find {
                        a_pnts.push(pnt); // `hxx:681-685`
                    }
                }
            }
        }

        // --- results (`hxx:690-802`) --------------------------------------
        let nb_sol = a_pnts.len();
        if nb_sol == 0 {
            self.done = false; // `hxx:691-695`
            return Ok(());
        }
        self.done = true; // `hxx:697`
        if nb_sol == 1 {
            self.points1.push(a_pnts[0].0); // `hxx:699-705`
            self.points2.push(a_pnts[0].1);
            return Ok(());
        }

        // `Extrema_GGenExtCC_comp` (`hxx:135-150`): by X, then by Y.
        a_pnts.sort_by(|a, b| {
            a.0.partial_cmp(&b.0)
                .unwrap_or(Ordering::Equal)
                .then(a.1.partial_cmp(&b.1).unwrap_or(Ordering::Equal))
        });

        let mut a_solutions: Vec<usize> = Vec::new(); // `hxx:709`
        let mut b_save_solution = true;
        let mut b_dirs_coinside = true;
        let mut b_different_solutions = false;
        let mut is_parallel = true;
        for an_idx in 0..nb_sol - 1 {
            let a_current = a_pnts[an_idx];
            let a_next = a_pnts[an_idx + 1];
            let mid = ((a_current.0 + a_next.0) * 0.5, (a_current.1 + a_next.1) * 0.5);
            let a_val = func.value(mid.0, mid.1).unwrap_or(f64::MAX); // `hxx:727`
            if (a_val - a_f).abs() < CONFUSION {
                if b_save_solution {
                    a_solutions.push(an_idx); // `hxx:728-735`
                    b_save_solution = false;
                }
            } else {
                is_parallel = false; // `hxx:736-741`
                a_solutions.push(an_idx);
                b_save_solution = true;
            }
            if !b_different_solutions {
                if a_next.0 > a_current.0 {
                    if a_next.1 > a_current.1 {
                        b_different_solutions = true;
                        b_dirs_coinside = true; // `hxx:743-758`
                    } else if a_next.1 < a_current.1 {
                        b_different_solutions = true;
                        b_dirs_coinside = false;
                    }
                }
            }
        }
        a_solutions.push(nb_sol - 1); // `hxx:760`

        if !b_different_solutions {
            is_parallel = false; // `hxx:762-763`
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
                let a_dist1 = proj_p_on_c(c2, &c1.d0(a_t1[i_t]), self.low_border.1, self.upp_border.1);
                let a_dist2 = proj_p_on_c(c1, &c2.d0(a_t2[i_t]), self.low_border.0, self.upp_border.0);
                is_parallel = (a_dist1.min(a_dist2) - a_f * a_f).abs() < CONFUSION; // `hxx:781`
            }
        }

        if is_parallel {
            self.points1.push(a_pnts[0].0); // `hxx:787-790`
            self.points2.push(a_pnts[0].1);
            self.parallel = true;
        } else {
            for idx in a_solutions {
                self.points1.push(a_pnts[idx].0); // `hxx:797-799`
                self.points2.push(a_pnts[idx].1);
            }
        }
        Ok(())
    }
}

/// The interval array of `Extrema_Curve2dTool::Intervals` for the continuity
/// `cont`: the curve's own break points clipped to the border box (the
/// endpoints of the box are always kept).
///
/// `Extrema_ECC2d(C1, C2)` (`Extrema_ECC2d.hxx:27-34`) is built from the two
/// adaptors, so `NbIntervals`/`Intervals` run over the adaptors' own range
/// (`Extrema_GGenExtCC.hxx:463-464`); this clip is what makes the engine usable
/// with explicit borders too (`SetParams`, `hxx:363-376`).
fn restricted_intervals(c: &dyn Curve2d, cont: u8, lo: f64, up: f64) -> Vec<f64> {
    let (lo, up) = if lo <= up { (lo, up) } else { (up, lo) };
    let mut out: Vec<f64> = Vec::new();
    for v in adaptor_intervals2d(c, cont) {
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

/// `Extrema_GGenExtCC_ChangeIntervals` (`Extrema_GGenExtCC.hxx:152-202`):
/// re-subdivide `the_ints` into `the_nb_ints` intervals.
///
/// OCCT allocates `NCollection_HArray1(1, theNbInts + 1)` up front (`hxx:157-158`)
/// and copies `aNbLast` old values into it; the port sizes the working vector
/// the same way, since the insertion loop writes one position past the current
/// last index.
fn change_intervals(the_ints: &mut Vec<f64>, the_nb_ints: i32) {
    let a_nb_ints = (the_ints.len() - 1) as i32;
    let mut a_nb_add = the_nb_ints - a_nb_ints;
    if a_nb_ints == 1 {
        let first = the_ints[0];
        let last = the_ints[the_ints.len() - 1];
        let dt = (last - first) / f64::from(the_nb_ints);
        let mut out = Vec::with_capacity(the_nb_ints as usize + 1);
        out.push(first); // `hxx:163`
        let mut t = first + dt;
        for _i in 2..=the_nb_ints {
            out.push(t); // `hxx:167-170`
            t += dt;
        }
        out.push(last); // `hxx:164`
        *the_ints = out;
        return;
    }
    // `aNewInts` has `theNbInts + 1` slots (`hxx:157-158`); `aNewInts(1..aNbLast)`
    // is the copy of `theInts` (`hxx:174-177`).
    let mut out: Vec<f64> = the_ints.clone();
    out.resize(the_nb_ints as usize + 1, f64::NAN);
    let mut a_nb_last = the_ints.len() as i32;
    while a_nb_add > 0 {
        let mut an_l_int_max = -1.0f64;
        let mut a_max_ind = -1i32;
        for i in 1..a_nb_last {
            let an_l = out[i as usize] - out[(i - 1) as usize]; // `hxx:182-190`
            if an_l > an_l_int_max {
                an_l_int_max = an_l;
                a_max_ind = i;
            }
        }
        let t = (out[a_max_ind as usize] + out[(a_max_ind - 1) as usize]) / 2.0; // `hxx:192`
        for i in (a_max_ind..a_nb_last).rev() {
            out[(i + 1) as usize] = out[i as usize]; // `hxx:193-196`
        }
        a_nb_last += 1;
        a_nb_add -= 1;
        out[a_max_ind as usize] = t; // `hxx:199`
    }
    the_ints.clear();
    the_ints.extend_from_slice(&out[..a_nb_last as usize]);
}

/// `Extrema_GGenExtCC_ProjPOnC` (`Extrema_GGenExtCC.hxx:248-263`): the smallest
/// **squared** distance of `the_p` to any range-restricted `Extrema_ExtPC2d`
/// solution, or `RealLast()` when there is none.
///
/// **UNPORTED (audit A16)**: OCCT calls `Extrema_ExtPC2d`
/// (`Extrema_ExtPC2d.hxx:31-38` = `Extrema_GGExtPC<...>`) initialized on the
/// restricted range (`Extrema_GGenExtCC.hxx:771-773`); the port's point/curve
/// extrema (`point_curve_extrema2d_all`, itself the A7 substitute for
/// `Extrema_ExtPC2d`) has no range-restricted form, so its candidates inside
/// `[u1, u2]` are used and an empty set stays `RealLast()`.
fn proj_p_on_c(c: &dyn Curve2d, p: &GpPnt2d, u1: f64, u2: f64) -> f64 {
    let mut a_dist = f64::MAX;
    for e in point_curve_extrema2d_all(c, p) {
        if e.u1 >= u1 - PCONFUSION && e.u1 <= u2 + PCONFUSION {
            a_dist = a_dist.min(e.p2.square_distance(p));
        }
    }
    a_dist
}