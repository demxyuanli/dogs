//! `Approx_FitAndDivide2d` - the degree-escalating divide-and-fit driver for a
//! [`AppCont_Function`](crate::app_cont::Function).
//!
//! Source: `ModelingData/TKGeomBase/Approx/Approx_ComputeCLine.gxx:1-486`
//! instantiated as `Approx_FitAndDivide2d` by
//! `Approx_FitAndDivide2d_0.cxx:18-22` (`MultiLine = AppCont_Function`,
//! `MultiLineTool` absent, `Approx_ComputeCLine` renamed). The class interface
//! is `Approx_FitAndDivide2d.hxx:30-136` (OCCT 8.0.0).
//!
//! The 3d twin `Approx_FitAndDivide` (`Approx_FitAndDivide_0.cxx`) is the same
//! `.gxx` with `MultiLine = AppDef_MultiLine`; only the 2d instantiation is
//! ported here because it is the one `ProjLib_ComputeApprox` uses
//! (`ProjLib_ComputeApprox.cxx:1323-1329`).
//!
//! Translation notes:
//! - `NCollection_Sequence<T>` becomes `Vec<T>`; the 1-based `Value(Index)` of
//!   the sequence is the `[index - 1]` element.
//! - `RealLast()` is `f64::MAX`, `RealEpsilon()` is
//!   [`REAL_EPSILON`](occt_core::precision::REAL_EPSILON).
//! - `AppCont_LeastSquare::Value()` is not a pure getter (it rebuilds
//!   `mySCU` from `myPoles`, `AppCont_LeastSquare.cxx:518-548`), so the port
//!   keeps the local least-square object mutable and clones the result.

use occt_core::precision::{CONFUSION, PAPPROXIMATION, REAL_EPSILON};

use crate::app_cont::function::Function;
use crate::app_cont::least_square::LeastSquare;
use crate::app_par_curves::constraint::Constraint;
use crate::app_par_curves::multi_curve::MultiCurve;

/// `MAXSEGM` (`gxx:27`).
const MAXSEGM: i32 = 1000;

/// `Approx_FitAndDivide2d` (`Approx_FitAndDivide2d.hxx:30-136`).
#[derive(Debug, Clone)]
pub struct FitAndDivide2d {
    /// `myMultiCurves` (`hxx:115`).
    my_multi_curves: Vec<MultiCurve>,
    /// `myfirstparam` (`hxx:116`).
    my_first_param: Vec<f64>,
    /// `mylastparam` (`hxx:117`).
    my_last_param: Vec<f64>,
    /// `TheMultiCurve` (`hxx:118`).
    the_multi_curve: MultiCurve,
    /// `alldone` (`hxx:119`).
    alldone: bool,
    /// `tolreached` (`hxx:120`).
    tolreached: bool,
    /// `Tolers3d` (`hxx:121`).
    tol_ers3d: Vec<f64>,
    /// `Tolers2d` (`hxx:122`).
    tol_ers2d: Vec<f64>,
    /// `mydegremin` (`hxx:123`).
    my_degremin: i32,
    /// `mydegremax` (`hxx:124`).
    my_degremax: i32,
    /// `mytol3d` (`hxx:125`).
    my_tol3d: f64,
    /// `mytol2d` (`hxx:126`).
    my_tol2d: f64,
    /// `currenttol3d` (`hxx:127`).
    current_tol3d: f64,
    /// `currenttol2d` (`hxx:128`).
    current_tol2d: f64,
    /// `mycut` (`hxx:129`).
    my_cut: bool,
    /// `myfirstC` (`hxx:130`).
    my_first_c: Constraint,
    /// `mylastC` (`hxx:131`).
    my_last_c: Constraint,
    /// `myMaxSegments` (`hxx:132`).
    my_max_segments: i32,
    /// `myInvOrder` (`hxx:133`).
    my_inv_order: bool,
    /// `myHangChecking` (`hxx:134`).
    my_hang_checking: bool,
}

impl FitAndDivide2d {
    /// `Approx_FitAndDivide2d(Line, degreemin, degreemax, Tolerance3d,
    /// Tolerance2d, cutting, FirstC, LastC)` (`hxx:39-47`, `gxx:37-61`).
    #[allow(clippy::too_many_arguments)]
    pub fn with_line<F: Function + ?Sized>(
        line: &F,
        degreemin: i32,
        degreemax: i32,
        tolerance3d: f64,
        tolerance2d: f64,
        cutting: bool,
        first_c: Constraint,
        last_c: Constraint,
    ) -> Self {
        let mut this = Self::from_settings(
            degreemin,
            degreemax,
            tolerance3d,
            tolerance2d,
            cutting,
            first_c,
            last_c,
        );
        this.perform(line);
        this
    }

    /// `Approx_FitAndDivide2d(degreemin, degreemax, Tolerance3d, Tolerance2d,
    /// cutting, FirstC, LastC)` (`hxx:50-58`, `gxx:67-90`): no `Perform`.
    #[allow(clippy::too_many_arguments)]
    pub fn from_settings(
        degreemin: i32,
        degreemax: i32,
        tolerance3d: f64,
        tolerance2d: f64,
        cutting: bool,
        first_c: Constraint,
        last_c: Constraint,
    ) -> Self {
        Self {
            my_multi_curves: Vec::new(),
            my_first_param: Vec::new(),
            my_last_param: Vec::new(),
            the_multi_curve: MultiCurve::new(),
            alldone: false,
            tolreached: false,
            tol_ers3d: Vec::new(),
            tol_ers2d: Vec::new(),
            my_degremin: degreemin,
            my_degremax: degreemax,
            my_tol3d: tolerance3d,
            my_tol2d: tolerance2d,
            current_tol3d: 0.0,
            current_tol2d: 0.0,
            my_cut: cutting,
            my_first_c: first_c,
            my_last_c: last_c,
            my_max_segments: MAXSEGM,
            my_inv_order: true,
            my_hang_checking: true,
        }
    }

    /// `Perform(Line)` (`gxx:96-229`).
    pub fn perform<F: Function + ?Sized>(&mut self, line: &F) {
        let u_first = line.first_parameter();
        let u_last = line.last_parameter();
        let mut begin = true;
        let mut ok = false;
        let mut thetol3d = CONFUSION;
        let mut thetol2d = CONFUSION;
        let tol_u = if self.my_hang_checking {
            ((u_last - u_first) * 1.0e-3).max(CONFUSION)
        } else {
            ((u_last - u_first) * 1.0e-5).max(PAPPROXIMATION)
        };
        let mut myfirst_u = u_first;
        let mut mylast_u = u_last;
        let mut a_max_segments = 0i32;
        let a_max_segments1 = self.my_max_segments - 1;
        let mut a_nb_cut = 0i32;
        let mut a_nb_imp = 0i32;
        let a_nb_comp = 10i32;

        if !self.my_cut {
            self.alldone = self.compute(line, u_first, u_last, &mut thetol3d, &mut thetol2d);
            if !self.alldone {
                self.tolreached = false;
                self.my_first_param.push(u_first);
                self.my_last_param.push(u_last);
                self.my_multi_curves.push(self.the_multi_curve.clone());
                self.tol_ers3d.push(self.current_tol3d);
                self.tol_ers2d.push(self.current_tol2d);
            }
        } else {
            // previous decision to be taken if we get worse with next cut (eap)
            let mut kept_multi_curve = MultiCurve::new();
            let mut kept_ufirst = 0.0f64;
            let mut kept_ulast = 0.0f64;
            let mut kept_t3d = f64::MAX;
            let mut kept_t2d = 0.0f64;

            // `while (!Finish)`: `Finish` is only ever set right before the
            // `return` (`gxx:153-156`), so the loop's single exit is that
            // return and the port keeps it as a plain `loop`.
            loop {
                // Management of the multiline splitting for approximation:
                if !begin {
                    if ok {
                        // Calcul de la partie a approximer.
                        myfirst_u = mylast_u;
                        mylast_u = u_last;
                        a_nb_cut = 0;
                        a_nb_imp = 0;
                        if (u_last - myfirst_u).abs() <= REAL_EPSILON
                            || a_max_segments >= self.my_max_segments
                        {
                            self.alldone = true;
                            return;
                        }
                        kept_t3d = f64::MAX;
                        kept_t2d = 0.0;
                        kept_ufirst = myfirst_u;
                        kept_ulast = mylast_u;
                    } else {
                        // keep best decision
                        if (thetol3d + thetol2d) < (kept_t3d + kept_t2d) {
                            kept_multi_curve = self.the_multi_curve.clone();
                            kept_ufirst = myfirst_u;
                            kept_ulast = mylast_u;
                            kept_t3d = thetol3d;
                            kept_t2d = thetol2d;
                            a_nb_imp += 1;
                        }

                        // cut an interval
                        mylast_u = (myfirst_u + mylast_u) / 2.0;
                        a_nb_cut += 1;
                    }
                }

                // Calculation of parameters on this new interval.
                ok = self.compute(line, myfirst_u, mylast_u, &mut thetol3d, &mut thetol2d);
                if ok {
                    a_max_segments += 1;
                }

                let mut a_stop_cutting = false;
                if self.my_hang_checking && a_nb_cut >= a_nb_comp {
                    if a_nb_cut > a_nb_imp + 1 {
                        a_stop_cutting = true;
                    }
                    a_nb_cut = 0;
                    a_nb_imp = 0;
                }
                // is new decision better?
                if !ok
                    && ((myfirst_u - mylast_u).abs() <= tol_u
                        || a_max_segments >= a_max_segments1
                        || a_stop_cutting)
                {
                    ok = true; // stop interval cutting, approx the rest part

                    if (thetol3d + thetol2d) < (kept_t3d + kept_t2d) {
                        kept_multi_curve = self.the_multi_curve.clone();
                        kept_ufirst = myfirst_u;
                        kept_ulast = mylast_u;
                        kept_t3d = thetol3d;
                        kept_t2d = thetol2d;
                    }

                    mylast_u = kept_ulast;

                    self.tolreached = false; // helas
                    self.my_multi_curves.push(kept_multi_curve.clone());
                    a_max_segments += 1;
                    self.tol_ers3d.push(kept_t3d);
                    self.tol_ers2d.push(kept_t2d);
                    self.my_first_param.push(kept_ufirst);
                    self.my_last_param.push(kept_ulast);
                }

                begin = false;
            }
        }
    }

    /// `Compute(Line, Ufirst, Ulast, TheTol3d, TheTol2d)` (`gxx:257-382`).
    fn compute<F: Function + ?Sized>(
        &mut self,
        line: &F,
        ufirst: f64,
        ulast: f64,
        the_tol3d: &mut f64,
        the_tol2d: &mut f64,
    ) -> bool {
        let nb_points_max = 24i32;
        let a_min_ratio = 0.05;
        let a_max_deg = 8i32;
        let mut _fv = 0.0f64;

        let mut a_prev_curve = MultiCurve::new();
        let mut a_prev_tol3d = f64::MAX;
        let mut a_prev_tol2d = f64::MAX;
        let mut a_prev_is_ok = false;
        let mut an_inv_order = self.my_inv_order;
        if an_inv_order && self.my_degremax > a_max_deg {
            if (ulast - ufirst) / (line.last_parameter() - line.first_parameter()) < a_min_ratio {
                an_inv_order = false;
            }
        }
        if an_inv_order {
            let mut deg = self.my_degremax;
            while deg >= self.my_degremin {
                let nb_points = (2 * deg + 1).min(nb_points_max);
                let mut l_square = LeastSquare::new(
                    line,
                    ufirst,
                    ulast,
                    self.my_first_c,
                    self.my_last_c,
                    deg,
                    nb_points,
                );
                let mydone = l_square.is_done();
                if mydone {
                    let (f, t3d, t2d) = l_square.error();
                    _fv = f;
                    *the_tol3d = t3d;
                    *the_tol2d = t2d;
                    if *the_tol3d <= self.my_tol3d && *the_tol2d <= self.my_tol2d {
                        if deg == self.my_degremin {
                            // Stockage de la multicurve approximee.
                            self.tolreached = true;
                            self.my_multi_curves.push(l_square.value().clone());
                            self.my_first_param.push(ufirst);
                            self.my_last_param.push(ulast);
                            self.tol_ers3d.push(*the_tol3d);
                            self.tol_ers2d.push(*the_tol2d);
                            return true;
                        }
                        a_prev_tol3d = *the_tol3d;
                        a_prev_tol2d = *the_tol2d;
                        a_prev_curve = l_square.value().clone();
                        a_prev_is_ok = true;
                        deg -= 1;
                        continue;
                    } else if a_prev_is_ok {
                        // Stockage de la multicurve approximee.
                        self.tolreached = true;
                        *the_tol3d = a_prev_tol3d;
                        *the_tol2d = a_prev_tol2d;
                        self.my_multi_curves.push(a_prev_curve.clone());
                        self.my_first_param.push(ufirst);
                        self.my_last_param.push(ulast);
                        self.tol_ers3d.push(a_prev_tol3d);
                        self.tol_ers2d.push(a_prev_tol2d);
                        return true;
                    }
                } else if a_prev_is_ok {
                    // Stockage de la multicurve approximee.
                    self.tolreached = true;
                    *the_tol3d = a_prev_tol3d;
                    *the_tol2d = a_prev_tol2d;
                    self.my_multi_curves.push(a_prev_curve.clone());
                    self.my_first_param.push(ufirst);
                    self.my_last_param.push(ulast);
                    self.tol_ers3d.push(a_prev_tol3d);
                    self.tol_ers2d.push(a_prev_tol2d);
                    return true;
                }
                if !a_prev_is_ok && deg == self.my_degremax {
                    self.the_multi_curve = l_square.value().clone();
                    self.current_tol3d = *the_tol3d;
                    self.current_tol2d = *the_tol2d;
                    // `aPrev*` / `aPrevCurve` are re-seeded here (`gxx:344-347`)
                    // but the loop breaks immediately after, so the writes are
                    // unobservable and the port drops them.
                    break;
                }
                deg -= 1;
            }
        } else {
            let mut deg = self.my_degremin;
            while deg <= self.my_degremax {
                let nb_points = (2 * deg + 1).min(nb_points_max);
                let mut l_square = LeastSquare::new(
                    line,
                    ufirst,
                    ulast,
                    self.my_first_c,
                    self.my_last_c,
                    deg,
                    nb_points,
                );
                let mydone = l_square.is_done();
                if mydone {
                    let (f, t3d, t2d) = l_square.error();
                    _fv = f;
                    *the_tol3d = t3d;
                    *the_tol2d = t2d;
                    if *the_tol3d <= self.my_tol3d && *the_tol2d <= self.my_tol2d {
                        // Stockage de la multicurve approximee.
                        self.tolreached = true;
                        self.my_multi_curves.push(l_square.value().clone());
                        self.my_first_param.push(ufirst);
                        self.my_last_param.push(ulast);
                        self.tol_ers3d.push(*the_tol3d);
                        self.tol_ers2d.push(*the_tol2d);
                        return true;
                    }
                }
                if deg == self.my_degremax {
                    self.the_multi_curve = l_square.value().clone();
                    self.current_tol3d = *the_tol3d;
                    self.current_tol2d = *the_tol2d;
                }
                deg += 1;
            }
        }
        false
    }

    /// `NbMultiCurves()` (`gxx:236-239`).
    pub fn nb_multi_curves(&self) -> usize {
        self.my_multi_curves.len()
    }

    /// `Value(Index)` (`gxx:246-250`).
    pub fn value(&self, index: i32) -> &MultiCurve {
        &self.my_multi_curves[(index - 1) as usize]
    }

    /// `Parameters(Index)` (`gxx:390-394`): `(firstpar, lastpar)`.
    pub fn parameters(&self, index: i32) -> (f64, f64) {
        let i = (index - 1) as usize;
        (self.my_first_param[i], self.my_last_param[i])
    }

    /// `SetDegrees` (`gxx:402-406`).
    pub fn set_degrees(&mut self, degreemin: i32, degreemax: i32) {
        self.my_degremin = degreemin;
        self.my_degremax = degreemax;
    }

    /// `SetTolerances` (`gxx:412-416`).
    pub fn set_tolerances(&mut self, tolerance3d: f64, tolerance2d: f64) {
        self.my_tol3d = tolerance3d;
        self.my_tol2d = tolerance2d;
    }

    /// `SetConstraints` (`gxx:423-428`).
    pub fn set_constraints(&mut self, first_c: Constraint, last_c: Constraint) {
        self.my_first_c = first_c;
        self.my_last_c = last_c;
    }

    /// `SetMaxSegments` (`gxx:435-438`).
    pub fn set_max_segments(&mut self, the_max_segments: i32) {
        self.my_max_segments = the_max_segments;
    }

    /// `SetInvOrder` (`gxx:443-446`).
    pub fn set_inv_order(&mut self, the_inv_order: bool) {
        self.my_inv_order = the_inv_order;
    }

    /// `SetHangChecking` (`gxx:451-454`).
    pub fn set_hang_checking(&mut self, the_hang_checking: bool) {
        self.my_hang_checking = the_hang_checking;
    }

    /// `IsAllApproximated` (`gxx:462-465`).
    pub fn is_all_approximated(&self) -> bool {
        self.alldone
    }

    /// `IsToleranceReached` (`gxx:472-475`).
    pub fn is_tolerance_reached(&self) -> bool {
        self.tolreached
    }

    /// `Error(Index, tol3d, tol2d)` (`gxx:482-486`).
    pub fn error(&self, index: i32) -> (f64, f64) {
        let i = (index - 1) as usize;
        (self.tol_ers3d[i], self.tol_ers2d[i])
    }
}
