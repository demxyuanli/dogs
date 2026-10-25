//! `Approx_ComputeLine::Perform` (`gxx:832-1219`): the non-cutting and
//! cutting loop that splits the line until every piece is approximated within
//! tolerance.
//!
//! Translation notes:
//! - `LineTool::MakeMLOneMorePoint` always returns `false` in the `AppDef`
//!   instantiation (`AppDef_MyLineTool.cxx:164-171`), so every
//!   `isOtherLineMade` branch is kept for control-flow parity but its body is
//!   unreachable; the out-parameter it never writes is `None` here.
//! - `LineTool::MakeMLBetween` returns its input (`cxx:157-162`), and
//!   `WhatStatus` returns `Approx_NoPointsAdded` (`cxx:173-176`); both are used
//!   verbatim so the `Approx_PointsAdded` arm stays as OCCT wrote it.

use occt_core::approx::{ApproxParametrizationType, ApproxStatus};
use occt_math::MathVector;

use crate::approx::line_tool::LineTool;
use crate::app_par_curves::constraint::ConstraintCouple;
use crate::app_par_curves::multi_curve::MultiCurve;
use crate::approx::compute_line::check_multi_curve::check_multi_curve;
use crate::approx::compute_line::ComputeLine;

impl ComputeLine {
    /// `Perform(Line)` (`gxx:832-1219`).
    ///
    /// `Finish` is set and then the routine returns on the same line in OCCT
    /// (`gxx:933-937`, `:949-953`, `:1055-1058`, `:1161-1165`), which Rust
    /// flags as a dead assignment; the assignments are kept for control-flow
    /// parity.
    #[allow(unused_assignments)]
    pub fn perform<L: LineTool>(&mut self, line: &L) {
        if !self.my_is_clear {
            self.my_multi_curves.clear();
            self.my_par.clear();
            self.tol_ers3d.clear();
            self.tol_ers2d.clear();
            self.my_multi_line_nb = 0;
        } else {
            self.my_is_clear = false;
        }

        let mut finish = false;
        let mut begin = true;
        let mut ok = false;
        let mut go_up = false;
        // `thetol3d` / `thetol2d` are `Compute`'s out-parameters (`gxx:886`,
        // `:1209`); `Perform` never reads them back, so they stay at
        // `Compute`'s discard value here.
        let thetol3d = 0.0f64;
        let thetol2d = 0.0f64;

        let thefirstpt = line.first_point();
        let thelastpt = line.last_point();
        let myfirstpt_init = thefirstpt;
        let mut myfirstpt = myfirstpt_init;
        let mut mylastpt = thelastpt;

        let mut my_couple1 = ConstraintCouple::from_parts(myfirstpt, self.my_first_c);
        let mut my_couple2 = ConstraintCouple::from_parts(mylastpt, self.my_last_c);
        self.my_constraints[0] = my_couple1.clone();
        self.my_constraints[1] = my_couple2.clone();

        let mut the_param = MathVector::new(thefirstpt as usize, thelastpt as usize);
        let mut oldlastpt;

        if !self.my_cut {
            if let Some(first_param) = self.my_first_param.clone() {
                self.load_initial_parameters(&mut the_param, &first_param, thefirstpt);
            } else {
                the_param = self.compute_parameters(line, thefirstpt, thelastpt);
            }
            self.the_multi_curve = MultiCurve::new();
            let (done, indbad) = self.compute(line, myfirstpt, mylastpt, &mut the_param);
            self.alldone = done;
            let _ = (thetol3d, thetol2d);
            // `gxx:886-896`: a rejected point (`indbad != 0`) is first offered
            // to `MakeMLOneMorePoint`, which may build a denser line to
            // approximate instead.
            let an_other_line0 = if indbad != 0 {
                line.make_ml_one_more_point(myfirstpt, mylastpt, indbad)
            } else {
                None
            };
            if let Some(an_other_line0) = an_other_line0 {
                self.my_is_clear = true;
                // `++myMultiLineNb;` is commented out in OCCT (`gxx:894`).
                self.perform(&an_other_line0);
                self.alldone = true;
            }
            if !self.alldone && self.the_multi_curve.nb_curves() > 0 {
                let mylen = mylastpt - myfirstpt + 1;
                self.push_current(myfirstpt, mylen);
            }
        } else {
            loop {
                if finish {
                    break;
                }
                oldlastpt = mylastpt;
                // Management of the multiline splitting for approximation.
                if !begin {
                    if !go_up {
                        if ok {
                            // Calculation of the part to approximate.
                            myfirstpt = mylastpt;
                            mylastpt = thelastpt;
                            if myfirstpt == thelastpt {
                                finish = true;
                                self.alldone = true;
                                return;
                            }
                        } else {
                            let nbp = mylastpt - myfirstpt + 1;
                            let my_status = line.what_status(myfirstpt, mylastpt);
                            if my_status == ApproxStatus::NoPointsAdded && nbp <= self.my_degremax + 1
                            {
                                let interpol = self.compute_curve(line, myfirstpt, mylastpt);
                                if interpol && mylastpt == thelastpt {
                                    finish = true;
                                    self.alldone = true;
                                    return;
                                }
                            }
                            mylastpt = (myfirstpt + mylastpt) / 2;
                        }
                    }
                    go_up = false;
                }

                // Verification of the number of remaining points against the
                // requested degree.
                let nbp = mylastpt - myfirstpt + 1;
                let my_status = line.what_status(myfirstpt, mylastpt);
                if nbp <= self.my_degremax + 5 {
                    go_up = false;
                    ok = true;
                    if my_status == ApproxStatus::PointsAdded {
                        // Recursive splitting.
                        go_up = true;
                        let an_other_line1 = line.make_ml_between(
                            myfirstpt,
                            mylastpt,
                            nbp - 1,
                        );
                        let nbpdsotherligne = an_other_line1.first_point()
                            - an_other_line1.last_point();

                        // If MakeML failed, return an empty line.
                        if nbpdsotherligne == 0 || self.my_multi_line_nb >= 3 {
                            if myfirstpt == mylastpt {
                                break;
                            }
                            my_couple1.set_index(myfirstpt);
                            my_couple2.set_index(mylastpt);
                            self.my_constraints[0] = my_couple1.clone();
                            self.my_constraints[1] = my_couple2.clone();

                            let save_par = self.par;
                            self.par = ApproxParametrizationType::IsoParametric;
                            let mut param = self.compute_parameters(line, myfirstpt, mylastpt);
                            self.the_multi_curve = MultiCurve::new();
                            let (done, indbad) =
                                self.compute(line, myfirstpt, mylastpt, &mut param);
                            ok = done;
                            // `gxx:1002-1013`.
                            let an_other_line2 = if indbad != 0 {
                                line.make_ml_one_more_point(myfirstpt, mylastpt, indbad)
                            } else {
                                None
                            };
                            if let Some(an_other_line2) = an_other_line2 {
                                self.my_is_clear = true;
                                // `++myMultiLineNb;` is commented out in OCCT
                                // (`gxx:1011`).
                                self.par = save_par;
                                self.perform(&an_other_line2);
                                ok = true;
                            }
                            if !ok {
                                let tt3d = self.current_tol3d;
                                let tt2d = self.current_tol2d;
                                let save_parameters = self.my_parameters.clone();
                                let save_multi_curve = self.the_multi_curve.clone();

                                if save_par != ApproxParametrizationType::IsoParametric {
                                    self.par = save_par;
                                } else {
                                    self.par = ApproxParametrizationType::ChordLength;
                                }

                                param = self.compute_parameters(line, myfirstpt, mylastpt);
                                let (done2, indbad2) =
                                    self.compute(line, myfirstpt, mylastpt, &mut param);
                                ok = done2;
                                // `gxx:1029-1043`: second `MakeMLOneMorePoint`
                                // attempt, this time on the chord-length
                                // parametrisation.
                                let an_other_line2 = if indbad2 != 0 {
                                    line.make_ml_one_more_point(myfirstpt, mylastpt, indbad2)
                                } else {
                                    None
                                };
                                if let Some(an_other_line2) = an_other_line2 {
                                    self.my_is_clear = true;
                                    // `++myMultiLineNb;` is commented out in
                                    // OCCT (`gxx:1040`).
                                    self.perform(&an_other_line2);
                                    ok = true;
                                }

                                if !ok
                                    && tt3d <= self.current_tol3d
                                    && tt2d <= self.current_tol2d
                                {
                                    self.current_tol3d = tt3d;
                                    self.current_tol2d = tt2d;
                                    self.my_parameters = save_parameters;
                                    self.the_multi_curve = save_multi_curve;
                                }
                            }
                            self.par = save_par;
                            if myfirstpt == thelastpt {
                                finish = true;
                                self.alldone = true;
                                return;
                            }

                            oldlastpt = mylastpt;
                            if !ok {
                                self.tolreached = false;
                                if self.the_multi_curve.nb_curves() == 0 {
                                    self.my_multi_curves.clear();
                                    return;
                                }
                                let mut indbad2 = 0;
                                let checked = check_multi_curve(
                                    &self.the_multi_curve.clone(),
                                    line,
                                    myfirstpt,
                                    mylastpt,
                                    &mut indbad2,
                                );
                                // `gxx:1074-1103`: a curve rejected by
                                // `CheckMultiCurve` is retried on a denser line.
                                let an_other_line3 = if !checked {
                                    line.make_ml_one_more_point(myfirstpt, mylastpt, indbad2)
                                } else {
                                    None
                                };
                                if let Some(an_other_line3) = an_other_line3 {
                                    self.my_is_clear = true;
                                    // `++myMultiLineNb;` is commented out in
                                    // OCCT (`gxx:1084`).
                                    self.perform(&an_other_line3);
                                    myfirstpt = mylastpt;
                                    mylastpt = thelastpt;
                                } else {
                                    let mylen = oldlastpt - myfirstpt + 1;
                                    self.push_current(myfirstpt, mylen);
                                }
                            }
                            myfirstpt = oldlastpt;
                            mylastpt = thelastpt;
                        } else {
                            self.my_is_clear = true;
                            self.my_multi_line_nb += 1;
                            self.perform(&an_other_line1);
                            myfirstpt = mylastpt;
                            mylastpt = thelastpt;
                        }
                    }

                    if my_status == ApproxStatus::NoPointsAdded && !begin {
                        // Return the best approximation obtained so far.
                        go_up = true;
                        self.tolreached = false;
                        if self.the_multi_curve.nb_curves() == 0 {
                            self.my_multi_curves.clear();
                            return;
                        }
                        let mylen = oldlastpt - myfirstpt + 1;
                        self.push_current(myfirstpt, mylen);

                        myfirstpt = oldlastpt;
                        mylastpt = thelastpt;
                    } else if my_status == ApproxStatus::NoApproximation {
                        // No approximation is done between myfirstpt and mylastpt.
                        go_up = true;
                        myfirstpt = mylastpt;
                        mylastpt = thelastpt;
                    }
                }

                if myfirstpt == thelastpt {
                    finish = true;
                    self.alldone = true;
                    return;
                }
                if !go_up {
                    if myfirstpt == mylastpt {
                        break;
                    }
                    my_couple1.set_index(myfirstpt);
                    my_couple2.set_index(mylastpt);
                    self.my_constraints[0] = my_couple1.clone();
                    self.my_constraints[1] = my_couple2.clone();

                    // Calculation of parameters on this new interval.
                    let mut param = MathVector::new(myfirstpt as usize, mylastpt as usize);
                    if begin {
                        if let Some(first_param) = self.my_first_param.clone() {
                            for i in first_param.0..=first_param.0 + first_param.1.len() as i32 - 1 {
                                param.set_value(i as usize, first_param.1[(i - first_param.0) as usize]);
                            }
                            self.my_first_param = None;
                        } else {
                            param = self.compute_parameters(line, myfirstpt, mylastpt);
                        }
                        the_param = param.clone();
                        begin = false;
                    } else {
                        let pfirst = the_param.value(myfirstpt as usize);
                        let plast = the_param.value(mylastpt as usize);
                        for i in myfirstpt..=mylastpt {
                            let v = the_param.value(i as usize);
                            param.set_value(i as usize, (v - pfirst) / (plast - pfirst));
                        }
                    }

                    self.the_multi_curve = MultiCurve::new();
                    let (done, _indbad) = self.compute(line, myfirstpt, mylastpt, &mut param);
                    ok = done;
                    if myfirstpt == thelastpt {
                        finish = true;
                        self.alldone = true;
                        return;
                    }
                }
            }
        }
        let _ = (thetol3d, thetol2d);
    }

    /// `myfirstParam` copy-in (`gxx:1186-1192`, `:871-879`).
    fn load_initial_parameters(
        &self,
        the_param: &mut MathVector,
        first_param: &(i32, Vec<f64>),
        thefirstpt: i32,
    ) {
        let (lo, vals) = first_param;
        for i in *lo..=*lo + vals.len() as i32 - 1 {
            the_param.set_value((i + thefirstpt - 1) as usize, vals[(i - lo) as usize]);
        }
    }

    /// Shared `myMultiCurves.Append(TheMultiCurve); Tolers3d.Append(...);
    /// myPar.Append(ThePar)` block (`gxx:904-914`, `:1092-1103`,
    /// `:1135-1143`).
    fn push_current(&mut self, myfirstpt: i32, mylen: i32) {
        self.my_multi_curves.push(self.the_multi_curve.clone());
        self.tol_ers3d.push(self.current_tol3d);
        self.tol_ers2d.push(self.current_tol2d);
        let my_par_len = self
            .my_parameters
            .as_ref()
            .map(|p| p.1.len() as i32)
            .unwrap_or(0);
        let a_len = if my_par_len > mylen { my_par_len } else { mylen };
        let (_, vals) = match self.my_parameters.clone() {
            Some(v) => v,
            None => (1, Vec::new()),
        };
        let mut the_par = vec![0.0f64; a_len.max(0) as usize];
        for i in 0..a_len {
            // `myParameters->Value(myParameters->Lower() + i)`
            let src = i as usize;
            the_par[i as usize] = if src < vals.len() { vals[src] } else { 0.0 };
        }
        self.my_par.push((myfirstpt, the_par));
    }
}






