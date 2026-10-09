//! Port of `IntCurve_ExactIntersectionPoint` instantiated as
//! `Geom2dInt_ExactIntersectionPointOfTheIntPCurvePCurveOfGInter`
//! (`IntCurve_ExactIntersectionPoint.gxx:26-271`).
//!
//! Newton refinement of a polygon-polygon interference point: it turns the
//! segment/parameter hint from `Intf_SectionPoint` into a bracket on both
//! curves, solves `C1(u) - C2(v) = 0` with `math_FunctionSetRoot`, and widens
//! the bracket segment by segment while no root is found.

use occt_core::intf::IntfPolygon2d;
use occt_math::{MathFunctionSetRoot, MathFunctionSetWithDerivatives, MathVector};

use crate::curve::Curve2d;
use super::curve_tool;
use super::dist_between_pcurves::DistBetweenPCurves;
use super::polygon2d::Geom2dIntPolygon2d;

/// Max iterations of the inner `math_FunctionSetRoot` (`gxx:244`).
const NB_ITERATIONS: i32 = 60;

/// `Geom2dInt_ExactIntersectionPointOfTheIntPCurvePCurveOfGInter`.
pub struct Geom2dIntExactIntersectionPoint<'a> {
    /// `done` (`hxx`).
    done: bool,
    /// `nbroots` (`hxx`).
    nbroots: i32,
    /// `myTol` = `Tol * Tol` (`gxx:30`).
    my_tol: f64,
    /// `FctDist` (`hxx`).
    fct_dist: DistBetweenPCurves<'a>,
    /// `ToleranceVector` (`gxx:32`).
    tolerance_vector: MathVector,
    /// `BInfVector` (`gxx:33`).
    b_inf_vector: MathVector,
    /// `BSupVector` (`gxx:34`).
    b_sup_vector: MathVector,
    /// `StartingPoint` (`gxx:35`).
    starting_point: MathVector,
    /// `Root` (`gxx:36`).
    root: MathVector,
    /// `anErrorOccurred` (`gxx:37`).
    an_error_occurred: bool,
}

impl<'a> Geom2dIntExactIntersectionPoint<'a> {
    /// `IntCurve_ExactIntersectionPoint(C1, C2, Tol)` (`gxx:26-41`).
    pub fn new(c1: &'a dyn Curve2d, c2: &'a dyn Curve2d, tol: f64) -> Self {
        let mut tolerance_vector = MathVector::new(1, 2);
        tolerance_vector.set_value(1, curve_tool::eps_x());
        let mut r = Self {
            done: false,
            nbroots: 0,
            my_tol: tol * tol,
            fct_dist: DistBetweenPCurves::new(c1, c2),
            tolerance_vector,
            b_inf_vector: MathVector::new(1, 2),
            b_sup_vector: MathVector::new(1, 2),
            starting_point: MathVector::new(1, 2),
            root: MathVector::new(1, 2),
            an_error_occurred: false,
        };
        r.tolerance_vector.set_value(2, curve_tool::eps_x());
        r
    }

    /// `Perform(Poly1, Poly2, NumSegOn1, NumSegOn2, ParamOnSeg1, ParamOnSeg2)`
    /// (`gxx:45-200`).
    pub fn perform(
        &mut self,
        poly1: &Geom2dIntPolygon2d,
        poly2: &Geom2dIntPolygon2d,
        num_seg_on1: &mut i32,
        num_seg_on2: &mut i32,
        param_on_seg1: &mut f64,
        param_on_seg2: &mut f64,
    ) {
        //----------------------------------------------------------------------
        //-- On prend comme bornes de recherches  :
        //--
        //--   Segment      :      i-1        i           i+1        i+2
        //--
        //--                  |---------|-----X-------|---------|----------|
        //--                Inf                                Sup
        //----------------------------------------------------------------------
        if *num_seg_on1 >= poly1.nb_segments() as i32 && *param_on_seg1 == 0.0 {
            *num_seg_on1 -= 1;
            *param_on_seg1 = 1.0;
        }
        if *num_seg_on2 >= poly2.nb_segments() as i32 && *param_on_seg2 == 0.0 {
            *num_seg_on2 -= 1;
            *param_on_seg2 = 1.0;
        }
        if *num_seg_on1 <= 0 {
            *num_seg_on1 = 1;
            *param_on_seg1 = 0.0;
        }
        if *num_seg_on2 <= 0 {
            *num_seg_on2 = 1;
            *param_on_seg2 = 0.0;
        }

        self.starting_point
            .set_value(1, poly1.approx_param_on_curve(*num_seg_on1 as usize, *param_on_seg1));
        if *num_seg_on1 <= 2 {
            self.b_inf_vector.set_value(1, poly1.inf_parameter());
        } else {
            self.b_inf_vector
                .set_value(1, poly1.approx_param_on_curve((*num_seg_on1 - 1) as usize, 0.0));
        }
        if *num_seg_on1 >= (poly1.nb_segments() as i32 - 2) {
            self.b_sup_vector.set_value(1, poly1.sup_parameter());
        } else {
            self.b_sup_vector
                .set_value(1, poly1.approx_param_on_curve((*num_seg_on1 + 2) as usize, 0.0));
        }

        self.starting_point
            .set_value(2, poly2.approx_param_on_curve(*num_seg_on2 as usize, *param_on_seg2));
        if *num_seg_on2 <= 2 {
            self.b_inf_vector.set_value(2, poly2.inf_parameter());
        } else {
            self.b_inf_vector
                .set_value(2, poly2.approx_param_on_curve((*num_seg_on2 - 1) as usize, 0.0));
        }
        if *num_seg_on2 >= (poly2.nb_segments() as i32 - 2) {
            self.b_sup_vector.set_value(2, poly2.sup_parameter());
        } else {
            self.b_sup_vector
                .set_value(2, poly2.approx_param_on_curve((*num_seg_on2 + 2) as usize, 0.0));
        }

        self.math_perform();
        if self.nbroots == 0 {
            //-- On risque de donner des bornes sur la courbe 1 trop etroites.
            let mut diff = 1i32;
            let an_binf_vector = self.b_inf_vector.value(1);
            let an_bsup_vector = self.b_sup_vector.value(1);
            //---------------- On elargit les bornes par la gauche --------------------
            loop {
                diff += 1;
                if (*num_seg_on1 - diff) <= 1 {
                    self.b_inf_vector.set_value(1, poly1.inf_parameter());
                    diff = 0;
                } else {
                    self.b_inf_vector
                        .set_value(1, poly1.approx_param_on_curve((*num_seg_on1 - diff) as usize, 0.0));
                }
                self.math_perform();
                //-- le 18 nov 97
                if diff > 3 {
                    diff += *num_seg_on1 / 2;
                }
                if !(self.nbroots == 0 && diff != 0) {
                    break;
                }
            }
            //---------------- On elargit les bornes par la droite --------------------
            if self.nbroots == 0 {
                self.b_inf_vector.set_value(1, an_binf_vector);
                diff = 1;
                loop {
                    diff += 1;
                    if (*num_seg_on1 + diff) >= (poly1.nb_segments() as i32 - 1) {
                        self.b_sup_vector.set_value(1, poly1.sup_parameter());
                        diff = 0;
                    } else {
                        self.b_sup_vector.set_value(
                            1,
                            poly1.approx_param_on_curve((*num_seg_on1 + 1 + diff) as usize, 0.0),
                        );
                    }
                    self.math_perform();
                    //-- le 18 nov 97
                    if diff > 3 {
                        diff += 1 + (poly1.nb_segments() as i32 - *num_seg_on1) / 2;
                    }
                    if !(self.nbroots == 0 && diff != 0) {
                        break;
                    }
                }
            }
            self.b_sup_vector.set_value(1, an_bsup_vector);

            if self.nbroots == 0 {
                //-- On risque de donner des bornes sur la courbe 1 trop etroites.
                let mut diff = 1i32;
                let an_binf_vector = self.b_inf_vector.value(2);
                let an_bsup_vector = self.b_sup_vector.value(2);
                //---------------- On elargit les bornes par la gauche --------------------
                loop {
                    diff += 1;
                    if (*num_seg_on2 - diff) <= 1 {
                        self.b_inf_vector.set_value(2, poly2.inf_parameter());
                        diff = 0;
                    } else {
                        self.b_inf_vector
                            .set_value(2, poly2.approx_param_on_curve((*num_seg_on2 - diff) as usize, 0.0));
                    }
                    self.math_perform();
                    //-- le 18 nov 97
                    if diff > 3 {
                        diff += *num_seg_on2 / 2;
                    }
                    if !(self.nbroots == 0 && diff != 0) {
                        break;
                    }
                }
                //---------------- On elargit les bornes par la droite --------------------
                if self.nbroots == 0 {
                    self.b_inf_vector.set_value(2, an_binf_vector);
                    diff = 1;
                    loop {
                        diff += 1;
                        if (*num_seg_on2 + diff) >= (poly2.nb_segments() as i32 - 1) {
                            self.b_sup_vector.set_value(2, poly2.sup_parameter());
                            diff = 0;
                        } else {
                            self.b_sup_vector.set_value(
                                2,
                                poly2.approx_param_on_curve((*num_seg_on2 + 1 + diff) as usize, 0.0),
                            );
                        }
                        self.math_perform();
                        //-- le 18 nov 97
                        if diff > 3 {
                            diff += 1 + (poly2.nb_segments() as i32 - *num_seg_on2) / 2;
                        }
                        if !(self.nbroots == 0 && diff != 0) {
                            break;
                        }
                    }
                }
                self.b_sup_vector.set_value(2, an_bsup_vector);
            }
        }
    }

    /// `Perform(Uo, Vo, UInf, VInf, USup, VSup)` (`gxx:205-221`).
    pub fn perform_bounds(&mut self, uo: f64, vo: f64, u_inf: f64, v_inf: f64, u_sup: f64, v_sup: f64) {
        self.done = true;
        self.b_inf_vector.set_value(1, u_inf);
        self.b_inf_vector.set_value(2, v_inf);
        self.b_sup_vector.set_value(1, u_sup);
        self.b_sup_vector.set_value(2, v_sup);
        self.starting_point.set_value(1, uo);
        self.starting_point.set_value(2, vo);

        self.math_perform();
    }

    /// `NbRoots()` (`gxx:227-230`).
    pub fn nb_roots(&self) -> i32 {
        self.nbroots
    }

    /// `Roots(U, V)` (`gxx:234-238`).
    pub fn roots(&self) -> (f64, f64) {
        (self.root.value(1), self.root.value(2))
    }

    /// `MathPerform()` (`gxx:242-263`).
    fn math_perform(&mut self) {
        let mut dist = DistBetweenPCurves::new(self.fct_dist.curve1, self.fct_dist.curve2);
        let mut fct = MathFunctionSetRoot::new(&self.fct_dist, &self.tolerance_vector, NB_ITERATIONS);
        fct.perform_with_bounds(
            &mut dist,
            &self.starting_point,
            &self.b_inf_vector,
            &self.b_sup_vector,
            false,
        );

        if fct.is_done() {
            self.root = fct.root();
            self.nbroots = 1;
            let mut xy = MathVector::new(1, 2);
            dist.value(&self.root, &mut xy);
            let dist2 = xy.value(1) * xy.value(1) + xy.value(2) * xy.value(2);
            if dist2 > self.my_tol {
                self.nbroots = 0;
            }
        } else {
            self.an_error_occurred = true;
            self.nbroots = 0;
        }
    }

    /// `AnErrorOccurred()` (`gxx:267-271`).
    pub fn an_error_occurred(&self) -> bool {
        self.an_error_occurred
    }

    /// `done` (`gxx:208`), only ever set by `Perform(Uo, Vo, ...)`. Kept for
    /// transcription fidelity with the OCCT member.
    pub fn is_done(&self) -> bool {
        self.done
    }
}
