//! `FairCurve_Batten` (`FairCurve_Batten.hxx:37-270`, `.cxx:35-632`,
//! `.lxx:27-113`): the batten (physical spline) fair-curve solver.
//!
//! Every entry point follows the `.cxx` body line by line; the OCCT
//! exceptions (`Standard_NullValue`, `Standard_NegativeValue`,
//! `Standard_ConstructionError`) become `Err(&'static str)` values or `None`
//! results, never a silent clamp.

use std::cell::RefCell;
use std::rc::Rc;

use occt_core::bspl::banded_interp::knot_sequence;
use occt_core::bspl::plib::coefficients_poles_2d;
use occt_core::gp::{GpPnt2d, GpVec2d, GpXY};
use occt_core::precision;
use occt_math::{hermite_coefficients, MathMultipleVarFunction, MathVector};

use crate::bspline_curve::Geom2dBSplineCurve;

use super::analysis_code::FairCurveAnalysisCode;
use super::batten_law::BattenLaw;
use super::energy_of_batten::EnergyOfBatten;
use super::newton::FairCurveNewton;

/// `FairCurve_Batten` (`FairCurve_Batten.hxx:37-270`).
#[derive(Debug, Clone)]
pub struct Batten {
    /// `myCode`.
    my_code: FairCurveAnalysisCode,
    /// `OldP1`.
    old_p1: GpPnt2d,
    /// `OldP2`.
    old_p2: GpPnt2d,
    /// `OldAngle1`.
    old_angle1: f64,
    /// `OldAngle2`.
    old_angle2: f64,
    /// `OldHeight`.
    old_height: f64,
    /// `OldSlope`.
    old_slope: f64,
    /// `OldSlidingFactor`.
    old_sliding_factor: f64,
    /// `OldFreeSliding`.
    old_free_sliding: bool,
    /// `OldConstraintOrder1`.
    old_constraint_order1: i32,
    /// `OldConstraintOrder2`.
    old_constraint_order2: i32,
    /// `NewP1`.
    new_p1: GpPnt2d,
    /// `NewP2`.
    new_p2: GpPnt2d,
    /// `NewAngle1`.
    new_angle1: f64,
    /// `NewAngle2`.
    new_angle2: f64,
    /// `NewHeight`.
    new_height: f64,
    /// `NewSlope`.
    new_slope: f64,
    /// `NewSlidingFactor`.
    new_sliding_factor: f64,
    /// `NewFreeSliding`.
    new_free_sliding: bool,
    /// `NewConstraintOrder1`.
    new_constraint_order1: i32,
    /// `NewConstraintOrder2`.
    new_constraint_order2: i32,
    /// `Degree`.
    degree: usize,
    /// `Poles`.
    poles: Vec<GpPnt2d>,
    /// `Flatknots`: the `BSplCLib::KnotSequence` of `Knots`/`Mults`/`Degree`.
    flatknots: Vec<f64>,
    /// `Knots`: the distinct knots.
    knots: Vec<f64>,
    /// `Mults`: the distinct multiplicities.
    mults: Vec<i32>,
}

impl Batten {
    /// `FairCurve_Batten(P1, P2, Height, Slope = 0)` (`.cxx:35-122`): the
    /// initial batten is the straight segment `P1P2` raised from degree 1 to
    /// `Degree = 9`.
    pub fn new(
        p1: GpPnt2d,
        p2: GpPnt2d,
        height: f64,
        slope: f64,
    ) -> Result<Self, &'static str> {
        // `.cxx:62-72`: the two throws.
        if p1.is_equal(&p2, precision::CONFUSION) {
            return Err("FairCurve : P1 and P2 are confused");
        }
        if height <= 0.0 {
            return Err("FairCurve : Height is not positive");
        }
        let degree = 9usize;
        // `.cxx:77-89`: `Knots = (0, 1)`, `Mults = (2, 2)`, `Poles = (P1, P2)`.
        let mut curve = Geom2dBSplineCurve::new(
            vec![p1.x(), p2.x()],
            vec![p1.y(), p2.y()],
            vec![0.0, 0.0, 1.0, 1.0],
            1,
        )?;
        // `.cxx:92-108`: `BSplCLib::IncreaseDegree(1, Degree, false, ...)`.
        curve.increase_degree(degree)?;
        let (knots, mults) = curve.distinct_knots_and_mults();
        Ok(Self {
            my_code: FairCurveAnalysisCode::Ok,
            old_p1: p1,
            old_p2: p2,
            old_angle1: 0.0,
            old_angle2: 0.0,
            old_height: height,
            old_slope: slope,
            old_sliding_factor: 1.0,
            old_free_sliding: false,
            old_constraint_order1: 1,
            old_constraint_order2: 1,
            new_p1: p1,
            new_p2: p2,
            new_angle1: 0.0,
            new_angle2: 0.0,
            new_height: height,
            new_slope: slope,
            new_sliding_factor: 1.0,
            new_free_sliding: false,
            new_constraint_order1: 1,
            new_constraint_order2: 1,
            degree,
            poles: curve
                .xs
                .iter()
                .zip(curve.ys.iter())
                .map(|(&x, &y)| GpPnt2d::new(x, y))
                .collect(),
            // `.cxx:115-122`: `Flatknots = BSplCLib::KnotSequence(Knots, Mults,
            // Degree, false)`.
            flatknots: curve.knots.clone(),
            knots,
            mults,
        })
    }

    /// `SetFreeSliding(FreeSliding)` (`.lxx:27-30`).
    pub fn set_free_sliding(&mut self, free_sliding: bool) {
        self.new_free_sliding = free_sliding;
    }

    /// `SetConstraintOrder1(ConstraintOrder)` (`.lxx:32-35`).
    pub fn set_constraint_order1(&mut self, constraint_order: i32) {
        self.new_constraint_order1 = constraint_order;
    }

    /// `SetConstraintOrder2(ConstraintOrder)` (`.lxx:37-40`).
    pub fn set_constraint_order2(&mut self, constraint_order: i32) {
        self.new_constraint_order2 = constraint_order;
    }

    /// `SetAngle1(Angle1)` (`.lxx:42-45`).
    pub fn set_angle1(&mut self, angle1: f64) {
        self.new_angle1 = angle1;
    }

    /// `SetAngle2(Angle2)` (`.lxx:47-50`).
    pub fn set_angle2(&mut self, angle2: f64) {
        self.new_angle2 = angle2;
    }

    /// `SetHeight(Height)` (`.lxx:52-57`).
    pub fn set_height(&mut self, height: f64) -> Result<(), &'static str> {
        if height <= 0.0 {
            return Err("FairCurve : Height is no positive");
        }
        self.new_height = height;
        Ok(())
    }

    /// `SetSlope(Slope)` (`.lxx:59-62`).
    pub fn set_slope(&mut self, slope: f64) {
        self.new_slope = slope;
    }

    /// `SetSlidingFactor(SlidingFactor)` (`.lxx:64-67`).
    pub fn set_sliding_factor(&mut self, sliding_factor: f64) {
        self.new_sliding_factor = sliding_factor;
    }

    /// `GetFreeSliding()` (`.lxx:69-72`).
    pub fn get_free_sliding(&self) -> bool {
        self.old_free_sliding
    }

    /// `GetConstraintOrder1()` (`.lxx:74-77`).
    pub fn get_constraint_order1(&self) -> i32 {
        self.old_constraint_order1
    }

    /// `GetConstraintOrder2()` (`.lxx:79-82`).
    pub fn get_constraint_order2(&self) -> i32 {
        self.old_constraint_order2
    }

    /// `GetP1()` (`.lxx:84-87`).
    pub fn get_p1(&self) -> &GpPnt2d {
        &self.old_p1
    }

    /// `GetP2()` (`.lxx:89-92`).
    pub fn get_p2(&self) -> &GpPnt2d {
        &self.old_p2
    }

    /// `GetAngle1()` (`.lxx:94-97`).
    pub fn get_angle1(&self) -> f64 {
        self.old_angle1
    }

    /// `GetAngle2()` (`.lxx:99-102`).
    pub fn get_angle2(&self) -> f64 {
        self.old_angle2
    }

    /// `GetHeight()` (`.lxx:104-107`).
    pub fn get_height(&self) -> f64 {
        self.old_height
    }

    /// `GetSlope()` (`.lxx:109-112`).
    pub fn get_slope(&self) -> f64 {
        self.old_slope
    }

    /// `GetSlidingFactor()` (`.lxx:114-117`).
    pub fn get_sliding_factor(&self) -> f64 {
        self.old_sliding_factor
    }

    /// `Angles(P1, P2)` (`.cxx:130-137`): rotate the stored end tangents with
    /// the segment direction change.
    fn angles(&mut self, p1: &GpPnt2d, p2: &GpPnt2d) {
        let old = GpVec2d::from_xy(self.new_p2.coord.subtracted(&self.new_p1.coord));
        let new = GpVec2d::from_xy(p2.coord.subtracted(&p1.coord));
        let dangle = old.angle(&new);
        self.new_angle1 -= dangle;
        self.new_angle2 += dangle;
    }

    /// `SetP1(P1)` (`.cxx:139-149`).
    pub fn set_p1(&mut self, p1: GpPnt2d) -> Result<(), &'static str> {
        if p1.is_equal(&self.new_p2, precision::CONFUSION) {
            return Err("FairCurve : P1 and P2 are confused");
        }
        let p2 = self.new_p2;
        self.angles(&p1, &p2);
        self.new_p1 = p1;
        Ok(())
    }

    /// `SetP2(P2)` (`.cxx:151-161`).
    pub fn set_p2(&mut self, p2: GpPnt2d) -> Result<(), &'static str> {
        if self.new_p1.is_equal(&p2, precision::CONFUSION) {
            return Err("FairCurve : P1 and P2 are confused");
        }
        let p1 = self.new_p1;
        self.angles(&p1, &p2);
        self.new_p2 = p2;
        Ok(())
    }

    /// `Compute(Code, NbIterations = 50, Tolerance = 1.0e-3)`
    /// (`.cxx:163-246`): the homotopy loop, which splits the requested
    /// deformation into steps small enough for the Newton minimisation.
    ///
    /// Returns `(Ok, Code)` in place of OCCT's `bool` return plus the
    /// `FairCurve_AnalysisCode&` out-parameter.
    pub fn compute(
        &mut self,
        nb_iterations: i32,
        tolerance: f64,
    ) -> (bool, FairCurveAnalysisCode) {
        let mut ok = true;
        let mut end = false;
        // `.cxx:170-171`.
        let angle_max = 0.7;
        let angle_min = 2.0 * std::f64::consts::PI / 100.0;
        let mut acode = FairCurveAnalysisCode::Ok;

        while ok && !end {
            let mut d_angle1 = self.new_angle1 - self.old_angle1;
            let mut d_angle2 = self.new_angle2 - self.old_angle2;
            let mut ratio = 1.0f64;
            if self.new_constraint_order1 > 0 {
                let fraction = d_angle1.abs()
                    / (angle_max * (-self.old_angle1.abs() / angle_max).exp() + angle_min);
                if fraction > 1.0 {
                    ratio = 1.0 / fraction;
                }
            }
            if self.new_constraint_order2 > 0 {
                let fraction = d_angle2.abs()
                    / (angle_max * (-self.old_angle2.abs() / angle_max).exp() + angle_min);
                if fraction > 1.0 {
                    ratio = if ratio < 1.0 / fraction {
                        ratio
                    } else {
                        1.0 / fraction
                    };
                }
            }

            let old_dist = self.old_p1.distance(&self.old_p2);
            let new_dist = self.new_p1.distance(&self.new_p2);
            let fraction = (old_dist - new_dist).abs() / (old_dist / 3.0);
            if fraction > 1.0 {
                ratio = if ratio < 1.0 / fraction {
                    ratio
                } else {
                    1.0 / fraction
                };
            }

            // `.cxx:210`: `gp_Vec2d DeltaP1(OldP1, NewP1)`.
            let mut delta_p1 = GpVec2d::from_xy(self.new_p1.coord.subtracted(&self.old_p1.coord));
            let mut delta_p2 = GpVec2d::from_xy(self.new_p2.coord.subtracted(&self.old_p2.coord));
            let toler;
            if ratio == 1.0 {
                end = true;
                toler = tolerance;
            } else {
                delta_p1.multiply_scalar(ratio);
                delta_p2.multiply_scalar(ratio);
                d_angle1 *= ratio;
                d_angle2 *= ratio;
                toler = 10.0 * tolerance;
            }

            ok = self.compute_delta(
                &delta_p1,
                &delta_p2,
                d_angle1,
                d_angle2,
                &mut acode,
                nb_iterations,
                toler,
            );

            if acode != FairCurveAnalysisCode::Ok {
                end = true;
            }
            if self.new_free_sliding {
                self.new_sliding_factor = self.old_sliding_factor;
            }
            if self.new_constraint_order1 == 0 {
                self.new_angle1 = self.old_angle1;
            }
            if self.new_constraint_order2 == 0 {
                self.new_angle2 = self.old_angle2;
            }
        }
        self.my_code = acode;
        (ok, acode)
    }

    /// `Compute(DeltaP1, DeltaP2, DeltaAngle1, DeltaAngle2, Code, NbIterations,
    /// Tolerance)` (`.cxx:247-521`): one homotopy step, a Hermite
    /// interpolation of the constraint increments added to the current batten
    /// and then minimised by [`FairCurveNewton`].
    #[allow(clippy::too_many_arguments)]
    fn compute_delta(
        &mut self,
        delta_p1: &GpVec2d,
        delta_p2: &GpVec2d,
        delta_angle1: f64,
        delta_angle2: f64,
        acode: &mut FairCurveAnalysisCode,
        nb_iterations: i32,
        tolerance: f64,
    ) -> bool {
        let mut ok_compute = true;
        *acode = FairCurveAnalysisCode::Ok;

        // `.cxx:260-262`.
        let l = 2 + self.new_constraint_order1 + self.new_constraint_order2;
        let l_usize = l as usize;

        // `.cxx:269-276`: `PLib::HermiteCoefficients(0, 1, ContrOrder1,
        // ContrOrder2, HermiteCoef)`; `false` propagates out of `Compute`.
        let hermite = match hermite_coefficients(
            0.0,
            1.0,
            self.new_constraint_order1,
            self.new_constraint_order2,
        ) {
            Some(matrix) => matrix,
            None => return false,
        };

        // `.cxx:280-282`.
        let old_v = GpVec2d::from_xy(self.old_p2.coord.subtracted(&self.old_p1.coord));
        let new_v = GpVec2d::from_xy(
            self.old_p2
                .coord
                .added(&delta_p2.coord)
                .subtracted(&self.old_p1.coord.added(&delta_p1.coord)),
        );
        let d_angle_ref = new_v.angle(&old_v);

        // `.cxx:284-308`: the constraint increments, `ADelta(1..L)`.
        let mut a_delta = vec![GpXY::zero(); l_usize];
        a_delta[0] = delta_p1.coord;
        let mut kk = 1usize;
        if self.new_constraint_order1 > 0 {
            let mut old_derive =
                GpVec2d::from_xy(self.poles[1].coord.subtracted(&self.poles[0].coord));
            let knot_gap1 = self.knots[1] - self.knots[0];
            if knot_gap1.abs() > precision::COMPUTATIONAL {
                old_derive.multiply_scalar(self.degree as f64 / knot_gap1);
            }
            a_delta[kk] = old_derive
                .rotated(delta_angle1 - d_angle_ref)
                .subtracted(&old_derive)
                .coord;
            kk += 1;
        }
        a_delta[kk] = delta_p2.coord;
        kk += 1;
        if self.new_constraint_order2 > 0 {
            let n = self.poles.len();
            let mut old_derive =
                GpVec2d::from_xy(self.poles[n - 1].coord.subtracted(&self.poles[n - 2].coord));
            let nk = self.knots.len();
            let knot_gap2 = self.knots[nk - 1] - self.knots[nk - 2];
            if knot_gap2.abs() > precision::COMPUTATIONAL {
                old_derive.multiply_scalar(self.degree as f64 / knot_gap2);
            }
            a_delta[kk] = old_derive
                .rotated(d_angle_ref - delta_angle2)
                .subtracted(&old_derive)
                .coord;
        }

        // `.cxx:310-322`: `Interpolation(ii) = sum_kk HermiteCoef(kk, ii) *
        // ADelta(kk)`.
        let mut interpolation = Vec::with_capacity(l_usize);
        for ii in 0..l_usize {
            let mut aux_xy = GpXY::zero();
            for (kk, delta) in a_delta.iter().enumerate() {
                aux_xy.add(&delta.multiplied_scalar(hermite[kk][ii]));
            }
            interpolation.push(GpPnt2d::from_xy(aux_xy));
        }
        // `.cxx:322`: `PLib::CoefficientsPoles(Interpolation, NoWeights(),
        // HermitePoles, NoWeights())`.
        let hermite_poles = coefficients_poles_2d(&interpolation);

        // `.cxx:324-327`: two knots of multiplicity `L`, degree `L - 1`.
        let mut delta_curve = match Geom2dBSplineCurve::new(
            hermite_poles.iter().map(|p| p.x()).collect(),
            hermite_poles.iter().map(|p| p.y()).collect(),
            vec![0.0; l_usize]
                .into_iter()
                .chain(vec![1.0; l_usize])
                .collect(),
            l_usize - 1,
        ) {
            Ok(curve) => curve,
            Err(_) => return false,
        };

        // `.cxx:329-333`: `IncreaseDegree(Degree)` then `InsertKnots(Knots,
        // Mults, 1.e-10)`.
        if delta_curve.increase_degree(self.degree).is_err() {
            return false;
        }
        if self.mults.len() > 2
            && delta_curve
                .insert_knots(&self.knots, Some(&self.mults), 1.0e-10, false)
                .is_err()
        {
            return false;
        }

        // `.cxx:335-341`: `NPoles = DeltaCurve->Poles()`, then add the poles of
        // the current batten. `IncreaseDegree` + `InsertKnots` reproduce
        // `Knots`/`Mults`, so both arrays have the same structure.
        let mut npoles: Vec<GpPnt2d> = delta_curve
            .xs
            .iter()
            .zip(delta_curve.ys.iter())
            .map(|(&x, &y)| GpPnt2d::new(x, y))
            .collect();
        if npoles.len() != self.poles.len() {
            return false;
        }
        for (kk, pole) in npoles.iter_mut().enumerate() {
            pole.coord.add(&self.poles[kk].coord);
        }

        // `.cxx:343-354`.
        let alph1 = self.old_angle1 + delta_angle1;
        let alph2 = self.old_angle2 + delta_angle2;
        let last = npoles.len() - 1;
        let dist = npoles[last].distance(&npoles[0]);
        let l_reference = self.sliding_of_reference(dist, alph1, alph2);
        let ox = GpVec2d::new(1.0, 0.0);
        let p1p2 = GpVec2d::from_xy(npoles[last].coord.subtracted(&npoles[0].coord));
        let angle1 = ox.angle(&p1p2) + alph1;
        let angle2 = -ox.angle(&p1p2) + alph2;

        // `.cxx:356-374`: the sliding length, imposed or initial.
        let sliding_length = if !self.new_free_sliding {
            self.new_sliding_factor * l_reference
        } else if self.old_free_sliding {
            self.old_sliding_factor * l_reference
        } else {
            self.sliding_of_reference(dist, alph1, alph2)
        };

        // `.cxx:376-387`.
        let law = BattenLaw::new(self.new_height, self.new_slope, sliding_length);
        let poles_handle: Rc<RefCell<Vec<GpPnt2d>>> = Rc::new(RefCell::new(npoles));
        let ebatten = EnergyOfBatten::new(
            (self.degree + 1) as i32,
            self.flatknots.clone(),
            Rc::clone(&poles_handle),
            self.new_constraint_order1,
            self.new_constraint_order2,
            &law,
            sliding_length,
            self.new_free_sliding,
            angle1,
            angle2,
        );

        // `.cxx:388-391`.
        let mut v_init = MathVector::new(1, ebatten.nb_variables());
        let mut v_convex = 0.01 * (self.new_height / sliding_length).powi(3);
        if v_convex < 1.0e-12 {
            v_convex = 1.0e-12;
        }
        let _ok = ebatten.variable(&mut v_init);

        // `.cxx:397-403`.
        let mut newton = FairCurveNewton::new(
            ebatten.nb_variables(),
            tolerance * p1p2.magnitude() / 10.0,
            tolerance,
            nb_iterations.max(0) as usize,
            v_convex,
            true,
        );
        newton.perform(&ebatten, &v_init);
        let ok = newton.is_done();

        if ok {
            // `.cxx:407-408`: `Poles = NPoles`: the shared array already holds
            // the poles of the last `ComputePoles`.
            self.poles = poles_handle.borrow().clone();
            let location = newton.location().clone();
            if self.new_free_sliding {
                self.old_sliding_factor = location.value(location.upper()) / l_reference;
            } else {
                self.old_sliding_factor = self.new_sliding_factor;
            }
            if self.new_constraint_order1 == 0 {
                let tangente =
                    GpVec2d::from_xy(self.poles[1].coord.subtracted(&self.poles[0].coord));
                self.old_angle1 = p1p2.angle(&tangente);
            } else {
                self.old_angle1 = alph1;
            }
            if self.new_constraint_order2 == 0 {
                let n = self.poles.len();
                let tangente =
                    GpVec2d::from_xy(self.poles[n - 1].coord.subtracted(&self.poles[n - 2].coord));
                self.old_angle2 = tangente.reversed().angle(&p1p2.reversed());
            } else {
                self.old_angle2 = alph2;
            }
            self.old_p1 = self.poles[0];
            self.old_p2 = self.poles[last];
            self.old_constraint_order1 = self.new_constraint_order1;
            self.old_constraint_order2 = self.new_constraint_order2;
            self.old_free_sliding = self.new_free_sliding;
            self.old_slope = self.new_slope;
            self.old_height = self.new_height;
        } else {
            // `.cxx:449-462`.
            *acode = ebatten.status();
            if law.value(0.0).is_none() || law.value(1.0).is_none() {
                *acode = FairCurveAnalysisCode::NullHeight;
            } else {
                ok_compute = false;
            }
            return ok_compute;
        }

        // `.cxx:464-466`: `EBatten.Variable(VInit)` refills the *initial*
        // variables; the free-sliding check below reads the auxiliary value the
        // energy settled on.
        let _ok = ebatten.variable(&mut v_init);

        // `.cxx:468-475`: non-convergence and infinite sliding.
        if !newton.is_converged() {
            *acode = FairCurveAnalysisCode::NotConverged;
        }
        if self.new_free_sliding && v_init.value(v_init.upper()) > 2.0 * l_reference {
            *acode = FairCurveAnalysisCode::InfiniteSliding;
        }

        // `.cxx:477-488`: eventual insertion of nodes.
        let mut new_knots = false;
        let mut nb_knots = self.knots.len() as i32;
        let val_angles = self.old_angle1.abs()
            + self.old_angle2.abs()
            + 2.0 * (self.old_angle2 - self.old_angle1).abs();
        while val_angles > ((2 * (nb_knots - 2) + 1) * (1 + 2 * nb_knots)) as f64 {
            new_knots = true;
            nb_knots += nb_knots - 1;
        }

        if new_knots {
            // `.cxx:490-510`.
            let flat = knot_sequence(&self.knots, &self.mults, self.degree as i32);
            let mut new_bs = match Geom2dBSplineCurve::new(
                self.poles.iter().map(|p| p.x()).collect(),
                self.poles.iter().map(|p| p.y()).collect(),
                flat,
                self.degree,
            ) {
                Ok(curve) => curve,
                Err(_) => return false,
            };
            let n_mults = vec![(self.degree as i32) - 3; nb_knots as usize];
            let n_knots: Vec<f64> = (1..=nb_knots as usize)
                .map(|ii| (ii - 1) as f64 / (nb_knots - 1) as f64)
                .collect();
            if new_bs
                .insert_knots(&n_knots, Some(&n_mults), 1.0e-10, false)
                .is_err()
            {
                return false;
            }
            let (knots, mults) = new_bs.distinct_knots_and_mults();
            self.poles = new_bs
                .xs
                .iter()
                .zip(new_bs.ys.iter())
                .map(|(&x, &y)| GpPnt2d::new(x, y))
                .collect();
            self.mults = mults;
            self.knots = knots;
            // `.cxx:509`: `Flatknots = NewBS->KnotSequence()`.
            self.flatknots = new_bs.knots.clone();
        }

        ok_compute
    }

    /// `SlidingOfReference()` (`.cxx:523-527`).
    pub fn sliding_of_reference_default(&self) -> f64 {
        self.sliding_of_reference(
            self.new_p1.distance(&self.new_p2),
            self.new_angle1,
            self.new_angle2,
        )
    }

    /// `SlidingOfReference(Dist, Angle1, Angle2)` (`.cxx:529-585`): the
    /// reference length of the batten for the new constraints.
    pub fn sliding_of_reference(&self, dist: f64, angle1: f64, angle2: f64) -> f64 {
        if self.new_constraint_order1 == 0 && self.new_constraint_order2 == 0 {
            return dist;
        }
        // `.cxx:542-559`. The second arm pairs `NewAngle1` with `Angle1`
        // exactly as `.cxx:553` does.
        let a1 = if self.new_constraint_order1 == 0 {
            if self.new_angle2.abs() < std::f64::consts::PI {
                (angle2 / 2.0).abs()
            } else {
                std::f64::consts::FRAC_PI_2
            }
        } else {
            angle1.abs()
        };
        let a2 = if self.new_constraint_order2 == 0 {
            if self.new_angle1.abs() < std::f64::consts::PI {
                (angle1 / 2.0).abs()
            } else {
                std::f64::consts::FRAC_PI_2
            }
        } else {
            angle2.abs()
        };

        if angle1 * angle2 >= 0.0 {
            // `.cxx:560-564`: angles of the same sign.
            self.compute_law_pair(dist, a1, a2)
        } else {
            // `.cxx:566-583`: angles of opposite sign.
            let a_sum = a1 + a2;
            if a_sum.abs() < precision::ANGULAR {
                return self.compute_law_pair(dist, a1, std::f64::consts::FRAC_PI_2);
            }
            let ratio = a1 / a_sum;
            let mut angle_milieu = (1.0 - ratio).powi(2) * a1 + ratio.powi(2) * a2;
            if angle_milieu > std::f64::consts::FRAC_PI_2 {
                angle_milieu = std::f64::consts::FRAC_PI_2;
            }
            ratio * self.compute_law_pair(dist, a1, angle_milieu)
                + (1.0 - ratio) * self.compute_law_pair(dist, a2, angle_milieu)
        }
    }

    /// `Compute(Dist, Angle1, Angle2)` (`.cxx:587-599`).
    fn compute_law_pair(&self, dist: f64, angle1: f64, angle2: f64) -> f64 {
        let mut l1 = self.compute_law(dist, angle1);
        let mut l2 = self.compute_law(dist, angle2);
        if l1 < l2 {
            std::mem::swap(&mut l1, &mut l2);
        }
        0.3 * l1 + 0.7 * l2
    }

    /// `Compute(Dist, Angle)` (`.cxx:601-621`).
    fn compute_law(&self, dist: f64, angle: f64) -> f64 {
        if angle < precision::ANGULAR {
            return dist;
        }
        if angle < std::f64::consts::FRAC_PI_2 {
            return angle * dist / angle.sin();
        }
        if angle > std::f64::consts::PI {
            return (angle * std::f64::consts::PI).sqrt() * dist;
        }
        angle * dist
    }

    /// `Curve()` (`.cxx:623-630`): the computed batten as a 2D B-spline.
    /// `None` where the stored knot structure is inconsistent (OCCT would
    /// raise `Standard_ConstructionError`).
    pub fn curve(&self) -> Option<Geom2dBSplineCurve> {
        let flat = knot_sequence(&self.knots, &self.mults, self.degree as i32);
        Geom2dBSplineCurve::new(
            self.poles.iter().map(|p| p.x()).collect(),
            self.poles.iter().map(|p| p.y()).collect(),
            flat,
            self.degree,
        )
        .ok()
    }

    /// `myCode`: the code of the last `Compute`.
    pub fn code(&self) -> FairCurveAnalysisCode {
        self.my_code
    }
}

// `FairCurve_Batten::Dump(Standard_OStream&)` (`.cxx:632-691`) is a debug
// stream helper with no Rust counterpart: UNPORTED.
