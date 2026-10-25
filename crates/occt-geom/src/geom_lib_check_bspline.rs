//! `GeomLib_CheckBSplineCurve` / `GeomLib_Check2dBSplineCurve` — detect and
//! repair reversed end tangents of a B-spline.
//!
//! Sources:
//! - `ModelingData/TKGeomBase/GeomLib/GeomLib_CheckBSplineCurve.hxx:26-81`,
//!   `.cxx:24-188`
//! - `ModelingData/TKGeomBase/GeomLib/GeomLib_Check2dBSplineCurve.hxx:26-66`,
//!   `.cxx:26-187`
//!
//! Both classes scan the poles from both ends: from pole 1 they look for the
//! first pole `ii` whose direction (relative to pole 1) is not parallel to the
//! first tangent *and* points backwards (`tangent.Dot(a_vector) < 0`), then
//! redistribute the poles between as a straight line. The same happens
//! mirrored from the last pole.
//!
//! Translation notes:
//! - The 3d check compares `CrossProd.SquareMagnitude()` with
//!   `AngularTolerance * AngularTolerance` (`cxx:46`, `:64`), the 2d check
//!   compares `abs(CrossProd)` with `AngularTolerance` (`2d.cxx:46`, `:66`).
//!   They are not the same test and both are kept as they stand.
//! - `tangent_normalized` / `avector_normalized` are default-constructed
//!   `gp_Vec(2d)` (zero) and only assigned when the corresponding magnitude
//!   passes `myTolerance`; they are also read only under the same guard, so the
//!   zero value is never used. The port computes them inside the guard.
//! - `myIndSecondPole` / `myIndPrelastPole` keep their 1-based OCCT values.
//! - The OCCT class holds a `handle<Geom_BSplineCurve>` to the caller's curve
//!   and `FixTangent` writes through it. The port takes the curve at the call
//!   (`fix_tangent(&mut self, curve: &mut ...)`), which is the same aliasing.
//! - `FixedTangent` (`cxx:132-142`) returns a null handle when neither flag
//!   applies; `None` here.

use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_geom2d::bspline_curve::Geom2dBSplineCurve;

use crate::bspline_curve::GeomBSplineCurve;

/// `gp_Vec(P1, P2)` for 2d points.
fn vec2d_between(p1: &GpPnt2d, p2: &GpPnt2d) -> GpVec2d {
    GpVec2d::from_xy(p2.xy().subtracted(p1.xy()))
}

/// `class GeomLib_CheckBSplineCurve` (`hxx:26-81`).
#[derive(Debug, Clone)]
pub struct CheckBSplineCurve {
    /// `myDone`
    done: bool,
    /// `myFixFirstTangent`
    fix_first_tangent: bool,
    /// `myFixLastTangent`
    fix_last_tangent: bool,
    /// `myAngularTolerance` — `std::abs(AngularTolerance)` (`cxx:32`)
    angular_tolerance: f64,
    /// `myTolerance` — `std::abs(Tolerance)` (`cxx:33`)
    tolerance: f64,
    /// `myIndSecondPole`
    ind_second_pole: i32,
    /// `myIndPrelastPole`
    ind_prelast_pole: i32,
}

impl CheckBSplineCurve {
    /// `GeomLib_CheckBSplineCurve(Curve, Tolerance, AngularTolerance)`
    /// (`cxx:24-120`).
    pub fn new(curve: &GeomBSplineCurve, tolerance: f64, angular_tolerance: f64) -> Self {
        let mut this = Self {
            done: false,
            fix_first_tangent: false,
            fix_last_tangent: false,
            angular_tolerance: angular_tolerance.abs(),
            tolerance: tolerance.abs(),
            ind_second_pole: -1,
            ind_prelast_pole: -1,
        };

        let num_poles = curve.nb_poles();
        if !curve.periodic() && num_poles >= 4 {
            let cross_prod_sq_tol = this.angular_tolerance * this.angular_tolerance;

            // Near first (`cxx:48-81`).
            let mut tangent = GpVec::from_pnts(curve.pole(0), curve.pole(1));
            let tangent_magnitude = tangent.magnitude();
            let mut tangent_normalized = GpVec::new(0.0, 0.0, 0.0);
            if tangent_magnitude > this.tolerance {
                tangent_normalized = tangent.divided(tangent_magnitude);
            }

            let mut ii = 3usize;
            while ii <= num_poles {
                let a_vector = GpVec::from_pnts(curve.pole(0), curve.pole(ii - 1));
                let vector_magnitude = a_vector.magnitude();

                if tangent_magnitude > this.tolerance && vector_magnitude > this.tolerance {
                    let avector_normalized = a_vector.divided(vector_magnitude);

                    let cross_prod = tangent_normalized.crossed(&avector_normalized);
                    if cross_prod.square_magnitude() > cross_prod_sq_tol {
                        break;
                    }

                    let value = tangent.dot(&a_vector);
                    if value < 0.0 {
                        this.fix_first_tangent = true;
                        this.ind_second_pole = ii as i32;
                        break;
                    }
                }
                ii += 1;
            }

            // Near last (`cxx:84-112`).
            tangent = GpVec::from_pnts(curve.pole(num_poles - 1), curve.pole(num_poles - 2));
            let tangent_magnitude = tangent.magnitude();
            let mut tangent_normalized = GpVec::new(0.0, 0.0, 0.0);
            if tangent_magnitude > this.tolerance {
                tangent_normalized = tangent.divided(tangent_magnitude);
            }

            let mut ii = num_poles as i32 - 2;
            while ii >= 1 {
                let a_vector = GpVec::from_pnts(curve.pole(num_poles - 1), curve.pole((ii - 1) as usize));
                let vector_magnitude = a_vector.magnitude();

                if tangent_magnitude > this.tolerance && vector_magnitude > this.tolerance {
                    let avector_normalized = a_vector.divided(vector_magnitude);

                    let cross_prod = tangent_normalized.crossed(&avector_normalized);
                    if cross_prod.square_magnitude() > cross_prod_sq_tol {
                        break;
                    }

                    let value = tangent.dot(&a_vector);
                    if value < 0.0 {
                        this.fix_last_tangent = true;
                        this.ind_prelast_pole = ii;
                        break;
                    }
                }
                ii -= 1;
            }
        } else {
            this.done = true;
        }

        this
    }

    /// `IsDone()` (`hxx:38`, `.lxx`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `NeedTangentFix(FirstFlag, LastFlag)` (`cxx:125-128`).
    pub fn need_tangent_fix(&self) -> (bool, bool) {
        (self.fix_first_tangent, self.fix_last_tangent)
    }

    /// `FixedTangent(FirstFlag, LastFlag)` (`cxx:132-142`): a copy of the curve
    /// with the flagged ends repaired, or `None` when neither flag applies.
    ///
    /// `&mut self` because OCCT's `FixedTangent` calls `FixTangentOnCurve`,
    /// which sets `myDone = true` on the checker itself (`cxx:187`).
    pub fn fixed_tangent(
        &mut self,
        curve: &GeomBSplineCurve,
        first_flag: bool,
        last_flag: bool,
    ) -> Option<GeomBSplineCurve> {
        if !((self.fix_first_tangent && first_flag) || (self.fix_last_tangent && last_flag)) {
            return None;
        }
        let mut new_curve = curve.clone();
        self.fix_tangent_on_curve(&mut new_curve, first_flag, last_flag);
        Some(new_curve)
    }

    /// `FixTangent(FirstFlag, LastFlag)` (`cxx:147-150`): writes the caller's
    /// curve in place.
    pub fn fix_tangent(
        &mut self,
        curve: &mut GeomBSplineCurve,
        first_flag: bool,
        last_flag: bool,
    ) {
        self.fix_tangent_on_curve(curve, first_flag, last_flag);
    }

    /// `FixTangentOnCurve(theCurve, FirstFlag, LastFlag)` (`cxx:154-188`).
    ///
    /// `myDone = true` is unconditional at the end (`cxx:187`), even when no
    /// flag applied.
    fn fix_tangent_on_curve(
        &mut self,
        curve: &mut GeomBSplineCurve,
        first_flag: bool,
        last_flag: bool,
    ) {
        if self.fix_first_tangent && first_flag {
            let xyz1 = *curve.pole(0).xyz();
            let xyz2 = *curve.pole((self.ind_second_pole - 1) as usize).xyz();
            let nb_samples = f64::from(self.ind_second_pole - 1);
            let mut i = 2usize;
            while i < self.ind_second_pole as usize {
                let ii = f64::from(i as i32 - 1);
                let a_new_pole = GpPnt::from_xyz(&xyz1.multiply_scalar(1.0 - ii / nb_samples).added(
                    &xyz2.multiply_scalar(ii / nb_samples),
                ));
                curve.set_pole(i - 1, a_new_pole);
                i += 1;
            }
        }

        if self.fix_last_tangent && last_flag {
            let num_poles = curve.nb_poles();

            let xyz1 = *curve.pole(num_poles - 1).xyz();
            let xyz2 = *curve.pole((self.ind_prelast_pole - 1) as usize).xyz();
            let nb_samples = f64::from(num_poles as i32 - self.ind_prelast_pole);
            let mut i = num_poles as i32 - 1;
            while i > self.ind_prelast_pole {
                let ii = f64::from(num_poles as i32 - i);
                let a_new_pole = GpPnt::from_xyz(&xyz1.multiply_scalar(1.0 - ii / nb_samples).added(
                    &xyz2.multiply_scalar(ii / nb_samples),
                ));
                curve.set_pole(i as usize - 1, a_new_pole);
                i -= 1;
            }
        }

        self.done = true;
    }
}

/// `class GeomLib_Check2dBSplineCurve` (`2d.hxx:26-66`).
#[derive(Debug, Clone)]
pub struct Check2dBSplineCurve {
    /// `myDone`
    done: bool,
    /// `myFixFirstTangent`
    fix_first_tangent: bool,
    /// `myFixLastTangent`
    fix_last_tangent: bool,
    /// `myAngularTolerance` — `std::abs(AngularTolerance)` (`2d.cxx:32`)
    angular_tolerance: f64,
    /// `myTolerance` — `std::abs(Tolerance)` (`2d.cxx:33`)
    tolerance: f64,
    /// `myIndSecondPole`
    ind_second_pole: i32,
    /// `myIndPrelastPole`
    ind_prelast_pole: i32,
}

impl Check2dBSplineCurve {
    /// `GeomLib_Check2dBSplineCurve(Curve, Tolerance, AngularTolerance)`
    /// (`2d.cxx:26-120`).
    pub fn new(curve: &Geom2dBSplineCurve, tolerance: f64, angular_tolerance: f64) -> Self {
        let mut this = Self {
            done: false,
            fix_first_tangent: false,
            fix_last_tangent: false,
            angular_tolerance: angular_tolerance.abs(),
            tolerance: tolerance.abs(),
            ind_second_pole: -1,
            ind_prelast_pole: -1,
        };

        let num_poles = curve.nb_poles();
        if !curve.periodic && num_poles >= 4 {
            let cross_prod_tol = this.angular_tolerance;

            // Near first (`2d.cxx:50-80`).
            let mut tangent = vec2d_between(&pnt2d(curve, 1), &pnt2d(curve, 2));
            let tangent_magnitude = tangent.magnitude();
            let mut tangent_normalized = GpVec2d::from_xy(occt_core::gp::GpXY::zero());
            if tangent_magnitude > this.tolerance {
                tangent_normalized = tangent.divided(tangent_magnitude);
            }

            let mut ii = 3usize;
            while ii <= num_poles {
                let a_vector = vec2d_between(&pnt2d(curve, 1), &pnt2d(curve, ii));
                let vector_magnitude = a_vector.magnitude();

                if tangent_magnitude > this.tolerance && vector_magnitude > this.tolerance {
                    let avector_normalized = a_vector.divided(vector_magnitude);

                    let cross_prod = tangent_normalized.crossed(&avector_normalized);
                    if cross_prod.abs() > cross_prod_tol {
                        break;
                    }

                    let value = tangent.dot(&a_vector);
                    if value < 0.0 {
                        this.fix_first_tangent = true;
                        this.ind_second_pole = ii as i32;
                        break;
                    }
                }
                ii += 1;
            }

            // Near last (`2d.cxx:90-112`).
            tangent = vec2d_between(&pnt2d(curve, num_poles), &pnt2d(curve, num_poles - 1));
            let tangent_magnitude = tangent.magnitude();
            let mut tangent_normalized = GpVec2d::from_xy(occt_core::gp::GpXY::zero());
            if tangent_magnitude > this.tolerance {
                tangent_normalized = tangent.divided(tangent_magnitude);
            }

            let mut ii = num_poles as i32 - 2;
            while ii >= 1 {
                let a_vector = vec2d_between(&pnt2d(curve, num_poles), &pnt2d(curve, ii as usize));
                let vector_magnitude = a_vector.magnitude();

                if tangent_magnitude > this.tolerance && vector_magnitude > this.tolerance {
                    let avector_normalized = a_vector.divided(vector_magnitude);

                    let cross_prod = tangent_normalized.crossed(&avector_normalized);
                    if cross_prod.abs() > cross_prod_tol {
                        break;
                    }

                    let value = tangent.dot(&a_vector);
                    if value < 0.0 {
                        this.fix_last_tangent = true;
                        this.ind_prelast_pole = ii;
                        break;
                    }
                }
                ii -= 1;
            }
        } else {
            this.done = true;
        }

        this
    }

    /// `IsDone()` (`2d.hxx:37`, `.lxx`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `NeedTangentFix(FirstFlag, SecondFlag)` (`2d.cxx:125-128`).
    pub fn need_tangent_fix(&self) -> (bool, bool) {
        (self.fix_first_tangent, self.fix_last_tangent)
    }

    /// `FixedTangent(FirstFlag, LastFlag)` (`2d.cxx:133-143`).
    ///
    /// `&mut self` because OCCT's `FixedTangent` calls `FixTangentOnCurve`,
    /// which sets `myDone = true` on the checker itself (`2d.cxx:186`).
    pub fn fixed_tangent(
        &mut self,
        curve: &Geom2dBSplineCurve,
        first_flag: bool,
        last_flag: bool,
    ) -> Option<Geom2dBSplineCurve> {
        if !((self.fix_first_tangent && first_flag) || (self.fix_last_tangent && last_flag)) {
            return None;
        }
        let mut new_curve = curve.clone();
        self.fix_tangent_on_curve(&mut new_curve, first_flag, last_flag);
        Some(new_curve)
    }

    /// `FixTangent(FirstFlag, LastFlag)` (`2d.cxx:147-150`).
    pub fn fix_tangent(
        &mut self,
        curve: &mut Geom2dBSplineCurve,
        first_flag: bool,
        last_flag: bool,
    ) {
        self.fix_tangent_on_curve(curve, first_flag, last_flag);
    }

    /// `FixTangentOnCurve(theCurve, FirstFlag, LastFlag)` (`2d.cxx:153-187`).
    fn fix_tangent_on_curve(
        &mut self,
        curve: &mut Geom2dBSplineCurve,
        first_flag: bool,
        last_flag: bool,
    ) {
        if self.fix_first_tangent && first_flag {
            let xy1 = *pnt2d(curve, 1).xy();
            let xy2 = *pnt2d(curve, self.ind_second_pole as usize).xy();
            let nb_samples = f64::from(self.ind_second_pole - 1);
            let mut i = 2usize;
            while i < self.ind_second_pole as usize {
                let ii = f64::from(i as i32 - 1);
                let a_new_pole = GpPnt2d::from_xy(
                    xy1.multiplied_scalar(1.0 - ii / nb_samples)
                        .added(&xy2.multiplied_scalar(ii / nb_samples)),
                );
                set_pnt2d(curve, i, a_new_pole);
                i += 1;
            }
        }

        if self.fix_last_tangent && last_flag {
            let num_poles = curve.nb_poles();

            let xy1 = *pnt2d(curve, num_poles).xy();
            let xy2 = *pnt2d(curve, self.ind_prelast_pole as usize).xy();
            let nb_samples = f64::from(num_poles as i32 - self.ind_prelast_pole);
            let mut i = num_poles as i32 - 1;
            while i > self.ind_prelast_pole {
                let ii = f64::from(num_poles as i32 - i);
                let a_new_pole = GpPnt2d::from_xy(
                    xy1.multiplied_scalar(1.0 - ii / nb_samples)
                        .added(&xy2.multiplied_scalar(ii / nb_samples)),
                );
                set_pnt2d(curve, i as usize, a_new_pole);
                i -= 1;
            }
        }

        self.done = true;
    }
}

/// `Geom2d_BSplineCurve::Pole(Index)` 1-based.
fn pnt2d(curve: &Geom2dBSplineCurve, index: usize) -> GpPnt2d {
    GpPnt2d::new(curve.xs[index - 1], curve.ys[index - 1])
}

/// `Geom2d_BSplineCurve::SetPole(Index, P)` 1-based.
fn set_pnt2d(curve: &mut Geom2dBSplineCurve, index: usize, p: GpPnt2d) {
    curve.xs[index - 1] = p.x();
    curve.ys[index - 1] = p.y();
}
