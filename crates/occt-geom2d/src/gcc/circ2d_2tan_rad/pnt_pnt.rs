//! `GccAna_Circ2d2TanRad(const gp_Pnt2d&, const gp_Pnt2d&, Radius, Tolerance)`
//! (`GccAna_Circ2d2TanRad_5.cxx:35-115`).
//!
//! `_5.cxx:60-63` reports `WellDone = true` with no solution when the two points
//! are further apart than `2 * Radius + Tol`. `_5.cxx:66-79` is the exact
//! two-diameter case: one circle on the midpoint. `_5.cxx:80-106` intersects the
//! two radius-`Radius` circles centred on the points.
//!
//! Translation notes:
//! - the two qualifiers are `GccEnt_noqualifier` in every branch and both
//!   `TheSame` flags stay `0`; the port writes them at every store as OCCT does,
//! - `_5.cxx:108-114` zeroes `pararg1`/`pararg2` **before** computing
//!   `par1sol`/`par2sol`, so both argument parameters are always `0`,
//! - there is no `TheSame`/`WellDone` initialisation loop in this constructor
//!   (`_5.cxx:52-54` only resets `NbrSol` and `WellDone`); [`super::GccAnaCirc2d2TanRad::empty`]
//!   covers the zeroed arrays.
//!
//! UNPORTED: none.

use occt_core::elib::clib2d::parameter_circ2d;
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::pnt2d::GpPnt2d;
use occt_core::intana2d::ana_intersection::IntAna2dAnaIntersection;

use super::{axis_at, GccAnaCirc2d2TanRad};
use crate::gcc::gcc_ent::GccEntPosition;

impl GccAnaCirc2d2TanRad {
    /// `GccAna_Circ2d2TanRad(const gp_Pnt2d&, const gp_Pnt2d&, double, double)`
    /// (`_5.cxx:35-115`).
    pub fn new_pnt_pnt(point1: &GpPnt2d, point2: &GpPnt2d, radius: f64, tolerance: f64) -> Self {
        let mut s = Self::empty();
        let tol = tolerance.abs();

        if radius < 0.0 {
            panic!("Standard_NegativeValue in GccAna_Circ2d2TanRad");
        }

        if point1.distance(point2) - radius * 2.0 > tol {
            // `_5.cxx:60-63`.
            s.well_done = true;
        } else if point1.distance(point2) - radius * 2.0 > 0.0 {
            // `_5.cxx:66-79`.
            s.well_done = true;
            s.nbr_sol = 1;
            let mid = GpPnt2d::new(
                0.5 * (point1.x() + point2.x()),
                0.5 * (point1.y() + point2.y()),
            );
            s.cirsol[0] = GpCirc2d::new(axis_at(mid), radius);
            s.qualifier1[0] = GccEntPosition::NoQualifier;
            s.qualifier2[0] = GccEntPosition::NoQualifier;
            s.the_same1[0] = 0;
            s.the_same2[0] = 0;
            s.pnttg1sol[0] = *point1;
            s.pnttg2sol[0] = *point2;
        } else {
            // `_5.cxx:80-106`.
            let c1 = GpCirc2d::new(axis_at(*point1), radius);
            let c2 = GpCirc2d::new(axis_at(*point2), radius);
            let mut intp = IntAna2dAnaIntersection::new();
            intp.perform_circ_circ(&c1, &c2);
            if intp.is_done() {
                if !intp.is_empty() {
                    for i in 1..=intp.nb_points() {
                        s.nbr_sol += 1;
                        let n = s.nbr_sol;
                        let center = *intp.point(i).value();
                        s.cirsol[n - 1] = GpCirc2d::new(axis_at(center), radius);
                        s.qualifier1[n - 1] = GccEntPosition::NoQualifier;
                        s.qualifier2[n - 1] = GccEntPosition::NoQualifier;
                        s.the_same1[n - 1] = 0;
                        s.the_same2[n - 1] = 0;
                        s.pnttg1sol[n - 1] = *point1;
                        s.pnttg2sol[n - 1] = *point2;
                    }
                }
                s.well_done = true;
            }
        }

        // `_5.cxx:108-114`.
        for i in 0..s.nbr_sol {
            s.pararg1[i] = 0.0;
            s.pararg2[i] = 0.0;
            s.par1sol[i] = parameter_circ2d(&s.cirsol[i], &s.pnttg1sol[i]);
            s.par2sol[i] = parameter_circ2d(&s.cirsol[i], &s.pnttg2sol[i]);
        }
        s
    }
}
