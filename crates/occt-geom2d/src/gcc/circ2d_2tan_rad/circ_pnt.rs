//! `GccAna_Circ2d2TanRad(GccEnt_QualifiedCirc&, const gp_Pnt2d&, Radius,
//! Tolerance)` (`GccAna_Circ2d2TanRad_2.cxx:43-291`).
//!
//! `_2.cxx:92-206` is the boundary table: three exact configurations collapse to
//! `nbsol = -1` (one circle placed on the line through `C1`'s centre and
//! `Point2`, `deport` away) and the rest build the auxiliary pair
//! `C(2j-1) = gp_Circ2d(C1.XAxis(), |Radius -/+ R1|)` /
//! `C(2j) = gp_Circ2d(gp_Ax2d(Point2, dirx), Radius)` for
//! `_2.cxx:207-262` to intersect.
//!
//! Translation notes:
//! - `dir1` starts as `gp_Dir2d(Point2 - C1.Location())` (`_2.cxx:83`), so a
//!   `Point2` coincident with the centre raises `Standard_ConstructionError` in
//!   OCCT; the port's `expect` panics with the same message,
//! - `dir1` is **reassigned** inside the intersection loop
//!   (`_2.cxx:240`, `Center - center1`), which changes the direction the
//!   `nbsol < 0` arm uses if that arm ran after a loop that never does; the port
//!   keeps the `mut` binding,
//! - `_2.cxx:246-256` picks the tangency point of the first argument through a
//!   three-arm test whose first and third arms are the same expression
//!   (`Center - Radius * dir1`), the middle one (`Center + Radius * dir1`) being
//!   taken only when `dist < R1 && Radius < R1`; kept as written,
//! - `_2.cxx:286-292` computes `pararg1` only when `TheSame1(i) == 0` and always
//!   zeroes `pararg2`; the port keeps both,
//! - `gp_Circ2d(C1.XAxis(), r)` (`_2.cxx:113`, `:136`, `:155`, `:194`, `:200`,
//!   `:203`) keeps `C1`'s X direction, unlike the `gp_Ax2d(center, dirx)` used
//!   for every solution circle; see [`super::circ_on_x_axis`].
//!
//! UNPORTED: none.

use occt_core::elib::clib2d::parameter_circ2d;
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::dir2d::GpDir2d;
use occt_core::gp::pnt2d::GpPnt2d;
use occt_core::intana2d::ana_intersection::IntAna2dAnaIntersection;

use super::{axis_at, circ_on_x_axis, GccAnaCirc2d2TanRad, MAX_SOL};
use crate::gcc::gcc_ent::{GccEntPosition, GccEntQualifiedCirc};

impl GccAnaCirc2d2TanRad {
    /// `GccAna_Circ2d2TanRad(const GccEnt_QualifiedCirc&, const gp_Pnt2d&,
    /// double, double)` (`_2.cxx:43-291`).
    pub fn new_circ_pnt(
        qualified1: &GccEntQualifiedCirc,
        point2: &GpPnt2d,
        radius: f64,
        tolerance: f64,
    ) -> Self {
        let mut s = Self::empty();
        let tol = tolerance.abs();
        // `_2.cxx:64-68`.
        let q1_ok = qualified1.is_enclosed()
            || qualified1.is_enclosing()
            || qualified1.is_outside()
            || qualified1.is_unqualified();
        if !q1_ok {
            panic!("GccEnt_BadQualifier in GccAna_Circ2d2TanRad");
        }
        // `for (i = 1; i <= 4; i++) { TheSame1(i) = 0; TheSame2(i) = 0; }`
        // (`_2.cxx:70-74`).
        for i in 0..MAX_SOL {
            s.the_same1[i] = 0;
            s.the_same2[i] = 0;
        }

        let mut deport = 0.0f64;
        let mut signe = 0i32;
        let mut nbsol = 0i32;

        let c1 = qualified1.qualified();
        // `NCollection_Array1<gp_Circ2d> C(1, 4)` (`_2.cxx:79`).
        let mut c: [GpCirc2d; MAX_SOL] = [c1; MAX_SOL];
        let r1 = c1.radius();
        let distance = c1.location().distance(point2);
        let dispc1 = c1.distance(point2);
        let center1 = c1.location();
        let mut dir1 = GpDir2d::new(point2.x() - center1.x(), point2.y() - center1.y())
            .expect("gp_Dir2d(Point2 - C1.Location())");

        if radius < 0.0 {
            panic!("Standard_NegativeValue in GccAna_Circ2d2TanRad");
        }

        if dispc1 - radius * 2.0 > tol {
            // `_2.cxx:92-95`.
            s.well_done = true;
        } else if qualified1.is_enclosed() {
            // `_2.cxx:96-116`.
            if distance - r1 > tol || radius - r1 > tol {
                s.well_done = true;
            } else if (distance - r1).abs() < tol {
                nbsol = -1;
                deport = r1 - radius;
                signe = 1;
            } else {
                c[0] = circ_on_x_axis(&c1, (radius - r1).abs());
                c[1] = GpCirc2d::new(axis_at(*point2), radius);
                nbsol = 1;
            }
        } else if qualified1.is_enclosing() {
            // `_2.cxx:117-140`.
            if tol < r1 - distance || tol < r1 - radius {
                s.well_done = true;
            } else if (distance - r1).abs() < tol {
                nbsol = -1;
                deport = r1 - radius;
                signe = 1;
            } else {
                c[0] = circ_on_x_axis(&c1, (radius - r1).abs());
                c[1] = GpCirc2d::new(axis_at(*point2), radius);
                nbsol = 1;
            }
        } else if qualified1.is_outside() {
            // `_2.cxx:141-160`.
            if tol < r1 - distance {
                s.well_done = true;
            } else if (distance - r1).abs() < tol || (dispc1 - radius * 2.0).abs() < tol {
                nbsol = -1;
                deport = r1 + radius;
                signe = -1;
            } else {
                c[0] = circ_on_x_axis(&c1, radius + r1);
                c[1] = GpCirc2d::new(axis_at(*point2), radius);
                nbsol = 1;
            }
        } else if qualified1.is_unqualified() {
            // `_2.cxx:161-206`.
            if (dispc1 - radius * 2.0).abs() < tol {
                s.well_done = true;
                let center = GpPnt2d::new(
                    center1.x() + (distance - radius) * dir1.x,
                    center1.y() + (distance - radius) * dir1.y,
                );
                s.cirsol[0] = GpCirc2d::new(axis_at(center), radius);
                if (center.distance(&center1) - r1).abs() < tol {
                    s.qualifier1[0] = GccEntPosition::Enclosed;
                } else {
                    s.qualifier1[0] = GccEntPosition::Outside;
                }
                s.qualifier2[0] = GccEntPosition::NoQualifier;
                s.pnttg1sol[0] = GpPnt2d::new(
                    center.x() - radius * dir1.x,
                    center.y() - radius * dir1.y,
                );
                s.pnttg2sol[0] = *point2;
                s.well_done = true;
                s.nbr_sol = 1;
            } else if (r1 - radius).abs() < tol && (distance - r1).abs() < tol {
                s.cirsol[0] = c1;
                s.qualifier1[0] = GccEntPosition::Unqualified;
                s.qualifier2[0] = GccEntPosition::NoQualifier;
                s.the_same1[0] = 1;
                s.pnttg2sol[0] = *point2;
                s.well_done = true;
                s.nbr_sol = 1;
                c[0] = circ_on_x_axis(&c1, radius + r1);
                c[1] = GpCirc2d::new(axis_at(*point2), radius);
                nbsol = 1;
            } else {
                c[0] = circ_on_x_axis(&c1, (radius - r1).abs());
                c[1] = GpCirc2d::new(axis_at(*point2), radius);
                c[2] = circ_on_x_axis(&c1, radius + r1);
                c[3] = GpCirc2d::new(axis_at(*point2), radius);
                nbsol = 2;
            }
        }

        if nbsol > 0 {
            // `_2.cxx:207-262`.
            for j in 1..=nbsol as usize {
                let mut intp = IntAna2dAnaIntersection::new();
                intp.perform_circ_circ(&c[2 * j - 2], &c[2 * j - 1]);
                if intp.is_done() {
                    if !intp.is_empty() {
                        for i in 1..=intp.nb_points() {
                            s.nbr_sol += 1;
                            let n = s.nbr_sol;
                            let center = *intp.point(i).value();
                            s.cirsol[n - 1] = GpCirc2d::new(axis_at(center), radius);
                            let distcc1 = center1.distance(&center);
                            if !qualified1.is_unqualified() {
                                s.qualifier1[n - 1] = qualified1.qualifier();
                            } else if (distcc1 + radius - r1).abs() < tol {
                                s.qualifier1[n - 1] = GccEntPosition::Enclosed;
                            } else if (distcc1 - r1 - radius).abs() < tol {
                                s.qualifier1[n - 1] = GccEntPosition::Outside;
                            } else {
                                s.qualifier1[n - 1] = GccEntPosition::Enclosing;
                            }
                            s.qualifier2[n - 1] = GccEntPosition::NoQualifier;
                            // `_2.cxx:240`: `dir1` is overwritten from the
                            // current centre, not from `Point2`.
                            dir1 = GpDir2d::new(center.x() - center1.x(), center.y() - center1.y())
                                .expect("gp_Dir2d(Center - center1)");
                            let dist_center = center.distance(&center1);
                            if dist_center > c1.radius() && radius < dist_center + c1.radius() {
                                s.pnttg1sol[n - 1] = GpPnt2d::new(
                                    center.x() - radius * dir1.x,
                                    center.y() - radius * dir1.y,
                                );
                            } else if dist_center < c1.radius() && radius < c1.radius() {
                                s.pnttg1sol[n - 1] = GpPnt2d::new(
                                    center.x() + radius * dir1.x,
                                    center.y() + radius * dir1.y,
                                );
                            } else {
                                s.pnttg1sol[n - 1] = GpPnt2d::new(
                                    center.x() - radius * dir1.x,
                                    center.y() - radius * dir1.y,
                                );
                            }
                            s.pnttg2sol[n - 1] = *point2;
                        }
                    }
                    s.well_done = true;
                }
            }
        } else if nbsol < 0 {
            // `_2.cxx:263-282`.
            let center = GpPnt2d::new(
                center1.x() + deport * dir1.x,
                center1.y() + deport * dir1.y,
            );
            s.cirsol[0] = GpCirc2d::new(axis_at(center), radius);
            s.qualifier1[0] = qualified1.qualifier();
            s.qualifier2[0] = GccEntPosition::NoQualifier;
            if deport.abs() <= tol && (radius - r1).abs() <= tol {
                s.the_same1[0] = 1;
            } else {
                s.pnttg1sol[0] = GpPnt2d::new(
                    center.x() + signe as f64 * radius * dir1.x,
                    center.y() + signe as f64 * radius * dir1.y,
                );
            }
            s.pnttg2sol[0] = *point2;
            s.well_done = true;
            s.nbr_sol = 1;
        }

        // `_2.cxx:284-293`.
        s.fill_solution_parameters();
        for i in 0..s.nbr_sol {
            if s.the_same1[i] == 0 {
                s.pararg1[i] = parameter_circ2d(&c1, &s.pnttg1sol[i]);
            }
            s.pararg2[i] = 0.0;
        }
        s
    }
}
