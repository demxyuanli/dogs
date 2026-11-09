//! `GccAna_Circ2d2TanRad(GccEnt_QualifiedLin&, const gp_Pnt2d&, Radius,
//! Tolerance)` (`GccAna_Circ2d2TanRad_3.cxx:46-271`).
//!
//! `_3.cxx:92-168` is the boundary table, keyed on the signed side
//! `d = -ydir * (cxloc - lxloc) + xdir * (cyloc - lyloc)` of the point with
//! respect to the line. It produces either `nbsol == 1` (one parallel of `L1`,
//! `nbcote` sides) or `nbsol == 2` (a single centre straight off `Point2`).
//! `_3.cxx:171-248` handles the two `nbsol == 1` sub-cases: `displ1 < 1.e-10`
//! builds the centres directly, otherwise the parallel is intersected with the
//! circle of centre `Point2`. `_3.cxx:250-262` is the `nbsol == 2` arm.
//!
//! Translation notes:
//! - `displ1 < 1.e-10` (`_3.cxx:173`) is a **hard-coded** literal, not
//!   `Tol`; kept as `1.0e-10`,
//! - `_3.cxx:110` presets `NbrSol = 1` in the `IsEnclosed` / `displ1 > 2*Radius`
//!   arm, but that arm is the `nbsol == 2` one, which sets `NbrSol = 1` again
//!   (`_3.cxx:253`); the preset is a dead store, kept for parity,
//! - `_3.cxx:250-261` (the `nbsol == 2` arm) assigns `qualifier2(1)` but never
//!   `qualifier1(1)`; OCCT leaves that slot at its default-constructed value and
//!   the port leaves [`GccEntPosition::Unqualified`],
//! - the `nbsol == 1` intersection loop stops at `NbrSol < 2` (`_3.cxx:217`),
//!   a per-solution clamp the port keeps,
//! - `_3.cxx:239` offsets the tangency point with `gp_XY(ydir, -xdir)`, the
//!   negative of the `gp_XY(-ydir, xdir)` used for the parallel itself
//!   (`_3.cxx:208-209`); both are preserved,
//! - `_3.cxx:263-270` computes `pararg1` against `L1` unconditionally, without
//!   the `TheSame1` guard the circle/circle and circle/point constructors use.
//!
//! UNPORTED: none.

use occt_core::elib::clib2d::parameter_lin2d;
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::dir2d::GpDir2d;
use occt_core::gp::lin2d::GpLin2d;
use occt_core::gp::pnt2d::GpPnt2d;
use occt_core::intana2d::ana_intersection::IntAna2dAnaIntersection;

use super::{axis_at, GccAnaCirc2d2TanRad};
use crate::gcc::gcc_ent::{GccEntPosition, GccEntQualifiedLin};

impl GccAnaCirc2d2TanRad {
    /// `GccAna_Circ2d2TanRad(const GccEnt_QualifiedLin&, const gp_Pnt2d&,
    /// double, double)` (`_3.cxx:46-271`).
    pub fn new_lin_pnt(
        qualified1: &GccEntQualifiedLin,
        point2: &GpPnt2d,
        radius: f64,
        tolerance: f64,
    ) -> Self {
        let mut s = Self::empty();
        let tol = tolerance.abs();
        // `_3.cxx:67-71`.
        let q1_ok = qualified1.is_enclosed()
            || qualified1.is_outside()
            || qualified1.is_unqualified();
        if !q1_ok {
            panic!("GccEnt_BadQualifier in GccAna_Circ2d2TanRad");
        }

        let mut nbsol = 0i32;
        let mut nbcote = 0usize;
        // `NCollection_Array1<double> cote(1, 2)` (`_3.cxx:74`).
        let mut cote = [0.0f64; 2];

        let l1 = qualified1.qualified();
        let displ1 = l1.distance(point2);
        let xdir = l1.direction().x;
        let ydir = l1.direction().y;
        let lxloc = l1.location().x();
        let lyloc = l1.location().y();
        let origin1 = GpPnt2d::new(lxloc, lyloc);
        let norm_l1 = GpDir2d::new(-ydir, xdir).expect("gp_Dir2d(-ydir, xdir)");
        let cxloc = point2.x();
        let cyloc = point2.y();
        // `_3.cxx:100`, `:123`, `:148`, `:154`.
        let side = -ydir * (cxloc - lxloc) + xdir * (cyloc - lyloc);

        if radius < 0.0 {
            panic!("Standard_NegativeValue in GccAna_Circ2d2TanRad");
        }

        if displ1 - radius * 2.0 > tol {
            // `_3.cxx:95-98`.
            s.well_done = true;
        } else if qualified1.is_enclosed() {
            // `_3.cxx:99-121`.
            if side < 0.0 {
                s.well_done = true;
            } else if displ1 - radius * 2.0 > 0.0 {
                nbsol = 2;
                s.nbr_sol = 1;
                nbcote = 1;
                cote[0] = 1.0;
            } else {
                nbsol = 1;
                nbcote = 1;
                cote[0] = 1.0;
            }
        } else if qualified1.is_outside() {
            // `_3.cxx:122-143`.
            if side > 0.0 {
                s.well_done = true;
            } else if displ1 - radius * 2.0 > 0.0 {
                nbsol = 2;
                nbcote = 1;
                cote[0] = -1.0;
            } else {
                nbsol = 1;
                nbcote = 1;
                cote[0] = -1.0;
            }
        } else if qualified1.is_unqualified() {
            // `_3.cxx:144-168`.
            if displ1 - radius * 2.0 > 0.0 {
                if side > 0.0 {
                    nbsol = 2;
                    nbcote = 1;
                    cote[0] = 1.0;
                } else if side < 0.0 {
                    nbsol = 2;
                    nbcote = 1;
                    cote[0] = -1.0;
                }
            } else {
                nbsol = 1;
                nbcote = 2;
                cote[0] = 1.0;
                cote[1] = -1.0;
            }
        }

        if nbsol == 1 {
            // `_3.cxx:171-248`.
            if displ1 < 1.0e-10 {
                // `_3.cxx:173-202`: `Point2` lies on the line, build the
                // centres directly.
                for jcote in 1..=nbcote {
                    s.nbr_sol += 1;
                    let n = s.nbr_sol;
                    let center = GpPnt2d::new(
                        cxloc - cote[jcote - 1] * ydir * radius,
                        cyloc + cote[jcote - 1] * xdir * radius,
                    );
                    s.cirsol[n - 1] = GpCirc2d::new(axis_at(center), radius);
                    s.qualifier2[n - 1] = GccEntPosition::NoQualifier;
                    if !qualified1.is_unqualified() {
                        s.qualifier1[n - 1] = qualified1.qualifier();
                    } else if cote[jcote - 1] > 0.0 {
                        s.qualifier1[n - 1] = GccEntPosition::Outside;
                    } else {
                        s.qualifier1[n - 1] = GccEntPosition::Enclosed;
                    }
                    s.the_same1[n - 1] = 0;
                    s.the_same2[n - 1] = 0;
                    s.pnttg1sol[n - 1] = *point2;
                    s.pnttg2sol[n - 1] = *point2;
                }
                s.well_done = true;
            } else {
                // `_3.cxx:204-247`.
                let cirint = GpCirc2d::new(axis_at(*point2), radius);
                for jcote in 1..=nbcote {
                    let linint = GpLin2d::from_pnt_dir(
                        GpPnt2d::new(
                            lxloc - cote[jcote - 1] * ydir * radius,
                            lyloc + cote[jcote - 1] * xdir * radius,
                        ),
                        *l1.direction(),
                    );
                    let mut intp = IntAna2dAnaIntersection::new();
                    intp.perform_lin_circ(&linint, &cirint);
                    if intp.is_done() {
                        if !intp.is_empty() {
                            for i in 1..=intp.nb_points() {
                                // `i <= Intp.NbPoints() && NbrSol < 2`.
                                if s.nbr_sol >= 2 {
                                    break;
                                }
                                s.nbr_sol += 1;
                                let n = s.nbr_sol;
                                let center = *intp.point(i).value();
                                s.cirsol[n - 1] = GpCirc2d::new(axis_at(center), radius);
                                let dc1 = GpDir2d::new(
                                    origin1.x() - center.x(),
                                    origin1.y() - center.y(),
                                )
                                .expect("gp_Dir2d(origin1 - Center)");
                                s.qualifier2[n - 1] = GccEntPosition::NoQualifier;
                                if !qualified1.is_unqualified() {
                                    s.qualifier1[n - 1] = qualified1.qualifier();
                                } else if dc1.x * norm_l1.x + dc1.y * norm_l1.y > 0.0 {
                                    s.qualifier1[n - 1] = GccEntPosition::Outside;
                                } else {
                                    s.qualifier1[n - 1] = GccEntPosition::Enclosed;
                                }
                                s.the_same1[n - 1] = 0;
                                s.the_same2[n - 1] = 0;
                                s.pnttg1sol[n - 1] = GpPnt2d::new(
                                    center.x() + cote[jcote - 1] * radius * ydir,
                                    center.y() + cote[jcote - 1] * radius * -xdir,
                                );
                                s.pnttg2sol[n - 1] = *point2;
                            }
                        }
                        s.well_done = true;
                    }
                }
            }
        } else if nbsol == 2 {
            // `_3.cxx:250-262`.
            let center = GpPnt2d::new(
                point2.x() + cote[0] * radius * -ydir,
                point2.y() + cote[0] * radius * xdir,
            );
            s.well_done = true;
            s.nbr_sol = 1;
            s.cirsol[0] = GpCirc2d::new(axis_at(center), radius);
            // `_3.cxx:256-258`: `qualifier1(1)` is never assigned here.
            s.qualifier2[0] = GccEntPosition::NoQualifier;
            s.the_same1[0] = 0;
            s.the_same2[0] = 0;
            s.pnttg1sol[0] = GpPnt2d::new(
                center.x() + cote[0] * radius * ydir,
                center.y() + cote[0] * radius * -xdir,
            );
            s.pnttg2sol[0] = *point2;
        }

        // `_3.cxx:263-270`.
        s.fill_solution_parameters();
        for i in 0..s.nbr_sol {
            s.pararg1[i] = parameter_lin2d(&l1, &s.pnttg1sol[i]);
            s.pararg2[i] = 0.0;
        }
        s
    }
}
