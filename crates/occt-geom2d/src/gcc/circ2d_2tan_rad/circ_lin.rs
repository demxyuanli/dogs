//! `GccAna_Circ2d2TanRad(GccEnt_QualifiedCirc&, GccEnt_QualifiedLin&,
//! Radius, Tolerance)` (`GccAna_Circ2d2TanRad_1.cxx:42-626`).
//!
//! The whole body turns on the signed side of the circle centre with respect
//! to the line, `d = -ydir * (cxloc - lxloc) + xdir * (cyloc - lyloc)`
//! (`_1.cxx:117` and the repeats at `:128`, `:153`, `:158`, ...). None of the
//! operands change inside the decision table, so the port evaluates it once in
//! `side` and reuses the value; the boundary cases and their `nbsol` codes
//! (`1..4`) are translated one for one.

use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::dir2d::GpDir2d;
use occt_core::gp::lin2d::GpLin2d;
use occt_core::gp::pnt2d::GpPnt2d;
use occt_core::intana2d::ana_intersection::IntAna2dAnaIntersection;

use super::{axis_at, GccAnaCirc2d2TanRad};
use crate::gcc::gcc_ent::{GccEntPosition, GccEntQualifiedCirc, GccEntQualifiedLin};

impl GccAnaCirc2d2TanRad {
    /// `GccAna_Circ2d2TanRad::GccAna_Circ2d2TanRad(Qualified1, Qualified2,
    /// Radius, Tolerance)` (`_1.cxx:42-626`).
    pub fn new_circ_lin(
        qualified1: &GccEntQualifiedCirc,
        qualified2: &GccEntQualifiedLin,
        radius: f64,
        tolerance: f64,
    ) -> Self {
        let mut s = Self::empty();
        let tol = tolerance.abs();
        // `_1.cxx:63-71`.
        let q1_ok = qualified1.is_enclosed()
            || qualified1.is_enclosing()
            || qualified1.is_outside()
            || qualified1.is_unqualified();
        let q2_ok = qualified2.is_enclosed()
            || qualified2.is_outside()
            || qualified2.is_unqualified();
        if !q1_ok || !q2_ok {
            panic!("GccEnt_BadQualifier in GccAna_Circ2d2TanRad");
        }
        // `for (i = 1; i <= 8; i++) TheSame2(i) = 0;` (`_1.cxx:72-75`).
        for i in 0..super::MAX_SOL {
            s.the_same2[i] = 0;
        }

        let mut ncote1 = 0usize;
        let mut ncote2 = 0usize;
        let mut nbsol = 0usize;
        let mut cote = 0.0f64;
        let mut ccote = 0.0f64;
        let mut cote1 = [0.0f64; 3];
        let mut cote2 = [0.0f64; 2];

        let c1 = qualified1.qualified();
        let l2 = qualified2.qualified();
        let xdir = l2.direction().x;
        let ydir = l2.direction().y;
        let norm_l2 = GpDir2d::new(-ydir, xdir).expect("gp_Dir2d(-ydir, xdir)");
        let lxloc = l2.location().x();
        let lyloc = l2.location().y();
        let cxloc = c1.location().x();
        let cyloc = c1.location().y();
        let r1 = c1.radius();
        let center1 = GpPnt2d::new(cxloc, cyloc);
        let origin2 = GpPnt2d::new(lxloc, lyloc);
        let dist = radius * 2.0 + r1;
        let distance = l2.distance(&center1);
        // `-ydir * (cxloc - lxloc) + xdir * (cyloc - lyloc)` (`_1.cxx:117`).
        let side = -ydir * (cxloc - lxloc) + xdir * (cyloc - lyloc);

        if radius < 0.0 {
            panic!("Standard_NegativeValue in GccAna_Circ2d2TanRad");
        }
        if distance - dist > tol {
            s.well_done = true;
        } else if qualified1.is_enclosed() {
            // `_1.cxx:108-220`.
            if qualified2.is_enclosed() {
                if distance - r1 > tol {
                    s.well_done = true;
                } else if ((side < 0.0) && (radius * 2.0 - r1 + distance > tol))
                    || ((side > 0.0) && (radius * 2.0 - r1 - distance > tol))
                    || (radius - r1 > tol)
                {
                    s.well_done = true;
                } else if ((side < 0.0) && (radius * 2.0 > r1 - distance))
                    || ((distance - r1).abs() < tol)
                    || ((side > 0.0) && (radius * 2.0 > r1 + distance + tol))
                {
                    cote = 1.0;
                    nbsol = 3;
                    ccote = 1.0;
                } else {
                    ncote1 = 1;
                    ncote2 = 1;
                    cote1[0] = r1 - radius;
                    cote2[0] = 1.0;
                    nbsol = 1;
                }
            } else if qualified2.is_outside() {
                if distance > r1 + tol {
                    s.well_done = true;
                } else if ((side > 0.0) && (radius * 2.0 - r1 + distance > tol))
                    || ((side < 0.0) && (radius * 2.0 - r1 - distance > tol))
                    || (radius - r1 > tol)
                {
                    s.well_done = true;
                } else if ((side > 0.0) && (radius * 2.0 > r1 - distance))
                    || ((r1 - radius).abs() < tol)
                    || ((side < 0.0) && (radius * 2.0 > r1 + distance))
                    || ((distance - r1).abs() < tol)
                {
                    cote = -1.0;
                    nbsol = 3;
                    ccote = 1.0;
                } else {
                    ncote1 = 1;
                    ncote2 = 1;
                    cote1[0] = r1 - radius;
                    cote2[0] = -1.0;
                    nbsol = 1;
                }
            } else if qualified2.is_unqualified() {
                if (distance - r1 > tol)
                    || (radius * 2.0 - r1 - distance > tol)
                    || (radius - r1 > tol)
                {
                    s.well_done = true;
                } else if distance > r1 {
                    if side > 0.0 {
                        cote = 1.0;
                    } else {
                        cote = -1.0;
                        // _1.cxx:202 -- nbsol = 3 is repeated at _1.cxx:204 right
                        // below, so the inner store is dead in OCCT too.
                    }
                    nbsol = 3;
                    ccote = 1.0;
                } else {
                    ncote1 = 1;
                    ncote2 = 2;
                    cote1[0] = r1 - radius;
                    cote2[0] = 1.0;
                    cote2[1] = -1.0;
                    nbsol = 1;
                }
            }
            // `if (nbsol == 3) { }` (`_1.cxx:217-219`): empty block.
        } else if qualified1.is_enclosing() {
            // `_1.cxx:221-302`.
            if qualified2.is_enclosed() {
                if (distance < r1 - tol)
                    || (radius < r1 - tol)
                    || (radius * 2.0 < distance + r1 - tol)
                    || (tol < r1 - radius)
                    || (side < 0.0)
                {
                    s.well_done = true;
                } else if (distance < r1) || (radius < r1) || (radius * 2.0 < distance + r1) {
                    cote = 1.0;
                    ccote = 1.0;
                    nbsol = 3;
                } else {
                    ncote1 = 1;
                    ncote2 = 1;
                    cote1[0] = radius - r1;
                    cote2[0] = 1.0;
                    nbsol = 1;
                }
            } else if qualified2.is_outside() {
                if (tol < r1 - distance)
                    || (tol < distance + r1 - radius * 2.0)
                    || (tol < r1 - radius)
                    || (side > 0.0)
                {
                    s.well_done = true;
                } else if (distance < r1) || (radius * 2.0 < distance + r1) || (radius < r1) {
                    cote = -1.0;
                    ccote = -1.0;
                    nbsol = 3;
                } else {
                    ncote1 = 1;
                    ncote2 = 1;
                    cote1[0] = radius - r1;
                    cote2[0] = -1.0;
                    nbsol = 1;
                }
            } else if qualified2.is_unqualified() {
                if (distance < r1 - tol)
                    || (radius * 2.0 < distance + r1 - tol)
                    || (radius < r1 - tol)
                {
                    s.well_done = true;
                } else if (distance < r1) || (radius * 2.0 < distance + r1) {
                    if side > 0.0 {
                        cote = 1.0;
                        ccote = 1.0;
                        nbsol = 3;
                    } else {
                        ccote = -1.0;
                        cote = -1.0;
                        nbsol = 3;
                    }
                } else {
                    ncote1 = 1;
                    ncote2 = 2;
                    cote1[0] = radius - r1;
                    cote2[0] = 1.0;
                    cote2[1] = -1.0;
                    nbsol = 1;
                }
            }
        } else if qualified1.is_outside() && qualified2.is_enclosed() {
            // `_1.cxx:303-332`.
            if ((side < 0.0) && (distance > r1 + tol)) || (distance > r1 + radius * 2.0 + tol) {
                s.well_done = true;
            } else if (side < 0.0) && (distance > r1) {
                cote = 1.0;
                nbsol = 2;
            } else if (side > 0.0) && (distance > r1 + 2.0 * radius) {
                cote = -1.0;
                nbsol = 2;
            } else {
                ncote1 = 1;
                ncote2 = 1;
                cote1[0] = radius + r1;
                cote2[0] = 1.0;
                nbsol = 1;
            }
        } else if qualified1.is_outside() && qualified2.is_outside() {
            // `_1.cxx:333-363`.
            if ((side > 0.0) && (distance > r1 + tol)) || (distance > r1 + radius * 2.0 + tol) {
                s.well_done = true;
            } else if (side > 0.0) && (distance > r1) {
                cote = -1.0;
                nbsol = 2;
            } else if (side < 0.0) && (distance > r1 + 2.0 * radius) {
                cote = 1.0;
                nbsol = 2;
            } else {
                ncote1 = 1;
                ncote2 = 1;
                cote1[0] = radius + r1;
                cote2[0] = -1.0;
                nbsol = 1;
            }
        } else if qualified1.is_outside() && qualified2.is_unqualified() {
            // `_1.cxx:364-392`.
            if distance > radius * 2.0 + r1 + tol {
                s.well_done = true;
            } else if distance > radius * 2.0 + r1 {
                if side > 0.0 {
                    cote = 1.0;
                } else {
                    cote = -1.0;
                }
                nbsol = 4;
            } else {
                ncote1 = 1;
                ncote2 = 2;
                cote1[0] = radius + r1;
                cote2[0] = 1.0;
                cote2[1] = -1.0;
                nbsol = 1;
            }
        } else if qualified1.is_unqualified() && qualified2.is_enclosed() {
            // `_1.cxx:393-425`.
            if (distance > r1 + radius * 2.0 + tol)
                || ((side > 0.0) && (distance > r1 + tol))
            {
                s.well_done = true;
            } else if (side > 0.0) && (distance > r1) {
                cote = 1.0;
                nbsol = 2;
            } else if (side < 0.0) && (distance > r1 + 2.0 * radius) {
                cote = -1.0;
                nbsol = 2;
            } else {
                ncote1 = 2;
                ncote2 = 1;
                cote1[0] = (radius - r1).abs();
                cote1[1] = radius + r1;
                cote2[0] = 1.0;
                nbsol = 1;
            }
        } else if qualified1.is_unqualified() && qualified2.is_outside() {
            // `_1.cxx:426-457`.
            if (distance > r1 + radius * 2.0 + tol)
                || ((side < 0.0) && (distance > r1 + tol))
            {
                s.well_done = true;
            } else if (side < 0.0) && (distance > r1) {
                cote = -1.0;
                nbsol = 2;
            } else if (side > 0.0) && (distance > r1 + 2.0 * radius) {
                cote = 1.0;
                nbsol = 2;
            } else {
                ncote1 = 2;
                ncote2 = 1;
                cote1[0] = (radius - r1).abs();
                cote1[1] = radius + r1;
                cote2[0] = -1.0;
                nbsol = 1;
            }
        } else if qualified1.is_unqualified() && qualified2.is_unqualified() {
            // `_1.cxx:458-487`.
            if distance > r1 + radius * 2.0 + tol {
                s.well_done = true;
            } else if distance > r1 + radius * 2.0 {
                if side > 0.0 {
                    cote = -1.0;
                } else {
                    cote = 1.0;
                }
                nbsol = 4;
            } else {
                ncote1 = 2;
                ncote2 = 2;
                cote1[0] = (radius - r1).abs();
                cote1[1] = radius + r1;
                cote2[0] = 1.0;
                cote2[1] = -1.0;
                nbsol = 1;
            }
        }

        // `_1.cxx:488-567`: the `nbsol == 1` intersection loop.
        if nbsol == 1 {
            for jcote1 in 1..=ncote1 {
                for jcote2 in 1..=ncote2 {
                    if cote1[jcote1 - 1] < tol {
                        continue;
                    }
                    let cirint = GpCirc2d::new(axis_at(center1), cote1[jcote1 - 1]);
                    let linint = GpLin2d::from_pnt_dir(
                        GpPnt2d::new(
                            lxloc - cote2[jcote2 - 1] * ydir * radius,
                            lyloc + cote2[jcote2 - 1] * xdir * radius,
                        ),
                        *l2.direction(),
                    );
                    let mut intp = IntAna2dAnaIntersection::new();
                    intp.perform_lin_circ(&linint, &cirint);
                    if intp.is_done() {
                        if !intp.is_empty() {
                            for i in 1..=intp.nb_points() {
                                s.nbr_sol += 1;
                                let n = s.nbr_sol;
                                let center = *intp.point(i).value();
                                s.cirsol[n - 1] = GpCirc2d::new(axis_at(center), radius);
                                let dc2 = GpDir2d::new(
                                    origin2.x() - center.x(),
                                    origin2.y() - center.y(),
                                )
                                .expect("gp_Dir2d(origin - center)");
                                let distcc1 = center.distance(&center1);
                                if !qualified1.is_unqualified() {
                                    s.qualifier1[n - 1] = qualified1.qualifier();
                                } else if (distcc1 + radius - r1).abs() < tol {
                                    s.qualifier1[n - 1] = GccEntPosition::Enclosed;
                                } else if (distcc1 - r1 - radius).abs() < tol {
                                    s.qualifier1[n - 1] = GccEntPosition::Outside;
                                } else {
                                    s.qualifier1[n - 1] = GccEntPosition::Enclosing;
                                }
                                if !qualified2.is_unqualified() {
                                    s.qualifier2[n - 1] = qualified2.qualifier();
                                } else if dc2.x * norm_l2.x + dc2.y * norm_l2.y > 0.0 {
                                    s.qualifier2[n - 1] = GccEntPosition::Outside;
                                } else {
                                    s.qualifier2[n - 1] = GccEntPosition::Enclosed;
                                }
                                s.the_same1[n - 1] = 0;
                                let dir_to_center = GpDir2d::new(
                                    center.x() - center1.x(),
                                    center.y() - center1.y(),
                                )
                                .expect("gp_Dir2d(center - center1)");
                                if (radius < r1) && center.distance(&center1) <= r1 {
                                    s.pnttg1sol[n - 1] = GpPnt2d::new(
                                        center.x() + radius * dir_to_center.x,
                                        center.y() + radius * dir_to_center.y,
                                    );
                                } else {
                                    s.pnttg1sol[n - 1] = GpPnt2d::new(
                                        center.x() - radius * dir_to_center.x,
                                        center.y() - radius * dir_to_center.y,
                                    );
                                }
                                s.pnttg2sol[n - 1] = GpPnt2d::new(
                                    center.x() + cote2[jcote2 - 1] * radius * ydir,
                                    center.y() + cote2[jcote2 - 1] * radius * -xdir,
                                );
                            }
                        }
                        s.well_done = true;
                    }
                }
            }
        } else if nbsol == 2 {
            // `_1.cxx:568-581`.
            let cen = GpPnt2d::new(
                center1.x() + cote * (r1 + radius) * -ydir,
                center1.y() + cote * (r1 + radius) * xdir,
            );
            s.cirsol[0] = GpCirc2d::new(axis_at(cen), radius);
            s.well_done = true;
            s.nbr_sol = 1;
            s.the_same1[0] = 0;
            s.qualifier1[0] = qualified1.qualifier();
            s.qualifier2[0] = qualified2.qualifier();
            s.pnttg1sol[0] =
                GpPnt2d::new(cen.x() + cote * radius * ydir, cen.y() + cote * radius * -xdir);
            s.pnttg2sol[0] = GpPnt2d::new(cen.x() + radius * ydir, cen.y() + radius * -xdir);
        } else if nbsol == 3 {
            // `_1.cxx:582-601`.
            s.well_done = true;
            s.nbr_sol = 1;
            let cen = GpPnt2d::new(
                center1.x() + cote * (r1 - radius) * ydir,
                center1.y() + cote * (r1 - radius) * -xdir,
            );
            s.cirsol[0] = GpCirc2d::new(axis_at(cen), radius);
            s.qualifier1[0] = qualified1.qualifier();
            s.qualifier2[0] = qualified2.qualifier();
            s.pnttg2sol[0] =
                GpPnt2d::new(cen.x() + cote * radius * ydir, cen.y() + cote * radius * -xdir);
            if (r1 - radius).abs() > 0.0 {
                s.pnttg1sol[0] =
                    GpPnt2d::new(cen.x() + ccote * radius * ydir, cen.y() + ccote * radius * -xdir);
            } else {
                s.the_same1[0] = 1;
            }
        } else if nbsol == 4 {
            // `_1.cxx:602-615`.
            let cent = GpPnt2d::new(
                center1.x() + cote * (r1 + radius) * -ydir,
                center1.y() + cote * (r1 + radius) * xdir,
            );
            s.cirsol[0] = GpCirc2d::new(axis_at(cent), radius);
            s.qualifier1[0] = qualified1.qualifier();
            s.qualifier2[0] = qualified2.qualifier();
            s.well_done = true;
            s.nbr_sol = 1;
            s.the_same1[0] = 0;
            s.pnttg1sol[0] = GpPnt2d::new(
                cent.x() + cote * radius * ydir,
                cent.y() + cote * radius * -xdir,
            );
            s.pnttg2sol[0] = GpPnt2d::new(
                cent.x() + cote * radius * -ydir,
                cent.y() + cote * radius * xdir,
            );
        }

        s.fill_solution_parameters();
        // `_1.cxx:616-625`: the argument parameter of the second tangency uses
        // the `gp_Lin2d` overload, the first the `gp_Circ2d` one.
        for i in 0..s.nbr_sol {
            if s.the_same1[i] == 0 {
                s.pararg1[i] = occt_core::elib::clib2d::parameter_circ2d(&c1, &s.pnttg1sol[i]);
            }
            s.pararg2[i] = occt_core::elib::clib2d::parameter_lin2d(&l2, &s.pnttg2sol[i]);
        }
        s
    }
}
