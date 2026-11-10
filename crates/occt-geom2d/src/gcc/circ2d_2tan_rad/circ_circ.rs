//! `GccAna_Circ2d2TanRad(GccEnt_QualifiedCirc&, GccEnt_QualifiedCirc&,
//! Radius, Tolerance)` (`GccAna_Circ2d2TanRad.cxx:43-802`).
//!
//! Every arm is the OCCT boundary-case table; the shared tail
//! (`cxx:692-801`) intersects the two concentric/offset auxiliary circles of
//! each surviving pair, then derives the qualifiers and tangency points from
//! the distance between the solution centre and each argument centre.

use occt_core::gp::circ2d::GpCirc2d;
use occt_core::intana2d::ana_intersection::IntAna2dAnaIntersection;
use occt_core::precision::CONFUSION;

use super::{axis_at, direction_between, offset_by, GccAnaCirc2d2TanRad, MAX_SOL};
use crate::gcc::gcc_ent::{GccEntPosition, GccEntQualifiedCirc};

impl GccAnaCirc2d2TanRad {
    /// `GccAna_Circ2d2TanRad::GccAna_Circ2d2TanRad(Qualified1, Qualified2,
    /// Radius, Tolerance)` (`cxx:43-802`).
    pub fn new_circ_circ(
        qualified1: &GccEntQualifiedCirc,
        qualified2: &GccEntQualifiedCirc,
        radius: f64,
        tolerance: f64,
    ) -> Self {
        let mut s = Self::empty();
        let tol = tolerance.abs();
        // `if (!(Q1.IsEnclosed() || Q1.IsEnclosing() || Q1.IsOutside() ||
        // Q1.IsUnqualified()) || !(...Q2...)) throw GccEnt_BadQualifier()`
        // (`cxx:63-71`).
        let q1_ok = qualified1.is_enclosed()
            || qualified1.is_enclosing()
            || qualified1.is_outside()
            || qualified1.is_unqualified();
        let q2_ok = qualified2.is_enclosed()
            || qualified2.is_enclosing()
            || qualified2.is_outside()
            || qualified2.is_unqualified();
        if !q1_ok || !q2_ok {
            panic!("GccEnt_BadQualifier in GccAna_Circ2d2TanRad");
        }

        let mut invers = false;
        // `gp_Circ2d C1 = Qualified1.Qualified(); gp_Circ2d C2 = ...`
        // (`cxx:72-73`). `C2` is reassigned in one branch of the enclosed /
        // enclosed arm (`cxx:128`), so both stay mutable.
        let c1 = qualified1.qualified();
        let mut c2 = qualified2.qualified();
        let mut r1 = c1.radius();
        let mut r2 = c2.radius();
        let rbid;
        let signe: i32;
        let dist;
        // `NCollection_Array1<gp_Circ2d> C(1, 8)` (`cxx:80`).
        let mut c = [c1; MAX_SOL];
        c[1] = c2;
        let mut nbsol = 0usize;
        let mut center1 = c1.location();
        let mut center2 = c2.location();
        let center3;

        if radius < 0.0 {
            panic!("Standard_NegativeValue in GccAna_Circ2d2TanRad");
        } else if c[1].location().is_equal(&c[0].location(), CONFUSION) {
            // `else if (C(2).Location().IsEqual(C(1).Location(),
            // Precision::Confusion())) WellDone = true;` (`cxx:92-95`).
            s.well_done = true;
        } else {
            let mut dir1 = direction_between(&c[0].location(), &c[1].location());
            dist = center2.distance(&center1);

            if qualified1.is_enclosed() && qualified2.is_enclosed() {
                // `cxx:99-155`.
                if radius * 2.0 - (r1 + r2 - dist).abs() > tol {
                    s.well_done = true;
                } else if (dist - r1 - r2 > tol)
                    || (tol < r1.max(r2) - dist - r1.min(r2))
                {
                    s.well_done = true;
                } else if (dist - r1 - r2).abs() <= tol
                    || (r1.max(r2) - dist - r1.min(r2)).abs() <= tol
                {
                    if (dist - r1 - r2).abs() <= tol {
                        rbid = r2 + (dist - r1 - r2) / 2.0;
                        signe = -1;
                    } else {
                        let r3;
                        if r1 > r2 {
                            r3 = r2;
                        } else {
                            c2 = c1;
                            r3 = r1;
                        }
                        rbid = -r3 + (r1.min(r2) + dist - r1.max(r2)) / 2.0;
                        signe = 1;
                    }
                    s.cirsol[0] = GpCirc2d::new(axis_at(offset_by(&center2, -rbid, &dir1)), radius);
                    s.qualifier1[0] = qualified1.qualifier();
                    s.qualifier2[0] = qualified2.qualifier();
                    s.pnttg1sol[0] = offset_by(&center2, signe as f64 * r2, &dir1);
                    s.pnttg2sol[0] = offset_by(&center1, r1, &dir1);
                    s.the_same1[0] = 0;
                    s.the_same2[0] = 0;
                    s.well_done = true;
                    s.nbr_sol = 1;
                } else {
                    c[0] = GpCirc2d::new(axis_at(c[0].location()), (radius - r1).abs());
                    c[1] = GpCirc2d::new(axis_at(c[1].location()), (radius - r2).abs());
                    nbsol = 1;
                }
            } else if (qualified1.is_enclosed() && qualified2.is_enclosing())
                || (qualified1.is_enclosing() && qualified2.is_enclosed())
            {
                // `cxx:156-236`.
                if qualified1.is_enclosing() {
                    let r3 = r1;
                    c[0] = c2;
                    r1 = r2;
                    c[1] = c1;
                    r2 = r3;
                    center3 = center1;
                    center1 = center2;
                    center2 = center3;
                    dir1.reverse();
                    invers = true;
                }
                if (r2 - radius > tol)
                    || (tol < radius - r1)
                    || (tol > r1 - dist - r2)
                    || (tol < radius * 2.0 - r1 - dist - r2)
                    || (tol < r1 + r2 - dist - radius * 2.0)
                {
                    s.well_done = true;
                } else if (r2 - radius).abs() <= tol || (r1 - radius).abs() <= tol {
                    if (r2 - radius).abs() <= tol {
                        c[1] = c2;
                        let r3 = r2;
                        s.the_same2[0] = 1;
                        s.the_same1[0] = 0;
                        s.pnttg1sol[0] = offset_by(&c[1].location(), r3, &dir1);
                    } else {
                        c[1] = c1;
                        let r3 = r1;
                        s.the_same1[0] = 1;
                        s.the_same2[0] = 0;
                        s.pnttg2sol[0] = offset_by(&c[1].location(), r3, &dir1);
                    }
                    s.well_done = true;
                    s.nbr_sol = 1;
                    s.cirsol[0] = GpCirc2d::new(axis_at(c[1].location()), radius);
                    s.qualifier1[0] = qualified1.qualifier();
                    s.qualifier2[0] = qualified2.qualifier();
                } else if (r2 + dist - r1).abs() <= tol
                    || (radius * 2.0 - r1 - dist - r2).abs() < tol
                {
                    if (r2 + dist - r1).abs() <= tol {
                        signe = 1;
                    } else {
                        signe = -1;
                    }
                    s.well_done = true;
                    s.nbr_sol = 1;
                    s.cirsol[0] = GpCirc2d::new(
                        axis_at(offset_by(&center1, signe as f64 * (r1 - radius), &dir1)),
                        radius,
                    );
                    s.qualifier1[0] = qualified1.qualifier();
                    s.qualifier2[0] = qualified2.qualifier();
                    s.pnttg1sol[0] = offset_by(&center1, signe as f64 * r1, &dir1);
                    s.pnttg2sol[0] = offset_by(&center2, r2, &dir1);
                    s.the_same1[0] = 0;
                    s.the_same2[0] = 0;
                } else {
                    c[0] = GpCirc2d::new(axis_at(c[0].location()), (radius - r1).abs());
                    c[1] = GpCirc2d::new(axis_at(c[1].location()), (radius - r2).abs());
                    nbsol = 1;
                }
            } else if (qualified1.is_enclosed() && qualified2.is_outside())
                || (qualified1.is_outside() && qualified2.is_enclosed())
            {
                // `cxx:237-306`.
                if qualified1.is_outside() {
                    c[1] = c1;
                    c[0] = c2;
                    let r3 = r1;
                    r1 = r2;
                    r2 = r3;
                    center3 = center1;
                    center1 = center2;
                    center2 = center3;
                    dir1.reverse();
                    invers = true;
                }
                if (radius - r1 > tol)
                    || (dist - r2 - r1 > tol)
                    || (r2 - dist - r1 + radius * 2.0 > tol)
                    || (r1 - r2 - dist - radius * 2.0 > tol)
                {
                    s.well_done = true;
                } else if ((radius - r1 > 0.0) && (dist - r1 - r2).abs() <= tol)
                    || (r1 - r2 + dist - radius * 2.0).abs() <= tol
                    || (r1 - r2 - dist - radius * 2.0).abs() <= tol
                    || (dist - r1 - r2).abs() <= tol
                {
                    if (r1 - r2 + dist - radius * 2.0).abs() <= tol {
                        signe = -1;
                    } else {
                        signe = 1;
                        if (radius - r1 > 0.0) && (dist - r1 - r2).abs() <= tol {
                            r2 = r1;
                        } else if (dist - r1 - r2).abs() <= tol {
                            r2 = r1;
                            if (r1 - radius).abs() <= tol {
                                s.the_same1[0] = 1;
                            }
                        }
                    }
                    s.well_done = true;
                    s.nbr_sol = 1;
                    s.cirsol[0] = GpCirc2d::new(
                        axis_at(offset_by(&center1, signe as f64 * (r1 - radius), &dir1)),
                        radius,
                    );
                    s.qualifier1[0] = qualified1.qualifier();
                    s.qualifier2[0] = qualified2.qualifier();
                    s.pnttg1sol[0] = offset_by(&center1, signe as f64 * r1, &dir1);
                    s.pnttg2sol[0] = offset_by(&center2, signe as f64 * r2, &dir1);
                    s.the_same1[0] = 0;
                    s.the_same2[0] = 0;
                } else {
                    c[0] = GpCirc2d::new(axis_at(c[0].location()), (radius - r1).abs());
                    c[1] = GpCirc2d::new(axis_at(c[1].location()), radius + r2);
                    nbsol = 1;
                }
            } else if (qualified1.is_enclosed() && qualified2.is_unqualified())
                || (qualified1.is_unqualified() && qualified2.is_enclosed())
            {
                // `cxx:308-362`.
                if qualified1.is_unqualified() {
                    c[0] = c2;
                    c[1] = c1;
                    let r3 = r1;
                    r1 = r2;
                    r2 = r3;
                    center3 = center1;
                    center1 = center2;
                    center2 = center3;
                    dir1.reverse();
                    invers = true;
                }
                if (radius - r1 > tol) || (dist - r2 - r1 > tol) {
                    s.well_done = true;
                } else if (dist - r2 - r1).abs() <= tol || (radius - r1).abs() <= tol {
                    s.well_done = true;
                    s.nbr_sol = 1;
                    s.cirsol[0] = GpCirc2d::new(
                        axis_at(offset_by(&center1, r1 - radius, &dir1)),
                        radius,
                    );
                    s.qualifier1[0] = qualified1.qualifier();
                    s.qualifier2[0] = qualified2.qualifier();
                    s.pnttg2sol[0] = offset_by(&center1, dist - r2, &dir1);
                    s.the_same2[0] = 0;
                    if (radius - r1).abs() > 0.0 {
                        s.the_same1[0] = 1;
                    } else {
                        s.the_same1[0] = 0;
                        s.pnttg1sol[0] = offset_by(&center1, r1, &dir1);
                    }
                } else {
                    c[2] = GpCirc2d::new(axis_at(c[0].location()), (radius - r1).abs());
                    c[3] = GpCirc2d::new(axis_at(c[1].location()), radius + r2);
                    c[0] = GpCirc2d::new(axis_at(c[0].location()), (radius - r1).abs());
                    c[1] = GpCirc2d::new(axis_at(c[1].location()), (radius - r2).abs());
                    nbsol = 2;
                }
            } else if qualified1.is_enclosing() && qualified2.is_enclosing() {
                // `cxx:362-429`.
                if (tol < r1.max(r2) - radius)
                    || (tol < r1.max(r2) - dist - r1.min(r2))
                    || (dist + r1 + r2 - radius * 2.0 > tol)
                {
                    s.well_done = true;
                } else if (dist + r1.min(r2) - r1.max(r2)).abs() <= tol
                    || (r1 + r2 + dist - 2.0 * radius).abs() <= tol
                {
                    if (dist + r1.min(r2) - r1.max(r2)).abs() <= tol {
                        signe = 1;
                    } else {
                        signe = -1;
                    }
                    s.well_done = true;
                    s.nbr_sol = 1;
                    s.cirsol[0] = GpCirc2d::new(
                        axis_at(offset_by(&center1, signe as f64 * (r1 - radius), &dir1)),
                        radius,
                    );
                    s.qualifier1[0] = qualified1.qualifier();
                    s.qualifier2[0] = qualified2.qualifier();
                    s.pnttg1sol[0] = offset_by(&center1, r1, &dir1);
                    s.pnttg2sol[0] = offset_by(&center2, signe as f64 * r2, &dir1);
                    s.the_same1[0] = 0;
                    s.the_same2[0] = 0;
                } else if (radius - r1.max(r2)).abs() <= tol {
                    s.well_done = true;
                    s.nbr_sol = 1;
                    if r1 > r2 {
                        c[0] = c1;
                        c[1] = c2;
                        s.the_same1[0] = 1;
                        s.the_same2[0] = 0;
                        s.pnttg2sol[0] = offset_by(&center1, dist + r2, &dir1);
                    } else {
                        c[0] = c2;
                        c[1] = c1;
                        s.the_same1[0] = 0;
                        s.the_same2[0] = 1;
                        s.pnttg1sol[0] = offset_by(&center1, -r1, &dir1);
                    }
                    s.cirsol[0] = GpCirc2d::new(axis_at(c[0].location()), radius);
                    s.qualifier1[0] = qualified1.qualifier();
                    s.qualifier2[0] = qualified2.qualifier();
                } else {
                    c[0] = GpCirc2d::new(axis_at(c[0].location()), (radius - r1).abs());
                    c[1] = GpCirc2d::new(axis_at(c[1].location()), (radius - r2).abs());
                    nbsol = 1;
                }
            } else if (qualified1.is_enclosing() && qualified2.is_outside())
                || (qualified1.is_outside() && qualified2.is_enclosing())
            {
                // `cxx:430-480`.
                if qualified1.is_outside() {
                    c[0] = c2;
                    c[1] = c1;
                    let r3 = r1;
                    r1 = r2;
                    r2 = r3;
                    center3 = center1;
                    center1 = center2;
                    center2 = center3;
                    dir1.reverse();
                    invers = true;
                }
                if (r1 - radius > tol)
                    || (tol < r1 + r2 - dist)
                    || (dist - r2 + r1 - radius * 2.0 > tol)
                {
                    s.well_done = true;
                } else if (r1 - radius).abs() <= tol || (r1 + r2 - dist).abs() <= tol
                    || (dist - r2 + r1 - radius * 2.0).abs() <= tol
                {
                    s.well_done = true;
                    s.nbr_sol = 1;
                    if (r1 - radius).abs() <= tol || (r1 + r2 - dist).abs() <= tol {
                        s.the_same1[0] = 1;
                    } else {
                        s.the_same1[0] = 0;
                    }
                    s.the_same2[0] = 0;
                    s.cirsol[0] = GpCirc2d::new(
                        axis_at(offset_by(&center1, r1 - radius, &dir1)),
                        radius,
                    );
                    s.qualifier1[0] = qualified1.qualifier();
                    s.qualifier2[0] = qualified2.qualifier();
                    s.pnttg1sol[0] = offset_by(&center1, -r1, &dir1);
                    s.pnttg2sol[0] = offset_by(&center1, dist - r2, &dir1);
                } else {
                    c[0] = GpCirc2d::new(axis_at(c[0].location()), (radius - r1).abs());
                    c[1] = GpCirc2d::new(axis_at(c[1].location()), radius + r2);
                    nbsol = 1;
                }
            } else if (qualified1.is_enclosing() && qualified2.is_unqualified())
                || (qualified1.is_unqualified() && qualified2.is_enclosing())
            {
                // `cxx:482-537`.
                if qualified1.is_unqualified() {
                    c[0] = c2;
                    c[1] = c1;
                    let r3 = r1;
                    r1 = r2;
                    r2 = r3;
                    center3 = center1;
                    center1 = center2;
                    center2 = center3;
                    invers = true;
                    dir1.reverse();
                }
                if (tol < r1 - dist - r2) || (r1 - radius > tol) {
                    s.well_done = true;
                } else if (r1 - radius).abs() <= tol || (r1 - dist - r2).abs() > 0.0 {
                    if (r1 - radius).abs() <= tol {
                        s.the_same1[0] = 1;
                        if (radius - r2).abs() <= tol && center1.distance(&center2) <= tol {
                            s.the_same2[0] = 1;
                        }
                    } else if (r1 - dist - r2).abs() > 0.0 {
                        s.the_same1[0] = 0;
                        s.the_same2[0] = 0;
                    }
                    s.well_done = true;
                    s.nbr_sol = 1;
                    s.cirsol[0] = GpCirc2d::new(
                        axis_at(offset_by(&center1, radius - r1, &dir1)),
                        radius,
                    );
                    s.qualifier1[0] = qualified1.qualifier();
                    s.qualifier2[0] = qualified2.qualifier();
                    s.pnttg1sol[0] = offset_by(&center1, r1, &dir1);
                    s.pnttg2sol[0] = offset_by(&center1, dist + r2, &dir1);
                } else {
                    c[2] = GpCirc2d::new(axis_at(c[0].location()), (radius - r1).abs());
                    c[3] = GpCirc2d::new(axis_at(c[1].location()), radius + r2);
                    c[0] = GpCirc2d::new(axis_at(c[0].location()), (radius - r1).abs());
                    c[1] = GpCirc2d::new(axis_at(c[1].location()), (radius - r2).abs());
                    nbsol = 2;
                }
            } else if qualified1.is_outside() && qualified2.is_outside() {
                // `cxx:538-594`.
                if tol < r1.max(r2) - dist - r1.min(r2) {
                    s.well_done = true;
                } else if dist - r1 - r2 - radius * 2.0 > tol {
                    s.well_done = true;
                } else if (dist + r1.min(r2) - r1.max(r2)).abs() <= tol {
                    s.well_done = true;
                    s.nbr_sol = 1;
                    if r1 < r2 {
                        signe = -1;
                    } else {
                        signe = 1;
                    }
                    s.cirsol[0] = GpCirc2d::new(
                        axis_at(offset_by(&center1, signe as f64 * (radius + r1), &dir1)),
                        radius,
                    );
                    s.qualifier1[0] = qualified1.qualifier();
                    s.qualifier2[0] = qualified2.qualifier();
                    s.pnttg1sol[0] = offset_by(&center1, signe as f64 * r1, &dir1);
                    s.pnttg2sol[0] = s.pnttg1sol[0];
                    s.the_same1[0] = 0;
                    s.the_same2[0] = 0;
                } else if (dist - r1 - r2 - radius * 2.0).abs() <= tol {
                    s.well_done = true;
                    s.nbr_sol = 1;
                    s.cirsol[0] = GpCirc2d::new(
                        axis_at(offset_by(&center1, r1 + radius, &dir1)),
                        radius,
                    );
                    s.qualifier1[0] = qualified1.qualifier();
                    s.qualifier2[0] = qualified2.qualifier();
                    s.pnttg1sol[0] = offset_by(&center1, r1, &dir1);
                    s.pnttg2sol[0] = offset_by(&center2, -r2, &dir1);
                    s.the_same1[0] = 0;
                    s.the_same2[0] = 0;
                } else {
                    c[0] = GpCirc2d::new(axis_at(c[0].location()), radius + r1);
                    c[1] = GpCirc2d::new(axis_at(c[1].location()), radius + r2);
                    nbsol = 1;
                }
            } else if (qualified1.is_outside() && qualified2.is_unqualified())
                || (qualified1.is_unqualified() && qualified2.is_outside())
            {
                // `cxx:596-656`.
                if qualified1.is_unqualified() {
                    c[0] = c2;
                    let r3 = r1;
                    r1 = r2;
                    c[1] = c1;
                    r2 = r3;
                    center3 = center1;
                    center1 = center2;
                    center2 = center3;
                    dir1.reverse();
                    invers = true;
                }
                if tol < r1 - dist - r2 {
                    s.well_done = true;
                } else if (tol < r2 - dist - r1) && (tol < radius * 2.0 - r2 - dist + r1) {
                    s.well_done = true;
                } else if (dist - r1 - r2 > tol) && (tol < dist - r1 - r2 - radius * 2.0) {
                    s.well_done = true;
                } else if (r1 - r2 - dist).abs() <= tol
                    || ((dist - r1 - r2).abs() <= tol
                        && (radius * 2.0 - dist + r1 + r2).abs() <= tol)
                    || ((dist + r1 - r2).abs() <= tol
                        && (r2 + dist - r1 - radius * 2.0).abs() <= tol)
                {
                    s.well_done = true;
                    s.nbr_sol = 1;
                    s.cirsol[0] = GpCirc2d::new(
                        axis_at(offset_by(&center1, radius + r1, &dir1)),
                        radius,
                    );
                    s.qualifier1[0] = qualified1.qualifier();
                    s.qualifier2[0] = qualified2.qualifier();
                    s.pnttg1sol[0] = offset_by(&center1, r1, &dir1);
                    s.pnttg2sol[0] = offset_by(&center1, dist + r2, &dir1);
                    s.the_same1[0] = 0;
                    s.the_same2[0] = 0;
                } else {
                    c[2] = GpCirc2d::new(axis_at(c[0].location()), radius + r1);
                    c[3] = GpCirc2d::new(axis_at(c[1].location()), radius + r2);
                    c[0] = GpCirc2d::new(axis_at(c[0].location()), radius + r1);
                    c[1] = GpCirc2d::new(axis_at(c[0].location()), (radius - r2).abs());
                    nbsol = 2;
                }
            } else if qualified1.is_unqualified() && qualified2.is_unqualified() {
                // `cxx:658-691`.
                if (dist - r1 - r2 > tol) && (tol < dist - r1 - r2 - radius * 2.0) {
                    s.well_done = true;
                } else if (r1.max(r2) - dist - r1.min(r2) > tol)
                    && (r1.max(r2) - dist - r1.min(r2) - radius * 2.0 > tol)
                {
                    s.well_done = true;
                } else {
                    // cxx:666-672 -- OCCT computes `p3` and assigns `dist` in that
                    // branch, but neither is read afterwards: only R1/R2 feed the
                    // solutions below, so both are omitted here.
                    let p4 = dist - r1 - r2;
                    let p5 = radius * 2.0 - dist + r1 + r2;
                    if p4 > 0.0 && p5 < 0.0 {
                        r1 = dist - r2 - radius * 2.0;
                    }
                    c[0] = GpCirc2d::new(axis_at(center1), (radius - r1).abs());
                    c[1] = GpCirc2d::new(axis_at(center2), (radius - r2).abs());
                    c[2] = GpCirc2d::new(axis_at(center1), (radius - r1).abs());
                    c[3] = GpCirc2d::new(axis_at(center2), radius + r2);
                    c[4] = GpCirc2d::new(axis_at(center1), radius + r1);
                    c[5] = GpCirc2d::new(axis_at(center2), (radius - r2).abs());
                    c[6] = GpCirc2d::new(axis_at(center1), radius + r1);
                    c[7] = GpCirc2d::new(axis_at(center2), radius + r2);
                    nbsol = 4;
                }
            }

            // `cxx:692-774`: the shared pair loop.
            if nbsol > 0 {
                for j in 1..=nbsol {
                    let mut intp = IntAna2dAnaIntersection::new();
                    intp.perform_circ_circ(&c[2 * j - 2], &c[2 * j - 1]);
                    if intp.is_done() {
                        if !intp.is_empty() {
                            for i in 1..=intp.nb_points() {
                                s.nbr_sol += 1;
                                let n = s.nbr_sol;
                                let center = *intp.point(i).value();
                                s.cirsol[n - 1] = GpCirc2d::new(axis_at(center), radius);
                                dir1 = direction_between(&center1, &center);
                                let dir2 = direction_between(&center2, &center);
                                let distcc1 = center.distance(&center1);
                                let distcc2 = center.distance(&center2);
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
                                } else if (distcc2 + radius - r2).abs() < tol {
                                    s.qualifier2[n - 1] = GccEntPosition::Enclosed;
                                } else if (distcc2 - r2 - radius).abs() < tol {
                                    s.qualifier2[n - 1] = GccEntPosition::Outside;
                                } else {
                                    s.qualifier2[n - 1] = GccEntPosition::Enclosing;
                                }
                                if (center.distance(&center1) > r1)
                                    && (radius < center.distance(&center1) + r1)
                                {
                                    s.pnttg1sol[n - 1] = offset_by(&center, -radius, &dir1);
                                } else if (center.distance(&center1) < r1) && (radius < r1) {
                                    s.pnttg1sol[n - 1] = offset_by(&center, radius, &dir1);
                                } else {
                                    s.pnttg1sol[n - 1] = offset_by(&center, -radius, &dir1);
                                }
                                if (center.distance(&center2) > r2)
                                    && (radius < center.distance(&center2) + r2)
                                {
                                    s.pnttg2sol[n - 1] = offset_by(&center, -radius, &dir2);
                                } else if (center.distance(&center2) < r2) && (radius < r2) {
                                    s.pnttg2sol[n - 1] = offset_by(&center, radius, &dir2);
                                } else {
                                    s.pnttg2sol[n - 1] = offset_by(&center, -radius, &dir2);
                                }
                                s.the_same1[n - 1] = 0;
                                s.the_same2[n - 1] = 0;
                            }
                        }
                        s.well_done = true;
                    }
                }
            }
        }

        // `cxx:777-786`: swap the tangency points when the arguments were
        // swapped.
        if invers {
            for i in 0..s.nbr_sol {
                let psav = s.pnttg1sol[i];
                s.pnttg1sol[i] = s.pnttg2sol[i];
                s.pnttg2sol[i] = psav;
            }
        }

        // `cxx:788-801`: the per-solution parameters, using the arguments as
        // they were handed to the constructor (`C1`, `C2`).
        for i in 0..s.nbr_sol {
            s.par1sol[i] =
                occt_core::elib::clib2d::parameter_circ2d(&s.cirsol[i], &s.pnttg1sol[i]);
            if s.the_same1[i] == 0 {
                s.pararg1[i] = occt_core::elib::clib2d::parameter_circ2d(&c1, &s.pnttg1sol[i]);
            }
            s.par2sol[i] =
                occt_core::elib::clib2d::parameter_circ2d(&s.cirsol[i], &s.pnttg2sol[i]);
            if s.the_same2[i] == 0 {
                s.pararg2[i] = occt_core::elib::clib2d::parameter_circ2d(&c2, &s.pnttg2sol[i]);
            }
        }

        s
    }
}
