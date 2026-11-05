//! `GccAna_Circ2d2TanRad(GccEnt_QualifiedLin&, GccEnt_QualifiedLin&,
//! Radius, Tolerance)` (`GccAna_Circ2d2TanRad_4.cxx:42-244`).
//!
//! Two offset lines per pair of sides, intersected with
//! `IntAna2d_AnaIntersection`; the qualifier of an unqualified argument comes
//! from which side of the line the centre lies on (`_4.cxx:200-222`).

use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::dir2d::GpDir2d;
use occt_core::gp::lin2d::GpLin2d;
use occt_core::gp::pnt2d::GpPnt2d;
use occt_core::intana2d::ana_intersection::IntAna2dAnaIntersection;
use occt_core::precision::ANGULAR;

use super::{axis_at, GccAnaCirc2d2TanRad};
use crate::gcc::gcc_ent::{GccEntPosition, GccEntQualifiedLin};

impl GccAnaCirc2d2TanRad {
    /// `GccAna_Circ2d2TanRad::GccAna_Circ2d2TanRad(Qualified1, Qualified2,
    /// Radius, Tolerance)` (`_4.cxx:42-244`).
    ///
    /// `_4.cxx:45` names the last parameter without a name: the line/line
    /// construction never consults the tolerance, so the port keeps the
    /// argument for signature parity and does not read it.
    pub fn new_lin_lin(
        qualified1: &GccEntQualifiedLin,
        qualified2: &GccEntQualifiedLin,
        radius: f64,
        _tolerance: f64,
    ) -> Self {
        let mut s = Self::empty();
        let mut cote1 = [0.0f64; 2];
        let mut cote2 = [0.0f64; 2];
        let mut nbrcote1 = 0usize;
        let mut nbrcote2 = 0usize;

        // `_4.cxx:68-72`.
        let q1_ok = qualified1.is_enclosed()
            || qualified1.is_outside()
            || qualified1.is_unqualified();
        let q2_ok = qualified2.is_enclosed()
            || qualified2.is_outside()
            || qualified2.is_unqualified();
        if !q1_ok || !q2_ok {
            panic!("GccEnt_BadQualifier in GccAna_Circ2d2TanRad");
        }

        let l1 = qualified1.qualified();
        let l2 = qualified2.qualified();
        let x1dir = l1.direction().x;
        let y1dir = l1.direction().y;
        let lx1loc = l1.location().x();
        let ly1loc = l1.location().y();
        let x2dir = l2.direction().x;
        let y2dir = l2.direction().y;
        let lx2loc = l2.location().x();
        let ly2loc = l2.location().y();
        let origin1 = GpPnt2d::new(lx1loc, ly1loc);
        let origin2 = GpPnt2d::new(lx2loc, ly2loc);
        // `_4.cxx:85-86` names these `normL1` / `normL2` but stores the line
        // directions, not their normals; the inside test at `:200-208` uses
        // them with `dc1`/`dc2` exactly as written.
        let norm_l1 = GpDir2d::new(x1dir, y1dir).expect("gp_Dir2d(L1.Direction())");
        let norm_l2 = GpDir2d::new(x2dir, y2dir).expect("gp_Dir2d(L2.Direction())");

        if radius < 0.0 {
            panic!("Standard_NegativeValue in GccAna_Circ2d2TanRad");
        }

        if l1.direction().is_parallel(l2.direction(), ANGULAR) {
            s.well_done = true;
        } else {
            // `_4.cxx:97-176`: the side table. The OCCT chain mixes `else if`
            // and separate `if`s (see the two `if`s at `:130` and `:138`, and
            // again at `:157` and `:165`); the assignments only ever add
            // entries, so the port keeps the same statement order.
            if qualified1.is_enclosed() && qualified2.is_enclosed() {
                nbrcote1 = 1;
                nbrcote2 = 1;
                cote1[0] = 1.0;
                cote2[0] = 1.0;
            } else if qualified1.is_enclosed() && qualified2.is_outside() {
                nbrcote1 = 1;
                nbrcote2 = 1;
                cote1[0] = 1.0;
                cote2[0] = -1.0;
            } else if qualified1.is_outside() && qualified2.is_enclosed() {
                nbrcote1 = 1;
                nbrcote2 = 1;
                cote1[0] = -1.0;
                cote2[0] = 1.0;
            } else if qualified1.is_outside() && qualified2.is_outside() {
                nbrcote1 = 1;
                nbrcote2 = 1;
                cote1[0] = -1.0;
                cote2[0] = -1.0;
            }
            if qualified1.is_enclosed() && qualified2.is_unqualified() {
                nbrcote1 = 1;
                nbrcote2 = 2;
                cote1[0] = 1.0;
                cote2[0] = 1.0;
                cote2[1] = -1.0;
            }
            if qualified1.is_unqualified() && qualified2.is_enclosed() {
                nbrcote1 = 2;
                nbrcote2 = 1;
                cote1[0] = 1.0;
                cote1[1] = -1.0;
                cote2[0] = 1.0;
            } else if qualified1.is_outside() && qualified2.is_unqualified() {
                nbrcote1 = 1;
                nbrcote2 = 2;
                cote1[0] = -1.0;
                cote2[0] = 1.0;
                cote2[1] = -1.0;
            }
            if qualified1.is_unqualified() && qualified2.is_outside() {
                nbrcote1 = 2;
                nbrcote2 = 1;
                cote1[0] = 1.0;
                cote1[1] = -1.0;
                cote2[0] = -1.0;
            } else if qualified1.is_unqualified() && qualified2.is_unqualified() {
                nbrcote1 = 2;
                nbrcote2 = 2;
                cote1[0] = 1.0;
                cote1[1] = -1.0;
                cote2[0] = 1.0;
                cote2[1] = -1.0;
            }

            // `_4.cxx:177-235`.
            for jcote1 in 1..=nbrcote1 {
                for jcote2 in 1..=nbrcote2 {
                    let linint1 = GpLin2d::from_pnt_dir(
                        GpPnt2d::new(
                            lx1loc - cote1[jcote1 - 1] * y1dir * radius,
                            ly1loc + cote1[jcote1 - 1] * x1dir * radius,
                        ),
                        *l1.direction(),
                    );
                    let linint2 = GpLin2d::from_pnt_dir(
                        GpPnt2d::new(
                            lx2loc - cote2[jcote2 - 1] * y2dir * radius,
                            ly2loc + cote2[jcote2 - 1] * x2dir * radius,
                        ),
                        *l2.direction(),
                    );
                    let mut intp = IntAna2dAnaIntersection::new();
                    intp.perform_lin_lin(&linint1, &linint2);
                    if intp.is_done() {
                        if !intp.is_empty() {
                            for i in 1..=intp.nb_points() {
                                s.nbr_sol += 1;
                                let n = s.nbr_sol;
                                let center = *intp.point(i).value();
                                s.cirsol[n - 1] = GpCirc2d::new(axis_at(center), radius);
                                let dc1 = GpDir2d::new(
                                    origin1.x() - center.x(),
                                    origin1.y() - center.y(),
                                )
                                .expect("gp_Dir2d(origin1 - center)");
                                let dc2 = GpDir2d::new(
                                    origin2.x() - center.x(),
                                    origin2.y() - center.y(),
                                )
                                .expect("gp_Dir2d(origin2 - center)");
                                if !qualified1.is_unqualified() {
                                    s.qualifier1[n - 1] = qualified1.qualifier();
                                } else if dc1.x * norm_l1.x + dc1.y * norm_l1.y > 0.0 {
                                    s.qualifier1[n - 1] = GccEntPosition::Outside;
                                } else {
                                    s.qualifier1[n - 1] = GccEntPosition::Enclosed;
                                }
                                if !qualified2.is_unqualified() {
                                    s.qualifier2[n - 1] = qualified2.qualifier();
                                } else if dc2.x * norm_l2.x + dc2.y * norm_l2.y > 0.0 {
                                    s.qualifier2[n - 1] = GccEntPosition::Outside;
                                } else {
                                    s.qualifier2[n - 1] = GccEntPosition::Enclosed;
                                }
                                s.the_same1[n - 1] = 0;
                                s.the_same2[n - 1] = 0;
                                s.pnttg1sol[n - 1] = GpPnt2d::new(
                                    center.x() + cote1[jcote1 - 1] * radius * y1dir,
                                    center.y() + cote1[jcote1 - 1] * radius * -x1dir,
                                );
                                s.pnttg2sol[n - 1] = GpPnt2d::new(
                                    center.x() + cote2[jcote2 - 1] * radius * y2dir,
                                    center.y() + cote2[jcote2 - 1] * radius * -x2dir,
                                );
                            }
                        }
                        s.well_done = true;
                    }
                }
            }
        }

        s.fill_solution_parameters();
        // `_4.cxx:237-243`.
        for i in 0..s.nbr_sol {
            s.pararg1[i] = occt_core::elib::clib2d::parameter_lin2d(&l1, &s.pnttg1sol[i]);
            s.pararg2[i] = occt_core::elib::clib2d::parameter_lin2d(&l2, &s.pnttg2sol[i]);
        }
        s
    }
}
