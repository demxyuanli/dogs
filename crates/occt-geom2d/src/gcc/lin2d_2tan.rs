//! `GccAna_Lin2d2Tan` (`GccAna_Lin2d2Tan.cxx:37-621`): lines through two
//! points, through a point and tangent to a circle, or tangent to two circles.
//!
//! | OCCT | lines |
//! |---|---|
//! | `GccAna_Lin2d2Tan(gp_Pnt2d, gp_Pnt2d, Tolerance)` | `cxx:37-71` |
//! | `GccAna_Lin2d2Tan(GccEnt_QualifiedCirc, gp_Pnt2d, Tolerance)` | `cxx:83-174` |
//! | `GccAna_Lin2d2Tan(GccEnt_QualifiedCirc, GccEnt_QualifiedCirc, Tolerance)` | `cxx:194-540` |
//! | `IsDone`..`Tangency2` | `cxx:542-621` |
//!
//! Translation notes:
//! - the qualifier decisions of the two-circle constructor are transcribed
//!   branch for branch, including the places where OCCT stores
//!   `Qualified1.Qualifier()` into `qualifier2` (`cxx:272`, `:303`, `:317`,
//!   `:343`, `:357`, `:395`, `:409`, `:424`, `:460`, `:474`, `:489`, `:511`,
//!   `:525`) and the `pnttg1sol(2) = P1` of `cxx:425`,
//! - `gp_Pnt2d::Rotated` is [`GpPnt2d::rotated`] with the same rotation
//!   transformation, so `P2.Location().Rotated(C1.Location(), angle)` is
//!   `c2.location().rotated(&c1.location(), angle)`,
//! - `NbSolutions` (`cxx:547-550`) and `ThisSolution` (`cxx:552-560`) do *not*
//!   test `WellDone` in OCCT; only the index is checked, and the port keeps
//!   that,
//! - `GccEnt_BadQualifier`, `StdFail_NotDone` and `Standard_OutOfRange` panic
//!   with the OCCT exception name.
//!
//! UNPORTED: none - all three constructors and all six accessors.

use occt_core::elib::clib2d::{parameter_circ2d, parameter_lin2d};
use occt_core::gp::dir2d::GpDir2d;
use occt_core::gp::lin2d::GpLin2d;
use occt_core::gp::pnt2d::GpPnt2d;
use occt_core::gp::xy::GpXY;

use super::gcc_ent::{GccEntPosition, GccEntQualifiedCirc};

/// `NCollection_Array1<...>(1, 4)`: the largest solution count, reached by the
/// two-circle constructor (`cxx:197-205`).
pub const MAX_SOL: usize = 4;

/// `GccAna_Lin2d2Tan` (`GccAna_Lin2d2Tan.hxx`).
#[derive(Debug, Clone, Copy)]
pub struct GccAnaLin2d2Tan {
    pub(crate) linsol: [GpLin2d; MAX_SOL],
    pub(crate) qualifier1: [GccEntPosition; MAX_SOL],
    pub(crate) qualifier2: [GccEntPosition; MAX_SOL],
    pub(crate) pnttg1sol: [GpPnt2d; MAX_SOL],
    pub(crate) pnttg2sol: [GpPnt2d; MAX_SOL],
    pub(crate) par1sol: [f64; MAX_SOL],
    pub(crate) par2sol: [f64; MAX_SOL],
    pub(crate) pararg1: [f64; MAX_SOL],
    pub(crate) pararg2: [f64; MAX_SOL],
    pub(crate) well_done: bool,
    pub(crate) nbr_sol: usize,
}

impl GccAnaLin2d2Tan {
    /// The `: linsol(1, 1), ...` / `(1, 2)` / `(1, 4)` initializer list every
    /// constructor shares (`cxx:40-48`, `cxx:86-94`, `cxx:197-205`).
    pub(crate) fn empty() -> Self {
        Self {
            linsol: [GpLin2d::default(); MAX_SOL],
            qualifier1: [GccEntPosition::Unqualified; MAX_SOL],
            qualifier2: [GccEntPosition::Unqualified; MAX_SOL],
            pnttg1sol: [GpPnt2d::zero(); MAX_SOL],
            pnttg2sol: [GpPnt2d::zero(); MAX_SOL],
            par1sol: [0.0; MAX_SOL],
            par2sol: [0.0; MAX_SOL],
            pararg1: [0.0; MAX_SOL],
            pararg2: [0.0; MAX_SOL],
            well_done: false,
            nbr_sol: 0,
        }
    }

    /// `GccAna_Lin2d2Tan(ThePoint1, ThePoint2, Tolerance)` (`cxx:37-71`).
    pub fn from_two_points(the_point1: &GpPnt2d, the_point2: &GpPnt2d, tolerance: f64) -> Self {
        let mut s = Self::empty();
        let tol = tolerance.abs();
        let dist = the_point1.distance(the_point2);
        s.qualifier1[0] = GccEntPosition::NoQualifier;
        s.qualifier2[0] = GccEntPosition::NoQualifier;
        if dist >= tol {
            let dir = GpDir2d::new(
                the_point2.x() - the_point1.x(),
                the_point2.y() - the_point1.y(),
            )
            .expect("gp_VectorWithNullMagnitude in GccAna_Lin2d2Tan");
            s.linsol[0] = GpLin2d::from_pnt_dir(*the_point1, dir);
            s.well_done = true;
            s.nbr_sol = 1;
            s.pnttg1sol[0] = *the_point1;
            s.pnttg2sol[0] = *the_point2;
            s.par1sol[0] = parameter_lin2d(&s.linsol[0], &s.pnttg1sol[0]);
            s.par2sol[0] = parameter_lin2d(&s.linsol[0], &s.pnttg2sol[0]);
            s.pararg1[0] = 0.0;
            s.pararg2[0] = 0.0;
        }
        s
    }

    /// `GccAna_Lin2d2Tan(Qualified1, ThePoint, Tolerance)` (`cxx:83-174`).
    pub fn from_circle_and_point(
        qualified1: &GccEntQualifiedCirc,
        the_point: &GpPnt2d,
        tolerance: f64,
    ) -> Self {
        let mut s = Self::empty();
        let tol = tolerance.abs();
        // `cxx:100-105`.
        if !(qualified1.is_enclosed()
            || qualified1.is_enclosing()
            || qualified1.is_outside()
            || qualified1.is_unqualified())
        {
            panic!("GccEnt_BadQualifier in GccAna_Lin2d2Tan");
        }
        let c1 = qualified1.qualified();
        let r1 = c1.radius();
        let center = c1.location();

        if qualified1.is_enclosed() {
            panic!("GccEnt_BadQualifier in GccAna_Lin2d2Tan");
        } else if tol < r1 - the_point.distance(&center) {
            // `cxx:114-117`: the point is inside the circle.
            s.well_done = true;
        } else if (the_point.distance(&center) - r1).abs() <= tol {
            // `cxx:118-129`: the point is on the circle.
            let dir = GpDir2d::new(
                the_point.x() - center.x(),
                the_point.y() - center.y(),
            )
            .expect("gp_VectorWithNullMagnitude in GccAna_Lin2d2Tan");
            s.linsol[0] = GpLin2d::from_pnt_dir(
                *the_point,
                GpDir2d::new(-dir.y, dir.x).expect("gp_Dir2d"),
            );
            s.qualifier1[0] = qualified1.qualifier();
            s.qualifier2[0] = GccEntPosition::NoQualifier;
            s.well_done = true;
            s.nbr_sol = 1;
            s.pnttg1sol[0] = *the_point;
            s.pnttg2sol[0] = *the_point;
        } else {
            // `cxx:130-166`.
            let mut signe = 1.0f64;
            let dist = the_point.distance(&center);
            let d = dist - (dist * dist - r1 * r1).sqrt();
            if qualified1.is_enclosing() {
                signe = 1.0;
                s.nbr_sol = 1;
            } else if qualified1.is_outside() {
                signe = -1.0;
                s.nbr_sol = 1;
            } else if qualified1.is_unqualified() {
                signe = 1.0;
                s.nbr_sol = 2;
            }
            for i in 0..s.nbr_sol {
                let mut p1 = center.rotated(the_point, (signe * r1 / dist).asin());
                let d1 = GpDir2d::new(the_point.x() - p1.x(), the_point.y() - p1.y())
                    .expect("gp_VectorWithNullMagnitude in GccAna_Lin2d2Tan");
                p1 = GpPnt2d::from_xy(
                    p1.xy().added(&GpXY::new(d1.x, d1.y).multiplied_scalar(d)),
                );
                let dir = GpDir2d::new(the_point.x() - p1.x(), the_point.y() - p1.y())
                    .expect("gp_VectorWithNullMagnitude in GccAna_Lin2d2Tan");
                s.linsol[i] = GpLin2d::from_pnt_dir(p1, dir);
                s.qualifier1[i] = qualified1.qualifier();
                s.qualifier2[i] = GccEntPosition::NoQualifier;
                s.pnttg1sol[i] = p1;
                s.pnttg2sol[i] = *the_point;
                signe = -signe;
            }
            s.well_done = true;
        }
        // `cxx:167-173`.
        for i in 0..s.nbr_sol {
            s.par1sol[i] = parameter_lin2d(&s.linsol[i], &s.pnttg1sol[i]);
            s.par2sol[i] = parameter_lin2d(&s.linsol[i], &s.pnttg2sol[i]);
            s.pararg1[i] = parameter_circ2d(&c1, &s.pnttg1sol[i]);
            s.pararg2[i] = 0.0;
        }
        s
    }

    /// `GccAna_Lin2d2Tan(Qualified1, Qualified2, Tolerance)` (`cxx:194-540`).
    pub fn from_two_circles(
        qualified1: &GccEntQualifiedCirc,
        qualified2: &GccEntQualifiedCirc,
        tolerance: f64,
    ) -> Self {
        let mut s = Self::empty();
        let tol = tolerance.abs();
        // `cxx:211-218`.
        let q1_ok = qualified1.is_enclosed()
            || qualified1.is_enclosing()
            || qualified1.is_outside()
            || qualified1.is_unqualified();
        let q2_ok = qualified2.is_enclosed()
            || qualified2.is_enclosing()
            || qualified2.is_outside()
            || qualified2.is_unqualified();
        if !q1_ok || !q2_ok {
            panic!("GccEnt_BadQualifier in GccAna_Lin2d2Tan");
        }
        let c1 = qualified1.qualified();
        let c2 = qualified2.qualified();
        if qualified1.is_enclosed() || qualified2.is_enclosed() {
            panic!("GccEnt_BadQualifier in GccAna_Lin2d2Tan");
        }

        let r1 = c1.radius();
        let r2 = c2.radius();
        let loc1 = c1.location();
        let loc2 = c2.location();
        let mut d1: GpDir2d;
        let mut signe = 1.0f64;
        let dist = loc1.distance(&loc2);
        let (rmin, rmax) = if r1 < r2 { (r1, r2) } else { (r2, r1) };

        if tol < rmax - dist - rmin {
            // `cxx:234-237`.
            s.well_done = true;
        } else if !qualified1.is_unqualified() || !qualified2.is_unqualified() {
            // `cxx:238-498`.
            if qualified1.is_enclosing() && qualified2.is_enclosing() {
                // `cxx:241-278`.
                if (dist + rmin - rmax).abs() <= tol && dist >= tol {
                    d1 = if r1 < r2 {
                        dir_between(&loc2, &loc1)
                    } else {
                        dir_between(&loc1, &loc2)
                    };
                    let p1 = GpPnt2d::from_xy(
                        loc1.xy().added(&GpXY::new(d1.x, d1.y).multiplied_scalar(r1)),
                    );
                    s.linsol[0] =
                        GpLin2d::from_pnt_dir(p1, GpDir2d::new(-d1.y, d1.x).expect("gp_Dir2d"));
                    s.qualifier1[0] = qualified1.qualifier();
                    s.qualifier2[0] = qualified2.qualifier();
                    s.pnttg1sol[0] = p1;
                    s.pnttg2sol[0] = p1;
                    s.well_done = true;
                    s.nbr_sol = 1;
                } else {
                    let p1 = loc2.rotated(&loc1, ((r1 - r2) / dist).asin());
                    d1 = dir_between(&loc1, &p1);
                    let p1 = GpPnt2d::from_xy(
                        loc1.xy().added(&GpXY::new(r1 * d1.y, -r1 * d1.x)),
                    );
                    s.linsol[0] = GpLin2d::from_pnt_dir(p1, d1);
                    s.qualifier1[0] = qualified1.qualifier();
                    s.qualifier2[0] = qualified1.qualifier();
                    s.pnttg1sol[0] = p1;
                    s.pnttg2sol[0] = GpPnt2d::from_xy(
                        loc2.xy().added(&GpXY::new(r2 * d1.y, -r2 * d1.x)),
                    );
                    s.well_done = true;
                    s.nbr_sol = 1;
                }
            } else if (qualified1.is_enclosing() && qualified2.is_outside())
                || (qualified2.is_enclosing() && qualified1.is_outside())
            {
                // `cxx:279-323`.
                signe = if qualified1.is_enclosing() && qualified2.is_outside() {
                    1.0
                } else {
                    -1.0
                };
                if r1 + r2 - dist > tol {
                    s.well_done = true;
                } else if (dist - r1 - r2).abs() < tol && dist > tol {
                    d1 = dir_between(&loc1, &loc2);
                    let p1 = GpPnt2d::from_xy(
                        loc1.xy().added(&GpXY::new(d1.x, d1.y).multiplied_scalar(r1)),
                    );
                    s.linsol[0] =
                        GpLin2d::from_pnt_dir(p1, GpDir2d::new(-d1.y, d1.x).expect("gp_Dir2d"));
                    s.qualifier1[0] = qualified1.qualifier();
                    s.qualifier2[0] = qualified1.qualifier();
                    s.pnttg1sol[0] = p1;
                    s.pnttg2sol[0] = p1;
                    s.well_done = true;
                    s.nbr_sol = 1;
                } else {
                    let p1 = loc2.rotated(&loc1, (signe * (r1 + r2) / dist).asin());
                    d1 = dir_between(&loc1, &p1);
                    let p1 = GpPnt2d::from_xy(
                        loc1.xy()
                            .added(&GpXY::new(r1 * d1.y, -r1 * d1.x).multiplied_scalar(signe)),
                    );
                    s.linsol[0] = GpLin2d::from_pnt_dir(p1, d1);
                    s.qualifier1[0] = qualified1.qualifier();
                    s.qualifier2[0] = qualified1.qualifier();
                    s.pnttg1sol[0] = p1;
                    s.pnttg2sol[0] = GpPnt2d::from_xy(loc2.xy().added(
                        &GpXY::new(-r2 * d1.y, r2 * d1.x).multiplied_scalar(signe),
                    ));
                    s.well_done = true;
                    s.nbr_sol = 1;
                }
            } else if qualified1.is_outside() && qualified2.is_outside() {
                // `cxx:324-363`.
                if (dist + rmin - rmax).abs() < tol && dist > tol {
                    d1 = if r1 < r2 {
                        dir_between(&loc2, &loc1)
                    } else {
                        dir_between(&loc1, &loc2)
                    };
                    let p1 = GpPnt2d::from_xy(
                        loc1.xy().added(&GpXY::new(d1.x, d1.y).multiplied_scalar(r1)),
                    );
                    s.linsol[0] =
                        GpLin2d::from_pnt_dir(p1, GpDir2d::new(d1.y, -d1.x).expect("gp_Dir2d"));
                    s.qualifier1[0] = qualified1.qualifier();
                    s.qualifier2[0] = qualified1.qualifier();
                    s.pnttg1sol[0] = p1;
                    s.pnttg2sol[0] = p1;
                    s.well_done = true;
                    s.nbr_sol = 1;
                } else {
                    let p1 = loc2.rotated(&loc1, ((r2 - r1) / dist).asin());
                    d1 = dir_between(&loc1, &p1);
                    let p1 = GpPnt2d::from_xy(
                        loc1.xy().added(&GpXY::new(-r1 * d1.y, r1 * d1.x)),
                    );
                    s.linsol[0] = GpLin2d::from_pnt_dir(p1, d1);
                    s.qualifier1[0] = qualified1.qualifier();
                    s.qualifier2[0] = qualified1.qualifier();
                    s.pnttg1sol[0] = p1;
                    s.pnttg2sol[0] = GpPnt2d::from_xy(
                        loc2.xy().added(&GpXY::new(-r2 * d1.y, r2 * d1.x)),
                    );
                    s.well_done = true;
                    s.nbr_sol = 1;
                }
            } else {
                // `cxx:364-497`.
                if (qualified1.is_unqualified() && qualified2.is_enclosing())
                    || (qualified2.is_unqualified() && qualified1.is_enclosing())
                {
                    // `cxx:366-430`.
                    signe = if qualified2.is_unqualified() { 1.0 } else { -1.0 };
                    if (dist + rmin - rmax).abs() < tol && dist > tol {
                        d1 = if r1 < r2 {
                            dir_between(&loc2, &loc1)
                        } else {
                            dir_between(&loc1, &loc2)
                        };
                        let p1 = GpPnt2d::from_xy(
                            loc1.xy().added(&GpXY::new(d1.x, d1.y).multiplied_scalar(r1)),
                        );
                        s.linsol[0] = GpLin2d::from_pnt_dir(
                            p1,
                            GpDir2d::new(-d1.y, d1.x).expect("gp_Dir2d"),
                        );
                        s.qualifier1[0] = qualified1.qualifier();
                        s.qualifier2[0] = qualified1.qualifier();
                        s.pnttg1sol[0] = p1;
                        s.pnttg2sol[0] = p1;
                        s.well_done = true;
                        s.nbr_sol = 1;
                    } else {
                        let p1 = loc2.rotated(&loc1, ((r1 - r2) / dist).asin());
                        d1 = dir_between(&loc1, &p1);
                        let p1 = GpPnt2d::from_xy(
                            loc1.xy().added(&GpXY::new(r1 * d1.y, -r1 * d1.x)),
                        );
                        s.linsol[0] = GpLin2d::from_pnt_dir(p1, d1);
                        s.qualifier1[0] = qualified1.qualifier();
                        s.qualifier2[0] = qualified1.qualifier();
                        s.pnttg1sol[0] = p1;
                        s.pnttg2sol[0] = GpPnt2d::from_xy(loc2.xy().added(
                            &GpXY::new(r2 * d1.y, -r2 * d1.x).multiplied_scalar(signe),
                        ));
                        s.well_done = true;
                        s.nbr_sol = 1;
                        if rmin + rmax < dist {
                            let p2 = loc2.rotated(&loc1, (signe * (r1 + r2) / dist).asin());
                            let d2 = dir_between(&loc1, &p2);
                            let p2 = GpPnt2d::from_xy(loc1.xy().added(
                                &GpXY::new(r1 * d2.y, -r1 * d2.x).multiplied_scalar(signe),
                            ));
                            s.linsol[1] = GpLin2d::from_pnt_dir(p2, d2);
                            s.qualifier1[0] = qualified1.qualifier();
                            s.qualifier2[0] = qualified1.qualifier();
                            // `cxx:425`: OCCT stores `P1`, not `P2`, here.
                            s.pnttg1sol[1] = p1;
                            s.pnttg2sol[1] = GpPnt2d::from_xy(
                                loc2.xy().added(&GpXY::new(-r2 * d2.y, r2 * d2.x)),
                            );
                            s.nbr_sol = 2;
                        }
                    }
                } else if (qualified1.is_unqualified() && qualified2.is_outside())
                    || (qualified2.is_unqualified() && qualified1.is_outside())
                {
                    // `cxx:431-496`.
                    signe = if qualified2.is_unqualified() { 1.0 } else { -1.0 };
                    if (dist + rmin - rmax).abs() <= tol && dist >= tol {
                        d1 = if r1 < r2 {
                            dir_between(&loc2, &loc1)
                        } else {
                            dir_between(&loc1, &loc2)
                        };
                        let p1 = GpPnt2d::from_xy(
                            loc1.xy().added(&GpXY::new(d1.x, d1.y).multiplied_scalar(r1)),
                        );
                        s.linsol[0] =
                            GpLin2d::from_pnt_dir(p1, GpDir2d::new(d1.y, -d1.x).expect("gp_Dir2d"));
                        s.qualifier1[0] = qualified1.qualifier();
                        s.qualifier2[0] = qualified1.qualifier();
                        s.pnttg1sol[0] = p1;
                        s.pnttg2sol[0] = p1;
                        s.well_done = true;
                        s.nbr_sol = 1;
                    } else {
                        let p1 =
                            loc2.rotated(&loc1, (signe * (r2 - r1) / dist).asin());
                        d1 = dir_between(&loc1, &p1);
                        let p1 = GpPnt2d::from_xy(
                            loc1.xy().added(&GpXY::new(-r1 * d1.y, r1 * d1.x)),
                        );
                        s.linsol[0] = GpLin2d::from_pnt_dir(p1, d1);
                        s.qualifier1[0] = qualified1.qualifier();
                        s.qualifier2[0] = qualified1.qualifier();
                        s.pnttg1sol[0] = p1;
                        s.pnttg2sol[0] = GpPnt2d::from_xy(loc2.xy().added(
                            &GpXY::new(-r2 * d1.y, r2 * d1.x).multiplied_scalar(signe),
                        ));
                        s.well_done = true;
                        s.nbr_sol = 1;
                        if rmin + rmax < dist {
                            let p2 = loc2
                                .rotated(&loc1, (signe * (-r2 - r1) / dist).asin());
                            let d2 = dir_between(&loc1, &p2);
                            let p2 = GpPnt2d::from_xy(loc1.xy().added(
                                &GpXY::new(-r1 * d2.y, r1 * d2.x).multiplied_scalar(signe),
                            ));
                            s.linsol[1] = GpLin2d::from_pnt_dir(p2, d2);
                            s.qualifier1[0] = qualified1.qualifier();
                            s.qualifier2[0] = qualified1.qualifier();
                            s.pnttg1sol[1] = p2;
                            s.pnttg2sol[1] = GpPnt2d::from_xy(loc2.xy().added(
                                &GpXY::new(r2 * d2.y, -r2 * d2.x).multiplied_scalar(signe),
                            ));
                            s.nbr_sol = 2;
                        }
                    }
                }
            }
        } else {
            // `cxx:499-531`: both arguments unqualified.
            for _ in 1..=2 {
                signe = -signe;
                s.nbr_sol += 1;
                let i = s.nbr_sol - 1;
                let p1 = loc2.rotated(&loc1, (signe * (r2 - r1) / dist).asin());
                d1 = dir_between(&loc1, &p1);
                let p1 = GpPnt2d::from_xy(loc1.xy().added(
                    &GpXY::new(-r1 * d1.y, r1 * d1.x).multiplied_scalar(signe),
                ));
                s.linsol[i] = GpLin2d::from_pnt_dir(p1, d1);
                s.qualifier1[i] = qualified1.qualifier();
                s.qualifier2[i] = qualified1.qualifier();
                s.pnttg1sol[i] = p1;
                s.pnttg2sol[i] = GpPnt2d::from_xy(loc2.xy().added(
                    &GpXY::new(-r2 * d1.y, r2 * d1.x).multiplied_scalar(signe),
                ));
                s.well_done = true;
                if rmin + rmax < dist {
                    let p2 = loc2.rotated(&loc1, (signe * (r2 + r1) / dist).asin());
                    let d2 = dir_between(&loc1, &p2);
                    let p2 = GpPnt2d::from_xy(loc1.xy().added(
                        &GpXY::new(r1 * d2.y, -r1 * d2.x).multiplied_scalar(signe),
                    ));
                    s.nbr_sol += 1;
                    let j = s.nbr_sol - 1;
                    s.linsol[j] = GpLin2d::from_pnt_dir(p2, d2);
                    s.qualifier1[j] = qualified1.qualifier();
                    s.qualifier2[j] = qualified1.qualifier();
                    s.pnttg1sol[j] = p2;
                    s.pnttg2sol[j] = GpPnt2d::from_xy(loc2.xy().added(
                        &GpXY::new(-r2 * d2.y, r2 * d2.x).multiplied_scalar(signe),
                    ));
                }
            }
        }

        // `cxx:533-539`.
        for i in 0..s.nbr_sol {
            s.par1sol[i] = parameter_lin2d(&s.linsol[i], &s.pnttg1sol[i]);
            s.par2sol[i] = parameter_lin2d(&s.linsol[i], &s.pnttg2sol[i]);
            s.pararg1[i] = parameter_circ2d(&c1, &s.pnttg1sol[i]);
            s.pararg2[i] = parameter_circ2d(&c2, &s.pnttg2sol[i]);
        }
        s
    }

    /// `IsDone()` (`cxx:542-545`).
    pub fn is_done(&self) -> bool {
        self.well_done
    }

    /// `NbSolutions()` (`cxx:547-550`): no `WellDone` test in OCCT.
    pub fn nb_solutions(&self) -> usize {
        self.nbr_sol
    }

    /// `ThisSolution(Index)` (`cxx:552-560`): only the index is checked.
    pub fn this_solution(&self, index: usize) -> GpLin2d {
        self.check_index(index);
        self.linsol[index - 1]
    }

    /// `WhichQualifier(Index, Qualif1, Qualif2)` (`cxx:562-579`).
    pub fn which_qualifier(&self, index: usize) -> (GccEntPosition, GccEntPosition) {
        self.check_done();
        self.check_index(index);
        (self.qualifier1[index - 1], self.qualifier2[index - 1])
    }

    /// `Tangency1(Index, ParSol, ParArg, PntSol)` (`cxx:581-600`).
    pub fn tangency1(&self, index: usize) -> (f64, f64, GpPnt2d) {
        self.check_done();
        self.check_index(index);
        (self.par1sol[index - 1], self.pararg1[index - 1], self.pnttg1sol[index - 1])
    }

    /// `Tangency2(Index, ParSol, ParArg, PntSol)` (`cxx:602-621`).
    pub fn tangency2(&self, index: usize) -> (f64, f64, GpPnt2d) {
        self.check_done();
        self.check_index(index);
        (self.par2sol[index - 1], self.pararg2[index - 1], self.pnttg2sol[index - 1])
    }

    /// `if (!WellDone) throw StdFail_NotDone()`.
    fn check_done(&self) {
        assert!(self.well_done, "StdFail_NotDone in GccAna_Lin2d2Tan");
    }

    /// `if (Index <= 0 || Index > NbrSol) throw Standard_OutOfRange()`.
    fn check_index(&self, index: usize) {
        assert!(
            index >= 1 && index <= self.nbr_sol,
            "Standard_OutOfRange in GccAna_Lin2d2Tan"
        );
    }
}

/// `gp_Dir2d(gp_Vec2d(from, to))`.
fn dir_between(from: &GpPnt2d, to: &GpPnt2d) -> GpDir2d {
    GpDir2d::new(to.x() - from.x(), to.y() - from.y())
        .expect("gp_VectorWithNullMagnitude in GccAna_Lin2d2Tan")
}
