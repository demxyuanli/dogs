//! `GeomFill_SectionPlacement::Perform(Path, Tol)` — the placement search
//! (`GeomFill_SectionPlacement.cxx:400-704`).
//!
//! `PonPath` is a scratch point in OCCT: the only reads are the
//! `myExt.Perform(PonPath)` / `PExt.Perform(PonSec)` inputs (`cxx:426`, `cxx:606`,
//! `cxx:738`) and the `V1(PonPath, TheAxe.Location())` plane distance
//! (`cxx:446`). Every write in the planar branch (`cxx:461`, `cxx:485`,
//! `cxx:668`, `cxx:694`) and the `PonSec` refresh next to it (`cxx:511-519`) is
//! dead in OCCT because those branches are exclusive with the only consumers.
//! The port keeps the live reads and drops the dead writes; the OCCT line
//! numbers are noted where they are dropped.

use std::f64::consts::FRAC_PI_2;

use occt_core::gp::{GpAx3, GpPnt, GpVec};
use occt_core::precision::{Precision, SQUARE_CONFUSION};

use crate::curve::Curve;
use crate::extrema_cc::curve_curve_extrema_all_range;
use crate::extrema_pc::{extrema_ext_pc_range_done, ExtPcSolution};
use crate::geom_fill::plane_curve_inter::plane_intersections;

use super::init_section::extrema_on_section;
use super::math_helpers::{choix, eval_angle, tangente};
use super::SectionPlacement;

impl SectionPlacement {
    /// `Perform(const occ::handle<Adaptor3d_Curve>& Path, Tol)` (`cxx:400-704`).
    pub(super) fn perform_path(&mut self, path: &dyn Curve, tol: f64) {
        let int_tol = 1.0e-5;
        let mut dist_center = Precision::INFINITE;

        if self.is_point {
            // `Extrema_ExtPC Projector(myPoint, *Path, Precision::Confusion())`
            let (sols, done) = extrema_ext_pc_range_done(
                path,
                &self.point,
                path.first_parameter(),
                path.last_parameter(),
            );
            dist_mini(
                path,
                &self.point,
                &sols,
                done,
                &mut self.dist,
                &mut self.path_param,
            );
            self.angle_max = FRAC_PI_2;
        } else {
            let section = self
                .section_curve
                .clone()
                .expect("GeomFill_SectionPlacement: null section");
            let sec = section.as_ref();

            self.path_param = path.first_parameter();
            self.sec_param = sec.first_parameter();

            let v_ref = GpVec::from_xyz(self.the_axe.direction().xyz());

            let (pon_path, dp1) = tangente(path, self.path_param);
            let pon_sec = sec.value(self.sec_param);
            self.dist = pon_path.distance(&pon_sec);
            if self.dist > tol {
                // On Cherche un meilleur point sur la section
                let (sols, done) = extrema_on_section(sec, &pon_path);
                if done {
                    dist_mini(
                        sec,
                        &pon_path,
                        &sols,
                        done,
                        &mut self.dist,
                        &mut self.sec_param,
                    );
                }
            }
            self.angle_max = eval_angle(&v_ref, &dp1);
            if self.is_plan {
                self.angle_max = FRAC_PI_2 - self.angle_max;
            }

            // OCCT declares `Trouve` before the `if (isplan)` block
            // (`cxx:439`) but only reads/writes it in the general branch.
            let mut trouve;

            // `PonSec` is only consumed by the `!isplan` fallback (`cxx:683`),
            // so it stays at its initial value here.
            let mut pon_sec = pon_sec;

            if self.is_plan {
                // (1.1) Distances Point-Plan
                let v1 = GpVec::from_pnts(&pon_path, self.the_axe.location());
                let dist_plan = v1.dot(&v_ref).abs();
                if dist_plan <= int_tol {
                    dist_center = v1.magnitude();
                }

                let plast = path.value(path.last_parameter());
                let v1 = GpVec::from_pnts(&plast, self.the_axe.location());
                let dist_plan = v1.dot(&v_ref).abs();
                if dist_plan <= int_tol {
                    let a_dist = v1.magnitude();
                    if a_dist < dist_center {
                        dist_center = a_dist;
                        // `PonPath = Plast;` (`cxx:461`) is dead, see the module
                        // comment.
                        self.path_param = path.last_parameter();
                    }
                }

                // (1.2) Intersection Plan-courbe
                let axe = GpAx3::from_point_dir(
                    *self.the_axe.location(),
                    *self.the_axe.direction(),
                )
                .expect("gp_Ax3: null direction");
                let axe_pnt = axe.location();
                let normal = GpVec::from_xyz(axe.direction().xyz());
                let pnts = plane_intersections(path, &axe_pnt, &normal);
                for (w, pw) in &pnts {
                    let a_dist = pw.distance(&axe_pnt);
                    if a_dist < dist_center {
                        dist_center = a_dist;
                        // `PonPath = P;` (`cxx:485`) is dead.
                        self.path_param = *w;
                    }
                }

                // `if (!Intersector.IsDone() || Intersector.NbPoints() == 0)`
                if pnts.is_empty() {
                    // Comparing the distances from the path's endpoints to the
                    // best matching plane of the profile.
                    let first_point = path.value(path.first_parameter());
                    let last_point = path.value(path.last_parameter());
                    let plane_square_distance = |p: &GpPnt| {
                        let d = GpVec::from_pnts(&axe_pnt, p).dot(&normal);
                        d * d
                    };
                    let first_distance = plane_square_distance(&first_point);
                    let last_distance = plane_square_distance(&last_point);

                    if ((first_distance.abs() < SQUARE_CONFUSION)
                        && (last_distance.abs() < SQUARE_CONFUSION))
                        || first_distance < last_distance
                    {
                        self.path_param = path.first_parameter();
                    } else {
                        self.path_param = path.last_parameter();
                        // `Tangente(*Path, PathParam, PonPath, dp1)` and the
                        // `PonSec`/`Dist` refresh (`cxx:510-519`) are dead.
                        let (_p, dp1) = tangente(path, self.path_param);
                        self.angle_max = eval_angle(&v_ref, &dp1);
                        self.angle_max = FRAC_PI_2 - self.angle_max;
                    }
                }
            }

            // Cas General
            if !self.is_plan {
                // (2.1) Distance avec les extremites ...
                let (sols, done) = extrema_on_section(sec, &pon_path);
                if done {
                    let (mut distaux, mut taux) = (0.0, 0.0);
                    dist_mini(sec, &pon_path, &sols, done, &mut distaux, &mut taux);
                    if distaux < self.dist {
                        self.dist = distaux;
                        self.sec_param = taux;
                    }
                }
                trouve = self.dist <= tol;
                if !trouve {
                    let (p_end, dp1) = tangente(path, path.last_parameter());
                    let alpha = eval_angle(&v_ref, &dp1);
                    let (sols, done) = extrema_on_section(sec, &p_end);
                    if done {
                        let (mut distaux, mut taux) = (0.0, 0.0);
                        dist_mini(sec, &p_end, &sols, done, &mut distaux, &mut taux);
                        if choix(distaux, alpha, self.dist, self.angle_max, self.gabarit) {
                            self.dist = distaux;
                            self.sec_param = taux;
                            self.angle_max = alpha;
                            // `PonPath = P;` (`cxx:635`) is dead.
                            self.path_param = path.last_parameter();
                        }
                    }
                    trouve = self.dist <= tol;
                }

                // (2.2) Distance courbe-courbe
                if !trouve {
                    // `Extrema_ExtCC Ext(*Path, myAdpSection, ...)`
                    //
                    // `Extrema_ExtCC::IsDone() && !IsParallel()` has no
                    // separate flag in the port: the parallel arm reports no
                    // solution, so an empty sequence is the same skip.
                    let pairs = curve_curve_extrema_all_range(
                        path,
                        sec,
                        path.first_parameter(),
                        path.last_parameter(),
                        sec.first_parameter(),
                        sec.last_parameter(),
                    );
                    for pair in pairs {
                        let distaux = pair.distance;
                        let (_, dp1) = tangente(path, pair.u1);
                        let alpha = eval_angle(&v_ref, &dp1);
                        if choix(distaux, alpha, self.dist, self.angle_max, self.gabarit) {
                            trouve = true;
                            self.dist = distaux;
                            self.path_param = pair.u1;
                            self.sec_param = pair.u2;
                            pon_sec = pair.p2;
                            // `PonPath = P;` (`cxx:668`) is dead.
                            self.angle_max = alpha;
                        }
                    }
                    if !trouve {
                        // Si l'on a toujours rien, on essai une distance
                        // point/path, c'est la derniere chance.
                        // `Extrema_ExtPC PExt; PExt.Initialize(*Path, ...)`
                        let (sols, done) = extrema_ext_pc_range_done(
                            path,
                            &pon_sec,
                            path.first_parameter(),
                            path.last_parameter(),
                        );
                        if done {
                            // OCC13595: `DistMini(PExt, *Path, ...)`, i.e. the
                            // endpoints come from the path, not the section.
                            let (mut distaux, mut taux) = (0.0, 0.0);
                            dist_mini(path, &pon_sec, &sols, done, &mut distaux, &mut taux);
                            let (_p, dp1) = tangente(path, taux);
                            let alpha = eval_angle(&v_ref, &dp1);
                            if choix(distaux, alpha, self.dist, self.angle_max, self.gabarit) {
                                self.dist = distaux;
                                // `PonPath = P;` (`cxx:694`) is dead.
                                self.angle_max = alpha;
                                self.path_param = taux;
                            }
                        }
                    }
                }
            }
        }

        self.done = true;
    }
}

/// `DistMini(Ext, C, Dist, Param)` (`cxx:114-153`): the endpoint distances of
/// `Extrema_ExtPC::TrimmedSquareDistances` (`Extrema_GGExtPC.hxx:145-156`)
/// first, then every stored solution. OCCT does not filter on `IsMin` here.
pub(super) fn dist_mini(
    curve: &dyn Curve,
    point: &GpPnt,
    sols: &[ExtPcSolution],
    done: bool,
    dist: &mut f64,
    param: &mut f64,
) {
    // `RealLast()`
    let real_last = f64::MAX;
    let mut dist2 = real_last;

    let u_inf = curve.first_parameter();
    let u_sup = curve.last_parameter();
    let dist1 = if Precision::is_infinite(u_inf) {
        real_last
    } else {
        point.square_distance(&curve.value(u_inf))
    };
    let dist2_end = if Precision::is_infinite(u_sup) {
        real_last
    } else {
        point.square_distance(&curve.value(u_sup))
    };

    if dist1 < dist2 || dist2_end < dist2 {
        if dist1 < dist2_end {
            dist2 = dist1;
            *param = u_inf;
        } else {
            dist2 = dist2_end;
            *param = u_sup;
        }
    }

    if done {
        for s in sols {
            if s.sq_dist < dist2 {
                dist2 = s.sq_dist;
                *param = s.u;
            }
        }
    }
    *dist = dist2.sqrt();
}
