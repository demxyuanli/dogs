//! `Bisector_BisecCC::Perform` (`Bisector_BisecCC.cxx:92-327`), plus the three
//! helpers that only the construction path uses: `SearchBound`
//! (`cxx:1741-1781`), `SupLastParameter` (`cxx:1485-1511`) and
//! `ComputePointEnd` (`cxx:1879-1927`).

use std::sync::Arc;

use occt_core::gp::GpPnt2d;
use occt_core::precision::{CONFUSION, INFINITE, PCONFUSION};

use crate::curve::Curve2d;

use super::super::{is_convex, BisectorPointOnBis};
use super::support::{curvature, point_by_int, proj_on_curve, test_extension};
use super::BisectorBisecCC;

/// `10 * Precision::Confusion()`, the `EpsMin` of `Perform` (`cxx:116`).
const EPS_MIN: f64 = 10.0 * CONFUSION;

impl BisectorBisecCC {
    /// `Perform(Cu1, Cu2, Side1, Side2, Origin, DistMax)` (`cxx:92-327`).
    pub fn perform(
        &mut self,
        cu1: &dyn Curve2d,
        cu2: &dyn Curve2d,
        side1: f64,
        side2: f64,
        origin: &GpPnt2d,
        dist_max: f64,
    ) {
        self.is_empty = false;
        self.dist_max = dist_max;

        self.curve1 = Arc::from(cu1.clone_dyn());
        self.curve2 = Arc::from(cu2.clone_dyn());

        self.sign1 = side1;
        self.sign2 = side2;

        self.is_convex1 = is_convex(&*self.curve1, self.sign1);
        self.is_convex2 = is_convex(&*self.curve2, self.sign2);

        let mut u = 0.0;
        let mut uc1 = 0.0;
        let mut uc2 = 0.0;
        let mut dist = 0.0;
        let mut nb_pnts: i32 = 21;
        let mut ya_poly = true;
        let mut ori_in_poly = false;

        // `cxx:120-127`: project the origin on the guide curve.
        let is_proj_done = proj_on_curve(origin, &*self.curve1, &mut u);
        if !is_proj_done {
            self.is_empty = true;
            return;
        }

        let mut p = self.value_by_int(u, &mut uc1, &mut uc2, &mut dist);
        if dist < CONFUSION {
            let a_p1 = self.curve1.value(uc1);
            let a_p2 = self.curve2.value(uc2);
            let dp = a_p1.distance(&p) + a_p2.distance(&p);
            let dorig = a_p1.distance(origin) + a_p2.distance(origin);
            if dp < dorig {
                self.is_empty = true;
                return;
            }
        }

        if dist < INFINITE {
            // `cxx:143-157`: the parameter of the origin point gives a point on
            // the polygon.
            self.my_polygon
                .append(BisectorPointOnBis::with_params(uc1, uc2, u, dist, p));
            self.start_intervals.push(u);
            if p.is_equal(origin, CONFUSION) {
                ori_in_poly = true;
            }
        } else {
            // `cxx:159-179`: the origin is on the extension; find the first
            // point of the polygon by dichotomy.
            let du = (self.curve1.last_parameter() - u) / ((nb_pnts - 1) as f64);
            u += du;
            for _ in 1..nb_pnts {
                p = self.value_by_int(u, &mut uc1, &mut uc2, &mut dist);
                if dist < INFINITE {
                    let u_sol = self.search_bound(u - du, u);
                    p = self.value_by_int(u_sol, &mut uc1, &mut uc2, &mut dist);
                    self.start_intervals.push(u_sol);
                    self.my_polygon
                        .append(BisectorPointOnBis::with_params(uc1, uc2, u_sol, dist, p));
                    break;
                }
                u += du;
            }
        }

        if self.my_polygon.length() != 0 {
            self.sup_last_parameter();
            // `cxx:184-225`: construction of the polygon of the bisector.
            u = self.first_parameter();
            let mut du = self.last_parameter() - u;
            if du < EPS_MIN {
                nb_pnts = 3;
            }
            du /= (nb_pnts - 1) as f64;

            u += du;
            let mut prev_pnt = p;
            for _ in 1..nb_pnts {
                p = self.value_by_int(u, &mut uc1, &mut uc2, &mut dist);
                if dist < INFINITE {
                    if p.distance(&prev_pnt) > CONFUSION {
                        self.my_polygon
                            .append(BisectorPointOnBis::with_params(uc1, uc2, u, dist, p));
                    }
                } else {
                    let u_sol = self.search_bound(u - du, u);
                    p = self.value_by_int(u_sol, &mut uc1, &mut uc2, &mut dist);
                    if let Some(first) = self.end_intervals.first_mut() {
                        *first = u_sol;
                    }
                    if p.distance(&prev_pnt) > CONFUSION {
                        self.my_polygon
                            .append(BisectorPointOnBis::with_params(uc1, uc2, u_sol, dist, p));
                    }
                    break;
                }
                u += du;
                prev_pnt = p;
            }
        } else {
            ya_poly = false;
        }

        self.extension_start = false;
        self.extension_end = false;
        self.point_start = *origin;

        if self.is_convex1 && self.is_convex2 {
            if ya_poly {
                self.point_end = self.my_polygon.last().point();
            }
        } else {
            // `cxx:246-320`: the curves are extended at the start and/or the end
            // when one of them is concave and the curves meet at that end with
            // an opening angle of PI.
            if ori_in_poly {
                self.extension_start = false;
            } else {
                self.extension_start = test_extension(&*self.curve1, &*self.curve2, 1);
            }
            self.extension_end = test_extension(&*self.curve1, &*self.curve2, 2);

            if self.extension_end {
                self.point_end = self.curve1.value(self.curve1.last_parameter());
            } else if ya_poly {
                self.point_end = self.my_polygon.last().point();
            } else {
                self.compute_point_end();
            }

            if ya_poly {
                if self.extension_start {
                    let p1 = self.my_polygon.first().point();
                    let u_first = self.start_first() - self.point_start.distance(&p1);
                    self.start_intervals.insert(0, u_first);
                    let second = self.start_value(2);
                    self.end_intervals.insert(0, second);
                }
                if self.extension_end {
                    let p1 = self.my_polygon.last().point();
                    let u_first = self.end_last();
                    let u_last = u_first + self.point_end.distance(&p1);
                    self.start_intervals.push(u_first);
                    self.end_intervals.push(u_last);
                }
            } else {
                // `cxx:308-311`: no polygon, the bisector is a segment.
                self.start_intervals.push(0.0);
                self.end_intervals
                    .push(self.point_end.distance(&self.point_start));
            }
        }
        if !ya_poly && !self.extension_start && !self.extension_end {
            self.is_empty = true;
        }
        if self.my_polygon.length() <= 2 {
            self.is_empty = true;
        }
    }

    /// `SearchBound(U1, U2)` (`cxx:1741-1781`): the smallest parameter of the
    /// bracket `[U1, U2]` whose `ValueByInt` distance leaves `Precision::Infinite`.
    pub fn search_bound(&self, u1: f64, u2: f64) -> f64 {
        const TOL_PNT: f64 = CONFUSION;
        const TOL_PAR: f64 = PCONFUSION;
        let mut u11 = u1;
        let mut u22 = u2;
        let mut uc1 = 0.0;
        let mut uc2 = 0.0;
        let mut dist1 = 0.0;
        let mut dist2 = 0.0;
        let mut p_bis_prec = self.value_by_int(u11, &mut uc1, &mut uc2, &mut dist1);
        let mut p_bis = self.value_by_int(u22, &mut uc1, &mut uc2, &mut dist2);

        while (u22 - u11) > TOL_PAR
            || (dist1 < INFINITE && dist2 < INFINITE && !p_bis.is_equal(&p_bis_prec, TOL_PNT))
        {
            p_bis_prec = p_bis;
            let u_mid = 0.5 * (u22 + u11);
            let mut dist_mid = 0.0;
            p_bis = self.value_by_int(u_mid, &mut uc1, &mut uc2, &mut dist_mid);
            if (dist1 < INFINITE) == (dist_mid < INFINITE) {
                u11 = u_mid;
                dist1 = dist_mid;
            } else {
                u22 = u_mid;
                dist2 = dist_mid;
            }
        }
        p_bis = self.value_by_int(u11, &mut uc1, &mut uc2, &mut dist1);
        let _ = p_bis;
        if dist1 < INFINITE {
            u11
        } else {
            u22
        }
    }

    /// `SupLastParameter()` (`cxx:1485-1511`).
    pub fn sup_last_parameter(&mut self) {
        self.end_intervals.push(self.curve1.last_parameter());
        // The parameter of curve1 associated with one of the extremities of
        // curve2 depends on the signs (`cxx:1492-1501`).
        let u_on_c2 = if self.sign1 == self.sign2 {
            self.curve2.first_parameter()
        } else {
            self.curve2.last_parameter()
        };
        let mut u_on_c1 = 0.0;
        let mut dist = 0.0;
        let ya_sol = point_by_int(
            &*self.curve2,
            &*self.curve1,
            self.sign2,
            self.sign1,
            u_on_c2,
            &mut u_on_c1,
            &mut dist,
        );
        if ya_sol && u_on_c1 > self.start_first() && u_on_c1 < self.end_last() {
            if let Some(first) = self.end_intervals.first_mut() {
                *first = u_on_c1;
            }
        }
    }

    /// `ComputePointEnd()` (`cxx:1879-1927`).
    pub fn compute_point_end(&mut self) {
        let u1 = self.curve1.first_parameter();
        let u2 = if self.sign1 == self.sign2 {
            self.curve2.last_parameter()
        } else {
            self.curve2.first_parameter()
        };
        let k1 = curvature(&*self.curve1, u1, CONFUSION);
        let k2 = curvature(&*self.curve2, u2, CONFUSION);
        let kc = if !self.is_convex1 && !self.is_convex2 {
            if k1 < k2 {
                k1
            } else {
                k2
            }
        } else if !self.is_convex1 {
            k1
        } else {
            k2
        };

        let (pf, mut tf) = self.curve1.d1(u1);
        let _ = tf.normalize();
        let rc = if kc != 0.0 { (1.0 / kc).abs() } else { INFINITE };
        self.point_end = GpPnt2d::new(
            pf.x() - self.sign1 * rc * tf.y(),
            pf.y() + self.sign1 * rc * tf.x(),
        );
    }
}
