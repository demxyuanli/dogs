//! `AppParCurves_LeastSquare::Init` (`AppParCurves_LeastSquare.gxx:275-508`).
//!
//! Resolves the effective constraints, flattens the input points into
//! `mypoints`, seeds the constrained poles and sizes the contact system.

use occt_core::gp::{GpPnt, GpPnt2d};

use super::*;

use super::affect::affect;

impl LeastSquare {
    /// `Init(SSP, FirstPoint, LastPoint)` (`gxx:275-508`).
    pub(super) fn init<L: LineTool>(&mut self, ml: &L, first_point: i32, last_point: i32) {
        // Variable de controle
        self.iscalculated = false;
        self.isready = true;

        self.myfirstp = first_point;
        self.mylastp = last_point;
        self.first_p = Self::the_first_point(self.first_constraint, self.myfirstp);
        self.last_p = Self::the_last_point(self.last_constraint, self.mylastp);

        self.nbp2d = ml.nb_p2d();
        self.nbp = ml.nb_p3d();
        let nb_p2d = self.nbp2d;
        let nb_p = self.nbp;
        let mut mynb_p2d = nb_p2d;
        let mut mynb_p = nb_p;
        if nb_p2d == 0 {
            mynb_p2d = 1;
        }
        if nb_p == 0 {
            mynb_p = 1;
        }

        self.deg = self.nbpoles - 1;

        if let Some(mults) = self.mymults.clone() {
            let mut sum = 0;
            for m in &mults {
                sum += *m;
            }
            self.deg = sum - self.nbpoles - 1;
            let knots = self.myknots.clone().unwrap_or_default();
            let mut k = 1;
            for i in 0..knots.len() {
                for _ in 1..=mults[i] {
                    self.vflatknots.set_value(k as usize, knots[i]);
                    k += 1;
                }
            }
        }

        let mut first_cons = self.first_constraint;
        let mut last_cons = self.last_constraint;
        affect(
            ml,
            first_point,
            &mut first_cons,
            &mut self.vec1t,
            &mut self.vec1c,
            nb_p,
            nb_p2d,
        );
        affect(
            ml,
            last_point,
            &mut last_cons,
            &mut self.vec2t,
            &mut self.vec2c,
            nb_p,
            nb_p2d,
        );
        self.first_constraint = first_cons;
        self.last_constraint = last_cons;

        for j in self.myfirstp..=self.mylastp {
            let mut i2 = 1;
            let (tab_p, tab_p2d) = if nb_p != 0 && nb_p2d != 0 {
                ml.value(j)
            } else if nb_p2d != 0 {
                (Vec::new(), ml.value_2d(j))
            } else {
                (ml.value_3d(j), Vec::new())
            };
            for i in 1..=nb_p {
                let p = tab_p[(i - 1) as usize];
                self.mypoints.set_value(j as usize, i2 as usize, p.x());
                self.mypoints.set_value(j as usize, (i2 + 1) as usize, p.y());
                self.mypoints.set_value(j as usize, (i2 + 2) as usize, p.z());
                i2 += 3;
            }
            for i in 1..=nb_p2d {
                let p = tab_p2d[(i - 1) as usize];
                self.mypoints.set_value(j as usize, i2 as usize, p.x());
                self.mypoints.set_value(j as usize, (i2 + 1) as usize, p.y());
                i2 += 2;
            }
        }

        let mut pole1 = MultiPoint::with_size(mynb_p.max(0) as usize, mynb_p2d.max(0) as usize);
        let mut pole_n = MultiPoint::with_size(mynb_p.max(0) as usize, mynb_p2d.max(0) as usize);

        if at_least_pass(self.first_constraint) {
            let mut i2 = 1;
            for i in 1..=nb_p {
                let poi = GpPnt::new(
                    self.mypoints.value(self.myfirstp as usize, i2 as usize),
                    self.mypoints.value(self.myfirstp as usize, (i2 + 1) as usize),
                    self.mypoints.value(self.myfirstp as usize, (i2 + 2) as usize),
                );
                pole1.set_point(i, &poi);
                i2 += 3;
            }
            for i in 1..=nb_p2d {
                let poi2d = GpPnt2d::new(
                    self.mypoints.value(self.myfirstp as usize, i2 as usize),
                    self.mypoints.value(self.myfirstp as usize, (i2 + 1) as usize),
                );
                pole1.set_point2d(i + nb_p, &poi2d);
                i2 += 2;
            }
            for i in 1..=self.mypoles.col_count() as i32 {
                let v = self.mypoints.value(self.myfirstp as usize, i as usize);
                self.mypoles.set_value(1, i as usize, v);
            }
        }

        // `gxx:376`: the third test reads `FirstConstraint` where the
        // surrounding arms use `LastConstraint`. Kept verbatim.
        if self.last_constraint == Constraint::PassPoint
            || self.last_constraint == Constraint::TangencyPoint
            || self.first_constraint == Constraint::CurvaturePoint
        {
            let mut i2 = 1;
            for i in 1..=nb_p {
                let poi = GpPnt::new(
                    self.mypoints.value(self.mylastp as usize, i2 as usize),
                    self.mypoints.value(self.mylastp as usize, (i2 + 1) as usize),
                    self.mypoints.value(self.mylastp as usize, (i2 + 2) as usize),
                );
                pole_n.set_point(i, &poi);
                i2 += 3;
            }
            for i in 1..=nb_p2d {
                let poi2d = GpPnt2d::new(
                    self.mypoints.value(self.mylastp as usize, i2 as usize),
                    self.mypoints.value(self.mylastp as usize, (i2 + 1) as usize),
                );
                pole_n.set_point2d(i + nb_p, &poi2d);
                i2 += 2;
            }
            for i in 1..=self.mypoles.col_count() as i32 {
                let v = self.mypoints.value(self.mylastp as usize, i as usize);
                self.mypoles.set_value(self.nbpoles as usize, i as usize, v);
            }
        }

        if self.first_constraint == Constraint::NoConstraint {
            self.resinit = 1;
            self.scu.base_mut().set_value(1, &pole1);
            if self.last_constraint == Constraint::NoConstraint {
                self.resfin = self.nbpoles;
            } else if self.last_constraint == Constraint::PassPoint {
                self.resfin = self.nbpoles - 1;
                self.scu.base_mut().set_value(self.nbpoles, &pole_n);
            } else if self.last_constraint == Constraint::TangencyPoint {
                self.resfin = self.nbpoles - 2;
                self.scu.base_mut().set_value(self.nbpoles, &pole_n);
            } else if self.last_constraint == Constraint::CurvaturePoint {
                self.resfin = self.nbpoles - 3;
                self.scu.base_mut().set_value(self.nbpoles, &pole_n);
            }
        } else if self.first_constraint == Constraint::PassPoint {
            self.resinit = 2;
            self.scu.base_mut().set_value(1, &pole1);
            if self.last_constraint == Constraint::NoConstraint {
                self.resfin = self.nbpoles;
            } else if self.last_constraint == Constraint::PassPoint {
                self.resfin = self.nbpoles - 1;
                self.scu.base_mut().set_value(self.nbpoles, &pole_n);
            } else if self.last_constraint == Constraint::TangencyPoint {
                self.resfin = self.nbpoles - 2;
                self.scu.base_mut().set_value(self.nbpoles, &pole_n);
            } else if self.last_constraint == Constraint::CurvaturePoint {
                self.resfin = self.nbpoles - 3;
                self.scu.base_mut().set_value(self.nbpoles, &pole_n);
            }
        } else if self.first_constraint == Constraint::TangencyPoint {
            self.resinit = 3;
            self.scu.base_mut().set_value(1, &pole1);
            if self.last_constraint == Constraint::NoConstraint {
                self.resfin = self.nbpoles;
            }
            if self.last_constraint == Constraint::PassPoint {
                self.resfin = self.nbpoles - 1;
                self.scu.base_mut().set_value(self.nbpoles, &pole_n);
            }
            if self.last_constraint == Constraint::TangencyPoint {
                self.resfin = self.nbpoles - 2;
                self.scu.base_mut().set_value(self.nbpoles, &pole_n);
            } else if self.last_constraint == Constraint::CurvaturePoint {
                self.resfin = self.nbpoles - 3;
                self.scu.base_mut().set_value(self.nbpoles, &pole_n);
            }
        } else if self.first_constraint == Constraint::CurvaturePoint {
            self.resinit = 4;
            self.scu.base_mut().set_value(1, &pole1);
            if self.last_constraint == Constraint::NoConstraint {
                self.resfin = self.nbpoles;
            }
            if self.last_constraint == Constraint::PassPoint {
                self.resfin = self.nbpoles - 1;
                self.scu.base_mut().set_value(self.nbpoles, &pole_n);
            }
            if self.last_constraint == Constraint::TangencyPoint {
                self.resfin = self.nbpoles - 2;
                self.scu.base_mut().set_value(self.nbpoles, &pole_n);
            } else if self.last_constraint == Constraint::CurvaturePoint {
                self.resfin = self.nbpoles - 3;
                self.scu.base_mut().set_value(self.nbpoles, &pole_n);
            }
        }

        let nincx = self.resfin - self.resinit + 1;
        if nincx < 1 {
            // Impossible d'aller plus loin
            self.isready = false;
            return;
        }
        let neq = self.last_p - self.first_p + 1;

        self.na = 3 * nb_p + 2 * nb_p2d;
        self.nlignes = self.na * neq;
        self.ninc = self.na * nincx;
        if at_least_tangency(self.first_constraint) {
            self.ninc += 1;
        }
        if at_least_tangency(self.last_constraint) {
            self.ninc += 1;
        }
    }
}



