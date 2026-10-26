//! `AppDef_Variational::EstTangent` / `EstSecnd` and the `Pnt` view they share.
//!
//! Sources: `AppDef_Variational.cxx:2401-2581` (`EstTangent`),
//! `cxx:2585-2679` (`EstSecnd`).

use occt_math::vector::MathVector;

use super::{Variational, CONFUSION};

impl Variational {
    /// `Pnt` view over `myTabPoints` starting at 1-based `adr`
    /// (`cxx:2422-2424`: `math_Vector(Pnt, 1, myDimension)`).
    pub(super) fn pnt_at(&self, adr: usize) -> MathVector {
        let dim = self.dimension as usize;
        MathVector::from_slice(&self.tab_points[adr..adr + dim])
    }

    /// `EstTangent(ipnt, VTang)` (`cxx:2401-2581`).
    pub(super) fn est_tangent(&self, ipnt: i32, vtang: &mut MathVector) {
        const EPS_NORM: f64 = 1.0e-9;
        let eps1 = CONFUSION * 0.01;
        let dim = self.dimension as i32;
        let mut wpnt = 1.0;

        if ipnt == self.first_point {
            if self.nb_points < 3 {
                wpnt = 0.0;
            } else {
                let adr1 = 1usize;
                let adr2 = adr1 + self.dimension as usize;
                let adr3 = adr2 + self.dimension as usize;
                let pnt1 = self.pnt_at(adr1);
                let pnt2 = self.pnt_at(adr2);
                let pnt3 = self.pnt_at(adr3);

                let v1 = pnt2.subtracted(&pnt1).norm();
                let mut v2 = 0.0;
                if v1 > eps1 {
                    v2 = pnt3.subtracted(&pnt2).norm();
                }
                if v2 > eps1 {
                    let mut d = v1 / (v1 + v2);
                    let d1 = 1.0 / (d * (1.0 - d));
                    d *= d;
                    let mut val = pnt1.multiplied_scalar(d - 1.0);
                    val.add(&pnt2);
                    val.subtract(&pnt3.multiplied_scalar(d));
                    val.multiply_scalar(d1);
                    *vtang = val;
                } else {
                    *vtang = pnt2.subtracted(&pnt1);
                }
            }
        } else if ipnt == self.last_point {
            if self.nb_points < 3 {
                wpnt = 0.0;
            } else {
                let adr1 = ((self.last_point - 3) * self.dimension + 1) as usize;
                let adr2 = adr1 + self.dimension as usize;
                let adr3 = adr2 + self.dimension as usize;
                let pnt1 = self.pnt_at(adr1);
                let pnt2 = self.pnt_at(adr2);
                let pnt3 = self.pnt_at(adr3);

                let v1 = pnt2.subtracted(&pnt1).norm();
                let mut v2 = 0.0;
                if v1 > eps1 {
                    v2 = pnt3.subtracted(&pnt2).norm();
                }
                if v2 > eps1 {
                    let mut d = v1 / (v1 + v2);
                    let d1 = 1.0 / (d * (1.0 - d));
                    d *= d - 2.0;
                    let mut val = pnt1.multiplied_scalar(d + 1.0);
                    val.subtract(&pnt2);
                    val.subtract(&pnt3.multiplied_scalar(d));
                    val.multiply_scalar(d1);
                    *vtang = val;
                } else {
                    *vtang = pnt3.subtracted(&pnt2);
                }
            }
        } else {
            let adr1 = ((ipnt - self.first_point - 1) * self.dimension + 1) as usize;
            let adr2 = adr1 + 2 * self.dimension as usize;
            let pnt1 = self.pnt_at(adr1);
            let pnt2 = self.pnt_at(adr2);
            *vtang = pnt2.subtracted(&pnt1);
        }

        let mut vnorm = vtang.norm();
        if vnorm <= EPS_NORM {
            vtang.init(0.0);
        } else {
            vtang.divide_scalar(vnorm);
        }

        // Estimation with constraints (`cxx:2514-2579`).
        let mut wcnt = 0.0;
        let mut idcnt = 1i32;
        let nb_constr = self.nb_pass_points + self.nb_tang_points + self.nb_curv_points;
        let mut vcnt = MathVector::with_init(1, dim as usize, 0.0);

        if nb_constr > 0 {
            while idcnt <= nb_constr
                && self.typ_constraints[(2 * idcnt - 1) as usize] < ipnt
            {
                idcnt += 1;
            }
            if idcnt <= nb_constr
                && self.typ_constraints[(2 * idcnt - 1) as usize] == ipnt
                && self.typ_constraints[(2 * idcnt) as usize] >= 1
            {
                wcnt = 1.0;
                let mut i0 = (2 * self.dimension * (idcnt - 1)) as usize;
                let mut k = 0usize;
                for _i in 1..=self.nb_p3d {
                    for _j in 1..=3 {
                        i0 += 1;
                        k += 1;
                        vcnt.set_value(k, self.tab_constraints[i0]);
                    }
                    i0 += 3;
                }
                for _i in 1..=self.nb_p2d {
                    for _j in 1..=2 {
                        i0 += 1;
                        k += 1;
                        vcnt.set_value(k, self.tab_constraints[i0]);
                    }
                    i0 += 2;
                }
            }
        }

        let denom = wpnt + wcnt;
        let denom = if denom == 0.0 { 1.0 } else { 1.0 / denom };

        let mut val = vtang.multiplied_scalar(wpnt);
        val.add(&vcnt.multiplied_scalar(wcnt));
        val.multiply_scalar(denom);
        *vtang = val;

        vnorm = vtang.norm();
        if vnorm <= EPS_NORM {
            vtang.init(0.0);
        } else {
            vtang.divide_scalar(vnorm);
        }
    }

    /// `EstSecnd(ipnt, VTang1, VTang2, Length, VScnd)` (`cxx:2585-2679`).
    pub(super) fn est_secnd(
        &self,
        ipnt: i32,
        vtang1: &MathVector,
        vtang2: &MathVector,
        length: f64,
        vscnd: &mut MathVector,
    ) {
        const EPS: f64 = 1.0e-9;
        let dim = self.dimension as i32;
        let wpnt = 1.0;

        let aux = if ipnt == self.first_point {
            self.param(ipnt + 1) - self.param(ipnt)
        } else if ipnt == self.last_point {
            self.param(ipnt) - self.param(ipnt - 1)
        } else {
            self.param(ipnt + 1) - self.param(ipnt - 1)
        };
        let aux = if aux <= EPS { 1.0 } else { 1.0 / aux };

        let mut val = vtang2.subtracted(vtang1);
        val.multiply_scalar(aux);
        *vscnd = val;

        // Estimation with constraints (`cxx:2622-2678`).
        let mut wcnt = 0.0;
        let mut idcnt = 1i32;
        let nb_constr = self.nb_pass_points + self.nb_tang_points + self.nb_curv_points;
        let mut vcnt = MathVector::with_init(1, dim as usize, 0.0);

        if nb_constr > 0 {
            while idcnt <= nb_constr
                && self.typ_constraints[(2 * idcnt - 1) as usize] < ipnt
            {
                idcnt += 1;
            }
            if idcnt <= nb_constr
                && self.typ_constraints[(2 * idcnt - 1) as usize] == ipnt
                && self.typ_constraints[(2 * idcnt) as usize] >= 2
            {
                wcnt = 1.0;
                let mut i0 = (2 * self.dimension * (idcnt - 1) + 3) as usize;
                let mut k = 0usize;
                for _i in 1..=self.nb_p3d {
                    for _j in 1..=3 {
                        i0 += 1;
                        k += 1;
                        vcnt.set_value(k, self.tab_constraints[i0]);
                    }
                    i0 += 3;
                }
                i0 -= 1;
                for _i in 1..=self.nb_p2d {
                    for _j in 1..=2 {
                        i0 += 1;
                        k += 1;
                        vcnt.set_value(k, self.tab_constraints[i0]);
                    }
                    i0 += 2;
                }
            }
        }

        let denom = wpnt + wcnt;
        let denom = if denom == 0.0 { 1.0 } else { 1.0 / denom };

        let mut val = vscnd.multiplied_scalar(wpnt);
        val.add(&vcnt.multiplied_scalar(wcnt * length));
        val.multiply_scalar(denom);
        *vscnd = val;
    }
}
