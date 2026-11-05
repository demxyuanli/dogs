//! `FairCurve_Energy` (`FairCurve_Energy.hxx:28-130`, `.cxx:28-530`,
//! `.lxx:17-34`): the variable packing, pole reconstruction and the
//! gradient/Hessian unpacking shared by the fair-curve objectives.
//!
//! Rust note: OCCT derives `FairCurve_EnergyOfBatten` from this class and
//! implements the pure virtual `Compute(DerivativeOrder, Result)`
//! (`.hxx:113`). The port keeps `FairCurve_Energy` concrete and passes that
//! `Compute` step in as the `compute` argument of [`Energy::value`] /
//! [`Energy::values`] / [`Energy::values_hessian`], so the arithmetic and the
//! call order of `.cxx:60-166` are unchanged.

use std::rc::Rc;

use occt_core::gp::pnt2d::GpPnt2d;
use occt_core::gp::xy::GpXY;
use occt_math::{MathMatrix, MathVector};

use super::distribution_of_energy::{pole, SharedPoles};

/// `Compute(DerivativeOrder, Result)` (`.hxx:113`): the sub-class step that
/// fills the packed energy/gradient/Hessian vector.
pub type EnergyCompute<'a> = &'a mut dyn FnMut(i32, &mut MathVector) -> bool;

/// `FairCurve_Energy` (`FairCurve_Energy.hxx:28-130`).
#[derive(Debug, Clone)]
pub struct Energy {
    /// `MyPoles` (`.hxx:116`): the same handle the distribution sub-objects
    /// read, so `ComputePoles` is visible to them (`.cxx:451-509`).
    poles: SharedPoles,
    /// `MyContrOrder1`.
    contr_order1: i32,
    /// `MyContrOrder2`.
    contr_order2: i32,
    /// `MyWithAuxValue`.
    with_aux_value: bool,
    /// `MyNbVar`.
    nb_var: i32,
    /// `MyNbValues`.
    nb_values: i32,
    /// `MyLinearForm` (`NCollection_Array1<gp_XY>` of bounds `(0, 1)`).
    linear_form: [GpXY; 2],
    /// `MyQuadForm` (`NCollection_Array1<gp_XY>` of bounds `(0, 1)`).
    quad_form: [GpXY; 2],
    /// `MyGradient` (`math_Vector MyGradient` of bounds `(0, MyNbValues)`).
    gradient: MathVector,
    /// `MyHessian` (bounds `(0, MyNbValues + MyNbValues*(MyNbValues+1)/2)`).
    hessian: MathVector,
}

impl Energy {
    /// `FairCurve_Energy(Poles, ContrOrder1, ContrOrder2, WithAuxValue = false,
    /// Angle1 = 0, Angle2 = 0, Degree = 2, Curvature1 = 0, Curvature2 = 0)`
    /// (`.cxx:28-58`).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        poles: SharedPoles,
        contr_order1: i32,
        contr_order2: i32,
        with_aux_value: bool,
        angle1: f64,
        angle2: f64,
        degree: i32,
        curvature1: f64,
        curvature2: f64,
    ) -> Self {
        let n_poles = poles.borrow().len() as i32;
        let aux = if with_aux_value { 1 } else { 0 };
        let nb_var = 2 * (n_poles - 2) - contr_order2 - contr_order1 + aux;
        let nb_values = 2 * n_poles + aux;

        // `.cxx:50-57`: `L0 = (cos(Angle1), sin(Angle1))`,
        // `L1 = (-cos(Angle2), sin(Angle2))`, and the quad forms scaled by
        // `Degree/(Degree-1) * Curvature`.
        let scale = degree as f64 / (degree - 1) as f64;
        let linear_form = [
            GpXY::new(angle1.cos(), angle1.sin()),
            GpXY::new(-angle2.cos(), angle2.sin()),
        ];
        let quad_form = [
            GpXY::new(-angle1.sin(), angle1.cos()).multiplied_scalar(scale * curvature1),
            GpXY::new(angle2.sin(), angle2.cos()).multiplied_scalar(scale * curvature2),
        ];

        Self {
            poles,
            contr_order1,
            contr_order2,
            with_aux_value,
            nb_var,
            nb_values,
            linear_form,
            quad_form,
            gradient: MathVector::new(0, nb_values as usize),
            hessian: MathVector::new(0, (nb_values + nb_values * (nb_values + 1) / 2) as usize),
        }
    }

    /// `NbVariables()` (`.lxx:17-20`).
    pub fn nb_variables(&self) -> i32 {
        self.nb_var
    }

    /// `Poles()` (`.lxx:22-25`).
    pub fn poles(&self) -> SharedPoles {
        Rc::clone(&self.poles)
    }

    /// `MyWithAuxValue` (`.hxx:119`).
    pub fn with_aux_value(&self) -> bool {
        self.with_aux_value
    }

    /// `Indice(i, j)` (`.lxx:27-30`): `MyNbValues + j + i*(i-1)/2`.
    pub fn indice(&self, i: i32, j: i32) -> i32 {
        self.nb_values + j + i * (i - 1) / 2
    }

    /// `Value(X, E)` (`.cxx:60-68`).
    pub fn value(&mut self, x: &MathVector, e: &mut f64, compute: EnergyCompute<'_>) -> bool {
        self.compute_poles(x);
        self.value_computed(e, compute)
    }

    /// The part of `Value(X, E)` after `ComputePoles(X)`, for sub-classes whose
    /// `ComputePoles` override (`FairCurve_EnergyOfBatten.cxx:66-73`) must be
    /// the one that runs.
    pub fn value_computed(&mut self, e: &mut f64, compute: EnergyCompute<'_>) -> bool {
        let mut energie = MathVector::new(0, 0);
        let is_done = compute(0, &mut energie);
        *e = energie.value(0);
        is_done
    }

    /// `Gradient(X, G)` (`.cxx:71-79`).
    pub fn gradient(
        &mut self,
        x: &MathVector,
        g: &mut MathVector,
        compute: EnergyCompute<'_>,
    ) -> bool {
        let mut e = 0.0;
        self.values(x, &mut e, g, compute)
    }

    /// `Values(X, E, G)` (`.cxx:135-149`).
    pub fn values(
        &mut self,
        x: &MathVector,
        e: &mut f64,
        g: &mut MathVector,
        compute: EnergyCompute<'_>,
    ) -> bool {
        self.compute_poles(x);
        self.values_computed(e, g, compute)
    }

    /// The part of `Values(X, E, G)` after `ComputePoles(X)`.
    pub fn values_computed(
        &mut self,
        e: &mut f64,
        g: &mut MathVector,
        compute: EnergyCompute<'_>,
    ) -> bool {
        let is_done = compute(1, &mut self.gradient);
        if is_done {
            *e = self.gradient.value(0);
            let vect = self.gradient.clone();
            self.gradient1(&vect, g);
        }
        is_done
    }

    /// `Values(X, E, G, H)` (`.cxx:152-165`).
    pub fn values_hessian(
        &mut self,
        x: &MathVector,
        e: &mut f64,
        g: &mut MathVector,
        h: &mut MathMatrix,
        compute: EnergyCompute<'_>,
    ) -> bool {
        self.compute_poles(x);
        self.values_hessian_computed(e, g, h, compute)
    }

    /// The part of `Values(X, E, G, H)` after `ComputePoles(X)`.
    pub fn values_hessian_computed(
        &mut self,
        e: &mut f64,
        g: &mut MathVector,
        h: &mut MathMatrix,
        compute: EnergyCompute<'_>,
    ) -> bool {
        let is_done = compute(2, &mut self.hessian);
        if is_done {
            *e = self.hessian.value(0);
            let vect = self.hessian.clone();
            self.gradient1(&vect, g);
            self.hessian1(&vect, h);
        }
        is_done
    }

    /// `Variable(X)` (`.cxx:403-448`).
    pub fn variable(&self, x: &mut MathVector) -> bool {
        let n_poles = self.poles.borrow().len() as i32;
        let aux = if self.with_aux_value { 1 } else { 0 };
        let mut index_deb1 = 2; // MyPoles->Lower() + 1
        let mut index_deb2 = x.lower() as i32;
        let mut index_fin1 = n_poles - 1; // MyPoles->Upper() - 1
        let mut index_fin2 = x.upper() as i32 - aux;

        {
            let poles = self.poles.borrow();

            // calculation of variables of constraints
            if self.contr_order1 >= 1 {
                x.set_value(
                    index_deb2 as usize,
                    pole(&poles, index_deb1).distance(&pole(&poles, index_deb1 + 1)),
                );
                index_deb1 += 1;
                index_deb2 += 1;
            }
            if self.contr_order1 == 2 {
                let b1b2 = pole(&poles, index_deb1 + 1)
                    .coord
                    .subtracted(&pole(&poles, index_deb1).coord);
                x.set_value(index_deb2 as usize, b1b2.dot(&self.linear_form[0]));
                index_deb1 += 1;
                index_deb2 += 1;
            }
            if self.contr_order2 == 2 {
                let bn2bn1 = pole(&poles, index_fin1)
                    .coord
                    .subtracted(&pole(&poles, index_fin1 - 1).coord);
                x.set_value(index_fin2 as usize, -bn2bn1.dot(&self.linear_form[1]));
                index_fin1 -= 1;
                index_fin2 -= 1;
            }
            if self.contr_order2 >= 1 {
                x.set_value(
                    index_fin2 as usize,
                    pole(&poles, n_poles).distance(&pole(&poles, n_poles - 1)),
                );
                index_fin1 -= 1;
            }

            // copy poles to variables
            let mut ii = index_deb1;
            while ii <= index_fin1 {
                x.set_value(index_deb2 as usize, pole(&poles, ii).x());
                x.set_value((index_deb2 + 1) as usize, pole(&poles, ii).y());
                index_deb2 += 2;
                ii += 1;
            }
        }
        true
    }

    /// `ComputePoles(X)` (`.cxx:451-509`).
    pub fn compute_poles(&mut self, x: &MathVector) {
        let n_poles = self.poles.borrow().len() as i32;
        let aux = if self.with_aux_value { 1 } else { 0 };
        let mut index_deb1 = 2;
        let mut index_deb2 = x.lower() as i32;
        let mut index_fin1 = n_poles - 1;
        let mut index_fin2 = x.upper() as i32 - aux;

        // calculation of pole constraints
        if self.contr_order1 >= 1 {
            index_deb1 += 1;
            index_deb2 += 1;
            let p1 = pole(&self.poles.borrow(), 1);
            let mut p2 = pole(&self.poles.borrow(), 2);
            self.compute_poles_g1(0, x.value(1), &p1, &mut p2);
            self.poles.borrow_mut()[1] = p2;
        }
        if self.contr_order1 == 2 {
            index_deb1 += 1;
            index_deb2 += 1;
            let p1 = pole(&self.poles.borrow(), 1);
            let mut p2 = pole(&self.poles.borrow(), 3);
            self.compute_poles_g2(0, x.value(1), x.value(2), &p1, &mut p2);
            self.poles.borrow_mut()[2] = p2;
        }
        if self.contr_order2 == 2 {
            index_fin1 -= 1;
            index_fin2 -= 1;
            let p1 = pole(&self.poles.borrow(), n_poles);
            let mut p2 = pole(&self.poles.borrow(), n_poles - 2);
            self.compute_poles_g2(
                1,
                x.value(index_fin2 as usize),
                x.value((index_fin2 + 1) as usize),
                &p1,
                &mut p2,
            );
            self.poles.borrow_mut()[(n_poles - 3) as usize] = p2;
        }
        if self.contr_order2 >= 1 {
            index_fin1 -= 1;
            let p1 = pole(&self.poles.borrow(), n_poles);
            let mut p2 = pole(&self.poles.borrow(), n_poles - 1);
            self.compute_poles_g1(1, x.value(index_fin2 as usize), &p1, &mut p2);
            self.poles.borrow_mut()[(n_poles - 2) as usize] = p2;
        }

        // recopy others
        {
            let mut poles = self.poles.borrow_mut();
            let mut ii = index_deb1;
            let mut j = index_deb2;
            while ii <= index_fin1 {
                poles[(ii - 1) as usize] = GpPnt2d::new(x.value(j as usize), x.value((j + 1) as usize));
                j += 2;
                ii += 1;
            }
        }
    }

    /// `ComputePolesG1(Side, Lambda, P1, P2)` (`.cxx:512-519`).
    pub fn compute_poles_g1(&self, side: usize, lambda: f64, p1: &GpPnt2d, p2: &mut GpPnt2d) {
        p2.coord = p1
            .coord
            .added(&self.linear_form[side].multiplied_scalar(lambda));
    }

    /// `ComputePolesG2(Side, Lambda, Rho, P1, P2)` (`.cxx:522-530`).
    pub fn compute_poles_g2(
        &self,
        side: usize,
        lambda: f64,
        rho: f64,
        p1: &GpPnt2d,
        p2: &mut GpPnt2d,
    ) {
        p2.coord = p1
            .coord
            .added(&self.linear_form[side].multiplied_scalar(lambda + rho))
            .added(&self.quad_form[side].multiplied_scalar(lambda * lambda));
    }

    /// `Gradient1(TheVector, G)` (`.cxx:83-132`).
    pub fn gradient1(&self, vect: &MathVector, grad: &mut MathVector) {
        let n_poles = self.poles.borrow().len() as i32;
        let mut deb_g = grad.lower() as i32;
        let mut fin_g = grad.upper() as i32;
        let mut vdeb = 3i32;
        let vfin = 2 * n_poles - 2;

        // .... by calculation
        if self.contr_order1 >= 1 {
            let d_pole = GpXY::new(
                vect.value(vdeb as usize),
                vect.value((vdeb + 1) as usize),
            );
            grad.set_value(deb_g as usize, self.linear_form[0].dot(&d_pole));
            vdeb += 2;
            deb_g += 1;
        }
        if self.contr_order1 == 2 {
            let lambda0 = {
                let poles = self.poles.borrow();
                pole(&poles, 1).distance(&pole(&poles, 2))
            };
            let d_pole = GpXY::new(
                vect.value(vdeb as usize),
                vect.value((vdeb + 1) as usize),
            );
            let laux = self.linear_form[0]
                .added(&self.quad_form[0].multiplied_scalar(2.0 * lambda0));
            let previous = grad.value((deb_g - 1) as usize);
            grad.set_value((deb_g - 1) as usize, previous + laux.dot(&d_pole));
            grad.set_value(deb_g as usize, self.linear_form[0].dot(&d_pole));
            vdeb += 2;
            deb_g += 1;
        }
        if self.with_aux_value {
            grad.set_value(fin_g as usize, vect.value((2 * n_poles + 1) as usize));
            fin_g -= 1;
        }
        if self.contr_order2 >= 1 {
            let d_pole = GpXY::new(
                vect.value((vfin - 1) as usize),
                vect.value(vfin as usize),
            );
            grad.set_value(fin_g as usize, self.linear_form[1].dot(&d_pole));
            fin_g -= 1;
        }
        if self.contr_order2 == 2 {
            let lambda1 = {
                let poles = self.poles.borrow();
                pole(&poles, n_poles).distance(&pole(&poles, n_poles - 1))
            };
            let d_pole = GpXY::new(
                vect.value((vfin - 3) as usize),
                vect.value((vfin - 2) as usize),
            );
            let laux = self.linear_form[1]
                .added(&self.quad_form[1].multiplied_scalar(2.0 * lambda1));
            grad.set_value(
                fin_g as usize,
                grad.value((fin_g + 1) as usize) + laux.dot(&d_pole),
            );
            grad.set_value((fin_g + 1) as usize, self.linear_form[1].dot(&d_pole));
            fin_g -= 1;
        }
        // ... by recopy
        let mut ii = deb_g;
        while ii <= fin_g {
            grad.set_value(ii as usize, vect.value(vdeb as usize));
            vdeb += 1;
            ii += 1;
        }
    }

    /// `Hessian1(TheVector, H)` (`.cxx:168-401`).
    pub fn hessian1(&self, vect: &MathVector, h: &mut MathMatrix) {
        let n_poles = self.poles.borrow().len() as i32;
        let aux = if self.with_aux_value { 1 } else { 0 };

        let vdeb = 3 + 2 * self.contr_order1;
        let vfin = 2 * n_poles - 2 * (self.contr_order2 + 1);
        let vup = 2 * n_poles + aux;
        let deb_h = 1 + self.contr_order1;
        let fin_h = self.nb_var - aux - self.contr_order2;

        let l0 = self.linear_form[0];
        let l1 = self.linear_form[1];
        let q0 = self.quad_form[0];
        let q1 = self.quad_form[1];
        let cos0 = l0.x * l0.x;
        let sin0 = l0.y * l0.y;
        let cos_sin0 = 2.0 * l0.x * l0.y;
        let cos1 = l1.x * l1.x;
        let sin1 = l1.y * l1.y;
        let cos_sin1 = 2.0 * l1.x * l1.y;
        let mut lambda0 = 0.0;
        let mut lambda1 = 0.0;

        if self.contr_order1 >= 1 {
            let poles = self.poles.borrow();
            lambda0 = pole(&poles, 1).distance(&pole(&poles, 2));
        }
        if self.contr_order2 >= 1 {
            let poles = self.poles.borrow();
            lambda1 = pole(&poles, n_poles).distance(&pole(&poles, n_poles - 1));
        }

        let v = |i: i32| vect.value(i as usize);

        if self.contr_order1 >= 1 {
            // calculate the left lambda column --------------------------------
            let mut jj = vdeb - 2 * self.contr_order1;
            let mut kk = self.indice(jj, jj); // X2X2
            let mut ii = self.indice(jj + 1, jj); // X2Y2
            h.set_value(1, 1, cos0 * v(kk) + cos_sin0 * v(ii) + sin0 * v(ii + 1));

            if self.contr_order1 >= 2 {
                let laux = l0.added(&q0.multiplied_scalar(2.0 * lambda0));
                jj = vdeb - 2 * (self.contr_order1 - 1);
                kk = self.indice(jj, jj - 2); // X1X2
                ii = self.indice(jj + 1, jj - 2); // X1Y2
                let aux_xy = GpXY::new(v(kk + 2), v(ii + 3));

                h.set_value(
                    1,
                    1,
                    h.value(1, 1)
                        + 2.0
                            * ((q0.x * v(5) + q0.y * v(6))
                                + (laux.x * (l0.x * v(kk) + l0.y * v(kk + 1))
                                    + laux.y * (l0.x * v(ii) + l0.y * v(ii + 1)))
                                + laux.x * laux.y * v(ii + 2))
                        + (laux.x * laux.x * v(kk + 2) + laux.y * laux.y * v(ii + 3)),
                );

                h.set_value(
                    2,
                    1,
                    (cos0 * v(kk) + cos_sin0 * (v(ii) + v(kk + 1)) / 2.0 + sin0 * v(ii + 1))
                        + laux.dot(&l0.multiplied_xy(&aux_xy))
                        + (laux.x * l0.y + laux.y * l0.x) * v(ii + 2),
                );
            }

            if self.with_aux_value {
                kk = self.indice(vup, vdeb - 2 * self.contr_order1);
                h.set_value(self.nb_var as usize, 1, l0.x * v(kk) + l0.y * v(kk + 1));
            }

            if self.contr_order2 >= 1 {
                h.set_value((self.nb_var - aux) as usize, 1, 0.0); // correct if there are less than 3 nodes
                if self.contr_order2 == 2 {
                    h.set_value((self.nb_var - aux - 1) as usize, 1, 0.0);
                }
            }

            let mut vk = vdeb;
            kk = self.indice(vk, vdeb - 2 * self.contr_order1);
            for ii in deb_h..=fin_h {
                h.set_value(ii as usize, 1, l0.x * v(kk) + l0.y * v(kk + 1));
                kk += vk;
                vk += 1;
            }
        }

        // calculate the left line mu ----------------------
        if self.contr_order1 >= 2 {
            let jj = vdeb - 2 * (self.contr_order1 - 1);
            let mut kk = self.indice(jj, jj); // X3X3
            let ii = self.indice(jj + 1, jj); // X3Y3
            h.set_value(2, 2, cos0 * v(kk) + cos_sin0 * v(ii) + sin0 * v(ii + 1));

            if self.with_aux_value {
                kk = self.indice(vup, vdeb - 2 * (self.contr_order1 - 1));
                let pole_xy = GpXY::new(v(kk), v(kk + 1));
                let previous = h.value(self.nb_var as usize, 1);
                h.set_value(
                    self.nb_var as usize,
                    1,
                    previous
                        + l0.added(&q0.multiplied_scalar(2.0 * lambda0)).dot(&pole_xy),
                );
                h.set_value(self.nb_var as usize, 2, l0.x * v(kk) + l0.y * v(kk + 1));
            }

            if self.contr_order2 >= 1 {
                h.set_value((self.nb_var - aux) as usize, 2, 0.0); // correct if there are less than 3 nodes
                if self.contr_order2 == 2 {
                    h.set_value((self.nb_var - aux - 1) as usize, 2, 0.0);
                }
            }
            let mut vk = vdeb;

            let laux = l0.added(&q0.multiplied_scalar(2.0 * lambda0));
            let xaux = laux.x;
            let yaux = laux.y;

            kk = self.indice(vk, vdeb - 2 * self.contr_order1 + 2);
            for ii in deb_h..=fin_h {
                h.set_value(ii as usize, 2, l0.x * v(kk) + l0.y * v(kk + 1));
                let previous = h.value(ii as usize, 1);
                h.set_value(ii as usize, 1, previous + xaux * v(kk) + yaux * v(kk + 1));
                kk += vk;
                vk += 1;
            }
        }

        // calculate the right lambda line -----------------------
        if self.contr_order2 >= 1 {
            let mut jj = fin_h + 1;
            let vk = vfin + 2 * self.contr_order2 - 1;
            let mut kk = self.indice(vk, vdeb);
            for ii in deb_h..=fin_h {
                h.set_value(jj as usize, ii as usize, l1.x * v(kk) + l1.y * v(kk + vk));
                kk += 1;
            }

            kk = self.indice(vk, vk);
            h.set_value(
                jj as usize,
                jj as usize,
                cos1 * v(kk) + cos_sin1 * v(kk + vk) + sin1 * v(kk + vk + 1),
            );

            if self.contr_order2 >= 2 {
                let laux = l1.added(&q1.multiplied_scalar(2.0 * lambda1));
                jj = vfin + 2 * self.contr_order2 - 3;
                let kk2 = self.indice(jj + 2, jj); // Xn-1Xn-2
                let ii2 = self.indice(jj + 3, jj); // Yn-1Xn-2
                let ll = self.indice(jj, jj);

                let previous = h.value((fin_h + 1) as usize, (fin_h + 1) as usize);
                let add = 2.0
                    * ((q1.x * v(jj) + q1.y * v(jj + 1))
                        + (laux.x * (l1.x * v(kk2) + l1.y * v(ii2))
                            + laux.y * (l1.x * v(kk2 + 1) + l1.y * v(ii2 + 1)))
                        + laux.x * laux.y * v(ll + jj))
                    + (laux.x * laux.x * v(ll) + laux.y * laux.y * v(ll + jj + 1));
                h.set_value((fin_h + 1) as usize, (fin_h + 1) as usize, previous + add);

                let mut value =
                    cos1 * v(kk2) + cos_sin1 * (v(ii2) + v(kk2 + 1)) / 2.0 + sin1 * v(ii2 + 1);
                let aux_xy = GpXY::new(v(ll), v(ll + jj + 1));
                value += laux.dot(&l1.multiplied_xy(&aux_xy))
                    + (laux.x * l1.y + laux.y * l1.x) * v(ll + jj);
                h.set_value((fin_h + 2) as usize, (fin_h + 1) as usize, value);
            }
        }

        // calculate the right line mu  -----------------------
        if self.contr_order2 >= 2 {
            let jj = fin_h + 2;
            let vk = vfin + 2 * self.contr_order2 - 3;
            let mut kk = self.indice(vk, vdeb);

            let laux = l1.added(&q1.multiplied_scalar(2.0 * lambda1));
            let xaux = laux.x;
            let yaux = laux.y;

            for ii in deb_h..=fin_h {
                h.set_value(jj as usize, ii as usize, l1.x * v(kk) + l1.y * v(kk + vk));
                // update the right line Lambda
                let previous = h.value((jj - 1) as usize, ii as usize);
                h.set_value(
                    (jj - 1) as usize,
                    ii as usize,
                    previous + xaux * v(kk) + yaux * v(kk + vk),
                );
                kk += 1;
            }
            kk = self.indice(vk, vk);
            let ii = self.indice(vk + 1, vk);
            h.set_value(
                jj as usize,
                jj as usize,
                cos1 * v(kk) + cos_sin1 * v(ii) + sin1 * v(ii + 1),
            );
        }

        // calculate the Auxiliary Variable line -----------------------
        if self.with_aux_value {
            let mut kk = self.indice(vup, vdeb);
            for ii in deb_h..=fin_h {
                h.set_value(self.nb_var as usize, ii as usize, v(kk));
                kk += 1;
            }

            if self.contr_order2 >= 1 {
                kk = self.indice(vup, vfin + 2 * self.contr_order2 - 1);
                h.set_value(
                    self.nb_var as usize,
                    (fin_h + 1) as usize,
                    l1.x * v(kk) + l1.y * v(kk + 1),
                );
            }
            if self.contr_order2 >= 2 {
                kk = self.indice(vup, vfin + 2 * self.contr_order2 - 3);
                let pole_xy = GpXY::new(v(kk), v(kk + 1));
                let previous = h.value(self.nb_var as usize, (fin_h + 1) as usize);
                h.set_value(
                    self.nb_var as usize,
                    (fin_h + 1) as usize,
                    previous + l1.added(&q1.multiplied_scalar(2.0 * lambda1)).dot(&pole_xy),
                );
                h.set_value(
                    self.nb_var as usize,
                    (fin_h + 2) as usize,
                    l1.dot(&pole_xy),
                );
            }
            kk = self.indice(vup, vup);
            let upper_row = h.upper_row();
            h.set_value(upper_row, upper_row, v(kk));
        }

        // recopy the internal block -----------------------------------
        let mut kk = self.indice(vdeb, vdeb);
        for ii in deb_h..=fin_h {
            let mut jj = deb_h;
            while jj <= ii {
                h.set_value(ii as usize, jj as usize, v(kk));
                kk += 1;
                jj += 1;
            }
            kk += vdeb - 1;
        }
        // symmetry
        let lower_row = h.lower_row();
        let upper_row = h.upper_row();
        for ii in lower_row..=upper_row {
            for jj in (ii + 1)..=upper_row {
                let value = h.value(jj, ii);
                h.set_value(ii, jj, value);
            }
        }
    }
}
