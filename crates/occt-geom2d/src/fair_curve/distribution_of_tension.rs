//! `FairCurve_DistributionOfTension`
//! (`FairCurve_DistributionOfTension.hxx:33-61`, `.cxx:28-217`): the local
//! tension term `Height * (|C'| - Sliding)^2 / Sliding` together with its
//! gradient (`.cxx:105-140`) and Hessian (`.cxx:141-214`) packed in OCCT's
//! triangular layout.

use occt_core::bspl::eval_basis::eval_bspline_basis;
use occt_core::gp::xy::GpXY;
use occt_math::{MathFunctionSet, MathVector};

use super::batten_law::BattenLaw;
use super::distribution_of_energy::{pole, DistributionOfEnergy, SharedPoles};

/// `FairCurve_DistributionOfTension`
/// (`FairCurve_DistributionOfTension.hxx:33-61`).
#[derive(Debug, Clone)]
pub struct DistributionOfTension {
    /// `FairCurve_DistributionOfEnergy` base state.
    pub(crate) base: DistributionOfEnergy,
    /// `MyLengthSliding`.
    my_length_sliding: f64,
    /// `MyLaw`.
    my_law: BattenLaw,
    /// `MyHeight` (non-zero only in the `Uniform` mode, `.cxx:39-47`).
    my_height: f64,
}

impl DistributionOfTension {
    /// `FairCurve_DistributionOfTension(BSplOrder, FlatKnots, Poles,
    /// DerivativeOrder, LengthSliding, Law, NbValAux = 0, Uniform = false)`
    /// (`.cxx:28-48`).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        bspl_order: i32,
        flat_knots: Vec<f64>,
        poles: SharedPoles,
        derivative_order: i32,
        length_sliding: f64,
        law: &BattenLaw,
        nb_val_aux: i32,
        uniform: bool,
    ) -> Self {
        let my_height = if uniform {
            // `MyLaw.Value(0.5, MyHeight)` (`.cxx:41`); the `Value` result is
            // discarded in OCCT, a rejected height leaves `MyHeight` at 0.
            law.value(0.5).unwrap_or(0.0)
        } else {
            0.0
        };
        Self {
            base: DistributionOfEnergy::new(
                bspl_order,
                flat_knots,
                poles,
                derivative_order,
                nb_val_aux,
            ),
            my_length_sliding: length_sliding,
            my_law: law.clone(),
            my_height,
        }
    }

    /// `SetLengthSliding(LengthSliding)`
    /// (`FairCurve_DistributionOfTension.lxx:18-21`).
    pub fn set_length_sliding(&mut self, length_sliding: f64) {
        self.my_length_sliding = length_sliding;
    }

    /// `SetDerivativeOrder(DerivativeOrder)` (`.hxx:43` of the base).
    pub fn set_derivative_order(&mut self, derivative_order: i32) {
        self.base.set_derivative_order(derivative_order);
    }
}

impl MathFunctionSet for DistributionOfTension {
    fn nb_variables(&self) -> usize {
        self.base.nb_variables()
    }

    fn nb_equations(&self) -> usize {
        self.base.nb_equations()
    }

    /// `Value(TParam, FTension)` (`.cxx:50-217`).
    fn value(&mut self, t_param: &MathVector, f_tension: &mut MathVector) -> bool {
        let ok = true;
        let bspl_order = self.base.bspl_order;
        let nb_val_aux = self.base.nb_val_aux;
        let derivative_order = self.base.derivative_order;

        // `(0.0) initialisations generales` (`.cxx:58-59`).
        f_tension.init(0.0);

        // `BSplCLib::EvalBsplineBasis(1, MyBSplOrder, MyFlatKnots->Array1(),
        // TParam(TParam.Lower()), FirstNonZero, Base)` (`.cxx:67-78`).
        let (first_non_zero0, base) = match eval_bspline_basis(
            1,
            bspl_order,
            &self.base.flat_knots,
            t_param.value(t_param.lower()),
            false,
        ) {
            Ok(value) => value,
            Err(_) => return false,
        };
        let last_zero = first_non_zero0 - 1;
        let first_non_zero = 2 * last_zero + 1;

        // `(0.1) evaluation de CPrim` (`.cxx:81-84`).
        let mut c_prim = GpXY::zero();
        {
            let poles = self.base.poles.borrow();
            for ii in 1..=bspl_order {
                let b = base[1][(ii - 1) as usize];
                c_prim.add(&pole(&poles, ii + last_zero).coord.multiplied_scalar(b));
            }
        }

        // `(1) Evaluation de la tension locale` (`.cxx:86-102`).
        let norme_c_prim = c_prim.modulus();
        let hauteur;
        if self.my_height > 0.0 {
            hauteur = self.my_height;
        } else {
            match self.my_law.value(t_param.value(t_param.lower())) {
                Some(value) => hauteur = value,
                // `if (!Ok) { return Ok; }` (`.cxx:97-101`).
                None => return false,
            }
        }
        let difference = norme_c_prim - self.my_length_sliding;
        f_tension.set_value(
            f_tension.lower(),
            hauteur * difference * difference / self.my_length_sliding,
        );

        if derivative_order < 1 {
            return ok;
        }

        // `(2) Evaluation du gradient de la tension locale` (`.cxx:107-140`).
        let n_poles = self.base.poles.borrow().len() as i32;
        let mut grad_difference = MathVector::new(1, (2 * bspl_order + nb_val_aux) as usize);
        let xaux = c_prim.x / norme_c_prim;
        let yaux = c_prim.y / norme_c_prim;
        let mut facteur = 2.0 * hauteur * difference / self.my_length_sliding;

        let mut kk = f_tension.lower() as i32 + first_non_zero;
        let mut jj = 1i32;
        for ii in 1..=bspl_order {
            let b = base[1][(ii - 1) as usize];
            grad_difference.set_value(jj as usize, b * xaux);
            f_tension.set_value(kk as usize, facteur * grad_difference.value(jj as usize));
            jj += 1;
            grad_difference.set_value(jj as usize, b * yaux);
            f_tension.set_value((kk + 1) as usize, facteur * grad_difference.value(jj as usize));
            jj += 1;
            kk += 2;
        }
        let last_gradient_index;
        if nb_val_aux == 1 {
            last_gradient_index = f_tension.lower() as i32 + 2 * n_poles + 1;
            let upper = grad_difference.upper();
            grad_difference.set_value(
                upper,
                1.0 - (norme_c_prim / self.my_length_sliding) * (norme_c_prim / self.my_length_sliding),
            );
            f_tension.set_value(
                last_gradient_index as usize,
                hauteur * grad_difference.value(upper),
            );
        } else {
            last_gradient_index = f_tension.lower() as i32 + 2 * n_poles;
        }

        if derivative_order < 2 {
            return ok;
        }

        // `(3) Evaluation du Hessien de la tension locale` (`.cxx:141-214`).
        let facteur_x = difference * (1.0 - xaux * xaux) / norme_c_prim;
        let facteur_y = difference * (1.0 - yaux * yaux) / norme_c_prim;
        let facteur_xy = -difference * xaux * yaux / norme_c_prim;
        facteur = 2.0 * hauteur / self.my_length_sliding;

        let mut kk = first_non_zero;
        let mut k2 = last_gradient_index + (kk - 1) * kk / 2;

        let mut ii = 2i32;
        while ii <= 2 * bspl_order {
            let mut k1 = k2 + first_non_zero;
            k2 = k1 + kk;
            kk += 2;
            let b_ii = base[1][(ii / 2 - 1) as usize];

            let mut jj = 2i32;
            while jj < ii {
                let produit = b_ii * base[1][(jj / 2 - 1) as usize];
                f_tension.set_value(
                    k1 as usize,
                    facteur
                        * (grad_difference.value((ii - 1) as usize)
                            * grad_difference.value((jj - 1) as usize)
                            + facteur_x * produit),
                );
                k1 += 1;
                f_tension.set_value(
                    k1 as usize,
                    facteur
                        * (grad_difference.value(ii as usize)
                            * grad_difference.value((jj - 1) as usize)
                            + facteur_xy * produit),
                );
                k1 += 1;
                f_tension.set_value(
                    k2 as usize,
                    facteur
                        * (grad_difference.value((ii - 1) as usize)
                            * grad_difference.value(jj as usize)
                            + facteur_xy * produit),
                );
                k2 += 1;
                f_tension.set_value(
                    k2 as usize,
                    facteur
                        * (grad_difference.value(ii as usize) * grad_difference.value(jj as usize)
                            + facteur_y * produit),
                );
                k2 += 1;
                jj += 2;
            }

            // case where jj = ii: triangular fill (`.cxx:182-193`).
            let produit = b_ii * b_ii;
            f_tension.set_value(
                k1 as usize,
                facteur
                    * (grad_difference.value((ii - 1) as usize)
                        * grad_difference.value((ii - 1) as usize)
                        + facteur_x * produit),
            );
            f_tension.set_value(
                k2 as usize,
                facteur
                    * (grad_difference.value(ii as usize) * grad_difference.value((ii - 1) as usize)
                        + facteur_xy * produit),
            );
            k2 += 1;
            f_tension.set_value(
                k2 as usize,
                facteur
                    * (grad_difference.value(ii as usize) * grad_difference.value(ii as usize)
                        + facteur_y * produit),
            );

            ii += 2;
        }

        if nb_val_aux == 1 {
            // `.cxx:195-212`.
            let facteur_x = -2.0 * c_prim.x * hauteur
                / (self.my_length_sliding * self.my_length_sliding);
            let facteur_y = -2.0 * c_prim.y * hauteur
                / (self.my_length_sliding * self.my_length_sliding);

            let ii = last_gradient_index - f_tension.lower() as i32;
            let mut kk = last_gradient_index + (ii - 1) * ii / 2 + first_non_zero;
            for ii in 1..=bspl_order {
                let b = base[1][(ii - 1) as usize];
                f_tension.set_value(kk as usize, facteur_x * b);
                kk += 1;
                f_tension.set_value(kk as usize, facteur_y * b);
                kk += 1;
            }
            f_tension.set_value(
                f_tension.upper(),
                2.0 * hauteur * (norme_c_prim / self.my_length_sliding)
                    * (norme_c_prim / self.my_length_sliding)
                    / self.my_length_sliding,
            );
        }

        ok
    }
}
