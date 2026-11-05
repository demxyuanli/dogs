//! `FairCurve_DistributionOfSagging`
//! (`FairCurve_DistributionOfSagging.hxx:33-56`, `.cxx:28-237`): the local
//! flexion term `(Height^3/12) * W^2` with `W = (C' x C'') / |C'|^2.5`, plus
//! its gradient (`.cxx:92-141`) and Hessian (`.cxx:143-233`).

use occt_core::bspl::eval_basis::eval_bspline_basis;
use occt_core::gp::xy::GpXY;
use occt_math::{MathFunctionSet, MathVector};

use super::batten_law::BattenLaw;
use super::distribution_of_energy::{pole, DistributionOfEnergy, SharedPoles};

/// `FairCurve_DistributionOfSagging`
/// (`FairCurve_DistributionOfSagging.hxx:33-56`).
#[derive(Debug, Clone)]
pub struct DistributionOfSagging {
    /// `FairCurve_DistributionOfEnergy` base state.
    pub(crate) base: DistributionOfEnergy,
    /// `MyLaw`.
    my_law: BattenLaw,
}

impl DistributionOfSagging {
    /// `FairCurve_DistributionOfSagging(BSplOrder, FlatKnots, Poles,
    /// DerivativeOrder, Law, NbValAux = 0)` (`.cxx:28-38`).
    pub fn new(
        bspl_order: i32,
        flat_knots: Vec<f64>,
        poles: SharedPoles,
        derivative_order: i32,
        law: &BattenLaw,
        nb_val_aux: i32,
    ) -> Self {
        Self {
            base: DistributionOfEnergy::new(
                bspl_order,
                flat_knots,
                poles,
                derivative_order,
                nb_val_aux,
            ),
            my_law: law.clone(),
        }
    }

    /// `SetDerivativeOrder(DerivativeOrder)` (`.hxx:43` of the base).
    pub fn set_derivative_order(&mut self, derivative_order: i32) {
        self.base.set_derivative_order(derivative_order);
    }
}

impl MathFunctionSet for DistributionOfSagging {
    fn nb_variables(&self) -> usize {
        self.base.nb_variables()
    }

    fn nb_equations(&self) -> usize {
        self.base.nb_equations()
    }

    /// `Value(TParam, Flexion)` (`.cxx:40-237`).
    fn value(&mut self, t_param: &MathVector, flexion: &mut MathVector) -> bool {
        let ok = true;
        let bspl_order = self.base.bspl_order;
        let nb_val_aux = self.base.nb_val_aux;
        let derivative_order = self.base.derivative_order;

        // `(0.0) initialisations generales` (`.cxx:47-48`).
        flexion.init(0.0);

        // `BSplCLib::EvalBsplineBasis(2, MyBSplOrder, MyFlatKnots->Array1(),
        // TParam(TParam.Lower()), FirstNonZero, Base)` (`.cxx:56-67`).
        let (first_non_zero0, base) = match eval_bspline_basis(
            2,
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

        // `(0.1) evaluation de CPrim et CScn` (`.cxx:69-73`).
        let mut c_prim = GpXY::zero();
        let mut c_secn = GpXY::zero();
        {
            let poles = self.base.poles.borrow();
            for ii in 1..=bspl_order {
                let coord = pole(&poles, ii + last_zero).coord;
                c_prim.add(&coord.multiplied_scalar(base[1][(ii - 1) as usize]));
                c_secn.add(&coord.multiplied_scalar(base[2][(ii - 1) as usize]));
            }
        }

        // `(1) Evaluation de la flexion locale = W*W` (`.cxx:75-92`).
        let norme_c_prim = c_prim.modulus();
        let inv_norme_c_prim = 1.0 / norme_c_prim;
        let numerateur = c_prim.crossed(&c_secn);
        let denominateur = norme_c_prim.powf(2.5);
        let hauteur = match self.my_law.value(t_param.value(t_param.lower())) {
            Some(value) => value,
            // `if (!Ok) { return Ok; }` (`.cxx:85-88`).
            None => return false,
        };
        let mesure = hauteur * hauteur * hauteur / 12.0;
        let w_val = numerateur / denominateur;
        flexion.set_value(flexion.lower(), mesure * w_val * w_val);

        if derivative_order < 1 {
            return ok;
        }

        // `(2) Evaluation du gradient de la flexion locale` (`.cxx:94-141`).
        let n_poles = self.base.poles.borrow().len() as i32;
        let n = (2 * bspl_order + nb_val_aux) as usize;
        let mut w_grad = MathVector::new(1, n);
        let mut num_grad = MathVector::new(1, n);
        let mut grad_norme_c_prim = MathVector::new(1, n);
        let mut numdu_grad = MathVector::new(1, n);

        let x_prim = c_prim.x;
        let y_prim = c_prim.y;
        let x_secn = c_secn.x;
        let y_secn = c_secn.y;
        let inv_denominateur = 1.0 / denominateur;
        let mut facteur = 2.0 * mesure * w_val;
        let aux0 = 2.5 * numerateur * inv_norme_c_prim;
        let mut kk = flexion.lower() as i32 + first_non_zero;

        let mut jj = 1i32;
        for ii in 1..=bspl_order {
            let b1 = base[1][(ii - 1) as usize];
            let b2 = base[2][(ii - 1) as usize];

            // `(2.1) Derivation en X` (`.cxx:115-121`).
            num_grad.set_value(jj as usize, y_secn * b1 - y_prim * b2);
            grad_norme_c_prim.set_value(jj as usize, x_prim * b1 * inv_norme_c_prim);
            numdu_grad.set_value(
                jj as usize,
                num_grad.value(jj as usize) - aux0 * grad_norme_c_prim.value(jj as usize),
            );
            w_grad.set_value(jj as usize, inv_denominateur * numdu_grad.value(jj as usize));
            flexion.set_value(kk as usize, facteur * w_grad.value(jj as usize));
            jj += 1;

            // `(2.2) Derivation en Y` (`.cxx:123-129`).
            num_grad.set_value(jj as usize, -x_secn * b1 + x_prim * b2);
            grad_norme_c_prim.set_value(jj as usize, y_prim * b1 * inv_norme_c_prim);
            numdu_grad.set_value(
                jj as usize,
                num_grad.value(jj as usize) - aux0 * grad_norme_c_prim.value(jj as usize),
            );
            w_grad.set_value(jj as usize, inv_denominateur * numdu_grad.value(jj as usize));
            flexion.set_value((kk + 1) as usize, facteur * w_grad.value(jj as usize));
            jj += 1;

            kk += 2;
        }

        let last_gradient_index;
        if nb_val_aux == 1 {
            // `(2.3) Gestion de la variable de glissement` (`.cxx:132-138`).
            last_gradient_index = flexion.lower() as i32 + 2 * n_poles + 1;
            let upper = w_grad.upper();
            w_grad.set_value(upper, 0.0);
        } else {
            last_gradient_index = flexion.lower() as i32 + 2 * n_poles;
        }

        if derivative_order < 2 {
            return ok;
        }

        // `(3) Evaluation du Hessien de la tension locale` (`.cxx:146-233`).
        let facteur_x = (1.0 - (x_prim * inv_norme_c_prim) * (x_prim * inv_norme_c_prim))
            * inv_norme_c_prim;
        let facteur_y = (1.0 - (y_prim * inv_norme_c_prim) * (y_prim * inv_norme_c_prim))
            * inv_norme_c_prim;
        let facteur_xy = -(x_prim * inv_norme_c_prim) * (y_prim * inv_norme_c_prim)
            * inv_norme_c_prim;
        let facteur_w = w_val * inv_norme_c_prim;

        let mut aux;
        let mut d_seconde;
        let mut v_intermed;

        facteur = 2.0 * mesure;

        let mut kk = first_non_zero;
        let mut k2 = last_gradient_index + (kk - 1) * kk / 2;

        let mut ii = 2i32;
        while ii <= 2 * bspl_order {
            let ii2 = ii / 2;
            let mut k1 = k2 + first_non_zero;
            k2 = k1 + kk;
            kk += 2;

            let b1_ii = base[1][(ii2 - 1) as usize];
            let b2_ii = base[2][(ii2 - 1) as usize];

            let mut jj = 2i32;
            while jj < ii {
                let jj2 = jj / 2;
                let b1_jj = base[1][(jj2 - 1) as usize];
                let b2_jj = base[2][(jj2 - 1) as usize];
                let produit = b1_ii * b1_jj;
                let n_seconde = b1_ii * b2_jj - b2_ii * b1_jj;

                // derivation en XiXj (`.cxx:174-181`).
                d_seconde = facteur_x * produit;
                aux = num_grad.value((ii - 1) as usize) * grad_norme_c_prim.value((jj - 1) as usize)
                    - 2.5
                        * (num_grad.value((jj - 1) as usize)
                            * grad_norme_c_prim.value((ii - 1) as usize)
                            + d_seconde * numerateur);
                v_intermed = inv_denominateur
                    * (aux
                        - 3.5 * grad_norme_c_prim.value((jj - 1) as usize)
                            * numdu_grad.value((ii - 1) as usize));
                flexion.set_value(
                    k1 as usize,
                    facteur
                        * (w_grad.value((ii - 1) as usize) * w_grad.value((jj - 1) as usize)
                            + facteur_w * v_intermed),
                );
                k1 += 1;

                // derivation en XiYj (`.cxx:183-190`).
                d_seconde = facteur_xy * produit;
                aux = norme_c_prim * n_seconde
                    + num_grad.value((ii - 1) as usize) * grad_norme_c_prim.value(jj as usize)
                    - 2.5
                        * (num_grad.value(jj as usize) * grad_norme_c_prim.value((ii - 1) as usize)
                            + d_seconde * numerateur);
                v_intermed = inv_denominateur
                    * (aux
                        - 3.5 * grad_norme_c_prim.value(jj as usize)
                            * numdu_grad.value((ii - 1) as usize));
                flexion.set_value(
                    k1 as usize,
                    facteur
                        * (w_grad.value((ii - 1) as usize) * w_grad.value(jj as usize)
                            + facteur_w * v_intermed),
                );
                k1 += 1;

                // derivation en YiXj (`.cxx:192-199`).
                aux = -norme_c_prim * n_seconde
                    + num_grad.value(ii as usize) * grad_norme_c_prim.value((jj - 1) as usize)
                    - 2.5
                        * (num_grad.value((jj - 1) as usize) * grad_norme_c_prim.value(ii as usize)
                            + d_seconde * numerateur);
                v_intermed = inv_denominateur
                    * (aux
                        - 3.5 * grad_norme_c_prim.value((jj - 1) as usize)
                            * numdu_grad.value(ii as usize));
                flexion.set_value(
                    k2 as usize,
                    facteur
                        * (w_grad.value(ii as usize) * w_grad.value((jj - 1) as usize)
                            + facteur_w * v_intermed),
                );
                k2 += 1;

                // derivation en YiYj (`.cxx:201-208`).
                d_seconde = facteur_y * produit;
                aux = num_grad.value(ii as usize) * grad_norme_c_prim.value(jj as usize)
                    - 2.5
                        * (num_grad.value(jj as usize) * grad_norme_c_prim.value(ii as usize)
                            + d_seconde * numerateur);
                v_intermed = inv_denominateur
                    * (aux
                        - 3.5 * grad_norme_c_prim.value(jj as usize)
                            * numdu_grad.value(ii as usize));
                flexion.set_value(
                    k2 as usize,
                    facteur
                        * (w_grad.value(ii as usize) * w_grad.value(jj as usize)
                            + facteur_w * v_intermed),
                );
                k2 += 1;

                jj += 2;
            }

            // case where jj = ii: triangular fill (`.cxx:210-232`).
            let produit = b1_ii * b1_ii;

            // derivation en XiXi.
            d_seconde = facteur_x * produit;
            aux = -1.5 * num_grad.value((ii - 1) as usize) * grad_norme_c_prim.value((ii - 1) as usize)
                - 2.5 * d_seconde * numerateur;
            v_intermed = inv_denominateur
                * (aux
                    - 3.5 * grad_norme_c_prim.value((ii - 1) as usize)
                        * numdu_grad.value((ii - 1) as usize));
            flexion.set_value(
                k1 as usize,
                facteur
                    * (w_grad.value((ii - 1) as usize) * w_grad.value((ii - 1) as usize)
                        + facteur_w * v_intermed),
            );

            // derivation en XiYi.
            d_seconde = facteur_xy * produit;
            aux = num_grad.value((ii - 1) as usize) * grad_norme_c_prim.value(ii as usize)
                - 2.5
                    * (num_grad.value(ii as usize) * grad_norme_c_prim.value((ii - 1) as usize)
                        + d_seconde * numerateur);
            v_intermed = inv_denominateur
                * (aux
                    - 3.5 * grad_norme_c_prim.value(ii as usize)
                        * numdu_grad.value((ii - 1) as usize));
            flexion.set_value(
                k2 as usize,
                facteur
                    * (w_grad.value(ii as usize) * w_grad.value((ii - 1) as usize)
                        + facteur_w * v_intermed),
            );
            k2 += 1;

            // derivation en YiYi.
            d_seconde = facteur_y * produit;
            aux = -1.5 * num_grad.value(ii as usize) * grad_norme_c_prim.value(ii as usize)
                - 2.5 * d_seconde * numerateur;
            v_intermed = inv_denominateur
                * (aux
                    - 3.5 * grad_norme_c_prim.value(ii as usize) * numdu_grad.value(ii as usize));
            flexion.set_value(
                k2 as usize,
                facteur
                    * (w_grad.value(ii as usize) * w_grad.value(ii as usize)
                        + facteur_w * v_intermed),
            );

            ii += 2;
        }

        ok
    }
}
