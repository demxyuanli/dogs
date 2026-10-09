//! `AdvApprox_SimpleApprox` — one-interval Jacobi approximation.
//! Source: `AdvApprox_SimpleApprox.cxx`.

use crate::bspl::plib_eval::{eval_poly0, hermite_interpolate};
use crate::bspl::plib_jacobi::JacobiPolynomial;
use crate::kernel::geomabs::Shape;

/// Evaluator: `param`, derivative request, writes `dimension` values, returns error code.
pub type EvalFn<'a> = &'a dyn Fn(f64, i32, &mut [f64]) -> i32;

pub struct SimpleApprox {
    total_num_ss: i32,
    total_dimension: usize,
    nb_gauss: i32,
    work_degree: i32,
    niv_constr: i32,
    jac: JacobiPolynomial,
    tab_points: Vec<f64>,
    tab_weights: Vec<Vec<f64>>,
    coeff: Vec<f64>,
    first_constr: Vec<Vec<f64>>,
    last_constr: Vec<Vec<f64>>,
    som: Vec<f64>,
    dif: Vec<f64>,
    degree: i32,
    max_error: Vec<f64>,
    average_error: Vec<f64>,
    done: bool,
}

impl SimpleApprox {
    pub fn new(
        total_dimension: usize,
        total_num_ss: i32,
        continuity: Shape,
        work_degree: i32,
        nb_gauss: i32,
        jac: JacobiPolynomial,
    ) -> Result<Self, ()> {
        let niv = match continuity {
            Shape::C0 => 0,
            Shape::C1 => 1,
            Shape::C2 => 2,
            _ => return Err(()),
        };
        let tab_points = jac.points(nb_gauss)?;
        let tab_weights = jac.weights(nb_gauss)?;
        let coeff_len = ((work_degree + 1) as usize) * total_dimension;
        let half = (nb_gauss / 2) as usize;
        Ok(Self {
            total_num_ss,
            total_dimension,
            nb_gauss,
            work_degree,
            niv_constr: niv,
            jac,
            tab_points,
            tab_weights,
            coeff: vec![0.0; coeff_len],
            first_constr: vec![vec![0.0; (niv + 1) as usize]; total_dimension],
            last_constr: vec![vec![0.0; (niv + 1) as usize]; total_dimension],
            som: vec![0.0; (half + 1) * total_dimension],
            dif: vec![0.0; (half + 1) * total_dimension],
            degree: 0,
            max_error: Vec::new(),
            average_error: Vec::new(),
            done: false,
        })
    }

    pub fn perform(
        &mut self,
        local_dimension: &[i32],
        local_tol: &[f64],
        first: f64,
        last: f64,
        max_degree: i32,
        eval: EvalFn<'_>,
    ) {
        self.done = false;
        let dim = self.total_dimension;
        let degree_r = 2 * self.niv_constr + 1;
        let degree_q = self.work_degree - 2 * (self.niv_constr + 1);
        let fact = (last - first) / 2.0;
        let mut result = vec![0.0; dim];
        for derive in (0..=self.niv_constr).rev() {
            if eval(first, derive, &mut result) != 0 {
                return;
            }
            if derive >= 1 {
                for v in &mut result {
                    *v *= fact;
                }
            }
            if derive == 2 {
                for v in &mut result {
                    *v *= fact;
                }
            }
            for idim in 0..dim {
                self.first_constr[idim][derive as usize] = result[idim];
            }
        }
        for derive in (0..=self.niv_constr).rev() {
            if eval(last, derive, &mut result) != 0 {
                return;
            }
            if derive >= 1 {
                for v in &mut result {
                    *v *= fact;
                }
            }
            if derive == 2 {
                for v in &mut result {
                    *v *= fact;
                }
            }
            for idim in 0..dim {
                self.last_constr[idim][derive as usize] = result[idim];
            }
        }
        if !hermite_interpolate(
            dim,
            -1.0,
            1.0,
            self.niv_constr,
            self.niv_constr,
            &self.first_constr,
            &self.last_constr,
            &mut self.coeff,
        ) {
            return;
        }
        let alin = (last - first) / 2.0;
        let blin = (last + first) / 2.0;
        let mut fti = vec![0.0; dim];
        let mut rpti = vec![0.0; dim];
        let mut rmti = vec![0.0; dim];
        let half = (self.nb_gauss / 2) as usize;
        let mut i_idim = dim;
        for i in 1..=half {
            let ti = self.tab_points[i];
            let tip = alin * ti + blin;
            if eval(tip, 0, &mut fti) != 0 {
                return;
            }
            for idim in 0..dim {
                self.som[i_idim] = fti[idim];
                self.dif[i_idim] = fti[idim];
                i_idim += 1;
            }
        }
        i_idim = dim;
        for i in 1..=half {
            let ti = self.tab_points[i];
            let tin = -alin * ti + blin;
            if eval(tin, 0, &mut fti) != 0 {
                return;
            }
            eval_poly0(&self.coeff, degree_r, dim, ti, &mut rpti);
            eval_poly0(&self.coeff, degree_r, dim, -ti, &mut rmti);
            for idim in 0..dim {
                self.som[i_idim] += fti[idim] - rpti[idim] - rmti[idim];
                self.dif[i_idim] -= fti[idim] + rpti[idim] - rmti[idim];
                i_idim += 1;
            }
        }
        if self.nb_gauss % 2 == 1 {
            let tip = blin;
            if eval(tip, 0, &mut fti) != 0 {
                return;
            }
            eval_poly0(&self.coeff, degree_r, dim, self.tab_points[0], &mut rpti);
            for idim in 0..dim {
                self.som[idim] = fti[idim] - rpti[idim];
                self.dif[idim] = fti[idim] - rpti[idim];
            }
        }
        for k in (0..=degree_q).step_by(2) {
            for idim in 0..dim {
                let mut sum = 0.0;
                for i in 1..=half {
                    sum += self.tab_weights[i][k as usize]
                        * self.som[i * dim + idim];
                }
                self.coeff[((k + degree_r + 1) as usize) * dim + idim] = sum;
            }
        }
        for k in (1..=degree_q).step_by(2) {
            for idim in 0..dim {
                let mut sum = 0.0;
                for i in 1..=half {
                    sum += self.tab_weights[i][k as usize]
                        * self.dif[i * dim + idim];
                }
                self.coeff[((k + degree_r + 1) as usize) * dim + idim] = sum;
            }
        }
        let mut jac_coeff = vec![0.0; dim * (self.work_degree as usize + 1)];
        let mut new_degree_max = 0i32;
        let mut rang_ss = 0usize;
        let mut rang_jac = 0usize;
        for numss in 0..self.total_num_ss as usize {
            let local_dim = local_dimension[numss] as usize;
            let mut rang_coeff = 0usize;
            let mut rang_dim = 0usize;
            for _k in 0..=self.work_degree as usize {
                for idim in 0..local_dim {
                    jac_coeff[rang_jac + rang_dim + idim] = self.coeff[rang_coeff + rang_ss + idim];
                }
                rang_dim += local_dim;
                rang_coeff += dim;
            }
            let (new_degree, _) = self.jac.reduce_degree(
                local_dimension[numss],
                max_degree,
                local_tol[numss],
                &jac_coeff[rang_jac..],
            );
            if new_degree > new_degree_max {
                new_degree_max = new_degree;
            }
            rang_ss += local_dim;
            rang_jac += (self.work_degree as usize + 1) * local_dim;
        }
        self.max_error = vec![0.0; self.total_num_ss as usize];
        self.average_error = vec![0.0; self.total_num_ss as usize];
        let mut _rang_ss = 0usize;
        let mut rang_jac = 0usize;
        for numss in 0..self.total_num_ss as usize {
            let local_dim = local_dimension[numss] as usize;
            self.max_error[numss] =
                self.jac
                    .max_error(local_dimension[numss], &jac_coeff[rang_jac..], new_degree_max);
            self.average_error[numss] = self.jac.average_error(
                local_dimension[numss],
                &jac_coeff[rang_jac..],
                new_degree_max,
            );
            _rang_ss += local_dim;
            rang_jac += (self.work_degree as usize + 1) * local_dim;
        }
        self.degree = new_degree_max;
        self.done = true;
    }

    pub fn is_done(&self) -> bool {
        self.done
    }

    pub fn degree(&self) -> i32 {
        self.degree
    }

    pub fn coefficients(&self) -> &[f64] {
        &self.coeff
    }

    pub fn max_error(&self, index: usize) -> f64 {
        self.max_error[index]
    }

    pub fn average_error(&self, index: usize) -> f64 {
        self.average_error[index]
    }
}
