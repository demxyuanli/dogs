//! Jacobi polynomial basis used by `AdvApprox_SimpleApprox`.
//! Source: `PLib_JacobiPolynomial.cxx`, `PLib.cxx` (`JacobiParameters`, `NivConstr`).

use std::sync::OnceLock;

use crate::kernel::geomabs::Shape;

const DATA: &str = include_str!("plib_jacobi_data.pxx");
const INVALID: f64 = -999.0;
const MAX_DEGREE: i32 = 30;
const NB: [i32; 9] = [8, 10, 15, 20, 25, 30, 40, 50, 61];

fn parse_array(name: &str) -> Vec<f64> {
    let key = format!("{name}");
    let Some(pos) = DATA.find(&key) else {
        return Vec::new();
    };
    let rest = &DATA[pos + key.len()..];
    let Some(brace) = rest.find('{') else {
        return Vec::new();
    };
    let rest = &rest[brace..];
    let Some(semi) = rest.find(';') else {
        return Vec::new();
    };
    let body = &rest[..semi];
    let mut out = Vec::new();
    let bytes = body.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        let c = bytes[i];
        let start_num = c == b'-' || c == b'+' || c.is_ascii_digit();
        if start_num {
            let start = i;
            i += 1;
            while i < bytes.len() {
                let d = bytes[i];
                if d.is_ascii_digit() || d == b'.' || d == b'e' || d == b'E' || d == b'+' || d == b'-'
                {
                    i += 1;
                } else {
                    break;
                }
            }
            if let Ok(v) = body[start..i].parse::<f64>() {
                out.push(v);
            }
        } else {
            i += 1;
        }
    }
    out
}

fn db() -> &'static JacobiDb {
    static DB: OnceLock<JacobiDb> = OnceLock::new();
    DB.get_or_init(|| JacobiDb {
        weights: [
            parse_array("WeightsDB_C0"),
            parse_array("WeightsDB_C1"),
            parse_array("WeightsDB_C2"),
        ],
        weights0: [
            parse_array("WeightsDB0_C0"),
            parse_array("WeightsDB0_C1"),
            parse_array("WeightsDB0_C2"),
        ],
        max_values: [
            parse_array("MaxValuesDB_C0"),
            parse_array("MaxValuesDB_C1"),
            parse_array("MaxValuesDB_C2"),
        ],
        trans: [
            parse_array("TransMatrix_C0"),
            parse_array("TransMatrix_C1"),
            parse_array("TransMatrix_C2"),
        ],
    })
}

struct JacobiDb {
    weights: [Vec<f64>; 3],
    weights0: [Vec<f64>; 3],
    max_values: [Vec<f64>; 3],
    trans: [Vec<f64>; 3],
}

/// `PLib::NivConstr`.
pub fn niv_constr(order: Shape) -> i32 {
    match order {
        Shape::C0 => 0,
        Shape::C1 => 1,
        Shape::C2 => 2,
        _ => 0,
    }
}

/// `PLib::JacobiParameters` (`PLib.cxx:2049-2189`).
pub fn jacobi_parameters(order: Shape, max_degree: i32, code: i32) -> Result<(i32, i32), ()> {
    let niv = niv_constr(order);
    if max_degree < 2 * niv + 1 {
        return Err(());
    }
    let work_degree = if code >= 1 {
        max_degree + 9
    } else {
        max_degree + 6
    };
    let ipmin = min_gauss(work_degree)?;
    let iwant = match code {
        -5 => 8,
        -4 => 10,
        -3 => 15,
        -2 => 20,
        -1 => 25,
        1 => 30,
        2 => 40,
        3 => 50,
        4 => 61,
        _ => return Err(()),
    };
    Ok((ipmin.max(iwant), work_degree))
}

fn min_gauss(work_degree: i32) -> Result<i32, ()> {
    for &n in &NB {
        if work_degree < n {
            return Ok(n);
        }
    }
    Err(())
}

/// Positive Legendre roots for Gauss quadrature on [-1, 1].
/// `PLib_JacobiPolynomial::Points` (`cxx:70-100`).
fn gauss_positive_half(n: i32) -> Vec<f64> {
    let (x, _) = gauss_legendre_m1p1(n as usize);
    let half = n as usize / 2;
    x[half..].to_vec()
}

fn gauss_legendre_m1p1(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut x = vec![0.0; n];
    let mut w = vec![0.0; n];
    let m = (n + 1) / 2;
    for i in 1..=m {
        let mut z = (std::f64::consts::PI * (i as f64 - 0.25) / (n as f64 + 0.5)).cos();
        let mut pp = 0.0;
        loop {
            let mut p1 = 1.0;
            let mut p2 = 0.0;
            for j in 1..=n {
                let p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j as f64 - 1.0) * z * p2 - (j as f64 - 1.0) * p3) / j as f64;
            }
            pp = n as f64 * (z * p1 - p2) / (z * z - 1.0);
            let z1 = z;
            z = z1 - p1 / pp;
            if (z - z1).abs() < 1e-15 {
                break;
            }
        }
        x[i - 1] = -z;
        x[n - i] = z;
        w[i - 1] = 2.0 / ((1.0 - z * z) * pp * pp);
        w[n - i] = w[i - 1];
    }
    (x, w)
}

/// `PLib_JacobiPolynomial`.
#[derive(Clone, Copy)]
pub struct JacobiPolynomial {
    work_degree: i32,
    niv_constr: i32,
    degree: i32,
}

impl JacobiPolynomial {
    pub fn new(work_degree: i32, order: Shape) -> Result<Self, ()> {
        let niv = niv_constr(order);
        let degree = work_degree - 2 * (niv + 1);
        if degree < 0 || degree > MAX_DEGREE {
            return Err(());
        }
        Ok(Self {
            work_degree,
            niv_constr: niv,
            degree,
        })
    }

    pub fn work_degree(&self) -> i32 {
        self.work_degree
    }

    pub fn niv_constr(&self) -> i32 {
        self.niv_constr
    }

    pub fn degree(&self) -> i32 {
        self.degree
    }

    /// `Points`: TabPoints[0] unused for even n; [1..=n/2] positive increasing roots.
    pub fn points(&self, nb_gauss: i32) -> Result<Vec<f64>, ()> {
        if !NB.contains(&nb_gauss) || nb_gauss <= self.degree {
            return Err(());
        }
        let pos = gauss_positive_half(nb_gauss);
        let half = (nb_gauss / 2) as usize;
        let mut tab = vec![INVALID; half + 1];
        if nb_gauss % 2 == 1 {
            tab[0] = 0.0;
        }
        for i in 0..half {
            tab[i + 1] = pos[i];
        }
        Ok(tab)
    }

    /// `Weights` (`cxx:104-186`). Row 0 is Gauss t=0 (odd n) or INVALID (even n).
    /// `tab[i][j]` = weight at positive root i for Jacobi degree j.
    pub fn weights(&self, nb_gauss: i32) -> Result<Vec<Vec<f64>>, ()> {
        if !NB.contains(&nb_gauss) || nb_gauss <= self.degree {
            return Err(());
        }
        let niv = self.niv_constr as usize;
        let min_degree = 2 * (self.niv_constr + 1);
        let mut ptr = 0usize;
        let bump = |n: i32, ptr: &mut usize| {
            if n % 2 == 0 {
                *ptr += (n * (n - min_degree) / 2) as usize;
            } else {
                *ptr += (((n - 1) / 2) * (n - min_degree)) as usize;
            }
        };
        if nb_gauss > 8 {
            bump(8, &mut ptr);
        }
        if nb_gauss > 10 {
            bump(10, &mut ptr);
        }
        if nb_gauss > 15 {
            bump(15, &mut ptr);
        }
        if nb_gauss > 20 {
            bump(20, &mut ptr);
        }
        if nb_gauss > 25 {
            bump(25, &mut ptr);
        }
        if nb_gauss > 30 {
            bump(30, &mut ptr);
        }
        if nb_gauss > 40 {
            bump(40, &mut ptr);
        }
        if nb_gauss > 50 {
            bump(50, &mut ptr);
        }
        let src = &db().weights[niv];
        let half = (nb_gauss / 2) as usize;
        let cols = (self.degree + 1) as usize;
        let mut tab = vec![vec![0.0; cols]; half + 1];
        let mut p = ptr;
        for j in 0..cols {
            for i in 1..=half {
                tab[i][j] = *src.get(p).unwrap_or(&0.0);
                p += 1;
            }
        }
        if nb_gauss % 2 == 1 {
            let mut p0 = 0usize;
            if nb_gauss > 15 {
                p0 += ((15 - 1 - min_degree) / 2 + 1) as usize;
            }
            if nb_gauss > 25 {
                p0 += ((25 - 1 - min_degree) / 2 + 1) as usize;
            }
            let s0 = &db().weights0[niv];
            let mut k = p0;
            let mut j = 0usize;
            while j < cols {
                tab[0][j] = *s0.get(k).unwrap_or(&0.0);
                k += 1;
                j += 2;
            }
        } else {
            for j in 0..cols {
                tab[0][j] = INVALID;
            }
        }
        Ok(tab)
    }

    pub fn max_value(&self) -> Vec<f64> {
        let src = &db().max_values[self.niv_constr as usize];
        let n = (self.degree + 2) as usize;
        let mut t = vec![0.0; n];
        for i in 0..n {
            t[i] = *src.get(i).unwrap_or(&0.0);
        }
        t
    }

    /// `MaxError` (`cxx:201-227`). `jac` is laid out degree-major, dimension-minor.
    pub fn max_error(&self, dimension: i32, jac: &[f64], new_degree: i32) -> f64 {
        let tab_max = self.max_value();
        let beg = 2 * (self.niv_constr + 1);
        let cut = beg.max(new_degree + 1);
        let dim = dimension as usize;
        let mut err = vec![0.0; dim];
        for d in 0..dim {
            for c in cut..=self.work_degree {
                let coeff = jac[(c as usize) * dim + d];
                let basis = tab_max[(c - beg) as usize];
                err[d] += coeff.abs() * basis;
            }
        }
        err.iter().map(|v| v * v).sum::<f64>().sqrt()
    }

    /// `ReduceDegree` (`cxx:231-289`).
    pub fn reduce_degree(
        &self,
        dimension: i32,
        max_degree: i32,
        tol: f64,
        jac: &[f64],
    ) -> (i32, f64) {
        let idx = 2 * (self.niv_constr + 1) - 1;
        let cut = idx + 1;
        let tab_max = self.max_value();
        let dim = dimension as usize;
        let mut err_dim = vec![0.0; dim];
        let mut new_degree = idx;
        let mut max_error = 0.0;
        for i in (cut..=self.work_degree).rev() {
            let off = (i as usize) * dim;
            for d in 0..dim {
                err_dim[d] += jac[off + d].abs() * tab_max[(i - cut) as usize];
            }
            let error = err_dim.iter().map(|v| v * v).sum::<f64>().sqrt();
            if error > tol && i <= max_degree {
                new_degree = i;
                break;
            }
            max_error = error;
        }
        if new_degree == idx {
            const EPS: f64 = 1.0e-9;
            new_degree = 0;
            for i in (1..=idx).rev() {
                let off = (i as usize) * dim;
                let mut bid = 0.0;
                for d in 0..dim {
                    bid += jac[off + d].abs();
                }
                if bid > EPS {
                    new_degree = i;
                    break;
                }
            }
        }
        (new_degree, max_error)
    }

    /// `AverageError` (`cxx:293-312`).
    pub fn average_error(&self, dimension: i32, jac: &[f64], new_degree: i32) -> f64 {
        let cut = (2 * (self.niv_constr + 1) + 1).max(new_degree + 1);
        let dim = dimension as usize;
        let mut avg = 0.0;
        for d in 0..dim {
            for i in cut..=self.degree {
                let c = jac[(i as usize) * dim + d];
                avg += c * c;
            }
        }
        (avg / 2.0).sqrt()
    }

    /// `ToCoefficients` (`cxx:316-372`).
    pub fn to_coefficients(&self, dimension: i32, degree: i32, jac: &[f64]) -> Vec<f64> {
        const MAX_M: i32 = MAX_DEGREE + 1;
        let dim = dimension as usize;
        let half = (degree / 2) as usize;
        let double_dim = 2 * dim;
        let tr = &db().trans[self.niv_constr as usize];
        let n = ((degree + 1) as usize) * dim;
        let mut coeff = vec![0.0; n.max(1)];
        for i in 0..=half {
            let ptr_idx = (i as i32 * MAX_M - ((i as i32 + 1) * i as i32) / 2) as usize;
            let coeff_off = double_dim * i;
            for d in 0..dim {
                let mut value = 0.0;
                for j in i..=half {
                    value += tr.get(ptr_idx + j).copied().unwrap_or(0.0)
                        * jac.get(double_dim * j + d).copied().unwrap_or(0.0);
                }
                if coeff_off + d < coeff.len() {
                    coeff[coeff_off + d] = value;
                }
            }
        }
        if degree == 0 {
            return coeff;
        }
        let tr_odd = MAX_M * (MAX_M + 1) / 2;
        let half_m1 = ((degree - 1) / 2) as usize;
        for i in 0..=half_m1 {
            let ptr_idx = (i as i32 * MAX_M - ((i as i32 + 1) * i as i32) / 2) as usize;
            let base = (2 * i + 1) * dim;
            let j_base = (2 * i + 1) * dim;
            for d in 0..dim {
                let mut value = 0.0;
                let mut jj = j_base + d;
                for j in i..=half_m1 {
                    value += tr
                        .get((tr_odd as usize) + ptr_idx + j)
                        .copied()
                        .unwrap_or(0.0)
                        * jac.get(jj).copied().unwrap_or(0.0);
                    jj += double_dim;
                    let _ = j;
                }
                if base + d < coeff.len() {
                    coeff[base + d] = value;
                }
            }
        }
        coeff
    }
}
