//! `BSplSLib::PrepareEval` and non-rational `BSplSLib::DN`.
//! Source: `BSplSLib.cxx:313-733` and `cxx:1519-1605`.
//! Pole/knot slices are 0-based; OCCT `Array` Lower is treated as 1.

use crate::bspl::bohm::bohm;
use crate::bspl::build_knots::build_knots;
use crate::bspl::knots::pole_index;
use crate::bspl::locate::locate_parameter;
use crate::bspl::rational_derivative::rational_derivative;
use crate::gp::{GpPnt, GpVec};
use crate::precision::epsilon;

/// Result of `PrepareEval`.
pub struct PreparedEval {
    pub u1: f64,
    pub u2: f64,
    pub d1: i32,
    pub d2: i32,
    pub rational: bool,
    pub ufirst: bool,
    pub poles: Vec<f64>,
    pub knots1: Vec<f64>,
    pub knots2: Vec<f64>,
}

fn pole_at(poles: &[Vec<GpPnt>], ip: i32, jp: i32) -> GpPnt {
    let nr = poles.len() as i32;
    let nc = poles.first().map(|r| r.len() as i32).unwrap_or(0);
    let mut i = ip;
    let mut j = jp;
    if i < 1 {
        i = nr;
    }
    if i > nr {
        i = 1;
    }
    if j < 1 {
        j = nc;
    }
    if j > nc {
        j = 1;
    }
    poles[(i - 1) as usize][(j - 1) as usize]
}

fn weight_at(weights: &[Vec<f64>], ip: i32, jp: i32) -> f64 {
    let nr = weights.len() as i32;
    let nc = weights.first().map(|r| r.len() as i32).unwrap_or(0);
    let mut i = ip;
    let mut j = jp;
    if i < 1 {
        i = nr;
    }
    if i > nr {
        i = 1;
    }
    if j < 1 {
        j = nc;
    }
    if j > nc {
        j = 1;
    }
    weights[(i - 1) as usize][(j - 1) as usize]
}

fn local_rational(weights: &[Vec<f64>], uindex: i32, vindex: i32, udeg: i32, vdeg: i32) -> bool {
    let p_lower_row = 1i32;
    let p_upper_row = weights.len() as i32;
    let p_lower_col = 1i32;
    let p_upper_col = weights.first().map(|r| r.len() as i32).unwrap_or(0);
    let mut ip = p_lower_row + uindex;
    let mut jp = p_lower_col + vindex;
    if ip < p_lower_row {
        ip = p_upper_row;
    }
    if jp < p_lower_col {
        jp = p_upper_col;
    }
    let w = weight_at(weights, ip, jp);
    let eps = epsilon(w);
    let mut rational = false;
    let mut i = 0;
    while i <= udeg && !rational {
        jp = p_lower_col + vindex;
        if jp < p_lower_col {
            jp = p_upper_col;
        }
        let mut j = 0;
        while j <= vdeg && !rational {
            let mut dw = weight_at(weights, ip, jp) - w;
            if dw < 0.0 {
                dw = -dw;
            }
            rational = dw > eps;
            jp += 1;
            if jp > p_upper_col {
                jp = p_lower_col;
            }
            j += 1;
        }
        ip += 1;
        if ip > p_upper_row {
            ip = p_lower_row;
        }
        i += 1;
    }
    rational
}

/// `PrepareEval` (`BSplSLib.cxx:313-733`). Returns `ufirst` as the C++ `bool`.
pub fn prepare_eval(
    u: f64,
    v: f64,
    u_index: i32,
    v_index: i32,
    u_degree: i32,
    v_degree: i32,
    u_rat: bool,
    v_rat: bool,
    u_per: bool,
    v_per: bool,
    poles: &[Vec<GpPnt>],
    weights: Option<&[Vec<f64>]>,
    u_knots: &[f64],
    v_knots: &[f64],
    u_mults: Option<&[i32]>,
    v_mults: Option<&[i32]>,
) -> PreparedEval {
    let mut rational = u_rat || v_rat;
    let mut uindex = u_index;
    let mut vindex = v_index;
    let uk_lower = 1i32;
    let vk_lower = 1i32;
    let p_lower_row = 1i32;
    let p_upper_row = poles.len() as i32;
    let p_lower_col = 1i32;
    let p_upper_col = poles.first().map(|r| r.len() as i32).unwrap_or(0);

    if u_degree <= v_degree {
        let (ui, u1) = if uindex < uk_lower || uindex > u_knots.len() as i32 {
            locate_parameter(u_degree, u_knots, u_mults, u, u_per, uindex)
        } else {
            (uindex, u)
        };
        uindex = ui;
        let (vi, u2) = if vindex < vk_lower || vindex > v_knots.len() as i32 {
            locate_parameter(v_degree, v_knots, v_mults, v, v_per, vindex)
        } else {
            (vindex, v)
        };
        vindex = vi;
        let d1 = u_degree;
        let d2 = v_degree;
        let knots1 = build_knots(u_degree, uindex, u_per, u_knots, u_mults);
        let knots2 = build_knots(v_degree, vindex, v_per, v_knots, v_mults);
        if let Some(m) = u_mults {
            uindex = pole_index(u_degree, uindex, u_per, m);
        } else {
            uindex -= uk_lower + u_degree;
        }
        if let Some(m) = v_mults {
            vindex = pole_index(v_degree, vindex, v_per, m);
        } else {
            vindex -= vk_lower + v_degree;
        }
        if rational {
            rational = match weights {
                Some(w) => local_rational(w, uindex, vindex, u_degree, v_degree),
                None => false,
            };
        }
        let dim = if rational { 4 } else { 3 };
        let mut pack = vec![0.0; (dim * (d1 + 1) * (d2 + 1)) as usize];
        let mut ip = p_lower_row + uindex;
        if ip < p_lower_row {
            ip = p_upper_row;
        }
        let mut slot = 0usize;
        for _i in 0..=d1 {
            let mut jp = p_lower_col + vindex;
            if jp < p_lower_col {
                jp = p_upper_col;
            }
            for _j in 0..=d2 {
                let p = pole_at(poles, ip, jp);
                if rational {
                    let w = weights.map(|ws| weight_at(ws, ip, jp)).unwrap_or(1.0);
                    pack[slot] = p.x() * w;
                    pack[slot + 1] = p.y() * w;
                    pack[slot + 2] = p.z() * w;
                    pack[slot + 3] = w;
                    slot += 4;
                } else {
                    pack[slot] = p.x();
                    pack[slot + 1] = p.y();
                    pack[slot + 2] = p.z();
                    slot += 3;
                }
                jp += 1;
                if jp > p_upper_col {
                    jp = p_lower_col;
                }
            }
            ip += 1;
            if ip > p_upper_row {
                ip = p_lower_row;
            }
        }
        PreparedEval {
            u1,
            u2,
            d1,
            d2,
            rational,
            ufirst: true,
            poles: pack,
            knots1,
            knots2,
        }
    } else {
        let (ui, u2) = if uindex < uk_lower || uindex > u_knots.len() as i32 {
            locate_parameter(u_degree, u_knots, u_mults, u, u_per, uindex)
        } else {
            (uindex, u)
        };
        uindex = ui;
        let (vi, u1) = if vindex < vk_lower || vindex > v_knots.len() as i32 {
            locate_parameter(v_degree, v_knots, v_mults, v, v_per, vindex)
        } else {
            (vindex, v)
        };
        vindex = vi;
        let d2 = u_degree;
        let d1 = v_degree;
        let knots2 = build_knots(u_degree, uindex, u_per, u_knots, u_mults);
        let knots1 = build_knots(v_degree, vindex, v_per, v_knots, v_mults);
        if let Some(m) = u_mults {
            uindex = pole_index(u_degree, uindex, u_per, m);
        } else {
            uindex -= uk_lower + u_degree;
        }
        if let Some(m) = v_mults {
            vindex = pole_index(v_degree, vindex, v_per, m);
        } else {
            vindex -= vk_lower + v_degree;
        }
        if rational {
            rational = match weights {
                Some(w) => local_rational(w, uindex, vindex, u_degree, v_degree),
                None => false,
            };
        }
        let dim = if rational { 4 } else { 3 };
        let mut pack = vec![0.0; (dim * (d1 + 1) * (d2 + 1)) as usize];
        let mut jp = p_lower_col + vindex;
        if jp < p_lower_col {
            jp = p_upper_col;
        }
        let mut slot = 0usize;
        for _i in 0..=d1 {
            let mut ip = p_lower_row + uindex;
            if ip < p_lower_row {
                ip = p_upper_row;
            }
            if !rational && ip > p_upper_row {
                ip = p_lower_row;
            }
            for _j in 0..=d2 {
                let p = pole_at(poles, ip, jp);
                if rational {
                    let w = weights.map(|ws| weight_at(ws, ip, jp)).unwrap_or(1.0);
                    pack[slot] = p.x() * w;
                    pack[slot + 1] = p.y() * w;
                    pack[slot + 2] = p.z() * w;
                    pack[slot + 3] = w;
                    slot += 4;
                } else {
                    pack[slot] = p.x();
                    pack[slot + 1] = p.y();
                    pack[slot + 2] = p.z();
                    slot += 3;
                }
                ip += 1;
                if ip > p_upper_row {
                    ip = p_lower_row;
                }
            }
            jp += 1;
            if jp > p_upper_col {
                jp = p_lower_col;
            }
        }
        PreparedEval {
            u1,
            u2,
            d1,
            d2,
            rational,
            ufirst: false,
            poles: pack,
            knots1,
            knots2,
        }
    }
}

/// `BSplSLib::DN` (`cxx:1519-1605`).
pub fn dn(
    u: f64,
    v: f64,
    nu: i32,
    nv: i32,
    u_index: i32,
    v_index: i32,
    poles: &[Vec<GpPnt>],
    weights: Option<&[Vec<f64>]>,
    u_knots: &[f64],
    v_knots: &[f64],
    u_mults: Option<&[i32]>,
    v_mults: Option<&[i32]>,
    u_degree: i32,
    v_degree: i32,
    u_rat: bool,
    v_rat: bool,
    u_per: bool,
    v_per: bool,
) -> GpVec {
    let mut prep = prepare_eval(
        u, v, u_index, v_index, u_degree, v_degree, u_rat, v_rat, u_per, v_per, poles, weights,
        u_knots, v_knots, u_mults, v_mults,
    );
    let dim = if prep.rational { 4 } else { 3 };
    if !prep.rational && (nu > u_degree || nv > v_degree) {
        return GpVec::new(0.0, 0.0, 0.0);
    }
    let n1 = if prep.ufirst { nu } else { nv };
    let n2 = if prep.ufirst { nv } else { nu };
    bohm(
        prep.u1,
        prep.d1,
        n1,
        &prep.knots1,
        dim * (prep.d2 + 1),
        &mut prep.poles,
    );
    let stride = (dim * (prep.d2 + 1)) as usize;
    let kmax = n1.min(prep.d1);
    for k in 0..=kmax {
        let start = (k as usize) * stride;
        bohm(
            prep.u2,
            prep.d2,
            n2,
            &prep.knots2,
            dim,
            &mut prep.poles[start..],
        );
    }
    if prep.rational {
        let mut ders = [0.0f64; 3];
        rational_derivative(prep.d1, prep.d2, n1, n2, &prep.poles, &mut ders, false);
        return GpVec::new(ders[0], ders[1], ders[2]);
    }
    let off = ((n1 * (prep.d2 + 1) + n2) * dim) as usize;
    GpVec::new(prep.poles[off], prep.poles[off + 1], prep.poles[off + 2])
}
