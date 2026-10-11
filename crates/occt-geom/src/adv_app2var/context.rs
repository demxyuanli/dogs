//! `AdvApp2Var_Context` (`AdvApp2Var_Context.cxx:23-290`).
//!
//! `lesparam`, the Jacobi maxima (`mma2jmx_`), the Legendre roots
//! (`mma2roo_`), and the C1 Gauss weights (`mmapptt_`) run here.
//! `Iso::MakeApprox` then calls `mma2fnc_` through `mma1cdi_`.

use super::gauss::gauss_weights;
use super::jacobi_max::jacobi_max;
use super::legendre_roots::legendre_roots;

/// `lesparam` (`AdvApp2Var_Context.cxx:22-78`).
fn lesparam(iordre: i32, ncflim: i32, icodeo: i32) -> Result<(i32, i32), String> {
    if icodeo < 0 {
        return Err("AdvApp2Var_Context".into());
    }
    let mut ndgjac = ncflim;
    if icodeo > 0 {
        ndgjac += 9 - (iordre + 1);
        ndgjac += (icodeo - 1) * 10;
    }
    let mut nbpnts = if ndgjac < 8 {
        8
    } else if ndgjac < 10 {
        10
    } else if ndgjac < 20 {
        20
    } else if ndgjac < 30 {
        30
    } else if ndgjac < 40 {
        40
    } else if ndgjac < 50 {
        50
    } else {
        50
    };
    if iordre > -1 {
        nbpnts += 2;
    }
    Ok((nbpnts, ndgjac))
}

/// Parameters gathered by `AdvApp2Var_ApproxAFunc2Var::Init`
/// (`AdvApp2Var_ApproxAFunc2Var.cxx:164-234`).
#[derive(Clone, Debug)]
#[allow(dead_code)] // Gauss tables are read by `mma1jak_` (`AdvApp2Var_ApproxF2var.cxx:1454`)
pub struct App2VarContext {
    pub fav: i32,
    pub ord_u: i32,
    pub ord_v: i32,
    pub lim_u: i32,
    pub lim_v: i32,
    pub nb_1d: i32,
    pub nb_2d: i32,
    pub nb_3d: i32,
    pub nb_u_root: i32,
    pub nb_v_root: i32,
    pub j_deg_u: i32,
    pub j_deg_v: i32,
    pub j_max_u: Vec<f64>,
    pub j_max_v: Vec<f64>,
    pub u_roots: Vec<f64>,
    pub v_roots: Vec<f64>,
    pub u_gauss: Vec<f64>,
    pub v_gauss: Vec<f64>,
    pub internal_tol: Vec<f64>,
    pub frontier_tol: Vec<[f64; 4]>,
    pub cutting_tol: Vec<[f64; 4]>,
    /// `my3DTolerances` (`AdvApp2Var_ApproxAFunc2Var.cxx:73`). Not the halved
    /// internal tolerance.
    pub tol_3d: Vec<f64>,
    /// `my3DTolOnFront` (`AdvApp2Var_ApproxAFunc2Var.cxx:76`).
    pub tol_3d_front: Vec<[f64; 4]>,
}

impl App2VarContext {
    /// `AdvApp2Var_Context::TotalDimension` (`AdvApp2Var_Context.cxx:295-298`).
    pub fn total_dimension(&self) -> i32 {
        self.nb_1d + 2 * self.nb_2d + 3 * self.nb_3d
    }

    /// `AdvApp2Var_Context::TotalNumberSSP` (`AdvApp2Var_Context.cxx:302-305`).
    pub fn total_ssp(&self) -> i32 {
        self.nb_1d + self.nb_2d + self.nb_3d
    }
}

fn h_max_factor(order: i32) -> f64 {
    match order {
        -1 => 0.0,
        0 => 1.0,
        1 => 1.5,
        _ => 1.75,
    }
}

fn append_internal(source: &[f64], count: i32, offset: i32, target: &mut [f64]) {
    for index in 1..=count {
        target[(offset + index - 1) as usize] = source[(index - 1) as usize];
    }
}

fn append_frontier(
    source: &[[f64; 4]],
    count: i32,
    offset: i32,
    frontier: &mut [[f64; 4]],
    cutting: &mut [[f64; 4]],
) {
    for index in 1..=count {
        let at = (offset + index - 1) as usize;
        frontier[at] = source[(index - 1) as usize];
        cutting[at] = [0.0; 4];
    }
}

/// `AdvApp2Var_Context::AdvApp2Var_Context` (`AdvApp2Var_Context.cxx:143-290`).
pub fn context_until_roots(
    fav: i32,
    iu: i32,
    iv: i32,
    nlimu: i32,
    nlimv: i32,
    iprecis: i32,
    nb1: i32,
    nb2: i32,
    nb3: i32,
    tol1: &[f64],
    tol2: &[f64],
    tol3: &[f64],
    tof1: &[[f64; 4]],
    tof2: &[[f64; 4]],
    tof3: &[[f64; 4]],
) -> Result<App2VarContext, String> {
    let mut coeff_u = nlimu;
    if coeff_u < 2 * iu + 2 {
        coeff_u = 2 * iu + 2;
    }
    let (nb_pnt_u, j_deg_u) = lesparam(iu, coeff_u, iprecis)?;
    let j_max_u = jacobi_max(j_deg_u, iu)?;

    let mut coeff_v = nlimv;
    if coeff_v < 2 * iv + 2 {
        coeff_v = 2 * iv + 2;
    }
    let (nb_pnt_v, j_deg_v) = lesparam(iv, coeff_v, iprecis)?;
    let j_max_v = jacobi_max(j_deg_v, iv)?;
    let internal_u = if iu > -1 { nb_pnt_u - 2 } else { nb_pnt_u };
    let internal_v = if iv > -1 { nb_pnt_v - 2 } else { nb_pnt_v };
    let filled_u = legendre_roots(internal_u)?;
    let filled_v = legendre_roots(internal_v)?;
    let mut u_roots = vec![0.0; nb_pnt_u.max(0) as usize];
    let mut v_roots = vec![0.0; nb_pnt_v.max(0) as usize];
    if filled_u.len() > u_roots.len() || filled_v.len() > v_roots.len() {
        return Err("AdvApp2Var_Context::mma2roo_: AdvApp2Var_Context.cxx:225".into());
    }
    u_roots[..filled_u.len()].copy_from_slice(&filled_u);
    v_roots[..filled_v.len()].copy_from_slice(&filled_v);

    let cols_u = j_deg_u - 2 * iu - 1;
    let cols_v = j_deg_v - 2 * iv - 1;
    let u_gauss = gauss_weights(j_deg_u, internal_u, iu, cols_u)?;
    let v_gauss = gauss_weights(j_deg_v, internal_v, iv, cols_v)?;

    let ssp = (nb1 + nb2 + nb3).max(0) as usize;
    let mut internal_tol = vec![0.0; ssp];
    append_internal(tol1, nb1, 0, &mut internal_tol);
    append_internal(tol2, nb2, nb1, &mut internal_tol);
    append_internal(tol3, nb3, nb1 + nb2, &mut internal_tol);
    if iu > -1 || iv > -1 {
        for value in &mut internal_tol {
            *value /= 2.0;
        }
    }
    let mut frontier_tol = vec![[0.0; 4]; ssp];
    let mut cutting_tol = vec![[0.0; 4]; ssp];
    append_frontier(tof1, nb1, 0, &mut frontier_tol, &mut cutting_tol);
    append_frontier(tof2, nb2, nb1, &mut frontier_tol, &mut cutting_tol);
    append_frontier(tof3, nb3, nb1 + nb2, &mut frontier_tol, &mut cutting_tol);
    if iu > -1 || iv > -1 {
        let weight = h_max_factor(iu) * h_max_factor(iv) + h_max_factor(iu) + h_max_factor(iv);
        for index in 0..ssp {
            for slot in 0..4 {
                let tol_min = internal_tol[index] / weight;
                if tol_min < frontier_tol[index][slot] {
                    frontier_tol[index][slot] = tol_min;
                }
                cutting_tol[index][slot] = tol_min;
            }
        }
    }

    Ok(App2VarContext {
        fav,
        ord_u: iu,
        ord_v: iv,
        lim_u: nlimu,
        lim_v: nlimv,
        nb_1d: nb1,
        nb_2d: nb2,
        nb_3d: nb3,
        nb_u_root: nb_pnt_u,
        nb_v_root: nb_pnt_v,
        j_deg_u,
        j_deg_v,
        j_max_u,
        j_max_v,
        u_roots,
        v_roots,
        u_gauss,
        v_gauss,
        internal_tol,
        frontier_tol,
        cutting_tol,
        tol_3d: tol3.to_vec(),
        tol_3d_front: tof3.to_vec(),
    })
}
