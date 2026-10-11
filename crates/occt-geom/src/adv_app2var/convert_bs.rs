//! `AdvApp2Var_ApproxAFunc2Var::ConvertBS` (`AdvApp2Var_ApproxAFunc2Var.cxx:889-990`).

use crate::bspline_surface::GeomBSplineSurface;

use super::context::App2VarContext;
use super::network::Network;

pub(super) fn convert_bs(net: &mut Network, ctx: &App2VarContext) -> Result<GeomBSplineSurface, String> {
    if ctx.total_ssp() != 1 || ctx.nb_3d != 1 {
        return Err("AdvApp2Var_Patch::Poles :  SSPIndex out of range".into());
    }
    net.same_degree(ctx.ord_u, ctx.ord_v);
    let nb_u = net.nb_patch_in_u();
    let nb_v = net.nb_patch_in_v();
    let size_eq = (ctx.lim_u)
        .checked_mul(ctx.lim_v)
        .and_then(|n| n.checked_mul(3))
        .ok_or_else(|| convert_bs_line(936))?;
    let n_patch = (nb_u)
        .checked_mul(nb_v)
        .ok_or_else(|| convert_bs_line(935))?;
    let poly_len = (n_patch)
        .checked_mul(size_eq)
        .ok_or_else(|| convert_bs_line(940))?;
    let poly_len = usize::try_from(poly_len).map_err(|_| convert_bs_line(940))?;
    let size_eq = usize::try_from(size_eq).map_err(|_| convert_bs_line(936))?;
    let mut poly = vec![0.0; poly_len];
    let mut counts = Vec::new();
    let mut icf = 0usize;
    if nb_v >= 1 && nb_u >= 1 {
        for j in 1..=nb_v {
            for i in 1..=nb_u {
                let (nu, nv) = net.coeff_counts(i, j)?;
                counts.push((nu, nv));
                let equation = net.equation(i, j)?;
                for ieq in 0..size_eq {
                    let value = *equation.get(ieq).ok_or_else(|| convert_bs_line(957))?;
                    let slot = poly.get_mut(icf).ok_or_else(|| convert_bs_line(957))?;
                    *slot = value;
                    icf += 1;
                }
            }
        }
    }
    let grid = super::grid_polynomial::grid_polynomial_to_poles(
        nb_u,
        nb_v,
        ctx.ord_u,
        ctx.ord_v,
        ctx.lim_u - 1,
        ctx.lim_v - 1,
        &counts,
        &poly,
        &[-1.0, 1.0],
        &[-1.0, 1.0],
        net.u_parameters(),
        net.v_parameters(),
    )?;
    GeomBSplineSurface::from_poles_knots_mults(
        grid.poles,
        grid.u_knots,
        grid.v_knots,
        grid.u_mults,
        grid.v_mults,
        grid.u_degree as usize,
        grid.v_degree as usize,
    )
}

fn convert_bs_line(line: i32) -> String {
    format!("AdvApp2Var_ApproxAFunc2Var::ConvertBS: AdvApp2Var_ApproxAFunc2Var.cxx:{line}")
}
