//! `AdvApp2Var_ApproxAFunc2Var::ConvertBS` (`AdvApp2Var_ApproxAFunc2Var.cxx:889-990`).

use crate::bspline_surface::GeomBSplineSurface;

use super::context::App2VarContext;
use super::network::Network;

pub(super) fn convert_bs(
    net: &mut Network,
    ctx: &App2VarContext,
) -> Result<Vec<GeomBSplineSurface>, String> {
    if ctx.nb_3d < 1 {
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
    let mut counts = Vec::new();
    if nb_v >= 1 && nb_u >= 1 {
        for j in 1..=nb_v {
            for i in 1..=nb_u {
                let (nu, nv) = net.coeff_counts(i, j)?;
                counts.push((nu, nv));
            }
        }
    }
    let first_3d_coordinate = ctx.nb_1d + 2 * ctx.nb_2d;
    let mut surfaces = Vec::with_capacity(ctx.nb_3d as usize);
    for surface_index in 0..ctx.nb_3d {
        let coordinate = first_3d_coordinate + 3 * surface_index;
        let mut poly = vec![0.0; poly_len];
        let mut output = 0usize;
        for j in 1..=nb_v {
            for i in 1..=nb_u {
                let equation = net.equation(i, j)?;
                for u in 0..ctx.lim_u {
                    for v in 0..ctx.lim_v {
                        let base = coordinate + (v + u * ctx.lim_v) * ctx.total_dimension();
                        for axis in 0..3 {
                            let value = *equation
                                .get(
                                    usize::try_from(base + axis)
                                        .map_err(|_| convert_bs_line(957))?,
                                )
                                .ok_or_else(|| convert_bs_line(957))?;
                            *poly.get_mut(output).ok_or_else(|| convert_bs_line(957))? = value;
                            output += 1;
                        }
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
        surfaces.push(GeomBSplineSurface::from_poles_knots_mults(
            grid.poles,
            grid.u_knots,
            grid.v_knots,
            grid.u_mults,
            grid.v_mults,
            grid.u_degree as usize,
            grid.v_degree as usize,
        )?);
    }
    Ok(surfaces)
}

fn convert_bs_line(line: i32) -> String {
    format!("AdvApp2Var_ApproxAFunc2Var::ConvertBS: AdvApp2Var_ApproxAFunc2Var.cxx:{line}")
}
