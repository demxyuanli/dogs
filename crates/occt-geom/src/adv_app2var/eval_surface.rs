//! `GeomConvert_ApproxSurface_Eval::Evaluate`
//! (`GeomConvert_ApproxSurface.cxx:55-310`).
//!
//! Orders 0, 1, and 2 use `D0` / `D1` / `D2`.
//! Order 3 stores `D3UUV` when `UOrder == 2` and `D3UVV` when `UOrder == 1`
//! (`cxx:193-205` and `:283-295`). Order 4 stores `DN` (`cxx:212` and `:302`).

use crate::Surface;

pub(super) fn evaluate(
    surf: &dyn Surface,
    dimension: i32,
    u_bounds: [f64; 2],
    v_bounds: [f64; 2],
    favor_iso: i32,
    const_param: f64,
    parameters: &[f64],
    u_order: i32,
    v_order: i32,
    result: &mut [f64],
) -> Result<i32, String> {
    // `UTrim` / `VTrim` (`:117-118`) restrict the adaptor to these bounds.
    // The trimmed derivatives forward to the basis at the same parameters.
    let _ = (u_bounds, v_bounds);
    let mut error = 0;
    if dimension != 3 {
        error = 1;
    }
    let order = u_order + v_order;
    let n = parameters.len() as i32;
    if favor_iso == 1 {
        let u_par = const_param;
        match order {
            0 => {
                for jpar in 1..=n {
                    let p = surf.d0(u_par, parameters[(jpar - 1) as usize]);
                    store3(result, jpar, dimension, p.x(), p.y(), p.z());
                }
            }
            1 => {
                for jpar in 1..=n {
                    let (_, d1u, d1v) = surf.d1(u_par, parameters[(jpar - 1) as usize]);
                    let v = if u_order == 1 { d1u } else { d1v };
                    store3(result, jpar, dimension, v.x(), v.y(), v.z());
                }
            }
            2 => {
                for jpar in 1..=n {
                    let (_, _, _, d2u, d2v, d2uv) = surf.d2(u_par, parameters[(jpar - 1) as usize]);
                    let v = if u_order == 2 {
                        d2u
                    } else if u_order == 1 {
                        d2uv
                    } else if u_order == 0 {
                        d2v
                    } else {
                        continue;
                    };
                    store3(result, jpar, dimension, v.x(), v.y(), v.z());
                }
            }
            3 => {
                for jpar in 1..=n {
                    store_d3(
                        surf,
                        u_par,
                        parameters[(jpar - 1) as usize],
                        u_order,
                        result,
                        jpar,
                        dimension,
                        "UNPORTED GeomConvert_ApproxSurface_Eval::D3: GeomConvert_ApproxSurface.cxx:193",
                    )?;
                }
            }
            4 => {
                for jpar in 1..=n {
                    store_dn(
                        surf,
                        u_par,
                        parameters[(jpar - 1) as usize],
                        u_order,
                        v_order,
                        result,
                        jpar,
                        dimension,
                        "UNPORTED GeomConvert_ApproxSurface_Eval::DN: GeomConvert_ApproxSurface.cxx:213",
                    )?;
                }
            }
            _ => {
                return Err(
                    "UNPORTED GeomConvert_ApproxSurface_Eval: GeomConvert_ApproxSurface.cxx:133"
                        .into(),
                );
            }
        }
    } else {
        let v_par = const_param;
        match order {
            0 => {
                for jpar in 1..=n {
                    let p = surf.d0(parameters[(jpar - 1) as usize], v_par);
                    store3(result, jpar, dimension, p.x(), p.y(), p.z());
                }
            }
            1 => {
                for jpar in 1..=n {
                    let (_, d1u, d1v) = surf.d1(parameters[(jpar - 1) as usize], v_par);
                    let v = if u_order == 1 { d1u } else { d1v };
                    store3(result, jpar, dimension, v.x(), v.y(), v.z());
                }
            }
            2 => {
                for jpar in 1..=n {
                    let (_, _, _, d2u, d2v, d2uv) = surf.d2(parameters[(jpar - 1) as usize], v_par);
                    let v = if u_order == 2 {
                        d2u
                    } else if u_order == 1 {
                        d2uv
                    } else if u_order == 0 {
                        d2v
                    } else {
                        continue;
                    };
                    store3(result, jpar, dimension, v.x(), v.y(), v.z());
                }
            }
            3 => {
                for jpar in 1..=n {
                    store_d3(
                        surf,
                        parameters[(jpar - 1) as usize],
                        v_par,
                        u_order,
                        result,
                        jpar,
                        dimension,
                        "UNPORTED GeomConvert_ApproxSurface_Eval::D3: GeomConvert_ApproxSurface.cxx:283",
                    )?;
                }
            }
            4 => {
                for jpar in 1..=n {
                    store_dn(
                        surf,
                        parameters[(jpar - 1) as usize],
                        v_par,
                        u_order,
                        v_order,
                        result,
                        jpar,
                        dimension,
                        "UNPORTED GeomConvert_ApproxSurface_Eval::DN: GeomConvert_ApproxSurface.cxx:302",
                    )?;
                }
            }
            _ => {
                return Err(
                    "UNPORTED GeomConvert_ApproxSurface_Eval: GeomConvert_ApproxSurface.cxx:223"
                        .into(),
                );
            }
        }
    }
    Ok(error)
}

fn map_missing(r: Result<occt_core::gp::GpVec, String>, cite: &str) -> Result<occt_core::gp::GpVec, String> {
    match r {
        Ok(v) => Ok(v),
        Err(e) if e == "UNPORTED Geom_Surface::EvalDN" => Err(cite.to_string()),
        Err(e) => Err(e),
    }
}

fn store_d3(
    surf: &dyn Surface,
    u: f64,
    v: f64,
    u_order: i32,
    result: &mut [f64],
    jpar: i32,
    dimension: i32,
    cite: &str,
) -> Result<(), String> {
    if u_order != 2 && u_order != 1 {
        return Ok(());
    }
    let vec = map_missing(surf.try_eval_d3_uv(u, v, u_order, 3 - u_order), cite)?;
    store3(result, jpar, dimension, vec.x(), vec.y(), vec.z());
    Ok(())
}

fn store_dn(
    surf: &dyn Surface,
    u: f64,
    v: f64,
    u_order: i32,
    v_order: i32,
    result: &mut [f64],
    jpar: i32,
    dimension: i32,
    cite: &str,
) -> Result<(), String> {
    let vec = map_missing(surf.try_eval_dn(u, v, u_order, v_order), cite)?;
    store3(result, jpar, dimension, vec.x(), vec.y(), vec.z());
    Ok(())
}

fn store3(result: &mut [f64], jpar: i32, dimension: i32, x: f64, y: f64, z: f64) {
    let base = ((jpar - 1) * dimension) as usize;
    if dimension >= 1 {
        result[base] = x;
    }
    if dimension >= 2 {
        result[base + 1] = y;
    }
    if dimension >= 3 {
        result[base + 2] = z;
    }
}
