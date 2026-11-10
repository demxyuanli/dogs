//! `ComputePatches` and the iso approximation loop (`cxx:340-476`).

use crate::Surface;

use super::super::context::App2VarContext;
use super::super::criterion::{Criterion, CriterionPatch, CriterionRepartition, CriterionType};

use super::*;

/// `InitGrid(1)` then `ComputeConstraints`, then the corner tables of
/// `Patch::Discretise`. Cuts call `UpdateInU` and `UpdateInV`. The next
/// `ConvertBS` returns the B-spline.
///
/// This is the criterion-free `ComputePatches`
/// (`AdvApp2Var_ApproxAFunc2Var.cxx:340-476`).
pub fn approximate_until_cdi(
    surf: &dyn Surface,
    u0: f64,
    u1: f64,
    v0: f64,
    v1: f64,
    ctx: &App2VarContext,
    max_patches: i32,
) -> Result<super::super::ApproxAFunc2VarResult, String> {
    compute_patches(surf, u0, u1, v0, v1, ctx, max_patches, None)
}

/// `ComputePatches(UChoice, VChoice, Func, Crit)`
/// (`AdvApp2Var_ApproxAFunc2Var.cxx:480-644`): the same loop as
/// [`approximate_until_cdi`], but every square asks the criterion for a
/// `CritValue` and `CutSense(Crit, NumDec)` decides, instead of
/// `myCutSense`.
pub fn approximate_until_cdi_with_criterion(
    surf: &dyn Surface,
    u0: f64,
    u1: f64,
    v0: f64,
    v1: f64,
    ctx: &App2VarContext,
    max_patches: i32,
    crit: &dyn Criterion,
) -> Result<super::super::ApproxAFunc2VarResult, String> {
    compute_patches(surf, u0, u1, v0, v1, ctx, max_patches, Some(crit))
}

/// `InitGrid(NbInt)` (`AdvApp2Var_ApproxAFunc2Var.cxx:239-331`): the seeded
/// grid is `InitGrid(1)`, then `NbInt - 1` regular cuts in `U` and `V`.
fn init_grid(
    u0: f64,
    u1: f64,
    v0: f64,
    v1: f64,
    ctx: &App2VarContext,
    nb_int: i32,
) -> Result<(super::super::framework::Framework, super::super::network::Network), String> {
    let mut frame = super::super::framework::Framework::init_grid(u0, u1, v0, v1, ctx.ord_u, ctx.ord_v);
    let mut net = super::super::network::Network::init(u0, u1, v0, v1, ctx.ord_u, ctx.ord_v);
    let deltu = (u1 - u0) / f64::from(nb_int);
    let deltv = (v1 - v0) / f64::from(nb_int);
    for iint in 1..nb_int {
        let cut_u = u0 + f64::from(iint) * deltu;
        let cut_v = v0 + f64::from(iint) * deltv;
        net.update_in_u(cut_u)?;
        frame.update_in_u(cut_u)?;
        net.update_in_v(cut_v)?;
        frame.update_in_v(cut_v)?;
    }
    Ok((frame, net))
}

/// `CutSense(Crit, NumDec)` (`AdvApp2Var_Patch.cxx:1043-1061`).
///
/// A relative criterion keeps the approximation's own answer while the square
/// is not approximated; otherwise the criterion decides between "no cut" (`0`)
/// and `NumDec`.
fn patch_cut_sense(
    crit: &dyn Criterion,
    patch: &CriterionPatch<'_>,
    num_dec: i32,
    cut_sense: i32,
) -> i32 {
    let relative = crit.criterion_type() == CriterionType::Relative;
    if relative && !patch.is_approximated {
        cut_sense
    } else if crit.is_satisfied(patch) {
        0
    } else {
        num_dec
    }
}

fn compute_patches(
    surf: &dyn Surface,
    u0: f64,
    u1: f64,
    v0: f64,
    v1: f64,
    ctx: &App2VarContext,
    max_patches: i32,
    crit: Option<&dyn Criterion>,
) -> Result<super::super::ApproxAFunc2VarResult, String> {
    // `GeomConvert_ApproxSurface.cxx:381-398`: C2 breaks are recommended,
    // C3 breaks are preferred, and `Weight` stays 5.
    let u_choice = crate::adv_approx::PrefAndRec::new(surf.u_intervals(4), surf.u_intervals(5), 5.0)
        .map_err(|_| "PrefAndRec : Weight is too small".to_string())?;
    let v_choice = crate::adv_approx::PrefAndRec::new(surf.v_intervals(4), surf.v_intervals(5), 5.0)
        .map_err(|_| "PrefAndRec : Weight is too small".to_string())?;
    // `AdvApp2Var_ApproxAFunc2Var.cxx:589-594` restarts the grid with
    // `InitGrid(++NbInt)`; `NbInt` starts at 1 (`:234`).
    let mut nb_int = 1;
    // `ComputeCritError` (`cxx:868-887`) takes the maximum `CritValue` over
    // the squares that survive; a restart drops every stored square, so the
    // running maximum restarts with it.
    let mut crit_max = 0.0f64;
    'grid: loop {
        let (mut frame, mut net) = init_grid(u0, u1, v0, v1, ctx, nb_int)?;
        loop {
            frame.compute_constraints(
                &mut net,
                max_patches,
                &u_choice,
                &v_choice,
                |iso, begin, end| approx_iso(surf, ctx, iso, begin, end),
            )?;
            let Some(index) = net.first_not_approx() else {
                let errors = net.compute_3d_errors(ctx.nb_3d, &ctx.tol_3d, &ctx.tol_3d_front)?;
                let surface = super::super::convert_bs::convert_bs(&mut net, ctx)?;
                return Ok(super::super::ApproxAFunc2VarResult {
                    surface,
                    is_done: errors.done,
                    max_error: errors.max_error,
                    average_error: errors.average_error,
                    u_front_error: errors.u_front_error,
                    v_front_error: errors.v_front_error,
                    criterion_error: crit_max,
                });
            };
            let (pu0, pu1, pv0, pv1) = net.domain(index)?;
            let corner = super::super::discretise::discretise_corners(
                &frame,
                ctx,
                surf,
                pu0,
                pu1,
                pv0,
                pv1,
                max_patches,
                net.nb_patch_in_u(),
                net.nb_patch_in_v(),
                crit.is_some_and(|c| c.criterion_type() == CriterionType::Absolute),
            )?;
            // `AdvApp2Var_ApproxAFunc2Var.cxx:574-577`.
            let mut num_dec = corner.num_dec;
            if num_dec >= 3 {
                num_dec -= 2;
            }
            // `:579-587`: the criterion's value on this square.
            let (decision, regular, crit_value) = match crit {
                Some(c) => {
                    let (equation, nb_coeff_u, nb_coeff_v) = match &corner.step.finished {
                        Some(done) => (done.equation.as_slice(), done.nb_coeff_u, done.nb_coeff_v),
                        None => (&[][..], 0, 0),
                    };
                    let mut patch = CriterionPatch {
                        u0: pu0,
                        u1: pu1,
                        v0: pv0,
                        v1: pv1,
                        nb_coeff_u,
                        nb_coeff_v,
                        equation,
                        crit_value: 0.0,
                        is_approximated: corner.step.finished.is_some(),
                    };
                    if corner.step.finished.is_some() {
                        c.value(&mut patch, ctx);
                    }
                    // `m1` (`:485`, `:581-585`) accumulates the maximum
                    // `CritValue` and is never read again in OCCT.
                    let decision = patch_cut_sense(c, &patch, num_dec, corner.step.sense);
                    (
                        decision,
                        c.repartition() == CriterionRepartition::Regular,
                        patch.crit_value,
                    )
                }
                None => (corner.step.sense, false, 0.0),
            };
            // `:588-594`: a regular criterion re-partitions the whole square
            // instead of splitting the current one.
            if regular && decision > 0 {
                nb_int += 1;
                crit_max = 0.0;
                continue 'grid;
            }
            // `:596-638`.
            if decision == 0 {
                if let Some(done) = corner.step.finished {
                    if crit.is_some() {
                        crit_max = crit_max.max(crit_value);
                    }
                    net.store_result(index, done)?;
                    continue;
                }
                return Err("AdvApp2Var_ApproxAFunc2Var : Surface Approximation Error".into());
            }
            match decision {
                1 => {
                    net.update_in_u(corner.u_dec)?;
                    frame.update_in_u(corner.u_dec)?;
                }
                2 => {
                    net.update_in_v(corner.v_dec)?;
                    frame.update_in_v(corner.v_dec)?;
                }
                3 => {
                    net.update_in_u(corner.u_dec)?;
                    frame.update_in_u(corner.u_dec)?;
                    net.update_in_v(corner.v_dec)?;
                    frame.update_in_v(corner.v_dec)?;
                }
                _ => {
                    return Err("AdvApp2Var_ApproxAFunc2Var : Surface Approximation Error".into());
                }
            }
        }
    }
}

fn approx_iso(
    surf: &dyn Surface,
    ctx: &App2VarContext,
    iso: &super::super::framework::GridIso,
    begin: &mut super::super::node::App2VarNode,
    end: &mut super::super::node::App2VarNode,
) -> Result<super::super::framework::IsoStatus, String> {
    let (extrem_order, deriv_order) = if iso.isofav == 2 {
        (ctx.ord_u, ctx.ord_v)
    } else {
        (ctx.ord_v, ctx.ord_u)
    };
    let mut errors = Vec::new();
    let mut somtab = Vec::new();
    let mut diftab = Vec::new();
    let mut polynom = Vec::new();
    let mut ncoeff = 0;
    let mut iercod = 0;
    for iapp in 0..=deriv_order {
        let cut = mma2fnc(
            surf,
            ctx,
            [iso.u0, iso.u1, iso.v0, iso.v1],
            iso.constant,
            iso.isofav,
            extrem_order,
            iapp,
            iso.position,
        )?;
        if cut.iercod > 0 {
            return Err("AdvApp2Var_ApproxAFunc2Var : Curve Approximation Error".into());
        }
        if cut.iercod == -1 {
            iercod = -1;
        }
        super::super::node::store_ends(
            begin,
            end,
            iapp,
            extrem_order,
            iso.isofav,
            &cut.contr1,
            &cut.contr2,
        )?;
        errors.push((cut.errmax, cut.errmoy));
        somtab.extend_from_slice(&cut.somtab);
        diftab.extend_from_slice(&cut.diftab);
        polynom.extend_from_slice(&cut.courbe);
        if ncoeff <= cut.ncoeff {
            ncoeff = cut.ncoeff;
        }
    }
    Ok(super::super::framework::IsoStatus {
        approximated: iercod == 0,
        has_result: iercod == 0 || iercod == -1,
        errors,
        somtab,
        diftab,
        polynom,
        nb_coeff: ncoeff,
    })
}
