//! `ChFi3d_FilBuilder::ExtentOneCorner` / `ExtentTwoCorner` /
//! `ExtentThreeCorner` — the pre-analysis extensions of the guidelines at a
//! vertex, fillet flavour.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_FilBuilder.cxx:2440-2591`
//! (OCCT 8.0.0). The `ChFi3d_ChBuilder` bodies (`ChFi3d_ChBuilder.cxx:1953-2110`)
//! are *not* the same code — the chamfer `ExtentOneCorner` has no
//! `ChFiDS_FreeBoundary` branch and its `ExtentTwoCorner` keeps a two-element
//! `Stripe[2]` / `Iedge[2]` state — so this module is the fillet-only
//! translation and does not serve the chamfer subclass.
//!
//! Translation notes:
//! - `Spine->SetFirstTgt` / `SetLastTgt` raise `Standard_Failure` on a periodic
//!   guideline (`ChFiDS_Spine.cxx:353-369`, `:373-392`); the port answers
//!   `None` and these bodies report `Err`, like `PerformSetOfKPart` does,
//! - `Spine->LastParameter(Spine->NbEdges())` answers `None` when the abscissa
//!   is not built yet (`ChFiDS_Spine.cxx:248-251`); OCCT would return an
//!   undefined `double` there, so the port stops instead of inventing `0.`,
//! - `ChFi3d_IndexOfSurfData` / `ChFi3d_EdgeFromV1` raise
//!   `Standard_ConstructionError` on a vertex that is on neither end of the
//!   spine; the port propagates their `Err`,
//! - `occ::handle` comparison `Stripe == itel.Value()` is OCCT handle equality,
//!   i.e. [`Arc::ptr_eq`] on the port's handles.

use std::sync::Arc;

use crate::chfi3d::builder_0_helpers::{ch_fi3d_edge_from_v1, ch_fi3d_index_of_surf_data};
use crate::chfi_ds::{SpineHandle, State as SpineState, StripeHandle, TypeOfConcavity};
use crate::shape::TopoShape;

/// `Coeff = 0.5` of `ExtentOneCorner` (`cxx:2447`).
const ONE_COEFF: f64 = 0.5;
/// `Coeff = 0.3` of `ExtentTwoCorner` (`cxx:2478`).
const TWO_COEFF: f64 = 0.3;
/// `Coeff = 0.1` of `ExtentThreeCorner` (`cxx:2553`).
const THREE_COEFF: f64 = 0.1;

/// The `ChFiDS_Spine` of `stripe`, the `S->Spine()` / `Stripe->Spine()` of the
/// three bodies.
fn spine_of(stripe: &StripeHandle) -> Result<SpineHandle, String> {
    let guard = stripe
        .read()
        .map_err(|_| "ExtentCorner: poisoned stripe".to_string())?;
    guard
        .spine()
        .cloned()
        .ok_or_else(|| "ExtentCorner: the stripe has no spine".to_string())
}

/// `Spine->LastParameter(Spine->NbEdges())`.
fn spine_end_parameter(spine: &SpineHandle) -> Result<f64, String> {
    let guard = spine
        .read()
        .map_err(|_| "ExtentCorner: poisoned spine".to_string())?;
    guard
        .last_parameter_of(guard.nb_edges() as i32)
        .ok_or_else(|| "ExtentCorner: ChFiDS_Spine::LastParameter(NbEdges)".to_string())
}

/// `ExtentOneCorner(V, S)` (`cxx:2440-2469`).
///
/// The guideline of `S` is extended by half its length towards `V`; a free
/// boundary doubles the coefficient so the walk goes to the end of the
/// support (`cxx:2454-2457`).
pub(crate) fn extent_one_corner(v: &TopoShape, s: &StripeHandle) -> Result<(), String> {
    // `int Sens = 0; double Coeff = 0.5;` (`cxx:2446-2447`).
    let mut sens = 0;
    let mut coeff = ONE_COEFF;
    // `occ::handle<ChFiDS_Spine> Spine = S->Spine();` (`cxx:2448`).
    let spine = spine_of(s)?;
    // `ChFi3d_IndexOfSurfData(V, S, Sens);` (`cxx:2449`).
    ch_fi3d_index_of_surf_data(v, s, &mut sens)?;
    // `double dU = Spine->LastParameter(Spine->NbEdges());` (`cxx:2450`).
    let d_u = spine_end_parameter(&spine)?;
    {
        let guard = spine
            .read()
            .map_err(|_| "ExtentOneCorner: poisoned spine".to_string())?;
        let is_first = sens == 1;
        // `if (Spine->IsTangencyExtremity((Sens == 1))) return;`
        // (`cxx:2451-2454`).
        if guard.is_tangency_extremity(is_first) {
            return Ok(());
        }
        // `if (Spine->Status((Sens == 1)) == ChFiDS_FreeBoundary) Coeff *= 2;`
        // (`cxx:2456-2459`).
        if guard.status(is_first) == SpineState::FreeBoundary {
            coeff *= 2.0;
        }
    }
    let mut guard = spine
        .write()
        .map_err(|_| "ExtentOneCorner: poisoned spine".to_string())?;
    if sens == 1 {
        // `Spine->SetFirstParameter(-dU * Coeff); Spine->SetFirstTgt(0.);`
        // (`cxx:2461-2465`).
        guard.set_first_parameter(-d_u * coeff);
        guard
            .set_first_tgt(0.0)
            .ok_or_else(|| "ExtentOneCorner: SetFirstTgt on a periodic spine".to_string())?;
    } else {
        // `Spine->SetLastParameter(dU * (1. + Coeff)); Spine->SetLastTgt(dU);`
        // (`cxx:2466-2470`).
        guard.set_last_parameter(d_u * (1.0 + coeff));
        guard
            .set_last_tgt(d_u)
            .ok_or_else(|| "ExtentOneCorner: SetLastTgt on a periodic spine".to_string())?;
    }
    Ok(())
}

/// `ExtentTwoCorner(V, LS)` (`cxx:2474-2543`).
///
/// A common extension `Eval` is computed from every stripe of `LS`
/// (`cxx:2486-2514`) and then applied to each of them on the side of `V`
/// (`cxx:2517-2542`).
pub(crate) fn extent_two_corner(v: &TopoShape, ls: &[StripeHandle]) -> Result<(), String> {
    // `double Coeff = 0.3, Eval = 0.0, dU, rad;` / `bool FF = true;`
    // (`cxx:2478-2482`).
    let coeff = TWO_COEFF;
    let mut eval = 0.0f64;
    // `occ::handle<ChFiDS_Stripe> Stripe;` (`cxx:2483`) survives the first loop
    // and is read by the `!FF && Stripe == itel.Value()` test of the second one
    // (`cxx:2520`), so it starts as the last stripe of `LS`.
    let mut stripe_prev: Option<StripeHandle> = None;

    // A value of symmetric extension is calculated (`cxx:2486-2514`).
    for stripe in ls {
        stripe_prev = Some(stripe.clone());
        let spine = spine_of(stripe)?;
        // `dU = Spine->LastParameter(Spine->NbEdges()) * Coeff;` (`cxx:2490`).
        let d_u = spine_end_parameter(&spine)? * coeff;
        // `occ::handle<ChFiDS_FilSpine> fsp = down_cast...;`
        // `if (fsp->IsConstant()) rad = fsp->Radius();` (`cxx:2491-2497`).
        let (is_constant, radius) = {
            let guard = spine
                .read()
                .map_err(|_| "ExtentTwoCorner: poisoned spine".to_string())?;
            (guard.is_fil_constant(), guard.fil_radius())
        };
        let rad = if is_constant {
            radius.ok_or_else(|| {
                "ExtentTwoCorner: ChFiDS_FilSpine::Radius raised Standard_DomainError".to_string()
            })?
        } else {
            // `TopoDS_Edge E = ChFi3d_EdgeFromV1(V, itel.Value(), Sens);`
            // `rad = MaxRad(fsp, E);` (`cxx:2498-2499`).
            let mut sens = 0;
            let e = ch_fi3d_edge_from_v1(v, stripe, &mut sens)?;
            super::statics::max_rad(&spine, &e)?
        };
        // `rad *= 1.5; if (rad > dU) dU = rad; if (dU > Eval) Eval = dU;`
        // (`cxx:2506-2513`).
        let d_u = (rad * 1.5).max(d_u);
        if d_u > eval {
            eval = d_u;
        }
    }

    // One applies (`cxx:2517-2542`).
    let mut ff = true;
    for stripe in ls {
        // `ChFi3d_IndexOfSurfData(V, itel.Value(), Sens);` (`cxx:2519`).
        let mut sens = 0;
        ch_fi3d_index_of_surf_data(v, stripe, &mut sens)?;
        // `if (!FF && Stripe == itel.Value()) Sens = -Sens;` (`cxx:2520-2523`).
        if !ff {
            if let Some(prev) = &stripe_prev {
                if Arc::ptr_eq(prev, stripe) {
                    sens = -sens;
                }
            }
        }
        // `Stripe = itel.Value(); Spine = Stripe->Spine();` (`cxx:2524-2525`).
        stripe_prev = Some(stripe.clone());
        let spine = spine_of(stripe)?;
        let is_first = sens == 1;
        let mut guard = spine
            .write()
            .map_err(|_| "ExtentTwoCorner: poisoned spine".to_string())?;
        // `if (!Spine->IsTangencyExtremity((Sens == 1)))` (`cxx:2526`).
        if !guard.is_tangency_extremity(is_first) {
            if sens == 1 {
                // `Spine->SetFirstParameter(-Eval); Spine->SetFirstTgt(0.);`
                // (`cxx:2529-2532`).
                guard.set_first_parameter(-eval);
                guard
                    .set_first_tgt(0.0)
                    .ok_or_else(|| "ExtentTwoCorner: SetFirstTgt on a periodic spine".to_string())?;
            } else {
                // `dU = Spine->LastParameter(Spine->NbEdges());`
                // `Spine->SetLastParameter(dU + Eval); Spine->SetLastTgt(dU);`
                // (`cxx:2533-2538`).
                let d_u = guard
                    .last_parameter_of(guard.nb_edges() as i32)
                    .ok_or_else(|| {
                        "ExtentTwoCorner: ChFiDS_Spine::LastParameter(NbEdges)".to_string()
                    })?;
                guard.set_last_parameter(d_u + eval);
                guard.set_last_tgt(d_u).ok_or_else(|| {
                    "ExtentTwoCorner: SetLastTgt on a periodic spine".to_string()
                })?;
            }
            // `FF = false;` (`cxx:2539`).
            ff = false;
        }
    }
    Ok(())
}

/// `ExtentThreeCorner(V, LS)` (`cxx:2548-2594`).
///
/// Each guideline is extended by a tenth of its length unless the spine meets
/// the vertex on the same face of a convex edge (`cxx:2573-2587`); the walk
/// stops at the first tangency extremity (`cxx:2569-2572`).
pub(crate) fn extent_three_corner(v: &TopoShape, ls: &[StripeHandle]) -> Result<(), String> {
    // `int Sens = 0; double Coeff = 0.1;` (`cxx:2552-2553`).
    let mut sens = 0;
    let coeff = THREE_COEFF;
    // `NCollection_List<occ::handle<ChFiDS_Stripe>> check;` (`cxx:2554`).
    let mut check: Vec<StripeHandle> = Vec::new();
    for stripe in ls {
        // `ChFi3d_IndexOfSurfData(V, Stripe, Sens);` (`cxx:2558`).
        ch_fi3d_index_of_surf_data(v, stripe, &mut sens)?;
        // `for (ich(check); ich.More(); ich.Next()) { if (Stripe ==
        // ich.Value()) { Sens = -Sens; break; } }` (`cxx:2559-2566`).
        for seen in &check {
            if Arc::ptr_eq(seen, stripe) {
                sens = -sens;
                break;
            }
        }
        // `occ::handle<ChFiDS_Spine> Spine = Stripe->Spine();` (`cxx:2567`).
        let spine = spine_of(stripe)?;
        let is_first = sens == 1;
        let mut guard = spine
            .write()
            .map_err(|_| "ExtentThreeCorner: poisoned spine".to_string())?;
        // `if (Spine->IsTangencyExtremity((Sens == 1))) return;`
        // (`cxx:2568-2571`).
        if guard.is_tangency_extremity(is_first) {
            return Ok(());
        }
        // `double dU = Spine->LastParameter(Spine->NbEdges());` (`cxx:2572`).
        let d_u = guard
            .last_parameter_of(guard.nb_edges() as i32)
            .ok_or_else(|| {
                "ExtentThreeCorner: ChFiDS_Spine::LastParameter(NbEdges)".to_string()
            })?;
        if sens == 1 {
            // `if (GetTypeOfConcavity() != ChFiDS_Convex || FirstStatus() !=
            // ChFiDS_OnSame)` (`cxx:2575-2579`).
            if guard.type_of_concavity() != TypeOfConcavity::Convex
                || guard.first_status() != SpineState::OnSame
            {
                guard.set_first_parameter(-d_u * coeff);
                guard.set_first_tgt(0.0).ok_or_else(|| {
                    "ExtentThreeCorner: SetFirstTgt on a periodic spine".to_string()
                })?;
            }
        } else {
            // `if (GetTypeOfConcavity() != ChFiDS_Convex || LastStatus() !=
            // ChFiDS_OnSame)` (`cxx:2582-2586`).
            if guard.type_of_concavity() != TypeOfConcavity::Convex
                || guard.last_status() != SpineState::OnSame
            {
                guard.set_last_parameter(d_u * (1.0 + coeff));
                guard.set_last_tgt(d_u).ok_or_else(|| {
                    "ExtentThreeCorner: SetLastTgt on a periodic spine".to_string()
                })?;
            }
        }
        // `check.Append(Stripe);` (`cxx:2589`).
        check.push(stripe.clone());
    }
    Ok(())
}
