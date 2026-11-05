//! `ChFi3d_ChBuilder::ExtentOneCorner` / `ExtentTwoCorner` /
//! `ExtentThreeCorner` — the pre-analysis extensions of the guidelines at a
//! vertex, chamfer flavour.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_ChBuilder.cxx:1953-2212`
//! (OCCT 8.0.0). These are **not** the fillet bodies
//! (`ChFi3d_FilBuilder.cxx:2440-2591`, ported in
//! [`crate::chfi3d::builder::fil_builder::extent_corner`]): the chamfer
//! `ExtentOneCorner` has no `ChFiDS_FreeBoundary` doubling (`cxx:1953-1976`),
//! and its two- / three-corner bodies extend by the *chamfer distances* read
//! off the two common faces (`ConexFaces` + `IsChamfer` / `Dist` / `Dists` /
//! `GetDistAngle`) instead of by a `MaxRad` multiple.
//!
//! Translation notes:
//! - `Spine->LastParameter(Spine->NbEdges())` answers `None` when the abscissa
//!   is not built yet (`ChFiDS_Spine.cxx:248-251`); OCCT would return an
//!   undefined `double`, the port reports `Err`,
//! - `ChFiDS_ChamfSpine::GetDist` / `Dists` / `GetDistAngle` raise
//!   `Standard_Failure` on a spine of another kind (`ChFiDS_ChamfSpine.cxx:47-97`);
//!   the port answers `None` and these bodies report `Err` — every call here is
//!   guarded by `IsChamfer() == <kind>`, exactly like OCCT,
//! - `occ::handle` comparison `Stripe == itel.Value()` is OCCT handle equality,
//!   i.e. [`Arc::ptr_eq`] on the port's handles,
//! - `ChFi3d_IndexOfSurfData` raises `Standard_ConstructionError` for a vertex
//!   on neither end of the spine; the port propagates its `Err`,
//! - `double dis[3][3]` (`cxx:2180`) is read by the last loop even when no face
//!   pair matched, i.e. OCCT reads uninitialised doubles there; the port starts
//!   it at `0.` (the same value `ExtentTwoCorner` initialises its own `dis` to)
//!   and does not add a guard.

use std::sync::Arc;

use crate::chfi3d::builder_0_helpers::ch_fi3d_index_of_surf_data;
use crate::chfi_ds::{ChamfMethod, ChFiDSMap, SpineClass, SpineHandle, State as SpineState, StripeHandle};
use crate::shape::{Face, TopoShape};
use crate::topo_tools::is_same;

use super::statics::{cone_x_faces, extent_spine_on_common_face};

/// `Coeff = 0.5` of `ExtentOneCorner` (`cxx:1957`).
const ONE_COEFF: f64 = 0.5;

/// `S->Spine()` of the three bodies.
fn spine_of(stripe: &StripeHandle) -> Result<SpineHandle, String> {
    let guard = stripe
        .read()
        .map_err(|_| "ChExtentCorner: poisoned stripe".to_string())?;
    guard
        .spine()
        .cloned()
        .ok_or_else(|| "ChExtentCorner: the stripe has no spine".to_string())
}

/// The two distances of one stripe: `cxx:2028-2043` and `cxx:2158-2175`.
///
/// `ChFiDS_Sym` answers `GetDist` twice, `ChFiDS_TwoDist` answers `Dists`, and
/// any other kind approximates the second distance with `tmd * tan(angle)`
/// (`cxx:2041`, `:2172`).
fn chamfer_dists(spine: &SpineHandle) -> Result<(f64, f64), String> {
    let guard = spine
        .read()
        .map_err(|_| "ChExtentCorner: poisoned spine".to_string())?;
    if guard.class() != SpineClass::Chamf {
        return Err("ChExtentCorner: the spine is not a ChFiDS_ChamfSpine".to_string());
    }
    match guard.is_chamfer() {
        // `chsp[i]->GetDist(d[j]); d[j + 1] = d[j];` (`cxx:2030-2032`).
        ChamfMethod::Sym => {
            let d = guard
                .get_dist()
                .ok_or_else(|| "ChExtentCorner: GetDist on a two-distance chamfer".to_string())?;
            Ok((d, d))
        }
        // `chsp[i]->Dists(d[j], d[j + 1]);` (`cxx:2034-2037`).
        ChamfMethod::TwoDist => guard
            .dists()
            .ok_or_else(|| "ChExtentCorner: Dists on a symmetric chamfer".to_string()),
        // `chsp[i]->GetDistAngle(tmd, tmpang); d[j] = tmd; d[j + 1] =
        // tmd * tan(tmpang);` (`cxx:2039-2043`).
        ChamfMethod::DistAngle => {
            let (tmd, tmpang) = guard
                .get_dist_angle()
                .ok_or_else(|| "ChExtentCorner: GetDistAngle on another kind".to_string())?;
            Ok((tmd, tmd * tmpang.tan()))
        }
    }
}

/// `ExtentOneCorner(V, S)` (`cxx:1953-1976`): the guideline of `S` is extended
/// by half its length towards `V`.
pub(super) fn extent_one_corner(v: &TopoShape, s: &StripeHandle) -> Result<(), String> {
    // `int Sens = 0; double Coeff = 0.5;` (`cxx:1955-1956`).
    let mut sens = 0;
    let coeff = ONE_COEFF;
    // `occ::handle<ChFiDS_Spine> Spine = S->Spine();`
    // `ChFi3d_IndexOfSurfData(V, S, Sens);` (`cxx:1958-1959`).
    let spine = spine_of(s)?;
    ch_fi3d_index_of_surf_data(v, s, &mut sens)?;
    let is_first = sens == 1;

    let mut guard = spine
        .write()
        .map_err(|_| "ChExtentOneCorner: poisoned spine".to_string())?;
    // `if (Spine->IsTangencyExtremity((Sens == 1))) return;` (`cxx:1960-1963`).
    if guard.is_tangency_extremity(is_first) {
        return Ok(());
    }
    // `double dU = Spine->LastParameter(Spine->NbEdges());` (`cxx:1964`).
    let d_u = guard
        .last_parameter_of(guard.nb_edges() as i32)
        .ok_or_else(|| "ChExtentOneCorner: ChFiDS_Spine::LastParameter(NbEdges)".to_string())?;
    if is_first {
        // `Spine->SetFirstParameter(-dU * Coeff); Spine->SetFirstTgt(0.);`
        // (`cxx:1966-1970`).
        guard.set_first_parameter(-d_u * coeff);
        guard
            .set_first_tgt(0.0)
            .ok_or_else(|| "ChExtentOneCorner: SetFirstTgt on a periodic spine".to_string())?;
    } else {
        // `Spine->SetLastParameter(dU * (1. + Coeff)); Spine->SetLastTgt(dU);`
        // (`cxx:1971-1975`).
        guard.set_last_parameter(d_u * (1.0 + coeff));
        guard
            .set_last_tgt(d_u)
            .ok_or_else(|| "ChExtentOneCorner: SetLastTgt on a periodic spine".to_string())?;
    }
    Ok(())
}

/// `ExtentTwoCorner(V, LS)` (`cxx:1982-2107`).
///
/// The distances are read off the two common faces of the two guidelines
/// (`cxx:2022-2045`), paired by `IsSame` (`cxx:2048-2067`), and then applied:
/// two `AllSame` end statuses extend both stripes with
/// [`extent_one_corner`], two `OnSame` ones with
/// [`extent_spine_on_common_face`] (`cxx:2070-2106`).
pub(super) fn extent_two_corner(
    ef_map: &ChFiDSMap,
    v: &TopoShape,
    ls: &[StripeHandle],
) -> Result<(), String> {
    // `bool FF = true; bool isfirst[2]; int Iedge[2] = {1, 1}; Stripe[2];
    // Spine[2];` (`cxx:1987-1995`).
    let mut ff = true;
    let mut isfirst = [false; 2];
    let mut iedge = [1usize, 1];
    let mut stripes: [Option<StripeHandle>; 2] = [None, None];
    let mut spines: [Option<SpineHandle>; 2] = [None, None];

    // `for (; itel.More(); itel.Next(), i++) { ... }` (`cxx:1997-2015`).
    for (i, stripe) in ls.iter().take(2).enumerate() {
        // `ChFi3d_IndexOfSurfData(V, itel.Value(), Sens);` (`cxx:1999`).
        let mut sens = 0;
        ch_fi3d_index_of_surf_data(v, stripe, &mut sens)?;
        // `if (!FF) { if (Stripe[1] == itel.Value()) Sens = -Sens; }`
        // (`cxx:2000-2006`). `Stripe[1]` is assigned *after* this test, so on
        // the second iteration it still holds the `Stripe[2]` initial null
        // handle and the negation never fires for a two-element list. The port
        // keeps that statement order.
        if !ff {
            if let Some(second) = &stripes[1] {
                if Arc::ptr_eq(second, stripe) {
                    sens = -sens;
                }
            }
        }
        // `Stripe[i] = itel.Value(); isfirst[i] = (Sens == 1);
        // Spine[i] = Stripe[i]->Spine(); if (!isfirst[i]) Iedge[i] =
        // Spine[i]->NbEdges(); FF = false;` (`cxx:2008-2014`).
        stripes[i] = Some(stripe.clone());
        isfirst[i] = sens == 1;
        let spine = spine_of(stripe)?;
        if !isfirst[i] {
            iedge[i] = spine
                .read()
                .map_err(|_| "ChExtentTwoCorner: poisoned spine".to_string())?
                .nb_edges();
        }
        spines[i] = Some(spine);
        ff = false;
    }

    // `double d[4], dis[2] = {0.0, 0.0}; TopoDS_Face F[4];` (`cxx:2017-2021`).
    let mut d = [0.0f64; 4];
    let mut dis = [0.0f64, 0.0];
    let mut faces: [Option<Face>; 4] = [None, None, None, None];

    // `for (i = 0, j = 0; i < 2; i++, j += 2) { ... }` (`cxx:2023-2045`).
    for i in 0..2 {
        let j = 2 * i;
        let Some(spine) = spines[i].clone() else {
            return Err("ChExtentTwoCorner: the stripe list has fewer than 2 stripes".to_string());
        };
        // `ConexFaces(Spine[i], Iedge[i], F[j], F[j + 1]);` (`cxx:2027`).
        let (f1, f2) = cone_x_faces(ef_map, &spine, iedge[i])?;
        faces[j] = Some(f1);
        faces[j + 1] = Some(f2);
        // The `IsChamfer()` switch (`cxx:2028-2043`).
        let (d0, d1) = chamfer_dists(&spine)?;
        d[j] = d0;
        d[j + 1] = d1;
    }

    // `bool notfound = true; i = 0; while (notfound && (i < 2)) { j = 0;
    // while (notfound && (j < 2)) { if (F[i].IsSame(F[j + 2])) { dis[0] =
    // d[i]; dis[1] = d[j + 2]; notfound = false; } j++; } i++; }`
    // (`cxx:2047-2068`).
    let mut notfound = true;
    for i in 0..2 {
        if !notfound {
            break;
        }
        for j in 0..2 {
            if !notfound {
                break;
            }
            let (Some(fi), Some(fj)) = (&faces[i], &faces[j + 2]) else {
                continue;
            };
            if is_same(&fi.0, &fj.0) {
                dis[0] = d[i];
                dis[1] = d[j + 2];
                notfound = false;
            }
        }
    }

    // `for (i = 0; i < 2; i++) State[i] = isfirst[i] ? Spine[i]->FirstStatus()
    // : Spine[i]->LastStatus();` (`cxx:2070-2081`).
    let mut state = [SpineState::OnSame; 2];
    for i in 0..2 {
        let Some(spine) = spines[i].clone() else {
            return Err("ChExtentTwoCorner: the stripe list has fewer than 2 stripes".to_string());
        };
        let guard = spine
            .read()
            .map_err(|_| "ChExtentTwoCorner: poisoned spine".to_string())?;
        state[i] = if isfirst[i] {
            guard.first_status()
        } else {
            guard.last_status()
        };
    }

    // `if (State[0] == ChFiDS_AllSame) { for (j = 0; j < 2; j++)
    // ExtentOneCorner(V, Stripe[j]); }` (`cxx:2082-2102`).
    if state[0] == SpineState::AllSame {
        for stripe in stripes.iter().flatten() {
            extent_one_corner(v, stripe)?;
        }
    }
    // `else if ((State[0] == ChFiDS_OnSame) && (State[1] == ChFiDS_OnSame))
    // ExtentSpineOnCommonFace(Spine[0], Spine[1], V, dis[0], dis[1],
    // isfirst[0], isfirst[1]);` (`cxx:2103-2106`).
    else if state[0] == SpineState::OnSame && state[1] == SpineState::OnSame {
        let (Some(spine0), Some(spine1)) = (spines[0].clone(), spines[1].clone()) else {
            return Err("ChExtentTwoCorner: the stripe list has fewer than 2 stripes".to_string());
        };
        extent_spine_on_common_face(
            &spine0, &spine1, v, dis[0], dis[1], isfirst[0], isfirst[1],
        )?;
    }
    Ok(())
}

/// `ExtentThreeCorner(V, LS)` (`cxx:2111-2212`).
///
/// Every guideline is paired with the next one (`j = (i + 1) % 3`) and both are
/// extended by the distance the pair shares on its common face
/// (`cxx:2182-2211`).
pub(super) fn extent_three_corner(
    ef_map: &ChFiDSMap,
    v: &TopoShape,
    ls: &[StripeHandle],
) -> Result<(), String> {
    // `NCollection_List<occ::handle<ChFiDS_Stripe>> check; bool isfirst[3];
    // int Iedge[3] = {1, 1, 1}; Spine[3];` (`cxx:2115-2122`).
    let mut check: Vec<StripeHandle> = Vec::new();
    let mut isfirst = [false; 3];
    let mut iedge = [1usize, 1, 1];
    let mut spines: [Option<SpineHandle>; 3] = [None, None, None];

    // `for (itel(LS); itel.More(); itel.Next(), i++) { ... }` (`cxx:2124-2147`).
    for (i, stripe) in ls.iter().take(3).enumerate() {
        // `ChFi3d_IndexOfSurfData(V, Stripe, Sens);` (`cxx:2128`).
        let mut sens = 0;
        ch_fi3d_index_of_surf_data(v, stripe, &mut sens)?;
        // `for (ich(check); ich.More(); ich.Next()) { if (Stripe == ich.Value())
        // { Sens = -Sens; break; } }` (`cxx:2129-2136`).
        for seen in &check {
            if Arc::ptr_eq(seen, stripe) {
                sens = -sens;
                break;
            }
        }
        // `isfirst[i] = (Sens == 1); Spine[i] = Stripe->Spine(); if
        // (!isfirst[i]) Iedge[i] = Spine[i]->NbEdges(); check.Append(Stripe);`
        // (`cxx:2138-2146`).
        isfirst[i] = sens == 1;
        let spine = spine_of(stripe)?;
        if !isfirst[i] {
            iedge[i] = spine
                .read()
                .map_err(|_| "ChExtentThreeCorner: poisoned spine".to_string())?
                .nb_edges();
        }
        spines[i] = Some(spine);
        check.push(stripe.clone());
    }

    // `double d[3][2], tmd, tmpangle; TopoDS_Face F[3][2]; chsp[3];`
    // (`cxx:2149-2152`).
    let mut d = [[0.0f64; 2]; 3];
    let mut faces: [[Option<Face>; 2]; 3] = [[None, None], [None, None], [None, None]];

    // `for (i = 0; i < 3; i++) { ... }` (`cxx:2154-2177`).
    for i in 0..3 {
        let Some(spine) = spines[i].clone() else {
            return Err(
                "ChExtentThreeCorner: the stripe list has fewer than 3 stripes".to_string(),
            );
        };
        // `ConexFaces(Spine[i], Iedge[i], F[i][0], F[i][1]);` (`cxx:2158`).
        let (f1, f2) = cone_x_faces(ef_map, &spine, iedge[i])?;
        faces[i][0] = Some(f1);
        faces[i][1] = Some(f2);
        // The `IsChamfer()` switch (`cxx:2159-2176`).
        let (d0, d1) = chamfer_dists(&spine)?;
        d[i][0] = d0;
        d[i][1] = d1;
    }

    // `double dis[3][3]; for (i = 0; i < 3; i++) { j = (i + 1) % 3;
    // bool notfound = true; k = 0; while (notfound && (k < 2)) { l = 0;
    // while (notfound && (l < 2)) { if (F[i][k].IsSame(F[j][l])) {
    // dis[i][j] = d[i][k]; dis[j][i] = d[j][l]; notfound = false; } l++; }
    // k++; } }` (`cxx:2180-2204`).
    let mut dis = [[0.0f64; 3]; 3];
    for i in 0..3 {
        let j = (i + 1) % 3;
        let mut notfound = true;
        let mut k = 0;
        while notfound && k < 2 {
            let mut l = 0;
            while notfound && l < 2 {
                let (Some(fik), Some(fjl)) = (&faces[i][k], &faces[j][l]) else {
                    l += 1;
                    continue;
                };
                if is_same(&fik.0, &fjl.0) {
                    dis[i][j] = d[i][k];
                    dis[j][i] = d[j][l];
                    notfound = false;
                }
                l += 1;
            }
            k += 1;
        }
    }

    // `for (i = 0; i < 3; i++) { j = (i + 1) % 3;
    // ExtentSpineOnCommonFace(Spine[i], Spine[j], V, dis[i][j], dis[j][i],
    // isfirst[i], isfirst[j]); }` (`cxx:2207-2211`).
    for i in 0..3 {
        let j = (i + 1) % 3;
        let (Some(spine_i), Some(spine_j)) = (spines[i].clone(), spines[j].clone()) else {
            return Err(
                "ChExtentThreeCorner: the stripe list has fewer than 3 stripes".to_string(),
            );
        };
        extent_spine_on_common_face(
            &spine_i,
            &spine_j,
            v,
            dis[i][j],
            dis[j][i],
            isfirst[i],
            isfirst[j],
        )?;
    }
    Ok(())
}
