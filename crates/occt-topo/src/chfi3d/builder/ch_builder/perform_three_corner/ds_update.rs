//! `ChFi3d_ChBuilder::PerformThreeCorner` - the DS update of the corner, and
//! the truncation of the three stripe sequences.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_ChBuilder_C3.cxx:888-1142`
//! (OCCT 8.0.0):
//!
//! - `cxx:889-1098` - the four stripes and the DS: [`write_ds`],
//! - `cxx:1101-1142` - the three `Remove` calls and the two appends:
//!   [`trim_and_append`].
//!
//! Differences from the `ChFi3d_FilBuilder` twin of the same block
//! (`ChFi3d_FilBuilder_C3.cxx:1016-1331`), all read off the `.cxx`:
//! - no `ChFi3d_EnlargeBox` and no `ChFi3d_SetPointTolerance`: the chamfer
//!   corner knows its p-curves exactly (`ComputeIntersection`), so the bounding
//!   boxes of the `SurfData` are not grown here,
//! - no `SameParameter` pass on the pivot p-curve (`cxx:1088-1089` assigns
//!   `fi.PCurveOnFace()` as it stands),
//! - the p-curves of the two stripe `CornerData` come from `debpc1` / `finpc1`
//!   in the `c1plan` case and from `ChFi3d_ComputePCurv(UV1, UV2, ..., rev)` -
//!   the 2D-only overload - otherwise,
//! - `Pl1` is a *reference* to `coin->ChangeVertexLastOnS1()` that is
//!   reassigned from `VertexFirstOnS1()` when `c1triangle` (`cxx:899-902`); the
//!   port copies the first common point into the last one and then reads it,
//!   which is what the C++ reference sees.
//!
//! Translation notes:
//! - `ChFiDS_SurfData& fddeb` / `fdfin` / `fdpiv` are reached through [`sd`]
//!   and [`with_sd`], which lock the stripe and take `Arc::make_mut` on the
//!   slot, so the aliasing OCCT gets from its handles (`fddeb == fdpiv` when
//!   the two stripes are one) is preserved,
//! - `CD[pivot]->InDS(isfirst)` is `in_ds(isfirst, 1)`, the default `Nb`.

use std::sync::Arc;

use occt_geom2d::curve::Curve2d;

use crate::abs::Orientation;
use crate::chfi3d::builder_0_helpers::{
    ch_fi3d_compute_arete, ch_fi3d_compute_pcurv_uv, ch_fi3d_index_point_in_ds,
    ComputeArete,
};
use crate::chfi_ds::{Regul, Stripe, SurfData};
use crate::topopebrep::ds::Curve as DsCurve;

use super::super::super::corner_sd::{data_surface, ds_curve, pcurve_on_surf, sd, with_sd};
use super::super::super::ChFi3dBuilder;
use super::{PlanCorner, ThreeChamfer};

/// `cxx:889-1098`, run when `done`.
pub(super) fn write_ds(
    base: &mut ChFi3dBuilder,
    c: &mut ThreeChamfer,
    corner: &mut Stripe,
    coin: &mut SurfData,
    plan: Option<&PlanCorner>,
) -> Result<(), String> {
    // `const ChFiDS_CommonPoint& Pf1 = coin->VertexFirstOnS1();
    // const ChFiDS_CommonPoint& Pf2 = coin->VertexFirstOnS2();
    // ChFiDS_CommonPoint& Pl1 = coin->ChangeVertexLastOnS1();
    // if (c1triangle) { Pl1 = coin->ChangeVertexFirstOnS1(); }
    // const ChFiDS_CommonPoint& Pl2 = coin->VertexLastOnS2();`
    // (`cxx:898-903`).
    let pf1 = coin.vertex(true, 1).clone();
    let pf2 = coin.vertex(true, 2).clone();
    if c.c1triangle {
        let first_on_s1 = coin.vertex(true, 1).clone();
        coin.vertex_mut(false, 1).clone_from(&first_on_s1);
    }
    let pl1 = coin.vertex(false, 1).clone();
    let pl2 = coin.vertex(false, 2).clone();

    // `ChFiDS_Regul regdeb, regfin;` (`cxx:907`).
    let mut regdeb = Regul::new();
    let mut regfin = Regul::new();

    // `If1 = ChFi3d_IndexPointInDS(Pf1, DStr); If2 = ...; if (c1triangle) {
    // Il1 = If1; } else { Il1 = ChFi3d_IndexPointInDS(Pl1, DStr); } Il2 = ...;`
    // (`cxx:908-918`).
    let if1 = ch_fi3d_index_point_in_ds(&pf1, base.my_ds.change_ds());
    let if2 = ch_fi3d_index_point_in_ds(&pf2, base.my_ds.change_ds());
    let il1 = if c.c1triangle {
        if1
    } else {
        ch_fi3d_index_point_in_ds(&pl1, base.my_ds.change_ds())
    };
    let il2 = ch_fi3d_index_point_in_ds(&pl2, base.my_ds.change_ds());

    // `coin->ChangeIndexOfS1(DStr.AddShape(face[pivot]));
    // coin->ChangeIndexOfS2(-fdpiv->Surf());` (`cxx:921-922`).
    let pivot_face = c.face[c.pivot]
        .clone()
        .ok_or_else(|| "PerformThreeCorner: the pivot has no face".to_string())?;
    let index_of_s1 = base.my_ds.change_ds().add_shape(&pivot_face.0);
    coin.set_index_of_s1(index_of_s1);
    let fdpiv_surf = sd(&c.stripes[c.pivot], c.i[c.pivot][c.deb])?.surf();
    coin.set_index_of_s2(-fdpiv_surf);

    // `DStr.Surface(coin->Surf()).Surface()` (`cxx:943`).
    let surf_coin = data_surface(base.my_ds.ds(), coin)?;

    // `// first points` (`cxx:924-956`).
    let icf: i32;
    let (p1deb, p2deb): (f64, f64);
    match plan {
        Some(plan) => {
            // `P1deb = DStr.Curve(Icf).Curve()->FirstParameter(); P2deb =
            // DStr.Curve(Icf).Curve()->LastParameter();` (`cxx:927-930`).
            let curve = ds_curve(base, plan.icf)?;
            p1deb = curve.first_parameter();
            p2deb = curve.last_parameter();
            icf = plan.icf;
        }
        None => {
            // `pp1 = coin->InterferenceOnS1().PCurveOnSurf()->Value(
            // coin->InterferenceOnS1().FirstParameter()); pp2 = ... S2 ...`
            // (`cxx:932-935`).
            let pp1 = pcurve_on_surf(coin, 1, coin.interference_on_s1().first_parameter())?;
            let pp2 = pcurve_on_surf(coin, 2, coin.interference_on_s2().first_parameter())?;
            // `ChFi3d_ComputeArete(Pf1, pp1, Pf2, pp2, DStr.Surface(coin->Surf())
            // .Surface(), C3d, corner->ChangeFirstPCurve(), P1deb, P2deb,
            // tolapp3d, tol2d, tolreached, 0);` (`cxx:937-952`).
            let mut first = ComputeArete::default();
            ch_fi3d_compute_arete(
                &pf1,
                &pp1,
                &pf2,
                &pp2,
                &surf_coin,
                &mut first,
                base.tolapp3d,
                base.tol2d,
                0,
            )?;
            p1deb = first.pardeb;
            p2deb = first.parfin;
            *corner.change_first_p_curve() = first.p_curve.clone();
            // `TopOpeBRepDS_Curve Tcurv(C3d, tolreached); Icf =
            // DStr.AddCurve(Tcurv);` (`cxx:953-955`).
            let c3d = first
                .c3d
                .clone()
                .ok_or_else(|| "PerformThreeCorner: ChFi3d_ComputeArete gave no C3d".to_string())?;
            let mut tcurv = DsCurve::with_curve(Some(c3d), first.tolreached, false);
            icf = base.my_ds.change_ds().add_curve(&mut tcurv);
        }
    }

    // `regdeb.SetCurve(Icf); regdeb.SetS1(coin->Surf(), false);
    // regdeb.SetS2(fddeb->Surf(), false); myRegul.Append(regdeb);`
    // (`cxx:957-960`).
    regdeb.set_curve(icf);
    regdeb.set_s1(coin.surf(), false);
    let fddeb_surf = sd(&c.stripes[c.deb], c.i[c.deb][c.pivot])?.surf();
    regdeb.set_s2(fddeb_surf, false);
    base.my_regul.push(regdeb);
    // `corner->ChangeFirstCurve(Icf); corner->ChangeFirstParameters(P1deb,
    // P2deb); corner->ChangeIndexFirstPointOnS1(If1);
    // corner->ChangeIndexFirstPointOnS2(If2);` (`cxx:961-964`).
    corner.change_first_curve(icf);
    corner.change_first_parameters(p1deb, p2deb);
    corner.change_index_first_point_on_s1(if1);
    corner.change_index_first_point_on_s2(if2);

    // `// last points` (`cxx:966-1002`).
    let icl: i32;
    let (p1fin, p2fin): (f64, f64);
    match plan {
        Some(plan) => {
            // `P1fin = DStr.Curve(Icl).Curve()->FirstParameter(); P2fin =
            // DStr.Curve(Icl).Curve()->LastParameter();` (`cxx:968-972`).
            let curve = ds_curve(base, plan.icl)?;
            p1fin = curve.first_parameter();
            p2fin = curve.last_parameter();
            icl = plan.icl;
        }
        None => {
            // `pp1 = coin->InterferenceOnS1().PCurveOnSurf()->Value(
            // coin->InterferenceOnS1().LastParameter()); pp2 = ... S2 ...`
            // (`cxx:974-977`).
            let pp1 = pcurve_on_surf(coin, 1, coin.interference_on_s1().last_parameter())?;
            let pp2 = pcurve_on_surf(coin, 2, coin.interference_on_s2().last_parameter())?;
            // `ChFi3d_ComputeArete(Pl1, pp1, Pl2, pp2, ..., C3d,
            // corner->ChangeLastPCurve(), P1fin, P2fin, tolapp3d, tol2d,
            // tolreached, 0);` (`cxx:979-993`).
            let mut last = ComputeArete::default();
            ch_fi3d_compute_arete(
                &pl1,
                &pp1,
                &pl2,
                &pp2,
                &surf_coin,
                &mut last,
                base.tolapp3d,
                base.tol2d,
                0,
            )?;
            p1fin = last.pardeb;
            p2fin = last.parfin;
            *corner.change_last_p_curve() = last.p_curve.clone();
            let c3d = last
                .c3d
                .clone()
                .ok_or_else(|| "PerformThreeCorner: ChFi3d_ComputeArete gave no C3d".to_string())?;
            let mut tcurv = DsCurve::with_curve(Some(c3d), last.tolreached, false);
            icl = base.my_ds.change_ds().add_curve(&mut tcurv);
        }
    }

    // `regfin.SetCurve(Icl); regfin.SetS1(coin->Surf(), false);
    // regfin.SetS2(fdfin->Surf(), false); myRegul.Append(regfin);`
    // (`cxx:995-998`).
    regfin.set_curve(icl);
    regfin.set_s1(coin.surf(), false);
    let fdfin_surf = sd(&c.stripes[c.fin], c.i[c.fin][c.pivot])?.surf();
    regfin.set_s2(fdfin_surf, false);
    base.my_regul.push(regfin);
    // `corner->ChangeLastCurve(Icl); corner->ChangeLastParameters(P1fin,
    // P2fin); corner->ChangeIndexLastPointOnS1(Il1);
    // corner->ChangeIndexLastPointOnS2(Il2);` (`cxx:999-1002`).
    corner.change_last_curve(icl);
    corner.change_last_parameters(p1fin, p2fin);
    corner.change_index_last_point_on_s1(il1);
    corner.change_index_last_point_on_s2(il2);

    // ---------------------------------------------------------------
    // `// puis la CornerData du debut,` (`cxx:1005-1036`).
    // ---------------------------------------------------------------
    // `bool isfirst = (sens[deb] == 1), rev = (jf[deb][fin] == 2); int isurf1
    // = 1, isurf2 = 2; double par = p[deb][pivot], par2 = p[deb][pivot]; if
    // (c1triangle) { par2 = p[deb][fin]; }` (`cxx:1008-1012`).
    let is_first = c.sens[c.deb] == 1;
    let rev = c.jf[c.deb][c.fin] == 2;
    let (mut isurf1, mut isurf2) = (1, 2);
    let par = c.p[c.deb][c.pivot];
    let par2 = if c.c1triangle {
        c.p[c.deb][c.fin]
    } else {
        c.p[c.deb][c.pivot]
    };
    if rev {
        isurf1 = 2;
        isurf2 = 1;
        c.stripes[c.deb]
            .write()
            .map_err(|_| "PerformThreeCorner: poisoned stripe".to_string())?
            .set_orientation(Orientation::Reversed, is_first);
    }
    // `CD[deb]->SetCurve(Icf, isfirst); CD[deb]->SetIndexPoint(If1, isfirst,
    // isurf1); CD[deb]->SetIndexPoint(If2, isfirst, isurf2);
    // CD[deb]->SetParameters(isfirst, P1deb, P2deb);` (`cxx:1019-1022`).
    {
        let mut guard = c.stripes[c.deb]
            .write()
            .map_err(|_| "PerformThreeCorner: poisoned stripe".to_string())?;
        guard.set_curve(icf, is_first);
        guard.set_index_point(if1, is_first, isurf1);
        guard.set_index_point(if2, is_first, isurf2);
        guard.set_parameters(is_first, p1deb, p2deb);
    }
    // `fddeb->ChangeVertex(isfirst, isurf1) = Pf1; fddeb->ChangeVertex(isfirst,
    // isurf2) = Pf2; fddeb->ChangeInterference(isurf1).SetParameter(par2,
    // isfirst); fddeb->ChangeInterference(isurf2).SetParameter(par, isfirst);`
    // (`cxx:1023-1026`).
    with_sd(&c.stripes[c.deb], c.i[c.deb][c.pivot], |data| {
        data.vertex_mut(is_first, isurf1).clone_from(&pf1);
        data.vertex_mut(is_first, isurf2).clone_from(&pf2);
        data.interference_mut(isurf1).set_parameter(par2, is_first);
        data.interference_mut(isurf2).set_parameter(par, is_first);
    })?;
    // `if (c1plan) { CD[deb]->ChangePCurve(isfirst) = debpc1; } else { pp1 =
    // fddeb->InterferenceOnS1().PCurveOnSurf()->Value(par); pp2 = ...;
    // ChFi3d_ComputePCurv(pp1, pp2, CD[deb]->ChangePCurve(isfirst), P1deb,
    // P2deb, rev); }` (`cxx:1027-1035`).
    let pc_deb: Option<Arc<dyn Curve2d>> = match plan {
        Some(plan) => plan.debpc1.clone(),
        None => {
            let fddeb = sd(&c.stripes[c.deb], c.i[c.deb][c.pivot])?;
            let dp1 = pcurve_on_surf(&fddeb, 1, par)?;
            let dp2 = pcurve_on_surf(&fddeb, 2, par)?;
            Some(ch_fi3d_compute_pcurv_uv(&dp1, &dp2, p1deb, p2deb, rev)?)
        }
    };
    {
        let mut guard = c.stripes[c.deb]
            .write()
            .map_err(|_| "PerformThreeCorner: poisoned stripe".to_string())?;
        *guard.change_p_curve(is_first) = pc_deb;
    }

    // ---------------------------------------------------------------
    // `// puis la CornerData de la fin,` (`cxx:1038-1075`).
    // ---------------------------------------------------------------
    // `isfirst = (sens[fin] == 1); rev = (jf[fin][deb] == 2); isurf1 = 1;
    // isurf2 = 2; par = p[fin][pivot]; par2 = p[fin][pivot]; if (c1triangle) {
    // par2 = p[fin][deb]; }` (`cxx:1040-1049`).
    let is_first = c.sens[c.fin] == 1;
    let rev = c.jf[c.fin][c.deb] == 2;
    let (mut isurf1, mut isurf2) = (1, 2);
    let par = c.p[c.fin][c.pivot];
    let par2 = if c.c1triangle {
        c.p[c.fin][c.deb]
    } else {
        c.p[c.fin][c.pivot]
    };
    if rev {
        isurf1 = 2;
        isurf2 = 1;
        c.stripes[c.fin]
            .write()
            .map_err(|_| "PerformThreeCorner: poisoned stripe".to_string())?
            .set_orientation(Orientation::Reversed, is_first);
    }
    // `CD[fin]->SetCurve(Icl, isfirst); CD[fin]->SetIndexPoint(Il1, isfirst,
    // isurf1); CD[fin]->SetIndexPoint(Il2, isfirst, isurf2);
    // CD[fin]->SetParameters(isfirst, P1fin, P2fin);` (`cxx:1056-1059`).
    {
        let mut guard = c.stripes[c.fin]
            .write()
            .map_err(|_| "PerformThreeCorner: poisoned stripe".to_string())?;
        guard.set_curve(icl, is_first);
        guard.set_index_point(il1, is_first, isurf1);
        guard.set_index_point(il2, is_first, isurf2);
        guard.set_parameters(is_first, p1fin, p2fin);
    }
    // `fdfin->ChangeVertex(isfirst, isurf1) = Pl1; fdfin->ChangeVertex(isfirst,
    // isurf2) = Pl2; fdfin->ChangeInterference(isurf1).SetParameter(par2,
    // isfirst); fdfin->ChangeInterference(isurf2).SetParameter(par, isfirst);`
    // (`cxx:1060-1063`).
    with_sd(&c.stripes[c.fin], c.i[c.fin][c.pivot], |data| {
        data.vertex_mut(is_first, isurf1).clone_from(&pl1);
        data.vertex_mut(is_first, isurf2).clone_from(&pl2);
        data.interference_mut(isurf1).set_parameter(par2, is_first);
        data.interference_mut(isurf2).set_parameter(par, is_first);
    })?;
    // `if (c1plan) { CD[fin]->ChangePCurve(isfirst) = finpc1; } else { ... }`
    // (`cxx:1064-1073`).
    let pc_fin: Option<Arc<dyn Curve2d>> = match plan {
        Some(plan) => plan.finpc1.clone(),
        None => {
            let fdfin = sd(&c.stripes[c.fin], c.i[c.fin][c.pivot])?;
            let fp1 = pcurve_on_surf(&fdfin, 1, par)?;
            let fp2 = pcurve_on_surf(&fdfin, 2, par)?;
            Some(ch_fi3d_compute_pcurv_uv(&fp1, &fp2, p1fin, p2fin, rev)?)
        }
    };
    {
        let mut guard = c.stripes[c.fin]
            .write()
            .map_err(|_| "PerformThreeCorner: poisoned stripe".to_string())?;
        *guard.change_p_curve(is_first) = pc_fin;
    }

    // ---------------------------------------------------------------
    // `// et enfin le pivot.` (`cxx:1077-1098`).
    // ---------------------------------------------------------------
    // `isfirst = (sens[pivot] == 1); rev = (jf[pivot][deb] == 2); isurf1 = 1;
    // isurf2 = 2; if (rev) { isurf1 = 2; isurf2 = 1;
    // CD[pivot]->SetOrientation(TopAbs_REVERSED, isfirst); }` (`cxx:1079-1087`).
    let is_first = c.sens[c.pivot] == 1;
    let rev = c.jf[c.pivot][c.deb] == 2;
    let (mut isurf1, mut isurf2) = (1, 2);
    if rev {
        isurf1 = 2;
        isurf2 = 1;
        c.stripes[c.pivot]
            .write()
            .map_err(|_| "PerformThreeCorner: poisoned stripe".to_string())?
            .set_orientation(Orientation::Reversed, is_first);
    }
    // `CD[pivot]->SetCurve(fi.LineIndex(), isfirst); CD[pivot]->ChangePCurve(
    // isfirst) = fi.PCurveOnFace(); CD[pivot]->SetParameters(isfirst,
    // fi.FirstParameter(), fi.LastParameter());` (`cxx:1088-1092`).
    let fi_line = coin.interference_on_s2().line_index();
    let fi_curve = coin.interference_on_s2().p_curve_on_face().cloned();
    let fi_first = coin.interference_on_s2().first_parameter();
    let fi_last = coin.interference_on_s2().last_parameter();
    {
        let mut guard = c.stripes[c.pivot]
            .write()
            .map_err(|_| "PerformThreeCorner: poisoned stripe".to_string())?;
        guard.set_curve(fi_line, is_first);
        *guard.change_p_curve(is_first) = fi_curve;
        guard.set_index_point(if2, is_first, isurf1);
        guard.set_index_point(il2, is_first, isurf2);
        guard.set_parameters(is_first, fi_first, fi_last);
        // `CD[pivot]->InDS(isfirst);` (`cxx:1097`).
        guard.in_ds(is_first, 1);
    }
    // `fdpiv->ChangeVertex(isfirst, isurf1) = Pf2; fdpiv->ChangeVertex(isfirst,
    // isurf2) = Pl2; fdpiv->ChangeInterference(isurf1).SetParameter(
    // p[pivot][deb], isfirst); fdpiv->ChangeInterference(isurf2).SetParameter(
    // p[pivot][fin], isfirst);` (`cxx:1093-1096`).
    let piv_deb = c.p[c.pivot][c.deb];
    let piv_fin = c.p[c.pivot][c.fin];
    with_sd(&c.stripes[c.pivot], c.i[c.pivot][c.deb], |data| {
        data.vertex_mut(is_first, isurf1).clone_from(&pf2);
        data.vertex_mut(is_first, isurf2).clone_from(&pl2);
        data.interference_mut(isurf1).set_parameter(piv_deb, is_first);
        data.interference_mut(isurf2).set_parameter(piv_fin, is_first);
    })?;

    Ok(())
}

/// `cxx:1101-1142`: the three `Remove` calls and the two appends.
pub(super) fn trim_and_append(
    base: &mut ChFi3dBuilder,
    c: &mut ThreeChamfer,
    corner: &mut Stripe,
    coin: &SurfData,
) -> Result<(), String> {
    // `if (i[deb][pivot] < Index[deb]) { CD[deb]->...Remove(i[deb][pivot] + 1,
    // Index[deb]); Index[deb] = i[deb][pivot]; } else if (i[deb][pivot] >
    // Index[deb]) { ...Remove(Index[deb], i[deb][pivot] - 1); i[deb][pivot] =
    // Index[deb]; }` (`cxx:1103-1112`), then the same for `fin`
    // (`cxx:1113-1122`) and `pivot` (`cxx:1124-1133`).
    // `NCollection_HSequence::Remove(From, To)` is 1-based and inclusive.
    fn trim_one(stripe: &mut Stripe, index: &mut i32, i_of: &mut i32) -> Result<(), String> {
        let from: i32;
        let to: i32;
        if *i_of < *index {
            from = *i_of + 1;
            to = *index;
            *index = *i_of;
        } else if *i_of > *index {
            from = *index;
            to = *i_of - 1;
            *i_of = *index;
        } else {
            return Ok(());
        }
        let first = usize::try_from(from - 1)
            .map_err(|_| "PerformThreeCorner: negative Remove range".to_string())?;
        let last = usize::try_from(to)
            .map_err(|_| "PerformThreeCorner: negative Remove range".to_string())?;
        let sequence = stripe.change_set_of_surf_data();
        if last > sequence.len() || first > last {
            return Err("PerformThreeCorner: Remove out of the SurfData range".to_string());
        }
        sequence.drain(first..last);
        Ok(())
    }

    {
        let mut guard = c.stripes[c.deb]
            .write()
            .map_err(|_| "PerformThreeCorner: poisoned stripe".to_string())?;
        trim_one(&mut guard, &mut c.index[c.deb], &mut c.i[c.deb][c.pivot])?;
    }
    {
        let mut guard = c.stripes[c.fin]
            .write()
            .map_err(|_| "PerformThreeCorner: poisoned stripe".to_string())?;
        trim_one(&mut guard, &mut c.index[c.fin], &mut c.i[c.fin][c.pivot])?;
    }
    {
        let mut guard = c.stripes[c.pivot]
            .write()
            .map_err(|_| "PerformThreeCorner: poisoned stripe".to_string())?;
        trim_one(
            &mut guard,
            &mut c.index[c.pivot],
            &mut c.i[c.pivot][c.deb],
        )?;
    }

    // `if (!myEVIMap.IsBound(Vtx)) { NCollection_List<int> li;
    // myEVIMap.Bind(Vtx, li); } myEVIMap.ChangeFind(Vtx).Append(coin->Surf());
    // corner->SetSolidIndex(CD[pivot]->SolidIndex());
    // myListStripe.Append(corner);` (`cxx:1134-1141`).
    base.change_evi(&c.vtx).push(coin.surf());
    let solid_index = c.stripes[c.pivot]
        .read()
        .map_err(|_| "PerformThreeCorner: poisoned stripe".to_string())?
        .solid_index();
    corner.set_solid_index(solid_index);
    // OCCT appended `coin` to `corner->ChangeSetOfSurfData()` when the corner
    // was created (`cxx:445-449`); the port kept it local, so the sequence
    // catches up here, before the stripe is shared.
    corner.change_set_of_surf_data().push(Arc::new(coin.clone()));
    base.my_list_stripe
        .push(Arc::new(std::sync::RwLock::new(corner.clone())));

    Ok(())
}
