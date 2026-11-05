//! The DS update of the `SameSide` arm: `ChFi3d_FilBuilder_C2.cxx:684-859`.
//!
//! `cxx:690-705` re-initializes the lost common point on the pivot arc, then
//! `cxx:706-859` rebuilds the two 3D edges on the corner surface
//! (`ChFi3d_ComputeArete`), computes their p-curves on the two stripe surfaces
//! (`ChFi3d_ComputePCurv`), stores both in the DS with a `ChFiDS_Regul` each
//! and rewrites the stripe / `SurfData` slots the corner replaced.

use std::sync::Arc;

use occt_geom::Surface;
use occt_geom2d::curve::Curve2d;

use crate::abs::Orientation;
use crate::chfi_ds::{Regul, Stripe, SurfData};
use crate::chfi3d::builder_0_helpers::{
    ch_fi3d_compute_arete, ch_fi3d_compute_pcurv, ch_fi3d_index_point_in_ds, ComputeArete,
};
use crate::chfi3d::builder::ChFi3dBuilder;
use crate::topopebrep::ds::Curve as DsCurve;

use super::super::CornerCtx;

/// `DStr.Surface(Fd->Surf()).Surface()` for one of the two `SurfData`.
fn data_surface(
    ds: &crate::topopebrep::ds::DataStructure,
    data: &SurfData,
) -> Result<Arc<dyn Surface>, String> {
    ds.surface(data.surf())
        .surface()
        .cloned()
        .ok_or_else(|| "PerformTwoCorner: DStr.Surface(Surf()).Surface() is null".to_string())
}

/// `cxx:690-859`, run when the corner surface is in place.
pub(super) fn update_ds(
    base: &mut ChFi3dBuilder,
    ctx: &mut CornerCtx,
    corner: &mut Stripe,
    coin: &mut SurfData,
) -> Result<(), String> {
    // `if (resetcp1) { gp_Pnt pjyl = CP1.Point(); double tolsav =
    // CP1.Tolerance(); CP1.Reset(); CP1.SetPoint(pjyl);
    // CP1.SetTolerance(tolsav); } else if (resetcp2) { ... }`
    // (`cxx:690-705`).
    if ctx.reset_cp1 {
        let pjyl = *ctx.cp1.point();
        let tolsav = ctx.cp1.tolerance();
        ctx.cp1.reset();
        ctx.cp1.set_point(&pjyl);
        ctx.cp1.set_tolerance(tolsav);
    } else if ctx.reset_cp2 {
        let pjyl = *ctx.cp2.point();
        let tolsav = ctx.cp2.tolerance();
        ctx.cp2.reset();
        ctx.cp2.set_point(&pjyl);
        ctx.cp2.set_tolerance(tolsav);
    }

    // `const ChFiDS_CommonPoint& Pf1 = coin->VertexFirstOnS1(); ChFiDS_CommonPoint&
    // Pf2 = coin->ChangeVertexFirstOnS2(); const ChFiDS_CommonPoint& Pl1 =
    // coin->VertexLastOnS1(); ChFiDS_CommonPoint& Pl2 =
    // coin->ChangeVertexLastOnS2(); Pf2 = CP1; Pl2 = CP2;`
    // (`cxx:708-712`).
    let pf1 = coin.vertex(true, 1).clone();
    let pl1 = coin.vertex(false, 1).clone();
    coin.vertex_mut(true, 2).clone_from(&ctx.cp1);
    coin.vertex_mut(false, 2).clone_from(&ctx.cp2);
    let pf2 = ctx.cp1.clone();
    let pl2 = ctx.cp2.clone();

    // `If1 = ChFi3d_IndexPointInDS(Pf1, DStr); If2 = ...; Il1 = ...; if
    // (sameparam) { Il2 = If2; } else { Il2 = ChFi3d_IndexPointInDS(Pl2, DStr);
    // }` (`cxx:718-728`).
    let if1 = ch_fi3d_index_point_in_ds(&pf1, base.my_ds.change_ds());
    let if2 = ch_fi3d_index_point_in_ds(&pf2, base.my_ds.change_ds());
    let il1 = ch_fi3d_index_point_in_ds(&pl1, base.my_ds.change_ds());
    let il2 = if ctx.same_param {
        if2
    } else {
        ch_fi3d_index_point_in_ds(&pl2, base.my_ds.change_ds())
    };

    // `ChFiDS_Regul regdeb, regfin;` (`cxx:716`).
    let mut regdeb = Regul::new();
    let mut regfin = Regul::new();

    // `pp1 = coin->InterferenceOnS1().PCurveOnSurf()->Value(
    // coin->InterferenceOnS1().FirstParameter());` and the same on S2
    // (`cxx:730-734`).
    let pp1 = coin
        .interference_on_s1()
        .p_curve_on_surf()
        .ok_or_else(|| "PerformTwoCorner: coin/1 has no p-curve".to_string())?
        .value(coin.interference_on_s1().first_parameter());
    let pp2 = coin
        .interference_on_s2()
        .p_curve_on_surf()
        .ok_or_else(|| "PerformTwoCorner: coin/2 has no p-curve".to_string())?
        .value(coin.interference_on_s2().first_parameter());

    // `ChFi3d_ComputeArete(Pf1, pp1, Pf2, pp2, DStr.Surface(coin->Surf())
    // .Surface(), C3d, corner->ChangeFirstPCurve(), P1deb, P2deb, tolapp3d,
    // tol2d, tolreached, 0);` (`cxx:735-749`).
    let surf_coin = data_surface(base.my_ds.ds(), coin)?;
    let mut arete = ComputeArete::default();
    ch_fi3d_compute_arete(
        &pf1,
        &pp1,
        &pf2,
        &pp2,
        &surf_coin,
        &mut arete,
        base.tolapp3d,
        base.tol2d,
        0,
    )?;
    let p1deb = arete.pardeb;
    let p2deb = arete.parfin;
    *corner.change_first_p_curve() = arete.p_curve.clone();

    // `double par1 = sd1->Interference(IFaArc1).Parameter(isfirst1); pp1 =
    // sd1->Interference(IFaCo1).PCurveOnSurf()->Value(par1); pp2 =
    // sd1->Interference(IFaArc1).PCurveOnSurf()->Value(par1);`
    // (`cxx:750-752`).
    let par1 = ctx.sd1.interference(ctx.i_fa_arc1).parameter(ctx.is_first1);
    let sp1 = ctx
        .sd1
        .interference(ctx.i_fa_co1)
        .p_curve_on_surf()
        .ok_or_else(|| "PerformTwoCorner: sd1/IFaCo1 has no p-curve".to_string())?
        .value(par1);
    let sp2 = ctx
        .sd1
        .interference(ctx.i_fa_arc1)
        .p_curve_on_surf()
        .ok_or_else(|| "PerformTwoCorner: sd1/IFaArc1 has no p-curve".to_string())?
        .value(par1);

    // `ChFi3d_ComputePCurv(C3d, pp1, pp2, st1->ChangePCurve(isfirst1),
    // DStr.Surface(sd1->Surf()).Surface(), P1deb, P2deb, tolapp3d, tolr1);
    // tolreached = std::max(tolreached, tolr1);`
    // (`cxx:753-764`).
    let surf1 = data_surface(base.my_ds.ds(), &ctx.sd1)?;
    let c3d = arete
        .c3d
        .clone()
        .ok_or_else(|| "PerformTwoCorner: ChFi3d_ComputeArete gave no C3d".to_string())?;
    let mut pc_st1: Option<Arc<dyn Curve2d>> = None;
    let tolr1 = ch_fi3d_compute_pcurv(
        c3d.as_ref(),
        &sp1,
        &sp2,
        &mut pc_st1,
        surf1.as_ref(),
        p1deb,
        p2deb,
        base.tolapp3d,
        false,
    )?;
    let mut tolreached = arete.tolreached.max(tolr1);
    {
        let mut guard = ctx
            .st1
            .write()
            .map_err(|_| "PerformTwoCorner: poisoned stripe".to_string())?;
        *guard.change_p_curve(ctx.is_first1) = pc_st1;
    }

    // `TopOpeBRepDS_Curve Tcurv1(C3d, tolreached); Icf =
    // DStr.AddCurve(Tcurv1); regdeb.SetCurve(Icf); regdeb.SetS1(coin->Surf(),
    // false); regdeb.SetS2(sd1->Surf(), false); myRegul.Append(regdeb);
    // corner->ChangeFirstCurve(Icf); corner->ChangeFirstParameters(P1deb,
    // P2deb); corner->ChangeIndexFirstPointOnS1(If1);
    // corner->ChangeIndexFirstPointOnS2(If2);`
    // (`cxx:765-774`).
    let mut tcurv1 = DsCurve::with_curve(Some(Arc::clone(&c3d)), tolreached, false);
    let icf = base.my_ds.change_ds().add_curve(&mut tcurv1);
    regdeb.set_curve(icf);
    regdeb.set_s1(coin.surf(), false);
    regdeb.set_s2(ctx.sd1.surf(), false);
    base.my_regul.push(regdeb);
    corner.change_first_curve(icf);
    corner.change_first_parameters(p1deb, p2deb);
    corner.change_index_first_point_on_s1(if1);
    corner.change_index_first_point_on_s2(if2);

    // `pp1 = coin->InterferenceOnS1().PCurveOnSurf()->Value(
    // coin->InterferenceOnS1().LastParameter());` and the same on S2
    // (`cxx:776-780`).
    let lp1 = coin
        .interference_on_s1()
        .p_curve_on_surf()
        .ok_or_else(|| "PerformTwoCorner: coin/1 has no p-curve".to_string())?
        .value(coin.interference_on_s1().last_parameter());
    let lp2 = coin
        .interference_on_s2()
        .p_curve_on_surf()
        .ok_or_else(|| "PerformTwoCorner: coin/2 has no p-curve".to_string())?
        .value(coin.interference_on_s2().last_parameter());

    // `ChFi3d_ComputeArete(Pl1, pp1, Pl2, pp2, ..., C3d,
    // corner->ChangeLastPCurve(), P1fin, P2fin, tolapp3d, tol2d, tolreached,
    // 0);` (`cxx:780-791`).
    let mut last_arete = ComputeArete::default();
    ch_fi3d_compute_arete(
        &pl1,
        &lp1,
        &pl2,
        &lp2,
        &surf_coin,
        &mut last_arete,
        base.tolapp3d,
        base.tol2d,
        0,
    )?;
    let p1fin = last_arete.pardeb;
    let p2fin = last_arete.parfin;
    *corner.change_last_p_curve() = last_arete.p_curve.clone();

    // `double par2 = sd2->Interference(IFaArc2).Parameter(isfirst2); ...`
    // (`cxx:792-794`).
    let par2 = ctx.sd2.interference(ctx.i_fa_arc2).parameter(ctx.is_first2);
    let sp1 = ctx
        .sd2
        .interference(ctx.i_fa_co2)
        .p_curve_on_surf()
        .ok_or_else(|| "PerformTwoCorner: sd2/IFaCo2 has no p-curve".to_string())?
        .value(par2);
    let sp2 = ctx
        .sd2
        .interference(ctx.i_fa_arc2)
        .p_curve_on_surf()
        .ok_or_else(|| "PerformTwoCorner: sd2/IFaArc2 has no p-curve".to_string())?
        .value(par2);

    // `ChFi3d_ComputePCurv(C3d, pp1, pp2, st2->ChangePCurve(isfirst2),
    // DStr.Surface(sd2->Surf()).Surface(), P1deb, P2deb, tolapp3d, tolr2);`
    // (`cxx:795-805`). OCCT passes the **first** pair of parameters here, not
    // `P1fin` / `P2fin`; the port keeps that.
    let surf2 = data_surface(base.my_ds.ds(), &ctx.sd2)?;
    let c3d = last_arete
        .c3d
        .clone()
        .ok_or_else(|| "PerformTwoCorner: ChFi3d_ComputeArete gave no C3d".to_string())?;
    let mut pc_st2: Option<Arc<dyn Curve2d>> = None;
    let tolr2 = ch_fi3d_compute_pcurv(
        c3d.as_ref(),
        &sp1,
        &sp2,
        &mut pc_st2,
        surf2.as_ref(),
        p1deb,
        p2deb,
        base.tolapp3d,
        false,
    )?;
    tolreached = tolreached.max(tolr2);
    {
        let mut guard = ctx
            .st2
            .write()
            .map_err(|_| "PerformTwoCorner: poisoned stripe".to_string())?;
        *guard.change_p_curve(ctx.is_first2) = pc_st2;
    }

    // `TopOpeBRepDS_Curve Tcurv2(C3d, tolreached); Icl =
    // DStr.AddCurve(Tcurv2); regfin.SetCurve(Icl); regfin.SetS1(coin->Surf(),
    // false); regfin.SetS2(sd2->Surf(), false); myRegul.Append(regfin);
    // corner->ChangeLastCurve(Icl); corner->ChangeLastParameters(P1fin,
    // P2fin); corner->ChangeIndexLastPointOnS1(Il1);
    // corner->ChangeIndexLastPointOnS2(Il2);`
    // (`cxx:806-815`).
    let mut tcurv2 = DsCurve::with_curve(Some(Arc::clone(&c3d)), tolreached, false);
    let icl = base.my_ds.change_ds().add_curve(&mut tcurv2);
    regfin.set_curve(icl);
    regfin.set_s1(coin.surf(), false);
    regfin.set_s2(ctx.sd2.surf(), false);
    base.my_regul.push(regfin);
    corner.change_last_curve(icl);
    corner.change_last_parameters(p1fin, p2fin);
    corner.change_index_last_point_on_s1(il1);
    corner.change_index_last_point_on_s2(il2);

    // `coin->ChangeIndexOfS1(DStr.AddShape(FaCo)); if (sameparam) {
    // coin->ChangeIndexOfS2(0); } else { coin->ChangeIndexOfS2(DStr.AddShape(
    // FaPiv)); coin->ChangeInterferenceOnS2().SetTransition(pctrans); }
    // corner->SetSolidIndex(st1->SolidIndex());`
    // (`cxx:817-827`).
    let fa_co = ctx.fa_co.0.clone();
    let i_fa_co = base.my_ds.change_ds().add_shape(&fa_co);
    coin.set_index_of_s1(i_fa_co);
    if ctx.same_param {
        coin.set_index_of_s2(0);
    } else {
        let fa_piv = ctx
            .fa_piv
            .as_ref()
            .ok_or_else(|| "PerformTwoCorner: null FaPiv".to_string())?
            .0
            .clone();
        let i_fa_piv = base.my_ds.change_ds().add_shape(&fa_piv);
        coin.set_index_of_s2(i_fa_piv);
        coin.interference_on_s2_mut().set_transition(ctx.pctrans);
    }
    let solid_index = ctx
        .st1
        .read()
        .map_err(|_| "PerformTwoCorner: poisoned stripe".to_string())?
        .solid_index();
    corner.set_solid_index(solid_index);

    // `st1->SetCurve(Icf, isfirst1); st1->SetIndexPoint(If1, isfirst1,
    // IFaCo1); st1->SetIndexPoint(If2, isfirst1, IFaArc1);
    // st1->SetParameters(isfirst1, P1deb, P2deb); sd1->ChangeVertex(isfirst1,
    // IFaCo1) = Pf1; sd1->ChangeVertex(isfirst1, IFaArc1) = Pf2;
    // sd1->ChangeInterference(IFaCo1).SetParameter(par1, isfirst1); if (IFaCo1
    // == 2) { st1->SetOrientation(TopAbs_REVERSED, isfirst1); }`
    // (`cxx:831-841`).
    {
        let mut guard = ctx
            .st1
            .write()
            .map_err(|_| "PerformTwoCorner: poisoned stripe".to_string())?;
        guard.set_curve(icf, ctx.is_first1);
        guard.set_index_point(if1, ctx.is_first1, ctx.i_fa_co1);
        guard.set_index_point(if2, ctx.is_first1, ctx.i_fa_arc1);
        guard.set_parameters(ctx.is_first1, p1deb, p2deb);
        if ctx.i_fa_co1 == 2 {
            guard.set_orientation(Orientation::Reversed, ctx.is_first1);
        }
    }
    ctx.sd1
        .vertex_mut(ctx.is_first1, ctx.i_fa_co1)
        .clone_from(&pf1);
    ctx.sd1
        .vertex_mut(ctx.is_first1, ctx.i_fa_arc1)
        .clone_from(&pf2);
    ctx.sd1
        .interference_mut(ctx.i_fa_co1)
        .set_parameter(par1, ctx.is_first1);

    // `st2->SetCurve(Icl, isfirst2); st2->SetIndexPoint(Il1, isfirst2,
    // IFaCo2); st2->SetIndexPoint(Il2, isfirst2, IFaArc2);
    // st2->SetParameters(isfirst2, P1fin, P2fin); sd2->ChangeVertex(isfirst2,
    // IFaCo2) = Pl1; sd2->ChangeVertex(isfirst2, IFaArc2) = Pl2;
    // sd2->ChangeInterference(IFaCo2).SetParameter(par2, isfirst2); if (IFaCo2
    // == 2) { st2->SetOrientation(TopAbs_REVERSED, isfirst2); }`
    // (`cxx:844-855`).
    {
        let mut guard = ctx
            .st2
            .write()
            .map_err(|_| "PerformTwoCorner: poisoned stripe".to_string())?;
        guard.set_curve(icl, ctx.is_first2);
        guard.set_index_point(il1, ctx.is_first2, ctx.i_fa_co2);
        guard.set_index_point(il2, ctx.is_first2, ctx.i_fa_arc2);
        guard.set_parameters(ctx.is_first2, p1fin, p2fin);
        if ctx.i_fa_co2 == 2 {
            guard.set_orientation(Orientation::Reversed, ctx.is_first2);
        }
    }
    ctx.sd2
        .vertex_mut(ctx.is_first2, ctx.i_fa_co2)
        .clone_from(&pl1);
    ctx.sd2
        .vertex_mut(ctx.is_first2, ctx.i_fa_arc2)
        .clone_from(&pl2);
    ctx.sd2
        .interference_mut(ctx.i_fa_co2)
        .set_parameter(par2, ctx.is_first2);

    Ok(())
}
