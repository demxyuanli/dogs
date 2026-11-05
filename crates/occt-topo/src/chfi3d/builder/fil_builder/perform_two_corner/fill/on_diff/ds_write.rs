//! The DS rewrite of the `OnSame` / `OnDiff` arm:
//! `ChFi3d_FilBuilder_C2.cxx:1006-1144`.
//!
//! `cxx:1010-1022` rebuilds the corner `SurfData`'s four common points from
//! `cpopdif`, `cxx:1032-1045` / `cxx:1076-1088` run the two
//! `ChFi3d_ComputeArete` passes, `cxx:1046-1053` / `cxx:1089-1095` the two
//! `ChFi3d_SameParameter` passes, and `cxx:1054-1141` writes the two DS curves,
//! the two `ChFiDS_Regul` and the two stripes.

use std::sync::Arc;

use occt_geom::curve::Curve;
use occt_geom::Surface;
use occt_geom2d::curve::Curve2d;

use crate::abs::Orientation;
use crate::chfi_ds::{Regul, Stripe, SurfData};
use crate::chfi3d::builder::ChFi3dBuilder;
use crate::chfi3d::builder_0_helpers::{
    ch_fi3d_compute_arete, ch_fi3d_fil_curve_in_ds, ch_fi3d_index_point_in_ds, ComputeArete,
};
use crate::chfi3d::builder_helpers::same_parameter::same_parameter;
use crate::topopebrep::ds::Curve as DsCurve;

use super::SortedCtx;

/// `DStr.Surface(coin->Surf()).Surface()`.
fn data_surface(
    ds: &crate::topopebrep::ds::DataStructure,
    data: &SurfData,
) -> Result<Arc<dyn Surface>, String> {
    ds.surface(data.surf())
        .surface()
        .cloned()
        .ok_or_else(|| "PerformTwoCorner: DStr.Surface(Surf()).Surface() is null".to_string())
}

/// One `ChFi3d_ComputeArete` + `ChFi3d_SameParameter` pair (`cxx:1032-1053`
/// for the first corner, `cxx:1076-1095` for the second).
struct AretePair {
    /// `C3d`
    c3d: Arc<dyn Curve>,
    /// `P1deb`, `P2deb` (or `P1fin`, `P2fin`)
    par_deb: f64,
    /// `P2deb` (or `P2fin`)
    par_fin: f64,
    /// `tolreached`, after the `SameParameter` maximum
    tol_reached: f64,
}

/// Runs `ChFi3d_ComputeArete` for one end of the corner and the
/// `ChFi3d_SameParameter` that follows it, writing the result into the corner
/// stripe's p-curve.
#[allow(clippy::too_many_arguments)]
fn arete_and_same_parameter(
    base: &mut ChFi3dBuilder,
    coin: &SurfData,
    first: bool,
    p1: &crate::chfi_ds::CommonPoint,
    p2: &crate::chfi_ds::CommonPoint,
    p_curve_slot: &mut Option<Arc<dyn Curve2d>>,
    pcurv: &mut Option<Arc<dyn Curve2d>>,
    surf_same: &dyn Surface,
) -> Result<AretePair, String> {
    // `pp1 = coin->InterferenceOnS1().PCurveOnSurf()->Value(
    // coin->InterferenceOnS1().{First,Last}Parameter()); pp2 = ... S2 ...`
    // (`cxx:1029-1031`, `cxx:1072-1075`).
    let u1 = if first {
        coin.interference_on_s1().first_parameter()
    } else {
        coin.interference_on_s1().last_parameter()
    };
    let u2 = if first {
        coin.interference_on_s2().first_parameter()
    } else {
        coin.interference_on_s2().last_parameter()
    };
    let pp1 = coin
        .interference_on_s1()
        .p_curve_on_surf()
        .ok_or_else(|| "PerformTwoCorner: coin/1 has no p-curve".to_string())?
        .value(u1);
    let pp2 = coin
        .interference_on_s2()
        .p_curve_on_surf()
        .ok_or_else(|| "PerformTwoCorner: coin/2 has no p-curve".to_string())?
        .value(u2);

    // `ChFi3d_ComputeArete(P1, pp1, P2, pp2, DStr.Surface(coin->Surf())
    // .Surface(), C3d, corner->Change{First,Last}PCurve(), P1deb, P2deb,
    // tolapp3d, tol2d, tolreached, 0);` (`cxx:1032-1045`, `cxx:1076-1088`).
    let surf_coin = data_surface(base.my_ds.ds(), coin)?;
    let mut arete = ComputeArete::default();
    ch_fi3d_compute_arete(
        p1,
        &pp1,
        p2,
        &pp2,
        &surf_coin,
        &mut arete,
        base.tolapp3d,
        base.tol2d,
        0,
    )?;
    let par_deb = arete.pardeb;
    let par_fin = arete.parfin;
    *p_curve_slot = arete.p_curve.clone();
    let c3d = arete
        .c3d
        .clone()
        .ok_or_else(|| "PerformTwoCorner: ChFi3d_ComputeArete gave no C3d".to_string())?;

    // `occ::handle<GeomAdaptor_Curve> HC3d = new GeomAdaptor_Curve(C3d);
    // ChFi3d_SameParameter(HC3d, pcurv, Hsurf, tolapp3d, tolr);
    // tolreached = std::max(tolreached, tolr);` (`cxx:1046-1050`,
    // `cxx:1089-1093`).
    //
    // `GeomAdaptor_Curve(C3d)` without a window keeps the curve's own range,
    // which is what `same_parameter` reads, so `C3d` goes straight through.
    let (_, tolr) = same_parameter(c3d.as_ref(), pcurv, surf_same, base.tolapp3d);
    let tol_reached = arete.tolreached.max(tolr);

    Ok(AretePair {
        c3d,
        par_deb,
        par_fin,
        tol_reached,
    })
}

/// `cxx:1006-1144`.
pub(super) fn write_ds(
    base: &mut ChFi3dBuilder,
    corner: &mut Stripe,
    coin: &mut SurfData,
    sorted: &mut SortedCtx,
) -> Result<(), String> {
    // `const ChFiDS_CommonPoint& Pf1 = coin->VertexFirstOnS1(); ChFiDS_CommonPoint&
    // Pf2 = coin->ChangeVertexFirstOnS2(); const ChFiDS_CommonPoint& Pl1 =
    // coin->VertexLastOnS1(); ChFiDS_CommonPoint& Pl2 =
    // coin->ChangeVertexLastOnS2(); Pf2 = Pl2 = cpopdif;` (`cxx:1013-1018`).
    let pf1 = coin.vertex(true, 1).clone();
    let pl1 = coin.vertex(false, 1).clone();
    coin.vertex_mut(true, 2).clone_from(&sorted.cp_op_dif);
    coin.vertex_mut(false, 2).clone_from(&sorted.cp_op_dif);
    let pf2 = sorted.cp_op_dif.clone();
    let pl2 = sorted.cp_op_dif.clone();

    // `If1 = ChFi3d_IndexPointInDS(Pf1, DStr); If2 = ...(Pf2); Il1 =
    // ...(Pl1); Il2 = If2;` (`cxx:1023-1026`).
    let if1 = ch_fi3d_index_point_in_ds(&pf1, base.my_ds.change_ds());
    let if2 = ch_fi3d_index_point_in_ds(&pf2, base.my_ds.change_ds());
    let il1 = ch_fi3d_index_point_in_ds(&pl1, base.my_ds.change_ds());
    let il2 = if2;

    // `ChFiDS_Regul regdeb, regfin;` (`cxx:1020`).
    let mut regdeb = Regul::new();
    let mut regfin = Regul::new();

    // `cxx:1032-1053`: the first corner.
    let hbr_fopsam_surface = sorted
        .hbr_fopsam
        .surface()
        .cloned()
        .ok_or_else(|| "ChFi3d_SameParameter: null Adaptor3d_Surface".to_string())?;
    let first = {
        let slot = corner.change_first_p_curve();
        let mut pc_fopsam = sorted.pc_fopsam.clone();
        let pair = arete_and_same_parameter(
            base,
            coin,
            true,
            &pf1,
            &pf2,
            slot,
            &mut pc_fopsam,
            hbr_fopsam_surface.as_ref(),
        )?;
        sorted.pc_fopsam = pc_fopsam;
        pair
    };

    // `TopOpeBRepDS_Curve Tcurv1(C3d, tolreached); Icf = DStr.AddCurve(Tcurv1);
    // OpcFopsam = sdsam->Interference(ifaopsam).Transition(); IFopsam =
    // sdsam->Index(ifaopsam); if (isfirstsam) OpcFopsam =
    // TopAbs::Reverse(OpcFopsam); interf = ChFi3d_FilCurveInDS(Icf, IFopsam,
    // pcFopsam, OpcFopsam); DStr.ChangeShapeInterferences(IFopsam).Append(
    // interf);` (`cxx:1051-1061`).
    let mut tcurv1 = DsCurve::with_curve(Some(Arc::clone(&first.c3d)), first.tol_reached, false);
    let icf = base.my_ds.change_ds().add_curve(&mut tcurv1);
    let i_fopsam = sorted.sd_sam.index(sorted.i_fa_op_sam);
    let mut o_pc_fopsam = sorted.sd_sam.interference(sorted.i_fa_op_sam).transition();
    if sorted.is_first_sam {
        o_pc_fopsam = o_pc_fopsam.reversed();
    }
    let interf = ch_fi3d_fil_curve_in_ds(icf, i_fopsam, sorted.pc_fopsam.clone(), o_pc_fopsam);
    base.my_ds
        .change_ds()
        .change_shape_interferences_at(i_fopsam)
        .push(interf);

    // `regdeb.SetCurve(Icf); regdeb.SetS1(coin->Surf(), false);
    // regdeb.SetS2(IFopsam, true); myRegul.Append(regdeb);
    // corner->ChangeFirstCurve(Icf); corner->ChangeFirstParameters(P1deb,
    // P2deb); corner->ChangeIndexFirstPointOnS1(If1);
    // corner->ChangeIndexFirstPointOnS2(If2);` (`cxx:1063-1070`).
    regdeb.set_curve(icf);
    regdeb.set_s1(coin.surf(), false);
    regdeb.set_s2(i_fopsam, true);
    base.my_regul.push(regdeb);
    corner.change_first_curve(icf);
    corner.change_first_parameters(first.par_deb, first.par_fin);
    corner.change_index_first_point_on_s1(if1);
    corner.change_index_first_point_on_s2(if2);

    // `cxx:1076-1095`: the second corner.
    let h_surf_dif = data_surface(base.my_ds.ds(), &sorted.sd_dif)?;
    let second = {
        let slot = corner.change_last_p_curve();
        let mut pc_surf_dif = sorted.pc_surf_dif.clone();
        let pair = arete_and_same_parameter(
            base,
            coin,
            false,
            &pl1,
            &pl2,
            slot,
            &mut pc_surf_dif,
            h_surf_dif.as_ref(),
        )?;
        sorted.pc_surf_dif = pc_surf_dif;
        pair
    };

    // `TopOpeBRepDS_Curve Tcurv2(C3d, tolreached); Icl = DStr.AddCurve(Tcurv2);
    // regfin.SetCurve(Icl); regfin.SetS1(coin->Surf(), false);
    // regfin.SetS2(sddif->Surf(), false); myRegul.Append(regfin);
    // corner->ChangeLastCurve(Icl); corner->ChangeLastParameters(P1fin,
    // P2fin); corner->ChangeIndexLastPointOnS1(Il1);
    // corner->ChangeIndexLastPointOnS2(Il2);` (`cxx:1094-1102`).
    let mut tcurv2 = DsCurve::with_curve(Some(Arc::clone(&second.c3d)), second.tol_reached, false);
    let icl = base.my_ds.change_ds().add_curve(&mut tcurv2);
    regfin.set_curve(icl);
    regfin.set_s1(coin.surf(), false);
    regfin.set_s2(sorted.sd_dif.surf(), false);
    base.my_regul.push(regfin);
    corner.change_last_curve(icl);
    corner.change_last_parameters(second.par_deb, second.par_fin);
    corner.change_index_last_point_on_s1(il1);
    corner.change_index_last_point_on_s2(il2);

    // `coin->ChangeIndexOfS1(-sdsam->Surf()); coin->ChangeIndexOfS2(0);`
    // (`cxx:1104-1105`).
    coin.set_index_of_s1(-sorted.sd_sam.surf());
    coin.set_index_of_s2(0);

    // `corner->SetSolidIndex(stsam->SolidIndex());` (`cxx:1107`).
    let solid_index = sorted
        .st_sam
        .read()
        .map_err(|_| "PerformTwoCorner: poisoned stripe".to_string())?
        .solid_index();
    corner.set_solid_index(solid_index);

    // `const ChFiDS_FaceInterference& intcoin1 =
    // coin->InterferenceOnS1(); stsam->SetCurve(intcoin1.LineIndex(),
    // isfirstsam); stsam->InDS(isfirstsam); stsam->ChangePCurve(isfirstsam) =
    // coin->InterferenceOnS1().PCurveOnFace(); stsam->SetIndexPoint(If1,
    // isfirstsam, ifaopsam); stsam->SetIndexPoint(Il1, isfirstsam, ifacosam);
    // stsam->SetParameters(isfirstsam, intcoin1.FirstParameter(),
    // intcoin1.LastParameter());` (`cxx:1111-1118`).
    let int_coin1_line = coin.interference_on_s1().line_index();
    let int_coin1_pc = coin.interference_on_s1().p_curve_on_face().cloned();
    let int_coin1_first = coin.interference_on_s1().first_parameter();
    let int_coin1_last = coin.interference_on_s1().last_parameter();
    {
        let mut guard = sorted
            .st_sam
            .write()
            .map_err(|_| "PerformTwoCorner: poisoned stripe".to_string())?;
        guard.set_curve(int_coin1_line, sorted.is_first_sam);
        guard.in_ds(sorted.is_first_sam, 1);
        *guard.change_p_curve(sorted.is_first_sam) = int_coin1_pc;
        guard.set_index_point(if1, sorted.is_first_sam, sorted.i_fa_op_sam);
        guard.set_index_point(il1, sorted.is_first_sam, sorted.i_fa_co_sam);
        guard.set_parameters(sorted.is_first_sam, int_coin1_first, int_coin1_last);
    }

    // `sdsam->ChangeVertex(isfirstsam, ifaopsam) = Pf1;
    // sdsam->ChangeVertex(isfirstsam, ifacosam) = Pl1;
    // sdsam->ChangeInterferenceOnS1().SetParameter(uintpcsam, isfirstsam);
    // sdsam->ChangeInterferenceOnS2().SetParameter(uintpcsam, isfirstsam); if
    // (ifaopsam == 2) stsam->SetOrientation(TopAbs_REVERSED, isfirstsam);`
    // (`cxx:1119-1124`).
    sorted
        .sd_sam
        .vertex_mut(sorted.is_first_sam, sorted.i_fa_op_sam)
        .clone_from(&pf1);
    sorted
        .sd_sam
        .vertex_mut(sorted.is_first_sam, sorted.i_fa_co_sam)
        .clone_from(&pl1);
    sorted
        .sd_sam
        .interference_on_s1_mut()
        .set_parameter(sorted.u_int_pc_sam, sorted.is_first_sam);
    sorted
        .sd_sam
        .interference_on_s2_mut()
        .set_parameter(sorted.u_int_pc_sam, sorted.is_first_sam);
    if sorted.i_fa_op_sam == 2 {
        sorted
            .st_sam
            .write()
            .map_err(|_| "PerformTwoCorner: poisoned stripe".to_string())?
            .set_orientation(Orientation::Reversed, sorted.is_first_sam);
    }

    // `stdif->SetCurve(Icl, isfirstdif); stdif->ChangePCurve(isfirstdif) =
    // pcsurfdif; stdif->SetIndexPoint(Il2, isfirstdif, ifaopdif);
    // stdif->SetIndexPoint(Il1, isfirstdif, ifacodif);
    // stdif->SetParameters(isfirstdif, P1fin, P2fin);` (`cxx:1129-1133`).
    {
        let mut guard = sorted
            .st_dif
            .write()
            .map_err(|_| "PerformTwoCorner: poisoned stripe".to_string())?;
        guard.set_curve(icl, sorted.is_first_dif);
        *guard.change_p_curve(sorted.is_first_dif) = sorted.pc_surf_dif.clone();
        guard.set_index_point(il2, sorted.is_first_dif, sorted.i_fa_op_dif);
        guard.set_index_point(il1, sorted.is_first_dif, sorted.i_fa_co_dif);
        guard.set_parameters(sorted.is_first_dif, second.par_deb, second.par_fin);
    }

    // `sddif->ChangeVertex(isfirstdif, ifaopdif) = Pl2;
    // sddif->ChangeVertex(isfirstdif, ifacodif) = Pl1;
    // sddif->ChangeInterference(ifacodif).SetParameter(uintpcdif, isfirstdif);
    // if (ifaopdif == 1) stdif->SetOrientation(TopAbs_REVERSED, isfirstdif);`
    // (`cxx:1134-1140`).
    sorted
        .sd_dif
        .vertex_mut(sorted.is_first_dif, sorted.i_fa_op_dif)
        .clone_from(&pl2);
    sorted
        .sd_dif
        .vertex_mut(sorted.is_first_dif, sorted.i_fa_co_dif)
        .clone_from(&pl1);
    sorted
        .sd_dif
        .interference_mut(sorted.i_fa_co_dif)
        .set_parameter(sorted.u_int_pc_dif, sorted.is_first_dif);
    if sorted.i_fa_op_dif == 1 {
        sorted
            .st_dif
            .write()
            .map_err(|_| "PerformTwoCorner: poisoned stripe".to_string())?
            .set_orientation(Orientation::Reversed, sorted.is_first_dif);
    }

    Ok(())
}
