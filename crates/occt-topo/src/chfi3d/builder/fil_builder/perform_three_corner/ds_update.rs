//! `ChFi3d_FilBuilder::PerformThreeCorner` - the DS update of the corner, and
//! the truncation of the three stripe sequences.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_FilBuilder_C3.cxx:1016-1331`
//! (OCCT 8.0.0):
//!
//! - `cxx:1016-1021` - the two `p[deb][fin]` / `p[fin][deb]` assignments of the
//!   `!c1pointu` case: [`super::perform_three_corner`],
//! - `cxx:1022-1288` - the four stripes and the DS: [`write_ds`],
//! - `cxx:1289-1331` - the three `Remove` calls and the two appends:
//!   [`trim_and_append`].
//!
//! Translation notes:
//! - `ChFiDS_SurfData& fddeb` / `fdfin` / `fdpiv` are reached through
//!   [`super::with_sd`] and [`super::sd`], which lock the stripe and take
//!   `Arc::make_mut` on the slot, so the aliasing OCCT gets from its handles
//!   (`fddeb == fdpiv` when the two stripes are one) is preserved,
//! - `Pl1` of `cxx:1036-1042` is a *reference* to `coin->ChangeVertexLastOnS1()`
//!   that is reassigned from `VertexFirstOnS1()` when `c1pointu`; the port
//!   copies the first common point into the last one and then reads it, which
//!   is what the C++ reference sees,
//! - `ChFi3d_EnlargeBox(DStr, Stripe, SurfData, box1, box2, isfirst)` is
//!   [`ch_fi3d_enlarge_box_surf_data`]; the four `Bnd_Box` are [`Boxes`], whose
//!   alias map covers `pbl1 = pbf1` (`cxx:1046-1051`),
//! - `ChFi3d_SameParameter(C3d, C2d, S, First, Last, Tol3d, TolReached)` is
//!   `same_parameter_range` of the helpers, the trimmed-curve overload,
//! - `CD[pivot]->InDS(isfirst)` is `in_ds(isfirst, 1)`, the default `Nb`.

use std::sync::Arc;

use occt_core::bnd::BndBox;
use occt_core::gp::{GpPnt, GpPnt2d};
use occt_geom::Surface;
use occt_geom2d::curve::Curve2d;

use crate::abs::Orientation;
use crate::chfi_ds::{Regul, Stripe, SurfData};
use crate::chfi3d::builder_0_helpers::{
    ch_fi3d_compute_arete, ch_fi3d_compute_pcurv, ch_fi3d_enlarge_box_surf_data,
    ch_fi3d_index_point_in_ds, ch_fi3d_set_point_tolerance, ComputeArete,
};
use crate::chfi3d::builder_helpers::same_parameter_range;
use crate::chfi3d::builder::ChFi3dBuilder;
use crate::topopebrep::ds::{Curve as DsCurve, DataStructure};

use super::ThreeCorner;
use super::super::super::corner_sd::{
    data_surface, ds_curve, pcurve_on_surf, raise_curve_tolerance, sd, with_sd,
};

/// The four `Bnd_Box` of `cxx:1046-1054` and the alias map `pbf1` / `pbl1` /
/// `pbf2` / `pbl2`, where `pbl1` is `pbf1` when `c1pointu`.
pub(super) struct Boxes {
    /// `bf1`, `bl1`, `bf2`, `bl2` in that order.
    pub b: [BndBox; 4],
    /// `pbf1`, `pbl1`, `pbf2`, `pbl2` as indices of [`Self::b`].
    pub p: [usize; 4],
}

impl Boxes {
    /// `pbf1 = &bf1, pbl1 = &bl1, ...; if (c1pointu) pbl1 = pbf1;`
    /// (`cxx:1046-1051`).
    pub fn new(c1pointu: bool) -> Self {
        Self {
            b: [BndBox::new(), BndBox::new(), BndBox::new(), BndBox::new()],
            p: [0, if c1pointu { 0 } else { 1 }, 2, 3],
        }
    }

    /// The two boxes of one `ChFi3d_EnlargeBox` call.
    pub fn pair(&mut self, a: usize, b: usize) -> (&mut BndBox, &mut BndBox) {
        if a < b {
            let (lo, hi) = self.b.split_at_mut(b);
            (&mut lo[a], &mut hi[0])
        } else {
            let (lo, hi) = self.b.split_at_mut(a);
            (&mut hi[0], &mut lo[b])
        }
    }

    /// `pbf1->Add(Pf1.Point()); pbf2->Add(Pf2.Point()); pbl1->Add(Pl1.Point());
    /// pbl2->Add(Pl2.Point());` (`cxx:1052-1055`).
    pub fn add_all(&mut self, pf1: &GpPnt, pf2: &GpPnt, pl1: &GpPnt, pl2: &GpPnt) {
        self.b[self.p[0]].add_point(pf1);
        self.b[self.p[2]].add_point(pf2);
        self.b[self.p[1]].add_point(pl1);
        self.b[self.p[3]].add_point(pl2);
    }
}

/// `cxx:1022-1288`, run when `done`.
pub(super) fn write_ds(
    base: &mut ChFi3dBuilder,
    c: &mut ThreeCorner,
    corner: &mut Stripe,
    coin: &mut SurfData,
) -> Result<(), String> {
    // `const ChFiDS_CommonPoint& Pf1 = coin->VertexFirstOnS1();
    // const ChFiDS_CommonPoint& Pf2 = coin->VertexFirstOnS2();
    // ChFiDS_CommonPoint& Pl1 = coin->ChangeVertexLastOnS1();
    // if (c1pointu) { Pl1 = coin->ChangeVertexFirstOnS1(); }
    // const ChFiDS_CommonPoint& Pl2 = coin->VertexLastOnS2();`
    // (`cxx:1035-1043`).
    let pf1 = coin.vertex(true, 1).clone();
    let pf2 = coin.vertex(true, 2).clone();
    if c.c1pointu {
        let first_on_s1 = coin.vertex(true, 1).clone();
        coin.vertex_mut(false, 1).clone_from(&first_on_s1);
    }
    let pl1 = coin.vertex(false, 1).clone();
    let pl2 = coin.vertex(false, 2).clone();

    // `Bnd_Box bf1, bl1, bf2, bl2; if (c1pointu) { pbl1 = pbf1; }
    // pbf1->Add(Pf1.Point()); pbf2->Add(Pf2.Point()); pbl1->Add(Pl1.Point());
    // pbl2->Add(Pl2.Point());` (`cxx:1045-1055`).
    let mut boxes = Boxes::new(c.c1pointu);
    boxes.add_all(pf1.point(), pf2.point(), pl1.point(), pl2.point());
    let (pbf1, pbl1, pbf2, pbl2) = (boxes.p[0], boxes.p[1], boxes.p[2], boxes.p[3]);

    // `ChFiDS_Regul regdeb, regfin;` (`cxx:1058`).
    let mut regdeb = Regul::new();
    let mut regfin = Regul::new();

    // `If1 = ChFi3d_IndexPointInDS(Pf1, DStr); If2 = ...; if (c1pointu) { Il1 =
    // If1; } else { Il1 = ChFi3d_IndexPointInDS(Pl1, DStr); } Il2 = ...;`
    // (`cxx:1059-1068`).
    let if1 = ch_fi3d_index_point_in_ds(&pf1, base.my_ds.change_ds());
    let if2 = ch_fi3d_index_point_in_ds(&pf2, base.my_ds.change_ds());
    let il1 = if c.c1pointu {
        if1
    } else {
        ch_fi3d_index_point_in_ds(&pl1, base.my_ds.change_ds())
    };
    let il2 = ch_fi3d_index_point_in_ds(&pl2, base.my_ds.change_ds());

    // `pp1 = coin->InterferenceOnS1().PCurveOnSurf()->Value(
    // coin->InterferenceOnS1().FirstParameter()); pp2 = ... S2 ...`
    // (`cxx:1070-1071`).
    let pp1 = pcurve_on_surf(coin, 1, coin.interference_on_s1().first_parameter())?;
    let pp2 = pcurve_on_surf(coin, 2, coin.interference_on_s2().first_parameter())?;

    // `if (c1pointu) { coin->ChangeIndexOfS1(0); } else {
    // coin->ChangeIndexOfS1(DStr.AddShape(face[pivot])); }
    // coin->ChangeIndexOfS2(-fdpiv->Surf());` (`cxx:1072-1080`).
    let index_of_s1 = if c.c1pointu {
        0
    } else {
        let pivot_face = c.face[c.pivot]
            .clone()
            .ok_or_else(|| "PerformThreeCorner: the pivot has no face".to_string())?
            .0
            .clone();
        base.my_ds.change_ds().add_shape(&pivot_face)
    };
    coin.set_index_of_s1(index_of_s1);
    let fdpiv_surf = sd(&c.stripes[c.pivot], c.i[c.pivot][c.deb])?.surf();
    coin.set_index_of_s2(-fdpiv_surf);

    // `ChFi3d_ComputeArete(Pf1, pp1, Pf2, pp2, DStr.Surface(coin->Surf())
    // .Surface(), C3d, corner->ChangeFirstPCurve(), P1deb, P2deb, tolapp3d,
    // tol2d, tolreached, 0);` (`cxx:1083-1096`).
    let surf_coin = data_surface(base.my_ds.ds(), coin)?;
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
    let p1deb = first.pardeb;
    let p2deb = first.parfin;
    *corner.change_first_p_curve() = first.p_curve.clone();

    // `TopOpeBRepDS_Curve Tcurv1(C3d, tolreached); Icf =
    // DStr.AddCurve(Tcurv1); regdeb.SetCurve(Icf); regdeb.SetS1(coin->Surf(),
    // false);
    let c3d_first = first
        .c3d
        .clone()
        .ok_or_else(|| "PerformThreeCorner: ChFi3d_ComputeArete gave no C3d".to_string())?;
    let tol_first = first.tolreached;
    let mut tcurv1 = DsCurve::with_curve(Some(Arc::clone(&c3d_first)), tol_first, false);
    let icf = base.my_ds.change_ds().add_curve(&mut tcurv1);
    regdeb.set_curve(icf);
    regdeb.set_s1(coin.surf(), false);
    let fddeb_surf = sd(&c.stripes[c.deb], c.i[c.deb][c.pivot])?.surf();
    regdeb.set_s2(fddeb_surf, false);
    base.my_regul.push(regdeb);
    // `corner->ChangeFirstCurve(Icf); corner->ChangeFirstParameters(P1deb,
    // P2deb); corner->ChangeIndexFirstPointOnS1(If1);
    // corner->ChangeIndexFirstPointOnS2(If2);` (`cxx:1104-1107`).
    corner.change_first_curve(icf);
    corner.change_first_parameters(p1deb, p2deb);
    corner.change_index_first_point_on_s1(if1);
    corner.change_index_first_point_on_s2(if2);

    // `ChFi3d_EnlargeBox(DStr, corner, coin, *pbf1, *pbf2, true);`
    // (`cxx:1108`).
    {
        let (b1, b2) = boxes.pair(pbf1, pbf2);
        ch_fi3d_enlarge_box_surf_data(
            base.my_ds.ds(),
            Some(&*corner),
            coin,
            b1,
            b2,
            true,
        )?;
    }

    // `pp1 = coin->InterferenceOnS1().PCurveOnSurf()->Value(
    // coin->InterferenceOnS1().LastParameter()); pp2 = ... S2 ...`
    // (`cxx:1110-1111`).
    let lp1 = pcurve_on_surf(coin, 1, coin.interference_on_s1().last_parameter())?;
    let lp2 = pcurve_on_surf(coin, 2, coin.interference_on_s2().last_parameter())?;

    // `ChFi3d_ComputeArete(Pl1, pp1, Pl2, pp2, ..., C3d,
    // corner->ChangeLastPCurve(), P1fin, P2fin, tolapp3d, tol2d, tolreached,
    // 0);` (`cxx:1112-1122`).
    let mut last = ComputeArete::default();
    ch_fi3d_compute_arete(
        &pl1,
        &lp1,
        &pl2,
        &lp2,
        &surf_coin,
        &mut last,
        base.tolapp3d,
        base.tol2d,
        0,
    )?;
    let p1fin = last.pardeb;
    let p2fin = last.parfin;
    *corner.change_last_p_curve() = last.p_curve.clone();

    // `TopOpeBRepDS_Curve Tcurv2(C3d, tolreached); Icl =
    // DStr.AddCurve(Tcurv2); regfin.SetCurve(Icl);
    // regfin.SetS1(coin->Surf(), false); regfin.SetS2(fdfin->Surf(), false);
    // myRegul.Append(regfin);` (`cxx:1123-1129`).
    let c3d_last = last
        .c3d
        .clone()
        .ok_or_else(|| "PerformThreeCorner: ChFi3d_ComputeArete gave no C3d".to_string())?;
    let tol_last = last.tolreached;
    let mut tcurv2 = DsCurve::with_curve(Some(Arc::clone(&c3d_last)), tol_last, false);
    let icl = base.my_ds.change_ds().add_curve(&mut tcurv2);
    regfin.set_curve(icl);
    regfin.set_s1(coin.surf(), false);
    let fdfin = sd(&c.stripes[c.fin], c.i[c.fin][c.pivot])?;
    regfin.set_s2(fdfin.surf(), false);
    base.my_regul.push(regfin);
    // `corner->ChangeLastCurve(Icl); corner->ChangeLastParameters(P1fin,
    // P2fin); corner->ChangeIndexLastPointOnS1(Il1);
    // corner->ChangeIndexLastPointOnS2(Il2);` (`cxx:1130-1133`).
    corner.change_last_curve(icl);
    corner.change_last_parameters(p1fin, p2fin);
    corner.change_index_last_point_on_s1(il1);
    corner.change_index_last_point_on_s2(il2);

    // `ChFi3d_EnlargeBox(DStr, corner, coin, *pbl1, *pbl2, false);`
    // (`cxx:1134`).
    {
        let (b1, b2) = boxes.pair(pbl1, pbl2);
        ch_fi3d_enlarge_box_surf_data(
            base.my_ds.ds(),
            Some(&*corner),
            coin,
            b1,
            b2,
            false,
        )?;
    }

    // ---------------------------------------------------------------
    // `// then CornerData of the beginning,` (`cxx:1136-1182`).
    // ---------------------------------------------------------------
    // `bool isfirst = (sens[deb] == 1), rev = (jf[deb][fin] == 2); int isurf1 =
    // 1, isurf2 = 2; parpp1 = p[deb][fin]; parpp2 = p[deb][pivot]; if (rev) {
    // isurf1 = 2; isurf2 = 1; parpp1 = p[deb][pivot]; parpp2 = p[deb][fin];
    // CD[deb]->SetOrientation(TopAbs_REVERSED, isfirst); }` (`cxx:1138-1149`).
    let is_first = c.sens[c.deb] == 1;
    let rev = c.jf[c.deb][c.fin] == 2;
    let (mut isurf1, mut isurf2) = (1, 2);
    let mut parpp1 = c.p[c.deb][c.fin];
    let mut parpp2 = c.p[c.deb][c.pivot];
    if rev {
        isurf1 = 2;
        isurf2 = 1;
        parpp1 = c.p[c.deb][c.pivot];
        parpp2 = c.p[c.deb][c.fin];
        c.stripes[c.deb]
            .write()
            .map_err(|_| "PerformThreeCorner: poisoned stripe".to_string())?
            .set_orientation(Orientation::Reversed, is_first);
    }
    let fddeb = sd(&c.stripes[c.deb], c.i[c.deb][c.pivot])?;

    // `pp1 = fddeb->InterferenceOnS1().PCurveOnSurf()->Value(parpp1); pp2 =
    // fddeb->InterferenceOnS2().PCurveOnSurf()->Value(parpp2);`
    // (`cxx:1150-1151`).
    let dp1 = pcurve_on_surf(&fddeb, 1, parpp1)?;
    let dp2 = pcurve_on_surf(&fddeb, 2, parpp2)?;

    // `CD[deb]->SetCurve(Icf, isfirst); CD[deb]->SetIndexPoint(If1, isfirst,
    // isurf1); CD[deb]->SetIndexPoint(If2, isfirst, isurf2);
    // CD[deb]->SetParameters(isfirst, P1deb, P2deb);` (`cxx:1152-1155`).
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
    // isurf2) = Pf2; fddeb->ChangeInterferenceOnS1().SetParameter(parpp1,
    // isfirst); fddeb->ChangeInterferenceOnS2().SetParameter(parpp2, isfirst);`
    // (`cxx:1156-1159`).
    with_sd(&c.stripes[c.deb], c.i[c.deb][c.pivot], |data| {
        data.vertex_mut(is_first, isurf1).clone_from(&pf1);
        data.vertex_mut(is_first, isurf2).clone_from(&pf2);
        data.interference_on_s1_mut().set_parameter(parpp1, is_first);
        data.interference_on_s2_mut().set_parameter(parpp2, is_first);
    })?;

    // `TopOpeBRepDS_Curve& tcdeb = DStr.ChangeCurve(Icf);
    // occ::handle<Geom_Curve> crefdeb = tcdeb.Curve(); ...
    // ChFi3d_ComputePCurv(crefdeb, pp1, pp2, CD[deb]->ChangePCurve(isfirst),
    // DStr.Surface(fddeb->Surf()).Surface(), P1deb, P2deb, tolapp3d, tolrdeb,
    // rev); tcdeb.Tolerance(std::max(tolrdeb, tcdeb.Tolerance()));`
    // (`cxx:1160-1174`).
    let crefdeb = ds_curve(base, icf)?;
    let fddeb = sd(&c.stripes[c.deb], c.i[c.deb][c.pivot])?;
    let surf_deb = data_surface(base.my_ds.ds(), &fddeb)?;
    let mut pc_deb: Option<Arc<dyn Curve2d>> = None;
    let tolrdeb = ch_fi3d_compute_pcurv(
        crefdeb.as_ref(),
        &dp1,
        &dp2,
        &mut pc_deb,
        surf_deb.as_ref(),
        p1deb,
        p2deb,
        base.tolapp3d,
        rev,
    )?;
    {
        let mut guard = c.stripes[c.deb]
            .write()
            .map_err(|_| "PerformThreeCorner: poisoned stripe".to_string())?;
        *guard.change_p_curve(is_first) = pc_deb;
    }
    raise_curve_tolerance(base, icf, tolrdeb);

    // `if (rev) { ChFi3d_EnlargeBox(DStr, CD[deb], fddeb, *pbf2, *pbf1,
    // isfirst); } else { ... *pbf1, *pbf2 ... }` (`cxx:1175-1182`).
    {
        let guard = c.stripes[c.deb]
            .read()
            .map_err(|_| "PerformThreeCorner: poisoned stripe".to_string())?;
        let data = guard
            .set_of_surf_data()
            .get(c.i[c.deb][c.pivot].saturating_sub(1) as usize)
            .ok_or_else(|| "PerformThreeCorner: SurfData is out of range".to_string())?;
        let (a, b) = if rev { (pbf2, pbf1) } else { (pbf1, pbf2) };
        let (b1, b2) = boxes.pair(a, b);
        ch_fi3d_enlarge_box_surf_data(
            base.my_ds.ds(),
            Some(&*guard),
            data.as_ref(),
            b1,
            b2,
            is_first,
        )?;
    }

    // ---------------------------------------------------------------
    // `// then the end CornerData,` (`cxx:1184-1231`).
    // ---------------------------------------------------------------
    // `isfirst = (sens[fin] == 1); rev = (jf[fin][deb] == 2); isurf1 = 1;
    // isurf2 = 2; parpp1 = p[fin][deb]; parpp2 = p[fin][pivot];`
    // (`cxx:1186-1191`).
    let is_first = c.sens[c.fin] == 1;
    let rev = c.jf[c.fin][c.deb] == 2;
    let (mut isurf1, mut isurf2) = (1, 2);
    let mut parpp1 = c.p[c.fin][c.deb];
    let mut parpp2 = c.p[c.fin][c.pivot];
    if rev {
        isurf1 = 2;
        isurf2 = 1;
        parpp1 = c.p[c.fin][c.pivot];
        parpp2 = c.p[c.fin][c.deb];
        c.stripes[c.fin]
            .write()
            .map_err(|_| "PerformThreeCorner: poisoned stripe".to_string())?
            .set_orientation(Orientation::Reversed, is_first);
    }
    let fdfin_data = sd(&c.stripes[c.fin], c.i[c.fin][c.pivot])?;

    // `pp1 = fdfin->InterferenceOnS1().PCurveOnSurf()->Value(parpp1); pp2 =
    // fdfin->InterferenceOnS2().PCurveOnSurf()->Value(parpp2);`
    // (`cxx:1199-1200`).
    let fp1 = pcurve_on_surf(&fdfin_data, 1, parpp1)?;
    let fp2 = pcurve_on_surf(&fdfin_data, 2, parpp2)?;

    // `CD[fin]->SetCurve(Icl, isfirst); CD[fin]->SetIndexPoint(Il1, isfirst,
    // isurf1); CD[fin]->SetIndexPoint(Il2, isfirst, isurf2);
    // CD[fin]->SetParameters(isfirst, P1fin, P2fin);` (`cxx:1201-1204`).
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
    // isurf2) = Pl2; fdfin->ChangeInterferenceOnS1().SetParameter(parpp1,
    // isfirst); fdfin->ChangeInterferenceOnS2().SetParameter(parpp2, isfirst);`
    // (`cxx:1205-1208`).
    with_sd(&c.stripes[c.fin], c.i[c.fin][c.pivot], |data| {
        data.vertex_mut(is_first, isurf1).clone_from(&pl1);
        data.vertex_mut(is_first, isurf2).clone_from(&pl2);
        data.interference_on_s1_mut().set_parameter(parpp1, is_first);
        data.interference_on_s2_mut().set_parameter(parpp2, is_first);
    })?;

    // `TopOpeBRepDS_Curve& tcfin = DStr.ChangeCurve(Icl); ...
    // ChFi3d_ComputePCurv(creffin, pp1, pp2, CD[fin]->ChangePCurve(isfirst),
    // DStr.Surface(fdfin->Surf()).Surface(), P1fin, P2fin, tolapp3d, tolrfin,
    // rev); tcfin.Tolerance(std::max(tolrfin, tcfin.Tolerance()));`
    // (`cxx:1209-1223`).
    let creffin = ds_curve(base, icl)?;
    let fdfin_data = sd(&c.stripes[c.fin], c.i[c.fin][c.pivot])?;
    let surf_fin = data_surface(base.my_ds.ds(), &fdfin_data)?;
    let mut pc_fin: Option<Arc<dyn Curve2d>> = None;
    let tolrfin = ch_fi3d_compute_pcurv(
        creffin.as_ref(),
        &fp1,
        &fp2,
        &mut pc_fin,
        surf_fin.as_ref(),
        p1fin,
        p2fin,
        base.tolapp3d,
        rev,
    )?;
    {
        let mut guard = c.stripes[c.fin]
            .write()
            .map_err(|_| "PerformThreeCorner: poisoned stripe".to_string())?;
        *guard.change_p_curve(is_first) = pc_fin;
    }
    raise_curve_tolerance(base, icl, tolrfin);

    // `if (rev) { ChFi3d_EnlargeBox(DStr, CD[fin], fdfin, *pbl2, *pbl1,
    // isfirst); } else { ... *pbl1, *pbl2 ... }` (`cxx:1224-1230`).
    {
        let guard = c.stripes[c.fin]
            .read()
            .map_err(|_| "PerformThreeCorner: poisoned stripe".to_string())?;
        let data = guard
            .set_of_surf_data()
            .get(c.i[c.fin][c.pivot].saturating_sub(1) as usize)
            .ok_or_else(|| "PerformThreeCorner: SurfData is out of range".to_string())?;
        let (a, b) = if rev { (pbl2, pbl1) } else { (pbl1, pbl2) };
        let (b1, b2) = boxes.pair(a, b);
        ch_fi3d_enlarge_box_surf_data(
            base.my_ds.ds(),
            Some(&*guard),
            data.as_ref(),
            b1,
            b2,
            is_first,
        )?;
    }

    // ---------------------------------------------------------------
    // `// anf finally the pivot.` (`cxx:1233-1278`).
    // ---------------------------------------------------------------
    // `isfirst = (sens[pivot] == 1); rev = (jf[pivot][deb] == 2); isurf1 = 1;
    // isurf2 = 2; if (rev) { isurf1 = 2; isurf2 = 1;
    // CD[pivot]->SetOrientation(TopAbs_REVERSED, isfirst); }`
    // (`cxx:1234-1245`).
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

    // `int ICcoinpiv = fi.LineIndex(); TopOpeBRepDS_Curve& TCcoinpiv =
    // DStr.ChangeCurve(ICcoinpiv); CD[pivot]->SetCurve(ICcoinpiv, isfirst);
    // Ccoinpiv = DStr.Curve(ICcoinpiv).Curve();`
    // (`cxx:1246-1249`).
    let ic_coinpiv = coin.interference_on_s2().line_index();
    c.stripes[c.pivot]
        .write()
        .map_err(|_| "PerformThreeCorner: poisoned stripe".to_string())?
        .set_curve(ic_coinpiv, is_first);
    let c_coinpiv = ds_curve(base, ic_coinpiv)?;

    // `occ::handle<Geom2d_Curve>& C2dOnPiv = fi.ChangePCurveOnFace();
    // Spiv = DStr.Surface(fdpiv->Surf()).Surface(); double tolr;
    // ChFi3d_SameParameter(Ccoinpiv, C2dOnPiv, Spiv, fi.FirstParameter(),
    // fi.LastParameter(), tolapp3d, tolr);
    // TCcoinpiv.Tolerance(std::max(TCcoinpiv.Tolerance(), tolr));`
    // (`cxx:1249-1259`).
    let mut c2d_on_piv = coin.interference_on_s2().p_curve_on_face().cloned();
    let fdpiv_data = sd(&c.stripes[c.pivot], c.i[c.pivot][c.deb])?;
    let spiv = data_surface(base.my_ds.ds(), &fdpiv_data)?;
    let fi_first = coin.interference_on_s2().first_parameter();
    let fi_last = coin.interference_on_s2().last_parameter();
    let (_, tolr) = same_parameter_range(
        c_coinpiv.as_ref(),
        &mut c2d_on_piv,
        spiv.as_ref(),
        fi_first,
        fi_last,
        base.tolapp3d,
    );
    raise_curve_tolerance(base, ic_coinpiv, tolr);

    // `CD[pivot]->ChangePCurve(isfirst) = C2dOnPiv; CD[pivot]->SetIndexPoint(
    // If2, isfirst, isurf1); CD[pivot]->SetIndexPoint(Il2, isfirst, isurf2);
    // CD[pivot]->SetParameters(isfirst, fi.FirstParameter(),
    // fi.LastParameter());` (`cxx:1260-1263`).
    coin
        .interference_on_s2_mut()
        .set_p_curve_on_face(c2d_on_piv.clone());
    {
        let mut guard = c.stripes[c.pivot]
            .write()
            .map_err(|_| "PerformThreeCorner: poisoned stripe".to_string())?;
        *guard.change_p_curve(is_first) = c2d_on_piv;
        guard.set_index_point(if2, is_first, isurf1);
        guard.set_index_point(il2, is_first, isurf2);
        guard.set_parameters(is_first, fi_first, fi_last);
        // `CD[pivot]->InDS(isfirst);` (`cxx:1268`).
        guard.in_ds(is_first, 1);
    }

    // `fdpiv->ChangeVertex(isfirst, isurf1) = Pf2; fdpiv->ChangeVertex(isfirst,
    // isurf2) = Pl2; fdpiv->ChangeInterference(isurf1).SetParameter(
    // p[pivot][deb], isfirst); fdpiv->ChangeInterference(isurf2).SetParameter(
    // p[pivot][fin], isfirst);` (`cxx:1264-1267`).
    let piv_deb = c.p[c.pivot][c.deb];
    let piv_fin = c.p[c.pivot][c.fin];
    with_sd(&c.stripes[c.pivot], c.i[c.pivot][c.deb], |data| {
        data.vertex_mut(is_first, isurf1).clone_from(&pf2);
        data.vertex_mut(is_first, isurf2).clone_from(&pl2);
        data.interference_mut(isurf1).set_parameter(piv_deb, is_first);
        data.interference_mut(isurf2).set_parameter(piv_fin, is_first);
    })?;

    // `if (rev) { ChFi3d_EnlargeBox(DStr, CD[pivot], fdpiv, *pbl2, *pbf2,
    // isfirst); } else { ... *pbf2, *pbl2 ... }` (`cxx:1269-1276`).
    {
        let guard = c.stripes[c.pivot]
            .read()
            .map_err(|_| "PerformThreeCorner: poisoned stripe".to_string())?;
        let data = guard
            .set_of_surf_data()
            .get(c.i[c.pivot][c.deb].saturating_sub(1) as usize)
            .ok_or_else(|| "PerformThreeCorner: SurfData is out of range".to_string())?;
        let (a, b) = if rev { (pbl2, pbf2) } else { (pbf2, pbl2) };
        let (b1, b2) = boxes.pair(a, b);
        ch_fi3d_enlarge_box_surf_data(
            base.my_ds.ds(),
            Some(&*guard),
            data.as_ref(),
            b1,
            b2,
            is_first,
        )?;
    }

    // `ChFi3d_SetPointTolerance(DStr, *pbf1, If1); ... *pbl2, Il2);`
    // (`cxx:1279-1283`).
    ch_fi3d_set_point_tolerance(base.my_ds.change_ds(), &boxes.b[pbf1], if1);
    ch_fi3d_set_point_tolerance(base.my_ds.change_ds(), &boxes.b[pbf2], if2);
    ch_fi3d_set_point_tolerance(base.my_ds.change_ds(), &boxes.b[pbl1], il1);
    ch_fi3d_set_point_tolerance(base.my_ds.change_ds(), &boxes.b[pbl2], il2);

    Ok(())
}

/// `cxx:1289-1331`: the three `Remove` calls and the two appends.
pub(super) fn trim_and_append(
    base: &mut ChFi3dBuilder,
    c: &mut ThreeCorner,
    corner: &mut Stripe,
    coin: &SurfData,
) -> Result<(), String> {
    // `if (i[deb][pivot] < Index[deb]) { CD[deb]->...Remove(i[deb][pivot] + 1,
    // Index[deb]); Index[deb] = i[deb][pivot]; } else if (i[deb][pivot] >
    // Index[deb]) { ...Remove(Index[deb], i[deb][pivot] - 1); i[deb][pivot] =
    // Index[deb]; }` (`cxx:1289-1298`), then the same for `fin`
    // (`cxx:1299-1308`) and `pivot` (`cxx:1309-1318`). `NCollection_HSequence
    // ::Remove(From, To)` is 1-based and inclusive.
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
    // myListStripe.Append(corner);` (`cxx:1319-1325`).
    base.change_evi(&c.vtx).push(coin.surf());
    let solid_index = c.stripes[c.pivot]
        .read()
        .map_err(|_| "PerformThreeCorner: poisoned stripe".to_string())?
        .solid_index();
    corner.set_solid_index(solid_index);
    // OCCT appended `coin` to `corner->ChangeSetOfSurfData()` when the corner
    // was created (`cxx:604-608`); the port kept it local, so the sequence
    // catches up here, before the stripe is shared.
    corner.change_set_of_surf_data().push(Arc::new(coin.clone()));
    base.my_list_stripe
        .push(Arc::new(std::sync::RwLock::new(corner.clone())));

    Ok(())
}
