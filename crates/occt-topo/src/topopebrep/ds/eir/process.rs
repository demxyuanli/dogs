//! `FUN_ProcessEdgeInterferences` and `CheckInterferenceIsValid` — the body of
//! `TopOpeBRepDS_EIR::ProcessEdgeInterferences(EIX)`.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_EIR.cxx` (OCCT 8.0.0):
//! - `FUN_ProcessEdgeInterferences` (`:827-958`)
//! - `CheckInterferenceIsValid` (`:1163-1221`)
//!
//! The two geometry kinds `K` and `G` of the C++ signature are unnamed in the
//! definition (`cxx:829-830`), so they are unused by the body.
//!
//! `FUN_unkeepEVIonGb1` is only applied when `EIX` is a section edge
//! (`cxx:953-956`); the two commented-out calls at `cxx:889-890` stay commented
//! out here.

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, TopoShape, Vertex};

use super::super::data_structure::DataStructure;
use super::super::interference::{Interference, InterferenceData};
use super::super::interference_iterator::InterferenceList;
use super::super::kind::Kind;
use super::super::process_interferences_tool::{
    fun_select_1d_i, fun_select_2d_i, fun_select_pure2d_i, fun_select_sk_interference,
    fun_select_trasha_interference, fun_select_3d_interference,
};
use super::pure2d::fun_reduce_pure2d_i;
use super::reducer_edge::{fun_reducer_edge, fun_reducer_sd_edge, fun_unkeep_evi_on_gb1};
use super::reducer_edge3d::fun_reducer_edge3d;
use super::{fun_keep_l3d_f, fun_recl_se};
use crate::topopebrep::tool::{fun_tool_proj_pon_e, fun_tool_value};

/// `FUN_ProcessEdgeInterferences(EIX, K, G, HDS, LI)` (`cxx:827-958`).
///
/// The reduction proper of one `(Kind, G)` group of the interferences attached
/// to the edge `EIX`.
pub fn fun_process_edge_interferences(
    eix: i32,
    _k: Kind,
    _g: i32,
    bds: &mut DataStructure,
    li: &mut InterferenceList,
) -> Result<(), String> {
    let e = bds
        .shape(eix, true)
        .ok_or_else(|| "FUN_ProcessEdgeInterferences: BDS.Shape(EIX)".to_string())?;

    // LI -> (lF + lFE) + lE + [LI]
    // lF  = {interference on edge <EIX> :  (T(face),G=POINT/VERTEX,S)
    // lFE = {interference on edge <EIX> :  (T(face),G=POINT/VERTEX,S=EDGE)
    // lE  = {interference on edge <EIX> :  (T(edge),G=POINT/VERTEX,S)
    let mut lf: InterferenceList = Vec::new();
    fun_select_trasha_interference(li, ShapeType::Face, &mut lf);
    let mut lfe: InterferenceList = Vec::new();
    fun_select_sk_interference(&mut lf, Kind::Edge, &mut lfe);
    let mut le: InterferenceList = Vec::new();
    fun_select_trasha_interference(li, ShapeType::Edge, &mut le);

    // xpu210798 : CTS21216 (e7,G=P3,S=E9,FTRASHA=f8,f10)
    //   EIX is section edge
    //   2dI1=(REVERSED(ftrasha1),G,ES), 2dI2=(FORWARD(ftrasha2),G,ES)
    //   for the compute of splitON(EIX), reduce 2dI1+2dI2 -> 2dR=(IN(ftrasha1,ftrasha2)
    let isse = bds.is_section_edge(&e, true);
    if isse {
        let mut li2dfe: InterferenceList = Vec::new();
        let mut lri2dfe: InterferenceList = Vec::new();
        // lF  -> lF
        // lFE -> lI2dFE + [lFE] / lI2dFE={FEI=(T(FTRASHA),G,SE) : no FFI=(T(FTRASHA),G,FTRASHA)}
        fun_select_pure2d_i(&lf, &mut lfe, &mut li2dfe);
        fun_reduce_pure2d_i(&mut li2dfe, &mut lri2dfe); // lI2dFE -> lRI2dFE + lI2dFE
        lfe.extend(li2dfe);
        lfe.extend(lri2dfe);
    }

    // I -> 3dI +2dI + 1dI:
    // -------------------
    // lFE -> l3dFE [+l3dFEresi] +l2dFE [+lFE (+lFEresi)]
    //    : l3dFE={I3d } to reduce; for I3dFE=(T(F),G,SE) there is I3dF=(T(F),G,F)
    //      l2dFE={I2d} to reduce
    //      lFEresi to remove
    // lF  -> l3dF [+lF] : l3dF to remove
    //
    // lE  -> l1dE +l1dEsd [+lE}
    //    : l1dE={I1d=(T(ETRA),G,ETRA)/same ETRA} to reduce
    //      l2dEsd={I1dsd=(T(ETRA),G,ETRA)/ E sdm ETRA}
    let mut lfe_resi: InterferenceList = Vec::new();
    let mut l3d_fe: InterferenceList = Vec::new();
    let mut l3d_f: InterferenceList = Vec::new();
    let mut l3d_fe_resi: InterferenceList = Vec::new();
    let mut l2d_fe: InterferenceList = Vec::new();
    fun_select_3d_interference(
        eix,
        bds,
        &mut lf,
        &mut l3d_f,
        &mut lfe,
        &mut lfe_resi,
        &mut l3d_fe,
        &mut l3d_fe_resi,
        &mut l2d_fe,
    )?;
    let mut l1d_e: InterferenceList = Vec::new();
    fun_select_2d_i(eix, bds, ShapeType::Edge, &mut le, &mut l1d_e)?;
    let mut l1d_esd: InterferenceList = Vec::new();
    fun_select_1d_i(eix, bds, &mut le, &mut l1d_esd)?;

    // reducer3d :
    // ----------
    // l3dFE -> lR3dFE [+l3dFE (non reduced 3dI)]
    let mut lr3d_fe: InterferenceList = Vec::new();
    fun_reducer_edge3d(eix, bds, &mut l3d_fe, &mut lr3d_fe)?;

    //  FUN_unkeepEVIonGb1(BDS,EIX,l1dE);  // filter :
    //  FUN_unkeepEVIonGb1(BDS,EIX,l2dFE); // filter :

    // no reduction on G : we keep l3dF(G)
    let mut l3d_f_keep: InterferenceList = Vec::new();
    fun_keep_l3d_f(eix, &l3d_f, &lr3d_fe, &mut l3d_f_keep);
    lf.extend(l3d_f_keep);

    // reducer2d :
    // ----------
    // l2dFE -> LR2dFE + l2dFE (non reduced 2dI)
    let mut lr2d_fe: InterferenceList = Vec::new();
    fun_reducer_edge(eix, bds, &mut l2d_fe, &mut lr2d_fe)?;

    // reducer1d :
    // ----------
    // xpu210498 : reduce interferences I1(T1(esd1),VG,esd1), I2(T2(esd2),VG,esd2)
    // if EIX sdm {esd1,esd2}
    // - t2_2 e32 has I1(in/ou(e32),v30,e20) && I2(ou/in(e20),v30,e20) -
    let mut lr1d_esd: InterferenceList = Vec::new();
    fun_reducer_sd_edge(eix, bds, &mut l1d_esd, &mut lr1d_esd)?;

    // xpu190298 : edge <EIX> same domain with closed support edge <SE>
    // at same G = vertex <VG>  = closing vertex of <SE>
    //  (I1 = (OU/IN(SE),VG,SE) && I2 = (IN/OU(SE),VG,SE)) -> Ir = (IN/IN(SE),VG,SE)
    let mut lr1d_closed_se: InterferenceList = Vec::new();
    fun_recl_se(eix, bds, &mut l1d_e, &mut lr1d_closed_se)?;

    // l1dE  -> LR1dE + l1dE (non reduced 2dI)
    let mut lr1d_e: InterferenceList = Vec::new();
    fun_reducer_edge(eix, bds, &mut l1d_e, &mut lr1d_e)?;

    //  attached to edge <EIX>,
    //  at same G : I  =(T,  G,S)   gives less information than
    //              Ir =(Tr,G,Sr) -reduced interference- with valid Tr
    //  -> unkeep I.
    // using reduced I : lR3dFE, LR2dFE, LR1dE
    let mut lri: InterferenceList = Vec::new();
    lri.extend(lr1d_esd);
    lri.extend(lr1d_e);
    lri.extend(lr1d_closed_se);
    lri.extend(lr2d_fe);
    lri.extend(lr3d_fe);

    lf.extend(lfe);
    lf.extend(l3d_fe);
    lf.extend(l2d_fe); // lF += LFE + l3dFE + l2dFE
    le.extend(l1d_e); // lE += l1dE
    le.extend(l1d_esd); // lE += l1dEsd xpu210498

    // xpu : 23-01-98 : attached to edge <EIX>,
    //  at same G : I  =(T, G, S) gives less information than
    //              Ir =(Tr,G,Sr) -reduced interference-
    //  -> I is added after Ir to help the edge builder.
    li.clear();
    li.extend(lri);
    li.extend(le);
    li.extend(lf);

    // xpu260698 : cto902A5, spOU(e6)
    //             cto801G1, spON(se26)
    if isse {
        fun_unkeep_evi_on_gb1(bds, eix, li); // filter :
    }
    Ok(())
}

/// `CheckInterferenceIsValid(I, theEdge, theSupportEdge, theVertex)`
/// (`cxx:1163-1221`).
///
/// The parameter of the interference on `theEdge` is projected on
/// `theSupportEdge`; the interference is valid when the two points are within
/// the largest of the three tolerances.
pub fn check_interference_is_valid(
    i: &Interference,
    the_edge: &TopoShape,
    the_support_edge: &TopoShape,
    the_vertex: Option<&TopoShape>,
) -> Result<bool, String> {
    let mut pref = 0.0;
    let mut ok = false;

    if let InterferenceData::CurvePoint { param } = i.data() {
        pref = *param;
        ok = true;
    }
    if !ok {
        if let InterferenceData::EdgeVertex(d) = i.data() {
            pref = d.param;
            ok = true;
        }
    }
    if !ok {
        return Ok(ok);
    }

    let Some(p3d1) = fun_tool_value(pref, the_edge) else {
        return Ok(false);
    };
    let Some((paron_support_e, _dist)) = fun_tool_proj_pon_e(&p3d1, the_support_edge) else {
        return Ok(false);
    };
    let Some(p3d2) = fun_tool_value(paron_support_e, the_support_edge) else {
        return Ok(false);
    };

    let tol_edge = BRepTool::edge_tolerance(&Edge(the_edge.clone()));
    let tol_support = BRepTool::edge_tolerance(&Edge(the_support_edge.clone()));
    let mut tolerance = if tol_edge > tol_support {
        tol_edge
    } else {
        tol_support
    };
    if let Some(v) = the_vertex {
        let tol_v = BRepTool::vertex_tolerance(&Vertex(v.clone()));
        tolerance = if tol_v > tolerance { tol_v } else { tolerance };
    }
    if p3d1.distance(&p3d2) > tolerance {
        ok = false;
    }
    Ok(ok)
}
