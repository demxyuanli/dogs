//! `TopOpeBRep_vpr.cxx` — the transition of a VPoint on a restriction edge, and
//! the interference-list selectors it needs.
//!
//! Sources (OCCT 8.0.0, `ModelingAlgorithms/TKBool/TopOpeBRep`):
//! - `FUN_LineRestF` (`TopOpeBRep_vpr.cxx:471-513`)
//! - `FUN_newtransEdge` (`vpr.cxx:514-682`)
//! - `FUN_ScanInterfList` (`vpr.cxx:684-708`)
//! - `FUN_selectTRAISHAinterference` (`vpr.cxx:710-738`)
//! - `FUN_selectGinterference` (`vpr.cxx:740-757`)
//! - `FUN_sameGsameS` (`vpr.cxx:759-780`)
//!
//! `FUN_newtransEdge` is the only caller of `TopOpeBRepTool_makeTransition` on
//! the `<edge>` / `<OOface>` side (`MKT` below); the transition of the
//! restriction edge against the same face is computed in
//! [`crate::topopebrep::faces_filler::vpclo`].

use occt_core::precision::ANGULAR;

use crate::abs::{Orientation, ShapeType};
use crate::shape::TopoShape;
use crate::topopebrep::ds::export::{fds_getupperlower, fds_loi_inf_sup};
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::ds::interference_iterator::InterferenceList;
use crate::topopebrep::ds::kind::Kind;
use crate::topopebrep::ds::point::Point;
use crate::topopebrep::ds::point_iterator::PointIterator;
use crate::topopebrep::ds::transition::Transition;
use crate::topopebrep::interf::v_point_inter::VPointInter;
use crate::topopebrep::line_inter::LineInter;
use crate::topopebrep::tool::make_transition::MakeTransition;
use crate::topopebrep::tool::topology::fun_tool_etg_f;
use crate::topopebrep::type_line_curve::TypeLineCurve;

/// `FUN_LineRestF(F, L, ERL, ER)` (`vpr.cxx:471-513`): `L` is on the
/// restriction `ER` of `F`, and `ER` is the only one of `ERL` it lies on.
pub fn fun_line_rest_f(
    f: &TopoShape,
    l: &mut LineInter,
    erl: &[TopoShape],
    global_tol_ff: f64,
) -> Result<Option<TopoShape>, String> {
    // `TopExp::MapShapes(F, TopAbs_EDGE, mapE)`
    let mut map_e: Vec<TopoShape> = Vec::new();
    {
        let mut ex = crate::topexp::Explorer::new(f, crate::abs::ShapeType::Edge);
        while ex.more() {
            map_e.push(ex.current().clone());
            ex.next();
        }
    }

    let mut erl_on_f: Vec<TopoShape> = Vec::new();
    for e in erl {
        if map_e.iter().any(|e0| crate::topo_tools::is_same(e0, e)) {
            erl_on_f.push(e.clone());
        }
    }

    let mut erl_on_f_on_l: Vec<TopoShape> = Vec::new();
    for e in erl_on_f.iter() {
        let el = [e.clone()];
        let is_on_l = super::process_gr::lsame_domain_erl(l, &el, global_tol_ff)?;
        if is_on_l {
            erl_on_f_on_l.push(e.clone());
        }
    }

    // `<L>` is on at most one edge restriction.
    if erl_on_f_on_l.len() != 1 {
        return Ok(None);
    }
    Ok(erl_on_f_on_l.first().cloned())
}

/// `FUN_ScanInterfList(PDS, HDS, loI, loIfound)` (`vpr.cxx:684-708`): every
/// interference of `loI` whose geometry is the point `PDS`.
pub fn fun_scan_interf_list(
    pds: &Point,
    hds: &HDataStructure,
    lo_i: &InterferenceList,
    lo_i_found: &mut InterferenceList,
) -> Result<(), String> {
    let mut it = PointIterator::new(lo_i.clone());
    while it.more() {
        let found = hds.scan_interf_list(&mut it, pds)?;
        if found {
            lo_i_found.push(it.value().clone());
            if it.more() {
                it.next();
            }
        } else {
            return Ok(());
        }
    }
    Ok(())
}

/// `FUN_selectTRAISHAinterference(lI, ITRASHA, lITRAonISHA)` (`vpr.cxx:710-738`).
///
/// `true` when at least one interference has `I->Transition().Index() ==
/// ITRASHA`; the `EXTERNAL` ones are skipped (`cxx:727-730`, `xpu030998`).
pub fn fun_select_traisha_interference(
    l_i: &InterferenceList,
    i_trasha: i32,
    l_it_ra_on_isha: &mut InterferenceList,
) -> Result<bool, String> {
    l_it_ra_on_isha.clear();
    for h in l_i {
        let (ori, index) = {
            let i = h.read().expect("interference lock");
            let t = i.transition();
            (
                t.orientation(crate::abs::State::In, ShapeType::Face),
                t.index(),
            )
        };
        if ori == Orientation::External {
            continue;
        }
        let Ok(index) = index else {
            // `T.Index()` raises when the two sides disagree; OCCT would throw.
            return Err("FUN_selectTRAISHAinterference : Transition::Index".to_string());
        };
        if index == i_trasha {
            l_it_ra_on_isha.push(h.clone());
        }
    }
    Ok(!l_it_ra_on_isha.is_empty())
}

/// `FUN_selectGinterference(lI, G, lIonG)` (`vpr.cxx:740-757`).
pub fn fun_select_g_interference(
    l_i: &InterferenceList,
    g: i32,
    l_i_on_g: &mut InterferenceList,
) -> bool {
    l_i_on_g.clear();
    for h in l_i {
        if h.read().expect("interference lock").geometry() == g {
            l_i_on_g.push(h.clone());
        }
    }
    !l_i_on_g.is_empty()
}

/// `FUN_sameGsameS(loI, G, S, loIfound)` (`vpr.cxx:759-780`): the
/// interferences of `loI` carrying geometry `G` and support `S`.
pub fn fun_same_g_same_s(
    lo_i: &InterferenceList,
    g: i32,
    s: i32,
    lo_i_found: &mut InterferenceList,
) -> bool {
    lo_i_found.clear();
    let mut pi = PointIterator::new(lo_i.clone());
    while pi.more() {
        let (gepi, sepi) = {
            let i = pi.current();
            (i.geometry(), i.support())
        };
        if gepi == g && sepi == s {
            lo_i_found.push(pi.value().clone());
        }
        pi.next();
    }
    !lo_i_found.is_empty()
}

/// `FUN_newtransEdge(HDS, FF, L, Lonrest, VP, PVKind, PVIndex, OOShapeIndex,
/// edge, ERL, T)` (`vpr.cxx:514-682`).
///
/// `FF` only provides `Face(OOShapeIndex)`, which the caller already has as
/// `oo_face`; `FF`'s tracing fields are debug-only.
#[allow(clippy::too_many_arguments)]
pub fn fun_newtrans_edge(
    hds: &mut HDataStructure,
    facint: &mut crate::topopebrep::faces_intersector::FacesIntersector,
    line_index: i32,
    lonrest: bool,
    vp: &VPointInter,
    pv_kind: Kind,
    pv_index: i32,
    oo_shape_index: i32,
    oo_face: &TopoShape,
    edge: &TopoShape,
    erl: &[TopoShape],
    t: &mut Transition,
) -> Result<bool, String> {
    t.before(crate::abs::State::Unknown, ShapeType::Face);
    t.after(crate::abs::State::Unknown, ShapeType::Face);

    // `cxx:531-538`: a face INTERNAL/EXTERNAL in its solid keeps its
    // orientation as the transition.
    let o_fie = oo_face.orientation();
    if o_fie == Orientation::Internal || o_fie == Orientation::External {
        t.set_orientation(o_fie);
        return Ok(true);
    }

    // `cxx:545-549`: the parameter of `VP` on `<edge>`.
    let mut paredge = 0.0;
    if !vp.paron_e(&crate::shape::Edge(edge.clone()), &mut paredge) {
        return Ok(false);
    }

    // `cxx:552-576`: the interval of `<edge>` around `paredge`.
    let (par1, par2) = if hds.ds().has_shape(edge, true) {
        if pv_index == 0 {
            let idx = hds.ds().shape_index(edge, true);
            let ds_edge = hds
                .ds()
                .shape(idx, true)
                .ok_or_else(|| "FUN_newtransEdge : no DS edge".to_string())?;
            fds_getupperlower(hds, &ds_edge, paredge)?
        } else {
            let loi = hds.ds().shape_interferences(edge, true).to_vec();
            let (pbef, paft, _isonper) = fds_loi_inf_sup(edge, paredge, pv_kind, pv_index, &loi)?;
            (pbef, paft)
        }
    } else {
        crate::topopebrep::tool::fun_tool_bounds(edge)
    };

    let uv = vp.surface_parameters(oo_shape_index)?;

    // `<Tr>` relative to the 3d `<OOface>` matter; the 2d `<OOface>` is only
    // taken into account when `<edge>` is normal to `<OOface>`.
    let tola = ANGULAR * 1.0e+4; // dealing with tolerances
    let etg_oof = fun_tool_etg_f(paredge, edge, &uv, oo_face, tola);
    let in_erl = super::vpr::fun_in_los(edge, erl);
    let isse = hds.ds().is_section_edge(edge, true);
    let rest = in_erl || isse;
    let mut interf2d = etg_oof && lonrest && rest;
    let interf3dtg = etg_oof && rest && !interf2d; // xpu260898 :cto902D6,(e15,p3,f9)

    let factor = 1.0e-2;
    let mut mkt = MakeTransition::new();
    if !mkt.initialize(edge, par1, par2, paredge, oo_face, &uv, factor) {
        return Ok(false);
    }
    let is_t2d = mkt.is_t2d();
    interf2d = interf2d && is_t2d;

    if interf2d {
        // `<tgLine>`:
        let on_oo_face;
        let ooer;
        let typ_l = facint.change_line(line_index).type_line_curve();
        if typ_l == TypeLineCurve::Restriction {
            on_oo_face = true;
            ooer = facint
                .change_line(line_index)
                .arc()
                .ok_or_else(|| "FUN_newtransEdge : no restriction arc".to_string())?
                .clone();
        } else {
            let tol = facint.global_tol_ff();
            let found =
                fun_line_rest_f(oo_face, facint.change_line(line_index), erl, tol)?;
            on_oo_face = found.is_some();
            ooer = found.unwrap_or_else(|| TopoShape::new(crate::abs::ShapeType::Edge));
        }
        if !on_oo_face {
            return Ok(false);
        }

        let mut oo_par = 0.0;
        let mut ok = vp.paron_e(&crate::shape::Edge(ooer.clone()), &mut oo_par);
        if !ok {
            if let Some(p) = crate::topopebrep::tool::project_curve::fun_tool_par_e(
                edge, paredge, &ooer,
            ) {
                oo_par = p;
                ok = true;
            }
        }
        if !ok {
            return Ok(false);
        }

        // xpu051098 : cto900L4 (edge18,OOface5)
        if !mkt.set_rest(&ooer, oo_par) {
            return Ok(false);
        }
    } else if interf3dtg {
        let absindex = vp.shape_index(); // 0,1,2,3
        let on2edges = absindex == 3;
        let has_on_edge = vp.state(oo_shape_index)? == crate::abs::State::On;
        let has_oo_edge = if on2edges { true } else { has_on_edge };

        if has_oo_edge {
            let (ooedge, oo_par) = if on2edges {
                (
                    vp.edge(oo_shape_index)
                        .ok_or_else(|| "FUN_newtransEdge : no OO edge".to_string())?
                        .clone(),
                    vp.edge_parameter(oo_shape_index),
                )
            } else {
                (
                    vp.edge_on(oo_shape_index)?.clone(),
                    vp.edge_on_parameter(oo_shape_index)?,
                )
            };
            if !mkt.set_rest(&ooedge, oo_par) {
                return Ok(false);
            }
        }
    }

    let mut stb = crate::abs::State::Unknown;
    let mut sta = crate::abs::State::Unknown;
    if !mkt.mk_ton_e(&mut stb, &mut sta) {
        return Ok(false);
    }
    t.before(stb, ShapeType::Face);
    t.after(sta, ShapeType::Face);
    Ok(true)
}