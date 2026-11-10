use crate::shape::TopoShape;

use super::super::data_structure::DataStructure;
use super::super::h_data_structure::HDataStructure;
use super::super::interference_iterator::InterferenceList;
use super::super::kind::Kind;
use super::super::tki::TKI;
use super::super::transition::Transition;

use super::*;

/// `FDS_LOIinfsup(BDS, E, pE, KDS, GDS, LOI, pbef, paft, isonboundper)`
/// (`cxx:3427-3566`): the parameters `pbef < pE < paft` bracketing `pE` on
/// `E`, taken from the interferences of `LOI` on a vertex/point other than
/// `(KDS, GDS)` whose transition is `(FACE, FACE)`.
///
/// The C++ tests `n == 0` *twice* (`cxx:3444`, `:3503`): the second one is
/// after the grouping loop and skips to the same answer. Both return `true`
/// (`cxx:3447`, `:3505`), so the port's `Err` only covers the
/// `FDS_Parameter` throw (`cxx:3508`, `:3555`).
pub fn fds_loi_inf_sup(
    e: &TopoShape,
    p_e: f64,
    k_ds: Kind,
    g_ds: i32,
    loi: &InterferenceList,
) -> Result<(f64, f64, bool), String> {
    let (f, l) = crate::topopebrep::tool::fun_tool_bounds(e);
    let mut pbef = f;
    let mut paft = l;
    if loi.is_empty() {
        return Ok((pbef, paft, false));
    }

    let mut loi_sans_gds: InterferenceList = Vec::new();
    let mut tki = TKI::new();
    tki.fill_on_geometry(loi);
    tki.init();
    while tki.more() {
        let group = tki.value().map(|(k, g, l)| (k, g, l.clone()));
        tki.next();
        let Some((k, g, list)) = group else {
            continue;
        };
        let pv = k == Kind::Point || k == Kind::Vertex;
        if !pv {
            continue;
        }
        let mk = k == k_ds;
        let mg = g == g_ds;
        if mk && mg {
            continue;
        }
        for h in list.iter() {
            let keep = {
                let i = h.read().expect("interference lock");
                let (tsb1, _isb1, tsa1, _isa1, _gt1, _g1, _st1, _s1) = fds_idata(&i);
                if tsb1 != crate::abs::ShapeType::Face {
                    false
                } else if tsa1 != crate::abs::ShapeType::Face {
                    false
                } else {
                    matches!(
                        i.data(),
                        super::super::interference::InterferenceData::CurvePoint { .. }
                            | super::super::interference::InterferenceData::EdgeVertex(_)
                    )
                }
            };
            if keep {
                loi_sans_gds.push(h.clone());
                break;
            }
        }
    }

    if loi_sans_gds.is_empty() {
        return Ok((pbef, paft, false));
    }

    let v = crate::topopebrep::tool::closed_e(e);
    let e_closed = v.is_some();
    let tole = crate::brep_tool::BRepTool::edge_tolerance(&crate::shape::Edge(e.clone()));
    let mut tol = occt_core::precision::Precision::parametric(tole, 1.0e2);
    let mut isonboundper = false;
    if e_closed {
        let v = v.expect("closed_e returned a vertex");
        let tolv = crate::brep_tool::BRepTool::vertex_tolerance(&crate::shape::Vertex(v));
        let tolv = occt_core::precision::Precision::parametric(tolv, 1.0e2);
        if tolv > tol {
            tol = tolv;
        }
        let p_e_is_ef = (p_e - f).abs() <= tol;
        let p_e_is_el = (p_e - l).abs() <= tol;
        isonboundper = p_e_is_ef || p_e_is_el;
    }

    for h in loi_sans_gds.iter() {
        let p = fds_parameter_or_fail(&h.read().expect("interference lock"))?;
        if isonboundper {
            // pbef > paft
            if p > pbef {
                pbef = p;
            }
            if p < paft {
                paft = p;
            }
        } else {
            if p > pbef && p < p_e {
                pbef = p;
            }
            if p < paft && p > p_e {
                paft = p;
            }
        }
    }
    Ok((pbef, paft, isonboundper))
}

/// `FDS_parbefaft(BDS, E, pE, pbef, paft, isonboundper, p1, p2)`
/// (`cxx:3572-3596`): the two interior parameters `p1` in `[pbef, pE[` and
/// `p2` in `]pE, paft]` with the fixed ratio `t = 0.3178947713`.
pub fn fds_parbefaft(e: &TopoShape, p_e: f64, pbef: f64, paft: f64, isonboundper: bool) -> (f64, f64) {
    let t = 0.3178947713;
    let (f, l) = crate::topopebrep::tool::fun_tool_bounds(e);
    if isonboundper {
        ((1.0 - t) * pbef + t * l, (1.0 - t) * f + t * paft)
    } else {
        ((1.0 - t) * pbef + t * p_e, (1.0 - t) * p_e + t * paft)
    }
}

/// `FDS_stateEwithF2d(BDS, E, pE, KDS, GDS, F1, TrmemeS)`
/// (`cxx:3598-3643`): the 2d transition of `E` on `F1` at `pE`, i.e. the
/// states of the two bracketing points.
///
/// `ON` is mapped to `IN` (`cxx:3629-3636`, "xpu190898 : cto014I2").
pub fn fds_state_e_with_f2d(
    bds: &DataStructure,
    e: &TopoShape,
    p_e: f64,
    k_ds: Kind,
    g_ds: i32,
    f1: &TopoShape,
    tr_meme_s: &mut Transition,
) -> Result<bool, String> {
    let loi = bds.shape_interferences(e, true).to_vec();
    let (pbef, paft, isonper) = fds_loi_inf_sup(e, p_e, k_ds, g_ds, &loi)?;
    let (t1, t2) = fds_parbefaft(e, p_e, pbef, paft, isonper);
    let (Some(p1), Some(p2)) = (
        crate::topopebrep::tool::fun_tool_value(t1, e),
        crate::topopebrep::tool::fun_tool_value(t2, e),
    ) else {
        return Ok(false);
    };

    let mut psc = crate::topopebrep::tool::fsc_get_psc(f1);
    let mut sta1 = crate::topopebrep::tool::fsc_state_pon_face(&p1, f1, &mut psc);
    let mut sta2 = crate::topopebrep::tool::fsc_state_pon_face(&p2, f1, &mut psc);

    // xpu190898 : cto014I2 (e5,fTRASHA14,vG7,eS10)
    if sta1 == crate::abs::State::On {
        sta1 = crate::abs::State::In;
    }
    if sta2 == crate::abs::State::On {
        sta2 = crate::abs::State::In;
    }

    tr_meme_s.before(sta1, crate::abs::ShapeType::Face);
    tr_meme_s.after(sta2, crate::abs::ShapeType::Face);
    Ok(true)
}

/// `FDS_getupperlower(HDS, edgeIndex, paredge, p1, p2)`
/// (`TopOpeBRepDS_EXPORT.cxx:617-643`): `p1`/`p2` start at the edge bounds and
/// are narrowed to the point-interference parameters bracketing `paredge`.
///
/// The C++ `EdgePoints(E)` walk yields only point/vertex interferences, whose
/// `Parameter()` never raises; the port still propagates the `Err` of the four
/// other interference types in case a list holds one.
pub fn fds_getupperlower(
    hds: &HDataStructure,
    e: &TopoShape,
    paredge: f64,
) -> Result<(f64, f64), String> {
    let (mut p1, mut p2) = crate::topopebrep::tool::fun_tool_bounds(e);
    let mut it = hds.edge_points(e);
    while it.more() {
        let par = it.parameter()?;
        it.next();
        let parsup1 = par > p1;
        let parinfe = par < paredge;
        let parinf2 = par < p2;
        let parsupe = par > paredge;
        if parsup1 && parinfe {
            p1 = par;
        }
        if parinf2 && parsupe {
            p2 = par;
        }
    }
    Ok((p1, p2))
}
