//! The 1d/2d edge reducers of `TopOpeBRepDS_EIR`.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_EIR.cxx` (OCCT 8.0.0):
//! - `FUN_ReducerEdge` (`:377-387`)
//! - `FUN_ReducerSDEdge` (`:389-505`)
//! - `FUN_unkeepEVIonGb1` (`:651-689`)
//!
//! `FUN_ReducerSDEdge` merges the two `(OU/IN(esd),VG,esd)` interferences of a
//! same-domain support edge into one `(IN/IN, VG, EDGE)` interference
//! (`cxx:465-494`). On the reduced arm `I1` is *replaced* by the new
//! interference: it is neither appended to `newLI` nor removed from `loi`,
//! which is discarded with the `TKI` (`cxx:491`).

use crate::abs::{Orientation, ShapeType, State};

use super::super::data_structure::DataStructure;
use super::super::export::fun_ds_sdm;
use super::super::interference::{Interference, InterferenceData};
use super::super::interference_iterator::InterferenceList;
use super::super::interference_tool::make_handle;
use super::super::kind::{Config, Kind};
use super::super::repvg::fds_repvg;
use super::super::tki::TKI;
use super::super::transition::Transition;
use super::fds_reduce_onface_interferences;
use crate::topopebrep::tool::fun_tool_orient_vin_e;

/// `M_FORWARD(O)` (`cxx:36`).
fn m_forward(o: Orientation) -> bool {
    o == Orientation::Forward
}

/// `M_REVERSED(O)` (`cxx:37`).
fn m_reversed(o: Orientation) -> bool {
    o == Orientation::Reversed
}

/// `FUN_ReducerEdge(SIX, BDS, LI, reducedLI)` (`cxx:377-387`).
///
/// The vertex pass, the `(ON, FACE)` filter and the point pass, in that order;
/// `reducedLI` accumulates across the two `FDS_repvg` calls.
pub fn fun_reducer_edge(
    six: i32,
    bds: &DataStructure,
    li: &mut InterferenceList,
    reduced_li: &mut InterferenceList,
) -> Result<(), String> {
    fds_repvg(bds, six, Kind::Vertex, li, reduced_li)?;
    fds_reduce_onface_interferences(li);
    fds_repvg(bds, six, Kind::Point, li, reduced_li)?;
    Ok(())
}

/// `FUN_ReducerSDEdge(SIX, BDS, LI, reducedLI)` (`cxx:389-505`).
pub fn fun_reducer_sd_edge(
    six: i32,
    bds: &DataStructure,
    li: &mut InterferenceList,
    reduced_li: &mut InterferenceList,
) -> Result<(), String> {
    reduced_li.clear();
    if li.len() <= 1 {
        return Ok(());
    }

    let mut new_li: InterferenceList = Vec::new();
    let eix = bds
        .shape(six, true)
        .ok_or_else(|| "FUN_ReducerSDEdge: BDS.Shape(SIX)".to_string())?;

    let mut tki = TKI::new();
    tki.fill_on_geometry(li);
    let mut order: Vec<(Kind, i32)> = Vec::new();
    tki.init();
    while tki.more() {
        if let Some((k, g, _)) = tki.value() {
            order.push((k, g));
        }
        tki.next();
    }

    for (k, g) in order {
        let mut loi: InterferenceList = tki.interferences(k, g).to_vec();

        let mut i1 = 0;
        while i1 < loi.len() {
            let h1 = loi[i1].clone();
            let (o1, ib1) = {
                let i = h1.read().expect("interference lock");
                let t1 = i.transition();
                (t1.orientation(State::In, ShapeType::Face), t1.index()?)
            };

            let cond1 = match bds.shape(ib1, true) {
                Some(s) => fun_ds_sdm(bds, &eix, &s),
                None => false,
            };
            if !cond1 {
                new_li.push(h1);
                i1 += 1;
                continue;
            }

            if i1 + 1 >= loi.len() {
                break;
            }

            let mut complex1d = false;
            let mut t = Transition::new(
                State::In,
                State::In,
                ShapeType::Edge,
                ShapeType::Edge,
            );

            let mut i2 = i1 + 1;
            while i2 < loi.len() {
                let h2 = loi[i2].clone();
                let (o2, ib2) = {
                    let i = h2.read().expect("interference lock");
                    let t2 = i.transition();
                    (t2.orientation(State::In, ShapeType::Face), t2.index()?)
                };

                let cond2 = match bds.shape(ib2, true) {
                    Some(s) => fun_ds_sdm(bds, &eix, &s),
                    None => false,
                };
                if !cond2 {
                    new_li.push(h2);
                    i2 += 1;
                    continue;
                }

                complex1d = (m_forward(o1) && m_reversed(o2))
                    || (m_forward(o2) && m_reversed(o1));
                if !complex1d {
                    new_li.push(h2);
                    i2 += 1;
                    continue;
                }

                let ib = if m_reversed(o1) { ib1 } else { ib2 };
                t.set_index_before(ib);
                let ia = if m_reversed(o1) { ib2 } else { ib1 };
                t.set_index_after(ia);

                loi.remove(i2);
                break; // no more than 2 interferences on sdmTRASHA at same G
            } // it2

            if complex1d {
                let (gbound, param) = {
                    let i = h1.read().expect("interference lock");
                    match i.as_edge_vertex() {
                        Some(d) => (d.shape_shape.g_bound, d.param),
                        None => {
                            return Err(
                                "FUN_ReducerSDEdge: down_cast<EdgeVertexInterference>".to_string()
                            )
                        }
                    }
                };
                let c_eix = bds.same_domain_ori(six);
                let c1 = bds.same_domain_ori(ib1);
                let conf = if c_eix == c1 {
                    Config::SameOriented
                } else {
                    Config::DiffOriented
                };
                let new_i = Interference::edge_vertex(
                    t.clone(),
                    Kind::Edge,
                    ib1,
                    g,
                    gbound,
                    conf,
                    param,
                );
                reduced_li.push(make_handle(new_i));
                i1 += 1;
            } else {
                new_li.push(h1);
                loi.remove(i1);
            }
        } // it1
    } // tki

    li.clear();
    li.extend(new_li);
    Ok(())
}

/// `FUN_unkeepEVIonGb1(BDS, EIX, LI)` (`cxx:651-689`).
///
/// Drops the `EdgeVertexInterference`s whose geometry is a vertex of `EIX`
/// (`FUN_tool_orientVinE != 0`).
pub fn fun_unkeep_evi_on_gb1(bds: &DataStructure, eix: i32, li: &mut InterferenceList) {
    let Some(e) = bds.shape(eix, true) else {
        return;
    };

    let mut i = 0;
    while i < li.len() {
        let (evi, gt, g) = {
            let x = li[i].read().expect("interference lock");
            (
                matches!(x.data(), InterferenceData::EdgeVertex(_)),
                x.geometry_type(),
                x.geometry(),
            )
        };
        if !evi {
            i += 1;
            continue;
        }
        if gt != Kind::Vertex {
            i += 1;
            continue;
        }

        let Some(v) = bds.shape(g, true) else {
            i += 1;
            continue;
        };
        let o = fun_tool_orient_vin_e(&v, &e);
        if o == 0 {
            i += 1;
            continue;
        }

        li.remove(i);
    }
}
