//! `TopOpeBRepDS_ProcessInterferencesTool` — the 3d selection filters.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_ProcessInterferencesTool.cxx`
//! (OCCT 8.0.0):
//! - `FUN_select3dI` (`:531-731`)
//! - `FUN_find3dISEsameISF` (`:733-771`)
//! - `FUN_select3dISEsameISF` (`:774-844`)
//! - `FUN_select3dinterference` (`:1147-1185`)
//!
//! Translation notes:
//! - `sameG` (`cxx:619`) is computed and never read: the C++ assigns it and
//!   never branches on it, so the port keeps the computation as `let _same_g`.
//! - `lFE.Remove(it)` at `cxx:719` removes through an iterator that belongs to
//!   `lFEresi`; `NCollection_List::Remove(Iterator&)` unlinks the node, so the
//!   element leaves `lFEresi`. The port spells that as `lfe_resi.remove(i)`.
//! - `TopoDS::Face(BDS.Shape(I))` raises when the shape is missing; those arms
//!   are `Err` here. `FUN_tool_Eshared` on a null `OOv` (`cxx:645`) answers
//!   `false`, which the port reproduces with `Option<TopoShape>`.

use crate::abs::{ShapeType, State};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, TopoShape};

use super::super::data_structure::DataStructure;
use super::super::export::{fds_data, fds_tdata, fun_ds_get_vsdm, fun_ds_sdm};
use super::super::interference::Interference;
use super::super::interference_iterator::InterferenceList;
use super::super::kind::Kind;
use super::fun_select_2d_i;
use super::fun_select_itrasha_interference;
use crate::topopebrep::tool::fun_tool_eshared;

/// `BDS.Shape(I)` followed by `TopoDS::Face` (`cxx:625-626`).
fn face_at(bds: &DataStructure, i: i32) -> Result<TopoShape, String> {
    let s = bds
        .shape(i, true)
        .ok_or_else(|| "FUN_select3dI: BDS.Shape".to_string())?;
    if s.shape_type() != ShapeType::Face {
        return Err("FUN_select3dI: TopoDS::Face".to_string());
    }
    Ok(s)
}

/// `FUN_select3dI(SIX, BDS, lFE, lFEresi, l3dFE)` (`cxx:531-731`).
///
/// Groups the `(T(face), G, S=edge)` interferences of an edge into the 3d
/// interferences that describe a complex transition:
/// - a pair on same-domain faces moves to `lFEresi` (`cxx:626-633`);
/// - a pair whose support edges are the same, or an edge shared by the two
///   faces with the bound vertex `OOv`, moves to `l3dFE` (`cxx:636-681`);
/// - a same-domain reference edge keeps the pair out (`cxx:661-670`).
///
/// The trailing pass (`cxx:694-729`) drops the `lFEresi` elements whose
/// transition index has no counterpart in `l3dFE`, or re-queues the whole
/// `lFEresi` when nothing was reduced.
pub fn fun_select_3d_i(
    six: i32,
    bds: &mut DataStructure,
    lfe: &mut InterferenceList,
    lfe_resi: &mut InterferenceList,
    l3d_fe: &mut InterferenceList,
) -> Result<i32, String> {
    l3d_fe.clear();
    lfe_resi.clear();
    if lfe.len() <= 1 {
        return Ok(0);
    }

    let e = bds
        .shape(six, true)
        .ok_or_else(|| "FUN_select3dI: BDS.Shape(SIX)".to_string())?;
    if e.shape_type() != ShapeType::Edge {
        return Err("FUN_select3dI: TopoDS::Edge".to_string());
    }
    let rank_e = bds.ancestor_rank_of(&e);
    let mut oov: Option<TopoShape> = None;
    let mut oog: i32 = 0;
    let mut gsta: i32 = 0;

    let mut i1 = 0;
    while i1 < lfe.len() {
        let h1 = lfe[i1].clone();
        let (gt1, g1, st1, s1, _sb1, ib1, _sa1, _ia1) = {
            let intf = h1.read().expect("interference lock");
            let (gt, g, st, s) = fds_data(&intf);
            let (sb, ib, sa, ia) = fds_tdata(&intf.transition());
            (gt, g, st, s, sb, ib, sa, ia)
        };
        let mut s1 = s1;
        let vertex1 = gt1 == Kind::Vertex;

        if i1 + 1 >= lfe.len() {
            break;
        }

        // <Gsta>, <OOv>
        if vertex1 {
            let vg1 = bds
                .shape(g1, true)
                .ok_or_else(|| "FUN_select3dI: BDS.Shape(G1)".to_string())?;
            let rank_vg1 = bds.ancestor_rank_of(&vg1);
            let g1hsd = fun_ds_get_vsdm(bds, g1);
            if let Some(sd) = g1hsd {
                oog = sd;
            }
            if rank_vg1 != rank_e {
                // vG1 not on E
                oov = Some(vg1);
                gsta = if g1hsd.is_some() { 3 } else { 2 };
            } else {
                // vG on E
                if g1hsd.is_some() {
                    oov = bds.shape(oog, true);
                }
                gsta = if g1hsd.is_some() { 3 } else { 1 };
            }
        }

        let mut complex3d = false;
        let mut i2 = i1 + 1;
        while i2 < lfe.len() {
            let h2 = lfe[i2].clone();
            let (gt2, g2, st2, s2i, _sb2, ib2, _sa2, _ia2) = {
                let intf = h2.read().expect("interference lock");
                let (gt, g, st, s) = fds_data(&intf);
                let (sb, ib, sa, ia) = fds_tdata(&intf.transition());
                (gt, g, st, s, sb, ib, sa, ia)
            };
            let mut s2 = s2i;

            if gt2 != gt1 {
                break;
            }
            // `cxx:619`: assigned, never read.
            let _same_g = g2 == g1 || g2 == oog;

            let f1 = face_at(bds, ib1)?;
            let f2 = face_at(bds, ib2)?;

            // same domain faces, -> 2d interferences
            let sdm_ft = ib1 == ib2 || fun_ds_sdm(bds, &f1, &f2);
            if sdm_ft {
                lfe_resi.push(h2);
                lfe.remove(i2);
                continue;
            }

            // a. if (S2 == S1) ->ok
            // b. else :
            //     if (G is POINT on Fi i = 1,2) ->continue
            //     else : look for E1(edge S1) on F1(IB1) shared by F2(IB2)
            //            with bound G.
            let same_s = st2 == st1 && s2 == s1;

            let mut same_s_or_shared = same_s;
            let has_oov = gsta > 1;
            if !same_s_or_shared && has_oov {
                let foundsh = match oov.as_ref() {
                    Some(v) => fun_tool_eshared(v, &f1, &f2),
                    None => None,
                };
                let Some(eshared) = foundsh else {
                    i2 += 1;
                    continue;
                };

                same_s_or_shared = true;
                if !bds.has_shape(&eshared, true) {
                    let oorank = match oov.as_ref() {
                        Some(v) => bds.ancestor_rank_of(v),
                        None => 0,
                    };
                    bds.add_shape_with_ancestor(&eshared, oorank);
                }
                s1 = bds.shape_index(&eshared, true);
                s2 = s1;
                // `cxx:659`: `S2 = S1;` — S2 is redeclared at every `it2`
                // iteration (`cxx:608`) and read only before this point, so the
                // store is dead in OCCT; kept for the transcription.
                let _ = s2;
            }

            if same_s_or_shared {
                // xpu : 09-03-98
                let eix = bds
                    .shape(six, true)
                    .ok_or_else(|| "FUN_select3dI: BDS.Shape(SIX)".to_string())?;
                if let Some(sd) = bds.shape(s1, true) {
                    if fun_ds_sdm(bds, &eix, &sd) {
                        i2 += 1;
                        continue;
                    }
                }
            }

            if same_s_or_shared {
                l3d_fe.push(h2);
                lfe.remove(i2);
                complex3d = true;
            } else {
                i2 += 1;
            }
        } // it2

        if complex3d {
            l3d_fe.push(h1);
            lfe.remove(i1);
        } else {
            i1 += 1;
        }
    } // it1

    // xpu170898 : cto009H1 (I1(T1(F),G,S1),I1'(T2(F),G,S2), I2)
    //             (I1,I2) -> to reduce => I1' to delete
    let n3d = l3d_fe.len() as i32;
    if lfe_resi.is_empty() {
        return Ok(n3d);
    }

    if n3d != 0 {
        let mut i = 0;
        while i < lfe_resi.len() {
            let tra = lfe_resi[i]
                .read()
                .expect("interference lock")
                .transition()
                .index()?;
            let mut lcopy = l3d_fe.clone();
            let mut lfound: InterferenceList = Vec::new();
            let nfound = fun_select_itrasha_interference(&mut lcopy, tra, &mut lfound);
            if nfound == 0 {
                lfe_resi.remove(i);
            } else {
                i += 1;
            }
        }
    } else {
        lfe.extend(lfe_resi.iter().cloned());
    }
    Ok(n3d)
}

/// `FUN_find3dISEsameISF(I1, it2)` (`cxx:733-771`): the first interference of
/// `l` that is `(T0, G0, S')` with the same transition orientation as `I1`.
///
/// The C++ leaves `it2` on the match; the port returns its index.
fn fun_find_3d_is_esame_isf(i1: &Interference, l: &InterferenceList) -> Option<usize> {
    let (o1, gt1, g1) = (
        i1.transition().orientation(State::In, ShapeType::Face),
        i1.geometry_type(),
        i1.geometry(),
    );
    let (sb1, ib1, _sa1, _ia1) = fds_tdata(&i1.transition());

    for (k, h) in l.iter().enumerate() {
        let i2 = h.read().expect("interference lock");
        let o2 = i2.transition().orientation(State::In, ShapeType::Face);
        let (gt2, g2, _st2, _s2) = fds_data(&i2);
        let (sb2, ib2, _sa2, _ia2) = fds_tdata(&i2.transition());
        let same_t = sb1 == sb2 && ib1 == ib2 && o1 == o2;
        let same_g = gt1 == gt2 && g1 == g2;
        if same_t && same_g {
            return Some(k);
        }
    }
    None
}

/// `FUN_select3dISEsameISF(lFE, l3dFE, l3dFEresi, lF, l3dF)` (`cxx:774-844`).
///
/// Matches every 3d interference on an edge with its `(T0, G0, S=face)`
/// counterpart: the counterpart moves from `lF` to `l3dF`. When no counterpart
/// is left in `lF` but one already went to `l3dF`, the edge interference is a
/// duplicate and moves to `l3dFEresi`; otherwise it goes back to `lFE`.
pub fn fun_select_3d_is_esame_isf(
    lfe: &mut InterferenceList,
    l3d_fe: &mut InterferenceList,
    l3d_fe_resi: &mut InterferenceList,
    lf: &mut InterferenceList,
    l3d_f: &mut InterferenceList,
) -> i32 {
    let mut i = 0;
    while i < l3d_fe.len() {
        let h1 = l3d_fe[i].clone();
        let found = {
            let intf = h1.read().expect("interference lock");
            fun_find_3d_is_esame_isf(&intf, lf)
        };

        if let Some(pos) = found {
            let moved = lf.remove(pos);
            l3d_f.push(moved);
            i += 1;
        } else {
            // xpu : 16-02-98 :
            // IFE1 = (T0,G0,Sedge1), IFE2 = (T0,G0,Sedge2)
            // IFF = (T0,G0,Sface)
            // => l3dFE += IFE1, l3dF += iFF,
            //    l3dFEresi += IFE2
            let found3 = {
                let intf = h1.read().expect("interference lock");
                fun_find_3d_is_esame_isf(&intf, l3d_f)
            };
            if found3.is_some() {
                l3d_fe_resi.push(h1);
            } else {
                lfe.push(h1);
            }
            l3d_fe.remove(i);
        }
    } // it1
    l3d_f.len() as i32
}

/// `FUN_select3dinterference(SIX, BDS, lF, l3dF, lFE, lFEresi, l3dFE,
/// l3dFEresi, l2dFE)` (`cxx:1147-1185`).
///
/// `lF = {I=(T(face),G=POINT/VERTEX,S=FACE)}` and `lFE` the same on
/// `S=EDGE`; on return `l3dFE`/`l3dF` hold the 3d interferences to reduce,
/// `l2dFE` the 2d ones, and `lFEresi` the leftovers to remove.
pub fn fun_select_3d_interference(
    six: i32,
    bds: &mut DataStructure,
    lf: &mut InterferenceList,
    l3d_f: &mut InterferenceList,
    lfe: &mut InterferenceList,
    lfe_resi: &mut InterferenceList,
    l3d_fe: &mut InterferenceList,
    l3d_fe_resi: &mut InterferenceList,
    l2d_fe: &mut InterferenceList,
) -> Result<(), String> {
    let e = bds
        .shape(six, true)
        .ok_or_else(|| "FUN_select3dinterference: BDS.Shape(SIX)".to_string())?;
    if BRepTool::is_degenerated(&Edge(e)) {
        return Ok(());
    }

    l3d_f.clear();
    l3d_fe.clear();

    fun_select_3d_i(six, bds, lfe, lfe_resi, l3d_fe)?;
    fun_select_3d_is_esame_isf(lfe, l3d_fe, l3d_fe_resi, lf, l3d_f);
    fun_select_2d_i(six, bds, ShapeType::Face, lfe, l2d_fe)?;
    Ok(())
}
