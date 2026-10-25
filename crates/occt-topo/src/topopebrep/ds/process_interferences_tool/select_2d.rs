//! `TopOpeBRepDS_ProcessInterferencesTool` — the 1d / 2d selection filters.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_ProcessInterferencesTool.cxx`
//! (OCCT 8.0.0):
//! - `FUN_select2dI` (`:847-977`)
//! - `FUN_selectpure2dI` (`:981-1025`)
//! - `FUN_select1dI` (`:1028-1145`)
//!
//! As in [`super`], a `FUN_select*` walks the input list and moves the selected
//! element to the output list; the return value is the *output* list extent
//! (`cxx:975`, `:1023`, `:1143`), not the number selected by this call.
//! `Transition::Index()` raises in OCCT when the two transition sides disagree,
//! which the port reports as `Err`.

use crate::abs::ShapeType;

use super::super::data_structure::DataStructure;
use super::super::export::{fds_copy, fds_data, fds_tdata, fun_ds_sdm};
use super::super::interference_iterator::InterferenceList;
use super::super::kind::Kind;
use super::super::tki::TKI;
use super::fun_select_sk_interference;

/// `FUN_select2dI(SIX, BDS, TRASHAk, lI, l2dI)` (`cxx:847-977`).
///
/// `TRASHAk == FACE` selects the 2d interferences on same-domain faces, `EDGE`
/// the ones on the same edge. The first walk (`cxx:866-891`, "xpu201098") drops
/// the interferences whose transition index is the support of a 3d interference
/// on a face of `SIX`: they are 3d face interferences, not 2d ones.
///
/// The second walk pairs two interferences on the same
/// `(SB, IB, IA)` transition and the same `(GT, G, ST)` geometry; the pair is a
/// 2d interference when they end on the same index `IB1 == IB2`, or (on faces)
/// when the two faces are same-domain.
pub fn fun_select_2d_i(
    six: i32,
    bds: &DataStructure,
    trasha_k: ShapeType,
    li: &mut InterferenceList,
    l2d_i: &mut InterferenceList,
) -> Result<i32, String> {
    l2d_i.clear();
    if li.len() <= 1 {
        return Ok(0);
    }

    let ton_face = trasha_k == ShapeType::Face;

    // xpu201098 : check FEI is not 3dFEI (cto904F6,SIX10,G6)
    let mut map_ftra: Vec<i32> = Vec::new();
    let mut lie: InterferenceList = Vec::new();
    fds_copy(
        &bds.shape_interferences_at(six, true).to_vec(),
        &mut lie,
    );
    let mut l3d_f: InterferenceList = Vec::new();
    fun_select_sk_interference(&mut lie, Kind::Face, &mut l3d_f);
    for h in l3d_f.iter() {
        map_ftra.push(h.read().expect("interference lock").support());
    }

    let mut lii: InterferenceList = Vec::new();
    let mut i = 0;
    while i < li.len() {
        let itrasha = li[i]
            .read()
            .expect("interference lock")
            .transition()
            .index()?;
        if map_ftra.contains(&itrasha) {
            i += 1;
        } else {
            lii.push(li[i].clone());
            li.remove(i);
        }
    }

    let mut i1 = 0;
    while i1 < lii.len() {
        let h1 = lii[i1].clone();
        let (gt1, g1, st1, _s1, sb1, ib1, _sa1, ia1) = {
            let intf = h1.read().expect("interference lock");
            let (gt, g, st, s) = fds_data(&intf);
            let (sb, ib, sa, ia) = fds_tdata(&intf.transition());
            (gt, g, st, s, sb, ib, sa, ia)
        };

        if sb1 != trasha_k {
            i1 += 1;
            continue;
        }

        if i1 + 1 >= lii.len() {
            break;
        }

        let mut complex2d = false;
        let mut i2 = i1 + 1;
        while i2 < lii.len() {
            let h2 = lii[i2].clone();
            let (gt2, g2, st2, _s2, sb2, ib2, _sa2, ia2) = {
                let intf = h2.read().expect("interference lock");
                let (gt, g, st, s) = fds_data(&intf);
                let (sb, ib, sa, ia) = fds_tdata(&intf.transition());
                (gt, g, st, s, sb, ib, sa, ia)
            };

            let cond = sb1 == sb2 && ib1 == ib2 && ia1 == ia2;
            if !cond {
                i2 += 1;
                continue;
            }

            let same_g = gt2 == gt1 && g2 == g1;
            let same_st = st2 == st1;
            if !(same_g && same_st) {
                i2 += 1;
                continue;
            }

            complex2d = ib1 == ib2;
            if ton_face {
                // `BDS.Shape(IB1)` / `BDS.Shape(IB2)` are faces of a valid DS.
                match (bds.shape(ib1, true), bds.shape(ib2, true)) {
                    (Some(f1), Some(f2)) => complex2d = complex2d || fun_ds_sdm(bds, &f1, &f2),
                    _ => {}
                }
            }
            if complex2d {
                l2d_i.push(h2);
                lii.remove(i2);
            } else {
                i2 += 1;
            }
        } // it2
        if complex2d {
            l2d_i.push(h1);
            lii.remove(i1);
        } else {
            i1 += 1;
        }
    } // it1
    li.extend(lii);
    Ok(l2d_i.len() as i32)
}

/// `FUN_selectpure2dI(lF, lFE, l2dFE)` (`cxx:981-1025`).
///
/// Keeps the `(T(face), G, S=edge)` interferences whose transition ends on the
/// same index (`IB == IA`) and for which `lF` holds no interference on that
/// index: no `(T(face), G, S=face)` counterpart exists, so the interference
/// cannot be part of a 3d reduction.
pub fn fun_select_pure2d_i(
    lf: &InterferenceList,
    lfe: &mut InterferenceList,
    l2d_fe: &mut InterferenceList,
) -> i32 {
    l2d_fe.clear();
    let mut i = 0;
    while i < lfe.len() {
        let (ib, ia) = {
            let intf = lfe[i].read().expect("interference lock");
            (
                intf.transition().index_before(),
                intf.transition().index_after(),
            )
        };
        if ib != ia {
            i += 1;
            continue;
        }

        let mut found_if = false;
        for h in lf.iter() {
            if h.read().expect("interference lock").support() == ib {
                found_if = true;
                break;
            }
        }
        if found_if {
            i += 1;
            continue;
        }

        let moved = lfe.remove(i);
        l2d_fe.push(moved);
    }
    l2d_fe.len() as i32
}

/// `FUN_select1dI(SIX, BDS, LI, l1dI)` (`cxx:1028-1145`).
///
/// `LI = {I = (T(ES), G, ES)}` on the edge `SIX`; the interferences grouped by
/// `(Kind, G)` through `TopOpeBRepDS_TKI` (`cxx:1047-1049`) are reduced on a
/// vertex geometry: a pair `(SB, IB) == (EDGE, IA)` with `EIX` same-domain with
/// the support edge `S1`/`S2` is a 1d interference.
///
/// As in OCCT, `LI` is rebuilt from the interferences that the walk kept, in
/// `TKI` group order, and the elements left unvisited by a `break`
/// (`cxx:1097-1101`) are dropped.
pub fn fun_select_1d_i(
    six: i32,
    bds: &DataStructure,
    li: &mut InterferenceList,
    l1d_i: &mut InterferenceList,
) -> Result<i32, String> {
    l1d_i.clear();
    if li.len() <= 1 {
        return Ok(0);
    }

    let mut new_li: InterferenceList = Vec::new();
    let eix = bds.shape(six, true);

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
        let loi: InterferenceList = tki.interferences(k, g).to_vec();
        let nloi = loi.len();

        let mut ok = k == Kind::Vertex;
        ok = ok && nloi > 1;
        if !ok {
            new_li.extend(loi.iter().cloned());
            continue;
        }

        let mut i1 = 0;
        while i1 < loi.len() {
            let h1 = loi[i1].clone();
            let (_gt1, _g1, _st1, s1, sb1, ib1, _sa1, ia1) = {
                let intf = h1.read().expect("interference lock");
                let (gt, g, st, s) = fds_data(&intf);
                let (sb, ib, sa, ia) = fds_tdata(&intf.transition());
                (gt, g, st, s, sb, ib, sa, ia)
            };

            let mut cond1 = sb1 == ShapeType::Edge && ib1 == ia1;
            if cond1 {
                cond1 = match (eix.as_ref(), bds.shape(s1, true)) {
                    (Some(a), Some(b)) => fun_ds_sdm(bds, a, &b),
                    _ => false,
                };
            }
            if !cond1 {
                new_li.push(h1);
                i1 += 1;
                continue;
            }

            if i1 + 1 >= loi.len() {
                break;
            }

            let mut complex1d = false;
            let mut i2 = i1 + 1;
            while i2 < loi.len() {
                let h2 = loi[i2].clone();
                let (_gt2, _g2, _st2, s2, sb2, ib2, _sa2, ia2) = {
                    let intf = h2.read().expect("interference lock");
                    let (gt, g, st, s) = fds_data(&intf);
                    let (sb, ib, sa, ia) = fds_tdata(&intf.transition());
                    (gt, g, st, s, sb, ib, sa, ia)
                };

                let cond2 = sb1 == sb2 && ib2 == ia2;
                if !cond2 {
                    new_li.push(h2);
                    i2 += 1;
                    continue;
                }
                complex1d = match (eix.as_ref(), bds.shape(s2, true)) {
                    (Some(a), Some(b)) => fun_ds_sdm(bds, a, &b),
                    _ => false,
                };
                if !complex1d {
                    new_li.push(h2);
                    i2 += 1;
                    continue;
                }
                l1d_i.push(h2);
                i2 += 1;
            } // it2
            if complex1d {
                l1d_i.push(h1);
            }
            i1 += 1;
        } // it1
    } // tki

    li.clear();
    li.extend(new_li);
    Ok(l1d_i.len() as i32)
}
