//! `FUN_reclSE` / `FUN_reclSE2` — the closed-support-edge 1d reduction of
//! `TopOpeBRepDS_EIR`.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_EIR.cxx` (OCCT 8.0.0):
//! - `FUN_reclSE2` (`:507-607`)
//! - `FUN_reclSE` (`:610-648`)
//!
//! On a reference edge same domain with a *closed* support edge, the two
//! interferences `(OU/IN(SE), VG, SE)` and `(IN/OU(SE), VG, SE)` at the closing
//! vertex reduce to `(IN/IN(SE), VG, SE)`: `I1`'s transition is set to
//! `INTERNAL` and `I1` moves to the reduced list (`cxx:595-601`).

use crate::abs::{Orientation, ShapeType, State};

use super::super::data_structure::DataStructure;
use super::super::export::{fds_data, fun_ds_sdm};
use super::super::interference_iterator::InterferenceList;
use super::super::kind::Kind;
use super::super::tki::TKI;
use crate::topopebrep::tool::closed::closed_e;

/// `FUN_reclSE2(SIX, BDS, LI, reducedLI)` (`cxx:507-607`).
///
/// `LI = {I = (T(edge), G0, edge)}`; `SIX` is same domain with the support edge
/// `SE` of `I`. `I1` is only reduced when its geometry is the closing vertex of
/// `SE` (`cxx:541-556`).
pub fn fun_recl_se2(
    six: i32,
    bds: &DataStructure,
    li: &mut InterferenceList,
    reduced_li: &mut InterferenceList,
) -> Result<(), String> {
    reduced_li.clear();

    let e = bds
        .shape(six, true)
        .ok_or_else(|| "FUN_reclSE2: BDS.Shape(SIX)".to_string())?;

    let mut i1 = 0;
    while i1 < li.len() {
        let h1 = li[i1].clone();
        let (_gt1, g1, _st1, s1) = {
            let i = h1.read().expect("interference lock");
            fds_data(&i)
        };
        let o1 = h1
            .read()
            .expect("interference lock")
            .transition()
            .orientation(State::In, ShapeType::Face);
        if o1 == Orientation::Internal || o1 == Orientation::External {
            i1 += 1;
            continue;
        }
        let co1 = o1.complement();

        let v1 = bds
            .shape(g1, true)
            .ok_or_else(|| "FUN_reclSE2: BDS.Shape(G1)".to_string())?;
        let e1 = bds
            .shape(s1, true)
            .ok_or_else(|| "FUN_reclSE2: BDS.Shape(S1)".to_string())?;
        let Some(vclo1) = closed_e(&e1) else {
            i1 += 1;
            continue;
        };
        if !crate::topo_tools::is_same(&vclo1, &v1) {
            i1 += 1;
            continue;
        }

        if !fun_ds_sdm(bds, &e, &e1) {
            i1 += 1;
            continue;
        }

        if i1 + 1 >= li.len() {
            i1 += 1;
            continue;
        }

        let mut has_co = false;
        let mut i2 = i1 + 1;
        while i2 < li.len() {
            let h2 = li[i2].clone();
            let (o2, _st2, s2) = {
                let i = h2.read().expect("interference lock");
                (
                    i.transition().orientation(State::In, ShapeType::Face),
                    i.support_type(),
                    i.support(),
                )
            };
            if s1 != s2 {
                i2 += 1;
                continue;
            }
            if o2 != co1 {
                i2 += 1;
                continue;
            }

            li.remove(i2);
            has_co = true;
            break;
        } // it2

        if has_co {
            h1.write()
                .expect("interference lock")
                .change_transition()
                .set_orientation(Orientation::Internal);
            reduced_li.push(h1);
            li.remove(i1);
        } else {
            i1 += 1;
        }
    } // it1
    Ok(())
}

/// `FUN_reclSE(EIX, BDS, LOI, RLOI)` (`cxx:610-648`).
///
/// Runs [`fun_recl_se2`] on every `(VERTEX, G)` group of two or more. Groups of
/// any other geometry kind are *dropped*: `cxx:626-629` rebuilds `LOI` from the
/// vertex groups only.
pub fn fun_recl_se(
    eix: i32,
    bds: &DataStructure,
    loi: &mut InterferenceList,
    rloi: &mut InterferenceList,
) -> Result<(), String> {
    let mut tki = TKI::new();
    tki.fill_on_geometry(loi);
    loi.clear();

    let mut order: Vec<(Kind, i32)> = Vec::new();
    tki.init();
    while tki.more() {
        if let Some((k, g, _)) = tki.value() {
            order.push((k, g));
        }
        tki.next();
    }

    for (k, g) in order {
        if k != Kind::Vertex {
            continue;
        }

        let mut loi_g: InterferenceList = tki.interferences(k, g).to_vec();
        let nloi = loi_g.len();
        if nloi == 0 {
            continue;
        } else if nloi == 1 {
            loi.extend(loi_g);
        } else {
            let mut rloi_g: InterferenceList = Vec::new();
            fun_recl_se2(eix, bds, &mut loi_g, &mut rloi_g)?;
            loi.extend(loi_g);
            rloi.extend(rloi_g);
        }
    }
    Ok(())
}
