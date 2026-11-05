//! `RemoveSD`, `CpOnEdge` and `RemoveSurfData` — the surf-data pruning of
//! `ChFi3d_Builder_CnCrn.cxx`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_CnCrn.cxx:222-238`
//! (`RemoveSD`), `:957-983` (`CpOnEdge`) and `:986-1080` (`RemoveSurfData`)
//! (OCCT 8.0.0).
//!
//! `RemoveSurfData` walks `myVDataMap(V1)` and drops the surf data of each
//! stripe that cannot reach one of the two edges adjacent to `V1`: from the
//! start of the sequence when the vertex is the spine's last one, from the end
//! when it is the first. The sequence is shortened in place
//! (`It.ChangeValue()`), so the port writes through the stripe handle.

use crate::chfi3d::builder_0_helpers::ch_fi3d_index_of_surf_data;
use crate::chfi3d::builder_helpers::{cherche_element, edge_common_faces};
use crate::chfi_ds::{ChFiDSMap, Stripe, StripeHandle, StripeMap};
use crate::shape::{Edge, Face, Vertex};
use crate::topo_tools::is_same;

use super::spine_edge_at;

/// `RemoveSD(Stripe, num1, num2)` (`cxx:222-238`): remove the 1-based inclusive
/// range `[num1, num2]` of the stripe's surf-data sequence.
///
/// An empty sequence returns immediately in OCCT; an out-of-range removal
/// raises from `NCollection_Sequence`, which the port answers by clamping.
pub fn remove_sd(stripe: &mut Stripe, num1: usize, num2: usize) {
    let seq = stripe.change_set_of_surf_data();
    if seq.is_empty() || num1 == 0 || num2 < num1 {
        return;
    }
    let start = num1 - 1;
    let end = num2.min(seq.len());
    if start >= end {
        return;
    }
    seq.drain(start..end);
}

/// `CpOnEdge(stripe, num, isfirst, Eadj1, Eadj2, compoint)` (`cxx:957-983`):
/// whether surf data `num` carries a common point on `Eadj1` or `Eadj2`.
pub fn cp_on_edge(
    stripe: &StripeHandle,
    num: usize,
    is_first: bool,
    eadj1: &Edge,
    eadj2: &Edge,
) -> Result<bool, String> {
    let guard = stripe
        .read()
        .map_err(|_| "CpOnEdge: poisoned stripe".to_string())?;
    let data = num
        .checked_sub(1)
        .and_then(|i| guard.set_of_surf_data().get(i))
        .ok_or_else(|| format!("CpOnEdge: SurfData {num} is out of range"))?;

    let mut compoint = false;
    // `ChFiDS_CommonPoint cp1 = ...ChangeVertex(isfirst, 1); cp2 = ...(isfirst, 2);`
    // (`cxx:968-969`).
    for on_s in [1, 2] {
        let cp = data.vertex(is_first, on_s);
        if cp.is_on_arc() {
            if let Some(arc) = cp.arc() {
                if is_same(arc, &eadj1.0) || is_same(arc, &eadj2.0) {
                    compoint = true;
                }
            }
        }
    }
    Ok(compoint)
}

/// `RemoveSurfData(myVDataMap, myEFMap, edgecouture, facecouture, V1)`
/// (`cxx:986-1080`).
pub fn remove_surf_data(
    v_data_map: &StripeMap,
    ef_map: &ChFiDSMap,
    edge_couture: Option<&Edge>,
    face_couture: Option<&Face>,
    v1: &Vertex,
) -> Result<(), String> {
    let stripes: Vec<StripeHandle> = v_data_map
        .find_from_key(&v1.0)
        .map(|l| l.to_vec())
        .unwrap_or_default();

    for stripe in &stripes {
        // `nbsurf = It.Value()->SetOfSurfData()->Length();` (`cxx:1005`).
        let nbsurf = {
            let guard = stripe
                .read()
                .map_err(|_| "RemoveSurfData: poisoned stripe".to_string())?;
            guard.set_of_surf_data().len()
        };
        if nbsurf == 1 {
            continue;
        }

        // `num = ChFi3d_IndexOfSurfData(V1, It.Value(), sense);` (`cxx:1009`).
        let mut sense = 0;
        let num = ch_fi3d_index_of_surf_data(&v1.0, stripe, &mut sense)?;

        // `if (sense == 1) Ecur = ...Edges(1); else Ecur = ...Edges(nbedge);`
        // (`cxx:1010-1016`).
        let ecur = spine_edge_at(stripe, sense == 1)?;

        // `ChFi3d_edge_common_faces(myEFMap(Ecur), F1, F2);` (`cxx:1018`).
        let (f1, f2) = edge_common_faces(ef_map.get_from_key(&ecur.0))
            .ok_or_else(|| "RemoveSurfData: the spine edge has no face".to_string())?;

        // `if (F1.IsSame(facecouture)) Eadj1 = edgecouture; else
        //  ChFi3d_cherche_element(V1, Ecur, F1, Eadj1, Vbid);` (`cxx:1019-1027`).
        let eadj1 = match face_couture {
            Some(fc) if is_same(&f1.0, &fc.0) => edge_couture.cloned(),
            _ => Some(cherche_element(v1, &ecur, &f1)?.0),
        };
        let eadj1 = eadj1.ok_or_else(|| "RemoveSurfData: Eadj1 is null".to_string())?;
        // `ChFi3d_edge_common_faces(myEFMap(Eadj1), Fg, Fd);` — the answer is
        // dropped by OCCT.
        let _ = edge_common_faces(ef_map.get_from_key(&eadj1.0));

        // The `F2` mirror (`cxx:1028-1036`).
        let eadj2 = match face_couture {
            Some(fc) if is_same(&f2.0, &fc.0) => edge_couture.cloned(),
            _ => Some(cherche_element(v1, &ecur, &f2)?.0),
        };
        let eadj2 = eadj2.ok_or_else(|| "RemoveSurfData: Eadj2 is null".to_string())?;
        let _ = edge_common_faces(ef_map.get_from_key(&eadj2.0));

        let is_first = sense == 1;
        if is_first {
            // Among the surf data find the greatest index `ind` that can reach
            // `Eadj1` / `Eadj2`, and remove 1 .. ind-1 (`cxx:1041-1058`).
            let mut ind = 0usize;
            for i in 1..=nbsurf {
                if cp_on_edge(stripe, i, is_first, &eadj1, &eadj2)? {
                    ind = i;
                }
            }
            if ind >= 2 {
                let mut guard = stripe
                    .write()
                    .map_err(|_| "RemoveSurfData: poisoned stripe".to_string())?;
                remove_sd(&mut guard, 1, ind - 1);
            }
        } else {
            // The `num` mirror: find the smallest index that can reach the two
            // edges, and remove ind+1 .. num (`cxx:1059-1075`).
            let mut ind = num;
            let mut i = num;
            while i >= 1 {
                if cp_on_edge(stripe, i as usize, is_first, &eadj1, &eadj2)? {
                    ind = i;
                }
                i -= 1;
            }
            if ind < num {
                let mut guard = stripe
                    .write()
                    .map_err(|_| "RemoveSurfData: poisoned stripe".to_string())?;
                remove_sd(&mut guard, (ind + 1) as usize, num as usize);
            }
        }
    }
    Ok(())
}
