//! ShapeExtend_WireData operations used by ComposeShell.
//! The OCCT class is represented by an ordered Vec of Edge throughout this
//! module, so these are plain vector/edge operations.

use crate::shape::Edge;
use crate::shape::Face;
use crate::tgeometry::GeometryRegistry;

/// ShapeExtend_WireData::ComputeSeams (ShapeExtend_WireData.cxx:163-222):
/// positions (1-based) of the seam occurrences of `edges`. The first forward
/// occurrence that IsSame as a reversed one becomes (seamF, seamR); further
/// pairs are appended to the returned list as (forward, reversed).
fn compute_seams(edges: &[Edge]) -> (i32, i32, Vec<i32>) {
    use std::collections::HashMap;
    let nb = edges.len();
    // First pass: map every REVERSED edge by TShape identity to its rank
    // (OCCT uses an IndexedMap; a later duplicate overwrites SE[num]).
    let mut reversed_rank: HashMap<usize, usize> = HashMap::new();
    for i in 1..=nb {
        let e = &edges[i - 1];
        if e.0.orientation().is_reversed() {
            reversed_rank.insert(GeometryRegistry::shape_key(&e.0), i);
        }
    }
    let mut seam_f = 0i32;
    let mut seam_r = 0i32;
    let mut seams: Vec<i32> = Vec::new();
    // Second pass: FORWARD edges that match a REVERSED one are seams.
    for i in 1..=nb {
        let e = &edges[i - 1];
        if e.0.orientation().is_reversed() {
            continue;
        }
        if let Some(&r) = reversed_rank.get(&GeometryRegistry::shape_key(&e.0)) {
            if seam_f == 0 {
                seam_f = i as i32;
                seam_r = r as i32;
            } else {
                seams.push(i as i32);
                seams.push(r as i32);
            }
        }
    }
    (seam_f, seam_r, seams)
}

/// ShapeExtend_WireData::SwapSeam (ShapeExtend_WireData.cxx:511-543): on the
/// FORWARD occurrence only, exchange PCurve/PCurve2 and restore the range of
/// the forward pcurve.
fn swap_seam(edge: &Edge, face: &Face) {
    if edge.0.orientation().is_reversed() {
        return; // cxx:518-521: le faire une seule fois
    }
    let face_key = GeometryRegistry::shape_key(&face.0);
    let reg = GeometryRegistry::global();
    let mut pcs = reg.edge_pcurves(&edge.0, face_key);
    if pcs.len() < 2 {
        return; // cxx:534-537: c2df or c2dr null
    }
    // cxx:531: uff/ulf from the FORWARD occurrence, before the swap.
    let (uff, ulf) = reg
        .pcurve_range(&edge.0, face_key)
        .unwrap_or((0.0, 0.0));
    pcs.swap(0, 1); // cxx:541: B.UpdateEdge(E, c2dr, c2df, face, 0.)
    reg.set_edge_pcurves(&edge.0, face_key, pcs);
    reg.set_pcurve_range(&edge.0, face_key, uff, ulf); // cxx:542
}

/// ShapeExtend_WireData::Reverse(face) (ShapeExtend_WireData.cxx:545-572):
/// Reverse() then ComputeSeams(true) and SwapSeam on mySeamF / mySeamR /
/// mySeams.
pub fn reverse_wire_data_on_face(edges: &mut [Edge], face: &Face) {
    reverse_wire_data(edges); // cxx:547
    let (seam_f, seam_r, seams) = compute_seams(edges); // cxx:557
    if seam_f > 0 {
        swap_seam(&edges[(seam_f - 1) as usize], face); // cxx:558-561
    }
    if seam_r > 0 {
        swap_seam(&edges[(seam_r - 1) as usize], face); // cxx:562-565
    }
    for s in seams {
        swap_seam(&edges[(s - 1) as usize], face); // cxx:566-570
    }
}

/// ShapeExtend_WireData::Reverse() (ShapeExtend_WireData.cxx:483-506): reverse
/// every edge and permute them, so the wire as a whole is reversed.
///
/// ShapeExtend_WireData::Reverse(face) additionally runs ComputeSeams(true)
/// and then SwapSeam on mySeamF / mySeamR / mySeams (cxx:508-572).
///
/// SwapSeam itself is already ported in this port as
/// crate::meshing::model_builder::swap_seam_pcurves (wire_builder.rs:448-459,
/// now pub(crate)), which does the FWD/REV pcurve exchange with the same
/// len() < 2 guard; a Reverse(face) wrapper still needs ComputeSeams,
/// which is UNPORTED (cxx:557).
pub fn reverse_wire_data(edges: &mut [Edge]) {
    let nb = edges.len();
    if nb == 0 {
        return;
    }
    for i in 1..=nb / 2 {
        // 1-based i: S1 = edges(i), S2 = edges(nb + 1 - i).
        let mut s1 = edges[i - 1].clone();
        s1.0.reverse();
        let mut s2 = edges[nb - i].clone();
        s2.0.reverse();
        edges[i - 1] = s2;
        edges[nb - i] = s1;
    }
    if nb % 2 == 1 {
        let i = (nb + 1) / 2;
        let mut si = edges[i - 1].clone();
        si.0.reverse();
        edges[i - 1] = si;
    }
}
