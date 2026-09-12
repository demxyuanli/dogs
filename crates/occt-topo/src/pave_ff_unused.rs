//! Unused section vertices collected by `PostTreatFF` (`_6.cxx:1203-1240`).
//!
//! For each F/F pair, `GetStickVertices` plus `RemoveUsedVertices` leaves the
//! stick vertices that were never put on a section curve. Those unused
//! vertices are fed into the nested PostTreatFF filler together with the
//! section edges, so they can fuse with section vertices.

use std::collections::HashSet;

use crate::pave_ff_misc::{get_stick_vertices, remove_used_vertices};
use crate::pave_filler::PaveFiller;
use crate::shape::TopoShape;

/// `PostTreatFF` unused-vertex harvest (`_6.cxx:1203`).
pub fn unused_stick_vertices(f: &PaveFiller) -> Vec<TopoShape> {
    let nb_ff = f.ds().interf_ff().len();
    let mut out: Vec<TopoShape> = Vec::new();
    let mut seen: HashSet<usize> = HashSet::new();
    for i in 0..nb_ff {
        let (n_f1, n_f2) = f.ds().interf_ff()[i].indices();
        let (mut mv, _mv_ef, _mi) = get_stick_vertices(f.ds(), n_f1, n_f2);
        remove_used_vertices(f.ds().interf_ff()[i].curves(), &mut mv);
        for n_v in mv {
            if !seen.insert(n_v) {
                continue;
            }
            if let Some(s) = f.ds().shape(n_v).cloned() {
                out.push(s);
            }
        }
    }
    out
}

/// True when a vertex index still has an ext-pave on any curve of `i`.
pub fn vertex_used_on_ff_curves(f: &PaveFiller, i: usize, n_v: usize) -> bool {
    for nc in f.ds().interf_ff()[i].curves() {
        for pb in nc.pave_blocks() {
            if pb.index1 == n_v || pb.index2 == n_v {
                return true;
            }
            if pb.ext_paves().iter().any(|p| p.index == n_v) {
                return true;
            }
        }
    }
    false
}
