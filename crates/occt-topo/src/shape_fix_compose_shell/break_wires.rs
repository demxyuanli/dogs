//! Port of ShapeFix_ComposeShell::BreakWires (ShapeFix_ComposeShell.cxx:2279-2385).
//!
//! myLoc is the identity here (port shapes carry no location), so the only
//! location-dependent code in this function (none) collapses.

use std::collections::HashSet;
use std::sync::Arc;

use crate::abs::Orientation;
use crate::shape::Vertex;

use super::shell::ComposeShell;
use super::wire_segment::WireSegment;

/// TopTools_ShapeMapHasher key: TShape identity (IsSame).
fn vkey(v: &Vertex) -> usize {
    Arc::as_ptr(&v.0.tshape) as usize
}

impl ComposeShell {
    /// ShapeFix_ComposeShell::BreakWires (cxx:2279-2385): split every wire
    /// segment at the vertices collected from EXTERNAL wires.
    pub fn break_wires(&mut self, seqw: &mut Vec<WireSegment>) {
        // cxx:2281-2306: first collect splitting vertices from the EXTERNAL /
        // INTERNAL wire segments.
        let mut split_vertices: HashSet<usize> = HashSet::new();
        for seg in seqw.iter() {
            let ori_wire = seg.orientation();
            if ori_wire != Orientation::External && ori_wire != Orientation::Internal {
                continue;
            }
            for e in seg.edges() {
                let ori_edge = if ori_wire == Orientation::External {
                    ori_wire
                } else {
                    e.0.orientation()
                };
                if ori_edge == Orientation::External {
                    if let Some(v) = crate::shhealing::first_vertex(e) {
                        split_vertices.insert(vkey(&v));
                    }
                    if let Some(v) = crate::shhealing::last_vertex(e) {
                        split_vertices.insert(vkey(&v));
                    }
                }
            }
        }

        // cxx:2310-2384: then split each wire. Each wire is supposed to be
        // connected (probably not closed).
        let mut i = 0usize; // 0-based position in seqw (OCCT uses 1-based i).
        while i < seqw.len() {
            let ori = seqw[i].orientation();
            if seqw[i].is_vertex() {
                i += 1;
                continue;
            }
            let wire = seqw[i].clone();
            let sbwd: Vec<crate::shape::Edge> = wire.edges().to_vec();
            let nb = sbwd.len();

            // cxx:2321-2329: find the first vertex for split.
            let mut j = 0usize;
            while j < nb {
                let found = crate::shhealing::first_vertex(&sbwd[j])
                    .map(|v| split_vertices.contains(&vkey(&v)))
                    .unwrap_or(false);
                if found {
                    break;
                }
                j += 1;
            }
            if j >= nb {
                i += 1;
                continue; // cxx:2330-2333: splitting not needed.
            }

            // cxx:2335-2345: if the first split of a closed edge is not its
            // start, make a permutation.
            let mut shift = 0usize;
            if j > 0 && !self.closed_mode && wire.is_closed() {
                let found = crate::shhealing::first_vertex(&sbwd[0])
                    .map(|v| split_vertices.contains(&vkey(&v)))
                    .unwrap_or(false);
                if !found {
                    shift = j;
                }
            }

            // cxx:2347-2383: perform the splitting.
            let mut nbnew = 0usize;
            let mut newwire = WireSegment::new();
            let mut cur_ori = ori;
            for ind in 0..nb {
                let j = (ind + shift) % nb;
                let mut edge = sbwd[j].clone();
                let at_first = crate::shhealing::first_vertex(&edge)
                    .map(|v| split_vertices.contains(&vkey(&v)))
                    .unwrap_or(false);
                if ind == 0 || at_first {
                    if newwire.nb_edges() != 0 {
                        newwire.set_orientation(cur_ori);
                        // cxx:2362: seqw.InsertBefore(i++, newwire).
                        seqw.insert(i, newwire.clone());
                        i += 1;
                        nbnew += 1;
                    }
                    newwire.clear();
                    cur_ori = ori;
                }
                let (iumin, iumax, ivmin, ivmax) =
                    wire.get_patch_index(j + 1).unwrap_or((0, 0, 0, 0));
                if ori == Orientation::Internal && edge.0.orientation() == Orientation::External {
                    cur_ori = Orientation::External;
                    edge.0.set_orientation(Orientation::Forward);
                    nbnew += 1;
                }
                newwire.add_edge_patch(0, edge, iumin, iumax, ivmin, ivmax);
            }
            if nbnew != 0 {
                newwire.set_orientation(cur_ori);
                // cxx:2382: seqw.SetValue(i, newwire).
                seqw[i] = newwire;
            }
            i += 1;
        }
    }
}
