//! `BRepFill_CompatibleWires::ComputeOrigin` (`cxx:1787-2126`).
//!
//! The closed-wire branch used when every section already carries the same
//! number of edges: it does not cut anything, it only re-origins and re-orients
//! every section so that the ruled surface between two sections does not twist.
//!
//! Each section is compared with the previous one over all `theLength`
//! rotations, forwards and backwards, by the sum of squared distances between
//! corresponding vertices (and, for sections of at most two edges, between
//! intermediate curve samples); the best rotation wins and the wire is
//! re-stacked in that order.
//!
//! Not ported, and why: `cxx:2128` opens an `#ifdef OCCT_DEBUG_EFV` block that
//! ends just before the method's closing brace. That block re-derives the
//! per-section origin with `SeqOfVertices` / `Transform` / an
//! `EdgesFromVertex`-based guess and is the only reader of the `polar`
//! argument (`const bool /*polar*/` in the signature, `cxx:1787`), so a default
//! build never compiles it. `polar` therefore stays unread here.
//!
//! Deviations, both OCCT assertion sites:
//! * `PrevSeq(n)` / `SeqVertices(k)` past the end assert in
//!   `NCollection_Sequence::Value`; a missing entry contributes no distance.
//! * a missing `BRepAdaptor_Curve` (`BRep_Tool::Curve` on a pcurve-only edge)
//!   would throw on `Value`; the sample contributes no distance.

use occt_core::gp::{GpPnt, GpVec};
use occt_core::precision::INFINITE;
use occt_geom::Curve;

use crate::abs::Orientation;
use crate::brep_fill_compatible_wires::statics::plane_of_wire_full;
use crate::brep_fill_compatible_wires::CompatibleWires;
use crate::brep_fill_generator::ThruSectionErrorStatus;
use crate::brep_tool::BRepTool;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Vertex, Wire};
use crate::topo_tools::{is_same, wire_vertices};

/// `ComputeOrigin(polar)` (`cxx:1787-2126`); `polar` is only read by the
/// `OCCT_DEBUG_EFV` block (`cxx:2215`) and is therefore unused.
pub(crate) fn compute_origin(cw: &mut CompatibleWires, _polar: bool) -> Result<(), String> {
    // `cxx:1799-1808`.
    let nb_sects = cw.my_work.len();
    let mut ideb = 1usize;
    let mut ifin = nb_sects;
    if cw.my_degen1 {
        ideb += 1;
    }
    if cw.my_degen2 {
        ifin -= 1;
    }
    // `cxx:1809`: only the `OCCT_DEBUG_EFV` block (`cxx:2135`) reads it.
    let _v_closed = (!cw.my_degen1)
        && (!cw.my_degen2)
        && ideb <= ifin
        && is_same(&cw.my_work[ideb - 1], &cw.my_work[ifin - 1]);

    // `cxx:1811-1838`: every section must be closed.
    let mut all_closed = true;
    if ideb <= ifin {
        for i in ideb..=ifin {
            let w = Wire(cw.my_work[i - 1].clone());
            let mut w_closed = w.0.is_closed();
            if !w_closed {
                let (v1, v2) = wire_vertices(&w);
                if let (Some(v1), Some(v2)) = (v1, v2) {
                    if is_same(&v1.0, &v2.0) {
                        w_closed = true;
                    }
                }
            }
            all_closed = all_closed && w_closed;
        }
    }
    if !all_closed {
        cw.my_status = ThruSectionErrorStatus::NotSameTopology;
        return Ok(());
    }

    // `cxx:1872-1881`: the first section's vertex / edge cycles.
    let wire = Wire(cw.my_work[ideb - 1].clone());
    let (mut prev_seq, mut prev_eseq) = wire_sequence(&wire);
    let the_length = prev_eseq.len();

    // `cxx:1883-1887`: extra samples only for the two-edge sections.
    let nb_samples = if the_length <= 2 { 4usize } else { 0usize };
    let (first_plane, _) = plane_of_wire_full(&wire);
    let mut prev_bary = first_plane.location();
    let ax1 = first_plane.axis();
    let d = ax1.direction();
    let normal = GpVec::new(d.x(), d.y(), d.z());

    // `cxx:1891-2124`.
    for i in (ideb + 1)..=ifin {
        let a_wire = Wire(cw.my_work[i - 1].clone());
        let (cur_plane, _) = plane_of_wire_full(&a_wire);
        let mut cur_bary = cur_plane.location();
        // `gp_Vec anOffsetProj = (aVec * NormalOfFirstPlane) * NormalOfFirstPlane`
        // (`cxx:1901-1902`): the projection of the barycenter step onto the
        // first plane's normal.
        let a_vec = GpVec::from_pnts(&prev_bary, &cur_bary);
        let an_offset_proj = normal.multiplied_scalar(a_vec.dot(&normal));
        cur_bary.translate_vec(&an_offset_proj.reversed());
        let offset = GpVec::from_pnts(&cur_bary, &prev_bary);

        let (seq_vertices, seq_edges) = wire_sequence(&a_wire);

        let mut min_sum_dist = INFINITE;
        let mut jmin = 1usize;
        let mut forward = false;
        if i == cw.my_work.len() && cw.my_degen2 {
            // `cxx:1921-1925`: the last section of a punctual pair is not
            // rotated.
            jmin = 1;
            forward = true;
        } else {
            for j in 1..=the_length {
                // `cxx:1930-1999`: forward.
                let mut sum = 0.0;
                for (k, n) in (j..=the_length).enumerate() {
                    add_pair(
                        &mut sum,
                        &prev_seq,
                        &prev_eseq,
                        &seq_vertices,
                        seq_edges.get(k),
                        k + 1,
                        n + 1,
                        &offset,
                        nb_samples,
                        false,
                    );
                }
                let mut n = the_length - j + 2;
                for k in 1..j {
                    add_pair(
                        &mut sum,
                        &prev_seq,
                        &prev_eseq,
                        &seq_vertices,
                        seq_edges.get(k - 1),
                        k,
                        n,
                        &offset,
                        nb_samples,
                        false,
                    );
                    n += 1;
                }
                if sum < min_sum_dist {
                    min_sum_dist = sum;
                    jmin = j;
                    forward = true;
                }

                // `cxx:2002-2073`: backward.
                let mut sum = 0.0;
                for (n, k) in (1..=j).zip((1..=j).rev()) {
                    // `cxx:2011-2016`: the paired edge is `SeqEdges(k-1)`,
                    // wrapped to the last edge when `k == 1`.
                    let cur_edge = if k >= 2 {
                        seq_edges.get(k - 2)
                    } else {
                        seq_edges.last()
                    };
                    add_pair(
                        &mut sum,
                        &prev_seq,
                        &prev_eseq,
                        &seq_vertices,
                        cur_edge,
                        k,
                        n,
                        &offset,
                        nb_samples,
                        true,
                    );
                }
                let mut n = j + 1;
                for k in ((j + 1)..=the_length).rev() {
                    add_pair(
                        &mut sum,
                        &prev_seq,
                        &prev_eseq,
                        &seq_vertices,
                        seq_edges.get(k - 2),
                        k,
                        n,
                        &offset,
                        nb_samples,
                        true,
                    );
                    n += 1;
                }
                if sum < min_sum_dist {
                    min_sum_dist = sum;
                    jmin = j;
                    forward = false;
                }
            }
        }

        // `cxx:2079-2120`: re-stack the section in the winning order.
        prev_seq.clear();
        prev_eseq.clear();
        let builder = TopoBuilder::new();
        let mut new_wire = Wire::new();
        if forward {
            for j in jmin..=the_length {
                builder.add(&mut new_wire.0, &seq_edges[j - 1].0);
                prev_seq.push(seq_vertices[j - 1].clone());
                prev_eseq.push(seq_edges[j - 1].clone());
            }
            for j in 1..jmin {
                builder.add(&mut new_wire.0, &seq_edges[j - 1].0);
                prev_seq.push(seq_vertices[j - 1].clone());
                prev_eseq.push(seq_edges[j - 1].clone());
            }
        } else {
            for j in (1..jmin).rev() {
                let mut e = seq_edges[j - 1].clone();
                e.0.reverse();
                builder.add(&mut new_wire.0, &e.0);
                prev_eseq.push(e);
            }
            for j in (jmin..=the_length).rev() {
                let mut e = seq_edges[j - 1].clone();
                e.0.reverse();
                builder.add(&mut new_wire.0, &e.0);
                prev_eseq.push(e);
            }
            // `cxx:2111-2118`: the vertices keep the un-reversed order.
            for j in (1..=jmin).rev() {
                prev_seq.push(seq_vertices[j - 1].clone());
            }
            for j in ((jmin + 1)..=the_length).rev() {
                prev_seq.push(seq_vertices[j - 1].clone());
            }
        }

        new_wire.0.set_closed(true);
        new_wire.0.set_orientation(Orientation::Forward);
        cw.my_work[i - 1] = new_wire.0;

        prev_bary = cur_bary;
    }

    Ok(())
}

/// `PrevSeq` / `PrevEseq` (`cxx:1876-1881`): one `WireExplorer` pass collecting
/// `CurrentVertex()` and `Current()` in step.
fn wire_sequence(w: &Wire) -> (Vec<Vertex>, Vec<Edge>) {
    let mut vertices = Vec::new();
    let mut edges = Vec::new();
    let mut exp = WireExplorer::new();
    exp.init(w, None);
    while exp.more() {
        if let Some(v) = exp.current_vertex() {
            vertices.push(v.clone());
        }
        if let Some(e) = exp.current_edge() {
            edges.push(e);
        }
        exp.next();
    }
    (vertices, edges)
}

/// One term of the `cxx:1930-2073` distance sums: the vertex distance of the
/// `n`-th previous sample to `SeqVertices(k)`, plus the `NbSamples - 1`
/// intermediate curve samples when `NbSamples > 0`.
///
/// `k` is the 1-based index into `SeqVertices`; the paired current edge is
/// passed separately because the backward loops sample `SeqEdges(k-1)` instead
/// of `SeqEdges(k)` (`cxx:2011-2016`, `cxx:2049`), wrapping `k = 1` to the last
/// edge. `cur_reverse` selects the backward sample direction (`cxx:2030`).
#[allow(clippy::too_many_arguments)]
fn add_pair(
    sum: &mut f64,
    prev_seq: &[Vertex],
    prev_eseq: &[Edge],
    seq_vertices: &[Vertex],
    cur_edge: Option<&Edge>,
    k: usize,
    n: usize,
    offset: &GpVec,
    nb_samples: usize,
    cur_reverse: bool,
) {
    let (Some(p_prev), Some(v_cur)) = (prev_seq.get(n - 1), seq_vertices.get(k - 1)) else {
        return;
    };
    let p_prev = BRepTool::vertex_point_world(p_prev);
    let p_cur = BRepTool::vertex_point_world(v_cur).translated_vec(offset);
    *sum += p_prev.distance(&p_cur);

    if nb_samples == 0 {
        return;
    }
    let (Some(prev_edge), Some(cur_edge)) = (prev_eseq.get(n - 1), cur_edge) else {
        return;
    };
    let (pf, pl) = BRepTool::edge_parameters(prev_edge);
    let (cf, cl) = BRepTool::edge_parameters(cur_edge);
    let sample_prev = (pl - pf) / nb_samples as f64;
    let sample_cur = (cl - cf) / nb_samples as f64;
    let (Some(prev_curve), Some(cur_curve)) = (
        BRepTool::edge_curve_adaptor_world(prev_edge),
        BRepTool::edge_curve_adaptor_world(cur_edge),
    ) else {
        return;
    };
    let prev_forward = prev_edge.orientation() == Orientation::Forward;
    let cur_forward = cur_edge.orientation() == Orientation::Forward;
    for nbs in 1..nb_samples {
        let n_f = nbs as f64;
        let m_f = (nb_samples - nbs) as f64;
        let par_prev = if prev_forward {
            pf + n_f * sample_prev
        } else {
            pf + m_f * sample_prev
        };
        // `cxx:1952-1957` / `cxx:2029-2034`.
        let par_cur = match (cur_reverse, cur_forward) {
            (false, true) => cf + n_f * sample_cur,
            (false, false) => cf + m_f * sample_cur,
            (true, true) => cf + m_f * sample_cur,
            (true, false) => cf + n_f * sample_cur,
        };
        let pon_prev = curve_value(&prev_curve, par_prev);
        let pon_cur = curve_value(&cur_curve, par_cur).translated_vec(offset);
        *sum += pon_prev.distance(&pon_cur);
    }
}

/// `BRepAdaptor_Curve::Value(U)`.
fn curve_value(curve: &std::sync::Arc<dyn Curve>, u: f64) -> GpPnt {
    curve.d0(u)
}
