//! Port of ShapeFix_ComposeShell::CollectWires (ShapeFix_ComposeShell.cxx:2512-2936).
//!
//! The accumulated wire is a raw ShapeExtend_WireData, represented here as a
//! Vec of Edge (exactly OCCT's ordered edge list), and only wrapped into a
//! WireSegment when it is appended to the result.

use occt_core::gp::{GpPnt2d, GpVec2d};
use occt_core::precision::ANGULAR;

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Vertex};

use super::helpers::*;
use super::shell::{ComposeShell, SHAPEEXTEND_FAIL2, SHAPEEXTEND_FAIL5};
use super::wire_data::reverse_wire_data_on_face;
use super::wire_segment::{WireSegment, MAXIND, MININD};

/// `TopoDS_Shape::IsSame` on two optional vertices (a null shape is not same
/// as anything).
fn same_v(a: &Option<Vertex>, b: &Option<Vertex>) -> bool {
    match (a, b) {
        (Some(x), Some(y)) => crate::topo_tools_full::is_same(&x.0, &y.0),
        _ => false,
    }
}

fn same_e(a: &Option<Edge>, b: &Option<Edge>) -> bool {
    match (a, b) {
        (Some(x), Some(y)) => crate::topo_tools_full::is_same(&x.0, &y.0),
        _ => false,
    }
}

impl ComposeShell {
    /// ShapeFix_ComposeShell::CollectWires (cxx:2512-2936): connect the wire
    /// segments of `seqw` into closed wires, appending them to `wires`.
    pub fn collect_wires(&mut self, wires: &mut Vec<WireSegment>, seqw: &mut Vec<WireSegment>) {
        let Some(face) = self.face.clone() else {
            return;
        };
        let nu = self.grid.nb_u_patches() as i32;
        let nv = self.grid.nb_v_patches() as i32;
        // cxx:2518: shorts is allocated on the initial seqw length; it is NOT
        // shifted when an element is removed (cxx:2524-2526), and the fill
        // order keeps it aligned with the shrunk sequence.
        let mut shorts: Vec<i32> = vec![0; seqw.len()];

        // cxx:2519-2549: move vertex / INTERNAL segments out, measure the rest.
        let mut i = 0usize;
        while i < seqw.len() {
            if seqw[i].is_vertex() || seqw[i].orientation() == Orientation::Internal {
                wires.push(seqw[i].clone());
                seqw.remove(i);
                continue; // cxx:2525-2526: i-- then i++ -> same i.
            }
            let isshort = is_short_segment(
                &seqw[i],
                &face,
                &self.grid,
                self.u_resolution,
                self.v_resolution,
            );
            shorts[i] = isshort;
            let one_degenerated = seqw[i].nb_edges() == 1
                && seqw[i].edge(1).map(BRepTool::is_degenerated).unwrap_or(false);
            if isshort > 0
                && (seqw[i].orientation() == Orientation::External || one_degenerated)
            {
                seqw[i].set_orientation(Orientation::Internal);
            }
            i += 1;
        }

        // cxx:2551-2567: connection state.
        let mut sbwd: Vec<Edge> = Vec::new(); // the accumulated WireData.
        let mut has_sbwd = false;
        let mut first_v: Option<Vertex> = None;
        let mut end_v: Option<Vertex> = None;
        let mut first_edge: Option<Edge> = None;
        let mut last_edge: Option<Edge> = None;
        let mut end_pnt = GpPnt2d::new(0.0, 0.0);
        let mut first_pnt = GpPnt2d::new(0.0, 0.0);
        let mut end_tan = GpVec2d::new(0.0, 0.0);
        let mut first_tan = GpVec2d::new(0.0, 0.0);
        let mut tol = 0.0f64;
        let (mut iumin, mut iumax, mut ivmin, mut ivmax) = (0i32, 0i32, 0i32, 0i32);
        let mut can_be_closed = false;
        let (mut dsu, mut dsv) = (0.0f64, 0.0f64);

        loop {
            // cxx:2562-2567: candidate search state.
            let mut index: Option<usize> = None;
            let mut misoriented = true;
            let mut samepatch = false;
            let mut reverse = false;
            let mut connected = false;
            let mut angle = -std::f64::consts::PI;
            let mut mindist = f64::MAX; // RealLast()
            let mut weigth = 0i32;
            let (mut shiftu, mut shiftv) = (0.0f64, 0.0f64);

            // cxx:2570-2724: find the next segment to connect.
            for i in 0..seqw.len() {
                let seg = &seqw[i];
                if seg.is_vertex() {
                    continue;
                }
                let an_or = seg.orientation();
                if an_or == Orientation::Internal {
                    continue;
                }
                if !has_sbwd {
                    // cxx:2584-2604: for the first segment, take any.
                    if shorts[i] > 0 || an_or == Orientation::External {
                        continue;
                    }
                    if an_or == Orientation::Forward {
                        reverse = true;
                    }
                    index = Some(i);
                    if let Some(p) = seg.get_patch_index(1) {
                        iumin = p.0;
                        iumax = p.1;
                        ivmin = p.2;
                        ivmax = p.3;
                    }
                    misoriented = false;
                    dsu = 0.0;
                    dsv = 0.0;
                    break;
                }
                // cxx:2607-2614: same-patch test and priority.
                let sp = is_same_patch(&seqw[i], nu, nv, &mut iumin, &mut iumax, &mut ivmin, &mut ivmax, false);
                if !sp && (can_be_closed || (index.is_some() && samepatch)) {
                    continue;
                }
                let wire_edges: Vec<Edge> = seg.edges().to_vec();
                for j in 0..2usize {
                    let candidate_v = if j == 1 { seg.last_vertex() } else { seg.first_vertex() };
                    if !same_v(&end_v, &candidate_v) {
                        continue; // cxx:2629-2632
                    }
                    let misor = an_or == if j == 1 { Orientation::Reversed } else { Orientation::Forward };
                    // cxx:2639-2652: returning by the same edge is lowest priority.
                    let back_edge = if j == 1 {
                        wire_edges[wire_edges.len() - 1].clone()
                    } else {
                        wire_edges[0].clone()
                    };
                    if same_e(&last_edge, &Some(back_edge)) {
                        if index.is_none() && !can_be_closed {
                            index = Some(i);
                            reverse = j != 0;
                            connected = true;
                            misoriented = misor;
                            samepatch = sp;
                            weigth = (if sp { 16 } else { 0 })
                                + (if connected { 8 } else { 0 })
                                + (if !misor { 4 } else { 0 });
                            dsu = 0.0;
                            dsv = 0.0;
                        }
                        continue;
                    }
                    // cxx:2654-2673: starting tangent.
                    let mut l_pnt = GpPnt2d::new(0.0, 0.0);
                    let mut l_vec = GpVec2d::new(0.0, 0.0);
                    let mut edge_tol = 0.0f64;
                    let mut k = 0usize;
                    while k < wire_edges.len() {
                        let e = if j == 1 {
                            let mut t = wire_edges[wire_edges.len() - 1 - k].clone();
                            t.0.reverse();
                            t
                        } else {
                            wire_edges[k].clone()
                        };
                        edge_tol = BRepTool::edge_tolerance(&e);
                        if get_end_tangent_2d(&e, &face, false, &mut l_pnt, &mut l_vec, 1e-3) {
                            break;
                        }
                        k += 1;
                    }
                    if k >= wire_edges.len() {
                        self.status |= SHAPEEXTEND_FAIL2; // cxx:2670-2673
                    }
                    // cxx:2675-2687: closed-mode period shift.
                    if self.closed_mode {
                        if self.u_closed {
                            shiftu = crate::shhealing::adjust_by_period(l_pnt.x(), end_pnt.x(), self.u_period);
                            l_pnt.set_x(l_pnt.x() + shiftu);
                        }
                        if self.v_closed {
                            shiftv = crate::shhealing::adjust_by_period(l_pnt.y(), end_pnt.y(), self.v_period);
                            l_pnt.set_y(l_pnt.y() + shiftv);
                        }
                    }
                    // cxx:2689-2694.
                    let mut ang = if shorts[i] > 0 { std::f64::consts::PI } else { end_tan.angle(&l_vec) };
                    if self.closed_mode && shorts[i] <= 0 && std::f64::consts::PI - ang < ANGULAR {
                        ang = 0.0;
                    }
                    // cxx:2696-2722.
                    let ctol = edge_tol.max(end_v.as_ref().map(BRepTool::vertex_tolerance).unwrap_or(0.0));
                    let conn = is_coincided(&end_pnt, &l_pnt, self.u_resolution, self.v_resolution, ctol);
                    let dist = end_pnt.square_distance(&l_pnt);
                    let w1 = (if sp { 16 } else { 0 })
                        + (if conn { 4 } else { 0 })
                        + (if !misor { 8 } else { 0 });
                    let tail1 = (if !conn && dist < mindist { 2 } else { 0 })
                        + (if ang > angle { 1 } else { 0 });
                    let tail2 = (if !connected && dist > mindist { 2 } else { 0 })
                        + (if ang < angle { 1 } else { 0 });
                    if w1 + tail1 <= weigth + tail2 {
                        continue;
                    }
                    index = Some(i);
                    reverse = j != 0;
                    angle = ang;
                    mindist = dist;
                    connected = conn;
                    misoriented = misor;
                    samepatch = sp;
                    weigth = w1;
                    dsu = shiftu;
                    dsv = shiftv;
                }
            }

            // cxx:2726-2779: connect the found segment.
            if let Some(idx) = index {
                if misoriented {
                    self.invert_edge_status = true; // cxx:2729-2732
                }
                let mut seg = seqw[idx].clone();
                if !has_sbwd {
                    sbwd.clear();
                    has_sbwd = true;
                    // cxx:2733-2737: sbwd = new ShapeExtend_WireData.
                } else if samepatch {
                    // cxx:2738-2748: extend the patch indices.
                    is_same_patch(&seg, nu, nv, &mut iumin, &mut iumax, &mut ivmin, &mut ivmax, true);
                }
                if self.closed_mode {
                    // cxx:2752-2755.
                    if let Some(p) = seg.get_patch_index(1) {
                        iumin = p.0;
                        iumax = p.1;
                        ivmin = p.2;
                        ivmax = p.3;
                    }
                }
                let seg_external = seg.orientation() == Orientation::External;
                if !reverse {
                    sbwd.extend(seg.edges().iter().cloned()); // cxx:2760
                } else {
                    // cxx:2762-2768: wire->Add then WireData::Reverse(myFace)
                    // (cxx:483-572: Reverse + ComputeSeams + SwapSeam).
                    let mut wire = seg.edges().to_vec();
                    reverse_wire_data_on_face(&mut wire, &face);
                    sbwd.extend(wire);
                }
                if seg_external {
                    seg.set_orientation(if reverse { Orientation::Reversed } else { Orientation::Forward });
                } else {
                    seg.set_orientation(Orientation::Internal);
                }
                seqw[idx] = seg; // cxx:2778
            } else if !has_sbwd {
                break; // cxx:2780-2783: no free segments left.
            }

            // cxx:2784-2791: for the first segment, remember the start point.
            if end_v.is_none() {
                first_edge = Some(sbwd[0].clone());
                first_v = crate::shhealing::first_vertex(&sbwd[0]);
                let _ = get_end_tangent_2d(&sbwd[0], &face, false, &mut first_pnt, &mut first_tan, 1e-3);
            }

            // cxx:2793-2822: update the last edge / vertex (not for short segs).
            let doupdate = index.is_some() && (shorts[index.unwrap()] <= 0 || end_v.is_none());
            if doupdate {
                let nbe = sbwd.len();
                last_edge = Some(sbwd[nbe - 1].clone());
                end_v = crate::shhealing::last_vertex(&sbwd[nbe - 1]);
                tol = end_v.as_ref().map(BRepTool::vertex_tolerance).unwrap_or(0.0);
                let mut k = nbe;
                let mut found = false;
                while k >= 1 {
                    if get_end_tangent_2d(&sbwd[k - 1], &face, true, &mut end_pnt, &mut end_tan, 1e-3) {
                        found = true;
                        break;
                    }
                    k -= 1;
                }
                if !found {
                    self.status |= SHAPEEXTEND_FAIL2; // cxx:2810-2813
                }
                if self.u_closed {
                    end_pnt.set_x(end_pnt.x() + dsu);
                }
                if self.v_closed {
                    end_pnt.set_y(end_pnt.y() + dsv);
                }
            }

            // cxx:2824-2846: if closed or no next segment found, append to wires.
            can_be_closed = same_v(&end_v, &first_v);
            let close = index.is_none()
                || (can_be_closed
                    && !same_e(&last_edge, &first_edge)
                    && is_coincided(&end_pnt, &first_pnt, self.u_resolution, self.v_resolution, 2.0 * tol));
            if close {
                let first_v_of_first = first_edge.as_ref().and_then(crate::shhealing::first_vertex);
                if !same_v(&end_v, &first_v_of_first) {
                    self.status |= SHAPEEXTEND_FAIL5; // cxx:2830-2836
                }
                let mut s = WireSegment::with_edges(sbwd.clone(), Orientation::Forward);
                s.define_iu_min(1, iumin);
                s.define_iu_max(1, iumax);
                s.define_iv_min(1, ivmin);
                s.define_iv_max(1, ivmax);
                wires.push(s);
                sbwd.clear();
                has_sbwd = false;
                end_v = None;
                can_be_closed = false;
            }
        }

        // cxx:2853-2935: merge short 3d segments into other wires.
        for i in 0..seqw.len() {
            if shorts[i] != 1
                || seqw[i].is_vertex()
                || seqw[i].orientation() == Orientation::Internal
                || seqw[i].orientation() == Orientation::External
            {
                continue;
            }
            let wd: Vec<Edge> = seqw[i].edges().to_vec();
            let v = seqw[i].first_vertex();
            let mut minj: Option<usize> = None;
            let mut mink = 0usize;
            let mut p2d = GpPnt2d::new(0.0, 0.0);
            let mut tangent = GpVec2d::new(0.0, 0.0);
            let mut mindist = 0.0f64;
            let mut samepatch = false;
            if let Some(p) = seqw[i].get_patch_index(1) {
                iumin = p.0;
                iumax = p.1;
                ivmin = p.2;
                ivmax = p.3;
            }
            let _ = get_end_tangent_2d(&wd[0], &face, false, &mut p2d, &mut tangent, 0.0);
            for j in 0..wires.len() {
                let cand: Vec<Edge> = wires[j].edges().to_vec();
                for k in 0..cand.len() {
                    let first_ok = same_v(&v, &crate::shhealing::first_vertex(&cand[k]));
                    if !first_ok {
                        continue; // cxx:2880
                    }
                    let sp = is_same_patch(&wires[j], nu, nv, &mut iumin, &mut iumax, &mut ivmin, &mut ivmax, false);
                    if samepatch && !sp {
                        continue;
                    }
                    let mut pp = GpPnt2d::new(0.0, 0.0);
                    let _ = get_end_tangent_2d(&cand[k], &face, false, &mut pp, &mut tangent, 0.0);
                    let dist = pp.square_distance(&p2d);
                    if sp && !samepatch {
                        minj = Some(j);
                        mink = k + 1;
                        mindist = dist;
                        samepatch = sp;
                    } else if minj.is_none() || mindist > dist {
                        minj = Some(j);
                        mink = k + 1;
                        mindist = dist;
                        samepatch = sp;
                    }
                }
            }
            let Some(target_idx) = minj else {
                // cxx:2914-2923: keep it as a separate wire.
                wires.push(WireSegment::with_edges(wd.clone(), Orientation::Forward));
                continue;
            };
            // cxx:2925-2931: sbwd->Add(wd->Edge(n), mink++).
            let target = &mut wires[target_idx];
            let mut m = mink;
            for e in wd.iter() {
                target.add_edge_patch(m, e.clone(), MININD, MAXIND, MININD, MAXIND);
                m += 1;
            }
        }
    }
}
