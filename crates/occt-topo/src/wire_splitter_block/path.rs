use super::prelude::*;
use super::SmartMap;
use super::angles::*;
use super::vertices::*;

pub(super) fn nothing_to_do(map: &SmartMap, edges: &[Edge]) -> bool {
    for infos in &map.infos {
        let mut cin = 0;
        let mut cout = 0;
        for ei in infos {
            if ei.is_in {
                cin += 1;
            } else {
                cout += 1;
            }
        }
        if cin != 1 || cout != 1 {
            return false;
        }
    }
    let mut n_by_tshape: HashMap<usize, usize> = HashMap::new();
    for e in edges {
        *n_by_tshape
            .entry(GeometryRegistry::shape_key(&e.0))
            .or_insert(0) += 1;
    }
    n_by_tshape.values().all(|&n| n == 1)
}

pub(super) fn path(
    face: &Face,
    vert_closed: &HashMap<usize, bool>,
    start_slot: usize,
    start_info: usize,
    map: &mut SmartMap,
    loops: &mut Vec<TopoShape>,
) {
    let mut ls: Vec<Edge> = Vec::new();
    let mut vert_va: Vec<TopoShape> = Vec::new();
    let mut coord_va: Vec<GpPnt2d> = Vec::new();
    let mut info_seq: Vec<(usize, usize)> = Vec::new();

    let mut va_slot = start_slot;
    let mut info_idx = start_info;
    // `BOPAlgo_WireSplitter_1.cxx:518` (`aEOuta = aLS.Last()`): after the
    // closing scan cuts the path back, the incoming edge is the truncated
    // path's last edge; `None` means "the edge appended in this iteration".
    let mut e_in: Option<Edge>;
    let eps = f64::EPSILON;
    loop {
        if ls.len() == 1 && is_same(&ls[0].0, &map.infos[va_slot][info_idx].edge.0) {
            return;
        }
        map.infos[va_slot][info_idx].passed = true;
        // Only the closing scan below may rewrite the incoming edge.
        e_in = None;
        let e_out = map.infos[va_slot][info_idx].edge.clone();
        let va = map.verts[va_slot].clone();
        ls.push(e_out.clone());
        vert_va.push(va.clone());
        info_seq.push((va_slot, info_idx));

        let mut p_va = va.clone();
        p_va.set_orientation(Orientation::Forward);
        let pa = coord2d(&p_va, &e_out, face);
        coord_va.push(pa);

        let vb = next_vertex(&p_va, &e_out);
        let pb = coord2d(&vb, &e_out, face);
        if trace() {
            eprintln!(
                "[sbp] stp{} e={} or={:?} va={} pa=({},{}) -> vb={} pb=({},{})",
                ls.len(),
                sk(&e_out.0),
                e_out.orientation(),
                sk(&va),
                pa.x(),
                pa.y(),
                sk(&vb),
                pb.x(),
                pb.y()
            );
        }
        let Some(&vb_slot) = map.index.get(&GeometryRegistry::shape_key(&vb)) else {
            return;
        };
        let tol2d = 2.0 * tolerance_2d(&Vertex(vb.clone()), face);
        let tol2d2 = tol2d * tol2d;
        let closed = vert_closed
            .get(&GeometryRegistry::shape_key(&vb))
            .copied()
            .unwrap_or(false);

        {
            let mut buf: Vec<Edge> = Vec::new();
            let mut has_edge = false;
            let nb = ls.len();
            let mut cut: Option<usize> = None;
            for i in (0..nb).rev() {
                let e_prev = &ls[i];
                buf.push(e_prev.clone());
                if !has_edge {
                    has_edge = !BRepTool::is_degenerated(e_prev);
                    if !has_edge {
                        continue;
                    }
                }
                let same_v = is_same(&vert_va[i], &vb);
                let mut same_v2d = same_v;
                if same_v && closed {
                    same_v2d = coord_va[i].square_distance(&pb) < tol2d2;
                    if same_v2d {
                        let ud = (coord_va[i].x() - pb.x()).abs();
                        let vd = (coord_va[i].y() - pb.y()).abs();
                        let (tu, tv) = uv_tolerance_2d(&Vertex(vb.clone()), face);
                        if ud > 2.0 * tu || vd > 2.0 * tv {
                            same_v2d = false;
                        }
                    }
                }
                if same_v && same_v2d {
                    if trace() {
                        eprintln!("[sbp]   close i={} d2={} tol2={}", i + 1, coord_va[i].square_distance(&pb), tol2d2);
                    }
                    let mut priz = true;
                    if buf.len() == 2 && is_same(&buf[0].0, &buf[1].0) {
                        priz = false;
                    }
                    if priz {
                        loops.push(make_wire(&buf));
                    }
                    if i < 1 {
                        return;
                    }
                    cut = Some(i);
                    break;
                }
            }
            if let Some(i) = cut {
                ls.truncate(i);
                vert_va.truncate(i);
                coord_va.truncate(i);
                info_seq.truncate(i);
                // `BOPAlgo_WireSplitter_1.cxx:518`: the incoming edge for the
                // next choice is the truncated path's last edge, not the edge
                // that just closed the loop.
                e_in = ls.last().cloned();
            }
        }

        let is_boundary = info_seq
            .last()
            .map(|&(s, ii)| !map.infos[s][ii].is_inside)
            .unwrap_or(true);
        let le = &map.infos[vb_slot];
        let angle_in = angle_in(e_in.as_ref().unwrap_or(&e_out), le);
        let i_cnt = nb_ways_out(le);
        let mut min_angle = 100.0;
        let mut chosen: Option<usize> = None;
        let mut only_inside: Option<usize> = None;
        let mut n_inside = 0usize;

        for (k, ei) in le.iter().enumerate() {
            if ei.is_in || ei.passed {
                continue;
            }
            if i_cnt == 0 {
                return;
            }
            if i_cnt == 1 {
                chosen = Some(k);
                break;
            }
            let ang = if is_same(&ei.edge.0, &e_out.0) {
                TWO_PI
            } else {
                if closed {
                    let p2 = coord2d_vf(&ei.edge, face);
                    if trace() {
                        eprintln!(
                            "[sbp]   vfcheck e={} or={:?} p2=({},{}) pb=({},{}) d2={} tol2={}",
                            sk(&ei.edge.0),
                            ei.edge.orientation(),
                            p2.x(),
                            p2.y(),
                            pb.x(),
                            pb.y(),
                            p2.square_distance(&pb),
                            tol2d2
                        );
                    }
                    if p2.square_distance(&pb) > tol2d2 {
                        continue;
                    }
                }
                clockwise_angle(angle_in, ei.angle)
            };
            if is_boundary && ei.is_inside {
                n_inside += 1;
                only_inside = Some(k);
            }
            if trace() {
                eprintln!(
                    "[sbp]   at {} cand e={} in={} ins={} ang={} cw={} cw_in={} ang_in={}",
                    sk(&vb),
                    sk(&ei.edge.0),
                    ei.is_in as u8,
                    ei.is_inside as u8,
                    ei.angle,
                    ang,
                    clockwise_angle(angle_in, ei.angle),
                    angle_in
                );
            }
            if ang < min_angle - eps {
                min_angle = ang;
                chosen = Some(k);
            }
        }
        if n_inside == 1 {
            chosen = only_inside;
        }
        let Some(k) = chosen else {
            return;
        };
        va_slot = vb_slot;
        info_idx = k;
    }
}

pub(super) fn make_wire(edges: &[Edge]) -> TopoShape {
    let w = TopoBuilder::new().make_wire(edges);
    w.0.set_closed(true);
    w.0
}

