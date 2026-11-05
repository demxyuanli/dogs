//! `BRepFill_CompatibleWires::SameNumberByPolarMethod` (`cxx:1008-1503`).
//!
//! The closed-wire branch of `Perform`: every wire is made to carry the same
//! number of edges by cutting the wires radially. Vertices of section `i` are
//! carried into the plane of section `i-1` (`Transform`), the ray from the
//! section centre to the carried point is intersected with the lower wire
//! (`EdgeIntersectOnWire`) and both sections walk forward to the last one
//! before the edges are re-stacked in `myWork(ideb)`'s order
//! (`BuildConnectedEdges` + the `scal1` / `scal2` orientation vote).
//!
//! Deviations, all of them OCCT throw / assertion sites that cannot be
//! reproduced as such:
//!
//! * `MapVLV(VRoot)` (`cxx:1278`) on an unbound key makes
//!   `NCollection_DataMap::operator()` throw `Standard_NoSuchObject`; the
//!   append is skipped here instead of aborting the whole `Perform`.
//! * `itF.Value()` / `itL.Value()` (`cxx:1316`) and `SeqV.Value(1)`
//!   (`cxx:1100`, `cxx:1107`) on an exhausted sequence assert in OCCT; here the
//!   search fails with `myStatus = Failed`.
//! * A null `Vsol` is not appended to `RMap` / `MapVLV`. OCCT appends the null
//!   `TopoDS_Vertex`; every reader is a `SearchRoot` / `SearchVertex` `IsSame`
//!   comparison, which a null shape never matches.

use std::collections::HashMap;

use occt_core::gp::{GpPnt, GpVec};
use occt_core::precision::CONFUSION;
use occt_geom::Curve;

use crate::brep_check::{wire_closed, Status};
use crate::brep_fill_compatible_wires::statics::{
    add_new_edge, build_connected_edges, edge_intersect_on_wire, plane_of_wire, search_root,
    search_vertex, seq_of_vertices, transform, vertex_map_bind, vertex_map_get,
    vertex_map_get_mut, VertexMap,
};
use crate::brep_fill_compatible_wires::CompatibleWires;
use crate::brep_fill_generator::ThruSectionErrorStatus;
use crate::brep_lib_make_wire::MakeWire;
use crate::brep_tool::BRepTool;
use crate::brep_tools_wire_explorer::{edge_vertices_cum_ori, WireExplorer};
use crate::builder::TopoBuilder;
use crate::shape::{Edge, TopoShape, Vertex, Wire};
use crate::topo_tools::{
    edges_of_wire, is_same, shape_map_key, ShapeMapKey,
};

/// `SameNumberByPolarMethod(WithRotation)` (`cxx:1008-1503`).
pub(crate) fn same_number_by_polar_method(
    cw: &mut CompatibleWires,
    with_rotation: bool,
) -> Result<(), String> {
    let nb_sects = cw.my_work.len();
    // `NCollection_DataMap<TopoDS_Shape, NCollection_Sequence<TopoDS_Shape>>
    // EdgeNewEdges` (`cxx:1014-1015`).
    let mut edge_new_edges: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();

    // `cxx:1019-1031`: every section must be a closed wire.
    let mut all_closed = true;
    for i in 1..=nb_sects {
        let w = Wire(cw.my_work[i - 1].clone());
        all_closed = all_closed && wire_closed(&w) == Status::NoError;
    }
    if !all_closed {
        cw.my_status = ThruSectionErrorStatus::NotSameTopology;
        return Ok(());
    }

    // `cxx:1033-1041`.
    let mut ideb = 1usize;
    let mut ifin = nb_sects;
    if cw.my_degen1 {
        ideb += 1;
    }
    if cw.my_degen2 {
        ifin -= 1;
    }
    // `cxx:1042`: looping sections.
    let v_closed = (!cw.my_degen1)
        && (!cw.my_degen2)
        && ideb <= ifin
        && is_same(&cw.my_work[ideb - 1], &cw.my_work[ifin - 1]);

    // `cxx:1043-1072`: remove the degenerated edges, in section order.
    for i in ideb..=ifin {
        let w = Wire(cw.my_work[i - 1].clone());
        let has_deg_edge = edges_of_wire(&w).iter().any(BRepTool::is_degenerated);
        if has_deg_edge {
            let mut new_wire = Wire::new();
            let builder = TopoBuilder::new();
            for e in edges_of_wire(&w) {
                if !BRepTool::is_degenerated(&e) {
                    builder.add(&mut new_wire.0, &e.0);
                }
            }
            cw.my_work[i - 1] = new_wire.0;
        }
    }

    // `cxx:1075-1081`: upper bound on the possible number of cuts.
    let mut nb_max_v = 0usize;
    for i in 1..=nb_sects {
        let w = Wire(cw.my_work[i - 1].clone());
        let mut exp = WireExplorer::new();
        exp.init(&w, None);
        while exp.more() {
            nb_max_v += 1;
            exp.next();
        }
    }

    // `cxx:1085-1109`: the position / axis of every wire's plane.
    let mut pos = vec![GpPnt::zero(); nb_sects + 1];
    let mut axe = vec![GpVec::new(0.0, 0.0, 0.0); nb_sects + 1];
    for i in ideb..=ifin {
        let w = Wire(cw.my_work[i - 1].clone());
        if let Some(p) = plane_of_wire(&w) {
            pos[i] = p.location();
            // `gp_Vec(P.Axis().Direction())` (`cxx:1093`).
            let ax1 = p.axis();
            let d = ax1.direction();
            axe[i] = GpVec::new(d.x(), d.y(), d.z());
        }
    }
    if cw.my_degen1 {
        // `cxx:1096-1102`.
        let w = Wire(cw.my_work[0].clone());
        if let Some(v) = seq_of_vertices(&w).first() {
            pos[1] = BRepTool::vertex_point_world(v);
        }
        axe[1] = axe[ideb];
    }
    if cw.my_degen2 {
        // `cxx:1103-1108`.
        let w = Wire(cw.my_work[nb_sects - 1].clone());
        if let Some(v) = seq_of_vertices(&w).first() {
            pos[nb_sects] = BRepTool::vertex_point_world(v);
        }
        axe[nb_sects] = axe[ifin];
    }

    // `cxx:1113-1176`: report the vertices of wire i onto wire i-1.
    let mut rmap: VertexMap = Vec::new();
    for i in (ideb + 1..=ifin).rev() {
        let wire1 = Wire(cw.my_work[i - 1].clone());
        let seq_v = seq_of_vertices(&wire1);
        if seq_v.len() > nb_max_v {
            cw.my_status = ThruSectionErrorStatus::Failed;
            return Ok(());
        }
        for vi in &seq_v {
            vertex_map_bind(&mut rmap, vi, Vec::new());
            let pi = BRepTool::vertex_point_world(vi);
            let pnew = transform(
                with_rotation,
                &pi,
                &pos[i],
                &axe[i],
                &pos[i - 1],
                &axe[i - 1],
            );
            if pnew.distance(&pos[i - 1]) > CONFUSION {
                // `myWork(i - 1)` is read by the intersection and written back
                // when it was cut.
                let wire_prev = Wire(cw.my_work[i - 2].clone());
                let res = edge_intersect_on_wire(
                    &pos[i - 1],
                    &pnew,
                    cw.my_percent,
                    &rmap,
                    &wire_prev,
                    &mut edge_new_edges,
                );
                if res.new_vertex {
                    if let Some(nw) = res.new_wire {
                        cw.my_work[i - 2] = nw.0;
                    }
                }
                if let Some(vsol) = &res.vsol {
                    if let Some(list) = vertex_map_get_mut(&mut rmap, vi) {
                        list.push(vsol.0.clone());
                    }
                }
            }
        }
    }

    // `cxx:1179-1209`: the vertex -> vertex report chains of the first section.
    let mut map_vlv: VertexMap = Vec::new();
    let wire_ideb = Wire(cw.my_work[ideb - 1].clone());
    let seq_v = seq_of_vertices(&wire_ideb);
    let size_map = seq_v.len();
    for vi in &seq_v {
        vertex_map_bind(&mut map_vlv, vi, vec![vi.0.clone()]);
        let mut nb_v = 1usize;
        let mut v0 = vi.clone();
        let mut v1 = search_root(&rmap, &v0);
        while let Some(cur) = v1.clone() {
            if let Some(list) = vertex_map_get_mut(&mut map_vlv, vi) {
                list.push(cur.0.clone());
            }
            nb_v += 1;
            // `cxx:1198-1206`: the guard that keeps a looping section finite.
            if is_same(&cur.0, &vi.0) || nb_v >= cw.my_work.len() {
                v1 = None;
            } else {
                v0 = cur;
                v1 = search_root(&rmap, &v0);
            }
        }
    }

    // `cxx:1211-1284`: report the vertices of wire i onto wire i+1.
    for i in ideb..ifin {
        let wire1 = Wire(cw.my_work[i - 1].clone());
        let seq_v = seq_of_vertices(&wire1);
        if seq_v.len() > nb_max_v || seq_v.len() > size_map {
            cw.my_status = ThruSectionErrorStatus::Failed;
            return Ok(());
        }
        let wire2 = Wire(cw.my_work[i].clone());
        for vi in &seq_v {
            let mut vroot: Option<Vertex> = None;
            let mut intersect = true;
            if let Some(root) = search_root(&map_vlv, vi) {
                let lvi = vertex_map_get(&map_vlv, &root).unwrap_or(&[]);
                intersect = search_vertex(lvi, &wire2).is_none();
                vroot = Some(root);
            }
            if !intersect {
                continue;
            }
            let pi = BRepTool::vertex_point_world(vi);
            let pnew = transform(
                with_rotation,
                &pi,
                &pos[i],
                &axe[i],
                &pos[i + 1],
                &axe[i + 1],
            );
            if pnew.distance(&pos[i + 1]) <= CONFUSION {
                continue;
            }
            let res = edge_intersect_on_wire(
                &pos[i + 1],
                &pnew,
                cw.my_percent,
                &map_vlv,
                &wire2,
                &mut edge_new_edges,
            );
            if let Some(vsol) = &res.vsol {
                // `MapVLV(VRoot)` throws `Standard_NoSuchObject` on an unbound
                // key; the append is skipped (see the module note).
                if let Some(root) = &vroot {
                    if let Some(list) = vertex_map_get_mut(&mut map_vlv, root) {
                        list.push(vsol.0.clone());
                    }
                }
            }
            if res.new_vertex {
                if let Some(nw) = res.new_wire {
                    cw.my_work[i] = nw.0;
                }
            }
        }
    }

    // `cxx:1286-1292`: re-stack every section in the first section's order.
    let wire = Wire(cw.my_work[ideb - 1].clone());
    let mut ibout = ifin;
    if v_closed {
        ibout -= 1;
    }
    // `cxx:1294`: a plane origin for the `Ns` / `N1` / `N2` directions.
    let p0 = pos[ideb];
    for i in (ideb + 1)..=ibout {
        let mut mw = MakeWire::new();
        let mut exp = WireExplorer::new();
        exp.init(&wire, None);
        let mut ecur = match exp.current_edge() {
            Some(e) => e,
            None => {
                cw.my_status = ThruSectionErrorStatus::ProfilesInconsistent;
                return Ok(());
            }
        };
        let (mut vf, mut vl) = edge_vertices_cum_ori(&ecur);
        let mut u1 = param_of(&ecur, vf.as_ref());
        let mut u2 = param_of(&ecur, vl.as_ref());
        let mut curve = BRepTool::edge_curve_adaptor_world(&ecur);
        let mut pps = curve_value(&curve, 0.1 * (u1 + 9.0 * u2));

        // `cxx:1308-1317`: the `i`-th counterpart of the first edge's vertices.
        let list_f = vf.as_ref().and_then(|v| vertex_map_get(&map_vlv, v));
        let list_l = vl.as_ref().and_then(|v| vertex_map_get(&map_vlv, v));
        let skip = i - ideb;
        let v1 = match list_f.and_then(|l| l.get(skip)) {
            Some(s) => Vertex(s.clone()),
            None => {
                cw.my_status = ThruSectionErrorStatus::Failed;
                return Ok(());
            }
        };
        let v2 = match list_l.and_then(|l| l.get(skip)) {
            Some(s) => Vertex(s.clone()),
            None => {
                cw.my_status = ThruSectionErrorStatus::Failed;
                return Ok(());
            }
        };

        // `cxx:1318-1380`: the edge of `myWork(i)` that joins (V1, V2) with the
        // best alignment against the first section's first edge.
        let mut esol: Option<Edge> = None;
        let mut scalmax = 0.0f64;
        let target = Wire(cw.my_work[i - 1].clone());
        for es in edges_of_wire(&target) {
            let (vvf, vvl) = edge_vertices_cum_ori(&es);
            let (Some(vvf), Some(vvl)) = (vvf, vvl) else {
                continue;
            };
            let straight = (is_same(&v1.0, &vvf.0) && is_same(&v2.0, &vvl.0))
                || (is_same(&v2.0, &vvf.0) && is_same(&v1.0, &vvl.0));
            if !straight {
                continue;
            }
            let u1param = BRepTool::parameter_of_vertex(&es.0, &vvf.0);
            let u2param = BRepTool::parameter_of_vertex(&es.0, &vvl.0);
            let curve_e = BRepTool::edge_curve_adaptor_world(&es);
            let mut pp1 = curve_value(&curve_e, 0.1 * (u1param + 9.0 * u2param));
            let mut pp2 = curve_value(&curve_e, 0.1 * (9.0 * u1param + u2param));
            for rang in (ideb + 1..=i).rev() {
                pp1 = transform(with_rotation, &pp1, &pos[rang], &axe[rang], &pos[rang - 1], &axe[rang - 1]);
                pp2 = transform(with_rotation, &pp2, &pos[rang], &axe[rang], &pos[rang - 1], &axe[rang - 1]);
            }
            let ns = GpVec::from_pnts(&p0, &pps).normalized();
            let n1 = GpVec::from_pnts(&p0, &pp1).normalized();
            let n2 = GpVec::from_pnts(&p0, &pp2).normalized();
            let scal1 = n1.dot(&ns);
            if scal1 > scalmax {
                scalmax = scal1;
                esol = Some(es.clone());
            }
            let scal2 = n2.dot(&ns);
            if scal2 > scalmax {
                scalmax = scal2;
                let mut er = es.clone();
                er.0.reverse();
                esol = Some(er);
            }
        }
        let Some(esol) = esol else {
            cw.my_status = ThruSectionErrorStatus::ProfilesInconsistent;
            return Ok(());
        };
        let _ = mw.add(&esol);

        // `cxx:1382-1438`: append the connected edges in the same order.
        let connected = build_connected_edges(&target, &esol, &v2);
        let mut it_ce = connected.into_iter();
        while exp.more() {
            let Some(next_ce) = it_ce.next() else { break };
            ecur = match exp.current_edge() {
                Some(e) => e,
                None => break,
            };
            let (f, l) = edge_vertices_cum_ori(&ecur);
            vf = f;
            vl = l;
            u1 = param_of(&ecur, vf.as_ref());
            u2 = param_of(&ecur, vl.as_ref());
            curve = BRepTool::edge_curve_adaptor_world(&ecur);
            pps = curve_value(&curve, 0.1 * (u1 + 9.0 * u2));

            let mut e = Edge(next_ce.clone());
            let (vvf, vvl) = edge_vertices_cum_ori(&e);
            let (Some(vvf), Some(vvl)) = (vvf, vvl) else {
                exp.next();
                continue;
            };
            u1 = BRepTool::parameter_of_vertex(&e.0, &vvf.0);
            u2 = BRepTool::parameter_of_vertex(&e.0, &vvl.0);
            // `Curve.Initialize(E)` (`cxx:1404`): the same adaptor is reused.
            curve = BRepTool::edge_curve_adaptor_world(&e);
            let mut pp1 = curve_value(&curve, 0.1 * (u1 + 9.0 * u2));
            let mut pp2 = curve_value(&curve, 0.1 * (9.0 * u1 + u2));
            for rang in (ideb + 1..=i).rev() {
                pp1 = transform(with_rotation, &pp1, &pos[rang], &axe[rang], &pos[rang - 1], &axe[rang - 1]);
                pp2 = transform(with_rotation, &pp2, &pos[rang], &axe[rang], &pos[rang - 1], &axe[rang - 1]);
            }
            let ns = GpVec::from_pnts(&p0, &pps).normalized();
            let n1 = GpVec::from_pnts(&p0, &pp1).normalized();
            let n2 = GpVec::from_pnts(&p0, &pp2).normalized();
            let scal1 = n1.dot(&ns);
            let scal2 = n2.dot(&ns);
            if scal2 > scal1 {
                e.0.reverse();
            }
            let _ = mw.add(&e);
            exp.next();
        }
        cw.my_work[i - 1] = mw.wire().0;
    }

    // `cxx:1441-1457`: the last section of a looping set is the first one.
    if v_closed {
        let last = cw.my_work.len();
        let tail = edges_of_wire(&Wire(cw.my_work[last - 1].clone()));
        let head = edges_of_wire(&Wire(cw.my_work[0].clone()));
        for (an_edge, a_new_edge) in tail.iter().zip(head.iter()) {
            if !is_same(&an_edge.0, &a_new_edge.0) {
                edge_new_edges.insert(shape_map_key(&an_edge.0), vec![a_new_edge.0.clone()]);
            }
        }
        cw.my_work[last - 1] = cw.my_work[0].clone();
    }

    // `cxx:1459-1485`: the debug pass that also rejects an unequal cut count.
    let mut nbmax = 0usize;
    let mut nbmin = 0usize;
    for i in ideb..=ifin {
        let w = Wire(cw.my_work[i - 1].clone());
        let mut exp = WireExplorer::new();
        exp.init(&w, None);
        let mut nb_edges = 0usize;
        while exp.more() {
            nb_edges += 1;
            exp.next();
        }
        if i == ideb {
            nbmin = nb_edges;
        }
        if nbmax < nb_edges {
            nbmax = nb_edges;
        }
        if nbmin > nb_edges {
            nbmin = nb_edges;
        }
    }
    if nbmin != nbmax {
        cw.my_status = ThruSectionErrorStatus::Failed;
        return Ok(());
    }

    // `cxx:1487-1501`: expand every `myMap` entry through `EdgeNewEdges`.
    let keys: Vec<ShapeMapKey> = cw.my_map.keys().cloned().collect();
    for key in keys {
        let Some(seed) = cw.my_map.get(&key).and_then(|v| v.first().cloned()) else {
            continue;
        };
        let mut list: Vec<TopoShape> = Vec::new();
        add_new_edge(&seed, &edge_new_edges, &mut list);
        cw.my_map.insert(key, list);
    }

    Ok(())
}

/// `BRep_Tool::Parameter(V, E)`; OCCT asserts on a null vertex.
fn param_of(e: &Edge, v: Option<&Vertex>) -> f64 {
    match v {
        Some(v) => BRepTool::parameter_of_vertex(&e.0, &v.0),
        None => 0.0,
    }
}

/// `BRepAdaptor_Curve::Value(U)`; OCCT dereferences the adaptor.
fn curve_value(curve: &Option<std::sync::Arc<dyn Curve>>, u: f64) -> GpPnt {
    match curve {
        Some(c) => c.d0(u),
        None => GpPnt::zero(),
    }
}
