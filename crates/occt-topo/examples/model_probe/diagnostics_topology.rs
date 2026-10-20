//! Topology diagnostics (--bndinf, --occ, --roots, --volroots, --wirehist).
use super::tname;
use super::Ctx;
use occt_topo::brep_bnd_lib;
use occt_topo::brep_tool::BRepTool;
use occt_topo::meshing::range_splitter::mesh_surface_type;
use occt_topo::tgeometry::GeometryRegistry;
use occt_topo::topo_tools::{edges_of_wire, faces_of, wires_of_face};

pub(crate) fn run(ctx: &Ctx) -> bool {
    let model = &ctx.model;
    let shape = ctx.shape.clone();
    if std::env::args().any(|a| a == "--bndinf") {
        // TEMP: which faces drive the document box open, and with what UV window.
        let doc = occt_topo::brep_exchange::document_shape(&model).expect("document shape");
        let mut nodes_open = 0usize;
        for (ni, s) in model.shapes.iter().enumerate() {
            let nb = brep_bnd_lib::shape_bnd_box(&s.shape);
            if !nb.is_finite() {
                nodes_open += 1;
                println!(
                    "BNDNODE n={ni} name={} open=({},{},{},{},{},{})",
                    s.name,
                    nb.is_open_xmin(),
                    nb.is_open_xmax(),
                    nb.is_open_ymin(),
                    nb.is_open_ymax(),
                    nb.is_open_zmin(),
                    nb.is_open_zmax()
                );
            }
        }
        println!("BNDNODE open_nodes={nodes_open} total={}", model.shapes.len());
        // TEMP: faces sharing one registered surface pointer (the pcurve
        // representation key), so a bad face's rep key can be traced to a
        // sibling face.
        let all_faces = faces_of(&doc);
        let mut surface_owners: std::collections::HashMap<usize, Vec<usize>> =
            std::collections::HashMap::new();
        for (i, f) in all_faces.iter().enumerate() {
            if let Some(s) = BRepTool::face_surface(f) {
                surface_owners
                    .entry(std::sync::Arc::as_ptr(&s) as *const () as usize)
                    .or_default()
                    .push(i);
            }
        }
        for (ptr, owners) in surface_owners.iter() {
            if owners.len() > 1 {
                println!("BNDSHARE surf={ptr:x} faces={owners:?}");
            }
        }
        println!("BNDSHARE shared_surfaces={}", surface_owners.iter().filter(|(_, v)| v.len() > 1).count());
        let mut bad = 0usize;
        for (i, f) in all_faces.iter().enumerate() {
            let bb = brep_bnd_lib::shape_bnd_box(&f.0);
            if bb.is_finite() {
                continue;
            }
            bad += 1;
            let surf = BRepTool::face_surface(f);
            let st = surf
                .as_ref()
                .map(|s| tname(mesh_surface_type(s.as_ref())))
                .unwrap_or("none");
            let (u1, u2, v1, v2) = occt_topo::brep_tools::uv_bounds(f);
            let (sur, svr) = surf
                .as_ref()
                .map(|s| (s.u_range(), s.v_range()))
                .unwrap_or(((0.0, 0.0), (0.0, 0.0)));
            let (b0, b1) = (bb.corner_min(), bb.corner_max());
            println!(
                "BNDINF face={} type={} uv=({:.6},{:.6},{:.6},{:.6}) ur={:?} vr={:?} open=({},{},{},{},{},{}) bbox=({:.6},{:.6},{:.6})-({:.6},{:.6},{:.6}) wires={}",
                i,
                st,
                u1,
                u2,
                v1,
                v2,
                sur,
                svr,
                bb.is_open_xmin(),
                bb.is_open_xmax(),
                bb.is_open_ymin(),
                bb.is_open_ymax(),
                bb.is_open_zmin(),
                bb.is_open_zmax(),
                b0.x(),
                b0.y(),
                b0.z(),
                b1.x(),
                b1.y(),
                b1.z(),
                wires_of_face(f).len()
            );
            let reg = GeometryRegistry::global();
            let key = GeometryRegistry::shape_key(&f.0);
            if let Some(s) = BRepTool::face_surface(f) {
                println!(
                    "BNDINF  facesurf={:x} kind={} ur={:?} vr={:?}",
                    std::sync::Arc::as_ptr(&s) as *const () as usize,
                    tname(mesh_surface_type(s.as_ref())),
                    s.u_range(),
                    s.v_range()
                );
            }
            for (wi, w) in wires_of_face(f).iter().enumerate() {
                let es = edges_of_wire(w);
                println!("BNDINF  wire[{wi}] ori={:?} nEdges={}", w.0.orientation(), es.len());
                for (ei, e) in es.iter().enumerate() {
                    let pcs = reg.edge_pcurves(&e.0, key);
                    let range = reg.pcurve_range(&e.0, key);
                    let reps = reg.edge_pcurve_reps(&e.0);
                    let rep_keys: Vec<String> = reps
                        .iter()
                        .map(|(s, _, _, _)| {
                            format!("{:x}", std::sync::Arc::as_ptr(s) as *const () as usize)
                        })
                        .collect();
                    println!(
                        "BNDINF   e[{ei}] ori={:?} deg={} deg_face={} npc={} rng={:?} onsurf={} reps={:?} facekey={key:x} reprkey={:x}",
                        e.0.orientation(),
                        BRepTool::is_degenerated(e),
                        BRepTool::is_closed_edge_face(e, f),
                        pcs.len(),
                        range,
                        occt_topo::boptools_2d::curve_on_surface_range(e, f).is_some(),
                        rep_keys,
                        reg.debug_repr_key(key)
                    );
                    for (rs, _, rf, rl) in &reps {
                        let ptr = std::sync::Arc::as_ptr(rs) as *const () as usize;
                        let owners = surface_owners
                            .get(&ptr)
                            .map(|v| format!("{v:?}"))
                            .unwrap_or_else(|| "none".to_string());
                        println!(
                            "BNDINF    rep surf={:x} kind={} ur={:?} vr={:?} rng=({rf:.6},{rl:.6}) owners={owners}",
                            ptr,
                            tname(mesh_surface_type(rs.as_ref())),
                            rs.u_range(),
                            rs.v_range()
                        );
                    }
                }
            }
        }
        let total = brep_bnd_lib::shape_bnd_box(&doc);
        match total.get() {
            Some((x0, x1, y0, y1, z0, z1)) => println!(
                "BNDINF get=[{x0:.6},{y0:.6},{z0:.6}]~[{x1:.6},{y1:.6},{z1:.6}]"
            ),
            None => println!("BNDINF get=none"),
        }
        println!("BNDINF bad_faces={bad} doc_open=({},{},{},{},{},{})", total.is_open_xmin(), total.is_open_xmax(), total.is_open_ymin(), total.is_open_ymax(), total.is_open_zmin(), total.is_open_zmax());
        return true;
    }
    if std::env::args().any(|a| a == "--occ") {
        use occt_topo::iterator::ShapeExplorer;
        let occ = ShapeExplorer::new(vec![shape.clone()], occt_topo::abs::ShapeType::Face).count();
        println!("OCC faces_distinct={} faces_occ={}", faces_of(&shape).len(), occ);
        return true;
    }
    if std::env::args().any(|a| a == "--roots") {
        println!("ROOTS count={}", model.shapes.len());
        for (i, s) in model.shapes.iter().enumerate() {
            let ptr = std::sync::Arc::as_ptr(&s.shape.tshape) as usize;
            println!("ROOT i={i} ptr={ptr:#x} type={:?} faces={}", s.shape.shape_type(), faces_of(&s.shape).len());
        }
        return true;
    }
    if std::env::args().any(|a| a == "--volroots") {
        println!("VOLROOT count={}", model.shapes.len());
        for (i, s) in model.shapes.iter().enumerate() {
            let ptr = std::sync::Arc::as_ptr(&s.shape.tshape) as usize;
            match occt_topo::brep_gprop_full::volume_properties(&s.shape) {
                Ok(props) => println!("VOLROOT i={i} ptr={ptr:#x} type={:?} faces={} mass={:.6}", s.shape.shape_type(), faces_of(&s.shape).len(), props.mass()),
                Err(e) => println!("VOLROOT i={i} ptr={ptr:#x} err={e}"),
            }
        }
        return true;
    }
    if std::env::args().any(|a| a == "--voltotal") {
        match occt_topo::brep_gprop_full::volume_properties(&shape) {
            Ok(props) => println!("VOLTOTAL mass={:.6}", props.mass()),
            Err(e) => println!("VOLTOTAL err={e}"),
        }
        return true;
    }
    let mode_facevol = std::env::args().any(|a| a == "--facevol");
    if mode_facevol {
        for (i, f) in faces_of(&shape).iter().enumerate() {
            match occt_topo::brep_gprop_full::face_volume_contribution(f) {
                Ok(v) => println!("FV face={i} contrib={v:.6}"),
                Err(e) => println!("FV face={i} err={e}"),
            }
        }
        return true;
    }
    if std::env::args().any(|a| a == "--wirehist") {
        use std::collections::BTreeMap;
        let faces = faces_of(&shape);
        let mut hist: BTreeMap<usize, usize> = BTreeMap::new();
        let mut n_wires = 0usize;
        let mut n2 = 0usize;
        for f in faces.iter() {
            let ws = wires_of_face(f);
            n_wires += ws.len();
            if ws.len() >= 2 {
                n2 += 1;
                let bb = brep_bnd_lib::shape_bnd_box(&f.0);
                let (mn, mx) = (bb.corner_min(), bb.corner_max());
                println!(
                    "MULTI bbox=({:.3},{:.3},{:.3})-({:.3},{:.3},{:.3}) wires={}",
                    mn.x(), mn.y(), mn.z(), mx.x(), mx.y(), mx.z(), ws.len()
                );
            }
            *hist.entry(ws.len()).or_insert(0) += 1;
        }
        let mut n1c = 0usize;
        let mut nf1c = 0usize;
        for f in faces.iter() {
            let mut hit = false;
            for w in wires_of_face(f).iter() {
                let es = edges_of_wire(w);
                if es.len() == 1 {
                    let fv = occt_topo::shhealing::first_vertex(&es[0]);
                    let lv = occt_topo::shhealing::last_vertex(&es[0]);
                    let closed = match (&fv, &lv) {
                        (Some(a), Some(b)) => occt_topo::topo_tools::is_same(&a.0, &b.0),
                        _ => false,
                    };
                    if closed {
                        n1c += 1;
                        hit = true;
                    }
                }
            }
            if hit {
                nf1c += 1;
            }
        }
        println!(
            "TOTAL faces={} wires={} faces_with_2plus_wires={} wires_1edge_closed={} faces_1edge_closed={}",
            faces.len(),
            n_wires,
            n2,
            n1c,
            nf1c
        );
        let s: Vec<String> = hist.iter().map(|(k, v)| format!("{k}:{v}")).collect();
        println!("WIREHIST {}", s.join(" "));
        return true;
    }
    false
}
