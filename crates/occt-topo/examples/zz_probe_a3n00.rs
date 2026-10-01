//! TEMPORARY diagnostic probe (delete after use).
use occt_topo::brep_bnd_lib;
use occt_topo::brep_exchange::prs3d_get_deflection;
use occt_topo::brep_tool::BRepTool;
use occt_topo::builder::TopoBuilder;
use occt_topo::meshing::incremental_mesh::IncrementalMesh;
use occt_topo::meshing::range_splitter::{classify_surface, SurfaceType};
use occt_topo::step::read_step_file;
use occt_topo::tgeometry::GeometryRegistry;
use occt_topo::topo_tools_full::{edge_vertices, edges_of, edges_of_wire, faces_of, vertices_of, wires_of_face};

fn tname(t: SurfaceType) -> &'static str {
    match t {
        SurfaceType::Plane => "Plane", SurfaceType::Cylinder => "Cylinder",
        SurfaceType::Cone => "Cone", SurfaceType::Sphere => "Sphere",
        SurfaceType::Torus => "Torus", SurfaceType::SurfaceOfRevolution => "Rev",
        SurfaceType::SurfaceOfExtrusion => "Extr", SurfaceType::BezierSurface => "Bezier",
        SurfaceType::BSplineSurface => "BSpline", SurfaceType::OffsetSurface => "Offset",
        SurfaceType::OtherSurface => "Other",
    }
}

fn main() {
    let path = std::env::args().nth(1).expect("usage: zz_probe_a3n00 <step> [lin_override|--edges]");
    let mode_edges = std::env::args().any(|a| a == "--edges");
    let mode_dev = std::env::args().any(|a| a == "--dev");
    let mode_low = std::env::args().any(|a| a == "--low");
    let mode_fixms = std::env::args().any(|a| a == "--fixms");
    // T-101: per-face vertex/triangle counts + bbox, so a face's mesh density can
    // be paired with OCCT's `occt_probe --mesh` `FACE k nodes=N triangles=T` line
    // by the six bbox coordinates (never by index).
    let mode_fstats = std::env::args().any(|a| a == "--fstats");
    let override_lin: Option<f64> = std::env::args().nth(2).and_then(|s| s.parse().ok());
    let model = read_step_file(&path).expect("read step");
    // Match the OCCT probe's `aReader.OneShape()`: the FIRST transfer root.
    // `model.shapes` can hold both a wrapper `Compound` and the `Solid` inside
    // it, so compounding every root would traverse the solid 2-3 times.
    let shape = model.shapes[0].shape.clone();
    if std::env::args().any(|a| a == "--occ") {
        use occt_topo::iterator::ShapeExplorer;
        let occ = ShapeExplorer::new(vec![shape.clone()], occt_topo::abs::ShapeType::Face).count();
        println!("OCC faces_distinct={} faces_occ={}", faces_of(&shape).len(), occ);
        return;
    }
    if std::env::args().any(|a| a == "--roots") {
        println!("ROOTS count={}", model.shapes.len());
        for (i, s) in model.shapes.iter().enumerate() {
            let ptr = std::sync::Arc::as_ptr(&s.shape.tshape) as usize;
            println!("ROOT i={i} ptr={ptr:#x} type={:?} faces={}", s.shape.shape_type(), faces_of(&s.shape).len());
        }
        return;
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
        return;
    }
    if std::env::args().any(|a| a == "--voltotal") {
        match occt_topo::brep_gprop_full::volume_properties(&shape) {
            Ok(props) => println!("VOLTOTAL mass={:.6}", props.mass()),
            Err(e) => println!("VOLTOTAL err={e}"),
        }
        return;
    }
    let mode_facevol = std::env::args().any(|a| a == "--facevol");
    if mode_facevol {
        for (i, f) in faces_of(&shape).iter().enumerate() {
            match occt_topo::brep_gprop_full::face_volume_contribution(f) {
                Ok(v) => println!("FV face={i} contrib={v:.6}"),
                Err(e) => println!("FV face={i} err={e}"),
            }
        }
        return;
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
                        (Some(a), Some(b)) => occt_topo::topo_tools_full::is_same(&a.0, &b.0),
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
        return;
    }
    if std::env::args().any(|a| a == "--fdump") {
        // `--fdump <idx>` dumps that face index; bare `--fdump` matches the
        // a3n00 F113 bbox (its original use).
        let want: Option<usize> = {
            let args: Vec<String> = std::env::args().collect();
            args.iter()
                .position(|a| a == "--fdump")
                .and_then(|i| args.get(i + 1))
                .and_then(|s| s.parse().ok())
        };
        for (i, f) in faces_of(&shape).iter().enumerate() {
            let bb = brep_bnd_lib::shape_bnd_box(&f.0);
            let (bmin, bmax) = (bb.corner_min(), bb.corner_max());
            match want {
                Some(w) => {
                    if i != w {
                        continue;
                    }
                }
                None => {
                    let bybox = (bmin.x() + 87.5).abs() < 1.0
                        && (bmax.x() - 87.5).abs() < 1.0
                        && (bmin.z() + 100.0).abs() < 1.0
                        && (bmax.z() + 32.0).abs() < 1.0;
                    if !bybox {
                        continue;
                    }
                }
            }
            println!(
                "FDUMP face={i} box=({:.5},{:.5},{:.5})-({:.5},{:.5},{:.5})",
                bmin.x(),
                bmin.y(),
                bmin.z(),
                bmax.x(),
                bmax.y(),
                bmax.z()
            );
            for (wi, w) in wires_of_face(f).iter().enumerate() {
                let es = edges_of_wire(w);
                println!("FDUMP  wire[{wi}] ori={:?} nEdges={}", w.0.orientation(), es.len());
                for (ei, e) in es.iter().enumerate() {
                    let fv = occt_topo::shhealing::first_vertex(e);
                    let lv = occt_topo::shhealing::last_vertex(e);
                    let p = |v: &Option<occt_topo::shape::Vertex>| match v {
                        Some(v) => {
                            let q = BRepTool::vertex_point(v);
                            format!("{:.6},{:.6},{:.6}", q.x(), q.y(), q.z())
                        }
                        None => "none".to_string(),
                    };
                    let ptr = |v: &Option<occt_topo::shape::Vertex>| match v {
                        Some(v) => format!("{:x}", std::sync::Arc::as_ptr(&v.0.tshape) as usize),
                        None => "0".to_string(),
                    };
                    println!(
                        "FDUMP   e[{ei}] ori={:?} first={}({}) last={}({}) deg={} closed={}",
                        e.0.orientation(),
                        ptr(&fv),
                        p(&fv),
                        ptr(&lv),
                        p(&lv),
                        BRepTool::is_degenerated(e),
                        BRepTool::is_closed_edge_face(e, f)
                    );
                }
            }
            return;
        }
        println!("FDUMP none");
        return;
    }
    let mode_cw = std::env::args().any(|a| a == "--checkwire");
    if mode_cw {
        use occt_topo::boptools_2d::curve_on_surface_range;
        for (i, f) in faces_of(&shape).iter().enumerate() {
            let Some(surf) = BRepTool::face_surface(f) else { continue };
            if !surf.is_v_periodic() && !surf.is_u_periodic() { continue; }
            let (uf, ul) = surf.u_range();
            let (vf, vl) = surf.v_range();
            println!("CW face={i} u=[{uf:.4},{ul:.4}] v=[{vf:.4},{vl:.4}] dU={:.4} dV={:.4}", (ul-uf).abs(), (vl-vf).abs());
            for (wi, w) in wires_of_face(f).iter().enumerate() {
                for (ei, e) in edges_of_wire(w).iter().enumerate() {
                    match curve_on_surface_range(e, f) {
                        Some((c, a, b)) => {
                            let pa = c.d0(a); let pb = c.d0(b);
                            let (e3f, e3l) = occt_topo::brep_tool::BRepTool::edge_parameters(e);
                            println!("  w{wi} e{ei} ptr={:#x} ori={:?} er=({e3f:.4},{e3l:.4}) f={a:.6} l={b:.6} df=({:.6},{:.6}) dl=({:.6},{:.6}) delta=({:.6},{:.6})",
                                std::sync::Arc::as_ptr(&e.0.tshape) as usize, e.0.orientation(), pa.x(), pa.y(), pb.x(), pb.y(), pb.x()-pa.x(), pb.y()-pa.y());
                        }
                        None => println!("  w{wi} e{ei} no pcurve"),
                    }
                }
            }
        }
        return;
    }
    if std::env::args().any(|a| a == "--raw") {
        for (i, f) in faces_of(&shape).iter().enumerate() {
            let wires = wires_of_face(f);
            if wires.len() != 2 { continue; }
            let surf = BRepTool::face_surface(f);
            let tn = surf.as_ref().map(|s| tname(classify_surface(s.as_ref()))).unwrap_or("none");
            if tn != "BSpline" { continue; }
            let face_children = f.0.tshape.read().unwrap().children.len();
            let mut s = format!("RAW face={i} {tn} face_children={face_children} wires={}", wires.len());
            for (wi, w) in wires.iter().enumerate() {
                let raw = w.0.tshape.read().unwrap().children.clone();
                s += &format!(" w{wi}.children={} [", raw.len());
                for c in &raw {
                    s += &format!("{:?}:{:#x}:{:?} ", c.shape_type(), std::sync::Arc::as_ptr(&c.tshape) as usize, c.orientation());
                }
                s += "]";
            }
            println!("{s}");
        }
        return;
    }
    if mode_fixms {
        for (i, f) in faces_of(&shape).iter().enumerate() {
            let Some(surf) = BRepTool::face_surface(f) else { continue };
            let uper = surf.is_u_periodic();
            let vper = surf.is_v_periodic();
            if !uper && !vper { continue; }
            let mut sff = occt_topo::shhealing::ShapeFixFace::with_face(f);
            let before_wires = wires_of_face(f).len();
            eprintln!("FIXMS-BEGIN face={i} uper={uper} vper={vper} wires={before_wires}");
            let r = sff.fix_missing_seam();
            let res = sff.result.clone();
            let desc = res.as_ref().map(|t| format!("{:?} faces={}", t.shape_type(), faces_of(t).len())).unwrap_or_else(|| "none".into());
            println!("FIXMS face={i} uper={uper} vper={vper} before_wires={before_wires} ret={r} result={desc}");
        }
        return;
    }
    if mode_fstats {
        let inc = IncrementalMesh::from_deflection(
            &shape,
            prs3d_get_deflection(&shape, 0.1),
            false,
            20.0_f64.to_radians(),
        );
        let stats = inc.face_stats();
        for (i, f) in faces_of(&shape).iter().enumerate() {
            let st = stats.get(i);
            let bb = occt_topo::brep_bnd_lib::shape_bnd_box(&f.0);
            let (b0, b1) = (bb.corner_min(), bb.corner_max());
            println!(
                "FSTAT face={i} mv={} mt={} bbox=({:.6},{:.6},{:.6})-({:.6},{:.6},{:.6})",
                st.map(|s| s.vertices).unwrap_or(0),
                st.map(|s| s.triangles).unwrap_or(0),
                b0.x(),
                b0.y(),
                b0.z(),
                b1.x(),
                b1.y(),
                b1.z()
            );
        }
        return;
    }
    if mode_low {
        let inc = IncrementalMesh::from_deflection(&shape, prs3d_get_deflection(&shape, 0.1), false, 20.0_f64.to_radians());
        let stats = inc.face_stats();
        for (i, f) in faces_of(&shape).iter().enumerate() {
            let st = stats.get(i);
            let mt = st.map(|s| s.triangles).unwrap_or(0);
            if mt >= 15 {
                continue;
            }
            let surf = BRepTool::face_surface(f);
            let tn = surf.as_ref().map(|s| tname(classify_surface(s.as_ref()))).unwrap_or("none");
            let (u0, u1, v0, v1) = surf.as_ref().map(|s| (s.u_range().0, s.u_range().1, s.v_range().0, s.v_range().1)).unwrap_or((0.0,0.0,0.0,0.0));
            let (uper, vper) = surf.as_ref().map(|s| (s.is_u_periodic(), s.is_v_periodic())).unwrap_or((false,false));
            let wires = wires_of_face(f);
            let we: Vec<usize> = wires.iter().map(|w| edges_of_wire(w).len()).collect();
            println!("LOW face={i} type={tn} mv={} mt={mt} wires={} edges={:?} u=[{u0:.4},{u1:.4}] v=[{v0:.4},{v1:.4}] uper={uper} vper={vper}",
                st.map(|s| s.vertices).unwrap_or(0), wires.len(), we);
        }
        return;
    }
    if mode_dev {
        use occt_topo::tgeometry::GeometryRegistry;
        let reg = GeometryRegistry::global();
        let mut printed = 0usize;
        for (fi, f) in faces_of(&shape).iter().enumerate() {
            let Some(surf) = BRepTool::face_surface(f) else { continue };
            let fk = GeometryRegistry::shape_key(&f.0);
            for w in wires_of_face(f) {
                for e in edges_of_wire(&w) {
                    let Some(c3d) = BRepTool::edge_curve(&e) else { continue };
                    let (first, last) = BRepTool::edge_parameters(&e);
                    let pcs = reg.edge_pcurves(&e.0, fk);
                    if pcs.is_empty() {
                        continue;
                    }
                    let (cf, cl) = reg.pcurve_range(&e.0, fk).unwrap_or((first, last));
                    let mut maxd: f64 = 0.0;
                    let mut samples: Vec<String> = Vec::new();
                    for k in 0..=8 {
                        let t = first + (last - first) * k as f64 / 8.0;
                        let pp = c3d.d0(t);
                        let uv = pcs[0].d0(t);
                        let q = surf.d0(uv.x(), uv.y());
                        let d = pp.distance(&q);
                        if d > maxd {
                            maxd = d;
                        }
                        if samples.len() < 9 {
                            samples.push(format!(
                                "t={t:.3} p=({:.3},{:.3},{:.3}) uv=({:.3},{:.3}) q=({:.3},{:.3},{:.3}) d={d:.3}",
                                pp.x(), pp.y(), pp.z(), uv.x(), uv.y(), q.x(), q.y(), q.z()
                            ));
                        }
                    }
                    if maxd > 5.0 && printed < 4 {
                        printed += 1;
                        let st = surf.as_ref();
                        let (u0, u1) = st.u_range();
                        let (v0, v1) = st.v_range();
                        println!(
                            "DEV face={fi} surf={} urange=[{u0:.4},{u1:.4}] vrange=[{v0:.4},{v1:.4}] uper={} vper={} range=[{first:.4},{last:.4}] prange=[{cf:.4},{cl:.4}] npc={} maxd={maxd:.4} c3d={} pc=[{:.4},{:.4}]",
                            tname(classify_surface(st)), st.is_u_periodic(), st.is_v_periodic(),
                            pcs.len(),
                            if c3d.is_line() { "Line" } else if c3d.gp_circ().is_some() { "Circle" } else if c3d.bspline_poles().is_some() { "BSpline" } else { "Other" },
                            pcs[0].first_parameter(), pcs[0].last_parameter()
                        );
                        for s in &samples {
                            println!("   {s}");
                        }
                    }
                }
            }
        }
        return;
    }
    if mode_edges {
        let mut n = 0usize; let mut bad = 0usize; let mut hist: std::collections::BTreeMap<String, (usize, usize)> = Default::default();
        for e in edges_of(&shape) {
            n += 1;
            let Some(curve) = BRepTool::edge_curve(&e) else { continue };
            let cn = if curve.is_line() { "Line" }
                else if curve.gp_circ().is_some() { "Circle" }
                else if curve.gp_ellipse().is_some() { "Ellipse" }
                else if curve.gp_hyperbola().is_some() { "Hyperbola" }
                else if curve.gp_parabola().is_some() { "Parabola" }
                else if curve.bspline_poles().is_some() { "BSpline" }
                else if curve.bezier_poles().is_some() { "Bezier" }
                else { "Other" };
            let (f, l) = BRepTool::edge_parameters(&e);
            let (v1, v2) = edge_vertices(&e);
            let (p1, p2) = match (&v1, &v2) {
                (Some(a), Some(b)) => (BRepTool::vertex_point(a), BRepTool::vertex_point(b)),
                _ => continue,
            };
            let d1 = curve.d0(f).distance(&p1);
            let d2 = curve.d0(l).distance(&p2);
            let d = d1.max(d2);
            let ent = hist.entry(cn.to_string()).or_insert((0, 0));
            ent.0 += 1;
            if d > 1e-3 { bad += 1; ent.1 += 1; 
                if d > 0.05 {
                    let (t1, t2) = (v1.as_ref().map(|v| BRepTool::vertex_tolerance(v)).unwrap_or(-1.0), v2.as_ref().map(|v| BRepTool::vertex_tolerance(v)).unwrap_or(-1.0));
                    let samples: String = (0..5).map(|k| {
                        let t = f + (l - f) * k as f64 / 4.0;
                        let q = curve.d0(t);
                        format!("({:.2},{:.2},{:.2})", q.x(), q.y(), q.z())
                    }).collect::<Vec<_>>().join(" ");
                    let info = format!(
                        "samples={samples} knots={:?} deg={:?} cr={:?} niv={} cont={} off={:?} trim={} basis={:?}",
                        curve.bspline_knots().map(|k| k.len()), curve.nurbs_degree(), curve.circle_radius(),
                        curve.nb_intervals(6), curve.continuity(), curve.offset_curve().map(|(_, o)| o),
                        curve.is_geom_trimmed(), curve.trimmed_basis_range());
                    println!("EDGE {n} type={cn} d1={d1:.6} d2={d2:.6} range=[{f:.6},{l:.6}] fp={:.6} lp={:.6} period={:.6} is_periodic={} closed={} {info} vtol=[{t1:.6},{t2:.6}] p1=({:.3},{:.3},{:.3}) p2=({:.3},{:.3},{:.3})",
                        curve.first_parameter(), curve.last_parameter(), curve.period(), curve.is_periodic(),
                        curve.d0(curve.first_parameter()).distance(&curve.d0(curve.last_parameter())) < 1e-9,
                        p1.x(), p1.y(), p1.z(), p2.x(), p2.y(), p2.z());
                }
            }
        }
        println!("EDGES total={n} residual>1e-3={bad}");
        for (k, (tot, b)) in &hist { println!("TYPE {k}: total={tot} bad={b}"); }
        return;
    }
    let b = brep_bnd_lib::shape_bnd_box(&shape);
    let (bmin, bmax) = (b.corner_min(), b.corner_max());
    println!("BBOX roots={} gap={:.6} min=({:.3},{:.3},{:.3}) max=({:.3},{:.3},{:.3})",
        model.shapes.len(), b.gap(), bmin.x(), bmin.y(), bmin.z(), bmax.x(), bmax.y(), bmax.z());
    let computed = prs3d_get_deflection(&shape, 0.1);
    let lin = override_lin.unwrap_or(computed);
    let angle = 20.0_f64.to_radians();
    let inc = IncrementalMesh::from_deflection(&shape, lin, false, angle);
    let stats = inc.face_stats();
    let (mv, mt) = inc.mesh().map(|m| (m.vertices.len(), m.triangles.len())).unwrap_or((0, 0));
    let faces = faces_of(&shape);
    println!("TOTAL faces={} computed_lin={:.6} used_lin={:.6} stats={} mesh_v={} mesh_t={}",
        faces.len(), computed, lin, stats.len(), mv, mt);
    let mut ftv: Vec<(f64, usize)> = faces.iter().enumerate().map(|(i, f)| (BRepTool::face_tolerance(f), i)).collect();
    ftv.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
    println!("TOP FACE TOLS: {}", ftv.iter().take(8).map(|(t, i)| format!("f{i}={t:.6}")).collect::<Vec<_>>().join(" "));
    let mut vtv: Vec<f64> = vertices_of(&shape).iter().map(|v| BRepTool::vertex_tolerance(v)).collect();
    vtv.sort_by(|a, b| b.partial_cmp(a).unwrap());
    println!("TOP VERT TOLS: {}", vtv.iter().take(8).map(|t| format!("{t:.6}")).collect::<Vec<_>>().join(" "));
    let mut matched = 0usize;
    let mut sum_mt = 0usize;
    let mut unmatched = 0usize;
    for f in faces.iter() {
        let key = GeometryRegistry::shape_key(&f.0);
        if let Some(s) = stats.iter().find(|s| s.shape_key == key) {
            matched += 1;
            sum_mt += s.triangles;
        } else {
            unmatched += 1;
        }
    }
    println!("STATMAP matched={matched} unmatched={unmatched} sum_mt={sum_mt} flat_mt={mt}");
    for (i, f) in faces.iter().enumerate() {
        let surf = BRepTool::face_surface(f);
        let tn = surf.as_ref().map(|s| tname(classify_surface(s.as_ref()))).unwrap_or("none");
        let wires = wires_of_face(f);
        let we: Vec<usize> = wires.iter().map(|w| edges_of_wire(w).len()).collect();
        // Pair the stat with the face by pointer identity: the mesh model skips
        // faces, so a positional lookup is off by the number of skipped faces.
        let key = GeometryRegistry::shape_key(&f.0);
        let st = stats.iter().find(|s| s.shape_key == key);
        let (sp_u, sp_v) = surf.as_ref().map(|s| { let (a, b) = s.u_range(); let (c, d) = s.v_range(); ((b - a).abs(), (d - c).abs()) }).unwrap_or((0.0, 0.0));
        let bb = brep_bnd_lib::shape_bnd_box(&f.0);
        let (bmin, bmax) = (bb.corner_min(), bb.corner_max());
        println!("F {} type={} tol={:.6} wires={} edges={:?} mv={} mt={} du={:.5} dv={:.5} stat_idx={:?} box=({:.5},{:.5},{:.5})-({:.5},{:.5},{:.5})",
            i, tn, BRepTool::face_tolerance(f), wires.len(), we,
            st.map(|s| s.vertices).unwrap_or(0), st.map(|s| s.triangles).unwrap_or(0),
            sp_u, sp_v, st.map(|s| s.index),
            bmin.x(), bmin.y(), bmin.z(), bmax.x(), bmax.y(), bmax.z());
    }
}