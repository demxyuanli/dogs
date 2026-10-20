//! Pcurve / wire census diagnostics (--pcensus, --ecensus, --checkwire).
use super::tname;
use super::Ctx;
use occt_topo::brep_bnd_lib;
use occt_topo::brep_tool::BRepTool;
use occt_topo::meshing::range_splitter::mesh_surface_type;
use occt_topo::tgeometry::GeometryRegistry;
use occt_topo::topo_tools::{edges_of_wire, faces_of, wires_of_face};

pub(crate) fn run(ctx: &Ctx) -> bool {
    let shape = ctx.shape.clone();
    if std::env::args().any(|a| a == "--pcensus") {
        // TEMP determinism probe: per (face, wire, edge, pcurve) the stored
        // pcurve kind, stored window and endpoint samples, so two runs can be
        // diffed to the exact edge whose pcurve differs.
        let reg = GeometryRegistry::global();
        let surfs: Vec<String> = faces_of(&shape)
            .iter()
            .enumerate()
            .map(|(i, f)| {
                let fk = GeometryRegistry::shape_key(&f.0);
                format!(
                    "SURF f={i} {}",
                    reg.debug_surface_desc(reg.debug_repr_key(fk))
                )
            })
            .collect();
        for s in surfs {
            println!("{s}");
        }
        for (i, f) in faces_of(&shape).iter().enumerate() {
            let fk = GeometryRegistry::shape_key(&f.0);
            for (wi, w) in wires_of_face(f).iter().enumerate() {
                for (ei, e) in edges_of_wire(w).iter().enumerate() {
                    let (f3, l3) = occt_topo::brep_tool::BRepTool::edge_parameters(e);
                    let c3 = occt_topo::brep_tool::BRepTool::edge_curve(e);
                    let c3k = match c3.as_ref() {
                        None => "none",
                        Some(c) if c.is_line() => "line",
                        Some(c) if c.gp_circ().is_some() => "circle",
                        Some(c) if c.nurbs_degree().is_some() => "bspline",
                        Some(_) => "other",
                    };
                    let pcs = reg.edge_pcurves(&e.0, fk);
                    let (rf, rl) = reg.pcurve_range(&e.0, fk).unwrap_or((f64::NAN, f64::NAN));
                    let mut desc = String::new();
                    for (k, pc) in pcs.iter().enumerate() {
                        let kind = if pc.is_geom2d_line() {
                            "line"
                        } else if pc.is_bspline2d() {
                            "bspline"
                        } else {
                            "other"
                        };
                        let (a, b) = (pc.first_parameter(), pc.last_parameter());
                        let (p1, p2) = if a.is_finite() && b.is_finite() {
                            let q1 = pc.d0(a);
                            let q2 = pc.d0(b);
                            (
                                format!("({:.15e},{:.15e})", q1.x(), q1.y()),
                                format!("({:.15e},{:.15e})", q2.x(), q2.y()),
                            )
                        } else {
                            ("inf".to_string(), "inf".to_string())
                        };
                        desc.push_str(&format!(
                            " pc{k}={kind}[{a:.15e},{b:.15e}]{p1}->{p2}"
                        ));
                    }
                    println!(
                        "PC f={i} w={wi} e={ei} c3={c3k} r3=[{f3:.15e},{l3:.15e}] rw=[{rf:.15e},{rl:.15e}] n={}{desc}",
                        pcs.len()
                    );
                }
            }
        }
        return true;
    }
    if std::env::args().any(|a| a == "--ecensus") {
        // T-102 read-only census: one line per face with `wires=` and the
        // per-wire edge counts, plus the geometric bbox and the face's
        // vertex-point signature, so the port's face can be paired with the
        // STEP entity table by GEOMETRY alone (never by face index).
        // No library call is added -- same traversal the default dump uses,
        // just without building the mesh.
        for (i, f) in faces_of(&shape).iter().enumerate() {
            let bb = brep_bnd_lib::shape_bnd_box(&f.0);
            let (mn, mx) = (bb.corner_min(), bb.corner_max());
            let ws = wires_of_face(f);
            let we: Vec<usize> = ws.iter().map(|w| edges_of_wire(w).len()).collect();
            let mut vs: Vec<(i64, i64, i64)> = Vec::new();
            for w in ws.iter() {
                for e in edges_of_wire(w) {
                    for v in [occt_topo::shhealing::first_vertex(&e),
                              occt_topo::shhealing::last_vertex(&e)]
                        .into_iter()
                        .flatten()
                    {
                        let p = BRepTool::vertex_point(&v);
                        vs.push((
                            (p.x() * 1e6).round() as i64,
                            (p.y() * 1e6).round() as i64,
                            (p.z() * 1e6).round() as i64,
                        ));
                    }
                }
            }
            vs.sort_unstable();
            vs.dedup();
            let (vmn, vmx) = if vs.is_empty() {
                ((0i64, 0i64, 0i64), (0i64, 0i64, 0i64))
            } else {
                (
                    (
                        vs.iter().map(|v| v.0).min().unwrap(),
                        vs.iter().map(|v| v.1).min().unwrap(),
                        vs.iter().map(|v| v.2).min().unwrap(),
                    ),
                    (
                        vs.iter().map(|v| v.0).max().unwrap(),
                        vs.iter().map(|v| v.1).max().unwrap(),
                        vs.iter().map(|v| v.2).max().unwrap(),
                    ),
                )
            };
            let vlist: Vec<String> = vs
                .iter()
                .map(|v| format!("{:.6},{:.6},{:.6}", v.0 as f64 / 1e6, v.1 as f64 / 1e6, v.2 as f64 / 1e6))
                .collect();
            println!(
                "CENSUS face={i} wires={} edges={:?} bbox=({:.6},{:.6},{:.6})-({:.6},{:.6},{:.6}) vbox=({:.6},{:.6},{:.6})-({:.6},{:.6},{:.6}) nv={} verts=[{}]",
                ws.len(),
                we,
                mn.x(), mn.y(), mn.z(), mx.x(), mx.y(), mx.z(),
                vmn.0 as f64 / 1e6, vmn.1 as f64 / 1e6, vmn.2 as f64 / 1e6,
                vmx.0 as f64 / 1e6, vmx.1 as f64 / 1e6, vmx.2 as f64 / 1e6,
                vs.len(),
                vlist.join(";")
            );
        }
        return true;
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
            if let Some(surf) = BRepTool::face_surface(f) {
                let (u0, u1) = surf.u_range();
                let (v0, v1) = surf.v_range();
                println!(
                    "FDUMP  surf={} u=[{u0:.6},{u1:.6}] v=[{v0:.6},{v1:.6}] uper={} vper={}",
                    tname(mesh_surface_type(surf.as_ref())),
                    surf.is_u_periodic(),
                    surf.is_v_periodic()
                );
            }
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
            return true;
        }
        println!("FDUMP none");
        return true;
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
        return true;
    }
    false
}
