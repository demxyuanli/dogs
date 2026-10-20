//! Default report: per-face summary or per-mode counters.
use super::tname;
use super::Ctx;
use occt_topo::brep_bnd_lib;
use occt_topo::brep_exchange::prs3d_get_deflection;
use occt_topo::brep_tool::BRepTool;
use occt_topo::meshing::incremental_mesh::IncrementalMesh;
use occt_topo::meshing::range_splitter::mesh_surface_type;
use occt_topo::tgeometry::GeometryRegistry;
use occt_topo::topo_tools::{edge_vertices, edges_of, edges_of_wire, faces_of, vertices_of, wires_of_face};

pub(crate) fn run(ctx: &Ctx) {
    let model = &ctx.model;
    let shape = ctx.shape.clone();
    let mode_low = ctx.mode_low;
    let mode_dev = ctx.mode_dev;
    let mode_edges = ctx.mode_edges;
    let override_lin = ctx.override_lin;
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
            let tn = surf.as_ref().map(|s| tname(mesh_surface_type(s.as_ref()))).unwrap_or("none");
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
                            tname(mesh_surface_type(st)), st.is_u_periodic(), st.is_v_periodic(),
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
        let tn = surf.as_ref().map(|s| tname(mesh_surface_type(s.as_ref()))).unwrap_or("none");
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
