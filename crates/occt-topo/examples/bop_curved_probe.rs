//! Scratch probe: analytic curved boolean through the faithful kernel.
//!
//! Mirrors `specs/occt_probe/occt_probe.cpp --fuse`: box[-1,1]^3 fused with
//! `BRepPrimAPI_MakeCylinder(0.4, 2.0)`; OCCT ground truth is
//! `FUSE volume=8.50265 faces=8` with the box top face split by the circle
//! (area 4 -> 3.49735) plus the cylinder lateral (2.51327) and top cap
//! (0.502655). Prints the same shape facts from the port.

use occt_core::gp::GpPnt;
use occt_topo::abs::{Orientation, ShapeType};
use occt_topo::bop_builder::{boolean, BoolOp};
use occt_topo::brep_gprop::volume as shape_volume;
use occt_topo::brep_surface::classify_surface;
use occt_topo::brep_tool::BRepTool;
use occt_topo::primitives::{BRepPrimBox, BRepPrimCylinder};
use occt_topo::shell_check::shell_invariants;
use occt_topo::topo_tools::{edges_of, faces_of, shapes_of};
use occt_topo::wireframe::wireframe_face_area;

fn report(tag: &str, a: &occt_topo::shape::TopoShape, b: &occt_topo::shape::TopoShape, op: BoolOp) {
    match boolean(a, b, op, 1e-6) {
        Ok(r) => {
            let nf = faces_of(&r.shape).len();
            let ns = shapes_of(&r.shape, ShapeType::Solid).len();
            let closed: Vec<bool> = r.shells.iter().map(|s| shell_invariants(s).is_valid_solid()).collect();
            let v = shape_volume(&r.shape, 0.02);
            let vg = occt_topo::brep_gprop_full::volume_properties(&r.shape)
                .map(|g| g.mass())
                .unwrap_or(f64::NAN);
            println!(
                "{tag}: type={:?} faces={nf} solids={ns} closed={closed:?} vol={v:.5} vol_gauss={vg:.5} warnings={:?}",
                r.shape.shape_type(),
                r.warnings
            );
            for (i, f) in faces_of(&r.shape).iter().enumerate() {
                let kind = BRepTool::face_surface(f).map(|s| classify_surface(s.as_ref()));
                let (u1, u2, v1, v2) = BRepTool::uv_bounds(f);
                let place = match BRepTool::face_surface(f).as_ref().and_then(|s| s.gp_pln()) {
                    Some(p) => {
                        let loc = p.pos.location();
                        let d = p.pos.direction();
                        format!("origin=({:.4},{:.4},{:.4}) n=({:.4},{:.4},{:.4})", loc.x(), loc.y(), loc.z(), d.x(), d.y(), d.z())
                    }
                    None => match BRepTool::face_surface(f).as_ref().and_then(|s| s.gp_cylinder()) {
                        Some(cy) => {
                            let loc = cy.pos.location();
                            format!("origin=({:.4},{:.4},{:.4}) r={:.4}", loc.x(), loc.y(), loc.z(), cy.radius())
                        }
                        None => "?".to_string(),
                    },
                };
                let wires: Vec<String> = occt_topo::topo_tools::wires_of_face(f)
                    .iter()
                    .map(|w| {
                        let es: Vec<String> = occt_topo::topo_tools::edges_of_wire(w)
                            .iter()
                            .map(|e| {
                                let r = BRepTool::edge_parameters(e);
                                format!("[{:.4},{:.4}]", r.0, r.1)
                            })
                            .collect();
                        format!("{}e {}", es.len(), es.join(" "))
                    })
                    .collect();
                println!(
                    "  RFACE {i} kind={kind:?} area={:.5} uv=({u1:.4},{u2:.4})x({v1:.4},{v2:.4}) {place} wires=[{}]",
                    wireframe_face_area(f),
                    wires.join(" | ")
                );
            }
        }
        Err(e) => println!("{tag}: ERR {e}"),
    }
}

/// Source-shape index map plus the face surface placement / UV bounds, so the
/// interference records below can be read against the actual geometry.
fn dump_shapes(a: &occt_topo::shape::TopoShape, b: &occt_topo::shape::TopoShape) {
    let mut idx = 0usize;
    for (tag, s) in [("A", a), ("B", b)] {
        for f in faces_of(s) {
            let surf = BRepTool::face_surface(&f);
            let kind = surf.as_ref().map(|s| classify_surface(s.as_ref()));
            let (u1, u2, v1, v2) = BRepTool::uv_bounds(&f);
            let place = match surf.as_ref().and_then(|s| s.gp_pln().map(|p| (p, 0))) {
                Some((p, _)) => {
                    let loc = p.pos.location();
                    let d = p.pos.direction();
                    format!("origin=({:.4},{:.4},{:.4}) n=({:.4},{:.4},{:.4})", loc.x(), loc.y(), loc.z(), d.x(), d.y(), d.z())
                }
                None => match surf.as_ref().and_then(|s| s.gp_cylinder()) {
                    Some(cy) => {
                        let loc = cy.pos.location();
                        let d = cy.pos.direction();
                        format!("origin=({:.4},{:.4},{:.4}) ax=({:.4},{:.4},{:.4}) r={:.4}", loc.x(), loc.y(), loc.z(), d.x(), d.y(), d.z(), cy.radius())
                    }
                    None => "?".to_string(),
                },
            };
            println!("SRC {tag}{idx}: Face kind={kind:?} uv=({u1:.4},{u2:.4})x({v1:.4},{v2:.4}) {place}");
            idx += 1;
        }
    }
}

fn dump_interf(ds: &occt_topo::bopds::BopdsDS) {
    println!(
        "INTERF: vv={} ve={} vf={} ee={} ef={} ff={}",
        ds.interf_vv().len(),
        ds.interf_ve().len(),
        ds.interf_vf().len(),
        ds.interf_ee().len(),
        ds.interf_ef().len(),
        ds.interf_ff().len()
    );
    for (tag, list) in [
        ("VV", ds.interf_vv()),
        ("VE", ds.interf_ve()),
        ("VF", ds.interf_vf()),
        ("EE", ds.interf_ee()),
        ("EF", ds.interf_ef()),
    ] {
        for r in list {
            let (i1, i2) = r.indices();
            let t1 = ds.shape(i1).map(|s| format!("{:?}", s.shape_type())).unwrap_or_default();
            let t2 = ds.shape(i2).map(|s| format!("{:?}", s.shape_type())).unwrap_or_default();
            println!("  {tag} ({i1}:{t1},{i2}:{t2}) new={} t=({:.5},{:.5})", r.index_new, r.common_first, r.common_last);
        }
    }
    for r in ds.interf_ff() {
        let (i1, i2) = r.indices();
        println!(
            "  FF ({i1},{i2}) tangent={} curves={} points={}",
            r.tangent_faces(),
            r.curves().len(),
            r.points().len(),
        );
    }
}

fn dump_edge_pcurves(ds: &occt_topo::bopds::BopdsDS, faces: &[(String, occt_topo::shape::TopoShape)]) {
    use occt_topo::tgeometry::GeometryRegistry;
    let reg = GeometryRegistry::global();
    for i in ds.nb_source_shapes()..ds.nb_shapes() {
        let Some(s) = ds.shape(i) else { continue };
        if s.shape_type() != occt_topo::abs::ShapeType::Edge {
            continue;
        }
        let e = match occt_topo::shape::Edge::wrap(s.clone()) {
            Some(e) => e,
            None => continue,
        };
        let r = BRepTool::edge_parameters(&e);
        print!("PCURVE edge {i} range=({:.5},{:.5})", r.0, r.1);
        if let Some(c3) = BRepTool::edge_curve(&e) {
            let a = c3.d0(r.0);
            let b = c3.d0(r.1);
            print!(
                " p3=( {:.4},{:.4},{:.4})->({:.4},{:.4},{:.4})",
                a.x(), a.y(), a.z(), b.x(), b.y(), b.z()
            );
        }
        for (tag, f) in faces {
            let fk = GeometryRegistry::shape_key(f);
            let rf = occt_topo::shape::Face::wrap(f.clone());
            let n = match &rf {
                Some(f2) => reg.edge_pcurves(&e.0, fk).len().max(usize::from(
                    occt_topo::boptools_2d::curve_on_surface(&e, f2).is_some(),
                )),
                None => 0,
            };
            if n > 0 {
                let mut info = String::new();
                if let Some(f2) = &rf {
                    if let Some((c2d, t1, t2)) = occt_topo::boptools_2d::curve_on_surface_range(&e, f2) {
                        let a = c2d.d0(t1);
                        let b = c2d.d0(t2);
                        let m = c2d.d0(0.5 * (t1 + t2));
                        info = format!(
                            " 2d=({:.4},{:.4})->({:.4},{:.4}) mid=({:.4},{:.4})",
                            a.x(), a.y(), b.x(), b.y(), m.x(), m.y()
                        );
                    }
                }
                print!(" | {tag}:n={n}{info}");
            }
        }
        println!();
    }
}

// `BOPAlgo_Builder::BuildSplitFaces` step 1 for one face, then the
// `BOPAlgo_WireSplitter::SplitBlock` walk: boundary images + section edges.
fn dump_split_block(
    ds: &occt_topo::bopds::BopdsDS,
    face_shape: &occt_topo::shape::TopoShape,
    label: &str,
    extra: &[usize],
) {
    use occt_topo::shape::Edge;
    let face = match occt_topo::shape::Face::wrap(face_shape.clone()) {
        Some(f) => f,
        None => return,
    };
    let idx_of = |s: &occt_topo::shape::TopoShape| -> Option<usize> {
        (0..ds.nb_shapes()).find(|&i| ds.shape(i).map(|x| x.same_tshape(s)).unwrap_or(false))
    };
    let mut le: Vec<Edge> = Vec::new();
    let face_edges: Vec<Edge> = {
        let mut v = Vec::new();
        for w in occt_topo::topo_tools::wires_of_face(&face) {
            v.extend(occt_topo::topo_tools::edges_of_wire(&w));
        }
        v
    };
    for be in face_edges {
        let ori = be.0.orientation();
        let mut replaced = false;
        let i = idx_of(&be.0);
        if let Some(i) = i {
            let pbs = ds.pave_blocks(i);
            if !pbs.is_empty() {
                for pb in pbs {
                    if let Some(e) = Edge::wrap(ds.shape(pb.edge()).cloned().unwrap_or_else(|| {
                        occt_topo::shape::TopoShape::new(occt_topo::abs::ShapeType::Edge)
                    })) {
                        let mut e = e;
                        e.0.set_orientation(ori);
                        le.push(e);
                    }
                }
                replaced = true;
            }
        }
        println!(
            "  BND {label} edge ds={:?} ori={ori:?} pbcount={:?} replaced={replaced} has_pc={}",
            i,
            i.map(|k| ds.pave_blocks(k).len()),
            occt_topo::boptools_2d::curve_on_surface(
                &be,
                &occt_topo::shape::Face::wrap(face_shape.clone()).unwrap()
            )
            .is_some(),
        );
        if !replaced {
            le.push(be);
        }
    }
    let nb_iface = le.len();
    for &i in extra {
        let Some(s) = ds.shape(i) else { continue };
        let Some(e) = Edge::wrap(s.clone()) else { continue };
        let mut f = e.clone();
        f.0.set_orientation(Orientation::Forward);
        le.push(f);
        let mut r = e.clone();
        r.0.set_orientation(Orientation::Reversed);
        le.push(r);
    }
    let n_curved = le.len();
    println!("SPLITBLOCK {label}: input={n_curved} (boundary={nb_iface} + section={})", n_curved - nb_iface);
    let mut le2 = le.clone();
    match occt_topo::wire_splitter::WireSplitter::split_block(&face, le.clone()) {
        Ok(wires) => {
            println!("  split_block -> {} wires", wires.len());
            for (i, w) in wires.iter().enumerate() {
                let es: Vec<String> = edges_of(w)
                    .iter()
                    .map(|e| match idx_of(&e.0) {
                        Some(k) => format!("{k}{}", if e.0.orientation() == Orientation::Reversed { "-" } else { "+" }),
                        None => "?".to_string(),
                    })
                    .collect();
                println!("    W{i} ({}) edges=[{}]", es.len(), es.join(" "));
            }
        }
        Err(e) => println!("  split_block ERR {e}"),
    }
    // Coordinate-only chaining, for contrast.
    match occt_topo::wire_splitter::WireSplitter::make_wire(&mut le2) {
        Ok(w) => println!("  greedy first loop: {} edges", edges_of(&w).len()),
        Err(e) => println!("  greedy ERR {e}"),
    }
}

fn dump_classifier(f: &occt_topo::shape::Face) {
    use occt_topo::abs::ShapeType as ST;
    let kind = BRepTool::face_surface(f).map(|s| classify_surface(s.as_ref()));
    println!("CLS face kind={kind:?}");
    match occt_topo::fclass2d::FClass2d::new(f, 1e-9) {
        Ok(cl) => {
            println!(
                "  outer={} holes={}",
                cl.outer_ring().map(|r| r.len()).unwrap_or(0),
                cl.hole_rings().len()
            );
            if let Some(r) = cl.outer_ring() {
                println!("    ring len={}", r.len());
                for p in r.iter() {
                    println!("    ring ({:.6},{:.6})", p.x(), p.y());
                }
            }
            for (i, h) in cl.hole_rings().iter().enumerate() {
                println!("  hole {i}: {} pts, first ({:.6},{:.6})", h.len(), h[0].x(), h[0].y());
            }
            for p in [(4.712389, 1.0), (0.0, 1.0), (1.5708, 1.0), (3.14159, 1.0)] {
                println!(
                    "  state ({:.5},{:.5}) = {:?}",
                    p.0,
                    p.1,
                    cl.perform(occt_core::gp::GpPnt2d::new(p.0, p.1))
                );
            }
        }
        Err(e) => println!("  CLS ERR {e}"),
    }
    let _ = ST::Face;
}

fn dump_pave_filler(a: &occt_topo::shape::TopoShape, b: &occt_topo::shape::TopoShape) {
    let mut pf = occt_topo::pave_filler::PaveFiller::new();
    pf.set_fuzzy_value(1e-7);
    pf.set_arguments(&[a.clone(), b.clone()]);
    match pf.perform() {
        Ok(()) => {}
        Err(e) => {
            println!("PAVEFILLER ERR: {e}");
            return;
        }
    }
    let ds = pf.ds();
    println!(
        "PAVEFILLER: nb_shapes={} nb_source={} nb_ranges={} n_interf_pairs={}",
        ds.nb_shapes(),
        ds.nb_source_shapes(),
        ds.nb_ranges(),
        ds.interferences().len()
    );
    let n_src = ds.nb_source_shapes();
    for i in 0..n_src {
        let Some(si) = ds.shape_info(i) else { continue };
        if si.shape_type() != occt_topo::abs::ShapeType::Face {
            continue;
        }
        if let Some(fi) = ds.face_info(i) {
            println!(
                "  FACE {i}: paves={} paves_in={} paves_on={} verts={} verts_sc={}",
                fi.paves().len(),
                fi.paves_in().len(),
                fi.paves_on().len(),
                fi.verts().len(),
                fi.verts_sc().len()
            );
        }
    }
    for i in n_src..ds.nb_shapes() {
        if let Some(s) = ds.shape(i) {
            let extra = if s.shape_type() == occt_topo::abs::ShapeType::Vertex {
                match occt_topo::shape::Vertex::wrap(s.clone()) {
                    Some(v) => {
                        let p = BRepTool::vertex_point(&v);
                        format!(" p=({:.5},{:.5},{:.5})", p.x(), p.y(), p.z())
                    }
                    None => String::new(),
                }
            } else if s.shape_type() == occt_topo::abs::ShapeType::Edge {
                match occt_topo::shape::Edge::wrap(s.clone()) {
                    Some(e) => {
                        let r = BRepTool::edge_parameters(&e);
                        format!(" range=({:.5},{:.5})", r.0, r.1)
                    }
                    None => String::new(),
                }
            } else {
                String::new()
            };
            println!("  NEW shape {i}: {:?}{extra}", s.shape_type());
        }
    }
    let ff = ds.interf_ff();
    println!("PAVEFILLER: interf_ff={} warnings={:?}", ff.len(), pf.warnings());
    dump_interf(ds);
    {
        let mut faces: Vec<(String, occt_topo::shape::TopoShape)> = Vec::new();
        for i in 0..ds.nb_source_shapes() {
            let Some(s) = ds.shape(i) else { continue };
            if s.shape_type() == occt_topo::abs::ShapeType::Face {
                faces.push((format!("F{i}"), s.clone()));
            }
        }
        dump_edge_pcurves(ds, &faces);
        for (tag, s) in &faces {
            if tag == "F11" || tag == "F36" {
                dump_split_block(ds, s, tag, &[42, 43]);
            }
        }
    }
    for rec in ff {
        let (i1, i2) = rec.indices();
        println!(
            "  FF ({i1},{i2}) tangent={} curves={} points={}",
            rec.tangent_faces(),
            rec.curves().len(),
            rec.points().len()
        );
        for c in rec.curves() {
            let (f, l) = c.range();
            let bx = c.bounding_box();
            let (mn, mx) = (bx.corner_min(), bx.corner_max());
            println!(
                "    curve kind range=({f:.5},{l:.5}) tol={:.3e} has_edge={} box=({:.3},{:.3},{:.3})-({:.3},{:.3},{:.3}) nblocks={}",
                c.tolerance(),
                c.has_edge(),
                mn.x(), mn.y(), mn.z(),
                mx.x(), mx.y(), mx.z(),
                c.pave_blocks().len(),
            );
            {
                let c3 = c.curve().expect("curve");
                let p0 = c3.d0(f);
                let pm = c3.d0(0.5 * (f + l));
                let p1 = c3.d0(l);
                println!(
                    "      c3d p0=({:.4},{:.4},{:.4}) pm=({:.4},{:.4},{:.4}) p1=({:.4},{:.4},{:.4})",
                    p0.x(), p0.y(), p0.z(), pm.x(), pm.y(), pm.z(), p1.x(), p1.y(), p1.z()
                );
            }
            for pb in c.pave_blocks() {
                let (a, b) = pb.indices();
                println!(
                    "      PB v=({a},{b}) t=({:.5},{:.5}) edge={} orig={}",
                    pb.first, pb.last, pb.edge(), pb.original_edge()
                );
            }
        }
    }
    for i in 0..ds.nb_shapes() {
        let Some(s) = ds.shape(i) else { continue };
        if s.shape_type() != occt_topo::abs::ShapeType::Edge {
            continue;
        }
        let pbs = ds.pave_blocks(i);
        if pbs.is_empty() {
            continue;
        }
        println!("  EDGE {i} pbs={}", pbs.len());
        for pb in pbs {
            let (a, b) = pb.indices();
            println!(
                "      PB v=({a},{b}) t=({:.5},{:.5}) edge={} orig={}",
                pb.first, pb.last, pb.edge(), pb.original_edge()
            );
        }
    }
}

fn main() {
    // Exact `--fuse` setup: box[-1,1]^3 and cylinder(r=0.4,h=2) at the origin.
    let box_s = BRepPrimBox::make_box_corner(&GpPnt::new(-1.0, -1.0, -1.0), &GpPnt::new(1.0, 1.0, 1.0)).solid;
    let cyl = BRepPrimCylinder::make_cylinder(0.4, 2.0).solid;
    let gauss = |s: &occt_topo::shape::TopoShape| {
        occt_topo::brep_gprop_full::volume_properties(s).map(|g| g.mass()).unwrap_or(f64::NAN)
    };
    println!(
        "box vol={:.5} vol_gauss={:.5} faces={}",
        shape_volume(&box_s.0, 0.02),
        gauss(&box_s.0),
        faces_of(&box_s.0).len()
    );
    println!(
        "cyl vol={:.5} vol_gauss={:.5} faces={}",
        shape_volume(&cyl.0, 0.02),
        gauss(&cyl.0),
        faces_of(&cyl.0).len()
    );
    for (r, h) in [(0.4, 1.0), (0.4, 2.0), (0.4, 4.0), (1.0, 1.0)] {
        let c = BRepPrimCylinder::make_cylinder(r, h).solid;
        println!(
            "cyl r={r} h={h} vol={:.5} vol_gauss={:.5} analytic={:.5}",
            shape_volume(&c.0, 0.02),
            gauss(&c.0),
            std::f64::consts::PI * r * r * h
        );
    }
    println!("OCCT GT: FUSE volume=8.50265 faces=8");

    report("box|cyl FUSE", &box_s.0, &cyl.0, BoolOp::Fuse);
    report("box|cyl CUT", &box_s.0, &cyl.0, BoolOp::Cut);
    report("box|cyl COMMON", &box_s.0, &cyl.0, BoolOp::Common);

    dump_shapes(&box_s.0, &cyl.0);
    // After the fuse has run, the pave-filler state is rebuilt on a fresh run so
    // the DS still holds the source shapes for the pcurve dump.
    dump_pave_filler(&box_s.0, &cyl.0);

    println!("--- classifiers ---");
    for (i, f) in faces_of(&cyl.0).iter().enumerate() {
        println!("cyl face {i}:");
        dump_classifier(f);
    }

    println!("--- FF box faces vs cylinder faces ---");
    let bf = faces_of(&box_s.0);
    let cf = faces_of(&cyl.0);
    for (i, f1) in bf.iter().enumerate() {
        let kind = BRepTool::face_surface(f1).map(|s| classify_surface(s.as_ref()));
        for (j, f2) in cf.iter().enumerate() {
            let mut ff = occt_topo::int_face_face::FaceFace::new();
            ff.set_face1(f1.clone());
            ff.set_face2(f2.clone());
            ff.set_tolerance(1e-6);
            if ff.perform().is_ok() {
                let r = ff.result();
                if r.nb_curves() > 0 {
                    let cs = r.curves();
                    let kinds: Vec<_> = cs.iter().map(|c| format!("{:?}", c.kind)).collect();
                    let tols: Vec<_> = cs.iter().map(|c| format!("{:.3e}", c.tolerance)).collect();
                    println!("box face {i} ({kind:?}) x cyl face {j}: n={} kinds={kinds:?} tols={tols:?}", cs.len());
                }
            }
        }
    }

    println!("--- direct analytic plane x cylinder ---");
    for f1 in &bf {
        let Some(pl) = BRepTool::face_surface(f1).and_then(|s| s.gp_pln()) else { continue };
        if (pl.pos.location().z() - 1.0).abs() > 1e-9 {
            continue;
        }
        for f2 in &cf {
            let Some(cy) = BRepTool::face_surface(f2).and_then(|s| s.gp_cylinder()) else { continue };
            if (cy.radius() - 0.4).abs() > 1e-9 {
                continue;
            }
            let cp = cy.pos;
            println!(
                "  cyl loc=({:.4},{:.4},{:.4}) zd=({:.4},{:.4},{:.4}) xd=({:.4},{:.4},{:.4}) yd=({:.4},{:.4},{:.4})",
                cp.location().x(), cp.location().y(), cp.location().z(),
                cp.direction().x(), cp.direction().y(), cp.direction().z(),
                cp.x_direction().x(), cp.x_direction().y(), cp.x_direction().z(),
                cp.y_direction().x(), cp.y_direction().y(), cp.y_direction().z()
            );
            let qi = occt_geom::intana::quadric_quadric_plane_cylinder(&pl, &cy, 1e-12, 1e-7, 0.0);
            match qi {
                occt_geom::intana::QuadricIntersection::Circle(c) => {
                    let p = c.position();
                    println!(
                        "  analytic Circle loc=({:.4},{:.4},{:.4}) xd=({:.4},{:.4},{:.4}) yd=({:.4},{:.4},{:.4}) r={:.4}",
                        p.location().x(), p.location().y(), p.location().z(),
                        p.x_direction().x(), p.x_direction().y(), p.x_direction().z(),
                        p.y_direction().x(), p.y_direction().y(), p.y_direction().z(),
                        c.radius()
                    );
                }
                other => println!("  analytic other={other:?}"),
            }
        }
    }
}
