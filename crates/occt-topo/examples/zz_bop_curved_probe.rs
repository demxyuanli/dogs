//! Scratch probe: analytic curved boolean through the faithful kernel.
//!
//! Mirrors `specs/occt_probe/occt_probe.cpp --fuse`: box[-1,1]^3 fused with
//! `BRepPrimAPI_MakeCylinder(0.4, 2.0)`; OCCT ground truth is
//! `FUSE volume=8.50265 faces=8` with the box top face split by the circle
//! (area 4 -> 3.49735) plus the cylinder lateral (2.51327) and top cap
//! (0.502655). Prints the same shape facts from the port.

use occt_core::gp::GpPnt;
use occt_topo::abs::ShapeType;
use occt_topo::bop_builder::{boolean, BoolOp};
use occt_topo::brep_gprop::volume as shape_volume;
use occt_topo::brep_surface::classify_surface;
use occt_topo::brep_tool::BRepTool;
use occt_topo::primitives::{BRepPrimBox, BRepPrimCylinder};
use occt_topo::shell_check::shell_invariants;
use occt_topo::topo_tools_full::{faces_of, shapes_of};
use occt_topo::wireframe::face_area;

fn report(tag: &str, a: &occt_topo::shape::TopoShape, b: &occt_topo::shape::TopoShape, op: BoolOp) {
    match boolean(a, b, op, 1e-6) {
        Ok(r) => {
            let nf = faces_of(&r.shape).len();
            let ns = shapes_of(&r.shape, ShapeType::Solid).len();
            let closed: Vec<bool> = r.shells.iter().map(|s| shell_invariants(s).is_valid_solid()).collect();
            let v = shape_volume(&r.shape, 0.02);
            println!(
                "{tag}: type={:?} faces={nf} solids={ns} closed={closed:?} vol={v:.5} warnings={:?}",
                r.shape.shape_type(),
                r.warnings
            );
            for (i, f) in faces_of(&r.shape).iter().enumerate() {
                let kind = BRepTool::face_surface(f).map(|s| classify_surface(s.as_ref()));
                println!("  RFACE {i} kind={kind:?} area={:.5}", face_area(f));
            }
        }
        Err(e) => println!("{tag}: ERR {e}"),
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
            println!("  NEW shape {i}: {:?}", s.shape_type());
        }
    }
    let ff = ds.interf_ff();
    println!("PAVEFILLER: interf_ff={} warnings={:?}", ff.len(), pf.warnings());
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
    println!("box vol={:.5} faces={}", shape_volume(&box_s.0, 0.02), faces_of(&box_s.0).len());
    println!("cyl vol={:.5} faces={}", shape_volume(&cyl.0, 0.02), faces_of(&cyl.0).len());
    println!("OCCT GT: FUSE volume=8.50265 faces=8");

    report("box|cyl FUSE", &box_s.0, &cyl.0, BoolOp::Fuse);
    report("box|cyl CUT", &box_s.0, &cyl.0, BoolOp::Cut);
    report("box|cyl COMMON", &box_s.0, &cyl.0, BoolOp::Common);

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
}
