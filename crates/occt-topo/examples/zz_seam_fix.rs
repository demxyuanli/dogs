//! TEMPORARY diagnostic probe (delete after T-93).
//!
//! Applies the port's `ShapeFixFace::fix_missing_seam` to one face of a STEP
//! file and reports:
//!   * the resulting wire/edge structure (to diff against OCCT's healed import),
//!   * whether the healed face produces a valid `RangeSplitter` range
//!     (`IsValid`, i.e. the gate that was rejecting it).
//!
//! Usage: `cargo run --example zz_seam_fix -- <file.step> <face-index>`
use occt_topo::brep_bnd_lib;
use occt_topo::brep_tool::BRepTool;
use occt_topo::builder::TopoBuilder;
use occt_topo::meshing::data_model::{MeshFace, MeshModel};
use occt_topo::meshing::edge_discret::EdgeDiscret;
use occt_topo::meshing::model_builder::ModelBuilder;
use occt_topo::meshing::parameters::MeshParameters;
use occt_topo::meshing::range_splitter::create_range_splitter;
use occt_topo::shhealing::ShapeFixFace;
use occt_topo::shape::Face;
use occt_topo::step::read_step_file;
use occt_topo::topo_tools_full::{edges_of_wire, faces_of, wires_of_face};

/// First / last vertex of a wire, in wire order (orientation-aware).
fn first_vertex_of(es: &[occt_topo::shape::Edge]) -> Option<occt_topo::shape::Vertex> {
    occt_topo::shhealing::first_vertex(es.first()?)
}

fn last_vertex_of(es: &[occt_topo::shape::Edge]) -> Option<occt_topo::shape::Vertex> {
    occt_topo::shhealing::last_vertex(es.last()?)
}

fn describe(tag: &str, face: &Face) {    let bb = brep_bnd_lib::shape_bnd_box(&face.0);
    let (b0, b1) = (bb.corner_min(), bb.corner_max());
    let surf = BRepTool::face_surface(face);
    let (ur, vr) = surf
        .as_ref()
        .map(|s| (s.u_range(), s.v_range()))
        .unwrap_or(((0.0, 0.0), (0.0, 0.0)));
    println!(
        "{tag} wires={} urange=[{:.9},{:.9}] vrange=[{:.9},{:.9}] bbox=({:.6},{:.6},{:.6})-({:.6},{:.6},{:.6})",
        wires_of_face(face).len(),
        ur.0, ur.1, vr.0, vr.1,
        b0.x(), b0.y(), b0.z(), b1.x(), b1.y(), b1.z()
    );
    for (wi, w) in wires_of_face(face).iter().enumerate() {
        let es = edges_of_wire(w);
        println!("{tag}  wire[{wi}] nEdges={}", es.len());
        for (ei, e) in es.iter().enumerate() {
            let (ef, el) = BRepTool::edge_parameters(e);
            println!(
                "{tag}   e[{ei}] ori={:?} deg={} par=[{:.15},{:.15}] dpar={:.15}",
                e.0.orientation(),
                BRepTool::is_degenerated(e),
                ef,
                el,
                el - ef
            );
        }
    }
}

/// Mesh just this one face through the real pipeline entry and report the
/// splitter verdict (`IsValid` is the gate that rejected F1691).
fn probe_splitter(face: &Face, label: &str) {
    let params = MeshParameters::default();
    let shape = face.0.clone();
    let mut model = match ModelBuilder::build_model(&shape, &params) {
        Ok(m) => m,
        Err(e) => {
            println!("SPLIT {label}: build_model failed: {e}");
            return;
        }
    };
    for i in 0..model.faces_nb() {
        let f = model.face(i).expect("face").clone();
        if let Some(s) = f.surface() {
            let (u0, u1) = s.u_range();
            let (v0, v1) = s.v_range();
            println!("SPLIT {label} modelface={i} wires={} urange=[{u0:.9},{u1:.9}] vrange=[{v0:.9},{v1:.9}]", f.wires().len());
        }
    }
}


fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = args.get(1).expect("usage: zz_seam_fix <file.step> <face-index>");
    let want: usize = args.get(2).expect("face index").parse().expect("index");
    let model = read_step_file(path).expect("read step");
    let shape = model.shapes[0].shape.clone();

    // `--all`: report how many faces change when the *Perform-sized* entry point
    // (`perform_fix_missing_seam`, which checks `fix_missing_seam_mode` and sets
    // myResult) is applied, as opposed to calling `fix_missing_seam` directly.
    if args.iter().any(|a| a == "--all") {
        let all = faces_of(&shape);
        let mut changed = 0usize;
        let mut grew: Vec<(usize, usize, usize)> = Vec::new();
        for (i, f) in all.iter().enumerate() {
            let before: usize = wires_of_face(f)
                .iter()
                .map(|w| edges_of_wire(w).len())
                .sum();
            let mut sff = ShapeFixFace::with_face(f);
            let res = sff.perform_fix_missing_seam();
            let after = res.as_ref().map(|r| {
                faces_of(r)
                    .iter()
                    .flat_map(|rf| wires_of_face(rf))
                    .map(|w| edges_of_wire(&w).len())
                    .sum::<usize>()
            });
            if let Some(a) = after {
                if a != before {
                    changed += 1;
                    grew.push((i, before, a));
                }
            }
        }
        println!(
            "PERFALL faces={} changed={}",
            all.len(),
            changed
        );
        for (i, b, a) in grew.iter() {
            println!("PERFALL  f={i} edges {b} -> {a}");
        }
        return;
    }

    // `--wireuv <f>`: for each wire of the named read face, print each edge's
    // pcurve evaluation at the parameters `CheckWire` uses, so a wire whose
    // summed displacement cancels can be diagnosed edge by edge.
    if args.iter().any(|a| a == "--wireuv") {
        use occt_topo::meshing::range_splitter::create_range_splitter;
        use occt_topo::shhealing::check_wire;
        let all = faces_of(&shape);
        let face = Face(all[want].0.clone());
        let (u_range, v_range) = match BRepTool::face_surface(&face) {
            Some(s) => {
                let (a, b) = s.u_range();
                let (c, d) = s.v_range();
                // Match FixMissingSeam's own range computation.
                let ur = (b - a).abs().min(1e100);
                let vr = (d - c).abs().min(1e100);
                (ur, vr)
            }
            None => (0.0, 0.0),
        };
        println!("WIREUV face={want} u_range={u_range:.9} v_range={v_range:.9}");
        {
            // T-99 checkpoint 1: raw inputs of the range guard (cxx:1781-1802).
            let s = BRepTool::face_surface(&face).unwrap();
            println!("WIREUV surf.u_range={:?}", s.u_range());
            println!("WIREUV surf.v_range={:?}", s.v_range());
            let (f_u1, f_u2, f_v1, f_v2) = occt_topo::brep_tools::uv_bounds(&face);
            println!("WIREUV uv_bounds face f_u=({f_u1},{f_u2}) f_v=({f_v1},{f_v2})");
        }
        let reg = occt_topo::tgeometry::GeometryRegistry::global();
        let fk = occt_topo::tgeometry::GeometryRegistry::shape_key(&face.0);
        for (wi, w) in wires_of_face(&face).iter().enumerate() {
            let es = edges_of_wire(w);
            println!("WIREUV wire[{wi}] nEdges={}", es.len());
            // T-99 decisive test: is this wire closed in 3D? If yes, the
            // per-edge UV displacement summing to (0,0) is CORRECT and the
            // problem is on the wire-PAIR SELECTION side, not the pcurve side.
            if let (Some(fv), Some(lv)) = (
                crate::first_vertex_of(&es),
                crate::last_vertex_of(&es),
            ) {
                let p0 = BRepTool::vertex_point(&fv);
                let p1 = BRepTool::vertex_point(&lv);
                let d = ((p1.x() - p0.x()).powi(2)
                    + (p1.y() - p0.y()).powi(2)
                    + (p1.z() - p0.z()).powi(2))
                .sqrt();
                println!(
                    "WIREUV wire[{wi}] first=({:.6},{:.6},{:.6}) last=({:.6},{:.6},{:.6}) gap={:.6} cl3d={}",
                    p0.x(), p0.y(), p0.z(), p1.x(), p1.y(), p1.z(), d,
                    if d < 1e-7 { "YES" } else { "NO" }
                );
                // T-99 decisive test 2: TShape identity of the wire's end
                // vertices. `ShapeFix_Wire` merges/reorders wires by vertex
                // identity (`TopoDS_Vertex::IsSame`), so two wires that meet
                // geometrically but do NOT share a TShape can never be joined.
                println!(
                    "WIREUV wire[{wi}] vkeys first={} last={}",
                    std::sync::Arc::as_ptr(&fv.0.tshape) as usize,
                    std::sync::Arc::as_ptr(&lv.0.tshape) as usize
                );
            }
            // Every wire's representative vertex TShapes, so cross-wire sharing
            // can be seen even when a wire has many edges.
            {
                let mut keys: Vec<usize> = Vec::new();
                for e in &es {
                    for v in [crate::first_vertex_of(&[e.clone()]), crate::last_vertex_of(&[e.clone()])]
                        .into_iter()
                        .flatten()
                    {
                        keys.push(std::sync::Arc::as_ptr(&v.0.tshape) as usize);
                    }
                }
                keys.sort_unstable();
                keys.dedup();
                println!("WIREUV wire[{wi}] all_vkeys({})={:?}", keys.len(), keys);
            }
            let mut sx = 0.0;
            let mut sy = 0.0;
            for (ei, e) in es.iter().enumerate() {
                let fkk = occt_topo::tgeometry::GeometryRegistry::shape_key(&face.0);
                if let Some((pc, f, l)) = occt_topo::boptools_2d::curve_on_surface_oriented(e, &face, true) {
                    let pf = pc.d0(f);
                    let pl = pc.d0(l);
                    sx += pl.x() - pf.x();
                    sy += pl.y() - pf.y();
                    println!(
                        "WIREUV   e[{ei}] ori={:?} par=[{f:.6},{l:.6}] pf=({:.6},{:.6}) pl=({:.6},{:.6}) d=({:.6},{:.6})",
                        e.0.orientation(),
                        pf.x(), pf.y(), pl.x(), pl.y(),
                        pl.x() - pf.x(), pl.y() - pf.y()
                    );
                } else {
                    println!("WIREUV   e[{ei}] no pcurve (fk={fk} fk2={fk} reg_n={})", reg.edge_pcurves(&e.0, fk).len());
                }
            }
            println!("WIREUV wire[{wi}] sum=({sx:.6},{sy:.6}) check_wire={:?}", check_wire(w, &face, u_range, v_range));
            let _ = create_range_splitter;
        }
        return;
    }

    let faces = faces_of(&shape);
    let face = Face(faces[want].0.clone());
    println!("--- BEFORE ---");
    describe("BEFORE", &face);
    if std::env::args().any(|a| a == "--ep") {
        for (wi, w) in occt_topo::topo_tools_full::wires_of_face(&face).iter().enumerate() {
            let es = occt_topo::topo_tools_full::edges_of_wire(w);
            let mut segs: Vec<String> = Vec::new();
            for e in &es {
                let f = occt_topo::shhealing::first_vertex(e)
                    .map(|v| { let q = BRepTool::vertex_point(&v); format!("({:.3},{:.3},{:.3})", q.x(), q.y(), q.z()) })
                    .unwrap_or_else(|| "none".to_string());
                let l = occt_topo::shhealing::last_vertex(e)
                    .map(|v| { let q = BRepTool::vertex_point(&v); format!("({:.3},{:.3},{:.3})", q.x(), q.y(), q.z()) })
                    .unwrap_or_else(|| "none".to_string());
                segs.push(format!("{f}->{l}"));
            }
            let zero = segs
                .iter()
                .filter(|t| { let p: Vec<&str> = t.split("->").collect(); p.len() == 2 && p[0] == p[1] })
                .count();
            println!("BEFOREP wire={wi} n={} zero={} {}", es.len(), zero, segs.join(" "));
        }
    }

    let mut sff = ShapeFixFace::with_face(&face);
    let ret = sff.fix_missing_seam();
    println!("fix_missing_seam returned {ret}");
    let Some(res) = sff.result.clone() else {
        println!("no result");
        return;
    };
    println!("result type={:?}", res.shape_type());
    let healed_faces = faces_of(&res);
    println!("result faces={}", healed_faces.len());
    for (i, hf) in healed_faces.iter().enumerate() {
        println!("--- HEALED face {i} ---");
        describe("HEALED", hf);
        // `--ep`: per-edge endpoints of the healed result (read-only instrument, §9.450).
        if std::env::args().any(|a| a == "--ep") {
            for (wi, w) in occt_topo::topo_tools_full::wires_of_face(hf).iter().enumerate() {
                let es = occt_topo::topo_tools_full::edges_of_wire(w);
                let segs: Vec<String> = es
                    .iter()
                    .map(|e| {
                        let f = occt_topo::shhealing::first_vertex(e)
                            .map(|v| {
                                let q = BRepTool::vertex_point(&v);
                                format!("({:.3},{:.3},{:.3})", q.x(), q.y(), q.z())
                            })
                            .unwrap_or_else(|| "none".to_string());
                        let l = occt_topo::shhealing::last_vertex(e)
                            .map(|v| {
                                let q = BRepTool::vertex_point(&v);
                                format!("({:.3},{:.3},{:.3})", q.x(), q.y(), q.z())
                            })
                            .unwrap_or_else(|| "none".to_string());
                        format!("{f}->{l}")
                    })
                    .collect();
                println!("HEALEDP face={i} wire={wi} n={} {}", es.len(), segs.join(" "));
            }
        }
        probe_splitter(hf, &format!("face{i}"));
    }
    let _ = (create_range_splitter, MeshFace::new, MeshModel::new, TopoBuilder::new);
}
