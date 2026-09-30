//! TEMPORARY diagnostic probe (delete after T-69).
//!
//! Runs the real `IncrementalMesh` pipeline over a STEP file and, for every
//! face, dumps the (u,v) point set that `NodeInsertionMeshAlgo::perform` feeds
//! into `RangeSplitter::add_point` (the port's `collectWirePoints` equivalent),
//! so it can be diffed against the OCCT probe's `--uv <face>` output.
//!
//! Usage: `cargo run --example zz_uv_feed -- <file.step> [--sum] [--dump <path>]`
//!
//! `--dump <path>` re-runs the pipeline a second time and writes the flat OBJ
//! triangle soup of that run to `path`, so the dump and the mesh always come
//! from the same code path.
use occt_topo::meshing::incremental_mesh::IncrementalMesh;
use occt_topo::step::read_step_file;
use std::io::Write;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = args.get(1).expect("usage: zz_uv_feed <file.step> [--sum] [--dump <path>]");
    let model = read_step_file(path).expect("read step");
    let shape = model.shapes[0].shape.clone();
    // T-99: `--lin <d>` overrides the port's deflection so the port's mesh can be
    // compared with the OCCT probe at the SAME deflection (`probe --mesh <d> <a>`).
    // The default (`prs3d_get_deflection(shape, 0.1)`) is what the OBJ gates use.
    let lin: f64 = args
        .iter()
        .position(|a| a == "--lin")
        .and_then(|i| args.get(i + 1))
        .and_then(|s| s.parse().ok())
        .unwrap_or_else(|| occt_topo::brep_exchange::prs3d_get_deflection(&shape, 0.1));

    let mode_ids = args.iter().any(|a| a == "--ids");
    std::env::set_var("OCCT_TOPO_DUMP_UV", "1");
    if !mode_ids && args.iter().any(|a| a == "--sum") {
        std::env::set_var("OCCT_TOPO_DUMP_UV_SUM", "1");
    }
    let inc = IncrementalMesh::from_deflection(&shape, lin, false, 20.0_f64.to_radians());
    eprintln!("zz_uv_feed done lin={lin:.6}");

    // `--fixall`: run the port's `ShapeFixFace::fix_missing_seam` over every face
    // of the shape (no reader change), rebuild the discrete model from the
    // healed shape and report how many faces now have a zero-width added-point
    // range. Measurement only -- the library is untouched.
    if args.iter().any(|a| a == "--fixall") {
        use occt_topo::brep_bnd_lib;
        use occt_topo::meshing::data_model::MeshModel;
        use occt_topo::meshing::edge_discret::EdgeDiscret;
        use occt_topo::meshing::incremental_mesh::IncrementalMesh;
        use occt_topo::meshing::model_builder::ModelBuilder;
        use occt_topo::meshing::parameters::MeshParameters;
        use occt_topo::shape::TopoShape;
        use occt_topo::shhealing::ShapeFixFace;
        use occt_topo::topo_tools_full::{faces_of, wires_of_face};

        // For each face, ask the port's seam fix whether it would change the
        // loop, and record the edge counts before/after. Read-only measurement.
        let mut n_fix = 0usize;
        let mut grown: Vec<(usize, usize, usize)> = Vec::new();
        for (i, f) in faces_of(&shape).iter().enumerate() {
            let before: usize = wires_of_face(f).iter().map(|w| occt_topo::topo_tools_full::edges_of_wire(w).len()).sum();
            let mut sff = ShapeFixFace::with_face(f);
            let ret = sff.fix_missing_seam();
            let after = sff.result.as_ref().map(|res| {
                if res.shape_type() == occt_topo::abs::ShapeType::Face {
                    let rf = occt_topo::shape::Face(res.clone());
                    wires_of_face(&rf)
                        .iter()
                        .map(|w| occt_topo::topo_tools_full::edges_of_wire(w).len())
                        .sum::<usize>()
                } else {
                    0
                }
            });
            if ret {
                n_fix += 1;
                if let Some(a) = after {
                    if a != before {
                        grown.push((i, before, a));
                    }
                }
            }
        }
        println!("FIXALL faces={} fix_returns_true={} edge_count_changed={}", faces_of(&shape).len(), n_fix, grown.len());
        for (i, b, a) in grown.iter() {
            println!("FIXALL  f={i} edges {b} -> {a}");
        }
        return;
    }

    // `--fixms <f>`: run the port's `ShapeFixFace::perform_fix_missing_seam` on
    // one read face and report the resulting wire/edge structure, so it can be
    // compared with OCCT's healed import for the same face.
    if args.iter().any(|a| a == "--fixms") {
        use occt_topo::brep_bnd_lib;
        use occt_topo::brep_tool::BRepTool;
        use occt_topo::shhealing::ShapeFixFace;
        use occt_topo::topo_tools_full::{edges_of_wire, faces_of, wires_of_face};
        let list: Vec<usize> = args
            .iter()
            .enumerate()
            .filter(|(_, a)| *a == "--fixms")
            .filter_map(|(i, _)| args.get(i + 1))
            .filter_map(|s| s.parse().ok())
            .collect();
        for (i, f) in faces_of(&shape).iter().enumerate() {
            if !list.is_empty() && !list.contains(&i) {
                continue;
            }
            let bb = brep_bnd_lib::shape_bnd_box(&f.0);
            let (b0, b1) = (bb.corner_min(), bb.corner_max());
            let before: Vec<usize> = wires_of_face(f).iter().map(|w| edges_of_wire(w).len()).collect();
            let mut sff = ShapeFixFace::with_face(f);
            let ret = sff.fix_missing_seam();
            let res = sff.result.clone();
            let after = res.as_ref().map(|t| {
                let fs = faces_of(t);
                let mut v = Vec::new();
                for ff in &fs {
                    for w in wires_of_face(ff) {
                        v.push(edges_of_wire(&w).len());
                    }
                }
                v
            });
            println!(
                "FIXMS f={i} ret={ret} before_wires={before:?} bbox=({:.6},{:.6},{:.6})-({:.6},{:.6},{:.6}) surf={} uper={} vper={} nfaces={} after_edges={after:?}",
                b0.x(), b0.y(), b0.z(), b1.x(), b1.y(), b1.z(),
                BRepTool::face_surface(f).map(|s| format!("{:?}", occt_topo::meshing::range_splitter::classify_surface(s.as_ref()))).unwrap_or("none".into()),
                BRepTool::face_surface(f).map(|s| s.is_u_periodic()).unwrap_or(false),
                BRepTool::face_surface(f).map(|s| s.is_v_periodic()).unwrap_or(false),
                res.as_ref().map(|t| faces_of(t).len()).unwrap_or(0),
            );
        }
        return;
    }

    // `--model`: rebuild the discrete model exactly as `IncrementalMesh` does
    // (ModelBuilder + EdgeDiscret) and dump per-face wires/edges so the port's
    // pre-mesh face structure can be compared with OCCT's raw / healed import.
    if args.iter().any(|a| a == "--model") {
        use occt_topo::brep_bnd_lib;
        use occt_topo::brep_tool::BRepTool;
        use occt_topo::meshing::data_model::MeshModel;
        use occt_topo::meshing::edge_discret::EdgeDiscret;
        use occt_topo::meshing::model_builder::ModelBuilder;
        use occt_topo::meshing::parameters::MeshParameters;
        let want: Option<usize> = args
            .iter()
            .position(|a| a == "--model")
            .and_then(|i| args.get(i + 1))
            .and_then(|s| s.parse().ok());
        let mut model = ModelBuilder::build_model(&shape, &MeshParameters::default())
            .expect("build_model");
        for i in 0..model.faces_nb() {
            if let Some(w) = want {
                if i != w {
                    continue;
                }
            }
            let face = model.face(i).expect("face").clone();
            let bb = brep_bnd_lib::shape_bnd_box(&face.face().0);
            let (b0, b1) = (bb.corner_min(), bb.corner_max());
            let (ur, vr) = face
                .surface()
                .map(|s| (s.u_range(), s.v_range()))
                .unwrap_or(((0.0, 0.0), (0.0, 0.0)));
            println!(
                "MODEL f={i} wires={} urange=[{:.17},{:.17}] vrange=[{:.17},{:.17}] bbox=({:.6},{:.6},{:.6})-({:.6},{:.6},{:.6})",
                face.wires().len(),
                ur.0, ur.1, vr.0, vr.1,
                b0.x(), b0.y(), b0.z(), b1.x(), b1.y(), b1.z()
            );
            for (wi, &wire_index) in face.wires().iter().enumerate() {
                let w = model.wire(wire_index).expect("wire");
                println!("MODEL  wire[{wi}] nEdges={}", w.edges_nb());
                for j in 0..w.edges_nb() {
                    let ei = w.edge(j).expect("edge index");
                    let ori = w.edge_orientation(j).expect("ori");
                    let e = model.edge(ei).expect("edge");
                    let deg = e.degenerated();
                    let (ef, el) = (e.first_parameter(), e.last_parameter());
                    let curve_desc = match e.curve() {
                        None => "none".to_string(),
                        Some(c) => {
                            if c.is_line() {
                                "Line".to_string()
                            } else if c.gp_circ().is_some() {
                                "Circle".to_string()
                            } else if c.bspline_poles().is_some() {
                                "BSpline".to_string()
                            } else {
                                "Other".to_string()
                            }
                        }
                    };
                    println!(
                        "MODEL   e[{j}] edge={ei} ori={ori:?} deg={deg} curve={curve_desc} par=[{ef:.15},{el:.15}] dpar={:.15} npc={} n3d={}",
                        el - ef,
                        e.pcurves_for(i).len(),
                        e.discretization().points().len()
                    );
                }
            }
        }
        let _ = BRepTool::face_tolerance;
        // `--model <f> --surf`: evaluate the face surface on a UV grid, so the
        // actual 3D patch the face stands on is visible (a full sphere vs a band).
        if args.iter().any(|a| a == "--surf") {
            if let Some(w) = want {
                let face = model.face(w).expect("face");
                if let Some(s) = face.surface() {
                    let (u0, u1) = s.u_range();
                    let (v0, v1) = s.v_range();
                    for i in 0..=4 {
                        let v = v0 + (v1 - v0) * (i as f64) / 4.0;
                        let p = s.d0(u0 + (u1 - u0) * 0.25, v);
                        println!(
                            "SURF v={v:.9} p=({:.9},{:.9},{:.9})",
                            p.x(), p.y(), p.z()
                        );
                    }
                    for i in 0..=4 {
                        let u = u0 + (u1 - u0) * (i as f64) / 4.0;
                        let p = s.d0(u, v0);
                        println!("SURF u={u:.9} v=v0 p=({:.9},{:.9},{:.9})", p.x(), p.y(), p.z());
                    }
                }
            }
        }
        return;
    }

    if args.iter().any(|a| a == "--sph") {
        use occt_topo::brep_bnd_lib;
        use occt_topo::brep_tool::BRepTool;
        use occt_topo::meshing::range_splitter::{classify_surface, SurfaceType};
        use occt_topo::topo_tools_full::{edges_of_wire, faces_of, wires_of_face};
        for (i, f) in faces_of(&shape).iter().enumerate() {
            let Some(surf) = BRepTool::face_surface(f) else { continue };
            if classify_surface(surf.as_ref()) != SurfaceType::Sphere {
                continue;
            }
            let bb = brep_bnd_lib::shape_bnd_box(&f.0);
            let (b0, b1) = (bb.corner_min(), bb.corner_max());
            let key = occt_topo::tgeometry::GeometryRegistry::shape_key(&f.0);
            let st = inc.face_stats().iter().find(|s| s.shape_key == key);
            let (u0, u1) = surf.u_range();
            let (v0, v1) = surf.v_range();
            println!(
                "SPH f={i} wires={} urange=[{u0:.17},{u1:.17}] vrange=[{v0:.17},{v1:.17}] uper={} vper={} mt={} bbox=({:.6},{:.6},{:.6})-({:.6},{:.6},{:.6})",
                wires_of_face(f).len(),
                surf.is_u_periodic(),
                surf.is_v_periodic(),
                st.map(|s| s.triangles).unwrap_or(usize::MAX),
                b0.x(), b0.y(), b0.z(), b1.x(), b1.y(), b1.z()
            );
            for (wi, w) in wires_of_face(f).iter().enumerate() {
                let es = edges_of_wire(w);
                println!("SPH  wire[{wi}] nEdges={}", es.len());
                for (ei, e) in es.iter().enumerate() {
                    let (ef, el) = BRepTool::edge_parameters(e);
                    let closed = BRepTool::is_closed_edge_face(e, f);
                    println!(
                        "SPH   e[{ei}] ori={:?} deg={} closed={} par=[{ef:.15},{el:.15}] dpar={:.15}",
                        e.0.orientation(),
                        BRepTool::is_degenerated(e),
                        closed,
                        el - ef
                    );
                }
            }
        }
        return;
    }

    if mode_ids {
        use occt_topo::brep_bnd_lib;
        use occt_topo::brep_tool::BRepTool;
        use occt_topo::topo_tools_full::{faces_of, wires_of_face};
        for (i, f) in faces_of(&shape).iter().enumerate() {
            let bb = brep_bnd_lib::shape_bnd_box(&f.0);
            let (b0, b1) = (bb.corner_min(), bb.corner_max());
            let key = occt_topo::tgeometry::GeometryRegistry::shape_key(&f.0);
            let st = inc.face_stats().iter().find(|s| s.shape_key == key);
            println!(
                "PORTID f={i} key={key} wires={} mt={} bbox=({:.6},{:.6},{:.6})-({:.6},{:.6},{:.6})",
                wires_of_face(f).len(),
                st.map(|s| s.triangles).unwrap_or(usize::MAX),
                b0.x(), b0.y(), b0.z(), b1.x(), b1.y(), b1.z()
            );
            let _ = BRepTool::face_tolerance(f);
        }
        return;
    }

    if let Some(pos) = args.iter().position(|a| a == "--dump") {
        let out_path = args.get(pos + 1).expect("--dump <path>");
        std::env::set_var("OCCT_TOPO_DUMP_UV_SUM", "1");
        let inc2 = IncrementalMesh::from_deflection(&shape, lin, false, 20.0_f64.to_radians());
        let mesh = inc2.mesh().expect("mesh");
        let mut f = std::fs::File::create(out_path).expect("create dump");
        for v in &mesh.vertices {
            writeln!(f, "v {} {} {}", v.x(), v.y(), v.z()).unwrap();
        }
        for t in &mesh.triangles {
            writeln!(f, "f {} {} {}", t.n0 + 1, t.n1 + 1, t.n2 + 1).unwrap();
        }
        eprintln!(
            "dumped {} verts {} tris (first run: {} verts)",
            mesh.vertices.len(),
            mesh.triangles.len(),
            inc.mesh().map(|m| m.vertices.len()).unwrap_or(0)
        );
    }
}
