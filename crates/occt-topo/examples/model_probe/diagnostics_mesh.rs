//! Mesh diagnostics (--raw, --fixms, --fstats, --eshare).
use super::tname;
use super::Ctx;
use occt_topo::brep_bnd_lib;
use occt_topo::brep_exchange::prs3d_get_deflection;
use occt_topo::brep_tool::BRepTool;
use occt_topo::meshing::incremental_mesh::IncrementalMesh;
use occt_topo::meshing::range_splitter::mesh_surface_type;
use occt_topo::tgeometry::GeometryRegistry;
use occt_topo::topo_tools::{faces_of, wires_of_face};

pub(crate) fn run(ctx: &Ctx) -> bool {
    let shape = ctx.shape.clone();
    let mode_fixms = ctx.mode_fixms;
    let mode_fstats = ctx.mode_fstats;
    if std::env::args().any(|a| a == "--raw") {
        for (i, f) in faces_of(&shape).iter().enumerate() {
            let wires = wires_of_face(f);
            if wires.len() != 2 { continue; }
            let surf = BRepTool::face_surface(f);
            let tn = surf.as_ref().map(|s| tname(mesh_surface_type(s.as_ref()))).unwrap_or("none");
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
        return true;
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
            sff.fix_periodic_degenerated();
            let r = sff.fix_missing_seam();
            let res = sff.result.clone();
            let desc = res.as_ref().map(|t| format!("{:?} faces={}", t.shape_type(), faces_of(t).len())).unwrap_or_else(|| "none".into());
            let bb = brep_bnd_lib::shape_bnd_box(&f.0);
            let (mn, mx) = (bb.corner_min(), bb.corner_max());
            println!("FIXMS face={i} uper={uper} vper={vper} before_wires={before_wires} ret={r} result={desc} bbox={:.4},{:.4},{:.4},{:.4},{:.4},{:.4}", mn.x(), mn.y(), mn.z(), mx.x(), mx.y(), mx.z());
        }
        return true;
    }
    if mode_fstats {
        // `FaceMeshStat` carries `index` (the mesh model's own face index) and
        // `shape_key` precisely because a positional lookup is off by the number
        // of skipped faces (`discret_root.rs:81-88`). Pair through `stat.index`
        // and resolve the bbox from the same pipeline's model, never from
        // `faces_of(&shape)[i]` (that ordering is a different traversal).
        let defl_override = std::env::var("FSTAT_DEFL").ok().and_then(|s| s.parse::<f64>().ok());
        let ang_override = std::env::var("FSTAT_ANGLE").ok().and_then(|s| s.parse::<f64>().ok());
        let params = occt_topo::meshing::parameters::MeshParameters {
            deflection: defl_override.unwrap_or_else(|| prs3d_get_deflection(&shape, 0.1)),
            angle: ang_override.unwrap_or_else(|| 20.0_f64.to_radians()),
            ..Default::default()
        };
        let inc = IncrementalMesh::from_deflection(&shape, params.deflection, false, params.angle);
        let stats = inc.face_stats().to_vec();
        let model = occt_topo::meshing::model_builder::ModelBuilder::build_model(&shape, &params)
            .expect("build_model");
        for st in &stats {
            let bb = model
                .face(st.index)
                .ok()
                .map(|f| occt_topo::brep_bnd_lib::shape_bnd_box(&f.face().0));
            let (b0, b1) = match bb {
                Some(b) => (b.corner_min(), b.corner_max()),
                None => (occt_core::gp::GpPnt::zero(), occt_core::gp::GpPnt::zero()),
            };
            println!(
                "FSTAT face={} key={} mv={} mt={} bbox=({:.6},{:.6},{:.6})-({:.6},{:.6},{:.6})",
                st.index,
                st.shape_key,
                st.vertices,
                st.triangles,
                b0.x(),
                b0.y(),
                b0.z(),
                b1.x(),
                b1.y(),
                b1.z()
            );
        }
        return true;
    }
    if std::env::args().any(|a| a == "--eshare") {
        // Edge-sharing census: per model edge, how many pcurves (faces) reference
        // it and how many pcurve samples it carries. An edge with one pcurve is
        // not shared between two faces even though a coincident edge may exist.
        let params = occt_topo::meshing::parameters::MeshParameters {
            deflection: prs3d_get_deflection(&shape, 0.1),
            angle: 20.0_f64.to_radians(),
            ..Default::default()
        };
        let model = occt_topo::meshing::model_builder::ModelBuilder::build_model(&shape, &params)
            .expect("build_model");
        let mut unshared = 0usize;
        let mut hist = std::collections::BTreeMap::new();
        for i in 0..model.edges_nb() {
            let e = model.edge(i).expect("edge");
            let npc = e.pcurves_nb();
            *hist.entry(npc).or_insert(0usize) += 1;
            if npc <= 1 {
                unshared += 1;
            }
            let bb = brep_bnd_lib::shape_bnd_box(e.edge());
            let (b0, b1) = (bb.corner_min(), bb.corner_max());
            let faces: Vec<usize> = (0..npc)
                .filter_map(|p| e.pcurve(p).ok().map(|pc| pc.face()))
                .collect();
            println!(
                "ESHARE edge={} pcurves={} faces={:?} bbox=({:.6},{:.6},{:.6})-({:.6},{:.6},{:.6})",
                i,
                npc,
                faces,
                b0.x(),
                b0.y(),
                b0.z(),
                b1.x(),
                b1.y(),
                b1.z()
            );
        }
        println!("ESHARE total={} unshared={} hist={:?}", model.edges_nb(), unshared, hist);
        return true;
    }
    if std::env::args().any(|a| a == "--fsurf") {
        let params = occt_topo::meshing::parameters::MeshParameters {
            deflection: prs3d_get_deflection(&shape, 0.1),
            angle: 20.0_f64.to_radians(),
            ..Default::default()
        };
        let inc = IncrementalMesh::from_deflection(&shape, params.deflection, false, params.angle);
        let stats = inc.face_stats().to_vec();
        let model = occt_topo::meshing::model_builder::ModelBuilder::build_model(&shape, &params)
            .expect("build_model");
        for st in &stats {
            let Ok(f) = model.face(st.index) else { continue };
            let surf = f.surface();
            let tn = surf.as_ref().map(|s| tname(mesh_surface_type(s.as_ref()))).unwrap_or("none");
            let (u0, u1, v0, v1) = surf
                .as_ref()
                .map(|s| (s.u_range().0, s.u_range().1, s.v_range().0, s.v_range().1))
                .unwrap_or((0.0, 0.0, 0.0, 0.0));
            let (uper, vper) = surf
                .as_ref()
                .map(|s| (s.is_u_periodic(), s.is_v_periodic()))
                .unwrap_or((false, false));
            let bb = brep_bnd_lib::shape_bnd_box(&f.face().0);
            let (b0, b1) = (bb.corner_min(), bb.corner_max());
            println!(
                "FSURF face={} key={} type={} mv={} mt={} wires={} urange=[{u0:.6},{u1:.6}] vrange=[{v0:.6},{v1:.6}] uper={uper} vper={vper} bbox=({:.6},{:.6},{:.6})-({:.6},{:.6},{:.6})",
                st.index,
                st.shape_key,
                tn,
                st.vertices,
                st.triangles,
                f.wires().len(),
                b0.x(),
                b0.y(),
                b0.z(),
                b1.x(),
                b1.y(),
                b1.z()
            );
        }
        return true;
    }
    if let Some(pos) = std::env::args().position(|a| a == "--tridump") {
        // `--tridump <TopExp faceIdx>` mirrors `occt_probe --facedump <modelIdx>`
        // for the same face: node UV/3D and triangle triples.
        let want: usize = std::env::args()
            .nth(pos + 1)
            .and_then(|s| s.parse().ok())
            .expect("--tridump <faceIdx>");
        let f = faces_of(&shape)[want].clone();
        let key = GeometryRegistry::shape_key(&f.0);
        let defl = std::env::var("FSTAT_DEFL")
            .ok()
            .and_then(|s| s.parse::<f64>().ok())
            .unwrap_or_else(|| prs3d_get_deflection(&shape, 0.1));
        let ang = std::env::var("FSTAT_ANGLE")
            .ok()
            .and_then(|s| s.parse::<f64>().ok())
            .unwrap_or_else(|| 20.0_f64.to_radians());
        let inc = IncrementalMesh::from_deflection(&shape, defl, false, ang);
        let Some(st) = inc.face_stats().iter().find(|s| s.shape_key == key) else {
            println!("TRIDUMP face={want} key={key} no-stats");
            return true;
        };
        let Some(ft) = inc
            .face_triangulations()
            .iter()
            .find(|t| t.face_index == st.index)
        else {
            println!("TRIDUMP face={want} key={key} model={} no-triangulation", st.index);
            return true;
        };
        println!(
            "TRIDUMP face={want} key={key} model={} nodes={} tris={}",
            st.index,
            ft.vertices.len(),
            ft.triangles.len()
        );
        for (i, (p, uv)) in ft.vertices.iter().zip(ft.uv.iter()).enumerate() {
            println!(
                "TRIDUMP n i={} uv=({:.9},{:.9}) p=({:.9},{:.9},{:.9})",
                i + 1,
                uv.x(),
                uv.y(),
                p.x(),
                p.y(),
                p.z()
            );
        }
        for (i, t) in ft.triangles.iter().enumerate() {
            println!(
                "TRIDUMP t i={} ({},{},{})",
                i + 1,
                t.n0 + 1,
                t.n1 + 1,
                t.n2 + 1
            );
        }
        return true;
    }
    false
}
