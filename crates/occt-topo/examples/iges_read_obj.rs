//! Read one IGES file and write an OBJ mesh.
//!
//! Args: iges path, obj path. The mesh call is `brep_to_obj(&shape, 0.1)`,
//! the same call as `examples/export_data_obj.rs`.

use std::fs;
use std::path::PathBuf;
use std::process::ExitCode;

use occt_topo::brep_bnd_lib::shape_bnd_box;
use occt_topo::brep_exchange::{brep_to_obj, prs3d_get_deflection};
use occt_topo::brep_tool::BRepTool;
use occt_topo::iges_read::read_iges_shape;
use occt_topo::meshing::data_model::MeshStatus;
use occt_topo::meshing::incremental_mesh::IncrementalMesh;
use occt_topo::meshing::range_splitter::{classify_surface, SurfaceType};
use occt_topo::tgeometry::GeometryRegistry;
use occt_topo::topo_tools_full::{edges_of_wire, faces_of, wires_of_face};

fn tname(t: SurfaceType) -> &'static str {
    match t {
        SurfaceType::Plane => "Plane",
        SurfaceType::Cylinder => "Cylinder",
        SurfaceType::Cone => "Cone",
        SurfaceType::Sphere => "Sphere",
        SurfaceType::Torus => "Torus",
        SurfaceType::SurfaceOfRevolution => "Rev",
        SurfaceType::SurfaceOfExtrusion => "Extr",
        SurfaceType::BezierSurface => "Bezier",
        SurfaceType::BSplineSurface => "BSpline",
        SurfaceType::OffsetSurface => "Offset",
        SurfaceType::OtherSurface => "Other",
    }
}

/// OCCT `GeomAbs_SurfaceType` code of a port surface type, so the `type=` field
/// of a `BOUND` line matches `occt_probe --boundary` (`0` plane … `10` other).
fn tcode(t: SurfaceType) -> i32 {
    match t {
        SurfaceType::Plane => 0,
        SurfaceType::Cylinder => 1,
        SurfaceType::Cone => 2,
        SurfaceType::Sphere => 3,
        SurfaceType::Torus => 4,
        SurfaceType::BezierSurface => 5,
        SurfaceType::BSplineSurface => 6,
        SurfaceType::SurfaceOfRevolution => 7,
        SurfaceType::SurfaceOfExtrusion => 8,
        SurfaceType::OffsetSurface => 9,
        SurfaceType::OtherSurface => 10,
    }
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let dump = args.iter().any(|a| a == "--fstats");
    let bdump = args.iter().any(|a| a == "--bdump");
    // `--fdump <faceIndex>`: per-face node/UV/triangle dump of the port's own
    // mesh, in the format `occt_probe --facedump <faceIndex> <defl> <angle>`
    // prints, so the two can be diffed line by line. `<faceIndex>` is the
    // `MeshModel::face` index (the same index `FSTAT face=` reports), which is
    // paired to the probe's index by bbox, not by position.
    let fdump: Option<usize> = args
        .iter()
        .position(|a| a == "--fdump")
        .and_then(|i| args.get(i + 1))
        .and_then(|s| s.parse().ok());
    let mut skip: Vec<String> = vec!["--fstats".to_string(), "--bdump".to_string()];
    if let Some(i) = args.iter().position(|a| a == "--fdump") {
        skip.push("--fdump".to_string());
        if let Some(v) = args.get(i + 1) {
            skip.push(v.clone());
        }
    }
    let paths: Vec<&String> =
        args.iter().filter(|a| !skip.iter().any(|s| s == *a)).collect();
    if paths.len() != 2 {
        eprintln!("usage: iges_read_obj <iges> <obj> [--fstats] [--bdump] [--fdump <faceIndex>]");
        return ExitCode::from(2);
    }
    let iges = paths[0];
    let obj_path = PathBuf::from(paths[1]);
    let shape = match read_iges_shape(iges) {
        Ok(shape) => shape,
        Err(err) => {
            eprintln!("{err}");
            return ExitCode::from(1);
        }
    };
    if dump {
        let lin = prs3d_get_deflection(&shape, 0.1);
        let ang = 20.0_f64.to_radians();
        let faces = faces_of(&shape);
        println!("FACES n={} lin={lin:.17} ang={ang:.6}", faces.len());
        for (i, face) in faces.iter().enumerate() {
            let surf = BRepTool::face_surface(face);
            let ty = surf
                .as_ref()
                .map(|s| tname(classify_surface(s.as_ref())))
                .unwrap_or("none");
            let (u0, u1, v0, v1) = surf
                .as_ref()
                .map(|s| (s.u_range().0, s.u_range().1, s.v_range().0, s.v_range().1))
                .unwrap_or((0.0, 0.0, 0.0, 0.0));
            let (uper, vper) = surf
                .as_ref()
                .map(|s| (s.is_u_periodic(), s.is_v_periodic()))
                .unwrap_or((false, false));
            let wires = wires_of_face(face);
            let edge_ns: Vec<usize> = wires.iter().map(|w| edges_of_wire(w).len()).collect();
            println!(
                "FACE i={i} type={ty} u=[{u0:.6},{u1:.6}] v=[{v0:.6},{v1:.6}] uper={uper} vper={vper} wires={} edges={edge_ns:?} bbox={}",
                wires.len(),
                {
                    let bb = shape_bnd_box(&face.0);
                    let (b0, b1) = (bb.corner_min(), bb.corner_max());
                    format!("({:.6},{:.6},{:.6})-({:.6},{:.6},{:.6})", b0.x(), b0.y(), b0.z(), b1.x(), b1.y(), b1.z())
                }
            );
            if let Some(surf) = surf.as_ref() {
                if let Some(axis) = surf.revolution_axis() {
                    let loc = axis.location();
                    let dir = axis.direction();
                    for k in 0..5 {
                        let t = k as f64 / 4.0;
                        let v = v0 + t * (v1 - v0);
                        let p = surf.d0(0.0, v);
                        let rel = occt_core::gp::GpVec::from_pnts(&loc, &p);
                        let along = dir.xyz().multiplied(rel.xyz().dot(dir.xyz()));
                        let rad = rel.xyz().subtracted(&along).modulus();
                        let iso_r = surf
                            .v_iso_curve(v)
                            .and_then(|c| c.circle_radius())
                            .unwrap_or(-1.0);
                        println!(
                            "  VISO v={v:.6} rad3d={rad:.6} iso_circle_r={iso_r:.6} p=({:.6},{:.6},{:.6})",
                            p.x(),
                            p.y(),
                            p.z()
                        );
                    }
                }
            }
            let fk = GeometryRegistry::shape_key(&face.0);
            let reg = GeometryRegistry::global();
            for (wi, w) in wires.iter().enumerate() {
                for (ei, e) in edges_of_wire(w).iter().enumerate() {
                    let (f3, l3) = BRepTool::edge_parameters(e);
                    let deg = BRepTool::is_degenerated(e);
                    let c3 = BRepTool::edge_curve(e)
                        .map(|c| format!("{}", std::any::type_name_of_val(c.as_ref())))
                        .unwrap_or_else(|| "None".into());
                    let pcs = reg.edge_pcurves(&e.0, fk);
                    let (pf, pl) = reg.pcurve_range(&e.0, fk).unwrap_or((f3, l3));
                    let pct = pcs
                        .first()
                        .map(|c| format!("{}", std::any::type_name_of_val(c.as_ref())))
                        .unwrap_or_else(|| "None".into());
                    println!(
                        "  E w={wi} e={ei} deg={deg} c3d={c3} range3=[{f3:.6},{l3:.6}] pc={pct} pcr=[{pf:.6},{pl:.6}] ori={:?}",
                        e.0.orientation()
                    );
                }
            }
        }
        let inc = IncrementalMesh::from_deflection(&shape, lin, false, ang);
        // Pair a face with the OCCT probe's `FSTAT_OCC` row by its bbox: the mesh
        // model's face order (`BRepLib::ReverseSortFaces`) is not the TopExp order.
        let params = occt_topo::meshing::parameters::MeshParameters {
            deflection: lin,
            angle: ang,
            ..Default::default()
        };
        let model = occt_topo::meshing::model_builder::ModelBuilder::build_model(&shape, &params)
            .expect("build_model");
        for st in inc.face_stats() {
            let bb = model
                .face(st.index)
                .ok()
                .map(|f| shape_bnd_box(&f.face().0));
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
    }
    if bdump {
        // Boundary stage only: `IMeshTools_Context` order `BuildModel` ->
        // `DiscretizeEdges` -> `HealModel` -> `PreProcessModel`, then dump what
        // `BRepMesh_BaseMeshAlgo::initDataStructure` (`cxx:87-125`) would read.
        // Output mirrors `occt_probe --boundary <dir>/boundary.txt` so the two
        // can be diffed line by line (faces paired by bbox, not by index).
        let lin = prs3d_get_deflection(&shape, 0.1);
        let ang = 20.0_f64.to_radians();
        let params = occt_topo::meshing::parameters::MeshParameters {
            deflection: lin,
            angle: ang,
            ..Default::default()
        };
        let model = match IncrementalMesh::build_boundary_model(&shape, &params) {
            Ok(m) => m,
            Err(err) => {
                eprintln!("bdump: {err}");
                return ExitCode::from(1);
            }
        };
        println!(
            "BOUNDARY defl={lin:.6} angle={ang:.6} faces={} edges={}",
            model.faces_nb(),
            model.edges_nb()
        );
        for f in 0..model.faces_nb() {
            let face = match model.face(f) {
                Ok(face) => face,
                Err(_) => continue,
            };
            let (b0, b1) = match shape_bnd_box(&face.face().0).get() {
                Some((x0, x1, y0, y1, z0, z1)) => (
                    occt_core::gp::GpPnt::new(x0, y0, z0),
                    occt_core::gp::GpPnt::new(x1, y1, z1),
                ),
                None => (occt_core::gp::GpPnt::zero(), occt_core::gp::GpPnt::zero()),
            };
            let ty = face
                .surface()
                .as_ref()
                .map(|s| tcode(classify_surface(s.as_ref())))
                .unwrap_or(10);
            let mut pts = 0usize;
            for w in 0..face.wires_nb() {
                let Ok(wi) = face.wire(w) else { continue };
                let Ok(dw) = model.wire(wi) else { continue };
                for e in 0..dw.edges_nb() {
                    let Ok(ei) = dw.edge(e) else { continue };
                    let Ok(de) = model.edge(ei) else { continue };
                    for pc in de.pcurves_for(f) {
                        pts += de.pcurve(pc).map(|p| p.parameters_nb()).unwrap_or(0);
                    }
                }
            }
            println!(
                "BOUND f={f} type={ty} wires={} pts={pts} bbox=({:.6},{:.6},{:.6})-({:.6},{:.6},{:.6})",
                face.wires_nb(),
                b0.x(),
                b0.y(),
                b0.z(),
                b1.x(),
                b1.y(),
                b1.z()
            );
            for w in 0..face.wires_nb() {
                let Ok(wi) = face.wire(w) else { continue };
                let Ok(dw) = model.wire(wi) else { continue };
                println!(
                    "W {w} edges={} selfint={} open={}",
                    dw.edges_nb(),
                    dw.is_status(MeshStatus::SELF_INTERSECTING_WIRE) as i32,
                    dw.is_status(MeshStatus::OPEN_WIRE) as i32
                );
                for e in 0..dw.edges_nb() {
                    let Ok(ei) = dw.edge(e) else { continue };
                    let Ok(de) = model.edge(ei) else { continue };
                    for pc in de.pcurves_for(f) {
                        let Ok(p) = de.pcurve(pc) else { continue };
                        println!(
                            "E w={w} e={e} pcId={pc} ori={} fwd={} n={}",
                            p.orientation() as i32,
                            p.is_forward() as i32,
                            p.parameters_nb()
                        );
                        for k in 0..p.parameters_nb() {
                            if let Ok(p2) = p.get_point(k) {
                                println!("P {:.17e} {:.17e}", p2.x(), p2.y());
                            }
                        }
                    }
                }
            }
        }
    }
    if let Some(want) = fdump {
        // Per-face dump of the port's own `FaceTriangulation`, in the probe's
        // `--facedump` format. Faces are ordered by `MeshModel::face` index, so
        // the probe row for the same face is found by bbox (FSTAT / FSTAT_OCC),
        // not by position.
        let lin = prs3d_get_deflection(&shape, 0.1);
        let ang = 20.0_f64.to_radians();
        let params = occt_topo::meshing::parameters::MeshParameters {
            deflection: lin,
            angle: ang,
            ..Default::default()
        };
        let model = match occt_topo::meshing::model_builder::ModelBuilder::build_model(&shape, &params)
        {
            Ok(m) => m,
            Err(err) => {
                eprintln!("fdump: {err}");
                return ExitCode::from(1);
            }
        };
        let inc = IncrementalMesh::from_deflection(&shape, lin, false, ang);
        println!("FDUMP_PORT defl={lin:.17} angle={ang:.6} faces={}", inc.face_triangulations().len());
        for ft in inc.face_triangulations() {
            if ft.face_index != want {
                continue;
            }
            let bb = model.face(ft.face_index).ok().map(|f| shape_bnd_box(&f.face().0));
            let (b0, b1) = match bb {
                Some(b) => (b.corner_min(), b.corner_max()),
                None => (occt_core::gp::GpPnt::zero(), occt_core::gp::GpPnt::zero()),
            };
            println!(
                "FDUMP face={} nodes={} tris={} bbox=({:.6},{:.6},{:.6})-({:.6},{:.6},{:.6})",
                ft.face_index,
                ft.vertices.len(),
                ft.triangles.len(),
                b0.x(),
                b0.y(),
                b0.z(),
                b1.x(),
                b1.y(),
                b1.z()
            );
            for (i, (uv, p)) in ft.uv.iter().zip(ft.vertices.iter()).enumerate() {
                println!(
                    "FDUMP n i={} uv=({:.9},{:.9}) p=({:.9},{:.9},{:.9})",
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
                    "FDUMP t i={} ({},{},{})",
                    i + 1,
                    t.n0 + 1,
                    t.n1 + 1,
                    t.n2 + 1
                );
            }
        }
    }
    let obj = brep_to_obj(&shape, 0.1);
    let nv = obj.lines().filter(|line| line.starts_with("v ")).count();
    let nf = obj.lines().filter(|line| line.starts_with("f ")).count();
    let bbox = match shape_bnd_box(&shape).get() {
        Some((x0, x1, y0, y1, z0, z1)) => {
            format!("exact=[{x0},{y0},{z0}]~[{x1},{y1},{z1}]")
        }
        None => "exact=void".to_string(),
    };
    if let Some(parent) = obj_path.parent() {
        if !parent.as_os_str().is_empty() {
            if let Err(err) = fs::create_dir_all(parent) {
                eprintln!("create {}: {err}", parent.display());
                return ExitCode::from(1);
            }
        }
    }
    if let Err(err) = fs::write(&obj_path, &obj) {
        eprintln!("write {}: {err}", obj_path.display());
        return ExitCode::from(1);
    }
    println!("v={nv} f={nf} {bbox}");
    if nv == 0 || nf == 0 {
        return ExitCode::from(1);
    }
    ExitCode::SUCCESS
}
