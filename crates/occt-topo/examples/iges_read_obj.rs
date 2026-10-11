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

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let dump = args.iter().any(|a| a == "--fstats");
    let paths: Vec<&String> = args.iter().filter(|a| a.as_str() != "--fstats").collect();
    if paths.len() != 2 {
        eprintln!("usage: iges_read_obj <iges> <obj> [--fstats]");
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
        println!("FACES n={} lin={lin:.6} ang={ang:.6}", faces.len());
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
                "FACE i={i} type={ty} u=[{u0:.6},{u1:.6}] v=[{v0:.6},{v1:.6}] uper={uper} vper={vper} wires={} edges={edge_ns:?}",
                wires.len()
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
        for st in inc.face_stats() {
            println!(
                "FSTAT face={} key={} mv={} mt={}",
                st.index, st.shape_key, st.vertices, st.triangles
            );
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
