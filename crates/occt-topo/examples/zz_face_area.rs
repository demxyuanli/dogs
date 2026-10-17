//! TEMPORARY probe: per-face `BRepGProp::SurfaceProperties` area, for
//! cross-checking `brep_gprop_full` against `specs/occt_probe/zz_face_area_probe.cpp`.
//!
//! Faces are emitted in `TopExp_Explorer(shape, TopAbs_FACE)` order (the
//! depth-first `solid -> shell -> face` descent), the same order the OCCT probe
//! prints, so the two dumps can be diffed line by line.
//!
//! Set `ZZFACE=<n>` to also dump the boundary arcs of that face.
//!
//! Usage: cargo run --example zz_face_area -- occ/T0M.stp [minArea maxArea]

use std::path::PathBuf;

use occt_geom::{Curve, Surface};
use occt_geom2d::Curve2d;
use occt_topo::abs::ShapeType;
use occt_topo::brep_gprop_full::surface_properties;
use occt_topo::brep_tool::BRepTool;
use occt_topo::shape::{Face, TopoShape};
use occt_topo::step::read_step_file;
use occt_topo::tgeometry::GeometryRegistry;
use occt_topo::topo_tools_full::{edges_of_wire, wires_of_face};

fn collect_faces(shape: &TopoShape, out: &mut Vec<TopoShape>) {
    let kids: Vec<TopoShape> = shape.tshape.read().unwrap().children.clone();
    for k in kids {
        match k.shape_type() {
            ShapeType::Face => out.push(k),
            ShapeType::Shell | ShapeType::Solid | ShapeType::CompSolid | ShapeType::Compound => {
                collect_faces(&k, out)
            }
            _ => {}
        }
    }
}

fn surf_name(s: &dyn Surface) -> &'static str {
    if s.gp_pln().is_some() {
        "Geom_Plane"
    } else if s.gp_cylinder().is_some() {
        "Geom_CylindricalSurface"
    } else if s.gp_cone().is_some() {
        "Geom_ConicalSurface"
    } else if s.gp_sphere().is_some() {
        "Geom_SphericalSurface"
    } else if s.gp_torus().is_some() {
        "Geom_ToroidalSurface"
    } else {
        "other"
    }
}

fn c2_name(c: &dyn Curve2d) -> &'static str {
    if c.trimmed_basis().is_some() {
        "Geom2d_TrimmedCurve"
    } else if c.gp_lin2d().is_some() {
        "Geom2d_Line"
    } else if c.gp_circ2d().is_some() {
        "Geom2d_Circle"
    } else if c.gp_elips2d().is_some() {
        "Geom2d_Ellipse"
    } else if c.bspline_poles2d().is_some() {
        "Geom2d_BSplineCurve"
    } else {
        "?"
    }
}

fn dump_face(face: &Face, index: usize) {
    let surf = BRepTool::face_surface(face).expect("surface");
    let fwd = Face(face.0.oriented(occt_topo::abs::Orientation::Forward));
    let face_key = GeometryRegistry::shape_key(&face.0);
    let reg = GeometryRegistry::global();
    println!("== dump FACE {index} {} ==", surf_name(surf.as_ref()));
    // Re-run this one face with the per-arc probe switched on so the arc dump
    // stays scoped to the requested face.
    std::env::set_var("PROBE_ARCS", "1");
    std::env::set_var("PROBE70", "1");
    let (_, own_area) = surface_properties(&face.0).expect("face area");
    std::env::remove_var("PROBE_ARCS");
    std::env::remove_var("PROBE70");
    println!(" own area = {own_area:.9}");
    for (wi, w) in wires_of_face(&fwd).iter().enumerate() {
        println!(
            " WIRE[{wi}] ori={:?}",
            w.0.orientation()
        );
        for e in edges_of_wire(w) {
            let pcs = reg.edge_pcurves(&e.0, face_key);
            let range = reg.pcurve_range(&e.0, face_key);
            let c3 = BRepTool::edge_curve(&e);
            print!(
                "  EDGE ori={:?} npc={} rng={:?}",
                e.0.orientation(),
                pcs.len(),
                range
            );
            if let Some(c) = pcs.first() {
                let (a, b) = range.unwrap_or((c.first_parameter(), c.last_parameter()));
                let p0 = c.d0(a);
                let p1 = c.d0(b);
                let pm = c.d0(0.5 * (a + b));
                print!(
                    " {} ({a:.6},{b:.6}) uv0=({:.6},{:.6}) uv1=({:.6},{:.6}) mid=({:.6},{:.6})",
                    c2_name(c.as_ref()),
                    p0.x(),
                    p0.y(),
                    p1.x(),
                    p1.y(),
                    pm.x(),
                    pm.y()
                );
            }
            if let Some(c) = c3 {
                let (a, b) = BRepTool::edge_parameters(&e);
                let p0 = c.d0(a);
                let p1 = c.d0(b);
                print!(
                    " c3 ({a:.6},{b:.6}) p0=({:.6},{:.6},{:.6}) p1=({:.6},{:.6},{:.6})",
                    p0.x(),
                    p0.y(),
                    p0.z(),
                    p1.x(),
                    p1.y(),
                    p1.z()
                );
            }
            println!();
        }
    }
}

fn main() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("..").join("..");
    let args: Vec<String> = std::env::args().skip(1).collect();
    let name = args.first().map_or("occ/T0M.stp", |s| s.as_str()).to_string();
    let min_area: f64 = args.get(1).map_or(-1e300, |s| s.parse().unwrap());
    let max_area: f64 = args.get(2).map_or(1e300, |s| s.parse().unwrap());
    let path = if name.contains('/') || name.contains('\\') {
        root.join(&name)
    } else {
        root.join("data").join(format!("{name}.step"))
    };
    let model = read_step_file(&path.to_string_lossy()).expect("read step");
    let shape = model.shapes[0].shape.clone();

    let (_, total) = surface_properties(&shape).expect("surface properties");
    println!("TOTAL area = {total:.9}");

    let mut faces: Vec<TopoShape> = Vec::new();
    if shape.shape_type() == ShapeType::Face {
        faces.push(shape.clone());
    } else {
        collect_faces(&shape, &mut faces);
    }
    println!("faces = {}", faces.len());
    let dump_index: Option<usize> = std::env::var("ZZFACE").ok().and_then(|v| v.parse().ok());
    if std::env::var("ZZNOPC").is_ok() {
        let reg = GeometryRegistry::global();
        for (fi, f) in faces.iter().enumerate() {
            let face = Face(f.clone());
            let fwd = Face(f.clone().oriented(occt_topo::abs::Orientation::Forward));
            let key = GeometryRegistry::shape_key(&f);
            let ws = wires_of_face(&fwd);
            let mut n_edges = 0usize;
            let mut n_pc = 0usize;
            for w in &ws {
                for e in edges_of_wire(w) {
                    n_edges += 1;
                    n_pc += reg.edge_pcurves(&e.0, key).len();
                }
            }
            if n_pc == 0 {
                let kind = BRepTool::face_surface(&face)
                    .map(|s| surf_name(s.as_ref()))
                    .unwrap_or("null");
                println!("NOPC FACE {} type={} wires={} edges={n_edges}", fi + 1, kind, ws.len());
            }
        }
    }
    let mut sum = 0.0;
    for (fi, f) in faces.iter().enumerate() {
        let face = Face(f.clone());
        let Ok((_, area)) = surface_properties(&f.clone()) else {
            continue;
        };
        sum += area;
        if dump_index == Some(fi + 1) {
            dump_face(&face, fi + 1);
        }
        if area < min_area || area > max_area {
            continue;
        }
        let kind = BRepTool::face_surface(&face)
            .map(|s| surf_name(s.as_ref()))
            .unwrap_or("null");
        let (bu1, bu2, bv1, bv2) = occt_topo::brep_tools::uv_bounds(&face);
        println!(
            "FACE {} area={:.9} type={} uv=({bu1:.12},{bu2:.12},{bv1:.12},{bv2:.12})",
            fi + 1,
            area,
            kind
        );
    }
    println!("SUM of per-face areas = {sum:.9}");
}
