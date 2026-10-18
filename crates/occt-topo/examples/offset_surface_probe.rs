//! TEMPORARY probe (T-126).
//!
//! Mirrors `specs/occt_probe/occt_probe.cpp --offsetsurf <dist>`: for every
//! face, print the basis surface type, the `BRepOffset::Surface`
//! (`BRepOffset.cxx:44-200`) result type, the `BRepOffset_Status` and the result
//! evaluated on a fixed `(u, v)` grid (independent of the face's UV window, so
//! this compares the surface rather than the port's `BRepTools::UVBounds`). The
//! face's trimmed window is unwrapped first, the way `BRepOffset_Offset::Init`
//! does (`BRepOffset_Offset.cxx:483-488`).
//!
//! Usage:
//!   cargo run --example offset_surface_probe -- data/Cone.step 0.5

use std::path::PathBuf;

use occt_geom::Surface;
use occt_topo::abs::ShapeType;
use occt_topo::brep_offset::brep_offset_surface;
use occt_topo::brep_tool::BRepTool;
use occt_topo::shape::{Face, TopoShape};
use occt_topo::step::read_step_file;

/// `DynamicType()->Name()` equivalent, for the types `BRepOffset::Surface`
/// dispatches on.
fn surf_kind(s: &dyn Surface) -> String {
    if s.rectangular_trimmed_basis().is_some() {
        "Geom_RectangularTrimmedSurface".into()
    } else if s.gp_pln().is_some() {
        "Geom_Plane".into()
    } else if s.gp_cylinder().is_some() {
        "Geom_CylindricalSurface".into()
    } else if s.gp_cone().is_some() {
        "Geom_ConicalSurface".into()
    } else if s.gp_sphere().is_some() {
        "Geom_SphericalSurface".into()
    } else if s.gp_torus().is_some() {
        "Geom_ToroidalSurface".into()
    } else if s.is_surface_of_revolution() {
        "Geom_SurfaceOfRevolution".into()
    } else if s.is_surface_of_linear_extrusion() {
        "Geom_SurfaceOfLinearExtrusion".into()
    } else if s.offset_basis_surface().is_some() {
        "Geom_OffsetSurface".into()
    } else if s.is_bspline_surface() {
        "Geom_BSplineSurface".into()
    } else if s.is_bezier_surface() {
        "Geom_BezierSurface".into()
    } else {
        "?".into()
    }
}

/// Faces of `shape` in `TopExp_Explorer(aShape, TopAbs_FACE)` order.
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

fn main() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("..").join("..");
    let args: Vec<String> = std::env::args().skip(1).collect();
    let name = args.first().expect("usage: offset_surface_probe <file.step> <dist>").clone();
    let dist: f64 = args.get(1).map_or(0.5, |s| s.parse().expect("distance"));
    let path = if name.contains('/') {
        root.join(&name)
    } else {
        root.join("data").join(format!("{name}.step"))
    };
    let model = read_step_file(&path.to_string_lossy()).expect("read step");
    let shape: TopoShape = model.shapes[0].shape.clone();
    let mut faces: Vec<TopoShape> = Vec::new();
    if shape.shape_type() == ShapeType::Face {
        faces.push(shape.clone());
    } else {
        collect_faces(&shape, &mut faces);
    }
    let us = [0.4, 1.1, 2.7];
    let vs = [-0.5, 0.3, 1.2];
    for (fi, f) in faces.iter().enumerate() {
        let face = Face(f.clone());
        let surf = BRepTool::face_surface(&face).expect("face surface");
        // `BRepOffset_Offset::Init` unwraps the trim before `BRepOffset::Surface`.
        let basis = surf
            .rectangular_trimmed_basis()
            .map_or_else(|| surf.clone(), |b| b);
        let (off, status) = brep_offset_surface(basis.as_ref(), dist, false).expect("offset");
        let mut line = format!(
            "OFFSURF f={fi} src={} res={} status={}",
            surf_kind(basis.as_ref()),
            surf_kind(off.as_ref()),
            status.code()
        );
        for i in 0..3 {
            let p = off.d0(us[i], vs[i]);
            line.push_str(&format!(" p{i}=({:.9},{:.9},{:.9})", p.x(), p.y(), p.z()));
        }
        println!("{line}");
    }
}
