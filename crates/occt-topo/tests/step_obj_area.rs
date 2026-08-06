//! STEP → OBJ **area** parity against the OCCT references in `data/`.
//!
//! The bbox parity test (`step_obj_parity.rs`) gates geometry *extent*; this
//! test gates the *surface area* of the exported mesh against the OCCT
//! reference OBJ. Both meshes are triangle approximations of the same analytic
//! shape, so the areas agree within the meshing tolerance: planar faces are
//! exact, curved faces differ by the chord/UV-grid sampling (< 4% at the 0.1
//! deflection the tests use).
//!
//! Pipeline: `step::read_step_file` → `brep_exchange::brep_to_obj` (0.1), both
//! areas measured with the same fan-triangle rule on [`ObjMesh::parse`].
use std::path::PathBuf;

use occt_core::io::obj::{ObjMesh, read_obj_file};
use occt_topo::brep_exchange::brep_to_obj;
use occt_topo::step::read_step_file;

/// `data/` directory (tests run with CWD = crate root of occt-topo).
fn data_dir() -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.push("..");
    p.push("..");
    p.push("data");
    p
}

/// Total triangle area of an ObjMesh (fan rule around the first vertex).
fn obj_area(obj: &ObjMesh) -> f64 {
    let mut area = 0.0;
    for f in &obj.faces {
        let vs: Vec<usize> = f.v.iter().map(|&i| i as usize).collect();
        if vs.len() < 3 {
            continue;
        }
        let a = &obj.vertices[vs[0]];
        for k in 1..vs.len() - 1 {
            let b = &obj.vertices[vs[k]];
            let c = &obj.vertices[vs[k + 1]];
            let ab = b.coord.subtracted(&a.coord);
            let ac = c.coord.subtracted(&a.coord);
            area += 0.5 * ab.crossed(&ac).modulus();
        }
    }
    area
}

/// Export `step` to OBJ, read the OCCT reference, and assert the area ratio is
/// within `tol` of 1. Returns the measured ratio (for the eprintln log).
fn check_area_parity(step: &str, occ: &str, tol: f64) -> f64 {
    let model = read_step_file(&data_dir().join(format!("{step}.step")).to_string_lossy())
        .unwrap_or_else(|e| panic!("read {step}.step: {e}"));
    assert!(!model.shapes.is_empty(), "{step}: parsed to no shapes");
    let our = obj_area(&ObjMesh::parse(&brep_to_obj(&model.shapes[0].shape, 0.1)).expect("parse our OBJ"));
    let occ_mesh =
        read_obj_file(&data_dir().join(occ).to_string_lossy()).unwrap_or_else(|e| panic!("read {occ}: {e}"));
    let occ_area = obj_area(&occ_mesh);
    let ratio = our / occ_area;
    eprintln!("[{step}] our={our:.2} occ={occ_area:.2} ratio={ratio:.4}");
    assert!(
        (ratio - 1.0).abs() < tol,
        "{step}: area ratio {ratio:.4} outside tolerance {tol} (ours {our:.2}, occ {occ_area:.2})"
    );
    ratio
}

// Planar / analytic shapes whose mesh areas agree essentially exactly.
const EXACT_AREA_TOL: f64 = 1e-3;
// Curved shapes: the OCCT reference and our UV grid sample the surface at
// different points; a chord/deflection difference of a few percent is expected.
const CURVED_AREA_TOL: f64 = 0.04;
// Complex B-spline / offset shapes: same sampling tolerance, slightly wider.
const COMPLEX_AREA_TOL: f64 = 0.05;

#[test]
fn cube_area_matches_occt() {
    check_area_parity("Cube", "occ-cube.obj", EXACT_AREA_TOL);
}

#[test]
fn offset_plane_hole_edge_area_matches_occt() {
    check_area_parity("OffsetPlaneHoleEdge", "occ-OffsetPlaneHoleEdge.obj", EXACT_AREA_TOL);
}

#[test]
fn cone_area_matches_occt() {
    check_area_parity("Cone", "occ-cone.obj", CURVED_AREA_TOL);
}

#[test]
fn cylinder_area_matches_occt() {
    check_area_parity("Cylinder", "occ-cylinder.obj", CURVED_AREA_TOL);
}

#[test]
fn torus_area_matches_occt() {
    check_area_parity("Torus", "occ-torus.obj", CURVED_AREA_TOL);
}

#[test]
fn sphere_area_matches_occt() {
    check_area_parity("Sphere", "occ-shpere.obj", CURVED_AREA_TOL);
}

#[test]
fn rev_area_matches_occt() {
    check_area_parity("rev", "occ-rev.obj", CURVED_AREA_TOL);
}

#[test]
fn holed_plate_area_matches_occt() {
    check_area_parity("HoledPlate", "occ-HoledPlate.obj", COMPLEX_AREA_TOL);
}

#[test]
fn shape2_area_matches_occt() {
    check_area_parity("Shape-2", "occ-shape-2.obj", COMPLEX_AREA_TOL);
}

#[test]
fn shape1_area_matches_occt() {
    // Shape-1.step: the analytic surfaces (cylinders/spheres/tori) are the
    // kettle's edge-trimming geometry — their face UV domains are bounded by the
    // pcurves, so the mesh covers only the trimmed patch. Bbox and area agree
    // with the OCCT reference once the pcurve clipping is honoured.
    check_area_parity("Shape-1", "occ-shape-1.obj", COMPLEX_AREA_TOL);
}

#[test]
fn shape_area_matches_occt() {
    check_area_parity("Shape", "occ-shape.obj", COMPLEX_AREA_TOL);
}
