//! Deep geometry validation for the STEP samples without an OCCT reference
//! OBJ (`Extrusion`, `Offset`).
//!
//! For these two, there is no `occ-*.obj` to compare against, so the validation
//! is geometric self-consistency:
//!
//!  1. **Surface area** against an analytic decomposition of the shape.
//!  2. **Mesh volume** via the divergence theorem (converges as the deflection
//!     refines) and via `shape_volume` (mesh with outward normals).
//!  3. **Topology**: the exported solid is watertight, a single closed shell
//!     with Euler characteristic 2 (a closed orientable genus-0 surface).
//!
//! `Extrusion` is a plain 20×5×15 box (6 planar faces) — area 950 and volume
//! 1500 are exact. `Offset` is a 10×10×10 cube with outward-rounded edges
//! (radius 2, twelve quarter-cylinders) and spherical corner bands: analytic
//! area 600 + 12·(π/4·2²·10) + 8·(spherical u90°×v90° band) ≈ 1047.
use std::path::PathBuf;

use occt_core::io::obj::{ObjMesh, read_obj_file};
use occt_topo::brep_exchange::brep_to_obj;
use occt_topo::shape_mesh::shape_volume;
use occt_topo::step::read_step_file;

fn data_dir() -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.push("..");
    p.push("..");
    p.push("data");
    p
}

/// Fan-triangle area of an OBJ mesh.
fn obj_area(obj: &ObjMesh) -> f64 {
    let mut area = 0.0;
    for f in &obj.faces {
        let vs: Vec<usize> = f.v.iter().map(|&i| i as usize).collect();
        if vs.len() < 3 { continue; }
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

/// Divergence-theorem volume from the signed triangles.
fn obj_volume(obj: &ObjMesh) -> f64 {
    let mut s = 0.0;
    for f in &obj.faces {
        let vs: Vec<usize> = f.v.iter().map(|&i| i as usize).collect();
        if vs.len() < 3 { continue; }
        let a = &obj.vertices[vs[0]];
        for k in 1..vs.len() - 1 {
            let b = &obj.vertices[vs[k]];
            let c = &obj.vertices[vs[k + 1]];
            let ab = b.coord.subtracted(&a.coord);
            let ac = c.coord.subtracted(&a.coord);
            s += a.coord.dot(&ab.crossed(&ac));
        }
    }
    s.abs() / 6.0
}

/// Assert the shape is a closed watertight solid with Euler characteristic 2.
fn assert_closed_solid(step: &str, sh: &occt_topo::shape::TopoShape) {
    assert!(occt_topo::shape_analysis::is_watertight(sh), "{step}: not watertight");
    let solid = occt_topo::shape::Solid(sh.clone());
    let shells = occt_topo::topo_tools_full::shapes_of(&solid.0, occt_topo::abs::ShapeType::Shell);
    assert_eq!(shells.len(), 1, "{step}: expected 1 shell, got {}", shells.len());
    for shl in &shells {
        let s2 = occt_topo::shape::Shell(shl.clone());
        assert!(occt_topo::shell_check::shell_is_closed(&s2), "{step}: shell open");
        let f = occt_topo::topo_tools_full::faces_of(&s2.0).len();
        let e = occt_topo::topo_tools_full::edges_of(&s2.0).len();
        let v = occt_topo::topo_tools_full::vertices_of(&s2.0).len();
        let chi = v as i64 - e as i64 + f as i64;
        assert_eq!(chi, 2, "{step}: Euler characteristic {chi}, expected 2");
    }
}

#[test]
fn extrusion_geometry_is_exact() {
    // Extrusion.step is a plain 20×5×15 box: exact area 950, volume 1500.
    let model = read_step_file(&data_dir().join("Extrusion.step").to_string_lossy())
        .unwrap_or_else(|e| panic!("read Extrusion.step: {e}"));
    let sh = &model.shapes[0].shape;
    let obj = ObjMesh::parse(&brep_to_obj(sh, 0.1)).expect("parse OBJ");
    let area = obj_area(&obj);
    let vol = obj_volume(&obj);
    let svol = shape_volume(sh, 0.1);
    assert!((area - 950.0).abs() < 1e-6, "Extrusion area {area}, expected 950");
    assert!((vol - 1500.0).abs() < 1e-6, "Extrusion volume {vol}, expected 1500");
    assert!((svol - 1500.0).abs() < 1e-6, "Extrusion shape_volume {svol}, expected 1500");
    assert_closed_solid("Extrusion", sh);
}

#[test]
fn offset_geometry_is_consistent() {
    // Offset.step is a 10×10×10 cube with outward-rounded edges (r=2) and
    // spherical corner bands. Analytic area:
    //   6 flat faces 10×10 = 600
    //   + 12 edge quarter-cylinders r=2 × length 10 = 12·(π/4·4·10) = 120π ≈ 377
    //   + 8 spherical corner bands (u90°×v90°) ≈ 8·8.88 ≈ 71
    //   = 1047. The mesh area must agree within the curved-face sampling tol;
    //   the volume must converge as the deflection refines.
    let model = read_step_file(&data_dir().join("Offset.step").to_string_lossy())
        .unwrap_or_else(|e| panic!("read Offset.step: {e}"));
    let sh = &model.shapes[0].shape;
    let obj = ObjMesh::parse(&brep_to_obj(sh, 0.02)).expect("parse OBJ");
    let area = obj_area(&obj);
    let analytic = 600.0 + 12.0 * std::f64::consts::PI * 4.0 * 10.0 / 4.0 + 8.0 * 8.88;
    assert!(
        (area - analytic).abs() / analytic < 0.03,
        "Offset area {area}, analytic ~{analytic:.1}"
    );

    // Volume: fine deflection must converge (both divergence and shape_volume).
    let v01 = obj_volume(&ObjMesh::parse(&brep_to_obj(sh, 0.01)).unwrap());
    let v005 = obj_volume(&ObjMesh::parse(&brep_to_obj(sh, 0.005)).unwrap());
    assert!(
        (v01 - v005).abs() / v005 < 0.01,
        "Offset volume not converging: 0.01→{v01:.1} 0.005→{v005:.1}"
    );
    let svol = shape_volume(sh, 0.01);
    assert!(
        (svol - v01).abs() / v01 < 0.02,
        "Offset divergence {v01:.1} vs shape_volume {svol:.1}"
    );
    assert_closed_solid("Offset", sh);
}

#[test]
fn offset_area_breaks_down_by_face_type() {
    // The 26 faces decompose into 6 planes, 12 quarter-cylinders, 8 spheres.
    // Mesh each face in isolation and check the per-type totals match the
    // analytic decomposition (600 / 376 / 71). The classifier recognises planes
    // and spheres; the twelve curved edges classify as `Other` (the classifier
    // has no cylinder branch yet), so they are counted as the remainder.
    let model = read_step_file(&data_dir().join("Offset.step").to_string_lossy())
        .unwrap_or_else(|e| panic!("read Offset.step: {e}"));
    let sh = &model.shapes[0].shape;
    let mut plane = 0.0f64;
    let mut sph = 0.0f64;
    let mut other = 0.0f64;
    let mut nplane = 0usize;
    let mut nsph = 0usize;
    let mut nother = 0usize;
    for f in occt_topo::topo_tools_full::faces_of(sh) {
        let surf = occt_topo::brep_tool::BRepTool::face_surface(&f).unwrap();
        let kind = occt_topo::brep_surface::classify_surface(surf.as_ref());
        let obj = ObjMesh::parse(&brep_to_obj(&f.0, 0.02)).unwrap();
        let a = obj_area(&obj);
        match kind {
            occt_topo::brep_surface::SurfaceKind::Plane => {
                plane += a;
                nplane += 1;
            }
            occt_topo::brep_surface::SurfaceKind::Sphere => {
                sph += a;
                nsph += 1;
            }
            _ => {
                other += a;
                nother += 1;
            }
        }
    }
    assert_eq!(nplane, 6, "expected 6 planar faces");
    assert_eq!(nsph, 8, "expected 8 spherical faces");
    assert_eq!(nother, 12, "expected 12 curved (cylinder) edge faces");
    assert!((plane - 600.0).abs() < 1.0, "planes area {plane}, expected 600");
    let cyl_analytic = 12.0 * std::f64::consts::PI * 4.0 * 10.0 / 4.0;
    assert!(
        (other - cyl_analytic).abs() / cyl_analytic < 0.03,
        "cylinder area {other}, expected {cyl_analytic:.1}"
    );
    assert!(sph > 40.0 && sph < 100.0, "sphere area {sph}, expected ~71");
}
