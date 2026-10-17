//! BRep modeling parity gates: the port's construction API vs the frozen OCCT
//! ground truth under `data/modeling/` (see `specs/modeling_probe`).
//!
//! Three layers, in increasing strictness of what they prove:
//! 1. `modeling_topology_matches_occt` - the shape census (`TopExp::MapShapes`
//!    per type) and `BRepCheck_Analyzer` validity. Exact equality; a differing
//!    face/edge/wire count is a topology gap, not a tolerance.
//! 2. `modeling_gprop_matches_occt` - `BRepGProp` volume / surface area /
//!    centroid against OCCT on the same nominal construction. Every G0 model
//!    also has a closed-form value (box 1000, cylinder 40*pi, cone 160*pi/3,
//!    sphere 500*pi/3, torus 80*pi^2, prism 24), so this layer does not depend
//!    on OCCT's own numerical choices.
//! 3. `modeling_mesh_matches_occt` - `BRepMesh_IncrementalMesh` at deflection
//!    0.1 exported to OBJ, compared with `data/modeling/<model>.obj`.
//!
//! Each layer iterates the whole model table and reports *every* mismatching
//! model, so one run carries the complete ladder status.
//!
//! Tolerances encode a measured gap; tightening one is a port improvement, not
//! a test change. The generation recipe and the current G0 reading (with its
//! open diffs) are in the `modeling_probe.cpp` header.

mod common;

use common::modeling::{build, oracle, oracle_obj, rust_counts, MODELS};
use common::{obj_area, obj_bbox};
use occt_core::gp::{GpPnt, GpVec};
use occt_core::io::obj::ObjMesh;
use occt_topo::brep_builder_api::{make_face_from_polygon, make_prism};
use occt_topo::brep_check::Analyzer;
use occt_topo::brep_exchange::brep_to_obj;
use occt_topo::brep_gprop_full::{surface_properties, volume_properties};

/// Volume / area relative tolerance. OCCT and the port both integrate the same
/// analytic surfaces, so this is a double-precision level gate.
const GPROP_REL: f64 = 1e-9;
/// Centroid absolute tolerance.
const CENTROID_ABS: f64 = 1e-7;
/// Planar models: the mesh equals the boundary exactly.
const BBOX_EXACT: f64 = 1e-6;
/// Curved models: the two meshers pick different chord points on the surface.
const BBOX_SAMPLING: f64 = 0.05;
/// Planar mesh area.
const AREA_EXACT: f64 = 1e-3;
/// Curved mesh area (chord approximation of the same deflection).
const AREA_CURVED: f64 = 0.05;

/// Per-model mesh tolerances (same split as the STEP -> OBJ gate).
fn mesh_tols(name: &str) -> (f64, f64) {
    match name {
        "prim-box" | "prim-prism" => (BBOX_EXACT, AREA_EXACT),
        _ => (BBOX_SAMPLING, AREA_CURVED),
    }
}

/// Layer 1: topology census + validity, exact.
#[test]
fn modeling_topology_matches_occt() {
    let mut fails = Vec::new();
    for name in MODELS {
        let o = oracle(name);
        let shape = build(name);
        let c = rust_counts(&shape);
        if c != o.counts {
            fails.push(format!(
                "{name}: census ours {c:?} occt {:?} (v{} e{} w{} f{} sh{} so{})",
                o.counts,
                c.v as i64 - o.counts.v as i64,
                c.e as i64 - o.counts.e as i64,
                c.w as i64 - o.counts.w as i64,
                c.f as i64 - o.counts.f as i64,
                c.sh as i64 - o.counts.sh as i64,
                c.so as i64 - o.counts.so as i64,
            ));
        }
        let analyzer = Analyzer::new(shape.clone());
        let valid = analyzer.is_valid();
        eprintln!(
            "[{name}] census ours {c:?} occt {:?} | valid ours {valid} occt {}",
            o.counts, o.valid
        );
        if valid != o.valid {
            fails.push(format!("{name}: BRepCheck validity ours {valid} occt {}", o.valid));
            for r in analyzer.results() {
                if !r.is_valid() {
                    eprintln!(
                        "  [{name}] {:?} {:?}: {:?}",
                        r.shape.shape_type(),
                        r.shape.tshape,
                        r.statuses
                    );
                }
            }
        }
    }
    assert!(fails.is_empty(), "topology gate:\n{}", fails.join("\n"));
}

/// Layer 2: mass properties, volume / area / centroid.
#[test]
fn modeling_gprop_matches_occt() {
    let mut fails = Vec::new();
    for name in MODELS {
        let o = oracle(name);
        let shape = build(name);
        let vol = match volume_properties(&shape) {
            Ok(v) => v,
            Err(e) => {
                fails.push(format!("{name}: volume_properties: {e}"));
                continue;
            }
        };
        let surf = match surface_properties(&shape) {
            Ok((s, _)) => s,
            Err(e) => {
                fails.push(format!("{name}: surface_properties: {e}"));
                continue;
            }
        };
        let (rv, ov) = (vol.mass(), o.volume);
        let (ra, oa) = (surf.mass(), o.area);
        eprintln!("[{name}] vol ours={rv:.9} occt={ov:.9} | area ours={ra:.9} occt={oa:.9}");
        for f in occt_topo::topo_tools_full::faces_of(&shape) {
            eprintln!(
                "  [{name}] face contrib = {:?}",
                occt_topo::brep_gprop_full::face_volume_contribution(&f)
            );
        }
        if *name == "prim-prism" {
            for (dx, dy, dz) in [(0.0, 0.0, 0.0), (2.0, 3.0, 0.0), (0.0, 0.0, 5.0), (2.0, 3.0, 5.0)] {
                let profile = [
                    GpPnt::new(dx, dy, dz),
                    GpPnt::new(2.0 + dx, dy, dz),
                    GpPnt::new(2.0 + dx, 3.0 + dy, dz),
                    GpPnt::new(dx, 3.0 + dy, dz),
                ];
                let f2 = make_face_from_polygon(&profile).expect("shifted face");
                let s2 = make_prism(&f2, &GpVec::new(0.0, 0.0, 1.0), 4.0)
                    .expect("shifted prism")
                    .solid
                    .0;
                let v2 = volume_properties(&s2).expect("shifted volume");
                let per: Vec<f64> = occt_topo::topo_tools_full::faces_of(&s2)
                    .iter()
                    .map(|f| occt_topo::brep_gprop_full::face_volume_contribution(f).unwrap_or(f64::NAN))
                    .collect();
                eprintln!(
                    "  [prism {:?}] vol = {:?} per-face = {:?} sum = {}",
                    (dx, dy, dz),
                    v2.mass(),
                    per,
                    per.iter().sum::<f64>()
                );
            }
        }
        if (rv - ov).abs() > GPROP_REL * ov.abs().max(1.0) {
            fails.push(format!("{name}: volume ours {rv:.9} occt {ov:.9}"));
        }
        if (ra - oa).abs() > GPROP_REL * oa.abs().max(1.0) {
            fails.push(format!("{name}: area ours {ra:.9} occt {oa:.9}"));
        }
        let c = vol.center();
        let ours = [c.x(), c.y(), c.z()];
        for (i, axis) in ["x", "y", "z"].iter().enumerate() {
            if (ours[i] - o.centroid[i]).abs() > CENTROID_ABS {
                fails.push(format!(
                    "{name}: centroid {axis} ours {:.9} occt {:.9}",
                    ours[i], o.centroid[i]
                ));
            }
        }
    }
    assert!(fails.is_empty(), "gprop gate:\n{}", fails.join("\n"));
}

/// Layer 3: the deflection-0.1 tessellation (area ratio and bounding box).
#[test]
fn modeling_mesh_matches_occt() {
    let mut fails = Vec::new();
    for name in MODELS {
        let (bbox_tol, area_tol) = mesh_tols(name);
        let shape = build(name);
        let ours = brep_to_obj(&shape, 0.1);
        if std::env::var("DUMP_OBJ").is_ok() {
            std::fs::write(format!("target/{name}-ours.obj"), &ours).unwrap();
            std::fs::write(format!("target/{name}-occt.obj"), oracle_obj(name)).unwrap();
        }
        let reference = oracle_obj(name);
        let our_mesh = match ObjMesh::parse(&ours) {
            Ok(m) => m,
            Err(e) => {
                fails.push(format!("{name}: our OBJ does not parse: {e} ({} bytes)", ours.len()));
                continue;
            }
        };
        let occ_mesh = ObjMesh::parse(&reference)
            .unwrap_or_else(|e| panic!("{name}: reference OBJ does not parse: {e}"));
        let our_area = obj_area(&our_mesh);
        let occ_area = obj_area(&occ_mesh);
        let ratio = our_area / occ_area;
        let (Some((our_min, our_max)), Some((occ_min, occ_max))) =
            (obj_bbox(&ours), obj_bbox(&reference))
        else {
            fails.push(format!(
                "{name}: OBJ without vertices (ours {} verts, occt {} verts)",
                our_mesh.vertices.len(),
                occ_mesh.vertices.len()
            ));
            continue;
        };
        eprintln!("[{name}] mesh area ours={our_area:.6} occt={occ_area:.6} ratio={ratio:.6}");
        if (ratio - 1.0).abs() >= area_tol {
            fails.push(format!(
                "{name}: mesh area ratio {ratio:.6} outside {area_tol} (ours {our_area:.6}, occt {occ_area:.6})"
            ));
        }
        for i in 0..3 {
            if (our_min[i] - occ_min[i]).abs() >= bbox_tol {
                fails.push(format!(
                    "{name}: bbox min[{i}] ours {} occt {} (tol {bbox_tol})",
                    our_min[i], occ_min[i]
                ));
            }
            if (our_max[i] - occ_max[i]).abs() >= bbox_tol {
                fails.push(format!(
                    "{name}: bbox max[{i}] ours {} occt {} (tol {bbox_tol})",
                    our_max[i], occ_max[i]
                ));
            }
        }
    }
    assert!(fails.is_empty(), "mesh gate:\n{}", fails.join("\n"));
}
