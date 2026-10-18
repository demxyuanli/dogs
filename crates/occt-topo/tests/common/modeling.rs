//! BRep modeling parity oracle: the frozen OCCT construction records under
//! `data/modeling/`, produced by `specs/modeling_probe` (see its header for the
//! generation recipe).
//!
//! Each `<name>.txt` carries the topology census (`TopExp::MapShapes` per shape
//! type), `BRepCheck_Analyzer` validity, `BRepGProp` volume / surface area /
//! centroid and the `Bnd_Box` of the shape built by the OCCT construction API;
//! each `<name>.obj` is that shape tessellated by `BRepMesh_IncrementalMesh` at
//! absolute deflection 0.1. `tests/modeling_gates.rs` rebuilds the same model
//! with the port's construction API and compares.
#![allow(dead_code, unused_imports)]

use occt_core::gp::{GpPnt, GpVec};
use occt_topo::abs::ShapeType;
use occt_topo::brep_builder_api::{make_face_from_polygon, make_prism};
use occt_topo::primitives::{
    BRepPrimBox, BRepPrimCone, BRepPrimCylinder, BRepPrimSphere, BRepPrimTorus,
};
use occt_topo::shape::TopoShape;
use occt_topo::topo_tools::shape_counts;

use super::data_dir;

/// Distinct-subshape census matching `TopExp::MapShapes` per type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Counts {
    pub v: usize,
    pub e: usize,
    pub w: usize,
    pub f: usize,
    pub sh: usize,
    pub so: usize,
}

/// One frozen OCCT construction record.
#[derive(Debug, Clone)]
pub struct Oracle {
    pub counts: Counts,
    pub valid: bool,
    pub volume: f64,
    pub area: f64,
    pub centroid: [f64; 3],
    pub bbox: [f64; 6],
}

/// Path of `data/modeling/`.
pub fn modeling_dir() -> std::path::PathBuf {
    data_dir().join("modeling")
}

/// Parse `data/modeling/<name>.txt`.
pub fn oracle(name: &str) -> Oracle {
    let path = modeling_dir().join(format!("{name}.txt"));
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("missing oracle {}: {e} (run specs/modeling_probe)", path.display()));
    let mut counts = None;
    let mut valid = None;
    let mut volume = None;
    let mut area = None;
    let mut centroid = None;
    let mut bbox = None;
    for line in text.lines() {
        let (key, rest) = match line.split_once(' ') {
            Some(kv) => kv,
            None => continue,
        };
        match key {
            "COUNTS" => {
                let map: std::collections::HashMap<&str, usize> = rest
                    .split_whitespace()
                    .filter_map(|t| t.split_once('='))
                    .map(|(k, v)| (k, v.parse().unwrap_or_else(|_| panic!("{name}: bad count {k}={v}"))))
                    .collect();
                counts = Some(Counts {
                    v: map["V"],
                    e: map["E"],
                    w: map["W"],
                    f: map["FA"],
                    sh: map["SH"],
                    so: map["SO"],
                });
            }
            "VALID" => valid = Some(rest.trim() == "1"),
            "VOLUME" => volume = Some(rest.trim().parse().unwrap_or_else(|_| panic!("{name}: {line}"))),
            "AREA" => area = Some(rest.trim().parse().unwrap_or_else(|_| panic!("{name}: {line}"))),
            "CENTROID" => {
                let v: Vec<f64> = rest.split_whitespace().map(|t| t.parse().unwrap()).collect();
                centroid = Some([v[0], v[1], v[2]]);
            }
            "BBOX" => {
                let v: Vec<f64> = rest.split_whitespace().map(|t| t.parse().unwrap()).collect();
                bbox = Some([v[0], v[1], v[2], v[3], v[4], v[5]]);
            }
            _ => {}
        }
    }
    Oracle {
        counts: counts.unwrap_or_else(|| panic!("{name}.txt: no COUNTS line")),
        valid: valid.unwrap_or_else(|| panic!("{name}.txt: no VALID line")),
        volume: volume.unwrap_or_else(|| panic!("{name}.txt: no VOLUME line")),
        area: area.unwrap_or_else(|| panic!("{name}.txt: no AREA line")),
        centroid: centroid.unwrap_or_else(|| panic!("{name}.txt: no CENTROID line")),
        bbox: bbox.unwrap_or_else(|| panic!("{name}.txt: no BBOX line")),
    }
}

/// The OCCT tessellation at deflection 0.1.
pub fn oracle_obj(name: &str) -> String {
    let path = modeling_dir().join(format!("{name}.obj"));
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("missing {}: {e}", path.display()))
}

/// Port-side census, same shape-keyed uniqueness as `TopExp::MapShapes`.
pub fn rust_counts(shape: &TopoShape) -> Counts {
    let m = shape_counts(shape);
    let n = |t: ShapeType| m.get(&t).copied().unwrap_or(0);
    Counts {
        v: n(ShapeType::Vertex),
        e: n(ShapeType::Edge),
        w: n(ShapeType::Wire),
        f: n(ShapeType::Face),
        sh: n(ShapeType::Shell),
        so: n(ShapeType::Solid),
    }
}

/// G0 primitive set: the six constructions the port already exposes.
/// Parameters equal `data/*.step` so the modeled solid also matches the
/// STEP-read solid (Cube 10^3, Cylinder r2 h10, Cone r4 h10, Sphere r5,
/// Torus R10 r2; the prism is a 2x3 rectangle extruded +Z by 4).
pub const MODELS: &[&str] = &[
    "prim-box",
    "prim-cylinder",
    "prim-cone",
    "prim-sphere",
    "prim-torus",
    "prim-prism",
];

/// Rebuild one modeling-gate model with the port's construction API.
pub fn build(name: &str) -> TopoShape {
    match name {
        "prim-box" => BRepPrimBox::make_box(10.0, 10.0, 10.0).solid.0,
        "prim-cylinder" => BRepPrimCylinder::make_cylinder(2.0, 10.0).solid.0,
        "prim-cone" => BRepPrimCone::make_cone(4.0, 10.0).solid.0,
        "prim-sphere" => BRepPrimSphere::make_sphere(5.0).solid.0,
        "prim-torus" => BRepPrimTorus::make_torus(10.0, 2.0).solid.0,
        "prim-prism" => {
            let profile = [
                GpPnt::new(0.0, 0.0, 0.0),
                GpPnt::new(2.0, 0.0, 0.0),
                GpPnt::new(2.0, 3.0, 0.0),
                GpPnt::new(0.0, 3.0, 0.0),
            ];
            let face = make_face_from_polygon(&profile).expect("prim-prism profile face");
            // The oracle builds this with `BRepPrimAPI_MakePrism`
            // (`specs/modeling_probe/modeling_probe.cpp`), so the gate calls the
            // port of that class: `brep_sweep::make_prism`
            // (`BRepPrimAPI_MakePrism.cxx:36-42` -> `BRepSweep_Prism`).
            occt_topo::brep_sweep::sweep_make_prism(&face.0, &GpVec::new(0.0, 0.0, 4.0))
                .expect("prim-prism sweep")
        }
        other => panic!("no port construction registered for modeling model {other}"),
    }
}
