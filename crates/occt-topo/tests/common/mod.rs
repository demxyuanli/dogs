//! Shared model table for the STEP -> OBJ gates.
//!
//! step_obj_parity (bounding box), step_obj_area (mesh surface area) and
//! step_to_obj (export validity) all iterate MODELS, so the three gates always
//! cover the same 23 STEP samples: the 15 data/*.step primitives plus the 8 OCCT
//! test models under data/occ/. Every entry has an OCCT WriteObj reference OBJ
//! produced from that same STEP file.
#![allow(dead_code)]

pub mod modeling;

use std::path::PathBuf;

use occt_core::io::obj::ObjMesh;
use occt_topo::brep_exchange::{document_shape, model_to_obj_caf};
use occt_topo::rwmesh::material_file_path;
use occt_topo::step::read_step_file;

/// data/ directory (tests run with CWD = the occt-topo crate root).
pub fn data_dir() -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.push("..");
    p.push("..");
    p.push("data");
    p
}

/// One STEP sample plus its OCCT reference OBJ and the gate tolerances.
pub struct Model {
    /// data/-relative STEP path.
    pub step: &'static str,
    /// data/-relative OCCT WriteObj reference OBJ.
    pub occ: &'static str,
    /// Bounding-box tolerance for step_obj_parity; None gates validity only.
    pub bbox_tol: Option<f64>,
    /// Mesh surface-area ratio tolerance for step_obj_area; None = not gated.
    pub area_tol: Option<f64>,
}

/// Faceted / planar shapes: the boundary vertices land exactly on the mesh.
const BBOX_EXACT: f64 = 1e-6;
/// Curved shapes: OCCT and the port sample the surface at different points, so
/// the bbox edge can differ by the chord sampling.
const BBOX_SAMPLING: f64 = 0.05;
/// Planar / analytic shapes whose mesh areas agree essentially exactly.
const AREA_EXACT: f64 = 1e-3;
/// Curved shapes: a few percent of chord/deflection difference is expected.
const AREA_CURVED: f64 = 0.04;
/// Complex B-spline / offset shapes: same sampling level, slightly wider.
const AREA_COMPLEX: f64 = 0.05;

/// The canonical model set, shared by all three STEP -> OBJ gates.
///
/// Tolerances encode the measured gap; tightening one is a port improvement,
/// not a test change. BBOX_SAMPLING / AREA_CURVED are the sampling-level gates;
/// the wide values are tracked STEP-import gaps recorded in
/// specs/_a3n00_gap_analysis.md.
pub const MODELS: &[Model] = &[
    // --- 15 data/*.step primitives / samples ---
    Model { step: "Cube.step", occ: "occ-cube.obj", bbox_tol: Some(BBOX_EXACT), area_tol: Some(AREA_EXACT) },
    Model { step: "Cone.step", occ: "occ-cone.obj", bbox_tol: Some(BBOX_EXACT), area_tol: Some(AREA_CURVED) },
    Model { step: "Cylinder.step", occ: "occ-cylinder.obj", bbox_tol: Some(BBOX_SAMPLING), area_tol: Some(AREA_CURVED) },
    Model { step: "Sphere.step", occ: "occ-sphere.obj", bbox_tol: Some(BBOX_SAMPLING), area_tol: Some(AREA_CURVED) },
    Model { step: "Torus.step", occ: "occ-torus.obj", bbox_tol: Some(BBOX_SAMPLING), area_tol: Some(AREA_CURVED) },
    Model { step: "HoledPlate.step", occ: "occ-HoledPlate.obj", bbox_tol: Some(BBOX_EXACT), area_tol: Some(AREA_COMPLEX) },
    Model { step: "OffsetPlaneHoleEdge.step", occ: "occ-OffsetPlaneHoleEdge.obj", bbox_tol: Some(BBOX_EXACT), area_tol: Some(AREA_EXACT) },
    Model { step: "rev.step", occ: "occ-rev.obj", bbox_tol: Some(BBOX_SAMPLING), area_tol: Some(AREA_CURVED) },
    Model { step: "Offset.step", occ: "occ-Offset.obj", bbox_tol: Some(BBOX_SAMPLING), area_tol: Some(AREA_COMPLEX) },
    // Shape / Shape-1 / Shape-2 are B-spline + offset-surface models: the bbox
    // tolerances are tracked STEP-import gaps (see the parity test notes).
    Model { step: "Shape.step", occ: "occ-shape.obj", bbox_tol: Some(0.1), area_tol: Some(AREA_COMPLEX) },
    Model { step: "Shape-1.step", occ: "occ-shape-1.obj", bbox_tol: Some(0.2), area_tol: Some(AREA_COMPLEX) },
    Model { step: "Shape-2.step", occ: "occ-shape-2.obj", bbox_tol: Some(BBOX_SAMPLING), area_tol: Some(AREA_COMPLEX) },
    Model { step: "Extrusion.step", occ: "occ-Extrusion.obj", bbox_tol: Some(BBOX_EXACT), area_tol: Some(AREA_EXACT) },
    Model { step: "linkrods.step", occ: "occ-linkrods.obj", bbox_tol: Some(BBOX_SAMPLING), area_tol: Some(AREA_CURVED) },
    Model { step: "screw.step", occ: "occ-screw.obj", bbox_tol: Some(BBOX_SAMPLING), area_tol: Some(AREA_CURVED) },
    // --- 8 OCCT test models under data/occ/ ---
    Model { step: "occ/ATU01038.step", occ: "occ/occ-ATU01038.obj", bbox_tol: Some(1e-4), area_tol: Some(0.01) },
    Model { step: "occ/bottom.step", occ: "occ/occ-bottom.obj", bbox_tol: Some(1e-3), area_tol: Some(0.05) },
    Model { step: "occ/motoc.step", occ: "occ/occ-motoc.obj", bbox_tol: Some(1e-3), area_tol: Some(0.05) },
    Model { step: "occ/top.step", occ: "occ/occ-top.obj", bbox_tol: Some(1e-3), area_tol: Some(0.01) },
    // T0M: tracked gap - the port still misses some face classes (density < 1);
    // the bbox/area tolerances encode that measured mismatch. The area tolerance
    // was widened 0.01 -> 0.025 by T-93 (a): the reader-side `FixMissingSeam`
    // bolt-on was removed because OCCT does not run `ShapeFix_Face::Perform` on
    // this path (specs/_a3n00_gap_analysis.md 9.219, proven by the link map), so
    // T0M's ratio moved 0.9995 -> 0.9786 (ours 188536.30 / occ 192658.64) while
    // a3n00's F113 started meshing (mv 0 -> 223). Re-measured at 9.625 (after
    // `ShapeFix_Face::FixOrientation`): 0.99998 (192654.62 / 192658.64), 67524
    // vs 67366 triangles - the tolerance is now slack, kept as-is. 9.626 wired
    // `ShapeFix_IntersectionTool::FixSelfIntersectWire`: 67526 vs 67366
    // (v 60641 vs 60548).
    // Re-measured 2026-10-12 after A13/A18/A19/A23 (fallback meshers deleted,
    // `face_discret` interior grid dropped, `A_ITERATIONS_NB = 11`): 73150 vs
    // 67366 (v 63522 vs 60548), area 192668.27 vs 192658.64 (ratio 1.0001).
    // The +8.6% triangle delta is **not** from the refinement cap: an A/B with
    // `A_ITERATIONS_NB = 4` gives 73162 vs 73150 (0.02%). It comes from faces
    // that used to fall back to the UV-grid mesher now going through Delaunay.
    // 2026-10-13: the missing `ShapeFix_Shape::Perform` tail of the STEP
    // `FixShape` pass (`ShapeFix_Shape.cxx:259-261` `SameParameter` +
    // `:267-291` `FixVertexTolerance`, wired in
    // `step/read_topology/resolver_context.rs::close_heal_context`) collapsed the delta: 67426 vs
    // 67366 (v 60576 vs 60548), area 192657.72 vs 192658.64 (ratio 1.0000). The
    // old +5784 triangles were one class: the port left `SameParameter=false`
    // on the torus face's u-iso BSpline edges, so `BRepMesh_EdgeDiscret::
    // CreateEdgeTessellator` (`cxx:57-59`) picked the `Adaptor3d_CurveOnSurface`
    // overload (`myType==GeomAbs_Circle` from `EvalKPart`) instead of the 3D
    // curve, and `GCPnts_TangentialDeflection::PerformCircular` sampled 17
    // uniform V parameters where OCCT's adaptive `PerformCurve` yields 11
    // (face 4: nvp 28->11, nv 17->11, nodes 525->315, mesh 629/1152 -> 407/720,
    // matching the OCCT probe row bit for bit).
    // 2026-10-13 (seam registration, see occ/TDB.stp): 67424 vs 67366
    // (v 60575 vs 60548).
    // 2026-10-16 (re-read on the frozen tranche, `.target-gate/t157/gates-raw.txt`):
    // 66897 vs 67366 (v 60199 vs 60548), area ratio 1.0017. The 2026-10-13 reading is kept as
    // history: T0M moved under the 10-15..10-16 packages that made the shared intersection /
    // extrema / root paths faithful. The tolerance below is the registered one and was not
    // widened; the move is away from OCCT, which matches the "port still misses some face
    // classes (density < 1)" gap this model tracks.
    // 2026-10-17 (`Extrema_ExtPElC.cxx:253-260`: the point-ellipse root equation
    // used `+B*Y*cos` where OCCT has `ko3 = -B*Y`; `poly_roots::ellipse_all` now
    // matches `ElCLib::EllipseValue`'s `+B*sin`): 67286 vs 67366 (v 60506 vs
    // 60548), area 192657.74 vs 192658.64 (ratio 1.0000). The old +321.92 came
    // from `ADVANCED_FACE #29213` (cylinder `#1514`, strip u in [-0.794,0.794]):
    // its 3.2e-7-long ELLIPSE edge `#15786` had both vertices projecting to the
    // mirrored root (+1.52e-4 instead of 2*pi-1.52e-4), so `UpdateParam3d`'s
    // `AdjustPeriodic` (`GeometricTool.cxx:274-281`) turned a 1.5e-4 arc into a
    // whole-period edge whose pcurve swept the full cylinder and widened
    // `BRepTools::UVBounds` to 1.25 periods (meshed face 507.11 vs 183.00).
    Model { step: "occ/T0M.stp", occ: "occ/occ-T0M.obj", bbox_tol: Some(0.9), area_tol: Some(0.025) },
    // acs10: widened 0.05 -> 0.12 by the same T-93 (a) change (ratio 0.9846 ->
    // 0.9025, ours 244773.63 / occ 271219.68). Re-measured at 9.625: 1.00004
    // (271230.61 / 271219.68), 46548 vs 46558 triangles. 9.626: 46544 vs 46558
    // (v 37287 vs 37296).
    // 2026-10-12 (A13/A18/A19/A23): 47196 vs 46558 (v 37606 vs 37296), area
    // 271232.13 vs 271219.68.
    // 2026-10-13 (seam registration, see occ/TDB.stp): 46538 vs 46558
    // (v 37280 vs 37296).
    Model { step: "occ/acs10.stp", occ: "occ/occ-acs10.obj", bbox_tol: Some(1e-3), area_tol: Some(0.12) },
    // TDB: 9.626 (FixSelfIntersectWire): 79912 vs 79977 faces
    // (v 73702 vs 73544).
    // 2026-10-12 (A13/A18/A19/A23): 81721 vs 79977 (v 74792 vs 73544), area
    // 235866.32 vs 234884.51 (ratio 1.0042).
    // 2026-10-13: the +1744 triangles were two seam-edge registration gaps.
    // `BRepMesh_ModelPreProcessor::splitEdge` evaluates the pcurve branch test
    // at `aPC1->FirstParameter()` (`BRepMesh_ModelPreProcessor.cxx:226-231`), so
    // a `Geom2d_Line` seam seed is evaluated at `-Precision::Infinite()` and the
    // two pcurves get cross-wired; the port used a finite midpoint and left the
    // seam edges ordered. `BRepMesh_BaseMeshAlgo::initDataStructure`
    // (`cxx:87-125`) then registers *every* pcurve of each wire edge instead of
    // the CheckOrder slot. Delta collapsed to 79926 vs 79977 (v 73516 vs
    // 73544), f-ratio 0.9994.
    // 2026-10-16 (re-read on the frozen tranche, `.target-gate/t157/gates-raw.txt`):
    // 79994 vs 79977 (v 73586 vs 73544), area ratio 1.0015. Same drift as T0M above (the
    // 2026-10-13 reading is kept as history); the registered tolerance was not widened.
    Model { step: "occ/TDB.stp", occ: "occ/occ-TDB.obj", bbox_tol: Some(1e-3), area_tol: Some(0.02) },
    // a3n00: the tracked T-59/T-69 face classes (density 0.72); the wide area
    // tolerance encodes the missing coverage, not a meshing approximation.
    // 2026-10-13 (seam registration, see occ/TDB.stp): 12442 vs 12466
    // (v 11138 vs 11145).
    Model { step: "occ/a3n00.stp", occ: "occ/occ-a3n00.obj", bbox_tol: Some(0.01), area_tol: Some(0.15) },
];

/// Read step, export the whole document to `(OBJ text, MTL text)` at
/// deflection 0.1.
///
/// This is the shipped exporter (`RWObj_CafWriter::Perform` via
/// [`model_to_obj_caf`]), not the plain [`brep_to_obj`] shape writer, so the
/// gates cover the same output `examples/export_data_obj.rs` produces: header,
/// `g <node>`, `mtllib` / `usemtl` / `newmtl`. The plain writer emits none of
/// those records, which is why the `.mtl` sidecar used to be untested.
///
/// Shape selection matches the exporter too: a single-root document exports its
/// root, a multi-root document is meshed as one compound of the roots.
pub fn export(step: &str) -> (String, String) {
    let model = read_step_file(&data_dir().join(step).to_string_lossy())
        .unwrap_or_else(|e| panic!("read {}: {e}", step));
    let whole = document_shape(&model)
        .unwrap_or_else(|| panic!("{}: parsed to no shapes", step));
    model_to_obj_caf(&model, &whole, &format!("{}.obj", stem(step)), 0.1)
}

/// `mtllib` record of an OBJ (the sidecar file name), if present.
pub fn obj_mtllib(obj: &str) -> Option<&str> {
    obj.lines().find_map(|l| l.strip_prefix("mtllib ")).map(str::trim)
}

/// `usemtl` names referenced by an OBJ, in file order (duplicates kept).
pub fn obj_usemtl(obj: &str) -> Vec<&str> {
    obj.lines().filter_map(|l| l.strip_prefix("usemtl ")).map(str::trim).collect()
}

/// `newmtl` names declared by an MTL library, in file order.
pub fn mtl_newmtl(mtl: &str) -> Vec<&str> {
    mtl.lines().filter_map(|l| l.strip_prefix("newmtl ")).map(str::trim).collect()
}

/// The OCCT reference OBJ text.
pub fn occ_text(occ: &str) -> String {
    std::fs::read_to_string(data_dir().join(occ)).unwrap_or_else(|_| panic!("missing {occ}"))
}

/// (vertices, faces) recorded in an OCCT WriteObj header.
pub fn occ_counts(obj: &str) -> (usize, usize) {
    let mut v = 0;
    let mut f = 0;
    for line in obj.lines() {
        if let Some(rest) = line.strip_prefix("#  Vertices:") {
            v = rest.trim().parse().unwrap_or(0);
        } else if let Some(rest) = line.strip_prefix("#     Faces:") {
            f = rest.trim().parse().unwrap_or(0);
        }
    }
    (v, f)
}

/// (vertices, faces) as written in the OBJ body.
pub fn obj_counts(obj: &str) -> (usize, usize) {
    (
        obj.lines().filter(|l| l.starts_with('v')).count(),
        obj.lines().filter(|l| l.starts_with('f')).count(),
    )
}

/// (min, max) of an OBJ v x y z lines.
pub fn obj_bbox(obj: &str) -> Option<([f64; 3], [f64; 3])> {
    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    let mut n = 0;
    for l in obj.lines() {
        if let Some(rest) = l.strip_prefix("v ") {
            let mut it = rest.split_whitespace();
            let x: f64 = it.next()?.parse().ok()?;
            let y: f64 = it.next()?.parse().ok()?;
            let z: f64 = it.next()?.parse().ok()?;
            min[0] = min[0].min(x);
            min[1] = min[1].min(y);
            min[2] = min[2].min(z);
            max[0] = max[0].max(x);
            max[1] = max[1].max(y);
            max[2] = max[2].max(z);
            n += 1;
        }
    }
    if n == 0 { None } else { Some((min, max)) }
}

/// Total triangle area (fan rule around the first vertex).
pub fn obj_area(m: &ObjMesh) -> f64 {
    let mut area = 0.0;
    for f in &m.faces {
        let vs: Vec<usize> = f.v.iter().map(|&i| i as usize).collect();
        if vs.len() < 3 {
            continue;
        }
        let a = &m.vertices[vs[0]];
        for k in 1..vs.len() - 1 {
            let b = &m.vertices[vs[k]];
            let c = &m.vertices[vs[k + 1]];
            let ab = b.coord.subtracted(&a.coord);
            let ac = c.coord.subtracted(&a.coord);
            area += 0.5 * ab.crossed(&ac).modulus();
        }
    }
    area
}

/// `data/output/gates/<stem>.obj` for a data/-relative STEP path.
///
/// The gates write to their own subdirectory: `examples/export_data_obj.rs`
/// owns `data/output/<stem>.obj` + `<stem>.mtl`, and both producers must not
/// overwrite each other's artifacts.
pub fn output_path(step: &str) -> PathBuf {
    data_dir().join("output").join("gates").join(format!("{}.obj", stem(step)))
}

/// data/-relative STEP path -> file stem.
pub fn stem(step: &str) -> String {
    std::path::Path::new(step)
        .file_stem()
        .expect("stem")
        .to_string_lossy()
        .to_string()
}

/// Write the exported OBJ + its `.mtl` sidecar to `data/output/gates/`
/// (deflection 0.1, for visual comparison against the OCCT reference). An empty
/// `mtl` removes a sidecar left by an earlier run instead of writing a stale
/// library.
pub fn write_output(step: &str, obj: &str, mtl: &str) {
    let path = output_path(step);
    std::fs::create_dir_all(path.parent().unwrap()).expect("create data/output/gates");
    std::fs::write(&path, obj).unwrap_or_else(|e| panic!("write {}: {e}", path.display()));
    let mat = material_file_path(&path.to_string_lossy());
    if mtl.is_empty() {
        let _ = std::fs::remove_file(&mat);
    } else {
        std::fs::write(&mat, mtl).unwrap_or_else(|e| panic!("write {mat}: {e}"));
    }
}
    // acs10: widened 0.05 -> 0.12 by the same T-93 (a) change (ratio 0.9846 ->
    // 0.9025, ours 244773.63 / occ 271219.68).
