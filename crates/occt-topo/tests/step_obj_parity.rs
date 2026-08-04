//! STEP → OBJ parity tests against the OCCT reference exports in `data/`.
//!
//! For every sample STEP that has an OCCT-generated `occ-*.obj` reference,
//! this test converts the parsed shape to OBJ and compares:
//!
//!  1. **Bounding box** — the geometry extent must match the OCCT reference
//!     to within a tolerance. This is the hard geometric gate: a correct STEP
//!     parser + surface evaluation gives the same extents as OCCT.
//!  2. **Mesh density** — the OCCT reference headers carry `# Vertices` /
//!     `# Faces`; we record the ratio but do not assert it (the Rust mesh
//!     tessellator is a UV grid, OCCT's is deflection-adaptive — the density
//!     gap is a known, tracked gap, not a geometry error).
//!
//! The pipeline under test: `step::read_step_file` → `brep_exchange::brep_to_obj`.

use std::path::PathBuf;

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

/// Parse the `# Vertices: N` / `# Faces: N` header lines of an OCCT OBJ.
struct OccHeader {
    vertices: usize,
    faces: usize,
}

fn parse_occ_header(obj: &str) -> OccHeader {
    let mut v = 0;
    let mut f = 0;
    for line in obj.lines() {
        if let Some(rest) = line.strip_prefix("#  Vertices:") {
            v = rest.trim().parse().unwrap_or(0);
        } else if let Some(rest) = line.strip_prefix("#     Faces:") {
            f = rest.trim().parse().unwrap_or(0);
        }
    }
    OccHeader { vertices: v, faces: f }
}

/// (min, max) bounding box of an OBJ's `v x y z` lines.
fn obj_bbox(obj: &str) -> Option<([f64; 3], [f64; 3])> {
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
    if n == 0 {
        None
    } else {
        Some((min, max))
    }
}

/// Compare our STEP→OBJ output to the OCCT reference: bbox must match
/// (within `tol`), density is recorded.
fn check_parity(step: &str, occ: &str, tol: f64) -> OccHeader {
    let model = read_step_file(&data_dir().join(format!("{step}.step")).to_string_lossy())
        .unwrap_or_else(|e| panic!("read {step}.step: {e}"));
    assert!(!model.shapes.is_empty(), "{step}: parsed to no shapes");
    let obj = brep_to_obj(&model.shapes[0].shape, 0.1);

    let (our_min, our_max) = obj_bbox(&obj).expect("our OBJ has vertices");

    let occ_path = data_dir().join(occ);
    let occ_obj = std::fs::read_to_string(&occ_path).unwrap_or_else(|_| panic!("missing {occ}"));
    let (occ_min, occ_max) = obj_bbox(&occ_obj).expect("occ OBJ has vertices");
    let header = parse_occ_header(&occ_obj);

    for i in 0..3 {
        assert!(
            (our_min[i] - occ_min[i]).abs() < tol,
            "{step}: bbox min[{i}] ours={} occ={}",
            our_min[i],
            occ_min[i]
        );
        assert!(
            (our_max[i] - occ_max[i]).abs() < tol,
            "{step}: bbox max[{i}] ours={} occ={}",
            our_max[i],
            occ_max[i]
        );
    }

    let nv = obj.lines().filter(|l| l.starts_with('v')).count();
    let nf = obj.lines().filter(|l| l.starts_with('f')).count();
    eprintln!(
        "[{step}] ours v={nv} f={nf} | occ v={} f={} | f-ratio {:.2}",
        header.vertices,
        header.faces,
        nf as f64 / header.faces as f64
    );
    header
}

// Bbox tolerance per shape. A tight tolerance (1e-6) means the parsed geometry
// matches OCCT exactly; a loose one is a *tracked gap* — the tolerance encodes
// the current mismatch, and porting the underlying STEP surface reconstruction
// should tighten it toward 1e-6. See the KNOWN_GAPS list in the test body.
// Bbox tolerance per shape. The *tight* cases (1e-6) are polyhedral/linear
// shapes whose boundary vertices land exactly on the mesh. Curved shapes get a
// small sampling tolerance (0.05): OCCT's deflection mesh puts vertices on the
// surface at its own sample points, our UV grid at ours — both are on-surface,
// so a few hundredths of a unit apart at the bbox edge is sampling, not error.
// The *wide* tolerances are tracked PARSE gaps (see each comment): the shape's
// parameterisation/coverage is wrong, and porting the STEP surface
// reconstruction should tighten them toward the sampling tolerance.
const CUBE_TOL: f64 = 1e-6;
const CONE_TOL: f64 = 1e-6;
const CYLINDER_TOL: f64 = 0.05; // fixed: pcurve-bounded UV domain (wave A)
const TORUS_TOL: f64 = 0.05; // on-surface sampling difference
const SPHERE_TOL: f64 = 0.05; // on-surface sampling difference

#[test]
fn cube_bbox_matches_occt() {
    check_parity("Cube", "occ-cube.obj", CUBE_TOL);
}

#[test]
fn cone_bbox_matches_occt() {
    check_parity("Cone", "occ-cone.obj", CONE_TOL);
}

#[test]
fn cylinder_bbox_matches_occt() {
    check_parity("Cylinder", "occ-cylinder.obj", CYLINDER_TOL);
}

#[test]
fn torus_bbox_matches_occt() {
    check_parity("Torus", "occ-torus.obj", TORUS_TOL);
}

#[test]
fn sphere_bbox_matches_occt() {
    check_parity("Sphere", "occ-shpere.obj", SPHERE_TOL);
}

#[test]
fn rev_export_is_valid() {
    // rev.step is a 6-face solid (two quarter-cylinder walls + caps) while
    // occ-rev.obj is a simpler 3-face quarter-cylinder — the reference is a
    // different/older shape, so no bbox parity is asserted. Assert a valid
    // export and record density.
    let model = read_step_file(&data_dir().join("rev.step").to_string_lossy())
        .unwrap_or_else(|e| panic!("read rev.step: {e}"));
    assert!(!model.shapes.is_empty(), "rev.step parsed to no shapes");
    let obj = brep_to_obj(&model.shapes[0].shape, 0.1);
    let nv = obj.lines().filter(|l| l.starts_with('v')).count();
    let nf = obj.lines().filter(|l| l.starts_with('f')).count();
    assert!(nv > 0 && nf > 0, "rev OBJ degenerate: v={nv} f={nf}");
    let occ = std::fs::read_to_string(data_dir().join("occ-rev.obj")).unwrap();
    let h = parse_occ_header(&occ);
    eprintln!(
        "[rev] ours v={nv} f={nf} | occ-reference(v={} f={}) | f-ratio {:.2}",
        h.vertices,
        h.faces,
        nf as f64 / h.faces as f64
    );
}

#[test]
fn offset_export_is_valid() {
    // Offset.step has no matching OCCT reference (occ-OffsetPlaneHoleEdge is a
    // different file), so no bbox parity is asserted — only that it reads and
    // exports a non-degenerate mesh. Density is recorded for reference.
    let model = read_step_file(&data_dir().join("Offset.step").to_string_lossy())
        .unwrap_or_else(|e| panic!("read Offset.step: {e}"));
    assert!(!model.shapes.is_empty(), "Offset.step parsed to no shapes");
    let obj = brep_to_obj(&model.shapes[0].shape, 0.1);
    let nv = obj.lines().filter(|l| l.starts_with('v')).count();
    let nf = obj.lines().filter(|l| l.starts_with('f')).count();
    assert!(nv > 0 && nf > 0, "Offset OBJ degenerate: v={nv} f={nf}");
    let occ = std::fs::read_to_string(data_dir().join("occ-OffsetPlaneHoleEdge.obj")).unwrap();
    let h = parse_occ_header(&occ);
    eprintln!(
        "[Offset] ours v={nv} f={nf} | occ-reference(v={} f={}) | f-ratio {:.2}",
        h.vertices,
        h.faces,
        nf as f64 / h.faces as f64
    );
}
