//! STEP -> OBJ gates: bbox parity, mesh-area parity and export validity over
//! the shared 23-model table (tests/common/mod.rs).
//!
//! All checks live in one test binary and share a OnceLock cache of the
//! read+exported OBJ texts, so the 23 STEP files are parsed and tessellated once
//! per run instead of once per gate.

mod common;

use std::sync::OnceLock;

use common::{
    data_dir, export, obj_area, obj_bbox, obj_counts, occ_counts, occ_text, output_path,
    write_output, MODELS,
};
use occt_core::io::obj::{read_obj_file, ObjMesh};

/// One model, read and exported once, plus its OCCT reference text.
struct Case {
    step: &'static str,
    occ: &'static str,
    bbox_tol: Option<f64>,
    area_tol: Option<f64>,
    obj: String,
    reference: String,
}

/// Read + export every model exactly once for the whole test binary.
fn cases() -> &'static [Case] {
    static CASES: OnceLock<Vec<Case>> = OnceLock::new();
    CASES.get_or_init(|| {
        MODELS
            .iter()
            .map(|m| Case {
                step: m.step,
                occ: m.occ,
                bbox_tol: m.bbox_tol,
                area_tol: m.area_tol,
                obj: export(m.step),
                reference: occ_text(m.occ),
            })
            .collect()
    })
}

/// step_obj_parity: bounding box of the exported mesh vs the OCCT reference.
#[test]
fn step_obj_parity_bboxes_match_occt() {
    for c in cases() {
        let (nv, nf) = obj_counts(&c.obj);
        assert!(nv > 0 && nf > 0, "{}: degenerate OBJ v={nv} f={nf}", c.step);
        let (our_min, our_max) = obj_bbox(&c.obj).expect("our OBJ has vertices");
        let (occ_min, occ_max) = obj_bbox(&c.reference).expect("reference OBJ has vertices");
        let (ov, of) = occ_counts(&c.reference);
        eprintln!(
            "[{}] ours v={nv} f={nf} | occ v={ov} f={of} | f-ratio {:.2}",
            c.step,
            nf as f64 / of as f64
        );
        if let Some(tol) = c.bbox_tol {
            for i in 0..3 {
                assert!(
                    (our_min[i] - occ_min[i]).abs() < tol,
                    "{}: bbox min[{i}] ours={} occ={} (tol {tol})",
                    c.step, our_min[i], occ_min[i]
                );
                assert!(
                    (our_max[i] - occ_max[i]).abs() < tol,
                    "{}: bbox max[{i}] ours={} occ={} (tol {tol})",
                    c.step, our_max[i], occ_max[i]
                );
            }
        }
    }
}

/// step_obj_area: mesh surface area vs the OCCT reference.
#[test]
fn step_obj_area_matches_occt() {
    for c in cases() {
        let our = obj_area(&ObjMesh::parse(&c.obj).expect("parse our OBJ"));
        let occ = ObjMesh::parse(&c.reference).expect("parse reference OBJ");
        let occ_area = obj_area(&occ);
        let ratio = our / occ_area;
        eprintln!("[{}] our={our:.2} occ={occ_area:.2} ratio={ratio:.4}", c.step);
        if let Some(tol) = c.area_tol {
            assert!(
                (ratio - 1.0).abs() < tol,
                "{}: area ratio {ratio:.4} outside tolerance {tol} (ours {our:.2}, occ {occ_area:.2})",
                c.step
            );
        }
    }
}

/// step_to_obj: every export is a non-degenerate, finite OBJ; write it out.
#[test]
fn step_to_obj_exports_valid_obj() {
    for c in cases() {
        let (nv, nf) = obj_counts(&c.obj);
        assert!(nv > 0 && nf > 0, "{}: degenerate OBJ v={nv} f={nf}", c.step);
        let (min, max) = obj_bbox(&c.obj).expect("our OBJ has vertices");
        assert!(
            min.iter().chain(max.iter()).all(|v| v.is_finite()),
            "{}: non-finite bbox {min:?} {max:?}",
            c.step
        );
        write_output(c.step, &c.obj);
    }
}

/// The written OBJ round-trips through the OBJ reader. Writes to its own
/// subdirectory so it never races the export gate's data/output/<stem>.obj.
#[test]
fn outputs_round_trip_through_obj_reader() {
    let dir = data_dir().join("output").join("roundtrip");
    std::fs::create_dir_all(&dir).expect("create data/output/roundtrip");
    for c in cases() {
        let stem = std::path::Path::new(c.step)
            .file_stem()
            .expect("stem")
            .to_string_lossy()
            .to_string();
        let out = dir.join(format!("{stem}.obj"));
        std::fs::write(&out, &c.obj).unwrap_or_else(|e| panic!("write {}: {e}", out.display()));
        let mesh = read_obj_file(&out.to_string_lossy()).expect("round-trip OBJ parse");
        assert!(!mesh.vertices.is_empty(), "{}: round-trip lost vertices", c.step);
        assert!(!mesh.faces.is_empty(), "{}: round-trip lost faces", c.step);
    }
}

/// Regression: a full-cylinder face's UV u-domain must span the whole period.
/// face_uv_bounds samples each boundary-edge pcurve; the full circle's pcurve is
/// unwrapped monotonically (u 0 -> -2pi) so the samples only reached 7pi/4,
/// dropping a pi/4 wedge and ~12.5% of the lateral area (133.9 vs 150.8).
#[test]
fn cylinder_step_mesh_covers_full_u_period() {
    let model = occt_topo::step::read_step_file(
        &data_dir().join("Cylinder.step").to_string_lossy(),
    )
    .unwrap_or_else(|e| panic!("read Cylinder.step: {e}"));
    let shape = &model.shapes[0].shape;
    let im = occt_topo::brepmesh::incremental_mesh(shape, 0.1).expect("mesh cylinder");
    let area = occt_topo::mesh::mesh_surface_area(&im.mesh);
    // Analytic: 2*pi*r^2 + 2*pi*r*h with r=2, h=10 -> 48pi ~ 150.80. The mesh is a
    // chord approximation of the curved lateral, so allow 2% under.
    let expect = 48.0 * std::f64::consts::PI;
    assert!(
        (area - expect).abs() / expect < 0.02,
        "cylinder mesh area {area} vs analytic {expect}"
    );
}
