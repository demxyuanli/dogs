use super::*;
use crate::primitives::{BRepPrimBox, BRepPrimSphere};


#[test]
fn box_iges_contains_expected_entities() {
    let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
    let iges = write_shape_iges(&b.solid.0);
    // The default `write.iges.brep.mode = 0` (Faces mode) writes the faces as
    // trimmed surfaces (144 over 142) grouped by 402 - not the 510/514/186
    // BRep-mode tree this port used to emit. Point entities (116) appear only
    // where OCCT references them (the location point of 192/194/196/198); the
    // port's earlier unconditional per-vertex 116 block was unreferenced output
    // and was removed in batch 71, so a box legitimately has none.
    for needle in ["110", "108", "142", "144", "402"] {
        assert!(iges.contains(needle), "missing entity {needle}");
    }
    for line in iges.lines() {
        assert_eq!(line.len(), 80, "line length {}: {:?}", line.len(), line);
    }
    let (geo, total) = iges_entity_counts(&iges);
    assert!(geo > 0 && total > 0, "entity counts ({geo}, {total})");
    let last = iges.lines().last().unwrap();
    assert!(last.contains("T0000001"), "terminate record: {last}");
}

#[test]
fn sphere_iges_has_arc_and_solid() {
    let s = BRepPrimSphere::make_sphere(2.0);
    let iges = write_shape_iges(&s.solid.0);
    assert!(
        iges.contains("100") || iges.contains("128"),
        "expected a circular arc or NURBS curve"
    );
    assert!(iges.contains("144"), "missing 144");
    for line in iges.lines() {
        assert_eq!(line.len(), 80, "line length {}", line.len());
    }
}

#[test]
fn iges_file_written() {
    let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
    let mut model = BRepModel::new();
    model.add("Box", b.solid.0.clone());
    let path = { let d = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("target").join("test-out"); let _ = std::fs::create_dir_all(&d); d }.join("occt_iges_test.igs");
    let p = path.to_str().unwrap();
    write_iges_file(p, &model).expect("write iges file");
    let content = std::fs::read_to_string(p).expect("read iges file");
    // The card format puts the data in columns 1-72, the section letter in
    // column 73 and the sequence number in columns 74-80
    // (`IGESData_IGESWriter.cxx:763-793` for the Start section).
    let first = content.lines().next().expect("first card");
    assert_eq!(first.len(), 80, "card length: {first:?}");
    assert_eq!(
        first.as_bytes().get(72),
        Some(&b'S'),
        "first card must be a Start card: {first:?}"
    );
    std::fs::remove_file(p).ok();
}

#[test]
fn section_terminators_present() {
    let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
    let iges = write_shape_iges(&b.solid.0);
    let last = iges.lines().last().unwrap();
    assert!(
        last.contains('S')
            && last.contains('G')
            && last.contains('D')
            && last.contains('P')
            && last.contains('T'),
        "terminate record missing section markers: {last}"
    );
}

#[test]
fn rec80_pads_and_truncates() {
    assert_eq!(rec80("abc").len(), 80);
    let long = "x".repeat(100);
    assert_eq!(rec80(&long).len(), 80);
}
