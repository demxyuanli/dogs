//! Phase 4 module: vrml — VRML 2.0 export.
//!
//! Port of OCCT's `VrmlAPI_Writer` output for the common case of an
//! `IndexedFaceSet`. Emits a single `Shape` with a `Coordinate` block and a
//! `coordIndex` triangle list (each face terminated by `-1`), optionally
//! nested inside a `Transform` when a shape carries a location translation.

use crate::mesh::ShapeMesh;
use crate::model::BRepModel;
use crate::shape::TopoShape;

/// Serialize `mesh` as a VRML 2.0 file with a single named `Shape`.
pub fn write_vrml(mesh: &ShapeMesh, name: &str) -> String {
    let mut s = String::new();
    s.push_str("#VRML V2.0 utf8\n");
    s.push_str(&format!("# {name}\n"));
    s.push_str(&vrml_shape_body(mesh, name, ""));
    s
}

/// Mesh `shape` and serialize it as VRML 2.0.
pub fn write_vrml_shape(shape: &TopoShape, name: &str, deflection: f64) -> String {
    write_vrml(&crate::shape_mesh::mesh_shape(shape, deflection), name)
}

/// Serialize a whole model: one `Shape` per entry, wrapped in a `Transform`
/// when the entry's shape carries a location transform.
pub fn write_vrml_scene(model: &BRepModel, deflection: f64) -> String {
    let mut s = String::new();
    s.push_str("#VRML V2.0 utf8\n");
    for entry in &model.shapes {
        let mesh = crate::shape_mesh::mesh_shape(&entry.shape, deflection);
        if entry.shape.location.is_identity() {
            s.push_str(&vrml_shape_body(&mesh, &entry.name, ""));
        } else {
            // Only the translation of the accumulated location is exported;
            // rotations/scales are ignored by this lightweight writer.
            let tr = entry.shape.location.transformation().translation_part();
            s.push_str("Transform {\n");
            s.push_str(&format!(
                "  translation {:.6} {:.6} {:.6}\n",
                tr.x, tr.y, tr.z
            ));
            s.push_str("  children [\n");
            s.push_str(&vrml_shape_body(&mesh, &entry.name, "    "));
            s.push_str("  ]\n");
            s.push_str("}\n");
        }
    }
    s
}

/// Write VRML text to a file.
pub fn write_vrml_file(path: &str, content: &str) -> std::io::Result<()> {
    std::fs::write(path, content)
}

/// Number of triangles in a VRML document, counting the `-1` `coordIndex`
/// face separators (test helper). Tokens are split on whitespace and commas,
/// so a coordinate like `-1.000000` is not counted.
pub fn vrml_triangle_count(content: &str) -> usize {
    content
        .split(|c: char| c.is_whitespace() || c == ',')
        .filter(|t| *t == "-1")
        .count()
}

/// The `Shape { geometry IndexedFaceSet { ... } }` block, indented by `indent`.
fn vrml_shape_body(mesh: &ShapeMesh, name: &str, indent: &str) -> String {
    let mut s = String::new();
    s.push_str(&format!("{indent}Shape {{\n"));
    s.push_str(&format!("{indent}  # {name}\n"));
    s.push_str(&format!("{indent}  geometry IndexedFaceSet {{\n"));
    s.push_str(&format!("{indent}    coord Coordinate {{\n"));
    s.push_str(&format!("{indent}      point [\n"));
    for v in &mesh.vertices {
        s.push_str(&format!(
            "{indent}        {:.6} {:.6} {:.6},\n",
            v.x(),
            v.y(),
            v.z()
        ));
    }
    s.push_str(&format!("{indent}      ]\n"));
    s.push_str(&format!("{indent}    }}\n"));
    s.push_str(&format!("{indent}    coordIndex [\n"));
    for t in &mesh.triangles {
        s.push_str(&format!("{indent}        {}, {}, {}, -1,\n", t.n0, t.n1, t.n2));
    }
    s.push_str(&format!("{indent}    ]\n"));
    s.push_str(&format!("{indent}    normalPerVertex FALSE\n"));
    s.push_str(&format!("{indent}    solid TRUE\n"));
    s.push_str(&format!("{indent}  }}\n"));
    s.push_str(&format!("{indent}}}\n"));
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use occt_core::gp::GpPnt;

    #[test]
    fn write_vrml_box_mesh() {
        let mesh = crate::mesh::mesh_box((GpPnt::zero(), GpPnt::new(1.0, 1.0, 1.0)));
        let s = write_vrml(&mesh, "box");
        assert!(s.starts_with("#VRML V2.0 utf8"));
        assert!(s.contains("IndexedFaceSet"));
        assert!(s.contains("coord"));
        assert!(s.contains("Coordinate"));
        assert_eq!(vrml_triangle_count(&s), 12);
        // The point section has no bare "-1" face separators.
        let point_section: Vec<&str> = s
            .lines()
            .skip_while(|l| !l.contains("point ["))
            .skip(1)
            .take_while(|l| !l.contains(']'))
            .collect();
        assert!(point_section.iter().all(|l| !l.contains(", -1,")));
    }

    #[test]
    fn write_vrml_shape_of_sphere_nonempty() {
        let sphere = crate::primitives::BRepPrimSphere::make_sphere(1.0).solid;
        let s = write_vrml_shape(&sphere, "sphere", 0.2);
        assert!(s.contains("#VRML V2.0 utf8"));
        assert!(s.contains("IndexedFaceSet"));
        assert!(vrml_triangle_count(&s) > 0);
    }

    #[test]
    fn scene_writes_one_shape_per_entry() {
        let mut model = BRepModel::new();
        let box_shape = crate::primitives::BRepPrimBox::make_box(1.0, 1.0, 1.0).solid;
        model.add("box1", box_shape.clone().into());
        let moved = crate::transform::translated(&box_shape, &occt_core::gp::GpVec::new(2.0, 0.0, 0.0));
        model.add("box2", moved);
        let s = write_vrml_scene(&model, 0.25);
        assert!(s.starts_with("#VRML V2.0 utf8"));
        assert_eq!(s.matches("IndexedFaceSet").count(), 2);
        assert!(s.contains("Transform {"));
        assert!(s.contains("translation 2.000000 0.000000 0.000000"));
        assert!(vrml_triangle_count(&s) >= 2);
    }

    #[test]
    fn vrml_file_roundtrip() {
        let mesh = crate::mesh::mesh_box((GpPnt::zero(), GpPnt::new(1.0, 1.0, 1.0)));
        let s = write_vrml(&mesh, "box");
        let path = std::env::temp_dir().join("occt_topo_test_vrml.wrl");
        let p = path.to_str().unwrap();
        write_vrml_file(p, &s).unwrap();
        let read = std::fs::read_to_string(p).unwrap();
        assert_eq!(read, s);
        std::fs::remove_file(p).ok();
    }
}
