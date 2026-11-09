//! Public API: the `IGESControl_Writer` entry points.
use super::*;

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Serialize the whole model to an IGES file (ANSI Y14.26M).
pub fn write_iges(model: &BRepModel) -> String {
    let mut w = IgesWriter::new();
    for ms in &model.shapes {
        w.emit_shape(&ms.shape);
    }
    w.finish()
}

/// Serialize the model to an IGES file on disk.
pub fn write_iges_file(path: &str, model: &BRepModel) -> std::io::Result<()> {
    std::fs::write(path, write_iges(model))
}

/// Serialize a single shape (wrapped in a one-shape model).
pub fn write_shape_iges(shape: &TopoShape) -> String {
    let mut model = BRepModel::new();
    model.add("Shape", shape.clone());
    write_iges(&model)
}

/// Count `(geometry entities, total records)` for testing.
///
/// Geometry entities are directory entries whose type number falls in the
/// geometric range (100–199: point/line/arc/plane/revolution/trimmed/MSBO, or
/// 500–599: B-rep face/shell). Total records is the number of directory
/// entries (D-section lines / 2).
pub fn iges_entity_counts(text: &str) -> (usize, usize) {
    // The section letter sits in column 73 and the entity type in columns 1-8
    // (`IGESData_IGESWriter.cxx:836-859`).
    let d_lines: Vec<&str> = text
        .lines()
        .filter(|l| l.as_bytes().get(72) == Some(&b'D'))
        .collect();
    let total = d_lines.len() / 2;
    let mut geo = 0usize;
    for (i, l) in d_lines.iter().enumerate() {
        if i % 2 == 0 {
            let ty: i32 = l.get(0..8).unwrap_or("0").trim().parse().unwrap_or(0);
            if (100..200).contains(&ty) || (500..600).contains(&ty) {
                geo += 1;
            }
        }
    }
    (geo, total)
}
