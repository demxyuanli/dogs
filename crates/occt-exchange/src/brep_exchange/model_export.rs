use crate::model::BRepModel;
use occt_topo::shape::TopoShape;

use super::caf_nodes::obj_caf_nodes;
use super::gltf::gltf_caf_nodes;

/// Export a document (`BRepModel`) through `RWGltf_CafWriter::Perform`
/// (`RWGltf_CafWriter.cxx:636-667`), returning the `.gltf`/`.glb` document and
/// the `.bin` body.
///
/// `gltf_path` supplies the `.bin` sidecar name in `buffers[0].uri`
/// (`cxx:262-264`); nothing is written. The coordinate converter in `options`
/// carries the glTF output system and the metre unit, and no input length unit
/// is set because the port's `BRepModel` has no
/// `XCAFDoc_DocumentTool::GetLengthUnit` counterpart — callers whose geometry is
/// in millimetres set it through `options.cs_trsf.set_input_length_unit(0.001)`,
/// as `RWGltf_CafReader` does on the read side (`rwmesh::gltf_file`).
pub fn model_to_gltf_caf(
    model: &BRepModel,
    whole: &TopoShape,
    gltf_path: &str,
    maximal_chordial_deviation: f64,
    options: &crate::rwmesh::GltfCafWriterOptions,
) -> Result<crate::rwmesh::GltfCafWriterOutput, String> {
    let nodes = gltf_caf_nodes(model, whole, maximal_chordial_deviation);
    // `RWGltf_CafWriter::Perform` (`cxx:645-649`): the document's length unit
    // feeds the converter before the writer runs; without one the converter
    // stays unscaled (`GetLengthUnit` FALSE).
    let mut options = options.clone();
    if let Some(length_unit) = model.length_unit {
        options.cs_trsf.set_input_length_unit(length_unit);
    }
    // `XSDRAWGLTF::WriteGltf` (`XSDRAWGLTF.cxx:307,505`): the Draw harness the
    // `data/occ-*.gltf` baselines come from declares the document frame as
    // Z-up, so the converter rotates Z-up -> glTF Y-up. `RWGltf_CafWriter`
    // leaves the input system unset; the caller supplies it.
    options
        .cs_trsf
        .set_input_coordinate_system(crate::rwmesh::RwMeshCoordinateSystem::ZUP);
    crate::rwmesh::perform_gltf(gltf_path, &nodes, &[], &options)
}

/// Write a document as `.gltf` plus its `.bin` sidecar, or as `.glb`
/// (`RWGltf_CafWriter::Perform`).
pub fn model_write_gltf_caf(
    path: &str,
    model: &BRepModel,
    whole: &TopoShape,
    maximal_chordial_deviation: f64,
    options: &crate::rwmesh::GltfCafWriterOptions,
) -> std::io::Result<()> {
    let out = model_to_gltf_caf(model, whole, path, maximal_chordial_deviation, options)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
    std::fs::write(path, &out.document)?;
    if !out.bin.is_empty() {
        std::fs::write(
            crate::rwmesh::bin_file_path(path, &out.bin_file_name_short),
            &out.bin,
        )?;
    }
    Ok(())
}

/// Export a document (`BRepModel`) through `RWPly_CafWriter::Perform`
/// (`RWPly_CafWriter.cxx:70-183`), returning the ASCII PLY text.
///
/// The node list is the one [`obj_caf_nodes`] builds for
/// `RWObj_CafWriter::Perform`: `XCAFPrs_DocumentExplorer(OnlyLeafNodes)` leaves
/// with the label-keyed styles, so both writers see the same faces, styles and
/// colors. `ply_path` supplies the context name only (the PLY header has no
/// file-name record).
pub fn model_to_ply_caf(
    model: &BRepModel,
    whole: &TopoShape,
    ply_path: &str,
    maximal_chordial_deviation: f64,
    options: &crate::rwmesh::PlyCafWriterOptions,
) -> Result<String, String> {
    let nodes = obj_caf_nodes(model, whole, maximal_chordial_deviation);
    crate::rwmesh::perform_ply(ply_path, &nodes, &[], options)
}

/// Write a document to an ASCII PLY file (`RWPly_CafWriter`).
pub fn model_write_ply_caf(
    path: &str,
    model: &BRepModel,
    whole: &TopoShape,
    maximal_chordial_deviation: f64,
) -> std::io::Result<()> {
    let ply = model_to_ply_caf(
        model,
        whole,
        path,
        maximal_chordial_deviation,
        &crate::rwmesh::PlyCafWriterOptions::default(),
    )
    .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
    std::fs::write(path, ply)
}

/// The single shape a document is meshed as, `XGetOneShape` in DRAW
/// (`data/_occ_ref_export.tcl`): the root itself for a single-root document,
/// a compound of the roots otherwise (`BRep_Builder::MakeCompound`), `None` for
/// an empty document.
///
/// Every read -> mesh -> export driver uses this (`examples/export_data_obj.rs`,
/// `examples/iges_check.rs`, the `step_obj_gates` model table), so a gate cannot
/// export a different part of the document than the shipped exporter.
pub fn document_shape(model: &BRepModel) -> Option<TopoShape> {
    match model.shapes.as_slice() {
        [] => None,
        [one] => Some(one.shape.clone()),
        many => {
            let parts: Vec<TopoShape> = many.iter().map(|s| s.shape.clone()).collect();
            Some(occt_topo::builder::TopoBuilder::new().make_compound_of(&parts).0)
        }
    }
}

/// Export a document (`BRepModel`) through `RWObj_CafWriter::Perform`
/// (`RWObj_CafWriter.cxx:82-230`), returning `(obj_text, mtl_text)`.
///
/// `obj_path` only supplies the `.mtl` sidecar name (`cxx:150-157`); nothing is
/// written. `mtl_text` is empty when no node carries a color
/// (`toCreateMatFile == false`), i.e. no `mtllib` record is emitted.
pub fn model_to_obj_caf(
    model: &BRepModel,
    whole: &TopoShape,
    obj_path: &str,
    maximal_chordial_deviation: f64,
) -> (String, String) {
    let nodes = obj_caf_nodes(model, whole, maximal_chordial_deviation);
    crate::rwmesh::write_obj(obj_path, &nodes)
}

/// Write a document to an OBJ file plus its `.mtl` sidecar when any face carries
/// a color (`RWObj_CafWriter::Perform` + `RWObj_ObjMaterialMap`).
pub fn model_write_obj_caf(
    path: &str,
    model: &BRepModel,
    whole: &TopoShape,
    maximal_chordial_deviation: f64,
) -> std::io::Result<()> {
    let (obj, mtl) = model_to_obj_caf(model, whole, path, maximal_chordial_deviation);
    std::fs::write(path, obj)?;
    if !mtl.is_empty() {
        std::fs::write(crate::rwmesh::material_file_path(path), mtl)?;
    }
    Ok(())
}

