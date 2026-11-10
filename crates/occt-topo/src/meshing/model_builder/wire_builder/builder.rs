use super::add_wire::add_wire;
use super::wire_edges::outer_of_wires;

use super::super::prelude::*;

/// Tool for building a discrete model from a topological shape.
/// Port of `BRepMesh_ModelBuilder` (+ the visiting logic of
/// `BRepMesh_ShapeVisitor`).
pub struct ModelBuilder;

impl ModelBuilder {
    /// Build the discrete model of `shape` under `params`.
    ///
    /// Mirrors `BRepMesh_ModelBuilder::performInternal` plus the
    /// `BRepMesh_ShapeVisitor` walk: every face gets its outer wire first
    /// (`ShapeAnalysis::OuterWire`), then its inner wires, each wire reordered
    /// into a connected chain (`ShapeAnalysis_Wire::CheckOrder`), with EXTERNAL
    /// edges skipped and shared edges deduplicated by TShape identity. Free edges
    /// (not bounding any face) are added as well.
    ///
    /// Returns `Err` when the shape is empty (void bounding box → `Message_Fail1`).
    pub fn build_model(shape: &TopoShape, params: &MeshParameters) -> Result<MeshModel, String> {
        let mut model = MeshModel::new(shape.clone());

        // Maximum size of the shape's bounding box (used by relative deflection).
        let bbox = shape_bbox(shape);
        if bbox.is_void() {
            return Err("BRepMesh_ModelBuilder::build_model: empty shape".to_string());
        }
        let max_size = if params.relative {
            box_max_dimension(&bbox).unwrap_or(0.0)
        } else {
            params.deflection.max(params.deflection_interior.max(0.0))
        };
        model.set_max_size(max_size);

        // Visit(Face): `IMeshTools_ShapeExplorer.cxx:101` stores only
        // `Oriented(FORWARD)` faces "to prevent inverse issue". A REVERSED
        // face would Compose hole wires to FORWARD and fill the hole.
        // Outer wire first; a failure on the outer wire fails the face, a
        // failure on an inner wire only marks it unoriented.
        let mut edge_index: HashMap<usize, usize> = HashMap::new();
        for f in faces_of(shape) {
            if crate::debug_flags::sweep_dbg() {
                let s = BRepTool::face_surface(&f);
                let sd = match &s {
                    Some(s) => {
                        let o = s.d0(0.0, 0.0);
                        let (_, du, dv) = s.d1(0.0, 0.0);
                        format!(
                            "o=({:.3},{:.3},{:.3}) du=({:.3},{:.3},{:.3}) dv=({:.3},{:.3},{:.3})",
                            o.x(), o.y(), o.z(), du.x(), du.y(), du.z(), dv.x(), dv.y(), dv.z()
                        )
                    }
                    None => "none".to_string(),
                };
                eprintln!(
                    "[model-dbg] face tsh={:?} loc={:?} surf={}",
                    Arc::as_ptr(&f.0.tshape) as usize,
                    f.0.location().transformation(),
                    sd
                );
            }
            let f = Face(f.0.oriented(Orientation::Forward));
            let face_index = model.add_face(f.clone());

            // `BRepMesh_ShapeVisitor::Visit` (`BRepMesh_ShapeVisitor.cxx:58-87`)
            // applies NO orientation fixup: it is `ShapeAnalysis::OuterWire` +
            // `TopExp_Explorer(face, TopAbs_WIRE)`, and `addWire` only runs
            // `ShapeAnalysis_Wire::CheckOrder`. `ShapeFix_Face::FixOrientation`
            // belongs to the Shape Healing pass the reader already ran
            // (`IGESControl_IGESBoundary` -> `ShapeFix_Face::Perform`,
            // port: `iges_read/iges_process_shape.rs`); repeating it inside the
            // mesher reversed the torus' `V in [-pi, pi]` wire and made the
            // FaceChecker report a self-intersection (`FaceChecker.cxx:150-158`)
            // so the face came out with zero triangles.
            let wires = wires_of_face(&f);
            let outer = outer_of_wires(&wires, &f);
            if let Some(outer_wire) = &outer {
                if !add_wire(&mut model, face_index, outer_wire, &mut edge_index) {
                    if let Ok(fm) = model.face_mut(face_index) {
                        fm.set_status(MeshStatus::FAILURE);
                    }
                    continue;
                }
            }

            for w in &wires {
                if let Some(outer_wire) = &outer {
                    if Arc::ptr_eq(&w.0.tshape, &outer_wire.0.tshape) {
                        continue;
                    }
                }
                if !add_wire(&mut model, face_index, w, &mut edge_index) {
                    if let Ok(fm) = model.face_mut(face_index) {
                        fm.set_status(MeshStatus::UNORIENTED_WIRE);
                    }
                }
            }
        }

        // Free edges that do not bound any face.
        for e in edges_of(shape) {
            let key = Arc::as_ptr(&e.0.tshape) as usize;
            if !edge_index.contains_key(&key) {
                model.add_edge(e.clone());
            }
        }

        Ok(model)
    }
}
