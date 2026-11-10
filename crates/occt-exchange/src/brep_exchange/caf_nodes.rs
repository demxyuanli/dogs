use std::collections::{HashMap, VecDeque};

use crate::model::BRepModel;
use crate::xcaf_prs::{dispatch_face_styles, XcafPrsStyle};
use occt_topo::shape::TopoShape;
use occt_topo::tgeometry::GeometryRegistry;
use occt_topo::topo_tools::faces_of;

use super::mesh_export::prs3d_get_deflection;
use super::style::{face_node, face_source};

/// Build the node list `RWObj_CafWriter::Perform` iterates for a document: one
/// node per `XCAFPrs_DocumentExplorer` leaf node (`BRepModel` shape) carrying
/// that node's `RWMesh_FaceIterator` faces with the style `RWMesh_ShapeIterator`
/// bound to each of them.
///
/// `whole` is the shape the single meshing pass runs on — DRAW's
/// `XGetOneShape s D` followed by `incmesh s <lin>` — so every node's faces are
/// read back from one mesh, exactly as `RWMesh_FaceIterator` reads the
/// `Poly_Triangulation` OCCT stored on the shared `TopoDS_Face`
/// (`data/_occ_ref_export.tcl:11-19`).
///
/// `RWMesh_FaceIterator` walks `TopExp_Explorer(shape, TopAbs_FACE)`, which this
/// crate's `faces_of` deduplicates by `TShape` identity (matching the mesh
/// model): a face shared by two nodes is written once, in the node that owns it
/// in mesh order.
///
/// A model filled from a document-style reader (`BRepModel::labels` non-empty)
/// is walked the way `RWObj_CafWriter::Perform` walks it
/// (`RWObj_CafWriter.cxx:186-221`): `XCAFPrs_DocumentExplorer` leaf nodes, the
/// Product label's name plus its parent names as the group name, and the
/// label-keyed `XCAFPrs::CollectStyleSettings` map for the per-face styles.
///
/// Locations are baked into the geometry here (`rwmesh::gltf_document`), so a
/// node's faces come from `DocumentNode::Label` — the instance label, whose
/// shape is the placed Product — instead of `RefLabel` plus
/// `DocumentNode::Location`; the styles and the group name still come from
/// `RefLabel`, as in `RWMesh_FaceIterator(aDocNode.RefLabel, aDocNode.Location,
/// ...)` (`RWObj_CafWriter.cxx:262-273`).
pub fn obj_caf_nodes(
    model: &BRepModel,
    whole: &TopoShape,
    maximal_chordial_deviation: f64,
) -> Vec<crate::rwmesh::ObjExportNode> {
    const PRS3D_DEV_ANGLE: f64 = 20.0 * std::f64::consts::PI / 180.0;
    let lin = prs3d_get_deflection(whole, maximal_chordial_deviation);
    if !model.labels.is_empty() && !model.roots.is_empty() {
        return obj_caf_nodes_from_document(model, whole, lin, PRS3D_DEV_ANGLE);
    }
    let inc = occt_topo::meshing::incremental_mesh::IncrementalMesh::from_deflection(
        whole,
        lin,
        false,
        PRS3D_DEV_ANGLE,
    );
    let stats = inc.face_stats();
    let tris = inc.face_triangulations();

    // Face shape key -> the mesh entries produced for it, in `MeshModel::face`
    // order. A face the mesher skipped (failure / reused) has no entry and is
    // skipped by the writer, like `toSkipFaceMesh` on a null triangulation
    // (`RWObj_CafWriter.cxx:74-77`).
    let mut by_key: HashMap<usize, VecDeque<usize>> = HashMap::new();
    for (i, st) in stats.iter().enumerate() {
        by_key.entry(st.shape_key).or_default().push_back(i);
    }

    let mut nodes = Vec::with_capacity(model.shapes.len());
    for ms in &model.shapes {
        // `XCAFPrs_DocumentNode::Style` (`XCAFPrs_DocumentExplorer.cxx:28-81`):
        // the node label's own color, which `RWMesh_ShapeIterator` keeps as
        // `myDefStyle`.
        let mut node_style = XcafPrsStyle::from_color(ms.color);
        // `XCAFPrs::CollectStyleSettings` (`XCAFPrs.cxx:174`, `cxx:251`) also
        // reads the label's visualization material:
        // `aStyle.SetMaterial(aMatTool->GetShapeMaterial(aLabel))`.
        node_style.set_material(model.shape_material(&ms.name).cloned());
        let face_styles = dispatch_face_styles(&ms.style_settings, &node_style);

        let mut faces = Vec::new();
        for face in faces_of(&ms.shape) {
            let key = GeometryRegistry::shape_key(&face.0);
            // `RWMesh_ShapeIterator::initShape` (`cxx:145-148`): the face's own
            // style, else the node style.
            let style = face_styles.get(&key).cloned().unwrap_or_else(|| node_style.clone());
            let Some(source) = face_source(&face, &mut by_key, tris) else { continue };
            faces.push(face_node(&face, source, style));
        }
        nodes.push(crate::rwmesh::ObjExportNode { name: ms.name.clone(), faces });
    }
    nodes
}

/// The export driver's tessellation stand-in for a document-style model: mesh
/// only the faces that carry no `BRep_Tool::Triangulation`, i.e. everything not
/// already filled by the reader (`readLateData`). Faces the document already
/// triangulated are never remeshed.
///
/// Returns the `MeshModel::face` entries keyed by face shape, plus the
/// triangulations they index into (`face_source`'s two inputs).
pub(super) fn mesh_missing_faces(
    leaf_faces: &[Vec<occt_topo::shape::Face>],
    lin: f64,
    ang: f64,
) -> (
    HashMap<usize, VecDeque<usize>>,
    Vec<occt_topo::meshing::triangulator::FaceTriangulation>,
) {
    let mut missing: Vec<TopoShape> = Vec::new();
    for faces in leaf_faces {
        for face in faces {
            if occt_topo::brep_tool::BRepTool::face_triangulation(face).is_some() {
                continue;
            }
            let key = GeometryRegistry::shape_key(&face.0);
            if missing.iter().any(|s| GeometryRegistry::shape_key(s) == key) {
                continue;
            }
            missing.push(face.0.clone());
        }
    }
    let mut by_key: HashMap<usize, VecDeque<usize>> = HashMap::new();
    let mut tris: Vec<occt_topo::meshing::triangulator::FaceTriangulation> = Vec::new();
    if !missing.is_empty() {
        let input = if missing.len() == 1 {
            missing[0].clone()
        } else {
            occt_topo::builder::TopoBuilder::new().make_compound_of(&missing).0
        };
        let inc = occt_topo::meshing::incremental_mesh::IncrementalMesh::from_deflection(&input, lin, false, ang);
        tris.extend_from_slice(inc.face_triangulations());
        for (i, st) in inc.face_stats().iter().enumerate() {
            by_key.entry(st.shape_key).or_default().push_back(i);
        }
    }
    (by_key, tris)
}

/// The document-style arm of [`obj_caf_nodes`]: `RWObj_CafWriter::Perform`
/// (`RWObj_CafWriter.cxx:186-221`) over `XCAFPrs_DocumentExplorer`.
fn obj_caf_nodes_from_document(
    model: &BRepModel,
    whole: &TopoShape,
    lin: f64,
    ang: f64,
) -> Vec<crate::rwmesh::ObjExportNode> {
    use crate::xcaf_prs_document_explorer::{DocumentExplorer, DocumentExplorerFlags};

    /// One explored leaf (`RWObj_CafWriter.cxx:192-208`).
    struct Leaf {
        name: String,
        /// `XCAFPrs_DocumentNode::Label`: the instance label. Its shape is the
        /// placed geometry, the port's stand-in for the node's
        /// `TopLoc_Location` (`rwmesh::gltf_document`, module docs).
        label: usize,
        /// `XCAFPrs_DocumentNode::RefLabel`: the Product label, the key of
        /// `XCAFPrs::CollectStyleSettings` and of `readNameAttribute`.
        ref_label: usize,
        style: XcafPrsStyle,
    }

    // `int aRootDepth = 0; if (theRootLabels.Size() == 1) ...` (`cxx:168-179`).
    let root_depth = if model.roots.len() == 1 {
        let ref_label = model.referred_label(model.roots[0]).unwrap_or(model.roots[0]);
        if model.label_name(ref_label).ends_with(".obj") {
            1
        } else {
            0
        }
    } else {
        0
    };

    let mut leaves: Vec<Leaf> = Vec::new();
    let mut explorer = DocumentExplorer::new(
        model,
        DocumentExplorerFlags::ONLY_LEAF_NODES,
        XcafPrsStyle::default(),
    );
    while explorer.more() {
        let node = explorer.current();
        // `aName = readNameAttribute(aDocNode.RefLabel)` plus the parent prefix
        // (`cxx:198-208`).
        let mut name = model.label_name(node.ref_label).to_string();
        let depth = explorer.current_depth();
        for parent in (root_depth..depth).rev() {
            let parent_name = explorer
                .current_at(parent)
                .map(|n| model.label_name(n.ref_label).to_string())
                .unwrap_or_default();
            if !parent_name.is_empty() {
                name = format!("{parent_name}/{name}");
            }
        }
        leaves.push(Leaf {
            name,
            label: node.label,
            ref_label: node.ref_label,
            style: node.style.clone(),
        });
        explorer.next();
    }

    // `RWMesh_FaceIterator` faces per node, in `faces_of` order.
    let leaf_faces: Vec<Vec<occt_topo::shape::Face>> = leaves
        .iter()
        .map(|leaf| model.label_shape(leaf.label).map(|s| faces_of(&s)).unwrap_or_default())
        .collect();
    let (mut by_key, tris) = mesh_missing_faces(&leaf_faces, lin, ang);
    let _ = whole;

    let mut nodes = Vec::with_capacity(leaves.len());
    for (leaf, faces) in leaves.iter().zip(leaf_faces) {
        // `RWMesh_ShapeIterator`'s style dispatch for the node's Product label
        // (`RWMesh_ShapeIterator.cxx:70-134`).
        let settings = crate::xcaf_prs::collect_style_settings(model, leaf.ref_label);
        let face_styles = dispatch_face_styles(&settings, &leaf.style);
        let mut out_faces = Vec::new();
        for face in &faces {
            let key = GeometryRegistry::shape_key(&face.0);
            let style = face_styles.get(&key).cloned().unwrap_or_else(|| leaf.style.clone());
            let Some(source) = face_source(face, &mut by_key, &tris) else { continue };
            out_faces.push(face_node(face, source, style));
        }
        nodes.push(crate::rwmesh::ObjExportNode { name: leaf.name.clone(), faces: out_faces });
    }
    nodes
}
