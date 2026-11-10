use crate::model::BRepModel;
use crate::xcaf_prs::{dispatch_face_styles, XcafPrsStyle};
use occt_topo::shape::TopoShape;
use occt_topo::tgeometry::GeometryRegistry;
use occt_topo::topo_tools::faces_of;

use super::caf_nodes::{mesh_missing_faces, obj_caf_nodes};
use super::mesh_export::prs3d_get_deflection;
use super::style::{face_node, face_source};

/// Build the node list `RWGltf_CafWriter::Perform` walks
/// (`RWGltf_CafWriter.cxx:664-666`): the whole `XCAFPrs_DocumentExplorer` tree
/// (flags `None`, `cxx:2286-2308`), assemblies included, each leaf carrying its
/// `RWMesh_FaceIterator` faces.
///
/// Unlike [`obj_caf_nodes`], a leaf keeps a face even when it has no
/// triangulation: `writeImages` / `writeMaterials` / `writeTextures` iterate the
/// face iterator without `toSkipShape` (`cxx:1948-1967`, `cxx:1997-2018`,
/// `cxx:2675-2690`), so an untriangulated face still contributes its style. The
/// empty node is what `RWMesh_FaceIterator::IsEmpty` means, and the mesh pass
/// skips it (`toSkipShape`).
///
/// A model without a label tree (no document-style reader filled it) becomes a
/// flat list of root nodes, one per `BRepModel::shapes` entry, the same shapes
/// [`obj_caf_nodes`] exports; its shape name stands for both `RWMesh::FormatName`
/// inputs because there is no instance/Product label to tell them apart.
pub fn gltf_caf_nodes(
    model: &BRepModel,
    whole: &TopoShape,
    maximal_chordial_deviation: f64,
) -> Vec<crate::rwmesh::GltfCafNode> {
    const PRS3D_DEV_ANGLE: f64 = 20.0 * std::f64::consts::PI / 180.0;
    let lin = prs3d_get_deflection(whole, maximal_chordial_deviation);
    if !model.labels.is_empty() && !model.roots.is_empty() {
        return gltf_caf_nodes_from_document(model, lin, PRS3D_DEV_ANGLE);
    }

    let nodes = obj_caf_nodes(model, whole, maximal_chordial_deviation);
    nodes
        .into_iter()
        .enumerate()
        .map(|(index, node)| {
            let name = model.shapes.get(index).map(|s| s.name.clone()).unwrap_or_default();
            let named_data = model.shapes.get(index).and_then(|s| s.named_data.clone());
            crate::rwmesh::GltfCafNode {
                id: crate::xcaf_prs_document_explorer::define_child_id(index, ""),
                label: index,
                ref_label: index,
                label_name: name.clone(),
                product_name: name,
                entry: index.to_string(),
                is_assembly: false,
                parent: None,
                local_trsf: occt_core::gp::GpTrsf::identity(),
                named_data,
                ref_named_data: None,
                faces: node.faces,
            }
        })
        .collect()
}

/// The document arm of [`gltf_caf_nodes`]: the full-tree walk of `writeNodes`
/// (`RWGltf_CafWriter.cxx:2286-2308`) plus the leaf faces of the binary/mesh
/// passes.
fn gltf_caf_nodes_from_document(
    model: &BRepModel,
    lin: f64,
    ang: f64,
) -> Vec<crate::rwmesh::GltfCafNode> {
    use crate::xcaf_prs_document_explorer::{DocumentExplorer, DocumentExplorerFlags};

    /// One explored node with the depth `CurrentDepth()` reports.
    struct ExploredNode {
        id: String,
        label: usize,
        ref_label: usize,
        /// `XCAFPrs_DocumentExplorer::CurrentDepth()` (`cxx:150-152`): an
        /// assembly's own depth, one past its parent.
        depth: usize,
        is_assembly: bool,
        /// `XCAFPrs_DocumentNode::Style`, the parent-merged style.
        style: XcafPrsStyle,
    }

    let mut explored: Vec<ExploredNode> = Vec::new();
    let mut explorer =
        DocumentExplorer::new(model, DocumentExplorerFlags::NONE, XcafPrsStyle::default());
    while explorer.more() {
        let node = explorer.current();
        explored.push(ExploredNode {
            id: node.id.clone(),
            label: node.label,
            ref_label: node.ref_label,
            depth: explorer.current_depth(),
            is_assembly: node.is_assembly,
            style: node.style.clone(),
        });
        explorer.next();
    }

    // `RWMesh_FaceIterator` faces per node, in `faces_of` order; assemblies own
    // no faces (`RWMesh_ShapeIterator` is only built for leaf nodes).
    let face_lists: Vec<Vec<occt_topo::shape::Face>> = explored
        .iter()
        .map(|node| {
            if node.is_assembly {
                return Vec::new();
            }
            model.label_shape(node.label).map(|s| faces_of(&s)).unwrap_or_default()
        })
        .collect();
    let (mut by_key, tris) = mesh_missing_faces(&face_lists, lin, ang);

    let mut nodes = Vec::with_capacity(explored.len());
    // Node index of the last node seen at each depth, for the parent link.
    let mut by_depth: Vec<usize> = Vec::new();
    for (index, (node, faces)) in explored.iter().zip(&face_lists).enumerate() {
        let parent = if node.depth == 0 {
            None
        } else {
            by_depth.get(node.depth - 1).copied()
        };
        by_depth.truncate(node.depth);
        by_depth.push(index);

        let mut out_faces = Vec::new();
        if !node.is_assembly {
            // `RWMesh_ShapeIterator`'s style dispatch for the node's Product
            // label (`RWMesh_ShapeIterator.cxx:70-134`).
            let settings = crate::xcaf_prs::collect_style_settings(model, node.ref_label);
            let face_styles = dispatch_face_styles(&settings, &node.style);
            for face in faces {
                let key = GeometryRegistry::shape_key(&face.0);
                let style = face_styles.get(&key).cloned().unwrap_or_else(|| node.style.clone());
                match face_source(face, &mut by_key, &tris) {
                    Some(source) => out_faces.push(face_node(face, source, style)),
                    None => out_faces.push(empty_face_node(face, style)),
                }
            }
        }

        nodes.push(crate::rwmesh::GltfCafNode {
            id: node.id.clone(),
            label: node.label,
            ref_label: node.ref_label,
            label_name: model.label_name(node.label).to_string(),
            product_name: model.label_name(node.ref_label).to_string(),
            entry: node.label.to_string(),
            is_assembly: node.is_assembly,
            parent,
            local_trsf: occt_core::gp::GpTrsf::identity(),
            named_data: model.label(node.label).and_then(|l| l.named_data.clone()),
            ref_named_data: model.label(node.ref_label).and_then(|l| l.named_data.clone()),
            faces: out_faces,
        });
    }
    nodes
}

/// `RWMesh_FaceIterator::IsEmpty` for a face the mesh has no triangulation for:
/// the style is still read by the glTF material passes, the geometry passes
/// skip it.
fn empty_face_node(face: &occt_topo::shape::Face, style: XcafPrsStyle) -> crate::rwmesh::ObjFaceNode {
    crate::rwmesh::ObjFaceNode {
        shape_key: GeometryRegistry::shape_key(&face.0),
        vertices: Vec::new(),
        triangles: Vec::new(),
        has_normals: false,
        has_tex_coords: false,
        normals: Vec::new(),
        tex_coords: Vec::new(),
        style,
    }
}
