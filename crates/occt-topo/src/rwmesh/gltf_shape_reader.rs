//! glTF read path, shape level: proxy faces/edges/vertices and document fill.
//!
//! Source: `RWGltf_GltfJsonParser::gltfParsePrimArray`
//! (`RWGltf_GltfJsonParser.cxx:1948-2035`) and
//! `RWGltf_CafReader::fillDocument` / `addShapeIntoDoc`
//! (`RWGltf_CafReader.cxx:435-688`).
//!
//! A glTF triangle primitive is *mesh data without a surface*: OCCT makes a
//! face whose triangulation is the primitive (`aBuilder.MakeFace(aFace,
//! aMeshData)`, `cxx:2002-2006`) and never attaches a `Geom_Surface`. Line
//! primitives become a polygon-backed edge and point primitives a compound of
//! vertices (`cxx:1972-1998`), both without a surface either. The
//! per-shape attributes (name + `XCAFPrs_Style` with the bound
//! `XCAFDoc_VisMaterial`) go into `myAttribMap` (`cxx:2011-2028`), which the
//! CAF reader later replays into the XDE document label by label.
//!
//! This crate has no XDE document; [`BRepModel`] stands in for it, keeping one
//! named shape per node and the material table with the per-shape links
//! (`XCAFDoc_VisMaterialTool::AddMaterial` / `SetShapeMaterial`). The label tree
//! `addShapeIntoDoc` builds lives in [`super::gltf_document`].
//!
//! The per-shape `extras` (`setShapeNamedData`, `cxx:658`) are kept on
//! [`crate::model::ModelShape::named_data`] as `TDataStd_NamedData`.

use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_core::poly::Triangulation;
use occt_core::poly::triangulation::Triangle;

use crate::builder::TopoBuilder;
use crate::model::BRepModel;
use crate::shape::{Compound, TopoShape};

use super::gltf::RwGltfPrimitiveMode;
use super::gltf_file::GltfReadOptions;
use super::scene::{MeshNode, MeshScene};

/// The `Poly_Triangulation` a glTF primitive turns into
/// (`RWGltf_GltfJsonParser.cxx:1948-1970`): the node's indexed triangle list,
/// with the per-node UV coordinates when the primitive carried TEXCOORD_0.
fn node_triangulation(node: &MeshNode) -> Triangulation {
    let tris = node
        .triangles
        .iter()
        .map(|&(a, b, c)| Triangle::new(a, b, c))
        .collect();
    let mut tri = Triangulation::new(node.vertices.clone(), tris);
    if let Some(uv) = &node.uv {
        if uv.len() == node.vertices.len() {
            tri = tri.with_uv(uv.iter().map(|p| (p.x(), p.y())).collect());
        }
    }
    // `RWGltf_TriangulationReader` stores the NORMAL accessor on the
    // triangulation (`setNodeNormal`, `RWGltf_TriangulationReader.cxx:330-364`),
    // which is what `Poly_Triangulation::HasNormals` / `Normal`
    // (`RWMesh_FaceIterator.cxx:59-67`) read back on export.
    if let Some(normals) = &node.normals {
        if normals.len() == node.vertices.len() {
            tri = tri.with_normals(normals.clone());
        }
    }
    tri
}

/// The shape `gltfParsePrimArray` builds for one primitive array
/// (`RWGltf_GltfJsonParser.cxx:1969-2013`):
///
/// * `Triangles` — a face whose only content is the node triangulation
///   (`cxx:2000-2007`);
/// * `Lines` — one edge over a `Poly_Polygon3D` built from the node indices
///   (`cxx:1985-1998`); the port uses [`TopoBuilder::make_edge_polygon`], the
///   same `BRep_Builder::MakeEdge(edge, polygon)` call;
/// * `Points` — a compound holding one vertex per node
///   (`BRep_Builder::MakeVertex(vertex, node, Precision::Confusion())`,
///   `cxx:1972-1983`).
///
/// This is the shape the reader binds attributes to and adds to the document
/// (`thePrimArrayShape`, `cxx:2031-2033`); [`node_proxy_shape`] wraps it in the
/// compound a scene node stands for.
pub fn prim_proxy_shape(node: &MeshNode) -> TopoShape {
    let b = TopoBuilder::new();
    match node.mode {
        RwGltfPrimitiveMode::Points => {
            let verts: Vec<TopoShape> = node
                .vertices
                .iter()
                .map(|p| {
                    // `cxx:1978-1981`: tolerance `Precision::Confusion()`.
                    b.make_vertex(p.clone(), occt_core::precision::CONFUSION).0
                })
                .collect();
            b.make_compound_of(&verts).0
        }
        RwGltfPrimitiveMode::Lines => {
            let nodes: Vec<GpPnt> = node
                .polyline
                .iter()
                .filter_map(|&i| node.vertices.get(i).cloned())
                .collect();
            b.make_edge_polygon(Arc::new(occt_core::poly::Polygon3D::new(nodes))).0
        }
        _ => b.make_face_from_triangulation(Arc::new(node_triangulation(node))).0,
    }
}

/// One node as a compound holding its proxy shape — the scene-level shape, and
/// what the flat [`BRepModel`] keeps per placed primitive array.
///
/// The node's vertices are already in world coordinates (the reader applies the
/// accumulated TRS while walking the graph), so the shape carries no location;
/// applying `node.transform` again would double-transform.
pub fn node_proxy_shape(node: &MeshNode) -> TopoShape {
    let b = TopoBuilder::new();
    let prim = prim_proxy_shape(node);
    let mut comp = Compound::new();
    b.add_compound(&mut comp, &prim);
    comp.0
}

/// One proxy shape per node of the scene, in node order.
pub fn scene_proxy_shapes(scene: &MeshScene) -> Vec<TopoShape> {
    scene.nodes.iter().map(node_proxy_shape).collect()
}

/// Read a glTF file into a [`BRepModel`]: the flat shape list (one proxy shape
/// per placed primitive array, named after the node) plus the document's
/// `XCAFDoc_ShapeTool` label tree.
///
/// The document side is `RWGltf_CafReader::fillDocument` (`cxx:435-467`) over the
/// per-shape attributes `RWGltf_GltfJsonParser::bindNamedShape` collected; see
/// [`super::gltf_document`] for what the label tree holds and which C++ steps it
/// does not reproduce.
pub fn read_gltf_model(path: &str) -> Result<BRepModel, String> {
    read_gltf_model_with_options(path, &GltfReadOptions::occt_defaults())
}

/// [`read_gltf_model`] with an explicit converter (`RWGltf_CafReader`'s
/// `SetCoordinateSystemConverter`, `RWGltf_CafReader.cxx:297`).
pub fn read_gltf_model_with_options(
    path: &str,
    opts: &GltfReadOptions,
) -> Result<BRepModel, String> {
    super::gltf_file::read_gltf_document_with_options(path, opts)
}

/// The flat document fill of [`read_gltf_scene`]'s result, reusable for scenes
/// that were built by other means: no label tree, one free shape per scene node.
pub fn scene_to_model(scene: &MeshScene) -> BRepModel {
    let mut model = BRepModel::new();
    // `myMaterials` is keyed by the material handle; clones of one
    // `XcafDocVisMaterial` share `id()`, so the table dedups exactly like
    // `XCAFDoc_VisMaterialTool::AddMaterial` does with `IsEqual`.
    let mut material_index: std::collections::HashMap<u64, usize> = std::collections::HashMap::new();

    for node in &scene.nodes {
        // OCCT can hold duplicate `TDataStd_Name` labels; `BRepModel` keys shapes
        // by name, so a repeated node name gets a positional suffix. This is a
        // container constraint, not a naming rule of the reader.
        let mut name = node.name.clone();
        let mut n = 2;
        while model.find(&name).is_some() {
            name = format!("{}_{n}", node.name);
            n += 1;
        }
        let n = model.len();
        model.add(&name, node_proxy_shape(node));
        if let Some(nd) = &node.named_data {
            model.shapes[n].named_data = Some(nd.clone());
        }
        if let Some(vis) = &node.vis_material {
            let idx = *material_index.entry(vis.id()).or_insert_with(|| {
                model.add_material(vis.clone(), vis.raw_name().unwrap_or(""))
            });
            model.set_shape_material(&name, Some(idx));
        }
    }
    model
}
