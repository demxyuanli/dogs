//! Module 2: brep_exchange — export BRep shapes to OBJ / STL / PLY.
//! Source: `RWObj`, `RWStl`, `RWPLY` (simplified OCCT writers).

use std::collections::{HashMap, VecDeque};

use crate::mesh::ShapeMesh;
use crate::model::BRepModel;
use crate::shape::TopoShape;
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::faces_of;
use crate::xcaf_prs::{dispatch_face_styles, XcafPrsStyle};
use occt_core::gp::GpPnt;
use occt_core::io::obj::{ObjFace, ObjMesh};
use occt_core::io::ply::PlyMesh;
use occt_core::io::stl::StlMesh;
use occt_core::poly::triangulation::Triangle;

/// Convert a `ShapeMesh` to an OBJ mesh (each triangle becomes an OBJ face).
fn shape_mesh_to_obj(mesh: &ShapeMesh) -> ObjMesh {
    ObjMesh {
        vertices: mesh.vertices.clone(),
        texcoords: Vec::new(),
        normals: Vec::new(),
        faces: mesh.triangles.iter().map(|t| ObjFace {
            v: vec![t.n0 as i32, t.n1 as i32, t.n2 as i32],
            vt: None,
            vn: None,
        }).collect(),
    }
}

/// Convert a `ShapeMesh` to an STL mesh (independent per-triangle vertices).
fn shape_mesh_to_stl(mesh: &ShapeMesh) -> StlMesh {
    StlMesh {
        triangles: mesh.triangles.iter().map(|t| [
            mesh.vertices[t.n0],
            mesh.vertices[t.n1],
            mesh.vertices[t.n2],
        ]).collect(),
        normals: Vec::new(),
    }
}

/// Convert a `ShapeMesh` to a PLY mesh.
fn shape_mesh_to_ply(mesh: &ShapeMesh) -> PlyMesh {
    PlyMesh {
        vertices: mesh.vertices.clone(),
        faces: mesh.triangles.iter().map(|t| vec![t.n0, t.n1, t.n2]).collect(),
    }
}

/// Mesh a shape for export using the OCCT `BRepMesh_IncrementalMesh` pipeline:
/// the Delaunay path (`meshing::incremental_mesh_to_shape_mesh`) shares each
/// edge's 3D polyline between its adjacent faces, so the output is watertight on
/// shared edges and periodic seams. A failed tessellation yields no triangles:
/// OCCT has no second mesher (audit A17/A18/A19).
///
/// ponytail: B-spline faces whose STEP `SURFACE_CURVE` pcurves are dropped
/// (`step.rs`) get a wrong boundary, so the Delaunay bbox drifts past `EXACT_TOL`
/// on `Shape.step`/`Shape-2.step` — a STEP-import gap, not a BRepMesh one.
/// `Prs3d::GetDeflection(shape, drawer)` (`Prs3d.hxx:82-103`):
/// `BRepBndLib::Add(shape, box, false)` then
/// `maxComp(diag) * DeviationCoefficient * 4`, coefficient default 0.001.
pub fn prs3d_get_deflection(shape: &TopoShape, maximal_chordial_deviation: f64) -> f64 {
    const DEVIATION_COEFFICIENT: f64 = 0.001;
    const CONFUSION: f64 = 1e-7;
    let b = crate::brep_bnd_lib::shape_bnd_box(shape);
    if b.is_void() {
        return maximal_chordial_deviation;
    }
    let b = if b.is_open() {
        if !b.has_finite_part() {
            return maximal_chordial_deviation;
        }
        b.finite_part()
    } else {
        b
    };
    let mn = b.corner_min();
    let mx = b.corner_max();
    let max_comp = (mx.x() - mn.x())
        .max(mx.y() - mn.y())
        .max(mx.z() - mn.z());
    (max_comp * DEVIATION_COEFFICIENT * 4.0).max(CONFUSION)
}

fn export_mesh(shape: &TopoShape, maximal_chordial_deviation: f64) -> ShapeMesh {
    // `RWObj_CafWriter` writes the vis triangulation
    // (`StdPrs_ToolTriangulatedShape::Tessellate` →
    // `BRepMesh_DiscretFactory::Discret(shape, GetDeflection, DeviationAngle)`).
    // Drawer defaults: 20 deg (`Prs3d_Drawer.cxx:96`) and relative
    // `GetDeflection` (`Prs3d.hxx:71`). The argument is the drawer's
    // `MaximalChordialDeviation` — it is **not** a linear deflection and only
    // matters for a void/unbounded bounding box (`Prs3d.hxx:82-103`).
    //
    // A failed tessellation yields an **empty** mesh: OCCT's OBJ writer reads
    // the triangulation already stored on the shape
    // (`RWMesh_FaceIterator.cxx:87`) and **skips** a face without one (`:89`) —
    // it never switches to another mesher. The two silent fallbacks that used to
    // sit here (legacy quadtree `brepmesh::incremental_mesh`, UV-grid
    // `shape_mesh::mesh_shape`) were non-OCCT algorithms and are removed
    // (audit A17); `brepmesh` itself is deleted as of A18/A19.
    const PRS3D_DEV_ANGLE: f64 = 20.0 * std::f64::consts::PI / 180.0;
    let lin = prs3d_get_deflection(shape, maximal_chordial_deviation);
    crate::meshing::incremental_mesh::IncrementalMesh::from_deflection(
        shape,
        lin,
        false,
        PRS3D_DEV_ANGLE,
    )
    .mesh()
    .cloned()
    .unwrap_or_else(|| ShapeMesh {
        vertices: Vec::new(),
        triangles: Vec::new(),
        source_shape: crate::abs::ShapeType::Shape,
    })
}

/// Export a shape to Wavefront OBJ text.
///
/// Vertices stay per face (`RWObj_CafWriter` / `RWMesh_FaceIterator`):
/// a box is 24 nodes / 12 triangles, not a globally welded 8-node mesh.
pub fn brep_to_obj(shape: &TopoShape, maximal_chordial_deviation: f64) -> String {
    let mesh = export_mesh(shape, maximal_chordial_deviation);
    occt_core::io::obj::obj_mesh_to_string(&shape_mesh_to_obj(&mesh))
}

/// Export a shape to ASCII STL text.
pub fn brep_to_stl_ascii(shape: &TopoShape, maximal_chordial_deviation: f64) -> String {
    let mesh = export_mesh(shape, maximal_chordial_deviation);
    occt_core::io::stl::write_ascii_stl(&shape_mesh_to_stl(&mesh))
}

/// Export a shape to binary STL bytes (vertices welded first).
pub fn brep_to_stl_binary(shape: &TopoShape, maximal_chordial_deviation: f64) -> Vec<u8> {
    let mut mesh = export_mesh(shape, maximal_chordial_deviation);
    crate::shape_mesh::weld_vertices(&mut mesh, 1e-9);
    occt_core::io::stl::write_binary_stl(&shape_mesh_to_stl(&mesh))
}

/// Export a shape to ASCII PLY text.
///
/// `RWPly_PlyWriterContext` writes each **face's own node block** with a running
/// `myVertOffset` (`RWPly_PlyWriterContext.cxx:276`, `:300`) — it never merges
/// nodes across faces, so this path must **not** weld: a box is 24 vertex
/// records (6 faces × 4), exactly like the OBJ path ([`brep_to_obj`]). The
/// previous body called `weld_vertices(1e-9)` first, producing a shared 8-node
/// mesh that OCCT never writes (audit A26).
pub fn brep_to_ply(shape: &TopoShape, maximal_chordial_deviation: f64) -> String {
    let mesh = export_mesh(shape, maximal_chordial_deviation);
    occt_core::io::ply::ply_mesh_to_string(&shape_mesh_to_ply(&mesh))
}

/// Write a shape to an OBJ file (per-face nodes, same as [`brep_to_obj`]).
pub fn brep_write_obj(path: &str, shape: &TopoShape, maximal_chordial_deviation: f64) -> std::io::Result<()> {
    let mesh = export_mesh(shape, maximal_chordial_deviation);
    occt_core::io::obj::write_obj_file(path, &shape_mesh_to_obj(&mesh))
}

// ---------------------------------------------------------------------------
// RWObj_CafWriter: OBJ + MTL for a document (`BRepModel`)
// ---------------------------------------------------------------------------

/// One face's mesh data as `RWMesh_FaceIterator::Next` obtains it
/// (`RWMesh_FaceIterator.cxx:83-101`): the face's own
/// `BRep_Tool::Triangulation`, else the triangulation the export driver meshed
/// into the face.
struct FaceSource {
    vertices: Vec<GpPnt>,
    triangles: Vec<Triangle>,
    uv: Vec<occt_core::gp::GpPnt2d>,
    /// `Poly_Triangulation::Normal` (`cxx:59-67`): the stored per-node normals.
    normals: Option<Vec<GpPnt>>,
}

/// `BRep_Tool::Triangulation(theFace, theLocation)` first (`cxx:87`), then the
/// mesh-pass entry for that face. `None` means the face has neither, which
/// `RWMesh_FaceIterator::Next` skips together with a zero-triangle
/// triangulation (`cxx:88-93`).
fn face_source(
    face: &crate::shape::Face,
    by_key: &mut HashMap<usize, VecDeque<usize>>,
    tris: &[crate::meshing::triangulator::FaceTriangulation],
) -> Option<FaceSource> {
    if let Some(tri) = crate::brep_tool::BRepTool::face_triangulation(face) {
        if tri.triangles.is_empty() {
            return None;
        }
        return Some(FaceSource {
            vertices: tri.nodes.clone(),
            triangles: tri.triangles.clone(),
            uv: tri
                .uv_nodes
                .clone()
                .unwrap_or_default()
                .into_iter()
                .map(|(u, v)| occt_core::gp::GpPnt2d::new(u, v))
                .collect(),
            normals: tri.normals.clone(),
        });
    }
    let key = GeometryRegistry::shape_key(&face.0);
    let entry = by_key.get_mut(&key)?.pop_front()?;
    let ft = &tris[entry];
    Some(FaceSource {
        vertices: ft.vertices.clone(),
        triangles: ft.triangles.clone(),
        uv: ft.uv.clone(),
        normals: None,
    })
}

/// `RWMesh_FaceIterator` for one face (`hxx:85-152`): orientation-applied
/// triangles, `HasNormals` / `HasTexCoords` and the transformed normals.
///
/// `BRepLProp_SLProps(1, 1e-12)` is the iterator's surface-properties tool
/// (`RWMesh_FaceIterator.cxx:35,47`); its resolution is the CSLib tolerance
/// `normal()` hands to the normal estimate.
const SL_PROPS_LIN_TOL: f64 = 1e-12;

fn face_node(face: &crate::shape::Face, source: FaceSource, style: XcafPrsStyle) -> crate::rwmesh::ObjFaceNode {
    // `RWMesh_FaceIterator::TriangleOriented` (`hxx:97-105`): a
    // `TopAbs_REVERSED` face occurrence swaps the last two nodes. The
    // `myIsMirrored` term needs a negative-determinant occurrence transform,
    // which this crate's shapes do not carry.
    let reversed = face.0.orientation().is_reversed();
    let triangles: Vec<Triangle> = source
        .triangles
        .iter()
        .map(|t| if reversed { Triangle::new(t.n0, t.n2, t.n1) } else { *t })
        .collect();

    // `RWMesh_FaceIterator::initFace` (`cxx:104-125`): `HasNormals()` is set
    // when the triangulation carries normals, or when it has UV nodes and the
    // face's surface (forward-oriented, unlocated, as in `cxx:115-121`) is
    // non-null. `HasTexCoords()` is `Poly_Triangulation::HasUVNodes`
    // (`hxx:112`).
    let has_tex_coords = !source.uv.is_empty() && source.uv.len() == source.vertices.len();
    let surface = if has_tex_coords {
        crate::brep_tool::BRepTool::face_surface(face)
    } else {
        None
    };
    let has_normals = source.normals.is_some() || surface.is_some();

    // `RWMesh_FaceIterator::NormalTransformed(i)` (`hxx:113-126`): the stored
    // normal, transformed by the face location (identity here) and reversed for
    // a `TopAbs_REVERSED` face. Without stored normals the value comes from
    // `normal()` (`cxx:59-77`): `BRepLProp_SLProps::Normal` at the node's UV,
    // and `gp::DZ()` when `IsNormalDefined()` is false (a singular point such as
    // a cone apex). `BRepLProp_SLProps(1, 1e-12)` is the iterator's tool
    // (`cxx:35,47`), so the CSLib tolerance is `1e-12`.
    let normals: Vec<GpPnt> = if let Some(mut ns) = source.normals {
        if reversed {
            for n in &mut ns {
                *n = GpPnt::new(-n.x(), -n.y(), -n.z());
            }
        }
        ns
    } else if has_normals {
        let mut ns = Vec::with_capacity(source.vertices.len());
        for uv in &source.uv {
            // `gp_Dir aNormal(gp::DZ());`
            let mut dir = GpPnt::new(0.0, 0.0, 1.0);
            if let Some(surf) = surface.as_deref() {
                let (_, du, dv) = surf.d1(uv.x(), uv.y());
                if let (occt_core::cslib::CSLibDerivativeStatus::Done, Some(n)) =
                    occt_core::cslib::normal_d1(&du, &dv, SL_PROPS_LIN_TOL)
                {
                    dir = GpPnt::new(n.x(), n.y(), n.z());
                }
            }
            if reversed {
                dir = GpPnt::new(-dir.x(), -dir.y(), -dir.z());
            }
            ns.push(dir);
        }
        ns
    } else {
        Vec::new()
    };

    crate::rwmesh::ObjFaceNode {
        shape_key: GeometryRegistry::shape_key(&face.0),
        vertices: source.vertices,
        triangles,
        has_normals,
        has_tex_coords,
        normals,
        tex_coords: source.uv,
        style,
    }
}

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
    let inc = crate::meshing::incremental_mesh::IncrementalMesh::from_deflection(
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
fn mesh_missing_faces(
    leaf_faces: &[Vec<crate::shape::Face>],
    lin: f64,
    ang: f64,
) -> (
    HashMap<usize, VecDeque<usize>>,
    Vec<crate::meshing::triangulator::FaceTriangulation>,
) {
    let mut missing: Vec<TopoShape> = Vec::new();
    for faces in leaf_faces {
        for face in faces {
            if crate::brep_tool::BRepTool::face_triangulation(face).is_some() {
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
    let mut tris: Vec<crate::meshing::triangulator::FaceTriangulation> = Vec::new();
    if !missing.is_empty() {
        let input = if missing.len() == 1 {
            missing[0].clone()
        } else {
            crate::builder::TopoBuilder::new().make_compound_of(&missing).0
        };
        let inc = crate::meshing::incremental_mesh::IncrementalMesh::from_deflection(&input, lin, false, ang);
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
    let leaf_faces: Vec<Vec<crate::shape::Face>> = leaves
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
    let face_lists: Vec<Vec<crate::shape::Face>> = explored
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
fn empty_face_node(face: &crate::shape::Face, style: XcafPrsStyle) -> crate::rwmesh::ObjFaceNode {
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
            Some(crate::builder::TopoBuilder::new().make_compound_of(&parts).0)
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

/// Write a shape to a (binary) STL file.
pub fn brep_write_stl(path: &str, shape: &TopoShape, maximal_chordial_deviation: f64) -> std::io::Result<()> {
    std::fs::write(path, brep_to_stl_binary(shape, maximal_chordial_deviation))
}

/// Write a shape to an ASCII PLY file.
pub fn brep_write_ply(path: &str, shape: &TopoShape, maximal_chordial_deviation: f64) -> std::io::Result<()> {
    std::fs::write(path, brep_to_ply(shape, maximal_chordial_deviation))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep_tool::BRepTool;
    use crate::builder::TopoBuilder;
    use crate::shape::Vertex;
    use occt_core::gp::{GpAx1, GpAx3, GpDir, GpPln, GpPnt};

    /// Build a unit cube [0,1]^3 as a Solid with real geometry.
    fn unit_box() -> TopoShape {
        let b = TopoBuilder::new();
        let mk = |x: f64, y: f64, z: f64| b.make_vertex(GpPnt::new(x, y, z), 0.0);
        let p0 = mk(0.,0.,0.); let p1 = mk(1.,0.,0.); let p2 = mk(1.,1.,0.); let p3 = mk(0.,1.,0.);
        let p4 = mk(0.,0.,1.); let p5 = mk(1.,0.,1.); let p6 = mk(1.,1.,1.); let p7 = mk(0.,1.,1.);

        let face = |c: [&Vertex; 4], normal: GpDir| {
            let e = [
                b.make_edge_segment(&BRepTool::vertex_point(c[0]), &BRepTool::vertex_point(c[1])),
                b.make_edge_segment(&BRepTool::vertex_point(c[1]), &BRepTool::vertex_point(c[2])),
                b.make_edge_segment(&BRepTool::vertex_point(c[2]), &BRepTool::vertex_point(c[3])),
                b.make_edge_segment(&BRepTool::vertex_point(c[3]), &BRepTool::vertex_point(c[0])),
            ];
            let wire = b.make_wire(&e);
            let pln = GpPln::new(GpAx3::from_ax1(&GpAx1::new(BRepTool::vertex_point(c[0]), normal)));
            let mut f = b.make_face_plane(&pln);
            b.add_wire(&mut f, &wire);
            f
        };

        let nz = GpDir::new(0.,0.,-1.).unwrap(); let pz = GpDir::new(0.,0.,1.).unwrap();
        let nx = GpDir::new(-1.,0.,0.).unwrap(); let px = GpDir::new(1.,0.,0.).unwrap();
        let ny = GpDir::new(0.,-1.,0.).unwrap(); let py = GpDir::new(0.,1.,0.).unwrap();

        let faces = vec![
            face([&p0, &p1, &p2, &p3], nz),
            face([&p4, &p5, &p6, &p7], pz),
            face([&p0, &p3, &p7, &p4], nx),
            face([&p1, &p5, &p6, &p2], px),
            face([&p0, &p1, &p5, &p4], ny),
            face([&p3, &p2, &p6, &p7], py),
        ];
        let shell = b.make_shell(&faces);
        let solid = b.make_solid(&[shell]);
        solid.into()
    }

    #[test]
    fn obj_export_has_vertices_and_faces() {
        let shape = unit_box();
        let obj = brep_to_obj(&shape, 0.1);
        assert!(obj.contains("v 0.000000000 0.000000000 0.000000000"), "origin corner present");
        assert!(obj.lines().filter(|l| l.starts_with("f ")).count() > 0);
    }

    #[test]
    fn stl_ascii_export() {
        let shape = unit_box();
        let stl = brep_to_stl_ascii(&shape, 0.1);
        assert!(stl.contains("solid"));
        assert!(stl.contains("endsolid"));
    }

    #[test]
    fn stl_binary_and_ply_export() {
        let shape = unit_box();
        let bytes = brep_to_stl_binary(&shape, 0.1);
        assert!(bytes.len() >= 84);
        assert_eq!((bytes.len() - 84) % 50, 0);
        let ply = brep_to_ply(&shape, 0.1);
        assert!(ply.starts_with("ply"));
    }

    #[test]
    fn write_obj_file_roundtrip() {
        let shape = unit_box();
        let path = { let d = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("target").join("test-out"); let _ = std::fs::create_dir_all(&d); d }.join(format!("occt_brep_obj_{}.obj", std::process::id()));
        let p = path.to_str().unwrap().to_string();
        brep_write_obj(&p, &shape, 0.1).unwrap();
        let m = occt_core::io::obj::read_obj_file(&p).unwrap();
        assert!(m.vertices.len() > 0);
        assert!(m.faces.len() > 0);
        let _ = std::fs::remove_file(&p);
    }
}
