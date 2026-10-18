//! Wavefront OBJ + MTL write side. Source: `RWObj` (`RWObj_CafWriter`,
//! `RWObj_ObjWriterContext`, `RWObj_ObjMaterialMap`) and `RWMesh_MaterialMap`.
//!
//! Faithful port of the OCCT export pipeline used by `DRAWEXE`'s `WriteObj`
//! command:
//!
//! * [`MeshMaterialMap`] — `RWMesh_MaterialMap` (`RWMesh_MaterialMap.cxx`):
//!   material registration, `mat_<n>` key / name generation, and the texture
//!   folder machinery (`CreateTextureFolder` / `CopyTexture` / `copyFileTo`).
//! * [`ObjMaterialMap`] — `RWObj_ObjMaterialMap`
//!   (`RWObj_ObjMaterialMap.cxx`): opens the `.mtl` lazily on the first
//!   material and writes `newmtl` / `Ka` / `Kd` / `Ks` / `Ns` / `Tr` / `map_Kd`.
//! * [`ObjWriterContext`] — `RWObj_ObjWriterContext`
//!   (`RWObj_ObjWriterContext.cxx`): the header (with `mtllib`), `usemtl`,
//!   `g`, `v`, `vn`, `vt` and `f` records.
//! * [`perform`] — `RWObj_CafWriter::Perform` (`RWObj_CafWriter.cxx:82-334`):
//!   the `toCreateMatFile` decision (`cxx:242-257`), the `.mtl` file naming
//!   (`cxx:150-157`) and the per-node `writeShape` loop.
//!
//! Not ported (no counterpart in this crate, documented rather than faked):
//! the `XCAFDoc_VisMaterial` common / PBR material properties
//! (`RWObj_ObjMaterialMap.cxx:78-92`). That arm is dead on the STEP -> OBJ path in
//! OCCT too (see the `xcaf_prs` module header); only the material's base-colour
//! texture reaches the style. Also not ported: the `WriteHeader` file-info map
//! (`RWObj_ObjWriterContext.cxx:112-135`).
//!
//! The `v` / `vn` / `vt` / `f` index shifts follow `WriteTriangle`
//! (`RWObj_ObjWriterContext.cxx:154-206`) and `FlushFace` (`cxx:303-312`)
//! verbatim: `myElemPosFirst` / `myElemNormFirst` / `myElemUVFirst` all advance
//! by the face's node count, so the three streams carry the same index.

use occt_core::gp::{GpPnt, GpPnt2d};
use occt_core::poly::triangulation::Triangle;

use crate::xcaf_prs::XcafPrsStyle;

mod material_map;
mod writer_context;
pub use material_map::*;
pub use writer_context::*;

// ---------------------------------------------------------------------------
// RWObj_CafWriter::Perform
// ---------------------------------------------------------------------------

/// One `RWMesh_FaceIterator` face: the face's own node block plus the style the
/// shape iterator bound to it (`RWMesh_ShapeIterator::initShape`, `cxx:139-162`).
pub struct ObjFaceNode {
    /// `RWMesh_FaceIterator`'s `TopoDS_Face`, as the `TShape` key
    /// (`GeometryRegistry::shape_key`). `RWGltf_StyledShape::Shape` compares
    /// faces by this identity (`RWGltf_CafWriter.hxx:491-495`); the OBJ writer
    /// only needs it to look the face up in the mesh.
    pub shape_key: usize,
    /// The face's `Poly_Triangulation` nodes (`RWMesh_FaceIterator::NbNodes`).
    pub vertices: Vec<GpPnt>,
    /// Triangles into [`Self::vertices`], 0-based (`Poly_Triangulation` indices
    /// minus one, `RWObj_CafWriter.cxx:407-408`).
    pub triangles: Vec<Triangle>,
    /// `RWMesh_FaceIterator::HasNormals` (`hxx:109-111` + `cxx:106-124`): the
    /// triangulation carries normals, or it has UV nodes and the face's surface
    /// is non-null (in which case [`Self::normals`] came from the surface).
    pub has_normals: bool,
    /// `RWMesh_FaceIterator::HasTexCoords` (`hxx:112`):
    /// `Poly_Triangulation::HasUVNodes`.
    pub has_tex_coords: bool,
    /// `RWMesh_FaceIterator::NormalTransformed(i)` per node (`hxx:113-126`),
    /// i.e. the surface normal at the node's UV, transformed by the face
    /// location and flipped for a `TopAbs_REVERSED` face. Only meaningful when
    /// [`Self::has_normals`] is set; `gp::DZ()` otherwise.
    pub normals: Vec<GpPnt>,
    /// `RWMesh_FaceIterator::NodeTexCoord(i)` per node (`hxx:139-142`), the
    /// `Poly_Triangulation` UV node.
    pub tex_coords: Vec<GpPnt2d>,
    /// `aFaceIter.FaceStyle()` (`cxx:302`) — the per-face style.
    pub style: XcafPrsStyle,
}

/// One document node for [`perform`]: one `RWObj_CafWriter::writeShape` call,
/// i.e. one `XCAFPrs_DocumentExplorer` leaf node and its faces
/// (`RWObj_CafWriter.cxx:192-221`).
pub struct ObjExportNode {
    /// `readNameAttribute(aDocNode.RefLabel)` plus the parent prefix
    /// (`cxx:197-207`) — the `g` group name.
    pub name: String,
    /// The node's faces in `RWMesh_FaceIterator` order.
    pub faces: Vec<ObjFaceNode>,
}

/// `RWObj_CafWriter::toSkipFaceMesh` (`RWObj_CafWriter.cxx:74-77`).
fn to_skip_face_mesh(face: &ObjFaceNode) -> bool {
    face.vertices.is_empty() || face.triangles.is_empty()
}

/// `RWObj_CafWriter::addFaceInfo` (`RWObj_CafWriter.cxx:235-257`).
///
/// The progress-step terms (`cxx:242-252`) only size the progress sentry, which
/// this port does not model, so only `theNbNodes` / `theNbElems` /
/// `theToCreateMatFile` are accumulated.
fn add_face_info(face: &ObjFaceNode, nb_nodes: &mut i32, nb_elems: &mut i32, to_create_mat_file: &mut bool) {
    *nb_nodes += face.vertices.len() as i32;
    *nb_elems += face.triangles.len() as i32;
    // `theToCreateMatFile || theFace.HasFaceColor() || (!BaseColorTexture()
    // .IsNull() && HasTexCoords())` (`cxx:254-256`).
    *to_create_mat_file = *to_create_mat_file
        || face.style.has_color()
        || (face.style.base_color_texture().is_some() && face.has_tex_coords);
}

/// `RWObj_CafWriter::writeShape` (`RWObj_CafWriter.cxx:262-334`), restricted to
/// the triangle-only path this crate meshes into.
fn write_shape(writer: &mut ObjWriterContext, mat_mgr: &mut ObjMaterialMap, node: &ObjExportNode) -> bool {
    let mut to_create_group = true;
    // `for (RWMesh_FaceIterator aFaceIter(theLabel, theParentTrsf, true,
    // theParentStyle); ...)` (`cxx:270-332`).
    for face in &node.faces {
        if to_skip_face_mesh(face) {
            continue;
        }

        writer.nb_faces += 1;
        {
            // `const bool hasNormals = aFaceIter.HasNormals(); const bool
            // hasTexCoords = aFaceIter.HasTexCoords(); if (theWriter.NbFaces !=
            // 1) { toCreateGroup = toCreateGroup || hasNormals !=
            // theWriter.HasNormals() || hasTexCoords !=
            // theWriter.HasTexCoords(); } theWriter.SetNormals(hasNormals);
            // theWriter.SetTexCoords(hasTexCoords);` (`cxx:280-291`). The
            // comparison reads the *previous* face's flags, so a node whose
            // faces alternate between having and not having normals splits into
            // several `g` records.
            let has_normals = face.has_normals;
            let has_tex_coords = face.has_tex_coords;
            if writer.nb_faces != 1 {
                to_create_group = to_create_group
                    || has_normals != writer.has_normals()
                    || has_tex_coords != writer.has_tex_coords();
            }
            writer.set_normals(has_normals);
            writer.set_tex_coords(has_tex_coords);
        }

        if to_create_group && !writer.write_group(&node.name) {
            return false;
        }
        to_create_group = false;

        // `TCollection_AsciiString aMatName; if (aFaceIter.HasFaceColor() ||
        // !aFaceIter.FaceStyle().BaseColorTexture().IsNull()) { aMatName =
        // theMatMgr.AddMaterial(aFaceIter.FaceStyle()); }` (`cxx:300-303`).
        let mut mat_name = String::new();
        if face.style.has_color() || face.style.base_color_texture().is_some() {
            mat_name = mat_mgr.add_material(&face.style);
        }
        if mat_name != writer.active_material() {
            // `theWriter.WriteActiveMaterial(aMatName)` (`cxx:304-307`).
            writer.write_active_material(&mat_name);
        }

        // write nodes (`cxx:311-314`, `writePositions` `cxx:337-353`).
        for v in &face.vertices {
            writer.write_vertex([v.x() as f32, v.y() as f32, v.z() as f32]);
        }

        // write normals (`cxx:317-320`, `writeNormals` `cxx:357-374`).
        if writer.has_normals() {
            for n in &face.normals {
                writer.write_normal([n.x() as f32, n.y() as f32, n.z() as f32]);
            }
        }

        // write texcoords (`cxx:322-326`, `writeTextCoords` `cxx:378-392`).
        if writer.has_tex_coords() {
            for uv in &face.tex_coords {
                writer.write_tex_coord([uv.x() as f32, uv.y() as f32]);
            }
        }

        // write indices (`cxx:328-331`, `writeIndices` `cxx:396-414`).
        // `WriteTriangle(Vec3(aTri(1),aTri(2),aTri(3)) - Vec3(anElemLower))`
        // (`cxx:407-408`): the `Poly_Triangle` values are 1-based and
        // `anElemLower == 1` for the port's per-face mesh, so the payload is the
        // 0-based index that `WriteTriangle` then shifts by `myElemPosFirst`.
        for t in &face.triangles {
            writer.write_triangle([t.n0 as i32, t.n1 as i32, t.n2 as i32]);
        }
        writer.flush_face(face.vertices.len() as i32);
    }
    true
}

/// `RWObj_CafWriter::Perform` (`RWObj_CafWriter.cxx:82-230`) over an in-memory
/// node list. `obj_path` supplies the base name for the `.mtl` sidecar.
///
/// Returns `(obj_text, mat_text)`; `mat_text` is empty when no node carries a
/// color (`toCreateMatFile == false`, `cxx:150-157`).
pub fn write_obj(obj_path: &str, nodes: &[ObjExportNode]) -> (String, String) {
    // `OSD_Path::FolderAndFileFromPath` / `FileNameAndExtension`
    // (`cxx:100-103`).
    let (folder, short_file_name_base) = split_path(obj_path);

    // First pass: `addFaceInfo` over every face of every node
    // (`cxx:119-144`).
    let mut nb_nodes_all = 0i32;
    let mut nb_elems_all = 0i32;
    let mut to_create_mat_file = false;
    for node in nodes {
        for face in &node.faces {
            if to_skip_face_mesh(face) {
                continue;
            }
            add_face_info(face, &mut nb_nodes_all, &mut nb_elems_all, &mut to_create_mat_file);
        }
    }

    // `TCollection_AsciiString aMatFileNameShort = aShortFileNameBase + ".mtl"`
    // (`cxx:150-157`); cleared when no material file is needed.
    let mut mat_file_name_short = format!("{short_file_name_base}.mtl");
    let mat_file_name_full = if !folder.is_empty() {
        format!("{folder}{mat_file_name_short}")
    } else {
        mat_file_name_short.clone()
    };
    if !to_create_mat_file {
        mat_file_name_short.clear();
    }

    let mut writer = ObjWriterContext::new(obj_path);
    let mut mat_mgr = ObjMaterialMap::new(&mat_file_name_full);
    // `aMatMgr.SetDefaultStyle(myDefaultStyle)` (`cxx:161`) — the port's default
    // style is empty, so this is a no-op.
    if !writer.write_header(nb_nodes_all, nb_elems_all, &mat_file_name_short, &[]) {
        return (writer.content().to_string(), mat_mgr.content().to_string());
    }

    // Second pass: `writeShape` per node (`cxx:184-221`). A node whose faces are
    // all empty writes no group but still runs the loop, exactly as OCCT's
    // `writeShape` skips the faces inside its own loop.
    for node in nodes {
        if !write_shape(&mut writer, &mut mat_mgr, node) {
            break;
        }
    }

    let _closed = writer.close();
    mat_mgr.close();
    (writer.content().to_string(), mat_mgr.content().to_string())
}

/// `OSD_Path::FolderAndFileFromPath` + `FileNameAndExtension` -> the folder
/// (with a trailing separator when non-empty) and the file stem.
fn split_path(path: &str) -> (String, String) {
    let normalized = path.replace('\\', "/");
    let (folder_part, file) = match normalized.rfind('/') {
        Some(i) => (&normalized[..i], &normalized[i + 1..]),
        None => ("", normalized.as_str()),
    };
    let stem = match file.rfind('.') {
        Some(i) => &file[..i],
        None => file,
    };
    let folder = if folder_part.is_empty() {
        String::new()
    } else {
        format!("{folder_part}/")
    };
    (folder, stem.to_string())
}

/// `OSD_Path::FolderAndFileFromPath`-based `.mtl` path for an OBJ path, exposed
/// so the caller can create the sidecar directory.
pub fn material_file_path(obj_path: &str) -> String {
    let (folder, stem) = split_path(obj_path);
    format!("{folder}{stem}.mtl")
}
