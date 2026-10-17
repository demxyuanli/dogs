//! glTF 2.0 / GLB write side. Source: `RWGltf_CafWriter`
//! (`RWGltf_CafWriter.cxx`), with `RWMesh::FormatName`
//! (`RWMesh.cxx:32-96`) and `RWGltf_GltfMaterialMap`
//! (`RWGltf_GltfMaterialMap.cxx`).
//!
//! Ported pipeline, in the C++'s order:
//!
//! 1. [`GltfCafWriter::write_bin_data`] — `writeBinData` (`cxx:811-1207`):
//!    collect one `RWGltf_GltfFace` per (face, style), then run the four array
//!    passes (`Position`, `Normal`, `TCoord0`, `Indices`) over that list,
//!    appending little-endian floats/uints to the binary body and recording one
//!    accessor per pass per face. A shape that appears in more than one document
//!    node is written once and its accessor copied (`aWrittenPrimData`,
//!    `cxx:971-998`).
//! 2. [`GltfCafWriter::write_json`] — `writeJson` (`cxx:1210-1426`): the
//!    `accessors` … `textures` members, then the GLB container when
//!    `is_binary` (12-byte header, JSON chunk padded with spaces to 4 bytes,
//!    BIN chunk).
//!
//! Notes carried over from the C++:
//!
//! * The *node* list is the document's labels in pre-order (the
//!   `XCAFPrs_DocumentExplorerFlags_None` walk `writeNodes` uses,
//!   `cxx:2287-2308`); the *mesh* list is the same list restricted to
//!   non-assembly nodes that carry geometry, in the same order (the
//!   `OnlyLeafNodes` walk, `cxx:1254-1302`). The two orders agree because both
//!   are depth-first pre-order.
//! * `myBinDataMap` is an `NCollection_IndexedDataMap` (`cxx:511-515`), i.e.
//!   insertion ordered: the binary layout, the accessor numbering and the
//!   `meshes` order all follow document order.
//! * Positions and normals go through `myCSTrsf` (`RWMesh_CoordinateSystemConverter`),
//!   which the writer's constructor sets to output unit 1.0 (metres) and the
//!   glTF Y-up system (`cxx:255-258`); `Perform` then feeds the document's
//!   length unit in (`cxx:645-649`).
//!
//! **UNPORTED** (no Rust counterpart, documented instead of faked):
//!
//! * `toMergeFaces` / `toSplitIndices16` (`cxx:670-753`, `cxx:2204-2246`) —
//!   merging several faces of one part into one primitive. The option would
//!   change both the `myBinDataMap` grouping and the per-face iterator
//!   offsets; the default (`myToMergeFaces = false`, `cxx:253`) is the ported
//!   path, so the flag is not exposed.
//! * `KHR_draco_mesh_compression` (`cxx:1094-1161`, `cxx:2092-2124`) — needs the
//!   Draco library; OCCT itself fails the export when it is unavailable
//!   (`cxx:816-822`). No branch is added.
//! * `RWMesh_EdgeIterator` / `RWMesh_VertexIterator` primitives (`Lines` /
//!   `Points` modes, `cxx:521-573`, `cxx:2079-2087`): the node lists this crate
//!   exports from come from shape documents whose leaves are triangulated
//!   faces, so only `TopAbs_FACE` reaches the writer. [`GltfPrimitiveMode`]
//!   keeps the enumerants.
//! * Per-node `LocalTrsf` (`cxx:2349-2422`) — the document readers here bake
//!   instance locations into the geometry (`rwmesh::gltf_document`, the OBJ
//!   writer's note at `brep_exchange::obj_caf_nodes`), so [`GltfCafNode`]
//!   carries an identity transform and `writeNodes` emits no TRS/matrix.
//!   `GltfWriterTrsfFormat` is kept for the day the readers keep locations.
//! * `OSD_Protection` copying in `RWMesh_MaterialMap::CreateTextureFolder`
//!   (`cxx:277-291`) — same gap the OBJ writer documents.

use std::collections::HashMap;

use occt_core::gp::GpTrsf;

use crate::image_texture::ImageTexture;
use crate::rwmesh::coordinate_system::RwMeshCoordinateSystemConverter;
use crate::rwmesh::mesh_name::MeshNameFormat;
use crate::rwmesh::ObjFaceNode;
use crate::tdata_std_named_data::TDataStdNamedData;
use crate::xcaf_prs::XcafPrsStyle;

pub(crate) mod bin_data;
mod json_writer;
mod material_map;
mod sections;

use bin_data::{GltfBufferView, GltfFaceData, GltfFaceGroup};
use json_writer::GltfJsonWriter;
use material_map::GltfMaterialMap;

/// `RWGltf_WriterTrsfFormat` (`RWGltf_WriterTrsfFormat.hxx:21-40`): how a
/// node's transformation is written.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum GltfWriterTrsfFormat {
    /// `RWGltf_WriterTrsfFormat_Compact` (the constructor default,
    /// `RWGltf_CafWriter.cxx:246`): `matrix` only when the transform carries
    /// rotation, scale and translation at once, `rotation`/`scale`/
    /// `translation` otherwise.
    #[default]
    Compact,
    /// `RWGltf_WriterTrsfFormat_Mat4`: always `matrix`.
    Mat4,
    /// `RWGltf_WriterTrsfFormat_TRS`: always `rotation`/`scale`/`translation`.
    Trs,
}

/// `RWGltf_CafWriter`'s setters (`RWGltf_CafWriter.hxx`) with the constructor
/// defaults (`RWGltf_CafWriter.cxx:244-258`).
#[derive(Debug, Clone)]
pub struct GltfCafWriterOptions {
    /// `SetBinaryFormat` / the constructor's `theIsBinary`: write a `.glb`
    /// container instead of `.gltf` + `.bin`.
    pub is_binary: bool,
    /// `SetTransformationFormat()`; `Compact` by default.
    pub trsf_format: GltfWriterTrsfFormat,
    /// `SetNodeNameFormat()`; `InstanceOrProduct` by default.
    pub node_name_format: MeshNameFormat,
    /// `SetMeshNameFormat()`; `Product` by default.
    pub mesh_name_format: MeshNameFormat,
    /// `SetForcedUVExport()`; FALSE by default.
    pub is_forced_uv_export: bool,
    /// `SetToEmbedTexturesInGlb()`; TRUE by default.
    pub to_embed_textures_in_glb: bool,
    /// `myCSTrsf`: the coordinate/length-unit converter, pre-set to the glTF
    /// output system and metres (`cxx:255-258`).
    pub cs_trsf: RwMeshCoordinateSystemConverter,
    /// `SetDefaultStyle()` (`cxx:652`): the material definition used when a
    /// style only carries colours.
    pub default_style: XcafPrsStyle,
}

impl Default for GltfCafWriterOptions {
    /// [`Self::new`] with `theIsBinary = false`.
    fn default() -> Self {
        Self::new(false)
    }
}

impl GltfCafWriterOptions {
    /// `RWGltf_CafWriter::RWGltf_CafWriter(theFile, theIsBinary)`
    /// (`RWGltf_CafWriter.cxx:244-266`), minus the file names (they depend on
    /// the output path [`perform`] receives).
    pub fn new(is_binary: bool) -> Self {
        let mut cs_trsf = RwMeshCoordinateSystemConverter::new();
        cs_trsf.set_output_length_unit(1.0); // meters
        cs_trsf.set_output_coordinate_system(crate::rwmesh::RwMeshCoordinateSystem::GLTF);
        Self {
            is_binary,
            trsf_format: GltfWriterTrsfFormat::Compact,
            node_name_format: MeshNameFormat::InstanceOrProduct,
            mesh_name_format: MeshNameFormat::Product,
            is_forced_uv_export: false,
            to_embed_textures_in_glb: true,
            cs_trsf,
            default_style: XcafPrsStyle::default(),
        }
    }
}

/// One `XCAFPrs_DocumentNode` of the exported document, plus the faces of the
/// node (`RWMesh_ShapeIterator`'s `RWMesh_FaceIterator` walk).
///
/// The names and the entry are the `TDF_Label` / `TDataStd_Name` values
/// `RWMesh::FormatName` reads; see [`format_name`](crate::rwmesh::format_name).
pub struct GltfCafNode {
    /// `XCAFPrs_DocumentNode::Id` — `DefineChildId` path, the node's identity.
    pub id: String,
    /// `XCAFPrs_DocumentNode::Label` — the instance label (index in the
    /// document).
    pub label: usize,
    /// `XCAFPrs_DocumentNode::RefLabel` — the Product label.
    pub ref_label: usize,
    /// `RWMesh::ReadNameAttribute(theLabel)`.
    pub label_name: String,
    /// `RWMesh::ReadNameAttribute(theRefLabel)`.
    pub product_name: String,
    /// `TDF_Tool::Entry(theLabel)`.
    pub entry: String,
    /// `XCAFPrs_DocumentNode::IsAssembly`.
    pub is_assembly: bool,
    /// Index of the parent node in the writer's node list; `None` at depth 0.
    pub parent: Option<usize>,
    /// `XCAFPrs_DocumentNode::LocalTrsf`.
    pub local_trsf: GpTrsf,
    /// `XCAFDoc_ShapeTool::GetNamedProperties(aDocNode.Label)`.
    pub named_data: Option<TDataStdNamedData>,
    /// `XCAFDoc_ShapeTool::GetNamedProperties(aDocNode.RefLabel)`.
    pub ref_named_data: Option<TDataStdNamedData>,
    /// The node's faces, in `RWMesh_FaceIterator` order.
    pub faces: Vec<ObjFaceNode>,
}

/// `RWGltf_CafWriter::toSkipShape` (`RWGltf_CafWriter.cxx:303-306`):
/// `RWMesh_ShapeIterator::IsEmpty`, i.e. no nodes or no triangles.
fn to_skip_shape(face: &ObjFaceNode) -> bool {
    face.vertices.is_empty() || face.triangles.is_empty()
}

/// The `RWGltf_StyledShape` identity (`RWGltf_CafWriter.hxx:491-495`):
/// `Shape.IsSame()` plus `Style.IsEqual()`.
fn same_styled_shape(group: &GltfFaceGroup, face: &ObjFaceNode) -> bool {
    group.shape_key == face.shape_key && group.style.is_equal(&face.style)
}

/// The output of [`perform`]: the document bytes and the binary body.
pub struct GltfCafWriterOutput {
    /// The `.gltf` JSON document, or the assembled `.glb` container.
    pub document: Vec<u8>,
    /// The `.bin` sidecar bytes; empty for a GLB (the body is inside
    /// [`Self::document`]).
    pub bin: Vec<u8>,
    /// The `.bin` file name the JSON `buffers[0].uri` points at
    /// (`<stem>.bin`; `<stem>.bin.tmp` internally for a GLB, `cxx:262-264`).
    pub bin_file_name_short: String,
}

/// `RWGltf_CafWriter` (`RWGltf_CafWriter.cxx`), writing into memory instead of
/// a `std::ostream`.
struct GltfCafWriter<'a> {
    options: &'a GltfCafWriterOptions,
    nodes: &'a [GltfCafNode],
    /// `myBinFileNameShort` (`cxx:262-264`).
    bin_file_name_short: String,
    /// The binary body (`myBinFileNameFull`'s contents).
    bin: Vec<u8>,
    /// `myBuffViewPos` / `myBuffViewNorm` / `myBuffViewTextCoord` /
    /// `myBuffViewInd`.
    buff_view_pos: GltfBufferView,
    buff_view_norm: GltfBufferView,
    buff_view_text_coord: GltfBufferView,
    buff_view_ind: GltfBufferView,
    /// `myBinDataMap` (`cxx:511-515`).
    bin_data_map: Vec<GltfFaceGroup>,
    /// `myBinDataLen64`.
    bin_data_len: i64,
    /// `aNbAccessors` handed to the four passes (`cxx:904`).
    accessor_nb: i32,
    /// `myMaterialMap` (`cxx:651-653`).
    material_map: GltfMaterialMap,
}

impl<'a> GltfCafWriter<'a> {
    /// `RWGltf_CafWriter(theFile, theIsBinary)` (`cxx:244-266`).
    fn new(file_name: &str, nodes: &'a [GltfCafNode], options: &'a GltfCafWriterOptions) -> Self {
        let (_folder, short_file_name) = ImageTexture::folder_and_file(file_name);
        let (short_file_name_base, _ext) = ImageTexture::file_name_and_extension(&short_file_name);
        let bin_file_name_short =
            format!("{short_file_name_base}.bin{}", if options.is_binary { ".tmp" } else { "" });
        Self {
            options,
            nodes,
            bin_file_name_short,
            bin: Vec::new(),
            buff_view_pos: GltfBufferView::new(0, 0),
            buff_view_norm: GltfBufferView::new(0, 0),
            buff_view_text_coord: GltfBufferView::new(0, 0),
            buff_view_ind: GltfBufferView::new(0, 0),
            bin_data_map: Vec::new(),
            bin_data_len: 0,
            accessor_nb: 0,
            material_map: GltfMaterialMap::new(file_name, 0),
        }
    }

    /// `RWGltf_CafWriter::writeBinData` (`cxx:811-1207`), without Draco.
    fn write_bin_data(&mut self) -> Result<(), String> {
        // `cxx:824-846`: the four buffer views' strides and targets.
        self.buff_view_pos = GltfBufferView::new(12, bin_data::TARGET_ARRAY_BUFFER);
        self.buff_view_norm = GltfBufferView::new(12, bin_data::TARGET_ARRAY_BUFFER);
        self.buff_view_text_coord = GltfBufferView::new(8, bin_data::TARGET_ARRAY_BUFFER);
        self.buff_view_ind = GltfBufferView::new(0, bin_data::TARGET_ELEMENT_ARRAY_BUFFER);
        self.bin_data_map.clear();
        self.bin_data_len = 0;
        self.bin.clear();

        // `cxx:874-901`: dispatchShapes over the leaf nodes.
        for (node_index, node) in self.nodes.iter().enumerate() {
            // `XCAFPrs_DocumentExplorerFlags_OnlyLeafNodes` never yields an
            // assembly (`XCAFPrs_DocumentExplorer.cxx:328-437`).
            if node.is_assembly {
                continue;
            }
            for (face_index, face) in node.faces.iter().enumerate() {
                if to_skip_shape(face) {
                    continue;
                }
                if self
                    .bin_data_map
                    .iter()
                    .any(|group| same_styled_shape(group, face))
                {
                    continue;
                }
                self.bin_data_map.push(GltfFaceGroup {
                    shape_key: face.shape_key,
                    style: face.style.clone(),
                    faces: vec![GltfFaceData::new(face.shape_key, node_index, face_index)],
                });
            }
        }

        // `cxx:905-1075`: the four array passes.
        for array_type in 0..4u8 {
            let mut buff_view = match array_type {
                0 => self.buff_view_pos,
                1 => self.buff_view_norm,
                2 => self.buff_view_text_coord,
                _ => self.buff_view_ind,
            };
            buff_view.byte_offset = self.bin.len() as i64;

            // `aWrittenPrimData` (`cxx:907-908`): the accessors of a shape that
            // was already written this pass, so a shape shared by several nodes
            // is serialized once. The C++ stores the `RWGltf_GltfFace` handle
            // and copies the accessor from its *current* state, so the entry
            // is a map index here rather than a snapshot.
            let mut written_prim_data: Vec<(usize, usize)> = Vec::new();
            for group_index in 0..self.bin_data_map.len() {
                let shape_key = self.bin_data_map[group_index].shape_key;
                let node_index = self.bin_data_map[group_index].faces[0].node;
                let face_index = self.bin_data_map[group_index].faces[0].face;
                let face = &self.nodes[node_index].faces[face_index];

                if let Some(position) = written_prim_data
                    .iter()
                    .position(|(key, _)| *key == shape_key)
                {
                    let old_group = written_prim_data[position].1;
                    let old = match array_type {
                        0 => self.bin_data_map[old_group].faces[0].node_pos.clone(),
                        1 => self.bin_data_map[old_group].faces[0].node_norm.clone(),
                        2 => self.bin_data_map[old_group].faces[0].node_uv.clone(),
                        _ => self.bin_data_map[old_group].faces[0].indices.clone(),
                    };
                    let data = &mut self.bin_data_map[group_index].faces[0];
                    match array_type {
                        0 => data.node_pos = old,
                        1 => data.node_norm = old,
                        2 => data.node_uv = old,
                        _ => data.indices = old,
                    }
                    continue;
                }
                written_prim_data.push((shape_key, group_index));

                let data = &mut self.bin_data_map[group_index].faces[0];
                let accessor_nb = &mut self.accessor_nb;
                match array_type {
                    0 => {
                        // `cxx:768-771`: reset before the indices pass.
                        data.nb_indexed_nodes = 0;
                        bin_data::save_nodes(
                            &mut self.bin,
                            data,
                            face,
                            accessor_nb,
                            &self.options.cs_trsf,
                            buff_view.byte_offset,
                        );
                    }
                    1 => bin_data::save_normals(
                        &mut self.bin,
                        data,
                        face,
                        accessor_nb,
                        &self.options.cs_trsf,
                        buff_view.byte_offset,
                    ),
                    2 => bin_data::save_text_coords(
                        &mut self.bin,
                        data,
                        face,
                        accessor_nb,
                        self.options.is_forced_uv_export,
                        buff_view.byte_offset,
                    ),
                    _ => bin_data::save_indices(
                        &mut self.bin,
                        data,
                        face,
                        accessor_nb,
                        buff_view.byte_offset,
                    ),
                }

                // `cxx:1049-1061`: align by 4 bytes (a uint16 index block can
                // leave a 2-byte tail).
                while self.bin.len() % 4 != 0 {
                    self.bin.push(b' ');
                }
            }

            buff_view.byte_length = self.bin.len() as i64 - buff_view.byte_offset;
            match array_type {
                0 => self.buff_view_pos = buff_view,
                1 => self.buff_view_norm = buff_view,
                2 => self.buff_view_text_coord = buff_view,
                _ => self.buff_view_ind = buff_view,
            }
        }

        // `cxx:1073-1092`: buffer view ids, in position/normal/uv/indices
        // order, skipping the empty ones.
        let mut buff_view_id = 0;
        for buff_view in [
            &mut self.buff_view_pos,
            &mut self.buff_view_norm,
            &mut self.buff_view_text_coord,
            &mut self.buff_view_ind,
        ] {
            if buff_view.byte_length > 0 {
                buff_view.id = buff_view_id;
                buff_view_id += 1;
            }
        }

        // `cxx:1163-1194`: embedded textures join the binary body in leaf-node
        // order.
        if self.options.is_binary && self.options.to_embed_textures_in_glb {
            for node in self.nodes {
                if node.is_assembly {
                    continue;
                }
                for face in &node.faces {
                    if to_skip_shape(face) {
                        continue;
                    }
                    self.material_map.add_glb_images(&mut self.bin, &face.style);
                }
            }
        }

        self.bin_data_len = self.bin.len() as i64;
        Ok(())
    }

    /// `RWGltf_CafWriter::writeJson` (`cxx:1210-1426`).
    fn write_json(&mut self, file_info: &[(String, String)]) -> Result<Vec<u8>, String> {
        let bin_data_buffer_id = 0;
        let default_scene_id = 0;

        // `cxx:1250-1302`: the scene node map, i.e. the leaf nodes that carry
        // geometry, in document order.
        let scene_nodes = self.scene_nodes();

        let mut writer = GltfJsonWriter::new();
        writer.start_object();
        self.write_accessors(&mut writer);
        self.write_asset(&mut writer, file_info);
        self.write_buffer_views(&mut writer, bin_data_buffer_id);
        self.write_buffers(&mut writer);
        self.write_images(&mut writer, &scene_nodes);
        self.write_materials(&mut writer, &scene_nodes);
        self.write_meshes(&mut writer, &scene_nodes);
        let scene_root_node_inds = self.write_nodes(&mut writer, &scene_nodes);
        self.write_samplers(&mut writer);
        self.write_scene(&mut writer, default_scene_id);
        self.write_scenes(&mut writer, &scene_root_node_inds);
        self.write_textures(&mut writer, &scene_nodes);
        writer.end_object();

        let mut json = writer.into_bytes();
        if !self.options.is_binary {
            return Ok(json);
        }

        // `cxx:1352-1360`: pad the JSON chunk to 4 bytes.
        while json.len() % 4 != 0 {
            json.push(b' ');
        }
        let bin_length = self.bin_data_len as u32;
        let content_length = json.len() as u32;
        let full_length = content_length as u64 + 20 + self.bin_data_len as u64 + 8;
        if full_length >= u32::MAX as u64 {
            return Err("glTF file content is too big for binary format".to_string());
        }

        let mut out = Vec::with_capacity(full_length as usize);
        // `cxx:1236-1248`: magic, version, total length, chunk length, chunk
        // type `JSON`.
        out.extend_from_slice(b"glTF");
        out.extend_from_slice(&2u32.to_le_bytes());
        out.extend_from_slice(&(full_length as u32).to_le_bytes());
        out.extend_from_slice(&content_length.to_le_bytes());
        out.extend_from_slice(&0x4E4F_534Au32.to_le_bytes());
        out.extend_from_slice(&json);
        // `cxx:1353-1358`: the BIN chunk header.
        out.extend_from_slice(&bin_length.to_le_bytes());
        out.extend_from_slice(&0x004E_4942u32.to_le_bytes());
        out.extend_from_slice(&self.bin);
        Ok(out)
    }

    /// `cxx:1254-1302`: `hasShapeData` per explored leaf node.
    fn scene_nodes(&self) -> Vec<usize> {
        self.nodes
            .iter()
            .enumerate()
            .filter(|(_, node)| !node.is_assembly && !node.faces.is_empty())
            .map(|(index, _)| index)
            .collect()
    }

    /// `aSceneNodeMapWithChildren`'s `children` lists (`cxx:2326-2343`):
    /// `TDF_ChildIterator(RefLabel)` is the document-order child list.
    fn node_children(&self) -> HashMap<usize, Vec<usize>> {
        let mut children: HashMap<usize, Vec<usize>> = HashMap::new();
        for (index, node) in self.nodes.iter().enumerate() {
            if let Some(parent) = node.parent {
                children.entry(parent).or_default().push(index);
            }
        }
        children
    }

    /// `cxx:2429-2433`: the node's mesh index, i.e. its position in the scene
    /// node map.
    fn mesh_index(scene_nodes: &[usize], node_index: usize) -> Option<i32> {
        scene_nodes
            .iter()
            .position(|index| *index == node_index)
            .map(|index| index as i32)
    }
}

/// `RWGltf_CafWriter::writeExtrasAttributes` (`cxx:2472-2584`): the
/// `TDataStd_NamedData` members in the C++'s order — integers, reals, strings,
/// bytes, arrays of integers, arrays of reals.
///
/// The bytes arm is UNPORTED: [`TDataStdNamedData`] carries no byte members
/// (the glTF extras reader has no `uint8` slot either,
/// `rwmesh::gltf_extras`).
pub(crate) fn write_extras_attributes(writer: &mut GltfJsonWriter, named_data: &TDataStdNamedData) {
    for (key, value) in named_data.integers() {
        writer.key(key);
        writer.int(value);
    }
    for (key, value) in named_data.reals() {
        writer.key(key);
        writer.double(value);
    }
    for (key, value) in named_data.strings() {
        writer.key(key);
        writer.string(value);
    }
    for (key, values) in named_data.arrays_of_integers() {
        writer.key(key);
        writer.start_array();
        for value in values {
            writer.int(*value);
        }
        writer.end_array();
    }
    for (key, values) in named_data.arrays_of_reals() {
        writer.key(key);
        writer.start_array();
        for value in values {
            writer.double(*value);
        }
        writer.end_array();
    }
}

/// `RWGltf_CafWriter::myBinFileNameFull` (`RWGltf_CafWriter.cxx:262-264`): the
/// `.bin` sidecar next to the `.gltf` path, i.e. `aFolder + myBinFileNameShort`.
pub fn bin_file_path(gltf_path: &str, bin_file_name_short: &str) -> String {
    match gltf_path.rfind(['/', '\\']) {
        Some(index) => format!("{}{}", &gltf_path[..=index], bin_file_name_short),
        None => bin_file_name_short.to_string(),
    }
}

/// `RWGltf_CafWriter::Perform` (`RWGltf_CafWriter.cxx:664-666`) over an
/// in-memory node list.
///
/// `path` supplies the base name of the `.bin` sidecar; nothing is written to
/// disk (the caller stores [`GltfCafWriterOutput::document`] and
/// [`GltfCafWriterOutput::bin`]), the same in-memory contract as the OBJ and
/// PLY ports. Texture files are copied by the material map, exactly like
/// `RWMesh_MaterialMap::CopyTexture` does for the OBJ writer.
pub fn perform(
    path: &str,
    nodes: &[GltfCafNode],
    file_info: &[(String, String)],
    options: &GltfCafWriterOptions,
) -> Result<GltfCafWriterOutput, String> {
    let mut writer = GltfCafWriter::new(path, nodes, options);
    writer.material_map.set_default_style(options.default_style.clone());
    writer.write_bin_data()?;
    let document = writer.write_json(file_info)?;
    let bin = if options.is_binary { Vec::new() } else { std::mem::take(&mut writer.bin) };
    Ok(GltfCafWriterOutput {
        document,
        bin,
        bin_file_name_short: writer.bin_file_name_short.clone(),
    })
}
