//! The member writers `RWGltf_CafWriter::writeJson` calls, in the C++'s order:
//! `writeAccessors` … `writeTextures` (`RWGltf_CafWriter.cxx:1434-2694`).
//!
//! Every member is written only when it has content, because glTF's validator
//! rejects empty arrays; the C++ spells that as an `isStarted` flag handed to
//! the material map (`writeImages` / `writeMaterials` / `writeTextures`) and as
//! an `INVALID_ID` check for the buffer views.

use crate::rwmesh::gltf_caf_writer::bin_data::{GltfBufferView, GltfFaceData, INVALID_ID};
use crate::rwmesh::gltf_caf_writer::json_writer::GltfJsonWriter;
use crate::rwmesh::gltf_caf_writer::{to_skip_shape, GltfCafWriter, GltfWriterTrsfFormat};
use crate::rwmesh::mesh_name::format_name;

/// `Precision::Confusion()` (`cxx:2361-2365`).
const CONFUSION: f64 = 1e-7;

impl GltfCafWriter<'_> {
    /// `writeAccessors` (`cxx:1434-1541`): one accessor per filled
    /// `RWGltf_GltfAccessor`, in id order.
    pub(crate) fn write_accessors(&self, writer: &mut GltfJsonWriter) {
        writer.key("accessors");
        writer.start_array();

        // `aWrittenIds` (`cxx:1447`) lives across the four passes.
        let mut written_ids: Vec<i32> = Vec::new();
        // `aWrittenFaces` (`cxx:1445-1446`) is per pass, but with
        // `toMergeFaces == false` every map entry owns a distinct list, so the
        // `Add` there always succeeds and the check has nothing to do.
        let mut nb_accessors = 0i32;
        for array_type in 0..4u8 {
            for group in &self.bin_data_map {
                for face in &group.faces {
                    let id = accessor_id(face, array_type);
                    // `cxx:1474-1478`.
                    if id == INVALID_ID || written_ids.contains(&id) {
                        continue;
                    }
                    written_ids.push(id);
                    // `cxx:1479-1483`: the ids are handed out sequentially by
                    // `writeBinData`, so a mismatch is an internal error there.
                    if id != nb_accessors {
                        continue;
                    }
                    nb_accessors += 1;
                    match array_type {
                        0 => self.write_positions(writer, face),
                        1 => self.write_normals(writer, face),
                        2 => self.write_text_coords(writer, face),
                        _ => self.write_indices(writer, face),
                    }
                }
            }
        }

        writer.end_array();
    }

    /// `writePositions` (`cxx:1547-1589`).
    fn write_positions(&self, writer: &mut GltfJsonWriter, face: &GltfFaceData) {
        if face.node_pos.id == INVALID_ID {
            return;
        }
        writer.start_object();
        writer.key("bufferView");
        writer.int(self.buff_view_pos.id);
        writer.key("byteOffset");
        writer.int64(face.node_pos.byte_offset);
        writer.key("componentType");
        writer.int(face.node_pos.comp_type.value());
        writer.key("count");
        writer.int64(face.node_pos.count);
        if let Some((min, max)) = face.node_pos.bnd_box {
            // `cxx:1569-1584`: `max` before `min`, from `Bnd_Box::CornerMax`
            // and `CornerMin`.
            writer.key("max");
            writer.start_array();
            for value in max {
                writer.double(value);
            }
            writer.end_array();
            writer.key("min");
            writer.start_array();
            for value in min {
                writer.double(value);
            }
            writer.end_array();
        }
        writer.key("type");
        writer.string(face.node_pos.layout.as_str());
        writer.end_object();
    }

    /// `writeNormals` (`cxx:1595-1642`): no `max`/`min` (`cxx:1617-1633`).
    fn write_normals(&self, writer: &mut GltfJsonWriter, face: &GltfFaceData) {
        if face.node_norm.id == INVALID_ID {
            return;
        }
        writer.start_object();
        writer.key("bufferView");
        writer.int(self.buff_view_norm.id);
        writer.key("byteOffset");
        writer.int64(face.node_norm.byte_offset);
        writer.key("componentType");
        writer.int(face.node_norm.comp_type.value());
        writer.key("count");
        writer.int64(face.node_norm.count);
        writer.key("type");
        writer.string(face.node_norm.layout.as_str());
        writer.end_object();
    }

    /// `writeTextCoords` (`cxx:1647-1694`).
    fn write_text_coords(&self, writer: &mut GltfJsonWriter, face: &GltfFaceData) {
        if face.node_uv.id == INVALID_ID {
            return;
        }
        writer.start_object();
        writer.key("bufferView");
        writer.int(self.buff_view_text_coord.id);
        writer.key("byteOffset");
        writer.int64(face.node_uv.byte_offset);
        writer.key("componentType");
        writer.int(face.node_uv.comp_type.value());
        writer.key("count");
        writer.int64(face.node_uv.count);
        writer.key("type");
        writer.string(face.node_uv.layout.as_str());
        writer.end_object();
    }

    /// `writeIndices` (`cxx:1699-1725`).
    fn write_indices(&self, writer: &mut GltfJsonWriter, face: &GltfFaceData) {
        if face.indices.id == INVALID_ID {
            return;
        }
        writer.start_object();
        writer.key("bufferView");
        writer.int(self.buff_view_ind.id);
        writer.key("byteOffset");
        writer.int64(face.indices.byte_offset);
        writer.key("componentType");
        writer.int(face.indices.comp_type.value());
        writer.key("count");
        writer.int64(face.indices.count);
        writer.key("type");
        writer.string(face.indices.layout.as_str());
        writer.end_object();
    }

    /// `writeAsset` (`cxx:1744-1780`): the generator string, the glTF version
    /// and the `theFileInfo` map as `extras`.
    pub(crate) fn write_asset(&self, writer: &mut GltfJsonWriter, file_info: &[(String, String)]) {
        writer.key("asset");
        writer.start_object();
        writer.key("generator");
        // `"Open CASCADE Technology " OCC_VERSION_STRING " [dev.opencascade.org]"`
        // (`cxx:1752`); the baselines in `data/occ-*.gltf` were written by the
        // 8.0 build, so the version string is spelled the way they carry it.
        writer.string("Open CASCADE Technology 8.0 [dev.opencascade.org]");
        writer.key("version");
        writer.string("2.0");
        if !file_info.is_empty() {
            writer.key("extras");
            writer.start_object();
            for (key, value) in file_info {
                writer.key(key);
                writer.string(value);
            }
            writer.end_object();
        }
        writer.end_object();
    }

    /// `writeBufferViews` (`cxx:1785-1880`) without the Draco views; the
    /// embedded image views follow through
    /// [`GltfMaterialMap::flush_glb_buffer_views`](super::material_map::GltfMaterialMap::flush_glb_buffer_views).
    pub(crate) fn write_buffer_views(&mut self, writer: &mut GltfJsonWriter, bin_data_buffer_id: i32) {
        let mut buff_view_id = 0;
        writer.key("bufferViews");
        writer.start_array();
        // `cxx:1806-1856`: the index view has no `byteStride`.
        let views = [
            (self.buff_view_pos, true),
            (self.buff_view_norm, true),
            (self.buff_view_text_coord, true),
            (self.buff_view_ind, false),
        ];
        for (view, has_stride) in views {
            if view.id == INVALID_ID {
                continue;
            }
            buff_view_id += 1;
            write_buffer_view(writer, bin_data_buffer_id, &view, has_stride);
        }
        self.material_map.flush_glb_buffer_views(writer, bin_data_buffer_id, &mut buff_view_id);
        writer.end_array();
    }

    /// `writeBuffers` (`cxx:1886-1908`): the single binary body, plus its
    /// `.bin` sidecar name for the non-binary flavour.
    pub(crate) fn write_buffers(&self, writer: &mut GltfJsonWriter) {
        writer.key("buffers");
        writer.start_array();
        writer.start_object();
        writer.key("byteLength");
        writer.int64(self.bin_data_len);
        if !self.options.is_binary {
            writer.key("uri");
            writer.string(&self.bin_file_name_short);
        }
        writer.end_object();
        writer.end_array();
    }

    /// `writeImages` (`cxx:1943-1977`): embedded buffer views for a GLB, copied
    /// texture files otherwise.
    pub(crate) fn write_images(&mut self, writer: &mut GltfJsonWriter, scene_nodes: &[usize]) {
        if self.options.is_binary && self.options.to_embed_textures_in_glb {
            self.material_map.flush_glb_images(writer);
            return;
        }
        let mut is_started = false;
        for &node_index in scene_nodes {
            for face in &self.nodes[node_index].faces {
                self.material_map.add_images(writer, &face.style, &mut is_started);
            }
        }
        if is_started {
            writer.end_array();
        }
    }

    /// `writeMaterials` (`cxx:1995-2022`): the face iterator is the only one
    /// with styles, so the edge/vertex arms have nothing to add.
    pub(crate) fn write_materials(&mut self, writer: &mut GltfJsonWriter, scene_nodes: &[usize]) {
        let mut is_started = false;
        for &node_index in scene_nodes {
            for face in &self.nodes[node_index].faces {
                self.material_map.add_material(writer, &face.style, &mut is_started);
            }
        }
        if is_started {
            writer.end_array();
        }
    }

    /// `writeMeshes` (`cxx:2185-2268`) without `toMergeFaces` / Draco.
    pub(crate) fn write_meshes(&self, writer: &mut GltfJsonWriter, scene_nodes: &[usize]) {
        writer.key("meshes");
        writer.start_array();
        for &node_index in scene_nodes {
            let node = &self.nodes[node_index];
            let name = format_name(
                self.options.mesh_name_format,
                &node.label_name,
                &node.product_name,
                &node.entry,
            );

            let mut to_start_prims = true;
            // `aWrittenShapes` (`cxx:2203`) is per node: a styled shape is
            // declared once in the node's `primitives`.
            let mut written_shapes: Vec<usize> = Vec::new();
            for face_index in 0..node.faces.len() {
                let face = &node.faces[face_index];
                if to_skip_shape(face) {
                    continue;
                }
                // `writeShapes` (`cxx:2147-2158`): the map lookup by
                // `RWGltf_StyledShape`.
                let Some(group_index) = self.bin_data_map.iter().position(|group| {
                    group.shape_key == face.shape_key && group.style.is_equal(&face.style)
                }) else {
                    continue;
                };
                if written_shapes.contains(&group_index) {
                    continue;
                }
                written_shapes.push(group_index);
                // `aGltfShapeList->First()` (`cxx:2161`).
                self.write_prim_array(writer, &self.bin_data_map[group_index].faces[0], &name, &mut to_start_prims);
            }
            if !to_start_prims {
                writer.end_array();
                writer.end_object();
            }
        }
        writer.end_array();
    }

    /// `writePrimArray` (`cxx:2027-2132`) without the Draco extension.
    fn write_prim_array(
        &self,
        writer: &mut GltfJsonWriter,
        face: &GltfFaceData,
        name: &str,
        to_start_prims: &mut bool,
    ) {
        if *to_start_prims {
            *to_start_prims = false;
            writer.start_object();
            if !name.is_empty() {
                writer.key("name");
                writer.string(name);
            }
            writer.key("primitives");
            writer.start_array();
        }

        let material_id = self
            .material_map
            .find_material(&self.node_face_style(face));
        writer.start_object();
        writer.key("attributes");
        writer.start_object();
        if face.node_norm.id != INVALID_ID {
            writer.key("NORMAL");
            writer.int(face.node_norm.id);
        }
        writer.key("POSITION");
        writer.int(face.node_pos.id);
        if face.node_uv.id != INVALID_ID {
            writer.key("TEXCOORD_0");
            writer.int(face.node_uv.id);
        }
        writer.end_object();

        writer.key("indices");
        writer.int(face.indices.id);
        if !material_id.is_empty() {
            // `writeMeshPrimitive` (`cxx:2050-2055`): `"material"` then the
            // decimal key `AsciiString::IntegerValue()`.
            writer.key("material");
            writer.int(material_id.parse::<i32>().unwrap_or(0));
        }
        writer.key("mode");
        // `getShapeType` is `TopAbs_FACE` for every document face (`cxx:2080-2089`).
        writer.int(super::bin_data::PRIMITIVE_MODE.value());
        writer.end_object();
    }

    /// The style of the document face a `GltfFaceData` was dispatched from.
    fn node_face_style(&self, face: &GltfFaceData) -> crate::xcaf_prs::XcafPrsStyle {
        self.nodes[face.node].faces[face.face].style.clone()
    }

    /// `writeNodes` (`cxx:2275-2464`): every document node, assemblies first in
    /// depth-first order. Returns the scene's root node indices.
    pub(crate) fn write_nodes(&self, writer: &mut GltfJsonWriter, scene_nodes: &[usize]) -> Vec<i32> {
        let children = self.node_children();

        writer.key("nodes");
        writer.start_array();

        let mut root_indices = Vec::new();
        for (node_index, node) in self.nodes.iter().enumerate() {
            // `cxx:2304-2308`: `CurrentDepth() == 0`.
            if node.parent.is_none() {
                root_indices.push(node_index as i32);
            }
            writer.start_object();

            // `cxx:2322-2345`: an assembly lists its children by document id.
            if node.is_assembly {
                writer.key("children");
                writer.start_array();
                if let Some(child_indices) = children.get(&node_index) {
                    for &child in child_indices {
                        writer.int(child as i32);
                    }
                }
                writer.end_array();
            }

            // `cxx:2349-2423`: `LocalTrsf`, skipped when the node carries none.
            self.write_node_transformation(writer, node);

            // `cxx:2425-2433`: the mesh index is the node's position in the
            // scene node map.
            if !node.is_assembly {
                if let Some(mesh_index) = Self::mesh_index(scene_nodes, node_index) {
                    writer.key("mesh");
                    writer.int(mesh_index);
                }
            }

            let name = format_name(
                self.options.node_name_format,
                &node.label_name,
                &node.product_name,
                &node.entry,
            );
            if !name.is_empty() {
                writer.key("name");
                writer.string(&name);
            }

            // `cxx:2445-2455`: both label and Product attributes, the label's
            // first.
            if node.named_data.is_some() || node.ref_named_data.is_some() {
                writer.key("extras");
                writer.start_object();
                if let Some(named_data) = &node.named_data {
                    super::write_extras_attributes(writer, named_data);
                }
                if let Some(named_data) = &node.ref_named_data {
                    super::write_extras_attributes(writer, named_data);
                }
                writer.end_object();
            }

            writer.end_object();
        }

        writer.end_array();
        root_indices
    }

    /// The transformation block of `writeNodes` (`cxx:2349-2423`).
    fn write_node_transformation(&self, writer: &mut GltfJsonWriter, node: &super::GltfCafNode) {
        if node.local_trsf.form() == occt_core::gp::TrsfForm::Identity {
            return;
        }
        let mut trsf = node.local_trsf.clone();
        self.options.cs_trsf.transform_transformation(&mut trsf);

        let quaternion = trsf.get_rotation();
        let resolution = 1e-15; // `gp::Resolution()`
        let has_rotation = quaternion.x().abs() > resolution
            || quaternion.y().abs() > resolution
            || quaternion.z().abs() > resolution
            || (quaternion.w() - 1.0).abs() > resolution;
        let scale_factor = trsf.scale_factor();
        let has_scale = (scale_factor - 1.0).abs() > CONFUSION;
        let translation = trsf.translation_part();
        let has_translation = translation.square_modulus() > resolution;

        let format = match self.options.trsf_format {
            GltfWriterTrsfFormat::Compact => {
                if has_rotation && has_scale && has_translation {
                    GltfWriterTrsfFormat::Mat4
                } else {
                    GltfWriterTrsfFormat::Trs
                }
            }
            other => other,
        };
        if format == GltfWriterTrsfFormat::Mat4 {
            // `cxx:2375-2389`: column-major 4x4.
            let matrix = trsf.get_mat4();
            let is_identity = (0..4).all(|col| {
                (0..4).all(|row| {
                    let expected = if row == col { 1.0 } else { 0.0 };
                    matrix[col][row] == expected
                })
            });
            if is_identity {
                return;
            }
            writer.key("matrix");
            writer.start_array();
            for col in 0..4 {
                for row in 0..4 {
                    writer.double(matrix[col][row]);
                }
            }
            writer.end_array();
            return;
        }
        if has_rotation {
            writer.key("rotation");
            writer.start_array();
            writer.double(quaternion.x());
            writer.double(quaternion.y());
            writer.double(quaternion.z());
            writer.double(quaternion.w());
            writer.end_array();
        }
        if has_scale {
            writer.key("scale");
            writer.start_array();
            for _ in 0..3 {
                writer.double(scale_factor);
            }
            writer.end_array();
        }
        if has_translation {
            writer.key("translation");
            writer.start_array();
            writer.double(translation.x());
            writer.double(translation.y());
            writer.double(translation.z());
            writer.end_array();
        }
    }

    /// `writeSamplers` (`cxx:2585-2608`): a single empty sampler, only when
    /// textures were added.
    pub(crate) fn write_samplers(&self, writer: &mut GltfJsonWriter) {
        if self.material_map.nb_images() == 0 {
            return;
        }
        writer.key("samplers");
        writer.start_array();
        writer.start_object();
        writer.end_object();
        writer.end_array();
    }

    /// `writeScene` (`cxx:2614-2624`).
    pub(crate) fn write_scene(&self, writer: &mut GltfJsonWriter, def_scene_id: i32) {
        writer.key("scene");
        writer.int(def_scene_id);
    }

    /// `writeScenes` (`cxx:2629-2653`).
    pub(crate) fn write_scenes(&self, writer: &mut GltfJsonWriter, root_indices: &[i32]) {
        writer.key("scenes");
        writer.start_array();
        writer.start_object();
        writer.key("nodes");
        writer.start_array();
        for root in root_indices {
            writer.int(*root);
        }
        writer.end_array();
        writer.end_object();
        writer.end_array();
    }

    /// `writeTextures` (`cxx:2667-2694`).
    pub(crate) fn write_textures(&mut self, writer: &mut GltfJsonWriter, scene_nodes: &[usize]) {
        let mut is_started = false;
        for &node_index in scene_nodes {
            for face in &self.nodes[node_index].faces {
                self.material_map.add_textures(writer, &face.style, &mut is_started);
            }
        }
        if is_started {
            writer.end_array();
        }
    }
}

/// The `RWGltf_GltfAccessor` a type pass writes (`switch` at `cxx:1468-1533`).
fn accessor_id(face: &GltfFaceData, array_type: u8) -> i32 {
    match array_type {
        0 => face.node_pos.id,
        1 => face.node_norm.id,
        2 => face.node_uv.id,
        _ => face.indices.id,
    }
}

/// One `bufferViews` entry (`cxx:1795-1855`).
fn write_buffer_view(
    writer: &mut GltfJsonWriter,
    bin_data_buffer_id: i32,
    view: &GltfBufferView,
    has_stride: bool,
) {
    writer.start_object();
    writer.key("buffer");
    writer.int(bin_data_buffer_id);
    writer.key("byteLength");
    writer.int64(view.byte_length);
    writer.key("byteOffset");
    writer.int64(view.byte_offset);
    if has_stride {
        writer.key("byteStride");
        writer.int64(view.byte_stride as i64);
    }
    writer.key("target");
    writer.int(view.target);
    writer.end_object();
}
