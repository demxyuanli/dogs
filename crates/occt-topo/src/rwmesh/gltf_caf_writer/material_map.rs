//! `RWGltf_GltfMaterialMap` — the `materials` / `images` / `textures` sections
//! of the glTF document.
//!
//! Source: `RWGltf_GltfMaterialMap.hxx` / `.cxx`, over `RWMesh_MaterialMap`
//! (`RWMesh_MaterialMap.cxx:29-232`), which
//! [`MeshMaterialMap`](crate::rwmesh::obj_caf_writer::MeshMaterialMap) shares
//! with the OBJ writer.
//!
//! Differences from the OBJ map that the glTF side adds:
//!
//! * `myMatNameAsKey` is FALSE (`RWGltf_GltfMaterialMap.cxx:58`), so a material
//!   key is the decimal material index and `name` is `mat_<index>`.
//! * `DefineMaterial` (`cxx:356-557`) turns the style into a
//!   `pbrMetallicRoughness` block: base colour, metallic/roughness factors, the
//!   double-sided flag, alpha mode/cut-off, emissive/normal/occlusion textures.
//! * Textures are copied next to the model (`AddImages`, `cxx:65-133`) or, for
//!   a GLB with `myToEmbedTexturesInGlb`, appended to the binary data as image
//!   buffer views (`AddGlbImages`, `cxx:86-174`).
//!
//! `myImageMap` is an `NCollection_IndexedDataMap`, i.e. insertion ordered, and
//! its index is the glTF `images` index (`cxx:320-322`).

use crate::image_texture::ImageTexture;
use crate::rwmesh::obj_caf_writer::MeshMaterialMap;
use crate::rwmesh::gltf_caf_writer::bin_data::GltfBufferView;
use crate::rwmesh::gltf_caf_writer::json_writer::GltfJsonWriter;
use crate::xcaf_doc_vis_material::{
    AlphaMode, LinearRgb, TypeOfBackfacingModel, XcafDocVisMaterial, XcafDocVisMaterialPbr,
};
use crate::xcaf_prs::XcafPrsStyle;

/// `equal_to<handle<Image_Texture>>` (`Image_Texture.hxx:105-118`) plus
/// `NCollection_IndexedDataMap`'s hash: two textures are the same entry when
/// their `TextureId()` matches.
fn same_texture(a: &ImageTexture, b: &ImageTexture) -> bool {
    a.texture_id() == b.texture_id()
}

/// `RWGltf_GltfMaterialMap::baseColorTexture` (`cxx:29-47`): the PBR base
/// colour texture, else the common material's diffuse texture, else none.
pub(crate) fn base_color_texture(material: Option<&XcafDocVisMaterial>) -> Option<&ImageTexture> {
    let material = material?;
    if material.has_pbr_material() {
        if let Some(texture) = material.pbr_material().base_color_texture.as_ref() {
            return Some(texture);
        }
    }
    if material.has_common_material() {
        if let Some(texture) = material.common_material().diffuse_texture.as_ref() {
            return Some(texture);
        }
    }
    None
}

/// `RWGltf_GltfMaterialMap` (`RWGltf_GltfMaterialMap.hxx:29-190`).
pub(crate) struct GltfMaterialMap {
    /// The `RWMesh_MaterialMap` base: key/name generation and texture copy.
    base: MeshMaterialMap,
    /// `myDefSamplerId` (`hxx:181`).
    def_sampler_id: i32,
    /// `myImageMap` (`hxx:185`): texture -> buffer view, insertion ordered.
    image_map: Vec<(ImageTexture, GltfBufferView)>,
    /// `myImageFailMap` (`hxx:186`): textures already known to fail.
    image_fail_map: Vec<ImageTexture>,
    /// `myTextureMap` (`hxx:187`): textures already declared in `textures`.
    texture_map: Vec<ImageTexture>,
    /// `myDefaultStyle` (`RWMesh_MaterialMap.hxx:86`).
    default_style: XcafPrsStyle,
}

impl GltfMaterialMap {
    /// `RWGltf_GltfMaterialMap(const AsciiString&, const int)`
    /// (`RWGltf_GltfMaterialMap.cxx:50-59`).
    pub(crate) fn new(file_name: &str, def_sampler_id: i32) -> Self {
        let mut base = MeshMaterialMap::new(file_name);
        base.set_mat_name_as_key(false);
        Self {
            base,
            def_sampler_id,
            image_map: Vec::new(),
            image_fail_map: Vec::new(),
            texture_map: Vec::new(),
            default_style: XcafPrsStyle::default(),
        }
    }

    /// `SetDefaultStyle` (`RWMesh_MaterialMap.hxx:71-74`).
    pub(crate) fn set_default_style(&mut self, style: XcafPrsStyle) {
        self.default_style = style;
    }

    /// `NbImages()` (`RWMesh_MaterialMap.hxx:51-54`): `myImageMap.Extent()`.
    pub(crate) fn nb_images(&self) -> usize {
        self.image_map.len()
    }

    /// `FindMaterial(style)` (`RWMesh_MaterialMap.hxx:41-48`): the material key
    /// of an already defined style (`""` when it has none).
    pub(crate) fn find_material(&self, style: &XcafPrsStyle) -> String {
        self.base.find_material(style).unwrap_or_default()
    }

    /// `RWGltf_GltfMaterialMap::AddImages` (`cxx:65-79`): the five texture
    /// slots of a material, copied next to the model.
    pub(crate) fn add_images(
        &mut self,
        writer: &mut GltfJsonWriter,
        style: &XcafPrsStyle,
        is_started: &mut bool,
    ) {
        let Some(material) = style.material().filter(|m| !m.is_empty()) else {
            return;
        };
        let pbr = material.pbr_material();
        let base_color = base_color_texture(Some(material)).cloned();
        self.add_image(writer, base_color.as_ref(), is_started);
        self.add_image(writer, pbr.metallic_roughness_texture.as_ref(), is_started);
        self.add_image(writer, pbr.normal_texture.as_ref(), is_started);
        self.add_image(writer, pbr.emissive_texture.as_ref(), is_started);
        self.add_image(writer, pbr.occlusion_texture.as_ref(), is_started);
    }

    /// `RWGltf_GltfMaterialMap::AddGlbImages` (`cxx:86-95`): the same five
    /// slots, written into the GLB binary body.
    pub(crate) fn add_glb_images(&mut self, bin: &mut Vec<u8>, style: &XcafPrsStyle) {
        let Some(material) = style.material().filter(|m| !m.is_empty()) else {
            return;
        };
        let pbr = material.pbr_material();
        let base_color = base_color_texture(Some(material)).cloned();
        self.add_glb_image(bin, base_color.as_ref());
        self.add_glb_image(bin, pbr.metallic_roughness_texture.as_ref());
        self.add_glb_image(bin, pbr.normal_texture.as_ref());
        self.add_glb_image(bin, pbr.emissive_texture.as_ref());
        self.add_glb_image(bin, pbr.occlusion_texture.as_ref());
    }

    /// `RWGltf_GltfMaterialMap::addImage` (`cxx:96-133`).
    fn add_image(
        &mut self,
        writer: &mut GltfJsonWriter,
        texture: Option<&ImageTexture>,
        is_started: &mut bool,
    ) {
        let Some(texture) = texture else { return };
        if self.image_map.iter().any(|(t, _)| same_texture(t, texture))
            || self.image_fail_map.iter().any(|t| same_texture(t, texture))
        {
            return;
        }

        // `const TCollection_AsciiString aGltfImgKey = myImageMap.Extent();`
        let gltf_img_key = self.image_map.len().to_string();
        // `if (!CopyTexture(aTextureUri, theTexture, aGltfImgKey))`
        let (texture_uri, ok) = self.base.copy_texture(texture, &gltf_img_key);
        if !ok {
            self.image_fail_map.push(texture.clone());
            return;
        }
        self.image_map.push((texture.clone(), GltfBufferView::new(0, 0)));

        if !*is_started {
            writer.key("images");
            writer.start_array();
            *is_started = true;
        }
        writer.start_object();
        writer.key("uri");
        writer.string(&texture_uri);
        writer.end_object();
    }

    /// `RWGltf_GltfMaterialMap::addGlbImage` (`cxx:140-173`).
    fn add_glb_image(&mut self, bin: &mut Vec<u8>, texture: Option<&ImageTexture>) {
        let Some(texture) = texture else { return };
        if self.image_map.iter().any(|(t, _)| same_texture(t, texture))
            || self.image_fail_map.iter().any(|t| same_texture(t, texture))
        {
            return;
        }

        let mut buff_image = GltfBufferView::new(0, 0);
        buff_image.byte_offset = bin.len() as i64;
        if !texture.write_image_stream(bin) {
            self.image_fail_map.push(texture.clone());
            return;
        }
        // alignment by 4 bytes
        while bin.len() % 4 != 0 {
            bin.push(b' ');
        }
        buff_image.byte_length = bin.len() as i64 - buff_image.byte_offset;
        if buff_image.byte_length <= 0 {
            self.image_fail_map.push(texture.clone());
            return;
        }
        self.image_map.push((texture.clone(), buff_image));
    }

    /// `RWGltf_GltfMaterialMap::FlushGlbBufferViews` (`cxx:177-207`): the image
    /// buffer views go at the end of the `bufferViews` array.
    pub(crate) fn flush_glb_buffer_views(
        &mut self,
        writer: &mut GltfJsonWriter,
        bin_data_buffer_id: i32,
        buff_view_id: &mut i32,
    ) {
        for (_, buff_view) in self.image_map.iter_mut() {
            if buff_view.byte_length <= 0 {
                continue;
            }
            buff_view.id = *buff_view_id;
            *buff_view_id += 1;
            writer.start_object();
            writer.key("buffer");
            writer.int(bin_data_buffer_id);
            writer.key("byteLength");
            writer.int64(buff_view.byte_length);
            writer.key("byteOffset");
            writer.int64(buff_view.byte_offset);
            writer.end_object();
        }
    }

    /// `RWGltf_GltfMaterialMap::FlushGlbImages` (`cxx:212-257`).
    pub(crate) fn flush_glb_images(&mut self, writer: &mut GltfJsonWriter) {
        let mut is_started = false;
        for (texture, buff_view) in &self.image_map {
            if buff_view.byte_length <= 0 {
                continue;
            }
            if !is_started {
                writer.key("images");
                writer.start_array();
                is_started = true;
            }
            writer.start_object();
            // `mimeType` first, then `bufferView` (`cxx:245-248`).
            writer.key("mimeType");
            writer.string(&texture.mime_type());
            writer.key("bufferView");
            writer.int(buff_view.id);
            writer.end_object();
        }
        if is_started {
            writer.end_array();
        }
    }

    /// `RWGltf_GltfMaterialMap::AddMaterial(writer, style, isStarted)`
    /// (`cxx:262-287`).
    pub(crate) fn add_material(
        &mut self,
        writer: &mut GltfJsonWriter,
        style: &XcafPrsStyle,
        is_started: &mut bool,
    ) {
        let has_material = style.material().map(|m| !m.is_empty()).unwrap_or(false);
        if !has_material && !style.is_set_color_surf() && !style.is_set_color_curv() {
            return;
        }
        if !*is_started {
            writer.key("materials");
            writer.start_array();
            *is_started = true;
        }
        let (_, name, is_new) = self.base.add_material(style);
        if is_new {
            self.define_material(writer, style, &name);
        }
    }

    /// `RWGltf_GltfMaterialMap::AddTextures` (`cxx:292-306`).
    pub(crate) fn add_textures(
        &mut self,
        writer: &mut GltfJsonWriter,
        style: &XcafPrsStyle,
        is_started: &mut bool,
    ) {
        let Some(material) = style.material().filter(|m| !m.is_empty()) else {
            return;
        };
        let pbr = material.pbr_material();
        let base_color = base_color_texture(Some(material)).cloned();
        self.add_texture(writer, base_color.as_ref(), is_started);
        self.add_texture(writer, pbr.metallic_roughness_texture.as_ref(), is_started);
        self.add_texture(writer, pbr.normal_texture.as_ref(), is_started);
        self.add_texture(writer, pbr.emissive_texture.as_ref(), is_started);
        self.add_texture(writer, pbr.occlusion_texture.as_ref(), is_started);
    }

    /// `RWGltf_GltfMaterialMap::addTexture` (`cxx:310-342`).
    fn add_texture(
        &mut self,
        writer: &mut GltfJsonWriter,
        texture: Option<&ImageTexture>,
        is_started: &mut bool,
    ) {
        let Some(texture) = texture else { return };
        if self.texture_map.iter().any(|t| same_texture(t, texture))
            || !self.image_map.iter().any(|(t, _)| same_texture(t, texture))
        {
            return;
        }
        // `myImageMap.FindIndex(theTexture) - 1` (`cxx:320-322`).
        let image_key = self
            .image_map
            .iter()
            .position(|(t, _)| same_texture(t, texture))
            .unwrap_or(0) as i32;
        self.texture_map.push(texture.clone());

        if !*is_started {
            writer.key("textures");
            writer.start_array();
            *is_started = true;
        }
        writer.start_object();
        writer.key("sampler");
        writer.int(self.def_sampler_id);
        writer.key("source");
        writer.int(image_key);
        writer.end_object();
    }

    /// `myImageMap.FindIndex(texture) - 1`, `None` when the texture is absent
    /// (`-1` in the C++).
    fn image_index(&self, texture: Option<&ImageTexture>) -> Option<i32> {
        let texture = texture?;
        self.image_map
            .iter()
            .position(|(t, _)| same_texture(t, texture))
            .map(|index| index as i32)
    }

    /// `RWGltf_GltfMaterialMap::DefineMaterial` (`cxx:356-557`).
    fn define_material(
        &mut self,
        writer: &mut GltfJsonWriter,
        style: &XcafPrsStyle,
        name: &str,
    ) {
        let has_material = style.material().map(|m| !m.is_empty()).unwrap_or(false);
        let mut pbr = if has_material {
            style
                .material()
                .map(|m| m.convert_to_pbr_material())
                .unwrap_or_else(XcafDocVisMaterialPbr::default)
        } else if let Some(default_material) =
            self.default_style.material().filter(|m| m.has_pbr_material())
        {
            default_material.pbr_material().clone()
        } else {
            XcafDocVisMaterialPbr::default()
        };

        if style.is_set_color_surf() {
            // `aPbrMat.BaseColor.SetRGB(theStyle.GetColorSurf())` — `SetRGB`
            // takes the colour's linear components and keeps the alpha.
            pbr.base_color.rgb = LinearRgb::from_srgb(style.get_color_surf());
            let rgba = style.get_color_surf_rgba();
            if rgba.a < 1.0 {
                pbr.metallic = 0.0;
                pbr.base_color.a = rgba.a;
            }
        } else if style.is_set_color_curv() {
            pbr.base_color.rgb = LinearRgb::from_srgb(style.get_color_curv());
        }

        writer.start_object();
        writer.key("name");
        writer.string(name);

        writer.key("pbrMetallicRoughness");
        writer.start_object();
        writer.key("baseColorFactor");
        writer.start_array();
        writer.double(pbr.base_color.rgb.r as f64);
        writer.double(pbr.base_color.rgb.g as f64);
        writer.double(pbr.base_color.rgb.b as f64);
        writer.double(pbr.base_color.a as f64);
        writer.end_array();

        let material = style.material();
        if let Some(base_texture) = base_color_texture(material) {
            if let Some(index) = self.image_index(Some(base_texture)) {
                writer.key("baseColorTexture");
                writer.start_object();
                writer.key("index");
                writer.int(index);
                writer.end_object();
            }
        }
        if has_material || pbr.metallic != 1.0 {
            writer.key("metallicFactor");
            writer.double(pbr.metallic as f64);
        }
        let metallic_roughness_index =
            self.image_index(pbr.metallic_roughness_texture.as_ref());
        if let Some(index) = metallic_roughness_index {
            writer.key("metallicRoughnessTexture");
            writer.start_object();
            writer.key("index");
            writer.int(index);
            writer.end_object();
        }
        if has_material || pbr.roughness != 1.0 {
            writer.key("roughnessFactor");
            writer.double(pbr.roughness as f64);
        }
        writer.end_object();

        // `cxx:452-464`: auto / double-sided / front-culled styles are written
        // as double-sided; front culling cannot be expressed in glTF.
        let double_sided = match material.map(|m| m.face_culling()) {
            None => true,
            Some(TypeOfBackfacingModel::Auto)
            | Some(TypeOfBackfacingModel::DoubleSided)
            | Some(TypeOfBackfacingModel::FrontCulled) => true,
            Some(TypeOfBackfacingModel::BackCulled) => false,
        };
        if double_sided {
            writer.key("doubleSided");
            writer.boolean(true);
        }

        let alpha_mode = match material {
            Some(m) => m.alpha_mode(),
            None => AlphaMode::BlendAuto,
        };
        match alpha_mode {
            AlphaMode::BlendAuto => {
                if pbr.base_color.a < 1.0 {
                    writer.key("alphaMode");
                    writer.string("BLEND");
                }
            }
            AlphaMode::Opaque => {}
            AlphaMode::Mask => {
                writer.key("alphaMode");
                writer.string("MASK");
            }
            AlphaMode::Blend | AlphaMode::MaskBlend => {
                writer.key("alphaMode");
                writer.string("BLEND");
            }
        }
        if let Some(material) = material {
            if material.alpha_cut_off() != 0.5 {
                writer.key("alphaCutoff");
                writer.double(material.alpha_cut_off() as f64);
            }
        }

        if pbr.emissive_factor != [0.0, 0.0, 0.0] {
            writer.key("emissiveFactor");
            writer.start_array();
            writer.double(pbr.emissive_factor[0] as f64);
            writer.double(pbr.emissive_factor[1] as f64);
            writer.double(pbr.emissive_factor[2] as f64);
            writer.end_array();
        }
        if let Some(index) = self.image_index(pbr.emissive_texture.as_ref()) {
            writer.key("emissiveTexture");
            writer.start_object();
            writer.key("index");
            writer.int(index);
            writer.end_object();
        }
        if let Some(index) = self.image_index(pbr.normal_texture.as_ref()) {
            writer.key("normalTexture");
            writer.start_object();
            writer.key("index");
            writer.int(index);
            writer.end_object();
        }
        if let Some(index) = self.image_index(pbr.occlusion_texture.as_ref()) {
            writer.key("occlusionTexture");
            writer.start_object();
            writer.key("index");
            writer.int(index);
            writer.end_object();
        }
        writer.end_object();
    }
}
