//! glTF material definitions and their binding to `XCAFDoc_VisMaterial`.
//!
//! Source: `RWGltf_MaterialCommon.hxx`, `RWGltf_MaterialMetallicRoughness.hxx`,
//! `RWGltf_GltfAlphaMode.hxx`, and `RWGltf_GltfJsonParser.cxx`:
//! `gltfParseMaterials` (`:793-863`), `gltfBindMaterial` (`:867-947`),
//! `gltfParseStdMaterial` (`:951-1017`), `gltfParsePbrMaterial` (`:1021-1137`)
//! and `gltfParseTexture` (`:1164-1311`).
//!
//! The two glTF material flavours coexist: glTF 2.0 defines a
//! metallic-roughness material, glTF 1.0 an obsolete common (Blinn-Phong) one.
//! OCCT parses both for every material entry and binds them into one
//! `XCAFDoc_VisMaterial`, whose PBR and common halves then feed the styles.
use std::path::Path;

use occt_core::quantity::Color;

use crate::image_texture::ImageTexture;
use crate::xcaf_doc_vis_material::{
    AlphaMode, LinearRgb, LinearRgba, TypeOfBackfacingModel, XcafDocVisMaterial,
    XcafDocVisMaterialCommon, XcafDocVisMaterialPbr, srgb_to_linear,
};

use super::gltf::{base64_decode, JVal};
use super::gltf_types::BufView;
use super::Material;

/// `RWGltf_GltfAlphaMode` (`RWGltf_GltfAlphaMode.hxx:21-27`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub(super) enum GltfAlphaMode {
    /// `RWGltf_GltfAlphaMode_Opaque`.
    #[default]
    Opaque,
    /// `RWGltf_GltfAlphaMode_Mask`.
    Mask,
    /// `RWGltf_GltfAlphaMode_Blend`.
    Blend,
}

/// `RWGltf_GltfParseAlphaMode` (`RWGltf_GltfAlphaMode.hxx:30-46`): unknown
/// strings fall back to `OPAQUE`.
pub(super) fn parse_alpha_mode(the_type: &str) -> GltfAlphaMode {
    if the_type.eq_ignore_ascii_case("OPAQUE") {
        GltfAlphaMode::Opaque
    } else if the_type.eq_ignore_ascii_case("MASK") {
        GltfAlphaMode::Mask
    } else if the_type.eq_ignore_ascii_case("BLEND") {
        GltfAlphaMode::Blend
    } else {
        GltfAlphaMode::Opaque
    }
}

/// `RWGltf_MaterialCommon` (`RWGltf_MaterialCommon.hxx:22-48`).
#[derive(Debug, Clone)]
pub(super) struct GltfMaterialCommon {
    pub ambient_texture: Option<ImageTexture>,
    pub diffuse_texture: Option<ImageTexture>,
    pub specular_texture: Option<ImageTexture>,
    pub id: String,
    pub name: String,
    pub ambient_color: LinearRgb,
    pub diffuse_color: LinearRgb,
    pub specular_color: LinearRgb,
    pub emissive_color: LinearRgb,
    pub shininess: f32,
    pub transparency: f32,
}

impl Default for GltfMaterialCommon {
    /// `RWGltf_MaterialCommon()` (`hxx:37-43`). The literals are sRGB and
    /// `Quantity_Color` stores linear components, so they are converted here.
    fn default() -> Self {
        Self {
            ambient_texture: None,
            diffuse_texture: None,
            specular_texture: None,
            id: String::new(),
            name: String::new(),
            ambient_color: LinearRgb::from_srgb(Color::new(0.1, 0.1, 0.1)),
            diffuse_color: LinearRgb::from_srgb(Color::new(0.8, 0.8, 0.8)),
            specular_color: LinearRgb::from_srgb(Color::new(0.2, 0.2, 0.2)),
            emissive_color: LinearRgb::new(0.0, 0.0, 0.0),
            shininess: 1.0,
            transparency: 0.0,
        }
    }
}

/// `RWGltf_MaterialMetallicRoughness` (`RWGltf_MaterialMetallicRoughness.hxx:34-64`).
#[derive(Debug, Clone)]
pub(super) struct GltfMaterialMetallicRoughness {
    pub base_color_texture: Option<ImageTexture>,
    pub metallic_roughness_texture: Option<ImageTexture>,
    pub emissive_texture: Option<ImageTexture>,
    pub occlusion_texture: Option<ImageTexture>,
    pub normal_texture: Option<ImageTexture>,
    pub id: String,
    pub name: String,
    pub base_color: LinearRgba,
    pub emissive_factor: [f32; 3],
    pub metallic: f32,
    pub roughness: f32,
    pub alpha_cut_off: f32,
    pub alpha_mode: GltfAlphaMode,
    pub is_double_sided: bool,
}

impl Default for GltfMaterialMetallicRoughness {
    /// `RWGltf_MaterialMetallicRoughness()` (`hxx:56-64`).
    fn default() -> Self {
        Self {
            base_color_texture: None,
            metallic_roughness_texture: None,
            emissive_texture: None,
            occlusion_texture: None,
            normal_texture: None,
            id: String::new(),
            name: String::new(),
            base_color: LinearRgba::WHITE,
            emissive_factor: [0.0, 0.0, 0.0],
            metallic: 1.0,
            roughness: 1.0,
            alpha_cut_off: 0.5,
            alpha_mode: GltfAlphaMode::Opaque,
            is_double_sided: false,
        }
    }
}

/// A `textures[]` entry (`texture.source`, `texture.target`).
pub(super) struct GltfTextureRef {
    pub(super) id: String,
    pub(super) source: Option<GltfKey>,
    pub(super) target: Option<f64>,
}

/// An `images[]` entry: a `uri` or an image slice inside a buffer.
pub(super) struct GltfImageDef {
    pub(super) id: String,
    pub(super) uri: Option<String>,
    pub(super) buffer_view: Option<usize>,
    /// glTF 1.0 `extensions.KHR_binary_glTF.bufferView`
    /// (`RWGltf_GltfJsonParser.cxx:1211-1220`).
    pub(super) khr_binary_buffer_view: Option<usize>,
}

/// A glTF reference that may be an index (2.0) or an id string (1.0).
#[derive(Debug, Clone, PartialEq)]
pub(super) enum GltfKey {
    Num(usize),
    Str(String),
}

impl GltfKey {
    pub(super) fn of(val: &JVal) -> Option<Self> {
        match val {
            JVal::Str(s) => Some(GltfKey::Str(s.clone())),
            _ => val.as_num().map(|n| GltfKey::Num(n as usize)),
        }
    }
}

/// Everything `gltfParseTexture` reads besides the texture reference itself.
pub(super) struct GltfTextureSource<'a> {
    pub(super) textures: &'a [GltfTextureRef],
    pub(super) images: &'a [GltfImageDef],
    pub(super) views: &'a [BufView],
    pub(super) buffers: &'a [Vec<u8>],
    /// `uri` of each buffer, parallel to `buffers` (empty for a GLB BIN chunk).
    pub(super) buffer_uris: &'a [String],
    pub(super) base_dir: &'a Path,
    pub(super) file_path: &'a str,
    /// `myBinBodyOffset` for a `.glb`, `None` for a plain `.gltf`.
    pub(super) bin_offset: Option<usize>,
    pub(super) is_gltf1: bool,
}

impl GltfTextureSource<'_> {
    fn texture(&self, key: &GltfKey) -> Option<&GltfTextureRef> {
        match key {
            GltfKey::Num(i) => self.textures.get(*i),
            GltfKey::Str(s) => self.textures.iter().find(|t| &t.id == s),
        }
    }

    fn image(&self, key: &GltfKey) -> Option<&GltfImageDef> {
        match key {
            GltfKey::Num(i) => self.images.get(*i),
            GltfKey::Str(s) => self.images.iter().find(|t| &t.id == s),
        }
    }

    /// `gltfParseTexturInGlbBuffer` (`cxx:1315-1356`): the image bytes sit in
    /// the GLB file's own BIN chunk, so they are addressed by file offset.
    fn texture_in_glb(&self, view_idx: usize) -> Option<ImageTexture> {
        let view = self.views.get(view_idx)?;
        if view.len == 0 {
            return None;
        }
        let off = self.bin_offset? + view.offset;
        Some(ImageTexture::new_file_part(self.file_path, off as i64, view.len as i64))
    }

    /// `gltfParseTextureInBufferView` (`cxx:1360-1455`): the image lives in a
    /// view of a regular buffer (data URI or external `.bin`).
    fn texture_in_buffer(&self, view_idx: usize, image_id: &str) -> Option<ImageTexture> {
        let view = self.views.get(view_idx)?;
        if view.len == 0 {
            return None;
        }
        let buf = self.buffers.get(view.buffer)?;
        if view.offset + view.len > buf.len() {
            return None;
        }
        let uri = self.buffer_uris.get(view.buffer).map(String::as_str).unwrap_or("");
        if let Some(rest) = uri.strip_prefix("data:") {
            let b64 = rest.split(";base64,").nth(1)?;
            let data = base64_decode(b64).ok()?;
            if view.offset + view.len > data.len() {
                return None;
            }
            return Some(ImageTexture::new_buffer(
                data[view.offset..view.offset + view.len].to_vec(),
                &format!("{}@{}", self.file_path, image_id),
            ));
        } else if uri.is_empty() {
            // A buffer without uri outside a GLB is rejected by the buffer pass.
            return None;
        }
        Some(ImageTexture::new_file_part(
            &self.base_dir.join(uri).to_string_lossy(),
            view.offset as i64,
            view.len as i64,
        ))
    }
}

/// `gltfParseTexture` (`RWGltf_GltfJsonParser.cxx:1164-1311`).
pub(super) fn parse_texture(key: &GltfKey, src: &GltfTextureSource) -> Option<ImageTexture> {
    let tex = src.texture(key)?;
    let source = tex.source.as_ref()?;
    // `cxx:1187-1191`: only GL_TEXTURE_2D is usable.
    if let Some(target) = tex.target {
        if target != 3553.0 {
            return None;
        }
    }
    let img = src.image(source)?;
    let image_id = match source {
        GltfKey::Num(i) => i.to_string(),
        GltfKey::Str(s) => s.clone(),
    };

    if src.bin_offset.is_some() {
        // `cxx:1204-1246`: 2.0 `bufferView`, else 1.0 `KHR_binary_glTF`.
        let mut view = img.buffer_view;
        if view.is_none() && src.is_gltf1 {
            view = img.khr_binary_buffer_view;
        }
        if let Some(v) = view {
            return src.texture_in_glb(v);
        }
    }

    let Some(uri) = img.uri.as_deref() else {
        // `cxx:1249-1269`: no uri, so the image must come from a bufferView.
        let view = img.buffer_view?;
        return src.texture_in_buffer(view, &image_id);
    };

    if let Some(rest) = uri.strip_prefix("data:") {
        // `cxx:1279-1303`: the base64 payload follows ";base64,".
        let b64 = rest.split(";base64,").nth(1)?;
        let data = base64_decode(b64).ok()?;
        return Some(ImageTexture::new_buffer(data, &format!("{}@{}", src.file_path, image_id)));
    }

    // `cxx:1306-1310`: a file next to the glTF.
    Some(ImageTexture::new_file(
        &src.base_dir.join(uri).to_string_lossy(),
    ))
}

/// `gltfReadVec4` (`RWGltf_GltfJsonParser.hxx:252-267`): exactly four numbers.
fn read_vec4(val: Option<&JVal>) -> Option<[f64; 4]> {
    let arr = val?.as_arr()?;
    if arr.len() != 4 {
        return None;
    }
    let mut out = [0.0; 4];
    for (i, v) in arr.iter().enumerate() {
        out[i] = v.as_num()?;
    }
    Some(out)
}

/// `gltfReadVec3` (`RWGltf_GltfJsonParser.hxx:278-294`): exactly three numbers.
fn read_vec3(val: Option<&JVal>) -> Option<[f64; 3]> {
    let arr = val?.as_arr()?;
    if arr.len() != 3 {
        return None;
    }
    let mut out = [0.0; 3];
    for (i, v) in arr.iter().enumerate() {
        out[i] = v.as_num()?;
    }
    Some(out)
}

/// `validateColor4` (`HXX:271-275`): every component within `[0, 1]`.
fn validate_color4(c: [f64; 4]) -> bool {
    c.iter().all(|v| (0.0..=1.0).contains(v))
}

/// `validateColor3` (`HXX:298-302`).
fn validate_color3(c: [f64; 3]) -> bool {
    c.iter().all(|v| (0.0..=1.0).contains(v))
}

/// `gltfParseStdMaterial` (`cxx:951-1017`): the glTF 1.0 `values` block.
pub(super) fn parse_std_material(
    node: &JVal,
    textures: &GltfTextureSource,
) -> Option<GltfMaterialCommon> {
    let values = node.get("values")?;
    let amb = values.get("ambient");
    let diff = values.get("diffuse");
    let emi = values.get("emission");
    let spec = values.get("specular");
    let shin = values.get("shininess");
    if amb.is_none() && diff.is_none() && emi.is_none() && spec.is_none() && shin.is_none() {
        return None;
    }

    let mut mat = GltfMaterialCommon::default();
    if let Some(v) = amb {
        if let Some(s) = v.as_str() {
            mat.ambient_texture = parse_texture(&GltfKey::Str(s.to_string()), textures);
        } else if let Some(c) = read_vec4(Some(v)).filter(|c| validate_color4(*c)) {
            mat.ambient_color = LinearRgb::from_srgb(Color::new(c[0] as f32, c[1] as f32, c[2] as f32));
        }
    }
    if let Some(v) = diff {
        if let Some(s) = v.as_str() {
            mat.diffuse_texture = parse_texture(&GltfKey::Str(s.to_string()), textures);
        } else if let Some(c) = read_vec4(Some(v)).filter(|c| validate_color4(*c)) {
            mat.diffuse_color = LinearRgb::from_srgb(Color::new(c[0] as f32, c[1] as f32, c[2] as f32));
            mat.transparency = (1.0 - c[3]) as f32;
        }
    }
    if let Some(c) = read_vec4(emi).filter(|c| validate_color4(*c)) {
        mat.emissive_color = LinearRgb::from_srgb(Color::new(c[0] as f32, c[1] as f32, c[2] as f32));
    }
    if let Some(v) = spec {
        if let Some(s) = v.as_str() {
            mat.specular_texture = parse_texture(&GltfKey::Str(s.to_string()), textures);
        }
    }
    if let Some(c) = read_vec4(spec).filter(|c| validate_color4(*c)) {
        mat.specular_color = LinearRgb::from_srgb(Color::new(c[0] as f32, c[1] as f32, c[2] as f32));
    }
    if let Some(n) = shin.and_then(|v| v.as_num()) {
        if n >= 0.0 {
            mat.shininess = (n / 1000.0).min(1.0) as f32;
        }
    }
    Some(mat)
}

/// `gltfParseCommonMaterial` (`cxx:1140-1158`): `extensions.KHR_materials_common`.
pub(super) fn parse_common_material(
    node: &JVal,
    textures: &GltfTextureSource,
) -> Option<GltfMaterialCommon> {
    let ext = node.get("extensions")?;
    let common = ext.get("KHR_materials_common")?;
    parse_std_material(common, textures)
}

/// `gltfParsePbrMaterial` (`cxx:1021-1137`). Returns `None` without a
/// `pbrMetallicRoughness` object — OCCT does not synthesise the spec defaults.
pub(super) fn parse_pbr_material(
    node: &JVal,
    textures: &GltfTextureSource,
) -> Option<GltfMaterialMetallicRoughness> {
    let mr = node.get("pbrMetallicRoughness")?;
    let mut mat = GltfMaterialMetallicRoughness::default();

    if let Some(v) = node.get("doubleSided").and_then(|v| v.as_bool()) {
        mat.is_double_sided = v;
    }
    if let Some(n) = node.get("alphaCutoff").and_then(|v| v.as_num()) {
        mat.alpha_cut_off = n as f32;
    }
    if let Some(s) = node.get("alphaMode").and_then(|v| v.as_str()) {
        mat.alpha_mode = parse_alpha_mode(s);
    }

    // `cxx:1074-1081`: baseColorTexture lives under pbrMetallicRoughness.
    if let Some(idx) = mr.get("baseColorTexture").and_then(|v| v.get("index")) {
        mat.base_color_texture = GltfKey::of(idx).and_then(|k| parse_texture(&k, textures));
    }
    if let Some(c) = read_vec4(mr.get("baseColorFactor")).filter(|c| validate_color4(*c)) {
        mat.base_color = LinearRgba::new(c[0] as f32, c[1] as f32, c[2] as f32, c[3] as f32);
    }
    if let Some(c) = read_vec3(node.get("emissiveFactor")).filter(|c| validate_color3(*c)) {
        mat.emissive_factor = [c[0] as f32, c[1] as f32, c[2] as f32];
    }
    // `cxx:1093-1100`.
    if let Some(idx) = mr.get("metallicRoughnessTexture").and_then(|v| v.get("index")) {
        mat.metallic_roughness_texture =
            GltfKey::of(idx).and_then(|k| parse_texture(&k, textures));
    }
    if let Some(n) = mr.get("metallicFactor").and_then(|v| v.as_num()) {
        mat.metallic = n as f32;
    }
    if let Some(n) = mr.get("roughnessFactor").and_then(|v| v.as_num()) {
        mat.roughness = n as f32;
    }
    // `cxx:1111-1135`.
    if let Some(idx) = node.get("normalTexture").and_then(|v| v.get("index")) {
        mat.normal_texture = GltfKey::of(idx).and_then(|k| parse_texture(&k, textures));
    }
    if let Some(idx) = node.get("emissiveTexture").and_then(|v| v.get("index")) {
        mat.emissive_texture = GltfKey::of(idx).and_then(|k| parse_texture(&k, textures));
    }
    if let Some(idx) = node.get("occlusionTexture").and_then(|v| v.get("index")) {
        mat.occlusion_texture = GltfKey::of(idx).and_then(|k| parse_texture(&k, textures));
    }
    Some(mat)
}

/// `gltfBindMaterial` (`cxx:867-947`): both halves land in one
/// `XCAFDoc_VisMaterial`, the PBR half providing the alpha mode and culling.
pub(super) fn bind_material(
    pbr: Option<GltfMaterialMetallicRoughness>,
    common: Option<GltfMaterialCommon>,
) -> Option<XcafDocVisMaterial> {
    if pbr.is_none() && common.is_none() {
        return None;
    }
    let mut mat = XcafDocVisMaterial::new();

    if let Some(c) = &common {
        let mut xde = XcafDocVisMaterialCommon::default();
        xde.is_defined = true;
        xde.ambient_color = c.ambient_color;
        xde.diffuse_color = c.diffuse_color;
        xde.specular_color = c.specular_color;
        xde.emissive_color = c.emissive_color;
        xde.shininess = c.shininess;
        xde.transparency = c.transparency;
        xde.diffuse_texture = c.diffuse_texture.clone();
        // `cxx:889-891`: fall back to the ambient texture.
        if xde.diffuse_texture.is_none() && c.ambient_texture.is_some() {
            xde.diffuse_texture = c.ambient_texture.clone();
        }
        mat.set_common_material(xde);
        if !c.name.is_empty() {
            mat.set_raw_name(Some(c.name.clone()));
        }
    }

    if let Some(p) = &pbr {
        let mut xde = XcafDocVisMaterialPbr::default();
        xde.is_defined = true;
        xde.metallic_roughness_texture = p.metallic_roughness_texture.clone();
        xde.base_color_texture = p.base_color_texture.clone();
        xde.emissive_texture = p.emissive_texture.clone();
        xde.occlusion_texture = p.occlusion_texture.clone();
        xde.normal_texture = p.normal_texture.clone();
        xde.base_color = p.base_color;
        xde.emissive_factor = p.emissive_factor;
        xde.metallic = p.metallic;
        xde.roughness = p.roughness;
        mat.set_pbr_material(xde);

        // `cxx:913-935`. The "transparent material with Opaque AlphaMode"
        // warning has no message stream in the port.
        let alpha = match p.alpha_mode {
            GltfAlphaMode::Opaque => AlphaMode::Opaque,
            GltfAlphaMode::Mask => AlphaMode::Mask,
            GltfAlphaMode::Blend => AlphaMode::Blend,
        };
        mat.set_alpha_mode(alpha, p.alpha_cut_off);
        // `cxx:936-939`: doubleSided -> Auto, else BackCulled.
        mat.set_face_culling(if p.is_double_sided {
            TypeOfBackfacingModel::Auto
        } else {
            TypeOfBackfacingModel::BackCulled
        });
        if !p.name.is_empty() {
            mat.set_raw_name(Some(p.name.clone()));
        }
    }
    Some(mat)
}

/// The `RWMesh` material view of a bound glTF material: the same values
/// `RWMesh_MaterialMap`/`RWObj_ObjMaterialMap` would read back out of the
/// `XCAFDoc_VisMaterial` (`XCAFDoc_VisMaterial::ConvertToCommonMaterial`,
/// `XCAFDoc_VisMaterial.cxx:158-186`).
pub(super) fn material_view(vis: &XcafDocVisMaterial) -> Material {
    let c = vis.convert_to_common_material();
    let (ar, ag, ab) = c.ambient_color.values_srgb();
    let (dr, dg, db) = c.diffuse_color.values_srgb();
    let (sr, sg, sb) = c.specular_color.values_srgb();
    Material {
        name: vis.raw_name().unwrap_or("default").to_string(),
        ambient: (ar, ag, ab),
        diffuse: (dr, dg, db),
        specular: (sr, sg, sb),
        shininess: c.shininess as f64,
        opacity: (1.0 - c.transparency) as f64,
        diffuse_texture: c.diffuse_texture.as_ref().map(|t| t.file_path().to_string()),
        specular_texture: None,
        bump_texture: None,
    }
}

/// `srgb_to_linear` re-export so the scene layer can build sRGB colours the
/// same way the material parser does.
pub(super) fn srgb_component(v: f64) -> f32 {
    srgb_to_linear(v)
}
