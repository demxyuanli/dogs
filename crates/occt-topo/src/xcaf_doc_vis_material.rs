//! `XCAFDoc_VisMaterial`, `XCAFDoc_VisMaterialCommon`, `XCAFDoc_VisMaterialPBR`
//! and the `Graphic3d_PBRMaterial` statics their conversions read.
//!
//! Source: `XCAFDoc_VisMaterial.hxx` / `.cxx`,
//! `XCAFDoc_VisMaterialCommon.hxx`, `XCAFDoc_VisMaterialPBR.hxx`,
//! `Graphic3d_PBRMaterial.hxx` / `.cxx`.
//!
//! `XCAFDoc_VisMaterial` is the attribute a document-level reader (glTF, XBF,
//! OBJ) writes for every material of the file; `XCAFPrs::CollectStyleSettings`
//! then hands it to the shapes through `XCAFPrs_Style::SetMaterial`, and the mesh
//! writers read it back through `XCAFPrs_Style::Material()`:
//!
//! * `RWObj_ObjMaterialMap::DefineMaterial` (`RWObj_ObjMaterialMap.cxx:85-92`)
//!   derives `Ka/Kd/Ks/Ns/Tr` from `ConvertToCommonMaterial()`.
//! * `RWMesh_MaterialMap::AddMaterial` (`RWMesh_MaterialMap.cxx:62-78`) uses the
//!   material's `TDataStd_Name` as the `newmtl` key.
//! * `RWMesh_ShapeIterator::initShape` (`RWMesh_ShapeIterator.cxx:145-152`) takes
//!   `BaseColor()` as the shape color.
//!
//! Colour storage: OCCT keeps `Quantity_Color` as a **linear** `myRgb` triple
//! (`Quantity_Color.hxx:36`), and every color these structs are constructed
//! with is linear (`Quantity_TOC_RGB`, `Quantity_Color.hxx:26-31`). The port
//! keeps the same linear triples in [`LinearRgb`] / [`LinearRgba`] and converts
//! to/from the crate's sRGB [`Color`] only at the API edge.
//!
//! **UNPORTED**: `FillMaterialAspect` / `FillAspect`
//! (`XCAFDoc_VisMaterial.cxx:222-350`) build `Graphic3d_MaterialAspect` and
//! `Graphic3d_Aspects` objects for the viewer; this crate has no visualization
//! aspects, and no mesh writer reads them. The attribute state they consume
//! (`AlphaMode`, `AlphaCutOff`, `FaceCulling`) is ported.

use std::sync::atomic::{AtomicU64, Ordering};

use occt_core::quantity::{Color, ColorRGBA};

use crate::image_texture::ImageTexture;

/// `Standard_ShortReal`'s `ShortRealEpsilon()` (`FLT_EPSILON`), the threshold
/// `XCAFDoc_VisMaterial.cxx:199` uses to decide whether a common material is
/// still opaque.
const SHORT_REAL_EPSILON: f32 = f32::EPSILON;

/// `Quantity_Color::Convert_sRGB_To_LinearRGB` — the **double** overload
/// (`Quantity_Color.hxx:311-316`); `SetValues(..., Quantity_TOC_sRGB)`
/// (`Quantity_Color.cxx:345-351`) truncates its result to the stored `float`.
pub(crate) fn srgb_to_linear(v: f64) -> f32 {
    let linear = if v <= 0.04045 {
        v / 12.92
    } else {
        ((v + 0.055) / 1.055).powf(2.4)
    };
    linear as f32
}

/// `Quantity_Color::Convert_LinearRGB_To_sRGB` — the **double** overload
/// (`Quantity_Color.hxx:297-300`), used by `Values(..., Quantity_TOC_sRGB)`
/// (`Quantity_Color.cxx:499-504`).
pub(crate) fn linear_to_srgb(v: f64) -> f64 {
    if v <= 0.0031308 {
        v * 12.92
    } else {
        v.powf(1.0 / 2.4) * 1.055 - 0.055
    }
}

/// `Quantity_Color`'s internal linear `myRgb` (`Quantity_Color.hxx:36`):
/// "other OCCT interfaces ... would expect them in linear space".
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LinearRgb {
    pub r: f32,
    pub g: f32,
    pub b: f32,
}

impl LinearRgb {
    /// `Quantity_Color(const NCollection_Vec3<float>&)` — "Define color from
    /// linear RGB values" (`Quantity_Color.hxx:63-66`).
    pub const fn new(r: f32, g: f32, b: f32) -> Self {
        Self { r, g, b }
    }

    /// `Quantity_Color::SetValues(..., Quantity_TOC_sRGB)`.
    pub fn from_srgb(color: Color) -> Self {
        Self::new(
            srgb_to_linear(color.r as f64),
            srgb_to_linear(color.g as f64),
            srgb_to_linear(color.b as f64),
        )
    }

    /// `Quantity_Color::Values(..., Quantity_TOC_sRGB)`.
    pub fn to_srgb(self) -> Color {
        let (r, g, b) = self.values_srgb();
        Color::new(r as f32, g as f32, b as f32)
    }

    /// The same conversion, kept in `f64` for the `.mtl` formatter.
    pub fn values_srgb(self) -> (f64, f64, f64) {
        (
            linear_to_srgb(self.r as f64),
            linear_to_srgb(self.g as f64),
            linear_to_srgb(self.b as f64),
        )
    }

    /// `(NCollection_Vec3<float>)theColor * theFactor`.
    pub fn scaled(self, factor: f32) -> Self {
        Self::new(self.r * factor, self.g * factor, self.b * factor)
    }

    /// `NCollection_Vec3<float>::maxComp()`.
    pub fn max_comp(self) -> f32 {
        self.r.max(self.g).max(self.b)
    }

    /// `NCollection_Vec3<float>::cwiseMin`.
    pub fn cwise_min(self, other: LinearRgb) -> Self {
        Self::new(self.r.min(other.r), self.g.min(other.g), self.b.min(other.b))
    }
}

/// `Quantity_ColorRGBA`'s internal storage: a linear RGB triple plus alpha.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LinearRgba {
    pub rgb: LinearRgb,
    pub a: f32,
}

impl LinearRgba {
    pub const fn new(r: f32, g: f32, b: f32, a: f32) -> Self {
        Self { rgb: LinearRgb::new(r, g, b), a }
    }

    /// `Quantity_NOC_WHITE` with alpha 1.
    pub const WHITE: Self = Self::new(1.0, 1.0, 1.0, 1.0);

    pub fn from_srgb(color: ColorRGBA) -> Self {
        Self { rgb: LinearRgb::from_srgb(Color::new(color.r, color.g, color.b)), a: color.a }
    }

    pub fn to_srgb(self) -> ColorRGBA {
        let c = self.rgb.to_srgb();
        ColorRGBA { r: c.r, g: c.g, b: c.b, a: self.a }
    }
}

/// `Graphic3d_AlphaMode` (`Graphic3d_AlphaMode.hxx:18-33`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlphaMode {
    /// `Graphic3d_AlphaMode_Opaque`.
    Opaque,
    /// `Graphic3d_AlphaMode_Mask`.
    Mask,
    /// `Graphic3d_AlphaMode_Blend`.
    Blend,
    /// `Graphic3d_AlphaMode_MaskBlend`.
    MaskBlend,
    /// `Graphic3d_AlphaMode_BlendAuto` — the attribute's default.
    BlendAuto,
}

/// `Graphic3d_TypeOfBackfacingModel`
/// (`Graphic3d_TypeOfBackfacingModel.hxx:21-36`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TypeOfBackfacingModel {
    /// `Graphic3d_TypeOfBackfacingModel_Auto` — the attribute's default.
    Auto,
    /// `Graphic3d_TypeOfBackfacingModel_DoubleSided`.
    DoubleSided,
    /// `Graphic3d_TypeOfBackfacingModel_BackCulled`.
    BackCulled,
    /// `Graphic3d_TypeOfBackfacingModel_FrontCulled`.
    FrontCulled,
}

/// `equal_to<handle<Image_Texture>>` (`Image_Texture.hxx:105-118`): two handles
/// are equal when their `TextureId()` strings match, not when their paths or
/// buffers do.
fn textures_equal(a: &Option<ImageTexture>, b: &Option<ImageTexture>) -> bool {
    match (a, b) {
        (None, None) => true,
        (Some(a), Some(b)) => a.texture_id() == b.texture_id(),
        _ => false,
    }
}

/// `XCAFDoc_VisMaterialCommon` (`XCAFDoc_VisMaterialCommon.hxx:29-104`): the
/// obsolete (OpenGL fixed-function) material definition.
#[derive(Debug, Clone)]
pub struct XcafDocVisMaterialCommon {
    /// `DiffuseTexture` — image defining the diffuse color.
    pub diffuse_texture: Option<ImageTexture>,
    /// `AmbientColor` (linear). Default `(0.1, 0.1, 0.1)`.
    pub ambient_color: LinearRgb,
    /// `DiffuseColor` (linear). Default `(0.8, 0.8, 0.8)`.
    pub diffuse_color: LinearRgb,
    /// `SpecularColor` (linear). Default `(0.2, 0.2, 0.2)`.
    pub specular_color: LinearRgb,
    /// `EmissiveColor` (linear). Default `(0.0, 0.0, 0.0)`.
    pub emissive_color: LinearRgb,
    /// `Shininess`. Default 1.
    pub shininess: f32,
    /// `Transparency` within `[0, 1]`, 0 meaning opaque.
    pub transparency: f32,
    /// `IsDefined` — TRUE by default here, FALSE for the attribute's own
    /// undefined state (`XCAFDoc_VisMaterial.cxx:45-46`).
    pub is_defined: bool,
}

impl Default for XcafDocVisMaterialCommon {
    /// `XCAFDoc_VisMaterialCommon::XCAFDoc_VisMaterialCommon()`
    /// (`XCAFDoc_VisMaterialCommon.hxx:48-57`). The constructor passes
    /// `Quantity_TOC_RGB`, i.e. the literals are linear.
    fn default() -> Self {
        Self {
            diffuse_texture: None,
            ambient_color: LinearRgb::new(0.1, 0.1, 0.1),
            diffuse_color: LinearRgb::new(0.8, 0.8, 0.8),
            specular_color: LinearRgb::new(0.2, 0.2, 0.2),
            emissive_color: LinearRgb::new(0.0, 0.0, 0.0),
            shininess: 1.0,
            transparency: 0.0,
            is_defined: true,
        }
    }
}

impl XcafDocVisMaterialCommon {
    /// The state `UnsetCommonMaterial` leaves behind
    /// (`XCAFDoc_VisMaterial.hxx:126-132`): default values, `IsDefined = false`.
    pub fn undefined() -> Self {
        Self { is_defined: false, ..Self::default() }
    }

    /// `XCAFDoc_VisMaterialCommon::IsEqual`
    /// (`XCAFDoc_VisMaterialCommon.hxx:70-96`).
    pub fn is_equal(&self, other: &Self) -> bool {
        if std::ptr::eq(self, other) {
            return true;
        } else if self.is_defined != other.is_defined {
            return false;
        } else if !self.is_defined {
            return true;
        }
        textures_equal(&self.diffuse_texture, &other.diffuse_texture)
            && self.ambient_color == other.ambient_color
            && self.diffuse_color == other.diffuse_color
            && self.specular_color == other.specular_color
            && self.emissive_color == other.emissive_color
            && self.shininess == other.shininess
            && self.transparency == other.transparency
    }
}

/// `XCAFDoc_VisMaterialPBR` (`XCAFDoc_VisMaterialPBR.hxx:29-110`): the
/// metallic-roughness material definition.
#[derive(Debug, Clone)]
pub struct XcafDocVisMaterialPbr {
    /// `BaseColorTexture` — RGB texture for the base color.
    pub base_color_texture: Option<ImageTexture>,
    /// `MetallicRoughnessTexture` — RG texture packing metallic and roughness.
    pub metallic_roughness_texture: Option<ImageTexture>,
    /// `EmissiveTexture`.
    pub emissive_texture: Option<ImageTexture>,
    /// `OcclusionTexture`.
    pub occlusion_texture: Option<ImageTexture>,
    /// `NormalTexture`.
    pub normal_texture: Option<ImageTexture>,
    /// `BaseColor` — base color (or scale factor for the texture); linear,
    /// `[1, 1, 1, 1]` by default.
    pub base_color: LinearRgba,
    /// `EmissiveFactor` (linear); `[0, 0, 0]` by default.
    pub emissive_factor: [f32; 3],
    /// `Metallic` in `[0, 1]`; 1 by default.
    pub metallic: f32,
    /// `Roughness` in `[0, 1]`; 1 by default.
    pub roughness: f32,
    /// `RefractionIndex` (IOR) in `[1, 3]`; 1.5 by default.
    pub refraction_index: f32,
    /// `IsDefined` — TRUE by default here, FALSE for the attribute's own
    /// undefined state (`XCAFDoc_VisMaterial.cxx:44`).
    pub is_defined: bool,
}

impl Default for XcafDocVisMaterialPbr {
    /// `XCAFDoc_VisMaterialPBR::XCAFDoc_VisMaterialPBR()`
    /// (`XCAFDoc_VisMaterialPBR.hxx:47-58`).
    fn default() -> Self {
        Self {
            base_color_texture: None,
            metallic_roughness_texture: None,
            emissive_texture: None,
            occlusion_texture: None,
            normal_texture: None,
            base_color: LinearRgba::WHITE,
            emissive_factor: [0.0, 0.0, 0.0],
            metallic: 1.0,
            roughness: 1.0,
            refraction_index: 1.5,
            is_defined: true,
        }
    }
}

impl XcafDocVisMaterialPbr {
    /// The state `UnsetPbrMaterial` leaves behind
    /// (`XCAFDoc_VisMaterial.hxx:114-120`): default values, `IsDefined = false`.
    pub fn undefined() -> Self {
        Self { is_defined: false, ..Self::default() }
    }

    /// `XCAFDoc_VisMaterialPBR::IsEqual`
    /// (`XCAFDoc_VisMaterialPBR.hxx:69-95`).
    pub fn is_equal(&self, other: &Self) -> bool {
        if std::ptr::eq(self, other) {
            return true;
        } else if self.is_defined != other.is_defined {
            return false;
        } else if !self.is_defined {
            return true;
        }
        textures_equal(&self.base_color_texture, &other.base_color_texture)
            && textures_equal(&self.metallic_roughness_texture, &other.metallic_roughness_texture)
            && textures_equal(&self.emissive_texture, &other.emissive_texture)
            && textures_equal(&self.occlusion_texture, &other.occlusion_texture)
            && textures_equal(&self.normal_texture, &other.normal_texture)
            && self.base_color == other.base_color
            && self.emissive_factor == other.emissive_factor
            && self.metallic == other.metallic
            && self.roughness == other.roughness
            && self.refraction_index == other.refraction_index
    }
}

/// `Graphic3d_PBRMaterial::MetallicFromSpecular`
/// (`Graphic3d_PBRMaterial.hxx:124-127`): the largest specular component.
fn metallic_from_specular(specular: LinearRgb) -> f32 {
    specular.max_comp()
}

/// `Graphic3d_PBRMaterial::RoughnessFromSpecular`
/// (`Graphic3d_PBRMaterial.cxx:23-36`).
fn roughness_from_specular(specular: LinearRgb, shininess: f64) -> f32 {
    let mut roughness_factor = 1.0 - shininess;
    // `theSpecular.Red() * 0.2125 + Green * 0.7154 + Blue * 0.0721` reads the
    // linear components, as `Red()/Green()/Blue()` do (`Quantity_Color.hxx`).
    let spec_intens = specular.r as f64 * 0.2125
        + specular.g as f64 * 0.7154
        + specular.b as f64 * 0.0721;
    if spec_intens < 0.1 {
        // low specular intensity should produce a rough material even if
        // shininess is high
        roughness_factor *= 1.0 - spec_intens;
    }
    roughness_factor as f32
}

/// Identity of the material handle. OCCT material references are
/// `handle<XCAFDoc_VisMaterial>`, and `XCAFPrs_Style::IsEqual`
/// (`XCAFPrs_Style.hxx:119`) compares those handles by pointer, so two
/// *separately created* materials with identical values must remain distinct
/// table entries. A clone keeps the id, exactly like a copied handle.
static NEXT_MATERIAL_ID: AtomicU64 = AtomicU64::new(1);

fn next_material_id() -> u64 {
    NEXT_MATERIAL_ID.fetch_add(1, Ordering::Relaxed)
}

/// `XCAFDoc_VisMaterial` (`XCAFDoc_VisMaterial.hxx:52-220`), minus the
/// `TDF_Attribute` plumbing (`Restore` / `NewEmpty` / `Paste` / `DumpJson`,
/// `cxx:79-121`, `cxx:355-380`), which has no document framework here.
#[derive(Debug, Clone)]
pub struct XcafDocVisMaterial {
    /// Handle identity — see `NEXT_MATERIAL_ID`.
    id: u64,
    /// `myRawName` — transient material name / tag, not stored in the document.
    raw_name: Option<String>,
    /// `TDataStd_Name` on the material's table label, which
    /// `RWMesh_MaterialMap::AddMaterial` reads as the `newmtl` key
    /// (`RWMesh_MaterialMap.cxx:67-77`). `None` when the material is not in a
    /// material table, or when it was added without a name.
    label_name: Option<String>,
    /// `myPbrMat`.
    pbr: XcafDocVisMaterialPbr,
    /// `myCommonMat`.
    common: XcafDocVisMaterialCommon,
    /// `myAlphaMode`.
    alpha_mode: AlphaMode,
    /// `myAlphaCutOff`.
    alpha_cut_off: f32,
    /// `myFaceCulling`.
    face_culling: TypeOfBackfacingModel,
}

impl Default for XcafDocVisMaterial {
    /// `XCAFDoc_VisMaterial::XCAFDoc_VisMaterial()`
    /// (`XCAFDoc_VisMaterial.cxx:38-47`): both definitions undefined, alpha
    /// `BlendAuto` with cutoff 0.5, face culling `Auto`.
    fn default() -> Self {
        Self {
            id: next_material_id(),
            raw_name: None,
            label_name: None,
            pbr: XcafDocVisMaterialPbr::undefined(),
            common: XcafDocVisMaterialCommon::undefined(),
            alpha_mode: AlphaMode::BlendAuto,
            alpha_cut_off: 0.5,
            face_culling: TypeOfBackfacingModel::Auto,
        }
    }
}

impl XcafDocVisMaterial {
    /// Empty constructor (`cxx:38-47`).
    pub fn new() -> Self {
        Self::default()
    }

    /// Handle identity for `XCAFPrs_Style::IsEqual` — clones share it.
    pub fn id(&self) -> u64 {
        self.id
    }

    /// `IsEmpty()` (`hxx:75`): neither definition is set.
    pub fn is_empty(&self) -> bool {
        !self.pbr.is_defined && !self.common.is_defined
    }

    /// `RawName()` (`hxx:180`).
    pub fn raw_name(&self) -> Option<&str> {
        self.raw_name.as_deref()
    }

    /// `SetRawName()` (`hxx:185`).
    pub fn set_raw_name(&mut self, name: Option<String>) {
        self.raw_name = name;
    }

    /// The `TDataStd_Name` of the material's table label
    /// (`XCAFDoc_VisMaterialTool::AddMaterial`, `cxx:83-86`).
    pub fn label_name(&self) -> Option<&str> {
        self.label_name.as_deref()
    }

    /// Assign the label name; called by the material table when the material is
    /// added (`XCAFDoc_VisMaterialTool.cxx:82-87`).
    pub(crate) fn set_label_name(&mut self, name: Option<String>) {
        self.label_name = name;
    }

    /// `HasPbrMaterial()` (`hxx:88`).
    pub fn has_pbr_material(&self) -> bool {
        self.pbr.is_defined
    }

    /// `PbrMaterial()` (`hxx:92`).
    pub fn pbr_material(&self) -> &XcafDocVisMaterialPbr {
        &self.pbr
    }

    /// `SetPbrMaterial()` (`cxx:55-59`).
    pub fn set_pbr_material(&mut self, material: XcafDocVisMaterialPbr) {
        self.pbr = material;
    }

    /// `UnsetPbrMaterial()` (`hxx:114-120`).
    pub fn unset_pbr_material(&mut self) {
        self.pbr = XcafDocVisMaterialPbr::undefined();
    }

    /// `HasCommonMaterial()` (`hxx:139`).
    pub fn has_common_material(&self) -> bool {
        self.common.is_defined
    }

    /// `CommonMaterial()` (`hxx:143`).
    pub fn common_material(&self) -> &XcafDocVisMaterialCommon {
        &self.common
    }

    /// `SetCommonMaterial()` (`cxx:63-67`).
    pub fn set_common_material(&mut self, material: XcafDocVisMaterialCommon) {
        self.common = material;
    }

    /// `UnsetCommonMaterial()` (`hxx:126-132`).
    pub fn unset_common_material(&mut self) {
        self.common = XcafDocVisMaterialCommon::undefined();
    }

    /// `BaseColor()` (`cxx:149-166`): the PBR base color, else the common
    /// diffuse color with `1 - Transparency` alpha, else opaque white. The
    /// result is linear, like `Quantity_ColorRGBA`'s internals.
    pub fn base_color(&self) -> LinearRgba {
        if self.pbr.is_defined {
            self.pbr.base_color
        } else if self.common.is_defined {
            LinearRgba { rgb: self.common.diffuse_color, a: 1.0 - self.common.transparency }
        } else {
            LinearRgba::WHITE
        }
    }

    /// `AlphaMode()` (`hxx:147`).
    pub fn alpha_mode(&self) -> AlphaMode {
        self.alpha_mode
    }

    /// `AlphaCutOff()` (`hxx:151`).
    pub fn alpha_cut_off(&self) -> f32 {
        self.alpha_cut_off
    }

    /// `SetAlphaMode()` (`cxx:69-74`).
    pub fn set_alpha_mode(&mut self, mode: AlphaMode, cut_off: f32) {
        self.alpha_mode = mode;
        self.alpha_cut_off = cut_off;
    }

    /// `FaceCulling()` (`hxx:160`).
    pub fn face_culling(&self) -> TypeOfBackfacingModel {
        self.face_culling
    }

    /// `SetFaceCulling()` (`cxx:76-80`).
    pub fn set_face_culling(&mut self, culling: TypeOfBackfacingModel) {
        self.face_culling = culling;
    }

    /// `IsEqual` (`hxx:196-208`): deep comparison by actual values.
    /// `XCAFPrs_Style` compares handles instead — see [`Self::id`].
    pub fn is_equal(&self, other: &Self) -> bool {
        if self.id == other.id {
            return true;
        }
        self.face_culling == other.face_culling
            && self.alpha_cut_off == other.alpha_cut_off
            && self.alpha_mode == other.alpha_mode
            && self.common.is_equal(&other.common)
            && self.pbr.is_equal(&other.pbr)
    }

    /// `ConvertToCommonMaterial()` (`cxx:158-186`).
    pub fn convert_to_common_material(&self) -> XcafDocVisMaterialCommon {
        if self.common.is_defined {
            return self.common.clone();
        } else if !self.pbr.is_defined {
            // `return XCAFDoc_VisMaterialCommon();` — the *defined* default.
            return XcafDocVisMaterialCommon::default();
        }

        // convert metal-roughness into common
        let rgb = self.pbr.base_color.rgb;
        let mut common = XcafDocVisMaterialCommon::default();
        common.is_defined = true;
        common.diffuse_texture = self.pbr.base_color_texture.clone();
        common.diffuse_color = rgb;
        common.specular_color = rgb.scaled(self.pbr.metallic);
        common.transparency = 1.0 - self.pbr.base_color.a;
        common.shininess = 1.0 - self.pbr.roughness;
        if self.pbr.emissive_texture.is_none() {
            let e = self.pbr.emissive_factor;
            common.emissive_color = LinearRgb::new(e[0], e[1], e[2])
                .cwise_min(LinearRgb::new(1.0, 1.0, 1.0));
        }
        common
    }

    /// `ConvertToPbrMaterial()` (`cxx:190-211`).
    pub fn convert_to_pbr_material(&self) -> XcafDocVisMaterialPbr {
        if self.pbr.is_defined {
            return self.pbr.clone();
        } else if !self.common.is_defined {
            // `return XCAFDoc_VisMaterialPBR();` — the *defined* default.
            return XcafDocVisMaterialPbr::default();
        }

        let mut pbr = XcafDocVisMaterialPbr::default();
        pbr.is_defined = true;
        pbr.base_color_texture = self.common.diffuse_texture.clone();
        pbr.base_color = LinearRgba {
            rgb: self.common.diffuse_color,
            a: 1.0 - self.common.transparency,
        };
        pbr.metallic = if self.common.transparency <= SHORT_REAL_EPSILON {
            metallic_from_specular(self.common.specular_color)
        } else {
            0.0
        };
        pbr.roughness =
            roughness_from_specular(self.common.specular_color, self.common.shininess as f64);
        pbr.emissive_factor =
            [self.common.emissive_color.r, self.common.emissive_color.g, self.common.emissive_color.b];
        pbr
    }
}
