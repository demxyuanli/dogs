//! Material sidecar writers: `RWMesh_MaterialMap` and `RWObj_ObjMaterialMap`,
//! plus the `Quantity_Color` colour-space conversions the `.mtl` numbers go
//! through.
//!
//! Source: `RWMesh_MaterialMap.cxx`, `RWObj_ObjMaterialMap.cxx`. Split out of
//! `obj_caf_writer/mod.rs` unchanged (2026-10-17); see that module for the export
//! pipeline and the "not ported" list.

use crate::image_texture::ImageTexture;
use crate::xcaf_doc_vis_material::{LinearRgb, XcafDocVisMaterialCommon};
use crate::xcaf_prs::XcafPrsStyle;

// ---------------------------------------------------------------------------
// Quantity_Color color-space conversions
// ---------------------------------------------------------------------------
//
// `linear_to_srgb` / `srgb_to_linear` and the linear `LinearRgb` triple now live
// in `xcaf_doc_vis_material` (`Quantity_Color`'s internal storage), because the
// material structs are constructed with linear values too.

// ---------------------------------------------------------------------------
// RWMesh_MaterialMap
// ---------------------------------------------------------------------------

/// `RWMesh_MaterialMap` (`RWMesh_MaterialMap.hxx:22-95`): registers materials,
/// hands out the `mat_<n>` keys used by `newmtl` / `usemtl`, and copies the
/// styles' textures next to the exported model.
///
/// Faithful to the base class: key / name generation (`AddMaterial`,
/// `cxx:53-116`), [Self::create_texture_folder] (`cxx:250-306`) and
/// [Self::copy_texture] (`cxx:200-232`). The abstract `DefineMaterial` lives on
/// [`ObjMaterialMap`].
pub struct MeshMaterialMap {
    /// `myFolder` (`hxx:79`) — the output model's folder, separator included
    /// (`OSD_Path::FolderAndFileFromPath` keeps it, `OSD_Path.cxx:1818`), or `"."`
    /// when the model path names no folder (`cxx:36-39`).
    folder: String,
    /// `myTexFolder` (`hxx:80`) — full path of the textures folder, empty until
    /// [Self::create_texture_folder].
    tex_folder: String,
    /// `myTexFolderShort` (`hxx:81`) — the same folder relative to the model.
    tex_folder_short: String,
    /// `myFileName` (`hxx:85`).
    file_name: String,
    /// `myShortFileNameBase` (`hxx:82`) — model file name without extension.
    short_file_name_base: String,
    /// `myKeyPrefix` — `"mat_"` (`RWMesh_MaterialMap.cxx:31`).
    key_prefix: String,
    /// `myNbMaterials` (`RWMesh_MaterialMap.hxx:93`).
    nb_materials: i32,
    /// `myStyles` — a `NCollection_DoubleMap<XCAFPrs_Style, AsciiString>` used
    /// as style -> key and key -> style in `AddMaterial` (`cxx:54-110`).
    styles: Vec<(XcafPrsStyle, String)>,
    /// `myImageMap` (`hxx:90`) — texture id -> already-copied relative path.
    image_map: Vec<(String, String)>,
    /// `myImageFailMap` (`hxx:91`) — texture ids that failed to copy, never
    /// retried (`RWObj_ObjMaterialMap.cxx:128`).
    image_fail_map: Vec<String>,
    /// `myIsFailed` (`hxx:94`).
    is_failed: bool,
    /// `myMatNameAsKey` (`hxx:95`), TRUE by default.
    mat_name_as_key: bool,
}

impl MeshMaterialMap {
    /// `RWMesh_MaterialMap::RWMesh_MaterialMap` (`RWMesh_MaterialMap.cxx:29-40`):
    /// split the model path into folder + stem, defaulting the folder to `"."`.
    pub(crate) fn new(file_name: &str) -> Self {
        let (mut folder, a_file_name) = crate::image_texture::ImageTexture::folder_and_file(file_name);
        let (short_file_name_base, _ext) =
            crate::image_texture::ImageTexture::file_name_and_extension(&a_file_name);
        if folder.is_empty() {
            // `if (myFolder.IsEmpty()) { myFolder = "."; }` (`cxx:36-39`).
            folder = ".".to_string();
        }
        Self {
            folder,
            tex_folder: String::new(),
            tex_folder_short: String::new(),
            file_name: file_name.to_string(),
            short_file_name_base,
            key_prefix: "mat_".to_string(),
            nb_materials: 0,
            styles: Vec::new(),
            image_map: Vec::new(),
            image_fail_map: Vec::new(),
            is_failed: false,
            mat_name_as_key: true,
        }
    }

    /// `RWMesh_MaterialMap::FindMaterial` (`RWMesh_MaterialMap.hxx:41-48`).
    /// The `NCollection_DoubleMap` lookup uses `XCAFPrs_Style::IsEqual`
    /// (`RWMesh_MaterialMap.cxx:103`), i.e. [`XcafPrsStyle::is_equal`].
    pub(crate) fn find_material(&self, style: &XcafPrsStyle) -> Option<String> {
        self.styles.iter().find(|(s, _)| s.is_equal(style)).map(|(_, k)| k.clone())
    }

    /// `RWMesh_MaterialMap::AddMaterial` (`RWMesh_MaterialMap.cxx:53-116`):
    /// the key / name generation, then the `DefineMaterial` call. The port
    /// turns the virtual callback into a return value the caller writes with.
    ///
    /// Returns `(key, name, is_new)`; `key` is the `newmtl` name, `name` the
    /// pretty name (unused by the OBJ writer, `RWObj_ObjMaterialMap.cxx:76`).
    pub(crate) fn add_material(&mut self, style: &XcafPrsStyle) -> (String, String, bool) {
        if let Some(key) = self.find_material(style) {
            return (key, String::new(), false);
        }

        let mut mat_key;
        let mat_name_suffix;
        let mut mat_name;
        let mut counter: i32 = 0;
        // `int* aCounterPtr = &myNbMaterials;` (`cxx:60`) — with `myMatNameAsKey`
        // and a material present the counter points at the local `aCounter`
        // (`cxx:64`).
        let mut counter_is_local = false;
        if self.mat_name_as_key {
            // `if (!theStyle.Material().IsNull() && !theStyle.Material()
            // ->IsEmpty())` (`cxx:62-78`): the material's own label name becomes
            // the key, falling back to `"mat"` when the table label carries no
            // `TDataStd_Name`.
            let material = style.material().filter(|m| !m.is_empty());
            if let Some(material) = material {
                counter_is_local = true;
                let material_name = material.label_name().unwrap_or("mat").to_string();
                mat_name_suffix = material_name.clone();
                mat_name = material_name;
            } else {
                self.nb_materials += 1;
                mat_name_suffix = self.key_prefix.clone();
                mat_name = format!("{}{}", mat_name_suffix, self.nb_materials);
            }
            mat_key = mat_name.clone();
        } else {
            // `else` branch (`cxx:85-91`) — unused with `myMatNameAsKey`.
            mat_key = self.nb_materials.to_string();
            self.nb_materials += 1;
            mat_name_suffix = self.key_prefix.clone();
            mat_name = format!("{}{}", mat_name_suffix, mat_key);
        }

        loop {
            // `for (;; ++(*aCounterPtr))` (`cxx:93-111`).
            let occupied = self.styles.iter().any(|(_, k)| *k == mat_key);
            if !occupied {
                break;
            }
            if self.mat_name_as_key {
                let n = if counter_is_local { counter } else { self.nb_materials };
                mat_name = format!("{mat_name_suffix}{n}");
                mat_key = mat_name.clone();
                if counter_is_local {
                    counter += 1;
                } else {
                    self.nb_materials += 1;
                }
            } else {
                mat_key = self.nb_materials.to_string();
                self.nb_materials += 1;
                mat_name = format!("{mat_name_suffix}{mat_key}");
            }
        }

        self.styles.push((style.clone(), mat_key.clone()));
        (mat_key, mat_name, true)
    }

    /// `RWMesh_MaterialMap::CreateTextureFolder` (`RWMesh_MaterialMap.cxx:250-306`).
    ///
    /// Lays the folder out as `<model folder>/<model stem>_textures/`
    /// (`cxx:257-259`) and creates it once. When the model folder itself does not
    /// exist the call fails (`cxx:270-275`); when creating the textures folder
    /// fails it falls back to the model folder and clears the short path
    /// (`cxx:293-300`).
    ///
    /// **UNPORTED by design**: the `OSD_Protection` copying (`cxx:277-291`) has
    /// no Rust counterpart and no effect on exported data — `OSD_Protection` is
    /// the read/write/execute permission triplet `OSD_Directory::Build` hands to
    /// the OS; `std::fs::create_dir` applies the process umask and exposes no
    /// equivalent setter. Reproducing it would mean hand-rolling `chmod`/ACL
    /// calls for a permission bit that neither the OBJ/MTL text nor the copied
    /// image bytes depend on.
    pub fn create_texture_folder(&mut self) -> bool {
        if !self.tex_folder.is_empty() {
            return true;
        }

        // `myTexFolderShort = myShortFileNameBase + "_textures/";` (`cxx:257`).
        self.tex_folder_short = format!("{}_textures/", self.short_file_name_base);
        // `myTexFolder = myFolder + "/" + myTexFolderShort;` (`cxx:258`).
        self.tex_folder = format!("{}/{}", self.folder, self.tex_folder_short);
        if crate::image_texture::ImageTexture::is_dir(&self.tex_folder) {
            return true;
        }

        // `if (!aResDir.Exists()) { Message::SendFail() << "Failed to create
        // textures folder '" << myFolder << "'"; return false; }` (`cxx:270-275`).
        if !crate::image_texture::ImageTexture::is_dir(&self.folder) {
            return false;
        }

        // `aTexDir.Build(aProt); if (aTexDir.Failed())` (`cxx:292-300`).
        if std::fs::create_dir(&self.tex_folder).is_err() {
            // fallback to the same folder as output model file
            self.tex_folder = self.folder.clone();
            self.tex_folder_short.clear();
        }
        true
    }

    /// `RWMesh_MaterialMap::CopyTexture` (`RWMesh_MaterialMap.cxx:200-232`):
    /// copy the image into the textures folder, renamed after `key` unless
    /// the source is a plain image file (then its own name is kept,
    /// `cxx:206-215`).
    ///
    /// Returns `(result texture path relative to the model, success)`. The path
    /// is assigned on every branch *before* the copy is attempted (as in the
    /// C++, `cxx:212-213` / `cxx:227-228`), so a failed copy still yields the
    /// name the caller would reference — `RWObj_ObjMaterialMap.cxx:137-140`
    /// writes `map_Kd` for any non-empty result.
    pub fn copy_texture(&mut self, texture: &ImageTexture, key: &str) -> (String, bool) {
        self.create_texture_folder();

        let mut tex_file_name;
        let texture_src = texture.file_path();
        if !texture_src.is_empty() && texture.file_offset() <= 0 && texture.file_length() <= 0 {
            // `OSD_Path::FolderAndFileFromPath(aTextureSrc, aSrcTexFolder,
            // aTexFileName)` (`cxx:209-210`).
            let (_src_folder, file_name) =
                crate::image_texture::ImageTexture::folder_and_file(texture_src);
            tex_file_name = file_name;
            let res_tex_file = format!("{}{}", self.tex_folder, tex_file_name);
            let res_texture = format!("{}{}", self.tex_folder_short, tex_file_name);
            let ok = Self::copy_file_to(texture_src, &res_tex_file);
            return (res_texture, ok);
        }

        // `TCollection_AsciiString anExt = theTexture->ProbeImageFileFormat();
        // if (anExt.IsEmpty()) { anExt = "bin"; }` (`cxx:217-221`).
        let mut ext = texture.probe_image_file_format();
        if ext.is_empty() {
            ext = "bin".to_string();
        }
        // `aTexFileName = theKey + "." + anExt;` (`cxx:222`).
        tex_file_name = format!("{key}.{ext}");

        let res_tex_file = format!("{}{}", self.tex_folder, tex_file_name);
        let res_texture = format!("{}{}", self.tex_folder_short, tex_file_name);
        // `return theTexture->WriteImage(aResTexFile);` (`cxx:231`).
        let ok = texture.write_image(&res_tex_file);
        (res_texture, ok)
    }

    /// `RWMesh_MaterialMap::copyFileTo` (`RWMesh_MaterialMap.cxx:120-148`): an
    /// empty side or a missing source fails, identical paths succeed without
    /// touching the file, otherwise the bytes are copied.
    fn copy_file_to(file_src: &str, file_dst: &str) -> bool {
        if file_src.is_empty() || file_dst.is_empty() {
            return false;
        } else if file_src == file_dst {
            return true;
        }
        // `OSD_File aFileSrc(aSrcPath); if (!aFileSrc.Exists()) ... return false;`
        // (`cxx:132-139`).
        if !std::path::Path::new(file_src).exists() {
            return false;
        }
        // `aFileSrc.Copy(aDstPath); return !aFileSrc.Failed();` (`cxx:140-141`).
        std::fs::copy(file_src, file_dst).is_ok()
    }

    /// `IsFailed()` (`hxx:90`).
    pub fn is_failed(&self) -> bool {
        self.is_failed
    }

    /// `myMatNameAsKey` (`RWMesh_MaterialMap.hxx:95`): the base class defaults
    /// it to TRUE; `RWGltf_GltfMaterialMap`'s constructor clears it
    /// (`RWGltf_GltfMaterialMap.cxx:58`).
    pub(crate) fn set_mat_name_as_key(&mut self, value: bool) {
        self.mat_name_as_key = value;
    }
}

// ---------------------------------------------------------------------------
// RWObj_ObjMaterialMap
// ---------------------------------------------------------------------------

/// `RWObj_ObjMaterialMap` (`RWObj_ObjMaterialMap.cxx`): the OBJ sidecar
/// material library writer. `AddMaterial` opens the file lazily on the first
/// material (`cxx:50-58`) and `DefineMaterial` (`cxx:72-141`) emits the body.
pub struct ObjMaterialMap {
    base: MeshMaterialMap,
    /// `myFileName` (`RWMesh_MaterialMap.hxx:85`).
    file_name: String,
    /// `myFile` (`RWObj_ObjMaterialMap.hxx`); `None` until the first material or
    /// after a failure.
    content: Option<String>,
}

impl ObjMaterialMap {
    /// `RWObj_ObjMaterialMap::RWObj_ObjMaterialMap` (`cxx:23-27`).
    pub fn new(file_name: &str) -> Self {
        Self {
            base: MeshMaterialMap::new(file_name),
            file_name: file_name.to_string(),
            content: None,
        }
    }

    /// `RWObj_ObjMaterialMap::AddMaterial` (`cxx:48-66`). Returns the material
    /// key (`""` when the file could not be created, `cxx:60-63`).
    pub fn add_material(&mut self, style: &XcafPrsStyle) -> String {
        if self.content.is_none() {
            // `myFile = OSD_OpenFile(myFileName, "wb")` (`cxx:53`): the
            // in-memory port always "opens" successfully.
            self.content = Some(String::from(
                "# Exported by Open CASCADE Technology [dev.opencascade.org]\n",
            ));
        }

        let (key, name, is_new) = self.base.add_material(style);
        if is_new {
            self.define_material(style, &key, &name);
        }
        key
    }

    /// `RWObj_ObjMaterialMap::DefineMaterial` (`RWObj_ObjMaterialMap.cxx:72-141`).
    fn define_material(&mut self, style: &XcafPrsStyle, key: &str, _name: &str) {
        if self.content.is_none() {
            return;
        }
        // The C++ writes record by record into `myFile`; the port accumulates the
        // same lines in `body` so the texture branch can borrow `self.base` (the
        // material map) while composing the output.
        let mut body = format!("newmtl {key}\n");

        let mut has_material = false;
        // `const XCAFDoc_VisMaterialCommon aDefMat = !myDefaultStyle.Material()
        // .IsNull() ? ... : XCAFDoc_VisMaterialCommon();` (`cxx:78-80`). The
        // writer's default style carries no material here (`Perform`'s
        // `theDefaultStyle` is empty on the DRAWEXE `WriteObj` path), so the
        // default common material values apply
        // (`XCAFDoc_VisMaterialCommon.hxx:48-56`).
        let default_material = XcafDocVisMaterialCommon::default();
        let mut amb_q = default_material.ambient_color;
        let mut diff_q = default_material.diffuse_color;
        let mut spec_q = default_material.specular_color;
        let mut transp = 0.0f32;
        // `aDefMat.Shininess * 1000.0f` (`cxx:84`); Shininess default 1.0.
        let mut specular = default_material.shininess * 1000.0;

        // `if (!theStyle.Material().IsNull() && !theStyle.Material()->IsEmpty())`
        // (`cxx:85-92`): the style's own material, converted to the common
        // model.
        if let Some(material) = style.material() {
            if !material.is_empty() {
                has_material = true;
                let com = material.convert_to_common_material();
                amb_q = com.ambient_color;
                diff_q = com.diffuse_color;
                spec_q = com.specular_color;
                transp = com.transparency;
                specular = com.shininess * 1000.0;
            }
        }
        if style.is_set_color_surf() {
            has_material = true;
            let c = style.get_color_surf();
            // `aDiffQ = theStyle.GetColorSurf()`; ambient is
            // `Vec3(color) * 0.25f` through the *linear* `Quantity_Color`
            // constructor (`cxx:96-99`).
            diff_q = LinearRgb::from_srgb(c);
            amb_q = diff_q.scaled(0.25);
            if style.get_color_surf_rgba().a < 1.0 {
                // `aTransp = 1.0f - theStyle.GetColorSurfRGBA().Alpha()`
                // (`cxx:100-103`).
                transp = 1.0 - style.get_color_surf_rgba().a;
            }
        }

        if has_material {
            // `Values(..., Quantity_TOC_sRGB)` (`cxx:108-110`) converts the
            // internally-linear colors back to sRGB for the file.
            let (ar, ag, ab) = amb_q.values_srgb();
            let (dr, dg, db) = diff_q.values_srgb();
            let (sr, sg, sb) = spec_q.values_srgb();
            body.push_str(&format!("Ka {ar:.6} {ag:.6} {ab:.6}\n"));
            body.push_str(&format!("Kd {dr:.6} {dg:.6} {db:.6}\n"));
            body.push_str(&format!("Ks {sr:.6} {sg:.6} {sb:.6}\n"));
            body.push_str(&format!("Ns {specular:.6}\n"));
            if transp >= 0.0001 {
                // `cxx:118-121`.
                body.push_str(&format!("Tr {transp:.6}\n"));
            }
        }

        // `if (const occ::handle<Image_Texture>& aBaseTexture =
        // theStyle.BaseColorTexture())` (`cxx:123-140`).
        if let Some(base_texture) = style.base_color_texture() {
            // `TCollection_AsciiString aTexture;` — assigned by the lookup or by
            // `CopyTexture` (`cxx:125-136`).
            let mut tex_out = String::new();
            let tex_id = base_texture.texture_id().to_string();
            if let Some((_, path)) = self.base.image_map.iter().find(|(t, _)| *t == tex_id) {
                tex_out = path.clone();
            } else if !self.base.image_fail_map.iter().any(|t| *t == tex_id) {
                // `CopyTexture(aTexture, aBaseTexture,
                // TCollection_AsciiString(myImageMap.Extent() + 1))`
                // (`cxx:128-129`): the key is the 1-based next index.
                let key_idx = format!("{}", self.base.image_map.len() + 1);
                let (res_texture, ok) = self.base.copy_texture(base_texture, &key_idx);
                if ok {
                    // `myImageMap.Bind(aBaseTexture, aTexture);` (`cxx:130-132`).
                    self.base.image_map.push((tex_id.clone(), res_texture.clone()));
                } else {
                    // `myImageFailMap.Add(aBaseTexture);` (`cxx:134-136`).
                    self.base.image_fail_map.push(tex_id.clone());
                }
                tex_out = res_texture;
            }
            // `if (!aTexture.IsEmpty()) { Fprintf(myFile, "map_Kd %s\n", ...); }`
            // (`cxx:137-140`). A failed copy still leaves `theResTexture` set
            // (`RWMesh_MaterialMap.cxx:212-213` / `cxx:227-228`), so the line is
            // written either way.
            if !tex_out.is_empty() {
                body.push_str(&format!("map_Kd {tex_out}\n"));
            }
        }

        if let Some(content) = self.content.as_mut() {
            content.push_str(&body);
        }
    }

    /// The accumulated `.mtl` text (`""` before any material).
    pub fn content(&self) -> &str {
        self.content.as_deref().unwrap_or("")
    }

    /// `RWObj_ObjMaterialMap::~RWObj_ObjMaterialMap` / `Close`
    /// (`cxx:30-38`). The in-memory port has nothing to flush.
    pub fn close(&mut self) {}

    /// `RWMesh_MaterialMap::FileName` equivalent (not exposed by OCCT; kept for
    /// the writer's diagnostics).
    pub fn file_name(&self) -> &str {
        &self.file_name
    }
}

