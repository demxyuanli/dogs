//! Wavefront OBJ + MTL write side. Source: `RWObj` (`RWObj_CafWriter`,
//! `RWObj_ObjWriterContext`, `RWObj_ObjMaterialMap`) and `RWMesh_MaterialMap`.
//!
//! Faithful port of the OCCT export pipeline used by `DRAWEXE`'s `WriteObj`
//! command:
//!
//! * [`MeshMaterialMap`] — `RWMesh_MaterialMap` (`RWMesh_MaterialMap.cxx`):
//!   material registration and `mat_<n>` key / name generation.
//! * [`ObjMaterialMap`] — `RWObj_ObjMaterialMap`
//!   (`RWObj_ObjMaterialMap.cxx`): opens the `.mtl` lazily on the first
//!   material and writes `newmtl` / `Ka` / `Kd` / `Ks` / `Ns` / `Tr`.
//! * [`ObjWriterContext`] — `RWObj_ObjWriterContext`
//!   (`RWObj_ObjWriterContext.cxx`): the header (with `mtllib`), `usemtl`,
//!   `g`, `v`, `vn`, `vt` and `f` records.
//! * [`perform`] — `RWObj_CafWriter::Perform` (`RWObj_CafWriter.cxx:82-334`):
//!   the `toCreateMatFile` decision (`cxx:242-257`), the `.mtl` file naming
//!   (`cxx:150-157`) and the per-node `writeShape` loop.
//!
//! Not ported (no counterpart in this crate, documented rather than faked):
//! the texture path (`RWObj_ObjMaterialMap::DefineMaterial` `map_Kd` +
//! `RWMesh_MaterialMap::CopyTexture` / `CreateTextureFolder`,
//! `RWObj_ObjMaterialMap.cxx:122-141`), the `XCAFDoc_VisMaterial` common /
//! PBR material properties (`cxx:78-92`) and the `WriteHeader` file-info map
//! (`RWObj_ObjWriterContext.cxx:112-135`).
//!
//! The `v` / `vn` / `vt` / `f` index shifts follow `WriteTriangle`
//! (`RWObj_ObjWriterContext.cxx:154-206`) and `FlushFace` (`cxx:303-312`)
//! verbatim: `myElemPosFirst` / `myElemNormFirst` / `myElemUVFirst` all advance
//! by the face's node count, so the three streams carry the same index.

use occt_core::gp::{GpPnt, GpPnt2d};
use occt_core::poly::triangulation::Triangle;
use occt_core::quantity::{Color, ColorRGBA};

use crate::xcaf_prs::XcafPrsStyle;

// ---------------------------------------------------------------------------
// Quantity_Color color-space conversions
// ---------------------------------------------------------------------------

/// `Quantity_Color::Convert_LinearRGB_To_sRGB` — the **double** overload
/// (`Quantity_Color.hxx:297-300`). `Values(..., Quantity_TOC_sRGB)`
/// (`Quantity_Color.cxx:499-504`) calls this one on `(double)myRgb`, so the
/// `.mtl` numbers are produced in double precision even though `myRgb` is
/// stored as `float`.
fn linear_to_srgb(v: f64) -> f64 {
    if v <= 0.0031308 {
        v * 12.92
    } else {
        v.powf(1.0 / 2.4) * 1.055 - 0.055
    }
}

/// `Quantity_Color::Convert_sRGB_To_LinearRGB` — the **double** overload
/// (`Quantity_Color.hxx:311-316`). `SetValues(..., Quantity_TOC_sRGB)`
/// (`Quantity_Color.cxx:345-351`) uses it on the incoming `double` and then
/// truncates to `float` for `myRgb`, which is what the STEP reader does
/// (`STEPConstruct_Styles.cxx:933-939`).
fn srgb_to_linear(v: f64) -> f32 {
    let linear = if v <= 0.04045 {
        v / 12.92
    } else {
        ((v + 0.055) / 1.055).powf(2.4)
    };
    linear as f32
}

/// A linearly-stored RGB triple, standing in for `Quantity_Color`'s internal
/// `myRgb` (`Quantity_Color.hxx:36` — "other OCCT interfaces ... would expect
/// them in linear space"). `Color` itself carries sRGB (STEP / OBJ use sRGB).
#[derive(Debug, Clone, Copy)]
struct LinearRgb {
    r: f32,
    g: f32,
    b: f32,
}

impl LinearRgb {
    /// `Quantity_Color(const NCollection_Vec3<float>&)` — "Define color from
    /// linear RGB values" (`Quantity_Color.hxx:63-66`, `Quantity_Color.cxx`).
    fn new(r: f32, g: f32, b: f32) -> Self {
        Self { r, g, b }
    }

    /// `Quantity_Color::Values(..., Quantity_TOC_sRGB)`
    /// (`Quantity_Color.cxx:499-504`): linear `myRgb` -> non-linear sRGB,
    /// computed in double precision.
    fn to_srgb(self) -> (f64, f64, f64) {
        (
            linear_to_srgb(self.r as f64),
            linear_to_srgb(self.g as f64),
            linear_to_srgb(self.b as f64),
        )
    }

    /// `(NCollection_Vec3<float>)theStyle.GetColorSurf() * 0.25f`
    /// (`RWObj_ObjMaterialMap.cxx:99`). The cast reads `myRgb` (linear).
    fn scaled(self, k: f32) -> Self {
        Self::new(self.r * k, self.g * k, self.b * k)
    }
}

// ---------------------------------------------------------------------------
// RWMesh_MaterialMap
// ---------------------------------------------------------------------------

/// `RWMesh_MaterialMap` (`RWMesh_MaterialMap.hxx:22-95`): registers materials
/// and hands out the `mat_<n>` keys used by `newmtl` / `usemtl`.
///
/// Only the base-class key / name generation is ported; the abstract
/// `DefineMaterial` and the texture machinery (`CreateTextureFolder` /
/// `CopyTexture`) are out of scope (see the module header).
pub struct MeshMaterialMap {
    /// `myKeyPrefix` — `"mat_"` (`RWMesh_MaterialMap.cxx:31`).
    key_prefix: String,
    /// `myNbMaterials` (`RWMesh_MaterialMap.hxx:93`).
    nb_materials: i32,
    /// `myMatNameAsKey` (`RWMesh_MaterialMap.hxx:95`), TRUE by default.
    mat_name_as_key: bool,
    /// `myStyles` — a `NCollection_DoubleMap<XCAFPrs_Style, AsciiString>` used
    /// as style -> key and key -> style in `AddMaterial` (`cxx:54-110`).
    styles: Vec<(XcafPrsStyle, String)>,
}

impl MeshMaterialMap {
    fn new() -> Self {
        Self {
            key_prefix: "mat_".to_string(),
            nb_materials: 0,
            mat_name_as_key: true,
            styles: Vec::new(),
        }
    }

    /// `RWMesh_MaterialMap::FindMaterial` (`RWMesh_MaterialMap.hxx:41-48`).
    /// The `NCollection_DoubleMap` lookup uses `XCAFPrs_Style::IsEqual`
    /// (`RWMesh_MaterialMap.cxx:103`), i.e. [`XcafPrsStyle::is_equal`].
    fn find_material(&self, style: &XcafPrsStyle) -> Option<String> {
        self.styles.iter().find(|(s, _)| s.is_equal(style)).map(|(_, k)| k.clone())
    }

    /// `RWMesh_MaterialMap::AddMaterial` (`RWMesh_MaterialMap.cxx:53-116`):
    /// the key / name generation, then the `DefineMaterial` call. The port
    /// turns the virtual callback into a return value the caller writes with.
    ///
    /// Returns `(key, name)`; `key` is the `newmtl` name, `name` the pretty
    /// name (unused by the OBJ writer, `RWObj_ObjMaterialMap.cxx:76`).
    fn add_material(&mut self, style: &XcafPrsStyle) -> (String, String, bool) {
        if let Some(key) = self.find_material(style) {
            return (key, String::new(), false);
        }

        let mut mat_key;
        let mat_name_suffix;
        let mut mat_name;
        let mut counter: i32 = 0;
        // `int* aCounterPtr = &myNbMaterials;` (`cxx:60`) — with `myMatNameAsKey`
        // and no material the counter points at the local `aCounter`.
        let mut counter_is_local = false;
        if self.mat_name_as_key {
            // `theStyle.Material()` is never non-null in this crate (STEP
            // carries only a surface color), so the material-name branch
            // (`cxx:62-78`) is not reachable here.
            self.nb_materials += 1;
            mat_name_suffix = self.key_prefix.clone();
            mat_name = format!("{}{}", mat_name_suffix, self.nb_materials);
            mat_key = mat_name.clone();
            counter_is_local = true;
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
            if counter_is_local {
                counter += 1;
                mat_name = format!("{}{}", mat_name_suffix, counter);
                mat_key = mat_name.clone();
            } else {
                counter += 1;
                mat_key = counter.to_string();
                mat_name = format!("{}{}", mat_name_suffix, mat_key);
            }
        }

        self.styles.push((style.clone(), mat_key.clone()));
        (mat_key, mat_name, true)
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
            base: MeshMaterialMap::new(),
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
        let Some(content) = self.content.as_mut() else { return };
        content.push_str(&format!("newmtl {key}\n"));

        let mut has_material = false;
        // `const XCAFDoc_VisMaterialCommon aDefMat = !myDefaultStyle.Material()
        // .IsNull() ? ... : XCAFDoc_VisMaterialCommon();` (`cxx:78-80`). The
        // port's default style never carries a material, so the default common
        // material values are used (`XCAFDoc_VisMaterialCommon.hxx:48-56`).
        let mut amb_q = LinearRgb::new(0.1, 0.1, 0.1);
        let mut diff_q = LinearRgb::new(0.8, 0.8, 0.8);
        let mut spec_q = LinearRgb::new(0.2, 0.2, 0.2);
        let mut transp = 0.0f32;
        // `aDefMat.Shininess * 1000.0f` (`cxx:84`); Shininess default 1.0.
        let specular = 1.0f32 * 1000.0;

        // `if (!theStyle.Material().IsNull() && !theStyle.Material()->IsEmpty())`
        // (`cxx:85-92`) — unreachable here: STEP styles carry only a color.
        if style.is_set_color_surf() {
            has_material = true;
            let c = style.get_color_surf();
            // `aDiffQ = theStyle.GetColorSurf()`; ambient is
            // `Vec3(color) * 0.25f` through the *linear* `Quantity_Color`
            // constructor (`cxx:96-99`).
            diff_q = LinearRgb::new(
                srgb_to_linear(c.r as f64),
                srgb_to_linear(c.g as f64),
                srgb_to_linear(c.b as f64),
            );
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
            let (ar, ag, ab) = amb_q.to_srgb();
            let (dr, dg, db) = diff_q.to_srgb();
            let (sr, sg, sb) = spec_q.to_srgb();
            content.push_str(&format!("Ka {ar:.6} {ag:.6} {ab:.6}\n"));
            content.push_str(&format!("Kd {dr:.6} {dg:.6} {db:.6}\n"));
            content.push_str(&format!("Ks {sr:.6} {sg:.6} {sb:.6}\n"));
            content.push_str(&format!("Ns {specular:.6}\n"));
            if transp >= 0.0001 {
                // `cxx:118-121`.
                content.push_str(&format!("Tr {transp:.6}\n"));
            }
        }

        // `if (const occ::handle<Image_Texture>& aBaseTexture = ...)`
        // (`cxx:123-140`) — texture path, unported (no `Image_Texture` here).
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

// ---------------------------------------------------------------------------
// RWObj_ObjWriterContext
// ---------------------------------------------------------------------------

/// `RWObj_ObjWriterContext` (`RWObj_ObjWriterContext.hxx:30-105`): the
/// low-level OBJ record writer. Writes into an in-memory buffer rather than a
/// `FILE*`, but keeps every OCCT control branch that this crate can reach.
pub struct ObjWriterContext {
    /// `NbFaces` (`RWObj_ObjWriterContext.hxx:92`).
    pub nb_faces: i32,
    /// `myFile` — the in-memory port always "opens" successfully.
    content: String,
    /// `myActiveMaterial` (`RWObj_ObjWriterContext.hxx:97`).
    active_material: String,
    /// `myElemPosFirst` / `myElemNormFirst` / `myElemUVFirst`
    /// (`RWObj_ObjWriterContext.hxx:99-101`) — `NCollection_Vec4<int>` with
    /// `xyz()` used as the index shift.
    elem_pos_first: [i32; 3],
    elem_norm_first: [i32; 3],
    elem_uv_first: [i32; 3],
    /// `myHasNormals` / `myHasTexCoords` (`RWObj_ObjWriterContext.hxx:103-104`).
    has_normals: bool,
    has_tex_coords: bool,
}

impl ObjWriterContext {
    /// `RWObj_ObjWriterContext::RWObj_ObjWriterContext`
    /// (`RWObj_ObjWriterContext.cxx:60-76`).
    pub fn new(_name: &str) -> Self {
        Self {
            nb_faces: 0,
            content: String::new(),
            active_material: String::new(),
            elem_pos_first: [1, 1, 1],
            elem_norm_first: [1, 1, 1],
            elem_uv_first: [1, 1, 1],
            has_normals: false,
            has_tex_coords: false,
        }
    }

    /// `RWObj_ObjWriterContext::Close` (`cxx:92-98`).
    pub fn close(&mut self) -> bool {
        true
    }

    /// `RWObj_ObjWriterContext::HasNormals` (`hxx:46`).
    pub fn has_normals(&self) -> bool {
        self.has_normals
    }

    /// `RWObj_ObjWriterContext::SetNormals` (`hxx:48`).
    pub fn set_normals(&mut self, has: bool) {
        self.has_normals = has;
    }

    /// `RWObj_ObjWriterContext::HasTexCoords` (`hxx:54`).
    pub fn has_tex_coords(&self) -> bool {
        self.has_tex_coords
    }

    /// `RWObj_ObjWriterContext::SetTexCoords` (`hxx:57`).
    pub fn set_tex_coords(&mut self, has: bool) {
        self.has_tex_coords = has;
    }

    /// `RWObj_ObjWriterContext::ActiveMaterial` (`hxx:66`).
    pub fn active_material(&self) -> &str {
        &self.active_material
    }

    /// `RWObj_ObjWriterContext::WriteHeader` (`cxx:100-141`). `file_info` is
    /// empty on the `RWObj_CafWriter` path (the DRAWEXE call passes none).
    pub fn write_header(
        &mut self,
        nb_nodes: i32,
        nb_elems: i32,
        mat_lib: &str,
        file_info: &[(String, String)],
    ) -> bool {
        self.content.push_str(
            "# Exported by Open CASCADE Technology [dev.opencascade.org]\n",
        );
        self.content.push_str(&format!("#  Vertices: {nb_nodes}\n"));
        self.content.push_str(&format!("#     Faces: {nb_elems}\n"));
        for (key, value) in file_info {
            // `splitLines` (`cxx:22-57`) + the key/value join (`cxx:112-136`).
            let key_lines: Vec<&str> = split_lines(key);
            let val_lines: Vec<&str> = split_lines(value);
            for (i, l) in key_lines.iter().enumerate() {
                self.content.push_str(if i > 0 { "\n# " } else { "# " });
                self.content.push_str(l);
            }
            self.content.push_str(if key_lines.is_empty() { "# " } else { ":" });
            for (i, l) in val_lines.iter().enumerate() {
                self.content.push_str(if i > 0 { "\n# " } else { " " });
                self.content.push_str(l);
            }
            self.content.push('\n');
        }
        if !mat_lib.is_empty() {
            self.content.push_str(&format!("mtllib {mat_lib}\n"));
        }
        true
    }

    /// `RWObj_ObjWriterContext::WriteActiveMaterial` (`cxx:145-150`).
    pub fn write_active_material(&mut self, material: &str) -> bool {
        self.active_material = material.to_string();
        if !material.is_empty() {
            self.content.push_str(&format!("usemtl {material}\n"));
        } else {
            self.content.push_str("usemtl\n");
        }
        true
    }

    /// `RWObj_ObjWriterContext::WriteTriangle` (`cxx:154-206`).
    pub fn write_triangle(&mut self, tri: [i32; 3]) -> bool {
        let pos = [tri[0] + self.elem_pos_first[0], tri[1] + self.elem_pos_first[1], tri[2] + self.elem_pos_first[2]];
        if self.has_normals {
            let norm = [
                tri[0] + self.elem_norm_first[0],
                tri[1] + self.elem_norm_first[1],
                tri[2] + self.elem_norm_first[2],
            ];
            if self.has_tex_coords {
                let uv = [
                    tri[0] + self.elem_uv_first[0],
                    tri[1] + self.elem_uv_first[1],
                    tri[2] + self.elem_uv_first[2],
                ];
                self.content.push_str(&format!(
                    "f {}/{}/{} {}/{}/{} {}/{}/{}\n",
                    pos[0], uv[0], norm[0], pos[1], uv[1], norm[1], pos[2], uv[2], norm[2]
                ));
            } else {
                self.content.push_str(&format!(
                    "f {}//{} {}//{} {}//{}\n",
                    pos[0], norm[0], pos[1], norm[1], pos[2], norm[2]
                ));
            }
            return true;
        }
        if self.has_tex_coords {
            let uv = [
                tri[0] + self.elem_uv_first[0],
                tri[1] + self.elem_uv_first[1],
                tri[2] + self.elem_uv_first[2],
            ];
            self.content.push_str(&format!(
                "f {}/{} {}/{} {}/{}\n",
                pos[0], uv[0], pos[1], uv[1], pos[2], uv[2]
            ));
            return true;
        }
        self.content.push_str(&format!("f {} {} {}\n", pos[0], pos[1], pos[2]));
        true
    }

    /// `RWObj_ObjWriterContext::WriteVertex` (`cxx:274-277`).
    pub fn write_vertex(&mut self, value: [f32; 3]) -> bool {
        self.content
            .push_str(&format!("v {:.*} {:.*} {:.*}\n", 6, value[0], 6, value[1], 6, value[2]));
        true
    }

    /// `RWObj_ObjWriterContext::WriteNormal` (`cxx:280-283`).
    pub fn write_normal(&mut self, value: [f32; 3]) -> bool {
        self.content
            .push_str(&format!("vn {:.*} {:.*} {:.*}\n", 6, value[0], 6, value[1], 6, value[2]));
        true
    }

    /// `RWObj_ObjWriterContext::WriteTexCoord` (`cxx:287-290`).
    pub fn write_tex_coord(&mut self, value: [f32; 2]) -> bool {
        self.content.push_str(&format!("vt {:.*} {:.*}\n", 6, value[0], 6, value[1]));
        true
    }

    /// `RWObj_ObjWriterContext::WriteGroup` (`cxx:296-298`).
    pub fn write_group(&mut self, value: &str) -> bool {
        if !value.is_empty() {
            self.content.push_str(&format!("g {value}\n"));
        } else {
            self.content.push_str("g\n");
        }
        true
    }

    /// `RWObj_ObjWriterContext::FlushFace` (`cxx:303-312`).
    pub fn flush_face(&mut self, nb_nodes: i32) {
        for i in 0..3 {
            self.elem_pos_first[i] += nb_nodes;
        }
        if self.has_normals {
            for i in 0..3 {
                self.elem_norm_first[i] += nb_nodes;
            }
        }
        if self.has_tex_coords {
            for i in 0..3 {
                self.elem_uv_first[i] += nb_nodes;
            }
        }
    }

    /// The accumulated OBJ text.
    pub fn content(&self) -> &str {
        &self.content
    }
}

/// `splitLines` (`RWObj_ObjWriterContext.cxx:22-57`): split on `\r` / `\n`,
/// dropping empty pieces, right-adjusting each line.
fn split_lines(text: &str) -> Vec<&str> {
    if text.is_empty() {
        return Vec::new();
    }
    text.split(['\r', '\n']).filter(|l| !l.is_empty()).collect()
}

// ---------------------------------------------------------------------------
// RWObj_CafWriter::Perform
// ---------------------------------------------------------------------------

/// One `RWMesh_FaceIterator` face: the face's own node block plus the style the
/// shape iterator bound to it (`RWMesh_ShapeIterator::initShape`, `cxx:139-162`).
pub struct ObjFaceNode {
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
    // .IsNull() && HasTexCoords())` (`cxx:254-256`); the texture term is
    // unreachable here, and `HasFaceColor()` is [`XcafPrsStyle::has_color`].
    *to_create_mat_file = *to_create_mat_file || face.style.has_color();
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
        if face.style.has_color() {
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
pub fn perform(obj_path: &str, nodes: &[ObjExportNode]) -> (String, String) {
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
