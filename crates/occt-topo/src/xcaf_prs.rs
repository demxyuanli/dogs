//! `XCAFPrs_Style` and the `XCAFPrs` style helpers. Source: `XCAFPrs_Style.hxx`,
//! `XCAFPrs.cxx`.
//!
//! `XCAFPrs_Style` is the carrier between the colour / material readers
//! (`STEPCAFControl_Reader::ReadColors`, the document-level format readers) and
//! the mesh writers (`RWMesh_ShapeIterator`, `RWMesh_FaceIterator`,
//! `RWObj_ObjMaterialMap`). The port carries the surface colour / curve colour /
//! material / visibility subset those callers read back:
//!
//! * [`XcafPrsStyle`] — `XCAFPrs_Style` (`XCAFPrs_Style.hxx:27-120`).
//! * [`dispatch_face_styles`] — `RWMesh_ShapeIterator::dispatchStyles`
//!   (`RWMesh_ShapeIterator.cxx:70-134`) restricted to `myShapeType ==
//!   TopAbs_FACE`, i.e. the function that turns the label-keyed style map of
//!   `XCAFPrs::CollectStyleSettings` into the per-face style the OBJ writer
//!   reads.
//! * [`collect_style_settings`] — `XCAFPrs::CollectStyleSettings`
//!   (`XCAFPrs.cxx:116-330`) over a [`BRepModel`] label tree.
//!
//! `myMaterial` is the [`XcafDocVisMaterial`] the document attaches to the
//! label (`XCAFPrs_Style.hxx:40-48`); `BaseColorTexture()` (`hxx:86-103`) reads
//! its PBR base-colour or common diffuse texture. In OCCT the material is
//! written into the document by a document-level reader (`RWGltf_CafReader`,
//! `RWObj_CafReader`, XBF) and read back by `XCAFPrs::CollectStyleSettings`
//! through `aMatTool->GetShapeMaterial(aLabel)` (`XCAFPrs.cxx:174`, `cxx:251`).
//! `STEPCAFControl_Reader` never writes one — `rg VisMaterial` over it is empty,
//! while `STEPCAFControl_Writer.cxx:209,1542,2242` only *reads* it back for
//! export — so on the STEP -> OBJ path `Material()` stays null in OCCT too and
//! `RWObj_ObjMaterialMap`'s `cxx:85-92` arm is equally dead there.
//!
//! The XCAF-document glTF read that is the production producer of materials is
//! ported (`rwmesh::gltf_document`): the reader builds the bare-triangulation
//! proxy faces of `BRep_Builder::MakeFace(aFace, aMeshData)`
//! (`RWGltf_GltfJsonParser.cxx:2002-2006`) through the `triangulations` slot of
//! `FaceGeomCore`, fills them in `readLateData`
//! (`RWGltf_CafReader.cxx:410-430`) and replays the attributes into the
//! `BRepModel` labels.
//!
//! Also UNPORTED: [`collect_style_settings`]'s layer (`XCAFDoc_LayerTool`) and
//! SHUO arms (`XCAFPrs.cxx:114-330`) — this crate's `BRepModel` has no layer
//! tree and no NAUO/SHUO graph; the reader-side flat analogue lives in
//! `step::read_colors`.

use std::collections::HashMap;

use occt_core::quantity::{Color, ColorRGBA};

use crate::abs::ShapeType;
use crate::image_texture::ImageTexture;
use crate::model::BRepModel;
use crate::shape::TopoShape;
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::faces_of;
use crate::xcaf_doc_vis_material::XcafDocVisMaterial;

/// `XCAFPrs_Style` (`XCAFPrs_Style.hxx:27-120`): a set of styling settings
/// applicable to a (sub)shape.
#[derive(Debug, Clone)]
pub struct XcafPrsStyle {
    /// `myHasColorSurf` / `myColorSurf`.
    color_surf: Option<ColorRGBA>,
    /// `myHasColorCurv` / `myColorCurv`.
    color_curv: Option<Color>,
    /// `myMaterial` (`XCAFPrs_Style.hxx:134`): the visualization material the
    /// document attaches to the label. `None` stands for a null handle.
    material: Option<XcafDocVisMaterial>,
    /// `myIsVisible`.
    is_visible: bool,
}

impl Default for XcafPrsStyle {
    fn default() -> Self {
        // "Empty constructor - colors are unset, visibility is TRUE"
        // (`XCAFPrs_Style.hxx:30-33`).
        Self {
            color_surf: None,
            color_curv: None,
            material: None,
            is_visible: true,
        }
    }
}

impl XcafPrsStyle {
    /// Build a style from a model shape's optional surface color. Used by the
    /// OBJ export for shapes whose colour comes from `ModelShape.color` rather
    /// than from a `XCAFDoc_ColorTool` label.
    pub fn from_color(color: Option<Color>) -> Self {
        Self {
            color_surf: color.map(|c| ColorRGBA { r: c.r, g: c.g, b: c.b, a: 1.0 }),
            ..Self::default()
        }
    }

    /// `XCAFPrs_Style::IsEmpty` (`hxx:35-38`): does not override any property.
    pub fn is_empty(&self) -> bool {
        self.color_surf.is_none()
            && self.color_curv.is_none()
            && self.material.is_none()
            && self.is_visible
    }

    /// `XCAFPrs_Style::Material()` (`hxx:40`).
    pub fn material(&self) -> Option<&XcafDocVisMaterial> {
        self.material.as_ref()
    }

    /// `XCAFPrs_Style::SetMaterial()` (`hxx:43-46`).
    pub fn set_material(&mut self, material: Option<XcafDocVisMaterial>) {
        self.material = material;
    }

    /// `XCAFPrs_Style::BaseColorTexture` (`hxx:86-103`): the PBR base-colour
    /// texture when the material defines PBR and carries one, else the common
    /// diffuse texture, else null.
    pub fn base_color_texture(&self) -> Option<&ImageTexture> {
        let material = self.material.as_ref()?;
        if material.has_pbr_material() {
            if let Some(tex) = material.pbr_material().base_color_texture.as_ref() {
                return Some(tex);
            }
        }
        if material.has_common_material() {
            if let Some(tex) = material.common_material().diffuse_texture.as_ref() {
                return Some(tex);
            }
        }
        None
    }

    /// `XCAFPrs_Style::IsSetColorSurf` (`hxx:50`).
    pub fn is_set_color_surf(&self) -> bool {
        self.color_surf.is_some()
    }

    /// `XCAFPrs_Style::SetColorSurf` (`hxx:62`) — the `Quantity_ColorRGBA`
    /// overload. `Color` carries sRGB in this port, so the stored triple is
    /// sRGB, the same value `Quantity_Color::Values(..., Quantity_TOC_sRGB)`
    /// returns for OCCT's internally-linear colour.
    pub fn set_color_surf(&mut self, color: ColorRGBA) {
        self.color_surf = Some(color);
    }

    /// `XCAFPrs_Style::UnSetColorSurf` (`hxx:68`).
    pub fn unset_color_surf(&mut self) {
        self.color_surf = None;
    }

    /// `XCAFPrs_Style::GetColorSurf` (`hxx:53`) — the sRGB `Quantity_Color`.
    pub fn color_surf(&self) -> Color {
        let c = self.color_surf.unwrap_or(ColorRGBA { r: 0.0, g: 0.0, b: 0.0, a: 1.0 });
        Color::new(c.r, c.g, c.b)
    }

    /// `XCAFPrs_Style::GetColorSurfRGBA` (`hxx:58`).
    pub fn color_surf_rgba(&self) -> ColorRGBA {
        self.color_surf.unwrap_or(ColorRGBA { r: 0.0, g: 0.0, b: 0.0, a: 1.0 })
    }

    /// `XCAFPrs_Style::IsSetColorCurv` (`hxx:66`).
    pub fn is_set_color_curv(&self) -> bool {
        self.color_curv.is_some()
    }

    /// `XCAFPrs_Style::SetColorCurv` (`hxx:76`).
    pub fn set_color_curv(&mut self, col: Color) {
        self.color_curv = Some(col);
    }

    /// `XCAFPrs_Style::GetColorCurv` (`hxx:72`).
    pub fn color_curv(&self) -> Color {
        self.color_curv.unwrap_or(Color::WHITE)
    }

    /// `XCAFPrs_Style::SetVisibility` (`hxx:80`).
    pub fn set_visibility(&mut self, visibility: bool) {
        self.is_visible = visibility;
    }

    /// `XCAFPrs_Style::IsVisible` (`hxx:83`).
    pub fn is_visible(&self) -> bool {
        self.is_visible
    }

    /// `XCAFPrs_Style::IsEqual` (`XCAFPrs_Style.hxx:107-122`) / `operator==`
    /// (`hxx:125`). Note the visibility short-circuit: an invisible style equals
    /// every other invisible style regardless of its unset colour slots
    /// (`hxx:112-115`). `myMaterial == theOther.myMaterial` (`hxx:119`) compares
    /// handles by pointer, which [`XcafDocVisMaterial::id`] reproduces: two
    /// separately created materials stay distinct even when their values match.
    pub fn is_equal(&self, other: &XcafPrsStyle) -> bool {
        if self.is_visible != other.is_visible {
            return false;
        } else if !self.is_visible {
            return true;
        }
        let mat_equal = match (&self.material, &other.material) {
            (None, None) => true,
            (Some(a), Some(b)) => a.id() == b.id(),
            _ => false,
        };
        self.color_surf.is_some() == other.color_surf.is_some()
            && self.color_curv.is_some() == other.color_curv.is_some()
            && mat_equal
            && (self.color_surf.is_none() || self.color_surf == other.color_surf)
            && (self.color_curv.is_none() || self.color_curv == other.color_curv)
    }

    /// `RWMesh_ShapeIterator::initShape`'s `myHasColor` flag
    /// (`RWMesh_ShapeIterator.cxx:140-163`): set by any material handle or by
    /// either colour slot.
    pub fn has_color(&self) -> bool {
        self.material.is_some() || self.is_set_color_surf() || self.is_set_color_curv()
    }

    /// `RWMesh_ShapeIterator::initShape`'s `myColor`
    /// (`RWMesh_ShapeIterator.cxx:149-162`): the material's base colour, else
    /// the surface colour, else the curve colour. The material arm returns
    /// `Quantity_ColorRGBA`'s externally-visible sRGB triple.
    pub fn color(&self) -> ColorRGBA {
        if let Some(material) = self.material.as_ref() {
            material.base_color().to_srgb()
        } else if self.is_set_color_surf() {
            self.color_surf_rgba()
        } else if self.is_set_color_curv() {
            let c = self.color_curv();
            ColorRGBA { r: c.r, g: c.g, b: c.b, a: 1.0 }
        } else {
            ColorRGBA { r: 0.0, g: 0.0, b: 0.0, a: 1.0 }
        }
    }

    /// The colour slots `RWMesh_ShapeIterator::dispatchStyles` fills from the
    /// parent style when the shape's own style leaves them unset
    /// (`RWMesh_ShapeIterator.cxx:103-112`), plus the material slot
    /// (`cxx:113-116`: `if (aCafStyle.Material().IsNull() && !theStyle.Material()
    /// .IsNull()) aCafStyle.SetMaterial(theStyle.Material())`).
    pub fn merge_missing_from(&mut self, other: &XcafPrsStyle) {
        if self.color_curv.is_none() && other.is_set_color_curv() {
            self.color_curv = Some(other.color_curv());
        }
        if self.color_surf.is_none() && other.is_set_color_surf() {
            self.color_surf = Some(other.color_surf_rgba());
        }
        if self.material.is_none() && other.material.is_some() {
            self.material = other.material.clone();
        }
    }
}

/// `RWMesh_ShapeIterator::dispatchStyles` (`RWMesh_ShapeIterator.cxx:70-134`)
/// for `myShapeType == TopAbs_FACE`.
///
/// The input is the shape-keyed style map `XCAFPrs::CollectStyleSettings`
/// produces for one document node (`theSettings`) plus that node's own style
/// (`theStyle`); the output is the face-keyed map
/// `RWMesh_ShapeIterator::initShape` looks up
/// (`RWMesh_ShapeIterator.cxx:145-147`), i.e. `face shape key -> style`.
///
/// Faithful to the `aTypeIter` loop order: FACE keys bind directly (last write
/// wins, `cxx:119`), then SHELL, SOLID, COMPSOLID and COMPOUND keys expand to
/// their faces with `if (!myStyles.IsBound(...))` first-write-wins
/// (`cxx:126-129`) — so the most specific key wins, exactly as in OCCT. Every
/// bound style takes the node style's unset colour slots
/// (`cxx:103-112`); a face that no key covers is bound to the node style by the
/// caller (`RWMesh_ShapeIterator::initShape`'s `myDefStyle` fallback,
/// `cxx:145-148`).
pub fn dispatch_face_styles(
    settings: &[(TopoShape, XcafPrsStyle)],
    node_style: &XcafPrsStyle,
) -> HashMap<usize, XcafPrsStyle> {
    let mut out: HashMap<usize, XcafPrsStyle> = HashMap::new();
    // `for (int aTypeIter = myShapeType; aTypeIter >= TopAbs_COMPOUND;
    // --aTypeIter)` (`cxx:79`).
    for ty in [
        ShapeType::Face,
        ShapeType::Shell,
        ShapeType::Solid,
        ShapeType::CompSolid,
        ShapeType::Compound,
    ] {
        for (key_shape, style) in settings {
            if key_shape.shape_type() != ty {
                continue;
            }
            let mut style = style.clone();
            style.merge_missing_from(node_style);
            if ty == ShapeType::Face {
                // `myStyles.Bind(aKeyShapeLocated, aCafStyle)` (`cxx:119`).
                out.insert(GeometryRegistry::shape_key(key_shape), style);
            } else {
                // `else { for (TopExp_Explorer aShapeIter(aKeyShapeLocated,
                // myShapeType); ...) if (!myStyles.IsBound(...)) myStyles.Bind(...) }`
                // (`cxx:122-131`).
                for f in faces_of(key_shape) {
                    out.entry(GeometryRegistry::shape_key(f.shape()))
                        .or_insert_with(|| style.clone());
                }
            }
        }
    }
    out
}

/// `XCAFPrs::CollectStyleSettings` (`XCAFPrs.cxx:116-330`) for one label of a
/// [`BRepModel`] document: the shape-keyed style map
/// `RWMesh_ShapeIterator::dispatchStyles` reads (`theSettings`).
///
/// The label walk is faithful (`cxx:120-330`): the referred shape first, then
/// the components of an assembly, then the sub-shape labels and the label
/// itself, each binding `XCAFPrs_Style` under the label's shape with last write
/// winning.
///
/// UNPORTED, with the `cxx` lines they come from: the `XCAFDoc_LayerTool` arms
/// (`cxx:125-142`, `:178-216`: hidden layers, colour-by-layer and the
/// `theLayerColor` fallback — `BRepModel` has no layer tree), the SHUO arm
/// (`cxx:222-303`, no NAUO/SHUO graph) and `XCAFDoc_ColorTool::IsColorByLayer`.
/// `IsVisible(aLabel)` is always true here: labels carry no visibility
/// attribute, so the style keeps the constructor's visible default
/// (`XCAFPrs_Style.hxx:30-33`).
pub fn collect_style_settings(model: &BRepModel, label: usize) -> Vec<(TopoShape, XcafPrsStyle)> {
    let mut out = Vec::new();
    collect_style_settings_into(model, label, &mut out);
    out
}

fn collect_style_settings_into(
    model: &BRepModel,
    label: usize,
    out: &mut Vec<(TopoShape, XcafPrsStyle)>,
) {
    if label >= model.labels.len() {
        return;
    }
    // "for references, first collect colors of referred shape" (`cxx:120-145`).
    if let Some(ref_label) = model.referred_label(label) {
        collect_style_settings_into(model, ref_label, out);
    }
    // "for assemblies, first collect colors defined in components"
    // (`cxx:147-158`).
    if model.is_assembly(label) {
        for &component in &model.labels[label].children {
            collect_style_settings_into(model, component, out);
        }
    }
    // "collect settings on subshapes" plus the label itself (`cxx:160-169`).
    let mut candidates: Vec<usize> = model.labels[label]
        .children
        .iter()
        .copied()
        .filter(|&child| matches!(model.labels.get(child).map(|l| l.role), Some(crate::model::LabelRole::SubShape)))
        .collect();
    candidates.push(label);
    for candidate in candidates {
        let mut style = XcafPrsStyle::default();
        // `aStyle.SetMaterial(aMatTool->GetShapeMaterial(aLabel))` (`cxx:176`).
        if let Some(material) = model.labels[candidate].material.and_then(|i| model.material(i)) {
            style.set_material(Some(material.clone()));
        }
        // `fillStyleColors` (`XCAFPrs.cxx:34-53`).
        if let Some(c) = model.labels[candidate].color {
            style.set_color_surf(ColorRGBA { r: c.r, g: c.g, b: c.b, a: 1.0 });
        }
        if let Some(c) = model.labels[candidate].color_curv {
            style.set_color_curv(c);
        }
        // `if (aStyle.IsEmpty()) continue;` (`cxx:305-309`).
        if style.is_empty() {
            continue;
        }
        let Some(shape) = model.labels[candidate].shape.clone() else { continue };
        // `if (aSubshape.ShapeType() == TopAbs_COMPOUND && aSubshape.NbChildren()
        // == 0) continue;` (`cxx:311-318`).
        if shape.shape_type() == ShapeType::Compound && shape.tshape.read().unwrap().nb_children() == 0 {
            continue;
        }
        bind_style(out, shape, style);
    }
}

/// `theSettings.ChangeSeek` / `Add` (`XCAFPrs.cxx:321-329`): replace the style
/// bound to the same shape (`TopTools_ShapeMapHasher` is `IsSame`).
fn bind_style(out: &mut Vec<(TopoShape, XcafPrsStyle)>, shape: TopoShape, style: XcafPrsStyle) {
    let key = GeometryRegistry::shape_key(&shape);
    if let Some(slot) = out.iter_mut().find(|(s, _)| GeometryRegistry::shape_key(s) == key) {
        slot.1 = style;
    } else {
        out.push((shape, style));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_style_is_visible() {
        let s = XcafPrsStyle::default();
        assert!(s.is_empty() && s.is_visible() && !s.is_set_color_surf());
    }

    #[test]
    fn surf_color_roundtrip() {
        let mut s = XcafPrsStyle::default();
        s.set_color_surf(ColorRGBA { r: 0.25, g: 0.5, b: 0.75, a: 0.5 });
        assert!(s.is_set_color_surf() && !s.is_empty());
        assert_eq!(s.color_surf(), Color::new(0.25, 0.5, 0.75));
        assert_eq!(s.color_surf_rgba().a, 0.5);
        assert!(s.has_color());
    }

    #[test]
    fn node_color_fills_missing_slot() {
        let mut face = XcafPrsStyle::default();
        face.set_color_curv(Color::RED);
        face.merge_missing_from(&XcafPrsStyle::from_color(Some(Color::GREEN)));
        assert!(face.is_set_color_surf());
        assert_eq!(face.color_surf(), Color::GREEN);
        assert_eq!(face.color_curv(), Color::RED);
    }
}
