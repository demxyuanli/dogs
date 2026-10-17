//! `XCAFPrs_Style` and the `XCAFPrs` style helpers. Source: `XCAFPrs_Style.hxx`,
//! `XCAFPrs.cxx`.
//!
//! `XCAFPrs_Style` is the carrier between the STEP colour reader
//! (`STEPCAFControl_Reader::ReadColors`) and the mesh writers
//! (`RWMesh_ShapeIterator`, `RWMesh_FaceIterator`, `RWObj_ObjMaterialMap`). The
//! port carries the shape / surface colour / curve colour / visibility subset
//! those callers read back:
//!
//! * [`XcafPrsStyle`] — `XCAFPrs_Style` (`XCAFPrs_Style.hxx:27-120`).
//! * [`dispatch_face_styles`] — `RWMesh_ShapeIterator::dispatchStyles`
//!   (`RWMesh_ShapeIterator.cxx:70-134`) restricted to `myShapeType ==
//!   TopAbs_FACE`, i.e. the function that turns the label-keyed style map of
//!   `XCAFPrs::CollectStyleSettings` into the per-face style the OBJ writer
//!   reads.
//!
//! **UNPORTED**: `Material()` is always unset — `XCAFDoc_VisMaterial` (common /
//! PBR material attributes, base colour texture) has no counterpart in this
//! crate, and the STEP read path never sets it either
//! (`STEPCAFControl_Reader::ReadColors` only calls `SetColor` / `SetVisibility`;
//! `ReadMaterials` is not ported). Every `Material()`-guarded branch of the OBJ
//! writer is therefore unreachable, as documented in
//! `rwmesh/obj_caf_writer.rs`. Also UNPORTED: `XCAFPrs::CollectStyleSettings`'s
//! layer (`XCAFDoc_LayerTool`) and SHUO arms (`XCAFPrs.cxx:114-330`) — this
//! crate's `BRepModel` has no document tree; the reader-side analogue lives in
//! `step::read_colors`.

use std::collections::HashMap;

use occt_core::quantity::{Color, ColorRGBA};

use crate::abs::ShapeType;
use crate::shape::TopoShape;
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::faces_of;

/// `XCAFPrs_Style` (`XCAFPrs_Style.hxx:27-120`): a set of styling settings
/// applicable to a (sub)shape.
#[derive(Debug, Clone)]
pub struct XcafPrsStyle {
    /// `myHasColorSurf` / `myColorSurf`.
    color_surf: Option<ColorRGBA>,
    /// `myHasColorCurv` / `myColorCurv`.
    color_curv: Option<Color>,
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
        self.color_surf.is_none() && self.color_curv.is_none() && self.is_visible
    }

    /// `XCAFPrs_Style::IsSetColorSurf` (`hxx:50`).
    pub fn is_set_color_surf(&self) -> bool {
        self.color_surf.is_some()
    }

    /// `XCAFPrs_Style::SetColorSurf` (`hxx:62`) — the `Quantity_ColorRGBA`
    /// overload. `Color` carries sRGB in this port, so the stored triple is
    /// sRGB, the same value `Quantity_Color::Values(..., Quantity_TOC_sRGB)`
    /// returns for OCCT's internally-linear colour.
    pub fn set_color_surf(&mut self, the_color: ColorRGBA) {
        self.color_surf = Some(the_color);
    }

    /// `XCAFPrs_Style::UnSetColorSurf` (`hxx:68`).
    pub fn unset_color_surf(&mut self) {
        self.color_surf = None;
    }

    /// `XCAFPrs_Style::GetColorSurf` (`hxx:53`) — the sRGB `Quantity_Color`.
    pub fn get_color_surf(&self) -> Color {
        let c = self.color_surf.unwrap_or(ColorRGBA { r: 0.0, g: 0.0, b: 0.0, a: 1.0 });
        Color::new(c.r, c.g, c.b)
    }

    /// `XCAFPrs_Style::GetColorSurfRGBA` (`hxx:58`).
    pub fn get_color_surf_rgba(&self) -> ColorRGBA {
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
    pub fn get_color_curv(&self) -> Color {
        self.color_curv.unwrap_or(Color::WHITE)
    }

    /// `XCAFPrs_Style::SetVisibility` (`hxx:80`).
    pub fn set_visibility(&mut self, the_visibility: bool) {
        self.is_visible = the_visibility;
    }

    /// `XCAFPrs_Style::IsVisible` (`hxx:83`).
    pub fn is_visible(&self) -> bool {
        self.is_visible
    }

    /// `XCAFPrs_Style::IsEqual` (`XCAFPrs_Style.hxx:107-122`) / `operator==`
    /// (`hxx:125`). Note the visibility short-circuit: an invisible style equals
    /// every other invisible style regardless of its unset colour slots
    /// (`hxx:112-115`). `myMaterial` is always null here (module header), so
    /// that term drops out.
    pub fn is_equal(&self, other: &XcafPrsStyle) -> bool {
        if self.is_visible != other.is_visible {
            return false;
        } else if !self.is_visible {
            return true;
        }
        self.color_surf.is_some() == other.color_surf.is_some()
            && self.color_curv.is_some() == other.color_curv.is_some()
            && (self.color_surf.is_none() || self.color_surf == other.color_surf)
            && (self.color_curv.is_none() || self.color_curv == other.color_curv)
    }

    /// `RWMesh_ShapeIterator::initShape`'s `myHasColor` flag
    /// (`RWMesh_ShapeIterator.cxx:145-162`): set by a material or by either
    /// colour slot. `Material()` is always null here (module header), so a
    /// surface or curve colour is what flags a face as coloured.
    pub fn has_color(&self) -> bool {
        self.is_set_color_surf() || self.is_set_color_curv()
    }

    /// `RWMesh_ShapeIterator::initShape`'s `myColor`
    /// (`RWMesh_ShapeIterator.cxx:149-162`): the material's base colour, else
    /// the surface colour, else the curve colour.
    pub fn color(&self) -> ColorRGBA {
        if self.is_set_color_surf() {
            self.get_color_surf_rgba()
        } else if self.is_set_color_curv() {
            let c = self.get_color_curv();
            ColorRGBA { r: c.r, g: c.g, b: c.b, a: 1.0 }
        } else {
            ColorRGBA { r: 0.0, g: 0.0, b: 0.0, a: 1.0 }
        }
    }

    /// The colour slots `RWMesh_ShapeIterator::dispatchStyles` fills from the
    /// parent style when the shape's own style leaves them unset
    /// (`RWMesh_ShapeIterator.cxx:103-112`).
    pub fn merge_missing_from(&mut self, other: &XcafPrsStyle) {
        if self.color_curv.is_none() && other.is_set_color_curv() {
            self.color_curv = Some(other.get_color_curv());
        }
        if self.color_surf.is_none() && other.is_set_color_surf() {
            self.color_surf = Some(other.get_color_surf_rgba());
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
        assert_eq!(s.get_color_surf(), Color::new(0.25, 0.5, 0.75));
        assert_eq!(s.get_color_surf_rgba().a, 0.5);
        assert!(s.has_color());
    }

    #[test]
    fn node_color_fills_missing_slot() {
        let mut face = XcafPrsStyle::default();
        face.set_color_curv(Color::RED);
        face.merge_missing_from(&XcafPrsStyle::from_color(Some(Color::GREEN)));
        assert!(face.is_set_color_surf());
        assert_eq!(face.get_color_surf(), Color::GREEN);
        assert_eq!(face.get_color_curv(), Color::RED);
    }
}
