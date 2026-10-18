//! `STEPConstruct_RenderingProperties`. Source:
//! `STEPConstruct_RenderingProperties.cxx` / `.hxx`.
//!
//! Reader side only: `Init(SurfaceStyleRenderingWithProperties)`
//! (`STEPConstruct_RenderingProperties.cxx:356-442`) and the accessors
//! `IsDefined` / `GetRGBAColor` / `SurfaceColor` / `Transparency` the STEP color
//! reader and the OBJ material writer use.
//!
//! **UNPORTED**: the writer-side constructors and queries —
//! `CreateRenderingProperties` (`cxx:194-282`), `CreateXCAFMaterial`
//! (`cxx:286-352`), `IsMaterialConvertible` (`cxx:652-656`) and the
//! `XCAFDoc_VisMaterialCommon` / `XCAFDoc_VisMaterial` / `Quantity_ColorRGBA`
//! `Init` overloads (`cxx:446-631`) — because this crate has no
//! `XCAFDoc_VisMaterial` and its STEP writer emits the plain `COLOUR_RGB`
//! `STYLED_ITEM` form only (`step/transfer.rs::write_step_model_with_colors`).
//! The reflectance fields are still carried, because `Init` sets `myIsDefined`
//! through them: a `SURFACE_STYLE_RENDERING_WITH_PROPERTIES` whose only property
//! is a reflectance makes the style defined even without a color.

use std::collections::HashMap;

use occt_core::quantity::{Color, ColorRGBA};

use super::construct_styles::decode_color;
use super::*;

/// `StepVisual_ShadingSurfaceMethod` (`StepVisual_ShadingSurfaceMethod.hxx`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShadingSurfaceMethod {
    ConstantShading,
    ColourShading,
    DotShading,
    NormalShading,
}

/// `STEPConstruct_RenderingProperties` (`STEPConstruct_RenderingProperties.hxx:83-224`).
#[derive(Debug, Clone)]
pub struct RenderingProperties {
    /// `mySurfaceColor` — a `Quantity_Color`, i.e. sRGB in this port.
    surface_color: Color,
    /// `myTransparency` (0.0 opaque - 1.0 fully transparent).
    transparency: f64,
    /// `myRenderingMethod`.
    rendering_method: ShadingSurfaceMethod,
    /// `myIsDefined`.
    is_defined: bool,
    /// `myAmbientReflectance` — `(value, is_set)`.
    ambient_reflectance: (f64, bool),
    /// `myDiffuseReflectance` — `(value, is_set)`.
    diffuse_reflectance: (f64, bool),
    /// `mySpecularReflectance` — `(value, is_set)`.
    specular_reflectance: (f64, bool),
    /// `mySpecularExponent` — `(value, is_set)`.
    specular_exponent: (f64, bool),
    /// `mySpecularColour` — `(value, is_set)`.
    specular_color: (Color, bool),
}

impl Default for RenderingProperties {
    /// `STEPConstruct_RenderingProperties::STEPConstruct_RenderingProperties`
    /// (`cxx:33-44`): `mySurfaceColor` is `Quantity_NOC_WHITE`.
    fn default() -> Self {
        Self {
            surface_color: Color::WHITE,
            transparency: 0.0,
            rendering_method: ShadingSurfaceMethod::NormalShading,
            is_defined: false,
            ambient_reflectance: (0.0, false),
            diffuse_reflectance: (0.0, false),
            specular_reflectance: (0.0, false),
            specular_exponent: (0.0, false),
            specular_color: (Color::WHITE, false),
        }
    }
}

impl RenderingProperties {
    /// `STEPConstruct_RenderingProperties::IsDefined` (`hxx:172`).
    pub fn is_defined(&self) -> bool {
        self.is_defined
    }

    /// `STEPConstruct_RenderingProperties::GetRGBAColor` (`hxx:146-149`).
    pub fn rgba_color(&self) -> ColorRGBA {
        ColorRGBA {
            r: self.surface_color.r,
            g: self.surface_color.g,
            b: self.surface_color.b,
            a: (1.0 - self.transparency) as f32,
        }
    }

    /// `STEPConstruct_RenderingProperties::SurfaceColor` (`hxx:153`).
    pub fn surface_color(&self) -> Color {
        self.surface_color
    }

    /// `STEPConstruct_RenderingProperties::Transparency` (`hxx:157`).
    pub fn transparency(&self) -> f64 {
        self.transparency
    }

    /// `STEPConstruct_RenderingProperties::RenderingMethod`
    /// (`STEPConstruct_RenderingProperties.hxx`).
    pub fn rendering_method(&self) -> ShadingSurfaceMethod {
        self.rendering_method
    }

    /// `STEPConstruct_RenderingProperties::Init(SurfaceStyleRenderingWithProperties)`
    /// (`cxx:356-442`).
    ///
    /// `rendering_properties` is the id of `@N =
    /// SURFACE_STYLE_RENDERING_WITH_PROPERTIES(rendering_method, surface_color,
    /// (properties))`
    /// (`RWStepVisual_RWSurfaceStyleRenderingWithProperties.cxx:35-108`).
    pub fn init(&mut self, records: &HashMap<usize, Record>, rendering_properties: usize) {
        self.surface_color = Color::WHITE;
        self.transparency = 0.0;
        self.rendering_method = ShadingSurfaceMethod::NormalShading;
        self.is_defined = false;
        self.ambient_reflectance = (0.0, false);
        self.diffuse_reflectance = (0.0, false);
        self.specular_reflectance = (0.0, false);
        self.specular_exponent = (0.0, false);
        self.specular_color = (Color::WHITE, false);

        let Some(rec) = records.get(&rendering_properties) else {
            return;
        };
        if rec.type_name != "SURFACE_STYLE_RENDERING_WITH_PROPERTIES" {
            return;
        }
        // `myRenderingMethod = theRenderingProperties->RenderingMethod()`
        // (`cxx:376`).
        self.rendering_method = match rec.args.first().map(|s| s.trim()) {
            Some(".CONSTANT_SHADING.") => ShadingSurfaceMethod::ConstantShading,
            Some(".COLOUR_SHADING.") => ShadingSurfaceMethod::ColourShading,
            Some(".DOT_SHADING.") => ShadingSurfaceMethod::DotShading,
            _ => ShadingSurfaceMethod::NormalShading,
        };

        // `occ::handle<StepVisual_Colour> aColor =
        // theRenderingProperties->SurfaceColour(); if (!aColor.IsNull()) { ... }`
        // (`cxx:377-386`).
        if let Some(color) = rec.args.get(1).and_then(|a| parse_ref(a)) {
            if let Some(decoded) = decode_color(records, color) {
                self.surface_color = decoded;
                self.is_defined = true;
            }
        }

        // `aProperties = theRenderingProperties->Properties()` (`cxx:388-442`).
        let props = rec.args.get(2).map(String::as_str).unwrap_or("()");
        for prop in parse_ref_list(props) {
            let Some(prec) = records.get(&prop) else { continue };
            let args = &prec.args;
            match prec.type_name.as_str() {
                "SURFACE_STYLE_TRANSPARENT" => {
                    // `cxx:396-403`.
                    if let Ok(t) = parse_f64(args.first().map(String::as_str).unwrap_or("$")) {
                        self.transparency = t;
                        self.is_defined = true;
                    }
                }
                // The `RenderingPropertiesSelect` downcasts in `cxx:405-424` are
                // on the *handle*, so the most derived member of the
                // `SURFACE_STYLE_REFLECTANCE_AMBIENT(_DIFFUSE(_SPECULAR))`
                // family matches first, exactly as below.
                "SURFACE_STYLE_REFLECTANCE_AMBIENT_DIFFUSE_SPECULAR" => {
                    // `RWStepVisual_RWSurfaceStyleReflectanceAmbientDiffuseSpecular.cxx:46-64`.
                    if let (Ok(amb), Ok(dif), Ok(spec), Ok(exp)) = (
                        parse_f64(args.first().map(String::as_str).unwrap_or("$")),
                        parse_f64(args.get(1).map(String::as_str).unwrap_or("$")),
                        parse_f64(args.get(2).map(String::as_str).unwrap_or("$")),
                        parse_f64(args.get(3).map(String::as_str).unwrap_or("$")),
                    ) {
                        self.is_defined = true;
                        self.ambient_reflectance = (amb, true);
                        self.diffuse_reflectance = (dif, true);
                        self.specular_reflectance = (spec, true);
                        self.specular_exponent = (exp, true);
                        self.specular_color = match args.get(4).and_then(|a| parse_ref(a)) {
                            Some(c) => match decode_color(records, c) {
                                Some(decoded) => (decoded, true),
                                None => (Color::WHITE, false),
                            },
                            None => (Color::WHITE, false),
                        };
                    }
                }
                "SURFACE_STYLE_REFLECTANCE_AMBIENT_DIFFUSE" => {
                    // `cxx:426-434`; attribute layout
                    // `RWStepVisual_RWSurfaceStyleReflectanceAmbientDiffuse.cxx:46-49`.
                    if let (Ok(amb), Ok(dif)) = (
                        parse_f64(args.first().map(String::as_str).unwrap_or("$")),
                        parse_f64(args.get(1).map(String::as_str).unwrap_or("$")),
                    ) {
                        self.ambient_reflectance = (amb, true);
                        self.diffuse_reflectance = (dif, true);
                        self.is_defined = true;
                    }
                }
                "SURFACE_STYLE_REFLECTANCE_AMBIENT" => {
                    // `cxx:436-440`; layout
                    // `RWStepVisual_RWSurfaceStyleReflectanceAmbient.cxx:46-49`.
                    if let Ok(amb) = parse_f64(args.first().map(String::as_str).unwrap_or("$")) {
                        self.ambient_reflectance = (amb, true);
                        self.is_defined = true;
                    }
                }
                _ => {}
            }
        }
    }
}
