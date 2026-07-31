//! Phase 5 module: xcaf — STEP assembly metadata (`XCAFDoc_ShapeTool`-lite).
//!
//! A minimal stand-in for OCCT's XCAF document: tracks the product list of an
//! assembly and per-product name / color / layer attributes, serializes them
//! alongside a STEP physical file, and reads them back. Because a fully valid
//! XCAF (STEP 214 `APPLICATION_PROTOCOL`) layer is out of scope, the metadata
//! is embedded as a comment block before the `DATA` section using a clear
//! marker, which keeps the STEP file valid and round-trippable.

use std::collections::HashMap;

use crate::model::BRepModel;

/// Assembly metadata: product list plus per-product name/color/layer.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct StepAssembly {
    pub products: Vec<String>,
    pub names: HashMap<String, String>,
    pub colors: HashMap<String, (f64, f64, f64)>,
    pub layers: HashMap<String, String>,
}

impl StepAssembly {
    /// Register a product (idempotent).
    pub fn add_product(&mut self, name: &str) {
        if !self.products.iter().any(|p| p == name) {
            self.products.push(name.to_string());
        }
        self.names.entry(name.to_string()).or_insert_with(|| name.to_string());
    }

    /// Set the RGB color of a product.
    pub fn set_color(&mut self, name: &str, rgb: (f64, f64, f64)) {
        self.add_product(name);
        self.colors.insert(name.to_string(), rgb);
    }

    /// Set the layer of a product.
    pub fn set_layer(&mut self, name: &str, layer: &str) {
        self.add_product(name);
        self.layers.insert(name.to_string(), layer.to_string());
    }

    /// Look up a product by name (returns its display name).
    pub fn find(&self, name: &str) -> Option<&str> {
        self.names.get(name).map(|s| s.as_str())
    }
}

/// Extract the assembly metadata (name/color/layer) from each model shape.
pub fn model_to_step_assembly(model: &BRepModel) -> StepAssembly {
    let mut asm = StepAssembly::default();
    for ms in &model.shapes {
        asm.add_product(&ms.name);
        if let Some(c) = ms.color {
            asm.set_color(&ms.name, (c.r as f64, c.g as f64, c.b as f64));
        }
        if !ms.layer.is_empty() {
            asm.set_layer(&ms.name, &ms.layer);
        }
    }
    asm
}

/// Marker lines delimiting the embedded assembly block.
const BLOCK_BEGIN: &str = "/* ASSEMBLY-BEGIN";
const BLOCK_END: &str = "ASSEMBLY-END";

/// The STEP text for `model` with the assembly metadata embedded as a comment
/// block before the `DATA` section.
///
/// The STEP body itself is produced by `crate::step::write_step`, which
/// already writes the product names into `ADVANCED_BREP_SHAPE_REPRESENTATION`
/// records. The comment block adds the color/layer attributes that plain STEP
/// does not carry. Lines are formatted `name|r,g,b|layer`; names or layers
/// containing a `|` are sanitized (replaced with `/`).
pub fn write_step_with_metadata(model: &BRepModel) -> String {
    let step = crate::step::write_step(model);
    let mut block = String::new();
    block.push_str(BLOCK_BEGIN);
    block.push('\n');
    for ms in &model.shapes {
        let name = sanitize(&ms.name);
        let color = match ms.color {
            Some(c) => format2(c),
            None => String::new(),
        };
        let layer = sanitize(&ms.layer);
        block.push_str(&format!("{name}|{color}|{layer}\n"));
    }
    block.push_str(BLOCK_END);
    block.push_str("\n*/\n");
    step.replacen("DATA;", &format!("{block}DATA;"), 1)
}

/// Format a color as `r,g,b` (rounded to 3 decimals, stable round-trip).
fn format2(c: occt_core::quantity::Color) -> String {
    let r = (c.r * 1000.0).round() / 1000.0;
    let g = (c.g * 1000.0).round() / 1000.0;
    let b = (c.b * 1000.0).round() / 1000.0;
    format!("{r},{g},{b}")
}

/// Replace characters that would break the `|`-delimited comment format.
fn sanitize(s: &str) -> String {
    s.replace(['|', '\n', '\r'], "/")
}

/// Parse the embedded assembly comment block back into a [`StepAssembly`].
///
/// Lines are `name|r,g,b|layer`; empty color/layer fields are skipped.
pub fn extract_metadata_from_step(step_text: &str) -> StepAssembly {
    let mut asm = StepAssembly::default();
    let Some(begin) = step_text.find(BLOCK_BEGIN) else {
        return asm;
    };
    let rest = &step_text[begin + BLOCK_BEGIN.len()..];
    let body = &rest[..rest.find(BLOCK_END).unwrap_or(rest.len())];
    for line in body.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('*') {
            continue;
        }
        let parts: Vec<&str> = line.split('|').collect();
        let name = parts[0].trim();
        if name.is_empty() {
            continue;
        }
        asm.add_product(name);
        if parts.len() >= 2 {
            let rgb: Vec<f64> = parts[1]
                .split(',')
                .filter_map(|s| s.trim().parse().ok())
                .collect();
            if rgb.len() == 3 {
                asm.colors.insert(name.to_string(), (rgb[0], rgb[1], rgb[2]));
            }
        }
        if parts.len() >= 3 && !parts[2].trim().is_empty() {
            asm.layers.insert(name.to_string(), parts[2].trim().to_string());
        }
    }
    asm
}

/// Write `model` with metadata and read it back, verifying the shape names
/// survive both the embedded block and the STEP representation records.
pub fn model_with_metadata_step_roundtrip(model: &BRepModel) -> Result<(), String> {
    let step = write_step_with_metadata(model);
    let asm = extract_metadata_from_step(&step);
    let got = crate::step::read_step(&step)?;
    for name in model.names() {
        if asm.find(name).is_none() {
            return Err(format!("assembly metadata lost product '{name}'"));
        }
        if !got.names().iter().any(|n| *n == name) {
            return Err(format!("STEP read-back lost shape name '{name}'"));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::{BRepPrimBox, BRepPrimSphere};
    use occt_core::quantity::Color;

    fn two_shape_model() -> BRepModel {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let s = BRepPrimSphere::make_sphere(0.5);
        let mut model = BRepModel::new();
        model.add_with_color("RedBox", b.solid.0.clone(), Color::RED);
        model.add("GreyBall", s.solid.0.clone());
        model.find_mut("GreyBall").unwrap().layer = "L1".into();
        model
    }

    #[test]
    fn model_to_assembly_has_both_products() {
        let model = two_shape_model();
        let asm = model_to_step_assembly(&model);
        assert!(asm.find("RedBox").is_some());
        assert!(asm.find("GreyBall").is_some());
        assert_eq!(asm.colors.get("RedBox"), Some(&(1.0, 0.0, 0.0)));
        assert_eq!(asm.layers.get("GreyBall").map(|s| s.as_str()), Some("L1"));
        assert_eq!(asm.products.len(), 2);
    }

    #[test]
    fn write_contains_names() {
        let model = two_shape_model();
        let step = write_step_with_metadata(&model);
        assert!(step.contains("RedBox"), "step:\n{step}");
        assert!(step.contains("GreyBall"));
        assert!(step.contains(BLOCK_BEGIN));
    }

    #[test]
    fn extract_roundtrips_names_and_colors() {
        let model = two_shape_model();
        let step = write_step_with_metadata(&model);
        let asm = extract_metadata_from_step(&step);
        assert!(asm.find("RedBox").is_some());
        assert!(asm.find("GreyBall").is_some());
        assert_eq!(asm.colors.get("RedBox"), Some(&(1.0, 0.0, 0.0)));
        assert_eq!(asm.layers.get("GreyBall").map(|s| s.as_str()), Some("L1"));
        assert_eq!(asm.products.len(), 2);
    }

    #[test]
    fn full_roundtrip_succeeds() {
        let model = two_shape_model();
        model_with_metadata_step_roundtrip(&model).expect("roundtrip ok");
    }
}
