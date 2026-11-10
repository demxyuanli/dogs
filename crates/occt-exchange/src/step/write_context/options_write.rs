use super::super::prelude::*;
use super::super::*;

use super::context::WriteCtx;
use super::options::emit_si_units;


/// Consolidated options for [`write_step_with_options`].
///
/// Combines the PRODUCT name, an optional color style, SI units and the
/// spline-capable writer into one call — the port of `STEPControl_Writer`'s
/// per-shape configuration (name attribute, presentation color and unit
/// system).
#[derive(Debug, Clone)]
pub struct StepWriteOptions {
    pub name: String,
    pub color: Option<(f64, f64, f64)>,
    pub units: StepUnits,
    pub splines: bool,
}

impl Default for StepWriteOptions {
    fn default() -> Self {
        Self {
            name: "Shape".into(),
            color: None,
            units: StepUnits::default(),
            splines: false,
        }
    }
}

/// Serialize a single shape with the full [`StepWriteOptions`] configuration.
///
/// This is the one-stop entry point: it writes the shape representation and
/// product scaffolding, then decorates it with a color style (when
/// `color` is set), an SI unit block, and (when `splines` is set) B-spline
/// output for non-analytic geometry.
pub fn write_step_with_options(shape: &TopoShape, opts: &StepWriteOptions) -> Result<String, String> {
    let mut ctx = WriteCtx::new();
    ctx.splines = opts.splines;
    let rep = ctx
        .write_shape_named(&opts.name, shape)
        .ok_or_else(|| "write_step_with_options: nothing to write".to_string())?;
    if let Some((r, g, b)) = opts.color {
        let color = ctx.w.emit(format!(
            "COLOUR_RGB('',{},{},{})",
            step_real(r),
            step_real(g),
            step_real(b)
        ));
        let fill = ctx.w.emit(format!("SURFACE_STYLE_FILL_AREA('',#{color})"));
        let usage = ctx.w.emit(format!("SURFACE_STYLE_USAGE('',#{fill})"));
        ctx.w.emit(format!("STYLED_ITEM('',(#{usage}),#{rep})"));
    }
    emit_si_units(&mut ctx.w, &opts.units);
    Ok(ctx.finish())
}

/// Serialize a single shape with an explicit name and a `COLOUR_RGB` style.
///
/// Combines [`write_step_with_name`] and [`write_step_with_color`]: the
/// product carries `name` and the shape representation is decorated with a
/// `SURFACE_STYLE_FILL_AREA` color.
pub fn write_step_with_name_and_color(
    shape: &TopoShape,
    name: &str,
    color: (f64, f64, f64),
) -> Result<String, String> {
    let mut ctx = WriteCtx::new();
    let rep = ctx
        .write_shape_named(name, shape)
        .ok_or_else(|| "write_step_with_name_and_color: nothing to write".to_string())?;
    let color = ctx.w.emit(format!(
        "COLOUR_RGB('',{},{},{})",
        step_real(color.0),
        step_real(color.1),
        step_real(color.2)
    ));
    let fill = ctx.w.emit(format!("SURFACE_STYLE_FILL_AREA('',#{color})"));
    let usage = ctx.w.emit(format!("SURFACE_STYLE_USAGE('',#{fill})"));
    ctx.w.emit(format!("STYLED_ITEM('',(#{usage}),#{rep})"));
    Ok(ctx.finish())
}
