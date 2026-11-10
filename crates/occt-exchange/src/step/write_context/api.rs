use super::super::prelude::*;
use super::super::*;

use super::context::WriteCtx;


/// Serialize the whole model to a STEP physical file.
pub fn write_step(model: &BRepModel) -> String {
    let mut ctx = WriteCtx::new();
    for ms in &model.shapes {
        ctx.write_shape_named(&ms.name, &ms.shape);
    }
    ctx.finish()
}

/// Serialize the model to a STEP file on disk.
pub fn write_step_file(path: &str, model: &BRepModel) -> std::io::Result<()> {
    std::fs::write(path, write_step(model))
}

/// Serialize a single shape (wrapped in a one-shape model).
pub fn write_shape_step(shape: &TopoShape) -> String {
    let mut model = BRepModel::new();
    model.add("Shape", shape.clone());
    write_step(&model)
}

/// Serialize a single shape, emitting B-spline geometry for non-analytic
/// curves and surfaces (full `STEPControl_Writer` spline coverage).
///
/// Analytic geometry — planes/cylinders/cones/spheres/tori and
/// lines/circles/ellipses/parabolas — is written exactly. Anything else is
/// approximated by an interpolating B-spline and written as
/// `B_SPLINE_CURVE_WITH_KNOTS` / `B_SPLINE_SURFACE_WITH_KNOTS`, which the
/// reader reconstructs as a valid B-spline face/edge.
pub fn write_step_with_splines(shape: &TopoShape) -> Result<String, String> {
    let mut ctx = WriteCtx::new();
    ctx.splines = true;
    ctx.write_shape_named("Shape", shape)
        .ok_or_else(|| "write_step_with_splines: nothing to write".to_string())?;
    Ok(ctx.finish())
}

/// Serialize a single shape, overriding the `PRODUCT` name attribute.
pub fn write_step_with_name(shape: &TopoShape, name: &str) -> Result<String, String> {
    let mut ctx = WriteCtx::new();
    ctx.write_shape_named(name, shape)
        .ok_or_else(|| "write_step_with_name: nothing to write".to_string())?;
    Ok(ctx.finish())
}

/// Serialize a single shape with a minimal STEP style block.
///
/// A `COLOUR_RGB` and the `SURFACE_STYLE_FILL_AREA` / `SURFACE_STYLE_USAGE`
/// style chain are attached to the shape's representation through a
/// `STYLED_ITEM` record. The style is decorative: STEP readers that do not
/// interpret presentation styles simply skip these records.
pub fn write_step_with_color(shape: &TopoShape, color: (f64, f64, f64)) -> Result<String, String> {
    let mut ctx = WriteCtx::new();
    let rep = ctx
        .write_shape_named("Shape", shape)
        .ok_or_else(|| "write_step_with_color: nothing to write".to_string())?;
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
