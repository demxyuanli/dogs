use super::super::prelude::*;
use super::super::*;

/// Shared assembly serialization; `splines` selects the geometry classifier.
pub(in crate::step) fn write_assembly_inner(a: &StepAssembly, splines: bool) -> Result<String, String> {
    let mut ctx = WriteCtx::new();
    ctx.splines = splines;
    let mut defs: HashMap<String, usize> = HashMap::new();

    // The assembly's own product definition is the root that children hang off.
    let top = esc_str(&a.name);
    let top_product = ctx.w.emit(format!("PRODUCT('{top}','{top}','',(#{}))", ctx.prod_ctx));
    let top_form = ctx.w.emit(format!("PRODUCT_DEFINITION_FORMATION('','',#{top_product})"));
    let top_def = ctx.w.emit(format!("PRODUCT_DEFINITION('','',#{top_form},#{})", ctx.def_ctx));
    defs.insert(a.name.clone(), top_def);

    for (name, shape) in &a.products {
        let n = esc_str(name);
        let product = ctx.w.emit(format!("PRODUCT('{n}','{n}','',(#{}))", ctx.prod_ctx));
        let formation = ctx.w.emit(format!("PRODUCT_DEFINITION_FORMATION('','',#{product})"));
        let def = ctx.w.emit(format!("PRODUCT_DEFINITION('','',#{formation},#{})", ctx.def_ctx));
        let items = ctx.emit_top(shape);
        if !items.is_empty() {
            let rep = ctx.w.emit(format!(
                "ADVANCED_BREP_SHAPE_REPRESENTATION('{n}',({}),#{})",
                join_refs(&items),
                ctx.geom_ctx
            ));
            let pds = ctx.w.emit(format!("PRODUCT_DEFINITION_SHAPE('','',#{def})"));
            ctx.w
                .emit(format!("SHAPE_DEFINITION_REPRESENTATION(#{pds},#{rep})"));
        }
        defs.insert(name.clone(), def);
    }

    for (parent, kids) in &a.children {
        let pdef = defs
            .get(parent)
            .copied()
            .ok_or_else(|| format!("write_step_assembly: unknown parent '{parent}'"))?;
        for kid in kids {
            let kdef = defs
                .get(kid)
                .copied()
                .ok_or_else(|| format!("write_step_assembly: unknown child '{kid}'"))?;
            let kn = esc_str(kid);
            ctx.w.emit(format!(
                "NEXT_ASSEMBLY_USAGE_OCCURRENCE('{kn}','{kn}','',#{pdef},#{kdef},'')"
            ));
        }
    }
    Ok(ctx.finish())
}

/// Serialize a model to STEP, emitting B-spline geometry for non-analytic
/// curves and surfaces (the spline-capable counterpart of [`write_step`]).
pub fn write_step_model_with_splines(model: &BRepModel) -> Result<String, String> {
    let mut ctx = WriteCtx::new();
    ctx.splines = true;
    for ms in &model.shapes {
        ctx.write_shape_named(&ms.name, &ms.shape);
    }
    Ok(ctx.finish())
}

/// Serialize a model to STEP, attaching a `COLOUR_RGB` style to each shape
/// whose `ModelShape.color` is set.
///
/// Shapes without a color are written plain. This is the model-level
/// counterpart of [`write_step_with_color`], reading the colors already
/// stored on the [`BRepModel`] (e.g. via `add_with_color`).
pub fn write_step_model_with_colors(model: &BRepModel) -> Result<String, String> {
    let mut ctx = WriteCtx::new();
    for ms in &model.shapes {
        let rep = ctx.write_shape_named(&ms.name, &ms.shape);
        if let (Some(rep), Some(c)) = (rep, ms.color) {
            let color = ctx.w.emit(format!(
                "COLOUR_RGB('',{},{},{})",
                step_real(c.r as f64),
                step_real(c.g as f64),
                step_real(c.b as f64)
            ));
            let fill = ctx.w.emit(format!("SURFACE_STYLE_FILL_AREA('',#{color})"));
            let usage = ctx.w.emit(format!("SURFACE_STYLE_USAGE('',#{fill})"));
            ctx.w.emit(format!("STYLED_ITEM('',(#{usage}),#{rep})"));
        }
    }
    Ok(ctx.finish())
}

/// Serialize a model with per-shape [`StepWriteOptions`].
///
/// Each shape in the model is written with the shared options; spline output
/// applies to every shape, and the color/units/name defaults apply as in
/// [`write_step_with_options`]. Shape names come from the model, so `opts.name`
/// is only used as a fallback for unnamed shapes.
pub fn write_step_model_with_options(model: &BRepModel, opts: &StepWriteOptions) -> Result<String, String> {
    let mut ctx = WriteCtx::new();
    ctx.splines = opts.splines;
    for ms in &model.shapes {
        let name = if ms.name.is_empty() { opts.name.clone() } else { ms.name.clone() };
        let rep = ctx.write_shape_named(&name, &ms.shape);
        if let (Some(rep), Some((r, g, b))) = (rep, opts.color) {
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
    }
    Ok(ctx.finish())
}

/// Read a STEP physical file from disk, also returning collected warnings for
/// skipped or unsupported entities.
pub fn read_step_file_with_warnings(path: &str) -> Result<(BRepModel, Vec<String>), String> {
    let content = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
    read_step_with_warnings(&content)
}

/// Serialize a list of named shapes into one STEP file.
///
/// Each `(name, shape)` pair becomes its own representation + product pair,
/// exactly as [`write_step`] does for a [`BRepModel`], without requiring the
/// caller to build a model object first. The reader returns each shape under
/// its `name`.
pub fn write_step_shapes(shapes: &[(String, TopoShape)]) -> String {
    let mut ctx = WriteCtx::new();
    for (name, shape) in shapes {
        ctx.write_shape_named(name, shape);
    }
    ctx.finish()
}

/// Spline-capable counterpart of [`write_step_shapes`]: non-analytic curves
/// and surfaces are emitted as B-splines.
pub fn write_step_shapes_with_splines(shapes: &[(String, TopoShape)]) -> Result<String, String> {
    let mut ctx = WriteCtx::new();
    ctx.splines = true;
    for (name, shape) in shapes {
        ctx.write_shape_named(name, shape);
    }
    Ok(ctx.finish())
}

/// Serialize a compound as its named children in one STEP file.
///
/// A [`TopoShape`] of type `Compound` is flattened: each child shape is written
/// as a named representation (`Child1`, `Child2`, ...), which the reader
/// reconstructs as separate model shapes. Non-compound shapes are written as a
/// single `"Shape"` representation.
pub fn write_step_compound(compound: &TopoShape) -> Result<String, String> {
    if !compound.is_compound() {
        return write_step_with_splines(compound);
    }
    let kids: Vec<TopoShape> = compound
        .tshape
        .read()
        .unwrap()
        .children
        .clone();
    if kids.is_empty() {
        return Err("write_step_compound: empty compound".into());
    }
    let named: Vec<(String, TopoShape)> = kids
        .into_iter()
        .enumerate()
        .map(|(i, k)| (format!("Child{}", i + 1), k))
        .collect();
    Ok(write_step_shapes_with_splines(&named)?)
}

