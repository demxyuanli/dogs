use super::super::prelude::*;
use super::super::*;



/// Write a shape to a STEP file on disk, using the spline-capable writer.
pub fn write_step_file_with_splines(path: &str, shape: &TopoShape) -> std::io::Result<()> {
    let content = write_step_with_splines(shape).map_err(|e| std::io::Error::other(e))?;
    std::fs::write(path, content)
}

/// Write an assembly to a STEP file on disk.
pub fn write_step_assembly_file(path: &str, a: &StepAssembly) -> std::io::Result<()> {
    let content = write_step_assembly(a).map_err(|e| std::io::Error::other(e))?;
    std::fs::write(path, content)
}

/// A lightweight assembly description: named parts plus a parent→children tree.
///
/// The top-level `name` is emitted as the assembly's own product definition,
/// so assembly trees may reference it as a parent in `children`.
#[derive(Debug, Clone)]
pub struct StepAssembly {
    pub name: String,
    pub products: Vec<(String, TopoShape)>,
    pub children: Vec<(String, Vec<String>)>,
}

/// Serialize an assembly: a `PRODUCT` / `PRODUCT_DEFINITION` pair per part
/// (with the part's shape representation), plus `NEXT_ASSEMBLY_USAGE_OCCURRENCE`
/// records linking the product definitions along the assembly tree.
///
/// Every name in `children` (both parents and children) must be either the
/// assembly `name` or a key of `products`; unknown names return an error.
pub fn write_step_assembly(a: &StepAssembly) -> Result<String, String> {
    write_assembly_inner(a, false)
}

/// Serialize an assembly with B-spline output for non-analytic part geometry.
///
/// Equivalent to [`write_step_assembly`] with the spline-capable geometry
/// classification enabled, so parts whose faces/edges are B-spline (or
/// otherwise non-analytic) round-trip through `B_SPLINE_SURFACE_WITH_KNOTS` /
/// `B_SPLINE_CURVE_WITH_KNOTS`.
pub fn write_step_assembly_with_splines(a: &StepAssembly) -> Result<String, String> {
    write_assembly_inner(a, true)
}
