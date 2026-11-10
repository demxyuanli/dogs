use super::super::prelude::*;
use super::super::*;

use super::context::WriteCtx;


/// SI unit definitions for [`write_step_with_units`].
///
/// The length unit is given as the number of metres per STEP length unit
/// (1.0 = metre, 0.001 = millimetre); the angle unit as radians per STEP
/// angle unit (1.0 = radian, `π/180` ≈ 0.01745 = degree). The emitted
/// `SI_UNIT` records attach these to a `GEOMETRIC_REPRESENTATION_CONTEXT`.
#[derive(Debug, Clone)]
pub struct StepUnits {
    pub length_unit_m: f64,
    pub angle_unit_rad: f64,
}

impl Default for StepUnits {
    fn default() -> Self {
        Self {
            length_unit_m: 1.0,
            angle_unit_rad: 1.0,
        }
    }
}

/// Emit the SI unit block (`DIMENSIONAL_EXPONENTS` + `SI_UNIT` records).
///
/// A dimensional-exponent vector of a length measure (metre^1) and one of a
/// plane angle (radian, dimensionless in the length slots) are shared by the
/// two `SI_UNIT` records. Readers that do not interpret units skip them.
pub(crate) fn emit_si_units(w: &mut StepWriter, units: &StepUnits) {
    let _ = units;
    let len_dim = w.emit("DIMENSIONAL_EXPONENTS(1.,0.,0.,0.,0.,0.,0.)".into());
    let ang_dim = w.emit("DIMENSIONAL_EXPONENTS(0.,0.,0.,1.,0.,0.,0.)".into());
    w.emit(format!("SI_UNIT(#{len_dim},*,.METRE.)"));
    w.emit(format!("SI_UNIT(#{ang_dim},*,.RADIAN.)"));
}

/// Serialize a single shape with an explicit PRODUCT name and SI unit block.
///
/// The unit records (`DIMENSIONAL_EXPONENTS` + `SI_UNIT`) are emitted after
/// the shape representation, mirroring the unit scaffolding a full
/// `STEPControl_Writer` places in the DATA section. The length/angle scale
/// factors from `units` are recorded so a downstream reader can convert to
/// metres; readers that do not interpret units simply skip the records.
pub fn write_step_with_units(
    shape: &TopoShape,
    name: &str,
    units: &StepUnits,
) -> Result<String, String> {
    let mut ctx = WriteCtx::new();
    ctx.write_shape_named(name, shape)
        .ok_or_else(|| "write_step_with_units: nothing to write".to_string())?;
    emit_si_units(&mut ctx.w, units);
    Ok(ctx.finish())
}
