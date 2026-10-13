//! `IGESToBRep_Actor::Recognize` and `Transfer` (`IGESToBRep_Actor.cxx:135-237`).
//!
//! After `TransferGeometry`, `ProcessShape` runs the ported `FixShape` steps
//! (`IGESToBRep_Actor.cxx:204-228`). See `iges_process_shape`.

use crate::shape::TopoShape;

use super::iges_data_iges_model::{IgesEntity, IgesModel};
use super::iges_process_shape::process_iges_shape;
use super::iges_to_brep_curve_and_surface::{is_curve_and_surface, transfer_geometry, Session};

/// `IGESToBRep_Actor::Recognize` (`IGESToBRep_Actor.cxx:116-131`). Needs the
/// model: `IsCurveAndSurface` inspects the entity (402 form 9's members).
pub(in crate::iges_read) fn actor_recognize(model: &IgesModel, index: usize) -> bool {
    // `cxx:120-123`: a null entity is not recognized.
    let Ok(ent) = model.entity(index) else {
        return false;
    };
    // `cxx:128-130`.
    if is_curve_and_surface(model, index) {
        return true;
    }
    if ent.type_number == 402 && matches!(ent.form_number, 1 | 7 | 14 | 15) {
        return true;
    }
    ent.type_number == 408 || ent.type_number == 308
}

pub(in crate::iges_read) fn actor_transfer(
    model: &IgesModel,
    index: usize,
) -> Result<Option<TopoShape>, String> {
    let ent = model.entity(index)?;
    if !actor_recognize(model, index) {
        return Ok(None);
    }
    // `read.precision.mode` default 0 uses `GlobalSection::Resolution`
    // (`IGESToBRep_Actor.cxx:171-175`). Not used as a length filter.
    let _eps = model.resolution;
    let mut session = Session::new(model);
    let shape = transfer_geometry(&mut session, index)?;
    Ok(shape.map(|s| process_iges_shape(&s, model)))
}
