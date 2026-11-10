

use occt_core::gp::{
    GpTrsf, TrsfForm,
};
use occt_core::toploc::TopLocLocation;

use occt_topo::shape::TopoShape;

use super::super::iges_to_brep_basic_curve::{
    convert_location, entity_location,
};
use super::super::iges_to_brep_curve_and_surface::Session;

/// `IGESToBRep_TopoSurface.cxx:1647-1658`: `SetEpsilon(1.E-04)` then
/// `IGESData_ToolLocation::ConvertLocation`. A failed conversion leaves the
/// default (identity) matrix in place and only warns.
pub(super) fn plane_trsf(session: &Session<'_>, index: usize) -> Result<GpTrsf, String> {
    if let Some(g) = entity_location(session.model, index)? {
        if let Some(tr) = convert_location(1.0e-4, &g, session.model.unit_factor) {
            return Ok(tr);
        }
    }
    Ok(GpTrsf::identity())
}

/// `TopLoc_Location loc(trsf); shape.Location(loc, false)` — the matrix
/// `TransferPlaneParts` computed (`IGESToBRep_TopoSurface.cxx:1511-1515`,
/// `:1597-1603`).
pub(super) fn apply_plane_trsf(mut shape: TopoShape, trsf: &GpTrsf) -> TopoShape {
    if trsf.form() != TrsfForm::Identity {
        shape.move_location(&TopLocLocation::composed(trsf, TopLocLocation::identity()));
    }
    shape
}
