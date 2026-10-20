//! `ComputeSRRWT` (`STEPConstruct_Styles` / `STEPControl_ActorRead`).

use super::super::prelude::*;
use super::*;

/// `STEPControl_ActorRead::ComputeSRRWT` (`STEPControl_ActorRead.cxx:2495-2546`).
///
/// Only `SHAPE_REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION` carries a
/// transformation: either a `CARTESIAN_TRANSFORMATION_OPERATOR_3D`
/// (`StepToGeom::MakeTransformation3d`, `cxx:2522`) or an
/// `ITEM_DEFINED_TRANSFORMATION` whose two `AXIS2_PLACEMENT_3D` items go through
/// `ComputeTransformation` (`cxx:2545`, axis-pair arm at `cxx:2484-2489`).
/// Both return "no transformation" when the result is the identity
/// (`cxx:2527`, `cxx:2489`).
///
/// **UNPORTED**: the `PrepareUnits(Rep2)` / restore pair (`cxx:2517-2526`) and
/// the `PrepareUnits(OrigContext/TargContext)` calls of `ComputeTransformation`
/// (`cxx:2467-2482`) — the per-representation unit context; the port uses the
/// file's global factors. Also UNPORTED: the axis-ownership / swapped-axis
/// warnings of `cxx:2440-2464`.
pub(super) fn compute_srrwt(
    resolver: &Resolver,
    records: &HashMap<usize, Record>,
    srr_id: usize,
) -> Option<occt_core::gp::GpTrsf> {
    let rec = records.get(&srr_id)?;
    if rec.type_name != "SHAPE_REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION" {
        return None;
    }
    let op_ref = parse_ref(rec.args.get(4)?)?;
    compute_srrwt_from_op(resolver, records, op_ref)
}

/// The `ComputeSRRWT` body above, starting from the operator entity: the same
/// operator reaches it through a `CDSR`'s complex relationship
/// (`read_assembly::nauo_target`), where the `SRRWT` record itself does not
/// exist as a standalone entity.
pub(in crate::step) fn compute_srrwt_from_op(
    resolver: &Resolver,
    records: &HashMap<usize, Record>,
    op_ref: usize,
) -> Option<occt_core::gp::GpTrsf> {
    let op = records.get(&op_ref)?;
    match op.type_name.as_str() {
        "CARTESIAN_TRANSFORMATION_OPERATOR_3D" => {
            let t = resolver.make_transformation3d(op_ref).ok()?;
            (t.form() != occt_core::gp::TrsfForm::Identity).then_some(t)
        }
        "ITEM_DEFINED_TRANSFORMATION" => {
            let ax1 = self_axis(resolver, op.args.get(2))?;
            let ax2 = self_axis(resolver, op.args.get(3))?;
            let t = resolver.compute_axis_transform(&ax1, &ax2);
            (t.form() != occt_core::gp::TrsfForm::Identity).then_some(t)
        }
        _ => None,
    }
}

fn self_axis(resolver: &Resolver, arg: Option<&String>) -> Option<GpAx2> {
    resolver.resolve_axis2(parse_ref(arg?)?).ok()
}
