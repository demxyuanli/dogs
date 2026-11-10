
use super::super::interference::Interference;
use super::super::interference_iterator::InterferenceList;



/// `FDS_Parameter(I, par)` (`cxx:464-475`): true only for an
/// `EdgeVertexInterference` or a `CurvePointInterference`.
pub fn fds_parameter(i: &Interference) -> Option<f64> {
    match i.data() {
        super::super::interference::InterferenceData::CurvePoint { .. }
        | super::super::interference::InterferenceData::EdgeVertex(_) => i.parameter(),
        _ => None,
    }
}

/// `FDS_Parameter(I)` (`cxx:478-495`): the `Standard_Failure` for any other
/// subclass is `None`.
pub fn fds_parameter_or_fail(i: &Interference) -> Result<f64, String> {
    fds_parameter(i).ok_or_else(|| "FDS_Parameter".to_string())
}

/// `FUN_selectTRAINTinterference(li, liINTERNAL)` (`cxx:696-712`): collect the
/// interferences whose `Orientation(TopAbs_IN)` is `INTERNAL`.
pub fn fun_select_train_tinterference(
    li: &InterferenceList,
    li_internal: &mut InterferenceList,
) -> bool {
    li_internal.clear();
    for i in li {
        let ori = i
            .read()
            .expect("interference lock")
            .transition()
            .orientation(crate::abs::State::In, crate::abs::ShapeType::Face);
        if ori == crate::abs::Orientation::Internal {
            li_internal.push(i.clone());
        }
    }
    !li_internal.is_empty()
}
