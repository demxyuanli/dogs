use super::super::prelude::*;
use super::super::*;


// ---------------------------------------------------------------------------
// Read precision: `STEPControl_ActorRead::myPrecision`
// ---------------------------------------------------------------------------

/// `Interface_StaticStandards.cxx:33-37`: `read.precision.mode` is initialised
/// `'e' ""` then `SetIVal(..., 0)`, so the default mode is 0.
const READ_PRECISION_MODE: i32 = 0;
/// `Interface_StaticStandards.cxx:39`: `read.precision.val 'r' "1.e-03"`.
const READ_PRECISION_VAL: f64 = 1.0e-3;
/// `StepData_Factors.cxx:24` `myCascadeUnit(1.)` - "length unit for current
/// transfer process (mm by default)".
const CASCADE_UNIT: f64 = 1.0;

/// `STEPConstruct_UnitContext::ConvertSiPrefix` (`STEPConstruct_UnitContext.cxx:210-250`).
/// An unknown prefix falls through to `1.` (`cxx:246-249`).
fn convert_si_prefix(prefix: &str) -> f64 {
    match prefix.trim().trim_matches('.').to_ascii_uppercase().as_str() {
        "EXA" => 1.0e18,
        "PETA" => 1.0e15,
        "TERA" => 1.0e12,
        "GIGA" => 1.0e9,
        "MEGA" => 1.0e6,
        "KILO" => 1.0e3,
        "HECTO" => 1.0e2,
        "DECA" => 1.0e1,
        "DECI" => 1.0e-1,
        "CENTI" => 1.0e-2,
        "MILLI" => 1.0e-3,
        "MICRO" => 1.0e-6,
        "NANO" => 1.0e-9,
        "PICO" => 1.0e-12,
        "FEMTO" => 1.0e-15,
        "ATTO" => 1.0e-18,
        _ => 1.0,
    }
}

/// `STEPConstruct_UnitContext::ComputeFactors` (`cxx:306-450`) for a unit whose
/// complex body is a `ConversionBasedUnitAndLengthUnit` or a
/// `SiUnitAndLengthUnit`: `lengthFactor = parameter * 1000. / aCascadeUnit`
/// (`cxx:424-428`, the `METER` build option is not defined, so the scaled
/// branch is the live one). Returns `None` when the members are not a length
/// unit.
///
/// `cxx:322-381` (conversion based, e.g. INCH): `theFactor = theSIPFactor *
/// theMVAL`, the conversion measure times the SI prefix of the unit it converts
/// to (`INCH` -> `LENGTH_MEASURE(25.4)` on `.MILLI.METRE.` -> `25.4`).
/// `cxx:384-417` (SI): `theFactor = theSIPFactor * theSIUNF` (the SI name
/// factor is `1.`). `cxx:344-348` / `cxx:370`: a missing conversion factor or a
/// non-SI target yields no parameter (`ComputeFactors` returns `-1` / `3` and
/// the pre-set default stays in place).
fn length_unit_factor(
    members: &[(String, Vec<String>)],
    records: &HashMap<usize, Record>,
) -> Option<f64> {
    if !members.iter().any(|(t, _)| t == "LENGTH_UNIT") {
        return None;
    }
    let parameter =
        if let Some((_, args)) = members.iter().find(|(t, _)| t == "CONVERSION_BASED_UNIT") {
            let conv = args.get(1).and_then(|a| parse_ref(a))?;
            let conv = records.get(&conv)?;
            let value = conv.args.first().and_then(|a| measure_value(a))?;
            let target = conv.args.get(1).and_then(|a| parse_ref(a))?;
            si_unit_member_scale(&records.get(&target)?.members)? * value
        } else {
            si_unit_member_scale(members)?
        };
    Some(parameter * 1000.0 / CASCADE_UNIT)
}

/// `STEPConstruct_UnitContext::ComputeFactors` (`cxx:384-417`) for the
/// `SI_UNIT(prefix, name)` member of a unit's complex body: `theFactor =
/// theSIPFactor * theSIUNF` (`cxx:405`). `SiUnitNameFactor` (`cxx:254-268`)
/// assigns `theSIUNFactor = 1.` before its switch, so `theSIUNF` is `1.` for
/// every name; an unrecognised name only sets `status = 11` at `cxx:397-400`
/// (`cxx:361-365` for the conversion-based arm) and computation continues, so
/// the prefix factor is returned for any name. `None` when the unit has no
/// `SI_UNIT` member at all, which is the `return 3` branch of `cxx:368-371`.
fn si_unit_member_scale(members: &[(String, Vec<String>)]) -> Option<f64> {
    let (_, args) = members.iter().find(|(t, _)| t == "SI_UNIT")?;
    let prefix = args.first().map(String::as_str).unwrap_or("$");
    Some(if prefix.trim() == "$" {
        1.0
    } else {
        convert_si_prefix(prefix)
    })
}

/// `STEPConstruct_UnitContext::ComputeFactors` (`cxx:306-450`) for a unit whose
/// complex body is a `ConversionBasedUnitAndPlaneAngleUnit` or a
/// `SiUnitAndPlaneAngleUnit` (`cxx:439-444`): the returned parameter is the
/// `planeAngleFactor`. `None` when the unit is not a plane angle unit.
///
/// `cxx:322-381` (conversion based): `theFactor = theSIPFactor * theMVAL`, the
/// conversion measure times the SI prefix of the unit it converts to
/// (`DEGREE` -> `PLANE_ANGLE_MEASURE(0.01745329252)` on `.RADIAN.`).
/// `cxx:344-348` / `cxx:370`: a missing conversion factor or a non-SI target
/// yields no parameter (`ComputeFactors` returns `-1` / `3` and the pre-set
/// default stays in place).
fn plane_angle_unit_factor(
    members: &[(String, Vec<String>)],
    records: &HashMap<usize, Record>,
) -> Option<f64> {
    if !members.iter().any(|(t, _)| t == "PLANE_ANGLE_UNIT") {
        return None;
    }
    if let Some((_, args)) = members.iter().find(|(t, _)| t == "CONVERSION_BASED_UNIT") {
        let conv = args.get(1).and_then(|a| parse_ref(a))?;
        let conv = records.get(&conv)?;
        let value = conv.args.first().and_then(|a| measure_value(a))?;
        let target = conv.args.get(1).and_then(|a| parse_ref(a))?;
        return Some(si_unit_member_scale(&records.get(&target)?.members)? * value);
    }
    si_unit_member_scale(members)
}

/// `STEPControl_ActorRead::PrepareUnits` -> `STEPConstruct_UnitContext::
/// ComputeFactors(theGUAC, ...)` (`STEPControl_ActorRead.cxx:2312-2355`): the
/// plane angle factor of the `GLOBAL_UNIT_ASSIGNED_CONTEXT`. It is the value
/// `StepData_Factors::PlaneAngleFactor()` hands to
/// `StepToGeom::MakeConicalSurface` (`StepToGeom.cxx:1316`) and the value
/// `FactorDegreeRadian()` hands to `GeomConvert_Units::DegreeToRadian`
/// (`StepToTopoDS_TranslateEdge.cxx:579`).
///
/// `ComputeFactors` pre-sets `planeAngleFactor = PI/180.` (`cxx:281`) and leaves
/// it when the context lists no plane angle unit, so a context that parses but
/// has no angle unit keeps `PI/180.`. A context whose unit list does not resolve
/// is the "Bad RepresentationContext, default unit taken" branch
/// (`STEPControl_ActorRead.cxx:2300-2305`), where `ResetUnits` (`cxx:2400`)
/// leaves every factor at `1.`. The port resolves one flat record map instead of
/// per-representation factors, so the first unit context in the file is used.
pub(in crate::step) fn context_plane_angle_factor(records: &HashMap<usize, Record>) -> f64 {
    let mut usable_context = false;
    for rec in records.values() {
        let Some((_, units)) = rec
            .members
            .iter()
            .find(|(t, _)| t == "GLOBAL_UNIT_ASSIGNED_CONTEXT")
        else {
            continue;
        };
        for unit_id in units.iter().flat_map(|a| parse_ref_list(a)) {
            let Some(unit) = records.get(&unit_id) else {
                continue;
            };
            usable_context = true;
            if let Some(f) = plane_angle_unit_factor(&unit.members, records) {
                return f;
            }
        }
    }
    if usable_context {
        PI / 180.0
    } else {
        1.0
    }
}

/// `STEPConstruct_UnitContext::ComputeFactors(theGUAC, ...)` (`cxx:272-302`):
/// the length factor comes from the length unit listed by the representation's
/// `GLOBAL_UNIT_ASSIGNED_CONTEXT`. `lengthFactor` is pre-set to `1.` (`cxx:280`).
pub(in crate::step) fn context_length_factor(records: &HashMap<usize, Record>) -> f64 {
    for rec in records.values() {
        let Some((_, units)) = rec
            .members
            .iter()
            .find(|(t, _)| t == "GLOBAL_UNIT_ASSIGNED_CONTEXT")
        else {
            continue;
        };
        for unit_id in units.iter().flat_map(|a| parse_ref_list(a)) {
            if let Some(unit) = records.get(&unit_id) {
                if let Some(f) = length_unit_factor(&unit.members, records) {
                    return f;
                }
            }
        }
    }
    1.0
}

/// `LENGTH_MEASURE(2.E-005)` -> `2.E-005`.
fn measure_value(arg: &str) -> Option<f64> {
    let open = arg.find('(')?;
    let close = arg.rfind(')')?;
    arg[open + 1..close].trim().parse().ok()
}

/// `STEPControl_ActorRead::myPrecision` (`STEPControl_ActorRead.cxx:2370-2384`),
/// the value `StepToTopoDS_TranslateEdgeLoop::Precision()` returns (`cxx:236`)
/// and that `ShapeFix_EdgeProjAux::Compute` (`cxx:844`) and
/// `XSAlgo_ShapeProcessor::CheckPCurve` (`cxx:875` -> `cxx:175`) receive.
///
///   if (ReadPrecisionMode == 1)        myPrecision = ReadPrecisionVal;
///   else if (myUnit.HasUncertainty())  myPrecision = Uncertainty() * LengthFactor();
///   else                               myPrecision = ReadPrecisionVal;
///
/// `myUnit` is the local `STEPConstruct_UnitContext` filled by
/// `STEPControl_ActorRead::PrepareUnits` (`cxx:2312-2384`): `ComputeFactors`
/// (`cxx:2347`) sets the length factor and `ComputeTolerance` (`cxx:2363`) the
/// uncertainty.
pub(in crate::step) fn step_precision(records: &HashMap<usize, Record>) -> f64 {
    if READ_PRECISION_MODE == 1 {
        return READ_PRECISION_VAL;
    }
    // `STEPConstruct_UnitContext::ComputeTolerance` (`cxx:480-541`): only a
    // `UNCERTAINTY_MEASURE_WITH_UNIT` whose unit component is an SI (`cxx:504`)
    // or conversion-based (`cxx:521`) length unit counts, and `theUncertainty`
    // (initialised to `RealLast`) is only ever lowered.
    let mut uncertainty = f64::MAX;
    let mut has_uncertainty = false;
    for rec in records.values() {
        if rec.type_name != "UNCERTAINTY_MEASURE_WITH_UNIT" {
            continue;
        }
        let Some(value) = rec.args.first().and_then(|a| measure_value(a)) else {
            continue;
        };
        let Some(unit_id) = rec.args.get(1).and_then(|a| parse_ref(a)) else {
            continue;
        };
        let Some(unit) = records.get(&unit_id) else {
            continue;
        };
        if !unit.members.iter().any(|(t, _)| t == "LENGTH_UNIT") {
            continue;
        }
        if uncertainty > value {
            uncertainty = value;
        }
        has_uncertainty = true;
    }
    if !has_uncertainty {
        return READ_PRECISION_VAL;
    }
    uncertainty * context_length_factor(records)
}
