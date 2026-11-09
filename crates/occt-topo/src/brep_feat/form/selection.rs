//! The selection of the tool parts to keep (`BRepFeat_Form.cxx:812-1211`).
//!
//! The parametric window comes from `BRepFeat::ParametricMinMax` over the
//! feature axis: the parts of the tool whose own interval overlaps that window
//! are kept. The `sens` of the intersection scan decides which end of the
//! window is open.

use std::sync::Arc;

use occt_core::elib::clib::in_period;
use occt_core::precision;
use occt_geom::Curve;

use crate::abs::ShapeType;
use crate::brep_builder_api::make_compound_shapes;
use crate::brep_feat::builder::BRepFeatBuilder;
use crate::brep_feat::perf_selection::PerfSelection;
use crate::brep_feat::statics::parametric_min_max;
use crate::brep_feat::status_error::StatusError;
use crate::locope::ShapeSet;
use crate::shape::TopoShape;
use crate::topo_tools::faces_of;

use super::descendants::descendants;
use super::{BRepFeatForm, FormCurves};

/// The parametric window plus the `flag` of the last `ParametricMinMax` call.
struct Window {
    prmin: f64,
    prmax: f64,
    pbmin: f64,
    pbmax: f64,
    flag: bool,
}

/// `cxx:812-1211`: `KeepPart` the parts of the tool inside the window.
///
/// `Ok(false)` is one of the three early `return`s of the reference (`cxx:827`,
/// `:876`, `:1173`), where `GlobalPerform` leaves without touching the result;
/// the caller performs that return.
pub(super) fn select_parts(
    form: &mut BRepFeatForm,
    curves: &mut dyn FormCurves,
    the_builder: &mut BRepFeatBuilder,
    lshape: &[TopoShape],
    sens_in: i32,
) -> Result<bool, String> {
    let mut sens = sens_in;

    // `cxx:819`.
    if lshape.is_empty() || form.my_perf_selection == PerfSelection::NoSelection {
        return Ok(true);
    }

    // `cxx:822-828`.
    let curve = match curves.baryc_curve() {
        Some(c) => c,
        None => {
            form.my_status_error = StatusError::EmptyBaryCurve;
            form.not_done();
            return Ok(false);
        }
    };

    let delta = precision::CONFUSION;
    let real_last = f64::MAX;
    let real_first = -f64::MAX;

    // `cxx:813-817`.
    let mut window = Window {
        prmin: real_last - 2.0 * delta,
        prmax: real_first + 2.0 * delta,
        pbmin: real_last,
        pbmax: real_first,
        flag: false,
    };

    // `cxx:830-928`: the window of each selection kind.
    let empty = TopoShape::new(ShapeType::Compound);
    if form.my_perf_selection == PerfSelection::SelectionSh {
        let base = form.my_sbase.as_ref().unwrap_or(&empty);
        let (a, b, c, d, f) = parametric_min_max(base, &curve, false);
        window.prmin = a;
        window.prmax = b;
        window.pbmin = c;
        window.pbmax = d;
        window.flag = f;
    } else if form.my_perf_selection == PerfSelection::SelectionFU {
        // `cxx:839-866`: both calls write the same `flag1`, so the second call
        // overwrites the first one's value.
        let sfrom = form.my_sfrom.as_ref().unwrap_or(&empty);
        let (prmin1, prmax1, prbmin1, prbmax1, _f1) = parametric_min_max(sfrom, &curve, false);
        let suntil = form.my_suntil.as_ref().unwrap_or(&empty);
        let (prmin2, prmax2, prbmin2, prbmax2, f2) = parametric_min_max(suntil, &curve, false);
        window.flag = f2;

        if curve.is_periodic() {
            let period = curve.period();
            window.prmax = prmax2;
            if window.flag {
                window.prmin = in_period(prmin1, window.prmax - period, window.prmax);
            } else {
                window.prmin = prmin1.min(prmin2);
            }
            window.pbmax = prbmax2;
            window.pbmin = in_period(prbmin1, window.pbmax - period, window.pbmax);
        } else {
            window.prmin = prmin1.min(prmin2);
            window.prmax = prmax1.max(prmax2);
            window.pbmin = prbmin1.min(prbmin2);
            window.pbmax = prbmax1.max(prbmax2);
        }
    } else if form.my_perf_selection == PerfSelection::SelectionShU {
        // `cxx:868-881`.
        if !form.my_just_feat && sens == 0 {
            sens = 1;
        }
        if sens == 0 {
            form.my_status_error = StatusError::IncDirection;
            form.not_done();
            return Ok(false);
        }

        let suntil = form.my_suntil.as_ref().unwrap_or(&empty);
        let (prmin1, prmax1, prbmin1, prbmax1, _f1) = parametric_min_max(suntil, &curve, false);
        let base = form.my_sbase.as_ref().unwrap_or(&empty);
        let (prmin2, prmax2, prbmin2, prbmax2, f2) = parametric_min_max(base, &curve, false);
        window.flag = f2;

        if sens == 1 {
            window.prmin = prmin2;
            window.prmax = prmax1;
            window.pbmin = prbmin2;
            window.pbmax = prbmax1;
        } else if sens == -1 {
            window.prmin = prmin1;
            window.prmax = prmax2;
            window.pbmin = prbmin1;
            window.pbmax = prbmax2;
        }
    } else if form.my_perf_selection == PerfSelection::SelectionU {
        // `cxx:901-910`.
        if sens == 0 {
            form.my_status_error = StatusError::IncDirection;
            form.not_done();
            return Ok(false);
        }

        let suntil = form.my_suntil.as_ref().unwrap_or(&empty);
        let (prmin1, prmax1, prbmin1, prbmax1, f1) = parametric_min_max(suntil, &curve, false);
        window.flag = f1;

        if sens == 1 {
            window.prmin = real_first;
            window.prmax = prmax1;
            window.pbmin = real_first;
            window.pbmax = prbmax1;
        } else if sens == -1 {
            window.prmin = prmin1;
            window.prmax = real_last;
            window.pbmin = prbmin1;
            window.pbmax = real_last;
        }
    }

    // `cxx:934-1112`: narrow the window where the tool meets From / Until.
    if form.my_perf_selection != PerfSelection::NoSelection {
        if let Some(suntil) = form.my_suntil.clone() {
            let mut map_funtil = ShapeSet::new();
            descendants(&suntil, the_builder, &mut map_funtil);
            if !map_funtil.is_empty() {
                for part in lshape {
                    for face in faces_of(part) {
                        if map_funtil.contains(&face.0) {
                            let (prmin1, prmax1, prbmin1, prbmax1, flag3) =
                                parametric_min_max(&face.0, &curve, false);
                            let (prmin2, prmax2, prbmin2, prbmax2, flag2) =
                                parametric_min_max(part, &curve, false);
                            narrow_until(
                                &mut window,
                                sens,
                                delta,
                                (prmin1, prmax1, prbmin1, prbmax1, flag3),
                                (prmin2, prmax2, prbmin2, prbmax2, flag2),
                            );
                            break;
                        }
                    }
                }
            }
        }
        if let Some(sfrom) = form.my_sfrom.clone() {
            let mut map_ffrom = ShapeSet::new();
            descendants(&sfrom, the_builder, &mut map_ffrom);
            if !map_ffrom.is_empty() {
                for part in lshape {
                    for face in faces_of(part) {
                        if map_ffrom.contains(&face.0) {
                            let (prmin1, prmax1, prbmin1, prbmax1, flag3) =
                                parametric_min_max(&face.0, &curve, false);
                            let (prmin2, prmax2, prbmin2, prbmax2, flag2) =
                                parametric_min_max(part, &curve, false);
                            narrow_from(
                                &mut window,
                                sens,
                                delta,
                                (prmin1, prmax1, prbmin1, prbmax1, flag3),
                                (prmin2, prmax2, prbmin2, prbmax2, flag2),
                            );
                            break;
                        }
                    }
                }
            }
        }
    }

    // `cxx:1114-1177`: keep the overlapping parts, or fail.
    if !form.my_just_feat {
        let mut keep_parts = false;
        for part in lshape {
            let (prmin1, prmax1, prbmin1, prbmax1, flag2) = part_window(part, &curve);
            let (pmin, pmax, min, max) = if !flag2 || !window.flag {
                (window.pbmin, window.pbmax, prbmin1, prbmax1)
            } else {
                (window.prmin, window.prmax, prmin1, prmax1)
            };
            // `cxx:1157`.
            if min <= pmax - delta && max >= pmin + delta {
                keep_parts = true;
                the_builder.keep_part(part);
            }
        }
        // `cxx:1163-1175`.
        if !keep_parts {
            form.my_status_error = StatusError::NoParts;
            form.not_done();
            return Ok(false);
        }
    } else {
        // `cxx:1178-1211`: every overlapping part is composed.
        let mut compo: Vec<TopoShape> = Vec::new();
        for part in lshape {
            let (prmin1, prmax1, prbmin1, prbmax1, flag2) = part_window(part, &curve);
            let (pmin, pmax, min, max) = if !flag2 || !window.flag {
                (window.pbmin, window.pbmax, prbmin1, prbmax1)
            } else {
                (window.prmin, window.prmax, prmin1, prmax1)
            };
            // `cxx:1203`.
            if min < pmax - delta && max > pmin + delta {
                compo.push(part.clone());
            }
        }
        form.my_shape = Some(make_compound_shapes(&compo).0);
    }

    Ok(true)
}

/// The interval of one part of the tool (`cxx:1123-1141`), with the periodic
/// arm folding both ends into one period.
fn part_window(part: &TopoShape, curve: &Arc<dyn Curve>) -> (f64, f64, f64, f64, bool) {
    if curve.is_periodic() {
        let period = curve.period();
        let (pr, prmax1, prb, prbmax1, flag2) = parametric_min_max(part, curve, true);
        let prmin1 = if flag2 {
            in_period(pr, prmax1 - period, prmax1)
        } else {
            pr
        };
        let prbmin1 = in_period(prb, prbmax1 - period, prbmax1);
        (prmin1, prmax1, prbmin1, prbmax1, flag2)
    } else {
        parametric_min_max(part, curve, false)
    }
}

/// `cxx:944-1027`: the `mySUntil` narrowing.
fn narrow_until(
    window: &mut Window,
    sens: i32,
    delta: f64,
    face: (f64, f64, f64, f64, bool),
    part: (f64, f64, f64, f64, bool),
) {
    let (prmin1, prmax1, prbmin1, prbmax1, flag3) = face;
    let (prmin2, prmax2, _, _, flag2) = part;

    if sens == 1 {
        // `cxx:970-994`.
        let mut test_ok = !flag2;
        if flag2 {
            test_ok = !window.flag;
            if window.flag && prmax2 > window.prmin + delta {
                test_ok = !flag3;
                if flag3 && prmax1 == prmax2 {
                    test_ok = true;
                }
            }
        }
        if prbmin1 < window.pbmax && test_ok {
            if flag2 {
                window.flag = flag2;
                window.prmax = prmax2;
            }
            window.pbmax = prbmin1;
        }
    } else if sens == -1 {
        // `cxx:996-1015`.
        let mut test_ok = !flag2;
        if flag2 {
            test_ok = !window.flag;
            if window.flag && prmin2 < window.prmax - delta {
                test_ok = !flag3;
                if flag3 && prmin1 == prmin2 {
                    test_ok = true;
                }
            }
        }
        if prbmax1 > window.pbmin && test_ok {
            if flag2 {
                window.flag = flag2;
                window.prmin = prmin2;
            }
            window.pbmin = prbmax1;
        }
    }
}

/// `cxx:1028-1111`: the `mySFrom` narrowing.
fn narrow_from(
    window: &mut Window,
    sens: i32,
    delta: f64,
    face: (f64, f64, f64, f64, bool),
    part: (f64, f64, f64, f64, bool),
) {
    let (prmin1, prmax1, prbmin1, prbmax1, flag3) = face;
    let (prmin2, prmax2, _, _, flag2) = part;

    if sens == 1 {
        // `cxx:1054-1083`.
        let mut test_ok = !flag2;
        if flag2 {
            test_ok = !window.flag;
            if window.flag && prmin2 < window.prmax - delta {
                test_ok = !flag3;
                if flag3 && prmin1 == prmin2 {
                    test_ok = true;
                }
            }
        }
        if prbmax1 > window.pbmin && test_ok {
            if flag2 {
                window.flag = flag2;
                window.prmin = prmin2;
            }
            window.pbmin = prbmax1;
        }
    } else if sens == -1 {
        // `cxx:1082-1107`.
        let mut test_ok = !flag2;
        if flag2 {
            test_ok = !window.flag;
            if window.flag && prmax2 > window.prmin + delta {
                test_ok = !flag3;
                if flag3 && prmax1 == prmax2 {
                    test_ok = true;
                }
            }
        }
        if prbmin1 < window.pbmax && test_ok {
            if flag2 {
                window.flag = flag2;
                window.prmax = prmax2;
            }
            window.pbmax = prbmin1;
        }
    }
}
