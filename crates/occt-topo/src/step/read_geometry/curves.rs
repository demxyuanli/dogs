//! Edge construction from a resolved 3D curve.

use super::super::prelude::*;
use super::*;

/// `StepToTopoDS_TranslateEdge::MakeFromCurve3D` (`cxx:442-449`):
/// `ShapeAnalysis_Curve::Project` + `UpdateParam3d`. Returns the (possibly
/// reversed) curve together with the projected `[U1, U2]` range.
pub(in crate::step) fn edge_from_curve3d(
    curve: Arc<dyn Curve>,
    p1: &GpPnt,
    p2: &GpPnt,
) -> (Arc<dyn Curve>, f64, f64) {
    const PRECI: f64 = 1e-3;
    let (first, last) = edge_params_for_curve(curve.as_ref(), p1, p2);
    let mut w1 = first;
    let mut w2 = last;
    let curve = update_param3d(curve, &mut w1, &mut w2, PRECI);
    let curve = shift_displaced_line(curve, p1, p2, w1, w2, PRECI);
    (curve, w1, w2)
}

/// `MakeFromCurve3D` (`TranslateEdge.cxx:451-472`): a line that keeps a
/// constant offset from both vertices and whose parameter span equals the
/// vertex distance is translated onto the vertices.
fn shift_displaced_line(
    curve: Arc<dyn Curve>,
    p1: &GpPnt,
    p2: &GpPnt,
    w1: f64,
    w2: f64,
    preci: f64,
) -> Arc<dyn Curve> {
    if !curve.is_line() {
        return curve;
    }
    let temp1 = curve.d0(w1).distance(p1);
    let temp2 = curve.d0(w2).distance(p2);
    if temp1 <= preci && temp2 <= preci {
        return curve;
    }
    if (temp1 - temp2).abs() >= preci {
        return curve;
    }
    if ((w2 - w1).abs() - p1.distance(p2)).abs() >= occt_core::precision::CONFUSION {
        return curve;
    }
    let d1 = curve.d1(0.0).1;
    let Ok(dir) = GpDir::from_vec(&d1) else {
        return curve;
    };
    let origin = p1.translated_vec(&GpVec::from_xyz(&dir.xyz().multiplied(-w1)));
    Arc::new(GeomLine::from_pnt_dir(origin, dir))
}

/// `StepToTopoDS_TranslateEdge::MakeFromCurve3D` (`cxx:442-444`): the edge range
/// is `ShapeAnalysis_Curve::Project` of both endpoint points with
/// `AdjustToEnds = false` — for **every** curve type, with no analytic-type
/// dispatch for the range itself and no endpoint-distance shortcut. `Project`
/// keeps its own endpoint early-exit for `Geom_BoundedCurve` (`cxx:161-182`) and
/// has no failure channel (`ShapeAnalysis_Curve.cxx:147-201` always produces a
/// parameter); `UpdateParam3d` (cxx:446) and the displaced-Line shift
/// (cxx:451-472) are applied by [`edge_from_curve3d`].
pub(in crate::step) fn edge_params_for_curve(curve: &dyn Curve, p1: &GpPnt, p2: &GpPnt) -> (f64, f64) {
    const PRECI: f64 = 1e-3;
    let (u1, _) = shape_analysis_project(curve, p1, PRECI);
    let (u2, _) = shape_analysis_project(curve, p2, PRECI);
    (u1, u2)
}

/// `StepToTopoDS_GeometricTool::UpdateParam3d` (`GeometricTool.cxx:227-409`).
fn update_param3d(
    curve: Arc<dyn Curve>,
    w1: &mut f64,
    w2: &mut f64,
    preci: f64,
) -> Arc<dyn Curve> {
    let cf = curve.first_parameter();
    let cl = curve.last_parameter();
    let start = curve.d0(cf);
    let end = curve.d0(cl);
    let is_closed = start.distance(&end) <= preci;

    // Bounded and not closed: clamp projections onto the natural domain
    // (`cxx:238-268`).
    if cf.is_finite() && cl.is_finite() && !curve.is_periodic() && !is_closed {
        if *w1 < cf {
            *w1 = cf;
        } else if *w1 > cl {
            *w1 = cl;
        }
        if *w2 < cf {
            *w2 = cf;
        } else if *w2 > cl {
            *w2 = cl;
        }
    }

    if *w1 < *w2 {
        return curve;
    }

    if curve.is_periodic() {
        crate::geom_bnd_lib_elclib2d::adjust_periodic(
            cf,
            cl,
            occt_core::precision::PCONFUSION,
            w1,
            w2,
        );
        return curve;
    }

    if is_closed {
        if (*w2 - cf).abs() < occt_core::precision::PCONFUSION {
            *w2 = cl;
        } else if (*w1 - cl).abs() < occt_core::precision::PCONFUSION {
            *w1 = cf;
        } else {
            if curve.d0(*w1).distance(&start) < preci {
                *w1 = cf;
            }
            if curve.d0(*w2).distance(&end) < preci {
                *w2 = cl;
            }
            if (*w2 - *w1).abs() < occt_core::precision::PCONFUSION {
                *w1 = cf;
                *w2 = cl;
            } else if *w1 > *w2 {
                std::mem::swap(w1, w2);
            }
        }
        return curve;
    }

    if curve.bspline_poles().is_some() && *w1 > *w2 {
        // Non-closed BSpline with reversed params: Reverse the per-edge copy
        // (`cxx:364-372`). The STEP curve cache is shared; do not mutate it.
        let rp = |w: f64| cf + cl - w;
        *w1 = rp(*w1);
        *w2 = rp(*w2);
        let rev: Arc<dyn Curve> = Arc::from(curve.reversed());
        if (*w1 - *w2).abs() < occt_core::precision::PCONFUSION {
            *w1 = rev.first_parameter();
            *w2 = rev.last_parameter();
        }
        return rev;
    }

    if *w1 > *w2 {
        std::mem::swap(w1, w2);
    }
    curve
}
