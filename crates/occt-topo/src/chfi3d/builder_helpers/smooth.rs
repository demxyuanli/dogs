//! `ChFi3d_IsSmooth` — a fillet guide is "smooth" when its curvature centre
//! does not swing by more than 60 degrees and its radius does not change by
//! more than a factor 2 between two samples.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:5855-5945`.
//!
//! UNPORTED: `ChFi3d_ApproxByC2` (`cxx:5835-5851`). It needs
//! `GeomAPI_PointsToBSpline(points, Approx_ChordLength, 3, 8, GeomAbs_C2,
//! 1.000001e-3)`; the port has global cubic interpolation
//! (`occt_core::geom::curve_interp3d`) but no chord-length approximation with a
//! degree range, so the branch is left out rather than faked.

use occt_core::gp::GpVec;
use occt_core::precision;
use occt_geom::geom_lprop::cl_props::ClProps;
use occt_geom::Curve;

use crate::brep_lprop::CN;

/// `ChFi3d_IsSmooth(C)` (`cxx:5855-5945`).
///
/// The curve is walked interval by interval at 30 samples; a sample with a
/// non-tangent point answers `false` (`cxx:5879-5881`), and from the first
/// sample with a curvature above `gp::Resolution()` onwards the vector from the
/// point to its curvature centre must keep its direction within `pi/3` and its
/// magnitude within a factor 2 (`cxx:5922-5938`). A curve with no curvature
/// anywhere is smooth (`cxx:5900-5903`).
///
/// The `Ratio > 2` test skips the very last sample of the very last interval
/// (`cxx:5934`), so the loops are driven by the interval index rather than by a
/// flat sample index.
///
/// `Err` covers the `LProp_NotDefined` throws OCCT lets escape from
/// `Curvature` / `CentreOfCurvature`; the guards make them unreachable for the
/// samples this function actually takes.
pub fn is_smooth(c: &dyn Curve) -> Result<bool, String> {
    let nb_intv = c.nb_intervals(CN);
    if nb_intv <= 0 {
        return Ok(true);
    }
    let ti = c.parameter_intervals(CN);
    // `NCollection_Array1<double> TI(1, nbIntv + 1);`: the adaptor always hands
    // back `NbIntervals + 1` values, so indexing `intv` and `intv + 1` is in
    // range; the accessor keeps a malformed adaptor from panicking.
    let knot = |i: i32| -> f64 {
        match usize::try_from(i).ok().and_then(|k| ti.get(k)) {
            Some(v) => *v,
            None => c.last_parameter(),
        }
    };
    let resolution = precision::RESOLUTION;
    let discretisation = 30;

    let mut lprop = ClProps::new(c, c.first_parameter(), 2, resolution)?;

    let mut prev_vec = GpVec::zero();
    let mut prev_vec_found = false;
    let mut intv_found = 0i32;
    for intv in 1..=nb_intv {
        let mut t = knot(intv);
        let step = (knot(intv + 1) - t) / discretisation as f64;
        for _ in 1..=discretisation {
            lprop.set_parameter(t);
            if !lprop.is_tangent_defined() {
                return Ok(false);
            }
            let curvature = lprop.curvature()?.abs();
            if curvature > resolution {
                let p1 = lprop.value();
                let p2 = lprop.centre_of_curvature()?;
                prev_vec = GpVec::from_pnts(&p1, &p2);
                prev_vec_found = true;
                break;
            }
            t += step;
        }
        if prev_vec_found {
            intv_found = intv;
            break;
        }
    }

    if !prev_vec_found {
        return Ok(true);
    }

    for intv in intv_found..=nb_intv {
        let mut t = knot(intv);
        let step = (knot(intv + 1) - t) / discretisation as f64;
        for ii in 1..=discretisation {
            lprop.set_parameter(t);
            if !lprop.is_tangent_defined() {
                return Ok(false);
            }
            let curvature = lprop.curvature()?.abs();
            if curvature > resolution {
                let p1 = lprop.value();
                let p2 = lprop.centre_of_curvature()?;
                let vec = GpVec::from_pnts(&p1, &p2);
                let angle = prev_vec.angle(&vec);
                if angle > std::f64::consts::PI / 3.0 {
                    return Ok(false);
                }
                let mut ratio = vec.magnitude() / prev_vec.magnitude();
                if ratio < 1.0 {
                    ratio = 1.0 / ratio;
                }
                if ratio > 2.0 && (intv != nb_intv || ii != discretisation) {
                    return Ok(false);
                }
                prev_vec = vec;
            }
            t += step;
        }
    }
    Ok(true)
}
