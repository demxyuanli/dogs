//! Shared 2D-curve sample counts, ported from the OCCT `Geom2d` stack.
//!
//! * [`nb_points`] ports the file-static `nbPoints(const Handle(Geom2d_Curve)&)`
//!   (`Geom2dAdaptor_Curve.cxx:1351-1389`), the body of
//!   `Geom2dAdaptor_Curve::NbSamples()` (`Geom2dAdaptor_Curve.cxx:1391-1394`).
//! * [`nb_samples`] ports
//!   `Geom2dInt_Geom2dCurveTool::NbSamples(const Adaptor2d_Curve2d&)`
//!   (`Geom2dInt_Geom2dCurveTool.cxx:73-91`).
//!
//! The `if (nbs > 2) nbs *= 4;` step is deliberately NOT in either function:
//! OCCT applies it at each call site (`BRepTopAdaptor_FClass2d.cxx:179-185`,
//! `ShapeAnalysis_Curve.cxx:1325-1331`).
//!
//! Moved here from `occt-topo::curve_sampling_2d` so the `Geom2dInt` port in
//! this crate can use it; `occt-topo` now re-exports it.

use crate::curve::Curve2d;

/// `nbPoints(theCurve)` (`Geom2dAdaptor_Curve.cxx:1351-1389`), i.e.
/// `Geom2dAdaptor_Curve::NbSamples()` (`Geom2dAdaptor_Curve.cxx:1391-1394`).
pub fn nb_points(pc: &dyn Curve2d) -> usize {
    // OCCT tests `IsKind(Geom2d_Line)` first (`cxx:1356-1359`), but
    // `Geom2dTrimmedCurve::is_line()` delegates to its basis, so the trimmed
    // arm must be tested first here; an OCCT `Geom2dAdaptor_Curve` never sees a
    // `Geom2d_TrimmedCurve` (`Geom2dAdaptor_Curve.cxx:285-288`).
    if let Some(basis) = pc.trimmed_basis() {
        return 20.max(nb_points(basis)); // cxx:1379-1383
    }
    if pc.is_line() {
        return 2; // cxx:1356-1359
    }
    if let Some(nb_poles) = pc.bezier_nb_poles() {
        return (3 + nb_poles).min(300); // cxx:1360-1363 + clamp 1384-1387
    }
    if let (Some(nb_knots), Some(degree)) = (pc.bspline_nb_knots(), pc.bspline_degree()) {
        return (nb_knots * degree).max(2).min(300); // cxx:1364-1372 + clamp
    }
    if let Some(basis) = pc.offset_basis() {
        return 20.max(nb_points(basis)); // cxx:1373-1377
    }
    20.min(300) // default nbs = 20 (cxx:1354) under the cxx:1384-1387 clamp
}

/// `Geom2dInt_Geom2dCurveTool::NbSamples(const Adaptor2d_Curve2d&)`
/// (`Geom2dInt_Geom2dCurveTool.cxx:73-91`).
pub fn nb_samples(pc: &dyn Curve2d, first: f64, last: f64) -> usize {
    let mut nbs = nb_points(pc);
    if let Some(circ) = pc.gp_circ2d() {
        if circ.radius() > 1.0 {
            let angl = 0.283079; // 2. * acos(1. - eps), eps = 0.01
            let n = ((last - first) / angl) as usize; // RealToInt truncates
            nbs = n.max(nbs);
        }
    }
    nbs
}

/// `Geom2dInt_Geom2dCurveTool::NbSamples(C, U0, U1)`
/// (`Geom2dInt_Geom2dCurveTool.cxx:23-70`).
pub fn nb_samples_range(pc: &dyn Curve2d, u0: f64, u1: f64) -> usize {
    let mut nbs = nb_points(pc);
    if let (Some(nb_knots), Some(degree)) = (pc.bspline_nb_knots(), pc.bspline_degree()) {
        let t = pc.last_parameter() - pc.first_parameter();
        if t > occt_core::precision::PCONFUSION {
            let mut t1 = u1 - u0;
            if t1 < 0.0 {
                t1 = -t1;
            }
            let anb = t1 / t * (nb_knots * degree) as f64;
            nbs = anb as usize;
            let min_pnt_nb = (degree + 1).max(4);
            if nbs < min_pnt_nb {
                nbs = min_pnt_nb;
            }
        }
    } else if let Some(circ) = pc.gp_circ2d() {
        if circ.radius() > 1.0 {
            let angl = 0.283079;
            let n = ((u1 - u0).abs() / angl) as usize;
            nbs = n.max(nbs);
        }
    }
    if nbs > 300 {
        nbs = 300;
    }
    nbs
}
