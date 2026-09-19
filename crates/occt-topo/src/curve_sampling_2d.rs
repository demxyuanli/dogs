//! Shared 2D-curve sample counts, ported from the OCCT `Geom2d` stack.
//!
//! Two OCCT functions live here, because two call sites in this crate need the
//! same numbers:
//!
//! * [`nb_points`] ports the file-static `nbPoints(const Handle(Geom2d_Curve)&)`
//!   (`Geom2dAdaptor_Curve.cxx:1351-1389`). That function is the body of
//!   `Geom2dAdaptor_Curve::NbSamples()` (`Geom2dAdaptor_Curve.cxx:1391-1394`).
//! * [`nb_samples`] ports
//!   `Geom2dInt_Geom2dCurveTool::NbSamples(const Adaptor2d_Curve2d&)`
//!   (`Geom2dInt_Geom2dCurveTool.cxx:73-91`).
//!
//! The `if (nbs > 2) nbs *= 4;` step is deliberately NOT in either function:
//! OCCT applies it at each call site, and at a different place in each caller:
//!
//! * `BRepTopAdaptor_FClass2d.cxx:179-185` (the face-classifier boundary ring)
//! * `ShapeAnalysis_Curve.cxx:1325-1331` (`GetSamplePoints` for a 2D curve)
//!
//! so each caller here re-applies it itself.
//!
//! All line numbers above refer to OCCT 8.0.0, the version this port targets and
//! the version the reference `DRAWEXE` and the reference `occ-*.obj` files come
//! from (`adm/cmake/version.cmake`: 8.0.0). A second, unrelated 7.9.3 tree exists
//! under the local vcpkg buildtrees; its line numbers differ and must NOT be
//! used as the authority for this project.

use occt_geom2d::curve::Curve2d;

/// `nbPoints(theCurve)` (`Geom2dAdaptor_Curve.cxx:1351-1389`), i.e.
/// `Geom2dAdaptor_Curve::NbSamples()` (`Geom2dAdaptor_Curve.cxx:1391-1394`).
pub(crate) fn nb_points(pc: &dyn Curve2d) -> usize {
    // Branch order note: OCCT tests `IsKind(Geom2d_Line)` first
    // (`Geom2dAdaptor_Curve.cxx:1356-1359`), but our `Geom2dTrimmedCurve::is_line()`
    // delegates to its basis, so the trimmed arm must be tested first here. An
    // OCCT `Geom2dAdaptor_Curve` never sees a `Geom2d_TrimmedCurve`: `load`
    // unwraps it onto its basis (`Geom2dAdaptor_Curve.cxx:285-288`). This only
    // corrects the delegation.
    if let Some(basis) = pc.trimmed_basis() {
        // `Geom2dAdaptor_Curve.cxx:1379-1383`: `max(nbs, nbPoints(basis))`,
        // returned without the 300 clamp in this arm (the recursive call
        // applies its own clamp).
        return 20.max(nb_points(basis));
    }
    if pc.is_line() {
        // `Geom2dAdaptor_Curve.cxx:1356-1359`.
        return 2;
    }
    if let Some(nb_poles) = pc.bezier_nb_poles() {
        // `Geom2dAdaptor_Curve.cxx:1360-1363`: `nbs = 3 + NbPoles`.
        return (3 + nb_poles).min(300); // clamp `cxx:1384-1387`
    }
    if let (Some(nb_knots), Some(degree)) = (pc.bspline_nb_knots(), pc.bspline_degree()) {
        // `Geom2dAdaptor_Curve.cxx:1364-1372`: `nbs = NbKnots * Degree`,
        // raised to 2 when below.
        return (nb_knots * degree).max(2).min(300); // clamp `cxx:1384-1387`
    }
    if let Some(basis) = pc.offset_basis() {
        // `Geom2dAdaptor_Curve.cxx:1373-1377`: `max(nbs, nbPoints(basis))`,
        // returned without the 300 clamp in this arm.
        return 20.max(nb_points(basis));
    }
    // Default `nbs = 20` (`Geom2dAdaptor_Curve.cxx:1354`) under the
    // `cxx:1384-1387` clamp.
    20.min(300)
}

/// `Geom2dInt_Geom2dCurveTool::NbSamples(const Adaptor2d_Curve2d&)`
/// (`Geom2dInt_Geom2dCurveTool.cxx:73-91`). `first`/`last` are the adaptor's
/// `FirstParameter()`/`LastParameter()`, i.e. the edge pcurve range on the face.
///
/// Note there is no `nbs *= 4` and no `> 300` clamp in this overload; both the
/// `*4` (`BRepTopAdaptor_FClass2d.cxx:182-184`) and the meshing site's own
/// steps stay with the caller. The `> 300` clamp does exist in the 3-argument
/// overload `Geom2dInt_Geom2dCurveTool.cxx:23-70` that meshing uses, which is
/// already applied inside [`nb_points`].
pub(crate) fn nb_samples(pc: &dyn Curve2d, first: f64, last: f64) -> usize {
    let mut nbs = nb_points(pc);
    if let Some(circ) = pc.gp_circ2d() {
        if circ.radius() > 1.0 {
            // `Geom2dInt_Geom2dCurveTool.cxx:79-87`: try to reach
            // deflection = eps * R with eps = 0.01. `RealToInt` truncates
            // toward zero; a negative range saturates to 0 and then loses the
            // `Max`, exactly as the OCCT `Max` with a negative `n` would.
            let angl = 0.283079; // 2. * acos(1. - eps)
            let n = ((last - first) / angl) as usize;
            nbs = n.max(nbs);
        }
    }
    nbs
}
