//! Curve classification and planar/normal analysis for the curve dispatch. Source: `GeomToIGES_GeomCurve::TransferCurve` (the `IsKind` chain) and `IGESConvGeom_GeomBuilder`.
use super::prelude::*;

// ---------------------------------------------------------------------------
// Curve classification
// ---------------------------------------------------------------------------

/// Curve kind for the IGES curve dispatch.
///
/// `GeomToIGES_GeomCurve::TransferCurve` (`GeomToIGES_GeomCurve.cxx:75-116`)
/// dispatches on the curve's **exact type** (`IsKind`), never on sampled
/// geometry: `Geom_BoundedCurve` (Bezier / B-spline / trimmed) → the 126/100
/// family, `Geom_Conic` → 104 (circle 100), `Geom_Line` → 110. A
/// `Geom_OffsetCurve` is not classified here: with the default
/// `write.iges.offset.mode = 0` it is approximated to a B-spline first
/// ([`IgesWriter::emit_offset_curve_default`]). The port has the equivalent
/// type queries on [`Curve`], so the previous
/// "6 samples of |C''| with a 2% spread" classifier (audit A26) is replaced by
/// this dispatch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum IgCurveKind {
    Line,
    Circle,
    /// Ellipse / hyperbola / parabola: OCCT emits conic-arc entity 104
    /// (`GeomToIGES_GeomCurve::TransferCurve(Geom_Conic)`,
    /// `GeomToIGES_GeomCurve.cxx:483-528` -> `:608-700`, `:707-773`,
    /// `:780-845`), see [`IgesWriter::emit_conic_arc`].
    Conic,
    /// Bezier / B-spline: OCCT emits entity 126
    /// (`GeomToIGES_GeomCurve::TransferBSplineCurve`), see
    /// [`IgesWriter::emit_bspline_curve`].
    Bounded,
    Other,
}

pub(super) fn iges_curve_kind(c: &dyn Curve) -> IgCurveKind {
    if c.is_line() {
        IgCurveKind::Line
    } else if c.gp_circ().is_some() {
        IgCurveKind::Circle
    } else if c.gp_ellipse().is_some() || c.gp_hyperbola().is_some() || c.gp_parabola().is_some() {
        IgCurveKind::Conic
    } else if c.bspline_poles().is_some() || c.bezier_poles().is_some() {
        IgCurveKind::Bounded
    } else {
        IgCurveKind::Other
    }
}

/// The curve `TransferCurve` approximates when `write.iges.offset.mode` is 0,
/// plus the range that approximation is then segmented to.
///
/// A trimmed curve is forwarded to its basis (`GeomToIGES_GeomCurve.cxx:467-475`).
/// `GeomTrimmedCurve` reports `FirstParameter() == 0` and `LastParameter() == 1`
/// while `untrimmed_basis` carries the basis interval, so `(a, b)` in that
/// unit interval is mapped back. `GeomTrimmedCurveBasis` already evaluates in
/// basis parameters (`first_parameter` is the trim start), and a bare
/// `Geom_OffsetCurve` is approximated as received.
pub(super) fn offset_curve_request(
    curve: &dyn Curve,
    a: f64,
    b: f64,
) -> Option<(Arc<dyn Curve>, f64, f64)> {
    if curve.is_geom_trimmed() {
        let (basis, t0, t1) = curve.untrimmed_basis()?;
        if basis.offset_curve().is_none() {
            return None;
        }
        let unit = curve.first_parameter().abs() <= f64::EPSILON
            && (curve.last_parameter() - 1.0).abs() <= f64::EPSILON;
        let basis_is_unit = t0.abs() <= 1e-12 && (t1 - 1.0).abs() <= 1e-12;
        let (u0, u1) = if unit && !basis_is_unit {
            let span = t1 - t0;
            (t0 + a * span, t0 + b * span)
        } else {
            (a, b)
        };
        return Some((basis, u0, u1));
    }
    if curve.offset_curve().is_some() {
        return Some((Arc::from(curve.clone_dyn()), a, b));
    }
    None
}

/// `mycurve->Copy()` for the `Geom_BSplineCurve` transfer
/// (`GeomToIGES_GeomCurve.cxx:301`, `:337`): the port cannot downcast
/// `&dyn Curve`, so the copy is rebuilt from the same
/// `bspline_poles`/`bspline_weights`/`bspline_knots`/`nurbs_degree`/
/// `is_periodic` queries (equal data, no re-approximation). `None` where OCCT's
/// `occ::down_cast` would yield a null handle.
pub(super) fn bspline_copy_of(curve: &dyn Curve) -> Option<occt_geom::bspline_curve::GeomBSplineCurve> {
    let poles = curve.bspline_poles()?.to_vec();
    let knots = curve.bspline_knots()?.to_vec();
    let degree = curve.nurbs_degree()?;
    Some(occt_geom::bspline_curve::GeomBSplineCurve::from_parts(
        poles,
        curve.bspline_weights().map(|w| w.to_vec()),
        knots,
        degree,
        curve.is_periodic(),
    ))
}

/// `ArePolesPlanar` (`GeomToIGES_GeomCurve.cxx:170-199`): the area vector
/// `P(n) x P(1) + sum_{i<n} P(i) x P(i+1)`, normalised; every pole must then sit
/// at the same distance from the plane through `P(1)`. Fewer than three poles is
/// planar by definition and gets an arbitrary normal perpendicular to the first
/// segment (`GetAnyNormal`).
pub(super) fn poles_planar_and_normal(poles: &[GpPnt]) -> (bool, GpVec) {
    let xyz = |p: &GpPnt| p.coord;
    if poles.len() < 3 {
        let d = xyz(&poles[0]).subtracted(&xyz(&poles[1]));
        return (true, any_normal(&d));
    }
    let mut n = xyz(&poles[poles.len() - 1]).crossed(&xyz(&poles[0]));
    for i in 0..poles.len() - 1 {
        n = n.added(&xyz(&poles[i]).crossed(&xyz(&poles[i + 1])));
    }
    let modulus = n.modulus();
    if modulus < occt_core::precision::CONFUSION {
        return (false, GpVec::new(0.0, 0.0, 1.0));
    }
    let normal = GpVec::from_xyz(&n.divided(modulus));
    let scl = xyz(&poles[0]).dot(&normal.xyz());
    for p in &poles[1..] {
        if (xyz(p).dot(&normal.xyz()) - scl).abs() > occt_core::precision::CONFUSION {
            return (false, normal);
        }
    }
    (true, normal)
}

/// `IsPlanar` (`GeomToIGES_GeomCurve.cxx:236-272`): the planar flag and plane
/// normal OCCT computes on the curve handed to `TransferCurve` (`cxx:414`,
/// `IsPlanar(start, Norm)`), i.e. before it is unperiodized or `Segment`-ed into
/// the copy whose poles are written out. The dispatch is OCCT's: `Geom_Line`
/// takes `GetAnyNormal` of its direction, the four conics their `Axis()`
/// direction, a `Geom_TrimmedCurve`/`Geom_OffsetCurve` recurses into the basis,
/// and `Geom_BSplineCurve`/`Geom_BezierCurve` fall through to
/// `ArePolesPlanar`. Anything else returns `false` with an untouched normal.
pub(super) fn curve_planar_and_normal(curve: &dyn Curve) -> (bool, GpVec) {
    if let Some(l) = curve.gp_line() {
        return (true, any_normal(&l.direction().xyz()));
    }
    let conic = curve
        .gp_circ()
        .map(|c| c.pos)
        .or_else(|| curve.gp_ellipse().map(|e| e.pos))
        .or_else(|| curve.gp_hyperbola().map(|h| h.pos))
        .or_else(|| curve.gp_parabola().map(|p| p.pos));
    if let Some(pos) = conic {
        return (true, GpVec::from_xyz(&pos.direction().xyz()));
    }
    if let Some((base, _, _)) = curve.untrimmed_basis() {
        return curve_planar_and_normal(base.as_ref());
    }
    if let Some((base, _)) = curve.offset_curve() {
        return curve_planar_and_normal(base.as_ref());
    }
    if let Some(p) = curve.bspline_poles() {
        return poles_planar_and_normal(p);
    }
    if let Some(p) = curve.bezier_poles() {
        return poles_planar_and_normal(p);
    }
    (false, GpVec::new(0.0, 0.0, 0.0))
}

/// `GetAnyNormal` (`GeomToIGES_GeomCurve.cxx:168-189`): `(0,0,1)` when `|Z|` is
/// below `Precision::Confusion()`, otherwise the unit vector `(Z, 0, -X)`.
pub(super) fn any_normal(v: &occt_core::gp::GpXyz) -> GpVec {
    if v.z().abs() < occt_core::precision::CONFUSION {
        return GpVec::new(0.0, 0.0, 1.0);
    }
    let n = occt_core::gp::GpXyz::new(v.z(), 0.0, -v.x());
    let m = n.modulus();
    if m < occt_core::precision::CONFUSION {
        GpVec::new(0.0, 0.0, 1.0)
    } else {
        GpVec::from_xyz(&n.divided(m))
    }
}

/// `Rational` (`Geom_BSplineCurve.cxx:98-108`): true when two consecutive
/// weights differ by more than `gp::Resolution()`.
pub(super) fn weights_are_rational(weights: &[f64]) -> bool {
    weights.windows(2).any(|pair| {
        (pair[0] - pair[1]).abs() > occt_core::precision::REAL_SMALL
    })
}

/// `(U1 + U2) / 2` over a surface's bounds, the mid point
/// `GeomToIGES_GeomSurface.cxx:1212-1213` evaluates the offset indicator at.
///
/// `Geom_Surface::Bounds` reports an unbounded direction as
/// `Precision::Infinite()` = `RealLast()`, so OCCT's mid point of a
/// doubly-infinite range is 0 and of a half-infinite range is `RealLast()/2`;
/// the port stores `f64::INFINITY` for those, and `INFINITY + -INFINITY` would be
/// NaN - a NaN parameter poisons the returned normal - so a non-finite range
/// collapses to 0 here. Elementary surfaces have constant derivatives in the
/// unbounded direction, so the indicator is unaffected either way.
pub(super) fn bounds_mid_point(a: f64, b: f64) -> f64 {
    if a.is_finite() && b.is_finite() {
        (a + b) / 2.0
    } else {
        0.0
    }
}

/// `CSLib::Normal(D1U, D1V, theSinTol, theStatus, theNormal)`
/// (`CSLib.cxx:45-80`), the overload reached from `GeomLProp_SLProps::Normal()`
/// through `LProp_SurfaceUtils::ComputeSurfNormal`
/// (`GeomLProp_SurfaceUtils.hxx:348-358`).
///
/// Returns the unit normal when the status is `CSLib_Done`. Every other status
/// (`D1IsNull`, `D1uIsNull`, `D1vIsNull`, `D1uIsParallelD1v`) makes OCCT's
/// `Normal()` throw `LProp_NotDefined`, and is reported here as `None`.
pub(super) fn surface_normal_from_d1(d1u: &GpVec, d1v: &GpVec, sin_tol: f64) -> Option<GpVec> {
    // `gp::Resolution()` is `RealSmall()` (`gp.hxx:60`), not `Precision`.
    const GP_RESOLUTION: f64 = occt_core::precision::REAL_SMALL;
    let mag_u = d1u.square_magnitude();
    let mag_v = d1v.square_magnitude();
    if mag_u <= GP_RESOLUTION || mag_v <= GP_RESOLUTION {
        return None;
    }
    let cross = d1u.crossed(d1v);
    let sin2 = cross.square_magnitude() / (mag_u * mag_v);
    if sin2 < sin_tol * sin_tol {
        return None;
    }
    // `gp_Dir(theD1UxD1V)`: OCCT normalises unconditionally here (the magnitude
    // is above `gp::Resolution()` because `sin2` is strictly positive), so this
    // divides explicitly instead of going through `GpVec::normalized`, which
    // keeps an unnormalised vector below `Precision::Resolution()`.
    let m = cross.magnitude();
    if m <= 0.0 {
        return None;
    }
    Some(cross.divided(m))
}

/// The `Geom_SweptSurface` family of `GeomToIGES_GeomSurface::TransferSurface`
/// (`GeomToIGES_GeomSurface.cxx:1000-1025`): a linear extrusion becomes entity
/// 122, a revolution entity 120.
#[derive(Clone, Copy, PartialEq, Eq)]
pub(super) enum SweptKind {
    Extrusion,
    Revolution,
}

/// Resolve the swept surface behind a face's surface the way
/// `GeomToIGES_GeomSurface::TransferSurface` dispatches: a
/// `Geom_RectangularTrimmedSurface` is peeled first (`cxx:492-515`, the
/// `Bounded` branch at `:177-181` recurses on `BasisSurface`), then the exact
/// `Geom_SweptSurface` type decides (`cxx:1013-1022`, extrusion before
/// revolution). Anything else - including a `Geom_OffsetSurface`, which is not a
/// `Geom_SweptSurface` - has no swept branch.
pub(super) fn swept_surface_kind(s: &dyn Surface) -> Option<SweptKind> {
    if let Some(b) = s.rectangular_trimmed_basis() {
        return swept_surface_kind(b.as_ref());
    }
    if s.is_surface_of_linear_extrusion() {
        Some(SweptKind::Extrusion)
    } else if s.is_surface_of_revolution() {
        Some(SweptKind::Revolution)
    } else {
        None
    }
}

/// `Geom_BSplineCurve::IsEqual` (`Geom_BSplineCurve_1.cxx:662-...`) between the
/// two iso curves that `Geom_BSplineSurface::IsUClosed`/`IsVClosed`
/// (`Geom_BSplineSurface_1.cxx:1350-1389`) compares: the pole lists must agree
/// **component-wise** within `Precision::Confusion()` (not by point distance), and
/// every other quantity `IsEqual` looks at (degree, knot count and values,
/// multiplicity list, rational flag, pole count) is shared by construction,
/// because both iso curves are built from the same surface
/// (`Geom_BSplineSurface::UIso`/`VIso`) and therefore use the same V/U knots,
/// multiplicities, degree and periodicity. Differing lengths fail the comparison,
/// as `IsEqual` does on a count mismatch.
pub(super) fn iso_rows_equal(a: &[GpPnt], b: &[GpPnt]) -> bool {
    let tol = occt_core::precision::CONFUSION;
    a.len() == b.len()
        && a.iter().zip(b).all(|(p, q)| {
            (p.x() - q.x()).abs() <= tol && (p.y() - q.y()).abs() <= tol && (p.z() - q.z()).abs() <= tol
        })
}

/// The rational arm of `Geom_BSplineCurve::IsEqual`
/// (`Geom_BSplineCurve_1.cxx`, `fabs(w1 - w2) > Epsilon(w1)`): weights are
/// compared against `Epsilon` of the **first** curve's weight
/// (`Standard_Real.hxx:242-246`).
pub(super) fn weight_equal(a: f64, b: f64) -> bool {
    (a - b).abs() <= occt_core::precision::epsilon(a)
}

/// The unit of the model space the Global section declares
/// (`IGESData_IGESWriter.cxx:38-45` writes `2Hmm` with unit flag 2 = millimetre,
/// so `GeomToIGES_GeomEntity::GetUnit()` is 1: distances are written as they are).
pub(super) const IGES_UNIT: f64 = 1.0;

/// `IGESConvGeom_GeomBuilder::IsIdentity` (`IGESConvGeom_GeomBuilder.cxx:159-175`,
/// `epsl = epsa = 1.E-10` at `:29-30`): the frame's 3x3 part is the identity
/// within `epsa` and its translation part vanishes within `epsl`. In OCCT `thepos`
/// is the frame's **local → global** transformation (`SetPosition(pos)` stores
/// `Trsf(pos, gp::XOY())`, `:137-143`) and `EvalXYZ` applies its inverse
/// (`:212-216`), so `thepos` is the identity exactly when the frame is the
/// absolute one.
pub(super) fn frame_is_identity(frame: &GpAx3) -> bool {
    let eps = 1e-10;
    let (x, y, z) = (
        *frame.x_direction().xyz(),
        *frame.y_direction().xyz(),
        *frame.direction().xyz(),
    );
    let o = frame.location().coord;
    let m = [
        [x.x(), y.x(), z.x()],
        [x.y(), y.y(), z.y()],
        [x.z(), y.z(), z.z()],
    ];
    for (i, row) in m.iter().enumerate() {
        for (j, v) in row.iter().enumerate() {
            let cons = if i == j { 1.0 } else { 0.0 };
            if (*v - cons).abs() > eps {
                return false;
            }
        }
    }
    o.x().abs() <= eps && o.y().abs() <= eps && o.z().abs() <= eps
}

/// `IGESConvGeom_GeomBuilder::EvalXYZ` (`IGESConvGeom_GeomBuilder.cxx:212-216`):
/// the point expressed in the frame's local coordinates (the frame's axes are
/// orthonormal, so this is `thepos.Inverted()` applied to the point).
pub(super) fn frame_local_point(frame: &GpAx3, p: &GpPnt) -> (f64, f64, f64) {
    let d = p.coord.subtracted(&frame.location().coord);
    (
        d.dot(frame.x_direction().xyz()),
        d.dot(frame.y_direction().xyz()),
        d.dot(&frame.direction().xyz()),
    )
}

/// `gp_Elips2d::Coefficients` (`gp_Elips2d.cxx:25-62`) evaluated on the frame
/// `gp_Ax22d(gp::Origin2d(), gp::DX2d(), gp::DY2d())` that
/// `GeomToIGES_GeomCurve.cxx:671` builds. With that frame the `gp_Trsf2d T` there
/// is the identity, so the general sums collapse to
/// `A = 1/DMaj, B = 1/DMin, C = D = E = 0, F = -1`; the degenerate arm
/// (`DMin <= gp::Resolution()`, `:35-45`) keeps its own form. The result is in
/// `gp_Elips2d::Coefficients`' **own** output order — the caller in
/// `GeomToIGES_GeomCurve.cxx:675` passes it in the shuffled order `(A, C, B, D,
/// E, F)`.
pub(super) fn gp_elips2d_coefficients(major: f64, minor: f64) -> (f64, f64, f64, f64, f64, f64) {
    let dmin = minor * minor;
    let dmaj = major * major;
    if dmin <= occt_core::precision::REAL_SMALL && dmaj <= occt_core::precision::REAL_SMALL {
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    } else if dmin <= occt_core::precision::REAL_SMALL {
        (1.0, 0.0, 0.0, 0.0, 0.0, -dmaj)
    } else {
        (1.0 / dmaj, 1.0 / dmin, 0.0, 0.0, 0.0, -1.0)
    }
}

/// `gp_Hypr2d::Coefficients` (`gp_Hypr2d.cxx:25-61`) on the same identity frame:
/// `A = 1/DMaj, B = -1/DMin, C = D = E = 0, F = -1`, with the same degenerate arm.
pub(super) fn gp_hypr2d_coefficients(major: f64, minor: f64) -> (f64, f64, f64, f64, f64, f64) {
    let dmin = minor * minor;
    let dmaj = major * major;
    if dmin <= occt_core::precision::REAL_SMALL && dmaj <= occt_core::precision::REAL_SMALL {
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    } else if dmin <= occt_core::precision::REAL_SMALL {
        (1.0, 0.0, 0.0, 0.0, 0.0, -dmaj)
    } else {
        (1.0 / dmaj, -1.0 / dmin, 0.0, 0.0, 0.0, -1.0)
    }
}

/// `gp_Parab2d::Coefficients` (`gp_Parab2d.cxx:44-60`) on the same identity frame,
/// with `P = 2 * focalLength` and `focalLength = start->Focal() * 2`
/// (`GeomToIGES_GeomCurve.cxx:818-819`): `A = C = 0`, `B = 1`, `D = -P`,
/// `E = F = 0`.
pub(super) fn gp_parab2d_coefficients(focal: f64) -> (f64, f64, f64, f64, f64, f64) {
    let p = 2.0 * (focal * 2.0);
    (0.0, 1.0, 0.0, -p, 0.0, 0.0)
}

/// `GeomToIGES_GeomCurve.cxx:722-729` (hyperbola) and `:795-802` (parabola): a
/// range end past `Precision::Infinite()` becomes ±`Precision::Infinite()`. The
/// ellipse branch has no such arm.
pub(super) fn occt_infinite_range(a: f64, b: f64) -> (f64, f64) {
    let u1 = if occt_core::precision::Precision::is_negative_infinite(a) {
        -occt_core::precision::INFINITE
    } else {
        a
    };
    let u2 = if occt_core::precision::Precision::is_positive_infinite(b) {
        occt_core::precision::INFINITE
    } else {
        b
    };
    (u1, u2)
}

/// `IGESGeom_ConicArc::ComputedFormNumber` (`IGESGeom_ConicArc.cxx:97-124`):
/// form 1 ellipse, 2 hyperbola, 3 parabola, 0 when the coefficients identify no
/// conic. `eps = 1.E-08` there, `eps2 = eps*eps`, `eps4 = eps2*eps2`.
pub(super) fn conic_form_number(a: f64, b: f64, c: f64, d: f64, e: f64, f: f64) -> i32 {
    let eps = 1e-8f64;
    let eps4 = (eps * eps) * (eps * eps);
    let q1 = a * (c * f - e * e / 4.0) + b / 2.0 * (e * d / 4.0 - b * f / 2.0)
        + d / 2.0 * (b * e / 4.0 - c * d / 2.0);
    let q2 = a * c - b * b / 4.0;
    let q3 = a + c;
    if q2 > eps4 && q1 * q3 < 0.0 {
        1
    } else if q2 < -eps4 && q1.abs() > eps4 {
        2
    } else if q2.abs() <= eps4 && q1.abs() > eps4 {
        3
    } else {
        0
    }
}

