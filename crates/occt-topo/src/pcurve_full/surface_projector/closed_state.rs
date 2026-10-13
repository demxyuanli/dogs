//! Closed/periodic surface state: `RestrictBounds` (`ShapeAnalysis_Surface.cxx:64-91`) and the `IsUClosed`/`IsVClosed` family (`cxx:660-863`, `cxx:865-1077`).

use super::*;

/// `RestrictBounds` (`ShapeAnalysis_Surface.cxx:64-91`): clamp an infinite
/// parameter range to a finite window - `[-1000, 1000]` when both ends are
/// infinite, otherwise a width of 2000 from the finite end.
pub(crate) fn restrict_bounds(first: &mut f64, last: &mut f64) {
    let is_first_inf = occt_core::precision::Precision::is_negative_infinite(*first);
    let is_last_inf = occt_core::precision::Precision::is_positive_infinite(*last);
    if is_first_inf || is_last_inf {
        if is_first_inf && is_last_inf {
            *first = -1000.0;
            *last = 1000.0;
        } else if is_first_inf {
            *first = *last - 2000.0;
        } else {
            *last = *first + 2000.0;
        }
    }
}

// ---------------------------------------------------------------------------
// ShapeAnalysis_Surface::IsUClosed / IsVClosed
// ---------------------------------------------------------------------------

/// `Geom_Surface::IsUClosed` union with `IsUPeriodic`.
///
/// `ShapeAnalysis_Surface::IsUClosed` (`ShapeAnalysis_Surface.cxx:674`) and
/// `IsVClosed` (`cxx:881`) start from the concrete surface's own flag. For a
/// B-spline that flag is `Geom_BSplineSurface::IsUClosed`
/// (`Geom_BSplineSurface_1.cxx:1350-1369`): `IsUPeriodic()` first, then the
/// boundary `UIso` curves compared with `IsEqual(.., Confusion)`. The port's
/// [`Surface::is_u_closed`] implements a pole-row approximation of the second
/// term (`crates/occt-geom/src/bspline_surface.rs`) but omits the periodic
/// term, so the periodic flag is added here. For every other surface type
/// `is_u_closed` already defaults to `is_u_periodic`, so this is a no-op.
pub(crate) fn geom_u_closed(s: &dyn Surface) -> bool {
    s.is_u_closed() || s.is_u_periodic()
}

/// `Geom_Surface::IsVClosed` union with `IsVPeriodic`; see [`geom_u_closed`].
pub(crate) fn geom_v_closed(s: &dyn Surface) -> bool {
    s.is_v_closed() || s.is_v_periodic()
}

/// `Geom_BSplineSurface::UMultiplicity(1)` / `(NbUKnots)` / the same for V:
/// the first / last unique-knot multiplicities from the flat knot vector.
pub(crate) fn end_mults(flat: &[f64]) -> (i32, i32) {
    let (_, mults) = occt_geom::bspline_surface::GeomBSplineSurface::unique_knots_mults(flat);
    let first = mults.first().copied().unwrap_or(0);
    let last = mults.last().copied().unwrap_or(0);
    (first, last)
}

/// `ShapeAnalysis_Surface::IsUClosed` (`ShapeAnalysis_Surface.cxx:660-863`).
///
/// Returns `(myUCloseVal, anUmidVal, myUDelt)`. `myUCloseVal` is the maximum
/// 3D distance between the two U-boundary curves (`cxx:853-854`), `anUmidVal`
/// the distance from the first boundary sample to the mid-U sample of the
/// surface (`-1` = OCCT's "not computed", so the symmetry test
/// `myUCloseVal > sqrt(anUmidVal)` is skipped), and `myUDelt` the
/// `Extrema_ExtPS` window extension `ValueOfUV` reads as
/// `du = min(myUDelt, UResolution(preci))` (`cxx:1346`): a twentieth of the U
/// range (`cxx:672`), reduced by half the U resolution at the minimum boundary
/// gap for the sampled arms (`cxx:764`, `:786`, `:818`, `:849`).
///
/// The `AN_ASSERT(myUCloseVal >= 0)` guard of `cxx:853` is the `myUCloseVal < 0`
/// cache test of `cxx:667`; the cache is not reproduced (each call recomputes -
/// the value is preci-independent, so the boolean is unchanged; see the module
/// notes on the stateless port).
pub(crate) fn sa_u_close_state(s: &dyn Surface) -> (f64, f64, f64) {
    use crate::geom_bnd_lib_surface3d::SurfaceKind as TypeKind;
    // `cxx:668-671`: `Bounds` + `RestrictBounds`.
    let (mut uf, mut ul) = s.u_range();
    let (vf, vl) = s.v_range();
    restrict_bounds(&mut uf, &mut ul);
    // `cxx:672`.
    let mut delt = (ul - uf).abs() / 20.0;
    // `cxx:674-681`.
    if geom_u_closed(s) {
        return (0.0, -1.0, 0.0);
    }
    // `cxx:684-689`: a `Geom_RectangularTrimmedSurface` forces the `default`
    // (101-point) arm even though the adaptor already sees the basis type.
    let surftype = if s.rectangular_trimmed_basis().is_some() {
        TypeKind::OtherSurface
    } else {
        crate::geom_bnd_lib_surface3d::surface_kind(s)
    };
    let mut close = f64::MAX; // `RealLast()`
    let mut mid = -1.0;
    match surftype {
        // `cxx:692-694`.
        TypeKind::Plane => {}
        // `cxx:695-717`.
        TypeKind::SurfaceOfExtrusion => {
            if let Some(crv) = s.extrusion_basis_curve() {
                let f = crv.first_parameter();
                let l = crv.last_parameter();
                if !occt_core::precision::Precision::is_infinite(f)
                    && !occt_core::precision::Precision::is_infinite(l)
                {
                    let p1 = crv.d0(f);
                    let p2 = crv.d0(l);
                    close = p1.square_distance(&p2);
                    let pm = crv.d0(0.5 * (f + l));
                    mid = p1.square_distance(&pm);
                }
            }
        }
        // `cxx:719-788`.
        TypeKind::BSplineSurface => {
            if let Some(bs) = s.osculating_bspline() {
                let nbup = bs.poles.len();
                let mut distmin = f64::MAX;
                let (m1, m2) = end_mults(&bs.knots_u);
                if bs.u_periodic {
                    // `cxx:724-728`.
                    close = 0.0;
                    delt = 0.0;
                } else if nbup < 3 {
                    // `cxx:729-731`.
                } else if bs.is_rational()
                    || m1 != bs.deg_u as i32 + 1
                    || m2 != bs.deg_u as i32 + 1
                {
                    // `cxx:732-765`: sample along the V knots when the boundary
                    // is not a simple pole row (rational, or a boundary knot
                    // multiplicity that is not `degree + 1`).
                    let (vknots, _) = occt_geom::bspline_surface::GeomBSplineSurface::unique_knots_mults(&bs.knots_v);
                    let nbvk = vknots.len();
                    let mut v = vknots[0];
                    let p1 = s.value(uf, v);
                    let p2 = s.value(ul, v);
                    close = p1.square_distance(&p2);
                    let pm = s.value(0.5 * (uf + ul), v);
                    mid = p1.square_distance(&pm);
                    distmin = close;
                    for i in 1..nbvk {
                        // `cxx:747-749`: `v = 0.5 * (VKnot(i-1) + VKnot(i))`.
                        v = 0.5 * (vknots[i - 1] + vknots[i]);
                        let p1 = s.value(uf, v);
                        let p2 = s.value(ul, v);
                        let a_dist = p1.square_distance(&p2);
                        if a_dist > close {
                            close = a_dist;
                            let pm = s.value(0.5 * (uf + ul), v);
                            mid = p1.square_distance(&pm);
                        } else {
                            distmin = distmin.min(a_dist);
                        }
                    }
                    // `cxx:763-764`.
                    delt = delt.min(0.5 * occt_geom::approx_same_parameter::u_resolution(s, distmin.sqrt()));
                } else {
                    // `cxx:767-787`: the boundary pole rows.
                    let nbvp = bs.nb_poles_v();
                    close = bs.poles[0][0].square_distance(&bs.poles[nbup - 1][0]);
                    mid = bs.poles[0][0].square_distance(&bs.poles[nbup / 2][0]);
                    distmin = close;
                    for i in 1..nbvp {
                        let a_dist =
                            bs.poles[0][i].square_distance(&bs.poles[nbup - 1][i]);
                        if a_dist > close {
                            close = a_dist;
                            mid = bs.poles[0][i].square_distance(&bs.poles[nbup / 2][i]);
                        } else {
                            distmin = distmin.min(a_dist);
                        }
                    }
                    // `cxx:785-786`.
                    delt = delt.min(0.5 * occt_geom::approx_same_parameter::u_resolution(s, distmin.sqrt()));
                }
            }
        }
        // `cxx:790-820`.
        TypeKind::BezierSurface => {
            if let Some(bz) = s.osculating_bspline() {
                let nbup = bz.poles.len();
                if nbup >= 3 {
                    let nbvp = bz.nb_poles_v();
                    close = bz.poles[0][0].square_distance(&bz.poles[nbup - 1][0]);
                    mid = bz.poles[0][0].square_distance(&bz.poles[nbup / 2][0]);
                    let mut distmin = close;
                    for i in 0..nbvp {
                        let a_dist =
                            bz.poles[0][i].square_distance(&bz.poles[nbup - 1][i]);
                        if a_dist > close {
                            close = a_dist;
                            mid = bz.poles[0][i].square_distance(&bz.poles[nbup / 2][i]);
                        } else {
                            distmin = distmin.min(a_dist);
                        }
                    }
                    // `cxx:817-818`.
                    delt = delt.min(0.5 * occt_geom::approx_same_parameter::u_resolution(s, distmin.sqrt()));
                }
            }
        }
        // `cxx:822-851`: `Geom_RectangularTrimmedSurface` and
        // `Geom_OffsetSurface`, sampled at 101 V parameters.
        _ => {
            let nbpoints = 101i32;
            let p1 = s.value(uf, vf);
            let p2 = s.value(ul, vf);
            close = p1.square_distance(&p2);
            let pm = s.value(0.5 * (uf + ul), vf);
            mid = p1.square_distance(&pm);
            let mut distmin = close;
            for i in 1..nbpoints {
                let vparam = vf + (vl - vf) * i as f64 / (nbpoints - 1) as f64;
                let p1 = s.value(uf, vparam);
                let p2 = s.value(ul, vparam);
                let a_dist = p1.square_distance(&p2);
                if a_dist > close {
                    close = a_dist;
                    let pm = s.value(0.5 * (uf + ul), vparam);
                    mid = p1.square_distance(&pm);
                } else {
                    distmin = distmin.min(a_dist);
                }
            }
            // `cxx:848-849`.
            delt = delt.min(0.5 * occt_geom::approx_same_parameter::u_resolution(s, distmin.sqrt()));
        }
    }
    // `cxx:853-854`.
    close = close.sqrt();
    // `cxx:857-861`: the mid-line symmetry test discards a "closed" verdict
    // when the boundary gap exceeds the distance to the mid-U line.
    if mid > 0.0 && close > mid.sqrt() {
        return (f64::MAX, -1.0, delt);
    }
    (close, mid, delt)
}

/// `ShapeAnalysis_Surface::IsVClosed` (`ShapeAnalysis_Surface.cxx:865-1077`).
/// Same contract as [`sa_u_close_state`] with the V switch; `myVDelt` is the
/// `dv` counterpart of `myUDelt` (`cxx:879`).
pub(crate) fn sa_v_close_state(s: &dyn Surface) -> (f64, f64, f64) {
    use crate::geom_bnd_lib_surface3d::SurfaceKind as TypeKind;
    let (uf, ul) = s.u_range();
    let (mut vf, mut vl) = s.v_range();
    restrict_bounds(&mut vf, &mut vl);
    // `cxx:879`.
    let mut delt = (vl - vf).abs() / 20.0;
    // `cxx:881-888`.
    if geom_v_closed(s) {
        return (0.0, -1.0, 0.0);
    }
    // `cxx:891-896`.
    let surftype = if s.rectangular_trimmed_basis().is_some() {
        TypeKind::OtherSurface
    } else {
        crate::geom_bnd_lib_surface3d::surface_kind(s)
    };
    let mut close = f64::MAX;
    let mut mid = -1.0;
    match surftype {
        // `cxx:900-906`.
        TypeKind::Plane
        | TypeKind::Cone
        | TypeKind::Cylinder
        | TypeKind::Sphere
        | TypeKind::SurfaceOfExtrusion => {}
        // `cxx:908-916`.
        TypeKind::SurfaceOfRevolution => {
            if let Some(crv) = s.revolution_basis_curve() {
                let p1 = crv.d0(crv.first_parameter());
                let p2 = crv.d0(crv.last_parameter());
                close = p1.square_distance(&p2);
            }
        }
        // `cxx:917-984`.
        TypeKind::BSplineSurface => {
            if let Some(bs) = s.osculating_bspline() {
                let nbvp = bs.nb_poles_v();
                let mut distmin = f64::MAX;
                let (m1, m2) = end_mults(&bs.knots_v);
                if bs.v_periodic {
                    close = 0.0;
                    delt = 0.0;
                } else if nbvp < 3 {
                    // `cxx:928-931`.
                } else if bs.is_rational()
                    || m1 != bs.deg_v as i32 + 1
                    || m2 != bs.deg_v as i32 + 1
                {
                    // `cxx:932-960`.
                    let (uknots, _) = occt_geom::bspline_surface::GeomBSplineSurface::unique_knots_mults(&bs.knots_u);
                    let nbuk = uknots.len();
                    let mut u = uknots[0];
                    let p1 = s.value(u, vf);
                    let p2 = s.value(u, vl);
                    close = p1.square_distance(&p2);
                    let pm = s.value(u, 0.5 * (vf + vl));
                    mid = p1.square_distance(&pm);
                    distmin = close;
                    for i in 1..nbuk {
                        u = 0.5 * (uknots[i - 1] + uknots[i]);
                        let p1 = s.value(u, vf);
                        let p2 = s.value(u, vl);
                        let a_dist = p1.square_distance(&p2);
                        if a_dist > close {
                            close = a_dist;
                            let pm = s.value(u, 0.5 * (vf + vl));
                            mid = p1.square_distance(&pm);
                        } else {
                            distmin = distmin.min(a_dist);
                        }
                    }
                    // `cxx:958-959`.
                    delt = delt.min(0.5 * occt_geom::approx_same_parameter::v_resolution(s, distmin.sqrt()));
                } else {
                    // `cxx:962-983`: the boundary pole columns.
                    let nbup = bs.poles.len();
                    close = bs.poles[0][0].square_distance(&bs.poles[0][nbvp - 1]);
                    mid = bs.poles[0][0].square_distance(&bs.poles[0][nbvp / 2]);
                    distmin = close;
                    for i in 1..nbup {
                        let a_dist = bs.poles[i][0].square_distance(&bs.poles[i][nbvp - 1]);
                        if a_dist > close {
                            close = a_dist;
                            mid = bs.poles[i][0].square_distance(&bs.poles[i][nbvp / 2]);
                        } else {
                            distmin = distmin.min(a_dist);
                        }
                    }
                    // `cxx:980-981`.
                    delt = delt.min(0.5 * occt_geom::approx_same_parameter::v_resolution(s, distmin.sqrt()));
                }
            }
        }
        // `cxx:985-1015`.
        TypeKind::BezierSurface => {
            if let Some(bz) = s.osculating_bspline() {
                let nbvp = bz.nb_poles_v();
                if nbvp >= 3 {
                    let nbup = bz.nb_poles_u();
                    close = bz.poles[0][0].square_distance(&bz.poles[0][nbvp - 1]);
                    mid = bz.poles[0][0].square_distance(&bz.poles[0][nbvp / 2]);
                    let mut distmin = close;
                    for i in 1..nbup {
                        let a_dist = bz.poles[i][0].square_distance(&bz.poles[i][nbvp - 1]);
                        if a_dist > close {
                            close = a_dist;
                            mid = bz.poles[i][0].square_distance(&bz.poles[i][nbvp / 2]);
                        } else {
                            distmin = distmin.min(a_dist);
                        }
                    }
                    // `cxx:1012-1013`.
                    delt = delt.min(0.5 * occt_geom::approx_same_parameter::v_resolution(s, distmin.sqrt()));
                }
            }
        }
        // `cxx:1017-1046`.
        _ => {
            let nbpoints = 101i32;
            let p1 = s.value(uf, vf);
            let p2 = s.value(uf, vl);
            let pm = s.value(uf, 0.5 * (vf + vl));
            close = p1.square_distance(&p2);
            mid = p1.square_distance(&pm);
            let mut distmin = close;
            for i in 1..nbpoints {
                let uparam = uf + (ul - uf) * i as f64 / (nbpoints - 1) as f64;
                let p1 = s.value(uparam, vf);
                let p2 = s.value(uparam, vl);
                let a_dist = p1.square_distance(&p2);
                if a_dist > close {
                    close = a_dist;
                    let pm = s.value(uparam, 0.5 * (vf + vl));
                    mid = p1.square_distance(&pm);
                } else {
                    distmin = distmin.min(a_dist);
                }
            }
            // `cxx:1043-1044`.
            delt = delt.min(0.5 * occt_geom::approx_same_parameter::v_resolution(s, distmin.sqrt()));
        }
    }
    // `cxx:1048-1049`.
    close = close.sqrt();
    // `cxx:1052-1056`.
    if mid > 0.0 && close > mid.sqrt() {
        return (f64::MAX, -1.0, delt);
    }
    (close, mid, delt)
}

/// `ShapeAnalysis_Surface::IsUClosed(preci)` (`ShapeAnalysis_Surface.cxx:660-863`):
/// `myUCloseVal <= max(preci, Precision::Confusion())` (`cxx:658`, `:863`).
pub(crate) fn sa_is_u_closed(s: &dyn Surface, preci: f64) -> bool {
    sa_u_close_state(s).0 <= preci.max(occt_core::precision::CONFUSION)
}

/// `ShapeAnalysis_Surface::IsVClosed(preci)` (`ShapeAnalysis_Surface.cxx:865-1077`).
pub(crate) fn sa_is_v_closed(s: &dyn Surface, preci: f64) -> bool {
    sa_v_close_state(s).0 <= preci.max(occt_core::precision::CONFUSION)
}

/// `myUDelt` as `ShapeAnalysis_Surface::ValueOfUV` leaves it before reading
/// `du = min(myUDelt, SurfAdapt.UResolution(preci))` (`cxx:1332-1346`): the
/// first `ValueOfUV` of a general-arm surface forces `IsUClosed()` when
/// `myUCloseVal < 0`, so the field always carries the state computed here.
pub(crate) fn sa_u_delt(s: &dyn Surface) -> f64 {
    sa_u_close_state(s).2
}

/// `myVDelt`; see [`sa_u_delt`] and `cxx:1336-1347`.
pub(crate) fn sa_v_delt(s: &dyn Surface) -> f64 {
    sa_v_close_state(s).2
}

