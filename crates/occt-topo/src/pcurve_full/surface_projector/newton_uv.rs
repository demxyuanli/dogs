//! `ShapeAnalysis_Surface::ValueOfUV` (`ShapeAnalysis_Surface.cxx:1245-1515`) with the `SurfaceNewton` polish (`cxx:1065-1156`) and `NextValueOfUV`.

use super::*;

/// `ShapeAnalysis_Surface::ValueOfUV` (`ShapeAnalysis_Surface.cxx:1245-1515`):
/// the closest `(u, v)` of a 3D point on a surface.
///
/// `cxx:1262-1293` is a type switch on `Adaptor3d_Surface::GetType()`: for the
/// five analytic surfaces the answer is `ElSLib::Parameters` directly, with
/// `ShapeAnalysis::AdjustByPeriod` applied to each periodic parameter so the
/// value lands next to the middle of the surface range (`cxx:1271` cylinder,
/// `cxx:1276` cone, `cxx:1281` sphere, `cxx:1287-1290` torus U and V). The
/// switch is on the *type*, so a `Geom_BSplineSurface` that is only
/// geometrically a cylinder must not take these arms.
///
/// `cxx:1294-1476` is the general arm: `Extrema_ExtPS` seeded over the bounds
/// (`cxx:1346-1352`, extended by one resolution) and then polished with
/// `SurfaceNewton` (`cxx:1390-1401`).
/// [`occt_geom::extrema_surf`] ports both (`Extrema_GenExtPS` grid seeds +
/// Newton on `F = ((S-P)·Su, (S-P)·Sv)`), replacing the removed coarse
/// grid + axis-descent projector, whose
/// coordinate descent stalls ~1e-3 short along a curved valley.
///
/// This is the `mpt[i] == 3`-independent numeric inversion the whole
/// `approxPCurve` population loop (`cxx:1310-1464`) and `getLine`
/// (`cxx:902-1102`) run on. `ShapeAnalysis_Surface::NextValueOfUV`
/// (`cxx:1164-1243`) is the same solve seeded from the previous `(u, v)`; the
/// extrema path is used here because it is global (no seed continuity needed
/// for the sampled point sequence the caller builds).
///
/// PORTED (this function): the window extension
/// `du = min(myUDelt, UResolution(preci))` /
/// `dv = min(myVDelt, VResolution(preci))` (`cxx:1340-1348`, both zero for
/// offset surfaces, `myUDelt = myVDelt = 0.01` at `cxx:106-107`), fed to
/// `Extrema_ExtPS` over `[uf - du, ul + du] x [vf - dv, vl + dv]`
/// (`cxx:1352`), and `RestrictBounds` (`cxx:1323`, `cxx:64-91`) before it.
/// The closest point is therefore allowed to leave the face's parameter range,
/// which is what OCCT does when the minimum sits on the domain border.
///
/// PORTED (this function): the residue branch after the extrema minimum
/// (`cxx:1381-1437`) - the extrema distance is recomputed from the actual
/// surface value (`cxx:1374-1378`, the OCC486 workaround), the local-minimum
/// normal test `possLockal` (`cxx:1402-1419`), and [`uv_from_iso`]
/// (`cxx:1423` / `cxx:1459`, both in place of a general closest-point search),
/// which wins the result when it lands closer (`cxx:1428-1437`).
///
/// NOT PORTED (reported):
/// - `myExtOK` caching (`cxx:1327`): OCCT builds the `Extrema_ExtPS` object
///   once per surface; [`occt_geom::extrema_surf::point_surface_extrema_box`]
///   re-seeds its grid on every call. The result is the same, only the cost
///   differs. `value_of_uv` is a free function, so no cache state exists here.
/// - `ComputeBoxes` (`cxx:1898-1922`) and the `anIsoBox->Distance(aPBox) >
///   theMin` early-outs (`cxx:1607-1610`, `cxx:1676-1678`): they skip a
///   candidate whose bounding box is already farther than the best distance, so
///   they cannot change the result - only the cost.
///
/// PORTED (this function): `ShapeAnalysis_Surface::SurfaceNewton`
/// (`ShapeAnalysis_Surface.cxx:1065-1156`) at both of its `ValueOfUV` call
/// sites - the conic-extrusion attempt (`cxx:1309`) and the residue polish
/// (`cxx:1390-1401`).
pub(crate) fn value_of_uv(s: &dyn Surface, p: &GpPnt, preci: f64) -> GpPnt2d {
    // `cxx:1254` `Bounds(uf, ul, vf, vl)` (`ShapeAnalysis_Surface.lxx:54-64`,
    // the surface's own parameter bounds buffered at construction).
    let (uf, ul) = s.u_range();
    let (vf, vl) = s.v_range();
    let two_pi = 2.0 * PI;
    // `cxx:1262-1293` reads `SurfAdapt.Cylinder()` / `Cone()` / `Sphere()` /
    // `Torus()`, and `GeomAdaptor_Surface::load` (`GeomAdaptor_Surface.cxx:423-430`)
    // answers those from the basis of a `Geom_RectangularTrimmedSurface`. The
    // period-adjustment target below stays the trimmed surface's own window
    // (`ShapeAnalysis_Surface.lxx:54-63`).
    let trimmed_basis = analytic_basis(s);
    let analytic: &dyn Surface = trimmed_basis.as_deref().unwrap_or(s);
    match classify_surface_kind(s) {
        // `cxx:1265-1268`: no periodic parameter, so no adjustment.
        SurfaceKind::Plane => {
            if let Some(pl) = analytic.gp_pln() {
                let (u, v) = occt_core::elib::slib::plane_parameters(&pl.pos, p);
                return GpPnt2d::new(u, v);
            }
        }
        // `cxx:1269-1275`.
        SurfaceKind::Cylinder => {
            if let Some(cy) = analytic.gp_cylinder() {
                let (mut u, v) = occt_core::elib::slib::cylinder_parameters(&cy.pos, p);
                u += crate::shhealing::adjust_by_period(u, 0.5 * (uf + ul), two_pi);
                return GpPnt2d::new(u, v);
            }
        }
        // `cxx:1276-1280`.
        SurfaceKind::Cone => {
            if let Some(co) = analytic.gp_cone() {
                let (mut u, v) =
                    occt_core::elib::slib::cone_parameters(&co.pos, co.radius, co.semi_angle, p);
                u += crate::shhealing::adjust_by_period(u, 0.5 * (uf + ul), two_pi);
                return GpPnt2d::new(u, v);
            }
        }
        // `cxx:1281-1285`.
        SurfaceKind::Sphere => {
            if let Some(sp) = analytic.gp_sphere() {
                let (mut u, v) = occt_core::elib::slib::sphere_parameters(&sp.pos, p);
                u += crate::shhealing::adjust_by_period(u, 0.5 * (uf + ul), two_pi);
                return GpPnt2d::new(u, v);
            }
        }
        // `cxx:1286-1293`, the only arm that adjusts V as well.
        SurfaceKind::Torus => {
            if let Some(to) = analytic.gp_torus() {
                let (mut u, mut v) = occt_core::elib::slib::torus_parameters(
                    &to.pos,
                    to.major_radius,
                    to.minor_radius,
                    p,
                );
                u += crate::shhealing::adjust_by_period(u, 0.5 * (uf + ul), two_pi);
                v += crate::shhealing::adjust_by_period(v, 0.5 * (vf + vl), two_pi);
                return GpPnt2d::new(u, v);
            }
        }
        SurfaceKind::Other => {}
    }
    // `cxx:1300-1301` `S = (uf + ul) / 2; T = (vf + vl) / 2;`: the general arm
    // starts from the middle of the window, which is also the answer it returns
    // when the extrema fail (`cxx:1454-1471`).
    let mut uf = uf;
    let mut ul = ul;
    let mut vf = vf;
    let mut vl = vl;
    let mut s_par = 0.5 * (uf + ul);
    let mut t_par = 0.5 * (vf + vl);
    // `cxx:1306-1321`: a `Geom_SurfaceOfExtrusion` with an infinite U range (a
    // conic extrusion - the basis curve is a line) is first solved by
    // `SurfaceNewton` from the middle of the window; on failure `uf, ul` become
    // `-500, 500` before `RestrictBounds`, which would otherwise take the
    // width-2000 window from the finite end.
    if s.is_surface_of_linear_extrusion()
        && occt_core::precision::Precision::is_infinite(uf)
        && occt_core::precision::Precision::is_infinite(ul)
    {
        // `cxx:1309-1313`: the Newton solution is returned outright.
        if let Some((solution, _)) =
            surface_newton(s, &GpPnt2d::new(s_par, t_par), p, preci)
        {
            return solution;
        }
        uf = -500.0;
        ul = 500.0;
    }
    // `cxx:1324` `RestrictBounds(uf, ul, vf, vl)`: an infinite range becomes
    // `[-1000, 1000]`, or a width-2000 window from the finite end.
    restrict_bounds(&mut uf, &mut ul);
    restrict_bounds(&mut vf, &mut vl);
    // `cxx:1332-1348`: the search window is the (restricted) bounds expanded by
    // `du = min(myUDelt, UResolution(preci))` / `dv = min(myVDelt, VResolution(preci))`;
    // the expansion is zero for an offset surface because evaluating one
    // outside its range can be undefined (id23943).
    //
    // `myUDelt` / `myVDelt` are NOT the constructor values
    // (`ShapeAnalysis_Surface.cxx:106-107`, `0.01`): the first general-arm
    // `ValueOfUV` forces `IsUClosed()` / `IsVClosed()` (`cxx:1332-1338`) which
    // overwrite them with `|ul - uf| / 20` (`cxx:672`), reduced by half the
    // resolution at the minimum boundary gap (`cxx:764`, `:786`, `:818`,
    // `:849`; `:959`, `:981`, `:1013`, `:1044`, `0` when the surface's own
    // `IsUClosed` already says closed).
    let mut du = 0.0;
    let mut dv = 0.0;
    if !s.is_offset_surface() {
        // `GeomAdaptor_Surface::UResolution` / `VResolution`
        // (`GeomAdaptor_Surface.cxx:1818`, `:1900`).
        du = sa_u_delt(s).min(occt_geom::approx_same_parameter::u_resolution(s, preci));
        dv = sa_v_delt(s).min(occt_geom::approx_same_parameter::v_resolution(s, preci));
    }
    let e = occt_geom::extrema_surf::point_surface_extrema_box(
        s,
        p,
        uf - du,
        ul + du,
        vf - dv,
        vl + dv,
    );
    let Some(v_min) = e.v2.filter(|v| e.u2.is_finite() && v.is_finite()) else {
        // `cxx:1448-1472`: no extrema solution, so `myGap = UVFromIso(P3D,
        // preci, UU, VV)` from the middle-of-window seed wins outright.
        let (_, iu, iv) = uv_from_iso(s, p, preci, s_par, t_par);
        return GpPnt2d::new(iu, iv);
    };
    s_par = e.u2;
    t_par = v_min;
    // `cxx:1376-1378`: the distance is recomputed from the actual surface value
    // (the OCC486 workaround - a surface of revolution can return a stale
    // distance from its extrema solution).
    let mut dis_surf = p.distance(&s.d0(s_par, t_par));
    if dis_surf > preci {
        // `cxx:1390-1401` PRO7226 #412920: the extrema foot is polished by
        // `SurfaceNewton` and the closer of the two keeps the answer.
        if let Some((pp, _)) = surface_newton(s, &GpPnt2d::new(s_par, t_par), p, preci) {
            let dist = p.distance(&s.value(pp.x(), pp.y()));
            if dist < dis_surf {
                dis_surf = dist;
                s_par = pp.x();
                t_par = pp.y();
            }
        }
        // `cxx:1402-1419`: a residue below `10 * preci` whose foot already lies
        // on the normal through `P3D` (`possLockal`) is a genuine local minimum,
        // and the iso search below is skipped for it.
        let mut poss_lockal = false;
        if dis_surf < 10.0 * preci && s.continuity() != 0 {
            // `GeomAbs_C0 == 0` (`GeomAbs_Shape`).
            let (pnt, d1u, d1v) = s.d1(s_par, t_par);
            let b = d1u.crossed(&d1v);
            let a = GpVec::from_pnts(&pnt, p);
            let ab = a.dot(&b);
            let nrm2 = b.square_magnitude();
            if nrm2 > 1e-10 {
                poss_lockal = a.square_magnitude() - (ab * ab) / nrm2
                    < occt_core::precision::CONFUSION * occt_core::precision::CONFUSION;
            }
        }
        if !poss_lockal {
            // `cxx:1423` (and `cxx:1459` on the failure path above).
            let (iso_dist, iu, iv) = uv_from_iso(s, p, preci, s_par, t_par);
            // `cxx:1428-1437`: the iso solution replaces the extrema one only
            // when it is closer.
            if dis_surf > iso_dist {
                s_par = iu;
                t_par = iv;
            }
        }
    }
    GpPnt2d::new(s_par, t_par)
}

/// [`value_of_uv`] plus the `myGap` that `ShapeAnalysis_Surface::ValueOfUV`
/// leaves in `Gap()`, which `ShapeConstruct_ProjectCurveOnSurface::getLine`
/// reads back at `cxx:973` / `:987` to grow the running `aTol2`.
///
/// `ValueOfUV` opens with `myGap = -1.` (`ShapeAnalysis_Surface.cxx:1249`) and
/// closes with `if (myGap <= 0) myGap = P3D.Distance(SurfAdapt.Value(S, T));`
/// (`cxx:1502-1513`). The five analytic arms (`cxx:1262-1293`) never assign it,
/// so the tail fills it with the distance between the sample and the surface
/// point of the returned parameters; the general arm stores the same quantity
/// as `min(disSurf, DistMinOnIso)` (`cxx:1432-1437`) or the `UVFromIso`
/// minimum (`cxx:1459`). The one path that returns with `myGap = -1.` is the
/// infinite-`U` `GeomAbs_SurfaceOfExtrusion` Newton arm (`cxx:1306-1313`),
/// which leaves before the tail - see [`conic_extrusion_newton_arm`].
pub(crate) fn value_of_uv_with_gap(s: &dyn Surface, p: &GpPnt, preci: f64) -> (GpPnt2d, f64) {
    let q = value_of_uv(s, p, preci);
    let gap = if conic_extrusion_newton_arm(s, p, preci) {
        -1.0
    } else {
        s.d0(q.x(), q.y()).distance(p)
    };
    (q, gap)
}

/// `ShapeAnalysis_Surface::ValueOfUV` (`cxx:1306-1313`): a
/// `Geom_SurfaceOfExtrusion` with an infinite `U` range (the basis curve is a
/// line - the conic case) is solved by `SurfaceNewton` from the middle of the
/// window and `return solution;` leaves the function before the `myGap` tail
/// (`cxx:1502-1513`), so `Gap()` keeps the `-1.` of `cxx:1249`.
/// On Newton failure the arm continues into the general path and the tail runs.
pub(crate) fn conic_extrusion_newton_arm(s: &dyn Surface, p: &GpPnt, preci: f64) -> bool {
    if !s.is_surface_of_linear_extrusion() {
        return false;
    }
    let (uf, ul) = s.u_range();
    if !(occt_core::precision::Precision::is_infinite(uf)
        && occt_core::precision::Precision::is_infinite(ul))
    {
        return false;
    }
    let (vf, vl) = s.v_range();
    surface_newton(
        s,
        &GpPnt2d::new(0.5 * (uf + ul), 0.5 * (vf + vl)),
        p,
        preci,
    )
    .is_some()
}

/// `ShapeAnalysis_Surface::SurfaceNewton` (`ShapeAnalysis_Surface.cxx:1065-1156`).
///
/// The Newton step of the stationary-point system
/// `F = ((S(u, v) - P3D) . Su, (S(u, v) - P3D) . Sv)` solved from the second
/// derivatives `SurfAdapt.D2` (`cxx:1088`). PORTED here: the `D2` normal test
/// (`cxx:1091-1095`), the discriminant and step (`cxx:1101-1115`), the window
/// rejection `[uf - du, ul + du] x [vf - dv, vl + dv]` (`cxx:1117-1120`, the
/// window is the one `Bounds` + `UResolution`/`VResolution` give, `cxx:1073-1077`),
/// the convergence test `|du| + |dv| <= max(1e-12, (U + V) * 10e-16)`
/// (`cxx:1126-1130`) and the two solution rejections (`cxx:1134-1137`
/// `rs2 > rsfirst`, `cxx:1140-1143` `rs2 - rsn^2 / nrm2 > Tol^2`).
///
/// Returns `(sol, res)`; `res == 2` is OCCT's "strange attractor" guard
/// (`cxx:1151` `nrm2 < 0.01 * ru2 * rv2`) which `NextValueOfUV` treats as
/// "never trust this solution without the iso fallback" (`cxx:1220`). `None` is
/// the C++ `false` return (0).
pub(crate) fn surface_newton(
    s: &dyn Surface,
    p2d_prev: &GpPnt2d,
    p3d: &GpPnt,
    preci: f64,
) -> Option<(GpPnt2d, i32)> {
    // `cxx:1072-1077`: `Bounds` widened by `SurfAdapt.UResolution(preci)` /
    // `VResolution(preci)` (`GeomAdaptor_Surface::UResolution`, the port in
    // `occt_geom::approx_same_parameter`). NOTE: the surface's own
    // `uv_resolution` is `None` for the analytic family (plane, cylinder, cone,
    // sphere, torus), which is why the resolution port is called directly here
    // instead of falling back to `Precision::Parametric(preci)`.
    let (uf, ul) = s.u_range();
    let (vf, vl) = s.v_range();
    let du_res = occt_geom::approx_same_parameter::u_resolution(s, preci);
    let dv_res = occt_geom::approx_same_parameter::v_resolution(s, preci);
    let u_f = uf - du_res;
    let u_l = ul + du_res;
    let v_f = vf - dv_res;
    let v_l = vl + dv_res;
    // `cxx:1079-1080`
    let tol = occt_core::precision::CONFUSION;
    let tol2 = tol * tol;
    let mut u = p2d_prev.x();
    let mut v = p2d_prev.y();
    // `cxx:1082` `rsfirst = P3D.XYZ() - Value(U, V).XYZ()`
    let rs_first = GpVec::from_pnts(&s.value(u, v), p3d);
    for _ in 0..25 {
        // `cxx:1088` `SurfAdapt.D2(U, V, pnt, ru, rv, ruu, rvv, ruv)`.
        let (pnt, ru, rv, ruu, rvv, ruv) = s.d2(u, v);
        // `cxx:1090-1096`.
        let ru2 = ru.dot(&ru);
        let rv2 = rv.dot(&rv);
        let n = ru.crossed(&rv);
        let nrm2 = n.square_magnitude();
        if nrm2 < 1e-10 || occt_core::precision::Precision::is_positive_infinite(nrm2) {
            break;
        }
        // `cxx:1098-1107`.
        let rs = GpVec::from_pnts(&pnt, p3d);
        let r_suu = rs.dot(&ruu);
        let r_svv = rs.dot(&rvv);
        let r_suv = rs.dot(&ruv);
        let d = -nrm2 + rv2 * r_suu + ru2 * r_svv - 2.0 * r_suv * ru.dot(&rv)
            + r_suv * r_suv
            - r_suu * r_svv;
        if d.abs() < 1e-10 {
            break;
        }
        // `cxx:1110-1115`.
        let fract = 1.0 / d;
        let step_u = rs.dot(
            &n.crossed(&rv)
                .added(&ru.multiplied_scalar(r_svv))
                .subtracted(&rv.multiplied_scalar(r_suv)),
        ) * fract;
        let step_v = rs.dot(
            &ru.crossed(&n)
                .added(&rv.multiplied_scalar(r_suu))
                .subtracted(&ru.multiplied_scalar(r_suv)),
        ) * fract;
        u += step_u;
        v += step_v;
        // `cxx:1117-1120`.
        if u < u_f || u > u_l || v < v_f || v > v_l {
            break;
        }
        // `cxx:1126-1130`.
        let a_resolution = 1e-12f64.max((u + v) * 10e-16);
        if step_u.abs() + step_v.abs() > a_resolution {
            continue;
        }
        // `cxx:1133-1137` PRO10109 4517.
        let rs2 = rs.square_magnitude();
        if rs2 > rs_first.square_magnitude() {
            break;
        }
        // `cxx:1140-1143`.
        let r_sn = rs.dot(&n);
        if rs2 - r_sn * r_sn / nrm2 > tol2 {
            break;
        }
        // `cxx:1147-1153`.
        let res = if nrm2 < 0.01 * ru2 * rv2 { 2 } else { 1 };
        return Some((GpPnt2d::new(u, v), res));
    }
    None
}

/// `GeomAdaptor_Surface::UContinuity` / `VContinuity` for a B-spline
/// (`GeomAdaptor_Surface.cxx:539-587` / `:591-639`) through
/// `LocalContinuity` (`cxx:110-165`): the direction is `GeomAbs_C0` when the
/// largest multiplicity of the unique knots strictly inside the window reaches
/// the degree (`cxx:139-153` wraps the end knots out of the range with the
/// `EpsKnot` nudges at `cxx:124-132`).
pub(crate) fn bspline_c0(direction_knots: &[f64], degree: usize, first: f64, last: f64) -> bool {
    let (knots, mults) = occt_geom::bspline_surface::GeomBSplineSurface::unique_knots_mults(
        direction_knots,
    );
    if knots.len() < 3 {
        return false;
    }
    let mut mult_max = 0i32;
    for (k, m) in knots.iter().zip(mults.iter()).skip(1).take(knots.len() - 2) {
        if *k > first && *k < last && *m > mult_max {
            mult_max = *m;
        }
    }
    (degree as i32) - mult_max <= 0
}

/// The `Geom_BSplineSurface` behind `s`, peeling `Geom_RectangularTrimmedSurface`
/// wrappers (`GeomAdaptor_Surface::load` dispatches a trimmed surface to its
/// basis type, so `SurfAdapt.BSpline()` is the basis, `Surface.cxx:423-430`).
pub(crate) fn bspline_behind(s: &dyn Surface, depth: usize) -> Option<occt_geom::bspline_surface::GeomBSplineSurface> {
    if let Some(basis) = s.rectangular_trimmed_basis() {
        if depth < 8 {
            return bspline_behind(basis.as_ref(), depth + 1);
        }
        return None;
    }
    s.osculating_bspline()
}

/// `ShapeAnalysis_Surface::NextValueOfUV` (`ShapeAnalysis_Surface.cxx:1164-1243`):
/// `SurfaceNewton` seeded from the previous `(u, v)` (the hint), used both by
/// `ShapeConstruct_ProjectCurveOnSurface::projectPoint` (`cxx:953`, `cxx:1410`)
/// and by the sample walk (`cxx:1358`, `cxx:1432`).
///
/// PORTED (this function): the type gate (`cxx:1174-1180`: `SurfaceNewton` is
/// only reachable for Bezier / B-spline / extrusion / revolution / offset
/// surfaces, every other type jumps straight to `ValueOfUV`, `cxx:1240`), the
/// C0-knot bail-out for a B-spline whose hint sits on a C0 knot
/// (`cxx:1182-1213`), and the `maxpreci` guard of the accepted solution
/// (`cxx:1216-1235`): `res == 2`, or a gap more than `Confusion` above
/// `maxpreci`, sends the result through `UVFromIso` and keeps it only when it is
/// at least as close as the iso minimum.
///
/// Returns `(uv, gap)` with the `myGap` OCCT leaves behind (`cxx:1226` / `:1232`
/// / through `ValueOfUV` at `cxx:1240`). The `ValueOfUV` fall-through returns
/// the distance recomputed at the returned `(u, v)`; see the note on
/// [`value_of_uv_with_gap`].
pub(crate) fn next_value_of_uv(
    s: &dyn Surface,
    p2d_prev: &GpPnt2d,
    p3d: &GpPnt,
    preci: f64,
    maxpreci: f64,
) -> (GpPnt2d, f64) {
    use crate::geom_bnd_lib_surface3d::SurfaceKind as TypeKind;
    let takes_newton = matches!(
        crate::geom_bnd_lib_surface3d::surface_kind(s),
        TypeKind::BezierSurface
            | TypeKind::BSplineSurface
            | TypeKind::SurfaceOfExtrusion
            | TypeKind::SurfaceOfRevolution
            | TypeKind::OffsetSurface
    );
    if takes_newton {
        // `cxx:1182-1213`: near a C0 knot the Newton step is unreliable, so the
        // global solve is used instead.
        if let Some(bs) = bspline_behind(s, 0) {
            let (uf, ul) = s.u_range();
            let (vf, vl) = s.v_range();
            let on_u_knot = bspline_c0(&bs.knots_u, bs.deg_u, uf, ul)
                && bs
                    .knots_u
                    .iter()
                    .any(|k| (k - p2d_prev.x()).abs() < occt_core::precision::CONFUSION);
            let on_v_knot = bspline_c0(&bs.knots_v, bs.deg_v, vf, vl)
                && bs
                    .knots_v
                    .iter()
                    .any(|k| (k - p2d_prev.y()).abs() < occt_core::precision::CONFUSION);
            if on_u_knot || on_v_knot {
                return value_of_uv_with_gap(s, p3d, preci);
            }
        }
        // `cxx:1215-1237`.
        if let Some((sol, res)) = surface_newton(s, p2d_prev, p3d, preci) {
            let gap = p3d.distance(&s.value(sol.x(), sol.y()));
            if res == 2 || (maxpreci > 0.0 && gap - maxpreci > occt_core::precision::CONFUSION) {
                // `cxx:1226-1230`: `UVFromIso` writes back into `U` / `V`.
                let (my_gap, u, v) = uv_from_iso(s, p3d, preci, sol.x(), sol.y());
                if gap >= my_gap {
                    return (GpPnt2d::new(u, v), my_gap);
                }
            }
            return (sol, gap);
        }
    }
    // `cxx:1240`.
    value_of_uv_with_gap(s, p3d, preci)
}

