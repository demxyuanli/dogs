//! Isoparametric-edge detection and projection: `isAnIsoparametric` plus `ShapeAnalysis_Curve::Project` on the iso

use super::*;

/// `ShapeAnalysis_Curve::Project(C3D, P3D, preci, proj, param, cf, cl, false)`
/// (`ShapeAnalysis_Curve.cxx:147-201`): an endpoint within
/// `Precision::Confusion()` is returned directly (`cxx:164-182`; the tolerance
/// is `Confusion` and not `preci` precisely because `AdjustToEnds` is false),
/// otherwise the range of a non-closed curve is widened by
/// `min(Resolution(preci), 0.1 * (umax - umin))` (`cxx:192-197`) and the minimum
/// over it is `ProjectAct` (`cxx:265-497`), the `Extrema_ExtPC` scan.
/// `ProjectAct(Adaptor3d_Curve, ..)` (`cxx:205-261`, taken by the offset arm at
/// `cxx:1806`) computes the same minimum, so one helper serves both arms and
/// `Adaptor3d_IsoCurve` is not needed.
///
/// Returns `(other, dist)` - OCCT's `param` / `distmin`.
///
/// UNPORTED (reported): `Geom_Curve::IsClosed()` is approximated by
/// `Curve::is_periodic()` plus coincident ends, because the trait exposes no
/// `IsClosed` (`Geom_Curve.hxx:155` is pure virtual and each curve implements
/// it, e.g. `Geom_BSplineCurve` through its first/last pole).
pub(crate) fn project_on_curve_range(
    iso: &dyn Curve,
    p: &GpPnt,
    preci: f64,
    cf: f64,
    cl: f64,
) -> Option<(f64, f64)> {
    let (u_min, u_max) = if cf < cl { (cf, cl) } else { (cl, cf) };
    let low = iso.d0(u_min);
    let d_low = low.distance(p);
    if d_low <= occt_core::precision::CONFUSION {
        return Some((u_min, d_low));
    }
    let high = iso.d0(u_max);
    let d_high = high.distance(p);
    if d_high <= occt_core::precision::CONFUSION {
        return Some((u_max, d_high));
    }
    let (mut a, mut b) = (u_min, u_max);
    let closed = iso.is_periodic() || low.distance(&high) <= occt_core::precision::CONFUSION;
    if !closed {
        let delta = iso.resolution(preci).min((u_max - u_min) * 0.1);
        a -= delta;
        b += delta;
    }
    let (t, _pt, d) = occt_geom::extrema_pc::extrema_ext_pc_min_in_range(iso, p, a, b)?;
    Some((t, d))
}

/// `ShapeAnalysis_Curve::NextProject` (`ShapeAnalysis_Curve.cxx:504-558`, the
/// `Geom_Curve` overload): the Newton walk `Extrema_LocateExtPC` from
/// `prev` (`cxx:572`), falling back to [`project_on_curve_range`] (OCCT
/// `Project(Adaptor3d_Curve,...)`, `cxx:205-261`) when the walk reports
/// `!IsDone`.
///
/// `AdjustToEnds` is false at the only call site (`cxx:2759`), so the bounded
/// endpoint snap uses `Precision::Confusion()` (`cxx:520-538`). As in OCCT,
/// a non-closed iso has its range widened by
/// `min(Resolution(preci), 0.1 * (umax - umin))` (`cxx:551-555`).
///
/// Returns `(param, dist)` - OCCT's `param` / the returned distance.
pub(crate) fn next_project_on_curve_range(
    iso: &dyn Curve,
    p: &GpPnt,
    preci: f64,
    prev: f64,
    cf: f64,
    cl: f64,
) -> Option<(f64, f64)> {
    let (mut u_min, mut u_max) = if cf < cl { (cf, cl) } else { (cl, cf) };
    if iso.first_parameter().is_finite() && iso.last_parameter().is_finite() {
        let low = iso.d0(u_min);
        let d_low = low.distance(p);
        if d_low <= occt_core::precision::CONFUSION {
            return Some((u_min, d_low));
        }
        let high = iso.d0(u_max);
        let d_high = high.distance(p);
        if d_high <= occt_core::precision::CONFUSION {
            return Some((u_max, d_high));
        }
    }
    let low = iso.d0(u_min);
    let high = iso.d0(u_max);
    let closed = iso.is_periodic() || low.distance(&high) <= occt_core::precision::CONFUSION;
    if !closed {
        let delta = iso.resolution(preci).min((u_max - u_min) * 0.1);
        u_min -= delta;
        u_max += delta;
    }
    // `Extrema_LocateExtPC(P3D, GAC, paramPrev, uMin, uMax, preci)` (`cxx:572`).
    if let Some((t, q)) =
        crate::int_tools_vertex_line::extrema_locate_ext_pc(iso, p, prev, u_min, u_max)
    {
        return Some((t, p.distance(&q)));
    }
    // `return Project(C3D, P3D, preci, proj, param, false)` (`cxx:578`), the
    // `Adaptor3d_Curve` overload on the already-widened `GAC` range.
    let d_low = iso.d0(u_min).distance(p);
    if d_low <= occt_core::precision::CONFUSION {
        return Some((u_min, d_low));
    }
    let d_high = iso.d0(u_max).distance(p);
    if d_high <= occt_core::precision::CONFUSION {
        return Some((u_max, d_high));
    }
    let (t, _pt, d) = occt_geom::extrema_pc::extrema_ext_pc_min_in_range(iso, p, u_min, u_max)?;
    if d < d_low + occt_core::precision::CONFUSION
        && d < d_high + occt_core::precision::CONFUSION
    {
        return Some((t, d));
    }
    if d_low < d_high {
        Some((u_min, d_low))
    } else {
        Some((u_max, d_high))
    }
}

/// `ShapeAnalysis_Surface::GetBoxUF/UL/VF/VL` (`ShapeAnalysis_Surface.cxx:1924-1946`)
/// reached through `ComputeBoxes` (`cxx:1897-1922`): the box of a bound iso
/// curve, `BndLib_Add3dCurve::Add(GeomAdaptor_Curve(iso), Precision::Confusion(),
/// B)` (`cxx:1904-1921`).
///
/// `BndLib_Add3dCurve::Add` (`BndLib_Add3dCurve.cxx:22-25`) uses the adaptor's
/// first/last parameter, i.e. the iso curve's own range; `box_curve` unwraps a
/// `Geom_TrimmedCurve` to its basis exactly like `GeomAdaptor_Curve::load`
/// (`GeomAdaptor_Curve.cxx:252-254`). The returned box carries `Gap = Confusion`
/// (`Bnd_Box.hxx:154` `Enlarge`), which is what `Bnd_Box::IsOut` tests
/// (`Bnd_Box.cxx`, `box3d.rs`).
pub(crate) fn iso_bound_box(iso: &dyn Curve) -> occt_core::bnd::BndBox {
    crate::geom_bnd_lib_curve3d::box_curve(
        iso,
        iso.first_parameter(),
        iso.last_parameter(),
        occt_core::precision::CONFUSION,
    )
}

/// Outputs of [`is_an_isoparametric`], the in-out reference parameters of
/// `ShapeConstruct_ProjectCurveOnSurface::isAnIsoparametric`
/// (`ShapeConstruct_ProjectCurveOnSurface.cxx:2461-2477`).
pub(crate) struct IsoParamOut {
    /// `theIsTypeU`: the selected iso is a `UIso` (fixed `U`).
    pub is_type_u: bool,
    /// `theP1OnIso` / `theValueP1`: the first sample lies on some bound iso.
    pub p1_on_iso: bool,
    pub value_p1: GpPnt2d,
    /// `theP2OnIso` / `theValueP2`: same for the last sample.
    pub p2_on_iso: bool,
    pub value_p2: GpPnt2d,
    /// `theIsoPar2d3d`: the 3D sample parameters equal the iso parameters, so
    /// no 2D projection is needed.
    pub iso_par2d3d: bool,
    /// `theCIso`: the selected bound iso curve.
    pub c_iso: Option<Arc<dyn Curve>>,
    /// `theT1` / `theT2`: the iso curve range corresponding to the surface
    /// window (`V1,V2` for a `UIso`, `U1,U2` for a `VIso`).
    pub t1: f64,
    pub t2: f64,
    /// `theParamsOut`: the iso parameters of the interior samples found by the
    /// `NextProject` distance walk (`cxx:2757-2768`).
    pub params_out: Vec<f64>,
}

impl IsoParamOut {
    pub(crate) fn new(n: usize) -> Self {
        Self {
            is_type_u: false,
            p1_on_iso: false,
            value_p1: GpPnt2d::zero(),
            p2_on_iso: false,
            value_p2: GpPnt2d::zero(),
            iso_par2d3d: false,
            c_iso: None,
            t1: 0.0,
            t2: 0.0,
            params_out: vec![0.0; n],
        }
    }
}

/// `ShapeConstruct_ProjectCurveOnSurface::isAnIsoparametric`
/// (`ShapeConstruct_ProjectCurveOnSurface.cxx:2461-2796`): test whether the
/// sampled 3D curve runs along one of the four boundary isos of the surface
/// parameter window, and, if so, the exact UV values that make its pcurve a
/// straight `(u, v)` line.
///
/// The four candidates are `UIso(U1)` / `UIso(U2)` / `VIso(V1)` / `VIso(V2)`
/// (`cxx:2512-2556`, `mySurf->Bounds` at `cxx:2485`; a
/// `Geom_RectangularTrimmedSurface` re-reads its own bounds at
/// `cxx:2487-2492`). For each, the first and last samples are matched against
/// the iso ends (`cxx:2588-2610`, modes 1/2) or projected on the iso by
/// `ShapeAnalysis_Curve::Project` (mode 3, `cxx:2631-2641`), with the bound
/// box `Bnd_Box::IsOut` rejection at `cxx:2619` and the spherical-`VIso` skip
/// at `cxx:2614`. The best pair (smallest sum of squared distances,
/// `cxx:2682-2696`) is kept.
///
/// When both ends land on the iso the samples are validated directly against
/// the iso 3D points (`cxx:2714-2731`, `theIsoPar2d3d`) or, failing that, by
/// the `NextProject` distance walk (`cxx:2741-2769`).
///
/// `mySurf->Surface()->Bounds(uf, ul, vf, vl)` (`cxx:1160-1161`) and
/// `mySurf->Bounds(U1, U2, V1, V2)` (`cxx:2485`) are both the port's
/// `Surface::u_range` / `v_range` (`ShapeAnalysis_Surface::Bounds`,
/// `ShapeAnalysis_Surface.lxx:54-62`, returns `myUF/myUL/myVF/myVL`). Those
/// members are copied from `mySurf->Bounds` by both constructors
/// (`ShapeAnalysis_Surface.cxx:118`, `:142`) and are written nowhere else:
/// `ShapeAnalysis_Surface::SetDomain` (`cxx:1887-1895`) has no caller anywhere
/// under `src/`, and every `RestrictBounds` call (`cxx:671`, `:878`, `:1323`,
/// `:1616`, `:1711`, `:1741`, `:1779`, `:1806`) takes local copies of the
/// members. The two bounds are therefore equal on every face, `Shape-1` (60
/// faces, 16 spherical) included.
pub(crate) fn is_an_isoparametric(
    surf: &dyn Surface,
    points: &[GpPnt],
    params: &[f64],
    out: &mut IsoParamOut,
) -> bool {
    // `cxx:2479`.
    let prec = occt_core::precision::CONFUSION;
    let nb_pnt = points.len();
    let mut iso_param = false;
    out.iso_par2d3d = false;
    let (u1, u2) = surf.u_range();
    let (v1, v2) = surf.v_range();
    // `cxx:2497-2502`.
    let mut mpt = [0i32; 2];
    let mut tpar = [0.0f64; 2];
    let mut iso_value = 0.0f64;
    let mut mindist2 = 4.0 * prec * prec;
    let mut mind2 = [mindist2; 2];
    out.p1_on_iso = false;
    out.p2_on_iso = false;
    // `cxx:2614`: `IsKind(STANDARD_TYPE(Geom_SphericalSurface))`, so a
    // `Geom_RectangularTrimmedSurface` wrapping a sphere does not qualify.
    let is_sphere = surf.gp_sphere().is_some();
    let mut chosen: Option<Arc<dyn Curve>> = None;

    for j in 0..4usize {
        // `cxx:2512-2556`.
        let (iso_u, iso_val, bound, tt1, tt2) = match j {
            0 => {
                if occt_core::precision::Precision::is_infinite(u1) {
                    continue;
                }
                (true, u1, u1, v1, v2)
            }
            1 => {
                if occt_core::precision::Precision::is_infinite(u2) {
                    continue;
                }
                (true, u2, u2, v1, v2)
            }
            2 => {
                if occt_core::precision::Precision::is_infinite(v1) {
                    continue;
                }
                (false, v1, v1, u1, u2)
            }
            _ => {
                if occt_core::precision::Precision::is_infinite(v2) {
                    continue;
                }
                (false, v2, v2, u1, u2)
            }
        };
        let ci = if iso_u {
            surf.u_iso_curve(bound)
        } else {
            surf.v_iso_curve(bound)
        };
        // `cxx:2558-2561`.
        let Some(ci) = ci else {
            continue;
        };
        let bnd_box = iso_bound_box(ci.as_ref());
        // `cxx:2578-2583`: a degenerate iso (all three probe points equal) is
        // skipped.
        let ext1 = ci.d0(tt1);
        let ext2 = ci.d0(tt2);
        let extmi = ci.d0(0.5 * (tt1 + tt2));
        if ext1.distance(&ext2) <= prec && ext1.distance(&extmi) <= prec {
            continue;
        }
        let mut pt_eq_ext1 = false;
        let mut pt_eq_ext2 = false;
        let mut currd2 = [0.0f64; 2];
        let mut tp = [0.0f64; 2];
        let mut mp = [0i32; 2];
        for i in 0..2usize {
            mp[i] = 0;
            // OCCT `k = (i == 0 ? 1 : theNbPnt)` (1-based) (`cxx:2594`).
            let k = if i == 0 { 0 } else { nb_pnt - 1 };
            currd2[i] = points[k].square_distance(&ext1);
            if currd2[i] <= prec * prec && !pt_eq_ext1 {
                mp[i] = 1;
                tp[i] = tt1;
                pt_eq_ext1 = true;
                continue;
            }
            currd2[i] = points[k].square_distance(&ext2);
            if currd2[i] <= prec * prec && !pt_eq_ext2 {
                mp[i] = 2;
                tp[i] = tt2;
                pt_eq_ext2 = true;
                continue;
            }
            // `cxx:2614-2617`: a spherical surface has no useful `VIso`
            // projection for its ends.
            if is_sphere && !iso_u {
                continue;
            }
            // `cxx:2619-2622`.
            if bnd_box.is_out(&points[k]) {
                continue;
            }
            let cf0 = ci.first_parameter();
            let cl0 = ci.last_parameter();
            let cf = if occt_core::precision::Precision::is_infinite(cf0) {
                -1000.0
            } else {
                cf0
            };
            let cl = if occt_core::precision::Precision::is_infinite(cl0) {
                1000.0
            } else {
                cl0
            };
            // `cxx:2635`: `sac.Project(cI, thePoints(k), prec, pt, t, Cf, Cl)`
            // with the default `AdjustToEnds = true`; since `prec` equals the
            // projection precision here the endpoint snap tolerance is the same
            // as the port's `AdjustToEnds = false` helper.
            if let Some((t, dist)) = project_on_curve_range(ci.as_ref(), &points[k], prec, cf, cl) {
                currd2[i] = dist * dist;
                if dist <= prec && t >= cf && t <= cl {
                    mp[i] = 3;
                    tp[i] = t;
                }
            }
        }
        // `cxx:2643-2647`.
        if mp[0] > 0 && mp[1] > 0 && (tp[0] - tp[1]).abs() < occt_core::precision::PCONFUSION {
            continue;
        }
        // `cxx:2649-2673`.
        if mp[0] > 0 && (!out.p1_on_iso || currd2[0] < mind2[0]) {
            out.p1_on_iso = true;
            mind2[0] = currd2[0];
            out.value_p1 = if iso_u {
                GpPnt2d::new(iso_val, tp[0])
            } else {
                GpPnt2d::new(tp[0], iso_val)
            };
        }
        if mp[1] > 0 && (!out.p2_on_iso || currd2[1] < mind2[1]) {
            out.p2_on_iso = true;
            mind2[1] = currd2[1];
            out.value_p2 = if iso_u {
                GpPnt2d::new(iso_val, tp[1])
            } else {
                GpPnt2d::new(tp[1], iso_val)
            };
        }
        // `cxx:2675-2696`.
        if mp[0] <= 0 || mp[1] <= 0 {
            continue;
        }
        let md2 = currd2[0] + currd2[1];
        if mindist2 <= md2 {
            continue;
        }
        mindist2 = md2;
        mpt[0] = mp[0];
        mpt[1] = mp[1];
        tpar[0] = tp[0];
        tpar[1] = tp[1];
        out.is_type_u = iso_u;
        iso_value = iso_val;
        out.t1 = tt1;
        out.t2 = tt2;
        chosen = Some(ci.clone());
    }

    // `cxx:2700-2771`.
    if mpt[0] > 0 && mpt[1] > 0 {
        let Some(ci) = chosen else {
            return false;
        };
        out.p1_on_iso = true;
        out.p2_on_iso = true;
        if out.is_type_u {
            out.value_p1 = GpPnt2d::new(iso_value, tpar[0]);
            out.value_p2 = GpPnt2d::new(iso_value, tpar[1]);
        } else {
            out.value_p1 = GpPnt2d::new(tpar[0], iso_value);
            out.value_p2 = GpPnt2d::new(tpar[1], iso_value);
        }
        if mpt[0] != 3 && mpt[1] != 3 {
            out.iso_par2d3d = true;
            for i in 1..nb_pnt.saturating_sub(1) {
                let t = if tpar[1] > tpar[0] {
                    params[i]
                } else {
                    out.t1 + out.t2 - params[i]
                };
                let pt = ci.d0(t);
                if points[i].distance(&pt) > prec {
                    out.iso_par2d3d = false;
                    break;
                }
            }
        }
        if out.iso_par2d3d {
            iso_param = true;
        } else {
            let mut prev_param = tpar[0];
            let cf0 = ci.first_parameter();
            let cl0 = ci.last_parameter();
            let cf = if occt_core::precision::Precision::is_infinite(cf0) {
                -1000.0
            } else {
                cf0
            };
            let cl = if occt_core::precision::Precision::is_infinite(cl0) {
                1000.0
            } else {
                cl0
            };
            let mut iso_by_distance = true;
            for i in 1..nb_pnt.saturating_sub(1) {
                match next_project_on_curve_range(ci.as_ref(), &points[i], prec, prev_param, cf, cl)
                {
                    Some((t, dist)) => {
                        prev_param = t;
                        out.params_out[i] = t;
                        if dist > prec || t < cf || t > cl {
                            iso_by_distance = false;
                            break;
                        }
                    }
                    None => {
                        iso_by_distance = false;
                        break;
                    }
                }
            }
            if iso_by_distance {
                iso_param = true;
            }
        }
        out.c_iso = Some(ci);
    }
    iso_param
}

/// `ShapeAnalysis_Surface::UVFromIso` (`ShapeAnalysis_Surface.cxx:1522-1842`):
/// the projection restricted to the two isos through the seed plus the four
/// window borders - candidates 0/1/2 are `UIso(myUF)` / `UIso(myUL)` / `UIso(U)`
/// and 3/4/5 `VIso(myVF)` / `VIso(myVL)` / `VIso(V)` (`cxx:1563-1603`) - then
/// refined by up to `MaxIters = 5` passes alternating the iso through the
/// running solution (`cxx:1688-1818`).
///
/// This is the "on the border" arm of `ValueOfUV`: when the extrema minimum
/// keeps a residue, the closest point of the face is often on one of the window
/// isos rather than at an interior stationary point.
///
/// Returns `(theMin, U, V)`; `u` / `v` are the seeds, matching OCCT where
/// `UU` / `VV` are in-out references and `theMin` is the returned residue.
pub(crate) fn uv_from_iso(s: &dyn Surface, p: &GpPnt, preci: f64, u: f64, v: f64) -> (f64, f64, f64) {
    // `cxx:1528-1537`: the seed itself is the first candidate.
    let mut uu = u;
    let mut vv = v;
    let mut min = s.d0(u, v).distance(p);
    if min < preci / 10.0 {
        // `cxx:1539-1542` "c etait deja OK"
        return (min, uu, vv);
    }
    // `cxx:1535-1538` reads `myUF` / `myUL` (the raw `Bounds`, not the
    // restricted window) and skips an infinite one (`cxx:1605`).
    let (uf, ul) = s.u_range();
    let (vf, vl) = s.v_range();
    let mut uv = true;
    for num in 0..6 {
        // `cxx:1563` `UV = (num < 3)`: an iso-U candidate fixes `par` to the U
        // of the iso and takes the projected parameter as V, and vice versa.
        uv = num < 3;
        let (par, iso) = match num {
            0 => (uf, s.u_iso_curve(uf)),
            1 => (ul, s.u_iso_curve(ul)),
            2 => (u, s.u_iso_curve(u)),
            3 => (vf, s.v_iso_curve(vf)),
            4 => (vl, s.v_iso_curve(vl)),
            _ => (v, s.v_iso_curve(v)),
        };
        if occt_core::precision::Precision::is_infinite(par) {
            continue;
        }
        // `cxx:1603` `iso.IsNull()`: no iso through that parameter.
        let Some(iso) = iso else {
            continue;
        };
        let (mut cf, mut cl) = (iso.first_parameter(), iso.last_parameter());
        restrict_bounds(&mut cf, &mut cl);
        if let Some((other, dist)) = project_on_curve_range(iso.as_ref(), p, preci, cf, cl) {
            if dist < min {
                min = dist;
                uu = if uv { par } else { other };
                vv = if uv { other } else { par };
            }
        }
    }
    // `cxx:1688-1818`: two projections per pass with opposite parities; OCCT
    // enters with `UV = false` (the last `UV = (num < 3)` has `num = 5`).
    let (mut prev_u, mut prev_v) = (u, v);
    let mut iters = 0;
    while (prev_u != uu || prev_v != vv) && iters < 5 && min > preci {
        prev_u = uu;
        prev_v = vv;
        for _ in 0..2 {
            let iso = if uv { s.u_iso_curve(uu) } else { s.v_iso_curve(vv) };
            if let Some(iso) = iso {
                let (mut cf, mut cl) = (iso.first_parameter(), iso.last_parameter());
                restrict_bounds(&mut cf, &mut cl);
                if let Some((other, dist)) = project_on_curve_range(iso.as_ref(), p, preci, cf, cl) {
                    if dist < min {
                        min = dist;
                        if uv {
                            vv = other;
                        } else {
                            uu = other;
                        }
                    }
                }
            }
            uv = !uv;
        }
        iters += 1;
    }
    (min, uu, vv)
}

/// `theIndCoord` accessor of `insertAdditionalPointOrAdjust` (`cxx:2039`):
/// 1 is U (the `x` coordinate), 2 is V (`thePoints2d.SetCoord(theIndCoord, ..)`).
pub(crate) fn iso_coord_of(p: &GpPnt2d, idx: usize) -> f64 {
    if idx == 1 { p.x() } else { p.y() }
}

pub(crate) fn set_iso_coord_of(p: &mut GpPnt2d, idx: usize, v: f64) {
    if idx == 1 { p.set_x(v) } else { p.set_y(v); }
}

/// `ShapeConstruct_ProjectCurveOnSurface::insertAdditionalPointOrAdjust`
/// (`ShapeConstruct_ProjectCurveOnSurface.cxx:2039-2147`): a period jump
/// between `thePrevCoord` and `theCurCoord` is either bridged by the projected
/// 3D midpoint (when the curve really crosses the seam) or removed by shifting
/// the current coordinate one period back onto the previous one.
/// `index` is 0-based here (OCCT's `theIndex` is 1-based); it is advanced past
/// the inserted midpoint. `params` / `pts2d` / `pts3d` grow together, as in
/// OCCT (`cxx:2111-2131` inserts the 3D midpoint in `thePoints` as well).
#[allow(clippy::too_many_arguments)]
pub(crate) fn insert_additional_point_or_adjust(
    to_adjust: &mut bool,
    idx_coord: usize,
    period: f64,
    tol_on_period: f64,
    cur_coord_raw: f64,
    cur_coord: &mut f64,
    prev_coord: f64,
    preci: f64,
    curve: &dyn Curve,
    index: &mut usize,
    params: &mut Vec<f64>,
    pts3d: &mut Vec<GpPnt>,
    pts2d: &mut Vec<GpPnt2d>,
    surf: &dyn Surface,
) {
    // `cxx:2053-2054`
    let corrected = occt_core::bspl::locate::in_period(
        cur_coord_raw,
        prev_coord - 0.5 * period,
        prev_coord + 0.5 * period,
    );

    if !*to_adjust {
        // `cxx:2057-2064`
        let cur_par = params[*index];
        let prev_par = params[*index - 1];
        let mut mid_par = 0.5 * (prev_par + cur_par);
        let mut mid_p3d = curve.d0(mid_par);
        let mut mid_p2d = value_of_uv(surf, &mid_p3d, preci);
        let mut mid_coord = occt_core::bspl::locate::in_period(
            iso_coord_of(&mid_p2d, idx_coord),
            prev_coord - 0.5 * period,
            prev_coord + 0.5 * period,
        );
        let mut first_coord = prev_coord;
        let mut last_coord = corrected;
        if last_coord < first_coord {
            std::mem::swap(&mut first_coord, &mut last_coord);
        }
        if last_coord - first_coord <= tol_on_period {
            // `cxx:2074-2077`
            *to_adjust = true;
        } else if first_coord <= mid_coord && mid_coord <= last_coord {
            // `cxx:2078-2081`
            *to_adjust = true;
        } else {
            // `cxx:2082-2141`: bisect in the 3D parameter until the midpoint
            // projection leaves the half-period window around `thePrevCoord`.
            let mut success = true;
            let mut first_t = prev_par;
            let mut last_t = cur_par;
            mid_coord = iso_coord_of(&mid_p2d, idx_coord);
            while (mid_coord - prev_coord).abs() >= 0.5 * period - tol_on_period
                || (cur_coord_raw - mid_coord).abs() >= 0.5 * period - tol_on_period
            {
                if mid_par - first_t <= occt_core::precision::PCONFUSION
                    || last_t - mid_par <= occt_core::precision::PCONFUSION
                {
                    success = false;
                    break;
                }
                if (mid_coord - prev_coord).abs() >= 0.5 * period - tol_on_period {
                    last_t = 0.5 * (first_t + last_t);
                } else {
                    first_t = 0.5 * (first_t + last_t);
                }
                mid_par = 0.5 * (first_t + last_t);
                mid_p3d = curve.d0(mid_par);
                mid_p2d = value_of_uv(surf, &mid_p3d, preci);
                mid_coord = iso_coord_of(&mid_p2d, idx_coord);
            }
            if success {
                // `cxx:2110-2135`: insert the extra sample at `theIndex` in all
                // three arrays.
                params.insert(*index, mid_par);
                pts3d.insert(*index, mid_p3d);
                pts2d.insert(*index, mid_p2d);
                // `cxx:2135` `theIndex++`: in OCCT the `for` loop then increments
                // `aPntIter` again (`cxx:1542`), so the next iteration lands back
                // on the shifted original sample. The Rust caller increments `i`
                // itself, so no advance is made here; that reproduces the same
                // revisit - and, with it, the same `minX` / `maxX` accumulation
                // of the original sample's coordinate.
            } else {
                *to_adjust = true;
            }
        }
    }
    if *to_adjust {
        // `cxx:2142-2146`
        *cur_coord = corrected;
        set_iso_coord_of(&mut pts2d[*index], idx_coord, corrected);
    }
}

