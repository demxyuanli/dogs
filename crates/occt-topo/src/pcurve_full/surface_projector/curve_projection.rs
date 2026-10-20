//! `ShapeConstruct_ProjectCurveOnSurface::Perform` (`cxx:707-774` onward) and the iso / analytic p-curve builders and entry points it dispatches to.

use super::*;

/// `ShapeConstruct_ProjectCurveOnSurface::Perform` (`cxx:707-774`):
/// `projectAnalytic` first (`cxx:731-736`), then `generateCurvePoints`
/// (`cxx:750-752`) into the `approxPCurve` chain (`cxx:760`). `preci` is
/// `myPreci` (`cxx:533` `myProjector->Init(sas, preci)`), which `approxPCurve`
/// forwards to `getLine` (`cxx:1119`). `tol_first` / `tol_last` are `Perform`'s
/// `theTolFirst` / `theTolLast`, the tolerances of the edge's end vertices
/// (`ShapeFix_Edge.cxx:521-531`); a negative value means `Precision::Confusion()`
/// (`cxx:1475-1476`).
///
/// The `myCache` endpoint cache and the B-spline corner cache are PORTED
/// (`pcurve_full/projection_cache.rs`): `getLine` seeds its four probes from them, anchors
/// `fixPeriodicityTroubles` on the cached point (`cxx:957-961`) and refills
/// `myCache` from the pcurve it returns (`cxx:1122-1142`); the fallback path
/// refills it at `cxx:1892-1903`. `cache` is `myCache`, cleared by `SetSurface`
/// (`cxx:676-683`) - the caller keeps one per reprojected edge chain.
///
/// The sample loop (`cxx:1311-1465`) is PORTED: interior samples walk from the
/// previous parameter point with `NextValueOfUV(p2d, p3d, myPreci,
/// Confusion + 1000 * gap)` (`cxx:1432`), the endpoints are re-projected from
/// `myCache` only when `getLine` fixed a period jump (`cxx:1401-1424`), and the
/// next seed is the extrapolation `2 * p2d - thePoints2d(prev)` (`cxx:1455-1463`).
///
/// PARK (unported, reported):
/// - `isBSplineCurveInvalid` / `PerformByProjLib` (`cxx:428`, `cxx:739-748`).
/// - the `isoParam` arms of the sample loop (`cxx:1317-1375`) and the
///   `p1OnIso` / `p2OnIso` endpoint overrides (`cxx:1379-1385`): both require
///   `isAnIsoparametric` (`cxx:1170`, impl `cxx:2461-2796`), an edge that lies
///   on a boundary iso of the surface parameter window and becomes an exact iso
///   pcurve (a straight line in `(u, v)`), with `theParamsOut` reparameterized
///   through `ShapeAnalysis_Curve::Project` / `NextProject`. Four things block
///   it:
///   1. `mySurf->UIso(U1)` / `UIso(U2)` / `VIso(V1)` / `VIso(V2)`
///      (`cxx:2519`, `:2530`, `:2540`, `:2549`). The trait's
///      `u_iso_curve` / `v_iso_curve` return `None` for the whole analytic
///      family (plane, cylinder, cone, sphere, torus), so those four arms would
///      silently never fire; implementing them means editing
///      `crates/occt-geom/src/**`, outside this port's allowed files.
///   2. `ShapeAnalysis_Surface::GetBoxUF/GetBoxUL/GetBoxVF/GetBoxVL`
///      (`cxx:2522`, `:2532`, `:2543`, `:2552`), used as the
///      `aBox->IsOut(thePoints(k))` rejection (`cxx:2619`).
///   3. `ShapeAnalysis_Curve::Project` (`ShapeAnalysis_Curve.cxx:126-202`,
///      `:205-261`, `ProjectAct` `:265-...`, extrema + Newton on a
///      `GeomAdaptor_Curve`) and `NextProject` (`:504-...`, `:561-...`); there
///      is no `ShapeAnalysis_Curve` module in this port.
///   4. `theCIso->D0` / `FirstParameter` / `LastParameter` are used on the
///      returned handle (`cxx:2578-2584`, `:2727`, `:2744`), which is available
///      through [`Curve`], but point 1 means the handle never exists.
///   Note that `theP1OnIso` / `theP2OnIso` can be set by `isAnIsoparametric`
///   even when it returns `isoParam == false` (`cxx:2647-2672`), so this port's
///   endpoints always take the `myCache` / `thePoints2d` arm. Measured on
///   ATU01038: the `getLine` arm above already collapses the iso-edge tolerance
///   band (top edge tolerance now equals OCCT's 0.007844513808 to 1e-13).
/// - `needResolveUJump` / `needResolveVJump` (`cxx:1436-1451`, consumed at
///   `cxx:1504` / `cxx:1596-1598`) are PORTED: latched in the sample loop above
///   and passed to [`resolve_closed_surface_period_jump`]. Their
///   `IsUClosed(myPreci)` / `IsVClosed(myPreci)` gates are the ported
///   [`sa_is_u_closed`] / [`sa_is_v_closed`]
///   (`ShapeAnalysis_Surface.cxx:660-863`, `:865-1077`).
/// - `theC2D = interpolatePCurve(...)` (`cxx:771`) is ported; `approximatePCurve`
///   (`cxx:2220`) has no caller in this version.
pub fn project_curve_on_surface_perform(
    curve: &dyn Curve,
    surf: &dyn Surface,
    first: f64,
    last: f64,
    preci: f64,
    tol_first: f64,
    tol_last: f64,
    cache: &mut super::super::projection_cache::ProjectorCache,
) -> Option<Arc<dyn Curve2d>> {
    cache.done4 = false;
    if !first.is_finite() || !last.is_finite() || last - first < 1e-15 {
        return None;
    }
    if let Some(c2d) = project_analytic(curve, surf) {
        return Some(c2d);
    }
    let n = generate_curve_point_count(curve, first, last);
    let mut t_vals: Vec<f64> = (0..n)
        .map(|i| {
            if i == 0 {
                first
            } else if i + 1 == n {
                last
            } else {
                first + (last - first) * (i as f64) / ((n - 1) as f64)
            }
        })
        .collect();
    let mut pts3d: Vec<GpPnt> = t_vals.iter().map(|&t| curve.d0(t)).collect();

    // TEMP PC3D: the 3-D sample points this p-curve is interpolated through, so
    // they can be diffed against the OCCT probe's `C3VAL` dump.
    if std::env::var("PC3D").is_ok()
        && (first - 1.1414320265529101).abs() < 1e-12
        && (last - 1.5714265148303701).abs() < 1e-12
    {
        for (i, p) in pts3d.iter().enumerate() {
            eprintln!(
                "PC3D i={i} t={:.17e} p=({:.17e},{:.17e},{:.17e})",
                t_vals[i],
                p.x(),
                p.y(),
                p.z()
            );
        }
    }
    // TEMP SURFD2: surface `D0` / `D2` at the Newton seed and at the point the
    // p-curve sample should land on, for a bitwise comparison against the OCCT
    // probe's `OCCTSRFD2` dump.
    if let Ok(spec) = std::env::var("SURFD2") {
        if std::env::var("PC3D").is_ok()
            && (first - 1.1414320265529101).abs() < 1e-12
            && (last - 1.5714265148303701).abs() < 1e-12
        {
            for pair in spec.split(';') {
                let mut it = pair.split(',').filter_map(|s| s.trim().parse::<f64>().ok());
                if let (Some(u), Some(v)) = (it.next(), it.next()) {
                    let p0 = surf.value(u, v);
                    let (pd, ru, rv, ruu, rvv, ruv) = surf.d2(u, v);
                    eprintln!(
                        "SURFD2 uv=({u:.17e},{v:.17e}) d0=({:.17e},{:.17e},{:.17e}) \
d2p=({:.17e},{:.17e},{:.17e}) ru=({:.17e},{:.17e},{:.17e}) rv=({:.17e},{:.17e},{:.17e}) \
ruu=({:.17e},{:.17e},{:.17e}) rvv=({:.17e},{:.17e},{:.17e}) ruv=({:.17e},{:.17e},{:.17e})",
                        p0.x(), p0.y(), p0.z(),
                        pd.x(), pd.y(), pd.z(),
                        ru.x(), ru.y(), ru.z(),
                        rv.x(), rv.y(), rv.z(),
                        ruu.x(), ruu.y(), ruu.z(),
                        rvv.x(), rvv.y(), rvv.z(),
                        ruv.x(), ruv.y(), ruv.z(),
                    );
                }
            }
        }
    }

    // TEMP PCUV2: `ShapeAnalysis_Surface::ValueOfUV` for every sample, i.e. the
    // seed-independent projection, so it can be diffed point by point against
    // the OCCT probe's `OCCTUV` dump (no cache / previous-sample chaining).
    if std::env::var("PCUV2").is_ok()
        && (first - 1.1414320265529101).abs() < 1e-12
        && (last - 1.5714265148303701).abs() < 1e-12
    {
        for (i, p) in pts3d.iter().enumerate() {
            let q = value_of_uv(surf, p, preci);
            let d = surf.value(q.x(), q.y()).distance(p);
            eprintln!(
                "PCUV2 i={i} t={:.17e} preci={preci:.17e} uv=({:.17e},{:.17e}) d={d:.17e}",
                t_vals[i],
                q.x(),
                q.y()
            );
        }
    }

    // `approxPCurve` (`cxx:1106-1143`) tries the straight-pcurve shortcut
    // `getLine` before everything else (`cxx:1119`); a hit gives the exact
    // `Geom2d_Line` / 2-pole pcurve OCCT stores for an isoparametric edge and
    // refills `myCache` from its two endpoints (`cxx:1122-1142`).
    let mut gl = super::super::projection_cache::GetLineOut::default();
    let trace_gl = crate::debug_flags::iges_trace();
    let hit = super::super::projection_cache::line_pcurve_from_points(surf, &pts3d, &t_vals, preci, cache, &mut gl);
    if trace_gl {
        let a = &pts3d[0];
        let b = &pts3d[pts3d.len() - 1];
        eprintln!(
            "TRACE getLine kind={:?} n={} cache_n={} hit={} is_recompute={} from_cache={} lin2d={} np2={} curve=({},{},{})-({},{},{}) probes=({:.4},{:.4})({:.4},{:.4})({:.4},{:.4})({:.4},{:.4})",
            classify_surface_kind(surf),
            pts3d.len(),
            cache.len(),
            hit.is_some(),
            gl.is_recompute,
            gl.is_from_cache,
            hit.as_ref().map_or(false, |c| c.gp_lin2d().is_some()),
            hit.as_ref().map_or(0, |c| c.bspline_poles2d().map_or(0, |(p, _)| p.len())),
            a.x(), a.y(), a.z(),
            b.x(), b.y(), b.z(),
            gl.probes2d[0].x(), gl.probes2d[0].y(),
            gl.probes2d[1].x(), gl.probes2d[1].y(),
            gl.probes2d[2].x(), gl.probes2d[2].y(),
            gl.probes2d[3].x(), gl.probes2d[3].y(),
        );
    }
    if std::env::var("IGES_TRACE_GL").is_ok() && surf.gp_sphere().is_some() {
        eprintln!(
            "GLSPH is_recompute={} hit={} up={} vp={} uprd={:?} vprd={:?} tval0={:.17e} tvaln={:.17e} p0=({:.17e},{:.17e}) p1=({:.17e},{:.17e}) p2=({:.17e},{:.17e}) p3=({:.17e},{:.17e})",
            gl.is_recompute,
            hit.is_some(),
            surf.is_u_periodic(),
            surf.is_v_periodic(),
            surf.u_period(),
            surf.v_period(),
            t_vals[0],
            t_vals[t_vals.len() - 1],
            gl.probes2d[0].x(), gl.probes2d[0].y(),
            gl.probes2d[1].x(), gl.probes2d[1].y(),
            gl.probes2d[2].x(), gl.probes2d[2].y(),
            gl.probes2d[3].x(), gl.probes2d[3].y(),
        );
    }
    if let Some(c2d) = hit {
        let change_cycle = cache.change_cycle(&pts3d[0], &pts3d[pts3d.len() - 1]);
        cache.store(&pts3d, &[gl.first2d, gl.last2d], change_cycle);
        return Some(c2d);
    }

    // `approxPCurve` (`cxx:1147-1287`): the isoparametric test and its pre-block.
    // `cxx:1153-1157`: `isAnalytic` is false only for the exact dynamic types
    // `Geom_BezierSurface` / `Geom_BSplineSurface`. `DynamicType()`, not
    // `IsKind`, is what OCCT tests, so a `Geom_RectangularTrimmedSurface`
    // wrapping either patch stays analytic.
    let is_analytic = !(surf.is_bspline_surface() || surf.osculating_bspline().is_some());
    let (uf, ul) = surf.u_range();
    let (vf, vl) = surf.v_range();
    let mut iso_out = IsoParamOut::new(n);
    let iso_param = is_an_isoparametric(surf, &pts3d, &t_vals, &mut iso_out);
    let mut iso_value = 0.0f64;
    let mut iso_par1 = 0.0f64;
    let mut iso_par2 = 0.0f64;
    let mut value_p1 = iso_out.value_p1;
    let mut value_p2 = iso_out.value_p2;
    let mut iso_closed = false;
    if iso_param {
        // `cxx:1190-1209`.
        let (parf, parl) = if iso_out.is_type_u {
            iso_value = value_p1.x();
            iso_par1 = value_p1.y();
            iso_par2 = value_p2.y();
            iso_closed = sa_is_v_closed(surf, preci);
            (vf, vl)
        } else {
            iso_value = value_p1.y();
            iso_par1 = value_p1.x();
            iso_par2 = value_p2.x();
            iso_closed = sa_is_u_closed(surf, preci);
            (uf, ul)
        };
        if !iso_out.iso_par2d3d && !is_analytic {
            // `cxx:1210-1244`: only for a non-analytic surface whose samples
            // are not already the iso points. OCCT computes `Cf` / `Cl` from
            // `cIso` here but never uses them (`cxx:1211-1222` is dead), so the
            // port omits that read.
            let tdeb = iso_out.params_out[1];
            if iso_closed && (iso_par1 == parf || iso_par1 == parl) {
                // `cxx:1224-1244`.
                if (tdeb - parf).abs() < (tdeb - parl).abs() {
                    iso_par1 = parf;
                } else {
                    iso_par1 = parl;
                }
                if iso_out.is_type_u {
                    value_p1.set_y(iso_par1);
                } else {
                    value_p1.set_x(iso_par1);
                }
            }
            if iso_closed && (iso_par2 == parf || iso_par2 == parl) {
                // `cxx:1245-1265`.
                let tfin = iso_out.params_out[n - 2];
                if (tfin - parf).abs() < (tfin - parl).abs() {
                    iso_par2 = parf;
                } else {
                    iso_par2 = parl;
                }
                if iso_out.is_type_u {
                    value_p2.set_y(iso_par2);
                } else {
                    value_p2.set_x(iso_par2);
                }
            }
            if !iso_closed {
                // `cxx:1267-1285`: swap the ends when the outer end is farther
                // from its sample.
                let tfin = iso_out.params_out[n - 2];
                if (tdeb - iso_par1).abs() > (tdeb - iso_par2).abs()
                    && (tfin - iso_par2).abs() > (tfin - iso_par1).abs()
                {
                    std::mem::swap(&mut value_p1, &mut value_p2);
                    if iso_out.is_type_u {
                        iso_value = value_p1.x();
                        iso_par1 = value_p1.y();
                        iso_par2 = value_p2.y();
                    } else {
                        iso_value = value_p1.y();
                        iso_par1 = value_p1.x();
                        iso_par2 = value_p2.x();
                    }
                }
            }
        }
    }

    // `approxPCurve` (`ShapeConstruct_ProjectCurveOnSurface.cxx:1290-1465`): the
    // point-sampling loop. Every interior sample is projected by
    // `NextValueOfUV` walking from the previous parameter point (`cxx:1432`)
    // instead of by an independent global projection. The two endpoints keep the
    // projections `getLine` left in `thePoints2d` (`cxx:1388-1399`) unless
    // `getLine` reported a fixed period jump (`cxx:1401-1424`), in which case
    // they are re-projected from `myCache` or globally. The `isoParam` arms
    // (`cxx:1317-1375`) and the `p1OnIso` / `p2OnIso` endpoint overrides
    // (`cxx:1379-1385`) come from `isAnIsoparametric`, ported above.
    //
    // `needResolveUJump` / `needResolveVJump` (`cxx:1436-1451`) ARE ported here
    // and consumed by `resolve_closed_surface_period_jump` (`cxx:1504`,
    // `cxx:1596-1598`).
    //
    // `cxx:1291`: the running seed tolerance, `myPreci` before the first sample.
    let mut gap = preci;
    // `cxx:1290-1291`: `Up` / `Vp` are the unrestricted surface bounds'
    // widths (`cxx:1160-1161`).
    let up = ul - uf;
    let vp = vl - vf;
    // `cxx:1296-1304`: the sample order reverses when the incoming cache is
    // anchored on the LAST point of this curve.
    let change_cycle_val = cache.change_cycle(&pts3d[0], &pts3d[n - 1]);
    let mut is_from_cache = false;
    let mut saved_point = GpPnt2d::zero();
    let mut pts: Vec<GpPnt2d> = vec![GpPnt2d::zero(); n];
    // `cxx:1020-1021`: `getLine` wrote `thePoints2d(1)` / `thePoints2d(aNb)`
    // before it failed; the `!isRecompute` endpoint arm reads them back.
    pts[0] = gl.first2d;
    pts[n - 1] = gl.last2d;
    // `cxx:1466-1467`: the singularities of the surface. OCCT triggers
    // `ComputeSingularities` lazily from `NbSingularities(myPreci)` inside the
    // loop below (`ShapeAnalysis_Surface.cxx:305-309`) and again at `cxx:1469`;
    // it is a pure function of the surface and its bounds, so it is computed
    // once here and shared by the latch, `correct_degenerated_points` and
    // `resolve_closed_surface_period_jump`.
    let singularities = super::super::singularities::compute_singularities(surf);
    // `cxx:1308-1310` `p2d`: loop-carried. OCCT leaves it uninitialized, but
    // every first-iteration arm assigns it before use.
    let mut p2d = GpPnt2d::zero();
    // `cxx:1308-1311`: the lattice of the jump test below.
    let mut need_resolve_u_jump = false;
    let mut need_resolve_v_jump = false;
    let mut prev_p3d = GpPnt::zero();
    let mut prev_p2d = GpPnt2d::zero();
    for ii in 1..=n {
        // `cxx:1313`: the sample order may be reversed.
        let pnt_index = if change_cycle_val { n - ii + 1 } else { ii };
        let p3d = pts3d[pnt_index - 1];
        if iso_param {
            // `cxx:1317-1375`: the isoparametric arms.
            let mut t_par = 0.0f64;
            if iso_out.iso_par2d3d {
                // `cxx:1319-1328`.
                t_par = if iso_par2 > iso_par1 {
                    t_vals[pnt_index - 1]
                } else {
                    iso_out.t1 + iso_out.t2 - t_vals[pnt_index - 1]
                };
            } else if !is_analytic {
                // `cxx:1330-1346`: the projected iso parameters `pout`.
                t_par = if pnt_index == 1 {
                    iso_par1
                } else if pnt_index == n {
                    iso_par2
                } else {
                    iso_out.params_out[pnt_index - 1]
                };
            }
            if !iso_out.iso_par2d3d && is_analytic {
                // `cxx:1347-1362`: an analytic surface whose samples are not
                // already the iso points still projects by `NextValueOfUV`;
                // only the two ends take the iso values.
                if pnt_index == 1 {
                    p2d = value_p1;
                } else if pnt_index == n {
                    p2d = value_p2;
                } else {
                    let (q, g) = next_value_of_uv(
                        surf,
                        &p2d,
                        &p3d,
                        preci,
                        occt_core::precision::CONFUSION + 1000.0 * gap,
                    );
                    p2d = q;
                    gap = g;
                }
            } else if iso_out.is_type_u {
                // `cxx:1364-1376`.
                p2d = GpPnt2d::new(iso_value, t_par);
            } else {
                p2d = GpPnt2d::new(t_par, iso_value);
            }
        } else if pnt_index == 1 && iso_out.p1_on_iso {
            // `cxx:1379-1381`: endpoint override from `isAnIsoparametric`,
            // which can be set even when `isoParam` is false.
            p2d = value_p1;
        } else if pnt_index == n && iso_out.p2_on_iso {
            // `cxx:1382-1384`.
            p2d = value_p2;
        } else if pnt_index == 1 || pnt_index == n {
            if !gl.is_recompute {
                // `cxx:1388-1399`: no period jump among the `getLine` probes, so
                // the endpoint keeps the seeded-Newton projection there.
                p2d = pts[pnt_index - 1];
                gap = gl.last_gap;
                if pnt_index == 1 {
                    is_from_cache = gl.is_from_cache;
                    saved_point = p2d;
                }
                continue;
            }
            // `cxx:1401-1424`: a period jump was fixed, so re-project the
            // endpoint from `myCache` (walking Newton) or, without a cache hit,
            // globally.
            match cache.find(&p3d, preci * preci) {
                Some(cache_pnt) => {
                    let (q, g) = next_value_of_uv(
                        surf,
                        &cache_pnt.second,
                        &p3d,
                        preci,
                        occt_core::precision::CONFUSION + gap,
                    );
                    p2d = q;
                    gap = g;
                    if pnt_index == 1 {
                        is_from_cache = true;
                        saved_point = cache_pnt.second;
                    }
                }
                None => {
                    let (q, g) = value_of_uv_with_gap(surf, &p3d, preci);
                    p2d = q;
                    gap = g;
                }
            }
        } else {
            // `cxx:1432`: the interior sample walks from the previous one.
            let (q, g) = next_value_of_uv(
                surf,
                &p2d,
                &p3d,
                preci,
                occt_core::precision::CONFUSION + 1000.0 * gap,
            );
            p2d = q;
            gap = g;
        }
        pts[pnt_index - 1] = p2d;
        // TEMP PCUV: every projected sample, so the port's `NextValueOfUV`
        // chain can be diffed sample by sample against OCCT's interpolant
        // knots (a clamped interpolant passes through its samples).
        if std::env::var("PCUV").is_ok()
            && (first - 1.1414320265529101).abs() < 1e-12
            && (last - 1.5714265148303701).abs() < 1e-12
        {
            eprintln!(
                "PCUV i={} t={:.17e} uv=({:.17e},{:.17e}) gap={:.17e}",
                pnt_index - 1,
                t_vals[pnt_index - 1],
                p2d.x(),
                p2d.y(),
                gap
            );
        }
        // `cxx:1436-1451`: `needResolveUJump` / `needResolveVJump`. A sample
        // whose U (V) parameter moved a full period while its 3D point stayed
        // within `myPreci` of the previous one is a seam crossing the point-wise
        // projection resolved the wrong way; the flags widen the consumers below
        // (`cxx:1504`, `cxx:1596-1598`) from "the surface is closed" to "this
        // pcurve jumped a period". Only an interior sample of a longer sequence
        // (`theNbPnt > 23`, `ii > 2`, `ii < theNbPnt`) that lands on a
        // `Geom_BSplineSurface` with a singularity and is not already
        // U(V)-closed qualifies.
        if n > 23 && ii > 2 && ii < n {
            let is_bspline = crate::geom_bnd_lib_surface3d::surface_kind(surf)
                == crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType::BSplineSurface;
            if is_bspline
                && prev_p3d.distance(&p3d) < preci
                && super::super::singularities::nb_singularities(&singularities, preci) > 0
            {
                if (p2d.x() - prev_p2d.x()).abs() > 0.95 * up
                    && !sa_is_u_closed(surf, preci)
                {
                    need_resolve_u_jump = true;
                }
                if (p2d.y() - prev_p2d.y()).abs() > 0.95 * vp
                    && !sa_is_v_closed(surf, preci)
                {
                    need_resolve_v_jump = true;
                }
            }
        }
        // `cxx:1452-1453`.
        prev_p3d = p3d;
        prev_p2d = p2d;
        // `cxx:1455-1463`: the seed for the next sample is the linear
        // extrapolation of this sample against the previously stored one.
        if ii > 1 {
            let base = if change_cycle_val {
                pts[pnt_index]
            } else {
                pts[pnt_index - 2]
            };
            p2d = GpPnt2d::new(2.0 * p2d.x() - base.x(), 2.0 * p2d.y() - base.y());
        }
    }
    // `approxPCurve` (`cxx:1466-1470`): the singularities of the surface
    // (`ShapeAnalysis_Surface::ComputeSingularities`) feed both the
    // degenerate-point projection of the sampled pcurve and the singularity
    // loop that corrects its ends. The degenerate-point projection is guarded
    // by `if (!isoPar2d3d)` at
    // `ShapeConstruct_ProjectCurveOnSurface.cxx:1467`.
    if !iso_out.iso_par2d3d {
        super::super::singularities::correct_degenerated_points(
            surf,
            curve,
            &singularities,
            preci,
            tol_first,
            tol_last,
            &pts3d,
            &t_vals,
            &mut pts,
        );
    }
    // `approxPCurve` (`cxx:1499-1890`): before any curve is built, remove the
    // full-period jumps that the point-wise projection left in the sampled UV
    // sequence of a closed surface, then run "Handle AdjustOverDegen".
    cache.done4 = resolve_closed_surface_period_jump(
        surf,
        curve,
        &singularities,
        preci,
        need_resolve_u_jump,
        need_resolve_v_jump,
        is_from_cache,
        saved_point,
        cache.adjust_over_degen_off,
        &mut t_vals,
        &mut pts3d,
        &mut pts,
    );
    // `Perform:cxx:771`: `theC2D = interpolatePCurve(aNbPini, aPoints2d, aParams)`,
    // the single other source of `theC2D` besides `getLine`. A `None` result is
    // OCCT's null `theC2D` (`cxx:773`, `ShapeExtend_FAIL1`).
    let c2d = super::super::projection_cache::interpolate_pcurve(&pts, &t_vals, preci)?;
    // `cxx:1892-1903`: `myCache` is refilled from the pcurve just built, so the
    // next pcurve on this surface is projected from these endpoints.
    let change_cycle = cache.change_cycle(&pts3d[0], &pts3d[pts3d.len() - 1]);
    cache.store(&pts3d, &pts, change_cycle);
    Some(c2d)
}

pub fn project_curve_on_surface(
    s: &dyn Surface,
    curve: &dyn Curve,
    t_vals: &[f64],
) -> Option<Vec<GpPnt2d>> {
    let kind = classify_surface_kind(s);
    if kind != SurfaceKind::Torus {
        return None;
    }
    let proj = Projector::from_surface(s, kind)?;
    Some(proj.project_seq(curve, t_vals))
}

/// Recover the plane's location and in-plane axes from its surface `d0`.
pub(crate) fn plane_axes(surf: &dyn Surface) -> Option<(GpPnt, GpVec, GpVec)> {
    let o = surf.d0(0.0, 0.0);
    let x = GpVec::from_pnts(&o, &surf.d0(1.0, 0.0));
    let y = GpVec::from_pnts(&o, &surf.d0(0.0, 1.0));
    if x.square_magnitude() < 1e-24 || y.square_magnitude() < 1e-24 {
        return None;
    }
    Some((o, x, y))
}

/// Recover a cylinder's axis point, unit axes and radius from its surface `d0`.
pub(crate) fn cylinder_axes(surf: &dyn Surface) -> Option<(GpPnt, GpVec, GpVec, GpVec)> {
    let a = surf.d0(0.0, 0.0);
    let b = surf.d0(PI, 0.0);
    let c = surf.d0(FRAC_PI_2, 0.0);
    let o = sampled_mid(&a, &b);
    let r = 0.5 * a.distance(&b);
    if r < 1e-12 {
        return None;
    }
    let x = GpVec::new((a.x() - o.x()) / r, (a.y() - o.y()) / r, (a.z() - o.z()) / r);
    let y = GpVec::new((c.x() - o.x()) / r, (c.y() - o.y()) / r, (c.z() - o.z()) / r);
    let z = GpVec::from_pnts(&a, &surf.d0(0.0, 1.0));
    if z.square_magnitude() < 1e-24 {
        return None;
    }
    Some((o, x, y, z))
}

/// Orthonormal radial frame of a cone: `x` is the `u = 0` radial direction
/// (the generatrix direction projected perpendicular to the axis).
pub(crate) fn cone_frame(s: &dyn Surface, _apex: &GpPnt, z: &GpVec) -> Option<(GpVec, GpVec, GpVec)> {
    let dz = GpVec::from_pnts(&s.d0(0.0, 0.0), &s.d0(0.0, 1.0));
    let x = dz.subtracted(&z.multiplied_scalar(dz.dot(z))).normalized();
    if x.square_magnitude() < 1e-20 {
        return None;
    }
    let y = z.crossed(&x).normalized();
    Some((x, y, *z))
}

/// The cone's placement point (the v=0 ring plane, `ElSLib::ConeValue` base).
///
/// `s.d0(0,0)` lies on the v=0 ring: `placement + RefRadius·x`. Subtracting
/// its radial component about the axis (from the apex) recovers the placement
/// point on the axis.
pub(crate) fn cone_placement(s: &dyn Surface, apex: &GpPnt, z: &GpVec) -> Option<GpPnt> {
    let p0 = s.d0(0.0, 0.0);
    let d = GpVec::from_pnts(apex, &p0);
    let radial = d.subtracted(&z.multiplied_scalar(d.dot(z)));
    Some(GpPnt::from_xyz(&p0.coord.subtracted(&radial.coord)))
}

// ---------------------------------------------------------------------------
// Analytic isoparametric pcurves
// ---------------------------------------------------------------------------

/// `Adaptor3d_Curve::GetType() == GeomAbs_Line`: a `Geom_Line` (or a trimmed
/// one). This is the *concrete class* test of
/// `GeomAdaptor_Curve::load` (`GeomAdaptor_Curve.cxx:252-311`), which unwraps a
/// `Geom_TrimmedCurve` to its basis (the port's trimmed wrappers forward the
/// query) and then tests `gp_line()`. It replaces the earlier
/// "unbounded parameter range" heuristic, which also matched an offset curve.
pub(crate) fn is_line(curve: &dyn Curve) -> bool {
    curve.gp_line().is_some()
}

/// `Adaptor3d_Curve::GetType() == GeomAbs_Circle`: a `Geom_Circle` (or a
/// trimmed one). Deliberately NOT a geometric test - a `Geom_BSplineCurve`
/// whose image is a circle stays `GeomAbs_BSplineCurve` in OCCT
/// (`GeomAdaptor_Curve::load`, `GeomAdaptor_Curve.cxx:252-311`; the same order
/// as `edge_edge/find_solutions/curve_kind.rs`), and
/// `ProjLib_ProjectedCurve::Project` (`ProjLib_ProjectedCurve.cxx:242-270`)
/// only has analytic overloads for `GeomAbs_Line` (`:247`), `GeomAbs_Circle`
/// (`:250`), `GeomAbs_Ellipse` (`:253`), `GeomAbs_Hyperbola` (`:256`) and
/// `GeomAbs_Parabola` (`:259`): `GeomAbs_BSplineCurve` / `BezierCurve` /
/// `OffsetCurve` / `OtherCurve` break out (`:262-266`) into the general
/// approximation (`ProjLib_ComputeApprox`, `ProjLib_ProjectedCurve.cxx:720`).
pub(crate) fn is_circle(curve: &dyn Curve) -> bool {
    curve.gp_circ().is_some()
}

/// Fit a clamped B-spline through `pts` with the exact parameter range [a, b].
/// Degree 1 is a polyline through the sample points (endpoints interpolated).
pub(crate) fn bspline_from_samples(pts: &[GpPnt2d], a: f64, b: f64, degree: usize) -> Result<Geom2dBSplineCurve, String> {
    if pts.len() < 2 {
        return Err("pcurve_full: need at least 2 points".into());
    }
    if degree >= pts.len() {
        return Err("pcurve_full: degree must be below pole count".into());
    }
    let mut knots = occt_core::bspl::knots::build_uniform_knots(pts.len(), degree);
    for k in knots.iter_mut() {
        *k = a + (b - a) * *k;
    }
    let xs: Vec<f64> = pts.iter().map(|p| p.x()).collect();
    let ys: Vec<f64> = pts.iter().map(|p| p.y()).collect();
    Geom2dBSplineCurve::new(xs, ys, knots, degree).map_err(|e| e.to_string())
}

/// Build a `Geom2dLine` (unit-speed) or a degree-1 B-spline (otherwise) through
/// the first and last projected points over the parameter range `[a, b]`.
pub(crate) fn line_from_pts(pts: &[GpPnt2d], a: f64, b: f64) -> Option<Arc<dyn Curve2d>> {
    let q0 = pts[0];
    let q1 = pts[pts.len() - 1];
    let d = GpVec2d::new(q1.x() - q0.x(), q1.y() - q0.y());
    let mag = d.magnitude();
    if mag < 1e-15 {
        return None;
    }
    let dir = GpDir2d::from_vec2d(&d).ok()?;
    let len = b - a;
    if (mag - len).abs() <= 1e-6 * len.max(1.0) {
        let loc = GpPnt2d::new(q0.x() - a * dir.x, q0.y() - a * dir.y);
        Some(Arc::new(Geom2dLine::new(GpAx2d::new(loc, dir))))
    } else {
        let bs = bspline_from_samples(pts, a, b, 1).ok()?;
        Some(Arc::new(bs))
    }
}

/// Analytic pcurve of a line edge on an analytic revolution surface: a
/// generatrix projects to a `u`-isoline (varying `v`) or `v`-isoline.
/// `None` when the projection is not an isoparametric line (the caller then
/// samples instead).
pub(crate) fn iso_line_pcurve(curve: &dyn Curve, proj: &Projector, a: f64, b: f64) -> Option<Arc<dyn Curve2d>> {
    if b - a < 1e-15 {
        return None;
    }
    let t_vals = [a, 0.5 * (a + b), b];
    let pts = proj.project_seq(curve, &t_vals);
    let tol = 1e-6 * (1.0 + pts[0].x().abs() + pts[0].y().abs());
    let u_const = pts.iter().all(|q| (q.x() - pts[0].x()).abs() < tol);
    let v_const = pts.iter().all(|q| (q.y() - pts[0].y()).abs() < tol);
    if !u_const && !v_const {
        return None;
    }
    line_from_pts(&pts, a, b)
}

/// Analytic pcurve of a circle edge on an analytic revolution surface: a
/// latitude (constant `v`) or meridian (constant `u`) projects to a straight
/// UV line; a U-periodic full circle unwraps across the seam. `None` when the
/// projection is not an isoparametric line.
pub(crate) fn iso_circle_pcurve(curve: &dyn Curve, proj: &Projector, a: f64, b: f64) -> Option<Arc<dyn Curve2d>> {
    if b - a < 1e-15 {
        return None;
    }
    let n = 9;
    let t_vals: Vec<f64> = (0..n).map(|i| a + (b - a) * i as f64 / (n - 1) as f64).collect();
    let pts = proj.project_seq(curve, &t_vals);
    let tol = 1e-6 * (1.0 + pts[0].x().abs() + pts[0].y().abs());
    let u_const = pts.iter().all(|q| (q.x() - pts[0].x()).abs() < tol);
    let v_const = pts.iter().all(|q| (q.y() - pts[0].y()).abs() < tol);
    if u_const == v_const {
        // Degenerate (both) or non-isoparametric (neither): let the sampling
        // path handle it.
        return None;
    }
    line_from_pts(&pts, a, b)
}

// ---------------------------------------------------------------------------
// Main entry points
// ---------------------------------------------------------------------------

