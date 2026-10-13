//! 2-D (UV) curve transfer: the p-curve is corrected, transformed and lifted to the plane Z=0. Source: `BRepToIGES_BRWire::TransferEdge` (`BRepToIGES_BRWire.cxx:340-585`) with `Geom2dToIGES_Geom2dCurve` and `ShapeBuild_Edge::TransformPCurve`.
use super::prelude::*;

// ---------------------------------------------------------------------------
// 2-D (UV) curve transfer
//
// 'BRepToIGES_BRWire::TransferEdge(edge, face, originMap, length, false)'
// ('BRepToIGES_BRWire.cxx:340-585') writes the edge's p-curve as an IGES curve
// entity: the p-curve is corrected per surface type, transformed by
// 'ShapeBuild_Edge::TransformPCurve', reversed for a REVERSED edge, and handed
// to 'Geom2dToIGES_Geom2dCurve::Transfer2dCurve', which lifts it into the plane
// Z=0 and reuses the 3-D 'GeomToIGES_GeomCurve::TransferCurve' dispatch.
//
// The port has no 'write.surfacecurve.mode' switch, so 'GetPCurveMode()' keeps
// OCCT's default 1 (On, 'Interface_StaticStandards.cxx:85'); every call here
// passes 'theIsBRepMode = false' ('BRepToIGES_BRShell.cxx:280', ':312', ':346'),
// which makes 'analyticMode = (GetConvertSurfaceMode() == 0 && false)' false
// ('BRepToIGES_BRWire.cxx:356').
// ---------------------------------------------------------------------------

/// 'gp_Dir2d::Transform' ('gp_Dir2d.cxx:80-105'). The port's
/// 'GpTrsf2d::transforms_xy' handles points; directions follow OCCT's own rule
/// (raw matrix, normalise, reverse for a negative scale factor).
pub(super) fn dir2d_transform(t: &GpTrsf2d, d: &GpDir2d) -> GpDir2d {
    match t.form() {
        TrsfForm::Identity | TrsfForm::Translation => *d,
        TrsfForm::PntMirror => GpDir2d::new(-d.x, -d.y).unwrap_or(*d),
        TrsfForm::Scale => {
            if t.scale_factor() < 0.0 {
                GpDir2d::new(-d.x, -d.y).unwrap_or(*d)
            } else {
                *d
            }
        }
        _ => {
            let mut xy = GpXY::new(d.x, d.y);
            xy.multiply_mat2d(t.h_vectorial_part());
            let _ = xy.normalize();
            if t.scale_factor() < 0.0 {
                xy.reverse();
            }
            GpDir2d::new(xy.x, xy.y).unwrap_or(*d)
        }
    }
}

/// 'gp_Pnt2d::Transform' ('gp_Pnt2d.hxx'), which is exactly what
/// 'GpTrsf2d::transforms_xy' applies.
pub(super) fn pnt2d_transform(t: &GpTrsf2d, p: &GpPnt2d) -> GpPnt2d {
    let mut c = p.coord;
    t.transforms_xy(&mut c);
    GpPnt2d::from_xy(c)
}

/// 'gp_Ax22d::Transform' ('gp_Ax22d.hxx:360-367'): location, X and Y directions
/// are transformed independently (a mirroring transform may leave the axis
/// indirect, exactly as OCCT's flag-free 'gp_Ax22d' does).
pub(super) fn ax22d_transform(t: &GpTrsf2d, a: &GpAx22d) -> GpAx22d {
    GpAx22d {
        point: pnt2d_transform(t, &a.point),
        vxdir: dir2d_transform(t, &a.vxdir),
        vydir: dir2d_transform(t, &a.vydir),
    }
}

/// 'gp_Trsf2d::SetMirror(const gp_Ax2d&)' ('gp_Trsf2d.cxx:31-46').
///
/// The port's 'GpTrsf2d::set_mirror_ax2d' stores the *positive* reflection
/// matrix together with 'scale = -1', so 'transforms_xy' applies '-R'; this
/// builds OCCT's own matrix (the negative reflection, 'scale = -1') so the
/// composite is the reflection R.
pub(super) fn mirror_ax2d_trsf(a: &GpAx2d) -> GpTrsf2d {
    let (vx, vy) = (a.direction().x, a.direction().y);
    let (x0, y0) = (a.location().x(), a.location().y());
    let matrix = GpMat2d::new(
        1.0 - 2.0 * vx * vx,
        -2.0 * vx * vy,
        -2.0 * vx * vy,
        1.0 - 2.0 * vy * vy,
    );
    let loc = GpXY::new(
        -2.0 * ((vx * vx - 1.0) * x0 + vx * vy * y0),
        -2.0 * (vx * vy * x0 + (vy * vy - 1.0) * y0),
    );
    GpTrsf2d {
        scale: -1.0,
        shape: TrsfForm::Ax1Mirror,
        matrix,
        loc,
    }
}

/// 'Mirror(gp_Ax2d(gp::Origin2d(), gp::Dir2d(1., 1.)))'
/// ('BRepToIGES_BRWire.cxx:463', ':471').
pub(super) fn mirror_origin_dir11() -> GpTrsf2d {
    let ax = GpAx2d::new(
        GpPnt2d::zero(),
        GpDir2d::new(1.0, 1.0).expect("gp::Dir2d(1,1) is non-null"),
    );
    mirror_ax2d_trsf(&ax)
}

/// 'Mirror(gp::OX2d())' ('BRepToIGES_BRWire.cxx:464', ':472').
pub(super) fn mirror_ox2d() -> GpTrsf2d {
    let ax = GpAx2d::new(
        GpPnt2d::zero(),
        GpDir2d::new(1.0, 0.0).expect("gp::DX2d is non-null"),
    );
    mirror_ax2d_trsf(&ax)
}

/// 'Curve2d->Translate(gp_Vec2d(dx, dy))'.
pub(super) fn translate2d(dx: f64, dy: f64) -> GpTrsf2d {
    let mut t = GpTrsf2d::identity();
    t.set_translation_vec(&GpVec2d::new(dx, dy));
    t
}

/// 'Geom2d_Curve::TransformedParameter' ('Geom2d_Curve.cxx:41-44' default;
/// 'Geom2d_Line.cxx:246-253', 'Geom2d_Parabola.cxx:249-256',
/// 'Geom2d_TrimmedCurve.cxx:297-300', 'Geom2d_OffsetCurve.cxx:423-426'): only a
/// line, a parabola and the trimmed/offset wrappers over them rescale the
/// parameter; every other type keeps it.
pub(super) fn transformed_parameter(curve: &dyn Curve2d, u: f64, t: &GpTrsf2d) -> f64 {
    if curve.gp_lin2d().is_some() || curve.gp_parab2d().is_some() {
        if occt_core::precision::Precision::is_infinite(u) {
            u
        } else {
            u * t.scale_factor().abs()
        }
    } else if let Some(b) = curve.trimmed_basis() {
        transformed_parameter(b, u, t)
    } else if let Some(b) = curve.offset_basis() {
        transformed_parameter(b, u, t)
    } else {
        u
    }
}

/// 'Geom2d_*::Transform' for the concrete types this port stores as p-curves -
/// 'Geom2d_Line', 'Geom2d_Circle', 'Geom2d_Ellipse', 'Geom2d_Hyperbola',
/// 'Geom2d_Parabola', 'Geom2d_BSplineCurve', 'Geom2d_TrimmedCurve'.
///
/// The port's 'Curve2d::transform' only moves the *location* of a line / circle /
/// conic ('GpLin2d::transform' etc. leave the axis directions untouched), so the
/// transform is rebuilt from the exact-type queries instead of delegating to it.
/// Any other curve (Bezier / Offset) keeps the port's own 'transform'.
pub(super) fn curve2d_transformed(curve: &Arc<dyn Curve2d>, t: &GpTrsf2d) -> Arc<dyn Curve2d> {
    // 'Geom2d_TrimmedCurve::Transform' ('Geom2d_TrimmedCurve.cxx:287-293'):
    // transform the basis, then set the trim to the transformed parameters.
    if let Some(b) = curve.trimmed_basis() {
        let nb = curve2d_transformed(&Arc::from(b.clone_dyn()), t);
        let u1 = transformed_parameter(b, curve.first_parameter(), t);
        let u2 = transformed_parameter(b, curve.last_parameter(), t);
        return Arc::new(Geom2dTrimmedCurve::new_sense(nb, u1, u2, true, false));
    }
    if let Some(l) = curve.gp_lin2d() {
        // 'Geom2d_Line::Transform' ('Geom2d_Line.cxx:239-242') ->
        // 'gp_Lin2d::Position().Transform' (location and direction).
        return Arc::new(Geom2dLine::new(GpAx2d::new(
            pnt2d_transform(t, &l.location()),
            dir2d_transform(t, l.direction()),
        )));
    }
    if let Some(c) = curve.gp_circ2d() {
        // 'gp_Circ2d::Transform': axis transformed, 'radius *= |ScaleFactor|'.
        return Arc::new(Geom2dCircle::new(GpCirc2d::new(
            ax22d_transform(t, c.position()),
            c.radius() * t.scale_factor().abs(),
        )));
    }
    if let Some(e) = curve.gp_elips2d() {
        // 'gp_Elips2d::Transform': axis transformed, both radii scaled.
        let s = t.scale_factor().abs();
        return Arc::new(Geom2dEllipse::new(GpElips2d::new(
            ax22d_transform(t, e.axis()),
            e.major_radius * s,
            e.minor_radius * s,
        )));
    }
    if let Some(h) = curve.gp_hypr2d() {
        let s = t.scale_factor().abs();
        return Arc::new(Geom2dHyperbola::new(GpHypr2d::new(
            ax22d_transform(t, h.axis()),
            h.major_radius * s,
            h.minor_radius * s,
        )));
    }
    if let Some(p) = curve.gp_parab2d() {
        return Arc::new(Geom2dParabola::new(GpParab2d::new(
            ax22d_transform(t, p.axis()),
            p.focal * t.scale_factor().abs(),
        )));
    }
    if let (Some((xs, ys)), Some(knots), Some(degree)) = (
        curve.bspline_poles2d(),
        curve.bspline_knots2d(),
        curve.bspline_degree(),
    ) {
        // 'Geom2d_BSplineCurve::Transform': the poles carry the affine map. The
        // port's 2-D B-spline is non-rational (no weights).
        let poles: Vec<GpPnt2d> = xs
            .iter()
            .zip(ys.iter())
            .map(|(x, y)| pnt2d_transform(t, &GpPnt2d::new(*x, *y)))
            .collect();
        let (out_xs, out_ys): (Vec<f64>, Vec<f64>) =
            poles.iter().map(|p| (p.x(), p.y())).unzip();
        return Arc::new(Geom2dBSplineCurve {
            xs: out_xs,
            ys: out_ys,
            weights: None,
            knots: knots.to_vec(),
            degree,
            periodic: curve.is_periodic(),
        });
    }
    Arc::from(curve.transformed(t))
}

/// 'ShapeBuild_Edge::TransformPCurve(pcurve, trans, uFact, aFirst, aLast)'
/// ('ShapeBuild_Edge.cxx:596-699').
///
/// Returns 'None' for the arms this port cannot reproduce - the caller then
/// writes no UV curve for that edge (the null handle of the OCCT transfer).
pub(super) fn transform_pcurve(
    pcurve: &Arc<dyn Curve2d>,
    trans: &GpTrsf2d,
    u_fact: f64,
    first: &mut f64,
    last: &mut f64,
) -> Option<Arc<dyn Curve2d>> {
    // 'cxx:602-608': a non-identity trans is applied and the range ends are
    // mapped through 'TransformedParameter'.
    let mut result: Arc<dyn Curve2d> = if trans.form() != TrsfForm::Identity {
        let r = curve2d_transformed(pcurve, trans);
        *first = transformed_parameter(r.as_ref(), *first, trans);
        *last = transformed_parameter(r.as_ref(), *last, trans);
        r
    } else {
        Arc::from(pcurve.clone_dyn())
    };
    if u_fact == 1.0 {
        return Some(result); // 'cxx:609-612'
    }
    // 'cxx:614-618': a trimmed curve is replaced by its basis.
    if let Some(b) = result.trimmed_basis() {
        result = Arc::from(b.clone_dyn());
    }
    // 'tMatu.SetAffinity(gp::OY2d(), uFact)' ('cxx:620-621') scales the X (U)
    // coordinate by uFact ('gp_GTrsf2d::SetAffinity', 'gp_GTrsf2d.cxx:24-38',
    // with the Y axis: 'matrix = (uFact, 0; 0, 1)').
    if result.gp_lin2d().is_some() {
        // 'cxx:624-641': scale the two range points, rebuild the line through
        // them and take the range from 'ElCLib::Parameter'.
        let pf = result.d0(*first);
        let pl = result.d0(*last);
        let pf = GpPnt2d::new(pf.x() * u_fact, pf.y());
        let pl = GpPnt2d::new(pl.x() * u_fact, pl.y());
        let v = GpVec2d::new(pl.x() - pf.x(), pl.y() - pf.y());
        let dir = GpDir2d::from_vec2d(&v).ok()?;
        let line = GpLin2d::from_pnt_dir(pf, dir);
        *first = occt_core::elib::clib2d::parameter_lin2d(&line, &pf);
        *last = occt_core::elib::clib2d::parameter_lin2d(&line, &pl);
        return Some(Arc::new(Geom2dLine::new(GpAx2d::new(pf, dir))));
    }
    // 'cxx:642-656': a 2-D Bezier scales every pole on X and keeps its own
    // range; the weights (rational Bezier) are carried over unchanged.
    if result.is_bezier2d() {
        let poles: Vec<GpPnt2d> = result
            .poles2d()?
            .iter()
            .map(|p| GpPnt2d::new(p.x() * u_fact, p.y()))
            .collect();
        let bez = match result.bezier_weights2d() {
            Some(w) => Geom2dBezierCurve::rational(poles, w.to_vec()),
            None => Geom2dBezierCurve::new(poles),
        }
        .ok()?;
        return Some(Arc::new(bez));
    }
    let is_conic = result.gp_circ2d().is_some()
        || result.gp_elips2d().is_some()
        || result.gp_hypr2d().is_some()
        || result.gp_parab2d().is_some();
    // 'cxx:658-698': the conic arm trims the curve to 'aFirst'/'aLast' and
    // re-runs it through 'Geom2dConvert_ApproxCurve' (falling back to
    // 'Geom2dConvert::CurveToBSplineCurve(thecurve, Convert_QuasiAngular)');
    // any remaining 2-D type goes through 'CurveToBSplineCurve(result,
    // Convert_QuasiAngular)'; a B-spline is used as it is. All three then have
    // their poles scaled on X.
    let mut bs = if is_conic {
        let tcurve: Arc<dyn Curve2d> =
            Arc::new(Geom2dTrimmedCurve::new(result.clone(), *first, *last));
        let approx = occt_geom2d::Geom2dConvertApproxCurve::new(
            tcurve.as_ref(),
            occt_core::precision::APPROXIMATION,
            occt_core::kernel::geomabs::Shape::C1,
            100,
            6,
        );
        let bs = match approx.curve() {
            Some(curve) => curve.clone(),
            None => occt_geom2d::geom2d_convert::curve_to_bspline_curve_bspl(
                tcurve.as_ref(),
                occt_core::convert::ParameterisationType::QuasiAngular,
            )?,
        };
        // 'cxx:673-674': 'aFirst'/'aLast' follow the approximation's range.
        *first = bs.first_parameter();
        *last = bs.last_parameter();
        bs
    } else if !result.is_bspline2d() {
        occt_geom2d::geom2d_convert::curve_to_bspline_curve_bspl(
            result.as_ref(),
            occt_core::convert::ParameterisationType::QuasiAngular,
        )?
    } else {
        result.bspline_copy2d()?
    };
    // 'cxx:688-698': transform the poles. 'aFirst'/'aLast' are only rewritten in
    // the conic arm above, so a plain B-spline keeps its range.
    bs.xs = bs.xs.iter().map(|x| x * u_fact).collect();
    Some(Arc::new(bs))
}

/// 'Adaptor3d_CurveOnSurface.cxx:72-78' (to3d of a 'gp_Ax22d' on
/// 'gp_Pln(0,0,1,0)'): 'gp_Ax2(P, VX.Crossed(VY), VX)' with 'P = (x, y, 0)'.
pub(super) fn promote_ax22d(a: &GpAx22d) -> Option<GpAx2> {
    let p = GpPnt::new(a.point.x(), a.point.y(), 0.0);
    let vx = GpDir::from_xyz(&GpXyz::new(a.vxdir.x, a.vxdir.y, 0.0)).ok()?;
    let vy = GpDir::from_xyz(&GpXyz::new(a.vydir.x, a.vydir.y, 0.0)).ok()?;
    let n = vx.crossed(&vy).ok()?;
    Some(GpAx2::new(p, n, vx).unwrap_or_else(|_| GpAx2::from_axis(p, n)))
}

/// 'Geom2dToIGES_Geom2dCurve::Transfer2dCurve' ('Geom2dToIGES_Geom2dCurve.cxx:52-68'):
/// 'GC.TransferCurve(GeomAPI::To3d(start, gp_Pln(0,0,1,0)), Udeb, Ufin)'.
///
/// 'GeomAPI::To3d' ('GeomAPI.cxx:56-65') puts the 2-D curve in an
/// 'Adaptor3d_CurveOnSurface' over the plane Z=0 and 'GeomAdaptor::MakeCurve'
/// ('GeomAdaptor.cxx:45-92') builds the concrete 'Geom_*' curve of the same
/// type. 'Geom2dAdaptor_Curve::load' unwraps a 'Geom2d_TrimmedCurve' and keeps
/// the basis, and 'GeomToIGES_GeomCurve::TransferCurve(Geom_TrimmedCurve)'
/// ('GeomToIGES_GeomCurve.cxx:456-477') transfers that basis with the requested
/// range - so promoting the basis with '[First, Last]' is equivalent.
pub(super) fn promote_curve2d(curve: &dyn Curve2d) -> Option<Box<dyn Curve>> {
    let basis: &dyn Curve2d = curve.trimmed_basis().unwrap_or(curve);
    if let Some(l) = basis.gp_lin2d() {
        let p = GpPnt::new(l.location().x(), l.location().y(), 0.0);
        let d = GpDir::from_xyz(&GpXyz::new(l.direction().x, l.direction().y, 0.0)).ok()?;
        return Some(Box::new(GeomLine::new(GpLin::from_pnt_dir(p, d))));
    }
    if let Some(c) = basis.gp_circ2d() {
        return Some(Box::new(GeomCircle::new(GpCirc::new(
            promote_ax22d(c.position())?,
            c.radius(),
        ))));
    }
    if let Some(e) = basis.gp_elips2d() {
        return Some(Box::new(GeomEllipse::new(GpElips::new(
            promote_ax22d(e.axis())?,
            e.major_radius,
            e.minor_radius,
        ))));
    }
    if let Some(h) = basis.gp_hypr2d() {
        return Some(Box::new(GeomHyperbola::new(GpHypr::new(
            promote_ax22d(h.axis())?,
            h.major_radius,
            h.minor_radius,
        ))));
    }
    if let Some(p) = basis.gp_parab2d() {
        return Some(Box::new(GeomParabola::new(GpParab::new(
            promote_ax22d(p.axis())?,
            p.focal,
        ))));
    }
    if let (Some((xs, ys)), Some(knots), Some(degree)) = (
        basis.bspline_poles2d(),
        basis.bspline_knots2d(),
        basis.bspline_degree(),
    ) {
        // `Adaptor3d_CurveOnSurface::BSpline()` (`cxx:1489-1520`) lifts every
        // pole with `to3d` (`cxx:57-60`): `ElSLib::Value` on the plane, not
        // `(x, y, 0)`. `GeomAPI::To3d` uses `gp_Pln(0, 0, 1, 0)`
        // (`Geom2dToIGES_Geom2dCurve.cxx:52-68`). The port's 2-D B-spline is
        // non-rational, so the `IsRational()` arm (weights array) has no
        // ported source.
        let plane = occt_core::gp::GpPln::from_coefficients(0.0, 0.0, 1.0, 0.0).ok()?;
        let poles: Vec<GpPnt> = xs
            .iter()
            .zip(ys.iter())
            .map(|(x, y)| occt_core::elib::slib::plane_value(&plane, *x, *y))
            .collect();
        return Some(Box::new(GeomBSplineCurve {
            poles,
            weights: None,
            knots: knots.to_vec(),
            degree,
            periodic: basis.is_periodic(),
        }));
    }
    // `Adaptor3d_CurveOnSurface::Bezier()` (`cxx:1458-1485`) lifts the 2-D
    // Bezier's poles to Z=0 and `GeomAdaptor::MakeCurve` (`GeomAdaptor.cxx:60-77`)
    // builds the matching `Geom_BezierCurve`.
    //
    // UNPORTED (`cxx:1474-1481`): the rational sub-arm
    // (`Bez2d->IsRational()` -> `Geom_BezierCurve(Poles, Weights)`) has no port
    // equivalent - `occt_geom::GeomBezierCurve` is non-rational, the same gap
    // noted in `occt_geom::convert_bspl::curve_to_bspline_curve` (`cxx:313-321`).
    if basis.is_bezier2d() {
        let plane = occt_core::gp::GpPln::from_coefficients(0.0, 0.0, 1.0, 0.0).ok()?;
        let poles: Vec<GpPnt> = basis
            .poles2d()?
            .iter()
            .map(|p| occt_core::elib::slib::plane_value(&plane, p.x(), p.y()))
            .collect();
        return Some(Box::new(GeomBezierCurve::new(poles).ok()?));
    }
    // UNPORTED: any other type has no 'GeomAdaptor::MakeCurve' arm either
    // ('GeomAdaptor.cxx:79-80' throws OtherCurve).
    None
}

