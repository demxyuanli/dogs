//! `StepToGeom::MakeCurve2d` and `TranslateEdge::MakePCurve` with `GeomConvert_Units::DegreeToRadian` / `SetScale`: 2D primitives, pcurves and their unit conversion.
use super::*;
use occt_geom2d::geom2d_convert::{curve_to_bspline_curve_bspl, DEFAULT_PARAMETERISATION};

impl<'a> Resolver<'a> {
    /// 2D `CARTESIAN_POINT` → `GpPnt2d` (STEP pcurves use two-component points).
    pub(in crate::step) fn resolve_point_2d(&self, id: usize) -> Result<GpPnt2d, String> {
        let rec = self.record(id)?;
        if rec.type_name != "CARTESIAN_POINT" {
            return Err(format!("expected CARTESIAN_POINT at #{id}"));
        }
        let (x, y) = parse_xy(&rec.args[1]).map_err(|e| format!("CARTESIAN_POINT #{id}: {e}"))?;
        Ok(GpPnt2d::new(x, y))
    }

    /// 2D `DIRECTION` → `GpDir2d`.
    pub(in crate::step) fn resolve_direction_2d(&self, id: usize) -> Result<GpDir2d, String> {
        let rec = self.record(id)?;
        if rec.type_name != "DIRECTION" {
            return Err(format!("expected DIRECTION at #{id}"));
        }
        let (x, y) = parse_xy(&rec.args[1]).map_err(|e| format!("DIRECTION #{id}: {e}"))?;
        GpDir2d::new(x, y).map_err(|e| format!("DIRECTION #{id}: {e}"))
    }

    /// 2D `VECTOR` → `GpVec2d`.
    pub(in crate::step) fn resolve_vector_2d(&self, id: usize) -> Result<GpVec2d, String> {
        let rec = self.record(id)?;
        if rec.type_name != "VECTOR" {
            return Err(format!("expected VECTOR at #{id}"));
        }
        let dir_ref = parse_ref(&rec.args[1]).ok_or("VECTOR: bad direction ref")?;
        let mag = parse_f64(&rec.args[2])?;
        let dir = self.resolve_direction_2d(dir_ref)?;
        Ok(GpVec2d::new(dir.x * mag, dir.y * mag))
    }

    /// 2D `AXIS2_PLACEMENT_2D` → `GpAx22d` (the ref_direction is X; Y is its
    /// counter-clockwise normal).
    pub(in crate::step) fn resolve_axis22d(&self, id: usize) -> Result<GpAx22d, String> {
        let rec = self.record(id)?;
        if rec.type_name != "AXIS2_PLACEMENT_2D" {
            return Err(format!("expected AXIS2_PLACEMENT_2D at #{id}"));
        }
        let loc = self.resolve_point_2d(parse_ref(&rec.args[1]).ok_or("AXIS2_2D: bad loc ref")?)?;
        let vx = match parse_ref(&rec.args[2]) {
            Some(r) => self.resolve_direction_2d(r)?,
            None => GpDir2d::new(1.0, 0.0).unwrap(),
        };
        let vy = GpDir2d::new(-vx.y, vx.x).map_err(|e| format!("AXIS2_2D: {e}"))?;
        GpAx22d::new(loc, vx, vy).map_err(|e| format!("AXIS2_PLACEMENT_2D: {e}"))
    }

    /// 2D curve → `Arc<dyn Curve2d>` (STEP pcurve geometry). Source:
    /// `StepToGeom::MakeCurve2d`.
    pub(in crate::step) fn resolve_curve_2d(&self, id: usize) -> Result<(Arc<dyn Curve2d>, (f64, f64)), String> {
        if let Some(c) = self.curve2d_cache.borrow().get(&id) {
            return Ok(c.clone());
        }
        let rec = self.record(id)?;
        let (curve, range): (Arc<dyn Curve2d>, (f64, f64)) = match rec.type_name.as_str() {
            "LINE" => {
                let pnt = parse_ref(&rec.args[1]).ok_or("LINE: bad point ref")?;
                let vec = parse_ref(&rec.args[2]).ok_or("LINE: bad vector ref")?;
                let p = self.resolve_point_2d(pnt)?;
                let v = self.resolve_vector_2d(vec)?;
                let d = GpDir2d::new(v.x(), v.y()).map_err(|e| format!("LINE: {e}"))?;
                let mag = (v.x() * v.x() + v.y() * v.y()).sqrt();
                (Arc::new(Geom2dLine::new(GpAx2d::new(p, d))), (0.0, mag))
            }
            "CIRCLE" => {
                let ax = parse_ref(&rec.args[1]).ok_or("CIRCLE: bad axis ref")?;
                let r = parse_f64(&rec.args[2])?;
                let ax22 = self.resolve_axis22d(ax)?;
                let c: Arc<dyn Curve2d> = Arc::new(Geom2dCircle::new(GpCirc2d::new(ax22, r)));
                let range = (c.first_parameter(), c.last_parameter());
                (c, range)
            }
            "ELLIPSE" => {
                let ax = parse_ref(&rec.args[1]).ok_or("ELLIPSE: bad axis ref")?;
                // `StepToGeom::MakeEllipse2d` (`StepToGeom.cxx:1571-1598`): no length
                // factor in the 2-D path (`cxx:1583-1584`), and when
                // `majorR - minorR < 0` the X direction is mirrored in Y —
                // `gp_Dir2d(X.X(), -X.Y())` (`cxx:1591-1593`) — with the radii swapped.
                let maj = parse_f64(&rec.args[2])?;
                let min = parse_f64(&rec.args[3])?;
                let mut ax22 = self.resolve_axis22d(ax)?;
                let (maj, min) = if maj - min >= 0.0 {
                    (maj, min)
                } else {
                    let x = ax22.x_direction();
                    if let Ok(mirror) = GpDir2d::new(x.x(), -x.y()) {
                        ax22.set_x_direction(mirror);
                    }
                    (min, maj)
                };
                let c: Arc<dyn Curve2d> =
                    Arc::new(Geom2dEllipse::new(GpElips2d::new(ax22, maj, min)));
                let range = (c.first_parameter(), c.last_parameter());
                (c, range)
            }
            "B_SPLINE_CURVE_WITH_KNOTS" => {
                let degree = parse_f64(&rec.args[1])? as usize;
                let pts: Vec<GpPnt2d> = parse_ref_list(&rec.args[2])
                    .into_iter()
                    .map(|r| self.resolve_point_2d(r))
                    .collect::<Result<Vec<_>, _>>()?;
                // `StepToGeom::MakeBSplineCurve2d` (`StepToGeom.cxx:952-963`) runs
                // the same `MakeBSplineCurveCommon`: duplicate knots merged,
                // multiplicities clamped, poles trimmed, periodic descriptors and
                // closed curves made periodic. The rational complex arm
                // (`cxx:896-913`) keeps the weights in slot 3 of the merged record
                // (`transfer/entity_parse.rs::merge_complex_body`).
                let weights = bspline_weights(rec, pts.len())?;
                let knot_mults =
                    parse_usize_list(rec.args.get(6).map(|s| s.as_str()).unwrap_or("()"));
                let knot_values =
                    parse_real_list(rec.args.get(7).map(|s| s.as_str()).unwrap_or("()"));
                let (xs, ys): (Vec<f64>, Vec<f64>) =
                    pts.iter().map(|p| (p.x(), p.y())).unzip();
                let c: Arc<dyn Curve2d> = Arc::new(make_bspline_curve_2d_with_knots(
                    degree,
                    xs,
                    ys,
                    weights,
                    &knot_mults,
                    &knot_values,
                    curve_record_closed(rec),
                )?);
                let range = (c.first_parameter(), c.last_parameter());
                (c, range)
            }
            "POLYLINE" => {
                let pts: Vec<GpPnt2d> = parse_ref_list(&rec.args[1])
                    .into_iter()
                    .map(|r| self.resolve_point_2d(r))
                    .collect::<Result<Vec<_>, _>>()?;
                if pts.len() < 2 {
                    return Err("POLYLINE: need at least 2 points".into());
                }
                let knots = uniform_knots_for(pts.len(), 1);
                let (xs, ys): (Vec<f64>, Vec<f64>) =
                    pts.iter().map(|p| (p.x(), p.y())).unzip();
                let c: Arc<dyn Curve2d> = Arc::new(
                    Geom2dBSplineCurve::new(xs, ys, knots, 1)
                        .map_err(|e| format!("POLYLINE: {e}"))?,
                );
                let range = (c.first_parameter(), c.last_parameter());
                (c, range)
            }
            "B_SPLINE_CURVE" => {
                // No OCCT branch: `MakeCurve2d` sends every `StepGeom_BoundedCurve`
                // to `MakeBoundedCurve2d` (`StepToGeom.cxx:1414-1417`), and that
                // only maps the with-knots / rational complex (`cxx:477-486`),
                // `TRIMMED_CURVE` (`cxx:487-490`) and `POLYLINE` (`cxx:491-494`);
                // the knotless base class falls through to an empty handle
                // (`cxx:495-496`). `MakeBSplineCurveCommon` likewise requires a
                // `StepGeom_BSplineCurveWithKnots` downcast (`cxx:759-775`).
                // Synthesising knots here would be a rule OCCT does not have.
                return Err(format!(
                    "B_SPLINE_CURVE(2d) without knots has no OCCT branch (#{id})"
                ));
            }
            "TRIMMED_CURVE" => {
                // `StepToGeom::MakeTrimmedCurve2d` (`StepToGeom.cxx:2503-2563`).
                let basis_ref = parse_ref(&rec.args[1]).ok_or("TRIMMED_CURVE: bad basis ref")?;
                let basis_rec = self.record(basis_ref)?;
                let basis_kind = basis_rec.type_name.clone();
                let (basis, basis_range) = self.resolve_curve_2d(basis_ref)?;
                // `cxx:2514-2517`: OCCT tests the **constructed** basis
                // (`theGeomBasis->IsKind(Geom2d_BSplineCurve)`), not the STEP entity
                // name, so a `POLYLINE` basis (whose `MakePolyline2d` also yields a
                // `Geom2d_BSplineCurve`) is returned **untrimmed** as well.
                if basis.is_bspline2d() {
                    return Ok((basis, basis_range));
                }
                let (params1, points1) = parse_trimming_select(&rec.args[2]);
                let (params2, points2) = parse_trimming_select(&rec.args[3]);
                // `cxx:2523-2525`: both trims must be a single parameter select.
                // Trimming by Cartesian points is left unimplemented in OCCT itself
                // (`cxx:2501-2502`).
                if params1.len() != 1
                    || params2.len() != 1
                    || !points1.is_empty()
                    || !points2.is_empty()
                {
                    return Err(format!("TRIMMED_CURVE(2d): trims are not single parameters (#{id})"));
                }
                let (u1, u2) = (params1[0], params2[0]);
                // `cxx:2529-2550`: Line → `Dir()->Magnitude()` (**no** LengthFactor in
                // the 2-D path); Circle/Ellipse → `PlaneAngleFactor`, plus a π/2 shift
                // for an ellipse with `SemiAxis1 - SemiAxis2 < 0`; parabola/hyperbola
                // are left as a TODO in OCCT itself (`cxx:2547-2550`).
                let (fact, shift) = match basis_kind.as_str() {
                    "LINE" => {
                        let vec_ref = parse_ref(&basis_rec.args[2]).ok_or("LINE: bad vector ref")?;
                        let vec_rec = self.record(vec_ref)?;
                        (parse_f64(&vec_rec.args[2])?, 0.0)
                    }
                    "CIRCLE" | "ELLIPSE" => {
                        let mut shift = 0.0;
                        if basis_kind == "ELLIPSE" {
                            let s1 = parse_f64(&basis_rec.args[2])?;
                            let s2 = parse_f64(&basis_rec.args[3])?;
                            if s1 - s2 < 0.0 {
                                shift = 0.5 * PI;
                            }
                        }
                        (self.plane_angle_factor, shift)
                    }
                    _ => (1.0, 0.0),
                };
                let new_u1 = shift + u1 * fact;
                let new_u2 = shift + u2 * fact;
                // `cxx:2555-2560`: `Geom2d_TrimmedCurve(theGeomBasis, newU1, newU2,
                // SenseAgreement(theTrimmed))` and then
                // `Geom2dConvert::CurveToBSplineCurve(theTrimmed)` with the default
                // `Convert_TgtThetaOver2` parameterisation. The conversion is where
                // OCCT produces a `Geom2d_BSplineCurve` for the pcurve.
                let sense = parse_logical(rec.args.get(4).map(|s| s.as_str()), true);
                let trimmed =
                    Arc::new(Geom2dTrimmedCurve::new_sense(basis, new_u1, new_u2, sense, true));
                let converted =
                    curve_to_bspline_curve_bspl(trimmed.as_ref(), DEFAULT_PARAMETERISATION)
                        .ok_or_else(|| {
                            format!("TRIMMED_CURVE(2d): CurveToBSplineCurve failed (#{id})")
                        })?;
                let range = (converted.first_parameter(), converted.last_parameter());
                (Arc::new(converted), range)
            }
            other => {
                self.warn(format!("unsupported 2D curve entity {other} (#{id})"));
                return Err(format!("unsupported 2D curve entity {other} (#{id})"));
            }
        };
        if let Some(spec) = crate::debug_flags::value("PCRAW") {
            self.dump_pcraw(id, curve.as_ref(), range, &spec);
        }
        self.curve2d_cache.borrow_mut().insert(id, (curve.clone(), range));
        Ok((curve, range))
    }

    /// T-186 `PCRAW`: the definition of the 2-D curve this reader builds
    /// from one STEP entity, dumped **before** `XSAlgo_ShapeProcessor::CheckPCurve`
    /// can replace it with a reprojection. The OCCT counterpart is
    /// `occt_probe --pcraw <first> <tol>`, which runs `StepToGeom::MakeCurve2d`
    /// (`StepToGeom.cxx` → `MakeBSplineCurve2d` / `MakeTrimmedCurve2d`) on the
    /// same `PCURVE` entity. `spec` is the wanted first parameter (empty dumps
    /// every curve). See `specs/occt_probe/README.md` for the live probe list.
    fn dump_pcraw(&self, id: usize, curve: &dyn Curve2d, range: (f64, f64), spec: &str) {
        let spec = spec.trim();
        if !spec.is_empty() {
            let Ok(want) = spec.parse::<f64>() else {
                eprintln!("PCRAW: <first> is not a number: {spec}");
                return;
            };
            if (range.0 - want).abs() > 1e-7 {
                return;
            }
        }
        let kind = if curve.is_bspline2d() {
            "bspline2d"
        } else if curve.is_geom2d_line() {
            "line"
        } else if curve.is_geom2d_circle() {
            "circle"
        } else {
            "other"
        };
        let Some(b) = curve.bspline_copy2d() else {
            eprintln!("PCRAW src={id} type={kind} pf={:.17e} pl={:.17e}", range.0, range.1);
            return;
        };
        let (knots, mults) = b.distinct_knots_and_mults();
        eprintln!(
            "PCRAW src={id} type={kind} pf={:.17e} pl={:.17e} deg={} npoles={} nknots={} periodic={} rational={} nflat={}",
            range.0,
            range.1,
            b.degree,
            b.xs.len(),
            knots.len(),
            b.periodic as u8,
            b.weights.is_some() as u8,
            b.knots.len()
        );
        for i in 0..b.xs.len() {
            match b.weights.as_ref() {
                Some(w) => eprintln!("PCRAW  pole[{}]={:.17e},{:.17e} w={:.17e}", i + 1, b.xs[i], b.ys[i], w[i]),
                None => eprintln!("PCRAW  pole[{}]={:.17e},{:.17e}", i + 1, b.xs[i], b.ys[i]),
            }
        }
        for (j, (k, m)) in knots.iter().zip(mults.iter()).enumerate() {
            eprintln!("PCRAW  knot[{}]={k:.17e} m={m}", j + 1);
        }
        let flat: Vec<String> = b.knots.iter().map(|k| format!("{k:.17e}")).collect();
        eprintln!("PCRAW  flat={}", flat.join(","));
    }

    /// `PCURVE(name, basis_surface, reference_to_curve)` → the pcurve's basis
    /// surface reference and its 2D curve (the first item of the
    /// DEFINITIONAL_REPRESENTATION). Source: `StepToTopoDS_TranslateEdge::MakePCurve`.
    pub(in crate::step) fn resolve_pcurve(&self, id: usize) -> Result<(usize, Arc<dyn Curve2d>), String> {
        let rec = self.record(id)?;
        if rec.type_name != "PCURVE" {
            return Err(format!("expected PCURVE at #{id}"));
        }
        let basis_surf = parse_ref(&rec.args[1]).ok_or("PCURVE: bad basis surface ref")?;
        let dri = parse_ref(&rec.args[2]).ok_or("PCURVE: bad reference_to_curve ref")?;
        let dri_rec = self.record(dri)?;
        if dri_rec.type_name != "DEFINITIONAL_REPRESENTATION" {
            return Err(format!("expected DEFINITIONAL_REPRESENTATION at #{dri}"));
        }
        let items = parse_ref_list(&dri_rec.args[1]);
        let curve_ref = *items.first().ok_or("PCURVE: empty reference_to_curve")?;
        let (c2d, _range) = self.resolve_curve_2d(curve_ref)?;
        // `StepToTopoDS_TranslateEdge::MakePCurve` then
        // `GeomConvert_Units::DegreeToRadian` (`GeomConvert_Units.cxx:173-302`).
        let surf = self.resolve_surface(basis_surf)?;
        let c2d = self.degree_to_radian(c2d, curve_ref, surf.as_ref());
        Ok((basis_surf, c2d))
    }

    /// `GeomConvert_Units::DegreeToRadian`. STEP stores cylinder/cone/sphere/torus
    /// pcurve `U` in the file angle unit and cone `V` as axis length; OCCT's
    /// `ElSLib` cone uses `V` along the generatrix (`Length / cos(semiAngle)`).
    /// The length factor is the file's `StepData_Factors::LengthFactor()`
    /// (`StepToTopoDS_TranslateEdge.cxx:573-580` passes it to
    /// `GeomConvert_Units::DegreeToRadian`); the angle factor is
    /// `theLocalFactors.FactorDegreeRadian()`.
    ///
    /// Line location/direction are scaled then the direction is renormalized
    /// (`GeomConvert_Units.cxx:260-276`), so a V-iso generatrix keeps parameter
    /// speed 1 while a U-iso circle at written `V = height` moves to ElSLib `V`.
    /// BSpline poles are scaled in place (`cxx:278-294`). Plane Circle/Ellipse
    /// use `SetScale(origin, LengthFact)` (`cxx:215-222`). Circle/Ellipse on a
    /// revolution surface go through `Geom2dConvert::CurveToBSplineCurve` first
    /// (`cxx:229-236`) — that convert is not in this crate; those entities keep
    /// the raw 2D curve. Parabola/Hyperbola return unchanged (`cxx:237-252`).
    pub(in crate::step) fn degree_to_radian(
        &self,
        c2d: Arc<dyn Curve2d>,
        curve_ref: usize,
        surf: &dyn Surface,
    ) -> Arc<dyn Curve2d> {
        let length_fact = self.length_factor;
        let angle_fact = self.plane_angle_factor;
        let is_plane = classify_surface(surf) == SurfaceKind::Plane;
        // `GeomConvert_Units.cxx:191-227`. Offset (and other non-analytic
        // kinds) return the pcurve unchanged. Revolution is U-angle only.
        let (u_fact, v_fact) = if let Some((_, alpha)) = surf.cone_ref() {
            (angle_fact, length_fact / alpha.cos())
        } else if surf.is_surface_of_revolution() {
            (angle_fact, length_fact)
        } else {
            match classify_surface(surf) {
                SurfaceKind::Sphere | SurfaceKind::Torus => (angle_fact, angle_fact),
                SurfaceKind::Cylinder => (angle_fact, length_fact),
                SurfaceKind::Plane => (length_fact, length_fact),
                SurfaceKind::Cone | SurfaceKind::Other => return c2d,
            }
        };
        let Ok(rec) = self.record(curve_ref) else {
            return c2d;
        };
        match rec.type_name.as_str() {
            "CIRCLE" | "ELLIPSE" if is_plane => {
                let mut t = occt_core::gp::GpTrsf2d::default();
                if t.set_scale(&GpPnt2d::new(0.0, 0.0), length_fact).is_err() {
                    return c2d;
                }
                let mut scaled = c2d.clone_dyn();
                scaled.transform(&t);
                Arc::from(scaled)
            }
            "PARABOLA" | "HYPERBOLA" => c2d,
            _ if (u_fact - 1.0).abs() <= 1e-15 && (v_fact - 1.0).abs() <= 1e-15 => c2d,
            "LINE" => {
                let loc = c2d.d0(0.0);
                let new_loc = GpPnt2d::new(loc.x() * u_fact, loc.y() * v_fact);
                let (_, tan) = c2d.d1(0.0);
                let Ok(new_dir) = GpDir2d::new(tan.x() * u_fact, tan.y() * v_fact) else {
                    return c2d;
                };
                Arc::new(Geom2dLine::from_pnt_dir(new_loc, new_dir))
            }
            "B_SPLINE_CURVE_WITH_KNOTS" | "B_SPLINE_CURVE" | "POLYLINE" => {
                match scale_pcurve_bspline_poles(c2d.as_ref(), u_fact, v_fact) {
                    Some(scaled) => scaled,
                    None => {
                        // `GeomConvert_Units.cxx:278-293` only ever runs
                        // `SetPole` on an existing `Geom2d_BSplineCurve`, so a
                        // failure here is a port defect, not an OCCT branch.
                        eprintln!("DegreeToRadian: BSpline pole affinity failed, curve kept unscaled");
                        c2d
                    }
                }
            }
            _ => c2d,
        }
    }
}

/// `GeomConvert_Units.cxx:278-294` — affinity on the poles of the already
/// constructed `Geom2d_BSplineCurve` (`aBSpline2d->SetPole(i, ...)`): knots,
/// multiplicities, weights and periodicity are kept. The construction cannot
/// fail for a curve that already is a B-spline, so `None` reports the port
/// defect described at the call site.
fn scale_pcurve_bspline_poles(
    c2d: &dyn Curve2d,
    u_fact: f64,
    v_fact: f64,
) -> Option<Arc<dyn Curve2d>> {
    let b = c2d.bspline_copy2d()?;
    let xs: Vec<f64> = b.xs.iter().map(|x| x * u_fact).collect();
    let ys: Vec<f64> = b.ys.iter().map(|y| y * v_fact).collect();
    Geom2dBSplineCurve::from_flat(xs, ys, b.weights, b.knots, b.degree, b.periodic)
        .ok()
        .map(|c| Arc::new(c) as Arc<dyn Curve2d>)
}
