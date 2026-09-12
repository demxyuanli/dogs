use super::prelude::*;
use super::*;

impl<'a> Resolver<'a> {

    pub(super) fn resolve_surface(&self, id: usize) -> Result<Arc<dyn Surface>, String> {
        if let Some(s) = self.surface_cache.borrow().get(&id) {
            return Ok(s.clone());
        }
        let rec = self.record(id)?;
        let surface: Arc<dyn Surface> = match rec.type_name.as_str() {
            "PLANE" => {
                let ax = parse_ref(&rec.args[1]).ok_or("PLANE: bad axis ref")?;
                Arc::new(GeomPlane::new(GpPln::new(self.resolve_axis2(ax)?.to_ax3())))
            }
            "CYLINDRICAL_SURFACE" => {
                let ax = parse_ref(&rec.args[1]).ok_or("CYLINDRICAL_SURFACE: bad axis ref")?;
                let r = parse_f64(&rec.args[2])?;
                Arc::new(GeomCylinder::new(
                    GpCylinder::new(self.resolve_axis2(ax)?.to_ax3(), r)
                        .map_err(|e| format!("CYLINDRICAL_SURFACE: {e}"))?,
                ))
            }
            "CONICAL_SURFACE" => {
                let ax = parse_ref(&rec.args[1]).ok_or("CONICAL_SURFACE: bad axis ref")?;
                let r = parse_f64(&rec.args[2])?;
                let a = parse_f64(&rec.args[3])?;
                Arc::new(GeomCone::new(
                    GpCone::new(self.resolve_axis2(ax)?.to_ax3(), r, a)
                        .map_err(|e| format!("CONICAL_SURFACE: {e}"))?,
                ))
            }
            "SPHERICAL_SURFACE" => {
                let ax = parse_ref(&rec.args[1]).ok_or("SPHERICAL_SURFACE: bad axis ref")?;
                let r = parse_f64(&rec.args[2])?;
                Arc::new(GeomSphere::new(
                    GpSphere::new(self.resolve_axis2(ax)?.to_ax3(), r)
                        .map_err(|e| format!("SPHERICAL_SURFACE: {e}"))?,
                ))
            }
            "TOROIDAL_SURFACE" => {
                let ax = parse_ref(&rec.args[1]).ok_or("TOROIDAL_SURFACE: bad axis ref")?;
                let maj = parse_f64(&rec.args[2])?;
                let min = parse_f64(&rec.args[3])?;
                Arc::new(GeomTorus::new(
                    GpTorus::new(self.resolve_axis2(ax)?.to_ax3(), maj, min)
                        .map_err(|e| format!("TOROIDAL_SURFACE: {e}"))?,
                ))
            }
            "SURFACE_OF_REVOLUTION" => {
                // SURFACE_OF_REVOLUTION(name, axis, generatrix) — the axis is an
                // AXIS1_PLACEMENT and the generatrix a curve. STEP order is
                // (axis, curve); FreeCAD emits (curve, axis). Detect which
                // argument holds the AXIS1_PLACEMENT by its record type.
                let r1 = parse_ref(&rec.args[1]);
                let r2 = parse_ref(&rec.args[2]);
                let is_axis1 = |r: usize| {
                    self.record(r)
                        .map(|rec| rec.type_name == "AXIS1_PLACEMENT")
                        .unwrap_or(false)
                };
                let (axis_ref, gen_ref) = match (r1, r2) {
                    (Some(a), Some(g)) if is_axis1(a) => (a, g),
                    (Some(g), Some(a)) if is_axis1(a) => (a, g),
                    (Some(a), Some(g)) => (a, g), // fallback: STEP order
                    (Some(a), None) => (a, a),
                    _ => return Err("SURFACE_OF_REVOLUTION: bad refs".into()),
                };
                let axis = self.resolve_axis1(axis_ref)?;
                let generatrix = self.resolve_curve(gen_ref)?;
                Arc::new(GeomSurfaceOfRevolution::new(generatrix, axis))
            }
            "B_SPLINE_SURFACE" => {
                // Plain B-spline surface (no explicit knots): clamped uniform
                // knot vectors derived from the pole grid and degrees.
                let deg_u = parse_f64(&rec.args[1])? as usize;
                let deg_v = parse_f64(&rec.args[2])? as usize;
                let poles: Vec<Vec<GpPnt>> = parse_nested_ref_list(&rec.args[3])
                    .into_iter()
                    .map(|row| {
                        row.into_iter()
                            .map(|r| self.resolve_point(r))
                            .collect::<Result<Vec<_>, _>>()
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                let nu = poles.len();
                let nv = poles.first().map_or(0, |r| r.len());
                let u_knots = uniform_knots_for(nu, deg_u);
                let v_knots = uniform_knots_for(nv, deg_v);
                Arc::new(
                    GeomBSplineSurface::new(poles, u_knots, v_knots, deg_u, deg_v)
                        .map_err(|e| format!("B_SPLINE_SURFACE: {e}"))?,
                )
            }
            "B_SPLINE_SURFACE_WITH_KNOTS" => {
                // Standard ISO 10303-42 layout (13 args): name, u_degree,
                // v_degree, control_points grid, surface_form, closed_u,
                // closed_v, self_intersect, u_multiplicities, v_multiplicities,
                // u_knots, v_knots, knot_spec. A merged rational surface carries
                // the weight grid as an optional 14th arg.
                let deg_u = parse_f64(&rec.args[1])? as usize;
                let deg_v = parse_f64(&rec.args[2])? as usize;
                let poles: Vec<Vec<GpPnt>> = parse_nested_ref_list(&rec.args[3])
                    .into_iter()
                    .map(|row| {
                        row.into_iter()
                            .map(|r| self.resolve_point(r))
                            .collect::<Result<Vec<_>, _>>()
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                let weights_arg = rec.args.get(13).map(|s| s.trim().to_string());
                let u_knots = expand_knots(
                    &parse_usize_list(rec.args.get(8).map(|s| s.as_str()).unwrap_or("()")),
                    &parse_real_list(rec.args.get(10).map(|s| s.as_str()).unwrap_or("()")),
                );
                let v_knots = expand_knots(
                    &parse_usize_list(rec.args.get(9).map(|s| s.as_str()).unwrap_or("()")),
                    &parse_real_list(rec.args.get(11).map(|s| s.as_str()).unwrap_or("()")),
                );
                let surface = match weights_arg.as_deref() {
                    None | Some("SELF") => {
                        GeomBSplineSurface::new(poles, u_knots, v_knots, deg_u, deg_v)
                    }
                    Some(w) => {
                        let wgrid = parse_nested_real_list(w);
                        GeomBSplineSurface::rational(poles, wgrid, u_knots, v_knots, deg_u, deg_v)
                    }
                };
                Arc::new(surface.map_err(|e| format!("B_SPLINE_SURFACE: {e}"))?)
            }
            "OFFSET_SURFACE" => {
                // Layout: (name, basis_surface, distance, self_intersect).
                let basis_ref = parse_ref(&rec.args[1]).ok_or("OFFSET_SURFACE: bad basis ref")?;
                let distance = parse_f64(&rec.args[2])?;
                let basis = self.resolve_surface(basis_ref)?;
                Arc::new(GeomOffsetSurface::new(basis, distance))
            }
            other => {
                self.warn(format!("unsupported surface entity {other} (#{id})"));
                return Err(format!("unsupported surface entity {other} (#{id})"));
            }
        };
        self.surface_cache.borrow_mut().insert(id, surface.clone());
        Ok(surface)
    }

    /// Resolve a representation record into (name, shapes).
    pub(super) fn resolve_representation(&self, id: usize) -> Result<(String, Vec<TopoShape>), String> {
        let rec = self.record(id)?;
        let name = parse_str(&rec.args[0]);
        let items = parse_ref_list(&rec.args[1]);
        // A representation's item list may mix the geometric/topological shape
        // entities (MANIFOLD_SOLID_BREP, …) with auxiliary placement / point /
        // direction entities (AXIS2_PLACEMENT_3D, CARTESIAN_POINT, DIRECTION,
        // VECTOR, …) that are referenced *by* the shapes but are not themselves
        // top-level shapes. OCCT's STEP reader skips these non-shape items; we
        // do the same so a placement next to the solid does not fail the whole
        // representation.
        let mut shapes = Vec::with_capacity(items.len());
        for &it in &items {
            match self.resolve_shape(it) {
                Ok(s) => shapes.push(s),
                Err(e) => {
                    if !self.is_auxiliary_entity(it) {
                        return Err(e);
                    }
                    // Auxiliary item (placement/point/direction/…): skip.
                }
            }
        }
        Ok((name, shapes))
    }

    /// Whether the record at `id` is an auxiliary geometric entity that a
    /// representation lists alongside its shapes but that is not a shape.
    pub(super) fn is_auxiliary_entity(&self, id: usize) -> bool {
        let Some(rec) = self.records.get(&id) else { return false };
        matches!(
            rec.type_name.as_str(),
            "AXIS2_PLACEMENT_3D"
                | "AXIS2_PLACEMENT_2D"
                | "AXIS1_PLACEMENT"
                | "CARTESIAN_POINT"
                | "DIRECTION"
                | "VECTOR"
                | "GEOMETRIC_REPRESENTATION_CONTEXT"
                | "REPRESENTATION_CONTEXT"
                | "PARAMETRIC_REPRESENTATION_CONTEXT"
                | "GLOBAL_UNIT_ASSIGNED_CONTEXT"
                | "APPLICATION_CONTEXT"
                | "PRODUCT_CONTEXT"
                | "PRODUCT_DEFINITION_CONTEXT"
                | "LENGTH_UNIT"
                | "PLANE_ANGLE_UNIT"
                | "SOLID_ANGLE_UNIT"
                | "SI_UNIT"
                | "NAMED_UNIT"
                | "UNCERTAINTY_MEASURE_WITH_UNIT"
        )
    }
}
/// `ShapeAnalysis_Curve::ProjectOnSegments` (`ShapeAnalysis_Curve.cxx:81-122`).
fn project_on_segments(
    curve: &dyn Curve,
    point: &GpPnt,
    segment_count: i32,
    start: &mut f64,
    end: &mut f64,
    proj_dist: &mut f64,
    proj_param: &mut f64,
) {
    if segment_count <= 0 {
        return;
    }
    let step = (*end - *start) / segment_count as f64;
    let mut min_sq = *proj_dist * *proj_dist;
    let mut changed = false;
    for i in 0..=segment_count {
        let u = *start + step * i as f64;
        let q = curve.d0(u);
        let sq = point.square_distance(&q);
        if sq < min_sq {
            min_sq = sq;
            *proj_param = u;
            changed = true;
        }
    }
    if changed {
        *proj_dist = min_sq.sqrt();
    }
    *end = (*end).min(*proj_param + step);
    *start = (*start).max(*proj_param - step);
}

/// `ElCLib::LineParameter` (`ElCLib.cxx:1192-1195`).
fn elclib_line_parameter(curve: &dyn Curve, p: &GpPnt) -> f64 {
    let origin = curve.d0(0.0);
    let dir = curve.d1(0.0).1;
    p.coord
        .subtracted(&origin.coord)
        .dot(dir.xyz())
}

/// `ShapeAnalysis_Curve::ProjectAct` (`cxx:265-497`).
fn shape_analysis_project_act(
    curve: &dyn Curve,
    point: &GpPnt,
    preci: f64,
    u_min: f64,
    u_max: f64,
) -> Option<(f64, f64)> {
    let mut ok = false;
    let mut proj_param = 0.0;
    let mut computed_param = 0.0;
    let mut computed_dist = f64::INFINITY;
    let mut have_old = false;
    let mut old_param = 0.0;
    let mut old_dist = f64::INFINITY;
    let mut mod_min = f64::INFINITY;
    if let Some((t, d)) = crate::int_tools_vertex_line::extrema_project_in_range(
        curve, point, u_min, u_max,
    ) {
        have_old = true;
        old_param = t;
        old_dist = d;
        computed_param = t;
        computed_dist = d;
        proj_param = t;
        mod_min = d;
        ok = d <= preci;
    }

    let mut closed = false;
    let mut period = 0.0;
    if ok {
        let lo = curve.d0(u_min);
        let hi = curve.d0(u_max);
        if lo.distance(&hi) <= occt_core::precision::CONFUSION || curve.is_periodic() {
            closed = true;
            period = u_max - u_min;
        }
    }

    if !ok {
        if curve.is_line() {
            proj_param = elclib_line_parameter(curve, point);
        } else if let Some(c) = curve.gp_circ() {
            let loc = c.position().location();
            if c.radius() <= occt_core::precision::RESOLUTION
                || point.square_distance(&loc) <= occt_core::precision::RESOLUTION
            {
                proj_param = u_min;
            } else {
                proj_param = occt_core::elib::clib::circle_parameter(&c.position(), point);
            }
            closed = true;
            period = 2.0 * std::f64::consts::PI;
        } else {
            let mut seg_lo = u_min;
            let mut seg_hi = u_max;
            let mut dist = f64::INFINITY;
            project_on_segments(
                curve,
                point,
                25,
                &mut seg_lo,
                &mut seg_hi,
                &mut dist,
                &mut proj_param,
            );
            if dist <= preci {
                return Some((proj_param, dist));
            }
            if let Some((t, q)) = crate::int_tools_vertex_line::extrema_locate_ext_pc(
                curve, point, proj_param, u_min, u_max,
            ) {
                let d_newton = point.distance(&q);
                if d_newton < mod_min {
                    return Some((t, d_newton));
                }
            }
            for n in [40, 20, 25, 40] {
                project_on_segments(
                    curve,
                    point,
                    n,
                    &mut seg_lo,
                    &mut seg_hi,
                    &mut dist,
                    &mut proj_param,
                );
                if dist <= preci {
                    return Some((proj_param, dist));
                }
            }
            if dist > mod_min && have_old {
                return Some((computed_param, computed_dist));
            }
            return Some((proj_param, dist));
        }
    }

    if closed && (proj_param < u_min || proj_param > u_max) && period.abs() > 1e-100 {
        // `ShapeAnalysis::AdjustByPeriod` (`ShapeAnalysis.cxx:48-62`).
        let to_val = 0.5 * (u_min + u_max);
        let diff = proj_param - to_val;
        let d = diff.abs();
        let p = period.abs();
        if d > 0.5 * p {
            let shift = if diff > 0.0 { -p } else { p } * (d / p + 0.5).floor();
            proj_param += shift;
        }
    }

    let q = curve.d0(proj_param);
    let new_dist = point.distance(&q);
    if have_old && old_dist * old_dist < new_dist * new_dist {
        return Some((old_param, old_dist));
    }
    Some((proj_param, new_dist))
}

/// `ShapeAnalysis_Curve::Project` (`cxx:147-201`) with `AdjustToEnds=false`
/// as `MakeFromCurve3D` (`TranslateEdge.cxx:443-444`).
fn shape_analysis_project(curve: &dyn Curve, p: &GpPnt, preci: f64) -> Option<(f64, f64)> {
    let (mut u_min, mut u_max) = (curve.first_parameter(), curve.last_parameter());
    if u_min > u_max {
        std::mem::swap(&mut u_min, &mut u_max);
    }
    if u_min.is_finite() && u_max.is_finite() {
        let low = curve.d0(u_min);
        let high = curve.d0(u_max);
        let dl = low.distance(p);
        if dl <= occt_core::precision::CONFUSION {
            return Some((u_min, dl));
        }
        let dh = high.distance(p);
        if dh <= occt_core::precision::CONFUSION {
            return Some((u_max, dh));
        }
        let closed = low.distance(&high) <= occt_core::precision::CONFUSION;
        if !closed {
            let delta = crate::int_tools_vertex_line::adaptor_resolution(curve, preci)
                .min((u_max - u_min) * 0.1);
            u_min -= delta;
            u_max += delta;
        }
    }
    shape_analysis_project_act(curve, p, preci, u_min, u_max)
}

/// `StepToTopoDS_TranslateEdge::MakeFromCurve3D` (`cxx:442-449`):
/// `ShapeAnalysis_Curve::Project` + `UpdateParam3d`. Returns the (possibly
/// reversed) curve together with the projected `[U1, U2]` range.
pub(super) fn edge_from_curve3d(
    curve: Arc<dyn Curve>,
    p1: &GpPnt,
    p2: &GpPnt,
) -> (Arc<dyn Curve>, f64, f64) {
    const PRECI: f64 = 1e-3;
    let (first, last) = edge_params_for_curve(curve.as_ref(), p1, p2);
    let mut w1 = first;
    let mut w2 = last;
    let curve = update_param3d(curve, &mut w1, &mut w2, PRECI);
    let curve = shift_displaced_line(curve, p1, p2, w1, w2, PRECI);
    (curve, w1, w2)
}

/// `MakeFromCurve3D` (`TranslateEdge.cxx:451-472`): a line that keeps a
/// constant offset from both vertices and whose parameter span equals the
/// vertex distance is translated onto the vertices.
fn shift_displaced_line(
    curve: Arc<dyn Curve>,
    p1: &GpPnt,
    p2: &GpPnt,
    w1: f64,
    w2: f64,
    preci: f64,
) -> Arc<dyn Curve> {
    if !curve.is_line() {
        return curve;
    }
    let temp1 = curve.d0(w1).distance(p1);
    let temp2 = curve.d0(w2).distance(p2);
    if temp1 <= preci && temp2 <= preci {
        return curve;
    }
    if (temp1 - temp2).abs() >= preci {
        return curve;
    }
    if ((w2 - w1).abs() - p1.distance(p2)).abs() >= occt_core::precision::CONFUSION {
        return curve;
    }
    let d1 = curve.d1(0.0).1;
    let Ok(dir) = GpDir::from_vec(&d1) else {
        return curve;
    };
    let origin = p1.translated_vec(&GpVec::from_xyz(&dir.xyz().multiplied(-w1)));
    Arc::new(GeomLine::from_pnt_dir(origin, dir))
}

/// Compute an edge's parameter range from its endpoint vertex points, based on
/// the reconstructed curve's analytic type.
pub(super) fn edge_params_for_curve(curve: &dyn Curve, p1: &GpPnt, p2: &GpPnt) -> (f64, f64) {
    let (f, l) = (curve.first_parameter(), curve.last_parameter());
    // Bounded non-periodic (BSpline / trimmed / Bezier): Project the vertices
    // (`MakeFromCurve3D` cxx:442-444). The previous natural-range early-return
    // used the whole knot domain for an EDGE_CURVE that only spans a portion,
    // so `ComputeDeflection` vertex-adjust and `MaxFaceTolerance` after
    // SameParameter inflated the face deflection (Shape-2: 18 vs Prs3d 0.6).
    if f.is_finite() && l.is_finite() && !curve.is_periodic() {
        const PRECI: f64 = 1e-3;
        // Vertices already near the natural ends: keep the knot domain.
        // Extrema Project on Shape.step's two [0,1] BSplines (nat=1.85)
        // still shrinks the 3D range to ~0.89 against a full-span pcurve
        // (6134 -> 6345). Shape-2 portion edges have nat 3..9, above this
        // floor. Leftover: CheckPCurves remap after Extrema-accurate ends.
        const ALIGNED: f64 = 2.0;
        if p1.distance(p2) < PRECI {
            return if l > f { (f, l) } else { (0.0, 1.0) };
        }
        let pf = curve.d0(f);
        let pl = curve.d0(l);
        if pf.distance(p1) <= ALIGNED && pl.distance(p2) <= ALIGNED {
            return (f, l);
        }
        if pf.distance(p2) <= ALIGNED && pl.distance(p1) <= ALIGNED {
            return (l, f);
        }
        if let (Some(a), Some(b)) = (
            shape_analysis_project(curve, p1, PRECI),
            shape_analysis_project(curve, p2, PRECI),
        ) {
            if (a.0 - b.0).abs() > occt_core::precision::PCONFUSION {
                return (a.0, b.0);
            }
        }
        return if l > f { (f, l) } else { (0.0, 1.0) };
    }
    // `StepToTopoDS_TranslateEdge::MakeFromCurve3D` (`cxx:442-446`):
    // `ShapeAnalysis_Curve::Project` then `UpdateParam3d`. For `Geom_Circle`
    // that is `ElCLib::CircleParameter` + periodic `AdjustPeriodic`.
    // Closed full-period edges (same vertex) keep the sampled atan2 branch:
    // `CircleParameter` normalizes into `[0, 2pi]`, so `(U, U+2pi)` can sit at
    // `[3pi/2, 7pi/2]` and this port's periodic UV box over-tessellates
    // `Torus.step` (2112/4032 vs occ 1369/2592). OCCT still uses that range;
    // the leftover is the UV-box interaction, not this Project call.
    if let Some(c) = curve.gp_circ() {
        if p1.distance(p2) >= 1e-9 {
            return circle_params_from_circ(&c, p1, p2);
        }
    }
    let (lo, hi) = if f.is_finite() && l.is_finite() && l > f {
        (f, l)
    } else {
        (0.0, 1.0)
    };
    match classify_curve(curve, lo, hi) {
        CurveKind::Line => {
            let origin = curve.d0(0.0);
            let d1 = curve.d1(0.0).1;
            let dir = GpVec::from_xyz(&d1.xyz().normalized());
            (
                GpVec::from_pnts(&origin, p1).dot(&dir),
                GpVec::from_pnts(&origin, p2).dot(&dir),
            )
        }
        CurveKind::Circle => circle_edge_params(curve, p1, p2, lo, hi),
        CurveKind::Ellipse => {
            let pa = curve.d0(lo);
            let pb = curve.d0(lo + (hi - lo) / 4.0);
            let pc = curve.d0(lo + (hi - lo) / 2.0);
            let center = midpoint(&pa, &pc);
            let a_major = center.distance(&pa).max(1e-30);
            let b_minor = center.distance(&pb).max(1e-30);
            let xdir = GpDir::from_vec(&GpVec::from_pnts(&center, &pa)).unwrap_or(dir_x());
            let ydir = GpDir::from_vec(&GpVec::from_pnts(&pb, &center)).unwrap_or(dir_y());
            let ang = |p: &GpPnt| {
                let v = GpVec::from_pnts(&center, p);
                (-v.xyz().dot(ydir.xyz()) / b_minor)
                    .atan2(v.xyz().dot(xdir.xyz()) / a_major)
            };
            let t1 = ang(p1);
            if p1.distance(p2) < 1e-9 {
                return (t1, t1 + 2.0 * PI);
            }
            let mut t2 = ang(p2);
            if t2 <= t1 {
                t2 += 2.0 * PI;
            }
            (t1, t2)
        }
        CurveKind::Parabola => {
            let vertex = curve.d0(0.0);
            let d1 = curve.d1(0.0).1;
            let ydir = GpDir::from_vec(&d1).unwrap_or(dir_y());
            let dir = GpVec::from_xyz(ydir.xyz());
            (
                GpVec::from_pnts(&vertex, p1).dot(&dir),
                GpVec::from_pnts(&vertex, p2).dot(&dir),
            )
        }
        CurveKind::Other => {
            if f.is_finite() && l.is_finite() && l > f {
                (f, l)
            } else {
                (0.0, 1.0)
            }
        }
    }
}

/// `StepToTopoDS_GeometricTool::UpdateParam3d` (`GeometricTool.cxx:227-409`).
fn update_param3d(
    curve: Arc<dyn Curve>,
    w1: &mut f64,
    w2: &mut f64,
    preci: f64,
) -> Arc<dyn Curve> {
    let cf = curve.first_parameter();
    let cl = curve.last_parameter();
    let start = curve.d0(cf);
    let end = curve.d0(cl);
    let is_closed = start.distance(&end) <= preci;

    // Bounded and not closed: clamp projections onto the natural domain
    // (`cxx:238-268`).
    if cf.is_finite() && cl.is_finite() && !curve.is_periodic() && !is_closed {
        if *w1 < cf {
            *w1 = cf;
        } else if *w1 > cl {
            *w1 = cl;
        }
        if *w2 < cf {
            *w2 = cf;
        } else if *w2 > cl {
            *w2 = cl;
        }
    }

    if *w1 < *w2 {
        return curve;
    }

    if curve.is_periodic() {
        crate::geom_bnd_lib_elclib2d::adjust_periodic(
            cf,
            cl,
            occt_core::precision::PCONFUSION,
            w1,
            w2,
        );
        return curve;
    }

    if is_closed {
        if (*w2 - cf).abs() < occt_core::precision::PCONFUSION {
            *w2 = cl;
        } else if (*w1 - cl).abs() < occt_core::precision::PCONFUSION {
            *w1 = cf;
        } else {
            if curve.d0(*w1).distance(&start) < preci {
                *w1 = cf;
            }
            if curve.d0(*w2).distance(&end) < preci {
                *w2 = cl;
            }
            if (*w2 - *w1).abs() < occt_core::precision::PCONFUSION {
                *w1 = cf;
                *w2 = cl;
            } else if *w1 > *w2 {
                std::mem::swap(w1, w2);
            }
        }
        return curve;
    }

    if curve.bspline_poles().is_some() && *w1 > *w2 {
        // Non-closed BSpline with reversed params: Reverse the per-edge copy
        // (`cxx:364-372`). The STEP curve cache is shared; do not mutate it.
        let rp = |w: f64| cf + cl - w;
        *w1 = rp(*w1);
        *w2 = rp(*w2);
        let rev: Arc<dyn Curve> = Arc::from(curve.reversed());
        if (*w1 - *w2).abs() < occt_core::precision::PCONFUSION {
            *w1 = rev.first_parameter();
            *w2 = rev.last_parameter();
        }
        return rev;
    }

    if *w1 > *w2 {
        std::mem::swap(w1, w2);
    }
    curve
}

/// `ElCLib::CircleParameter` + `StepToTopoDS_GeometricTool::UpdateParam3d`
/// periodic arm (`GeometricTool.cxx:270-279`).
fn circle_params_from_circ(c: &GpCirc, p1: &GpPnt, p2: &GpPnt) -> (f64, f64) {
    let pos = c.position();
    let mut w1 = occt_core::elib::clib::circle_parameter(&pos, p1);
    let mut w2 = occt_core::elib::clib::circle_parameter(&pos, p2);
    // Closed circle: `BRep_Tool::Range` starts at the vertex parameter,
    // not always 0. `[0, 2pi]` makes `Value(first)` the opposite point on
    // the circle, so `ComputeDeflection`'s vertex-adjust floor becomes
    // the diameter and tessellation stops at 2 points (`Torus.step` minor
    // seam). Source: `BRep_Tool::Range` + `BRepMesh_Deflection.cxx:86-94`.
    if p1.distance(p2) < 1e-9 {
        return (w1, w1 + 2.0 * PI);
    }
    if w1 < w2 {
        return (w1, w2);
    }
    crate::geom_bnd_lib_elclib2d::adjust_periodic(
        0.0,
        2.0 * PI,
        occt_core::precision::PCONFUSION,
        &mut w1,
        &mut w2,
    );
    (w1, w2)
}

/// Sampled-circle fallback when the 3D curve is not a `Geom_Circle`.
fn circle_edge_params(
    curve: &dyn Curve,
    p1: &GpPnt,
    p2: &GpPnt,
    lo: f64,
    hi: f64,
) -> (f64, f64) {
    let pa = curve.d0(lo);
    let pb = curve.d0(lo + (hi - lo) / 4.0);
    let pc = curve.d0(lo + (hi - lo) / 2.0);
    let center = circle_center3(&pa, &pb, &pc).unwrap_or_else(GpPnt::zero);
    let xdir = GpDir::from_vec(&GpVec::from_pnts(&center, &pa)).unwrap_or(dir_x());
    let ydir = GpDir::from_vec(&GpVec::from_pnts(&center, &pb)).unwrap_or(dir_y());
    let ang = |p: &GpPnt| {
        let v = GpVec::from_pnts(&center, p);
        v.xyz().dot(ydir.xyz()).atan2(v.xyz().dot(xdir.xyz()))
    };
    let t1 = ang(p1);
    if p1.distance(p2) < 1e-9 {
        return (t1, t1 + 2.0 * PI);
    }
    let mut t2 = ang(p2);
    if t2 <= t1 {
        t2 += 2.0 * PI;
    }
    (t1, t2)
}

pub(super) fn parse_usize_list(s: &str) -> Vec<usize> {
    let s = s.trim();
    if !(s.starts_with('(') && s.ends_with(')')) {
        return Vec::new();
    }
    split_top(&s[1..s.len() - 1])
        .into_iter()
        .filter_map(|a| a.trim().parse().ok())
        .collect()
}

pub(super) fn parse_real_list(s: &str) -> Vec<f64> {
    let s = s.trim();
    if !(s.starts_with('(') && s.ends_with(')')) {
        return Vec::new();
    }
    split_top(&s[1..s.len() - 1])
        .into_iter()
        .filter_map(|a| a.trim().parse().ok())
        .collect()
}

/// Split a nested `((...),(...),...)` argument into its inner lists.
pub(super) fn parse_nested_lists(s: &str) -> Vec<Vec<String>> {
    let s = s.trim();
    if !(s.starts_with('(') && s.ends_with(')')) {
        return Vec::new();
    }
    split_top(&s[1..s.len() - 1])
        .into_iter()
        .filter_map(|row| {
            let row = row.trim();
            if row.starts_with('(') && row.ends_with(')') {
                Some(split_top(&row[1..row.len() - 1]))
            } else {
                None
            }
        })
        .collect()
}

/// Parse a nested list of entity references (e.g. a B-spline surface pole grid).
pub(super) fn parse_nested_ref_list(s: &str) -> Vec<Vec<usize>> {
    parse_nested_lists(s)
        .into_iter()
        .map(|row| row.into_iter().filter_map(|a| parse_ref(&a)).collect())
        .collect()
}

/// Parse a nested list of real values (e.g. a B-spline surface weight grid).
pub(super) fn parse_nested_real_list(s: &str) -> Vec<Vec<f64>> {
    parse_nested_lists(s)
        .into_iter()
        .map(|row| row.into_iter().filter_map(|a| a.trim().parse().ok()).collect())
        .collect()
}

/// Parse a STEP physical file, returning the model plus collected warnings.
/// Extract the `DATA` section of a physical file (between `DATA;` and the
/// closing `ENDSEC;`), validating the file terminator.
pub(super) fn data_section(content: &str) -> Result<&str, String> {
    if !content.contains("END-ISO-10303-21") {
        return Err("STEP file: missing END-ISO-10303-21 terminator".into());
    }
    let data_start = content
        .find("DATA;")
        .ok_or("STEP file: missing DATA section")?;
    let after_data = &content[data_start + 5..];
    let data_end = after_data
        .find("ENDSEC;")
        .ok_or("STEP file: DATA section not closed by ENDSEC")?;
    Ok(&after_data[..data_end])
}

pub(super) fn read_step_impl(content: &str) -> Result<(BRepModel, Vec<String>), String> {
    let data = data_section(content)?;
    let records = parse_records(data)?;
    let resolver = Resolver::new(&records);

    let mut rep_ids: Vec<usize> = records
        .iter()
        .filter(|(_, r)| {
            r.type_name == "ADVANCED_BREP_SHAPE_REPRESENTATION" || r.type_name == "SHAPE_REPRESENTATION"
        })
        .map(|(id, _)| *id)
        .collect();
    rep_ids.sort_unstable();

    let mut model = BRepModel::new();
    if !rep_ids.is_empty() {
        for id in rep_ids {
            match resolver.resolve_representation(id) {
                Ok((name, shapes)) => {
                    if shapes.is_empty() {
                        continue;
                    }
                    let shape = if shapes.len() == 1 {
                        shapes.into_iter().next().unwrap()
                    } else {
                        let b = TopoBuilder::new();
                        b.make_compound_of(&shapes).0
                    };
                    model.add(&name, shape);
                }
                Err(e) => resolver.warn(format!("representation #{id}: {e}")),
            }
        }
    } else {
        // Fallback for bare files without a representation layer: treat any
        // MANIFOLD_SOLID_BREP as a root shape.
        let mut roots: Vec<usize> = records
            .iter()
            .filter(|(_, r)| r.type_name == "MANIFOLD_SOLID_BREP")
            .map(|(id, _)| *id)
            .collect();
        roots.sort_unstable();
        for id in roots {
            match resolver.resolve_shape(id) {
                Ok(s) => {
                    model.add("", s);
                }
                Err(e) => resolver.warn(format!("shape #{id}: {e}")),
            }
        }
    }

    let warnings = resolver.warnings.into_inner();
    Ok((model, warnings))
}

/// Parse a STEP physical file into a `BRepModel`.
pub fn read_step(content: &str) -> Result<BRepModel, String> {
    Ok(read_step_impl(content)?.0)
}

/// Parse a STEP physical file, also returning collected warnings for skipped
/// or unsupported entities.
pub fn read_step_with_warnings(content: &str) -> Result<(BRepModel, Vec<String>), String> {
    read_step_impl(content)
}

/// Read a STEP physical file from disk.
pub fn read_step_file(path: &str) -> Result<BRepModel, String> {
    let content = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
    read_step(&content)
}

/// Read a STEP physical file back into a [`StepAssembly`].
///
/// Reconstructs the product tree from `PRODUCT` / `PRODUCT_DEFINITION` /
/// `NEXT_ASSEMBLY_USAGE_OCCURRENCE` records. The root name is the product that
/// is never referenced as a child; part shapes are re-read from the part
/// representations when present (parts without a representation keep an empty
/// `TopoShape`).
pub fn read_step_assembly(content: &str) -> Result<StepAssembly, String> {
    let data = data_section(content)?;
    let records = parse_records(data)?;
    let resolver = Resolver::new(&records);

    // PRODUCT id -> product name.
    let mut product_names: HashMap<usize, String> = HashMap::new();
    for (id, rec) in &records {
        if rec.type_name == "PRODUCT" {
            product_names.insert(*id, parse_str(&rec.args[0]));
        }
    }

    // PRODUCT_DEFINITION_FORMATION id -> owning PRODUCT id.
    let mut form_product: HashMap<usize, usize> = HashMap::new();
    for (id, rec) in &records {
        if rec.type_name == "PRODUCT_DEFINITION_FORMATION" {
            if let Some(p) = parse_ref(&rec.args[2]) {
                form_product.insert(*id, p);
            }
        }
    }

    // PRODUCT_DEFINITION id -> product name (via its formation).
    let mut def_name: HashMap<usize, String> = HashMap::new();
    for (id, rec) in &records {
        if rec.type_name == "PRODUCT_DEFINITION" {
            let name = parse_ref(&rec.args[3])
                .and_then(|f| form_product.get(&f).copied())
                .and_then(|p| product_names.get(&p).cloned())
                .unwrap_or_default();
            def_name.insert(*id, name);
        }
    }

    // PRODUCT_DEFINITION_SHAPE id -> PRODUCT id (the shape's owner).
    let mut pds_product: HashMap<usize, usize> = HashMap::new();
    for (id, rec) in &records {
        if rec.type_name == "PRODUCT_DEFINITION_SHAPE" {
            if let Some(p) = parse_ref(&rec.args[2]) {
                pds_product.insert(*id, p);
            }
        }
    }

    // representation id -> PRODUCT id, via the SHAPE_REPRESENTATION link.
    let mut rep_product: HashMap<usize, usize> = HashMap::new();
    for (id, rec) in &records {
        if rec.type_name == "PRODUCT_DEFINITION_SHAPE_REPRESENTATION" {
            if let (Some(pds), Some(rep)) = (parse_ref(&rec.args[1]), parse_ref(&rec.args[2])) {
                if let Some(p) = pds_product.get(&pds) {
                    rep_product.insert(rep, *p);
                }
            }
        }
    }

    // Parts: resolve each representation owned by a product.
    let mut products: Vec<(String, TopoShape)> = Vec::new();
    let mut seen: HashSet<String> = HashSet::new();
    let mut reps: Vec<usize> = records
        .iter()
        .filter(|(_, r)| {
            r.type_name == "ADVANCED_BREP_SHAPE_REPRESENTATION" || r.type_name == "SHAPE_REPRESENTATION"
        })
        .map(|(id, _)| *id)
        .collect();
    reps.sort_unstable();
    for id in reps {
        let Some(pid) = rep_product.get(&id) else { continue };
        let Some(name) = product_names.get(pid) else { continue };
        let name = name.clone();
        if seen.contains(&name) {
            continue;
        }
        match resolver.resolve_representation(id) {
            Ok((_, shapes)) => {
                if shapes.is_empty() {
                    continue;
                }
                let shape = if shapes.len() == 1 {
                    shapes.into_iter().next().unwrap()
                } else {
                    let b = TopoBuilder::new();
                    b.make_compound_of(&shapes).0
                };
                seen.insert(name.clone());
                products.push((name, shape));
            }
            Err(e) => resolver.warn(format!("representation #{id}: {e}")),
        }
    }

    // Assembly tree from NEXT_ASSEMBLY_USAGE_OCCURRENCE records.
    let mut children: Vec<(String, Vec<String>)> = Vec::new();
    let mut has_parent: HashSet<String> = HashSet::new();
    for (_, rec) in &records {
        if rec.type_name == "NEXT_ASSEMBLY_USAGE_OCCURRENCE" {
            if let (Some(rel), Some(red)) = (parse_ref(&rec.args[3]), parse_ref(&rec.args[4])) {
                if let (Some(pn), Some(cn)) = (def_name.get(&rel), def_name.get(&red)) {
                    if let Some(entry) = children.iter_mut().find(|(p, _)| p == pn) {
                        entry.1.push(cn.clone());
                    } else {
                        children.push((pn.clone(), vec![cn.clone()]));
                    }
                    has_parent.insert(cn.clone());
                }
            }
        }
    }

    // The assembly root is the product never referenced as a child.
    let name = product_names
        .values()
        .find(|n| !has_parent.contains(*n))
        .cloned()
        .unwrap_or_else(|| "Assembly".to_string());

    Ok(StepAssembly {
        name,
        products,
        children,
    })
}

/// Read a STEP assembly physical file from disk.
pub fn read_step_assembly_file(path: &str) -> Result<StepAssembly, String> {
    let content = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
    read_step_assembly(&content)
}

/// Round-trip helper: write `shape` to STEP, read it back, and count the
/// distinct vertices / edges / faces in the reconstructed model.
pub fn step_roundtrip_counts(shape: &TopoShape) -> Result<(usize, usize, usize), String> {
    let step = write_shape_step(shape);
    let model = read_step(&step)?;
    let mut nv = 0usize;
    let mut ne = 0usize;
    let mut nf = 0usize;
    for ms in &model.shapes {
        let c = crate::topo_tools_full::shape_counts(&ms.shape);
        nv += c.get(&ShapeType::Vertex).copied().unwrap_or(0);
        ne += c.get(&ShapeType::Edge).copied().unwrap_or(0);
        nf += c.get(&ShapeType::Face).copied().unwrap_or(0);
    }
    Ok((nv, ne, nf))
}
