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
                // `StepToGeom::MakeCylindricalSurface` (`StepToGeom.cxx:1444-1453`):
                // `SS->Radius() * LengthFactor()` (`cxx:1452`).
                let r = parse_f64(&rec.args[2])? * self.length_factor;
                Arc::new(GeomCylinder::new(
                    GpCylinder::new(self.resolve_axis2(ax)?.to_ax3(), r)
                        .map_err(|e| format!("CYLINDRICAL_SURFACE: {e}"))?,
                ))
            }
            "CONICAL_SURFACE" => {
                let ax = parse_ref(&rec.args[1]).ok_or("CONICAL_SURFACE: bad axis ref")?;
                // `StepToGeom::MakeConicalSurface` (`StepToGeom.cxx:1307-1321`):
                // `R = SS->Radius() * theLocalFactors.LengthFactor()` (`cxx:1315`),
                // `Ang = SS->SemiAngle() * theLocalFactors.PlaneAngleFactor()`
                // (`cxx:1316`). OCCT then floors the angle at
                // `Precision::Angular()` (`cxx:1319`); that floor is unported here.
                let r = parse_f64(&rec.args[2])? * self.length_factor;
                let a = parse_f64(&rec.args[3])? * self.plane_angle_factor;
                Arc::new(GeomCone::new(
                    GpCone::new(self.resolve_axis2(ax)?.to_ax3(), r, a)
                        .map_err(|e| format!("CONICAL_SURFACE: {e}"))?,
                ))
            }
            "SPHERICAL_SURFACE" => {
                let ax = parse_ref(&rec.args[1]).ok_or("SPHERICAL_SURFACE: bad axis ref")?;
                // `StepToGeom::MakeSphericalSurface` (`StepToGeom.cxx:1893-1901`):
                // `SS->Radius() * LengthFactor()` (`cxx:1899`).
                let r = parse_f64(&rec.args[2])? * self.length_factor;
                Arc::new(GeomSphere::new(
                    GpSphere::new(self.resolve_axis2(ax)?.to_ax3(), r)
                        .map_err(|e| format!("SPHERICAL_SURFACE: {e}"))?,
                ))
            }
            "TOROIDAL_SURFACE" => {
                let ax = parse_ref(&rec.args[1]).ok_or("TOROIDAL_SURFACE: bad axis ref")?;
                // `StepToGeom::MakeToroidalSurface` (`StepToGeom.cxx:2104-2115`):
                // `std::abs(MajorRadius * LF)`, `std::abs(MinorRadius * LF)`
                // (`cxx:2111-2113`). OCCT takes the absolute value; this port
                // keeps the raw sign because `step_surface_is_reversed`
                // (`StepToTopoDS_TranslateFace.cxx:479-491`) already handles the
                // negative-major-radius face orientation.
                let maj = parse_f64(&rec.args[2])? * self.length_factor;
                let min = parse_f64(&rec.args[3])? * self.length_factor;
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
            "SURFACE_OF_LINEAR_EXTRUSION" => {
                // SURFACE_OF_LINEAR_EXTRUSION(name, swept_curve, extrusion_axis).
                // `RWStepGeom_RWSurfaceOfLinearExtrusion.cxx:44-55` reads
                // `swept_curve` from attribute 2 and `extrusion_axis` from
                // attribute 3, so the order is positional.
                let curve_ref = parse_ref(&rec.args[1])
                    .ok_or("SURFACE_OF_LINEAR_EXTRUSION: bad curve ref")?;
                let axis_ref = parse_ref(&rec.args[2])
                    .ok_or("SURFACE_OF_LINEAR_EXTRUSION: bad axis ref")?;
                // `StepToGeom::MakeSurfaceOfLinearExtrusion`
                // (`StepToGeom.cxx:2007-2029`): `gp_Dir D(V->Vec())`
                // (`cxx:2020`) drops the `VECTOR` magnitude, which
                // `resolve_vector` has already scaled by `LengthFactor`.
                let axis_vec = self.resolve_vector(axis_ref)?;
                let direction = GpDir::from_vec(&axis_vec)
                    .map_err(|e| format!("SURFACE_OF_LINEAR_EXTRUSION: {e}"))?;
                let basis = self.resolve_curve(curve_ref)?;
                // `cxx:2021-2024`: a `Geom_Line` basis parallel to the extrusion
                // direction is degenerate and OCCT returns a null surface. The
                // port has no `Geom_Line` downcast, so `is_line()` (which a
                // `Geom_TrimmedCurve` would forward) plus the `!is_geom_trimmed()`
                // test stands in for it, and the constant tangent of the line
                // stands in for `Lin().Direction()`.
                if basis.is_line() && !basis.is_geom_trimmed() {
                    if let Ok(line_dir) = GpDir::from_vec(&basis.d1(0.0).1) {
                        if line_dir.is_parallel_tol(&direction, occt_core::precision::ANGULAR) {
                            return Err(
                                "SURFACE_OF_LINEAR_EXTRUSION: degenerate line basis".into()
                            );
                        }
                    }
                }
                Arc::new(GeomSurfaceOfLinearExtrusion::new(basis, direction))
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
                // `StepToGeom::MakeSurface` Offset arm (`StepToGeom.cxx:1940-1962`):
                // `anOffset = OS->Distance() * LengthFactor()` (`cxx:1947`).
                let distance = parse_f64(&rec.args[2])? * self.length_factor;
                let basis = self.resolve_surface(basis_ref)?;
                Arc::new(GeomOffsetSurface::new(basis, distance))
            }
            "RECTANGULAR_TRIMMED_SURFACE" => {
                // Layout: (name, basis_surface, u1, u2, v1, v2, usense, vsense).
                // `StepToGeom::MakeRectangularTrimmedSurface`
                // (`StepToGeom.cxx:1834-1885`): the basis surface is built first
                // (`:1838`), then the trim window is scaled by `uFact`/`vFact` chosen
                // from the basis type (`:1845-1875`) before
                // `Geom_RectangularTrimmedSurface` is constructed (`:1882`).
                let basis_ref = parse_ref(&rec.args[1])
                    .ok_or("RECTANGULAR_TRIMMED_SURFACE: bad basis ref")?;
                let basis = self.resolve_surface(basis_ref)?;
                let (u_fact, v_fact) = if basis.gp_sphere().is_some() || basis.gp_torus().is_some() {
                    (self.plane_angle_factor, self.plane_angle_factor)
                } else if basis.gp_cylinder().is_some() {
                    (self.plane_angle_factor, self.length_factor)
                } else if basis.is_surface_of_revolution() {
                    (self.plane_angle_factor, 1.0)
                } else if let Some(co) = basis.gp_cone() {
                    // `cxx:1866-1871`: `vFact = LengthFact / cos(SemiAngle)`.
                    (self.plane_angle_factor, self.length_factor / co.semi_angle().cos())
                } else if basis.gp_pln().is_some() {
                    (self.length_factor, self.length_factor)
                } else {
                    // `cxx:1845-1846`: `uFact = vFact = 1.` for any other basis.
                    (1.0, 1.0)
                };
                let u1 = parse_f64(&rec.args[2])? * u_fact;
                let u2 = parse_f64(&rec.args[3])? * u_fact;
                let v1 = parse_f64(&rec.args[4])? * v_fact;
                let v2 = parse_f64(&rec.args[5])? * v_fact;
                // UNPORTED: the `Usense`/`Vsense` flags (`rec.args[6]`/`[7]`) are the
                // reversal arguments of `Geom_RectangularTrimmedSurface`; the port's
                // `rectangular_trimmed::uv` carries no sense state and performs no
                // `SetTrim` normalisation (`Geom_RectangularTrimmedSurface.cxx`
                // `SetTrim`, including `ElCLib::AdjustPeriodic` for a periodic basis).
                Arc::new(GeomRectangularTrimmedSurface::uv(basis, u1, u2, v1, v2))
            }
            "BEZIER_SURFACE" | "UNIFORM_SURFACE" | "QUASI_UNIFORM_SURFACE" => {
                // `StepToGeom::MakeSurface` (`StepToGeom.cxx:522-743`): STEP Bezier,
                // uniform and quasi-uniform surfaces are each converted into a
                // `BSplineSurfaceWithKnots` before being mapped onto Geom:
                //   Bezier       -> knots {0,1}, multiplicities degree+1 (`:537-548`)
                //   Uniform      -> n_poles + degree + 1 knots i-1, mults 1 (`:570-593`)
                //   QuasiUniform -> n_poles - degree + 1 knots i-1, mults 1 except
                //                   both ends = degree+1 (`:612-638`)
                // The `_AND_RATIONAL_B_SPLINE_SURFACE` complex forms carry the
                // weight grid (`:644-687`, `:692-739`); a rational **Bezier**
                // surface is the exception — OCCT's Bezier arm never reads
                // `WeightsData` (`:524-555`), so its weights are dropped just as
                // OCCT drops them.
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
                let (u_knots, v_knots) = match rec.type_name.as_str() {
                    "BEZIER_SURFACE" => (bezier_knots(deg_u), bezier_knots(deg_v)),
                    "UNIFORM_SURFACE" => (
                        uniform_open_knots(nu, deg_u),
                        uniform_open_knots(nv, deg_v),
                    ),
                    _ => (
                        quasi_uniform_knots(nu, deg_u),
                        quasi_uniform_knots(nv, deg_v),
                    ),
                };
                // A merged rational complex carries the weight grid at index 4; a
                // plain entity has `surface_form` there, which never is a list.
                let weights_arg = rec.args.get(4).map(|s| s.trim().to_string());
                let weights = match weights_arg.as_deref() {
                    Some(w) if w.starts_with('(') && rec.type_name != "BEZIER_SURFACE" => {
                        Some(parse_nested_real_list(w))
                    }
                    _ => None,
                };
                let surface = match weights {
                    Some(wgrid) => GeomBSplineSurface::rational(
                        poles, wgrid, u_knots, v_knots, deg_u, deg_v,
                    ),
                    None => GeomBSplineSurface::new(poles, u_knots, v_knots, deg_u, deg_v),
                };
                Arc::new(surface.map_err(|e| format!("{}: {e}", rec.type_name))?)
            }
            "SURFACE_REPLICA" => {
                // `StepToGeom::MakeSurface` SurfaceReplica arm (`StepToGeom.cxx:1967-1986`):
                // `S1 = MakeSurface(ParentSurface)`, then `S1->Transform(T1)` with
                // `T1 = MakeTransformation3d(Transformation)`; the guard
                // `!T.IsNull() && PS != SS` (`cxx:1973`) rejects a cyclic replica.
                let parent_ref = parse_ref(&rec.args[1])
                    .ok_or("SURFACE_REPLICA: bad parent ref")?;
                let trsf_ref = parse_ref(&rec.args[2])
                    .ok_or("SURFACE_REPLICA: bad transformation ref")?;
                if parent_ref == id {
                    return Err(format!("SURFACE_REPLICA: cyclic parent (#{id})"));
                }
                let parent = self.resolve_surface(parent_ref)?;
                let t = self.make_transformation3d(trsf_ref)?;
                Arc::from(parent.transformed(&t))
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
///
/// OCCT's `ProjectAct` returns a parameter on every path (`ShapeAnalysis_Curve.cxx`
/// has no failure channel), so this port has no `Option` either.
fn shape_analysis_project_act(
    curve: &dyn Curve,
    point: &GpPnt,
    preci: f64,
    u_min: f64,
    u_max: f64,
) -> (f64, f64) {
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
        // `ShapeAnalysis_Curve.cxx:340-344`: `theCurve.IsClosed()` only.
        // Do not treat `is_periodic` alone as closed (invent vs cxx) — that
        // forced AdjustByPeriod on Extrema hits and could shrink / wrap
        // edge ranges after always-Project (`TranslateEdge.cxx:442-444`).
        let lo = curve.d0(u_min);
        let hi = curve.d0(u_max);
        if lo.distance(&hi) <= occt_core::precision::CONFUSION {
            closed = true;
            period = u_max - u_min;
        }
    }

    if !ok {
        // `ShapeAnalysis_Curve.cxx:355-477`: `switch (theCurve.GetType())` with a
        // precise `ElCLib::Parameter` arm for every analytic type, and the
        // segmented `ProjectOnSegments` search only for the `default:` case
        // (B-spline / Bezier / other).
        if let Some(c) = curve.gp_circ() {
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
        } else if let Some(h) = curve.gp_hyperbola() {
            proj_param = occt_core::elib::clib::parameter_hypr(&h, point);
        } else if let Some(prb) = curve.gp_parabola() {
            proj_param = occt_core::elib::clib::parameter_parab(&prb, point);
        } else if curve.is_line() {
            proj_param = elclib_line_parameter(curve, point);
        } else if let Some(el) = curve.gp_ellipse() {
            proj_param = occt_core::elib::clib::parameter_elips(&el, point);
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
                return (proj_param, dist);
            }
            if let Some((t, q)) = crate::int_tools_vertex_line::extrema_locate_ext_pc(
                curve, point, proj_param, u_min, u_max,
            ) {
                let d_newton = point.distance(&q);
                if d_newton < mod_min {
                    return (t, d_newton);
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
                    return (proj_param, dist);
                }
            }
            if dist > mod_min && have_old {
                return (computed_param, computed_dist);
            }
            return (proj_param, dist);
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
        return (old_param, old_dist);
    }
    (proj_param, new_dist)
}

/// `ShapeAnalysis_Curve::Project` (`cxx:147-201`) with `AdjustToEnds=false`
/// as `MakeFromCurve3D` (`TranslateEdge.cxx:443-444`).
fn shape_analysis_project(curve: &dyn Curve, p: &GpPnt, preci: f64) -> (f64, f64) {
    let (mut u_min, mut u_max) = (curve.first_parameter(), curve.last_parameter());
    if u_min > u_max {
        std::mem::swap(&mut u_min, &mut u_max);
    }
    if u_min.is_finite() && u_max.is_finite() {
        let low = curve.d0(u_min);
        let high = curve.d0(u_max);
        let dl = low.distance(p);
        if dl <= occt_core::precision::CONFUSION {
            return (u_min, dl);
        }
        let dh = high.distance(p);
        if dh <= occt_core::precision::CONFUSION {
            return (u_max, dh);
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

/// `StepToTopoDS_TranslateEdge::MakeFromCurve3D` (`cxx:442-444`): the edge range
/// is `ShapeAnalysis_Curve::Project` of both endpoint points with
/// `AdjustToEnds = false` — for **every** curve type, with no analytic-type
/// dispatch for the range itself and no endpoint-distance shortcut. `Project`
/// keeps its own endpoint early-exit for `Geom_BoundedCurve` (`cxx:161-182`) and
/// has no failure channel (`ShapeAnalysis_Curve.cxx:147-201` always produces a
/// parameter); `UpdateParam3d` (cxx:446) and the displaced-Line shift
/// (cxx:451-472) are applied by [`edge_from_curve3d`].
pub(super) fn edge_params_for_curve(curve: &dyn Curve, p1: &GpPnt, p2: &GpPnt) -> (f64, f64) {
    const PRECI: f64 = 1e-3;
    let (u1, _) = shape_analysis_project(curve, p1, PRECI);
    let (u2, _) = shape_analysis_project(curve, p2, PRECI);
    (u1, u2)
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
    // UNPORTED: `STEPControl_ActorRead::TransferEntity` composes each mapped
    // shape representation with the shapes related to it by
    // `SHAPE_REPRESENTATION_RELATIONSHIP` / `..._WITH_TRANSFORMATION`
    // (`TransferRelatedSRR`, `STEPControl_ActorRead.cxx:2061-2091`, and
    // `ComputeSRRWT`, `cxx:2493-2545`), and `read.step.shape.relationship`
    // defaults to true (`DESTEP_Parameters.hxx:170`). Here every related
    // representation is emitted as its own root shape instead of being composed
    // into the product shape. No sample in `data/` carries such a relation.

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

    // PRODUCT_DEFINITION_SHAPE id -> product name. Its `definition` field
    // (`characterized_definition`) is a PRODUCT_DEFINITION in a standard file
    // (`SHAPE_DEFINITION_REPRESENTATION`, as OCCT writes) and was a PRODUCT id
    // in this port's earlier files; accept both.
    let mut pds_name: HashMap<usize, String> = HashMap::new();
    for (id, rec) in &records {
        if rec.type_name == "PRODUCT_DEFINITION_SHAPE" {
            if let Some(d) = parse_ref(&rec.args[2]) {
                if let Some(n) = def_name.get(&d).or_else(|| product_names.get(&d)) {
                    pds_name.insert(*id, n.clone());
                }
            }
        }
    }

    // representation id -> product name, via the shape-definition link.
    //  * `SHAPE_DEFINITION_REPRESENTATION(definition, used_representation)` is
    //    the ISO 10303-42 entity OCCT's writer emits.
    //  * `PRODUCT_DEFINITION_SHAPE_REPRESENTATION` is this port's legacy record.
    let mut rep_name: HashMap<usize, String> = HashMap::new();
    for (_, rec) in &records {
        let (pds, rep) = match rec.type_name.as_str() {
            "SHAPE_DEFINITION_REPRESENTATION" => {
                (parse_ref(&rec.args[0]), parse_ref(&rec.args[1]))
            }
            "PRODUCT_DEFINITION_SHAPE_REPRESENTATION" => {
                (parse_ref(&rec.args[1]), parse_ref(&rec.args[2]))
            }
            _ => (None, None),
        };
        if let (Some(pds), Some(rep)) = (pds, rep) {
            if let Some(n) = pds_name.get(&pds) {
                rep_name.insert(rep, n.clone());
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
        let Some(name) = rep_name.get(&id) else { continue };
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
