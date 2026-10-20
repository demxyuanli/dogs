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
                if crate::debug_flags::iges_trace() {
                    eprintln!(
                        "TRACE surf#{id}: nu={} nv={} ku={} kv={} deg=({deg_u},{deg_v})",
                        surface.as_ref().map_or(0, |s| s.nb_poles_u()),
                        surface.as_ref().map_or(0, |s| s.nb_poles_v()),
                        surface.as_ref().map_or(0, |s| s.bspline_surface_uknots().map_or(0, |k| k.len())),
                        surface.as_ref().map_or(0, |s| s.bspline_surface_vknots().map_or(0, |k| k.len())),
                    );
                }
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
                // `SS->Usense()` / `SS->Vsense()` (`cxx:1882`): the reversal
                // arguments of the patch, used only when the basis is periodic
                // in that direction (see `SetTrim`).
                let u_sense = parse_logical(rec.args.get(6).map(String::as_str), true);
                let v_sense = parse_logical(rec.args.get(7).map(String::as_str), true);
                Arc::new(
                    GeomRectangularTrimmedSurface::uv_sense(
                        basis, u1, u2, v1, v2, u_sense, v_sense,
                    )
                    .map_err(|e| format!("RECTANGULAR_TRIMMED_SURFACE: {e}"))?,
                )
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

mod assembly;
mod curves;
mod parse;
mod roundtrip;
mod shape_analysis;
mod srr_composer;
mod srrwt;
mod step_file;

pub use self::assembly::*;
pub use self::curves::*;
pub use self::parse::*;
pub use self::roundtrip::*;
pub use self::shape_analysis::*;
pub use self::srr_composer::*;
pub use self::srrwt::*;
pub use self::step_file::*;
