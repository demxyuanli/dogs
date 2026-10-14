//! `StepToGeom` placement and axis derivations: `MakeCartesianPoint`, `MakeDirection`, `MakeVectorWithMagnitude`, the axis placements, `ExtractParameter` (`cxx:2221-2317`) and `StepGeom_Axis2Placement3d::HasRefDirection`.
use super::*;

impl<'a> Resolver<'a> {
    pub(in crate::step) fn resolve_point(&self, id: usize) -> Result<GpPnt, String> {
        if let Some(p) = self.point_cache.borrow().get(&id) {
            return Ok(*p);
        }
        let rec = self.record(id)?;
        if rec.type_name != "CARTESIAN_POINT" {
            self.warn(format!("expected CARTESIAN_POINT, got {} (#{id})", rec.type_name));
            return Err(format!("expected CARTESIAN_POINT at #{id}"));
        }
        let v = parse_xyz(&rec.args[1]).map_err(|e| format!("CARTESIAN_POINT #{id}: {e}"))?;
        // `StepToGeom::MakeCartesianPoint` (`StepToGeom.cxx:1173-1186`): every
        // coordinate is multiplied by `LF = theLocalFactors.LengthFactor()`
        // (`cxx:1179`). The 2D reader (`MakeCartesianPoint2d`, `cxx:1191-1204`)
        // deliberately does NOT scale, so `resolve_point_2d` stays unchanged.
        let p = GpPnt::from_xyz(&v.multiplied(self.length_factor));
        self.point_cache.borrow_mut().insert(id, p);
        Ok(p)
    }

    pub(in crate::step) fn resolve_direction(&self, id: usize) -> Result<GpDir, String> {
        if let Some(d) = self.dir_cache.borrow().get(&id) {
            return Ok(*d);
        }
        let rec = self.record(id)?;
        if rec.type_name != "DIRECTION" {
            self.warn(format!("expected DIRECTION, got {} (#{id})", rec.type_name));
            return Err(format!("expected DIRECTION at #{id}"));
        }
        let v = parse_xyz(&rec.args[1]).map_err(|e| format!("DIRECTION #{id}: {e}"))?;
        let d = GpDir::new(v.x, v.y, v.z).map_err(|e| format!("DIRECTION #{id}: {e}"))?;
        self.dir_cache.borrow_mut().insert(id, d);
        Ok(d)
    }

    pub(in crate::step) fn resolve_vector(&self, id: usize) -> Result<GpVec, String> {
        let rec = self.record(id)?;
        if rec.type_name != "VECTOR" {
            self.warn(format!("expected VECTOR, got {} (#{id})", rec.type_name));
            return Err(format!("expected VECTOR at #{id}"));
        }
        let dir_ref = parse_ref(&rec.args[1]).ok_or("VECTOR: bad direction ref")?;
        let mag = parse_f64(&rec.args[2])?;
        let dir = self.resolve_direction(dir_ref)?;
        // `StepToGeom::MakeVectorWithMagnitude` (`StepToGeom.cxx:2569-2581`):
        // `V = D->Dir().XYZ() * SV->Magnitude() * LengthFactor` (`cxx:2577`).
        // `MakeVectorWithMagnitude2d` (`cxx:2586-2597`) does not scale, so
        // `resolve_vector_2d` stays unchanged.
        Ok(GpVec::from_xyz(
            &dir.xyz().multiplied(mag * self.length_factor),
        ))
    }

    /// `AXIS1_PLACEMENT(name, location, axis_direction)` → `GpAx1`.
    pub(in crate::step) fn resolve_axis1(&self, id: usize) -> Result<GpAx1, String> {
        let rec = self.record(id)?;
        if rec.type_name != "AXIS1_PLACEMENT" {
            self.warn(format!("expected AXIS1_PLACEMENT, got {} (#{id})", rec.type_name));
            return Err(format!("expected AXIS1_PLACEMENT at #{id}"));
        }
        let loc_ref = parse_ref(&rec.args[1]).ok_or("AXIS1: bad location ref")?;
        let dir_ref = parse_ref(&rec.args[2]).ok_or("AXIS1: bad direction ref")?;
        let loc = self.resolve_point(loc_ref)?;
        let dir = self.resolve_direction(dir_ref)?;
        Ok(GpAx1::new(loc, dir))
    }

    pub(in crate::step) fn resolve_axis2(&self, id: usize) -> Result<GpAx2, String> {
        if let Some(a) = self.axis_cache.borrow().get(&id) {
            return Ok(*a);
        }
        let rec = self.record(id)?;
        if rec.type_name != "AXIS2_PLACEMENT_3D" {
            self.warn(format!(
                "expected AXIS2_PLACEMENT_3D, got {} (#{id})",
                rec.type_name
            ));
            return Err(format!("expected AXIS2_PLACEMENT_3D at #{id}"));
        }
        let loc = self
            .resolve_point(parse_ref(&rec.args[1]).ok_or("AXIS2: bad location ref")?)?;
        let axis = match parse_ref(&rec.args[2]) {
            Some(r) => self.resolve_direction(r)?,
            None => dir_z(),
        };
        let refd = match parse_ref(&rec.args[3]) {
            Some(r) => self.resolve_direction(r)?,
            None => dir_x(),
        };
        let ax2 = GpAx2::new(loc, axis, refd).map_err(|e| format!("AXIS2_PLACEMENT_3D: {e}"))?;
        self.axis_cache.borrow_mut().insert(id, ax2);
        Ok(ax2)
    }

    /// One side of `StepToGeom::ExtractParameter` (`StepToGeom.cxx:2221-2317`):
    /// the master representation decides which select is authoritative; a
    /// parameter select becomes `Shift + Factor * value` (`cxx:2238`, `:2286`),
    /// a point select is projected onto the curve (`cxx:2249-2251`,
    /// `ShapeAnalysis_Curve::Project`; the port uses its Extrema-based
    /// projection, `int_tools_vertex_line::project_point_on_curve_param`).
    pub(super) fn extract_parameter(
        &self,
        curve: &dyn Curve,
        params: &[f64],
        points: &[usize],
        master_rep: i32,
        fact: f64,
        shift: f64,
    ) -> Result<Option<f64>, String> {
        if master_rep == 2 {
            if let Some(p) = params.first() {
                return Ok(Some(shift + fact * p));
            }
        } else if master_rep == 1 {
            if let Some(&r) = points.first() {
                let p = self.resolve_point(r)?;
                return Ok(crate::int_tools_vertex_line::project_point_on_curve_param(curve, &p));
            }
        }
        // `cxx:2278-2314`: with an unspecified master representation the parameter
        // is preferred, and the point is projected only when no parameter exists.
        if let Some(p) = params.first() {
            return Ok(Some(shift + fact * p));
        }
        if let Some(&r) = points.first() {
            let p = self.resolve_point(r)?;
            return Ok(crate::int_tools_vertex_line::project_point_on_curve_param(curve, &p));
        }
        Ok(None)
    }

    /// `StepGeom_Axis2Placement3d::HasRefDirection` (`cxx:2404`): the STEP record
    /// carries `$` in the reference-direction slot when it is absent.
    pub(in crate::step) fn axis2_has_ref_direction(&self, id: usize) -> Result<bool, String> {
        let rec = self.record(id)?;
        if rec.type_name != "AXIS2_PLACEMENT_3D" {
            return Ok(true);
        }
        Ok(match rec.args.get(3).map(|s| s.trim()) {
            None | Some("") | Some("$") | Some(".F.") => false,
            _ => true,
        })
    }

}
