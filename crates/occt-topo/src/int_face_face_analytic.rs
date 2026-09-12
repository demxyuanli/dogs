//! Analytic FaceFace pair intersections (IntAna / intpatch closed forms).
use std::sync::Arc;

use occt_core::gp::{GpAx1, GpDir, GpLin, GpPnt, GpPnt2d, GpVec};
use occt_core::precision::ANGULAR;
use occt_geom::intana::{
    quadric_quadric_cone_cone, quadric_quadric_cylinder_cone, quadric_quadric_cylinder_cylinder,
    quadric_quadric_sphere_cone, quadric_quadric_plane_cone, quadric_quadric_cylinder_sphere,
    QuadricIntersection,
};
use occt_geom::{
    Curve, GeomCircle, GeomEllipse, GeomHyperbola, GeomLine, GeomParabola, Surface,
};

use crate::brep_surface::{classify_surface, SurfaceKind};
use crate::intpatch;
use crate::inttools_data::{CurveKind, IntRange};
use crate::pcurve_full::torus_params;
use crate::shape::Face;

use super::{
    cone_from_surface, curve_range, cylinder_from_surface, cylinder_params, lin2d_through,
    line_in_uv_rect, plane_cylinder_kind, plane_torus_circles, sphere_from_surface, FaceFace,
    FaceFaceCurve,
};

impl FaceFace {
    /// Plane ∩ plane — an exact intersection line trimmed to the overlap of
    /// the two faces' UV domains (mirrors `PerformPlanes`).
    pub(crate) fn plane_plane(
        &mut self,
        sa: &dyn Surface,
        sb: &dyn Surface,
        fa: &Face,
        fb: &Face,
        tol: f64,
    ) -> Result<Vec<FaceFaceCurve>, String> {
        let pa = intpatch::plane_from_surface(sa).ok_or("FaceFace: plane extraction (a)")?;
        let pb = intpatch::plane_from_surface(sb).ok_or("FaceFace: plane extraction (b)")?;
        let n1 = GpVec::from_xyz(pa.axis().direction().xyz());
        let n2 = GpVec::from_xyz(pb.axis().direction().xyz());
        if n1.cross_magnitude(&n2) <= tol {
            // Parallel planes: coincident → tangent, otherwise no intersection.
            let dd = n1
                .normalized()
                .dot(&GpVec::from_pnts(&pa.location(), &pb.location()))
                .abs();
            if dd <= tol {
                self.tangent_faces = true;
            }
            return Ok(Vec::new());
        }
        let (origin, dir) = crate::face_face::plane_plane_intersection(&pa, &pb)
            .ok_or("FaceFace: plane-plane intersection failed")?;
        match self.trim_plane_plane_line(origin, dir, sa, sb, fa, fb, tol)? {
            Some(c) => Ok(vec![c]),
            None => Ok(Vec::new()),
        }
    }

    /// Trim the infinite plane/plane intersection line to the segment lying on
    /// both bounded faces, by projecting the line into each face's UV domain
    /// and clipping against the face UV bounds.
    pub(crate) fn trim_plane_plane_line(
        &self,
        origin: GpPnt,
        dir: GpVec,
        sa: &dyn Surface,
        sb: &dyn Surface,
        fa: &Face,
        fb: &Face,
        tol: f64,
    ) -> Result<Option<FaceFaceCurve>, String> {
        let dir_u = dir.normalized();
        let p0 = origin;
        let p1 = origin.translated_vec(&dir_u);
        let proj_a0 = intpatch::project_params(sa, &p0);
        let proj_a1 = intpatch::project_params(sa, &p1);
        let proj_b0 = intpatch::project_params(sb, &p0);
        let proj_b1 = intpatch::project_params(sb, &p1);
        let a0 = GpPnt2d::new(proj_a0.0, proj_a0.1);
        let a1 = GpPnt2d::new(proj_a1.0, proj_a1.1);
        let b0 = GpPnt2d::new(proj_b0.0, proj_b0.1);
        let b1 = GpPnt2d::new(proj_b1.0, proj_b1.1);
        let ba = crate::wireframe::face_uv_bounds(fa, sa);
        let bb = crate::wireframe::face_uv_bounds(fb, sb);

        // For each face, clip the projected line against the face UV rectangle.
        // The 3D line is p(t) = origin + t·dir (unit dir); the projected 2D
        // line maps p(t) → q0 + t·(q1 − q0), so a 2D parameter interval
        // [s0, s1] corresponds to the 3D interval [s0/k, s1/k], k = |q1 − q0|.
        let mut t_lo = f64::NEG_INFINITY;
        let mut t_hi = f64::INFINITY;
        for (q0, q1, bounds) in [(a0, a1, ba), (b0, b1, bb)] {
            let Some(lin) = lin2d_through(q0, q1) else { continue };
            let k = q0.distance(&q1);
            if k < 1e-12 {
                continue;
            }
            match line_in_uv_rect(&lin, bounds) {
                Some((s0, s1)) => {
                    t_lo = t_lo.max(s0 / k);
                    t_hi = t_hi.min(s1 / k);
                }
                None => return Ok(None),
            }
        }
        if !t_lo.is_finite() || !t_hi.is_finite() || t_hi - t_lo <= tol {
            return Ok(None);
        }

        let lin = GpLin::from_pnt_dir(origin, GpDir::from_vec(&dir_u).map_err(|e| e.to_string())?);
        let curve: Arc<dyn Curve> = Arc::new(GeomLine::new(lin));
        let range = IntRange::new_unchecked(t_lo, t_hi);
        Ok(Some(self.finish_curve(curve, CurveKind::Line, range, fa, fb)))
    }

    /// Plane ∩ sphere — an exact circle (or no curve when the plane misses the
    /// sphere), falling back to the tracer when the closed form is empty.
    pub(crate) fn plane_sphere(
        &mut self,
        sa: &dyn Surface,
        sb: &dyn Surface,
        fa: &Face,
        fb: &Face,
        tol: f64,
    ) -> Result<Vec<FaceFaceCurve>, String> {
        let is_plane_a = classify_surface(sa) == SurfaceKind::Plane;
        let (pln, (c, r)) = if is_plane_a {
            (
                intpatch::plane_from_surface(sa).ok_or("FaceFace: plane extraction")?,
                intpatch::sphere_params(sb).ok_or("FaceFace: sphere extraction")?,
            )
        } else {
            (
                intpatch::plane_from_surface(sb).ok_or("FaceFace: plane extraction")?,
                intpatch::sphere_params(sa).ok_or("FaceFace: sphere extraction")?,
            )
        };
        let mut out = Vec::new();
        if let Some(ic) = intpatch::intersect_plane_sphere(&pln, c, r) {
            out.push(self.curve_from_ic(ic, CurveKind::Circle, fa, fb));
        }
        if out.is_empty() {
            out = self.general(sa, sb, fa, fb, tol)?;
        }
        Ok(out)
    }

    /// Sphere ∩ sphere — an exact circle (or no curve when disjoint/nested),
    /// falling back to the tracer when the closed form is empty.
    pub(crate) fn sphere_sphere(
        &mut self,
        sa: &dyn Surface,
        sb: &dyn Surface,
        fa: &Face,
        fb: &Face,
        tol: f64,
    ) -> Result<Vec<FaceFaceCurve>, String> {
        let (c1, r1) = intpatch::sphere_params(sa).ok_or("FaceFace: sphere extraction (a)")?;
        let (c2, r2) = intpatch::sphere_params(sb).ok_or("FaceFace: sphere extraction (b)")?;
        let mut out = Vec::new();
        if let Some(ic) = intpatch::intersect_sphere_sphere(c1, r1, c2, r2) {
            out.push(self.curve_from_ic(ic, CurveKind::Circle, fa, fb));
        }
        if out.is_empty() {
            out = self.general(sa, sb, fa, fb, tol)?;
        }
        Ok(out)
    }

    /// Plane ∩ cylinder — an exact circle / generatrix line / ellipse, falling
    /// back to the tracer when the closed form is empty.
    pub(crate) fn plane_cylinder(
        &mut self,
        sa: &dyn Surface,
        sb: &dyn Surface,
        fa: &Face,
        fb: &Face,
        tol: f64,
    ) -> Result<Vec<FaceFaceCurve>, String> {
        let is_plane_a = classify_surface(sa) == SurfaceKind::Plane;
        let (pln, cyl) = if is_plane_a {
            (
                intpatch::plane_from_surface(sa).ok_or("FaceFace: plane extraction")?,
                sb,
            )
        } else {
            (
                intpatch::plane_from_surface(sb).ok_or("FaceFace: plane extraction")?,
                sa,
            )
        };
        let (center, ax, rad) = cylinder_params(cyl).ok_or("FaceFace: cylinder extraction")?;
        let ax1 = GpAx1::new(center, GpDir::from_vec(&ax).map_err(|e| e.to_string())?);
        let kind = plane_cylinder_kind(&pln, &ax1);
        let mut out = Vec::new();
        if let Some(ic) = intpatch::intersect_plane_cylinder(&pln, &ax1, rad) {
            out.push(self.curve_from_ic(ic, kind, fa, fb));
        }
        if out.is_empty() {
            out = self.general(sa, sb, fa, fb, tol)?;
        }
        Ok(out)
    }

    /// Plane ∩ cone — the exact conic section (circle / ellipse / parabola /
    /// hyperbola / two generatrix lines) via `IntAna_QuadQuadGeo::Perform
    /// (gp_Pln, gp_Cone)` (intana::quadric_quadric_plane_cone), falling back to
    /// the tracer when the closed form yields nothing (tangent apex point,
    /// degenerate or disjoint).
    pub(crate) fn plane_cone(
        &mut self,
        sa: &dyn Surface,
        sb: &dyn Surface,
        fa: &Face,
        fb: &Face,
        tol: f64,
    ) -> Result<Vec<FaceFaceCurve>, String> {
        let is_plane_a = classify_surface(sa) == SurfaceKind::Plane;
        let (pln, cone) = if is_plane_a {
            (
                intpatch::plane_from_surface(sa).ok_or("FaceFace: plane extraction")?,
                cone_from_surface(sb).ok_or("FaceFace: cone extraction")?,
            )
        } else {
            (
                intpatch::plane_from_surface(sb).ok_or("FaceFace: plane extraction")?,
                cone_from_surface(sa).ok_or("FaceFace: cone extraction")?,
            )
        };
        let qi = quadric_quadric_plane_cone(&pln, &cone, 1e-12, 1e-7);
        let out = self.conics_to_curves(qi, fa, fb);
        if out.is_empty() {
            self.general(sa, sb, fa, fb, tol)
        } else {
            Ok(out)
        }
    }

    /// Plane ∩ torus — the exact circles of `IntAna_QuadQuadGeo::Perform
    /// (gp_Pln, gp_Torus)` (up to two: a perpendicular cut at `major ± dt`, or
    /// the two `minor` circles when the axis lies in the plane), falling back
    /// to the tracer for an oblique section.
    pub(crate) fn plane_torus(
        &mut self,
        sa: &dyn Surface,
        sb: &dyn Surface,
        fa: &Face,
        fb: &Face,
        tol: f64,
    ) -> Result<Vec<FaceFaceCurve>, String> {
        let is_plane_a = classify_surface(sa) == SurfaceKind::Plane;
        let (pln, (center, ax, major, minor)) = if is_plane_a {
            (
                intpatch::plane_from_surface(sa).ok_or("FaceFace: plane extraction")?,
                torus_params(sb).ok_or("FaceFace: torus extraction")?,
            )
        } else {
            (
                intpatch::plane_from_surface(sb).ok_or("FaceFace: plane extraction")?,
                torus_params(sa).ok_or("FaceFace: torus extraction")?,
            )
        };
        let mut out = Vec::new();
        if let Some(circs) = plane_torus_circles(&pln, center, &ax, major, minor, tol) {
            for c in circs {
                out.push(self.conic_curve(Arc::new(GeomCircle::new(c)), CurveKind::Circle, fa, fb));
            }
        }
        if out.is_empty() {
            self.general(sa, sb, fa, fb, tol)
        } else {
            Ok(out)
        }
    }

    /// Cylinder ∩ cylinder — generatrix lines / `Same` / two ellipses /
    /// tangent via `IntAna_QuadQuadGeo::Perform(gp_Cylinder, gp_Cylinder)`,
    /// falling back to the tracer when there is no analytic closed form.
    pub(crate) fn cylinder_cylinder(
        &mut self,
        sa: &dyn Surface,
        sb: &dyn Surface,
        fa: &Face,
        fb: &Face,
        tol: f64,
    ) -> Result<Vec<FaceFaceCurve>, String> {
        let c1 = cylinder_from_surface(sa).ok_or("FaceFace: cylinder extraction (a)")?;
        let c2 = cylinder_from_surface(sb).ok_or("FaceFace: cylinder extraction (b)")?;
        let qi = quadric_quadric_cylinder_cylinder(&c1, &c2, tol);
        let out = self.conics_to_curves(qi, fa, fb);
        if out.is_empty() {
            self.general(sa, sb, fa, fb, tol)
        } else {
            Ok(out)
        }
    }

    /// Cylinder ∩ sphere — one or two circles when the sphere center lies on
    /// the cylinder axis (via `Perform(gp_Cylinder, gp_Sphere)`), else tracer.
    pub(crate) fn cylinder_sphere(
        &mut self,
        sa: &dyn Surface,
        sb: &dyn Surface,
        fa: &Face,
        fb: &Face,
        tol: f64,
    ) -> Result<Vec<FaceFaceCurve>, String> {
        // Sorted order: face1 (sa) is the sphere, face2 (sb) the cylinder.
        let sph = sphere_from_surface(sa).ok_or("FaceFace: sphere extraction")?;
        let cyl = cylinder_from_surface(sb).ok_or("FaceFace: cylinder extraction")?;
        let qi = quadric_quadric_cylinder_sphere(&cyl, &sph, tol);
        let out = self.conics_to_curves(qi, fa, fb);
        if out.is_empty() {
            self.general(sa, sb, fa, fb, tol)
        } else {
            Ok(out)
        }
    }

    /// Sphere ∩ cone — one or two circles when the sphere center lies on the
    /// cone axis (via `Perform(gp_Sphere, gp_Cone)`), else tracer.
    pub(crate) fn sphere_cone(
        &mut self,
        sa: &dyn Surface,
        sb: &dyn Surface,
        fa: &Face,
        fb: &Face,
        tol: f64,
    ) -> Result<Vec<FaceFaceCurve>, String> {
        // Sorted order: face1 (sa) is the sphere, face2 (sb) the cone.
        let sph = sphere_from_surface(sa).ok_or("FaceFace: sphere extraction")?;
        let cone = cone_from_surface(sb).ok_or("FaceFace: cone extraction")?;
        let qi = quadric_quadric_sphere_cone(&sph, &cone, tol);
        let out = self.conics_to_curves(qi, fa, fb);
        if out.is_empty() {
            self.general(sa, sb, fa, fb, tol)
        } else {
            Ok(out)
        }
    }

    /// Cone ∩ cylinder — two circles for coincident axes (via `Perform
    /// (gp_Cylinder, gp_Cone)`), else tracer.
    pub(crate) fn cylinder_cone(
        &mut self,
        sa: &dyn Surface,
        sb: &dyn Surface,
        fa: &Face,
        fb: &Face,
        tol: f64,
    ) -> Result<Vec<FaceFaceCurve>, String> {
        // Sorted order: face1 (sa) is the cone, face2 (sb) the cylinder.
        let cone = cone_from_surface(sa).ok_or("FaceFace: cone extraction")?;
        let cyl = cylinder_from_surface(sb).ok_or("FaceFace: cylinder extraction")?;
        let qi = quadric_quadric_cylinder_cone(&cyl, &cone, tol);
        let out = self.conics_to_curves(qi, fa, fb);
        if out.is_empty() {
            self.general(sa, sb, fa, fb, tol)
        } else {
            Ok(out)
        }
    }

    /// Cone ∩ cone — coincident-axis circles, parallel equal-angle conic, or
    /// shared-apex generatrices (via `Perform(gp_Cone, gp_Cone)`), else tracer.
    pub(crate) fn cone_cone(
        &mut self,
        sa: &dyn Surface,
        sb: &dyn Surface,
        fa: &Face,
        fb: &Face,
        tol: f64,
    ) -> Result<Vec<FaceFaceCurve>, String> {
        let c1 = cone_from_surface(sa).ok_or("FaceFace: cone extraction (a)")?;
        let c2 = cone_from_surface(sb).ok_or("FaceFace: cone extraction (b)")?;
        let qi = quadric_quadric_cone_cone(&c1, &c2, ANGULAR, tol);
        let out = self.conics_to_curves(qi, fa, fb);
        if out.is_empty() {
            self.general(sa, sb, fa, fb, tol)
        } else {
            Ok(out)
        }
    }

    /// Wrap an `intana` conic (analytic section) into a `FaceFaceCurve`.
    pub(crate) fn conic_curve(
        &self,
        curve: Arc<dyn Curve>,
        kind: CurveKind,
        fa: &Face,
        fb: &Face,
    ) -> FaceFaceCurve {
        let range = curve_range(curve.as_ref());
        self.finish_curve(curve, kind, range, fa, fb)
    }

    /// Convert a `QuadricIntersection` (exact conic section) into section
    /// curves. A tangent `Point`, `Same` and `None` produce no curve — the
    /// caller falls back to the tracer.
    pub(crate) fn conics_to_curves(&self, qi: QuadricIntersection, fa: &Face, fb: &Face) -> Vec<FaceFaceCurve> {
        use QuadricIntersection::*;
        match qi {
            Line(l) => vec![self.conic_curve(Arc::new(GeomLine::new(l)), CurveKind::Line, fa, fb)],
            TwoLines(l1, l2) => vec![
                self.conic_curve(Arc::new(GeomLine::new(l1)), CurveKind::Line, fa, fb),
                self.conic_curve(Arc::new(GeomLine::new(l2)), CurveKind::Line, fa, fb),
            ],
            Circle(c) => vec![self.conic_curve(Arc::new(GeomCircle::new(c)), CurveKind::Circle, fa, fb)],
            TwoCircles(c1, c2) => vec![
                self.conic_curve(Arc::new(GeomCircle::new(c1)), CurveKind::Circle, fa, fb),
                self.conic_curve(Arc::new(GeomCircle::new(c2)), CurveKind::Circle, fa, fb),
            ],
            Ellipse(e) => vec![self.conic_curve(Arc::new(GeomEllipse::new(e)), CurveKind::Ellipse, fa, fb)],
            TwoEllipses(e1, e2) => vec![
                self.conic_curve(Arc::new(GeomEllipse::new(e1)), CurveKind::Ellipse, fa, fb),
                self.conic_curve(Arc::new(GeomEllipse::new(e2)), CurveKind::Ellipse, fa, fb),
            ],
            Parabola(p) => vec![self.conic_curve(Arc::new(GeomParabola::new(p)), CurveKind::Parabola, fa, fb)],
            Hyperbola(h) => vec![self.conic_curve(Arc::new(GeomHyperbola::new(h)), CurveKind::Hyperbola, fa, fb)],
            Point(_) | Same | None => Vec::new(),
        }
    }
}
