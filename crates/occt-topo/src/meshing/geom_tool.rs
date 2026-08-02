//! Port of OCCT `BRepMesh_GeomTool` — geometric evaluation for tessellation.
//!
//! Source: `src/ModelingAlgorithms/TKMesh/BRepMesh/BRepMesh_GeomTool.{hxx,cxx}`.
//!
//! The OCCT class discretizes a geometric curve (or an iso-curve of a surface)
//! under linear + angular deflection, holds the resulting polyline, and exposes
//! static helpers for surface point/tangent/normal evaluation, first/second
//! derivatives, and 2-D segment intersection tests. This port keeps the same
//! split: [`GeomTool`] holds a discretized polyline, and the associated static
//! methods evaluate surface/curve geometry. Surfaces/curves are passed as
//! `&dyn Surface` / `&dyn Curve`, which accepts the crate's `Arc<dyn Surface>` /
//! `Arc<dyn Curve>` handles via deref coercion.

use occt_core::gcpnts::{CurveDeriv, CurveSample, TangentialDeflection};
use occt_core::gp::{GpDir, GpPnt, GpPnt2d, GpVec, GpXY};
use occt_core::precision::{ANGULAR, PCONFUSION, RESOLUTION, SQUARE_CONFUSION};
use occt_geom::{Curve, Surface};

/// Iso-curve type. Source: `GeomAbs_IsoType`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IsoType {
    /// Fixed `U`, varying `V`.
    U,
    /// Fixed `V`, varying `U`.
    V,
}

/// Status of a 2-D segment intersection check. Source: `BRepMesh_GeomTool::IntFlag`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntFlag {
    NoIntersection,
    Cross,
    EndPointTouch,
    PointOnSegment,
    Glued,
    Same,
}

/// Tool accumulating common geometrical functions for tessellation.
///
/// Instance methods discretize a curve or surface iso-curve; static methods
/// evaluate points, tangents, normals and derivatives on geometry.
pub struct GeomTool {
    params: Vec<f64>,
    points: Vec<GpPnt>,
    iso: Option<(IsoType, f64)>,
}

impl GeomTool {
    /// Discretizes a curve on `[first, last]` under linear (`lin_deflection`)
    /// and angular (`ang_deflection`) tolerances, with at least `min_points`
    /// points. Returns `(parameter, point)` pairs.
    pub fn new_curve(
        curve: &dyn Curve,
        first: f64,
        last: f64,
        lin_deflection: f64,
        ang_deflection: f64,
        min_points: usize,
    ) -> Result<Self, String> {
        if !(first.is_finite() && last.is_finite() && last >= first) {
            return Err("GeomTool::new_curve: invalid parameter range".to_string());
        }
        let sampled = Self::discretize_curve(
            curve,
            first,
            last,
            lin_deflection,
            ang_deflection,
            min_points,
        );
        let params = sampled.iter().map(|&(t, _)| t).collect();
        let points = sampled.into_iter().map(|(_, p)| p).collect();
        Ok(Self { params, points, iso: None })
    }

    /// Discretizes an iso-curve of a surface on `[first, last]`.
    /// For [`IsoType::U`] `first`/`last` are `V` bounds; for [`IsoType::V`] they
    /// are `U` bounds. Returns `(parameter, point)` pairs plus the reconstructed
    /// `(U, V)` surface parameter for every point.
    pub fn new_iso_curve(
        surface: &dyn Surface,
        iso_type: IsoType,
        iso_param: f64,
        first: f64,
        last: f64,
        lin_deflection: f64,
        ang_deflection: f64,
        min_points: usize,
    ) -> Result<Self, String> {
        if !(first.is_finite() && last.is_finite() && last >= first) {
            return Err("GeomTool::new_iso_curve: invalid parameter range".to_string());
        }
        let sampled = Self::discretize_iso_curve(
            surface,
            iso_type,
            iso_param,
            first,
            last,
            lin_deflection,
            ang_deflection,
            min_points,
        );
        let params = sampled.iter().map(|&(t, _, _)| t).collect();
        let points = sampled.into_iter().map(|(_, p, _)| p).collect();
        Ok(Self { params, points, iso: Some((iso_type, iso_param)) })
    }

    /// Number of discretization points.
    pub fn nb_points(&self) -> usize {
        self.params.len()
    }

    /// Parameter of the `index`-th point (0-based).
    pub fn parameter(&self, index: usize) -> Option<f64> {
        self.params.get(index).copied()
    }

    /// Point of the `index`-th discretization point (0-based).
    pub fn point(&self, index: usize) -> Option<GpPnt> {
        self.points.get(index).copied()
    }

    /// `(parameter, point)` of the `index`-th point (0-based).
    pub fn value(&self, index: usize) -> Option<(f64, GpPnt)> {
        match (self.params.get(index), self.points.get(index)) {
            (Some(&t), Some(&p)) => Some((t, p)),
            _ => None,
        }
    }

    /// `(parameter, point, uv)` of the `index`-th point (0-based).
    ///
    /// The UV pair is only available for iso-curve discretizations; for a
    /// plain curve this returns `None`.
    pub fn value_uv(&self, index: usize) -> Option<(f64, GpPnt, GpPnt2d)> {
        let (t, p) = self.value(index)?;
        match self.iso {
            Some((IsoType::U, iso_param)) => Some((t, p, GpPnt2d::new(iso_param, t))),
            Some((IsoType::V, iso_param)) => Some((t, p, GpPnt2d::new(t, iso_param))),
            None => None,
        }
    }

    /// Adds a point at parameter `param`, keeping the parameter list sorted.
    /// When `replace` is true and a point already exists within parametric
    /// tolerance, that point is replaced instead. Returns the 0-based index of
    /// the new/replaced point.
    pub fn add_point(&mut self, point: GpPnt, param: f64, replace: bool) -> usize {
        let prec = PCONFUSION;
        let mut idx = 0;
        while idx < self.params.len() && self.params[idx] < param {
            idx += 1;
        }
        if replace {
            for j in [idx.wrapping_sub(1), idx] {
                if let Some(&p) = self.params.get(j) {
                    if (p - param).abs() < prec {
                        self.params[j] = param;
                        self.points[j] = point;
                        return j;
                    }
                }
            }
        }
        self.params.insert(idx, param);
        self.points.insert(idx, point);
        idx
    }

    // ------------------------------------------------------------------
    // Static geometry evaluation
    // ------------------------------------------------------------------

    /// Point of `surface` at parameters `(u, v)`. Source: `ValueOnSurface`
    pub fn value_on_surface(surface: &dyn Surface, u: f64, v: f64) -> GpPnt {
        surface.d0(u, v)
    }

    /// Point and both first partials `(du, dv)` of `surface` at `(u, v)`.
    ///
    /// Uses `Surface::d1` when it is non-degenerate; otherwise falls back to
    /// central finite differences of `d0` (the elementary-surface
    /// implementations in `occt-geom` leave `d1` zero for cylinder/cone/sphere/
    /// torus). Source: `TangentOnSurface` / `ValueAndTangents`
    pub fn tangent_on_surface(surface: &dyn Surface, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
        let (p, du, dv) = surface.d1(u, v);
        let (du, dv) = if du.magnitude() <= 1e-12 || dv.magnitude() <= 1e-12 {
            let (fdu, fdv) = fd_partials(surface, u, v);
            let du = if du.magnitude() <= 1e-12 { fdu } else { du };
            let dv = if dv.magnitude() <= 1e-12 { fdv } else { dv };
            (du, dv)
        } else {
            (du, dv)
        };
        (p, du, dv)
    }

    /// Point and both first partials of `surface` at `(u, v)`.
    /// Alias of [`GeomTool::tangent_on_surface`] matching the OCCT name.
    pub fn value_and_tangents(surface: &dyn Surface, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
        Self::tangent_on_surface(surface, u, v)
    }

    /// Point and outward-oriented unit normal of `surface` at `(u, v)`.
    ///
    /// Degenerate points (parallel first partials, e.g. a sphere pole) fall back
    /// to a second-derivative candidate search; if no stable normal is found an
    /// error is returned, mirroring OCCT's "FALSE if the normal cannot be
    /// computed". Source: `Normal`
    pub fn normal_on_surface(surface: &dyn Surface, u: f64, v: f64) -> Result<(GpPnt, GpDir), String> {
        let (p, du, dv) = Self::tangent_on_surface(surface, u, v);
        let n = du.crossed(&dv);
        let tol = ANGULAR * du.magnitude().max(1e-12) * dv.magnitude().max(1e-12);
        if n.magnitude() > tol {
            return Ok((p, GpDir::from_vec(&n).map_err(|e| e.to_string())?));
        }

        // Degenerate: try combinations involving second derivatives.
        let (duu, dvv, duv) = fd_d2(surface, u, v);
        let candidates = [
            du.crossed(&dvv),
            duu.crossed(&dv),
            duu.crossed(&dvv),
            duu.crossed(&duv),
            duv.crossed(&dv),
        ];
        let mut best = GpVec::zero();
        for c in candidates {
            if c.square_magnitude() > best.square_magnitude() {
                best = c;
            }
        }
        if best.magnitude() > tol {
            return Ok((p, GpDir::from_vec(&best).map_err(|e| e.to_string())?));
        }
        Err("GeomTool::normal_on_surface: cannot compute surface normal (degenerate point)".to_string())
    }

    /// Point and first derivative of `curve` at `u`. Source: `EvaluateCubicalDerivative`
    pub fn evaluate_curve_d1(curve: &dyn Curve, u: f64) -> (GpPnt, GpVec) {
        curve.d1(u)
    }

    /// Point, first and second derivatives of `curve` at `u`.
    pub fn evaluate_curve_d2(curve: &dyn Curve, u: f64) -> (GpPnt, GpVec, GpVec) {
        curve.d2(u)
    }

    /// Point and first partials of `surface` at `(u, v)`.
    pub fn evaluate_surface_d1(surface: &dyn Surface, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
        Self::tangent_on_surface(surface, u, v)
    }

    /// Point, first and second partials `(du, dv, duu, dvv, duv)` of `surface`
    /// at `(u, v)`. Second partials come from finite differences of `d0` since
    /// `Surface` exposes no `d2`.
    pub fn evaluate_surface_d2(
        surface: &dyn Surface,
        u: f64,
        v: f64,
    ) -> (GpPnt, GpVec, GpVec, GpVec, GpVec, GpVec) {
        let (p, du, dv) = Self::tangent_on_surface(surface, u, v);
        let (duu, dvv, duv) = fd_d2(surface, u, v);
        (p, du, dv, duu, dvv, duv)
    }

    /// Square distance of `mid` to the chord `first..last`. When the chord is
    /// degenerate the square distance to `first` is returned. Source:
    /// `SquareDeflectionOfSegment`
    pub fn square_deflection_of_segment(first: &GpPnt, last: &GpPnt, mid: &GpPnt) -> f64 {
        if first.square_distance(last) > SQUARE_CONFUSION {
            let v = GpVec::from_pnts(first, last);
            let vm = GpVec::from_pnts(first, mid);
            let t = vm.dot(&v) / v.square_magnitude();
            let proj = GpPnt::new(
                first.x() + t * v.x(),
                first.y() + t * v.y(),
                first.z() + t * v.z(),
            );
            proj.square_distance(mid)
        } else {
            first.square_distance(mid)
        }
    }

    /// Intersection of two infinite lines (each given by two points), returns
    /// the status, the intersection point and the parameters on both lines.
    /// Source: `IntLinLin`
    pub fn int_lin_lin(
        start1: &GpXY,
        end1: &GpXY,
        start2: &GpXY,
        end2: &GpXY,
    ) -> (IntFlag, GpXY, [f64; 2]) {
        let v1 = end1.subtracted(start1);
        let v2 = end2.subtracted(start2);
        let vo1o2 = start2.subtracted(start1);

        let cross_d1d2 = v1.crossed(&v2);
        let cross_d1d3 = vo1o2.crossed(&v2);

        let prec = RESOLUTION;
        if cross_d1d2.abs() < prec {
            return if cross_d1d3.abs() < prec {
                (IntFlag::Same, GpXY::zero(), [0.0; 2])
            } else {
                (IntFlag::NoIntersection, GpXY::zero(), [0.0; 2])
            };
        }

        let mut param = [0.0f64; 2];
        param[0] = cross_d1d3 / cross_d1d2;
        let int_pnt = start1.added(&v1.multiplied(param[0]));

        let cross_d2d3 = vo1o2.crossed(&v1);
        param[1] = cross_d2d3 / cross_d1d2;

        (IntFlag::Cross, int_pnt, param)
    }

    /// Intersection of two segments, honoring the endpoint-touch and
    /// point-on-segment considerations. Source: `IntSegSeg`
    pub fn int_seg_seg(
        start1: &GpXY,
        end1: &GpXY,
        start2: &GpXY,
        end2: &GpXY,
        consider_end_point_touch: bool,
        consider_point_on_segment: bool,
    ) -> (IntFlag, GpXY) {
        let hash = [
            classify_point(start1, end1, start2),
            classify_point(start1, end1, end2),
            classify_point(start2, end2, start1),
            classify_point(start2, end2, end1),
        ];
        let pos_hash = hash.iter().sum::<i32>();

        // Shared vertex case.
        if hash[0] < 0 || hash[1] < 0 {
            if pos_hash == -1 {
                return (IntFlag::Glued, GpXY::zero());
            }
            return if consider_end_point_touch {
                (IntFlag::EndPointTouch, GpXY::zero())
            } else {
                (IntFlag::NoIntersection, GpXY::zero())
            };
        }

        if pos_hash == 1 {
            if consider_point_on_segment {
                let pnt = if hash[0] == 1 {
                    *start1
                } else if hash[1] == 1 {
                    *end1
                } else if hash[2] == 1 {
                    *start2
                } else {
                    *end2
                };
                return (IntFlag::PointOnSegment, pnt);
            }
            return (IntFlag::NoIntersection, GpXY::zero());
        } else if pos_hash == 2 {
            return (IntFlag::Glued, GpXY::zero());
        }

        let (flag, pnt, param) = Self::int_lin_lin(start1, end1, start2, end2);
        if flag == IntFlag::NoIntersection {
            return (IntFlag::NoIntersection, GpXY::zero());
        }
        if flag == IntFlag::Same {
            if pos_hash < -2 {
                return (IntFlag::Same, GpXY::zero());
            } else if pos_hash == -1 {
                return (IntFlag::Glued, GpXY::zero());
            }
            return (IntFlag::NoIntersection, GpXY::zero());
        }

        // Cross: intersection must lie within both segments.
        let prec = PCONFUSION;
        let end_prec = 1.0 - prec;
        if param.iter().any(|&p| p < prec || p > end_prec) {
            return (IntFlag::NoIntersection, GpXY::zero());
        }
        (IntFlag::Cross, pnt)
    }

    /// Discretizes `curve` on `[first, last]` with adaptive tangential
    /// deflection. Returns `(parameter, point)` pairs.
    pub fn discretize_curve(
        curve: &dyn Curve,
        first: f64,
        last: f64,
        lin_deflection: f64,
        ang_deflection: f64,
        min_points: usize,
    ) -> Vec<(f64, GpPnt)> {
        let adapter = CurveAdapter { curve };
        let td = TangentialDeflection::from_curve_with_deriv(
            &adapter,
            first,
            last,
            lin_deflection,
            ang_deflection,
        );
        enforce_min_points(&td.params, &td.points, min_points, &adapter)
    }

    /// Discretizes an iso-curve of `surface`. Returns `(parameter, point, uv)`.
    pub fn discretize_iso_curve(
        surface: &dyn Surface,
        iso_type: IsoType,
        iso_param: f64,
        first: f64,
        last: f64,
        lin_deflection: f64,
        ang_deflection: f64,
        min_points: usize,
    ) -> Vec<(f64, GpPnt, GpPnt2d)> {
        let adapter = IsoCurveAdapter { surface, iso_type, iso_param };
        let td = TangentialDeflection::from_curve_with_deriv(
            &adapter,
            first,
            last,
            lin_deflection,
            ang_deflection,
        );
        let pairs = enforce_min_points(&td.params, &td.points, min_points, &adapter);
        pairs
            .into_iter()
            .map(|(t, p)| {
                let uv = match iso_type {
                    IsoType::U => GpPnt2d::new(iso_param, t),
                    IsoType::V => GpPnt2d::new(t, iso_param),
                };
                (t, p, uv)
            })
            .collect()
    }
}

/// Enforces the minimum-points contract: if the adaptive sampler produced fewer
/// than `min_points`, re-samples uniformly. Otherwise keeps the adaptive points.
fn enforce_min_points<A: CurveSample>(
    params: &[f64],
    points: &[GpPnt],
    min_points: usize,
    adapter: &A,
) -> Vec<(f64, GpPnt)> {
    if params.len() >= min_points.max(2) {
        return params.iter().zip(points.iter()).map(|(&t, &p)| (t, p)).collect();
    }
    let n = min_points.max(2);
    let (a, b) = (params[0], *params.last().unwrap());
    (0..n)
        .map(|i| {
            let t = a + (b - a) * i as f64 / (n - 1) as f64;
            (t, adapter.point(t))
        })
        .collect()
}

/// `CurveSample`/`CurveDeriv` adapter over a `dyn Curve` for the gcpnts samplers.
struct CurveAdapter<'a> {
    curve: &'a dyn Curve,
}

impl CurveSample for CurveAdapter<'_> {
    fn point(&self, u: f64) -> GpPnt {
        self.curve.d0(u)
    }
}

impl CurveDeriv for CurveAdapter<'_> {
    fn tangent(&self, u: f64) -> GpVec {
        self.curve.d1(u).1
    }
}

/// `CurveSample`/`CurveDeriv` adapter over a surface iso-curve.
struct IsoCurveAdapter<'a> {
    surface: &'a dyn Surface,
    iso_type: IsoType,
    iso_param: f64,
}

impl IsoCurveAdapter<'_> {
    fn params(&self, t: f64) -> (f64, f64) {
        match self.iso_type {
            IsoType::U => (self.iso_param, t),
            IsoType::V => (t, self.iso_param),
        }
    }
}

impl CurveSample for IsoCurveAdapter<'_> {
    fn point(&self, t: f64) -> GpPnt {
        let (u, v) = self.params(t);
        self.surface.d0(u, v)
    }
}

impl CurveDeriv for IsoCurveAdapter<'_> {
    fn tangent(&self, t: f64) -> GpVec {
        let (u, v) = self.params(t);
        let (_, du, dv) = GeomTool::tangent_on_surface(self.surface, u, v);
        match self.iso_type {
            IsoType::U => dv,
            IsoType::V => du,
        }
    }
}

/// Finite-difference step for a surface, scaled to its parametric extent.
fn fd_step(surface: &dyn Surface) -> f64 {
    let (u0, u1) = surface.u_range();
    let (v0, v1) = surface.v_range();
    let du = u1 - u0;
    let dv = v1 - v0;
    let range = if du.is_finite() && dv.is_finite() {
        du.abs().max(dv.abs())
    } else {
        1.0
    };
    (range * 1e-5).max(1e-7)
}

/// Central-difference first partials of a surface from `d0` only.
fn fd_partials(surface: &dyn Surface, u: f64, v: f64) -> (GpVec, GpVec) {
    let h = fd_step(surface);
    let pu = surface.d0(u + h, v);
    let pm = surface.d0(u - h, v);
    let pv = surface.d0(u, v + h);
    let pn = surface.d0(u, v - h);
    let du = GpVec::from_xyz(&pu.coord.subtracted(&pm.coord).divided(2.0 * h));
    let dv = GpVec::from_xyz(&pv.coord.subtracted(&pn.coord).divided(2.0 * h));
    (du, dv)
}

/// Central-difference second partials `(duu, dvv, duv)` of a surface from `d0`.
fn fd_d2(surface: &dyn Surface, u: f64, v: f64) -> (GpVec, GpVec, GpVec) {
    let h = fd_step(surface);
    let h2 = h * h;
    let p = surface.d0(u, v);
    let pu = surface.d0(u + h, v);
    let pmu = surface.d0(u - h, v);
    let pv = surface.d0(u, v + h);
    let pmv = surface.d0(u, v - h);
    let puv = surface.d0(u + h, v + h);
    let pumv = surface.d0(u + h, v - h);
    let pmuv = surface.d0(u - h, v + h);
    let pmumv = surface.d0(u - h, v - h);

    let duu = pu.coord.subtracted(&p.coord.multiplied(2.0)).added(&pmu.coord).divided(h2);
    let dvv = pv.coord.subtracted(&p.coord.multiplied(2.0)).added(&pmv.coord).divided(h2);
    let duv = puv
        .coord
        .subtracted(&pumv.coord)
        .subtracted(&pmuv.coord)
        .added(&pmumv.coord)
        .divided(4.0 * h2);
    (GpVec::from_xyz(&duu), GpVec::from_xyz(&dvv), GpVec::from_xyz(&duv))
}

/// Classifies `check` against segment `p1..p2`: `1` interior, `-1` coincident
/// with an endpoint, `0` outside. Source: `BRepMesh_GeomTool::classifyPoint`
fn classify_point(p1: &GpXY, p2: &GpXY, check: &GpXY) -> i32 {
    let ap1 = p2.subtracted(p1);
    let ap2 = check.subtracted(p1);

    let prec = PCONFUSION;
    let sq_prec = prec * prec;
    let mut dist = ap1.crossed(&ap2).abs();
    if dist > prec {
        dist = (dist * dist) / ap1.square_modulus();
        if dist > sq_prec {
            return 0;
        }
    }

    let mult = ap1.multiplied_xy(&ap2);
    if mult.x() < 0.0 || mult.y() < 0.0 {
        return 0;
    }
    if ap1.square_modulus() < ap2.square_modulus() {
        return 0;
    }
    if check.is_equal(p1, prec) || check.is_equal(p2, prec) {
        return -1;
    }
    1
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;
    use std::sync::Arc;

    use occt_core::gp::{GpAx2, GpAx3, GpCirc, GpCylinder, GpPln, GpSphere};
    use occt_geom::{GeomCircle, GeomCylinder, GeomPlane, GeomSphere};

    fn plane() -> Arc<dyn Surface> {
        Arc::new(GeomPlane::new(GpPln::new(GpAx3::standard())))
    }

    fn cylinder(r: f64) -> Arc<dyn Surface> {
        Arc::new(GeomCylinder::new(GpCylinder::new(GpAx3::standard(), r).unwrap()))
    }

    fn sphere(r: f64) -> Arc<dyn Surface> {
        Arc::new(GeomSphere::new(GpSphere::new(GpAx3::standard(), r).unwrap()))
    }

    #[test]
    fn plane_point_and_tangents_are_exact() {
        let s = plane();
        // d0(u, v) = origin + u*X + v*Y = (u, v, 0)
        let p = GeomTool::value_on_surface(s.as_ref(), 1.0, 2.0);
        assert!((p.x() - 1.0).abs() < 1e-12);
        assert!((p.y() - 2.0).abs() < 1e-12);
        assert!(p.z().abs() < 1e-12);

        let (pt, du, dv) = GeomTool::tangent_on_surface(s.as_ref(), 0.5, -1.0);
        assert!((pt.x() - 0.5).abs() < 1e-12 && (pt.y() + 1.0).abs() < 1e-12);
        assert!((du.x() - 1.0).abs() < 1e-12 && du.y().abs() < 1e-12 && du.z().abs() < 1e-12);
        assert!(dv.x().abs() < 1e-12 && (dv.y() - 1.0).abs() < 1e-12 && dv.z().abs() < 1e-12);
    }

    #[test]
    fn cylinder_point_and_tangents() {
        let s = cylinder(2.0);
        // (r cos u, r sin u, v) at (0, 3) = (2, 0, 3)
        let p = GeomTool::value_on_surface(s.as_ref(), 0.0, 3.0);
        assert!((p.x() - 2.0).abs() < 1e-12);
        assert!(p.y().abs() < 1e-12);
        assert!((p.z() - 3.0).abs() < 1e-12);

        // du = (-r sin u, r cos u, 0) = (0, 2, 0); dv = (0, 0, 1)
        let (_, du, dv) = GeomTool::tangent_on_surface(s.as_ref(), 0.0, 0.0);
        assert!(du.x().abs() < 1e-4 && (du.y() - 2.0).abs() < 1e-4 && du.z().abs() < 1e-4);
        assert!(dv.x().abs() < 1e-4 && dv.y().abs() < 1e-4 && (dv.z() - 1.0).abs() < 1e-4);
    }

    #[test]
    fn sphere_point_and_tangents() {
        let s = sphere(3.0);
        // (r cos v cos u, r cos v sin u, r sin v) at (0, 0) = (3, 0, 0)
        let p = GeomTool::value_on_surface(s.as_ref(), 0.0, 0.0);
        assert!((p.x() - 3.0).abs() < 1e-12);
        assert!(p.y().abs() < 1e-12);
        assert!(p.z().abs() < 1e-12);

        // du = (0, r, 0) = (0, 3, 0); dv = (0, 0, r) = (0, 0, 3)
        let (_, du, dv) = GeomTool::tangent_on_surface(s.as_ref(), 0.0, 0.0);
        assert!(du.x().abs() < 1e-4 && (du.y() - 3.0).abs() < 1e-4 && du.z().abs() < 1e-4);
        assert!(dv.x().abs() < 1e-4 && dv.y().abs() < 1e-4 && (dv.z() - 3.0).abs() < 1e-4);
    }

    #[test]
    fn normal_on_plane_and_cylinder() {
        let s = plane();
        let (_, n) = GeomTool::normal_on_surface(s.as_ref(), 0.0, 0.0).unwrap();
        assert!((n.z() - 1.0).abs() < 1e-12);

        let c = cylinder(2.0);
        let (_, n) = GeomTool::normal_on_surface(c.as_ref(), 0.0, 0.0).unwrap();
        assert!((n.x() - 1.0).abs() < 1e-4);
        assert!(n.y().abs() < 1e-4);
        assert!(n.z().abs() < 1e-4);
    }

    #[test]
    fn circle_curve_derivatives() {
        let c: Arc<dyn Curve> = Arc::new(GeomCircle::new(GpCirc::new(GpAx2::standard(), 2.0)));
        let (p, d1) = GeomTool::evaluate_curve_d1(c.as_ref(), 0.0);
        assert!((p.x() - 2.0).abs() < 1e-12 && p.y().abs() < 1e-12);
        assert!(d1.x().abs() < 1e-12 && (d1.y() - 2.0).abs() < 1e-12);

        let (_, _, d2) = GeomTool::evaluate_curve_d2(c.as_ref(), 0.0);
        assert!((d2.x() + 2.0).abs() < 1e-12 && d2.y().abs() < 1e-12);
    }

    #[test]
    fn square_deflection_of_segment_matches_occt() {
        let a = GpPnt::new(0.0, 0.0, 0.0);
        let b = GpPnt::new(10.0, 0.0, 0.0);
        assert!(GeomTool::square_deflection_of_segment(&a, &b, &GpPnt::new(5.0, 0.0, 0.0)) < 1e-20);
        assert!((GeomTool::square_deflection_of_segment(&a, &b, &GpPnt::new(5.0, 1.0, 0.0)) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn int_seg_seg_crossing() {
        // Horizontal (0,0)-(2,0) crossed by vertical (1,-1)-(1,1) at (1, 0).
        let (flag, pnt) = GeomTool::int_seg_seg(
            &GpXY::new(0.0, 0.0),
            &GpXY::new(2.0, 0.0),
            &GpXY::new(1.0, -1.0),
            &GpXY::new(1.0, 1.0),
            true,
            true,
        );
        assert_eq!(flag, IntFlag::Cross);
        assert!((pnt.x() - 1.0).abs() < 1e-12 && (pnt.y() - 0.0).abs() < 1e-12);
    }

    #[test]
    fn int_seg_seg_disjoint() {
        let (flag, _) = GeomTool::int_seg_seg(
            &GpXY::new(0.0, 0.0),
            &GpXY::new(1.0, 0.0),
            &GpXY::new(2.0, 0.0),
            &GpXY::new(3.0, 0.0),
            true,
            true,
        );
        assert_eq!(flag, IntFlag::NoIntersection);
    }

    #[test]
    fn discretize_curve_respects_deflection() {
        let c: Arc<dyn Curve> = Arc::new(GeomCircle::new(GpCirc::new(GpAx2::standard(), 2.0)));
        let sampled = GeomTool::discretize_curve(c.as_ref(), 0.0, PI, 0.02, 0.2, 4);
        assert!(sampled.len() >= 2);
        let points: Vec<GpPnt> = sampled.iter().map(|&(_, p)| p).collect();
        let params: Vec<f64> = sampled.iter().map(|&(t, _)| t).collect();
        let err = occt_core::gcpnts::total_chord_error(&points, &|u| c.d0(u), &params);
        assert!(err < 0.03, "chord error {err}");
    }

    #[test]
    fn discretize_iso_curve_uv_reconstruction() {
        let s = cylinder(1.0);
        // V-iso at v = 2.0 is a circle of radius 1 in plane z = 2.
        let sampled = GeomTool::discretize_iso_curve(s.as_ref(), IsoType::V, 2.0, 0.0, PI, 0.01, 0.1, 8);
        assert!(sampled.len() >= 2);
        for &(t, p, uv) in sampled.iter() {
            assert!((uv.x() - t).abs() < 1e-12 && (uv.y() - 2.0).abs() < 1e-12);
            assert!((p.z() - 2.0).abs() < 1e-9);
        }
    }
}
