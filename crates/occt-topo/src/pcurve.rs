//! Edge→face UV curve (pcurve) construction.
//!
//! Projects an edge's 3D curve into the parameter domain `(u, v)` of a face's
//! surface. Analytic cases get exact 2D curves: a line edge on a plane or a
//! cylinder generatrix becomes a `Geom2dLine`; a circle edge in a plane
//! parallel to the face becomes a `Geom2dCircle`. Every other case samples
//! `brep_surface::edge_pcurve_on_face` and fits a polyline B-spline whose
//! parameter range matches the edge range.
//!
//! Source: `BOPTools_AlgoTools2D::MakePCurveOnFace` (TKBO) — the OCCT routine
//! that builds `Geom2d_Curve` pcurves by projection onto the face surface.

use std::sync::Arc;

use occt_core::gp::{
    GpAx1, GpAx2d, GpAx22d, GpAx3, GpCirc2d, GpDir, GpDir2d, GpElips2d, GpPln, GpPnt, GpPnt2d, GpVec,
    GpVec2d,
};
use occt_geom::projlib;
use occt_geom::{Curve, Surface};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::{Geom2dBSplineCurve, Geom2dCircle, Geom2dEllipse, Geom2dLine};

use crate::brep_surface::{edge_pcurve_on_face, is_planar};
use crate::shape::{Edge, Face};
use crate::tgeometry::GeometryRegistry;

/// `GeomProjLib::ProjectOnPlane` KeepParam=true as a 2d curve: UV of the
/// normal projection of `C3D(t)`. Used when `C3D` is `Geom_TrimmedCurve`
/// (`GeomProjLib.cxx:339-343`) so same-t matches the remapped `[0,1]` domain
/// instead of the basis `gp_Circ` / BSpline poles.
struct PlaneKeepParam2d {
    curve: Arc<dyn Curve>,
    pln: GpPln,
}

impl PlaneKeepParam2d {
    fn uv_of(&self, p: &GpPnt) -> GpPnt2d {
        let q = projlib::project_pnt_on_plane(&self.pln, p);
        projlib::eval_pln_pnt2d(&self.pln, &q)
    }

    fn duv_of(&self, v: &GpVec) -> GpVec2d {
        let z = self.pln.pos.direction();
        let vz = v.x() * z.x() + v.y() * z.y() + v.z() * z.z();
        let tx = v.x() - vz * z.x();
        let ty = v.y() - vz * z.y();
        let tz = v.z() - vz * z.z();
        let x = self.pln.pos.x_direction();
        let y = self.pln.pos.y_direction();
        GpVec2d::new(
            tx * x.x() + ty * x.y() + tz * x.z(),
            tx * y.x() + ty * y.y() + tz * y.z(),
        )
    }
}

impl Curve2d for PlaneKeepParam2d {
    fn d0(&self, u: f64) -> GpPnt2d {
        self.uv_of(&self.curve.d0(u))
    }
    fn d1(&self, u: f64) -> (GpPnt2d, GpVec2d) {
        let (p, v) = self.curve.d1(u);
        (self.uv_of(&p), self.duv_of(&v))
    }
    fn d2(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
        let (p, v1, v2) = self.curve.d2(u);
        (self.uv_of(&p), self.duv_of(&v1), self.duv_of(&v2))
    }
    fn first_parameter(&self) -> f64 {
        self.curve.first_parameter()
    }
    fn last_parameter(&self) -> f64 {
        self.curve.last_parameter()
    }
    fn continuity(&self) -> u8 {
        self.curve.continuity()
    }
    fn transform(&mut self, _t: &occt_core::gp::GpTrsf2d) {}
    fn reverse(&mut self) {
        let mut c = self.curve.clone_dyn();
        c.reverse();
        self.curve = Arc::from(c);
    }
    // Kind delegation to the projected 3D curve. `GeomProjLib::ProjectOnPlane`
    // (`GeomProjLib.cxx:313-366`) returns the projected curve in the same
    // representation as the source (Line / Circle / Ellipse / Bezier / BSpline
    // arms, `cxx:327-354`) and re-wraps a trimmed source in a `Geom_TrimmedCurve`
    // (`cxx:358-363`); `BRep_Tool::CurveOnPlane` then runs
    // `ProjLib_ProjectedCurve` (`BRep_Tool.cxx:399-410`), which preserves the
    // curve kind. So the sample count of this wrapper
    // (`Geom2dAdaptor_Curve::NbSamples`, `Geom2dAdaptor_Curve.cxx:1351-1394`) is
    // the source curve's count, not the trait default 20. A bare default makes
    // a trimmed/Bezier/BSpline edge projected onto a plane face under-sample.
    //
    // `gp_circ2d` is NOT delegated: the projection of a circle is a circle only
    // when the plane is parallel to the circle's plane, otherwise an ellipse
    // (`ProjLib_Plane::Project(gp_Circ)`, `ProjLib_Plane.cxx:110-123`), and that
    // detection is not ported here. Reporting `None` keeps the circle arm of
    // `Geom2dInt_Geom2dCurveTool::NbSamples` off instead of guessing.
    fn is_line(&self) -> bool {
        self.curve.is_line()
    }
    fn bezier_nb_poles(&self) -> Option<usize> {
        self.curve.bezier_poles().map(|poles| poles.len())
    }
    fn bspline_nb_knots(&self) -> Option<usize> {
        self.curve.bspline_knots().map(|knots| knots.len())
    }
    fn bspline_degree(&self) -> Option<usize> {
        self.curve.nurbs_degree()
    }
    fn clone_dyn(&self) -> Box<dyn Curve2d> {
        Box::new(Self {
            curve: Arc::from(self.curve.clone_dyn()),
            pln: self.pln.clone(),
        })
    }
}

/// `BRep_Tool::CurveOnPlane` (`BRep_Tool.cxx:379-449`):
/// `GeomProjLib::ProjectOnPlane` (KeepParam=true) then `ProjLib_ProjectedCurve`
/// on the plane. Unwraps a 2d TrimmedCurve to its basis (`cxx:443-447`).
fn curve_on_plane(curve: &dyn Curve, surf: &dyn Surface) -> Option<Arc<dyn Curve2d>> {
    let pln = pln_from_surface(surf)?;
    project_curve_on_plane(curve, &pln)
}

/// `GeomProjLib::ProjectOnPlane` (KeepParam=true) then `ProjLib_ProjectedCurve`
/// on the plane given by its own `gp_Pln`: the normal projection of `curve`
/// expressed in that plane's `(u, v)`. Split out of [`curve_on_plane`] so the
/// caller can pass a plane that is not the surface's own frame - the
/// `ShapeConstruct_ProjectCurveOnSurface::projectAnalytic` arm projects onto the
/// *basis* plane behind a trimmed / offset surface wrapper
/// (`ShapeConstruct_ProjectCurveOnSurface.cxx:848-897`).
pub(crate) fn project_curve_on_plane(curve: &dyn Curve, pln: &GpPln) -> Option<Arc<dyn Curve2d>> {
    // `GeomProjLib.cxx:339-343`: a trimmed 3D curve stays KeepParam on its
    // own `[First, Last]` (our STEP trim is remapped to `[0, 1]`).
    if curve.is_geom_trimmed()
        || curve.bspline_poles().is_some()
        || curve.bezier_poles().is_some()
    {
        // Trimmed / BSpline / Bezier: KeepParam UV of C3D(t). Pole-copy 2d
        // BSpline is non-rational (`Geom2dBSplineCurve`) and same-t of a
        // rational 3D BSpline was 1.322 on Shape-2 (cxx KeepParam + weights).
        return Some(Arc::new(PlaneKeepParam2d {
            curve: Arc::from(curve.clone_dyn()),
            pln: pln.clone(),
        }));
    }
    if let Some(c) = curve.gp_circ() {
        // `ProjLib_Plane::Project(gp_Circ)` (`cxx:110-123`).
        let p2d = projlib::eval_pln_pnt2d(&pln, &c.location());
        let (xx, xy) = projlib::eval_pln_dir2d(&pln, c.position().x_direction());
        let (yx, yy) = projlib::eval_pln_dir2d(&pln, c.position().y_direction());
        let vx = GpDir2d::new(xx, xy).ok()?;
        let vy = GpDir2d::new(yx, yy).ok()?;
        let ax = GpAx22d::new(p2d, vx, vy).ok()?;
        return Some(Arc::new(Geom2dCircle::new(GpCirc2d::new(ax, c.radius()))));
    }
    if let Some(e) = curve.gp_ellipse() {
        // `ProjLib_Plane::Project(gp_Elips)` (`ProjLib_Plane.cxx:126-139`):
        // project the ellipse axes and keep the radii, so the 2D parameter
        // `u` stays the 3D ellipse parameter.
        let p2d = projlib::eval_pln_pnt2d(&pln, &e.location());
        let (xx, xy) = projlib::eval_pln_dir2d(&pln, e.position().x_direction());
        let (yx, yy) = projlib::eval_pln_dir2d(&pln, e.position().y_direction());
        let vx = GpDir2d::new(xx, xy).ok()?;
        let vy = GpDir2d::new(yx, yy).ok()?;
        let ax = GpAx22d::new(p2d, vx, vy).ok()?;
        return Some(Arc::new(Geom2dEllipse::new(GpElips2d::new(
            ax,
            e.major_radius(),
            e.minor_radius(),
        ))));
    }
    if curve.is_line() {
        // `ProjLib_Plane::Project(gp_Lin)` (`cxx:101-106`).
        let o = curve.d0(0.0);
        let dir3 = curve.d1(0.0).1;
        let d = GpDir::from_vec(&dir3).ok()?;
        let loc = projlib::eval_pln_pnt2d(&pln, &o);
        let (dx, dy) = projlib::eval_pln_dir2d(&pln, &d);
        let dir2 = GpDir2d::new(dx, dy).ok()?;
        return Some(Arc::new(Geom2dLine::new(GpAx2d::new(loc, dir2))));
    }
    None
}

fn pln_from_surface(surf: &dyn Surface) -> Option<GpPln> {
    let o = surf.d0(0.0, 0.0);
    let (_, du, dv) = surf.d1(0.0, 0.0);
    let z = GpDir::from_vec(&du.crossed(&dv)).ok()?;
    if let Ok(x) = GpDir::from_vec(&du) {
        if let Ok(ax) = GpAx3::new(o, z, &x) {
            return Some(GpPln::new(ax));
        }
    }
    Some(GpPln::new(GpAx3::from_ax1(&GpAx1::new(o, z))))
}

/// Analytic kind of a pcurve, deduced from geometric invariants of the
/// `Curve2d` trait (no downcasting available on `Arc<dyn Curve2d>`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CurveKind {
    Line,
    Circle,
    BSpline,
    Other,
}

/// Classify a 2D curve: periodic 2π ⇒ Circle, unbounded range ⇒ Line,
/// otherwise a finite-range B-spline.
pub fn pc_curve_kind(curve: &dyn Curve2d) -> CurveKind {
    if curve.is_periodic() && (curve.period() - 2.0 * std::f64::consts::PI).abs() < 1e-9 {
        CurveKind::Circle
    } else if !curve.first_parameter().is_finite() || !curve.last_parameter().is_finite() {
        CurveKind::Line
    } else {
        CurveKind::BSpline
    }
}

/// Construct the pcurve of `edge` on `face`: the edge's 3D curve projected
/// into the face surface's `(u, v)` parameter domain.
///
/// The returned pcurve is parameterized over the edge's range `[first, last]`
/// so `d0(first)`/`d0(last)` land on the projected endpoints.
pub fn make_pcurve_on_face(edge: &Edge, face: &Face) -> Result<Arc<dyn Curve2d>, String> {
    let Some(curve) = GeometryRegistry::global().edge_curve(&edge.0) else {
        return Err("make_pcurve_on_face: edge has no 3D curve".into());
    };
    let (a, b) = GeometryRegistry::global().edge_parameters(&edge.0);
    if !a.is_finite() || !b.is_finite() || b - a < 1e-15 {
        return Err("make_pcurve_on_face: edge range is empty or unbounded".into());
    }

    let Some(surf) = GeometryRegistry::global().face_surface(&face.0) else {
        return Err("make_pcurve_on_face: face has no surface".into());
    };
    let planar = is_planar(surf.as_ref(), 8, 8, 1e-6);
    if planar {
        if let Some(pc) = curve_on_plane(curve.as_ref(), surf.as_ref()) {
            return Ok(pc);
        }
    }

    // Analytic cases on plane/cylinder faces.
    if let Ok(proj) = face_projector(face) {
        if is_line(curve.as_ref()) {
            return line_pcurve_on_face(curve.as_ref(), &proj, a, b);
        }
        if is_circle(curve.as_ref()) {
            if let Some(c) = circle_pcurve_on_face(curve.as_ref(), &proj, a, b) {
                return Ok(c);
            }
            // Phase mismatch or non-conformal projection: fall through to the
            // sampling path, which still reproduces the endpoint UVs.
        }
    }

    // General: sample the projected points and fit a polyline B-spline.
    let pts = edge_pcurve_on_face(edge, face, 32);
    if pts.len() < 2 {
        return Err("make_pcurve_on_face: too few projected samples".into());
    }
    let bs = bspline_from_samples(&pts, a, b, 1)?;
    Ok(Arc::new(bs))
}

/// `Adaptor3d_Curve::GetType() == GeomAbs_Line`: a `Geom_Line` (or a trimmed
/// one). This is the *concrete class* test of `GeomAdaptor_Curve::load`
/// (`GeomAdaptor_Curve.cxx:252-311`), which unwraps a `Geom_TrimmedCurve` to
/// its basis and then tests `gp_line()`. It replaces the earlier
/// "unbounded parameter range" heuristic, which also matched an offset curve.
fn is_line(curve: &dyn Curve) -> bool {
    curve.gp_line().is_some()
}

/// `Adaptor3d_Curve::GetType() == GeomAbs_Circle`: a `Geom_Circle` (or a
/// trimmed one). Deliberately NOT a geometric test - a `Geom_BSplineCurve`
/// whose image is a circle stays `GeomAbs_BSplineCurve` in OCCT
/// (`GeomAdaptor_Curve.cxx:252-311`), and `ProjLib_ProjectedCurve::Project`
/// only has analytic overloads for `GeomAbs_Line/Circle/Ellipse/Hyperbola/
/// Parabola` (`ProjLib_ProjectedCurve.cxx:247-259`): `GeomAbs_BSplineCurve`
/// and the rest break out (`:262-266`) into the general approximation. The
/// earlier "periodic with period 2π" heuristic classified a B-spline circle
/// as an analytic one, which produced an open `Geom2d_Line` spanning a full
/// period for a closed edge (T-69).
fn is_circle(curve: &dyn Curve) -> bool {
    curve.gp_circ().is_some()
}

/// UV projection for an analytic face (plane or cylinder). Recovered from the
/// face surface's own parameterization so the results agree with
/// `brep_surface::edge_pcurve_on_face`.
enum Projector {
    Plane { o: GpPnt, x: GpVec, y: GpVec },
    Cylinder { o: GpPnt, x: GpVec, y: GpVec, z: GpVec },
}

impl Projector {
    /// Project a 3D point onto the face's `(u, v)` parameter domain.
    fn project(&self, p: &GpPnt) -> GpPnt2d {
        match self {
            Projector::Plane { o, x, y } => {
                let d = p.coord.subtracted(&o.coord);
                // Solve [X·X X·Y; X·Y Y·Y] [u; v] = [d·X; d·Y] for an affine
                // (possibly non-orthonormal) plane parameterization.
                let xx = x.xyz().dot(x.xyz());
                let yy = y.xyz().dot(y.xyz());
                let xy = x.xyz().dot(y.xyz());
                let dx = d.dot(x.xyz());
                let dy = d.dot(y.xyz());
                let det = xx * yy - xy * xy;
                if det.abs() < 1e-24 {
                    GpPnt2d::new(dx, dy)
                } else {
                    GpPnt2d::new((dx * yy - dy * xy) / det, (dy * xx - dx * xy) / det)
                }
            }
            Projector::Cylinder { o, x, y, z } => {
                // Surface: P(u,v) = O + X·r·cos u + Y·r·sin u + Z·v, so
                // u = atan2(d·Y, d·X), v = d·Z with d = p − O.
                let d = p.coord.subtracted(&o.coord);
                let u = d.dot(y.xyz()).atan2(d.dot(x.xyz()));
                let v = d.dot(z.xyz());
                GpPnt2d::new(u, v)
            }
        }
    }
}

/// Recover the plane's location and in-plane axes from its surface `d0`.
fn plane_axes(surf: &dyn Surface) -> Option<(GpPnt, GpVec, GpVec)> {
    let o = surf.d0(0.0, 0.0);
    let x = GpVec::from_pnts(&o, &surf.d0(1.0, 0.0));
    let y = GpVec::from_pnts(&o, &surf.d0(0.0, 1.0));
    if x.xyz().square_modulus() < 1e-24 || y.xyz().square_modulus() < 1e-24 {
        return None;
    }
    Some((o, x, y))
}

/// Recover a cylinder's axis point (at v = 0), unit axes and radius from its
/// surface `d0`. A U-periodic surface is assumed to be a cylinder.
fn cylinder_axes(surf: &dyn Surface) -> Option<(GpPnt, GpVec, GpVec, GpVec)> {
    let a = surf.d0(0.0, 0.0);
    let b = surf.d0(std::f64::consts::PI, 0.0);
    let c = surf.d0(std::f64::consts::FRAC_PI_2, 0.0);
    let o = midpoint(&a, &b);
    let r = 0.5 * a.distance(&b);
    if r < 1e-12 {
        return None;
    }
    let x = GpVec::new((a.x() - o.x()) / r, (a.y() - o.y()) / r, (a.z() - o.z()) / r);
    let y = GpVec::new((c.x() - o.x()) / r, (c.y() - o.y()) / r, (c.z() - o.z()) / r);
    let z = GpVec::from_pnts(&a, &surf.d0(0.0, 1.0));
    if z.xyz().square_modulus() < 1e-24 {
        return None;
    }
    Some((o, x, y, z))
}

/// Build a UV projector for a face whose surface is a plane or a cylinder.
fn face_projector(face: &Face) -> Result<Projector, String> {
    let Some(surf) = GeometryRegistry::global().face_surface(&face.0) else {
        return Err("face has no surface".into());
    };
    if is_planar(surf.as_ref(), 8, 8, 1e-6) {
        let (o, x, y) = plane_axes(surf.as_ref()).ok_or("plane axes not recoverable")?;
        Ok(Projector::Plane { o, x, y })
    } else if surf.is_u_periodic() {
        let (o, x, y, z) = cylinder_axes(surf.as_ref()).ok_or("cylinder axes not recoverable")?;
        Ok(Projector::Cylinder { o, x, y, z })
    } else {
        Err("face surface is neither plane nor cylinder".into())
    }
}

/// Project a line edge onto an analytic face. Returns a `Geom2dLine` when the
/// projection is unit-speed (the common arc-length parameterized case); a
/// degree-1 B-spline otherwise preserves the exact endpoint UVs over [a, b].
fn line_pcurve_on_face(
    curve: &dyn Curve,
    proj: &Projector,
    a: f64,
    b: f64,
) -> Result<Arc<dyn Curve2d>, String> {
    let q0 = proj.project(&curve.d0(a));
    let q1 = proj.project(&curve.d0(b));
    let d = GpVec2d::new(q1.x() - q0.x(), q1.y() - q0.y());
    let mag = d.magnitude();
    if mag < 1e-15 {
        // Degenerate UV extent — a point pcurve.
        let bs = bspline_from_samples(&[q0, q1], a, b, 1)?;
        return Ok(Arc::new(bs));
    }
    let dir = GpDir2d::from_vec2d(&d).map_err(|e| e.to_string())?;
    let len = b - a;
    if (mag - len).abs() <= 1e-6 * len.max(1.0) {
        // Unit-speed: Geom2dLine positioned so d0(a) = q0 and d0(b) = q1.
        let loc = GpPnt2d::new(q0.x() - a * dir.x, q0.y() - a * dir.y);
        Ok(Arc::new(Geom2dLine::new(GpAx2d::new(loc, dir))))
    } else {
        // Non-unit-speed: exact polyline over [a, b].
        let bs = bspline_from_samples(&[q0, q1], a, b, 1)?;
        Ok(Arc::new(bs))
    }
}

/// Project a circle edge onto an analytic face. Returns a `Geom2dCircle` when
/// the projection is an (unrotated) circle — the case of a circle in a plane
/// parallel to the face; otherwise `None` lets the caller sample instead.
fn circle_pcurve_on_face(
    curve: &dyn Curve,
    proj: &Projector,
    a: f64,
    b: f64,
) -> Option<Arc<dyn Curve2d>> {
    // 3D circle centre from two opposite samples (valid for any trim range).
    let center3 = midpoint(&curve.d0(0.0), &curve.d0(std::f64::consts::PI));
    let center_uv = proj.project(&center3);
    let q0 = proj.project(&curve.d0(a));
    let q1 = proj.project(&curve.d0(b));
    let r_uv = center_uv.distance(&q0);
    if r_uv < 1e-12 {
        return None;
    }
    // Phase of the projected start relative to the 2D circle parameter.
    let phi = (q0.y() - center_uv.y()).atan2(q0.x() - center_uv.x());
    if normalize_angle(phi - a).abs() > 1e-6 {
        return None; // rotated projection → sampling fallback
    }
    let ax = GpAx22d::from_xdir(center_uv, GpDir2d::default());
    let circ = Geom2dCircle::new(GpCirc2d::new(ax, r_uv));
    let tol = 1e-6 * r_uv.max(1.0);
    if circ.d0(a).distance(&q0) <= tol && circ.d0(b).distance(&q1) <= tol {
        Some(Arc::new(circ))
    } else {
        None
    }
}

/// Fit a clamped B-spline through `pts` with the exact parameter range [a, b].
/// Degree 1 is a polyline through the sample points (endpoints interpolated).
fn bspline_from_samples(
    pts: &[GpPnt2d],
    a: f64,
    b: f64,
    degree: usize,
) -> Result<Geom2dBSplineCurve, String> {
    if pts.len() < 2 {
        return Err("bspline_from_samples: need at least 2 points".into());
    }
    if degree >= pts.len() {
        return Err("bspline_from_samples: degree must be below pole count".into());
    }
    let mut knots = occt_core::bspl::knots::build_uniform_knots(pts.len(), degree);
    for k in knots.iter_mut() {
        *k = a + (b - a) * *k;
    }
    let xs: Vec<f64> = pts.iter().map(|p| p.x()).collect();
    let ys: Vec<f64> = pts.iter().map(|p| p.y()).collect();
    Geom2dBSplineCurve::new(xs, ys, knots, degree).map_err(|e| e.to_string())
}

fn midpoint(a: &GpPnt, b: &GpPnt) -> GpPnt {
    GpPnt::new(0.5 * (a.x() + b.x()), 0.5 * (a.y() + b.y()), 0.5 * (a.z() + b.z()))
}

fn normalize_angle(x: f64) -> f64 {
    let two_pi = 2.0 * std::f64::consts::PI;
    let mut r = x % two_pi;
    if r > std::f64::consts::PI {
        r -= two_pi;
    } else if r < -std::f64::consts::PI {
        r += two_pi;
    }
    r
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep_extrema::test_box::unit_box;
    use crate::brep_surface::edge_pcurve_on_face;
    use crate::builder::TopoBuilder;
    use occt_core::gp::{GpAx2, GpAx3, GpCylinder, GpPln, GpPnt};
    use occt_geom::{GeomBSplineCurve, GeomCylinder};

    #[test]
    fn box_bottom_line_edge_projects_to_line_pcurve() {
        let b = unit_box();
        let edge = &b.edges[0]; // (0,0,0) → (1,0,0)
        let face = &b.faces[0]; // bottom (z = 0)
        let pc = make_pcurve_on_face(edge, face).expect("pcurve");
        assert_eq!(pc_curve_kind(pc.as_ref()), CurveKind::Line);
        let (a, z) = GeometryRegistry::global().edge_parameters(&edge.0);
        let samples = edge_pcurve_on_face(edge, face, 8);
        let p0 = pc.d0(a);
        let p1 = pc.d0(z);
        assert!((p0.x() - samples[0].x()).abs() < 1e-6, "start u: {} vs {}", p0.x(), samples[0].x());
        assert!((p0.y() - samples[0].y()).abs() < 1e-6);
        assert!((p1.x() - samples[7].x()).abs() < 1e-6);
        assert!((p1.y() - samples[7].y()).abs() < 1e-6);
    }

    #[test]
    fn box_top_edge_hits_face_boundary() {
        let b = unit_box();
        let edge = &b.edges[4]; // (0,0,1) → (1,0,1)
        let face = &b.faces[1]; // top (z = 1)
        let pc = make_pcurve_on_face(edge, face).expect("pcurve");
        assert_eq!(pc_curve_kind(pc.as_ref()), CurveKind::Line);
        // The top face's UV domain is [0,1]²; both projected corners sit on
        // the v = 0 boundary.
        let p0 = pc.d0(0.0);
        let p1 = pc.d0(1.0);
        assert!((p0.x() - 0.0).abs() < 1e-6 && (p0.y() - 0.0).abs() < 1e-6, "start {:?}", p0);
        assert!((p1.x() - 1.0).abs() < 1e-6 && (p1.y() - 0.0).abs() < 1e-6, "end {:?}", p1);
    }

    #[test]
    fn circle_edge_on_plane_face_is_circle_pcurve() {
        let b = TopoBuilder::new();
        let face = b.make_face_plane(&GpPln::new(GpAx3::standard()));
        let ax2 = GpAx2::standard();
        let edge = b.make_edge_circle(&ax2, 1.0, 0.0, 2.0 * std::f64::consts::PI);
        let pc = make_pcurve_on_face(&edge, &face).expect("pcurve");
        assert_eq!(pc_curve_kind(pc.as_ref()), CurveKind::Circle);
        // Projection of (cos t, sin t, 0) on the XY plane is (cos t, sin t).
        let p0 = pc.d0(0.0);
        let pm = pc.d0(std::f64::consts::FRAC_PI_2);
        assert!((p0.x() - 1.0).abs() < 1e-6 && (p0.y() - 0.0).abs() < 1e-6, "start {:?}", p0);
        assert!((pm.x() - 0.0).abs() < 1e-6 && (pm.y() - 1.0).abs() < 1e-6, "mid {:?}", pm);
    }

    #[test]
    fn cylinder_generatrix_projects_to_line_pcurve() {
        let b = TopoBuilder::new();
        let cyl = GpCylinder::new(GpAx3::standard(), 1.0).unwrap();
        let face = b.make_face(Arc::new(GeomCylinder::new(cyl)), &[]);
        let edge = b.make_edge_segment(&GpPnt::new(1.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 1.0));
        let pc = make_pcurve_on_face(&edge, &face).expect("pcurve");
        assert_eq!(pc_curve_kind(pc.as_ref()), CurveKind::Line);
        // Generatrix at angle 0: u = atan2(0, 1) = 0, v = axial coordinate.
        let p0 = pc.d0(0.0);
        let p1 = pc.d0(1.0);
        assert!((p0.x() - 0.0).abs() < 1e-6 && (p0.y() - 0.0).abs() < 1e-6, "start {:?}", p0);
        assert!((p1.x() - 0.0).abs() < 1e-6 && (p1.y() - 1.0).abs() < 1e-6, "end {:?}", p1);
    }

    #[test]
    fn general_bspline_edge_samples_hit_face_boundary() {
        let b = TopoBuilder::new();
        let face = b.make_face_plane(&GpPln::new(GpAx3::standard()));
        let poles = vec![
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(0.5, 0.4, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
        ];
        let curve = Arc::new(
            GeomBSplineCurve::new(poles, vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 2).unwrap(),
        );
        let edge = b.make_edge(curve, 0.0, 1.0);
        let pc = make_pcurve_on_face(&edge, &face).expect("pcurve");
        assert_eq!(pc_curve_kind(pc.as_ref()), CurveKind::BSpline);
        // Clamped endpoints interpolate p0 / p2, which project onto the face
        // boundary UV (0,0) and (1,0).
        let p0 = pc.d0(0.0);
        let p1 = pc.d0(1.0);
        assert!((p0.x() - 0.0).abs() < 1e-6 && (p0.y() - 0.0).abs() < 1e-6, "start {:?}", p0);
        assert!((p1.x() - 1.0).abs() < 1e-6 && (p1.y() - 0.0).abs() < 1e-6, "end {:?}", p1);
    }

    #[test]
    fn registry_pcurve_roundtrip_and_default_empty() {
        let reg = GeometryRegistry::global();
        let b = TopoBuilder::new();
        let edge = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 0.0));
        let face = b.make_face_plane(&GpPln::new(GpAx3::standard()));
        let face_key = GeometryRegistry::shape_key(&face.0);
        // Default: no pcurve attached.
        assert!(reg.edge_pcurve(&edge.0, face_key).is_none());
        let pc: Arc<dyn Curve2d> = Arc::new(Geom2dLine::new(GpAx2d::new(
            GpPnt2d::zero(),
            GpDir2d::default(),
        )));
        reg.set_edge_pcurve(&edge.0, face_key, pc.clone());
        let got = reg.edge_pcurve(&edge.0, face_key).expect("pcurve present");
        assert!((got.d0(3.5).x() - 3.5).abs() < 1e-12);
        // Clearing the edge drops the attached pcurves too.
        reg.clear_shape(&edge.0);
        assert!(reg.edge_pcurve(&edge.0, face_key).is_none());
    }
}
