//! `GeomProjLib::Curve2d` — the 2D projection of a 3D curve onto a surface.
//!
//! Source: `ModelingData/TKGeomBase/GeomProjLib/GeomProjLib.cxx:54-203`
//! (OCCT 8.0.0).
//!
//! The body builds a `GeomAdaptor_Curve AC(C, First, Last)` and a
//! `GeomAdaptor_Surface AS(S, UDeb, UFin, VDeb, VFin)`, runs
//! `ProjLib_ProjectedCurve Proj(HS, HC, Tolerance)` and maps `Proj.GetType()` to
//! the matching `Geom2d_*` curve; a type outside the seven analytic cases
//! returns a null handle (`cxx:108-111`). The tolerance is raised to
//! `Precision::PConfusion()` first (`cxx:64`) and overwritten with
//! `Proj.GetTolerance()` on the way out (`cxx:114`, `:130`).
//!
//! The two adaptors only carry the parameter range into the projector. In the
//! port the projected curve reads the range from the curve and the UV bounds
//! from the surface itself ([`crate::projlib::ProjectedCurve::perform`]), and
//! every caller of this function passes `C->FirstParameter()` /
//! `C->LastParameter()` with `S->Bounds(...)` — the adaptors' own defaults — so
//! [`curve_2d_range`] keeps the ranges as the trim clamp of `cxx:119-128` and
//! lets the projector see the same range through the curve.
//!
//! The trailing `Geom2d_TrimmedCurve` (`cxx:119-128`) is re-applied on a
//! `Geom_TrimmedCurve` input only, and only when the projected curve is not
//! periodic; the port asks [`Curve::is_geom_trimmed`] for the `IsKind` test.

use std::sync::Arc;

use occt_core::gp::{GpDir, GpPln, GpPnt, GpTrsf, GpVec};
use occt_core::kernel::geomabs::CurveType;
use occt_core::precision;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;
use occt_geom2d::{
    Curve2d, Geom2dCircle, Geom2dEllipse, Geom2dHyperbola, Geom2dLine, Geom2dParabola,
};

use crate::curve::Curve;
use crate::projlib::ProjectedCurve;
use crate::surface::Surface;

/// `GeomProjLib::Curve2d(C, S)` (`cxx:185-192`): the curve's own range,
/// `Tol = Precision::PConfusion()`.
pub fn curve_2d(c: &Arc<dyn Curve>, s: &Arc<dyn Surface>) -> Option<Arc<dyn Curve2d>> {
    let mut tolerance = precision::CONFUSION;
    curve_2d_with_tolerance(c, s, &mut tolerance)
}

/// `GeomProjLib::Curve2d(C, S, Tolerance)` (`cxx:196-203`).
pub fn curve_2d_with_tolerance(
    c: &Arc<dyn Curve>,
    s: &Arc<dyn Surface>,
    tolerance: &mut f64,
) -> Option<Arc<dyn Curve2d>> {
    curve_2d_range(c, c.first_parameter(), c.last_parameter(), s, tolerance)
}

/// `GeomProjLib::Curve2d(C, First, Last, S, UDeb, UFin, VDeb, VFin,
/// Tolerance)` (`cxx:54-134`).
///
/// `udeb` / `ufin` / `vdeb` / `vfin` are the `GeomAdaptor_Surface` bounds of
/// `cxx:68`; the port's projector takes them from the surface, which is the
/// same box for every caller of this overload, so they are not passed here.
pub fn curve_2d_range(
    c: &Arc<dyn Curve>,
    first: f64,
    last: f64,
    s: &Arc<dyn Surface>,
    tolerance: &mut f64,
) -> Option<Arc<dyn Curve2d>> {
    // `Tolerance = std::max(Precision::PConfusion(), Tolerance);` (`cxx:64`).
    *tolerance = tolerance.max(precision::CONFUSION);
    // `ProjLib_ProjectedCurve Proj(HS, HC, Tolerance);` (`cxx:74`).
    let proj = ProjectedCurve::with_surface_curve_tol(Arc::clone(s), Arc::clone(c), *tolerance);

    // `switch (Proj.GetType())` (`cxx:78-111`).
    let g2d_c: Option<Arc<dyn Curve2d>> = match proj.get_type() {
        CurveType::Line => proj
            .line()
            .ok()
            .map(|l| Arc::new(Geom2dLine::new(*l.position())) as Arc<dyn Curve2d>),
        CurveType::Circle => proj
            .circle()
            .ok()
            .map(|l| Arc::new(Geom2dCircle::new(l)) as Arc<dyn Curve2d>),
        CurveType::Ellipse => proj
            .ellipse()
            .ok()
            .map(|l| Arc::new(Geom2dEllipse::new(l)) as Arc<dyn Curve2d>),
        CurveType::Parabola => proj
            .parabola()
            .ok()
            .map(|l| Arc::new(Geom2dParabola::new(l)) as Arc<dyn Curve2d>),
        CurveType::Hyperbola => proj
            .hyperbola()
            .ok()
            .map(|l| Arc::new(Geom2dHyperbola::new(l)) as Arc<dyn Curve2d>),
        CurveType::BezierCurve => proj.bezier().map(|b| b as Arc<dyn Curve2d>),
        CurveType::BSplineCurve => proj.bspline().map(|b| b as Arc<dyn Curve2d>),
        _ => return None,
    };

    // `if (G2dC.IsNull()) { Tolerance = Proj.GetTolerance(); return G2dC; }`
    // (`cxx:113-117`).
    let Some(g2d_c) = g2d_c else {
        *tolerance = proj.get_tolerance();
        return None;
    };

    // `if (C->IsKind(STANDARD_TYPE(Geom_TrimmedCurve))) { ... }` (`cxx:119-128`).
    let g2d_c: Arc<dyn Curve2d> = if c.is_geom_trimmed() {
        let mut u1 = first;
        let mut u2 = last;
        if !g2d_c.is_periodic() {
            u1 = u1.max(g2d_c.first_parameter());
            u2 = u2.min(g2d_c.last_parameter());
        }
        Arc::new(Geom2dTrimmedCurve::new(g2d_c, u1, u2))
    } else {
        g2d_c
    };

    // `Tolerance = Proj.GetTolerance();` (`cxx:130`).
    *tolerance = proj.get_tolerance();
    Some(g2d_c)
}

/// `GeomProjLib::ProjectOnPlane(Curve, Plane, Dir, KeepParametrization)`
/// (`GeomProjLib.cxx:293-345`).
///
/// The reference projects every point of `Curve` along `Dir` onto the plane
/// (`ProjLib_ProjectOnPlane.cxx` `OnPlane_Value` / `OnPlane_D1` / `OnPlane_D2`):
///
/// ```text
/// Proj(U) = P(U) + alpha * D,   alpha = ((O - P(U)) . Z) / (D . Z)
/// ```
///
/// with `O` the plane location and `Z` its direction, and the derivatives
/// `d/du Proj = d/du P - ((d/du P . Z) / (D . Z)) * D` (linear, so the same for
/// every order). `GeomProjLib.cxx:307-335` then re-wraps the result in the
/// analytic type `Proj.GetType()` reports.
///
/// UNPORTED: `ProjLib_ProjectOnPlane` (`ProjLib_ProjectOnPlane.cxx:1-1567`) is
/// not ported, so there is no analytic dispatch (Line / Circle / Ellipse /
/// Parabola / Hyperbola / Bezier / BSpline, `GeomProjLib.cxx:307-335`) and no
/// approximation (`Approx_FitAndDivide`, `ProjLib_ProjectOnPlane.cxx`
/// `Project()`'s `!myKeepParam` arm) here. The port answers a generic curve in
/// world coordinates on the plane, at the source parameter of `Curve`, and the
/// points and derivatives are the reference's `OnPlane_Value` / `OnPlane_D1` /
/// `OnPlane_D2` exactly. The trailing `Geom_TrimmedCurve` re-wrap
/// (`GeomProjLib.cxx:339-343`) is not needed because the source curve already
/// carries the trim interval in that arm.
pub fn project_on_plane(c: &Arc<dyn Curve>, pln: &GpPln, dir: &GpDir) -> Arc<dyn Curve> {
    Arc::new(PlaneProjected3d {
        curve: Arc::clone(c),
        pln: pln.clone(),
        dir: *dir,
    })
}

/// The `OnPlane_*` family of `ProjLib_ProjectOnPlane.cxx`.
struct PlaneProjected3d {
    curve: Arc<dyn Curve>,
    pln: GpPln,
    dir: GpDir,
}

impl PlaneProjected3d {
    /// `(D . Z)` (`ProjLib_ProjectOnPlane.cxx`, `OnPlane_Value`). A zero
    /// denominator (projection direction in the plane) raises in the reference;
    /// the port keeps a finite value so the curve stays evaluable.
    fn denom(&self) -> f64 {
        let z = self.pln.position().direction();
        self.dir.dot(&z)
    }

    fn plane_z(&self) -> GpVec {
        let z = self.pln.position().direction();
        GpVec::new(z.x(), z.y(), z.z())
    }

    fn dir_vec(&self) -> GpVec {
        GpVec::new(self.dir.x(), self.dir.y(), self.dir.z())
    }

    fn project_pnt(&self, p: &GpPnt) -> GpPnt {
        let z = self.plane_z();
        let po = GpVec::from_pnts(p, &self.pln.location());
        let alpha = po.dot(&z) / self.denom();
        let d = self.dir_vec();
        GpPnt::new(
            p.x() + alpha * d.x(),
            p.y() + alpha * d.y(),
            p.z() + alpha * d.z(),
        )
    }

    fn project_vec(&self, v: &GpVec) -> GpVec {
        let z = self.plane_z();
        let k = v.dot(&z) / self.denom();
        let d = self.dir_vec();
        GpVec::new(v.x() - k * d.x(), v.y() - k * d.y(), v.z() - k * d.z())
    }
}

impl Curve for PlaneProjected3d {
    fn d0(&self, u: f64) -> GpPnt {
        self.project_pnt(&self.curve.d0(u))
    }

    fn d1(&self, u: f64) -> (GpPnt, GpVec) {
        let (p, v) = self.curve.d1(u);
        (self.project_pnt(&p), self.project_vec(&v))
    }

    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        let (p, v1, v2) = self.curve.d2(u);
        (
            self.project_pnt(&p),
            self.project_vec(&v1),
            self.project_vec(&v2),
        )
    }

    fn d3(&self, u: f64) -> (GpPnt, GpVec, GpVec, GpVec) {
        let (p, v1, v2, v3) = self.curve.d3(u);
        (
            self.project_pnt(&p),
            self.project_vec(&v1),
            self.project_vec(&v2),
            self.project_vec(&v3),
        )
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

    fn parameter_intervals(&self, continuity: u8) -> Vec<f64> {
        self.curve.parameter_intervals(continuity)
    }

    /// The reference has no `Transform` on `ProjLib_ProjectOnPlane`; the port
    /// carries the transform into both the source curve and the plane, which
    /// keeps `Proj` of the moved points for a rigid `t`.
    fn transform(&mut self, t: &GpTrsf) {
        let mut c = self.curve.clone_dyn();
        c.transform(t);
        self.curve = Arc::from(c);
        self.pln.transform(t);
        let mut xyz = *self.dir.xyz();
        t.transforms_xyz(&mut xyz);
        if let Ok(d) = GpDir::from_xyz(&xyz) {
            self.dir = d;
        }
    }

    fn reverse(&mut self) {
        let mut c = self.curve.clone_dyn();
        c.reverse();
        self.curve = Arc::from(c);
    }

    fn clone_dyn(&self) -> Box<dyn Curve> {
        Box::new(PlaneProjected3d {
            curve: Arc::clone(&self.curve),
            pln: self.pln.clone(),
            dir: self.dir,
        })
    }
}
