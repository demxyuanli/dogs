//! `GeomProjLib::ProjectOnPlane` (KeepParam=true) and `ProjLib_ProjectedCurve`
//! on a plane -- the geometry half of `BRep_Tool::CurveOnPlane`
//! (`BRep_Tool.cxx:379-449`).
//!
//! Pure geometry: a 3D curve in, a `Geom2d_Curve` in the plane's `(u, v)` out.
//! `occt-topo`'s `pcurve` module re-exports it, and `occt-topology`'s
//! `brep_tool` uses it for the `BRep_Tool::CurveOnSurface` plane fallback
//! (`BRep_Tool.cxx:367-372`).

use std::sync::Arc;

use occt_core::gp::{
    GpAx1, GpAx2d, GpAx22d, GpAx3, GpCirc2d, GpDir, GpDir2d, GpElips2d, GpPln, GpPnt, GpPnt2d,
    GpVec, GpVec2d,
};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::{
    Geom2dBezierCurve, Geom2dBSplineCurve, Geom2dCircle, Geom2dEllipse, Geom2dHyperbola, Geom2dLine,
    Geom2dParabola,
};

use crate::projlib;
use crate::{Curve, Surface};
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
    /// `Geom2dAdaptor_Curve::Resolution` over `ProjLib_ProjectedCurve`'s result.
    ///
    /// `GeomProjLib::ProjectOnPlane` is called with `KeepParam = true`
    /// (`BRep_Tool.cxx:379-449`), so this wrapper shares the source curve's
    /// parameter space (`first_parameter`/`last_parameter` are delegated) and
    /// `ProjLib_ProjectedCurve` preserves the curve kind. Like the other kind
    /// queries on this wrapper the answer therefore comes from the source curve.
    /// The projection contracts lengths, so the faithful projected-curve
    /// resolution is at least this value; returning the source's keeps the step
    /// on the conservative side instead of falling back to the trait default
    /// (`Precision::Parametric`), which would under-step by orders of magnitude.
    fn resolution(&self, ruv: f64) -> f64 {
        self.curve.resolution(ruv)
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
pub fn curve_on_plane(curve: &dyn Curve, surf: &dyn Surface) -> Option<Arc<dyn Curve2d>> {
    let pln = pln_from_surface(surf)?;
    project_curve_on_plane(curve, &pln)
}

/// `ProjLib_ComputeApprox.cxx:1197-1252`: project the poles of a B-spline /
/// Bezier onto the plane, keeping the knot vector, multiplicities and weights.
///
/// The plane's `(u, v)` map is affine, so the pole-wise image is the same curve
/// as the projected curve expressed in `(u, v)` -- no approximation is needed,
/// which is exactly why OCCT short-circuits `ProjLib_ComputeApprox` for
/// `GeomAbs_Plane`.
fn projected_pole_curve(curve: &dyn Curve, pln: &GpPln) -> Option<Arc<dyn Curve2d>> {
    if let (Some(poles), Some(knots), Some(degree)) = (
        curve.bspline_poles(),
        curve.bspline_knots(),
        curve.nurbs_degree(),
    ) {
        let uv: Vec<GpPnt2d> = poles
            .iter()
            .map(|p| projlib::eval_pln_pnt2d(pln, p))
            .collect();
        let xs: Vec<f64> = uv.iter().map(|p| p.x()).collect();
        let ys: Vec<f64> = uv.iter().map(|p| p.y()).collect();
        let mut bs = match curve.bspline_weights() {
            Some(w) if w.len() == poles.len() => Geom2dBSplineCurve::rational(
                xs,
                ys,
                w.to_vec(),
                knots.to_vec(),
                degree,
            )
            .ok()?,
            _ => Geom2dBSplineCurve::new(xs, ys, knots.to_vec(), degree).ok()?,
        };
        // `Geom2d_BSplineCurve(Poles, ..., BS->IsPeriodic())` (`cxx:1219-1225`).
        if curve.is_periodic() {
            bs.set_periodic();
        }
        return Some(Arc::new(bs));
    }
    if let Some(poles) = curve.bezier_poles() {
        let uv: Vec<GpPnt2d> = poles
            .iter()
            .map(|p| projlib::eval_pln_pnt2d(pln, p))
            .collect();
        let bez = match curve.bezier_weights() {
            Some(w) if w.len() == poles.len() => {
                Geom2dBezierCurve::rational(uv, w.to_vec()).ok()?
            }
            _ => Geom2dBezierCurve::new(uv).ok()?,
        };
        return Some(Arc::new(bez));
    }
    None
}

/// `GeomProjLib::ProjectOnPlane` (KeepParam=true) then `ProjLib_ProjectedCurve`
/// on the plane given by its own `gp_Pln`: the normal projection of `curve`
/// expressed in that plane's `(u, v)`. Split out of [`curve_on_plane`] so the
/// caller can pass a plane that is not the surface's own frame - the
/// `ShapeConstruct_ProjectCurveOnSurface::projectAnalytic` arm projects onto the
/// *basis* plane behind a trimmed / offset surface wrapper
/// (`ShapeConstruct_ProjectCurveOnSurface.cxx:848-897`).
pub fn project_curve_on_plane(curve: &dyn Curve, pln: &GpPln) -> Option<Arc<dyn Curve2d>> {
    // `BRep_Tool::CurveOnPlane` wraps `Geom_TrimmedCurve(C,f,l)` then
    // `GeomProjLib::ProjectOnPlane`. OCCT's trim keeps basis parameters; unwrap
    // to the basis for the pole / analytic arms so a basis-parameter trim does
    // not force the KeepParam fallback.
    let basis = curve
        .untrimmed_basis()
        .map(|(c, _, _)| c)
        .unwrap_or_else(|| Arc::from(curve.clone_dyn()));
    let basis = basis.as_ref();

    // `ProjLib_ComputeApprox.cxx:1197-1252`: for a B-spline or Bezier 3D curve
    // on a plane the projection is not approximated -- every pole is projected
    // onto the plane's `(u, v)` and the knots / multiplicities / weights are
    // kept, so the stored pcurve is a concrete `Geom2d_BSplineCurve` /
    // `Geom2d_BezierCurve` in the source curve's own parameter space. Keeping
    // the concrete class matters for the consumers that branch on it, e.g.
    // `GeomBndLib_Curve2d` picks `GeomBndLib_BSplineCurve2d::Box` (pole hull,
    // `GeomBndLib_BSplineCurve2d.cxx:52-57`) instead of the sampling box.
    if let Some(pc) = projected_pole_curve(basis, pln) {
        return Some(pc);
    }
    // `GeomProjLib.cxx:339-343`: a trimmed 3D curve stays KeepParam on its
    // own `[First, Last]` when the basis is not a concrete pole curve.
    if curve.is_geom_trimmed()
        || basis.bspline_poles().is_some()
        || basis.bezier_poles().is_some()
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
    if let Some(p) = curve.gp_parabola() {
        // `ProjLib_Plane::Project(gp_Parab)` (`cxx:143-154`): project the
        // position, keep `Focal()`.
        let p2d = projlib::eval_pln_pnt2d(&pln, &p.location());
        let (xx, xy) = projlib::eval_pln_dir2d(&pln, p.position().x_direction());
        let (yx, yy) = projlib::eval_pln_dir2d(&pln, p.position().y_direction());
        let vx = GpDir2d::new(xx, xy).ok()?;
        let vy = GpDir2d::new(yx, yy).ok()?;
        let ax = GpAx22d::new(p2d, vx, vy).ok()?;
        return Some(Arc::new(Geom2dParabola::from_axes(ax, p.focal)));
    }
    if let Some(h) = curve.gp_hyperbola() {
        // `ProjLib_Plane::Project(gp_Hypr)` (`cxx:158-169`): project the
        // position, keep both radii.
        let p2d = projlib::eval_pln_pnt2d(&pln, &h.location());
        let (xx, xy) = projlib::eval_pln_dir2d(&pln, h.position().x_direction());
        let (yx, yy) = projlib::eval_pln_dir2d(&pln, h.position().y_direction());
        let vx = GpDir2d::new(xx, xy).ok()?;
        let vy = GpDir2d::new(yx, yy).ok()?;
        let ax = GpAx22d::new(p2d, vx, vy).ok()?;
        return Some(Arc::new(Geom2dHyperbola::from_axes(
            ax,
            h.major_radius,
            h.minor_radius,
        )));
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

