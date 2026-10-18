//! OCCT port: GeomAPI projection, intersection and extremum helpers.
//! Source: `GeomAPI_ProjectPointOnCurve.hxx`, `GeomAPI_ProjectPointOnSurf.hxx`,
//! `GeomAPI_IntCS.hxx`, `GeomAPI_ExtremaCurveCurve.hxx`,
//! `GeomAPI_ExtremaSurfaceSurface.hxx` (see [`ExtremaSurfaceSurface`]).
//!
//! **Ported (faithful)**: [`project_point_on_curve`] = `Extrema_ExtPC` via
//! [`crate::extrema_pc`] with `GeomAPI_ProjectPointOnCurve`'s own
//! `IsDone`/minimum rules (`GeomAPI_ProjectPointOnCurve.cxx:135-155`), and
//! [`project_point_on_surface`] = `Extrema_ExtPS` via [`crate::extrema_surf::ExtPs`]
//! with `GeomAPI_ProjectPointOnSurf`'s rules
//! (`GeomAPI_ProjectPointOnSurf.cxx:214-246`).
//!
//! **A16 求交 half — 已解决（2026-09-21 批 81）**: the three port-internal
//! samplers that used to stand in for the intersection APIs
//! (`curve_surface_intersections`, `curve_curve_intersections`,
//! `curve_curve_distance`) are **deleted** together with their private helpers
//! and unit tests. Their production callers in `occt-topo` now use the faithful
//! engines: `crate::intcurvesurface::perform_curve_surface`
//! (= `IntCurveSurface_HInter`, `IntTools_EdgeFace.cxx:426-445`) and
//! `crate::edge_edge::EdgeEdge`
//! (= `IntTools_EdgeEdge`, `IntTools_EdgeEdge.cxx:185-243`). `occt-geom` cannot
//! host those (it is the lower crate), so nothing replaces them here.

use std::sync::Arc;

use occt_core::gp::{
    GpAx2d, GpAx22d, GpDir2d, GpPln, GpPnt, GpPnt2d, GpVec, GpVec2d,
};
use occt_geom2d::bezier_curve::Geom2dBezierCurve;
use occt_geom2d::bspline_curve::Geom2dBSplineCurve;
use occt_geom2d::circle::Geom2dCircle;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::ellipse::Geom2dEllipse;
use occt_geom2d::hyperbola::Geom2dHyperbola;
use occt_geom2d::line::Geom2dLine;
use occt_geom2d::parabola::Geom2dParabola;

use crate::curve::Curve;
use crate::extrema_ss::ExtremaExtSS;
use crate::projlib;
use crate::surface::Surface;

/// Orthogonal projection of a point onto a curve (nearest solution).
pub struct PointOnCurve {
    pub parameter: f64,
    pub point: GpPnt,
    pub distance: f64,
}

/// Orthogonal projection of a point onto a surface (nearest solution).
pub struct PointOnSurface {
    pub u: f64,
    pub v: f64,
    pub point: GpPnt,
    pub distance: f64,
}

/// Nearest point on `c` to `p`.
///
/// Port of `GeomAPI_ProjectPointOnCurve::Perform`
/// (`GeomAPI_ProjectPointOnCurve.cxx:135-155`): the extrema engine is
/// `Extrema_ExtPC` over the curve's own range (`myExtPC.Initialize(myC,
/// myC.FirstParameter(), myC.LastParameter())`, `cxx:58`), `IsDone()` requires
/// `IsDone() && NbExt() > 0` (`cxx:139`) and the answer is the **minimum**
/// (`cxx:143-152`). OCCT performs **no** endpoint adjustment here — that is
/// `ShapeAnalysis_Curve::Project`'s job (`ShapeAnalysis_Curve.cxx:161-182`).
/// The faithful `Extrema_ExtPC` port is [`crate::extrema_pc::point_curve_extrema`]
/// (audit A7/T-43: `ExtPElC` analytic arms + `Extrema_GGExtPC`'s `default:` arm,
/// with `None` standing for OCCT's `IsDone() == false`).
///
/// The previous body was an invented coarse scan + golden-section refinement
/// over an expanding window (audit A16).
pub fn project_point_on_curve(c: &dyn Curve, p: &GpPnt, _tol: f64) -> Option<PointOnCurve> {
    let e = crate::extrema_pc::point_curve_extrema(c, p)?;
    Some(PointOnCurve { parameter: e.u1, point: e.p2, distance: e.distance })
}

/// Nearest (u, v) parameters of `p` on `s` (the minimum of the extrema).
///
/// Port of `GeomAPI_ProjectPointOnSurf::Perform`
/// (`GeomAPI_ProjectPointOnSurf.cxx:214-246`): `Extrema_ExtPS` over the
/// surface's own range with `Tolerance` used for both parameters, `IsDone()`
/// requiring `IsDone() && NbExt() > 0` (`cxx:83`) and the answer being the
/// smallest `SquareDistance` (`cxx:88-100`). The engine is
/// [`crate::extrema_surf::ExtPs`] (audit A1/T-67 step 2): elementary surfaces go
/// to the exact `Extrema_ExtPElS` arms, extrusion/revolution surfaces to
/// `Extrema_ExtPExtS` / `Extrema_ExtPRevS`, everything else to the ported
/// `Extrema_GenExtPS` (`extrema_surf/gen_ext_ps.rs`, `Extrema_ExtPS.cxx:346`).
/// The grid+Newton `point_surface_newton_all*` helpers are not on this path.
///
/// The previous body was an invented 16×16 grid + hill-climb + golden-section
/// over an expanding window (audit A16).
pub fn project_point_on_surface(s: &dyn Surface, p: &GpPnt, tol: f64) -> Option<PointOnSurface> {
    let tol = if tol > 0.0 { tol } else { occt_core::precision::PCONFUSION };
    let ex = crate::extrema_surf::ExtPs::with_surface(p, s, tol, tol);
    if !ex.is_done() {
        return None;
    }
    let n = ex.nb_ext();
    if n == 0 {
        return None;
    }
    let mut best = 1usize;
    for i in 2..=n {
        if ex.square_distance(i) < ex.square_distance(best) {
            best = i;
        }
    }
    let (u, v, point) = ex.point(best);
    Some(PointOnSurface { u, v, point, distance: point.distance(p) })
}

/// `GeomAPI_ExtremaSurfaceSurface`
/// (`GeomAPI_ExtremaSurfaceSurface.hxx:49-182`,
/// `GeomAPI_ExtremaSurfaceSurface.cxx:26-243`).
///
/// `Init` runs `Extrema_ExtSS(TheSurface1, TheSurface2, Tol, Tol)` with
/// `Tol = Precision::PConfusion()` (`cxx:60-69`) and then keeps the index of
/// the smallest `SquareDistance` in `myIndex` (`cxx:75-89`). `NearestPoints`
/// (`cxx:202-208`), `LowerDistanceParameters` (`cxx:210-221`) and
/// `LowerDistance` (`cxx:223-231`) all read that index; `Points`, `Parameters`
/// and `Distance` answer for the index their caller passes (`cxx:162-200`).
///
/// UNPORTED overloads: the 10-argument `Init` (`cxx:93-145`) and its
/// constructor (`cxx:43-56`), which take explicit parameter windows instead of
/// the adaptor's own; no consumer in this port uses them.
pub struct ExtremaSurfaceSurface {
    /// `myExtSS`.
    ext_ss: ExtremaExtSS,
    /// `myIsDone`: `myExtSS.IsDone() && (myExtSS.NbExt() > 0)` (`cxx:68`).
    is_done: bool,
    /// `myIndex`: the index of the smallest `SquareDistance` (`cxx:79-88`);
    /// `0` while `!myIsDone`, as in the cxx default constructor (`cxx:29-30`).
    index: i32,
}

impl ExtremaSurfaceSurface {
    /// `GeomAPI_ExtremaSurfaceSurface(S1, S2)` (`cxx:37-40`) -> `Init`.
    pub fn new(s1: &dyn Surface, s2: &dyn Surface) -> Self {
        let ext_ss = ExtremaExtSS::new(
            s1,
            s2,
            occt_core::precision::PCONFUSION,
            occt_core::precision::PCONFUSION,
        );
        let is_done = ext_ss.is_done() && (ext_ss.nb_ext() > 0);
        let mut index = 0;
        if is_done {
            // evaluate the lower distance and its index (cxx:75-89)
            let mut dist2_min = ext_ss.square_distance(1);
            index = 1;
            for i in 2..=ext_ss.nb_ext() as i32 {
                let dist2 = ext_ss.square_distance(i);
                if dist2 < dist2_min {
                    dist2_min = dist2;
                    index = i;
                }
            }
        }
        ExtremaSurfaceSurface {
            ext_ss,
            is_done,
            index,
        }
    }

    /// `NbExtrema()` (`cxx:147-157`): `myExtSS.NbExt()` while done, else `0`.
    pub fn nb_extrema(&self) -> i32 {
        if self.is_done {
            self.ext_ss.nb_ext() as i32
        } else {
            0
        }
    }

    /// `IsParallel()` (`hxx:150`).
    pub fn is_parallel(&self) -> bool {
        self.ext_ss.is_parallel()
    }

    /// `Points(Index, P1, P2)` (`cxx:162-172`).
    pub fn points(&self, index: i32) -> (GpPnt, GpPnt) {
        assert!(
            index >= 1 && index <= self.nb_extrema(),
            "GeomAPI_ExtremaSurfaceSurface::Points: Standard_OutOfRange"
        );
        let (ps1, ps2) = self.ext_ss.points(index);
        (ps1.value(), ps2.value())
    }

    /// `Parameters(Index, U1, V1, U2, V2)` (`cxx:174-190`).
    pub fn parameters(&self, index: i32) -> (f64, f64, f64, f64) {
        assert!(
            index >= 1 && index <= self.nb_extrema(),
            "GeomAPI_ExtremaSurfaceSurface::Parameters: Standard_OutOfRange"
        );
        let (ps1, ps2) = self.ext_ss.points(index);
        let (u1, v1) = ps1.parameter();
        let (u2, v2) = ps2.parameter();
        (u1, v1, u2, v2)
    }

    /// `Distance(Index)` (`cxx:192-200`):
    /// `sqrt(myExtSS.SquareDistance(Index))`.
    pub fn distance(&self, index: i32) -> f64 {
        assert!(
            index >= 1 && index <= self.nb_extrema(),
            "GeomAPI_ExtremaSurfaceSurface::Distance: Standard_OutOfRange"
        );
        self.ext_ss.square_distance(index).sqrt()
    }

    /// `NearestPoints(P1, P2)` (`cxx:202-208`): `Points(myIndex, P1, P2)`.
    ///
    /// `StdFail_NotDone` when `!myIsDone` (`cxx:206`); on the parallel-plane
    /// arm the recorded index has no points and `Extrema_ExtSS::Points`
    /// (`Extrema_ExtSS.cxx:292-300`) raises `Standard_OutOfRange`.
    pub fn nearest_points(&self) -> (GpPnt, GpPnt) {
        assert!(
            self.is_done,
            "GeomAPI_ExtremaSurfaceSurface::NearestPoints: StdFail_NotDone"
        );
        self.points(self.index)
    }

    /// `LowerDistanceParameters(U1, V1, U2, V2)` (`cxx:210-221`):
    /// `Parameters(myIndex, ...)`.
    pub fn lower_distance_parameters(&self) -> (f64, f64, f64, f64) {
        assert!(
            self.is_done,
            "GeomAPI_ExtremaSurfaceSurface::LowerDistanceParameters: StdFail_NotDone"
        );
        self.parameters(self.index)
    }

    /// `LowerDistance()` (`cxx:223-231`):
    /// `sqrt(myExtSS.SquareDistance(myIndex))`.
    pub fn lower_distance(&self) -> f64 {
        assert!(
            self.is_done,
            "GeomAPI_ExtremaSurfaceSurface::LowerDistance: StdFail_NotDone"
        );
        self.ext_ss.square_distance(self.index).sqrt()
    }

    /// `Extrema()` (`hxx:176`): the algorithmic object from `Extrema`.
    pub fn extrema(&self) -> &ExtremaExtSS {
        &self.ext_ss
    }
}

/// Distance from point `p` to the nearest point on curve `c`.
pub fn distance_point_curve(c: &dyn Curve, p: &GpPnt) -> f64 {
    project_point_on_curve(c, p, 1e-7).map_or(f64::INFINITY, |pc| pc.distance)
}

/// Distance from point `p` to the nearest point on surface `s`.
pub fn distance_point_surface(s: &dyn Surface, p: &GpPnt) -> f64 {
    project_point_on_surface(s, p, 1e-7).map_or(f64::INFINITY, |ps| ps.distance)
}

/// Unit tangent vector of `c` at `u` (from `d1`, normalized).
pub fn tangent_at(c: &dyn Curve, u: f64) -> GpVec {
    let (_, v) = c.d1(u);
    v.normalized()
}

/// Approximate bounding box of a curve by sampling `(pmin, pmax)`.
/// Unbounded curves are sampled over a clamped `[-1, 1]` proxy window.
pub fn curve_bounding_box(c: &dyn Curve, samples: usize) -> (GpPnt, GpPnt) {
    let (a, b) = finite_range(c.first_parameter(), c.last_parameter());
    let n = samples.max(1);
    let mut min = GpPnt::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut max = GpPnt::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for i in 0..=n {
        let p = c.d0(a + (b - a) * i as f64 / n as f64);
        min.coord.x = min.coord.x.min(p.x());
        min.coord.y = min.coord.y.min(p.y());
        min.coord.z = min.coord.z.min(p.z());
        max.coord.x = max.coord.x.max(p.x());
        max.coord.y = max.coord.y.max(p.y());
        max.coord.z = max.coord.z.max(p.z());
    }
    (min, max)
}

/// Approximate p-curve of `c` on `s`: projection of each curve sample onto the
/// surface, sampled at `samples` parameter values.
pub fn pcurve_of_curve_on_surface(c: &dyn Curve, s: &dyn Surface, samples: usize) -> Vec<GpPnt2d> {
    let (a, b) = finite_range(c.first_parameter(), c.last_parameter());
    let n = samples.max(2);
    let mut out: Vec<GpPnt2d> = Vec::with_capacity(n);
    for i in 0..n {
        let u = a + (b - a) * i as f64 / (n - 1) as f64;
        match project_point_on_surface(s, &c.d0(u), occt_core::precision::PCONFUSION) {
            Some(ps) => out.push(GpPnt2d::new(ps.u, ps.v)),
            // `GeomAPI_ProjectPointOnSurf` reports not-done for this sample;
            // repeat the previous parameter rather than dropping the point.
            // (OCCT's projector retries with `ShapeAnalysis_Surface::ValueOfUV`,
            // `ShapeAnalysis_Surface.cxx:1449-1459`, which the port lacks.)
            None => out.push(out.last().copied().unwrap_or(GpPnt2d::new(0.0, 0.0))),
        }
    }
    out
}

// --- internals -------------------------------------------------------------

/// `ProjLib_Plane::Project` overloads (`ProjLib_Plane.cxx:99-169`): the exact
/// projection of an elementary curve onto a plane. `None` is the arm the
/// OCCT dispatch skips (`Project(ProjLib_Projector&, C)` breaking for
/// `GeomAbs_BSplineCurve` / `GeomAbs_BezierCurve` / `GeomAbs_OtherCurve`,
/// `ProjLib_ProjectedCurve.cxx:262-268`).
///
/// `ProjLib_Plane::EvalPnt2d` / `EvalDir2d` are
/// [`projlib::eval_pln_pnt2d`] / [`projlib::eval_pln_dir2d`].
fn proj_lib_plane_project(c: &dyn Curve, pln: &GpPln) -> Option<Arc<dyn Curve2d>> {
    // `void ProjLib_Plane::Project(const gp_Lin& L)` (`cxx:101-106`).
    if let Some(l) = c.gp_line() {
        let p = projlib::eval_pln_pnt2d(pln, &l.location());
        let (dx, dy) = projlib::eval_pln_dir2d(pln, &l.direction());
        let d = GpDir2d::new(dx, dy).ok()?;
        return Some(Arc::new(Geom2dLine::new(GpAx2d::new(p, d))));
    }
    // `Project(const gp_Circ& C)` (`cxx:110-123`).
    if let Some(cc) = c.gp_circ() {
        let p = projlib::eval_pln_pnt2d(pln, &cc.location());
        let ax = ax22d_on_plane(pln, &cc.position())?;
        return Some(Arc::new(Geom2dCircle::new(
            occt_core::gp::GpCirc2d::new(ax, cc.radius()),
        )));
    }
    // `Project(const gp_Elips& E)` (`cxx:127-139`).
    if let Some(e) = c.gp_ellipse() {
        let ax = ax22d_on_plane(pln, e.position())?;
        return Some(Arc::new(Geom2dEllipse::new(
            occt_core::gp::GpElips2d::new(ax, e.major_radius, e.minor_radius),
        )));
    }
    // `Project(const gp_Parab& P)` (`cxx:143-153`).
    if let Some(p) = c.gp_parabola() {
        let ax = ax22d_on_plane(pln, p.position())?;
        return Some(Arc::new(Geom2dParabola::new(occt_core::gp::GpParab2d::new(
            ax, p.focal,
        ))));
    }
    // `Project(const gp_Hypr& H)` (`cxx:157-168`).
    if let Some(h) = c.gp_hyperbola() {
        let ax = ax22d_on_plane(pln, h.position())?;
        return Some(Arc::new(Geom2dHyperbola::new(
            occt_core::gp::GpHypr2d::new(ax, h.major_radius, h.minor_radius),
        )));
    }
    None
}

/// `gp_Ax22d Ax(EvalPnt2d(C.Location(), myPlane),
/// EvalDir2d(C.Position().XDirection(), myPlane),
/// EvalDir2d(C.Position().YDirection(), myPlane))` (`ProjLib_Plane.cxx:116-119`).
fn ax22d_on_plane(pln: &GpPln, pos: &occt_core::gp::GpAx2) -> Option<GpAx22d> {
    let p = projlib::eval_pln_pnt2d(pln, &pos.location());
    let (xx, xy) = projlib::eval_pln_dir2d(pln, pos.x_direction());
    let (yx, yy) = projlib::eval_pln_dir2d(pln, pos.y_direction());
    let x2d = GpDir2d::new(xx, xy).ok()?;
    let y2d = GpDir2d::new(yx, yy).ok()?;
    GpAx22d::new(p, x2d, y2d).ok()
}

/// `GeomAPI::To2d(const Handle(Geom_Curve)& C, const gp_Pln& P)`
/// (`GeomAPI.cxx:37-52`).
///
/// The projection runs through `ProjLib_ProjectedCurve(S = Geom_Plane(P), C)`:
/// `Perform` dispatches `GeomAbs_Plane` to `ProjLib_Plane` + `Project`
/// (`ProjLib_ProjectedCurve.cxx:391-396`); when that leaves the result not done
/// for a B-spline / Bezier curve, the analytic arm
/// `ProjLib_ComputeApprox::Perform` projects the poles exactly
/// (`ProjLib_ComputeApprox.cxx:1197-1253`). GeomAPI keeps the result only for
/// `GetType() != GeomAbs_OffsetCurve && != GeomAbs_OtherCurve` (`cxx:46-49`).
///
/// UNPORTED arm: the general approximation
/// (`ProjLib_ComputeApprox::Perform`'s `ProjLib_Function` branch,
/// `ProjLib_ComputeApprox.cxx:1254+`) for a non-elementary curve on a plane.
/// Not reachable from `ShapeFix_Wire::RemoveLoop`, whose `To2d` input is always
/// the B-spline that `GeomConvert_CompCurveToBSplineCurve` produced.
pub fn to_2d(c: &dyn Curve, pln: &GpPln) -> Option<Arc<dyn Curve2d>> {
    if let Some(c2d) = proj_lib_plane_project(c, pln) {
        return Some(c2d);
    }
    // `if (CType == GeomAbs_BSplineCurve && SType == GeomAbs_Plane)`
    // (`ProjLib_ComputeApprox.cxx:1197-1227`): project every pole, keep the
    // knots, multiplicities, degree, weights and the periodicity flag.
    if let (Some(poles), Some(degree), Some(knots)) =
        (c.bspline_poles(), c.nurbs_degree(), c.bspline_knots())
    {
        let xs: Vec<f64> = poles.iter().map(|p| projlib::eval_pln_pnt2d(pln, p).x()).collect();
        let ys: Vec<f64> = poles.iter().map(|p| projlib::eval_pln_pnt2d(pln, p).y()).collect();
        let weights = c.bspline_weights().map(|w| w.to_vec());
        let bs = Geom2dBSplineCurve::from_flat(
            xs,
            ys,
            weights,
            knots.to_vec(),
            degree,
            // `myResult.SetPeriodic()` (`ProjLib_ProjectedCurve.cxx:781-784`).
            c.is_periodic(),
        )
        .ok()?;
        return Some(Arc::new(bs));
    }
    // `else if (CType == GeomAbs_BezierCurve && SType == GeomAbs_Plane)`
    // (`ProjLib_ComputeApprox.cxx:1229-1253`).
    if let Some(poles) = c.bezier_poles() {
        let poles2d: Vec<GpPnt2d> = poles.iter().map(|p| projlib::eval_pln_pnt2d(pln, p)).collect();
        let bez = Geom2dBezierCurve::new(poles2d).ok()?;
        return Some(Arc::new(bez));
    }
    None
}

/// Finite, sane sampling bounds; unbounded ranges clamp to `[-1, 1]`.
fn finite_range(a: f64, b: f64) -> (f64, f64) {
    if a.is_finite() && b.is_finite() && b > a {
        (a, b)
    } else {
        (-1.0, 1.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::circle::GeomCircle;
    use crate::line::GeomLine;
    use crate::plane::GeomPlane;
    use occt_core::gp::{GpAx2, GpAx3, GpCirc, GpDir, GpPln};

    #[test]
    fn project_point_onto_circle() {
        let circle = GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0));
        let pc = project_point_on_curve(&circle as &dyn Curve, &GpPnt::new(0.0, 2.0, 0.0), 1e-6).unwrap();
        assert!((pc.distance - 1.0).abs() < 1e-6, "distance={}", pc.distance);
        assert!((pc.parameter - std::f64::consts::FRAC_PI_2).abs() < 1e-4, "parameter={}", pc.parameter);
        assert!(pc.point.distance(&GpPnt::new(0.0, 1.0, 0.0)) < 1e-6);
    }

    #[test]
    fn project_point_onto_plane() {
        let plane = GeomPlane::new(GpPln::new(GpAx3::standard()));
        let ps = project_point_on_surface(&plane as &dyn Surface, &GpPnt::new(0.5, 0.5, 1.0), 1e-6).unwrap();
        assert!((ps.distance - 1.0).abs() < 1e-6, "distance={}", ps.distance);
        assert!(ps.point.distance(&GpPnt::new(0.5, 0.5, 0.0)) < 1e-6);
        // Far point on the unbounded plane forces window expansion.
        let ps2 = project_point_on_surface(&plane as &dyn Surface, &GpPnt::new(10.0, 10.0, 1.0), 1e-6).unwrap();
        assert!((ps2.distance - 1.0).abs() < 1e-6, "distance={}", ps2.distance);
        assert!(ps2.point.distance(&GpPnt::new(10.0, 10.0, 0.0)) < 1e-6);
    }

    #[test]
    fn distance_and_tangent() {
        let circle = GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0));
        let d = distance_point_curve(&circle as &dyn Curve, &GpPnt::new(3.0, 0.0, 0.0));
        assert!((d - 2.0).abs() < 1e-6, "d={d}");
        let t = tangent_at(&circle as &dyn Curve, 0.0);
        assert!(t.subtracted(&GpVec::new(0.0, 1.0, 0.0)).magnitude() < 1e-6, "tangent={t:?}");
    }

    #[test]
    fn bbox_and_pcurve() {
        let circle = GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0));
        let (pmin, pmax) = curve_bounding_box(&circle as &dyn Curve, 64);
        assert!((pmin.x() + 1.0).abs() < 1e-6 && (pmax.x() - 1.0).abs() < 1e-6);

        let plane = GeomPlane::new(GpPln::new(GpAx3::standard()));
        let line = GeomLine::from_pnt_dir(GpPnt::new(0.5, 0.5, 0.0), GpDir::new(1.0, 0.0, 0.0).unwrap());
        let pc = pcurve_of_curve_on_surface(&line as &dyn Curve, &plane as &dyn Surface, 4);
        assert_eq!(pc.len(), 4);
        // Unbounded line sampled over the [-1, 1] proxy window: first/last
        // samples land at (0.5, 0.5, 0) ± (1, 0, 0) → (u, v) = (x, y).
        assert!((pc[0].x() - -0.5).abs() < 1e-6 && (pc[0].y() - 0.5).abs() < 1e-6);
        assert!((pc[3].x() - 1.5).abs() < 1e-6 && (pc[3].y() - 0.5).abs() < 1e-6);
    }
}
