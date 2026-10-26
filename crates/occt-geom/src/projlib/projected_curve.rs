//! `ProjLib_ProjectedCurve` — project a 3D curve into a surface's `(u, v)`
//! parameter space.
//!
//! Source: `ModelingData/TKGeomBase/ProjLib/ProjLib_ProjectedCurve.hxx:27-232`,
//! `.cxx:25-1178` (OCCT 8.0.0).
//!
//! `Perform` (`cxx:369-886`) dispatches on `Adaptor3d_Surface::GetType()`:
//! the five elementary surfaces go through the `ProjLib_Projector` family and
//! keep an analytic 2D result (line / circle / ellipse / parabola /
//! hyperbola); every other type goes through an approximator.
//!
//! **UNPORTED** (approximators, same file):
//! - `ProjLib_ComputeApproxOnPolarSurface` (`ProjLib_ComputeApproxOnPolarSurface.cxx`,
//!   2214 lines) — the `GeomAbs_BezierSurface` / `GeomAbs_BSplineSurface` arm
//!   (`cxx:496-536`). The pole trimming of `:466-493` and the `ExtendC2d`
//!   tail (`:507-535`) are ported around it.
//! - `Approx_CurveOnSurface` + `ProjLib_HCompProjectedCurve`
//!   (`ProjLib_CompProjectedCurve.cxx`, 2069 lines) — the `default:` arm
//!   (`cxx:539-715`), including the `Extrema_ExtPC` search for the singular
//!   parameter on a surface of revolution's basis curve (`cxx:558-610`).
//! - `ExtendC2d` (`cxx:162-238`) and the `GeomLib::SameRange` remake of a
//!   trimmed approximation result (`cxx:740-766`); the `ProjLib_ComputeApprox`
//!   fallback itself is ported ([`crate::projlib::ComputeApprox`]).

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d};
use occt_core::kernel::geomabs::{CurveType, SurfaceType};
use occt_core::precision;
use occt_geom2d::bezier_curve::Geom2dBezierCurve;
use occt_geom2d::bspline_curve::Geom2dBSplineCurve;
use occt_geom2d::curve::Curve2d;

use super::{cone, cylinder, plane, projector::ProjectorState, sphere, torus, Projector};
use crate::curve::Curve;
use crate::surface::Surface;
use crate::trimmed::GeomTrimmedCurveBasis;

/// `AppParCurves_Constraint` (`AppParCurves_Constraint.hxx`): the boundary
/// constraint handed to the approximators. Only read by the UNPORTED
/// approximation arms; the value is kept so the accessors stay faithful.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Constraint {
    NoConstraint,
    PassPoint,
    TangencyPoint,
}

/// `Adaptor3d_Surface::GetType()` (`GeomAdaptor_Surface.cxx:480-...`): the
/// dynamic-type test order of `GeomAdaptor_Surface::Load`. The elementary
/// surfaces are recognized through the trait accessors that already answer
/// `Adaptor3d_Surface::Plane()` / `Cylinder()` / ... , then the parametric
/// patches, the swept surfaces and the offset.
pub fn surface_type(s: &dyn Surface) -> SurfaceType {
    if s.gp_pln().is_some() {
        SurfaceType::Plane
    } else if s.gp_cylinder().is_some() {
        SurfaceType::Cylinder
    } else if s.gp_cone().is_some() {
        SurfaceType::Cone
    } else if s.gp_sphere().is_some() {
        SurfaceType::Sphere
    } else if s.gp_torus().is_some() {
        SurfaceType::Torus
    } else if s.is_bezier_surface() {
        SurfaceType::BezierSurface
    } else if s.is_bspline_surface() {
        SurfaceType::BSplineSurface
    } else if s.is_surface_of_revolution() {
        SurfaceType::SurfaceOfRevolution
    } else if s.is_surface_of_linear_extrusion() {
        SurfaceType::SurfaceOfExtrusion
    } else if s.is_offset_surface() {
        SurfaceType::OffsetSurface
    } else {
        SurfaceType::OtherSurface
    }
}

/// `Adaptor3d_Curve::GetType()` (`GeomAdaptor_Curve::GetType`,
/// `GeomAdaptor_Curve.cxx:314-360`): the analytic overloads tested before the
/// "parametric" answer, in the `GeomAdaptor_Curve::load` order.
pub fn curve_type(c: &dyn Curve) -> CurveType {
    if c.gp_line().is_some() {
        CurveType::Line
    } else if c.gp_circ().is_some() {
        CurveType::Circle
    } else if c.gp_ellipse().is_some() {
        CurveType::Ellipse
    } else if c.gp_hyperbola().is_some() {
        CurveType::Hyperbola
    } else if c.gp_parabola().is_some() {
        CurveType::Parabola
    } else if c.bezier_poles().is_some() {
        CurveType::BezierCurve
    } else if c.bspline_poles().is_some() {
        CurveType::BSplineCurve
    } else if c.offset_curve().is_some() {
        CurveType::OffsetCurve
    } else {
        CurveType::OtherCurve
    }
}

/// `ComputeTolU` (`cxx:53-61`): the surface's U resolution, capped at a
/// hundredth of the U period on a periodic surface.
pub fn compute_tol_u(s: &dyn Surface, tolerance: f64) -> Option<f64> {
    let (mut tol_u, _) = s.uv_resolution(tolerance)?;
    if s.is_u_periodic() {
        tol_u = tol_u.min(0.01 * s.u_period());
    }
    Some(tol_u)
}

/// `ComputeTolV` (`cxx:64-74`).
pub fn compute_tol_v(s: &dyn Surface, tolerance: f64) -> Option<f64> {
    let (_, mut tol_v) = s.uv_resolution(tolerance)?;
    if s.is_v_periodic() {
        tol_v = tol_v.min(0.01 * s.v_period());
    }
    Some(tol_v)
}

/// `IsoIsDeg(S, Param, IT, TolMin, TolMax)` (`cxx:80-126`): an iso of the
/// surface is degenerate when the norm of the *other* first derivative stays
/// inside `[TolMin, TolMax]` over the whole range, sampled in ten steps.
pub fn iso_is_deg(s: &dyn Surface, param: f64, iso_v: bool, tol_min: f64, tol_max: f64) -> bool {
    let (u1, u2) = s.u_range();
    let (v1, v2) = s.v_range();
    let mut d1_norm_max = 0.0f64;
    if iso_v {
        // `GeomAbs_IsoV`: vary U at constant V, watch |D1U|.
        let step = (u2 - u1) / 10.0;
        let mut t = u1;
        while t <= u2 {
            let (_, du, _) = s.d1(t, param);
            d1_norm_max = d1_norm_max.max(du.magnitude());
            t += step;
        }
    } else {
        // `GeomAbs_IsoU`: vary V at constant U, watch |D1V|.
        let step = (v2 - v1) / 10.0;
        let mut t = v1;
        while t <= v2 {
            let (_, _, dv) = s.d1(param, t);
            d1_norm_max = d1_norm_max.max(dv.magnitude());
            t += step;
        }
    }
    !(d1_norm_max > tol_max || d1_norm_max < tol_min)
}

/// `TrimC3d(myCurve, IsTrimmed, dt, Pole, SingularCase, NumberOfSingularCase,
/// TolConf)` (`cxx:129-159`): shorten the curve by `dt` at an end that sits on
/// the degenerate pole, recording which end and which singular case.
pub fn trim_c3d(
    curve: &mut Arc<dyn Curve>,
    is_trimmed: &mut [bool; 2],
    dt: f64,
    pole: &GpPnt,
    singular_case: &mut [i32; 2],
    number_of_singular_case: i32,
    tol_conf: f64,
) {
    let mut f = curve.first_parameter();
    let mut l = curve.last_parameter();

    let p = curve.d0(f);
    if p.distance(pole) <= tol_conf {
        is_trimmed[0] = true;
        f += dt;
        *curve = Arc::new(GeomTrimmedCurveBasis::new(
            curve.clone(),
            f,
            l,
        ));
        singular_case[0] = number_of_singular_case;
    }

    let p = curve.d0(l);
    if p.distance(pole) <= tol_conf {
        is_trimmed[1] = true;
        l -= dt;
        *curve = Arc::new(GeomTrimmedCurveBasis::new(
            curve.clone(),
            f,
            l,
        ));
        singular_case[1] = number_of_singular_case;
    }
}

/// `ExtendC2d` (`cxx:162-238`) — **UNPORTED**: it is only reached from the
/// three approximation arms above (singular-pole extensions of the fitted
/// B-spline). It needs `Geom2dConvert_CompCurveToBSplineCurve` wiring into
/// `GeomLib::SameRange`, neither of which is reachable today.

/// `Project(ProjLib_Projector&, Adaptor3d_Curve&)` (`cxx:241-270`): the
/// double-dispatch on `Adaptor3d_Curve::GetType()`. Parametric / offset /
/// other curves fall through without touching the projector, which OCCT
/// does by `break` (`cxx:262-265`) so the caller keeps an undone projector
/// and takes its approximation branch.
pub fn project<P: Projector + ?Sized>(p: &mut P, c: &dyn Curve) {
    match curve_type(c) {
        CurveType::Line => p.project_lin(&c.gp_line().expect("GeomAbs_Line")),
        CurveType::Circle => p.project_circ(&c.gp_circ().expect("GeomAbs_Circle")),
        CurveType::Ellipse => p.project_elips(&c.gp_ellipse().expect("GeomAbs_Ellipse")),
        CurveType::Hyperbola => p.project_hypr(&c.gp_hyperbola().expect("GeomAbs_Hyperbola")),
        CurveType::Parabola => p.project_parab(&c.gp_parabola().expect("GeomAbs_Parabola")),
        CurveType::BSplineCurve
        | CurveType::BezierCurve
        | CurveType::OffsetCurve
        | CurveType::OtherCurve => {}
    }
}

/// `class ProjLib_ProjectedCurve` (`hxx:27-232`).
///
/// The C++ members map to:
/// `mySurface` → [`ProjectedCurve::surface`], `myCurve` →
/// [`ProjectedCurve::curve`] (replaced by `Trim` in `Perform`), `myResult` →
/// [`ProjectedCurve::result`] (the `ProjLib_Projector` base state),
/// `myTolerance` → [`ProjectedCurve::tolerance`], and
/// `myDegMin` / `myDegMax` / `myMaxSegments` / `myMaxDist` / `myBndPnt` are
/// the approximation settings kept by the `Set*` accessors.
pub struct ProjectedCurve {
    /// `mySurface` (`hxx:212`).
    pub surface: Option<Arc<dyn Surface>>,
    /// `myCurve` (`hxx:213`), after any `TrimC3d`.
    pub curve: Option<Arc<dyn Curve>>,
    /// `myResult` (`hxx:214`), the `ProjLib_Projector` state.
    pub result: ProjectorState,
    /// `myTolerance` (`hxx:215`).
    tolerance: f64,
    /// `myDegMin` / `myDegMax` (`hxx:216-217`).
    deg_min: i32,
    deg_max: i32,
    /// `myMaxSegments` (`hxx:218`).
    max_segments: i32,
    /// `myMaxDist` (`hxx:219`).
    max_dist: f64,
    /// `myBndPnt` (`hxx:220`).
    bnd_pnt: Constraint,
}

impl ProjectedCurve {
    /// `ProjLib_ProjectedCurve()` (`cxx:274-282`).
    pub fn new() -> Self {
        Self {
            surface: None,
            curve: None,
            result: ProjectorState::new(),
            tolerance: precision::CONFUSION,
            deg_min: -1,
            deg_max: -1,
            max_segments: -1,
            max_dist: -1.0,
            bnd_pnt: Constraint::TangencyPoint,
        }
    }

    /// `ProjLib_ProjectedCurve(S)` (`cxx:286-295`).
    pub fn with_surface(s: Arc<dyn Surface>) -> Self {
        let mut p = Self::new();
        p.load_surface(s);
        p
    }

    /// `ProjLib_ProjectedCurve(S, C)` (`cxx:299-310`).
    pub fn with_surface_curve(s: Arc<dyn Surface>, c: Arc<dyn Curve>) -> Self {
        let mut p = Self::with_surface(s);
        p.perform(c);
        p
    }

    /// `ProjLib_ProjectedCurve(S, C, Tol)` (`cxx:314-327`): the tolerance is
    /// raised to `Precision::Confusion()`.
    pub fn with_surface_curve_tol(s: Arc<dyn Surface>, c: Arc<dyn Curve>, tol: f64) -> Self {
        let mut p = Self::with_surface(s);
        p.load_tolerance(tol.max(precision::CONFUSION));
        p.perform(c);
        p
    }

    /// `Load(S)` (`cxx:355-358`).
    pub fn load_surface(&mut self, s: Arc<dyn Surface>) {
        self.surface = Some(s);
    }

    /// `Load(theTol)` (`cxx:362-365`).
    pub fn load_tolerance(&mut self, tol: f64) {
        self.tolerance = tol;
    }

    /// `SetDegree(theDegMin, theDegMax)` (`cxx:890-895`).
    pub fn set_degree(&mut self, deg_min: i32, deg_max: i32) {
        self.deg_min = deg_min;
        self.deg_max = deg_max;
    }

    /// `SetMaxSegments(theMaxSegments)` (`cxx:898-902`).
    pub fn set_max_segments(&mut self, max_segments: i32) {
        self.max_segments = max_segments;
    }

    /// `SetBndPnt(theBndPnt)` (`cxx:905-908`).
    pub fn set_bnd_pnt(&mut self, bnd_pnt: Constraint) {
        self.bnd_pnt = bnd_pnt;
    }

    /// `SetMaxDist(theMaxDist)` (`cxx:912-915`).
    pub fn set_max_dist(&mut self, max_dist: f64) {
        self.max_dist = max_dist;
    }

    /// `GetTolerance()` (`cxx:933-936`).
    pub fn get_tolerance(&self) -> f64 {
        self.tolerance
    }

    /// `GetType()` (`cxx:1047-1051`): the type of the projector result.
    pub fn get_type(&self) -> CurveType {
        self.result.curve_type
    }

    /// `Line()` (`cxx:1054-1057`) delegating to `ProjLib_Projector::Line()`
    /// (`ProjLib_Projector.cxx:88-96`): OCCT raises `Standard_NoSuchObject`
    /// when the projector type is not `GeomAbs_Line`.
    pub fn line(&self) -> Result<occt_core::gp::GpLin2d, String> {
        if self.result.curve_type != CurveType::Line {
            return Err("Standard_NoSuchObject: ProjLib_Projector::Line".into());
        }
        Ok(self.result.lin)
    }

    /// `Bezier()` (`cxx:1157-1160`) delegating to
    /// `ProjLib_Projector::Bezier()` (`ProjLib_Projector.cxx:234-237`), which
    /// returns `myBezier` with no type check — hence `Option`, not `Result`.
    pub fn bezier(&self) -> Option<Arc<Geom2dBezierCurve>> {
        self.result.bezier.clone()
    }

    /// `BSpline()` (`cxx:1164-1167`) delegating to
    /// `ProjLib_Projector::BSpline()` (`ProjLib_Projector.cxx:245-251`), also
    /// an unguarded `myBSpline` reader.
    pub fn bspline(&self) -> Option<Arc<Geom2dBSplineCurve>> {
        self.result.bspline.clone()
    }

    /// `IsDone()` of the `myResult` projector — the port's stand-in for
    /// `!myResult.IsNull()` / `myResult.GetType()`, since C++ has no
    /// `IsDone` on the base. `false` until an arm sets a result.
    pub fn is_done(&self) -> bool {
        self.result.done
    }

    /// `FirstParameter()` (`cxx:940-943`).
    pub fn first_parameter(&self) -> f64 {
        self.curve
            .as_ref()
            .map(|c| c.first_parameter())
            .unwrap_or(0.0)
    }

    /// `LastParameter()` (`cxx:947-950`).
    pub fn last_parameter(&self) -> f64 {
        self.curve
            .as_ref()
            .map(|c| c.last_parameter())
            .unwrap_or(0.0)
    }

    /// `Perform(C)` (`cxx:369-886`).
    pub fn perform(&mut self, c: Arc<dyn Curve>) {
        self.tolerance = self.tolerance.max(precision::CONFUSION);
        let Some(surface) = self.surface.clone() else {
            return;
        };
        let mut curve = c.clone();
        self.curve = Some(c);

        let first_par = curve.first_parameter();
        let last_par = curve.last_parameter();
        let s_type = surface_type(surface.as_ref());
        // `CType` (`cxx:376`), read by the `ProjLib_ComputeApprox` fallback
        // (`cxx:769`, `:781`).
        let c_type = curve_type(curve.as_ref());
        let mut is_analytical_surf = true;
        let mut is_trimmed = [false, false];
        let mut singular_case = [0i32, 0i32];
        let eps = 0.01;
        let mut tol_conf = precision::CONFUSION;
        let mut dt = (last_par - first_par) * eps;
        let (u1, u2) = surface.u_range();
        let (v1, v2) = surface.v_range();

        match s_type {
            SurfaceType::Plane => {
                let mut p = plane::Plane::from_plane(&surface.gp_pln().expect("GeomAbs_Plane"));
                project(&mut p, curve.as_ref());
                self.result = p.state().clone();
            }
            SurfaceType::Cylinder => {
                let mut p =
                    cylinder::Cylinder::from_cylinder(&surface.gp_cylinder().expect("GeomAbs_Cylinder"));
                project(&mut p, curve.as_ref());
                self.result = p.state().clone();
            }
            SurfaceType::Cone => {
                let mut p = cone::Cone::from_cone(&surface.gp_cone().expect("GeomAbs_Cone"));
                project(&mut p, curve.as_ref());
                self.result = p.state().clone();
            }
            SurfaceType::Sphere => {
                let mut p = sphere::Sphere::from_sphere(&surface.gp_sphere().expect("GeomAbs_Sphere"));
                project(&mut p, curve.as_ref());
                if p.is_done() {
                    // Place into the pseudo-period (the sphere is not periodic
                    // in V), `cxx:415-419`.
                    p.set_in_bounds(curve.first_parameter());
                } else {
                    // `cxx:421-443`: the circle is not a meridian/parallel, so
                    // trim the ends that sit on the pole before approximating.
                    let vmax = std::f64::consts::FRAC_PI_2;
                    let vmin = -vmax;
                    let minang = 1.0e-5 * std::f64::consts::PI;
                    let sph = surface.gp_sphere().expect("GeomAbs_Sphere");
                    let r = sph.radius;
                    let f = curve.first_parameter();
                    let l = curve.last_parameter();
                    let pf = curve.d0(f);
                    let pl = curve.d0(l);
                    let a_loc = sph.pos.location();
                    let maxdist = pf.distance(&a_loc).max(pl.distance(&a_loc));
                    tol_conf = (r * minang).max((r - maxdist).abs());

                    let mut pole = surface.d0(u1, vmin);
                    trim_c3d(
                        &mut curve,
                        &mut is_trimmed,
                        dt,
                        &pole,
                        &mut singular_case,
                        3,
                        tol_conf,
                    );
                    pole = surface.d0(u1, vmax);
                    trim_c3d(
                        &mut curve,
                        &mut is_trimmed,
                        dt,
                        &pole,
                        &mut singular_case,
                        4,
                        tol_conf,
                    );
                    self.curve = Some(curve.clone());
                }
                self.result = p.state().clone();
            }
            SurfaceType::Torus => {
                let mut p = torus::Torus::from_torus(&surface.gp_torus().expect("GeomAbs_Torus"));
                project(&mut p, curve.as_ref());
                self.result = p.state().clone();
            }
            SurfaceType::BezierSurface | SurfaceType::BSplineSurface => {
                is_analytical_surf = false;
                let f = curve.first_parameter();
                let l = curve.last_parameter();
                dt = (l - f) * eps;

                // `cxx:469-493`: trim the curve where an iso edge is
                // degenerate (a pole of the patch).
                if iso_is_deg(surface.as_ref(), u1, false, 0.0, self.tolerance) {
                    let pole = surface.d0(u1, v1);
                    trim_c3d(
                        &mut curve,
                        &mut is_trimmed,
                        dt,
                        &pole,
                        &mut singular_case,
                        1,
                        tol_conf,
                    );
                }
                if iso_is_deg(surface.as_ref(), u2, false, 0.0, self.tolerance) {
                    let pole = surface.d0(u2, v1);
                    trim_c3d(
                        &mut curve,
                        &mut is_trimmed,
                        dt,
                        &pole,
                        &mut singular_case,
                        2,
                        tol_conf,
                    );
                }
                if iso_is_deg(surface.as_ref(), v1, true, 0.0, self.tolerance) {
                    let pole = surface.d0(u1, v1);
                    trim_c3d(
                        &mut curve,
                        &mut is_trimmed,
                        dt,
                        &pole,
                        &mut singular_case,
                        3,
                        tol_conf,
                    );
                }
                if iso_is_deg(surface.as_ref(), v2, true, 0.0, self.tolerance) {
                    let pole = surface.d0(u1, v2);
                    trim_c3d(
                        &mut curve,
                        &mut is_trimmed,
                        dt,
                        &pole,
                        &mut singular_case,
                        4,
                        tol_conf,
                    );
                }
                self.curve = Some(curve.clone());

                // UNPORTED: `ProjLib_ComputeApproxOnPolarSurface`
                // (`ProjLib_ComputeApproxOnPolarSurface.cxx`) — `cxx:496-536`.
                // OCCT leaves `myResult` undone when the fit produces no
                // B-spline, which is what this arm returns.
            }
            _ => {
                // UNPORTED: `cxx:539-715`. The singular-parameter search on a
                // surface of revolution (`Extrema_ExtPC` over
                // `mySurface->BasisCurve()`, `cxx:558-610`), the
                // `ProjLib_HCompProjectedCurve` boundary search and the
                // `Approx_CurveOnSurface` fit (`cxx:613-713`).
                self.curve = Some(curve);
                return;
            }
        }

        if !self.result.done && is_analytical_surf {
            // `cxx:717-785`: use the advanced analytical projector when the
            // base analytical projection failed.
            let mut comp = crate::projlib::ComputeApprox::new();
            comp.set_tolerance(self.tolerance);
            comp.set_degree(self.deg_min, self.deg_max);
            comp.set_max_segments(self.max_segments);
            comp.set_bnd_pnt(self.bnd_pnt);
            comp.perform(curve.as_ref(), surface.as_ref());
            if comp.bezier().is_none() && comp.bspline().is_none() {
                // the advanced projector has failed too (`cxx:727-730`).
                self.curve = Some(curve);
                return;
            }
            if is_trimmed[0] || is_trimmed[1] {
                // UNPORTED: `cxx:740-766` — the `ExtendC2d` segment on the
                // trimmed end, then `GeomLib::SameRange` back onto
                // `[FirstPar, LastPar]`. Both run on a `Geom2d_BSplineCurve`
                // through `Geom2dConvert_CompCurveToBSplineCurve`; the port
                // leaves the projector undone, so the caller sees no pcurve
                // for a curve trimmed at a singular pole.
                self.curve = Some(curve);
                return;
            }
            self.result.done = true;
            if s_type == SurfaceType::Plane && c_type == CurveType::BezierCurve {
                // `cxx:769-774`.
                self.result.curve_type = CurveType::BezierCurve;
                self.result.bezier = comp.bezier();
            } else {
                // `cxx:775-778`.
                self.result.curve_type = CurveType::BSplineCurve;
                self.result.bspline = comp.bspline();
            }
            // `cxx:781-784`: set the periodicity flag.
            if s_type == SurfaceType::Plane && c_type == CurveType::BSplineCurve && curve.is_periodic()
            {
                self.result.periodic = true;
            }
            self.tolerance = comp.tolerance();
        }

        // `cxx:787-886`: keep a periodic result inside the surface's parameter
        // range. `GeomAbs_Line` results move through `UFrame` / `VFrame`;
        // B-spline results are translated by the dominant period multiple.
        self.clamp_periodic_result();
        self.curve = Some(curve);
    }

    /// The `cxx:787-886` tail of `Perform`.
    fn clamp_periodic_result(&mut self) {
        let Some(surface) = self.surface.clone() else {
            return;
        };
        let is_periodic = [surface.is_u_periodic(), surface.is_v_periodic()];
        if !self.result.done || !(is_periodic[0] || is_periodic[1]) {
            return;
        }
        let (uf, vf) = (surface.u_range().0, surface.v_range().0);
        let surf_first_par = [uf, vf];
        let mut surf_period = [0.0f64, 0.0];
        if is_periodic[0] {
            surf_period[0] = surface.u_period();
        }
        if is_periodic[1] {
            surf_period[1] = surface.v_period();
        }

        for an_idx in 1..=2usize {
            if !is_periodic[an_idx - 1] {
                continue;
            }
            if self.result.curve_type == CurveType::BSplineCurve {
                let Some(a_res) = self.result.bspline.clone() else {
                    continue;
                };
                let a_deg = a_res.degree() as i32;
                let mut map: std::collections::BTreeMap<i32, i32> = std::collections::BTreeMap::new();
                let mut knot_idx = a_res.first_u_knot_index();
                let last_knot = a_res.last_u_knot_index();
                while knot_idx < last_knot {
                    let a_first = a_res.knot(knot_idx as usize);
                    let a_last = a_res.knot(knot_idx as usize + 1);
                    for int_idx in 0..=a_deg {
                        let curr = a_first
                            + (a_last - a_first) * (int_idx as f64) / (a_deg as f64 + 1.0);
                        let p = a_res.d0(curr);
                        let coord = if an_idx == 1 { p.x() } else { p.y() };
                        let mut key = ((coord - surf_first_par[an_idx - 1])
                            / surf_period[an_idx - 1])
                            as i32;
                        if coord - surf_first_par[an_idx - 1] < 0.0 {
                            key -= 1;
                        }
                        *map.entry(key).or_insert(0) += 1;
                    }
                    knot_idx += 1;
                }
                let mut max_points = 0;
                let mut max_idx = 0;
                for (k, v) in map {
                    if v > max_points {
                        max_points = v;
                        max_idx = k;
                    }
                }
                if max_idx != 0 {
                    let first_pnt = a_res.d0(a_res.first_parameter());
                    let mut second = first_pnt;
                    if an_idx == 1 {
                        second.set_x(first_pnt.x() - surf_period[0] * max_idx as f64);
                    } else {
                        second.set_y(first_pnt.y() - surf_period[1] * max_idx as f64);
                    }
                    let mut trsf = occt_core::gp::GpTrsf2d::identity();
                    trsf.set_translation_vec(&occt_core::gp::GpVec2d::new(
                        second.x() - first_pnt.x(),
                        second.y() - first_pnt.y(),
                    ));
                    let mut moved = (*a_res).clone();
                    moved.transform(&trsf);
                    self.result.bspline = Some(Arc::new(moved));
                }
            }

            if self.result.curve_type == CurveType::Line {
                let t1 = self.first_parameter();
                let t2 = self.last_parameter();
                if an_idx == 1 {
                    self.result.u_frame(t1, t2, surf_first_par[0], surf_period[0]);
                } else {
                    self.result.v_frame(t1, t2, surf_first_par[1], surf_period[1]);
                }
            }
        }
    }
}

impl Default for ProjectedCurve {
    fn default() -> Self {
        Self::new()
    }
}

/// `ProjLib_Projector` result sampling helper: `D0(t)` of the analytic 2D
/// result the projector produced, using the projection `Perform` stored.
///
/// This mirrors `ProjLib_Projector`'s `Line()` / `Circle()` / ... readers and
/// `ProjLib_ProjectedCurve::Value` (`cxx:998-1003`), which evaluate the
/// analytic form when one exists. `None` when no analytic result is stored.
pub fn projector_value(state: &ProjectorState, t: f64) -> Option<GpPnt2d> {
    use occt_core::elib::clib;
    match state.curve_type {
        CurveType::Line => Some(clib::line2d_value(&state.lin, t)),
        CurveType::Circle => Some(clib::circle2d_value(&state.circ, t)),
        _ => None,
    }
}

/// `Geom2d_BSplineCurve` reference kept for the `ExtendC2d` tail; unused while
/// that helper is UNPORTED.
#[allow(dead_code)]
fn _keep(bs: &Geom2dBSplineCurve) -> usize {
    bs.nb_poles()
}

/// `Arc<dyn Curve2d>` clone helper for the `C2D->Copy()` uses.
#[allow(dead_code)]
pub(crate) fn copy_curve2d(c: &dyn Curve2d) -> Arc<dyn Curve2d> {
    Arc::from(c.clone_dyn())
}
