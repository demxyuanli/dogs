//! `ProjLib_ComputeApprox` — the advanced projector that fits a 2d B-spline to
//! the parameter-space image of a 3d curve on an analytic surface.
//!
//! Source: `ModelingData/TKGeomBase/ProjLib/ProjLib_ComputeApprox.cxx` and
//! `.hxx` (OCCT 8.0.0). `Perform` (`cxx:1171-1493`) has three arms:
//!
//! - a `GeomAbs_BSplineCurve` on a `GeomAbs_Plane` (`cxx:1197-1227`) projects
//!   the poles directly and rebuilds the 2d B-spline;
//! - a `GeomAbs_BezierCurve` on a `GeomAbs_Plane` (`cxx:1228-1252`) does the
//!   same for a 2d Bezier;
//! - everything else (`cxx:1253-1492`) feeds [`ProjLibFunction`] to
//!   [`FitAndDivide2d`] (`Approx_ComputeCLine.gxx` instantiated for
//!   `AppCont_Function`), glues the segments with
//!   [`CompBezierCurves2d`](occt_core::convert::CompBezierCurves2d), reparametrizes,
//!   and finally translates / mirrors the result so its midpoint matches the
//!   direct projection of the 3d midpoint ("return curve home").
//!
//! The file statics of the `.cxx` — `PlaneProjector`, `IsEqual`,
//! `Function_Value`, `Function_D1`, `Function_ComputeStep` — live in
//! [`function`]; `Function_SetUVBounds` in [`uv_bounds`]; the `AppCont_Function`
//! adapter in [`projlib_function`].
//!
//! `myTolerance` is rescaled by `aNewTol2d / aTol2d` at `cxx:1417`, so a caller
//! must read [`ComputeApprox::tolerance`] after `Perform` (that is what
//! `ProjLib_ProjectedCurve` does at `ProjLib_ProjectedCurve.cxx:786`).

mod function;
mod projlib_function;
mod uv_bounds;

pub use projlib_function::ProjLibFunction;
pub use uv_bounds::UvBounds;

use std::sync::Arc;

use occt_core::bspl::{banded_interp, knots as bspl_knots};
use occt_core::elib::{clib, slib};
use occt_core::gp::{GpAx2d, GpDir2d, GpPnt2d, GpTrsf2d, GpVec2d};
use occt_core::kernel::geomabs::{CurveType, SurfaceType};
use occt_core::precision::{CONFUSION, PAPPROXIMATION, PCONFUSION};
use occt_core::convert::CompBezierCurves2d;
use occt_geom2d::bezier_curve::Geom2dBezierCurve;
use occt_geom2d::bspline_curve::Geom2dBSplineCurve;
use occt_geom2d::curve::Curve2d;

use crate::app_cont::function::Function;
use crate::app_par_curves::constraint::Constraint;
use crate::approx::fit_and_divide::FitAndDivide2d;
use crate::curve::Curve;
use crate::projlib::projected_curve::{compute_tol_u, compute_tol_v, curve_type, surface_type, Constraint as ProjConstraint};
use crate::projlib::is_ana_surf;
use crate::surface::Surface;

use function::PlaneProjector;

const PI: f64 = std::f64::consts::PI;

/// `class ProjLib_ComputeApprox` (`ProjLib_ComputeApprox.hxx`).
pub struct ComputeApprox {
    /// `myBSpline` (`hxx`).
    bspline: Option<Arc<Geom2dBSplineCurve>>,
    /// `myBezier` (`hxx`).
    bezier: Option<Arc<Geom2dBezierCurve>>,
    /// `myTolerance` (`hxx`).
    tolerance: f64,
    /// `myDegMin` / `myDegMax` (`hxx`).
    deg_min: i32,
    deg_max: i32,
    /// `myMaxSegments` (`hxx`).
    max_segments: i32,
    /// `myBndPnt` (`hxx`).
    bnd_pnt: ProjConstraint,
}

impl ComputeApprox {
    /// `ProjLib_ComputeApprox()` (`cxx:1146-1154`).
    pub fn new() -> Self {
        Self {
            bspline: None,
            bezier: None,
            tolerance: PAPPROXIMATION,
            deg_min: -1,
            deg_max: -1,
            max_segments: -1,
            bnd_pnt: ProjConstraint::TangencyPoint,
        }
    }

    /// `ProjLib_ComputeApprox(C, S, Tol)` (`cxx:1157-1168`): the tolerance is
    /// raised to `Precision::PApproximation()`.
    pub fn with_surface_curve_tol(c: &dyn Curve, s: &dyn Surface, tol: f64) -> Self {
        let mut this = Self::new();
        this.tolerance = tol.max(PAPPROXIMATION);
        this.perform(c, s);
        this
    }

    /// `Perform(C, S)` (`cxx:1171-1493`).
    pub fn perform(&mut self, c: &dyn Curve, s: &dyn Surface) {
        let c_type = curve_type(c);
        let s_type = surface_type(s);

        let surf_is_anal = is_ana_surf(s);
        let curv_is_anal = c_type != CurveType::BSplineCurve
            && c_type != CurveType::BezierCurve
            && c_type != CurveType::OffsetCurve
            && c_type != CurveType::OtherCurve;

        let mut simplecase = surf_is_anal && curv_is_anal;
        if c_type == CurveType::BSplineCurve || c_type == CurveType::BezierCurve {
            let mut a_nb_knots = 1;
            if c_type == CurveType::BSplineCurve {
                a_nb_knots = c
                    .as_bspline_curve()
                    .map(|b| b.distinct_knots().len() as i32)
                    .unwrap_or(1);
            }
            simplecase =
                simplecase && c.nurbs_degree().unwrap_or(0) as i32 <= 2 && a_nb_knots <= 2;
        }

        if c_type == CurveType::BSplineCurve && s_type == SurfaceType::Plane {
            // `cxx:1197-1227`: project the poles, keep knots and weights.
            if let Some(bs) = c.as_bspline_curve() {
                let proj = PlaneProjector::new(&s.gp_pln().expect("GeomAbs_Plane").pos);
                let poles2d: Vec<GpPnt2d> = bs.poles.iter().map(|p| proj.project(p)).collect();
                self.bspline = build_bspline(
                    &poles2d,
                    bs.weights.clone(),
                    bs.knots.clone(),
                    bs.degree,
                    bs.periodic(),
                )
                .map(Arc::new);
            }
        } else if c_type == CurveType::BezierCurve && s_type == SurfaceType::Plane {
            // `cxx:1228-1252`.
            if let Some(p3d) = c.bezier_poles() {
                let proj = PlaneProjector::new(&s.gp_pln().expect("GeomAbs_Plane").pos);
                let poles2d: Vec<GpPnt2d> = p3d.iter().map(|p| proj.project(p)).collect();
                self.bezier = match c.bezier_weights() {
                    Some(w) => Geom2dBezierCurve::rational(poles2d, w.to_vec()),
                    None => Geom2dBezierCurve::new(poles2d),
                }
                .ok()
                .map(Arc::new);
            }
        } else {
            self.fit_and_divide(c, s, s_type, simplecase);
        }
    }

    /// `cxx:1253-1492`: the `FitAndDivide2d` arm.
    fn fit_and_divide(&mut self, c: &dyn Curve, s: &dyn Surface, s_type: SurfaceType, simplecase: bool) {
        let f = ProjLibFunction::new(c, s);

        // `cxx:1287-1305`.
        let mut deg1 = 5;
        let deg2 = if simplecase { 8 } else { 10 };
        if self.deg_min > 0 {
            deg1 = self.deg_min;
        }
        let deg2 = if self.deg_max > 0 { self.deg_max } else { deg2 };

        let mut a_max_segments = 1000;
        if self.max_segments > 0 {
            a_max_segments = self.max_segments;
        }

        let mut a_fist_c = Constraint::TangencyPoint;
        let mut a_last_c = Constraint::TangencyPoint;
        if self.bnd_pnt != ProjConstraint::TangencyPoint {
            let mapped = match self.bnd_pnt {
                ProjConstraint::NoConstraint => Constraint::NoConstraint,
                ProjConstraint::PassPoint => Constraint::PassPoint,
                ProjConstraint::TangencyPoint => Constraint::TangencyPoint,
            };
            a_fist_c = mapped;
            a_last_c = mapped;
        }

        // `cxx:1318-1321`.
        let a_tol_u = compute_tol_u(s, self.tolerance);
        let a_tol_v = compute_tol_v(s, self.tolerance);
        let a_tol2d = (a_tol_u * a_tol_u + a_tol_v * a_tol_v).sqrt().max(PCONFUSION);

        // `cxx:1323-1329`.
        let mut fit = FitAndDivide2d::from_settings(
            deg1,
            deg2,
            self.tolerance,
            a_tol2d,
            true,
            a_fist_c,
            a_last_c,
        );
        fit.set_max_segments(a_max_segments);
        if simplecase {
            fit.set_hang_checking(false);
        }
        fit.perform(&f);

        let mut a_new_tol2d = 0.0f64;
        if fit.is_all_approximated() {
            // `cxx:1333-1397`.
            let nb_curves = fit.nb_multi_curves() as i32;
            let mut conv = CompBezierCurves2d::new(1.0e-4);
            for i in 1..=nb_curves {
                let (_, tol2d) = fit.error(i);
                a_new_tol2d = a_new_tol2d.max(tol2d);
                let mc = fit.value(i);
                let poles2d = mc.curve2d(1);
                conv.add_curve(&poles2d);
            }
            conv.perform();
            let nb_poles = conv.nb_poles();
            let nb_knots = conv.nb_knots();
            if nb_poles == 0 || nb_poles > 100000 {
                return;
            }
            if nb_knots == 0 || nb_knots > 100000 {
                return;
            }

            let (knots, mults) = conv.knots_and_mults();
            let mults = mults.to_vec();
            let mut new_knots = knots.to_vec();
            let degree = conv.degree() as usize;
            let new_poles: Vec<GpPnt2d> = conv.poles().to_vec();

            bspl_knots::reparametrize(c.first_parameter(), c.last_parameter(), &mut new_knots);
            // Set NewKnots(NbKnots) exactly C->LastParameter() to avoid
            // problems if trim is used (`cxx:1376-1377`).
            let last = new_knots.len() - 1;
            new_knots[last] = c.last_parameter();

            let flat = banded_interp::knot_sequence(&new_knots, &mults, degree as i32);
            let mut built = build_bspline(&new_poles, None, flat, degree, false);

            if a_fist_c == Constraint::PassPoint || a_last_c == Constraint::PassPoint {
                // `cxx:1385-1397`: try to smooth the curve to `GeomAbs_C1`.
                // `OK` short-circuits the loop in OCCT, so a failed removal
                // stops the remaining attempts.
                if let Some(bs) = built.as_mut() {
                    let a_deg = bs.degree() as i32;
                    let a_smooth_tol = CONFUSION.max(a_new_tol2d);
                    let mut ok = true;
                    for ij in 2..nb_knots as i32 {
                        ok = ok && bs.remove_knot(ij, a_deg - 1, a_smooth_tol).unwrap_or(false);
                    }
                }
            }
            self.bspline = built.map(Arc::new);
        } else {
            // `cxx:1398-1408`.
            let nb_curves = fit.nb_multi_curves() as i32;
            if nb_curves != 0 {
                let (_, tol2d) = fit.error(nb_curves);
                a_new_tol2d = tol2d;
            }
        }

        // restore tolerance 3d from 2d (`cxx:1417`).
        self.tolerance *= a_new_tol2d / a_tol2d;

        // Return curve home (`cxx:1419-1491`).
        let u_first = f.first_parameter();
        let u_last = f.last_parameter();
        let u_mid = (u_first + u_last) / 2.0;
        let p3d = c.d0(u_mid);
        let (u, mut v) = match s_type {
            SurfaceType::Plane => slib::plane_parameters(&s.gp_pln().expect("Plane").pos, &p3d),
            SurfaceType::Cylinder => {
                let cy = s.gp_cylinder().expect("Cylinder");
                slib::cylinder_parameters(&cy.pos, &p3d)
            }
            SurfaceType::Cone => {
                let co = s.gp_cone().expect("Cone");
                slib::cone_parameters(&co.pos, co.radius, co.semi_angle, &p3d)
            }
            SurfaceType::Sphere => slib::sphere_parameters(&s.gp_sphere().expect("Sphere").pos, &p3d),
            SurfaceType::Torus => {
                let to = s.gp_torus().expect("Torus");
                slib::torus_parameters(&to.pos, to.major_radius, to.minor_radius, &p3d)
            }
            _ => panic!("Standard_NoSuchObject: ProjLib_ComputeApprox::Value"),
        };

        let b = f.bounds();
        let mut to_mirror = false;
        let mut du = 0.0f64;
        let mut dv = 0.0f64;
        if b.v_couture {
            if s_type == SurfaceType::Sphere && (u - b.u1).abs() > PI {
                to_mirror = true;
                dv = -PI;
                v = PI - v;
            }
            let new_v = clib::in_period(v, b.v1, b.v2);
            let number = ((new_v - v) / (b.v2 - b.v1)).floor() as i32;
            dv -= number as f64 * (b.v2 - b.v1);
        }
        if b.u_couture || (b.v_couture && s_type == SurfaceType::Sphere) {
            let p2d = f.value_2d(u_mid);
            du = u - p2d.x();
            du = if du < 0.0 { du - PCONFUSION } else { du + PCONFUSION };
            let number = (du / PI).trunc() as i32;
            du = number as f64 * PI;
        }

        if let Some(bs) = self.bspline.as_ref() {
            let mut moved = (**bs).clone();
            if du != 0.0 || dv != 0.0 {
                let mut trsf = GpTrsf2d::identity();
                trsf.set_translation_vec(&GpVec2d::new(du, dv));
                moved.transform(&trsf);
            }
            if to_mirror {
                let mut trsf = GpTrsf2d::identity();
                trsf.set_mirror_ax2d(&GpAx2d::new(
                    GpPnt2d::new(0.0, 0.0),
                    GpDir2d::new(1.0, 0.0).expect("gp::DX2d"),
                ));
                moved.transform(&trsf);
            }
            self.bspline = Some(Arc::new(moved));
        }
    }

    /// `SetTolerance(theTolerance)` (`cxx:1496-1499`).
    pub fn set_tolerance(&mut self, the_tolerance: f64) {
        self.tolerance = the_tolerance;
    }

    /// `SetDegree(theDegMin, theDegMax)` (`cxx:1503-1507`).
    pub fn set_degree(&mut self, the_deg_min: i32, the_deg_max: i32) {
        self.deg_min = the_deg_min;
        self.deg_max = the_deg_max;
    }

    /// `SetMaxSegments(theMaxSegments)` (`cxx:1511-1514`).
    pub fn set_max_segments(&mut self, the_max_segments: i32) {
        self.max_segments = the_max_segments;
    }

    /// `SetBndPnt(theBndPnt)` (`cxx:1518-1521`).
    pub fn set_bnd_pnt(&mut self, the_bnd_pnt: ProjConstraint) {
        self.bnd_pnt = the_bnd_pnt;
    }

    /// `BSpline()` (`cxx:1526-1529`).
    pub fn bspline(&self) -> Option<Arc<Geom2dBSplineCurve>> {
        self.bspline.clone()
    }

    /// `Bezier()` (`cxx:1534-1537`).
    pub fn bezier(&self) -> Option<Arc<Geom2dBezierCurve>> {
        self.bezier.clone()
    }

    /// `Tolerance()` (`cxx:1542-1545`).
    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }
}

impl Default for ComputeApprox {
    fn default() -> Self {
        Self::new()
    }
}

/// `new Geom2d_BSplineCurve(Poles[, Weights], Knots, Mults, Degree, IsPeriodic)`:
/// the poles are split into the two coordinate arrays the port's 2d B-spline
/// stores, and the flat knot vector is handed over as is.
fn build_bspline(
    poles: &[GpPnt2d],
    weights: Option<Vec<f64>>,
    knots: Vec<f64>,
    degree: usize,
    periodic: bool,
) -> Option<Geom2dBSplineCurve> {
    let xs: Vec<f64> = poles.iter().map(|p| p.x()).collect();
    let ys: Vec<f64> = poles.iter().map(|p| p.y()).collect();
    Geom2dBSplineCurve::from_flat(xs, ys, weights, knots, degree, periodic).ok()
}
