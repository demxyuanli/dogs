//! `GeomFill` — the package-level static helpers, split by the OCCT
//! overload family: this file keeps the ruled/special surface through two
//! curves, the section shape / knots / mults of a
//! `Convert_ParameterisationType` and the minimal weights; [`circle`] keeps
//! the three `GetCircle` overloads.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill.cxx:40-340` (OCCT 8.0.0) and
//! `GeomFill.hxx`.
//!
//! Translation notes:
//! - `NCollection_Array1<gp_Pnt>` / `<double>` / `<int>` out-parameters become
//!   slices; `GeomFill::GetShape` sizes them, exactly like OCCT
//!   (`cxx:186-225`).
//! - `Convert_ParameterisationType&` of `GetShape` (`cxx:189`) is
//!   `&mut ParameterisationType` because the default arm rewrites it
//!   (`cxx:214-224`).
//! - `Geom_TrimmedCurve::BasisCurve()` / `FirstParameter()` / `LastParameter()`
//!   (`cxx:60-83`) map to [`Curve::untrimmed_basis`], which reports the basis
//!   parameters of the trim.
//! - `occ::handle<Geom_Surface>` is `Arc<dyn Surface>`; the
//!   `GeomFill_Generator` fallback (`cxx:173-178`) is
//!   [`generator_surface`].

mod circle;

use std::sync::Arc;

use occt_core::convert::ParameterisationType;
use occt_core::gp::{GpAx2, GpAx3, GpCirc, GpCone, GpCylinder, GpDir, GpPln, GpVec};
use occt_core::precision::{ANGULAR, CONFUSION, PCONFUSION};

use crate::circle::GeomCircle;
use crate::cone::GeomCone;
use crate::convert_bspl::curve_to_bspline_curve;
use crate::curve::Curve;
use crate::cylinder::GeomCylinder;
use crate::geom_fill::generator_surface;
use crate::plane::GeomPlane;
use crate::rectangular_trimmed::GeomRectangularTrimmedSurface;
use crate::surface::Surface;
use crate::trimmed::GeomTrimmedCurveBasis;

/// `gp_Vec` view of a `gp_Dir`.
fn dir_vec(d: &GpDir) -> GpVec {
    GpVec::new(d.x(), d.y(), d.z())
}

/// Unit circle in the XY plane, `gp_Circ(gp_Ax2(gp_Pnt(0,0,0), gp_Dir(0,0,1)), 1)`
/// (`cxx:242-244`).
fn unit_circle() -> Arc<dyn Curve> {
    let circ = GpCirc::new(GpAx2::standard(), 1.0);
    Arc::new(GeomCircle::new(circ))
}

/// `GeomFill` (`GeomFill.hxx`), the static utility class.
pub struct GeomFill;

impl GeomFill {
    /// `GeomFill::Surface(Curve1, Curve2)` (`cxx:40-181`).
    ///
    /// A plane / cylinder / cone when the two curves are two trimmed lines or
    /// two circles, otherwise the `GeomFill_Generator` ruled surface.
    pub fn surface(
        curve1: &Arc<dyn Curve>,
        curve2: &Arc<dyn Curve>,
    ) -> Result<Arc<dyn Surface>, String> {
        let mut a1 = 0.0;
        let mut a2 = 0.0;
        let mut b1 = 0.0;
        let mut b2 = 0.0;
        let mut trim1 = false;
        let mut trim2 = false;

        // `cxx:60-83`: a trimmed curve is unwrapped to its basis and its trim
        // parameters; any other curve is copied as-is.
        let the_curve1: Arc<dyn Curve> = if curve1.is_geom_trimmed() {
            let (basis, first, last) = curve1
                .untrimmed_basis()
                .ok_or_else(|| "GeomFill::Surface: trimmed curve without basis".to_string())?;
            a1 = first;
            b1 = last;
            trim1 = true;
            basis
        } else {
            Arc::from(curve1.clone_dyn())
        };
        let the_curve2: Arc<dyn Curve> = if curve2.is_geom_trimmed() {
            let (basis, first, last) = curve2
                .untrimmed_basis()
                .ok_or_else(|| "GeomFill::Surface: trimmed curve without basis".to_string())?;
            a2 = first;
            b2 = last;
            trim2 = true;
            basis
        } else {
            Arc::from(curve2.clone_dyn())
        };

        let mut surf: Option<Arc<dyn Surface>> = None;

        let line1 = the_curve1.gp_line();
        let line2 = the_curve2.gp_line();
        if let (Some(l1), Some(l2)) = (line1, line2) {
            if trim1 && trim2 {
                let d1 = l1.direction();
                let d2 = l2.direction();
                if d1.is_parallel_tol(&d2, ANGULAR) {
                    let p1p2 = GpVec::from_pnts(&l1.location(), &l2.location());
                    let proj = p1p2.dot(&dir_vec(&d1));
                    if d1.is_equal_tol(&d2, ANGULAR) {
                        if (a1 - proj - a2).abs() <= CONFUSION
                            && (b1 - proj - b2).abs() <= CONFUSION
                        {
                            let normal = GpDir::from_vec(&dir_vec(&d1).crossed(&p1p2))?;
                            let ax = GpAx3::new(l1.location(), normal, &d1)?;
                            let v = p1p2.dot(&dir_vec(ax.y_direction()));
                            let plane: Arc<dyn Surface> = Arc::new(GeomPlane::new(GpPln::new(ax)));
                            surf = Some(Arc::new(GeomRectangularTrimmedSurface::uv(
                                plane,
                                a1,
                                b1,
                                0.0f64.min(v),
                                0.0f64.max(v),
                            )?));
                        }
                    }
                    if d1.is_opposite_ang(&d2, ANGULAR) {
                        if (a1 - proj + b2).abs() <= CONFUSION
                            && (b1 - proj + a2).abs() <= CONFUSION
                        {
                            let normal = GpDir::from_vec(&dir_vec(&d1).crossed(&p1p2))?;
                            let ax = GpAx3::new(l1.location(), normal, &d1)?;
                            let v = p1p2.dot(&dir_vec(ax.y_direction()));
                            let plane: Arc<dyn Surface> = Arc::new(GeomPlane::new(GpPln::new(ax)));
                            surf = Some(Arc::new(GeomRectangularTrimmedSurface::uv(
                                plane,
                                a1,
                                b1,
                                0.0f64.min(v),
                                0.0f64.max(v),
                            )?));
                        }
                    }
                }
            }
        } else if let (Some(c1), Some(c2)) = (the_curve1.gp_circ(), the_curve2.gp_circ()) {
            // `cxx:128-171`: two coaxial circles give a cylinder or a cone.
            let a1_ax = c1.position();
            let a2_ax = c2.position();
            if a1_ax
                .axis()
                .is_coaxial_tol(&a2_ax.axis(), ANGULAR, CONFUSION)
            {
                let mut v = GpVec::from_pnts(&a1_ax.location(), &a2_ax.location())
                    .dot(&dir_vec(&a1_ax.direction()));
                if !trim1 && !trim2 {
                    let mut pos = a1_ax.to_ax3();
                    if (c1.radius() - c2.radius()).abs() < CONFUSION {
                        let cyl = GeomCylinder::new(GpCylinder::new(pos, c1.radius())?);
                        let cyl: Arc<dyn Surface> = Arc::new(cyl);
                        surf = Some(Arc::new(GeomRectangularTrimmedSurface::one_param(
                            cyl,
                            0.0f64.min(v),
                            0.0f64.max(v),
                            false,
                        )?));
                    } else {
                        let rad = c2.radius() - c1.radius();
                        let mut ang = (rad / v).atan();
                        if ang < 0.0 {
                            pos.z_reverse();
                            v = -v;
                            ang = -ang;
                        }
                        let cone = GeomCone::new(GpCone::new(pos, c1.radius(), ang)?);
                        v /= ang.cos();
                        let cone: Arc<dyn Surface> = Arc::new(cone);
                        surf = Some(Arc::new(GeomRectangularTrimmedSurface::one_param(
                            cone,
                            0.0f64.min(v),
                            0.0f64.max(v),
                            false,
                        )?));
                    }
                }
            }
        }

        // `cxx:173-178`: no special surface, fall back to the generator.
        match surf {
            Some(s) => Ok(s),
            None => generator_surface(&[curve1.clone(), curve2.clone()]),
        }
    }

    /// `GeomFill::GetShape(MaxAng, NbPoles, NbKnots, Degree, TConv)`
    /// (`cxx:186-225`): the section shape of a parameterisation.
    pub fn get_shape(
        max_ang: f64,
        t_conv: &mut ParameterisationType,
    ) -> (i32, i32, i32) {
        match *t_conv {
            ParameterisationType::QuasiAngular => (7, 2, 6),
            ParameterisationType::Polynomial => (8, 2, 7),
            _ => {
                let nb_span = (3.0 * max_ang.abs() / 2.0 / std::f64::consts::PI).ceil() as i32;
                let nb_poles = 2 * nb_span + 1;
                let nb_knots = nb_span + 1;
                match nb_span {
                    1 => *t_conv = ParameterisationType::TgtThetaOver2_1,
                    2 => *t_conv = ParameterisationType::TgtThetaOver2_2,
                    3 => *t_conv = ParameterisationType::TgtThetaOver2_3,
                    _ => {}
                }
                (nb_poles, nb_knots, 2)
            }
        }
    }

    /// `GeomFill::GetMinimalWeights(TConv, MinAng, MaxAng, Weights)`
    /// (`cxx:235-264`).
    ///
    /// `Weights` must already have the length `GetShape(MaxAng, ...)` gives.
    pub fn get_minimal_weights(
        t_conv: ParameterisationType,
        min_ang: f64,
        max_ang: f64,
        weights: &mut [f64],
    ) -> Result<(), String> {
        if t_conv == ParameterisationType::Polynomial {
            weights.fill(1.0);
            return Ok(());
        }
        let circ = unit_circle();
        let sect1 = GeomTrimmedCurveBasis::new(circ.clone(), 0.0, max_ang);
        let c1 = curve_to_bspline_curve(&sect1, t_conv).map_err(|e| format!("{e:?}"))?;
        if c1.nb_poles() != weights.len() {
            return Err("Standard_OutOfRange: GeomFill::GetMinimalWeights".into());
        }
        for (i, w) in weights.iter_mut().enumerate() {
            *w = c1.weight(i as i32 + 1);
        }
        let angle_min = PCONFUSION.max(min_ang);
        let sect2 = GeomTrimmedCurveBasis::new(circ, 0.0, angle_min);
        let c2 = curve_to_bspline_curve(&sect2, t_conv).map_err(|e| format!("{e:?}"))?;
        if c2.nb_poles() != weights.len() {
            return Err("Standard_OutOfRange: GeomFill::GetMinimalWeights".into());
        }
        for (i, w) in weights.iter_mut().enumerate() {
            let p = c2.weight(i as i32 + 1);
            if p < *w {
                *w = p;
            }
        }
        Ok(())
    }

    /// `GeomFill::Knots(TConv, TKnots)` (`cxx:268-285`).
    pub fn knots(t_conv: ParameterisationType, t_knots: &mut [f64]) {
        if t_conv != ParameterisationType::QuasiAngular
            && t_conv != ParameterisationType::Polynomial
        {
            for (i, k) in t_knots.iter_mut().enumerate() {
                *k = i as f64;
            }
        } else {
            if t_knots.len() < 2 {
                return;
            }
            t_knots[0] = 0.0;
            t_knots[1] = 1.0;
        }
    }

    /// `GeomFill::Mults(TConv, TMults)` (`cxx:288-315`).
    pub fn mults(t_conv: ParameterisationType, t_mults: &mut [i32]) {
        match t_conv {
            ParameterisationType::QuasiAngular => {
                if t_mults.len() >= 2 {
                    t_mults[0] = 7;
                    t_mults[1] = 7;
                }
            }
            ParameterisationType::Polynomial => {
                if t_mults.len() >= 2 {
                    t_mults[0] = 8;
                    t_mults[1] = 8;
                }
            }
            _ => {
                if t_mults.is_empty() {
                    return;
                }
                let last = t_mults.len() - 1;
                t_mults[0] = 3;
                for m in t_mults.iter_mut().take(last).skip(1) {
                    *m = 2;
                }
                t_mults[last] = 3;
            }
        }
    }

    /// `GeomFill::GetTolerance(TConv, AngleMin, Radius, AngularTol, SpatialTol)`
    /// (`cxx:325-340`).
    pub fn get_tolerance(
        t_conv: ParameterisationType,
        angle_min: f64,
        radius: f64,
        angular_tol: f64,
        spatial_tol: f64,
    ) -> Result<f64, String> {
        let circ = GpCirc::new(GpAx2::standard(), radius);
        let circle = Arc::new(GeomCircle::new(circ)) as Arc<dyn Curve>;
        // `0.02` is about one degree below which the tangency is not
        // considered, so the approximation tolerance does not go to zero.
        let sect = GeomTrimmedCurveBasis::new(circle, 0.0, angle_min.max(0.02));
        let c = curve_to_bspline_curve(&sect, t_conv).map_err(|e| format!("{e:?}"))?;
        let dist = c.pole(0).distance(c.pole(1)) + spatial_tol;
        Ok(dist * angular_tol / 2.0)
    }
}
