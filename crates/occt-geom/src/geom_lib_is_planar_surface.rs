//! `GeomLib_IsPlanarSurface` (`GeomLib_IsPlanarSurface.cxx:1-260`,
//! `GeomLib_IsPlanarSurface.hxx:20-60`): tests whether a surface is planar
//! within a tolerance and, when it is, returns the plane that supports it.
//!
//! Control flow follows the `.cxx` by surface type:
//!
//! * `GeomAbs_Plane` answers `Adaptor3d_Surface::Plane()` directly.
//! * Cylinder / cone / sphere / torus are never planar.
//! * A surface of revolution is tested by the angle between its axis and the
//!   normal at the mid-parameter, then by `Controle` on its `UIso(Umin)`.
//! * A surface of extrusion is tested by `gp_Dir::IsNormal(AS.Direction(),
//!   Tol / (Vmax - Vmin))`, then by `Controle` on its `VIso` mid-value.
//! * Anything else is sampled on `NbU * NbV` points and passed to the
//!   inertia-based `Controle`.
//!
//! Deviations: `StdFail_NotDone` out of `Plan()` answers `Err` (the port has
//! no exception), and `GeomAdaptor_Surface` is replaced by the `Surface`
//! trait accessors that carry the same elementary geometries
//! (`projlib::projected_curve::surface_type` uses the same dispatch order).

use std::sync::Arc;

use occt_core::gp::{GpAx3, GpDir, GpPln, GpPnt};
use occt_core::kernel::geomabs::{CurveType, Shape, SurfaceType};
use occt_core::precision::{CONFUSION, RESOLUTION};

use crate::curve::Curve;
use crate::geom_lib::geom_lib_inertia;
use crate::projlib::{curve_type, surface_type};
use crate::surface::Surface;

/// `GeomLib_IsPlanarSurface` (`hxx:20-60`): the `IsPlan` flag and `myPlan`.
pub struct IsPlanarSurface {
    is_plan: bool,
    plan: GpPln,
}

impl IsPlanarSurface {
    /// `GeomLib_IsPlanarSurface::GeomLib_IsPlanarSurface(S, Tol)`
    /// (`cxx:155-260`), the whole `switch (AS.GetType())`.
    pub fn new(s: &Arc<dyn Surface>, tol: f64) -> Self {
        let mut this = Self {
            is_plan: false,
            plan: GpPln::new(GpAx3::standard()),
        };
        match surface_type(s.as_ref()) {
            SurfaceType::Plane => {
                // `cxx:163-168`: trust the elementary plane.
                this.is_plan = true;
                if let Some(pln) = s.gp_pln() {
                    this.plan = pln;
                }
            }
            SurfaceType::Cylinder
            | SurfaceType::Cone
            | SurfaceType::Sphere
            | SurfaceType::Torus => {
                // `cxx:169-174`.
                this.is_plan = false;
            }
            SurfaceType::SurfaceOfRevolution => {
                this.revolution_case(s, tol);
            }
            SurfaceType::SurfaceOfExtrusion => {
                this.extrusion_case(s, tol);
            }
            _ => {
                this.sampled_case(s, tol);
            }
        }
        this
    }

    /// `cxx:176-234`: the `GeomAbs_SurfaceOfRevolution` arm.
    fn revolution_case(&mut self, s: &Arc<dyn Surface>, tol: f64) {
        let mut essai = true;
        let (umin, umax) = s.u_range();
        let (vmin, vmax) = s.v_range();
        let (_, mut du, mut dv) = s.d1((umin + umax) / 2.0, (vmin + vmax) / 2.0);
        if du.magnitude() <= RESOLUTION || dv.magnitude() <= RESOLUTION {
            let new_u = (umin + umax) / 2.0 + (umax - umin) * 0.1;
            let new_v = (vmin + vmax) / 2.0 + (vmax - vmin) * 0.1;
            let (_, du2, dv2) = s.d1(new_u, new_v);
            du = du2;
            dv = dv2;
        }
        let dn = du.crossed(&dv);
        if dn.magnitude() > 1.0e-7 {
            let Some(axis) = s.revolution_axis() else {
                self.is_plan = false;
                return;
            };
            let mut dir = axis.direction().clone();
            let mut angle = dir.angle(&match GpDir::from_vec(&dn) {
                Ok(d) => d,
                Err(_) => {
                    self.is_plan = false;
                    return;
                }
            });
            if angle > std::f64::consts::PI / 2.0 {
                angle = std::f64::consts::PI - angle;
                dir.reverse();
            }
            essai = angle < 0.1;
            if essai {
                // `cxx:203-215`: the plane through the mid point, `Dir` as
                // normal, `DU` as X direction.
                let p = s.d0((umin + umax) / 2.0, (vmin + vmax) / 2.0);
                let xdir = match GpDir::from_vec(&du) {
                    Ok(d) => d,
                    Err(_) => {
                        self.is_plan = false;
                        return;
                    }
                };
                if let Ok(ax3) = GpAx3::new(p.clone(), dir, &xdir) {
                    self.plan = GpPln::new(ax3);
                    self.plan.set_location(&p);
                    if let Some(c) = s.u_iso_curve(umin) {
                        self.is_plan = controle_curve(&c, &self.plan, tol);
                    } else {
                        self.is_plan = false;
                    }
                } else {
                    self.is_plan = false;
                }
            } else {
                self.is_plan = false;
            }
        } else {
            // `cxx:196-202`: `Essai` keeps its `true` initial value, so the
            // plane is still built and checked.
            if essai {
                let p = s.d0((umin + umax) / 2.0, (vmin + vmax) / 2.0);
                let dir = match GpDir::from_vec(&dn) {
                    Ok(d) => d,
                    Err(_) => {
                        self.is_plan = false;
                        return;
                    }
                };
                let xdir = match GpDir::from_vec(&du) {
                    Ok(d) => d,
                    Err(_) => {
                        self.is_plan = false;
                        return;
                    }
                };
                if let Ok(ax3) = GpAx3::new(p.clone(), dir, &xdir) {
                    self.plan = GpPln::new(ax3);
                    self.plan.set_location(&p);
                    if let Some(c) = s.u_iso_curve(umin) {
                        self.is_plan = controle_curve(&c, &self.plan, tol);
                    } else {
                        self.is_plan = false;
                    }
                } else {
                    self.is_plan = false;
                }
            }
        }
    }

    /// `cxx:235-281`: the `GeomAbs_SurfaceOfExtrusion` arm.
    fn extrusion_case(&mut self, s: &Arc<dyn Surface>, tol: f64) {
        let (umin, umax) = s.u_range();
        let (vmin, vmax) = s.v_range();
        let (p, mut du, mut dv) = s.d1((umin + umax) / 2.0, (vmin + vmax) / 2.0);
        if du.magnitude() <= RESOLUTION || dv.magnitude() <= RESOLUTION {
            let new_u = (umin + umax) / 2.0 + (umax - umin) * 0.1;
            let new_v = (vmin + vmax) / 2.0 + (vmax - vmin) * 0.1;
            let (_, du2, dv2) = s.d1(new_u, new_v);
            du = du2;
            dv = dv2;
        }
        let mut dn = du.crossed(&dv);
        let norm = dn.magnitude();
        let mut essai = false;
        if norm > 1.0e-15 {
            dn = dn.divided(norm);
            let angmax = tol / (vmax - vmin);
            if let (Ok(d), Some(dir)) = (GpDir::from_vec(&dn), s.extrusion_direction()) {
                essai = d.is_normal_tol(&dir, angmax);
            }
        }
        if essai {
            let zdir = match GpDir::from_vec(&dn) {
                Ok(d) => d,
                Err(_) => {
                    self.is_plan = false;
                    return;
                }
            };
            let xdir = match GpDir::from_vec(&du) {
                Ok(d) => d,
                Err(_) => {
                    self.is_plan = false;
                    return;
                }
            };
            if let Ok(ax3) = GpAx3::new(p.clone(), zdir, &xdir) {
                self.plan = GpPln::new(ax3);
                self.plan.set_location(&p);
                if let Some(c) = s.v_iso_curve((vmin + vmax) / 2.0) {
                    self.is_plan = controle_curve(&c, &self.plan, tol);
                } else {
                    self.is_plan = false;
                }
            } else {
                self.is_plan = false;
            }
        } else {
            self.is_plan = false;
        }
    }

    /// `cxx:282-315`: the `default` arm, `NbU * NbV` samples.
    fn sampled_case(&mut self, s: &Arc<dyn Surface>, tol: f64) {
        let nbu = 8 + 3 * s.nb_u_intervals(Shape::CN as u8);
        let nbv = 8 + 3 * s.nb_v_intervals(Shape::CN as u8);
        let (umin, umax) = s.u_range();
        let (vmin, vmax) = s.v_range();
        let du = (umax - umin) / (nbu as f64 - 1.0);
        let dv = (vmax - vmin) / (nbv as f64 - 1.0);
        let mut pnts = Vec::with_capacity((nbu * nbv) as usize);
        for ii in 0..nbu {
            let u = umin + du * ii as f64;
            for jj in 0..nbv {
                let v = vmin + dv * jj as f64;
                pnts.push(s.d0(u, v));
            }
        }
        let mut plan = GpPln::new(GpAx3::standard());
        self.is_plan = controle_points(&pnts, tol, s, &mut plan);
        if self.is_plan {
            self.plan = plan;
        }
    }

    /// `GeomLib_IsPlanarSurface::IsPlanar` (`cxx:317-320`).
    pub fn is_planar(&self) -> bool {
        self.is_plan
    }

    /// `GeomLib_IsPlanarSurface::Plan` (`cxx:322-330`): `StdFail_NotDone` when
    /// the surface is not planar.
    pub fn plan(&self) -> Result<GpPln, String> {
        if !self.is_plan {
            return Err(" GeomLib_IsPlanarSurface".to_string());
        }
        Ok(self.plan.clone())
    }
}

/// `Controle(Poles, Tol, S, Plan)` (`cxx:32-98`): the inertia of the point
/// cloud gives the two minor axes, and `Plan` keeps `DX` when it is closer to
/// `DU` than `DY` is.
fn controle_points(
    poles: &[GpPnt],
    tol: f64,
    s: &Arc<dyn Surface>,
    plan: &mut GpPln,
) -> bool {
    let inertia = geom_lib_inertia(poles);
    let (gx, gy, gz) = (inertia.x_gap, inertia.y_gap, inertia.z_gap);
    let _ = gx;
    let mut dx = inertia.x_dir.clone();
    let mut dy = inertia.y_dir.clone();
    if !(gz < tol && gy > CONFUSION) {
        return false;
    }
    let (umin, umax) = s.u_range();
    let (vmin, vmax) = s.v_range();
    let (_, du_v, dv_v) = s.d1((umin + umax) / 2.0, (vmin + vmax) / 2.0);
    if !(du_v.square_magnitude() > RESOLUTION && dv_v.square_magnitude() > RESOLUTION) {
        return false;
    }
    let mut du = match GpDir::from_vec(&du_v) {
        Ok(d) => d,
        Err(_) => return false,
    };
    let mut angle1 = du.angle(&dx);
    let mut angle2 = du.angle(&dy);
    if angle1 > std::f64::consts::PI / 2.0 {
        angle1 = std::f64::consts::PI - angle1;
    }
    if angle2 > std::f64::consts::PI / 2.0 {
        angle2 = std::f64::consts::PI - angle2;
    }
    if angle2 < angle1 {
        // `cxx:76-80`: `du = DY; DY = DX; DX = du;`
        let saved = du;
        du = dy.clone();
        dy = dx.clone();
        dx = saved;
    }
    if dx.angle(&du) > std::f64::consts::PI / 2.0 {
        dx.reverse();
    }
    if let Ok(dv) = GpDir::from_vec(&dv_v) {
        if dy.angle(&dv) > std::f64::consts::PI / 2.0 {
            dy.reverse();
        }
    }
    let normal = match dx.crossed(&dy) {
        Ok(d) => d,
        Err(_) => return false,
    };
    let bary = inertia.bary.clone();
    let Ok(ax3) = GpAx3::new(bary.clone(), normal, &dx) else {
        return false;
    };
    plan.set_position(&ax3);
    plan.set_location(&bary);
    true
}

/// `Controle(C, Plan, Tol)` (`cxx:100-153`): sample count by curve type, then
/// all samples must lie closer than `Tol` to the plane.
fn controle_curve(c: &Arc<dyn Curve>, plan: &GpPln, tol: f64) -> bool {
    let nb = match curve_type(c.as_ref()) {
        CurveType::Line => 2,
        CurveType::Circle => 3,
        CurveType::Ellipse | CurveType::Hyperbola | CurveType::Parabola => 5,
        CurveType::BezierCurve => c.bezier_poles().map_or(16, |p| p.len() as i32),
        CurveType::BSplineCurve => c.bspline_poles().map_or(16, |p| p.len() as i32),
        _ => 8 + 3 * c.nb_intervals(Shape::CN as u8),
    };
    let f = c.first_parameter();
    let l = c.last_parameter();
    let du = (l - f) / (nb as f64 - 1.0);
    let mut b = true;
    for ii in 1..=nb {
        let u = (ii - 1) as f64 * du + f;
        let p = c.value(u);
        let d = plan.signed_distance_pnt(&p).abs();
        b = d < tol;
        if !b {
            break;
        }
    }
    b
}
