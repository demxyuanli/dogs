//! `IntSurf_Quadric`. Source: `IntSurf_Quadric.hxx`, `.lxx`, `.cxx`.
//!
//! Implicit representation of the five elementary quadrics (plane, cylinder,
//! cone, sphere, torus) used by the `IntSurf`/`IntCurveSurface`/`IntPatch`
//! analytic intersection arms.
//!
//! OCCT overloads the constructor and `SetValue` on the five `gp` quadrics;
//! Rust cannot, so each overload is spelled out (`from_plane`, `set_value_plane`
//! ... `from_torus`, `set_value_torus`).

use occt_core::elib::{clib, slib, surface_eval};
use occt_core::gp::{
    GpAx3, GpCone, GpCylinder, GpDir, GpLin, GpPln, GpPnt, GpSphere, GpTorus, GpVec,
};
use occt_core::kernel::geomabs::SurfaceType;
use occt_core::precision::REAL_SMALL;

/// `hasMagnitudeForNormalization` (`IntSurf_Quadric.cxx:30-34`):
/// `SquareMagnitude() > gp::Resolution()^2`.
fn has_magnitude_for_normalization(v: &GpVec) -> bool {
    v.square_magnitude() > REAL_SMALL * REAL_SMALL
}

/// `IntSurf_Quadric` (`IntSurf_Quadric.hxx:40-110`).
#[derive(Debug, Clone)]
pub struct IntSurfQuadric {
    ax3: GpAx3,
    lin: GpLin,
    typ: SurfaceType,
    prm1: f64,
    prm2: f64,
    prm3: f64,
    prm4: f64,
    ax3direc: bool,
}

impl Default for IntSurfQuadric {
    /// `IntSurf_Quadric()` (`cxx:37-46`).
    fn default() -> Self {
        Self {
            ax3: GpAx3::default(),
            lin: GpLin::default(),
            typ: SurfaceType::OtherSurface,
            prm1: 0.0,
            prm2: 0.0,
            prm3: 0.0,
            prm4: 0.0,
            ax3direc: false,
        }
    }
}

impl IntSurfQuadric {
    /// `IntSurf_Quadric()` (`cxx:37-46`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `IntSurf_Quadric(const gp_Pln&)` (`cxx:49-54`).
    pub fn from_plane(p: &GpPln) -> Self {
        let mut r = Self {
            ax3: p.position(),
            typ: SurfaceType::Plane,
            ..Self::default()
        };
        r.ax3direc = r.ax3.is_direct();
        let (a, b, c, d) = p.coefficients();
        r.prm1 = a;
        r.prm2 = b;
        r.prm3 = c;
        r.prm4 = d;
        r
    }

    /// `IntSurf_Quadric(const gp_Cylinder&)` (`cxx:57-68`).
    pub fn from_cylinder(c: &GpCylinder) -> Self {
        let mut r = Self {
            ax3: c.position().clone(),
            typ: SurfaceType::Cylinder,
            ..Self::default()
        };
        r.lin = GpLin::new(r.ax3.axis().clone());
        r.prm2 = 0.0;
        r.prm3 = 0.0;
        r.prm4 = 0.0;
        r.ax3direc = r.ax3.is_direct();
        r.prm1 = c.radius;
        r
    }

    /// `IntSurf_Quadric(const gp_Sphere&)` (`cxx:71-81`).
    pub fn from_sphere(s: &GpSphere) -> Self {
        let mut r = Self {
            ax3: s.position().clone(),
            typ: SurfaceType::Sphere,
            ..Self::default()
        };
        r.lin = GpLin::new(r.ax3.axis().clone());
        r.prm2 = 0.0;
        r.prm3 = 0.0;
        r.prm4 = 0.0;
        r.ax3direc = r.ax3.is_direct();
        r.prm1 = s.radius;
        r
    }

    /// `IntSurf_Quadric(const gp_Cone&)` (`cxx:84-96`).
    pub fn from_cone(c: &GpCone) -> Self {
        let mut r = Self {
            ax3: c.position().clone(),
            typ: SurfaceType::Cone,
            ..Self::default()
        };
        r.ax3direc = r.ax3.is_direct();
        r.lin.set_position(r.ax3.axis().clone());
        r.prm1 = c.radius;
        r.prm2 = c.semi_angle;
        r.prm3 = r.prm2.cos();
        r.prm4 = 0.0;
        r
    }

    /// `IntSurf_Quadric(const gp_Torus&)` (`cxx:99-111`).
    pub fn from_torus(t: &GpTorus) -> Self {
        let mut r = Self {
            ax3: t.position().clone(),
            typ: SurfaceType::Torus,
            ..Self::default()
        };
        r.ax3direc = r.ax3.is_direct();
        r.lin.set_position(r.ax3.axis().clone());
        r.prm1 = t.major_radius;
        r.prm2 = t.minor_radius;
        r.prm3 = 0.0;
        r.prm4 = 0.0;
        r
    }

    /// `SetValue(const gp_Pln&)` (`cxx:115-120`).
    pub fn set_value_plane(&mut self, p: &GpPln) {
        self.typ = SurfaceType::Plane;
        self.ax3 = p.position();
        self.ax3direc = self.ax3.is_direct();
        let (a, b, c, d) = p.coefficients();
        self.prm1 = a;
        self.prm2 = b;
        self.prm3 = c;
        self.prm4 = d;
    }

    /// `SetValue(const gp_Cylinder&)` (`cxx:123-130`).
    pub fn set_value_cylinder(&mut self, c: &GpCylinder) {
        self.typ = SurfaceType::Cylinder;
        self.ax3 = c.position().clone();
        self.ax3direc = self.ax3.is_direct();
        self.lin.set_position(self.ax3.axis().clone());
        self.prm1 = c.radius;
        self.prm2 = 0.0;
        self.prm3 = 0.0;
        self.prm4 = 0.0;
    }

    /// `SetValue(const gp_Sphere&)` (`cxx:133-141`).
    pub fn set_value_sphere(&mut self, s: &GpSphere) {
        self.typ = SurfaceType::Sphere;
        self.ax3 = s.position().clone();
        self.ax3direc = self.ax3.is_direct();
        self.lin.set_position(self.ax3.axis().clone());
        self.prm1 = s.radius;
        self.prm2 = 0.0;
        self.prm3 = 0.0;
        self.prm4 = 0.0;
    }

    /// `SetValue(const gp_Cone&)` (`cxx:144-153`).
    pub fn set_value_cone(&mut self, c: &GpCone) {
        self.typ = SurfaceType::Cone;
        self.ax3 = c.position().clone();
        self.ax3direc = self.ax3.is_direct();
        self.lin.set_position(self.ax3.axis().clone());
        self.prm1 = c.radius;
        self.prm2 = c.semi_angle;
        self.prm3 = self.prm2.cos();
        self.prm4 = 0.0;
    }

    /// `SetValue(const gp_Torus&)` (`cxx:156-167`).
    pub fn set_value_torus(&mut self, t: &GpTorus) {
        self.typ = SurfaceType::Torus;
        self.ax3 = t.position().clone();
        self.ax3direc = self.ax3.is_direct();
        self.lin.set_position(self.ax3.axis().clone());
        self.prm1 = t.major_radius;
        self.prm2 = t.minor_radius;
        self.prm3 = 0.0;
        self.prm4 = 0.0;
    }

    /// `TypeQuadric()` (`lxx:21-24`).
    pub fn type_quadric(&self) -> SurfaceType {
        self.typ
    }

    /// `Plane()` (`lxx:28-31`).
    pub fn plane(&self) -> GpPln {
        GpPln::new(self.ax3.clone())
    }

    /// `Sphere()` (`lxx:35-38`).
    pub fn sphere(&self) -> GpSphere {
        GpSphere::new(self.ax3.clone(), self.prm1)
            .expect("IntSurf_Quadric::Sphere: negative radius")
    }

    /// `Cylinder()` (`lxx:42-45`).
    pub fn cylinder(&self) -> GpCylinder {
        GpCylinder::new(self.ax3.clone(), self.prm1)
            .expect("IntSurf_Quadric::Cylinder: negative radius")
    }

    /// `Cone()` (`lxx:48-51`).
    pub fn cone(&self) -> GpCone {
        GpCone::new(self.ax3.clone(), self.prm1, self.prm2)
            .expect("IntSurf_Quadric::Cone: negative radius")
    }

    /// `Torus()` (`lxx:53-56`).
    pub fn torus(&self) -> GpTorus {
        GpTorus::new(self.ax3.clone(), self.prm1, self.prm2)
            .expect("IntSurf_Quadric::Torus: negative radius")
    }

    /// `Distance(const gp_Pnt&)` (`cxx:172-210`).
    pub fn distance(&self, p: &GpPnt) -> f64 {
        match self.typ {
            SurfaceType::Plane => self.prm1 * p.x() + self.prm2 * p.y() + self.prm3 * p.z() + self.prm4,
            SurfaceType::Cylinder => self.lin.distance(p) - self.prm1,
            SurfaceType::Sphere => self.lin.pos.loc.distance(p) - self.prm1,
            SurfaceType::Cone => {
                let dist = self.lin.distance(p);
                let (u, v) = slib::cone_parameters(&self.ax3, self.prm1, self.prm2, p);
                let pp = slib::cone_value(&self.cone(), u, v);
                let distp = self.lin.distance(&pp);
                (dist - distp) / self.prm3
            }
            SurfaceType::Torus => {
                let o = self.ax3.location();
                let oz = self.ax3.direction();
                let oz_v = GpVec::from_xyz(oz.xyz());
                let dot = GpVec::from_pnts(&o, p).dot(&oz_v);
                let pp = p.translated_vec(&oz_v.multiplied_scalar(-dot));
                let dop = if o.square_distance(&pp) < 1e-14 {
                    self.ax3.x_direction().clone()
                } else {
                    GpDir::from_vec(&GpVec::from_pnts(&o, &pp)).expect("torus: degenerate axis")
                };
                let pt = GpPnt::from_xyz(&o.coord.added(&dop.xyz().multiplied(self.prm1)));
                p.distance(&pt) - self.prm2
            }
            _ => 0.0,
        }
    }

    /// `Gradient(const gp_Pnt&)` (`cxx:213-303`).
    pub fn gradient(&self, p: &GpPnt) -> GpVec {
        match self.typ {
            SurfaceType::Plane => GpVec::new(self.prm1, self.prm2, self.prm3),
            SurfaceType::Cylinder => {
                let param = clib::line_parameter(&self.lin, p);
                let pp = self
                    .lin
                    .pos
                    .loc
                    .coord
                    .added(&self.lin.pos.vdir.xyz().multiplied(param));
                let mut grad = GpVec::from_xyz(&p.coord.subtracted(&pp));
                let n = grad.magnitude();
                if n > 1e-14 {
                    grad = grad.divide(n);
                } else {
                    grad = GpVec::new(0.0, 0.0, 0.0);
                }
                grad
            }
            SurfaceType::Sphere => {
                let mut grad = GpVec::from_xyz(&p.coord.subtracted(&self.lin.pos.loc.coord));
                let n = grad.magnitude();
                if n > 1e-14 {
                    grad = grad.divide(n);
                } else {
                    grad = GpVec::new(0.0, 0.0, 0.0);
                }
                grad
            }
            SurfaceType::Cone => {
                let (u, v) = slib::cone_parameters(&self.ax3, self.prm1, self.prm2, p);
                let (_, d1u, d1v) = surface_eval::cone_d1(&self.cone(), u, v);
                let mut grad = d1u.crossed(&d1v);
                if !self.ax3direc {
                    grad = grad.reversed();
                }
                if has_magnitude_for_normalization(&grad) {
                    grad = grad.normalized();
                } else {
                    grad = GpVec::new(0.0, 0.0, 0.0);
                }
                grad
            }
            SurfaceType::Torus => {
                let o = self.ax3.location();
                let oz = self.ax3.direction();
                let oz_v = GpVec::from_xyz(oz.xyz());
                let dot = GpVec::from_pnts(&o, p).dot(&oz_v);
                let pp = p.translated_vec(&oz_v.multiplied_scalar(-dot));
                let dop = if o.square_distance(&pp) < 1e-14 {
                    self.ax3.x_direction().clone()
                } else {
                    GpDir::from_vec(&GpVec::from_pnts(&o, &pp)).expect("torus: degenerate axis")
                };
                let pt = GpPnt::from_xyz(&o.coord.added(&dop.xyz().multiplied(self.prm1)));
                let mut grad = GpVec::from_xyz(&p.coord.subtracted(&pt.coord));
                let n = grad.magnitude();
                if n > 1e-14 {
                    grad = grad.divide(n);
                } else {
                    grad = GpVec::new(0.0, 0.0, 0.0);
                }
                grad
            }
            _ => GpVec::new(0.0, 0.0, 0.0),
        }
    }

    /// `ValAndGrad(const gp_Pnt&, double&, gp_Vec&)` (`cxx:306-404`).
    pub fn val_and_grad(&self, p: &GpPnt) -> (f64, GpVec) {
        match self.typ {
            SurfaceType::Plane => (
                self.prm1 * p.x() + self.prm2 * p.y() + self.prm3 * p.z() + self.prm4,
                GpVec::new(self.prm1, self.prm2, self.prm3),
            ),
            SurfaceType::Cylinder => {
                let dist = self.lin.distance(p) - self.prm1;
                let param = clib::line_parameter(&self.lin, p);
                let pp = self
                    .lin
                    .pos
                    .loc
                    .coord
                    .added(&self.lin.pos.vdir.xyz().multiplied(param));
                let mut grad = GpVec::from_xyz(&p.coord.subtracted(&pp));
                let n = grad.magnitude();
                if n > 1e-14 {
                    grad = grad.divide(n);
                } else {
                    grad = GpVec::new(0.0, 0.0, 0.0);
                }
                (dist, grad)
            }
            SurfaceType::Sphere => {
                let dist = self.lin.pos.loc.distance(p) - self.prm1;
                let mut grad = GpVec::from_xyz(&p.coord.subtracted(&self.lin.pos.loc.coord));
                let n = grad.magnitude();
                if n > 1e-14 {
                    grad = grad.divide(n);
                } else {
                    grad = GpVec::new(0.0, 0.0, 0.0);
                }
                (dist, grad)
            }
            SurfaceType::Cone => {
                let dist = self.lin.distance(p);
                let (u, v) = slib::cone_parameters(&self.ax3, self.prm1, self.prm2, p);
                let (pp, d1u, d1v) = surface_eval::cone_d1(&self.cone(), u, v);
                let distp = self.lin.distance(&pp);
                let dist = (dist - distp) / self.prm3;
                let mut grad = d1u.crossed(&d1v);
                if !self.ax3direc {
                    grad = grad.reversed();
                }
                if has_magnitude_for_normalization(&grad) {
                    grad = grad.normalized();
                } else {
                    grad = GpVec::new(0.0, 0.0, 0.0);
                }
                (dist, grad)
            }
            SurfaceType::Torus => {
                let o = self.ax3.location();
                let oz = self.ax3.direction();
                let oz_v = GpVec::from_xyz(oz.xyz());
                let dot = GpVec::from_pnts(&o, p).dot(&oz_v);
                let pp = p.translated_vec(&oz_v.multiplied_scalar(-dot));
                let dop = if o.square_distance(&pp) < 1e-14 {
                    self.ax3.x_direction().clone()
                } else {
                    GpDir::from_vec(&GpVec::from_pnts(&o, &pp)).expect("torus: degenerate axis")
                };
                let pt = GpPnt::from_xyz(&o.coord.added(&dop.xyz().multiplied(self.prm1)));
                let dist = p.distance(&pt) - self.prm2;
                let mut grad = GpVec::from_xyz(&p.coord.subtracted(&pt.coord));
                let n = grad.magnitude();
                if n > 1e-14 {
                    grad = grad.divide(n);
                } else {
                    grad = GpVec::new(0.0, 0.0, 0.0);
                }
                (dist, grad)
            }
            _ => (0.0, GpVec::new(0.0, 0.0, 0.0)),
        }
    }

    /// `Value(double, double)` (`cxx:407-428`).
    pub fn value(&self, u: f64, v: f64) -> GpPnt {
        match self.typ {
            SurfaceType::Plane => slib::plane_value(&self.plane(), u, v),
            SurfaceType::Cylinder => slib::cylinder_value(&self.cylinder(), u, v),
            SurfaceType::Sphere => slib::sphere_value(&self.sphere(), u, v),
            SurfaceType::Cone => slib::cone_value(&self.cone(), u, v),
            SurfaceType::Torus => slib::torus_value(&self.torus(), u, v),
            _ => GpPnt::new(0.0, 0.0, 0.0),
        }
    }

    /// `D1(double, double, gp_Pnt&, gp_Vec&, gp_Vec&)` (`cxx:432-456`).
    pub fn d1(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
        match self.typ {
            SurfaceType::Plane => surface_eval::plane_d1(&self.plane(), u, v),
            SurfaceType::Cylinder => surface_eval::cylinder_d1(&self.cylinder(), u, v),
            SurfaceType::Sphere => surface_eval::sphere_d1(&self.sphere(), u, v),
            SurfaceType::Cone => surface_eval::cone_d1(&self.cone(), u, v),
            SurfaceType::Torus => surface_eval::torus_d1(&self.torus(), u, v),
            _ => (
                GpPnt::new(0.0, 0.0, 0.0),
                GpVec::new(0.0, 0.0, 0.0),
                GpVec::new(0.0, 0.0, 0.0),
            ),
        }
    }

    /// `DN(double, double, int, int)` (`cxx:459-481`).
    pub fn dn(&self, u: f64, v: f64, nu: i32, nv: i32) -> GpVec {
        match self.typ {
            SurfaceType::Plane => surface_eval::plane_dn(u, v, &self.ax3, nu, nv),
            SurfaceType::Cylinder => {
                surface_eval::cylinder_dn(u, v, &self.ax3, self.prm1, nu, nv)
            }
            SurfaceType::Sphere => surface_eval::sphere_dn(u, v, &self.ax3, self.prm1, nu, nv),
            SurfaceType::Cone => {
                surface_eval::cone_dn(u, v, &self.ax3, self.prm1, self.prm2, nu, nv)
            }
            SurfaceType::Torus => {
                surface_eval::torus_dn(u, v, &self.ax3, self.prm1, self.prm2, nu, nv)
            }
            _ => GpVec::new(0.0, 0.0, 0.0),
        }
    }

    /// `Normale(double, double)` (`cxx:484-521`).
    pub fn normale_uv(&self, u: f64, v: f64) -> GpVec {
        match self.typ {
            SurfaceType::Plane => {
                let d = self.ax3.direction();
                if self.ax3direc {
                    GpVec::from_xyz(d.xyz())
                } else {
                    GpVec::from_xyz(d.xyz()).reversed()
                }
            }
            SurfaceType::Cylinder => self.normale(&self.value(u, v)),
            SurfaceType::Sphere => self.normale(&self.value(u, v)),
            SurfaceType::Cone => {
                let (_, d1u, d1v) = surface_eval::cone_d1(&self.cone(), u, v);
                if d1u.magnitude() < 0.0000001 {
                    return GpVec::new(0.0, 0.0, 0.0);
                }
                d1u.crossed(&d1v)
            }
            SurfaceType::Torus => self.normale(&self.value(u, v)),
            _ => GpVec::new(0.0, 0.0, 0.0),
        }
    }

    /// `Normale(const gp_Pnt&)` (`cxx:524-587`).
    pub fn normale(&self, p: &GpPnt) -> GpVec {
        match self.typ {
            SurfaceType::Plane => {
                let d = self.ax3.direction();
                if self.ax3direc {
                    GpVec::from_xyz(d.xyz())
                } else {
                    GpVec::from_xyz(d.xyz()).reversed()
                }
            }
            SurfaceType::Cylinder => {
                let d = self.lin.normal(p).pos.vdir;
                if self.ax3direc {
                    GpVec::from_xyz(d.xyz())
                } else {
                    GpVec::from_xyz(d.xyz()).reversed()
                }
            }
            SurfaceType::Sphere => {
                let loc = self.ax3.location();
                if self.ax3direc {
                    let v = GpVec::from_pnts(&loc, p);
                    GpVec::from_xyz(v.xyz())
                } else {
                    let v = GpVec::from_pnts(p, &loc);
                    GpVec::from_xyz(v.xyz())
                }
            }
            SurfaceType::Cone => {
                let (u, v) = slib::cone_parameters(&self.ax3, self.prm1, self.prm2, p);
                self.normale_uv(u, v)
            }
            SurfaceType::Torus => {
                let o = self.ax3.location();
                let oz = self.ax3.direction();
                let oz_v = GpVec::from_xyz(oz.xyz());
                let dot = GpVec::from_pnts(&o, p).dot(&oz_v);
                let pp = p.translated_vec(&oz_v.multiplied_scalar(-dot));
                let dop = if o.square_distance(&pp) < 1e-14 {
                    self.ax3.x_direction().clone()
                } else {
                    GpDir::from_vec(&GpVec::from_pnts(&o, &pp)).expect("torus: degenerate axis")
                };
                let pt = GpPnt::from_xyz(&o.coord.added(&dop.xyz().multiplied(self.prm1)));
                if pt.square_distance(p) < 1e-14 {
                    return GpVec::from_xyz(oz_v.xyz());
                }
                if self.ax3direc {
                    GpVec::from_pnts(&pt, p)
                } else {
                    GpVec::from_pnts(p, &pt)
                }
            }
            _ => GpVec::new(0.0, 0.0, 0.0),
        }
    }

    /// `Parameters(const gp_Pnt&, double&, double&)` (`cxx:590-612`).
    pub fn parameters(&self, p: &GpPnt) -> (f64, f64) {
        match self.typ {
            SurfaceType::Plane => slib::plane_parameters(&self.ax3, p),
            SurfaceType::Cylinder => slib::cylinder_parameters(&self.ax3, p),
            SurfaceType::Sphere => slib::sphere_parameters(&self.ax3, p),
            SurfaceType::Cone => slib::cone_parameters(&self.ax3, self.prm1, self.prm2, p),
            SurfaceType::Torus => slib::torus_parameters(&self.ax3, self.prm1, self.prm2, p),
            _ => (0.0, 0.0),
        }
    }
}
