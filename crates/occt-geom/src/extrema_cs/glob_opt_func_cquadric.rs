//! `Extrema_GlobOptFuncCQuadric` — the squared distance between a curve and a
//! quadric surface, parameterized by the curve parameter `ct` alone.
//! Source: `Extrema_GlobOptFuncCQuadric.cxx:22-274`, `.hxx:23-73`.
//!
//! For every `ct` the curve point `C(ct)` is projected onto the quadric's
//! analytic parameters (`ElSLib::Parameters`), then the distance is taken
//! against that projection, against the four trimmed corners `myPTrim[0..3]`,
//! and the smallest wins (`cxx:63-73`). `QuadricParameters` (`cxx:207-273`)
//! repeats the same minimization but returns the *parameters*, with the same
//! tie-breaking order as `value`.
//!
//! This is the 1-variable reduction `GlobMinCQuadric`
//! (`Extrema_GenExtCS.cxx:711-823`) runs when the surface is a quadric
//! (`Extrema_GenExtCS.cxx:325-329`).

use occt_core::precision::INFINITE;
use occt_math::{MathMultipleVarFunction, MathVector};

use super::prelude::*;

/// `Precision::RealLast()`.
const REAL_LAST: f64 = f64::MAX;

/// `GeomAbs_SurfaceType` restricted to the quadrics `LoadQuad` stores
/// (`Extrema_GlobOptFuncCQuadric.cxx:161-181`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum QuadricType {
    Plane,
    Cylinder,
    Cone,
    Sphere,
    Torus,
    /// The `default:` arm of both switches (`cxx:178-180`, `:234-238`).
    Other,
}

/// `Extrema_GlobOptFuncCQuadric` (`Extrema_GlobOptFuncCQuadric.hxx:23-73`).
pub struct GlobOptFuncCQuadric<'a> {
    /// `myC`.
    c: &'a dyn Curve,
    /// `myS` — set by [`GlobOptFuncCQuadric::load_quad`].
    s: Option<&'a dyn Surface>,
    /// `myTf` / `myTl`.
    t_f: f64,
    t_l: f64,
    /// `myUf` / `myUl` / `myVf` / `myVl`.
    u_f: f64,
    u_l: f64,
    v_f: f64,
    v_l: f64,
    /// `myPTrim[4]`.
    p_trim: [GpPnt; 4],
    /// `mySType`.
    s_type: QuadricType,
    /// `myPln`.
    pln: Option<GpPln>,
    /// `myCylinder`.
    cylinder: Option<GpCylinder>,
    /// `myCone`.
    cone: Option<GpCone>,
    /// `mySphere`.
    sphere: Option<GpSphere>,
    /// `myTorus`.
    torus: Option<GpTorus>,
}

impl<'a> GlobOptFuncCQuadric<'a> {
    /// `Extrema_GlobOptFuncCQuadric(C, S)` (`cxx:88-98`): the curve range plus
    /// `LoadQuad(S, ...)` over the surface's own range.
    pub fn with_surface(c: &'a dyn Curve, s: &'a dyn Surface) -> Self {
        let mut this = GlobOptFuncCQuadric::with_curve(c);
        let (u_f, u_l) = s.u_range();
        let (v_f, v_l) = s.v_range();
        this.load_quad(s, u_f, u_l, v_f, v_l);
        this
    }

    /// `Extrema_GlobOptFuncCQuadric(C)` (`cxx:102-106`).
    pub fn with_curve(c: &'a dyn Curve) -> Self {
        GlobOptFuncCQuadric {
            c,
            s: None,
            t_f: c.first_parameter(),
            t_l: c.last_parameter(),
            u_f: 0.0,
            u_l: 0.0,
            v_f: 0.0,
            v_l: 0.0,
            p_trim: [
                GpPnt::new(0.0, 0.0, 0.0),
                GpPnt::new(0.0, 0.0, 0.0),
                GpPnt::new(0.0, 0.0, 0.0),
                GpPnt::new(0.0, 0.0, 0.0),
            ],
            s_type: QuadricType::Other,
            pln: None,
            cylinder: None,
            cone: None,
            sphere: None,
            torus: None,
        }
    }

    /// `Extrema_GlobOptFuncCQuadric(C, theTf, theTl)` (`cxx:111-117`).
    pub fn with_range(c: &'a dyn Curve, the_tf: f64, the_tl: f64) -> Self {
        let mut this = GlobOptFuncCQuadric::with_curve(c);
        this.t_f = the_tf;
        this.t_l = the_tl;
        this
    }

    /// `checkInputData(X, ct)` (`cxx:79-84`): `None` is the `false` return.
    fn check_input_data(&self, x: &MathVector) -> Option<f64> {
        let ct = x.value(1);
        if ct >= self.t_f && ct <= self.t_l {
            Some(ct)
        } else {
            None
        }
    }

    /// The `ElSLib::Parameters` switch shared by `value` (`cxx:26-46`) and
    /// `QuadricParameters` (`cxx:218-238`). `None` is the `default:` arm.
    fn parameters_of(&self, p: &GpPnt) -> Option<(f64, f64)> {
        match self.s_type {
            QuadricType::Plane => self.pln.as_ref().map(|pl| slib::plane_parameters(&pl.pos, p)),
            QuadricType::Cylinder => self
                .cylinder
                .as_ref()
                .map(|cy| slib::cylinder_parameters(&cy.pos, p)),
            QuadricType::Cone => self
                .cone
                .as_ref()
                .map(|co| slib::cone_parameters(&co.pos, co.radius, co.semi_angle, p)),
            QuadricType::Sphere => self
                .sphere
                .as_ref()
                .map(|sp| slib::sphere_parameters(&sp.pos, p)),
            QuadricType::Torus => self
                .torus
                .as_ref()
                .map(|to| slib::torus_parameters(&to.pos, to.major_radius, to.minor_radius, p)),
            QuadricType::Other => None,
        }
    }

    /// The `mySType != GeomAbs_Plane` / `== GeomAbs_Torus` period fold of
    /// `cxx:48-62` and `cxx:240-252`.
    fn fold_period(&self, u: &mut f64, v: &mut f64) {
        if self.s_type != QuadricType::Plane
            && self.u_l > 2.0 * std::f64::consts::PI + PCONFUSION
        {
            *u += 2.0 * std::f64::consts::PI;
        }
        if self.s_type == QuadricType::Torus
            && self.v_l > 2.0 * std::f64::consts::PI + PCONFUSION
        {
            *v += 2.0 * std::f64::consts::PI;
        }
    }

    /// `value(ct, F)` (`cxx:22-74`).
    fn value_at(&self, ct: f64) -> f64 {
        let a_cp = self.c.d0(ct);
        let (mut u, mut v) = match self.parameters_of(&a_cp) {
            Some(uv) => uv,
            None => return INFINITE,
        };
        self.fold_period(&mut u, &mut v);

        let mut f = REAL_LAST;
        if u >= self.u_f && u <= self.u_l && v >= self.v_f && v <= self.v_l {
            if let Some(s) = self.s {
                let a_ps = s.d0(u, v);
                f = f.min(a_cp.square_distance(&a_ps));
            }
        }
        for i in 0..4 {
            f = f.min(a_cp.square_distance(&self.p_trim[i]));
        }
        f
    }

    /// `NbVariables()` (`cxx:186-189`): 1.
    pub fn nb_variables(&self) -> usize {
        1
    }

    /// `Value(X, F)` (`cxx:193-203`): `false` when `X` is out of range or `F`
    /// is infinite.
    pub fn value(&self, x: &MathVector, f: &mut f64) -> bool {
        match self.check_input_data(x) {
            Some(ct) => {
                let val = self.value_at(ct);
                *f = val;
                !occt_core::precision::Precision::is_infinite(val)
            }
            None => false,
        }
    }

    /// `QuadricParameters(theCT, theUV)` (`cxx:207-273`).
    pub fn quadric_parameters(&self, ct: f64) -> (f64, f64) {
        // `cxx:212-214`: the corner parameters paired with `myPTrim[0..3]`.
        let uext = [self.u_f, self.u_l, self.u_l, self.u_f];
        let vext = [self.v_f, self.v_f, self.v_l, self.v_l];
        let a_cp = self.c.d0(ct);
        let (mut u, mut v) = match self.parameters_of(&a_cp) {
            Some(uv) => uv,
            // `cxx:234-238`: `theUV(1) = myUf; theUV(2) = myUl; return;`
            None => return (self.u_f, self.u_l),
        };
        self.fold_period(&mut u, &mut v);

        let mut f = REAL_LAST;
        if u >= self.u_f && u <= self.u_l && v >= self.v_f && v <= self.v_l {
            if let Some(s) = self.s {
                let a_ps = s.d0(u, v);
                f = a_cp.square_distance(&a_ps);
            }
        }
        for i in 0..4 {
            let fi = a_cp.square_distance(&self.p_trim[i]);
            if fi < f {
                f = fi;
                u = uext[i];
                v = vext[i];
            }
        }
        (u, v)
    }

    /// `LoadQuad(S, theUf, theUl, theVf, theVl)` (`cxx:123-182`).
    pub fn load_quad(
        &mut self,
        s: &'a dyn Surface,
        the_uf: f64,
        the_ul: f64,
        the_vf: f64,
        the_vl: f64,
    ) {
        self.s = Some(s);
        self.u_f = the_uf;
        self.u_l = the_ul;
        self.v_f = the_vf;
        self.v_l = the_vl;
        let a_t_max = 2.0 * std::f64::consts::PI + PCONFUSION;
        if s.is_u_periodic() {
            if self.u_f > a_t_max
                || self.u_f < -PCONFUSION
                || (self.u_l - self.u_f).abs() > a_t_max
            {
                let preci = ((self.u_l - self.u_f).abs() / 2.0).min(PCONFUSION);
                let mut u1 = self.u_f;
                let mut u2 = self.u_l;
                clib::adjust_periodic(0.0, 2.0 * std::f64::consts::PI, preci, &mut u1, &mut u2);
                self.u_f = u1;
                self.u_l = u2;
            }
        }
        if s.is_v_periodic() {
            if self.v_f > a_t_max
                || self.v_f < -PCONFUSION
                || (self.v_l - self.v_f).abs() > a_t_max
            {
                let preci = ((self.v_l - self.v_f).abs() / 2.0).min(PCONFUSION);
                let mut v1 = self.v_f;
                let mut v2 = self.v_l;
                clib::adjust_periodic(0.0, 2.0 * std::f64::consts::PI, preci, &mut v1, &mut v2);
                self.v_f = v1;
                self.v_l = v2;
            }
        }
        self.p_trim[0] = s.d0(self.u_f, self.v_f);
        self.p_trim[1] = s.d0(self.u_l, self.v_f);
        self.p_trim[2] = s.d0(self.u_l, self.v_l);
        self.p_trim[3] = s.d0(self.u_f, self.v_l);

        if let Some(pl) = s.gp_pln() {
            self.s_type = QuadricType::Plane;
            self.pln = Some(pl);
        } else if let Some(cy) = s.gp_cylinder() {
            self.s_type = QuadricType::Cylinder;
            self.cylinder = Some(cy);
        } else if let Some(co) = s.gp_cone() {
            self.s_type = QuadricType::Cone;
            self.cone = Some(co);
        } else if let Some(sp) = s.gp_sphere() {
            self.s_type = QuadricType::Sphere;
            self.sphere = Some(sp);
        } else if let Some(to) = s.gp_torus() {
            self.s_type = QuadricType::Torus;
            self.torus = Some(to);
        } else {
            self.s_type = QuadricType::Other;
        }
    }
}

impl MathMultipleVarFunction for GlobOptFuncCQuadric<'_> {
    fn nb_variables(&self) -> usize {
        GlobOptFuncCQuadric::nb_variables(self)
    }

    fn value(&self, x: &MathVector, f: &mut f64) -> bool {
        GlobOptFuncCQuadric::value(self, x, f)
    }
}
