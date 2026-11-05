//! `Extrema_GlobOptFuncConicS` — the squared distance between a conic and a
//! surface, parameterized by the two surface parameters `(su, sv)`.
//! Source: `Extrema_GlobOptFuncConicS.cxx:22-239`, `.hxx:24-74`.
//!
//! ```text
//! F(su,sv) = min over the conic of |S(su,sv) - C(t)|^2
//! ```
//!
//! For every `(su, sv)` the point `S(su,sv)` is projected onto the conic's
//! analytic parameter (`ElCLib::Parameter`), then the distance is taken against
//! that projection, against the trimmed end points `myCPf` / `myCPl`, and the
//! smallest wins (`cxx:52-64`). The reduced search dimension is what makes
//! `GlobMinConicS` (`Extrema_GenExtCS.cxx:523-708`) a 2-variable PSO.
//!
//! `ConicParameter` (`cxx:188-238`) repeats the same minimization but returns
//! the *parameter*, with the same tie-breaking order as `value`.

use occt_core::precision::INFINITE;
use occt_math::{MathMultipleVarFunction, MathVector};

use super::prelude::*;

/// `Precision::RealLast()` — the initial `F` of `value` / `ConicParameter`.
const REAL_LAST: f64 = f64::MAX;

/// `GeomAbs_CurveType` restricted to the conics `LoadConic` stores
/// (`cxx:143-161`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ConicType {
    /// `GeomAbs_Line`.
    Line,
    /// `GeomAbs_Circle`.
    Circle,
    /// `GeomAbs_Ellipse`.
    Ellipse,
    /// `GeomAbs_Hyperbola`.
    Hyperbola,
    /// `GeomAbs_Parabola`.
    Parabola,
    /// The `default:` arm of both switches (`cxx:158-160`, `:207-210`).
    Other,
}

/// `Extrema_GlobOptFuncConicS` (`Extrema_GlobOptFuncConicS.hxx:24-74`).
pub struct GlobOptFuncConicS<'a> {
    /// `myS`.
    s: &'a dyn Surface,
    /// `myUf` / `myUl` / `myVf` / `myVl`.
    u_f: f64,
    u_l: f64,
    v_f: f64,
    v_l: f64,
    /// `myC` — set by [`GlobOptFuncConicS::load_conic`] only.
    c: Option<&'a dyn Curve>,
    /// `myTf` / `myTl`.
    t_f: f64,
    t_l: f64,
    /// `myCPf` / `myCPl`.
    c_pf: GpPnt,
    c_pl: GpPnt,
    /// `myCType`.
    c_type: ConicType,
    /// `myLin`.
    lin: Option<GpLin>,
    /// `myCirc`.
    circ: Option<GpCirc>,
    /// `myElips`.
    elips: Option<GpElips>,
    /// `myHypr`.
    hypr: Option<GpHypr>,
    /// `myParab`.
    parab: Option<GpParab>,
}

/// `GpCirc` is not part of the shared prelude; the circle arm needs it.
use occt_core::gp::GpCirc;

impl<'a> GlobOptFuncConicS<'a> {
    /// `Extrema_GlobOptFuncConicS(S, theUf, theUl, theVf, theVl)`
    /// (`cxx:78-89`).
    pub fn new(s: &'a dyn Surface, u_f: f64, u_l: f64, v_f: f64, v_l: f64) -> Self {
        GlobOptFuncConicS {
            s,
            u_f,
            u_l,
            v_f,
            v_l,
            c: None,
            t_f: 0.0,
            t_l: 0.0,
            c_pf: GpPnt::new(0.0, 0.0, 0.0),
            c_pl: GpPnt::new(0.0, 0.0, 0.0),
            c_type: ConicType::Other,
            lin: None,
            circ: None,
            elips: None,
            hypr: None,
            parab: None,
        }
    }

    /// `Extrema_GlobOptFuncConicS(S)` (`cxx:93-101`): the surface's own range.
    pub fn with_surface(s: &'a dyn Surface) -> Self {
        let (u_f, u_l) = s.u_range();
        let (v_f, v_l) = s.v_range();
        GlobOptFuncConicS::new(s, u_f, u_l, v_f, v_l)
    }

    /// `Extrema_GlobOptFuncConicS(C, S)` (`cxx:105-117`).
    pub fn with_conic(c: &'a dyn Curve, s: &'a dyn Surface) -> Self {
        let (u_f, u_l) = s.u_range();
        let (v_f, v_l) = s.v_range();
        let mut this = GlobOptFuncConicS::new(s, u_f, u_l, v_f, v_l);
        let a_ctf = c.first_parameter();
        let a_ctl = c.last_parameter();
        this.load_conic(c, a_ctf, a_ctl);
        this
    }

    /// `checkInputData(X, su, sv)` (`cxx:68-74`): `None` is the `false`
    /// return.
    fn check_input_data(&self, x: &MathVector) -> Option<(f64, f64)> {
        let su = x.value(1);
        let sv = x.value(2);
        if su >= self.u_f && su <= self.u_l && sv >= self.v_f && sv <= self.v_l {
            Some((su, sv))
        } else {
            None
        }
    }

    /// The conic parameter of `p` under `myCType` (`cxx:26-47` and
    /// `cxx:192-211` share this switch). `None` is the `default:` arm.
    fn conic_parameter_of(&self, p: &GpPnt) -> Option<f64> {
        match self.c_type {
            ConicType::Line => self.lin.as_ref().map(|l| clib::line_parameter(l, p)),
            ConicType::Circle => self
                .circ
                .as_ref()
                .map(|c| clib::circle_parameter(&c.pos, p)),
            ConicType::Ellipse => self.elips.as_ref().map(|e| {
                clib::ellipse_parameter(&e.pos, e.major_radius, e.minor_radius, p)
            }),
            ConicType::Hyperbola => self.hypr.as_ref().map(|h| {
                clib::hyperbola_parameter(&h.pos, h.major_radius, h.minor_radius, p)
            }),
            ConicType::Parabola => self.parab.as_ref().map(|pr| clib::parabola_parameter(&pr.pos, p)),
            ConicType::Other => None,
        }
    }

    /// `value(su, sv, F)` (`cxx:22-64`).
    fn value_at(&self, su: f64, sv: f64) -> f64 {
        let a_ps = self.s.d0(su, sv);
        let mut ct = match self.conic_parameter_of(&a_ps) {
            Some(ct) => ct,
            None => return INFINITE,
        };
        if self.c_type == ConicType::Circle || self.c_type == ConicType::Ellipse {
            if self.t_l > 2.0 * std::f64::consts::PI + PCONFUSION {
                ct += 2.0 * std::f64::consts::PI;
            }
        }
        let mut f = REAL_LAST;
        if ct >= self.t_f && ct <= self.t_l {
            if let Some(c) = self.c {
                let a_pc = c.d0(ct);
                f = f.min(a_ps.square_distance(&a_pc));
            }
        }
        f = f.min(a_ps.square_distance(&self.c_pf));
        f = f.min(a_ps.square_distance(&self.c_pl));
        f
    }

    /// `NbVariables()` (`cxx:164-167`): 2.
    pub fn nb_variables(&self) -> usize {
        2
    }

    /// `Value(X, F)` (`cxx:172-184`): `false` when `X` is out of range or `F`
    /// is infinite.
    pub fn value(&self, x: &MathVector, f: &mut f64) -> bool {
        match self.check_input_data(x) {
            Some((su, sv)) => {
                let v = self.value_at(su, sv);
                *f = v;
                !occt_core::precision::Precision::is_infinite(v)
            }
            None => false,
        }
    }

    /// `ConicParameter(theUV)` (`cxx:188-238`).
    pub fn conic_parameter(&self, uv: &MathVector) -> f64 {
        let a_ps = self.s.d0(uv.value(1), uv.value(2));
        let mut ct = match self.conic_parameter_of(&a_ps) {
            Some(ct) => ct,
            // `cxx:207-210`: `ct = myTf; return ct;`
            None => return self.t_f,
        };
        if self.c_type == ConicType::Circle || self.c_type == ConicType::Ellipse {
            if self.t_l > 2.0 * std::f64::consts::PI + PCONFUSION {
                ct += 2.0 * std::f64::consts::PI;
            }
        }
        let mut f = REAL_LAST;
        if ct >= self.t_f && ct <= self.t_l {
            if let Some(c) = self.c {
                let a_pc = c.d0(ct);
                f = f.min(a_ps.square_distance(&a_pc));
            }
        }
        let fext = a_ps.square_distance(&self.c_pf);
        if fext < f {
            f = fext;
            ct = self.t_f;
        }
        let fext = a_ps.square_distance(&self.c_pl);
        if fext < f {
            ct = self.t_l;
        }
        ct
    }

    /// `LoadConic(C, theTf, theTl)` (`cxx:120-162`).
    pub fn load_conic(&mut self, c: &'a dyn Curve, the_tf: f64, the_tl: f64) {
        self.c = Some(c);
        self.t_f = the_tf;
        self.t_l = the_tl;
        if c.is_periodic() {
            let a_t_max = 2.0 * std::f64::consts::PI + PCONFUSION;
            if self.t_f > a_t_max
                || self.t_f < -PCONFUSION
                || (self.t_l - self.t_f).abs() > a_t_max
            {
                let preci = ((self.t_l - self.t_f).abs() / 2.0).min(PCONFUSION);
                let mut u1 = self.t_f;
                let mut u2 = self.t_l;
                clib::adjust_periodic(0.0, 2.0 * std::f64::consts::PI, preci, &mut u1, &mut u2);
                self.t_f = u1;
                self.t_l = u2;
            }
        }
        self.c_pf = c.d0(self.t_f);
        self.c_pl = c.d0(self.t_l);
        if let Some(l) = c.gp_line() {
            self.c_type = ConicType::Line;
            self.lin = Some(l);
        } else if let Some(ci) = c.gp_circ() {
            self.c_type = ConicType::Circle;
            self.circ = Some(ci);
        } else if let Some(e) = c.gp_ellipse() {
            self.c_type = ConicType::Ellipse;
            self.elips = Some(e);
        } else if let Some(h) = c.gp_hyperbola() {
            self.c_type = ConicType::Hyperbola;
            self.hypr = Some(h);
        } else if let Some(pr) = c.gp_parabola() {
            self.c_type = ConicType::Parabola;
            self.parab = Some(pr);
        } else {
            self.c_type = ConicType::Other;
        }
    }
}

impl MathMultipleVarFunction for GlobOptFuncConicS<'_> {
    fn nb_variables(&self) -> usize {
        GlobOptFuncConicS::nb_variables(self)
    }

    fn value(&self, x: &MathVector, f: &mut f64) -> bool {
        GlobOptFuncConicS::value(self, x, f)
    }
}
