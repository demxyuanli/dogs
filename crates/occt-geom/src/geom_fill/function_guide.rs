//! `GeomFill_FunctionGuide` — the residual `G(w) - S(teta, v)` of the
//! "guide curve against surface of revolution" solve.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_FunctionGuide.cxx:36-272` and
//! `GeomFill_FunctionGuide.hxx:33-133` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `Adaptor3d_Curve` / `Geom_Curve` are [`Curve`]; the revolution surface is
//!   [`GeomSurfaceOfRevolution`] (`Geom_SurfaceOfRevolution`).
//! - `TheCurve` is built by `SetParam` from the section law: a trimmed copy of
//!   `TheConst` when the law is constant (`cxx:92-96`), otherwise the law
//!   section at `TheUonS` as a `Geom_BSplineCurve` (`cxx:97-116`). The trim is
//!   [`GeomTrimmedCurveBasis`], the parameter-preserving `Geom_TrimmedCurve`
//!   view (`cxx:212-215` evaluates in basis parameters), so `TheSurface`'s `V`
//!   range is `First..Last` exactly as in OCCT.
//! - `SetParam` raises through `gp_Ax3` when `XDir` is not perpendicular to
//!   `Dir`, so it returns `Result<(), String>` here.
//! - `Deriv2T` (`hxx:93-96`) is declared in OCCT but has no definition: the
//!   whole block `cxx:280-364` is commented out. It stays UNPORTED.
//! - `DSDT` (`cxx:231-272`) leaves the `Q` derivative at zero for a
//!   non-constant law; OCCT prints `"Not implemented"` (`cxx:248-251`).

use std::sync::Arc;

use occt_core::gp::{GpAx1, GpAx3, GpDir, GpPnt, GpTrsf, GpVec, GpXyz};
use occt_core::precision::CONFUSION;
use occt_math::function_set_root::MathFunctionSetWithDerivatives;
use occt_math::{MathMatrix, MathVector};

use crate::bspline_curve::GeomBSplineCurve;
use crate::curve::Curve;
use crate::geom_fill::section_law::SectionLaw;
use crate::surface::Surface;
use crate::surface_of_revolution::GeomSurfaceOfRevolution;
use crate::trimmed::GeomTrimmedCurveBasis;

/// `gp_Vec::Coord(i)` with OCCT's 1-based index.
fn comp(v: &GpVec, i: usize) -> f64 {
    match i {
        1 => v.coord.x,
        2 => v.coord.y,
        _ => v.coord.z,
    }
}

/// `gp_Pnt::Coord(i)` with OCCT's 1-based index.
fn comp_p(p: &GpPnt, i: usize) -> f64 {
    match i {
        1 => p.coord.x,
        2 => p.coord.y,
        _ => p.coord.z,
    }
}

/// `GeomFill_FunctionGuide` (`GeomFill_FunctionGuide.hxx:33-133`).
pub struct FunctionGuide {
    /// `TheGuide` (`hxx`): the guide curve.
    guide: Arc<dyn Curve>,
    /// `TheLaw` (`hxx`): the section law.
    law: Arc<dyn SectionLaw>,
    /// `isconst` (`hxx`): `TheLaw->IsConstant(Confusion)` (`cxx:54-65`).
    isconst: bool,
    /// `TheCurve` (`hxx`): the transformed section curve.
    the_curve: Option<Arc<dyn Curve>>,
    /// `TheConst` (`hxx`): `TheLaw->ConstantSection()` when constant.
    the_const: Option<Arc<dyn Curve>>,
    /// `TheSurface` (`hxx`): revolution of `TheCurve` about `(Centre, Dir)`.
    the_surface: Option<GeomSurfaceOfRevolution>,
    /// `First` (`hxx`).
    first: f64,
    /// `Last` (`hxx`).
    last: f64,
    /// `TheUonS` (`hxx`): parameter on the law.
    the_uon_s: f64,
    /// `Centre` (`hxx`).
    centre: GpXyz,
    /// `Dir` (`hxx`).
    dir: GpXyz,
}

impl FunctionGuide {
    /// `GeomFill_FunctionGuide(S, C, Param)` (`cxx:43-69`).
    pub fn new(law: Arc<dyn SectionLaw>, guide: Arc<dyn Curve>, param: f64) -> Self {
        let mut tol = CONFUSION;
        let mut me = Self {
            guide,
            law,
            isconst: false,
            the_curve: None,
            the_const: None,
            the_surface: None,
            first: 0.0,
            last: 0.0,
            the_uon_s: param,
            centre: GpXyz::zero(),
            dir: GpXyz::zero(),
        };
        if me.law.is_constant(&mut tol) {
            me.isconst = true;
            // `IsConstant` succeeded, so `ConstantSection` is available.
            let the_const = me.law.constant_section().ok();
            if let Some(c) = &the_const {
                me.first = c.first_parameter();
                me.last = c.last_parameter();
            }
            me.the_const = the_const;
        } else {
            me.isconst = false;
            me.the_const = None;
        }
        me
    }

    /// `SetParam(Param, C, D, DX)` (`cxx:76-121`): rebuild `TheCurve` and its
    /// revolution surface for the frame `(C, D, DX)`.
    pub fn set_param(&mut self, c: &GpPnt, d: &GpXyz, dx: &GpXyz) -> Result<(), String> {
        self.centre = c.coord;
        self.dir = *d;

        // Fixed frame `gp::XOY()`.
        let rep = GpAx3::new(
            GpPnt::zero(),
            GpDir::new(0.0, 0.0, 1.0)?,
            &GpDir::new(1.0, 0.0, 0.0)?,
        )?;

        // Transformation between the trihedron and Oxyz (`cxx:86-91`).
        let z = GpDir::from_xyz(d)?;
        let x = GpDir::from_xyz(dx)?;
        let rep_triedre = GpAx3::new(*c, z, &x)?;
        let mut transfo = GpTrsf::identity();
        transfo.set_transformation_from_to(&rep_triedre, &rep);

        let mut curve: Box<dyn Curve> = if self.isconst {
            let base = self
                .the_const
                .as_ref()
                .ok_or_else(|| "GeomFill_FunctionGuide: TheConst null".to_string())?;
            // `Geom_TrimmedCurve(TheConst->Copy(), First, Last)` (`cxx:94`).
            let copy: Arc<dyn Curve> = Arc::from(base.as_ref().clone_dyn());
            Box::new(GeomTrimmedCurveBasis::new(copy, self.first, self.last))
        } else {
            let (nb_poles, nb_knots, deg) = self.law.section_shape();
            let mults = self.law.mults();
            let knots = self.law.knots();
            let (poles, weights) = self.law.d0(self.the_uon_s)?;
            if nb_poles as usize != poles.len() || nb_knots as usize != knots.len() {
                return Err("GeomFill_FunctionGuide: inconsistent section shape".into());
            }
            let flat = occt_core::bspl::banded_interp::knot_sequence(&knots, &mults, deg);
            let periodic = self.law.is_u_periodic();
            let bs = if self.law.is_rational() {
                GeomBSplineCurve::from_parts(poles, Some(weights), flat, deg.max(0) as usize, periodic)
            } else {
                GeomBSplineCurve::from_parts(poles, None, flat, deg.max(0) as usize, periodic)
            };
            Box::new(bs)
        };

        let axe = GpAx1::new(*c, GpDir::from_xyz(d)?);
        curve.transform(&transfo);
        let generatrix: Arc<dyn Curve> = Arc::from(curve);
        self.the_curve = Some(generatrix.clone());
        self.the_surface = Some(GeomSurfaceOfRevolution::new(generatrix, axe));
        Ok(())
    }

    /// `NbVariables` (`cxx:128-131`): `(w, u, v)`.
    pub fn nb_variables(&self) -> usize {
        3
    }

    /// `NbEquations` (`cxx:137-140`).
    pub fn nb_equations(&self) -> usize {
        3
    }

    /// `Value(X, F)` (`cxx:145-155`): `G(X(1)) - S(X(2), X(3))`.
    pub fn value(&self, x: &MathVector, f: &mut MathVector) -> bool {
        let Some(surface) = &self.the_surface else {
            return false;
        };
        let p = self.guide.d0(x.value(1));
        let p1 = surface.d0(x.value(2), x.value(3));
        for i in 1..=3 {
            f.set_value(i, comp_p(&p, i) - comp_p(&p1, i));
        }
        true
    }

    /// `Derivatives(X, D)` (`cxx:162-179`).
    pub fn derivatives(&self, x: &MathVector, d: &mut MathMatrix) -> bool {
        let Some(surface) = &self.the_surface else {
            return false;
        };
        let (_, dp) = self.guide.d1(x.value(1));
        let (_, dp1u, dp1v) = surface.d1(x.value(2), x.value(3));
        for i in 1..=3 {
            d.set_value(i, 1, comp(&dp, i));
            d.set_value(i, 2, -comp(&dp1u, i));
            d.set_value(i, 3, -comp(&dp1v, i));
        }
        true
    }

    /// `Values(X, F, D)` (`cxx:186-203`).
    pub fn values(&self, x: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        let Some(surface) = &self.the_surface else {
            return false;
        };
        let (p, dp) = self.guide.d1(x.value(1));
        let (p1, dp1u, dp1v) = surface.d1(x.value(2), x.value(3));
        for i in 1..=3 {
            f.set_value(i, comp_p(&p, i) - comp_p(&p1, i));
            d.set_value(i, 1, comp(&dp, i));
            d.set_value(i, 2, -comp(&dp1u, i));
            d.set_value(i, 3, -comp(&dp1v, i));
        }
        true
    }

    /// `DerivT(X, DCentre, DDir, F)` (`cxx:208-226`).
    pub fn deriv_t(&self, x: &MathVector, dcentre: &GpXyz, ddir: &GpXyz, f: &mut MathVector) -> bool {
        let Ok(ds) = self.dsdt(x.value(2), x.value(3), dcentre, ddir) else {
            return false;
        };
        let Some(curve) = &self.the_curve else {
            return false;
        };
        let p = curve.d0(x.value(1));
        for i in 1..=3 {
            f.set_value(i, comp_p(&p, i) - comp(&ds, i));
        }
        true
    }

    /// `DSDT(U, V, DC, DDir, DS)` (`cxx:231-272`): derivative of the
    /// revolution surface with respect to the revolution parameter `t`.
    fn dsdt(&self, u: f64, v: f64, dc: &GpXyz, ddir: &GpXyz) -> Result<GpVec, String> {
        let Some(curve) = &self.the_curve else {
            return Err("GeomFill_FunctionGuide: TheCurve null".into());
        };
        // `Q(v)` (`cxx:244`).
        let mut q = curve.d0(v).coord;
        let mut dq = GpXyz::zero();
        if !self.isconst {
            // `cxx:248-251`: OCCT prints "Not implemented"; `DQ` stays zero.
        }

        q = q.subtract(&self.centre); // CQ
        dq = dq.subtract(dc);

        // `DDir x CQ + Dir x DQ`, scaled by sin(U) (`cxx:254-259`).
        let dv_cross_cq = ddir.crossed(&q).add(&self.dir.crossed(&dq)).multiplied(u.sin());

        // `(DDir.Q + Dir.DQ) * Dir + (Dir.Q) * DDir` (`cxx:261-266`).
        let mut dv_dot_cq = GpXyz::zero();
        let s1 = ddir.dot(&q) + self.dir.dot(&dq);
        let s2 = self.dir.dot(&q);
        dv_dot_cq.set_linear_form(s1, &self.dir, s2, ddir);
        dv_dot_cq = dv_dot_cq.add(&dv_cross_cq);

        dq = dq.multiplied(u.cos()).add(&dv_dot_cq).add(dc);
        Ok(GpVec::from_xyz(&dq))
    }
}

impl MathFunctionSetWithDerivatives for FunctionGuide {
    fn nb_variables(&self) -> usize {
        FunctionGuide::nb_variables(self)
    }

    fn nb_equations(&self) -> usize {
        FunctionGuide::nb_equations(self)
    }

    fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        FunctionGuide::value(self, x, f)
    }

    fn derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        FunctionGuide::derivatives(self, x, d)
    }

    fn values(&mut self, x: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        FunctionGuide::values(self, x, f, d)
    }
}
