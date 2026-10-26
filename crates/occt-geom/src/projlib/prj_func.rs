//! `ProjLib_PrjFunc` - the normal-projection function set.
//!
//! Source: `ModelingData/TKGeomBase/ProjLib/ProjLib_PrjFunc.cxx:25-157` plus
//! `ProjLib_PrjFunc.hxx:22-60` (OCCT 8.0.0).
//!
//! The two equations are `F1 = v . S_u * Norm` and `F2 = v . S_v * Norm`
//! where `v = S(u, v) - C(t)` (`cxx:105-109`). `myFix` decides which of
//! `t` / `u` / `v` is held constant and which two are the variables:
//!
//! - `1`: `t` fixed, variables `(u, v)` - a normal projection at a given `t`,
//! - `2`: `u` fixed, variables `(t, v)` - solve for the `u` iso,
//! - `3`: `v` fixed, variables `(t, u)` - solve for the `v` iso.
//!
//! Translation notes:
//! - `myNorm = min(1., min(UResolution(1.), VResolution(1.)))`
//!   (`cxx:37`) is [`crate::approx_same_parameter::u_resolution`] /
//!   `v_resolution`, i.e. `GeomAdaptor_Surface::UResolution`
//!   (`GeomAdaptor_Surface.cxx:1818-1896`), not `Precision::Parametric`.
//! - `Values` writes the last evaluated `(t, u, v)` back into `myT` / `myU` /
//!   `myV` (`cxx:138-140`); [`PrjFunc::solution`] reads them back.
//! - `Value` / `Derivatives` both go through `Values` (`cxx:67-75`), so the
//!   Jacobian and the residuals are always one evaluation.

use occt_core::gp::GpVec;
use occt_math::{MathFunctionSetWithDerivatives, MathMatrix, MathVector};

use crate::approx_same_parameter::{u_resolution, v_resolution};
use crate::{Curve, Surface};

/// `ProjLib_PrjFunc` (`ProjLib_PrjFunc.hxx:22-60`).
pub struct PrjFunc<'a> {
    curve: &'a dyn Curve,
    surface: &'a dyn Surface,
    /// `myNorm` (`hxx:46`).
    norm: f64,
    /// `myt` (`hxx:47`).
    my_t: f64,
    /// `myU` (`hxx:48`).
    my_u: f64,
    /// `myV` (`hxx:49`).
    my_v: f64,
    /// `myFix` (`hxx:50`).
    fix: i32,
}

impl<'a> PrjFunc<'a> {
    /// `ProjLib_PrjFunc(C, FixVal, S, Fix)` (`cxx:25-54`).
    ///
    /// `myFix` outside `1..=3` raises `Standard_ConstructionError`
    /// (`cxx:49-51`); the port answers `Err` here.
    pub fn new(
        curve: &'a dyn Curve,
        fix_val: f64,
        surface: &'a dyn Surface,
        fix: i32,
    ) -> Result<Self, String> {
        let norm = 1.0_f64.min(u_resolution(surface, 1.0).min(v_resolution(surface, 1.0)));
        let mut this = Self {
            curve,
            surface,
            norm,
            my_t: 0.0,
            my_u: 0.0,
            my_v: 0.0,
            fix,
        };
        match fix {
            1 => this.my_t = fix_val,
            2 => this.my_u = fix_val,
            3 => this.my_v = fix_val,
            _ => return Err("ProjLib_PrjFunc: myFix must be in 1..=3 (cxx:49-51)".into()),
        }
        Ok(this)
    }

    /// `ProjLib_PrjFunc::Solution()` (`cxx:145-157`).
    pub fn solution(&self) -> (f64, f64) {
        match self.fix {
            1 => (self.my_u, self.my_v),
            2 => (self.my_t, self.my_v),
            3 => (self.my_t, self.my_u),
            _ => (0.0, 0.0),
        }
    }

    /// `Values(X, F, D)` (`cxx:77-143`).
    fn values(&mut self, x: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        let (mut u, mut v, mut t) = (0.0_f64, 0.0_f64, 0.0_f64);
        match self.fix {
            1 => {
                t = self.my_t;
                u = x.value(1);
                v = x.value(2);
            }
            2 => {
                t = x.value(1);
                u = self.my_u;
                v = x.value(2);
            }
            3 => {
                t = x.value(1);
                u = x.value(2);
                v = self.my_v;
            }
            _ => {}
        }

        let (c, dc1_t) = self.curve.d1(t);
        let (s, ds1_u, ds1_v, ds2_u, ds2_v, ds2_uv) = self.surface.d2(u, v);
        // `gp_Vec V(C, S)`
        let w = GpVec::from_pnts(&c, &s);

        f.set_value(1, w.dot(&ds1_u) * self.norm);
        f.set_value(2, w.dot(&ds1_v) * self.norm);

        match self.fix {
            1 => {
                d.set_value(1, 1, (ds1_u.square_magnitude() + w.dot(&ds2_u)) * self.norm);
                d.set_value(1, 2, (ds1_v.dot(&ds1_u) + w.dot(&ds2_uv)) * self.norm);
                d.set_value(2, 1, d.value(1, 2));
                d.set_value(2, 2, (ds1_v.square_magnitude() + w.dot(&ds2_v)) * self.norm);
            }
            2 => {
                d.set_value(1, 1, (-dc1_t.dot(&ds1_u)) * self.norm);
                d.set_value(1, 2, (ds1_v.dot(&ds1_u) + w.dot(&ds2_uv)) * self.norm);
                d.set_value(2, 1, (-dc1_t.dot(&ds1_v)) * self.norm);
                d.set_value(2, 2, (ds1_v.square_magnitude() + w.dot(&ds2_v)) * self.norm);
            }
            3 => {
                d.set_value(1, 1, -dc1_t.dot(&ds1_u) * self.norm);
                d.set_value(1, 2, (ds1_u.square_magnitude() + w.dot(&ds2_u)) * self.norm);
                d.set_value(2, 1, -dc1_t.dot(&ds1_v) * self.norm);
                d.set_value(2, 2, (ds1_v.dot(&ds1_u) + w.dot(&ds2_uv)) * self.norm);
            }
            _ => {}
        }

        self.my_u = u;
        self.my_v = v;
        self.my_t = t;
        true
    }
}

impl MathFunctionSetWithDerivatives for PrjFunc<'_> {
    fn nb_variables(&self) -> usize {
        2
    }

    fn nb_equations(&self) -> usize {
        2
    }

    fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        let mut d = MathMatrix::with_init(1, 2, 1, 2, 0.0);
        self.values(x, f, &mut d)
    }

    fn derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        let mut f = MathVector::with_init(1, 2, 0.0);
        self.values(x, &mut f, d)
    }

    fn values(&mut self, x: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        self.values(x, f, d)
    }
}
