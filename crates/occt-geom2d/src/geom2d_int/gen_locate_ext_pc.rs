//! Port of `Extrema_GenLocateExtPC` (`Extrema_GenLocateExtPC.hxx:37-167`)
//! monomorphised onto the `Geom2dInt` instantiation
//! (`Geom2dInt_TheLocateExtPCOfTheProjPCurOfGInter.hxx:28-33`):
//! `TheCurve=Curve2d`, `TheTool=Geom2dInt_Geom2dCurveTool`, `ThePOnC=POnCurv2d`,
//! `ThePnt=GpPnt2d`, `ThePCLocF=GFuncExtPC2d`.

use occt_core::gp::GpPnt2d;
use occt_math::{MathFunctionRoot, FunctionWithDerivative};

use crate::curve::Curve2d;
use super::gfunc_ext_pc::{GFuncExtPC2d, POnCurv2d};

/// `Extrema_GenLocateExtPC` (`hxx:42-167`), 2D instantiation.
pub struct GenLocateExtPC2d<'a> {
    done: bool,
    tol_u: f64,
    u_min: f64,
    u_sup: f64,
    f: GFuncExtPC2d<'a>,
}

impl<'a> Default for GenLocateExtPC2d<'a> {
    /// `Extrema_GenLocateExtPC()` (`hxx:48-54`).
    fn default() -> Self {
        Self {
            done: false,
            tol_u: 0.0,
            u_min: 0.0,
            u_sup: 0.0,
            f: GFuncExtPC2d::default(),
        }
    }
}

impl<'a> GenLocateExtPC2d<'a> {
    /// `Extrema_GenLocateExtPC(P, C, U0, TolU)` (`hxx:64-71`): the search
    /// window is the curve's own `[FirstParameter, LastParameter]`.
    pub fn with_params(p: &GpPnt2d, c: &'a dyn Curve2d, u0: f64, tol_u: f64) -> Self {
        let mut r = Self::default();
        r.initialize(c, c.first_parameter(), c.last_parameter(), tol_u);
        r.perform(p, u0);
        r
    }

    /// `Initialize(theC, Umin, Usup, TolU)` (`hxx:94-104`).
    pub fn initialize(&mut self, c: &'a dyn Curve2d, u_min: f64, u_sup: f64, tol_u: f64) {
        self.done = false;
        self.f.initialize(c);
        self.u_min = u_min;
        self.u_sup = u_sup;
        self.tol_u = tol_u;
    }

    /// `Perform(theP, U0)` (`hxx:108-126`).
    pub fn perform(&mut self, p: &GpPnt2d, u0: f64) {
        self.f.set_point(p);
        let s = MathFunctionRoot::new_with_bounds(
            &mut self.f,
            u0,
            self.tol_u,
            self.u_min,
            self.u_sup,
            100,
        );
        self.done = s.is_done();
        if self.done {
            let uu = self.point(1).parameter();
            let mut ff = 0.0f64;
            if self.f.value(uu, &mut ff) {
                if ff.abs() >= 1.0e-07 {
                    self.done = false;
                }
            } else {
                self.done = false;
            }
        }
    }

    /// `IsDone` (`hxx:129`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `SquareDistance` (`hxx:132-139`).
    pub fn square_distance(&self) -> f64 {
        self.f.square_distance(1)
    }

    /// `IsMin` (`hxx:142-149`).
    pub fn is_min(&self) -> bool {
        self.f.is_min(1)
    }

    /// `Point` (`hxx:152-159`), 1-based.
    pub fn point(&self, n: usize) -> POnCurv2d {
        self.f.point(n)
    }
}
