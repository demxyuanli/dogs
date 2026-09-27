//! Port of `Extrema_GenLocateExtPC` (`Extrema_GenLocateExtPC.hxx:37-167`)
//! monomorphised onto the `Geom2dInt` instantiation
//! (`Geom2dInt_TheLocateExtPCOfTheProjPCurOfGInter.hxx:28-33`):
//! `TheCurve=Curve2d`, `TheTool=Geom2dInt_Geom2dCurveTool`, `ThePOnC=POnCurv2d`,
//! `ThePnt=GpPnt2d`, `ThePCLocF=GFuncExtPC`.

use occt_core::gp::GpPnt2d;
use occt_math::{MathFunctionRoot, MathFunctionWithDerivative};

use crate::curve::Curve2d;
use super::gfunc_ext_pc::{GFuncExtPC, POnCurv2d};

/// `Extrema_GenLocateExtPC` (`hxx:42-167`), 2D instantiation.
pub struct GenLocateExtPC<'a> {
    my_done: bool,
    my_tol_u: f64,
    my_u_min: f64,
    my_u_sup: f64,
    my_f: GFuncExtPC<'a>,
}

impl<'a> Default for GenLocateExtPC<'a> {
    /// `Extrema_GenLocateExtPC()` (`hxx:48-54`).
    fn default() -> Self {
        Self {
            my_done: false,
            my_tol_u: 0.0,
            my_u_min: 0.0,
            my_u_sup: 0.0,
            my_f: GFuncExtPC::default(),
        }
    }
}

impl<'a> GenLocateExtPC<'a> {
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
        self.my_done = false;
        self.my_f.initialize(c);
        self.my_u_min = u_min;
        self.my_u_sup = u_sup;
        self.my_tol_u = tol_u;
    }

    /// `Perform(theP, U0)` (`hxx:108-126`).
    pub fn perform(&mut self, p: &GpPnt2d, u0: f64) {
        self.my_f.set_point(p);
        let s = MathFunctionRoot::new_with_bounds(
            &mut self.my_f,
            u0,
            self.my_tol_u,
            self.my_u_min,
            self.my_u_sup,
            100,
        );
        self.my_done = s.is_done();
        if self.my_done {
            let uu = self.point(1).parameter();
            let mut ff = 0.0f64;
            if self.my_f.value(uu, &mut ff) {
                if ff.abs() >= 1.0e-07 {
                    self.my_done = false;
                }
            } else {
                self.my_done = false;
            }
        }
    }

    /// `IsDone` (`hxx:129`).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// `SquareDistance` (`hxx:132-139`).
    pub fn square_distance(&self) -> f64 {
        self.my_f.square_distance(1)
    }

    /// `IsMin` (`hxx:142-149`).
    pub fn is_min(&self) -> bool {
        self.my_f.is_min(1)
    }

    /// `Point` (`hxx:152-159`), 1-based.
    pub fn point(&self, n: usize) -> POnCurv2d {
        self.my_f.point(n)
    }
}
