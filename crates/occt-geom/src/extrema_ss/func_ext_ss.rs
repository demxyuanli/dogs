//! `Extrema_FuncExtSS` and the `Extrema_FuncDistSS` objective of
//! `Extrema_GenExtSS`. Source: `Extrema_FuncExtSS.cxx` / `.hxx`,
//! `Extrema_GenExtSS.cxx:29-90`.
//!
//! `Extrema_FuncExtSS` is the 4-equation / 4-variable orthogonality system
//! solved by `math_FunctionSetRoot`:
//!
//! ```text
//! { F1 = (S1-S2)·Du1s1 } { F2 = (S1-S2)·Dv1s1 }
//! { F3 = (S1-S2)·Du2s2 } { F4 = (S1-S2)·Dv2s2 }
//! ```
//!
//! `Extrema_FuncDistSS` is the plain squared-distance objective (with its
//! analytic gradient) handed to `math_BFGS`.

use occt_core::gp::{GpPnt, GpVec};
use occt_math::{
    MathFunctionSetWithDerivatives, MathMatrix, MathMultipleVarFunction,
    MathMultipleVarFunctionWithGradient, MathVector,
};

use crate::extrema_surf::ExtremaPOnSurf;
use crate::surface::Surface;

/// `Extrema_FuncDistSS` (`Extrema_GenExtSS.cxx:29-90`).
pub struct ExtremaFuncDistSS<'a> {
    s1: &'a dyn Surface,
    s2: &'a dyn Surface,
}

impl<'a> ExtremaFuncDistSS<'a> {
    /// `Extrema_FuncDistSS(S1, S2)` (`Extrema_GenExtSS.cxx:33-38`).
    pub fn new(s1: &'a dyn Surface, s2: &'a dyn Surface) -> Self {
        ExtremaFuncDistSS { s1, s2 }
    }

    /// `NbVariables()` (`Extrema_GenExtSS.cxx:40`).
    pub fn nb_variables(&self) -> usize {
        4
    }

    /// `Value(X, F)` (`Extrema_GenExtSS.cxx:42-47`).
    pub fn value(&self, x: &MathVector) -> f64 {
        self.s1
            .d0(x.value(1), x.value(2))
            .square_distance(&self.s2.d0(x.value(3), x.value(4)))
    }

    /// `Gradient(X, G)` (`Extrema_GenExtSS.cxx:49-66`).
    ///
    /// `G(1) = P1P2·Du1s1`, `G(2) = P1P2·Dv1s1`, `G(3) = -P1P2·Du2s2`,
    /// `G(4) = -P1P2·Dv2s2`, with `P1P2 = gp_Vec(P2, P1)`.
    pub fn gradient(&self, x: &MathVector) -> MathVector {
        let (p1, du1, dv1) = self.s1.d1(x.value(1), x.value(2));
        let (p2, du2, dv2) = self.s2.d1(x.value(3), x.value(4));
        let p1p2 = GpVec::from_pnts(&p2, &p1);
        let mut g = MathVector::new(1, 4);
        g.set_value(1, p1p2.dot(&du1));
        g.set_value(2, p1p2.dot(&dv1));
        g.set_value(3, -p1p2.dot(&du2));
        g.set_value(4, -p1p2.dot(&dv2));
        g
    }

    /// `Values(X, F, G)` (`Extrema_GenExtSS.cxx:68-85`) — the same two
    /// expressions as `Value` + `Gradient`, computed together.
    pub fn values(&self, x: &MathVector) -> (f64, MathVector) {
        (self.value(x), self.gradient(x))
    }
}

impl MathMultipleVarFunction for ExtremaFuncDistSS<'_> {
    fn nb_variables(&self) -> usize {
        ExtremaFuncDistSS::nb_variables(self)
    }

    fn value(&self, x: &MathVector, f: &mut f64) -> bool {
        *f = ExtremaFuncDistSS::value(self, x);
        true
    }
}

impl MathMultipleVarFunctionWithGradient for ExtremaFuncDistSS<'_> {
    fn gradient(&self, x: &MathVector, g: &mut MathVector) -> bool {
        *g = ExtremaFuncDistSS::gradient(self, x);
        true
    }

    fn values(&self, x: &MathVector, f: &mut f64, g: &mut MathVector) -> bool {
        let (v, grad) = ExtremaFuncDistSS::values(self, x);
        *f = v;
        *g = grad;
        true
    }
}

/// `Extrema_FuncExtSS` (`Extrema_FuncExtSS.hxx:28-100`).
///
/// Holds the last evaluated `(U, V, P)` pair, exactly as the OCCT class does
/// (`myU1`/`myV1`/`myP1`, `myU2`/`myV2`/`myP2`), plus the three result
/// sequences `mySqDist`, `myPoint1`, `myPoint2` filled by `GetStateNumber`.
pub struct ExtremaFuncExtSS<'a> {
    s1: Option<&'a dyn Surface>,
    s2: Option<&'a dyn Surface>,
    u1: f64,
    v1: f64,
    u2: f64,
    v2: f64,
    p1: GpPnt,
    p2: GpPnt,
    /// `mySqDist`.
    sq_dist: Vec<f64>,
    /// `myPoint1`.
    point1: Vec<ExtremaPOnSurf>,
    /// `myPoint2`.
    point2: Vec<ExtremaPOnSurf>,
}

impl<'a> ExtremaFuncExtSS<'a> {
    /// `Extrema_FuncExtSS()` (`Extrema_FuncExtSS.cxx:52-64`): both surfaces
    /// null, `myS1init` / `myS2init` false.
    pub fn new() -> Self {
        ExtremaFuncExtSS {
            s1: None,
            s2: None,
            u1: 0.0,
            v1: 0.0,
            u2: 0.0,
            v2: 0.0,
            p1: GpPnt::new(0.0, 0.0, 0.0),
            p2: GpPnt::new(0.0, 0.0, 0.0),
            sq_dist: Vec::new(),
            point1: Vec::new(),
            point2: Vec::new(),
        }
    }

    /// `Extrema_FuncExtSS(S1, S2)` (`Extrema_FuncExtSS.cxx:66-78`).
    pub fn with_surfaces(s1: &'a dyn Surface, s2: &'a dyn Surface) -> Self {
        ExtremaFuncExtSS {
            s1: Some(s1),
            s2: Some(s2),
            ..ExtremaFuncExtSS::new()
        }
    }

    /// `Extrema_FuncExtSS::Initialize(S1, S2)` (`Extrema_FuncExtSS.cxx:82-91`):
    /// rebinds the surfaces and clears the three sequences.
    pub fn initialize(&mut self, s1: &'a dyn Surface, s2: &'a dyn Surface) {
        self.s1 = Some(s1);
        self.s2 = Some(s2);
        self.point1.clear();
        self.point2.clear();
        self.sq_dist.clear();
    }

    fn surfaces(&self) -> (&'a dyn Surface, &'a dyn Surface) {
        // `Extrema_FuncExtSS.cxx:112-115` throws `Standard_TypeMismatch` when
        // either surface is unset.
        (
            self.s1.expect("Extrema_FuncExtSS: Standard_TypeMismatch"),
            self.s2.expect("Extrema_FuncExtSS: Standard_TypeMismatch"),
        )
    }

    /// `NbExt()` (`Extrema_FuncExtSS.cxx:207-210`): `mySqDist.Length()`.
    pub fn nb_ext(&self) -> usize {
        self.sq_dist.len()
    }

    /// `SquareDistance(N)` (`Extrema_FuncExtSS.cxx:214-220`).
    pub fn square_distance(&self, n: i32) -> f64 {
        self.sq_dist[(n - 1) as usize]
    }

    /// `PointOnS1(N)` (`Extrema_FuncExtSS.cxx:224-231`).
    pub fn point_on_s1(&self, n: i32) -> ExtremaPOnSurf {
        self.point1[(n - 1) as usize]
    }

    /// `PointOnS2(N)` (`Extrema_FuncExtSS.cxx:234-241`).
    pub fn point_on_s2(&self, n: i32) -> ExtremaPOnSurf {
        self.point2[(n - 1) as usize]
    }
}

impl Default for ExtremaFuncExtSS<'_> {
    fn default() -> Self {
        ExtremaFuncExtSS::new()
    }
}

impl MathFunctionSetWithDerivatives for ExtremaFuncExtSS<'_> {
    /// `NbVariables()` (`Extrema_FuncExtSS.cxx:96-99`).
    fn nb_variables(&self) -> usize {
        4
    }

    /// `NbEquations()` (`Extrema_FuncExtSS.cxx:103-106`).
    fn nb_equations(&self) -> usize {
        4
    }

    /// `Extrema_FuncExtSS::Value` (`Extrema_FuncExtSS.cxx:109-133`).
    ///
    /// `F(1) = P1P2·Du1s1`, `F(2) = P1P2·Dv1s1`, `F(3) = P1P2·Du2s2`,
    /// `F(4) = P1P2·Dv2s2` with `P1P2 = gp_Vec(P2, P1)`. (Note the sign
    /// difference against `Extrema_FuncDistSS::Gradient`, which minimises the
    /// squared distance instead of driving these to zero.)
    fn value(&mut self, uv: &MathVector, f: &mut MathVector) -> bool {
        let (s1, s2) = self.surfaces();
        self.u1 = uv.value(1);
        self.v1 = uv.value(2);
        self.u2 = uv.value(3);
        self.v2 = uv.value(4);
        let (p1, du1, dv1) = s1.d1(self.u1, self.v1);
        let (p2, du2, dv2) = s2.d1(self.u2, self.v2);
        self.p1 = p1;
        self.p2 = p2;
        let p1p2 = GpVec::from_pnts(&p2, &p1);
        f.set_value(1, p1p2.dot(&du1));
        f.set_value(2, p1p2.dot(&dv1));
        f.set_value(3, p1p2.dot(&du2));
        f.set_value(4, p1p2.dot(&dv2));
        true
    }

    /// `Extrema_FuncExtSS::Derivatives` (`Extrema_FuncExtSS.cxx:137-140`):
    /// computes `Values` and discards `F`.
    fn derivatives(&mut self, uv: &MathVector, d: &mut MathMatrix) -> bool {
        let mut f = MathVector::new(1, 4);
        self.values(uv, &mut f, d)
    }

    /// `Extrema_FuncExtSS::Values` (`Extrema_FuncExtSS.cxx:145-187`).
    ///
    /// The Jacobian uses the second derivatives `D2` of both surfaces; the
    /// entries `(2,1)`, `(3,1)`, `(3,2)`, `(4,1)`, `(4,2)`, `(4,3)` are the
    /// signed duplicates the `.cxx` writes directly.
    fn values(&mut self, uv: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        let (s1, s2) = self.surfaces();
        self.u1 = uv.value(1);
        self.v1 = uv.value(2);
        self.u2 = uv.value(3);
        self.v2 = uv.value(4);
        let (p1, du1, dv1, du1u1, dv1v1, du1v1) = s1.d2(self.u1, self.v1);
        let (p2, du2, dv2, du2u2, dv2v2, du2v2) = s2.d2(self.u2, self.v2);
        self.p1 = p1;
        self.p2 = p2;
        let p1p2 = GpVec::from_pnts(&p2, &p1);

        f.set_value(1, p1p2.dot(&du1));
        f.set_value(2, p1p2.dot(&dv1));
        f.set_value(3, p1p2.dot(&du2));
        f.set_value(4, p1p2.dot(&dv2));

        d.set_value(1, 1, du1.square_magnitude() + p1p2.dot(&du1u1));
        d.set_value(1, 2, dv1.dot(&du1) + p1p2.dot(&du1v1));
        d.set_value(1, 3, -du2.dot(&du1));
        d.set_value(1, 4, -dv2.dot(&du1));

        d.set_value(2, 1, d.value(1, 2));
        d.set_value(2, 2, dv1.square_magnitude() + p1p2.dot(&dv1v1));
        d.set_value(2, 3, -du2.dot(&dv1));
        d.set_value(2, 4, -dv2.dot(&dv1));

        d.set_value(3, 1, -d.value(1, 3));
        d.set_value(3, 2, -d.value(2, 3));
        d.set_value(3, 3, -du2.square_magnitude() + p1p2.dot(&du2u2));
        d.set_value(3, 4, -dv2.dot(&du2) + p1p2.dot(&du2v2));

        d.set_value(4, 1, -d.value(1, 4));
        d.set_value(4, 2, -d.value(2, 4));
        d.set_value(4, 3, d.value(3, 4));
        d.set_value(4, 4, -dv2.square_magnitude() + p1p2.dot(&dv2v2));

        true
    }

    /// `Extrema_FuncExtSS::GetStateNumber` (`Extrema_FuncExtSS.cxx:192-202`):
    /// appends `myP1.SquareDistance(myP2)`, `(myU1, myV1, myP1)` and
    /// `(myU2, myV2, myP2)`.
    fn get_state_number(&mut self) -> i32 {
        let (_, _) = self.surfaces();
        self.sq_dist.push(self.p1.square_distance(&self.p2));
        self.point1.push(ExtremaPOnSurf::new(self.u1, self.v1, self.p1));
        self.point2.push(ExtremaPOnSurf::new(self.u2, self.v2, self.p2));
        0
    }
}
