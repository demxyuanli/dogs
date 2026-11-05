//! `Extrema_FuncExtCS` — the 3-equation orthogonality system between a curve
//! and a surface. Source: `Extrema_FuncExtCS.cxx:45-238`, `.hxx:35-100`.
//!
//! ```text
//! { F1(t,u,v) = (C(t)-S(u,v)).Dtc(t) }
//! { F2(t,u,v) = (C(t)-S(u,v)).Dus(u,v) }
//! { F3(t,u,v) = (C(t)-S(u,v)).Dvs(u,v) }
//! ```
//!
//! `Values` additionally fills the 3x3 Jacobian exactly as `cxx:150-165`
//! writes it, including the signed duplicates for the lower triangle.
//!
//! `GetStateNumber` (`cxx:172-197`) appends one solution per accepted root; the
//! duplicate test compares the curve parameter of the stored solutions against
//! the current `myt` within `Precision::SquarePConfusion()`.

use occt_math::{MathFunctionSetWithDerivatives, MathMatrix, MathVector};

use super::prelude::*;

/// `Precision::SquarePConfusion()` (`Precision.cxx`), used by
/// `GetStateNumber` (`cxx:178`).
const SQUARE_PCONFUSION: f64 = PCONFUSION * PCONFUSION;

/// `Extrema_FuncExtCS` (`Extrema_FuncExtCS.hxx:35-100`).
///
/// Holds the last evaluated `(t, u, v)` triple and its two points (`myt`,
/// `myU`, `myV`, `myP1`, `myP2`), plus the three result sequences
/// `mySqDist` / `myPoint1` / `myPoint2` that `GetStateNumber` fills.
pub struct FuncExtCS<'a> {
    /// `myC`.
    c: Option<&'a dyn Curve>,
    /// `myS`.
    s: Option<&'a dyn Surface>,
    /// `myP1`.
    p1: GpPnt,
    /// `myP2`.
    p2: GpPnt,
    /// `myt`.
    t: f64,
    /// `myU`.
    u: f64,
    /// `myV`.
    v: f64,
    /// `mySqDist`.
    sq_dist: Vec<f64>,
    /// `myPoint1`.
    point1: Vec<ExtremaPOnCurv>,
    /// `myPoint2`.
    point2: Vec<ExtremaPOnSurf>,
}

impl<'a> FuncExtCS<'a> {
    /// `Extrema_FuncExtCS()` (`cxx:47-57`): `myC` / `myS` null, `myt` / `myU`
    /// / `myV` zero, `myCinit` / `mySinit` false.
    pub fn new() -> Self {
        FuncExtCS {
            c: None,
            s: None,
            p1: GpPnt::new(0.0, 0.0, 0.0),
            p2: GpPnt::new(0.0, 0.0, 0.0),
            t: 0.0,
            u: 0.0,
            v: 0.0,
            sq_dist: Vec::new(),
            point1: Vec::new(),
            point2: Vec::new(),
        }
    }

    /// `Extrema_FuncExtCS(C, S)` (`cxx:61-65`): `Initialize(C, S)`.
    pub fn with_curve_surface(c: &'a dyn Curve, s: &'a dyn Surface) -> Self {
        let mut this = FuncExtCS::new();
        this.initialize(c, s);
        this
    }

    /// `Initialize(C, S)` (`cxx:69-78`): rebind both handles and clear the
    /// three sequences.
    pub fn initialize(&mut self, c: &'a dyn Curve, s: &'a dyn Surface) {
        self.c = Some(c);
        self.s = Some(s);
        self.point1.clear();
        self.point2.clear();
        self.sq_dist.clear();
    }

    fn handles(&self) -> (&'a dyn Curve, &'a dyn Surface) {
        // `cxx:96-100` / `:135-139` throw `Standard_TypeMismatch` when either
        // handle is unset.
        (
            self.c.expect("Extrema_FuncExtCS: Standard_TypeMismatch"),
            self.s.expect("Extrema_FuncExtCS: Standard_TypeMismatch"),
        )
    }

    /// `NbExt()` (`cxx:202-205`): `mySqDist.Length()`.
    pub fn nb_ext(&self) -> usize {
        self.sq_dist.len()
    }

    /// `SquareDistance(N)` (`cxx:209-216`).
    pub fn square_distance(&self, n: usize) -> f64 {
        self.sq_dist[n - 1]
    }

    /// `PointOnCurve(N)` (`cxx:220-227`).
    pub fn point_on_curve(&self, n: usize) -> ExtremaPOnCurv {
        self.point1[n - 1]
    }

    /// `PointOnSurface(N)` (`cxx:231-238`).
    pub fn point_on_surface(&self, n: usize) -> ExtremaPOnSurf {
        self.point2[n - 1]
    }

    /// Mutable access to `mySqDist` (`hxx:77-78`), which
    /// `Extrema_GenExtCS::Perform` rewrites in the false-extrema filter
    /// (`Extrema_GenExtCS.cxx:390-420`).
    pub fn square_distances_mut(&mut self) -> &mut Vec<f64> {
        &mut self.sq_dist
    }

    /// Mutable access to `myPoint1` (`hxx:80-81`).
    pub fn points_on_curve_mut(&mut self) -> &mut Vec<ExtremaPOnCurv> {
        &mut self.point1
    }

    /// Mutable access to `myPoint2` (`hxx:83-84`).
    pub fn points_on_surf_mut(&mut self) -> &mut Vec<ExtremaPOnSurf> {
        &mut self.point2
    }
}

impl Default for FuncExtCS<'_> {
    fn default() -> Self {
        FuncExtCS::new()
    }
}

impl MathFunctionSetWithDerivatives for FuncExtCS<'_> {
    /// `NbVariables()` (`cxx:82-85`): 3.
    fn nb_variables(&self) -> usize {
        3
    }

    /// `NbEquations()` (`cxx:89-92`): 3.
    fn nb_equations(&self) -> usize {
        3
    }

    /// `Value(UV, F)` (`cxx:96-121`): sets `myt` / `myU` / `myV` / `myP1` /
    /// `myP2`, then `F = (P2->P1)·(Dtc, Dus, Dvs)`.
    fn value(&mut self, uv: &MathVector, f: &mut MathVector) -> bool {
        let (c, s) = self.handles();
        self.t = uv.value(1);
        self.u = uv.value(2);
        self.v = uv.value(3);
        let (p1, dtc) = c.d1(self.t);
        let (p2, dus, dvs) = s.d1(self.u, self.v);
        self.p1 = p1;
        self.p2 = p2;
        let p1p2 = GpVec::from_pnts(&p2, &p1);
        f.set_value(1, p1p2.dot(&dtc));
        f.set_value(2, p1p2.dot(&dus));
        f.set_value(3, p1p2.dot(&dvs));
        true
    }

    /// `Derivatives(UV, DF)` (`cxx:125-129`): `Values` with a scratch `F`.
    fn derivatives(&mut self, uv: &MathVector, d: &mut MathMatrix) -> bool {
        let mut f = MathVector::new(1, 3);
        self.values(uv, &mut f, d)
    }

    /// `Values(UV, F, Df)` (`cxx:133-168`).
    fn values(&mut self, uv: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        let (c, s) = self.handles();
        self.t = uv.value(1);
        self.u = uv.value(2);
        self.v = uv.value(3);
        let (p1, dtc, dttc) = c.d2(self.t);
        let (p2, dus, dvs, duus, dvvs, duvs) = s.d2(self.u, self.v);
        self.p1 = p1;
        self.p2 = p2;
        let p1p2 = GpVec::from_pnts(&p2, &p1);

        f.set_value(1, p1p2.dot(&dtc));
        f.set_value(2, p1p2.dot(&dus));
        f.set_value(3, p1p2.dot(&dvs));

        d.set_value(1, 1, dtc.square_magnitude() + p1p2.dot(&dttc));
        d.set_value(1, 2, -dus.dot(&dtc));
        d.set_value(1, 3, -dvs.dot(&dtc));

        d.set_value(2, 1, -d.value(1, 2));
        d.set_value(2, 2, -dus.square_magnitude() + p1p2.dot(&duus));
        d.set_value(2, 3, -dvs.dot(&dus) + p1p2.dot(&duvs));

        d.set_value(3, 1, -d.value(1, 3));
        d.set_value(3, 2, d.value(2, 3));
        d.set_value(3, 3, -dvs.square_magnitude() + p1p2.dot(&dvvs));

        true
    }

    /// `GetStateNumber()` (`cxx:172-197`): skip the root when some stored
    /// solution already has the same curve parameter within
    /// `Precision::SquarePConfusion()`; otherwise append the squared distance
    /// and both points.
    fn state_number(&mut self) -> i32 {
        let tol2d = SQUARE_PCONFUSION;
        let nb_sol = self.sq_dist.len();
        let mut i = 0;
        while i < nb_sol {
            let mut a_t = self.point1[i].parameter();
            a_t -= self.t;
            a_t *= a_t;
            if a_t <= tol2d {
                break;
            }
            i += 1;
        }
        if i < nb_sol {
            return 0;
        }
        self.sq_dist.push(self.p1.square_distance(&self.p2));
        self.point1.push(ExtremaPOnCurv::new(self.t, self.p1));
        self.point2.push(ExtremaPOnSurf::new(self.u, self.v, self.p2));
        0
    }
}
