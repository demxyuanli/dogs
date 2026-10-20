use occt_math::{MathFunctionSetWithDerivatives, MathMatrix, MathVector};

use super::super::prelude::*;
use super::point_params::ExtremaPOnSurf;

// ---------------------------------------------------------------------------
// Extrema_FuncPSNorm.hxx / .cxx (194 lines)
// ---------------------------------------------------------------------------

/// Extrema_FuncPSNorm (Extrema_FuncPSNorm.hxx:54-104,
/// Extrema_FuncPSNorm.cxx): the math_FunctionSetWithDerivatives whose two
/// equations are F1 = (S-P) dot Su, F2 = (S-P) dot Sv (cxx:83-135) and which
/// accumulates the deduplicated extrema in GetStateNumber (cxx:139-165).
pub struct ExtremaFuncPsNorm<'a> {
    /// myP.
    p: GpPnt,
    /// myS.
    s: Option<&'a dyn Surface>,
    /// myU / myV (cxx:89-90, cxx:118-119).
    u: f64,
    v: f64,
    /// myPs.
    ps: GpPnt,
    /// mySqDist (Extrema_FuncPSNorm.cxx:100).
    sq_dist: Vec<f64>,
    /// myPoint (Extrema_FuncPSNorm.cxx:101).
    point: Vec<ExtremaPOnSurf>,
    /// myPinit / mySinit (Extrema_FuncPSNorm.cxx:102-103).
    p_init: bool,
    s_init: bool,
}

impl<'a> Default for ExtremaFuncPsNorm<'a> {
    /// Extrema_FuncPSNorm() (Extrema_FuncPSNorm.cxx:26-33).
    fn default() -> Self {
        ExtremaFuncPsNorm {
            p: GpPnt::zero(),
            s: None,
            u: 0.0,
            v: 0.0,
            ps: GpPnt::zero(),
            sq_dist: Vec::new(),
            point: Vec::new(),
            p_init: false,
            s_init: false,
        }
    }
}

impl<'a> ExtremaFuncPsNorm<'a> {
    /// Extrema_FuncPSNorm(P, S) (Extrema_FuncPSNorm.cxx:36-44).
    pub fn new(p: GpPnt, s: &'a dyn Surface) -> Self {
        ExtremaFuncPsNorm {
            p,
            s: Some(s),
            u: 0.0,
            v: 0.0,
            ps: GpPnt::zero(),
            sq_dist: Vec::new(),
            point: Vec::new(),
            p_init: true,
            s_init: true,
        }
    }

    /// Extrema_FuncPSNorm::Initialize(S) (Extrema_FuncPSNorm.cxx:47-53).
    pub fn initialize(&mut self, s: &'a dyn Surface) {
        self.s = Some(s);
        self.s_init = true;
        self.point.clear();
        self.sq_dist.clear();
    }

    /// Extrema_FuncPSNorm::SetPoint(P) (Extrema_FuncPSNorm.cxx:57-63).
    pub fn set_point(&mut self, p: GpPnt) {
        self.p = p;
        self.p_init = true;
        self.point.clear();
        self.sq_dist.clear();
    }

    /// Extrema_FuncPSNorm::NbExt (Extrema_FuncPSNorm.cxx:169-172).
    pub fn nb_ext(&self) -> usize {
        self.sq_dist.len()
    }

    /// Extrema_FuncPSNorm::SquareDistance(N) (Extrema_FuncPSNorm.cxx:176-183).
    /// N is 1-based; OCCT raises Standard_OutOfRange on a bad index.
    pub fn square_distance(&self, n: usize) -> f64 {
        self.sq_dist[n - 1]
    }

    /// Extrema_FuncPSNorm::Point(N) (Extrema_FuncPSNorm.cxx:187-194).
    pub fn point(&self, n: usize) -> ExtremaPOnSurf {
        self.point[n - 1]
    }
}

impl<'a> MathFunctionSetWithDerivatives for ExtremaFuncPsNorm<'a> {
    /// Extrema_FuncPSNorm::NbVariables (Extrema_FuncPSNorm.cxx:69-72).
    fn nb_variables(&self) -> usize {
        2
    }

    /// Extrema_FuncPSNorm::NbEquations (Extrema_FuncPSNorm.cxx:76-79).
    fn nb_equations(&self) -> usize {
        2
    }

    /// Extrema_FuncPSNorm::Value(UV, F) (Extrema_FuncPSNorm.cxx:83-100).
    /// The Standard_TypeMismatch guard is cxx:85-88; the Rust port reports the
    /// same failure as false (the solver's "no value" answer) because the trait
    /// has no exception channel.
    fn value(&mut self, uv: &MathVector, f: &mut MathVector) -> bool {
        let Some(s) = self.s else {
            return false;
        };
        if !self.p_init || !self.s_init {
            return false;
        }
        self.u = uv.value(1);
        self.v = uv.value(2);
        let (ps, dus, dvs) = s.d1(self.u, self.v);
        self.ps = ps;
        let pps = GpVec::from_pnts(&self.p, &ps);
        f.set_value(1, pps.dot(&dus));
        f.set_value(2, pps.dot(&dvs));
        true
    }

    /// Extrema_FuncPSNorm::Derivatives(UV, Df) (Extrema_FuncPSNorm.cxx:104-108).
    fn derivatives(&mut self, uv: &MathVector, d: &mut MathMatrix) -> bool {
        let mut f = MathVector::new(1, 2);
        self.values(uv, &mut f, d)
    }

    /// Extrema_FuncPSNorm::Values(UV, F, Df) (Extrema_FuncPSNorm.cxx:112-135).
    fn values(&mut self, uv: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        let Some(s) = self.s else {
            return false;
        };
        if !self.p_init || !self.s_init {
            return false;
        }
        self.u = uv.value(1);
        self.v = uv.value(2);
        let (ps, dus, dvs, duus, dvvs, duvs) = s.d2(self.u, self.v);
        self.ps = ps;
        let pps = GpVec::from_pnts(&self.p, &ps);
        let d11 = dus.square_magnitude() + pps.dot(&duus);
        let d12 = dvs.dot(&dus) + pps.dot(&duvs);
        let d22 = dvs.square_magnitude() + pps.dot(&dvvs);
        d.set_value(1, 1, d11);
        d.set_value(1, 2, d12);
        d.set_value(2, 1, d12);
        d.set_value(2, 2, d22);
        f.set_value(1, pps.dot(&dus));
        f.set_value(2, pps.dot(&dvs));
        true
    }

    /// Extrema_FuncPSNorm::GetStateNumber (Extrema_FuncPSNorm.cxx:139-165):
    /// append the current stationary point unless it is within
    /// Precision::PConfusion()^2 of an already stored one.
    fn state_number(&mut self) -> i32 {
        let tol2d = PCONFUSION * PCONFUSION;
        let nb_sol = self.sq_dist.len();
        let mut i = 0usize;
        while i < nb_sol {
            let (au, av) = self.point[i].parameter();
            if (self.u - au) * (self.u - au) + (self.v - av) * (self.v - av) <= tol2d {
                break;
            }
            i += 1;
        }
        if i < nb_sol {
            // cxx:158-161: already stored.
            return 0;
        }
        self.sq_dist.push(self.ps.square_distance(&self.p));
        self.point.push(ExtremaPOnSurf::new(self.u, self.v, self.ps));
        0
    }
}
