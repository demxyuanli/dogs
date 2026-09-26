//! Extrema_GenExtPS - general point/surface extrema engine.
//!
//! Source: Extrema_GenExtPS.cxx (1195 lines) + Extrema_GenExtPS.hxx (181),
//! src/ModelingData/TKGeomBase/Extrema/. Line references below are to the
//! .cxx unless another file is named. The helper classes it needs live here too:
//!
//! * [ExtremaPOnSurf] - Extrema_POnSurf.hxx:26-71.
//! * [ExtremaPOnSurfParams] - Extrema_POnSurfParams.hxx:29-88.
//! * [ExtremaElementType] - Extrema_ElementType.hxx:20-26.
//! * [ExtremaFuncPsNorm] - the myF functional, Extrema_FuncPSNorm.hxx /
//!   .cxx (194 lines); this is the math_FunctionSetWithDerivatives the Newton
//!   refinement runs on.
//!
//! The control flow is ported branch for branch:
//! Initialize (cxx:275-319), GetGridPoints + fillParams (cxx:321-454),
//! ComputeEdgeParameters (cxx:460-526), BuildGrid (cxx:528-753), LengthOfIso /
//! CorrectNbSamples (cxx:755-854), BuildTree (cxx:856-931), FindSolution
//! (cxx:933-952), SetFlag / SetAlgo (cxx:954-966), Perform (cxx:968-1151) and
//! the accessors (cxx:1155-1193).
//!
//! Two implementation devices are **not** control-flow differences:
//!
//! * GeomGridEval_Surface::EvaluateGrid (cxx:569-572, cxx:914-917) is an
//!   evaluator with span caching; the grid it produces is exactly
//!   S.Value(U_i, V_j) for every pair, so the port calls
//!   [crate::surface::Surface::d0] node by node.
//! * The NCollection_UBTree built by BuildTree (cxx:907-930) is an acceleration
//!   structure. Bnd_SphereUBTreeSelectorMin/Max (cxx:77-169) selects the grid
//!   sphere minimizing / maximizing Bnd_Sphere::Distance(P); the port keeps
//!   OCCT's mySphereArray (cxx:911) and performs the equivalent scan in
//!   [ExtremaGenExtPs::perform] instead of maintaining the tree. Bnd_Sphere
//!   itself is occt_core::bnd::BndSphere.
//!
//! Extrema_ExtPExtS / Extrema_ExtPRevS are ported in
//! extrusion_point_extrema.rs / revolution_point_extrema.rs and dispatched from
//! point_surface_extrema.rs (Extrema_ExtPS.cxx:292-343); they drive this engine
//! whenever their basis curve is not analytically computable.

use super::prelude::*;
use super::*;

use occt_core::bnd::BndSphere;
use occt_math::{MathFunctionSetRoot, MathFunctionSetWithDerivatives, MathMatrix, MathVector};

/// Precision::Infinite() (Precision.hxx:371) - the default of Extrema_POnSurf
/// (Extrema_POnSurf.hxx:33-36).
const PRECISION_INFINITE: f64 = occt_core::precision::INFINITE;

// ---------------------------------------------------------------------------
// Extrema_ElementType.hxx:20-26
// ---------------------------------------------------------------------------

/// Extrema_ElementType (Extrema_ElementType.hxx:20-26): which grid element a
/// [ExtremaPOnSurfParams] sits on.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub enum ExtremaElementType {
    #[default]
    Node,
    UIsoEdge,
    VIsoEdge,
    Face,
}

// ---------------------------------------------------------------------------
// Extrema_POnSurf.hxx:26-71
// ---------------------------------------------------------------------------

/// Extrema_POnSurf (Extrema_POnSurf.hxx:26-71): a 3d point and its (U, V)
/// parameters.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExtremaPOnSurf {
    u: f64,
    v: f64,
    p: GpPnt,
}

impl ExtremaPOnSurf {
    /// Extrema_POnSurf(theU, theV, theP) (Extrema_POnSurf.hxx:41-46).
    pub fn new(u: f64, v: f64, p: GpPnt) -> Self {
        ExtremaPOnSurf { u, v, p }
    }

    /// Value() (Extrema_POnSurf.hxx:49).
    pub fn value(&self) -> GpPnt {
        self.p
    }

    /// SetParameters (Extrema_POnSurf.hxx:53-58).
    pub fn set_parameters(&mut self, u: f64, v: f64, p: GpPnt) {
        self.u = u;
        self.v = v;
        self.p = p;
    }

    /// Parameter(U, V) (Extrema_POnSurf.hxx:61-65).
    pub fn parameter(&self) -> (f64, f64) {
        (self.u, self.v)
    }
}

impl Default for ExtremaPOnSurf {
    /// Extrema_POnSurf() (Extrema_POnSurf.hxx:32-37): (Infinite, Infinite) and
    /// P(Infinite, Infinite, Infinite).
    fn default() -> Self {
        ExtremaPOnSurf {
            u: PRECISION_INFINITE,
            v: PRECISION_INFINITE,
            p: GpPnt::new(PRECISION_INFINITE, PRECISION_INFINITE, PRECISION_INFINITE),
        }
    }
}

// ---------------------------------------------------------------------------
// Extrema_POnSurfParams.hxx:29-88
// ---------------------------------------------------------------------------

/// Extrema_POnSurfParams (Extrema_POnSurfParams.hxx:29-88): an [ExtremaPOnSurf]
/// plus the square distance to the projected point, the element type and the
/// grid indices it belongs to.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExtremaPOnSurfParams {
    point: ExtremaPOnSurf,
    sqr_distance: f64,
    element_type: ExtremaElementType,
    index_u: i32,
    index_v: i32,
}

impl ExtremaPOnSurfParams {
    /// Extrema_POnSurfParams() (Extrema_POnSurfParams.hxx:35-41).
    pub fn new_default() -> Self {
        ExtremaPOnSurfParams {
            point: ExtremaPOnSurf::default(),
            sqr_distance: 0.0,
            element_type: ExtremaElementType::Node,
            index_u: 0,
            index_v: 0,
        }
    }

    /// Extrema_POnSurfParams(theU, theV, thePnt)
    /// (Extrema_POnSurfParams.hxx:45-52).
    pub fn new(u: f64, v: f64, p: GpPnt) -> Self {
        ExtremaPOnSurfParams {
            point: ExtremaPOnSurf::new(u, v, p),
            sqr_distance: 0.0,
            element_type: ExtremaElementType::Node,
            index_u: 0,
            index_v: 0,
        }
    }

    /// SetSqrDistance (Extrema_POnSurfParams.hxx:56).
    pub fn set_sqr_distance(&mut self, d: f64) {
        self.sqr_distance = d;
    }

    /// GetSqrDistance (Extrema_POnSurfParams.hxx:59).
    pub fn get_sqr_distance(&self) -> f64 {
        self.sqr_distance
    }

    /// SetElementType (Extrema_POnSurfParams.hxx:62).
    pub fn set_element_type(&mut self, t: ExtremaElementType) {
        self.element_type = t;
    }

    /// GetElementType (Extrema_POnSurfParams.hxx:65).
    pub fn get_element_type(&self) -> ExtremaElementType {
        self.element_type
    }

    /// SetIndices (Extrema_POnSurfParams.hxx:69-73).
    pub fn set_indices(&mut self, iu: i32, iv: i32) {
        self.index_u = iu;
        self.index_v = iv;
    }

    /// GetIndices (Extrema_POnSurfParams.hxx:77-81).
    pub fn get_indices(&self) -> (i32, i32) {
        (self.index_u, self.index_v)
    }

    /// Extrema_POnSurf::Value() through the inheritance.
    pub fn value(&self) -> GpPnt {
        self.point.value()
    }

    /// Extrema_POnSurf::SetParameters() through the inheritance.
    pub fn set_parameters(&mut self, u: f64, v: f64, p: GpPnt) {
        self.point.set_parameters(u, v, p);
    }

    /// Extrema_POnSurf::Parameter() through the inheritance.
    pub fn parameter(&self) -> (f64, f64) {
        self.point.parameter()
    }
}

impl Default for ExtremaPOnSurfParams {
    fn default() -> Self {
        Self::new_default()
    }
}

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
    fn get_state_number(&mut self) -> i32 {
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

// ---------------------------------------------------------------------------
// Helpers (Extrema_GenExtPS.cxx:321-369, :755-854)
// ---------------------------------------------------------------------------

/// The distinct knot values of a flat knot vector - Geom_BSplineSurface::UKnots()
/// / Geom_BSplineCurve::Knots(), which OCCT stores as distinct knots plus
/// multiplicities (Geom_BSplineSurface.cxx:1148-1165). The port's surfaces keep
/// the flat sequence, so collapsing equal neighbours recovers the OCCT array.
fn distinct_knots(flat: &[f64]) -> Vec<f64> {
    let mut out: Vec<f64> = Vec::new();
    for &k in flat {
        if out.last() != Some(&k) {
            out.push(k);
        }
    }
    out
}

/// fillParams (Extrema_GenExtPS.cxx:321-369): parametric samples derived from
/// the knot spans of a B-spline / Bezier geometry. Returns None where OCCT
/// leaves theParams null (too few points, cxx:360-363).
fn fill_params(
    knots: &[f64],
    degree: i32,
    par_min: f64,
    par_max: f64,
    sample: i32,
) -> Option<Vec<f64>> {
    let mut params: Vec<f64> = Vec::new();
    let mut prev_par = par_min;
    params.push(prev_par);
    let mut i = 1usize;
    while i < knots.len() && knots[i - 1] < (par_max - PCONFUSION) {
        if knots[i] < par_min + PCONFUSION {
            i += 1;
            continue;
        }
        let step = (knots[i] - knots[i - 1]) / degree.max(2) as f64;
        let mut k = 1;
        while k <= degree {
            let par = knots[i - 1] + k as f64 * step;
            if par > par_max - PCONFUSION {
                break;
            }
            if par > prev_par + PCONFUSION {
                params.push(par);
                prev_par = par;
            }
            k += 1;
        }
        i += 1;
    }
    params.push(par_max);
    let nb_par = params.len();
    if nb_par < sample as usize {
        return None;
    }
    Some(params)
}

/// LengthOfIso (Extrema_GenExtPS.cxx:755-791): the polyline length of the iso
/// line theIso at parameter thePar, sampled in theNbPnts points from thePar1 to
/// thePar2.
fn length_of_iso(
    s: &dyn Surface,
    iso_u: bool,
    par1: f64,
    par2: f64,
    nb_pnts: i32,
    par: f64,
) -> f64 {
    let mut len = 0.0;
    let dpar = (par2 - par1) / (nb_pnts - 1) as f64;
    let mut a_par = par1 + dpar;
    let mut p1 = if iso_u { s.d0(par, par1) } else { s.d0(par1, par) };
    for _i in 2..=nb_pnts {
        let p2 = if iso_u { s.d0(par, a_par) } else { s.d0(a_par, par) };
        len += p1.distance(&p2);
        p1 = p2;
        a_par += dpar;
    }
    len
}

/// CorrectNbSamples (Extrema_GenExtPS.cxx:793-854). Faithful, including the
/// cxx:846-853 branch which scales theNbV under the aRatio < 0.1 test as the
/// OCCT source does.
fn correct_nb_samples(
    s: &dyn Surface,
    u1: f64,
    u2: f64,
    nb_u: &mut i32,
    v1: f64,
    v2: f64,
    nb_v: &mut i32,
) {
    let min_len = 1.0e-3;
    let nbp = 23.min(*nb_v);
    let mut len_u1 = length_of_iso(s, true, v1, v2, nbp, u1);
    if len_u1 <= min_len {
        let l = length_of_iso(s, true, v1, v2, nbp, 0.7 * u1 + 0.3 * u2);
        len_u1 = l.max(len_u1);
    }
    let mut len_u2 = length_of_iso(s, true, v1, v2, nbp, u2);
    if len_u2 <= min_len {
        let l = length_of_iso(s, true, v1, v2, nbp, 0.3 * u1 + 0.7 * u2);
        len_u2 = l.max(len_u2);
    }
    let nbp = 23.min(*nb_v);
    let mut len_v1 = length_of_iso(s, false, u1, u2, nbp, v1);
    if len_v1 <= min_len {
        let l = length_of_iso(s, false, u1, u2, nbp, 0.7 * v1 + 0.3 * v2);
        len_v1 = l.max(len_v1);
    }
    let mut len_v2 = length_of_iso(s, false, u1, u2, nbp, v2);
    if len_v2 <= min_len {
        let l = length_of_iso(s, false, u1, u2, nbp, 0.3 * v1 + 0.7 * v2);
        len_v2 = l.max(len_v2);
    }

    let step_v1 = len_u1 / *nb_v as f64;
    let step_v2 = len_u2 / *nb_v as f64;
    let step_u1 = len_v1 / *nb_u as f64;
    let step_u2 = len_v2 / *nb_u as f64;

    let max_step_v = step_v1.max(step_v2);
    let max_step_u = step_u1.max(step_u2);

    let ratio = max_step_v / max_step_u;
    if ratio > 10.0 {
        let mult = real_to_int(ratio.ln());
        if mult > 1 {
            *nb_v *= mult;
        }
    } else if ratio < 0.1 {
        let mult = real_to_int(-ratio.ln());
        if mult > 1 {
            *nb_v *= mult;
        }
    }
}

/// RealToInt (Standard_Real.hxx:319-328): truncation toward zero with the int
/// range clamp.
fn real_to_int(v: f64) -> i32 {
    if v < i32::MIN as f64 {
        i32::MIN
    } else if v > i32::MAX as f64 {
        i32::MAX
    } else {
        v as i32
    }
}

// ---------------------------------------------------------------------------
// Extrema_GenExtPS
// ---------------------------------------------------------------------------

/// Extrema_GenExtPS (Extrema_GenExtPS.hxx:36-179).
pub struct ExtremaGenExtPs<'a> {
    /// myDone.
    done: bool,
    /// myInit.
    init: bool,
    /// myumin / myusup / myvmin / myvsup.
    umin: f64,
    usup: f64,
    vmin: f64,
    vsup: f64,
    /// myusample / myvsample.
    usample: i32,
    vsample: i32,
    /// mytolu / mytolv.
    tolu: f64,
    tolv: f64,
    /// myPoints (cxx:166), bounds 0..=usample+1 x 0..=vsample+1.
    points: Vec<Vec<ExtremaPOnSurfParams>>,
    /// mySphereArray (cxx:168); index 0..usample*vsample.
    sphere_array: Vec<BndSphere>,
    /// The (NoU, NoV) each mySphereArray entry was built for - the U/V carried
    /// by Bnd_Sphere (cxx:924, Bnd_Sphere.hxx), which the port's BndSphere does
    /// not store.
    sphere_uv: Vec<(i32, i32)>,
    /// !mySphereUBTree.IsNull() (cxx:859). The tree itself is an acceleration
    /// structure; see the module header.
    tree_built: bool,
    /// myF (cxx:169).
    f: ExtremaFuncPsNorm<'a>,
    /// myS (cxx:170).
    s: Option<&'a dyn Surface>,
    /// myFlag / myAlgo.
    flag: ExtremaExtFlag,
    algo: ExtremaExtAlgo,
    /// myUParams / myVParams (cxx:173-174). Vec index i is OCCT index i + 1.
    u_params: Option<Vec<f64>>,
    v_params: Option<Vec<f64>>,
    /// myFacePntParams (cxx:175), bounds 0..=usample x 0..=vsample.
    face_pnt_params: Vec<Vec<ExtremaPOnSurfParams>>,
    /// myUEdgePntParams (cxx:176), OCCT bounds 1..=usample-1 x 1..=vsample.
    u_edge_pnt_params: Vec<Vec<ExtremaPOnSurfParams>>,
    /// myVEdgePntParams (cxx:177), OCCT bounds 1..=usample x 1..=vsample-1.
    v_edge_pnt_params: Vec<Vec<ExtremaPOnSurfParams>>,
    /// myGridParam (cxx:178).
    grid_param: ExtremaPOnSurfParams,
}

impl<'a> Default for ExtremaGenExtPs<'a> {
    fn default() -> Self {
        Self::new()
    }
}

impl<'a> ExtremaGenExtPs<'a> {
    /// Extrema_GenExtPS() (cxx:216-231).
    pub fn new() -> Self {
        ExtremaGenExtPs {
            done: false,
            init: false,
            umin: 0.0,
            usup: 0.0,
            vmin: 0.0,
            vsup: 0.0,
            usample: 0,
            vsample: 0,
            tolu: 0.0,
            tolv: 0.0,
            points: Vec::new(),
            sphere_array: Vec::new(),
            sphere_uv: Vec::new(),
            tree_built: false,
            f: ExtremaFuncPsNorm::default(),
            s: None,
            flag: ExtremaExtFlag::MinMax,
            algo: ExtremaExtAlgo::Grad,
            u_params: None,
            v_params: None,
            face_pnt_params: Vec::new(),
            u_edge_pnt_params: Vec::new(),
            v_edge_pnt_params: Vec::new(),
            grid_param: ExtremaPOnSurfParams::new_default(),
        }
    }

    /// Extrema_GenExtPS(P, S, NbU, NbV, TolU, TolV, F, A) (cxx:239-253):
    /// Initialize then Perform.
    #[allow(clippy::too_many_arguments)]
    pub fn with_surface(
        p: &GpPnt,
        s: &'a dyn Surface,
        nb_u: i32,
        nb_v: i32,
        tol_u: f64,
        tol_v: f64,
        flag: ExtremaExtFlag,
        algo: ExtremaExtAlgo,
    ) -> Self {
        let mut e = ExtremaGenExtPs::new();
        e.flag = flag;
        e.algo = algo;
        e.initialize(s, nb_u, nb_v, tol_u, tol_v);
        e.perform(p);
        e
    }

    /// Extrema_GenExtPS(P, S, NbU, NbV, Umin, Usup, Vmin, Vsup, TolU, TolV, F, A)
    /// (cxx:255-273).
    #[allow(clippy::too_many_arguments)]
    pub fn with_window(
        p: &GpPnt,
        s: &'a dyn Surface,
        nb_u: i32,
        nb_v: i32,
        umin: f64,
        usup: f64,
        vmin: f64,
        vsup: f64,
        tol_u: f64,
        tol_v: f64,
        flag: ExtremaExtFlag,
        algo: ExtremaExtAlgo,
    ) -> Self {
        let mut e = ExtremaGenExtPs::new();
        e.flag = flag;
        e.algo = algo;
        e.initialize_window(s, nb_u, nb_v, umin, usup, vmin, vsup, tol_u, tol_v);
        e.perform(p);
        e
    }

    /// Extrema_GenExtPS::Initialize(S, NbU, NbV, TolU, TolV) (cxx:275-286).
    pub fn initialize(&mut self, s: &'a dyn Surface, nb_u: i32, nb_v: i32, tol_u: f64, tol_v: f64) {
        let (umin, usup) = s.u_range();
        let (vmin, vsup) = s.v_range();
        self.initialize_window(s, nb_u, nb_v, umin, usup, vmin, vsup, tol_u, tol_v);
    }

    /// Extrema_GenExtPS::Initialize(S, NbU, NbV, Umin, Usup, Vmin, Vsup, TolU,
    /// TolV) (cxx:288-319).
    #[allow(clippy::too_many_arguments)]
    pub fn initialize_window(
        &mut self,
        s: &'a dyn Surface,
        nb_u: i32,
        nb_v: i32,
        umin: f64,
        usup: f64,
        vmin: f64,
        vsup: f64,
        tol_u: f64,
        tol_v: f64,
    ) {
        self.s = Some(s);
        self.usample = nb_u;
        self.vsample = nb_v;
        self.tolu = tol_u;
        self.tolv = tol_v;
        self.umin = umin;
        self.usup = usup;
        self.vmin = vmin;
        self.vsup = vsup;

        if self.usample < 2 || self.vsample < 2 {
            // cxx:308-311 throws Standard_OutOfRange.
            panic!("Extrema_GenExtPS::Initialize: nb samples < 2");
        }

        self.f.initialize(s);

        // cxx:315-317: mySphereUBTree.Nullify(); myUParams.Nullify();
        // myVParams.Nullify();
        self.tree_built = false;
        self.sphere_array.clear();
        self.sphere_uv.clear();
        self.u_params = None;
        self.v_params = None;
        self.init = false;
    }

    /// Extrema_GenExtPS::GetGridPoints (cxx:371-454).
    fn get_grid_points(&mut self, s: &dyn Surface) {
        if s.is_offset_surface() {
            // cxx:375-378.
            if let Some(basis) = s.offset_basis_surface() {
                self.get_grid_points(&*basis);
            }
        } else if s.is_bspline_surface() {
            // cxx:380-390.
            if let (Some(uk), Some(vk)) = (s.bspline_surface_uknots(), s.bspline_surface_vknots()) {
                let ud = distinct_knots(uk);
                let vd = distinct_knots(vk);
                self.u_params = fill_params(&ud, s.u_degree(), self.umin, self.usup, self.usample);
                self.v_params = fill_params(&vd, s.v_degree(), self.vmin, self.vsup, self.vsample);
            }
        } else if s.is_bezier_surface() {
            // cxx:392-405.
            let (u1, u2) = s.u_range();
            let (v1, v2) = s.v_range();
            self.u_params = fill_params(&[u1, u2], s.u_degree(), self.umin, self.usup, self.usample);
            self.v_params = fill_params(&[v1, v2], s.v_degree(), self.vmin, self.vsup, self.vsample);
        } else if s.is_surface_of_revolution() || s.is_surface_of_linear_extrusion() {
            // cxx:407-444: parametric points follow the generating curve when it
            // is a B-spline or a Bezier.
            let basis = if s.is_surface_of_revolution() {
                s.revolution_basis_curve()
            } else {
                s.extrusion_basis_curve()
            };
            let Some(curve) = basis else {
                return;
            };
            let mut arr: Option<Vec<f64>> = None;
            let mut degree: i32 = 0;
            if let Some(knots) = curve.bspline_knots() {
                // BasisCurve()->GetType() == GeomAbs_BSplineCurve (cxx:412-420).
                arr = Some(distinct_knots(knots));
                degree = curve.nurbs_degree().unwrap_or(1) as i32;
            }
            if arr.is_none() && curve.bezier_poles().is_some() {
                // BasisCurve()->GetType() == GeomAbs_BezierCurve (cxx:421-431).
                arr = Some(vec![curve.first_parameter(), curve.last_parameter()]);
                degree = curve.nurbs_degree().unwrap_or(1) as i32;
            }
            let Some(arr) = arr else {
                return;
            };
            if s.is_surface_of_revolution() {
                self.v_params = fill_params(&arr, degree, self.vmin, self.vsup, self.vsample);
            } else {
                self.u_params = fill_params(&arr, degree, self.umin, self.usup, self.usample);
            }
        }
        // cxx:445-453: update the number of points in sample.
        if let Some(u) = &self.u_params {
            self.usample = u.len() as i32;
        }
        if let Some(v) = &self.v_params {
            self.vsample = v.len() as i32;
        }
    }

    /// Extrema_GenExtPS::ComputeEdgeParameters (cxx:460-526). OCCT returns a
    /// reference to myGridParam (or to one of the inputs); the port returns the
    /// value, which is what every caller stores.
    fn compute_edge_parameters(
        &mut self,
        is_u_edge: bool,
        param0: &ExtremaPOnSurfParams,
        param1: &ExtremaPOnSurfParams,
        the_point: &GpPnt,
        diff_tol: f64,
    ) -> ExtremaPOnSurfParams {
        let s = self.s.expect("Extrema_GenExtPS: surface not set");
        let sqr_dist01 = param0.value().square_distance(&param1.value());

        if sqr_dist01 <= diff_tol {
            // cxx:469-473: the points are confused.
            return *param0;
        }
        let diff_dist = (param0.get_sqr_distance() - param1.get_sqr_distance()).abs();
        if diff_dist >= sqr_dist01 - diff_tol {
            // cxx:476-491: the shortest distance is one of the nodes.
            return if param0.get_sqr_distance() > param1.get_sqr_distance() {
                *param1
            } else {
                *param0
            };
        }
        // cxx:493-524: the shortest distance is inside the edge.
        let pop = GpVec::from_pnts(&param0.value(), the_point);
        let pop1 = GpVec::from_pnts(&param0.value(), &param1.value());
        let ratio = pop.dot(&pop1) / sqr_dist01;
        let (au0, av0) = param0.parameter();
        let (au1, av1) = param1.parameter();
        let (mut u_par, mut v_par) = (au0, av0);
        if is_u_edge {
            u_par += ratio * (au1 - au0);
        } else {
            v_par += ratio * (av1 - av0);
        }
        self.grid_param.set_parameters(u_par, v_par, s.d0(u_par, v_par));
        let (i0, j0) = param0.get_indices();
        self.grid_param.set_element_type(if is_u_edge {
            ExtremaElementType::UIsoEdge
        } else {
            ExtremaElementType::VIsoEdge
        });
        self.grid_param
            .set_sqr_distance(the_point.square_distance(&self.grid_param.value()));
        self.grid_param.set_indices(i0, j0);
        self.grid_param
    }

    /// Extrema_GenExtPS::BuildGrid (cxx:528-753).
    fn build_grid(&mut self, the_point: &GpPnt) {
        let s = self.s.expect("Extrema_GenExtPS: surface not set");

        if !self.init {
            // cxx:534: build parametric grid for complex geometry.
            self.get_grid_points(s);
            let usample = self.usample as usize;
            let vsample = self.vsample as usize;

            // cxx:537-549: build grid in other cases (U).
            if self.u_params.is_none() {
                let mut pas_u = self.usup - self.umin;
                let mut u0 = pas_u / self.usample as f64 / 100.0;
                pas_u = (pas_u - u0) / (self.usample as f64 - 1.0);
                u0 = u0 / 2.0 + self.umin;
                let mut params = Vec::with_capacity(usample);
                let mut u = u0;
                for _ in 1..=usample {
                    params.push(u);
                    u += pas_u;
                }
                self.u_params = Some(params);
            }
            // cxx:551-564: build grid in other cases (V).
            if self.v_params.is_none() {
                let mut pas_v = self.vsup - self.vmin;
                let mut v0 = pas_v / self.vsample as f64 / 100.0;
                pas_v = (pas_v - v0) / (self.vsample as f64 - 1.0);
                v0 = v0 / 2.0 + self.vmin;
                let mut params = Vec::with_capacity(vsample);
                let mut v = v0;
                for _ in 1..=vsample {
                    params.push(v);
                    v += pas_v;
                }
                self.v_params = Some(params);
            }

            // cxx:566-567: myPoints.Resize(0, usample+1, 0, vsample+1, false).
            let nu = self.usample as usize + 2;
            let nv = self.vsample as usize + 2;
            self.points = vec![vec![ExtremaPOnSurfParams::new_default(); nv]; nu];

            // cxx:569-572: GeomGridEval_Surface::EvaluateGrid over the two
            // parameter arrays - the same points as Value(U_i, V_j).
            let u_params = self.u_params.clone().expect("u params");
            let v_params = self.v_params.clone().expect("v params");
            for no_u in 1..=self.usample as usize {
                for no_v in 1..=self.vsample as usize {
                    let p1 = s.d0(u_params[no_u - 1], v_params[no_v - 1]);
                    let mut a_param =
                        ExtremaPOnSurfParams::new(u_params[no_u - 1], v_params[no_v - 1], p1);
                    a_param.set_element_type(ExtremaElementType::Node);
                    a_param.set_indices(no_u as i32, no_v as i32);
                    self.points[no_u][no_v] = a_param;
                }
            }

            // cxx:591-593.
            let fu = self.usample as usize + 1;
            let fv = self.vsample as usize + 1;
            self.face_pnt_params = vec![vec![ExtremaPOnSurfParams::new_default(); fv]; fu];
            self.u_edge_pnt_params = vec![vec![ExtremaPOnSurfParams::new_default(); fv]; fu];
            self.v_edge_pnt_params = vec![vec![ExtremaPOnSurfParams::new_default(); fv]; fu];

            // cxx:595-607: fill the boundary with negative square distance (used
            // for the maximum search).
            for no_v in 0..=self.vsample as usize {
                self.points[0][no_v].set_sqr_distance(-1.0);
                self.points[self.usample as usize + 1][no_v].set_sqr_distance(-1.0);
            }
            for no_u in 1..=self.usample as usize {
                self.points[no_u][0].set_sqr_distance(-1.0);
                self.points[no_u][self.vsample as usize + 1].set_sqr_distance(-1.0);
            }

            self.init = true;
        }

        let usample = self.usample as usize;
        let vsample = self.vsample as usize;

        // cxx:612-622: step 1, distances to the nodes.
        for no_u in 1..=usample {
            for no_v in 1..=vsample {
                let d = the_point.square_distance(&self.points[no_u][no_v].value());
                self.points[no_u][no_v].set_sqr_distance(d);
            }
        }

        // cxx:624-752: for the minimum, compute distances to the mesh (edges
        // then faces).
        if self.flag == ExtremaExtFlag::Min || self.flag == ExtremaExtFlag::MinMax {
            let diff_tol = self.tolu + self.tolv;

            // cxx:631-660: step 2, distances to the edges.
            for no_u in 1..=usample {
                for no_v in 1..=vsample {
                    let param0 = self.points[no_u][no_v];
                    if no_u < usample {
                        let param1 = self.points[no_u + 1][no_v];
                        let edge =
                            self.compute_edge_parameters(true, &param0, &param1, the_point, diff_tol);
                        self.u_edge_pnt_params[no_u][no_v] = edge;
                    }
                    if no_v < vsample {
                        let param1 = self.points[no_u][no_v + 1];
                        let edge = self.compute_edge_parameters(
                            false, &param0, &param1, the_point, diff_tol,
                        );
                        self.v_edge_pnt_params[no_u][no_v] = edge;
                    }
                }
            }

            // cxx:662-738: step 3, distances to the faces.
            for no_u in 1..usample {
                for no_v in 1..vsample {
                    let ue0 = self.u_edge_pnt_params[no_u][no_v];
                    let ue1 = self.u_edge_pnt_params[no_u][no_v + 1];
                    let ve0 = self.v_edge_pnt_params[no_u][no_v];
                    let ve1 = self.v_edge_pnt_params[no_u + 1][no_v];

                    let mut sqr_dist01 = ue0.value().square_distance(&ue1.value());
                    let mut diff_dist = (ue0.get_sqr_distance() - ue1.get_sqr_distance()).abs();
                    let mut is_out = false;
                    if diff_dist >= sqr_dist01 - diff_tol {
                        // cxx:684-688: projection outside the face.
                        is_out = true;
                    } else {
                        sqr_dist01 = ve0.value().square_distance(&ve1.value());
                        diff_dist = (ve0.get_sqr_distance() - ve1.get_sqr_distance()).abs();
                        if diff_dist >= sqr_dist01 - diff_tol {
                            is_out = true;
                        }
                    }

                    if is_out {
                        // cxx:701-712: closest point on an edge.
                        let ue_min = if ue0.get_sqr_distance() < ue1.get_sqr_distance() {
                            ue0
                        } else {
                            ue1
                        };
                        let ve_min = if ve0.get_sqr_distance() < ve1.get_sqr_distance() {
                            ve0
                        } else {
                            ve1
                        };
                        let e_min = if ue_min.get_sqr_distance() < ve_min.get_sqr_distance() {
                            ue_min
                        } else {
                            ve_min
                        };
                        self.face_pnt_params[no_u][no_v] = e_min;
                    } else {
                        // cxx:713-736: closest point inside the face.
                        let (au0, _av0) = ue0.parameter();
                        let (au1, _av1) = ue1.parameter();
                        let u_par = 0.5 * (au0 + au1);
                        let (_bu0, bv0) = ve0.parameter();
                        let (_bu1, bv1) = ve1.parameter();
                        let v_par = 0.5 * (bv0 + bv1);

                        let mut a_param =
                            ExtremaPOnSurfParams::new(u_par, v_par, s.d0(u_par, v_par));
                        a_param.set_element_type(ExtremaElementType::Face);
                        a_param.set_sqr_distance(the_point.square_distance(&a_param.value()));
                        a_param.set_indices(no_u as i32, no_v as i32);
                        self.face_pnt_params[no_u][no_v] = a_param;
                    }
                }
            }

            // cxx:740-751: fill the boundary with RealLast.
            for no_v in 0..=self.vsample as usize {
                self.face_pnt_params[0][no_v].set_sqr_distance(f64::MAX);
                self.face_pnt_params[self.usample as usize][no_v].set_sqr_distance(f64::MAX);
            }
            for no_u in 1..self.usample as usize {
                self.face_pnt_params[no_u][0].set_sqr_distance(f64::MAX);
                self.face_pnt_params[no_u][self.vsample as usize].set_sqr_distance(f64::MAX);
            }
        }
    }

    /// Extrema_GenExtPS::BuildTree (cxx:856-931). The grid spheres are the
    /// selection input (mySphereArray, cxx:911-929); the NCollection_UBTree is
    /// not reproduced - see the module header.
    fn build_tree(&mut self) {
        if self.tree_built {
            // cxx:858-862.
            return;
        }
        let s = self.s.expect("Extrema_GenExtPS: surface not set");

        // cxx:864-878: raise the sampling for B-spline surfaces.
        if s.is_bspline_surface() {
            let uk = s
                .bspline_surface_uknots()
                .map(distinct_knots)
                .unwrap_or_default();
            let vk = s
                .bspline_surface_vknots()
                .map(distinct_knots)
                .unwrap_or_default();
            let a_u_value = s.u_degree() * uk.len() as i32;
            let a_v_value = s.v_degree() * vk.len() as i32;
            if a_u_value > self.usample {
                self.usample = a_u_value.min(300);
            }
            if a_v_value > self.vsample {
                self.vsample = a_v_value.min(300);
            }
        }
        // cxx:880.
        correct_nb_samples(
            s,
            self.umin,
            self.usup,
            &mut self.usample,
            self.vmin,
            self.vsup,
            &mut self.vsample,
        );

        // cxx:882-903: uniform parametric grid.
        let mut pas_u = self.usup - self.umin;
        let mut pas_v = self.vsup - self.vmin;
        let mut u0 = pas_u / self.usample as f64 / 100.0;
        let mut v0 = pas_v / self.vsample as f64 / 100.0;
        pas_u = (pas_u - u0) / (self.usample as f64 - 1.0);
        pas_v = (pas_v - v0) / (self.vsample as f64 - 1.0);
        u0 = u0 / 2.0 + self.umin;
        v0 = v0 / 2.0 + self.vmin;

        let mut u_params = Vec::with_capacity(self.usample as usize);
        let mut v_params = Vec::with_capacity(self.vsample as usize);
        let mut u = u0;
        for _ in 1..=self.usample {
            u_params.push(u);
            u += pas_u;
        }
        let mut v = v0;
        for _ in 1..=self.vsample {
            v_params.push(v);
            v += pas_v;
        }
        self.u_params = Some(u_params);
        self.v_params = Some(v_params);

        // cxx:905-929: build the sphere array.
        self.sphere_array.clear();
        self.sphere_uv.clear();
        let u_params = self.u_params.clone().unwrap();
        let v_params = self.v_params.clone().unwrap();
        for no_u in 1..=self.usample as usize {
            for no_v in 1..=self.vsample as usize {
                let p1 = s.d0(u_params[no_u - 1], v_params[no_v - 1]);
                // Bnd_Sphere(aP1.XYZ(), 0, NoU, NoV) (cxx:924).
                self.sphere_array.push(BndSphere::from_center_radius(p1, 0.0));
                self.sphere_uv.push((no_u as i32, no_v as i32));
            }
        }
        // cxx:930: aFiller.Fill().
        self.tree_built = true;
    }

    /// Extrema_GenExtPS::FindSolution (cxx:933-952).
    fn find_solution(&mut self, _p: &GpPnt, the_params: &ExtremaPOnSurfParams) {
        let mut tol = MathVector::new(1, 2);
        tol.set_value(1, self.tolu);
        tol.set_value(2, self.tolv);

        let mut uv = MathVector::new(1, 2);
        let (u, v) = the_params.parameter();
        uv.set_value(1, u);
        uv.set_value(2, v);

        let mut uv_inf = MathVector::new(1, 2);
        let mut uv_sup = MathVector::new(1, 2);
        uv_inf.set_value(1, self.umin);
        uv_inf.set_value(2, self.vmin);
        uv_sup.set_value(1, self.usup);
        uv_sup.set_value(2, self.vsup);

        // math_FunctionSetRoot S(myF, Tol); S.Perform(myF, UV, UVinf, UVsup);
        let mut solver = MathFunctionSetRoot::new(&self.f, &tol, 100);
        solver.perform_with_bounds(&mut self.f, &uv, &uv_inf, &uv_sup, false);

        self.done = true;
    }

    /// Extrema_GenExtPS::SetFlag (cxx:954-957).
    pub fn set_flag(&mut self, f: ExtremaExtFlag) {
        self.flag = f;
    }

    /// Extrema_GenExtPS::SetAlgo (cxx:959-966).
    pub fn set_algo(&mut self, a: ExtremaExtAlgo) {
        if self.algo != a {
            self.init = false;
        }
        self.algo = a;
    }

    /// The flag set by [ExtremaGenExtPs::set_flag].
    pub fn flag(&self) -> ExtremaExtFlag {
        self.flag
    }

    /// The algorithm set by [ExtremaGenExtPs::set_algo].
    pub fn algo(&self) -> ExtremaExtAlgo {
        self.algo
    }

    /// Extrema_GenExtPS::Perform (cxx:968-1151).
    pub fn perform(&mut self, p: &GpPnt) {
        self.done = false;
        self.f.set_point(*p);

        if self.algo == ExtremaExtAlgo::Grad {
            self.build_grid(p);
            let usample = self.usample;
            let vsample = self.vsample;

            if self.flag == ExtremaExtFlag::Min || self.flag == ExtremaExtFlag::MinMax {
                // cxx:978-1079 - minimums over the faces.
                for no_u in 1..usample as usize {
                    for no_v in 1..vsample as usize {
                        let a_param = self.face_pnt_params[no_u][no_v];
                        let mut is_min = false;
                        let an_elem_type = a_param.get_element_type();

                        if an_elem_type == ExtremaElementType::Face {
                            is_min = true;
                        } else {
                            let (i_u, i_v) = a_param.get_indices();
                            if an_elem_type == ExtremaElementType::UIsoEdge {
                                is_min = i_v == 1 || i_v == vsample;
                            } else if an_elem_type == ExtremaElementType::VIsoEdge {
                                is_min = i_u == 1 || i_u == usample;
                            } else if an_elem_type == ExtremaElementType::Node {
                                is_min =
                                    (i_u == 1 || i_u == usample) && (i_v == 1 || i_v == vsample);
                            }

                            if !is_min {
                                if an_elem_type == ExtremaElementType::UIsoEdge
                                    || (an_elem_type == ExtremaElementType::Node
                                        && (i_u == 1 || i_u == usample))
                                {
                                    // cxx:1025-1033: check the down face.
                                    let down = self.face_pnt_params[no_u][no_v - 1];
                                    if down.get_element_type() == an_elem_type {
                                        let (i_u2, i_v2) = down.get_indices();
                                        is_min = i_u == i_u2 && i_v == i_v2;
                                    }
                                } else if an_elem_type == ExtremaElementType::VIsoEdge
                                    || (an_elem_type == ExtremaElementType::Node
                                        && (i_v == 1 || i_v == vsample))
                                {
                                    // cxx:1034-1045: check the right face.
                                    let right = self.face_pnt_params[no_u - 1][no_v];
                                    if right.get_element_type() == an_elem_type {
                                        let (i_u2, i_v2) = right.get_indices();
                                        is_min = i_u == i_u2 && i_v == i_v2;
                                    }
                                } else if i_u == no_u as i32 && i_v == no_v as i32 {
                                    // cxx:1046-1069: the lower-left node.
                                    is_min = true;
                                    let others = [
                                        self.face_pnt_params[no_u][no_v - 1],
                                        self.face_pnt_params[no_u - 1][no_v - 1],
                                        self.face_pnt_params[no_u - 1][no_v],
                                    ];
                                    for other in others.iter() {
                                        if !is_min {
                                            break;
                                        }
                                        if other.get_element_type() == ExtremaElementType::Node {
                                            let (i_u2, i_v2) = other.get_indices();
                                            is_min = i_u == i_u2 && i_v == i_v2;
                                        } else {
                                            is_min = false;
                                        }
                                    }
                                }
                            }
                        }

                        if is_min {
                            self.find_solution(p, &a_param);
                        }
                    }
                }
            }

            if self.flag == ExtremaExtFlag::Max || self.flag == ExtremaExtFlag::MinMax {
                // cxx:1081-1111 - maximums over the 3x3 node neighbourhood.
                for no_u in 1..=usample as usize {
                    for no_v in 1..=vsample as usize {
                        let main = self.points[no_u][no_v];
                        let p1 = self.points[no_u - 1][no_v - 1];
                        let p2 = self.points[no_u - 1][no_v];
                        let p3 = self.points[no_u - 1][no_v + 1];
                        let p4 = self.points[no_u][no_v - 1];
                        let p5 = self.points[no_u][no_v + 1];
                        let p6 = self.points[no_u + 1][no_v - 1];
                        let p7 = self.points[no_u + 1][no_v];
                        let p8 = self.points[no_u + 1][no_v + 1];

                        let dist = main.get_sqr_distance();
                        if p1.get_sqr_distance() <= dist
                            && p2.get_sqr_distance() <= dist
                            && p3.get_sqr_distance() <= dist
                            && p4.get_sqr_distance() <= dist
                            && p5.get_sqr_distance() <= dist
                            && p6.get_sqr_distance() <= dist
                            && p7.get_sqr_distance() <= dist
                            && p8.get_sqr_distance() <= dist
                        {
                            self.find_solution(p, &main);
                        }
                    }
                }
            }
        } else {
            // cxx:1113-1150 - the tree algorithm.
            self.build_tree();

            if self.flag == ExtremaExtFlag::Min || self.flag == ExtremaExtFlag::MinMax {
                let idx = self.select_sphere(p, true);
                let (nu, nv) = self.sphere_uv[idx];
                let a_u = self.u_params.as_ref().unwrap()[(nu - 1) as usize];
                let a_v = self.v_params.as_ref().unwrap()[(nv - 1) as usize];
                let s = self.s.expect("Extrema_GenExtPS: surface not set");
                let mut a_params = ExtremaPOnSurfParams::new(a_u, a_v, s.d0(a_u, a_v));
                a_params.set_sqr_distance(p.square_distance(&a_params.value()));
                a_params.set_indices(nu, nv);
                self.find_solution(p, &a_params);
            }
            if self.flag == ExtremaExtFlag::Max || self.flag == ExtremaExtFlag::MinMax {
                let idx = self.select_sphere(p, false);
                let (nu, nv) = self.sphere_uv[idx];
                let a_u = self.u_params.as_ref().unwrap()[(nu - 1) as usize];
                let a_v = self.v_params.as_ref().unwrap()[(nv - 1) as usize];
                let s = self.s.expect("Extrema_GenExtPS: surface not set");
                let mut a_params = ExtremaPOnSurfParams::new(a_u, a_v, s.d0(a_u, a_v));
                a_params.set_sqr_distance(p.square_distance(&a_params.value()));
                a_params.set_indices(nu, nv);
                self.find_solution(p, &a_params);
            }
        }
    }

    /// Bnd_SphereUBTreeSelectorMin/Max::Accept (cxx:104-169) over the whole
    /// mySphereArray: select the sphere minimizing (min == true) or maximizing
    /// Bnd_Sphere::Distance(theCheckPoint). The UBTree only prunes candidates, so
    /// the selected sphere is the same.
    fn select_sphere(&self, p: &GpPnt, min: bool) -> usize {
        let mut best = 0usize;
        let mut best_dist = self.sphere_array[0].distance(p);
        for i in 1..self.sphere_array.len() {
            let d = self.sphere_array[i].distance(p);
            if (min && d < best_dist) || (!min && d > best_dist) {
                best = i;
                best_dist = d;
            }
        }
        best
    }

    /// Extrema_GenExtPS::IsDone (cxx:1155-1158).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// Extrema_GenExtPS::NbExt (cxx:1162-1169). OCCT throws StdFail_NotDone when
    /// the search failed; the port reports the raw count.
    pub fn nb_ext(&self) -> usize {
        self.f.nb_ext()
    }

    /// Extrema_GenExtPS::SquareDistance(N) (cxx:1173-1181).
    pub fn square_distance(&self, n: usize) -> f64 {
        self.f.square_distance(n)
    }

    /// Extrema_GenExtPS::Point(N) (cxx:1185-1193).
    pub fn point(&self, n: usize) -> (f64, f64, GpPnt) {
        let ps = self.f.point(n);
        let (u, v) = ps.parameter();
        (u, v, ps.value())
    }
}
