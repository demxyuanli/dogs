//! `Extrema_GenExtSS` (TKGeomBase). Source: `Extrema_GenExtSS.cxx` / `.hxx`.
//!
//! The general surface-surface extrema engine that `Extrema_ExtSS` falls back
//! to for every type pair except `Plane` x `Plane` (`Extrema_ExtSS.cxx:131-180`
//! and `cxx:180-221`, both with `NbU = NbV = 20`):
//!
//! 1. sample a `NbU` x `NbV` grid on S2 (`Initialize`, `cxx:169-227`) and on S1
//!    (`Perform`, `cxx:260-306`), inset by `Pas/nsample/100` at each end;
//! 2. scan the `Nsample^2` x `Nsample^2` grid pairs for the minimum and the
//!    maximum squared distance (`cxx:336-370`);
//! 3. `math_BFGS(4)` from the grid minimum over `Extrema_FuncDistSS`
//!    (`cxx:379-387`); on failure restore the seed and solve with
//!    `math_FunctionSetRoot` over `Extrema_FuncExtSS`
//!    (`cxx:390-400`);
//! 4. solve the maximum with `math_FunctionSetRoot` from the grid maximum
//!    (`cxx:406-414`);
//! 5. `myDone = true` (`cxx:417`).
//!
//! `GeomGridEval_Surface::EvaluateGrid` (`GeomGridEval_Surface.cxx:340-364` /
//! `GeomGridEval_OtherSurface.cxx:112-138`) is a batched `S.Value(u, v)` over
//! the two parameter arrays; the port evaluates `d0` directly at the same
//! `(u, v)`.

use occt_math::{BFGS, MathFunctionSetRoot, MathFunctionSetWithDerivatives, MathVector};
use occt_core::gp::GpPnt;
use crate::extrema_surf::ExtremaPOnSurf;
use crate::surface::Surface;

use super::func_ext_ss::{ExtremaFuncDistSS, ExtremaFuncExtSS};

/// `NCollection_HArray2<gp_Pnt>(0, NbU + 1, 0, NbV + 1)`
/// (`Extrema_GenExtSS.cxx:179-180`): a `(NbU + 2)` x `(NbV + 2)` grid indexed
/// `1..=NbU` x `1..=NbV`.
struct PntGrid {
    data: Vec<GpPnt>,
    stride: usize,
}

impl PntGrid {
    fn new(nb_u: usize, nb_v: usize) -> Self {
        PntGrid {
            data: vec![GpPnt::new(0.0, 0.0, 0.0); (nb_u + 2) * (nb_v + 2)],
            stride: nb_v + 2,
        }
    }

    fn set(&mut self, no_u: usize, no_v: usize, p: GpPnt) {
        self.data[no_u * self.stride + no_v] = p;
    }

    fn get(&self, no_u: usize, no_v: usize) -> GpPnt {
        self.data[no_u * self.stride + no_v]
    }
}

/// `Extrema_GenExtSS` (`Extrema_GenExtSS.hxx:31-120`).
pub struct ExtremaGenExtSS<'a> {
    /// `myF`.
    f: ExtremaFuncExtSS<'a>,
    /// `myS2`.
    s2: Option<&'a dyn Surface>,
    u1min: f64,
    u1sup: f64,
    v1min: f64,
    v1sup: f64,
    u2min: f64,
    u2sup: f64,
    v2min: f64,
    v2sup: f64,
    /// `myusample`.
    usample: usize,
    /// `myvsample`.
    vsample: usize,
    /// `mytol1`.
    tol1: f64,
    /// `mytol2`.
    tol2: f64,
    /// `myDone`.
    done: bool,
    /// `myInit` (`Extrema_GenExtSS.hxx:133`). OCCT writes it once in the
    /// default constructor (`cxx:110`) and never reads it; kept so the member
    /// set stays complete.
    #[allow(dead_code)]
    init: bool,
    /// `mypoints1`.
    points1: PntGrid,
    /// `mypoints2`.
    points2: PntGrid,
}

impl<'a> ExtremaGenExtSS<'a> {
    /// `Extrema_GenExtSS()` (`Extrema_GenExtSS.cxx:93-111`).
    pub fn new() -> Self {
        ExtremaGenExtSS {
            f: ExtremaFuncExtSS::new(),
            s2: None,
            u1min: 0.0,
            u1sup: 0.0,
            v1min: 0.0,
            v1sup: 0.0,
            u2min: 0.0,
            u2sup: 0.0,
            v2min: 0.0,
            v2sup: 0.0,
            usample: 0,
            vsample: 0,
            tol1: 0.0,
            tol2: 0.0,
            done: false,
            init: false,
            points1: PntGrid::new(0, 0),
            points2: PntGrid::new(0, 0),
        }
    }

    /// `Extrema_GenExtSS(S1, S2, NbU, NbV, Tol1, Tol2)`
    /// (`Extrema_GenExtSS.cxx:119-129`): `Initialize(S2, NbU, NbV, Tol2)` then
    /// `Perform(S1, Tol1)`.
    pub fn with_surfaces(
        s1: &'a dyn Surface,
        s2: &'a dyn Surface,
        nb_u: usize,
        nb_v: usize,
        tol1: f64,
        tol2: f64,
    ) -> Self {
        let mut me = ExtremaGenExtSS::new();
        let (u2min, u2sup) = super::occt_range(s2.u_range());
        let (v2min, v2sup) = super::occt_range(s2.v_range());
        me.initialize(s2, nb_u, nb_v, u2min, u2sup, v2min, v2sup, tol2);
        me.perform(s1, tol1);
        me
    }

    /// `Initialize(S2, NbU, NbV, Tol2)` (`Extrema_GenExtSS.cxx:156-167`):
    /// take S2's own parameter bounds.
    pub fn initialize_default(&mut self, s2: &'a dyn Surface, nb_u: usize, nb_v: usize, tol2: f64) {
        let (u2min, u2sup) = super::occt_range(s2.u_range());
        let (v2min, v2sup) = super::occt_range(s2.v_range());
        self.initialize(s2, nb_u, nb_v, u2min, u2sup, v2min, v2sup, tol2);
    }

    /// The 8-argument `Initialize` (`Extrema_GenExtSS.cxx:169-227`).
    #[allow(clippy::too_many_arguments)]
    pub fn initialize(
        &mut self,
        s2: &'a dyn Surface,
        nb_u: usize,
        nb_v: usize,
        u2min: f64,
        u2sup: f64,
        v2min: f64,
        v2sup: f64,
        tol2: f64,
    ) {
        self.s2 = Some(s2);
        self.points1 = PntGrid::new(nb_u, nb_v);
        self.points2 = PntGrid::new(nb_u, nb_v);
        self.usample = nb_u;
        self.vsample = nb_v;
        self.u2min = u2min;
        self.u2sup = u2sup;
        self.v2min = v2min;
        self.v2sup = v2sup;
        self.tol2 = tol2;

        // Parametrage de l echantillon sur S2 — cxx:190-197.
        let (u0, pas_u, v0, pas_v) = Self::sample_axes(u2min, u2sup, v2min, v2sup, nb_u, nb_v);

        let mut u_params = Vec::with_capacity(nb_u);
        let mut u = u0;
        for _ in 1..=nb_u {
            u_params.push(u);
            u += pas_u;
        }
        let mut v_params = Vec::with_capacity(nb_v);
        let mut v = v0;
        for _ in 1..=nb_v {
            v_params.push(v);
            v += pas_v;
        }

        for (iu, &uu) in u_params.iter().enumerate() {
            for (iv, &vv) in v_params.iter().enumerate() {
                self.points2.set(iu + 1, iv + 1, s2.d0(uu, vv));
            }
        }
    }

    /// `Perform(S1, Tol1)` (`Extrema_GenExtSS.cxx:232-238`).
    pub fn perform(&mut self, s1: &'a dyn Surface, tol1: f64) {
        let (u1min, u1sup) = super::occt_range(s1.u_range());
        let (v1min, v1sup) = super::occt_range(s1.v_range());
        self.perform_range(s1, u1min, u1sup, v1min, v1sup, tol1);
    }

    /// `Perform(S1, U1min, U1sup, V1min, V1sup, Tol1)`
    /// (`Extrema_GenExtSS.cxx:242-417`).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_range(
        &mut self,
        s1: &'a dyn Surface,
        u1min: f64,
        u1sup: f64,
        v1min: f64,
        v1sup: f64,
        tol1: f64,
    ) {
        let s2 = self.s2.expect("Extrema_GenExtSS: no S2");
        self.f.initialize(s1, s2);
        self.u1min = u1min;
        self.u1sup = u1sup;
        self.v1min = v1min;
        self.v1sup = v1sup;
        self.tol1 = tol1;

        let (usample, vsample) = (self.usample, self.vsample);

        // Parametrage de l echantillon sur S1 — cxx:260-268.
        let (u10, pas_u1, v10, pas_v1) =
            Self::sample_axes(u1min, u1sup, v1min, v1sup, usample, vsample);
        // Parametrage de l echantillon sur S2 — cxx:270-277.
        let (u20, pas_u2, v20, pas_v2) = Self::sample_axes(
            self.u2min, self.u2sup, self.v2min, self.v2sup, usample, vsample,
        );

        // Batch evaluation of S1 over the sample grid — cxx:279-306.
        for no_u in 1..=usample {
            for no_v in 1..=vsample {
                let u = u10 + (no_u - 1) as f64 * pas_u1;
                let v = v10 + (no_v - 1) as f64 * pas_v1;
                self.points1.set(no_u, no_v, s1.d0(u, v));
            }
        }

        // b- Calcul des minima — cxx:309-370.
        let mut distmin = f64::MAX; // RealLast()
        let mut distmax = 0.0;
        let (mut n1u_min, mut n1v_min, mut n2u_min, mut n2v_min) = (0usize, 0usize, 0usize, 0usize);
        let (mut n1u_max, mut n1v_max, mut n2u_max, mut n2v_max) = (0usize, 0usize, 0usize, 0usize);
        for no_u1 in 1..=usample {
            for no_v1 in 1..=vsample {
                let p1 = self.points1.get(no_u1, no_v1);
                for no_u2 in 1..=usample {
                    for no_v2 in 1..=vsample {
                        let p2 = self.points2.get(no_u2, no_v2);
                        let the_dist = p1.square_distance(&p2);
                        if the_dist < distmin {
                            distmin = the_dist;
                            n1u_min = no_u1;
                            n1v_min = no_v1;
                            n2u_min = no_u2;
                            n2v_min = no_v2;
                        }
                        if the_dist > distmax {
                            distmax = the_dist;
                            n1u_max = no_u1;
                            n1v_max = no_v1;
                            n2u_max = no_u2;
                            n2v_max = no_v2;
                        }
                    }
                }
            }
        }

        let mut tol = MathVector::new(1, 4);
        tol.set_value(1, self.tol1);
        tol.set_value(2, self.tol1);
        tol.set_value(3, self.tol2);
        tol.set_value(4, self.tol2);
        let mut uv_inf = MathVector::new(1, 4);
        uv_inf.set_value(1, self.u1min);
        uv_inf.set_value(2, self.v1min);
        uv_inf.set_value(3, self.u2min);
        uv_inf.set_value(4, self.v2min);
        let mut uv_sup = MathVector::new(1, 4);
        uv_sup.set_value(1, self.u1sup);
        uv_sup.set_value(2, self.v1sup);
        uv_sup.set_value(3, self.u2sup);
        uv_sup.set_value(4, self.v2sup);

        // UV = the grid minimum — cxx:372-375.
        let mut uv = MathVector::new(1, 4);
        uv.set_value(1, u10 + (n1u_min - 1) as f64 * pas_u1);
        uv.set_value(2, v10 + (n1v_min - 1) as f64 * pas_v1);
        uv.set_value(3, u20 + (n2u_min - 1) as f64 * pas_u2);
        uv.set_value(4, v20 + (n2v_min - 1) as f64 * pas_v2);

        // Extrema_FuncDistSS aGFSS(S1, *myS2); math_BFGS aBFGSSolver(4);
        // aBFGSSolver.Perform(aGFSS, UV); — cxx:377-379.
        let fdist = ExtremaFuncDistSS::new(s1, s2);
        let bfgs = BFGS::new();
        let bfgs_result = bfgs.minimize(|x| fdist.value(x), |x| fdist.gradient(x), &uv);
        match bfgs_result {
            Ok(location) => {
                // aBFGSSolver.Location(UV); myF.Value(UV, UV);
                // myF.GetStateNumber(); — cxx:381-387.
                uv = location;
                let mut scratch = MathVector::new(1, 4);
                self.f.value(&uv, &mut scratch);
                self.f.get_state_number();
            }
            Err(_) => {
                // If optimum is not computed successfully then compute by old
                // approach. Restore initial point. — cxx:388-400.
                let mut sr1 = MathFunctionSetRoot::new(&self.f, &tol, 100);
                sr1.perform_with_bounds(&mut self.f, &uv, &uv_inf, &uv_sup, false);
            }
        }

        // UV = the grid maximum — cxx:402-408.
        uv.set_value(1, u10 + (n1u_max - 1) as f64 * pas_u1);
        uv.set_value(2, v10 + (n1v_max - 1) as f64 * pas_v1);
        uv.set_value(3, u20 + (n2u_max - 1) as f64 * pas_u2);
        uv.set_value(4, v20 + (n2v_max - 1) as f64 * pas_v2);

        // It is impossible to compute max distance in the same manner, since
        // for the distance functional for max have bad definition. So, for max
        // computation old approach is used. — cxx:410-414.
        let mut sr2 = MathFunctionSetRoot::new(&self.f, &tol, 100);
        sr2.perform_with_bounds(&mut self.f, &uv, &uv_inf, &uv_sup, false);

        self.done = true;
    }

    /// The shared "Parametrage de l echantillon" of `Initialize` /
    /// `Perform` (`Extrema_GenExtSS.cxx:190-197` and `cxx:262-267`): inset the
    /// span by `span/nsample/100`, then distribute `nsample` points over what
    /// remains.
    fn sample_axes(
        u_min: f64,
        u_sup: f64,
        v_min: f64,
        v_sup: f64,
        nb_u: usize,
        nb_v: usize,
    ) -> (f64, f64, f64, f64) {
        let mut pas_u = u_sup - u_min;
        let mut pas_v = v_sup - v_min;
        let mut u0 = pas_u / nb_u as f64 / 100.;
        let mut v0 = pas_v / nb_v as f64 / 100.;
        pas_u = (pas_u - u0) / (nb_u - 1) as f64;
        pas_v = (pas_v - v0) / (nb_v - 1) as f64;
        u0 = u_min + u0 / 2.;
        v0 = v_min + v0 / 2.;
        (u0, pas_u, v0, pas_v)
    }

    /// `IsDone()` (`Extrema_GenExtSS.cxx:421-424`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `NbExt()` (`Extrema_GenExtSS.cxx:428-434`): `myF.NbExt()`.
    pub fn nb_ext(&self) -> usize {
        assert!(self.is_done(), "Extrema_GenExtSS::NbExt: StdFail_NotDone");
        self.f.nb_ext()
    }

    /// `SquareDistance(N)` (`Extrema_GenExtSS.cxx:438-446`).
    pub fn square_distance(&self, n: i32) -> f64 {
        assert!(
            n >= 1 && n as usize <= self.nb_ext(),
            "Extrema_GenExtSS::SquareDistance: Standard_OutOfRange"
        );
        self.f.square_distance(n)
    }

    /// `PointOnS1(N)` (`Extrema_GenExtSS.cxx:450-458`).
    pub fn point_on_s1(&self, n: i32) -> ExtremaPOnSurf {
        assert!(
            n >= 1 && n as usize <= self.nb_ext(),
            "Extrema_GenExtSS::PointOnS1: Standard_OutOfRange"
        );
        self.f.point_on_s1(n)
    }

    /// `PointOnS2(N)` (`Extrema_GenExtSS.cxx:461-470`).
    pub fn point_on_s2(&self, n: i32) -> ExtremaPOnSurf {
        assert!(
            n >= 1 && n as usize <= self.nb_ext(),
            "Extrema_GenExtSS::PointOnS2: Standard_OutOfRange"
        );
        self.f.point_on_s2(n)
    }
}

impl Default for ExtremaGenExtSS<'_> {
    fn default() -> Self {
        ExtremaGenExtSS::new()
    }
}
