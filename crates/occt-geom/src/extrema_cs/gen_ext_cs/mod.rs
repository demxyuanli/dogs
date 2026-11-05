//! `Extrema_GenExtCS` — the general curve/surface extrema engine.
//! Source: `Extrema_GenExtCS.cxx:136-877`, `.hxx:31-163`.
//!
//! ```text
//! Initialize -> sample a NbU x NbV grid on S (cxx:195-287)
//! Perform    -> per curve interval, GlobMin* seeds math_PSO from the grid
//!               plus NbT curve samples, then math_FunctionSetRoot polishes
//!               over Extrema_FuncExtCS (cxx:298-424)
//! ```
//!
//! The search dimension depends on the type pair (`cxx:325-336`):
//!
//! | surface | curve | dim | arm |
//! |---|---|---|---|
//! | quadric (`IsQuadric`, `cxx:44-66`) | any | 1 | [`GenExtCS::glob_min_c_quadric`] |
//! | other | conic (`IsConic`, `cxx:68-91`) | 2 | [`GenExtCS::glob_min_conic_s`] |
//! | other | other | 3 | [`GenExtCS::glob_min_gen_cs`] |
//!
//! `GeomGridEval_Surface::EvaluateGrid` (`GeomGridEval_Surface.cxx:340-364`,
//! `GeomGridEval_OtherSurface.cxx:112-138`) is a batched `S.Value(u, v)` over
//! the two parameter arrays; the port evaluates `d0` at the same `(u, v)`.
//!
//! `GetSurfMaxParamVals` / `GetCurvMaxParamVal` (`cxx:94-133`) cap the
//! infinite parameter spans at `1e10`, except for a hyperbola (or a surface
//! swept from one) which stops at `23` to keep `sinh`/`cosh` finite.
//!
//! The three `GlobMin*` arms live in their own files (`glob_min_gen_cs.rs`,
//! `glob_min_conic_s.rs`, `glob_min_c_quadric.rs`); their bodies are the
//! `Extrema_GenExtCS.cxx` bodies unchanged, only relocated to keep this
//! directory's files under the 800-line gate.

mod glob_min_c_quadric;
mod glob_min_conic_s;
mod glob_min_gen_cs;

use occt_core::precision::Precision;
use occt_math::{
    MathFunctionSetRoot, MathFunctionSetWithDerivatives, MathVector, Pso, PsoParticlesPool,
};

use super::func_ext_cs::FuncExtCS;
use super::glob_opt_func_conic_s::GlobOptFuncConicS;
use super::glob_opt_func_cquadric::GlobOptFuncCQuadric;
use super::glob_opt_func_cs::GlobOptFuncCS;
use super::prelude::*;
use crate::extrema_surf::GenLocateExtPs;

/// `MaxParamVal` (`Extrema_GenExtCS.cxx:39`).
const MAX_PARAM_VAL: f64 = 1.0e10;
/// `aBorderDivisor` (`cxx:40`).
const A_BORDER_DIVISOR: f64 = 1.0e4;
/// `HyperbolaLimit = ln(MaxParamVal)` (`cxx:41`).
const HYPERBOLA_LIMIT: f64 = 23.0;

/// `IsQuadric` (`cxx:44-66`).
fn is_quadric(the_s_type: &dyn Surface) -> bool {
    the_s_type.gp_pln().is_some()
        || the_s_type.gp_cylinder().is_some()
        || the_s_type.gp_cone().is_some()
        || the_s_type.gp_sphere().is_some()
        || the_s_type.gp_torus().is_some()
}

/// `IsConic` (`cxx:68-91`).
fn is_conic(the_c_type: &dyn Curve) -> bool {
    the_c_type.gp_line().is_some()
        || the_c_type.gp_circ().is_some()
        || the_c_type.gp_ellipse().is_some()
        || the_c_type.gp_hyperbola().is_some()
        || the_c_type.gp_parabola().is_some()
}

/// `GetCurvMaxParamVal` (`cxx:94-113`): restrict the maximal parameter of a
/// hyperbola (including one behind a `Geom_OffsetCurve`, optionally through a
/// `Geom_TrimmedCurve`) to avoid FPE.
fn curv_max_param_val(the_c: &dyn Curve) -> f64 {
    if the_c.gp_hyperbola().is_some() {
        return HYPERBOLA_LIMIT;
    }
    if let Some((basis, _)) = the_c.offset_curve() {
        let a_bc = if basis.is_geom_trimmed() {
            basis.trim_basis().unwrap_or(basis)
        } else {
            basis
        };
        if a_bc.gp_hyperbola().is_some() {
            return HYPERBOLA_LIMIT;
        }
    }
    MAX_PARAM_VAL
}

/// `GetSurfMaxParamVals` (`cxx:116-133`).
fn surf_max_param_vals(the_s: &dyn Surface) -> (f64, f64) {
    let mut the_u_max = MAX_PARAM_VAL;
    let mut the_v_max = MAX_PARAM_VAL;

    if the_s.is_surface_of_linear_extrusion() {
        if let Some(bc) = the_s.extrusion_basis_curve() {
            the_u_max = curv_max_param_val(bc.as_ref());
        }
    } else if the_s.is_surface_of_revolution() {
        if let Some(bc) = the_s.revolution_basis_curve() {
            the_v_max = curv_max_param_val(bc.as_ref());
        }
    } else if the_s.is_offset_surface() {
        if let Some(basis) = the_s.offset_basis_surface() {
            let (u, v) = surf_max_param_vals(basis.as_ref());
            the_u_max = u;
            the_v_max = v;
        }
    }
    (the_u_max, the_v_max)
}

/// `RealToInt` (`Standard_Real.hxx`): round half away from zero.
fn real_to_int(v: f64) -> usize {
    let r = if v >= 0.0 { v + 0.5 } else { v - 0.5 };
    r as usize
}

/// `Extrema_GenExtCS` (`Extrema_GenExtCS.hxx:31-163`).
pub struct GenExtCS<'a> {
    /// `myF`.
    f: FuncExtCS<'a>,
    /// `myS`.
    s: Option<&'a dyn Surface>,
    /// `mytmin` / `mytsup`.
    tmin: f64,
    tsup: f64,
    /// `myumin` / `myusup` / `myvmin` / `myvsup`.
    umin: f64,
    usup: f64,
    vmin: f64,
    vsup: f64,
    /// `mytsample` / `myusample` / `myvsample`.
    tsample: usize,
    usample: usize,
    vsample: usize,
    /// `mytol1` / `mytol2`.
    tol1: f64,
    tol2: f64,
    /// `myDone`.
    done: bool,
    /// `mySurfPnts` (`NCollection_HArray2<gp_Pnt>(0, NbU, 0, NbV)`, `cxx:357`):
    /// a `(usample + 1)` x `(vsample + 1)` grid stored row-major.
    surf_pnts: Vec<GpPnt>,
}

impl<'a> GenExtCS<'a> {
    /// `Extrema_GenExtCS()` (`cxx:136-155`).
    pub fn new() -> Self {
        GenExtCS {
            f: FuncExtCS::new(),
            s: None,
            tmin: 0.0,
            tsup: 0.0,
            umin: 0.0,
            usup: 0.0,
            vmin: 0.0,
            vsup: 0.0,
            tsample: 0,
            usample: 0,
            vsample: 0,
            tol1: 0.0,
            tol2: 0.0,
            done: false,
            surf_pnts: Vec::new(),
        }
    }

    /// `Extrema_GenExtCS(C, S, NbT, NbU, NbV, Tol1, Tol2)` (`cxx:160-170`):
    /// `Initialize(S, NbU, NbV, Tol2)` then `Perform(C, NbT, Tol1)`.
    pub fn with_curve_surface(
        c: &'a dyn Curve,
        s: &'a dyn Surface,
        nb_t: usize,
        nb_u: usize,
        nb_v: usize,
        tol1: f64,
        tol2: f64,
    ) -> Self {
        let mut me = GenExtCS::new();
        me.initialize_default(s, nb_u, nb_v, tol2);
        me.perform(c, nb_t, tol1);
        me
    }

    /// `Extrema_GenExtCS(C, S, NbT, NbU, NbV, tmin, tsup, Umin, Usup, Vmin,
    /// Vsup, Tol1, Tol2)` (`cxx:174-191`).
    #[allow(clippy::too_many_arguments)]
    pub fn with_curve_surface_range(
        c: &'a dyn Curve,
        s: &'a dyn Surface,
        nb_t: usize,
        nb_u: usize,
        nb_v: usize,
        tmin: f64,
        tsup: f64,
        u_min: f64,
        u_sup: f64,
        v_min: f64,
        v_sup: f64,
        tol1: f64,
        tol2: f64,
    ) -> Self {
        let mut me = GenExtCS::new();
        me.initialize(s, nb_u, nb_v, u_min, u_sup, v_min, v_sup, tol2);
        me.perform_range(c, nb_t, tmin, tsup, tol1);
        me
    }

    /// `Initialize(S, NbU, NbV, Tol2)` (`cxx:195-205`): take S's own range.
    pub fn initialize_default(&mut self, s: &'a dyn Surface, nb_u: usize, nb_v: usize, tol2: f64) {
        let (umin, usup) = s.u_range();
        let (vmin, vsup) = s.v_range();
        self.initialize(s, nb_u, nb_v, umin, usup, vmin, vsup, tol2);
    }

    /// `Initialize(S, NbU, NbV, Umin, Usup, Vmin, Vsup, Tol2)`
    /// (`cxx:209-287`).
    #[allow(clippy::too_many_arguments)]
    pub fn initialize(
        &mut self,
        s: &'a dyn Surface,
        nb_u: usize,
        nb_v: usize,
        u_min: f64,
        u_sup: f64,
        v_min: f64,
        v_sup: f64,
        tol2: f64,
    ) {
        self.s = Some(s);
        self.usample = nb_u;
        self.vsample = nb_v;
        self.umin = u_min;
        self.usup = u_sup;
        self.vmin = v_min;
        self.vsup = v_sup;
        self.tol2 = tol2;

        let (umaxpar, vmaxpar) = surf_max_param_vals(s);
        if Precision::is_infinite(self.usup) {
            self.usup = umaxpar;
        }
        if Precision::is_infinite(self.umin) {
            self.umin = -umaxpar;
        }
        if Precision::is_infinite(self.vsup) {
            self.vsup = vmaxpar;
        }
        if Precision::is_infinite(self.vmin) {
            self.vmin = -vmaxpar;
        }

        let du = (self.usup - self.umin) / A_BORDER_DIVISOR;
        let dv = (self.vsup - self.vmin) / A_BORDER_DIVISOR;
        let a_min_u = self.umin + du;
        let a_min_v = self.vmin + dv;
        let a_max_u = self.usup - du;
        let a_max_v = self.vsup - dv;

        let a_step_s_u = (a_max_u - a_min_u) / self.usample as f64;
        let a_step_s_v = (a_max_v - a_min_v) / self.vsample as f64;

        self.surf_pnts = vec![GpPnt::new(0.0, 0.0, 0.0); (nb_u + 1) * (nb_v + 1)];
        let stride = self.vsample + 1;
        let mut a_su = a_min_u;
        for a_sui in 0..=self.usample {
            let mut a_sv = a_min_v;
            for a_svi in 0..=self.vsample {
                self.surf_pnts[a_sui * stride + a_svi] = s.d0(a_su, a_sv);
                a_sv += a_step_s_v;
            }
            a_su += a_step_s_u;
        }
    }

    /// `Perform(C, NbT, Tol1)` (`cxx:290-297`): the curve's own range.
    pub fn perform(&mut self, c: &'a dyn Curve, nb_t: usize, tol1: f64) {
        self.tmin = c.first_parameter();
        self.tsup = c.last_parameter();
        self.perform_range(c, nb_t, self.tmin, self.tsup, tol1);
    }

    /// `Perform(C, NbT, tmin, tsup, Tol1)` (`cxx:301-424`).
    pub fn perform_range(&mut self, c: &'a dyn Curve, nb_t: usize, tmin: f64, tsup: f64, tol1: f64) {
        let s = self.s.expect("Extrema_GenExtCS: no surface");
        self.done = false;
        self.f.initialize(c, s);
        self.tmin = tmin;
        self.tsup = tsup;
        self.tol1 = tol1;
        self.tsample = nb_t;

        // `cxx:312`: keep the un-trimmed surface bounds for the PSO borders.
        let trim_usup = self.usup;
        let trim_umin = self.umin;
        let trim_vsup = self.vsup;
        let trim_vmin = self.vmin;

        let a_c_max_val = curv_max_param_val(c);
        if Precision::is_infinite(self.tsup) {
            self.tsup = a_c_max_val;
        }
        if Precision::is_infinite(self.tmin) {
            self.tmin = -a_c_max_val;
        }

        // `cxx:325-336`: reduce the search dimension.
        let mut a_nb_var = 3;
        if is_quadric(s) {
            a_nb_var = 1;
        } else if is_conic(c) {
            a_nb_var = 2;
        }

        let mut tol = MathVector::new(1, 3);
        tol.set_value(1, self.tol1);
        tol.set_value(2, self.tol2);
        tol.set_value(3, self.tol2);

        const A_NB_PARTICLES: usize = 48;

        let mut a_nb_int_c = 1;
        if c.is_closed() || c.is_periodic() {
            let a_period = c.period();
            if c.last_parameter() - c.first_parameter() > 2.0 * a_period / 3.0 {
                a_nb_int_c = 2;
            }
        }

        let d_t = (self.tsup - self.tmin) / a_nb_int_c as f64;
        for an_int in 1..=a_nb_int_c {
            let mut tuvinf = MathVector::new(1, 3);
            tuvinf.set_value(1, self.tmin + (an_int - 1) as f64 * d_t);
            tuvinf.set_value(2, trim_umin);
            tuvinf.set_value(3, trim_vmin);

            let mut tuvsup = MathVector::new(1, 3);
            tuvsup.set_value(1, tuvinf.value(1) + d_t);
            tuvsup.set_value(2, trim_usup);
            tuvsup.set_value(3, trim_vsup);

            let mut tuv = MathVector::new(1, 3);
            if a_nb_var == 3 {
                self.glob_min_gen_cs(c, A_NB_PARTICLES, &tuvinf, &tuvsup, &mut tuv);
            } else if a_nb_var == 2 {
                self.glob_min_conic_s(c, A_NB_PARTICLES, &tuvinf, &tuvsup, &mut tuv);
            } else {
                self.glob_min_c_quadric(c, A_NB_PARTICLES, &tuvinf, &tuvsup, &mut tuv);
            }

            // Find min approximation — cxx:385-386.
            let mut root = MathFunctionSetRoot::new(&self.f, &tol, 100);
            root.perform_with_bounds(&mut self.f, &mut tuv, &tuvinf, &tuvsup, false);
        }

        // Try to remove "false" extrema caused by dividing curve interval —
        // cxx:388-421.
        if a_nb_int_c > 1 && self.f.nb_ext() > 1 {
            let a_sq_dists1 = self.f.square_distances_mut().clone();
            let a_pnts_on_crv1 = self.f.points_on_curve_mut().clone();
            let a_pnts_on_surf1 = self.f.points_on_surf_mut().clone();

            let mut a_min_dist = a_sq_dists1[0];
            for &d in a_sq_dists1.iter().skip(1) {
                if d < a_min_dist {
                    a_min_dist = d;
                }
            }

            self.f.square_distances_mut().clear();
            self.f.points_on_curve_mut().clear();
            self.f.points_on_surf_mut().clear();

            let a_tol = SQUARE_CONFUSION;
            for (i, &d) in a_sq_dists1.iter().enumerate() {
                if (d - a_min_dist).abs() <= a_tol {
                    self.f.square_distances_mut().push(d);
                    self.f.points_on_curve_mut().push(a_pnts_on_crv1[i]);
                    self.f.points_on_surf_mut().push(a_pnts_on_surf1[i]);
                }
            }
        }

        self.done = true;
    }

    /// `IsDone()` (`cxx:827-830`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `NbExt()` (`cxx:834-841`): OCCT throws `StdFail_NotDone`.
    pub fn nb_ext(&self) -> usize {
        assert!(self.is_done(), "Extrema_GenExtCS::NbExt: StdFail_NotDone");
        self.f.nb_ext()
    }

    /// `SquareDistance(N)` (`cxx:845-853`).
    pub fn square_distance(&self, n: usize) -> f64 {
        assert!(
            n >= 1 && n <= self.nb_ext(),
            "Extrema_GenExtCS::SquareDistance: Standard_OutOfRange"
        );
        self.f.square_distance(n)
    }

    /// `PointOnCurve(N)` (`cxx:857-866`).
    pub fn point_on_curve(&self, n: usize) -> ExtremaPOnCurv {
        assert!(
            n >= 1 && n <= self.nb_ext(),
            "Extrema_GenExtCS::PointOnCurve: Standard_OutOfRange"
        );
        self.f.point_on_curve(n)
    }

    /// `PointOnSurface(N)` (`cxx:870-877`).
    pub fn point_on_surface(&self, n: usize) -> ExtremaPOnSurf {
        assert!(
            n >= 1 && n <= self.nb_ext(),
            "Extrema_GenExtCS::PointOnSurface: Standard_OutOfRange"
        );
        self.f.point_on_surface(n)
    }
}

impl Default for GenExtCS<'_> {
    fn default() -> Self {
        GenExtCS::new()
    }
}
