//! `Extrema_ExtPS` — point/surface extrema with the type dispatch, the
//! parameter window and the degenerate-iso sampling rule.
//!
//! Source: `Extrema_ExtPS.cxx` (428 lines) + `Extrema_ExtPS.hxx` (140 lines),
//! `src/ModelingData/TKGeomBase/Extrema/`. Line references below are to the
//! `.cxx`.
//!
//! This module is the faithful **dispatch layer**. Its three engines are:
//!
//! * `Extrema_ExtPElS` — ported ([`super::point_plane_extrema`] and friends,
//!   `Extrema_ExtPElS.cxx`), used for Plane/Cylinder/Cone/Sphere/Torus
//!   (`cxx:276-290`).
//! * `Extrema_ExtPExtS` / `Extrema_ExtPRevS` — **UNPORTED**. OCCT builds them
//!   for SurfaceOfExtrusion / SurfaceOfRevolution (`cxx:292-343`); the port has
//!   neither, so those two types fall through to the general arm (recorded at
//!   the `match` below).
//! * `Extrema_GenExtPS` — **UNPORTED** (`cxx:346`, 1195 lines in 8.0.0:
//!   `GetGridPoints` + `BuildGrid` + `FindSolution` over a `Bnd_Sphere` UBTree
//!   with `math_FunctionSetRoot`). The port substitutes the 24×24 grid +
//!   numeric-Jacobian Newton of `numeric_extrema.rs`; see the module header.
//!
//! What *is* ported here and was missing before: the ±1e10 clamp of the window
//! (`cxx:216-231`), the sampling counts `nbU/nbV` = 44 (B-spline/Bezier) else
//! 32 with 300 on a degenerate iso (`cxx:237-259`), `IsoIsDeg` (`cxx:32-93`),
//! `TreatSolution`'s periodic normalization + window test (`cxx:97-135`), the
//! engine `IsDone()`/`NbExt()` semantics, and `TrimmedSquareDistances`.

use super::prelude::*;
use super::*;

use occt_core::precision;

/// `Extrema_ExtFlag` (`Extrema_ExtFlag.hxx`).
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub enum ExtremaExtFlag {
    Min,
    Max,
    #[default]
    MinMax,
}

/// `Extrema_ExtAlgo` (`Extrema_ExtAlgo.hxx`).
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub enum ExtremaExtAlgo {
    #[default]
    Grad,
    Tree,
}

/// The port's `Adaptor3d_Surface::GetType()` (`GeomAbs_SurfaceType`), limited to
/// the values `Extrema_ExtPS` distinguishes.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum ExtPsSurfaceType {
    Plane,
    Cylinder,
    Cone,
    Sphere,
    Torus,
    SurfaceOfExtrusion,
    SurfaceOfRevolution,
    BSpline,
    Bezier,
    Other,
}

/// `Adaptor3d_Surface::GetType()` taken from the `dyn Surface` type queries
/// (`Surface::gp_pln` etc.). The elementary queries come first because OCCT's
/// adaptor reports the *basis* type for a trimmed surface, and the port's
/// trimmed surfaces forward the same queries.
pub fn ext_ps_surface_type(s: &dyn Surface) -> ExtPsSurfaceType {
    if s.gp_pln().is_some() {
        return ExtPsSurfaceType::Plane;
    }
    if s.gp_cylinder().is_some() {
        return ExtPsSurfaceType::Cylinder;
    }
    if s.gp_cone().is_some() {
        return ExtPsSurfaceType::Cone;
    }
    if s.gp_sphere().is_some() {
        return ExtPsSurfaceType::Sphere;
    }
    if s.gp_torus().is_some() {
        return ExtPsSurfaceType::Torus;
    }
    if s.is_surface_of_linear_extrusion() {
        return ExtPsSurfaceType::SurfaceOfExtrusion;
    }
    if s.is_surface_of_revolution() {
        return ExtPsSurfaceType::SurfaceOfRevolution;
    }
    if s.is_bspline_surface() {
        return ExtPsSurfaceType::BSpline;
    }
    if s.is_bezier_surface() {
        return ExtPsSurfaceType::Bezier;
    }
    ExtPsSurfaceType::Other
}

/// `IsoIsDeg(S, Param, IT, TolMin, TolMax)` (`Extrema_ExtPS.cxx:32-93`).
///
/// Walks the 11 samples of the *other* parameter (`(U2-U1)/10`, aborting on a
/// step below `PConfusion`) and reports whether the first derivative along
/// `Param`'s direction stays inside `[TolMin, TolMax]` everywhere. `iso_u` is
/// `IT == GeomAbs_IsoU` (vary `V`, probe `D1V`); `false` is `GeomAbs_IsoV`
/// (vary `U`, probe `D1U`).
pub fn iso_is_deg(
    s: &dyn Surface,
    param: f64,
    iso_u: bool,
    tol_min: f64,
    tol_max: f64,
) -> bool {
    let (u1, u2) = s.u_range();
    let (v1, v2) = s.v_range();
    let mut along = true;
    if !iso_u {
        if !precision::Precision::is_infinite(u1) && !precision::Precision::is_infinite(u2) {
            let step = (u2 - u1) / 10.0;
            if step < precision::PCONFUSION {
                return false;
            }
            let mut d1_norm_max = 0.0f64;
            let mut t = u1;
            while t <= u2 {
                let (_, d1u, _) = s.d1(t, param);
                d1_norm_max = d1_norm_max.max(d1u.magnitude());
                t += step;
            }
            if d1_norm_max > tol_max || d1_norm_max < tol_min {
                along = false;
            }
        }
    } else if !precision::Precision::is_infinite(v1) && !precision::Precision::is_infinite(v2) {
        let step = (v2 - v1) / 10.0;
        if step < precision::PCONFUSION {
            return false;
        }
        let mut d1_norm_max = 0.0f64;
        let mut t = v1;
        while t <= v2 {
            let (_, _, d1v) = s.d1(param, t);
            d1_norm_max = d1_norm_max.max(d1v.magnitude());
            t += step;
        }
        if d1_norm_max > tol_max || d1_norm_max < tol_min {
            along = false;
        }
    }
    along
}

/// Port of `Extrema_ExtPS` (`Extrema_ExtPS.hxx:30-120`): the window/type layer
/// that `ShapeAnalysis_Surface::ValueOfUV` (`ShapeAnalysis_Surface.cxx:1350`)
/// and the projection call sites drive.
pub struct ExtPs<'a> {
    s: Option<&'a dyn Surface>,
    /// `myuinf`/`myusup`/`myvinf`/`myvsup` — the window, after the infinite
    /// clamp of `Initialize` (`cxx:216-231`).
    uinf: f64,
    usup: f64,
    vinf: f64,
    vsup: f64,
    tolu: f64,
    tolv: f64,
    my_type: ExtPsSurfaceType,
    /// `nbU`/`nbV` handed to `myExtPS.Initialize` (`cxx:237-261`).
    nb_u: i32,
    nb_v: i32,
    b_u_iso_deg: bool,
    b_v_iso_deg: bool,
    done: bool,
    /// `myPoints` — treated solutions only (`TreatSolution`, `cxx:129-134`).
    points: Vec<(f64, f64, GpPnt)>,
    /// `mySqDist`, parallel to `points`.
    sq_dist: Vec<f64>,
    /// `d11..d22`/`P11..P22` (`TrimmedSquareDistances`, `cxx:401-418`). In 8.0.0
    /// nothing ever writes them (they are 0 / default in the constructor,
    /// `cxx:148-152`); the port keeps them to mirror the API.
    d11: f64,
    d12: f64,
    d21: f64,
    d22: f64,
    p11: GpPnt,
    p12: GpPnt,
    p21: GpPnt,
    p22: GpPnt,
    flag: ExtremaExtFlag,
    algo: ExtremaExtAlgo,
}

impl<'a> Default for ExtPs<'a> {
    fn default() -> Self {
        Self::new()
    }
}

impl<'a> ExtPs<'a> {
    /// `Extrema_ExtPS::Extrema_ExtPS()` (`cxx:139-154`).
    pub fn new() -> Self {
        ExtPs {
            s: None,
            uinf: 0.0,
            usup: 0.0,
            vinf: 0.0,
            vsup: 0.0,
            tolu: 0.0,
            tolv: 0.0,
            my_type: ExtPsSurfaceType::Other,
            nb_u: 0,
            nb_v: 0,
            b_u_iso_deg: false,
            b_v_iso_deg: false,
            done: false,
            points: Vec::new(),
            sq_dist: Vec::new(),
            d11: 0.0,
            d12: 0.0,
            d21: 0.0,
            d22: 0.0,
            p11: GpPnt::zero(),
            p12: GpPnt::zero(),
            p21: GpPnt::zero(),
            p22: GpPnt::zero(),
            flag: ExtremaExtFlag::MinMax,
            algo: ExtremaExtAlgo::Grad,
        }
    }

    /// `Extrema_ExtPS(P, S, Uinf, Usup, Vinf, Vsup, TolU, TolV, F, A)`
    /// (`cxx:181-198`): `Initialize` then `Perform`.
    pub fn with_window(
        p: &GpPnt,
        s: &'a dyn Surface,
        uinf: f64,
        usup: f64,
        vinf: f64,
        vsup: f64,
        tolu: f64,
        tolv: f64,
    ) -> Self {
        let mut e = ExtPs::new();
        e.initialize(s, uinf, usup, vinf, vsup, tolu, tolv);
        e.perform(p);
        e
    }

    /// `Extrema_ExtPS(P, S, TolU, TolV, F, A)` (`cxx:158-177`): the natural
    /// parameter range of the surface.
    pub fn with_surface(p: &GpPnt, s: &'a dyn Surface, tolu: f64, tolv: f64) -> Self {
        let (u0, u1) = s.u_range();
        let (v0, v1) = s.v_range();
        ExtPs::with_window(p, s, u0, u1, v0, v1, tolu, tolv)
    }

    /// `Extrema_ExtPS::Initialize` (`cxx:202-265`).
    pub fn initialize(
        &mut self,
        s: &'a dyn Surface,
        uinf: f64,
        usup: f64,
        vinf: f64,
        vsup: f64,
        tolu: f64,
        tolv: f64,
    ) {
        self.s = Some(s);
        self.uinf = uinf;
        self.usup = usup;
        self.vinf = vinf;
        self.vsup = vsup;

        if precision::Precision::is_negative_infinite(self.uinf) {
            self.uinf = -1e10;
        }
        if precision::Precision::is_positive_infinite(self.usup) {
            self.usup = 1e10;
        }
        if precision::Precision::is_negative_infinite(self.vinf) {
            self.vinf = -1e10;
        }
        if precision::Precision::is_positive_infinite(self.vsup) {
            self.vsup = 1e10;
        }

        self.tolu = tolu;
        self.tolv = tolv;
        self.my_type = ext_ps_surface_type(s);

        let is_b = matches!(
            self.my_type,
            ExtPsSurfaceType::BSpline | ExtPsSurfaceType::Bezier
        );
        let mut nb_u = if is_b { 44 } else { 32 };
        let mut nb_v = if is_b { 44 } else { 32 };

        let mut b_u_iso_deg = false;
        let mut b_v_iso_deg = false;
        if self.my_type != ExtPsSurfaceType::Plane {
            b_u_iso_deg = iso_is_deg(s, self.uinf, true, 0.0, 1.0e-9)
                || iso_is_deg(s, self.usup, true, 0.0, 1.0e-9);
            b_v_iso_deg = iso_is_deg(s, self.vinf, false, 0.0, 1.0e-9)
                || iso_is_deg(s, self.vsup, false, 0.0, 1.0e-9);
        }
        if b_u_iso_deg {
            nb_u = 300;
        }
        if b_v_iso_deg {
            nb_v = 300;
        }
        self.nb_u = nb_u;
        self.nb_v = nb_v;
        self.b_u_iso_deg = b_u_iso_deg;
        self.b_v_iso_deg = b_v_iso_deg;

        // UNPORTED (T-67): `myExtPS.Initialize(*myS, nbU, nbV, myuinf, myusup,
        // myvinf, myvsup, mytolu, mytolv)` (`cxx:261`) — `Extrema_GenExtPS` is
        // not ported; the substitute of `numeric_extrema.rs` has no separate initialization,
        // it builds its grid in `Perform`'s general arm. `cxx:263-264`
        // (`myExtPExtS.Nullify(); myExtPRevS.Nullify();`) has no Rust
        // counterpart: the two engines do not exist here.
    }

    /// `Extrema_ExtPS::Perform` (`cxx:269-367`).
    pub fn perform(&mut self, p: &GpPnt) {
        self.points.clear();
        self.sq_dist.clear();
        let Some(s) = self.s else {
            self.done = false;
            return;
        };

        match self.my_type {
            ExtPsSurfaceType::Plane => {
                let pl = s.gp_pln().unwrap();
                let e = point_plane_extrema(&pl, p);
                self.done = true; // `Extrema_ExtPElS` plane arm: 1 solution, done
                self.treat_solution(s, e.u2, e.v2.unwrap_or(0.0), &e.p2, e.distance * e.distance);
                return;
            }
            ExtPsSurfaceType::Cylinder => {
                let cy = s.gp_cylinder().unwrap();
                let sols = point_cylinder_extrema(&cy, p);
                self.done = !sols.is_empty();
                self.treat_all(s, sols);
                return;
            }
            ExtPsSurfaceType::Cone => {
                let co = s.gp_cone().unwrap();
                let sols = point_cone_extrema(&co, p);
                self.done = !sols.is_empty();
                self.treat_all(s, sols);
                return;
            }
            ExtPsSurfaceType::Sphere => {
                let sp = s.gp_sphere().unwrap();
                let sols = point_sphere_extrema(&sp, p);
                self.done = !sols.is_empty();
                self.treat_all(s, sols);
                return;
            }
            ExtPsSurfaceType::Torus => {
                let to = s.gp_torus().unwrap();
                let sols = point_torus_extrema(&to, p);
                self.done = !sols.is_empty();
                self.treat_all(s, sols);
                return;
            }
            ExtPsSurfaceType::SurfaceOfExtrusion | ExtPsSurfaceType::SurfaceOfRevolution => {
                // UNPORTED (T-67): OCCT builds `Extrema_ExtPExtS` / `Extrema_ExtPRevS`
                // here (`cxx:292-343`) — the extrusion/revolution engines are not
                // ported. The general arm below stands in for them, which is how
                // the port behaved before this dispatch existed; it is *not* what
                // OCCT runs for these two types.
                self.perform_general(s, p);
                return;
            }
            ExtPsSurfaceType::BSpline
            | ExtPsSurfaceType::Bezier
            | ExtPsSurfaceType::Other => {
                self.perform_general(s, p);
                return;
            }
        }
    }

    /// `default:` arm of `Extrema_ExtPS::Perform` (`cxx:345-356`) —
    /// `myExtPS.Perform(thePoint)` then `TreatSolution` over its results.
    ///
    /// UNPORTED (T-67): the engine is `Extrema_GenExtPS` (`GenExtPS.cxx:968`),
    /// 1195 lines over `GeomGridEval_Surface` + a `Bnd_Sphere` UBTree +
    /// `math_FunctionSetRoot`; the port substitutes `numeric_extrema.rs`'s 24×24 grid and
    /// numeric-Jacobian Newton. The substitute cannot sample an unbounded range,
    /// so a window left at the `Initialize` ±1e10 clamp is replaced by the
    /// natural-bounds clamp of `analytic_solvers::surf_bound_{u,v}` — a substitute-only
    /// concession, not an OCCT branch.
    fn perform_general(&mut self, s: &dyn Surface, p: &GpPnt) {
        let u_clamped = self.uinf <= -1e9 || self.usup >= 1e9;
        let v_clamped = self.vinf <= -1e9 || self.vsup >= 1e9;
        let sols = if u_clamped && v_clamped {
            point_surface_newton_all(s, p)
        } else {
            let (gu0, gu1) = if u_clamped {
                surf_bound_u(s)
            } else {
                (self.uinf, self.usup)
            };
            let (gv0, gv1) = if v_clamped {
                surf_bound_v(s)
            } else {
                (self.vinf, self.vsup)
            };
            point_surface_newton_all_box(s, p, gu0, gu1, gv0, gv1)
        };
        self.done = !sols.is_empty();
        self.treat_all(s, sols);
    }

    fn treat_all(&mut self, s: &dyn Surface, sols: Vec<ExtremaPair>) {
        for e in sols {
            let v = e.v2.unwrap_or(0.0);
            self.treat_solution(s, e.u2, v, &e.p2, e.distance * e.distance);
        }
    }

    /// `Extrema_ExtPS::TreatSolution` (`cxx:97-135`): normalize periodic
    /// parameters into the window (allowing one period of overshoot for a
    /// trimmed surface, `cxx:105-113`/`:118-127`) and keep the solution only
    /// inside `[uinf - tolu, usup + tolu] × [vinf - tolv, vsup + tolv]`.
    fn treat_solution(&mut self, s: &dyn Surface, mut u: f64, mut v: f64, q: &GpPnt, val: f64) {
        if s.is_u_periodic() {
            u = clib::in_period(u, self.uinf, self.uinf + s.u_period());
            if u > self.usup + self.tolu {
                u -= s.u_period();
            }
            if u < self.uinf - self.tolu {
                u += s.u_period();
            }
        }
        if s.is_v_periodic() {
            v = clib::in_period(v, self.vinf, self.vinf + s.v_period());
            if v > self.vsup + self.tolv {
                v -= s.v_period();
            }
            if v < self.vinf - self.tolv {
                v += s.v_period();
            }
        }
        if (self.uinf - u) <= self.tolu
            && (u - self.usup) <= self.tolu
            && (self.vinf - v) <= self.tolv
            && (v - self.vsup) <= self.tolv
        {
            self.points.push((u, v, *q));
            self.sq_dist.push(val);
        }
    }

    /// `Extrema_ExtPS::IsDone` (`cxx:369-372`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `Extrema_ExtPS::NbExt` (`cxx:383-390`). OCCT throws `StdFail_NotDone`
    /// when the search failed; the port returns 0.
    pub fn nb_ext(&self) -> usize {
        if !self.is_done() {
            return 0;
        }
        self.sq_dist.len()
    }

    /// `Extrema_ExtPS::SquareDistance` (`cxx:374-381`).
    pub fn square_distance(&self, n: usize) -> f64 {
        self.sq_dist[n - 1]
    }

    /// `Extrema_ExtPS::Point` (`cxx:392-399`): `(U, V, P)`.
    pub fn point(&self, n: usize) -> (f64, f64, GpPnt) {
        self.points[n - 1]
    }

    /// `Extrema_ExtPS::TrimmedSquareDistances` (`cxx:401-418`).
    pub fn trimmed_square_distances(&self) -> (f64, f64, f64, f64, GpPnt, GpPnt, GpPnt, GpPnt) {
        (
            self.d11, self.d12, self.d21, self.d22, self.p11, self.p12, self.p21, self.p22,
        )
    }

    /// `Extrema_ExtPS::SetFlag` (`cxx:420-423`).
    ///
    /// UNPORTED (T-67): OCCT forwards the flag to `Extrema_GenExtPS`, which uses
    /// it to keep only minima / only maxima / both; with the gen engine
    /// unported the flag is stored and has no effect (the analytic arms always
    /// return every extremum, exactly as `Extrema_ExtPElS` does).
    pub fn set_flag(&mut self, f: ExtremaExtFlag) {
        self.flag = f;
    }

    /// `Extrema_ExtPS::SetAlgo` (`cxx:425-428`). UNPORTED (T-67) — see
    /// [`ExtPs::set_flag`].
    pub fn set_algo(&mut self, a: ExtremaExtAlgo) {
        self.algo = a;
    }

    /// The flag last passed to [`ExtPs::set_flag`].
    pub fn flag(&self) -> ExtremaExtFlag {
        self.flag
    }

    /// The algorithm last passed to [`ExtPs::set_algo`].
    pub fn algo(&self) -> ExtremaExtAlgo {
        self.algo
    }

    /// `nbU` as handed to `myExtPS.Initialize` (`cxx:237-259`): 44 for
    /// B-spline/Bezier, 32 otherwise, 300 on a degenerate iso.
    pub fn sample_counts(&self) -> (i32, i32) {
        (self.nb_u, self.nb_v)
    }

    /// `bUIsoIsDeg` / `bVIsoIsDeg` (`cxx:242-259`).
    pub fn iso_degenerate(&self) -> (bool, bool) {
        (self.b_u_iso_deg, self.b_v_iso_deg)
    }

    /// The surface type this run dispatched on (`mytype`, `cxx:235`).
    pub fn surface_type(&self) -> ExtPsSurfaceType {
        self.my_type
    }
}
