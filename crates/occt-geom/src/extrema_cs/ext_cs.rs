//! `Extrema_ExtCS` — the curve/surface extrema driver.
//! Source: `Extrema_ExtCS.cxx:37-630`, `.hxx:53-130`.
//!
//! ```text
//! Perform(curve, uinf, usup)                       cxx:108-527
//!   switch curve type / surface type               cxx:124-297
//!     GeomAbs_Line     -> Extrema_ExtElCS (Sphere / Cylinder / Plane)
//!                         Plane not parallel falls through to the generic arm
//!     GeomAbs_Circle   -> Extrema_ExtElCS (Cylinder / Plane / Sphere)
//!     GeomAbs_Hyperbola-> Extrema_ExtElCS (Plane)
//!     anything else    -> generic arm
//!   analytic arm                                    cxx:297-415
//!     IsDone            -> AddSolution per point, else try the curve ends
//!   generic arm                                     cxx:416-527
//!     Extrema_GenExtCS over a NbT=12 / NbU=NbV=10 grid, NbT=13 for a circle
//!     or a B-spline and a [-20, 20] hyperbola window, then the C1 sharp
//!     points projected by Extrema_ExtPS, cut down to min + max
//! ```
//!
//! **UNPORTED** — `cxx:155-197`: for a `GeomAbs_Line` that falls through to the
//! generic arm, OCCT shrinks `[cfirst, clast]` to the span of the curve's line
//! over the surface's bounding box (`Bnd_Box` + `BndLib_AddSurface::Add`,
//! `cxx:165-197`). The `BndLib_AddSurface` port lives in `occt-topo`
//! (`crates/occt-topo/src/geom_bnd_lib_surface3d.rs:769`), which depends on
//! `occt-geom`, so this crate cannot call it. The window therefore stays the
//! caller's `[uinf, usup]` — a superset of OCCT's, so the generic engine
//! searches at least the same interval.

use crate::extrema_pc::{extrema_curve_type, GeomAbsCurveType};
use crate::extrema_surf::ExtPs;
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::Precision;

use super::ext_el_cs::ExtElCS;
use super::gen_ext_cs::GenExtCS;
use super::prelude::*;

/// `GeomAdaptor_Surface::GetType()` split for the `cxx:131-295` switch. The
/// elementary queries come from the same type predicates
/// `extrema_surf::ext_ps_surface_type` uses.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SurKind {
    Plane,
    Cylinder,
    Cone,
    Sphere,
    Torus,
    Bezier,
    BSpline,
    Revolution,
    Extrusion,
    Offset,
    Other,
}

fn surface_kind(s: &dyn Surface) -> SurKind {
    if s.gp_pln().is_some() {
        return SurKind::Plane;
    }
    if s.gp_cylinder().is_some() {
        return SurKind::Cylinder;
    }
    if s.gp_cone().is_some() {
        return SurKind::Cone;
    }
    if s.gp_sphere().is_some() {
        return SurKind::Sphere;
    }
    if s.gp_torus().is_some() {
        return SurKind::Torus;
    }
    if s.is_surface_of_linear_extrusion() {
        return SurKind::Extrusion;
    }
    if s.is_surface_of_revolution() {
        return SurKind::Revolution;
    }
    if s.is_bspline_surface() {
        return SurKind::BSpline;
    }
    if s.is_bezier_surface() {
        return SurKind::Bezier;
    }
    if s.is_offset_surface() {
        return SurKind::Offset;
    }
    SurKind::Other
}

/// `Extrema_ExtCS` (`Extrema_ExtCS.hxx:53-130`).
pub struct ExtCS<'a> {
    /// `myS` (`hxx:113`).
    s: Option<&'a dyn Surface>,
    /// `myDone` (`hxx:114`).
    done: bool,
    /// `myIsPar` (`hxx:115`).
    is_par: bool,
    /// `myuinf` / `myusup` / `myvinf` / `myvsup` (`hxx:116-119`).
    uinf: f64,
    usup: f64,
    vinf: f64,
    vsup: f64,
    /// `mytolC` / `mytolS` (`hxx:120-121`).
    tol_c: f64,
    tol_s: f64,
    /// `myucinf` / `myucsup` (`hxx:122-123`).
    ucinf: f64,
    ucsup: f64,
    /// `myStype` (`hxx:124`).
    stype: SurKind,
    /// `mySqDist` / `myPOnC` / `myPOnS` (`hxx:107-109`).
    sq_dist: Vec<f64>,
    p_on_c: Vec<ExtremaPOnCurv>,
    p_on_s: Vec<ExtremaPOnSurf>,
    /// `myExtElCS` (`hxx:112`).
    ext_el_cs: ExtElCS,
}

impl<'a> ExtCS<'a> {
    /// `Extrema_ExtCS()` (`cxx:37-51`).
    pub fn new() -> Self {
        ExtCS {
            s: None,
            done: false,
            is_par: false,
            uinf: 0.0,
            usup: 0.0,
            vinf: 0.0,
            vsup: 0.0,
            tol_c: 0.0,
            tol_s: 0.0,
            ucinf: 0.0,
            ucsup: 0.0,
            stype: SurKind::Other,
            sq_dist: Vec::new(),
            p_on_c: Vec::new(),
            p_on_s: Vec::new(),
            ext_el_cs: ExtElCS::new(),
        }
    }

    /// `Extrema_ExtCS(C, S, TolC, TolS)` (`cxx:52-61`).
    pub fn with_curve_surface(c: &dyn Curve, s: &'a dyn Surface, tol_c: f64, tol_s: f64) -> Self {
        let mut me = ExtCS::new();
        me.initialize(s, tol_c, tol_s);
        me.perform(c, c.first_parameter(), c.last_parameter());
        me
    }

    /// `Extrema_ExtCS(C, S, UCinf, UCsup, Uinf, Usup, Vinf, Vsup, TolC, TolS)`
    /// (`cxx:63-76`).
    #[allow(clippy::too_many_arguments)]
    pub fn with_curve_surface_range(
        c: &dyn Curve,
        s: &'a dyn Surface,
        uc_inf: f64,
        uc_sup: f64,
        u_inf: f64,
        u_sup: f64,
        v_inf: f64,
        v_sup: f64,
        tol_c: f64,
        tol_s: f64,
    ) -> Self {
        let mut me = ExtCS::new();
        me.initialize_range(s, u_inf, u_sup, v_inf, v_sup, tol_c, tol_s);
        me.perform(c, uc_inf, uc_sup);
        me
    }

    /// `Initialize(S, TolC, TolS)` (`cxx:78-88`).
    pub fn initialize(&mut self, s: &'a dyn Surface, tol_c: f64, tol_s: f64) {
        let (uinf, usup) = s.u_range();
        let (vinf, vsup) = s.v_range();
        self.initialize_range(s, uinf, usup, vinf, vsup, tol_c, tol_s);
    }

    /// `Initialize(S, Uinf, Usup, Vinf, Vsup, TolC, TolS)` (`cxx:90-106`).
    #[allow(clippy::too_many_arguments)]
    pub fn initialize_range(
        &mut self,
        s: &'a dyn Surface,
        u_inf: f64,
        u_sup: f64,
        v_inf: f64,
        v_sup: f64,
        tol_c: f64,
        tol_s: f64,
    ) {
        self.s = Some(s);
        self.is_par = false;
        self.uinf = u_inf;
        self.usup = u_sup;
        self.vinf = v_inf;
        self.vsup = v_sup;
        self.tol_c = tol_c;
        self.tol_s = tol_s;
        self.stype = surface_kind(s);
    }

    /// `Perform(C, Uinf, Usup)` (`cxx:108-527`).
    pub fn perform(&mut self, c: &dyn Curve, uinf: f64, usup: f64) {
        let Some(s) = self.s else {
            return;
        };
        self.ucinf = uinf;
        self.ucsup = usup;
        self.p_on_s.clear();
        self.p_on_c.clear();
        self.sq_dist.clear();

        let mut nb_t = 12usize;
        let mut nb_u = 10usize;
        let mut nb_v = 10usize;
        let ctype = extrema_curve_type(c);

        self.done = false;
        let mut is_compute_analytic = true;

        match ctype {
            GeomAbsCurveType::Line => {
                let line = c.gp_line().expect("GeomAbs_Line without gp_Lin");
                let to_generic = match self.stype {
                    SurKind::Sphere => {
                        let sph = s.gp_sphere().expect("GeomAbs_Sphere");
                        self.ext_el_cs.perform_line_sphere(&line, &sph);
                        false
                    }
                    SurKind::Cylinder => {
                        let cyl = s.gp_cylinder().expect("GeomAbs_Cylinder");
                        self.ext_el_cs.perform_line_cylinder(&line, &cyl);
                        false
                    }
                    SurKind::Plane => {
                        let pln = s.gp_pln().expect("GeomAbs_Plane");
                        self.ext_el_cs.perform_line_plane(&line, &pln);
                        // `cxx:141-146`: a parallel line keeps the analytic
                        // answer, otherwise the whole `cxx:148-259` block runs.
                        !self.ext_el_cs.is_parallel()
                    }
                    _ => true,
                };
                if to_generic {
                    self.perform_line_generic(c, &mut nb_t, &mut nb_u, &mut nb_v);
                    return;
                }
            }
            GeomAbsCurveType::Circle => {
                if self.stype == SurKind::Cylinder {
                    let circ = c.gp_circ().expect("GeomAbs_Circle");
                    let cyl = s.gp_cylinder().expect("GeomAbs_Cylinder");
                    self.ext_el_cs.perform_circle_cylinder(&circ, &cyl);
                } else if self.stype == SurKind::Plane {
                    let circ = c.gp_circ().expect("GeomAbs_Circle");
                    let pln = s.gp_pln().expect("GeomAbs_Plane");
                    self.ext_el_cs.perform_circle_plane(&circ, &pln);
                } else if self.stype == SurKind::Sphere {
                    let circ = c.gp_circ().expect("GeomAbs_Circle");
                    let sph = s.gp_sphere().expect("GeomAbs_Sphere");
                    self.ext_el_cs.perform_circle_sphere(&circ, &sph);
                } else {
                    // `cxx:280-295`: [[fallthrough]] to `default`.
                    is_compute_analytic = false;
                }
            }
            GeomAbsCurveType::Hyperbola => {
                if self.stype == SurKind::Plane {
                    let hypr = c.gp_hyperbola().expect("GeomAbs_Hyperbola");
                    let pln = s.gp_pln().expect("GeomAbs_Plane");
                    self.ext_el_cs.perform_hyperbola_plane(&hypr, &pln);
                } else {
                    is_compute_analytic = false;
                }
            }
            _ => {
                is_compute_analytic = false;
            }
        }

        if is_compute_analytic {
            if self.ext_el_cs.is_done() {
                self.done = true;
                self.is_par = self.ext_el_cs.is_parallel();
                if self.is_par {
                    let d = self.ext_el_cs.square_distance(1);
                    self.sq_dist.push(d);
                } else {
                    let nb_ext = self.ext_el_cs.nb_ext();
                    for i in 1..=nb_ext {
                        let (pc, ps) = self.ext_el_cs.points(i);
                        let u_curve = pc.parameter();
                        let (u, v) = ps.parameter();
                        let d = self.ext_el_cs.square_distance(i);
                        self.add_solution(c, u_curve, u, v, *pc.value(), ps.value(), d);
                    }

                    if self.sq_dist.is_empty() && nb_ext > 0 {
                        // `cxx:322-414`: the analytic extrema lie outside the
                        // curve/surface window — try the curve ends.
                        self.add_curve_ends(c, s);
                    }
                }
                return;
            }
        }

        // `cxx:416-449`: the generic engine.
        let mut ext = GenExtCS::new();
        ext.initialize_default(s, nb_u, nb_v, self.tol_s);
        if ctype == GeomAbsCurveType::Hyperbola {
            let tmin = (-20.0f64).max(c.first_parameter());
            let tmax = 20.0f64.min(c.last_parameter());
            ext.perform_range(c, nb_t, tmin, tmax, self.tol_c);
        } else {
            if (ctype == GeomAbsCurveType::Circle || ctype == GeomAbsCurveType::BSplineCurve)
                && nb_t < 13
            {
                nb_t = 13;
            }
            ext.perform(c, nb_t, self.tol_c);
        }

        self.done = ext.is_done();
        if !self.done {
            return;
        }

        let nb_ext = ext.nb_ext();
        for i in 1..=nb_ext {
            let pc = ext.point_on_curve(i);
            let ps = ext.point_on_surface(i);
            let t = pc.parameter();
            let (u, v) = ps.parameter();
            self.add_solution(c, t, u, v, *pc.value(), ps.value(), ext.square_distance(i));
        }

        // `cxx:451-527`: add the curve's C1 sharp points, then cut the extra
        // ones down to the min and the max.
        let sol_number = self.sq_dist.len();
        let sharp_points = crate::extrema_pc::extrema_intervals(c, Shape::C1 as u8);
        let mut proj = ExtPs::new();
        proj.initialize(
            s,
            s.u_range().0,
            s.u_range().1,
            s.v_range().0,
            s.v_range().1,
            self.tol_s,
            self.tol_s,
        );
        // `cxx:469`: `for (i = 2; i < SharpPoints.Upper(); ++i)` — the interior
        // interval bounds of the 1-based array.
        if sharp_points.len() > 2 {
            for idx in 1..sharp_points.len() - 1 {
                let t = sharp_points[idx];
                let a_pnt = c.value(t);
                proj.perform(&a_pnt);
                if !proj.is_done() {
                    continue;
                }
                let nb_proj = proj.nb_ext();
                let mut jmin = 0usize;
                let mut min_sq_dist = f64::MAX;
                for j in 1..=nb_proj {
                    let d = proj.square_distance(j);
                    if d < min_sq_dist {
                        min_sq_dist = d;
                        jmin = j;
                    }
                }
                if jmin != 0 {
                    let (u, v, p) = proj.point(jmin);
                    self.add_solution(c, t, u, v, a_pnt, p, min_sq_dist);
                }
            }
        }

        let total = self.sq_dist.len();
        let mut imin = sol_number;
        let mut imax = total.wrapping_sub(1);
        for i in sol_number..total {
            if self.sq_dist[i] < self.sq_dist[imin] {
                imin = i;
            }
            if self.sq_dist[i] > self.sq_dist[imax] {
                imax = i;
            }
        }
        if total > sol_number + 2 {
            let min_sq_dist = self.sq_dist[imin];
            let max_sq_dist = self.sq_dist[imax];
            let min_pc = self.p_on_c[imin];
            let max_pc = self.p_on_c[imax];
            let min_ps = self.p_on_s[imin];
            let max_ps = self.p_on_s[imax];
            self.sq_dist.truncate(sol_number);
            self.p_on_c.truncate(sol_number);
            self.p_on_s.truncate(sol_number);
            self.sq_dist.push(min_sq_dist);
            self.p_on_c.push(min_pc);
            self.p_on_s.push(min_ps);
            self.sq_dist.push(max_sq_dist);
            self.p_on_c.push(max_pc);
            self.p_on_s.push(max_ps);
        }
    }

    /// `cxx:148-259`: the `GeomAbs_Line` fall-through into the generic engine.
    fn perform_line_generic(
        &mut self,
        c: &dyn Curve,
        nb_t: &mut usize,
        nb_u: &mut usize,
        nb_v: &mut usize,
    ) {
        let Some(s) = self.s else {
            return;
        };
        // `cxx:158-197` would shrink `[cfirst, clast]` here (UNPORTED, see the
        // module note), so the pair stays immutable.
        let cfirst = self.ucinf;
        let clast = self.ucsup;
        let (ufirst, ulast) = s.u_range();
        let (vfirst, vlast) = s.v_range();

        // `cxx:158-197`: UNPORTED bounding-box clip of `[cfirst, clast]` (see
        // the module note). The test itself is kept so the sampling counts and
        // the degenerate-interval branch stay on the same path as OCCT.
        let _surf_unbounded = Precision::is_infinite(ufirst)
            || Precision::is_infinite(ulast)
            || Precision::is_infinite(vfirst)
            || Precision::is_infinite(vlast);

        if s.is_u_periodic() {
            *nb_u = 13;
        }
        if s.is_v_periodic() {
            *nb_v = 13;
        }

        if clast - cfirst <= CONFUSION {
            // `cxx:206-225`: a degenerate curve interval — project its middle
            // point instead. NOTE: OCCT passes a default-constructed
            // `Extrema_POnCurv PC` to `AddSolution`, whose `Value()` is the
            // origin, so the stored curve point is `(0,0,0)`.
            let a_cpar = 0.5 * (cfirst + clast);
            let a_pm = c.value(a_cpar);
            let an_ext_ps = ExtPs::with_window(
                &a_pm,
                s,
                ufirst,
                ulast,
                vfirst,
                vlast,
                self.tol_s,
                self.tol_s,
            );
            self.done = an_ext_ps.is_done();
            if self.done {
                let nb = an_ext_ps.nb_ext();
                for i in 1..=nb {
                    let (u, v, p) = an_ext_ps.point(i);
                    let d = an_ext_ps.square_distance(i);
                    self.add_solution(c, a_cpar, u, v, GpPnt::zero(), p, d);
                }
            }
            return;
        }

        // `cxx:228-257`.
        let ext = GenExtCS::with_curve_surface_range(
            c,
            s,
            *nb_t,
            *nb_u,
            *nb_v,
            cfirst,
            clast,
            ufirst,
            ulast,
            vfirst,
            vlast,
            self.tol_c,
            self.tol_s,
        );
        self.done = ext.is_done();
        if self.done {
            let nb = ext.nb_ext();
            for i in 1..=nb {
                let pc = ext.point_on_curve(i);
                let ps = ext.point_on_surface(i);
                let t = pc.parameter();
                let (u, v) = ps.parameter();
                self.add_solution(c, t, u, v, *pc.value(), ps.value(), ext.square_distance(i));
            }
        }
    }

    /// `cxx:322-414`: the curve-end fallback of the analytic arm.
    fn add_curve_ends(&mut self, c: &dyn Curve, s: &dyn Surface) {
        let a_t = [self.ucinf, self.ucsup];
        let mut a_p_on_c = [GpPnt::zero(); 2];
        let mut a_p_on_s = [GpPnt::zero(); 2];
        let mut u = [0.0f64; 2];
        let mut v = [0.0f64; 2];
        let mut a_dist = [-1.0f64; 2];

        for i in 0..2 {
            if Precision::is_infinite(a_t[i]) {
                continue;
            }
            a_p_on_c[i] = c.value(a_t[i]);
            match self.stype {
                SurKind::Plane => {
                    let pln = s.gp_pln().expect("GeomAbs_Plane");
                    let (u_, v_) = slib::plane_parameters(&pln.pos, &a_p_on_c[i]);
                    u[i] = u_;
                    v[i] = v_;
                    a_p_on_s[i] = slib::plane_value(&pln, u_, v_);
                }
                SurKind::Sphere => {
                    let sph = s.gp_sphere().expect("GeomAbs_Sphere");
                    let (u_, v_) = slib::sphere_parameters(&sph.position(), &a_p_on_c[i]);
                    u[i] = u_;
                    v[i] = v_;
                    a_p_on_s[i] = slib::sphere_value(&sph, u_, v_);
                }
                SurKind::Cylinder => {
                    let cyl = s.gp_cylinder().expect("GeomAbs_Cylinder");
                    let (u_, v_) = slib::cylinder_parameters(&cyl.position(), &a_p_on_c[i]);
                    u[i] = u_;
                    v[i] = v_;
                    a_p_on_s[i] = slib::cylinder_value(&cyl, u_, v_);
                }
                SurKind::Torus => {
                    let tor = s.gp_torus().expect("GeomAbs_Torus");
                    let (u_, v_) = slib::torus_parameters(
                        &tor.pos,
                        tor.major_radius,
                        tor.minor_radius,
                        &a_p_on_c[i],
                    );
                    u[i] = u_;
                    v[i] = v_;
                    a_p_on_s[i] = slib::torus_value(&tor, u_, v_);
                }
                SurKind::Cone => {
                    let con = s.gp_cone().expect("GeomAbs_Cone");
                    let (u_, v_) = slib::cone_parameters(
                        &con.position(),
                        con.radius(),
                        con.semi_angle(),
                        &a_p_on_c[i],
                    );
                    u[i] = u_;
                    v[i] = v_;
                    a_p_on_s[i] = slib::cone_value(&con, u_, v_);
                }
                _ => continue,
            }
            a_dist[i] = a_p_on_c[i].square_distance(&a_p_on_s[i]);
        }

        // `cxx:370-403`: pick the closer end, or both when they tie.
        let mut b_add = [false; 2];
        if a_dist[0] >= 0.0 && a_dist[1] >= 0.0 {
            let a_diff = a_dist[0] - a_dist[1];
            if a_diff.abs() < CONFUSION {
                b_add[0] = true;
                b_add[1] = true;
            } else if a_diff < 0.0 {
                b_add[0] = true;
            } else {
                b_add[1] = true;
            }
        } else if a_dist[0] >= 0.0 {
            b_add[0] = true;
        } else if a_dist[1] >= 0.0 {
            b_add[1] = true;
        }

        for i in 0..2 {
            if b_add[i] {
                self.add_solution(
                    c,
                    a_t[i],
                    u[i],
                    v[i],
                    a_p_on_c[i],
                    a_p_on_s[i],
                    a_dist[i],
                );
            }
        }
    }

    /// `AddSolution` (`cxx:577-630`): fold the parameters into the windows and
    /// append the solution unless one with the same parameters exists.
    #[allow(clippy::too_many_arguments)]
    fn add_solution(
        &mut self,
        c: &dyn Curve,
        a_t: f64,
        a_u: f64,
        a_v: f64,
        point_on_curve: GpPnt,
        point_on_surf: GpPnt,
        square_dist: f64,
    ) -> bool {
        let Some(s) = self.s else {
            return false;
        };
        let mut t = a_t;
        let mut u = a_u;
        let mut v = a_v;

        if c.is_periodic() {
            t = clib::in_period(t, self.ucinf, self.ucinf + c.period());
        }
        if s.is_u_periodic() {
            u = clib::in_period(u, self.uinf, self.uinf + s.u_period());
        }
        if s.is_v_periodic() {
            v = clib::in_period(v, self.vinf, self.vinf + s.v_period());
        }

        if (self.ucinf - t) <= self.tol_c
            && (t - self.ucsup) <= self.tol_c
            && (self.uinf - u) <= self.tol_s
            && (u - self.usup) <= self.tol_s
            && (self.vinf - v) <= self.tol_s
            && (v - self.vsup) <= self.tol_s
        {
            let mut is_new_solution = true;
            for j in 0..self.sq_dist.len() {
                let tj = self.p_on_c[j].parameter();
                let (uj, vj) = self.p_on_s[j].parameter();
                if (t - tj).abs() <= self.tol_c
                    && (u - uj).abs() <= self.tol_s
                    && (v - vj).abs() <= self.tol_s
                {
                    is_new_solution = false;
                    break;
                }
            }
            if is_new_solution {
                self.sq_dist.push(square_dist);
                self.p_on_c.push(ExtremaPOnCurv::new(t, point_on_curve));
                self.p_on_s.push(ExtremaPOnSurf::new(u, v, point_on_surf));
                return true;
            }
        }
        false
    }

    /// `IsDone()` (`cxx:529-532`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `IsParallel()` (`cxx:534-544`).
    pub fn is_parallel(&self) -> bool {
        assert!(self.is_done(), "Extrema_ExtCS::IsParallel: StdFail_NotDone");
        self.is_par
    }

    /// `NbExt()` (`cxx:556-562`).
    pub fn nb_ext(&self) -> usize {
        assert!(self.is_done(), "Extrema_ExtCS::NbExt: StdFail_NotDone");
        self.sq_dist.len()
    }

    /// `SquareDistance(N)` (`cxx:546-554`), 1-based.
    pub fn square_distance(&self, n: usize) -> f64 {
        assert!(
            n >= 1 && n <= self.nb_ext(),
            "Extrema_ExtCS::SquareDistance: Standard_OutOfRange"
        );
        self.sq_dist[n - 1]
    }

    /// `Points(N, P1, P2)` (`cxx:564-573`), 1-based.
    pub fn points(&self, n: usize) -> (ExtremaPOnCurv, ExtremaPOnSurf) {
        assert!(
            n >= 1 && n <= self.nb_ext(),
            "Extrema_ExtCS::Points: Standard_OutOfRange"
        );
        (self.p_on_c[n - 1], self.p_on_s[n - 1])
    }
}

impl Default for ExtCS<'_> {
    fn default() -> Self {
        ExtCS::new()
    }
}
