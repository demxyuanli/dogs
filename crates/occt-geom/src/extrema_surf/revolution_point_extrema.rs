//! `Extrema_ExtPRevS` - point / surface-of-revolution extrema.
//!
//! Source: `Extrema_ExtPRevS.cxx` (599 lines) + `Extrema_ExtPRevS.hxx`
//! (90 lines), `src/ModelingData/TKGeomBase/Extrema/`. Line references below
//! are to the `.cxx`.
//!
//! `Extrema_ExtPS::Perform` builds this engine for a
//! `GeomAbs_SurfaceOfRevolution` and merges its solutions through
//! `TreatSolution` (`Extrema_ExtPS.cxx:319-343`). Two analytic halves are
//! tried: the point is rotated back by the meridian angle `U` derived from
//! its meridian plane (`cxx:314-363`), the point/curve extrema of the
//! generatrix are read (`cxx:369-453`), and the same is repeated after a
//! rotation by `M_PI` (`cxx:455-544`). Everything that is not an analytic
//! basis runs the general `Extrema_GenExtPS` engine instead (`cxx:305-312`).
//!
//! Control flow ported branch for branch: `Initialize` (`cxx:259-296`),
//! `Perform` (`cxx:300-545`), `GetPosition` (`cxx:38-79`),
//! `HasSingularity` (`cxx:83-102`), `PerformExtPElC` (`cxx:106-131`),
//! `IsCaseAnalyticallyComputable` (`cxx:135-164`), `IsOriginalPnt`
//! (`cxx:168-178`), `IsExtremum` (`cxx:182-205`) and the accessors
//! (`cxx:549-599`).
//!
//! The revolution adaptor `Extrema_ExtPS` builds (`cxx:322-323`) is the
//! surface the port already holds, so the basis curve and the axis are read
//! back from it.

use std::sync::Arc;

use super::elementary_curve_extrema::{
    ext_pelc_perform, is_elementary_curve, is_original_pnt, pnt_equal, plane_square_distance,
    vec_angle_with_ref, CurveMap, CurvePlane,
};
use super::prelude::*;
use super::*;

/// `Extrema_ExtPRevS` (`Extrema_ExtPRevS.hxx:28-88`).
pub struct ExtremaExtPRevS<'a> {
    /// `myS` (`hxx:77`).
    my_s: &'a dyn Surface,
    /// `myS->BasisCurve()`, the generatrix the surface stores.
    my_c: Arc<dyn Curve>,
    /// `myS->AxeOfRevolution()`.
    my_axis: GpAx1,
    /// `myPosition` (`hxx:81`) - only `Location()` / `Direction()` are read.
    my_position: CurvePlane,
    /// `myvinf` / `myvsup` / `mytolv` (`hxx:78-80`).
    my_vinf: f64,
    my_vsup: f64,
    my_tolv: f64,
    /// `myExtPS` (`hxx:82`).
    my_ext_ps: ExtremaGenExtPs<'a>,
    /// `myIsAnalyticallyComputable` (`hxx:83`).
    my_is_analytically_computable: bool,
    /// `myDone` (`hxx:84`).
    my_done: bool,
    /// `myNbExt` (`hxx:85`).
    my_nb_ext: usize,
    /// `mySqDist[8]` (`hxx:86`), initialised to `RealLast()`
    /// (`cxx:217-220`).
    my_sq_dist: [f64; 8],
    /// `myPoint[8]` (`hxx:87`).
    my_point: [ExtremaPOnSurf; 8],
}

impl<'a> ExtremaExtPRevS<'a> {
    /// `Extrema_ExtPRevS(P, S, Umin, Usup, Vmin, Vsup, TolU, TolV)`
    /// (`cxx:225-237`): `Initialize` then `Perform`.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        p: &GpPnt,
        s: &'a dyn Surface,
        umin: f64,
        usup: f64,
        vmin: f64,
        vsup: f64,
        tolu: f64,
        tolv: f64,
    ) -> Self {
        let mut e = ExtremaExtPRevS {
            my_s: s,
            my_c: s
                .revolution_basis_curve()
                .expect("Extrema_ExtPRevS: revolution surface without a basis curve"),
            my_axis: s
                .revolution_axis()
                .expect("Extrema_ExtPRevS: revolution surface without an axis"),
            my_position: CurvePlane::default(),
            my_vinf: 0.0,
            my_vsup: 0.0,
            my_tolv: 0.0,
            my_ext_ps: ExtremaGenExtPs::new(),
            my_is_analytically_computable: false,
            my_done: false,
            my_nb_ext: 0,
            my_sq_dist: [f64::MAX; 8],
            my_point: [ExtremaPOnSurf::default(); 8],
        };
        e.initialize(s, umin, usup, vmin, vsup, tolu, tolv);
        e.perform(p);
        e
    }

    /// `Extrema_ExtPRevS::Initialize` (`cxx:259-296`).
    ///
    /// OCCT only recomputes `myPosition` / `myIsAnalyticallyComputable` when
    /// `myS != theS` (`cxx:277-283`); this port builds the engine once per
    /// surface (`Extrema_ExtPS::Perform`, `Extrema_ExtPS.cxx:319-331`), so
    /// the block always runs.
    #[allow(clippy::too_many_arguments)]
    pub fn initialize(
        &mut self,
        s: &'a dyn Surface,
        umin: f64,
        usup: f64,
        vmin: f64,
        vsup: f64,
        tolu: f64,
        tolv: f64,
    ) {
        self.my_vinf = vmin;
        self.my_vsup = vsup;
        self.my_tolv = tolv;

        self.my_done = false;
        self.my_nb_ext = 0;
        self.my_is_analytically_computable = false;

        let an_a_curve = s
            .revolution_basis_curve()
            .expect("Extrema_ExtPRevS: revolution surface without a basis curve");

        self.my_s = s;
        self.my_c = an_a_curve;
        self.my_axis = s
            .revolution_axis()
            .expect("Extrema_ExtPRevS: revolution surface without an axis");
        self.my_position = get_position(s);
        self.my_is_analytically_computable =
            is_case_analytically_computable(self.my_c.as_ref(), &self.my_position, &self.my_axis);

        if !self.my_is_analytically_computable {
            // `int aNbu = 32, aNbv = 32; if (HasSingularity(*theS)) aNbv = 100;`
            // (`cxx:287-292`).
            let a_nbu = 32;
            let a_nbv = if has_singularity(s) { 100 } else { 32 };
            self.my_ext_ps
                .initialize_window(s, a_nbu, a_nbv, umin, usup, vmin, vsup, tolu, tolv);
        }
    }

    /// `Extrema_ExtPRevS::Perform` (`cxx:300-545`).
    pub fn perform(&mut self, p: &GpPnt) {
        self.my_done = false;
        self.my_nb_ext = 0;

        if !self.my_is_analytically_computable {
            self.my_ext_ps.perform(p);
            self.my_done = self.my_ext_ps.is_done();
            self.my_nb_ext = self.my_ext_ps.nb_ext();
            return;
        }

        let ax = self.my_axis;
        let dir = GpVec::from_xyz(ax.direction().xyz());
        let z = self.my_position.dir;
        let o = ax.location();

        // `Pp = P translated by -((O,P).Dir) * Dir; return when P is on the
        // axis of revolution` (`cxx:320-325`).
        let op_dir = GpVec::from_pnts(&o, p).dot(&dir);
        let pp = p.translated_vec(&dir.multiplied_scalar(-op_dir));
        if pnt_equal(&o, &pp) {
            return;
        }

        let mut u;
        let ppp;
        let op_pz = GpVec::from_pnts(&o, &pp).dot(&GpVec::from_xyz(z.xyz()));
        if op_pz.abs() <= RESOLUTION {
            ppp = pp;
            u = 0.0;
        } else {
            ppp = pp.translated_vec(&GpVec::from_xyz(z.xyz()).multiplied_scalar(-op_pz));
            if pnt_equal(&o, &ppp) {
                u = std::f64::consts::PI / 2.0;
            } else {
                // `U = gp_Vec(O, Ppp).AngleWithRef(gp_Vec(O, Pp), Dir);`
                // (`cxx:344`).
                u = vec_angle_with_ref(
                    &GpVec::from_pnts(&o, &ppp),
                    &GpVec::from_pnts(&o, &pp),
                    &dir,
                );
            }
        }

        // `gp_Vec OPpp(O, Ppp), OPq(O, myS->Value(M_PI / 2, 0));`
        // (`cxx:348-359`).
        let op_pp = GpVec::from_pnts(&o, &ppp);
        let mut op_q = GpVec::from_pnts(&o, &self.my_s.d0(std::f64::consts::PI / 2.0, 0.0));
        if u != std::f64::consts::PI / 2.0 {
            if op_q.magnitude() <= RESOLUTION {
                let last = self.my_c.last_parameter();
                op_q = GpVec::from_pnts(&o, &self.my_s.d0(std::f64::consts::PI / 2.0, last / 10.0));
            }
            if vec_angle_with_ref(&op_pp, &op_q, &dir) < 0.0 {
                u += std::f64::consts::PI;
            }
        }

        // `gp_Trsf T; T.SetRotation(Ax, -U); P1 = P.Transformed(T);`
        // (`cxx:361-363`).
        let mut t = GpTrsf::identity();
        let _ = t.set_rotation_ax1(&ax, -u);
        let mut p1 = p.transformed(&t);

        let map = CurveMap::of(&self.my_c);
        let an_ext = ext_pelc_perform(map.analytic.as_ref(), &p1, self.my_tolv);
        if an_ext.done {
            self.my_done = true;
            self.add_solutions(p, &an_ext.sols, &map, u);
        }

        // `T.SetRotation(Ax, M_PI); P1.Transform(T);` (`cxx:455-456`).
        let _ = t.set_rotation_ax1(&ax, std::f64::consts::PI);
        p1.transform(&t);

        let an_ext2 = ext_pelc_perform(map.analytic.as_ref(), &p1, self.my_tolv);
        if an_ext2.done {
            self.my_done = true;
            self.add_solutions(p, &an_ext2.sols, &map, u + std::f64::consts::PI);
        }
    }

    /// One pass of `Extrema_ExtPRevS::Perform`'s solution loop
    /// (`cxx:375-453` for the first half, `cxx:465-543` for the one rotated by
    /// `M_PI`; both use the same `IsVSup` argument per branch): fold `V` into
    /// the window, keep it only when `IsExtremum` accepts it, and store the
    /// deduplicated point.
    fn add_solutions(
        &mut self,
        p: &GpPnt,
        sols: &[super::elementary_curve_extrema::ExtPelcSolution],
        map: &CurveMap,
        u: f64,
    ) {
        // `(anACurve->GetType() == GeomAbs_Circle) || (... GeomAbs_Ellipse)`
        // (`cxx:385`, `:417`).
        let periodic_v = self.my_c.gp_circ().is_some() || self.my_c.gp_ellipse().is_some();

        for sol in sols {
            let mut v = map.to_curve(sol.u);
            let mut e = GpPnt::zero();
            let mut dist2 = 0.0;
            if v > self.my_vsup {
                let mut new_v = self.my_vsup;
                if periodic_v {
                    new_v = clib::in_period(v, self.my_vinf, self.my_vinf + 2.0 * std::f64::consts::PI);
                    if new_v > self.my_vsup {
                        new_v -= 2.0 * std::f64::consts::PI;
                        if new_v + self.my_tolv < self.my_vinf {
                            new_v = self.my_vsup;
                        } else if new_v < self.my_vinf {
                            new_v = self.my_vinf;
                        }
                    }
                }
                v = new_v;
                if !self.is_extremum(u, v, p, true, sol.is_min, &mut e, &mut dist2) {
                    continue;
                }
            } else if v < self.my_vinf {
                let mut new_v = self.my_vinf;
                if periodic_v {
                    // `ElCLib::InPeriod(V, myvsup - 2*M_PI, myvsup);`
                    // (`cxx:419`).
                    new_v = clib::in_period(v, self.my_vsup - 2.0 * std::f64::consts::PI, self.my_vsup);
                    if new_v < self.my_vinf {
                        new_v += 2.0 * std::f64::consts::PI;
                        if new_v - self.my_tolv > self.my_vsup {
                            new_v = self.my_vinf;
                        } else if new_v > self.my_vsup {
                            new_v = self.my_vsup;
                        }
                    }
                }
                v = new_v;
                if !self.is_extremum(u, v, p, false, sol.is_min, &mut e, &mut dist2) {
                    continue;
                }
            } else {
                e = self.my_s.d0(u, v);
                dist2 = p.square_distance(&e);
            }
            if is_original_pnt(&e, &self.my_point, self.my_nb_ext) {
                self.my_point[self.my_nb_ext] = ExtremaPOnSurf::new(u, v, e);
                self.my_sq_dist[self.my_nb_ext] = dist2;
                self.my_nb_ext += 1;
            }
        }
    }

    /// `IsExtremum` (`Extrema_ExtPRevS.cxx:182-205`).
    #[allow(clippy::too_many_arguments)]
    fn is_extremum(
        &self,
        u: f64,
        v: f64,
        p: &GpPnt,
        is_vsup: bool,
        is_min: bool,
        e: &mut GpPnt,
        dist2: &mut f64,
    ) -> bool {
        *e = self.my_s.d0(u, v);
        *dist2 = p.square_distance(e);
        let v_step = if is_vsup { v - 1.0 } else { v + 1.0 };
        let up = p.square_distance(&self.my_s.d0(u + 1.0, v));
        let down = p.square_distance(&self.my_s.d0(u - 1.0, v));
        let side = p.square_distance(&self.my_s.d0(u, v_step));
        if is_min {
            *dist2 < up && *dist2 < down && *dist2 < side
        } else {
            *dist2 > up && *dist2 > down && *dist2 > side
        }
    }

    /// `Extrema_ExtPRevS::IsDone()` (`cxx:549-552`).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// `Extrema_ExtPRevS::NbExt()` (`cxx:556-563`).
    pub fn nb_ext(&self) -> usize {
        if !self.is_done() {
            return 0;
        }
        if self.my_is_analytically_computable {
            self.my_nb_ext
        } else {
            self.my_ext_ps.nb_ext()
        }
    }

    /// `Extrema_ExtPRevS::SquareDistance(N)` (`cxx:567-581`).
    pub fn square_distance(&self, n: usize) -> f64 {
        if self.my_is_analytically_computable {
            self.my_sq_dist[n - 1]
        } else {
            self.my_ext_ps.square_distance(n)
        }
    }

    /// `Extrema_ExtPRevS::Point(N)` (`cxx:585-599`) as `(U, V, P)`.
    pub fn point(&self, n: usize) -> (f64, f64, GpPnt) {
        if self.my_is_analytically_computable {
            let (u, v) = self.my_point[n - 1].parameter();
            (u, v, self.my_point[n - 1].value())
        } else {
            self.my_ext_ps.point(n)
        }
    }
}

/// `GetPosition(S)` (`Extrema_ExtPRevS.cxx:38-79`): the supporting plane of
/// the generatrix, as `(location, normal)`.
fn get_position(s: &dyn Surface) -> CurvePlane {
    let c = s
        .revolution_basis_curve()
        .expect("Extrema_ExtPRevS: revolution surface without a basis curve");
    let ax = s
        .revolution_axis()
        .expect("Extrema_ExtPRevS: revolution surface without an axis");

    if let Some(l) = c.gp_line() {
        let mut n = *ax.direction();
        if n.is_parallel_tol(&l.direction(), ANGULAR) {
            let mut oo = GpVec::from_pnts(&l.location(), &ax.location());
            if oo.magnitude() <= RESOLUTION {
                oo = GpVec::from_pnts(&l.location(), &clib::line_value(&l, 100.0));
                if let Ok(d) = GpDir::from_vec(&oo) {
                    if n.is_parallel_tol(&d, ANGULAR) {
                        return CurvePlane::default(); // line and axis coincide
                    }
                }
            }
            // `N ^= OO;` (`cxx:60`) - gp_Dir Cross with the normalised OO.
            let cross = GpDir::from_vec(&oo).and_then(|d| n.crossed(&d));
            n = cross.unwrap_or(n);
        } else {
            // `N ^= L.Direction();` (`cxx:64`).
            n = n.crossed(&l.direction()).unwrap_or(n);
        }
        return CurvePlane {
            loc: l.location(),
            dir: n,
        };
    }
    if let Some(gc) = c.gp_circ() {
        let p = gc.position();
        return CurvePlane {
            loc: p.location(),
            dir: p.direction(),
        };
    }
    if let Some(e) = c.gp_ellipse() {
        let p = *e.position();
        return CurvePlane {
            loc: p.location(),
            dir: p.direction(),
        };
    }
    if let Some(h) = c.gp_hyperbola() {
        let p = *h.position();
        return CurvePlane {
            loc: p.location(),
            dir: p.direction(),
        };
    }
    if let Some(pa) = c.gp_parabola() {
        let p = *pa.position();
        return CurvePlane {
            loc: p.location(),
            dir: p.direction(),
        };
    }
    CurvePlane::default()
}

/// `HasSingularity` (`Extrema_ExtPRevS.cxx:83-102`).
fn has_singularity(s: &dyn Surface) -> bool {
    let c = s
        .revolution_basis_curve()
        .expect("Extrema_ExtPRevS: revolution surface without a basis curve");
    let axis = s
        .revolution_axis()
        .expect("Extrema_ExtPRevS: revolution surface without an axis");
    let l = GpLin::from_pnt_dir(*axis.location(), *axis.direction());
    let p1 = c.d0(c.first_parameter());
    if l.square_distance(&p1) < SQUARE_CONFUSION {
        return true;
    }
    let p2 = c.d0(c.last_parameter());
    l.square_distance(&p2) < SQUARE_CONFUSION
}

/// `IsCaseAnalyticallyComputable` (`Extrema_ExtPRevS.cxx:135-164`): the five
/// elementary types, and the axis of revolution inside the generatrix's plane.
fn is_case_analytically_computable(c: &dyn Curve, plane: &CurvePlane, axis: &GpAx1) -> bool {
    if !is_elementary_curve(c) {
        return false;
    }
    // `double dist = 100., aThreshold = Angular()^2 * dist^2; p2 =
    // Axe.Location + dist * Axe.Direction;` (`cxx:154-156`).
    let dist = 100.0;
    let a_threshold = ANGULAR * ANGULAR * dist * dist;
    let p1 = axis.location();
    let p2 = GpPnt::from_xyz(&axis.location().xyz().added(&axis.direction().xyz().multiplied(dist)));
    plane_square_distance(&plane.loc, &plane.dir, &p1) < a_threshold
        && plane_square_distance(&plane.loc, &plane.dir, &p2) < a_threshold
}
