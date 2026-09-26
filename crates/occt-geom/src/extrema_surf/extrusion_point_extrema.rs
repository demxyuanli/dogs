//! `Extrema_ExtPExtS` - point / surface-of-linear-extrusion extrema.
//!
//! Source: `Extrema_ExtPExtS.cxx` (630 lines) + `Extrema_ExtPExtS.hxx`
//! (104 lines), `src/ModelingData/TKGeomBase/Extrema/`. Line references below
//! are to the `.cxx`.
//!
//! `Extrema_ExtPS::Perform` builds this engine for a
//! `GeomAbs_SurfaceOfExtrusion` and merges its solutions through
//! `TreatSolution` (`Extrema_ExtPS.cxx:292-317`). The engine reduces the
//! problem to a point/elementary-curve extrema on the basis curve whenever the
//! basis is one of the five analytic types supported by `Extrema_ExtPElC` and
//! its supporting plane is not parallel to the extrusion direction
//! (`IsCaseAnalyticallyComputable`, `cxx:567-585`); the remaining cases run
//! the general `Extrema_GenExtPS` engine on the surface itself
//! (`cxx:276-284`).
//!
//! Control flow ported branch for branch: `Initialize` (`cxx:230-264`),
//! `Perform` (`cxx:268-438`), `MakePreciser` (`cxx:81-145`),
//! `GetPosition` (`cxx:509-534`), `PerformExtPElC` (`cxx:538-563`),
//! `IsCaseAnalyticallyComputable` (`cxx:567-585`), `GetValue`
//! (`cxx:589-606`) and the accessors (`cxx:442-505`).
//!
//! Two implementation devices are not control-flow differences:
//!
//! * `GetValue(U, C)` evaluates the analytic primitive (`cxx:589-606`); the
//!   port evaluates the curve object the surface itself stores, which is the
//!   same point for every untrimmed basis and, for this port's trimmed view, is
//!   the parameterisation the surface uses.
//! * The extrusion adaptor `Extrema_ExtPS` builds
//!   (`cxx:295-297`) is the surface the port already holds, so the basis
//!   curve and direction are read back from it.

use std::sync::Arc;

use super::elementary_curve_extrema::{
    ext_pelc_perform, is_elementary_curve, is_original_pnt, project_pnt, CurveMap, CurvePlane,
};
use super::prelude::*;
use super::*;

use occt_core::elib::clib2d;
use occt_math::{MathFunctionSetRoot, MathFunctionSetWithDerivatives, MathVector};

/// `Extrema_ExtPExtS` (`Extrema_ExtPExtS.hxx:30-102`).
pub struct ExtremaExtPExtS<'a> {
    /// `myuinf` / `myusup` / `mytolu` (`hxx:85-87`).
    my_uinf: f64,
    my_usup: f64,
    my_tolu: f64,
    /// `myvinf` / `myvsup` / `mytolv` (`hxx:88-90`).
    my_vinf: f64,
    my_vsup: f64,
    my_tolv: f64,
    /// `myF` - `Extrema_FuncPSNorm` on the extrusion surface (`hxx:91`).
    my_f: ExtremaFuncPsNorm<'a>,
    /// `myC` (`hxx:92`) - the basis curve exactly as the surface stores it.
    my_c: Arc<dyn Curve>,
    /// `myS` (`hxx:93`).
    my_s: &'a dyn Surface,
    /// `myDirection` (`hxx:94`).
    my_direction: GpDir,
    /// `myPosition` (`hxx:95`) - only `Location()` / `Direction()` are read.
    my_position: CurvePlane,
    /// `myExtPS` (`hxx:96`).
    my_ext_ps: ExtremaGenExtPs<'a>,
    /// `myIsAnalyticallyComputable` (`hxx:97`).
    my_is_analytically_computable: bool,
    /// `myDone` (`hxx:98`).
    my_done: bool,
    /// `myNbExt` (`hxx:99`) - the 0-based count into `my_point` /
    /// `my_sq_dist` (`cxx:326-328`).
    my_nb_ext: usize,
    /// `mySqDist[4]` (`hxx:100`), initialised to `RealLast()`
    /// (`cxx:160-163`).
    my_sq_dist: [f64; 4],
    /// `myPoint[4]` (`hxx:101`).
    my_point: [ExtremaPOnSurf; 4],
}

impl<'a> ExtremaExtPExtS<'a> {
    /// `Extrema_ExtPExtS(P, S, Umin, Usup, Vmin, Vsup, TolU, TolV)`
    /// (`cxx:168-194`): `Initialize` then `Perform`.
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
        let mut e = ExtremaExtPExtS {
            my_uinf: 0.0,
            my_usup: 0.0,
            my_tolu: 0.0,
            my_vinf: 0.0,
            my_vsup: 0.0,
            my_tolv: 0.0,
            my_f: ExtremaFuncPsNorm::default(),
            my_c: s
                .extrusion_basis_curve()
                .expect("Extrema_ExtPExtS: extrusion surface without a basis curve"),
            my_s: s,
            my_direction: s
                .extrusion_direction()
                .expect("Extrema_ExtPExtS: extrusion surface without a direction"),
            my_position: CurvePlane::default(),
            my_ext_ps: ExtremaGenExtPs::new(),
            my_is_analytically_computable: false,
            my_done: false,
            my_nb_ext: 0,
            my_sq_dist: [f64::MAX; 4],
            my_point: [ExtremaPOnSurf::default(); 4],
        };
        e.initialize(s, umin, usup, vmin, vsup, tolu, tolv);
        e.perform(p);
        e
    }

    /// `Extrema_ExtPExtS::Initialize` (`cxx:230-264`).
    #[allow(clippy::too_many_arguments)]
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
        self.my_uinf = uinf;
        self.my_usup = usup;
        self.my_tolu = tolu;

        self.my_vinf = vinf;
        self.my_vsup = vsup;
        self.my_tolv = tolv;

        self.my_is_analytically_computable = false;
        self.my_done = false;
        self.my_nb_ext = 0;

        let an_a_curve = s
            .extrusion_basis_curve()
            .expect("Extrema_ExtPExtS: extrusion surface without a basis curve");

        // `myF.Initialize(*theS);` (`cxx:252`).
        self.my_f.initialize(s);
        self.my_c = an_a_curve;
        self.my_s = s;
        self.my_position = get_position(self.my_c.as_ref());
        self.my_direction = s
            .extrusion_direction()
            .expect("Extrema_ExtPExtS: extrusion surface without a direction");
        self.my_is_analytically_computable =
            is_case_analytically_computable(self.my_c.as_ref(), &self.my_position, &self.my_direction);

        if !self.my_is_analytically_computable {
            // `myExtPS.Initialize(*theS, 32, 32, ...)` (`cxx:262`).
            self.my_ext_ps
                .initialize_window(s, 32, 32, uinf, usup, vinf, vsup, tolu, tolv);
        }
    }

    /// `Extrema_ExtPExtS::Perform` (`cxx:268-438`).
    pub fn perform(&mut self, p: &GpPnt) {
        /// `const int NbExtMax = 4` - dimension of `myPoint` / `mySqDist`
        /// (`cxx:270`).
        const NB_EXT_MAX: usize = 4;

        self.my_done = false;
        self.my_nb_ext = 0;

        if !self.my_is_analytically_computable {
            self.my_ext_ps.perform(p);
            self.my_done = self.my_ext_ps.is_done();
            self.my_nb_ext = self.my_ext_ps.nb_ext();
            return;
        }

        let dir = GpVec::from_xyz(self.my_direction.xyz());
        let pp = project_pnt(&self.my_position.loc, &self.my_position.dir, &dir, p);
        // `Extrema_ExtPElC anExt; PerformExtPElC(anExt, Pp, myC, mytolu);`
        // (`cxx:287-288`).
        let map = CurveMap::of(&self.my_c);
        let an_ext = ext_pelc_perform(map.analytic.as_ref(), &pp, self.my_tolu);
        if !an_ext.done {
            return;
        }

        // `bool isSimpleCase = myDirection.IsParallel(myPosition.Direction(),
        // Precision::Angular());` (`cxx:296`).
        let is_simple_case = dir.is_parallel_ang(&GpVec::from_xyz(self.my_position.dir.xyz()), ANGULAR);

        let mut uv = MathVector::new(1, 2);
        let mut tol = MathVector::new(1, 2);
        let mut uv_inf = MathVector::new(1, 2);
        let mut uv_sup = MathVector::new(1, 2);
        tol.set_value(1, self.my_tolu);
        tol.set_value(2, self.my_tolv);
        uv_inf.set_value(1, self.my_uinf);
        uv_inf.set_value(2, self.my_vinf);
        uv_sup.set_value(1, self.my_usup);
        uv_sup.set_value(2, self.my_vsup);

        // `math_FunctionSetRoot aFSR(myF, Tol);` (`cxx:345`).
        let mut a_fsr = MathFunctionSetRoot::new(&self.my_f, &tol, 100);

        for sol in &an_ext.sols {
            let mut u = map.to_curve(sol.u);
            // modified by jgv, 23.12.2008 for OCC17194 (`cxx:310-316`).
            if self.my_c.is_periodic() {
                let mut u2 = u;
                clib2d::adjust_periodic(
                    self.my_uinf,
                    self.my_uinf + 2.0 * std::f64::consts::PI,
                    PCONFUSION,
                    &mut u,
                    &mut u2,
                );
            }
            let e = sol.point;
            let pe = project_pnt(p, &self.my_direction, &dir, &e);

            if is_simple_case {
                // `V = gp_Vec(E, Pe) * gp_Vec(myDirection);` (`cxx:322`).
                let v = GpVec::from_pnts(&e, &pe).dot(&dir);
                self.my_point[self.my_nb_ext] = ExtremaPOnSurf::new(u, v, pe);
                self.my_sq_dist[self.my_nb_ext] = sol.sq_dist;
                self.my_nb_ext += 1;
                if self.my_nb_ext == NB_EXT_MAX {
                    break;
                }
            } else {
                self.my_f.set_point(*p);
                // `isMin = anExt.IsMin(i);` (`cxx:338`).
                let is_min = sol.is_min;

                self.make_preciser(&mut u, p, is_min, p, &self.my_direction);
                let e2 = self.my_c.d0(u);
                let pe2 = project_pnt(p, &self.my_direction, &dir, &e2);
                let v = GpVec::from_pnts(&e2, &pe2).dot(&dir);
                uv.set_value(1, u);
                uv.set_value(2, v);
                a_fsr.perform_with_bounds(&mut self.my_f, &uv, &uv_inf, &uv_sup, false);

                for k in 1..=self.my_f.nb_ext() {
                    let pt = self.my_f.point(k);
                    if is_original_pnt(&pt.value(), &self.my_point, self.my_nb_ext) {
                        self.my_point[self.my_nb_ext] = pt;
                        self.my_sq_dist[self.my_nb_ext] = self.my_f.square_distance(k);
                        self.my_nb_ext += 1;
                        if self.my_nb_ext == NB_EXT_MAX {
                            break;
                        }
                    }
                }
                if self.my_nb_ext == NB_EXT_MAX {
                    break;
                }

                // try symmetric point (`cxx:370-379`).
                self.my_f.set_point(*p); // To clear previous solutions
                u *= -1.0;
                self.make_preciser(&mut u, p, is_min, p, &self.my_direction);
                let e3 = self.my_c.d0(u);
                let pe3 = project_pnt(p, &self.my_direction, &dir, &e3);
                let v = GpVec::from_pnts(&e3, &pe3).dot(&dir);
                uv.set_value(1, u);
                uv.set_value(2, v);

                a_fsr.perform_with_bounds(&mut self.my_f, &uv, &uv_inf, &uv_sup, false);

                for k in 1..=self.my_f.nb_ext() {
                    if self.my_f.square_distance(k) > CONFUSION * CONFUSION {
                        // Additional checking solution: FSR sometimes is wrong
                        // when starting point is far from solution
                        // (`cxx:383-415`).
                        let dist = self.my_f.square_distance(k).sqrt();
                        let mut vals = MathVector::new(1, 2);
                        let pon_s = self.my_f.point(k);
                        let (u2, v2) = pon_s.parameter();
                        uv.set_value(1, u2);
                        uv.set_value(2, v2);
                        let _ = MathFunctionSetWithDerivatives::value(&mut self.my_f, &uv, &mut vals);
                        let (_, du, dv) = self.my_s.d1(u2, v2);
                        let mdu = du.magnitude();
                        let mdv = dv.magnitude();
                        let du_abs = vals.value(1).abs();
                        let dv_abs = vals.value(2).abs();
                        if mdu > PCONFUSION && du_abs / dist / mdu > PCONFUSION {
                            continue;
                        }
                        if mdv > PCONFUSION && dv_abs / dist / mdv > PCONFUSION {
                            continue;
                        }
                    }
                    let pt = self.my_f.point(k);
                    if is_original_pnt(&pt.value(), &self.my_point, self.my_nb_ext) {
                        self.my_point[self.my_nb_ext] = pt;
                        self.my_sq_dist[self.my_nb_ext] = self.my_f.square_distance(k);
                        self.my_nb_ext += 1;
                        if self.my_nb_ext == NB_EXT_MAX {
                            break;
                        }
                    }
                }
                if self.my_nb_ext == NB_EXT_MAX {
                    break;
                }
            }
        }
        self.my_done = true;
    }

    /// `Extrema_ExtPExtS::MakePreciser` (`cxx:81-145`). `OrtogSection` is
    /// `gp_Ax2(P, myDirection)` (`cxx:294`).
    ///
    /// The `pnext = pprev` assignment of `cxx:116` is never read afterwards
    /// (the loop recomputes `pnext` before its next use); it is kept for
    /// fidelity with the OCCT source.
    #[allow(unused_assignments)]
    fn make_preciser(
        &self,
        u: &mut f64,
        p: &GpPnt,
        is_min: bool,
        ortog_loc: &GpPnt,
        ortog_dir: &GpDir,
    ) {
        if *u > self.my_usup {
            *u = self.my_usup;
        } else if *u < self.my_uinf {
            *u = self.my_uinf;
        } else {
            let mut step = (self.my_usup - self.my_uinf) / 30.0;
            let dir = GpVec::from_xyz(self.my_direction.xyz());
            let pe = project_pnt(ortog_loc, ortog_dir, &dir, &self.my_c.d0(*u));
            let pprev = project_pnt(ortog_loc, ortog_dir, &dir, &self.my_c.d0(*u - step));
            let mut pnext = project_pnt(ortog_loc, ortog_dir, &dir, &self.my_c.d0(*u + step));
            let mut d2e = p.square_distance(&pe);
            let mut d2next = p.square_distance(&pnext);
            let d2prev = p.square_distance(&pprev);
            let mut not_found = if is_min {
                d2e > d2prev || d2e > d2next
            } else {
                d2e < d2prev || d2e < d2next
            };

            if not_found && (d2e < d2next && is_min) {
                step = -step;
                d2next = d2prev;
                pnext = pprev;
            }
            while not_found {
                *u += step;
                if *u > self.my_usup {
                    *u = self.my_usup;
                    break;
                }
                if *u < self.my_uinf {
                    *u = self.my_uinf;
                    break;
                }
                d2e = d2next;
                pnext = project_pnt(ortog_loc, ortog_dir, &dir, &self.my_c.d0(*u + step));
                d2next = p.square_distance(&pnext);
                not_found = if is_min { d2e > d2next } else { d2e < d2next };
            }
        }
    }

    /// `Extrema_ExtPExtS::IsDone()` (`cxx:442-445`).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// `Extrema_ExtPExtS::NbExt()` (`cxx:449-463`).
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

    /// `Extrema_ExtPExtS::SquareDistance(N)` (`cxx:467-484`).
    pub fn square_distance(&self, n: usize) -> f64 {
        if self.my_is_analytically_computable {
            self.my_sq_dist[n - 1]
        } else {
            self.my_ext_ps.square_distance(n)
        }
    }

    /// `Extrema_ExtPExtS::Point(N)` (`cxx:488-505`) as `(U, V, P)`.
    pub fn point(&self, n: usize) -> (f64, f64, GpPnt) {
        if self.my_is_analytically_computable {
            let (u, v) = self.my_point[n - 1].parameter();
            (u, v, self.my_point[n - 1].value())
        } else {
            self.my_ext_ps.point(n)
        }
    }
}

/// `GetPosition(C)` (`Extrema_ExtPExtS.cxx:509-534`): the supporting plane
/// of the elementary basis, as `(location, normal)`.
fn get_position(c: &dyn Curve) -> CurvePlane {
    if let Some(l) = c.gp_line() {
        // gp_Lin L = C->Line(); gp_Pln Pln(L.Location(), L.Direction());
        // gp_Ax2 Pos(Pln.Location(), Pln.Position().Direction(), ...);
        return CurvePlane {
            loc: l.location(),
            dir: l.direction(),
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

/// `IsCaseAnalyticallyComputable` (`Extrema_ExtPExtS.cxx:567-585`): the five
/// elementary types, and the curve's supporting plane not parallel to the
/// extrusion direction.
fn is_case_analytically_computable(
    c: &dyn Curve,
    curve_pos: &CurvePlane,
    surface_direction: &GpDir,
) -> bool {
    if !is_elementary_curve(c) {
        return false;
    }
    (curve_pos.dir.dot(surface_direction)).abs() > RESOLUTION
}
