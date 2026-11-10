//! `BRepBlend_Walking::InternalPerform` and its `evalpinit` helper.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_Walking.cxx:1848-2605`
//! (OCCT 8.0.0): `evalpinit` at `:1848-1911`, `InternalPerform` at `:1913-2605`.
//!
//! `InternalPerform` marches from `parameter` towards `Bound` one section at a
//! time. Each iteration checks the 3d angular deflection of the guide
//! (`CosRef3D = 0.88`, `cxx:1820`), solves `Blend_Function` with the previous
//! solution extrapolated along the last 2d tangents (`evalpinit`), classifies
//! the result on both domains, recadres the sections that left a domain onto
//! the nearest restriction (`Recadre`), corrects extremities that land on a
//! corner (`CorrectExtremityOnOneRst`), decides with `TestArret` and updates the
//! line. `stepw` is halved on a too large step and grown by 1.5 on a too small
//! one, between the guide tolerance and `pasmax`.
//!
//! Translation notes:
//! - `Ext1.SetValue(P, U, V, Tol)` is the four-argument curve overload of
//!   `BRepBlend_Extremity` (`hxx:81-84`), i.e. `set_on_curve` here.
//! - `Precision::IsInfinite` is [`occt_core::precision::Precision::is_infinite`].

use occt_core::gp::{GpDir, GpVec};
use occt_core::precision::{Precision, REAL_SMALL};
use occt_geom::Curve;
use occt_math::MathVector;

use crate::blend::{BlendFuncInv, BlendFunction, BlendPoint, BlendStatus};
use crate::brep_blend::Extremity;
use crate::brep_top_adaptor_hvertex::HVertex;
use crate::fclass2d::FaceState;

use super::{arc_curve, classify_domain, solve_root_set, Walking};

/// `static const double CosRef3D = 0.88;` (`cxx:1820`): the guide deflection
/// limit of the walking loop.
const COS_REF3D: f64 = 0.88;
/// `const double TolProd = 1.e-5` of the same-direction test (`cxx:2191`).
const TOL_PROD: f64 = 1.0e-5;

impl Walking {
    /// `InternalPerform(Func, FuncInv, Bound)` (`cxx:1913-2605`).
    pub(super) fn internal_perform(
        &mut self,
        func: &mut dyn BlendFunction,
        func_inv: &mut dyn BlendFuncInv,
        bound: f64,
    ) {
        let sens = self.sens;
        let tolgui = self.tolgui;
        let mut stepw = self.pasmax;

        let nbp = self.line.as_ref().map(|l| l.nb_points()).unwrap_or(0);
        if nbp >= 2 {
            // Take the last step back when it is not too small (`cxx:1921-1931`).
            let line = self.line.as_ref().expect("line");
            stepw = if sens < 0.0 {
                line.point(2).map(|p| p.parameter()).unwrap_or(0.0)
                    - line.point(1).map(|p| p.parameter()).unwrap_or(0.0)
            } else {
                line.point(nbp).map(|p| p.parameter()).unwrap_or(0.0)
                    - line.point(nbp - 1).map(|p| p.parameter()).unwrap_or(0.0)
            };
            stepw = stepw.max(100.0 * tolgui);
        }

        let mut parprec = self.param;
        let (_pt_on_guide, mut prev_tg_on_guide) = self.hguide.d1(parprec);

        if sens * (parprec - bound) >= -tolgui {
            return;
        }

        let mut tolerance = MathVector::with_init(1, 4, 0.0);
        let mut infbound = MathVector::with_init(1, 4, 0.0);
        let mut supbound = MathVector::with_init(1, 4, 0.0);
        let mut parinit = self.sol.clone();
        func.get_tolerance(&mut tolerance, self.tolpoint3d);
        func.get_bounds(&mut infbound, &mut supbound);
        // `math_FunctionSetRoot rsnld(Func, tolerance, 30)` (`cxx:1960`); the
        // solver copies the vector, so every iteration reuses this one.
        let solver_toler = tolerance.clone();

        let mut arrive = false;
        let mut state = BlendStatus::OnRst12;
        let mut situ1 = FaceState::In;
        let mut situ2 = FaceState::In;
        let mut w1;
        let mut w2;
        let mut index1 = 0usize;
        let mut index2 = 0usize;
        let mut recad1;
        let mut recad2;
        let mut is_vtx1 = false;
        let mut is_vtx2 = false;
        let mut vtx1: Option<HVertex> = None;
        let mut vtx2: Option<HVertex> = None;
        let mut solrst1 = MathVector::with_init(1, 4, 0.0);
        let mut solrst2 = MathVector::with_init(1, 4, 0.0);
        let mut ext1 = Extremity::new();
        let mut ext2 = Extremity::new();

        self.param = parprec + sens * stepw;
        if sens * (self.param - bound) > 0.0 {
            stepw = sens * (bound - parprec) * 0.5;
            self.param = parprec + sens * stepw;
        }
        evalpinit(
            &mut parinit,
            &self.previous_p,
            parprec,
            self.param,
            &infbound,
            &supbound,
            self.clasons1,
            self.clasons2,
        );

        while !arrive {
            let (_pt, tg_on_guide) = self.hguide.d1(self.param);
            // Check deflection on guide (`cxx:1981-2007`).
            let cosi = prev_tg_on_guide.dot(&tg_on_guide);
            let cosi2 = if cosi < REAL_SMALL {
                0.0
            } else {
                cosi * cosi / prev_tg_on_guide.square_magnitude() / tg_on_guide.square_magnitude()
            };
            if cosi2 < COS_REF3D {
                // Angle 3d too great.
                state = BlendStatus::StepTooLarge;
                stepw /= 2.0;
                self.param = parprec + sens * stepw; // Cannot pass Bound.
                if stepw.abs() < tolgui {
                    self.set_both_extremities(&mut ext1, &mut ext2);
                    arrive = true;
                }
                continue;
            }
            prev_tg_on_guide = tg_on_guide;

            let mut bonpoint = true;
            func.set(self.param);
            match solve_root_set(func, &parinit, &solver_toler, &infbound, &supbound, 30) {
                None => {
                    state = BlendStatus::StepTooLarge;
                    bonpoint = false;
                }
                Some(root) => {
                    self.sol = root;
                    situ1 = if self.clasons1 {
                        classify_domain(
                            &self.domain1,
                            self.sol.value(1),
                            self.sol.value(2),
                            tolerance.value(1).min(tolerance.value(2)),
                            false,
                        )
                    } else {
                        FaceState::In
                    };
                    situ2 = if self.clasons2 {
                        classify_domain(
                            &self.domain2,
                            self.sol.value(3),
                            self.sol.value(4),
                            tolerance.value(3).min(tolerance.value(4)),
                            false,
                        )
                    } else {
                        FaceState::In
                    };
                }
            }

            let nb_points = self.line.as_ref().map(|l| l.nb_points()).unwrap_or(0);
            if bonpoint && nb_points == 1 && (situ1 != FaceState::In || situ2 != FaceState::In) {
                state = BlendStatus::StepTooLarge;
                bonpoint = false;
            }

            if bonpoint {
                w1 = bound;
                w2 = bound;
                recad1 = false;
                recad2 = false;
                let mut echecrecad = false;
                let mut control = false;

                if situ1 == FaceState::Out || situ1 == FaceState::On {
                    // Inverse problem on the first surface (`cxx:2058-2088`).
                    let snapshot = self.sol.clone();
                    match self.recadre(func_inv, true, &snapshot, 0.0) {
                        Some(r) => {
                            recad1 = true;
                            solrst1 = r.solrst;
                            index1 = r.index;
                            is_vtx1 = r.is_vtx;
                            vtx1 = r.vtx;
                            let wtemp = solrst1.value(2);
                            if (self.param - wtemp) / sens >= -10.0 * tolgui {
                                w1 = solrst1.value(2);
                                control = true;
                            } else {
                                // The wrong solution was taken.
                                echecrecad = true;
                                recad1 = false;
                                // cxx:1951-1952 -- State/bonpoint are written here but
                                // both are overwritten at cxx:2216-2229 before being read.
                                stepw /= 2.0;
                            }
                        }
                        None => echecrecad = true,
                    }
                }
                if situ2 == FaceState::Out || situ2 == FaceState::On {
                    // Inverse problem on the second surface (`cxx:2090-2116`).
                    let snapshot = self.sol.clone();
                    match self.recadre(func_inv, false, &snapshot, 0.0) {
                        Some(r) => {
                            recad2 = true;
                            solrst2 = r.solrst;
                            index2 = r.index;
                            is_vtx2 = r.is_vtx;
                            vtx2 = r.vtx;
                            let wtemp = solrst2.value(2);
                            if (self.param - wtemp) / sens >= -10.0 * tolgui {
                                w2 = solrst2.value(2);
                                control = true;
                            } else {
                                echecrecad = true;
                                recad2 = false;
                                // cxx:1981-1982 -- same dead store as cxx:1951-1952.
                                stepw /= 2.0;
                            }
                        }
                        None => echecrecad = true,
                    }
                }

                // Which recadrage has to be checked (`cxx:2119-2146`).
                if recad1 && recad2 {
                    if (w1 - w2).abs() <= 10.0 * tolgui {
                        // No need to check the recadrages: the check can fail
                        // (cf model blend10).
                        control = false;
                    } else if sens * (w1 - w2) < 0.0 {
                        recad2 = false;
                    } else {
                        recad1 = false;
                    }
                }

                // Effective check of the recadrages (`cxx:2143-2165`).
                if control {
                    if recad1 && self.clasons2 {
                        let situ = classify_domain(
                            &self.recdomain2,
                            solrst1.value(3),
                            solrst1.value(4),
                            tolerance.value(3).min(tolerance.value(4)),
                            true,
                        );
                        if situ == FaceState::Out {
                            recad1 = false;
                            echecrecad = true;
                        }
                    } else if recad2 && self.clasons1 {
                        let situ = classify_domain(
                            &self.recdomain1,
                            solrst2.value(3),
                            solrst2.value(4),
                            tolerance.value(1).min(tolerance.value(1)),
                            true,
                        );
                        if situ == FaceState::Out {
                            recad2 = false;
                            echecrecad = true;
                        }
                    }
                }

                if recad1 || recad2 {
                    echecrecad = false;
                }

                if !echecrecad {
                    if recad1 && recad2 {
                        state = self.on_both_restrictions(&solrst1, &solrst2, w1, w2);
                    } else if recad1 {
                        state = BlendStatus::OnRst1;
                        self.param = w1;
                        if let Some(p2d) = self.arc_point_on(true, solrst1.value(1)) {
                            self.sol.set_value(1, p2d.x());
                            self.sol.set_value(2, p2d.y());
                        }
                        self.sol.set_value(3, solrst1.value(3));
                        self.sol.set_value(4, solrst1.value(4));
                        let the_pnt_on_rst = self.surf1.d0(self.sol.value(1), self.sol.value(2));
                        if let Some((u, v, _pnt, param)) = self.correct_extremity_on_one_rst(
                            1,
                            self.sol.value(3),
                            self.sol.value(4),
                            self.param,
                            &the_pnt_on_rst,
                        ) {
                            self.param = param;
                            self.sol.set_value(3, u);
                            self.sol.set_value(4, v);
                        }
                    } else if recad2 {
                        state = BlendStatus::OnRst2;
                        self.param = w2;
                        self.sol.set_value(1, solrst2.value(3));
                        self.sol.set_value(2, solrst2.value(4));
                        if let Some(p2d) = self.arc_point_on(false, solrst2.value(1)) {
                            self.sol.set_value(3, p2d.x());
                            self.sol.set_value(4, p2d.y());
                        }
                        let the_pnt_on_rst = self.surf2.d0(self.sol.value(3), self.sol.value(4));
                        if let Some((u, v, _pnt, param)) = self.correct_extremity_on_one_rst(
                            2,
                            self.sol.value(1),
                            self.sol.value(2),
                            self.param,
                            &the_pnt_on_rst,
                        ) {
                            self.param = param;
                            self.sol.set_value(1, u);
                            self.sol.set_value(2, v);
                        }
                    } else {
                        state = BlendStatus::Ok;
                    }

                    // A non-orthodox step is better than no recadrage (PMN).
                    if recad1 || recad2 {
                        func.set(self.param);
                        state = self.test_arret(
                            func,
                            state,
                            stepw.abs() > 3.0 * tolgui,
                            false,
                            true,
                        );
                    } else {
                        state = self.test_arret(func, state, true, true, false);
                    }
                } else if stepw > 2.0 * tolgui {
                    // The maximal step may be badly set: divide.
                    state = BlendStatus::StepTooLarge;
                } else {
                    // Recadrage failure: leave with PointsConfondus.
                    state = BlendStatus::SamePoints;
                }
            }

            match state {
                BlendStatus::Ok => {
                    self.record_point(sens);
                    parprec = self.param;
                    if self.param == bound {
                        arrive = true;
                        // End on Bound.
                        self.set_both_extremities(&mut ext1, &mut ext2);
                    } else {
                        self.param += sens * stepw;
                        if sens * (self.param - bound) > -tolgui {
                            self.param = bound;
                        }
                    }
                    evalpinit(
                        &mut parinit,
                        &self.previous_p,
                        parprec,
                        self.param,
                        &infbound,
                        &supbound,
                        self.clasons1,
                        self.clasons2,
                    );
                }
                BlendStatus::StepTooLarge => {
                    stepw /= 2.0;
                    if stepw.abs() < tolgui {
                        self.set_both_extremities(&mut ext1, &mut ext2);
                        arrive = true;
                    } else {
                        self.param = parprec + sens * stepw;
                        evalpinit(
                            &mut parinit,
                            &self.previous_p,
                            parprec,
                            self.param,
                            &infbound,
                            &supbound,
                            self.clasons1,
                            self.clasons2,
                        );
                    }
                }
                BlendStatus::StepTooSmall => {
                    self.record_point(sens);
                    parprec = self.param;
                    stepw = (1.5 * stepw).min(self.pasmax);
                    if self.param == bound {
                        arrive = true;
                        // End on Bound.
                        self.set_both_extremities(&mut ext1, &mut ext2);
                    } else {
                        self.param += sens * stepw;
                        if sens * (self.param - bound) > -tolgui {
                            self.param = bound;
                        }
                    }
                    evalpinit(
                        &mut parinit,
                        &self.previous_p,
                        parprec,
                        self.param,
                        &infbound,
                        &supbound,
                        self.clasons1,
                        self.clasons2,
                    );
                }
                BlendStatus::OnRst1 => {
                    self.record_point(sens);
                    self.make_extremity(
                        &mut ext1,
                        true,
                        index1,
                        solrst1.value(1),
                        is_vtx1,
                        vtx1.as_ref(),
                    );
                    // The singular case where one of the recadrages failed
                    // (`cxx:2499-2512`).
                    if self
                        .previous_p
                        .point_on_s1()
                        .is_equal_with_tolerance(self.previous_p.point_on_s2(), 2.0 * self.tolpoint3d)
                    {
                        ext2.set_on_curve(
                            *self.previous_p.point_on_s1(),
                            self.sol.value(3),
                            self.sol.value(4),
                            self.tolpoint3d,
                        );
                        if is_vtx1 {
                            if let Some(v) = vtx1.as_ref() {
                                self.make_singular_extremity(&mut ext2, false, v);
                            }
                        }
                    } else {
                        ext2.set_on_surface(
                            *self.previous_p.point_on_s2(),
                            self.sol.value(3),
                            self.sol.value(4),
                            self.previous_p.parameter(),
                            self.tolpoint3d,
                        );
                    }
                    arrive = true;
                }
                BlendStatus::OnRst2 => {
                    self.record_point(sens);
                    if self
                        .previous_p
                        .point_on_s1()
                        .is_equal_with_tolerance(self.previous_p.point_on_s2(), 2.0 * self.tolpoint3d)
                    {
                        ext1.set_on_curve(
                            *self.previous_p.point_on_s2(),
                            self.sol.value(1),
                            self.sol.value(2),
                            self.tolpoint3d,
                        );
                        if is_vtx2 {
                            if let Some(v) = vtx2.as_ref() {
                                self.make_singular_extremity(&mut ext1, true, v);
                            }
                        }
                    } else {
                        ext1.set_on_surface(
                            *self.previous_p.point_on_s1(),
                            self.sol.value(1),
                            self.sol.value(2),
                            self.previous_p.parameter(),
                            self.tolpoint3d,
                        );
                    }
                    self.make_extremity(
                        &mut ext2,
                        false,
                        index2,
                        solrst2.value(1),
                        is_vtx2,
                        vtx2.as_ref(),
                    );
                    arrive = true;
                }
                BlendStatus::OnRst12 => {
                    self.record_point(sens);
                    if is_vtx1 != is_vtx2
                        && self
                            .previous_p
                            .point_on_s1()
                            .is_equal_with_tolerance(
                                self.previous_p.point_on_s2(),
                                2.0 * self.tolpoint3d,
                            )
                    {
                        // Only one recadrage recognised as a vertex
                        // (`cxx:2557-2569`).
                        if is_vtx1 {
                            is_vtx2 = true;
                            vtx2 = vtx1.clone();
                        } else {
                            is_vtx1 = true;
                            vtx1 = vtx2.clone();
                        }
                    }
                    self.make_extremity(
                        &mut ext1,
                        true,
                        index1,
                        solrst1.value(1),
                        is_vtx1,
                        vtx1.as_ref(),
                    );
                    self.make_extremity(
                        &mut ext2,
                        false,
                        index2,
                        solrst2.value(1),
                        is_vtx2,
                        vtx2.as_ref(),
                    );
                    arrive = true;
                }
                BlendStatus::SamePoints => {
                    // Points confused during the walking.
                    self.set_both_extremities(&mut ext1, &mut ext2);
                    arrive = true;
                }
                BlendStatus::Backward => {}
            }

            if arrive {
                if let Some(line) = self.line.as_mut() {
                    if sens > 0.0 {
                        line.set_end_points(ext1.clone(), ext2.clone());
                    } else {
                        line.set_start_points(ext1.clone(), ext2.clone());
                    }
                }
            }
        }
    }

    /// The `recad1 && recad2` arm of `InternalPerform` (`cxx:2176-2264`): the
    /// solution lies on both restrictions, possibly on a corner of the guide.
    fn on_both_restrictions(
        &mut self,
        solrst1: &MathVector,
        solrst2: &MathVector,
        w1: f64,
        w2: f64,
    ) -> BlendStatus {
        self.param = (w1 + w2) / 2.0;
        // Through the arcs, to avoid the periodic surfaces' problems.
        if let Some(p2d) = self.arc_point_on(true, solrst1.value(1)) {
            self.sol.set_value(1, p2d.x());
            self.sol.set_value(2, p2d.y());
        }
        let pnt1 = self.surf1.d0(self.sol.value(1), self.sol.value(2));
        if let Some(p2d) = self.arc_point_on(false, solrst2.value(1)) {
            self.sol.set_value(3, p2d.x());
            self.sol.set_value(4, p2d.y());
        }
        let pnt2 = self.surf2.d0(self.sol.value(3), self.sol.value(4));

        // Check the original first and last parameters of the guide curve for
        // equality with <param>: the tangent of the guide and the normal of the
        // plane built on the guide point and the two restriction points have to
        // be parallel (`cxx:2192-2229`).
        let saved = [
            self.hguide.get_saved_first_parameter(),
            self.hguide.get_saved_last_parameter(),
        ];
        let mut same_dirs = [false, false];
        for (ind, &saved_param) in saved.iter().enumerate() {
            if Precision::is_infinite(saved_param) {
                continue;
            }
            let (pnt0, dir0) = self.hguide.d1(saved_param);
            let length = dir0.magnitude();
            if length <= REAL_SMALL {
                continue;
            }
            let Ok(dir0) = GpDir::from_vec(&dir0.divided(length)) else {
                continue;
            };
            let Some(plane) = occt_geom::intana::plane_through_three(&pnt0, &pnt1, &pnt2) else {
                continue;
            };
            let plane_axis = plane.axis();
            let dir_plane = plane_axis.direction();
            let the_prod =
                GpVec::from_xyz(dir0.xyz()).crossed(&GpVec::from_xyz(dir_plane.xyz()));
            let prod_mod = the_prod.magnitude();
            if prod_mod <= TOL_PROD {
                same_dirs[ind] = true;
            }
        }

        // Choose the closest parameter (`cxx:2230-2263`).
        let mut the_param = Precision::INFINITE;
        if same_dirs[0] && same_dirs[1] {
            the_param = if (self.param - saved[0]).abs() < (self.param - saved[1]).abs() {
                saved[0]
            } else {
                saved[1]
            };
        } else if same_dirs[0] {
            the_param = saved[0];
        } else if same_dirs[1] {
            the_param = saved[1];
        }

        if let Some((_u, _v, _pnt, new_param)) = self.correct_extremity_on_one_rst(
            1,
            self.sol.value(3),
            self.sol.value(4),
            self.param,
            &pnt1,
        ) {
            if (self.param - new_param).abs() < (self.param - the_param).abs() {
                the_param = new_param;
            }
        }

        if !Precision::is_infinite(the_param) {
            self.param = the_param;
        }

        BlendStatus::OnRst12
    }

    /// The four-argument `Extremity::SetValue(P, U, V, Tol)` pairs of
    /// `cxx:1999-2000`, `:2400-2405`, `:2426-2431`, `:2583-2587`: both
    /// extremities are put back on the previous section.
    fn set_both_extremities(&self, ext1: &mut Extremity, ext2: &mut Extremity) {
        ext1.set_on_surface(
            *self.previous_p.point_on_s1(),
            self.sol.value(1),
            self.sol.value(2),
            self.previous_p.parameter(),
            self.tolpoint3d,
        );
        ext2.set_on_surface(
            *self.previous_p.point_on_s2(),
            self.sol.value(3),
            self.sol.value(4),
            self.previous_p.parameter(),
            self.tolpoint3d,
        );
        if let Some(tg) = self.previous_tangent(true) {
            ext1.set_tangent(tg);
        }
        if let Some(tg) = self.previous_tangent(false) {
            ext2.set_tangent(tg);
        }
    }

    /// `line->Append` / `line->Prepend` of `previousP` (`cxx:2384-2392`).
    fn record_point(&mut self, sens: f64) {
        if let Some(line) = self.line.as_mut() {
            if sens > 0.0 {
                line.append(self.previous_p);
            } else {
                line.prepend(self.previous_p);
            }
        }
    }

    /// `p2d = BRepBlend_HCurve2dTool::Value(recdomain->Value(), u)`: the arc the
    /// iterator was left on, without an `Init()` (`cxx:2182`, `:2187`, `:2276`,
    /// `:2313`).
    fn arc_point_on(&self, on_first: bool, arc_param: f64) -> Option<occt_core::gp::GpPnt2d> {
        let tool = if on_first {
            &self.recdomain1
        } else {
            &self.recdomain2
        };
        let arc = Self::current_arc(tool)?;
        let curve = arc_curve(&arc)?;
        Some(curve.d0(arc_param))
    }
}

/// `evalpinit(parinit, previousP, parprec, param, infbound, supbound, classonS1,
/// classonS2)` (`cxx:1848-1911`): extrapolate the previous solution along the
/// 2d tangents, falling back to the previous solution when the guess leaves the
/// domain it is still classified on.
#[allow(clippy::too_many_arguments)]
fn evalpinit(
    parinit: &mut MathVector,
    previous_p: &BlendPoint,
    parprec: f64,
    param: f64,
    infbound: &MathVector,
    supbound: &MathVector,
    classon_s1: bool,
    classon_s2: bool,
) {
    if previous_p.is_tangency_point() {
        if let Ok((u, v)) = previous_p.parameters_on_s1() {
            parinit.set_value(1, u);
            parinit.set_value(2, v);
        }
        if let Ok((u, v)) = previous_p.parameters_on_s2() {
            parinit.set_value(3, u);
            parinit.set_value(4, v);
        }
        return;
    }

    let mut inside = true;
    let (mut u1, mut v1) = previous_p.parameters_on_s1().unwrap_or((0.0, 0.0));
    let (mut u2, mut v2) = previous_p.parameters_on_s2().unwrap_or((0.0, 0.0));
    // `Tangent2dOnS1/2()` raise on a tangency point, which the branch above
    // already took.
    let (du1, dv1) = match previous_p.tangent2d_on_s1() {
        Ok(t) => (t.coord.x, t.coord.y),
        Err(_) => (0.0, 0.0),
    };
    let (du2, dv2) = match previous_p.tangent2d_on_s2() {
        Ok(t) => (t.coord.x, t.coord.y),
        Err(_) => (0.0, 0.0),
    };

    let step = param - parprec;
    u1 += step * du1;
    v1 += step * dv1;
    if classon_s1 {
        if u1 < infbound.value(1) || u1 > supbound.value(1) {
            inside = false;
        }
        if v1 < infbound.value(2) || v1 > supbound.value(2) {
            inside = false;
        }
    }
    u2 += step * du2;
    v2 += step * dv2;
    if classon_s2 {
        if u2 < infbound.value(3) || u2 > supbound.value(3) {
            inside = false;
        }
        if v2 < infbound.value(4) || v2 > supbound.value(4) {
            inside = false;
        }
    }

    if inside {
        parinit.set_value(1, u1);
        parinit.set_value(2, v1);
        parinit.set_value(3, u2);
        parinit.set_value(4, v2);
    } else if let Ok((u, v)) = previous_p.parameters_on_s1() {
        // Do not play the smart one (`cxx:1904-1909`).
        parinit.set_value(1, u);
        parinit.set_value(2, v);
    }
    if !inside {
        if let Ok((u, v)) = previous_p.parameters_on_s2() {
            parinit.set_value(3, u);
            parinit.set_value(4, v);
        }
    }
}
