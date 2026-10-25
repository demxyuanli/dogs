//! `BRepBlend_Walking::Transition`, `MakeExtremity`, `MakeSingularExtremity`
//! and `CorrectExtremityOnOneRst`.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_Walking.cxx:1612-1817`
//! (transitions and extremities) and `:2608-2769` (the corner-plane correction
//! of an extremity that landed on a restriction), OCCT 8.0.0.
//!
//! `Transition` (`cxx:1612-1726`) rebuilds the line/arc transitions of one
//! extremity from the mixed product of the chord and the arc tangent with the
//! surface normal, using `CSLib::Normal` on the first derivatives and falling
//! back to the order-2 expansion of `CSLib` when they are degenerate.
//!
//! `CorrectExtremityOnOneRst` (`cxx:2608-2769`) is the corner fix: the ends of
//! the arc are projected on the guide, the closest spine vertex gives the corner
//! plane, and the partner point is projected on that plane and back on its
//! surface. `TolAng = 0.001` is the "bug OCC25701" value of `cxx:2619`, not a
//! port threshold.

use std::rc::Rc;

use occt_core::cslib::{normal_d1_mag, normal_max_order, CSLibNormalStatus};
use occt_core::gp::{GpAx1, GpPnt, GpVec};
use occt_core::precision::INFINITE;
use occt_geom::extrema_pc::point_curve_extrema_all;
use occt_geom::extrema_surf::point_surface_extrema;
use occt_geom::intana::plane_through_three;

use crate::blend::BlendPoint;
use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_blend::{BlendTool, Extremity};
use crate::brep_top_adaptor_hvertex::HVertex;
use crate::intsurf_transition::{make_transition, IntSurfTransition};

use super::Walking;

/// `CSLib::Normal` sine tolerance of `Transition` (`cxx:1673`, `:1706`).
const SIN_TOL: f64 = 1.0e-9;
/// `TolAng` of `CorrectExtremityOnOneRst` (`cxx:2619`, bug OCC25701).
const TOL_ANG: f64 = 0.001;

impl Walking {
    /// `Transition(OnFirst, A, Param, TLine, TArc)` (`cxx:1612-1726`).
    ///
    /// The transitions are left UNDECIDED when OCCT would raise: a tangency
    /// point with fewer than two points in the line returns early (`cxx:1622-1626`),
    /// and a tangent or a normal that cannot be evaluated is not recoverable.
    pub(super) fn transition(
        &self,
        on_first: bool,
        arc: &Curve2dAdaptor,
        param: f64,
        tline: &mut IntSurfTransition,
        tarc: &mut IntSurfTransition,
    ) {
        let mut compute_with_chord = false;
        let mut prevprev: Option<BlendPoint> = None;

        if self.previous_p.is_tangency_point() {
            let Some(line) = self.line.as_ref() else {
                return;
            };
            if line.nb_points() < 2 {
                return;
            }
            compute_with_chord = true;
            prevprev = if self.sens < 0.0 {
                line.point(2).copied()
            } else {
                line.point(line.nb_points() - 1).copied()
            };
        }

        let Some((p2d, dp2d)) = arc.d1(param) else {
            return;
        };

        let surf = if on_first { &self.surf1 } else { &self.surf2 };
        let (_pbid, d1u, d1v) = surf.d1(p2d.x(), p2d.y());

        let tgline = if !compute_with_chord {
            match self.previous_tangent(on_first) {
                Some(v) => v,
                None => return,
            }
        } else {
            let Some(prevprev) = prevprev else {
                return;
            };
            if on_first {
                GpVec::from_pnts(prevprev.point_on_s1(), self.previous_p.point_on_s1())
            } else {
                GpVec::from_pnts(prevprev.point_on_s2(), self.previous_p.point_on_s2())
            }
        };

        // `tgrst.SetLinearForm(dp2d.X(), d1u, dp2d.Y(), d1v)` (`cxx:1672`).
        let tgrst = GpVec::zero().set_linear_form_2(dp2d.x(), &d1u, dp2d.y(), &d1v);

        // `CSLib::Normal(d1u, d1v, 1.e-9, stat, thenormal)` (`cxx:1674`) is the
        // magnitude-guarded overload (`CSLib.cxx:156-179`).
        let (stat, thenormal) = normal_d1_mag(&d1u, &d1v, SIN_TOL);
        let normale = match thenormal {
            Some(n) if stat == CSLibNormalStatus::Defined => n,
            _ => {
                // Fallback to the order-2 expansion (`cxx:1679-1724`).
                let (ux, vx) = (p2d.x(), p2d.y());
                let (_p, d1u, d1v, d2u, d2v, d2uv) = surf.d2(ux, vx);
                let mut der = vec![vec![GpVec::zero(); 3]; 3];
                der[1][0] = d1u;
                der[0][1] = d1v;
                der[2][0] = d2u;
                der[0][2] = d2v;
                der[1][1] = d2uv;
                // `Adaptor3d_HSurfaceTool::DN` (`cxx:1702-1704`) — the surface's
                // `EvalDN`. A surface that does not implement it raises there, so
                // the port leaves `Transition` altogether.
                let dn = |nu: i32, nv: i32| surf.try_eval_dn(ux, vx, nu, nv).ok();
                let (Some(d21), Some(d12), Some(d22)) = (dn(2, 1), dn(1, 2), dn(2, 2)) else {
                    return;
                };
                der[2][1] = d21;
                der[1][2] = d12;
                der[2][2] = d22;
                let (u_min, u_max) = surf.u_range();
                let (v_min, v_max) = surf.v_range();
                let (_st, nrm, _iu, _iv) = normal_max_order(
                    2, &der, SIN_TOL, ux, vx, u_min, u_max, v_min, v_max,
                );
                match nrm {
                    Some(n) => n,
                    // `gp_Dir Norm(Normale)` on a null vector raises in OCCT.
                    None => return,
                }
            }
        };

        make_transition(&tgline, &tgrst, &normale, tline, tarc);
    }

    /// `MakeExtremity(Extrem, OnFirst, Index, Param, IsVtx, Vtx)`
    /// (`cxx:1728-1773`).
    ///
    /// `IsVtx` is only ever true together with a vertex found by `Recadre`
    /// (`cxx:1585-1600`), so the null handle OCCT would store cannot occur.
    pub(super) fn make_extremity(
        &self,
        extrem: &mut Extremity,
        on_first: bool,
        index: usize,
        param: f64,
        is_vtx: bool,
        vtx: Option<&HVertex>,
    ) {
        let mut tline = IntSurfTransition::new();
        let mut tarc = IntSurfTransition::new();

        if on_first {
            extrem.set_on_surface(
                *self.previous_p.point_on_s1(),
                self.sol.value(1),
                self.sol.value(2),
                self.previous_p.parameter(),
                self.tolpoint3d,
            );
            if let Some(tg) = self.previous_tangent(true) {
                extrem.set_tangent(tg);
            }
        } else {
            extrem.set_on_surface(
                *self.previous_p.point_on_s2(),
                self.sol.value(3),
                self.sol.value(4),
                self.previous_p.parameter(),
                self.tolpoint3d,
            );
            if let Some(tg) = self.previous_tangent(false) {
                extrem.set_tangent(tg);
            }
        }

        let tool = if on_first {
            Rc::clone(&self.recdomain1)
        } else {
            Rc::clone(&self.recdomain2)
        };
        let Some(arc) = Self::domain_arc(&tool, index) else {
            return;
        };

        self.transition(on_first, &arc, param, &mut tline, &mut tarc);
        if let Some(curve) = super::arc_curve(&arc) {
            extrem.add_arc(curve, param, tline, tarc);
        }
        if is_vtx {
            if let Some(v) = vtx {
                extrem.set_vertex(v.clone());
            }
        }
    }

    /// `MakeSingularExtremity(Extrem, OnFirst, Vtx)` (`cxx:1775-1817`): the
    /// extremity of the partner surface at a vertex where both restrictions
    /// meet, built from every arc carrying that vertex.
    pub(super) fn make_singular_extremity(
        &self,
        extrem: &mut Extremity,
        on_first: bool,
        vtx: &HVertex,
    ) {
        let tool = if on_first {
            Rc::clone(&self.recdomain1)
        } else {
            Rc::clone(&self.recdomain2)
        };

        if let Some(tg) = self.previous_tangent(on_first) {
            extrem.set_tangent(tg);
        }

        let mut iter = tool.borrow_mut();
        iter.init();
        extrem.set_vertex(vtx.clone());
        while iter.more() {
            let Some(arc) = iter.value().cloned() else {
                break;
            };
            iter.initialize_curve(&arc);
            iter.init_vertex_iterator();
            while iter.more_vertex() {
                let Some(current) = iter.vertex() else {
                    break;
                };
                if current.is_same(vtx) {
                    let prm = BlendTool::parameter(vtx, &arc);
                    let mut tline = IntSurfTransition::new();
                    let mut tarc = IntSurfTransition::new();
                    self.transition(on_first, &arc, prm, &mut tline, &mut tarc);
                    if let Some(curve) = super::arc_curve(&arc) {
                        extrem.add_arc(curve, prm, tline, tarc);
                    }
                }
                iter.next_vertex();
            }
            iter.next();
        }
    }

    /// `CorrectExtremityOnOneRst(IndexOfRst, theU, theV, theParam, thePntOnRst,
    /// NewU, NewV, NewPoint, NewParam)` (`cxx:2608-2769`): `None` where OCCT
    /// answers `false`.
    pub(super) fn correct_extremity_on_one_rst(
        &self,
        index_of_rst: i32,
        the_u: f64,
        the_v: f64,
        the_param: f64,
        the_pnt_on_rst: &GpPnt,
    ) -> Option<(f64, f64, GpPnt, f64)> {
        if self.hguide.nb_vertices() == 0 {
            return None;
        }

        let (domain_of_rst, surf_of_rst, another_surf) = if index_of_rst == 1 {
            (
                Rc::clone(&self.recdomain1),
                std::sync::Arc::clone(&self.surf1),
                std::sync::Arc::clone(&self.surf2),
            )
        } else {
            (
                Rc::clone(&self.recdomain2),
                std::sync::Arc::clone(&self.surf2),
                std::sync::Arc::clone(&self.surf1),
            )
        };

        // `DomainOfRst->Value()`, read without an `Init()` (`cxx:2632`): the
        // iterator is left on the arc `Recadre` selected.
        let thearc = domain_of_rst.borrow().value().cloned()?;
        let thecurve = super::arc_curve(&thearc)?;
        let ends = [thecurve.first_parameter(), thecurve.last_parameter()];

        let mut global_min_sq_dist = INFINITE;
        let mut param_on_guide = 0.0f64;
        let mut point_on_guide = GpPnt::new(0.0, 0.0, 0.0);
        for &end in ends.iter() {
            let p2d_on_end = thecurve.d0(end);
            let pnt_on_end = surf_of_rst.d0(p2d_on_end.x(), p2d_on_end.y());
            // `Extrema_ExtPC projoncurv(PntOnEnd, theElSpine)` (`cxx:2645`);
            // the minimum `SquareDistance` is the first of the sorted extrema.
            let Some(best) = point_curve_extrema_all(self.hguide.as_ref(), &pnt_on_end)
                .into_iter()
                .next()
            else {
                continue;
            };
            let min_sq_dist = best.distance * best.distance;
            if min_sq_dist < global_min_sq_dist {
                global_min_sq_dist = min_sq_dist;
                param_on_guide = best.u1;
                point_on_guide = best.p2;
            }
        }

        let mut new_param = param_on_guide;
        if self.hguide.is_periodic() {
            let period = self.hguide.period().unwrap_or(0.0);
            let sign = if new_param < the_param { 1.0 } else { -1.0 };
            while (new_param - the_param).abs() > period / 2.0 {
                new_param += sign * period;
            }
        }

        // The closest spine vertex gives the corner axis (`cxx:2680-2691`).
        let mut global_min_sq_dist = INFINITE;
        let mut the_ax1: Option<GpAx1> = None;
        for ind in 1..=self.hguide.nb_vertices() {
            let Some(an_ax1) = self.hguide.vertex_with_tangent(ind) else {
                continue;
            };
            let a_sq_dist = point_on_guide.square_distance(an_ax1.location());
            if a_sq_dist < global_min_sq_dist {
                global_min_sq_dist = a_sq_dist;
                the_ax1 = Some(*an_ax1);
            }
        }
        let the_ax1 = the_ax1?;
        let pnt0 = *the_ax1.location();
        let dir0 = *the_ax1.direction();

        // Check the new point: is it a real solution? (`cxx:2695-2712`).
        let old_p_on_guide = self.hguide.eval_d0(the_param)?;
        let pnt_on_surf2 = another_surf.d0(the_u, the_v);
        let old_plane = plane_through_three(the_pnt_on_rst, &old_p_on_guide, &pnt_on_surf2)?;
        let old_axis = old_plane.axis();
        let old_dir = old_axis.direction();
        let mut angle = old_dir.angle(&dir0);
        if angle > std::f64::consts::PI / 2.0 {
            angle = std::f64::consts::PI - angle;
        }
        if angle > TOL_ANG {
            return None;
        }

        // Project (theU, theV) on the plane (Pnt0, Dir0) (`cxx:2714-2717`).
        let a_vec = GpVec::from_pnts(&pnt0, &pnt_on_surf2);
        let dir0_vec = GpVec::from_xyz(dir0.xyz());
        let a_translation =
            GpVec::from_xyz(&a_vec.xyz().multiply_xyz(dir0_vec.xyz()).multiply_xyz(dir0_vec.xyz()));
        let pnt_on_plane = pnt_on_surf2.translated_vec(&a_translation.reversed());

        // Check the new point again (`cxx:2719-2733`).
        let new_plane = plane_through_three(the_pnt_on_rst, &pnt0, &pnt_on_plane)?;
        let new_axis = new_plane.axis();
        let dir_of_new_plane = new_axis.direction();
        angle = dir0.angle(&dir_of_new_plane);
        if angle > std::f64::consts::PI / 2.0 {
            angle = std::f64::consts::PI - angle;
        }
        if angle > TOL_ANG {
            return None;
        }

        // Project PntOnPlane on the partner surface (`cxx:2736-2766`).
        let new_p_on_surf2 = point_surface_extrema(another_surf.as_ref(), &pnt_on_plane)?;
        let new_point = new_p_on_surf2.p2;
        let mut new_u = new_p_on_surf2.u2;
        let mut new_v = new_p_on_surf2.v2?;
        let uperiod = if another_surf.is_u_periodic() {
            another_surf.u_period()
        } else {
            0.0
        };
        let vperiod = if another_surf.is_v_periodic() {
            another_surf.v_period()
        } else {
            0.0
        };
        recadre_if_periodic(&mut new_u, &mut new_v, the_u, the_v, uperiod, vperiod);
        Some((new_u, new_v, new_point, new_param))
    }
}

/// `RecadreIfPeriodic(NewU, NewV, OldU, OldV, UPeriod, VPeriod)`
/// (`cxx:1822-1846`): bring the new parameters to the period of the old ones.
pub(super) fn recadre_if_periodic(
    new_u: &mut f64,
    new_v: &mut f64,
    old_u: f64,
    old_v: f64,
    u_period: f64,
    v_period: f64,
) {
    if u_period > 0.0 {
        let sign = if *new_u < old_u { 1.0 } else { -1.0 };
        while (*new_u - old_u).abs() > u_period / 2.0 {
            *new_u += sign * u_period;
        }
    }
    if v_period > 0.0 {
        let sign = if *new_v < old_v { 1.0 } else { -1.0 };
        while (*new_v - old_v).abs() > v_period / 2.0 {
            *new_v += sign * v_period;
        }
    }
}
