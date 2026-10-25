//! `TopOpeBRep_FFTransitionTool` — the face/face transition factory.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRep/TopOpeBRep_FFTransitionTool.hxx:30-69`,
//! `.cxx:31-435`.
//!
//! Translation notes:
//! - The file-local helpers stay module functions: `TransitionToOrientation`
//!   (`cxx:32-74`), `FUN_ProjectPoint` (`cxx:304-327`) and `FUN_GeomTrans`
//!   (`cxx:335-388`).
//! - `TT.Set(O)` is [`Transition::set_orientation`]
//!   (`TopOpeBRepDS_Transition.cxx:205-229`); `TT.Set(TopAbs_UNKNOWN,
//!   TopAbs_UNKNOWN)` is the four-argument `Set` with its `TopAbs_FACE`
//!   defaults (`hxx:42-45`), so the shapes are written to `Face` too.
//! - `TT.Before(S)` / `TT.After(S)` take OCCT's default `TopAbs_FACE` shape
//!   (`hxx:54-58`).
//! - `TopoDS::Edge(RR)` / `TopoDS::Face(FF)` are `static_cast`-style
//!   downcasts; the port takes the typed `Edge` / `Face` parameters instead.
//! - `Index` outside `{1, 2}` leaves OCCT's local `IntSurf_Transition T`
//!   default-constructed (UNDECIDED), which `TransitionToOrientation` reports
//!   as undefined; the port keeps the same default value.
//! - `P.Transition*()` / `P.ParametersOnS*()` dereference the point handle
//!   without a null test in OCCT; the port's accessors answer `Err` when the
//!   point is absent, and the transition then keeps the OCCT default value
//!   instead of crashing.
//! - `ProcessEdgeONTransition` reads the face surface unconditionally in OCCT
//!   (`FUN_GeomTrans` has no null test on `S1`); an unregistered surface makes
//!   the port leave the transition default rather than dereference nothing.
//! - `ProcessLineTransition(P, LI)` indexes `VPoint(1)` / `VPoint(nbv)`; with no
//!   `VPoint` at all OCCT reads out of bounds, so the port answers the default
//!   transition.

use occt_core::gp::GpPnt;
use occt_geom::curve::Curve;

use crate::abs::{Orientation, ShapeType, State};
use crate::brep_tool::BRepTool;
use crate::intsurf_transition::{IntSurfSituation, IntSurfTransition, IntSurfTypeTrans};
use crate::topopebrep::ds::Transition;
use crate::topopebrep::interf::VPointInter;
use crate::topopebrep::line_inter::LineInter;
use crate::{Edge, Face};

/// `TransitionToOrientation` (`cxx:32-74`): the transition of an
/// `IntSurf_Transition` as a `TopAbs_Orientation`, plus OCCT's "defined" flag.
fn transition_to_orientation(t: &IntSurfTransition) -> (Orientation, bool) {
    let mut defined = true;
    let mut result = Orientation::Forward;
    match t.transition_type() {
        IntSurfTypeTrans::In => result = Orientation::Forward,
        IntSurfTypeTrans::Out => result = Orientation::Reversed,
        IntSurfTypeTrans::Touch => {
            // `T.Situation()` raises `Standard_DomainError` unless the
            // transition is a TOUCH one, which it just tested.
            let situa = t.situation().unwrap_or(IntSurfSituation::Unknown);
            match situa {
                IntSurfSituation::Inside => result = Orientation::Internal,
                IntSurfSituation::Outside => result = Orientation::External,
                IntSurfSituation::Unknown => defined = false,
            }
        }
        IntSurfTypeTrans::Undecided => defined = false,
    }
    (result, defined)
}

/// `FUN_ProjectPoint` (`cxx:304-327`): parameter of `p1` on `c2` restricted to
/// `[fc2, lc2]` through `GeomAPI_ProjectPointOnCurve`.
fn fun_project_point(
    p1: &GpPnt,
    c2: Option<&dyn Curve>,
    fc2: f64,
    lc2: f64,
) -> Option<f64> {
    let c2 = c2?;
    occt_geom::extrema_pc::point_curve_extrema_range(c2, p1, fc2, lc2).map(|e| e.u1)
}

/// `FUN_GeomTrans` (`cxx:335-388`): transition on `c1` at `t1` when crossing
/// `c2`, where `S1(u1, v1) == c1(t1)`.
#[allow(clippy::too_many_arguments)]
fn fun_geom_trans(
    s1: &dyn occt_geom::surface::Surface,
    u1: f64,
    v1: f64,
    c1: Option<&dyn Curve>,
    t1: f64,
    c2: Option<&dyn Curve>,
    fc2: f64,
    lc2: f64,
    trans: &mut Transition,
) -> bool {
    let (Some(c1), Some(c2)) = (c1, c2) else {
        return false;
    };

    // P1 = C1(T1), D1_C1 = D1(C1(T1))
    let (p1, d1_c1) = c1.d1(t1);

    let Some(t2) = fun_project_point(&p1, Some(c2), fc2, lc2) else {
        return false;
    };
    // D1_C2 = D1(C2(P1))
    let (_, d1_c2) = c2.d1(t2);

    // N1 = D1(S1(U1,V1))
    let (_, mut d1u, mut d1v) = s1.d1(u1, v1);
    d1u.normalize();
    d1v.normalize();
    let mut n1 = d1u.crossed(&d1v);
    n1.normalize();

    let n1d1_c1 = n1.crossed(&d1_c1);
    let dot = n1d1_c1.dot(&d1_c2);
    if dot > 0.0 {
        trans.before(State::Out, ShapeType::Face);
        trans.after(State::In, ShapeType::Face);
    } else {
        trans.before(State::In, ShapeType::Face);
        trans.after(State::Out, ShapeType::Face);
    }

    true
}

/// `class TopOpeBRep_FFTransitionTool` (`hxx:30-69`).
pub struct FFTransitionTool;

impl FFTransitionTool {
    /// `ProcessLineTransition(P, Index, EdgeOrientation)` (`cxx:86-126`):
    /// transition of the intersection point `P` on the intersected shape of
    /// index `index` (1 or 2) for a line crossing an edge.
    pub fn process_line_transition_with_index(
        p: &VPointInter,
        index: i32,
        edge_orientation: Orientation,
    ) -> Transition {
        let mut tt = Transition::default();

        if edge_orientation == Orientation::Internal || edge_orientation == Orientation::External
        {
            tt.set_orientation(edge_orientation);
        } else {
            let mut t = IntSurfTransition::new();
            match index {
                1 => {
                    if let Ok(v) = p.transition_line_arc1() {
                        t = v;
                    }
                }
                2 => {
                    if let Ok(v) = p.transition_line_arc2() {
                        t = v;
                    }
                }
                _ => {}
            }
            let (mut o, o_defined) = transition_to_orientation(&t);
            if o_defined {
                if edge_orientation == Orientation::Reversed {
                    o = o.complement();
                }
                tt.set_orientation(o);
            } else {
                tt.set(State::Unknown, State::Unknown, ShapeType::Face, ShapeType::Face);
            }
        }
        tt
    }

    /// `ProcessLineTransition(P, LI)` (`cxx:135-165`): transition of `P` on
    /// line `li`, `P` lying on neither of the intersecting shapes.
    pub fn process_line_transition_on_line(p: &VPointInter, li: &LineInter) -> Transition {
        let mut tt = Transition::default();

        // P.IsOnDomS1() and P.IsOnDomS2() are both false.
        let nbv = li.nb_v_point();
        let (Some(p1), Some(pn)) = (li.v_point(1), li.v_point(nbv)) else {
            return tt;
        };
        let par1 = p1.parameter_on_line().unwrap_or(f64::NAN);
        let parn = pn.parameter_on_line().unwrap_or(f64::NAN);

        let par = p.parameter_on_line().unwrap_or(f64::NAN);
        let result = if par == par1 {
            Orientation::Forward
        } else if par == parn {
            Orientation::Reversed
        } else {
            Orientation::Internal
        };

        tt.set_orientation(result);
        tt
    }

    /// `ProcessEdgeTransition(P, Index, LineOrientation)` (`cxx:173-213`): the
    /// transition of `P` on the intersected shape of index `index` (1 or 2).
    pub fn process_edge_transition(
        p: &VPointInter,
        index: i32,
        face_transition: Orientation,
    ) -> Transition {
        let mut tt = Transition::default();

        if face_transition == Orientation::Internal || face_transition == Orientation::External {
            tt.set_orientation(face_transition);
        } else {
            let mut t = IntSurfTransition::new();
            if index == 1 {
                if let Ok(v) = p.transition_on_s1() {
                    t = v;
                }
            } else if index == 2 {
                if let Ok(v) = p.transition_on_s2() {
                    t = v;
                }
            }

            let (mut o, defined) = transition_to_orientation(&t);
            if defined {
                if face_transition == Orientation::Reversed {
                    o = o.complement();
                }
                tt.set_orientation(o);
            } else {
                tt.set(State::Unknown, State::Unknown, ShapeType::Face, ShapeType::Face);
            }
        }

        tt
    }

    /// `ProcessFaceTransition(L, Index, FaceOrientation)` (`cxx:220-297`): the
    /// transition read from the line itself. `index == 1` is the first shape.
    pub fn process_face_transition(
        l: &LineInter,
        index: i32,
        face_orientation: Orientation,
    ) -> Transition {
        let mut tt = Transition::default();

        if face_orientation == Orientation::Internal || face_orientation == Orientation::External {
            tt.set_orientation(face_orientation);
        } else {
            let mut o_defined = true;
            let mut o = Orientation::Forward;

            let trans = if index == 1 {
                l.transition_on_s1()
            } else {
                l.transition_on_s2()
            };

            match trans {
                IntSurfTypeTrans::In => o = Orientation::Forward,
                IntSurfTypeTrans::Out => o = Orientation::Reversed,
                IntSurfTypeTrans::Touch => {
                    let situa = if index == 1 {
                        l.situation_s1()
                    } else {
                        l.situation_s2()
                    };
                    match situa {
                        IntSurfSituation::Inside => o = Orientation::Internal,
                        IntSurfSituation::Outside => o = Orientation::External,
                        IntSurfSituation::Unknown => o_defined = false,
                    }
                }
                IntSurfTypeTrans::Undecided => o_defined = false,
            }

            if o_defined {
                if face_orientation == Orientation::Reversed {
                    o = o.complement();
                }
                tt.set_orientation(o);
            } else {
                tt.set(State::Unknown, State::Unknown, ShapeType::Face, ShapeType::Face);
            }
        }

        tt
    }

    /// `ProcessEdgeONTransition(VP, ShapeIndex, R, E, F)` (`cxx:392-435`):
    /// transition on the `IntPatch_Restriction` edge `r` when crossing edge `e`
    /// of face `f` at point `vp`, `vp` lying ON edge `r`.
    pub fn process_edge_on_transition(
        vp: &VPointInter,
        shape_index: i32,
        r: &Edge,
        e: &Edge,
        f: &Face,
    ) -> Transition {
        let ori_e = e.orientation();

        let s = BRepTool::face_surface(f);
        let mut u = 0.0;
        let mut v = 0.0;
        if shape_index == 1 {
            if let Ok(p) = vp.parameters_on_s1() {
                u = p.0;
                v = p.1;
            }
        } else if shape_index == 2 {
            if let Ok(p) = vp.parameters_on_s2() {
                u = p.0;
                v = p.1;
            }
        }

        let (fe, le) = BRepTool::edge_parameters(e);
        let ce = BRepTool::edge_curve(e);
        let te = vp.edge_parameter(shape_index);

        let (fr, lr) = BRepTool::edge_parameters(r);
        let cr = BRepTool::edge_curve(r);

        let mut trans = Transition::default();
        let trans_ok = match s.as_deref() {
            Some(s) => fun_geom_trans(
                s,
                u,
                v,
                ce.as_deref(),
                te,
                cr.as_deref(),
                fr,
                lr,
                &mut trans,
            ),
            None => false,
        };
        let _ = (fe, le);
        if trans_ok {
            // Trans: transition on R when crossing edge E oriented in face F.
            if ori_e == Orientation::Reversed {
                trans = trans.complement();
            }
        }

        trans
    }
}
