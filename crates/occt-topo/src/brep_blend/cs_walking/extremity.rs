//! `BRepBlend_CSWalking::Transition` and `MakeExtremity`.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_CSWalking.cxx:741-822`
//! (OCCT 8.0.0).
//!
//! `Transition` (`cxx:741-759`) rebuilds the line/arc transitions of one
//! extremity from the mixed product of the surface tangent, the arc tangent and
//! the surface normal. Unlike `BRepBlend_Walking::Transition` it calls
//! `CSLib` nowhere: the normal is the plain cross product of the first
//! derivatives (`cxx:757`), and `IntSurf::MakeTransition` (`cxx:759`) does the
//! classification.
//!
//! `MakeExtremity` (`cxx:761-822`) stores the surface point and the guide
//! parameter, then either attaches the `Index`-th arc of the domain or, at a
//! vertex, every arc of the domain that carries that vertex. That vertex form is
//! only reachable through the recadrage call OCCT comments out (`cxx:894`), so
//! the port keeps it for the day that call comes back.

use super::CsWalking;
use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_blend::{BlendTool, Extremity};
use crate::brep_top_adaptor_hvertex::HVertex;
use crate::intsurf_transition::{make_transition, IntSurfTransition};
use occt_core::gp::GpDir;
use occt_core::gp::GpVec;

impl CsWalking {
    /// `Transition(A, Param, TLine, TArc)` (`cxx:741-759`).
    ///
    /// A null arc tangent or a null surface normal makes OCCT raise in
    /// `gp_Dir`; the port leaves the transitions undecided instead, the same
    /// way [`super::super::walking::Walking::transition`] does.
    pub(super) fn transition(
        &self,
        arc: &Curve2dAdaptor,
        param: f64,
        tline: &mut IntSurfTransition,
        tarc: &mut IntSurfTransition,
    ) {
        let Some((p2d, dp2d)) = arc.d1(param) else {
            return;
        };
        let (_pbid, d1u, d1v) = self.surf.d1(p2d.x(), p2d.y());

        // `tgrst.SetLinearForm(dp2d.X(), d1u, dp2d.Y(), d1v)` (`cxx:756`).
        let tgrst = GpVec::zero().set_linear_form_2(dp2d.x(), &d1u, dp2d.y(), &d1v);
        // `normale = d1u.Crossed(d1v)` (`cxx:757`); `gp_Dir` normalises it.
        let normale = match GpDir::from_vec(&d1u.crossed(&d1v)) {
            Ok(n) => n,
            Err(_) => return,
        };

        let tgline = self.previous_tangent();
        make_transition(&tgline, &tgrst, &normale, tline, tarc);
    }

    /// `MakeExtremity(Extrem, Index, Param, IsVtx, Vtx)` (`cxx:761-822`).
    pub(super) fn make_extremity(
        &self,
        extrem: &mut Extremity,
        index: usize,
        param: f64,
        is_vtx: bool,
        vtx: Option<&HVertex>,
    ) {
        let mut tline = IntSurfTransition::new();
        let mut tarc = IntSurfTransition::new();

        let (u, v) = self.previous_parameters_on_s();
        extrem.set_on_surface(
            *self.previous_p.point_on_s(),
            u,
            v,
            self.previous_p.parameter(),
            self.tolpoint3d,
        );

        let mut iter = self.domain.borrow_mut();
        iter.init();
        let mut nbarc = 1usize;

        if !is_vtx {
            while nbarc < index {
                if !iter.more() {
                    return;
                }
                nbarc += 1;
                iter.next();
            }
            let Some(arc) = iter.value().cloned() else {
                return;
            };
            self.transition(&arc, param, &mut tline, &mut tarc);
            if arc.curve().is_some() {
                extrem.add_arc(&arc, param, tline, tarc);
            }
        } else {
            let Some(vertex) = vtx else {
                return;
            };
            extrem.set_vertex(vertex.clone());
            while iter.more() {
                let Some(arc) = iter.value().cloned() else {
                    return;
                };
                if nbarc != index {
                    iter.initialize_curve(&arc);
                    iter.init_vertex_iterator();
                    while iter.more_vertex() {
                        let identical = match iter.vertex() {
                            Some(v) => crate::brep_top_adaptor_topol_tool::TopolTool::identical(
                                vertex, &v,
                            ),
                            None => false,
                        };
                        if identical {
                            let prm = BlendTool::parameter(vertex, &arc);
                            self.transition(&arc, prm, &mut tline, &mut tarc);
                            if arc.curve().is_some() {
                                extrem.add_arc(&arc, prm, tline, tarc);
                            }
                        }
                        iter.next_vertex();
                    }
                } else {
                    self.transition(&arc, param, &mut tline, &mut tarc);
                    if arc.curve().is_some() {
                        extrem.add_arc(&arc, param, tline, tarc);
                    }
                }
                nbarc += 1;
                iter.next();
            }
        }
    }
}
