//! `BRepBlend_RstRstLineBuilder::Transition` and `MakeExtremity`.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_RstRstLineBuilder.cxx:1452-1586`
//! (OCCT 8.0.0).
//!
//! `Transition` (`cxx:1452-1516`) is the surface/restriction sibling of the one
//! in [`super::super::surf_rst_line_builder`]: the arc tangent is read as
//! `Arc->D1`, the surface normal is `d1u ^ d1v`, and the line tangent is either
//! `previousP.TangentOnC1/2()` or, at a tangency point, the chord from the point
//! before it (`cxx:1462-1477`).
//!
//! `MakeExtremity` (`cxx:1522-1586`) stores the point on `surf1` when `OnFirst`
//! and on `surf2` otherwise, both with `sol` as the curve parameter, and then
//! attaches the given arc or, at a vertex, every arc of the domain carrying that
//! vertex.
//!
//! Translation notes:
//! - The `OnFirst == false` branch sets the extremity tangent from
//!   `previousP.TangentOnC1()` (`cxx:1545`), not `TangentOnC2`; the surrounding
//!   `OnFirst` branch right above uses `TangentOnC1` as well (`cxx:1536`). OCCT
//!   8.0.0 does that, so the port keeps it.
//! - `Arc->D1` failing and a null surface normal make OCCT raise in `gp_Dir`;
//!   the port leaves the transitions undecided and attaches the arc anyway, as
//!   the sibling builders do.

use occt_core::gp::GpDir;
use occt_core::gp::GpVec;

use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_blend::{BlendTool, Extremity};
use crate::brep_top_adaptor_hvertex::HVertex;
use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::intsurf_transition::{make_transition, IntSurfTransition};

use super::RstRstLineBuilder;

impl RstRstLineBuilder {
    /// `Transition(OnFirst, Arc, Param, TLine, TArc)` (`cxx:1452-1516`).
    pub(super) fn transition(
        &self,
        on_first: bool,
        arc: &Curve2dAdaptor,
        param: f64,
        tline: &mut IntSurfTransition,
        tarc: &mut IntSurfTransition,
    ) {
        let mut chord = false;
        let mut prevprev = None;

        if self.previous_p.is_tangency_point() {
            if self.line_nb_points() < 2 {
                return;
            }
            chord = true;
            prevprev = if self.sens < 0.0 {
                self.line_point(2)
            } else {
                self.line_point(self.line_nb_points() - 1)
            };
        }

        let Some((p2d, dp2d)) = arc.d1(param) else {
            return;
        };

        let (d1u, d1v) = if on_first {
            let (_p, d1u, d1v) = self.surf1.d1(p2d.x(), p2d.y());
            (d1u, d1v)
        } else {
            let (_p, d1u, d1v) = self.surf2.d1(p2d.x(), p2d.y());
            (d1u, d1v)
        };

        // `tgrst.SetLinearForm(dp2d.X(), d1u, dp2d.Y(), d1v)` (`cxx:1510`).
        let tgrst = GpVec::zero().set_linear_form_2(dp2d.x(), &d1u, dp2d.y(), &d1v);
        // `normale = d1u.Crossed(d1v)` (`cxx:1511`); `gp_Dir` normalises it.
        let Ok(normale) = GpDir::from_vec(&d1u.crossed(&d1v)) else {
            return;
        };

        let tgline = if !chord {
            self.previous_tangent_of(on_first)
        } else {
            let Some(prevprev) = prevprev else {
                return;
            };
            if on_first {
                GpVec::from_pnts(prevprev.point_on_c1(), self.previous_p.point_on_c1())
            } else {
                GpVec::from_pnts(prevprev.point_on_c2(), self.previous_p.point_on_c2())
            }
        };

        make_transition(&tgline, &tgrst, &normale, tline, tarc);
    }

    /// `MakeExtremity(Extrem, OnFirst, Arc, Param, IsVtx, Vtx)`
    /// (`cxx:1522-1586`).
    pub(super) fn make_extremity(
        &self,
        extrem: &mut Extremity,
        on_first: bool,
        arc: &Curve2dAdaptor,
        param: f64,
        is_vtx: bool,
        vtx: Option<&HVertex>,
    ) {
        let mut tline = IntSurfTransition::new();
        let mut tarc = IntSurfTransition::new();

        let iter = if on_first {
            self.domain1.clone()
        } else {
            self.domain2.clone()
        };

        if on_first {
            extrem.set_on_curve(
                *self.previous_p.point_on_c1(),
                self.sol.value(1),
                self.previous_p.parameter(),
                self.tolpoint3d,
            );
            if !self.previous_p.is_tangency_point() {
                if let Ok(tg) = self.previous_p.tangent_on_c1() {
                    extrem.set_tangent(*tg);
                }
            }
        } else {
            extrem.set_on_curve(
                *self.previous_p.point_on_c2(),
                self.sol.value(2),
                self.previous_p.parameter(),
                self.tolpoint3d,
            );
            if !self.previous_p.is_tangency_point() {
                // `TangentOnC1()` in the `else` branch (`cxx:1545`).
                if let Ok(tg) = self.previous_p.tangent_on_c1() {
                    extrem.set_tangent(*tg);
                }
            }
        }

        iter.borrow_mut().init();
        if !is_vtx {
            self.transition(on_first, arc, param, &mut tline, &mut tarc);
            if let Some(curve) = super::super::walking::arc_curve(arc) {
                extrem.add_arc(curve, param, tline, tarc);
            }
            return;
        }

        // `Extrem.SetVertex(Vtx)` is called unguarded in OCCT; the port can only
        // be here with a vertex handle.
        let Some(vtx) = vtx else {
            return;
        };
        extrem.set_vertex(vtx.clone());
        loop {
            if !iter.borrow().more() {
                break;
            }
            let Some(other) = iter.borrow().value().cloned() else {
                break;
            };
            if !Self::same_arc(&other, arc) {
                iter.borrow_mut().initialize_curve(&other);
                iter.borrow_mut().init_vertex_iterator();
                while iter.borrow().more_vertex() {
                    let iv = iter.borrow().vertex();
                    if let Some(iv) = iv {
                        if TopolTool::identical(vtx, &iv) {
                            let prm = BlendTool::parameter(vtx, &other);
                            self.transition(on_first, &other, prm, &mut tline, &mut tarc);
                            if let Some(curve) = super::super::walking::arc_curve(&other) {
                                extrem.add_arc(curve, prm, tline, tarc);
                            }
                        }
                    }
                    iter.borrow_mut().next_vertex();
                }
            } else {
                self.transition(on_first, &other, param, &mut tline, &mut tarc);
                if let Some(curve) = super::super::walking::arc_curve(&other) {
                    extrem.add_arc(curve, param, tline, tarc);
                }
            }
            iter.borrow_mut().next();
        }
    }
}
