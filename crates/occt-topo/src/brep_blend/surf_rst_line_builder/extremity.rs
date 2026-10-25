//! `BRepBlend_SurfRstLineBuilder::Transition` and `MakeExtremity`.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_SurfRstLineBuilder.cxx:1293-1423`
//! (OCCT 8.0.0).
//!
//! `Transition` (`cxx:1293-1355`) rebuilds the line/arc transitions of one
//! extremity from the mixed product of the surface tangent, the arc tangent and
//! the surface normal. For a tangency point it falls back to the chord between
//! `previousP` and the point before it (`cxx:1310-1318`), and it gives up when
//! the line does not have the two points that chord needs (`cxx:1304-1307`).
//!
//! `MakeExtremity` (`cxx:1359-1423`) stores the point -- on `surf1` when
//! `OnFirst`, on `rst` otherwise -- and its parameter, then either attaches the
//! arc it was given or, at a vertex, every arc of the domain that carries that
//! vertex.

use occt_core::gp::GpDir;
use occt_core::gp::GpVec;

use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_blend::{BlendTool, Extremity};
use crate::brep_top_adaptor_hvertex::HVertex;
use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::intsurf_transition::{make_transition, IntSurfTransition};

use super::SurfRstLineBuilder;

impl SurfRstLineBuilder {
    /// `Transition(OnFirst, Arc, Param, TLine, TArc)` (`cxx:1293-1355`).
    ///
    /// A null arc tangent or a null surface normal makes OCCT raise in
    /// `gp_Dir`; as in [`super::super::cs_walking`] the port leaves the
    /// transitions undecided instead.
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

        // `tgrst.SetLinearForm(dp2d.X(), d1u, dp2d.Y(), d1v)` (`cxx:1351`).
        let tgrst = GpVec::zero().set_linear_form_2(dp2d.x(), &d1u, dp2d.y(), &d1v);
        // `normale = d1u.Crossed(d1v)` (`cxx:1352`); `gp_Dir` normalises it.
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
                GpVec::from_pnts(prevprev.point_on_s1(), self.previous_p.point_on_s1())
            } else {
                GpVec::from_pnts(prevprev.point_on_s2(), self.previous_p.point_on_s2())
            }
        };

        make_transition(&tgline, &tgrst, &normale, tline, tarc);
    }

    /// `MakeExtremity(Extrem, OnFirst, Arc, Param, IsVtx, Vtx)`
    /// (`cxx:1359-1423`).
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
            extrem.set_on_surface(
                *self.previous_p.point_on_s(),
                self.sol.value(1),
                self.sol.value(2),
                self.previous_p.parameter(),
                self.tolpoint3d,
            );
            if !self.previous_p.is_tangency_point() {
                if let Ok(tg) = self.previous_p.tangent_on_s() {
                    extrem.set_tangent(*tg);
                }
            }
        } else {
            extrem.set_on_curve(
                *self.previous_p.point_on_c(),
                self.sol.value(3),
                self.previous_p.parameter(),
                self.tolpoint3d,
            );
            if !self.previous_p.is_tangency_point() {
                if let Ok(tg) = self.previous_p.tangent_on_c() {
                    extrem.set_tangent(*tg);
                }
            }
        }

        iter.borrow_mut().init();
        if !is_vtx {
            self.transition(on_first, arc, param, &mut tline, &mut tarc);
            if arc.curve().is_some() {
                extrem.add_arc(arc, param, tline, tarc);
            }
            return;
        }

        // `Extrem.SetVertex(Vtx)` is called unguarded in OCCT; the port can
        // only be here with a vertex handle.
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
                            if other.curve().is_some() {
                                extrem.add_arc(&other, prm, tline, tarc);
                            }
                        }
                    }
                    iter.borrow_mut().next_vertex();
                }
            } else {
                self.transition(on_first, &other, param, &mut tline, &mut tarc);
                if other.curve().is_some() {
                    extrem.add_arc(&other, param, tline, tarc);
                }
            }
            iter.borrow_mut().next();
        }
    }
}
