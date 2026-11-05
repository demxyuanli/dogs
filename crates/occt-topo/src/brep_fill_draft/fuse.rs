//! `BRepFill_Draft::Fuse` (`BRepFill_Draft.cxx:538-824`): the Boolean between
//! the skin and the stopping shape.
//!
//! The reference runs one `BOPAlgo_PaveFiller`, hands it to
//! `BRepAlgoAPI_Section` (`cxx:599`), to `BOPAlgo_Builder::PerformWithFiller`
//! (`cxx:694`) and to the gluer (`cxx:750-757`). Only public builder entry
//! points exist in the port, and each of them runs the intersection phase
//! itself; since `PaveFiller::perform` depends on the arguments and options
//! alone, every run fills the same data structure. `SetGlue(GlueShift)`
//! (`cxx:752`) is forwarded through
//! [`BopBuilder::set_glue`](crate::bop_builder2::BopBuilder::set_glue), so the
//! gluing itself is the ported `GlueEnum::Shift` path.
//!
//! `aHistory->Merge(aLO, aBuilder)` (`cxx:747`) restricts the cut history to
//! the first argument; the port merges the whole cut history
//! ([`BopHistory::merge_algorithm`]), which differs only for shapes the code
//! never asks about.

use occt_core::precision::Precision;
use occt_geom::geom_lprop::sl_props::SlProps;

use crate::abs::ShapeType;
use crate::bop_builder2::{BoolOp2, BopBuilder};
use crate::bop_hist::BopHistory;
use crate::brep_extrema_dist_shape_shape::DistShapeShape;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::fclass2d::FaceState;
use crate::pave_filler::GlueEnum;
use crate::shape::{Face, Shell, Solid, TopoShape};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::explorer_shapes;

use super::BrepFillDraft;

impl BrepFillDraft {
    /// `Fuse(StopShape, KeepOutSide)` (`cxx:538-824`).
    ///
    /// `Ok(false)` is the reference's `return false`; the `Err` arm is the
    /// `TopoDS::Solid` / `TopoDS::Shell` / `TopoDS::Face` cast failures, which
    /// the reference raises as `Standard_TypeMismatch`.
    pub(super) fn fuse(
        &mut self,
        stop_shape: &TopoShape,
        keep_out_side: bool,
    ) -> Result<bool, String> {
        let b = TopoBuilder::new();

        // `cxx:540-556`.
        let mut is_solid = false;
        let sol1: Solid = match self.my_shape.clone() {
            None => return Ok(false),
            Some(shape) if shape.shape_type() == ShapeType::Solid => {
                is_solid = true;
                Solid::wrap(shape).ok_or("Standard_TypeMismatch: TopoDS::Solid")?
            }
            Some(shape) => {
                // `B.MakeSolid(Sol1); B.Add(Sol1, myShape);`
                let shell = Shell::wrap(shape).ok_or("Standard_TypeMismatch: TopoDS::Shell")?;
                b.make_solid(&[shell])
            }
        };

        // `cxx:556-584`.
        let sol2: Solid = match stop_shape.shape_type() {
            ShapeType::Compound => {
                // `TopoDS_Iterator It(StopShape); return Fuse(It.Value(), ...)`
                return match crate::bop_occt_util::iter_children(stop_shape).first().cloned() {
                    Some(first) => self.fuse(&first, keep_out_side),
                    None => Ok(false),
                };
            }
            ShapeType::Solid => {
                Solid::wrap(stop_shape.clone()).ok_or("Standard_TypeMismatch: TopoDS::Solid")?
            }
            ShapeType::Shell => {
                let shell =
                    Shell::wrap(stop_shape.clone()).ok_or("Standard_TypeMismatch: TopoDS::Shell")?;
                b.make_solid(&[shell])
            }
            ShapeType::Face => {
                // `cxx:572-579`.
                let face = Face::wrap(stop_shape.clone()).ok_or("Standard_TypeMismatch: TopoDS::Face")?;
                let shell = b.make_shell(&[face]);
                shell.0.set_closed(BRepTool::is_closed(&shell.0));
                b.make_solid(&[shell])
            }
            // `cxx:581-583`: impossible to do.
            _ => return Ok(false),
        };

        // `cxx:586-597`: `aPF.SetArguments(anArgs); aPF.Perform();`
        // `cxx:599-601`: `BRepAlgoAPI_Section aSec(Sol1, Sol2, aPF)`.
        let mut section = BopBuilder::new();
        section.set_arguments(&[sol1.0.clone(), sol2.0.clone()]);
        let a_section = match section.perform_section(&[sol1.0.clone()], &[sol2.0.clone()]) {
            Ok(s) => s,
            Err(_) => return Ok(false),
        };

        // `cxx:602-607`.
        let section_edges = explorer_shapes(&a_section, ShapeType::Edge);
        if section_edges.is_empty() {
            return Ok(false);
        }

        // `cxx:542-543`: both states start `TopAbs_OUT`.
        let mut state1 = FaceState::Out;
        let mut state2 = FaceState::Out;

        // `cxx:609-674`: choose the state by the geometry when the stop shape is
        // not a solid.
        if stop_shape.shape_type() != ShapeType::Solid {
            let mut a_se_min: Option<crate::shape::Edge> = None;
            let mut d_min = Precision::INFINITE;
            let mut dist = DistShapeShape::default();
            dist.load_s1(&self.my_wire.0);

            for a_se in &section_edges {
                dist.load_s2(a_se);
                dist.perform();
                if dist.is_done() {
                    let d = dist.value();
                    if d < d_min {
                        d_min = d;
                        a_se_min = crate::shape::Edge::wrap(a_se.clone());
                        if d_min < Precision::CONFUSION {
                            break;
                        }
                    }
                }
            }

            if let Some(a_se_min) = a_se_min {
                // `cxx:641-673`.
                // `BRep_Tool::CurveOnSurface(aSEMin, C2d, S, L, f, l, 2)`:
                // the second `BRep_CurveOnSurface` of the edge.
                let reps = GeometryRegistry::global().edge_pcurve_reps(&a_se_min.0);
                if let Some((surface, c2d, f, l)) = reps.get(1).cloned() {
                    let p2d = c2d.d0((f + l) / 2.0);
                    let mut sp = SlProps::new(surface.as_ref(), p2d.x(), p2d.y(), 1, 1.0e-12)?;
                    if !sp.is_normal_defined(surface.as_ref()) {
                        let p2d = c2d.d0((3.0 * f + l) / 4.0);
                        sp.set_parameters(surface.as_ref(), p2d.x(), p2d.y());
                        if !sp.is_normal_defined(surface.as_ref()) {
                            let p2d = c2d.d0((f + 3.0 * l) / 4.0);
                            sp.set_parameters(surface.as_ref(), p2d.x(), p2d.y());
                        }
                    }

                    if sp.is_normal_defined(surface.as_ref()) {
                        let normal = sp.normal(surface.as_ref())?;
                        // `cxx:665-671`: "Subtract State1".
                        state1 = if self.my_dir.angle(&normal) < std::f64::consts::PI / 2.0 {
                            FaceState::In
                        } else {
                            FaceState::Out
                        };
                    }
                }
            }
        }

        // `cxx:678-688`: "Invert State2".
        if !keep_out_side {
            state2 = if state2 == FaceState::In {
                FaceState::Out
            } else {
                FaceState::In
            };
        }

        // `cxx:691-698`: `aBuilder.AddArgument(Sol1/Sol2); PerformWithFiller(aPF)`.
        let mut bop = BopBuilder::new();
        bop.set_arguments(&[sol1.0.clone(), sol2.0.clone()]);
        if bop.perform().is_err() {
            return Ok(false);
        }
        // The arguments the filler registered, i.e. the shapes the data
        // structure knows (`BOPAlgo_Builder::BuildBOP`'s "unknown shape" check
        // reads the same list).
        let args = bop.arguments().to_vec();
        if args.len() < 2 {
            return Ok(false);
        }

        let mut history = BopHistory::new();
        let mut result: Option<TopoShape> = None;
        let mut single_op_needed = true;

        // `cxx:703-769`: "To get rid of the unnecessary parts of first solid
        // make the cutting first".
        if state1 == FaceState::Out {
            if bop
                .build_bop_operation(&args[0..1], &args[1..2], BoolOp2::Cut)
                .is_ok()
            {
                // `cxx:713-741`: the cut solid closest to `myWire`.
                let solids = explorer_shapes(bop.result(), ShapeType::Solid);
                let mut a_cut_min: Option<Solid> = None;
                if let Some(first) = solids.first() {
                    a_cut_min = Solid::wrap(first.clone());
                    if solids.len() > 1 {
                        let mut d_min = Precision::INFINITE;
                        let mut dist = DistShapeShape::default();
                        dist.load_s1(&self.my_wire.0);
                        for cut in &solids {
                            dist.load_s2(cut);
                            dist.perform();
                            if dist.is_done() {
                                let d = dist.value();
                                if d < d_min {
                                    d_min = d;
                                    a_cut_min = Solid::wrap(cut.clone());
                                }
                            }
                        }
                    }
                }

                if let Some(a_cut_min) = a_cut_min {
                    // `cxx:745-747`: "Save history for first argument only".
                    history.merge_algorithm(&[sol1.0.clone()], bop.history());

                    // `cxx:749-757`.
                    let mut gluer = BopBuilder::new();
                    gluer.set_arguments(&[a_cut_min.0.clone(), sol2.0.clone()]);
                    gluer.set_glue(GlueEnum::Shift);
                    let _ = gluer.perform();
                    let gargs = gluer.arguments().to_vec();
                    let glued = if gargs.len() < 2 {
                        false
                    } else {
                        gluer
                            .build_bop(&gargs[0..1], state1, &gargs[1..2], state2)
                            .is_ok()
                    };

                    // `cxx:759-768`.
                    if glued {
                        history.merge(gluer.history());
                        let r = gluer.result().clone();
                        single_op_needed = explorer_shapes(&r, ShapeType::Solid).is_empty();
                        result = Some(r);
                    }
                }
            }
        }

        // `cxx:771-787`.
        if single_op_needed {
            history = BopHistory::new();
            if bop.build_bop(&args[0..1], state1, &args[1..2], state2).is_err() {
                return Ok(false);
            }
            history.merge(bop.history());
            result = Some(bop.result().clone());
        }

        let Some(result) = result else {
            return Ok(false);
        };

        // `cxx:789-800`.
        if is_solid {
            self.my_shape = Some(result);
        } else {
            // `cxx:793-800`: keep the first shell of the result.
            if let Some(shell) = explorer_shapes(&result, ShapeType::Shell).first() {
                self.my_shape = Some(shell.clone());
            }
        }

        // `cxx:803-820`: "Update the History".
        let nb_law = self.nb_law();
        for ii in 1..=nb_law {
            let old = self.my_faces.value(1, ii).cloned();
            if let Some(old) = old {
                if let Some(first) = history.modified(&old).and_then(|l| l.first().cloned()) {
                    self.my_faces.set_value(1, ii, Some(first));
                }
            }
        }
        for ii in 1..=nb_law + 1 {
            let old = self.my_sections.value(1, ii).cloned();
            if let Some(old) = old {
                if let Some(first) = history.modified(&old).and_then(|l| l.first().cloned()) {
                    self.my_sections.set_value(1, ii, Some(first));
                }
            }
        }

        Ok(true)
    }
}
