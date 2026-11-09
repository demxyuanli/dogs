//! `BRepFeat_MakePrism::PerformUntilEnd` (`cxx:659-700`),
//! `PerformFromEnd` (`cxx:705-844`) and `PerformThruAll` (`cxx:849-897`).
//!
//! The three semi-infinite / through-all variants. A cut feature answers the
//! Boolean of the basis with the primitive; a fuse feature always goes through
//! `GlobalPerform`, which selects the parts of the tool that stay.

use std::sync::Arc;

use occt_core::precision::PCONFUSION;
use occt_geom::Curve;

use crate::abs::Orientation;
use crate::brep_feat::algo_api::cut;
use crate::brep_feat::form::update_descendants_bop;
use crate::brep_feat::perf_selection::PerfSelection;
use crate::brep_feat::status_error::StatusError;
use crate::brep_feat::statics::tool;
use crate::brep_lib_make_face;
use crate::brep_tool::BRepTool;
use crate::locope::{CSIntersector, Prism};
use crate::shape::{Face, TopoShape};
use crate::topo_tools::{faces_of, is_same};

use super::statics;
use super::MakePrism;

impl MakePrism {
    /// `PerformUntilEnd()` (`cxx:659-700`): a semi-infinite prism limited only
    /// by the basis, so the sweep is a plain two-height extrusion.
    pub fn perform_until_end(&mut self) -> Result<(), String> {
        // `cxx:665-670`.
        self.form.my_perf_selection = PerfSelection::SelectionSh;
        self.form.perf_selection_valid();
        self.form.my_glued_f.clear();
        self.form.my_suntil = None;
        self.form.shape_until_valid();
        self.form.my_sfrom = None;
        self.form.shape_from_valid();

        // `cxx:671-678`.
        let pbase = self
            .my_pbase
            .clone()
            .ok_or_else(|| "StdFail_NotDone".to_string())?;
        let height = statics::height_max(
            self.form.my_sbase.as_ref(),
            self.form.my_skface.as_ref(),
            self.form.my_sfrom.as_ref(),
            self.form.my_suntil.as_ref(),
        );
        let v = self.scaled_dir(2.0 * height);
        let the_prism = Prism::with_base(&pbase, &v);
        let vraiprism = the_prism
            .shape()
            .cloned()
            .ok_or_else(|| "StdFail_NotDone".to_string())?;

        // `cxx:679-682`.
        statics::maj_map(
            &pbase,
            &the_prism,
            &mut self.form.my_map,
            &mut self.form.my_f_shape,
            &mut self.form.my_l_shape,
        );
        self.form.my_gshape = Some(vraiprism);
        self.form.generated_shape_valid();
        self.form.glued_faces_valid();

        // `cxx:684-699`.
        if !self.form.my_fuse {
            if let Some(sbase) = self.form.my_sbase.clone() {
                if let Some(gshape) = self.form.my_gshape.clone() {
                    // `cxx:687-695`: the result is kept only when the Boolean
                    // succeeded; a failure leaves the feature not done.
                    let mut c = cut(&sbase, &gshape);
                    if c.build().is_ok() {
                        if let Some(shape) = c.shape().cloned() {
                            self.form.my_shape = Some(shape.clone());
                            update_descendants_bop(&mut self.form, &c, &shape, false)?;
                            self.form.done();
                        }
                    }
                }
            }
        } else {
            self.take_prism_curves(&the_prism);
            self.global_perform()?;
        }
        Ok(())
    }

    /// `PerformFromEnd(Until)` (`cxx:705-844`): a semi-infinite prism running
    /// away from the basis and cut at `Until`.
    pub fn perform_from_end(&mut self, until: &TopoShape) -> Result<(), String> {
        // `cxx:711-716`.
        if until.is_null() {
            return Err("Standard_ConstructionError".to_string());
        }
        if let Some(sk) = self.form.my_skface.clone() {
            if is_same(until, &sk.0) {
                // `cxx:717-721`.
                self.reverse_dir();
                return self.perform_until_end();
            }
        }
        if faces_of(until).is_empty() {
            return Err("Standard_ConstructionError".to_string());
        }

        // `cxx:727-733`.
        self.form.my_perf_selection = PerfSelection::SelectionShU;
        self.form.perf_selection_valid();
        self.form.my_sfrom = None;
        self.form.shape_from_valid();
        self.form.my_suntil = Some(until.clone());
        let trf = self.form.transform_shape_fu(1);
        self.form.shape_until_valid();

        // `cxx:734-741`.
        let pbase = self
            .my_pbase
            .clone()
            .ok_or_else(|| "StdFail_NotDone".to_string())?;
        let c = statics::test_curve(&pbase, self.dir());
        let sens = match &c {
            Some(c) => statics::sens_of_prism(c, until),
            None => 1,
        };
        let height = statics::height_max(
            self.form.my_sbase.as_ref(),
            self.form.my_skface.as_ref(),
            self.form.my_sfrom.as_ref(),
            self.form.my_suntil.as_ref(),
        );
        let vtra = self.scaled_dir(-3.0 * height * sens as f64 / 2.0);
        let vect = self.scaled_dir(3.0 * sens as f64 * height);
        let the_prism = Prism::with_base_and_translation(&pbase, &vect, &vtra);
        let vraiprism = the_prism
            .shape()
            .cloned()
            .ok_or_else(|| "StdFail_NotDone".to_string())?;

        if !trf {
            // `cxx:743-754`.
            statics::maj_map(
                &pbase,
                &the_prism,
                &mut self.form.my_map,
                &mut self.form.my_f_shape,
                &mut self.form.my_l_shape,
            );
            self.form.my_gshape = Some(vraiprism);
            self.form.generated_shape_valid();
            self.form.my_glued_f.clear();
            self.form.glued_faces_valid();
            self.take_prism_curves(&the_prism);
            self.global_perform()?;
        } else {
            // `cxx:755-843`.
            statics::maj_map(
                &pbase,
                &the_prism,
                &mut self.form.my_map,
                &mut self.form.my_f_shape,
                &mut self.form.my_l_shape,
            );
            let scur: Vec<Option<Arc<dyn Curve>>> = match &c {
                Some(c) => vec![Some(c.clone())],
                None => vec![None],
            };
            let sbase = self.form.my_sbase.clone();
            let mut asi1 = CSIntersector::with_shape(until);
            let mut asi2 = match &sbase {
                Some(sb) => CSIntersector::with_shape(sb),
                None => return Ok(()),
            };
            let _ = asi1.perform_curves(&scur);
            let _ = asi2.perform_curves(&scur);

            // `cxx:764-775`. The two orientations start at FORWARD.
            let mut or_u = Orientation::Forward;
            let mut or_f = Orientation::Forward;
            let mut f_until: Option<Face> = None;
            let mut f_from: Option<Face> = None;
            let nb1 = if asi1.is_done() {
                asi1.nb_points(1).unwrap_or(0)
            } else {
                0
            };
            if nb1 >= 1 {
                or_u = asi1
                    .point(1, 1)
                    .map(|p| p.orientation())
                    .unwrap_or(Orientation::Forward);
                if sens == -1 {
                    or_u = or_u.reversed();
                }
                f_until = asi1.point(1, 1).ok().map(|p| p.face().clone());
            }
            let nb2 = if asi2.is_done() {
                asi2.nb_points(1).unwrap_or(0)
            } else {
                0
            };
            if nb2 >= 1 {
                or_f = asi2
                    .point(1, 1)
                    .map(|p| p.orientation())
                    .unwrap_or(Orientation::Forward);
                f_from = asi2.point(1, 1).ok().map(|p| p.face().clone());
                // `cxx:779-789`: the support face is replaced by an unbounded
                // face on its basis surface, then trimmed by `TransformShapeFU`.
                if let Some(ff) = &f_from {
                    if let Some(s) = BRepTool::face_surface(ff) {
                        let s = s.rectangular_trimmed_basis().unwrap_or(s);
                        let fac = brep_lib_make_face::make_face_from_surface(s, PCONFUSION);
                        self.form.my_sfrom = Some(fac.0.clone());
                        let _ = self.form.transform_shape_fu(0);
                        f_from = Face::wrap(fac.0);
                    }
                }
            }

            // `cxx:791-815`.
            let mut tools: Vec<TopoShape> = Vec::new();
            let tool_u = match (until, f_until) {
                (u, Some(f)) => tool(u, &f, or_u),
                _ => None,
            };
            match tool_u {
                Some(s) => tools.push(s.0),
                None => {
                    self.set_error(StatusError::NullToolU);
                    return Ok(());
                }
            }
            let tool_f = match (&self.form.my_sfrom, f_from) {
                (Some(sf), Some(f)) => tool(sf, &f, or_f),
                _ => None,
            };
            match tool_f {
                Some(s) => tools.push(s.0),
                None => {
                    self.set_error(StatusError::NullToolF);
                    return Ok(());
                }
            }

            // `cxx:817-822`.
            let mut tr_p = crate::brep_feat::algo_api::BooleanOperation::new(
                crate::bop_builder2::BoolOp2::Cut,
            );
            tr_p.set_arguments(&[vraiprism.clone()]);
            tr_p.set_tools(&tools);
            let (tr_p, tr_shape) = statics::boolean_result(tr_p)?;
            update_descendants_bop(&mut self.form, &tr_p, &tr_shape, false)?;

            // `cxx:824-842`.
            self.finish_with_boolean_on(&tr_shape)?;
        }
        Ok(())
    }

    /// `PerformThruAll()` (`cxx:849-897`): the prism is swept through the whole
    /// bounding box, three heights in both directions around the profile.
    pub fn perform_thru_all(&mut self) -> Result<(), String> {
        // `cxx:855-870`.
        self.form.my_suntil = None;
        self.form.shape_until_valid();
        self.form.my_sfrom = None;
        self.form.shape_from_valid();
        if !self.form.my_fuse {
            self.form.my_perf_selection = PerfSelection::NoSelection;
        } else {
            self.form.my_perf_selection = PerfSelection::SelectionSh;
        }
        self.form.perf_selection_valid();
        self.form.my_glued_f.clear();
        self.form.glued_faces_valid();

        // `cxx:872-877`.
        let pbase = self
            .my_pbase
            .clone()
            .ok_or_else(|| "StdFail_NotDone".to_string())?;
        let height = statics::height_max(
            self.form.my_sbase.as_ref(),
            self.form.my_skface.as_ref(),
            self.form.my_sfrom.as_ref(),
            self.form.my_suntil.as_ref(),
        );
        let v = self.scaled_dir(3.0 * height);
        let vtra = self.scaled_dir(-3.0 * height / 2.0);
        let the_prism = Prism::with_base_and_translation(&pbase, &v, &vtra);
        let vraiprism = the_prism
            .shape()
            .cloned()
            .ok_or_else(|| "StdFail_NotDone".to_string())?;

        // `cxx:878-880`.
        statics::maj_map(
            &pbase,
            &the_prism,
            &mut self.form.my_map,
            &mut self.form.my_f_shape,
            &mut self.form.my_l_shape,
        );
        self.form.my_gshape = Some(vraiprism);
        self.form.generated_shape_valid();

        // `cxx:882-896`.
        if !self.form.my_fuse {
            if let Some(sbase) = self.form.my_sbase.clone() {
                if let Some(gshape) = self.form.my_gshape.clone() {
                    let mut c = cut(&sbase, &gshape);
                    if c.build().is_ok() {
                        if let Some(shape) = c.shape().cloned() {
                            self.form.my_shape = Some(shape.clone());
                            update_descendants_bop(&mut self.form, &c, &shape, false)?;
                            self.form.done();
                        }
                    }
                }
            }
        } else {
            self.take_prism_curves(&the_prism);
            self.global_perform()?;
        }
        Ok(())
    }

    /// `myDir.Reverse()` (`cxx:719`), kept next to its only call site.
    fn reverse_dir(&mut self) {
        let mut d = *self.dir();
        d.reverse();
        self.set_dir(d);
    }
}
