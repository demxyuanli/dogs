//! `BRepFeat_MakeDPrism::PerformUntilEnd` (`cxx:620-648`),
//! `PerformFromEnd` (`cxx:655-816`) and `PerformThruAll` (`cxx:824-875`).
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
use crate::builder::TopoBuilder;
use crate::locope::{CSIntersector, DPrism};
use crate::shape::{Face, TopoShape};
use crate::topo_tools::{faces_of, is_same};

use super::statics;
use super::MakeDPrism;

impl MakeDPrism {
    /// `PerformUntilEnd()` (`cxx:620-648`): a semi-infinite draft prism limited
    /// only by the basis, so the sweep is the two-height extrusion of the
    /// bounding-box extent.
    pub fn perform_until_end(&mut self) -> Result<(), String> {
        // `cxx:626-634`.
        self.form.my_perf_selection = PerfSelection::SelectionSh;
        self.form.perf_selection_valid();
        self.form.my_glued_f.clear();
        self.form.my_suntil = None;
        self.form.shape_until_valid();
        self.form.my_sfrom = None;
        self.form.shape_from_valid();

        // `cxx:634-647`.
        let pbase = self.pbase_face()?;
        let height = statics::height_max(
            self.form.my_sbase.as_ref(),
            self.form.my_skface.as_ref(),
            self.form.my_sfrom.as_ref(),
            self.form.my_suntil.as_ref(),
        );
        // `cxx:637-638`: the single-height constructor.
        let the_dprism = DPrism::with_height(&pbase, height, self.angle());
        let vraiprism = the_dprism
            .shape()
            .cloned()
            .ok_or_else(|| "StdFail_NotDone".to_string())?;

        statics::maj_map(
            &pbase.0,
            &the_dprism,
            &mut self.form.my_map,
            &mut self.form.my_f_shape,
            &mut self.form.my_l_shape,
        );
        self.form.my_gshape = Some(vraiprism);
        self.form.generated_shape_valid();
        self.form.glued_faces_valid();
        self.take_dprism_curves(&the_dprism);
        self.global_perform()
    }

    /// `PerformFromEnd(Until)` (`cxx:655-816`): a semi-infinite draft prism
    /// running away from the basis and cut at `Until`.
    pub fn perform_from_end(&mut self, until: &TopoShape) -> Result<(), String> {
        // `cxx:661-672`.
        if until.is_null() {
            return Err("Standard_ConstructionError".to_string());
        }
        if let Some(sk) = self.form.my_skface.clone() {
            if is_same(until, &sk.0) {
                return self.perform_until_end();
            }
        }
        if faces_of(until).is_empty() {
            return Err("Standard_ConstructionError".to_string());
        }

        // `cxx:678-686`.
        self.form.my_perf_selection = PerfSelection::SelectionShU;
        self.form.perf_selection_valid();
        self.form.my_sfrom = None;
        self.form.shape_from_valid();
        self.form.my_suntil = Some(until.clone());
        let trf = self.form.transform_shape_fu(1);
        self.form.shape_until_valid();

        // `cxx:687-692`.
        let pbase = self.pbase_face()?;
        let c = statics::test_curve(&pbase);
        let sens = match &c {
            Some(c) => statics::sens_of_prism(c, until),
            // `TestCurve` answers a null handle on a non-planar profile; the
            // direction is undefined there, and the reference would raise.
            None => 1,
        };
        let height = statics::height_max(
            self.form.my_sbase.as_ref(),
            self.form.my_skface.as_ref(),
            self.form.my_sfrom.as_ref(),
            self.form.my_suntil.as_ref(),
        );
        let the_dprism =
            DPrism::with_two_heights(&pbase, sens as f64 * height, sens as f64 * height, self.angle());
        let vraiprism = the_dprism
            .shape()
            .cloned()
            .ok_or_else(|| "StdFail_NotDone".to_string())?;
        // `cxx:692-696`: `Shape()` already raised on an unfinished prism, so the
        // null check answers the empty result of a done-but-null `myRes`.
        if vraiprism.is_null() {
            self.set_error(StatusError::NullRealTool);
            return Ok(());
        }

        // `cxx:699` and `cxx:711` both start with the same `MajMap`.
        statics::maj_map(
            &pbase.0,
            &the_dprism,
            &mut self.form.my_map,
            &mut self.form.my_f_shape,
            &mut self.form.my_l_shape,
        );

        if !trf {
            // `cxx:698-708`: the limit is a face of the basis.
            self.form.my_gshape = Some(vraiprism);
            self.form.generated_shape_valid();
            self.form.my_glued_f.clear();
            self.form.glued_faces_valid();
            self.take_dprism_curves(&the_dprism);
            self.global_perform()?;
        } else {
            // `cxx:709-815`: the support was trimmed, so both the limit and the
            // basis become tool solids.
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

            // `cxx:727-738`. Unlike the other overloads, a miss against `Until`
            // does not stop the computation; the orientation stays FORWARD.
            let mut or_u = Orientation::Forward;
            let mut f_until: Option<Face> = None;
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
                let prm = asi1.point(1, 1).map(|p| p.parameter()).unwrap_or(0.0);
                if prm < 0.0 {
                    or_u = or_u.reversed();
                }
                f_until = asi1.point(1, 1).ok().map(|p| p.face().clone());
            }

            // `cxx:740-766`: the basis point with the smallest parameter wins.
            let mut or_f = Orientation::Forward;
            let mut f_from: Option<Face> = None;
            let nb2 = if asi2.is_done() {
                asi2.nb_points(1).unwrap_or(0)
            } else {
                0
            };
            if nb2 >= 1 {
                let mut prm = asi2.point(1, 1).map(|p| p.parameter()).unwrap_or(0.0);
                f_from = asi2.point(1, 1).ok().map(|p| p.face().clone());
                or_f = asi2
                    .point(1, 1)
                    .map(|p| p.orientation().reversed())
                    .unwrap_or(Orientation::Forward);
                for iii in 1..=nb2 {
                    if asi2.point(1, iii).map(|p| p.parameter()).unwrap_or(0.0) < prm {
                        prm = asi2.point(1, iii).map(|p| p.parameter()).unwrap_or(0.0);
                        f_from = asi2.point(1, iii).ok().map(|p| p.face().clone());
                        or_f = asi2
                            .point(1, iii)
                            .map(|p| p.orientation().reversed())
                            .unwrap_or(Orientation::Forward);
                    }
                }
                // `cxx:758-765`: the support face is replaced by an unbounded
                // face on its basis surface.
                if let Some(ff) = &f_from {
                    if let Some(s) = BRepTool::face_surface(ff) {
                        let s = s.rectangular_trimmed_basis().unwrap_or(s);
                        let fac = brep_lib_make_face::make_face_from_surface(s, PCONFUSION);
                        self.form.my_sfrom = Some(fac.0.clone());
                        // `Trf` is reassigned by the reference here and never
                        // read again (`cxx:764-765`).
                        let _ = self.form.transform_shape_fu(0);
                    }
                }
            }

            // `cxx:768-792`.
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

            // `cxx:794-815`.
            let comp = TopoBuilder::new().make_compound_of(&tools);
            let (_tr_p, tr_shape) = statics::boolean_result(cut(&vraiprism, &comp.0))?;
            self.finish_with_boolean_on(&tr_shape)?;
        }
        Ok(())
    }

    /// `PerformThruAll()` (`cxx:824-875`): the draft prism is swept through the
    /// whole bounding box, two heights in both directions around the profile.
    pub fn perform_thru_all(&mut self) -> Result<(), String> {
        // `cxx:831-844`.
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

        // `cxx:848-853`.
        let pbase = self.pbase_face()?;
        let height = statics::height_max(
            self.form.my_sbase.as_ref(),
            self.form.my_skface.as_ref(),
            self.form.my_sfrom.as_ref(),
            self.form.my_suntil.as_ref(),
        );
        let the_dprism = DPrism::with_two_heights(&pbase, height, height, self.angle());
        let vraiprism = the_dprism
            .shape()
            .cloned()
            .ok_or_else(|| "StdFail_NotDone".to_string())?;

        statics::maj_map(
            &pbase.0,
            &the_dprism,
            &mut self.form.my_map,
            &mut self.form.my_f_shape,
            &mut self.form.my_l_shape,
        );
        self.form.my_gshape = Some(vraiprism);
        self.form.generated_shape_valid();

        // `cxx:857-874`.
        if !self.form.my_fuse {
            if let Some(sbase) = self.form.my_sbase.clone() {
                if let Some(gshape) = self.form.my_gshape.clone() {
                    // `cxx:858-865`: the result is kept only when the Boolean
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
            self.take_dprism_curves(&the_dprism);
            self.global_perform()?;
        }
        Ok(())
    }
}
