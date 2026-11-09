//! `BRepFeat_MakePrism::Perform(Length)` (`cxx:210-319`) and
//! `PerformUntilHeight(Until, Length)` (`cxx:903-1005`).
//!
//! Both sweep the profile by a vector of known length: `Perform(Length)` takes
//! it directly, `PerformUntilHeight` multiplies the length by the direction
//! `SensOfPrism` finds. When the top of the prism is not tangent to the basis
//! they answer the bare primitive, otherwise the geometry goes through
//! `GlobalPerform` (gluer) or a plain Boolean on the primitive.

use std::sync::Arc;

use occt_geom::Curve;

use crate::abs::{Orientation, ShapeType};
use crate::brep_feat::algo_api::cut;
use crate::brep_feat::form::update_descendants_bop;
use crate::brep_feat::perf_selection::PerfSelection;
use crate::brep_feat::statics::{is_inside, tool};
use crate::builder::TopoBuilder;
use crate::locope::{CSIntersector, Prism};
use crate::shape::{Face, TopoShape};
use crate::topo_tools::{explorer_shapes, faces_of, is_same};

use super::statics;
use super::MakePrism;

impl MakePrism {
    /// `Perform(Length)` (`cxx:210-319`): an extrusion of `Length` along
    /// `myDir`, then the glue-or-Boolean decision (`cxx:284-318`).
    pub fn perform_length(&mut self, length: f64) -> Result<(), String> {
        // `cxx:218-224`.
        self.form.my_sfrom = None;
        self.form.shape_from_valid();
        self.form.my_suntil = None;
        self.form.shape_until_valid();
        self.form.my_glued_f.clear();
        self.form.my_perf_selection = PerfSelection::NoSelection;
        self.form.perf_selection_valid();

        // `cxx:225-228`.
        let v = self.scaled_dir(length);
        let pbase = self
            .my_pbase
            .clone()
            .ok_or_else(|| "StdFail_NotDone".to_string())?;
        let the_prism = Prism::with_base(&pbase, &v);
        let vraiprism = the_prism
            .shape()
            .cloned()
            .ok_or_else(|| "StdFail_NotDone".to_string())?;

        // `cxx:232`: `MajMap(myPbase, thePrism, myMap, myFShape, myLShape)`.
        statics::maj_map(
            &pbase,
            &the_prism,
            &mut self.form.my_map,
            &mut self.form.my_f_shape,
            &mut self.form.my_l_shape,
        );

        // `cxx:234-236`.
        self.form.my_gshape = Some(vraiprism.clone());
        self.form.generated_shape_valid();

        // `cxx:238-281`: the reference looks for the face of the prism whose
        // border wire is `myLShape` and then walks the basis faces testing
        // `ToFuse` / `IsInside`. Neither `FFace` nor `found` is read after the
        // block, so it has no effect on the result; it is kept branch-for-branch
        // because the reference's traversal is part of the control flow.
        let mut fface: Option<Face> = None;
        if self.form.my_skface.is_some() || !self.my_slface.is_empty() {
            if let Some(lshape) = self.form.my_l_shape.clone() {
                if lshape.shape_type() == ShapeType::Wire {
                    'faces: for f in explorer_shapes(&vraiprism, ShapeType::Face) {
                        for w in explorer_shapes(&f, ShapeType::Wire) {
                            if is_same(&w, &lshape) {
                                fface = Face::wrap(f.clone());
                                break 'faces;
                            }
                        }
                    }
                }
            }
            if let Some(sbase) = self.form.my_sbase.clone() {
                for ff in faces_of(&sbase) {
                    if !statics::to_fuse(Some(&ff), fface.as_ref()) {
                        continue;
                    }
                    // `cxx:277-280`. `ToFuse` is false on a null `FFace`, so the
                    // face is present here.
                    if let Some(fc) = fface.as_ref() {
                        if !is_same(&fc.0, &pbase) && is_inside(&ff, fc) {
                            break;
                        }
                    }
                }
            }
        }

        // `cxx:284-297`.
        self.form.glued_faces_valid();
        if !self.form.my_glued_f.is_empty() {
            self.form.my_just_gluer = true;
            self.take_prism_curves(&the_prism);
            self.global_perform()?;
        }

        // `cxx:297-318`.
        if !self.form.my_just_gluer {
            self.finish_with_boolean_on(&vraiprism)?;
        }
        Ok(())
    }

    /// `PerformUntilHeight(Until, Length)` (`cxx:903-1005`).
    ///
    /// The reference calls `Perform(Length)` on a null `Until` and
    /// `Perform(Until)` on a zero length *without returning*, so the selected
    /// overload runs and the body continues; that control flow is kept.
    pub fn perform_until_height(&mut self, until: &TopoShape, length: f64) -> Result<(), String> {
        if until.is_null() {
            self.perform_length(length)?;
        }
        if length == 0.0 {
            self.perform_until(until)?;
        }
        if faces_of(until).is_empty() {
            return Err("Standard_ConstructionError".to_string());
        }

        // `cxx:922-929`.
        self.form.my_glued_f.clear();
        self.form.my_perf_selection = PerfSelection::NoSelection;
        self.form.perf_selection_valid();
        self.form.my_sfrom = None;
        self.form.shape_from_valid();
        self.form.my_suntil = Some(until.clone());
        let trf = self.form.transform_shape_fu(1);
        self.form.shape_until_valid();

        // `cxx:930-933`.
        let pbase = self
            .my_pbase
            .clone()
            .ok_or_else(|| "StdFail_NotDone".to_string())?;
        let c = statics::test_curve(&pbase, self.dir());
        let sens = match &c {
            Some(c) => statics::sens_of_prism(c, until),
            // `TestCurve` always answers a line in the reference; the port only
            // loses it on an empty sample, where the parameter is undefined.
            None => 1,
        };
        let v = self.scaled_dir(sens as f64 * length);
        let the_prism = Prism::with_base(&pbase, &v);
        let vraiprism = the_prism
            .shape()
            .cloned()
            .ok_or_else(|| "StdFail_NotDone".to_string())?;

        // `cxx:937-946` / `cxx:948-1001`.
        statics::maj_map(
            &pbase,
            &the_prism,
            &mut self.form.my_map,
            &mut self.form.my_f_shape,
            &mut self.form.my_l_shape,
        );

        if !trf {
            self.form.my_gshape = Some(vraiprism);
            self.form.generated_shape_valid();
            self.form.glued_faces_valid();
            self.take_prism_curves(&the_prism);
            self.global_perform()?;
        } else {
            let scur: Vec<Option<Arc<dyn Curve>>> = match &c {
                Some(c) => vec![Some(c.clone())],
                None => vec![None],
            };
            let mut asi = CSIntersector::with_shape(until);
            let _ = asi.perform_curves(&scur);
            let nb = if asi.is_done() {
                asi.nb_points(1).unwrap_or(0)
            } else {
                0
            };
            if nb >= 1 {
                // `cxx:957-970`.
                let mut or = if self.form.my_fuse {
                    asi.point(1, 1).ok().map(|p| p.orientation())
                } else {
                    asi.point(1, nb).ok().map(|p| p.orientation())
                }
                .unwrap_or(Orientation::Forward);
                if sens == -1 {
                    or = or.reversed();
                }
                // `cxx:971-981`: the tool solid goes into a one-element compound.
                let mut tools: Vec<TopoShape> = Vec::new();
                if let Some(f_until) = asi.point(1, 1).ok().map(|p| p.face().clone()) {
                    if let Some(s) = tool(until, &f_until, or) {
                        tools.push(s.0);
                    }
                }
                let comp = TopoBuilder::new().make_compound_of(&tools);
                // `cxx:982-983`.
                let (tr_p, tr_shape) = statics::boolean_result(cut(&vraiprism, &comp.0))?;
                update_descendants_bop(&mut self.form, &tr_p, &tr_shape, false)?;
                // `cxx:984-1001`.
                self.finish_with_boolean_on(&tr_shape)?;
            }
        }
        Ok(())
    }
}
