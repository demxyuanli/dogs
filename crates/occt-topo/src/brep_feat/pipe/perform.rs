//! The three `BRepFeat_MakePipe::Perform` overloads (`cxx:177-336`).
//!
//! All three build the same primitive - [`crate::locope::Pipe`] over `myPbase`
//! along `mySpine` - record its first/last wire through `MajMap`, and then
//! split: `Perform()` fuses or cuts the primitive directly when no glued face
//! exists, while the two limited overloads always go through `GlobalPerform`.

use crate::brep_feat::algo_api::{cut, fuse};
use crate::brep_feat::form::update_descendants_bop;
use crate::brep_feat::perf_selection::PerfSelection;
use crate::locope::{sample_edges, Pipe};
use crate::shape::TopoShape;
use crate::topo_tools::{faces_of, is_same};

use super::statics;
use super::MakePipe;

impl MakePipe {
    /// `Perform()` (`cxx:177-228`).
    pub fn perform(&mut self) -> Result<(), String> {
        // `cxx:182-188`.
        self.form.my_sfrom = None;
        self.form.shape_from_valid();
        self.form.my_suntil = None;
        self.form.shape_until_valid();
        self.form.my_glued_f.clear();
        self.form.my_perf_selection = PerfSelection::NoSelection;
        self.form.perf_selection_valid();

        // `cxx:189-194`.
        let mut the_pipe = self.sweep()?;
        self.form.glued_faces_valid();

        if self.form.my_glued_f.is_empty() {
            // `cxx:196-218`.
            let sbase = self
                .form
                .my_sbase
                .clone()
                .ok_or_else(|| "StdFail_NotDone".to_string())?;
            let gshape = self
                .form
                .my_gshape
                .clone()
                .ok_or_else(|| "StdFail_NotDone".to_string())?;
            if self.form.my_fuse {
                // `cxx:198-204`.
                let (f, shape) = statics::boolean_result(fuse(&sbase, &gshape))?;
                self.form.my_shape = Some(shape.clone());
                update_descendants_bop(&mut self.form, &f, &shape, false)?;
                self.form.done();
            } else {
                // `cxx:205-217`. The reference's third arm (`myFuse == 1` false
                // and `myFuse == 0` false) cannot be reached by a `bool`, and
                // `myFuse == 0` is the cut arm.
                let (c, shape) = statics::boolean_result(cut(&sbase, &gshape))?;
                self.form.my_shape = Some(shape.clone());
                update_descendants_bop(&mut self.form, &c, &shape, false)?;
                self.form.done();
            }
        } else {
            // `cxx:219-227`.
            self.curves_from_pipe(&mut the_pipe)?;
            self.global_perform()?;
        }
        Ok(())
    }

    /// `Perform(Until)` (`cxx:234-270`).
    pub fn perform_until(&mut self, until: &TopoShape) -> Result<(), String> {
        // `cxx:239-247`.
        if until.is_null() {
            return Err("Standard_ConstructionError".to_string());
        }
        if faces_of(until).is_empty() {
            return Err("Standard_ConstructionError".to_string());
        }

        // `cxx:248-255`.
        self.form.my_glued_f.clear();
        self.form.my_perf_selection = PerfSelection::SelectionU;
        self.form.perf_selection_valid();
        self.form.my_sfrom = None;
        self.form.shape_from_valid();
        self.form.my_suntil = Some(until.clone());
        let _ = self.form.transform_shape_fu(1);
        self.form.shape_until_valid();

        // `cxx:256-268`.
        let mut the_pipe = self.sweep()?;
        self.form.glued_faces_valid();
        self.curves_from_pipe(&mut the_pipe)?;

        // `cxx:269`.
        self.global_perform()
    }

    /// `Perform(From, Until)` (`cxx:276-336`).
    pub fn perform_from_until(&mut self, from: &TopoShape, until: &TopoShape) -> Result<(), String> {
        // `cxx:281-284`.
        if from.is_null() || until.is_null() {
            return Err("Standard_ConstructionError".to_string());
        }
        // `cxx:285-295`: a limit equal to the sketch face degenerates to the
        // one-limit overload.
        if let Some(sk) = self.my_skface.clone() {
            if is_same(from, &sk.0) {
                return self.perform_until(until);
            } else if is_same(until, &sk.0) {
                return self.perform_until(from);
            }
        }

        // `cxx:296-317`.
        self.form.my_glued_f.clear();
        self.form.my_perf_selection = PerfSelection::SelectionFU;
        self.form.perf_selection_valid();
        if faces_of(from).is_empty() {
            return Err("Standard_ConstructionError".to_string());
        }
        if faces_of(until).is_empty() {
            return Err("Standard_ConstructionError".to_string());
        }
        self.form.my_sfrom = Some(from.clone());
        let _ = self.form.transform_shape_fu(0);
        self.form.shape_from_valid();
        self.form.my_suntil = Some(until.clone());
        let _ = self.form.transform_shape_fu(1);
        self.form.shape_until_valid();

        // `cxx:318-330`.
        let mut the_pipe = self.sweep()?;
        self.form.glued_faces_valid();
        self.curves_from_pipe(&mut the_pipe)?;

        // `cxx:335`.
        self.global_perform()
    }

    /// `cxx:189-194`, `cxx:256-259` and `cxx:318-321`: build the pipe, run
    /// `MajMap` and publish the primitive as `myGShape`.
    fn sweep(&mut self) -> Result<Pipe, String> {
        let pbase = self
            .my_pbase
            .clone()
            .ok_or_else(|| "StdFail_NotDone".to_string())?;
        let spine = self
            .my_spine
            .clone()
            .ok_or_else(|| "StdFail_NotDone".to_string())?;

        // `LocOpe_Pipe thePipe(mySpine, theBase)`; `Shape()` raises
        // `StdFail_NotDone` when the sweep did not complete
        // (`LocOpe_Pipe.cxx:279-282`).
        let mut the_pipe = Pipe::new(&spine, &pbase)?;
        let vrai_pipe = the_pipe
            .shape()
            .cloned()
            .ok_or_else(|| "StdFail_NotDone".to_string())?;

        statics::maj_map(
            &pbase,
            &mut the_pipe,
            &mut self.form.my_map,
            &mut self.form.my_f_shape,
            &mut self.form.my_l_shape,
        )?;
        self.form.my_gshape = Some(vrai_pipe);
        self.form.generated_shape_valid();
        Ok(the_pipe)
    }

    /// `cxx:221-225` (and the same three lines in the two limited overloads):
    /// the top of the form is the pipe's first shape, and the axes are sampled
    /// from it.
    fn curves_from_pipe(&mut self, the_pipe: &mut Pipe) -> Result<(), String> {
        self.form.my_f_shape = the_pipe.first_shape().cloned();
        // `LocOpe::SampleEdges(myFShape, spt)` clears the sequence; a null
        // `myFShape` samples nothing.
        let spt = match self.form.my_f_shape.as_ref() {
            Some(shape) => sample_edges(shape),
            None => Vec::new(),
        };
        self.my_curves = the_pipe.curves(&spt)?.to_vec();
        self.my_bcurve = the_pipe.baryc_curve()?;
        Ok(())
    }
}
