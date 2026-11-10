use super::super::prelude::*;
use super::super::*;

impl BopBuilder {
    // -----------------------------------------------------------------------
    // Performing the operation
    // -----------------------------------------------------------------------

    /// Performs the General Fuse operation: runs the intersection phase
    /// (`PaveFiller::perform`) followed by FillImages + BuildResult + post-treat.
    ///
    /// The result of the General Fuse is a compound containing all split parts
    /// of the arguments. Object/tool IN/OUT filtering is not applied here
    /// (`BOPAlgo_Builder::PerformInternal1`).
    pub fn perform(&mut self) -> Result<TopoShape, String> {
        let args = self.arguments.clone();
        self.perform_internal(&args, FaceState::Out, &[], FaceState::Out)
    }

    /// `BOPAlgo_BOP::Perform` (`BOPAlgo_BOP.cxx:364-409`): the intersection
    /// phase receives `objects ∪ tools`, then `PerformInternal1` runs the
    /// Boolean building phase with the given IN/OUT states.
    ///
    /// The Rust [`Self::perform`] is `BOPAlgo_Builder::Perform` (General Fuse,
    /// no tool group); this is the `BOPAlgo_BOP` override that
    /// `BRepFeat_Builder::Perform` inherits.
    pub fn perform_bop(
        &mut self,
        objects: &[TopoShape],
        obj_state: FaceState,
        tools: &[TopoShape],
        tools_state: FaceState,
    ) -> Result<TopoShape, String> {
        self.perform_internal(objects, obj_state, tools, tools_state)
    }

    /// `BOPAlgo_Builder::PerformWithFiller` (`BOPAlgo_Builder.cxx:198-208`):
    /// runs the building phase on top of the filler the caller already
    /// performed.
    ///
    /// The option copy of the reference (`myNonDestructive` / `myFuzzyValue` /
    /// `myGlue` / `myUseOBB` taken from the filler) is free here: the port keeps
    /// those options on the [`PaveFiller`] itself, which is exactly where
    /// [`Self::non_destructive`] and [`Self::fuzzy_value`] read them.
    pub fn perform_with_filler(&mut self) -> Result<TopoShape, String> {
        self.errors.clear();
        self.warnings.clear();
        // `CheckData`: no argument means nothing to build.
        if self.arguments.is_empty() {
            let msg = "BOPAlgo_Builder: too few arguments".to_string();
            self.add_error(msg.clone());
            return Err(msg);
        }
        if self.filler.has_errors() {
            let msgs: Vec<String> = self.filler.errors().to_vec();
            for m in msgs {
                self.add_error(m);
            }
            return Err(self.errors.first().cloned().unwrap_or_default());
        }
        // `PerformInternal1`: FillImages → BuildResult → PrepareHistory →
        // PostTreat. `myShape` starts as an empty compound (`Prepare`).
        self.reset_run();
        self.fill_images()?;
        self.build_result()?;
        self.prepare_history();
        self.post_treat()?;
        Ok(self.result_shape.clone())
    }

    /// Shared pipeline of [`BopBuilder::perform`] and `BOPAlgo_BOP`.
    ///
    /// Mirrors `PerformInternal1`: FillImages (GF, no obj/tool filter) then
    /// [`Self::build_result`] (GF assembly of every argument image). When tools
    /// are present this is a Boolean run, so `BOPAlgo_BOP::BuildShape` applies
    /// Fuse/Cut/Common afterwards (`BOPAlgo_BOP.cxx:563`).
    pub(crate) fn perform_internal(
        &mut self,
        objects: &[TopoShape],
        obj_state: FaceState,
        tools: &[TopoShape],
        tools_state: FaceState,
    ) -> Result<TopoShape, String> {
        self.reset_run();
        if self.arguments.len() < 2 {
            let msg = "BOPAlgo_Builder: too few arguments".to_string();
            self.add_error(msg.clone());
            return Err(msg);
        }

        // Intersection phase.
        self.filler.set_arguments(&self.arguments);
        if let Err(e) = self.filler.perform() {
            let msgs: Vec<String> = self.filler.errors().to_vec();
            for m in msgs {
                self.add_error(m);
            }
            if self.errors.is_empty() {
                self.add_error(e.clone());
            }
            return Err(self.errors.first().cloned().unwrap_or_else(|| e));
        }
        if self.filler.has_errors() {
            let msgs: Vec<String> = self.filler.errors().to_vec();
            for m in msgs {
                self.add_error(m);
            }
            return Err(self.errors.first().cloned().unwrap_or_default());
        }

        // GF building phase: record the BOP groups for BuildShape, then fill
        // images without IN/OUT selection (`FillImages*` does not read states).
        self.objects = objects.to_vec();
        self.tools = tools.to_vec();
        self.obj_state = obj_state;
        self.tools_state = tools_state;
        self.fill_images()?;
        self.build_result()?;
        if !self.tools.is_empty() {
            crate::bop_bop::build_shape(self)?;
        }
        self.prepare_history();
        self.post_treat()?;
        Ok(self.result_shape.clone())
    }

    /// Resets the per-run state (report, history, origins, same-domain map,
    /// result). The arguments and the filler options survive.
    pub(crate) fn reset_run(&mut self) {
        self.errors.clear();
        self.warnings.clear();
        self.history.clear();
        self.origins.clear();
        self.shapes_sd.clear();
        self.objects.clear();
        self.tools.clear();
        self.obj_state = FaceState::Out;
        self.tools_state = FaceState::Out;
        self.result_shape = TopoShape::new(ShapeType::Compound);
        self.in_parts.clear();
    }

    /// Clears the content of the algorithm: the report, the history, the
    /// arguments and the filler (`BOPAlgo_Builder::Clear`).
    pub fn clear(&mut self) {
        self.errors.clear();
        self.warnings.clear();
        self.history.clear();
        self.origins.clear();
        self.shapes_sd.clear();
        self.arguments.clear();
        self.objects.clear();
        self.tools.clear();
        self.result_shape = TopoShape::new(ShapeType::Compound);
        self.in_parts.clear();
    }
}
