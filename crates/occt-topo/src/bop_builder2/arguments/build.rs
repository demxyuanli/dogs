use super::super::prelude::*;
use super::super::*;

impl BopBuilder {
    // -----------------------------------------------------------------------
    // Building the result
    // -----------------------------------------------------------------------

    /// Records object/tool groups and fills GF images.
    ///
    /// This is **not** `BOPAlgo_Builder::BuildBOP` (that is
    /// [`BopBuilder::build_bop`]) but the FillImages sequence of
    /// `PerformInternal1`. IN/OUT is validated here only so a direct caller
    /// cannot bind an illegal BOP state; FillImages itself does not read the
    /// states. Fuse/Cut/Common filtering is `crate::bop_bop::build_shape`.
    pub fn fill_images_for_groups(
        &mut self,
        objects: &[TopoShape],
        obj_state: FaceState,
        tools: &[TopoShape],
        tools_state: FaceState,
    ) -> Result<(), String> {
        if self.has_errors() {
            return Err(self.errors.first().cloned().unwrap_or_default());
        }
        if !matches!(obj_state, FaceState::In | FaceState::Out)
            || !matches!(tools_state, FaceState::In | FaceState::Out)
        {
            let msg = "BOPAlgo_Builder: invalid state for the operation".to_string();
            self.add_error(msg.clone());
            return Err(msg);
        }
        if objects.is_empty() && tools.is_empty() {
            let msg = "BOPAlgo_Builder: too few arguments".to_string();
            self.add_error(msg.clone());
            return Err(msg);
        }
        for s in objects.iter().chain(tools.iter()) {
            if self.filler.ds().index(s).is_none() {
                let msg = "BOPAlgo_Builder: unknown shape for the operation".to_string();
                self.add_error(msg.clone());
                return Err(msg);
            }
        }
        self.objects = objects.to_vec();
        self.tools = tools.to_vec();
        self.obj_state = obj_state;
        self.tools_state = tools_state;
        self.fill_images()
    }

    /// Builds the result of the BOP from the faces of the two groups selected
    /// by their IN/OUT states.
    ///
    /// Port of `BOPAlgo_Builder::BuildBOP(theObjects, theObjState, theTools,
    /// theToolsState, ...)` (`BOPAlgo_Builder.cxx:483-841`): input validation,
    /// face selection ([`crate::bop_build_bop::build_bop`]), unused-face shell
    /// assembly, internal parts, result compound and `PrepareHistory`.
    ///
    /// UNPORTED: `myShape.IsNull()` (`cxx:496`) — [`TopoShape::is_null`] is
    /// always false in the port, so the `BOPAlgo_AlertBuilderFailed` arm can
    /// only be reached through `myArguments.IsEmpty()`. The
    /// `SetContext(myContext)` / `SetRunParallel` calls on `BOPAlgo_BuilderSolid`
    /// (`cxx:723-724`) have no counterpart: the port's `BuilderSolid` builds its
    /// own `IntToolsContext`.
    pub fn build_bop(
        &mut self,
        objects: &[TopoShape],
        obj_state: FaceState,
        tools: &[TopoShape],
        tools_state: FaceState,
    ) -> Result<(), String> {
        // `cxx:489-490`: a failed Builder leaves the result untouched.
        if self.has_errors() {
            return Err(self.errors.first().cloned().unwrap_or_default());
        }
        // `cxx:496-501`.
        if self.arguments.is_empty() {
            let msg = "BOPAlgo_AlertBuilderFailed".to_string();
            self.add_error(msg.clone());
            return Err(msg);
        }
        // `cxx:502-508`.
        if !matches!(obj_state, FaceState::In | FaceState::Out)
            || !matches!(tools_state, FaceState::In | FaceState::Out)
        {
            let msg = "BOPAlgo_AlertBOPNotSet".to_string();
            self.add_error(msg.clone());
            return Err(msg);
        }
        // `cxx:510-517`.
        if objects.is_empty() && tools.is_empty() {
            let msg = "BOPAlgo_AlertTooFewArguments".to_string();
            self.add_error(msg.clone());
            return Err(msg);
        }
        // `cxx:518-552`: every given shape must belong to the arguments and be
        // a SOLID or a compound/comp-solid made of solids only.
        for group in [objects, tools] {
            for s in group {
                if self.filler.ds().index(s).is_none() {
                    let msg = "BOPAlgo_AlertUnknownShape".to_string();
                    self.add_error(msg.clone());
                    return Err(msg);
                }
                if s.shape_type() == ShapeType::Solid {
                    continue;
                }
                let mut treated: Vec<TopoShape> = Vec::new();
                let mut fence: HashSet<ShapeMapKey> = HashSet::new();
                crate::bop_occt_util::treat_compound(s, &mut treated, &mut fence);
                for sx in &treated {
                    if sx.shape_type() != ShapeType::Solid
                        && sx.shape_type() != ShapeType::CompSolid
                    {
                        let msg = "BOPAlgo_AlertUnsupportedType".to_string();
                        self.add_error(msg.clone());
                        return Err(msg);
                    }
                }
            }
        }

        self.objects = objects.to_vec();
        self.tools = tools.to_vec();
        self.obj_state = obj_state;
        self.tools_state = tools_state;
        crate::bop_build_bop::build_bop(self)?;
        // `cxx:840`: the result history is rebuilt from the new result shape.
        self.prepare_history();
        Ok(())
    }

    /// Builds the result of the Boolean operation `op` from the two groups.
    ///
    /// Port of the inline `BOPAlgo_Builder::BuildBOP(theObjects, theTools,
    /// theOperation, ...)` overload (`BOPAlgo_Builder.hxx:213-247`): the
    /// operation type is converted into the object/tool states and forwarded to
    /// [`BopBuilder::build_bop`].
    ///
    /// UNPORTED: the `BOPAlgo_CUT21` (`TopAbs_IN`, `TopAbs_OUT`) and the
    /// `default` (`TopAbs_UNKNOWN`, `TopAbs_UNKNOWN`) arms — [`BoolOp2`] has no
    /// `Cut21`/`Unknown` variant (`crates/occt-topo/src/bop_builder2/data_access.rs:105`).
    pub fn build_bop_operation(
        &mut self,
        objects: &[TopoShape],
        tools: &[TopoShape],
        op: BoolOp2,
    ) -> Result<(), String> {
        let (obj_state, tools_state) = op.states();
        self.build_bop(objects, obj_state, tools, tools_state)
    }

    /// Assembles the result compound from the images of the arguments.
    ///
    /// Port of `BOPAlgo_Builder::BuildResult` (`Builder_1.cxx:130`) and
    /// `BOPAlgo_BOP::BuildResult` (`BOP.cxx:323`): for every argument of a
    /// given type, add its images, or the argument itself if it has none.
    /// Fence by `TShape` identity. No object/tool IN/OUT skip — Cut/Common
    /// membership is `BOPAlgo_BOP::BuildShape` / `BuildRC`.
    pub(crate) fn build_result(&mut self) -> Result<TopoShape, String> {
        let b = TopoBuilder::new();
        let mut result = b.make_compound_of(&[]);
        let mut fence: Vec<TopoShape> = Vec::new();
        for t in RESULT_TYPES {
            for arg in &self.arguments {
                if arg.shape_type() != t {
                    continue;
                }
                let imgs: Vec<TopoShape> = match self.history.image(arg) {
                    Some(list) => list.to_vec(),
                    None => vec![arg.clone()],
                };
                for img in imgs {
                    if !fence.iter().any(|f| is_same(f, &img)) {
                        fence.push(img.clone());
                        b.add(&mut result.0, &img);
                    }
                }
            }
        }
        self.result_shape = result.0.clone();
        Ok(result.0)
    }
}
