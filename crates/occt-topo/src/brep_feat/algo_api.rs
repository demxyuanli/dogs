//! `BRepAlgoAPI_BuilderAlgo` / `BRepAlgoAPI_BooleanOperation` with history.
//!
//! Source: `ModelingAlgorithms/TKBO/BRepAlgoAPI/BRepAlgoAPI_BuilderAlgo.cxx:48-333`,
//! `BRepAlgoAPI_BooleanOperation.cxx:34-200`, `BRepAlgoAPI_Fuse.cxx`,
//! `BRepAlgoAPI_Cut.cxx`.
//!
//! The `BRepFeat` feature classes keep the Boolean operation alive after
//! `Build()` because `UpdateDescendants` queries its history
//! (`BRepFeat_Form.cxx:470`, `BRepFeat_MakePrism.cxx:796`), so this module is a
//! facade over [`crate::bop_builder2`] that answers `Shape` / `Modified` /
//! `IsDeleted` from the run's [`crate::bop_builder2::BopBuilder`] instead of
//! throwing the algorithm away.
//!
//! Deviations, all forced by the language:
//! - the constructors are [`BooleanOperation::new`] (operation only, the two
//!   groups follow through `SetArguments` / `SetTools`) and the free functions
//!   [`cut`] / [`fuse`], which stand for the two-shape constructors of
//!   `BRepAlgoAPI_Cut` / `_Fuse`.
//! - `Build()` raises in the reference (`StdFail_NotDone` through `Shape()`);
//!   here the failure comes back as the `Err` of [`BooleanOperation::build`],
//!   and `HasErrors` stays observable through [`BooleanOperation::has_errors`].

use crate::bop_builder2::{builder_bop_algo_with_history, BoolOp2, BopBuilder};
use crate::shape::TopoShape;

/// `BRepAlgoAPI_BooleanOperation` (`BRepAlgoAPI_BooleanOperation.hxx`).
#[derive(Debug, Clone)]
pub struct BooleanOperation {
    /// `myObjects` - the object group (`SetArguments`).
    objects: Vec<TopoShape>,
    /// `myTools` - the tool group (`SetTools`).
    tools: Vec<TopoShape>,
    /// `myOperation` (`BRepAlgoAPI_BooleanOperation::SetOperation`).
    op: BoolOp2,
    /// The result of the last successful `Build()` (`myShape`).
    result: Option<TopoShape>,
    /// The performed `BOPAlgo_BOP`, kept for the history queries. The reference
    /// keeps `myBuilder` (`BRepAlgoAPI_BuilderAlgo.hxx:157`) for the same
    /// reason.
    run: Option<BopBuilder>,
    /// `myReport` errors of the last run.
    errors: Vec<String>,
}

impl BooleanOperation {
    /// `BRepAlgoAPI_BooleanOperation()` plus `SetOperation`: the reference
    /// builds the operation on the `BOPAlgo_BOP` that the caller configured.
    pub fn new(op: BoolOp2) -> Self {
        Self {
            objects: Vec::new(),
            tools: Vec::new(),
            op,
            result: None,
            run: None,
            errors: Vec::new(),
        }
    }

    /// `SetOperation` (`BRepAlgoAPI_BooleanOperation::SetOperation`).
    pub fn set_operation(&mut self, op: BoolOp2) {
        self.op = op;
    }

    /// `SetArguments` (`BRepAlgoAPI_BuilderAlgo::SetArguments`): the object
    /// group.
    pub fn set_arguments(&mut self, objects: &[TopoShape]) {
        self.objects = objects.to_vec();
    }

    /// `SetTools` (`BRepAlgoAPI_BooleanOperation::SetTools`): the tool group.
    pub fn set_tools(&mut self, tools: &[TopoShape]) {
        self.tools = tools.to_vec();
    }

    /// `Build()` (`BRepAlgoAPI_BuilderAlgo::Build` /
    /// `BRepAlgoAPI_BooleanOperation::Build`).
    ///
    /// The reference raises `StdFail_NotDone` from `Shape()` when the run
    /// reports errors; the port answers the first error instead, and leaves
    /// `myShape` unset exactly as the reference does.
    pub fn build(&mut self) -> Result<(), String> {
        self.errors.clear();
        self.result = None;
        // `BRepAlgoAPI_BuilderAlgo::BuildResult`: `SetToFillHistory(true)` on
        // the builder, then the run. `BOPAlgo_BOP::Perform` is the inherited
        // perform that the two-shape Booleans use.
        let run = builder_bop_algo_with_history(&self.objects, &self.tools, self.op, true);
        if run.has_errors() {
            self.errors = run.errors().to_vec();
            self.run = Some(run);
            return Err(self
                .errors
                .first()
                .cloned()
                .unwrap_or_else(|| "StdFail_NotDone".to_string()));
        }
        self.result = Some(run.result().clone());
        self.run = Some(run);
        Ok(())
    }

    /// `Shape()`: the result of the last successful `Build()`.
    pub fn shape(&self) -> Option<&TopoShape> {
        self.result.as_ref()
    }

    /// `Modified(S)` (`BRepAlgoAPI_BuilderAlgo::Modified`): the images of `S`
    /// in the result, from the run's history.
    pub fn modified(&self, s: &TopoShape) -> Vec<TopoShape> {
        match &self.run {
            Some(run) => run
                .history()
                .modified(s)
                .map(|v| v.to_vec())
                .unwrap_or_default(),
            None => Vec::new(),
        }
    }

    /// `IsDeleted(S)` (`BRepAlgoAPI_BuilderAlgo::IsDeleted`).
    pub fn is_deleted(&self, s: &TopoShape) -> bool {
        match &self.run {
            Some(run) => run.history().is_deleted(s),
            None => false,
        }
    }

    /// `HasErrors` of the last run.
    pub fn has_errors(&self) -> bool {
        !self.errors.is_empty()
    }

    /// The errors of the last run.
    pub fn errors(&self) -> &[String] {
        &self.errors
    }
}

/// `BRepAlgoAPI_Cut(object, tool)` (`BRepAlgoAPI_Cut.cxx`).
pub fn cut(object: &TopoShape, tool: &TopoShape) -> BooleanOperation {
    let mut bop = BooleanOperation::new(BoolOp2::Cut);
    bop.set_arguments(&[object.clone()]);
    bop.set_tools(&[tool.clone()]);
    bop
}

/// `BRepAlgoAPI_Fuse(object, tool)` (`BRepAlgoAPI_Fuse.cxx`).
pub fn fuse(object: &TopoShape, tool: &TopoShape) -> BooleanOperation {
    let mut bop = BooleanOperation::new(BoolOp2::Fuse);
    bop.set_arguments(&[object.clone()]);
    bop.set_tools(&[tool.clone()]);
    bop
}

/// `BRepAlgoAPI_Common(object, tool)` (`BRepAlgoAPI_Common.cxx`).
pub fn common(object: &TopoShape, tool: &TopoShape) -> BooleanOperation {
    let mut bop = BooleanOperation::new(BoolOp2::Common);
    bop.set_arguments(&[object.clone()]);
    bop.set_tools(&[tool.clone()]);
    bop
}
