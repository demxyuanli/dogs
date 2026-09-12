//! `BRepAlgoAPI` — Fuse / Cut / Common / BooleanOperation.
//!
//! Source: `BRepAlgoAPI_BooleanOperation.cxx`, `BRepAlgoAPI_Fuse.cxx`,
//! `BRepAlgoAPI_Cut.cxx`, `BRepAlgoAPI_Common.cxx`. Each class is a thin
//! wrapper over `BOPAlgo_BOP`: set objects + tools, set the operation, run.

use crate::bop_builder2::{builder_bop_with_fuzzy, BoolOp2};
use crate::shape::TopoShape;

/// `BRepAlgoAPI_BooleanOperation` — two-group Boolean with an explicit op.
#[derive(Debug, Clone)]
pub struct BooleanOperation {
    objects: Vec<TopoShape>,
    tools: Vec<TopoShape>,
    op: BoolOp2,
    fuzzy: f64,
    result: Option<TopoShape>,
    errors: Vec<String>,
}

impl BooleanOperation {
    /// Empty algorithm (`BRepAlgoAPI_BooleanOperation()`).
    pub fn new() -> Self {
        Self {
            objects: Vec::new(),
            tools: Vec::new(),
            op: BoolOp2::Fuse,
            fuzzy: 1e-7,
            result: None,
            errors: Vec::new(),
        }
    }

    /// Two operands and an operation (`BRepAlgoAPI_*` two-shape constructors).
    pub fn from_shapes(object: TopoShape, tool: TopoShape, op: BoolOp2) -> Self {
        Self {
            objects: vec![object],
            tools: vec![tool],
            op,
            fuzzy: 1e-7,
            result: None,
            errors: Vec::new(),
        }
    }

    /// `SetOperation`.
    pub fn set_operation(&mut self, op: BoolOp2) {
        self.op = op;
    }

    /// `SetFuzzyValue`.
    pub fn set_fuzzy_value(&mut self, fuzzy: f64) {
        self.fuzzy = fuzzy.max(1e-9);
    }

    /// `SetTools` / extra objects.
    pub fn set_tools(&mut self, tools: &[TopoShape]) {
        self.tools = tools.to_vec();
    }

    /// `SetArguments` for the object group.
    pub fn set_objects(&mut self, objects: &[TopoShape]) {
        self.objects = objects.to_vec();
    }

    /// `Build` / `Perform`.
    pub fn build(&mut self) -> Result<(), String> {
        self.errors.clear();
        match builder_bop_with_fuzzy(&self.objects, &self.tools, self.op, self.fuzzy) {
            Ok(s) => {
                self.result = Some(s);
                Ok(())
            }
            Err(e) => {
                self.errors.push(e.clone());
                Err(e)
            }
        }
    }

    /// `Shape()`.
    pub fn shape(&self) -> Option<&TopoShape> {
        self.result.as_ref()
    }

    /// `HasErrors`.
    pub fn has_errors(&self) -> bool {
        !self.errors.is_empty()
    }
}

impl Default for BooleanOperation {
    fn default() -> Self {
        Self::new()
    }
}

/// `BRepAlgoAPI_Fuse`.
pub fn fuse(object: &TopoShape, tool: &TopoShape) -> Result<TopoShape, String> {
    builder_bop_with_fuzzy(&[object.clone()], &[tool.clone()], BoolOp2::Fuse, 1e-7)
}

/// `BRepAlgoAPI_Cut`.
pub fn cut(object: &TopoShape, tool: &TopoShape) -> Result<TopoShape, String> {
    builder_bop_with_fuzzy(&[object.clone()], &[tool.clone()], BoolOp2::Cut, 1e-7)
}

/// `BRepAlgoAPI_Common`.
pub fn common(object: &TopoShape, tool: &TopoShape) -> Result<TopoShape, String> {
    builder_bop_with_fuzzy(&[object.clone()], &[tool.clone()], BoolOp2::Common, 1e-7)
}

/// `BRepAlgoAPI_Fuse` with an explicit fuzzy value.
pub fn fuse_with_fuzzy(object: &TopoShape, tool: &TopoShape, fuzzy: f64) -> Result<TopoShape, String> {
    builder_bop_with_fuzzy(&[object.clone()], &[tool.clone()], BoolOp2::Fuse, fuzzy)
}

/// `BRepAlgoAPI_Cut` with an explicit fuzzy value.
pub fn cut_with_fuzzy(object: &TopoShape, tool: &TopoShape, fuzzy: f64) -> Result<TopoShape, String> {
    builder_bop_with_fuzzy(&[object.clone()], &[tool.clone()], BoolOp2::Cut, fuzzy)
}

/// `BRepAlgoAPI_Common` with an explicit fuzzy value.
pub fn common_with_fuzzy(object: &TopoShape, tool: &TopoShape, fuzzy: f64) -> Result<TopoShape, String> {
    builder_bop_with_fuzzy(&[object.clone()], &[tool.clone()], BoolOp2::Common, fuzzy)
}
