use super::super::prelude::*;
use super::super::*;

impl BopBuilder {
    // -----------------------------------------------------------------------
    // Arguments
    // -----------------------------------------------------------------------

    /// Sets the list of arguments for the operation, deduplicated by `TShape`
    /// identity (mirrors `BOPAlgo_Builder::SetArguments` + `AddArgument` fence).
    pub fn set_arguments(&mut self, shapes: &[TopoShape]) {
        self.arguments.clear();
        for s in shapes {
            let flat = flatten_location(s);
            if !self.arguments.iter().any(|a| is_same(a, &flat)) {
                self.arguments.push(flat);
            }
        }
    }

    /// Adds an argument to the operation, skipping duplicates.
    pub fn add_argument(&mut self, s: &TopoShape) {
        let flat = flatten_location(s);
        if !self.arguments.iter().any(|a| is_same(a, &flat)) {
            self.arguments.push(flat);
        }
    }

    /// `BOPAlgo_ToolsProvider::AddTool` (`BOPAlgo_ToolsProvider.cxx:47-53`):
    /// appends to the tool group when the shape is not already there.
    pub fn add_tool(&mut self, s: &TopoShape) {
        let flat = flatten_location(s);
        if !self.tools.iter().any(|t| is_same(t, &flat)) {
            self.tools.push(flat);
        }
    }

    /// `BOPAlgo_ToolsProvider::SetTools` (`BOPAlgo_ToolsProvider.cxx:57-65`).
    pub fn set_tools(&mut self, shapes: &[TopoShape]) {
        self.tools.clear();
        for s in shapes {
            self.add_tool(s);
        }
    }

    /// Returns the arguments of the operation.
    pub fn arguments(&self) -> &[TopoShape] {
        &self.arguments
    }

    /// Returns the objects group of the current run.
    pub fn objects(&self) -> &[TopoShape] {
        &self.objects
    }

    /// Returns the tools group of the current run.
    pub fn tools(&self) -> &[TopoShape] {
        &self.tools
    }

    /// Returns the state the object faces must have to pass into the result.
    pub fn obj_state(&self) -> FaceState {
        self.obj_state
    }

    /// Returns the state the tool faces must have to pass into the result.
    pub fn tools_state(&self) -> FaceState {
        self.tools_state
    }

    /// Stores the operation's IN/OUT states without running anything.
    ///
    /// `BOPAlgo_BOP::SetOperation` keeps `myOperation`; the port carries it as
    /// the state pair `BOPAlgo_BOP::BuildShape` / `BuildRC` read, so a subclass
    /// that re-runs the building stages (`BRepFeat_Builder::PerformResult`,
    /// `BRepFeat_Builder.cxx:194`) must be able to set it.
    pub fn set_bop_states(&mut self, obj_state: FaceState, tools_state: FaceState) {
        self.obj_state = obj_state;
        self.tools_state = tools_state;
    }

    /// `CheckArgsForOpenSolid` (`BOPAlgo_BOP::BuildShape`, `cxx:886`): false
    /// when the subclass suppresses the alternative `BuildBOP` path.
    pub fn check_args_for_open_solid(&self) -> bool {
        self.check_open_solid
    }

    /// Sets the `CheckArgsForOpenSolid` override (`BRepFeat_Builder.hxx:119-121`).
    pub fn set_check_args_for_open_solid(&mut self, flag: bool) {
        self.check_open_solid = flag;
    }
}
