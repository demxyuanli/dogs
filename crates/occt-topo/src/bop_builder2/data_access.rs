use super::prelude::*;
use super::*;

/// The face-reconstruction host contract (`crate::bop_build_faces`).
///
/// Every accessor delegates to the corresponding inherent [`BopBuilder`]
/// method; the trait is what lets the sibling Phase-20 modules
/// (`crate::bop_build_faces::fill_images_faces` & friends) run against this
/// builder.
impl BopBuilderLike for BopBuilder {
    fn ds(&self) -> &BopdsDS {
        self.ds()
    }
    fn ds_mut(&mut self) -> &mut BopdsDS {
        self.ds_mut()
    }
    fn history(&self) -> &BopHistory {
        self.history()
    }
    fn history_mut(&mut self) -> &mut BopHistory {
        self.history_mut()
    }
    fn has_errors(&self) -> bool {
        self.has_errors()
    }
    fn add_error(&mut self, msg: String) {
        self.add_error(msg);
    }
    fn add_warning(&mut self, msg: String) {
        self.add_warning(msg);
    }
    fn non_destructive(&self) -> bool {
        self.non_destructive()
    }
    fn fuzzy_value(&self) -> f64 {
        self.fuzzy_value()
    }
    fn bind_shapes_sd(&mut self, shape: TopoShape, sd: TopoShape) {
        self.bind_shapes_sd(shape, sd);
    }
    fn seek_shapes_sd(&self, shape: &TopoShape) -> Option<TopoShape> {
        self.seek_shapes_sd(shape)
    }
    fn origins_mut(&mut self) -> &mut HashMap<usize, Vec<TopoShape>> {
        self.origins_mut()
    }
}

/// The container/internal-shape host contract (`crate::bop_build_common`).
///
/// All accessors delegate to inherent [`BopBuilder`] methods; the trait lets
/// the solid stage (`crate::bop_build_solids`) and the container stages
/// (`crate::bop_build_common`) run against this builder.
impl BopBuildOps for BopBuilder {
    fn ds(&self) -> &BopdsDS {
        self.ds()
    }
    fn history(&self) -> &BopHistory {
        self.history()
    }
    fn history_mut(&mut self) -> &mut BopHistory {
        self.history_mut()
    }
    fn fuzzy_value(&self) -> f64 {
        self.fuzzy_value()
    }
    fn arguments(&self) -> &[TopoShape] {
        self.arguments()
    }
    fn origins_mut(&mut self) -> &mut HashMap<usize, Vec<TopoShape>> {
        self.origins_mut()
    }
}

impl BopSolidHost for BopBuilder {
    fn ds(&self) -> &BopdsDS {
        self.ds()
    }
    fn history(&self) -> &BopHistory {
        self.history()
    }
    fn history_mut(&mut self) -> &mut BopHistory {
        self.history_mut()
    }
    fn fuzzy_value(&self) -> f64 {
        self.fuzzy_value()
    }
    fn arguments(&self) -> &[TopoShape] {
        self.arguments()
    }
    fn origins_mut(&mut self) -> &mut HashMap<usize, Vec<TopoShape>> {
        self.origins_mut()
    }
    fn add_warning(&mut self, msg: String) {
        self.add_warning(msg);
    }
    fn seek_shapes_sd(&self, shape: &TopoShape) -> Option<TopoShape> {
        self.seek_shapes_sd(shape)
    }
    fn bind_shapes_sd(&mut self, shape: TopoShape, sd: TopoShape) {
        self.bind_shapes_sd(shape, sd);
    }
}

/// The Boolean operation to apply. Source: `BOPAlgo_Operation`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoolOp2 {
    /// A ∪ B.
    Fuse,
    /// A − B.
    Cut,
    /// A ∩ B.
    Common,
}

impl BoolOp2 {
    /// Converts the operation into the states of the object/tool faces that
    /// pass into the result — the `BOPAlgo_Builder::BuildBOP` conversion:
    ///
    /// - [`BoolOp2::Fuse`] — faces of the objects OUT of the tools + faces of
    ///   the tools OUT of the objects → the union. States `(Out, Out)`.
    /// - [`BoolOp2::Cut`] — faces of the objects OUT of the tools + faces of
    ///   the tools IN the objects → object minus tools. States `(Out, In)`.
    /// - [`BoolOp2::Common`] — faces of the objects IN the tools + faces of
    ///   the tools IN the objects → the intersection. States `(In, In)`.
    pub fn states(&self) -> (FaceState, FaceState) {
        match self {
            BoolOp2::Fuse => (FaceState::Out, FaceState::Out),
            BoolOp2::Cut => (FaceState::Out, FaceState::In),
            BoolOp2::Common => (FaceState::In, FaceState::In),
        }
    }
}

/// Performs the Boolean operation of type `op` on the `objects` and `tools`
/// groups.
///
/// Port of `BOPAlgo_BOP::Perform`: both groups are intersected together (GF
/// FillImages + BuildResult), then the operation type is converted into
/// object/tool states for `BOPAlgo_BOP::BuildShape`.
pub fn builder_bop(
    objects: &[TopoShape],
    tools: &[TopoShape],
    op: BoolOp2,
) -> Result<TopoShape, String> {
    builder_bop_with_fuzzy(objects, tools, op, 1e-7)
}

/// [`builder_bop`] with an explicit fuzzy value (`BOPAlgo_PaveFiller::SetFuzzyValue`).
pub fn builder_bop_with_fuzzy(
    objects: &[TopoShape],
    tools: &[TopoShape],
    op: BoolOp2,
    fuzzy: f64,
) -> Result<TopoShape, String> {
    let objects: Vec<TopoShape> = objects.iter().map(flatten_location).collect();
    let tools: Vec<TopoShape> = tools.iter().map(flatten_location).collect();
    let mut all: Vec<TopoShape> = objects.clone();
    for t in &tools {
        if !all.iter().any(|a| a.same_tshape(t)) {
            all.push(t.clone());
        }
    }
    let mut b = BopBuilder::new();
    b.set_fuzzy_value(fuzzy);
    b.set_arguments(&all);
    let (os, ts) = op.states();
    b.perform_internal(&objects, os, &tools, ts)
}

/// Union of two groups of shapes. Source: `BOPAlgo_BOP` with `BOPAlgo_FUSE`.
pub fn fuse(objects: &[TopoShape], tools: &[TopoShape]) -> Result<TopoShape, String> {
    builder_bop(objects, tools, BoolOp2::Fuse)
}

/// Subtraction of the tools from the objects. Source: `BOPAlgo_CUT`.
pub fn cut(objects: &[TopoShape], tools: &[TopoShape]) -> Result<TopoShape, String> {
    builder_bop(objects, tools, BoolOp2::Cut)
}

/// Intersection of the two groups. Source: `BOPAlgo_COMMON`.
pub fn common(objects: &[TopoShape], tools: &[TopoShape]) -> Result<TopoShape, String> {
    builder_bop(objects, tools, BoolOp2::Common)
}

/// Stable identity key of a shape — the address of its shared `TShape`
/// (the same identity `TopTools_ShapeMapHasher` uses).
pub(super) fn shape_key(s: &TopoShape) -> usize {
    Arc::as_ptr(&s.tshape) as usize
}

/// Bakes a shape's `TopoShape` location into its geometry, returning a copy at
/// the world position with an identity location. Shapes without a location are
/// returned unchanged (the shared `TShape` is preserved).
///
/// The Phase 19/20 intersection and rebuild stages read geometry from the
/// registry without applying the `TopoShape` location, so a located argument
/// (e.g. `crate::transform::translated`) would otherwise be intersected in its
/// local frame. Flattening the location at the builder boundary keeps the rest
/// of the pipeline location-agnostic.
pub(super) fn flatten_location(s: &TopoShape) -> TopoShape {
    if s.location.is_identity() {
        return s.clone();
    }
    let t = s.location.transformation();
    crate::shape_ops::transformed_copy(s, &t).unwrap_or_else(|_| s.clone())
}

/// True if the pave block `pb` lies on a common block of the data structure:
/// its edge (or original edge) is an edge of a common block and its parameter
/// range matches one of the common-block ranges.
pub(super) fn is_common_block_on_edge(ds: &BopdsDS, pb: &BopdsPaveBlock) -> bool {
    ds.common_blocks().iter().any(|cb| {
        let on_edge = cb.contains_index(pb.edge()) || cb.contains_index(pb.original_edge());
        on_edge && cb.contains_range(pb.first, pb.last, 1e-7)
    })
}
