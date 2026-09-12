use super::prelude::*;
use super::*;

/// Result-shape assembly order — matches the per-type `BuildResult` calls of
/// OCCT `BOPAlgo_Builder::PerformInternal1` (lowest type first).

pub(super) const RESULT_TYPES: [ShapeType; 8] = [
    ShapeType::Vertex,
    ShapeType::Edge,
    ShapeType::Wire,
    ShapeType::Face,
    ShapeType::Shell,
    ShapeType::Solid,
    ShapeType::CompSolid,
    ShapeType::Compound,
];

/// The General Fuse algorithm — base algorithm of the Boolean Component.
///
/// Source: `BOPAlgo_Builder`.
///
/// The class owns the [`PaveFiller`] (intersection phase), the
/// [`BopHistory`] images table (each old sub-shape → its split pieces), the
/// origins back-map, the same-domain shapes map and the alert report of the
/// run. FillImages / BuildResult are General Fuse: they do not read object
/// vs tool IN/OUT. Those states are stored for `BOPAlgo_BOP::BuildShape`.
#[derive(Debug, Clone)]
pub struct BopBuilder {
    /// Pave Filler — the intersection phase of the algorithm.
    pub(super) filler: PaveFiller,
    /// Images naming table + history: each old shape → the pieces it was
    /// split/expanded into during the intersection (`myImages`).
    pub(super) history: BopHistory,
    /// Origins — the back map of the images table: split shape → the source
    /// shapes that expand into it, keyed by the split shape's `TShape` address
    /// (`myOrigins`).
    pub(super) origins: HashMap<usize, Vec<TopoShape>>,
    /// Same-domain shapes: source index → the coincident shape it maps to
    /// (`myShapesSD`).
    pub(super) shapes_sd: HashMap<usize, usize>,
    /// Same-domain bindings keyed by TShape identity, covering reconstructed
    /// faces that are not in the DS (`myShapesSD` on shape handles).
    pub(super) sd_by_key: HashMap<usize, TopoShape>,
    /// Arguments of the operation (deduplicated by TShape identity).
    pub(super) arguments: Vec<TopoShape>,
    /// Objects group of the current BOP run.
    pub(super) objects: Vec<TopoShape>,
    /// Tools group of the current BOP run.
    pub(super) tools: Vec<TopoShape>,
    /// Object IN/OUT for `BOPAlgo_BOP::BuildShape` / `BuildBOP` only.
    /// FillImagesSolids and BuildResult do not read this.
    pub(super) obj_state: FaceState,
    /// Tool IN/OUT for `BOPAlgo_BOP::BuildShape` / `BuildBOP` only.
    pub(super) tools_state: FaceState,
    /// Fatal alerts — a non-empty list means the algorithm has failed.
    pub(super) errors: Vec<String>,
    /// Non-fatal alerts.
    pub(super) warnings: Vec<String>,
    /// The shape produced by the last run.
    pub(super) result_shape: TopoShape,
    /// `myInParts`: source solid TShape key → IN + INTERNAL faces from
    /// FillIn3DParts. Consumed by `BOPAlgo_Builder::BuildBOP`.
    pub(super) in_parts: HashMap<usize, Vec<TopoShape>>,
}

impl Default for BopBuilder {
    fn default() -> Self {
        Self::new()
    }
}
