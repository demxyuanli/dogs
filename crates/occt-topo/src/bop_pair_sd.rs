//! Parallel-algo payloads used by FillSameDomainFaces and FillInternalVertices.
//!
//! Source: `BOPAlgo_Builder_2.cxx`:
//! * `BOPAlgo_PairOfShapeBoolean` at `:63-113` (`Perform` at `:94`);
//! * `BOPAlgo_VFI` at `:151-208` (`Perform` at `:188`).
//!
//! OCCT runs these through `BOPTools_Parallel::Perform`. The port executes
//! them sequentially with a shared [`IntToolsContext`], which is the
//! `myRunParallel == false` path of that helper.

use crate::algo_tools_face::are_faces_same_domain;
use crate::int_tools_full::IntToolsContext;
use crate::shape::{Face, TopoShape, Vertex};

/// `BOPAlgo_PairOfShapeBoolean` — two faces plus the same-domain flag.
///
/// `Perform` (`_2.cxx:94-105`) casts both shapes to `TopoDS_Face` and stores
/// `BOPTools_AlgoTools::AreFacesSameDomain(aFj, aFk, myContext, myFuzzyValue)`
/// in `myFlag`.
#[derive(Debug, Clone)]
pub struct PairOfShapeBoolean {
    shape1: TopoShape,
    shape2: TopoShape,
    flag: bool,
    fuzzy: f64,
}

impl PairOfShapeBoolean {
    /// Empty constructor (`myFlag = false`).
    pub fn new() -> Self {
        Self {
            shape1: TopoShape::new(crate::abs::ShapeType::Face),
            shape2: TopoShape::new(crate::abs::ShapeType::Face),
            flag: false,
            fuzzy: 0.0,
        }
    }

    /// `Shape1()`.
    pub fn shape1(&self) -> &TopoShape {
        &self.shape1
    }

    /// `Shape2()`.
    pub fn shape2(&self) -> &TopoShape {
        &self.shape2
    }

    /// `Shape1()` mutable (`aPSB.Shape1() = aF1`).
    pub fn shape1_mut(&mut self) -> &mut TopoShape {
        &mut self.shape1
    }

    /// `Shape2()` mutable.
    pub fn shape2_mut(&mut self) -> &mut TopoShape {
        &mut self.shape2
    }

    /// `Flag()` — true when the pair is same-domain.
    pub fn flag(&self) -> bool {
        self.flag
    }

    /// `SetFuzzyValue` (inherited from `BOPAlgo_Options`, `_2.cxx:705`).
    pub fn set_fuzzy_value(&mut self, fuzzy: f64) {
        self.fuzzy = fuzzy;
    }

    /// `BOPAlgo_PairOfShapeBoolean::Perform` (`_2.cxx:94-105`).
    ///
    /// Progress-scope / UserBreak is omitted; the geometric test is the
    /// whole body.
    pub fn perform(&mut self, ctx: &mut IntToolsContext) {
        let fj = Face(self.shape1.clone());
        let fk = Face(self.shape2.clone());
        self.flag = are_faces_same_domain(&fj, &fk, ctx, self.fuzzy);
    }
}

impl Default for PairOfShapeBoolean {
    fn default() -> Self {
        Self::new()
    }
}

/// `BOPAlgo_VectorOfPairOfShapeBoolean` — run every pair.
pub fn perform_pair_vector(pairs: &mut [PairOfShapeBoolean], ctx: &mut IntToolsContext) {
    for p in pairs.iter_mut() {
        p.perform(ctx);
    }
}

/// `BOPAlgo_VFI` — classify one vertex as internal to one face.
///
/// `Perform` (`_2.cxx:188-200`) calls `IntTools_Context::ComputeVF` and
/// stores `myIsInternal = (iFlag == 0)`.
#[derive(Debug, Clone)]
pub struct VertexFaceInternal {
    vertex: Vertex,
    face: Face,
    is_internal: bool,
    fuzzy: f64,
}

impl VertexFaceInternal {
    /// Empty constructor (`myIsInternal = false`).
    pub fn new() -> Self {
        Self {
            vertex: Vertex(TopoShape::new(crate::abs::ShapeType::Vertex)),
            face: Face(TopoShape::new(crate::abs::ShapeType::Face)),
            is_internal: false,
            fuzzy: 0.0,
        }
    }

    /// `SetVertex`.
    pub fn set_vertex(&mut self, v: Vertex) {
        self.vertex = v;
    }

    /// `Vertex()` — OCCT returns a non-const `TopoDS_Vertex&`
    /// (`BOPAlgo_Builder_2.cxx:168`) and `FillInternalVertices` binds it as
    /// `TopoDS_Vertex& aV = aVFI.Vertex()` (`_2.cxx:904`).
    pub fn vertex_mut(&mut self) -> &mut Vertex {
        &mut self.vertex
    }

    /// `SetFace`.
    pub fn set_face(&mut self, f: Face) {
        self.face = f;
    }

    /// `Face()` — non-const `TopoDS_Face&` (`_2.cxx:174`), bound as
    /// `TopoDS_Face& aF = aVFI.Face()` at `_2.cxx:905`.
    pub fn face_mut(&mut self) -> &mut Face {
        &mut self.face
    }

    /// `IsInternal()`.
    pub fn is_internal(&self) -> bool {
        self.is_internal
    }

    /// `SetFuzzyValue`.
    pub fn set_fuzzy_value(&mut self, fuzzy: f64) {
        self.fuzzy = fuzzy;
    }

    /// `BOPAlgo_VFI::Perform` (`_2.cxx:188-200`).
    ///
    /// `ComputeVF` returns 0 when the vertex projects strictly inside the
    /// face within the summed tolerances.
    pub fn perform(&mut self, ctx: &mut IntToolsContext) {
        let i_flag = ctx.compute_vf(&self.vertex, &self.face, self.fuzzy);
        self.is_internal = i_flag == 0;
    }
}

impl Default for VertexFaceInternal {
    fn default() -> Self {
        Self::new()
    }
}

/// `BOPAlgo_VectorOfVFI` — run every vertex/face pair.
pub fn perform_vfi_vector(pairs: &mut [VertexFaceInternal], ctx: &mut IntToolsContext) {
    for p in pairs.iter_mut() {
        p.perform(ctx);
    }
}
