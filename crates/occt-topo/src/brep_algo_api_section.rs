//! `BRepAlgoAPI_Section` — the SECTION high-level API.
//!
//! Source: `BRepAlgoAPI_Section.cxx:48-333` (TKBO/BRepAlgoAPI). The class is a
//! thin front end: it stores `Objects` / `Tools`, pins the operation to
//! `BOPAlgo_SECTION` with `myApprox` / `myComputePCurve1` / `myComputePCurve2`
//! off, delegates the run to `BOPAlgo_Section` ([`BopBuilder::perform_section`])
//! and keeps the filler for the `HasAncestorFaceOn1/2` queries.
//!
//! UNPORTED: `BRepAlgoAPI_Section.cxx:85-121` / `:150-186` / `:324-333` — the
//! `gp_Pln` / `Geom_Surface` overloads and the file-local `MakeShape` need
//! `BRepBuilderAPI_MakeFace(Geom_Surface, tol)` / `BRepBuilderAPI_MakeShell`
//! on a bare surface, which this crate does not have; every caller in the
//! current tree passes faces.
//!
//! UNPORTED: `BOPAlgo_PaveFiller_6.cxx:506` — `IntTools_FaceFace::SetParameters`
//! is never called by the port's F/F stage for *any* operation, so the
//! `BOPAlgo_SectionAttribute(false, false, false)` that `SetAttributes`
//! (`BRepAlgoAPI_Section.cxx:209-213`) installs has no effect here. OCCT builds
//! neither the 3-D approximation nor the pcurves from that flag; the port always
//! keeps the intersector's pcurves, which is the *builder* default.
use crate::abs::ShapeType;
use crate::bop_builder2::BopBuilder;
use crate::shape::TopoShape;
use crate::topo_tools::is_same;

/// Result of the SECTION operation plus the filler needed for the ancestor
/// queries. Source: `BRepAlgoAPI_Section` + `BRepAlgoAPI_BuilderAlgo`.
pub struct SectionAlgo {
    /// `myArguments` — the first operand group (`Init1`).
    arguments: Vec<TopoShape>,
    /// `myTools` — the second operand group (`Init2`).
    tools: Vec<TopoShape>,
    /// The builder that ran (or will run) the section.
    builder: BopBuilder,
    /// `myShape` — the section result.
    shape: TopoShape,
    /// `myIsDone`.
    done: bool,
}

impl Default for SectionAlgo {
    fn default() -> Self {
        Self::new()
    }
}

impl SectionAlgo {
    /// `BRepAlgoAPI_Section()` (`cxx:48-51`): `Init(false)`.
    pub fn new() -> Self {
        Self {
            arguments: Vec::new(),
            tools: Vec::new(),
            builder: BopBuilder::new(),
            shape: TopoShape::new(ShapeType::Compound),
            done: false,
        }
    }

    /// `BRepAlgoAPI_Section(Sh1, Sh2, PerformNow)` (`cxx:63-71`).
    pub fn with_shapes(s1: &TopoShape, s2: &TopoShape, perform_now: bool) -> Self {
        let mut s = Self::new();
        s.init1(s1);
        s.init2(s2);
        s.init(perform_now);
        s
    }

    /// `BRepAlgoAPI_Section::Init` (`cxx:128-140`): pin the operation to
    /// `BOPAlgo_SECTION` with the three section-attribute flags off.
    fn init(&mut self, perform_now: bool) {
        if perform_now {
            let _ = self.build();
        }
    }

    /// `Init1` (`cxx:144-166`): the first operand group replaces `myArguments`.
    pub fn init1(&mut self, s: &TopoShape) {
        self.arguments.clear();
        self.arguments.push(s.clone());
    }

    /// `Init2` (`cxx:168-184`): the second operand group replaces `myTools`.
    pub fn init2(&mut self, s: &TopoShape) {
        self.tools.clear();
        self.tools.push(s.clone());
    }

    /// True when the last run produced a result (`myIsDone`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `myShape` — a compound of the section vertices and edges.
    pub fn shape(&self) -> &TopoShape {
        &self.shape
    }

    /// `BRepAlgoAPI_Section::Build` (`cxx:217-220`) on top of
    /// `BRepAlgoAPI_BooleanOperation::Build` (`cxx:125-230`).
    ///
    /// The intersection runs over `Objects + Tools` together
    /// (`BRepAlgoAPI_BooleanOperation.cxx:176-183`), then `BOPAlgo_Section`
    /// gets that combined list as its arguments
    /// (`BRepAlgoAPI_BooleanOperation.cxx:197-200`).
    pub fn build(&mut self) -> Result<&TopoShape, String> {
        self.done = false;
        if self.arguments.is_empty() || self.tools.is_empty() {
            return Err("BRepAlgoAPI_Section: too few arguments".to_string());
        }
        let mut all: Vec<TopoShape> = Vec::new();
        for s in self.arguments.iter().chain(self.tools.iter()) {
            if !all.iter().any(|a| is_same(a, s)) {
                all.push(s.clone());
            }
        }
        self.builder.set_arguments(&all);
        let objects = self.arguments.clone();
        let tools = self.tools.clone();
        let result = self.builder.perform_section(&objects, &tools)?;
        self.shape = result;
        self.done = true;
        Ok(&self.shape)
    }

    /// `BRepAlgoAPI_Section::HasAncestorFaceOn1` (`cxx:225-231`).
    pub fn has_ancestor_face_on1(&self, e: &TopoShape) -> Option<TopoShape> {
        self.has_ancestor_face(1, e)
    }

    /// `BRepAlgoAPI_Section::HasAncestorFaceOn2` (`cxx:236-242`).
    pub fn has_ancestor_face_on2(&self, e: &TopoShape) -> Option<TopoShape> {
        self.has_ancestor_face(2, e)
    }

    /// `HasAncestorFace` (`cxx:246-272`).
    fn has_ancestor_face(&self, index: i32, e: &TopoShape) -> Option<TopoShape> {
        if e.shape_type() != ShapeType::Edge {
            return None;
        }
        let (f1, f2) = self.has_ancestor_faces(e)?;
        Some(if index == 1 { f1 } else { f2 })
    }

    /// `HasAncestorFaces` (`cxx:276-324`): the section edge `aEx` is reached
    /// through the `PaveBlocks` of the `InterfFF` curves; the two faces of that
    /// record are the answer.
    fn has_ancestor_faces(&self, e: &TopoShape) -> Option<(TopoShape, TopoShape)> {
        let ds = self.builder.ds();
        for ff in ds.interf_ff() {
            let (n_f1, n_f2) = ff.indices();
            for c in ff.curves() {
                for pb in c.pave_blocks() {
                    let n_e = pb.edge();
                    if n_e == 0 {
                        continue;
                    }
                    let Some(shape) = ds.shape(n_e) else { continue };
                    if is_same(shape, e) {
                        let f1 = ds.shape(n_f1)?.clone();
                        let f2 = ds.shape(n_f2)?.clone();
                        return Some((f1, f2));
                    }
                }
            }
        }
        None
    }

    /// Shared edges of the section result, in result order. Convenience for the
    /// `TopExp_Explorer(shape, TopAbs_EDGE)` walk the callers of
    /// `BRepAlgoAPI_Section` do.
    pub fn result_edges(&self) -> Vec<TopoShape> {
        crate::topo_tools::explorer_shapes(&self.shape, ShapeType::Edge)
    }
}

/// Runs the SECTION operation on two shapes and answers the result shape.
///
/// `BRepAlgoAPI_Section(Sh1, Sh2, true)` (`cxx:63-71`) followed by `Shape()`.
pub fn section(s1: &TopoShape, s2: &TopoShape) -> Result<TopoShape, String> {
    section_with_fuzzy(s1, s2, 1e-7)
}

/// [`section`] with an explicit fuzzy value (`BOPAlgo_PaveFiller::SetFuzzyValue`).
pub fn section_with_fuzzy(
    s1: &TopoShape,
    s2: &TopoShape,
    fuzzy: f64,
) -> Result<TopoShape, String> {
    let mut algo = SectionAlgo::with_shapes(s1, s2, false);
    algo.builder.set_fuzzy_value(fuzzy);
    algo.build()?;
    Ok(algo.shape().clone())
}
