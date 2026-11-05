//! `BRepBuilderAPI_Sewing` (`BRepBuilderAPI_Sewing.cxx`, 5,946 lines;
//! `.hxx:68-420`; `.lxx`).
//!
//! Line citations are on `D:\source\OCCT-src\src` (OCCT 8.0.0). Function
//! ranges are re-baselined against that tree; the statement-level ranges
//! inside a body still carry the older (7.9) offset, so locate a statement by
//! its function name first.
//!
//! Sewing pipeline (`Perform`, `cxx:2187-2343`):
//! `FaceAnalysis` (drop too-small edges/faces) -> `FindFreeBoundaries`
//! (edge -> faces map, node map) -> `VerticesAssembling` (glue near vertices)
//! -> `Cutting` (split bounds at other bounds' nodes) -> `Merging`
//! (same-parameter the contiguous bounds) -> `EdgeProcessing` /
//! `CreateSewedShape` / `EdgeRegularity` / `SameParameterShape` ->
//! `CreateOutputInformations`.
//!
//! The shape substitution context is `BRepTools_ReShape`: the port uses
//! [`MapReShape`], which now carries `Remove` (`BRepTools_ReShape.cxx:137-142`)
//! so a face deleted by `FaceAnalysis` disappears from its parents.
//!
//! **UNPORTED**: `Message_ProgressRange` / `Message_ProgressScope`
//! (`cxx:2191`, `:2610`, `:3343`, `:3553`, `:3794`, `:4468`, `:4917`, `:5009`)
//! — there is no progress port, so every `aPS.More()`
//! early return is unreachable and the full algorithm always runs. The
//! `OCCT_DEBUG` timers and `Dump`'s `std::cout` output (`cxx:2540-2584`) are
//! ported as a returned string.

mod candidates;
mod closed_surface;
mod cutting;
mod edge_processing;
mod face_analysis;
mod merging;
mod output;
mod same_parameter;
mod shape_maps;
mod vertices;

use std::sync::Arc;

use occt_core::precision::{CONFUSION, INFINITE, REAL_EPSILON};
use occt_geom::Surface;

use crate::abs::ShapeType;
use crate::shape::{Edge, Face, TopoShape};
use crate::shape_fix_compose_shell::MapReShape;

use shape_maps::*;

/// `BRepBuilderAPI_Sewing` (`BRepBuilderAPI_Sewing.hxx:76-427`).
pub struct Sewing {
    // Options (`hxx:87-98`).
    tolerance: f64,
    sewing: bool,
    analysis: bool,
    cutting: bool,
    nonmanifold: bool,
    face_mode: bool,
    floating_edges_mode: bool,
    same_parameter_mode: bool,
    local_tolerance_mode: bool,
    min_tolerance: f64,
    max_tolerance: f64,
    // Shape context (`hxx:395-412`).
    old_shapes: ShapeKeyedMap<TopoShape>,
    sewed_shape: Option<TopoShape>,
    degenerated: ShapeSet,
    free_edges: ShapeSet,
    multiple_edges: ShapeSet,
    contigous_edges: ShapeKeyedMap<Vec<TopoShape>>,
    contig_sec_bound: ShapeKeyedMap<TopoShape>,
    bound_faces: ShapeKeyedMap<Vec<TopoShape>>,
    bound_sections: ShapeKeyedMap<Vec<TopoShape>>,
    section_bound: ShapeKeyedMap<TopoShape>,
    vertex_node: ShapeKeyedMap<TopoShape>,
    vertex_node_free: ShapeKeyedMap<TopoShape>,
    node_sections: ShapeKeyedMap<Vec<TopoShape>>,
    cutting_node: ShapeKeyedMap<Vec<TopoShape>>,
    little_face: ShapeSet,
    shape: Option<TopoShape>,
    reshape: MapReShape,
    // Counters.
    nb_shapes: usize,
    nb_vertices: usize,
    #[allow(dead_code)]
    nb_edges: usize,
    // `myMergedEdges` (`hxx:424`).
    merged_edges: ShapeSet,
}

impl Default for Sewing {
    fn default() -> Self {
        Self::new(1.0e-6, true, true, true, false)
    }
}

impl Sewing {
    /// `BRepBuilderAPI_Sewing(tolerance, option1, option2, option3, option4)`
    /// (`cxx:2090-2103`).
    pub fn new(
        tolerance: f64,
        option_sewing: bool,
        option_analysis: bool,
        option_cutting: bool,
        option_nonmanifold: bool,
    ) -> Self {
        let mut s = Self {
            tolerance: 0.0,
            sewing: false,
            analysis: false,
            cutting: false,
            nonmanifold: false,
            face_mode: true,
            floating_edges_mode: false,
            same_parameter_mode: true,
            local_tolerance_mode: false,
            min_tolerance: 0.0,
            max_tolerance: INFINITE,
            old_shapes: ShapeKeyedMap::new(),
            sewed_shape: None,
            degenerated: ShapeSet::new(),
            free_edges: ShapeSet::new(),
            multiple_edges: ShapeSet::new(),
            contigous_edges: ShapeKeyedMap::new(),
            contig_sec_bound: ShapeKeyedMap::new(),
            bound_faces: ShapeKeyedMap::new(),
            bound_sections: ShapeKeyedMap::new(),
            section_bound: ShapeKeyedMap::new(),
            vertex_node: ShapeKeyedMap::new(),
            vertex_node_free: ShapeKeyedMap::new(),
            node_sections: ShapeKeyedMap::new(),
            cutting_node: ShapeKeyedMap::new(),
            little_face: ShapeSet::new(),
            shape: None,
            reshape: MapReShape::new(),
            nb_shapes: 0,
            nb_vertices: 0,
            nb_edges: 0,
            merged_edges: ShapeSet::new(),
        };
        s.init(
            tolerance,
            option_sewing,
            option_analysis,
            option_cutting,
            option_nonmanifold,
        );
        s
    }

    /// `Init(tolerance, option1, option2, option3, option4)`
    /// (`cxx:2104-2133`).
    pub fn init(
        &mut self,
        tolerance: f64,
        option_sewing: bool,
        option_analysis: bool,
        option_cutting: bool,
        option_nonmanifold: bool,
    ) {
        self.tolerance = tolerance.max(CONFUSION);
        self.sewing = option_sewing;
        self.analysis = option_analysis;
        self.cutting = option_cutting;
        self.nonmanifold = option_nonmanifold;
        self.min_tolerance = self.tolerance * 1e-4;
        if self.min_tolerance < CONFUSION {
            self.min_tolerance = CONFUSION;
        }
        self.max_tolerance = INFINITE;
        self.face_mode = true;
        self.floating_edges_mode = false;
        self.same_parameter_mode = true;
        self.local_tolerance_mode = false;
        self.sewed_shape = None;
        self.load(None);
    }

    /// `Load(shape)` (`cxx:2136-2167`).
    pub fn load(&mut self, shape: Option<&TopoShape>) {
        self.reshape = MapReShape::new();
        self.shape = match shape {
            None => None,
            Some(s) => self.reshape.apply_option(s),
        };
        self.sewed_shape = None;
        self.nb_shapes = 0;
        self.nb_edges = 0;
        self.nb_vertices = 0;
        self.old_shapes.clear();
        self.degenerated.clear();
        self.free_edges.clear();
        self.multiple_edges.clear();
        self.contigous_edges.clear();
        self.contig_sec_bound.clear();
        self.bound_faces.clear();
        self.bound_sections.clear();
        self.vertex_node.clear();
        self.vertex_node_free.clear();
        self.node_sections.clear();
        self.cutting_node.clear();
        self.section_bound.clear();
        self.little_face.clear();
    }

    /// `Add(shape)` (`cxx:2170-2180`).
    pub fn add(&mut self, shape: &TopoShape) {
        if shape.is_null() {
            return;
        }
        if let Some(o) = self.reshape.apply_option(shape) {
            self.old_shapes.bind(shape, o);
        }
        self.nb_shapes = self.old_shapes.extent();
    }

    /// `Perform` (`cxx:2187-2343`).
    ///
    /// Without a progress indicator (`Message_ProgressScope`) the stage order
    /// is the only control flow; every `aPS.More()` early exit is unreachable.
    pub fn perform(&mut self) -> Result<(), String> {
        if self.analysis {
            self.face_analysis();
        }
        if self.nb_shapes > 0 || self.shape.is_some() {
            self.find_free_boundaries();
            if self.bound_faces.extent() > 0 {
                self.vertices_assembling();
                if self.cutting {
                    self.cutting_bounds();
                }
                self.merging(true);
            }
            if self.sewing {
                self.edge_processing();
                self.create_sewed_shape();
                self.edge_regularity();
                if self.same_parameter_mode && self.face_mode {
                    self.same_parameter_shape();
                }
            }
            self.create_output_informations();
        }
        Ok(())
    }

    /// `SewedShape()` (`cxx:2352-2356`).
    pub fn sewed_shape(&self) -> Option<&TopoShape> {
        self.sewed_shape.as_ref()
    }

    /// `GetContext()` (`cxx:5935-5939`).
    pub fn context(&self) -> &MapReShape {
        &self.reshape
    }

    /// `SetContext(theContext)` (`cxx:5942-5946`).
    pub fn set_context(&mut self, context: MapReShape) {
        self.reshape = context;
    }

    /// `Tolerance()` / `SetTolerance(theToler)` (`.lxx:32-42`).
    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }
    pub fn set_tolerance(&mut self, t: f64) {
        self.tolerance = t;
    }

    /// `MinTolerance()` / `SetMinTolerance(theMinToler)` (`.lxx:46-57`).
    pub fn min_tolerance(&self) -> f64 {
        self.min_tolerance
    }
    pub fn set_min_tolerance(&mut self, t: f64) {
        self.min_tolerance = t;
    }

    /// `MaxTolerance()` / `SetMaxTolerance(theMaxToler)` (`.lxx:18-27`).
    pub fn max_tolerance(&self) -> f64 {
        self.max_tolerance
    }
    pub fn set_max_tolerance(&mut self, t: f64) {
        self.max_tolerance = t;
    }

    /// `FaceMode()` / `SetFaceMode(theFaceMode)` (`.lxx:60-71`).
    pub fn face_mode(&self) -> bool {
        self.face_mode
    }
    pub fn set_face_mode(&mut self, on: bool) {
        self.face_mode = on;
    }

    /// `FloatingEdgesMode()` / `SetFloatingEdgesMode(...)` (`.lxx:75-87`).
    pub fn floating_edges_mode(&self) -> bool {
        self.floating_edges_mode
    }
    pub fn set_floating_edges_mode(&mut self, on: bool) {
        self.floating_edges_mode = on;
    }

    /// `SameParameterMode()` / `SetSameParameterMode(...)` (`.lxx:105-116`).
    pub fn same_parameter_mode(&self) -> bool {
        self.same_parameter_mode
    }
    pub fn set_same_parameter_mode(&mut self, on: bool) {
        self.same_parameter_mode = on;
    }

    /// `LocalTolerancesMode()` / `SetLocalTolerancesMode(...)` (`.lxx:119-130`).
    pub fn local_tolerances_mode(&self) -> bool {
        self.local_tolerance_mode
    }
    pub fn set_local_tolerances_mode(&mut self, on: bool) {
        self.local_tolerance_mode = on;
    }

    /// `NonManifoldMode()` / `SetNonManifoldMode(...)` (`.lxx:133-143`).
    pub fn non_manifold_mode(&self) -> bool {
        self.nonmanifold
    }
    pub fn set_non_manifold_mode(&mut self, on: bool) {
        self.nonmanifold = on;
    }

    /// `NbFreeEdges()` / `FreeEdge(index)` (`cxx:2359-2363`).
    pub fn nb_free_edges(&self) -> usize {
        self.free_edges.extent()
    }
    pub fn free_edge(&self, index: usize) -> Option<Edge> {
        let s = self.free_edges.key_at(index)?;
        Edge::wrap(s.clone())
    }

    /// `NbMultipleEdges()` / `MultipleEdge(index)` (`cxx:2375-2388`).
    pub fn nb_multiple_edges(&self) -> usize {
        self.multiple_edges.extent()
    }
    pub fn multiple_edge(&self, index: usize) -> Option<Edge> {
        let s = self.multiple_edges.key_at(index)?;
        Edge::wrap(s.clone())
    }

    /// `NbContigousEdges()` / `ContigousEdge(index)` (`cxx:2391-2404`).
    pub fn nb_contigous_edges(&self) -> usize {
        self.contigous_edges.extent()
    }
    pub fn contigous_edge(&self, index: usize) -> Option<Edge> {
        let s = self.contigous_edges.key_at(index)?;
        Edge::wrap(s.clone())
    }

    /// `ContigousEdgeCouple(index)` (`cxx:2407-2414`).
    pub fn contigous_edge_couple(&self, index: usize) -> Option<&Vec<TopoShape>> {
        self.contigous_edges.value_at(index)
    }

    /// `IsSectionBound(section)` (`cxx:2417-2421`).
    pub fn is_section_bound(&self, section: &TopoShape) -> bool {
        self.contig_sec_bound.contains(section)
    }

    /// `SectionToBoundary(section)` (`cxx:2424-2430`).
    pub fn section_to_boundary(&self, section: &TopoShape) -> Option<Edge> {
        let b = self.contig_sec_bound.get(section)?;
        Edge::wrap(b)
    }

    /// `NbDegeneratedShapes()` / `DegeneratedShape(index)` (`cxx:2449-2462`).
    pub fn nb_degenerated_shapes(&self) -> usize {
        self.degenerated.extent()
    }
    pub fn degenerated_shape(&self, index: usize) -> Option<TopoShape> {
        self.degenerated.key_at(index).cloned()
    }

    /// `IsDegenerated(shape)` (`cxx:2465-2494`).
    pub fn is_degenerated(&mut self, shape: &TopoShape) -> bool {
        let new_shape = self.reshape.apply_option(shape);
        // Degenerated face: the replacement is null.
        if shape.shape_type() == ShapeType::Face {
            return new_shape.is_none();
        }
        let Some(new_shape) = new_shape else {
            return false;
        };
        if new_shape.shape_type() == ShapeType::Edge {
            return Edge::wrap(new_shape)
                .map(|e| crate::brep_tool::BRepTool::is_degenerated(&e))
                .unwrap_or(false);
        }
        if new_shape.shape_type() == ShapeType::Wire {
            let mut is_deg = true;
            for child in crate::iterator::cumulated_children(&new_shape) {
                if !is_deg {
                    break;
                }
                is_deg = match Edge::wrap(child) {
                    Some(e) => crate::brep_tool::BRepTool::is_degenerated(&e),
                    None => false,
                };
            }
            return is_deg;
        }
        false
    }

    /// `IsModified(shape)` (`cxx:2497-2510`).
    pub fn is_modified(&self, shape: &TopoShape) -> bool {
        let new_shape = self.old_shapes.get(shape).unwrap_or_else(|| shape.clone());
        !crate::topo_tools::is_same(&new_shape, shape)
    }

    /// `Modified(shape)` (`cxx:2513-2522`).
    pub fn modified<'a>(&'a self, shape: &'a TopoShape) -> &'a TopoShape {
        match self.old_shapes.find(shape) {
            Some(s) => s,
            None => shape,
        }
    }

    /// `IsModifiedSubShape(shape)` (`cxx:2525-2530`).
    pub fn is_modified_sub_shape(&mut self, shape: &TopoShape) -> bool {
        match self.reshape.apply_option(shape) {
            Some(s) => !crate::topo_tools::is_same(&s, shape),
            None => true,
        }
    }

    /// `ModifiedSubShape(shape)` (`cxx:2533-2537`).
    pub fn modified_sub_shape(&mut self, shape: &TopoShape) -> Option<TopoShape> {
        self.reshape.apply_option(shape)
    }

    /// `NbDeletedFaces()` / `DeletedFace(index)` (`cxx:2433-2446`).
    pub fn nb_deleted_faces(&self) -> usize {
        self.little_face.extent()
    }
    pub fn deleted_face(&self, index: usize) -> Option<Face> {
        let s = self.little_face.key_at(index)?;
        Face::wrap(s.clone())
    }

    /// `WhichFace(edg, index)` (`cxx:145-166`): the `index`-th face the bound
    /// of `theEdg` is attached to.
    pub fn which_face(&self, edge: &TopoShape, index: usize) -> Option<Face> {
        let mut bound = edge.clone();
        if let Some(b) = self.section_bound.get(edge) {
            bound = b;
        }
        let faces = self.bound_faces.find(&bound)?;
        let f = faces.get(index.checked_sub(1)?)?;
        Face::wrap(f.clone())
    }

    /// `Dump()` (`cxx:2540-2584`) as a string.
    pub fn dump(&mut self) -> String {
        let nb_bounds = self.bound_faces.extent();
        let mut nb_sections = 0usize;
        let mut edges = ShapeSet::new();
        let mut vertices = ShapeSet::new();
        for i in 1..=nb_bounds {
            let bound = self.bound_faces.key_at(i).cloned().unwrap();
            match self.bound_sections.find(&bound) {
                Some(l) => nb_sections += l.len(),
                None => nb_sections += 1,
            }
            let applied = self.reshape.apply_option(&bound);
            if let Some(ap) = applied {
                for edge in crate::topo_tools::explorer_shapes(&ap, ShapeType::Edge) {
                    edges.add(&edge);
                    let e = Edge::wrap(edge);
                    if let Some(e) = e {
                        let (v1, v2) = crate::topo_tools::edge_vertices(&e);
                        if let Some(v) = v1 {
                            vertices.add(&v.0);
                        }
                        if let Some(v) = v2 {
                            vertices.add(&v.0);
                        }
                    }
                }
            }
        }
        let mut out = String::new();
        out.push_str("\n                        Information                         \n");
        out.push_str(" ===========================================================\n");
        out.push_str(&format!(
            " Number of input shapes      : {}\n",
            self.old_shapes.extent()
        ));
        out.push_str(&format!(" Number of actual shapes     : {}\n", self.nb_shapes));
        out.push_str(&format!(" Number of Bounds            : {}\n", nb_bounds));
        out.push_str(&format!(" Number of Sections          : {}\n", nb_sections));
        out.push_str(&format!(" Number of Edges             : {}\n", edges.extent()));
        out.push_str(&format!(" Number of Vertices          : {}\n", self.nb_vertices));
        out.push_str(&format!(" Number of Nodes             : {}\n", vertices.extent()));
        out.push_str(&format!(
            " Number of Free Edges        : {}\n",
            self.free_edges.extent()
        ));
        out.push_str(&format!(
            " Number of Contiguous Edges  : {}\n",
            self.contigous_edges.extent()
        ));
        out.push_str(&format!(
            " Number of Multiple Edges    : {}\n",
            self.multiple_edges.extent()
        ));
        out.push_str(&format!(
            " Number of Degenerated Edges : {}\n",
            self.degenerated.extent()
        ));
        out.push_str(" ===========================================================\n");
        out
    }

    // --- helpers shared by the stage modules -----------------------------

    /// `BRep_Tool::Surface(face, loc)` with the location composed in
    /// (`BRep_Tool.cxx:430-460`): the surface the sewing stages sample.
    pub(crate) fn face_surface_located(&self, face: &Face) -> Option<Arc<dyn Surface>> {
        let surf = crate::brep_tool::BRepTool::face_surface(face)?;
        let loc = face.0.location();
        if loc.is_identity() {
            Some(surf)
        } else {
            Some(Arc::from(surf.transformed(&loc.transformation())))
        }
    }

    /// `TopoDS_Iterator(shape, /*cumOri=*/false)` children: the stored
    /// components with the orientation and location they were added with.
    pub(crate) fn raw_children(shape: &TopoShape) -> Vec<TopoShape> {
        shape.tshape.read().expect("poisoned TShape lock").children.clone()
    }

    /// `BRep_Builder`'s `RealEpsilon()` use in `ComputeToleranceVertex`
    /// (`cxx:522`).
    pub(crate) const EPS: f64 = REAL_EPSILON;
}



