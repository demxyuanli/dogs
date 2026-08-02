//! Meshing context — port of `IMeshTools_Context` plus the algorithm
//! interfaces it drives: `IMeshTools_MeshAlgo`, `IMeshTools_ModelAlgo`,
//! `IMeshTools_ModelBuilder` and `IMeshTools_ShapeExplorer`.
//!
//! A [`MeshContext`] caches the discrete [`MeshModel`] built for a shape and
//! the tools that process it. The pipeline it orchestrates mirrors
//! `BRepMesh_IncrementalMesh`:
//!
//! 1. [`MeshBuilder`] builds the discrete model (shape → model).
//! 2. [`ModelAlgo`]s discretize edges, heal, pre/post-process and mesh faces.
//! 3. [`MeshAlgo`] is the per-face triangulation interface used inside the
//!    face discretizer.
//! 4. [`ShapeExplorer`] walks a shape's faces and free edges and feeds them to
//!    a [`ShapeVisitor`](super::shape_tool::ShapeVisitor).

use std::collections::HashSet;
use std::sync::Arc;

use crate::abs::Orientation;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools_full::{edges_of, faces_of};

use super::data_model::{MeshFace, MeshModel};
use super::parameters::MeshParameters;
use super::shape_tool::ShapeVisitor;

/// Identity key of a `TShape` (its heap address). Used to deduplicate shared
/// sub-shapes that are visited more than once.
fn tshape_ptr(s: &TopoShape) -> usize {
    Arc::as_ptr(&s.tshape) as usize
}

/// Algorithm that builds a triangle mesh for a single discrete face.
/// Source: `IMeshTools_MeshAlgo::Perform`.
pub trait MeshAlgo {
    /// Meshes the given discrete face.
    fn perform(&mut self, face: &mut MeshFace, parameters: &MeshParameters) -> Result<(), String>;
}

/// Algorithm that updates or modifies a whole discrete model.
/// Source: `IMeshTools_ModelAlgo::Perform`.
pub trait ModelAlgo {
    /// Processes the given discrete model.
    fn perform(&mut self, model: &mut MeshModel, parameters: &MeshParameters) -> Result<(), String>;
}

/// Tool that builds a discrete model from a topological shape.
/// Source: `IMeshTools_ModelBuilder::Perform`.
pub trait MeshBuilder {
    /// Builds the discrete model for the given shape.
    fn build_model(&self, shape: &TopoShape, parameters: &MeshParameters)
        -> Result<MeshModel, String>;
}

/// Explores a shape for parts to be meshed — faces and free edges — visiting
/// each one through a [`ShapeVisitor`].
/// Source: `IMeshTools_ShapeExplorer::Accept`.
pub trait ShapeExplorer {
    /// Visits all faces and edges of the explored shape.
    fn accept(&self, visitor: &mut dyn ShapeVisitor);
}

/// Default [`ShapeExplorer`]: walks the shape's distinct faces (and the edges
/// belonging to each), then visits the faces themselves. Edges that are not
/// part of any face are visited first as free edges, mirroring
/// `IMeshTools_ShapeExplorer::Accept`.
pub struct TopoShapeExplorer {
    shape: TopoShape,
}

impl TopoShapeExplorer {
    /// Creates an explorer for the given shape.
    pub fn new(shape: TopoShape) -> Self {
        Self { shape }
    }

    /// The shape being explored.
    pub fn shape(&self) -> &TopoShape {
        &self.shape
    }
}

impl ShapeExplorer for TopoShapeExplorer {
    fn accept(&self, visitor: &mut dyn ShapeVisitor) {
        let shape = &self.shape;

        // Edges that appear in at least one face.
        let face_edge_ptrs: HashSet<usize> = faces_of(shape)
            .iter()
            .flat_map(|f| edges_of(&f.0))
            .map(|e| tshape_ptr(&e.0))
            .collect();

        // Free edges first (not bound to any face), visited forward.
        for e in edges_of(shape) {
            if !face_edge_ptrs.contains(&tshape_ptr(&e.0)) {
                visitor.visit_edge(&Edge(e.0.oriented(Orientation::Forward)));
            }
        }

        // Faces with their edges, forward-oriented. Edges shared between faces
        // are visited more than once; the visitor deduplicates them.
        for f in faces_of(shape) {
            for e in edges_of(&f.0) {
                visitor.visit_edge(&e);
            }
            visitor.visit_face(&Face(f.0.oriented(Orientation::Forward)));
        }
    }
}

/// Identifies which model-processing slot of a [`MeshContext`] is invoked.
#[derive(Debug, Clone, Copy)]
enum AlgoSlot {
    EdgeDiscret,
    ModelHealer,
    PreProcessor,
    FaceDiscret,
    PostProcessor,
}

/// Context of the BRepMesh algorithm: caches the discrete model and the
/// instances of the tools used to build and process it.
/// Source: `IMeshTools_Context`.
pub struct MeshContext {
    shape: TopoShape,
    parameters: MeshParameters,
    model: Option<MeshModel>,
    model_builder: Option<Box<dyn MeshBuilder>>,
    mesh_algo: Option<Box<dyn MeshAlgo>>,
    edge_discret: Option<Box<dyn ModelAlgo>>,
    model_healer: Option<Box<dyn ModelAlgo>>,
    pre_processor: Option<Box<dyn ModelAlgo>>,
    face_discret: Option<Box<dyn ModelAlgo>>,
    post_processor: Option<Box<dyn ModelAlgo>>,
}

impl MeshContext {
    /// Creates a context for the given shape and meshing parameters.
    pub fn new(shape: TopoShape, parameters: MeshParameters) -> Self {
        Self {
            shape,
            parameters,
            model: None,
            model_builder: None,
            mesh_algo: None,
            edge_discret: None,
            model_healer: None,
            pre_processor: None,
            face_discret: None,
            post_processor: None,
        }
    }

    /// The shape being meshed.
    pub fn shape(&self) -> &TopoShape {
        &self.shape
    }

    /// Meshing parameters.
    pub fn parameters(&self) -> &MeshParameters {
        &self.parameters
    }

    /// Mutable reference to the meshing parameters.
    pub fn change_parameters(&mut self) -> &mut MeshParameters {
        &mut self.parameters
    }

    /// The discrete model built for the shape, if any.
    pub fn model(&self) -> Option<&MeshModel> {
        self.model.as_ref()
    }

    /// Replaces the discrete model.
    pub fn set_model(&mut self, model: MeshModel) {
        self.model = Some(model);
    }

    /// Takes the discrete model out of the context.
    pub fn take_model(&mut self) -> Option<MeshModel> {
        self.model.take()
    }

    /// The assigned model builder, if any.
    pub fn model_builder(&self) -> Option<&dyn MeshBuilder> {
        self.model_builder.as_deref()
    }

    /// Assigns the model builder.
    pub fn set_model_builder(&mut self, builder: Box<dyn MeshBuilder>) {
        self.model_builder = Some(builder);
    }

    /// The assigned per-face mesh algorithm, if any.
    pub fn mesh_algo(&self) -> Option<&dyn MeshAlgo> {
        self.mesh_algo.as_deref()
    }

    /// Assigns the per-face mesh algorithm.
    pub fn set_mesh_algo(&mut self, algo: Box<dyn MeshAlgo>) {
        self.mesh_algo = Some(algo);
    }

    /// The assigned edge discretizer, if any.
    pub fn edge_discret(&self) -> Option<&dyn ModelAlgo> {
        self.edge_discret.as_deref()
    }

    /// Assigns the edge discretizer.
    pub fn set_edge_discret(&mut self, algo: Box<dyn ModelAlgo>) {
        self.edge_discret = Some(algo);
    }

    /// The assigned model healer, if any.
    pub fn model_healer(&self) -> Option<&dyn ModelAlgo> {
        self.model_healer.as_deref()
    }

    /// Assigns the model healer.
    pub fn set_model_healer(&mut self, algo: Box<dyn ModelAlgo>) {
        self.model_healer = Some(algo);
    }

    /// The assigned pre-processing algorithm, if any.
    pub fn pre_processor(&self) -> Option<&dyn ModelAlgo> {
        self.pre_processor.as_deref()
    }

    /// Assigns the pre-processing algorithm.
    pub fn set_pre_processor(&mut self, algo: Box<dyn ModelAlgo>) {
        self.pre_processor = Some(algo);
    }

    /// The assigned face discretizer, if any.
    pub fn face_discret(&self) -> Option<&dyn ModelAlgo> {
        self.face_discret.as_deref()
    }

    /// Assigns the face discretizer.
    pub fn set_face_discret(&mut self, algo: Box<dyn ModelAlgo>) {
        self.face_discret = Some(algo);
    }

    /// The assigned post-processing algorithm, if any.
    pub fn post_processor(&self) -> Option<&dyn ModelAlgo> {
        self.post_processor.as_deref()
    }

    /// Assigns the post-processing algorithm.
    pub fn set_post_processor(&mut self, algo: Box<dyn ModelAlgo>) {
        self.post_processor = Some(algo);
    }

    /// Builds the discrete model using the assigned model builder. Fails when
    /// no builder has been assigned.
    pub fn build_model(&mut self) -> Result<(), String> {
        let builder = self
            .model_builder
            .as_deref()
            .ok_or("MeshContext::build_model: no model builder assigned")?;
        let model = builder.build_model(&self.shape, &self.parameters)?;
        self.model = Some(model);
        Ok(())
    }

    /// Discretizes the edges of the model using the assigned edge discretizer.
    pub fn discretize_edges(&mut self) -> Result<(), String> {
        self.run(AlgoSlot::EdgeDiscret)
    }

    /// Heals the discrete model using the assigned healer.
    pub fn heal_model(&mut self) -> Result<(), String> {
        self.run(AlgoSlot::ModelHealer)
    }

    /// Pre-processes the discrete model (e.g. cleaning old triangulation).
    pub fn pre_process_model(&mut self) -> Result<(), String> {
        self.run(AlgoSlot::PreProcessor)
    }

    /// Meshes the faces of the discrete model using the assigned face discretizer.
    pub fn discretize_faces(&mut self) -> Result<(), String> {
        self.run(AlgoSlot::FaceDiscret)
    }

    /// Post-processes the discrete model using the assigned post-processor.
    pub fn post_process_model(&mut self) -> Result<(), String> {
        self.run(AlgoSlot::PostProcessor)
    }

    /// Cleans temporary context data: drops the cached model when
    /// `clean_model` is requested.
    pub fn clean(&mut self) {
        if self.parameters.clean_model {
            self.model = None;
        }
    }

    /// Runs the model algorithm in the given slot against the cached model.
    fn run(&mut self, slot: AlgoSlot) -> Result<(), String> {
        let algo = match slot {
            AlgoSlot::EdgeDiscret => self
                .edge_discret
                .as_mut()
                .ok_or("MeshContext: edge discretizer not assigned")?,
            AlgoSlot::ModelHealer => self
                .model_healer
                .as_mut()
                .ok_or("MeshContext: model healer not assigned")?,
            AlgoSlot::PreProcessor => self
                .pre_processor
                .as_mut()
                .ok_or("MeshContext: pre-processor not assigned")?,
            AlgoSlot::FaceDiscret => self
                .face_discret
                .as_mut()
                .ok_or("MeshContext: face discretizer not assigned")?,
            AlgoSlot::PostProcessor => self
                .post_processor
                .as_mut()
                .ok_or("MeshContext: post-processor not assigned")?,
        };
        let model = self.model.as_mut().ok_or("MeshContext: model not built")?;
        algo.perform(model, &self.parameters)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::BRepPrimBox;

    /// Test builder: creates an empty model for the shape.
    struct TestModelBuilder;

    impl MeshBuilder for TestModelBuilder {
        fn build_model(
            &self,
            shape: &TopoShape,
            _parameters: &MeshParameters,
        ) -> Result<MeshModel, String> {
            Ok(MeshModel::new(shape.clone()))
        }
    }

    /// No-op per-face mesh algorithm.
    struct NoopMeshAlgo;

    impl MeshAlgo for NoopMeshAlgo {
        fn perform(
            &mut self,
            _face: &mut MeshFace,
            _parameters: &MeshParameters,
        ) -> Result<(), String> {
            Ok(())
        }
    }

    /// No-op whole-model algorithm.
    struct NoopModelAlgo;

    impl ModelAlgo for NoopModelAlgo {
        fn perform(
            &mut self,
            _model: &mut MeshModel,
            _parameters: &MeshParameters,
        ) -> Result<(), String> {
            Ok(())
        }
    }

    /// Visitor that counts distinct visited edges and faces by `TShape` id.
    #[derive(Default)]
    struct CountingVisitor {
        edges: HashSet<usize>,
        faces: HashSet<usize>,
    }

    impl ShapeVisitor for CountingVisitor {
        fn visit_edge(&mut self, edge: &Edge) {
            self.edges.insert(tshape_ptr(&edge.0));
        }

        fn visit_face(&mut self, face: &Face) {
            self.faces.insert(tshape_ptr(&face.0));
        }
    }

    #[test]
    fn context_assembles_and_builds_model() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let mut ctx = MeshContext::new(b.solid.0, MeshParameters::default());

        // Defaults mirror OCCT's IMeshTools_Parameters.
        assert_eq!(ctx.parameters().angle, 0.5);
        assert_eq!(ctx.parameters().deflection, 0.001);
        assert!(ctx.parameters().clean_model);

        // Nothing assigned yet; build fails gracefully.
        assert!(ctx.model_builder().is_none());
        assert!(ctx.model().is_none());
        assert!(ctx.build_model().is_err());

        // Assemble the tools.
        ctx.set_model_builder(Box::new(TestModelBuilder));
        ctx.set_mesh_algo(Box::new(NoopMeshAlgo));
        ctx.set_edge_discret(Box::new(NoopModelAlgo));
        ctx.set_model_healer(Box::new(NoopModelAlgo));
        ctx.set_pre_processor(Box::new(NoopModelAlgo));
        ctx.set_face_discret(Box::new(NoopModelAlgo));
        ctx.set_post_processor(Box::new(NoopModelAlgo));
        assert!(ctx.model_builder().is_some());
        assert!(ctx.mesh_algo().is_some());

        // Build → cache a model, run a full pipeline pass.
        ctx.build_model().expect("build model");
        let model = ctx.model().expect("model cached");
        assert_eq!(model.faces_nb(), 0, "test builder creates no faces");

        // Model algorithms require a built model.
        ctx.discretize_edges().expect("edge discret");
        ctx.heal_model().expect("heal");
        ctx.pre_process_model().expect("pre-process");
        ctx.discretize_faces().expect("face discret");
        ctx.post_process_model().expect("post-process");

        // Mutating parameters through change_parameters.
        ctx.change_parameters().deflection = 0.01;
        assert_eq!(ctx.parameters().deflection, 0.01);

        // clean() drops the cached model because clean_model defaults to true.
        ctx.clean();
        assert!(ctx.model().is_none());
    }

    #[test]
    fn shape_explorer_visits_box_edges_and_faces() {
        let b = BRepPrimBox::make_box(1.0, 2.0, 3.0);
        let explorer = TopoShapeExplorer::new(b.solid.0);
        assert!(explorer.shape().is_solid());

        let mut visitor = CountingVisitor::default();
        explorer.accept(&mut visitor);

        assert_eq!(visitor.edges.len(), 12, "box has 12 distinct edges");
        assert_eq!(visitor.faces.len(), 6, "box has 6 faces");
    }

    #[test]
    fn shape_explorer_accepts_through_trait_object() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let explorer: Box<dyn ShapeExplorer> = Box::new(TopoShapeExplorer::new(b.solid.0));
        let mut visitor = CountingVisitor::default();
        explorer.accept(&mut visitor);
        assert_eq!(visitor.edges.len(), 12);
        assert_eq!(visitor.faces.len(), 6);
    }
}
