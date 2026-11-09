//! `BRepMesh_ModelPostProcessor` (`BRepMesh_ModelPostProcessor.cxx:29-238`,
//! `.hxx:31-63`).
//!
//! The model post-processing step `IMeshTools_MeshBuilder::Perform`
//! (`IMeshTools_MeshBuilder.cxx:63`) runs after `DiscretizeFaces`: the face
//! triangulations are already installed on the `BRep_TFace`s
//! (`BRepMesh_BaseMeshAlgo::commitSurfaceTriangulation`,
//! `BRepMesh_BaseMeshAlgo.cxx:229-241`), and this algorithm
//!
//! 1. commits the edge polygons - `Poly_Polygon3D` for free edges
//!    (`PolygonCommitter::commitPolygon3D`, `cxx:64-81`) and
//!    `Poly_PolygonOnTriangulation` for edges shared with faces
//!    (`commitPolygons`, `cxx:84-139`), the producer
//!    `HLRBRep_PolyAlgo::InitBiPoints`'s connexity step reads back through
//!    `BRep_Tool::PolygonOnTriangulation`; then
//! 2. recomputes each face's triangulation deflection
//!    (`DeflectionEstimator`, `cxx:166-202`), which needs both the polygons and
//!    the triangulation to exist (`BRepLib::UpdateDeflection`,
//!    `BRepLib.cxx:2793-2890`) and then records the parameters it was built for.
//!
//! The two `OSD_Parallel::For` loops (`cxx:227-234`) force the polygon pass to a
//! single thread (the comment at `cxx:226` says the edges share `TShape`s) and
//! run the deflection pass single-threaded when `InParallel` is off; both are
//! sequential here, which is that same order.
//!
//! UNPORTED / deviations:
//! * `performInternal`'s null-model answer (`cxx:220-224`): the model is a
//!   reference, so there is no null case to report.
//! * `Poly_PolygonOnTriangulation` node indices are stored 1-based exactly as
//!   OCCT stores them (`collectPolygon`, `cxx:148` copies
//!   `IMeshData_PCurve::Index` verbatim); every reader that indexes the crate's
//!   0-based `Triangulation::nodes` has to subtract one, see
//!   [`collect_polygon`].

use std::sync::Arc;

use occt_core::poly::{Polygon3D, PolyPolygonOnTriangulation, TriangulationParameters};

use crate::brep_lib_update_deflection::update_deflection;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::Face;

use super::context::ModelAlgo;
use super::data_model::{MeshEdge, MeshModel, MeshPCurve, MeshStatus};
use super::parameters::MeshParameters;

/// Commits 3D polygons and polygons on triangulations for corresponding edges.
/// Source: the file-local `PolygonCommitter` (`cxx:33-161`).
pub struct PolygonCommitter<'a> {
    model: &'a MeshModel,
}

impl<'a> PolygonCommitter<'a> {
    /// `PolygonCommitter(theModel)` (`cxx:36-39`).
    pub fn new(model: &'a MeshModel) -> Self {
        Self { model }
    }

    /// `operator()(theEdgeIndex)` (`cxx:42-61`).
    pub fn commit(&self, edge_index: usize) -> Result<(), String> {
        let d_edge = self.model.edge(edge_index)?;
        // `cxx:44-47`: an edge without a 3D discretization carries nothing.
        if d_edge.discretization().parameters_nb() == 0 {
            return Ok(());
        }
        if d_edge.is_free() {
            // `cxx:50-56`: a reused edge keeps the polygon it already has.
            if !d_edge.is_status(MeshStatus::REUSED) {
                self.commit_polygon3d(d_edge)?;
            }
        } else {
            self.commit_polygons(d_edge)?;
        }
        Ok(())
    }

    /// `commitPolygon3D(theDEdge)` (`cxx:64-81`): `Poly_Polygon3D` with the
    /// 3D curve's points and parameters, stored with the edge's deflection.
    fn commit_polygon3d(&self, d_edge: &MeshEdge) -> Result<(), String> {
        let curve = d_edge.discretization();
        let nb = curve.parameters_nb();
        let mut nodes = Vec::with_capacity(nb);
        let mut params = Vec::with_capacity(nb);
        for i in 0..nb {
            nodes.push(curve.point(i)?);
            params.push(curve.parameter(i)?);
        }
        let mut polygon = Polygon3D::new(nodes).with_params(params);
        // `cxx:79`: `aPoly3D->Deflection(theDEdge->GetDeflection())`.
        polygon.deflection = d_edge.deflection();
        // `cxx:81`: `BRepMesh_ShapeTool::UpdateEdge(theDEdge->GetEdge(), aPoly3D)`
        // -> `BRep_Builder::UpdateEdge(E, P, L)` (`BRep_Builder.cxx:751-786`).
        TopoBuilder::new().update_edge_polygon3d(d_edge.edge(), Some(Arc::new(polygon)));
        Ok(())
    }

    /// `commitPolygons(theDEdge)` (`cxx:84-139`): group the edge's pcurves by
    /// face, then store one polygon on each face's triangulation (two when the
    /// edge has a pair of pcurves on it - a seam).
    fn commit_polygons(&self, d_edge: &MeshEdge) -> Result<(), String> {
        // `cxx:87-105`: collect pcurves per face, skipping the faces that
        // failed or were reused.
        let mut groups: Vec<(usize, Vec<&MeshPCurve>)> = Vec::new();
        for pc_index in 0..d_edge.pcurves_nb() {
            let pcurve = d_edge.pcurve(pc_index)?;
            let face_index = pcurve.face();
            let face = self.model.face(face_index)?;
            if face.is_status(MeshStatus::FAILURE) || face.is_status(MeshStatus::REUSED) {
                continue;
            }
            match groups.iter_mut().find(|(index, _)| *index == face_index) {
                Some((_, pcurves)) => pcurves.push(pcurve),
                None => groups.push((face_index, vec![pcurve])),
            }
        }

        // `cxx:107-138`: commit the polygons of each face's triangulation.
        for (face_index, pcurves) in groups {
            let face = self.model.face(face_index)?;
            let Some(triangulation) = BRepTool::face_triangulation(face.face()) else {
                continue;
            };
            let Some(first) = pcurves.first() else {
                continue;
            };
            let deflection = d_edge.deflection();
            if pcurves.len() == 2 {
                // `cxx:120-129`: the closed-triangulation form stores both
                // pcurves (`BRep_Builder::UpdateEdge(E, P1, P2, T, L)`,
                // `BRep_Builder.cxx:837-883`).
                let polygon1 = collect_polygon(first, deflection)?;
                let polygon2 = collect_polygon(pcurves[pcurves.len() - 1], deflection)?;
                TopoBuilder::new().update_edge_polygons_on_triangulation(
                    d_edge.edge(),
                    Arc::new(polygon1),
                    Some(Arc::new(polygon2)),
                    triangulation,
                );
            } else {
                // `cxx:130-136`: the single-polygon form
                // (`BRep_Builder::UpdateEdge(E, P, T, L)`,
                // `BRep_Builder.cxx:790-832`).
                let polygon = collect_polygon(first, deflection)?;
                TopoBuilder::new().update_edge_polygon_on_triangulation(
                    d_edge.edge(),
                    Arc::new(polygon),
                    triangulation,
                );
            }
        }
        Ok(())
    }
}

/// `collectPolygon(thePCurve, theDeflection)` (`cxx:144-157`): the pcurve's
/// node indices and parameters as a `Poly_PolygonOnTriangulation`.
fn collect_polygon(
    pcurve: &MeshPCurve,
    deflection: f64,
) -> Result<PolyPolygonOnTriangulation, String> {
    let nb = pcurve.parameters_nb();
    let mut nodes = Vec::with_capacity(nb);
    let mut params = Vec::with_capacity(nb);
    for i in 0..nb {
        // `cxx:148`: `aNodes(i) = thePCurve->GetIndex(i - 1)`, stored verbatim
        // as the 1-based triangulation node index OCCT consumers expect.
        //
        // A pcurve whose points were never registered keeps the 0 that
        // `BRepMeshData_PCurve` initializes its indices with
        // (`BRepMeshData_PCurve.cxx:47`, `:56`) - OCCT stores it verbatim too;
        // the situation only arises for the wires `BRepMesh_BaseMeshAlgo::
        // initDataStructure` skips (`BRepMesh_BaseMeshAlgo.cxx:80-85`).
        nodes.push(pcurve.index(i)?.max(0) as usize);
        params.push(pcurve.parameter(i)?);
    }
    let mut polygon = PolyPolygonOnTriangulation::from_nodes_params(nodes, params)?;
    // `cxx:155`: `aPolygon->Deflection(theDeflection)`.
    polygon.set_deflection(deflection);
    Ok(polygon)
}

/// Estimates and updates deflection of triangulations for corresponding faces.
/// Source: the file-local `DeflectionEstimator` (`cxx:164-203`).
pub struct DeflectionEstimator<'a> {
    model: &'a MeshModel,
    /// The parameter set the triangulations were built for: OCCT keeps a
    /// `Poly_TriangulationParameters(theParams.Deflection, theParams.Angle,
    /// theParams.MinSize)` (`cxx:170-176`).
    parameters: TriangulationParameters,
}

impl<'a> DeflectionEstimator<'a> {
    /// `DeflectionEstimator(theModel, theParams)` (`cxx:169-177`).
    pub fn new(model: &'a MeshModel, params: &MeshParameters) -> Self {
        Self {
            model,
            parameters: TriangulationParameters::new(
                params.deflection,
                params.angle,
                params.min_size,
            ),
        }
    }

    /// `operator()(theFaceIndex)` (`cxx:179-197`).
    pub fn estimate(&self, face_index: usize) -> Result<(), String> {
        let d_face = self.model.face(face_index)?;
        if d_face.is_status(MeshStatus::FAILURE) || d_face.is_status(MeshStatus::REUSED) {
            return Ok(());
        }

        // `cxx:188`: `BRepLib::UpdateDeflection(aDFace->GetFace())`.
        update_deflection(&d_face.face().0);

        // `cxx:190-196`: record the parameters on the installed triangulation
        // (`Poly_Triangulation::Parameters(theParams)`).
        if BRepTool::face_triangulation(d_face.face()).is_some() {
            install_parameters(d_face.face(), self.parameters);
        }
        Ok(())
    }
}

/// `aTriangulation->Parameters(theParams)` (`BRepMesh_ModelPostProcessor.cxx:196`).
///
/// `Poly_Triangulation` is a shared object whose handle is the key of the
/// `BRep_PolygonOnTriangulation` representations written by
/// [`PolygonCommitter`] (`BRep_PolygonOnTriangulation.cxx:36-49`), so the
/// parameter set is rewritten in place through [`Triangulation::set_parameters`]
/// and the handle itself never changes.
fn install_parameters(face: &Face, parameters: TriangulationParameters) {
    if let Some(triangulation) = BRepTool::face_triangulation(face) {
        triangulation.set_parameters(Some(parameters));
    }
}

/// `BRepMesh_ModelPostProcessor` (`BRepMesh_ModelPostProcessor.hxx:31-63`).
///
/// `Perform` (`cxx:215-235`) never touches the discrete model itself: it writes
/// onto the topological shapes through `BRep_Builder`, so the `&mut MeshModel`
/// of [`ModelAlgo::perform`] is only borrowed to read the edge/face tables.
pub struct BRepMeshModelPostProcessor;

impl Default for BRepMeshModelPostProcessor {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepMeshModelPostProcessor {
    /// Default constructor (`cxx:216`).
    pub fn new() -> Self {
        Self
    }
}

impl ModelAlgo for BRepMeshModelPostProcessor {
    /// `performInternal(theModel, theParameters, theRange)`
    /// (`cxx:215-235`).
    fn perform(&mut self, model: &mut MeshModel, parameters: &MeshParameters) -> Result<(), String> {
        // `cxx:227-229`: PolygonCommitter over every edge, forced single
        // threaded because edges share `TShape`s (`cxx:226`).
        let committer = PolygonCommitter::new(model);
        for edge_index in 0..model.edges_nb() {
            committer.commit(edge_index)?;
        }

        // `cxx:231-234`: DeflectionEstimator over every face.
        let estimator = DeflectionEstimator::new(model, parameters);
        for face_index in 0..model.faces_nb() {
            estimator.estimate(face_index)?;
        }
        Ok(())
    }
}
