//! BRepMesh toolkit depth port (OCCT `src/ModelingAlgorithms/TKMesh`).
//!
//! Modern incremental meshing pipeline, class-by-class port of the 59-class
//! `BRepMesh` toolkit + `BRepMeshData`(6) + `IMeshData`(17) + `IMeshTools`(14).
//! Reference: `D:\source\occt-src\src\ModelingAlgorithms\TKMesh\`.
//!
//! Pipeline: `IncrementalMesh::perform` → `ModelBuilder`/`ModelPreProcessor`
//! (shape → `MeshModel` with edges/faces/curves/pcurves) → edge/face
//! discretization with deflection (`EdgeDiscret`/`FaceDiscret`/`GeomTool`) →
//! 2D UV-space Delaunay triangulation (`Delaun`) → deflection-controlled
//! refinement (`DelaunayDeflectionControlMeshAlgo`) → healing
//! (`ModelHealer`/`ModelPostProcessor`).

// Wave 1 — framework + data + geometry + pipeline entry.
pub mod parameters;      // IMeshTools_Parameters
pub mod data_model;      // BRepMeshData::{Model,Edge,Face,Wire,Curve,PCurve} + IMeshData::{Status,Types}
pub mod context;         // IMeshTools_Context + MeshAlgo/ModelAlgo/MeshBuilder/ShapeExplorer traits
pub mod geom_tool;       // BRepMesh_GeomTool
pub mod deflection;      // BRepMesh_Deflection
pub mod shape_tool;      // BRepMesh_ShapeTool + ShapeVisitor
pub mod edge_discret;    // BRepMesh_EdgeDiscret + EdgeParameterProvider + CurveTessellator
pub mod face_discret;    // BRepMesh_FaceDiscret + FaceChecker + Classifier
pub mod incremental_mesh;// BRepMesh_IncrementalMesh + DiscretRoot
pub mod model_builder;   // BRepMesh_ModelBuilder + ModelPreProcessor
