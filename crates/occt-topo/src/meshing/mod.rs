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
pub mod geomlib_norm;    // GeomLib::NormEstim
pub mod deflection;      // BRepMesh_Deflection
pub mod shape_tool;      // BRepMesh_ShapeTool + ShapeVisitor
pub mod edge_discret;    // BRepMesh_EdgeDiscret + EdgeParameterProvider + CurveTessellator
pub mod face_discret;    // BRepMesh_FaceDiscret + FaceChecker + Classifier
pub mod incremental_mesh;// BRepMesh_IncrementalMesh + DiscretRoot
pub mod model_builder;   // BRepMesh_ModelBuilder + ModelPreProcessor
pub mod wire_order;      // ShapeAnalysis_WireOrder (2D pcurve mode)

// Wave 2 — Delaunay core.
pub mod delaun_types;    // BRepMesh_{Vertex,Triangle,Circle,Edge} + OrientedEdge + PairOfIndex
pub mod delaun_data;     // BRepMesh_DataStructureOfDelaun + SelectorOfDataStructureOfDelaun
pub mod delaun;          // BRepMesh_Delaun (2D UV-space incremental Delaunay)
pub mod delaun_index;    // BRepMesh_{VertexTool,VertexInspector,CircleTool,CircleInspector}
pub mod mesh_tool;       // BRepMesh_MeshTool

// Wave 3 — mesh algorithms + deflection control.
pub mod mesh_algo;         // BRepMesh_{BaseMeshAlgo,ConstrainedBaseMeshAlgo,CustomBaseMeshAlgo,CustomDelaunayBaseMeshAlgo,DelaunayBaseMeshAlgo}
pub mod node_insertion;    // BRepMesh_{NodeInsertionMeshAlgo,DelaunayNodeInsertionMeshAlgo}
pub mod deflection_control;// BRepMesh_DelaunayDeflectionControlMeshAlgo
pub mod triangulator;      // BRepMesh_Triangulator
pub mod fast_discret;      // BRepMesh_FastDiscret

// Wave 4 — range splitters, healing, factories, Delabella + pipeline integration.
pub mod range_splitter;    // BRepMesh_{DefaultRangeSplitter,UVParamRangeSplitter,UndefinedRangeSplitter,BoundaryParamsRangeSplitter, CylinderRangeSplitter,ConeRangeSplitter,SphereRangeSplitter,TorusRangeSplitter,NURBSRangeSplitter,ExtrusionRangeSplitter}
pub mod model_healer;      // BRepMesh_{ModelHealer,ModelPostProcessor}
pub mod factories;         // BRepMesh_{DiscretFactory,DiscretAlgoFactory,MeshAlgoFactory,IncrementalMeshFactory}
pub mod delabella;         // BRepMesh_{DelabellaBaseMeshAlgo,DelabellaMeshAlgoFactory}
pub mod degree_of_freedom; // BRepMesh_DegreeOfFreedom
