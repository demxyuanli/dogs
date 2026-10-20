// Shared include set and cross-file declarations for the OCCT ground-truth probe.
// The probe used to be one translation unit; it is split by command group, so every
// declaration below is one that `main` or another group calls. Code bodies are
// unchanged, only moved.
#pragma once

// OCCT ground-truth probe (T-69/T-80/T-87/T-05 oracle).
// Reads a STEP file with OCCT 8.0.0 and prints what the port must reproduce:
//   * BRepCheck_Analyzer verdict
//   * BRepGProp::VolumeProperties (the OCCT volume, analytic by default)
//   * per-face surface type, orientation, composed normal and outward sign
//   * cylinder/spline frames, so the port's reconstruction can be compared
#include <BRepAdaptor_Surface.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
// T-99: `BRepMesh_DefaultRangeSplitter` exposes `GetRangeU/GetRangeV/GetDelta/
// GetToleranceUV` publicly, which is all the diagnostic below needs. Its
// `computeLengthU/V` are declared `private` in the shipped header and their symbols
// are emitted with private access, so they cannot be called from outside without
// either editing the installed header or recompiling the library - the same length
// values are already reported by the `--uvsum` dump, so they are not repeated here.
#include <BRepMesh_GeomTool.hxx>
#include <BRepMesh_DefaultRangeSplitter.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepMesh_ModelBuilder.hxx>
#include <BRepMesh_EdgeDiscret.hxx>
#include <BRepMesh_ModelHealer.hxx>
#include <ShapeAnalysis_Wire.hxx>
#include <ShapeAnalysis_WireOrder.hxx>
#include <ShapeExtend_WireData.hxx>
#include <ShapeFix_Wire.hxx>
#include <ShapeFix_Face.hxx>
#include <ShapeBuild_ReShape.hxx>
#include <ShapeAnalysis_Surface.hxx>
#include <ShapeExtend.hxx>
#include <IMeshData_Model.hxx>
#include <IMeshData_Face.hxx>
#include <IMeshData_Wire.hxx>
#include <IMeshData_Edge.hxx>
#include <IMeshData_PCurve.hxx>
#include <IMeshTools_Parameters.hxx>
#include <Message_ProgressRange.hxx>
// T-104: OCCT's own node-insertion algo, subclassed to dump the structure it
// hands to `BRepMesh_Delaun` (`--structdump`). Kept after the `IMeshData_*`
// includes because `BRepMesh_NodeInsertionMeshAlgo.hxx` needs those types.
#include <BRepMesh_DelaunayNodeInsertionMeshAlgo.hxx>
#include <BRepMesh_DelaunayBaseMeshAlgo.hxx>
// T-123: the `--algospy` replay subclasses the real node-insertion algo, so the face
// discretiser that owns the algo factory, the factory itself, the algorithm base and the
// Delaunay structure the dump reports have to be complete types.
#include <BRepMesh_FaceDiscret.hxx>
#include <BRepMesh_MeshAlgoFactory.hxx>
#include <BRepMesh_DataStructureOfDelaun.hxx>
// TEMP DCPASS: the per-pass instrumentation below needs the private loop state
// (`myControlNodes` / `myCouplesMap` / `myMaxSqDeflection`) of the template
// header, so it is included with access control disabled. The class layout is
// untouched (access specifiers do not affect it) and the instantiation is local.
#define private public
#include <BRepMesh_DelaunayDeflectionControlMeshAlgo.hxx>
#undef private
#include <BRepMesh_DefaultRangeSplitter.hxx>
#include <BRepMesh_BoundaryParamsRangeSplitter.hxx>
#include <BRepMesh_ExtrusionRangeSplitter.hxx>
#include <BRepMesh_NURBSRangeSplitter.hxx>
#include <BRepMesh_UndefinedRangeSplitter.hxx>
#include <IGESControl_Controller.hxx>
#include <IGESControl_Reader.hxx>
#include <IGESControl_Writer.hxx>
#include <BRepMesh_TorusRangeSplitter.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <BRepMesh_SphereRangeSplitter.hxx>
#include <BRepMesh_CylinderRangeSplitter.hxx>
#include <BRepMesh_ConeRangeSplitter.hxx>
#include <Precision.hxx>
#include <Poly_Triangulation.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <cstdlib>
#include <filesystem>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <BRepGProp_Face.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <STEPControl_Reader.hxx>
#include <XSControl_WorkSession.hxx>
#include <XSControl_TransferReader.hxx>
#include <Interface_InterfaceModel.hxx>
#include <Interface_Static.hxx>
#include <ShapeProcess.hxx>
#include <ShapeAnalysis.hxx>
#include <ShapeAnalysis_Edge.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <ShapeAnalysis_Surface.hxx>
#include <ShapeConstruct_ProjectCurveOnSurface.hxx>
#include <DE_ShapeFixParameters.hxx>
#include <StepData_StepModel.hxx>
#include <BRepGProp_Domain.hxx>
#include <BRepGProp_Vinert.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BOPAlgo_PaveFiller.hxx>
#include <BOPDS_DS.hxx>
#include <BOPAlgo_Builder.hxx>
#include <BOPTools_AlgoTools.hxx>
#include <BOPTools_CoupleOfShape.hxx>
#include <IntTools_Context.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BOPTools_AlgoTools3D.hxx>
#include <BOPTools_AlgoTools2D.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <Geom_Plane.hxx>
#include <BRepOffset.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <BRep_Tool.hxx>
#include <BRepTools.hxx>
#include <TopTools_ListOfShape.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <StepShape_FaceSurface.hxx>
#include <string>
#include <TopExp_Explorer.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopExp.hxx>
#include <TopoDS_Vertex.hxx>
// T-191: the `--pcdef` dump needs the concrete 2-D curve types behind a
// `BRep_Tool::CurveOnSurface` handle.
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom_Curve.hxx>
#include <Geom_BSplineCurve.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
// T-191: the `--spcheck` replay of `XSAlgo_ShapeProcessor::CheckPCurve` needs
// the pcurve-fixer and the temporary-edge builder.
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <ShapeBuild_Edge.hxx>
#include <ShapeFix_Edge.hxx>
#include <BRep_Builder.hxx>

#include <Geom2d_Curve.hxx>
#include <Geom2d_Line.hxx>
#include <Bnd_Box2d.hxx>
// TEMP FCLS2D: `BRepTopAdaptor_FClass2d` boundary ring / state ground truth.
#include <BRepTopAdaptor_FClass2d.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2dAdaptor_Curve.hxx>
#include <gp_Pnt2d.hxx>
#include <algorithm>
#include <cmath>
#include <map>
#include <unordered_map>
#include <vector>

#include <cstdio>
#include <cstdlib>
#include <iostream>

// --- TEMP T-99 feasibility probe (`--delauncheck`) -------------------------------
// §9.365 established that OCCT meshes all 226 faces (min 2 triangles) while the port's
// `Delaun::add_vertices` yields 0 domain elements on 16 of them, so the next
// measurement must be OCCT's *Delaunay state*, not its product. The node-insertion
// algorithm is templated and not directly hookable, so the clean route is to feed
// `BRepMesh_Delaun` the port's own registered points. This block only proves that the
// class is reachable/linkable from this probe's link set; it constructs a Delaun over
// a tiny synthetic vertex array and reports the resulting structure size.
#include <BRepMesh_CircleTool.hxx>
#include <BRepMesh_Delaun.hxx>
#include <BRepMesh_Vertex.hxx>
#include <BRepMesh_DataStructureOfDelaun.hxx>
#include <BRepMesh_VertexTool.hxx>
#include <BRepMesh_Context.hxx>
#include <IMeshTools_Context.hxx>
#include <NCollection_IncAllocator.hxx>
#include <IMeshData_Types.hxx>
#include <array>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <optional>

// Argument/shape context handed to the `argv[2]` dispatch groups below. The shape is
// the imported root (`OneShape`); the reader is the STEP reader that produced it
// (`--faceids` / `--entities` read its transfer model).
struct ProbeCtx
{
  int                 argc;
  char**              argv;
  const TopoDS_Shape& shape;
  STEPControl_Reader& reader;
  bool                frames;
};

// `argv[2]`-dispatched groups. `std::nullopt` means `theCtx` does not ask for any
// mode of that group; otherwise the returned value is the process exit code.
std::optional<int> runMeshOracleModes(const ProbeCtx& theCtx);
std::optional<int> runPcurveModes(const ProbeCtx& theCtx);
std::optional<int> runUvModes(const ProbeCtx& theCtx);
std::optional<int> runFaceIdModes(const ProbeCtx& theCtx);
std::optional<int> runEdgeFaceModes(const ProbeCtx& theCtx);
std::optional<int> runFixFaceModes(const ProbeCtx& theCtx);
std::optional<int> runTopologyReportModes(const ProbeCtx& theCtx);
std::optional<int> runBooleanModes(const ProbeCtx& theCtx);
std::optional<int> runClassifyModes(const ProbeCtx& theCtx);
// Fallback report (`BRepCheck`, volume, per-face type/orientation/edges).
int                runDefaultReport(const ProbeCtx& theCtx);

// Structural dumps and per-stage replays called from `main`'s `run*` loop. The
// definition of each lives in the group file named after the command it serves.
int runDelaunCheck();
int runDelaunFeed(const std::string& theDir);
int runDelaunStruct(const std::string& theDir);
int runDelaunSteps(const std::string& theDir, const int theFace);
int runDelaunWhy(const std::string& theDir, const int theFace, const int theStep);
int runBoundary(const TopoDS_Shape& theShape, const std::string& theDir);
int runLenSamples(const TopoDS_Shape& theShape, const int theFaceIdx,
                  const double u0, const double u1, const double v0, const double v1);
int runStructDump(const TopoDS_Shape& theShape, const int theFaceIdx,
                  const bool theFactoryAlgo = false, const bool theSrfNodes = false,
                  const double theDefl = 0.0, const double theAngle = 0.0);
int runSrfParams(const TopoDS_Shape& theShape, const int theFaceIdx,
                 const double theDeflection, const double theAngle);
int runFaceStats(const TopoDS_Shape& theShape, const double theDeflection,
                 const double theAngle);
int runFaceDump(const TopoDS_Shape& theShape, const int theFaceIdx,
                const double theDeflection, const double theAngle);
int runEdgeDu(const TopoDS_Shape& theShape, const int theFaceIdx,
              const double theDeflection, const double theAngle);
int runTorusDbg(const TopoDS_Shape& theShape, const double theDeflection,
                const double theAngle);
int runHealerDbg(const TopoDS_Shape& theShape,
                 const double        theDeflection,
                 const double        theAngle,
                 const double*       theBBox);
int runCtxDbg(const TopoDS_Shape& theShape,
              const double        theDeflection,
              const double        theAngle,
              const double*       theBBox);
int runRangeReplay(const TopoDS_Shape& theShape,
                   const int           theFaceIdx,
                   const double        theDeflection,
                   const double        theAngle);
int runAlgoSpy(const TopoDS_Shape& theShape,
               const double        theDeflection,
               const double        theAngle,
               const int*          theFaces,
               const int           theFacesNb);
int runOrderDbg(const TopoDS_Shape& theShape, const double* theBBox);
int runWireDbg(const TopoDS_Shape& theShape,
               const double        theDeflection,
               const double        theAngle,
               const double*       theBBox,
               const int           theFaceIdx = -1);
int runIges(const TopoDS_Shape& theShape, const char* thePath);
int runIgesPc(const TopoDS_Shape& theShape, const char* thePath);
int runOffsetSurf(const TopoDS_Shape& theShape, const double theOffset);
int runExtSS();
