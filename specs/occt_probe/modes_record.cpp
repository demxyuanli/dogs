// Record mode: dump a shape as the RECORD v1 text the port's `record_dump`
// example also emits, so a parity run is a field-by-field text diff.
//
// Both halves of the comparison must start from the same *source-level*
// pipeline, never from a model file produced by an external application:
//
//   occt_probe <in.step> --record <out.txt> [defl] [angle]
//       read `in.step` with `STEPControl_Reader` and dump the record.
//   occt_probe <in.step> --record-rt iges <mid.iges> <out.txt> [defl] [angle]
//       read, write with `IGESControl_Writer`, read the written file back with
//       `IGESControl_Reader`, dump the record of the re-read shape. This is the
//       OCCT half of the write-side round trip: the writer and the reader are
//       both compiled from the OCCT sources in this tree's link set.
//   `--record-rt step` is the same with `STEPControl_Writer` /
//   `STEPControl_Reader`.
//
// The port half is `crates/occt-topo/examples/record_dump.rs`, which performs
// the identical pipeline with the port's reader/writer.
//
// Fields, one per line, in this order (both sides):
//   RECORD v1                     format tag
//   STAGE read|rt-iges|rt-step    which pipeline produced the shape
//   INPUT <path>
//   TYPE <Vertex|Edge|...>        root shape type
//   COUNTS V E W F SH SO CS C     `TopExp::MapShapes` census per type
//   VALID 0|1                     `BRepCheck_Analyzer`
//   VOLUME / AREA / CENTROID      `BRepGProp` on the whole shape
//   BBOX x0 y0 z0 x1 y1 z1        `BRepBndLib::Add(shape, box, false)`
//   FACE type= area= wires= edges= bbox=   one per face, sorted by text so the
//                                 traversal order does not matter
//   MESH nv nf area defl angle    `BRepMesh_IncrementalMesh` at the given
//                                 parameters, summed over all faces
//   END
#include "occt_probe_common.h"

#include <STEPControl_Writer.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <BRep_TEdge.hxx>
#include <BRep_CurveRepresentation.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <gp_Trsf.hxx>
#include <algorithm>
#include <iomanip>

namespace {

std::string shape_type_name(const TopAbs_ShapeEnum theType)
{
  switch (theType)
  {
    case TopAbs_VERTEX: return "Vertex";
    case TopAbs_EDGE: return "Edge";
    case TopAbs_WIRE: return "Wire";
    case TopAbs_FACE: return "Face";
    case TopAbs_SHELL: return "Shell";
    case TopAbs_SOLID: return "Solid";
    case TopAbs_COMPSOLID: return "CompSolid";
    case TopAbs_COMPOUND: return "Compound";
    default: return "Shape";
  }
}

// `GeomAbs_SurfaceType` name, spelled the way the port's
// `meshing::range_splitter::MeshSurfaceType` is printed by `tname`, so a face
// line can be compared as text.
const char* surface_type_name(const GeomAbs_SurfaceType theType)
{
  switch (theType)
  {
    case GeomAbs_Plane: return "Plane";
    case GeomAbs_Cylinder: return "Cylinder";
    case GeomAbs_Cone: return "Cone";
    case GeomAbs_Sphere: return "Sphere";
    case GeomAbs_Torus: return "Torus";
    case GeomAbs_BezierSurface: return "Bezier";
    case GeomAbs_BSplineSurface: return "BSpline";
    case GeomAbs_SurfaceOfRevolution: return "Rev";
    case GeomAbs_SurfaceOfExtrusion: return "Extr";
    case GeomAbs_OffsetSurface: return "Offset";
    case GeomAbs_OtherSurface: return "Other";
    default: return "Other";
  }
}

int count_of_type(const TopoDS_Shape& theShape, const TopAbs_ShapeEnum theType)
{
  TopTools_IndexedMapOfShape aMap;
  TopExp::MapShapes(theShape, theType, aMap);
  return aMap.Extent();
}

std::string fmt_g17(const double theValue)
{
  char aBuf[64];
  std::snprintf(aBuf, sizeof(aBuf), "%.17g", theValue);
  return std::string(aBuf);
}

std::string fmt_f9(const double theValue)
{
  char aBuf[64];
  std::snprintf(aBuf, sizeof(aBuf), "%.9f", theValue);
  return std::string(aBuf);
}

struct FaceRow
{
  std::string type;
  double      area;
  int         wires;
  int         edges;
  double      bbox[6];
  std::string text;
};

void bbox_of(const TopoDS_Shape& theShape, double* theOut)
{
  Bnd_Box aBox;
  BRepBndLib::Add(theShape, aBox, false);
  for (int k = 0; k < 6; ++k)
  {
    theOut[k] = 0.0;
  }
  if (!aBox.IsVoid())
  {
    aBox.Get(theOut[0], theOut[1], theOut[2], theOut[3], theOut[4], theOut[5]);
  }
}

// Sum of the 3D triangle areas of every face triangulation.
double mesh_area(const TopoDS_Shape& theShape, int& theNodes, int& theTriangles)
{
  theNodes = 0;
  theTriangles = 0;
  double anArea = 0.0;
  for (TopExp_Explorer anExp(theShape, TopAbs_FACE); anExp.More(); anExp.Next())
  {
    const TopoDS_Face& aFace = TopoDS::Face(anExp.Current());
    TopLoc_Location    aLoc;
    const occ::handle<Poly_Triangulation>& aTri = BRep_Tool::Triangulation(aFace, aLoc);
    if (aTri.IsNull())
    {
      continue;
    }
    theNodes += aTri->NbNodes();
    theTriangles += aTri->NbTriangles();
    const gp_Trsf& aTrsf = aLoc.Transformation();
    for (int t = 1; t <= aTri->NbTriangles(); ++t)
    {
      Standard_Integer n1 = 0, n2 = 0, n3 = 0;
      aTri->Triangle(t).Get(n1, n2, n3);
      const gp_Pnt p1 = aTri->Node(n1).Transformed(aTrsf);
      const gp_Pnt p2 = aTri->Node(n2).Transformed(aTrsf);
      const gp_Pnt p3 = aTri->Node(n3).Transformed(aTrsf);
      anArea += 0.5 * gp_Vec(p1, p2).Crossed(gp_Vec(p1, p3)).Magnitude();
    }
  }
  return anArea;
}

int write_record(const TopoDS_Shape& theShape,
                 const char*         theStage,
                 const char*         theInput,
                 const char*         theOutPath,
                 const double        theDeflection,
                 const double        theAngle)
{
  std::ofstream anOut(theOutPath);
  if (!anOut)
  {
    std::cerr << "cannot write " << theOutPath << "\n";
    return 6;
  }
  anOut << "RECORD v1\n";
  anOut << "STAGE " << theStage << "\n";
  anOut << "INPUT " << theInput << "\n";
  anOut << "TYPE " << shape_type_name(theShape.ShapeType()) << "\n";
  anOut << "COUNTS"
        << " V=" << count_of_type(theShape, TopAbs_VERTEX)
        << " E=" << count_of_type(theShape, TopAbs_EDGE)
        << " W=" << count_of_type(theShape, TopAbs_WIRE)
        << " F=" << count_of_type(theShape, TopAbs_FACE)
        << " SH=" << count_of_type(theShape, TopAbs_SHELL)
        << " SO=" << count_of_type(theShape, TopAbs_SOLID)
        << " CS=" << count_of_type(theShape, TopAbs_COMPSOLID)
        << " C=" << count_of_type(theShape, TopAbs_COMPOUND) << "\n";
  anOut << "VALID " << (BRepCheck_Analyzer(theShape).IsValid() ? 1 : 0) << "\n";

  {
    GProp_GProps aVol;
    BRepGProp::VolumeProperties(theShape, aVol);
    const gp_Pnt aC = aVol.CentreOfMass();
    anOut << "VOLUME " << fmt_g17(aVol.Mass()) << "\n";
    anOut << "CENTROID " << fmt_g17(aC.X()) << " " << fmt_g17(aC.Y()) << " " << fmt_g17(aC.Z())
          << "\n";
    GProp_GProps aSurf;
    BRepGProp::SurfaceProperties(theShape, aSurf);
    anOut << "AREA " << fmt_g17(aSurf.Mass()) << "\n";
  }

  double aBox[6];
  bbox_of(theShape, aBox);
  anOut << "BBOX";
  for (int k = 0; k < 6; ++k)
  {
    anOut << " " << fmt_f9(aBox[k]);
  }
  anOut << "\n";

  std::vector<FaceRow> aRows;
  for (TopExp_Explorer anExp(theShape, TopAbs_FACE); anExp.More(); anExp.Next())
  {
    const TopoDS_Face& aFace = TopoDS::Face(anExp.Current());
    FaceRow            aRow;
    const BRepAdaptor_Surface aAdaptor(aFace, Standard_False);
    aRow.type = surface_type_name(aAdaptor.GetType());
    GProp_GProps aProps;
    BRepGProp::SurfaceProperties(aFace, aProps);
    aRow.area = aProps.Mass();
    int aWireCount = 0;
    int anEdgeCount = 0;
    for (TopExp_Explorer aWireExp(aFace, TopAbs_WIRE); aWireExp.More(); aWireExp.Next())
    {
      ++aWireCount;
      for (TopExp_Explorer anEdgeExp(aWireExp.Current(), TopAbs_EDGE); anEdgeExp.More();
           anEdgeExp.Next())
      {
        ++anEdgeCount;
      }
    }
    aRow.wires = aWireCount;
    aRow.edges = anEdgeCount;
    bbox_of(aFace, aRow.bbox);
    std::ostringstream aText;
    aText << "FACE type=" << aRow.type << " area=" << fmt_g17(aRow.area)
          << " wires=" << aRow.wires << " edges=" << aRow.edges << " bbox=";
    for (int k = 0; k < 6; ++k)
    {
      aText << (k == 0 ? "" : " ") << fmt_f9(aRow.bbox[k]);
    }
    aRow.text = aText.str();
    aRows.push_back(aRow);
  }
  // Face traversal order differs between the two sides, so the comparison is a
  // multiset comparison: sort the rendered lines.
  std::sort(aRows.begin(), aRows.end(),
            [](const FaceRow& a, const FaceRow& b) { return a.text < b.text; });
  for (const FaceRow& aRow : aRows)
  {
    anOut << aRow.text << "\n";
  }

  {
    BRepMesh_IncrementalMesh aMesher(theShape, theDeflection, Standard_False, theAngle);
    int anNodes = 0;
    int aTris = 0;
    const double anArea = mesh_area(theShape, anNodes, aTris);
    anOut << "MESH nv=" << anNodes << " nf=" << aTris << " area=" << fmt_g17(anArea)
          << " defl=" << fmt_g17(theDeflection) << " angle=" << fmt_g17(theAngle) << "\n";
  }
  anOut << "END\n";
  anOut.close();
  std::cout << "RECORD " << theStage << " faces=" << aRows.size() << " out=" << theOutPath
            << std::endl;
  return 0;
}

bool read_step_shape(const char* thePath, TopoDS_Shape& theShape)
{
  STEPControl_Reader aReader;
  if (aReader.ReadFile(thePath) != IFSelect_RetDone)
  {
    std::cerr << "read failed: " << thePath << "\n";
    return false;
  }
  aReader.TransferRoots();
  theShape = aReader.OneShape();
  return !theShape.IsNull();
}

bool read_iges_shape(const char* thePath, TopoDS_Shape& theShape)
{
  IGESControl_Controller::Init();
  IGESControl_Reader aReader;
  if (aReader.ReadFile(thePath) != IFSelect_RetDone)
  {
    std::cerr << "read failed: " << thePath << "\n";
    return false;
  }
  aReader.TransferRoots();
  theShape = aReader.OneShape();
  return !theShape.IsNull();
}

bool write_iges_shape(const TopoDS_Shape& theShape, const char* thePath)
{
  IGESControl_Controller::Init();
  IGESControl_Writer aWriter("MM", 0);
  aWriter.AddShape(theShape);
  aWriter.ComputeModel();
  return aWriter.Write(thePath) != Standard_False;
}

bool write_step_shape(const TopoDS_Shape& theShape, const char* thePath)
{
  STEPControl_Writer aWriter;
  if (aWriter.Transfer(theShape, STEPControl_AsIs) != IFSelect_RetDone)
  {
    return false;
  }
  return aWriter.Write(thePath) == IFSelect_RetDone;
}

} // namespace

// `--edgerep`: per-edge representation census. `BRepCheck_Edge::Minimum`
// (`BRepCheck_Edge.cxx:74-148`) decides `No3DCurve` / `InvalidDegeneratedFlag`
// purely from the `BRep_TEdge` representation list, so the port's validity
// divergence on degenerated edges can only be read from this table:
//   EDGEREP n=<k> deg=<0|1> spar=<0|1> srange=<0|1> c3d=<null|nonnull>
//           cos=<number of CurveOnSurface reps> range=[f,l]
int runEdgeRep(const TopoDS_Shape& theShape)
{
  TopTools_IndexedMapOfShape anEdges;
  TopExp::MapShapes(theShape, TopAbs_EDGE, anEdges);
  int n = 0;
  for (int i = 1; i <= anEdges.Extent(); ++i)
  {
    const TopoDS_Edge& anEdge = TopoDS::Edge(anEdges(i));
    const occ::handle<BRep_TEdge>& aTE =
      *((occ::handle<BRep_TEdge>*)&anEdge.TShape());
    int  aCos = 0;
    bool aHasCurve3DRep = false;
    bool aCurve3DNull = true;
    for (NCollection_List<occ::handle<BRep_CurveRepresentation>>::Iterator itcr(aTE->Curves());
         itcr.More(); itcr.Next())
    {
      const occ::handle<BRep_CurveRepresentation>& aCr = itcr.Value();
      if (aCr->IsCurve3D())
      {
        aHasCurve3DRep = true;
        if (!aCr->Curve3D().IsNull())
        {
          aCurve3DNull = false;
        }
      }
      else if (aCr->IsCurveOnSurface())
      {
        ++aCos;
      }
    }
    std::cout << "EDGEREP n=" << (++n) << " deg=" << (aTE->Degenerated() ? 1 : 0)
              << " spar=" << (aTE->SameParameter() ? 1 : 0)
              << " srange=" << (aTE->SameRange() ? 1 : 0)
              << " c3drep=" << (aHasCurve3DRep ? 1 : 0)
              << " c3dnull=" << (aCurve3DNull ? 1 : 0)
              << " cos=" << aCos << std::endl;
  }
  return 0;
}

// `--vertexrep`: vertex census. `TopExp::MapShapes` (`TopTools_ShapeMapHasher`
// = `IsSame`) counts one entry per vertex TShape, so a junction one side shares
// and the other leaves split can only be read off this table:
//   VERTEX n=<k> p=[x,y,z] tol=<tol>
//   WIRE n=<k> edges=[edge n,...]
//   EDGE n=<k> first=<vertex n|-1> last=<vertex n|-1>
int runVertexRep(const TopoDS_Shape& theShape)
{
  TopTools_IndexedMapOfShape aVerts;
  TopExp::MapShapes(theShape, TopAbs_VERTEX, aVerts);
  for (int i = 1; i <= aVerts.Extent(); ++i)
  {
    const TopoDS_Vertex& aV = TopoDS::Vertex(aVerts(i));
    const gp_Pnt         aP = BRep_Tool::Pnt(aV);
    std::cout << "VERTEX n=" << i << " p=[" << aP.X() << "," << aP.Y() << "," << aP.Z()
              << "] tol=" << BRep_Tool::Tolerance(aV) << std::endl;
  }
  TopTools_IndexedMapOfShape aWires;
  TopExp::MapShapes(theShape, TopAbs_WIRE, aWires);
  TopTools_IndexedMapOfShape anEdges;
  TopExp::MapShapes(theShape, TopAbs_EDGE, anEdges);
  for (int i = 1; i <= aWires.Extent(); ++i)
  {
    std::cout << "WIRE n=" << i << " edges=[";
    bool aFirst = true;
    for (TopoDS_Iterator it(aWires(i), false, true); it.More(); it.Next())
    {
      if (!aFirst)
      {
        std::cout << ",";
      }
      aFirst = false;
      const int aIdx = anEdges.FindIndex(it.Value());
      if (aIdx == 0)
      {
        std::cout << "?" << static_cast<int>(it.Value().ShapeType());
      }
      else
      {
        std::cout << aIdx;
      }
    }
    std::cout << "]" << std::endl;
  }
  for (int i = 1; i <= anEdges.Extent(); ++i)
  {
    TopoDS_Vertex aVf, aVl;
    TopExp::Vertices(TopoDS::Edge(anEdges(i)), aVf, aVl);
    const int aIf = aVf.IsNull() ? -1 : aVerts.FindIndex(aVf);
    const int aIl = aVl.IsNull() ? -1 : aVerts.FindIndex(aVl);
    std::cout << "EDGE n=" << i << " first=" << aIf << " last=" << aIl << std::endl;
  }
  return 0;
}

// `--record` / `--record-rt`. The shape is the one `main` already read from
// `argv[1]`, so the read half is identical to every other mode of this probe.
int runRecord(const TopoDS_Shape& theShape, const char* theInput, const char* theOutPath,
              const double theDeflection, const double theAngle)
{
  return write_record(theShape, "read", theInput, theOutPath, theDeflection, theAngle);
}

int runRecordRoundTrip(const TopoDS_Shape& theShape, const char* theInput,
                       const char* theFormat, const char* theMidPath, const char* theOutPath,
                       const double theDeflection, const double theAngle)
{
  const std::string aFormat(theFormat);
  TopoDS_Shape      aReadBack;
  if (aFormat == "iges")
  {
    if (!write_iges_shape(theShape, theMidPath))
    {
      std::cerr << "iges write failed: " << theMidPath << "\n";
      return 7;
    }
    if (!read_iges_shape(theMidPath, aReadBack))
    {
      return 3;
    }
    return write_record(aReadBack, "rt-iges", theInput, theOutPath, theDeflection,
                        theAngle);
  }
  if (aFormat == "step")
  {
    if (!write_step_shape(theShape, theMidPath))
    {
      std::cerr << "step write failed: " << theMidPath << "\n";
      return 7;
    }
    if (!read_step_shape(theMidPath, aReadBack))
    {
      return 3;
    }
    return write_record(aReadBack, "rt-step", theInput, theOutPath, theDeflection,
                        theAngle);
  }
  std::cerr << "unknown round-trip format: " << theFormat << "\n";
  return 5;
}

namespace {

const char* status_name(const BRepCheck_Status theStatus)
{
  switch (theStatus)
  {
    case BRepCheck_NoError: return "NoError";
    case BRepCheck_InvalidPointOnCurve: return "InvalidPointOnCurve";
    case BRepCheck_InvalidPointOnCurveOnSurface: return "InvalidPointOnCurveOnSurface";
    case BRepCheck_InvalidPointOnSurface: return "InvalidPointOnSurface";
    case BRepCheck_No3DCurve: return "No3DCurve";
    case BRepCheck_Multiple3DCurve: return "Multiple3DCurve";
    case BRepCheck_Invalid3DCurve: return "Invalid3DCurve";
    case BRepCheck_NoCurveOnSurface: return "NoCurveOnSurface";
    case BRepCheck_InvalidCurveOnSurface: return "InvalidCurveOnSurface";
    case BRepCheck_InvalidCurveOnClosedSurface: return "InvalidCurveOnClosedSurface";
    case BRepCheck_InvalidSameRangeFlag: return "InvalidSameRangeFlag";
    case BRepCheck_InvalidSameParameterFlag: return "InvalidSameParameterFlag";
    case BRepCheck_InvalidDegeneratedFlag: return "InvalidDegeneratedFlag";
    case BRepCheck_FreeEdge: return "FreeEdge";
    case BRepCheck_InvalidMultiConnexity: return "InvalidMultiConnexity";
    case BRepCheck_InvalidRange: return "InvalidRange";
    case BRepCheck_EmptyWire: return "EmptyWire";
    case BRepCheck_RedundantEdge: return "RedundantEdge";
    case BRepCheck_SelfIntersectingWire: return "SelfIntersectingWire";
    case BRepCheck_NoSurface: return "NoSurface";
    case BRepCheck_InvalidWire: return "InvalidWire";
    case BRepCheck_RedundantWire: return "RedundantWire";
    case BRepCheck_IntersectingWires: return "IntersectingWires";
    case BRepCheck_InvalidImbricationOfWires: return "InvalidImbricationOfWires";
    case BRepCheck_EmptyShell: return "EmptyShell";
    case BRepCheck_RedundantFace: return "RedundantFace";
    case BRepCheck_InvalidImbricationOfShells: return "InvalidImbricationOfShells";
    case BRepCheck_UnorientableShape: return "UnorientableShape";
    case BRepCheck_NotClosed: return "NotClosed";
    case BRepCheck_NotConnected: return "NotConnected";
    case BRepCheck_SubshapeNotInShape: return "SubshapeNotInShape";
    case BRepCheck_BadOrientation: return "BadOrientation";
    case BRepCheck_BadOrientationOfSubshape: return "BadOrientationOfSubshape";
    case BRepCheck_InvalidPolygonOnTriangulation: return "InvalidPolygonOnTriangulation";
    case BRepCheck_InvalidToleranceValue: return "InvalidToleranceValue";
    case BRepCheck_EnclosedRegion: return "EnclosedRegion";
    case BRepCheck_CheckFail: return "CheckFail";
  }
  return "Unknown";
}

void print_status_list(const NCollection_List<BRepCheck_Status>& theList)
{
  bool aFirst = true;
  for (NCollection_List<BRepCheck_Status>::Iterator it(theList); it.More(); it.Next())
  {
    if (!aFirst)
    {
      std::cout << ",";
    }
    aFirst = false;
    std::cout << status_name(it.Value());
  }
  if (aFirst)
  {
    std::cout << "-";
  }
}

} // namespace

// `--statusdump`: `BRepCheck_Analyzer` verdict per sub-shape, the OCCT half of
// the port's `record_dump --statuses`. `BRepCheck_Result` keeps one status list
// for the shape itself (`Status()`) and one per context shape
// (`StatusOnShape(ctx)`); the face-context list is where
// `BRepCheck_Edge::InContext` / `BRepCheck_Wire::InContext` record
// `InvalidCurveOnSurface` / `SelfIntersectingWire`, so both are printed:
//   VALID <0|1>
//   STATUS type=<..> n=<k> min=<list> [ctx=<type>:<list> ...]
int runStatusDump(const TopoDS_Shape& theShape)
{
  BRepCheck_Analyzer aCheck(theShape);
  std::cout << "VALID " << (aCheck.IsValid() ? 1 : 0) << std::endl;
  TopTools_IndexedMapOfShape aMap;
  TopExp::MapShapes(theShape, aMap);
  for (int i = 1; i <= aMap.Extent(); ++i)
  {
    const TopoDS_Shape&               aSub    = aMap(i);
    const occ::handle<BRepCheck_Result>& aResult = aCheck.Result(aSub);
    if (aResult.IsNull())
    {
      continue;
    }
    std::cout << "STATUS type=" << shape_type_name(aSub.ShapeType()) << " n=" << i
              << " min=";
    print_status_list(aResult->Status());
    aResult->InitContextIterator();
    for (; aResult->MoreShapeInContext(); aResult->NextShapeInContext())
    {
      std::cout << " ctx=" << shape_type_name(aResult->ContextualShape().ShapeType()) << ":";
      print_status_list(aResult->StatusOnShape());
    }
    std::cout << std::endl;
  }
  return 0;
}
