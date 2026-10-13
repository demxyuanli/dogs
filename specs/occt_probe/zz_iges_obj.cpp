// Minimal "standard OCCT tool" IGES -> mesh, used to answer one question:
// when a plain OCCT chain (IGESControl_Controller::Init + IGESControl_Reader +
// TransferRoots + OneShape + BRepMesh_IncrementalMesh) reports a face that
// meshes only a fragment, does any of the usual post-processing steps
// (ShapeFix_Shape, ShapeUpgrade_UnifySameDomain, BRepBuilderAPI_Sewing)
// change the face/wire domain so the whole surface gets meshed?
//
// Each mode re-reads the file so the post-processing steps cannot contaminate
// each other. Output is one line per mode: face count, mesh nodes, mesh
// triangles and the mesh bbox, i.e. exactly the quantities that separate a
// "wrong domain" result (bbox too small) from a "too coarse" one (bbox right,
// counts low).
#include <IGESControl_Controller.hxx>
#include <IGESControl_Reader.hxx>
#include <IGESControl_Writer.hxx>
#include <STEPControl_Reader.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <Geom_Surface.hxx>
#include <Standard_Type.hxx>
#include <IGESGeom_TrimmedSurface.hxx>
#include <IGESGeom_CurveOnSurface.hxx>
#include <Interface_InterfaceModel.hxx>
#include <IGESData_IGESEntity.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Curve2d.hxx>
#include <Geom2d_Curve.hxx>
#include <TopoDS_Shape.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRep_Tool.hxx>
#include <Poly_Triangulation.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopLoc_Location.hxx>
#include <ShapeFix_Shape.hxx>
#include <ShapeProcess.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>

static int countFaces(const TopoDS_Shape& theShape)
{
  int aN = 0;
  for (TopExp_Explorer anEx(theShape, TopAbs_FACE); anEx.More(); anEx.Next())
  {
    ++aN;
  }
  return aN;
}

static int countWires(const TopoDS_Shape& theShape)
{
  int aN = 0;
  for (TopExp_Explorer anEx(theShape, TopAbs_WIRE); anEx.More(); anEx.Next())
  {
    ++aN;
  }
  return aN;
}

static int countEdges(const TopoDS_Shape& theShape)
{
  int aN = 0;
  for (TopExp_Explorer anEx(theShape, TopAbs_EDGE); anEx.More(); anEx.Next())
  {
    ++aN;
  }
  return aN;
}

// Print what the reader actually produced: surface kind, periodicity, the
// wire-driven UV domain and each boundary loop. This is what separates an
// analytic single-loop face (whole torus) from a derived-surface face cut by
// several loops (fragment).
static void dumpTopology(const TopoDS_Shape& theShape)
{
  for (TopExp_Explorer anEx(theShape, TopAbs_FACE); anEx.More(); anEx.Next())
  {
    const TopoDS_Face aFace = TopoDS::Face(anEx.Current());
    double aU1 = 0., aU2 = 0., aV1 = 0., aV2 = 0.;
    BRepTools::UVBounds(aFace, aU1, aU2, aV1, aV2);
    const occ::handle<Geom_Surface>& aSurf = BRep_Tool::Surface(aFace);
    std::cout << "face surface=" << aSurf->DynamicType()->Name()
              << " uperiodic=" << aSurf->IsUPeriodic() << "(u=" << aSurf->UPeriod()
              << ") vperiodic=" << aSurf->IsVPeriodic() << "(v=" << aSurf->VPeriod() << ")"
              << std::fixed << std::setprecision(7) << " uv=[" << aU1 << "," << aU2 << "]x[" << aV1 << ","
              << aV2 << "]\n";
    int aW = 0;
    for (TopExp_Explorer aWx(aFace, TopAbs_WIRE); aWx.More(); aWx.Next(), ++aW)
    {
      int    aE   = 0;
      double aLen = 0.;
      for (TopExp_Explorer aEx(aWx.Current(), TopAbs_EDGE); aEx.More(); aEx.Next())
      {
        ++aE;
        double aF = 0., aL = 0.;
        BRep_Tool::Range(TopoDS::Edge(aEx.Current()), aF, aL);
        aLen += aL - aF;
      }
      std::cout << "  wire " << aW << " edges=" << aE << std::fixed << std::setprecision(4)
                << " total_range=" << aLen << "\n";
      int aI = 0;
      for (TopExp_Explorer aEx(aWx.Current(), TopAbs_EDGE); aEx.More(); aEx.Next(), ++aI)
      {
        const TopoDS_Edge& anEdge = TopoDS::Edge(aEx.Current());
        double             aF = 0., aL = 0.;
        BRepAdaptor_Curve  aC(anEdge);
        const gp_Pnt       aP0 = aC.Value(aC.FirstParameter());
        const gp_Pnt       aP1 = aC.Value(aC.LastParameter());
        std::cout << "    edge " << aI << " 3dtype=" << aC.GetType() << std::fixed
                  << std::setprecision(4) << " p0=(" << aP0.X() << "," << aP0.Y() << "," << aP0.Z()
                  << ") p1=(" << aP1.X() << "," << aP1.Y() << "," << aP1.Z() << ")";
        if (BRep_Tool::IsClosed(anEdge))
        {
          std::cout << " closed";
        }
        // Distinguish which representation built this edge: a 3D curve alone, a
        // pcurve alone, or a 3D curve carrying a transferred pcurve.
        TopLoc_Location     aL3d;
        double              aT1 = 0., aT2 = 0.;
        const occ::handle<Geom_Curve>& aC3d = BRep_Tool::Curve(anEdge, aL3d, aT1, aT2);
        std::cout << " has3d=" << (!aC3d.IsNull());
        occ::handle<Geom2d_Curve> aPc;
        double                    aPF = 0., aPL = 0.;
        aPc = BRep_Tool::CurveOnSurface(anEdge, aFace, aPF, aPL);
        std::cout << " haspcurve=" << (!aPc.IsNull());
        if (!aPc.IsNull())
        {
          gp_Pnt2d aU0, aU1;
          aPc->D0(aPF, aU0);
          aPc->D0(aPL, aU1);
          std::cout << " pcurve_uv=(" << aU0.X() << "," << aU0.Y() << ")-(" << aU1.X() << ","
                    << aU1.Y() << ")";
        }
        std::cout << "\n";
      }
    }
  }
}

// Ask the IGES model itself what it declares: how many boundary loops the
// trimmed surface carries and, for the curve-on-surface, which of the three
// representations (CRV3D / parameter-space / model-space) are present. This
// separates "the file asks for two loops" from "the reader invented a second".
static void dumpIgesModel(const char* theFile)
{
  IGESControl_Controller::Init();
  IGESControl_Reader aReader;
  if (aReader.ReadFile(theFile) != IFSelect_RetDone)
  {
    return;
  }
  const occ::handle<Interface_InterfaceModel>& aModel = aReader.Model();
  for (int i = 1; i <= aModel->NbEntities(); ++i)
  {
    const occ::handle<Standard_Transient>& anEnt = aModel->Value(i);
    if (anEnt.IsNull())
    {
      continue;
    }
    if (anEnt->IsKind(STANDARD_TYPE(IGESGeom_TrimmedSurface)))
    {
      const occ::handle<IGESGeom_TrimmedSurface>& aTS =
        occ::down_cast<IGESGeom_TrimmedSurface>(anEnt);
      std::cout << "IGES 144 #" << i << " outer=" << aTS->HasOuterContour()
                << " nInner=" << aTS->NbInnerContours();
      if (aTS->HasOuterContour())
      {
        std::cout << " outerType=" << aTS->OuterContour()->DynamicType()->Name();
      }
      std::cout << "\n";
    }
    else if (anEnt->IsKind(STANDARD_TYPE(IGESGeom_CurveOnSurface)))
    {
      const occ::handle<IGESGeom_CurveOnSurface>& aCS = occ::down_cast<IGESGeom_CurveOnSurface>(anEnt);
      std::cout << "IGES 142 #" << i << " pref=" << aCS->PreferenceMode()
                << " crv3d=" << (aCS->Curve3D().IsNull() ? "null" : "set")
                << " uv=" << (aCS->CurveUV().IsNull() ? "null" : "set") << "\n";
    }
  }
}

static void meshStats(const TopoDS_Shape& theShape,
                      const double      theDeflection,
                      const double      theAngle,
                      const char*       theTag)
{
  BRepMesh_IncrementalMesh aMesher(theShape, theDeflection, false, theAngle);
  int    aNodes = 0, aTris = 0, aFaces = 0;
  double aMin[3] = {1e300, 1e300, 1e300};
  double aMax[3] = {-1e300, -1e300, -1e300};
  for (TopExp_Explorer anEx(theShape, TopAbs_FACE); anEx.More(); anEx.Next())
  {
    ++aFaces;
    TopLoc_Location aLoc;
    const occ::handle<Poly_Triangulation>& aTri =
      BRep_Tool::Triangulation(TopoDS::Face(anEx.Current()), aLoc);
    if (aTri.IsNull())
    {
      continue;
    }
    const gp_Trsf& aT = aLoc.Transformation();
    for (int i = 1; i <= aTri->NbNodes(); ++i)
    {
      const gp_Pnt p = aTri->Node(i).Transformed(aT);
      aMin[0]        = std::min(aMin[0], p.X());
      aMin[1]        = std::min(aMin[1], p.Y());
      aMin[2]        = std::min(aMin[2], p.Z());
      aMax[0]        = std::max(aMax[0], p.X());
      aMax[1]        = std::max(aMax[1], p.Y());
      aMax[2]        = std::max(aMax[2], p.Z());
    }
    aNodes += aTri->NbNodes();
    aTris += aTri->NbTriangles();
  }
  std::cout << theTag << " faces=" << aFaces << " wires=" << countWires(theShape)
            << " nodes=" << aNodes << " triangles=" << aTris;
  if (aNodes > 0)
  {
    std::cout << std::fixed << std::setprecision(4) << " bbox=(" << aMin[0] << "," << aMin[1] << ","
              << aMin[2] << ")-(" << aMax[0] << "," << aMax[1] << "," << aMax[2] << ")";
  }
  else
  {
    std::cout << " bbox=(empty)";
  }
  std::cout << "\n";
}

static TopoDS_Shape readIges(const char* theFile, const bool theNoFix)
{
  // `.brep` inputs let a foreign reader (FreeCAD) hand over the shape it built,
  // so the mesher can be compared on identical topology.
  const std::string aName(theFile);
  if (aName.size() >= 5 && aName.compare(aName.size() - 5, 5, ".brep") == 0)
  {
    TopoDS_Shape aShape;
    BRep_Builder aBuilder;
    BRepTools::Read(aShape, theFile, aBuilder);
    if (aShape.IsNull())
    {
      std::cerr << "brep read failed: " << theFile << "\n";
    }
    return aShape;
  }
  // STEP is the other oracle for the same solids: an analytic surface with a
  // single seam-bounded loop, i.e. the shape the IGES is supposed to describe.
  const bool isStep = aName.size() >= 5
                      && (aName.compare(aName.size() - 5, 5, ".step") == 0
                          || aName.compare(aName.size() - 4, 4, ".stp") == 0);
  if (isStep)
  {
    STEPControl_Reader aStepReader;
    if (aStepReader.ReadFile(theFile) != IFSelect_RetDone)
    {
      std::cerr << "step read failed: " << theFile << "\n";
      return TopoDS_Shape();
    }
    aStepReader.TransferRoots();
    return aStepReader.OneShape();
  }
  IGESControl_Controller::Init();
  IGESControl_Reader aReader;
  if (theNoFix)
  {
    ShapeProcess::OperationsFlags aFlags;
    aFlags.reset();
    aReader.SetShapeProcessFlags(aFlags);
  }
  if (aReader.ReadFile(theFile) != IFSelect_RetDone)
  {
    std::cerr << "read failed: " << theFile << "\n";
    return TopoDS_Shape();
  }
  aReader.TransferRoots();
  return aReader.OneShape();
}

int main(int argc, char** argv)
{
  if (argc < 2)
  {
    std::cerr << "usage: zz_iges_obj <file.iges|file.brep|file.step> [defl] [angle] [nofix]\n";
    std::cerr << "       zz_iges_obj data/x.step --iges <out.iges> [brepmode 0|1]\n";
    return 2;
  }
  if (std::string(argv[1]) == "--iges")
  {
    // Faces mode (mode 0, what the port's oracle uses) maps an analytic torus
    // to a 120 SurfaceOfRevolution; BRep mode (mode 1) reaches
    // GeomToIGES_GeomSurface::TransferToroidalSurface and writes a 198.
    // usage: zz_iges_obj --iges <in.step> <out.iges> <mode>
    const char*        aIn   = argv[2];
    const char*        aOut  = argv[3];
    const int          aMode = (argc > 4) ? std::atoi(argv[4]) : 0;
    const TopoDS_Shape aShape = readIges(aIn, false);
    IGESControl_Controller::Init();
    IGESControl_Writer aWriter("MM", aMode);
    aWriter.AddShape(aShape);
    aWriter.ComputeModel();
    std::cout << "iges-write in=" << aIn << " out=" << aOut << " mode=" << aMode
              << " ok=" << (int)aWriter.Write(aOut) << "\n";
    return 0;
  }
  const char*  aFile       = argv[1];
  const double aDeflection = (argc > 2) ? std::atof(argv[2]) : 0.096;
  const double anAngle     = (argc > 3) ? std::atof(argv[3]) : 0.3490658503988659;
  const bool   aNoFix      = (argc > 4) && std::string(argv[4]) == "nofix";

  {
    const TopoDS_Shape aShape = readIges(aFile, aNoFix);
    std::cout << "read faces=" << countFaces(aShape) << " wires=" << countWires(aShape)
              << " edges=" << countEdges(aShape) << " null=" << aShape.IsNull()
              << " nofix=" << aNoFix << "\n";
    dumpIgesModel(aFile);
    dumpTopology(aShape);
    meshStats(aShape, aDeflection, anAngle, "plain           ");
  }
  {
    TopoDS_Shape    aShape = readIges(aFile, aNoFix);
    ShapeFix_Shape aFix(aShape);
    aFix.Perform();
    meshStats(aFix.Shape(), aDeflection, anAngle, "+fixshape       ");
  }
  {
    TopoDS_Shape                     aShape = readIges(aFile, aNoFix);
    ShapeUpgrade_UnifySameDomain aUnify(aShape, false, true, false);
    aUnify.Build();
    meshStats(aUnify.Shape(), aDeflection, anAngle, "+unifysamedomain");
  }
  {
    TopoDS_Shape           aShape = readIges(aFile, aNoFix);
    BRepBuilderAPI_Sewing aSewing;
    aSewing.Add(aShape);
    aSewing.Perform();
    meshStats(aSewing.SewedShape(), aDeflection, anAngle, "+sewing         ");
  }
  {
    TopoDS_Shape      aShape = readIges(aFile, aNoFix);
    ShapeFix_Shape aFix(aShape);
    aFix.Perform();
    ShapeUpgrade_UnifySameDomain aUnify(aFix.Shape(), false, true, false);
    aUnify.Build();
    meshStats(aUnify.Shape(), aDeflection, anAngle, "+fix+unify      ");
  }
  return 0;
}
