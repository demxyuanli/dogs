// OCCT ground-truth generator for the polygonal HLR parity gate.
//
// Reads a STEP model, meshes it with BRepMesh_IncrementalMesh, runs
// HLRBRep_PolyAlgo -> HLRBRep_PolyHLRToShape over the same projection the Rust
// `hlr_poly_probe` example uses, and dumps, in the same text format:
//   - the mesh census (faces / triangles / committed face triangulations),
//   - the visible / hidden bi-point segments of `HLRBRep_PolyHLRToShape::Update`
//     (HLRBRep_PolyHLRToShape.cxx:37-79),
//   - the edge count of each of the eight compounds (lxx:35-143).
//
// Usage: hlr_probe <step-file> <deflection> <angular-rad> <view: xoy|xoz|yoz>
//
// The coordinates printed are the doubles OCCT computes; the Rust probe prints
// them with the same "%.9f" formatting so the two transcripts diff.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRep_Tool.hxx>
#include <HLRBRep_Curve.hxx>
#include <HLRAlgo_BiPoint.hxx>
#include <HLRAlgo_EdgeIterator.hxx>
#include <HLRAlgo_EdgeStatus.hxx>
#include <HLRAlgo_PolyAlgo.hxx>
#include <HLRAlgo_PolyData.hxx>
#include <HLRAlgo_PolyShellData.hxx>
#include <HLRAlgo_TriangleData.hxx>
#include <HLRBRep_PolyAlgo.hxx>
#include <HLRBRep_PolyHLRToShape.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_List.hxx>
#include <Poly_Triangulation.hxx>
#include <STEPControl_Reader.hxx>
#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <gp_Ax2.hxx>
#include <gp_Dir.hxx>
#include <gp_Dir2d.hxx>
#include <gp_Lin2d.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Trsf.hxx>
#include <gp_Vec2d.hxx>
#include <gp_XY.hxx>
#include <gp_XYZ.hxx>

static const char* boolText(Standard_Boolean theFlag)
{
  return theFlag ? "true" : "false";
}

static int countEdges(const TopoDS_Shape& theShape)
{
  if (theShape.IsNull())
  {
    return -1;
  }
  int aNb = 0;
  for (TopExp_Explorer anExp(theShape, TopAbs_EDGE); anExp.More(); anExp.Next())
  {
    ++aNb;
  }
  return aNb;
}

static const char* baseName(const char* thePath)
{
  const char* aSlash = std::strrchr(thePath, '/');
  const char* aBack  = std::strrchr(thePath, '\\');
  const char* aLast  = (aBack != nullptr && (aSlash == nullptr || aBack > aSlash)) ? aBack : aSlash;
  return (aLast == nullptr) ? thePath : aLast + 1;
}

static void printCompound(const char* theName, const TopoDS_Shape& theShape)
{
  if (theShape.IsNull())
  {
    std::printf("%s: null\n", theName);
    return;
  }
  std::printf("%s: edges=%d\n", theName, countEdges(theShape));
}

int main(int argc, char** argv)
{
  if (argc < 2)
  {
    std::fprintf(stderr, "usage: hlr_probe <step> [deflection] [angle-rad] [xoy|xoz|yoz]\n");
    return 2;
  }
  const char*       aPath  = argv[1];
  const double      aLin   = (argc > 2) ? std::atof(argv[2]) : 0.1;
  const double      aAngle = (argc > 3) ? std::atof(argv[3]) : 0.3490658503988659;
  const std::string aView  = (argc > 4) ? argv[4] : "xoy";
  const std::string aMode  = (argc > 5) ? argv[5] : "";

  STEPControl_Reader aReader;
  if (aReader.ReadFile(aPath) != IFSelect_RetDone)
  {
    std::fprintf(stderr, "read failed: %s\n", aPath);
    return 1;
  }
  aReader.TransferRoots();
  TopoDS_Shape aShape = aReader.OneShape();
  if (aShape.IsNull())
  {
    std::fprintf(stderr, "no shape\n");
    return 1;
  }

  BRepMesh_IncrementalMesh aMesher(aShape, aLin, Standard_False, aAngle, Standard_True);
  aMesher.Perform();

  int aNbFaces = 0, aNbTris = 0, aCommitted = 0;
  for (TopExp_Explorer anExp(aShape, TopAbs_FACE); anExp.More(); anExp.Next())
  {
    ++aNbFaces;
    TopLoc_Location            aLoc;
    Handle(Poly_Triangulation) aTri =
      BRep_Tool::Triangulation(TopoDS::Face(anExp.Current()), aLoc);
    if (!aTri.IsNull())
    {
      ++aCommitted;
      aNbTris += aTri->NbTriangles();
    }
  }
  std::printf("PROBE %s lin=%g view=%s\n", baseName(aPath), aLin, aView.c_str());
  std::printf("mesh faces=%d triangles=%d\n", aNbFaces, aNbTris);
  std::printf("committed_faces=%d uncommitted_faces=%d\n", aCommitted, aNbFaces - aCommitted);

  gp_Dir aDir(0.0, 0.0, 1.0);
  gp_Dir aVx(1.0, 0.0, 0.0);
  if (aView == "xoz")
  {
    aDir = gp_Dir(0.0, -1.0, 0.0);
  }
  else if (aView == "yoz")
  {
    aDir = gp_Dir(-1.0, 0.0, 0.0);
    aVx  = gp_Dir(0.0, 1.0, 0.0);
  }
  gp_Ax2            aPlane(gp_Pnt(0.0, 0.0, 0.0), aDir, aVx);
  HLRAlgo_Projector aProjector(aPlane);

  // Exact-HLR curve layer ground truth: the `HLRBRep_Curve` state after
  // `Update` plus the 2D / 3D evaluation at five parameters, for every edge of
  // `TopExp::MapShapes(shape, TopAbs_EDGE)`. Mirrored by the Rust example
  // `hlrbrep_curve_probe`.
  if (aMode == "curve")
  {
    TopTools_IndexedMapOfShape anEMap;
    TopExp::MapShapes(aShape, TopAbs_EDGE, anEMap);
    std::printf("emap=%d\n", anEMap.Extent());
    for (int ei = 1; ei <= anEMap.Extent(); ++ei)
    {
      const TopoDS_Edge& anEdge = TopoDS::Edge(anEMap.FindKey(ei));
      HLRBRep_Curve aCurve;
      aCurve.Projector(&aProjector);
      aCurve.Curve(anEdge);
      double aMin[16], aMax[16];
      for (int k = 0; k < 16; ++k)
      {
        aMin[k] = RealLast();
        aMax[k] = -RealLast();
      }
      const double aTol = aCurve.Update(aMin, aMax);
      std::printf("EDGE %d type=%d tol=%.9f first=%.9f last=%.9f\n",
                  ei,
                  (int)aCurve.GetType(),
                  aTol,
                  aCurve.FirstParameter(),
                  aCurve.LastParameter());
      std::printf(" MIN");
      for (int k = 0; k < 16; ++k)
      {
        std::printf(" %.9f", aMin[k]);
      }
      std::printf("\n MAX");
      for (int k = 0; k < 16; ++k)
      {
        std::printf(" %.9f", aMax[k]);
      }
      std::printf("\n");

      BRepAdaptor_Curve& a3d = aCurve.Curve();
      const double       u0  = a3d.FirstParameter();
      const double       u1  = a3d.LastParameter();
      for (int s = 0; s <= 4; ++s)
      {
        const double     u   = u0 + (u1 - u0) * s / 4.0;
        const gp_Pnt     p3  = aCurve.Value3D(u);
        gp_Pnt2d         pd0;
        aCurve.D0(u, pd0);
        gp_Pnt2d pd1;
        gp_Vec2d vd1;
        aCurve.D1(u, pd1, vd1);
        gp_Pnt2d pd2;
        gp_Vec2d vd21, vd22;
        aCurve.D2(u, pd2, vd21, vd22);
        const double az    = aCurve.Z(u);
        const double apar2 = aCurve.Parameter2d(u);
        const double apar3 = aCurve.Parameter3d(apar2);
        std::printf(
          " S%d u=%.9f p3=(%.9f,%.9f,%.9f) d0=(%.9f,%.9f) d1p=(%.9f,%.9f) d1v=(%.9f,%.9f)"
          " d2p=(%.9f,%.9f) d2v1=(%.9f,%.9f) d2v2=(%.9f,%.9f) z=%.9f par2=%.9f par3=%.9f\n",
          s,
          u,
          p3.X(),
          p3.Y(),
          p3.Z(),
          pd0.X(),
          pd0.Y(),
          pd1.X(),
          pd1.Y(),
          vd1.X(),
          vd1.Y(),
          pd2.X(),
          pd2.Y(),
          vd21.X(),
          vd21.Y(),
          vd22.X(),
          vd22.Y(),
          az,
          apar2,
          apar3);
      }
      for (int at = 0; at <= 1; ++at)
      {
        gp_Pnt2d TP;
        gp_Dir2d TD;
        try
        {
          aCurve.Tangent(at == 0, TP, TD);
          std::printf(" TAN%d p=(%.9f,%.9f) d=(%.9f,%.9f)\n", at + 1, TP.X(), TP.Y(), TD.X(), TD.Y());
        }
        catch (...)
        {
          std::printf(" TAN%d undef\n", at + 1);
        }
      }
      if (aCurve.GetType() == GeomAbs_Line)
      {
        const gp_Lin2d aL = aCurve.Line();
        std::printf(" LINE p=(%.9f,%.9f) d=(%.9f,%.9f)\n",
                    aL.Location().X(),
                    aL.Location().Y(),
                    aL.Direction().X(),
                    aL.Direction().Y());
      }
    }
    std::printf("PROBE END\n");
    return 0;
  }

  occ::handle<HLRBRep_PolyAlgo> aAlgo = new HLRBRep_PolyAlgo(aShape);
  aAlgo->Projector(aProjector);
  aAlgo->Update();
  std::printf("poly_algo more_hide=%s\n", boolText(aAlgo->MoreHide()));

  // Intermediate dump of every shell / face / triangle / bi-point, so the port
  // can be diffed stage by stage (PolyShellData::Edges/PolyData are public).
  {
    const HLRAlgo_PolyAlgo&                          aPA     = *aAlgo->Algo();
    const NCollection_Array1<occ::handle<HLRAlgo_PolyShellData>>& aShells = aPA.PolyShell();
    std::printf("shells=%d\n", aShells.Length());
    for (int s = aShells.Lower(); s <= aShells.Upper(); ++s)
    {
      const occ::handle<HLRAlgo_PolyShellData>& aShell = aShells.Value(s);
      if (aShell.IsNull())
      {
        std::printf("SHELL %d null\n", s);
        continue;
      }
      std::printf("SHELL %d hiding=%d npoly=%d nhidpoly=%d\n",
                  s,
                  aShell->Hiding() ? 1 : 0,
                  aShell->PolyData().Length(),
                  aShell->HidingPolyData().Length());
      for (int f = aShell->PolyData().Lower(); f <= aShell->PolyData().Upper(); ++f)
      {
        const occ::handle<HLRAlgo_PolyData>& aPD = aShell->PolyData().Value(f);
        if (aPD.IsNull())
        {
          std::printf(" FACE %d null\n", f);
          continue;
        }
        const NCollection_Array1<gp_XYZ>&               aNodes = aPD->Nodes();
        const NCollection_Array1<HLRAlgo_TriangleData>& aTD    = aPD->TData();
        std::printf(" FACE %d nodes=%d tris=%d hid=%d\n",
                    aPD->FaceIndex(),
                    aNodes.Length(),
                    aTD.Length(),
                    aPD->Hiding() ? 1 : 0);
        for (int n = aNodes.Lower(); n <= aNodes.Upper(); ++n)
        {
          const gp_XYZ& aP = aNodes.Value(n);
          std::printf("  NODE %d (%.9f,%.9f,%.9f)\n", n, aP.X(), aP.Y(), aP.Z());
        }
        for (int t = aTD.Lower(); t <= aTD.Upper(); ++t)
        {
          const HLRAlgo_TriangleData& aTri = aTD.Value(t);
          std::printf("  TRI %d n=(%d,%d,%d) flags=0x%x\n",
                      t,
                      aTri.Node1,
                      aTri.Node2,
                      aTri.Node3,
                      static_cast<unsigned>(aTri.Flags));
        }
        if (aPD->Hiding())
        {
          const NCollection_Array1<HLRAlgo_PolyHidingData>& aPH = aPD->PHDat();
          for (int h = aPH.Lower(); h <= aPH.Upper(); ++h)
          {
            HLRAlgo_PolyHidingData       aItem = aPH.Value(h);
            const HLRAlgo_PolyHidingData::PlaneT& aPl  = aItem.Plane();
            std::printf("  HIDE %d tri=%d min=%d max=%d plane=(%.9f,%.9f,%.9f,%.9f)\n",
                        h,
                        aItem.Indices().Index,
                        aItem.Indices().Min,
                        aItem.Indices().Max,
                        aPl.Normal.X(),
                        aPl.Normal.Y(),
                        aPl.Normal.Z(),
                        aPl.D);
          }
        }
      }
      int aNbBp = 0;
      for (NCollection_List<HLRAlgo_BiPoint>::Iterator anIt(aShell->Edges()); anIt.More();
           anIt.Next())
      {
        const HLRAlgo_BiPoint& aBP = anIt.Value();
        const HLRAlgo_BiPoint::IndicesT aInd = const_cast<HLRAlgo_BiPoint&>(aBP).Indices();
        const HLRAlgo_BiPoint::PointsT  aPts = const_cast<HLRAlgo_BiPoint&>(aBP).Points();
        std::printf(
          " BP %d p1=(%.9f,%.9f,%.9f) p2=(%.9f,%.9f,%.9f) t1=(%.9f,%.9f,%.9f) t2=(%.9f,%.9f,%.9f)"
          " idx=%d c1=%d f1=(%d,%d) c2=%d f2=(%d,%d) minseg=%d maxseg=%d hid=%d rg1=%d rgn=%d"
          " outl=%d intl=%d\n",
          aNbBp,
          aPts.Pnt1.X(),
          aPts.Pnt1.Y(),
          aPts.Pnt1.Z(),
          aPts.Pnt2.X(),
          aPts.Pnt2.Y(),
          aPts.Pnt2.Z(),
          aPts.PntP1.X(),
          aPts.PntP1.Y(),
          aPts.PntP1.Z(),
          aPts.PntP2.X(),
          aPts.PntP2.Y(),
          aPts.PntP2.Z(),
          aInd.ShapeIndex,
          aInd.FaceConex1,
          aInd.Face1Pt1,
          aInd.Face1Pt2,
          aInd.FaceConex2,
          aInd.Face2Pt1,
          aInd.Face2Pt2,
          aInd.MinSeg,
          aInd.MaxSeg,
          aBP.Hidden() ? 1 : 0,
          aBP.Rg1Line() ? 1 : 0,
          aBP.RgNLine() ? 1 : 0,
          aBP.OutLine() ? 1 : 0,
          aBP.IntLine() ? 1 : 0);
        ++aNbBp;
      }
    }
  }

  // Replay of HLRBRep_PolyHLRToShape::Update (cxx:37-79) so every segment and
  // its flags can be printed.
  const gp_Trsf&     T = aProjector.Transformation();
  HLRAlgo_EdgeStatus status;
  HLRAlgo_EdgeIterator anIt;
  double               sta = 0.0, end = 0.0;
  float                tolsta = 0.0f, tolend = 0.0f;
  std::string          aVis, aHid;
  int                  aNbVis = 0, aNbHid = 0;
  char                 aLine[256];

  for (aAlgo->InitHide(); aAlgo->MoreHide(); aAlgo->NextHide())
  {
    TopoDS_Shape aSegShape;
    Standard_Boolean reg1 = Standard_False, regn = Standard_False, outl = Standard_False,
                     intl = Standard_False;
    const HLRAlgo_BiPoint::PointsT& aPoints =
      aAlgo->Hide(status, aSegShape, reg1, regn, outl, intl);
    gp_XYZ aSta = aPoints.Pnt1;
    gp_XYZ aEnd = aPoints.Pnt2;
    T.Transforms(aSta);
    T.Transforms(aEnd);
    const gp_XY aSta2D(aSta.X(), aSta.Y());
    const gp_XY aEnd2D(aEnd.X(), aEnd.Y());
    const gp_XY aD = aEnd2D - aSta2D;
    if (aD.Modulus() <= 1.0e-10)
    {
      continue;
    }
    for (anIt.InitVisible(status); anIt.MoreVisible(); anIt.NextVisible())
    {
      anIt.Visible(sta, tolsta, end, tolend);
      const gp_XY aP1 = aSta2D + sta * aD;
      const gp_XY aP2 = aSta2D + end * aD;
      std::snprintf(aLine,
                    sizeof(aLine),
                    "VIS (%.9f,%.9f)-(%.9f,%.9f) rg1=%s rgn=%s outl=%s intl=%s\n",
                    aP1.X(),
                    aP1.Y(),
                    aP2.X(),
                    aP2.Y(),
                    boolText(reg1),
                    boolText(regn),
                    boolText(outl),
                    boolText(intl));
      aVis += aLine;
      ++aNbVis;
    }
    for (anIt.InitHidden(status); anIt.MoreHidden(); anIt.NextHidden())
    {
      anIt.Hidden(sta, tolsta, end, tolend);
      const gp_XY aP1 = aSta2D + sta * aD;
      const gp_XY aP2 = aSta2D + end * aD;
      std::snprintf(aLine,
                    sizeof(aLine),
                    "HID (%.9f,%.9f)-(%.9f,%.9f) rg1=%s rgn=%s outl=%s intl=%s\n",
                    aP1.X(),
                    aP1.Y(),
                    aP2.X(),
                    aP2.Y(),
                    boolText(reg1),
                    boolText(regn),
                    boolText(outl),
                    boolText(intl));
      aHid += aLine;
      ++aNbHid;
    }
  }

  HLRBRep_PolyHLRToShape aToShape;
  aToShape.Update(aAlgo);

  // Edge polygon-on-triangulation dump: what `BRepMesh_ModelPostProcessor`'s
  // PolygonCommitter committed, read back exactly as `InitBiPointsWithConnexity`
  // reads it (HLRBRep_PolyAlgo.cxx:900-910).
  {
    TopTools_IndexedMapOfShape anEMap, aFMap;
    TopExp::MapShapes(aShape, TopAbs_EDGE, anEMap);
    TopExp::MapShapes(aShape, TopAbs_FACE, aFMap);
    std::printf("emap=%d fmap=%d\n", anEMap.Extent(), aFMap.Extent());
    for (int fi = 1; fi <= aFMap.Extent(); ++fi)
    {
      const TopoDS_Face& aFace = TopoDS::Face(aFMap.FindKey(fi));
      TopLoc_Location    aLoc;
      Handle(Poly_Triangulation) aTr = BRep_Tool::Triangulation(aFace, aLoc);
      for (TopExp_Explorer anExp(aFace, TopAbs_EDGE); anExp.More(); anExp.Next())
      {
        const TopoDS_Edge& aEdge = TopoDS::Edge(anExp.Current());
        const int          ei    = anEMap.FindIndex(aEdge);
        Handle(Poly_PolygonOnTriangulation) aPol =
          BRep_Tool::PolygonOnTriangulation(aEdge, aTr, aLoc);
        if (aPol.IsNull())
        {
          std::printf("POLY face=%d edge=%d null\n", fi, ei);
          continue;
        }
        std::printf("POLY face=%d edge=%d n=%d nodes=", fi, ei, aPol->NbNodes());
        for (int n = aPol->Nodes().Lower(); n <= aPol->Nodes().Upper(); ++n)
        {
          std::printf("%s%d", (n == aPol->Nodes().Lower()) ? "" : ",", aPol->Nodes().Value(n));
        }
        std::printf(" params=");
        if (!aPol->Parameters().IsNull())
        {
          for (int n = aPol->Parameters()->Lower(); n <= aPol->Parameters()->Upper(); ++n)
          {
            std::printf("%s%.9f", (n == aPol->Parameters()->Lower()) ? "" : ",", aPol->Parameters()->Value(n));
          }
        }
        std::printf("\n");
      }
    }
  }

  std::printf("bi_pnt_vis=%d bi_pnt_hid=%d\n", aNbVis, aNbHid);
  printCompound("VCompound", aToShape.VCompound());
  printCompound("Rg1LineVCompound", aToShape.Rg1LineVCompound());
  printCompound("RgNLineVCompound", aToShape.RgNLineVCompound());
  printCompound("OutLineVCompound", aToShape.OutLineVCompound());
  printCompound("HCompound", aToShape.HCompound());
  printCompound("Rg1LineHCompound", aToShape.Rg1LineHCompound());
  printCompound("RgNLineHCompound", aToShape.RgNLineHCompound());
  printCompound("OutLineHCompound", aToShape.OutLineHCompound());
  std::printf("%s", aVis.c_str());
  std::printf("%s", aHid.c_str());
  std::printf("PROBE END\n");
  return 0;
}
