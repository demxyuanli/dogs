// Scratch ground-truth probe: box with a Z through hole cut by a perpendicular
// Y cylinder. Dumps the OCCT result (volume, per-face area/wires/UV) and the
// per-occurrence table of `BOPAlgo_WireSplitter::SplitBlock` plus its loops for
// the Z cylinder face, mirroring `BOPAlgo_Builder::BuildSplitFaces`
// (`BOPAlgo_Builder_2.cxx:359-494`).
#include <BOPAlgo_PaveFiller.hxx>
#include <BOPAlgo_WireEdgeSet.hxx>
#include <BOPAlgo_WireSplitter.hxx>
#include <BOPDS_DS.hxx>
#include <BOPDS_FaceInfo.hxx>
#include <BOPDS_PaveBlock.hxx>
#include <BOPTools_AlgoTools.hxx>
#include <BOPTools_AlgoTools2D.hxx>
#include <BOPTools_AlgoTools3D.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBndLib.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box2d.hxx>
#include <GProp_GProps.hxx>
#include <Geom2dAdaptor_Curve.hxx>
#include <Geom2d_Curve.hxx>
#include <GeomLib.hxx>
#include <NCollection_DataMap.hxx>
#include <NCollection_Sequence.hxx>
#include <Precision.hxx>
#include <Geom_Surface.hxx>
#include <IntTools_Context.hxx>
#include <Message_ProgressRange.hxx>
#include <Message_Report.hxx>
#include <NCollection_List.hxx>
#include <NCollection_Map.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Iterator.hxx>
#include <gp_Pnt2d.hxx>
#include <cmath>
#include <iostream>

namespace
{

const double THE_TWO_PI = 2.0 * M_PI;

const char* SurfName(GeomAbs_SurfaceType theT)
{
  switch (theT)
  {
    case GeomAbs_Plane: return "Plane";
    case GeomAbs_Cylinder: return "Cylinder";
    case GeomAbs_Cone: return "Cone";
    case GeomAbs_Sphere: return "Sphere";
    case GeomAbs_Torus: return "Torus";
    default: return "Other";
  }
}

gp_Pnt2d Coord2dOf(const TopoDS_Vertex& theV,
                   const TopoDS_Edge&   theE,
                   const TopoDS_Face&   theF,
                   double&              theT)
{
  double                    aFirst = 0., aLast = 0.;
  occ::handle<Geom2d_Curve> aC     = BRep_Tool::CurveOnSurface(theE, theF, aFirst, aLast);
  theT                             = BRep_Tool::Parameter(theV, theE, theF);
  gp_Pnt2d aP(99., 99.);
  if (!aC.IsNull())
  {
    aC->D0(theT, aP);
  }
  return aP;
}

void DumpEdgeReps(const TopoDS_Edge& theE, const TopoDS_Face& theF)
{
  double                    f = 0., l = 0.;
  occ::handle<Geom2d_Curve> aC = BRep_Tool::CurveOnSurface(theE, theF, f, l);
  std::cout << " pcRange=[" << f << "," << l << "]";
  if (!aC.IsNull())
  {
    gp_Pnt2d a0, a1;
    aC->D0(f, a0);
    aC->D0(l, a1);
    std::cout << " pc=(" << a0.X() << "," << a0.Y() << ")->(" << a1.X() << "," << a1.Y() << ")";
  }
}

// Fill `LE` exactly as `BOPAlgo_Builder_2.cxx:359-494`.
void CollectLE(BOPAlgo_PaveFiller&  thePF,
               const int                  theFace,
               NCollection_List<TopoDS_Shape>& theLE)
{
  const BOPDS_DS&         aDS  = thePF.DS();
  const TopoDS_Face       aF   = TopoDS::Face(aDS.Shape(theFace));
  const BOPDS_FaceInfo&   aFI  = aDS.FaceInfo(theFace);
  TopoDS_Face             aFF  = aF;
  aFF.Orientation(TopAbs_FORWARD);

  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aMFence;
  bool isChecked = false, isUClosed = false, isVClosed = false;

  for (TopExp_Explorer anExp(aFF, TopAbs_EDGE); anExp.More(); anExp.Next())
  {
    const TopoDS_Edge        aE     = TopoDS::Edge(anExp.Current());
    const TopAbs_Orientation anOriE = aE.Orientation();
    const int                nE     = aDS.Index(aE);
    if (nE <= 0 || !aDS.HasPaveBlocks(nE))
    {
      if (anOriE == TopAbs_INTERNAL)
      {
        TopoDS_Shape aEE = aE;
        aEE.Orientation(TopAbs_FORWARD);
        theLE.Append(aEE);
        aEE.Orientation(TopAbs_REVERSED);
        theLE.Append(aEE);
      }
      else
      {
        theLE.Append(aE);
      }
      continue;
    }
    if (!isChecked)
    {
      const occ::handle<Geom_Surface> aSurf = BRep_Tool::Surface(aFF);
      GeomLib::IsClosed(aSurf, BRep_Tool::Tolerance(aE), isUClosed, isVClosed);
      isChecked = true;
    }
    bool bIsClosed = false;
    if ((isUClosed || isVClosed) && BRep_Tool::IsClosed(aE, aFF))
    {
      bool isUIso = false, isVIso = false;
      BOPTools_AlgoTools2D::IsEdgeIsoline(aE, aFF, isUIso, isVIso);
      bIsClosed = ((isUClosed && isUIso) || (isVClosed && isVIso));
    }
    const bool                                           bIsDegenerated = BRep_Tool::Degenerated(aE);
    const NCollection_List<occ::handle<BOPDS_PaveBlock>>& aLPB          = aDS.PaveBlocks(nE);
    for (NCollection_List<occ::handle<BOPDS_PaveBlock>>::Iterator anIt(aLPB); anIt.More();
         anIt.Next())
    {
      TopoDS_Shape aSp = aDS.Shape(anIt.Value()->Edge());
      if (aSp.IsNull())
      {
        continue;
      }
      if (bIsDegenerated)
      {
        aSp.Orientation(anOriE);
        theLE.Append(aSp);
        continue;
      }
      if (anOriE == TopAbs_INTERNAL)
      {
        aSp.Orientation(TopAbs_FORWARD);
        theLE.Append(aSp);
        aSp.Orientation(TopAbs_REVERSED);
        theLE.Append(aSp);
        continue;
      }
      if (bIsClosed)
      {
        if (aMFence.Add(aSp))
        {
          if (!BRep_Tool::IsClosed(TopoDS::Edge(aSp), aFF))
          {
            if (!BOPTools_AlgoTools3D::DoSplitSEAMOnFace(TopoDS::Edge(aSp), aFF))
            {
              BOPTools_AlgoTools3D::DoSplitSEAMOnFace(aE, TopoDS::Edge(aSp), aFF);
            }
          }
          aSp.Orientation(TopAbs_FORWARD);
          theLE.Append(aSp);
          aSp.Orientation(TopAbs_REVERSED);
          theLE.Append(aSp);
        }
        continue;
      }
      aSp.Orientation(anOriE);
      occ::handle<Message_Report> aRep = new Message_Report;
      if (BOPTools_AlgoTools::IsSplitToReverseWithWarn(aSp, aE, nullptr, aRep))
      {
        aSp.Reverse();
      }
      theLE.Append(aSp);
    }
  }
  for (int aPass = 0; aPass < 2; ++aPass)
  {
    const NCollection_IndexedMap<occ::handle<BOPDS_PaveBlock>>& aMPB =
      aPass == 0 ? aFI.PaveBlocksIn() : aFI.PaveBlocksSc();
    for (int j = 1; j <= aMPB.Extent(); ++j)
    {
      TopoDS_Shape aSp = aDS.Shape(aMPB(j)->Edge());
      if (aSp.IsNull())
      {
        continue;
      }
      aSp.Orientation(TopAbs_FORWARD);
      theLE.Append(aSp);
      aSp.Orientation(TopAbs_REVERSED);
      theLE.Append(aSp);
    }
  }
}

void DumpFaceSplit(BOPAlgo_PaveFiller& thePF, const int theFace)
{
  const BOPDS_DS&       aDS = thePF.DS();
  const TopoDS_Face     aF  = TopoDS::Face(aDS.Shape(theFace));
  TopoDS_Face           aFF = aF;
  aFF.Orientation(TopAbs_FORWARD);
  BRepAdaptor_Surface aBAS(aFF);
  const BOPDS_FaceInfo& aFI = aDS.FaceInfo(theFace);
  std::cout << "\n== FACE ds=" << theFace << " " << SurfName(aBAS.GetType())
            << " orient=" << (int)aF.Orientation() << " in=" << aFI.PaveBlocksIn().Extent()
            << " on=" << aFI.PaveBlocksOn().Extent() << " sc=" << aFI.PaveBlocksSc().Extent()
            << " uclosed=" << (aBAS.IsUClosed() ? 1 : 0)
            << " vclosed=" << (aBAS.IsVClosed() ? 1 : 0) << "\n";

  NCollection_List<TopoDS_Shape> aLE;
  CollectLE(thePF, theFace, aLE);
  std::cout << " LE=" << aLE.Extent() << "\n";

  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aMS;
  for (NCollection_List<TopoDS_Shape>::Iterator anIt(aLE); anIt.More(); anIt.Next())
  {
    const TopoDS_Edge& aE        = TopoDS::Edge(anIt.Value());
    const bool         bIsClosed = BRep_Tool::Degenerated(aE) || BRep_Tool::IsClosed(aE, aFF);
    if (!aMS.Add(aE) && !bIsClosed)
    {
      aMS.Remove(aE);
    }
  }
  for (NCollection_List<TopoDS_Shape>::Iterator anIt(aLE); anIt.More(); anIt.Next())
  {
    const TopoDS_Edge& aE = TopoDS::Edge(anIt.Value());
    std::cout << "  E ds=" << aDS.Index(aE) << " ori=" << (int)aE.Orientation()
              << " inside=" << (aMS.Contains(aE) ? 0 : 1)
              << " deg=" << (BRep_Tool::Degenerated(aE) ? 1 : 0);
    DumpEdgeReps(aE, aFF);
    int anIv = 0;
    for (TopoDS_Iterator aVI(aE); aVI.More(); aVI.Next(), ++anIv)
    {
      TopoDS_Vertex aV = TopoDS::Vertex(aVI.Value());
      double        aT = 0.;
      const gp_Pnt2d aP = Coord2dOf(aV, aE, aFF, aT);
      std::cout << "  V" << anIv << " ds=" << aDS.Index(aV)
                << " ori=" << (int)aV.Orientation() << " uv=(" << aP.X() << "," << aP.Y() << ")";
    }
    std::cout << "\n";
  }

  // OCCT's own loop decomposition.
  occ::handle<IntTools_Context> aCtx = new IntTools_Context;
  BOPAlgo_WireEdgeSet aWES;
  aWES.SetFace(aFF);
  for (NCollection_List<TopoDS_Shape>::Iterator anIt(aLE); anIt.More(); anIt.Next())
  {
    aWES.AddStartElement(anIt.Value());
  }
  BOPAlgo_WireSplitter aWS;
  aWS.SetWES(aWES);
  aWS.SetContext(aCtx);
  aWS.SetRunParallel(false);
  aWS.Perform(Message_ProgressRange());
  if (aWS.HasErrors())
  {
    std::cout << " splitter errors\n";
    aWS.DumpErrors(std::cout);
  }
  int aNw = 0;
  for (NCollection_List<TopoDS_Shape>::Iterator anIt(aWES.Shapes()); anIt.More(); anIt.Next())
  {
    std::cout << " LOOP " << aNw++ << ":";
    for (TopoDS_Iterator aEI(anIt.Value()); aEI.More(); aEI.Next())
    {
      const TopoDS_Edge& aE = TopoDS::Edge(aEI.Value());
      std::cout << " " << aDS.Index(aE) << (aE.Orientation() == TopAbs_REVERSED ? "-" : "+");
    }
    std::cout << "\n";
  }
}

void DumpResult(const TopoDS_Shape& theS, const char* theTag)
{
  GProp_GProps aV;
  BRepGProp::VolumeProperties(theS, aV);
  std::cout << "\n=== " << theTag << ": volume=" << aV.Mass() << " faces="
            << [&] {
                 int n = 0;
                 for (TopExp_Explorer anEx(theS, TopAbs_FACE); anEx.More(); anEx.Next())
                 {
                   ++n;
                 }
                 return n;
               }()
            << "\n";
  int aI = 0;
  for (TopExp_Explorer anEx(theS, TopAbs_FACE); anEx.More(); anEx.Next(), ++aI)
  {
    const TopoDS_Face aF = TopoDS::Face(anEx.Current());
    BRepAdaptor_Surface aS(aF);
    GProp_GProps aA;
    BRepGProp::SurfaceProperties(aF, aA);
    int  aNw = 0;
    int  aNe = 0;
    for (TopoDS_Iterator aWI(aF); aWI.More(); aWI.Next(), ++aNw)
    {
      for (TopoDS_Iterator aEI(aWI.Value()); aEI.More(); aEI.Next())
      {
        ++aNe;
      }
    }
    Bnd_Box2d aUV;
    BRepTools::AddUVBounds(aF, aUV);
    double u1, u2, v1, v2;
    aUV.Get(u1, v1, u2, v2);
    std::cout << "  F" << aI << " " << SurfName(aS.GetType()) << " area=" << aA.Mass()
              << " wires=" << aNw << " edges=" << aNe << " uv=(" << u1 << "," << u2 << ","
              << v1 << "," << v2 << ")\n";
  }
}

// ============================================================================
// `BOPAlgo_WireSplitter` internals, copied from `BOPAlgo_WireSplitter_1.cxx`
// with prints: the smart map, the angle table (`SplitBlock` steps 1-3, `:136-321`),
// `RefineAngles` (`:903-1029`) and `Path` (`:355-617`).
// ============================================================================

const double THE_EPS = 2.2204460492503131e-16;

struct EInfo
{
  TopoDS_Edge edge;
  bool        passed   = false;
  bool        isIn     = false;
  bool        isInside = false;
  double      angle    = -1.;
};

double AngleOf2dDir(const gp_Dir2d& theDir)
{
  gp_Dir2d aRef(gp_Dir2d::D::X);
  double   anAngle = aRef.Angle(theDir);
  if (anAngle < 0.)
  {
    anAngle += M_PI + M_PI;
  }
  return anAngle;
}

double Tol2D(const TopoDS_Vertex& theV, const BRepAdaptor_Surface& theGAS)
{
  const double aTolV3D = BRep_Tool::Tolerance(theV);
  const double anUr    = theGAS.UResolution(aTolV3D);
  const double aVr     = theGAS.VResolution(aTolV3D);
  double       aTol2D  = (aVr > anUr) ? aVr : anUr;
  if (aTol2D < aTolV3D)
  {
    aTol2D = aTolV3D;
  }
  if (theGAS.GetType() == GeomAbs_BSplineSurface)
  {
    aTol2D = 1.1 * aTol2D;
  }
  return aTol2D;
}

double Angle2DOf(const TopoDS_Vertex&       theV,
                 const TopoDS_Edge&         theE,
                 const TopoDS_Face&         theF,
                 const BRepAdaptor_Surface& theGAS,
                 const bool                 theIsIn,
                 const occ::handle<IntTools_Context>& theCtx,
                 double&                    theTV,
                 gp_Pnt2d&                  thePV)
{
  double                    aFirst = 0., aLast = 0., aToler = 0.;
  occ::handle<Geom2d_Curve> aC2D;
  theTV = BRep_Tool::Parameter(theV, theE, theF);
  if (Precision::IsInfinite(theTV))
  {
    return 0.;
  }
  BOPTools_AlgoTools2D::CurveOnSurface(theE, theF, aC2D, aFirst, aLast, aToler, theCtx);
  std::cout << "  a2d f=" << aFirst << " l=" << aLast << " tv=" << theTV
            << " null=" << (aC2D.IsNull() ? 1 : 0) << "\n";
  const double tol2d = 2. * Tol2D(theV, theGAS);
  Geom2dAdaptor_Curve aGAC2D(aC2D);
  double              dt = std::max(aGAC2D.Resolution(tol2d), Precision::PConfusion());
  std::cout << "  a2d tol2d=" << tol2d << " dt=" << dt << " type=" << (int)aGAC2D.GetType()
            << "\n";
  if (aGAC2D.GetType() != GeomAbs_Line)
  {
    gp_Pnt2d aP;
    gp_Vec2d aD1, aD2;
    aC2D->D2(theTV, aP, aD1, aD2);
    const double aN1 = aD1.Magnitude();
    if (aN1 > Precision::PConfusion())
    {
      double R = std::abs(aD1.Crossed(aD2)) / (aN1 * aN1 * aN1);
      if (R > Precision::PConfusion())
      {
        R             = 1. / R;
        const double cosphi = R / (R + tol2d);
        dt            = std::max(dt, std::acos(cosphi));
      }
    }
  }
  double aTX = 0.05 * (aLast - aFirst);
  if (aTX < 5.e-5)
  {
    aTX = std::min(5.e-5, (aLast - aFirst) / 2.);
  }
  if (dt > aTX)
  {
    dt = aTX;
  }
  const double aTV1 = (std::abs(theTV - aFirst) < std::abs(theTV - aLast)) ? theTV + dt : theTV - dt;
  gp_Pnt2d     aPV1;
  aGAC2D.D0(aTV1, aPV1);
  aGAC2D.D0(theTV, thePV);
  gp_Vec2d aV2D = theIsIn ? gp_Vec2d(aPV1, thePV) : gp_Vec2d(thePV, aPV1);
  return AngleOf2dDir(gp_Dir2d(aV2D));
}

gp_Pnt2d Coord2dOf2(const TopoDS_Vertex& theV, const TopoDS_Edge& theE, const TopoDS_Face& theF)
{
  double                    aFirst = 0., aLast = 0.;
  occ::handle<Geom2d_Curve> aC2D   = BRep_Tool::CurveOnSurface(theE, theF, aFirst, aLast);
  const double              aT     = BRep_Tool::Parameter(theV, theE, theF);
  gp_Pnt2d                  aP(99., 99.);
  if (!aC2D.IsNull())
  {
    aC2D->D0(aT, aP);
  }
  return aP;
}

gp_Pnt2d Coord2dVfOf(const TopoDS_Edge& theE, const TopoDS_Face& theF)
{
  for (TopoDS_Iterator anIt(theE); anIt.More(); anIt.Next())
  {
    if (anIt.Value().Orientation() == TopAbs_FORWARD)
    {
      return Coord2dOf2(TopoDS::Vertex(anIt.Value()), theE, theF);
    }
  }
  return gp_Pnt2d(99., 99.);
}

void GetNextVert(const TopoDS_Vertex& theV, const TopoDS_Edge& theE, TopoDS_Vertex& theV1)
{
  for (TopoDS_Iterator anIt(theE); anIt.More(); anIt.Next())
  {
    if (!anIt.Value().IsEqual(theV))
    {
      theV1 = TopoDS::Vertex(anIt.Value());
      return;
    }
  }
  theV1 = theV;
}

double ClockAngle(const double theAngleIn, const double theAngleOut)
{
  double aIn = theAngleIn, aOut = theAngleOut;
  if (aIn >= THE_TWO_PI)
  {
    aIn -= THE_TWO_PI;
  }
  if (aOut >= THE_TWO_PI)
  {
    aOut -= THE_TWO_PI;
  }
  double a1 = aIn + M_PI;
  if (a1 >= THE_TWO_PI)
  {
    a1 -= THE_TWO_PI;
  }
  double dA = a1 - aOut;
  if (dA <= 0.)
  {
    dA = THE_TWO_PI + dA;
  }
  else if (dA <= 1.e-14)
  {
    dA = THE_TWO_PI;
  }
  return dA;
}

double AngIn(const TopoDS_Edge& theEIn, const NCollection_List<EInfo>& theLE)
{
  for (NCollection_List<EInfo>::Iterator anIt(theLE); anIt.More(); anIt.Next())
  {
    if (anIt.Value().isIn && anIt.Value().edge == theEIn)
    {
      return anIt.Value().angle;
    }
  }
  return 0.;
}

int NbsWaysOut(const NCollection_List<EInfo>& theLE)
{
  int aCnt = 0;
  for (NCollection_List<EInfo>::Iterator anIt(theLE); anIt.More(); anIt.Next())
  {
    if (!anIt.Value().isIn && !anIt.Value().passed)
    {
      ++aCnt;
    }
  }
  return aCnt;
}

typedef NCollection_IndexedDataMap<TopoDS_Shape,
                                   NCollection_List<EInfo>,
                                   TopTools_ShapeMapHasher>
  MySmartMap;
typedef NCollection_DataMap<TopoDS_Shape, bool, TopTools_ShapeMapHasher> MyVertMap;

int DsIdx(const BOPDS_DS& theDS, const TopoDS_Shape& theS) { return theDS.Index(theS); }

void PrintPnt(const BOPDS_DS& theDS, const TopoDS_Vertex& theV)
{
  const gp_Pnt aP = BRep_Tool::Pnt(theV);
  std::cout << "V" << DsIdx(theDS, theV) << "(" << aP.X() << "," << aP.Y() << "," << aP.Z() << ")";
}

void PrintOcc(const BOPDS_DS& theDS, const EInfo& theEI)
{
  std::cout << "e" << DsIdx(theDS, theEI.edge)
            << (theEI.edge.Orientation() == TopAbs_REVERSED ? "-" : "+")
            << " in=" << (theEI.isIn ? 1 : 0) << " ins=" << (theEI.isInside ? 1 : 0)
            << " pass=" << (theEI.passed ? 1 : 0) << " ang=" << theEI.angle;
}

void PathTrace(const BRepAdaptor_Surface&                      theGAS,
               const TopoDS_Face&                              theFace,
               const MyVertMap&                                theVertMap,
               const BOPDS_DS&                                 theDS,
               const TopoDS_Vertex&                            theVFirst,
               const TopoDS_Edge&                              theEFirst,
               EInfo&                                          theEIFirst,
               MySmartMap&                                     theMap,
               int                                             theRun)
{
  NCollection_Sequence<TopoDS_Shape> aLS;
  NCollection_Sequence<TopoDS_Shape> aVertVa;
  NCollection_Sequence<gp_Pnt2d>     aCoordVa;
  NCollection_Sequence<EInfo*>       anInfoSeq;

  TopoDS_Vertex aVa = theVFirst;
  TopoDS_Edge   aEOuta = theEFirst;
  EInfo*        anEdgeInfo = &theEIFirst;
  const double  eps = THE_EPS;

  for (;;)
  {
    int aNb = aLS.Length();
    if (aNb == 1 && aLS(aNb).IsSame(aEOuta))
    {
      std::cout << "  run" << theRun << " ABANDON (came back to the first edge)\n";
      return;
    }
    anEdgeInfo->passed = true;
    aLS.Append(aEOuta);
    aVertVa.Append(aVa);
    anInfoSeq.Append(anEdgeInfo);

    TopoDS_Vertex pVa = aVa;
    pVa.Orientation(TopAbs_FORWARD);
    const gp_Pnt2d aPa = Coord2dOf2(pVa, aEOuta, theFace);
    aCoordVa.Append(aPa);
    std::cout << "  run" << theRun << " stp" << aLS.Length() << " " << DsIdx(theDS, aEOuta)
              << (aEOuta.Orientation() == TopAbs_REVERSED ? "-" : "+") << " V"
              << DsIdx(theDS, aVa) << " pa=(" << aPa.X() << "," << aPa.Y() << ")";

    TopoDS_Vertex aVb;
    GetNextVert(pVa, aEOuta, aVb);
    const gp_Pnt2d aPb = Coord2dOf2(aVb, aEOuta, theFace);
    std::cout << " -> V" << DsIdx(theDS, aVb) << " pb=(" << aPb.X() << "," << aPb.Y() << ")";

    const NCollection_List<EInfo>& aLEInfo = theMap.FindFromKey(aVb);
    double aTol2D  = 2. * Tol2D(aVb, theGAS);
    double aTol2D2 = aTol2D * aTol2D;
    const bool bIsClosed = theVertMap.IsBound(aVb) ? theVertMap.Find(aVb) : false;
    std::cout << " closed=" << (bIsClosed ? 1 : 0) << "\n";
    {
      NCollection_List<TopoDS_Shape> aBuf;
      bool                           bHasEdge = false;
      aNb = aLS.Length();
      for (int i = aNb; i > 0; --i)
      {
        const TopoDS_Shape& aVPrev  = aVertVa(i);
        const gp_Pnt2d&     aPaPrev = aCoordVa(i);
        const TopoDS_Edge&  aEPrev  = TopoDS::Edge(aLS(i));
        aBuf.Append(aEPrev);
        if (!bHasEdge)
        {
          bHasEdge = !BRep_Tool::Degenerated(aEPrev);
          if (!bHasEdge)
          {
            continue;
          }
        }
        bool anIsSameV   = aVPrev.IsSame(aVb);
        bool anIsSameV2d = anIsSameV;
        if (anIsSameV && bIsClosed)
        {
          const double aD2 = aPaPrev.SquareDistance(aPb);
          anIsSameV2d      = aD2 < aTol2D2;
          std::cout << "   close? i=" << i << " d2=" << aD2 << " tol2=" << aTol2D2
                    << " -> " << (anIsSameV2d ? 1 : 0) << "\n";
          if (anIsSameV2d)
          {
            const double udist = std::abs(aPaPrev.X() - aPb.X());
            const double vdist = std::abs(aPaPrev.Y() - aPb.Y());
            const double aTolU = 2. * theGAS.UResolution(BRep_Tool::Tolerance(aVb));
            const double aTolV = 2. * theGAS.VResolution(BRep_Tool::Tolerance(aVb));
            if ((udist > aTolU) || (vdist > aTolV))
            {
              anIsSameV2d = false;
            }
            std::cout << "   uvclose i=" << i << " du=" << udist << " dv=" << vdist
                      << " -> " << (anIsSameV2d ? 1 : 0) << "\n";
          }
        }
        else if (anIsSameV)
        {
          std::cout << "   close? i=" << i << " (not closed on face) -> 1\n";
        }
        if (anIsSameV && anIsSameV2d)
        {
          int iPriz = 1;
          if (aBuf.Extent() == 2 && aBuf.First().IsSame(aBuf.Last()))
          {
            iPriz = 0;
          }
          if (iPriz)
          {
            std::cout << "   LOOP:";
            for (NCollection_List<TopoDS_Shape>::Iterator aBI(aBuf); aBI.More(); aBI.Next())
            {
              const TopoDS_Edge& aBE = TopoDS::Edge(aBI.Value());
              std::cout << " " << DsIdx(theDS, aBE)
                        << (aBE.Orientation() == TopAbs_REVERSED ? "-" : "+");
            }
            std::cout << "\n";
          }
          const int aNbj = i - 1;
          if (aNbj < 1)
          {
            aLS.Clear();
            aVertVa.Clear();
            aCoordVa.Clear();
            return;
          }
          NCollection_Sequence<TopoDS_Shape> aLSt, aVertVat;
          NCollection_Sequence<gp_Pnt2d>     aCoordVat;
          NCollection_Sequence<EInfo*>       anInfoSeqTmp;
          aVb = TopoDS::Vertex(aVertVa(i));
          for (int j = 1; j <= aNbj; ++j)
          {
            aLSt.Append(aLS(j));
            aVertVat.Append(aVertVa(j));
            aCoordVat.Append(aCoordVa(j));
            anInfoSeqTmp.Append(anInfoSeq(j));
          }
          aLS       = aLSt;
          aVertVa   = aVertVat;
          aCoordVa  = aCoordVat;
          anInfoSeq = anInfoSeqTmp;
          aEOuta     = TopoDS::Edge(aLS.Last());
          anEdgeInfo = anInfoSeq.Last();
          std::cout << "   cut back to " << aLS.Length() << " edges\n";
          break;
        }
      }
    }

    EInfo* pEdgeInfo = nullptr;
    const double anAngleIn = AngIn(aEOuta, aLEInfo);
    double       aMinAngle = 100.;
    int          iCnt      = NbsWaysOut(aLEInfo);
    const bool   isBoundary = !anEdgeInfo->isInside;
    int          aNbWaysInside = 0;
    EInfo*       pOnlyWayIn = nullptr;
    std::cout << "   at V" << DsIdx(theDS, aVb) << " angleIn=" << anAngleIn << " nout=" << iCnt
              << " bnd=" << (isBoundary ? 1 : 0) << "\n";
    for (NCollection_List<EInfo>::Iterator anIt(aLEInfo); anIt.More(); anIt.Next())
    {
      EInfo&       anEI   = anIt.ChangeValue();
      const bool   anIsOut = !anEI.isIn;
      const bool   anIsNotPassed = !anEI.passed;
      if (!(anIsOut && anIsNotPassed))
      {
        continue;
      }
      std::cout << "    cand ";
      PrintOcc(theDS, anEI);
      if (!iCnt)
      {
        std::cout << " -> no way out\n";
        return;
      }
      if (iCnt == 1)
      {
        pEdgeInfo = &anEI;
        std::cout << " -> only\n";
        break;
      }
      double anAngle;
      if (anEI.edge.IsSame(aEOuta))
      {
        anAngle = THE_TWO_PI;
      }
      else
      {
        if (bIsClosed)
        {
          const gp_Pnt2d aP2Dx = Coord2dVfOf(anEI.edge, theFace);
          if (aP2Dx.SquareDistance(aPb) > aTol2D2)
          {
            std::cout << " -> skip (vf d2)\n";
            continue;
          }
        }
        anAngle = ClockAngle(anAngleIn, anEI.angle);
      }
      if (isBoundary && anEI.isInside)
      {
        ++aNbWaysInside;
        pOnlyWayIn = &anEI;
      }
      std::cout << " cw=" << anAngle << "\n";
      if (anAngle < aMinAngle - eps)
      {
        aMinAngle = anAngle;
        pEdgeInfo = &anEI;
      }
    }
    if (aNbWaysInside == 1)
    {
      pEdgeInfo = pOnlyWayIn;
      std::cout << "   only inside\n";
    }
    if (!pEdgeInfo)
    {
      std::cout << "  run" << theRun << " ABANDON (no way to go)\n";
      return;
    }
    std::cout << "   chose e" << DsIdx(theDS, pEdgeInfo->edge)
              << (pEdgeInfo->edge.Orientation() == TopAbs_REVERSED ? "-" : "+") << "\n";
    aVa        = aVb;
    aEOuta     = pEdgeInfo->edge;
    anEdgeInfo = pEdgeInfo;
  }
}

void DumpTrace(BOPAlgo_PaveFiller& thePF, const int theFace)
{
  const BOPDS_DS& theDS = thePF.DS();
  TopoDS_Face     aFF   = TopoDS::Face(theDS.Shape(theFace));
  aFF.Orientation(TopAbs_FORWARD);
  BRepAdaptor_Surface aBAS(aFF);
  occ::handle<IntTools_Context> aCtx = new IntTools_Context;

  NCollection_List<TopoDS_Shape> aLE;
  CollectLE(thePF, theFace, aLE);
  std::cout << " TRACE ds=" << theFace << " LE=" << aLE.Extent() << "\n";

  MySmartMap  aMap;
  MyVertMap   aVertMap;
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aMS;
  for (NCollection_List<TopoDS_Shape>::Iterator anIt(aLE); anIt.More(); anIt.Next())
  {
    const TopoDS_Edge& aE = TopoDS::Edge(anIt.Value());
    std::cout << "  m ds=" << theDS.Index(aE) << " ori=" << (int)aE.Orientation() << "\n";
    if (!BOPTools_AlgoTools2D::HasCurveOnSurface(aE, aFF))
    {
      continue;
    }
    bool bIsClosed = BRep_Tool::Degenerated(aE) || BRep_Tool::IsClosed(aE, aFF);
    if (!aMS.Add(aE) && !bIsClosed)
    {
      aMS.Remove(aE);
    }
    TopoDS_Shape aV1;
    for (TopoDS_Iterator aItS(aE); aItS.More(); aItS.Next())
    {
      const TopoDS_Shape& aV  = aItS.Value();
      const int           aIx = aMap.FindIndex(aV);
      if (!aIx)
      {
        NCollection_List<EInfo> aLEIx;
        aMap.Add(aV, aLEIx);
      }
      NCollection_List<EInfo>& aLEI = aMap(aMap.FindIndex(aV));
      EInfo                    aEI;
      aEI.edge   = aE;
      aEI.isIn   = (aV.Orientation() == TopAbs_REVERSED);
      aLEI.Append(aEI);
      if (aV1.IsNull())
      {
        aV1 = aV;
      }
      else
      {
        bIsClosed = bIsClosed || aV1.IsSame(aV);
      }
      if (aVertMap.IsBound(aV))
      {
        if (bIsClosed)
        {
          aVertMap.ChangeFind(aV) = bIsClosed;
        }
      }
      else
      {
        aVertMap.Bind(aV, bIsClosed);
      }
    }
  }

  for (int i = 1; i <= aMap.Extent(); ++i)
  {
    for (NCollection_List<EInfo>::Iterator anIt(aMap(i)); anIt.More(); anIt.Next())
    {
      EInfo& anEI = anIt.ChangeValue();
      anEI.isInside = !aMS.Contains(anEI.edge);
    }
  }

  // step 3: angles
  for (int i = 1; i <= aMap.Extent(); ++i)
  {
    TopoDS_Vertex aV = TopoDS::Vertex(aMap.FindKey(i));
    for (NCollection_List<EInfo>::Iterator anIt(aMap(i)); anIt.More(); anIt.Next())
    {
      EInfo& anEI = anIt.ChangeValue();
      TopoDS_Vertex aVV = aV;
      aVV.Orientation(anEI.isIn ? TopAbs_REVERSED : TopAbs_FORWARD);
      double   aTV = 0.;
      gp_Pnt2d aPV;
      std::cout << " T V" << DsIdx(theDS, aV) << " e" << DsIdx(theDS, anEI.edge) << " ins="
                << (anEI.isInside ? 1 : 0) << " in=" << (anEI.isIn ? 1 : 0) << "\n";
      anEI.angle = Angle2DOf(aVV, anEI.edge, aFF, aBAS, anEI.isIn, aCtx, aTV, aPV);
      std::cout << " T2 V" << DsIdx(theDS, aV) << " ";
      PrintOcc(theDS, anEI);
      std::cout << " tv=" << aTV << " pv=(" << aPV.X() << "," << aPV.Y() << ")\n";
    }
  }

  // RefineAngles
  for (int i = 1; i <= aMap.Extent(); ++i)
  {
    TopoDS_Vertex              aV = TopoDS::Vertex(aMap.FindKey(i));
    NCollection_List<EInfo>&   aLEI = aMap(i);
    double aA1 = 0., aA2 = 0.;
    int    iCntBnd = 0, iCntInt = 0;
    for (NCollection_List<EInfo>::Iterator anIt(aLEI); anIt.More(); anIt.Next())
    {
      const EInfo& aEI = anIt.Value();
      if (!aEI.isInside)
      {
        ++iCntBnd;
        if (!aEI.isIn)
        {
          aA1 = aEI.angle;
        }
        else
        {
          aA2 = aEI.angle;
        }
      }
      else
      {
        ++iCntInt;
      }
    }
    std::cout << " R V" << DsIdx(theDS, aV) << " bnd=" << iCntBnd << " int=" << iCntInt
              << " a1=" << aA1 << " a2=" << aA2;
    if (iCntBnd != 2)
    {
      std::cout << " (skip)\n";
      continue;
    }
    std::cout << " delta=" << ClockAngle(aA2, aA1) << "\n";
  }

  // Path
  int aRun = 0;
  for (int i = 1; i <= aMap.Extent(); ++i)
  {
    TopoDS_Vertex aVa = TopoDS::Vertex(aMap.FindKey(i));
    for (NCollection_List<EInfo>::Iterator anIt(aMap(i)); anIt.More(); anIt.Next())
    {
      EInfo& anEI = anIt.ChangeValue();
      if (!anEI.isIn && !anEI.passed)
      {
        std::cout << " RUN" << aRun << " start V" << DsIdx(theDS, aVa) << " e"
                  << DsIdx(theDS, anEI.edge) << "\n";
        PathTrace(aBAS, aFF, aVertMap, theDS, aVa, anEI.edge, anEI, aMap, aRun);
        ++aRun;
      }
    }
  }
  std::cout << " RUNS=" << aRun << "\n";
}

} // namespace

int main()
{
  std::cout.setf(std::ios::unitbuf);
  const TopoDS_Shape aBox =
    BRepPrimAPI_MakeBox(gp_Pnt(0., 0., 0.), gp_Pnt(2., 2., 1.)).Shape();
  const TopoDS_Shape aCylZ =
    BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(1., 1., -0.5), gp_Dir(0., 0., 1.)), 0.3, 2.0).Shape();
  const TopoDS_Shape aCylY =
    BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(1., -2., 0.5), gp_Dir(0., 1., 0.)), 0.3, 4.0).Shape();

  BRepAlgoAPI_Cut aCut1(aBox, aCylZ);
  aCut1.Build();
  if (!aCut1.IsDone())
  {
    std::cout << "cut1 failed\n";
    return 1;
  }
  DumpResult(aCut1.Shape(), "box - Zcyl");

  BRepAlgoAPI_Cut aCut2(aCut1.Shape(), aCylY);
  aCut2.Build();
  if (!aCut2.IsDone())
  {
    std::cout << "cut2 failed\n";
    return 1;
  }
  DumpResult(aCut2.Shape(), "box - Zcyl - Ycyl");

  BOPAlgo_PaveFiller   aPF;
  TopTools_ListOfShape aArgs;
  aArgs.Append(aCut1.Shape());
  aArgs.Append(aCylY);
  aPF.SetArguments(aArgs);
  aPF.SetRunParallel(false);
  aPF.Perform(Message_ProgressRange());
  if (aPF.HasErrors())
  {
    std::cout << "pave filler errors\n";
    aPF.DumpErrors(std::cout);
    return 1;
  }
  const BOPDS_DS& aDS = aPF.DS();
  for (int i = 0; i < aDS.NbShapes(); ++i)
  {
    if (aDS.Shape(i).IsNull() || aDS.Shape(i).ShapeType() != TopAbs_FACE)
    {
      continue;
    }
    TopoDS_Face         aF = TopoDS::Face(aDS.Shape(i));
    BRepAdaptor_Surface aS(aF);
    if (aS.GetType() != GeomAbs_Cylinder)
    {
      continue;
    }
    if (!aDS.HasFaceInfo(i))
    {
      continue;
    }
    DumpTrace(aPF, i);
    DumpFaceSplit(aPF, i);
  }
  return 0;
}
