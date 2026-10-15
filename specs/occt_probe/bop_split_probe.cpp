// T-79 ground truth: run OCCT's own BOPAlgo_WireSplitter on the cylinder lateral
// face of the box[-1,1]^3 | cylinder(r=0.4,h=2) fuse and dump the edge set `LE`
// exactly as `BOPAlgo_Builder::BuildSplitFaces` (`BOPAlgo_Builder_2.cxx:359-494`)
// collects it, together with the per-occurrence data `SplitBlock`
// (`BOPAlgo_WireSplitter_1.cxx:113-354`) derives from it.
//
// Everything here mirrors the OCCT source named in the comments; nothing is
// re-derived from the Rust port.
#include <BOPAlgo_PaveFiller.hxx>
#include <BOPAlgo_WireSplitter.hxx>
#include <BOPAlgo_WireEdgeSet.hxx>
#include <BOPDS_DS.hxx>
#include <BOPDS_FaceInfo.hxx>
#include <BOPDS_PaveBlock.hxx>
#include <BOPTools_AlgoTools.hxx>
#include <BOPTools_AlgoTools2D.hxx>
#include <BOPTools_AlgoTools3D.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRep_Builder.hxx>
#include <BRep_CurveOnSurface.hxx>
#include <BRep_GCurve.hxx>
#include <BRep_TEdge.hxx>
#include <BRep_Tool.hxx>
#include <Geom2dAdaptor_Curve.hxx>
#include <Geom2d_Curve.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <Geom_Circle.hxx>
#include <Geom_Curve.hxx>
#include <GeomLib.hxx>
#include <GeomLProp_CLProps.hxx>
#include <Geom_Surface.hxx>
#include <IntTools_Context.hxx>
#include <Message_ProgressRange.hxx>
#include <Message_Report.hxx>
#include <NCollection_List.hxx>
#include <NCollection_Map.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Iterator.hxx>
#include <gp_Dir2d.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Vec2d.hxx>
#include <algorithm>
#include <cmath>
#include <iostream>

namespace
{

const double THE_TWO_PI = 2.0 * M_PI;

double Tol2D(const TopoDS_Vertex& theV, const BRepAdaptor_Surface& theS)
{
  const double aTolV3D = BRep_Tool::Tolerance(theV);
  const double anUr    = theS.UResolution(aTolV3D);
  const double aVr     = theS.VResolution(aTolV3D);
  double       aTol2D  = (aVr > anUr) ? aVr : anUr;
  if (aTol2D < aTolV3D)
  {
    aTol2D = aTolV3D;
  }
  return aTol2D;
}

double AngleOfDir(const gp_Dir2d& theDir)
{
  const gp_Dir2d aRef(gp_Dir2d::D::X);
  double         anAngle = aRef.Angle(theDir);
  if (anAngle < 0.)
  {
    anAngle += THE_TWO_PI;
  }
  return anAngle;
}

// `BOPAlgo_WireSplitter_1.cxx:768-840`.
double Angle2D(const TopoDS_Vertex&                 theV,
               const TopoDS_Edge&                   theE,
               const TopoDS_Face&                   theF,
               BRepAdaptor_Surface&                 theS,
               const bool                           theIsIn,
               const occ::handle<IntTools_Context>& theCtx)
{
  const double aT = BRep_Tool::Parameter(theV, theE, theF);
  if (Precision::IsInfinite(aT))
  {
    return 0.;
  }
  occ::handle<Geom2d_Curve> aC2D;
  double                    aFirst = 0., aLast = 0., aToler = 0.;
  BOPTools_AlgoTools2D::CurveOnSurface(theE, theF, aC2D, aFirst, aLast, aToler, theCtx);
  if (aC2D.IsNull())
  {
    return 0.;
  }
  const double        aTol2d = 2. * Tol2D(theV, theS);
  Geom2dAdaptor_Curve aGAC2D(aC2D);
  double              dt = std::max(aGAC2D.Resolution(aTol2d), Precision::PConfusion());
  if (aGAC2D.GetType() != GeomAbs_Line)
  {
    GeomLProp_CLProps2d aLProp(aC2D, aT, 2, Precision::PConfusion());
    if (aLProp.IsTangentDefined())
    {
      double aR = aLProp.Curvature();
      if (aR > Precision::PConfusion())
      {
        aR = 1. / aR;
        const double aCosPhi = aR / (aR + aTol2d);
        dt = std::max(dt, std::acos(aCosPhi));
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
  const double aT1 = (std::fabs(aT - aFirst) < std::fabs(aT - aLast)) ? aT + dt : aT - dt;
  gp_Pnt2d aP1, aP;
  aGAC2D.D0(aT1, aP1);
  aGAC2D.D0(aT, aP);
  const gp_Vec2d aV2D = theIsIn ? gp_Vec2d(aP1, aP) : gp_Vec2d(aP, aP1);
  return AngleOfDir(gp_Dir2d(aV2D));
}

// `Coord2d` (`_1.cxx:663-674`): `BRep_Tool::Parameter(V, E, F)` + the p-curve
// selected by `BRep_Tool::CurveOnSurface(E, F, first, last)`.
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

// 3D curve type of the edge plus, for a circle, its `gp_Ax2` — which fixes
// where parameter 0 sits on the surface.
void DumpCurve3d(const TopoDS_Edge& theE)
{
  double                  f = 0., l = 0.;
  occ::handle<Geom_Curve> aC = BRep_Tool::Curve(theE, f, l);
  if (aC.IsNull())
  {
    std::cout << " c3d=null";
    return;
  }
  GeomAdaptor_Curve aGA(aC);
  std::cout << " c3d=" << (int)aGA.GetType();
  if (aGA.GetType() == GeomAbs_Circle)
  {
    const gp_Circ aCi  = aGA.Circle();
    const gp_Ax2 anAx  = aCi.Position();
    std::cout << " R=" << aCi.Radius() << " loc=(" << anAx.Location().X() << ","
              << anAx.Location().Y() << "," << anAx.Location().Z() << ")"
              << " xd=(" << anAx.XDirection().X() << "," << anAx.XDirection().Y() << ","
              << anAx.XDirection().Z() << ")"
              << " yd=(" << anAx.YDirection().X() << "," << anAx.YDirection().Y() << ","
              << anAx.YDirection().Z() << ")";
  }
  const gp_Pnt aP0 = aC->Value(f);
  const gp_Pnt aP1 = aC->Value(l);
  std::cout << " p(f)=(" << aP0.X() << "," << aP0.Y() << "," << aP0.Z() << ")"
            << " p(l)=(" << aP1.X() << "," << aP1.Y() << "," << aP1.Z() << ")";
}

// Number of p-curves on this surface and both of them, when OCCT keeps two
// (`BRep_Tool::IsClosed(edge, face)` == the representation is a
// `BRep_CurveOnClosedSurface`).
void DumpRepresentations(const TopoDS_Edge& theE, const TopoDS_Face& theF)
{
  TopLoc_Location                  aLoc;
  const occ::handle<Geom_Surface>& aS    = BRep_Tool::Surface(theF, aLoc);
  const TopLoc_Location            aLocE = aLoc.Predivided(theE.Location());
  const BRep_TEdge*                aTE   = static_cast<const BRep_TEdge*>(theE.TShape().get());
  int                              aNbU  = 0;
  const BRep_CurveOnSurface*       aRep  = nullptr;
  for (NCollection_List<occ::handle<BRep_CurveRepresentation>>::Iterator anIt(aTE->Curves());
       anIt.More();
       anIt.Next())
  {
    const occ::handle<BRep_CurveRepresentation>& aCr = anIt.Value();
    if (aCr->IsCurveOnSurface(aS, aLocE))
    {
      ++aNbU;
      aRep = dynamic_cast<const BRep_CurveOnSurface*>(aCr.get());
    }
  }
  std::cout << " reps=" << aNbU << " IsClosed(E,F)=" << (BRep_Tool::IsClosed(theE, theF) ? 1 : 0);
  if (aRep == nullptr)
  {
    return;
  }
  const BRep_GCurve* aGC = static_cast<const BRep_GCurve*>(aRep);
  double             f = 0., l = 0.;
  aGC->Range(f, l);
  gp_Pnt2d aP1, aP2;
  aRep->PCurve()->D0(f, aP1);
  aRep->PCurve()->D0(l, aP2);
  std::cout << " pc1[" << f << "," << l << "]=(" << aP1.X() << "," << aP1.Y() << ")->(" << aP2.X()
            << "," << aP2.Y() << ")";
  if (aRep->IsCurveOnClosedSurface())
  {
    occ::handle<Geom2d_Curve> aPc2 = aRep->PCurve2();
    aPc2->D0(f, aP1);
    aPc2->D0(l, aP2);
    std::cout << " pc2=(" << aP1.X() << "," << aP1.Y() << ")->(" << aP2.X() << "," << aP2.Y() << ")";
  }
}

} // namespace

int main()
{
  const TopoDS_Shape aBox =
    BRepPrimAPI_MakeBox(gp_Pnt(-1., -1., -1.), gp_Pnt(1., 1., 1.)).Shape();
  const TopoDS_Shape aCyl = BRepPrimAPI_MakeCylinder(0.4, 2.0).Shape();

  BOPAlgo_PaveFiller   aPF;
  TopTools_ListOfShape aArgs;
  aArgs.Append(aBox);
  aArgs.Append(aCyl);
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

  // The cylinder lateral face: the only U-closed cylinder in the DS.
  int         nFace = -1;
  TopoDS_Face aFace;
  for (int i = 0; i < aDS.NbShapes(); ++i)
  {
    if (aDS.Shape(i).IsNull() || aDS.Shape(i).ShapeType() != TopAbs_FACE)
    {
      continue;
    }
    TopoDS_Face         aF = TopoDS::Face(aDS.Shape(i));
    BRepAdaptor_Surface aS(aF);
    if (aS.GetType() == GeomAbs_Cylinder && aS.IsUClosed())
    {
      nFace = i;
      aFace = aF;
      break;
    }
  }
  if (nFace < 0)
  {
    std::cout << "cylinder lateral face not found\n";
    return 1;
  }
  std::cout << "FACE ds=" << nFace << " orient=" << (int)aFace.Orientation() << "\n";

  const BOPDS_FaceInfo& aFI = aDS.FaceInfo(nFace);
  std::cout << "  pavesIn=" << aFI.PaveBlocksIn().Extent()
            << " pavesOn=" << aFI.PaveBlocksOn().Extent()
            << " pavesSc=" << aFI.PaveBlocksSc().Extent() << "\n";

  // 1. Fill LE exactly as `BOPAlgo_Builder_2.cxx:359-494`.
  TopoDS_Face aFF = aFace;
  aFF.Orientation(TopAbs_FORWARD);

  NCollection_List<TopoDS_Shape>                        aLE;
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aMFence;
  bool                                                   isChecked = false, isUClosed = false,
                                                                      isVClosed = false;
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
        aLE.Append(aEE);
        aEE.Orientation(TopAbs_REVERSED);
        aLE.Append(aEE);
      }
      else
      {
        aLE.Append(aE);
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
    const NCollection_List<occ::handle<BOPDS_PaveBlock>>& aLPB           = aDS.PaveBlocks(nE);
    std::cout << "  BND ds=" << nE << " ori=" << (int)anOriE << " uclosed=" << (isUClosed ? 1 : 0)
              << " vclosed=" << (isVClosed ? 1 : 0) << " closed=" << (bIsClosed ? 1 : 0)
              << " degenerated=" << (bIsDegenerated ? 1 : 0) << " pbs=" << aLPB.Extent() << "\n";
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
        aLE.Append(aSp);
        continue;
      }
      if (anOriE == TopAbs_INTERNAL)
      {
        aSp.Orientation(TopAbs_FORWARD);
        aLE.Append(aSp);
        aSp.Orientation(TopAbs_REVERSED);
        aLE.Append(aSp);
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
              if (!BOPTools_AlgoTools3D::DoSplitSEAMOnFace(aE, TopoDS::Edge(aSp), aFF))
              {
                std::cout << "    DoSplitSEAMOnFace FAILED for split ds=" << aDS.Index(aSp)
                          << "\n";
              }
            }
          }
          aSp.Orientation(TopAbs_FORWARD);
          aLE.Append(aSp);
          aSp.Orientation(TopAbs_REVERSED);
          aLE.Append(aSp);
        }
        continue;
      }
      aSp.Orientation(anOriE);
      occ::handle<Message_Report> aRep = new Message_Report;
      if (BOPTools_AlgoTools::IsSplitToReverseWithWarn(aSp, aE, nullptr, aRep))
      {
        aSp.Reverse();
      }
      aLE.Append(aSp);
    }
  }
  // 1.2 In edges and 1.3 Section edges (`_2.cxx:467-494`): both doubled.
  for (int aPass = 0; aPass < 2; ++aPass)
  {
    const NCollection_IndexedMap<occ::handle<BOPDS_PaveBlock>>& aMPB =
      aPass == 0 ? aFI.PaveBlocksIn() : aFI.PaveBlocksSc();
    std::cout << ((aPass == 0) ? "  IN" : "  SC") << " pbs=" << aMPB.Extent() << "\n";
    for (int j = 1; j <= aMPB.Extent(); ++j)
    {
      TopoDS_Shape aSp = aDS.Shape(aMPB(j)->Edge());
      if (aSp.IsNull())
      {
        continue;
      }
      aSp.Orientation(TopAbs_FORWARD);
      aLE.Append(aSp);
      aSp.Orientation(TopAbs_REVERSED);
      aLE.Append(aSp);
    }
  }

  // 2. Per-occurrence table of `SplitBlock` (`_1.cxx:133-317`).
  occ::handle<IntTools_Context> aCtx = new IntTools_Context;
  BRepAdaptor_Surface&          aBAS = aCtx->SurfaceAdaptor(aFF);
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
  std::cout << "LE edges=" << aLE.Extent() << "\n";
  int aK = -1;
  for (NCollection_List<TopoDS_Shape>::Iterator anIt(aLE); anIt.More(); anIt.Next(), ++aK)
  {
    const TopoDS_Edge& aE = TopoDS::Edge(anIt.Value());
    double             f = 0., l = 0.;
    BRep_Tool::Range(aE, f, l);
    std::cout << " E" << aK << " ds=" << aDS.Index(aE) << " ori=" << (int)aE.Orientation()
              << " range=[" << f << "," << l << "] inside=" << (aMS.Contains(aE) ? 0 : 1)
              << " degenerated=" << (BRep_Tool::Degenerated(aE) ? 1 : 0);
    DumpRepresentations(aE, aFF);
    DumpCurve3d(aE);
    std::cout << "\n";
    int anIv = 0;
    for (TopoDS_Iterator aVI(aE); aVI.More(); aVI.Next(), ++anIv)
    {
      TopoDS_Vertex aV    = TopoDS::Vertex(aVI.Value());
      const bool    bIsIN = (aV.Orientation() == TopAbs_REVERSED);
      double        aT    = 0.;
      const gp_Pnt2d aP   = Coord2dOf(aV, aE, aFF, aT);
      std::cout << "    V" << anIv << " ori=" << (int)aV.Orientation() << " in=" << (bIsIN ? 1 : 0)
                << " t=" << aT << " uv=(" << aP.X() << "," << aP.Y() << ")"
                << " angle=" << Angle2D(aV, aE, aFF, aBAS, bIsIN, aCtx) << "\n";
    }
  }

  // 2b. Smart map exactly as `BOPAlgo_WireSplitter::SplitBlock` builds it
  //     (`_1.cxx:133-317`): keyed by vertex TShape (`TopTools_ShapeMapHasher`
  //     uses `IsSame`), one EdgeInfo per (edge occurrence, vertex occurrence).
  {
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
    struct Info
    {
      TopoDS_Edge      e;
      TopAbs_Orientation vor;
      bool             isIn;
      bool             isInside;
      double           angle;
    };
    NCollection_IndexedDataMap<TopoDS_Shape, int, TopTools_ShapeMapHasher>     aVSlot;
    std::vector<std::vector<Info>>                                             aVInfos;
    NCollection_DataMap<TopoDS_Shape, bool, TopTools_ShapeMapHasher> aVertMap;
    int aEdgeOrd = 0;
    for (NCollection_List<TopoDS_Shape>::Iterator anIt(aLE); anIt.More(); anIt.Next(), ++aEdgeOrd)
    {
      const TopoDS_Edge& aE = TopoDS::Edge(anIt.Value());
      bool               bIsClosed = BRep_Tool::Degenerated(aE) || BRep_Tool::IsClosed(aE, aFF);
      TopoDS_Iterator    aItS(aE);
      TopoDS_Shape       aV1;
      for (int i = 0; aItS.More(); aItS.Next(), ++i)
      {
        const TopoDS_Shape& aV   = aItS.Value();
        int                 slot = aVSlot.FindIndex(aV);
        if (!slot)
        {
          aVInfos.emplace_back();
          slot = aVSlot.Add(aV, int(aVInfos.size()) - 1);
        }
        Info anInfo;
        anInfo.e        = aE;
        anInfo.vor      = aV.Orientation();
        anInfo.isIn     = (aV.Orientation() == TopAbs_REVERSED);
        anInfo.isInside = !aMS.Contains(aE);
        anInfo.angle    = 0.;
        aVInfos[size_t(slot - 1)].push_back(anInfo);
        if (!i)
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
    for (int i = 1; i <= aVSlot.Extent(); ++i)
    {
      const TopoDS_Shape& aV = aVSlot.FindKey(i);
      const int           aL = aVSlot(i);
      // 3. Angles in mySmartMap (`_1.cxx:299-317`).
      for (auto& anInfo : aVInfos[size_t(aL)])
      {
        TopoDS_Vertex aVV = TopoDS::Vertex(aV);
        aVV.Orientation(anInfo.isIn ? TopAbs_REVERSED : TopAbs_FORWARD);
        anInfo.angle = Angle2D(aVV, anInfo.e, aFF, aBAS, anInfo.isIn, aCtx);
      }
    }
    for (int i = 1; i <= aVSlot.Extent(); ++i)
    {
      const TopoDS_Shape& aV = aVSlot.FindKey(i);
      const int           aL = aVSlot(i);
      std::cout << " SMV" << (i - 1) << " ds=" << aDS.Index(aV)
                << " closed=" << (aVertMap.Find(aV) ? 1 : 0) << " infos=" << aVInfos[size_t(aL)].size()
                << "\n";
      for (const auto& anInfo : aVInfos[size_t(aL)])
      {
        TopoDS_Vertex aVV = TopoDS::Vertex(aV);
        aVV.Orientation(TopAbs_FORWARD);
        double        aT = 0.;
        const gp_Pnt2d aP = Coord2dOf(aVV, anInfo.e, aFF, aT);
        std::cout << "   e ds=" << aDS.Index(anInfo.e) << " ori=" << (int)anInfo.e.Orientation()
                  << " vori=" << (int)anInfo.vor << " in=" << (anInfo.isIn ? 1 : 0)
                  << " inside=" << (anInfo.isInside ? 1 : 0) << " t=" << aT << " uv=(" << aP.X()
                  << "," << aP.Y() << ")"
                  << " angle=" << anInfo.angle;
        TopoDS_Vertex aVF = TopoDS::Vertex(aV);
        aVF.Orientation(TopAbs_FORWARD);
        TopoDS_Iterator aEI(anInfo.e);
        for (; aEI.More(); aEI.Next())
        {
          if (aEI.Value().Orientation() == TopAbs_FORWARD)
          {
            double         aTf = 0.;
            const gp_Pnt2d aPf = Coord2dOf(TopoDS::Vertex(aEI.Value()), anInfo.e, aFF, aTf);
            std::cout << " vf=(" << aPf.X() << "," << aPf.Y() << ")";
            break;
          }
        }
        std::cout << "\n";
      }
    }
    // The section edge's own vertices: same TShape (periodic closed edge)?
    for (NCollection_List<TopoDS_Shape>::Iterator anIt(aLE); anIt.More(); anIt.Next())
    {
      const TopoDS_Edge& aE = TopoDS::Edge(anIt.Value());
      if (aE.Orientation() != TopAbs_FORWARD || !aMS.Contains(aE))
      {
        continue;
      }
      TopoDS_Vertex aVf, aVl;
      TopExp::Vertices(aE, aVf, aVl);
      std::cout << " EDGE ds=" << aDS.Index(aE) << " V0.IsSame(V1)=" << (aVf.IsSame(aVl) ? 1 : 0)
                << " V0.IsEqual(V1)=" << (aVf.IsEqual(aVl) ? 1 : 0) << "\n";
    }
  }

  // 3. OCCT's own loop decomposition (`BOPAlgo_WireSplitter::Perform`).
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
    std::cout << "wire splitter errors\n";
    aWS.DumpErrors(std::cout);
  }
  int aNw = 0;
  for (NCollection_List<TopoDS_Shape>::Iterator anIt(aWES.Shapes()); anIt.More(); anIt.Next())
  {
    const TopoDS_Shape& aW = anIt.Value();
    std::cout << " LOOP " << aNw++ << ":";
    for (TopoDS_Iterator aEI(aW); aEI.More(); aEI.Next())
    {
      const TopoDS_Edge& aE = TopoDS::Edge(aEI.Value());
      std::cout << " " << aDS.Index(aE) << (aE.Orientation() == TopAbs_REVERSED ? "-" : "+");
    }
    std::cout << "\n";
  }
  std::cout << "loops=" << aNw << "\n";
  return 0;
}
