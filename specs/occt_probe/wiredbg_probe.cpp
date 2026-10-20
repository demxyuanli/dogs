// Wire ordering dumps: `--orderdbg`, `--wiredbg`, `--wiredbgi`.
#include "occt_probe_common.h"

// --- TEMP T-103 (`--orderdbg <defl> <angle> <x1> <y1> <z1> <x2> <y2> <z2>`) -----------
// Dumps, for every face whose bbox matches, the `BRepMesh_ShapeVisitor::addWire`
// ordering stage: `ShapeExtend_WireData` list + `ShapeAnalysis_Wire::CheckOrder`
// (`ShapeAnalysis_WireOrder::Perform`) result, i.e. the exact `Ordered(i)` chain
// the port's `WireOrder` must reproduce.
int runOrderDbg(const TopoDS_Shape& theShape, const double* theBBox)
{
  int aFound = 0;
  TopExp_Explorer aEx(theShape, TopAbs_FACE);
  for (; aEx.More(); aEx.Next())
  {
    const TopoDS_Face& aFace = TopoDS::Face(aEx.Current());
    Bnd_Box            aBox;
    BRepBndLib::Add(aFace, aBox, false);
    if (aBox.IsVoid())
    {
      continue;
    }
    double v[6] = {0, 0, 0, 0, 0, 0};
    aBox.Get(v[0], v[1], v[2], v[3], v[4], v[5]);
    if (!(std::abs(v[0] - theBBox[0]) < 1e-3 && std::abs(v[1] - theBBox[1]) < 1e-3
          && std::abs(v[2] - theBBox[2]) < 1e-3 && std::abs(v[3] - theBBox[3]) < 1e-3
          && std::abs(v[4] - theBBox[4]) < 1e-3 && std::abs(v[5] - theBBox[5]) < 1e-3))
    {
      continue;
    }
    ++aFound;
    // TEMP: print the occurrence face orientation, and repeat the whole dump for
    // the oriented-FORWARD face the real `BRepMesh_ShapeVisitor::Visit` uses.
    std::printf("ODBG faceori=%d\n", (int)aFace.Orientation());
    const TopoDS_Face aFaceFwd = TopoDS::Face(aFace.Oriented(TopAbs_FORWARD));
    for (int aPass = 0; aPass < 2; ++aPass)
    {
    const TopoDS_Face& aSrc = (aPass == 0) ? aFace : aFaceFwd;
    std::printf("ODBG pass=%d\n", aPass);
    int aWireIt = 0;
    for (TopExp_Explorer aW(aSrc, TopAbs_WIRE); aW.More(); aW.Next(), ++aWireIt)
    {
      const TopoDS_Wire& aWire     = TopoDS::Wire(aW.Current());
      occ::handle<ShapeExtend_WireData> aWireData = new ShapeExtend_WireData(aWire, true, false);
      std::printf("ODBG wire=%d stored=%d ori=%d\n", aWireIt, aWireData->NbEdges(),
                  (int)aWire.Orientation());
      for (int i = 1; i <= aWireData->NbEdges(); ++i)
      {
        const TopoDS_Edge&   anEdge = aWireData->Edge(i);
        BRepAdaptor_Curve    aC3d(TopoDS::Edge(anEdge.Oriented(TopAbs_FORWARD)));
        double               f = 0., l = 0.;
        const occ::handle<Geom2d_Curve> aC2d = BRep_Tool::CurveOnSurface(
          TopoDS::Edge(anEdge.Oriented(TopAbs_FORWARD)), aFace, f, l);
        gp_Pnt2d aB(0., 0.), aE(0., 0.);
        if (!aC2d.IsNull())
        {
          const bool aRev = (anEdge.Orientation() == TopAbs_REVERSED);
          aB = aC2d->Value(aRev ? l : f);
          aE = aC2d->Value(aRev ? f : l);
        }
        const gp_Pnt aP0 = aC3d.Value(aC3d.FirstParameter());
        const gp_Pnt aP1 = aC3d.Value(aC3d.LastParameter());
        const char*  aKind = aC2d.IsNull() ? "null" : aC2d->DynamicType()->Name();
        const bool   aPer  = !aC2d.IsNull() && aC2d->IsPeriodic();
        std::printf("ODBG   i=%d edge=%p ori=%d closed=%d deg=%d uv=(%.9f,%.9f)->(%.9f,%.9f) "
                    "p3d=(%.6f,%.6f,%.6f)->(%.6f,%.6f,%.6f)\n",
                    i,
                    (void*)anEdge.TShape().get(),
                    (int)anEdge.Orientation(),
                    BRep_Tool::IsClosed(anEdge, aFace) ? 1 : 0,
                    BRep_Tool::Degenerated(anEdge) ? 1 : 0,
                    aB.X(),
                    aB.Y(),
                    aE.X(),
                    aE.Y(),
                    aP0.X(),
                    aP0.Y(),
                    aP0.Z(),
                    aP1.X(),
                    aP1.Y(),
                    aP1.Z());
        std::printf("ODBG     c2d=%s per=%d range2d=(%.9f,%.9f) range3d=(%.9f,%.9f) sameParam=%d "
                    "sameRange=%d v(f)=(%.9f,%.9f) v(l)=(%.9f,%.9f)\n",
                    aKind,
                    aPer ? 1 : 0,
                    f,
                    l,
                    aC3d.FirstParameter(),
                    aC3d.LastParameter(),
                    BRep_Tool::SameParameter(TopoDS::Edge(anEdge.Oriented(TopAbs_FORWARD))) ? 1 : 0,
                    BRep_Tool::SameRange(TopoDS::Edge(anEdge.Oriented(TopAbs_FORWARD))) ? 1 : 0,
                    aC2d.IsNull() ? 0.0 : aC2d->Value(f).X(),
                    aC2d.IsNull() ? 0.0 : aC2d->Value(f).Y(),
                    aC2d.IsNull() ? 0.0 : aC2d->Value(l).X(),
                    aC2d.IsNull() ? 0.0 : aC2d->Value(l).Y());
      }
      ShapeAnalysis_Wire       aWireTool(aWireData, aFace, Precision::Confusion());
      ShapeAnalysis_WireOrder  aOrderTool;
      aWireTool.CheckOrder(aOrderTool, true, false);
      std::printf("ODBG wire=%d status: fail=%d done3=%d nb=%d ordered=[",
                  aWireIt,
                  aWireTool.LastCheckStatus(ShapeExtend_FAIL) ? 1 : 0,
                  aWireTool.LastCheckStatus(ShapeExtend_DONE3) ? 1 : 0,
                  aOrderTool.NbEdges());
      for (int i = 1; i <= aOrderTool.NbEdges(); ++i)
      {
        std::printf("%s%d", i > 1 ? " " : "", aOrderTool.Ordered(i));
      }
      std::printf("]\n");
    }
    }
    break;
  }
  std::printf("ODBG found=%d\n", aFound);
  return 0;
}

// --- TEMP T-103 (`--wiredbg <defl> <angle> <x1> <y1> <z1> <x2> <y2> <z2>`) ------------
// Ground truth for the port's `CWDBG`/`CBWIN`. Replays the real `IMeshTools_MeshBuilder`
// prefix (`BuildModel` -> `DiscretizeEdges` -> `HealModel`) and dumps, for every wire
// slot of the target face, the `IMeshData_Wire` order: the slot orientation, the
// pcurve orientation / `IsForward`, `SameParam` / `SameRange`, both parameter ranges
// and the full 3D + 2D point arrays in `GetPoint(i)` order (the arrays the classifier
// polygon and the Delaunay frontier are built from).
int runWireDbg(const TopoDS_Shape& theShape,
                      const double        theDeflection,
                      const double        theAngle,
                      const double*       theBBox,
                      const int           theFaceIdx)
{
  IMeshTools_Parameters aParams;
  aParams.Deflection    = theDeflection;
  aParams.Angle         = theAngle;
  aParams.InParallel    = false;
  aParams.Relative      = false;
  aParams.MinSize       = (std::max)(IMeshTools_Parameters::RelMinSize() * theDeflection,
                                     Precision::Confusion());
  aParams.AdjustMinSize = false;

  occ::handle<BRepMesh_ModelBuilder> aBuilder = new BRepMesh_ModelBuilder;
  const occ::handle<IMeshData_Model> aModel   = aBuilder->Perform(theShape, aParams);
  if (aModel.IsNull())
  {
    std::printf("WDBG model-null\n");
    return 0;
  }
  // `IMeshTools_MeshBuilder::Perform` (`cxx:52-56`): DiscretizeEdges runs BEFORE
  // HealModel, so the pre-heal arrays below are the EdgeDiscret output.
  occ::handle<BRepMesh_EdgeDiscret> anEdgeDiscret = new BRepMesh_EdgeDiscret;
  anEdgeDiscret->Perform(aModel, aParams, Message_ProgressRange());
  occ::handle<BRepMesh_ModelHealer> aHealer = new BRepMesh_ModelHealer;
  aHealer->Perform(aModel, aParams, Message_ProgressRange());

  for (int f = 0; f < aModel->FacesNb(); ++f)
  {
    const IMeshData::IFaceHandle& aDFace = aModel->GetFace(f);
    const TopoDS_Face&            aFace  = aDFace->GetFace();
    Bnd_Box                       aBox;
    BRepBndLib::Add(aFace, aBox, false);
    if (aBox.IsVoid())
    {
      continue;
    }
    double v[6] = {0, 0, 0, 0, 0, 0};
    aBox.Get(v[0], v[1], v[2], v[3], v[4], v[5]);
    // Selection: by `IMeshData_Model::GetFace` index when the caller gave one
    // (the port's own bbox is tolerance-inflated, so the six-coordinate form
    // cannot pair the two sides), otherwise by the exact geometry box.
    if (theFaceIdx >= 0)
    {
      if (f != theFaceIdx)
      {
        continue;
      }
    }
    else if (!(std::abs(v[0] - theBBox[0]) < 1e-3 && std::abs(v[1] - theBBox[1]) < 1e-3
               && std::abs(v[2] - theBBox[2]) < 1e-3 && std::abs(v[3] - theBBox[3]) < 1e-3
               && std::abs(v[4] - theBBox[4]) < 1e-3 && std::abs(v[5] - theBBox[5]) < 1e-3))
    {
      continue;
    }
    std::printf("WDBG face=%d bbox=(%.9f,%.9f,%.9f)-(%.9f,%.9f,%.9f)\n", f, v[0], v[1], v[2], v[3],
                v[4], v[5]);
    std::printf("WDBG face=%d surfType=%d wires=%d\n", f, (int)aDFace->GetSurface()->GetType(),
                aDFace->WiresNb());
    for (int w = 0; w < aDFace->WiresNb(); ++w)
    {
      const IMeshData::IWireHandle& aDWire = aDFace->GetWire(w);
      std::printf("WDBG  w=%d edges=%d\n", w, aDWire->EdgesNb());
      for (int e = 0; e < aDWire->EdgesNb(); ++e)
      {
        const IMeshData::IEdgeHandle& aDEdge  = aDWire->GetEdge(e);
        const TopAbs_Orientation      aSlotOri = aDWire->GetEdgeOrientation(e);
        const IMeshData::IPCurveHandle& aPCurve = aDEdge->GetPCurve(aDFace.get(), aSlotOri);
        const TopoDS_Edge&              aEdge   = aDEdge->GetEdge();
        double                          r3f = 0., r3l = 0.;
        BRep_Tool::Range(TopoDS::Edge(aEdge.Oriented(TopAbs_FORWARD)), r3f, r3l);
        std::printf("WDBG   e=%d ptr=%p slotOri=%d pcsNb=%d pcOri=%d fwd=%d sameParam=%d "
                    "sameRange=%d deg=%d range3d=(%.9f,%.9f) nb3d=%d nb2d=%d\n",
                    e,
                    (void*)aEdge.TShape().get(),
                    (int)aSlotOri,
                    aDEdge->PCurvesNb(),
                    (int)aPCurve->GetOrientation(),
                    aPCurve->IsForward() ? 1 : 0,
                    aDEdge->GetSameParam() ? 1 : 0,
                    aDEdge->GetSameRange() ? 1 : 0,
                    aDEdge->GetDegenerated() ? 1 : 0,
                    r3f,
                    r3l,
                    aDEdge->GetCurve()->ParametersNb(),
                    aPCurve->ParametersNb());
        std::printf("WDBG   p2d=[");
        for (int k = 0; k < aPCurve->ParametersNb(); ++k)
        {
          std::printf("%s(%.17g,%.17g)", k ? " " : "", aPCurve->GetPoint(k).X(),
                      aPCurve->GetPoint(k).Y());
        }
        std::printf("]\n");
        std::printf("WDBG   p3d=[");
        for (int k = 0; k < aDEdge->GetCurve()->ParametersNb(); ++k)
        {
          const gp_Pnt& p = aDEdge->GetCurve()->GetPoint(k);
          std::printf("%s(%.6f,%.6f,%.6f)", k ? " " : "", p.X(), p.Y(), p.Z());
        }
        std::printf("]\n");
      }
    }
  }
  return 0;
}
