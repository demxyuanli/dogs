// Per-stage mesh replays: `--torusdbg`, `--healerdbg`, `--ctxdbg`, `--rangereplay`,
// `--algospy`.
#include "occt_probe_common.h"

// --- TEMP T-103 (`--torusdbg <defl> <angle>`) --------------------------------------
// Replays OCCT's own `BRepMesh_TorusRangeSplitter` on every torus face of the
// discrete model built with the given parameters, printing the exact interior
// node count it produces. `_dbg/ZZ_TorusRangeSplitter.cxx` supplies both the
// byte-identical splitter copy and `ZZTorDbgFace`.
void ZZTorDbgFace(const IMeshData::IFaceHandle& theDFace,
                  const IMeshTools_Parameters&  theParameters,
                  const int                     theFaceIdx);

int runTorusDbg(const TopoDS_Shape& theShape, const double theDeflection,
                       const double theAngle)
{
  IMeshTools_Parameters aParams;
  aParams.Deflection    = theDeflection;
  aParams.Angle         = theAngle;
  aParams.InParallel    = false;
  aParams.Relative      = false;
  // Mirror `BRepMesh_IncrementalMesh::initParameters` (`hxx:91-96`): the real
  // mesher runs with `MinSize = max(RelMinSize * min(defl, deflInterior), Confusion)`.
  // The earlier Confusion setting made `BRepMesh_EdgeDiscret` emit finer pcurves
  // than the real run, so the splitter inputs were not comparable.
  aParams.MinSize       = (std::max)(IMeshTools_Parameters::RelMinSize() * theDeflection,
                                     Precision::Confusion());
  aParams.AdjustMinSize = false;

  occ::handle<BRepMesh_ModelBuilder> aBuilder = new BRepMesh_ModelBuilder;
  const occ::handle<IMeshData_Model> aModel   = aBuilder->Perform(theShape, aParams);
  if (aModel.IsNull())
  {
    std::cout << "ZZTORDBG model-null\n";
    return 0;
  }
  // The builder only enumerates the faces; the pcurves that
  // `collectWirePoints` walks are filled by `BRepMesh_EdgeDiscret`
  // (`BRepMesh_ModelBuilder` + `BRepMesh_ModelHealer` + `BRepMesh_EdgeDiscret`
  // is the `IMeshTools_Context` chain). Without this every replay reports
  // `ParametersNb() == 0` and `AdjustRange` leaves the splitter invalid.
  occ::handle<BRepMesh_EdgeDiscret> anEdgeDiscret = new BRepMesh_EdgeDiscret;
  anEdgeDiscret->Perform(aModel, aParams, Message_ProgressRange());
  for (int f = 0; f < aModel->FacesNb(); ++f)
  {
    const IMeshData::IFaceHandle& aDFace = aModel->GetFace(f);
    if (aDFace->GetSurface()->GetType() != GeomAbs_Torus)
    {
      continue;
    }
    ZZTorDbgFace(aDFace, aParams, f);
  }
  return 0;
}

// --- TEMP T-103 (`--healerdbg <defl> <angle> <x1> <y1> <z1> <x2> <y2> <z2>`) ---------
// Dumps the per-wire-edge pcurve state of every discrete face whose bbox matches
// the given target, BEFORE and AFTER `BRepMesh_ModelHealer::Perform`. This is the
// oracle for the port's `ModelHealer::fix_face_boundaries` snap: the pre-heal
// columns must be bit-equal to the port's EdgeDiscret output, and the post-heal
// columns show what OCCT's own `connectClosestPoints`/`adjustSamePoints` write.
static void dumpFacePcurves(const IMeshData::IFaceHandle& theDFace,
                            const int                     theFaceIdx,
                            const char*                   theTag)
{
  for (int aWireIt = 0; aWireIt < theDFace->WiresNb(); ++aWireIt)
  {
    const IMeshData::IWireHandle& aDWire = theDFace->GetWire(aWireIt);
    for (int aEdgeIt = 0; aEdgeIt < aDWire->EdgesNb(); ++aEdgeIt)
    {
      const IMeshData::IEdgeHandle&   aDEdge  = aDWire->GetEdge(aEdgeIt);
      const IMeshData::IPCurveHandle& aPCurve = aDEdge->GetPCurve(
        theDFace.get(), aDWire->GetEdgeOrientation(aEdgeIt));
      const int aNb = aPCurve->ParametersNb();
      std::printf("HDBG %s face=%d w=%d e=%d nb=%d fwd=%d ori=%d same=%d free=%d int=%d "
                  "first=(%.9f,%.9f) last=(%.9f,%.9f)\n",
                  theTag,
                  theFaceIdx,
                  aWireIt,
                  aEdgeIt,
                  aNb,
                  aPCurve->IsForward() ? 1 : 0,
                  (int)aPCurve->GetOrientation(),
                  aDEdge->GetSameParam() ? 1 : 0,
                  aDEdge->IsFree() ? 1 : 0,
                  aPCurve->IsInternal() ? 1 : 0,
                  aNb > 0 ? aPCurve->GetPoint(0).X() : 0.0,
                  aNb > 0 ? aPCurve->GetPoint(0).Y() : 0.0,
                  aNb > 0 ? aPCurve->GetPoint(aNb - 1).X() : 0.0,
                  aNb > 0 ? aPCurve->GetPoint(aNb - 1).Y() : 0.0);
      if (aNb > 0 && aNb <= 64)
      {
        std::printf("HDBG %s face=%d e=%d pts=[", theTag, theFaceIdx, aEdgeIt);
        for (int k = 0; k < aNb; ++k)
        {
          std::printf("%s(%.6f,%.6f)", k ? " " : "", aPCurve->GetPoint(k).X(), aPCurve->GetPoint(k).Y());
        }
        std::printf("]\n");
      }
    }
  }
}

int runHealerDbg(const TopoDS_Shape& theShape,
                        const double        theDeflection,
                        const double        theAngle,
                        const double*       theBBox)
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
    std::cout << "HDBG model-null" << std::endl;
    return 0;
  }
  occ::handle<BRepMesh_EdgeDiscret> anEdgeDiscret = new BRepMesh_EdgeDiscret;
  anEdgeDiscret->Perform(aModel, aParams, Message_ProgressRange());

  std::vector<int> aHits;
  for (int f = 0; f < aModel->FacesNb(); ++f)
  {
    const TopoDS_Face& aFace = aModel->GetFace(f)->GetFace();
    Bnd_Box            aBox;
    BRepBndLib::Add(aFace, aBox, false);
    if (aBox.IsVoid())
    {
      continue;
    }
    double v[6] = {0, 0, 0, 0, 0, 0};
    aBox.Get(v[0], v[1], v[2], v[3], v[4], v[5]);
    if (std::abs(v[0] - theBBox[0]) < 1e-3 && std::abs(v[1] - theBBox[1]) < 1e-3
        && std::abs(v[2] - theBBox[2]) < 1e-3 && std::abs(v[3] - theBBox[3]) < 1e-3
        && std::abs(v[4] - theBBox[4]) < 1e-3 && std::abs(v[5] - theBBox[5]) < 1e-3)
    {
      aHits.push_back(f);
    }
  }
  std::cout << "HDBG hits=" << aHits.size() << std::endl;
  for (size_t i = 0; i < aHits.size(); ++i)
  {
    dumpFacePcurves(aModel->GetFace(aHits[i]), aHits[i], "pre");
  }

  occ::handle<BRepMesh_ModelHealer> aHealer = new BRepMesh_ModelHealer;
  aHealer->Perform(aModel, aParams, Message_ProgressRange());

  for (size_t i = 0; i < aHits.size(); ++i)
  {
    dumpFacePcurves(aModel->GetFace(aHits[i]), aHits[i], "post");
  }

  // TEMP T-122: the face flags the healer left behind (`BRepMesh_ModelHealer::
  // Perform`, `cxx:156-166`: a FaceChecker failure marks SelfIntersectingWire +
  // Failure, and `BRepMesh_FaceDiscret::process` (`cxx:102`) then skips the face
  // entirely, so it ends up with no triangulation at all).
  for (size_t i = 0; i < aHits.size(); ++i)
  {
    const IMeshData::IFaceHandle& aDFace = aModel->GetFace(aHits[i]);
    int aOpenWire = 0;
    for (int w = 0; w < aDFace->WiresNb(); ++w)
    {
      if (aDFace->GetWire(w)->IsSet(IMeshData_OpenWire))
      {
        aOpenWire = 1;
      }
    }
    std::cout << "HDBGFLAG f=" << aHits[i] << " wires=" << aDFace->WiresNb()
              << " fail=" << (aDFace->IsSet(IMeshData_Failure) ? 1 : 0)
              << " selfint=" << (aDFace->IsSet(IMeshData_SelfIntersectingWire) ? 1 : 0)
              << " open=" << aOpenWire
              << " outdated=" << (aDFace->IsSet(IMeshData_Outdated) ? 1 : 0)
              << " unoriented=" << (aDFace->IsSet(IMeshData_UnorientedWire) ? 1 : 0)
              << " defl=" << std::setprecision(17) << aDFace->GetDeflection() << std::endl;
  }
  return 0;
}

// --- TEMP T-122 (`--ctxdbg <defl> <angle> <x1> <y1> <z1> <x2> <y2> <z2>`) -------------
// Runs the real `IMeshTools_MeshBuilder` stage sequence (`IMeshTools_MeshBuilder.cxx:
// 52-62`: BuildModel -> DiscretizeEdges -> HealModel -> PreProcessModel ->
// DiscretizeFaces -> PostProcessModel) on the bbox-selected faces and dumps their
// `IMeshData_*` status flags plus the resulting `Poly_Triangulation` size after EVERY
// stage. The stage whose dump flips a flag - or whose triangulation is still missing -
// is the stage that leaves the face unmeshed, so the behaviour gets ported branch by
// branch instead of guessed from the final state. The parameter block mirrors
// `BRepMesh_IncrementalMesh::initParameters` (`BRepMesh_IncrementalMesh.hxx:79-108`);
// the per-stage totals as well as the final `nodes=`/`tris=` must match
// `--mesh <defl> <angle>` on the same shape, which is this replay's self-check.
namespace
{
class ZZStageCtx : public BRepMesh_Context
{
public:
  ZZStageCtx(const double* theBBox)
  {
    for (int i = 0; i < 6; ++i)
    {
      myBBox[i] = theBBox[i];
    }
  }

  bool BuildModel() override
  {
    const bool aRes = BRepMesh_Context::BuildModel();
    findFaces();
    dumpStage("build", aRes);
    return aRes;
  }

  bool DiscretizeEdges() override
  {
    const bool aRes = BRepMesh_Context::DiscretizeEdges();
    dumpStage("disc", aRes);
    return aRes;
  }

  bool HealModel() override
  {
    const bool aRes = BRepMesh_Context::HealModel();
    dumpStage("heal", aRes);
    return aRes;
  }

  bool PreProcessModel() override
  {
    const bool aRes = BRepMesh_Context::PreProcessModel();
    dumpStage("preproc", aRes);
    return aRes;
  }

  bool DiscretizeFaces(const Message_ProgressRange& theRange) override
  {
    const bool aRes = BRepMesh_Context::DiscretizeFaces(theRange);
    dumpStage("faces", aRes);
    return aRes;
  }

  bool PostProcessModel() override
  {
    const bool aRes = BRepMesh_Context::PostProcessModel();
    dumpStage("postproc", aRes);
    return aRes;
  }

private:
  void findFaces()
  {
    const occ::handle<IMeshData_Model>& aModel = GetModel();
    if (aModel.IsNull())
    {
      return;
    }
    for (int f = 0; f < aModel->FacesNb(); ++f)
    {
      Bnd_Box aBox;
      BRepBndLib::Add(aModel->GetFace(f)->GetFace(), aBox, false);
      if (aBox.IsVoid())
      {
        continue;
      }
      double v[6] = {0, 0, 0, 0, 0, 0};
      aBox.Get(v[0], v[1], v[2], v[3], v[4], v[5]);
      bool aHit = true;
      for (int k = 0; k < 6; ++k)
      {
        if (std::abs(v[k] - myBBox[k]) >= 1e-3)
        {
          aHit = false;
        }
      }
      if (aHit)
      {
        myHits.push_back(f);
      }
    }
  }

  void dumpStage(const char* theTag, const bool theRes) const
  {
    const occ::handle<IMeshData_Model>& aModel = GetModel();
    std::cout << "CTXDBG " << theTag << " res=" << (theRes ? 1 : 0)
              << " model_null=" << (aModel.IsNull() ? 1 : 0);
    if (aModel.IsNull())
    {
      std::cout << std::endl;
      return;
    }
    int aNodes = 0, aTris = 0;
    for (int f = 0; f < aModel->FacesNb(); ++f)
    {
      TopLoc_Location                        aLoc;
      const occ::handle<Poly_Triangulation>& aTri =
        BRep_Tool::Triangulation(aModel->GetFace(f)->GetFace(), aLoc);
      if (!aTri.IsNull())
      {
        aNodes += aTri->NbNodes();
        aTris += aTri->NbTriangles();
      }
    }
    std::cout << " faces=" << aModel->FacesNb() << " nodes=" << aNodes << " tris=" << aTris
              << std::endl;
    for (size_t i = 0; i < myHits.size(); ++i)
    {
      const IMeshData::IFaceHandle& aDFace = aModel->GetFace(myHits[i]);
      int                           aOpenWire = 0;
      for (int w = 0; w < aDFace->WiresNb(); ++w)
      {
        if (aDFace->GetWire(w)->IsSet(IMeshData_OpenWire))
        {
          aOpenWire = 1;
        }
      }
      TopLoc_Location                        aLoc;
      const occ::handle<Poly_Triangulation>& aTri =
        BRep_Tool::Triangulation(aDFace->GetFace(), aLoc);
      std::cout << "CTXFLAG " << theTag << " f=" << myHits[i]
                << " stype=" << static_cast<int>(aDFace->GetSurface()->GetType())
                << " fail=" << (aDFace->IsSet(IMeshData_Failure) ? 1 : 0)
                << " reused=" << (aDFace->IsSet(IMeshData_Reused) ? 1 : 0)
                << " remesh=" << (aDFace->IsSet(IMeshData_ReMesh) ? 1 : 0)
                << " selfint=" << (aDFace->IsSet(IMeshData_SelfIntersectingWire) ? 1 : 0)
                << " unoriented=" << (aDFace->IsSet(IMeshData_UnorientedWire) ? 1 : 0)
                << " toofew=" << (aDFace->IsSet(IMeshData_TooFewPoints) ? 1 : 0)
                << " outdated=" << (aDFace->IsSet(IMeshData_Outdated) ? 1 : 0)
                << " open=" << aOpenWire << " tri=" << (aTri.IsNull() ? -1 : aTri->NbNodes())
                << "/" << (aTri.IsNull() ? -1 : aTri->NbTriangles()) << std::endl;
      // TEMP T-123: the wire statuses are what `initDataStructure` branches on
      // (`BRepMesh_BaseMeshAlgo.cxx:80`, `BRepMesh_NodeInsertionMeshAlgo.hxx:70`), and
      // they are invisible in the face-level flags above.
      for (int w = 0; w < aDFace->WiresNb(); ++w)
      {
        const IMeshData::IWireHandle& aDWire = aDFace->GetWire(w);
        std::cout << "CTXWIRE " << theTag << " f=" << myHits[i] << " w=" << w
                  << " edges=" << aDWire->EdgesNb()
                  << " selfint=" << (aDWire->IsSet(IMeshData_SelfIntersectingWire) ? 1 : 0)
                  << " open=" << (aDWire->IsSet(IMeshData_OpenWire) ? 1 : 0) << std::endl;
      }
    }
  }

  double           myBBox[6];
  std::vector<int> myHits;
};
} // namespace

int runCtxDbg(const TopoDS_Shape& theShape,
                     const double        theDeflection,
                     const double        theAngle,
                     const double*       theBBox)
{
  IMeshTools_Parameters aParams;
  aParams.Deflection         = theDeflection;
  aParams.DeflectionInterior = theDeflection;
  aParams.Angle              = theAngle;
  aParams.AngleInterior      = 2.0 * theAngle;
  aParams.MinSize =
    (std::max)(IMeshTools_Parameters::RelMinSize() * theDeflection, Precision::Confusion());
  aParams.AdjustMinSize = false;
  aParams.InParallel    = false;
  aParams.Relative      = false;
  aParams.CleanModel    = false;

  std::cout << "CTXDBG params defl=" << std::setprecision(17) << aParams.Deflection
            << " deflint=" << aParams.DeflectionInterior << " angle=" << aParams.Angle
            << " angint=" << aParams.AngleInterior << " minsize=" << aParams.MinSize
            << " ivm=" << (aParams.InternalVerticesMode ? 1 : 0)
            << " ecsda=" << (aParams.EnableControlSurfaceDeflectionAllSurfaces ? 1 : 0)
            << std::endl;

  occ::handle<ZZStageCtx> aCtx = new ZZStageCtx(theBBox);
  aCtx->SetShape(theShape);
  aCtx->ChangeParameters() = aParams;

  if (!aCtx->BuildModel() || !aCtx->DiscretizeEdges() || !aCtx->HealModel()
      || !aCtx->PreProcessModel() || !aCtx->DiscretizeFaces(Message_ProgressRange()))
  {
    return 0;
  }
  aCtx->PostProcessModel();
  return 0;
}

// --- TEMP T-122 (`--rangereplay <faceIdx> <defl> <angle>`) ---------------------------
// Replays the splitter half of `BRepMesh_NodeInsertionMeshAlgo::initDataStructure` for
// one model face: `Reset` -> `AddPoint` for every point of every pcurve of wire 0 in the
// wire's edge orientation (`BRepMesh_NodeInsertionMeshAlgo.hxx:143-183`, i.e. the same
// `collectWirePoints` loop the `SeamEdgeAmplifier` uses to size the cone seam step) ->
// `AdjustRange` (`BRepMesh_DefaultRangeSplitter.cxx:45-81`) -> `IsValid` (`cxx:85-88`).
// The splitter class follows `BRepMesh_MeshAlgoFactory::GetAlgo`
// (`BRepMesh_MeshAlgoFactory.cxx:60-118`) for the face's surface type, so `valid=0` here
// is exactly the branch that makes OCCT drop a face: `initDataStructure`
// (`BRepMesh_NodeInsertionMeshAlgo.hxx:79-83`) sets `IMeshData_Failure` and returns
// false, `BRepMesh_BaseMeshAlgo::Perform` (`cxx:47-61`) then never runs
// `generateMesh`/`commitSurfaceTriangulation`, and the face keeps a null
// `Poly_Triangulation`. The pipeline stops after `PreProcessModel` because the
// `SeamEdgeAmplifier` inserts the seam points before `DiscretizeFaces` runs, so those
// points have to be in the replayed pcurves.
template <class Splitter>
static void replayRangeOn(const IMeshData::IFaceHandle& theDFace,
                          const IMeshTools_Parameters&  theParams,
                          const char*                   theName)
{
  Splitter aSplitter;
  aSplitter.Reset(theDFace, theParams);

  const IMeshData::IWireHandle& aDWire = theDFace->GetWire(0);
  double aMinU = 1e100, aMaxU = -1e100, aMinV = 1e100, aMaxV = -1e100;
  int    aPoints = 0;
  for (int e = 0; e < aDWire->EdgesNb(); ++e)
  {
    const IMeshData::IEdgeHandle&   aDEdge = aDWire->GetEdge(e);
    const IMeshData::IPCurveHandle& aPCurve =
      aDEdge->GetPCurve(theDFace.get(), aDWire->GetEdgeOrientation(e));
    for (int p = 0; p < aPCurve->ParametersNb(); ++p)
    {
      const gp_Pnt2d& aPnt = aPCurve->GetPoint(p);
      aSplitter.AddPoint(aPnt);
      aMinU = (std::min)(aMinU, aPnt.X());
      aMaxU = (std::max)(aMaxU, aPnt.X());
      aMinV = (std::min)(aMinV, aPnt.Y());
      aMaxV = (std::max)(aMaxV, aPnt.Y());
      ++aPoints;
    }
  }
  aSplitter.AdjustRange();

  const occ::handle<BRepAdaptor_Surface>& aSurf  = theDFace->GetSurface();
  const std::pair<double, double>&        aRU    = aSplitter.GetRangeU();
  const std::pair<double, double>&        aRV    = aSplitter.GetRangeV();
  const std::pair<double, double>&        aDelta = aSplitter.GetDelta();
  const std::pair<double, double>&        aTol   = aSplitter.GetToleranceUV();
  std::printf("RANGEREPLAY %s face_type=%d valid=%d pts=%d "
              "disc_u=[%.17g,%.17g] disc_v=[%.17g,%.17g] "
              "ru=[%.17g,%.17g] rv=[%.17g,%.17g] "
              "geom_u=[%.17g,%.17g] pu=%d geom_v=[%.17g,%.17g] pv=%d "
              "delta=[%.17g,%.17g] tol=[%.17g,%.17g]\n",
              theName,
              static_cast<int>(aSurf->GetType()),
              aSplitter.IsValid() ? 1 : 0,
              aPoints,
              aMinU,
              aMaxU,
              aMinV,
              aMaxV,
              aRU.first,
              aRU.second,
              aRV.first,
              aRV.second,
              aSurf->FirstUParameter(),
              aSurf->LastUParameter(),
              aSurf->IsUPeriodic() ? 1 : 0,
              aSurf->FirstVParameter(),
              aSurf->LastVParameter(),
              aSurf->IsVPeriodic() ? 1 : 0,
              aDelta.first,
              aDelta.second,
              aTol.first,
              aTol.second);
}

int runRangeReplay(const TopoDS_Shape& theShape,
                          const int           theFaceIdx,
                          const double        theDeflection,
                          const double        theAngle)
{
  IMeshTools_Parameters aParams;
  aParams.Deflection         = theDeflection;
  aParams.DeflectionInterior = theDeflection;
  aParams.Angle              = theAngle;
  aParams.AngleInterior      = 2.0 * theAngle;
  aParams.MinSize =
    (std::max)(IMeshTools_Parameters::RelMinSize() * theDeflection, Precision::Confusion());
  aParams.AdjustMinSize = false;
  aParams.InParallel    = false;
  aParams.Relative      = false;
  aParams.CleanModel    = false;

  occ::handle<BRepMesh_Context> aCtx = new BRepMesh_Context;
  aCtx->SetShape(theShape);
  aCtx->ChangeParameters() = aParams;
  if (!aCtx->BuildModel() || !aCtx->DiscretizeEdges() || !aCtx->HealModel()
      || !aCtx->PreProcessModel())
  {
    return 0;
  }
  const occ::handle<IMeshData_Model>& aModel = aCtx->GetModel();
  if (aModel.IsNull() || theFaceIdx < 0 || theFaceIdx >= aModel->FacesNb())
  {
    std::cout << "RANGEREPLAY model-null-or-range\n";
    return 0;
  }
  const IMeshData::IFaceHandle& aDFace = aModel->GetFace(theFaceIdx);
  switch (aDFace->GetSurface()->GetType())
  {
    case GeomAbs_Cone:
      replayRangeOn<BRepMesh_ConeRangeSplitter>(aDFace, aParams, "cone");
      break;
    case GeomAbs_Cylinder:
      replayRangeOn<BRepMesh_CylinderRangeSplitter>(aDFace, aParams, "cylinder");
      break;
    case GeomAbs_Sphere:
      replayRangeOn<BRepMesh_SphereRangeSplitter>(aDFace, aParams, "sphere");
      break;
    case GeomAbs_Torus:
      replayRangeOn<BRepMesh_TorusRangeSplitter>(aDFace, aParams, "torus");
      break;
    case GeomAbs_BezierSurface:
    case GeomAbs_BSplineSurface:
      replayRangeOn<BRepMesh_NURBSRangeSplitter>(aDFace, aParams, "nurbs");
      break;
    default:
      replayRangeOn<BRepMesh_DefaultRangeSplitter>(aDFace, aParams, "default");
      break;
  }
  return 0;
}

// --- TEMP T-123 (`--algospy <defl> <angle> <faceIdx> [<faceIdx> ...]`) ---------------
// Replays `BRepMesh_FaceDiscret::process` for the given model faces with the real
// `BRepMesh_DelaunayNodeInsertionMeshAlgo<...RangeSplitter>` that
// `BRepMesh_MeshAlgoFactory::GetAlgo` (`BRepMesh_MeshAlgoFactory.cxx:60-118`) hands out
// for a cone, subclassed so the internals become observable. Two and only two branches
// raise `IMeshData_Failure` for a face whose status flags are clean before face
// discretization (`--ctxdbg` shows that): `BRepMesh_NodeInsertionMeshAlgo::
// initDataStructure` (`BRepMesh_NodeInsertionMeshAlgo.hxx:61-111`) returning false after
// `AdjustRange`/`IsValid`, and `BRepMesh_BaseMeshAlgo::commitSurfaceTriangulation`
// (`BRepMesh_BaseMeshAlgo.cxx:228-241`) finding `ElementsOfDomain()` empty - which
// `collectTriangles` (`cxx:245-252`) reports as a null triangulation. The `init` line
// tells them apart: no `init` line at all = the first branch, `elems=0` after `gen` = the
// second, and `elems` collapsing between `gen_pre`/`gen_post`/`ins_pre`/`ins_post`
// locates the loop (`BRepMesh_Delaun` base pass vs `insertNodes`) that empties it.
// `elems` is the `collectTriangles` triangle count, `links` the contour link count
// `initDataStructure` built, so the same lines also prove whether the wire loops of
// `initDataStructure` (skipped for `SelfIntersectingWire`, or `OpenWire` on non-outer
// wires) contributed anything.
namespace
{
template <class Splitter>
class ZZAlgoSpy
    : public BRepMesh_DelaunayNodeInsertionMeshAlgo<Splitter, BRepMesh_DelaunayBaseMeshAlgo>
{
  using Base =
    BRepMesh_DelaunayNodeInsertionMeshAlgo<Splitter, BRepMesh_DelaunayBaseMeshAlgo>;

protected:
  bool initDataStructure() override
  {
    const bool aRes = Base::initDataStructure();
    dump("init", aRes ? 1 : 0);
    return aRes;
  }

  void generateMesh(const Message_ProgressRange& theRange) override
  {
    dump("gen_pre", -1);
    dumpMeshInput();
    Base::generateMesh(theRange);
    dump("gen_post", -1);
  }

  void postProcessMesh(BRepMesh_Delaun& theMesher, const Message_ProgressRange& theRange) override
  {
    dump("ins_pre", -1);
    Base::postProcessMesh(theMesher, theRange);
    dump("ins_post", -1);
  }

private:
  void dump(const char* theTag, const int theRes)
  {
    const IMeshData::IFaceHandle&                      aDFace  = this->getDFace();
    const Splitter&                                    aSplit  = this->getRangeSplitter();
    const occ::handle<BRepMesh_DataStructureOfDelaun>& aStruct = this->getStructure();
    const std::pair<double, double>&                   aRU     = aSplit.GetRangeU();
    const std::pair<double, double>&                   aRV     = aSplit.GetRangeV();
    const std::pair<double, double>&                   aDelta  = aSplit.GetDelta();
    std::printf("ALGOSPY %s res=%d valid=%d nodes=%d links=%d elems=%d linkdom=%d "
                "ru=[%.17g,%.17g] rv=[%.17g,%.17g] delta=[%.17g,%.17g]\n",
                theTag,
                theRes,
                const_cast<Splitter&>(aSplit).IsValid() ? 1 : 0,
                aStruct.IsNull() ? -1 : aStruct->NbNodes(),
                aStruct.IsNull() ? -1 : aStruct->NbLinks(),
                aStruct.IsNull() ? -1 : aStruct->ElementsOfDomain().Extent(),
                aStruct.IsNull() ? -1 : aStruct->LinksOfDomain().Extent(),
                aRU.first,
                aRU.second,
                aRV.first,
                aRV.second,
                aDelta.first,
                aDelta.second);
    std::fflush(stdout);
    if (std::string(theTag) == "gen_pre")
    {
      dumpPolygon(theTag);
    }
  }

  // Dumps the two inputs `BRepMesh_DelaunayBaseMeshAlgo::generateMesh` (`cxx:34-45`)
  // derives from the structure: the vertex list length (`aVerticesOrder` = one entry per
  // registered 3D point, `cxx:38-41`) and the acceleration-grid cell counts
  // (`getCellsCount` -> `BRepMesh_GeomTool::CellsCount`, `cxx:42`), plus the node tool
  // tolerance/cell size `BRepMesh_NodeInsertionMeshAlgo::initDataStructure`
  // (`hxx:86-92`) set. The cell size decides which triangles `BRepMesh_CircleTool::Select`
  // offers to `createTrianglesOnNewVertices` for a candidate vertex, so a difference here
  // is enough to change a degenerate triangulation's output.
  void dumpMeshInput() const
  {
    const occ::handle<BRepMesh_DataStructureOfDelaun>& aStruct = this->getStructure();
    if (aStruct.IsNull())
    {
      return;
    }
    double aTolU = 0.0, aTolV = 0.0;
    aStruct->Data()->GetTolerance(aTolU, aTolV);
    const std::pair<int, int> aCells =
      const_cast<ZZAlgoSpy*>(this)->getCellsCount(this->getNodesMap()->Length());
    std::printf("ALGOSPYM nodes=%d vorder=%d cells=[%d,%d] tol=[%.17g,%.17g]\n",
                aStruct->NbNodes(),
                this->getNodesMap()->Length(),
                aCells.first,
                aCells.second,
                aTolU,
                aTolV);
    std::fflush(stdout);
  }

  // Dumps the structure `BRepMesh_Delaun` is constructed from in
  // `BRepMesh_DelaunayBaseMeshAlgo::generateMesh` (`cxx:36-47`): every node in the scaled
  // basis (`BRepMesh_DefaultRangeSplitter::Scale(thePoint, true)`) plus the links of the
  // domain, which is the boundary polygon the mesher has to triangulate.
  void dumpPolygon(const char* theTag) const
  {
    const occ::handle<BRepMesh_DataStructureOfDelaun>& aStruct = this->getStructure();
    if (aStruct.IsNull())
    {
      return;
    }
    for (int i = 1; i <= aStruct->NbNodes(); ++i)
    {
      const BRepMesh_Vertex& aV = aStruct->GetNode(i);
      std::printf("ALGOSPYN %s %d %.17g %.17g %d\n",
                  theTag,
                  i,
                  aV.Coord().X(),
                  aV.Coord().Y(),
                  static_cast<int>(aV.Movability()));
    }
    for (IMeshData::IteratorOfMapOfInteger aLinkIt(aStruct->LinksOfDomain()); aLinkIt.More();
         aLinkIt.Next())
    {
      const BRepMesh_Edge& aL = aStruct->GetLink(aLinkIt.Key());
      std::printf("ALGOSPYP %s %d %d %d %d\n",
                  theTag,
                  aLinkIt.Key(),
                  aL.FirstNode(),
                  aL.LastNode(),
                  static_cast<int>(aL.Movability()));
    }
    std::fflush(stdout);
  }
};
} // namespace

// Maps a surface type to the algo `BRepMesh_MeshAlgoFactory::GetAlgo`
// (`BRepMesh_MeshAlgoFactory.cxx:60-118`) picks for the default parameters used here
// (`InternalVerticesMode = false`, `EnableControlSurfaceDeflectionAllSurfaces = false`),
// i.e. the node-insertion algo for every type this replay supports.
static occ::handle<IMeshTools_MeshAlgo> makeSpyAlgo(const GeomAbs_SurfaceType theType)
{
  switch (theType)
  {
    case GeomAbs_Cone:
      return new ZZAlgoSpy<BRepMesh_ConeRangeSplitter>;
    default:
      return occ::handle<IMeshTools_MeshAlgo>();
  }
}

int runAlgoSpy(const TopoDS_Shape& theShape,
                      const double        theDeflection,
                      const double        theAngle,
                      const int*          theFaces,
                      const int           theFacesNb)
{
  IMeshTools_Parameters aParams;
  aParams.Deflection         = theDeflection;
  aParams.DeflectionInterior = theDeflection;
  aParams.Angle              = theAngle;
  aParams.AngleInterior      = 2.0 * theAngle;
  aParams.MinSize =
    (std::max)(IMeshTools_Parameters::RelMinSize() * theDeflection, Precision::Confusion());
  aParams.AdjustMinSize = false;
  aParams.InParallel    = false;
  aParams.Relative      = false;
  aParams.CleanModel    = false;

  occ::handle<BRepMesh_Context> aCtx = new BRepMesh_Context;
  aCtx->SetShape(theShape);
  aCtx->ChangeParameters() = aParams;

  if (!aCtx->BuildModel() || !aCtx->DiscretizeEdges() || !aCtx->HealModel()
      || !aCtx->PreProcessModel())
  {
    std::cout << "ALGOSPY stage-failed\n";
    return 0;
  }

  const occ::handle<IMeshData_Model>& aModel = aCtx->GetModel();
  if (aModel.IsNull())
  {
    std::cout << "ALGOSPY model-null\n";
    return 0;
  }

  for (int i = 0; i < theFacesNb; ++i)
  {
    if (theFaces[i] < 0 || theFaces[i] >= aModel->FacesNb())
    {
      continue;
    }
    const IMeshData::IFaceHandle&   aDFace = aModel->GetFace(theFaces[i]);
    const GeomAbs_SurfaceType       aType  = aDFace->GetSurface()->GetType();
    std::printf("ALGOSPY face f=%d type=%d wires=%d fail=%d reusable=%d\n",
                theFaces[i],
                static_cast<int>(aType),
                aDFace->WiresNb(),
                aDFace->IsSet(IMeshData_Failure) ? 1 : 0,
                aDFace->IsSet(IMeshData_Reused) ? 1 : 0);
    for (int w = 0; w < aDFace->WiresNb(); ++w)
    {
      const IMeshData::IWireHandle& aDWire = aDFace->GetWire(w);
      std::printf("ALGOSPY wire f=%d w=%d edges=%d selfint=%d open=%d\n",
                  theFaces[i],
                  w,
                  aDWire->EdgesNb(),
                  aDWire->IsSet(IMeshData_SelfIntersectingWire) ? 1 : 0,
                  aDWire->IsSet(IMeshData_OpenWire) ? 1 : 0);
    }
    std::fflush(stdout);

    // `BRepMesh_FaceDiscret::process` (`BRepMesh_FaceDiscret.cxx:99-131`).
    if (aDFace->IsSet(IMeshData_Failure) || aDFace->IsSet(IMeshData_Reused))
    {
      std::printf("ALGOSPY skip f=%d\n", theFaces[i]);
      continue;
    }
    occ::handle<IMeshTools_MeshAlgo> aAlgo = makeSpyAlgo(aType);
    if (aAlgo.IsNull())
    {
      std::printf("ALGOSPY unsupported f=%d type=%d\n", theFaces[i], static_cast<int>(aType));
      continue;
    }
    aAlgo->Perform(aDFace, aParams, Message_ProgressRange());

    TopLoc_Location                        aLoc;
    const occ::handle<Poly_Triangulation>& aTri =
      BRep_Tool::Triangulation(aDFace->GetFace(), aLoc);
    std::printf("ALGOSPY final f=%d fail=%d tri=%d/%d\n",
                theFaces[i],
                aDFace->IsSet(IMeshData_Failure) ? 1 : 0,
                aTri.IsNull() ? -1 : aTri->NbNodes(),
                aTri.IsNull() ? -1 : aTri->NbTriangles());
    std::fflush(stdout);
  }
  return 0;
}
