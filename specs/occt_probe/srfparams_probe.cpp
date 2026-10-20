// Splitter-parameter dumps: `--srparams`, `--facestats`, `--facedump`, `--edgedu`.
#include "occt_probe_common.h"

// --- TEMP T-101 (`--facestats <defl> <angle>`) -------------------------------------
// Per-model-face `nodes/triangles/bbox` with the *same* mesher parameters the port
// uses, in the model's own face order (`IMeshData_Model::GetFace`), so a face can be
// paired with the port's `zz_probe_a3n00 --fstats` by the six bbox coordinates.
// `--uvsum` cannot be used for this: that branch meshes with a hard-coded 0.5 rad
// angle (its per-face `triangles=` sums to 8054 on a3n00 while the 20-degree mesh has
// 12324), so its counts are not comparable with the port's.
// TEMP T-132 (`--srparams <faceIdx> <defl> <angle>`): prints the raw U/V parameter
// sets a factory-selected UV range splitter hands to its grain filter, together with
// the grid those parameters produce (`BRepMesh_NURBSRangeSplitter::
// GenerateSurfaceNodes`, BRepMesh_NURBSRangeSplitter.cxx:314-400), so the port's
// `NURBSGRID` dump can be compared entry by entry instead of by node totals.
static void dumpRawParams(const char* theLabel, const IMeshData::IMapOfReal& theParams)
{
  std::vector<double> aVals;
  aVals.reserve(theParams.Extent());
  for (int i = 1; i <= theParams.Extent(); ++i)
  {
    aVals.push_back(theParams(i));
  }
  std::sort(aVals.begin(), aVals.end());
  std::printf("OCCTSRP param %s n=%d vals=", theLabel, static_cast<int>(aVals.size()));
  for (size_t i = 0; i < aVals.size(); ++i)
  {
    std::printf("%s%.17g", (i == 0) ? "" : ",", aVals[i]);
  }
  std::printf("\n");
}

// TEMP T-132: same splitter/algo pairing as `ZZSrfNodeCountAlgo`, but the dump is the
// splitter's own parameter state (`GetParametersU/V`, `BRepMesh_UVParamRangeSplitter.
// hxx:46-58`) plus the generated grid. The parameters are filled during
// `initDataStructure` -> `collectWirePoints` -> `AddPoint`
// (`BRepMesh_NodeInsertionMeshAlgo.hxx:143-183`, `BRepMesh_BoundaryParamsRangeSplitter.
// hxx:33-38`) and never change afterwards, so reading them in `postProcessMesh` is the
// same state `insertNodes` filters.
template <class Splitter>
class ZZSrfParamsAlgo
    : public BRepMesh_DelaunayNodeInsertionMeshAlgo<Splitter, BRepMesh_DelaunayBaseMeshAlgo>
{
  typedef BRepMesh_DelaunayNodeInsertionMeshAlgo<Splitter, BRepMesh_DelaunayBaseMeshAlgo> Base;

protected:
  void postProcessMesh(BRepMesh_Delaun& theMesher, const Message_ProgressRange& theRange) override
  {
    const Splitter&                 aSplitter = this->getRangeSplitter();
    const IMeshData::IFaceHandle&   aDFace    = this->getDFace();
    const std::pair<double, double>& aRangeU   = aSplitter.GetRangeU();
    const std::pair<double, double>& aRangeV   = aSplitter.GetRangeV();
    const std::pair<double, double>& aDelta    = aSplitter.GetDelta();
    const std::pair<double, double>& aTol      = aSplitter.GetToleranceUV();
    const IMeshTools_Parameters&    aParams   = this->getParameters();
    const double                    aDefl     = aDFace->GetDeflection();
    std::printf("OCCTSRP face type=%d valid=%d ru=(%.17g,%.17g) rv=(%.17g,%.17g) delta=(%.17g,%.17g) "
                "tol=(%.17g,%.17g) ures=%.17g vres=%.17g facedefl=%.17g minsize=%.17g "
                "angint=%.17g angle=%.17g defl=%.17g\n",
                static_cast<int>(aDFace->GetSurface()->GetType()),
                const_cast<Splitter&>(aSplitter).IsValid() ? 1 : 0,
                aRangeU.first,
                aRangeU.second,
                aRangeV.first,
                aRangeV.second,
                aDelta.first,
                aDelta.second,
                aTol.first,
                aTol.second,
                aSplitter.GetSurface()->UResolution(aDefl),
                aSplitter.GetSurface()->VResolution(aDefl),
                aDefl,
                aParams.MinSize,
                aParams.AngleInterior,
                aParams.Angle,
                aParams.Deflection);
    dumpRawParams("u", aSplitter.GetParametersU());
    dumpRawParams("v", aSplitter.GetParametersV());
    {
      const occ::handle<BRepAdaptor_Surface>& aSurf = aSplitter.GetSurface();
      const int aNbU = aSurf->NbUIntervals(GeomAbs_CN);
      const int aNbV = aSurf->NbVIntervals(GeomAbs_CN);
      TColStd_Array1OfReal aIvU(1, (std::max)(aNbU + 1, 2));
      TColStd_Array1OfReal aIvV(1, (std::max)(aNbV + 1, 2));
      aSurf->UIntervals(aIvU, GeomAbs_CN);
      aSurf->VIntervals(aIvV, GeomAbs_CN);
      std::printf("OCCTSRP uint nb=%d vals=", aNbU);
      for (int i = aIvU.Lower(); i <= aIvU.Upper(); ++i)
      {
        std::printf("%s%.17g", (i == aIvU.Lower()) ? "" : ",", aIvU.Value(i));
      }
      std::printf("\nOCCTSRP vint nb=%d vals=", aNbV);
      for (int i = aIvV.Lower(); i <= aIvV.Upper(); ++i)
      {
        std::printf("%s%.17g", (i == aIvV.Lower()) ? "" : ",", aIvV.Value(i));
      }
      std::printf("\n");
    }

    // `GenerateSurfaceNodes` fills the parameter maps on its first call for the
    // splitters that derive them lazily (`BRepMesh_NURBSRangeSplitter::
    // initParameters`, cxx:452-488, writes through the const_casts at :469/:476),
    // so the grid is built first and both maps are dumped a second time after it.
    const Handle(IMeshData::ListOfPnt2d) aNodes = aSplitter.GenerateSurfaceNodes(aParams);
    dumpRawParams("u_after", aSplitter.GetParametersU());
    dumpRawParams("v_after", aSplitter.GetParametersV());
    std::printf("OCCTSRP nodes=%d\n", aNodes.IsNull() ? 0 : static_cast<int>(aNodes->Size()));
    std::printf("OCCTSRP meshpre nodes=%d tris=%d\n",
                theMesher.Result()->NbNodes(),
                theMesher.Result()->NbElements());
    if (!aNodes.IsNull())
    {
      for (IMeshData::ListOfPnt2d::Iterator anIt(*aNodes); anIt.More(); anIt.Next())
      {
        std::printf("OCCTSRPN u=%.17g v=%.17g\n", anIt.Value().X(), anIt.Value().Y());
      }
    }
    Base::postProcessMesh(theMesher, theRange);
    std::printf("OCCTSRP meshpost nodes=%d tris=%d\n",
                theMesher.Result()->NbNodes(),
                theMesher.Result()->NbElements());
  }
};

int runSrfParams(const TopoDS_Shape& theShape, const int theFaceIdx,
                        const double theDeflection, const double theAngle)
{
  IMeshTools_Parameters aParams;
  aParams.Deflection    = theDeflection;
  aParams.Angle         = theAngle;
  aParams.InParallel    = false;
  aParams.Relative      = false;
  aParams.MinSize       = 0.0;
  aParams.AdjustMinSize = false;
  // `BRepMesh_IncrementalMesh::initParameters` (BRepMesh_IncrementalMesh.hxx:79-108),
  // i.e. the parameter state the real mesher runs with.
  if (aParams.DeflectionInterior < Precision::Confusion())
  {
    aParams.DeflectionInterior = aParams.Deflection;
  }
  if (aParams.MinSize < Precision::Confusion())
  {
    aParams.MinSize =
      (std::max)(IMeshTools_Parameters::RelMinSize()
                   * (std::min)(aParams.Deflection, aParams.DeflectionInterior),
                 Precision::Confusion());
  }
  if (aParams.AngleInterior < Precision::Angular())
  {
    aParams.AngleInterior = 2.0 * aParams.Angle;
  }
  occ::handle<BRepMesh_Context> aCtx = new BRepMesh_Context;
  aCtx->SetShape(theShape);
  aCtx->ChangeParameters()            = aParams;
  aCtx->ChangeParameters().CleanModel = false;
  const bool aBuilt   = aCtx->BuildModel();
  const bool aDiscret = aBuilt && aCtx->DiscretizeEdges();
  const bool aHealed  = aDiscret && aCtx->HealModel();
  const bool aPreproc = aHealed && aCtx->PreProcessModel();
  const occ::handle<IMeshData_Model>& aModel = aCtx->GetModel();
  if (aModel.IsNull() || theFaceIdx < 0 || theFaceIdx >= aModel->FacesNb())
  {
    std::cout << "SRPARAM model-null-or-range\n";
    return 0;
  }
  const IMeshData::IFaceHandle& aDFace = aModel->GetFace(theFaceIdx);
  std::printf("OCCTSRP shape defl=%.17g angle=%.17g angint=%.17g minsize=%.17g built=%d discret=%d "
              "healed=%d preproc=%d faces=%d f=%d\n",
              aParams.Deflection,
              aParams.Angle,
              aParams.AngleInterior,
              aParams.MinSize,
              aBuilt ? 1 : 0,
              aDiscret ? 1 : 0,
              aHealed ? 1 : 0,
              aPreproc ? 1 : 0,
              aModel->FacesNb(),
              theFaceIdx);
  // TEMP: post-heal model pcurve nodes, so a UV end value that differs from the
  // port can be read from the same data the range splitter consumes.
  for (int w = 0; w < aDFace->WiresNb(); ++w)
  {
    const IMeshData::IWireHandle& aW = aDFace->GetWire(w);
    for (int e = 0; e < aW->EdgesNb(); ++e)
    {
      const IMeshData::IEdgePtr& aE = aW->GetEdge(e);
      for (int p = 0; p < aE->PCurvesNb(); ++p)
      {
        const IMeshData::IPCurveHandle& aP = aE->GetPCurve(p);
        if (aP->GetFace() != aDFace)
        {
          continue;
        }
        if (std::getenv("SRPCDEF") != nullptr)
        {
          // Geometry actually evaluated by `BRepMesh_EdgeParameterProvider`:
          // `BRepAdaptor_Curve::Initialize(oriented edge, face)` ->
          // `BRep_Tool::CurveOnSurface(E, F, pf, pl)` ->
          // `Geom2dAdaptor_Curve::Load(PC, pf, pl)`.
          const TopoDS_Edge anOriented = TopoDS::Edge(aE->GetEdge().Oriented(aW->GetEdgeOrientation(e)));
          double f = 0., l = 0.;
          occ::handle<Geom2d_Curve> aPC = BRep_Tool::CurveOnSurface(anOriented, aDFace->GetFace(), f, l);
          if (!aPC.IsNull())
          {
            occ::handle<Geom2d_Curve> aBasis = aPC;
            if (aPC->IsKind(STANDARD_TYPE(Geom2d_TrimmedCurve)))
            {
              aBasis = occ::down_cast<Geom2d_TrimmedCurve>(aPC)->BasisCurve();
            }
            std::printf("OCCTSRPDEF w=%d e=%d p=%d type=%s pf=%.17g pl=%.17g",
                        w,
                        e,
                        p,
                        aBasis->DynamicType()->Name(),
                        f,
                        l);
            if (aBasis->IsKind(STANDARD_TYPE(Geom2d_BSplineCurve)))
            {
              const occ::handle<Geom2d_BSplineCurve> aBS = occ::down_cast<Geom2d_BSplineCurve>(aBasis);
              std::printf(" deg=%d npoles=%d nknots=%d rational=%d\n",
                          aBS->Degree(),
                          aBS->NbPoles(),
                          aBS->NbKnots(),
                          aBS->IsRational() ? 1 : 0);
              for (int k = 1; k <= aBS->NbPoles(); ++k)
              {
                const gp_Pnt2d aPole = aBS->Pole(k);
                std::printf("OCCTSRPDEF  pole[%d]=%.17g,%.17g\n", k, aPole.X(), aPole.Y());
              }
              // TEMP: knot values and the interpolant's values there, i.e. the
              // projected UV samples the curve was built through (a clamped
              // interpolant passes through them), matching the port's `knotval`.
              for (int k = 1; k <= aBS->NbKnots(); ++k)
              {
                const gp_Pnt2d aVal = aBS->Value(aBS->Knot(k));
                std::printf("OCCTSRPDEF  knot[%d]=%.17g m=%d val=%.17g,%.17g\n",
                            k,
                            aBS->Knot(k),
                            aBS->Multiplicity(k),
                            aVal.X(),
                            aVal.Y());
              }
            }
            else
            {
              std::printf("\n");
            }
          }
          // TEMP: the pcurve stored for the OPPOSITE edge orientation. On a seam
          // edge `BRep_Tool::CurveOnSurface` returns the other of the two
          // `CurveOnClosedSurface` branches, so this prints the pcurve the wire
          // edge at hand does not use (`OCCTSRP pc` above lists both point sets,
          // which is why the per-orientation poles are needed to tell which
          // branch carries the `+-2*PI` rounding).
          {
            const TopoDS_Edge anOpp = TopoDS::Edge(aE->GetEdge().Oriented(
              aW->GetEdgeOrientation(e) == TopAbs_FORWARD ? TopAbs_REVERSED : TopAbs_FORWARD));
            double of = 0., ol = 0.;
            occ::handle<Geom2d_Curve> aOppC = BRep_Tool::CurveOnSurface(anOpp, aDFace->GetFace(), of, ol);
            if (!aOppC.IsNull())
            {
              occ::handle<Geom2d_Curve> aOppBasis = aOppC;
              if (aOppC->IsKind(STANDARD_TYPE(Geom2d_TrimmedCurve)))
              {
                aOppBasis = occ::down_cast<Geom2d_TrimmedCurve>(aOppC)->BasisCurve();
              }
              std::printf("OCCTSRPDEF2 w=%d e=%d type=%s pf=%.17g pl=%.17g\n",
                          w,
                          e,
                          aOppBasis->DynamicType()->Name(),
                          of,
                          ol);
              if (aOppBasis->IsKind(STANDARD_TYPE(Geom2d_BSplineCurve)))
              {
                const occ::handle<Geom2d_BSplineCurve> aOBS =
                  occ::down_cast<Geom2d_BSplineCurve>(aOppBasis);
                std::printf("OCCTSRPDEF2  deg=%d npoles=%d rational=%d\n",
                            aOBS->Degree(),
                            aOBS->NbPoles(),
                            aOBS->IsRational() ? 1 : 0);
                for (int k = 1; k <= aOBS->NbPoles(); ++k)
                {
                  const gp_Pnt2d aPole = aOBS->Pole(k);
                  std::printf("OCCTSRPDEF2  pole[%d]=%.17g,%.17g\n", k, aPole.X(), aPole.Y());
                }
              }
              else if (aOppBasis->IsKind(STANDARD_TYPE(Geom2d_Line)))
              {
                const occ::handle<Geom2d_Line> aOLine = occ::down_cast<Geom2d_Line>(aOppBasis);
                std::printf("OCCTSRPDEF2  line p=(%.17g,%.17g) d=(%.17g,%.17g)\n",
                            aOLine->Location().X(),
                            aOLine->Location().Y(),
                            aOLine->Direction().X(),
                            aOLine->Direction().Y());
              }
              else
              {
                std::printf("\n");
              }
            }
          }
        }
        if (p == 0 && std::getenv("SRPCC2") != nullptr)
        {
          // TEMP: both STORED 2D representations of the edge, i.e. what
          // `BRep_Tool::CurveOnClosedSurface` used to hand out for a seam edge
          // (`CurveOnSurface(E, C, S, L, f, l, Index)` in 8.0). The port keeps
          // the same pair in `GeometryRegistry::edge_pcurves`.
          occ::handle<Geom2d_Curve> aC1, aC2;
          occ::handle<Geom_Surface> aS1, aS2;
          TopLoc_Location           aL1, aL2;
          double                    aF1 = 0., aL1v = 0., aF2 = 0., aL2v = 0.;
          BRep_Tool::CurveOnSurface(aE->GetEdge(), aC1, aS1, aL1, aF1, aL1v, 1);
          BRep_Tool::CurveOnSurface(aE->GetEdge(), aC2, aS2, aL2, aF2, aL2v, 2);
          std::printf("OCCTPCC2 w=%d e=%d c1null=%d c2null=%d same=%d f1=%.17g l1=%.17g "
                      "f2=%.17g l2=%.17g\n",
                      w,
                      e,
                      aC1.IsNull() ? 1 : 0,
                      aC2.IsNull() ? 1 : 0,
                      (!aC1.IsNull() && !aC2.IsNull() && aC1 == aC2) ? 1 : 0,
                      aF1,
                      aL1v,
                      aF2,
                      aL2v);
          for (int q = 0; q < 2; ++q)
          {
            const occ::handle<Geom2d_Curve>& aC = (q == 0) ? aC1 : aC2;
            if (aC.IsNull())
            {
              continue;
            }
            std::printf("OCCTPCC2  q=%d type=%s\n", q, aC->DynamicType()->Name());
            if (aC->IsKind(STANDARD_TYPE(Geom2d_BSplineCurve)))
            {
              const occ::handle<Geom2d_BSplineCurve> aB = occ::down_cast<Geom2d_BSplineCurve>(aC);
              std::printf("OCCTPCC2  q=%d deg=%d np=%d nk=%d\n",
                          q,
                          aB->Degree(),
                          aB->NbPoles(),
                          aB->NbKnots());
              for (int k = 1; k <= aB->NbPoles(); ++k)
              {
                const gp_Pnt2d aPole = aB->Pole(k);
                std::printf("OCCTPCC2  q=%d pole[%d]=%.17g,%.17g\n",
                            q,
                            k,
                            aPole.X(),
                            aPole.Y());
              }
            }
            else if (aC->IsKind(STANDARD_TYPE(Geom2d_Line)))
            {
              const occ::handle<Geom2d_Line> aL = occ::down_cast<Geom2d_Line>(aC);
              std::printf("OCCTPCC2  q=%d line p=%.17g,%.17g d=%.17g,%.17g\n",
                          q,
                          aL->Location().X(),
                          aL->Location().Y(),
                          aL->Direction().X(),
                          aL->Direction().Y());
            }
          }
        }
        std::printf("OCCTSRP pc w=%d e=%d p=%d ori=%d int=%d n=%d",
                    w,
                    e,
                    p,
                    static_cast<int>(aW->GetEdgeOrientation(e)),
                    aP->IsInternal() ? 1 : 0,
                    aP->ParametersNb());
        for (int k = 0; k < aP->ParametersNb(); ++k)
        {
          std::printf(" [%d]%.17g,%.17g,%.17g",
                      k,
                      aP->GetPoint(k).X(),
                      aP->GetPoint(k).Y(),
                      aP->GetParameter(k));
        }
        std::printf("\n");
      }
    }
  }

  occ::handle<IMeshTools_MeshAlgo> aAlgo;
  switch (aDFace->GetSurface()->GetType())
  {
    case GeomAbs_SurfaceOfRevolution:
      aAlgo = new ZZSrfParamsAlgo<BRepMesh_BoundaryParamsRangeSplitter>();
      break;
    case GeomAbs_Torus:
      aAlgo = new ZZSrfParamsAlgo<BRepMesh_TorusRangeSplitter>();
      break;
    case GeomAbs_SurfaceOfExtrusion:
      aAlgo = new ZZSrfParamsAlgo<BRepMesh_ExtrusionRangeSplitter>();
      break;
    case GeomAbs_BezierSurface:
    case GeomAbs_BSplineSurface:
      aAlgo = new ZZSrfParamsAlgo<BRepMesh_NURBSRangeSplitter>();
      break;
    default:
      std::printf("OCCTSRP unsupported stype=%d\n",
                  static_cast<int>(aDFace->GetSurface()->GetType()));
      return 0;
  }
  aAlgo->Perform(aDFace, aParams, Message_ProgressRange());
  return 0;
}

int runFaceStats(const TopoDS_Shape& theShape, const double theDeflection,
                        const double theAngle)
{
  IMeshTools_Parameters aParams;
  aParams.Deflection    = theDeflection;
  aParams.Angle         = theAngle;
  aParams.InParallel    = false;
  aParams.Relative      = false;
  aParams.MinSize       = Precision::Confusion();
  aParams.AdjustMinSize = false;

  occ::handle<BRepMesh_ModelBuilder> aBuilder = new BRepMesh_ModelBuilder;
  const occ::handle<IMeshData_Model> aModel   = aBuilder->Perform(theShape, aParams);
  if (aModel.IsNull())
  {
    std::cout << "FSTAT_OCC model-null\n";
    return 0;
  }
  {
    BRepMesh_IncrementalMesh aMesher(theShape, theDeflection, false, theAngle);
  }
  std::cout << "FSTAT_OCC defl=" << theDeflection << " angle=" << theAngle
            << " faces=" << aModel->FacesNb() << std::endl;
  for (int f = 0; f < aModel->FacesNb(); ++f)
  {
    const IMeshData::IFaceHandle& aDFace = aModel->GetFace(f);
    const TopoDS_Face&            aFace  = aDFace->GetFace();
    Bnd_Box                       aBox;
    BRepBndLib::Add(aFace, aBox, false);
    double v[6] = {0, 0, 0, 0, 0, 0};
    if (!aBox.IsVoid())
    {
      aBox.Get(v[0], v[1], v[2], v[3], v[4], v[5]);
    }
    TopLoc_Location                     aLoc;
    const occ::handle<Poly_Triangulation>& aTri = BRep_Tool::Triangulation(aFace, aLoc);
    const int aNodes = aTri.IsNull() ? 0 : aTri->NbNodes();
    const int aTris  = aTri.IsNull() ? 0 : aTri->NbTriangles();
    std::cout << std::fixed << std::setprecision(6) << "FSTAT_OCC face=" << f
              << " nodes=" << aNodes << " triangles=" << aTris << " bbox=(" << v[0] << "," << v[1]
              << "," << v[2] << ")-(" << v[3] << "," << v[4] << "," << v[5] << ")" << std::endl;
  }
  return 0;
}

// --- TEMP T-121 (`--edgedu <faceIdx> <defl> <angle>`) -------------------------------
// TEMP T-101b (`--facedump <faceIdx> <defl> <angle>`): the node and triangle arrays of
// one *model* face (`IMeshData_Model::GetFace`, same indexing as `--facestats`), so a
// per-face count divergence can be read as geometry instead of a total.
int runFaceDump(const TopoDS_Shape& theShape, const int theFaceIdx,
                       const double theDeflection, const double theAngle)
{
  {
    BRepMesh_IncrementalMesh aMesher(theShape, theDeflection, false, theAngle);
  }
  IMeshTools_Parameters aParams;
  aParams.Deflection    = theDeflection;
  aParams.Angle         = theAngle;
  aParams.InParallel    = false;
  aParams.Relative      = false;
  aParams.MinSize       = Precision::Confusion();
  aParams.AdjustMinSize = false;
  occ::handle<BRepMesh_ModelBuilder> aBuilder = new BRepMesh_ModelBuilder;
  const occ::handle<IMeshData_Model> aModel   = aBuilder->Perform(theShape, aParams);
  if (aModel.IsNull() || theFaceIdx < 0 || theFaceIdx >= aModel->FacesNb())
  {
    std::cout << "FDUMP model-null-or-range\n";
    return 0;
  }
  const TopoDS_Face& aFace = aModel->GetFace(theFaceIdx)->GetFace();
  TopLoc_Location    aLoc;
  const occ::handle<Poly_Triangulation>& aTri = BRep_Tool::Triangulation(aFace, aLoc);
  std::cout << std::fixed << std::setprecision(9) << "FDUMP face=" << theFaceIdx
            << " nodes=" << (aTri.IsNull() ? 0 : aTri->NbNodes())
            << " tris=" << (aTri.IsNull() ? 0 : aTri->NbTriangles()) << std::endl;
  if (aTri.IsNull())
    return 0;
  for (int i = 1; i <= aTri->NbNodes(); ++i)
  {
    const gp_Pnt aP = aTri->Node(i).Transformed(aLoc.Transformation());
    std::cout << "FDUMP n i=" << i;
    if (aTri->HasUVNodes())
    {
      const gp_Pnt2d aUV = aTri->UVNode(i);
      std::cout << " uv=(" << aUV.X() << "," << aUV.Y() << ")";
    }
    std::cout << " p=(" << aP.X() << "," << aP.Y() << "," << aP.Z() << ")" << std::endl;
  }
  for (int i = 1; i <= aTri->NbTriangles(); ++i)
  {
    int a1, a2, a3;
    aTri->Triangle(i).Get(a1, a2, a3);
    std::cout << "FDUMP t i=" << i << " (" << a1 << "," << a2 << "," << a3 << ")" << std::endl;
  }
  return 0;
}

// --- TEMP T-121 (`--edgedu <faceIdx> <defl> <angle>`) -------------------------------// Prints, for every edge occurrence of one face, the `BRepAdaptor_Curve` type and the
// exact inputs of `GCPnts_TangentialDeflection::EvaluateDu` (`cxx:302-322`) at the
// edge's first parameter: `|T|`, `|N|`, `|N x T|` and `Ln = |N x T| / |T|`. The
// resulting `du_init = sqrt(8 * max(curvatureDeflection, minLen) / Ln)` is printed for
// both the face deflection and half of it, so a divergence in the first tessellation
// step can be attributed to the curve (`Ln`) or to the deflection (`D`).
int runEdgeDu(const TopoDS_Shape& theShape, const int theFaceIdx,
                     const double theDeflection, const double theAngle)
{
  {
    BRepMesh_IncrementalMesh aMesher(theShape, theDeflection, false, theAngle);
  }
  int aFaceIdx = 0;
  for (TopExp_Explorer aFaceEx(theShape, TopAbs_FACE); aFaceEx.More(); aFaceEx.Next(), ++aFaceIdx)
  {
    if (aFaceIdx != theFaceIdx)
      continue;
    const TopoDS_Face aFace = TopoDS::Face(aFaceEx.Current());
    int               anEdgeIdx = 0;
    for (TopExp_Explorer anEdgeEx(aFace, TopAbs_EDGE); anEdgeEx.More();
         anEdgeEx.Next(), ++anEdgeIdx)
    {
      const TopoDS_Edge anEdge = TopoDS::Edge(anEdgeEx.Current());
      std::printf("EDGEDU face=%d i=%d degenerated=%d\n",
                  aFaceIdx,
                  anEdgeIdx,
                  BRep_Tool::Degenerated(anEdge) ? 1 : 0);
      try
      {
        BRepAdaptor_Curve  ac(anEdge);
        const GeomAbs_CurveType aType = ac.GetType();
        const double       aF = ac.FirstParameter();
        const double       aL = ac.LastParameter();
        gp_Pnt             aP;
        gp_Vec             aT, aN;
        ac.D2(aF, aP, aT, aN);
        const double aLt    = aT.Magnitude();
        const double aNMag  = aN.Magnitude();
        const double aCross = aN.CrossMagnitude(aT);
        const double aLn    = (aLt > Precision::Confusion()) ? aCross / aLt : 0.0;
        std::printf("EDGEDU face=%d i=%d type=%d f=%.17g l=%.17g |T|=%.17g |N|=%.17g "
                    "|NxT|=%.17g Ln=%.17g p0=(%.17g,%.17g,%.17g) du(D)=%.17g du(D/2)=%.17g\n",
                    aFaceIdx,
                    anEdgeIdx,
                    static_cast<int>(aType),
                    aF,
                    aL,
                    aLt,
                    aNMag,
                    aCross,
                    aLn,
                    aP.X(),
                    aP.Y(),
                    aP.Z(),
                    (aLn > Precision::Confusion()) ? std::sqrt(8.0 * theDeflection / aLn) : 0.0,
                    (aLn > Precision::Confusion()) ? std::sqrt(8.0 * (theDeflection / 2.0) / aLn)
                                                   : 0.0);
        if (aType == GeomAbs_Circle)
        {
          std::printf("EDGEDU face=%d i=%d circle R=%.17g\n",
                      aFaceIdx,
                      anEdgeIdx,
                      ac.Circle().Radius());
        }
      }
      catch (const Standard_Failure& anErr)
      {
        std::printf("EDGEDU face=%d i=%d thrown=%s\n",
                    aFaceIdx,
                    anEdgeIdx,
                    anErr.GetMessageString() ? anErr.GetMessageString() : "?");
      }
    }
  }
  return 0;
}
