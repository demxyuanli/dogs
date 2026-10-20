// Delaunay input dumps: `--lens`, `--structdump`, `--dcstruct`, `--srfnodes`.
// The `ZZ*Algo` subclasses replay OCCT own node-insertion algo and print the
// structure it hands to `BRepMesh_Delaun`.
#include "occt_probe_common.h"

// --- TEMP T-104 (`--structdump <faceIdx>`) -----------------------------------------
// Prints the exact `BRepMesh_DataStructureOfDelaun` OCCT's own node-insertion algo
// hands to `BRepMesh_Delaun` for one discrete face, in the port dump's `N`/`L` line
// format (`crates/occt-topo/src/meshing/node_insertion.rs::dump_delaun_structure`),
// so the two files diff line by line. The model is built exactly like `--boundary`
// (same deflection/angle as the port's `prs3d_get_deflection(shape, 0.1)`), and the
// face index is the `IMeshData_Model` order - pair it with the port's face by the
// printed bbox, not by index.
template <class Splitter>
class ZZStructDumpAlgo
    : public BRepMesh_DelaunayNodeInsertionMeshAlgo<Splitter, BRepMesh_DelaunayBaseMeshAlgo>
{protected:
  // TEMP T-105: dumps, in the port's own `N`/`L`/`TOL`/`CELLS`/`V` format, the exact
  // structure the *real* algo passes to `BRepMesh_Delaun` (the moment `generateMesh`
  // is entered). Diffing this against the port's `delaun_in_f*.txt` says whether the
  // two sides really start from the same input.
  void generateMesh(const Message_ProgressRange& theRange) override
  {
    typedef BRepMesh_DelaunayNodeInsertionMeshAlgo<Splitter, BRepMesh_DelaunayBaseMeshAlgo> Base;
    const char* aPath = std::getenv("ZZDELAUNIN");
    if (aPath != nullptr)
    {
      const occ::handle<BRepMesh_DataStructureOfDelaun>& aSt = this->getStructure();
      const int aNb = this->getNodesMap()->Length();
      double    aTolX = 0.0, aTolY = 0.0;
      aSt->Data()->GetTolerance(aTolX, aTolY);
      const std::pair<int, int> aCells = this->getCellsCount(aNb);
      std::ofstream aFile(aPath);
      aFile << "TOL " << std::setprecision(18) << aTolX << " " << aTolY << "\n";
      aFile << "CELLS " << aCells.first << " " << aCells.second << "\n";
      for (int i = 1; i <= aSt->NbNodes(); ++i)
      {
        const BRepMesh_Vertex& aV = aSt->GetNode(i);
        aFile << "N " << i << " " << std::setprecision(18) << aV.Coord().X() << " "
              << aV.Coord().Y() << " " << static_cast<int>(aV.Movability()) << "\n";
      }
      for (int i = 1; i <= aSt->NbLinks(); ++i)
      {
        const BRepMesh_Edge& aL = aSt->GetLink(i);
        aFile << "L " << i << " " << aL.FirstNode() << " " << aL.LastNode() << " "
              << static_cast<int>(aL.Movability()) << "\n";
      }
      aFile << "V";
      for (int i = 1; i <= aNb; ++i)
      {
        aFile << " " << i;
      }
      aFile << "\n";
      aFile.close();
      std::printf("OCCTPOST DELAUNIN wrote %s nodes=%d links=%d vorder=%d tol=(%.17e,%.17e) "
                  "cells=(%d,%d)\n",
                  aPath,
                  aSt->NbNodes(),
                  aSt->NbLinks(),
                  aNb,
                  aTolX,
                  aTolY,
                  aCells.first,
                  aCells.second);
    }
    Base::generateMesh(theRange);
  }

  // TEMP T-105: the real factory algo (`BRepMesh_MeshAlgoFactory`, torus ->
  // `NodeInsertionMeshAlgo<TorusRangeSplitter>`) is a template, so the only way to
  // observe the state it hands `BRepMesh_Delaun` and the state it leaves behind is to
  // subclass it. `postProcessMesh` is the seam between the base Delaunay
  // (`generateMesh`) and the surface-node insertion (`insertNodes` -> `AddVertices`),
  // which is exactly where the port's torus faces lose every triangle.
  void postProcessMesh(BRepMesh_Delaun& theMesher, const Message_ProgressRange& theRange) override
  {
    typedef BRepMesh_DelaunayNodeInsertionMeshAlgo<Splitter, BRepMesh_DelaunayBaseMeshAlgo> Base;
    dumpAlgoState("POST_BASEMESH", theMesher);
    {
      double aTolX = 0.0, aTolY = 0.0;
      const occ::handle<BRepMesh_DataStructureOfDelaun>& aSt = this->getStructure();
      aSt->Data()->GetTolerance(aTolX, aTolY);
      const std::pair<int, int>          aCells  = this->getCellsCount(aSt->NbNodes() - 3);
      const std::pair<double, double>&   aDelta  = this->getRangeSplitter().GetDelta();
      const std::pair<double, double>&   aTolUV  = this->getRangeSplitter().GetToleranceUV();
      const std::pair<double, double>&   aRangeU = this->getRangeSplitter().GetRangeU();
      const std::pair<double, double>&   aRangeV = this->getRangeSplitter().GetRangeV();
      std::printf("OCCTPOST PARAMS defl=%.9f tol=(%.17e,%.17e) cells=(%d,%d) delta=(%.17e,%.17e) "
                  "toluv=(%.17e,%.17e) rangeU=(%.9f,%.9f) rangeV=(%.9f,%.9f)\n",
                  this->getDFace()->GetDeflection(),
                  aTolX,
                  aTolY,
                  aCells.first,
                  aCells.second,
                  aDelta.first,
                  aDelta.second,
                  aTolUV.first,
                  aTolUV.second,
                  aRangeU.first,
                  aRangeU.second,
                  aRangeV.first,
                  aRangeV.second);
    }
    {
      const Handle(IMeshData::ListOfPnt2d) aSurfaceNodes =
        this->getRangeSplitter().GenerateSurfaceNodes(this->getParameters());
      int aIn = 0;
      const int aTotal = aSurfaceNodes.IsNull() ? 0 : aSurfaceNodes->Size();
      if (!aSurfaceNodes.IsNull())
      {
        IMeshData::ListOfPnt2d::Iterator aNodesIt(*aSurfaceNodes);
        for (; aNodesIt.More(); aNodesIt.Next())
        {
          if (this->getClassifier()->Perform(aNodesIt.Value()) == TopAbs_IN)
          {
            ++aIn;
          }
        }
      }
      std::printf("OCCTPOST surface_nodes total=%d in=%d\n", aTotal, aIn);
    }
    Base::postProcessMesh(theMesher, theRange);
    dumpAlgoState("POST_INSERT", theMesher);
  }

  void dumpAlgoState(const char* theTag, BRepMesh_Delaun& theMesher) const
  {
    const occ::handle<BRepMesh_DataStructureOfDelaun>& aSt = this->getStructure();
    const occ::handle<BRepMesh_DataStructureOfDelaun>& aRes = theMesher.Result();
    int aFrontier = 0, aFixed = 0, aFree = 0, aOther = 0;
    for (int e = 1; e <= aRes->NbLinks(); ++e)
    {
      switch (aRes->GetLink(e).Movability())
      {
        case BRepMesh_Frontier: ++aFrontier; break;
        case BRepMesh_Fixed:    ++aFixed;    break;
        case BRepMesh_Free:     ++aFree;     break;
        default:                ++aOther;    break;
      }
    }
    std::printf("OCCTPOST %s struct_nodes=%d struct_links=%d nodes=%d links=%d domain=%d "
                "frontier=%d fixed=%d free=%d other=%d\n",
                theTag,
                aSt.IsNull() ? -1 : aSt->NbNodes(),
                aSt.IsNull() ? -1 : aSt->NbLinks(),
                aRes->NbNodes(),
                aRes->NbLinks(),
                aRes->ElementsOfDomain().Extent(),
                aFrontier,
                aFixed,
                aFree,
                aOther);
    IMeshData::IteratorOfMapOfInteger aTriIt(aRes->ElementsOfDomain());
    for (; aTriIt.More(); aTriIt.Next())
    {
      int aNodes[3];
      aRes->ElementNodes(aRes->GetElement(aTriIt.Key()), aNodes);
      std::printf("OCCTPOST %s T %d %d %d %d\n",
                  theTag,
                  aTriIt.Key(),
                  aNodes[0],
                  aNodes[1],
                  aNodes[2]);
    }
    // TEMP T-105b: the full element sequence, deleted slots included, so the port's
    // `add_element` trace (`DELADDT`) can be diffed call by call.
    for (int i = 1; i <= aRes->NbElements(); ++i)
    {
      const BRepMesh_Triangle& anElt = aRes->GetElement(i);
      int                      aN[3];
      aRes->ElementNodes(anElt, aN);
      std::printf("OCCTEL %s E %d %d %d %d %d\n",
                  theTag,
                  i,
                  aN[0],
                  aN[1],
                  aN[2],
                  static_cast<int>(anElt.Movability()));
    }
  }

  bool initDataStructure() override
  {
    typedef BRepMesh_DelaunayNodeInsertionMeshAlgo<Splitter, BRepMesh_DelaunayBaseMeshAlgo> Base;
    const bool isOk = Base::initDataStructure();

    const occ::handle<BRepMesh_DataStructureOfDelaun>& aSt = this->getStructure();
    if (!aSt.IsNull())
    {
      std::printf("OCCTSTRUCT nodes=%d links=%d\n", aSt->NbNodes(), aSt->NbLinks());
      for (int i = 1; i <= aSt->NbNodes(); ++i)
      {
        const BRepMesh_Vertex& aV = aSt->GetNode(i);
        std::printf("N %d %.17e %.17e %d\n",
                    i,
                    aV.Coord().X(),
                    aV.Coord().Y(),
                    static_cast<int>(aV.Movability()));
      }
      for (int i = 1; i <= aSt->NbLinks(); ++i)
      {
        const BRepMesh_Edge& aL = aSt->GetLink(i);
        std::printf("L %d %d %d %d\n",
                    i,
                    aL.FirstNode(),
                    aL.LastNode(),
                    static_cast<int>(aL.Movability()));
      }
    }
    return isOk;
  }
};

// TEMP T-119: `BRepMesh_MeshAlgoFactory::GetAlgo` (`cxx:100-101`) pairs
// `GeomAbs_SurfaceOfRevolution` with `DeflectionControlMeshAlgo<
// BRepMesh_BoundaryParamsRangeSplitter>` (`cxx:42-46`). The DC pass is where that
// face's node count is decided, so hook `optimizeMesh` to dump the structure it
// receives (`POST_INSERT`, = base surface-node insertion) and the one it leaves
// (`POST_DC`, = the mesh `Compute` returns).
template <class Splitter>
class ZZStructDumpDCAlgo
    : public BRepMesh_DelaunayDeflectionControlMeshAlgo<Splitter, BRepMesh_DelaunayBaseMeshAlgo>
{
  typedef BRepMesh_DelaunayDeflectionControlMeshAlgo<Splitter, BRepMesh_DelaunayBaseMeshAlgo> BaseDC;

protected:
  void postProcessMesh(BRepMesh_Delaun& theMesher, const Message_ProgressRange& theRange) override
  {
    // TEMP POST_BASE: state right after the base mesher ran and before any
    // surface node is registered (`BRepMesh_DelaunayNodeInsertionMeshAlgo::
    // postProcessMesh` -> `InsertionBaseClass::postProcessMesh`).
    this->dumpState("POST_BASE", theMesher);
    BaseDC::postProcessMesh(theMesher, theRange);
  }

  void optimizeMesh(BRepMesh_Delaun& theMesher, const Message_ProgressRange& theRange) override
  {
    this->dumpState("POST_INSERT", theMesher);
    // TEMP DCPASS: verbatim copy of `BaseDC::optimizeMesh`'s loop with a per-pass
    // print, so the port's `DELSTAGE after_vertex_loop` chain can be lined up
    // pass by pass against OCCT.
    occ::handle<NCollection_IncAllocator> aTmpAlloc =
      new NCollection_IncAllocator(IMeshData::MEMORY_BLOCK_SIZE_HUGE);

    this->mySqMinSize = this->getParameters().MinSize * this->getParameters().MinSize;
    this->myCouplesMap =
      new IMeshData::MapOfOrientedEdges(3 * this->getStructure()->ElementsOfDomain().Extent(),
                                        aTmpAlloc);
    this->myControlNodes = new IMeshData::ListOfPnt2d(aTmpAlloc);
    this->myCircles      = &theMesher.Circles();

    const int aIterationsNb = 11;
    bool      isInserted    = true;
    Message_ProgressScope aPS(theRange, "Iteration", aIterationsNb);
    for (int aPass = 1; aPass <= aIterationsNb && isInserted && !this->myIsAllDegenerated; ++aPass)
    {
      if (!aPS.More())
      {
        return;
      }
      this->myMaxSqDeflection  = -1.;
      this->myIsAllDegenerated = true;
      this->myControlNodes->Clear();

      const int aDomainBefore = this->getStructure()->ElementsOfDomain().Extent();
      if (aDomainBefore < 1)
      {
        break;
      }
      IMeshData::IteratorOfMapOfInteger aTriangleIt(this->getStructure()->ElementsOfDomain());
      for (; aTriangleIt.More(); aTriangleIt.Next())
      {
        const BRepMesh_Triangle& aTriangle = this->getStructure()->GetElement(aTriangleIt.Key());
        this->splitTriangleGeometry(aTriangle);
      }

      std::printf("OCCTDCPASS pass=%d domain_before=%d ctrl=%d maxsq=%.17g alldeg=%d\n",
                  aPass,
                  aDomainBefore,
                  (int)this->myControlNodes->Size(),
                  this->myMaxSqDeflection,
                  this->myIsAllDegenerated ? 1 : 0);
      std::fflush(stdout);
      // TEMP CTRLDUMP: the exact control node list this pass hands to
      // `insertNodes`, in face-basis coordinates.
      {
        int aIdx = 0;
        for (IMeshData::ListOfPnt2d::Iterator aCtrlIt(*this->myControlNodes); aCtrlIt.More();
             aCtrlIt.Next())
        {
          std::printf("OCCTCTRL P%d %d %.17g %.17g\n",
                      aPass,
                      ++aIdx,
                      aCtrlIt.Value().X(),
                      aCtrlIt.Value().Y());
        }
        std::fflush(stdout);
      }

      isInserted = this->insertNodes(this->myControlNodes, theMesher, aPS.Next());

      std::printf("OCCTDCPASS pass=%d domain_after=%d nodes=%d inserted=%d\n",
                  aPass,
                  this->getStructure()->ElementsOfDomain().Extent(),
                  this->getStructure()->NbNodes(),
                  isInserted ? 1 : 0);
      std::fflush(stdout);
      {
        char aTag[16];
        std::snprintf(aTag, sizeof(aTag), "P%d", aPass);
        this->dumpState(aTag, theMesher);
      }
    }

    this->myCouplesMap.Nullify();
    this->myControlNodes.Nullify();

    if (!(this->myMaxSqDeflection < 0.))
    {
      this->getDFace()->SetDeflection(std::sqrt(this->myMaxSqDeflection));
    }
    this->dumpState("POST_DC", theMesher);
  }

  void dumpState(const char* theTag, BRepMesh_Delaun& theMesher) const
  {
    const occ::handle<BRepMesh_DataStructureOfDelaun>& aSt  = this->getStructure();
    const occ::handle<BRepMesh_DataStructureOfDelaun>& aRes = theMesher.Result();
    int aFrontier = 0, aFixed = 0, aFree = 0, aOther = 0;
    for (int e = 1; e <= aRes->NbLinks(); ++e)
    {
      switch (aRes->GetLink(e).Movability())
      {
        case BRepMesh_Frontier: ++aFrontier; break;
        case BRepMesh_Fixed:    ++aFixed;    break;
        case BRepMesh_Free:     ++aFree;     break;
        default:                ++aOther;    break;
      }
    }
    std::printf("OCCTDC %s struct_nodes=%d struct_links=%d nodes=%d links=%d domain=%d "
                "frontier=%d fixed=%d free=%d other=%d\n",
                theTag,
                aSt.IsNull() ? -1 : aSt->NbNodes(),
                aSt.IsNull() ? -1 : aSt->NbLinks(),
                aRes->NbNodes(),
                aRes->NbLinks(),
                aRes->ElementsOfDomain().Extent(),
                aFrontier,
                aFixed,
                aFree,
                aOther);
    IMeshData::IteratorOfMapOfInteger aTriIt(aRes->ElementsOfDomain());
    for (; aTriIt.More(); aTriIt.Next())
    {
      int aNodes[3];
      aRes->ElementNodes(aRes->GetElement(aTriIt.Key()), aNodes);
        std::printf("OCCTDC %s T %d %d %d %d\n",
                    theTag,
                    aTriIt.Key(),
                    aNodes[0],
                    aNodes[1],
                    aNodes[2]);
    }
    // TEMP T-119: node coordinates of the mesher's structure, so the port's
    // registered node list can be compared entry by entry (same ids).
    for (int i = 1; i <= aRes->NbNodes(); ++i)
    {
      const BRepMesh_Vertex& aV = aRes->GetNode(i);
      std::printf("OCCTDCN %s %d %.17g %.17g %d\n",
                  theTag,
                  i,
                  aV.Coord().X(),
                  aV.Coord().Y(),
                  static_cast<int>(aV.Movability()));
    }
    for (int i = 1; i <= aRes->NbElements(); ++i)
    {
      const BRepMesh_Triangle& anElt = aRes->GetElement(i);
      int                      aN[3];
      aRes->ElementNodes(anElt, aN);
      std::printf("OCCTEL %s E %d %d %d %d %d\n",
                  theTag,
                  i,
                  aN[0],
                  aN[1],
                  aN[2],
                  static_cast<int>(anElt.Movability()));
    }
  }
};

// TEMP T-120 (`--srfnodes <faceIdx>`): prints, for one discrete face, the surface
// nodes the factory-selected splitter produces and the verdict
// `BRepMesh_Classifier::Perform` gives each of them, i.e. exactly the filter
// `BRepMesh_DelaunayNodeInsertionMeshAlgo::insertNodes` applies
// (`BRepMesh_DelaunayNodeInsertionMeshAlgo.hxx:115-123`). Ground truth for the
// port's `list_surface_nodes` (`node_insertion.rs`) drop (Sphere.iges: 323
// generated, 255 IN, 68 ON). The classifier filter is identical in the DC variant,
// so the non-DC node-insertion algo is used here.
template <class Splitter>
class ZZSrfNodeCountAlgo
    : public BRepMesh_DelaunayNodeInsertionMeshAlgo<Splitter, BRepMesh_DelaunayBaseMeshAlgo>
{
  typedef BRepMesh_DelaunayNodeInsertionMeshAlgo<Splitter, BRepMesh_DelaunayBaseMeshAlgo> Base;

protected:
  void postProcessMesh(BRepMesh_Delaun& theMesher, const Message_ProgressRange& theRange) override
  {
    const Handle(IMeshData::ListOfPnt2d) aNodes =
      this->getRangeSplitter().GenerateSurfaceNodes(this->getParameters());
    int aIn = 0, aOut = 0, aOn = 0;
    const int aTotal = aNodes.IsNull() ? 0 : aNodes->Size();
    if (!aNodes.IsNull())
    {
      for (IMeshData::ListOfPnt2d::Iterator anIt(*aNodes); anIt.More(); anIt.Next())
      {
        switch (this->getClassifier()->Perform(anIt.Value()))
        {
          case TopAbs_IN:
            ++aIn;
            break;
          case TopAbs_OUT:
            ++aOut;
            std::printf("OCCTSRF OUT u=%.17g v=%.17g\n", anIt.Value().X(), anIt.Value().Y());
            break;
          default:
            ++aOn;
            break;
        }
      }
    }
    std::printf("OCCTSRF face total=%d in=%d out=%d on=%d stype=%d\n",
                aTotal,
                aIn,
                aOut,
                aOn,
                static_cast<int>(this->getDFace()->GetSurface()->GetType()));
    Base::postProcessMesh(theMesher, theRange);
  }
};

// --- TEMP T-133 (`--lens <faceIdx> <u0> <u1> <v0> <v1>`) ---------------------------
// Replays `BRepMesh_DefaultRangeSplitter::computeLengthU` / `computeLengthV`
// (`BRepMesh_DefaultRangeSplitter.cxx:144-198`) on one discrete face with an
// explicitly given range and prints every `BRepAdaptor_Surface::D0` sample at full
// precision plus the two lengths. The range is passed in because it is already the
// same on both sides (`--srparams` beats the port's NURBSGRID print), so any sample
// difference can only come from the surface evaluator, and the length difference
// from there into `computeDelta` (`cxx:129-137`).
int runLenSamples(const TopoDS_Shape& theShape, const int theFaceIdx,
                         const double u0, const double u1, const double v0, const double v1)
{
  occ::handle<BRepMesh_Context> aCtx = new BRepMesh_Context;
  aCtx->SetShape(theShape);
  IMeshTools_Parameters aParams;
  aParams.InParallel           = false;
  aParams.Relative             = false;
  aCtx->ChangeParameters()     = aParams;
  aCtx->ChangeParameters().CleanModel = false;
  if (!aCtx->BuildModel())
  {
    std::printf("LENS built=0\n");
    return 0;
  }
  const occ::handle<IMeshData_Model>& aModel = aCtx->GetModel();
  if (aModel.IsNull() || theFaceIdx < 0 || theFaceIdx >= aModel->FacesNb())
  {
    std::printf("LENS bad face\n");
    return 0;
  }
  const occ::handle<BRepAdaptor_Surface>& aSurf = aModel->GetFace(theFaceIdx)->GetSurface();
  std::printf("LENS f=%d type=%d ru=(%.17g,%.17g) rv=(%.17g,%.17g)\n",
              theFaceIdx,
              static_cast<int>(aSurf->GetType()),
              u0,
              u1,
              v0,
              v1);
  {
    double longv = 0.0;
    gp_Pnt P11, P12, P21, P22, P31, P32;
    double dv     = 0.05 * (v1 - v0);
    double dfuave = 0.5 * (u1 + u0);
    double dfvcur;
    aSurf->D0(u0, v0, P11);
    aSurf->D0(dfuave, v0, P21);
    aSurf->D0(u1, v0, P31);
    std::printf("LENS V k=0 (%.17g,%.17g,%.17g) (%.17g,%.17g,%.17g) (%.17g,%.17g,%.17g)\n",
                P11.X(), P11.Y(), P11.Z(), P21.X(), P21.Y(), P21.Z(), P31.X(), P31.Y(), P31.Z());
    double aV = v0 + dv;
    for (int i1 = 1; i1 <= 20; i1++, aV += dv)
    {
      aSurf->D0(u0, aV, P12);
      aSurf->D0(dfuave, aV, P22);
      aSurf->D0(u1, aV, P32);
      longv += (P11.Distance(P12) + P21.Distance(P22) + P31.Distance(P32));
      P11 = P12;
      P21 = P22;
      P31 = P32;
      std::printf("LENS V k=%d (%.17g,%.17g,%.17g) (%.17g,%.17g,%.17g) (%.17g,%.17g,%.17g)\n",
                  i1, P11.X(), P11.Y(), P11.Z(), P21.X(), P21.Y(), P21.Z(), P31.X(), P31.Y(), P31.Z());
    }
    std::printf("LENS V longv=%.17g\n", longv / 3.);
  }
  {
    double longu = 0.0;
    gp_Pnt P11, P12, P21, P22, P31, P32;
    double du     = 0.05 * (u1 - u0);
    double dfvave = 0.5 * (v1 + v0);
    aSurf->D0(u0, v0, P11);
    aSurf->D0(u0, dfvave, P21);
    aSurf->D0(u0, v1, P31);
    std::printf("LENS U k=0 (%.17g,%.17g,%.17g) (%.17g,%.17g,%.17g) (%.17g,%.17g,%.17g)\n",
                P11.X(), P11.Y(), P11.Z(), P21.X(), P21.Y(), P21.Z(), P31.X(), P31.Y(), P31.Z());
    double aU = u0 + du;
    for (int i1 = 1; i1 <= 20; i1++, aU += du)
    {
      aSurf->D0(aU, v0, P12);
      aSurf->D0(aU, dfvave, P22);
      aSurf->D0(aU, v1, P32);
      longu += (P11.Distance(P12) + P21.Distance(P22) + P31.Distance(P32));
      P11 = P12;
      P21 = P22;
      P31 = P32;
      std::printf("LENS U k=%d (%.17g,%.17g,%.17g) (%.17g,%.17g,%.17g) (%.17g,%.17g,%.17g)\n",
                  i1, P11.X(), P11.Y(), P11.Z(), P21.X(), P21.Y(), P21.Z(), P31.X(), P31.Y(), P31.Z());
    }
    std::printf("LENS U longu=%.17g\n", longu / 3.);
  }
  // TEMP T-133: the same two replays through the raw `Geom_Surface` (the port's
  // `Surface::d0`, i.e. `Geom_BSplineSurface::EvalD0` / `BSplSLib::D0`), so the
  // cache arm of `GeomAdaptor_Surface::EvalD0` (`GeomAdaptor_Surface.cxx:1063-1073`
  // -> `BSplSLib_Cache::D0`) can be told apart from the raw evaluator.
  {
    TopLoc_Location                 aLoc;
    const occ::handle<Geom_Surface> aRaw =
      BRep_Tool::Surface(aModel->GetFace(theFaceIdx)->GetFace(), aLoc);
    const bool     isIdentity = aLoc.IsIdentity();
    const gp_Trsf& aTrsf      = aLoc.Transformation();
    auto           aD0        = [&](double u, double v) {
      gp_Pnt aP = aRaw->EvalD0(u, v);
      if (!isIdentity)
      {
        aP.Transform(aTrsf);
      }
      return aP;
    };
    {
      double       longv = 0.0;
      const double dv    = 0.05 * (v1 - v0);
      const double uave  = 0.5 * (u1 + u0);
      gp_Pnt       P11 = aD0(u0, v0), P21 = aD0(uave, v0), P31 = aD0(u1, v0);
      std::printf("LENRAW V k=0 (%.17g,%.17g,%.17g) (%.17g,%.17g,%.17g) (%.17g,%.17g,%.17g)\n",
                  P11.X(), P11.Y(), P11.Z(), P21.X(), P21.Y(), P21.Z(), P31.X(), P31.Y(), P31.Z());
      double aV = v0 + dv;
      for (int i1 = 1; i1 <= 20; i1++, aV += dv)
      {
        gp_Pnt P12 = aD0(u0, aV), P22 = aD0(uave, aV), P32 = aD0(u1, aV);
        longv += (P11.Distance(P12) + P21.Distance(P22) + P31.Distance(P32));
        P11 = P12;
        P21 = P22;
        P31 = P32;
        std::printf("LENRAW V k=%d (%.17g,%.17g,%.17g) (%.17g,%.17g,%.17g) (%.17g,%.17g,%.17g)\n",
                    i1, P11.X(), P11.Y(), P11.Z(), P21.X(), P21.Y(), P21.Z(), P31.X(), P31.Y(), P31.Z());
      }
      std::printf("LENRAW V longv=%.17g\n", longv / 3.);
    }
    {
      double       longu = 0.0;
      const double du    = 0.05 * (u1 - u0);
      const double vave  = 0.5 * (v1 + v0);
      gp_Pnt       P11 = aD0(u0, v0), P21 = aD0(u0, vave), P31 = aD0(u0, v1);
      std::printf("LENRAW U k=0 (%.17g,%.17g,%.17g) (%.17g,%.17g,%.17g) (%.17g,%.17g,%.17g)\n",
                  P11.X(), P11.Y(), P11.Z(), P21.X(), P21.Y(), P21.Z(), P31.X(), P31.Y(), P31.Z());
      double aU = u0 + du;
      for (int i1 = 1; i1 <= 20; i1++, aU += du)
      {
        gp_Pnt P12 = aD0(aU, v0), P22 = aD0(aU, vave), P32 = aD0(aU, v1);
        longu += (P11.Distance(P12) + P21.Distance(P22) + P31.Distance(P32));
        P11 = P12;
        P21 = P22;
        P31 = P32;
        std::printf("LENRAW U k=%d (%.17g,%.17g,%.17g) (%.17g,%.17g,%.17g) (%.17g,%.17g,%.17g)\n",
                    i1, P11.X(), P11.Y(), P11.Z(), P21.X(), P21.Y(), P21.Z(), P31.X(), P31.Y(), P31.Z());
      }
      std::printf("LENRAW U longu=%.17g\n", longu / 3.);
    }
  }
  return 0;
}

int runStructDump(const TopoDS_Shape& theShape, const int theFaceIdx,
                         const bool theFactoryAlgo, const bool theSrfNodes,
                         const double theDefl, const double theAngle)
{
  Bnd_Box aShapeBox;
  BRepBndLib::Add(theShape, aShapeBox, false);
  double aBox[6] = {0, 0, 0, 0, 0, 0};
  if (!aShapeBox.IsVoid())
  {
    aShapeBox.Get(aBox[0], aBox[1], aBox[2], aBox[3], aBox[4], aBox[5]);
  }
  const double aMaxComp =
    (std::max)((std::max)(aBox[3] - aBox[0], aBox[4] - aBox[1]), aBox[5] - aBox[2]);

  IMeshTools_Parameters aParams;
  // TEMP T-121: the port derives its deflection from the *reader's* own bbox, so
  // the derived value here (`aMaxComp * 0.004`) differs from the port's in the
  // 8th digit, which is enough to change the NURBS/SoR grain filter and therefore
  // the whole surface-node grid. When the caller passes the port's exact numbers,
  // use them so the two DC-pass structures are comparable at all.
  aParams.Deflection    = (theDefl > 0.0) ? theDefl : aMaxComp * 0.001 * 4.0;
  aParams.Angle         = (theAngle > 0.0) ? theAngle : 20.0 * M_PI / 180.0;
  aParams.InParallel    = false;
  aParams.Relative      = false;
  aParams.MinSize       = (std::max)(IMeshTools_Parameters::RelMinSize() * aParams.Deflection,
                                     Precision::Confusion());
  aParams.AdjustMinSize = false;
  // TEMP T-119: `BRepMesh_IncrementalMesh::initParameters` (`hxx:79-107`) fills
  // `DeflectionInterior`/`AngleInterior` when they are left at their sentinel
  // (`-1`). Skipping it feeds the NURBS/SoR splitter a negative interior angle and
  // deflection, which defeats its grain filter and inflates the generated node
  // grid (measured: a 95x95 sphere grid instead of the real one).
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
  std::printf("OCCTSTRUCT defl=%.9f angle=%.9f angint=%.9f minsize=%.9f built=%d discret=%d healed=%d "
              "preproc=%d faces=%d\n",
              aParams.Deflection,
              aParams.Angle,
              aParams.AngleInterior,
              aParams.MinSize,
              aBuilt ? 1 : 0,
              aDiscret ? 1 : 0,
              aHealed ? 1 : 0,
              aPreproc ? 1 : 0,
              aModel.IsNull() ? -1 : aModel->FacesNb());
  if (aModel.IsNull() || theFaceIdx < 0 || theFaceIdx >= aModel->FacesNb())
  {
    return 0;
  }

  const IMeshData::IFaceHandle& aDFace = aModel->GetFace(theFaceIdx);
  Bnd_Box                       aFaceBox;
  BRepBndLib::Add(aDFace->GetFace(), aFaceBox, false);
  double v[6] = {0, 0, 0, 0, 0, 0};
  if (!aFaceBox.IsVoid())
  {
    aFaceBox.Get(v[0], v[1], v[2], v[3], v[4], v[5]);
  }
  std::printf("OCCTSTRUCT f=%d type=%d wires=%d bbox=(%.6f,%.6f,%.6f)-(%.6f,%.6f,%.6f)\n",
              theFaceIdx,
              static_cast<int>(aDFace->GetSurface()->GetType()),
              aDFace->WiresNb(),
              v[0],
              v[1],
              v[2],
              v[3],
              v[4],
              v[5]);
  for (int w = 0; w < aDFace->WiresNb(); ++w)
  {
    const IMeshData::IWireHandle& aDWire = aDFace->GetWire(w);
    std::printf("OCCTSTRUCT W %d edges=%d\n", w, aDWire->EdgesNb());
    for (int e = 0; e < aDWire->EdgesNb(); ++e)
    {
      const IMeshData::IEdgeHandle&   aDEdge  = aDWire->GetEdge(e);
      const IMeshData::ListOfInteger& aList   = aDEdge->GetPCurves(aDFace.get());
      for (IMeshData::ListOfInteger::Iterator it(aList); it.More(); it.Next())
      {
        const IMeshData::IPCurveHandle& aPCurve = aDEdge->GetPCurve(it.Value());
        std::printf("OCCTSTRUCT E w=%d e=%d pcId=%d ori=%d fwd=%d n=%d first=(%.9f,%.9f) "
                    "last=(%.9f,%.9f)\n",
                    w,
                    e,
                    it.Value(),
                    static_cast<int>(aPCurve->GetOrientation()),
                    aPCurve->IsForward() ? 1 : 0,
                    aPCurve->ParametersNb(),
                    aPCurve->GetPoint(0).X(),
                    aPCurve->GetPoint(0).Y(),
                    aPCurve->GetPoint(aPCurve->ParametersNb() - 1).X(),
                    aPCurve->GetPoint(aPCurve->ParametersNb() - 1).Y());
        // TEMP EPC: full-precision pcurve points + parameters, so the port's
        // `EPROBE uv=/tt=` list can be compared entry by entry.
        for (int k = 0; k < aPCurve->ParametersNb(); ++k)
        {
          const gp_Pnt2d aP = aPCurve->GetPoint(k);
          std::printf("OCCTEPC w=%d e=%d pcId=%d k=%d u=%.17g v=%.17g t=%.17g\n",
                      w,
                      e,
                      it.Value(),
                      k,
                      aP.X(),
                      aP.Y(),
                      aPCurve->GetParameter(k));
        }
      }
      // TEMP PC2D: the face's own pcurve of this edge
      // (`BRep_Tool::CurveOnSurface`), i.e. what `BRepMesh_ModelBuilder`
      // copies into `IMeshData::IPCurve`, with its dynamic type and a few
      // full-precision evaluations.
      {
        TopLoc_Location aPcLoc;
        (void)aPcLoc;
        double aPcF = 0., aPcL = 0.;
        const occ::handle<Geom2d_Curve>& aPc2d =
          BRep_Tool::CurveOnSurface(aDEdge->GetEdge(), aDFace->GetFace(), aPcF, aPcL);
        if (!aPc2d.IsNull())
        {
          const double aPcMid = 0.5 * (aPcF + aPcL);
          const gp_Pnt2d aPa = aPc2d->Value(aPcF);
          const gp_Pnt2d aPb = aPc2d->Value(aPcMid);
          const gp_Pnt2d aPc = aPc2d->Value(aPcL);
          std::printf("OCCTPC2D w=%d e=%d type=%s f=%.17g l=%.17g pf=(%.17g,%.17g) "
                      "pm=(%.17g,%.17g) pl=(%.17g,%.17g)\n",
                      w,
                      e,
                      aPc2d->DynamicType()->Name(),
                      aPcF,
                      aPcL,
                      aPa.X(),
                      aPa.Y(),
                      aPb.X(),
                      aPb.Y(),
                      aPc.X(),
                      aPc.Y());
        }
        else
        {
          std::printf("OCCTPC2D w=%d e=%d type=NULL\n", w, e);
        }
      }
      // TEMP PC2DORI: the pcurve of THIS edge on the face taken through both
      // edge orientations, i.e. exactly the two objects
      // `BRepMesh_ModelPreProcessor::splitEdge` (`cxx:213-219`) picks between:
      // `aE.Orientation(FORWARD); aPC1 = BRep_Tool::CurveOnSurface(aE, aF, af, al)`
      // and the same with REVERSED. Printing the type, the own range and the
      // value at first/mid/last shows which of the two carries the +-2*pi
      // (seam) branch, i.e. how `splitEdge` can insert a point on the opposite
      // branch of the same physical seam line.
      {
        const TopoDS_Edge aERaw = aDEdge->GetEdge();
        for (int aOriIt = 0; aOriIt < 2; ++aOriIt)
        {
          TopoDS_Edge aEOri = aERaw;
          aEOri.Orientation(aOriIt == 0 ? TopAbs_FORWARD : TopAbs_REVERSED);
          double aOF = 0., aOL = 0.;
          const occ::handle<Geom2d_Curve>& aOC =
            BRep_Tool::CurveOnSurface(aEOri, aDFace->GetFace(), aOF, aOL);
          if (aOC.IsNull())
          {
            std::printf("OCCTPC2DORI w=%d e=%d ori=%d type=NULL\n", w, e, aOriIt);
            continue;
          }
          const double   aOFirst = aOC->FirstParameter();
          const double   aOLast  = aOC->LastParameter();
          const gp_Pnt2d aOA     = aOC->Value(aOFirst);
          const gp_Pnt2d aOB     = aOC->Value(0.5 * (aOFirst + aOLast));
          const gp_Pnt2d aOCc    = aOC->Value(aOLast);
          std::printf("OCCTPC2DORI w=%d e=%d ori=%d type=%s range=%.17g %.17g evf=%.17g evl=%.17g "
                      "pf=(%.17g,%.17g) pm=(%.17g,%.17g) pl=(%.17g,%.17g)\n",
                      w,
                      e,
                      aOriIt,
                      aOC->DynamicType()->Name(),
                      aOFirst,
                      aOLast,
                      aOF,
                      aOL,
                      aOA.X(),
                      aOA.Y(),
                      aOB.X(),
                      aOB.Y(),
                      aOCc.X(),
                      aOCc.Y());
        }
        std::printf("OCCTPC2DFLAGS w=%d e=%d sameparam=%d samerange=%d pcnb=%d\n",
                    w,
                    e,
                    BRep_Tool::SameParameter(aERaw) ? 1 : 0,
                    BRep_Tool::SameRange(aERaw) ? 1 : 0,
                    aDEdge->PCurvesNb());
      }
    }
  }

  occ::handle<IMeshTools_MeshAlgo> aAlgo;
  if (theSrfNodes)
  {
    // TEMP T-120: factory-faithful splitter per `BRepMesh_MeshAlgoFactory::GetAlgo`
    // (`cxx:60-118`), node-insertion algo (the classifier filter under test is the
    // same in `DelaunayNodeInsertionMeshAlgo` and `DelaunayDeflectionControlMeshAlgo`).
    switch (aDFace->GetSurface()->GetType())
    {
      case GeomAbs_Plane:
        aAlgo = new ZZSrfNodeCountAlgo<BRepMesh_DefaultRangeSplitter>();
        break;
      case GeomAbs_Sphere:
        aAlgo = new ZZSrfNodeCountAlgo<BRepMesh_SphereRangeSplitter>();
        break;
      case GeomAbs_Cylinder:
        aAlgo = new ZZSrfNodeCountAlgo<BRepMesh_CylinderRangeSplitter>();
        break;
      case GeomAbs_Cone: aAlgo = new ZZSrfNodeCountAlgo<BRepMesh_ConeRangeSplitter>(); break;
      case GeomAbs_Torus: aAlgo = new ZZSrfNodeCountAlgo<BRepMesh_TorusRangeSplitter>(); break;
      case GeomAbs_SurfaceOfRevolution:
        aAlgo = new ZZSrfNodeCountAlgo<BRepMesh_BoundaryParamsRangeSplitter>();
        break;
      case GeomAbs_SurfaceOfExtrusion:
        aAlgo = new ZZSrfNodeCountAlgo<BRepMesh_ExtrusionRangeSplitter>();
        break;
      case GeomAbs_BezierSurface:
      case GeomAbs_BSplineSurface:
        aAlgo = new ZZSrfNodeCountAlgo<BRepMesh_NURBSRangeSplitter>();
        break;
      default: aAlgo = new ZZSrfNodeCountAlgo<BRepMesh_UndefinedRangeSplitter>(); break;
    }
    aAlgo->Perform(aDFace, aParams, Message_ProgressRange());
    return 0;
  }
  if (theFactoryAlgo)
  {
    // TEMP T-119: mirror `BRepMesh_MeshAlgoFactory::GetAlgo` (`cxx:60-118`) with
    // the probe's parameters (`EnableControlSurfaceDeflectionAllSurfaces` false,
    // `InternalVerticesMode` true). Only the DC-capable pairs are listed; the
    // `BaseMeshAlgo::Type` fallback for plane/cylinder is not installed here
    // because this probe always wants a dumpable structure.
    switch (aDFace->GetSurface()->GetType())
    {
      case GeomAbs_Sphere:
        aAlgo = new ZZStructDumpDCAlgo<BRepMesh_SphereRangeSplitter>();
        break;
      case GeomAbs_Cone: aAlgo = new ZZStructDumpDCAlgo<BRepMesh_ConeRangeSplitter>(); break;
      case GeomAbs_Torus: aAlgo = new ZZStructDumpDCAlgo<BRepMesh_TorusRangeSplitter>(); break;
      case GeomAbs_SurfaceOfRevolution:
        aAlgo = new ZZStructDumpDCAlgo<BRepMesh_BoundaryParamsRangeSplitter>();
        break;
      case GeomAbs_SurfaceOfExtrusion:
        aAlgo = new ZZStructDumpDCAlgo<BRepMesh_ExtrusionRangeSplitter>();
        break;
      case GeomAbs_BezierSurface:
      case GeomAbs_BSplineSurface:
        aAlgo = new ZZStructDumpDCAlgo<BRepMesh_NURBSRangeSplitter>();
        break;
      case GeomAbs_Cylinder:
        aAlgo = new ZZStructDumpDCAlgo<BRepMesh_CylinderRangeSplitter>();
        break;
      case GeomAbs_Plane:
        aAlgo = new ZZStructDumpDCAlgo<BRepMesh_DefaultRangeSplitter>();
        break;
      default: aAlgo = new ZZStructDumpDCAlgo<BRepMesh_UndefinedRangeSplitter>(); break;
    }
    aAlgo->Perform(aDFace, aParams, Message_ProgressRange());
    return 0;
  }

  switch (aDFace->GetSurface()->GetType())
  {
    case GeomAbs_Plane: aAlgo = new ZZStructDumpAlgo<BRepMesh_DefaultRangeSplitter>(); break;
    case GeomAbs_Sphere: aAlgo = new ZZStructDumpAlgo<BRepMesh_SphereRangeSplitter>(); break;
    case GeomAbs_Cylinder:
      aAlgo = new ZZStructDumpAlgo<BRepMesh_CylinderRangeSplitter>();
      break;
    case GeomAbs_Cone: aAlgo = new ZZStructDumpAlgo<BRepMesh_ConeRangeSplitter>(); break;
    case GeomAbs_Torus: aAlgo = new ZZStructDumpAlgo<BRepMesh_TorusRangeSplitter>(); break;
    default: aAlgo = new ZZStructDumpAlgo<BRepMesh_DefaultRangeSplitter>(); break;
  }
  aAlgo->Perform(aDFace, aParams, Message_ProgressRange());
  return 0;
}
