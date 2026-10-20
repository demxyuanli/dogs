// `argv[2]` modes: `--uv` / `--uvsum` / `--nofix` (node-insertion UV inputs).
#include "occt_probe_common.h"

std::optional<int> runUvModes(const ProbeCtx& theCtx)
{
  const int           argc = theCtx.argc;
  char** const        argv = theCtx.argv;
  const TopoDS_Shape& aShape = theCtx.shape;

  // T-69 oracle: the (u,v) point set that BRepMesh_NodeInsertionMeshAlgo::
  // collectWirePoints feeds into BRepMesh_DefaultRangeSplitter::AddPoint for one
  // face, plus the min/max range those points produce. This is the exact input
  // whose spread decides `computeLengthU/V` and therefore `IsValid`
  // (`BRepMesh_DefaultRangeSplitter.cxx:35-41` / `:72-80`).
  if (argc > 2 && (std::string(argv[2]) == "--uv" || std::string(argv[2]) == "--uvsum"
                   || std::string(argv[2]) == "--nofix"))
  {
    int  aWant = -1;
    bool doSum = false;
    double aDefl = 0.1, anAngle = 0.349066;
    for (int i = 2; i < argc; ++i)
    {
      const std::string a = argv[i];
      if (a == "--uv" && i + 1 < argc)
      {
        aWant = std::atoi(argv[i + 1]);
      }
      else if (a == "--uvdef" && i + 2 < argc)
      {
        aDefl   = std::atof(argv[i + 1]);
        anAngle = std::atof(argv[i + 2]);
        i += 2;
      }
      else if (a == "--uvsum")
      {
        doSum = true;
      }
    }
    IMeshTools_Parameters aParams;
    aParams.Deflection    = aDefl;
    aParams.Angle         = anAngle;
    aParams.InParallel    = false;
    aParams.Relative      = false;
    aParams.MinSize       = Precision::Confusion();
    aParams.AdjustMinSize = false;

    occ::handle<BRepMesh_ModelBuilder> aBuilder = new BRepMesh_ModelBuilder;
    const occ::handle<IMeshData_Model> aModel   = aBuilder->Perform(aShape, aParams);
    if (aModel.IsNull())
    {
      std::cout << "UV model-null\n";
      return 0;
    }
    occ::handle<BRepMesh_EdgeDiscret> anEdgeDiscret = new BRepMesh_EdgeDiscret;
    const bool aDiscretDone = anEdgeDiscret->Perform(aModel, aParams, Message_ProgressRange());
    // T-99: this branch only builds the discrete model and discretizes edges - it
    // never triangulates, so every face's `Poly_Triangulation` is null and the
    // summary's `triangles=` field came out 0 (measured). Triangulate the shape
    // here so that field reports OCCT's real per-face triangle count. The
    // deflection is taken from argv ("--uvsum <d>", default 0.1) so the caller can
    // match whatever deflection the ground-truth OBJ was produced with.
    {
      double aUvDefl = 0.1;
      for (int i = 2; i + 1 < argc; ++i)
      {
        if (std::string(argv[i]) == "--uvsum")
        {
          aUvDefl = std::atof(argv[i + 1]);
        }
      }
      BRepMesh_IncrementalMesh aSumMesher(aShape, aUvDefl, false, 0.5);
    }
    {
      int aPCurvesTotal = 0, aFreeEdges = 0, aFilled = 0, aEmptyPCurves = 0;
      for (int e = 0; e < aModel->EdgesNb(); ++e)
      {
        const IMeshData::IEdgeHandle& aE = aModel->GetEdge(e);
        if (aE->IsFree())
        {
          ++aFreeEdges;
        }
        for (int p = 0; p < aE->PCurvesNb(); ++p)
        {
          ++aPCurvesTotal;
          const int aNb = aE->GetPCurve(p)->ParametersNb();
          if (aNb > 0)
          {
            ++aFilled;
          }
          else
          {
            ++aEmptyPCurves;
          }
        }
      }
      std::cout << "UV model edges=" << aModel->EdgesNb() << " faces=" << aModel->FacesNb()
                << " pcurves=" << aPCurvesTotal << " free=" << aFreeEdges << " filled=" << aFilled
                << " empty=" << aEmptyPCurves << " discret_done=" << (aDiscretDone ? 1 : 0)
                << " maxsize=" << aModel->GetMaxSize() << "\n";
    }

    if (aWant >= 0)
    {
    const IMeshData::IFaceHandle& aDFace = aModel->GetFace(aWant);
    const TopoDS_Face&            aFace  = aDFace->GetFace();
    {
      Bnd_Box aBox;
      BRepBndLib::Add(aFace, aBox);
      double v[6] = {0, 0, 0, 0, 0, 0};
      if (!aBox.IsVoid())
      {
        aBox.Get(v[0], v[1], v[2], v[3], v[4], v[5]);
      }
      std::cout << std::fixed << std::setprecision(6) << "UV face=" << aWant << " wires=" << aDFace->WiresNb()
                << " u=[" << aDFace->GetSurface()->FirstUParameter() << ","
                << aDFace->GetSurface()->LastUParameter() << "] v=["
                << aDFace->GetSurface()->FirstVParameter() << ","
                << aDFace->GetSurface()->LastVParameter() << "]"
                << " uper=" << (aDFace->GetSurface()->IsUPeriodic() ? 1 : 0)
                << " vper=" << (aDFace->GetSurface()->IsVPeriodic() ? 1 : 0)
                << " bbox=(" << v[0] << "," << v[1] << "," << v[2] << ")-(" << v[3] << "," << v[4] << ","
                << v[5] << ")\n";
    }

    double uMin = 1e100, uMax = -1e100, vMin = 1e100, vMax = -1e100;
    int    nPts = 0, onUMin = 0, onUMax = 0, onVMin = 0, onVMax = 0;
    for (int aWireIt = 0; aWireIt < aDFace->WiresNb(); ++aWireIt)
    {
      const IMeshData::IWireHandle& aDWire = aDFace->GetWire(aWireIt);
      for (int aEdgeIt = 0; aEdgeIt < aDWire->EdgesNb(); ++aEdgeIt)
      {
        const IMeshData::IEdgeHandle& aDEdge = aDWire->GetEdge(aEdgeIt);
        const TopAbs_Orientation      anOri = aDWire->GetEdgeOrientation(aEdgeIt);
        const IMeshData::IPCurveHandle& aPCurve = aDEdge->GetPCurve(aDFace.get(), anOri);
        std::cout << "UV  w" << aWireIt << " e" << aEdgeIt << " ori=" << (int)anOri
                  << " prefix=" << (anOri == TopAbs_REVERSED ? "REV" : "FWD")
                  << " nat=" << (aPCurve->IsForward() ? 1 : 0)
                  << " npc=" << aPCurve->ParametersNb();
        if (aPCurve->ParametersNb() > 0)
        {
          const gp_Pnt2d aF = aPCurve->GetPoint(0);
          const gp_Pnt2d aL = aPCurve->GetPoint(aPCurve->ParametersNb() - 1);
          std::cout << " first=(" << aF.X() << "," << aF.Y() << ") last=(" << aL.X() << "," << aL.Y() << ")"
                    << " p0=" << aPCurve->GetParameter(0) << " pn="
                    << aPCurve->GetParameter(aPCurve->ParametersNb() - 1);
        }
        std::cout << "\n";
        for (int aPt = 0; aPt < aPCurve->ParametersNb(); ++aPt)
        {
          const gp_Pnt2d& aP = aPCurve->GetPoint(aPt);
          std::cout << std::setprecision(17) << "UV   p" << aPt << " uv=(" << aP.X() << "," << aP.Y()
                    << ") t=" << aPCurve->GetParameter(aPt) << "\n";
          ++nPts;
          uMin = std::min(uMin, aP.X());
          uMax = std::max(uMax, aP.X());
          vMin = std::min(vMin, aP.Y());
          vMax = std::max(vMax, aP.Y());
        }
      }
    }
    // Count how many collected points sit exactly on the extremal coordinates:
    // a zero-width direction shows every point on both the min and the max.
    for (int aWireIt = 0; aWireIt < aDFace->WiresNb(); ++aWireIt)
    {
      const IMeshData::IWireHandle& aDWire = aDFace->GetWire(aWireIt);
      for (int aEdgeIt = 0; aEdgeIt < aDWire->EdgesNb(); ++aEdgeIt)
      {
        const IMeshData::IPCurveHandle& aPCurve =
          aDWire->GetEdge(aEdgeIt)->GetPCurve(aDFace.get(), aDWire->GetEdgeOrientation(aEdgeIt));
        for (int aPt = 0; aPt < aPCurve->ParametersNb(); ++aPt)
        {
          const gp_Pnt2d& aP = aPCurve->GetPoint(aPt);
          if (aP.X() == uMin) ++onUMin;
          if (aP.X() == uMax) ++onUMax;
          if (aP.Y() == vMin) ++onVMin;
          if (aP.Y() == vMax) ++onVMax;
        }
      }
    }
    // Reproduce `BRepMesh_NodeInsertionMeshAlgo::collectWirePoints`
    // (`BRepMesh_NodeInsertionMeshAlgo.hxx:142-184`): every pcurve sample except
    // the traversal's last one.
    double rUMin = 1e100, rUMax = -1e100, rVMin = 1e100, rVMax = -1e100;
    for (int aWireIt = 0; aWireIt < aDFace->WiresNb(); ++aWireIt)
    {
      const IMeshData::IWireHandle& aDWire = aDFace->GetWire(aWireIt);
      if (aDWire->IsSet(IMeshData_SelfIntersectingWire)
          || (aDWire->IsSet(IMeshData_OpenWire) && aWireIt != 0))
      {
        continue;
      }
      for (int aEdgeIt = 0; aEdgeIt < aDWire->EdgesNb(); ++aEdgeIt)
      {
        const IMeshData::IPCurveHandle& aPCurve =
          aDWire->GetEdge(aEdgeIt)->GetPCurve(aDFace.get(), aDWire->GetEdgeOrientation(aEdgeIt));
        int aIt, aEnd, aInc;
        if (aPCurve->IsForward())
        {
          aEnd = aPCurve->ParametersNb() - 1;
          aIt  = (std::min)(0, aEnd);
          aInc = 1;
        }
        else
        {
          aIt  = aPCurve->ParametersNb() - 1;
          aEnd = (std::min)(0, aIt);
          aInc = -1;
        }
        for (; aIt != aEnd; aIt += aInc)
        {
          const gp_Pnt2d& aP = aPCurve->GetPoint(aIt);
          rUMin = std::min(rUMin, aP.X());
          rUMax = std::max(rUMax, aP.X());
          rVMin = std::min(rVMin, aP.Y());
          rVMax = std::max(rVMax, aP.Y());
        }
      }
    }
    std::cout << std::setprecision(17) << "UV RANGE pts=" << nPts << " u=[" << uMin << "," << uMax
              << "] v=[" << vMin << "," << vMax << "] dU=" << (uMax - uMin) << " dV=" << (vMax - vMin)
              << " onUMin=" << onUMin << " onUMax=" << onUMax << " onVMin=" << onVMin
              << " onVMax=" << onVMax << "\n";
    std::cout << std::setprecision(17) << "UV ADDEDRANGE u=[" << rUMin << "," << rUMax
              << "] v=[" << rVMin << "," << rVMax << "] dU=" << (rUMax - rUMin)
              << " dV=" << (rVMax - rVMin) << "\n";
    // The predicate `AdjustRange` uses to decide `IsValid`
    // (`BRepMesh_DefaultRangeSplitter.cxx:72-80`) on the added-point range.
    const double aDu = 0.05 * (rUMax - rUMin);
    const double aDv = 0.05 * (rVMax - rVMin);
    const occ::handle<BRepAdaptor_Surface>& aSurf = aDFace->GetSurface();
    double      aLenU = 0., aLenV = 0.;
    gp_Pnt      P11, P12, P21, P22, P31, P32;
    const double aVave = 0.5 * (rVMax + rVMin);
    aSurf->D0(rUMin, rVMin, P11);
    aSurf->D0(rUMin, aVave, P21);
    aSurf->D0(rUMin, rVMax, P31);
    for (int i1 = 1; i1 <= 20; ++i1)
    {
      const double aU = rUMin + aDu * i1;
      aSurf->D0(aU, rVMin, P12);
      aSurf->D0(aU, aVave, P22);
      aSurf->D0(aU, rVMax, P32);
      aLenU += (P11.Distance(P12) + P21.Distance(P22) + P31.Distance(P32));
      P11 = P12;
      P21 = P22;
      P31 = P32;
    }
    const double aUave = 0.5 * (rUMax + rUMin);
    aSurf->D0(rUMin, rVMin, P11);
    aSurf->D0(aUave, rVMin, P21);
    aSurf->D0(rUMax, rVMin, P31);
    for (int i1 = 1; i1 <= 20; ++i1)
    {
      const double aV = rVMin + aDv * i1;
      aSurf->D0(rUMin, aV, P12);
      aSurf->D0(aUave, aV, P22);
      aSurf->D0(rUMax, aV, P32);
      aLenV += (P11.Distance(P12) + P21.Distance(P22) + P31.Distance(P32));
      P11 = P12;
      P21 = P22;
      P31 = P32;
    }
    aLenU /= 3.;
    aLenV /= 3.;
    std::cout << std::setprecision(17) << "UV LEN lenU=" << aLenU << " lenV=" << aLenV
              << " pconf=" << Precision::PConfusion() << " valid="
              << ((aLenU > Precision::PConfusion() && aLenV > Precision::PConfusion()) ? 1 : 0) << "\n";
    if (!doSum)
    {
      return 0;
    }
    }  // end single-face (`--uv <k>`) block

    // Summary: one line per face, so the whole shape can be paired by geometry.
    for (int f = 0; f < aModel->FacesNb(); ++f)
    {
      const IMeshData::IFaceHandle& aDF = aModel->GetFace(f);
      Bnd_Box                       aBox;
      BRepBndLib::Add(aDF->GetFace(), aBox);
      double w[6] = {0, 0, 0, 0, 0, 0};
      if (!aBox.IsVoid())
      {
        aBox.Get(w[0], w[1], w[2], w[3], w[4], w[5]);
      }
      // T-99: this face's own triangulation, so a per-face triangle count can be
      // attached to the SAME face index as this dump. `--mesh` cannot be used for
      // that: it enumerates `TopExp_Explorer(aShape, TopAbs_FACE)` while this loop
      // enumerates `IMeshData_Model::GetFace(f)`, and the two orderings were
      // measured to agree at only 1 of 226 indices. Taking the count HERE removes
      // every cross-dump index assumption. `BRep_Tool::Triangulation` is the
      // OCCT API for exactly this (see `--mesh` above); a face may legitimately
      // carry no triangulation, which prints as 0.
      int aFaceTris = 0;
      {
        TopLoc_Location                        aLoc;
        const occ::handle<Poly_Triangulation>& aTri =
          BRep_Tool::Triangulation(TopoDS::Face(aDF->GetFace()), aLoc);
        if (!aTri.IsNull())
        {
          aFaceTris = aTri->NbTriangles();
        }
      }
      double sUMin = 1e100, sUMax = -1e100, sVMin = 1e100, sVMax = -1e100;
      for (int aWireIt = 0; aWireIt < aDF->WiresNb(); ++aWireIt)
      {
        const IMeshData::IWireHandle& aDWire = aDF->GetWire(aWireIt);
        if (aDWire->IsSet(IMeshData_SelfIntersectingWire)
            || (aDWire->IsSet(IMeshData_OpenWire) && aWireIt != 0))
        {
          continue;
        }
        for (int aEdgeIt = 0; aEdgeIt < aDWire->EdgesNb(); ++aEdgeIt)
        {
          const IMeshData::IPCurveHandle& aPC =
            aDWire->GetEdge(aEdgeIt)->GetPCurve(aDF.get(), aDWire->GetEdgeOrientation(aEdgeIt));
          int aIt, aEnd, aInc;
          if (aPC->IsForward())
          {
            aEnd = aPC->ParametersNb() - 1;
            aIt  = (std::min)(0, aEnd);
            aInc = 1;
          }
          else
          {
            aIt  = aPC->ParametersNb() - 1;
            aEnd = (std::min)(0, aIt);
            aInc = -1;
          }
          for (; aIt != aEnd; aIt += aInc)
          {
            const gp_Pnt2d& aP = aPC->GetPoint(aIt);
            sUMin = std::min(sUMin, aP.X());
            sUMax = std::max(sUMax, aP.X());
            sVMin = std::min(sVMin, aP.Y());
            sVMax = std::max(sVMax, aP.Y());
          }
        }
      }
      const int    sType = (int)aDF->GetSurface()->GetType();
      const double sDU = sUMax - sUMin, sDV = sVMax - sVMin;
      // T-99: reproduce the range splitter's own intermediate values at the same
      // pipeline point the port reports (`adjust_range` -> `generate_surface_nodes`),
      // so the two sides can be diffed value by value. Every formula on this path was
      // verified line-for-line against the .cxx, so a divergence has to come from
      // these INPUTS. The splitter is built here exactly as the port does it:
      // Reset(face) + AddPoint(each pcurve sample) + AdjustRange().
      if (std::getenv("OCCT_TOPO_TRACE_SPLITTER") != nullptr)
      {
        BRepMesh_DefaultRangeSplitter aSplitter;
        aSplitter.Reset(aDF, aParams);
        for (int aWi = 0; aWi < aDF->WiresNb(); ++aWi)
        {
          const IMeshData::IWireHandle& aDW = aDF->GetWire(aWi);
          for (int aEi = 0; aEi < aDW->EdgesNb(); ++aEi)
          {
            const IMeshData::IPCurveHandle& aPC =
              aDW->GetEdge(aEi)->GetPCurve(aDF.get(), aDW->GetEdgeOrientation(aEi));
            for (int aPt = 0; aPt < aPC->ParametersNb(); ++aPt)
            {
              aSplitter.AddPoint(aPC->GetPoint(aPt));
            }
          }
        }
        aSplitter.AdjustRange();
        const std::pair<double, double>& aRu = aSplitter.GetRangeU();
        const std::pair<double, double>& aRv = aSplitter.GetRangeV();
        const std::pair<double, double>& aDe = aSplitter.GetDelta();
        const std::pair<double, double>& aTo = aSplitter.GetToleranceUV();
        std::cout << "SPLITTER f=" << f << " preset=1 ru=[" << std::setprecision(15)
                  << aRu.first << "," << aRu.second << "] rv=[" << aRv.first << "," << aRv.second
                  << "] delta=[" << aDe.first << "," << aDe.second << "]";
        std::cout << std::scientific << std::setprecision(6) << " tol=[" << aTo.first << ","
                  << aTo.second << "]";
        std::cout << std::fixed << " vnb=0 defl=" << aParams.Deflection;
        // T-99: the cells count this face would get. This is the decision number -
        // range, delta and every formula on the path already matched between the port
        // and OCCT (see specs/_a3n00_gap_analysis.md §9.343/§9.344/§9.346), so if the
        // cell counts agree the divergence must be downstream of it. Printed twice
        // because OCCT passes `aStructure->NbNodes()` (Deleted included) while the
        // port's `node_insertion.rs:362` excludes Deleted nodes.
        const std::pair<int, int> aCellsZero =
          BRepMesh_GeomTool::CellsCount(aDF->GetSurface(), 0, aParams.Deflection, &aSplitter);
        std::cout << " cells_vnb0=[" << aCellsZero.first << "," << aCellsZero.second << "]";
        std::cout << std::setprecision(6) << " bbox=(" << w[0] << "," << w[1] << "," << w[2]
                  << ")-(" << w[3] << "," << w[4] << "," << w[5] << ")\n";
      }
      const double sDu = 0.05 * sDU, sDv = 0.05 * sDV;
      const occ::handle<BRepAdaptor_Surface>& aS2 = aDF->GetSurface();
      double lU = 0., lV = 0.;
      gp_Pnt A11, A12, A21, A22, A31, A32;
      const double vA = 0.5 * (sVMax + sVMin);
      aS2->D0(sUMin, sVMin, A11);
      aS2->D0(sUMin, vA, A21);
      aS2->D0(sUMin, sVMax, A31);
      for (int i1 = 1; i1 <= 20; ++i1)
      {
        const double aU = sUMin + sDu * i1;
        aS2->D0(aU, sVMin, A12);
        aS2->D0(aU, vA, A22);
        aS2->D0(aU, sVMax, A32);
        lU += (A11.Distance(A12) + A21.Distance(A22) + A31.Distance(A32));
        A11 = A12;
        A21 = A22;
        A31 = A32;
      }
      const double uA = 0.5 * (sUMax + sUMin);
      aS2->D0(sUMin, sVMin, A11);
      aS2->D0(uA, sVMin, A21);
      aS2->D0(sUMax, sVMin, A31);
      for (int i1 = 1; i1 <= 20; ++i1)
      {
        const double aV = sVMin + sDv * i1;
        aS2->D0(sUMin, aV, A12);
        aS2->D0(uA, aV, A22);
        aS2->D0(sUMax, aV, A32);
        lV += (A11.Distance(A12) + A21.Distance(A22) + A31.Distance(A32));
        A11 = A12;
        A21 = A22;
        A31 = A32;
      }
      lU /= 3.;
      lV /= 3.;
      const bool sValid = (lU > Precision::PConfusion() && lV > Precision::PConfusion());
      std::cout << std::fixed << std::setprecision(9) << "UVSUM f=" << f << " type=" << sType
                << " wires=" << aDF->WiresNb() << " triangles=" << aFaceTris
                << " bbox=(" << w[0] << "," << w[1] << "," << w[2]
                << ")-(" << w[3] << "," << w[4] << "," << w[5] << ")";
      std::cout << std::setprecision(17) << " urange=[" << sUMin << "," << sUMax << "] vrange=["
                << sVMin << "," << sVMax << "] dU=" << sDU << " dV=" << sDV << " lenU=" << lU
                << " lenV=" << lV << " valid=" << (sValid ? 1 : 0) << "\n";
    }
    return 0;
  }
  return std::nullopt;
}
