// Real-pipeline boundary UV dump (`--boundary <dir>`).
#include "occt_probe_common.h"

int runBoundary(const TopoDS_Shape& theShape, const std::string& theDir)
{
  Bnd_Box aShapeBox;
  BRepBndLib::Add(theShape, aShapeBox, false);
  double aBox[6] = {0, 0, 0, 0, 0, 0};
  if (!aShapeBox.IsVoid())
  {
    aShapeBox.Get(aBox[0], aBox[1], aBox[2], aBox[3], aBox[4], aBox[5]);
  }
  const double aMaxComp =
    std::max(std::max(aBox[3] - aBox[0], aBox[4] - aBox[1]), aBox[5] - aBox[2]);

  IMeshTools_Parameters aParams;
  // The port's a3n00 gate meshes with `prs3d_get_deflection(shape, 0.1)`
  // (`maxComp * 0.001 * 4`) and 20 degrees, so match that exactly.
  aParams.Deflection    = aMaxComp * 0.001 * 4.0;
  aParams.Angle         = 20.0 * M_PI / 180.0;
  aParams.InParallel    = false;
  aParams.Relative      = false;
  // Mirror `BRepMesh_IncrementalMesh::initParameters` (`hxx:91-96`): the real
  // mesher runs with `MinSize = max(RelMinSize * min(defl, deflInterior), Confusion)`.
  // The earlier `Precision::Confusion()` made `BRepMesh_CurveTessellator::init` hand
  // `GCPnts_TangentialDeflection` a `MinimumLength` of 1e-7, so every arc got the full
  // 10-degree (half-angle) step and the dump was NOT the structure the mesher builds.
  aParams.MinSize       = (std::max)(IMeshTools_Parameters::RelMinSize() * aParams.Deflection,
                                     Precision::Confusion());
  aParams.AdjustMinSize = false;

  occ::handle<BRepMesh_Context> aCtx = new BRepMesh_Context;
  aCtx->SetShape(theShape);
  aCtx->ChangeParameters()            = aParams;
  aCtx->ChangeParameters().CleanModel = false;
  const bool aBuilt     = aCtx->BuildModel();
  const bool aDiscret   = aBuilt && aCtx->DiscretizeEdges();
  const bool aHealed    = aDiscret && aCtx->HealModel();
  const bool aPreproc   = aHealed && aCtx->PreProcessModel();
  const occ::handle<IMeshData_Model>& aModel = aCtx->GetModel();
  std::cout << "BOUNDARY defl=" << aParams.Deflection << " angle=" << aParams.Angle
            << " built=" << (aBuilt ? 1 : 0) << " discret=" << (aDiscret ? 1 : 0)
            << " healed=" << (aHealed ? 1 : 0) << " preproc=" << (aPreproc ? 1 : 0)
            << " model_null=" << (aModel.IsNull() ? 1 : 0) << std::endl;
  if (aModel.IsNull())
  {
    return 0;
  }

  const int aFaces[] = {0,   208, 209, 169, 174, 189, 192, 193, 194,
                        195, 196, 197, 198, 199, 204, 206, 210, 212};
  (void)aFaces;

  std::ostringstream anOut;
  for (int f = 0; f < aModel->FacesNb(); ++f)
  {
    const IMeshData::IFaceHandle& aDFace = aModel->GetFace(f);
    Bnd_Box                       aFaceBox;
    BRepBndLib::Add(aDFace->GetFace(), aFaceBox, false);
    double v[6] = {0, 0, 0, 0, 0, 0};
    if (!aFaceBox.IsVoid())
    {
      aFaceBox.Get(v[0], v[1], v[2], v[3], v[4], v[5]);
    }
    int aPts = 0;
    for (int w = 0; w < aDFace->WiresNb(); ++w)
    {
      const IMeshData::IWireHandle& aDWire = aDFace->GetWire(w);
      for (int e = 0; e < aDWire->EdgesNb(); ++e)
      {
        const IMeshData::ListOfInteger& aList = aDWire->GetEdge(e)->GetPCurves(aDFace.get());
        for (IMeshData::ListOfInteger::Iterator it(aList); it.More(); it.Next())
        {
          aPts += aDWire->GetEdge(e)->GetPCurve(it.Value())->ParametersNb();
        }
      }
    }
    anOut << std::fixed << "BOUND f=" << f << " type=" << (int)aDFace->GetSurface()->GetType()
          << " wires=" << aDFace->WiresNb() << " pts=" << aPts << " bbox=(" << v[0] << "," << v[1]
          << "," << v[2] << ")-(" << v[3] << "," << v[4] << "," << v[5] << ")\n";
    for (int w = 0; w < aDFace->WiresNb(); ++w)
    {
      const IMeshData::IWireHandle& aDWire = aDFace->GetWire(w);
      anOut << "W " << w << " edges=" << aDWire->EdgesNb()
            << " selfint=" << (aDWire->IsSet(IMeshData_SelfIntersectingWire) ? 1 : 0)
            << " open=" << (aDWire->IsSet(IMeshData_OpenWire) ? 1 : 0) << "\n";
      for (int e = 0; e < aDWire->EdgesNb(); ++e)
      {
        const IMeshData::IEdgeHandle& aDEdge = aDWire->GetEdge(e);
        if (true)
        {
          // `BRepMesh_CurveTessellator::init` inputs (`cxx:60-120`): the edge's stored
          // deflections and the curve type decide how many points
          // `GCPnts_TangentialDeflection` emits, so dump them next to the points.
          BRepAdaptor_Curve aBRepCurve(aDEdge->GetEdge());
          anOut << "  EINFO w=" << w << " e=" << e << " ctype=" << (int)aBRepCurve.GetType()
                << " defl=" << aDEdge->GetDeflection()
                << " ang=" << aDEdge->GetAngularDeflection()
                << " sameparam=" << (aDEdge->GetSameParam() ? 1 : 0)
                << " samerange=" << (aDEdge->GetSameRange() ? 1 : 0)
                << " degen=" << (aDEdge->GetDegenerated() ? 1 : 0)
                << " free=" << (aDEdge->IsFree() ? 1 : 0)
                << " first=" << aBRepCurve.FirstParameter()
                << " last=" << aBRepCurve.LastParameter();
          if (aBRepCurve.GetType() == GeomAbs_Circle)
          {
            anOut << " radius=" << aBRepCurve.Circle().Radius();
          }
          anOut << "\n";
        }
        // Exactly `BRepMesh_BaseMeshAlgo.cxx:91-96`: the pcurves of this edge on
        // this face, in the structure's own order.
        const IMeshData::ListOfInteger& aListOfPCurves = aDEdge->GetPCurves(aDFace.get());
        for (IMeshData::ListOfInteger::Iterator aPCurveIt(aListOfPCurves); aPCurveIt.More();
             aPCurveIt.Next())
        {
          const IMeshData::IPCurveHandle& aPCurve = aDEdge->GetPCurve(aPCurveIt.Value());
          anOut << "E w=" << w << " e=" << e << " pcId=" << aPCurveIt.Value()
                << " ori=" << static_cast<int>(aPCurve->GetOrientation())
                << " fwd=" << (aPCurve->IsForward() ? 1 : 0)
                << " n=" << aPCurve->ParametersNb() << "\n";
          for (int p = 0; p < aPCurve->ParametersNb(); ++p)
          {
            const gp_Pnt2d& aP2 = aPCurve->GetPoint(p);
            char            aBuf[128];
            std::snprintf(aBuf, sizeof(aBuf), "P %.17e %.17e\n", aP2.X(), aP2.Y());
            anOut << aBuf;
          }
        }
      }
    }
  }

  const std::string aPath = theDir + "/boundary.txt";
  std::ofstream     aFile(aPath.c_str());
  if (!aFile)
  {
    std::cout << "BOUNDARY cannot write " << aPath << std::endl;
    return 0;
  }
  aFile << anOut.str();
  aFile.close();
  std::cout << "BOUNDARY wrote " << aPath << std::endl;
  return 0;
}
