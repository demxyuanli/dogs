// OCCT ground-truth probe (T-69/T-80/T-87/T-05 oracle).
// Reads a STEP file with OCCT 8.0.0 and prints what the port must reproduce:
//   * BRepCheck_Analyzer verdict
//   * BRepGProp::VolumeProperties (the OCCT volume, analytic by default)
//   * per-face surface type, orientation, composed normal and outward sign
//   * cylinder/spline frames, so the port's reconstruction can be compared
#include <BRepAdaptor_Surface.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <Poly_Triangulation.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <cstdlib>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <BRepGProp_Face.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <STEPControl_Reader.hxx>
#include <ShapeProcess.hxx>
#include <DE_ShapeFixParameters.hxx>
#include <StepData_StepModel.hxx>
#include <BRepGProp_Domain.hxx>
#include <BRepGProp_Vinert.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BOPAlgo_PaveFiller.hxx>
#include <BOPDS_DS.hxx>
#include <BOPAlgo_Builder.hxx>
#include <BOPTools_AlgoTools.hxx>
#include <BOPTools_CoupleOfShape.hxx>
#include <IntTools_Context.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BOPTools_AlgoTools3D.hxx>
#include <BOPTools_AlgoTools2D.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <Geom_Plane.hxx>
#include <BRep_Tool.hxx>
#include <BRepTools.hxx>
#include <TopTools_ListOfShape.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <StepShape_FaceSurface.hxx>
#include <string>
#include <TopExp_Explorer.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopExp.hxx>
#include <TopoDS_Vertex.hxx>
#include <Geom_Curve.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

#include <cstdio>
#include <cstdlib>
#include <iostream>

int main(int argc, char** argv)
{
  if (argc < 2)
  {
    std::cerr << "usage: occt_probe <file.step> [--frames] [--nofix]\n";
    return 2;
  }
  const bool frames = (argc > 2);
  bool       nofix  = false;
  bool       faces0 = false;
  for (int i = 2; i < argc; ++i)
  {
    if (std::string(argv[i]) == "--nofix")
      nofix = true;
    if (std::string(argv[i]) == "--faces0")
      faces0 = true;
  }

  STEPControl_Reader aReader;
  if (nofix)
  {
    ShapeProcess::OperationsFlags aFlags;
    aFlags.reset();
    aReader.SetShapeProcessFlags(aFlags);
    std::cout << "SHAPE-PROCESS disabled\n";
  }
  if (faces0)
  {
    DE_ShapeFixParameters aFix;
    aFix.FixFaceOrientationMode = static_cast<DE_ShapeFixParameters::FixMode>(0); // NeedFix(0) => do not fix
    aReader.SetShapeFixParameters(aFix);
    std::cout << "SHP-FIX face-orientation disabled\n";
  }
  if (aReader.ReadFile(argv[1]) != IFSelect_RetDone)
  {
    std::cerr << "read failed: " << argv[1] << "\n";
    return 3;
  }
  if (nofix)
  {
    ShapeProcess::OperationsFlags aFlags;
    aFlags.reset();
    aReader.SetShapeProcessFlags(aFlags);
    DE_ShapeFixParameters aFix;
    aFix.FixShellOrientationMode = static_cast<DE_ShapeFixParameters::FixMode>(0);
    aFix.FixFaceOrientationMode = static_cast<DE_ShapeFixParameters::FixMode>(0);
    aFix.FixOrientationMode = static_cast<DE_ShapeFixParameters::FixMode>(0);
    aReader.SetShapeFixParameters(aFix);
    std::cout << "PROC-OFF (after read, before transfer)\n";
  }
  aReader.TransferRoots();
  const TopoDS_Shape aShape = aReader.OneShape();
  if (aShape.IsNull())
  {
    std::cerr << "null shape\n";
    return 4;
  }

  // T-54 oracle: OCCT's triangulation counts at a given deflection
  // (`BRepMesh_IncrementalMesh`), per face and in total. The port's
  // `export_data_obj` gate compares its own `v=`/`f=` counts, so this is
  // what a faithful constrained-Delaunay port must reproduce.
  if (argc > 2 && std::string(argv[2]) == "--mesh")
  {
    const double aDeflection = (argc > 3) ? std::atof(argv[3]) : 0.1;
    BRepMesh_IncrementalMesh aMesher(aShape, aDeflection);
    TopExp_Explorer            anEx(aShape, TopAbs_FACE);
    int                        anIdx = 0, aNodes = 0, aTris = 0;
    for (; anEx.More(); anEx.Next(), ++anIdx)
    {
      TopLoc_Location                     aLoc;
      const occ::handle<Poly_Triangulation>& aTri =
        BRep_Tool::Triangulation(TopoDS::Face(anEx.Current()), aLoc);
      const int aN = aTri.IsNull() ? 0 : aTri->NbNodes();
      const int aT = aTri.IsNull() ? 0 : aTri->NbTriangles();
      std::cout << "FACE " << anIdx << " nodes=" << aN << " triangles=" << aT << "\n";
      aNodes += aN;
      aTris += aT;
    }
    std::cout << "TOTAL faces=" << anIdx << " nodes=" << aNodes << " triangles=" << aTris
              << " deflection=" << aDeflection << "\n";
    return 0;
  }

  if (argc > 2 && std::string(argv[2]) == "--entities")
  {
    occ::handle<StepData_StepModel> aModel = aReader.StepModel();
    int aSameT = 0, aSameF = 0;
    for (int i = 1; i <= aModel->NbEntities(); ++i)
    {
      occ::handle<StepShape_FaceSurface> aFS =
        occ::down_cast<StepShape_FaceSurface>(aModel->Value(i));
      if (!aFS.IsNull())
      {
        if (aFS->SameSense()) ++aSameT; else ++aSameF;
        std::cout << "FS ent=" << i << " sameSense=" << (aFS->SameSense() ? 1 : 0) << "\n";
      }
    }
    std::cout << "FS total sameSense T=" << aSameT << " F=" << aSameF << "\n";
    return 0;
  }
  bool perface = false;
  for (int i = 1; i < argc; ++i) { if (std::string(argv[i]) == "--perface") perface = true; }
  if (perface)
  {
    GProp_GProps aP;
    BRepGProp::VolumeProperties(aShape, aP);
    const gp_Pnt aC(0.0, 0.0, 0.0);
    std::cout << "PERFACE total=" << aP.Mass() << " ref=(" << aC.X() << "," << aC.Y() << "," << aC.Z() << ") cm=(" << aP.CentreOfMass().X() << "," << aP.CentreOfMass().Y() << "," << aP.CentreOfMass().Z() << ")\n";
    int k = 0;
    double aSum = 0.0;
    for (TopExp_Explorer ex(aShape, TopAbs_FACE); ex.More(); ex.Next(), ++k)
    {
      const TopoDS_Face aF = TopoDS::Face(ex.Current());
      BRepGProp_Face   aBF;
      aBF.Load(aF);
      BRepGProp_Domain aBD;
      aBD.Init(aF);
      BRepGProp_Vinert aG;
      aG.SetLocation(aC);
      BRepAdaptor_Surface aBAS2(aF);
      aG.Perform(aBF, aBD);
      aSum += aG.Mass();
      std::cout << "PF face " << k << " type=" << (int)aBAS2.GetType() << " orient=" << (int)aF.Orientation() << " contrib=" << aG.Mass();
      if (aBAS2.GetType() == GeomAbs_Sphere) { const gp_Ax3 ax = aBAS2.Sphere().Position(); std::cout << " sphc=(" << ax.Location().X() << "," << ax.Location().Y() << "," << ax.Location().Z() << ")"; }
      std::cout << "\n";
    }
    std::cout << "PF sum=" << aSum << "\n";
    return 0;
  }
  if (argc > 2 && std::string(argv[2]) == "--ds")
  {
    auto report = [](const char* tag, double ax) {
      TopoDS_Shape a = BRepPrimAPI_MakeBox(gp_Pnt(ax, 0., 0.), 1., 1., 1.).Shape();
      TopoDS_Shape b = BRepPrimAPI_MakeBox(gp_Pnt(0.5, 0., 0.), 1., 1., 1.).Shape();
      (void)b;
      TopTools_ListOfShape aLS;
      aLS.Append(a);
      TopoDS_Shape b2 = BRepPrimAPI_MakeBox(gp_Pnt(ax + 3.0, 0., 0.), 1., 1., 1.).Shape();
      aLS.Append(b2);
      BOPAlgo_PaveFiller aPF;
      aPF.SetArguments(aLS);
      aPF.Perform();
      std::cout << "DS " << tag << " NbShapes=" << aPF.DS().NbShapes()
                << " NbSourceShapes=" << aPF.DS().NbSourceShapes() << "\n";
    };
    report("disjoint", 0.0);
    {
      TopoDS_Shape a = BRepPrimAPI_MakeBox(gp_Pnt(0., 0., 0.), 1., 1., 1.).Shape();
      TopoDS_Shape b = BRepPrimAPI_MakeBox(gp_Pnt(0.5, 0., 0.), 1., 1., 1.).Shape();
      TopTools_ListOfShape aLS; aLS.Append(a); aLS.Append(b);
      BOPAlgo_PaveFiller aPF; aPF.SetArguments(aLS); aPF.Perform();
      std::cout << "DS overlap NbShapes=" << aPF.DS().NbShapes()
                << " NbSourceShapes=" << aPF.DS().NbSourceShapes() << "\n";
    }
    return 0;
  }
  if (argc > 2 && std::string(argv[2]) == "--fuse")
  {
    TopoDS_Shape aBox = BRepPrimAPI_MakeBox(gp_Pnt(-1., -1., -1.), gp_Pnt(1., 1., 1.)).Shape();
    TopoDS_Shape aCyl = BRepPrimAPI_MakeCylinder(0.4, 2.0).Shape();
    BRepAlgoAPI_Fuse aFuse(aBox, aCyl);
    aFuse.Build();
    if (!aFuse.IsDone()) { std::cout << "FUSE failed\n"; return 1; }
    TopoDS_Shape aRes = aFuse.Shape();
    GProp_GProps aP; BRepGProp::VolumeProperties(aRes, aP);
    std::cout << "FUSE volume=" << aP.Mass() << " faces=";
    int nf = 0; for (TopExp_Explorer ex(aRes, TopAbs_FACE); ex.More(); ex.Next()) ++nf;
    std::cout << nf << "\n";
    const gp_Pnt aOrig(0., 0., 0.);
    int i = 0;
    for (TopExp_Explorer ex(aRes, TopAbs_FACE); ex.More(); ex.Next(), ++i)
    {
      const TopoDS_Face aF = TopoDS::Face(ex.Current());
      BRepAdaptor_Surface aBAS(aF);
      GProp_GProps aFA; BRepGProp::SurfaceProperties(aF, aFA);
      BRepGProp_Vinert aV; aV.SetLocation(aOrig);
      BRepGProp_Face aBF; aBF.Load(aF); BRepGProp_Domain aBD; aBD.Init(aF); aV.Perform(aBF, aBD);
      std::cout << "FFACE " << i << " type=" << (int)aBAS.GetType() << " orient=" << (int)aF.Orientation()
                << " area=" << aFA.Mass() << " contrib=" << aV.Mass() << "\n";
    }
    return 0;
  }
  if (argc > 2 && std::string(argv[2]) == "--dir2")
  {
    // Faithful re-implementation of BOPTools_AlgoTools.cxx's file-local
    // FindPointInFace (:2160-2231) so the final aDB of GetFaceDir (:2110-2152)
    // can be observed.
    struct FindRes { bool found; gp_Pnt pOut; gp_Dir db; };
    static bool zzVerbose = false;
    auto FindPointInFace = [](const TopoDS_Face& aF, const gp_Pnt& aP, gp_Dir& aDB, gp_Pnt& aPOut,
                              const occ::handle<IntTools_Context>& ctx,
                              GeomAPI_ProjectPointOnSurf& aProjPL, double aDt, double aTolE) -> FindRes {
      double aDTol = Precision::Angular();
      const double aPM = aP.XYZ().Modulus();
      if (aPM > 1000.) aDTol = 5.e-16 * aPM;
      int aNbItMax = 15;
      const double anEps = Precision::SquareConfusion();
      GeomAPI_ProjectPointOnSurf& aProj = ctx->ProjPS(aF);
      gp_Pnt aPS = aP;
      aProj.Perform(aPS);
      if (!aProj.IsDone()) { if (zzVerbose) std::cout << "      early#1 not done\n"; return {false, aPOut, aDB}; }
      aPS = aProj.NearestPoint();
      aProjPL.Perform(aPS);
      aPS = aProjPL.NearestPoint();
      aPS.SetXYZ(aPS.XYZ() + 2. * aTolE * aDB.XYZ());
      aProj.Perform(aPS);
      if (!aProj.IsDone()) { if (zzVerbose) std::cout << "      early#2 not done\n"; return {false, aPOut, aDB}; }
      aPS = aProj.NearestPoint();
      aProjPL.Perform(aPS);
      aPS = aProjPL.NearestPoint();
      bool bRet = false;
      double aDist = 0.;
      do {
        const gp_Pnt aP1(aPS.XYZ() + aDt * aDB.XYZ());
        aProj.Perform(aP1);
        if (!aProj.IsDone()) { if (zzVerbose) std::cout << "      early#3 not done\n"; return {false, aPOut, aDB}; }
        aPOut = aProj.NearestPoint();
        aDist = aProj.LowerDistance();
        aProjPL.Perform(aPOut);
        aPOut = aProjPL.NearestPoint();
        const gp_Vec aV(aPS, aPOut);
        if (zzVerbose)
          std::cout << "      it aDist=" << aDist << " |aV|=" << aV.Magnitude()
                    << " anEps=" << anEps << " dbNow=(" << aDB.X() << "," << aDB.Y() << "," << aDB.Z() << ")\n";
        if (aV.SquareMagnitude() < anEps) return {false, aPOut, aDB};
        aDB.SetXYZ(aV.XYZ());
        aPS = aPOut;
      } while (aDist > aDTol && --aNbItMax);
      bRet = aDist < aDTol;
      return {bRet, aPOut, aDB};
    };

    TopoDS_Shape aBox = BRepPrimAPI_MakeBox(gp_Pnt(-1., -1., -1.), gp_Pnt(1., 1., 1.)).Shape();
    TopoDS_Shape aCyl = BRepPrimAPI_MakeCylinder(0.4, 2.0).Shape();
    TopTools_ListOfShape aArgs;
    aArgs.Append(aBox);
    aArgs.Append(aCyl);
    BOPAlgo_Builder aGF;
    aGF.SetArguments(aArgs);
    aGF.Perform();
    TopoDS_Shape aRes = aGF.Shape();
    TopoDS_Edge aSec;
    for (TopExp_Explorer ex(aRes, TopAbs_EDGE); ex.More(); ex.Next())
    {
      const TopoDS_Edge e = TopoDS::Edge(ex.Current());
      BRepAdaptor_Curve ac(e);
      if (ac.GetType() != GeomAbs_Circle) continue;
      const gp_Circ c = ac.Circle();
      if (std::abs(c.Radius() - 0.4) < 1e-9 && std::abs(c.Location().Z() - 1.0) < 1e-9) { aSec = e; break; }
    }
    if (aSec.IsNull()) { std::cout << "no sect\n"; return 1; }
    NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher> aMEF;
    TopExp::MapShapesAndAncestors(aRes, TopAbs_EDGE, TopAbs_FACE, aMEF);
    const NCollection_List<TopoDS_Shape>& aLF = aMEF.FindFromKey(aSec);
    occ::handle<IntTools_Context> aCtx = new IntTools_Context;
    // current face = first cylinder
    TopoDS_Face f1;
    TopoDS_Edge e1;
    for (NCollection_List<TopoDS_Shape>::Iterator it(aLF); it.More(); it.Next())
    {
      const TopoDS_Face f = TopoDS::Face(it.Value());
      BRepAdaptor_Surface b(f);
      if (b.GetType() == GeomAbs_Cylinder) { f1 = f; BOPTools_AlgoTools::GetEdgeOnFace(aSec, f, e1); break; }
    }
    double t1, t2;
    BRep_Tool::Curve(e1, t1, t2);
    const double t = 0.5 * (t1 + t2);
    gp_Pnt aPx;
    BRepAdaptor_Curve ac1(e1);
    ac1.D0(t, aPx);
    gp_Vec tau1;
    BOPTools_AlgoTools2D::EdgeTangent(e1, t, tau1);
    const gp_Dir aDTgt(tau1);
    gp_Dir aDN1;
    BOPTools_AlgoTools3D::GetNormalToFaceOnEdge(e1, f1, t, aDN1, aCtx);
    gp_Dir aDBF(aDN1.XYZ().Crossed(aDTgt.XYZ()));
    // aProjPL: plane through aPx with normal aDTgt
    occ::handle<Geom_Plane> aPL = new Geom_Plane(aPx, aDTgt);
    double u1, u2, v1, v2;
    aPL->Bounds(u1, u2, v1, v2);
    GeomAPI_ProjectPointOnSurf aProjPL;
    aProjPL.Init(aPL, u1, u2, v1, v2);
    const double tolE = BRep_Tool::Tolerance(e1);
    // Faithful MinStep3D (BOPTools_AlgoTools.cxx:2235-2346).
    double dtMax = -1.;
    double dtMin = 5.e-6;
    for (NCollection_List<TopoDS_Shape>::Iterator it(aLF); it.More(); it.Next())
    {
      const TopoDS_Face f = TopoDS::Face(it.Value());
      dtMax = std::max(dtMax, 2. * (tolE + BRep_Tool::Tolerance(f)));
      BRepAdaptor_Surface aBAS(f);
      double aR = 0.;
      switch (aBAS.GetType())
      {
        case GeomAbs_Cylinder: aR = aBAS.Cylinder().Radius(); break;
        case GeomAbs_Cone: { gp_Lin aL(aBAS.Cone().Axis()); aR = aL.Distance(aPx); } break;
        case GeomAbs_Sphere: dtMin = std::max(dtMin, 5.e-4); aR = aBAS.Sphere().Radius(); break;
        case GeomAbs_Torus: aR = aBAS.Torus().MajorRadius(); break;
        default: dtMin = std::max(dtMin, 5.e-4); break;
      }
      if (aR > 100.)
      {
        constexpr double d = 10 * Precision::PConfusion();
        dtMin = std::max(dtMin, std::sqrt(d * d + 2 * d * aR));
      }
    }
    if (dtMax < dtMin) dtMax = dtMin;
    const double aDt = dtMax;
    gp_Pnt p1 = aPx;
    gp_Dir db1 = aDBF;
    zzVerbose = true;
    FindRes r1 = FindPointInFace(f1, aPx, db1, p1, aCtx, aProjPL, aDt, tolE);
    zzVerbose = false;
    if (r1.found) {
      aDBF = db1;
    } else {
      // GetFaceDir's fallback (BOPTools_AlgoTools.cxx:2139-2149)
      gp_Pnt aPx2;
      gp_Dir aDN2f;
      BOPTools_AlgoTools3D::GetApproxNormalToFaceOnEdge(e1, f1, t, aDt, aPx2, aDN2f, aCtx);
      aProjPL.Perform(aPx2);
      aPx2 = aProjPL.NearestPoint();
      const gp_Vec aVec(aPx, aPx2);
      aDBF.SetXYZ(aVec.XYZ());
      aDN1 = aDN2f;
    }
    const gp_Dir aDTF(aDN1.XYZ().Crossed(aDBF.XYZ()));
    BRepAdaptor_Surface b1(f1);
    std::cout << "CUR type=" << (int)b1.GetType() << " eo=" << (int)e1.Orientation()
              << " dt=" << aDt << " found1=" << (r1.found ? 1 : 0)
              << " dbf=(" << aDBF.X() << "," << aDBF.Y() << "," << aDBF.Z() << ")"
              << " dtf=(" << aDTF.X() << "," << aDTF.Y() << "," << aDTF.Z() << ")\n";
    for (NCollection_List<TopoDS_Shape>::Iterator it(aLF); it.More(); it.Next())
    {
      const TopoDS_Face f2 = TopoDS::Face(it.Value());
      if (f2.IsSame(f1)) continue;
      TopoDS_Edge e2;
      if (!BOPTools_AlgoTools::GetEdgeOnFace(aSec, f2, e2)) continue;
      const gp_Dir aDTgt2 = (e2.Orientation() == e1.Orientation()) ? aDTgt : aDTgt.Reversed();
      gp_Dir aDN2;
      BOPTools_AlgoTools3D::GetNormalToFaceOnEdge(e2, f2, t, aDN2, aCtx);
      gp_Dir db2(aDN2.XYZ().Crossed(aDTgt2.XYZ()));
      const gp_Dir db2_init = db2;
      gp_Pnt p2 = aPx;
      FindRes r2 = FindPointInFace(f2, aPx, db2, p2, aCtx, aProjPL, aDt, tolE);
      const double aAngle = aDBF.AngleWithRef(db2, aDTF) * 180. / M_PI;
      BRepAdaptor_Surface b2(f2);
      std::cout << "  CAND type=" << (int)b2.GetType() << " eo=" << (int)e2.Orientation()
                << " found=" << (r2.found ? 1 : 0)
                << " db_init=(" << db2_init.X() << "," << db2_init.Y() << "," << db2_init.Z() << ")"
                << " db_fin=(" << db2.X() << "," << db2.Y() << "," << db2.Z() << ")"
                << " angle=" << aAngle << "\n";
    }
    return 0;
  }
  if (argc > 2 && std::string(argv[2]) == "--dir")
  {
    TopoDS_Shape aBox = BRepPrimAPI_MakeBox(gp_Pnt(-1., -1., -1.), gp_Pnt(1., 1., 1.)).Shape();
    TopoDS_Shape aCyl = BRepPrimAPI_MakeCylinder(0.4, 2.0).Shape();
    TopTools_ListOfShape aArgs;
    aArgs.Append(aBox);
    aArgs.Append(aCyl);
    BOPAlgo_Builder aGF;
    aGF.SetArguments(aArgs);
    aGF.Perform();
    TopoDS_Shape aRes = aGF.Shape();
    TopoDS_Edge aSec;
    for (TopExp_Explorer ex(aRes, TopAbs_EDGE); ex.More(); ex.Next())
    {
      const TopoDS_Edge e = TopoDS::Edge(ex.Current());
      BRepAdaptor_Curve ac(e);
      if (ac.GetType() != GeomAbs_Circle) continue;
      const gp_Circ c = ac.Circle();
      if (std::abs(c.Radius() - 0.4) < 1e-9 && std::abs(c.Location().Z() - 1.0) < 1e-9) { aSec = e; break; }
    }
    if (aSec.IsNull()) { std::cout << "no sect\n"; return 1; }
    NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher> aMEF;
    TopExp::MapShapesAndAncestors(aRes, TopAbs_EDGE, TopAbs_FACE, aMEF);
    const NCollection_List<TopoDS_Shape>& aLF = aMEF.FindFromKey(aSec);
    occ::handle<IntTools_Context> aCtx = new IntTools_Context;
    int k = 0;
    for (NCollection_List<TopoDS_Shape>::Iterator it(aLF); it.More(); it.Next(), ++k)
    {
      const TopoDS_Face f = TopoDS::Face(it.Value());
      TopoDS_Edge eo;
      if (!BOPTools_AlgoTools::GetEdgeOnFace(aSec, f, eo)) continue;
      double t1, t2;
      BRep_Tool::Curve(eo, t1, t2);
      const double t = 0.5 * (t1 + t2);
      gp_Dir aD;
      BOPTools_AlgoTools3D::GetNormalToFaceOnEdge(eo, f, t, aD, aCtx);
      gp_Vec tau;
      BOPTools_AlgoTools2D::EdgeTangent(eo, t, tau);
      gp_Dir aDTgt(tau);
      const gp_Vec db = aD.XYZ().Crossed(aDTgt.XYZ());
      BRepAdaptor_Surface b(f);
      std::cout << "DIR f" << k << " type=" << (int)b.GetType() << " orient=" << (int)f.Orientation()
                << " eo=" << (int)eo.Orientation()
                << " dN=(" << aD.X() << "," << aD.Y() << "," << aD.Z() << ")"
                << " dTgt=(" << aDTgt.X() << "," << aDTgt.Y() << "," << aDTgt.Z() << ")"
                << " dB=(" << db.X() << "," << db.Y() << "," << db.Z() << ")\n";
    }
    return 0;
  }
  if (argc > 2 && std::string(argv[2]) == "--faceoff2")
  {
    // Same as --faceoff, but builds the candidate couples exactly the way
    // BOPAlgo_ShellSplitter::SplitBlock does: BOPTools_AlgoTools::GetEdgeOff
    // (BOPTools_AlgoTools.cxx:1099-1127) insists on the edge view whose
    // orientation is the OPPOSITE of the current one, and skips the candidate
    // otherwise. Compare the two modes to see whether the couple construction
    // decides the answer.
    TopoDS_Shape aBox = BRepPrimAPI_MakeBox(gp_Pnt(-1., -1., -1.), gp_Pnt(1., 1., 1.)).Shape();
    TopoDS_Shape aCyl = BRepPrimAPI_MakeCylinder(0.4, 2.0).Shape();
    TopTools_ListOfShape aArgs;
    aArgs.Append(aBox);
    aArgs.Append(aCyl);
    BOPAlgo_Builder aGF;
    aGF.SetArguments(aArgs);
    aGF.Perform();
    TopoDS_Shape aRes = aGF.Shape();
    TopoDS_Edge aSec;
    for (TopExp_Explorer ex(aRes, TopAbs_EDGE); ex.More(); ex.Next())
    {
      const TopoDS_Edge e = TopoDS::Edge(ex.Current());
      BRepAdaptor_Curve ac(e);
      if (ac.GetType() != GeomAbs_Circle) continue;
      const gp_Circ c = ac.Circle();
      if (std::abs(c.Radius() - 0.4) < 1e-9 && std::abs(c.Location().Z() - 1.0) < 1e-9) { aSec = e; break; }
    }
    if (aSec.IsNull()) { std::cout << "no section circle\n"; return 1; }
    NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher> aMEF;
    TopExp::MapShapesAndAncestors(aRes, TopAbs_EDGE, TopAbs_FACE, aMEF);
    const NCollection_List<TopoDS_Shape>& aLF = aMEF.FindFromKey(aSec);
    occ::handle<IntTools_Context> aCtx = new IntTools_Context;
    int k = 0;
    for (NCollection_List<TopoDS_Shape>::Iterator it(aLF); it.More(); it.Next(), ++k)
    {
      const TopoDS_Face f1 = TopoDS::Face(it.Value());
      TopoDS_Edge e1;
      if (!BOPTools_AlgoTools::GetEdgeOnFace(aSec, f1, e1)) continue;
      NCollection_List<BOPTools_CoupleOfShape> aLC;
      int skipped = 0, kept = 0;
      for (NCollection_List<TopoDS_Shape>::Iterator jt(aLF); jt.More(); jt.Next())
      {
        const TopoDS_Face f2 = TopoDS::Face(jt.Value());
        if (f2.IsSame(f1)) continue;
        TopoDS_Edge e2;
        if (!BOPTools_AlgoTools::GetEdgeOff(e1, f2, e2)) { ++skipped; continue; }
        BOPTools_CoupleOfShape cs;
        cs.SetShape1(e2);
        cs.SetShape2(f2);
        aLC.Append(cs);
        ++kept;
      }
      TopoDS_Face aOff;
      const bool done = BOPTools_AlgoTools::GetFaceOff(e1, f1, aLC, aOff, aCtx);
      BRepAdaptor_Surface b1(f1);
      if (!aOff.IsNull())
      {
        BRepAdaptor_Surface b2(aOff);
        std::cout << "O2 from f" << k << "(type=" << (int)b1.GetType() << " eo=" << (int)e1.Orientation()
                  << ") couples=" << kept << "/skip=" << skipped << " -> type=" << (int)b2.GetType()
                  << " orient=" << (int)aOff.Orientation() << " done=" << (done ? 1 : 0) << "\n";
      }
      else
      {
        std::cout << "O2 from f" << k << "(type=" << (int)b1.GetType() << " eo=" << (int)e1.Orientation()
                  << ") couples=" << kept << "/skip=" << skipped << " -> null done=" << (done ? 1 : 0) << "\n";
      }
    }
    return 0;
  }
  if (argc > 2 && std::string(argv[2]) == "--faceoff")
  {
    TopoDS_Shape aBox = BRepPrimAPI_MakeBox(gp_Pnt(-1., -1., -1.), gp_Pnt(1., 1., 1.)).Shape();
    TopoDS_Shape aCyl = BRepPrimAPI_MakeCylinder(0.4, 2.0).Shape();
    TopTools_ListOfShape aArgs;
    aArgs.Append(aBox);
    aArgs.Append(aCyl);
    BOPAlgo_Builder aGF;
    aGF.SetArguments(aArgs);
    aGF.Perform();
    if (aGF.HasErrors()) { std::cout << "GF failed\n"; return 1; }
    TopoDS_Shape aRes = aGF.Shape();
    // section circle: radius 0.4 at z=1
    TopoDS_Edge aSec;
    for (TopExp_Explorer ex(aRes, TopAbs_EDGE); ex.More(); ex.Next())
    {
      const TopoDS_Edge e = TopoDS::Edge(ex.Current());
      BRepAdaptor_Curve ac(e);
      if (ac.GetType() != GeomAbs_Circle) continue;
      const gp_Circ c = ac.Circle();
      if (std::abs(c.Radius() - 0.4) < 1e-9 && std::abs(c.Location().Z() - 1.0) < 1e-9)
      {
        aSec = e;
        break;
      }
    }
    if (aSec.IsNull()) { std::cout << "no section circle\n"; return 1; }
    NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher> aMEF;
    TopExp::MapShapesAndAncestors(aRes, TopAbs_EDGE, TopAbs_FACE, aMEF);
    if (!aMEF.Contains(aSec)) { std::cout << "edge not in map\n"; return 1; }
    const NCollection_List<TopoDS_Shape>& aLF = aMEF.FindFromKey(aSec);
    std::cout << "SEC eo=" << (int)aSec.Orientation() << " adj_faces=" << aLF.Extent() << "\n";
    NCollection_List<BOPTools_CoupleOfShape> aLCS;
    for (NCollection_List<TopoDS_Shape>::Iterator it(aLF); it.More(); it.Next())
    {
      const TopoDS_Face f = TopoDS::Face(it.Value());
      BRepAdaptor_Surface bas(f);
      TopoDS_Edge eo;
      const bool ok = BOPTools_AlgoTools::GetEdgeOnFace(aSec, f, eo);
      std::cout << "  face type=" << (int)bas.GetType() << " orient=" << (int)f.Orientation()
                << " eo=" << (ok ? (int)eo.Orientation() : -1) << " ok=" << (ok ? 1 : 0) << "\n";
      if (ok)
      {
        BOPTools_CoupleOfShape cs;
        cs.SetShape1(eo);
        cs.SetShape2(f);
        aLCS.Append(cs);
      }
    }
    // call GetFaceOff from each adjacent face
    occ::handle<IntTools_Context> aCtx = new IntTools_Context;
    int k = 0;
    for (NCollection_List<TopoDS_Shape>::Iterator it(aLF); it.More(); it.Next(), ++k)
    {
      const TopoDS_Face f1 = TopoDS::Face(it.Value());
      TopoDS_Edge e1;
      if (!BOPTools_AlgoTools::GetEdgeOnFace(aSec, f1, e1)) continue;
      NCollection_List<BOPTools_CoupleOfShape> aLC;
      for (NCollection_List<BOPTools_CoupleOfShape>::Iterator jt(aLCS); jt.More(); jt.Next())
      {
        if (jt.Value().Shape2().IsSame(f1)) continue;
        BOPTools_CoupleOfShape cs;
        cs.SetShape1(jt.Value().Shape1());
        cs.SetShape2(jt.Value().Shape2());
        aLC.Append(cs);
      }
      TopoDS_Face aOff;
      const bool done = BOPTools_AlgoTools::GetFaceOff(e1, f1, aLC, aOff, aCtx);
      BRepAdaptor_Surface b1(f1);
      if (!aOff.IsNull())
      {
        BRepAdaptor_Surface b2(aOff);
        std::cout << "OFF from f" << k << "(type=" << (int)b1.GetType() << ") -> type="
                  << (int)b2.GetType() << " orient=" << (int)aOff.Orientation()
                  << " done=" << (done ? 1 : 0) << "\n";
      }
      else
      {
        std::cout << "OFF from f" << k << "(type=" << (int)b1.GetType() << ") -> null done="
                  << (done ? 1 : 0) << "\n";
      }
    }
    return 0;
  }
  if (argc > 2 && std::string(argv[2]) == "--cylall")
  {
    TopoDS_Shape aPrim = BRepPrimAPI_MakeCylinder(0.4, 2.0).Shape();
    int ie = 0;
    for (TopExp_Explorer ex(aPrim, TopAbs_FACE); ex.More(); ex.Next(), ++ie)
    {
      const TopoDS_Face aF = TopoDS::Face(ex.Current());
      BRepAdaptor_Surface aBAS(aF);
      std::cout << "FACE " << ie << " type=" << (int)aBAS.GetType()
                << " orient=" << (int)aF.Orientation() << "\n";
      for (TopExp_Explorer exe(aF, TopAbs_EDGE); exe.More(); exe.Next())
      {
        const TopoDS_Edge aE = TopoDS::Edge(exe.Current());
        double f = 0., l = 0.;
        occ::handle<Geom2d_Curve> aC = BRep_Tool::CurveOnSurface(aE, aF, f, l);
        gp_Pnt2d pm(99., 99.);
        if (!aC.IsNull()) aC->D0(0.5 * (f + l), pm);
        std::cout << "  E ori=" << (int)aE.Orientation() << " uv_mid=(" << pm.X() << "," << pm.Y()
                  << ")\n";
      }
    }
    return 0;
  }
  if (argc > 2 && std::string(argv[2]) == "--cylface")
  {
    TopoDS_Shape aPrim = BRepPrimAPI_MakeCylinder(0.4, 2.0).Shape();
    for (TopExp_Explorer ex(aPrim, TopAbs_FACE); ex.More(); ex.Next())
    {
      const TopoDS_Face aF = TopoDS::Face(ex.Current());
      BRepAdaptor_Surface aBAS(aF);
      if (aBAS.GetType() != GeomAbs_Cylinder) continue;
      std::cout << "CYLFACE orient=" << (int)aF.Orientation() << "\n";
      int ie = 0;
      for (TopExp_Explorer exe(aF, TopAbs_EDGE); exe.More(); exe.Next(), ++ie)
      {
        const TopoDS_Edge aE = TopoDS::Edge(exe.Current());
        double f = 0., l = 0.;
        occ::handle<Geom2d_Curve> aC = BRep_Tool::CurveOnSurface(aE, aF, f, l);
        std::cout << "  E" << ie << " ori=" << (int)aE.Orientation()
                  << " closed=" << (BRep_Tool::IsClosed(aE, aF) ? 1 : 0)
                  << " range=[" << f << "," << l << "] null=" << (aC.IsNull() ? 1 : 0);
        if (!aC.IsNull())
        {
          gp_Pnt2d pm;
          aC->D0(0.5 * (f + l), pm);
          std::cout << " mid=(" << pm.X() << "," << pm.Y() << ")";
        }
        std::cout << "\n";
        TopoDS_Iterator it(aE);
        int iv = 0;
        for (; it.More(); it.Next(), ++iv)
        {
          const TopoDS_Vertex aV = TopoDS::Vertex(it.Value());
          const double        t  = BRep_Tool::Parameter(aV, aE, aF);
          gp_Pnt2d            pv(99., 99.);
          if (!aC.IsNull()) aC->D0(t, pv);
          std::cout << "    V" << iv << " ori=" << (int)aV.Orientation() << " param=" << t
                    << " uv=(" << pv.X() << "," << pv.Y() << ")\n";
        }
        // Static replay of BOPAlgo_WireSplitter::GetNextVertex for each child
        // forced FORWARD (the Path call passes pVa.Orientation(FORWARD)).
        for (int k = 0; k < 2; ++k)
        {
          TopoDS_Iterator itk(aE);
          int             m = 0;
          for (; itk.More() && m < k; itk.Next(), ++m) {}
          if (!itk.More()) continue;
          TopoDS_Vertex aV = TopoDS::Vertex(itk.Value());
          aV.Orientation(TopAbs_FORWARD);
          TopoDS_Iterator itv(aE);
          TopoDS_Vertex   aV1;
          for (; itv.More(); itv.Next())
          {
            const TopoDS_Shape& aVx = itv.Value();
            if (!aVx.IsEqual(aV))
            {
              aV1 = TopoDS::Vertex(aVx);
              break;
            }
          }
          double   tn = BRep_Tool::Parameter(aV1, aE, aF);
          gp_Pnt2d pn(99., 99.);
          if (!aC.IsNull()) aC->D0(tn, pn);
          std::cout << "    NEXT(from V" << k << " forced F) -> ori=" << (int)aV1.Orientation()
                    << " param=" << tn << " uv=(" << pn.X() << "," << pn.Y() << ")\n";
        }
      }
    }
    return 0;
  }
  if (argc > 2 && (std::string(argv[2]) == "--cyl" || std::string(argv[2]) == "--sph"))
  {
    TopoDS_Shape aPrim;
    if (std::string(argv[2]) == "--cyl")
      aPrim = BRepPrimAPI_MakeCylinder(0.4, 2.0).Shape();
    else
      aPrim = BRepPrimAPI_MakeSphere(1.0).Shape();
    GProp_GProps aP;
    BRepGProp::VolumeProperties(aPrim, aP);
    std::cout << "PRIM volume=" << aP.Mass() << " faces=";
    int nf = 0;
    for (TopExp_Explorer ex(aPrim, TopAbs_FACE); ex.More(); ex.Next()) ++nf;
    std::cout << nf << "\n";
    int i = 0;
    for (TopExp_Explorer ex(aPrim, TopAbs_FACE); ex.More(); ex.Next(), ++i)
    {
      const TopoDS_Face aF = TopoDS::Face(ex.Current());
      BRepAdaptor_Surface aBAS(aF);
      BRepGProp_Face aBF; aBF.Load(aF);
      double u1,u2,v1,v2; aBF.Bounds(u1,u2,v1,v2);
      std::cout << "PRIMFACE " << i << " type=" << (int)aBAS.GetType()
                << " orient=" << (int)aF.Orientation() << " uv=[" << u1 << "," << u2 << "]x["
                << v1 << "," << v2 << "]\n";
      int j = 0;
      for (TopExp_Explorer exe(aF, TopAbs_EDGE); exe.More(); exe.Next(), ++j)
      {
        const TopoDS_Edge aE = TopoDS::Edge(exe.Current());
        if (!aBF.Load(aE)) { std::cout << "  E" << j << " load=false\n"; continue; }
        const double lp1 = aBF.FirstParameter(), lp2 = aBF.LastParameter();
        const double lm = 0.5 * (lp1 + lp2);
        gp_Pnt2d Puv; gp_Vec2d Vuv;
        aBF.D12d(lm, Puv, Vuv);
        std::cout << "  E" << j << " eo=" << (int)aE.Orientation() << " range=[" << lp1 << ","
                  << lp2 << "] midU=" << Puv.X() << " midV=" << Puv.Y() << " du=" << Vuv.X()
                  << " dv=" << Vuv.Y() << "\n";
      }
    }
    return 0;
  }
  if (argc > 3 && std::string(argv[2]) == "--facek")
  {
    const int k = std::atoi(argv[3]);
    int i = 0;
    for (TopExp_Explorer ex(aShape, TopAbs_FACE); ex.More(); ex.Next(), ++i)
    {
      if (i != k) continue;
      const TopoDS_Face aF = TopoDS::Face(ex.Current());
      BRepAdaptor_Surface aBAS(aF);
      BRepGProp_Face aBF; aBF.Load(aF);
      double u1,u2,v1,v2; aBF.Bounds(u1,u2,v1,v2);
      std::cout << "FK face=" << k << " type=" << (int)aBAS.GetType() << " orient=" << (int)aF.Orientation()
                << " uv=[" << u1 << "," << u2 << "]x[" << v1 << "," << v2 << "]\n";
      int j = 0;
      for (TopExp_Explorer exe(aF, TopAbs_EDGE); exe.More(); exe.Next(), ++j)
      {
        const TopoDS_Edge aE = TopoDS::Edge(exe.Current());
        if (!aBF.Load(aE)) { std::cout << "  FKE" << j << " load=false\n"; continue; }
        const double lp1 = aBF.FirstParameter(), lp2 = aBF.LastParameter();
        const double lm = 0.5 * (lp1 + lp2);
        gp_Pnt2d Puv; gp_Vec2d Vuv; aBF.D12d(lm, Puv, Vuv);
        std::cout << "  FKE" << j << " eo=" << (int)aE.Orientation() << " range=[" << lp1 << "," << lp2
                  << "] midU=" << Puv.X() << " midV=" << Puv.Y() << " du=" << Vuv.X() << " dv=" << Vuv.Y() << "\n";
      }
      return 0;
    }
    return 0;
  }
  int aNbFaces = 0;
  for (TopExp_Explorer ex(aShape, TopAbs_FACE); ex.More(); ex.Next())
    ++aNbFaces;
  std::cout << "SHAPE null=0 faces=" << aNbFaces << " rootOrient=" << (int)aShape.Orientation()
            << " rootType=" << (int)aShape.ShapeType() << "\n";
  for (TopExp_Explorer ex(aShape, TopAbs_SOLID); ex.More(); ex.Next())
    std::cout << "SOLID orient=" << (int)ex.Current().Orientation() << "\n";
  for (TopExp_Explorer ex(aShape, TopAbs_SHELL); ex.More(); ex.Next())
    std::cout << "SHELL orient=" << (int)ex.Current().Orientation() << "\n";
  {
    int j = 0;
    for (TopExp_Explorer exs(aShape, TopAbs_SHELL); exs.More(); exs.Next())
    {
      std::cout << "SHELLCHILD shellOrient=" << (int)exs.Current().Orientation() << "\n";
      // cumOri=false: the *raw* stored orientation of each face, without
      // composing the shell's own orientation.
      for (TopoDS_Iterator it(exs.Current(), Standard_False, Standard_False); it.More();
           it.Next(), ++j)
      {
        std::cout << "RAW ptr=" << (const void*)it.Value().TShape().get()
                  << " orient=" << (int)it.Value().Orientation() << "\n";
      }
    }
    for (TopExp_Explorer exf(aShape, TopAbs_FACE); exf.More(); exf.Next())
    {
      std::cout << "EXPLORER ptr=" << (const void*)exf.Current().TShape().get()
                << " orient=" << (int)exf.Current().Orientation() << "\n";
    }
  }

  BRepCheck_Analyzer aCheck(aShape);
  std::cout << "BRepCheck valid=" << (aCheck.IsValid() ? 1 : 0) << "\n";

  GProp_GProps aProps;
  BRepGProp::VolumeProperties(aShape, aProps);
  std::cout.setf(std::ios::fixed);
  std::cout.precision(6);
  std::cout << "BRepGProp volume=" << aProps.Mass() << "\n";
  const gp_Pnt aCenter = aProps.CentreOfMass();
  std::cout << "BRepGProp centre=(" << aCenter.X() << "," << aCenter.Y() << "," << aCenter.Z()
            << ")\n";

  int i = 0;
  for (TopExp_Explorer ex(aShape, TopAbs_FACE); ex.More(); ex.Next(), ++i)
  {
    const TopoDS_Face aFace = TopoDS::Face(ex.Current());
    BRepAdaptor_Surface aBAS(aFace);
    const GeomAbs_SurfaceType aType = aBAS.GetType();

    BRepGProp_Face aBF;
    aBF.Load(aFace);
    double u1 = 0., u2 = 0., v1 = 0., v2 = 0.;
    aBF.Bounds(u1, u2, v1, v2);
    gp_Pnt aP;
    gp_Vec aN;
    aBF.Normal(0.5 * (u1 + u2), 0.5 * (v1 + v2), aP, aN);
    const double anOutward = gp_Vec(aCenter, aP).Dot(aN);
    // Raw (unflipped) natural normal: D1U ^ D1V at the same parameters.
    gp_Pnt aP2;
    gp_Vec aDU, aDV;
    aBAS.D1(0.5 * (u1 + u2), 0.5 * (v1 + v2), aP2, aDU, aDV);
    const gp_Vec aRaw = aDU.Crossed(aDV);
    const double anOutwardRaw = gp_Vec(aCenter, aP2).Dot(aRaw);

    const char* aOri = "?";
    switch (aFace.Orientation())
    {
      case TopAbs_FORWARD: aOri = "FWD"; break;
      case TopAbs_REVERSED: aOri = "REV"; break;
      case TopAbs_INTERNAL: aOri = "INT"; break;
      case TopAbs_EXTERNAL: aOri = "EXT"; break;
    }
    std::cout << "face " << i << " type=" << (int)aType << " orient=" << aOri
              << " outw=" << anOutward << " raw=" << anOutwardRaw << " uv=[" << u1 << "," << u2 << "]x[" << v1 << "," << v2
              << "]";
    if (frames)
    {
      if (aType == GeomAbs_Cylinder)
      {
        const gp_Cylinder aCyl = aBAS.Cylinder();
        const gp_Ax3 aAx = aCyl.Position();
        std::cout << " cyl loc=(" << aAx.Location().X() << "," << aAx.Location().Y() << ","
                  << aAx.Location().Z() << ") Z=(" << aAx.Direction().X() << ","
                  << aAx.Direction().Y() << "," << aAx.Direction().Z() << ") X=("
                  << aAx.XDirection().X() << "," << aAx.XDirection().Y() << ","
                  << aAx.XDirection().Z() << ") Y=(" << aAx.YDirection().X() << ","
                  << aAx.YDirection().Y() << "," << aAx.YDirection().Z()
                  << ") R=" << aCyl.Radius();
      }
      else if (aType == GeomAbs_Plane)
      {
        const gp_Pln aPln = aBAS.Plane();
        const gp_Ax3 aAx = aPln.Position();
        std::cout << " pln loc=(" << aAx.Location().X() << "," << aAx.Location().Y() << ","
                  << aAx.Location().Z() << ") Z=(" << aAx.Direction().X() << ","
                  << aAx.Direction().Y() << "," << aAx.Direction().Z() << ") X=("
                  << aAx.XDirection().X() << "," << aAx.XDirection().Y() << ","
                  << aAx.XDirection().Z() << ")";
      }
      else if (aType == GeomAbs_Sphere)
      {
        const gp_Sphere aSph = aBAS.Sphere();
        const gp_Ax3   aAx  = aSph.Position();
        std::cout << " sph loc=(" << aAx.Location().X() << "," << aAx.Location().Y() << ","
                  << aAx.Location().Z() << ") Z=(" << aAx.Direction().X() << ","
                  << aAx.Direction().Y() << "," << aAx.Direction().Z() << ") R=" << aSph.Radius();
      }
    }
    std::cout << "\n";
    // Per-face wire edges with their composed orientation and endpoints, so the
    // port's edge-orientation mapping can be diffed edge by edge.
    {
      int ei = 0;
      for (TopExp_Explorer exw(aFace, TopAbs_WIRE); exw.More(); exw.Next())
      {
        for (TopExp_Explorer exe(exw.Current(), TopAbs_EDGE); exe.More(); exe.Next(), ++ei)
        {
          const TopoDS_Edge anE = TopoDS::Edge(exe.Current());
          double f = 0., l = 0.;
          const occ::handle<Geom_Curve>& aC = BRep_Tool::Curve(anE, f, l);
          TopoDS_Vertex aV1, aV2;
          TopExp::Vertices(anE, aV1, aV2);
          const gp_Pnt p1 = BRep_Tool::Pnt(aV1);
          const gp_Pnt p2 = BRep_Tool::Pnt(aV2);
          std::cout << "EDGE face=" << i << " n=" << ei
                    << " orient=" << (int)anE.Orientation() << " null3d=" << (aC.IsNull() ? 1 : 0)
                    << " p1=(" << p1.X() << "," << p1.Y() << "," << p1.Z() << ")"
                    << " p2=(" << p2.X() << "," << p2.Y() << "," << p2.Z() << ")\n";
        }
      }
    }
  }
  return 0;
}
