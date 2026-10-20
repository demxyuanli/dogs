// Boolean ground truth: `--ds`, `--fuse`, `--gf`, `--cut` / `--common`,
// `--dir` / `--dir2`, `--faceoff` / `--faceoff2`.
#include "occt_probe_common.h"

std::optional<int> runBooleanModes(const ProbeCtx& theCtx)
{
  const int           argc = theCtx.argc;
  char** const        argv = theCtx.argv;
  const TopoDS_Shape& aShape = theCtx.shape;
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
  if (argc > 2 && std::string(argv[2]) == "--gf")
  {
    // General Fuse of the same operands. The GF result carries the split
    // images of both arguments, so the solid images (the BOPAlgo_SplitSolid
    // areas of BuildSplitSolids) can be read off directly.
    TopoDS_Shape aBox = BRepPrimAPI_MakeBox(gp_Pnt(-1., -1., -1.), gp_Pnt(1., 1., 1.)).Shape();
    TopoDS_Shape aCyl = BRepPrimAPI_MakeCylinder(0.4, 2.0).Shape();
    NCollection_List<TopoDS_Shape> aLS;
    aLS.Append(aBox);
    aLS.Append(aCyl);
    BOPAlgo_Builder aGF;
    aGF.SetArguments(aLS);
    aGF.Perform();
    if (aGF.HasErrors()) { std::cout << "GF failed\n"; return 1; }
    const TopoDS_Shape aRes = aGF.Shape();
    int i = 0;
    for (TopExp_Explorer ex(aRes, TopAbs_SOLID); ex.More(); ex.Next(), ++i)
    {
      const TopoDS_Solid aSol = TopoDS::Solid(ex.Current());
      GProp_GProps aP; BRepGProp::VolumeProperties(aSol, aP);
      int nf = 0; for (TopExp_Explorer exf(aSol, TopAbs_FACE); exf.More(); exf.Next()) ++nf;
      int nc = 0;
      for (TopoDS_Iterator it(aSol); it.More(); it.Next())
      {
        if (it.Value().ShapeType() == TopAbs_SHELL) ++nc;
      }
      std::cout << "GFSOL " << i << " orient=" << (int)aSol.Orientation()
                << " volume=" << aP.Mass() << " faces=" << nf << " subshapes=" << nc << "\n";
      int j = 0;
      for (TopExp_Explorer exf(aSol, TopAbs_FACE); exf.More(); exf.Next(), ++j)
      {
        const TopoDS_Face aF = TopoDS::Face(exf.Current());
        BRepAdaptor_Surface aBAS(aF);
        GProp_GProps aFA; BRepGProp::SurfaceProperties(aF, aFA);
        std::cout << "  GFSF " << j << " type=" << (int)aBAS.GetType()
                  << " orient=" << (int)aF.Orientation() << " area=" << aFA.Mass() << "\n";
      }
    }
    std::cout << "GF faces=" ;
    int nall = 0; for (TopExp_Explorer ex(aRes, TopAbs_FACE); ex.More(); ex.Next()) ++nall;
    std::cout << nall << "\n";
    return 0;
  }
  if (argc > 2 && (std::string(argv[2]) == "--cut" || std::string(argv[2]) == "--common"))
  {
    // Same operands as --fuse: box[-1,1]^3 minus / intersected with
    // BRepPrimAPI_MakeCylinder(0.4, 2.0).
    const bool bCut = (std::string(argv[2]) == "--cut");
    TopoDS_Shape aBox = BRepPrimAPI_MakeBox(gp_Pnt(-1., -1., -1.), gp_Pnt(1., 1., 1.)).Shape();
    TopoDS_Shape aCyl = BRepPrimAPI_MakeCylinder(0.4, 2.0).Shape();
    TopoDS_Shape aRes;
    if (bCut)
    {
      BRepAlgoAPI_Cut anOp(aBox, aCyl);
      anOp.Build();
      if (!anOp.IsDone()) { std::cout << "CUT failed\n"; return 1; }
      aRes = anOp.Shape();
    }
    else
    {
      BRepAlgoAPI_Common anOp(aBox, aCyl);
      anOp.Build();
      if (!anOp.IsDone()) { std::cout << "COMMON failed\n"; return 1; }
      aRes = anOp.Shape();
    }
    GProp_GProps aP; BRepGProp::VolumeProperties(aRes, aP);
    std::cout << (bCut ? "CUT" : "COMMON") << " volume=" << aP.Mass() << " faces=";
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
      std::cout << "AFFACE " << i << " type=" << (int)aBAS.GetType() << " orient=" << (int)aF.Orientation()
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
  return std::nullopt;
}
