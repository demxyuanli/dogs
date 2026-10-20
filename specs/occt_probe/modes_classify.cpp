// `argv[2]` modes: `--fcls2d`, `--cylall`, `--cylface`, `--cyl` / `--sph`,
// `--facek`.
#include "occt_probe_common.h"

std::optional<int> runClassifyModes(const ProbeCtx& theCtx)
{
  const int           argc = theCtx.argc;
  char** const        argv = theCtx.argv;
  const TopoDS_Shape& aShape = theCtx.shape;
  if (argc > 2 && std::string(argv[2]) == "--fcls2d")
  {
    // Ground truth for `BRepTopAdaptor_FClass2d` on the cylinder lateral face of
    // `BRepPrimAPI_MakeCylinder(0.4, 2.0)`: the `BRepTools_WireExplorer` order
    // (which pcurve each seam occurrence walks), the classifier state over a u
    // grid at v=1, and the `IntTools_Context::IsPointInOnFace` path for the
    // pave point (0,-0.4,1).
    TopoDS_Shape aPrim = BRepPrimAPI_MakeCylinder(0.4, 2.0).Shape();
    for (TopExp_Explorer ex(aPrim, TopAbs_FACE); ex.More(); ex.Next())
    {
      TopoDS_Face aF = TopoDS::Face(ex.Current());
      BRepAdaptor_Surface aBAS(aF, false);
      if (aBAS.GetType() != GeomAbs_Cylinder) continue;
      aF.Orientation(TopAbs_FORWARD);
      std::cout << "FCLS2D cylinder face\n";
      for (TopExp_Explorer exw(aF, TopAbs_WIRE); exw.More(); exw.Next())
      {
        const TopoDS_Wire       aW = TopoDS::Wire(exw.Current());
        BRepTools_WireExplorer  aWE;
        int                     k = 0;
        for (aWE.Init(aW, aF); aWE.More(); aWE.Next(), ++k)
        {
          const TopoDS_Edge         aE = aWE.Current();
          double                    f = 0., l = 0.;
          occ::handle<Geom2d_Curve> aC = BRep_Tool::CurveOnSurface(aE, aF, f, l);
          gp_Pnt2d                  p1(99., 99.), p2(99., 99.), pv(99., 99.);
          if (!aC.IsNull())
          {
            aC->D0(f, p1);
            aC->D0(l, p2);
          }
          const TopoDS_Vertex aV = aWE.CurrentVertex();
          if (!aV.IsNull() && !aC.IsNull())
          {
            aC->D0(BRep_Tool::Parameter(aV, aE, aF), pv);
          }
          std::cout << "  WE" << k << " ori=" << (int)aWE.Orientation()
                    << " closed=" << (BRep_Tool::IsClosed(aE, aF) ? 1 : 0) << " range=[" << f << ","
                    << l << "] pc_f=(" << p1.X() << "," << p1.Y() << ") pc_l=(" << p2.X() << ","
                    << p2.Y() << ") curV=(" << pv.X() << "," << pv.Y() << ")\n";
        }
      }
      BRepTopAdaptor_FClass2d aFc(aF, Precision::Confusion());
      BRepClass_FaceClassifier aCls;
      for (int i = 0; i <= 16; ++i)
      {
        const gp_Pnt2d aP(2.0 * M_PI * i / 16.0, 1.0);
        aCls.Perform(aF, aP, Precision::Confusion());
        std::cout << "  P u=" << aP.X() << " fclass=" << (int)aFc.Perform(aP)
                  << " faceclass=" << (int)aCls.State() << "\n";
      }
      GeomAPI_ProjectPointOnSurf aProj(gp_Pnt(0.0, -0.4, 1.0), BRep_Tool::Surface(aF));
      if (aProj.IsDone())
      {
        double u = 0., v = 0.;
        aProj.LowerDistanceParameters(u, v);
        const gp_Pnt2d aP2d(u, v);
        aCls.Perform(aF, aP2d, Precision::Confusion());
        std::cout << "  ISPOINT uv=(" << u << "," << v << ") fclass=" << (int)aFc.Perform(aP2d)
                  << " faceclass=" << (int)aCls.State() << "\n";
      }
      std::cout << "  InfPoint fclass=" << (int)aFc.PerformInfinitePoint() << "\n";
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
  return std::nullopt;
}
