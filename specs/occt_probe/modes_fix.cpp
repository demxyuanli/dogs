// `argv[2]` modes: `--nbdump`, `--obdump`, `--vmap`, `--degchk` / `--degchkraw` /
// `--degchksph`, `--facefix`, `--wirefix` / `--wirefixraw`.
#include "occt_probe_common.h"

std::optional<int> runFixFaceModes(const ProbeCtx& theCtx)
{
  const int           argc = theCtx.argc;
  char** const        argv = theCtx.argv;
  const TopoDS_Shape& aShape = theCtx.shape;
  // TEMP T-121: `--nbdump` prints, per face, the values `ShapeFix_Face::
  // isNeedAddNaturalBound` (`ShapeFix_Face.cxx:1120-1160`) reads: the surface
  // U/V bounds and periodicity, the wire count, `ShapeAnalysis::IsOuterBound`
  // (`ShapeAnalysis.cxx:203-230`) and, per wire, `ShapeAnalysis::TotCross2D`
  // (`ShapeAnalysis.cxx:115-150`) with every edge orientation and p-curve
  // range, so the port's decision can be compared step by step.
  if (argc > 2 && std::string(argv[2]) == "--nbdump")
  {
    int k = 0;
    for (TopExp_Explorer fx(aShape, TopAbs_FACE); fx.More(); fx.Next(), ++k)
    {
      const TopoDS_Face aF = TopoDS::Face(fx.Current());
      BRepAdaptor_Surface aS(aF);
      int aNW = 0;
      for (TopExp_Explorer wx(aF, TopAbs_WIRE); wx.More(); wx.Next())
      {
        ++aNW;
      }
      std::cout << "NB face " << k << " type=" << (int)aS.GetType()
                << " uper=" << (aS.IsUPeriodic() ? 1 : 0) << " vper=" << (aS.IsVPeriodic() ? 1 : 0)
                << " u=[" << aS.FirstUParameter() << "," << aS.LastUParameter() << "]"
                << " v=[" << aS.FirstVParameter() << "," << aS.LastVParameter() << "]"
                << " wires=" << aNW << " outer=" << (ShapeAnalysis::IsOuterBound(aF) ? 1 : 0) << "\n";
      int aWI = 0;
      for (TopExp_Explorer wx(aF, TopAbs_WIRE); wx.More(); wx.Next(), ++aWI)
      {
        TopoDS_Wire                       aW = TopoDS::Wire(wx.Current());
        occ::handle<ShapeExtend_WireData> aWD = new ShapeExtend_WireData(aW);
        std::cout << "NB  w" << aWI << " ori=" << (int)aW.Orientation()
                  << " n=" << aWD->NbEdges() << " cross=" << ShapeAnalysis::TotCross2D(aWD, aF)
                  << "\n";
        for (int i = 1; i <= aWD->NbEdges(); ++i)
        {
          const TopoDS_Edge&        e = aWD->Edge(i);
          double                    f = 0., l = 0.;
          occ::handle<Geom2d_Curve> c2 = BRep_Tool::CurveOnSurface(e, aF, f, l);
          gp_Pnt2d                  b(0., 0.), en(0., 0.);
          if (!c2.IsNull())
          {
            b  = c2->Value(f);
            en = c2->Value(l);
          }
          std::cout << "NB   e" << i << " ori=" << (int)e.Orientation()
                    << " closed=" << (BRep_Tool::IsClosed(e, aF) ? 1 : 0)
                    << " deg=" << (BRep_Tool::Degenerated(e) ? 1 : 0) << " uv=(" << b.X() << ","
                    << b.Y() << ")->(" << en.X() << "," << en.Y() << ")\n";
        }
      }
    }
    return 0;
  }
  // TEMP T-121: `--obdump <face>` replays `ShapeAnalysis::IsOuterBound`
  // (`ShapeAnalysis.cxx:203-230`) frame for frame: `F.Orientation(FORWARD)`,
  // the wire handed to `ShapeExtend_WireData`, its edge order/orientations and
  // the `TotCross2D` walk, so the port's own wire data can be compared against
  // the exact sequence OCCT sums.
  if (argc > 3 && std::string(argv[2]) == "--obdump")
  {
    int  kTarget = std::atoi(argv[3]);
    int  k = 0;
    TopoDS_Face aFound;
    for (TopExp_Explorer fx(aShape, TopAbs_FACE); fx.More(); fx.Next(), ++k)
    {
      if (k == kTarget)
      {
        aFound = TopoDS::Face(fx.Current());
        break;
      }
    }
    if (aFound.IsNull())
    {
      std::cout << "OB face " << kTarget << " not found\n";
      return 0;
    }
    TopoDS_Face F = aFound;
    std::cout << "OB face " << kTarget << " ori=" << (int)aFound.Orientation()
              << " outer=" << (ShapeAnalysis::IsOuterBound(aFound) ? 1 : 0) << "\n";
    F.Orientation(TopAbs_FORWARD);
    TopoDS_Wire aW;
    for (TopExp_Explorer wx(F, TopAbs_WIRE); wx.More(); wx.Next())
    {
      aW = TopoDS::Wire(wx.Current());
    }
    std::cout << "OB wire ori=" << (int)aW.Orientation() << "\n";
    occ::handle<ShapeExtend_WireData> aWD = new ShapeExtend_WireData(aW);
    std::cout << "OB n=" << aWD->NbEdges() << " cross=" << ShapeAnalysis::TotCross2D(aWD, F)
              << "\n";
    for (int i = 1; i <= aWD->NbEdges(); ++i)
    {
      const TopoDS_Edge&        e = aWD->Edge(i);
      double                    f = 0., l = 0.;
      occ::handle<Geom2d_Curve> c2 = BRep_Tool::CurveOnSurface(e, F, f, l);
      gp_Pnt2d                  b(0., 0.), en(0., 0.);
      if (!c2.IsNull())
      {
        b  = c2->Value(f);
        en = c2->Value(l);
      }
      std::cout << "OB   e" << i << " ori=" << (int)e.Orientation()
                << " closed=" << (BRep_Tool::IsClosed(e, F) ? 1 : 0) << " t=(" << f << "," << l
                << ") uv=(" << b.X() << "," << b.Y() << ")->(" << en.X() << "," << en.Y()
                << ")\n";
      // PCurve1 / PCurve2 as the reader baked them (`BRep_CurveOnClosedSurface`):
      // a FORWARD copy yields PCurve1, a REVERSED copy yields PCurve2.
      if (BRep_Tool::IsClosed(e, F))
      {
        const char* aTag[2] = {"PCurve1", "PCurve2"};
        for (int s = 0; s < 2; ++s)
        {
          TopoDS_Edge ec = e;
          ec.Orientation(s == 0 ? TopAbs_FORWARD : TopAbs_REVERSED);
          double                    cf = 0., cl = 0.;
          occ::handle<Geom2d_Curve> cc = BRep_Tool::CurveOnSurface(ec, F, cf, cl);
          if (cc.IsNull())
            continue;
          const gp_Pnt2d p0 = cc->Value(cf);
          const gp_Pnt2d p1 = cc->Value(cl);
          std::cout << "OB    " << aTag[s] << " t=(" << cf << "," << cl << ") uv=(" << p0.X()
                    << "," << p0.Y() << ")->(" << p1.X() << "," << p1.Y() << ")\n";
        }
      }
    }
    return 0;
  }
  // TEMP T-121: `--vmap` walks every wire of every face and prints, per edge,
  // the identity and point of `ShapeAnalysis_Edge::FirstVertex` / `LastVertex`
  // (`ShapeAnalysis_Edge.cxx:228-257`), so edge closure (`v1.IsSame(v2)`) and
  // vertex sharing between adjacent edges can be read off. The identity is the
  // raw `TShape` address, stable inside one process.
  if (argc > 2 && std::string(argv[2]) == "--vmap")
  {
    int k = 0;
    for (TopExp_Explorer fx(aShape, TopAbs_FACE); fx.More(); fx.Next(), ++k)
    {
      std::cout << "VMAP face " << k << "\n";
      int aWI = 0;
      for (TopExp_Explorer wx(fx.Current(), TopAbs_WIRE); wx.More(); wx.Next(), ++aWI)
      {
        ShapeAnalysis_Edge                aEA;
        ShapeExtend_WireData              aWD(TopoDS::Wire(wx.Current()));
        std::cout << "VMAP  w" << aWI << " n=" << aWD.NbEdges() << "\n";
        for (int i = 1; i <= aWD.NbEdges(); ++i)
        {
          const TopoDS_Edge   e  = aWD.Edge(i);
          const TopoDS_Vertex v1 = aEA.FirstVertex(e);
          const TopoDS_Vertex v2 = aEA.LastVertex(e);
          std::cout << "VMAP   e" << i << " ori=" << (int)e.Orientation()
                    << " v1=" << (v1.IsNull() ? 0L : (long)v1.TShape().get())
                    << " v2=" << (v2.IsNull() ? 0L : (long)v2.TShape().get())
                    << " same=" << ((!v1.IsNull() && !v2.IsNull() && v1.IsSame(v2)) ? 1 : 0);
          if (!v1.IsNull())
          {
            const gp_Pnt p = BRep_Tool::Pnt(v1);
            std::cout << " p1=(" << p.X() << "," << p.Y() << "," << p.Z() << ")";
          }
          if (!v2.IsNull())
          {
            const gp_Pnt p = BRep_Tool::Pnt(v2);
            std::cout << " p2=(" << p.X() << "," << p.Y() << "," << p.Z() << ")";
          }
          std::cout << "\n";
        }
      }
    }
    return 0;
  }
  // TEMP: `--degchk <faceIdx>` replays the reader-time wire passes up to
  // `ShapeFix_Wire::FixEdgeCurves` and then prints
  // `ShapeAnalysis_Wire::CheckDegenerated` per edge (endpoints, singularity
  // list, DONE/FAIL bits) plus the edge count after each `FixDegenerated(i)`.
  int  aDegChkFace = -1;
  bool aDegChkSph  = false;
  for (int i = 2; i + 1 < argc; ++i)
  {
    if (std::string(argv[i]) == "--degchk" || std::string(argv[i]) == "--degchkraw")
    {
      aDegChkFace = std::atoi(argv[i + 1]);
      break;
    }
  }
  for (int i = 2; i < argc; ++i)
  {
    if (std::string(argv[i]) == "--degchksph")
    {
      aDegChkSph = true;
    }
  }
  if (aDegChkFace >= 0 || aDegChkSph)
  {
    TopoDS_Face aFace;
    int         aFaceIdx = aDegChkFace;
    if (aDegChkSph)
    {
      int k = 0;
      for (TopExp_Explorer fx(aShape, TopAbs_FACE); fx.More(); fx.Next(), ++k)
      {
        const TopoDS_Face        fc = TopoDS::Face(fx.Current());
        const occ::handle<Geom_Surface> sf = BRep_Tool::Surface(fc);
        if (!sf.IsNull() && (new ShapeAnalysis_Surface(sf))->NbSingularities(Precision::Confusion()) > 0)
        {
          const TopoDS_Wire w0 = TopoDS::Wire(TopExp_Explorer(fc, TopAbs_WIRE).Current());
          ShapeFix_Wire     s0;
          s0.SetFace(fc, new ShapeAnalysis_Surface(sf));
          s0.Load(w0);
          s0.FixReorder();
          s0.FixSmall(0);
          s0.FixConnected();
          s0.FixEdgeCurves();
          std::cout << "SPH face=" << k << " n=" << s0.NbEdges();
          occ::handle<ShapeExtend_WireData> wd = new ShapeExtend_WireData(s0.Wire());
          ShapeAnalysis_Edge                se;
          for (int j = 1; j <= wd->NbEdges(); ++j)
          {
            const TopoDS_Edge&        e = wd->Edge(j);
            double                    f = 0., l = 0.;
            occ::handle<Geom2d_Curve> c2 = BRep_Tool::CurveOnSurface(e, fc, f, l);
            std::cout << " [" << (int)e.Orientation() << (BRep_Tool::Degenerated(e) ? "d" : "");
            if (!c2.IsNull())
            {
              const gp_Pnt2d b = c2->Value(f), en = c2->Value(l);
              std::cout << " (" << b.X() << "," << b.Y() << ")->(" << en.X() << "," << en.Y()
                        << ")";
            }
            std::cout << "]";
          }
          std::cout << "\n";
        }
      }
      return 0;
    }
    {
      int k = 0;
      for (TopExp_Explorer fx(aShape, TopAbs_FACE); fx.More(); fx.Next(), ++k)
      {
        if (k == aFaceIdx)
        {
          aFace = TopoDS::Face(fx.Current());
          break;
        }
      }
    }
    if (aFace.IsNull())
    {
      std::cout << "DC no-face\n";
      return 0;
    }
    occ::handle<Geom_Surface>          aSurf = BRep_Tool::Surface(aFace);
    occ::handle<ShapeAnalysis_Surface> aSAS  = new ShapeAnalysis_Surface(aSurf);
    const TopoDS_Wire aW = TopoDS::Wire(TopExp_Explorer(aFace, TopAbs_WIRE).Current());
    ShapeFix_Wire aSFW;
    aSFW.SetFace(aFace, aSAS);
    aSFW.Load(aW);
    aSFW.FixReorder();
    aSFW.FixSmall(0);
    aSFW.FixConnected();
    aSFW.FixEdgeCurves();
    occ::handle<ShapeExtend_WireData> aWD = new ShapeExtend_WireData(aSFW.Wire());
    ShapeAnalysis_Wire                aWA(aWD, aFace, Precision::Confusion());
    std::cout << "DC n=" << aWD->NbEdges() << "\n";
    for (int j = 1; j <= aSAS->NbSingularities(Precision::Confusion()); ++j)
    {
      double    preci, par1, par2;
      gp_Pnt2d  pd1, pd2;
      gp_Pnt    p3;
      bool      uiso = false;
      aSAS->Singularity(j, preci, p3, pd1, pd2, par1, par2, uiso);
      std::cout << "DC sing" << j << " preci=" << preci << " p3=(" << p3.X() << "," << p3.Y()
                << "," << p3.Z() << ") pd1=(" << pd1.X() << "," << pd1.Y() << ") pd2=("
                << pd2.X() << "," << pd2.Y() << ") uiso=" << (uiso ? 1 : 0) << "\n";
    }
    ShapeAnalysis_Edge sae;
    for (int i = 1; i <= aWD->NbEdges(); ++i)
    {
      const int   n2 = i;
      const int   n1 = (n2 > 1 ? n2 - 1 : aWD->NbEdges());
      TopoDS_Edge E1 = aWD->Edge(n1);
      TopoDS_Edge E2 = aWD->Edge(n2);
      gp_Pnt2d    pd1, pd2;
      aWA.CheckDegenerated(i, pd1, pd2);
      std::cout << "DC i=" << i << " n1=" << n1
                << " done=" << (aWA.LastCheckStatus(ShapeExtend_DONE) ? 1 : 0)
                << " done1=" << (aWA.LastCheckStatus(ShapeExtend_DONE1) ? 1 : 0)
                << " done2=" << (aWA.LastCheckStatus(ShapeExtend_DONE2) ? 1 : 0)
                << " fail1=" << (aWA.LastCheckStatus(ShapeExtend_FAIL1) ? 1 : 0)
                << " fail2=" << (aWA.LastCheckStatus(ShapeExtend_FAIL2) ? 1 : 0)
                << " pd1=(" << pd1.X() << "," << pd1.Y() << ") pd2=(" << pd2.X() << ","
                << pd2.Y() << ")";
      TopoDS_Vertex vp = sae.FirstVertex(E1), v0 = sae.LastVertex(E1);
      TopoDS_Vertex v1 = sae.FirstVertex(E2), v2 = sae.LastVertex(E2);
      if (!vp.IsNull() && !v0.IsNull() && !v1.IsNull() && !v2.IsNull())
      {
        gp_Pnt pp = BRep_Tool::Pnt(vp), p0 = BRep_Tool::Pnt(v0);
        gp_Pnt p1 = BRep_Tool::Pnt(v1), p2 = BRep_Tool::Pnt(v2);
        std::cout << " pp=(" << pp.X() << "," << pp.Y() << "," << pp.Z() << ") p0=(" << p0.X()
                  << "," << p0.Y() << "," << p0.Z() << ") p1=(" << p1.X() << "," << p1.Y()
                  << "," << p1.Z() << ") p2=(" << p2.X() << "," << p2.Y() << "," << p2.Z()
                  << ")";
      }
      std::cout << "\n";
    }
    for (int i = aSFW.NbEdges(); i >= 1; --i)
    {
      const int  aBefore = aSFW.NbEdges();
      const bool aRes    = aSFW.FixDegenerated(i);
      std::cout << "DC fix i=" << i << " r=" << (aRes ? 1 : 0) << " n=" << aSFW.NbEdges()
                << " was=" << aBefore << "\n";
    }
    occ::handle<ShapeExtend_WireData> aEnd = new ShapeExtend_WireData(aSFW.Wire());
    for (int i = 1; i <= aEnd->NbEdges(); ++i)
    {
      const TopoDS_Edge&        e  = aEnd->Edge(i);
      double                    f  = 0., l = 0.;
      occ::handle<Geom2d_Curve> c2 = BRep_Tool::CurveOnSurface(e, aFace, f, l);
      std::cout << "DC e" << i << " ori=" << (int)e.Orientation() << " deg="
                << (BRep_Tool::Degenerated(e) ? 1 : 0);
      if (!c2.IsNull())
      {
        const gp_Pnt2d b = c2->Value(f), en = c2->Value(l);
        std::cout << " uv=(" << b.X() << "," << b.Y() << ")->(" << en.X() << "," << en.Y() << ")";
      }
      std::cout << "\n";
    }
    return 0;
  }
  // TEMP T-123: `--facefix <faceIdx>` runs the real `ShapeFix_Face::Perform`
  // (all modes default) on one face of the loaded shape and prints, before and
  // after, every wire's edge count plus each junction `CheckLacking` verdict
  // (`ShapeAnalysis_Wire::CheckLacking`), so the port's second round
  // (`perform_wire_second_round` -> `fix_lacking_all`) can be compared against
  // what the same C++ call decides on the same wire.
  if (argc > 3 && std::string(argv[2]) == "--facefix")
  {
    const int aFaceIdx = std::atoi(argv[3]);
    TopoDS_Face aFace;
    int         k = 0;
    for (TopExp_Explorer fx(aShape, TopAbs_FACE); fx.More(); fx.Next(), ++k)
    {
      if (k == aFaceIdx)
      {
        aFace = TopoDS::Face(fx.Current());
        break;
      }
    }
    if (aFace.IsNull())
    {
      std::cout << "FF no-face\n";
      return 0;
    }

    auto report = [&](const char* tag, const TopoDS_Face& f) {
      int nw = 0;
      for (TopExp_Explorer wx(f, TopAbs_WIRE); wx.More(); wx.Next())
      {
        ++nw;
        occ::handle<ShapeExtend_WireData> wd = new ShapeExtend_WireData(TopoDS::Wire(wx.Current()));
        std::cout << "FF " << tag << " w" << nw << " ori=" << (int)wx.Current().Orientation()
                  << " n=" << wd->NbEdges() << " cross=" << ShapeAnalysis::TotCross2D(wd, f)
                  << "\n";
        ShapeAnalysis_Wire aWA(wd, f, Precision::Confusion());
        for (int i = 1; i <= wd->NbEdges(); ++i)
        {
          const TopoDS_Edge&        e  = wd->Edge(i);
          double                    f1 = 0., l1 = 0.;
          occ::handle<Geom2d_Curve> c2 = BRep_Tool::CurveOnSurface(e, f, f1, l1);
          gp_Pnt2d                  b(0., 0.), en(0., 0.);
          if (!c2.IsNull())
          {
            b  = c2->Value(f1);
            en = c2->Value(l1);
          }
          gp_Pnt2d p1(0., 0.), p2(0., 0.);
          const bool lack = aWA.CheckLacking(i, 0., p1, p2);
          std::cout << "FF   e" << i << " ori=" << (int)e.Orientation()
                    << " closed=" << (BRep_Tool::IsClosed(e, f) ? 1 : 0) << " uv=(" << b.X() << ","
                    << b.Y() << ")->(" << en.X() << "," << en.Y() << ")";
          if (lack)
          {
            std::cout << " lack=1 d2=" << (aWA.LastCheckStatus(ShapeExtend_DONE2) ? 1 : 0)
                      << " max2d=" << aWA.MaxDistance2d() << " max3d=" << aWA.MaxDistance3d()
                      << " p2d1=(" << p1.X() << "," << p1.Y() << ") p2d2=(" << p2.X() << ","
                      << p2.Y() << ")";
          }
          std::cout << "\n";
        }
      }
    };

    report("pre", aFace);

    occ::handle<ShapeBuild_ReShape> aCtx = new ShapeBuild_ReShape;
    ShapeFix_Face                   aSFF(aFace);
    aSFF.SetContext(aCtx);
    aSFF.SetPrecision(Precision::Confusion());
    const bool aDone = aSFF.Perform(Message_ProgressRange());
    std::cout << "FF perform=" << (aDone ? 1 : 0) << "\n";
    if (!aSFF.FixWireTool().IsNull())
    {
      occ::handle<ShapeFix_Wire> aFW = aSFF.FixWireTool();
      std::cout << "FF status n=" << aFW->NbEdges()
                << " lacking=" << (aFW->StatusLacking(ShapeExtend_DONE) ? 1 : 0)
                << " reorder=" << (aFW->StatusReorder(ShapeExtend_DONE) ? 1 : 0)
                << " closed=" << (aFW->StatusClosed(ShapeExtend_DONE) ? 1 : 0) << "\n";
    }
    const TopoDS_Shape aRes = aSFF.Result();
    if (aRes.IsNull() || aRes.IsSame(aFace))
    {
      std::cout << "FF result-same\n";
      report("post", aFace);
    }
    else
    {
      report("post", TopoDS::Face(aRes));
    }
    return 0;
  }
  // TEMP T-121: `--wirefix <faceIdx>` replays `ShapeFix_Face::Perform`'s wire
  // pass (`ShapeFix_Face.cxx:400-455`) on one face of the loaded (usually
  // `--nofix`) shape: for the full pass and then for every pass with one
  // `ShapeFix_Wire` mode disabled, it prints `ShapeAnalysis::TotCross2D`
  // (`ShapeAnalysis.cxx:115-150`) and the resulting edge orientations, so the
  // port's `perform_wire_first_round` can be compared stage by stage.
  // `--wirefixraw` is the same replay on the wire as `ShapeFix_Face::Perform`
  // sees it: `TopoDS_Iterator(face, false)` does not compose the face
  // orientation (`ShapeFix_Face.cxx:401`), while `TopExp_Explorer` does, so
  // `ShapeExtend_WireData` prepends here and appends there.
  if (argc > 3
      && (std::string(argv[2]) == "--wirefix" || std::string(argv[2]) == "--wirefixraw"))
  {
    const bool aRawOri = (std::string(argv[2]) == "--wirefixraw");
    const int  aFaceIdx = std::atoi(argv[3]);
    TopoDS_Face aFace;
    int         k = 0;
    for (TopExp_Explorer fx(aShape, TopAbs_FACE); fx.More(); fx.Next(), ++k)
    {
      if (k == aFaceIdx)
      {
        aFace = TopoDS::Face(fx.Current());
        break;
      }
    }
    if (aFace.IsNull())
    {
      std::cout << "WF no-face\n";
      return 0;
    }
    occ::handle<Geom_Surface>          aSurf = BRep_Tool::Surface(aFace);
    occ::handle<ShapeAnalysis_Surface> aSAS  = new ShapeAnalysis_Surface(aSurf);
    TopoDS_Wire                        aW;
    if (aRawOri)
    {
      for (TopoDS_Iterator aIt(aFace, false); aIt.More(); aIt.Next())
      {
        if (aIt.Value().ShapeType() == TopAbs_WIRE)
        {
          aW = TopoDS::Wire(aIt.Value());
          break;
        }
      }
      std::cout << "WF face-ori=" << (int)aFace.Orientation()
                << " wire-ori=" << (int)aW.Orientation() << "\n";
      for (TopoDS_Iterator aE(aW, false); aE.More(); aE.Next())
      {
        std::cout << "WF stored-ori=" << (int)aE.Value().Orientation() << "\n";
      }
    }
    else
    {
      aW = TopoDS::Wire(TopExp_Explorer(aFace, TopAbs_WIRE).Current());
    }

    auto report = [&](const char* tag, const TopoDS_Wire& w) {
      occ::handle<ShapeExtend_WireData> aWD = new ShapeExtend_WireData(w);
      std::cout << "WF " << tag << " n=" << aWD->NbEdges()
                << " cross=" << ShapeAnalysis::TotCross2D(aWD, aFace) << "\n";
      for (int i = 1; i <= aWD->NbEdges(); ++i)
      {
        const TopoDS_Edge&        e = aWD->Edge(i);
        double                    f = 0., l = 0.;
        occ::handle<Geom2d_Curve> c2 = BRep_Tool::CurveOnSurface(e, aFace, f, l);
        gp_Pnt2d                  b(0., 0.), en(0., 0.);
        if (!c2.IsNull())
        {
          b  = c2->Value(f);
          en = c2->Value(l);
        }
        std::cout << "WF   e" << i << " ori=" << (int)e.Orientation() << " uv=(" << b.X() << ","
                  << b.Y() << ")->(" << en.X() << "," << en.Y() << ")\n";
      }
    };

    report("raw", aW);

    {
      ShapeAnalysis_Edge aEA;
      occ::handle<ShapeExtend_WireData> aWD = new ShapeExtend_WireData(aW);
      for (int i = 1; i <= aWD->NbEdges(); ++i)
      {
        const TopoDS_Edge   e = aWD->Edge(i);
        const TopoDS_Vertex v1 = aEA.FirstVertex(e);
        const TopoDS_Vertex v2 = aEA.LastVertex(e);
        int                 aNbV = 0;
        for (TopoDS_Iterator aVi(e, false); aVi.More(); aVi.Next())
        {
          ++aNbV;
        }
        std::cout << "WF v" << i << " children=" << aNbV
                  << " same=" << ((!v1.IsNull() && !v2.IsNull() && v1.IsSame(v2)) ? 1 : 0)
                  << " null=" << (v1.IsNull() || v2.IsNull() ? 1 : 0);
        if (!v1.IsNull() && !v2.IsNull())
        {
          const gp_Pnt p1 = BRep_Tool::Pnt(v1);
          const gp_Pnt p2 = BRep_Tool::Pnt(v2);
          std::cout << " p1=(" << p1.X() << "," << p1.Y() << "," << p1.Z() << ") p2=(" << p2.X()
                    << "," << p2.Y() << "," << p2.Z() << ")";
        }
        std::cout << "\n";
      }
    }

    {
      ShapeAnalysis_Wire      aWT(new ShapeExtend_WireData(aW), aFace, Precision::Confusion());
      ShapeAnalysis_WireOrder aO3;
      aWT.CheckOrder(aO3, true, true);
      std::cout << "WF order3d nb=" << aO3.NbEdges() << " status=" << aO3.Status() << " ordered=[";
      for (int i = 1; i <= aO3.NbEdges(); ++i)
      {
        std::cout << (i > 1 ? " " : "") << aO3.Ordered(i);
      }
      std::cout << "]\n";
      ShapeAnalysis_WireOrder aO2;
      aWT.CheckOrder(aO2, true, false);
      std::cout << "WF order2d nb=" << aO2.NbEdges() << " status=" << aO2.Status() << " ordered=[";
      for (int i = 1; i <= aO2.NbEdges(); ++i)
      {
        std::cout << (i > 1 ? " " : "") << aO2.Ordered(i);
      }
      std::cout << "]\n";
    }

    const char* aModes[] = {"ALL",          "FixReorder",  "FixSmall",     "FixConnected",
                            "FixEdgeCurves", "FixDegenerated", "FixTails",
                            "FixShifted",   "FixAddPCurve", "FixAddCurve3d", "FixReversed2d",
                            "FixSeam",      "ROff+COn",     "ROff+SmallOn"};
    const int   aNbModes = (int)(sizeof(aModes) / sizeof(aModes[0]));
    for (int m = 0; m < aNbModes; ++m)
    {
      ShapeFix_Wire aSFW;
      aSFW.SetFace(aFace, aSAS);
      aSFW.Load(aW);
      aSFW.FixLackingMode()          = false;
      aSFW.FixSelfIntersectionMode() = false;
      const std::string mn = aModes[m];
      if (mn == "FixReorder")
        aSFW.FixReorderMode() = 0;
      else if (mn == "FixSmall")
        aSFW.FixSmallMode() = 0;
      else if (mn == "FixConnected")
        aSFW.FixConnectedMode() = 0;
      else if (mn == "FixEdgeCurves")
        aSFW.FixEdgeCurvesMode() = 0;
      else if (mn == "FixDegenerated")
        aSFW.FixDegeneratedMode() = 0;
      else if (mn == "FixTails")
        aSFW.FixTailMode() = 0;
      else if (mn == "FixShifted")
        aSFW.FixShiftedMode() = 0;
      else if (mn == "FixAddPCurve")
        aSFW.FixAddPCurveMode() = 0;
      else if (mn == "FixAddCurve3d")
        aSFW.FixAddCurve3dMode() = 0;
      else if (mn == "FixReversed2d")
        aSFW.FixReversed2dMode() = 0;
      else if (mn == "FixSeam")
        aSFW.FixSeamMode() = 0;
      else if (mn == "ROff+COn")
      {
        aSFW.FixReorderMode()   = 0;
        aSFW.FixConnectedMode() = 1;
      }
      else if (mn == "ROff+SmallOn")
      {
        aSFW.FixReorderMode() = 0;
        aSFW.FixSmallMode()   = 1;
      }
      const bool aRes = aSFW.Perform(Message_ProgressRange());
      report(mn.c_str(), aSFW.Wire());
      std::cout << "WF   res=" << (aRes ? 1 : 0) << "\n";
    }
    return 0;
  }
  return std::nullopt;
}
