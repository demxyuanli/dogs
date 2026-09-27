// Read-only probe: per-face wire census + UV-degenerate-wire census on OCCT's shape.
// A wire is "UV-degenerate" when its pcurve UV bounding box has ~zero area, which is
// the signature the port's mesh pipeline chokes on (zero-area loops -> open chains).
#include <BRepAdaptor_Surface.hxx>
#include <BRepBndLib.hxx>
#include <Bnd_Box.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRep_Tool.hxx>
#include <Geom2d_Curve.hxx>
#include <STEPControl_Reader.hxx>
#include <Interface_Static.hxx>
#include <ShapeProcess.hxx>
#include <ShapeFix_Face.hxx>
#include <ShapeFix_Shape.hxx>
#include <ShapeBuild_ReShape.hxx>
#include <ShapeFix.hxx>
#include <XSAlgo_ShapeProcessor.hxx>
#include <TCollection_AsciiString.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Wire.hxx>
#include <Poly_Triangulation.hxx>
#include <TopLoc_Location.hxx>
#include <cstdlib>
#include <iostream>
#include <map>
#include <cmath>

int main(int argc, char** argv)
{
  if (argc < 2)
  {
    std::cout << "usage: wires_probe <file.step> [faceIndex]\n";
    return 2;
  }
  // "raw <...>": disable the FromSTEP FixShape operator to inspect the shape as
  // transferred, before any healing (oracle for the port's 4-vs-2 wire question).
  if (argc > 2 && std::string(argv[2]) == "raw")
  {
    Interface_Static::SetCVal("FromSTEP.exec.op", "");
  }

  STEPControl_Reader aReader;
  if (aReader.ReadFile(argv[1]) != IFSelect_RetDone)
  {
    std::cout << "read failed\n";
    return 1;
  }
  // "setp <key> <value>": set one entry of the actor's ShapeFix parameter map
  // (bypasses the text eval that Interface_Static::SetCVal breaks on).
  if (argc > 4 && std::string(argv[2]) == "setp")
  {
    XSAlgo_ShapeProcessor::ParameterMap aMap = aReader.GetShapeFixParameters();
    XSAlgo_ShapeProcessor::SetParameter(argv[3], TCollection_AsciiString(argv[4]), true, aMap);
    aReader.SetShapeFixParameters(std::move(aMap));
  }

  // "set"/"setc <static> <value>": tweak one Interface_Static AFTER ReadFile (the
  // controller registers the FromSTEP.FixShape.* statics as type 't' lazily; setting
  // them before the reader exists corrupts them and crashes the probe).
  if (argc > 4 && std::string(argv[2]) == "set")
  {
    Interface_Static::SetCVal(argv[3], argv[4]);
  }
  if (argc > 4 && std::string(argv[2]) == "setc")
  {
    Interface_Static::SetCVal(argv[3], argv[4]);
  }

  // "noop": clear the actor's ShapeProcess operations, so the per-entity
  // ProcessShape (STEPControl_ActorRead.cxx:2150) becomes a no-op and the
  // transferred shape is seen before any healing. Must run after ReadFile
  // (which creates the workspace/actor), before TransferRoots.
  if ((argc > 2 && std::string(argv[2]) == "noop")
      || (argc > 4 && std::string(argv[4]) == "noop"))
  {
    aReader.SetShapeProcessFlags(ShapeProcess::OperationsFlags());
  }
  // "oneop <index>": enable exactly one ShapeProcess operation, to identify
  // which operator merges the multi-bound wires.
  if (argc > 3 && std::string(argv[2]) == "oneop")
  {
    ShapeProcess::OperationsFlags f;
    f.set(std::atoi(argv[3]));
    aReader.SetShapeProcessFlags(f);
  }
  aReader.TransferRoots();
  const TopoDS_Shape aShape = aReader.OneShape();
  const int aWant = (argc > 2) ? std::atoi(argv[2]) : -1;

  // "wdump <faceIndex>": per-wire edge census with each edge's 3D endpoints, to
  // see how the transfer's wires and the healed wires relate.
  if (argc > 3 && std::string(argv[2]) == "wdump")
  {
    const int idx = std::atoi(argv[3]);
    TopExp_Explorer ex(aShape, TopAbs_FACE);
    int k = 0;
    for (; ex.More(); ex.Next(), ++k)
    {
      if (k != idx) continue;
      TopoDS_Face F = TopoDS::Face(ex.Current());
      int wi = 0;
      for (TopoDS_Iterator it(F, false); it.More(); it.Next(), ++wi)
      {
        if (it.Value().ShapeType() != TopAbs_WIRE) continue;
        std::cout << "W" << wi << " edges=";
        int ei = 0;
        for (TopoDS_Iterator et(it.Value(), false); et.More(); et.Next(), ++ei)
        {
          TopoDS_Edge E = TopoDS::Edge(et.Value());
          double a, b;
          occ::handle<Geom_Curve> c = BRep_Tool::Curve(E, a, b);
          if (c.IsNull()) { std::cout << "[?]"; continue; }
          gp_Pnt p0 = c->Value(a), p1 = c->Value(b);
          std::cout << " [" << p0.X() << "," << p0.Y() << "," << p0.Z() << "->"
                    << p1.X() << "," << p1.Y() << "," << p1.Z() << "]";
        }
        std::cout << "\n";
      }
      break;
    }
    return 0;
  }

  // "modes <faceIndex> [noop]": compare the mode fields of a standalone
  // ShapeFix_Face against the one owned by a ShapeFix_Shape.
  if (argc > 3 && std::string(argv[2]) == "modes")
  {
    const int idx = std::atoi(argv[3]);
    TopExp_Explorer ex(aShape, TopAbs_FACE);
    int k = 0;
    auto dump = [](const char* tag, const occ::handle<ShapeFix_Face>& f) {
      std::cout << tag
        << " wire=" << f->FixWireMode()
        << " miss=" << f->FixMissingSeamMode()
        << " orient=" << f->FixOrientationMode()
        << " natural=" << f->FixAddNaturalBoundMode()
        << " small=" << f->FixSmallAreaWireMode()
        << " intersect=" << f->FixIntersectingWiresMode()
        << " loop=" << f->FixLoopWiresMode()
        << " split=" << f->FixSplitFaceMode()
        << " wtopo=" << f->FixWireTool()->ModifyTopologyMode()
        << " wclosed=" << f->FixWireTool()->ClosedWireMode()
        << " wsmall=" << f->FixWireTool()->FixSmallMode()
        << " wconn=" << f->FixWireTool()->FixConnectedMode()
        << " wlacking=" << f->FixWireTool()->FixLackingMode()
        << " wself=" << f->FixWireTool()->FixSelfIntersectionMode()
        << "\n";
    };
    for (; ex.More(); ex.Next(), ++k)
    {
      if (k != idx) continue;
      TopoDS_Face F = TopoDS::Face(ex.Current());
      ShapeFix_Shape a;
      occ::handle<ShapeFix_Face> fa = a.FixFaceTool();
      ShapeFix_Shape b;
      b.Init(F);
      occ::handle<ShapeFix_Face> fb = b.FixFaceTool();
      dump("MODE standalone ", fa);
      dump("MODE shapefix   ", fb);
      break;
    }
    return 0;
  }

  // "emulm <faceIndex> [noop] <method>": run one ShapeFix_Face method in the
  // ShapeFix_Shape-owned tool and report the resulting wire count.
  if (argc > 3 && std::string(argv[2]) == "emulm")
  {
    const int idx = std::atoi(argv[3]);
    std::string m = (argc > 5) ? argv[5] : "none";
    TopExp_Explorer ex(aShape, TopAbs_FACE);
    int k = 0;
    auto cnt = [](const TopoDS_Face& f) {
      int w = 0;
      for (TopoDS_Iterator it(f, false); it.More(); it.Next())
        if (it.Value().ShapeType() == TopAbs_WIRE) ++w;
      return w;
    };
    for (; ex.More(); ex.Next(), ++k)
    {
      if (k != idx) continue;
      TopoDS_Face F = TopoDS::Face(ex.Current());
      ShapeFix_Shape sfs;
      sfs.Init(F);
      occ::handle<ShapeFix_Face> sff = sfs.FixFaceTool();
      sff->Init(TopoDS::Face(sfs.Context()->Apply(F)));
      sff->SetContext(sfs.Context());
      bool r = false;
      if (m == "miss") r = sff->FixMissingSeam();
      else if (m == "orient") r = sff->FixOrientation();
      else if (m == "natural") r = sff->FixAddNaturalBound();
      else if (m == "small") r = sff->FixSmallAreaWire(true);
      else if (m == "coinc") r = sff->FixWiresTwoCoincEdges();
      else if (m == "intersect") r = sff->FixIntersectingWires();
      else if (m == "periodic") r = sff->FixPeriodicDegenerated();
      std::cout << "EMULM " << m << " ret=" << r << " face_wires=" << cnt(sff->Face())
                << " result_wires=" << (sff->Result().IsNull() ? -1 : (int)sff->Result().NbChildren()) << "\n";
      break;
    }
    return 0;
  }

  // "emul <faceIndex> [noop]": emulate ShapeFix_Shape::Perform step by step.
  if (argc > 3 && std::string(argv[2]) == "emul")
  {
    const int idx = std::atoi(argv[3]);
    TopExp_Explorer ex(aShape, TopAbs_FACE);
    int k = 0;
    auto cnt = [](const TopoDS_Shape& s) {
      int w = 0;
      if (!s.IsNull())
        for (TopExp_Explorer e(s, TopAbs_WIRE); e.More(); e.Next()) ++w;
      return w;
    };
    for (; ex.More(); ex.Next(), ++k)
    {
      if (k != idx) continue;
      TopoDS_Face F = TopoDS::Face(ex.Current());
      ShapeFix_Shape sfs;
      sfs.Init(F);
      TopoDS_Shape S = sfs.Context()->Apply(F);
      std::cout << "EMUL init wires=" << cnt(S) << " same=" << S.IsSame(F)
                << " sori=" << (int)S.Orientation() << " fori=" << (int)F.Orientation()
                << " schild=" << S.NbChildren() << " fchild=" << F.NbChildren() << "\n";
      occ::handle<ShapeFix_Face> sff = sfs.FixFaceTool();
      sff->Init(TopoDS::Face(S));
      sff->SetContext(sfs.Context());
      sff->FixWireTool()->ModifyTopologyMode() = true;
      if (argc > 4) {
        std::string m4 = argv[4];
        if (m4 == "nowire") sff->FixWireMode() = 0;
        if (m4 == "nomiss") sff->FixMissingSeamMode() = 0;
        if (m4 == "noorient") sff->FixOrientationMode() = 0;
        if (m4 == "noloop") sff->FixLoopWiresMode() = 0;
        if (m4 == "nointersect") sff->FixIntersectingWiresMode() = 0;
        if (m4 == "nosplit") sff->FixSplitFaceMode() = 0;
        if (m4 == "nosmall") sff->FixSmallAreaWireMode() = 0;
        if (m4 == "nonatural") sff->FixAddNaturalBoundMode() = 0;
      }
      if (argc > 4 && std::string(argv[4]) == "fresh") sff->SetContext(new ShapeBuild_ReShape);
      if (argc > 4 && std::string(argv[4]) == "notopo") sff->FixWireTool()->ModifyTopologyMode() = false;
      if (argc > 4 && std::string(argv[4]) == "noctx") sff->SetContext(occ::handle<ShapeBuild_ReShape>());
      if (argc > 4 && std::string(argv[4]) == "local")
      {
        ShapeFix_Face sl;
        sl.Init(TopoDS::Face(S));
        sl.SetContext(sfs.Context());
        sl.FixWireTool()->ModifyTopologyMode() = true;
        sl.Perform();
        int ff = 0;
        for (TopExp_Explorer e2(sl.Result(), TopAbs_FACE); e2.More(); e2.Next()) ++ff;
        std::cout << "EMUL local faces=" << ff << " wires=" << cnt(sl.Result()) << "\n";
        break;
      }
      sff->Perform();
      {
        int ff = 0;
        for (TopExp_Explorer e2(sff->Result(), TopAbs_FACE); e2.More(); e2.Next()) ++ff;
        std::cout << "EMUL face faces=" << ff << " wires=" << cnt(sff->Result()) << "\n";
      }
      TopoDS_Shape mr = sfs.Context()->Apply(S);
      std::cout << "EMUL apply wires=" << cnt(mr) << "\n";
      ShapeFix::SameParameter(mr, false, 0.0);
      std::cout << "EMUL same wires=" << cnt(mr) << "\n";
      break;
    }
    return 0;
  }

  // "shp <faceIndex> [noop]": run ShapeFix_Shape::Perform on one face.
  if (argc > 3 && std::string(argv[2]) == "shp")
  {
    const int idx = std::atoi(argv[3]);
    TopExp_Explorer ex(aShape, TopAbs_FACE);
    int k = 0;
    for (; ex.More(); ex.Next(), ++k)
    {
      if (k != idx) continue;
      TopoDS_Face F = TopoDS::Face(ex.Current());
      int w0 = 0;
      for (TopoDS_Iterator it(F, false); it.More(); it.Next())
        if (it.Value().ShapeType() == TopAbs_WIRE) ++w0;
      ShapeFix_Shape sfs;
      sfs.Init(F);
      if (argc > 5 && std::string(argv[5]) == "nosame") sfs.FixSameParameterMode() = 0;
      if (argc > 5 && std::string(argv[5]) == "novtol") sfs.FixVertexTolMode() = 0;
      bool r = sfs.Perform();
      TopoDS_Shape res = sfs.Shape();
      int faces = 0, wires = 0;
      if (!res.IsNull())
      {
        for (TopExp_Explorer e2(res, TopAbs_FACE); e2.More(); e2.Next()) ++faces;
        for (TopExp_Explorer e2(res, TopAbs_WIRE); e2.More(); e2.Next()) ++wires;
      }
      std::cout << "SHP idx=" << idx << " in_wires=" << w0 << " ret=" << r
                << " null=" << res.IsNull() << " faces=" << faces << " wires=" << wires << "\n";
      break;
    }
    return 0;
  }

  // "perf <faceIndex> [noop]": run the full ShapeFix_Face::Perform on one face.
  if (argc > 3 && std::string(argv[2]) == "perf")
  {
    const int idx = std::atoi(argv[3]);
    TopExp_Explorer ex(aShape, TopAbs_FACE);
    int k = 0;
    for (; ex.More(); ex.Next(), ++k)
    {
      if (k != idx) continue;
      TopoDS_Face F = TopoDS::Face(ex.Current());
      int w0 = 0;
      for (TopoDS_Iterator it(F, false); it.More(); it.Next())
        if (it.Value().ShapeType() == TopAbs_WIRE) ++w0;
      ShapeFix_Shape tmpS;
      if (argc > 4 && std::string(argv[4]) == "init") tmpS.Init(F);
      occ::handle<ShapeFix_Face> sffH = tmpS.FixFaceTool();
      ShapeFix_Face* sffp = sffH.operator->();
      if (argc > 4 && std::string(argv[4]) == "apply")
      {
        tmpS.Init(F);
        sffp->Init(TopoDS::Face(tmpS.Context()->Apply(F)));
      }
      if (argc > 4 && std::string(argv[4]) == "applyboth")
      {
        tmpS.Init(F);
        sffp->Init(TopoDS::Face(tmpS.Context()->Apply(F)));
        sffp->SetContext(tmpS.Context());
        sffp->FixWireTool()->ModifyTopologyMode() = true;
      }
      else
      {
        sffp->Init(F);
      }
      if (argc > 4 && std::string(argv[4]) == "topo") sffp->FixWireTool()->ModifyTopologyMode() = true;
      if (argc > 4 && std::string(argv[4]) == "ctx") sffp->SetContext(new ShapeBuild_ReShape);
      if (argc > 4 && std::string(argv[4]) == "both") { sffp->SetContext(new ShapeBuild_ReShape); sffp->FixWireTool()->ModifyTopologyMode() = true; }
      if (argc > 4 && std::string(argv[4]) == "mloc") { occ::handle<ShapeBuild_ReShape> c = new ShapeBuild_ReShape; c->ModeConsiderLocation() = true; sffp->SetContext(c); sffp->FixWireTool()->ModifyTopologyMode() = true; }
      if (argc > 4 && std::string(argv[4]) == "mloc2") { occ::handle<ShapeBuild_ReShape> c = new ShapeBuild_ReShape; c->ModeConsiderLocation() = true; sffp->SetContext(c); }
      bool r = sffp->Perform();
      TopoDS_Shape res = sffp->Result();
      int faces = 0, wires = 0;
      if (!res.IsNull())
      {
        for (TopExp_Explorer e2(res, TopAbs_FACE); e2.More(); e2.Next()) ++faces;
        for (TopExp_Explorer e2(res, TopAbs_WIRE); e2.More(); e2.Next()) ++wires;
      }
      std::cout << "PERF idx=" << idx << " in_wires=" << w0 << " ret=" << r
                << " null=" << res.IsNull() << " faces=" << faces << " wires=" << wires << "\n";
      break;
    }
    return 0;
  }

  // "fms <faceIndex>": run ShapeFix_Face::FixMissingSeam on one face and report
  // the result's type / face / wire counts (OCCT oracle for the port's Shell).
  if (argc > 3 && std::string(argv[2]) == "fms")
  {
    const int idx = std::atoi(argv[3]);
    TopExp_Explorer ex(aShape, TopAbs_FACE);
    int k = 0;
    for (; ex.More(); ex.Next(), ++k)
    {
      if (k != idx) continue;
      TopoDS_Face F = TopoDS::Face(ex.Current());
      if (BRep_Tool::Surface(F).IsNull()) { std::cout << "FMS no surface\n"; return 0; }
      ShapeFix_Face sff;
      sff.Init(F);
      bool r = sff.FixMissingSeam();
      TopoDS_Shape res = sff.Result();
      int faces = 0, wires = 0;
      if (!res.IsNull())
      {
        for (TopExp_Explorer e2(res, TopAbs_FACE); e2.More(); e2.Next()) ++faces;
        for (TopExp_Explorer e2(res, TopAbs_WIRE); e2.More(); e2.Next()) ++wires;
      }
      std::cout << "FMS idx=" << idx << " ret=" << r << " null=" << res.IsNull()
                << " type=" << (res.IsNull() ? -1 : (int)res.ShapeType())
                << " faces=" << faces << " wires=" << wires << "\n";
      break;
    }
    return 0;
  }



  TopExp_Explorer anEx(aShape, TopAbs_FACE);
  int aIdx = 0, aMul = 0, aDegen = 0, aTotalWires = 0;
  std::map<int, int> aHist;
  for (; anEx.More(); anEx.Next(), ++aIdx)
  {
    const TopoDS_Face& aFace = TopoDS::Face(anEx.Current());
    int aWires = 0;
    std::string aPerWire;
    std::string aPerWireSpan;
    bool aFaceDegen = false;
    for (TopExp_Explorer aW(aFace, TopAbs_WIRE); aW.More(); aW.Next())
    {
      ++aWires;
      int    aEdges = 0;
      double u0 = 1e100, u1 = -1e100, v0 = 1e100, v1 = -1e100;
      for (TopExp_Explorer aE(aW.Current(), TopAbs_EDGE); aE.More(); aE.Next())
      {
        ++aEdges;
        const TopoDS_Edge& anEdge = TopoDS::Edge(aE.Current());
        double             f = 0., l = 0.;
        const occ::handle<Geom2d_Curve> aC2d = BRep_Tool::CurveOnSurface(anEdge, aFace, f, l);
        if (aC2d.IsNull()) continue;
        if (aIdx == aWant)
        {
          const gp_Pnt2d aP0 = aC2d->Value(f);
          const gp_Pnt2d aP1 = aC2d->Value(l);
          std::cout << "  EDGE f=" << f << " l=" << l
                    << " p0=(" << aP0.X() << "," << aP0.Y() << ")"
                    << " p1=(" << aP1.X() << "," << aP1.Y() << ")"
                    << " ori=" << anEdge.Orientation()
                    << " deg=" << (BRep_Tool::Degenerated(anEdge) ? 1 : 0)
                    << " closedOnFace=" << (BRep_Tool::IsClosed(anEdge, aFace) ? 1 : 0) << "\n";
        }
        if (BRep_Tool::IsClosed(anEdge, aFace))
        {
          double f2 = 0., l2 = 0.;
          const TopoDS_Edge aRev = TopoDS::Edge(anEdge.Reversed());
          const occ::handle<Geom2d_Curve> aC2d2 = BRep_Tool::CurveOnSurface(aRev, aFace, f2, l2);
          const gp_Pnt2d a = aC2d->Value(f);
          const gp_Pnt2d b = aC2d->Value(l);
          std::cout << "SEAM face=" << aIdx << " c1=(" << f << "," << l << ")"
                    << " p1=(" << a.X() << "," << a.Y() << ")->(" << b.X() << "," << b.Y() << ")";
          if (!aC2d2.IsNull())
          {
            const gp_Pnt2d c = aC2d2->Value(f2);
            const gp_Pnt2d d = aC2d2->Value(l2);
            std::cout << " c2=(" << f2 << "," << l2 << ")"
                      << " p2=(" << c.X() << "," << c.Y() << ")->(" << d.X() << "," << d.Y() << ")";
          }
          std::cout << "\n";
        }
        for (int k = 0; k <= 8; ++k)
        {
          const gp_Pnt2d aP = aC2d->Value(f + (l - f) * k / 8.0);
          u0 = std::min(u0, aP.X()); u1 = std::max(u1, aP.X());
          v0 = std::min(v0, aP.Y()); v1 = std::max(v1, aP.Y());
        }
      }
      aPerWire += " " + std::to_string(aEdges);
      const double aSpanU = (u1 > u0) ? (u1 - u0) : 0.0;
      const double aSpanV = (v1 > v0) ? (v1 - v0) : 0.0;
      char aBuf[96];
      std::snprintf(aBuf, sizeof(aBuf), " u=%.4f v=%.4f", aSpanU, aSpanV);
      aPerWireSpan += aBuf;
      if (aSpanU < 1e-9 || aSpanV < 1e-9) aFaceDegen = true;
    }
    aTotalWires += aWires;
    aHist[aWires]++;
    if (aWires >= 2) ++aMul;
    if (aFaceDegen) ++aDegen;
    if (aIdx == aWant || true)
    {
      BRepAdaptor_Surface aSurf(aFace);
      Bnd_Box aB; BRepBndLib::Add(aFace, aB);
      double bx0=0,by0=0,bz0=0,bx1=0,by1=0,bz1=0; aB.Get(bx0,by0,bz0,bx1,by1,bz1);
      std::cout << "FACE " << aIdx << " type=" << aSurf.GetType() << " wires=" << aWires
                << " edges_per_wire:" << aPerWire << " spans:" << aPerWireSpan
                << " box=(" << bx0 << "," << by0 << "," << bz0 << ")-(" << bx1 << "," << by1 << "," << bz1 << ")"
                << (aFaceDegen ? "  [UV-DEGENERATE-WIRE]" : "") << "\n";
    }
  }
  std::cout << "TOTAL faces=" << aIdx << " wires=" << aTotalWires << " multiwire_faces=" << aMul
            << " uv_degenerate_faces=" << aDegen << " hist:";
  for (const auto& kv : aHist) std::cout << " " << kv.first << "->" << kv.second;
  std::cout << "\n";

  // Per-face OCCT mesh census at the port's linear deflection (default
  // 1.076007 for a3n00 via prs3d_get_deflection(shape, 0.1)).
  const double aLin = (argc > 3) ? std::atof(argv[3]) : 1.076007;
  BRepMesh_IncrementalMesh aMesh(aShape, aLin, false, 20. * M_PI / 180., true);
  TopExp_Explorer aMF(aShape, TopAbs_FACE);
  int aMI = 0;
  for (; aMF.More(); aMF.Next(), ++aMI)
  {
    const TopoDS_Face& aFace = TopoDS::Face(aMF.Current());
    TopLoc_Location aLoc;
    const occ::handle<Poly_Triangulation>& aTri = BRep_Tool::Triangulation(aFace, aLoc);
    if (aTri.IsNull())
    {
      std::cout << "MF " << aMI << " null\n";
      continue;
    }
    BRepAdaptor_Surface aS(aFace);
    Bnd_Box             aBox;
    BRepBndLib::Add(aFace, aBox);
    double bx0 = 0, by0 = 0, bz0 = 0, bx1 = 0, by1 = 0, bz1 = 0;
    if (!aBox.IsVoid())
      aBox.Get(bx0, by0, bz0, bx1, by1, bz1);
    std::cout << "MF " << aMI << " type=" << aS.GetType() << " nodes=" << aTri->NbNodes()
              << " tris=" << aTri->NbTriangles()
              << " box=(" << bx0 << "," << by0 << "," << bz0 << ")-(" << bx1 << "," << by1 << "," << bz1 << ")\n";
  }
  return 0;
}
