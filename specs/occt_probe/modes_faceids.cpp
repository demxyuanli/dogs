// `argv[2]` modes: `--faceids`, `--fid` (STEP entity ids and per-face detail).
#include "occt_probe_common.h"

std::optional<int> runFaceIdModes(const ProbeCtx& theCtx)
{
  const int           argc = theCtx.argc;
  char** const        argv = theCtx.argv;
  const TopoDS_Shape& aShape = theCtx.shape;
  STEPControl_Reader& aReader = theCtx.reader;

  // T-69 oracle: wire/edge structure of the *imported* shape. Answers "how many
  // wires does OCCT's import produce for a face whose port twin has two wires of
  // coincident closed edges". Read-only.
  if (argc > 2 && std::string(argv[2]) == "--faceids")
  {
    occ::handle<XSControl_TransferReader>  tr = aReader.WS()->TransferReader();
    occ::handle<Interface_InterfaceModel>  model = aReader.Model();
    int k = 0;
    for (TopExp_Explorer ex(aShape, TopAbs_FACE); ex.More(); ex.Next(), ++k)
    {
      occ::handle<Standard_Transient> ent = tr->EntityFromShapeResult(ex.Current(), 1);
      int id = 0;
      if (!ent.IsNull())
      {
        id = model->Number(ent);
      }
      int nw = 0;
      for (TopExp_Explorer we(ex.Current(), TopAbs_WIRE); we.More(); we.Next()) ++nw;
      // TEMP T-102: attach the face bbox so an id=0 (created by healing) face can
      // be paired with the port's `--ecensus` bbox.
      Bnd_Box aBox;
      BRepBndLib::Add(TopoDS::Face(ex.Current()), aBox, false);
      double v[6] = {0, 0, 0, 0, 0, 0};
      if (!aBox.IsVoid())
      {
        aBox.Get(v[0], v[1], v[2], v[3], v[4], v[5]);
      }
      std::cout << std::fixed << std::setprecision(6) << "STEPFACE face=" << k << " id=" << id
                << " wires=" << nw << " bbox=(" << v[0] << "," << v[1] << "," << v[2] << ")-("
                << v[3] << "," << v[4] << "," << v[5] << ")\n";
    }
    return 0;
  }

  // TEMP T-104: `--fid <i0>[:<i1>]` dumps the same per-face detail as `--fbox`
  // but selects by INDEX. `--fbox` accepts a match within 0.5 of the query box,
  // which is ambiguous when two faces share a box (the failure mode behind the
  // "cone vs cylinder" reading), and the transferred entity id printed here is
  // the only key that survives healing renumbering.
  if (argc > 3 && std::string(argv[2]) == "--fid")
  {
    int i0 = 0, i1 = 0;
    {
      const std::string s(argv[3]);
      const size_t      c = s.find(':');
      i0 = std::atoi(s.c_str());
      i1 = (c == std::string::npos) ? i0 : std::atoi(s.c_str() + c + 1);
    }
    occ::handle<XSControl_TransferReader> tr    = aReader.WS()->TransferReader();
    occ::handle<Interface_InterfaceModel> model = aReader.Model();
    int                                   k     = 0;
    for (TopExp_Explorer ex(aShape, TopAbs_FACE); ex.More(); ex.Next(), ++k)
    {
      if (k < i0 || k > i1)
        continue;
      const TopoDS_Face& aF = TopoDS::Face(ex.Current());
      BRepAdaptor_Surface anAS(aF);
      TopLoc_Location     aLoc;
      const occ::handle<Geom_Surface>& aRaw = BRep_Tool::Surface(aF, aLoc);
      int                              aNs = -1, aUc = 0, aVc = 0;
      if (!aRaw.IsNull())
      {
        occ::handle<ShapeAnalysis_Surface> aSAS = new ShapeAnalysis_Surface(aRaw);
        aNs = (int)aSAS->NbSingularities(1.e-7);
        aUc = aSAS->IsUClosed(1.e-7) ? 1 : 0;
        aVc = aSAS->IsVClosed(1.e-7) ? 1 : 0;
      }
      int                            id = 0;
      occ::handle<Standard_Transient> ent = tr->EntityFromShapeResult(aF, 1);
      if (!ent.IsNull())
        id = model->Number(ent);
      Bnd_Box aBox;
      BRepBndLib::Add(aF, aBox);
      double v[6] = {0, 0, 0, 0, 0, 0};
      if (!aBox.IsVoid())
        aBox.Get(v[0], v[1], v[2], v[3], v[4], v[5]);
      std::cout << std::fixed << std::setprecision(6) << "FID face=" << k << " id=" << id
                << " surf=" << (int)anAS.GetType()
                << " raw=" << (aRaw.IsNull() ? "null" : aRaw->DynamicType()->Name())
                << " nsing=" << aNs << " uclosed=" << aUc << " vclosed=" << aVc << " bbox=(" << v[0]
                << "," << v[1] << "," << v[2] << ")-(" << v[3] << "," << v[4] << "," << v[5] << ")\n";
      int wi = 0;
      for (TopExp_Explorer we(aF, TopAbs_WIRE); we.More(); we.Next(), ++wi)
      {
        int ne = 0;
        for (TopExp_Explorer ee(we.Current(), TopAbs_EDGE); ee.More(); ee.Next())
          ++ne;
        std::cout << "FID  wire[" << wi << "] nEdges=" << ne << "\n";
        int ei = 0;
        for (TopExp_Explorer ee(we.Current(), TopAbs_EDGE); ee.More(); ee.Next(), ++ei)
        {
          const TopoDS_Edge E = TopoDS::Edge(ee.Current());
          TopoDS_Vertex     vf, vl;
          TopExp::Vertices(E, vf, vl);
          std::cout << std::setprecision(6) << "FID   e[" << ei << "] ori=" << (int)E.Orientation()
                    << " deg=" << (BRep_Tool::Degenerated(E) ? 1 : 0);
          if (!vf.IsNull())
          {
            const gp_Pnt pf = BRep_Tool::Pnt(vf), pl = BRep_Tool::Pnt(vl);
            std::cout << " first=(" << pf.X() << "," << pf.Y() << "," << pf.Z() << ") last=("
                      << pl.X() << "," << pl.Y() << "," << pl.Z() << ")";
          }
          else
          {
            std::cout << " first=null last=null";
          }
          double                        f = 0, l = 0;
          const occ::handle<Geom_Curve>& aC3 = BRep_Tool::Curve(E, f, l);
          std::cout << " c3=" << (aC3.IsNull() ? "null" : aC3->DynamicType()->Name())
                    << " r3=" << f << ".." << l;
          double                            f2 = 0, l2 = 0;
          const occ::handle<Geom2d_Curve>& aC2 = BRep_Tool::CurveOnSurface(E, aF, f2, l2);
          std::cout << " c2=" << (aC2.IsNull() ? "null" : aC2->DynamicType()->Name())
                    << " r2=" << f2 << ".." << l2;
          if (!aC2.IsNull())
          {
            // TEMP T-104: pcurve own range + control-net hull (live and segmented).
            std::cout << " c2own=" << std::setprecision(9) << aC2->FirstParameter() << ".."
                      << aC2->LastParameter();
            occ::handle<Geom2d_BSplineCurve> aBS = occ::down_cast<Geom2d_BSplineCurve>(aC2);
            if (!aBS.IsNull())
            {
              Bnd_Box2d pb;
              for (int q = 1; q <= aBS->NbPoles(); ++q)
              {
                pb.Add(aBS->Pole(q));
              }
              double x1, y1, x2, y2;
              pb.Get(x1, y1, x2, y2);
              std::cout << " nbpoles=" << aBS->NbPoles() << " polebox=(" << x1 << "," << y1 << ")-("
                        << x2 << "," << y2 << ")";
              occ::handle<Geom2d_BSplineCurve> aCopy = occ::down_cast<Geom2d_BSplineCurve>(aBS->Copy());
              try
              {
                aCopy->Segment(f2, l2);
                Bnd_Box2d sb;
                for (int q = 1; q <= aCopy->NbPoles(); ++q)
                {
                  sb.Add(aCopy->Pole(q));
                }
                double a1, b1, a2, b2;
                sb.Get(a1, b1, a2, b2);
                std::cout << " nbpoles2=" << aCopy->NbPoles() << " segbox=(" << a1 << "," << b1 << ")-("
                          << a2 << "," << b2 << ")";
              }
              catch (const Standard_Failure&)
              {
                std::cout << " segthrow";
              }
            }
          }
          std::cout << "\n";
        }
      }
    }
    return 0;
  }
  return std::nullopt;
}
