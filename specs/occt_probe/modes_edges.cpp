// `argv[2]` modes: `--edev`, `--etol`, `--fbox`, `--wires`, `--edges`, `--topo`.
#include "occt_probe_common.h"

std::optional<int> runEdgeFaceModes(const ProbeCtx& theCtx)
{
  const int           argc = theCtx.argc;
  char** const        argv = theCtx.argv;
  const TopoDS_Shape& aShape = theCtx.shape;

  if (argc > 3 && std::string(argv[2]) == "--edev")
  {
    // TEMP: per-edge `ShapeAnalysis_Edge::CheckSameParameter` deviation on every
    // ancestor face of the edges of one face (`ShapeFix_Edge.cxx:869-928`).
    const int want = std::atoi(argv[3]);
    TopTools_IndexedDataMapOfShapeListOfShape aEF;
    TopExp::MapShapesAndAncestors(aShape, TopAbs_EDGE, TopAbs_FACE, aEF);
    std::vector<TopoDS_Face> aFaces;
    for (TopExp_Explorer ex(aShape, TopAbs_FACE); ex.More(); ex.Next())
    {
      aFaces.push_back(TopoDS::Face(ex.Current()));
    }
    int k = 0;
    for (TopExp_Explorer ex(aShape, TopAbs_FACE); ex.More(); ex.Next(), ++k)
    {
      if (k != want)
      {
        continue;
      }
      int ei = 0;
      for (TopExp_Explorer ee(ex.Current(), TopAbs_EDGE); ee.More(); ee.Next(), ++ei)
      {
        const TopoDS_Edge& anE = TopoDS::Edge(ee.Current());
        std::cout << std::fixed << std::setprecision(9) << "EDEV face=" << k << " e[" << ei
                  << "] etol=" << BRep_Tool::Tolerance(anE) << "\n";
        if (!aEF.Contains(anE))
        {
          std::cout << "  (no ancestor face)\n";
          continue;
        }
        for (TopTools_ListIteratorOfListOfShape it(aEF.FindFromKey(anE)); it.More(); it.Next())
        {
          const TopoDS_Face& aF = TopoDS::Face(it.Value());
          double             aDev = 0.0;
          ShapeAnalysis_Edge sae;
          const bool         st = sae.CheckSameParameter(anE, aF, aDev);
          int                fi = -1;
          for (int j = 0; j < (int)aFaces.size(); ++j)
          {
            if (aFaces[j].IsSame(aF))
            {
              fi = j;
              break;
            }
          }
          BRepAdaptor_Surface aAS(aF);
          double                cf = 0.0, cl = 0.0, f3 = 0.0, l3 = 0.0;
          BRep_Tool::Curve(anE, f3, l3);
          occ::handle<Geom2d_Curve> aPC = BRep_Tool::CurveOnSurface(anE, aF, cf, cl);
          std::cout << "  -> face=" << fi << " type=" << (int)aAS.GetType()
                    << " dev=" << std::setprecision(9) << aDev
                    << " ok=" << st << " r3=[" << f3 << "," << l3 << "]"
                    << " pc=" << (aPC.IsNull() ? "null" : aPC->DynamicType()->Name())
                    << " rp=[" << cf << "," << cl << "]"
                    << " etol=" << BRep_Tool::Tolerance(anE) << "\n";
        }
      }
      return 0;
    }
    std::cout << "EDEV none\n";
    return 0;
  }

  if (argc > 3 && std::string(argv[2]) == "--etol")
  {
    // TEMP: per-edge tolerance and the `ShapeFix::SameParameter` planar
    // deviation (`ShapeFix.cxx:224-243`) for one face.
    const int want = std::atoi(argv[3]);
    int       k    = 0;
    for (TopExp_Explorer ex(aShape, TopAbs_FACE); ex.More(); ex.Next(), ++k)
    {
      if (k != want)
      {
        continue;
      }
      const TopoDS_Face&          aF = TopoDS::Face(ex.Current());
      BRepAdaptor_Surface         aAS(aF);
      const occ::handle<Geom_Surface>& aSurf = BRep_Tool::Surface(aF);
      std::cout << std::fixed << std::setprecision(9) << "ETOL face=" << k
                << " ftol=" << BRep_Tool::Tolerance(aF)
                << " surf=" << (aSurf.IsNull() ? "null" : aSurf->DynamicType()->Name()) << "\n";
      int ei = 0;
      for (TopExp_Explorer ee(aF, TopAbs_EDGE); ee.More(); ee.Next(), ++ei)
      {
        TopoDS_Edge              E = TopoDS::Edge(ee.Current());
        double                   f3, l3;
        occ::handle<Geom_Curve>  c3d = BRep_Tool::Curve(E, f3, l3);
        double                   f = f3, l = l3;
        occ::handle<Geom2d_Curve> c2d = BRep_Tool::CurveOnSurface(E, aF, f, l);
        std::cout << std::setprecision(9) << "  e[" << ei << "] ori=" << (int)E.Orientation()
                  << " etol=" << BRep_Tool::Tolerance(E)
                  << " c3d=" << (c3d.IsNull() ? "null" : c3d->DynamicType()->Name())
                  << " r3=" << f3 << "," << l3;
        if (c2d.IsNull())
        {
          std::cout << " pc=null\n";
          continue;
        }
        std::cout << " pc=" << c2d->DynamicType()->Name() << " rp=" << f << "," << l;
        double maxd = 0.0;
        for (int i = 0; i < 23; i++)
        {
          double par = (f * (22 - i) + l * i) / 22.0;
          gp_Pnt p3  = c3d->Value(par);
          gp_Pnt2d p2 = c2d->Value(par);
          gp_Pnt pq  = aAS.Value(p2.X(), p2.Y());
          double d   = p3.Distance(pq);
          if (d > maxd)
          {
            maxd = d;
          }
        }
        std::cout << " sp_maxd=" << std::setprecision(9) << maxd
                  << " tol_formula=" << 1.00005 * maxd << "\n";
      }
      return 0;
    }
    std::cout << "ETOL none\n";
    return 0;
  }

  if (argc > 8 && std::string(argv[2]) == "--fbox")
  {
    const double q[6] = {atof(argv[3]), atof(argv[4]), atof(argv[5]),
                         atof(argv[6]), atof(argv[7]), atof(argv[8])};
    int k = 0;
    for (TopExp_Explorer ex(aShape, TopAbs_FACE); ex.More(); ex.Next(), ++k)
    {
      Bnd_Box b;
      BRepBndLib::Add(ex.Current(), b);
      if (b.IsVoid()) continue;
      double v[6];
      b.Get(v[0], v[1], v[2], v[3], v[4], v[5]);
      bool ok = true;
      for (int i = 0; i < 6; ++i)
        if (std::abs(v[i] - q[i]) > 0.5) ok = false;
      if (!ok) continue;
      const TopoDS_Face& aF = TopoDS::Face(ex.Current());
      BRepAdaptor_Surface anAS(aF);
      TopLoc_Location aLoc;
      const occ::handle<Geom_Surface>& aRaw = BRep_Tool::Surface(aF, aLoc);
      occ::handle<ShapeAnalysis_Surface> aSAS = new ShapeAnalysis_Surface(aRaw);
      std::cout << std::fixed << std::setprecision(5) << "FBOX face=" << k << " bbox (" << v[0] << "," << v[1]
                << "," << v[2] << ")-(" << v[3] << "," << v[4] << "," << v[5] << ")"
                << " ftol=" << std::setprecision(9) << BRep_Tool::Tolerance(aF)
                << std::setprecision(5)
                << " surf=" << (int)anAS.GetType()
                << " raw=" << (aRaw.IsNull() ? "null" : aRaw->DynamicType()->Name())
                << " astype=" << (int)aSAS->Adaptor3d()->GetType()
                << " nsing=" << aSAS->NbSingularities(1.e-7)
                << " uclosed=" << aSAS->IsUClosed(1.e-7)
                << " vclosed=" << aSAS->IsVClosed(1.e-7) << "\n";
      int wi = 0;
      for (TopExp_Explorer we(ex.Current(), TopAbs_WIRE); we.More(); we.Next(), ++wi)
      {
        int ne = 0;
        for (TopExp_Explorer ee(we.Current(), TopAbs_EDGE); ee.More(); ee.Next()) ++ne;
        std::cout << "FBOX  wire[" << wi << "] nEdges=" << ne << "\n";
        int ei = 0;
        for (TopExp_Explorer ee(we.Current(), TopAbs_EDGE); ee.More(); ee.Next(), ++ei)
        {
          TopoDS_Edge   E = TopoDS::Edge(ee.Current());
          TopoDS_Vertex vf, vl;
          TopExp::Vertices(E, vf, vl);
          gp_Pnt pf = BRep_Tool::Pnt(vf), pl = BRep_Tool::Pnt(vl);
          std::cout << std::fixed << std::setprecision(6) << "FBOX   e[" << ei
                    << "] ori=" << (int)E.Orientation() << " first=(" << pf.X() << "," << pf.Y() << ","
                    << pf.Z() << ") last=(" << pl.X() << "," << pl.Y() << "," << pl.Z()
                    << ") deg=" << BRep_Tool::Degenerated(E)
                    << " etol=" << std::setprecision(9) << BRep_Tool::Tolerance(E)
                    << " v1tol=" << BRep_Tool::Tolerance(vf)
                    << " v2tol=" << BRep_Tool::Tolerance(vl) << std::setprecision(6) << "\n";
        }
      }
      return 0;
    }
    std::cout << "FBOX none\n";
    return 0;
  }

  if (argc > 2 && std::string(argv[2]) == "--wires")
  {
    struct EW
    {
      bool          deg = false, closed = false, hasPC = false;
      TopoDS_Vertex v1, v2;
      double        u1 = 0, w1 = 0, u2 = 0, w2 = 0;
      int           idx = 0;
    };

    int                nFaces = 0, nWires = 0, nFaces2 = 0, nFacesDup = 0;
    int                nWires1Closed = 0, nWiresDupAny = 0, nFaces1Closed = 0, nFacesDupAny = 0;
    std::map<int, int> aHist;
    TopExp_Explorer    aFaceEx(aShape, TopAbs_FACE);
    for (; aFaceEx.More(); aFaceEx.Next())
    {
      const TopoDS_Face aFace = TopoDS::Face(aFaceEx.Current());
      ++nFaces;
      {
        int zzW = 0;
        for (TopExp_Explorer cwe(aFace, TopAbs_WIRE); cwe.More(); cwe.Next()) ++zzW;
        if (zzW >= 2)
        {
          Bnd_Box zb;
          BRepBndLib::Add(aFace, zb);
          if (!zb.IsVoid())
          {
            double x1, y1, z1, x2, y2, z2;
            zb.Get(x1, y1, z1, x2, y2, z2);
            std::cout << std::fixed << std::setprecision(3) << "MULTI bbox=(" << x1 << "," << y1 << ","
                      << z1 << ")-(" << x2 << "," << y2 << "," << z2 << ") wires=" << zzW << "\n";
          }
        }
      }
      int             aWireNb = 0;
      bool            isDup    = false;
      bool            isDupAny = false;
      std::string     aDetail;
      TopExp_Explorer aWireEx(aFace, TopAbs_WIRE);
      for (; aWireEx.More(); aWireEx.Next(), ++aWireNb)
      {
        ++nWires;
        std::vector<EW> aEdges;
        int             anEdgeIdx = 0;
        TopExp_Explorer anEdgeEx(aWireEx.Current(), TopAbs_EDGE);
        for (; anEdgeEx.More(); anEdgeEx.Next(), ++anEdgeIdx)
        {
          const TopoDS_Edge anEdge = TopoDS::Edge(anEdgeEx.Current());
          EW                anE;
          anE.idx = anEdgeIdx;
          anE.deg = BRep_Tool::Degenerated(anEdge);
          TopExp::Vertices(anEdge, anE.v1, anE.v2);
          anE.closed = (!anE.v1.IsNull() && !anE.v2.IsNull() && anE.v1.IsSame(anE.v2));
          double                    aF = 0, aL = 0;
          occ::handle<Geom2d_Curve> aC2d = BRep_Tool::CurveOnSurface(anEdge, aFace, aF, aL);
          if (!aC2d.IsNull())
          {
            const gp_Pnt2d aP1 = aC2d->Value(aF);
            const gp_Pnt2d aP2 = aC2d->Value(aL);
            anE.hasPC          = true;
            anE.u1 = aP1.X(); anE.w1 = aP1.Y();
            anE.u2 = aP2.X(); anE.w2 = aP2.Y();
          }
          aEdges.push_back(anE);
        }
        if (aEdges.size() == 1 && aEdges[0].closed)
        {
          ++nWires1Closed;
          ++nFaces1Closed;
          if (nWires1Closed <= 400)
          {
            const EW& e0 = aEdges[0];
            BRepAdaptor_Surface aSurf0(aFace);
            const double aD2 = std::sqrt((e0.u2 - e0.u1) * (e0.u2 - e0.u1)
                                         + (e0.w2 - e0.w1) * (e0.w2 - e0.w1));
            int                      aType = -1;
            double                   aF = 0, aL = 0, aCF = 0, aCL = 0;
            TCollection_AsciiString  aTN("none");
            {
              TopExp_Explorer ex0(aWireEx.Current(), TopAbs_EDGE);
              for (; ex0.More(); ex0.Next())
              {
                const TopoDS_Edge         e0x = TopoDS::Edge(ex0.Current());
                double                    f0 = 0, l0 = 0;
                occ::handle<Geom2d_Curve> c0 = BRep_Tool::CurveOnSurface(e0x, aFace, f0, l0);
                if (c0.IsNull())
                  continue;
                Geom2dAdaptor_Curve ad0(c0);
                aType = (int)ad0.GetType();
                aF    = f0; aL = l0;
                aCF   = c0->FirstParameter(); aCL = c0->LastParameter();
                aTN   = c0->DynamicType()->Name();
                break;
              }
            }
            std::cout << "C1 face=" << nFaces << " surf=" << (int)aSurf0.GetType()
                      << " d2d=" << aD2 << " hasPC=" << (e0.hasPC ? 1 : 0)
                      << " pctype=" << aType << " tn=" << aTN.ToCString()
                      << " range=" << aF << ".." << aL << " cfl=" << aCF << ".." << aCL
                      << " u=" << e0.u1 << ".." << e0.u2
                      << " v=" << e0.w1 << ".." << e0.w2 << "\n";
          }
        }
        for (size_t a = 0; a < aEdges.size(); ++a)
        {
          for (size_t b = a + 1; b < aEdges.size(); ++b)
          {
            const EW& ea = aEdges[a];
            const EW& eb = aEdges[b];
            if (!ea.hasPC || !eb.hasPC)
              continue;
            {
              const double aF2 = std::abs(ea.u1 - eb.u1) + std::abs(ea.w1 - eb.w1)
                                 + std::abs(ea.u2 - eb.u2) + std::abs(ea.w2 - eb.w2);
              const double aR2 = std::abs(ea.u1 - eb.u2) + std::abs(ea.w1 - eb.w2)
                                 + std::abs(ea.u2 - eb.u1) + std::abs(ea.w2 - eb.w1);
              if (std::min(aF2, aR2) < 1e-7)
                isDupAny = true;
            }
            if (isDup)
              continue;
            if (!ea.closed || !eb.closed)
              continue;
            if (!ea.v1.IsSame(eb.v1) || !ea.v2.IsSame(eb.v2))
              continue;
            const double aFwd = std::abs(ea.u1 - eb.u1) + std::abs(ea.w1 - eb.w1)
                                + std::abs(ea.u2 - eb.u2) + std::abs(ea.w2 - eb.w2);
            const double aRev = std::abs(ea.u1 - eb.u2) + std::abs(ea.w1 - eb.w2)
                                + std::abs(ea.u2 - eb.u1) + std::abs(ea.w2 - eb.w1);
            if (std::min(aFwd, aRev) < 1e-7)
              isDup = true;
          }
        }
        char aBuf[256];
        std::snprintf(aBuf, sizeof(aBuf), " w%d(e=%d)", aWireNb, (int)aEdges.size());
        aDetail += aBuf;
        for (const EW& e : aEdges)
        {
          std::snprintf(aBuf, sizeof(aBuf), " [cl=%d dg=%d u=%.4f..%.4f v=%.4f..%.4f]",
                        e.closed ? 1 : 0, e.deg ? 1 : 0, e.u1, e.u2, e.w1, e.w2);
          aDetail += aBuf;
        }
      }
      if (isDupAny)
      {
        ++nWiresDupAny;
        ++nFacesDupAny;
      }
      aHist[aWireNb]++;
      if (aWireNb >= 2)
        ++nFaces2;
      if (isDup)
      {
        ++nFacesDup;
        BRepAdaptor_Surface aSurf(aFace);
        std::cout << "DUP face=" << nFaces << " surf=" << (int)aSurf.GetType()
                  << " wires=" << aWireNb << aDetail << "\n";
      }
    }
    std::cout << "TOTAL faces=" << nFaces << " wires=" << nWires
              << " faces_with_2plus_wires=" << nFaces2
              << " faces_with_coincident_closed_edges=" << nFacesDup
              << " wires_1edge_closed=" << nWires1Closed << " faces_1edge_closed=" << nFaces1Closed
              << " wires_dup_pc_any=" << nWiresDupAny << " faces_dup_pc_any=" << nFacesDupAny << "\n";
    std::cout << "WIREHIST";
    for (const auto& aPair : aHist)
      std::cout << " " << aPair.first << ":" << aPair.second;
    std::cout << "\n";
    return 0;
  }

  // T-54 step-1 oracle: OCCT's **per-edge** discretization, i.e. how many
  // nodes `BRepMesh` put on each edge's polygon-on-triangulation. The port has
  // already ported the `compute_nb_samples*` family but never wired it, so this
  // is the table that family must reproduce.
  if (argc > 2 && std::string(argv[2]) == "--edges")
  {
    const double aDeflection = (argc > 3) ? std::atof(argv[3]) : 0.1;
    const double anAngle      = (argc > 4) ? std::atof(argv[4]) : 0.5;
    BRepMesh_IncrementalMesh aMesher(aShape, aDeflection, false, anAngle);
    TopExp_Explorer aFaceEx(aShape, TopAbs_FACE);
    int             aFaceIdx = 0, anEdgeIdx = 0, aTotal = 0, aNoPoly = 0;
    for (; aFaceEx.More(); aFaceEx.Next(), ++aFaceIdx)
    {
      TopLoc_Location                     aFaceLoc;
      const occ::handle<Poly_Triangulation>& aFaceTri =
        BRep_Tool::Triangulation(TopoDS::Face(aFaceEx.Current()), aFaceLoc);
      TopExp_Explorer anEdgeEx(aFaceEx.Current(), TopAbs_EDGE);
      for (; anEdgeEx.More(); anEdgeEx.Next(), ++anEdgeIdx)
      {
        const TopoDS_Edge& anEdge = TopoDS::Edge(anEdgeEx.Current());
        int                aN      = 0;
        if (!aFaceTri.IsNull())
        {
          const occ::handle<Poly_PolygonOnTriangulation>& aPoly =
            BRep_Tool::PolygonOnTriangulation(anEdge, aFaceTri, aFaceLoc);
          if (aPoly.IsNull())
          {
            ++aNoPoly;
          }
          else
          {
            aN = aPoly->NbNodes();
          }
        }
        aTotal += aN;
        std::cout << "EDGE face=" << aFaceIdx << " i=" << anEdgeIdx << " nodes=" << aN
                  << " degenerated=" << (BRep_Tool::Degenerated(anEdge) ? 1 : 0) << "\n";
      }
    }
    std::cout << "TOTAL faces=" << aFaceIdx << " edge_occurrences=" << anEdgeIdx
              << " edge_nodes=" << aTotal << " without_polygon=" << aNoPoly
              << " deflection=" << aDeflection << "\n";
    return 0;
  }

  // `Prs3d::GetDeflection` needs `maxComp(BRepBndLib::Add(shape, box, false))`; print
  // it so a GT run can use the same parameters the port's export path uses.
  if (argc > 2 && std::string(argv[2]) == "--topo")
  {
    int k = 0;
    for (TopExp_Explorer fx(aShape, TopAbs_FACE); fx.More(); fx.Next(), ++k)
    {
      const TopoDS_Face aF = TopoDS::Face(fx.Current());
      BRepAdaptor_Surface aS(aF);
      std::cout << "TOPO face " << k << " type=" << (int)aS.GetType()
                << " uper=" << (aS.IsUPeriodic() ? 1 : 0) << " vper=" << (aS.IsVPeriodic() ? 1 : 0)
                << " uw=[" << aS.FirstUParameter() << "," << aS.LastUParameter() << "]"
                << " vw=[" << aS.FirstVParameter() << "," << aS.LastVParameter() << "]";
      int aNW = 0;
      for (TopExp_Explorer wx(aF, TopAbs_WIRE); wx.More(); wx.Next(), ++aNW)
      {
        int anE = 0, aSeam = 0, anOriented = 0;
        for (TopExp_Explorer ex(wx.Current(), TopAbs_EDGE); ex.More(); ex.Next(), ++anE)
        {
          if (BRep_Tool::IsClosed(TopoDS::Edge(ex.Current()), aF))
            ++aSeam;
        }
        std::cout << " w" << aNW << "={e" << anE << ",seam" << aSeam << "}";
      }
      std::cout << "\n";
    }
    return 0;
  }
  return std::nullopt;
}
