// `argv[2]` modes: `--pcdef`, `--pcf`, `--spcheck` (p-curve definitions and
// re-projection).
#include "occt_probe_common.h"

std::optional<int> runPcurveModes(const ProbeCtx& theCtx)
{
  const int           argc = theCtx.argc;
  char** const        argv = theCtx.argv;
  const TopoDS_Shape& aShape = theCtx.shape;

  // T-191 oracle: full definition of a face's p-curve, so the port's stored STEP
  // 2-D curve can be diffed pole by pole against OCCT's object. Selection is by
  // the p-curve's own first parameter (`--pcdef <first> [<tol>]`) because the
  // face/edge indices of the port and of `TopExp_Explorer` are independent.
  {
    double aWantFirst = 0.;
    double aWantTol   = 1.e-7;
    bool   hasPcDef   = false;
    for (int i = 2; i < argc; ++i)
    {
      const std::string a = argv[i];
      if (a == "--pcdef" && i + 1 < argc)
      {
        aWantFirst = std::atof(argv[i + 1]);
        hasPcDef   = true;
        if (i + 2 < argc && argv[i + 2][0] != '-')
        {
          aWantTol = std::atof(argv[i + 2]);
        }
      }
    }
    if (hasPcDef)
    {
      int          fi = 0;
      int          nHit = 0;
      TopExp_Explorer fe(aShape, TopAbs_FACE);
      for (; fe.More(); fe.Next(), ++fi)
      {
        const TopoDS_Face& aFace = TopoDS::Face(fe.Current());
        int                wi = 0;
        for (TopExp_Explorer we(aFace, TopAbs_WIRE); we.More(); we.Next(), ++wi)
        {
          int ei = 0;
          for (TopExp_Explorer ee(we.Current(), TopAbs_EDGE); ee.More(); ee.Next(), ++ei)
          {
            const TopoDS_Edge anE = TopoDS::Edge(ee.Current());
            double            pf = 0., pl = 0.;
            occ::handle<Geom2d_Curve> aPC = BRep_Tool::CurveOnSurface(anE, aFace, pf, pl);
            if (aPC.IsNull() || std::abs(pf - aWantFirst) > aWantTol)
            {
              continue;
            }
            occ::handle<Geom2d_Curve> aBasis = aPC;
            if (aPC->IsKind(STANDARD_TYPE(Geom2d_TrimmedCurve)))
            {
              aBasis = occ::down_cast<Geom2d_TrimmedCurve>(aPC)->BasisCurve();
            }
            ++nHit;
            std::cout << std::setprecision(17) << "PCDEF face=" << fi << " wire=" << wi
                      << " edge=" << ei << " type=" << aBasis->DynamicType()->Name()
                      << " pf=" << pf << " pl=" << pl << "\n";
            if (std::getenv("SURFD2") != nullptr)
            {
              // TEMP: surface `D0` / `D2` at the same (u, v) pairs the port's
              // `SURFD2` dump uses, from the adaptor the mesher itself evaluates.
              const BRepAdaptor_Surface aS(aFace);
              const std::string         aSpec = std::getenv("SURFD2");
              size_t                    aPos  = 0;
              while (aPos <= aSpec.size())
              {
                const size_t aEnd = aSpec.find(';', aPos);
                const std::string aPair =
                  aSpec.substr(aPos, aEnd == std::string::npos ? std::string::npos : aEnd - aPos);
                double aU = 0., aV = 0.;
                if (std::sscanf(aPair.c_str(), "%lf,%lf", &aU, &aV) == 2)
                {
                  gp_Pnt aP0, aPD;
                  gp_Vec aRU, aRV, aRUU, aRVV, aRUV;
                  aS.D0(aU, aV, aP0);
                  aS.D2(aU, aV, aPD, aRU, aRV, aRUU, aRVV, aRUV);
                  std::cout << std::setprecision(17) << "OCCTSRFD2 uv=(" << aU << "," << aV
                            << ") d0=(" << aP0.X() << "," << aP0.Y() << "," << aP0.Z()
                            << ") d2p=(" << aPD.X() << "," << aPD.Y() << "," << aPD.Z()
                            << ") ru=(" << aRU.X() << "," << aRU.Y() << "," << aRU.Z()
                            << ") rv=(" << aRV.X() << "," << aRV.Y() << "," << aRV.Z()
                            << ") ruu=(" << aRUU.X() << "," << aRUU.Y() << "," << aRUU.Z()
                            << ") rvv=(" << aRVV.X() << "," << aRVV.Y() << "," << aRVV.Z()
                            << ") ruv=(" << aRUV.X() << "," << aRUV.Y() << "," << aRUV.Z() << ")\n";
                }
                if (aEnd == std::string::npos)
                {
                  break;
                }
                aPos = aEnd + 1;
              }
            }
            if (std::getenv("OCCTUV") != nullptr)
            {
              // TEMP: `ShapeAnalysis_Surface::ValueOfUV` for the same 23 samples
              // the port's `PCUV2` dump projects, so the projectors can be
              // compared without the `NextValueOfUV` chain or `getLine` seed.
              const char* aPrecEnv = std::getenv("OCCTUV_PRECI");
              const double aPreci = (aPrecEnv != nullptr) ? std::atof(aPrecEnv) : 1.e-7;
              ShapeAnalysis_Surface aSAS(BRep_Tool::Surface(aFace));
              int aN = 23;
              double f3 = 0., l3 = 0.;
              const occ::handle<Geom_Curve> aC3 = BRep_Tool::Curve(anE, f3, l3);
              if (!aC3.IsNull())
              {
                for (int k = 0; k < aN; ++k)
                {
                  const double aT = pf + (pl - pf) * (double)k / (double)(aN - 1);
                  const gp_Pnt   aP3 = aC3->Value(aT);
                  const gp_Pnt2d aUV = aSAS.ValueOfUV(aP3, aPreci);
                  const double   aD  = aSAS.Value(aUV.X(), aUV.Y()).Distance(aP3);
                  std::cout << "OCCTUV i=" << k << " preci=" << aPreci << " uv=(" << aUV.X()
                            << "," << aUV.Y() << ") d=" << aD << "\n";
                }
              }
            }
            if (std::getenv("C3VAL") != nullptr)
            {
              // TEMP: the 3-D curve points `ShapeConstruct_ProjectCurveOnSurface`
              // samples before projecting them, matching the port's `PC3D` dump.
              double f3 = 0., l3 = 0.;
              const occ::handle<Geom_Curve> aC3 = BRep_Tool::Curve(anE, f3, l3);
              if (!aC3.IsNull())
              {
                int aN = 23;
                for (int k = 0; k < aN; ++k)
                {
                  const double aT = pf + (pl - pf) * (double)k / (double)(aN - 1);
                  const gp_Pnt aP3 = aC3->Value(aT);
                  std::cout << "C3VAL i=" << k << " t=" << aT << " p=(" << aP3.X() << ","
                            << aP3.Y() << "," << aP3.Z() << ")\n";
                }
              }
            }
            if (aBasis->IsKind(STANDARD_TYPE(Geom2d_BSplineCurve)))
            {
              occ::handle<Geom2d_BSplineCurve> aBS = occ::down_cast<Geom2d_BSplineCurve>(aBasis);
              const bool aRat = aBS->IsRational();
              std::cout << "PCDEF  degree=" << aBS->Degree() << " npoles=" << aBS->NbPoles()
                        << " nknots=" << aBS->NbKnots() << " periodic=" << (aBS->IsPeriodic() ? 1 : 0)
                        << " rational=" << (aRat ? 1 : 0) << "\n";
              for (int i = 1; i <= aBS->NbPoles(); ++i)
              {
                const gp_Pnt2d aP = aBS->Pole(i);
                std::cout << "PCDEF  pole[" << i << "]=" << aP.X() << "," << aP.Y();
                if (aRat)
                {
                  std::cout << " w=" << aBS->Weight(i);
                }
                std::cout << "\n";
              }
              for (int i = 1; i <= aBS->NbKnots(); ++i)
              {
                std::cout << "PCDEF  knot[" << i << "]=" << aBS->Knot(i)
                          << " m=" << aBS->Multiplicity(i) << "\n";
              }
              std::cout << "PCDEF  knotval";
              for (int i = 1; i <= aBS->NbKnots(); ++i)
              {
                const gp_Pnt2d aV = aBS->Value(aBS->Knot(i));
                std::cout << " [" << i << "]" << aV.X() << "," << aV.Y();
              }
              std::cout << "\n";
              std::cout << "PCDEF  val";
              for (int k = 0; k <= 4; ++k)
              {
                const double aU = pf + (pl - pf) * (double)k / 4.0;
                const gp_Pnt2d aV = aBS->Value(aU);
                std::cout << " [" << k << "]" << aV.X() << "," << aV.Y();
              }
              {
                const gp_Pnt2d aV = aBS->Value(-0.19481982426506805);
                std::cout << " [t]" << aV.X() << "," << aV.Y();
              }
              std::cout << "\n";
              std::cout << "PCDEF  flat=";
              {
                const TColStd_Array1OfReal& aFlat = aBS->Knots();
                for (int i = aFlat.Lower(); i <= aFlat.Upper(); ++i)
                {
                  std::cout << aFlat.Value(i) << (i < aFlat.Upper() ? "," : "\n");
                }
              }
            }
          }
        }
      }
      std::cout << "PCDEF hits=" << nHit << "\n";
      return 0;
    }
  }

  // T-192 oracle: per-edge pcurve dump for one face index (`--pcf <face>`),
  // sorted into the same order as the port's `PCDEF f=<face>` dump so the two
  // can be diffed edge by edge.
  {
    int  aWantFace = -1;
    bool hasPcF    = false;
    for (int i = 2; i < argc; ++i)
    {
      const std::string a = argv[i];
      if (a == "--pcf" && i + 1 < argc)
      {
        aWantFace = std::atoi(argv[i + 1]);
        hasPcF    = true;
      }
    }
    if (hasPcF)
    {
      int fi = 0;
      for (TopExp_Explorer fe(aShape, TopAbs_FACE); fe.More(); fe.Next(), ++fi)
      {
        if (fi != aWantFace)
        {
          continue;
        }
        const TopoDS_Face& aFace = TopoDS::Face(fe.Current());
        int                wi = 0;
        for (TopExp_Explorer we(aFace, TopAbs_WIRE); we.More(); we.Next(), ++wi)
        {
          int ei = 0;
          for (TopExp_Explorer ee(we.Current(), TopAbs_EDGE); ee.More(); ee.Next(), ++ei)
          {
            const TopoDS_Edge anE = TopoDS::Edge(ee.Current());
            double            pf = 0., pl = 0., e3f = 0., e3l = 0.;
            occ::handle<Geom2d_Curve> aPC = BRep_Tool::CurveOnSurface(anE, aFace, pf, pl);
            occ::handle<Geom_Curve>   aC3 = BRep_Tool::Curve(anE, e3f, e3l);
            occ::handle<Geom2d_Curve> aBasis = aPC;
            if (!aPC.IsNull() && aPC->IsKind(STANDARD_TYPE(Geom2d_TrimmedCurve)))
            {
              aBasis = occ::down_cast<Geom2d_TrimmedCurve>(aPC)->BasisCurve();
            }
            int np = -1, nk = -1, dg = -1;
            if (!aBasis.IsNull() && aBasis->IsKind(STANDARD_TYPE(Geom2d_BSplineCurve)))
            {
              const occ::handle<Geom2d_BSplineCurve> aBS =
                occ::down_cast<Geom2d_BSplineCurve>(aBasis);
              np = aBS->NbPoles();
              nk = aBS->NbKnots();
              dg = aBS->Degree();
            }
            std::cout << std::setprecision(17) << "PCF face=" << fi << " wire=" << wi
                      << " edge=" << ei << " c3f=" << e3f << " c3l=" << e3l
                      << " pf=" << pf << " pl=" << pl
                      << " pctype=" << (aBasis.IsNull() ? "none" : aBasis->DynamicType()->Name())
                      << " np=" << np << " nk=" << nk << " deg=" << dg;
            if (!aPC.IsNull())
            {
              for (int k = 0; k <= 4; ++k)
              {
                const double   aU = pf + (pl - pf) * (double)k / 4.0;
                const gp_Pnt2d aV = aPC->Value(aU);
                std::cout << " [" << k << "]" << aV.X() << "," << aV.Y();
              }
            }
            std::cout << "\n";
          }
        }
      }
      return 0;
    }
  }

  // T-191 oracle #2: replay `XSAlgo_ShapeProcessor::CheckPCurve`
  // (`XSAlgo_ShapeProcessor.cxx:344-505`) on the imported shape. The tolerance
  // of the temporary edge after `ShapeFix_Edge::FixSameParameter` (`cxx:444`)
  // decides whether OCCT re-projects the pcurve at all (`cxx:449`), so it is
  // the number to compare with the port's
  // `shhealing::pcurve_ranges::check_pcurve` (`a_tol`). Selection is by the
  // pcurve's own first parameter, like `--pcdef`.
  {
    double aWantFirst = 0.;
    double aWantTol   = 1.e-7;
    double aPreci     = Interface_Static::RVal("read.precision.val");
    bool   hasSpChk   = false;
    for (int i = 2; i < argc; ++i)
    {
      const std::string a = argv[i];
      if (a == "--spcheck" && i + 1 < argc)
      {
        aWantFirst = std::atof(argv[i + 1]);
        hasSpChk   = true;
        if (i + 2 < argc && argv[i + 2][0] != '-')
        {
          aWantTol = std::atof(argv[i + 2]);
          if (i + 3 < argc && argv[i + 3][0] != '-')
          {
            aPreci = std::atof(argv[i + 3]);
          }
        }
      }
    }
    if (hasSpChk)
    {
      // `XSAlgo_ShapeProcessor.cxx:317-340` `MakeEdgeOnCurve`.
      auto makeEdgeOnCurve = [](const TopoDS_Edge& theEdge) -> TopoDS_Edge {
        ShapeAnalysis_Edge anAnalyzer;
        occ::handle<Geom_Curve> aCurve;
        double f = 0., l = 0.;
        if (!anAnalyzer.Curve3d(theEdge, aCurve, f, l, false) || aCurve.IsNull())
        {
          return TopoDS_Edge();
        }
        BRepBuilderAPI_MakeEdge aMaker(aCurve, aCurve->Value(f), aCurve->Value(l), f, l);
        ShapeBuild_Edge         aSBE;
        aSBE.SetRange3d(aMaker, f, l);
        return aMaker.Edge();
      };
      auto poleCount = [](const occ::handle<Geom2d_Curve>& thePC) -> int {
        occ::handle<Geom2d_Curve> aBasis = thePC;
        if (thePC->IsKind(STANDARD_TYPE(Geom2d_TrimmedCurve)))
        {
          aBasis = occ::down_cast<Geom2d_TrimmedCurve>(thePC)->BasisCurve();
        }
        if (aBasis->IsKind(STANDARD_TYPE(Geom2d_BSplineCurve)))
        {
          return occ::down_cast<Geom2d_BSplineCurve>(aBasis)->NbPoles();
        }
        return -1;
      };

      int fi = 0, nHit = 0;
      std::cout << std::setprecision(17) << "SPCHK preci=" << aPreci << "\n";
      for (TopExp_Explorer fe(aShape, TopAbs_FACE); fe.More(); fe.Next(), ++fi)
      {
        const TopoDS_Face& aFace = TopoDS::Face(fe.Current());
        int                wi = 0;
        for (TopExp_Explorer we(aFace, TopAbs_WIRE); we.More(); we.Next(), ++wi)
        {
          int ei = 0;
          for (TopExp_Explorer ee(we.Current(), TopAbs_EDGE); ee.More(); ee.Next(), ++ei)
          {
            const TopoDS_Edge anE = TopoDS::Edge(ee.Current());
            ShapeAnalysis_Edge anAnalyzer;
            occ::handle<Geom2d_Curve> aPC;
            double pf = 0., pl = 0.;
            if (!anAnalyzer.PCurve(anE, aFace, aPC, pf, pl, false)
                || std::abs(pf - aWantFirst) > aWantTol)
            {
              continue;
            }
            ++nHit;
            const int npSrc = poleCount(aPC);
            // `cxx:411-446`: temporary edge + pcurve, SameRange(false),
            // FixSameParameter, then read the flags back.
            TopoDS_Edge aTmp = makeEdgeOnCurve(anE);
            if (aTmp.IsNull())
            {
              std::cout << "SPCHK face=" << fi << " wire=" << wi << " edge=" << ei
                        << " no 3d curve\n";
              continue;
            }
            BRep_Builder B;
            B.UpdateEdge(aTmp, aPC, aFace, 0.);
            B.Range(aTmp, aFace, pf, pl);
            B.SameRange(aTmp, false);
            if (Interface_Static::IVal("read.stdsameparameter.mode"))
            {
              B.SameParameter(aTmp, false);
            }
            occ::handle<ShapeFix_Edge> aFix = new ShapeFix_Edge;
            aFix->FixSameParameter(aTmp);
            const double aTol1 = BRep_Tool::Tolerance(aTmp);
            const bool   aSR1  = BRep_Tool::SameRange(aTmp);
            const bool   aSP1  = BRep_Tool::SameParameter(aTmp);
            std::cout << "SPCHK face=" << fi << " wire=" << wi << " edge=" << ei
                      << " pf=" << pf << " pl=" << pl << " npoles_src=" << npSrc
                      << " tol=" << aTol1 << " same_range=" << (aSR1 ? 1 : 0)
                      << " same_param=" << (aSP1 ? 1 : 0);
            const double aLimit = std::min(1., 2. * aPreci);
            if (aTol1 > aLimit || !aSR1)
            {
              // `cxx:451-465`: fresh temporary edge, FixAddPCurve (projection)
              // then FixSameParameter; take it only when it is better.
              TopoDS_Edge aPr = makeEdgeOnCurve(anE);
              aFix->FixAddPCurve(aPr, aFace, false, aPreci);
              occ::handle<Geom2d_Curve> aPC1b;
              double                    f1b = 0., l1b = 0.;
              if (!anAnalyzer.PCurve(aPr, aFace, aPC1b, f1b, l1b, false))
              {
                std::cout << " reproj_pc=none";
              }
              else
              {
                // T-191: which projection branch produced the pcurve.
                occ::handle<Geom2d_Curve> aSrcBasis = aPC;
                if (aSrcBasis->IsKind(STANDARD_TYPE(Geom2d_TrimmedCurve)))
                {
                  aSrcBasis = occ::down_cast<Geom2d_TrimmedCurve>(aSrcBasis)->BasisCurve();
                }
                int    aSrcNk = -1;
                int    aNewNk = -1;
                int    aNewDeg = -1;
                if (aSrcBasis->IsKind(STANDARD_TYPE(Geom2d_BSplineCurve)))
                {
                  aSrcNk = occ::down_cast<Geom2d_BSplineCurve>(aSrcBasis)->NbKnots();
                }
                occ::handle<Geom2d_Curve> aNewBasis = aPC1b;
                if (aNewBasis->IsKind(STANDARD_TYPE(Geom2d_TrimmedCurve)))
                {
                  aNewBasis = occ::down_cast<Geom2d_TrimmedCurve>(aNewBasis)->BasisCurve();
                }
                if (aNewBasis->IsKind(STANDARD_TYPE(Geom2d_BSplineCurve)))
                {
                  aNewNk  = occ::down_cast<Geom2d_BSplineCurve>(aNewBasis)->NbKnots();
                  aNewDeg = occ::down_cast<Geom2d_BSplineCurve>(aNewBasis)->Degree();
                }
                double aMaxUv = 0.;
                for (int k = 0; k <= 10; ++k)
                {
                  const double aT = k / 10.;
                  const gp_Pnt2d aA =
                    aPC->Value(pf + (pl - pf) * aT);
                  const gp_Pnt2d aB =
                    aPC1b->Value(f1b + (l1b - f1b) * aT);
                  aMaxUv = std::max(aMaxUv, aA.Distance(aB));
                }
                occ::handle<Geom_Surface> aSurf = BRep_Tool::Surface(aFace);
                int                        aUsedKnots = -1;
                int                        aMinPnt    = -1;
                occ::handle<Geom_Curve>    aC3d;
                double                     a3f = 0., a3l = 0.;
                anAnalyzer.Curve3d(anE, aC3d, a3f, a3l, false);
                if (!aC3d.IsNull() && aC3d->IsKind(STANDARD_TYPE(Geom_BSplineCurve)))
                {
                  const occ::handle<Geom_BSplineCurve> a3bs =
                    occ::down_cast<Geom_BSplineCurve>(aC3d);
                  aUsedKnots = 0;
                  for (int ki = 1; ki < a3bs->NbKnots(); ++ki)
                  {
                    if (a3bs->Knot(ki + 1) > a3f && a3bs->Knot(ki) < a3l)
                    {
                      aUsedKnots++;
                    }
                  }
                  aMinPnt = aUsedKnots * (((const Geom_BSplineCurve&)*a3bs).Degree() + 1);
                }
                std::cout << " reproj_type=" << aPC1b->DynamicType()->Name()
                          << " reproj_nknots=" << aNewNk << " reproj_deg=" << aNewDeg
                          << " src_npoles=" << npSrc << " src_nknots=" << aSrcNk
                          << " uv_dev=" << aMaxUv
                          << " surf=" << aSurf->DynamicType()->Name()
                          << " used_knots=" << aUsedKnots << " min_pnt=" << aMinPnt;
              }

              aFix->FixSameParameter(aPr);
              const double aTol2 = BRep_Tool::Tolerance(aPr);
              occ::handle<Geom2d_Curve> aPC2;
              double                    f2 = 0., l2 = 0.;
              anAnalyzer.PCurve(aPr, aFace, aPC2, f2, l2, false);
              std::cout << " reproj_tol=" << aTol2 << " reproj_npoles="
                        << (aPC2.IsNull() ? -1 : poleCount(aPC2))
                        << " pre_fix_npoles=" << (aPC1b.IsNull() ? -1 : poleCount(aPC1b))
                        << " take=" << ((aTol2 < aTol1 || !aSR1) ? 1 : 0);
            }
            else
            {
              std::cout << " reproj=skipped";
            }
            std::cout << "\n";
          }
        }
      }
      std::cout << "SPCHK hits=" << nHit << "\n";
      return 0;
    }
  }
  return std::nullopt;
}
