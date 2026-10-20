// IGES writer oracle: `--iges` / `--igespc`.
#include "occt_probe_common.h"

// --- T-118 IGES oracle -------------------------------------------------------------
// --iges <out.iges> writes the transferred shape with OCCT's own
// IGESControl_Writer in the settings the Rust writer ports: IGESControl_Controller
// initialised, Faces mode (write.iges.brep.mode = 0, IGESData.cxx:90-94) and
// millimetres. The result is the reference GeomToIGES_* / BRepToIGES_BRShell
// output, so the port's entity set and every parameter block can be compared
// against it instead of against a previous Rust run.
int runIges(const TopoDS_Shape& theShape, const char* thePath)
{
  IGESControl_Controller::Init();
  IGESControl_Writer aWriter("MM", 0);
  aWriter.AddShape(theShape);
  aWriter.ComputeModel();
  if (!aWriter.Write(thePath))
  {
    std::cerr << "IGES write failed: " << thePath << "\n";
    return 5;
  }
  std::cout << "IGES written: " << thePath << "\n";
  return 0;
}

// --- end T-118 ---------------------------------------------------------------------

// --- T-111 oracle: raw face/edge pcurve dump ---------------------------------------
// `--igespc <out.txt>` lists, in TopExp face order, the surface type and U/V bounds
// of every face together with the *raw* edge pcurves (`BRep_Tool::CurveOnSurface`)
// exactly as `BRepToIGES_BRShell::TransferFace` sees them. It answers whether a
// reversed spherical face carries a Geom2d_Line or a projected BSpline pcurve in
// the shape the IGES writer is handed, without guessing from the writer output.
int runIgesPc(const TopoDS_Shape& theShape, const char* thePath)
{
  std::ofstream aOut(thePath);
  int aFaceIdx = 0;
  for (TopExp_Explorer aFx(theShape, TopAbs_FACE); aFx.More(); aFx.Next(), ++aFaceIdx)
  {
    const TopoDS_Face aFace = TopoDS::Face(aFx.Current());
    TopLoc_Location aLoc;
    const Handle(Geom_Surface) aSurf = BRep_Tool::Surface(aFace, aLoc);
    double aU1 = 0, aU2 = 0, aV1 = 0, aV2 = 0;
    aSurf->Bounds(aU1, aU2, aV1, aV2);
    aOut << "FACE[" << aFaceIdx << "] "
         << (aFace.Orientation() == TopAbs_REVERSED ? "REV" : "FWD")
         << " surf=" << aSurf->DynamicType()->Name()
         << " u=[" << aU1 << "," << aU2 << "] v=[" << aV1 << "," << aV2 << "]"
         << " uclosed=" << (aSurf->IsUClosed() ? 1 : 0)
         << " vclosed=" << (aSurf->IsVClosed() ? 1 : 0) << "\n";
    int aEdgeIdx = 0;
    for (TopExp_Explorer aEx(aFace, TopAbs_EDGE); aEx.More(); aEx.Next(), ++aEdgeIdx)
    {
      const TopoDS_Edge aEdge = TopoDS::Edge(aEx.Current());
      // The pcurve range and the 3D range are read into SEPARATE variables:
      // `BRep_Tool::Curve` writes its own range into its two by-reference
      // arguments (`BRep_Tool.cxx:172-196`), so sharing them would print the 3D
      // range in `pcrange` and evaluate `aPc` at the 3D parameters. On an edge
      // whose 3D representation holds a null curve (`BRep_Builder::Degenerated`,
      // `BRep_Builder.cxx:1073-1085`) that range is the `BRep_GCurve` default
      // `[-RealLast, RealLast]` and the pcurve numbers become meaningless.
      double                aPCF = 0, aPCL = 0;
      Handle(Geom2d_Curve)  aPc  = BRep_Tool::CurveOnSurface(aEdge, aFace, aPCF, aPCL);
      TopLoc_Location       aELoc;
      double                aF = 0, aL = 0;
      const Handle(Geom_Curve) aC3 = BRep_Tool::Curve(aEdge, aELoc, aF, aL);
      aOut << "  EDGE[" << aEdgeIdx << "] "
           << (aEdge.Orientation() == TopAbs_REVERSED ? "REV" : "FWD")
           << " c3=" << (aC3.IsNull() ? "null" : aC3->DynamicType()->Name())
           << " pc=" << (aPc.IsNull() ? "null" : aPc->DynamicType()->Name())
           << " pcrange=[" << aPCF << "," << aPCL << "]"
           << " null3d=" << (aC3.IsNull() ? 1 : 0);
      {
        // The edge's own range, the value `BRep_Tool::Range` resolves over the
        // representation list (`BRep_Tool.cxx:931-960`): it skips a `BRep_Curve3D`
        // whose curve is null and falls through to the curve-on-surface one.
        double aEF = 0, aEL = 0;
        BRep_Tool::Range(aEdge, aEF, aEL);
        aOut << " erange=[" << aEF << "," << aEL << "]";
      }
      if (!aPc.IsNull())
      {
        gp_Pnt2d aP1, aP2;
        aPc->D0(aPCF, aP1);
        aPc->D0(aPCL, aP2);
        aOut << " p1=(" << aP1.X() << "," << aP1.Y() << ")"
             << " p2=(" << aP2.X() << "," << aP2.Y() << ")";
        if (aPc->DynamicType() == STANDARD_TYPE(Geom2d_Line))
        {
          const gp_Dir2d aD = Handle(Geom2d_Line)::DownCast(aPc)->Lin2d().Direction();
          aOut << " dir=(" << aD.X() << "," << aD.Y() << ")";
        }
        if (aPc->DynamicType() == STANDARD_TYPE(Geom2d_BSplineCurve))
        {
          const Handle(Geom2d_BSplineCurve) aBs = Handle(Geom2d_BSplineCurve)::DownCast(aPc);
          aOut << " deg=" << aBs->Degree() << " npol=" << aBs->NbPoles()
               << " nknot=" << aBs->NbKnots() << " rational=" << (aBs->IsRational() ? 1 : 0)
               << " periodic=" << (aBs->IsPeriodic() ? 1 : 0) << " knots=[";
          for (int k = 1; k <= aBs->NbKnots(); ++k)
          {
            aOut << (k > 1 ? "," : "") << aBs->Knot(k);
          }
          aOut << "]";
        }
        // T-111: 17-digit form of the same numbers, plus the pcurve's own
        // parameter range, so the port's `getLine` inputs (the 3D curve range
        // `theParams(1)` / `theParams(n)`, the two endpoint projections and the
        // `|aLLength - dPar| <= PConfusion` decision) can be compared exactly.
        {
          aOut << "\n    PREC";
          aOut << std::setprecision(17);
          aOut << " c3f=" << aF << " c3l=" << aL;
          double aPF = 0, aPL = 0;
          Handle(Geom2d_Curve) aPC2 = BRep_Tool::CurveOnSurface(aEdge, aFace, aPF, aPL);
          aOut << " pcr=[" << aPF << "," << aPL << "]";
          // `q1` / `q2` / `dPar` / `aLLength` sample the pcurve at the 3D curve's
          // own range, which is the parameter set `getLine` uses. Without a 3D
          // curve that range is the null-representation default `[-RealLast,
          // RealLast]`, so these four numbers are only printed when `aC3` exists.
          if (!aC3.IsNull())
          {
            gp_Pnt2d aQ1, aQ2;
            aPc->D0(aF, aQ1);
            aPc->D0(aL, aQ2);
            aOut << " q1=(" << aQ1.X() << "," << aQ1.Y() << ")"
                 << " q2=(" << aQ2.X() << "," << aQ2.Y() << ")";
            aOut << " dPar=" << (aL - aF)
                 << " aLLength=" << aQ1.Distance(aQ2);
          }
          aOut << std::setprecision(6);
        }
        // T-111: replay `ShapeFix_Edge::FixAddPCurve`'s projection call on this
        // very edge/face, with a fresh projector (empty `myCache`), and report
        // what `ShapeConstruct_ProjectCurveOnSurface::Perform` returns. This
        // separates "the port's `getLine` verdict differs" from "the stored
        // pcurve came from somewhere else in the healing chain".
        if (!aC3.IsNull())
        {
          Handle(ShapeAnalysis_Surface) aSAS = new ShapeAnalysis_Surface(aSurf);
          Handle(ShapeConstruct_ProjectCurveOnSurface) aProj =
            new ShapeConstruct_ProjectCurveOnSurface();
          const double aPrec = BRep_Tool::Tolerance(aEdge);
          aProj->Init(aSAS, aPrec);
          double aTolFirst = -1, aTolLast = -1;
          TopoDS_Vertex aV1, aV2;
          TopExp::Vertices(aEdge, aV1, aV2);
          if (!aV1.IsNull())
          {
            aTolFirst = BRep_Tool::Tolerance(aV1);
          }
          if (!aV2.IsNull())
          {
            aTolLast = BRep_Tool::Tolerance(aV2);
          }
          Handle(Geom2d_Curve) aRep;
          aProj->Perform(aC3, aF, aL, aRep, aTolFirst, aTolLast);
          aOut << "\n    REPROJ prec=" << aPrec << " tf=" << aTolFirst << " tl=" << aTolLast
               << " type=" << (aRep.IsNull() ? "null" : aRep->DynamicType()->Name());
          if (!aRep.IsNull())
          {
            aOut << std::setprecision(17);
            if (aRep->DynamicType() == STANDARD_TYPE(Geom2d_BSplineCurve))
            {
              const Handle(Geom2d_BSplineCurve) aRB = Handle(Geom2d_BSplineCurve)::DownCast(aRep);
              aOut << " rdeg=" << aRB->Degree() << " rnpol=" << aRB->NbPoles()
                   << " rnknot=" << aRB->NbKnots();
              gp_Pnt2d aR1, aR2;
              aRB->D0(aRB->FirstParameter(), aR1);
              aRB->D0(aRB->LastParameter(), aR2);
              aOut << " rrng=[" << aRB->FirstParameter() << "," << aRB->LastParameter() << "]"
                   << " rq1=(" << aR1.X() << "," << aR1.Y() << ") rq2=(" << aR2.X() << "," << aR2.Y()
                   << ")";
            }
            else
            {
              gp_Pnt2d aR1, aR2;
              aRep->D0(aF, aR1);
              aRep->D0(aL, aR2);
              aOut << " rrng=[" << aRep->FirstParameter() << "," << aRep->LastParameter() << "]"
                   << " rq1=(" << aR1.X() << "," << aR1.Y() << ") rq2=(" << aR2.X() << "," << aR2.Y()
                   << ")";
            }
            aOut << std::setprecision(6);
          }
        }
      }
      aOut << "\n";
    }
  }
  aOut.close();
  std::cout << "igespc written: " << thePath << "\n";
  return 0;
}

// --- end T-111 oracle --------------------------------------------------------------
