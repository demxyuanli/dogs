// Offset and surface/surface extrema probes: `--offsetsurf`, `--ssx`.
#include "occt_probe_common.h"

// T-126: `--offsetsurf <dist>` -- for every face, the `BRepOffset::Surface`
// result that the port's `brep_offset::brep_offset_surface`
// (`crates/occt-topo/src/brep_offset/brep_offset_surface.rs`) must reproduce:
// the basis surface's dynamic type, the result's dynamic type, the
// `BRepOffset_Status`, and the result evaluated on a fixed `(u, v)` grid. The
// grid is deliberately independent of the face's UV window, so this compares
// the *surface* and not the port's `BRepTools::UVBounds`. The face's trimmed
// window is unwrapped first, the way `BRepOffset_Offset::Init` does
// (`BRepOffset_Offset.cxx:483-488`).
int runOffsetSurf(const TopoDS_Shape& theShape, const double theOffset)
{
  std::cout << std::fixed << std::setprecision(9);
  const double aUs[3] = {0.4, 1.1, 2.7};
  const double aVs[3] = {-0.5, 0.3, 1.2};
  int          fi     = -1;
  for (TopExp_Explorer exp(theShape, TopAbs_FACE); exp.More(); exp.Next())
  {
    ++fi;
    const TopoDS_Face&        aFace = TopoDS::Face(exp.Current());
    TopLoc_Location           L;
    occ::handle<Geom_Surface> S = BRep_Tool::Surface(aFace, L);
    if (occ::handle<Geom_RectangularTrimmedSurface> RT =
          occ::down_cast<Geom_RectangularTrimmedSurface>(S))
    {
      S = RT->BasisSurface();
    }
    BRepOffset_Status         aStatus;
    occ::handle<Geom_Surface> anOff = BRepOffset::Surface(S, theOffset, aStatus);
    std::cout << "OFFSURF f=" << fi << " src=" << S->DynamicType()->Name()
              << " res=" << anOff->DynamicType()->Name()
              << " status=" << static_cast<int>(aStatus);
    for (int i = 0; i < 3; ++i)
    {
      const gp_Pnt aP = anOff->Value(aUs[i], aVs[i]);
      std::cout << " p" << i << "=(" << aP.X() << "," << aP.Y() << "," << aP.Z() << ")";
    }
    std::cout << std::endl;
  }
  return 0;
}

// --- TEMP B6d (`--ssx`): Extrema_ExtSS / Extrema_GenExtSS ground truth --------
// Prints what OCCT itself returns for the surface pairs the port's
// `occt_geom::extrema_ss::surface_surface_extrema_all` dispatches on, so the
// GenExtSS arm of `Extrema_ExtSS::Perform` can be compared number for number.
#include <Extrema_ExtSS.hxx>
#include <Extrema_GenExtSS.hxx>
#include <Extrema_POnSurf.hxx>
#include <Geom_SphericalSurface.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <GeomAPI_ExtremaSurfaceSurface.hxx>
#include <gp_Ax3.hxx>

static void dumpExtSS(const char* theTag, const Extrema_ExtSS& theExt)
{
  std::cout << "SSX " << theTag << " done=" << (theExt.IsDone() ? 1 : 0);
  if (!theExt.IsDone())
  {
    std::cout << std::endl;
    return;
  }
  std::cout << " isPar=" << (theExt.IsParallel() ? 1 : 0) << " nbExt=" << theExt.NbExt()
            << std::endl;
  for (int i = 1; i <= theExt.NbExt(); ++i)
  {
    std::cout << "SSX  " << theTag << "[" << i << "] sqDist=" << std::setprecision(17)
              << theExt.SquareDistance(i);
    if (!theExt.IsParallel())
    {
      double u1, v1, u2, v2;
      Extrema_POnSurf aP1, aP2;
      theExt.Points(i, aP1, aP2);
      aP1.Parameter(u1, v1);
      aP2.Parameter(u2, v2);
      const gp_Pnt aQ1 = aP1.Value();
      const gp_Pnt aQ2 = aP2.Value();
      std::cout << " uv1=(" << u1 << "," << v1 << ") uv2=(" << u2 << "," << v2 << ")"
                << " q1=(" << aQ1.X() << "," << aQ1.Y() << "," << aQ1.Z() << ")"
                << " q2=(" << aQ2.X() << "," << aQ2.Y() << "," << aQ2.Z() << ")";
    }
    std::cout << std::endl;
  }
}

int runExtSS()
{
  const double aTol = Precision::PConfusion();
  gp_Ax3 aAx;
  occ::handle<Geom_SphericalSurface> aSph1 = new Geom_SphericalSurface(aAx, 1.0);
  occ::handle<Geom_SphericalSurface> aSph2 = new Geom_SphericalSurface(aAx, 1.0);
  aSph2->Translate(gp_Vec(4.0, 0.0, 0.0));
  occ::handle<Geom_Plane>   aPln = new Geom_Plane(aAx);
  aPln->Translate(gp_Vec(0.0, 0.0, 5.0));

  GeomAdaptor_Surface anA1(aSph1);
  GeomAdaptor_Surface anA2(aSph2);
  GeomAdaptor_Surface anAp(aPln);

  std::cout << std::setprecision(17);
  std::cout << "SSX rng sph u=[" << anA1.FirstUParameter() << "," << anA1.LastUParameter()
            << "] v=[" << anA1.FirstVParameter() << "," << anA1.LastVParameter() << "]"
            << " pln u=[" << anAp.FirstUParameter() << "," << anAp.LastUParameter() << "] v=["
            << anAp.FirstVParameter() << "," << anAp.LastVParameter() << "]" << std::endl;
  std::cout << "SSX types s1=" << static_cast<int>(anA1.GetType())
            << " s2=" << static_cast<int>(anA2.GetType())
            << " pl=" << static_cast<int>(anAp.GetType()) << std::endl;

  dumpExtSS("sph-sph", Extrema_ExtSS(anA1, anA2, aTol, aTol));
  dumpExtSS("sph-pln", Extrema_ExtSS(anA1, anAp, aTol, aTol));
  dumpExtSS("pln-sph", Extrema_ExtSS(anAp, anA1, aTol, aTol));

  // The 12-argument constructor: explicit windows instead of the adaptor's own.
  dumpExtSS("sph-pln-win",
            Extrema_ExtSS(anA1,
                          anAp,
                          0.0,
                          2.0 * M_PI,
                          -M_PI / 2.0,
                          M_PI / 2.0,
                          -10.0,
                          10.0,
                          -10.0,
                          10.0,
                          aTol,
                          aTol));

  // GeomAPI_ExtremaSurfaceSurface (GeomAPI_ExtremaSurfaceSurface.cxx:60-69).
  GeomAPI_ExtremaSurfaceSurface aApi(aSph1, aPln);
  std::cout << "SSX api-sph-pln nbExt=" << aApi.NbExtrema();
  if (aApi.NbExtrema() > 0)
  {
    gp_Pnt aQ1, aQ2;
    aApi.NearestPoints(aQ1, aQ2);
    std::cout << " dist=" << std::setprecision(17) << aApi.LowerDistance()
              << " q1=(" << aQ1.X() << "," << aQ1.Y() << "," << aQ1.Z() << ")"
              << " q2=(" << aQ2.X() << "," << aQ2.Y() << "," << aQ2.Z() << ")";
  }
  std::cout << std::endl;

  // Raw Extrema_GenExtSS with the same NbU/NbV the Extrema_ExtSS arm uses.
  Extrema_GenExtSS aGen(anA1, anAp, 20, 20, aTol, aTol);
  std::cout << "SSX gen-sph-pln done=" << (aGen.IsDone() ? 1 : 0) << " nbExt=" << aGen.NbExt();
  for (int i = 1; i <= aGen.NbExt(); ++i)
  {
    std::cout << " [" << i << "] sqDist=" << std::setprecision(17) << aGen.SquareDistance(i);
  }
  std::cout << std::endl;
  return 0;
}
// --- end TEMP B6d -------------------------------------------------------------
