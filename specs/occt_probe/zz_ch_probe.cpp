// Scratch ground truth for `BRepFeat_MakeCylindricalHole` on OCCT 8.0.0:
// the five entry points (`Perform(Radius)`, `Perform(Radius, PFrom, PTo, Cont)`,
// `PerformThruNext`, `PerformUntilEnd`, `PerformBlind`) plus `Build` on a 2x2x1
// box, dumping status, volume, and the per-face area/wire/edge table.
#include <BRepAdaptor_Surface.hxx>
#include <BRepFeat_MakeCylindricalHole.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box2d.hxx>
#include <GProp_GProps.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Iterator.hxx>
#include <gp_Ax1.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <iostream>

namespace
{

const char* SurfName(GeomAbs_SurfaceType theT)
{
  switch (theT)
  {
    case GeomAbs_Plane: return "Plane";
    case GeomAbs_Cylinder: return "Cylinder";
    case GeomAbs_Cone: return "Cone";
    case GeomAbs_Sphere: return "Sphere";
    case GeomAbs_Torus: return "Torus";
    default: return "Other";
  }
}

void DumpResult(const TopoDS_Shape& theS, const char* theTag, const int theStatus)
{
  if (theS.IsNull())
  {
    std::cout << "=== " << theTag << ": status=" << theStatus << " shape=null\n";
    return;
  }
  GProp_GProps aV;
  BRepGProp::VolumeProperties(theS, aV);
  int aNbFaces = 0;
  for (TopExp_Explorer anEx(theS, TopAbs_FACE); anEx.More(); anEx.Next())
  {
    ++aNbFaces;
  }
  int aNbShells = 0;
  for (TopExp_Explorer anEx(theS, TopAbs_SHELL); anEx.More(); anEx.Next())
  {
    ++aNbShells;
  }
  int aNbSolids = 0;
  for (TopExp_Explorer anEx(theS, TopAbs_SOLID); anEx.More(); anEx.Next())
  {
    ++aNbSolids;
  }
  std::cout << "\n=== " << theTag << ": status=" << theStatus << " volume=" << aV.Mass()
            << " faces=" << aNbFaces << " shells=" << aNbShells << " solids=" << aNbSolids
            << "\n";
  int aI = 0;
  for (TopExp_Explorer anEx(theS, TopAbs_FACE); anEx.More(); anEx.Next(), ++aI)
  {
    const TopoDS_Face aF = TopoDS::Face(anEx.Current());
    BRepAdaptor_Surface aS(aF);
    GProp_GProps aA;
    BRepGProp::SurfaceProperties(aF, aA);
    int aNw = 0;
    int aNe = 0;
    for (TopoDS_Iterator aWI(aF); aWI.More(); aWI.Next(), ++aNw)
    {
      for (TopoDS_Iterator aEI(aWI.Value()); aEI.More(); aEI.Next())
      {
        ++aNe;
      }
    }
    Bnd_Box2d aUV;
    BRepTools::AddUVBounds(aF, aUV);
    double u1, u2, v1, v2;
    aUV.Get(u1, v1, u2, v2);
    std::cout << "  F" << aI << " " << SurfName(aS.GetType())
              << " ori=" << (int)aF.Orientation() << " area=" << aA.Mass() << " wires=" << aNw
              << " edges=" << aNe << " uv=(" << u1 << "," << u2 << "," << v1 << "," << v2
              << ")\n";
  }
}

} // namespace

int main()
{
  std::cout.setf(std::ios::unitbuf);
  const TopoDS_Shape aBox =
    BRepPrimAPI_MakeBox(gp_Pnt(0., 0., 0.), gp_Pnt(2., 2., 1.)).Shape();
  const double aR = 0.3;

  // Axis from below the box, pointing up: parameters 0.5 (z=0) .. 1.5 (z=1).
  const gp_Ax1 anAxisUp(gp_Pnt(1., 1., -0.5), gp_Dir(0., 0., 1.));
  // Axis from above the box, pointing down: parameters 0.2 (z=1) .. 1.2 (z=0).
  const gp_Ax1 anAxisDown(gp_Pnt(1., 1., 1.2), gp_Dir(0., 0., -1.));

  {
    BRepFeat_MakeCylindricalHole aHole;
    aHole.Init(aBox, anAxisUp);
    aHole.Perform(aR);
    aHole.Build();
    DumpResult(aHole.Shape(), "Perform(Radius)", (int)aHole.Status());
  }
  {
    BRepFeat_MakeCylindricalHole aHole;
    aHole.Init(aBox, anAxisUp);
    aHole.Perform(aR, 0.5, 1.0, true);
    aHole.Build();
    DumpResult(aHole.Shape(), "Perform(Radius,0.5,1.0)", (int)aHole.Status());
  }
  {
    BRepFeat_MakeCylindricalHole aHole;
    aHole.Init(aBox, anAxisUp);
    aHole.Perform(aR, 1.0, 0.5, true);
    aHole.Build();
    DumpResult(aHole.Shape(), "Perform(Radius,1.0,0.5)", (int)aHole.Status());
  }
  {
    BRepFeat_MakeCylindricalHole aHole;
    aHole.Init(aBox, anAxisUp);
    aHole.PerformThruNext(aR, true);
    aHole.Build();
    DumpResult(aHole.Shape(), "PerformThruNext", (int)aHole.Status());
  }
  {
    BRepFeat_MakeCylindricalHole aHole;
    aHole.Init(aBox, anAxisUp);
    aHole.PerformUntilEnd(aR, true);
    aHole.Build();
    DumpResult(aHole.Shape(), "PerformUntilEnd", (int)aHole.Status());
  }
  {
    BRepFeat_MakeCylindricalHole aHole;
    aHole.Init(aBox, anAxisDown);
    aHole.PerformBlind(aR, 0.8, true);
    aHole.Build();
    DumpResult(aHole.Shape(), "PerformBlind(0.8)", (int)aHole.Status());
  }
  {
    // `Length` beyond the far side: the reference answers `HoleTooLong`.
    BRepFeat_MakeCylindricalHole aHole;
    aHole.Init(aBox, anAxisDown);
    aHole.PerformBlind(aR, 2.0, true);
    aHole.Build();
    DumpResult(aHole.Shape(), "PerformBlind(2.0)", (int)aHole.Status());
  }
  return 0;
}
