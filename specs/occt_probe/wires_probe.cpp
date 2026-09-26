// Read-only probe: per-face wire census + UV-degenerate-wire census on OCCT's shape.
// A wire is "UV-degenerate" when its pcurve UV bounding box has ~zero area, which is
// the signature the port's mesh pipeline chokes on (zero-area loops -> open chains).
#include <BRepAdaptor_Surface.hxx>
#include <BRep_Tool.hxx>
#include <Geom2d_Curve.hxx>
#include <STEPControl_Reader.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Wire.hxx>
#include <cstdlib>
#include <iostream>
#include <map>

int main(int argc, char** argv)
{
  if (argc < 2)
  {
    std::cout << "usage: wires_probe <file.step> [faceIndex]\n";
    return 2;
  }
  STEPControl_Reader aReader;
  if (aReader.ReadFile(argv[1]) != IFSelect_RetDone)
  {
    std::cout << "read failed\n";
    return 1;
  }
  aReader.TransferRoots();
  const TopoDS_Shape aShape = aReader.OneShape();
  const int aWant = (argc > 2) ? std::atoi(argv[2]) : -1;

  TopExp_Explorer anEx(aShape, TopAbs_FACE);
  int aIdx = 0, aMul = 0, aDegen = 0, aTotalWires = 0;
  std::map<int, int> aHist;
  for (; anEx.More(); anEx.Next(), ++aIdx)
  {
    const TopoDS_Face& aFace = TopoDS::Face(anEx.Current());
    int aWires = 0;
    std::string aPerWire;
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
      if (aSpanU < 1e-9 || aSpanV < 1e-9) aFaceDegen = true;
    }
    aTotalWires += aWires;
    aHist[aWires]++;
    if (aWires >= 2) ++aMul;
    if (aFaceDegen) ++aDegen;
    if (aIdx == aWant || (aWires >= 2 && aFaceDegen))
    {
      std::cout << "FACE " << aIdx << " wires=" << aWires << " edges_per_wire:" << aPerWire
                << (aFaceDegen ? "  [UV-DEGENERATE-WIRE]" : "") << "\n";
    }
  }
  std::cout << "TOTAL faces=" << aIdx << " wires=" << aTotalWires << " multiwire_faces=" << aMul
            << " uv_degenerate_faces=" << aDegen << " hist:";
  for (const auto& kv : aHist) std::cout << " " << kv.first << "->" << kv.second;
  std::cout << "\n";
  return 0;
}
