// Separate read-only probe: per-face wire/edge census on OCCT's imported shape.
// Built independently of occt_probe.cpp so concurrent edits cannot clash.
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
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
    std::cout << "usage: wires_probe <file.step> [faceIndex] [--nofix]\n";
    return 2;
  }
  bool aNoFix = false;
  for (int i = 1; i < argc; ++i)
  {
    if (std::string(argv[i]) == "--nofix") aNoFix = true;
  }
  STEPControl_Reader aReader;
  if (aReader.ReadFile(argv[1]) != IFSelect_RetDone)
  {
    std::cout << "read failed\n";
    return 1;
  }
  aReader.TransferRoots();
  TopoDS_Shape aShape = aReader.OneShape();
  if (!aNoFix)
  {
    // default: whatever the reader produced (ShapeFix already ran inside Transfer)
  }

  const int aWant = (argc > 2 && argv[2][0] != '-') ? std::atoi(argv[2]) : -1;
  TopExp_Explorer anEx(aShape, TopAbs_FACE);
  int aIdx = 0, aTotalWires = 0;
  std::map<int, int> aHist;
  for (; anEx.More(); anEx.Next(), ++aIdx)
  {
    const TopoDS_Face& aFace = TopoDS::Face(anEx.Current());
    int aWires = 0;
    std::string aPerWire;
    TopExp_Explorer aW(aFace, TopAbs_WIRE);
    for (; aW.More(); aW.Next())
    {
      ++aWires;
      int aEdges = 0;
      TopExp_Explorer aE(aW.Current(), TopAbs_EDGE);
      for (; aE.More(); aE.Next()) ++aEdges;
      aPerWire += " " + std::to_string(aEdges);
    }
    aTotalWires += aWires;
    aHist[aWires]++;
    if (aIdx == aWant || aWires >= 2)
    {
      std::cout << "FACE " << aIdx << " wires=" << aWires << " edges_per_wire:" << aPerWire << "\n";
    }
  }
  std::cout << "TOTAL faces=" << aIdx << " wires=" << aTotalWires << " hist(wires->faces):";
  for (const auto& kv : aHist) std::cout << " " << kv.first << "->" << kv.second;
  std::cout << "\n";
  return 0;
}
