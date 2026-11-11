// Scratch ground-truth probe: General Fuse split of (box - Z hole) and the
// perpendicular Y cylinder, mirroring the state `BOPAlgo_Builder::BuildSplitSolids`
// hands to `BOPAlgo_BOP::BuildRC`.
//
// Dumps, for every argument, the solids `Modified()` reports (the `myImages`
// entries filled by `BuildSplitSolids`, `BOPAlgo_Builder_3.cxx:413-618`) with
// their faces, so the Rust split can be compared face by face.
#include <BOPAlgo_BOP.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_BuilderAlgo.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBndLib.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopoDS.hxx>
#include <iostream>

namespace
{

const char* SurfName(const GeomAbs_SurfaceType theType)
{
  switch (theType)
  {
    case GeomAbs_Plane: return "Plane";
    case GeomAbs_Cylinder: return "Cylinder";
    case GeomAbs_Cone: return "Cone";
    case GeomAbs_Sphere: return "Sphere";
    case GeomAbs_Torus: return "Torus";
    default: return "Other";
  }
}

void DumpSolid(const TopoDS_Shape& theSolid, const char* theTag)
{
  GProp_GProps aProps;
  BRepGProp::VolumeProperties(theSolid, aProps);
  int aNbFaces = 0;
  for (TopExp_Explorer anExp(theSolid, TopAbs_FACE); anExp.More(); anExp.Next())
  {
    ++aNbFaces;
  }
  std::cout << theTag << ": type=" << theSolid.ShapeType() << " faces=" << aNbFaces
            << " vol=" << aProps.Mass() << "\n";
  for (TopExp_Explorer anExp(theSolid, TopAbs_FACE); anExp.More(); anExp.Next())
  {
    const TopoDS_Face& aF = TopoDS::Face(anExp.Current());
    GProp_GProps anArea;
    BRepGProp::SurfaceProperties(aF, anArea);
    BRepAdaptor_Surface anAdaptor(aF);
    int aNbWires = 0;
    for (TopExp_Explorer aW(aF, TopAbs_WIRE); aW.More(); aW.Next())
    {
      ++aNbWires;
    }
    int aNbEdges = 0;
    for (TopExp_Explorer anE(aF, TopAbs_EDGE); anE.More(); anE.Next())
    {
      ++aNbEdges;
    }
    Bnd_Box aBox;
    BRepBndLib::Add(aF, aBox);
    double aXmin = 0., aYmin = 0., aZmin = 0., aXmax = 0., aYmax = 0., aZmax = 0.;
    if (!aBox.IsVoid())
    {
      aBox.Get(aXmin, aYmin, aZmin, aXmax, aYmax, aZmax);
    }
    const char* anOri = "F";
    switch (anExp.Current().Orientation())
    {
      case TopAbs_REVERSED: anOri = "R"; break;
      case TopAbs_INTERNAL: anOri = "I"; break;
      case TopAbs_EXTERNAL: anOri = "E"; break;
      default: break;
    }
    std::cout << "   " << SurfName(anAdaptor.GetType()) << " ori=" << anOri
              << " area=" << anArea.Mass()
              << " wires=" << aNbWires << " edges=" << aNbEdges << " bb=(" << aXmin << ","
              << aYmin << "," << aZmin << ")-(" << aXmax << "," << aYmax << "," << aZmax
              << ")\n";
  }
}

void DumpImages(BRepAlgoAPI_BuilderAlgo& theGF,
                const TopoDS_Shape&            theSource,
                const char*                    theTag)
{
  const NCollection_List<TopoDS_Shape>& anImages = theGF.Modified(theSource);
  std::cout << theTag << ": images=" << anImages.Extent() << "\n";
  int anIdx = 0;
  for (NCollection_List<TopoDS_Shape>::Iterator anIt(anImages); anIt.More(); anIt.Next())
  {
    DumpSolid(anIt.Value(), ("  " + std::to_string(anIdx++)).c_str());
  }
}

void DumpImageEntry(const NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>& theMap,
                    const TopoDS_Shape&                        theSource,
                    const char*                                theTag)
{
  if (!theMap.IsBound(theSource))
  {
    std::cout << theTag << ": not bound\n";
    return;
  }
  const NCollection_List<TopoDS_Shape>& aList = theMap.Find(theSource);
  std::cout << theTag << ": images=" << aList.Extent() << "\n";
  int anIdx = 0;
  for (NCollection_List<TopoDS_Shape>::Iterator anIt(aList); anIt.More(); anIt.Next())
  {
    DumpSolid(anIt.Value(), ("  " + std::to_string(anIdx++)).c_str());
  }
}

} // namespace

int main()
{  std::cout.setf(std::ios::unitbuf);
  const TopoDS_Shape aBox =
    BRepPrimAPI_MakeBox(gp_Pnt(0., 0., 0.), gp_Pnt(2., 2., 1.)).Shape();
  const TopoDS_Shape aCylZ =
    BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(1., 1., -0.5), gp_Dir(0., 0., 1.)), 0.3, 2.0).Shape();
  const TopoDS_Shape aCylY =
    BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(1., -2., 0.5), gp_Dir(0., 1., 0.)), 0.3, 4.0).Shape();

  BRepAlgoAPI_Cut aCut1(aBox, aCylZ);
  aCut1.SetRunParallel(false);
  aCut1.Build();
  if (!aCut1.IsDone())
  {
    std::cout << "cut1 failed\n";
    return 1;
  }

  BRepAlgoAPI_BuilderAlgo aGF;
  TopTools_ListOfShape    aArgs;
  aArgs.Append(aCut1.Shape());
  aArgs.Append(aCylY);
  aGF.SetArguments(aArgs);
  aGF.SetRunParallel(false);
  aGF.SetToFillHistory(true);
  aGF.Build();
  if (aGF.HasErrors())
  {
    std::cout << "GF failed\n";
    aGF.DumpErrors(std::cout);
    return 1;
  }
  DumpSolid(aGF.Shape(), "GF shape");
  DumpImages(aGF, aCut1.Shape(), "arg (box - Zcyl)");
  DumpImages(aGF, aCylY, "tool (Ycyl)");

  BRepAlgoAPI_Cut aCut2(aCut1.Shape(), aCylY);
  aCut2.SetRunParallel(false);
  aCut2.SetToFillHistory(true);
  aCut2.Build();
  if (!aCut2.IsDone())
  {
    std::cout << "cut2 failed\n";
    return 1;
  }
  DumpSolid(aCut2.Shape(), "cut2 shape");
  DumpImages(aCut2, aCut1.Shape(), "cut arg images");
  DumpImages(aCut2, aCylY, "cut tool images");

  // myImages of the building phase (`BuildSplitSolids` output) for the Cut.
  {
    BOPAlgo_BOP aBOP;
    TopTools_ListOfShape aObj;
    aObj.Append(aCut1.Shape());
    aBOP.SetArguments(aObj);
    aBOP.AddTool(aCylY);
    aBOP.SetOperation(BOPAlgo_CUT);
    aBOP.SetRunParallel(false);
    aBOP.Perform(Message_ProgressRange());
    if (aBOP.HasErrors())
    {
      std::cout << "BOP errors\n";
      aBOP.DumpErrors(std::cout);
      return 1;
    }
    DumpSolid(aBOP.Shape(), "BOP shape");
    const NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>& aMap = aBOP.Images();
    std::cout << "myImages extent=" << aMap.Extent() << "\n";
    DumpImageEntry(aMap, aCut1.Shape(), "myImages[arg]");
    DumpImageEntry(aMap, aCylY, "myImages[tool]");
    for (NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>::Iterator anIt(aMap); anIt.More(); anIt.Next())
    {
      if (anIt.Key().ShapeType() != TopAbs_SOLID)
      {
        continue;
      }
      DumpImageEntry(aMap, anIt.Key(), "myImages[solid]");
    }
  }
  // Case 2: a box with an internal cavity (Cut of a fully contained box).
  {
    std::cout << "\n=== internal box cut ===\n";
    const TopoDS_Shape anOuter =
      BRepPrimAPI_MakeBox(gp_Pnt(0., 0., 0.), gp_Pnt(1., 1., 1.)).Shape();
    const TopoDS_Shape anInner =
      BRepPrimAPI_MakeBox(gp_Pnt(0.2, 0.2, 0.2), gp_Pnt(0.8, 0.8, 0.8)).Shape();
    BOPAlgo_BOP aBOP;
    TopTools_ListOfShape aObj;
    aObj.Append(anOuter);
    aBOP.SetArguments(aObj);
    aBOP.AddTool(anInner);
    aBOP.SetOperation(BOPAlgo_CUT);
    aBOP.SetRunParallel(false);
    aBOP.Perform(Message_ProgressRange());
    if (aBOP.HasErrors())
    {
      std::cout << "internal BOP errors\n";
      aBOP.DumpErrors(std::cout);
    }
    else
    {
      DumpSolid(aBOP.Shape(), "internal BOP shape");
      const NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>&
        aMap = aBOP.Images();
      DumpImageEntry(aMap, anOuter, "internal myImages[arg]");
      DumpImageEntry(aMap, anInner, "internal myImages[tool]");
      // Shell structure of every image of the argument.
      if (aMap.IsBound(anOuter))
      {
        int anIdx = 0;
        for (NCollection_List<TopoDS_Shape>::Iterator anIt(aMap.Find(anOuter)); anIt.More();
             anIt.Next(), ++anIdx)
        {
          const TopoDS_Shape& aSolid = anIt.Value();
          std::cout << "image " << anIdx << " shells:\n";
          for (TopExp_Explorer aShExp(aSolid, TopAbs_SHELL); aShExp.More(); aShExp.Next())
          {
            const TopoDS_Shape& aShell = aShExp.Current();
            int                 aNbF  = 0;
            for (TopExp_Explorer aFExp(aShell, TopAbs_FACE); aFExp.More(); aFExp.Next())
            {
              ++aNbF;
            }
            std::cout << "   shell ori=" << (int)aShell.Orientation() << " faces=" << aNbF
                      << "\n";
            for (TopExp_Explorer aFExp(aShell, TopAbs_FACE); aFExp.More(); aFExp.Next())
            {
              const TopoDS_Shape& aF = aFExp.Current();
              Bnd_Box             aBox;
              BRepBndLib::Add(aF, aBox);
              double aXmin = 0., aYmin = 0., aZmin = 0., aXmax = 0., aYmax = 0., aZmax = 0.;
              if (!aBox.IsVoid())
              {
                aBox.Get(aXmin, aYmin, aZmin, aXmax, aYmax, aZmax);
              }
              std::cout << "      face ori=" << (int)aF.Orientation() << " bb=(" << aXmin << ","
                        << aYmin << "," << aZmin << ")-(" << aXmax << "," << aYmax << ","
                        << aZmax << ")\n";
            }
          }
        }
      }
    }
  }
  return 0;
}

