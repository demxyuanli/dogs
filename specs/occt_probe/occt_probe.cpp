// OCCT ground-truth probe (T-69/T-80/T-87/T-05 oracle).
// Reads a STEP file with OCCT 8.0.0 and prints what the port must reproduce:
//   * BRepCheck_Analyzer verdict
//   * BRepGProp::VolumeProperties (the OCCT volume, analytic by default)
//   * per-face surface type, orientation, composed normal and outward sign
//   * cylinder/spline frames, so the port's reconstruction can be compared
#include <BRepAdaptor_Surface.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
// T-99: `BRepMesh_DefaultRangeSplitter` exposes `GetRangeU/GetRangeV/GetDelta/
// GetToleranceUV` publicly, which is all the diagnostic below needs. Its
// `computeLengthU/V` are declared `private` in the shipped header and their symbols
// are emitted with private access, so they cannot be called from outside without
// either editing the installed header or recompiling the library - the same length
// values are already reported by the `--uvsum` dump, so they are not repeated here.
#include <BRepMesh_GeomTool.hxx>
#include <BRepMesh_DefaultRangeSplitter.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepMesh_ModelBuilder.hxx>
#include <BRepMesh_EdgeDiscret.hxx>
#include <BRepMesh_ModelHealer.hxx>
#include <ShapeAnalysis_Wire.hxx>
#include <ShapeAnalysis_WireOrder.hxx>
#include <ShapeExtend_WireData.hxx>
#include <ShapeFix_Wire.hxx>
#include <ShapeExtend.hxx>
#include <IMeshData_Model.hxx>
#include <IMeshData_Face.hxx>
#include <IMeshData_Wire.hxx>
#include <IMeshData_Edge.hxx>
#include <IMeshData_PCurve.hxx>
#include <IMeshTools_Parameters.hxx>
#include <Message_ProgressRange.hxx>
// T-104: OCCT's own node-insertion algo, subclassed to dump the structure it
// hands to `BRepMesh_Delaun` (`--structdump`). Kept after the `IMeshData_*`
// includes because `BRepMesh_NodeInsertionMeshAlgo.hxx` needs those types.
#include <BRepMesh_DelaunayNodeInsertionMeshAlgo.hxx>
#include <BRepMesh_DelaunayBaseMeshAlgo.hxx>
// T-123: the `--algospy` replay subclasses the real node-insertion algo, so the face
// discretiser that owns the algo factory, the factory itself, the algorithm base and the
// Delaunay structure the dump reports have to be complete types.
#include <BRepMesh_FaceDiscret.hxx>
#include <BRepMesh_MeshAlgoFactory.hxx>
#include <BRepMesh_DataStructureOfDelaun.hxx>
// TEMP DCPASS: the per-pass instrumentation below needs the private loop state
// (`myControlNodes` / `myCouplesMap` / `myMaxSqDeflection`) of the template
// header, so it is included with access control disabled. The class layout is
// untouched (access specifiers do not affect it) and the instantiation is local.
#define private public
#include <BRepMesh_DelaunayDeflectionControlMeshAlgo.hxx>
#undef private
#include <BRepMesh_DefaultRangeSplitter.hxx>
#include <BRepMesh_BoundaryParamsRangeSplitter.hxx>
#include <BRepMesh_ExtrusionRangeSplitter.hxx>
#include <BRepMesh_NURBSRangeSplitter.hxx>
#include <BRepMesh_UndefinedRangeSplitter.hxx>
#include <IGESControl_Controller.hxx>
#include <IGESControl_Reader.hxx>
#include <IGESControl_Writer.hxx>
#include <BRepMesh_TorusRangeSplitter.hxx>
#include <BRepMesh_SphereRangeSplitter.hxx>
#include <BRepMesh_CylinderRangeSplitter.hxx>
#include <BRepMesh_ConeRangeSplitter.hxx>
#include <Precision.hxx>
#include <Poly_Triangulation.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <cstdlib>
#include <filesystem>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <BRepGProp_Face.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <STEPControl_Reader.hxx>
#include <XSControl_WorkSession.hxx>
#include <XSControl_TransferReader.hxx>
#include <Interface_InterfaceModel.hxx>
#include <Interface_Static.hxx>
#include <ShapeProcess.hxx>
#include <ShapeAnalysis.hxx>
#include <ShapeAnalysis_Edge.hxx>
#include <ShapeAnalysis_Surface.hxx>
#include <ShapeConstruct_ProjectCurveOnSurface.hxx>
#include <DE_ShapeFixParameters.hxx>
#include <StepData_StepModel.hxx>
#include <BRepGProp_Domain.hxx>
#include <BRepGProp_Vinert.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BOPAlgo_PaveFiller.hxx>
#include <BOPDS_DS.hxx>
#include <BOPAlgo_Builder.hxx>
#include <BOPTools_AlgoTools.hxx>
#include <BOPTools_CoupleOfShape.hxx>
#include <IntTools_Context.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BOPTools_AlgoTools3D.hxx>
#include <BOPTools_AlgoTools2D.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <Geom_Plane.hxx>
#include <BRep_Tool.hxx>
#include <BRepTools.hxx>
#include <TopTools_ListOfShape.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <StepShape_FaceSurface.hxx>
#include <string>
#include <TopExp_Explorer.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopExp.hxx>
#include <TopoDS_Vertex.hxx>
#include <Geom_Curve.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

#include <Geom2d_Curve.hxx>
#include <Geom2d_Line.hxx>
#include <Bnd_Box2d.hxx>
// TEMP FCLS2D: `BRepTopAdaptor_FClass2d` boundary ring / state ground truth.
#include <BRepTopAdaptor_FClass2d.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2dAdaptor_Curve.hxx>
#include <gp_Pnt2d.hxx>
#include <algorithm>
#include <cmath>
#include <map>
#include <unordered_map>
#include <vector>

#include <cstdio>
#include <cstdlib>
#include <iostream>

// --- TEMP T-99 feasibility probe (`--delauncheck`) -------------------------------
// §9.365 established that OCCT meshes all 226 faces (min 2 triangles) while the port's
// `Delaun::add_vertices` yields 0 domain elements on 16 of them, so the next
// measurement must be OCCT's *Delaunay state*, not its product. The node-insertion
// algorithm is templated and not directly hookable, so the clean route is to feed
// `BRepMesh_Delaun` the port's own registered points. This block only proves that the
// class is reachable/linkable from this probe's link set; it constructs a Delaun over
// a tiny synthetic vertex array and reports the resulting structure size.
#include <BRepMesh_CircleTool.hxx>
#include <BRepMesh_Delaun.hxx>
#include <BRepMesh_Vertex.hxx>
#include <BRepMesh_DataStructureOfDelaun.hxx>
#include <BRepMesh_VertexTool.hxx>
#include <BRepMesh_Context.hxx>
#include <IMeshTools_Context.hxx>
#include <NCollection_IncAllocator.hxx>
#include <IMeshData_Types.hxx>
#include <array>
#include <fstream>
#include <iomanip>
#include <sstream>

static int runDelaunCheck()
{
  IMeshData::Array1OfVertexOfDelaun aVerts(1, 5);
  aVerts(1) = BRepMesh_Vertex(gp_XY(0.0, 0.0), 1, BRepMesh_Frontier);
  aVerts(2) = BRepMesh_Vertex(gp_XY(10.0, 0.0), 2, BRepMesh_Frontier);
  aVerts(3) = BRepMesh_Vertex(gp_XY(10.0, 10.0), 3, BRepMesh_Frontier);
  aVerts(4) = BRepMesh_Vertex(gp_XY(0.0, 10.0), 4, BRepMesh_Frontier);
  aVerts(5) = BRepMesh_Vertex(gp_XY(5.0, 5.0), 5, BRepMesh_Free);

  BRepMesh_Delaun aDelaun(aVerts);
  const occ::handle<BRepMesh_DataStructureOfDelaun>& aSt = aDelaun.Result();
  // NOTE: `aDelaun.Frontier()` cannot be used from here - it calls the private
  // `getEdgesByType`, whose access level mangles the symbol and breaks the link.
  // `Result()` and the `BRepMesh_DataStructureOfDelaun` accessors are public.
  std::cout << "DELAUNCHECK nodes=" << aSt->NbNodes() << " links=" << aSt->NbLinks()
            << " domain=" << aSt->ElementsOfDomain().Extent() << std::endl;
  return 0;
}
// --- end TEMP T-99 ---------------------------------------------------------------

// --- TEMP T-99 (`--delaunfeed <dir>`) ---------------------------------------------
// SUPERSEDED by `--delaunstruct` below (§9.368). Feeding only the points leaves the
// structure with **zero** frontier links, and `cleanupMesh` (`BRepMesh_Delaun.cxx:1028`,
// which skips only Frontier links at `:832-835`) then deletes every triangle, so the
// `domain=0` this mode prints is an artifact of the missing constraints - it is **not**
// evidence that the point set cannot be triangulated. Kept only as a record of the
// §9.367 run; use `--delaunstruct` for comparisons.
// Feeds the port's own registered Delaun points (dumped by `OCCT_TOPO_DUMP_DELAUN`)
// straight into OCCT's `BRepMesh_Delaun`, bypassing the reader and the whole
// `IMeshData` pipeline. §9.365 showed OCCT meshes all 226 faces (min 2 triangles)
// while the port's Delaun yields 0 domain elements on 16 of them; this decides
// whether that is a data problem or a code problem.
static int runDelaunFeed(const std::string& theDir)
{
  for (int f : {169, 174, 189, 192, 193, 194, 195, 196, 197, 198, 199, 204, 206, 209, 210, 212})
  {
    char aName[512];
    std::snprintf(aName, sizeof(aName), "%s/delaun_f%d.txt", theDir.c_str(), f);
    std::ifstream aIn(aName);
    if (!aIn)
    {
      std::cout << "DELAUNFEED f=" << f << " (missing " << aName << ")" << std::endl;
      continue;
    }
    std::vector<std::array<double, 3>> aPts;
    std::string aLine;
    while (std::getline(aIn, aLine))
    {
      if (aLine.empty() || aLine[0] == '#')
      {
        continue;
      }
      std::istringstream aSS(aLine);
      double u = 0.0, v = 0.0; int m = 0;
      if (aSS >> u >> v >> m)
      {
        aPts.push_back({u, v, static_cast<double>(m)});
      }
    }
    if (aPts.empty())
    {
      std::cout << "DELAUNFEED f=" << f << " (no points)" << std::endl;
      continue;
    }

    IMeshData::Array1OfVertexOfDelaun aVerts(1, static_cast<int>(aPts.size()));
    for (size_t i = 0; i < aPts.size(); ++i)
    {
      // movability: 5 == BRepMesh_Frontier (see VertexState ordering in the port);
      // anything else is treated as a free/interior vertex.
      const BRepMesh_DegreeOfFreedom aMov =
        (static_cast<int>(aPts[i][2]) == 5) ? BRepMesh_Frontier : BRepMesh_Free;
      aVerts(static_cast<int>(i) + 1) =
        BRepMesh_Vertex(gp_XY(aPts[i][0], aPts[i][1]), static_cast<int>(i) + 1, aMov);
    }

    // `ProcessConstraints()` cannot be called from here: it is an inlined header
    // method whose body references the private `insertInternalEdges()` /
    // `frontierAdjust()`, so the access level mangles the symbols and the link fails
    // (same trap as `Frontier()`, see §9.366). The constructor already runs
    // `Init` -> `perform` -> `createTrianglesOnNewVertices`, which is the comparable
    // step; constraint handling would additionally need the port's link set.
    BRepMesh_Delaun aDelaun(aVerts);
    const occ::handle<BRepMesh_DataStructureOfDelaun>& aSt = aDelaun.Result();

    // Link states, counted from the public accessors (Frontier()/FreeEdges() would
    // call the private getEdgesByType and fail to link, see §9.366).
    int aFrontier = 0, aFixed = 0, aFree = 0, aDeleted = 0;
    for (int e = 1; e <= aSt->NbLinks(); ++e)
    {
      switch (aSt->GetLink(e).Movability())
      {
        case BRepMesh_Frontier: ++aFrontier; break;
        case BRepMesh_Fixed:    ++aFixed;    break;
        case BRepMesh_Free:     ++aFree;     break;
        default:                ++aDeleted;  break;
      }
    }

    std::cout << "DELAUNFEED f=" << f << " pts=" << aPts.size()
              << " nodes=" << aSt->NbNodes() << " links=" << aSt->NbLinks()
              << " domain=" << aSt->ElementsOfDomain().Extent()
              << " frontier=" << aFrontier << " fixed=" << aFixed
              << " free=" << aFree << " deleted=" << aDeleted << std::endl;
  }
  return 0;
}
// --- TEMP T-99 (`--delaunstruct <dir>`) -------------------------------------------
// Feeds OCCT's `BRepMesh_Delaun` the *identical* structure the port builds: same
// nodes (face-basis coordinates), same constraint links, same tolerance/cell size
// and the same vertex order (dumped by `OCCT_TOPO_DUMP_DELAUN`). §9.367 could only
// feed the bare points (0 constraints), and with no frontier links `cleanupMesh`
// deletes every triangle, so that run's `domain=0` said nothing about the points.
// This closes the gap: constraints are inserted through the public
// `AddLink` + `BRepMesh_Delaun(structure, indices, cellsU, cellsV)` constructor,
// which reaches `ProcessConstraints()` as a member call (§9.366's LNK2019 trap is
// avoided entirely).
static int runDelaunStruct(const std::string& theDir)
{
  // TEMP T-123: the face list is scanned from the directory (`delaun_in_f<face>.txt`)
  // instead of being hard-coded, so replaying a newly tagged `DELSTAGE=<face>` run needs
  // no probe edit (T-121's fixed list only covered Shape-2 / T-119 faces).
  std::vector<int> aFaces;
  {
    const std::string aPrefix  = "delaun_in_f";
    const std::string aSuffix  = ".txt";
    std::error_code   anErr;
    for (const std::filesystem::directory_entry& anEntry :
         std::filesystem::directory_iterator(theDir, anErr))
    {
      const std::string aFileName = anEntry.path().filename().string();
      if (aFileName.size() <= aPrefix.size() + aSuffix.size()
          || aFileName.compare(0, aPrefix.size(), aPrefix) != 0
          || aFileName.compare(aFileName.size() - aSuffix.size(), aSuffix.size(), aSuffix) != 0)
      {
        continue;
      }
      aFaces.push_back(std::atoi(
        aFileName.substr(aPrefix.size(), aFileName.size() - aPrefix.size() - aSuffix.size())
          .c_str()));
    }
  }
  std::sort(aFaces.begin(), aFaces.end());
  for (int f : aFaces)
  {
    char aName[512];
    std::snprintf(aName, sizeof(aName), "%s/delaun_in_f%d.txt", theDir.c_str(), f);
    std::ifstream aIn(aName);
    if (!aIn)
    {
      std::cout << "DELAUNSTRUCT f=" << f << " (missing " << aName << ")" << std::endl;
      continue;
    }

    std::vector<std::array<double, 3>> aNodes; // x, y, movability
    std::vector<std::array<int, 3>>    aLinks; // first, last, movability
    std::vector<int>                   aOrder;
    std::vector<int>                   aOrder1;
    std::vector<int>                   aOrder2;
    double                             aTolU = 0.0, aTolV = 0.0;
    int                                aCellsU = -1, aCellsV = -1;
    std::string                        aLine;
    while (std::getline(aIn, aLine))
    {
      if (aLine.empty() || aLine[0] == '#')
      {
        continue;
      }
      std::istringstream aSS(aLine);
      std::string        aTag;
      aSS >> aTag;
      if (aTag == "TOL")
      {
        aSS >> aTolU >> aTolV;
      }
      else if (aTag == "CELLS")
      {
        aSS >> aCellsU >> aCellsV;
      }
      else if (aTag == "N")
      {
        int i, m;
        double x, y;
        aSS >> i >> x >> y >> m;
        aNodes.push_back({x, y, static_cast<double>(m)});
      }
      else if (aTag == "L")
      {
        int i, a, b, m;
        aSS >> i >> a >> b >> m;
        aLinks.push_back({a, b, m});
      }
      else if (aTag == "V")
      {
        int v;
        while (aSS >> v)
        {
          aOrder.push_back(v);
        }
      }
      else if (aTag == "V1")
      {
        int v;
        while (aSS >> v)
        {
          aOrder1.push_back(v);
        }
      }
      else if (aTag == "V2")
      {
        int v;
        while (aSS >> v)
        {
          aOrder2.push_back(v);
        }
      }
      else if (aTag == "AFTER_CTOR" || aTag == "FINAL")
      {
        // Echo what the port itself measured on this very structure.
        std::cout << "PORT   f=" << f << " " << aTag;
        std::string aTok;
        while (aSS >> aTok)
        {
          const size_t aEq = aTok.find('=');
          if (aEq != std::string::npos)
          {
            std::cout << " " << aTok;
          }
        }
        std::cout << std::endl;
      }
    }

    if (aNodes.empty())
    {
      std::cout << "DELAUNSTRUCT f=" << f << " (no nodes)" << std::endl;
      continue;
    }
    std::cout << "TRACE f=" << f << " parsed nodes=" << aNodes.size()
              << " links=" << aLinks.size() << " v1=" << aOrder1.size()
              << " v2=" << aOrder2.size() << std::endl;

    occ::handle<NCollection_IncAllocator>       anAlloc = new NCollection_IncAllocator(IMeshData::MEMORY_BLOCK_SIZE_HUGE);
    occ::handle<BRepMesh_DataStructureOfDelaun> aStruct =
      new BRepMesh_DataStructureOfDelaun(anAlloc, static_cast<int>(aNodes.size()));
    aStruct->Data()->SetTolerance(aTolU, aTolV);
    aStruct->Data()->SetCellSize(14.0 * aTolU, 14.0 * aTolV);

    std::vector<int> aMap(aNodes.size() + 1, 0);
    const std::vector<int>& aOrder1Ref = aOrder1;
    const std::vector<int>& aOrder2Ref = aOrder2;

    if (!aOrder2Ref.empty())
    {
      // Two-phase replay of `BRepMesh_DelaunayBaseMeshAlgo::generateMesh` +
      // `BRepMesh_DelaunayNodeInsertionMeshAlgo::postProcessMesh`: the base
      // mesher is built over the phase-1 vertices, `EraseFreeLinks` runs, and
      // only then are the surface nodes registered and passed to
      // `AddVertices`. Same structure, same order, same cells -> any triangle
      // difference is a port bug, not an input difference.
      for (int v : aOrder1Ref)
      {
        const BRepMesh_Vertex aVertex(gp_XY(aNodes[v - 1][0], aNodes[v - 1][1]),
                                      v - 1,
                                      static_cast<BRepMesh_DegreeOfFreedom>(
                                        static_cast<int>(aNodes[v - 1][2])));
        aMap[v] = aStruct->AddNode(aVertex);
      }
      std::cout << "TRACE f=" << f << " phase1 nodes added, links=" << aStruct->NbLinks()
                << std::endl;
      for (const std::array<int, 3>& aL : aLinks)
      {
        // Only constraint links are fed here: the auxiliary links dumped from the
        // port (`Deleted`, still holding the recycled super-vertex ids) are
        // re-created by this replay's own base-mesher run.
        if (aL[2] == static_cast<int>(BRepMesh_Deleted) || aMap[aL[0]] == 0 || aMap[aL[1]] == 0)
        {
          continue;
        }
        aStruct->AddLink(BRepMesh_Edge(aMap[aL[0]],
                                       aMap[aL[1]],
                                       static_cast<BRepMesh_DegreeOfFreedom>(aL[2])));
      }
      std::cout << "TRACE f=" << f << " links added=" << aStruct->NbLinks() << std::endl;

      IMeshData::VectorOfInteger aIndices1;
      for (int v : aOrder1Ref)
      {
        aIndices1.Append(aMap[v]);
      }

      std::cout << "TRACE f=" << f << " entering base mesher" << std::endl;
      BRepMesh_Delaun aDelaun(aStruct, aIndices1, aCellsU, aCellsV);
      std::cout << "TRACE f=" << f << " base mesher done" << std::endl;
      std::cout << "OCCTBASE f=" << f << " nodes=" << aStruct->NbNodes()
                << " links=" << aStruct->NbLinks()
                << " domain=" << aStruct->ElementsOfDomain().Extent() << std::endl;

      for (int v : aOrder2Ref)
      {
        const BRepMesh_Vertex aVertex(gp_XY(aNodes[v - 1][0], aNodes[v - 1][1]),
                                      v - 1,
                                      static_cast<BRepMesh_DegreeOfFreedom>(
                                        static_cast<int>(aNodes[v - 1][2])));
        aMap[v] = aStruct->AddNode(aVertex);
      }
      std::cout << "TRACE f=" << f << " phase2 nodes added=" << aStruct->NbNodes() << std::endl;

      IMeshData::VectorOfInteger aIndices2;
      for (int v : aOrder2Ref)
      {
        aIndices2.Append(aMap[v]);
      }
      std::cout << "TRACE f=" << f << " calling AddVertices n=" << aIndices2.Length() << std::endl;
      aDelaun.AddVertices(aIndices2);
      std::cout << "TRACE f=" << f << " AddVertices done" << std::endl;
      // TEMP T-119: `AddVertices` heap sorts the index vector in place
      // (`cxx:2415-2417`), so the vector now holds the exact order
      // `createTrianglesOnNewVertices` walks. Dumping it lets the port's
      // `heap_sort_by_xy_sum` be diffed index by index instead of only through
      // the resulting triangles.
      std::printf("OCCTSORT f=%d n=%d", f, aIndices2.Length());
      for (int anIdx = aIndices2.Lower(); anIdx <= aIndices2.Upper(); ++anIdx)
      {
        std::printf(" %d", aIndices2(anIdx));
      }
      std::printf("\n");
      std::fflush(stdout);

      std::cout << "OCCT   f=" << f << " nodes=" << aStruct->NbNodes()
                << " links=" << aStruct->NbLinks()
                << " domain=" << aStruct->ElementsOfDomain().Extent()
                << " phase1=" << aOrder1Ref.size() << " phase2=" << aOrder2Ref.size()
                << std::endl;

      std::unordered_map<int, int> anInv;
      for (int v = 1; v <= static_cast<int>(aNodes.size()); ++v)
      {
        if (aMap[v] > 0)
        {
          anInv[aMap[v]] = v;
        }
      }
      std::vector<std::array<int, 3>> aTriangles;
      for (IMeshData::IteratorOfMapOfInteger aTriIt(aStruct->ElementsOfDomain());
           aTriIt.More();
           aTriIt.Next())
      {
        int aNode3[3];
        aStruct->ElementNodes(aStruct->GetElement(aTriIt.Key()), aNode3);
        std::array<int, 3> aTri = {anInv[aNode3[0]], anInv[aNode3[1]], anInv[aNode3[2]]};
        std::sort(aTri.begin(), aTri.end());
        aTriangles.push_back(aTri);
      }
      std::sort(aTriangles.begin(), aTriangles.end());
      char aTriName[512];
      std::snprintf(aTriName,
                    sizeof(aTriName),
                    "%s/delaun_tris_occt_f%d.txt",
                    theDir.c_str(),
                    f);
      std::ofstream aOut(aTriName);
      aOut << "TRIS " << aTriangles.size() << "\n";
      for (const std::array<int, 3>& aTri : aTriangles)
      {
        aOut << "T " << aTri[0] << " " << aTri[1] << " " << aTri[2] << "\n";
      }
      continue;
    }

    int              aIdxMismatch = 0;
    for (size_t i = 0; i < aNodes.size(); ++i)
    {
      const BRepMesh_Vertex aVertex(gp_XY(aNodes[i][0], aNodes[i][1]),
                                    static_cast<int>(i),
                                    static_cast<BRepMesh_DegreeOfFreedom>(static_cast<int>(aNodes[i][2])));
      const int aId = aStruct->AddNode(aVertex);
      aMap[i + 1]     = aId;
      if (aId != static_cast<int>(i) + 1)
      {
        ++aIdxMismatch;
      }
    }
    for (const std::array<int, 3>& aL : aLinks)
    {
      aStruct->AddLink(BRepMesh_Edge(aMap[aL[0]],
                                     aMap[aL[1]],
                                     static_cast<BRepMesh_DegreeOfFreedom>(aL[2])));
    }
    const int aLinksDumped = static_cast<int>(aLinks.size());
    const int aLinksInStruct = aStruct->NbLinks();

    IMeshData::VectorOfInteger aIndices;
    for (int v : aOrder)
    {
      aIndices.Append(aMap[v]);
    }

    BRepMesh_Delaun aDelaun(aStruct, aIndices, aCellsU, aCellsV);
    const occ::handle<BRepMesh_DataStructureOfDelaun>& aSt = aDelaun.Result();

    int aFrontier = 0, aFixed = 0, aFree = 0, aDeleted = 0;
    for (int e = 1; e <= aSt->NbLinks(); ++e)
    {
      switch (aSt->GetLink(e).Movability())
      {
        case BRepMesh_Frontier: ++aFrontier; break;
        case BRepMesh_Fixed:    ++aFixed;    break;
        case BRepMesh_Free:     ++aFree;     break;
        default:                ++aDeleted;  break;
      }
    }

    std::cout << "OCCT   f=" << f << " nodes=" << aSt->NbNodes() << " links=" << aSt->NbLinks()
              << " domain=" << aSt->ElementsOfDomain().Extent() << " frontier=" << aFrontier
              << " fixed=" << aFixed << " free=" << aFree << " deleted=" << aDeleted
              << " idx_mismatch=" << aIdxMismatch << " links_dumped=" << aLinksDumped
              << " links_in_struct=" << aLinksInStruct << std::endl;
  }
  return 0;
}

// --- TEMP T-99 (`--delaunsteps <dir> <face>`) -------------------------------------
// Step-by-step variant of `--delaunstruct`: same two-phase replay, but the phase-2
// nodes are fed to `BRepMesh_Delaun::AddVertices` ONE AT A TIME (each call is the
// exact `AddVertices({v})` control flow: a one-element heap sort is a no-op, the
// insertion loop runs for that one vertex, then `ProcessConstraints()`), and the
// domain triangle set is written after every step. Diffed against the port's
// `DELONESTEP` dump, the first divergent insertion index identifies the step whose
// cavity/selection differs, instead of only seeing the end state.
static int runDelaunSteps(const std::string& theDir, const int theFace)
{
  char aName[512];
  std::snprintf(aName, sizeof(aName), "%s/delaun_in_f%d.txt", theDir.c_str(), theFace);
  std::ifstream aIn(aName);
  if (!aIn)
  {
    std::cout << "DELAUNSTEPS f=" << theFace << " (missing " << aName << ")" << std::endl;
    return 1;
  }

  std::vector<std::array<double, 3>> aNodes; // x, y, movability
  std::vector<std::array<int, 3>>    aLinks; // first, last, movability
  std::vector<int>                   aOrder1;
  std::vector<int>                   aOrder2;
  double                             aTolU = 0.0, aTolV = 0.0;
  int                                aCellsU = -1, aCellsV = -1;
  std::string                        aLine;
  while (std::getline(aIn, aLine))
  {
    if (aLine.empty() || aLine[0] == '#')
    {
      continue;
    }
    std::istringstream aSS(aLine);
    std::string        aTag;
    aSS >> aTag;
    if (aTag == "TOL")
    {
      aSS >> aTolU >> aTolV;
    }
    else if (aTag == "CELLS")
    {
      aSS >> aCellsU >> aCellsV;
    }
    else if (aTag == "N")
    {
      int i, m;
      double x, y;
      aSS >> i >> x >> y >> m;
      aNodes.push_back({x, y, static_cast<double>(m)});
    }
    else if (aTag == "L")
    {
      int i, a, b, m;
      aSS >> i >> a >> b >> m;
      aLinks.push_back({a, b, m});
    }
    else if (aTag == "V1")
    {
      int v;
      while (aSS >> v)
      {
        aOrder1.push_back(v);
      }
    }
    else if (aTag == "V2")
    {
      int v;
      while (aSS >> v)
      {
        aOrder2.push_back(v);
      }
    }
  }

  if (aNodes.empty() || aOrder2.empty())
  {
    std::cout << "DELAUNSTEPS f=" << theFace << " (no nodes / no phase2)" << std::endl;
    return 1;
  }

  occ::handle<NCollection_IncAllocator> anAlloc =
    new NCollection_IncAllocator(IMeshData::MEMORY_BLOCK_SIZE_HUGE);
  occ::handle<BRepMesh_DataStructureOfDelaun> aStruct =
    new BRepMesh_DataStructureOfDelaun(anAlloc, static_cast<int>(aNodes.size()));
  aStruct->Data()->SetTolerance(aTolU, aTolV);
  aStruct->Data()->SetCellSize(14.0 * aTolU, 14.0 * aTolV);

  std::vector<int> aMap(aNodes.size() + 1, 0);
  for (int v : aOrder1)
  {
    const BRepMesh_Vertex aVertex(gp_XY(aNodes[v - 1][0], aNodes[v - 1][1]),
                                  v - 1,
                                  static_cast<BRepMesh_DegreeOfFreedom>(
                                    static_cast<int>(aNodes[v - 1][2])));
    aMap[v] = aStruct->AddNode(aVertex);
  }
  for (const std::array<int, 3>& aL : aLinks)
  {
    if (aL[2] == static_cast<int>(BRepMesh_Deleted) || aMap[aL[0]] == 0 || aMap[aL[1]] == 0)
    {
      continue;
    }
    aStruct->AddLink(BRepMesh_Edge(aMap[aL[0]],
                                   aMap[aL[1]],
                                   static_cast<BRepMesh_DegreeOfFreedom>(aL[2])));
  }

  IMeshData::VectorOfInteger aIndices1;
  for (int v : aOrder1)
  {
    aIndices1.Append(aMap[v]);
  }
  BRepMesh_Delaun aDelaun(aStruct, aIndices1, aCellsU, aCellsV);
  std::cout << "DELAUNSTEPS f=" << theFace << " base nodes=" << aStruct->NbNodes()
            << " links=" << aStruct->NbLinks()
            << " domain=" << aStruct->ElementsOfDomain().Extent()
            << " phase2=" << aOrder2.size() << std::endl;

  for (int v : aOrder2)
  {
    const BRepMesh_Vertex aVertex(gp_XY(aNodes[v - 1][0], aNodes[v - 1][1]),
                                  v - 1,
                                  static_cast<BRepMesh_DegreeOfFreedom>(
                                    static_cast<int>(aNodes[v - 1][2])));
    aMap[v] = aStruct->AddNode(aVertex);
  }

  std::unordered_map<int, int> anInv;
  for (int v = 1; v <= static_cast<int>(aNodes.size()); ++v)
  {
    if (aMap[v] > 0)
    {
      anInv[aMap[v]] = v;
    }
  }

  char aOutName[512];
  std::snprintf(aOutName,
                sizeof(aOutName),
                "%s/delaun_steps_occt_f%d.txt",
                theDir.c_str(),
                theFace);
  std::ofstream aOut(aOutName);
  if (!aOut)
  {
    std::cout << "DELAUNSTEPS f=" << theFace << " (cannot write " << aOutName << ")" << std::endl;
    return 1;
  }

  // TEMP T-119: the batched `BRepMesh_Delaun::AddVertices` (`cxx:2415-2417`)
  // heap sorts the index vector before `createTrianglesOnNewVertices`, so the
  // real path walks the sorted order. `OCCT_DELAUNSTEPS_SORT=1` replays that
  // same order one vertex at a time, to be diffed against the port's
  // `DELONESTEP` dump of its own batched order.
  std::vector<int> aOrder(aOrder2);
  if (std::getenv("OCCT_DELAUNSTEPS_SORT") != nullptr)
  {
    auto aLess = [&aNodes](int theLeft, int theRight) {
      return aNodes[theLeft - 1][0] + aNodes[theLeft - 1][1]
             < aNodes[theRight - 1][0] + aNodes[theRight - 1][1];
    };
    std::make_heap(aOrder.begin(), aOrder.end(), aLess);
    std::sort_heap(aOrder.begin(), aOrder.end(), aLess);
  }

  for (size_t k = 0; k < aOrder.size(); ++k)
  {
    IMeshData::VectorOfInteger aOne;
    aOne.Append(aMap[aOrder[k]]);
    aDelaun.AddVertices(aOne);

    std::vector<std::array<int, 3>> aTriangles;
    for (IMeshData::IteratorOfMapOfInteger aTriIt(aStruct->ElementsOfDomain());
         aTriIt.More();
         aTriIt.Next())
    {
      int aNode3[3];
      aStruct->ElementNodes(aStruct->GetElement(aTriIt.Key()), aNode3);
      std::array<int, 3> aTri = {anInv[aNode3[0]], anInv[aNode3[1]], anInv[aNode3[2]]};
      std::sort(aTri.begin(), aTri.end());
      aTriangles.push_back(aTri);
    }
    std::sort(aTriangles.begin(), aTriangles.end());
    aOut << "STEP " << k << "\n";
    for (const std::array<int, 3>& aTri : aTriangles)
    {
      aOut << "T " << aTri[0] << " " << aTri[1] << " " << aTri[2] << "\n";
    }
  }

  std::cout << "DELAUNSTEPS f=" << theFace << " steps=" << aOrder2.size() << " written" << std::endl;
  return 0;
}

// --- TEMP T-121 (`--delaunwhy <dir> <face> [step]`) -------------------------------
// Per-insertion verdict dump: for the requested step (all steps when <step> < 0) it
// prints the exact inputs of ONE `createTrianglesOnNewVertices` iteration
// (`BRepMesh_Delaun.cxx:629-700`):
//   * the vertex that is being inserted (structure id + the port's own node id),
//   * `myCircles.Select(coord)` - the candidate triangle list, in list order,
//   * `Contains()` per candidate (`BRepMesh_Delaun.cxx:2578-2637`) together with the
//     `theEdgeOn` link id and that link's movability,
//   * the triangle the caller's own skip rule (`cxx:653-660`) would select.
// The port prints the same lines under `DCWHY`, so the first field that differs
// names the missing `.cxx` arm instead of only showing the divergent triangle set.
static int runDelaunWhy(const std::string& theDir, const int theFace, const int theStep)
{
  char aName[512];
  std::snprintf(aName, sizeof(aName), "%s/delaun_in_f%d.txt", theDir.c_str(), theFace);
  std::ifstream aIn(aName);
  if (!aIn)
  {
    std::cout << "DELAUNWHY f=" << theFace << " (missing " << aName << ")" << std::endl;
    return 1;
  }

  std::vector<std::array<double, 3>> aNodes; // x, y, movability
  std::vector<std::array<int, 3>>    aLinks; // first, last, movability
  std::vector<int>                   aOrder1;
  std::vector<int>                   aOrder2;
  double                             aTolU = 0.0, aTolV = 0.0;
  int                                aCellsU = -1, aCellsV = -1;
  std::string                        aLine;
  while (std::getline(aIn, aLine))
  {
    if (aLine.empty() || aLine[0] == '#')
    {
      continue;
    }
    std::istringstream aSS(aLine);
    std::string        aTag;
    aSS >> aTag;
    if (aTag == "TOL")
    {
      aSS >> aTolU >> aTolV;
    }
    else if (aTag == "CELLS")
    {
      aSS >> aCellsU >> aCellsV;
    }
    else if (aTag == "N")
    {
      int i, m;
      double x, y;
      aSS >> i >> x >> y >> m;
      aNodes.push_back({x, y, static_cast<double>(m)});
    }
    else if (aTag == "L")
    {
      int i, a, b, m;
      aSS >> i >> a >> b >> m;
      aLinks.push_back({a, b, m});
    }
    else if (aTag == "V1")
    {
      int v;
      while (aSS >> v)
      {
        aOrder1.push_back(v);
      }
    }
    else if (aTag == "V2")
    {
      int v;
      while (aSS >> v)
      {
        aOrder2.push_back(v);
      }
    }
  }

  if (aNodes.empty() || aOrder2.empty())
  {
    std::cout << "DELAUNWHY f=" << theFace << " (no nodes / no phase2)" << std::endl;
    return 1;
  }

  occ::handle<NCollection_IncAllocator> anAlloc =
    new NCollection_IncAllocator(IMeshData::MEMORY_BLOCK_SIZE_HUGE);
  occ::handle<BRepMesh_DataStructureOfDelaun> aStruct =
    new BRepMesh_DataStructureOfDelaun(anAlloc, static_cast<int>(aNodes.size()));
  aStruct->Data()->SetTolerance(aTolU, aTolV);
  aStruct->Data()->SetCellSize(14.0 * aTolU, 14.0 * aTolV);

  std::vector<int> aMap(aNodes.size() + 1, 0);
  for (int v : aOrder1)
  {
    const BRepMesh_Vertex aVertex(gp_XY(aNodes[v - 1][0], aNodes[v - 1][1]),
                                  v - 1,
                                  static_cast<BRepMesh_DegreeOfFreedom>(
                                    static_cast<int>(aNodes[v - 1][2])));
    aMap[v] = aStruct->AddNode(aVertex);
  }
  for (const std::array<int, 3>& aL : aLinks)
  {
    if (aL[2] == static_cast<int>(BRepMesh_Deleted) || aMap[aL[0]] == 0 || aMap[aL[1]] == 0)
    {
      continue;
    }
    aStruct->AddLink(BRepMesh_Edge(aMap[aL[0]],
                                   aMap[aL[1]],
                                   static_cast<BRepMesh_DegreeOfFreedom>(aL[2])));
  }

  IMeshData::VectorOfInteger aIndices1;
  for (int v : aOrder1)
  {
    aIndices1.Append(aMap[v]);
  }
  BRepMesh_Delaun aDelaun(aStruct, aIndices1, aCellsU, aCellsV);
  std::cout << "DELAUNWHY f=" << theFace << " base nodes=" << aStruct->NbNodes()
            << " links=" << aStruct->NbLinks()
            << " domain=" << aStruct->ElementsOfDomain().Extent()
            << " phase2=" << aOrder2.size() << " step=" << theStep << std::endl;

  for (int v : aOrder2)
  {
    const BRepMesh_Vertex aVertex(gp_XY(aNodes[v - 1][0], aNodes[v - 1][1]),
                                  v - 1,
                                  static_cast<BRepMesh_DegreeOfFreedom>(
                                    static_cast<int>(aNodes[v - 1][2])));
    aMap[v] = aStruct->AddNode(aVertex);
  }

  // The batched path heap sorts first (`AddVertices`, `cxx:2415-2417`).
  std::vector<int> aOrder(aOrder2);
  auto             aLess = [&aNodes](int theLeft, int theRight) {
    return aNodes[theLeft - 1][0] + aNodes[theLeft - 1][1] < aNodes[theRight - 1][0] + aNodes[theRight - 1][1];
  };
  std::make_heap(aOrder.begin(), aOrder.end(), aLess);
  std::sort_heap(aOrder.begin(), aOrder.end(), aLess);

  const double aSqTol = aTolU * aTolU + aTolV * aTolV;

  for (size_t k = 0; k < aOrder.size(); ++k)
  {
    const int              aPortV = aOrder[k];
    const int              aIdx   = aMap[aPortV];
    const BRepMesh_Vertex& aV     = aDelaun.GetVertex(aIdx);
    if (theStep < 0 || static_cast<int>(k) == theStep)
    {
      std::printf("WHY step=%d port_v=%d struct=%d uv=%.17e %.17e mov=%d sqtol=%.17e ntri=%d",
                  static_cast<int>(k),
                  aPortV,
                  aIdx,
                  aV.Coord().X(),
                  aV.Coord().Y(),
                  static_cast<int>(aV.Movability()),
                  aSqTol,
                  aStruct->ElementsOfDomain().Extent());
      // `Circles()` is const (the mesh exposes the tool read-only) while `Select`
      // is not: it only appends to the inspector's own candidate list, so casting
      // the (non-const) `aDelaun`'s tool is the same call OCCT itself makes at
      // `BRepMesh_Delaun.cxx:642`.
      BRepMesh_CircleTool&      aTool  = const_cast<BRepMesh_CircleTool&>(aDelaun.Circles());
      IMeshData::ListOfInteger& aCands = aTool.Select(aV.Coord());
      std::printf(" ncand=%d\n", aCands.Extent());
      for (IMeshData::ListOfInteger::Iterator aCIt(aCands); aCIt.More(); aCIt.Next())
      {
        std::printf("  CAND t=%d\n", aCIt.Value());
      }
      int aChosen = 0;
      for (IMeshData::ListOfInteger::Iterator aCIt(aCands); aCIt.More(); aCIt.Next())
      {
        int  anEdgeOn = 0;
        bool aCont    = aDelaun.Contains(aCIt.Value(), aV, aSqTol, anEdgeOn);
        int  anEdgeMov = anEdgeOn != 0 ? static_cast<int>(aDelaun.GetEdge(anEdgeOn).Movability()) : -1;
        std::printf("  CTEST t=%d contains=%d edge=%d edgemov=%d\n",
                    aCIt.Value(),
                    aCont ? 1 : 0,
                    anEdgeOn,
                    anEdgeMov);
        if (aChosen != 0 || !aCont)
        {
          continue;
        }
        if (anEdgeOn != 0 && aDelaun.GetEdge(anEdgeOn).Movability() != BRepMesh_Free)
        {
          if (aV.Movability() == BRepMesh_Free)
          {
            continue;
          }
        }
        aChosen = aCIt.Value();
      }
      std::printf("  CHOSEN t=%d\n", aChosen);
    }

    IMeshData::VectorOfInteger aOne;
    aOne.Append(aIdx);
    aDelaun.AddVertices(aOne);
    if (theStep < 0 || static_cast<int>(k) == theStep)
    {
      std::printf("  AFTER tris=%d\n", aStruct->ElementsOfDomain().Extent());
    }
  }

  std::cout << "DELAUNWHY f=" << theFace << " done tris=" << aStruct->ElementsOfDomain().Extent()
            << std::endl;
  return 0;
}

// --- TEMP T-99 (`--boundary <dir>`) ----------------------------------------------
// Dumps, for the faces of interest, the pcurve points OCCT's real pipeline feeds
// into `BRepMesh_BaseMeshAlgo::initDataStructure` (wire -> edge -> pcurve ->
// `GetPoint(0..ParametersNb()-1)`, i.e. exactly the loop at
// `BRepMesh_BaseMeshAlgo.cxx:87-125`). The model is built through
// `BRepMesh_Context` in `IMeshTools_MeshBuilder`'s order (BuildModel ->
// DiscretizeEdges -> HealModel -> PreProcessModel) so the points are the ones the
// real mesher sees, and each face's bbox is printed so the port's face can be
// paired by geometry rather than by index (D21).
static int runBoundary(const TopoDS_Shape& theShape, const std::string& theDir)
{
  Bnd_Box aShapeBox;
  BRepBndLib::Add(theShape, aShapeBox, false);
  double aBox[6] = {0, 0, 0, 0, 0, 0};
  if (!aShapeBox.IsVoid())
  {
    aShapeBox.Get(aBox[0], aBox[1], aBox[2], aBox[3], aBox[4], aBox[5]);
  }
  const double aMaxComp =
    std::max(std::max(aBox[3] - aBox[0], aBox[4] - aBox[1]), aBox[5] - aBox[2]);

  IMeshTools_Parameters aParams;
  // The port's a3n00 gate meshes with `prs3d_get_deflection(shape, 0.1)`
  // (`maxComp * 0.001 * 4`) and 20 degrees, so match that exactly.
  aParams.Deflection    = aMaxComp * 0.001 * 4.0;
  aParams.Angle         = 20.0 * M_PI / 180.0;
  aParams.InParallel    = false;
  aParams.Relative      = false;
  // Mirror `BRepMesh_IncrementalMesh::initParameters` (`hxx:91-96`): the real
  // mesher runs with `MinSize = max(RelMinSize * min(defl, deflInterior), Confusion)`.
  // The earlier `Precision::Confusion()` made `BRepMesh_CurveTessellator::init` hand
  // `GCPnts_TangentialDeflection` a `MinimumLength` of 1e-7, so every arc got the full
  // 10-degree (half-angle) step and the dump was NOT the structure the mesher builds.
  aParams.MinSize       = (std::max)(IMeshTools_Parameters::RelMinSize() * aParams.Deflection,
                                     Precision::Confusion());
  aParams.AdjustMinSize = false;

  occ::handle<BRepMesh_Context> aCtx = new BRepMesh_Context;
  aCtx->SetShape(theShape);
  aCtx->ChangeParameters()            = aParams;
  aCtx->ChangeParameters().CleanModel = false;
  const bool aBuilt     = aCtx->BuildModel();
  const bool aDiscret   = aBuilt && aCtx->DiscretizeEdges();
  const bool aHealed    = aDiscret && aCtx->HealModel();
  const bool aPreproc   = aHealed && aCtx->PreProcessModel();
  const occ::handle<IMeshData_Model>& aModel = aCtx->GetModel();
  std::cout << "BOUNDARY defl=" << aParams.Deflection << " angle=" << aParams.Angle
            << " built=" << (aBuilt ? 1 : 0) << " discret=" << (aDiscret ? 1 : 0)
            << " healed=" << (aHealed ? 1 : 0) << " preproc=" << (aPreproc ? 1 : 0)
            << " model_null=" << (aModel.IsNull() ? 1 : 0) << std::endl;
  if (aModel.IsNull())
  {
    return 0;
  }

  const int aFaces[] = {0,   208, 209, 169, 174, 189, 192, 193, 194,
                        195, 196, 197, 198, 199, 204, 206, 210, 212};
  (void)aFaces;

  std::ostringstream anOut;
  for (int f = 0; f < aModel->FacesNb(); ++f)
  {
    const IMeshData::IFaceHandle& aDFace = aModel->GetFace(f);
    Bnd_Box                       aFaceBox;
    BRepBndLib::Add(aDFace->GetFace(), aFaceBox, false);
    double v[6] = {0, 0, 0, 0, 0, 0};
    if (!aFaceBox.IsVoid())
    {
      aFaceBox.Get(v[0], v[1], v[2], v[3], v[4], v[5]);
    }
    int aPts = 0;
    for (int w = 0; w < aDFace->WiresNb(); ++w)
    {
      const IMeshData::IWireHandle& aDWire = aDFace->GetWire(w);
      for (int e = 0; e < aDWire->EdgesNb(); ++e)
      {
        const IMeshData::ListOfInteger& aList = aDWire->GetEdge(e)->GetPCurves(aDFace.get());
        for (IMeshData::ListOfInteger::Iterator it(aList); it.More(); it.Next())
        {
          aPts += aDWire->GetEdge(e)->GetPCurve(it.Value())->ParametersNb();
        }
      }
    }
    anOut << std::fixed << "BOUND f=" << f << " type=" << (int)aDFace->GetSurface()->GetType()
          << " wires=" << aDFace->WiresNb() << " pts=" << aPts << " bbox=(" << v[0] << "," << v[1]
          << "," << v[2] << ")-(" << v[3] << "," << v[4] << "," << v[5] << ")\n";
    for (int w = 0; w < aDFace->WiresNb(); ++w)
    {
      const IMeshData::IWireHandle& aDWire = aDFace->GetWire(w);
      anOut << "W " << w << " edges=" << aDWire->EdgesNb()
            << " selfint=" << (aDWire->IsSet(IMeshData_SelfIntersectingWire) ? 1 : 0)
            << " open=" << (aDWire->IsSet(IMeshData_OpenWire) ? 1 : 0) << "\n";
      for (int e = 0; e < aDWire->EdgesNb(); ++e)
      {
        const IMeshData::IEdgeHandle& aDEdge = aDWire->GetEdge(e);
        if (true)
        {
          // `BRepMesh_CurveTessellator::init` inputs (`cxx:60-120`): the edge's stored
          // deflections and the curve type decide how many points
          // `GCPnts_TangentialDeflection` emits, so dump them next to the points.
          BRepAdaptor_Curve aBRepCurve(aDEdge->GetEdge());
          anOut << "  EINFO w=" << w << " e=" << e << " ctype=" << (int)aBRepCurve.GetType()
                << " defl=" << aDEdge->GetDeflection()
                << " ang=" << aDEdge->GetAngularDeflection()
                << " sameparam=" << (aDEdge->GetSameParam() ? 1 : 0)
                << " samerange=" << (aDEdge->GetSameRange() ? 1 : 0)
                << " degen=" << (aDEdge->GetDegenerated() ? 1 : 0)
                << " free=" << (aDEdge->IsFree() ? 1 : 0)
                << " first=" << aBRepCurve.FirstParameter()
                << " last=" << aBRepCurve.LastParameter();
          if (aBRepCurve.GetType() == GeomAbs_Circle)
          {
            anOut << " radius=" << aBRepCurve.Circle().Radius();
          }
          anOut << "\n";
        }
        // Exactly `BRepMesh_BaseMeshAlgo.cxx:91-96`: the pcurves of this edge on
        // this face, in the structure's own order.
        const IMeshData::ListOfInteger& aListOfPCurves = aDEdge->GetPCurves(aDFace.get());
        for (IMeshData::ListOfInteger::Iterator aPCurveIt(aListOfPCurves); aPCurveIt.More();
             aPCurveIt.Next())
        {
          const IMeshData::IPCurveHandle& aPCurve = aDEdge->GetPCurve(aPCurveIt.Value());
          anOut << "E w=" << w << " e=" << e << " pcId=" << aPCurveIt.Value()
                << " ori=" << static_cast<int>(aPCurve->GetOrientation())
                << " fwd=" << (aPCurve->IsForward() ? 1 : 0)
                << " n=" << aPCurve->ParametersNb() << "\n";
          for (int p = 0; p < aPCurve->ParametersNb(); ++p)
          {
            const gp_Pnt2d& aP2 = aPCurve->GetPoint(p);
            char            aBuf[128];
            std::snprintf(aBuf, sizeof(aBuf), "P %.17e %.17e\n", aP2.X(), aP2.Y());
            anOut << aBuf;
          }
        }
      }
    }
  }

  const std::string aPath = theDir + "/boundary.txt";
  std::ofstream     aFile(aPath.c_str());
  if (!aFile)
  {
    std::cout << "BOUNDARY cannot write " << aPath << std::endl;
    return 0;
  }
  aFile << anOut.str();
  aFile.close();
  std::cout << "BOUNDARY wrote " << aPath << std::endl;
  return 0;
}

// --- TEMP T-104 (`--structdump <faceIdx>`) -----------------------------------------
// Prints the exact `BRepMesh_DataStructureOfDelaun` OCCT's own node-insertion algo
// hands to `BRepMesh_Delaun` for one discrete face, in the port dump's `N`/`L` line
// format (`crates/occt-topo/src/meshing/node_insertion.rs::dump_delaun_structure`),
// so the two files diff line by line. The model is built exactly like `--boundary`
// (same deflection/angle as the port's `prs3d_get_deflection(shape, 0.1)`), and the
// face index is the `IMeshData_Model` order - pair it with the port's face by the
// printed bbox, not by index.
template <class Splitter>
class ZZStructDumpAlgo
    : public BRepMesh_DelaunayNodeInsertionMeshAlgo<Splitter, BRepMesh_DelaunayBaseMeshAlgo>
{protected:
  // TEMP T-105: dumps, in the port's own `N`/`L`/`TOL`/`CELLS`/`V` format, the exact
  // structure the *real* algo passes to `BRepMesh_Delaun` (the moment `generateMesh`
  // is entered). Diffing this against the port's `delaun_in_f*.txt` says whether the
  // two sides really start from the same input.
  void generateMesh(const Message_ProgressRange& theRange) override
  {
    typedef BRepMesh_DelaunayNodeInsertionMeshAlgo<Splitter, BRepMesh_DelaunayBaseMeshAlgo> Base;
    const char* aPath = std::getenv("ZZDELAUNIN");
    if (aPath != nullptr)
    {
      const occ::handle<BRepMesh_DataStructureOfDelaun>& aSt = this->getStructure();
      const int aNb = this->getNodesMap()->Length();
      double    aTolX = 0.0, aTolY = 0.0;
      aSt->Data()->GetTolerance(aTolX, aTolY);
      const std::pair<int, int> aCells = this->getCellsCount(aNb);
      std::ofstream aFile(aPath);
      aFile << "TOL " << std::setprecision(18) << aTolX << " " << aTolY << "\n";
      aFile << "CELLS " << aCells.first << " " << aCells.second << "\n";
      for (int i = 1; i <= aSt->NbNodes(); ++i)
      {
        const BRepMesh_Vertex& aV = aSt->GetNode(i);
        aFile << "N " << i << " " << std::setprecision(18) << aV.Coord().X() << " "
              << aV.Coord().Y() << " " << static_cast<int>(aV.Movability()) << "\n";
      }
      for (int i = 1; i <= aSt->NbLinks(); ++i)
      {
        const BRepMesh_Edge& aL = aSt->GetLink(i);
        aFile << "L " << i << " " << aL.FirstNode() << " " << aL.LastNode() << " "
              << static_cast<int>(aL.Movability()) << "\n";
      }
      aFile << "V";
      for (int i = 1; i <= aNb; ++i)
      {
        aFile << " " << i;
      }
      aFile << "\n";
      aFile.close();
      std::printf("OCCTPOST DELAUNIN wrote %s nodes=%d links=%d vorder=%d tol=(%.17e,%.17e) "
                  "cells=(%d,%d)\n",
                  aPath,
                  aSt->NbNodes(),
                  aSt->NbLinks(),
                  aNb,
                  aTolX,
                  aTolY,
                  aCells.first,
                  aCells.second);
    }
    Base::generateMesh(theRange);
  }

  // TEMP T-105: the real factory algo (`BRepMesh_MeshAlgoFactory`, torus ->
  // `NodeInsertionMeshAlgo<TorusRangeSplitter>`) is a template, so the only way to
  // observe the state it hands `BRepMesh_Delaun` and the state it leaves behind is to
  // subclass it. `postProcessMesh` is the seam between the base Delaunay
  // (`generateMesh`) and the surface-node insertion (`insertNodes` -> `AddVertices`),
  // which is exactly where the port's torus faces lose every triangle.
  void postProcessMesh(BRepMesh_Delaun& theMesher, const Message_ProgressRange& theRange) override
  {
    typedef BRepMesh_DelaunayNodeInsertionMeshAlgo<Splitter, BRepMesh_DelaunayBaseMeshAlgo> Base;
    dumpAlgoState("POST_BASEMESH", theMesher);
    {
      double aTolX = 0.0, aTolY = 0.0;
      const occ::handle<BRepMesh_DataStructureOfDelaun>& aSt = this->getStructure();
      aSt->Data()->GetTolerance(aTolX, aTolY);
      const std::pair<int, int>          aCells  = this->getCellsCount(aSt->NbNodes() - 3);
      const std::pair<double, double>&   aDelta  = this->getRangeSplitter().GetDelta();
      const std::pair<double, double>&   aTolUV  = this->getRangeSplitter().GetToleranceUV();
      const std::pair<double, double>&   aRangeU = this->getRangeSplitter().GetRangeU();
      const std::pair<double, double>&   aRangeV = this->getRangeSplitter().GetRangeV();
      std::printf("OCCTPOST PARAMS defl=%.9f tol=(%.17e,%.17e) cells=(%d,%d) delta=(%.17e,%.17e) "
                  "toluv=(%.17e,%.17e) rangeU=(%.9f,%.9f) rangeV=(%.9f,%.9f)\n",
                  this->getDFace()->GetDeflection(),
                  aTolX,
                  aTolY,
                  aCells.first,
                  aCells.second,
                  aDelta.first,
                  aDelta.second,
                  aTolUV.first,
                  aTolUV.second,
                  aRangeU.first,
                  aRangeU.second,
                  aRangeV.first,
                  aRangeV.second);
    }
    {
      const Handle(IMeshData::ListOfPnt2d) aSurfaceNodes =
        this->getRangeSplitter().GenerateSurfaceNodes(this->getParameters());
      int aIn = 0;
      const int aTotal = aSurfaceNodes.IsNull() ? 0 : aSurfaceNodes->Size();
      if (!aSurfaceNodes.IsNull())
      {
        IMeshData::ListOfPnt2d::Iterator aNodesIt(*aSurfaceNodes);
        for (; aNodesIt.More(); aNodesIt.Next())
        {
          if (this->getClassifier()->Perform(aNodesIt.Value()) == TopAbs_IN)
          {
            ++aIn;
          }
        }
      }
      std::printf("OCCTPOST surface_nodes total=%d in=%d\n", aTotal, aIn);
    }
    Base::postProcessMesh(theMesher, theRange);
    dumpAlgoState("POST_INSERT", theMesher);
  }

  void dumpAlgoState(const char* theTag, BRepMesh_Delaun& theMesher) const
  {
    const occ::handle<BRepMesh_DataStructureOfDelaun>& aSt = this->getStructure();
    const occ::handle<BRepMesh_DataStructureOfDelaun>& aRes = theMesher.Result();
    int aFrontier = 0, aFixed = 0, aFree = 0, aOther = 0;
    for (int e = 1; e <= aRes->NbLinks(); ++e)
    {
      switch (aRes->GetLink(e).Movability())
      {
        case BRepMesh_Frontier: ++aFrontier; break;
        case BRepMesh_Fixed:    ++aFixed;    break;
        case BRepMesh_Free:     ++aFree;     break;
        default:                ++aOther;    break;
      }
    }
    std::printf("OCCTPOST %s struct_nodes=%d struct_links=%d nodes=%d links=%d domain=%d "
                "frontier=%d fixed=%d free=%d other=%d\n",
                theTag,
                aSt.IsNull() ? -1 : aSt->NbNodes(),
                aSt.IsNull() ? -1 : aSt->NbLinks(),
                aRes->NbNodes(),
                aRes->NbLinks(),
                aRes->ElementsOfDomain().Extent(),
                aFrontier,
                aFixed,
                aFree,
                aOther);
    IMeshData::IteratorOfMapOfInteger aTriIt(aRes->ElementsOfDomain());
    for (; aTriIt.More(); aTriIt.Next())
    {
      int aNodes[3];
      aRes->ElementNodes(aRes->GetElement(aTriIt.Key()), aNodes);
      std::printf("OCCTPOST %s T %d %d %d %d\n",
                  theTag,
                  aTriIt.Key(),
                  aNodes[0],
                  aNodes[1],
                  aNodes[2]);
    }
    // TEMP T-105b: the full element sequence, deleted slots included, so the port's
    // `add_element` trace (`DELADDT`) can be diffed call by call.
    for (int i = 1; i <= aRes->NbElements(); ++i)
    {
      const BRepMesh_Triangle& anElt = aRes->GetElement(i);
      int                      aN[3];
      aRes->ElementNodes(anElt, aN);
      std::printf("OCCTEL %s E %d %d %d %d %d\n",
                  theTag,
                  i,
                  aN[0],
                  aN[1],
                  aN[2],
                  static_cast<int>(anElt.Movability()));
    }
  }

  bool initDataStructure() override
  {
    typedef BRepMesh_DelaunayNodeInsertionMeshAlgo<Splitter, BRepMesh_DelaunayBaseMeshAlgo> Base;
    const bool isOk = Base::initDataStructure();

    const occ::handle<BRepMesh_DataStructureOfDelaun>& aSt = this->getStructure();
    if (!aSt.IsNull())
    {
      std::printf("OCCTSTRUCT nodes=%d links=%d\n", aSt->NbNodes(), aSt->NbLinks());
      for (int i = 1; i <= aSt->NbNodes(); ++i)
      {
        const BRepMesh_Vertex& aV = aSt->GetNode(i);
        std::printf("N %d %.17e %.17e %d\n",
                    i,
                    aV.Coord().X(),
                    aV.Coord().Y(),
                    static_cast<int>(aV.Movability()));
      }
      for (int i = 1; i <= aSt->NbLinks(); ++i)
      {
        const BRepMesh_Edge& aL = aSt->GetLink(i);
        std::printf("L %d %d %d %d\n",
                    i,
                    aL.FirstNode(),
                    aL.LastNode(),
                    static_cast<int>(aL.Movability()));
      }
    }
    return isOk;
  }
};

// TEMP T-119: `BRepMesh_MeshAlgoFactory::GetAlgo` (`cxx:100-101`) pairs
// `GeomAbs_SurfaceOfRevolution` with `DeflectionControlMeshAlgo<
// BRepMesh_BoundaryParamsRangeSplitter>` (`cxx:42-46`). The DC pass is where that
// face's node count is decided, so hook `optimizeMesh` to dump the structure it
// receives (`POST_INSERT`, = base surface-node insertion) and the one it leaves
// (`POST_DC`, = the mesh `Compute` returns).
template <class Splitter>
class ZZStructDumpDCAlgo
    : public BRepMesh_DelaunayDeflectionControlMeshAlgo<Splitter, BRepMesh_DelaunayBaseMeshAlgo>
{
  typedef BRepMesh_DelaunayDeflectionControlMeshAlgo<Splitter, BRepMesh_DelaunayBaseMeshAlgo> BaseDC;

protected:
  void postProcessMesh(BRepMesh_Delaun& theMesher, const Message_ProgressRange& theRange) override
  {
    // TEMP POST_BASE: state right after the base mesher ran and before any
    // surface node is registered (`BRepMesh_DelaunayNodeInsertionMeshAlgo::
    // postProcessMesh` -> `InsertionBaseClass::postProcessMesh`).
    this->dumpState("POST_BASE", theMesher);
    BaseDC::postProcessMesh(theMesher, theRange);
  }

  void optimizeMesh(BRepMesh_Delaun& theMesher, const Message_ProgressRange& theRange) override
  {
    this->dumpState("POST_INSERT", theMesher);
    // TEMP DCPASS: verbatim copy of `BaseDC::optimizeMesh`'s loop with a per-pass
    // print, so the port's `DELSTAGE after_vertex_loop` chain can be lined up
    // pass by pass against OCCT.
    occ::handle<NCollection_IncAllocator> aTmpAlloc =
      new NCollection_IncAllocator(IMeshData::MEMORY_BLOCK_SIZE_HUGE);

    this->mySqMinSize = this->getParameters().MinSize * this->getParameters().MinSize;
    this->myCouplesMap =
      new IMeshData::MapOfOrientedEdges(3 * this->getStructure()->ElementsOfDomain().Extent(),
                                        aTmpAlloc);
    this->myControlNodes = new IMeshData::ListOfPnt2d(aTmpAlloc);
    this->myCircles      = &theMesher.Circles();

    const int aIterationsNb = 11;
    bool      isInserted    = true;
    Message_ProgressScope aPS(theRange, "Iteration", aIterationsNb);
    for (int aPass = 1; aPass <= aIterationsNb && isInserted && !this->myIsAllDegenerated; ++aPass)
    {
      if (!aPS.More())
      {
        return;
      }
      this->myMaxSqDeflection  = -1.;
      this->myIsAllDegenerated = true;
      this->myControlNodes->Clear();

      const int aDomainBefore = this->getStructure()->ElementsOfDomain().Extent();
      if (aDomainBefore < 1)
      {
        break;
      }
      IMeshData::IteratorOfMapOfInteger aTriangleIt(this->getStructure()->ElementsOfDomain());
      for (; aTriangleIt.More(); aTriangleIt.Next())
      {
        const BRepMesh_Triangle& aTriangle = this->getStructure()->GetElement(aTriangleIt.Key());
        this->splitTriangleGeometry(aTriangle);
      }

      std::printf("OCCTDCPASS pass=%d domain_before=%d ctrl=%d maxsq=%.17g alldeg=%d\n",
                  aPass,
                  aDomainBefore,
                  (int)this->myControlNodes->Size(),
                  this->myMaxSqDeflection,
                  this->myIsAllDegenerated ? 1 : 0);
      std::fflush(stdout);
      // TEMP CTRLDUMP: the exact control node list this pass hands to
      // `insertNodes`, in face-basis coordinates.
      {
        int aIdx = 0;
        for (IMeshData::ListOfPnt2d::Iterator aCtrlIt(*this->myControlNodes); aCtrlIt.More();
             aCtrlIt.Next())
        {
          std::printf("OCCTCTRL P%d %d %.17g %.17g\n",
                      aPass,
                      ++aIdx,
                      aCtrlIt.Value().X(),
                      aCtrlIt.Value().Y());
        }
        std::fflush(stdout);
      }

      isInserted = this->insertNodes(this->myControlNodes, theMesher, aPS.Next());

      std::printf("OCCTDCPASS pass=%d domain_after=%d nodes=%d inserted=%d\n",
                  aPass,
                  this->getStructure()->ElementsOfDomain().Extent(),
                  this->getStructure()->NbNodes(),
                  isInserted ? 1 : 0);
      std::fflush(stdout);
      {
        char aTag[16];
        std::snprintf(aTag, sizeof(aTag), "P%d", aPass);
        this->dumpState(aTag, theMesher);
      }
    }

    this->myCouplesMap.Nullify();
    this->myControlNodes.Nullify();

    if (!(this->myMaxSqDeflection < 0.))
    {
      this->getDFace()->SetDeflection(std::sqrt(this->myMaxSqDeflection));
    }
    this->dumpState("POST_DC", theMesher);
  }

  void dumpState(const char* theTag, BRepMesh_Delaun& theMesher) const
  {
    const occ::handle<BRepMesh_DataStructureOfDelaun>& aSt  = this->getStructure();
    const occ::handle<BRepMesh_DataStructureOfDelaun>& aRes = theMesher.Result();
    int aFrontier = 0, aFixed = 0, aFree = 0, aOther = 0;
    for (int e = 1; e <= aRes->NbLinks(); ++e)
    {
      switch (aRes->GetLink(e).Movability())
      {
        case BRepMesh_Frontier: ++aFrontier; break;
        case BRepMesh_Fixed:    ++aFixed;    break;
        case BRepMesh_Free:     ++aFree;     break;
        default:                ++aOther;    break;
      }
    }
    std::printf("OCCTDC %s struct_nodes=%d struct_links=%d nodes=%d links=%d domain=%d "
                "frontier=%d fixed=%d free=%d other=%d\n",
                theTag,
                aSt.IsNull() ? -1 : aSt->NbNodes(),
                aSt.IsNull() ? -1 : aSt->NbLinks(),
                aRes->NbNodes(),
                aRes->NbLinks(),
                aRes->ElementsOfDomain().Extent(),
                aFrontier,
                aFixed,
                aFree,
                aOther);
    IMeshData::IteratorOfMapOfInteger aTriIt(aRes->ElementsOfDomain());
    for (; aTriIt.More(); aTriIt.Next())
    {
      int aNodes[3];
      aRes->ElementNodes(aRes->GetElement(aTriIt.Key()), aNodes);
        std::printf("OCCTDC %s T %d %d %d %d\n",
                    theTag,
                    aTriIt.Key(),
                    aNodes[0],
                    aNodes[1],
                    aNodes[2]);
    }
    // TEMP T-119: node coordinates of the mesher's structure, so the port's
    // registered node list can be compared entry by entry (same ids).
    for (int i = 1; i <= aRes->NbNodes(); ++i)
    {
      const BRepMesh_Vertex& aV = aRes->GetNode(i);
      std::printf("OCCTDCN %s %d %.17g %.17g %d\n",
                  theTag,
                  i,
                  aV.Coord().X(),
                  aV.Coord().Y(),
                  static_cast<int>(aV.Movability()));
    }
    for (int i = 1; i <= aRes->NbElements(); ++i)
    {
      const BRepMesh_Triangle& anElt = aRes->GetElement(i);
      int                      aN[3];
      aRes->ElementNodes(anElt, aN);
      std::printf("OCCTEL %s E %d %d %d %d %d\n",
                  theTag,
                  i,
                  aN[0],
                  aN[1],
                  aN[2],
                  static_cast<int>(anElt.Movability()));
    }
  }
};

// TEMP T-120 (`--srfnodes <faceIdx>`): prints, for one discrete face, the surface
// nodes the factory-selected splitter produces and the verdict
// `BRepMesh_Classifier::Perform` gives each of them, i.e. exactly the filter
// `BRepMesh_DelaunayNodeInsertionMeshAlgo::insertNodes` applies
// (`BRepMesh_DelaunayNodeInsertionMeshAlgo.hxx:115-123`). Ground truth for the
// port's `list_surface_nodes` (`node_insertion.rs`) drop (Sphere.iges: 323
// generated, 255 IN, 68 ON). The classifier filter is identical in the DC variant,
// so the non-DC node-insertion algo is used here.
template <class Splitter>
class ZZSrfNodeCountAlgo
    : public BRepMesh_DelaunayNodeInsertionMeshAlgo<Splitter, BRepMesh_DelaunayBaseMeshAlgo>
{
  typedef BRepMesh_DelaunayNodeInsertionMeshAlgo<Splitter, BRepMesh_DelaunayBaseMeshAlgo> Base;

protected:
  void postProcessMesh(BRepMesh_Delaun& theMesher, const Message_ProgressRange& theRange) override
  {
    const Handle(IMeshData::ListOfPnt2d) aNodes =
      this->getRangeSplitter().GenerateSurfaceNodes(this->getParameters());
    int aIn = 0, aOut = 0, aOn = 0;
    const int aTotal = aNodes.IsNull() ? 0 : aNodes->Size();
    if (!aNodes.IsNull())
    {
      for (IMeshData::ListOfPnt2d::Iterator anIt(*aNodes); anIt.More(); anIt.Next())
      {
        switch (this->getClassifier()->Perform(anIt.Value()))
        {
          case TopAbs_IN:
            ++aIn;
            break;
          case TopAbs_OUT:
            ++aOut;
            std::printf("OCCTSRF OUT u=%.17g v=%.17g\n", anIt.Value().X(), anIt.Value().Y());
            break;
          default:
            ++aOn;
            break;
        }
      }
    }
    std::printf("OCCTSRF face total=%d in=%d out=%d on=%d stype=%d\n",
                aTotal,
                aIn,
                aOut,
                aOn,
                static_cast<int>(this->getDFace()->GetSurface()->GetType()));
    Base::postProcessMesh(theMesher, theRange);
  }
};

// --- TEMP T-133 (`--lens <faceIdx> <u0> <u1> <v0> <v1>`) ---------------------------
// Replays `BRepMesh_DefaultRangeSplitter::computeLengthU` / `computeLengthV`
// (`BRepMesh_DefaultRangeSplitter.cxx:144-198`) on one discrete face with an
// explicitly given range and prints every `BRepAdaptor_Surface::D0` sample at full
// precision plus the two lengths. The range is passed in because it is already the
// same on both sides (`--srparams` beats the port's NURBSGRID print), so any sample
// difference can only come from the surface evaluator, and the length difference
// from there into `computeDelta` (`cxx:129-137`).
static int runLenSamples(const TopoDS_Shape& theShape, const int theFaceIdx,
                         const double u0, const double u1, const double v0, const double v1)
{
  occ::handle<BRepMesh_Context> aCtx = new BRepMesh_Context;
  aCtx->SetShape(theShape);
  IMeshTools_Parameters aParams;
  aParams.InParallel           = false;
  aParams.Relative             = false;
  aCtx->ChangeParameters()     = aParams;
  aCtx->ChangeParameters().CleanModel = false;
  if (!aCtx->BuildModel())
  {
    std::printf("LENS built=0\n");
    return 0;
  }
  const occ::handle<IMeshData_Model>& aModel = aCtx->GetModel();
  if (aModel.IsNull() || theFaceIdx < 0 || theFaceIdx >= aModel->FacesNb())
  {
    std::printf("LENS bad face\n");
    return 0;
  }
  const occ::handle<BRepAdaptor_Surface>& aSurf = aModel->GetFace(theFaceIdx)->GetSurface();
  std::printf("LENS f=%d type=%d ru=(%.17g,%.17g) rv=(%.17g,%.17g)\n",
              theFaceIdx,
              static_cast<int>(aSurf->GetType()),
              u0,
              u1,
              v0,
              v1);
  {
    double longv = 0.0;
    gp_Pnt P11, P12, P21, P22, P31, P32;
    double dv     = 0.05 * (v1 - v0);
    double dfuave = 0.5 * (u1 + u0);
    double dfvcur;
    aSurf->D0(u0, v0, P11);
    aSurf->D0(dfuave, v0, P21);
    aSurf->D0(u1, v0, P31);
    std::printf("LENS V k=0 (%.17g,%.17g,%.17g) (%.17g,%.17g,%.17g) (%.17g,%.17g,%.17g)\n",
                P11.X(), P11.Y(), P11.Z(), P21.X(), P21.Y(), P21.Z(), P31.X(), P31.Y(), P31.Z());
    double aV = v0 + dv;
    for (int i1 = 1; i1 <= 20; i1++, aV += dv)
    {
      aSurf->D0(u0, aV, P12);
      aSurf->D0(dfuave, aV, P22);
      aSurf->D0(u1, aV, P32);
      longv += (P11.Distance(P12) + P21.Distance(P22) + P31.Distance(P32));
      P11 = P12;
      P21 = P22;
      P31 = P32;
      std::printf("LENS V k=%d (%.17g,%.17g,%.17g) (%.17g,%.17g,%.17g) (%.17g,%.17g,%.17g)\n",
                  i1, P11.X(), P11.Y(), P11.Z(), P21.X(), P21.Y(), P21.Z(), P31.X(), P31.Y(), P31.Z());
    }
    std::printf("LENS V longv=%.17g\n", longv / 3.);
  }
  {
    double longu = 0.0;
    gp_Pnt P11, P12, P21, P22, P31, P32;
    double du     = 0.05 * (u1 - u0);
    double dfvave = 0.5 * (v1 + v0);
    aSurf->D0(u0, v0, P11);
    aSurf->D0(u0, dfvave, P21);
    aSurf->D0(u0, v1, P31);
    std::printf("LENS U k=0 (%.17g,%.17g,%.17g) (%.17g,%.17g,%.17g) (%.17g,%.17g,%.17g)\n",
                P11.X(), P11.Y(), P11.Z(), P21.X(), P21.Y(), P21.Z(), P31.X(), P31.Y(), P31.Z());
    double aU = u0 + du;
    for (int i1 = 1; i1 <= 20; i1++, aU += du)
    {
      aSurf->D0(aU, v0, P12);
      aSurf->D0(aU, dfvave, P22);
      aSurf->D0(aU, v1, P32);
      longu += (P11.Distance(P12) + P21.Distance(P22) + P31.Distance(P32));
      P11 = P12;
      P21 = P22;
      P31 = P32;
      std::printf("LENS U k=%d (%.17g,%.17g,%.17g) (%.17g,%.17g,%.17g) (%.17g,%.17g,%.17g)\n",
                  i1, P11.X(), P11.Y(), P11.Z(), P21.X(), P21.Y(), P21.Z(), P31.X(), P31.Y(), P31.Z());
    }
    std::printf("LENS U longu=%.17g\n", longu / 3.);
  }
  // TEMP T-133: the same two replays through the raw `Geom_Surface` (the port's
  // `Surface::d0`, i.e. `Geom_BSplineSurface::EvalD0` / `BSplSLib::D0`), so the
  // cache arm of `GeomAdaptor_Surface::EvalD0` (`GeomAdaptor_Surface.cxx:1063-1073`
  // -> `BSplSLib_Cache::D0`) can be told apart from the raw evaluator.
  {
    TopLoc_Location                 aLoc;
    const occ::handle<Geom_Surface> aRaw =
      BRep_Tool::Surface(aModel->GetFace(theFaceIdx)->GetFace(), aLoc);
    const bool     isIdentity = aLoc.IsIdentity();
    const gp_Trsf& aTrsf      = aLoc.Transformation();
    auto           aD0        = [&](double u, double v) {
      gp_Pnt aP = aRaw->EvalD0(u, v);
      if (!isIdentity)
      {
        aP.Transform(aTrsf);
      }
      return aP;
    };
    {
      double       longv = 0.0;
      const double dv    = 0.05 * (v1 - v0);
      const double uave  = 0.5 * (u1 + u0);
      gp_Pnt       P11 = aD0(u0, v0), P21 = aD0(uave, v0), P31 = aD0(u1, v0);
      std::printf("LENRAW V k=0 (%.17g,%.17g,%.17g) (%.17g,%.17g,%.17g) (%.17g,%.17g,%.17g)\n",
                  P11.X(), P11.Y(), P11.Z(), P21.X(), P21.Y(), P21.Z(), P31.X(), P31.Y(), P31.Z());
      double aV = v0 + dv;
      for (int i1 = 1; i1 <= 20; i1++, aV += dv)
      {
        gp_Pnt P12 = aD0(u0, aV), P22 = aD0(uave, aV), P32 = aD0(u1, aV);
        longv += (P11.Distance(P12) + P21.Distance(P22) + P31.Distance(P32));
        P11 = P12;
        P21 = P22;
        P31 = P32;
        std::printf("LENRAW V k=%d (%.17g,%.17g,%.17g) (%.17g,%.17g,%.17g) (%.17g,%.17g,%.17g)\n",
                    i1, P11.X(), P11.Y(), P11.Z(), P21.X(), P21.Y(), P21.Z(), P31.X(), P31.Y(), P31.Z());
      }
      std::printf("LENRAW V longv=%.17g\n", longv / 3.);
    }
    {
      double       longu = 0.0;
      const double du    = 0.05 * (u1 - u0);
      const double vave  = 0.5 * (v1 + v0);
      gp_Pnt       P11 = aD0(u0, v0), P21 = aD0(u0, vave), P31 = aD0(u0, v1);
      std::printf("LENRAW U k=0 (%.17g,%.17g,%.17g) (%.17g,%.17g,%.17g) (%.17g,%.17g,%.17g)\n",
                  P11.X(), P11.Y(), P11.Z(), P21.X(), P21.Y(), P21.Z(), P31.X(), P31.Y(), P31.Z());
      double aU = u0 + du;
      for (int i1 = 1; i1 <= 20; i1++, aU += du)
      {
        gp_Pnt P12 = aD0(aU, v0), P22 = aD0(aU, vave), P32 = aD0(aU, v1);
        longu += (P11.Distance(P12) + P21.Distance(P22) + P31.Distance(P32));
        P11 = P12;
        P21 = P22;
        P31 = P32;
        std::printf("LENRAW U k=%d (%.17g,%.17g,%.17g) (%.17g,%.17g,%.17g) (%.17g,%.17g,%.17g)\n",
                    i1, P11.X(), P11.Y(), P11.Z(), P21.X(), P21.Y(), P21.Z(), P31.X(), P31.Y(), P31.Z());
      }
      std::printf("LENRAW U longu=%.17g\n", longu / 3.);
    }
  }
  return 0;
}

static int runStructDump(const TopoDS_Shape& theShape, const int theFaceIdx,
                         const bool theFactoryAlgo = false, const bool theSrfNodes = false,
                         const double theDefl = 0.0, const double theAngle = 0.0)
{
  Bnd_Box aShapeBox;
  BRepBndLib::Add(theShape, aShapeBox, false);
  double aBox[6] = {0, 0, 0, 0, 0, 0};
  if (!aShapeBox.IsVoid())
  {
    aShapeBox.Get(aBox[0], aBox[1], aBox[2], aBox[3], aBox[4], aBox[5]);
  }
  const double aMaxComp =
    (std::max)((std::max)(aBox[3] - aBox[0], aBox[4] - aBox[1]), aBox[5] - aBox[2]);

  IMeshTools_Parameters aParams;
  // TEMP T-121: the port derives its deflection from the *reader's* own bbox, so
  // the derived value here (`aMaxComp * 0.004`) differs from the port's in the
  // 8th digit, which is enough to change the NURBS/SoR grain filter and therefore
  // the whole surface-node grid. When the caller passes the port's exact numbers,
  // use them so the two DC-pass structures are comparable at all.
  aParams.Deflection    = (theDefl > 0.0) ? theDefl : aMaxComp * 0.001 * 4.0;
  aParams.Angle         = (theAngle > 0.0) ? theAngle : 20.0 * M_PI / 180.0;
  aParams.InParallel    = false;
  aParams.Relative      = false;
  aParams.MinSize       = (std::max)(IMeshTools_Parameters::RelMinSize() * aParams.Deflection,
                                     Precision::Confusion());
  aParams.AdjustMinSize = false;
  // TEMP T-119: `BRepMesh_IncrementalMesh::initParameters` (`hxx:79-107`) fills
  // `DeflectionInterior`/`AngleInterior` when they are left at their sentinel
  // (`-1`). Skipping it feeds the NURBS/SoR splitter a negative interior angle and
  // deflection, which defeats its grain filter and inflates the generated node
  // grid (measured: a 95x95 sphere grid instead of the real one).
  if (aParams.DeflectionInterior < Precision::Confusion())
  {
    aParams.DeflectionInterior = aParams.Deflection;
  }
  if (aParams.MinSize < Precision::Confusion())
  {
    aParams.MinSize =
      (std::max)(IMeshTools_Parameters::RelMinSize()
                   * (std::min)(aParams.Deflection, aParams.DeflectionInterior),
                 Precision::Confusion());
  }
  if (aParams.AngleInterior < Precision::Angular())
  {
    aParams.AngleInterior = 2.0 * aParams.Angle;
  }

  occ::handle<BRepMesh_Context> aCtx = new BRepMesh_Context;
  aCtx->SetShape(theShape);
  aCtx->ChangeParameters()            = aParams;
  aCtx->ChangeParameters().CleanModel = false;
  const bool aBuilt   = aCtx->BuildModel();
  const bool aDiscret = aBuilt && aCtx->DiscretizeEdges();
  const bool aHealed  = aDiscret && aCtx->HealModel();
  const bool aPreproc = aHealed && aCtx->PreProcessModel();
  const occ::handle<IMeshData_Model>& aModel = aCtx->GetModel();
  std::printf("OCCTSTRUCT defl=%.9f angle=%.9f angint=%.9f minsize=%.9f built=%d discret=%d healed=%d "
              "preproc=%d faces=%d\n",
              aParams.Deflection,
              aParams.Angle,
              aParams.AngleInterior,
              aParams.MinSize,
              aBuilt ? 1 : 0,
              aDiscret ? 1 : 0,
              aHealed ? 1 : 0,
              aPreproc ? 1 : 0,
              aModel.IsNull() ? -1 : aModel->FacesNb());
  if (aModel.IsNull() || theFaceIdx < 0 || theFaceIdx >= aModel->FacesNb())
  {
    return 0;
  }

  const IMeshData::IFaceHandle& aDFace = aModel->GetFace(theFaceIdx);
  Bnd_Box                       aFaceBox;
  BRepBndLib::Add(aDFace->GetFace(), aFaceBox, false);
  double v[6] = {0, 0, 0, 0, 0, 0};
  if (!aFaceBox.IsVoid())
  {
    aFaceBox.Get(v[0], v[1], v[2], v[3], v[4], v[5]);
  }
  std::printf("OCCTSTRUCT f=%d type=%d wires=%d bbox=(%.6f,%.6f,%.6f)-(%.6f,%.6f,%.6f)\n",
              theFaceIdx,
              static_cast<int>(aDFace->GetSurface()->GetType()),
              aDFace->WiresNb(),
              v[0],
              v[1],
              v[2],
              v[3],
              v[4],
              v[5]);
  for (int w = 0; w < aDFace->WiresNb(); ++w)
  {
    const IMeshData::IWireHandle& aDWire = aDFace->GetWire(w);
    std::printf("OCCTSTRUCT W %d edges=%d\n", w, aDWire->EdgesNb());
    for (int e = 0; e < aDWire->EdgesNb(); ++e)
    {
      const IMeshData::IEdgeHandle&   aDEdge  = aDWire->GetEdge(e);
      const IMeshData::ListOfInteger& aList   = aDEdge->GetPCurves(aDFace.get());
      for (IMeshData::ListOfInteger::Iterator it(aList); it.More(); it.Next())
      {
        const IMeshData::IPCurveHandle& aPCurve = aDEdge->GetPCurve(it.Value());
        std::printf("OCCTSTRUCT E w=%d e=%d pcId=%d ori=%d fwd=%d n=%d first=(%.9f,%.9f) "
                    "last=(%.9f,%.9f)\n",
                    w,
                    e,
                    it.Value(),
                    static_cast<int>(aPCurve->GetOrientation()),
                    aPCurve->IsForward() ? 1 : 0,
                    aPCurve->ParametersNb(),
                    aPCurve->GetPoint(0).X(),
                    aPCurve->GetPoint(0).Y(),
                    aPCurve->GetPoint(aPCurve->ParametersNb() - 1).X(),
                    aPCurve->GetPoint(aPCurve->ParametersNb() - 1).Y());
        // TEMP EPC: full-precision pcurve points + parameters, so the port's
        // `EPROBE uv=/tt=` list can be compared entry by entry.
        for (int k = 0; k < aPCurve->ParametersNb(); ++k)
        {
          const gp_Pnt2d aP = aPCurve->GetPoint(k);
          std::printf("OCCTEPC w=%d e=%d pcId=%d k=%d u=%.17g v=%.17g t=%.17g\n",
                      w,
                      e,
                      it.Value(),
                      k,
                      aP.X(),
                      aP.Y(),
                      aPCurve->GetParameter(k));
        }
      }
      // TEMP PC2D: the face's own pcurve of this edge
      // (`BRep_Tool::CurveOnSurface`), i.e. what `BRepMesh_ModelBuilder`
      // copies into `IMeshData::IPCurve`, with its dynamic type and a few
      // full-precision evaluations.
      {
        TopLoc_Location aPcLoc;
        (void)aPcLoc;
        double aPcF = 0., aPcL = 0.;
        const occ::handle<Geom2d_Curve>& aPc2d =
          BRep_Tool::CurveOnSurface(aDEdge->GetEdge(), aDFace->GetFace(), aPcF, aPcL);
        if (!aPc2d.IsNull())
        {
          const double aPcMid = 0.5 * (aPcF + aPcL);
          const gp_Pnt2d aPa = aPc2d->Value(aPcF);
          const gp_Pnt2d aPb = aPc2d->Value(aPcMid);
          const gp_Pnt2d aPc = aPc2d->Value(aPcL);
          std::printf("OCCTPC2D w=%d e=%d type=%s f=%.17g l=%.17g pf=(%.17g,%.17g) "
                      "pm=(%.17g,%.17g) pl=(%.17g,%.17g)\n",
                      w,
                      e,
                      aPc2d->DynamicType()->Name(),
                      aPcF,
                      aPcL,
                      aPa.X(),
                      aPa.Y(),
                      aPb.X(),
                      aPb.Y(),
                      aPc.X(),
                      aPc.Y());
        }
        else
        {
          std::printf("OCCTPC2D w=%d e=%d type=NULL\n", w, e);
        }
      }
      // TEMP PC2DORI: the pcurve of THIS edge on the face taken through both
      // edge orientations, i.e. exactly the two objects
      // `BRepMesh_ModelPreProcessor::splitEdge` (`cxx:213-219`) picks between:
      // `aE.Orientation(FORWARD); aPC1 = BRep_Tool::CurveOnSurface(aE, aF, af, al)`
      // and the same with REVERSED. Printing the type, the own range and the
      // value at first/mid/last shows which of the two carries the +-2*pi
      // (seam) branch, i.e. how `splitEdge` can insert a point on the opposite
      // branch of the same physical seam line.
      {
        const TopoDS_Edge aERaw = aDEdge->GetEdge();
        for (int aOriIt = 0; aOriIt < 2; ++aOriIt)
        {
          TopoDS_Edge aEOri = aERaw;
          aEOri.Orientation(aOriIt == 0 ? TopAbs_FORWARD : TopAbs_REVERSED);
          double aOF = 0., aOL = 0.;
          const occ::handle<Geom2d_Curve>& aOC =
            BRep_Tool::CurveOnSurface(aEOri, aDFace->GetFace(), aOF, aOL);
          if (aOC.IsNull())
          {
            std::printf("OCCTPC2DORI w=%d e=%d ori=%d type=NULL\n", w, e, aOriIt);
            continue;
          }
          const double   aOFirst = aOC->FirstParameter();
          const double   aOLast  = aOC->LastParameter();
          const gp_Pnt2d aOA     = aOC->Value(aOFirst);
          const gp_Pnt2d aOB     = aOC->Value(0.5 * (aOFirst + aOLast));
          const gp_Pnt2d aOCc    = aOC->Value(aOLast);
          std::printf("OCCTPC2DORI w=%d e=%d ori=%d type=%s range=%.17g %.17g evf=%.17g evl=%.17g "
                      "pf=(%.17g,%.17g) pm=(%.17g,%.17g) pl=(%.17g,%.17g)\n",
                      w,
                      e,
                      aOriIt,
                      aOC->DynamicType()->Name(),
                      aOFirst,
                      aOLast,
                      aOF,
                      aOL,
                      aOA.X(),
                      aOA.Y(),
                      aOB.X(),
                      aOB.Y(),
                      aOCc.X(),
                      aOCc.Y());
        }
        std::printf("OCCTPC2DFLAGS w=%d e=%d sameparam=%d samerange=%d pcnb=%d\n",
                    w,
                    e,
                    BRep_Tool::SameParameter(aERaw) ? 1 : 0,
                    BRep_Tool::SameRange(aERaw) ? 1 : 0,
                    aDEdge->PCurvesNb());
      }
    }
  }

  occ::handle<IMeshTools_MeshAlgo> aAlgo;
  if (theSrfNodes)
  {
    // TEMP T-120: factory-faithful splitter per `BRepMesh_MeshAlgoFactory::GetAlgo`
    // (`cxx:60-118`), node-insertion algo (the classifier filter under test is the
    // same in `DelaunayNodeInsertionMeshAlgo` and `DelaunayDeflectionControlMeshAlgo`).
    switch (aDFace->GetSurface()->GetType())
    {
      case GeomAbs_Plane:
        aAlgo = new ZZSrfNodeCountAlgo<BRepMesh_DefaultRangeSplitter>();
        break;
      case GeomAbs_Sphere:
        aAlgo = new ZZSrfNodeCountAlgo<BRepMesh_SphereRangeSplitter>();
        break;
      case GeomAbs_Cylinder:
        aAlgo = new ZZSrfNodeCountAlgo<BRepMesh_CylinderRangeSplitter>();
        break;
      case GeomAbs_Cone: aAlgo = new ZZSrfNodeCountAlgo<BRepMesh_ConeRangeSplitter>(); break;
      case GeomAbs_Torus: aAlgo = new ZZSrfNodeCountAlgo<BRepMesh_TorusRangeSplitter>(); break;
      case GeomAbs_SurfaceOfRevolution:
        aAlgo = new ZZSrfNodeCountAlgo<BRepMesh_BoundaryParamsRangeSplitter>();
        break;
      case GeomAbs_SurfaceOfExtrusion:
        aAlgo = new ZZSrfNodeCountAlgo<BRepMesh_ExtrusionRangeSplitter>();
        break;
      case GeomAbs_BezierSurface:
      case GeomAbs_BSplineSurface:
        aAlgo = new ZZSrfNodeCountAlgo<BRepMesh_NURBSRangeSplitter>();
        break;
      default: aAlgo = new ZZSrfNodeCountAlgo<BRepMesh_UndefinedRangeSplitter>(); break;
    }
    aAlgo->Perform(aDFace, aParams, Message_ProgressRange());
    return 0;
  }
  if (theFactoryAlgo)
  {
    // TEMP T-119: mirror `BRepMesh_MeshAlgoFactory::GetAlgo` (`cxx:60-118`) with
    // the probe's parameters (`EnableControlSurfaceDeflectionAllSurfaces` false,
    // `InternalVerticesMode` true). Only the DC-capable pairs are listed; the
    // `BaseMeshAlgo::Type` fallback for plane/cylinder is not installed here
    // because this probe always wants a dumpable structure.
    switch (aDFace->GetSurface()->GetType())
    {
      case GeomAbs_Sphere:
        aAlgo = new ZZStructDumpDCAlgo<BRepMesh_SphereRangeSplitter>();
        break;
      case GeomAbs_Cone: aAlgo = new ZZStructDumpDCAlgo<BRepMesh_ConeRangeSplitter>(); break;
      case GeomAbs_Torus: aAlgo = new ZZStructDumpDCAlgo<BRepMesh_TorusRangeSplitter>(); break;
      case GeomAbs_SurfaceOfRevolution:
        aAlgo = new ZZStructDumpDCAlgo<BRepMesh_BoundaryParamsRangeSplitter>();
        break;
      case GeomAbs_SurfaceOfExtrusion:
        aAlgo = new ZZStructDumpDCAlgo<BRepMesh_ExtrusionRangeSplitter>();
        break;
      case GeomAbs_BezierSurface:
      case GeomAbs_BSplineSurface:
        aAlgo = new ZZStructDumpDCAlgo<BRepMesh_NURBSRangeSplitter>();
        break;
      case GeomAbs_Cylinder:
        aAlgo = new ZZStructDumpDCAlgo<BRepMesh_CylinderRangeSplitter>();
        break;
      case GeomAbs_Plane:
        aAlgo = new ZZStructDumpDCAlgo<BRepMesh_DefaultRangeSplitter>();
        break;
      default: aAlgo = new ZZStructDumpDCAlgo<BRepMesh_UndefinedRangeSplitter>(); break;
    }
    aAlgo->Perform(aDFace, aParams, Message_ProgressRange());
    return 0;
  }

  switch (aDFace->GetSurface()->GetType())
  {
    case GeomAbs_Plane: aAlgo = new ZZStructDumpAlgo<BRepMesh_DefaultRangeSplitter>(); break;
    case GeomAbs_Sphere: aAlgo = new ZZStructDumpAlgo<BRepMesh_SphereRangeSplitter>(); break;
    case GeomAbs_Cylinder:
      aAlgo = new ZZStructDumpAlgo<BRepMesh_CylinderRangeSplitter>();
      break;
    case GeomAbs_Cone: aAlgo = new ZZStructDumpAlgo<BRepMesh_ConeRangeSplitter>(); break;
    case GeomAbs_Torus: aAlgo = new ZZStructDumpAlgo<BRepMesh_TorusRangeSplitter>(); break;
    default: aAlgo = new ZZStructDumpAlgo<BRepMesh_DefaultRangeSplitter>(); break;
  }
  aAlgo->Perform(aDFace, aParams, Message_ProgressRange());
  return 0;
}

// --- TEMP T-101 (`--facestats <defl> <angle>`) -------------------------------------
// Per-model-face `nodes/triangles/bbox` with the *same* mesher parameters the port
// uses, in the model's own face order (`IMeshData_Model::GetFace`), so a face can be
// paired with the port's `zz_probe_a3n00 --fstats` by the six bbox coordinates.
// `--uvsum` cannot be used for this: that branch meshes with a hard-coded 0.5 rad
// angle (its per-face `triangles=` sums to 8054 on a3n00 while the 20-degree mesh has
// 12324), so its counts are not comparable with the port's.
// TEMP T-132 (`--srparams <faceIdx> <defl> <angle>`): prints the raw U/V parameter
// sets a factory-selected UV range splitter hands to its grain filter, together with
// the grid those parameters produce (`BRepMesh_NURBSRangeSplitter::
// GenerateSurfaceNodes`, BRepMesh_NURBSRangeSplitter.cxx:314-400), so the port's
// `NURBSGRID` dump can be compared entry by entry instead of by node totals.
static void dumpRawParams(const char* theLabel, const IMeshData::IMapOfReal& theParams)
{
  std::vector<double> aVals;
  aVals.reserve(theParams.Extent());
  for (int i = 1; i <= theParams.Extent(); ++i)
  {
    aVals.push_back(theParams(i));
  }
  std::sort(aVals.begin(), aVals.end());
  std::printf("OCCTSRP param %s n=%d vals=", theLabel, static_cast<int>(aVals.size()));
  for (size_t i = 0; i < aVals.size(); ++i)
  {
    std::printf("%s%.17g", (i == 0) ? "" : ",", aVals[i]);
  }
  std::printf("\n");
}

// TEMP T-132: same splitter/algo pairing as `ZZSrfNodeCountAlgo`, but the dump is the
// splitter's own parameter state (`GetParametersU/V`, `BRepMesh_UVParamRangeSplitter.
// hxx:46-58`) plus the generated grid. The parameters are filled during
// `initDataStructure` -> `collectWirePoints` -> `AddPoint`
// (`BRepMesh_NodeInsertionMeshAlgo.hxx:143-183`, `BRepMesh_BoundaryParamsRangeSplitter.
// hxx:33-38`) and never change afterwards, so reading them in `postProcessMesh` is the
// same state `insertNodes` filters.
template <class Splitter>
class ZZSrfParamsAlgo
    : public BRepMesh_DelaunayNodeInsertionMeshAlgo<Splitter, BRepMesh_DelaunayBaseMeshAlgo>
{
  typedef BRepMesh_DelaunayNodeInsertionMeshAlgo<Splitter, BRepMesh_DelaunayBaseMeshAlgo> Base;

protected:
  void postProcessMesh(BRepMesh_Delaun& theMesher, const Message_ProgressRange& theRange) override
  {
    const Splitter&                 aSplitter = this->getRangeSplitter();
    const IMeshData::IFaceHandle&   aDFace    = this->getDFace();
    const std::pair<double, double>& aRangeU   = aSplitter.GetRangeU();
    const std::pair<double, double>& aRangeV   = aSplitter.GetRangeV();
    const std::pair<double, double>& aDelta    = aSplitter.GetDelta();
    const std::pair<double, double>& aTol      = aSplitter.GetToleranceUV();
    const IMeshTools_Parameters&    aParams   = this->getParameters();
    const double                    aDefl     = aDFace->GetDeflection();
    std::printf("OCCTSRP face type=%d valid=%d ru=(%.17g,%.17g) rv=(%.17g,%.17g) delta=(%.17g,%.17g) "
                "tol=(%.17g,%.17g) ures=%.17g vres=%.17g facedefl=%.17g minsize=%.17g "
                "angint=%.17g angle=%.17g defl=%.17g\n",
                static_cast<int>(aDFace->GetSurface()->GetType()),
                const_cast<Splitter&>(aSplitter).IsValid() ? 1 : 0,
                aRangeU.first,
                aRangeU.second,
                aRangeV.first,
                aRangeV.second,
                aDelta.first,
                aDelta.second,
                aTol.first,
                aTol.second,
                aSplitter.GetSurface()->UResolution(aDefl),
                aSplitter.GetSurface()->VResolution(aDefl),
                aDefl,
                aParams.MinSize,
                aParams.AngleInterior,
                aParams.Angle,
                aParams.Deflection);
    dumpRawParams("u", aSplitter.GetParametersU());
    dumpRawParams("v", aSplitter.GetParametersV());

    // `GenerateSurfaceNodes` fills the parameter maps on its first call for the
    // splitters that derive them lazily (`BRepMesh_NURBSRangeSplitter::
    // initParameters`, cxx:452-488, writes through the const_casts at :469/:476),
    // so the grid is built first and both maps are dumped a second time after it.
    const Handle(IMeshData::ListOfPnt2d) aNodes = aSplitter.GenerateSurfaceNodes(aParams);
    dumpRawParams("u_after", aSplitter.GetParametersU());
    dumpRawParams("v_after", aSplitter.GetParametersV());
    std::printf("OCCTSRP nodes=%d\n", aNodes.IsNull() ? 0 : static_cast<int>(aNodes->Size()));
    if (!aNodes.IsNull())
    {
      for (IMeshData::ListOfPnt2d::Iterator anIt(*aNodes); anIt.More(); anIt.Next())
      {
        std::printf("OCCTSRPN u=%.17g v=%.17g\n", anIt.Value().X(), anIt.Value().Y());
      }
    }
    Base::postProcessMesh(theMesher, theRange);
  }
};

static int runSrfParams(const TopoDS_Shape& theShape, const int theFaceIdx,
                        const double theDeflection, const double theAngle)
{
  IMeshTools_Parameters aParams;
  aParams.Deflection    = theDeflection;
  aParams.Angle         = theAngle;
  aParams.InParallel    = false;
  aParams.Relative      = false;
  aParams.MinSize       = 0.0;
  aParams.AdjustMinSize = false;
  // `BRepMesh_IncrementalMesh::initParameters` (BRepMesh_IncrementalMesh.hxx:79-108),
  // i.e. the parameter state the real mesher runs with.
  if (aParams.DeflectionInterior < Precision::Confusion())
  {
    aParams.DeflectionInterior = aParams.Deflection;
  }
  if (aParams.MinSize < Precision::Confusion())
  {
    aParams.MinSize =
      (std::max)(IMeshTools_Parameters::RelMinSize()
                   * (std::min)(aParams.Deflection, aParams.DeflectionInterior),
                 Precision::Confusion());
  }
  if (aParams.AngleInterior < Precision::Angular())
  {
    aParams.AngleInterior = 2.0 * aParams.Angle;
  }
  occ::handle<BRepMesh_Context> aCtx = new BRepMesh_Context;
  aCtx->SetShape(theShape);
  aCtx->ChangeParameters()            = aParams;
  aCtx->ChangeParameters().CleanModel = false;
  const bool aBuilt   = aCtx->BuildModel();
  const bool aDiscret = aBuilt && aCtx->DiscretizeEdges();
  const bool aHealed  = aDiscret && aCtx->HealModel();
  const bool aPreproc = aHealed && aCtx->PreProcessModel();
  const occ::handle<IMeshData_Model>& aModel = aCtx->GetModel();
  if (aModel.IsNull() || theFaceIdx < 0 || theFaceIdx >= aModel->FacesNb())
  {
    std::cout << "SRPARAM model-null-or-range\n";
    return 0;
  }
  const IMeshData::IFaceHandle& aDFace = aModel->GetFace(theFaceIdx);
  std::printf("OCCTSRP shape defl=%.17g angle=%.17g angint=%.17g minsize=%.17g built=%d discret=%d "
              "healed=%d preproc=%d faces=%d f=%d\n",
              aParams.Deflection,
              aParams.Angle,
              aParams.AngleInterior,
              aParams.MinSize,
              aBuilt ? 1 : 0,
              aDiscret ? 1 : 0,
              aHealed ? 1 : 0,
              aPreproc ? 1 : 0,
              aModel->FacesNb(),
              theFaceIdx);
  occ::handle<IMeshTools_MeshAlgo> aAlgo;
  switch (aDFace->GetSurface()->GetType())
  {
    case GeomAbs_SurfaceOfRevolution:
      aAlgo = new ZZSrfParamsAlgo<BRepMesh_BoundaryParamsRangeSplitter>();
      break;
    case GeomAbs_SurfaceOfExtrusion:
      aAlgo = new ZZSrfParamsAlgo<BRepMesh_ExtrusionRangeSplitter>();
      break;
    case GeomAbs_BezierSurface:
    case GeomAbs_BSplineSurface:
      aAlgo = new ZZSrfParamsAlgo<BRepMesh_NURBSRangeSplitter>();
      break;
    default:
      std::printf("OCCTSRP unsupported stype=%d\n",
                  static_cast<int>(aDFace->GetSurface()->GetType()));
      return 0;
  }
  aAlgo->Perform(aDFace, aParams, Message_ProgressRange());
  return 0;
}

static int runFaceStats(const TopoDS_Shape& theShape, const double theDeflection,
                        const double theAngle)
{
  IMeshTools_Parameters aParams;
  aParams.Deflection    = theDeflection;
  aParams.Angle         = theAngle;
  aParams.InParallel    = false;
  aParams.Relative      = false;
  aParams.MinSize       = Precision::Confusion();
  aParams.AdjustMinSize = false;

  occ::handle<BRepMesh_ModelBuilder> aBuilder = new BRepMesh_ModelBuilder;
  const occ::handle<IMeshData_Model> aModel   = aBuilder->Perform(theShape, aParams);
  if (aModel.IsNull())
  {
    std::cout << "FSTAT_OCC model-null\n";
    return 0;
  }
  {
    BRepMesh_IncrementalMesh aMesher(theShape, theDeflection, false, theAngle);
  }
  std::cout << "FSTAT_OCC defl=" << theDeflection << " angle=" << theAngle
            << " faces=" << aModel->FacesNb() << std::endl;
  for (int f = 0; f < aModel->FacesNb(); ++f)
  {
    const IMeshData::IFaceHandle& aDFace = aModel->GetFace(f);
    const TopoDS_Face&            aFace  = aDFace->GetFace();
    Bnd_Box                       aBox;
    BRepBndLib::Add(aFace, aBox, false);
    double v[6] = {0, 0, 0, 0, 0, 0};
    if (!aBox.IsVoid())
    {
      aBox.Get(v[0], v[1], v[2], v[3], v[4], v[5]);
    }
    TopLoc_Location                     aLoc;
    const occ::handle<Poly_Triangulation>& aTri = BRep_Tool::Triangulation(aFace, aLoc);
    const int aNodes = aTri.IsNull() ? 0 : aTri->NbNodes();
    const int aTris  = aTri.IsNull() ? 0 : aTri->NbTriangles();
    std::cout << std::fixed << std::setprecision(6) << "FSTAT_OCC face=" << f
              << " nodes=" << aNodes << " triangles=" << aTris << " bbox=(" << v[0] << "," << v[1]
              << "," << v[2] << ")-(" << v[3] << "," << v[4] << "," << v[5] << ")" << std::endl;
  }
  return 0;
}

// --- TEMP T-121 (`--edgedu <faceIdx> <defl> <angle>`) -------------------------------
// TEMP T-101b (`--facedump <faceIdx> <defl> <angle>`): the node and triangle arrays of
// one *model* face (`IMeshData_Model::GetFace`, same indexing as `--facestats`), so a
// per-face count divergence can be read as geometry instead of a total.
static int runFaceDump(const TopoDS_Shape& theShape, const int theFaceIdx,
                       const double theDeflection, const double theAngle)
{
  {
    BRepMesh_IncrementalMesh aMesher(theShape, theDeflection, false, theAngle);
  }
  IMeshTools_Parameters aParams;
  aParams.Deflection    = theDeflection;
  aParams.Angle         = theAngle;
  aParams.InParallel    = false;
  aParams.Relative      = false;
  aParams.MinSize       = Precision::Confusion();
  aParams.AdjustMinSize = false;
  occ::handle<BRepMesh_ModelBuilder> aBuilder = new BRepMesh_ModelBuilder;
  const occ::handle<IMeshData_Model> aModel   = aBuilder->Perform(theShape, aParams);
  if (aModel.IsNull() || theFaceIdx < 0 || theFaceIdx >= aModel->FacesNb())
  {
    std::cout << "FDUMP model-null-or-range\n";
    return 0;
  }
  const TopoDS_Face& aFace = aModel->GetFace(theFaceIdx)->GetFace();
  TopLoc_Location    aLoc;
  const occ::handle<Poly_Triangulation>& aTri = BRep_Tool::Triangulation(aFace, aLoc);
  std::cout << std::fixed << std::setprecision(9) << "FDUMP face=" << theFaceIdx
            << " nodes=" << (aTri.IsNull() ? 0 : aTri->NbNodes())
            << " tris=" << (aTri.IsNull() ? 0 : aTri->NbTriangles()) << std::endl;
  if (aTri.IsNull())
    return 0;
  for (int i = 1; i <= aTri->NbNodes(); ++i)
  {
    const gp_Pnt aP = aTri->Node(i).Transformed(aLoc.Transformation());
    std::cout << "FDUMP n i=" << i;
    if (aTri->HasUVNodes())
    {
      const gp_Pnt2d aUV = aTri->UVNode(i);
      std::cout << " uv=(" << aUV.X() << "," << aUV.Y() << ")";
    }
    std::cout << " p=(" << aP.X() << "," << aP.Y() << "," << aP.Z() << ")" << std::endl;
  }
  for (int i = 1; i <= aTri->NbTriangles(); ++i)
  {
    int a1, a2, a3;
    aTri->Triangle(i).Get(a1, a2, a3);
    std::cout << "FDUMP t i=" << i << " (" << a1 << "," << a2 << "," << a3 << ")" << std::endl;
  }
  return 0;
}

// --- TEMP T-121 (`--edgedu <faceIdx> <defl> <angle>`) -------------------------------// Prints, for every edge occurrence of one face, the `BRepAdaptor_Curve` type and the
// exact inputs of `GCPnts_TangentialDeflection::EvaluateDu` (`cxx:302-322`) at the
// edge's first parameter: `|T|`, `|N|`, `|N x T|` and `Ln = |N x T| / |T|`. The
// resulting `du_init = sqrt(8 * max(curvatureDeflection, minLen) / Ln)` is printed for
// both the face deflection and half of it, so a divergence in the first tessellation
// step can be attributed to the curve (`Ln`) or to the deflection (`D`).
static int runEdgeDu(const TopoDS_Shape& theShape, const int theFaceIdx,
                     const double theDeflection, const double theAngle)
{
  {
    BRepMesh_IncrementalMesh aMesher(theShape, theDeflection, false, theAngle);
  }
  int aFaceIdx = 0;
  for (TopExp_Explorer aFaceEx(theShape, TopAbs_FACE); aFaceEx.More(); aFaceEx.Next(), ++aFaceIdx)
  {
    if (aFaceIdx != theFaceIdx)
      continue;
    const TopoDS_Face aFace = TopoDS::Face(aFaceEx.Current());
    int               anEdgeIdx = 0;
    for (TopExp_Explorer anEdgeEx(aFace, TopAbs_EDGE); anEdgeEx.More();
         anEdgeEx.Next(), ++anEdgeIdx)
    {
      const TopoDS_Edge anEdge = TopoDS::Edge(anEdgeEx.Current());
      std::printf("EDGEDU face=%d i=%d degenerated=%d\n",
                  aFaceIdx,
                  anEdgeIdx,
                  BRep_Tool::Degenerated(anEdge) ? 1 : 0);
      try
      {
        BRepAdaptor_Curve  ac(anEdge);
        const GeomAbs_CurveType aType = ac.GetType();
        const double       aF = ac.FirstParameter();
        const double       aL = ac.LastParameter();
        gp_Pnt             aP;
        gp_Vec             aT, aN;
        ac.D2(aF, aP, aT, aN);
        const double aLt    = aT.Magnitude();
        const double aNMag  = aN.Magnitude();
        const double aCross = aN.CrossMagnitude(aT);
        const double aLn    = (aLt > Precision::Confusion()) ? aCross / aLt : 0.0;
        std::printf("EDGEDU face=%d i=%d type=%d f=%.17g l=%.17g |T|=%.17g |N|=%.17g "
                    "|NxT|=%.17g Ln=%.17g p0=(%.17g,%.17g,%.17g) du(D)=%.17g du(D/2)=%.17g\n",
                    aFaceIdx,
                    anEdgeIdx,
                    static_cast<int>(aType),
                    aF,
                    aL,
                    aLt,
                    aNMag,
                    aCross,
                    aLn,
                    aP.X(),
                    aP.Y(),
                    aP.Z(),
                    (aLn > Precision::Confusion()) ? std::sqrt(8.0 * theDeflection / aLn) : 0.0,
                    (aLn > Precision::Confusion()) ? std::sqrt(8.0 * (theDeflection / 2.0) / aLn)
                                                   : 0.0);
        if (aType == GeomAbs_Circle)
        {
          std::printf("EDGEDU face=%d i=%d circle R=%.17g\n",
                      aFaceIdx,
                      anEdgeIdx,
                      ac.Circle().Radius());
        }
      }
      catch (const Standard_Failure& anErr)
      {
        std::printf("EDGEDU face=%d i=%d thrown=%s\n",
                    aFaceIdx,
                    anEdgeIdx,
                    anErr.GetMessageString() ? anErr.GetMessageString() : "?");
      }
    }
  }
  return 0;
}

// --- TEMP T-103 (`--torusdbg <defl> <angle>`) --------------------------------------
// Replays OCCT's own `BRepMesh_TorusRangeSplitter` on every torus face of the
// discrete model built with the given parameters, printing the exact interior
// node count it produces. `_dbg/ZZ_TorusRangeSplitter.cxx` supplies both the
// byte-identical splitter copy and `ZZTorDbgFace`.
void ZZTorDbgFace(const IMeshData::IFaceHandle& theDFace,
                  const IMeshTools_Parameters&  theParameters,
                  const int                     theFaceIdx);

static int runTorusDbg(const TopoDS_Shape& theShape, const double theDeflection,
                       const double theAngle)
{
  IMeshTools_Parameters aParams;
  aParams.Deflection    = theDeflection;
  aParams.Angle         = theAngle;
  aParams.InParallel    = false;
  aParams.Relative      = false;
  // Mirror `BRepMesh_IncrementalMesh::initParameters` (`hxx:91-96`): the real
  // mesher runs with `MinSize = max(RelMinSize * min(defl, deflInterior), Confusion)`.
  // The earlier Confusion setting made `BRepMesh_EdgeDiscret` emit finer pcurves
  // than the real run, so the splitter inputs were not comparable.
  aParams.MinSize       = (std::max)(IMeshTools_Parameters::RelMinSize() * theDeflection,
                                     Precision::Confusion());
  aParams.AdjustMinSize = false;

  occ::handle<BRepMesh_ModelBuilder> aBuilder = new BRepMesh_ModelBuilder;
  const occ::handle<IMeshData_Model> aModel   = aBuilder->Perform(theShape, aParams);
  if (aModel.IsNull())
  {
    std::cout << "ZZTORDBG model-null\n";
    return 0;
  }
  // The builder only enumerates the faces; the pcurves that
  // `collectWirePoints` walks are filled by `BRepMesh_EdgeDiscret`
  // (`BRepMesh_ModelBuilder` + `BRepMesh_ModelHealer` + `BRepMesh_EdgeDiscret`
  // is the `IMeshTools_Context` chain). Without this every replay reports
  // `ParametersNb() == 0` and `AdjustRange` leaves the splitter invalid.
  occ::handle<BRepMesh_EdgeDiscret> anEdgeDiscret = new BRepMesh_EdgeDiscret;
  anEdgeDiscret->Perform(aModel, aParams, Message_ProgressRange());
  for (int f = 0; f < aModel->FacesNb(); ++f)
  {
    const IMeshData::IFaceHandle& aDFace = aModel->GetFace(f);
    if (aDFace->GetSurface()->GetType() != GeomAbs_Torus)
    {
      continue;
    }
    ZZTorDbgFace(aDFace, aParams, f);
  }
  return 0;
}

// --- TEMP T-103 (`--healerdbg <defl> <angle> <x1> <y1> <z1> <x2> <y2> <z2>`) ---------
// Dumps the per-wire-edge pcurve state of every discrete face whose bbox matches
// the given target, BEFORE and AFTER `BRepMesh_ModelHealer::Perform`. This is the
// oracle for the port's `ModelHealer::fix_face_boundaries` snap: the pre-heal
// columns must be bit-equal to the port's EdgeDiscret output, and the post-heal
// columns show what OCCT's own `connectClosestPoints`/`adjustSamePoints` write.
static void dumpFacePcurves(const IMeshData::IFaceHandle& theDFace,
                            const int                     theFaceIdx,
                            const char*                   theTag)
{
  for (int aWireIt = 0; aWireIt < theDFace->WiresNb(); ++aWireIt)
  {
    const IMeshData::IWireHandle& aDWire = theDFace->GetWire(aWireIt);
    for (int aEdgeIt = 0; aEdgeIt < aDWire->EdgesNb(); ++aEdgeIt)
    {
      const IMeshData::IEdgeHandle&   aDEdge  = aDWire->GetEdge(aEdgeIt);
      const IMeshData::IPCurveHandle& aPCurve = aDEdge->GetPCurve(
        theDFace.get(), aDWire->GetEdgeOrientation(aEdgeIt));
      const int aNb = aPCurve->ParametersNb();
      std::printf("HDBG %s face=%d w=%d e=%d nb=%d fwd=%d ori=%d same=%d free=%d int=%d "
                  "first=(%.9f,%.9f) last=(%.9f,%.9f)\n",
                  theTag,
                  theFaceIdx,
                  aWireIt,
                  aEdgeIt,
                  aNb,
                  aPCurve->IsForward() ? 1 : 0,
                  (int)aPCurve->GetOrientation(),
                  aDEdge->GetSameParam() ? 1 : 0,
                  aDEdge->IsFree() ? 1 : 0,
                  aPCurve->IsInternal() ? 1 : 0,
                  aNb > 0 ? aPCurve->GetPoint(0).X() : 0.0,
                  aNb > 0 ? aPCurve->GetPoint(0).Y() : 0.0,
                  aNb > 0 ? aPCurve->GetPoint(aNb - 1).X() : 0.0,
                  aNb > 0 ? aPCurve->GetPoint(aNb - 1).Y() : 0.0);
      if (aNb > 0 && aNb <= 64)
      {
        std::printf("HDBG %s face=%d e=%d pts=[", theTag, theFaceIdx, aEdgeIt);
        for (int k = 0; k < aNb; ++k)
        {
          std::printf("%s(%.6f,%.6f)", k ? " " : "", aPCurve->GetPoint(k).X(), aPCurve->GetPoint(k).Y());
        }
        std::printf("]\n");
      }
    }
  }
}

static int runHealerDbg(const TopoDS_Shape& theShape,
                        const double        theDeflection,
                        const double        theAngle,
                        const double*       theBBox)
{
  IMeshTools_Parameters aParams;
  aParams.Deflection    = theDeflection;
  aParams.Angle         = theAngle;
  aParams.InParallel    = false;
  aParams.Relative      = false;
  aParams.MinSize       = (std::max)(IMeshTools_Parameters::RelMinSize() * theDeflection,
                                     Precision::Confusion());
  aParams.AdjustMinSize = false;

  occ::handle<BRepMesh_ModelBuilder> aBuilder = new BRepMesh_ModelBuilder;
  const occ::handle<IMeshData_Model> aModel   = aBuilder->Perform(theShape, aParams);
  if (aModel.IsNull())
  {
    std::cout << "HDBG model-null" << std::endl;
    return 0;
  }
  occ::handle<BRepMesh_EdgeDiscret> anEdgeDiscret = new BRepMesh_EdgeDiscret;
  anEdgeDiscret->Perform(aModel, aParams, Message_ProgressRange());

  std::vector<int> aHits;
  for (int f = 0; f < aModel->FacesNb(); ++f)
  {
    const TopoDS_Face& aFace = aModel->GetFace(f)->GetFace();
    Bnd_Box            aBox;
    BRepBndLib::Add(aFace, aBox, false);
    if (aBox.IsVoid())
    {
      continue;
    }
    double v[6] = {0, 0, 0, 0, 0, 0};
    aBox.Get(v[0], v[1], v[2], v[3], v[4], v[5]);
    if (std::abs(v[0] - theBBox[0]) < 1e-3 && std::abs(v[1] - theBBox[1]) < 1e-3
        && std::abs(v[2] - theBBox[2]) < 1e-3 && std::abs(v[3] - theBBox[3]) < 1e-3
        && std::abs(v[4] - theBBox[4]) < 1e-3 && std::abs(v[5] - theBBox[5]) < 1e-3)
    {
      aHits.push_back(f);
    }
  }
  std::cout << "HDBG hits=" << aHits.size() << std::endl;
  for (size_t i = 0; i < aHits.size(); ++i)
  {
    dumpFacePcurves(aModel->GetFace(aHits[i]), aHits[i], "pre");
  }

  occ::handle<BRepMesh_ModelHealer> aHealer = new BRepMesh_ModelHealer;
  aHealer->Perform(aModel, aParams, Message_ProgressRange());

  for (size_t i = 0; i < aHits.size(); ++i)
  {
    dumpFacePcurves(aModel->GetFace(aHits[i]), aHits[i], "post");
  }

  // TEMP T-122: the face flags the healer left behind (`BRepMesh_ModelHealer::
  // Perform`, `cxx:156-166`: a FaceChecker failure marks SelfIntersectingWire +
  // Failure, and `BRepMesh_FaceDiscret::process` (`cxx:102`) then skips the face
  // entirely, so it ends up with no triangulation at all).
  for (size_t i = 0; i < aHits.size(); ++i)
  {
    const IMeshData::IFaceHandle& aDFace = aModel->GetFace(aHits[i]);
    int aOpenWire = 0;
    for (int w = 0; w < aDFace->WiresNb(); ++w)
    {
      if (aDFace->GetWire(w)->IsSet(IMeshData_OpenWire))
      {
        aOpenWire = 1;
      }
    }
    std::cout << "HDBGFLAG f=" << aHits[i] << " wires=" << aDFace->WiresNb()
              << " fail=" << (aDFace->IsSet(IMeshData_Failure) ? 1 : 0)
              << " selfint=" << (aDFace->IsSet(IMeshData_SelfIntersectingWire) ? 1 : 0)
              << " open=" << aOpenWire
              << " outdated=" << (aDFace->IsSet(IMeshData_Outdated) ? 1 : 0)
              << " unoriented=" << (aDFace->IsSet(IMeshData_UnorientedWire) ? 1 : 0)
              << " defl=" << std::setprecision(17) << aDFace->GetDeflection() << std::endl;
  }
  return 0;
}

// --- TEMP T-122 (`--ctxdbg <defl> <angle> <x1> <y1> <z1> <x2> <y2> <z2>`) -------------
// Runs the real `IMeshTools_MeshBuilder` stage sequence (`IMeshTools_MeshBuilder.cxx:
// 52-62`: BuildModel -> DiscretizeEdges -> HealModel -> PreProcessModel ->
// DiscretizeFaces -> PostProcessModel) on the bbox-selected faces and dumps their
// `IMeshData_*` status flags plus the resulting `Poly_Triangulation` size after EVERY
// stage. The stage whose dump flips a flag - or whose triangulation is still missing -
// is the stage that leaves the face unmeshed, so the behaviour gets ported branch by
// branch instead of guessed from the final state. The parameter block mirrors
// `BRepMesh_IncrementalMesh::initParameters` (`BRepMesh_IncrementalMesh.hxx:79-108`);
// the per-stage totals as well as the final `nodes=`/`tris=` must match
// `--mesh <defl> <angle>` on the same shape, which is this replay's self-check.
namespace
{
class ZZStageCtx : public BRepMesh_Context
{
public:
  ZZStageCtx(const double* theBBox)
  {
    for (int i = 0; i < 6; ++i)
    {
      myBBox[i] = theBBox[i];
    }
  }

  bool BuildModel() override
  {
    const bool aRes = BRepMesh_Context::BuildModel();
    findFaces();
    dumpStage("build", aRes);
    return aRes;
  }

  bool DiscretizeEdges() override
  {
    const bool aRes = BRepMesh_Context::DiscretizeEdges();
    dumpStage("disc", aRes);
    return aRes;
  }

  bool HealModel() override
  {
    const bool aRes = BRepMesh_Context::HealModel();
    dumpStage("heal", aRes);
    return aRes;
  }

  bool PreProcessModel() override
  {
    const bool aRes = BRepMesh_Context::PreProcessModel();
    dumpStage("preproc", aRes);
    return aRes;
  }

  bool DiscretizeFaces(const Message_ProgressRange& theRange) override
  {
    const bool aRes = BRepMesh_Context::DiscretizeFaces(theRange);
    dumpStage("faces", aRes);
    return aRes;
  }

  bool PostProcessModel() override
  {
    const bool aRes = BRepMesh_Context::PostProcessModel();
    dumpStage("postproc", aRes);
    return aRes;
  }

private:
  void findFaces()
  {
    const occ::handle<IMeshData_Model>& aModel = GetModel();
    if (aModel.IsNull())
    {
      return;
    }
    for (int f = 0; f < aModel->FacesNb(); ++f)
    {
      Bnd_Box aBox;
      BRepBndLib::Add(aModel->GetFace(f)->GetFace(), aBox, false);
      if (aBox.IsVoid())
      {
        continue;
      }
      double v[6] = {0, 0, 0, 0, 0, 0};
      aBox.Get(v[0], v[1], v[2], v[3], v[4], v[5]);
      bool aHit = true;
      for (int k = 0; k < 6; ++k)
      {
        if (std::abs(v[k] - myBBox[k]) >= 1e-3)
        {
          aHit = false;
        }
      }
      if (aHit)
      {
        myHits.push_back(f);
      }
    }
  }

  void dumpStage(const char* theTag, const bool theRes) const
  {
    const occ::handle<IMeshData_Model>& aModel = GetModel();
    std::cout << "CTXDBG " << theTag << " res=" << (theRes ? 1 : 0)
              << " model_null=" << (aModel.IsNull() ? 1 : 0);
    if (aModel.IsNull())
    {
      std::cout << std::endl;
      return;
    }
    int aNodes = 0, aTris = 0;
    for (int f = 0; f < aModel->FacesNb(); ++f)
    {
      TopLoc_Location                        aLoc;
      const occ::handle<Poly_Triangulation>& aTri =
        BRep_Tool::Triangulation(aModel->GetFace(f)->GetFace(), aLoc);
      if (!aTri.IsNull())
      {
        aNodes += aTri->NbNodes();
        aTris += aTri->NbTriangles();
      }
    }
    std::cout << " faces=" << aModel->FacesNb() << " nodes=" << aNodes << " tris=" << aTris
              << std::endl;
    for (size_t i = 0; i < myHits.size(); ++i)
    {
      const IMeshData::IFaceHandle& aDFace = aModel->GetFace(myHits[i]);
      int                           aOpenWire = 0;
      for (int w = 0; w < aDFace->WiresNb(); ++w)
      {
        if (aDFace->GetWire(w)->IsSet(IMeshData_OpenWire))
        {
          aOpenWire = 1;
        }
      }
      TopLoc_Location                        aLoc;
      const occ::handle<Poly_Triangulation>& aTri =
        BRep_Tool::Triangulation(aDFace->GetFace(), aLoc);
      std::cout << "CTXFLAG " << theTag << " f=" << myHits[i]
                << " stype=" << static_cast<int>(aDFace->GetSurface()->GetType())
                << " fail=" << (aDFace->IsSet(IMeshData_Failure) ? 1 : 0)
                << " reused=" << (aDFace->IsSet(IMeshData_Reused) ? 1 : 0)
                << " remesh=" << (aDFace->IsSet(IMeshData_ReMesh) ? 1 : 0)
                << " selfint=" << (aDFace->IsSet(IMeshData_SelfIntersectingWire) ? 1 : 0)
                << " unoriented=" << (aDFace->IsSet(IMeshData_UnorientedWire) ? 1 : 0)
                << " toofew=" << (aDFace->IsSet(IMeshData_TooFewPoints) ? 1 : 0)
                << " outdated=" << (aDFace->IsSet(IMeshData_Outdated) ? 1 : 0)
                << " open=" << aOpenWire << " tri=" << (aTri.IsNull() ? -1 : aTri->NbNodes())
                << "/" << (aTri.IsNull() ? -1 : aTri->NbTriangles()) << std::endl;
      // TEMP T-123: the wire statuses are what `initDataStructure` branches on
      // (`BRepMesh_BaseMeshAlgo.cxx:80`, `BRepMesh_NodeInsertionMeshAlgo.hxx:70`), and
      // they are invisible in the face-level flags above.
      for (int w = 0; w < aDFace->WiresNb(); ++w)
      {
        const IMeshData::IWireHandle& aDWire = aDFace->GetWire(w);
        std::cout << "CTXWIRE " << theTag << " f=" << myHits[i] << " w=" << w
                  << " edges=" << aDWire->EdgesNb()
                  << " selfint=" << (aDWire->IsSet(IMeshData_SelfIntersectingWire) ? 1 : 0)
                  << " open=" << (aDWire->IsSet(IMeshData_OpenWire) ? 1 : 0) << std::endl;
      }
    }
  }

  double           myBBox[6];
  std::vector<int> myHits;
};
} // namespace

static int runCtxDbg(const TopoDS_Shape& theShape,
                     const double        theDeflection,
                     const double        theAngle,
                     const double*       theBBox)
{
  IMeshTools_Parameters aParams;
  aParams.Deflection         = theDeflection;
  aParams.DeflectionInterior = theDeflection;
  aParams.Angle              = theAngle;
  aParams.AngleInterior      = 2.0 * theAngle;
  aParams.MinSize =
    (std::max)(IMeshTools_Parameters::RelMinSize() * theDeflection, Precision::Confusion());
  aParams.AdjustMinSize = false;
  aParams.InParallel    = false;
  aParams.Relative      = false;
  aParams.CleanModel    = false;

  std::cout << "CTXDBG params defl=" << std::setprecision(17) << aParams.Deflection
            << " deflint=" << aParams.DeflectionInterior << " angle=" << aParams.Angle
            << " angint=" << aParams.AngleInterior << " minsize=" << aParams.MinSize
            << " ivm=" << (aParams.InternalVerticesMode ? 1 : 0)
            << " ecsda=" << (aParams.EnableControlSurfaceDeflectionAllSurfaces ? 1 : 0)
            << std::endl;

  occ::handle<ZZStageCtx> aCtx = new ZZStageCtx(theBBox);
  aCtx->SetShape(theShape);
  aCtx->ChangeParameters() = aParams;

  if (!aCtx->BuildModel() || !aCtx->DiscretizeEdges() || !aCtx->HealModel()
      || !aCtx->PreProcessModel() || !aCtx->DiscretizeFaces(Message_ProgressRange()))
  {
    return 0;
  }
  aCtx->PostProcessModel();
  return 0;
}

// --- TEMP T-122 (`--rangereplay <faceIdx> <defl> <angle>`) ---------------------------
// Replays the splitter half of `BRepMesh_NodeInsertionMeshAlgo::initDataStructure` for
// one model face: `Reset` -> `AddPoint` for every point of every pcurve of wire 0 in the
// wire's edge orientation (`BRepMesh_NodeInsertionMeshAlgo.hxx:143-183`, i.e. the same
// `collectWirePoints` loop the `SeamEdgeAmplifier` uses to size the cone seam step) ->
// `AdjustRange` (`BRepMesh_DefaultRangeSplitter.cxx:45-81`) -> `IsValid` (`cxx:85-88`).
// The splitter class follows `BRepMesh_MeshAlgoFactory::GetAlgo`
// (`BRepMesh_MeshAlgoFactory.cxx:60-118`) for the face's surface type, so `valid=0` here
// is exactly the branch that makes OCCT drop a face: `initDataStructure`
// (`BRepMesh_NodeInsertionMeshAlgo.hxx:79-83`) sets `IMeshData_Failure` and returns
// false, `BRepMesh_BaseMeshAlgo::Perform` (`cxx:47-61`) then never runs
// `generateMesh`/`commitSurfaceTriangulation`, and the face keeps a null
// `Poly_Triangulation`. The pipeline stops after `PreProcessModel` because the
// `SeamEdgeAmplifier` inserts the seam points before `DiscretizeFaces` runs, so those
// points have to be in the replayed pcurves.
template <class Splitter>
static void replayRangeOn(const IMeshData::IFaceHandle& theDFace,
                          const IMeshTools_Parameters&  theParams,
                          const char*                   theName)
{
  Splitter aSplitter;
  aSplitter.Reset(theDFace, theParams);

  const IMeshData::IWireHandle& aDWire = theDFace->GetWire(0);
  double aMinU = 1e100, aMaxU = -1e100, aMinV = 1e100, aMaxV = -1e100;
  int    aPoints = 0;
  for (int e = 0; e < aDWire->EdgesNb(); ++e)
  {
    const IMeshData::IEdgeHandle&   aDEdge = aDWire->GetEdge(e);
    const IMeshData::IPCurveHandle& aPCurve =
      aDEdge->GetPCurve(theDFace.get(), aDWire->GetEdgeOrientation(e));
    for (int p = 0; p < aPCurve->ParametersNb(); ++p)
    {
      const gp_Pnt2d& aPnt = aPCurve->GetPoint(p);
      aSplitter.AddPoint(aPnt);
      aMinU = (std::min)(aMinU, aPnt.X());
      aMaxU = (std::max)(aMaxU, aPnt.X());
      aMinV = (std::min)(aMinV, aPnt.Y());
      aMaxV = (std::max)(aMaxV, aPnt.Y());
      ++aPoints;
    }
  }
  aSplitter.AdjustRange();

  const occ::handle<BRepAdaptor_Surface>& aSurf  = theDFace->GetSurface();
  const std::pair<double, double>&        aRU    = aSplitter.GetRangeU();
  const std::pair<double, double>&        aRV    = aSplitter.GetRangeV();
  const std::pair<double, double>&        aDelta = aSplitter.GetDelta();
  const std::pair<double, double>&        aTol   = aSplitter.GetToleranceUV();
  std::printf("RANGEREPLAY %s face_type=%d valid=%d pts=%d "
              "disc_u=[%.17g,%.17g] disc_v=[%.17g,%.17g] "
              "ru=[%.17g,%.17g] rv=[%.17g,%.17g] "
              "geom_u=[%.17g,%.17g] pu=%d geom_v=[%.17g,%.17g] pv=%d "
              "delta=[%.17g,%.17g] tol=[%.17g,%.17g]\n",
              theName,
              static_cast<int>(aSurf->GetType()),
              aSplitter.IsValid() ? 1 : 0,
              aPoints,
              aMinU,
              aMaxU,
              aMinV,
              aMaxV,
              aRU.first,
              aRU.second,
              aRV.first,
              aRV.second,
              aSurf->FirstUParameter(),
              aSurf->LastUParameter(),
              aSurf->IsUPeriodic() ? 1 : 0,
              aSurf->FirstVParameter(),
              aSurf->LastVParameter(),
              aSurf->IsVPeriodic() ? 1 : 0,
              aDelta.first,
              aDelta.second,
              aTol.first,
              aTol.second);
}

static int runRangeReplay(const TopoDS_Shape& theShape,
                          const int           theFaceIdx,
                          const double        theDeflection,
                          const double        theAngle)
{
  IMeshTools_Parameters aParams;
  aParams.Deflection         = theDeflection;
  aParams.DeflectionInterior = theDeflection;
  aParams.Angle              = theAngle;
  aParams.AngleInterior      = 2.0 * theAngle;
  aParams.MinSize =
    (std::max)(IMeshTools_Parameters::RelMinSize() * theDeflection, Precision::Confusion());
  aParams.AdjustMinSize = false;
  aParams.InParallel    = false;
  aParams.Relative      = false;
  aParams.CleanModel    = false;

  occ::handle<BRepMesh_Context> aCtx = new BRepMesh_Context;
  aCtx->SetShape(theShape);
  aCtx->ChangeParameters() = aParams;
  if (!aCtx->BuildModel() || !aCtx->DiscretizeEdges() || !aCtx->HealModel()
      || !aCtx->PreProcessModel())
  {
    return 0;
  }
  const occ::handle<IMeshData_Model>& aModel = aCtx->GetModel();
  if (aModel.IsNull() || theFaceIdx < 0 || theFaceIdx >= aModel->FacesNb())
  {
    std::cout << "RANGEREPLAY model-null-or-range\n";
    return 0;
  }
  const IMeshData::IFaceHandle& aDFace = aModel->GetFace(theFaceIdx);
  switch (aDFace->GetSurface()->GetType())
  {
    case GeomAbs_Cone:
      replayRangeOn<BRepMesh_ConeRangeSplitter>(aDFace, aParams, "cone");
      break;
    case GeomAbs_Cylinder:
      replayRangeOn<BRepMesh_CylinderRangeSplitter>(aDFace, aParams, "cylinder");
      break;
    case GeomAbs_Sphere:
      replayRangeOn<BRepMesh_SphereRangeSplitter>(aDFace, aParams, "sphere");
      break;
    case GeomAbs_Torus:
      replayRangeOn<BRepMesh_TorusRangeSplitter>(aDFace, aParams, "torus");
      break;
    case GeomAbs_BezierSurface:
    case GeomAbs_BSplineSurface:
      replayRangeOn<BRepMesh_NURBSRangeSplitter>(aDFace, aParams, "nurbs");
      break;
    default:
      replayRangeOn<BRepMesh_DefaultRangeSplitter>(aDFace, aParams, "default");
      break;
  }
  return 0;
}

// --- TEMP T-123 (`--algospy <defl> <angle> <faceIdx> [<faceIdx> ...]`) ---------------
// Replays `BRepMesh_FaceDiscret::process` for the given model faces with the real
// `BRepMesh_DelaunayNodeInsertionMeshAlgo<...RangeSplitter>` that
// `BRepMesh_MeshAlgoFactory::GetAlgo` (`BRepMesh_MeshAlgoFactory.cxx:60-118`) hands out
// for a cone, subclassed so the internals become observable. Two and only two branches
// raise `IMeshData_Failure` for a face whose status flags are clean before face
// discretization (`--ctxdbg` shows that): `BRepMesh_NodeInsertionMeshAlgo::
// initDataStructure` (`BRepMesh_NodeInsertionMeshAlgo.hxx:61-111`) returning false after
// `AdjustRange`/`IsValid`, and `BRepMesh_BaseMeshAlgo::commitSurfaceTriangulation`
// (`BRepMesh_BaseMeshAlgo.cxx:228-241`) finding `ElementsOfDomain()` empty - which
// `collectTriangles` (`cxx:245-252`) reports as a null triangulation. The `init` line
// tells them apart: no `init` line at all = the first branch, `elems=0` after `gen` = the
// second, and `elems` collapsing between `gen_pre`/`gen_post`/`ins_pre`/`ins_post`
// locates the loop (`BRepMesh_Delaun` base pass vs `insertNodes`) that empties it.
// `elems` is the `collectTriangles` triangle count, `links` the contour link count
// `initDataStructure` built, so the same lines also prove whether the wire loops of
// `initDataStructure` (skipped for `SelfIntersectingWire`, or `OpenWire` on non-outer
// wires) contributed anything.
namespace
{
template <class Splitter>
class ZZAlgoSpy
    : public BRepMesh_DelaunayNodeInsertionMeshAlgo<Splitter, BRepMesh_DelaunayBaseMeshAlgo>
{
  using Base =
    BRepMesh_DelaunayNodeInsertionMeshAlgo<Splitter, BRepMesh_DelaunayBaseMeshAlgo>;

protected:
  bool initDataStructure() override
  {
    const bool aRes = Base::initDataStructure();
    dump("init", aRes ? 1 : 0);
    return aRes;
  }

  void generateMesh(const Message_ProgressRange& theRange) override
  {
    dump("gen_pre", -1);
    dumpMeshInput();
    Base::generateMesh(theRange);
    dump("gen_post", -1);
  }

  void postProcessMesh(BRepMesh_Delaun& theMesher, const Message_ProgressRange& theRange) override
  {
    dump("ins_pre", -1);
    Base::postProcessMesh(theMesher, theRange);
    dump("ins_post", -1);
  }

private:
  void dump(const char* theTag, const int theRes)
  {
    const IMeshData::IFaceHandle&                      aDFace  = this->getDFace();
    const Splitter&                                    aSplit  = this->getRangeSplitter();
    const occ::handle<BRepMesh_DataStructureOfDelaun>& aStruct = this->getStructure();
    const std::pair<double, double>&                   aRU     = aSplit.GetRangeU();
    const std::pair<double, double>&                   aRV     = aSplit.GetRangeV();
    const std::pair<double, double>&                   aDelta  = aSplit.GetDelta();
    std::printf("ALGOSPY %s res=%d valid=%d nodes=%d links=%d elems=%d linkdom=%d "
                "ru=[%.17g,%.17g] rv=[%.17g,%.17g] delta=[%.17g,%.17g]\n",
                theTag,
                theRes,
                const_cast<Splitter&>(aSplit).IsValid() ? 1 : 0,
                aStruct.IsNull() ? -1 : aStruct->NbNodes(),
                aStruct.IsNull() ? -1 : aStruct->NbLinks(),
                aStruct.IsNull() ? -1 : aStruct->ElementsOfDomain().Extent(),
                aStruct.IsNull() ? -1 : aStruct->LinksOfDomain().Extent(),
                aRU.first,
                aRU.second,
                aRV.first,
                aRV.second,
                aDelta.first,
                aDelta.second);
    std::fflush(stdout);
    if (std::string(theTag) == "gen_pre")
    {
      dumpPolygon(theTag);
    }
  }

  // Dumps the two inputs `BRepMesh_DelaunayBaseMeshAlgo::generateMesh` (`cxx:34-45`)
  // derives from the structure: the vertex list length (`aVerticesOrder` = one entry per
  // registered 3D point, `cxx:38-41`) and the acceleration-grid cell counts
  // (`getCellsCount` -> `BRepMesh_GeomTool::CellsCount`, `cxx:42`), plus the node tool
  // tolerance/cell size `BRepMesh_NodeInsertionMeshAlgo::initDataStructure`
  // (`hxx:86-92`) set. The cell size decides which triangles `BRepMesh_CircleTool::Select`
  // offers to `createTrianglesOnNewVertices` for a candidate vertex, so a difference here
  // is enough to change a degenerate triangulation's output.
  void dumpMeshInput() const
  {
    const occ::handle<BRepMesh_DataStructureOfDelaun>& aStruct = this->getStructure();
    if (aStruct.IsNull())
    {
      return;
    }
    double aTolU = 0.0, aTolV = 0.0;
    aStruct->Data()->GetTolerance(aTolU, aTolV);
    const std::pair<int, int> aCells =
      const_cast<ZZAlgoSpy*>(this)->getCellsCount(this->getNodesMap()->Length());
    std::printf("ALGOSPYM nodes=%d vorder=%d cells=[%d,%d] tol=[%.17g,%.17g]\n",
                aStruct->NbNodes(),
                this->getNodesMap()->Length(),
                aCells.first,
                aCells.second,
                aTolU,
                aTolV);
    std::fflush(stdout);
  }

  // Dumps the structure `BRepMesh_Delaun` is constructed from in
  // `BRepMesh_DelaunayBaseMeshAlgo::generateMesh` (`cxx:36-47`): every node in the scaled
  // basis (`BRepMesh_DefaultRangeSplitter::Scale(thePoint, true)`) plus the links of the
  // domain, which is the boundary polygon the mesher has to triangulate.
  void dumpPolygon(const char* theTag) const
  {
    const occ::handle<BRepMesh_DataStructureOfDelaun>& aStruct = this->getStructure();
    if (aStruct.IsNull())
    {
      return;
    }
    for (int i = 1; i <= aStruct->NbNodes(); ++i)
    {
      const BRepMesh_Vertex& aV = aStruct->GetNode(i);
      std::printf("ALGOSPYN %s %d %.17g %.17g %d\n",
                  theTag,
                  i,
                  aV.Coord().X(),
                  aV.Coord().Y(),
                  static_cast<int>(aV.Movability()));
    }
    for (IMeshData::IteratorOfMapOfInteger aLinkIt(aStruct->LinksOfDomain()); aLinkIt.More();
         aLinkIt.Next())
    {
      const BRepMesh_Edge& aL = aStruct->GetLink(aLinkIt.Key());
      std::printf("ALGOSPYP %s %d %d %d %d\n",
                  theTag,
                  aLinkIt.Key(),
                  aL.FirstNode(),
                  aL.LastNode(),
                  static_cast<int>(aL.Movability()));
    }
    std::fflush(stdout);
  }
};
} // namespace

// Maps a surface type to the algo `BRepMesh_MeshAlgoFactory::GetAlgo`
// (`BRepMesh_MeshAlgoFactory.cxx:60-118`) picks for the default parameters used here
// (`InternalVerticesMode = false`, `EnableControlSurfaceDeflectionAllSurfaces = false`),
// i.e. the node-insertion algo for every type this replay supports.
static occ::handle<IMeshTools_MeshAlgo> makeSpyAlgo(const GeomAbs_SurfaceType theType)
{
  switch (theType)
  {
    case GeomAbs_Cone:
      return new ZZAlgoSpy<BRepMesh_ConeRangeSplitter>;
    default:
      return occ::handle<IMeshTools_MeshAlgo>();
  }
}

static int runAlgoSpy(const TopoDS_Shape& theShape,
                      const double        theDeflection,
                      const double        theAngle,
                      const int*          theFaces,
                      const int           theFacesNb)
{
  IMeshTools_Parameters aParams;
  aParams.Deflection         = theDeflection;
  aParams.DeflectionInterior = theDeflection;
  aParams.Angle              = theAngle;
  aParams.AngleInterior      = 2.0 * theAngle;
  aParams.MinSize =
    (std::max)(IMeshTools_Parameters::RelMinSize() * theDeflection, Precision::Confusion());
  aParams.AdjustMinSize = false;
  aParams.InParallel    = false;
  aParams.Relative      = false;
  aParams.CleanModel    = false;

  occ::handle<BRepMesh_Context> aCtx = new BRepMesh_Context;
  aCtx->SetShape(theShape);
  aCtx->ChangeParameters() = aParams;

  if (!aCtx->BuildModel() || !aCtx->DiscretizeEdges() || !aCtx->HealModel()
      || !aCtx->PreProcessModel())
  {
    std::cout << "ALGOSPY stage-failed\n";
    return 0;
  }

  const occ::handle<IMeshData_Model>& aModel = aCtx->GetModel();
  if (aModel.IsNull())
  {
    std::cout << "ALGOSPY model-null\n";
    return 0;
  }

  for (int i = 0; i < theFacesNb; ++i)
  {
    if (theFaces[i] < 0 || theFaces[i] >= aModel->FacesNb())
    {
      continue;
    }
    const IMeshData::IFaceHandle&   aDFace = aModel->GetFace(theFaces[i]);
    const GeomAbs_SurfaceType       aType  = aDFace->GetSurface()->GetType();
    std::printf("ALGOSPY face f=%d type=%d wires=%d fail=%d reusable=%d\n",
                theFaces[i],
                static_cast<int>(aType),
                aDFace->WiresNb(),
                aDFace->IsSet(IMeshData_Failure) ? 1 : 0,
                aDFace->IsSet(IMeshData_Reused) ? 1 : 0);
    for (int w = 0; w < aDFace->WiresNb(); ++w)
    {
      const IMeshData::IWireHandle& aDWire = aDFace->GetWire(w);
      std::printf("ALGOSPY wire f=%d w=%d edges=%d selfint=%d open=%d\n",
                  theFaces[i],
                  w,
                  aDWire->EdgesNb(),
                  aDWire->IsSet(IMeshData_SelfIntersectingWire) ? 1 : 0,
                  aDWire->IsSet(IMeshData_OpenWire) ? 1 : 0);
    }
    std::fflush(stdout);

    // `BRepMesh_FaceDiscret::process` (`BRepMesh_FaceDiscret.cxx:99-131`).
    if (aDFace->IsSet(IMeshData_Failure) || aDFace->IsSet(IMeshData_Reused))
    {
      std::printf("ALGOSPY skip f=%d\n", theFaces[i]);
      continue;
    }
    occ::handle<IMeshTools_MeshAlgo> aAlgo = makeSpyAlgo(aType);
    if (aAlgo.IsNull())
    {
      std::printf("ALGOSPY unsupported f=%d type=%d\n", theFaces[i], static_cast<int>(aType));
      continue;
    }
    aAlgo->Perform(aDFace, aParams, Message_ProgressRange());

    TopLoc_Location                        aLoc;
    const occ::handle<Poly_Triangulation>& aTri =
      BRep_Tool::Triangulation(aDFace->GetFace(), aLoc);
    std::printf("ALGOSPY final f=%d fail=%d tri=%d/%d\n",
                theFaces[i],
                aDFace->IsSet(IMeshData_Failure) ? 1 : 0,
                aTri.IsNull() ? -1 : aTri->NbNodes(),
                aTri.IsNull() ? -1 : aTri->NbTriangles());
    std::fflush(stdout);
  }
  return 0;
}

// --- TEMP T-103 (`--orderdbg <defl> <angle> <x1> <y1> <z1> <x2> <y2> <z2>`) -----------
// Dumps, for every face whose bbox matches, the `BRepMesh_ShapeVisitor::addWire`
// ordering stage: `ShapeExtend_WireData` list + `ShapeAnalysis_Wire::CheckOrder`
// (`ShapeAnalysis_WireOrder::Perform`) result, i.e. the exact `Ordered(i)` chain
// the port's `WireOrder` must reproduce.
static int runOrderDbg(const TopoDS_Shape& theShape, const double* theBBox)
{
  int aFound = 0;
  TopExp_Explorer aEx(theShape, TopAbs_FACE);
  for (; aEx.More(); aEx.Next())
  {
    const TopoDS_Face& aFace = TopoDS::Face(aEx.Current());
    Bnd_Box            aBox;
    BRepBndLib::Add(aFace, aBox, false);
    if (aBox.IsVoid())
    {
      continue;
    }
    double v[6] = {0, 0, 0, 0, 0, 0};
    aBox.Get(v[0], v[1], v[2], v[3], v[4], v[5]);
    if (!(std::abs(v[0] - theBBox[0]) < 1e-3 && std::abs(v[1] - theBBox[1]) < 1e-3
          && std::abs(v[2] - theBBox[2]) < 1e-3 && std::abs(v[3] - theBBox[3]) < 1e-3
          && std::abs(v[4] - theBBox[4]) < 1e-3 && std::abs(v[5] - theBBox[5]) < 1e-3))
    {
      continue;
    }
    ++aFound;
    // TEMP: print the occurrence face orientation, and repeat the whole dump for
    // the oriented-FORWARD face the real `BRepMesh_ShapeVisitor::Visit` uses.
    std::printf("ODBG faceori=%d\n", (int)aFace.Orientation());
    const TopoDS_Face aFaceFwd = TopoDS::Face(aFace.Oriented(TopAbs_FORWARD));
    for (int aPass = 0; aPass < 2; ++aPass)
    {
    const TopoDS_Face& aSrc = (aPass == 0) ? aFace : aFaceFwd;
    std::printf("ODBG pass=%d\n", aPass);
    int aWireIt = 0;
    for (TopExp_Explorer aW(aSrc, TopAbs_WIRE); aW.More(); aW.Next(), ++aWireIt)
    {
      const TopoDS_Wire& aWire     = TopoDS::Wire(aW.Current());
      occ::handle<ShapeExtend_WireData> aWireData = new ShapeExtend_WireData(aWire, true, false);
      std::printf("ODBG wire=%d stored=%d ori=%d\n", aWireIt, aWireData->NbEdges(),
                  (int)aWire.Orientation());
      for (int i = 1; i <= aWireData->NbEdges(); ++i)
      {
        const TopoDS_Edge&   anEdge = aWireData->Edge(i);
        BRepAdaptor_Curve    aC3d(TopoDS::Edge(anEdge.Oriented(TopAbs_FORWARD)));
        double               f = 0., l = 0.;
        const occ::handle<Geom2d_Curve> aC2d = BRep_Tool::CurveOnSurface(
          TopoDS::Edge(anEdge.Oriented(TopAbs_FORWARD)), aFace, f, l);
        gp_Pnt2d aB(0., 0.), aE(0., 0.);
        if (!aC2d.IsNull())
        {
          const bool aRev = (anEdge.Orientation() == TopAbs_REVERSED);
          aB = aC2d->Value(aRev ? l : f);
          aE = aC2d->Value(aRev ? f : l);
        }
        const gp_Pnt aP0 = aC3d.Value(aC3d.FirstParameter());
        const gp_Pnt aP1 = aC3d.Value(aC3d.LastParameter());
        const char*  aKind = aC2d.IsNull() ? "null" : aC2d->DynamicType()->Name();
        const bool   aPer  = !aC2d.IsNull() && aC2d->IsPeriodic();
        std::printf("ODBG   i=%d edge=%p ori=%d closed=%d deg=%d uv=(%.9f,%.9f)->(%.9f,%.9f) "
                    "p3d=(%.6f,%.6f,%.6f)->(%.6f,%.6f,%.6f)\n",
                    i,
                    (void*)anEdge.TShape().get(),
                    (int)anEdge.Orientation(),
                    BRep_Tool::IsClosed(anEdge, aFace) ? 1 : 0,
                    BRep_Tool::Degenerated(anEdge) ? 1 : 0,
                    aB.X(),
                    aB.Y(),
                    aE.X(),
                    aE.Y(),
                    aP0.X(),
                    aP0.Y(),
                    aP0.Z(),
                    aP1.X(),
                    aP1.Y(),
                    aP1.Z());
        std::printf("ODBG     c2d=%s per=%d range2d=(%.9f,%.9f) range3d=(%.9f,%.9f) sameParam=%d "
                    "sameRange=%d v(f)=(%.9f,%.9f) v(l)=(%.9f,%.9f)\n",
                    aKind,
                    aPer ? 1 : 0,
                    f,
                    l,
                    aC3d.FirstParameter(),
                    aC3d.LastParameter(),
                    BRep_Tool::SameParameter(TopoDS::Edge(anEdge.Oriented(TopAbs_FORWARD))) ? 1 : 0,
                    BRep_Tool::SameRange(TopoDS::Edge(anEdge.Oriented(TopAbs_FORWARD))) ? 1 : 0,
                    aC2d.IsNull() ? 0.0 : aC2d->Value(f).X(),
                    aC2d.IsNull() ? 0.0 : aC2d->Value(f).Y(),
                    aC2d.IsNull() ? 0.0 : aC2d->Value(l).X(),
                    aC2d.IsNull() ? 0.0 : aC2d->Value(l).Y());
      }
      ShapeAnalysis_Wire       aWireTool(aWireData, aFace, Precision::Confusion());
      ShapeAnalysis_WireOrder  aOrderTool;
      aWireTool.CheckOrder(aOrderTool, true, false);
      std::printf("ODBG wire=%d status: fail=%d done3=%d nb=%d ordered=[",
                  aWireIt,
                  aWireTool.LastCheckStatus(ShapeExtend_FAIL) ? 1 : 0,
                  aWireTool.LastCheckStatus(ShapeExtend_DONE3) ? 1 : 0,
                  aOrderTool.NbEdges());
      for (int i = 1; i <= aOrderTool.NbEdges(); ++i)
      {
        std::printf("%s%d", i > 1 ? " " : "", aOrderTool.Ordered(i));
      }
      std::printf("]\n");
    }
    }
    break;
  }
  std::printf("ODBG found=%d\n", aFound);
  return 0;
}

// --- TEMP T-103 (`--wiredbg <defl> <angle> <x1> <y1> <z1> <x2> <y2> <z2>`) ------------
// Ground truth for the port's `CWDBG`/`CBWIN`. Replays the real `IMeshTools_MeshBuilder`
// prefix (`BuildModel` -> `DiscretizeEdges` -> `HealModel`) and dumps, for every wire
// slot of the target face, the `IMeshData_Wire` order: the slot orientation, the
// pcurve orientation / `IsForward`, `SameParam` / `SameRange`, both parameter ranges
// and the full 3D + 2D point arrays in `GetPoint(i)` order (the arrays the classifier
// polygon and the Delaunay frontier are built from).
static int runWireDbg(const TopoDS_Shape& theShape,
                      const double        theDeflection,
                      const double        theAngle,
                      const double*       theBBox,
                      const int           theFaceIdx = -1)
{
  IMeshTools_Parameters aParams;
  aParams.Deflection    = theDeflection;
  aParams.Angle         = theAngle;
  aParams.InParallel    = false;
  aParams.Relative      = false;
  aParams.MinSize       = (std::max)(IMeshTools_Parameters::RelMinSize() * theDeflection,
                                     Precision::Confusion());
  aParams.AdjustMinSize = false;

  occ::handle<BRepMesh_ModelBuilder> aBuilder = new BRepMesh_ModelBuilder;
  const occ::handle<IMeshData_Model> aModel   = aBuilder->Perform(theShape, aParams);
  if (aModel.IsNull())
  {
    std::printf("WDBG model-null\n");
    return 0;
  }
  // `IMeshTools_MeshBuilder::Perform` (`cxx:52-56`): DiscretizeEdges runs BEFORE
  // HealModel, so the pre-heal arrays below are the EdgeDiscret output.
  occ::handle<BRepMesh_EdgeDiscret> anEdgeDiscret = new BRepMesh_EdgeDiscret;
  anEdgeDiscret->Perform(aModel, aParams, Message_ProgressRange());
  occ::handle<BRepMesh_ModelHealer> aHealer = new BRepMesh_ModelHealer;
  aHealer->Perform(aModel, aParams, Message_ProgressRange());

  for (int f = 0; f < aModel->FacesNb(); ++f)
  {
    const IMeshData::IFaceHandle& aDFace = aModel->GetFace(f);
    const TopoDS_Face&            aFace  = aDFace->GetFace();
    Bnd_Box                       aBox;
    BRepBndLib::Add(aFace, aBox, false);
    if (aBox.IsVoid())
    {
      continue;
    }
    double v[6] = {0, 0, 0, 0, 0, 0};
    aBox.Get(v[0], v[1], v[2], v[3], v[4], v[5]);
    // Selection: by `IMeshData_Model::GetFace` index when the caller gave one
    // (the port's own bbox is tolerance-inflated, so the six-coordinate form
    // cannot pair the two sides), otherwise by the exact geometry box.
    if (theFaceIdx >= 0)
    {
      if (f != theFaceIdx)
      {
        continue;
      }
    }
    else if (!(std::abs(v[0] - theBBox[0]) < 1e-3 && std::abs(v[1] - theBBox[1]) < 1e-3
               && std::abs(v[2] - theBBox[2]) < 1e-3 && std::abs(v[3] - theBBox[3]) < 1e-3
               && std::abs(v[4] - theBBox[4]) < 1e-3 && std::abs(v[5] - theBBox[5]) < 1e-3))
    {
      continue;
    }
    std::printf("WDBG face=%d bbox=(%.9f,%.9f,%.9f)-(%.9f,%.9f,%.9f)\n", f, v[0], v[1], v[2], v[3],
                v[4], v[5]);
    std::printf("WDBG face=%d surfType=%d wires=%d\n", f, (int)aDFace->GetSurface()->GetType(),
                aDFace->WiresNb());
    for (int w = 0; w < aDFace->WiresNb(); ++w)
    {
      const IMeshData::IWireHandle& aDWire = aDFace->GetWire(w);
      std::printf("WDBG  w=%d edges=%d\n", w, aDWire->EdgesNb());
      for (int e = 0; e < aDWire->EdgesNb(); ++e)
      {
        const IMeshData::IEdgeHandle& aDEdge  = aDWire->GetEdge(e);
        const TopAbs_Orientation      aSlotOri = aDWire->GetEdgeOrientation(e);
        const IMeshData::IPCurveHandle& aPCurve = aDEdge->GetPCurve(aDFace.get(), aSlotOri);
        const TopoDS_Edge&              aEdge   = aDEdge->GetEdge();
        double                          r3f = 0., r3l = 0.;
        BRep_Tool::Range(TopoDS::Edge(aEdge.Oriented(TopAbs_FORWARD)), r3f, r3l);
        std::printf("WDBG   e=%d ptr=%p slotOri=%d pcsNb=%d pcOri=%d fwd=%d sameParam=%d "
                    "sameRange=%d deg=%d range3d=(%.9f,%.9f) nb3d=%d nb2d=%d\n",
                    e,
                    (void*)aEdge.TShape().get(),
                    (int)aSlotOri,
                    aDEdge->PCurvesNb(),
                    (int)aPCurve->GetOrientation(),
                    aPCurve->IsForward() ? 1 : 0,
                    aDEdge->GetSameParam() ? 1 : 0,
                    aDEdge->GetSameRange() ? 1 : 0,
                    aDEdge->GetDegenerated() ? 1 : 0,
                    r3f,
                    r3l,
                    aDEdge->GetCurve()->ParametersNb(),
                    aPCurve->ParametersNb());
        std::printf("WDBG   p2d=[");
        for (int k = 0; k < aPCurve->ParametersNb(); ++k)
        {
          std::printf("%s(%.17g,%.17g)", k ? " " : "", aPCurve->GetPoint(k).X(),
                      aPCurve->GetPoint(k).Y());
        }
        std::printf("]\n");
        std::printf("WDBG   p3d=[");
        for (int k = 0; k < aDEdge->GetCurve()->ParametersNb(); ++k)
        {
          const gp_Pnt& p = aDEdge->GetCurve()->GetPoint(k);
          std::printf("%s(%.6f,%.6f,%.6f)", k ? " " : "", p.X(), p.Y(), p.Z());
        }
        std::printf("]\n");
      }
    }
  }
  return 0;
}

// --- T-118 IGES oracle -------------------------------------------------------------
// --iges <out.iges> writes the transferred shape with OCCT's own
// IGESControl_Writer in the settings the Rust writer ports: IGESControl_Controller
// initialised, Faces mode (write.iges.brep.mode = 0, IGESData.cxx:90-94) and
// millimetres. The result is the reference GeomToIGES_* / BRepToIGES_BRShell
// output, so the port's entity set and every parameter block can be compared
// against it instead of against a previous Rust run.
static int runIges(const TopoDS_Shape& theShape, const char* thePath)
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
static int runIgesPc(const TopoDS_Shape& theShape, const char* thePath)
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

int main(int argc, char** argv)
{
  if (argc < 2)
  {
    std::cerr << "usage: occt_probe <file.step> [--frames] [--nofix]\n";
    return 2;
  }
  const bool frames = (argc > 2);
  // TEMP T-104: face index for `--structdump` (dispatched once the shape is read).
  int aStructDumpFace = -1;
  // TEMP T-121: optional exact `<defl> <angle>` override for `--structdump` /
  // `--dcstruct`, so both sides can be driven with the port's own numbers.
  double aStructDumpDefl = 0.0;
  double aStructDumpAngle = 0.0;
  // TEMP T-119: face index for `--dcstruct` (factory-faithful splitter + DC pass).
  int aDCStructFace = -1;
  // TEMP T-120: face index for `--srfnodes` (factory-faithful splitter + classifier).
  int aSrfNodeFace = -1;
  // TEMP T-99: `--delauncheck` proves BRepMesh_Delaun is reachable from this probe.
  for (int i = 2; i < argc; ++i)
  {
    if (std::string(argv[i]) == "--delauncheck")
    {
      return runDelaunCheck();
    }
    if (std::string(argv[i]) == "--delaunfeed" && i + 1 < argc)
    {
      return runDelaunFeed(argv[i + 1]);
    }
    if (std::string(argv[i]) == "--structdump" && i + 1 < argc)
    {
      aStructDumpFace = std::atoi(argv[i + 1]);
      if (i + 3 < argc)
      {
        aStructDumpDefl  = std::atof(argv[i + 2]);
        aStructDumpAngle = std::atof(argv[i + 3]) * M_PI / 180.0;
      }
    }
    // TEMP T-119: same dump, but through the splitter/algo pair the real
    // `BRepMesh_MeshAlgoFactory` picks (so a `GeomAbs_SurfaceOfRevolution` face is
    // meshed by `BoundaryParamsRangeSplitter` + the DC pass, not `DefaultRangeSplitter`).
    if (std::string(argv[i]) == "--dcstruct" && i + 1 < argc)
    {
      aDCStructFace = std::atoi(argv[i + 1]);
      if (i + 3 < argc)
      {
        aStructDumpDefl  = std::atof(argv[i + 2]);
        aStructDumpAngle = std::atof(argv[i + 3]) * M_PI / 180.0;
      }
    }
    // TEMP T-120: `--srfnodes <faceIdx>` prints OCCT's surface-node classifier
    // verdicts for the factory-selected splitter (see `ZZSrfNodeCountAlgo`).
    if (std::string(argv[i]) == "--srfnodes" && i + 1 < argc)
    {
      aSrfNodeFace = std::atoi(argv[i + 1]);
    }
    if (std::string(argv[i]) == "--delaunstruct" && i + 1 < argc)
    {
      return runDelaunStruct(argv[i + 1]);
    }
    // TEMP T-99: per-vertex variant for step-by-step diffs (`--delaunsteps <dir> <face>`).
    if (std::string(argv[i]) == "--delaunsteps" && i + 2 < argc)
    {
      return runDelaunSteps(argv[i + 1], std::atoi(argv[i + 2]));
    }
    // TEMP T-121: per-insertion verdict dump
    // (`--delaunwhy <dir> <face> [step]`, `<step>` < 0 = every step).
    if (std::string(argv[i]) == "--delaunwhy" && i + 2 < argc)
    {
      const int aStep = (i + 3 < argc) ? std::atoi(argv[i + 3]) : -1;
      return runDelaunWhy(argv[i + 1], std::atoi(argv[i + 2]), aStep);
    }
  }
  bool       nofix  = false;
  bool       faces0 = false;
  for (int i = 2; i < argc; ++i)
  {
    if (std::string(argv[i]) == "--nofix")
      nofix = true;
    if (std::string(argv[i]) == "--faces0")
      faces0 = true;
  }

  // TEMP T-102: `--fix <field> <int>` overrides one `DE_ShapeFixParameters`
  // field (applied after `ReadFile`, see below). `--set <static> <ival>`
  // overrides one `Interface_Static` (also after `ReadFile`: the controller
  // registers the `FromSTEP.FixShape.*` statics lazily as type 't', and the
  // `ShapeProcess_ShapeContext` only reads them while the FixShape operator runs
  // inside TransferRoots). See also `wires_probe.cpp` (`set`/`setc`/`setp`).
  DE_ShapeFixParameters aFixParams;
  bool                 anyFixParam = false;
  {
    auto setFixMode = [&](const std::string& n, int v) -> bool {
      const DE_ShapeFixParameters::FixMode m =
        static_cast<DE_ShapeFixParameters::FixMode>(v);
      if (false)
      {
      }
#define FIXMODE(N) else if (n == #N) { aFixParams.N = m; }
      FIXMODE(FixFreeShellMode)
      FIXMODE(FixFreeFaceMode)
      FIXMODE(FixFreeWireMode)
      FIXMODE(FixSameParameterMode)
      FIXMODE(FixSolidMode)
      FIXMODE(FixShellOrientationMode)
      FIXMODE(CreateOpenSolidMode)
      FIXMODE(FixShellMode)
      FIXMODE(FixFaceOrientationMode)
      FIXMODE(FixFaceMode)
      FIXMODE(FixWireMode)
      FIXMODE(FixOrientationMode)
      FIXMODE(FixAddNaturalBoundMode)
      FIXMODE(FixMissingSeamMode)
      FIXMODE(FixSmallAreaWireMode)
      FIXMODE(RemoveSmallAreaFaceMode)
      FIXMODE(FixIntersectingWiresMode)
      FIXMODE(FixLoopWiresMode)
      FIXMODE(FixSplitFaceMode)
      FIXMODE(AutoCorrectPrecisionMode)
      FIXMODE(ModifyTopologyMode)
      FIXMODE(ModifyGeometryMode)
      FIXMODE(ClosedWireMode)
      FIXMODE(PreferencePCurveMode)
      FIXMODE(FixReorderMode)
      FIXMODE(FixSmallMode)
      FIXMODE(FixConnectedMode)
      FIXMODE(FixEdgeCurvesMode)
      FIXMODE(FixDegeneratedMode)
      FIXMODE(FixLackingMode)
      FIXMODE(FixSelfIntersectionMode)
      FIXMODE(RemoveLoopMode)
      FIXMODE(FixReversed2dMode)
      FIXMODE(FixRemovePCurveMode)
      FIXMODE(FixRemoveCurve3dMode)
      FIXMODE(FixAddPCurveMode)
      FIXMODE(FixAddCurve3dMode)
      FIXMODE(FixSeamMode)
      FIXMODE(FixShiftedMode)
      FIXMODE(FixEdgeSameParameterMode)
      FIXMODE(FixNotchedEdgesMode)
      FIXMODE(FixTailMode)
      FIXMODE(MaxTailAngle)
      FIXMODE(MaxTailWidth)
      FIXMODE(FixSelfIntersectingEdgeMode)
      FIXMODE(FixIntersectingEdgesMode)
      FIXMODE(FixNonAdjacentIntersectingEdgesMode)
      FIXMODE(FixVertexPositionMode)
      FIXMODE(FixVertexToleranceMode)
#undef FIXMODE
      else
      {
        return false;
      }
      return true;
    };
    for (int i = 2; i + 2 < argc; ++i)
    {
      if (std::string(argv[i]) == "--fix")
      {
        if (!setFixMode(argv[i + 1], std::atoi(argv[i + 2])))
        {
          std::cerr << "unknown fix parameter: " << argv[i + 1] << "\n";
          return 5;
        }
        anyFixParam = true;
        std::cout << "FIX " << argv[i + 1] << " = " << argv[i + 2] << "\n";
      }
    }
  }

  // TEMP T-119: `.iges`/`.igs` inputs are read through `IGESControl_Reader`
  // (`IGESControl_Controller::Init` + `ReadFile` + `TransferRoots` + `OneShape`,
  // the same entry chain the port's `iges_read` mirrors), so the pipeline modes
  // below (`--structdump`, `--mesh`, ...) can be pointed at an IGES oracle.
  const std::string aReadPath(argv[1]);
  const bool         isIges =
    (aReadPath.size() >= 5 && aReadPath.compare(aReadPath.size() - 5, 5, ".iges") == 0)
    || (aReadPath.size() >= 4 && aReadPath.compare(aReadPath.size() - 4, 4, ".igs") == 0);

  TopoDS_Shape aShape;
  // Declared for the STEP-only modes below (`--faceids`, `--entities`, ...); the
  // IGES branch leaves it default-constructed, so those modes remain STEP-only.
  STEPControl_Reader aReader;
  if (isIges)
  {
    IGESControl_Controller::Init();
    IGESControl_Reader aIgesReader;
    if (faces0)
    {
      aFixParams.FixFaceOrientationMode =
        static_cast<DE_ShapeFixParameters::FixMode>(0); // NeedFix(0) => do not fix
      anyFixParam = true;
      std::cout << "SHP-FIX face-orientation disabled\n";
    }
    if (aIgesReader.ReadFile(argv[1]) != IFSelect_RetDone)
    {
      std::cerr << "read failed: " << argv[1] << "\n";
      return 3;
    }
    // Same delegation caveat as the STEP branch below: the reader's transfer
    // actor only exists after `ReadFile`.
    for (int i = 2; i + 2 < argc; ++i)
    {
      if (std::string(argv[i]) == "--set")
      {
        Interface_Static::SetIVal(argv[i + 1], std::atoi(argv[i + 2]));
        std::cout << "SET " << argv[i + 1] << " = " << Interface_Static::IVal(argv[i + 1]) << "\n";
      }
    }
    if (anyFixParam)
    {
      aIgesReader.SetShapeFixParameters(aFixParams);
    }
    if (nofix)
    {
      ShapeProcess::OperationsFlags aFlags;
      aFlags.reset();
      aIgesReader.SetShapeProcessFlags(aFlags);
      std::cout << "SHAPE-PROCESS disabled\n";
    }
    aIgesReader.TransferRoots();
    aShape = aIgesReader.OneShape();
  }
  else
  {
    if (nofix)
    {
      ShapeProcess::OperationsFlags aFlags;
      aFlags.reset();
      aReader.SetShapeProcessFlags(aFlags);
      std::cout << "SHAPE-PROCESS disabled\n";
    }
    if (faces0)
    {
      aFixParams.FixFaceOrientationMode =
        static_cast<DE_ShapeFixParameters::FixMode>(0); // NeedFix(0) => do not fix
      anyFixParam = true;
      std::cout << "SHP-FIX face-orientation disabled\n";
    }
    if (aReader.ReadFile(argv[1]) != IFSelect_RetDone)
    {
      std::cerr << "read failed: " << argv[1] << "\n";
      return 3;
    }
    // NOTE: `SetShapeProcessFlags` / `SetShapeFixParameters` delegate to the
    // transfer actor, which only exists after `ReadFile`; calling them earlier is
    // a no-op (see `XSControl_Reader::GetActor`), and the same holds for the
    // `FromSTEP.FixShape.*` statics.
    for (int i = 2; i + 2 < argc; ++i)
    {
      if (std::string(argv[i]) == "--set")
      {
        Interface_Static::SetIVal(argv[i + 1], std::atoi(argv[i + 2]));
        std::cout << "SET " << argv[i + 1] << " = " << Interface_Static::IVal(argv[i + 1]) << "\n";
      }
    }
    if (anyFixParam)
    {
      aReader.SetShapeFixParameters(aFixParams);
    }
    if (nofix)
    {
      ShapeProcess::OperationsFlags aFlags;
      aFlags.reset();
      aReader.SetShapeProcessFlags(aFlags);
      DE_ShapeFixParameters aFix;
      aFix.FixShellOrientationMode = static_cast<DE_ShapeFixParameters::FixMode>(0);
      aFix.FixFaceOrientationMode = static_cast<DE_ShapeFixParameters::FixMode>(0);
      aFix.FixOrientationMode = static_cast<DE_ShapeFixParameters::FixMode>(0);
      aReader.SetShapeFixParameters(aFix);
      std::cout << "PROC-OFF (after read, before transfer)\n";
    }
    aReader.TransferRoots();
    aShape = aReader.OneShape();
  }
  if (aShape.IsNull())
  {
    std::cerr << "null shape\n";
    return 4;
  }
  // TEMP T-99: `--boundary <dir>` dumps the real pipeline's Delaunay boundary UV.
  for (int i = 2; i + 1 < argc; ++i)
  {
    // T-118: `--iges <out.iges>` writes the shape with OCCT's own IGESControl_Writer.
    if (std::string(argv[i]) == "--iges")
    {
      return runIges(aShape, argv[i + 1]);
    }
    // T-111: `--igespc <out.txt>` dumps the raw face/edge pcurve state seen by the
    // IGES writer, for comparing against the port's `curve_on_surface` results.
    if (std::string(argv[i]) == "--igespc")
    {
      return runIgesPc(aShape, argv[i + 1]);
    }
    if (std::string(argv[i]) == "--boundary")
    {
      return runBoundary(aShape, argv[i + 1]);
    }
    // TEMP T-133: `--lens <faceIdx> <u0> <u1> <v0> <v1>` replays the range
    // splitter's length sampling (see `runLenSamples`).
    if (std::string(argv[i]) == "--lens" && i + 5 < argc)
    {
      return runLenSamples(aShape,
                           std::atoi(argv[i + 1]),
                           std::atof(argv[i + 2]),
                           std::atof(argv[i + 3]),
                           std::atof(argv[i + 4]),
                           std::atof(argv[i + 5]));
    }
    // TEMP T-101: `--facestats <defl> <angle>` -- per-model-face nodes/triangles/bbox.
    if (std::string(argv[i]) == "--facestats" && i + 2 < argc)
    {
      return runFaceStats(aShape, std::atof(argv[i + 1]), std::atof(argv[i + 2]));
    }
    // TEMP T-101b: `--facedump <faceIdx> <defl> <angle>` -- one model face's mesh arrays.
    if (std::string(argv[i]) == "--facedump" && i + 3 < argc)
    {
      return runFaceDump(aShape,
                         std::atoi(argv[i + 1]),
                         std::atof(argv[i + 2]),
                         std::atof(argv[i + 3]));
    }
    // TEMP T-132: `--srparams <faceIdx> <defl> <angle>` -- the raw splitter parameters
    // behind one model face's generated grid.
    if (std::string(argv[i]) == "--srparams" && i + 3 < argc)
    {
      return runSrfParams(aShape,
                          std::atoi(argv[i + 1]),
                          std::atof(argv[i + 2]),
                          std::atof(argv[i + 3]));
    }
    // TEMP T-121: `--edgedu <faceIdx> <defl> <angle>` -- `EvaluateDu` inputs per edge.
    if (std::string(argv[i]) == "--edgedu" && i + 3 < argc)
    {
      return runEdgeDu(aShape,
                       std::atoi(argv[i + 1]),
                       std::atof(argv[i + 2]),
                       std::atof(argv[i + 3]));
    }
    // TEMP T-103: `--torusdbg <defl> <angle>` -- replay OCCT's torus range
    // splitter per torus face of the discrete model (see _dbg/ZZ_TorusRangeSplitter.cxx).
    if (std::string(argv[i]) == "--torusdbg" && i + 2 < argc)
    {
      return runTorusDbg(aShape, std::atof(argv[i + 1]), std::atof(argv[i + 2]));
    }
    // TEMP T-103: `--healerdbg <defl> <angle> <x1> <y1> <z1> <x2> <y2> <z2>` -- dump the
    // per-edge pcurve state of the matching face(s) before/after the model healer.
    if (std::string(argv[i]) == "--healerdbg" && i + 8 < argc)
    {
      double aBB[6];
      for (int k = 0; k < 6; ++k)
      {
        aBB[k] = std::atof(argv[i + 3 + k]);
      }
      return runHealerDbg(aShape, std::atof(argv[i + 1]), std::atof(argv[i + 2]), aBB);
    }
    // TEMP T-122: `--ctxdbg <defl> <angle> <x1> <y1> <z1> <x2> <y2> <z2>` -- per-stage
    // status flags of the matching face(s) through the real mesh pipeline.
    if (std::string(argv[i]) == "--ctxdbg" && i + 8 < argc)
    {
      double aBB[6];
      for (int k = 0; k < 6; ++k)
      {
        aBB[k] = std::atof(argv[i + 3 + k]);
      }
      return runCtxDbg(aShape, std::atof(argv[i + 1]), std::atof(argv[i + 2]), aBB);
    }
    // TEMP T-122: `--rangereplay <faceIdx> <defl> <angle>` -- replay the splitter half of
    // `initDataStructure` (AdjustRange/IsValid) for one model face.
    if (std::string(argv[i]) == "--rangereplay" && i + 3 < argc)
    {
      return runRangeReplay(aShape, std::atoi(argv[i + 1]), std::atof(argv[i + 2]),
                            std::atof(argv[i + 3]));
    }
    // TEMP T-123: `--algospy <defl> <angle> <faceIdx> [...]` -- replay the real
    // cone face-discretization algo (initDataStructure/Delaunay/insertNodes).
    if (std::string(argv[i]) == "--algospy" && i + 3 < argc)
    {
      std::vector<int> aFaces;
      for (int k = i + 3; k < argc; ++k)
      {
        aFaces.push_back(std::atoi(argv[k]));
      }
      return runAlgoSpy(aShape, std::atof(argv[i + 1]), std::atof(argv[i + 2]), aFaces.data(),
                        static_cast<int>(aFaces.size()));
    }
    // TEMP T-103: `--orderdbg <x1> <y1> <z1> <x2> <y2> <z2>` -- dump the wire ordering
    // (`ShapeAnalysis_Wire::CheckOrder`) for the matching face.
    if (std::string(argv[i]) == "--orderdbg" && i + 6 < argc)
    {
      double aBB[6];
      for (int k = 0; k < 6; ++k)
      {
        aBB[k] = std::atof(argv[i + 1 + k]);
      }
      return runOrderDbg(aShape, aBB);
    }
    // TEMP T-103: `--wiredbg <defl> <angle> <x1> <y1> <z1> <x2> <y2> <z2>` -- ground
    // truth for the port's `CWDBG`/`CBWIN` (per-slot pcurve point arrays after
    // DiscretizeEdges + HealModel).
    if (std::string(argv[i]) == "--wiredbg" && i + 8 < argc)
    {
      double aBB[6];
      for (int k = 0; k < 6; ++k)
      {
        aBB[k] = std::atof(argv[i + 3 + k]);
      }
      return runWireDbg(aShape, std::atof(argv[i + 1]), std::atof(argv[i + 2]), aBB);
    }
    // TEMP: `--wiredbgi <faceIdx> <defl> <angle>` -- same dump as `--wiredbg` but
    // selects the face by `IMeshData_Model::GetFace` index, because the port's own
    // bbox is tolerance-inflated and cannot be matched against the geometry box.
    if (std::string(argv[i]) == "--wiredbgi" && i + 3 < argc)
    {
      const double aNoBox[6] = {0, 0, 0, 0, 0, 0};
      return runWireDbg(aShape, std::atof(argv[i + 2]), std::atof(argv[i + 3]), aNoBox,
                        std::atoi(argv[i + 1]));
    }
  }

  // T-54 oracle: OCCT's triangulation counts at a given deflection
  // (`BRepMesh_IncrementalMesh`), per face and in total. The port's
  // `export_data_obj` gate compares its own `v=`/`f=` counts, so this is
  // what a faithful constrained-Delaunay port must reproduce.
  if (argc > 2 && std::string(argv[2]) == "--mesh")
  {
    const double aDeflection = (argc > 3 && argv[3][0] != '-') ? std::atof(argv[3]) : 0.1;
    const double anAngle      = (argc > 4 && argv[4][0] != '-') ? std::atof(argv[4]) : 0.5;
    BRepMesh_IncrementalMesh aMesher(aShape, aDeflection, false, anAngle);
    TopExp_Explorer            anEx(aShape, TopAbs_FACE);
    int                        anIdx = 0, aNodes = 0, aTris = 0, aNegTris = 0;
    double                     aMeshVol = 0.0;
    for (; anEx.More(); anEx.Next(), ++anIdx)
    {
      TopLoc_Location                     aLoc;
      const occ::handle<Poly_Triangulation>& aTri =
        BRep_Tool::Triangulation(TopoDS::Face(anEx.Current()), aLoc);
      const int aN = aTri.IsNull() ? 0 : aTri->NbNodes();
      const int aT = aTri.IsNull() ? 0 : aTri->NbTriangles();
      // T-99: also emit the face's bbox so this per-face triangle count can be
      // paired with the port's `PORTID` line by the six bbox coordinates (D21;
      // bbox-centre nearest neighbour is known to mis-pair). Same helper as the
      // `UVSUM` path below (`BRepBndLib::Add`).
      Bnd_Box     aFaceBox;
      BRepBndLib::Add(TopoDS::Face(anEx.Current()), aFaceBox);
      Standard_Real fw[6];
      aFaceBox.Get(fw[0], fw[1], fw[2], fw[3], fw[4], fw[5]);
      // T-99: `--mesh` (TopExp_Explorer) and `--uvsum` (IMeshData_Model::GetFace)
      // index faces in DIFFERENT orders - measured, their bboxes agree at only
      // 1 of 226 indices - so a triangle count from here cannot be attached to a
      // face by index. Print the face's `TShape` address instead: it is unique
      // per face instance and is exactly what the port's `PORTID` prints as
      // `key=` (round-trip check: the OCCT run and the port run are separate
      // processes, so compare the two address SETS only if both are stable; the
      // bbox remains the portable cross-process key).
      std::cout << "FACE " << anIdx << " nodes=" << aN << " triangles=" << aT
                << " tshape=" << TopoDS::Face(anEx.Current()).TShape().get() << std::fixed
                << std::setprecision(9) << " bbox=(" << fw[0] << "," << fw[1] << "," << fw[2]
                << ")-(" << fw[3] << "," << fw[4] << "," << fw[5] << ")\n";
      aNodes += aN;
      aTris += aT;
      // Signed mesh volume of this face's triangles w.r.t. the origin -- the
      // quantity T-87 compares (divergence theorem over the triangulation).
      if (!aTri.IsNull())
      {
        for (int t = 1; t <= aTri->NbTriangles(); ++t)
        {
          int n1, n2, n3;
          aTri->Triangle(t).Get(n1, n2, n3);
          const gp_Pnt a = aTri->Node(n1).Transformed(aLoc.Transformation());
          const gp_Pnt b = aTri->Node(n2).Transformed(aLoc.Transformation());
          const gp_Pnt c = aTri->Node(n3).Transformed(aLoc.Transformation());
          const double v = a.XYZ().Dot(b.XYZ().Crossed(c.XYZ())) / 6.0;
          aMeshVol += v;
          if (v < 0.0)
          {
            ++aNegTris;
          }
        }
      }
    }
    std::cout << "TOTAL faces=" << anIdx << " nodes=" << aNodes << " triangles=" << aTris
              << " meshvol=" << aMeshVol << " neg_triangles=" << aNegTris
              << " deflection=" << aDeflection << " angle=" << anAngle << "\n";
    return 0;
  }

  // Decisive mesh dump (`--dumpmesh <out.obj> [defl] [angle]`): write OCCT's own
  // `Poly_Triangulation` after `BRepMesh_IncrementalMesh` with the port's
  // parameters, so a port-side OBJ can be compared node by node instead of by
  // counts and scalars alone. `--mesh` prints the same numbers but no geometry:
  // a box that includes the triangulation (`useTriangulation = true`, the
  // default) can disagree with one that does not (`BRepBndLib::Add(face, box,
  // false)`), and that disagreement has to be attributed to the nodes rather
  // than guessed at.
  if (argc > 3 && std::string(argv[2]) == "--dumpmesh")
  {
    const double aDeflection = (argc > 4 && argv[4][0] != '-') ? std::atof(argv[4]) : 0.1;
    const double anAngle      = (argc > 5 && argv[5][0] != '-') ? std::atof(argv[5]) : 0.5;
    BRepMesh_IncrementalMesh aMesher(aShape, aDeflection, false, anAngle);
    std::ofstream aOut(argv[3]);
    int           aBase = 0, aFaces = 0, aAllNodes = 0, aAllTris = 0;
    double        aMin[3] = {1e300, 1e300, 1e300};
    double        aMax[3] = {-1e300, -1e300, -1e300};
    int           aIdx = 0;
    for (TopExp_Explorer anEx(aShape, TopAbs_FACE); anEx.More(); anEx.Next(), ++aIdx)
    {
      TopLoc_Location aLoc;
      const occ::handle<Poly_Triangulation>& aTri =
        BRep_Tool::Triangulation(TopoDS::Face(anEx.Current()), aLoc);
      ++aFaces;
      if (aTri.IsNull())
      {
        std::cout << "DUMPFACE i=" << aIdx << " nodes=0 triangles=0 null\n";
        continue;
      }
      const gp_Trsf& aT = aLoc.Transformation();
      for (int i = 1; i <= aTri->NbNodes(); ++i)
      {
        const gp_Pnt p = aTri->Node(i).Transformed(aT);
        aOut << "v " << std::setprecision(9) << p.X() << " " << p.Y() << " " << p.Z() << "\n";
        aMin[0]       = std::min(aMin[0], p.X());
        aMin[1]       = std::min(aMin[1], p.Y());
        aMin[2]       = std::min(aMin[2], p.Z());
        aMax[0]       = std::max(aMax[0], p.X());
        aMax[1]       = std::max(aMax[1], p.Y());
        aMax[2]       = std::max(aMax[2], p.Z());
      }
      for (int t = 1; t <= aTri->NbTriangles(); ++t)
      {
        int n1, n2, n3;
        aTri->Triangle(t).Get(n1, n2, n3);
        aOut << "f " << aBase + n1 << " " << aBase + n2 << " " << aBase + n3 << "\n";
      }
      std::cout << "DUMPFACE i=" << aIdx << " nodes=" << aTri->NbNodes()
                << " triangles=" << aTri->NbTriangles() << " base=" << aBase << "\n";
      aBase += aTri->NbNodes();
      aAllNodes += aTri->NbNodes();
      aAllTris += aTri->NbTriangles();
    }
    std::cout << std::fixed << std::setprecision(6) << "DUMPMESH faces=" << aFaces
              << " nodes=" << aAllNodes << " triangles=" << aAllTris << " bbox=(" << aMin[0] << ","
              << aMin[1] << "," << aMin[2] << ")-(" << aMax[0] << "," << aMax[1] << "," << aMax[2]
              << ") file=" << argv[3] << "\n";
    return 0;
  }


  // T-93 oracle: per-face wire/edge inventory of the *imported* shape, in the
  // same enumeration order the port uses (TopExp_Explorer over faces). Used to
  // diff the port's healed loop structure against OCCT's face by face.
  if (argc > 2 && std::string(argv[2]) == "--wireinv")
  {
    int k = 0;
    for (TopExp_Explorer ex(aShape, TopAbs_FACE); ex.More(); ex.Next(), ++k)
    {
      const TopoDS_Face aFace = TopoDS::Face(ex.Current());
      Bnd_Box           aBox;
      BRepBndLib::Add(aFace, aBox);
      double v[6] = {0, 0, 0, 0, 0, 0};
      if (!aBox.IsVoid())
      {
        aBox.Get(v[0], v[1], v[2], v[3], v[4], v[5]);
      }
      BRepAdaptor_Surface aS(aFace, false);
      int                 nw = 0, ne = 0;
      for (TopExp_Explorer we(aFace, TopAbs_WIRE); we.More(); we.Next())
      {
        ++nw;
        for (TopExp_Explorer ee(we.Current(), TopAbs_EDGE); ee.More(); ee.Next())
        {
          ++ne;
        }
      }
      std::cout << std::fixed << std::setprecision(9) << "WIREINV f=" << k << " type=" << (int)aS.GetType()
                << " wires=" << nw << " edges=" << ne << " bbox=(" << v[0] << "," << v[1] << "," << v[2]
                << ")-(" << v[3] << "," << v[4] << "," << v[5] << ")\n";
    }
    return 0;
  }

  // T-93 oracle: list every face of a given `GeomAbs_SurfaceType` with its bbox
  // and wire/edge counts. Used to locate a band/side face by surface identity
  // when the face order between imports differs.
  if (argc > 3 && std::string(argv[2]) == "--typefaces")
  {
    const int aType = std::atoi(argv[3]);
    const double aR  = (argc > 4) ? std::atof(argv[4]) : -1.0;
    int          k   = 0;
    for (TopExp_Explorer ex(aShape, TopAbs_FACE); ex.More(); ex.Next(), ++k)
    {
      const TopoDS_Face       aFace = TopoDS::Face(ex.Current());
      BRepAdaptor_Surface     aS(aFace, false);
      if ((int)aS.GetType() != aType)
      {
        continue;
      }
      if (aR >= 0.0)
      {
        double r = -1.0;
        if (aS.GetType() == GeomAbs_Sphere)
        {
          r = aS.Sphere().Radius();
        }
        else if (aS.GetType() == GeomAbs_Cylinder)
        {
          r = aS.Cylinder().Radius();
        }
        if (std::abs(r - aR) > 1e-6)
        {
          continue;
        }
      }
      Bnd_Box aBox;
      BRepBndLib::Add(aFace, aBox);
      double v[6] = {0, 0, 0, 0, 0, 0};
      if (!aBox.IsVoid())
      {
        aBox.Get(v[0], v[1], v[2], v[3], v[4], v[5]);
      }
      int nw = 0, ne = 0;
      for (TopExp_Explorer we(aFace, TopAbs_WIRE); we.More(); we.Next())
      {
        ++nw;
        for (TopExp_Explorer ee(we.Current(), TopAbs_EDGE); ee.More(); ee.Next())
        {
          ++ne;
        }
      }
      std::cout << std::fixed << std::setprecision(6) << "TF f=" << k << " type=" << aType
                << " wires=" << nw << " edges=" << ne << " bbox=(" << v[0] << "," << v[1] << "," << v[2]
                << ")-(" << v[3] << "," << v[4] << "," << v[5] << ")";
      std::cout << std::setprecision(17) << " urange=[" << aS.FirstUParameter() << ","
                << aS.LastUParameter() << "] vrange=[" << aS.FirstVParameter() << ","
                << aS.LastVParameter() << "] uper=" << (aS.IsUPeriodic() ? 1 : 0)
                << " uclosed=" << (aS.IsUClosed() ? 1 : 0) << " vper=" << (aS.IsVPeriodic() ? 1 : 0)
                << " vclosed=" << (aS.IsVClosed() ? 1 : 0);
      if (aS.GetType() == GeomAbs_Sphere)
      {
        const gp_Sphere aSph = aS.Sphere();
        std::cout << " sphR=" << aSph.Radius() << " sphLoc=(" << aSph.Position().Location().X() << ","
                  << aSph.Position().Location().Y() << "," << aSph.Position().Location().Z() << ")";
        std::cout << " sphAxis=(" << aSph.Position().Direction().X() << ","
                  << aSph.Position().Direction().Y() << "," << aSph.Position().Direction().Z() << ")";
      }
      std::cout << "\n";
      // Each edge with its pcurve parameter range, so a single-edge loop is
      // visible as "the whole boundary is one curve".
      int ei = 0;
      for (TopExp_Explorer we(aFace, TopAbs_WIRE); we.More(); we.Next())
      {
        for (TopExp_Explorer ee(we.Current(), TopAbs_EDGE); ee.More(); ee.Next(), ++ei)
        {
          const TopoDS_Edge anE = TopoDS::Edge(ee.Current());
          double            f = 0., l = 0., pf = 0., pl = 0.;
          const occ::handle<Geom_Curve>&   aC  = BRep_Tool::Curve(anE, f, l);
          const occ::handle<Geom2d_Curve>& aPC = BRep_Tool::CurveOnSurface(anE, aFace, pf, pl);
          std::cout << "TF   e" << ei << " ori=" << (int)anE.Orientation()
                    << " deg=" << (BRep_Tool::Degenerated(anE) ? 1 : 0) << " par=[" << f << "," << l
                    << "]";
          if (!aPC.IsNull())
          {
            const gp_Pnt2d aF = aPC->Value(pf);
            const gp_Pnt2d aL = aPC->Value(pl);
            std::cout << " pcpar=[" << pf << "," << pl << "] pc_first=(" << aF.X() << "," << aF.Y()
                      << ") pc_last=(" << aL.X() << "," << aL.Y() << ")";
          }
          std::cout << "\n";
        }
      }
    }
    return 0;
  }

  // T-69/T-93 oracle: the *imported* shape's own boundary for one face, read
  // straight off `TopoDS_Face` without the discrete model. Answers "how many
  // edges does OCCT's import put in this face loop, and what are they" so the
  // port's wire can be diffed edge by edge.
  //
  // Selection is by INDEX (`--adv <k>`) or by BBOX (`--advbox x0 y0 z0 x1 y1 z1`,
  // tolerance 0.5) because the face order of `top_exp_face` and the order the
  // port enumerates faces are independent -- pairing must be geometric.
  {
    // Mode may appear anywhere from argv[2] on (so `--nofix` etc. can be
    // appended without shifting the geometry arguments).
    int    aWant = -1;
    bool   byBox = false;
    double q[6] = {0, 0, 0, 0, 0, 0};
    for (int i = 2; i < argc; ++i)
    {
      const std::string a = argv[i];
      if (a == "--adv" && i + 1 < argc)
      {
        aWant = std::atoi(argv[i + 1]);
      }
      else if (a == "--advbox" && i + 6 < argc)
      {
        byBox = true;
        for (int j = 0; j < 6; ++j)
        {
          q[j] = std::atof(argv[i + 1 + j]);
        }
      }
    }
    if (aWant >= 0 || byBox)
    {
      int k = 0;
      for (TopExp_Explorer ex(aShape, TopAbs_FACE); ex.More(); ex.Next(), ++k)
      {
        const TopoDS_Face aFace = TopoDS::Face(ex.Current());
        Bnd_Box           aBox;
        BRepBndLib::Add(aFace, aBox);
        double v[6] = {0, 0, 0, 0, 0, 0};
        if (aBox.IsVoid())
        {
          continue;
        }
        aBox.Get(v[0], v[1], v[2], v[3], v[4], v[5]);
        if (byBox)
        {
          bool ok = true;
          for (int i = 0; i < 6; ++i)
          {
            if (std::abs(v[i] - q[i]) > 0.5)
            {
              ok = false;
            }
          }
          if (!ok)
          {
            continue;
          }
        }
        else if (k != aWant)
        {
          continue;
        }
        BRepAdaptor_Surface aS(aFace, false);
        std::cout << std::fixed << std::setprecision(6) << "ADV face=" << k
                  << " type=" << (int)aS.GetType() << " wires=";
        int nw = 0;
        for (TopExp_Explorer we(aFace, TopAbs_WIRE); we.More(); we.Next())
        {
          ++nw;
        }
        std::cout << nw << " bbox=(" << v[0] << "," << v[1] << "," << v[2] << ")-(" << v[3] << ","
                  << v[4] << "," << v[5] << ")";
        std::cout << std::setprecision(17) << " u=[" << aS.FirstUParameter() << ","
                  << aS.LastUParameter() << "] v=[" << aS.FirstVParameter() << ","
                  << aS.LastVParameter() << "] uper=" << (aS.IsUPeriodic() ? 1 : 0)
                  << " vper=" << (aS.IsVPeriodic() ? 1 : 0);
        if (aS.GetType() == GeomAbs_Sphere)
        {
          const gp_Sphere aSph = aS.Sphere();
          std::cout << " sphR=" << aSph.Radius() << " sphLoc=(" << aSph.Position().Location().X()
                    << "," << aSph.Position().Location().Y() << ","
                    << aSph.Position().Location().Z() << ")";
        }
        std::cout << "\n";
        int wi = 0;
        for (TopExp_Explorer we(aFace, TopAbs_WIRE); we.More(); we.Next(), ++wi)
        {
          int ne = 0;
          for (TopExp_Explorer ee(we.Current(), TopAbs_EDGE); ee.More(); ee.Next())
          {
            ++ne;
          }
          std::cout << "ADV  wire[" << wi << "] nEdges=" << ne << "\n";
          int ei = 0;
          for (TopExp_Explorer ee(we.Current(), TopAbs_EDGE); ee.More(); ee.Next(), ++ei)
          {
            const TopoDS_Edge anE = TopoDS::Edge(ee.Current());
            double            f = 0., l = 0.;
            const occ::handle<Geom_Curve>& aC = BRep_Tool::Curve(anE, f, l);
            double            pf = 0., pl = 0.;
            const occ::handle<Geom2d_Curve>& aPC = BRep_Tool::CurveOnSurface(anE, aFace, pf, pl);
            TopoDS_Vertex v1, v2;
            TopExp::Vertices(anE, v1, v2);
            const gp_Pnt p1 = BRep_Tool::Pnt(v1);
            const gp_Pnt p2 = BRep_Tool::Pnt(v2);
            std::cout << "ADV   e[" << ei << "] ori=" << (int)anE.Orientation()
                      << " deg=" << (BRep_Tool::Degenerated(anE) ? 1 : 0)
                      << " null3d=" << (aC.IsNull() ? 1 : 0) << " nullpc=" << (aPC.IsNull() ? 1 : 0);
            std::cout << std::setprecision(9) << " p1=(" << p1.X() << "," << p1.Y() << "," << p1.Z()
                      << ") p2=(" << p2.X() << "," << p2.Y() << "," << p2.Z() << ")";
            std::cout << std::setprecision(15) << " par=[" << f << "," << l << "]";
            if (!aPC.IsNull())
            {
              const gp_Pnt2d aF = aPC->Value(pf);
              const gp_Pnt2d aL = aPC->Value(pl);
              std::cout << " pcpar=[" << pf << "," << pl << "] pc_first=(" << aF.X() << "," << aF.Y()
                        << ") pc_last=(" << aL.X() << "," << aL.Y() << ")";
            }
            bool closed = false;
            if (!v1.IsNull() && !v2.IsNull())
            {
              closed = v1.IsSame(v2);
            }
            std::cout << " closedV=" << (closed ? 1 : 0) << "\n";
          }
        }
        if (!byBox)
        {
          return 0;
        }
      }
      std::cout << "ADV none\n";
      return 0;
    }
  }

  // TEMP T-104: `--structdump <faceIdx>` (dispatched here so the imported shape
  // is available; the index is the `IMeshData_Model` face order).
  if (aStructDumpFace >= 0)
  {
    return runStructDump(aShape, aStructDumpFace, false, false, aStructDumpDefl, aStructDumpAngle);
  }
  // TEMP T-119: factory-faithful splitter + DC pass dump.
  if (aDCStructFace >= 0)
  {
    return runStructDump(aShape, aDCStructFace, true, false, aStructDumpDefl, aStructDumpAngle);
  }
  // TEMP T-120: factory-faithful splitter + surface-node classifier verdicts.
  if (aSrfNodeFace >= 0)
  {
    return runStructDump(aShape, aSrfNodeFace, false, true);
  }

  // T-69 oracle: the (u,v) point set that BRepMesh_NodeInsertionMeshAlgo::
  // collectWirePoints feeds into BRepMesh_DefaultRangeSplitter::AddPoint for one
  // face, plus the min/max range those points produce. This is the exact input
  // whose spread decides `computeLengthU/V` and therefore `IsValid`
  // (`BRepMesh_DefaultRangeSplitter.cxx:35-41` / `:72-80`).
  if (argc > 2 && (std::string(argv[2]) == "--uv" || std::string(argv[2]) == "--uvsum"
                   || std::string(argv[2]) == "--nofix"))
  {
    int  aWant = -1;
    bool doSum = false;
    for (int i = 2; i < argc; ++i)
    {
      const std::string a = argv[i];
      if (a == "--uv" && i + 1 < argc)
      {
        aWant = std::atoi(argv[i + 1]);
      }
      else if (a == "--uvsum")
      {
        doSum = true;
      }
    }
    IMeshTools_Parameters aParams;
    aParams.Deflection    = 0.1;
    aParams.Angle         = 0.349066;
    aParams.InParallel    = false;
    aParams.Relative      = false;
    aParams.MinSize       = Precision::Confusion();
    aParams.AdjustMinSize = false;

    occ::handle<BRepMesh_ModelBuilder> aBuilder = new BRepMesh_ModelBuilder;
    const occ::handle<IMeshData_Model> aModel   = aBuilder->Perform(aShape, aParams);
    if (aModel.IsNull())
    {
      std::cout << "UV model-null\n";
      return 0;
    }
    occ::handle<BRepMesh_EdgeDiscret> anEdgeDiscret = new BRepMesh_EdgeDiscret;
    const bool aDiscretDone = anEdgeDiscret->Perform(aModel, aParams, Message_ProgressRange());
    // T-99: this branch only builds the discrete model and discretizes edges - it
    // never triangulates, so every face's `Poly_Triangulation` is null and the
    // summary's `triangles=` field came out 0 (measured). Triangulate the shape
    // here so that field reports OCCT's real per-face triangle count. The
    // deflection is taken from argv ("--uvsum <d>", default 0.1) so the caller can
    // match whatever deflection the ground-truth OBJ was produced with.
    {
      double aUvDefl = 0.1;
      for (int i = 2; i + 1 < argc; ++i)
      {
        if (std::string(argv[i]) == "--uvsum")
        {
          aUvDefl = std::atof(argv[i + 1]);
        }
      }
      BRepMesh_IncrementalMesh aSumMesher(aShape, aUvDefl, false, 0.5);
    }
    {
      int aPCurvesTotal = 0, aFreeEdges = 0, aFilled = 0, aEmptyPCurves = 0;
      for (int e = 0; e < aModel->EdgesNb(); ++e)
      {
        const IMeshData::IEdgeHandle& aE = aModel->GetEdge(e);
        if (aE->IsFree())
        {
          ++aFreeEdges;
        }
        for (int p = 0; p < aE->PCurvesNb(); ++p)
        {
          ++aPCurvesTotal;
          const int aNb = aE->GetPCurve(p)->ParametersNb();
          if (aNb > 0)
          {
            ++aFilled;
          }
          else
          {
            ++aEmptyPCurves;
          }
        }
      }
      std::cout << "UV model edges=" << aModel->EdgesNb() << " faces=" << aModel->FacesNb()
                << " pcurves=" << aPCurvesTotal << " free=" << aFreeEdges << " filled=" << aFilled
                << " empty=" << aEmptyPCurves << " discret_done=" << (aDiscretDone ? 1 : 0)
                << " maxsize=" << aModel->GetMaxSize() << "\n";
    }

    if (aWant >= 0)
    {
    const IMeshData::IFaceHandle& aDFace = aModel->GetFace(aWant);
    const TopoDS_Face&            aFace  = aDFace->GetFace();
    {
      Bnd_Box aBox;
      BRepBndLib::Add(aFace, aBox);
      double v[6] = {0, 0, 0, 0, 0, 0};
      if (!aBox.IsVoid())
      {
        aBox.Get(v[0], v[1], v[2], v[3], v[4], v[5]);
      }
      std::cout << std::fixed << std::setprecision(6) << "UV face=" << aWant << " wires=" << aDFace->WiresNb()
                << " u=[" << aDFace->GetSurface()->FirstUParameter() << ","
                << aDFace->GetSurface()->LastUParameter() << "] v=["
                << aDFace->GetSurface()->FirstVParameter() << ","
                << aDFace->GetSurface()->LastVParameter() << "]"
                << " uper=" << (aDFace->GetSurface()->IsUPeriodic() ? 1 : 0)
                << " vper=" << (aDFace->GetSurface()->IsVPeriodic() ? 1 : 0)
                << " bbox=(" << v[0] << "," << v[1] << "," << v[2] << ")-(" << v[3] << "," << v[4] << ","
                << v[5] << ")\n";
    }

    double uMin = 1e100, uMax = -1e100, vMin = 1e100, vMax = -1e100;
    int    nPts = 0, onUMin = 0, onUMax = 0, onVMin = 0, onVMax = 0;
    for (int aWireIt = 0; aWireIt < aDFace->WiresNb(); ++aWireIt)
    {
      const IMeshData::IWireHandle& aDWire = aDFace->GetWire(aWireIt);
      for (int aEdgeIt = 0; aEdgeIt < aDWire->EdgesNb(); ++aEdgeIt)
      {
        const IMeshData::IEdgeHandle& aDEdge = aDWire->GetEdge(aEdgeIt);
        const TopAbs_Orientation      anOri = aDWire->GetEdgeOrientation(aEdgeIt);
        const IMeshData::IPCurveHandle& aPCurve = aDEdge->GetPCurve(aDFace.get(), anOri);
        std::cout << "UV  w" << aWireIt << " e" << aEdgeIt << " ori=" << (int)anOri
                  << " prefix=" << (anOri == TopAbs_REVERSED ? "REV" : "FWD")
                  << " nat=" << (aPCurve->IsForward() ? 1 : 0)
                  << " npc=" << aPCurve->ParametersNb();
        if (aPCurve->ParametersNb() > 0)
        {
          const gp_Pnt2d aF = aPCurve->GetPoint(0);
          const gp_Pnt2d aL = aPCurve->GetPoint(aPCurve->ParametersNb() - 1);
          std::cout << " first=(" << aF.X() << "," << aF.Y() << ") last=(" << aL.X() << "," << aL.Y() << ")"
                    << " p0=" << aPCurve->GetParameter(0) << " pn="
                    << aPCurve->GetParameter(aPCurve->ParametersNb() - 1);
        }
        std::cout << "\n";
        for (int aPt = 0; aPt < aPCurve->ParametersNb(); ++aPt)
        {
          const gp_Pnt2d& aP = aPCurve->GetPoint(aPt);
          std::cout << std::setprecision(17) << "UV   p" << aPt << " uv=(" << aP.X() << "," << aP.Y()
                    << ") t=" << aPCurve->GetParameter(aPt) << "\n";
          ++nPts;
          uMin = std::min(uMin, aP.X());
          uMax = std::max(uMax, aP.X());
          vMin = std::min(vMin, aP.Y());
          vMax = std::max(vMax, aP.Y());
        }
      }
    }
    // Count how many collected points sit exactly on the extremal coordinates:
    // a zero-width direction shows every point on both the min and the max.
    for (int aWireIt = 0; aWireIt < aDFace->WiresNb(); ++aWireIt)
    {
      const IMeshData::IWireHandle& aDWire = aDFace->GetWire(aWireIt);
      for (int aEdgeIt = 0; aEdgeIt < aDWire->EdgesNb(); ++aEdgeIt)
      {
        const IMeshData::IPCurveHandle& aPCurve =
          aDWire->GetEdge(aEdgeIt)->GetPCurve(aDFace.get(), aDWire->GetEdgeOrientation(aEdgeIt));
        for (int aPt = 0; aPt < aPCurve->ParametersNb(); ++aPt)
        {
          const gp_Pnt2d& aP = aPCurve->GetPoint(aPt);
          if (aP.X() == uMin) ++onUMin;
          if (aP.X() == uMax) ++onUMax;
          if (aP.Y() == vMin) ++onVMin;
          if (aP.Y() == vMax) ++onVMax;
        }
      }
    }
    // Reproduce `BRepMesh_NodeInsertionMeshAlgo::collectWirePoints`
    // (`BRepMesh_NodeInsertionMeshAlgo.hxx:142-184`): every pcurve sample except
    // the traversal's last one.
    double rUMin = 1e100, rUMax = -1e100, rVMin = 1e100, rVMax = -1e100;
    for (int aWireIt = 0; aWireIt < aDFace->WiresNb(); ++aWireIt)
    {
      const IMeshData::IWireHandle& aDWire = aDFace->GetWire(aWireIt);
      if (aDWire->IsSet(IMeshData_SelfIntersectingWire)
          || (aDWire->IsSet(IMeshData_OpenWire) && aWireIt != 0))
      {
        continue;
      }
      for (int aEdgeIt = 0; aEdgeIt < aDWire->EdgesNb(); ++aEdgeIt)
      {
        const IMeshData::IPCurveHandle& aPCurve =
          aDWire->GetEdge(aEdgeIt)->GetPCurve(aDFace.get(), aDWire->GetEdgeOrientation(aEdgeIt));
        int aIt, aEnd, aInc;
        if (aPCurve->IsForward())
        {
          aEnd = aPCurve->ParametersNb() - 1;
          aIt  = (std::min)(0, aEnd);
          aInc = 1;
        }
        else
        {
          aIt  = aPCurve->ParametersNb() - 1;
          aEnd = (std::min)(0, aIt);
          aInc = -1;
        }
        for (; aIt != aEnd; aIt += aInc)
        {
          const gp_Pnt2d& aP = aPCurve->GetPoint(aIt);
          rUMin = std::min(rUMin, aP.X());
          rUMax = std::max(rUMax, aP.X());
          rVMin = std::min(rVMin, aP.Y());
          rVMax = std::max(rVMax, aP.Y());
        }
      }
    }
    std::cout << std::setprecision(17) << "UV RANGE pts=" << nPts << " u=[" << uMin << "," << uMax
              << "] v=[" << vMin << "," << vMax << "] dU=" << (uMax - uMin) << " dV=" << (vMax - vMin)
              << " onUMin=" << onUMin << " onUMax=" << onUMax << " onVMin=" << onVMin
              << " onVMax=" << onVMax << "\n";
    std::cout << std::setprecision(17) << "UV ADDEDRANGE u=[" << rUMin << "," << rUMax
              << "] v=[" << rVMin << "," << rVMax << "] dU=" << (rUMax - rUMin)
              << " dV=" << (rVMax - rVMin) << "\n";
    // The predicate `AdjustRange` uses to decide `IsValid`
    // (`BRepMesh_DefaultRangeSplitter.cxx:72-80`) on the added-point range.
    const double aDu = 0.05 * (rUMax - rUMin);
    const double aDv = 0.05 * (rVMax - rVMin);
    const occ::handle<BRepAdaptor_Surface>& aSurf = aDFace->GetSurface();
    double      aLenU = 0., aLenV = 0.;
    gp_Pnt      P11, P12, P21, P22, P31, P32;
    const double aVave = 0.5 * (rVMax + rVMin);
    aSurf->D0(rUMin, rVMin, P11);
    aSurf->D0(rUMin, aVave, P21);
    aSurf->D0(rUMin, rVMax, P31);
    for (int i1 = 1; i1 <= 20; ++i1)
    {
      const double aU = rUMin + aDu * i1;
      aSurf->D0(aU, rVMin, P12);
      aSurf->D0(aU, aVave, P22);
      aSurf->D0(aU, rVMax, P32);
      aLenU += (P11.Distance(P12) + P21.Distance(P22) + P31.Distance(P32));
      P11 = P12;
      P21 = P22;
      P31 = P32;
    }
    const double aUave = 0.5 * (rUMax + rUMin);
    aSurf->D0(rUMin, rVMin, P11);
    aSurf->D0(aUave, rVMin, P21);
    aSurf->D0(rUMax, rVMin, P31);
    for (int i1 = 1; i1 <= 20; ++i1)
    {
      const double aV = rVMin + aDv * i1;
      aSurf->D0(rUMin, aV, P12);
      aSurf->D0(aUave, aV, P22);
      aSurf->D0(rUMax, aV, P32);
      aLenV += (P11.Distance(P12) + P21.Distance(P22) + P31.Distance(P32));
      P11 = P12;
      P21 = P22;
      P31 = P32;
    }
    aLenU /= 3.;
    aLenV /= 3.;
    std::cout << std::setprecision(17) << "UV LEN lenU=" << aLenU << " lenV=" << aLenV
              << " pconf=" << Precision::PConfusion() << " valid="
              << ((aLenU > Precision::PConfusion() && aLenV > Precision::PConfusion()) ? 1 : 0) << "\n";
    if (!doSum)
    {
      return 0;
    }
    }  // end single-face (`--uv <k>`) block

    // Summary: one line per face, so the whole shape can be paired by geometry.
    for (int f = 0; f < aModel->FacesNb(); ++f)
    {
      const IMeshData::IFaceHandle& aDF = aModel->GetFace(f);
      Bnd_Box                       aBox;
      BRepBndLib::Add(aDF->GetFace(), aBox);
      double w[6] = {0, 0, 0, 0, 0, 0};
      if (!aBox.IsVoid())
      {
        aBox.Get(w[0], w[1], w[2], w[3], w[4], w[5]);
      }
      // T-99: this face's own triangulation, so a per-face triangle count can be
      // attached to the SAME face index as this dump. `--mesh` cannot be used for
      // that: it enumerates `TopExp_Explorer(aShape, TopAbs_FACE)` while this loop
      // enumerates `IMeshData_Model::GetFace(f)`, and the two orderings were
      // measured to agree at only 1 of 226 indices. Taking the count HERE removes
      // every cross-dump index assumption. `BRep_Tool::Triangulation` is the
      // OCCT API for exactly this (see `--mesh` above); a face may legitimately
      // carry no triangulation, which prints as 0.
      int aFaceTris = 0;
      {
        TopLoc_Location                        aLoc;
        const occ::handle<Poly_Triangulation>& aTri =
          BRep_Tool::Triangulation(TopoDS::Face(aDF->GetFace()), aLoc);
        if (!aTri.IsNull())
        {
          aFaceTris = aTri->NbTriangles();
        }
      }
      double sUMin = 1e100, sUMax = -1e100, sVMin = 1e100, sVMax = -1e100;
      for (int aWireIt = 0; aWireIt < aDF->WiresNb(); ++aWireIt)
      {
        const IMeshData::IWireHandle& aDWire = aDF->GetWire(aWireIt);
        if (aDWire->IsSet(IMeshData_SelfIntersectingWire)
            || (aDWire->IsSet(IMeshData_OpenWire) && aWireIt != 0))
        {
          continue;
        }
        for (int aEdgeIt = 0; aEdgeIt < aDWire->EdgesNb(); ++aEdgeIt)
        {
          const IMeshData::IPCurveHandle& aPC =
            aDWire->GetEdge(aEdgeIt)->GetPCurve(aDF.get(), aDWire->GetEdgeOrientation(aEdgeIt));
          int aIt, aEnd, aInc;
          if (aPC->IsForward())
          {
            aEnd = aPC->ParametersNb() - 1;
            aIt  = (std::min)(0, aEnd);
            aInc = 1;
          }
          else
          {
            aIt  = aPC->ParametersNb() - 1;
            aEnd = (std::min)(0, aIt);
            aInc = -1;
          }
          for (; aIt != aEnd; aIt += aInc)
          {
            const gp_Pnt2d& aP = aPC->GetPoint(aIt);
            sUMin = std::min(sUMin, aP.X());
            sUMax = std::max(sUMax, aP.X());
            sVMin = std::min(sVMin, aP.Y());
            sVMax = std::max(sVMax, aP.Y());
          }
        }
      }
      const int    sType = (int)aDF->GetSurface()->GetType();
      const double sDU = sUMax - sUMin, sDV = sVMax - sVMin;
      // T-99: reproduce the range splitter's own intermediate values at the same
      // pipeline point the port reports (`adjust_range` -> `generate_surface_nodes`),
      // so the two sides can be diffed value by value. Every formula on this path was
      // verified line-for-line against the .cxx, so a divergence has to come from
      // these INPUTS. The splitter is built here exactly as the port does it:
      // Reset(face) + AddPoint(each pcurve sample) + AdjustRange().
      if (std::getenv("OCCT_TOPO_TRACE_SPLITTER") != nullptr)
      {
        BRepMesh_DefaultRangeSplitter aSplitter;
        aSplitter.Reset(aDF, aParams);
        for (int aWi = 0; aWi < aDF->WiresNb(); ++aWi)
        {
          const IMeshData::IWireHandle& aDW = aDF->GetWire(aWi);
          for (int aEi = 0; aEi < aDW->EdgesNb(); ++aEi)
          {
            const IMeshData::IPCurveHandle& aPC =
              aDW->GetEdge(aEi)->GetPCurve(aDF.get(), aDW->GetEdgeOrientation(aEi));
            for (int aPt = 0; aPt < aPC->ParametersNb(); ++aPt)
            {
              aSplitter.AddPoint(aPC->GetPoint(aPt));
            }
          }
        }
        aSplitter.AdjustRange();
        const std::pair<double, double>& aRu = aSplitter.GetRangeU();
        const std::pair<double, double>& aRv = aSplitter.GetRangeV();
        const std::pair<double, double>& aDe = aSplitter.GetDelta();
        const std::pair<double, double>& aTo = aSplitter.GetToleranceUV();
        std::cout << "SPLITTER f=" << f << " preset=1 ru=[" << std::setprecision(15)
                  << aRu.first << "," << aRu.second << "] rv=[" << aRv.first << "," << aRv.second
                  << "] delta=[" << aDe.first << "," << aDe.second << "]";
        std::cout << std::scientific << std::setprecision(6) << " tol=[" << aTo.first << ","
                  << aTo.second << "]";
        std::cout << std::fixed << " vnb=0 defl=" << aParams.Deflection;
        // T-99: the cells count this face would get. This is the decision number -
        // range, delta and every formula on the path already matched between the port
        // and OCCT (see specs/_a3n00_gap_analysis.md §9.343/§9.344/§9.346), so if the
        // cell counts agree the divergence must be downstream of it. Printed twice
        // because OCCT passes `aStructure->NbNodes()` (Deleted included) while the
        // port's `node_insertion.rs:362` excludes Deleted nodes.
        const std::pair<int, int> aCellsZero =
          BRepMesh_GeomTool::CellsCount(aDF->GetSurface(), 0, aParams.Deflection, &aSplitter);
        std::cout << " cells_vnb0=[" << aCellsZero.first << "," << aCellsZero.second << "]";
        std::cout << std::setprecision(6) << " bbox=(" << w[0] << "," << w[1] << "," << w[2]
                  << ")-(" << w[3] << "," << w[4] << "," << w[5] << ")\n";
      }
      const double sDu = 0.05 * sDU, sDv = 0.05 * sDV;
      const occ::handle<BRepAdaptor_Surface>& aS2 = aDF->GetSurface();
      double lU = 0., lV = 0.;
      gp_Pnt A11, A12, A21, A22, A31, A32;
      const double vA = 0.5 * (sVMax + sVMin);
      aS2->D0(sUMin, sVMin, A11);
      aS2->D0(sUMin, vA, A21);
      aS2->D0(sUMin, sVMax, A31);
      for (int i1 = 1; i1 <= 20; ++i1)
      {
        const double aU = sUMin + sDu * i1;
        aS2->D0(aU, sVMin, A12);
        aS2->D0(aU, vA, A22);
        aS2->D0(aU, sVMax, A32);
        lU += (A11.Distance(A12) + A21.Distance(A22) + A31.Distance(A32));
        A11 = A12;
        A21 = A22;
        A31 = A32;
      }
      const double uA = 0.5 * (sUMax + sUMin);
      aS2->D0(sUMin, sVMin, A11);
      aS2->D0(uA, sVMin, A21);
      aS2->D0(sUMax, sVMin, A31);
      for (int i1 = 1; i1 <= 20; ++i1)
      {
        const double aV = sVMin + sDv * i1;
        aS2->D0(sUMin, aV, A12);
        aS2->D0(uA, aV, A22);
        aS2->D0(sUMax, aV, A32);
        lV += (A11.Distance(A12) + A21.Distance(A22) + A31.Distance(A32));
        A11 = A12;
        A21 = A22;
        A31 = A32;
      }
      lU /= 3.;
      lV /= 3.;
      const bool sValid = (lU > Precision::PConfusion() && lV > Precision::PConfusion());
      std::cout << std::fixed << std::setprecision(9) << "UVSUM f=" << f << " type=" << sType
                << " wires=" << aDF->WiresNb() << " triangles=" << aFaceTris
                << " bbox=(" << w[0] << "," << w[1] << "," << w[2]
                << ")-(" << w[3] << "," << w[4] << "," << w[5] << ")";
      std::cout << std::setprecision(17) << " urange=[" << sUMin << "," << sUMax << "] vrange=["
                << sVMin << "," << sVMax << "] dU=" << sDU << " dV=" << sDV << " lenU=" << lU
                << " lenV=" << lV << " valid=" << (sValid ? 1 : 0) << "\n";
    }
    return 0;
  }

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
                    << ") deg=" << BRep_Tool::Degenerated(E) << "\n";
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
  // TEMP T-121: `--nbdump` prints, per face, the values `ShapeFix_Face::
  // isNeedAddNaturalBound` (`ShapeFix_Face.cxx:1120-1160`) reads: the surface
  // U/V bounds and periodicity, the wire count, `ShapeAnalysis::IsOuterBound`
  // (`ShapeAnalysis.cxx:203-230`) and, per wire, `ShapeAnalysis::TotCross2D`
  // (`ShapeAnalysis.cxx:115-150`) with every edge orientation and p-curve
  // range, so the port's decision can be compared step by step.
  if (argc > 2 && std::string(argv[2]) == "--nbdump")
  {
    int k = 0;
    for (TopExp_Explorer fx(aShape, TopAbs_FACE); fx.More(); fx.Next(), ++k)
    {
      const TopoDS_Face aF = TopoDS::Face(fx.Current());
      BRepAdaptor_Surface aS(aF);
      int aNW = 0;
      for (TopExp_Explorer wx(aF, TopAbs_WIRE); wx.More(); wx.Next())
      {
        ++aNW;
      }
      std::cout << "NB face " << k << " type=" << (int)aS.GetType()
                << " uper=" << (aS.IsUPeriodic() ? 1 : 0) << " vper=" << (aS.IsVPeriodic() ? 1 : 0)
                << " u=[" << aS.FirstUParameter() << "," << aS.LastUParameter() << "]"
                << " v=[" << aS.FirstVParameter() << "," << aS.LastVParameter() << "]"
                << " wires=" << aNW << " outer=" << (ShapeAnalysis::IsOuterBound(aF) ? 1 : 0) << "\n";
      int aWI = 0;
      for (TopExp_Explorer wx(aF, TopAbs_WIRE); wx.More(); wx.Next(), ++aWI)
      {
        TopoDS_Wire                       aW = TopoDS::Wire(wx.Current());
        occ::handle<ShapeExtend_WireData> aWD = new ShapeExtend_WireData(aW);
        std::cout << "NB  w" << aWI << " ori=" << (int)aW.Orientation()
                  << " n=" << aWD->NbEdges() << " cross=" << ShapeAnalysis::TotCross2D(aWD, aF)
                  << "\n";
        for (int i = 1; i <= aWD->NbEdges(); ++i)
        {
          const TopoDS_Edge&        e = aWD->Edge(i);
          double                    f = 0., l = 0.;
          occ::handle<Geom2d_Curve> c2 = BRep_Tool::CurveOnSurface(e, aF, f, l);
          gp_Pnt2d                  b(0., 0.), en(0., 0.);
          if (!c2.IsNull())
          {
            b  = c2->Value(f);
            en = c2->Value(l);
          }
          std::cout << "NB   e" << i << " ori=" << (int)e.Orientation()
                    << " closed=" << (BRep_Tool::IsClosed(e, aF) ? 1 : 0)
                    << " deg=" << (BRep_Tool::Degenerated(e) ? 1 : 0) << " uv=(" << b.X() << ","
                    << b.Y() << ")->(" << en.X() << "," << en.Y() << ")\n";
        }
      }
    }
    return 0;
  }
  // TEMP T-121: `--vmap` walks every wire of every face and prints, per edge,
  // the identity and point of `ShapeAnalysis_Edge::FirstVertex` / `LastVertex`
  // (`ShapeAnalysis_Edge.cxx:228-257`), so edge closure (`v1.IsSame(v2)`) and
  // vertex sharing between adjacent edges can be read off. The identity is the
  // raw `TShape` address, stable inside one process.
  if (argc > 2 && std::string(argv[2]) == "--vmap")
  {
    int k = 0;
    for (TopExp_Explorer fx(aShape, TopAbs_FACE); fx.More(); fx.Next(), ++k)
    {
      std::cout << "VMAP face " << k << "\n";
      int aWI = 0;
      for (TopExp_Explorer wx(fx.Current(), TopAbs_WIRE); wx.More(); wx.Next(), ++aWI)
      {
        ShapeAnalysis_Edge                aEA;
        ShapeExtend_WireData              aWD(TopoDS::Wire(wx.Current()));
        std::cout << "VMAP  w" << aWI << " n=" << aWD.NbEdges() << "\n";
        for (int i = 1; i <= aWD.NbEdges(); ++i)
        {
          const TopoDS_Edge   e  = aWD.Edge(i);
          const TopoDS_Vertex v1 = aEA.FirstVertex(e);
          const TopoDS_Vertex v2 = aEA.LastVertex(e);
          std::cout << "VMAP   e" << i << " ori=" << (int)e.Orientation()
                    << " v1=" << (v1.IsNull() ? 0L : (long)v1.TShape().get())
                    << " v2=" << (v2.IsNull() ? 0L : (long)v2.TShape().get())
                    << " same=" << ((!v1.IsNull() && !v2.IsNull() && v1.IsSame(v2)) ? 1 : 0);
          if (!v1.IsNull())
          {
            const gp_Pnt p = BRep_Tool::Pnt(v1);
            std::cout << " p1=(" << p.X() << "," << p.Y() << "," << p.Z() << ")";
          }
          if (!v2.IsNull())
          {
            const gp_Pnt p = BRep_Tool::Pnt(v2);
            std::cout << " p2=(" << p.X() << "," << p.Y() << "," << p.Z() << ")";
          }
          std::cout << "\n";
        }
      }
    }
    return 0;
  }
  // TEMP T-121: `--wirefix <faceIdx>` replays `ShapeFix_Face::Perform`'s wire
  // pass (`ShapeFix_Face.cxx:400-455`) on one face of the loaded (usually
  // `--nofix`) shape: for the full pass and then for every pass with one
  // `ShapeFix_Wire` mode disabled, it prints `ShapeAnalysis::TotCross2D`
  // (`ShapeAnalysis.cxx:115-150`) and the resulting edge orientations, so the
  // port's `perform_wire_first_round` can be compared stage by stage.
  if (argc > 3 && std::string(argv[2]) == "--wirefix")
  {
    const int aFaceIdx = std::atoi(argv[3]);
    TopoDS_Face aFace;
    int         k = 0;
    for (TopExp_Explorer fx(aShape, TopAbs_FACE); fx.More(); fx.Next(), ++k)
    {
      if (k == aFaceIdx)
      {
        aFace = TopoDS::Face(fx.Current());
        break;
      }
    }
    if (aFace.IsNull())
    {
      std::cout << "WF no-face\n";
      return 0;
    }
    occ::handle<Geom_Surface>          aSurf = BRep_Tool::Surface(aFace);
    occ::handle<ShapeAnalysis_Surface> aSAS  = new ShapeAnalysis_Surface(aSurf);
    const TopoDS_Wire aW = TopoDS::Wire(TopExp_Explorer(aFace, TopAbs_WIRE).Current());

    auto report = [&](const char* tag, const TopoDS_Wire& w) {
      occ::handle<ShapeExtend_WireData> aWD = new ShapeExtend_WireData(w);
      std::cout << "WF " << tag << " n=" << aWD->NbEdges()
                << " cross=" << ShapeAnalysis::TotCross2D(aWD, aFace) << "\n";
      for (int i = 1; i <= aWD->NbEdges(); ++i)
      {
        const TopoDS_Edge&        e = aWD->Edge(i);
        double                    f = 0., l = 0.;
        occ::handle<Geom2d_Curve> c2 = BRep_Tool::CurveOnSurface(e, aFace, f, l);
        gp_Pnt2d                  b(0., 0.), en(0., 0.);
        if (!c2.IsNull())
        {
          b  = c2->Value(f);
          en = c2->Value(l);
        }
        std::cout << "WF   e" << i << " ori=" << (int)e.Orientation() << " uv=(" << b.X() << ","
                  << b.Y() << ")->(" << en.X() << "," << en.Y() << ")\n";
      }
    };

    report("raw", aW);

    {
      ShapeAnalysis_Edge aEA;
      occ::handle<ShapeExtend_WireData> aWD = new ShapeExtend_WireData(aW);
      for (int i = 1; i <= aWD->NbEdges(); ++i)
      {
        const TopoDS_Edge   e = aWD->Edge(i);
        const TopoDS_Vertex v1 = aEA.FirstVertex(e);
        const TopoDS_Vertex v2 = aEA.LastVertex(e);
        int                 aNbV = 0;
        for (TopoDS_Iterator aVi(e, false); aVi.More(); aVi.Next())
        {
          ++aNbV;
        }
        std::cout << "WF v" << i << " children=" << aNbV
                  << " same=" << ((!v1.IsNull() && !v2.IsNull() && v1.IsSame(v2)) ? 1 : 0)
                  << " null=" << (v1.IsNull() || v2.IsNull() ? 1 : 0);
        if (!v1.IsNull() && !v2.IsNull())
        {
          const gp_Pnt p1 = BRep_Tool::Pnt(v1);
          const gp_Pnt p2 = BRep_Tool::Pnt(v2);
          std::cout << " p1=(" << p1.X() << "," << p1.Y() << "," << p1.Z() << ") p2=(" << p2.X()
                    << "," << p2.Y() << "," << p2.Z() << ")";
        }
        std::cout << "\n";
      }
    }

    {
      ShapeAnalysis_Wire      aWT(new ShapeExtend_WireData(aW), aFace, Precision::Confusion());
      ShapeAnalysis_WireOrder aO3;
      aWT.CheckOrder(aO3, true, true);
      std::cout << "WF order3d nb=" << aO3.NbEdges() << " status=" << aO3.Status() << " ordered=[";
      for (int i = 1; i <= aO3.NbEdges(); ++i)
      {
        std::cout << (i > 1 ? " " : "") << aO3.Ordered(i);
      }
      std::cout << "]\n";
      ShapeAnalysis_WireOrder aO2;
      aWT.CheckOrder(aO2, true, false);
      std::cout << "WF order2d nb=" << aO2.NbEdges() << " status=" << aO2.Status() << " ordered=[";
      for (int i = 1; i <= aO2.NbEdges(); ++i)
      {
        std::cout << (i > 1 ? " " : "") << aO2.Ordered(i);
      }
      std::cout << "]\n";
    }

    const char* aModes[] = {"ALL",          "FixReorder",  "FixSmall",     "FixConnected",
                            "FixEdgeCurves", "FixDegenerated", "FixTails",
                            "FixShifted",   "FixAddPCurve", "FixAddCurve3d", "FixReversed2d",
                            "FixSeam",      "ROff+COn",     "ROff+SmallOn"};
    const int   aNbModes = (int)(sizeof(aModes) / sizeof(aModes[0]));
    for (int m = 0; m < aNbModes; ++m)
    {
      ShapeFix_Wire aSFW;
      aSFW.SetFace(aFace, aSAS);
      aSFW.Load(aW);
      aSFW.FixLackingMode()          = false;
      aSFW.FixSelfIntersectionMode() = false;
      const std::string mn = aModes[m];
      if (mn == "FixReorder")
        aSFW.FixReorderMode() = 0;
      else if (mn == "FixSmall")
        aSFW.FixSmallMode() = 0;
      else if (mn == "FixConnected")
        aSFW.FixConnectedMode() = 0;
      else if (mn == "FixEdgeCurves")
        aSFW.FixEdgeCurvesMode() = 0;
      else if (mn == "FixDegenerated")
        aSFW.FixDegeneratedMode() = 0;
      else if (mn == "FixTails")
        aSFW.FixTailMode() = 0;
      else if (mn == "FixShifted")
        aSFW.FixShiftedMode() = 0;
      else if (mn == "FixAddPCurve")
        aSFW.FixAddPCurveMode() = 0;
      else if (mn == "FixAddCurve3d")
        aSFW.FixAddCurve3dMode() = 0;
      else if (mn == "FixReversed2d")
        aSFW.FixReversed2dMode() = 0;
      else if (mn == "FixSeam")
        aSFW.FixSeamMode() = 0;
      else if (mn == "ROff+COn")
      {
        aSFW.FixReorderMode()   = 0;
        aSFW.FixConnectedMode() = 1;
      }
      else if (mn == "ROff+SmallOn")
      {
        aSFW.FixReorderMode() = 0;
        aSFW.FixSmallMode()   = 1;
      }
      const bool aRes = aSFW.Perform(Message_ProgressRange());
      report(mn.c_str(), aSFW.Wire());
      std::cout << "WF   res=" << (aRes ? 1 : 0) << "\n";
    }
    return 0;
  }
  if (argc > 2 && std::string(argv[2]) == "--period")
  {
    int k = 0;
    for (TopExp_Explorer ex(aShape, TopAbs_FACE); ex.More(); ex.Next(), ++k)
    {
      const TopoDS_Face aF = TopoDS::Face(ex.Current());
      BRepAdaptor_Surface aS(aF);
      std::cout << "PER face " << k << " type=" << (int)aS.GetType()
                << " uper=" << (aS.IsUPeriodic() ? 1 : 0) << " vper=" << (aS.IsVPeriodic() ? 1 : 0)
                << " up=" << aS.UPeriod() << " vp=" << aS.VPeriod()
                << " uclo=" << (aS.IsUClosed() ? 1 : 0) << " vclo=" << (aS.IsVClosed() ? 1 : 0)
                << " u=[" << aS.FirstUParameter() << "," << aS.LastUParameter() << "]"
                << " v=[" << aS.FirstVParameter() << "," << aS.LastVParameter() << "]\n";
    }
    return 0;
  }
  if (argc > 2 && std::string(argv[2]) == "--bbox")
  {
    Bnd_Box aBox;
    BRepBndLib::Add(aShape, aBox, false);
    const double aMaxComp = aBox.IsVoid() || aBox.IsOpen()
                              ? 0.0
                              : std::max({aBox.CornerMax().X() - aBox.CornerMin().X(),
                                          aBox.CornerMax().Y() - aBox.CornerMin().Y(),
                                          aBox.CornerMax().Z() - aBox.CornerMin().Z()});
    std::cout << "BBOX maxcomp=" << aMaxComp << " lin=" << aMaxComp * 0.001 * 4.0 << "\n";
    std::printf("BBOX17 maxcomp=%.17e lin=%.17e gap=%.17e\n",
                aMaxComp,
                aMaxComp * 0.001 * 4.0,
                aBox.GetGap());
    double gb[6] = {0, 0, 0, 0, 0, 0};
    aBox.Get(gb[0], gb[1], gb[2], gb[3], gb[4], gb[5]);
    std::printf("BBOX17 get=%.17e %.17e %.17e %.17e %.17e %.17e\n",
                gb[0], gb[1], gb[2], gb[3], gb[4], gb[5]);
    int fi = 0;
    for (TopExp_Explorer fx(aShape, TopAbs_FACE); fx.More(); fx.Next(), ++fi)
    {
      const TopoDS_Face F = TopoDS::Face(fx.Current());
      std::printf("BBOXTOL face=%d ftol=%.17e type=%d\n",
                  fi,
                  BRep_Tool::Tolerance(F),
                  (int)BRepAdaptor_Surface(F).GetType());
      int ei = 0;
      for (TopExp_Explorer ex(F, TopAbs_EDGE); ex.More(); ex.Next(), ++ei)
        std::printf("BBOXTOL   e%d etol=%.17e\n",
                    ei,
                    BRep_Tool::Tolerance(TopoDS::Edge(ex.Current())));
    }
    return 0;
  }

  if (argc > 2 && std::string(argv[2]) == "--entities")
  {
    occ::handle<StepData_StepModel> aModel = aReader.StepModel();
    int aSameT = 0, aSameF = 0;
    for (int i = 1; i <= aModel->NbEntities(); ++i)
    {
      occ::handle<StepShape_FaceSurface> aFS =
        occ::down_cast<StepShape_FaceSurface>(aModel->Value(i));
      if (!aFS.IsNull())
      {
        if (aFS->SameSense()) ++aSameT; else ++aSameF;
        std::cout << "FS ent=" << i << " sameSense=" << (aFS->SameSense() ? 1 : 0) << "\n";
      }
    }
    std::cout << "FS total sameSense T=" << aSameT << " F=" << aSameF << "\n";
    return 0;
  }
  bool perface = false;
  for (int i = 1; i < argc; ++i) { if (std::string(argv[i]) == "--perface") perface = true; }
  if (perface)
  {
    GProp_GProps aP;
    BRepGProp::VolumeProperties(aShape, aP);
    const gp_Pnt aC(0.0, 0.0, 0.0);
    std::cout << "PERFACE total=" << aP.Mass() << " ref=(" << aC.X() << "," << aC.Y() << "," << aC.Z() << ") cm=(" << aP.CentreOfMass().X() << "," << aP.CentreOfMass().Y() << "," << aP.CentreOfMass().Z() << ")\n";
    int k = 0;
    double aSum = 0.0;
    for (TopExp_Explorer ex(aShape, TopAbs_FACE); ex.More(); ex.Next(), ++k)
    {
      const TopoDS_Face aF = TopoDS::Face(ex.Current());
      BRepGProp_Face   aBF;
      aBF.Load(aF);
      BRepGProp_Domain aBD;
      aBD.Init(aF);
      BRepGProp_Vinert aG;
      aG.SetLocation(aC);
      BRepAdaptor_Surface aBAS2(aF);
      aG.Perform(aBF, aBD);
      aSum += aG.Mass();
      std::cout << "PF face " << k << " type=" << (int)aBAS2.GetType() << " orient=" << (int)aF.Orientation() << " contrib=" << aG.Mass();
      if (aBAS2.GetType() == GeomAbs_Sphere) { const gp_Ax3 ax = aBAS2.Sphere().Position(); std::cout << " sphc=(" << ax.Location().X() << "," << ax.Location().Y() << "," << ax.Location().Z() << ")"; }
      std::cout << "\n";
    }
    std::cout << "PF sum=" << aSum << "\n";
    return 0;
  }
  if (argc > 2 && std::string(argv[2]) == "--ds")
  {
    auto report = [](const char* tag, double ax) {
      TopoDS_Shape a = BRepPrimAPI_MakeBox(gp_Pnt(ax, 0., 0.), 1., 1., 1.).Shape();
      TopoDS_Shape b = BRepPrimAPI_MakeBox(gp_Pnt(0.5, 0., 0.), 1., 1., 1.).Shape();
      (void)b;
      TopTools_ListOfShape aLS;
      aLS.Append(a);
      TopoDS_Shape b2 = BRepPrimAPI_MakeBox(gp_Pnt(ax + 3.0, 0., 0.), 1., 1., 1.).Shape();
      aLS.Append(b2);
      BOPAlgo_PaveFiller aPF;
      aPF.SetArguments(aLS);
      aPF.Perform();
      std::cout << "DS " << tag << " NbShapes=" << aPF.DS().NbShapes()
                << " NbSourceShapes=" << aPF.DS().NbSourceShapes() << "\n";
    };
    report("disjoint", 0.0);
    {
      TopoDS_Shape a = BRepPrimAPI_MakeBox(gp_Pnt(0., 0., 0.), 1., 1., 1.).Shape();
      TopoDS_Shape b = BRepPrimAPI_MakeBox(gp_Pnt(0.5, 0., 0.), 1., 1., 1.).Shape();
      TopTools_ListOfShape aLS; aLS.Append(a); aLS.Append(b);
      BOPAlgo_PaveFiller aPF; aPF.SetArguments(aLS); aPF.Perform();
      std::cout << "DS overlap NbShapes=" << aPF.DS().NbShapes()
                << " NbSourceShapes=" << aPF.DS().NbSourceShapes() << "\n";
    }
    return 0;
  }
  if (argc > 2 && std::string(argv[2]) == "--fuse")
  {
    TopoDS_Shape aBox = BRepPrimAPI_MakeBox(gp_Pnt(-1., -1., -1.), gp_Pnt(1., 1., 1.)).Shape();
    TopoDS_Shape aCyl = BRepPrimAPI_MakeCylinder(0.4, 2.0).Shape();
    BRepAlgoAPI_Fuse aFuse(aBox, aCyl);
    aFuse.Build();
    if (!aFuse.IsDone()) { std::cout << "FUSE failed\n"; return 1; }
    TopoDS_Shape aRes = aFuse.Shape();
    GProp_GProps aP; BRepGProp::VolumeProperties(aRes, aP);
    std::cout << "FUSE volume=" << aP.Mass() << " faces=";
    int nf = 0; for (TopExp_Explorer ex(aRes, TopAbs_FACE); ex.More(); ex.Next()) ++nf;
    std::cout << nf << "\n";
    const gp_Pnt aOrig(0., 0., 0.);
    int i = 0;
    for (TopExp_Explorer ex(aRes, TopAbs_FACE); ex.More(); ex.Next(), ++i)
    {
      const TopoDS_Face aF = TopoDS::Face(ex.Current());
      BRepAdaptor_Surface aBAS(aF);
      GProp_GProps aFA; BRepGProp::SurfaceProperties(aF, aFA);
      BRepGProp_Vinert aV; aV.SetLocation(aOrig);
      BRepGProp_Face aBF; aBF.Load(aF); BRepGProp_Domain aBD; aBD.Init(aF); aV.Perform(aBF, aBD);
      std::cout << "FFACE " << i << " type=" << (int)aBAS.GetType() << " orient=" << (int)aF.Orientation()
                << " area=" << aFA.Mass() << " contrib=" << aV.Mass() << "\n";
    }
    return 0;
  }
  if (argc > 2 && std::string(argv[2]) == "--dir2")
  {
    // Faithful re-implementation of BOPTools_AlgoTools.cxx's file-local
    // FindPointInFace (:2160-2231) so the final aDB of GetFaceDir (:2110-2152)
    // can be observed.
    struct FindRes { bool found; gp_Pnt pOut; gp_Dir db; };
    static bool zzVerbose = false;
    auto FindPointInFace = [](const TopoDS_Face& aF, const gp_Pnt& aP, gp_Dir& aDB, gp_Pnt& aPOut,
                              const occ::handle<IntTools_Context>& ctx,
                              GeomAPI_ProjectPointOnSurf& aProjPL, double aDt, double aTolE) -> FindRes {
      double aDTol = Precision::Angular();
      const double aPM = aP.XYZ().Modulus();
      if (aPM > 1000.) aDTol = 5.e-16 * aPM;
      int aNbItMax = 15;
      const double anEps = Precision::SquareConfusion();
      GeomAPI_ProjectPointOnSurf& aProj = ctx->ProjPS(aF);
      gp_Pnt aPS = aP;
      aProj.Perform(aPS);
      if (!aProj.IsDone()) { if (zzVerbose) std::cout << "      early#1 not done\n"; return {false, aPOut, aDB}; }
      aPS = aProj.NearestPoint();
      aProjPL.Perform(aPS);
      aPS = aProjPL.NearestPoint();
      aPS.SetXYZ(aPS.XYZ() + 2. * aTolE * aDB.XYZ());
      aProj.Perform(aPS);
      if (!aProj.IsDone()) { if (zzVerbose) std::cout << "      early#2 not done\n"; return {false, aPOut, aDB}; }
      aPS = aProj.NearestPoint();
      aProjPL.Perform(aPS);
      aPS = aProjPL.NearestPoint();
      bool bRet = false;
      double aDist = 0.;
      do {
        const gp_Pnt aP1(aPS.XYZ() + aDt * aDB.XYZ());
        aProj.Perform(aP1);
        if (!aProj.IsDone()) { if (zzVerbose) std::cout << "      early#3 not done\n"; return {false, aPOut, aDB}; }
        aPOut = aProj.NearestPoint();
        aDist = aProj.LowerDistance();
        aProjPL.Perform(aPOut);
        aPOut = aProjPL.NearestPoint();
        const gp_Vec aV(aPS, aPOut);
        if (zzVerbose)
          std::cout << "      it aDist=" << aDist << " |aV|=" << aV.Magnitude()
                    << " anEps=" << anEps << " dbNow=(" << aDB.X() << "," << aDB.Y() << "," << aDB.Z() << ")\n";
        if (aV.SquareMagnitude() < anEps) return {false, aPOut, aDB};
        aDB.SetXYZ(aV.XYZ());
        aPS = aPOut;
      } while (aDist > aDTol && --aNbItMax);
      bRet = aDist < aDTol;
      return {bRet, aPOut, aDB};
    };

    TopoDS_Shape aBox = BRepPrimAPI_MakeBox(gp_Pnt(-1., -1., -1.), gp_Pnt(1., 1., 1.)).Shape();
    TopoDS_Shape aCyl = BRepPrimAPI_MakeCylinder(0.4, 2.0).Shape();
    TopTools_ListOfShape aArgs;
    aArgs.Append(aBox);
    aArgs.Append(aCyl);
    BOPAlgo_Builder aGF;
    aGF.SetArguments(aArgs);
    aGF.Perform();
    TopoDS_Shape aRes = aGF.Shape();
    TopoDS_Edge aSec;
    for (TopExp_Explorer ex(aRes, TopAbs_EDGE); ex.More(); ex.Next())
    {
      const TopoDS_Edge e = TopoDS::Edge(ex.Current());
      BRepAdaptor_Curve ac(e);
      if (ac.GetType() != GeomAbs_Circle) continue;
      const gp_Circ c = ac.Circle();
      if (std::abs(c.Radius() - 0.4) < 1e-9 && std::abs(c.Location().Z() - 1.0) < 1e-9) { aSec = e; break; }
    }
    if (aSec.IsNull()) { std::cout << "no sect\n"; return 1; }
    NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher> aMEF;
    TopExp::MapShapesAndAncestors(aRes, TopAbs_EDGE, TopAbs_FACE, aMEF);
    const NCollection_List<TopoDS_Shape>& aLF = aMEF.FindFromKey(aSec);
    occ::handle<IntTools_Context> aCtx = new IntTools_Context;
    // current face = first cylinder
    TopoDS_Face f1;
    TopoDS_Edge e1;
    for (NCollection_List<TopoDS_Shape>::Iterator it(aLF); it.More(); it.Next())
    {
      const TopoDS_Face f = TopoDS::Face(it.Value());
      BRepAdaptor_Surface b(f);
      if (b.GetType() == GeomAbs_Cylinder) { f1 = f; BOPTools_AlgoTools::GetEdgeOnFace(aSec, f, e1); break; }
    }
    double t1, t2;
    BRep_Tool::Curve(e1, t1, t2);
    const double t = 0.5 * (t1 + t2);
    gp_Pnt aPx;
    BRepAdaptor_Curve ac1(e1);
    ac1.D0(t, aPx);
    gp_Vec tau1;
    BOPTools_AlgoTools2D::EdgeTangent(e1, t, tau1);
    const gp_Dir aDTgt(tau1);
    gp_Dir aDN1;
    BOPTools_AlgoTools3D::GetNormalToFaceOnEdge(e1, f1, t, aDN1, aCtx);
    gp_Dir aDBF(aDN1.XYZ().Crossed(aDTgt.XYZ()));
    // aProjPL: plane through aPx with normal aDTgt
    occ::handle<Geom_Plane> aPL = new Geom_Plane(aPx, aDTgt);
    double u1, u2, v1, v2;
    aPL->Bounds(u1, u2, v1, v2);
    GeomAPI_ProjectPointOnSurf aProjPL;
    aProjPL.Init(aPL, u1, u2, v1, v2);
    const double tolE = BRep_Tool::Tolerance(e1);
    // Faithful MinStep3D (BOPTools_AlgoTools.cxx:2235-2346).
    double dtMax = -1.;
    double dtMin = 5.e-6;
    for (NCollection_List<TopoDS_Shape>::Iterator it(aLF); it.More(); it.Next())
    {
      const TopoDS_Face f = TopoDS::Face(it.Value());
      dtMax = std::max(dtMax, 2. * (tolE + BRep_Tool::Tolerance(f)));
      BRepAdaptor_Surface aBAS(f);
      double aR = 0.;
      switch (aBAS.GetType())
      {
        case GeomAbs_Cylinder: aR = aBAS.Cylinder().Radius(); break;
        case GeomAbs_Cone: { gp_Lin aL(aBAS.Cone().Axis()); aR = aL.Distance(aPx); } break;
        case GeomAbs_Sphere: dtMin = std::max(dtMin, 5.e-4); aR = aBAS.Sphere().Radius(); break;
        case GeomAbs_Torus: aR = aBAS.Torus().MajorRadius(); break;
        default: dtMin = std::max(dtMin, 5.e-4); break;
      }
      if (aR > 100.)
      {
        constexpr double d = 10 * Precision::PConfusion();
        dtMin = std::max(dtMin, std::sqrt(d * d + 2 * d * aR));
      }
    }
    if (dtMax < dtMin) dtMax = dtMin;
    const double aDt = dtMax;
    gp_Pnt p1 = aPx;
    gp_Dir db1 = aDBF;
    zzVerbose = true;
    FindRes r1 = FindPointInFace(f1, aPx, db1, p1, aCtx, aProjPL, aDt, tolE);
    zzVerbose = false;
    if (r1.found) {
      aDBF = db1;
    } else {
      // GetFaceDir's fallback (BOPTools_AlgoTools.cxx:2139-2149)
      gp_Pnt aPx2;
      gp_Dir aDN2f;
      BOPTools_AlgoTools3D::GetApproxNormalToFaceOnEdge(e1, f1, t, aDt, aPx2, aDN2f, aCtx);
      aProjPL.Perform(aPx2);
      aPx2 = aProjPL.NearestPoint();
      const gp_Vec aVec(aPx, aPx2);
      aDBF.SetXYZ(aVec.XYZ());
      aDN1 = aDN2f;
    }
    const gp_Dir aDTF(aDN1.XYZ().Crossed(aDBF.XYZ()));
    BRepAdaptor_Surface b1(f1);
    std::cout << "CUR type=" << (int)b1.GetType() << " eo=" << (int)e1.Orientation()
              << " dt=" << aDt << " found1=" << (r1.found ? 1 : 0)
              << " dbf=(" << aDBF.X() << "," << aDBF.Y() << "," << aDBF.Z() << ")"
              << " dtf=(" << aDTF.X() << "," << aDTF.Y() << "," << aDTF.Z() << ")\n";
    for (NCollection_List<TopoDS_Shape>::Iterator it(aLF); it.More(); it.Next())
    {
      const TopoDS_Face f2 = TopoDS::Face(it.Value());
      if (f2.IsSame(f1)) continue;
      TopoDS_Edge e2;
      if (!BOPTools_AlgoTools::GetEdgeOnFace(aSec, f2, e2)) continue;
      const gp_Dir aDTgt2 = (e2.Orientation() == e1.Orientation()) ? aDTgt : aDTgt.Reversed();
      gp_Dir aDN2;
      BOPTools_AlgoTools3D::GetNormalToFaceOnEdge(e2, f2, t, aDN2, aCtx);
      gp_Dir db2(aDN2.XYZ().Crossed(aDTgt2.XYZ()));
      const gp_Dir db2_init = db2;
      gp_Pnt p2 = aPx;
      FindRes r2 = FindPointInFace(f2, aPx, db2, p2, aCtx, aProjPL, aDt, tolE);
      const double aAngle = aDBF.AngleWithRef(db2, aDTF) * 180. / M_PI;
      BRepAdaptor_Surface b2(f2);
      std::cout << "  CAND type=" << (int)b2.GetType() << " eo=" << (int)e2.Orientation()
                << " found=" << (r2.found ? 1 : 0)
                << " db_init=(" << db2_init.X() << "," << db2_init.Y() << "," << db2_init.Z() << ")"
                << " db_fin=(" << db2.X() << "," << db2.Y() << "," << db2.Z() << ")"
                << " angle=" << aAngle << "\n";
    }
    return 0;
  }
  if (argc > 2 && std::string(argv[2]) == "--dir")
  {
    TopoDS_Shape aBox = BRepPrimAPI_MakeBox(gp_Pnt(-1., -1., -1.), gp_Pnt(1., 1., 1.)).Shape();
    TopoDS_Shape aCyl = BRepPrimAPI_MakeCylinder(0.4, 2.0).Shape();
    TopTools_ListOfShape aArgs;
    aArgs.Append(aBox);
    aArgs.Append(aCyl);
    BOPAlgo_Builder aGF;
    aGF.SetArguments(aArgs);
    aGF.Perform();
    TopoDS_Shape aRes = aGF.Shape();
    TopoDS_Edge aSec;
    for (TopExp_Explorer ex(aRes, TopAbs_EDGE); ex.More(); ex.Next())
    {
      const TopoDS_Edge e = TopoDS::Edge(ex.Current());
      BRepAdaptor_Curve ac(e);
      if (ac.GetType() != GeomAbs_Circle) continue;
      const gp_Circ c = ac.Circle();
      if (std::abs(c.Radius() - 0.4) < 1e-9 && std::abs(c.Location().Z() - 1.0) < 1e-9) { aSec = e; break; }
    }
    if (aSec.IsNull()) { std::cout << "no sect\n"; return 1; }
    NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher> aMEF;
    TopExp::MapShapesAndAncestors(aRes, TopAbs_EDGE, TopAbs_FACE, aMEF);
    const NCollection_List<TopoDS_Shape>& aLF = aMEF.FindFromKey(aSec);
    occ::handle<IntTools_Context> aCtx = new IntTools_Context;
    int k = 0;
    for (NCollection_List<TopoDS_Shape>::Iterator it(aLF); it.More(); it.Next(), ++k)
    {
      const TopoDS_Face f = TopoDS::Face(it.Value());
      TopoDS_Edge eo;
      if (!BOPTools_AlgoTools::GetEdgeOnFace(aSec, f, eo)) continue;
      double t1, t2;
      BRep_Tool::Curve(eo, t1, t2);
      const double t = 0.5 * (t1 + t2);
      gp_Dir aD;
      BOPTools_AlgoTools3D::GetNormalToFaceOnEdge(eo, f, t, aD, aCtx);
      gp_Vec tau;
      BOPTools_AlgoTools2D::EdgeTangent(eo, t, tau);
      gp_Dir aDTgt(tau);
      const gp_Vec db = aD.XYZ().Crossed(aDTgt.XYZ());
      BRepAdaptor_Surface b(f);
      std::cout << "DIR f" << k << " type=" << (int)b.GetType() << " orient=" << (int)f.Orientation()
                << " eo=" << (int)eo.Orientation()
                << " dN=(" << aD.X() << "," << aD.Y() << "," << aD.Z() << ")"
                << " dTgt=(" << aDTgt.X() << "," << aDTgt.Y() << "," << aDTgt.Z() << ")"
                << " dB=(" << db.X() << "," << db.Y() << "," << db.Z() << ")\n";
    }
    return 0;
  }
  if (argc > 2 && std::string(argv[2]) == "--faceoff2")
  {
    // Same as --faceoff, but builds the candidate couples exactly the way
    // BOPAlgo_ShellSplitter::SplitBlock does: BOPTools_AlgoTools::GetEdgeOff
    // (BOPTools_AlgoTools.cxx:1099-1127) insists on the edge view whose
    // orientation is the OPPOSITE of the current one, and skips the candidate
    // otherwise. Compare the two modes to see whether the couple construction
    // decides the answer.
    TopoDS_Shape aBox = BRepPrimAPI_MakeBox(gp_Pnt(-1., -1., -1.), gp_Pnt(1., 1., 1.)).Shape();
    TopoDS_Shape aCyl = BRepPrimAPI_MakeCylinder(0.4, 2.0).Shape();
    TopTools_ListOfShape aArgs;
    aArgs.Append(aBox);
    aArgs.Append(aCyl);
    BOPAlgo_Builder aGF;
    aGF.SetArguments(aArgs);
    aGF.Perform();
    TopoDS_Shape aRes = aGF.Shape();
    TopoDS_Edge aSec;
    for (TopExp_Explorer ex(aRes, TopAbs_EDGE); ex.More(); ex.Next())
    {
      const TopoDS_Edge e = TopoDS::Edge(ex.Current());
      BRepAdaptor_Curve ac(e);
      if (ac.GetType() != GeomAbs_Circle) continue;
      const gp_Circ c = ac.Circle();
      if (std::abs(c.Radius() - 0.4) < 1e-9 && std::abs(c.Location().Z() - 1.0) < 1e-9) { aSec = e; break; }
    }
    if (aSec.IsNull()) { std::cout << "no section circle\n"; return 1; }
    NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher> aMEF;
    TopExp::MapShapesAndAncestors(aRes, TopAbs_EDGE, TopAbs_FACE, aMEF);
    const NCollection_List<TopoDS_Shape>& aLF = aMEF.FindFromKey(aSec);
    occ::handle<IntTools_Context> aCtx = new IntTools_Context;
    int k = 0;
    for (NCollection_List<TopoDS_Shape>::Iterator it(aLF); it.More(); it.Next(), ++k)
    {
      const TopoDS_Face f1 = TopoDS::Face(it.Value());
      TopoDS_Edge e1;
      if (!BOPTools_AlgoTools::GetEdgeOnFace(aSec, f1, e1)) continue;
      NCollection_List<BOPTools_CoupleOfShape> aLC;
      int skipped = 0, kept = 0;
      for (NCollection_List<TopoDS_Shape>::Iterator jt(aLF); jt.More(); jt.Next())
      {
        const TopoDS_Face f2 = TopoDS::Face(jt.Value());
        if (f2.IsSame(f1)) continue;
        TopoDS_Edge e2;
        if (!BOPTools_AlgoTools::GetEdgeOff(e1, f2, e2)) { ++skipped; continue; }
        BOPTools_CoupleOfShape cs;
        cs.SetShape1(e2);
        cs.SetShape2(f2);
        aLC.Append(cs);
        ++kept;
      }
      TopoDS_Face aOff;
      const bool done = BOPTools_AlgoTools::GetFaceOff(e1, f1, aLC, aOff, aCtx);
      BRepAdaptor_Surface b1(f1);
      if (!aOff.IsNull())
      {
        BRepAdaptor_Surface b2(aOff);
        std::cout << "O2 from f" << k << "(type=" << (int)b1.GetType() << " eo=" << (int)e1.Orientation()
                  << ") couples=" << kept << "/skip=" << skipped << " -> type=" << (int)b2.GetType()
                  << " orient=" << (int)aOff.Orientation() << " done=" << (done ? 1 : 0) << "\n";
      }
      else
      {
        std::cout << "O2 from f" << k << "(type=" << (int)b1.GetType() << " eo=" << (int)e1.Orientation()
                  << ") couples=" << kept << "/skip=" << skipped << " -> null done=" << (done ? 1 : 0) << "\n";
      }
    }
    return 0;
  }
  if (argc > 2 && std::string(argv[2]) == "--faceoff")
  {
    TopoDS_Shape aBox = BRepPrimAPI_MakeBox(gp_Pnt(-1., -1., -1.), gp_Pnt(1., 1., 1.)).Shape();
    TopoDS_Shape aCyl = BRepPrimAPI_MakeCylinder(0.4, 2.0).Shape();
    TopTools_ListOfShape aArgs;
    aArgs.Append(aBox);
    aArgs.Append(aCyl);
    BOPAlgo_Builder aGF;
    aGF.SetArguments(aArgs);
    aGF.Perform();
    if (aGF.HasErrors()) { std::cout << "GF failed\n"; return 1; }
    TopoDS_Shape aRes = aGF.Shape();
    // section circle: radius 0.4 at z=1
    TopoDS_Edge aSec;
    for (TopExp_Explorer ex(aRes, TopAbs_EDGE); ex.More(); ex.Next())
    {
      const TopoDS_Edge e = TopoDS::Edge(ex.Current());
      BRepAdaptor_Curve ac(e);
      if (ac.GetType() != GeomAbs_Circle) continue;
      const gp_Circ c = ac.Circle();
      if (std::abs(c.Radius() - 0.4) < 1e-9 && std::abs(c.Location().Z() - 1.0) < 1e-9)
      {
        aSec = e;
        break;
      }
    }
    if (aSec.IsNull()) { std::cout << "no section circle\n"; return 1; }
    NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher> aMEF;
    TopExp::MapShapesAndAncestors(aRes, TopAbs_EDGE, TopAbs_FACE, aMEF);
    if (!aMEF.Contains(aSec)) { std::cout << "edge not in map\n"; return 1; }
    const NCollection_List<TopoDS_Shape>& aLF = aMEF.FindFromKey(aSec);
    std::cout << "SEC eo=" << (int)aSec.Orientation() << " adj_faces=" << aLF.Extent() << "\n";
    NCollection_List<BOPTools_CoupleOfShape> aLCS;
    for (NCollection_List<TopoDS_Shape>::Iterator it(aLF); it.More(); it.Next())
    {
      const TopoDS_Face f = TopoDS::Face(it.Value());
      BRepAdaptor_Surface bas(f);
      TopoDS_Edge eo;
      const bool ok = BOPTools_AlgoTools::GetEdgeOnFace(aSec, f, eo);
      std::cout << "  face type=" << (int)bas.GetType() << " orient=" << (int)f.Orientation()
                << " eo=" << (ok ? (int)eo.Orientation() : -1) << " ok=" << (ok ? 1 : 0) << "\n";
      if (ok)
      {
        BOPTools_CoupleOfShape cs;
        cs.SetShape1(eo);
        cs.SetShape2(f);
        aLCS.Append(cs);
      }
    }
    // call GetFaceOff from each adjacent face
    occ::handle<IntTools_Context> aCtx = new IntTools_Context;
    int k = 0;
    for (NCollection_List<TopoDS_Shape>::Iterator it(aLF); it.More(); it.Next(), ++k)
    {
      const TopoDS_Face f1 = TopoDS::Face(it.Value());
      TopoDS_Edge e1;
      if (!BOPTools_AlgoTools::GetEdgeOnFace(aSec, f1, e1)) continue;
      NCollection_List<BOPTools_CoupleOfShape> aLC;
      for (NCollection_List<BOPTools_CoupleOfShape>::Iterator jt(aLCS); jt.More(); jt.Next())
      {
        if (jt.Value().Shape2().IsSame(f1)) continue;
        BOPTools_CoupleOfShape cs;
        cs.SetShape1(jt.Value().Shape1());
        cs.SetShape2(jt.Value().Shape2());
        aLC.Append(cs);
      }
      TopoDS_Face aOff;
      const bool done = BOPTools_AlgoTools::GetFaceOff(e1, f1, aLC, aOff, aCtx);
      BRepAdaptor_Surface b1(f1);
      if (!aOff.IsNull())
      {
        BRepAdaptor_Surface b2(aOff);
        std::cout << "OFF from f" << k << "(type=" << (int)b1.GetType() << ") -> type="
                  << (int)b2.GetType() << " orient=" << (int)aOff.Orientation()
                  << " done=" << (done ? 1 : 0) << "\n";
      }
      else
      {
        std::cout << "OFF from f" << k << "(type=" << (int)b1.GetType() << ") -> null done="
                  << (done ? 1 : 0) << "\n";
      }
    }
    return 0;
  }
  if (argc > 2 && std::string(argv[2]) == "--fcls2d")
  {
    // Ground truth for `BRepTopAdaptor_FClass2d` on the cylinder lateral face of
    // `BRepPrimAPI_MakeCylinder(0.4, 2.0)`: the `BRepTools_WireExplorer` order
    // (which pcurve each seam occurrence walks), the classifier state over a u
    // grid at v=1, and the `IntTools_Context::IsPointInOnFace` path for the
    // pave point (0,-0.4,1).
    TopoDS_Shape aPrim = BRepPrimAPI_MakeCylinder(0.4, 2.0).Shape();
    for (TopExp_Explorer ex(aPrim, TopAbs_FACE); ex.More(); ex.Next())
    {
      TopoDS_Face aF = TopoDS::Face(ex.Current());
      BRepAdaptor_Surface aBAS(aF, false);
      if (aBAS.GetType() != GeomAbs_Cylinder) continue;
      aF.Orientation(TopAbs_FORWARD);
      std::cout << "FCLS2D cylinder face\n";
      for (TopExp_Explorer exw(aF, TopAbs_WIRE); exw.More(); exw.Next())
      {
        const TopoDS_Wire       aW = TopoDS::Wire(exw.Current());
        BRepTools_WireExplorer  aWE;
        int                     k = 0;
        for (aWE.Init(aW, aF); aWE.More(); aWE.Next(), ++k)
        {
          const TopoDS_Edge         aE = aWE.Current();
          double                    f = 0., l = 0.;
          occ::handle<Geom2d_Curve> aC = BRep_Tool::CurveOnSurface(aE, aF, f, l);
          gp_Pnt2d                  p1(99., 99.), p2(99., 99.), pv(99., 99.);
          if (!aC.IsNull())
          {
            aC->D0(f, p1);
            aC->D0(l, p2);
          }
          const TopoDS_Vertex aV = aWE.CurrentVertex();
          if (!aV.IsNull() && !aC.IsNull())
          {
            aC->D0(BRep_Tool::Parameter(aV, aE, aF), pv);
          }
          std::cout << "  WE" << k << " ori=" << (int)aWE.Orientation()
                    << " closed=" << (BRep_Tool::IsClosed(aE, aF) ? 1 : 0) << " range=[" << f << ","
                    << l << "] pc_f=(" << p1.X() << "," << p1.Y() << ") pc_l=(" << p2.X() << ","
                    << p2.Y() << ") curV=(" << pv.X() << "," << pv.Y() << ")\n";
        }
      }
      BRepTopAdaptor_FClass2d aFc(aF, Precision::Confusion());
      BRepClass_FaceClassifier aCls;
      for (int i = 0; i <= 16; ++i)
      {
        const gp_Pnt2d aP(2.0 * M_PI * i / 16.0, 1.0);
        aCls.Perform(aF, aP, Precision::Confusion());
        std::cout << "  P u=" << aP.X() << " fclass=" << (int)aFc.Perform(aP)
                  << " faceclass=" << (int)aCls.State() << "\n";
      }
      GeomAPI_ProjectPointOnSurf aProj(gp_Pnt(0.0, -0.4, 1.0), BRep_Tool::Surface(aF));
      if (aProj.IsDone())
      {
        double u = 0., v = 0.;
        aProj.LowerDistanceParameters(u, v);
        const gp_Pnt2d aP2d(u, v);
        aCls.Perform(aF, aP2d, Precision::Confusion());
        std::cout << "  ISPOINT uv=(" << u << "," << v << ") fclass=" << (int)aFc.Perform(aP2d)
                  << " faceclass=" << (int)aCls.State() << "\n";
      }
      std::cout << "  InfPoint fclass=" << (int)aFc.PerformInfinitePoint() << "\n";
    }
    return 0;
  }
  if (argc > 2 && std::string(argv[2]) == "--cylall")
  {
    TopoDS_Shape aPrim = BRepPrimAPI_MakeCylinder(0.4, 2.0).Shape();
    int ie = 0;
    for (TopExp_Explorer ex(aPrim, TopAbs_FACE); ex.More(); ex.Next(), ++ie)
    {
      const TopoDS_Face aF = TopoDS::Face(ex.Current());
      BRepAdaptor_Surface aBAS(aF);
      std::cout << "FACE " << ie << " type=" << (int)aBAS.GetType()
                << " orient=" << (int)aF.Orientation() << "\n";
      for (TopExp_Explorer exe(aF, TopAbs_EDGE); exe.More(); exe.Next())
      {
        const TopoDS_Edge aE = TopoDS::Edge(exe.Current());
        double f = 0., l = 0.;
        occ::handle<Geom2d_Curve> aC = BRep_Tool::CurveOnSurface(aE, aF, f, l);
        gp_Pnt2d pm(99., 99.);
        if (!aC.IsNull()) aC->D0(0.5 * (f + l), pm);
        std::cout << "  E ori=" << (int)aE.Orientation() << " uv_mid=(" << pm.X() << "," << pm.Y()
                  << ")\n";
      }
    }
    return 0;
  }
  if (argc > 2 && std::string(argv[2]) == "--cylface")
  {
    TopoDS_Shape aPrim = BRepPrimAPI_MakeCylinder(0.4, 2.0).Shape();
    for (TopExp_Explorer ex(aPrim, TopAbs_FACE); ex.More(); ex.Next())
    {
      const TopoDS_Face aF = TopoDS::Face(ex.Current());
      BRepAdaptor_Surface aBAS(aF);
      if (aBAS.GetType() != GeomAbs_Cylinder) continue;
      std::cout << "CYLFACE orient=" << (int)aF.Orientation() << "\n";
      int ie = 0;
      for (TopExp_Explorer exe(aF, TopAbs_EDGE); exe.More(); exe.Next(), ++ie)
      {
        const TopoDS_Edge aE = TopoDS::Edge(exe.Current());
        double f = 0., l = 0.;
        occ::handle<Geom2d_Curve> aC = BRep_Tool::CurveOnSurface(aE, aF, f, l);
        std::cout << "  E" << ie << " ori=" << (int)aE.Orientation()
                  << " closed=" << (BRep_Tool::IsClosed(aE, aF) ? 1 : 0)
                  << " range=[" << f << "," << l << "] null=" << (aC.IsNull() ? 1 : 0);
        if (!aC.IsNull())
        {
          gp_Pnt2d pm;
          aC->D0(0.5 * (f + l), pm);
          std::cout << " mid=(" << pm.X() << "," << pm.Y() << ")";
        }
        std::cout << "\n";
        TopoDS_Iterator it(aE);
        int iv = 0;
        for (; it.More(); it.Next(), ++iv)
        {
          const TopoDS_Vertex aV = TopoDS::Vertex(it.Value());
          const double        t  = BRep_Tool::Parameter(aV, aE, aF);
          gp_Pnt2d            pv(99., 99.);
          if (!aC.IsNull()) aC->D0(t, pv);
          std::cout << "    V" << iv << " ori=" << (int)aV.Orientation() << " param=" << t
                    << " uv=(" << pv.X() << "," << pv.Y() << ")\n";
        }
        // Static replay of BOPAlgo_WireSplitter::GetNextVertex for each child
        // forced FORWARD (the Path call passes pVa.Orientation(FORWARD)).
        for (int k = 0; k < 2; ++k)
        {
          TopoDS_Iterator itk(aE);
          int             m = 0;
          for (; itk.More() && m < k; itk.Next(), ++m) {}
          if (!itk.More()) continue;
          TopoDS_Vertex aV = TopoDS::Vertex(itk.Value());
          aV.Orientation(TopAbs_FORWARD);
          TopoDS_Iterator itv(aE);
          TopoDS_Vertex   aV1;
          for (; itv.More(); itv.Next())
          {
            const TopoDS_Shape& aVx = itv.Value();
            if (!aVx.IsEqual(aV))
            {
              aV1 = TopoDS::Vertex(aVx);
              break;
            }
          }
          double   tn = BRep_Tool::Parameter(aV1, aE, aF);
          gp_Pnt2d pn(99., 99.);
          if (!aC.IsNull()) aC->D0(tn, pn);
          std::cout << "    NEXT(from V" << k << " forced F) -> ori=" << (int)aV1.Orientation()
                    << " param=" << tn << " uv=(" << pn.X() << "," << pn.Y() << ")\n";
        }
      }
    }
    return 0;
  }
  if (argc > 2 && (std::string(argv[2]) == "--cyl" || std::string(argv[2]) == "--sph"))
  {
    TopoDS_Shape aPrim;
    if (std::string(argv[2]) == "--cyl")
      aPrim = BRepPrimAPI_MakeCylinder(0.4, 2.0).Shape();
    else
      aPrim = BRepPrimAPI_MakeSphere(1.0).Shape();
    GProp_GProps aP;
    BRepGProp::VolumeProperties(aPrim, aP);
    std::cout << "PRIM volume=" << aP.Mass() << " faces=";
    int nf = 0;
    for (TopExp_Explorer ex(aPrim, TopAbs_FACE); ex.More(); ex.Next()) ++nf;
    std::cout << nf << "\n";
    int i = 0;
    for (TopExp_Explorer ex(aPrim, TopAbs_FACE); ex.More(); ex.Next(), ++i)
    {
      const TopoDS_Face aF = TopoDS::Face(ex.Current());
      BRepAdaptor_Surface aBAS(aF);
      BRepGProp_Face aBF; aBF.Load(aF);
      double u1,u2,v1,v2; aBF.Bounds(u1,u2,v1,v2);
      std::cout << "PRIMFACE " << i << " type=" << (int)aBAS.GetType()
                << " orient=" << (int)aF.Orientation() << " uv=[" << u1 << "," << u2 << "]x["
                << v1 << "," << v2 << "]\n";
      int j = 0;
      for (TopExp_Explorer exe(aF, TopAbs_EDGE); exe.More(); exe.Next(), ++j)
      {
        const TopoDS_Edge aE = TopoDS::Edge(exe.Current());
        if (!aBF.Load(aE)) { std::cout << "  E" << j << " load=false\n"; continue; }
        const double lp1 = aBF.FirstParameter(), lp2 = aBF.LastParameter();
        const double lm = 0.5 * (lp1 + lp2);
        gp_Pnt2d Puv; gp_Vec2d Vuv;
        aBF.D12d(lm, Puv, Vuv);
        std::cout << "  E" << j << " eo=" << (int)aE.Orientation() << " range=[" << lp1 << ","
                  << lp2 << "] midU=" << Puv.X() << " midV=" << Puv.Y() << " du=" << Vuv.X()
                  << " dv=" << Vuv.Y() << "\n";
      }
    }
    return 0;
  }
  if (argc > 3 && std::string(argv[2]) == "--facek")
  {
    const int k = std::atoi(argv[3]);
    int i = 0;
    for (TopExp_Explorer ex(aShape, TopAbs_FACE); ex.More(); ex.Next(), ++i)
    {
      if (i != k) continue;
      const TopoDS_Face aF = TopoDS::Face(ex.Current());
      BRepAdaptor_Surface aBAS(aF);
      BRepGProp_Face aBF; aBF.Load(aF);
      double u1,u2,v1,v2; aBF.Bounds(u1,u2,v1,v2);
      std::cout << "FK face=" << k << " type=" << (int)aBAS.GetType() << " orient=" << (int)aF.Orientation()
                << " uv=[" << u1 << "," << u2 << "]x[" << v1 << "," << v2 << "]\n";
      int j = 0;
      for (TopExp_Explorer exe(aF, TopAbs_EDGE); exe.More(); exe.Next(), ++j)
      {
        const TopoDS_Edge aE = TopoDS::Edge(exe.Current());
        if (!aBF.Load(aE)) { std::cout << "  FKE" << j << " load=false\n"; continue; }
        const double lp1 = aBF.FirstParameter(), lp2 = aBF.LastParameter();
        const double lm = 0.5 * (lp1 + lp2);
        gp_Pnt2d Puv; gp_Vec2d Vuv; aBF.D12d(lm, Puv, Vuv);
        std::cout << "  FKE" << j << " eo=" << (int)aE.Orientation() << " range=[" << lp1 << "," << lp2
                  << "] midU=" << Puv.X() << " midV=" << Puv.Y() << " du=" << Vuv.X() << " dv=" << Vuv.Y() << "\n";
      }
      return 0;
    }
    return 0;
  }
  int aNbFaces = 0;
  for (TopExp_Explorer ex(aShape, TopAbs_FACE); ex.More(); ex.Next())
    ++aNbFaces;
  std::cout << "SHAPE null=0 faces=" << aNbFaces << " rootOrient=" << (int)aShape.Orientation()
            << " rootType=" << (int)aShape.ShapeType() << "\n";
  for (TopExp_Explorer ex(aShape, TopAbs_SOLID); ex.More(); ex.Next())
    std::cout << "SOLID orient=" << (int)ex.Current().Orientation() << "\n";
  for (TopExp_Explorer ex(aShape, TopAbs_SHELL); ex.More(); ex.Next())
    std::cout << "SHELL orient=" << (int)ex.Current().Orientation() << "\n";
  {
    int j = 0;
    for (TopExp_Explorer exs(aShape, TopAbs_SHELL); exs.More(); exs.Next())
    {
      std::cout << "SHELLCHILD shellOrient=" << (int)exs.Current().Orientation() << "\n";
      // cumOri=false: the *raw* stored orientation of each face, without
      // composing the shell's own orientation.
      for (TopoDS_Iterator it(exs.Current(), Standard_False, Standard_False); it.More();
           it.Next(), ++j)
      {
        std::cout << "RAW ptr=" << (const void*)it.Value().TShape().get()
                  << " orient=" << (int)it.Value().Orientation() << "\n";
      }
    }
    for (TopExp_Explorer exf(aShape, TopAbs_FACE); exf.More(); exf.Next())
    {
      std::cout << "EXPLORER ptr=" << (const void*)exf.Current().TShape().get()
                << " orient=" << (int)exf.Current().Orientation() << "\n";
    }
  }

  BRepCheck_Analyzer aCheck(aShape);
  std::cout << "BRepCheck valid=" << (aCheck.IsValid() ? 1 : 0) << "\n";

  GProp_GProps aProps;
  BRepGProp::VolumeProperties(aShape, aProps);
  std::cout.setf(std::ios::fixed);
  std::cout.precision(6);
  std::cout << "BRepGProp volume=" << aProps.Mass() << "\n";
  const gp_Pnt aCenter = aProps.CentreOfMass();
  std::cout << "BRepGProp centre=(" << aCenter.X() << "," << aCenter.Y() << "," << aCenter.Z()
            << ")\n";

  int i = 0;
  for (TopExp_Explorer ex(aShape, TopAbs_FACE); ex.More(); ex.Next(), ++i)
  {
    const TopoDS_Face aFace = TopoDS::Face(ex.Current());
    BRepAdaptor_Surface aBAS(aFace);
    const GeomAbs_SurfaceType aType = aBAS.GetType();

    BRepGProp_Face aBF;
    aBF.Load(aFace);
    double u1 = 0., u2 = 0., v1 = 0., v2 = 0.;
    aBF.Bounds(u1, u2, v1, v2);
    gp_Pnt aP;
    gp_Vec aN;
    aBF.Normal(0.5 * (u1 + u2), 0.5 * (v1 + v2), aP, aN);
    const double anOutward = gp_Vec(aCenter, aP).Dot(aN);
    // Raw (unflipped) natural normal: D1U ^ D1V at the same parameters.
    gp_Pnt aP2;
    gp_Vec aDU, aDV;
    aBAS.D1(0.5 * (u1 + u2), 0.5 * (v1 + v2), aP2, aDU, aDV);
    const gp_Vec aRaw = aDU.Crossed(aDV);
    const double anOutwardRaw = gp_Vec(aCenter, aP2).Dot(aRaw);

    const char* aOri = "?";
    switch (aFace.Orientation())
    {
      case TopAbs_FORWARD: aOri = "FWD"; break;
      case TopAbs_REVERSED: aOri = "REV"; break;
      case TopAbs_INTERNAL: aOri = "INT"; break;
      case TopAbs_EXTERNAL: aOri = "EXT"; break;
    }
    std::cout << "face " << i << " type=" << (int)aType << " orient=" << aOri
              << " outw=" << anOutward << " raw=" << anOutwardRaw << " uv=[" << u1 << "," << u2 << "]x[" << v1 << "," << v2
              << "]";
    if (frames)
    {
      if (aType == GeomAbs_Cylinder)
      {
        const gp_Cylinder aCyl = aBAS.Cylinder();
        const gp_Ax3 aAx = aCyl.Position();
        std::cout << " cyl loc=(" << aAx.Location().X() << "," << aAx.Location().Y() << ","
                  << aAx.Location().Z() << ") Z=(" << aAx.Direction().X() << ","
                  << aAx.Direction().Y() << "," << aAx.Direction().Z() << ") X=("
                  << aAx.XDirection().X() << "," << aAx.XDirection().Y() << ","
                  << aAx.XDirection().Z() << ") Y=(" << aAx.YDirection().X() << ","
                  << aAx.YDirection().Y() << "," << aAx.YDirection().Z()
                  << ") R=" << aCyl.Radius();
      }
      else if (aType == GeomAbs_Plane)
      {
        const gp_Pln aPln = aBAS.Plane();
        const gp_Ax3 aAx = aPln.Position();
        std::cout << " pln loc=(" << aAx.Location().X() << "," << aAx.Location().Y() << ","
                  << aAx.Location().Z() << ") Z=(" << aAx.Direction().X() << ","
                  << aAx.Direction().Y() << "," << aAx.Direction().Z() << ") X=("
                  << aAx.XDirection().X() << "," << aAx.XDirection().Y() << ","
                  << aAx.XDirection().Z() << ")";
      }
      else if (aType == GeomAbs_Sphere)
      {
        const gp_Sphere aSph = aBAS.Sphere();
        const gp_Ax3   aAx  = aSph.Position();
        std::cout << " sph loc=(" << aAx.Location().X() << "," << aAx.Location().Y() << ","
                  << aAx.Location().Z() << ") Z=(" << aAx.Direction().X() << ","
                  << aAx.Direction().Y() << "," << aAx.Direction().Z() << ") R=" << aSph.Radius();
      }
    }
    std::cout << "\n";
    // Per-face wire edges with their composed orientation and endpoints, so the
    // port's edge-orientation mapping can be diffed edge by edge.
    {
      int ei = 0;
      for (TopExp_Explorer exw(aFace, TopAbs_WIRE); exw.More(); exw.Next())
      {
        for (TopExp_Explorer exe(exw.Current(), TopAbs_EDGE); exe.More(); exe.Next(), ++ei)
        {
          const TopoDS_Edge anE = TopoDS::Edge(exe.Current());
          double f = 0., l = 0.;
          const occ::handle<Geom_Curve>& aC = BRep_Tool::Curve(anE, f, l);
          TopoDS_Vertex aV1, aV2;
          TopExp::Vertices(anE, aV1, aV2);
          const gp_Pnt p1 = BRep_Tool::Pnt(aV1);
          const gp_Pnt p2 = BRep_Tool::Pnt(aV2);
          std::cout << "EDGE face=" << i << " n=" << ei
                    << " orient=" << (int)anE.Orientation() << " null3d=" << (aC.IsNull() ? 1 : 0)
                    << " p1=(" << p1.X() << "," << p1.Y() << "," << p1.Z() << ")"
                    << " p2=(" << p2.X() << "," << p2.Y() << "," << p2.Z() << ")\n";
        }
      }
    }
  }
  return 0;
}
