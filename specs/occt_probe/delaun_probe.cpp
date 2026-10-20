// Delaunay state probes: `--delauncheck`, `--delaunfeed`, `--delaunstruct`,
// `--delaunsteps`, `--delaunwhy`.
#include "occt_probe_common.h"

int runDelaunCheck()
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
int runDelaunFeed(const std::string& theDir)
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
int runDelaunStruct(const std::string& theDir)
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
int runDelaunSteps(const std::string& theDir, const int theFace)
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
int runDelaunWhy(const std::string& theDir, const int theFace, const int theStep)
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
