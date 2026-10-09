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
#include <IGESControl_Controller.hxx>
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
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2dAdaptor_Curve.hxx>
#include <gp_Pnt2d.hxx>
#include <algorithm>
#include <cmath>
#include <map>
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
  const int aFaces[] = {0,   208, 209, 169, 174, 189, 192, 193, 194,
                        195, 196, 197, 198, 199, 204, 206, 210, 212, 433,
                        51,  52,  53,  54,  321};
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

    occ::handle<NCollection_IncAllocator>       anAlloc = new NCollection_IncAllocator(IMeshData::MEMORY_BLOCK_SIZE_HUGE);
    occ::handle<BRepMesh_DataStructureOfDelaun> aStruct =
      new BRepMesh_DataStructureOfDelaun(anAlloc, static_cast<int>(aNodes.size()));
    aStruct->Data()->SetTolerance(aTolU, aTolV);
    aStruct->Data()->SetCellSize(14.0 * aTolU, 14.0 * aTolV);

    std::vector<int> aMap(aNodes.size() + 1, 0);
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
        if (aDWire->EdgesNb() == 16)
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
{
protected:
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

static int runStructDump(const TopoDS_Shape& theShape, const int theFaceIdx)
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
  aParams.Deflection    = aMaxComp * 0.001 * 4.0;
  aParams.Angle         = 20.0 * M_PI / 180.0;
  aParams.InParallel    = false;
  aParams.Relative      = false;
  aParams.MinSize       = (std::max)(IMeshTools_Parameters::RelMinSize() * aParams.Deflection,
                                     Precision::Confusion());
  aParams.AdjustMinSize = false;

  occ::handle<BRepMesh_Context> aCtx = new BRepMesh_Context;
  aCtx->SetShape(theShape);
  aCtx->ChangeParameters()            = aParams;
  aCtx->ChangeParameters().CleanModel = false;
  const bool aBuilt   = aCtx->BuildModel();
  const bool aDiscret = aBuilt && aCtx->DiscretizeEdges();
  const bool aHealed  = aDiscret && aCtx->HealModel();
  const bool aPreproc = aHealed && aCtx->PreProcessModel();
  const occ::handle<IMeshData_Model>& aModel = aCtx->GetModel();
  std::printf("OCCTSTRUCT defl=%.9f angle=%.9f built=%d discret=%d healed=%d preproc=%d faces=%d\n",
              aParams.Deflection,
              aParams.Angle,
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
      }
    }
  }

  occ::handle<IMeshTools_MeshAlgo> aAlgo;
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
    if (!(std::abs(v[0] - theBBox[0]) < 1e-3 && std::abs(v[1] - theBBox[1]) < 1e-3
          && std::abs(v[2] - theBBox[2]) < 1e-3 && std::abs(v[3] - theBBox[3]) < 1e-3
          && std::abs(v[4] - theBBox[4]) < 1e-3 && std::abs(v[5] - theBBox[5]) < 1e-3))
    {
      continue;
    }
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
          std::printf("%s(%.9f,%.9f)", k ? " " : "", aPCurve->GetPoint(k).X(),
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
      double aF = 0, aL = 0;
      Handle(Geom2d_Curve) aPc = BRep_Tool::CurveOnSurface(aEdge, aFace, aF, aL);
      TopLoc_Location aELoc;
      const Handle(Geom_Curve) aC3 = BRep_Tool::Curve(aEdge, aELoc, aF, aL);
      aOut << "  EDGE[" << aEdgeIdx << "] "
           << (aEdge.Orientation() == TopAbs_REVERSED ? "REV" : "FWD")
           << " c3=" << (aC3.IsNull() ? "null" : aC3->DynamicType()->Name())
           << " pc=" << (aPc.IsNull() ? "null" : aPc->DynamicType()->Name())
           << " pcrange=[" << aF << "," << aL << "]";
      if (!aPc.IsNull())
      {
        gp_Pnt2d aP1, aP2;
        aPc->D0(aF, aP1);
        aPc->D0(aL, aP2);
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
          gp_Pnt2d aQ1, aQ2;
          aPc->D0(aF, aQ1);
          aPc->D0(aL, aQ2);
          aOut << " q1=(" << aQ1.X() << "," << aQ1.Y() << ")"
               << " q2=(" << aQ2.X() << "," << aQ2.Y() << ")";
          aOut << " dPar=" << (aL - aF)
               << " aLLength=" << aQ1.Distance(aQ2);
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
    }
    if (std::string(argv[i]) == "--delaunstruct" && i + 1 < argc)
    {
      return runDelaunStruct(argv[i + 1]);
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

  STEPControl_Reader aReader;
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
  const TopoDS_Shape aShape = aReader.OneShape();
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
    // TEMP T-101: `--facestats <defl> <angle>` — per-model-face nodes/triangles/bbox.
    if (std::string(argv[i]) == "--facestats" && i + 2 < argc)
    {
      return runFaceStats(aShape, std::atof(argv[i + 1]), std::atof(argv[i + 2]));
    }
    // TEMP T-103: `--torusdbg <defl> <angle>` — replay OCCT's torus range
    // splitter per torus face of the discrete model (see _dbg/ZZ_TorusRangeSplitter.cxx).
    if (std::string(argv[i]) == "--torusdbg" && i + 2 < argc)
    {
      return runTorusDbg(aShape, std::atof(argv[i + 1]), std::atof(argv[i + 2]));
    }
    // TEMP T-103: `--healerdbg <defl> <angle> <x1> <y1> <z1> <x2> <y2> <z2>` — dump the
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
    // TEMP T-103: `--orderdbg <x1> <y1> <z1> <x2> <y2> <z2>` — dump the wire ordering
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
    // TEMP T-103: `--wiredbg <defl> <angle> <x1> <y1> <z1> <x2> <y2> <z2>` — ground
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
    return runStructDump(aShape, aStructDumpFace);
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
