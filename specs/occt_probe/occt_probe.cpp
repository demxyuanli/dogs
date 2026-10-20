// OCCT ground-truth probe (T-69/T-80/T-87/T-05 oracle) - entry point and `run*` loop.
// The `argv[2]` command groups live in the `modes_*.cpp` files; the structural dumps
// they share live in the matching `*_probe.cpp` file.
#include "occt_probe_common.h"

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
    if (std::string(argv[i]) == "--ssx")
    {
      return runExtSS();
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
    if (std::string(argv[i]) == "--nofix" || std::string(argv[i]) == "--degchkraw")
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
    // T-126: `--offsetsurf <dist>` dumps `BRepOffset::Surface` per face.
    if (std::string(argv[i]) == "--offsetsurf" && i + 1 < argc)
    {
      return runOffsetSurf(aShape, std::atof(argv[i + 1]));
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

  const ProbeCtx aCtx = {argc, argv, aShape, aReader, frames};
  if (std::optional<int> aRc = runMeshOracleModes(aCtx)) { return *aRc; }
  if (std::optional<int> aRc = runPcurveModes(aCtx)) { return *aRc; }
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
  if (std::optional<int> aRc = runUvModes(aCtx)) { return *aRc; }
  if (std::optional<int> aRc = runFaceIdModes(aCtx)) { return *aRc; }
  if (std::optional<int> aRc = runEdgeFaceModes(aCtx)) { return *aRc; }
  if (std::optional<int> aRc = runFixFaceModes(aCtx)) { return *aRc; }
  if (std::optional<int> aRc = runTopologyReportModes(aCtx)) { return *aRc; }
  if (std::optional<int> aRc = runBooleanModes(aCtx)) { return *aRc; }
  if (std::optional<int> aRc = runClassifyModes(aCtx)) { return *aRc; }
  return runDefaultReport(aCtx);
}
