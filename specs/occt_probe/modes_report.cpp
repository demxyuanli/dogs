// `argv[2]` modes: `--period`, `--bbox`, `--entities`, `--perface`.
#include "occt_probe_common.h"

std::optional<int> runTopologyReportModes(const ProbeCtx& theCtx)
{
  const int           argc = theCtx.argc;
  char** const        argv = theCtx.argv;
  const TopoDS_Shape& aShape = theCtx.shape;
  STEPControl_Reader& aReader = theCtx.reader;
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
  return std::nullopt;
}
