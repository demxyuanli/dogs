// Fallback report printed when no `argv[2]` mode matched: `BRepCheck` verdict,
// `BRepGProp` volume/centre, and one line per face (type, orientation, outward sign,
// UV window, surface frame under `--frames`, wire edges).
#include "occt_probe_common.h"

int runDefaultReport(const ProbeCtx& theCtx)
{
  const TopoDS_Shape& aShape = theCtx.shape;
  const bool          frames = theCtx.frames;
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
