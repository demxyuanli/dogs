// `argv[2]` modes: `--mesh`, `--dumpmesh`, `--wireinv`, `--typefaces`,
// `--adv` / `--advbox`.
#include "occt_probe_common.h"

std::optional<int> runMeshOracleModes(const ProbeCtx& theCtx)
{
  const int           argc = theCtx.argc;
  char** const        argv = theCtx.argv;
  const TopoDS_Shape& aShape = theCtx.shape;
  STEPControl_Reader& aReader = theCtx.reader;

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
  return std::nullopt;
}
