// One-off OCCT ground-truth probe: per-face BRepGProp area plus the face's
// pcurve U branch, used to cross-check brep_gprop_full against BRepGProp.
//
// usage: zz_face_area_probe <step-file> [minArea maxArea]
#include <BRepGProp.hxx>
#include <BRepGProp_Face.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <GProp_GProps.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom_Surface.hxx>
#include <STEPControl_Reader.hxx>
#include <TopExp_Explorer.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Shape.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <cstdio>
#include <cstdlib>

int main(int argc, char** argv)
{
  if (argc < 2)
  {
    std::printf("usage: %s <step-file> [minArea maxArea]\n", argv[0]);
    return 1;
  }
  const double aMin = (argc > 3) ? std::atof(argv[2]) : -1e300;
  const double aMax = (argc > 3) ? std::atof(argv[3]) : 1e300;

  STEPControl_Reader aReader;
  if (aReader.ReadFile(argv[1]) != IFSelect_RetDone)
  {
    std::printf("read failed\n");
    return 1;
  }
  aReader.TransferRoots();
  TopoDS_Shape aShape = aReader.OneShape();

  GProp_GProps aTotal;
  BRepGProp::SurfaceProperties(aShape, aTotal);
  std::printf("TOTAL area = %.9f\n", aTotal.Mass());

  TopTools_IndexedMapOfShape aFaces;
  TopExp::MapShapes(aShape, TopAbs_FACE, aFaces);
  std::printf("faces = %d\n", aFaces.Extent());

  for (int i = 1; i <= aFaces.Extent(); ++i)
  {
    const TopoDS_Face& aFace = TopoDS::Face(aFaces.FindKey(i));
    GProp_GProps       aProps;
    BRepGProp::SurfaceProperties(aFace, aProps);
    const double anArea = aProps.Mass();
    if (anArea < aMin || anArea > aMax)
      continue;

    BRepAdaptor_Surface anAdaptor(aFace);
    double               anU1, anU2, aV1, aV2;
    BRepTools::UVBounds(aFace, anU1, anU2, aV1, aV2);
    std::printf("FACE %d area=%.9f type=%d uv=(%.12f,%.12f,%.12f,%.12f) nchild=%d\n",
                i,
                anArea,
                (int)anAdaptor.GetType(),
                anU1,
                anU2,
                aV1,
                aV2,
                aFace.NbChildren());

    for (TopExp_Explorer aW(aFace, TopAbs_WIRE); aW.More(); aW.Next())
    {
      for (BRepTools_WireExplorer aWE(TopoDS::Wire(aW.Current()), aFace); aWE.More(); aWE.Next())
      {
        const TopoDS_Edge&                 anEdge = aWE.Current();
        double                             a, b;
        TopLoc_Location                    aLoc;
        const occ::handle<Geom2d_Curve> aC2d = BRep_Tool::CurveOnSurface(anEdge, aFace, a, b);
        if (aC2d.IsNull())
        {
          std::printf("  EDGE orient=%d nopcurve\n", (int)anEdge.Orientation());
          continue;
        }
        const gp_Pnt2d aP0 = aC2d->Value(a);
        const gp_Pnt2d aP1 = aC2d->Value(b);
        std::printf(
          "  EDGE orient=%d t=(%.15g,%.15g) uv0=(%.12f,%.12f) uv1=(%.12f,%.12f) mid=(%.12f,%.12f)\n",
          (int)anEdge.Orientation(),
          a,
          b,
          aP0.X(),
          aP0.Y(),
          aP1.X(),
          aP1.Y(),
          aC2d->Value(0.5 * (a + b)).X(),
          aC2d->Value(0.5 * (a + b)).Y());
      }
    }
  }
  return 0;
}
