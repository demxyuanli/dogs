// OCCT ground-truth generator for the BRep modeling parity gate.
//
// Constructs the modeling-gate models with the OCCT 8.0.0 construction APIs
// (BRepPrimAPI_* / BRepBuilderAPI_* / BRepPrimAPI_MakePrism) and freezes, for
// each model:
//   - the topology census (TopExp::MapShapes per shape type),
//   - BRepCheck_Analyzer validity,
//   - BRepGProp volume / surface area / centroid,
//   - a BRepMesh_IncrementalMesh(0.1) tessellation written as OBJ.
//
// The frozen records (data/modeling/<model>.txt + <model>.obj) are the oracle
// of crates/occt-topo/tests/modeling_gates.rs. They are produced once and
// committed; the gate never rebuilds OCCT.
//
// Usage: modeling_probe <output-dir>
//
// Model parameters match the existing data/*.step samples so the modeled solid
// can also be compared against the STEP-read solid (Cube 10^3, Cylinder r2 h10,
// Cone r4 h10, Sphere r5, Torus R10 r2).
//
// G0 reading (2026-10-16, first run of tests/modeling_gates.rs against these
// records): prim-box, prim-cylinder and prim-torus pass all three layers.
// Open port diffs, with the OCCT control flow each one is missing:
//   prim-cone   topology: E ours 2 vs occt 3. `BRepPrim_OneAxis::LateralWire`
//               (BRepPrim_OneAxis.cxx:660-683) is [TopEdge, EndEdge, BottomEdge,
//               StartEdge]; with a finite VMax whose meridian sits on the axis
//               `TopEdge` is a degenerate edge at the apex
//               (BRepPrim_OneAxis.cxx:1205-1213, MakeDegeneratedEdge(ETOP)) and
//               `StartEdge` aliases `EndEdge` (cxx:976-979, `!HasSides()`), so
//               OCCT has 3 edges = base circle + seam + apex degenerate. The
//               port has no apex degenerate edge.
//               gprop: volume 1091.476 vs 167.552, area 576.404 vs 185.609,
//               centroid (-0.394, 0, 5.240) vs (0, 0, 2.5). The port builds the
//               lateral surface with `semi_angle = +atan(radius/height)`
//               (primitives.rs make_cone), whereas OCCT uses
//               `atan((R2 - R1) / H)` (BRepPrim_Cone.cxx:170-173), negative for
//               this cone, and parameterizes V along the generatrix
//               (`ElSLib::ConeD0`, ElSLib.cxx:572-590: R = Radius + V*sin(SAngle),
//               height = V*cos(SAngle)), with VMax = sqrt(H^2 + (R2-R1)^2) = the
//               slant length. With the port's sign the apex is off the surface,
//               so the lateral face carries a wire that does not lie on it.
//               mesh: our OBJ has 0 vertices (occt 103) - same malformed face.
//   prim-sphere topology: E ours 4 vs occt 3, and BRepCheck ours false vs occt
//               true. The port builds two distinct meridian edges (start, end),
//               where OCCT aliases them (`StartEdge`/`EndEdge`,
//               BRepPrim_OneAxis.cxx:976-979, `!HasSides()` is true for a full
//               revolution); OCCT's sphere wire is seam + one degenerate edge
//               per pole (cxx:1212 / cxx:1264).
//   prim-prism  gprop volume ours 16 vs occt 24. Area (52), centroid and the
//               deflection-0.1 mesh all match exactly, and the census matches,
//               so the shape is right and the volume integral loses one face
//               contribution (the boundary sum is 48 instead of 72 = 3*V).
// These are port gaps, not oracle noise: each occt value equals the closed form
// (cone 160*pi/3, prism 24, sphere 500*pi/3).

#include <algorithm>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <utility>
#include <vector>

#include <Bnd_Box.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBndLib.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepPrimAPI_MakeTorus.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <Poly_Triangulation.hxx>
#include <TopAbs.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Ax2.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

namespace {

struct Model {
    std::string name;
    TopoDS_Shape shape;
};

int counts_of(const TopoDS_Shape& s, TopAbs_ShapeEnum t) {
    TopTools_IndexedMapOfShape m;
    TopExp::MapShapes(s, t, m);
    return m.Extent();
}

std::string shape_type_name(const TopoDS_Shape& s) {
    switch (s.ShapeType()) {
        case TopAbs_VERTEX: return "Vertex";
        case TopAbs_EDGE: return "Edge";
        case TopAbs_WIRE: return "Wire";
        case TopAbs_FACE: return "Face";
        case TopAbs_SHELL: return "Shell";
        case TopAbs_SOLID: return "Solid";
        case TopAbs_COMPSOLID: return "CompSolid";
        case TopAbs_COMPOUND: return "Compound";
        default: return "Shape";
    }
}

// Appends the face's triangulation vertices (in world coordinates) to `v` and
// the 1-based triangles to `f`. No vn/vt records: the gate counts `v `/`f `
// lines, so extra record types would corrupt the census.
void collect_face(const TopoDS_Face& face, std::vector<gp_Pnt>& v,
                  std::vector<std::array<int, 3> >& f) {
    TopLoc_Location loc;
    Handle(Poly_Triangulation) tri = BRep_Tool::Triangulation(face, loc);
    if (tri.IsNull()) return;
    const gp_Trsf& tr = loc.Transformation();
    const int base = static_cast<int>(v.size());
    for (int i = 1; i <= tri->NbNodes(); ++i) {
        gp_Pnt p = tri->Node(i);
        if (!loc.IsIdentity()) p.Transform(tr);
        v.push_back(p);
    }
    for (int i = 1; i <= tri->NbTriangles(); ++i) {
        Standard_Integer n1, n2, n3;
        tri->Triangle(i).Get(n1, n2, n3);
        if (face.Orientation() == TopAbs_REVERSED) std::swap(n2, n3);
        std::array<int, 3> t = {base + n1, base + n2, base + n3};
        f.push_back(t);
    }
}

void write_obj(const std::string& path, const std::string& name,
               const TopoDS_Shape& shape) {
    std::vector<gp_Pnt> v;
    std::vector<std::array<int, 3> > f;
    for (TopExp_Explorer ex(shape, TopAbs_FACE); ex.More(); ex.Next()) {
        collect_face(TopoDS::Face(ex.Current()), v, f);
    }
    FILE* fp = fopen(path.c_str(), "w");
    if (fp == NULL) {
        fprintf(stderr, "cannot write %s\n", path.c_str());
        exit(2);
    }
    fprintf(fp, "# Exported by Open CASCADE Technology [dev.opencascade.org]\n");
    fprintf(fp, "#  Vertices: %d\n", static_cast<int>(v.size()));
    fprintf(fp, "#     Faces: %d\n", static_cast<int>(f.size()));
    fprintf(fp, "g %s\n", name.c_str());
    for (size_t i = 0; i < v.size(); ++i) {
        fprintf(fp, "v %.6f %.6f %.6f\n", v[i].X(), v[i].Y(), v[i].Z());
    }
    for (size_t i = 0; i < f.size(); ++i) {
        fprintf(fp, "f %d %d %d\n", f[i][0], f[i][1], f[i][2]);
    }
    fclose(fp);
}

void emit(const std::string& dir, const Model& m) {
    const TopoDS_Shape& s = m.shape;
    const std::string txt = dir + "/" + m.name + ".txt";
    const std::string obj = dir + "/" + m.name + ".obj";

    BRepCheck_Analyzer ana(s);
    GProp_GProps vp, sp;
    BRepGProp::VolumeProperties(s, vp);
    BRepGProp::SurfaceProperties(s, sp);
    const gp_Pnt c = vp.CentreOfMass();

    Bnd_Box box;
    BRepBndLib::Add(s, box);
    Standard_Real xmin = 0, ymin = 0, zmin = 0, xmax = 0, ymax = 0, zmax = 0;
    box.Get(xmin, ymin, zmin, xmax, ymax, zmax);

    BRepMesh_IncrementalMesh mesh(s, 0.1);
    write_obj(obj, m.name, s);

    FILE* fp = fopen(txt.c_str(), "w");
    if (fp == NULL) {
        fprintf(stderr, "cannot write %s\n", txt.c_str());
        exit(2);
    }
    fprintf(fp, "MODEL %s\n", m.name.c_str());
    fprintf(fp, "TYPE %s\n", shape_type_name(s).c_str());
    fprintf(fp, "COUNTS V=%d E=%d W=%d FA=%d SH=%d SO=%d CO=%d CP=%d\n",
            counts_of(s, TopAbs_VERTEX), counts_of(s, TopAbs_EDGE),
            counts_of(s, TopAbs_WIRE), counts_of(s, TopAbs_FACE),
            counts_of(s, TopAbs_SHELL), counts_of(s, TopAbs_SOLID),
            counts_of(s, TopAbs_COMPSOLID), counts_of(s, TopAbs_COMPOUND));
    fprintf(fp, "VALID %d\n", ana.IsValid() ? 1 : 0);
    fprintf(fp, "VOLUME %.17g\n", vp.Mass());
    fprintf(fp, "AREA %.17g\n", sp.Mass());
    fprintf(fp, "CENTROID %.17g %.17g %.17g\n", c.X(), c.Y(), c.Z());
    fprintf(fp, "BBOX %.17g %.17g %.17g %.17g %.17g %.17g\n", xmin, ymin, zmin,
            xmax, ymax, zmax);
    fprintf(fp, "END\n");
    fclose(fp);
    printf("MODEL %-16s faces=%d valid=%d vol=%.6f area=%.6f\n", m.name.c_str(),
           counts_of(s, TopAbs_FACE), ana.IsValid() ? 1 : 0, vp.Mass(), sp.Mass());
}

std::vector<Model> build_models() {
    std::vector<Model> models;

    // G0 primitives: parameters equal to the data/*.step samples.
    models.push_back(Model{"prim-box", BRepPrimAPI_MakeBox(10.0, 10.0, 10.0).Shape()});
    models.push_back(Model{"prim-cylinder", BRepPrimAPI_MakeCylinder(2.0, 10.0).Shape()});
    models.push_back(Model{"prim-cone", BRepPrimAPI_MakeCone(4.0, 0.0, 10.0).Shape()});
    models.push_back(Model{"prim-sphere", BRepPrimAPI_MakeSphere(5.0).Shape()});
    models.push_back(Model{"prim-torus", BRepPrimAPI_MakeTorus(10.0, 2.0).Shape()});

    // G0 prism: a 2x3 rectangle in XY extruded by +Z 4 (BRepBuilderAPI profile
    // + BRepPrimAPI_MakePrism).
    {
        BRepBuilderAPI_MakePolygon poly;
        poly.Add(gp_Pnt(0.0, 0.0, 0.0));
        poly.Add(gp_Pnt(2.0, 0.0, 0.0));
        poly.Add(gp_Pnt(2.0, 3.0, 0.0));
        poly.Add(gp_Pnt(0.0, 3.0, 0.0));
        poly.Close();
        BRepBuilderAPI_MakeFace mkFace(poly.Wire());
        models.push_back(
            Model{"prim-prism", BRepPrimAPI_MakePrism(mkFace.Face(), gp_Vec(0.0, 0.0, 4.0)).Shape()});
    }

    return models;
}

}  // namespace

int main(int argc, char** argv) {
    const std::string dir = (argc > 1) ? argv[1] : "data/modeling";
    const std::vector<Model> models = build_models();
    for (size_t i = 0; i < models.size(); ++i) emit(dir, models[i]);
    printf("DONE %d models -> %s\n", static_cast<int>(models.size()), dir.c_str());
    return 0;
}
