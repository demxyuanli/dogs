import sys, traceback

P = r"D:\source\repos\dogs\_igesdump\oracle\Torus.iges"

import FreeCAD, Part, Import

try:
    Import.open(P)
    doc = FreeCAD.ActiveDocument
    objs = [o for o in doc.Objects if hasattr(o, "Shape")]
    for o in objs:
        s = o.Shape
        print("obj %s: faces=%d shells=%d solids=%d" % (o.Name, len(s.Faces), len(s.Shells), len(s.Solids)))
        for i, f in enumerate(s.Faces):
            print("  face %d type=%s area=%.6f" % (i, f.Surface.TypeId, f.Area))
            for j, w in enumerate(f.Wires):
                print("    wire %d edges=%d area=%.6f" % (j, len(w.Edges), w.Area))
            # parametric bounds of the surface as FreeCAD sees it, plus the
            # face's own (wire-driven) parameter range
            try:
                print("    surf.ParameterRange=%s" % (f.Surface.ParameterRange,))
            except Exception as e:
                print("    surf.ParameterRange error %s" % e)
            try:
                print("    face.ParameterRange=%s" % (f.ParameterRange,))
            except Exception as e:
                print("    face.ParameterRange error %s" % e)
            try:
                print("    uvs=%s" % (Part.Face(f).Surface.parameter(f.CenterOfMass),))
            except Exception as e:
                print("    parameter error %s" % e)
except Exception:
    traceback.print_exc()

sys.stdout.flush()
