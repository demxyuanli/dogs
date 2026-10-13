import sys, traceback

STEP = r"D:\source\repos\dogs\data\Torus.step"
OUT = r"D:\source\repos\dogs\_tmp\t024mesh\freecad_torus_out.iges"

import FreeCAD, Part, Import

try:
    Import.open(STEP)
    doc = FreeCAD.ActiveDocument
    objs = [o for o in doc.Objects if hasattr(o, "Shape")]
    for o in objs:
        s = o.Shape
        print("obj %s: faces=%d" % (o.Name, len(s.Faces)))
        for i, f in enumerate(s.Faces):
            print("  face %d type=%s area=%.6f nwires=%d" % (
                i, f.Surface.TypeId, f.Area, len(f.Wires)))
    Part.export(objs, OUT)
    print("wrote", OUT)
except Exception:
    traceback.print_exc()

sys.stdout.flush()
