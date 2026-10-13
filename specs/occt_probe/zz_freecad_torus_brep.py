import sys, traceback

P = r"D:\source\repos\dogs\_igesdump\oracle\Torus.iges"
BREP = r"D:\source\repos\dogs\_tmp\t024mesh\freecad_torus.brep"

import FreeCAD, Part, Import

try:
    Import.open(P)
    doc = FreeCAD.ActiveDocument
    objs = [o for o in doc.Objects if hasattr(o, "Shape")]
    for o in objs:
        s = o.Shape
        print("obj %s: faces=%d wires_total=%d edges_total=%d shells=%d solids=%d" % (
            o.Name, len(s.Faces), len(s.Wires), len(s.Edges), len(s.Shells), len(s.Solids)))
        for i, f in enumerate(s.Faces):
            ws = list(f.Wires)
            print("  face %d type=%s area=%.6f nwires=%d" % (i, f.Surface.TypeId, f.Area, len(ws)))
            for j, w in enumerate(ws):
                print("    wire %d edges=%d area=%.6f closed=%s" % (
                    j, len(w.Edges), w.Area, w.isClosed()))
        s.exportBrep(BREP)
        print("wrote", BREP)
except Exception:
    traceback.print_exc()

sys.stdout.flush()
