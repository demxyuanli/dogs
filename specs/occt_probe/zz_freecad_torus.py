import sys, traceback

P = r"D:\source\repos\dogs\_igesdump\oracle\Torus.iges"
OUT = r"D:\source\repos\dogs\_tmp\t024mesh\freecad_torus.obj"
DEFL = 0.096
ANG = 0.3490658503988659

import FreeCAD, Part, MeshPart

def report(tag, shape):
    m = MeshPart.meshFromShape(
        Shape=shape,
        LinearDeflection=DEFL,
        AngularDeflection=ANG,
        Relative=False,
    )
    bb = m.BoundBox
    print("%s faces=%d points=%d facets=%d bbox=(%.6f,%.6f,%.6f)-(%.6f,%.6f,%.6f)" % (
        tag, len(shape.Faces), m.CountPoints, m.CountFacets,
        bb.XMin, bb.YMin, bb.ZMin, bb.XMax, bb.YMax, bb.ZMax))
    return m

try:
    import Import
    Import.open(P)
    doc = FreeCAD.ActiveDocument
    objs = [o for o in doc.Objects if hasattr(o, "Shape")]
    print("objects:", [(o.Name, o.TypeId) for o in objs])
    for o in objs:
        s = o.Shape
        print("  obj %s faces=%d solids=%d shells=%d" % (
            o.Name, len(s.Faces), len(s.Solids), len(s.Shells)))
        m = report("freecad:" + o.Name, s)
        m.write(OUT)
        print("  wrote", OUT)
except Exception:
    traceback.print_exc()

sys.stdout.flush()
