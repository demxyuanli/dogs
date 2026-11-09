//! File-local statics of `LocOpe_Gluer.cxx`.

use crate::shape::TopoShape;
use crate::topo_tools::is_same;

/// `Contains(const NCollection_List<TopoDS_Shape>& L, const TopoDS_Shape& S)`
/// (`cxx:456-468`): `Perform` uses it to decide whether a face of the glued
/// shape came back from `LocOpe_GluedShape::OrientedFaces` (`cxx:222-225`).
pub(super) fn contains(list: &[TopoShape], s: &TopoShape) -> bool {
    list.iter().any(|v| is_same(v, s))
}

/// `AddEdges()` (`cxx:471-555`) is **not ported**: it has no output.
///
/// The reference maps the faces of `mySn` (`cxx:483`), then for every face of
/// `myRes` that is not one of them builds `MapV` (its vertices,
/// `cxx:488-489`) and accumulates `MapE` (its edges, `cxx:490`, never cleared
/// across faces), and for every edge of `mySn` outside `MapE` with a vertex in
/// `MapV` runs `BRepExtrema_ExtPF` against that face (`cxx:492-548`). The
/// resulting `flag` only guards the empty block `if (flag == 1) { }`
/// (`cxx:549-551`); `myRes`, `myEdges` and `myTgtEdges` are never touched, and
/// no other member is read afterwards.
///
/// The call site (`cxx:235`) is kept so `Perform`'s statement order still
/// matches the reference.
pub(super) fn add_edges() {}
