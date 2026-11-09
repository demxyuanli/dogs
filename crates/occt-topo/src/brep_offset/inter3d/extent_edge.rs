//! File-local `ExtentEdge` of `BRepOffset_Inter3d.cxx:62-88` (a different
//! function from `BRepOffset_Tool::ExtentEdge`): a copy of `E` whose 3D range
//! is stretched by 100 lengths on both sides, with two fresh vertices built at
//! the stretched ends.
//!
//! Callers in the reference are `ContextIntByInt` (`:1137`) and
//! `ContextIntByArc` (`:1326`).

use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face};
use crate::tgeometry::GeometryRegistry;

/// `ExtentEdge(F, E, NE)` (`cxx:62-88`). The face argument is unused in the
/// reference (`/*F*/`) and stays unused here.
pub(crate) fn extent_edge(e: &Edge, _f: &Face) -> Edge {
    // `NE = TopoDS::Edge(E.EmptyCopied()); NE.Orientation(FORWARD)`.
    let mut ne = e.0.empty_copied();
    ne.set_orientation(crate::abs::Orientation::Forward);

    // `BRep_Tool::Range(E, f, l); length = l - f; f -= 100*length;
    // l += 100*length; B.Range(NE, f, l)`.
    let (f, l) = BRepTool::edge_parameters(e);
    let length = l - f;
    let new_f = f - 100.0 * length;
    let new_l = l + 100.0 * length;
    GeometryRegistry::global().set_ranges_all(&ne, new_f, new_l);

    // `BRepAdaptor_Curve CE(E)` is the edge's own 3D curve (no surface), so
    // `CE.Value(f)` extrapolates it at the stretched parameters.
    if let Some(curve) = BRepTool::edge_curve_world(e) {
        let v1 = TopoBuilder::new().make_vertex(curve.d0(new_f), occt_core::precision::CONFUSION);
        let v2 = TopoBuilder::new().make_vertex(curve.d0(new_l), occt_core::precision::CONFUSION);
        // `B.Add(NE, V1.Oriented(FORWARD)); B.Add(NE, V2.Oriented(REVERSED))`.
        TopoBuilder::new().add(&mut ne, &v1.0.oriented(crate::abs::Orientation::Forward));
        TopoBuilder::new().add(&mut ne, &v2.0.oriented(crate::abs::Orientation::Reversed));
    }

    // `NE.Orientation(E.Orientation())`.
    ne.set_orientation(e.0.orientation());
    Edge(ne)
}
