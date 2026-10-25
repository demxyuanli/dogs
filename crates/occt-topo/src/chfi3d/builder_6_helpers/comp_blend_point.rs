//! `CompBlendPoint` — the `Blend_Point` of a tangency on a vertex.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_6.cxx:304-339`
//! (OCCT 8.0.0, a file-static helper of `ChFi3d_Builder`).
//!
//! The point is the vertex itself on both surfaces, its 1D parameter on the
//! guide is `W`, and the two 2D parameters are read on the p-curves of the edge
//! on the two faces:
//!
//! ```text
//! P3d   = Pnt(V);
//! param = Parameter(V, E, F1);  P1 = pc1->Value(param);
//! param = Parameter(V, E, F2);  P2 = pc2->Value(param);
//! BP.SetValue(P3d, P3d, W, P1.X(), P1.Y(), P2.X(), P2.Y());
//! ```
//!
//! Divergences that are not heuristics:
//! - `CurveOnSurface(E, F, f, l)` is
//!   [`boptools_2d::curve_on_surface_range`], which keeps the
//!   `BRep_Tool.cxx:367-372` fallback to `CurveOnPlane` when no
//!   `BRep_CurveOnSurface` representation is stored,
//! - `BRep_Tool::Parameter(V, E, F)` is
//!   [`BRepTool::parameter_of_vertex_on_face`]; OCCT raises
//!   `Standard_NoSuchObject` when the vertex is not a child of the edge, which
//!   cannot happen on the `ChFiDS_Spine` edges this helper is called on
//!   (`_6.cxx:1783-1786`, `:1801-1804`: `Spine->FirstVertex()` /
//!   `Edges(1)` and `Spine->LastVertex()` / `Edges(Spine->NbEdges())`),
//! - the `pc.IsNull()` guard becomes the `None` of `curve_on_surface_range`,
//!   and a null p-curve answers `None` exactly like OCCT's `false`,
//! - `Blend_Point::SetValue(P1, P2, Param, U1, V1, U2, V2)`
//!   (`Blend_Point.cxx:111-131`) is [`BlendPoint::set_on_two_surfaces`].

use crate::blend::point::BlendPoint;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, Vertex};

/// `CompBlendPoint(V, E, W, F1, F2, BP)` (`cxx:310-339`): the blending point of
/// a tangency at vertex `v` of edge `e`, with guide parameter `w`.
///
/// `None` is OCCT's `false`: at least one of the two faces carries no p-curve
/// for `e`.
pub fn comp_blend_point(
    v: &Vertex,
    e: &Edge,
    w: f64,
    f1: &Face,
    f2: &Face,
) -> Option<BlendPoint> {
    // `P3d = BRep_Tool::Pnt(V)` (`cxx:322`).
    let p3d = BRepTool::vertex_point(v);

    // `param = BRep_Tool::Parameter(V, E, F1); pc = CurveOnSurface(E, F1, f, l)`
    // (`cxx:323-324`).
    let param = BRepTool::parameter_of_vertex_on_face(&v.0, e, f1);
    let (pc, _, _) = crate::boptools_2d::curve_on_surface_range(e, f1)?;
    // `P1 = pc->Value(param)` (`cxx:329`).
    let p1 = pc.value(param);

    // `param = BRep_Tool::Parameter(V, E, F2); pc = CurveOnSurface(E, F2, f, l)`
    // (`cxx:330-331`).
    let param = BRepTool::parameter_of_vertex_on_face(&v.0, e, f2);
    let (pc, _, _) = crate::boptools_2d::curve_on_surface_range(e, f2)?;
    // `P2 = pc->Value(param)` (`cxx:336`).
    let p2 = pc.value(param);

    // `BP.SetValue(P3d, P3d, W, P1.X(), P1.Y(), P2.X(), P2.Y())`
    // (`cxx:337`).
    let mut bp = BlendPoint::new();
    bp.set_on_two_surfaces(p3d, p3d, w, p1.x(), p1.y(), p2.x(), p2.y());
    Some(bp)
}
