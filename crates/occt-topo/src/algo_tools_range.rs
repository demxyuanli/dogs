//! Remainder of `BOPTools_AlgoTools_2.cxx` and `IntTools_Tools.cxx` helpers
//! used by `PerformEE` / `PerformEF` / `MakeSplitEdges`.
//!
//! Source:
//! - `BOPTools_AlgoTools_2.cxx` (`UpdateVertex` at 36/57/80, `MakeSectEdge`
//!   at 102, `CopyEdge` at 124, `MakeSplitEdge` at 138, `MakeNewVertex` at
//!   187/224/254, `CorrectRange` at 284/364, `Dimensions` at 467)
//! - `IntTools_Tools.cxx` (`VertexParameters` at 593, `VertexParameter` at
//!   615, `IsOnPave1` at 627, `IsInRange` at 650)
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

use occt_core::gp::GpPnt;
use occt_core::precision::{PCONFUSION, RESOLUTION};
use occt_geom::Curve;

use crate::abs::{Orientation, ShapeType};
use crate::algo_tools::{AlgoTools, D_TOLERANCE};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::inttools_data::{CommonPrt, IntRange};
use crate::iterator::ShapeIterator;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;

/// `BOPTools_AlgoTools::PointOnEdge`.
pub fn point_on_edge(edge: &Edge, t: f64) -> Option<GpPnt> {
    BRepTool::edge_curve(edge).map(|c| c.d0(t))
}

/// `BOPTools_AlgoTools::UpdateVertex` (vertex / vertex).
pub fn update_vertex_vv(vf: &Vertex, a_new: &Vertex) {
    let pvf = BRepTool::vertex_point(vf);
    let pnew = BRepTool::vertex_point(a_new);
    let tol_vf = BRepTool::vertex_tolerance(vf);
    let tol_new = BRepTool::vertex_tolerance(a_new);
    let dist = pvf.distance(&pnew);
    let new_tol = dist + tol_new;
    if new_tol > tol_vf {
        vf.set_tolerance(new_tol + D_TOLERANCE);
    }
}

/// `BOPTools_AlgoTools::UpdateVertex` (edge, parameter, vertex).
pub fn update_vertex_et(a_e: &Edge, a_t: f64, v: &Vertex) {
    let pv = BRepTool::vertex_point(v);
    let tol_v = BRepTool::vertex_tolerance(v);
    let Some(pc) = point_on_edge(a_e, a_t) else {
        return;
    };
    let dist = pv.distance(&pc);
    if dist > tol_v {
        v.set_tolerance(dist + D_TOLERANCE);
    }
}

/// `BOPTools_AlgoTools::UpdateVertex` (curve, parameter, vertex).
pub fn update_vertex_ct(curve: &dyn Curve, a_t: f64, v: &Vertex) {
    let pv = BRepTool::vertex_point(v);
    let tol_v = BRepTool::vertex_tolerance(v);
    let pc = curve.d0(a_t);
    let dist = pv.distance(&pc);
    if dist > tol_v {
        v.set_tolerance(dist + D_TOLERANCE);
    }
}

/// `BOPTools_AlgoTools::MakeNewVertex` (two vertices).
pub fn make_new_vertex_vv(a_v1: &Vertex, v2: &Vertex) -> Result<TopoShape, String> {
    let pnt1 = BRepTool::vertex_point(a_v1);
    let tol1 = BRepTool::vertex_tolerance(a_v1);
    let pnt2 = BRepTool::vertex_point(v2);
    let tol2 = BRepTool::vertex_tolerance(v2);
    let dist = pnt1.distance(&pnt2);
    let max_tol = tol1.max(tol2) + 0.5 * dist;
    let a_new = GpPnt::new(
        0.5 * (pnt1.x() + pnt2.x()),
        0.5 * (pnt1.y() + pnt2.y()),
        0.5 * (pnt1.z() + pnt2.z()),
    );
    AlgoTools::make_new_vertex(&a_new, max_tol)
}

/// `BOPTools_AlgoTools::MakeNewVertex` (two edges and parameters).
pub fn make_new_vertex_ee(
    e1: &Edge,
    parm1: f64,
    e2: &Edge,
    parm2: f64,
) -> Result<TopoShape, String> {
    let pnt1 = point_on_edge(e1, parm1).ok_or("make_new_vertex_ee: no curve on edge 1")?;
    let pnt2 = point_on_edge(e2, parm2).ok_or("make_new_vertex_ee: no curve on edge 2")?;
    let tol1 = BRepTool::edge_tolerance(e1);
    let tol2 = BRepTool::edge_tolerance(e2);
    let dist = pnt1.distance(&pnt2);
    let max_tol = tol1.max(tol2) + 0.5 * dist;
    let a_new = GpPnt::new(
        0.5 * (pnt1.x() + pnt2.x()),
        0.5 * (pnt1.y() + pnt2.y()),
        0.5 * (pnt1.z() + pnt2.z()),
    );
    AlgoTools::make_new_vertex(&a_new, max_tol)
}

/// `BOPTools_AlgoTools::MakeNewVertex` (edge, parameter, face).
pub fn make_new_vertex_ef(e1: &Edge, parm1: f64, f1: &Face) -> Result<TopoShape, String> {
    let pnt = point_on_edge(e1, parm1).ok_or("make_new_vertex_ef: no curve on edge")?;
    let tol1 = BRepTool::edge_tolerance(e1);
    let tol2 = BRepTool::face_tolerance(f1);
    let max_tol = tol1 + tol2 + D_TOLERANCE;
    AlgoTools::make_new_vertex(&pnt, max_tol)
}

/// `BOPTools_AlgoTools::MakeVertex` — fuse a chain of vertices into one.
pub fn make_vertex_from_list(lvsd: &[TopoShape]) -> Result<TopoShape, String> {
    if lvsd.is_empty() {
        return Err("make_vertex_from_list: empty chain".into());
    }
    let mut cx = 0.0;
    let mut cy = 0.0;
    let mut cz = 0.0;
    let mut tol: f64 = 0.0;
    for s in lvsd {
        let v = Vertex(s.clone());
        let p = BRepTool::vertex_point(&v);
        cx += p.x();
        cy += p.y();
        cz += p.z();
        tol = tol.max(BRepTool::vertex_tolerance(&v));
    }
    let n = lvsd.len() as f64;
    let p = GpPnt::new(cx / n, cy / n, cz / n);
    for s in lvsd {
        let v = Vertex(s.clone());
        let d = BRepTool::vertex_point(&v).distance(&p);
        tol = tol.max(d + BRepTool::vertex_tolerance(&v));
    }
    AlgoTools::make_new_vertex(&p, tol + D_TOLERANCE)
}

/// `aCT == GeomAbs_BezierCurve || aCT == GeomAbs_BSplineCurve ||
/// aCT == GeomAbs_OffsetCurve || aCT == GeomAbs_OtherCurve`
/// (`BOPTools_AlgoTools_2.cxx:318-319`, `:392-393`).
///
/// This is `GeomAdaptor_Curve::load`'s type tag (`GeomAdaptor_Curve.cxx:250-311`):
/// a `Geom_TrimmedCurve` recurses to its basis, Bezier/B-spline are reported by
/// [`Curve::nurbs_degree`], a `Geom_OffsetCurve` by [`Curve::offset_curve`], and
/// everything that is none of the seven concrete classes - the `else` arm
/// `myTypeCurve = GeomAbs_OtherCurve` - is `OtherCurve`.
///
/// The previous body tested `!curve_is_line` (a 3-point sampled collinearity
/// test), which also swept the circle, ellipse, parabola and hyperbola classes
/// into the derivative path; OCCT sends those through `aBC.Resolution` instead.
fn curve_is_freeform(c: &dyn Curve) -> bool {
    c.nurbs_degree().is_some()
        || c.offset_curve().is_some()
        || !(c.is_line()
            || c.gp_circ().is_some()
            || c.gp_ellipse().is_some()
            || c.gp_parabola().is_some()
            || c.gp_hyperbola().is_some())
}

/// `BOPTools_AlgoTools::CorrectRange` (edge / edge).
pub fn correct_range_ee(e1: &Edge, e2: &Edge, sr: IntRange) -> IntRange {
    let mut new_sr = sr;
    let Some(c) = BRepTool::edge_curve(e1) else {
        return new_sr;
    };
    if c.is_line() {
        return new_sr;
    }
    let d_t = PCONFUSION;
    let tf = sr.first;
    let tl = sr.last;
    let tol_e1 = BRepTool::edge_tolerance(e1);
    let tol_e2 = BRepTool::edge_tolerance(e2);
    let freeform = curve_is_freeform(c.as_ref());
    for i in 0..2 {
        let mut res = 2.0 * (tol_e1 + tol_e2);
        if freeform {
            let t = if i == 0 { tf } else { tl };
            let (_, der) = c.d1(t);
            let mgn = der.magnitude();
            if mgn > 1.0e-12 {
                res /= mgn;
            } else {
                res = c.resolution(res);
            }
        } else {
            res = c.resolution(res);
        }
        if i == 0 {
            new_sr.first = tf + res;
        } else {
            new_sr.last = tl - res;
        }
        if (new_sr.last - new_sr.first) < d_t {
            new_sr = sr;
        }
    }
    new_sr
}

/// `BOPTools_AlgoTools::CorrectRange` (edge / face).
pub fn correct_range_ef(a_e: &Edge, a_f: &Face, sr: IntRange) -> IntRange {
    let mut new_sr = sr;
    let Some(c) = BRepTool::edge_curve(a_e) else {
        return new_sr;
    };
    let d_t = PCONFUSION;
    let tf = sr.first;
    let tl = sr.last;
    let tol_f = BRepTool::face_tolerance(a_f);
    let freeform = curve_is_freeform(c.as_ref());
    for i in 0..2 {
        let mut res = tol_f;
        if freeform {
            let t = if i == 0 { tf } else { tl };
            let (_, der) = c.d1(t);
            let mgn = der.magnitude();
            if mgn > 1.0e-12 {
                res /= mgn;
            } else {
                res = c.resolution(res);
            }
        } else {
            res = c.resolution(res);
        }
        if i == 0 {
            new_sr.first = tf + res;
        } else {
            new_sr.last = tl - res;
        }
        if (new_sr.last - new_sr.first) < d_t {
            new_sr = sr;
        }
    }
    new_sr
}

/// `BOPTools_AlgoTools::CopyEdge`.
pub fn copy_edge(edge: &Edge) -> Result<Edge, String> {
    let (f, l) = BRepTool::edge_parameters(edge);
    let verts: Vec<TopoShape> = ShapeIterator::of_shape(&edge.0)
        .filter(|s| s.shape_type() == ShapeType::Vertex)
        .collect();
    let v1 = verts.first();
    let v2 = verts.get(1).or(verts.first());
    let mut sp = AlgoTools::make_split_edge(edge, v1, f, v2, l)?;
    sp.0.set_orientation(edge.0.orientation());
    Ok(sp)
}

/// Elementary-shape dimension (`BOPTools_AlgoTools_2.cxx` anonymous
/// `dimension`).
pub fn dimension(the_s: &TopoShape) -> i32 {
    match the_s.shape_type() {
        ShapeType::Vertex => 0,
        ShapeType::Edge | ShapeType::Wire => 1,
        ShapeType::Face | ShapeType::Shell => 2,
        ShapeType::Solid | ShapeType::CompSolid => 3,
        _ => -1,
    }
}

/// `BOPTools_AlgoTools::Dimensions`.
pub fn dimensions(the_s: &TopoShape) -> (i32, i32) {
    let d = dimension(the_s);
    if d >= 0 {
        return (d, d);
    }
    let mut d_min = 4;
    let mut d_max = -1;
    for sub in ShapeIterator::of_shape(the_s) {
        let di = dimension(&sub);
        if di < 0 {
            continue;
        }
        if di < d_min {
            d_min = di;
        }
        if di > d_max {
            d_max = di;
        }
    }
    if d_max < 0 {
        (-1, -1)
    } else {
        (d_min, d_max)
    }
}

/// `IntTools_Tools::VertexParameter`.
pub fn algo_tools_vertex_parameter(c_part: &CommonPrt, param: Option<f64>) -> f64 {
    let r = c_part.range;
    let mut a_t = 0.5 * (r.first + r.last);
    if let Some(vp) = param {
        if vp >= r.first && vp <= r.last {
            a_t = vp;
        }
    }
    a_t
}

/// `IntTools_Tools::VertexParameters`.
pub fn algo_tools_vertex_parameters(c_part: &CommonPrt, t1: Option<f64>, t2: Option<f64>) -> (f64, f64) {
    let a_t1 = algo_tools_vertex_parameter(c_part, t1);
    let r2 = c_part.range;
    let mut a_t2 = 0.5 * (r2.first + r2.last);
    if let Some(vp) = t2 {
        if vp >= r2.first && vp <= r2.last {
            a_t2 = vp;
        }
    }
    (a_t1, a_t2)
}

/// `IntTools_Tools::IsOnPave1`.
pub fn algo_tools_is_on_pave1(tr: f64, cp_first: f64, cp_last: f64, tolerance: f64) -> bool {
    if tr >= cp_first && tr <= cp_last {
        return true;
    }
    let d_t1 = (tr - cp_first).abs();
    let d_t2 = (tr - cp_last).abs();
    d_t1 <= tolerance || d_t2 <= tolerance
}

/// `IntTools_Tools::IsInRange`.
pub fn algo_tools_is_in_range(r_ref: IntRange, r: IntRange, tolerance: f64) -> bool {
    let a_t1 = r.first;
    let a_t2 = r.last;
    let t_ref1 = r_ref.first - tolerance;
    let t_ref2 = r_ref.last + tolerance;
    (a_t1 >= t_ref1 && a_t1 <= t_ref2) || (a_t2 >= t_ref1 && a_t2 <= t_ref2)
}

/// `BOPTools_AlgoTools::MakeSectEdge` — wrap [`AlgoTools::make_edge`].
pub fn make_sect_edge(
    curve: std::sync::Arc<dyn Curve>,
    v1: Option<&TopoShape>,
    t1: f64,
    v2: Option<&TopoShape>,
    t2: f64,
    tol: f64,
) -> Result<TopoShape, String> {
    AlgoTools::make_edge(curve, v1, t1, v2, t2, tol)
}

/// Raise a vertex tolerance if `dist` exceeds the current value.
pub fn raise_vertex_tolerance(v: &Vertex, dist: f64) {
    let a_tol = BRepTool::vertex_tolerance(v);
    if dist > a_tol {
        v.set_tolerance(dist + D_TOLERANCE);
    }
}

/// Copy the registered 3-D curve of `edge` (used by `CopyEdge` callers that
/// only need geometry).
pub fn edge_curve_clone(edge: &Edge) -> Option<std::sync::Arc<dyn Curve>> {
    GeometryRegistry::global()
        .edge_geom(&edge.0)
        .map(|g| g.curve.clone())
}

/// `BRep_Builder::MakeWire` + `Add` for a sequence of edges.
pub fn make_wire_of(edges: &[Edge]) -> crate::shape::Wire {
    TopoBuilder::new().make_wire(edges)
}

/// Unused import guard for RESOLUTION (curve D1 near-zero).
pub fn curve_der_ok(c: &dyn Curve, t: f64) -> bool {
    c.d1(t).1.square_magnitude() > RESOLUTION
}

/// Orientation helper used when assembling split edges (`MakeSplitEdge`).
pub fn vertex_ori_for_split(t1: f64, t2: f64, first: bool) -> Orientation {
    if first {
        if t1 < t2 {
            Orientation::Forward
        } else {
            Orientation::Reversed
        }
    } else if t1 < t2 {
        Orientation::Reversed
    } else {
        Orientation::Forward
    }
}
