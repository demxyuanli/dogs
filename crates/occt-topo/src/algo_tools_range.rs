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
pub fn update_vertex_vv(a_vf: &Vertex, a_new: &Vertex) {
    let a_pvf = BRepTool::vertex_point(a_vf);
    let a_pnew = BRepTool::vertex_point(a_new);
    let a_tol_vf = BRepTool::vertex_tolerance(a_vf);
    let a_tol_new = BRepTool::vertex_tolerance(a_new);
    let a_dist = a_pvf.distance(&a_pnew);
    let a_new_tol = a_dist + a_tol_new;
    if a_new_tol > a_tol_vf {
        a_vf.set_tolerance(a_new_tol + D_TOLERANCE);
    }
}

/// `BOPTools_AlgoTools::UpdateVertex` (edge, parameter, vertex).
pub fn update_vertex_et(a_e: &Edge, a_t: f64, a_v: &Vertex) {
    let a_pv = BRepTool::vertex_point(a_v);
    let a_tol_v = BRepTool::vertex_tolerance(a_v);
    let Some(a_pc) = point_on_edge(a_e, a_t) else {
        return;
    };
    let a_dist = a_pv.distance(&a_pc);
    if a_dist > a_tol_v {
        a_v.set_tolerance(a_dist + D_TOLERANCE);
    }
}

/// `BOPTools_AlgoTools::UpdateVertex` (curve, parameter, vertex).
pub fn update_vertex_ct(curve: &dyn Curve, a_t: f64, a_v: &Vertex) {
    let a_pv = BRepTool::vertex_point(a_v);
    let a_tol_v = BRepTool::vertex_tolerance(a_v);
    let a_pc = curve.d0(a_t);
    let a_dist = a_pv.distance(&a_pc);
    if a_dist > a_tol_v {
        a_v.set_tolerance(a_dist + D_TOLERANCE);
    }
}

/// `BOPTools_AlgoTools::MakeNewVertex` (two vertices).
pub fn make_new_vertex_vv(a_v1: &Vertex, a_v2: &Vertex) -> Result<TopoShape, String> {
    let a_pnt1 = BRepTool::vertex_point(a_v1);
    let a_tol1 = BRepTool::vertex_tolerance(a_v1);
    let a_pnt2 = BRepTool::vertex_point(a_v2);
    let a_tol2 = BRepTool::vertex_tolerance(a_v2);
    let a_dist = a_pnt1.distance(&a_pnt2);
    let a_max_tol = a_tol1.max(a_tol2) + 0.5 * a_dist;
    let a_new = GpPnt::new(
        0.5 * (a_pnt1.x() + a_pnt2.x()),
        0.5 * (a_pnt1.y() + a_pnt2.y()),
        0.5 * (a_pnt1.z() + a_pnt2.z()),
    );
    AlgoTools::make_new_vertex(&a_new, a_max_tol)
}

/// `BOPTools_AlgoTools::MakeNewVertex` (two edges and parameters).
pub fn make_new_vertex_ee(
    a_e1: &Edge,
    a_parm1: f64,
    a_e2: &Edge,
    a_parm2: f64,
) -> Result<TopoShape, String> {
    let a_pnt1 = point_on_edge(a_e1, a_parm1).ok_or("make_new_vertex_ee: no curve on edge 1")?;
    let a_pnt2 = point_on_edge(a_e2, a_parm2).ok_or("make_new_vertex_ee: no curve on edge 2")?;
    let a_tol1 = BRepTool::edge_tolerance(a_e1);
    let a_tol2 = BRepTool::edge_tolerance(a_e2);
    let a_dist = a_pnt1.distance(&a_pnt2);
    let a_max_tol = a_tol1.max(a_tol2) + 0.5 * a_dist;
    let a_new = GpPnt::new(
        0.5 * (a_pnt1.x() + a_pnt2.x()),
        0.5 * (a_pnt1.y() + a_pnt2.y()),
        0.5 * (a_pnt1.z() + a_pnt2.z()),
    );
    AlgoTools::make_new_vertex(&a_new, a_max_tol)
}

/// `BOPTools_AlgoTools::MakeNewVertex` (edge, parameter, face).
pub fn make_new_vertex_ef(a_e1: &Edge, a_parm1: f64, a_f1: &Face) -> Result<TopoShape, String> {
    let a_pnt = point_on_edge(a_e1, a_parm1).ok_or("make_new_vertex_ef: no curve on edge")?;
    let a_tol1 = BRepTool::edge_tolerance(a_e1);
    let a_tol2 = BRepTool::face_tolerance(a_f1);
    let a_max_tol = a_tol1 + a_tol2 + D_TOLERANCE;
    AlgoTools::make_new_vertex(&a_pnt, a_max_tol)
}

/// `BOPTools_AlgoTools::MakeVertex` — fuse a chain of vertices into one.
pub fn make_vertex_from_list(a_lvsd: &[TopoShape]) -> Result<TopoShape, String> {
    if a_lvsd.is_empty() {
        return Err("make_vertex_from_list: empty chain".into());
    }
    let mut cx = 0.0;
    let mut cy = 0.0;
    let mut cz = 0.0;
    let mut tol: f64 = 0.0;
    for s in a_lvsd {
        let v = Vertex(s.clone());
        let p = BRepTool::vertex_point(&v);
        cx += p.x();
        cy += p.y();
        cz += p.z();
        tol = tol.max(BRepTool::vertex_tolerance(&v));
    }
    let n = a_lvsd.len() as f64;
    let p = GpPnt::new(cx / n, cy / n, cz / n);
    for s in a_lvsd {
        let v = Vertex(s.clone());
        let d = BRepTool::vertex_point(&v).distance(&p);
        tol = tol.max(d + BRepTool::vertex_tolerance(&v));
    }
    AlgoTools::make_new_vertex(&p, tol + D_TOLERANCE)
}

/// True when the curve is a straight line (sampled collinearity).
fn curve_is_line(c: &dyn Curve) -> bool {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !a.is_finite() || !b.is_finite() || (b - a).abs() <= 1e-15 {
        return false;
    }
    let p0 = c.d0(a);
    let p1 = c.d0(0.5 * (a + b));
    let p2 = c.d0(b);
    let v1 = occt_core::gp::GpVec::from_pnts(&p0, &p1);
    let v2 = occt_core::gp::GpVec::from_pnts(&p0, &p2);
    let m = v1.crossed(&v2).magnitude();
    let scale = v1.magnitude() * v2.magnitude();
    m < 1e-6 * scale.max(1e-12)
}

/// True when the curve is a free-form type that uses D1 for range correction.
fn curve_is_freeform(c: &dyn Curve) -> bool {
    !curve_is_line(c)
}

/// Approximate `BRepAdaptor_Curve::Resolution(res)` as `res / |D1|`.
fn curve_resolution(c: &dyn Curve, t: f64, res: f64) -> f64 {
    let (_, der) = c.d1(t);
    let mgn = der.magnitude();
    if mgn > 1.0e-12 {
        res / mgn
    } else {
        res
    }
}

/// `BOPTools_AlgoTools::CorrectRange` (edge / edge).
pub fn correct_range_ee(a_e1: &Edge, a_e2: &Edge, a_sr: IntRange) -> IntRange {
    let mut a_new_sr = a_sr;
    let Some(c) = BRepTool::edge_curve(a_e1) else {
        return a_new_sr;
    };
    if curve_is_line(c.as_ref()) {
        return a_new_sr;
    }
    let d_t = PCONFUSION;
    let a_tf = a_sr.first;
    let a_tl = a_sr.last;
    let a_tol_e1 = BRepTool::edge_tolerance(a_e1);
    let a_tol_e2 = BRepTool::edge_tolerance(a_e2);
    let freeform = curve_is_freeform(c.as_ref());
    for i in 0..2 {
        let mut a_res = 2.0 * (a_tol_e1 + a_tol_e2);
        if freeform {
            let t = if i == 0 { a_tf } else { a_tl };
            let (_, der) = c.d1(t);
            let a_mgn = der.magnitude();
            if a_mgn > 1.0e-12 {
                a_res /= a_mgn;
            } else {
                a_res = curve_resolution(c.as_ref(), t, a_res);
            }
        } else {
            let t = if i == 0 { a_tf } else { a_tl };
            a_res = curve_resolution(c.as_ref(), t, a_res);
        }
        if i == 0 {
            a_new_sr.first = a_tf + a_res;
        } else {
            a_new_sr.last = a_tl - a_res;
        }
        if (a_new_sr.last - a_new_sr.first) < d_t {
            a_new_sr = a_sr;
        }
    }
    a_new_sr
}

/// `BOPTools_AlgoTools::CorrectRange` (edge / face).
pub fn correct_range_ef(a_e: &Edge, a_f: &Face, a_sr: IntRange) -> IntRange {
    let mut a_new_sr = a_sr;
    let Some(c) = BRepTool::edge_curve(a_e) else {
        return a_new_sr;
    };
    let d_t = PCONFUSION;
    let a_tf = a_sr.first;
    let a_tl = a_sr.last;
    let a_tol_f = BRepTool::face_tolerance(a_f);
    let freeform = curve_is_freeform(c.as_ref()) && !curve_is_line(c.as_ref());
    for i in 0..2 {
        let mut a_res = a_tol_f;
        if freeform {
            let t = if i == 0 { a_tf } else { a_tl };
            let (_, der) = c.d1(t);
            let a_mgn = der.magnitude();
            if a_mgn > 1.0e-12 {
                a_res /= a_mgn;
            } else {
                a_res = curve_resolution(c.as_ref(), t, a_res);
            }
        } else {
            let t = if i == 0 { a_tf } else { a_tl };
            a_res = curve_resolution(c.as_ref(), t, a_res);
        }
        if i == 0 {
            a_new_sr.first = a_tf + a_res;
        } else {
            a_new_sr.last = a_tl - a_res;
        }
        if (a_new_sr.last - a_new_sr.first) < d_t {
            a_new_sr = a_sr;
        }
    }
    a_new_sr
}

/// `BOPTools_AlgoTools::CopyEdge`.
pub fn copy_edge(the_edge: &Edge) -> Result<Edge, String> {
    let (f, l) = BRepTool::edge_parameters(the_edge);
    let verts: Vec<TopoShape> = ShapeIterator::of_shape(&the_edge.0)
        .filter(|s| s.shape_type() == ShapeType::Vertex)
        .collect();
    let v1 = verts.first();
    let v2 = verts.get(1).or(verts.first());
    let mut sp = AlgoTools::make_split_edge(the_edge, v1, f, v2, l)?;
    sp.0.set_orientation(the_edge.0.orientation());
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
pub fn vertex_parameter(a_c_part: &CommonPrt, param: Option<f64>) -> f64 {
    let a_r = a_c_part.range;
    let mut a_t = 0.5 * (a_r.first + a_r.last);
    if let Some(vp) = param {
        if vp >= a_r.first && vp <= a_r.last {
            a_t = vp;
        }
    }
    a_t
}

/// `IntTools_Tools::VertexParameters`.
pub fn vertex_parameters(a_c_part: &CommonPrt, t1: Option<f64>, t2: Option<f64>) -> (f64, f64) {
    let a_t1 = vertex_parameter(a_c_part, t1);
    let a_r2 = a_c_part.range;
    let mut a_t2 = 0.5 * (a_r2.first + a_r2.last);
    if let Some(vp) = t2 {
        if vp >= a_r2.first && vp <= a_r2.last {
            a_t2 = vp;
        }
    }
    (a_t1, a_t2)
}

/// `IntTools_Tools::IsOnPave1`.
pub fn is_on_pave1(a_tr: f64, a_cp_first: f64, a_cp_last: f64, a_tolerance: f64) -> bool {
    if a_tr >= a_cp_first && a_tr <= a_cp_last {
        return true;
    }
    let d_t1 = (a_tr - a_cp_first).abs();
    let d_t2 = (a_tr - a_cp_last).abs();
    d_t1 <= a_tolerance || d_t2 <= a_tolerance
}

/// `IntTools_Tools::IsInRange`.
pub fn is_in_range(a_r_ref: IntRange, a_r: IntRange, a_tolerance: f64) -> bool {
    let a_t1 = a_r.first;
    let a_t2 = a_r.last;
    let a_t_ref1 = a_r_ref.first - a_tolerance;
    let a_t_ref2 = a_r_ref.last + a_tolerance;
    (a_t1 >= a_t_ref1 && a_t1 <= a_t_ref2) || (a_t2 >= a_t_ref1 && a_t2 <= a_t_ref2)
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
