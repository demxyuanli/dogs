//! The closure / orientation-in-parent helpers of `TopOpeBRepTool_TOOL`.
//!
//! Source: `TKBool/TopOpeBRepTool/TopOpeBRepTool_TOOL.cxx:94-436`,
//! `TopOpeBRepTool_TOOL.hxx:54-110`.
//!
//! `OriinSor` and friends return `int` codes rather than a `TopAbs_Orientation`
//! because the `#define`s at `cxx:79-83` add two values the enum does not have
//! (`FORWARD (1)`, `REVERSED (2)`, `INTERNAL (3)`, `EXTERNAL (4)`,
//! `CLOSING (5)`), and `0` means "not found". [`ori_code`] maps the enum onto
//! the four shared codes, leaving `0` and `CLOSING` unambiguous.

use std::sync::Arc;

use occt_geom::surface::Surface;

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::iterator::cumulated_children;
use crate::shape::{Edge, Face, TopoShape};
use crate::topexp::Explorer;
use crate::topo_tools;

use super::adaptor::CurveAdaptor;

/// `FORWARD (1)` (`cxx:79`).
pub const FORWARD: i32 = 1;
/// `REVERSED (2)` (`cxx:80`).
pub const REVERSED: i32 = 2;
/// `INTERNAL (3)` (`cxx:81`).
pub const INTERNAL: i32 = 3;
/// `EXTERNAL (4)` (`cxx:82`).
pub const EXTERNAL: i32 = 4;
/// `CLOSING (5)` (`cxx:83`).
pub const CLOSING: i32 = 5;

/// `TopAbs_Orientation` to its `FORWARD..EXTERNAL` code.
pub fn ori_code(o: Orientation) -> i32 {
    match o {
        Orientation::Forward => FORWARD,
        Orientation::Reversed => REVERSED,
        Orientation::Internal => INTERNAL,
        Orientation::External => EXTERNAL,
    }
}

/// `static bool CheckEdgeLength(E)` (`cxx:96-136`).
///
/// Only a single-vertex edge is measured: its ten samples must span more than
/// `Precision::Confusion()`. Every other edge returns `true` unmeasured.
pub fn check_edge_length(e: &TopoShape) -> bool {
    if unique_vertex_count(e) != 1 {
        return true;
    }
    let Some(curve) = BRepTool::edge_curve_adaptor_world(&Edge(e.clone())) else {
        return true;
    };
    let (f, l) = BRepTool::edge_parameters(&Edge(e.clone()));
    let n = 10;
    let dt = (l - f) / (n as f64);
    let mut ln = 0.0;
    let mut p1 = curve.d0(f);
    for i in 1..=n {
        let p2 = if i == n {
            curve.d0(l)
        } else {
            curve.d0(f + (i as f64) * dt)
        };
        ln += p1.distance(&p2);
        p1 = p2;
    }
    ln > occt_core::precision::CONFUSION
}

/// `TopExp::MapShapes(E, TopAbs_VERTEX, M).Extent()` — occurrences deduplicated
/// by `IsSame`.
fn unique_vertex_count(e: &TopoShape) -> usize {
    let mut seen: Vec<TopoShape> = Vec::new();
    let mut ex = Explorer::new(e, ShapeType::Vertex);
    while ex.more() {
        let v = ex.current().clone();
        if !seen.iter().any(|s| topo_tools::is_same(s, &v)) {
            seen.push(v);
        }
        ex.next();
    }
    seen.len()
}

/// `TopOpeBRepTool_TOOL::OriinSor(sub, S, checkclo)` (`cxx:142-212`).
///
/// With `check_clo`, a shape closed in itself answers `CLOSING` for its closing
/// vertex (`S` an edge) or its closing edge (`S` a face) before the occurrence
/// scan runs.
pub fn ori_in_sor(sub: &TopoShape, s: &TopoShape, check_clo: bool) -> i32 {
    if check_clo {
        match s.shape_type() {
            ShapeType::Edge => {
                if sub.shape_type() != ShapeType::Vertex {
                    return 0;
                }
                if let Some(vclo) = closed_e(s) {
                    if topo_tools::is_same(sub, &vclo) {
                        return CLOSING;
                    }
                }
            }
            ShapeType::Face => {
                if sub.shape_type() != ShapeType::Edge {
                    return 0;
                }
                if closed_s(s) && is_closing_e(sub, s) {
                    return CLOSING;
                }
            }
            _ => {}
        }
    }

    let mut ex = Explorer::new(s, sub.shape_type());
    while ex.more() {
        let ssub = ex.current();
        if topo_tools::is_same(ssub, sub) {
            return ori_code(ssub.orientation());
        }
        ex.next();
    }
    0
}

/// `TopOpeBRepTool_TOOL::OriinSorclosed(sub, S)` (`cxx:216-253`): the same over
/// `TopoDS_Iterator`, i.e. direct children only, compared with `IsEqual`
/// (location and orientation both matter).
pub fn ori_in_sor_closed(sub: &TopoShape, s: &TopoShape) -> i32 {
    match s.shape_type() {
        ShapeType::Edge => {
            if sub.shape_type() != ShapeType::Vertex {
                return 0;
            }
        }
        ShapeType::Face => {
            if sub.shape_type() != ShapeType::Edge {
                return 0;
            }
        }
        _ => {}
    }
    for ssub in cumulated_children(s) {
        if !topo_tools::is_equal(&ssub, sub) {
            continue;
        }
        if ssub.orientation() == Orientation::Forward {
            return FORWARD;
        }
        if ssub.orientation() == Orientation::Reversed {
            return REVERSED;
        }
    }
    0
}

/// `TopOpeBRepTool_TOOL::ClosedE(E, vclo)` (`cxx:256-286`): the closing vertex of
/// `E`, or `None`.
///
/// `INTERNAL` vertices are skipped, and a degenerated edge has no closing
/// vertex. The first non-`INTERNAL` vertex is the one the later occurrences are
/// compared against.
pub fn closed_e(e: &TopoShape) -> Option<TopoShape> {
    if BRepTool::is_degenerated(&Edge(e.clone())) {
        return None;
    }
    let mut vv: Option<TopoShape> = None;
    let mut ex = Explorer::new(e, ShapeType::Vertex);
    while ex.more() {
        let v = ex.current().clone();
        ex.next();
        if v.orientation() == Orientation::Internal {
            continue;
        }
        match &vv {
            None => vv = Some(v),
            Some(first) => {
                if topo_tools::is_same(&v, first) {
                    return Some(first.clone());
                }
            }
        }
    }
    None
}

/// `TopOpeBRepTool_ShapeTool::BASISSURFACE(S)` (`ShapeTool.cxx:122-141`): unwrap
/// `Geom_OffsetSurface` / `Geom_RectangularTrimmedSurface` down to the surface
/// whose closure flags are read.
pub fn basis_surface(s: &Arc<dyn Surface>) -> Arc<dyn Surface> {
    if let Some(basis) = s.offset_basis_surface() {
        return basis_surface(&basis);
    }
    if let Some(basis) = s.rectangular_trimmed_basis() {
        return basis_surface(&basis);
    }
    s.clone()
}

/// `TopOpeBRepTool_TOOL::ClosedS(F)` (`cxx:291-309`).
///
/// A closed direction must also be periodic: `IsUClosed` alone is not enough.
pub fn closed_s(f: &TopoShape) -> bool {
    let Some(surface) = BRepTool::face_surface_world(&Face(f.clone())) else {
        return false;
    };
    let basis = basis_surface(&surface);
    let mut uclosed = basis.is_u_closed();
    if uclosed {
        uclosed = basis.is_u_periodic();
    }
    let mut vclosed = basis.is_v_closed();
    if vclosed {
        vclosed = basis.is_v_periodic();
    }
    uclosed || vclosed
}

/// `TopOpeBRepTool_TOOL::IsClosingE(E, F)` (`cxx:313-328`).
pub fn is_closing_e(e: &TopoShape, f: &TopoShape) -> bool {
    closing_e_in(e, f, f)
}

/// `TopOpeBRepTool_TOOL::IsClosingE(E, W, F)` (`cxx:333-351`): the occurrences
/// are counted in `W`, the closure test still uses `F`.
pub fn is_closing_e_in_wire(e: &TopoShape, w: &TopoShape, f: &TopoShape) -> bool {
    closing_e_in(e, w, f)
}

/// `FUN_tool_IsClosingE(E, S, F)` (`TopOpeBRepTool_TOPOLOGY.cxx:779-797`): the
/// occurrences of `E` are counted in the arbitrary sub-shape `S` (a wire or a
/// face), the closure test uses `F`.
pub fn is_closing_e_in_shape(e: &TopoShape, s: &TopoShape, f: &TopoShape) -> bool {
    closing_e_in(e, s, f)
}

fn closing_e_in(e: &TopoShape, w: &TopoShape, f: &TopoShape) -> bool {
    let mut nbocc = 0;
    let mut ex = Explorer::new(w, ShapeType::Edge);
    while ex.more() {
        if topo_tools::is_same(ex.current(), e) {
            nbocc += 1;
        }
        ex.next();
    }
    if nbocc != 2 {
        return false;
    }
    BRepTool::is_closed_edge_face(&Edge(e.clone()), &Face(f.clone()))
}

/// `TopOpeBRepTool_TOOL::Vertices(E, Vces)` (`cxx:355-383`): `Vces[1]` holds the
/// vertex at the smaller parameter, `Vces[2]` the other one, with the pair
/// swapped again for a `REVERSED` edge handle.
pub fn vertices(e: &TopoShape) -> (Option<TopoShape>, Option<TopoShape>) {
    let ori_e = e.orientation();
    let (v1, v2) = topo_tools::edge_vertices(&Edge(e.clone()));
    let v1 = v1.map(|v| v.0);
    let v2 = v2.map(|v| v.0);
    let par = |v: &Option<TopoShape>| match v {
        Some(v) => BRepTool::parameter_of_vertex(e, v),
        None => f64::NAN,
    };
    let (par1, par2) = (par(&v1), par(&v2));

    let mut iv_par_sma = if par1 < par2 { FORWARD } else { REVERSED };
    let mut iv_par_sup = if par1 < par2 { REVERSED } else { FORWARD };
    if ori_e == Orientation::Reversed {
        iv_par_sma = if iv_par_sma == FORWARD { REVERSED } else { FORWARD };
        iv_par_sup = if iv_par_sup == REVERSED { FORWARD } else { REVERSED };
    }
    let mut vces: [Option<TopoShape>; 2] = [None, None];
    vces[(iv_par_sma - 1) as usize] = v1;
    vces[(iv_par_sup - 1) as usize] = v2;
    let [a, b] = vces;
    (a, b)
}

/// `TopOpeBRepTool_TOOL::Vertex(Iv, E)` (`cxx:387-394`).
pub fn vertex(iv: i32, e: &TopoShape) -> Option<TopoShape> {
    let (a, b) = vertices(e);
    if iv == 1 {
        a
    } else {
        b
    }
}

/// `TopOpeBRepTool_TOOL::ParE(Iv, E)` (`cxx:397-402`).
pub fn par_e(iv: i32, e: &TopoShape) -> Option<f64> {
    let v = vertex(iv, e)?;
    Some(BRepTool::parameter_of_vertex(e, &v))
}

/// `TopOpeBRepTool_TOOL::OnBoundary(par, e)` (`cxx:406-433`).
pub fn on_boundary(par: f64, e: &TopoShape) -> i32 {
    let adaptor = CurveAdaptor::of_edge(&Edge(e.clone()));
    let closed = adaptor.is_closed();
    let first = adaptor.first_parameter();
    let last = adaptor.last_parameter();
    let tolp = adaptor.resolution(BRepTool::edge_tolerance(&Edge(e.clone())));

    let onf = (par - first).abs() < tolp;
    let onl = (par - last).abs() < tolp;
    if (onf || onl) && closed {
        return CLOSING;
    }
    if onf {
        return FORWARD;
    }
    if onl {
        return REVERSED;
    }
    if first < par && par < last {
        return INTERNAL;
    }
    EXTERNAL
}
