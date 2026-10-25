//! The `ChFiDS` fillet-side DS interference statics of
//! `ChFi3d_Builder_0.cxx`: the three interference builders, the "already
//! present" queries, the vertex-on-edge cut and the end-point search.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx` (OCCT
//! 8.0.0):
//!
//! - `ChFi3d_FilCurveInDS` (`cxx:2156-2170`),
//! - `ChFi3d_IndexPointInDS` is in [`super::point_in_ds`], not here,
//! - `ChFi3d_FilPointInDS` (`cxx:2341-2367`),
//! - `ChFi3d_FilVertexInDS` (`cxx:2370-2385`),
//! - `ChFi3d_Orientation` (`cxx:2394-2436`),
//! - `ChFi3d_Contains` (`cxx:2440-2453`),
//! - `QueryAddVertexInEdge` (`cxx:2457-2483`),
//! - `CutEdge` (`cxx:2488-2514`),
//! - `findIndexPoint` (`cxx:2523-2561`).
//!
//! Translation notes:
//! - `occ::handle<TopOpeBRepDS_*Interference>` is
//!   [`InterferenceHandle`]; the subclass of a built handle is the
//!   `InterferenceData` variant (`SurfaceCurve` for
//!   `TopOpeBRepDS_SurfaceCurveInterference`, `CurvePoint` for
//!   `TopOpeBRepDS_CurvePointInterference`), which is what the subclasses
//!   `occ::down_cast` to,
//! - `TopOpeBRepDS_Transition(Et)` is
//!   [`Transition::from_orientation`], and `Transition().Orientation(TopAbs_IN)`
//!   is [`Transition::orientation`] with the `Transition.hxx` default
//!   `TopAbs_FACE` as the second argument,
//! - `DStr.ShapeInterferences(I)` / `DStr.CurveInterferences(I)` /
//!   `DStr.SurfaceInterferences(I)` answer the empty list for an unbound
//!   index, the same way [`DataStructure`] does,
//! - `DStr.ChangeShapeInterferences(I)` answers a throw-away list for an
//!   unbound index; every caller here has bound the shape first with
//!   `DStr.AddShape`, exactly like OCCT,
//! - `BRep_Tool::Parameter(V, E)` is
//!   [`BRepTool::parameter_of_vertex`] and `TopExp_Explorer(E,
//!   TopAbs_VERTEX)` is [`cumulated_children`] over the edge's stored vertex
//!   children, so the orientation of a visited vertex is the composed one OCCT
//!   reports,
//! - `TopoDS::Edge(DStr.Shape(IC))` raises in OCCT when the bound shape is not
//!   an edge; here that is the `Err` of [`cut_edge`].

use std::sync::{Arc, RwLock};

use occt_geom2d::curve::Curve2d;

use crate::abs::{Orientation, ShapeType, State};
use crate::brep_tool::BRepTool;
use crate::chfi_ds::{CommonPoint, SurfData};
use crate::iterator::cumulated_children;
use crate::shape::TopoShape;
use crate::topo_tools::is_same;
use crate::topopebrep::ds::{
    DataStructure, Interference, InterferenceHandle, Kind, Transition,
};

/// `ChFi3d_FilCurveInDS(Icurv, Isurf, Pcurv, Et)` (`cxx:2156-2170`).
///
/// `TopOpeBRepDS_SurfaceCurveInterference(Transition(Et), SURFACE, Isurf,
/// CURVE, Icurv, Pcurv)`: the support is the surface, the geometry the curve.
pub fn ch_fi3d_fil_curve_in_ds(
    icurv: i32,
    isurf: i32,
    p_curv: Option<Arc<dyn Curve2d>>,
    et: Orientation,
) -> InterferenceHandle {
    Arc::new(RwLock::new(Interference::surface_curve(
        Transition::from_orientation(et),
        Kind::Surface,
        isurf,
        Kind::Curve,
        icurv,
        p_curv,
    )))
}

/// `ChFi3d_FilPointInDS(Et, Ic, Ip, Par, IsVertex)` (`cxx:2341-2367`).
///
/// `TopOpeBRepDS_CurvePointInterference(Transition(Et), CURVE, Ic, VERTEX or
/// POINT, Ip, Par)`: the point/vertex kind is the `IsVertex` flag, so a vertex
/// and a point that share the DS index stay distinguishable.
pub fn ch_fi3d_fil_point_in_ds(
    et: Orientation,
    ic: i32,
    ip: i32,
    par: f64,
    is_vertex: bool,
) -> InterferenceHandle {
    let geometry_kind = if is_vertex { Kind::Vertex } else { Kind::Point };
    Arc::new(RwLock::new(Interference::curve_point(
        Transition::from_orientation(et),
        Kind::Curve,
        ic,
        geometry_kind,
        ip,
        par,
    )))
}

/// `ChFi3d_FilVertexInDS(Et, Ic, Ip, Par)` (`cxx:2370-2385`): the same
/// interference as [`ch_fi3d_fil_point_in_ds`] with the geometry kind fixed to
/// `VERTEX`.
pub fn ch_fi3d_fil_vertex_in_ds(
    et: Orientation,
    ic: i32,
    ip: i32,
    par: f64,
) -> InterferenceHandle {
    ch_fi3d_fil_point_in_ds(et, ic, ip, par, true)
}

/// `ChFi3d_Orientation(LI, igros, ipetit, Or, isvertex, aprendre)`
/// (`cxx:2394-2436`).
///
/// Answers the `TopAbs_IN` orientation of the first interference of `li` whose
/// support is `igros` and whose geometry is `ipetit`; `aprendre` additionally
/// requires the geometry kind to be `typepetit`, which is the "is this the
/// point or the vertex?" disambiguation of the comment at `cxx:2400-2401`.
///
/// OCCT reports the answer through the `Or` out-parameter and a `bool`; the
/// port folds both into `Option<Orientation>`.
pub fn ch_fi3d_orientation(
    li: &[InterferenceHandle],
    igros: i32,
    ipetit: i32,
    is_vertex: bool,
    aprendre: bool,
) -> Option<Orientation> {
    let type_petit = if is_vertex { Kind::Vertex } else { Kind::Point };
    for cur in li {
        let Ok(i) = cur.read() else {
            continue;
        };
        let (gk, g, _sk, s) = i.gkgsks();
        let hit = if aprendre {
            s == igros && g == ipetit && gk == type_petit
        } else {
            s == igros && g == ipetit
        };
        if hit {
            return Some(i.transition().orientation(State::In, ShapeType::Face));
        }
    }
    None
}

/// `ChFi3d_Contains(LI, igros, ipetit, isvertex, aprendre)` (`cxx:2440-2453`):
/// "does the interference already exist?", asked through
/// [`ch_fi3d_orientation`] with a throw-away orientation.
pub fn ch_fi3d_contains(
    li: &[InterferenceHandle],
    igros: i32,
    ipetit: i32,
    is_vertex: bool,
    aprendre: bool,
) -> bool {
    ch_fi3d_orientation(li, igros, ipetit, is_vertex, aprendre).is_some()
}

/// `QueryAddVertexInEdge(LI, IC, IV, par, Or)` (`cxx:2457-2483`).
///
/// Appends the `CurvePointInterference` of vertex `IV` at parameter `par` on
/// curve `IC` unless an identical one (same vertex, kind, orientation and a
/// parameter within `1.0e-10`) is already in the list.
pub fn query_add_vertex_in_edge(
    li: &mut Vec<InterferenceHandle>,
    ic: i32,
    iv: i32,
    par: f64,
    or: Orientation,
) {
    for cur in li.iter() {
        let Ok(cpi) = cur.read() else {
            continue;
        };
        if !cpi.is_curve_point() {
            continue;
        }
        let new_iv = cpi.geometry();
        let kv = cpi.geometry_type();
        let new_or = cpi.transition().orientation(State::In, ShapeType::Face);
        let Some(new_par) = cpi.parameter() else {
            continue;
        };
        if iv == new_iv && kv == Kind::Vertex && or == new_or && (par - new_par).abs() < 1.0e-10 {
            return;
        }
    }
    li.push(ch_fi3d_fil_vertex_in_ds(or, ic, iv, par));
}

/// `CutEdge(V, SD, DStr, <unused>, ons)` (`cxx:2488-2514`).
///
/// Injects the vertex `v` into the interference list of the fillet edge stored
/// at `sd`'s `IndexOfC(ons)`, one `CurvePointInterference` per occurrence of
/// `v` among that edge's vertices. The fourth OCCT parameter is unnamed and
/// unused (`cxx:2492`).
pub fn cut_edge(
    v: &TopoShape,
    sd: &SurfData,
    ds: &mut DataStructure,
    ons: i32,
) -> Result<(), String> {
    // `if (!SD->IsOnCurve(ons)) return;` (`cxx:2493-2496`).
    if !sd.is_on_curve(ons) {
        return Ok(());
    }
    // `int IC = SD->IndexOfC(ons);` (`cxx:2497`).
    let ic = sd
        .index_of_c(ons)
        .ok_or_else(|| "CutEdge: the SurfData is not on a curve".to_string())?;
    // `int IV = DStr.AddShape(V);` (`cxx:2498`).
    let iv = ds.add_shape(v);
    // `TopoDS_Edge E = TopoDS::Edge(DStr.Shape(IC)); E.Orientation(FORWARD);`
    // (`cxx:2500-2502`).
    let e = ds
        .shape(ic, false)
        .ok_or_else(|| "CutEdge: the edge index is not bound in the DS".to_string())?;
    if e.shape_type() != ShapeType::Edge {
        return Err("CutEdge: TopoDS::Edge on a non-edge shape".to_string());
    }
    let e = e.oriented(Orientation::Forward);
    // `for (ex.Init(E, TopAbs_VERTEX); ...)` (`cxx:2506-2514`).
    for vv in cumulated_children(&e) {
        if vv.shape_type() != ShapeType::Vertex || !is_same(&vv, v) {
            continue;
        }
        let or = vv.orientation().reversed();
        let par = BRepTool::parameter_of_vertex(&e, &vv);
        let li = ds.change_shape_interferences_at(ic);
        query_add_vertex_in_edge(li, ic, iv, par, or);
    }
    Ok(())
}

/// `findIndexPoint(DStr, Fd, OnS, ipoin)` (`cxx:2523-2561`).
///
/// Searches the `CurvePointInterference`s of the surface interferences of
/// `fd`'s fillet surface for the DS point that coincides with the last common
/// point of `fd` on side `on_s`.
///
/// OCCT raises `Standard_ProgramError` from `DStr.Point(iPoint)` for an
/// out-of-range index; the port answers that as `Err`.
pub fn find_index_point(
    ds: &DataStructure,
    fd: &SurfData,
    on_s: i32,
) -> Result<Option<i32>, String> {
    let p = fd.vertex(false, on_s).point();
    for sci_handle in ds.surface_interferences(fd.surf()) {
        let Ok(sci) = sci_handle.read() else {
            continue;
        };
        // `occ::down_cast<SurfaceCurveInterference>` (`cxx:2537-2538`).
        if sci.as_surface_curve().is_none() {
            continue;
        }
        for cpi_handle in ds.curve_interferences(sci.geometry()) {
            let Ok(cpi) = cpi_handle.read() else {
                continue;
            };
            if !cpi.is_curve_point() {
                continue;
            }
            let i_point = cpi.geometry();
            let tp = ds.point(i_point).ok_or_else(|| {
                format!("findIndexPoint: DS point {i_point} is out of range")
            })?;
            if p.is_equal_with_tolerance(tp.point(), tp.tolerance()) {
                return Ok(Some(i_point));
            }
        }
    }
    Ok(None)
}

/// `ChFiDS_CommonPoint` pair swap of the degenerated-interference pass
/// (`cxx:2645-2649`):
///
/// ```text
/// CP1.Reset(); CP1.SetPoint(CP2.Point());
/// CP2.Reset(); CP2.SetPoint(CP1.Point());
/// ```
///
/// Both points end up with `CP2`'s original coordinates because `CP1` is
/// overwritten before it is read. The Rust version takes the source point by
/// value so the two `&mut` accesses of two different `SurfData` stay
/// sequential.
pub fn reset_both_points(cp1: &mut CommonPoint, cp2: &mut CommonPoint) {
    let p = *cp2.point();
    cp1.reset();
    cp1.set_point(&p);
    cp2.reset();
    cp2.set_point(&p);
}
