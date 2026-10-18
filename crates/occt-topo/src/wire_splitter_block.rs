//! `BOPAlgo_WireSplitter::SplitBlock` (`BOPAlgo_WireSplitter_1.cxx`).
//!
//! Irregular connexity blocks (a vertex with other than one IN and one OUT,
//! or TShape-coincident edges) are walked by `Path`: unused OUT edges, the
//! smallest `ClockWiseAngle` from the incoming IN angle, then `MakeWire` of
//! each closed buffer. Regular 1-in/1-out blocks without TShape coincidence
//! fall through to a single `MakeWire`.

use std::collections::HashMap;
use std::f64::consts::PI;

use occt_core::gp::{GpDir2d, GpPnt2d, GpVec2d};
use occt_core::precision::{Precision, ANGULAR, PCONFUSION};
use occt_geom2d::cl_props2d::GeomLPropCLProps2d;
use occt_geom2d::line::Geom2dLine;

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d::make_2d;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::{is_equal, is_same};

const TWO_PI: f64 = PI + PI;

/// `BOPAlgo_EdgeInfo`.
#[derive(Clone)]
struct EdgeInfo {
    edge: Edge,
    passed: bool,
    is_in: bool,
    is_inside: bool,
    angle: f64,
}

struct SmartMap {
    verts: Vec<TopoShape>,
    infos: Vec<Vec<EdgeInfo>>,
    index: HashMap<usize, usize>,
}

impl SmartMap {
    fn new() -> Self {
        Self {
            verts: Vec::new(),
            infos: Vec::new(),
            index: HashMap::new(),
        }
    }

    fn slot(&mut self, v: &TopoShape) -> usize {
        let k = GeometryRegistry::shape_key(v);
        if let Some(&i) = self.index.get(&k) {
            return i;
        }
        let i = self.verts.len();
        self.index.insert(k, i);
        self.verts.push(v.clone());
        self.infos.push(Vec::new());
        i
    }
}

/// `BOPAlgo_WireSplitter::SplitBlock`.
pub(crate) fn split_wire_block(face: &Face, edges: &[Edge]) -> Result<Vec<TopoShape>, String> {
    let mut map = SmartMap::new();
    let mut unique: HashMap<usize, usize> = HashMap::new();
    let mut vert_closed: HashMap<usize, bool> = HashMap::new();
    let fk = GeometryRegistry::shape_key(&face.0);

    for e in edges {
        if make_2d(e, face).is_err() {
            continue;
        }
        let ek = GeometryRegistry::shape_key(&e.0);
        let closed = edge_closed_on_face(e, face, fk);
        let n = unique.entry(ek).or_insert(0);
        *n += 1;
        if *n > 1 && !closed {
            unique.insert(ek, 0);
        }
        let kids = oriented_vertices(e);
        for (i, v) in kids.iter().enumerate() {
            let slot = map.slot(v);
            let is_in = v.orientation() == Orientation::Reversed;
            map.infos[slot].push(EdgeInfo {
                edge: e.clone(),
                passed: false,
                is_in,
                is_inside: false,
                angle: -1.0,
            });
            let vk = GeometryRegistry::shape_key(v);
            let cyc = closed || (kids.len() >= 2 && is_same(&kids[0], &kids[1]));
            let e0 = vert_closed.entry(vk).or_insert(false);
            if cyc || i > 0 && is_same(&kids[0], v) {
                *e0 = true;
            }
            if cyc {
                *e0 = true;
            }
        }
    }

    if map.verts.is_empty() {
        return Err("SplitBlock: no edges with a pcurve on the face".into());
    }

    for infos in &mut map.infos {
        for ei in infos.iter_mut() {
            let ek = GeometryRegistry::shape_key(&ei.edge.0);
            ei.is_inside = unique.get(&ek).copied().unwrap_or(0) == 0;
        }
    }

    if nothing_to_do(&map, edges) {
        return Ok(vec![make_wire(edges)]);
    }

    for slot in 0..map.verts.len() {
        let v = map.verts[slot].clone();
        let n = map.infos[slot].len();
        for i in 0..n {
            let e = map.infos[slot][i].edge.clone();
            let is_in = map.infos[slot][i].is_in;
            let mut vv = v.clone();
            vv.set_orientation(if is_in {
                Orientation::Reversed
            } else {
                Orientation::Forward
            });
            map.infos[slot][i].angle = angle_2d(&vv, &e, face, is_in);
        }
    }
    refine_angles(face, &mut map);

    let mut loops: Vec<TopoShape> = Vec::new();
    let nverts = map.verts.len();
    for slot in 0..nverts {
        let ninfo = map.infos[slot].len();
        for i in 0..ninfo {
            if map.infos[slot][i].is_in || map.infos[slot][i].passed {
                continue;
            }
            path(
                face,
                &vert_closed,
                slot,
                i,
                &mut map,
                &mut loops,
            );
        }
    }
    if loops.is_empty() {
        return Err("SplitBlock: Path produced no wires".into());
    }
    Ok(loops)
}

fn nothing_to_do(map: &SmartMap, edges: &[Edge]) -> bool {
    for infos in &map.infos {
        let mut cin = 0;
        let mut cout = 0;
        for ei in infos {
            if ei.is_in {
                cin += 1;
            } else {
                cout += 1;
            }
        }
        if cin != 1 || cout != 1 {
            return false;
        }
    }
    let mut n_by_tshape: HashMap<usize, usize> = HashMap::new();
    for e in edges {
        *n_by_tshape
            .entry(GeometryRegistry::shape_key(&e.0))
            .or_insert(0) += 1;
    }
    n_by_tshape.values().all(|&n| n == 1)
}

fn path(
    face: &Face,
    vert_closed: &HashMap<usize, bool>,
    start_slot: usize,
    start_info: usize,
    map: &mut SmartMap,
    loops: &mut Vec<TopoShape>,
) {
    let mut ls: Vec<Edge> = Vec::new();
    let mut vert_va: Vec<TopoShape> = Vec::new();
    let mut coord_va: Vec<GpPnt2d> = Vec::new();
    let mut info_seq: Vec<(usize, usize)> = Vec::new();

    let mut va_slot = start_slot;
    let mut info_idx = start_info;
    // `BOPAlgo_WireSplitter_1.cxx:518` (`aEOuta = aLS.Last()`): after the
    // closing scan cuts the path back, the incoming edge is the truncated
    // path's last edge; `None` means "the edge appended in this iteration".
    let mut e_in: Option<Edge> = None;
    let eps = f64::EPSILON;
    loop {
        if ls.len() == 1 && is_same(&ls[0].0, &map.infos[va_slot][info_idx].edge.0) {
            return;
        }
        map.infos[va_slot][info_idx].passed = true;
        // Only the closing scan below may rewrite the incoming edge.
        e_in = None;
        let e_out = map.infos[va_slot][info_idx].edge.clone();
        let va = map.verts[va_slot].clone();
        ls.push(e_out.clone());
        vert_va.push(va.clone());
        info_seq.push((va_slot, info_idx));

        let mut p_va = va.clone();
        p_va.set_orientation(Orientation::Forward);
        let pa = coord2d(&p_va, &e_out, face);
        coord_va.push(pa);

        let vb = next_vertex(&p_va, &e_out);
        let pb = coord2d(&vb, &e_out, face);
        let Some(&vb_slot) = map.index.get(&GeometryRegistry::shape_key(&vb)) else {
            return;
        };
        let tol2d = 2.0 * tolerance_2d(&Vertex(vb.clone()), face);
        let tol2d2 = tol2d * tol2d;
        let closed = vert_closed
            .get(&GeometryRegistry::shape_key(&vb))
            .copied()
            .unwrap_or(false);

        {
            let mut buf: Vec<Edge> = Vec::new();
            let mut has_edge = false;
            let nb = ls.len();
            let mut cut: Option<usize> = None;
            for i in (0..nb).rev() {
                let e_prev = &ls[i];
                buf.push(e_prev.clone());
                if !has_edge {
                    has_edge = !BRepTool::is_degenerated(e_prev);
                    if !has_edge {
                        continue;
                    }
                }
                let same_v = is_same(&vert_va[i], &vb);
                let mut same_v2d = same_v;
                if same_v && closed {
                    same_v2d = coord_va[i].square_distance(&pb) < tol2d2;
                    if same_v2d {
                        let ud = (coord_va[i].x() - pb.x()).abs();
                        let vd = (coord_va[i].y() - pb.y()).abs();
                        let (tu, tv) = uv_tolerance_2d(&Vertex(vb.clone()), face);
                        if ud > 2.0 * tu || vd > 2.0 * tv {
                            same_v2d = false;
                        }
                    }
                }
                if same_v && same_v2d {
                    let mut priz = true;
                    if buf.len() == 2 && is_same(&buf[0].0, &buf[1].0) {
                        priz = false;
                    }
                    if priz {
                        loops.push(make_wire(&buf));
                    }
                    if i < 1 {
                        return;
                    }
                    cut = Some(i);
                    break;
                }
            }
            if let Some(i) = cut {
                ls.truncate(i);
                vert_va.truncate(i);
                coord_va.truncate(i);
                info_seq.truncate(i);
                // `BOPAlgo_WireSplitter_1.cxx:518`: the incoming edge for the
                // next choice is the truncated path's last edge, not the edge
                // that just closed the loop.
                e_in = ls.last().cloned();
            }
        }

        let is_boundary = info_seq
            .last()
            .map(|&(s, ii)| !map.infos[s][ii].is_inside)
            .unwrap_or(true);
        let le = &map.infos[vb_slot];
        let angle_in = angle_in(e_in.as_ref().unwrap_or(&e_out), le);
        let i_cnt = nb_ways_out(le);
        let mut min_angle = 100.0;
        let mut chosen: Option<usize> = None;
        let mut only_inside: Option<usize> = None;
        let mut n_inside = 0usize;

        for (k, ei) in le.iter().enumerate() {
            if ei.is_in || ei.passed {
                continue;
            }
            if i_cnt == 0 {
                return;
            }
            if i_cnt == 1 {
                chosen = Some(k);
                break;
            }
            let ang = if is_same(&ei.edge.0, &e_out.0) {
                TWO_PI
            } else {
                if closed {
                    let p2 = coord2d_vf(&ei.edge, face);
                    if p2.square_distance(&pb) > tol2d2 {
                        continue;
                    }
                }
                clockwise_angle(angle_in, ei.angle)
            };
            if is_boundary && ei.is_inside {
                n_inside += 1;
                only_inside = Some(k);
            }
            if ang < min_angle - eps {
                min_angle = ang;
                chosen = Some(k);
            }
        }
        if n_inside == 1 {
            chosen = only_inside;
        }
        let Some(k) = chosen else {
            return;
        };
        va_slot = vb_slot;
        info_idx = k;
    }
}

fn make_wire(edges: &[Edge]) -> TopoShape {
    let w = TopoBuilder::new().make_wire(edges);
    w.0.set_closed(true);
    w.0
}

/// `ClockWiseAngle`.
fn clockwise_angle(angle_in: f64, angle_out: f64) -> f64 {
    let mut in_value = angle_in;
    let mut out = angle_out;
    if in_value >= TWO_PI {
        in_value -= TWO_PI;
    }
    if out >= TWO_PI {
        out -= TWO_PI;
    }
    let mut a1 = in_value + PI;
    if a1 >= TWO_PI {
        a1 -= TWO_PI;
    }
    let mut da = a1 - out;
    if da <= 0.0 {
        da += TWO_PI;
    } else if da <= 1.0e-14 {
        da = TWO_PI;
    }
    da
}

/// `gp_Dir2d::Angle` from +X, folded into `[0, 2pi)`.
fn dir2d_angle(d: &GpDir2d) -> f64 {
    let a = d.y().atan2(d.x());
    if a < 0.0 {
        a + TWO_PI
    } else {
        a
    }
}

fn coord2d(v: &TopoShape, e: &Edge, face: &Face) -> GpPnt2d {
    let t = vertex_parameter_on_face(v, e);
    match make_2d(e, face) {
        Ok(c) => c.d0(t),
        Err(_) => GpPnt2d::new(99.0, 99.0),
    }
}

/// `BRep_Tool::Parameter(V, E, S, L)` (`BRep_Tool.cxx:301-352`): the parameter of
/// the vertex `v` on the edge's pcurve.
///
/// The occurrence of `v` among the edge's vertices is searched with
/// `TopoDS_Iterator(E.Oriented(FORWARD))`, i.e. over the *stored* vertex
/// orientations; a second occurrence of the same vertex (a closed edge whose two
/// ends are one shape, e.g. a ring on a periodic surface) records
/// `rev = (E.Orientation() == REVERSED)` and replaces the occurrence only when
/// its orientation equals `v`'s. The chosen occurrence's orientation then picks
/// the end: `FORWARD -> first` (flipped by `rev`), `REVERSED -> last` (flipped by
/// `rev`). Without this, the two ends of a closed edge collapse onto the same UV
/// and `BOPAlgo_WireSplitter::Path` closes a bogus one-edge loop at the seam.
///
/// The `INTERNAL` / not-found branch of OCCT consults the vertex's
/// `BRep_PointRepresentation`s on the pcurve (`BRep_TVertex::Points`), which this
/// port does not model; it falls back to the 3D `BRep_Tool::Parameter(V, E)`.
fn vertex_parameter_on_face(v: &TopoShape, e: &Edge) -> f64 {
    let (first, last) = BRepTool::edge_parameters(e);
    let mut vf: Option<Orientation> = None;
    let mut rev = false;
    for k in stored_vertices_oriented(e) {
        if !is_same(&k, v) {
            continue;
        }
        match vf {
            None => vf = Some(k.orientation()),
            Some(_) => {
                rev = e.0.orientation() == Orientation::Reversed;
                if k.orientation() == v.orientation() {
                    vf = Some(k.orientation());
                }
            }
        }
    }
    match vf {
        Some(Orientation::Forward) => {
            if rev {
                last
            } else {
                first
            }
        }
        Some(Orientation::Reversed) => {
            if rev {
                first
            } else {
                last
            }
        }
        _ => vertex_parameter(v, e),
    }
}

fn coord2d_vf(e: &Edge, face: &Face) -> GpPnt2d {
    for v in oriented_vertices(e) {
        if v.orientation() == Orientation::Forward {
            return coord2d(&v, e, face);
        }
    }
    GpPnt2d::new(99.0, 99.0)
}

fn vertex_parameter(v: &TopoShape, e: &Edge) -> f64 {
    let (first, last) = BRepTool::edge_parameters(e);
    let kids = stored_vertices(e);
    if kids.first().is_some_and(|k| is_same(k, v)) {
        first
    } else {
        last
    }
}

fn next_vertex(v: &TopoShape, e: &Edge) -> TopoShape {
    for vx in oriented_vertices(e) {
        if !is_equal(&vx, v) {
            return vx;
        }
    }
    v.clone()
}

fn nb_ways_out(le: &[EdgeInfo]) -> usize {
    le.iter().filter(|ei| !ei.is_in && !ei.passed).count()
}

fn angle_in(e_in: &Edge, le: &[EdgeInfo]) -> f64 {
    for ei in le {
        if ei.is_in && is_equal(&ei.edge.0, &e_in.0) {
            return ei.angle;
        }
    }
    0.0
}

fn angle_2d(v: &TopoShape, e: &Edge, face: &Face, is_in: bool) -> f64 {
    let tv = vertex_parameter_on_face(v, e);
    if Precision::is_infinite(tv) {
        return 0.0;
    }
    let Ok(c) = make_2d(e, face) else {
        return 0.0;
    };
    let (first, last) = BRepTool::edge_parameters(e);
    let tol2d = 2.0 * tolerance_2d(&Vertex(v.clone()), face);
    // `dt = std::max(aGAC2D.Resolution(tol2d), Precision::PConfusion())`
    // (`BOPAlgo_WireSplitter_1.cxx:792`), with `aType = aGAC2D.GetType()`
    // (`:795`) driving the non-line curvature correction (`:796-807`).
    let mut dt = c.resolution(tol2d).max(PCONFUSION);
    if !c.is_line() {
        // `GeomLProp_CLProps2d LProp(aC2D, aTV, 2, Precision::PConfusion())`
        // then `IsTangentDefined()` / `Curvature()` (`_1.cxx:797-807`).
        let mut lprop = GeomLPropCLProps2d::new(c.as_ref(), tv, 2, PCONFUSION);
        if lprop.is_tangent_defined() {
            let mut r = lprop.curvature();
            if r > PCONFUSION {
                r = 1.0 / r;
                let cosphi = r / (r + tol2d);
                dt = dt.max(cosphi.acos());
            }
        }
    }
    let mut tx = 0.05 * (last - first);
    if tx < 5.0e-5 {
        tx = 5.0e-5_f64.min((last - first) / 2.0);
    }
    if dt > tx {
        dt = tx;
    }
    let tv1 = if (tv - first).abs() < (tv - last).abs() {
        tv + dt
    } else {
        tv - dt
    };
    let pv = c.d0(tv);
    let pv1 = c.d0(tv1);
    let v2 = if is_in {
        GpVec2d::new(pv.x() - pv1.x(), pv.y() - pv1.y())
    } else {
        GpVec2d::new(pv1.x() - pv.x(), pv1.y() - pv.y())
    };
    match GpDir2d::from_vec2d(&v2) {
        Ok(d) => dir2d_angle(&d),
        Err(_) => 0.0,
    }
}

/// `BOPAlgo_WireSplitter::Tolerance2D` (`BOPAlgo_WireSplitter_1.cxx:859-881`):
/// `max(UResolution, VResolution, tolerance)` with the 1.1 factor of
/// `GeomAbs_BSplineSurface`.
fn tolerance_2d(v: &Vertex, face: &Face) -> f64 {
    let t3 = BRepTool::vertex_tolerance(v);
    let (u, vv) = uv_tolerance_2d(v, face);
    let mut t2 = u.max(vv).max(t3);
    if let Some(s) = BRepTool::face_surface(face) {
        if crate::geom_bnd_lib_surface3d::surface_kind(s.as_ref())
            == crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType::BSplineSurface
        {
            // `cxx:875-878`.
            t2 *= 1.1;
        }
    }
    t2
}

/// `BOPAlgo_WireSplitter::UTolerance2D` / `VTolerance2D` (`_1.cxx:885-901`):
/// `BRepAdaptor_Surface::U/VResolution` (`GeomAdaptor_Surface.cxx:1818-1958`).
///
/// OCCT has one implementation of that switch and the port keeps one too, in
/// [`occt_geom::approx_same_parameter`]; this forwards to it rather than
/// re-deriving the per-class arms (the former hand-rolled version left the
/// revolution / extrusion arms on a finite-difference estimate, missed the
/// `Geom_RectangularTrimmedSurface` unwrap of `cxx:423-430` and the
/// `GeomAbs_OffsetSurface` arms of `cxx:1882-1885` / `:1944-1947`, and had no
/// unbounded-cone guard for `cxx:1855-1858`).
fn uv_tolerance_2d(v: &Vertex, face: &Face) -> (f64, f64) {
    let t3 = BRepTool::vertex_tolerance(v);
    let Some(surf) = BRepTool::face_surface(face) else {
        return (t3, t3);
    };
    let s = surf.as_ref();
    (
        occt_geom::approx_same_parameter::u_resolution(s, t3),
        occt_geom::approx_same_parameter::v_resolution(s, t3),
    )
}

/// `BOPAlgo_WireSplitter::RefineAngles` (`BOPAlgo_WireSplitter_1.cxx:905-918`,
/// `:925-1029`): for every vertex, when exactly **two** boundary edges meet,
/// re-angle the interior OUT edges that do not already lie inside the boundary
/// wedge — through [`refine_angle_2d`], or, when the vertex carries exactly two
/// interior edges and the 2-D intersection finds nothing, by nudging the angle
/// just inside the wedge (`:996-1000`). The updated angle is looked up by edge
/// (`aDMSR` is a map of edge → angle), so both occurrences of a shared section
/// edge at the vertex get it, the IN one shifted by π (`:1009-1027`).
fn refine_angles(face: &Face, map: &mut SmartMap) {
    for slot in 0..map.verts.len() {
        let v = map.verts[slot].clone();
        refine_angles_at_vertex(face, &v, &mut map.infos[slot]);
    }
}

fn refine_angles_at_vertex(face: &Face, v: &TopoShape, infos: &mut [EdgeInfo]) {
    let mut a1 = 0.0_f64; // angle of the outgoing boundary edge
    let mut a2 = 0.0_f64; // angle of the incoming boundary edge
    let mut i_cnt_bnd = 0usize;
    let mut i_cnt_int = 0usize;
    for ei in infos.iter() {
        if !ei.is_inside {
            i_cnt_bnd += 1;
            if ei.is_in {
                a2 = ei.angle;
            } else {
                a1 = ei.angle;
            }
        } else {
            i_cnt_int += 1;
        }
    }
    if i_cnt_bnd != 2 {
        return;
    }
    let delta = clockwise_angle(a2, a1);
    let mut refined: HashMap<usize, f64> = HashMap::new();
    for ei in infos.iter() {
        if !ei.is_inside || ei.is_in {
            continue;
        }
        if clockwise_angle(a2, ei.angle) < delta {
            continue; // already inside the wedge
        }
        let ek = GeometryRegistry::shape_key(&ei.edge.0);
        match refine_angle_2d(face, v, &ei.edge, a1, a2, delta) {
            Some(new_val) => {
                refined.insert(ek, new_val);
            }
            None if i_cnt_int == 2 => {
                let new_val = if ei.angle <= a1 {
                    a1 + ANGULAR
                } else {
                    a2 - ANGULAR
                };
                refined.insert(ek, new_val);
            }
            None => {}
        }
    }
    if refined.is_empty() {
        return;
    }
    for ei in infos.iter_mut() {
        if let Some(&a) = refined.get(&GeometryRegistry::shape_key(&ei.edge.0)) {
            ei.angle = if ei.is_in { a + PI } else { a };
        }
    }
}

/// `RefineAngle2D` (`BOPAlgo_WireSplitter_1.cxx:1033-1125`): intersect the
/// edge's p-curve with the ray through the vertex along each boundary angle and
/// take the angle of the point where the curve leaves the wedge. Returns the
/// refined angle, or `None` when neither ray produces a usable point.
/// `BOPAlgo_WireSplitter_1.cxx::RefineAngle2D` (`:1033-1125`).
///
/// Substitutions registered for audit §13 / task T-33:
/// * OCCT intersects through `Geom2dInt_GInter::Perform(aGAC1, aDomain1,
///   aGAC2, aDomain2, aTolInt, aTolInt)` with an explicit `IntRes2d_Domain`
///   (`cxx:1080`); this port calls
///   [`occt_geom2d::geom2d_api::intersect_curves`] and applies the
///   `[t1, t2]` restriction itself (`:1084-1100`). The two agree when the
///   intersector returns every intersection point of that domain.
fn refine_angle_2d(
    face: &Face,
    v: &TopoShape,
    e: &Edge,
    a1: f64,
    a2: f64,
    delta: f64,
) -> Option<f64> {
    let cf = 0.01;
    let tol_int = 1.0e-10;
    let c = make_2d(e, face).ok()?;
    let (t1, t2) = BRepTool::edge_parameters(e);
    // `aTV = BRep_Tool::Parameter(aV, aE, myFace)` (`cxx:1060`) — the vertex
    // parameter on the *pcurve*, not on the 3D edge.
    let tv = vertex_parameter_on_face(v, e);
    let pv = c.d0(tv);
    let t_op = if (tv - t1).abs() < (tv - t2).abs() {
        t2
    } else {
        t1
    };
    // `aGAC1.Load(aC2D, aT1, aT2)`: the intersection is restricted to the
    // CurveOnSurface range, which for a bounded p-curve is its own domain.
    let max_dt = 0.3 * (t2 - t1);
    for i in 0..2 {
        let ai = if i == 0 { a1 } else { a2 + PI };
        let dir = GpDir2d::from_vec2d(&GpVec2d::new(ai.cos(), ai.sin())).ok()?;
        let line = Geom2dLine::from_pnt_dir(pv, dir);
        let points = occt_geom2d::geom2d_api::intersect_curves(c.as_ref(), &line, tol_int);
        let mut t1max = tv;
        let mut t2max = -1.0_f64;
        for p in points {
            if p.u1 < t1 - tol_int || p.u1 > t2 + tol_int {
                continue;
            }
            if p.u2 > t2max && (p.u1 - tv).abs() < max_dt {
                t2max = p.u2;
                t1max = p.u1;
            }
        }
        if t2max <= 0.0 {
            continue;
        }
        let d_t = t_op - t1max;
        if d_t.abs() < tol_int {
            continue;
        }
        let t = t1max + cf * d_t;
        let p_val = c.d0(t);
        let vec2 = GpVec2d::new(p_val.x() - pv.x(), p_val.y() - pv.y());
        let Ok(dir2) = GpDir2d::from_vec2d(&vec2) else {
            continue;
        };
        let angle = dir2d_angle(&dir2);
        if clockwise_angle(a2, angle) < delta {
            return Some(angle);
        }
    }
    None
}

fn edge_closed_on_face(e: &Edge, _face: &Face, face_key: usize) -> bool {
    if BRepTool::is_degenerated(e) {
        return true;
    }
    GeometryRegistry::global().edge_pcurves(&e.0, face_key).len() >= 2
}

fn stored_vertices(e: &Edge) -> Vec<TopoShape> {
    e.0.tshape
        .read()
        .expect("poisoned TShape lock")
        .children
        .iter()
        .filter(|s| s.shape_type() == ShapeType::Vertex)
        .cloned()
        .collect()
}

/// The edge's stored vertices with the `(FORWARD, REVERSED)` orientation pair
/// OCCT guarantees (`BRep_Builder::Add`); the port also builds edges whose two
/// children are stored Forward, so the last child is normalised to Reversed
/// there. No composition with the edge's own orientation.
fn stored_vertices_oriented(e: &Edge) -> Vec<TopoShape> {
    let stored = stored_vertices(e);
    let mut oris: Vec<Orientation> = stored.iter().map(|v| v.orientation()).collect();
    if oris.len() >= 2
        && oris[0] == Orientation::Forward
        && oris[1] == Orientation::Forward
    {
        oris[1] = Orientation::Reversed;
    }
    stored
        .into_iter()
        .zip(oris)
        .map(|(mut v, o)| {
            v.set_orientation(o);
            v
        })
        .collect()
}

/// Edge vertices with `TopoDS_Iterator(cumOri=true)` orientation.
///
/// When both children are stored Forward (`make_edge_segment` in this port),
/// the last child is treated as Reversed, matching `BRepLib_MakeEdge`.
fn oriented_vertices(e: &Edge) -> Vec<TopoShape> {
    stored_vertices_oriented(e)
        .into_iter()
        .map(|mut v| {
            v.set_orientation(Orientation::compose(e.0.orientation(), v.orientation()));
            v
        })
        .collect()
}
