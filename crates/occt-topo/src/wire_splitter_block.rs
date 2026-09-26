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
use occt_geom2d::curve::Curve2d;
use occt_geom2d::line::Geom2dLine;

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d::make_2d;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::pcurve::{pc_curve_kind, CurveKind};
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;

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
pub(crate) fn split_block(face: &Face, edges: &[Edge]) -> Result<Vec<TopoShape>, String> {
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
            let cyc = closed || (kids.len() >= 2 && kids[0].same_tshape(&kids[1]));
            let e0 = vert_closed.entry(vk).or_insert(false);
            if cyc || i > 0 && kids[0].same_tshape(v) {
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
        if ls.len() == 1 && ls[0].0.same_tshape(&map.infos[va_slot][info_idx].edge.0) {
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

        let vb = get_next_vertex(&p_va, &e_out);
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
            let a_nb = ls.len();
            let mut cut: Option<usize> = None;
            for i in (0..a_nb).rev() {
                let e_prev = &ls[i];
                buf.push(e_prev.clone());
                if !has_edge {
                    has_edge = !BRepTool::is_degenerated(e_prev);
                    if !has_edge {
                        continue;
                    }
                }
                let same_v = vert_va[i].same_tshape(&vb);
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
                    if buf.len() == 2 && buf[0].0.same_tshape(&buf[1].0) {
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
        let a_le = &map.infos[vb_slot];
        let angle_in = angle_in(e_in.as_ref().unwrap_or(&e_out), a_le);
        let i_cnt = nb_ways_out(a_le);
        let mut min_angle = 100.0;
        let mut chosen: Option<usize> = None;
        let mut only_inside: Option<usize> = None;
        let mut n_inside = 0usize;

        for (k, ei) in a_le.iter().enumerate() {
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
            let ang = if ei.edge.0.same_tshape(&e_out.0) {
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
    let mut a_in = angle_in;
    let mut a_out = angle_out;
    if a_in >= TWO_PI {
        a_in -= TWO_PI;
    }
    if a_out >= TWO_PI {
        a_out -= TWO_PI;
    }
    let mut a1 = a_in + PI;
    if a1 >= TWO_PI {
        a1 -= TWO_PI;
    }
    let mut da = a1 - a_out;
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
        if !k.same_tshape(v) {
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
    if kids.first().is_some_and(|k| k.same_tshape(v)) {
        first
    } else {
        last
    }
}

fn get_next_vertex(v: &TopoShape, e: &Edge) -> TopoShape {
    for vx in oriented_vertices(e) {
        if !(vx.same_tshape(v) && vx.orientation() == v.orientation()) {
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
        if ei.is_in
            && ei.edge.0.same_tshape(&e_in.0)
            && ei.edge.0.orientation() == e_in.0.orientation()
        {
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
    let mut dt = curve_resolution(c.as_ref(), tol2d).max(PCONFUSION);
    if pc_curve_kind(c.as_ref()) != CurveKind::Line {
        let r = curve2d_radius(c.as_ref(), tv);
        if r > PCONFUSION {
            let cosphi = r / (r + tol2d);
            if cosphi.abs() <= 1.0 {
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

fn curve_resolution(c: &dyn Curve2d, tol2d: f64) -> f64 {
    let mid = 0.5 * (c.first_parameter() + c.last_parameter());
    let t = if mid.is_finite() {
        mid
    } else {
        0.0
    };
    let (_, d1) = c.d1(t);
    let mag = d1.magnitude();
    if mag > PCONFUSION {
        tol2d / mag
    } else {
        tol2d
    }
}

fn curve2d_radius(c: &dyn Curve2d, t: f64) -> f64 {
    let (_, d1, d2) = c.d2(t);
    let n2 = d1.square_magnitude();
    if n2 <= PCONFUSION {
        return 0.0;
    }
    let k = d1.crossed(&d2).abs() / n2.powf(1.5);
    if k > PCONFUSION {
        1.0 / k
    } else {
        0.0
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
            == crate::geom_bnd_lib_surface3d::SurfaceKind::BSplineSurface
        {
            // `cxx:875-878`.
            t2 *= 1.1;
        }
    }
    t2
}

/// `BOPAlgo_WireSplitter::UTolerance2D` / `VTolerance2D` (`_1.cxx:885-901`):
/// `BRepAdaptor_Surface::U/VResolution` (`GeomAdaptor_Surface.cxx:1818-1945`).
///
/// UNPORTED: the revolution (V) and extrusion (U) arms need the *basis curve's*
/// `Resolution`; those two keep the former finite-difference estimate.
fn uv_tolerance_2d(v: &Vertex, face: &Face) -> (f64, f64) {
    let t3 = BRepTool::vertex_tolerance(v);
    let Some(surf) = BRepTool::face_surface(face) else {
        return (t3, t3);
    };
    let s = surf.as_ref();
    use crate::geom_bnd_lib_surface3d::{surface_kind, SurfaceKind};
    // `if (Res <= 1.) return 2*asin(Res); return 2*pi;` (`cxx:1913-1916`, `:1943-1945`).
    let angular = |res: f64| {
        if res <= 1.0 {
            2.0 * res.asin()
        } else {
            2.0 * PI
        }
    };
    let exact_uv = || s.uv_resolution(t3);
    match surface_kind(s) {
        // `case GeomAbs_Torus: Res = R3d / (2*(Major+Minor))` (U), `Minor` (V).
        SurfaceKind::Torus => {
            let (major, minor) = s
                .gp_torus()
                .map(|t| (t.major_radius(), t.minor_radius()))
                .unwrap_or((0.0, 0.0));
            let ru = if major + minor > occt_core::precision::CONFUSION {
                t3 / (2.0 * (major + minor))
            } else {
                0.0
            };
            let rv = if minor > occt_core::precision::CONFUSION { t3 / (2.0 * minor) } else { 0.0 };
            (angular(ru), angular(rv))
        }
        // Sphere: both directions `R3d / (2*R)`.
        SurfaceKind::Sphere => {
            let r = s.gp_sphere().map(|x| x.radius()).unwrap_or(0.0);
            let res = if r > occt_core::precision::CONFUSION { t3 / (2.0 * r) } else { 0.0 };
            (angular(res), angular(res))
        }
        // Cylinder: U `R3d/(2*R)`, V `R3d`.
        SurfaceKind::Cylinder => {
            let r = s.gp_cylinder().map(|x| x.radius()).unwrap_or(0.0);
            let res = if r > occt_core::precision::CONFUSION { t3 / (2.0 * r) } else { 0.0 };
            (angular(res), t3)
        }
        // Cone: U uses the largest `VIso` radius over the V range, V is `R3d`.
        SurfaceKind::Cone => {
            let (v0, v1) = s.v_range();
            let radius_at = |vv: f64| s.d0(0.0, vv).distance(&s.d0(PI, vv)) * 0.5;
            let r = if v0.is_finite() && v1.is_finite() {
                radius_at(v0).max(radius_at(v1))
            } else {
                0.0
            };
            let ru = if r > occt_core::precision::CONFUSION { t3 / r } else { 0.0 };
            (ru, t3)
        }
        // Plane / extrusion (V): both `R3d`.
        SurfaceKind::Plane => (t3, t3),
        // Bezier / BSpline / offset: the surface's own `Resolution`.
        SurfaceKind::BezierSurface | SurfaceKind::BSplineSurface | SurfaceKind::OffsetSurface => {
            match exact_uv() {
                Some((ur, vr)) => (ur, vr),
                None => fd_uv_tolerance(v, s, t3),
            }
        }
        // Revolution (V) and extrusion (U) need the basis curve's `Resolution`
        // (UNPORTED — see the doc comment); default arm is
        // `Precision::Parametric(R3d)` = `R3d * 0.01` (`Precision.hxx:328`).
        SurfaceKind::SurfaceOfRevolution | SurfaceKind::SurfaceOfExtrusion => {
            fd_uv_tolerance(v, s, t3)
        }
        _ => (t3 * 0.01, t3 * 0.01),
    }
}

/// Former finite-difference estimate `t3 / |dS/du|` — kept only for the
/// UNPORTED revolution/extrusion arms (and as a fallback when a surface
/// reports no `Resolution`).
fn fd_uv_tolerance(v: &Vertex, s: &dyn occt_geom::surface::Surface, t3: f64) -> (f64, f64) {
    let t3 = t3.max(PCONFUSION);
    let p = BRepTool::vertex_point(v);
    // UNPORTED: `fd_uv_tolerance` stands in for the unported
    // revolution/extrusion `Resolution` arms (see the caller's doc comment);
    // the grid supplies the (u, v) at which the finite-difference derivative is
    // evaluated.
    let (u, vv) = crate::brep_surface::surface_closest_params(s, &p, 8, 8);
    let (_, du, dv) = s.d1(u, vv);
    let ur = if du.magnitude() > PCONFUSION {
        t3 / du.magnitude()
    } else {
        t3
    };
    let vr = if dv.magnitude() > PCONFUSION {
        t3 / dv.magnitude()
    } else {
        t3
    };
    (ur, vr)
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
    let mut a_a1 = 0.0_f64; // angle of the outgoing boundary edge
    let mut a_a2 = 0.0_f64; // angle of the incoming boundary edge
    let mut i_cnt_bnd = 0usize;
    let mut i_cnt_int = 0usize;
    for ei in infos.iter() {
        if !ei.is_inside {
            i_cnt_bnd += 1;
            if ei.is_in {
                a_a2 = ei.angle;
            } else {
                a_a1 = ei.angle;
            }
        } else {
            i_cnt_int += 1;
        }
    }
    if i_cnt_bnd != 2 {
        return;
    }
    let a_delta = clockwise_angle(a_a2, a_a1);
    let mut refined: HashMap<usize, f64> = HashMap::new();
    for ei in infos.iter() {
        if !ei.is_inside || ei.is_in {
            continue;
        }
        if clockwise_angle(a_a2, ei.angle) < a_delta {
            continue; // already inside the wedge
        }
        let ek = GeometryRegistry::shape_key(&ei.edge.0);
        match refine_angle_2d(face, v, &ei.edge, a_a1, a_a2, a_delta) {
            Some(a_new) => {
                refined.insert(ek, a_new);
            }
            None if i_cnt_int == 2 => {
                let a_new = if ei.angle <= a_a1 {
                    a_a1 + ANGULAR
                } else {
                    a_a2 - ANGULAR
                };
                refined.insert(ek, a_new);
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
///   `[a_t1, a_t2]` restriction itself (`:1084-1100`). The two agree when the
///   intersector returns every intersection point of that domain.
/// * **UNPORTED**: OCCT takes the vertex parameter from the *pcurve*
///   (`aTV = BRep_Tool::Parameter(aV, aE, myFace)`, `cxx:1060`); this port uses
///   the **3D edge parameter** ([`vertex_parameter`]). That is only valid when
///   the pcurve and the 3D curve share a parameterization (SameParameter): the
///   degenerate "pcurve parameterization ≠ 3D edge parameterization" case is
///   handled silently as if it were SameParameter. `IntTools_Curve` and
///   `MakeSplitEdge` keep that property for the BOP inputs, but this is not
///   line-by-line equivalent.
fn refine_angle_2d(
    face: &Face,
    v: &TopoShape,
    e: &Edge,
    a_a1: f64,
    a_a2: f64,
    a_delta: f64,
) -> Option<f64> {
    let a_cf = 0.01;
    let a_tol_int = 1.0e-10;
    let c = make_2d(e, face).ok()?;
    let (a_t1, a_t2) = BRepTool::edge_parameters(e);
    let a_tv = vertex_parameter(v, e);
    let a_pv = c.d0(a_tv);
    let a_t_op = if (a_tv - a_t1).abs() < (a_tv - a_t2).abs() {
        a_t2
    } else {
        a_t1
    };
    // `aGAC1.Load(aC2D, aT1, aT2)`: the intersection is restricted to the
    // CurveOnSurface range, which for a bounded p-curve is its own domain.
    let max_dt = 0.3 * (a_t2 - a_t1);
    for i in 0..2 {
        let a_ai = if i == 0 { a_a1 } else { a_a2 + PI };
        let dir = GpDir2d::from_vec2d(&GpVec2d::new(a_ai.cos(), a_ai.sin())).ok()?;
        let line = Geom2dLine::from_pnt_dir(a_pv, dir);
        let points = occt_geom2d::geom2d_api::intersect_curves(c.as_ref(), &line, a_tol_int);
        let mut a_t1max = a_tv;
        let mut a_t2max = -1.0_f64;
        for p in points {
            if p.u1 < a_t1 - a_tol_int || p.u1 > a_t2 + a_tol_int {
                continue;
            }
            if p.u2 > a_t2max && (p.u1 - a_tv).abs() < max_dt {
                a_t2max = p.u2;
                a_t1max = p.u1;
            }
        }
        if a_t2max <= 0.0 {
            continue;
        }
        let d_t = a_t_op - a_t1max;
        if d_t.abs() < a_tol_int {
            continue;
        }
        let a_t = a_t1max + a_cf * d_t;
        let a_p = c.d0(a_t);
        let vec2 = GpVec2d::new(a_p.x() - a_pv.x(), a_p.y() - a_pv.y());
        let Ok(dir2) = GpDir2d::from_vec2d(&vec2) else {
            continue;
        };
        let a_angle = dir2d_angle(&dir2);
        if clockwise_angle(a_a2, a_angle) < a_delta {
            return Some(a_angle);
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
