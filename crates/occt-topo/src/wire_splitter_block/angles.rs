use super::prelude::*;
use super::EdgeInfo;
use super::vertices::*;

/// `ClockWiseAngle`.
pub(super) fn clockwise_angle(angle_in: f64, angle_out: f64) -> f64 {
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
pub(super) fn dir2d_angle(d: &GpDir2d) -> f64 {
    let a = d.y().atan2(d.x());
    if a < 0.0 {
        a + TWO_PI
    } else {
        a
    }
}

pub(super) fn coord2d(v: &TopoShape, e: &Edge, face: &Face) -> GpPnt2d {
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
pub(super) fn vertex_parameter_on_face(v: &TopoShape, e: &Edge) -> f64 {
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

pub(super) fn coord2d_vf(e: &Edge, face: &Face) -> GpPnt2d {
    for v in oriented_vertices(e) {
        if v.orientation() == Orientation::Forward {
            return coord2d(&v, e, face);
        }
    }
    GpPnt2d::new(99.0, 99.0)
}

pub(super) fn vertex_parameter(v: &TopoShape, e: &Edge) -> f64 {
    let (first, last) = BRepTool::edge_parameters(e);
    let kids = stored_vertices(e);
    if kids.first().is_some_and(|k| is_same(k, v)) {
        first
    } else {
        last
    }
}

pub(super) fn next_vertex(v: &TopoShape, e: &Edge) -> TopoShape {
    for vx in oriented_vertices(e) {
        if !is_equal(&vx, v) {
            return vx;
        }
    }
    v.clone()
}

pub(super) fn nb_ways_out(le: &[EdgeInfo]) -> usize {
    le.iter().filter(|ei| !ei.is_in && !ei.passed).count()
}

pub(super) fn angle_in(e_in: &Edge, le: &[EdgeInfo]) -> f64 {
    for ei in le {
        if ei.is_in && is_equal(&ei.edge.0, &e_in.0) {
            return ei.angle;
        }
    }
    0.0
}

pub(super) fn angle_2d(v: &TopoShape, e: &Edge, face: &Face, is_in: bool) -> f64 {
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
pub(super) fn tolerance_2d(v: &Vertex, face: &Face) -> f64 {
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
pub(super) fn uv_tolerance_2d(v: &Vertex, face: &Face) -> (f64, f64) {
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

