//! The static geometry predicates of `BRepFill_Evolved.cxx:98-116`, `:136-146`,
//! `:172-186`, `:219-302`, `:626-651`, `:1524-1541`, `:1822-1840`,
//! `:2949-3010`, `:3046-3097`, `:3104-3119`.
//!
//! These decide which of the three elementary constructions (`Elementary`,
//! `Planar`, `Vertical`) a piece of the profile needs (`Side`), which end of a
//! parallel edge points the same way as the bisector edge (`IsInversed`), and
//! whether a parallel edge may be closed (`ConcaveSide`).

use std::sync::Arc;

use occt_core::gp::{GpDir, GpLin2d, GpPln, GpPnt};
use occt_geom::curve::Curve;
use occt_core::kernel::geomabs::CurveType;
use occt_geom2d::adaptor_curve::Geom2dAdaptorCurve;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d;
use crate::brep_lib_make_wire::MakeWire;
use crate::brep_tool::BRepTool;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};
use crate::topo_tools::{edge_vertices, edges_of, is_same};

/// `BRepFill_Confusion()` (`cxx:99-103`): the `1.e-6` tolerance the whole
/// class works with.
pub const FILL_CONFUSION: f64 = 1.0e-6;

/// `EdgeVertices(E, V1, V2)` (`cxx:172-186`): `TopExp::Vertices` followed by
/// the orientation swap of a `TopAbs_REVERSED` edge.
///
/// The reference has no null check (`TopoDS::Vertex` of an edge without two
/// vertices raises), so a missing vertex answers the null handle as `None`.
pub fn edge_vertices_oriented(e: &Edge) -> (Option<Vertex>, Option<Vertex>) {
    let (v1, v2) = edge_vertices(e);
    if e.0.orientation() == Orientation::Reversed {
        (v2, v1)
    } else {
        (v1, v2)
    }
}

/// `IsVertical(E)` (`cxx:219-237`): the two endpoint points share `Y` and the
/// edge's 3D curve is a `Geom_Line`.
pub fn is_vertical(e: &Edge) -> bool {
    let (Some(v1), Some(v2)) = edge_vertices(e) else {
        return false;
    };
    let p1 = BRepTool::vertex_point(&v1);
    let p2 = BRepTool::vertex_point(&v2);
    if (p1.y() - p2.y()).abs() < FILL_CONFUSION {
        // `GC->DynamicType() == STANDARD_TYPE(Geom_Line)` (`cxx:232`).
        if let Some(curve) = BRepTool::edge_curve(e) {
            return curve.gp_line().is_some();
        }
    }
    false
}

/// `IsPlanar(E)` (`cxx:242-260`): the two endpoint points share `Z` and the
/// edge's 3D curve is a `Geom_Line`.
pub fn is_planar(e: &Edge) -> bool {
    let (Some(v1), Some(v2)) = edge_vertices(e) else {
        return false;
    };
    let p1 = BRepTool::vertex_point(&v1);
    let p2 = BRepTool::vertex_point(&v2);
    if (p1.z() - p2.z()).abs() < FILL_CONFUSION {
        if let Some(curve) = BRepTool::edge_curve(e) {
            return curve.gp_line().is_some();
        }
    }
    false
}

/// `Side(Profil, Tol)` (`cxx:276-303`): only the first edge of the wire is
/// tested.
///
/// 1 / 2 / 3 = left (`MAT_Left`), left + planar, left + vertical;
/// 4 / 5 / 6 = right, right + planar, right + vertical.
pub fn side(profil: &Wire, tol: f64) -> i32 {
    let Some(first) = edges_of(&profil.0).into_iter().next() else {
        // `TopExp_Explorer` on a wire without edges answers a null
        // `Current()`; `TopoDS::Edge` on it raises `Standard_NullObject`.
        return 1;
    };
    let (v1, v2) = edge_vertices(&first);
    let (Some(v1), Some(v2)) = (v1, v2) else {
        return 1;
    };
    let p1 = BRepTool::vertex_point(&v1);
    let p2 = BRepTool::vertex_point(&v2);

    let mut the_side = if p1.y() < -tol || p2.y() < -tol { 4 } else { 1 };
    if is_vertical(&first) {
        the_side += 2;
    } else if is_planar(&first) {
        the_side += 1;
    }
    the_side
}

/// `DistanceToOZ(V)` (`cxx:3087-3091`).
pub fn distance_to_oz(v: &Vertex) -> f64 {
    BRepTool::vertex_point(v).y().abs()
}

/// `Altitud(V)` (`cxx:3095-3099`).
pub fn altitud(v: &Vertex) -> f64 {
    BRepTool::vertex_point(v).z()
}

/// `PosOnFace(d1, d2, d3)` (`cxx:3017-3043`).
///
/// 0 outside, 1 on the `d2` edge, 2 inside, 3 on the `d3` edge.
pub fn pos_on_face(d1: f64, d2: f64, d3: f64) -> i32 {
    if (d1 - d2).abs() <= FILL_CONFUSION {
        return 1;
    }
    if (d1 - d3).abs() <= FILL_CONFUSION {
        return 3;
    }
    if d2 < d3 {
        if d1 > d2 + FILL_CONFUSION && d1 < d3 - FILL_CONFUSION {
            return 2;
        }
    } else if d1 > d3 + FILL_CONFUSION && d1 < d2 - FILL_CONFUSION {
        return 2;
    }
    0
}

/// `DoubleOrNotInFace(EC, V)` (`cxx:3052-3084`): true when `V` appears zero or
/// two times among the endpoints of the control edges.
pub fn double_or_not_in_face(ec: &[TopoShape], v: &Vertex) -> bool {
    let mut vu = false;
    for candidate in ec {
        let Some(e) = Edge::wrap(candidate.clone()) else {
            continue;
        };
        let (v1, v2) = edge_vertices(&e);
        for vertex in [v1, v2].into_iter().flatten() {
            if is_same(&vertex.0, &v.0) {
                if vu {
                    return true;
                }
                vu = true;
            }
        }
    }
    !vu
}

/// `Bubble(Seq)` (`cxx:1527-1542`): the growth sort of the whole sequence.
pub fn bubble(seq: &mut [f64]) {
    let nb = seq.len();
    if nb < 2 {
        return;
    }
    let mut invert = true;
    while invert {
        invert = false;
        for i in 0..nb - 1 {
            if seq[i + 1] < seq[i] {
                seq.swap(i, i + 1);
                invert = true;
            }
        }
    }
}

/// `Compare(E1, E2)` (`cxx:1825-1840`): `REVERSED` when the first vertex of
/// `E1` is closer to the last vertex of `E2` than to the first one.
pub fn compare(e1: &Edge, e2: &Edge) -> Orientation {
    let (v11, _) = edge_vertices(e1);
    let (v21, v22) = edge_vertices(e2);
    let (Some(v11), Some(v21), Some(v22)) = (v11, v21, v22) else {
        return Orientation::Forward;
    };
    let p1 = BRepTool::vertex_point(&v11);
    let p2 = BRepTool::vertex_point(&v21);
    let p3 = BRepTool::vertex_point(&v22);
    if p1.distance(&p3) < p1.distance(&p2) {
        Orientation::Reversed
    } else {
        Orientation::Forward
    }
}

/// `DS` of `IsInversed` (`cxx:560-572`): the tangent of `S` at the end its
/// orientation points at.
fn end_tangent(shape: &TopoShape) -> Option<occt_core::gp::GpVec> {
    let e = Edge::wrap(shape.clone())?;
    let curve = BRepTool::edge_curve(&e)?;
    let (first, last) = BRepTool::edge_parameters(&e);
    let u = if e.0.orientation() == Orientation::Forward {
        first
    } else {
        last
    };
    let (_, mut d) = curve.d1(u);
    if e.0.orientation() != Orientation::Forward {
        d.reverse();
    }
    Some(d)
}

/// `IsInversed(S, E1, E2, Inverse)` (`cxx:547-610`).
///
/// `Inverse[i]` is set for both slots even when `S` is not an edge: the
/// reference clears them first and returns.
pub fn is_inversed(s: &TopoShape, e1: &Edge, e2: &Edge, inverse: &mut [bool; 2]) {
    inverse[0] = false;
    inverse[1] = false;
    if s.shape_type() != ShapeType::Edge {
        return;
    }
    let Some(ds) = end_tangent(s) else {
        return;
    };

    if BRepTool::is_degenerated(e1) {
        // `cxx:587-590`.
        inverse[0] = true;
    } else if let Some(dc1) = end_tangent(&e1.0) {
        inverse[0] = ds.dot(&dc1) < 0.0;
    }

    if BRepTool::is_degenerated(e2) {
        // `cxx:606-609`.
        inverse[1] = true;
    } else if let Some(dc2) = end_tangent(&e2.0) {
        inverse[1] = ds.dot(&dc2) < 0.0;
    }
}

/// `ConcaveSide(S, F)` (`cxx:628-652`): true when the 2D circle of `S` on `F`
/// is direct, flipped by a `TopAbs_REVERSED` orientation.
///
/// The reference is documented as "not finished, done only for circles": every
/// non-circle answers `false`.
pub fn concave_side(s: &TopoShape, f: &Face) -> bool {
    if s.shape_type() == ShapeType::Vertex {
        return false;
    }
    if s.shape_type() != ShapeType::Edge {
        return false;
    }
    let Some(edge) = Edge::wrap(s.clone()) else {
        return false;
    };
    let Some((curve, first, last)) = boptools_2d::curve_on_surface_range(&edge, f) else {
        return false;
    };
    let adaptor = Geom2dAdaptorCurve::with_range(curve, first, last);
    if adaptor.get_type() != Some(CurveType::Circle) {
        return false;
    }
    let Some(circle) = adaptor.circle() else {
        return false;
    };
    let mut direct = circle.is_direct();
    if s.orientation() == Orientation::Reversed {
        direct = !direct;
    }
    direct
}

/// `Relative(W1, W2, V, Commun)` (`cxx:2952-3002`): whether the two wires
/// share `V` (`Commun`) and which orientation the shared part has.
pub fn relative(w1: &Wire, w2: &Wire, v: &Vertex) -> (Orientation, bool) {
    let find_edge = |w: &Wire| -> Option<Edge> {
        for e in edges_of(&w.0) {
            let (v1, v2) = edge_vertices(&e);
            for vertex in [v1, v2].into_iter().flatten() {
                if is_same(&vertex.0, &v.0) {
                    return Some(e);
                }
            }
        }
        None
    };
    let (Some(e1), Some(e2)) = (find_edge(w1), find_edge(w2)) else {
        return (Orientation::Forward, false);
    };
    let make_wire = |e: &Edge| -> Wire {
        let mut mw = MakeWire::new();
        let _ = mw.add(e);
        mw.wire()
    };
    let ww1 = make_wire(&e1);
    let ww2 = make_wire(&e2);
    let tol = FILL_CONFUSION;
    if side(&ww1, tol) < 4 && side(&ww2, tol) < 4 {
        return (Orientation::Forward, true);
    }
    if side(&ww1, tol) > 4 && side(&ww2, tol) > 4 {
        return (Orientation::Forward, true);
    }
    (Orientation::Reversed, true)
}

/// `SimpleExpression(B, Bis)` (`cxx:3105-3119`).
///
/// The reference unwraps a `Geom2d_TrimmedCurve` over a `Bisector_BisecAna`
/// down to the analytic curve of the bisector, re-trimmed to the same
/// parameters. The port's [`BisectorBisec::value`] already is that trimmed
/// curve and the port's `Geom2dTrimmedCurve` evaluates its basis at the same
/// parameters (`occt-geom2d/src/trimmed.rs`), so the unwrapping is the
/// identity and only the `Value()` step is kept.
///
/// [`BisectorBisec::value`]: occt_geom2d::bisector::BisectorBisec::value
pub fn simple_expression(bis: &Geom2dTrimmedCurve) -> Arc<dyn Curve2d> {
    Arc::new(bis.clone())
}

/// The 2D line `gp::OY2d()` of `PrepareProfile` (`cxx:1566`): the `Y` axis of
/// the `xOy` plane.
pub fn oy2d_line() -> GpLin2d {
    GpLin2d::from_pnt_dir(
        occt_core::gp::GpPnt2d::new(0.0, 0.0),
        occt_core::gp::GpDir2d::new(0.0, 1.0).expect("gp::OY2d"),
    )
}

/// `gp_Pln(P, V)` for the cap faces of `AddTopAndBottom`: the plane through
/// `P = (0, 0, z)` with normal `+Z`. `gp_Pln(const gp_Pnt&, const gp_Dir&)`
/// keeps the direction of `gp_Ax3(P, V)` as `+X` whenever `V` is not `+Z`; for
/// `V = +Z` the X direction is arbitrary and the port fixes `gp::DX()`.
pub fn plane_at_z(z: f64) -> GpPln {
    let ax3 = occt_core::gp::GpAx3::new(
        GpPnt::new(0.0, 0.0, z),
        GpDir::new(0.0, 0.0, 1.0).expect("gp::DZ"),
        &GpDir::new(1.0, 0.0, 0.0).expect("gp::DX"),
    )
    .expect("gp_Ax3: DX is normal to DZ");
    GpPln::new(ax3)
}

/// `gp_Ax3(gp::YOZ())` (`cxx:1565`): the `yOz` plane, normal `+X`, `X` along
/// `+Y`.
pub fn yoz_plane() -> GpPln {
    GpPln::new(
        occt_core::gp::GpAx3::new(
            GpPnt::new(0.0, 0.0, 0.0),
            occt_core::gp::GpDir::new(1.0, 0.0, 0.0).expect("gp::DX"),
            &occt_core::gp::GpDir::new(0.0, 1.0, 0.0).expect("gp::DY"),
        )
        .expect("gp::YOZ"),
    )
}

/// `TopAbs::Complement` (`TopAbs.cxx`): `FORWARD` <-> `REVERSED`.
pub fn complement(o: Orientation) -> Orientation {
    match o {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        other => other,
    }
}

/// The first edge of `wire` as `BRepTools_WireExplorer` sees it: the explorer
/// reorders the edges into a connected sequence starting at the wire's first
/// vertex (`BRepTools_WireExplorer.cxx:60-140`).
pub fn wire_first_edge(wire: &Wire) -> Option<Edge> {
    let exp = WireExplorer::of_wire(wire);
    exp.current_edge()
}

/// `Curve` of an edge as `BRepAdaptor_Curve` sees it: the stored curve, or the
/// pcurve lifted to 3D when the edge has no 3D representation.
pub fn adaptor_curve(e: &Edge) -> Option<Arc<dyn Curve>> {
    BRepTool::edge_curve_adaptor_world(e)
}
