//! The file-local statics of `BRepOffsetAPI_MiddlePath`.
//!
//! Source: `ModelingAlgorithms/TKOffset/BRepOffsetAPI/BRepOffsetAPI_MiddlePath.cxx`
//! (OCCT 8.0.0):
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `IsLinear(anEdge, aLine)` | `cxx:60-99` | [`is_linear`] |
//! | `TypeOfEdge(anEdge)` | `cxx:101-111` | [`type_of_edge`] |
//! | `TangentOfEdge(aShape, OnFirst)` | `cxx:113-139` | [`tangent_of_edge`] |
//! | `IsValidEdge(theEdge, theFace)` | `cxx:141-177` | [`is_valid_edge`] |
//! | `GetUnifiedWire(theWire, theUnifier)` | `cxx:266-293` | [`get_unified_wire`] |
//!
//! `IsLinear` is the `Geom_Line` / 2-pole `Bezier` / 2-pole `BSpline` test; the
//! `gce_MakeLin(P1, P2)` of the two pole arms (`cxx:83`, `:94`) raises
//! `gce_ConfusedPoints` for coincident poles and the port answers `None`.
//! `TypeOfEdge` falls back to `BRepAdaptor_Curve::GetType()` (`cxx:109-110`),
//! i.e. the adaptor's `GeomAbs_CurveType`, which is
//! [`curve_kind`](crate::geom_abs_curve_type::curve_kind) here; an edge with no
//! 3D curve raises in the reference and answers `GeomAbs_OtherCurve` here.
//! `IsValidEdge` runs `BRepExtrema_DistShapeShape(theEdge, anEdge)` against
//! every edge of the face with the default flag / algorithm
//! (`Extrema_ExtFlag_MIN`, `Extrema_ExtAlgo_Grad`) and rejects a solution that
//! lands *on* the other edge or on a vertex that is neither end of `theEdge`.
//! `GetUnifiedWire` replays the `ShapeUpgrade_UnifySameDomain` history.

use std::sync::Arc;

use occt_core::gp::{GpDir, GpLin, GpPnt, GpVec};
use occt_core::precision::CONFUSION;
use occt_geom::extrema_surf::{ExtremaExtAlgo, ExtremaExtFlag};
use occt_geom::{Curve, Surface};
use occt_geom2d::Curve2d;

use crate::abs::{Orientation, ShapeType};
use crate::bop_hist::BopHistory;
use crate::brep_extrema_dist_shape_shape::DistShapeShape;
use crate::brep_extrema_solution_elem::SupportType;
use crate::brep_lib_build_curves3d::build_curve3d;
use crate::brep_lib_make_edge;
use crate::brep_lib_make_wire::MakeWire;
use crate::brep_tool::BRepTool;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::builder::TopoBuilder;
use crate::geom_abs_curve_type::{curve_kind, GeomAbsCurveType};
use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{edge_vertices, is_same, shape_map_key};

/// `static bool IsLinear(const TopoDS_Edge& anEdge, gp_Lin& aLine)`
/// (`cxx:60-99`).
pub(super) fn is_linear(an_edge: &Edge) -> Option<GpLin> {
    // `cxx:64-68`: `BRep_Tool::Curve` then `IsInstance(Geom_TrimmedCurve)`.
    let mut a_curve = BRepTool::edge_curve(an_edge)?;
    if let Some((basis, _, _)) = a_curve.untrimmed_basis() {
        a_curve = basis;
    }
    // `cxx:70-75`.
    if let Some(lin) = a_curve.gp_line() {
        return Some(lin);
    }
    // `cxx:76-86`: a 2-pole Bezier.
    if let Some(poles) = a_curve.bezier_poles() {
        if poles.len() == 2 {
            return line_through(&poles[0], &poles[1]);
        }
    }
    // `cxx:87-95`: a 2-pole BSpline.
    if let Some(poles) = a_curve.bspline_poles() {
        if poles.len() == 2 {
            return line_through(&poles[0], &poles[1]);
        }
    }
    // `cxx:97`.
    None
}

/// `gce_MakeLin(Pnt1, Pnt2)` (`cxx:83`, `:94`): `None` for coincident points
/// (`gce_ConfusedPoints`).
fn line_through(p1: &GpPnt, p2: &GpPnt) -> Option<GpLin> {
    if p1.distance(p2) <= f64::EPSILON {
        return None;
    }
    let dir = GpDir::from_vec(&GpVec::from_pnts(p1, p2)).ok()?;
    Some(GpLin::from_pnt_dir(*p1, dir))
}

/// `static GeomAbs_CurveType TypeOfEdge(const TopoDS_Edge& anEdge)`
/// (`cxx:101-111`).
pub(super) fn type_of_edge(an_edge: &Edge) -> GeomAbsCurveType {
    // `cxx:105-108`.
    if is_linear(an_edge).is_some() {
        return GeomAbsCurveType::Line;
    }
    // `cxx:109-110`: `BRepAdaptor_Curve(anEdge).GetType()`.
    match BRepTool::edge_curve(an_edge) {
        Some(c) => curve_kind(c.as_ref()),
        None => GeomAbsCurveType::OtherCurve,
    }
}

/// `static gp_Vec TangentOfEdge(const TopoDS_Shape& aShape, const bool
/// OnFirst)` (`cxx:113-139`).
pub(super) fn tangent_of_edge(a_shape: &TopoShape, on_first: bool) -> GpVec {
    // `cxx:115-116`: `TopoDS::Edge(aShape)`.
    let Some(an_edge) = Edge::wrap(a_shape.clone()) else {
        return GpVec::zero();
    };
    let an_or = a_shape.orientation();
    // `cxx:119-128`.
    let (fpar, lpar) = BRepTool::edge_parameters(&an_edge);
    let the_par = if on_first {
        if an_or == Orientation::Forward {
            fpar
        } else {
            lpar
        }
    } else if an_or == Orientation::Forward {
        lpar
    } else {
        fpar
    };
    // `cxx:130-137`.
    let Some(curve) = BRepTool::edge_curve(&an_edge) else {
        return GpVec::zero();
    };
    let (_, mut tangent) = curve.d1(the_par);
    if an_or == Orientation::Reversed {
        tangent.reverse();
    }
    tangent
}

/// `static bool IsValidEdge(const TopoDS_Edge& theEdge, const TopoDS_Face&
/// theFace)` (`cxx:141-177`).
pub(super) fn is_valid_edge(the_edge: &Edge, the_face: &Face) -> bool {
    // `cxx:144-145`: `TopExp::Vertices(theEdge, V1, V2)`, i.e. CumOri = false.
    let (v1, v2) = edge_vertices(the_edge);
    // `cxx:147`.
    let tol = CONFUSION;
    // `cxx:149-150`: `TopExp_Explorer(theFace, TopAbs_EDGE)`.
    for an_edge in crate::topo_tools::explorer_shapes(&the_face.0, ShapeType::Edge) {
        let an_edge = Edge(an_edge);
        // `cxx:154-155`: the default flag / algorithm of the constructor.
        let dist = DistShapeShape::new(
            &the_edge.0,
            &an_edge.0,
            CONFUSION,
            ExtremaExtFlag::MinMax,
            ExtremaExtAlgo::Grad,
        );
        if dist.value() > tol {
            continue;
        }
        // `cxx:157-170`.
        for i in 1..=dist.nb_solution() {
            if dist.support_type_shape2(i) == SupportType::IsOnEdge {
                return false;
            }
            // theType is "IsVertex".
            let a_vertex = dist.support_on_shape2(i);
            let Some(a_vertex) = a_vertex else {
                return false;
            };
            let same_v1 = v1.as_ref().map_or(false, |v| is_same(&a_vertex, &v.0));
            let same_v2 = v2.as_ref().map_or(false, |v| is_same(&a_vertex, &v.0));
            if !(same_v1 || same_v2) {
                return false;
            }
        }
    }
    // `cxx:175`.
    true
}

/// `static TopoDS_Wire GetUnifiedWire(const TopoDS_Wire& theWire,
/// ShapeUpgrade_UnifySameDomain& theUnifier)` (`cxx:266-293`).
///
pub(super) fn get_unified_wire(the_wire: &Wire, history: &BopHistory) -> Result<Wire, String> {
    let mut a_w_maker = MakeWire::new();
    let mut generated_edges = std::collections::HashSet::new();
    let mut wexp = WireExplorer::of_wire(the_wire);
    // `cxx:269-292`.
    while wexp.more() {
        let Some(an_edge) = wexp.current_edge() else {
            break;
        };
        if let Some(modified) = history.modified(&an_edge.0) {
            for shape in modified {
                if !generated_edges.insert(shape_map_key(shape)) {
                    continue;
                }
                let edge = Edge::wrap(shape.clone())
                    .ok_or_else(|| "Standard_TypeMismatch: TopoDS::Edge".to_string())?;
                a_w_maker
                    .add(&edge)
                    .map_err(|e| format!("BRepLib_MakeWire: {e:?}"))?;
            }
        } else {
            a_w_maker
                .add(&an_edge)
                .map_err(|e| format!("BRepLib_MakeWire: {e:?}"))?;
        }
        wexp.next();
    }
    // `cxx:293`.
    Ok(a_w_maker.wire())
}

/// `BRepLib::BuildCurve3d(NewEdge)` (`cxx:698`, `:737`, `:745`, `:810`) with
/// the `BRepLib.hxx` defaults (`1.0e-5`, `GeomAbs_C1`, `14`, `0`).
pub(super) fn build_curve3d_default(edge: &Edge) {
    let _ = build_curve3d(edge, 1.0e-5, occt_core::kernel::geomabs::Shape::C1, 14, 0);
}

/// `GC_MakeLine2d(P1, P2)` (`cxx:692`, `:732`, `:738`, `:739`): the 2D line
/// through two points as a `Geom2d_Line`. `None` is `gce_ConfusedPoints`
/// (`gce_MakeLin2d`).
pub(super) fn make_line_2d(
    p1: &occt_core::gp::GpPnt2d,
    p2: &occt_core::gp::GpPnt2d,
) -> Option<Arc<dyn Curve2d>> {
    let v = occt_core::gp::GpVec2d::new(p2.x() - p1.x(), p2.y() - p1.y());
    let dir = occt_core::gp::GpDir2d::from_vec2d(&v).ok()?;
    let line = occt_geom2d::line::Geom2dLine::from_pnt_dir(*p1, dir);
    Some(Arc::new(line) as Arc<dyn Curve2d>)
}

/// `BRep_Tool::CurveOnSurface(E, F, f, l)` (`BRep_Tool.cxx:310-360`): the
/// stored p-curve of `e` on `f`, the surface handle and the p-curve range.
pub(super) fn pcurve_on_face(
    e: &Edge,
    f: &Face,
) -> Option<(Arc<dyn Curve2d>, Arc<dyn Surface>, f64, f64)> {
    let s = BRepTool::face_surface(f)?;
    let key = Arc::as_ptr(&s) as *const () as usize;
    let reps: Vec<(Arc<dyn Curve2d>, f64, f64)> = GeometryRegistry::global()
        .edge_pcurve_reps(&e.0)
        .into_iter()
        .filter(|(surf, _, _, _)| Arc::as_ptr(surf) as *const () as usize == key)
        .map(|(_, pc, a, b)| (pc, a, b))
        .collect();
    // `BRep_Tool::CurveOnSurface` asks for `PCurve2()` when the representation
    // is a closed surface and the edge is REVERSED, `PCurve()` otherwise.
    if reps.len() >= 2 && e.0.orientation() == Orientation::Reversed {
        let (pc, a, b) = reps.into_iter().nth(1)?;
        return Some((pc, s, a, b));
    }
    let (pc, a, b) = reps.into_iter().next()?;
    Some((pc, s, a, b))
}

/// `BRepLib_MakeEdge(C2d, S, V1, V2, f, l)` (`cxx:695`, `:734`, `:742`, `:748`).
pub(super) fn make_edge_2d(
    pc: &Arc<dyn Curve2d>,
    surf: &Arc<dyn Surface>,
    face: &Face,
    v1: Option<&Vertex>,
    v2: Option<&Vertex>,
    p1: f64,
    p2: f64,
) -> Option<Edge> {
    let build = brep_lib_make_edge::init_2d_vertices(pc, surf, v1, v2, p1, p2)?;
    let mut e = Edge::new();
    let reg = GeometryRegistry::global();
    let face_key = GeometryRegistry::shape_key(&face.0);
    reg.set_edge_pcurve(&e.0, face_key, pc.clone());
    reg.set_ranges_all(&e.0, build.first, build.last);
    reg.raise_edge_tolerance(&e.0, CONFUSION);
    brep_lib_make_edge::add_endpoint_vertices(&TopoBuilder::new(), &mut e, &build);
    if build.reverse {
        e.0.set_orientation(Orientation::Reversed);
    }
    Some(e)
}

/// `BRepLib_MakeEdge(const Handle(Geom_Curve)& C, const Standard_Real U1,
/// const Standard_Real U2)` (`cxx:1012`): `None` is the `!IsDone()` arm of
/// `BRepLib_MakeEdge::Build` (`myError != BRepLib_EdgeDone`), which the
/// reference turns into a null `MidEdges(i)` (`cxx:1014`).
pub(super) fn make_edge_range(curve: &Arc<dyn Curve>, first: f64, last: f64) -> Option<Edge> {
    let build = brep_lib_make_edge::init_3d(curve, first, last)?;
    let b = TopoBuilder::new();
    let mut e = b.make_edge(curve.clone(), first, last);
    brep_lib_make_edge::add_endpoint_vertices(&b, &mut e, &build);
    if build.closed {
        e.0.set_closed(true);
    }
    Some(e)
}

/// `BRepLib_MakeEdge(const Handle(Geom_Curve)& C)` (`cxx:1100`): the whole
/// curve's parameter range.
pub(super) fn make_edge_over_curve(curve: &Arc<dyn Curve>) -> Option<Edge> {
    make_edge_range(curve, curve.first_parameter(), curve.last_parameter())
}
