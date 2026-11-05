//! `BRepAlgo::ConvertWire` (`BRepAlgo.cxx:80-219`) and
//! `BRepAlgo::ConvertFace` (`cxx:223-242`).
//!
//! Every edge's p-curve is approximated by arcs of circles / line segments
//! (`Geom2dConvert_ApproxArcsSegments`), the pieces are lifted back onto the
//! face's surface and re-assembled into a wire of `BRepLib_MakeEdge(C2d, S, V1,
//! V2, p1, p2)` edges. `BRepLib::BuildCurves3d` is **not** called here - OCCT
//! leaves it to the two callers (`BRepOffsetAPI_MakeOffset.cxx:81`, `:370`,
//! `:379`).

use std::sync::Arc;

use occt_core::precision::CONFUSION;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::geom2d_convert_approx_arcs_segments::ApproxArcsSegments;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use crate::abs::{Orientation, ShapeType};
use crate::brep_lib_make_edge;
use crate::brep_lib_make_face;
use crate::brep_lib_make_wire::MakeWire;
use crate::brep_tool::BRepTool;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, Vertex, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::explorer_shapes;

/// `MINIMAL_TOLERANCE` (`cxx:56`): the floor on the tolerance an edge's p-curve
/// approximation runs with.
const MINIMAL_TOLERANCE: f64 = 0.0001;

/// `struct OrientedCurve` (`cxx:61-74`): one piece of the approximation, plus
/// the flag saying whether the edge it came from is stored REVERSED.
struct OrientedCurve {
    /// `Curve` (`cxx:63`): the `Geom2d_TrimmedCurve` `GetResult()` handed back.
    curve: Arc<dyn Curve2d>,
    /// `IsReverse` (`cxx:64`).
    is_reverse: bool,
}

impl OrientedCurve {
    /// `Point(isEnd)` (`cxx:67-72`): `StartPoint()` when `isEnd == IsReverse`,
    /// `EndPoint()` otherwise.
    fn point(&self, is_end: bool) -> occt_core::gp::GpPnt2d {
        let u = if is_end == self.is_reverse {
            self.curve.first_parameter()
        } else {
            self.curve.last_parameter()
        };
        self.curve.d0(u)
    }
}

/// `cxx:914-922` of `BRepLib_MakeEdge::Init(C2d, ...)`: the trim-kill OCCT
/// applies before it stores the p-curve, so the edge's `[p1, p2]` range is in
/// basis parameters (`Geom2d_TrimmedCurve::First/LastParameter` already are).
fn untrimmed_2d(curve: &Arc<dyn Curve2d>) -> Arc<dyn Curve2d> {
    let mut c = curve.clone();
    loop {
        let next = c.trimmed_basis().map(|b| Arc::from(b.clone_dyn()));
        match next {
            Some(n) => c = n,
            None => return c,
        }
    }
}

/// `BRepBuilderAPI_MakeEdge(C2d, S, V1, V2, p1, p2)`
/// (`BRepBuilderAPI_MakeEdge.cxx:483-496`) over
/// `BRepLib_MakeEdge::Init(C2d, S, VV1, VV2, pp1, pp2)` (`cxx:905-1094`).
///
/// The edge carries the p-curve on `face` and **no** 3D curve, exactly as OCCT
/// leaves it: `B.MakeEdge(E)` (`cxx:1075`) + `B.UpdateEdge(E, C, S,
/// TopLoc_Location(), preci)` (`cxx:1077`) only append a CurveOnSurface
/// representation; `BRepLib::BuildCurve3d` materialises the 3D curve later
/// (`BRepLib.cxx:301-455`). `B.Range(E, p1, p2)` (`cxx:1087`) writes both
/// representations' range, and `Add(E, V1)` / `Add(E, V2)` happen before it
/// (`cxx:1079-1086`).
#[allow(clippy::too_many_arguments)]
fn make_edge_2d(
    builder: &TopoBuilder,
    curve: &Arc<dyn Curve2d>,
    surface: &Arc<dyn occt_geom::Surface>,
    face: &Face,
    v1: &Vertex,
    v2: &Vertex,
    p1: f64,
    p2: f64,
) -> Option<Edge> {
    let build = brep_lib_make_edge::init_2d_vertices(curve, surface, Some(v1), Some(v2), p1, p2)?;
    let mut edge = Edge::new();
    let reg = GeometryRegistry::global();
    let face_key = GeometryRegistry::shape_key(&face.0);
    reg.set_edge_pcurve(&edge.0, face_key, untrimmed_2d(curve));
    // `cxx:1087`: `B.Range(E, p1, p2)`.
    reg.set_ranges_all(&edge.0, build.first, build.last);
    // `cxx:1077`: `TE->UpdateTolerance(BRepLib::Precision())`.
    reg.raise_edge_tolerance(&edge.0, CONFUSION);
    brep_lib_make_edge::add_endpoint_vertices(builder, &mut edge, &build);
    if build.reverse {
        // `cxx:1089-1092`.
        edge.0.set_orientation(Orientation::Reversed);
    }
    Some(edge)
}

/// `BRepAlgo::ConvertWire(theWire, theAngleTol, theFace)` (`cxx:80-219`).
///
/// Each edge of `theWire` is walked in connection order by
/// `BRepTools_WireExplorer(W, F)` (`cxx:89`), its p-curve approximated by
/// `Geom2dConvert_ApproxArcsSegments(aCurve, max(Tolerance(E),
/// MINIMAL_TOLERANCE), AngleTol)` (`cxx:91-102`) and the pieces collected in
/// edge order, reversed edge by reversed edge (`cxx:108-133`). The chain is
/// then rebuilt: two consecutive pieces whose ends are closer than
/// `aMaxTol + Precision::Confusion()` share one vertex placed at their midpoint
/// (`cxx:143-149`, `:180-187`), and when the wire's closing point is that close
/// to the first vertex the first vertex is reused (`cxx:167-178`). Every new
/// edge is `BRepBuilderAPI_MakeEdge(C2d, S, V1, V2, p1, p2)`.
///
/// `None` is the null `TopoDS_Wire` of `cxx:134` / `:212-215`. The single
/// divergence is a face without a surface: OCCT dereferences the null handle at
/// `cxx:103` (undefined behaviour), the port answers `None`.
pub fn convert_wire(w: &Wire, angle_tol: f64, face: &Face) -> Option<Wire> {
    use occt_core::gp::GpPnt2d;

    // `cxx:86`: `const Handle(Geom_Surface) aSurf = BRep_Tool::Surface(theFace)`.
    let surface = BRepTool::face_surface(face)?;
    let mut a_max_tol = 0.0f64;
    let mut vec_curve: Vec<OrientedCurve> = Vec::new();

    // `cxx:89-133`: `BRepTools_WireExplorer anExpE(theWire, theFace)`.
    let mut wexp = WireExplorer::of_wire_on_face(w, face);
    while wexp.more() {
        let edge = wexp.current_edge()?;
        // `cxx:94-99`: `aTol = max(BRep_Tool::Tolerance(anEdge), MINIMAL_TOLERANCE)`.
        let a_tol = BRepTool::edge_tolerance(&edge).max(MINIMAL_TOLERANCE);
        if a_tol > a_max_tol {
            a_max_tol = a_tol;
        }
        // `cxx:93`: `BRepAdaptor_Curve2d aCurve(anEdge, theFace)`. The adaptor
        // is the stored p-curve plus the CurveOnSurface range
        // (`Geom2dAdaptor_Curve::Load`); `Geom2dTrimmedCurve` is the port's
        // carrier for exactly that pair (basis + `[First, Last]`) and forwards
        // the type queries the approximation branches on.
        let (pc, first, last) = crate::boptools_2d::curve_on_surface_range(&edge, face)?;
        let adaptor: Arc<dyn Curve2d> =
            Arc::new(Geom2dTrimmedCurve::new_sense(pc, first, last, true, false));
        // `cxx:102`: `Geom2dConvert_ApproxArcsSegments anAlgo(aCurve, aTol, theAngleTol)`.
        let algo = ApproxArcsSegments::new(adaptor.as_ref(), a_tol, angle_tol);
        let approx = algo.result();

        // `cxx:108-133`: the pieces in edge order, or in reverse edge order
        // when the edge is REVERSED.
        if edge.0.orientation() == Orientation::Reversed {
            for crv in approx.iter().rev() {
                if let Some(crv) = crv {
                    vec_curve.push(OrientedCurve { curve: crv.clone(), is_reverse: true });
                }
            }
        } else {
            for crv in approx.iter() {
                if let Some(crv) = crv {
                    vec_curve.push(OrientedCurve { curve: crv.clone(), is_reverse: false });
                }
            }
        }
        wexp.next();
    }

    // `cxx:136`: nothing approximated -> null wire.
    if vec_curve.is_empty() {
        return None;
    }

    let builder = TopoBuilder::new();

    // `cxx:140-154`: the first vertex, shared by the wire's two ends when the
    // mismatch is within `aMaxTol + Precision::Confusion()`.
    let mut a_pnt = [vec_curve[0].point(false), vec_curve[vec_curve.len() - 1].point(true)];
    let mut a_dist = a_pnt[0].distance(&a_pnt[1]);
    if a_dist > a_max_tol + CONFUSION {
        // `cxx:144-146`: tolerances do not match -> a plain `Confusion` vertex.
        a_dist = CONFUSION;
    } else {
        // `cxx:147-151`.
        a_dist = 0.5 * a_dist + CONFUSION;
        a_pnt[0] = GpPnt2d::from_xy(a_pnt[0].coord.added(&a_pnt[1].coord).multiplied_scalar(0.5));
    }
    // `cxx:152-155`: `aSurf->D0(aPnt[0].X(), aPnt[0].Y(), aPnt3d)`.
    let p3d = surface.d0(a_pnt[0].x(), a_pnt[0].y());
    let first_vertex = builder.make_vertex(p3d, a_dist);

    // `cxx:158-211`: the edges.
    let mut mk_wire = MakeWire::new();
    let mut vertex = first_vertex.clone();
    for i in 0..vec_curve.len() {
        let ocurve = &vec_curve[i];
        a_pnt[0] = ocurve.point(true);
        let next_vertex;
        if i == vec_curve.len() - 1 {
            // `cxx:167-178`: the closing point.
            a_pnt[1] = vec_curve[0].point(false);
            a_dist = a_pnt[0].distance(&a_pnt[1]);
            if a_dist > a_max_tol + CONFUSION {
                let p = surface.d0(a_pnt[0].x(), a_pnt[0].y());
                next_vertex = builder.make_vertex(p, CONFUSION);
            } else {
                next_vertex = first_vertex.clone();
            }
        } else {
            // `cxx:180-187`: the junction with the next piece.
            a_pnt[1] = vec_curve[i + 1].point(false);
            a_dist = 0.5 * a_pnt[0].distance(&a_pnt[1]) + CONFUSION;
            a_pnt[0] =
                GpPnt2d::from_xy(a_pnt[0].coord.added(&a_pnt[1].coord).multiplied_scalar(0.5));
            let p = surface.d0(a_pnt[0].x(), a_pnt[0].y());
            next_vertex = builder.make_vertex(p, a_dist);
        }
        // `cxx:188`: the parameter range of the piece.
        let (param0, param1) =
            (ocurve.curve.first_parameter(), ocurve.curve.last_parameter());
        let (v1, v2) = if ocurve.is_reverse {
            // `cxx:190-198`.
            (&next_vertex, &vertex)
        } else {
            // `cxx:200-207`.
            (&vertex, &next_vertex)
        };
        let mut edge = make_edge_2d(&builder, &ocurve.curve, &surface, face, v1, v2, param0, param1)?;
        if ocurve.is_reverse {
            // `cxx:197`: `anEdgeRes.Orientation(TopAbs_REVERSED)`.
            edge.0.set_orientation(Orientation::Reversed);
        }
        vertex = next_vertex;
        // `cxx:210`: `aMkWire.Add(anEdgeRes)`.
        let _ = mk_wire.add(&edge);
    }

    // `cxx:213-216`: `if (aMkWire.IsDone()) aResult = aMkWire.Wire();`.
    if mk_wire.is_done() {
        Some(mk_wire.wire())
    } else {
        None
    }
}

/// `BRepAlgo::ConvertFace(theFace, theAngleTolerance)` (`cxx:223-242`).
///
/// `BRepBuilderAPI_MakeFace(aSurf, Precision::Confusion())` (`cxx:229`) is
/// `BRepLib_MakeFace(S, TolDegen)` -> `Init(S, Standard_True, TolDegen)`
/// (`BRepLib_MakeFace.cxx:110-113`), i.e. a face bounded by the surface's
/// natural UV rectangle, natural-restriction set; every converted wire is then
/// added (`cxx:232-235`) and the face's orientation is restored (`cxx:235`).
///
/// `None` is the null `TopoDS_Face` when no surface is available; OCCT builds a
/// surface-less face there and dereferences the null handle in `ConvertWire`.
pub fn convert_face(face: &Face, angle_tol: f64) -> Option<Face> {
    // `cxx:228`.
    let surface = BRepTool::face_surface(face)?;
    // `cxx:229`.
    let mut result = brep_lib_make_face::make_face_from_surface(surface, CONFUSION);
    let reg = GeometryRegistry::global();
    // `BRepLib_MakeFace::Init(S, Standard_True, TolDegen)`
    // (`BRepLib_MakeFace.cxx:366-386`): `B.NaturalRestriction(face, true)`.
    reg.set_natural_restriction(&result.0, true);
    let builder = TopoBuilder::new();
    // `cxx:68` of `BRepOffsetAPI_MakeOffset::ConvertFace` is the other caller's
    // FORWARD step; here `cxx:230-233` passes the wires as `TopExp_Explorer`
    // yields them, and `ConvertWire` uses the face's surface either way.
    // `cxx:229-235`: every direct wire of the face.
    for sub in explorer_shapes(&face.0, ShapeType::Wire) {
        let wire = Wire(sub);
        let converted = convert_wire(&wire, angle_tol, face);
        // `BRepLib_MakeFace::Add` (`BRepLib_MakeFace.cxx:873-879`) adds the
        // wire and clears `NaturalRestriction` even when the wire is null.
        reg.set_natural_restriction(&result.0, false);
        // `cxx:232-234`: a null converted wire makes `BRep_Builder::Add` a
        // no-op, so it is skipped here.
        if let Some(new_wire) = converted {
            builder.add_wire(&mut result, &new_wire);
        }
    }
    // `cxx:237`: `aNewFace.Orientation(anOr)`.
    result.0.set_orientation(face.0.orientation());
    // `cxx:237-241`: `aMkFace.IsDone()` is `FaceDone` as soon as `Init` ran.
    Some(result)
}
