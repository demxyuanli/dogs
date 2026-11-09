//! The static services of `HLRBRep` (`ModelingAlgorithms/TKHLR/HLRBRep/HLRBRep.hxx`,
//! `HLRBRep.cxx:37-258`): turning a projected `HLRBRep_Curve` segment back into a
//! `TopoDS_Edge`.
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `MakeEdge(ec, U1, U2)` | `cxx:37-183` | [`make_edge`] |
//! | `MakeEdge3d(ec, U1, U2)` | `cxx:187-233` | [`make_edge_3d`] |
//! | `PolyHLRAngleAndDeflection` | `cxx:239-258` | [`poly_hlr_angle_and_deflection`] |
//!
//! # Translation notes
//!
//! * `MakeEdge` is `BRepLib_MakeEdge2d(C, sta, end)` in every arm. The port
//!   builds the same 2D curve object and runs it through
//!   `BRepLib_MakeEdge::Init(C2d, Plane, V1, V2, p1, p2)`
//!   (`BRepLib_MakeEdge.cxx:905-1094`), then stores it as a p-curve on
//!   `BRepLib::Plane()` with `BRep_Builder::UpdateEdge`
//!   (`BRepLib_MakeEdge2d.cxx:619-633`) - see [`make_edge2d`].
//!   The three `IsDone()`-guarded arms (`cxx:79-84`, `:148-153`, `:177-181`)
//!   answer a null edge on failure, i.e. `None` here.
//! * `MakeEdge3d` copies the source edge (`TopoDS::Edge(anEdge.EmptyCopied())`,
//!   `cxx:202`) and only rewrites the range and the two vertices. The port's
//!   geometry lives in `GeometryRegistry`, so the copy re-registers the source
//!   `EdgeGeom` - the same `BRep_TEdge::EmptyCopy` content (`BRep_TEdge.cxx:105-133`).
//! * `BRepAdaptor_Curve BAcurve(anEdge)` (`cxx:197`) is
//!   [`occt_topo::brep_tool::BRepTool::edge_curve_adaptor_world`] plus the edge's
//!   registered range, which is what the adaptor evaluates.
//! * `PolyHLRAngleAndDeflection(InAngl, OutAngl, OutDefl)` returns the two
//!   outputs as a tuple; the three static bound variables (`cxx:241-243`) are
//!   the literals of `cxx:241-243` in radians.

use std::sync::Arc;

use occt_core::bspl::banded_interp::knot_sequence;
use occt_core::bspl::knots::knot_sequence_periodic;
use occt_core::gp::{GpLin2d, GpPnt2d};
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom2d::bezier_curve::Geom2dBezierCurve;
use occt_geom2d::bspline_curve::Geom2dBSplineCurve;
use occt_geom2d::circle::Geom2dCircle;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::ellipse::Geom2dEllipse;
use occt_geom2d::hyperbola::Geom2dHyperbola;
use occt_geom2d::line::Geom2dLine;
use occt_geom2d::parabola::Geom2dParabola;

use occt_topo::abs::Orientation;
use occt_topo::brep_lib_make_edge;
use occt_topo::brep_tool::BRepTool;
use occt_topo::builder::TopoBuilder;
use occt_topo::geom_abs_curve_type::GeomAbsCurveType;
use crate::hlrbrep_curve::HlrBRepCurve;
use crate::hlrbrep_poly_hlr_to_shape::brep_lib_plane;
use occt_topo::shape::{Edge, Vertex};
use occt_topo::tgeometry::GeometryRegistry;

/// `BRepLib_MakeEdge2d(const Handle(Geom2d_Curve)& C, U1, U2)`
/// (`BRepLib_MakeEdge2d.cxx:387-403` then `:465-633`) written the way
/// [`crate::hlrbrep_poly_hlr_to_shape::make_edge_2d`] writes the two-point
/// overload: `Init(C2d, Plane, V1=null, V2=null, p1, p2)`
/// (`BRepLib_MakeEdge.cxx:905-1094`) stores the p-curve on `BRepLib::Plane()`
/// and the vertices at `BRepLib::Precision()`, then `B.Range(E, p1, p2)`
/// (`cxx:630`) writes the range to every representation.
///
/// `None` is the `!IsDone()` of `ParameterOutOfRange`
/// (`BRepLib_MakeEdge.cxx:960-964`), which `MakeEdge` reports as a null edge.
fn make_edge2d(curve: Arc<dyn Curve2d>, first: f64, last: f64) -> Option<Edge> {
    let plane = brep_lib_plane();
    let build = brep_lib_make_edge::init_2d(&curve, plane, first, last)?;
    let mut edge = Edge::new();
    let reg = GeometryRegistry::global();
    reg.set_ranges_all(&edge.0, build.first, build.last);
    // `B.UpdateEdge(E, C, BRepLib::Plane(), TopLoc_Location(), preci)`
    // (`BRepLib_MakeEdge2d.cxx:621`): the edge carries the p-curve only.
    reg.update_edge_pcurve_on_surface(&edge.0, plane, Some(curve), CONFUSION);
    reg.set_pcurve_range_on_surface(&edge.0, plane, build.first, build.last);
    let b = TopoBuilder::new();
    brep_lib_make_edge::add_endpoint_vertices(&b, &mut edge, &build);
    if build.reverse {
        edge.0.set_orientation(Orientation::Reversed);
    }
    Some(edge)
}

/// `BSplCLib::KnotSequence` / `BSplCLib::KnotSequence(Knots, Mults, Deg,
/// Periodic)`: the flat knot vector the `Geom2d_BSplineCurve` constructor
/// derives in `updateKnots()`.
fn flat_knots(knots: &[f64], mults: &[i32], degree: usize, periodic: bool) -> Vec<f64> {
    if periodic {
        knot_sequence_periodic(knots, mults, degree as i32)
    } else {
        knot_sequence(knots, mults, degree as i32)
    }
}

/// `Geom2d_BSplineCurve(Poles, Knots, Mults, Degree, IsPeriodic)`
/// (`Geom2d_BSplineCurve.cxx:71-118`), rational and non-rational arms.
fn bspline2d(
    poles: Vec<GpPnt2d>,
    weights: Option<Vec<f64>>,
    knots: &[f64],
    mults: &[i32],
    degree: usize,
    periodic: bool,
) -> Option<Geom2dBSplineCurve> {
    let xs: Vec<f64> = poles.iter().map(|p| p.x()).collect();
    let ys: Vec<f64> = poles.iter().map(|p| p.y()).collect();
    let flat = flat_knots(knots, mults, degree, periodic);
    Geom2dBSplineCurve::from_flat(xs, ys, weights, flat, degree, periodic).ok()
}

/// `HLRBRep::MakeEdge(const HLRBRep_Curve& ec, const double U1, const double U2)`
/// (`HLRBRep.cxx:37-183`): the 2D edge of the projected curve between two
/// parameters.
///
/// `None` is the null `TopoDS_Edge` the source returns when the
/// `BRepLib_MakeEdge2d` construction fails (`cxx:79-84`, `:148-153`, `:177-181`)
/// or when the curve is not one of the handled types (returns `Edg` untouched,
/// `cxx:182`).
pub fn make_edge(ec: &HlrBRepCurve, u1: f64, u2: f64) -> Option<Edge> {
    // `cxx:40-41`: the two 2D parameters of the segment ends.
    let sta = ec.parameter2d(u1);
    let end = ec.parameter2d(u2);

    match ec.get_type() {
        // `cxx:44-47`: `BRepLib_MakeEdge2d(ec.Line(), sta, end)`.
        GeomAbsCurveType::Line => {
            let lin: GpLin2d = ec.line()?;
            make_edge2d(Arc::new(Geom2dLine::new(lin.pos)), sta, end)
        }
        // `cxx:49-51`.
        GeomAbsCurveType::Circle => {
            let c = ec.circle()?;
            make_edge2d(Arc::new(Geom2dCircle::new(c)), sta, end)
        }
        // `cxx:53-55`.
        GeomAbsCurveType::Ellipse => {
            let e = ec.ellipse()?;
            make_edge2d(Arc::new(Geom2dEllipse::new(e)), sta, end)
        }
        // `cxx:57-59`.
        GeomAbsCurveType::Hyperbola => {
            let h = ec.hyperbola();
            make_edge2d(Arc::new(Geom2dHyperbola::new(h)), sta, end)
        }
        // `cxx:61-63`.
        GeomAbsCurveType::Parabola => {
            let p = ec.parabola();
            make_edge2d(Arc::new(Geom2dParabola::new(p)), sta, end)
        }
        // `cxx:65-85`: `Geom2d_BezierCurve(Poles[, Weights])`.
        GeomAbsCurveType::BezierCurve => {
            let curve = if ec.is_rational() {
                let (poles, weights) = ec.poles_and_weights();
                Geom2dBezierCurve::rational(poles, weights).ok()?
            } else {
                Geom2dBezierCurve::new(ec.poles()).ok()?
            };
            make_edge2d(Arc::new(curve), sta, end)
        }
        // `cxx:87-154`.
        GeomAbsCurveType::BSplineCurve => {
            let curve = bspline_arm(ec, sta, end)?;
            make_edge2d(Arc::new(curve), sta, end)
        }
        // `cxx:155-181`: 15 sampled points through a degree-1 B-spline.
        _ => {
            const NB_PNT: usize = 15;
            let mut poles: Vec<GpPnt2d> = Vec::with_capacity(NB_PNT);
            let mut knots: Vec<f64> = Vec::with_capacity(NB_PNT);
            let mut mults: Vec<i32> = vec![1; NB_PNT];
            // `cxx:161-162`.
            mults[0] = 2;
            mults[NB_PNT - 1] = 2;
            // `cxx:163-171`.
            let step = (u2 - u1) / (NB_PNT as f64 - 1.0);
            let mut par3d = u1;
            for _ in 0..(NB_PNT - 1) {
                poles.push(ec.value(par3d));
                knots.push(par3d);
                par3d += step;
            }
            // `cxx:172-173`.
            poles.push(ec.value(u2));
            knots.push(u2);
            // `cxx:175`: `Geom2d_BSplineCurve(Poles, knots, mults, 1)`.
            let curve = bspline2d(poles, None, &knots, &mults, 1, false)?;
            make_edge2d(Arc::new(curve), sta, end)
        }
    }
}

/// `HLRBRep.cxx:88-153`: the B-spline arm of [`make_edge`].
///
/// The periodic sub-branch (`cxx:100-126`) segments a **copy of the 3D
/// B-spline** and reads the post-segment poles, knots, multiplicities and
/// degree; the other one (`cxx:127-147`) reads them straight from
/// `HLRBRep_Curve`. `None` is a failed `BRep_Tool::Curve` / `down_cast`
/// (`cxx:91-97`), which OCCT answers with a dereferenced null handle.
fn bspline_arm(ec: &HlrBRepCurve, sta: f64, end: f64) -> Option<Geom2dBSplineCurve> {
    // `cxx:89-90`: `ec.GetCurve().Curve()` and `ec.GetCurve().Edge()`.
    let an_edge = ec.edge()?.clone();
    // `cxx:92-94`: `BRep_Tool::Curve(anEdge, fpar, lpar)`.
    let mut a_curve = BRepTool::edge_curve(&an_edge)?;
    // `cxx:95-98`: a `Geom_TrimmedCurve` is replaced by its basis curve.
    if a_curve.is_geom_trimmed() {
        a_curve = a_curve.trim_basis()?;
    }
    // `cxx:98-99`: `occ::down_cast<Geom_BSplineCurve>(aCurve)`.
    let bspl = a_curve.as_bspline_curve().cloned()?;
    // `cxx:99-100`: `BSplCurve->Copy()`.
    let mut the_curve = bspl;

    // `cxx:100`: `theCurve->IsPeriodic() && !GAcurve.IsClosed()`.
    if the_curve.periodic() && !ec.is_closed() {
        // `cxx:102`: `theCurve->Segment(sta, end)`.
        the_curve.segment(sta, end, PCONFUSION).ok()?;
        // `cxx:103-107`.
        let nb_poles = the_curve.nb_poles();
        let poles: Vec<GpPnt2d> = the_curve_poles(ec, &the_curve)
            .into_iter()
            .take(nb_poles)
            .collect();
        let knots = the_curve.distinct_knots().to_vec();
        let mults = the_curve.multiplicities().to_vec();
        let weights = if the_curve.is_rational() {
            Some(the_curve_weights(&the_curve).into_iter().take(nb_poles).collect())
        } else {
            None
        };
        // `cxx:109-126`: both constructors pass `theCurve->IsPeriodic()`,
        // which `Segment` has already reset to `false` (`Geom_BSplineCurve.cxx:617`).
        return bspline2d(
            poles,
            weights,
            &knots,
            &mults,
            the_curve.degree(),
            the_curve.periodic(),
        );
    }

    // `cxx:127-147`: the poles/knots/mults of `HLRBRep_Curve` itself.
    let knots = ec.knots();
    let mults = ec.multiplicities();
    let degree = ec.degree().max(0) as usize;
    let (poles, weights) = if ec.is_rational() {
        let (p, w) = ec.poles_and_weights();
        (p, Some(w))
    } else {
        (ec.poles(), None)
    };
    bspline2d(poles, weights, &knots, &mults, degree, ec.is_periodic())
}

/// The 2D poles of the given 3D B-spline, projected and reduced to `X`/`Y`
/// (`HLRBRep_Curve::Poles(C, TP)`, `HLRBRep_Curve.cxx:511-528`).
fn the_curve_poles(ec: &HlrBRepCurve, c: &occt_geom::GeomBSplineCurve) -> Vec<GpPnt2d> {
    let mut out = Vec::with_capacity(c.nb_poles());
    for i in 0..c.nb_poles() {
        let p = *c.pole(i);
        out.push(ec.project_to_2d(&p));
    }
    out
}

/// `HLRBRep_Curve::PolesAndWeights(C, TP, TW)` (`HLRBRep_Curve.cxx:579-602`):
/// `Geom_BSplineCurve::WeightsArray()`, i.e. the stored weights.
fn the_curve_weights(c: &occt_geom::GeomBSplineCurve) -> Vec<f64> {
    c.weights
        .clone()
        .unwrap_or_else(|| vec![1.0; c.nb_poles()])
}

/// All vertices of `an_edge`, split into the last FORWARD and the last REVERSED
/// ones, as `TopExp::Vertices(E, V1, V2)` picks them (`TopExp.cxx:214-251`).
fn edge_vertices(an_edge: &Edge) -> (Option<Vertex>, Option<Vertex>) {
    let mut v_first = None;
    let mut v_last = None;
    for v in occt_topo::bop_occt_util::vertices_of_edge_direct(&an_edge.0) {
        match v.0.orientation() {
            Orientation::Forward => v_first = Some(v),
            Orientation::Reversed => v_last = Some(v),
            _ => {}
        }
    }
    (v_first, v_last)
}

/// `HLRBRep::MakeEdge3d(const HLRBRep_Curve& ec, const double U1, const double U2)`
/// (`HLRBRep.cxx:187-233`): the 3D edge of the underlying topological edge,
/// restricted to `[U1, U2]` with vertices shared when the parameter matches the
/// edge's own end.
///
/// `None` is a missing source edge or 3D curve (`cxx:193` dereferences both).
pub fn make_edge_3d(ec: &HlrBRepCurve, u1: f64, u2: f64) -> Option<Edge> {
    // `cxx:193`: `ec.GetCurve().Edge()`.
    let an_edge = ec.edge()?.clone();
    // `cxx:196-199`: `BRepAdaptor_Curve BAcurve(anEdge)`.
    let ba_curve = BRepTool::edge_curve_adaptor_world(&an_edge)?;
    let (fpar, lpar) = BRepTool::edge_parameters(&an_edge);

    // `cxx:202-205`: `TopoDS::Edge(anEdge.EmptyCopied())`,
    // `Edg.Orientation(TopAbs_FORWARD)`, `BB.Range(Edg, U1, U2)`.
    let mut edg = Edge::wrap(an_edge.0.empty_copied())?;
    if let Some(geom) = GeometryRegistry::global().edge_geom(&an_edge.0) {
        GeometryRegistry::global().set_edge(&edg.0, geom);
    }
    edg.0.set_orientation(Orientation::Forward);
    GeometryRegistry::global().set_ranges_all(&edg.0, u1, u2);

    // `cxx:208-231`: share the endpoints when the range end is unchanged,
    // otherwise build a fresh vertex at the evaluated point.
    let (v1, v2) = edge_vertices(&an_edge);
    let b = TopoBuilder::new();
    let v1_new = if (fpar - u1).abs() <= PCONFUSION {
        v1
    } else {
        Some(b.make_vertex(ba_curve.d0(u1), CONFUSION))
    };
    let v2_new = if (lpar - u2).abs() <= PCONFUSION {
        v2
    } else {
        Some(b.make_vertex(ba_curve.d0(u2), CONFUSION))
    };

    // `cxx:230-232`: `V1new.Orientation(FORWARD)`, `V2new.Orientation(REVERSED)`,
    // `BB.Add(Edg, V1new)`, `BB.Add(Edg, V2new)`.
    if let Some(v) = &v1_new {
        b.add(&mut edg.0, &v.0.oriented(Orientation::Forward));
    }
    if let Some(v) = &v2_new {
        b.add(&mut edg.0, &v.0.oriented(Orientation::Reversed));
    }
    Some(edg)
}

/// `HLRBRep::PolyHLRAngleAndDeflection(InAngl, OutAngl, OutDefl)`
/// (`HLRBRep.cxx:239-258`): the polygonal HLR angle the caller asks for,
/// clamped to `[1 deg, 35 deg]` and mapped onto `[5 deg, 35 deg]`, and the
/// deflection `OutAngl^2 / 2` derived from it.
pub fn poly_hlr_angle_and_deflection(in_angl: f64) -> (f64, f64) {
    // `cxx:241-243`: the three static bounds, in radians.
    let h_ang_min = 1.0 * std::f64::consts::PI / 180.0;
    let h_ang_lim = 5.0 * std::f64::consts::PI / 180.0;
    let h_ang_max = 35.0 * std::f64::consts::PI / 180.0;

    let mut out_angl = in_angl;
    // `cxx:246-252`.
    if out_angl < h_ang_min {
        out_angl = h_ang_min;
    }
    if out_angl > h_ang_max {
        out_angl = h_ang_max;
    }
    // `cxx:253-256`.
    out_angl = h_ang_lim
        + ((out_angl - h_ang_min) * (h_ang_max - h_ang_lim) * (h_ang_max - h_ang_lim)
            / (h_ang_max - h_ang_min))
            .sqrt();
    // `cxx:257`.
    let out_defl = out_angl * out_angl * 0.5;
    (out_angl, out_defl)
}
