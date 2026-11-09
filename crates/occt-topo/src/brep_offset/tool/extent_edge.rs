//! File-local `ExtentEdge` (`BRepOffset_Tool.cxx:2072-2157`),
//! `ProjectVertexOnEdge` (`:2161-2249`), `SelectEdge` (`:2536-2576`) and
//! `TryParameter` (`:3638-3681`), the edge-level helpers `ExtentFace`
//! (`:3793-4333`) and `Inter2d` (`:2253-2532`) call.
//!
//! `ExtentEdge` extends an edge's pcurve on one face in tangency until it
//! reaches the surface's UV box (clamped to `+-1.e2`), then rebuilds the edge
//! with that extended pcurve and no 3D curve. It answers an *empty* edge for
//! the analytic curve types (`cxx:2085-2089`) and for a `Geom_BoundedCurve`
//! down-cast failure (`cxx:2139-2142`) - the caller's null test is on the
//! `TopoDS_Edge` handle, which stays non-null in both arms, so the port
//! returns the empty edge as well.
//!
//! `ExtentEdge` reads the surface of `F` (`cxx:2092`) but attaches the new
//! pcurve to `EF` (`cxx:2154`): the parameter space is shared because `EF` is
//! an enlargement of `F` over the same surface. The two arguments are not
//! interchangeable and both are used, each at its own call.

use std::sync::Arc;

use occt_core::convert::ParameterisationType;
use occt_core::gp::{GpPln, GpPnt};
use occt_core::precision::{Precision, CONFUSION};
use occt_geom::convert_bspl::curve_to_bspline_curve;
use occt_geom::extrema_pc::extrema_ext_pc_range;
use occt_geom::geom_api::to_2d;
use occt_geom::geom_lib::to_3d;
use occt_geom::geom_lib_extend_curve::extend_curve_to_point;
use occt_geom::Curve;
use occt_geom2d::Curve2d;

use crate::abs::{Orientation, ShapeType};
use crate::brep_builder_update_vertex::update_vertex_on_edge;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::geom_abs_curve_type::{curve_kind, GeomAbsCurveType};
use crate::shape::{Edge, Face, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topopebrep::tool::adaptor::GeomAdaptorCurve;

use super::BRepOffsetTool;

impl BRepOffsetTool {
    /// `ExtentEdge(F, EF, E, NE)` (`cxx:2072-2157`).
    ///
    /// `f` is the surface whose bounds bound the extension, `ef` the face that
    /// receives the new pcurve. A `None` is the null `BRep_Tool::Surface(F)`
    /// handle OCCT dereferences at `cxx:2095`; a missing `CurveOnSurface(E, F)`
    /// is the null `C2d` dereferenced at `cxx:2107`.
    pub fn extent_edge(f: &Face, ef: &Face, e: &Edge) -> Option<Edge> {
        let builder = TopoBuilder::new();
        // `TopoDS_Shape aLocalEdge = E.EmptyCopied(); NE = TopoDS::Edge(...)`
        // (`cxx:2078-2080`): an edge with a fresh empty TShape. `B.MakeEdge(NE)`
        // at `cxx:2152` re-makes it anyway, so the port starts from an empty edge.
        let empty_edge = || Edge(builder.make_shape(ShapeType::Edge));

        let c = BRepTool::edge_curve_world(e)?;
        // `BRepAdaptor_Curve CE(E); CE.GetType()` (`cxx:2076-2077`).
        let kind = curve_kind(c.as_ref());
        if matches!(
            kind,
            GeomAbsCurveType::Line
                | GeomAbsCurveType::Circle
                | GeomAbsCurveType::Ellipse
                | GeomAbsCurveType::Hyperbola
                | GeomAbsCurveType::Parabola
        ) {
            // `cxx:2085-2089`.
            return Some(empty_edge());
        }

        // `cxx:2089-2098`: `PMax = 1.e2` and `S->Bounds` clamped to it.
        const P_MAX: f64 = 1.0e2;
        let s = BRepTool::face_surface(f)?;
        let (mut umin, mut umax) = s.u_range();
        let (mut vmin, mut vmax) = s.v_range();
        umin = umin.max(-P_MAX);
        vmin = vmin.max(-P_MAX);
        umax = umax.min(P_MAX);
        vmax = vmax.min(P_MAX);

        // `cxx:2100-2101`: the pcurve and its COS range.
        let (c2d, _cf, _cl) = crate::boptools_2d::curve_on_surface_range(e, f)?;
        // `BRepAdaptor_Curve CE(E)` bounds: `CE.FirstParameter()` / `LastParameter()`.
        let (par_first, par_last) = BRepTool::edge_parameters(e);

        // `cxx:2103-2120`: the tangent extension at the first parameter, aimed at
        // the UV box. `tx` / `ty` start at `Precision::Infinite()` and only the
        // finite-direction case is considered.
        let (p, tang) = c2d.d1(par_first);
        let mut tx = Precision::INFINITE;
        let mut ty = Precision::INFINITE;
        if tang.x().abs() > CONFUSION {
            tx = ((umax - p.x()) / tang.x())
                .abs()
                .min(((umin - p.x()) / tang.x()).abs());
        }
        if tang.y().abs() > CONFUSION {
            ty = ((vmax - p.y()) / tang.y())
                .abs()
                .min(((vmin - p.y()) / tang.y()).abs());
        }
        let tmin = tx.min(ty);
        // `Tang = tmin * Tang; PF2d(P.X() - Tang.X(), P.Y() - Tang.Y())`.
        let pf2d = (p.x() - tmin * tang.x(), p.y() - tmin * tang.y());

        // `cxx:2122-2135`: the same at the last parameter, but with a `+`.
        let (p, tang) = c2d.d1(par_last);
        let mut tx = Precision::INFINITE;
        let mut ty = Precision::INFINITE;
        if tang.x().abs() > CONFUSION {
            tx = ((umax - p.x()) / tang.x())
                .abs()
                .min(((umin - p.x()) / tang.x()).abs());
        }
        if tang.y().abs() > CONFUSION {
            ty = ((vmax - p.y()) / tang.y())
                .abs()
                .min(((vmin - p.y()) / tang.y()).abs());
        }
        let tmin = tx.min(ty);
        let pl2d = (p.x() + tmin * tang.x(), p.y() + tmin * tang.y());

        // `cxx:2137-2142`: `GeomAPI::To3d(C2d, gp_Pln(gp::XOY()))` then the
        // `Geom_BoundedCurve` down-cast. `Geom_BoundedCurve` covers B-spline,
        // Bezier and trimmed curves only, so anything else answers the empty edge.
        let xoy_pos = occt_core::gp::GpAx2::standard();
        let cc = to_3d(&xoy_pos, c2d.as_ref())?;
        let is_bounded = matches!(
            curve_kind(cc.as_ref()),
            GeomAbsCurveType::BSplineCurve | GeomAbsCurveType::BezierCurve
        ) || cc.untrimmed_basis().is_some();
        if !is_bounded {
            return Some(empty_edge());
        }

        // `cxx:2145-2147`: `GeomLib::ExtendCurveToPoint(ExtC, PF, 1, false)` then
        // `(ExtC, PL, 1, true)`. OCCT's helper internally converts through
        // `GeomConvert_CompCurveToBSplineCurve(Convert_QuasiAngular)`
        // (`GeomLib.cxx:1286`) and reassigns the handle to that B-spline
        // (`cxx:1412`), which is what the port does explicitly.
        let mut bs =
            curve_to_bspline_curve(cc.as_ref(), ParameterisationType::QuasiAngular).ok()?;
        extend_curve_to_point(&mut bs, &GpPnt::new(pf2d.0, pf2d.1, 0.0), 1, false).ok()?;
        extend_curve_to_point(&mut bs, &GpPnt::new(pl2d.0, pl2d.1, 0.0), 1, true).ok()?;
        let ext_c: Arc<dyn Curve> = Arc::new(bs);

        // `cxx:2149`: `GeomAPI::To2d(ExtC, gp_Pln(gp::XOY()))`.
        let xoy = GpPln::default();
        let cne2d: Arc<dyn Curve2d> = to_2d(ext_c.as_ref(), &xoy)?;

        // `cxx:2151-2156`.
        let mut ne = builder.make_shape(ShapeType::Edge);
        let ne_edge = Edge(ne.clone());
        let reg = GeometryRegistry::global();
        let ef_key = GeometryRegistry::shape_key(&ef.0);
        let tol = BRepTool::edge_tolerance(e);
        reg.set_edge_pcurve(&ne_edge.0, ef_key, cne2d.clone());
        reg.set_ranges_all(&ne_edge.0, cne2d.first_parameter(), cne2d.last_parameter());
        reg.set_edge_tolerance(&ne_edge.0, tol);
        // `NE.Orientation(E.Orientation())`.
        ne.set_orientation(e.0.orientation());
        Some(Edge(ne))
    }

    /// `ProjectVertexOnEdge(V, E, TolConf)` (`cxx:2161-2249`): locate `V` on `E`.
    ///
    /// Order of the three arms, all from the reference:
    /// 1. a `FORWARD` vertex whose point is the edge's first point (`cxx:2175-2184`),
    /// 2. a `REVERSED` vertex whose point is the edge's last point (`cxx:2186-2198`),
    /// 3. the closest `Extrema_ExtPC` extremum (`cxx:2199-2215`).
    ///
    /// On a hit `B.UpdateVertex(V.Oriented(INTERNAL), U, E.Oriented(FORWARD), tol)`
    /// is applied (`cxx:2240-2247`). `v` is taken by reference and mutated in
    /// place, like OCCT's `TopoDS_Vertex&`.
    ///
    /// The `OCCT_DEBUG` prints (`cxx:2218-2239`) are not ported; they have no
    /// effect on the returned flag.
    pub fn project_vertex_on_edge(v: &Vertex, e: &Edge, tol_conf: f64) -> bool {
        let p = BRepTool::vertex_point(v);

        // `BRepAdaptor_Curve C(E); f = C.FirstParameter(); l = C.LastParameter()`
        // (`cxx:2170-2172`).
        let Some(curve) = BRepTool::edge_curve_world(e) else {
            // Port guard: `BRepAdaptor_Curve` over an edge without a 3D curve has
            // nothing to evaluate, while OCCT raises from `GeomAdaptor_Curve`.
            return false;
        };
        let adaptor = GeomAdaptorCurve::of_curve(curve.clone());
        let f = adaptor.first_parameter();
        let l = adaptor.last_parameter();

        let mut u = 0.0;
        let mut found = false;

        if v.0.orientation() == Orientation::Forward && !Precision::is_infinite(f) {
            let pf = adaptor.d0(f);
            if pf.is_equal_with_tolerance(&p, tol_conf) {
                u = f;
                found = true;
            }
        }
        if v.0.orientation() == Orientation::Reversed && !found && !Precision::is_infinite(l) {
            let pl = adaptor.d0(l);
            if pl.is_equal_with_tolerance(&p, tol_conf) {
                u = l;
                found = true;
            }
        }
        if !found {
            // `Extrema_ExtPC Proj(P, C)`, then the extremum with the smallest
            // `SquareDistance` (`cxx:2200-2215`).
            let sols = extrema_ext_pc_range(&adaptor, &p, f, l);
            let mut best: Option<(f64, f64)> = None;
            for sol in &sols {
                let mut d = 0.0;
                for (a, b) in [
                    (sol.point.x(), p.x()),
                    (sol.point.y(), p.y()),
                    (sol.point.z(), p.z()),
                ] {
                    d += (a - b) * (a - b);
                }
                if best.map(|(bd, _)| d < bd).unwrap_or(true) {
                    best = Some((d, sol.u));
                }
            }
            if let Some((_, par)) = best {
                u = par;
                found = true;
            }
        }

        if found {
            // `cxx:2242-2247`.
            let ee = e.0.oriented(Orientation::Forward);
            let vv = v.0.oriented(Orientation::Internal);
            let tol = BRepTool::edge_tolerance(e);
            let _ = update_vertex_on_edge(&vv, u, &ee, tol);
        }
        found
    }

    /// `SelectEdge(F, EF, E, LInt)` (`cxx:2536-2576`, the two faces unused):
    /// of the section edges in `lint`, keep the single one whose endpoints are
    /// closest to `e`'s endpoints. `lint` is rewritten in place, like OCCT's
    /// `NCollection_List&`.
    ///
    /// When every candidate is farther than the running `dU = 1.e100` the
    /// list is still cleared and `GE` stays null in OCCT (`cxx:2573-2575`);
    /// the port then leaves the list empty, which is what the caller's
    /// `LInt1.First()` would have faulted on.
    pub fn select_edge(e: &Edge, lint: &mut Vec<Edge>) {
        let mut du = 1.0e100;
        let mut ge: Option<Edge> = None;

        let (fst, lst) = BRepTool::edge_parameters(e);
        let Some(c1) = BRepTool::edge_curve_world(e) else {
            lint.clear();
            return;
        };
        let p_first = c1.d0(fst);
        let p_last = c1.d0(lst);

        for ei in lint.iter() {
            let (fst, lst) = BRepTool::edge_parameters(ei);
            let Some(c2) = BRepTool::edge_curve_world(ei) else {
                continue;
            };
            let p1 = c2.d0(fst);
            let p2 = c2.d0(lst);
            let tmp = p1.distance(&p_first) + p2.distance(&p_last);
            if tmp <= du {
                du = tmp;
                ge = Some(ei.clone());
            }
        }
        lint.clear();
        if let Some(ge) = ge {
            lint.push(ge);
        }
    }

    /// `TryParameter(OE, V, NE, TolConf)` (`cxx:3638-3681`): when a vertex of
    /// `oe` also lies on `ne` at the same parameter and that parameter is
    /// strictly inside `ne`'s range, bind the vertex to it with
    /// `B.UpdateVertex(V.Oriented(INTERNAL), U, NE, tolerance(NE))`.
    ///
    /// The two arms are not exclusive: the last claim wins (`cxx:3653-3668`),
    /// and `BRep_Tool::Tolerance(NE)` (not the vertex's) is the tolerance
    /// passed on (`cxx:3676`).
    pub fn try_parameter(oe: &Edge, v: &Vertex, ne: &Edge, tol_conf: f64) -> bool {
        let Some(oc) = BRepTool::edge_curve_world(oe) else {
            return false;
        };
        let Some(nc) = BRepTool::edge_curve_world(ne) else {
            return false;
        };
        let (of, ol) = BRepTool::edge_parameters(oe);
        let (nf, nl) = BRepTool::edge_parameters(ne);
        let mut u = 0.0;
        let p = BRepTool::vertex_point(v);
        let mut ok = false;

        if p.distance(&oc.d0(of)) < tol_conf
            && of > nf
            && of < nl
            && p.distance(&nc.d0(of)) < tol_conf
        {
            ok = true;
            u = of;
        }
        if p.distance(&oc.d0(ol)) < tol_conf
            && ol > nf
            && ol < nl
            && p.distance(&nc.d0(ol)) < tol_conf
        {
            ok = true;
            u = ol;
        }
        if ok {
            let vv = v.0.oriented(Orientation::Internal);
            let tol = BRepTool::edge_tolerance(ne);
            let _ = update_vertex_on_edge(&vv, u, &ne.0, tol);
        }
        ok
    }
}
