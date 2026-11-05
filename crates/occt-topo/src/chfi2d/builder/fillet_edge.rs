//! `BuildFilletEdge` (`ChFi2d_Builder.cxx:764-1230`): solve the 2D circle
//! tangent to the two adjacent edges with the requested radius
//! (`Geom2dGcc_Circ2d2TanRad`), pick the solution whose tangency lies in the
//! admissible parameter window, and build the fillet edge over it.
//!
//! The parameters the whole function works with are the **p-curve** parameters
//! of `newFace`: `BRep_Tool::CurveOnSurface(E, newFace, f, l)` for the ranges,
//! `BRep_Tool::Parameter(V, E)` for the vertex parameters. On a planar face the
//! two parameterisations coincide (`ChFi2d` relies on that through
//! `BRepLib::SameParameter`), which is why the C++ compares them directly.

use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::elib::{clib, clib2d};
use occt_core::gp::{GpCirc2d, GpPln, GpPnt2d, GpVec2d};
use occt_core::kernel::geomabs::CurveType;
use occt_core::precision::{ANGULAR, CONFUSION, PCONFUSION};
use occt_geom::{Curve, GeomCircle};
use occt_geom2d::adaptor_curve::{basis2d, type2d};
use occt_geom2d::gcc::{GccEntPosition, Geom2dGccCirc2d2TanRad, Geom2dGccQualifiedCurve};

use crate::abs::Orientation;
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_lib_make_edge::{add_endpoint_vertices, init_3d_vertices};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::chfi2d::ConstructionError;
use crate::shape::{Edge, TopoShape, Vertex};
use crate::topo_tools::{edge_vertices, is_same};

use super::ChFi2dBuilder;

/// The 2D circle solution plus what `BuildFilletEdge` needs to know about the
/// first argument: whether it is a `Geom2d_Circle` (so `Sens1` comes from
/// `gp_Circ2d::IsDirect()`), and its `ElCLib::D1` evaluation.
fn d1_on_argument(
    curve: &dyn occt_geom2d::curve::Curve2d,
    param: f64,
) -> Option<(GpPnt2d, GpVec2d, bool)> {
    match type2d(curve) {
        CurveType::Circle => {
            let circ = basis2d(curve).gp_circ2d()?;
            let (p, v) = clib2d::circle_d1_ax22d(param, circ.position(), circ.radius());
            Some((p, v, circ.is_direct()))
        }
        CurveType::Line => {
            let lin = basis2d(curve).gp_lin2d()?;
            let (p, v) = clib2d::line_d1_ax2d(param, lin.position());
            Some((p, v, true))
        }
        // `cxx:865-887` ends with an unchecked `down_cast<Geom2d_Line>`; a
        // curve that is neither is already refused by `IsLineOrCircle` at
        // `AddFillet` (`cxx:259-263`).
        _ => None,
    }
}

impl ChFi2dBuilder {
    /// `BuildFilletEdge(V, AdjEdge1, AdjEdge2, Radius, NewExtr1, NewExtr2)`
    /// (`cxx:764-1230`).
    ///
    /// Returns the fillet edge and the two new extremities of the adjacent
    /// edges. `None` is one of the OCCT early returns — `TangencyError`
    /// (`cxx:912`), `ComputationError` (`cxx:1021`, `:1096`, `:1114`) — or the
    /// null-vertex / solve-failure arms that throw in OCCT.
    pub(super) fn build_fillet_edge(
        &mut self,
        v: &TopoShape,
        adj_edge1: &Edge,
        adj_edge2: &Edge,
        radius: f64,
    ) -> Option<(Edge, TopoShape, TopoShape)> {
        let face = self.new_face.clone()?;
        let plane = self.ref_plane()?;
        let surface = self.ref_surface()?;

        // `cxx:774-780`: the four endpoint vertices, taken **before** the two
        // edges are forced FORWARD (`cxx:788-790`).
        let (v1, v2) = edge_vertices(adj_edge1);
        let (v3, v4) = edge_vertices(adj_edge2);
        let (v1, v2, v3, v4) = (v1?, v2?, v3?, v4?);

        let mut e1 = adj_edge1.clone();
        let mut e2 = adj_edge2.clone();
        // `cxx:785-790`.
        let oe1 = e1.0.orientation();
        e1.0.set_orientation(Orientation::Forward);
        e2.0.set_orientation(Orientation::Forward);

        // `cxx:806-823`: the parameters of `V` and of the opposite vertex on
        // each edge, in the order that puts the vertex `V` first.
        let o1;
        let (param1, param2) = if is_same(&v1.0, v) {
            o1 = v2.0.orientation();
            (
                BRepTool::parameter_of_vertex(&e1.0, &v1.0),
                BRepTool::parameter_of_vertex(&e1.0, &v2.0),
            )
        } else {
            o1 = v1.0.orientation();
            (
                BRepTool::parameter_of_vertex(&e1.0, &v2.0),
                BRepTool::parameter_of_vertex(&e1.0, &v1.0),
            )
        };
        let (param3, param4) = if is_same(&v3.0, v) {
            (
                BRepTool::parameter_of_vertex(&e2.0, &v3.0),
                BRepTool::parameter_of_vertex(&e2.0, &v4.0),
            )
        } else {
            (
                BRepTool::parameter_of_vertex(&e2.0, &v4.0),
                BRepTool::parameter_of_vertex(&e2.0, &v3.0),
            )
        };

        // `cxx:834-835`: `BRep_Tool::CurveOnSurface(E1, newFace, ufirst1, ulast1)`.
        let (c1, ufirst1, ulast1) = curve_on_surface_range(&e1, &face)?;
        let (c2, ufirst2, ulast2) = curve_on_surface_range(&e2, &face)?;
        let c1_is_circle = type2d(c1.as_ref()) == CurveType::Circle;
        let c2_is_circle = type2d(c2.as_ref()) == CurveType::Circle;

        // `cxx:864-888`: point and tangent at `V` on each argument, plus the
        // direction sense of a circle argument.
        let (p, ve1, sens1) = d1_on_argument(c1.as_ref(), param1)?;
        let (_p2, ve2, sens2) = d1_on_argument(c2.as_ref(), param3)?;

        // `cxx:892-916`: the tangency / downcast-point processing.
        let mut cross = ve1.crossed(&ve2);
        let mut ve3 = ve1;
        let mut ve4 = ve2;
        if ve1.is_parallel(&ve2, ANGULAR) {
            cross = 0.0;
            if param1 < param2 {
                ve3 = ve1.reversed();
            }
            if param3 > param4 {
                ve4 = ve2.reversed();
            }
            if !ve4.is_opposite(&ve3, ANGULAR) {
                // `cxx:910-914`: a true tangency point.
                self.status = ConstructionError::TangencyError;
                return None;
            }
        }

        // `cxx:918-1012`: the qualifier of each argument from the sign of the
        // cross product and the two parameter orders.
        let (qual1, qual2) = if cross < 0.0 {
            let qual1 = if param3 > param4 {
                if sens1 { GccEntPosition::Enclosed } else { GccEntPosition::Outside }
            } else if sens1 {
                GccEntPosition::Outside
            } else {
                GccEntPosition::Enclosed
            };
            let qual2 = if param1 > param2 {
                if sens2 { GccEntPosition::Outside } else { GccEntPosition::Enclosed }
            } else if sens2 {
                GccEntPosition::Enclosed
            } else {
                GccEntPosition::Outside
            };
            (qual1, qual2)
        } else {
            let qual1 = if param3 > param4 {
                if sens1 { GccEntPosition::Outside } else { GccEntPosition::Enclosed }
            } else if sens1 {
                GccEntPosition::Enclosed
            } else {
                GccEntPosition::Outside
            };
            let qual2 = if param1 > param2 {
                if sens2 { GccEntPosition::Enclosed } else { GccEntPosition::Outside }
            } else if sens2 {
                GccEntPosition::Outside
            } else {
                GccEntPosition::Enclosed
            };
            (qual1, qual2)
        };

        // `cxx:1013-1018`.
        let tol = CONFUSION;
        let fillet = Geom2dGccCirc2d2TanRad::new(
            &Geom2dGccQualifiedCurve::new(c1.clone(), qual1),
            &Geom2dGccQualifiedCurve::new(c2.clone(), qual2),
            radius,
            tol,
        );
        // `cxx:1019-1023`.
        if !fillet.is_done() || fillet.nb_solutions() == 0 {
            self.status = ConstructionError::ComputationError;
            return None;
        }

        // `cxx:1024-1071`: keep the solution whose tangency point on the first
        // argument is closest to `V` and inside that argument's span.
        let mut numsol = 1usize;
        let mut dist1 = 1.0e40;
        let mut nsol = 1usize;
        while nsol <= fillet.nb_solutions() {
            // `Fillet.Tangency1(nsol, PU1, PU2, Ptg1)` throws when
            // `IsTheSame1(nsol)`; the port skips that solution.
            if let Some((_pu1, pu2, ptg1)) = fillet.tangency1(nsol) {
                let mut dist = ptg1.distance(&p);
                if !c1_is_circle {
                    // `cxx:1039-1046`.
                    let inside = (pu2 < param1 && pu2 > param2) || (pu2 < param2 && pu2 > param1);
                    if inside && dist < dist1 {
                        numsol = nsol;
                        dist1 = dist;
                    }
                } else if let Some((_ppu1, ppu2, ptg2)) = fillet.tangency2(nsol) {
                    // `cxx:1048-1067`.
                    dist = ptg2.distance(&p);
                    let mut inside =
                        (ppu2 < param3 && ppu2 > param4) || (ppu2 < param4 && ppu2 > param3);
                    // Case of an arc of circle passing on the sewing.
                    if c2_is_circle
                        && ((2.0 * PI < param3 && 2.0 * PI > param4)
                            || (2.0 * PI < param4 && 2.0 * PI > param3))
                    {
                        // case param3 < param4
                        inside =
                            (param3 < ppu2 && ppu2 < 2.0 * PI) || (0.0 <= ppu2 && ppu2 < param4 - 2.0 * PI);
                        // case param4 < param3
                        inside = inside
                            || (param4 < ppu2 && ppu2 < 2.0 * PI)
                            || (0.0 <= ppu2 && ppu2 < param3 - 2.0 * PI);
                    }
                    if inside && dist < dist1 {
                        numsol = nsol;
                        dist1 = dist;
                    }
                }
            }
            nsol += 1;
        }

        // `cxx:1071-1079`.
        let cir: GpCirc2d = fillet.this_solution(numsol);
        let (mut u1, u2, ptg1) = fillet.tangency1(numsol)?;
        let (mut vv1, vv2, ptg2) = fillet.tangency2(numsol)?;

        // `cxx:1080-1098`: the tangency parameter on the first argument must
        // be inside its span (the `>=` / `<=` form of bug 0022695).
        let mut inside = (u2 < param1 && u2 >= param2) || (u2 <= param2 && u2 > param1);
        if c1_is_circle
            && ((2.0 * PI < param1 && 2.0 * PI > param2) || (2.0 * PI < param2 && 2.0 * PI > param1))
        {
            // arc of circle containing the circle origin
            inside = (param1 < u2 && u2 < 2.0 * PI) || (0.0 <= u2 && u2 < param2 - 2.0 * PI);
            inside = inside || (param2 < u2 && u2 < 2.0 * PI) || (0.0 <= u2 && u2 < param1 - 2.0 * PI);
        }
        if !inside {
            self.status = ConstructionError::ComputationError;
            return None;
        }

        // `cxx:1100-1116`: the same check on the second argument.
        inside = (vv2 < param3 && vv2 >= param4) || (vv2 <= param4 && vv2 > param3);
        if c2_is_circle
            && ((2.0 * PI < param3 && 2.0 * PI > param4) || (2.0 * PI < param4 && 2.0 * PI > param3))
        {
            inside = (param3 < vv2 && vv2 < 2.0 * PI) || (0.0 <= vv2 && vv2 < param4 - 2.0 * PI);
            inside = inside || (param4 < vv2 && vv2 < 2.0 * PI) || (0.0 <= vv2 && vv2 < param3 - 2.0 * PI);
        }
        if !inside {
            self.status = ConstructionError::ComputationError;
            return None;
        }

        // `cxx:1118-1142`: the new extremities; a tangency that lands on a
        // span end reuses that endpoint vertex.
        let b = TopoBuilder::new();
        let p1 = surface.d0(ptg1.x(), ptg1.y());
        let p2 = surface.d0(ptg2.x(), ptg2.y());
        let vertex1 = b.make_vertex(p1, tol);
        let mut new_extr1 = vertex1.0.clone();
        if (u2 - ufirst1).abs() <= PCONFUSION {
            new_extr1 = v1.0.clone();
        }
        if (u2 - ulast1).abs() <= PCONFUSION {
            new_extr1 = v2.0.clone();
        }
        let vertex2 = b.make_vertex(p2, tol);
        let mut new_extr2 = vertex2.0.clone();
        if (vv2 - ufirst2).abs() <= PCONFUSION {
            new_extr2 = v3.0.clone();
        }
        if (vv2 - ulast2).abs() <= PCONFUSION {
            new_extr2 = v4.0.clone();
        }

        // `cxx:1146-1165`: the vertex `V` and the opposite vertex of `E1`, in
        // 2D.
        let pntbid = BRepTool::vertex_point(&Vertex(v.clone()));
        let sommet = GpPnt2d::new(pntbid.x(), pntbid.y());
        let sombid = if is_same(&v1.0, v) {
            let q = BRepTool::vertex_point(&v2);
            GpPnt2d::new(q.x(), q.y())
        } else {
            let q = BRepTool::vertex_point(&v1);
            GpPnt2d::new(q.x(), q.y())
        };

        // `cxx:1167-1188`: the fillet's tangent at `U1`, and the tangent of the
        // first argument at `V`.
        let vec = clib2d::circle_d1_ax22d(u1, cir.position(), cir.radius()).1;
        let mut vec1 = match type2d(c1.as_ref()) {
            CurveType::Circle => {
                let circ = basis2d(c1.as_ref()).gp_circ2d()?;
                let par = clib2d::parameter_circ2d(&circ, &ptg1);
                clib2d::circle_d1_ax22d(par, circ.position(), circ.radius()).1
            }
            CurveType::Line => {
                let lin = basis2d(c1.as_ref()).gp_lin2d()?;
                let par = clib2d::parameter_lin2d(&lin, &sommet);
                // `cxx:1181-1185` assigns `gp_Vec2d(sommet - somBid)` and then
                // overwrites it with the `ElCLib::D1` derivative of the same
                // point, so the assignment is dead in OCCT as well.
                let _unused = GpVec2d::new(sommet.x() - sombid.x(), sommet.y() - sombid.y());
                clib2d::line_d1_ax2d(par, lin.position()).1
            }
            _ => return None,
        };
        if oe1 == Orientation::Reversed {
            vec1.reverse();
        }
        let sense = vec1.dot(&vec) > 0.0;

        // `cxx:1191-1195`.
        if u1 > vv1 && u1 > 2.0 * PI {
            clib::adjust_periodic(0.0, 2.0 * PI, CONFUSION, &mut u1, &mut vv1);
        }

        // `cxx:1196-1226`: the edge order follows `O1` (the orientation the
        // vertex `V` carries inside `E1`) against `OE1` (`E1`'s own
        // orientation); `Sense` then flips the edge orientation.
        let forward_pair = (o1 == Orientation::Forward && oe1 == Orientation::Forward)
            || (o1 == Orientation::Reversed && oe1 == Orientation::Reversed);
        let mut fillet_edge = if forward_pair {
            make_edge_from_2d_circle(&cir, &plane, &new_extr1, &new_extr2, u1, vv1)?
        } else {
            make_edge_from_2d_circle(&cir, &plane, &new_extr2, &new_extr1, vv1, u1)?
        };
        if !sense {
            let s1 = fillet_edge.0.orientation();
            fillet_edge = if forward_pair {
                make_edge_from_2d_circle(&cir, &plane, &new_extr2, &new_extr1, vv1, u1)?
            } else {
                make_edge_from_2d_circle(&cir, &plane, &new_extr1, &new_extr2, u1, vv1)?
            };
            // `cxx:1216-1222`.
            fillet_edge.0.set_orientation(if s1 == Orientation::Forward {
                Orientation::Reversed
            } else {
                Orientation::Forward
            });
        }

        // `cxx:1228`: `BRepLib::BuildCurves3d(filletEdge)` is the `has a 3D
        // curve already` early return of `BRepLib::BuildCurve3d`
        // (`BRepLib.cxx:321-325`); the edge above was built with its 3D curve.

        self.status = ConstructionError::IsDone;
        Some((fillet_edge, new_extr1, new_extr2))
    }
}

/// `BRepLib_MakeEdge(circle, refSurf, V1, V2, p1, p2)`
/// (`Builder.cxx:1197`, `:1205`, `:1210`, `:1215`) followed by the
/// `BRepLib::BuildCurves3d` that lifts the p-curve to 3D
/// (`BRepLib.cxx:361-365`, `GeomLib::To3d`).
///
/// The 3D curve `GeomLib::To3d` builds for a `Geom2d_Circle` on a plane is the
/// circle `ElCLib::To3d(plane.position().ax2(), circle2d)`, and its parameters
/// are the p-curve parameters, which is what makes the single-step build
/// equivalent.
fn make_edge_from_2d_circle(
    circle: &GpCirc2d,
    plane: &GpPln,
    v1: &TopoShape,
    v2: &TopoShape,
    p1: f64,
    p2: f64,
) -> Option<Edge> {
    let circ3d = clib::to_3d_circ(&plane.position().ax2(), circle).ok()?;
    let curve: Arc<dyn Curve> = Arc::new(GeomCircle::new(circ3d));
    // `BRepLib_MakeEdge::Init(C, V1, V2, p1, p2)` with both vertices given
    // (`BRepLib_MakeEdge.cxx:603-796`); a periodic curve only adjusts the
    // period (`cxx:628-631`), it never reverses.
    let build = init_3d_vertices(
        &curve,
        Some(&Vertex(v1.clone())),
        Some(&Vertex(v2.clone())),
        p1,
        p2,
    )?;
    let b = TopoBuilder::new();
    let mut edge = b.make_edge(curve, build.first, build.last);
    add_endpoint_vertices(&b, &mut edge, &build);
    Some(edge)
}
