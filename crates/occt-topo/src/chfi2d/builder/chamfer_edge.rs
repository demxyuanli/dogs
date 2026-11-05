//! `BuildChamferEdge` (both overloads), `ComputePoint` and `OrientChamfer`
//! (`ChFi2d_Builder_0.cxx:524-865`).
//!
//! Translation notes:
//! - `BRepAdaptor_Curve(E)` is [`crate::topopebrep::tool::adaptor::CurveAdaptor`]
//!   only for the range; the basis accessors (`c.Line()` / `c.Circle()`) OCCT
//!   calls are `GeomAdaptor_Curve::Line` / `::Circle`, i.e. the **local** basis
//!   curve, while the transform is applied inside `D0`/`D1`. On a `ChFi2d` face
//!   the locations are identity, so the local curve is also what `D0` returns
//!   and [`BRepTool::edge_curve`] is the same handle,
//! - `BRepAdaptor_Curve(E, F)` (`_0.cxx:617`) evaluates the **p-curve** of `E`
//!   on `F`. On a planar `F` the plane parameterisation is an isometry, so the
//!   3D tangent is `Pu * dx + Pv * dy` of the p-curve derivative; that is the
//!   lift used here (`_0.cxx:618-626` only needs the tangent direction),
//! - `GCPnts_AbscissaPoint(cc, D, first)` is
//!   [`occt_geom::gcpnts::abscissa_point`]; `Parameter()` on a `!IsDone`
//!   instance throws in OCCT and is `None` here,
//! - the line of the `Geom2dInt_GInter` helper goes through
//!   `GeomAPI::To2d(newLine->Lin(), refSurf->Pln())`
//!   ([`occt_geom::geom_api::to_2d`]) and the p-curve of `E` on `F`
//!   ([`crate::boptools_2d::curve_on_surface_range`]).

use std::sync::Arc;

use occt_core::elib::clib;
use occt_core::gp::{GpAx1, GpDir, GpLin, GpPnt, GpVec};
use occt_core::precision::{CONFUSION, PINTERSECTION, REAL_SMALL};
use occt_geom::curve::Curve;
use occt_geom::gcpnts::abscissa_point;
use occt_geom::line::GeomLine;
use occt_geom::Surface;
use occt_geom2d::geom2d_int::Geom2dIntGInter;

use crate::abs::Orientation;
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::chfi2d::ConstructionError;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{edge_vertices, is_same};

use super::ChFi2dBuilder;

/// The basis curve behind a `Geom_TrimmedCurve` / `Geom_OffsetCurve`, i.e. what
/// `GeomAdaptor_Curve` keeps in `myCurve` for the `Line()` / `Circle()`
/// accessors.
fn basis_of(curve: &Arc<dyn Curve>) -> Arc<dyn Curve> {
    curve.trim_basis().unwrap_or_else(|| curve.clone())
}

/// `ComputePoint(V, E, D, Param)` (`_0.cxx:702-797`): the point at distance
/// `D` from `V` along `E`, in the direction that leaves the edge.
///
/// `None` is one of OCCT's unchecked arms — a null 3D curve, a zero-length
/// `gp_Vec` normalisation, or a `GCPnts_AbscissaPoint` that is not done.
fn compute_point_on_edge(v: &TopoShape, e: &Edge, d: f64) -> Option<(GpPnt, f64)> {
    // `BRepAdaptor_Curve c(E); first = c.FirstParameter(); ...` (`cxx:705-708`).
    let curve = BRepTool::edge_curve(e)?;
    let (first, last) = BRepTool::edge_parameters(e);
    let (v1, v2) = edge_vertices(e);
    let (v1, v2) = (v1?, v2?);
    let basis = basis_of(&curve);

    // `if (c.GetType() == GeomAbs_Line)` (`cxx:711-736`).
    if let Some(lin) = basis.gp_line() {
        let p1 = BRepTool::vertex_point(&v1);
        let p2 = BRepTool::vertex_point(&v2);
        // `myVec = gp_Vec(p1, p2); myVec.Normalize(); myVec *= D;`
        // (`cxx:718-720`). `Normalize()` raises on a null magnitude
        // (`gp_Vec.hxx:314-321`); [`GpVec::normalize`] answers the magnitude, so
        // that raise is the `REAL_SMALL` test here.
        let mut my_vec = GpVec::from_pnts(&p1, &p2);
        if my_vec.normalize() <= REAL_SMALL {
            return None;
        }
        my_vec.multiply_scalar(d);
        let mut the_point;
        if is_same(&v2.0, v) {
            // `myVec *= -1; thePoint = p2.Translated(myVec);` (`cxx:722-725`).
            my_vec.reverse();
            the_point = p2.translated_vec(&my_vec);
        } else {
            // `thePoint = p1.Translated(myVec);` (`cxx:727-729`).
            the_point = p1.translated_vec(&my_vec);
        }
        // `Param = ElCLib::Parameter(c.Line(), thePoint); c.D0(Param, thePoint);`
        // (`cxx:731-734`, the `OCC20823` re-projection).
        let param = clib::line_parameter(&lin, &the_point);
        the_point = curve.d0(param);
        return Some((the_point, param));
    }

    // `if (c.GetType() == GeomAbs_Circle)` (`cxx:738-767`).
    if let Some(cir) = basis.gp_circ() {
        let radius = cir.radius();
        // `cxx:744-756`: the parameter of `V` first, then of the other vertex.
        let (param1, param2) = if is_same(&v1.0, v) {
            (
                BRepTool::parameter_of_vertex(&e.0, &v1.0),
                BRepTool::parameter_of_vertex(&e.0, &v2.0),
            )
        } else {
            (
                BRepTool::parameter_of_vertex(&e.0, &v2.0),
                BRepTool::parameter_of_vertex(&e.0, &v1.0),
            )
        };
        // `deltaAlpha = D / radius; Param = param1 -/+ deltaAlpha`
        // (`cxx:757-764`).
        let delta_alpha = d / radius;
        let param = if param1 > param2 {
            param1 - delta_alpha
        } else {
            param1 + delta_alpha
        };
        return Some((curve.d0(param), param));
    }

    // `else` (`cxx:769-796`): `GCPnts_AbscissaPoint` on the adaptor curve.
    let p = if is_same(&v1.0, v) {
        BRepTool::vertex_point(&v1)
    } else {
        BRepTool::vertex_point(&v2)
    };
    let param = if p.distance(&curve.d0(first)) <= CONFUSION {
        abscissa_point(curve.as_ref(), d, first).ok()?
    } else {
        abscissa_point(curve.as_ref(), d, last).ok()?
    };
    Some((curve.d0(param), param))
}

/// The `BRepAdaptor_Curve(E, F)` evaluation inside `BuildChamferEdge`'s
/// `(D, Ang)` overload (`_0.cxx:617-626`): the 3D point and tangent of the
/// p-curve `c2d` of `E` on `F` at parameter `u`. See the module notes for the
/// lift of the 2D derivative to 3D.
fn lift_pcurve_point_and_tangent(
    c2d: &dyn occt_geom2d::curve::Curve2d,
    surface: &Arc<dyn Surface>,
    u: f64,
) -> Option<(GpPnt, GpVec)> {
    let (p2d, v2d) = c2d.d1(u);
    let (_p, du, dv) = surface.d1(p2d.x(), p2d.y());
    // `du * dx + dv * dy`: the plane Jacobian applied to the p-curve
    // derivative.
    let tan = du
        .multiplied_scalar(v2d.x())
        .added(&dv.multiplied_scalar(v2d.y()));
    Some((surface.d0(p2d.x(), p2d.y()), tan))
}

/// `ComputePoint(F, L, E, Param)` (`_0.cxx:802-838`): the intersection of the
/// 3D line `L` (projected onto `F`'s plane) with `E`'s p-curve on `F`, taking
/// the point with the smallest parameter on the line.
///
/// `None` where OCCT reads the uninitialised `p2d` of a failed intersection
/// (`cxx:819-836`: `p2d` is only written inside the `NbPoints` loop).
fn compute_point_on_face(line: &GpLin, e: &Edge, f: &Face) -> Option<(GpPnt, f64)> {
    let surface = GeometryRegistry::global().face_surface(&f.0)?;
    let pln = surface.gp_pln()?;
    // `GeomAPI::To2d(L, refSurf->Pln())` (`cxx:807`).
    let line_curve: Arc<dyn Curve> = Arc::new(GeomLine::new(*line));
    let lin2d = occt_geom::geom_api::to_2d(line_curve.as_ref(), &pln)?;
    // `c2d = BRep_Tool::CurveOnSurface(E, F, first, last)` (`cxx:810-811`).
    let (c2d, _first, _last) = curve_on_surface_range(e, f)?;

    // `Geom2dInt_GInter Intersection(adaptorL, adaptorC, PIntersection,
    //  PIntersection)` (`cxx:814-817`).
    let tol = PINTERSECTION;
    let mut intersection = Geom2dIntGInter::new();
    intersection.perform_curves(lin2d.as_ref(), c2d.as_ref(), tol, tol);
    // `if (Intersection.IsDone())` (`cxx:820`).
    if !intersection.is_done() {
        return None;
    }
    // `cxx:821-835`: the point with the smallest `ParamOnFirst`, whose
    // `ParamOnSecond` is the answer's `Param`.
    let mut param_on_line = 1.0e300;
    let mut chosen: Option<(f64, f64, f64)> = None;
    let mut i = 1usize;
    while i <= intersection.nb_points() {
        let ip = intersection.point(i);
        if ip.param_on_first() < param_on_line {
            param_on_line = ip.param_on_first();
            chosen = Some((ip.value().x(), ip.value().y(), ip.param_on_second()));
        }
        i += 1;
    }
    let (x, y, param) = chosen?;
    // `thePoint = Adaptor3dSurface.Value(p2d.X(), p2d.Y())` (`cxx:836-837`).
    Some((surface.d0(x, y), param))
}

/// `OrientChamfer(chamfer, E, V)` (`_0.cxx:842-865`): `FORWARD` when `E`'s
/// orientation and the orientation of `E`'s vertex opposite `V` agree.
fn orient_chamfer(chamfer: &mut Edge, e: &Edge, v: &TopoShape) {
    let orient = e.0.orientation();
    let (v1, v2) = edge_vertices(e);
    let v_orient = match (v1, v2) {
        (Some(first), Some(last)) => {
            if is_same(&first, v) {
                last.orientation()
            } else {
                first.orientation()
            }
        }
        // OCCT reads a null `TopoDS_Vertex` here and compares its orientation
        // with `TopAbs_FORWARD`: a default `TopoDS_Shape` is `FORWARD`
        // (`TopoDS_Shape.cxx:30-36`).
        _ => Orientation::Forward,
    };
    let forward = (orient == Orientation::Forward && v_orient == Orientation::Forward)
        || (orient == Orientation::Reversed && v_orient == Orientation::Reversed);
    chamfer
        .0
        .set_orientation(if forward { Orientation::Forward } else { Orientation::Reversed });
}

impl ChFi2dBuilder {
    /// `BuildChamferEdge(V, AdjEdge1, AdjEdge2, D1, D2, NewExtr1, NewExtr2)`
    /// (`_0.cxx:524-592`).
    pub(super) fn build_chamfer_edge_d1d2(
        &mut self,
        v: &TopoShape,
        adj_edge1: &Edge,
        adj_edge2: &Edge,
        d1: f64,
        d2: f64,
    ) -> Option<(Edge, TopoShape, TopoShape)> {
        // `if (D1 <= 0 || D2 <= 0)` (`cxx:532-536`).
        if d1 <= 0.0 || d2 <= 0.0 {
            self.status = ConstructionError::ParametersError;
            return None;
        }
        // `cxx:538-540`.
        let (p1, param1) = compute_point_on_edge(v, adj_edge1, d1)?;
        let (p2, param2) = compute_point_on_edge(v, adj_edge2, d2)?;

        let tol = CONFUSION;
        let b = TopoBuilder::new();
        // `B.MakeVertex(NewExtr1, p1, tol); B.MakeVertex(NewExtr2, p2, tol);
        //  NewExtr1.Orientation(FORWARD); NewExtr2.Orientation(REVERSED);`
        // (`cxx:542-548`).
        let mut new_extr1 = b.make_vertex(p1, tol).0;
        let mut new_extr2 = b.make_vertex(p2, tol).0;
        new_extr1.set_orientation(Orientation::Forward);
        new_extr2.set_orientation(Orientation::Reversed);

        // `gp_Vec myVec(p1, p2); gp_Dir myDir(myVec);` (`cxx:554-555`); a zero
        // vector makes `gp_Dir` raise `Standard_ConstructionError`.
        let my_vec = GpVec::from_pnts(&p1, &p2);
        let my_dir = GpDir::from_vec(&my_vec).ok()?;
        let new_line = GpLin::from_pnt_dir(p1, my_dir);
        // `B.MakeEdge(chamfer, newLine, tol); B.Range(chamfer, 0., param);`
        // (`cxx:557-558`), `param = ElCLib::Parameter(newLine->Lin(), p2)`.
        let param = clib::line_parameter(&new_line, &p2);
        let curve: Arc<dyn Curve> = Arc::new(GeomLine::new(new_line));
        let mut chamfer = b.make_edge(curve, 0.0, param);
        // `B.Add(chamfer, NewExtr1); B.UpdateVertex(NewExtr1, 0., chamfer, tol);
        //  B.Add(chamfer, NewExtr2); B.UpdateVertex(NewExtr2, param, chamfer,
        //  tol);` (`cxx:559-562`).
        b.add_edge_vertices(
            &mut chamfer,
            &Vertex(new_extr1.clone()),
            &Vertex(new_extr2.clone()),
        );
        let reg = GeometryRegistry::global();
        reg.set_vertex_parameter_on_edge(&chamfer.0, &new_extr1, 0.0);
        reg.set_vertex_parameter_on_edge(&chamfer.0, &new_extr2, param);
        // `cxx:563`.
        orient_chamfer(&mut chamfer, adj_edge1, v);

        // `cxx:566-587`: the orientation each new extremity carries inside the
        // adjacent edge it belongs to.
        set_extremity_orientation(&mut new_extr1, adj_edge1, v);
        set_extremity_orientation(&mut new_extr2, adj_edge2, v);
        // `B.UpdateVertex(NewExtr1, param1, AdjEdge1, tol);
        //  B.UpdateVertex(NewExtr2, param2, AdjEdge2, tol);` (`cxx:588-589`).
        reg.set_vertex_parameter_on_edge(&adj_edge1.0, &new_extr1, param1);
        reg.set_vertex_parameter_on_edge(&adj_edge2.0, &new_extr2, param2);

        // `status = ChFi2d_IsDone; return chamfer;` (`cxx:591-592`).
        self.status = ConstructionError::IsDone;
        Some((chamfer, new_extr1, new_extr2))
    }

    /// `BuildChamferEdge(V, AdjEdge1, D, Ang, AdjEdge2, NewExtr1, NewExtr2)`
    /// (`_0.cxx:596-698`): the chamfer defined by a distance on `AdjEdge1` and
    /// an angle with `AdjEdge2`.
    pub(super) fn build_chamfer_edge_da(
        &mut self,
        v: &TopoShape,
        adj_edge1: &Edge,
        d: f64,
        ang: f64,
        adj_edge2: &Edge,
    ) -> Option<(Edge, TopoShape, TopoShape)> {
        let face = self.new_face.clone()?;
        // `if (D <= 0 || Ang <= 0)` (`cxx:606-610`).
        if d <= 0.0 || ang <= 0.0 {
            self.status = ConstructionError::ParametersError;
            return None;
        }
        // `cxx:611-612`.
        let (p1, param1) = compute_point_on_edge(v, adj_edge1, d)?;
        let p = BRepTool::vertex_point(&Vertex(v.clone()));
        let my_vec = GpVec::from_pnts(&p1, &p);

        // `BRepAdaptor_Curve c(AdjEdge2, refFace); first/last; D1(first)`
        // (`cxx:617-626`): the tangent at the end that touches `V`.
        let (c2d, first, last) = curve_on_surface_range(adj_edge2, &face)?;
        let surface = GeometryRegistry::global().face_surface(&face.0)?;
        let (a_point, mut tan) = lift_pcurve_point_and_tangent(c2d.as_ref(), &surface, first)?;
        if a_point.distance(&p) > CONFUSION {
            let (_p, t) = lift_pcurve_point_and_tangent(c2d.as_ref(), &surface, last)?;
            tan = t;
        }
        // `cxx:629-644`: the tangent orientation follows the orientation of the
        // `AdjEdge2` vertex at `V`.
        let (v1, v2) = edge_vertices(adj_edge2);
        let orient = match (v1, v2) {
            (Some(first_v), Some(last_v)) => {
                if is_same(&first_v, v) {
                    first_v.orientation()
                } else {
                    last_v.orientation()
                }
            }
            _ => Orientation::Forward,
        };
        if orient == Orientation::Reversed {
            tan.reverse();
        }

        // `gp_Ax1 RotAxe(p1, tan ^ myVec); gp_Vec vecLin = myVec.Rotated(RotAxe,
        //  -Ang); gp_Dir myDir(vecLin); newLine = Geom_Line(p1, myDir)`
        // (`cxx:646-649`).
        let axis_dir = GpDir::from_vec(&tan.crossed(&my_vec)).ok()?;
        let rot_axe = GpAx1::new(p1, axis_dir);
        let vec_lin = my_vec.rotated(&rot_axe, -ang);
        let my_dir = GpDir::from_vec(&vec_lin).ok()?;
        let new_line = GpLin::from_pnt_dir(p1, my_dir);

        // `gp_Pnt p2 = ComputePoint(refFace, newLine, AdjEdge2, param2);`
        // (`cxx:651-652`).
        let (p2, param2) = compute_point_on_face(&new_line, adj_edge2, &face)?;

        let tol = CONFUSION;
        let b = TopoBuilder::new();
        // `cxx:654-659`.
        let mut new_extr1 = b.make_vertex(p1, tol).0;
        let mut new_extr2 = b.make_vertex(p2, tol).0;
        new_extr1.set_orientation(Orientation::Forward);
        new_extr2.set_orientation(Orientation::Reversed);

        // `param = ElCLib::Parameter(newLine->Lin(), p2); B.MakeEdge(chamfer,
        //  newLine, tol); B.Range(chamfer, 0., param); B.Add ... UpdateVertex`
        // (`cxx:662-668`).
        let param = clib::line_parameter(&new_line, &p2);
        let curve: Arc<dyn Curve> = Arc::new(GeomLine::new(new_line));
        let mut chamfer = b.make_edge(curve, 0.0, param);
        b.add_edge_vertices(
            &mut chamfer,
            &Vertex(new_extr1.clone()),
            &Vertex(new_extr2.clone()),
        );
        let reg = GeometryRegistry::global();
        reg.set_vertex_parameter_on_edge(&chamfer.0, &new_extr1, 0.0);
        reg.set_vertex_parameter_on_edge(&chamfer.0, &new_extr2, param);
        // `cxx:669`.
        orient_chamfer(&mut chamfer, adj_edge1, v);

        // `cxx:672-693`.
        set_extremity_orientation(&mut new_extr1, adj_edge1, v);
        set_extremity_orientation(&mut new_extr2, adj_edge2, v);
        // `cxx:694-695`.
        reg.set_vertex_parameter_on_edge(&adj_edge1.0, &new_extr1, param1);
        reg.set_vertex_parameter_on_edge(&adj_edge2.0, &new_extr2, param2);

        self.status = ConstructionError::IsDone;
        Some((chamfer, new_extr1, new_extr2))
    }
}

/// `cxx:566-587` / `cxx:672-693`: `NewExtr.Orientation(V1.Orientation())` where
/// `V1` is the `E` vertex at `V`, otherwise the other vertex's orientation.
fn set_extremity_orientation(new_extr: &mut TopoShape, e: &Edge, v: &TopoShape) {
    let (v1, v2) = edge_vertices(e);
    let orient = match (v1, v2) {
        (Some(first), Some(last)) => {
            if is_same(&first, v) {
                first.orientation()
            } else {
                last.orientation()
            }
        }
        // See `orient_chamfer`: a null vertex is a default `TopoDS_Shape`,
        // whose orientation is `FORWARD`.
        _ => Orientation::Forward,
    };
    new_extr.set_orientation(orient);
}
