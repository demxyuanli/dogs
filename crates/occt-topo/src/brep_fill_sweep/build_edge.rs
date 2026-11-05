//! The two `BuildEdge` helpers of `BRepFill_Sweep`:
//! `cxx:766-862` (3D curve + pcurve for the four borders of a `BuildFace`) and
//! `cxx:1476-1655` (an iso curve of the swept surface).
//!
//! Translation notes:
//! - `BRep_Builder::MakeEdge(E)` + `UpdateEdge(E, C2d, S, L, Tol)` + `Range` +
//!   `Degenerated(E, true)` is [`degenerated_pcurve_edge`]: the port models the
//!   "null 3D curve but present representation" state as
//!   `EdgeGeom::curve3d_null` (see [`crate::tgeometry::EdgeGeom`]), whose `curve`
//!   carries the `Adaptor3d_CurveOnSurface` the 2D representation implies.
//! - `BRepLib_MakeEdge::Init(C, V1, V2, p1, p2)` is
//!   [`crate::brep_lib_make_edge::init_3d_vertices`] plus
//!   [`crate::brep_lib_make_edge::add_endpoint_vertices`].
//! - `BOPTools_AlgoTools::IsMicroEdge(E, nullCtx)` is
//!   [`crate::pave_ff_se::is_micro_edge`] with `check_splittable = false`.
//! - `Standard_ConstructionError` is the `Err` of
//!   [`build_edge_iso`], carrying OCCT's message.
//! - `BRep_Tool::Pnt` applies the vertex location, hence
//!   [`crate::brep_tool::BRepTool::vertex_point_world`].

use std::sync::Arc;

use occt_core::gp::{GpDir2d, GpPnt2d};
use occt_core::precision::CONFUSION;
use occt_geom::curve::Curve;
use occt_geom::curve_on_surface::CurveOnSurface;
use occt_geom::gcpnts::curve_length_range;
use occt_geom::surface::Surface;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::line::Geom2dLine;

use crate::brep_lib_make_edge::{add_endpoint_vertices, init_3d_vertices};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Vertex};
use crate::tgeometry::{EdgeGeom, GeometryRegistry};
use crate::topo_tools::is_same;

use super::helpers;
use super::same_parameter::check_same_parameter;

/// `BuildEdge(C3d, C2d, S, VF, VL, f, l, Tol3d)` (`cxx:766-862`).
///
/// `None` is the reference's null edge: `BRepLib_MakeEdge` failed at `cxx:839`
/// and the `if (!MkE.IsDone())` arm is empty, so OCCT continues with a null
/// `E`.
pub fn build_edge_3d(
    c3d: &Arc<dyn Curve>,
    c2d: &Arc<dyn Curve2d>,
    s: &Arc<dyn Surface>,
    vf: &Vertex,
    vl: &Vertex,
    f: f64,
    l: f64,
    tol3d: f64,
) -> Option<Edge> {
    let reg = GeometryRegistry::global();
    let b = TopoBuilder::new();

    // `cxx:782-790`.
    let p1 = BRepTool::vertex_point_world(vf);
    let tol1 = BRepTool::vertex_tolerance(vf);
    let p2 = BRepTool::vertex_point_world(vl);
    let tol2 = BRepTool::vertex_tolerance(vl);
    let mut tol = tol1.max(tol2);

    if is_same(&vf.0, &vl.0) || p1.distance(&p2) < tol {
        // Degenerated case (`cxx:794-821`).
        let uv = c2d.d0(f);
        let p = s.d0(uv.x(), uv.y());
        let d = p1.distance(&p);
        if d > tol {
            tol = d;
        }
        let uv = c2d.d0(l);
        let p = s.d0(uv.x(), uv.y());
        let d = p2.distance(&p);
        if d > tol {
            tol = d;
        }

        reg.raise_vertex_tolerance(&vf.0, tol);
        reg.raise_vertex_tolerance(&vl.0, tol);

        return Some(degenerated_pcurve_edge(c2d, s, vf, vl, f, l, tol, &b));
    }

    // `cxx:824-837`.
    let p = c3d.d0(f);
    let d = p1.distance(&p);
    if d > tol1 {
        reg.raise_vertex_tolerance(&vf.0, d);
    }
    let p = c3d.d0(l);
    let d = p2.distance(&p);
    if d > tol2 {
        reg.raise_vertex_tolerance(&vl.0, d);
    }

    // `BRepLib_MakeEdge MkE(C3d, VF, VL, f, l)` (`cxx:839-845`).
    let build = init_3d_vertices(c3d, Some(vf), Some(vl), f, l)?;
    let mut e = b.make_edge(c3d.clone(), build.first, build.last);
    add_endpoint_vertices(&b, &mut e, &build);

    // `B.UpdateEdge(E, C2d, S, Loc, Tol3d)` (`cxx:846-847`).
    reg.update_edge_pcurve_on_surface(&e.0, s, Some(c2d.clone()), tol3d);

    // `cxx:848-858`.
    if crate::pave_ff_se::is_micro_edge(&e, false) {
        reg.raise_vertex_tolerance(&vf.0, p1.distance(&p2));
        return Some(degenerated_pcurve_edge(c2d, s, vf, vf, f, l, tol, &b));
    }

    Some(e)
}

/// `BRep_Builder::MakeEdge(E)` + `UpdateEdge(E, C2d, S, Loc, Tol)` +
/// `Add(E, V1)` + `Add(E, V2)` + `Range(E, f, l)` + `Degenerated(E, true)`
/// (`cxx:815-820` and `cxx:852-858`).
#[allow(clippy::too_many_arguments)]
fn degenerated_pcurve_edge(
    c2d: &Arc<dyn Curve2d>,
    s: &Arc<dyn Surface>,
    vf: &Vertex,
    vl: &Vertex,
    f: f64,
    l: f64,
    tol: f64,
    b: &TopoBuilder,
) -> Edge {
    let reg = GeometryRegistry::global();
    let mut e = Edge::new();
    let cos = CurveOnSurface::new(c2d.clone(), s.clone());
    let mut g = EdgeGeom::new(Arc::new(cos), f, l);
    g.tolerance = tol;
    g.degenerated = true;
    g.curve3d_null = true;
    reg.set_edge(&e.0, g);
    reg.update_edge_pcurve_on_surface(&e.0, s, Some(c2d.clone()), tol);
    b.add(&mut e.0, &vf.0);
    b.add(&mut e.0, &vl.0);
    e
}

/// `BuildEdge(S, isUiso, ValIso, VFirst, VLast, Tol)` (`cxx:1476-1655`): the
/// iso curve of the swept surface as an edge, with its `Geom2d_Line` pcurve.
///
/// `Err` is the `Standard_ConstructionError` of `cxx:1612-1615`, carrying the
/// reference's message.
pub fn build_edge_iso(
    s: &Arc<dyn Surface>,
    is_u_iso: bool,
    val_iso: f64,
    v_first: &Vertex,
    v_last: &Vertex,
    tol: f64,
) -> Result<Edge, String> {
    let reg = GeometryRegistry::global();
    let b = TopoBuilder::new();

    // `cxx:1487-1495`: `S->UIso(ValIso)` / `S->VIso(ValIso)`.
    let iso: Arc<dyn Curve> = if is_u_iso {
        s.u_iso_curve(val_iso)
            .unwrap_or_else(|| occt_geom::surface_ops::iso_v(s.as_ref(), val_iso))
    } else {
        s.v_iso_curve(val_iso)
            .unwrap_or_else(|| occt_geom::surface_ops::iso_u(s.as_ref(), val_iso))
    };

    // `cxx:1497-1514`: singular case.
    let mut sing = false;
    if is_same(&v_first.0, &v_last.0) {
        let mut tol_v = BRepTool::vertex_tolerance(v_first);
        if tol > tol_v {
            tol_v = tol;
        }
        let mid = (iso.first_parameter() + iso.last_parameter()) / 2.0;
        let p = iso.d0(mid);
        if p.distance(&BRepTool::vertex_point_world(v_first)) < tol_v {
            let len = curve_length_range(
                iso.as_ref(),
                iso.first_parameter(),
                iso.last_parameter(),
                tol_v / 4.0,
            );
            sing = len < tol_v;
        }
    }

    let mut e;
    if sing {
        // `cxx:1517-1526`.
        let mut v = v_first.clone();
        e = helpers::null_edge(&mut v);
        reg.set_degenerated(&e.0, true);
    } else {
        // `cxx:1529-1617`.
        let p1 = iso.first_parameter();
        let p2 = iso.last_parameter();
        let pa = iso.d0(p1);
        let pb = iso.d0(p2);

        let t1 = BRepTool::vertex_tolerance(v_first);
        let t2 = BRepTool::vertex_tolerance(v_last);

        let p11 = pa.distance(&BRepTool::vertex_point_world(v_first));
        let p22 = pb.distance(&BRepTool::vertex_point_world(v_last));
        let p12 = pa.distance(&BRepTool::vertex_point_world(v_last));
        let p21 = pb.distance(&BRepTool::vertex_point_world(v_first));

        let fwd = p11 < p12 && p22 < p21;

        if fwd {
            if p11 >= t1 {
                reg.raise_vertex_tolerance(&v_first.0, 1.01 * p11);
            }
            if p22 >= t2 {
                reg.raise_vertex_tolerance(&v_last.0, 1.01 * p22);
            }
        } else {
            if p12 >= t2 {
                reg.raise_vertex_tolerance(&v_last.0, 1.01 * p12);
            }
            if p21 >= t1 {
                reg.raise_vertex_tolerance(&v_first.0, 1.01 * p21);
            }
        }

        // `MkE.Init(Iso, VF, VL, first, last)` (`cxx:1585-1603`).
        let build = if fwd {
            init_3d_vertices(&iso, Some(v_first), Some(v_last), p1, p2)
        } else {
            init_3d_vertices(&iso, Some(v_last), Some(v_first), p1, p2)
        };
        let Some(build) = build else {
            // `throw Standard_ConstructionError("BRepFill_Sweep::BuildEdge")`
            // (`cxx:1612-1615`).
            return Err("BRepFill_Sweep::BuildEdge".to_string());
        };
        e = b.make_edge(iso.clone(), build.first, build.last);
        add_endpoint_vertices(&b, &mut e, &build);
    }

    // `cxx:1619-1638`: associate the 2D line of the iso line.
    let (u_min, u_max) = s.u_range();
    let (v_min, v_max) = s.v_range();
    let iso_first = iso.first_parameter();
    let line = if is_u_iso {
        Geom2dLine::from_pnt_dir(
            GpPnt2d::new(val_iso, v_min - iso_first),
            GpDir2d::new(0.0, 1.0).expect("unit direction"),
        )
    } else {
        Geom2dLine::from_pnt_dir(
            GpPnt2d::new(u_min - iso_first, val_iso),
            GpDir2d::new(1.0, 0.0).expect("unit direction"),
        )
    };
    let _ = (u_max, v_max);
    let line: Arc<dyn Curve2d> = Arc::new(line);
    // `B.UpdateEdge(E, L, S, Loc, Precision::Confusion())` (`cxx:1639`).
    reg.update_edge_pcurve_on_surface(&e.0, s, Some(line.clone()), CONFUSION);
    if sing {
        // `B.Range(E, S, Loc, Iso->FirstParameter(), Iso->LastParameter())`
        // (`cxx:1640-1643`).
        reg.set_pcurve_range_on_surface(&e.0, s, iso.first_parameter(), iso.last_parameter());
    }

    // `cxx:1645-1652`.
    let max_tol = 1.0e-4;
    let (_, the_tol) = check_same_parameter(
        iso.as_ref(),
        iso.first_parameter(),
        iso.last_parameter(),
        line.as_ref(),
        s.as_ref(),
        max_tol,
    );
    // `B.UpdateEdge(E, theTol)` (`cxx:1652`).
    reg.raise_edge_tolerance(&e.0, the_tol);

    Ok(e)
}
