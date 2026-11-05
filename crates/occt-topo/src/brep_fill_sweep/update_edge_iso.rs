//! `BRepFill_Sweep`'s static `UpdateEdge` (`BRepFill_Sweep.cxx:1659-1820`):
//! re-attach the swept surface's iso curve as the edge's pcurve, keeping the
//! edge range, direction and vertex tolerances consistent.
//!
//! Translation notes:
//! - `GeomLib::SameRange` (`cxx:1766`) is
//!   [`occt_geom::geom_lib_same_range_2d::same_range_2d`].
//! - `Couture(E, S, Loc)`'s `Loc` out-parameter is unused by the caller; the
//!   port's [`super::helpers::couture`] returns the pcurve alone.
//! - `BRep_Builder::UpdateEdge(E, nullCurve, 0.)` + `Degenerated(E, true)` +
//!   `Range(E, F2d, L2d)` (`cxx:1757-1761`) maps onto
//!   `remove_edge_curve3d` / `set_degenerated` / `set_edge_range`.
//! - The trailing `B.Range(E, S, L, First, Last)` (`cxx:1818`) writes the range
//!   on the only `BRep_CurveOnSurface` representation this edge can carry at
//!   that point, so it is `set_ranges_all`.
//! - A null iso curve (OCCT would dereference it at `cxx:1688`) is `Err`.

use std::sync::Arc;

use occt_core::gp::{GpDir2d, GpPnt2d, GpVec2d};
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom::curve::Curve;
use occt_geom::gcpnts::curve_length_range;
use occt_geom::geom_lib_same_range_2d::same_range_2d;
use occt_geom::surface::Surface;
use occt_geom2d::line::Geom2dLine;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;
use occt_geom2d::Curve2d;

use crate::brep_tool::BRepTool;
use crate::shape::Edge;
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::edge_vertices;

use super::helpers;
use super::same_parameter::same_parameter;

/// `UpdateEdge(E, S, isUiso, ValIso)` (`cxx:1659-1820`).
pub fn update_edge_iso(
    e: &mut Edge,
    s: &Arc<dyn Surface>,
    is_u_iso: bool,
    val_iso: f64,
) -> Result<(), String> {
    let reg = GeometryRegistry::global();
    let (u_first, u_last) = s.u_range();
    let (v_first, v_last) = s.v_range();

    // `cxx:1671-1681`.
    let iso: Arc<dyn Curve> = if is_u_iso {
        s.u_iso_curve(val_iso)
    } else {
        s.v_iso_curve(val_iso)
    }
    .ok_or_else(|| "UpdateEdge: null iso curve".to_string())?;

    let (v_f, v_l) = edge_vertices(e);

    // `cxx:1683-1695`: the singular case.
    let mut sing = false;
    if let (Some(vf), Some(vl)) = (&v_f, &v_l) {
        if crate::topo_tools::is_same(&vf.0, &vl.0) {
            let tol = BRepTool::vertex_tolerance(vf);
            let pmid = iso.d0((iso.first_parameter() + iso.last_parameter()) / 2.0);
            if pmid.distance(&BRepTool::vertex_point(vf)) < tol {
                sing = curve_length_range(
                    iso.as_ref(),
                    iso.first_parameter(),
                    iso.last_parameter(),
                    tol / 4.0,
                ) < tol;
            }
        }
    }

    // `cxx:1697-1711`: the 2D line and its default range.
    let (line, mut f2d, mut l2d) = if is_u_iso {
        let p = GpPnt2d::new(val_iso, 0.0);
        let v = GpVec2d::new(0.0, 1.0);
        let dir = GpDir2d::from_vec2d(&v).map_err(|e| format!("gp_Dir2d: {e}"))?;
        (
            Geom2dLine::new(occt_core::gp::ax2d::GpAx2d::new(p, dir)),
            v_first,
            v_last,
        )
    } else {
        let p = GpPnt2d::new(0.0, val_iso);
        let v = GpVec2d::new(1.0, 0.0);
        let dir = GpDir2d::from_vec2d(&v).map_err(|e| format!("gp_Dir2d: {e}"))?;
        (
            Geom2dLine::new(occt_core::gp::ax2d::GpAx2d::new(p, dir)),
            u_first,
            u_last,
        )
    };

    let line_arc: Arc<dyn Curve2d> = Arc::new(line);
    let mut cl: Arc<dyn Curve2d> =
        Arc::new(Geom2dTrimmedCurve::new(line_arc.clone(), f2d, l2d));

    // `cxx:1716-1725`.
    let tol = 1.0e-4;
    let (mut first, mut last) = BRepTool::edge_parameters(e);
    let mut reverse = false;

    if !sing_case(&v_f, &v_l) {
        // `cxx:1728-1731`.
        let p2d = cl.d0(f2d);
        let p_on_s = s.value(p2d.x(), p2d.y());
        let (Some(vf), Some(vl)) = (&v_f, &v_l) else {
            return Err("UpdateEdge: missing vertex".into());
        };
        reverse = p_on_s.distance(&BRepTool::vertex_point(vl))
            < p_on_s.distance(&BRepTool::vertex_point(vf));
    } else if !sing {
        // `cxx:1733-1747`: compare the 3D tangent with the surface tangent of
        // the 2D line.
        let Some(c3d) = reg.edge_curve(&e.0) else {
            return Err("UpdateEdge: null 3D curve".into());
        };
        let (_, d_c3d) = c3d.d1(first);
        let (p2d, v2d) = cl.d1(f2d);
        let (_, du, dv) = s.d1(p2d.x(), p2d.y());
        let v3d = du.multiplied_scalar(v2d.x()).added(&dv.multiplied_scalar(v2d.y()));
        reverse = d_c3d.angle(&v3d) > tol;
    }

    if reverse {
        // `cxx:1748-1754`.
        let mut c = Geom2dTrimmedCurve::new(line_arc.clone(), f2d, l2d);
        c.reverse();
        f2d = c.first_parameter();
        l2d = c.last_parameter();
        cl = Arc::new(c);
    }

    // `cxx:1756-1763`.
    if sing {
        reg.remove_edge_curve3d(&e.0);
        reg.set_degenerated(&e.0, true);
        reg.set_edge_range(&e.0, f2d, l2d);
        first = f2d;
        last = l2d;
    }

    // `cxx:1765-1770`.
    if first != f2d || last != l2d {
        let c2d = same_range_2d(PCONFUSION, &cl, f2d, l2d, first, last)?;
        cl = Arc::new(Geom2dTrimmedCurve::new(c2d, first, last));
    }

    // `cxx:1772-1789`: raise the endpoint vertex tolerances.
    let p2d = cl.d0(first);
    let p_on_s = s.value(p2d.x(), p2d.y());
    let (v_first, v_last) = edge_vertices(e);
    if let Some(vf) = &v_first {
        reg.raise_vertex_tolerance(&vf.0, p_on_s.distance(&BRepTool::vertex_point(vf)));
    }
    let p2d = cl.d0(last);
    let p_on_s = s.value(p2d.x(), p2d.y());
    if let Some(vl) = &v_last {
        reg.raise_vertex_tolerance(&vl.0, p_on_s.distance(&BRepTool::vertex_point(vl)));
    }

    // `cxx:1791-1794`: `SameParameter(E, CL, S, Tol, R)` reads `Tol` and
    // answers the error through `R`, so the port needs no separate seed.
    if !sing {
        let (ok, err) = same_parameter(e, &mut cl, s, tol);
        if ok {
            reg.raise_edge_tolerance(&e.0, err);
        }
    }

    // `cxx:1796-1813`.
    if helpers::couture(e, s).is_none() {
        reg.update_edge_pcurve_on_surface(&e.0, s, Some(cl.clone()), CONFUSION);
    } else {
        // Sewing edge: keep both pcurves, ordered by the edge orientation.
        let mut probe = e.clone();
        helpers::oriente(s, &mut probe);
        let curves = if probe.0.orientation() == crate::abs::Orientation::Reversed {
            vec![cl.clone(), pcurve_of(e, s)?]
        } else {
            vec![pcurve_of(e, s)?, cl.clone()]
        };
        reg.set_edge_pcurves(&e.0, Arc::as_ptr(s) as *const () as usize, curves);
    }

    // `cxx:1815-1819`.
    if !reg.same_range(&e.0) {
        reg.set_ranges_all(&e.0, first, last);
    }

    Ok(())
}

/// `Vf.IsSame(Vl)` (`cxx:1685`, `cxx:1727`).
fn sing_case(v_f: &Option<crate::shape::Vertex>, v_l: &Option<crate::shape::Vertex>) -> bool {
    match (v_f, v_l) {
        (Some(a), Some(b)) => crate::topo_tools::is_same(&a.0, &b.0),
        _ => true,
    }
}

/// `PCurve` of the sewing-edge arm: the edge's current pcurve on `s`.
fn pcurve_of(e: &Edge, s: &Arc<dyn Surface>) -> Result<Arc<dyn Curve2d>, String> {
    helpers::couture(e, s).ok_or_else(|| "UpdateEdge: missing pcurve".to_string())
}
