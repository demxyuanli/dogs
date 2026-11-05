//! Surface-closure predicates and the two closure tests of
//! `BRepBuilderAPI_Sewing`: `IsClosedShape` (`cxx:170-213`),
//! `IsClosedByIsos` (`cxx:217-257`), `IsUClosedSurface` (`cxx:260-291`),
//! `IsVClosedSurface` (`cxx:294-325`) and `IsMergedClosed` (`cxx:1525-1615`).

use std::sync::Arc;

use occt_core::bnd::BndBox2d;
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom::gcpnts::curve_length_range;
use occt_geom::{Curve, Surface};
use occt_geom2d::curve::Curve2d;

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape};
use crate::tgeometry::GeometryRegistry;

use super::Sewing;

/// `BRep_Tool::CurveOnSurface(edge, surface, loc, first, last)` for the
/// surface-based lookup: the sewing stages pass the surface they already hold
/// (`cxx:264-267`, `:292-296`), so the representation is found by the surface
/// handle the pcurves are keyed by.
pub(crate) fn pcurve_on_surface(
    edge: &TopoShape,
    surf: &Arc<dyn Surface>,
) -> Option<(Arc<dyn Curve2d>, f64, f64)> {
    let want = Arc::as_ptr(surf) as *const () as usize;
    for (s, c, f, l) in GeometryRegistry::global().edge_pcurve_reps(edge) {
        if Arc::as_ptr(&s) as *const () as usize == want {
            return Some((c, f, l));
        }
    }
    None
}

/// `IsClosedShape(shape, v1, v2)` (`cxx:170-213`): the total length of the
/// edges of `shape` that touch `v1` or `v2`, compared with the vertex
/// distance.
pub(crate) fn is_closed_shape(shape: &TopoShape, v1: &TopoShape, v2: &TopoShape) -> bool {
    let mut tot_length = 0.0;
    for edge in crate::topo_tools::explorer_shapes(shape, ShapeType::Edge) {
        let Some(e) = Edge::wrap(edge) else { continue };
        let (ve1, ve2) = crate::topo_tools::edge_vertices(&e);
        let (Some(ve1), Some(ve2)) = (ve1, ve2) else {
            continue;
        };
        if !crate::topo_tools::is_same(&ve1.0, v1) && !crate::topo_tools::is_same(&ve1.0, v2) {
            continue;
        }
        if BRepTool::is_degenerated(&e) {
            continue;
        }
        let Some(c3d) = BRepTool::edge_curve(&e) else {
            continue;
        };
        let (first, last) = BRepTool::edge_parameters(&e);
        let length = curve_length_range(c3d.as_ref(), first, last, CONFUSION);
        tot_length += length;
        if crate::topo_tools::is_same(&ve2.0, v1) || crate::topo_tools::is_same(&ve2.0, v2) {
            break;
        }
    }
    if tot_length > 0.0 {
        let p1 = BRepTool::vertex_point_world(&crate::shape::Vertex(v1.clone()));
        let p2 = BRepTool::vertex_point_world(&crate::shape::Vertex(v2.clone()));
        return p1.distance(&p2) < tot_length / (1.2 * std::f64::consts::PI);
    }
    false
}

/// `IsClosedByIsos(surf, acrv2d, f2d, l2d, isUIsos)` (`cxx:217-257`).
///
/// `UIso` / `VIso` answer `None` for surfaces whose port has no iso builder;
/// OCCT dereferences the handle it always gets, so a `None` answers "not
/// closed" here (documented deviation, no invented geometry).
fn is_closed_by_isos(
    surf: &Arc<dyn Surface>,
    acrv2d: &Arc<dyn Curve2d>,
    f2d: f64,
    l2d: f64,
    is_u_isos: bool,
) -> bool {
    let p1 = if acrv2d.is_periodic() {
        acrv2d.d0(f2d)
    } else {
        acrv2d.d0(f2d.max(acrv2d.first_parameter()))
    };
    let p2 = if acrv2d.is_periodic() {
        acrv2d.d0(l2d)
    } else {
        acrv2d.d0(l2d.min(acrv2d.last_parameter()))
    };
    let (crv1, crv2) = if is_u_isos {
        (surf.u_iso_curve(p1.x()), surf.u_iso_curve(p2.x()))
    } else {
        (surf.v_iso_curve(p1.y()), surf.v_iso_curve(p2.y()))
    };
    let (Some(crv1), Some(crv2)) = (crv1, crv2) else {
        return false;
    };
    // `cxx:230-241`: the iso curves are 3D; `gp_Xyz` norms are compared.
    let sample = |c: &Arc<dyn Curve>| -> (occt_core::gp::GpPnt, occt_core::gp::GpPnt, occt_core::gp::GpPnt) {
        let af = c.first_parameter();
        let al = c.last_parameter();
        (c.d0(af), c.d0((af + al) * 0.5), c.d0(al))
    };
    let (p11, p1m, p12) = sample(&crv1);
    let (p21, p2m, p22) = sample(&crv2);
    (p11.distance(&p12) < p11.distance(&p1m) - CONFUSION)
        && (p21.distance(&p22) < p21.distance(&p2m) - CONFUSION)
}

impl Sewing {
    /// `IsUClosedSurface(surf, theEdge, theloc)` (`cxx:260-291`).
    pub(crate) fn is_u_closed_surface(
        &self,
        surf: &Arc<dyn Surface>,
        edge: &TopoShape,
    ) -> bool {
        if let Some(basis) = surf.rectangular_trimmed_basis() {
            return self.is_u_closed_surface(&basis, edge);
        }
        if let Some(basis) = surf.offset_basis_surface() {
            return self.is_u_closed_surface(&basis, edge);
        }
        let mut is_closed = surf.is_u_closed();
        if !is_closed {
            if let Some((c2d, f2d, l2d)) = pcurve_on_surface(edge, surf) {
                // `cxx:269` passes `Standard_False`.
                is_closed = is_closed_by_isos(surf, &c2d, f2d, l2d, false);
            }
        }
        is_closed
    }

    /// `IsVClosedSurface(surf, theEdge, theloc)` (`cxx:294-325`).
    pub(crate) fn is_v_closed_surface(
        &self,
        surf: &Arc<dyn Surface>,
        edge: &TopoShape,
    ) -> bool {
        if let Some(basis) = surf.rectangular_trimmed_basis() {
            return self.is_v_closed_surface(&basis, edge);
        }
        if let Some(basis) = surf.offset_basis_surface() {
            return self.is_v_closed_surface(&basis, edge);
        }
        let mut is_closed = surf.is_v_closed();
        if !is_closed {
            if let Some((c2d, f2d, l2d)) = pcurve_on_surface(edge, surf) {
                // `cxx:297` passes `Standard_True`.
                is_closed = is_closed_by_isos(surf, &c2d, f2d, l2d, true);
            }
        }
        is_closed
    }

    /// `IsMergedClosed(Edge1, Edge2, face)` (`cxx:1525-1615`).
    ///
    /// The two pcurve boxes come from the dispatched `BndLib_Add2dCurve::Add`
    /// port ([`crate::bnd_lib_add2d::add_geom2d_range`]) with
    /// `Precision::PConfusion()` (`cxx:1571-1572`), same branch as OCCT.
    pub(crate) fn is_merged_closed(
        &self,
        edge1: &Edge,
        edge2: &Edge,
        face: &Face,
    ) -> bool {
        let Some(surf) = BRepTool::face_surface(face) else {
            return false;
        };
        let is_u_closed = self.is_u_closed_surface(&surf, &edge1.0);
        let is_v_closed = self.is_v_closed_surface(&surf, &edge1.0);
        if !is_u_closed && !is_v_closed {
            return false;
        }
        let Some(c2d1) = crate::boptools_2d::curve_on_surface_range(edge1, face) else {
            return false;
        };
        let Some(c2d2) = crate::boptools_2d::curve_on_surface_range(edge2, face) else {
            return false;
        };
        let b1 = box_2d_curve(&c2d1.0, c2d1.1, c2d1.2);
        let b2 = box_2d_curve(&c2d2.0, c2d2.1, c2d2.2);
        let (c1u_min, c1v_min, c1u_max, c1v_max) = b1;
        let (c2u_min, c2v_min, c2u_max, c2v_max) = b2;
        let du = c1u_max - c1u_min;
        let dv = c1v_max - c1v_min;
        let is_u_long_c1 = dv <= du;
        let is_v_long_c1 = du <= dv;
        let du = c2u_max - c2u_min;
        let dv = c2v_max - c2v_min;
        let is_u_long_c2 = dv <= du;
        let is_v_long_c2 = du <= dv;
        let (su_min, su_max) = surf.u_range();
        let (sv_min, sv_max) = surf.v_range();
        if is_u_closed && is_v_long_c1 && is_v_long_c2 {
            let dist = (c2v_min - c1v_max).max(c1v_min - c2v_max);
            if dist < 0.0 {
                let dist_inner = (c2u_min - c1u_max).max(c1u_min - c2u_max);
                let dist_outer =
                    (su_max - su_min) - (c2u_max - c1u_min).max(c1u_max - c2u_min);
                if dist_outer <= dist_inner {
                    return true;
                }
            }
        }
        if is_v_closed && is_u_long_c1 && is_u_long_c2 {
            let dist = (c2u_min - c1u_max).max(c1u_min - c2u_max);
            if dist < 0.0 {
                let dist_inner = (c2v_min - c1v_max).max(c1v_min - c2v_max);
                let dist_outer =
                    (sv_max - sv_min) - (c2v_max - c1v_min).max(c1v_max - c2v_min);
                if dist_outer <= dist_inner {
                    return true;
                }
            }
        }
        false
    }
}

/// `BndLib_Add2dCurve::Add(C2d, first, last, Precision::PConfusion(), Box)`
/// (`BndLib_Add2dCurve.cxx:76-83`): the dispatched [`crate::bnd_lib_add2d`]
/// port, i.e. pole boxes for Bezier/B-spline pcurves and the 33-point sampling
/// only for the curve types OCCT also falls back on.
fn box_2d_curve(c: &Arc<dyn Curve2d>, first: f64, last: f64) -> (f64, f64, f64, f64) {
    let mut b = BndBox2d::new();
    crate::bnd_lib_add2d::add_geom2d_range(c.as_ref(), first, last, PCONFUSION, &mut b);
    // `cxx:1571-1573` reads `Bnd_Box2d::Get` directly; a void box leaves the
    // four doubles at `0.0` in this port.
    let (xmin, ymin, xmax, ymax) = b.get().unwrap_or((0., 0., 0., 0.));
    (xmin, ymin, xmax, ymax)
}

