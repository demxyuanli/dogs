use super::prelude::*;
use super::*;

use occt_core::bnd::BndBox2d;
use occt_core::gp::{
    GpAx3, GpDir2d, GpLin, GpPln, GpPnt, GpPnt2d, GpTrsf2d, GpVec, GpVec2d,
};
use occt_core::precision::{CONFUSION, PCONFUSION, REAL_SMALL};
use occt_geom::{Curve, GeomPlane, Surface};
use occt_geom2d::bspline_curve::Geom2dBSplineCurve;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::line::Geom2dLine;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use crate::abs::Orientation;
use crate::boptools_2d::{curve_on_surface_oriented, replace_pcurve};
use crate::brep_surface::SurfaceKind;
use crate::brep_tool::BRepTool;
use crate::meshing::wire_order::{WireOrder, WireOrderStatus};
use crate::pcurve_full::{
    classify_surface_kind, project_curve_on_surface_perform, reparam_curve2d,
    SurfaceSingularities,
};
use super::shape_analysis_curve::{project_adaptor, Projection};
use super::transfer_params::TransferParametersProj;

mod shifted_wire;
mod pcurve_range;
mod same_parameter;
mod reorder_connected;
mod small;
mod degenerated;
mod lacking;
mod notched;
mod singularity;
mod add_curve3d;

pub use shifted_wire::*;
pub use pcurve_range::*;
pub use same_parameter::*;
pub use reorder_connected::*;
pub use small::*;
pub use degenerated::*;
pub use lacking::*;
pub use notched::*;
pub use add_curve3d::*;


/// `ShapeFix_Wire::MaxTolerance()` at read time. `ShapeFix_Shape::SetMaxTolerance`
/// propagates down to the wire tool (`ShapeFix_Solid.cxx:745-749`,
/// `ShapeFix_Shell.cxx:1709-1713`, `ShapeFix_Face.cxx:173-177`) and
/// `ShapeProcess_OperLibrary.cxx:807` sets it from
/// `FromSTEP.FixShape.MaxTolerance3d` (`STEPControl_Controller.cxx:204-206`),
/// which `STEPControl_ActorRead.cxx:2384` fills with `myMaxTol` =
/// `max(myPrecision, ReadMaxPrecisionVal)` (default `1.0`,
/// `DE_ShapeFixParameters.hxx:32`).
pub(super) const SHAPE_FIX_MAX_TOLERANCE: f64 = 1.0;

/// `ShapeAnalysis::AdjustByPeriod` (`ShapeAnalysis.cxx:48-62`).
pub fn adjust_by_period(val: f64, to_val: f64, period: f64) -> f64 {
    let diff = val - to_val;
    let d = diff.abs();
    let p = period.abs();
    if d <= 0.5 * p {
        return 0.0;
    }
    if p < 1e-100 {
        return diff;
    }
    (if diff > 0.0 { -p } else { p }) * (d / p + 0.5).floor()
}

pub fn first_vertex(edge: &Edge) -> Option<Vertex> {
    let (f, l) = edge_vertices(edge);
    if edge.0.orientation().is_reversed() {
        l
    } else {
        f
    }
}

pub fn last_vertex(edge: &Edge) -> Option<Vertex> {
    let (f, l) = edge_vertices(edge);
    if edge.0.orientation().is_reversed() {
        f
    } else {
        l
    }
}

fn pcurve_at(
    edge: &Edge,
    face: &Face,
) -> Option<(Arc<dyn Curve2d>, f64, f64, GpPnt2d, GpPnt2d)> {
    // `ShapeAnalysis_Edge::PCurve(edge, face, c2d, cf, cl, orient=true)`
    // (`ShapeAnalysis_Edge.cxx:171-188`): takes the face's *surface* plus its
    // location and stays on the surface-level `BRep_Tool::CurveOnSurface`
    // overload. The face orientation is deliberately not composed (see the
    // "COMMENTED BACK" note at `ShapeAnalysis_Edge.cxx:178-184`).
    let (c2, a, b) = curve_on_surface_oriented(edge, face, true)?;
    Some((c2.clone(), a, b, c2.d0(a), c2.d0(b)))
}

fn xy_cross(a: &GpPnt2d, b: &GpPnt2d, x: &GpVec2d) -> f64 {
    let dx = a.x() - b.x();
    let dy = a.y() - b.y();
    dx * x.y() - dy * x.x()
}

fn xy_dot_delta(a: &GpPnt2d, b: &GpPnt2d, x: &GpVec2d) -> f64 {
    (a.x() - b.x()) * x.x() + (a.y() - b.y()) * x.y()
}

fn xy_dot_axis(p: &GpPnt2d, x: &GpVec2d) -> f64 {
    p.x() * x.x() + p.y() * x.y()
}

/// `ShapeAnalysis_Surface::DegeneratedValues(P3d, preci, firstP2d, lastP2d,
/// firstPar, lastPar)` (`ShapeAnalysis_Surface.cxx:373-411`): the singularity
/// whose 3D point is closest to `p3d`, walked in ascending `preci` order over
/// the `ComputeSingularities` arrays (`SortSingularities`, `cxx:1845-1883`,
/// called from `cxx:293`).
///
/// This used to keep its own copy of the cone / torus / sphere / corner arms
/// without the sort, so the `myPreci[i] <= preci` loop condition at `cxx:390`
/// stopped at whichever corner came first instead of the one with the smallest
/// precision. `ShapeAnalysis_Surface::ComputeSingularities` is a pure function
/// of the surface and its bounds, so the port keeps one copy:
/// `pcurve_full::SurfaceSingularities` (whose `min_gap` is this function).
pub(super) fn degenerated_values(
    surf: &dyn Surface,
    p3d: &GpPnt,
    preci: f64,
) -> Option<(GpPnt2d, GpPnt2d, f64, f64)> {
    SurfaceSingularities::compute(surf).min_gap(p3d, preci)
}

/// `ShapeAnalysis_Surface::IsDegenerated(p2d1, p2d2, tol, ratio)`
/// (`ShapeAnalysis_Surface.cxx:547-573`).
fn is_degenerated_2d(surf: &dyn Surface, p1: GpPnt2d, p2: GpPnt2d, tol: f64, ratio: f64) -> bool {
    let a = surf.d0(p1.x(), p1.y());
    let b = surf.d0(p2.x(), p2.y());
    let m = surf.d0(0.5 * (p1.x() + p2.x()), 0.5 * (p1.y() + p2.y()));
    let mut max3d = a.distance(&b).max(m.distance(&a)).max(m.distance(&b));
    if max3d > tol {
        return false;
    }
    // `cxx:562-569`: the parametric deltas are divided by
    // `GeomAdaptor_Surface::UResolution(1.)` / `VResolution(1.)`; a resolution
    // below `Precision::PConfusion()` aborts.
    let ru = occt_geom::approx_same_parameter::u_resolution(surf, 1.0);
    let rv = occt_geom::approx_same_parameter::v_resolution(surf, 1.0);
    if ru < PCONFUSION || rv < PCONFUSION {
        return false;
    }
    let du = (p1.x() - p2.x()).abs() / ru;
    let dv = (p1.y() - p2.y()).abs() / rv;
    max3d *= ratio;
    du * du + dv * dv > max3d * max3d
}

/// `ShapeExtend_WireData::Remove(num)` (`ShapeExtend_WireData.cxx:446-452`).
pub fn remove_wire_edge(wire: &mut Wire, at: usize) {
    degenerated::wire_remove_edge(wire, at);
}

/// One junction of `FixConnected(num, prec)` (`ShapeFix_Wire.cxx:1476-1485`).
pub fn fix_connected_at(wire: &mut Wire, precision: f64, prec: f64, num: usize) -> bool {
    reorder_connected::fix_connected_one(wire, precision, prec, num)
}

/// `ShapeFix_Wire::FixConnected` (`ShapeFix_Wire.cxx:557-569`). `closed` false
/// stops at junction 1 (`aStop = 1`), which is the free-wire arm when
/// `ClosedWireMode` was cleared (`ShapeFix_Shape.cxx:221-223`).
pub fn fix_connected_closed(
    wire: &mut Wire,
    precision: f64,
    max_tolerance: f64,
    closed: bool,
) -> bool {
    let nb = degenerated::wire_edges_nb(wire);
    if nb == 0 {
        return false;
    }
    let stop = if closed { 0 } else { 1 };
    let mut done = false;
    let mut index = nb;
    while index > stop {
        done |= reorder_connected::fix_connected_one(wire, precision, max_tolerance, index);
        index -= 1;
    }
    done
}
