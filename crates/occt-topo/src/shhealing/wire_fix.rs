use super::prelude::*;
use super::*;

use occt_core::bnd::BndBox2d;
use occt_core::bspl::curve_tools::reparameterize;
use occt_core::gp::{
    GpAx3, GpDir, GpDir2d, GpLin, GpPln, GpPnt, GpPnt2d, GpTrsf2d, GpVec, GpVec2d,
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
};
use super::shape_analysis_curve::{project_adaptor, Projection};
use super::transfer_params::TransferParametersProj;

/// `ShapeFix_Wire::MaxTolerance()` at read time. `ShapeFix_Shape::SetMaxTolerance`
/// propagates down to the wire tool (`ShapeFix_Solid.cxx:745-749`,
/// `ShapeFix_Shell.cxx:1709-1713`, `ShapeFix_Face.cxx:173-177`) and
/// `ShapeProcess_OperLibrary.cxx:807` sets it from
/// `FromSTEP.FixShape.MaxTolerance3d` (`STEPControl_Controller.cxx:204-206`),
/// which `STEPControl_ActorRead.cxx:2384` fills with `myMaxTol` =
/// `max(myPrecision, ReadMaxPrecisionVal)` (default `1.0`,
/// `DE_ShapeFixParameters.hxx:32`).
const SHAPE_FIX_MAX_TOLERANCE: f64 = 1.0;

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

pub(super) fn first_vertex(edge: &Edge) -> Option<Vertex> {
    let (f, l) = edge_vertices(edge);
    if edge.0.orientation().is_reversed() {
        l
    } else {
        f
    }
}

pub(super) fn last_vertex(edge: &Edge) -> Option<Vertex> {
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

/// Analytic singularities of `ShapeAnalysis_Surface::ComputeSingularities`
/// for cone / torus / sphere (`cxx:201-246`).
pub(super) fn degenerated_values(
    surf: &dyn Surface,
    p3d: &GpPnt,
    preci: f64,
) -> Option<(GpPnt2d, GpPnt2d, f64, f64)> {
    let (su1, su2) = surf.u_range();
    let (sv1, sv2) = surf.v_range();
    let mut hits: Vec<(f64, GpPnt, GpPnt2d, GpPnt2d)> = Vec::new();
    if let Some((radius, alpha)) = surf.cone_ref() {
        let sin_a = alpha.sin();
        if sin_a.abs() > 1e-16 {
            let v_apex = -radius / sin_a;
            let apex = surf.d0(0.0, v_apex);
            hits.push((
                0.0,
                apex,
                GpPnt2d::new(su1, v_apex),
                GpPnt2d::new(su2, v_apex),
            ));
        }
    } else if let Some(tor) = surf.gp_torus() {
        let minor = tor.minor_radius();
        let major = tor.major_radius();
        let ang = (major / minor).min(1.0).acos();
        let pre = (major - minor).max(0.0);
        hits.push((
            pre,
            surf.d0(0.0, std::f64::consts::PI - ang),
            GpPnt2d::new(su1, std::f64::consts::PI - ang),
            GpPnt2d::new(su2, std::f64::consts::PI - ang),
        ));
        if major <= minor {
            hits.push((
                pre,
                surf.d0(0.0, std::f64::consts::PI + ang),
                GpPnt2d::new(su2, std::f64::consts::PI + ang),
                GpPnt2d::new(su1, std::f64::consts::PI + ang),
            ));
        }
    } else if surf.gp_sphere().is_some() {
        hits.push((
            0.0,
            surf.d0(su1, sv2),
            GpPnt2d::new(su2, sv2),
            GpPnt2d::new(su1, sv2),
        ));
        hits.push((
            0.0,
            surf.d0(su1, sv1),
            GpPnt2d::new(su1, sv1),
            GpPnt2d::new(su2, sv1),
        ));
    }
    let mut best: Option<(f64, GpPnt2d, GpPnt2d)> = None;
    for (pre, q, a, b) in hits {
        if pre > preci {
            continue;
        }
        let gap = q.distance(p3d);
        if gap <= preci && best.as_ref().is_none_or(|(g, _, _)| gap < *g) {
            best = Some((gap, a, b));
        }
    }
    best.map(|(_, a, b)| (a, b, su1, su2))
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

/// `ShapeFix_Wire::FixShifted` (`ShapeFix_Wire.cxx:1661-2126`).
pub fn fix_shifted_wire(wire: &Wire, face: &Face) -> bool {
    let Some(surf) = BRepTool::face_surface(face) else {
        return false;
    };
    // `ShapeFix_Wire.cxx:1671-1673`: `surf->IsUClosed(Precision())` /
    // `IsVClosed(Precision())` are the `ShapeAnalysis_Surface` versions, not the
    // bare `Geom_Surface` flags, with `Precision() = ShapeFix_Root::myPrecision`
    // (`ShapeFix_Root.lxx:34`, `ShapeFix_Root.cxx:26` = `Precision::Confusion()`).
    // The extra `Geom_SphericalSurface` arm is OCCT's own, because a sphere is
    // closed in V without being V-periodic.
    let mut u_closed = crate::pcurve_full::sa_is_u_closed(surf.as_ref(), CONFUSION);
    let mut v_closed =
        crate::pcurve_full::sa_is_v_closed(surf.as_ref(), CONFUSION) || surf.gp_sphere().is_some();
    let mut v_range = 1.0;
    let mut v_crv_closed = false;
    if surf.is_surface_of_revolution() {
        // `ShapeFix_Wire.cxx:1682-1707`: for a `Geom_SurfaceOfRevolution` the
        // basis curve decides `vclosed` when it is periodic (CTS18546-2: a 2d
        // contour shifted by 2*PI against a V range of length 2*PI). The
        // `Geom_TrimmedCurve` unwrap is `cxx:1692-1696`; the
        // `Geom_OffsetCurve` unwrap (`cxx:1688-1691`) is UNPORTED, this port
        // has no 3d offset curve type.
        let mut basis = surf.revolution_basis_curve();
        loop {
            let Some(cur) = basis.as_ref() else {
                break;
            };
            if !cur.is_geom_trimmed() {
                break;
            }
            let cur = cur.clone();
            let Some((base, _, _)) = cur.untrimmed_basis() else {
                break;
            };
            basis = Some(base);
        }
        if let Some(basis) = basis {
            if basis.is_periodic() {
                v_closed = true;
                v_range = basis.period();
                v_crv_closed = true;
            }
        }
    }
    if !u_closed && !v_closed {
        return false;
    }
    let (suf, sul) = surf.u_range();
    let (svf, svl) = surf.v_range();
    let su_mid = 0.5 * (suf + sul);
    let sv_mid = 0.5 * (svf + svl);
    let u_range = if u_closed {
        (sul - suf).abs()
    } else {
        f64::MAX
    };
    if !v_crv_closed {
        v_range = if v_closed {
            (svl - svf).abs()
        } else {
            f64::MAX
        };
    }
    if !u_range.is_finite() && !v_range.is_finite() {
        return false;
    }
    let u_tol = 0.2 * u_range;
    let v_tol = 0.2 * v_range;
    let edges: Vec<Edge> = edges_of_wire(wire)
        .into_iter()
        .filter(|e| !(BRepTool::is_degenerated(e) && pcurve_at(e, face).is_none()))
        .collect();
    let nb = edges.len();
    if nb == 0 {
        return false;
    }
    let mut done = false;
    let mut stop = nb;
    let mut ended = nb == 0;
    let mut degstop = false;
    let mut degn2 = 0usize;
    let mut pdeg = GpPnt::new(0.0, 0.0, 0.0);
    let mut n2 = 0usize;
    let mut n1 = nb;
    // `cxx:1764` / `cxx:2015`: first-pass box is junction `p2d1` only.
    let mut box2 = BndBox2d::new();
    while !ended {
        n2 += 1;
        if n2 > nb {
            n2 = 1;
        }
        if n2 == stop {
            ended = true;
        }
        let e1 = &edges[n1 - 1];
        let e2 = &edges[n2 - 1];
        n1 = n2;
        if BRepTool::is_degenerated(e1) || BRepTool::is_degenerated(e2) {
            if !degstop {
                stop = n2;
                degstop = true;
            }
            continue;
        }
        let Some(v) = first_vertex(e2) else {
            continue;
        };
        let p = BRepTool::vertex_point(&v);
        let mut is_deg = 0i32;
        let preci = CONFUSION.max(BRepTool::vertex_tolerance(&v));
        if let Some((deg_p1, deg_p2, _, _)) = degenerated_values(surf.as_ref(), &p, preci) {
            is_deg = if (deg_p1.x() - deg_p2.x()).abs() > (deg_p1.y() - deg_p2.y()).abs() {
                1
            } else {
                2
            };
        }
        const MAX_TOL: f64 = SHAPE_FIX_MAX_TOLERANCE;
        if surf.is_surface_of_revolution() {
            if is_deg == 0 && !v_closed {
                if let Some((_, _, _, _, p_b1)) = pcurve_at(e1, face) {
                    let q1 = GpPnt2d::new(suf, p_b1.y());
                    let q2 = GpPnt2d::new(sul, p_b1.y());
                    if let Some((_, a1, b1, pa, pb)) = pcurve_at(e1, face) {
                        let _ = (a1, b1, pa);
                        if is_degenerated_2d(surf.as_ref(), q1, q2, MAX_TOL, 10.0)
                            && !is_degenerated_2d(surf.as_ref(), pa, pb, MAX_TOL, 10.0)
                        {
                            is_deg = 1;
                        }
                    }
                }
            }
            if is_deg == 0 && !u_closed {
                if let Some((_, _, _, _, p_b1)) = pcurve_at(e1, face) {
                    let q1 = GpPnt2d::new(p_b1.x(), svf);
                    let q2 = GpPnt2d::new(p_b1.x(), svl);
                    if let Some((_, _, _, pa, pb)) = pcurve_at(e1, face) {
                        if is_degenerated_2d(surf.as_ref(), q1, q2, MAX_TOL, 10.0)
                            && !is_degenerated_2d(surf.as_ref(), pa, pb, MAX_TOL, 10.0)
                        {
                            is_deg = 2;
                        }
                    }
                }
            }
        }
        if is_deg != 0 {
            if !degstop {
                stop = n2;
                degstop = true;
            }
            if degn2 == 0 {
                degn2 = n2;
                pdeg = p;
            } else if pdeg.square_distance(&p) < CONFUSION * CONFUSION {
                degn2 = n2;
            } else if try_bi_meridian(
                &edges,
                face,
                degn2,
                n2,
                nb,
                e1,
                e2,
                u_closed,
                if u_closed { u_range } else { v_range },
            ) {
                done = true;
                continue;
            }
        }
        let Some((_, _, _, _, p2d1)) = pcurve_at(e1, face) else {
            continue;
        };
        let Some((c2, _, _, p2d2, _)) = pcurve_at(e2, face) else {
            continue;
        };
        box2.add_point(&p2d1);
        let mut du = 0.0;
        let mut dv = 0.0;
        if u_closed && is_deg != 1 {
            let dx = (p2d2.x() - p2d1.x()).abs();
            if dx > u_range - u_tol {
                du = adjust_by_period(p2d2.x(), p2d1.x(), u_range);
            } else if dx > u_tol && stop == nb {
                stop = n2;
            }
        }
        if v_closed && is_deg != 2 {
            let dy = (p2d2.y() - p2d1.y()).abs();
            if dy > v_range - v_tol {
                dv = adjust_by_period(p2d2.y(), p2d1.y(), v_range);
            } else if dy > v_tol && stop == nb {
                stop = n2;
            }
        }
        if du != 0.0 || dv != 0.0 {
            let mut shift = GpTrsf2d::default();
            shift.set_translation_vec(&GpVec2d::new(du, dv));
            replace_pcurve(e2, face, Arc::from(c2.transformed(&shift)));
            done = true;
        }
    }

    // `cxx:2056-2067`: early-out uses the first-pass `p2d1` box, not `[a,mid]`.
    if box2.is_void() {
        return false;
    }
    let (umin, vmin, umax, vmax) = box2.get().expect("non-void box");
    if (umin + umax - suf - sul).abs() < u_range
        && (vmin + vmax - svf - svl).abs() < v_range
        && !done
    {
        return false;
    }
    // `cxx:2069-2082`: rebuild from `Value(a)` and `Value((a+b)/2)`.
    box2.set_void();
    for e in &edges {
        let Some((c2, a, b, _, _)) = pcurve_at(e, face) else {
            continue;
        };
        box2.add_point(&c2.d0(a));
        box2.add_point(&c2.d0(0.5 * (a + b)));
    }
    let Some((umin, vmin, umax, vmax)) = box2.get() else {
        return done;
    };
    let mut du = 0.0;
    let mut dv = 0.0;
    if u_closed {
        du = adjust_by_period(0.5 * (umin + umax), su_mid, u_range);
    }
    if v_closed {
        dv = adjust_by_period(0.5 * (vmin + vmax), sv_mid, v_range);
    }
    if du == 0.0 && dv == 0.0 {
        return done;
    }
    for e in edges_of_wire(wire) {
        let Some((c2, _, _, _, _)) = pcurve_at(&e, face) else {
            continue;
        };
        let mut shift = GpTrsf2d::default();
        shift.set_translation_vec(&GpVec2d::new(du, dv));
        replace_pcurve(&e, face, Arc::from(c2.transformed(&shift)));
    }
    true
}

fn try_bi_meridian(
    edges: &[Edge],
    face: &Face,
    degn2: usize,
    n2: usize,
    nb: usize,
    e1: &Edge,
    e2: &Edge,
    u_closed: bool,
    u_range: f64,
) -> bool {
    let prev = if degn2 > 1 { degn2 - 1 } else { nb };
    let Some((_, _, _, _, pn1)) = pcurve_at(e1, face) else {
        return false;
    };
    let Some((_, _, _, pn2, _)) = pcurve_at(e2, face) else {
        return false;
    };
    let Some((_, _, _, _, pd1)) = pcurve_at(&edges[prev - 1], face) else {
        return false;
    };
    let Some((_, _, _, pd2, _)) = pcurve_at(&edges[degn2 - 1], face) else {
        return false;
    };
    let (x, period) = if u_closed {
        (GpVec2d::new(1.0, 0.0), u_range)
    } else {
        (GpVec2d::new(0.0, 1.0), u_range)
    };
    let rot1 = xy_cross(&pn1, &pd2, &x);
    let rot2 = xy_cross(&pd1, &pn2, &x);
    let scld = xy_dot_delta(&pd2, &pd1, &x);
    let scln = xy_dot_delta(&pn2, &pn1, &x);
    if !(rot1 * rot2 < -PCONFUSION
        && scld * scln < -PCONFUSION
        && scln.abs() > 0.1 * period
        && scld.abs() > 0.1 * period
        && rot1 * scld > PCONFUSION
        && rot2 * scln > PCONFUSION)
    {
        return false;
    }
    let sign = if rot2 > 0.0 { 1.0 } else { -1.0 };
    let Some((_, a2, b2, _, pb2)) = pcurve_at(e2, face) else {
        return false;
    };
    let Some((_, ax1, bx1, _, _)) = pcurve_at(&edges[prev - 1], face) else {
        return false;
    };
    let Some((cx1, _, _, _, _)) = pcurve_at(&edges[prev - 1], face) else {
        return false;
    };
    let Some((c2d2, _, _, _, _)) = pcurve_at(e2, face) else {
        return false;
    };
    let Some((c2d1, a1, b1, _, _)) = pcurve_at(e1, face) else {
        return false;
    };
    let Some((cx2, ax2, bx2, _, _)) = pcurve_at(&edges[degn2 - 1], face) else {
        return false;
    };
    let deep1 = [
        sign * xy_dot_axis(&pn2, &x),
        sign * xy_dot_axis(&pd1, &x),
        sign * xy_dot_axis(&pb2, &x),
        sign * xy_dot_axis(&cx1.d0(ax1), &x),
        sign * xy_dot_axis(&c2d2.d0(0.5 * (a2 + b2)), &x),
        sign * xy_dot_axis(&cx1.d0(0.5 * (ax1 + bx1)), &x),
    ]
    .into_iter()
    .fold(f64::INFINITY, f64::min);
    let deep2 = [
        sign * xy_dot_axis(&pn1, &x),
        sign * xy_dot_axis(&pd2, &x),
        sign * xy_dot_axis(&c2d1.d0(a1), &x),
        sign * xy_dot_axis(&cx2.d0(bx2), &x),
        sign * xy_dot_axis(&c2d1.d0(0.5 * (a1 + b1)), &x),
        sign * xy_dot_axis(&cx2.d0(0.5 * (ax2 + bx2)), &x),
    ]
    .into_iter()
    .fold(f64::NEG_INFINITY, f64::max);
    let deep = deep2 - deep1;
    let dx = adjust_by_period(deep, 0.5 * (PCONFUSION + period + PCONFUSION), period);
    let scale = if scld > 0.0 { -dx } else { dx };
    let mut k = degn2;
    loop {
        if k > nb {
            k = 1;
        }
        if k == n2 {
            break;
        }
        if let Some((cx, _, _, _, _)) = pcurve_at(&edges[k - 1], face) {
            let mut shift = GpTrsf2d::default();
            shift.set_translation_vec(&GpVec2d::new(x.x() * scale, x.y() * scale));
            replace_pcurve(&edges[k - 1], face, Arc::from(cx.transformed(&shift)));
        }
        k += 1;
    }
    true
}

/// `GeomLib::SameRange` (`GeomLib.cxx:842-922`): remap a 2d curve from
/// `[first_on, last_on]` onto `[req_first, req_last]`. Linear reparam matches
/// the BSpline knot `BSplCLib::Reparametrize` evaluation.
fn is_geom2d_line_curve(c: &dyn Curve2d) -> bool {
    !c.first_parameter().is_finite()
        && !c.last_parameter().is_finite()
        && c.d2(0.0).2.square_magnitude() < 1e-30
}

fn same_range_is_line(c: &dyn Curve2d) -> bool {
    is_geom2d_line_curve(c) || c.is_line()
}

/// `Geom2dConvert::CurveToBSplineCurve` Line arm (`Geom2dConvert.cxx:211-225`)
/// then `BSplCLib::Reparametrize` (`GeomLib.cxx:962-968`).
fn line_trim_to_reparam_bspline(
    c2d: Arc<dyn Curve2d>,
    udeb: f64,
    ufin: f64,
    req_first: f64,
    req_last: f64,
) -> Option<Arc<dyn Curve2d>> {
    if (ufin - udeb).abs() <= PCONFUSION {
        return None;
    }
    // Intermediate TrimmedCurve so Start/End and First/Last match Convert.
    let tc = Geom2dTrimmedCurve::new(c2d, udeb, ufin);
    let k0 = tc.first_parameter();
    let k1 = tc.last_parameter();
    if (k1 - k0).abs() <= PCONFUSION {
        return None;
    }
    let p0 = tc.d0(k0);
    let p1 = tc.d0(k1);
    // Flat knots for Mults=[2,2], Degree=1 unique Knots=[k0,k1].
    let mut knots = vec![k0, k0, k1, k1];
    // `BSplCLib::Reparametrize(RequestedFirst, RequestedLast, Knots)`.
    let u_first = req_first.min(req_last);
    let u_last = req_first.max(req_last);
    reparameterize(&mut knots, k0, k1, u_first, u_last);
    let bs = Geom2dBSplineCurve::new(vec![p0.x(), p1.x()], vec![p0.y(), p1.y()], knots, 1).ok()?;
    Some(Arc::new(bs))
}

/// `GeomLib.cxx:924-969` unequal-span segment for a Line: TrimmedCurve bounds
/// then CurveToBSplineCurve + knot Reparametrize (not a ReparamCurve2d wrapper).
fn same_range_unequal_line(
    c2d: Arc<dyn Curve2d>,
    first_on: f64,
    last_on: f64,
    req_first: f64,
    req_last: f64,
) -> Option<Arc<dyn Curve2d>> {
    let mut a_check: &dyn Curve2d = c2d.as_ref();
    if let Some(basis) = c2d.trimmed_basis() {
        a_check = basis;
    }
    let (udeb, ufin) = if a_check.is_periodic() {
        // `cxx:934-944`: periodic basis — trim to [FirstOn, LastOn] when span
        // is non-degenerate, else the curve's own parameter window.
        if (last_on - first_on).abs() > PCONFUSION {
            (first_on, last_on)
        } else {
            (c2d.first_parameter(), c2d.last_parameter())
        }
    } else {
        // `cxx:946-958`.
        let cf = c2d.first_parameter();
        let cl = c2d.last_parameter();
        let udeb = if cf.is_finite() {
            cf.max(first_on)
        } else {
            first_on
        };
        let ufin = if cl.is_finite() {
            cl.min(last_on)
        } else {
            last_on
        };
        if (ufin - udeb).abs() > PCONFUSION {
            (udeb, ufin)
        } else {
            (cf, cl)
        }
    };
    if !udeb.is_finite() || !ufin.is_finite() {
        return None;
    }
    line_trim_to_reparam_bspline(c2d, udeb, ufin, req_first, req_last)
}

/// `GeomLib::SameRange` (`GeomLib.cxx:842-969`) for equal-span Line/reparam
/// and unequal-span Line Trimmed+Convert.
pub(crate) fn geom_lib_same_range(
    c2d: Arc<dyn Curve2d>,
    first_on: f64,
    last_on: f64,
    req_first: f64,
    req_last: f64,
) -> Arc<dyn Curve2d> {
    if (last_on - req_last).abs() <= PCONFUSION && (first_on - req_first).abs() <= PCONFUSION {
        return c2d;
    }
    if !first_on.is_finite() || !last_on.is_finite() {
        return c2d;
    }
    let equal_span =
        ((last_on - first_on) - (req_last - req_first)).abs() <= PCONFUSION;
    if equal_span {
        // `GeomLib.cxx:864-870` `IsKind(Geom2d_Line)` — exact infinite Line only
        // (not Trimmed/Reparam wrappers that forward `is_line`).
        if is_geom2d_line_curve(c2d.as_ref()) {
            let du = first_on - req_first;
            let (_, dir) = c2d.d1(0.0);
            let mut trsf = GpTrsf2d::identity();
            trsf.set_translation_vec(&GpVec2d::new(dir.x() * du, dir.y() * du));
            return Arc::from(c2d.transformed(&trsf));
        }
        // `GeomLib.cxx:890-900`: TrimmedCurve recurses on the basis then
        // re-trims to RequestedFirst/Last (period-shifted Init2d windows).
        if c2d.trimmed_basis().is_some() {
            let basis = Arc::from(c2d.trimmed_basis().unwrap().clone_dyn());
            let new_basis =
                geom_lib_same_range(basis, first_on, last_on, req_first, req_last);
            return Arc::new(Geom2dTrimmedCurve::new(new_basis, req_first, req_last));
        }
        // Equal-span Circle (`cxx:872-888`) / other: Reparam stand-in until
        // Convert arms are ported (`cxx:908-921`).
        return reparam_curve2d(c2d, req_first, req_last, first_on, last_on);
    }
    // `GeomLib.cxx:924-969`: unequal span — Line uses Trimmed+CurveToBSpline
    // + Reparametrize. Other kinds remain on the Reparam wrapper until Convert
    // is ported (`Geom2dConvert.cxx` non-Line arms).
    if same_range_is_line(c2d.as_ref()) {
        if let Some(bs) =
            same_range_unequal_line(c2d.clone(), first_on, last_on, req_first, req_last)
        {
            return bs;
        }
    }
    reparam_curve2d(c2d, req_first, req_last, first_on, last_on)
}

/// `ShapeAnalysis_Edge::CheckPCurveRange` (`ShapeAnalysis_Edge.cxx:999-1032`).
fn check_pcurve_range(first: f64, last: f64, pc: &dyn Curve2d) -> bool {
    let eps = PCONFUSION;
    let (mut fp, mut lp) = (pc.first_parameter(), pc.last_parameter());
    let mut is_periodic = pc.is_periodic();
    let mut period = if is_periodic { pc.period() } else { f64::MAX };
    if let Some(basis) = pc.trimmed_basis() {
        fp = basis.first_parameter();
        lp = basis.last_parameter();
        is_periodic = basis.is_periodic();
        if is_periodic {
            period = basis.period();
        }
    }
    if is_periodic {
        if last - first > period + eps {
            return false;
        }
    } else if first < fp - eps || last > lp + eps {
        return false;
    }
    true
}

/// `ShapeFix_Edge.cxx:171-325` `TranslatePCurve` for a both-closed seam.
fn translate_pcurve_seam(surf: &dyn Surface, c2d: &dyn Curve2d, tol: f64) -> Arc<dyn Curve2d> {
    let (uf, ul) = surf.u_range();
    let (vf, vl) = surf.v_range();
    let a = c2d.first_parameter();
    let b = c2d.last_parameter();
    let (p0, p1) = if a.is_finite() && b.is_finite() {
        (c2d.d0(a), c2d.d0(b))
    } else {
        let (_, tan) = c2d.d1(0.0);
        return translate_pcurve_dir(surf, c2d, &c2d.d0(0.0), &tan, tol);
    };
    let vec = GpVec2d::new(p1.x() - p0.x(), p1.y() - p0.y());
    let iso_u = GpVec2d::new(0.0, vl - vf);
    let iso_v = GpVec2d::new(ul - uf, 0.0);
    let mut shift = GpTrsf2d::default();
    if vec.is_parallel(&iso_u, tol) {
        if (p0.x() - uf).abs() < (p0.x() - ul).abs() {
            shift.set_translation_vec(&GpVec2d::new(ul - uf, 0.0));
        } else {
            shift.set_translation_vec(&GpVec2d::new(uf - ul, 0.0));
        }
        return Arc::from(c2d.transformed(&shift));
    }
    if vec.is_parallel(&iso_v, tol) {
        if (p0.y() - vf).abs() < (p0.y() - vl).abs() {
            shift.set_translation_vec(&GpVec2d::new(0.0, vl - vf));
        } else {
            shift.set_translation_vec(&GpVec2d::new(0.0, vf - vl));
        }
        return Arc::from(c2d.transformed(&shift));
    }
    Arc::from(c2d.clone_dyn())
}

fn translate_pcurve_dir(
    surf: &dyn Surface,
    c2d: &dyn Curve2d,
    loc: &GpPnt2d,
    dir: &GpVec2d,
    tol: f64,
) -> Arc<dyn Curve2d> {
    let (uf, ul) = surf.u_range();
    let (vf, vl) = surf.v_range();
    let mut shift = GpTrsf2d::default();
    if dir.x().abs() <= tol && dir.y().abs() >= tol {
        if (loc.x() - uf).abs() < (loc.x() - ul).abs() {
            shift.set_translation_vec(&GpVec2d::new(ul - uf, 0.0));
        } else {
            shift.set_translation_vec(&GpVec2d::new(uf - ul, 0.0));
        }
        return Arc::from(c2d.transformed(&shift));
    }
    if dir.x().abs() >= tol && dir.y().abs() <= tol {
        if (loc.y() - vf).abs() < (loc.y() - vl).abs() {
            shift.set_translation_vec(&GpVec2d::new(0.0, vl - vf));
        } else {
            shift.set_translation_vec(&GpVec2d::new(0.0, vf - vl));
        }
        return Arc::from(c2d.transformed(&shift));
    }
    Arc::from(c2d.clone_dyn())
}

/// `ShapeFix_Edge::FixAddPCurve` (`ShapeFix_Edge.cxx:470-614`). `prec` is the
/// caller's tolerance; `ShapeFix_Edge.cxx:499` falls back to the edge tolerance
/// when it is not positive, and that value becomes `myPreci` of the projector
/// (`cxx:533`).
pub(super) fn fix_add_pcurve(
    edge: &Edge,
    face: &Face,
    is_seam: bool,
    prec: f64,
    cache: &mut crate::pcurve_full::ProjectorCache,
) -> bool {
    let face_key = GeometryRegistry::shape_key(&face.0);
    // `ShapeFix_Edge.cxx:475-485`: nothing to do when the edge already carries
    // the pcurve (`ShapeAnalysis_Edge::HasPCurve`) or both seam pcurves
    // (`ShapeAnalysis_Edge::IsSeam`).
    let existing = GeometryRegistry::global().edge_pcurves(&edge.0, face_key);
    if (!is_seam && !existing.is_empty()) || (is_seam && existing.len() >= 2) {
        return false;
    }
    let Some(surf) = BRepTool::face_surface(face) else {
        return false;
    };
    // `ShapeFix_Edge.cxx:487-489`: a pcurve on a plane is not computed. The
    // C++ tests `surf->IsKind(STANDARD_TYPE(Geom_Plane))`, i.e. the *concrete*
    // `Geom_Plane`; a trimmed / offset surface built on a plane is not one and
    // goes to the projector below (`ShapeConstruct_ProjectCurveOnSurface`,
    // whose `projectAnalytic` resolves the wrapper's basis plane).
    if surf.gp_pln().is_some() {
        return false;
    }
    let Some(c3d) = GeometryRegistry::global().edge_curve(&edge.0) else {
        return false;
    };
    let (first, last) = GeometryRegistry::global().edge_parameters(&edge.0);
    // `ShapeFix_Edge.cxx:499`: `preci = (prec > 0. ? prec : BRep_Tool::Tolerance(edge))`.
    let preci = if prec > 0.0 {
        prec
    } else {
        BRepTool::edge_tolerance(edge)
    };
    // `ShapeFix_Edge.cxx:521-531`: `TolFirst` / `TolLast` are the tolerances of
    // the edge's end vertices, `-1` when a vertex is null.
    let tol_first = first_vertex(edge).map_or(-1.0, |v| BRepTool::vertex_tolerance(&v));
    let tol_last = last_vertex(edge).map_or(-1.0, |v| BRepTool::vertex_tolerance(&v));
    let Some(c2d) = project_curve_on_surface_perform(
        c3d.as_ref(),
        surf.as_ref(),
        first,
        last,
        preci,
        tol_first,
        tol_last,
        cache,
    ) else {
        return false;
    };
    if !is_seam {
        GeometryRegistry::global().set_edge_pcurve(&edge.0, face_key, c2d);
        return true;
    }
    let prec = GeometryRegistry::global().edge_tolerance(&edge.0).max(CONFUSION);
    let (uf, ul) = surf.u_range();
    let (vf, vl) = surf.v_range();
    let mut shift = GpTrsf2d::default();
    // `ShapeFix_Edge.cxx:563-579`: which of the two seam pcurves moves by one
    // period is decided by `sas->IsUClosed(prec)` / `IsVClosed(prec)` (the same
    // `preci` as `cxx:499`), falling back to the double-closed
    // `TranslatePCurve` for `sas->IsUClosed()` / `IsVClosed()` with the default
    // `preci = -1`, i.e. `Precision::Confusion()` (`cxx:578`).
    let u_closed = crate::pcurve_full::sa_is_u_closed(surf.as_ref(), preci);
    let v_closed = crate::pcurve_full::sa_is_v_closed(surf.as_ref(), preci);
    let c2d2 = if u_closed && !v_closed {
        shift.set_translation_vec(&GpVec2d::new(ul - uf, 0.0));
        Arc::from(c2d.transformed(&shift))
    } else if v_closed && !u_closed {
        shift.set_translation_vec(&GpVec2d::new(0.0, vl - vf));
        Arc::from(c2d.transformed(&shift))
    } else if crate::pcurve_full::sa_is_u_closed(surf.as_ref(), -1.0)
        && crate::pcurve_full::sa_is_v_closed(surf.as_ref(), -1.0)
    {
        translate_pcurve_seam(surf.as_ref(), c2d.as_ref(), prec)
    } else {
        Arc::from(c2d.clone_dyn())
    };
    GeometryRegistry::global().set_edge_pcurves(&edge.0, face_key, vec![c2d, c2d2]);
    true
}

pub(super) fn is_seam_use(edges: &[Edge], i: usize) -> bool {
    let key = GeometryRegistry::shape_key(&edges[i].0);
    edges
        .iter()
        .filter(|e| GeometryRegistry::shape_key(&e.0) == key)
        .count()
        >= 2
}

/// `StepToTopoDS_TranslateEdgeLoop::CheckPCurves`
/// (`StepToTopoDS_TranslateEdgeLoop.cxx:106-177`): the planar early return,
/// then the per-edge 2D parameter checks (`cxx:126-172`) and the
/// `XSAlgo_ShapeProcessor::CheckPCurve` advanced check (`cxx:175`).
fn check_pcurve_rep_range(wire: &Wire, face: &Face, preci: f64) {
    let Some(surf) = BRepTool::face_surface(face) else {
        return;
    };
    // `cxx:110-114`: planar faces drop STEP pcurves; meshing uses
    // `BRep_Tool::CurveOnPlane` via `make_pcurve_on_face`.
    if classify_surface_kind(surf.as_ref()) == SurfaceKind::Plane {
        let reg = GeometryRegistry::global();
        for e in edges_of_wire(wire) {
            reg.remove_pcurves_on_surface(&e.0, &face.0);
        }
        return;
    }
    let face_key = GeometryRegistry::shape_key(&face.0);
    let reg = GeometryRegistry::global();
    for e in edges_of_wire(wire) {
        let Some((the_pc, mut w1, mut w2)) = curve_on_surface_oriented(&e, face, false) else {
            continue;
        };
        let cf = the_pc.first_parameter();
        let cl = the_pc.last_parameter();
        if w1 == w2 {
            reg.remove_pcurves_on_surface(&e.0, &face.0);
            continue;
        }
        if !the_pc.is_periodic() {
            if w1 < cf {
                w1 = cf;
                reg.set_pcurve_range(&e.0, face_key, w1, w2);
            }
            if w2 > cl {
                w2 = cl;
                reg.set_pcurve_range(&e.0, face_key, w1, w2);
            }
        }
        if w1 > w2 && surf.is_u_periodic() {
            let (u1, u2) = surf.u_range();
            let preci_p = ((w2 - w1).abs() / 2.0).min(PCONFUSION);
            crate::geom_bnd_lib_elclib2d::adjust_periodic(u1, u2, preci_p, &mut w1, &mut w2);
            reg.set_pcurve_range(&e.0, face_key, w1, w2);
        }
    }
    // `cxx:175` ends the `cxx:120-176` loop body with
    // `XSAlgo_ShapeProcessor::CheckPCurve(myEdge, aFace, preci, sbwd->IsSeam(i))`
    // (`XSAlgo_ShapeProcessor.cxx:344-505`), reached from
    // `StepToTopoDS_TranslateEdgeLoop.cxx:875` `CheckPCurves` and therefore from
    // `step/pcurve_ranges.rs` where this function is invoked. That call is the missing
    // piece behind linkrods FACE 35: the gate at `cxx:392-401`
    // (`aDist11/aDist22 > thePrecision` -> `RemovePCurve` + return) is what
    // leaves the 2.172 edge's tolerance at the imported 1.2432146551310009e-07
    // instead of the 3.415648804e-03 our `fix_same_parameter` writes. Measured
    // on that edge (instrumented `check_pcurve`): d1=1.230581e-04,
    // d2=8.621578e-05. `preci` is `StepToTopoDS_Root::Precision()` set from
    // `STEPControl_ActorRead` (`STEPControl_ActorRead.cxx:2370-2384`): the
    // file's `UNCERTAINTY_MEASURE_WITH_UNIT` length measure times the length
    // factor (linkrods.step declares 2.E-005, so preci = 2e-05), else
    // `read.precision.val` (`Interface_StaticStandards.cxx:39`, default
    // 1.e-03); `step/wire_fix::step_precision` parses that.
    //
    // ENABLED (t336). The two blockers recorded here against enabling it are
    // both gone:
    //  - the `Offset` blow-up (708/876 -> 30264/57316) was not caused by this
    //    gate: it was `pcurve_on_reversed` using `Oriented(TopAbs_REVERSED)`
    //    (`TopoDS_Shape::Oriented`, which *sets* the orientation) where OCCT
    //    uses `TopoDS::Edge(theEdge.Reversed())` (`TopoDS_Shape::Reversed()`,
    //    which *flips* FORWARD <-> REVERSED). After that fix the gate on/off
    //    measures `Offset` 712/892 both ways.
    //  - the 53 dropped linkrods pcurves were a symptom of the same wrong
    //    pcurve selection on reversed edges, not of `preci`.
    // Measured on/off over all 16 models (`export_data_obj`, preci from the
    // STEP file): only `linkrods` 3461/5007 -> 3485/5060 (OCCT reference
    // 3494/5078) and `Shape-1` 3327/4304 -> 3291/4230 (OCCT 3343/4336) move;
    // every other model, including all seven hard holds, is unchanged. Cost:
    // the full run goes 44.8 s -> 279 s.
    //
    // The re-projection tail `cxx:403-505` is separately gated by
    // `CHECK_PCURVE_REPROJECT` in `shhealing/pcurve_ranges.rs`; measured inert for all 16
    // models, and still parked for the reason recorded there.
    if T312_WIRE_XSALGO_CHECKPCURVE {
        check_pcurves_xsalgo(wire, face, preci);
    }
}

/// True when `check_pcurve_rep_range` runs the `cxx:175`
/// `XSAlgo_ShapeProcessor::CheckPCurve` call. ENABLED (t336); see the call site
/// for the on/off measurements that replaced the earlier block on it.
const T312_WIRE_XSALGO_CHECKPCURVE: bool = true;

/// `ShapeFix_Edge.cxx:335-464` `TempSameRange`.
/// Remap each COS pcurve whose representation range differs from the 3D
/// range onto that 3D range, then `B.Range` all representations and
/// `B.SameRange(true)`.
fn temp_same_range(edge: &Edge) {
    let reg = GeometryRegistry::global();
    let Some(g) = reg.edge_geom(&edge.0) else {
        return;
    };
    let (current_first, current_last) = (g.first, g.last);
    if !current_first.is_finite() || !current_last.is_finite() {
        return;
    }
    let face_keys: Vec<usize> = g.pcurves.keys().copied().collect();
    for fk in &face_keys {
        let Some((first, last)) = g.pcurve_ranges.get(fk).copied() else {
            continue;
        };
        if (first - current_first).abs() <= PCONFUSION
            && (last - current_last).abs() <= PCONFUSION
        {
            continue;
        }
        let pcs = g.pcurves.get(fk).cloned().unwrap_or_default();
        if pcs.is_empty() {
            continue;
        }
        let new_pcs: Vec<Arc<dyn Curve2d>> = pcs
            .into_iter()
            .map(|pc| geom_lib_same_range(pc, first, last, current_first, current_last))
            .collect();
        if new_pcs.len() == 1 {
            reg.set_edge_pcurve(&edge.0, *fk, new_pcs.into_iter().next().unwrap());
        } else {
            reg.set_edge_pcurves(&edge.0, *fk, new_pcs);
        }
    }
    for fk in face_keys {
        reg.set_pcurve_range(&edge.0, fk, current_first, current_last);
    }
    reg.set_same_range(&edge.0, true);
}

/// `Extrema_LocateExtPC` Newton walk (`int_tools_vertex_line` / ValidateEdge).
/// Returns square distance at the local extremum, or `None` when `!IsDone`.
fn locate_ext_pc_sq<D0, D1>(
    d0: D0,
    d1: D1,
    p: &GpPnt,
    t_seed: f64,
    first: f64,
    last: f64,
) -> Option<f64>
where
    D0: Fn(f64) -> GpPnt,
    D1: Fn(f64) -> (GpPnt, GpVec),
{
    let lo = first.min(last);
    let hi = first.max(last);
    let mut t = t_seed.clamp(lo, hi);
    let mut last_good = d0(t);
    for _ in 0..12 {
        let (q, der) = d1(t);
        let speed2 = der.square_magnitude();
        if speed2 <= 1e-20 {
            return None;
        }
        let num = (q.x() - p.x()) * der.x() + (q.y() - p.y()) * der.y() + (q.z() - p.z()) * der.z();
        let next = t - num / speed2;
        if !next.is_finite() {
            return None;
        }
        let next = next.clamp(lo, hi);
        let qn = d0(next);
        if qn.distance(&last_good) < 1e-10 {
            return Some(qn.square_distance(p));
        }
        t = next;
        last_good = qn;
    }
    Some(last_good.square_distance(p))
}

fn cos_d0(pc: &dyn Curve2d, surf: &dyn Surface, t: f64) -> GpPnt {
    let uv = pc.d0(t);
    surf.d0(uv.x(), uv.y())
}

fn cos_d1(pc: &dyn Curve2d, surf: &dyn Surface, t: f64) -> (GpPnt, GpVec) {
    let (uv, duv) = pc.d1(t);
    let (p, su, sv) = surf.d1(uv.x(), uv.y());
    let tan = su
        .multiplied_scalar(duv.x())
        .added(&sv.multiplied_scalar(duv.y()));
    (p, tan)
}

/// `BRepLib_ValidateEdge::processApprox` (`cxx:103-229`) plus
/// `UpdateTolerance` (`cxx:57-63`).
/// `ShapeAnalysis_Edge::CheckSameParameter` (`cxx:704-838`).
fn check_same_parameter(edge: &Edge, face: &Face) -> f64 {
    if BRepTool::is_degenerated(edge) {
        return 0.0;
    }
    let Some(c3d) = BRepTool::edge_curve(edge) else {
        return 0.0;
    };
    let (first, last) = BRepTool::edge_parameters(edge);
    if !first.is_finite() || !last.is_finite() {
        return 0.0;
    }
    let Some(surf) = BRepTool::face_surface(face) else {
        return 0.0;
    };
    let same_p = BRepTool::same_parameter(edge);
    let face_key = GeometryRegistry::shape_key(&face.0);
    let reg = GeometryRegistry::global();
    let is_plane = classify_surface_kind(surf.as_ref()) == SurfaceKind::Plane;
    // After CheckPCurves RemovePCurves (`cxx:110-114`) a plane has no stored
    // pcurve. Do not same-t leftover STEP pcs (Shape-2 1.322 inflation).
    let mut pcs = if is_plane {
        Vec::new()
    } else {
        reg.edge_pcurves(&edge.0, face_key)
    };
    let (mut cf, mut cl) = reg
        .pcurve_range(&edge.0, face_key)
        .unwrap_or((first, last));
    // `ShapeAnalysis_Edge.cxx:801-826`: `CurveOnPlane` then
    // `Geom2dAdaptor_Curve(aPC, aFirst, aLast)` + `CurveOnSurface` +
    // `BRepLib_ValidateEdge` (`aFirst/aLast` = 3D edge range, cxx:418-419).
    if pcs.is_empty() && is_plane {
        if let Ok(pc) = crate::pcurve::make_pcurve_on_face(edge, face) {
            pcs.push(pc);
            cf = first;
            cl = last;
        }
    }
    let mut max_sq: f64 = 0.0;
    // `SetControlPointsNumber(NbControl - 1)` with default NbControl=23.
    const NB: i32 = 22;
    for pc in &pcs {
        let projection = !same_p
            || (cf - first).abs() > PCONFUSION
            || (cl - last).abs() > PCONFUSION;

        if projection {
            // `BRepLib_ValidateEdge.cxx:141-227`.
            let saved = max_sq;
            let p_ref0 = c3d.d0(first);
            let p_oth0 = cos_d0(pc.as_ref(), surf.as_ref(), cf);
            max_sq = max_sq.max(p_ref0.square_distance(&p_oth0));
            let p_ref1 = c3d.d0(last);
            let p_oth1 = cos_d0(pc.as_ref(), surf.as_ref(), cl);
            max_sq = max_sq.max(p_ref1.square_distance(&p_oth1));
            let mut aborted = false;
            for i in 1..NB {
                let t_ref = ((NB - i) as f64 * first + i as f64 * last) / NB as f64;
                let t_oth = ((NB - i) as f64 * cf + i as f64 * cl) / NB as f64;
                let p_ref = c3d.d0(t_ref);
                let p_oth = cos_d0(pc.as_ref(), surf.as_ref(), t_oth);
                let Some(sq_ref) = locate_ext_pc_sq(
                    |t| c3d.d0(t),
                    |t| c3d.d1(t),
                    &p_oth,
                    t_ref,
                    first,
                    last,
                ) else {
                    aborted = true;
                    break;
                };
                max_sq = max_sq.max(sq_ref);
                let Some(sq_oth) = locate_ext_pc_sq(
                    |t| cos_d0(pc.as_ref(), surf.as_ref(), t),
                    |t| cos_d1(pc.as_ref(), surf.as_ref(), t),
                    &p_ref,
                    t_oth,
                    cf,
                    cl,
                ) else {
                    aborted = true;
                    break;
                };
                max_sq = max_sq.max(sq_oth);
            }
            if aborted {
                // cxx:202-205 / 222-225: !IsDone leaves CalculatedDistance unset
                // for this pair; earlier pairs already contributed to maxdev.
                max_sq = saved;
            }
            continue;
        }
        for i in 0..=NB {
            let t = ((NB - i) as f64 * first + i as f64 * last) / NB as f64;
            let uv = pc.d0(t);
            let sq = c3d.d0(t).square_distance(&surf.d0(uv.x(), uv.y()));
            max_sq = max_sq.max(sq);
        }
    }
    let maxdev = max_sq.sqrt() * 1.00001;
    maxdev
}

/// `ShapeAnalysis_Edge::CheckVertexTolerance` with `checkAll=true`
/// (`cxx:568-667`) then `B.UpdateVertex`.
fn fix_vertex_tolerance_edge(edge: &Edge) {
    let (Some(v1), Some(v2)) = (first_vertex(edge), last_vertex(edge)) else {
        return;
    };
    let p1 = BRepTool::vertex_point(&v1);
    let p2 = BRepTool::vertex_point(&v2);
    let mut t1: f64 = 0.0;
    let mut t2: f64 = 0.0;
    if let Some(c3d) = BRepTool::edge_curve(edge) {
        let (mut a, mut b) = BRepTool::edge_parameters(edge);
        if a.is_finite() && b.is_finite() {
            if edge.0.orientation().is_reversed() {
                std::mem::swap(&mut a, &mut b);
            }
            t1 = p1.square_distance(&c3d.d0(a));
            t2 = p2.square_distance(&c3d.d0(b));
        }
    }
    let reg = GeometryRegistry::global();
    for (surf, pc, mut a, mut b) in reg.edge_pcurve_reps(&edge.0) {
        if !(a.is_finite() && b.is_finite()) {
            continue;
        }
        if edge.0.orientation().is_reversed() {
            std::mem::swap(&mut a, &mut b);
        }
        let uv1 = pc.d0(a);
        let uv2 = pc.d0(b);
        t1 = t1.max(p1.square_distance(&surf.d0(uv1.x(), uv1.y())));
        t2 = t2.max(p2.square_distance(&surf.d0(uv2.x(), uv2.y())));
    }
    let tole = BRepTool::edge_tolerance(edge);
    v1.set_tolerance((1.0000001 * t1.sqrt()).max(tole));
    v2.set_tolerance((1.0000001 * t2.sqrt()).max(tole));
}

/// `ShapeBuild_Edge::CopyReplaceVertices` (`ShapeBuild_Edge.cxx:59-140`): a new
/// edge over the same 3D curve and location, carrying the source vertices
/// (`cxx:99-119`, first FORWARD and last REVERSED) and the source
/// tolerance/flags/ranges/pcurves (`BRep_TEdge::EmptyCopy`, `BRep_TEdge.cxx:104-129`,
/// appends a `Copy()` of every `BRep_GCurve` including its pcurve).
/// `EmptyCopied` also keeps the source orientation (`TopoDS_Shape.hxx:294-302`)
/// and `CopyRanges` (`cxx:125`) re-ranges every representation.
fn copy_replace_vertices(edge: &Edge) -> Edge {
    copy_replace_vertices_with(edge, None, None)
}

/// `ShapeBuild_Edge::CopyReplaceVertices(edge, V1, V2)`
/// (`ShapeBuild_Edge.cxx:59-140`). A null argument keeps the edge's own
/// endpoint (`cxx:65-96`), which is the `TopExp` FORWARD / REVERSED vertex
/// (`cxx:78-93`) and, in this port, `first_vertex` / `last_vertex`.
fn copy_replace_vertices_with(edge: &Edge, v1: Option<&Vertex>, v2: Option<&Vertex>) -> Edge {
    let reg = GeometryRegistry::global();
    let Some(mut geom) = reg.edge_geom(&edge.0) else {
        return edge.clone();
    };
    let builder = TopoBuilder::new();
    let mut copy = builder.make_edge(geom.curve.clone(), geom.first, geom.last);
    let fv = v1.cloned().or_else(|| first_vertex(edge));
    let lv = v2.cloned().or_else(|| last_vertex(edge));
    if let (Some(v1), Some(v2)) = (fv, lv) {
        builder.add_edge_vertices(&mut copy, &v1, &v2);
    }
    // `TopoDS_Shape::EmptyCopied` (`TopoDS_Shape.hxx:294-302`) keeps the source
    // orientation, and `BRep_TEdge::EmptyCopy` (`BRep_TEdge.cxx:104-129`) keeps
    // the tolerance, the `SameParameter` / `SameRange` / `Degenerated` flags and
    // a `Copy()` of every curve representation - the 3D curve and the pcurves.
    // Registering the source `geom` reproduces all of that (see `EdgeGeom`).
    copy.0.set_orientation(edge.0.orientation());
    copy.0.set_location(edge.0.location());
    reg.set_edge(&copy.0, geom);
    copy
}

/// `BRepTools::Compare(V1, V2)` (`BRepTools.cxx:527-545`): same vertex, or a
/// gap no larger than either vertex tolerance.
fn vertices_coincide(a: &Vertex, b: &Vertex) -> bool {
    if is_same(&a.0, &b.0) {
        return true;
    }
    let l = BRepTool::vertex_point(a).distance(&BRepTool::vertex_point(b));
    l <= BRepTool::vertex_tolerance(a) || l <= BRepTool::vertex_tolerance(b)
}

/// `ShapeBuild_Vertex::CombineVertex(V1, V2, tolFactor)`
/// (`ShapeBuild_Vertex.cxx:26-71`).
pub(super) fn combine_vertex(v1: &Vertex, v2: &Vertex, tol_factor: f64) -> Vertex {
    let p1 = BRepTool::vertex_point(v1);
    let p2 = BRepTool::vertex_point(v2);
    let tol1 = BRepTool::vertex_tolerance(v1);
    let tol2 = BRepTool::vertex_tolerance(v2);
    let dist = p1.distance(&p2);
    let (pos, tol) = if dist + tol2 <= tol1 {
        (p1, tol1)
    } else if dist + tol1 <= tol2 {
        (p2, tol2)
    } else {
        let tol = 0.5 * (dist + tol1 + tol2);
        // `ShapeBuild_Vertex.cxx:64`: guarded against a zero distance.
        let s = if dist > 0.0 { (tol2 - tol1) / dist } else { 0.0 };
        (
            GpPnt::new(
                0.5 * ((1.0 - s) * p1.x() + (1.0 + s) * p2.x()),
                0.5 * ((1.0 - s) * p1.y() + (1.0 + s) * p2.y()),
                0.5 * ((1.0 - s) * p1.z() + (1.0 + s) * p2.z()),
            ),
            tol,
        )
    };
    TopoBuilder::new().make_vertex(pos, tol_factor * tol)
}

/// `ShapeBuild_Edge::CopyPCurves` (`ShapeBuild_Edge.cxx:360-413`): copy every
/// `CurveOnSurface` representation of `from` onto `to` (replacing the
/// representation on the same surface, `cxx:376-398`) together with its
/// pcurve `Copy()` and range (`cxx:400-411`).
fn copy_pcurves(to: &Edge, from: &Edge) {
    let reg = GeometryRegistry::global();
    let Some(geom) = reg.edge_geom(&from.0) else {
        return;
    };
    for (face_key, pcurves) in &geom.pcurves {
        let copies: Vec<Arc<dyn Curve2d>> =
            pcurves.iter().map(|pc| Arc::from(pc.clone_dyn())).collect();
        reg.set_edge_pcurves(&to.0, *face_key, copies);
        if let Some(&(first, last)) = geom.pcurve_ranges.get(face_key) {
            reg.set_pcurve_range(&to.0, *face_key, first, last);
        }
    }
}

/// `ShapeFix_Edge::FixSameParameter` (`ShapeFix_Edge.cxx:798-935`).
///
/// Both arms are ported: `wasSP` (`cxx:872-928`) and `!wasSP`
/// (`cxx:839-858` copy + `BRepLib::SameParameter(copy, tol)` with its
/// `Approx_SameParameter` re-projection at `BRepLib.cxx:1631`, then the
/// pick-best block `cxx:889-912` with `ShapeBuild_Edge::CopyPCurves`).
///
/// UNPORTED, inside the `!wasSP` arm: the cxx widens the deviation check to
/// *all* pcurves by passing an empty face (`cxx:877-881`) and
/// `BRepLib::SameParameter(edge, tol)` walks every `CurveOnSurface`
/// representation (`BRepLib.cxx:1263-1291`); our `check_same_parameter` and
/// `same_parameter_edge` are face based, so they cover the call-site face only.
///
/// Also unmodelled in both arms: the `ShapeExtend` status writes of the method
/// (`cxx:853-855` FAIL2, `cxx:883-886` FAIL1, `cxx:894` DONE3, `cxx:902-905`,
/// `cxx:910` DONE5, `cxx:931-934` DONE2); this port has no `myStatus`
/// accumulator.
pub(super) fn fix_same_parameter(edge: &Edge, face: &Face) {
    let reg = GeometryRegistry::global();
    let same_range = reg.edge_geom(&edge.0).map(|g| g.same_range).unwrap_or(true);
    if reg.is_degenerated_edge(&edge.0) {
        if !same_range {
            temp_same_range(edge);
        }
        reg.set_same_parameter(&edge.0, true);
        return;
    }
    let tol = BRepTool::edge_tolerance(edge);
    let (first_3d, last_3d) = BRepTool::edge_parameters(edge);
    // `cxx:820-823`: `V1`/`V2` and their tolerances are read *before*
    // `BRepLib::SameParameter` runs. The copy made by the `!wasSP` arm carries
    // these same two vertices (`ShapeBuild_Edge.cxx:105-115`), and that call
    // widens them (`BRepLib.cxx:1220-1232`), so the restore below needs the
    // values captured here.
    let v1 = first_vertex(edge);
    let v2 = last_vertex(edge);
    let tol_fv = v1.as_ref().map(BRepTool::vertex_tolerance).unwrap_or(0.0);
    let tol_lv = v2.as_ref().map(BRepTool::vertex_tolerance).unwrap_or(0.0);
    // `cxx:826`: `bool wasSP = BRep_Tool::SameParameter(edge), SP = false`.
    let was_sp = BRepTool::same_parameter(edge);
    let mut copy_edge: Option<Edge> = None;
    let mut sp = false;
    if !same_range {
        temp_same_range(edge);
    }
    if !was_sp {
        // `cxx:839-858`: heal on a copy of the edge, keeping the original.
        let copy = copy_replace_vertices(edge); // cxx:841
        copy_pcurves(&copy, edge); // cxx:841 `Copy(edge, false)`
        reg.set_same_parameter(&copy.0, false); // cxx:842
        // `cxx:844-849`: `Copy` may shift a periodic 3D range; enforce the
        // original range so the copy's pcurves are not rewritten under it.
        if let Some(mut geom) = reg.edge_geom(&copy.0) {
            geom.first = first_3d;
            geom.last = last_3d;
            reg.set_edge(&copy.0, geom);
        }
        // `cxx:850`: the two argument overload passes the default
        // `tolerance = 0` (`ShapeFix_Edge.hxx:193`), so `tol` is selected.
        crate::brep_lib_same_parameter::same_parameter_edge_inplace(&copy, face, tol);
        sp = BRepTool::same_parameter(&copy); // cxx:851
        copy_edge = Some(copy);
    }
    reg.set_same_parameter(&edge.0, true); // cxx:872
    // `cxx:883`: deviation on the original pcurves.
    let mut maxdev = check_same_parameter(edge, face);
    if sp {
        if let Some(copy) = copy_edge.as_ref() {
            // `cxx:889-912`: compare and select the best variant.
            let mut brl_tol = BRepTool::edge_tolerance(copy); // cxx:891
            let brl_dev = check_same_parameter(copy, face); // cxx:892
            if brl_tol < brl_dev {
                brl_tol = brl_dev; // cxx:895-898
            }
            if brl_tol < maxdev {
                copy_pcurves(edge, copy); // cxx:906
                maxdev = brl_tol; // cxx:907
                reg.set_edge_tolerance(&edge.0, brl_tol); // cxx:908
            }
        }
    }
    // `cxx:915-922`: restore the vertex tolerances captured at `cxx:820-823`,
    // because `BRepLib::SameParameter` may have widened them through the copy.
    // `ShapeFix_ShapeTolerance::SetTolerance` writes the value *exactly*
    // (`ShapeFix_ShapeTolerance.cxx:142-164`) and ignores a non-positive
    // precision (`cxx:149-152`).
    for (v, tol_v) in [(v1.as_ref(), tol_fv), (v2.as_ref(), tol_lv)] {
        if let Some(v) = v {
            let value = maxdev.max(tol_v);
            if value > 0.0 {
                v.set_tolerance(value);
            }
        }
    }
    if maxdev > tol {
        reg.set_edge_tolerance(&edge.0, maxdev);
        fix_vertex_tolerance_edge(edge);
    }
}

/// `ShapeFix_Wire::FixReorder` (`cxx:487-534`) then `FixReorder(sawo)`
/// (`cxx:1351-1399`). `Perform` calls `FixReorder()` with `theModeBoth=false`,
/// so CheckOrder is 3D (`cxx:504`) for every surface type. Status 0 leaves the
/// wire unchanged.
///
/// Returns the `ReorderOK` flag of `Perform` (`cxx:317-325`): false only when
/// the checked order could not be applied (`ShapeExtend_FAIL`, which then
/// disables FixShifted at `cxx:368-374`). A wire that is already ordered
/// returns true (`FixReorder(sawo)` returns false for status 0 without setting
/// FAIL).
pub fn fix_reorder_wire(wire: &Wire, face: &Face) -> bool {
    let stored: Vec<Edge> = edges_of_wire(wire)
        .into_iter()
        .filter(|e| {
            let o = e.0.orientation();
            o == Orientation::Forward || o == Orientation::Reversed
        })
        .collect();
    if stored.len() < 2 {
        return true;
    }
    let mut order = WireOrder::new();
    for e in &stored {
        // `ShapeAnalysis_Wire::CheckOrder` (`ShapeAnalysis_Wire.cxx:617-621`)
        // sets FAIL2 and leaves the order flag unset, so `Perform` sees
        // `ReorderOK == true` and does not even attempt FixReorder.
        let (Some(v1), Some(v2)) = (first_vertex(e), last_vertex(e)) else {
            return true;
        };
        order.add_edge_xyz(BRepTool::vertex_point(&v1), BRepTool::vertex_point(&v2));
    }
    order.perform();
    if order.status() == WireOrderStatus::Same {
        return true;
    }
    if order.nb_edges() != stored.len() {
        return false;
    }
    let mut new_edges = Vec::with_capacity(stored.len());
    for i in 1..=stored.len() {
        let signed = order.ordered(i);
        if signed == 0 {
            return false;
        }
        let mut e = stored[signed.unsigned_abs() as usize - 1].clone();
        if signed < 0 {
            e.0.reverse();
        }
        new_edges.push(e);
    }
    let extras: Vec<TopoShape> = {
        let ts = wire.0.tshape.read().expect("poisoned TShape lock");
        ts.children
            .iter()
            .filter(|c| {
                !c.is_edge()
                    || (c.orientation() != Orientation::Forward
                        && c.orientation() != Orientation::Reversed)
            })
            .cloned()
            .collect()
    };
    // `ShapeExtend_WireData` stores the edges with the wire orientation already
    // composed (`ShapeExtend_WireData::Init`, `ShapeExtend_WireData.cxx:80-120`)
    // and `ShapeExtend_WireData::Wire()` rebuilds a fresh FORWARD wire from that
    // list (`ShapeExtend_WireData.cxx:651-685`). Our children stay under the
    // original wire, so undo the wire's own orientation before storing, keeping
    // the composed view identical to OCCT's WireData.
    let wire_reversed = wire.0.orientation() == Orientation::Reversed;
    {
        let mut ts = wire.0.tshape.write().expect("poisoned TShape lock");
        ts.children.clear();
        for e in &new_edges {
            let mut child = e.0.clone();
            if wire_reversed {
                child.reverse();
            }
            ts.children.push(child);
        }
        ts.children.extend(extras);
    }
    true
}

/// `ShapeAnalysis_Wire::CheckConnected(num, prec)` status
/// (`ShapeAnalysis_Wire.cxx:693-762`). The three `ShapeExtend_DONE` cases are
/// the only ones `FixConnected` acts on; `FAIL1` / `FAIL2` both return false
/// before any fix (`cxx:1487-1493`).
#[derive(Clone, Copy, PartialEq, Eq)]
enum ConnectedCheck {
    None,
    /// `cxx:727`: `myMin3d <= gp::Resolution()`, "absolutely confused".
    Done1,
    /// `cxx:730`: `myMin3d <= myPrecision`.
    Done2,
    /// `cxx:733`: `myMin3d <= prec`.
    Done3,
}

/// `ShapeAnalysis_Wire::CheckConnected(num, prec)`
/// (`ShapeAnalysis_Wire.cxx:693-762`). `my_precision` is the analyzer's
/// `myPrecision` (`ShapeAnalysis_Wire.cxx:203`), the value
/// `ShapeFix_Wire::Perform` was driven with.
fn check_connected(wire: &Wire, my_precision: f64, prec: f64, num: usize) -> ConnectedCheck {
    let nb = wire_edges_nb(wire);
    if nb < 1 {
        return ConnectedCheck::None;
    }
    let n2 = if num > 0 { num } else { nb };
    let n1 = if n2 > 1 { n2 - 1 } else { nb };
    let edges = edges_of_wire(wire);
    let (Some(v1), Some(v2)) = (last_vertex(&edges[n1 - 1]), first_vertex(&edges[n2 - 1])) else {
        return ConnectedCheck::None;
    };
    if is_same(&v1.0, &v2.0) {
        return ConnectedCheck::None;
    }
    let p1 = BRepTool::vertex_point(&v1);
    let p2 = BRepTool::vertex_point(&v2);
    let mut min3d = p1.distance(&p2);
    if min3d <= occt_core::precision::REAL_SMALL {
        return ConnectedCheck::Done1;
    }
    if min3d <= my_precision {
        return ConnectedCheck::Done2;
    }
    if min3d <= prec {
        return ConnectedCheck::Done3;
    }
    // `cxx:737-758`: the reverse test on the last vertex of E2.
    if n1 == n2 {
        return ConnectedCheck::None;
    }
    let Some(v2b) = last_vertex(&edges[n2 - 1]) else {
        return ConnectedCheck::None;
    };
    let p2b = BRepTool::vertex_point(&v2b);
    let dist = p1.distance(&p2b);
    if dist > my_precision {
        return ConnectedCheck::None;
    }
    min3d = dist;
    let _ = min3d;
    ConnectedCheck::None
}

/// `ShapeFix_Wire::FixConnected(num, prec, theUpdateWire)`
/// (`ShapeFix_Wire.cxx:1476-1608`), at read time (`theUpdateWire = false`,
/// `cxx:569`). `Context()` is not modelled: the two `CopyReplaceVertices`
/// writes to `WireData` (`cxx:1585-1601`) are applied directly to the wire's
/// child list, the same convention `fix_degenerated` uses.
fn fix_connected_one(wire: &mut Wire, my_precision: f64, prec: f64, num: usize) -> bool {
    let check = check_connected(wire, my_precision, prec, num);
    if check == ConnectedCheck::None {
        return false;
    }
    let nb = wire_edges_nb(wire);
    if nb == 0 {
        return false;
    }
    let n2 = if num > 0 { num } else { nb };
    let n1 = if n2 > 1 { n2 - 1 } else { nb };
    let edges = edges_of_wire(wire);
    let e1 = edges[n1 - 1].clone();
    let e2 = edges[n2 - 1].clone();
    let v1 = last_vertex(&e1);
    let v2 = first_vertex(&e2);
    let (Some(v1), Some(v2)) = (v1, v2) else {
        return false;
    };
    // `cxx:1505-1542`: the merged vertex.
    let combined: Vertex = if check == ConnectedCheck::Done1 {
        let v2_is_last = last_vertex(&e2)
            .map(|l| is_same(&v2.0, &l.0))
            .unwrap_or(false);
        if v2_is_last {
            v2.clone()
        } else {
            v1.clone()
        }
    } else {
        combine_vertex(&v1, &v2, 1.0001)
    };

    // `cxx:1544-1601`: this port's edges are frozen by `TopoBuilder::add`, so
    // the `E1.Free() && E2.Free() && myTopoMode` arm never runs; the
    // `CopyReplaceVertices` arm does.
    if nb < 2 {
        // `cxx:1559-1565`: `sbe.CopyReplaceVertices(E2, V, V)` + `Set(tmpE, n2)`.
        let tmp_e = copy_replace_vertices_with(&e2, Some(&combined), Some(&combined));
        wire_set_edge_composed(wire, n2, &tmp_e);
    } else {
        // `cxx:1584-1599`: `CopyReplaceVertices(E2, V, null)` + `Set(tmpE2, n2)`,
        // then the same on E1 unless DONE1 with a non-closed E2.
        let tmp_e2 = copy_replace_vertices_with(&e2, Some(&combined), None);
        wire_set_edge_composed(wire, n2, &tmp_e2);
        let first_is_last = first_vertex(&e2)
            .zip(last_vertex(&e2))
            .map(|(f, l)| is_same(&f.0, &l.0))
            .unwrap_or(false);
        if check != ConnectedCheck::Done1 || first_is_last {
            let tmp_e1 = copy_replace_vertices_with(&e1, None, Some(&combined));
            wire_set_edge_composed(wire, n1, &tmp_e1);
        }
    }
    true
}

/// `ShapeFix_Wire::FixConnected(prec)` (`ShapeFix_Wire.cxx:557-596`), the API
/// driver `Perform` calls at `cxx:353-357` under `myFixConnectedMode`
/// (`FromSTEP.FixShape.FixConnectedMode`, `STEPControl_Controller.cxx:235`,
/// default `-1` -> `NeedFix(-1, ReorderOK)` = ReorderOK). `myClosedMode` is
/// true (`ShapeFix_Wire.cxx:842`), so `aStop = 0` and every junction is
/// visited, from `NbEdges()` down to 1. `Perform` passes the default
/// `prec = -1.`, which `cxx:1485` turns into `MaxTolerance()`.
/// `ShapeFix_Shape::SetMaxTolerance` propagates through
/// `ShapeFix_Solid`/`ShapeFix_Shell`/`ShapeFix_Face` to the wire tool
/// (`ShapeFix_Solid.cxx:745-749`, `ShapeFix_Shell.cxx:1709-1713`,
/// `ShapeFix_Face.cxx:173-177`), so this is `FromSTEP.FixShape.MaxTolerance3d`
/// (`STEPControl_Controller.cxx:204-206`, set from
/// `STEPControl_ActorRead::myMaxTol` at `cxx:1146`).
fn fix_connected_all(wire: &mut Wire, my_precision: f64, max_tolerance: f64) -> bool {
    let nb = wire_edges_nb(wire);
    if nb == 0 {
        return false;
    }
    let mut done = false;
    for a_i in (1..=nb).rev() {
        done |= fix_connected_one(wire, my_precision, max_tolerance, a_i);
    }
    done
}

/// `ShapeFix_Root::MinTolerance()` at read time: `FromSTEP.FixShape.MinTolerance3d`
/// is `"1.e-7"` (`STEPControl_Controller.cxx:207`), `ShapeProcess_OperLibrary.cxx:807`
/// feeds it to `ShapeFix_Shape::SetMinTolerance`, which propagates through
/// `ShapeFix_Solid`/`_Shell`/`_Face` down to the wire tool
/// (`ShapeFix_Shape.cxx:339-342`, `ShapeFix_Solid.cxx:737-740`,
/// `ShapeFix_Shell.cxx:1701-1704`, `ShapeFix_Face.cxx:165-168`).
const SHAPE_FIX_MIN_TOLERANCE: f64 = CONFUSION;

/// Status of `ShapeAnalysis_Wire::CheckSmall(num, precsmall)`
/// (`ShapeAnalysis_Wire.cxx:765-835`): `ok` is that function's own return value,
/// `fail` and `done2` are the `ShapeExtend` bits it leaves in `myStatus`
/// (`DONE1` is `ok && !done2`, `DONE2` is `ok && done2`).
struct SmallCheck {
    ok: bool,
    /// `ShapeExtend_FAIL1` / `FAIL2`. OCCT only ORs these into
    /// `ShapeFix_Wire::myLastFixStatus` (`ShapeFix_Wire.cxx:1420-1423`,
    /// `1434-1443`), a status no caller of this port reads.
    #[allow(dead_code)]
    fail: bool,
    done2: bool,
}

/// `ShapeAnalysis_Wire::CheckSmall(num, precsmall)`
/// (`ShapeAnalysis_Wire.cxx:765-835`). `num == 0` means the last edge (`cxx:774`).
fn check_small(wire: &Wire, face: &Face, precsmall: f64, num: usize) -> SmallCheck {
    let nb = wire_edges_nb(wire);
    // `cxx:768-771`: `!IsLoaded() || NbEdges() <= 1`.
    if nb <= 1 {
        return SmallCheck { ok: false, fail: false, done2: false };
    }
    let edges = edges_of_wire(wire);
    let n = if num > 0 { num } else { nb };
    let Some(e) = edges.get(n - 1) else {
        return SmallCheck { ok: false, fail: false, done2: false };
    };
    // `cxx:777-786`: a degenerated edge that still carries a pcurve on the face
    // is left alone; one without a pcurve is a `FAIL1` and falls through.
    let mut fail = false;
    if BRepTool::is_degenerated(e) {
        if curve_on_surface_oriented(e, face, false).is_some() {
            return SmallCheck { ok: false, fail: false, done2: false };
        }
        fail = true;
    }
    // `cxx:787-794`: missing vertices are a `FAIL2`.
    let (Some(v1), Some(v2)) = (first_vertex(e), last_vertex(e)) else {
        return SmallCheck { ok: false, fail: true, done2: false };
    };
    let p1 = BRepTool::vertex_point(&v1);
    let p2 = BRepTool::vertex_point(&v2);
    // `cxx:796`: `Min(myPrecision, precsmall)` is commented out in OCCT, so the
    // caller's `precsmall` is used as-is.
    let prec = precsmall;
    // `cxx:797-800`: `dist > prec` is not small enough (any `FAIL1` survives).
    if p1.distance(&p2) > prec {
        return SmallCheck { ok: false, fail, done2: false };
    }
    // `cxx:806-826`: the midpoint of the 3d curve, else of the pcurve lifted by
    // the surface; with neither, `FAIL1` and the first vertex is used.
    let mid = match BRepTool::edge_curve(e) {
        Some(c3d) => {
            let (cf, cl) = BRepTool::edge_parameters(e);
            c3d.d0(0.5 * (cf + cl))
        }
        None => {
            let lifted = BRepTool::face_surface(face).and_then(|surf| {
                curve_on_surface_oriented(e, face, false).map(|(c2d, cf, cl)| {
                    let p2m = c2d.d0(0.5 * (cf + cl));
                    surf.d0(p2m.x(), p2m.y())
                })
            });
            match lifted {
                Some(m) => m,
                None => {
                    fail = true;
                    p1.clone()
                }
            }
        }
    };
    // `cxx:827-830`.
    if mid.distance(&p1) > prec || mid.distance(&p2) > prec {
        return SmallCheck { ok: false, fail, done2: false };
    }
    // `cxx:832`: `DONE1` when both vertices are the same shape, else `DONE2`.
    SmallCheck { ok: true, fail, done2: !is_same(&v1.0, &v2.0) }
}

/// `ShapeFix_Wire::FixSmall(num, lockvtx, precsmall)`
/// (`ShapeFix_Wire.cxx:1404-1472`).
fn fix_small_one(
    wire: &mut Wire,
    face: &Face,
    my_precision: f64,
    precsmall: f64,
    lockvtx: bool,
    num: usize,
) -> bool {
    let nb = wire_edges_nb(wire);
    // `cxx:1407-1417`: `!IsLoaded() || NbEdges() <= 1`, then the null analyzer
    // guard (this port always has the face, so the analyzer is never null).
    if nb <= 1 {
        return false;
    }
    let n = if num > 0 { num } else { nb };
    let check = check_small(wire, face, precsmall, n);
    // `cxx:1420-1428`: a `FAIL` is recorded but not fatal, and without a `DONE`
    // bit there is nothing to fix.
    if !check.ok {
        return false;
    }
    // `cxx:1433-1444`: a small edge whose vertices are not the same is only
    // removed when `myTopoMode` is on and the vertices are not locked.
    if check.done2 && lockvtx {
        return false;
    }
    // `cxx:1449-1456`: `Context()->Remove(WireData()->Edge(n))` then
    // `WireData()->Remove(n)`, plus the `FixAdvWire.FixSmall.MSG0` warning.
    // UNPORTED: this port has no `ShapeBuild_ReShape` context at read time, so
    // only the wire removal is applied (same convention as `fix_degenerated`);
    // `SendWarning` needs `Message_Msg`, which is not ported either.
    wire_remove_edge(wire, n);
    // `cxx:1459-1468`: `FixConnected(n <= NbEdges() ? n : 1, precsmall)` when the
    // removed edge had two distinct vertices. A `FAIL` there only sets `FAIL3`
    // in `myLastFixStatus`, which the driver ignores, so the result is dropped.
    if check.done2 {
        let nb2 = wire_edges_nb(wire);
        let cnum = if n <= nb2 { n } else { 1 };
        let _ = fix_connected_one(wire, my_precision, precsmall, cnum);
    }
    true
}

/// `ShapeFix_Wire::FixSmall(lockvtx, precsmall)` (`ShapeFix_Wire.cxx:538-553`),
/// the API driver `Perform` calls at `cxx:334` under
/// `NeedFix(myFixSmallMode, myTopoMode)`. `FromSTEP.FixShape.FixSmallMode` is
/// "-1" (`STEPControl_Controller.cxx:233`), so with `myTopoMode` true on the
/// face path (`ShapeFix_Shape.cxx:200`) the pass runs on every read-in wire,
/// with `lockvtx = !myTopoMode || !ReorderOK` and `precsmall = MinTolerance()`.
fn fix_small_all(
    wire: &mut Wire,
    face: &Face,
    my_precision: f64,
    precsmall: f64,
    lockvtx: bool,
) -> bool {
    // `cxx:546-550`: `for (int i = NbEdges(); i > 0; i--)`. `i` is initialised
    // only once, so after a removal the following indices address the shrunken
    // list exactly as OCCT's do.
    let mut done = false;
    let mut i = wire_edges_nb(wire);
    while i > 0 {
        done |= fix_small_one(wire, face, my_precision, precsmall, lockvtx, i);
        i -= 1;
    }
    done
}

/// `ShapeFix_Wire::FixDegenerated` (`ShapeFix_Wire.cxx:2130-2205`), the pass
/// `ShapeFix_Wire::Perform` runs at `cxx:386-392` under `myFixDegeneratedMode`
/// (`FromSTEP.FixShape.FixDegeneratedMode`, `STEPControl_Controller.cxx:236`,
/// default `-1`). The read-time chain is `ShapeProcess_OperLibrary.cxx:846`
/// (`sfw->FixDegeneratedMode() = ctx->IntegerVal("FixDegeneratedMode", -1)`) ->
/// `ShapeFix_Face::Perform` -> its first wire loop (`ShapeFix_Face.cxx:423`
/// `theAdvFixWire->Perform`, `StatusDegenerated` read at `cxx:431`; the second
/// loop disables the mode again at `cxx:520`) -> `ShapeFix_Wire::Perform`.
/// `NeedFix` maps `-1` to true (`ShapeFix_Root.lxx:101-104` with the default
/// `def = true`, `ShapeFix_Root.hxx:112`), so OCCT does run this pass on a STEP
/// read.
///
/// `Shape-1`'s 16 `SPHERICAL_SURFACE` faces are the read-shape effect of this
/// pass: the file has 116 `EDGE_CURVE`s while the OCCT read shape has 132
/// edges, the 16 extra being one degenerated edge per spherical surface. Each
/// adds one vertex reference forward and the same one reversed (`+301 0 -301 0
/// *` in the OCCT read shape), so the read shape keeps the file's 61
/// `VERTEX_POINT`s with no new vertex -- `FixMissingSeam` would have created one
/// per edge (`ShapeFix_Face.cxx:1984` `B.MakeVertex(V, mySurf->Value(...))`) and
/// given 77. All 16 new edges carry the range `0 1.01893653209489`, which is
/// `vect2d.Magnitude()` for a sphere at `-22.5 -6.129 2.5` r=1.9 (the face's own
/// u-span; all 16 spheres are the same shape, so it is one value), not `2*pi`,
/// again ruling out `FixMissingSeam`'s sphere branch (`ShapeFix_Face.cxx:1920-1925`,
/// `aRange = 2*pi`; that pass is enabled at read time too,
/// `ShapeProcess_OperLibrary.cxx:830`, so it is the range, not the mode, that
/// excludes it). Wire and vertex counts are unchanged (63 wires, 61 vertices in
/// file and reference).
///
/// `ShapeAnalysis_Wire::CheckDegenerated` (`ShapeAnalysis_Wire.cxx:896-1116`).
///
/// `num` is the 1-based wire position of the edge to test (`0` = the last edge,
/// `cxx:903`). The singularity side reuses `SurfaceSingularities`
/// (`pcurve_full/pcurve_ranges.rs`), the port of `ShapeAnalysis_Surface::ComputeSingularities`
/// / `DegeneratedValues` / `IsDegenerated`.
enum DegeneratedCheck {
    /// Returned false without `ShapeExtend_FAIL2`: nothing to fix.
    None,
    /// `ShapeExtend_FAIL2` (`cxx:930-931`, `cxx:945-947`, `cxx:1047-1048`,
    /// `cxx:1093-1097`): `ShapeFix_Wire::FixDegenerated` drops the edge
    /// (`cxx:2148-2152`).
    Remove,
    /// `ShapeExtend_DONE` (`cxx:1110`): insert (`lack` = `ShapeExtend_DONE1`,
    /// `ShapeExtend_WireData::Add`) or replace (`ShapeExtend_DONE2`,
    /// `ShapeExtend_WireData::Set`) a degenerated edge spanning `p2d1`..`p2d2`.
    Found { p2d1: GpPnt2d, p2d2: GpPnt2d, lack: bool },
}

/// `ShapeAnalysis_Edge::HasPCurve(edge, face)` (`ShapeAnalysis_Edge.cxx:183-191`):
/// a `BRep_Tool::CurveOnSurface` representation exists.
fn has_pcurve(edge: &Edge, face: &Face) -> bool {
    curve_on_surface_oriented(edge, face, true).is_some()
}

fn check_degenerated(wire: &Wire, face: &Face, prec: f64, num: usize) -> DegeneratedCheck {
    let Some(surf) = BRepTool::face_surface(face) else {
        return DegeneratedCheck::None;
    };
    let edges = edges_of_wire(wire);
    let nb = edges.len();
    if nb < 1 || num > nb {
        return DegeneratedCheck::None;
    }
    // `cxx:903-905`.
    let n2 = if num > 0 { num } else { nb };
    let n1 = if n2 > 1 { n2 - 1 } else { nb };
    let n3 = if n2 < nb { n2 + 1 } else { 1 };
    let e1 = edges[n1 - 1].clone();
    let e2 = edges[n2 - 1].clone();
    let e3 = edges[n3 - 1].clone();

    // `cxx:913-934`: an edge already flagged degenerated *and* carrying a
    // pcurve is only re-checked for a pcurve that no longer spans the
    // singularity (OCC7630) -- that one is removable.
    if BRepTool::is_degenerated(&e2) && has_pcurve(&e2, face) {
        if has_pcurve(&e1, face) && has_pcurve(&e3, face) {
            let (Some(pc2), Some(pc1), Some(pc3)) = (
                pcurve_at(&e2, face),
                pcurve_at(&e1, face),
                pcurve_at(&e3, face),
            ) else {
                return DegeneratedCheck::None;
            };
            let p21 = pc2.0.d0(pc2.1);
            let p22 = pc2.0.d0(pc2.2);
            let p12 = pc1.4;
            let p31 = pc3.3;
            if (p12.distance(&p31) - p21.distance(&p22)).abs() > 2.0 * PCONFUSION {
                return DegeneratedCheck::Remove;
            }
        }
        return DegeneratedCheck::None;
    }

    // `cxx:938-948`: `n1 != n2` (OCC320, two sequences of degenerated edges on
    // separate surface bounds) and the previous edge is degenerated without a
    // pcurve.
    if n1 != n2 && BRepTool::is_degenerated(&e1) && !has_pcurve(&e1, face) {
        if BRepTool::is_degenerated(&e2) {
            return DegeneratedCheck::Remove;
        }
        return DegeneratedCheck::None;
    }

    // `cxx:950-970`.
    let (Some(vp), Some(v0), Some(v1v), Some(v2v)) = (
        first_vertex(&e1),
        last_vertex(&e1),
        first_vertex(&e2),
        last_vertex(&e2),
    ) else {
        return DegeneratedCheck::None;
    };
    let pp = BRepTool::vertex_point(&vp);
    let p0 = BRepTool::vertex_point(&v0);
    let p1 = BRepTool::vertex_point(&v1v);
    let p2 = BRepTool::vertex_point(&v2v);
    let tol1 = BRepTool::vertex_tolerance(&v1v);
    let prec_first = prec.min(tol1);
    let prec_fin = prec.max(tol1);
    let prec_vtx = if prec < tol1 { 2.0 * prec_fin } else { prec_fin };
    // OCCT's `forward` (`cxx:959`) is unused by
    // `ShapeAnalysis_Surface::DegeneratedValues` (`cxx:373-411`,
    // `const bool /*forward*/`).

    let sing = crate::pcurve_full::SurfaceSingularities::compute(surf.as_ref());
    let mut p2d1 = GpPnt2d::new(0.0, 0.0);
    let mut p2d2 = GpPnt2d::new(0.0, 0.0);
    let mut dgnr = false;
    // `cxx:974-991`: the edge's own wire is already closed on the singular
    // point.
    if p1.distance(&p2) <= prec_first {
        if let Some((a, b, _, _)) = sing.min_gap(&p1, prec_vtx) {
            p2d1 = a;
            p2d2 = b;
            dgnr = true;
            // `cxx:979-990`: do not turn a closed edge whose mid-point is away
            // from the singular point into a degenerated one.
            let (a3, b3) = BRepTool::edge_parameters(&e2);
            if a3.is_finite() && b3.is_finite() {
                if let Some(c3d) = BRepTool::edge_curve(&e2) {
                    let pm = c3d.d0(0.5 * (a3 + b3));
                    if pm.square_distance(&p1) > prec_vtx * prec_vtx {
                        dgnr = false;
                    }
                }
            }
        }
    }
    let mut lack = false;
    if !dgnr {
        // `cxx:995-999`.
        if n1 != n2
            && p1.distance(&pp) <= prec_first
            && sing.is_degenerated(&pp, prec_first)
            && !BRepTool::is_degenerated(&e1)
        {
            return DegeneratedCheck::None;
        }
        // `cxx:1007-1037`: a missing degenerated edge, found by the same
        // minimum-gap singularity search as `DegeneratedValues`.
        if p0.distance(&p1) <= prec_fin {
            if let Some((a, b, _, _)) = sing.min_gap(&p1, prec_vtx) {
                p2d1 = a;
                p2d2 = b;
                lack = true;
            }
        }
    }
    if !lack && !dgnr {
        if BRepTool::is_degenerated(&e2) && !has_pcurve(&e2, face) {
            return DegeneratedCheck::Remove;
        }
        return DegeneratedCheck::None;
    }

    // `cxx:1064-1086`: parametrize the new pcurve exactly from the end of the
    // previous edge's pcurve to the start of the next one.
    if lack || n1 != n2 {
        if let Some((_, _, _, _, pb)) = pcurve_at(&e1, face) {
            p2d1 = pb;
        }
        let target = if dgnr { &e3 } else { &e2 };
        if let Some((_, _, _, a, _)) = pcurve_at(target, face) {
            p2d2 = a;
        }
    }

    // `cxx:1089-1098`: the fix is postponed to `ShapeFix_Wire::FixLacking` when
    // the pcurve is no longer degenerate.
    if !is_degenerated_2d(surf.as_ref(), p2d1, p2d2, prec_vtx, 10.0) {
        if BRepTool::is_degenerated(&e2) {
            return DegeneratedCheck::Remove;
        }
        return DegeneratedCheck::None;
    }

    // `cxx:1101-1107`: the parametric space is already closed. `gp::Resolution()`
    // is `REAL_SMALL`.
    let ads_max = occt_geom::approx_same_parameter::u_resolution(surf.as_ref(), prec)
        .max(occt_geom::approx_same_parameter::v_resolution(surf.as_ref(), prec));
    if p2d1.distance(&p2d2) <= ads_max + occt_core::precision::REAL_SMALL {
        return DegeneratedCheck::None;
    }
    DegeneratedCheck::Found { p2d1, p2d2, lack }
}

/// `ShapeExtend_WireData::Add(edge, atnum)` (`ShapeExtend_WireData.cxx:259-281`):
/// `myEdges->InsertBefore(atnum, edge)` with a 1-based `atnum` (`0` appends).
/// `TShape::children` is the port's `myEdges`, so the new edge is inserted at
/// the same position instead of through `TopoBuilder::add_edge` (which is
/// guarded by `TShape::free()` and would skip a wire already bound to its
/// face, `builder.rs:32-52`).
fn wire_insert_edge_before(wire: &mut Wire, at: usize, edge: &Edge) {
    let mut t = wire.0.tshape.write().expect("poisoned TShape lock");
    let pos = if at == 0 {
        t.children.len()
    } else {
        (at - 1).min(t.children.len())
    };
    t.children.insert(pos, edge.0.clone());
    t.set_modified(true);
}

/// `ShapeExtend_WireData::Set(edge, num)` (`ShapeExtend_WireData.cxx:456-476`):
/// `myEdges->SetValue(num, edge)`, which ignores an out-of-range `num`.
fn wire_set_edge(wire: &mut Wire, at: usize, edge: &Edge) {
    let mut t = wire.0.tshape.write().expect("poisoned TShape lock");
    let len = t.children.len();
    let pos = if at == 0 { len.saturating_sub(1) } else { at - 1 };
    if pos < len {
        t.children[pos] = edge.0.clone();
        t.set_modified(true);
    }
}

/// `ShapeExtend_WireData::Set(edge, num)` for an edge read back through
/// [`edges_of_wire`]. `ShapeExtend_WireData::Init` stores the edges with the
/// wire's own orientation already composed (`ShapeExtend_WireData.cxx:114-121`)
/// and `ShapeExtend_WireData::Wire()` rebuilds a fresh FORWARD wire out of that
/// list (`ShapeExtend_WireData.cxx:651-685`). Our children stay under the
/// original wire, so undo the wire's own orientation before storing, exactly as
/// `fix_reorder_wire` does; `edges_of_wire` then re-composes to the value OCCT
/// would have.
fn wire_set_edge_composed(wire: &mut Wire, at: usize, edge: &Edge) {
    let mut e = edge.clone();
    if wire.0.orientation() == Orientation::Reversed {
        e.0.reverse();
    }
    wire_set_edge(wire, at, &e);
}

/// `ShapeExtend_WireData::Remove(num)` (`ShapeExtend_WireData.cxx:446-452`).
fn wire_remove_edge(wire: &mut Wire, at: usize) {
    let mut t = wire.0.tshape.write().expect("poisoned TShape lock");
    let len = t.children.len();
    let pos = if at == 0 { len.saturating_sub(1) } else { at - 1 };
    if pos < len {
        t.children.remove(pos);
        t.set_modified(true);
    }
}

/// Number of edges `ShapeExtend_WireData` holds (`ShapeExtend_WireData.cxx:576-579`).
fn wire_edges_nb(wire: &Wire) -> usize {
    edges_of_wire(wire).len()
}

/// `ShapeFix_Wire::FixDegenerated(const int num)`
/// (`ShapeFix_Wire.cxx:2130-2205`).
///
/// Returns `(ShapeExtend_DONE, ShapeExtend_DONE2)`; the `DONE2` bit is what the
/// no-argument driver uses to drop duplicated degenerated edges (`cxx:1050`).
fn fix_degenerated(wire: &mut Wire, face: &Face, prec: f64, num: usize) -> (bool, bool) {
    let nb = wire_edges_nb(wire);
    if nb < 1 {
        return (false, false);
    }
    match check_degenerated(wire, face, prec, num) {
        DegeneratedCheck::None => (false, false),
        // `cxx:2148-2152`: FAIL2 -> `WireData()->Remove(num)` + DONE3.
        DegeneratedCheck::Remove => {
            wire_remove_edge(wire, num);
            (true, false)
        }
        DegeneratedCheck::Found { p2d1, p2d2, lack } => {
            let Some(surf) = BRepTool::face_surface(face) else {
                return (false, false);
            };
            // `cxx:2160-2163`.
            let vect2d = GpVec2d::new(p2d2.x() - p2d1.x(), p2d2.y() - p2d1.y());
            let Ok(dir2d) = GpDir2d::from_vec2d(&vect2d) else {
                return (false, false);
            };
            let line2d: Arc<dyn Curve2d> = Arc::new(Geom2dLine::from_pnt_dir(p2d1, dir2d));
            let mag = vect2d.magnitude();

            // `cxx:2165-2168`: `B.MakeEdge(degEdge)`, `B.Degenerated(degEdge,
            // true)`, `B.UpdateEdge(degEdge, line2d, Face(), Precision::Confusion())`,
            // `B.Range(degEdge, Face(), 0., vect2d.Magnitude())`.
            //
            // `BRep_Builder::Degenerated` clears the 3D curve (`cxx:1082-1084`),
            // and this port's `EdgeGeom` cannot hold a null one, so the edge is
            // built over the curve `BRepAdaptor_Curve(edge, face)` produces once
            // the 3D curve is gone: the pcurve image on the surface
            // (`Adaptor3d_CurveOnSurface`, `meshing/edge_discret.rs`).
            let deg_curve: Arc<dyn Curve> = Arc::new(
                crate::meshing::edge_discret::CurveOnSurface::new(
                    line2d.clone(),
                    surf.clone(),
                    0.0,
                    mag,
                ),
            );
            let b = TopoBuilder::new();
            let mut deg = b.make_edge(deg_curve, 0.0, mag);
            let reg = GeometryRegistry::global();
            let face_key = GeometryRegistry::shape_key(&face.0);
            reg.set_degenerated(&deg.0, true);
            reg.set_edge_pcurve(&deg.0, face_key, line2d);
            reg.set_pcurve_range(&deg.0, face_key, 0.0, mag);

            // `cxx:2171-2186`.
            let n2 = if num > 0 { num } else { nb };
            let n1 = if n2 > 1 { n2 - 1 } else { nb };
            let n3 = if lack { n2 } else if n2 < nb { n2 + 1 } else { 1 };
            let edges = edges_of_wire(wire);
            if let Some(v1) = last_vertex(&edges[n1 - 1]) {
                let mut v = v1.0.clone();
                v.set_orientation(Orientation::Forward);
                b.add(&mut deg.0, &v);
            }
            if let Some(v2) = first_vertex(&edges[n3 - 1]) {
                let mut v = v2.0.clone();
                v.set_orientation(Orientation::Reversed);
                b.add(&mut deg.0, &v);
            }
            deg.0.set_orientation(Orientation::Forward);

            // `cxx:2188-2198`: `Add` (DONE1) inserts, `Set` (DONE2) replaces.
            if lack {
                wire_insert_edge_before(wire, n2, &deg);
            } else {
                wire_set_edge(wire, n2, &deg);
            }
            (true, !lack)
        }
    }
}

/// `ShapeFix_Wire::FixDegenerated()` (`ShapeFix_Wire.cxx:1034-1076`), the
/// driver `ShapeFix_Wire::Perform` calls at `cxx:386-392`.
fn fix_degenerated_all(wire: &mut Wire, face: &Face, prec: f64) -> bool {
    let mut done = false;
    let mut last_coded = -1i32;
    let mut prev_coded = 0i32;
    // `myClosedMode` is true (`ShapeFix_Wire.cxx:171`) so `stop = 0`
    // (`cxx:1045`) and the scan runs down to index 1.
    let mut i = wire_edges_nb(wire) as isize;
    while i > 0 {
        let idx = i as usize;
        let (d, coded2) = fix_degenerated(wire, face, prec, idx);
        done |= d;
        let coded = i32::from(coded2);
        if last_coded == -1 {
            last_coded = coded;
        }
        // `cxx:1050-1071`: PRO7226 -- drop a duplicated degenerated edge and
        // clear the flag of the edge that follows it.
        if coded == 1
            && (prev_coded == 1 || (i == 1 && last_coded == 1))
            && wire_edges_nb(wire) > 1
        {
            wire_remove_edge(wire, idx);
            if prev_coded == 0 {
                i = wire_edges_nb(wire) as isize;
            }
            if let Some(e) = edges_of_wire(wire).get(i as usize - 1) {
                GeometryRegistry::global().set_degenerated(&e.0, false);
            }
            // `B.Degenerated(sbwd->Edge(i++), false)` -- the `i++` is undone by
            // the loop's `i--`.
            i += 1;
            prev_coded = 0;
        } else {
            prev_coded = coded;
        }
        i -= 1;
    }
    done
}

/// `ShapeAnalysis_Wire::CheckLacking(num, Tolerance, p2d1, p2d2)`
/// (`ShapeAnalysis_Wire.cxx:1711-1793`).
///
/// `None` covers both `FAIL` (missing vertex / pcurve, `cxx:1732-1750`) and the
/// "no gap" exit `myMax2d < tol2d * tol2d` (`cxx:1776-1780`) - the two states
/// `FixLacking` treats the same way (`cxx:3629-3635`).
struct LackingCheck {
    p2d1: GpPnt2d,
    p2d2: GpPnt2d,
    /// `MaxDistance2d()` (`cxx:1782`).
    max2d: f64,
    /// `MaxDistance3d()` (`cxx:1783`).
    max3d: f64,
    /// `LastCheckStatus(ShapeExtend_DONE2)` (`cxx:1786-1791`).
    done2: bool,
}

fn check_lacking(wire: &Wire, face: &Face, num: usize, tolerance: f64) -> Option<LackingCheck> {
    let nb = wire_edges_nb(wire);
    if nb < 1 {
        return None;
    }
    // `cxx:1723-1726`.
    let n2 = if num > 0 { num } else { nb };
    let n1 = if n2 > 1 { n2 - 1 } else { nb };
    let edges = edges_of_wire(wire);
    let e1 = edges[n1 - 1].clone();
    let e2 = edges[n2 - 1].clone();
    let surf = BRepTool::face_surface(face)?;
    // `cxx:1729-1741`: `LastVertex(E1)` / `FirstVertex(E2)` through
    // `BRepTools::Compare`.
    let (v1, v2) = (last_vertex(&e1)?, first_vertex(&e2)?);
    if !vertices_coincide(&v1, &v2) {
        return None;
    }
    // `cxx:1746-1767`: `sae.PCurve(..., orient=true)` orders the range as
    // FirstVertex -> LastVertex, so `b` is the junction parameter of E1 and `a`
    // the junction parameter of E2 (`ShapeAnalysis_Edge.cxx:218-227`).
    let (c2d1, _a1, b1) = curve_on_surface_oriented(&e1, face, true)?;
    let (p2d1, mut v1t) = c2d1.d1(b1);
    if e1.0.orientation().is_reversed() {
        v1t.reverse();
    }
    let (c2d2, a2, _b2) = curve_on_surface_oriented(&e2, face, true)?;
    let (p2d2, mut v2t) = c2d2.d1(a2);
    if e2.0.orientation().is_reversed() {
        v2t.reverse();
    }
    // `cxx:1768-1780`.
    let v12 = GpVec2d::new(p2d2.x() - p2d1.x(), p2d2.y() - p2d1.y());
    let mut max2d = v12.square_magnitude();
    let mut tol = BRepTool::vertex_tolerance(&v1).max(BRepTool::vertex_tolerance(&v2));
    if tolerance > occt_core::precision::REAL_SMALL && tolerance < tol {
        tol = tolerance; // `cxx:1773`
    }
    let tol2d = 2.0
        * occt_geom::approx_same_parameter::u_resolution(surf.as_ref(), tol)
            .max(occt_geom::approx_same_parameter::v_resolution(surf.as_ref(), tol));
    if max2d < tol2d * tol2d {
        return None;
    }
    max2d = max2d.sqrt();
    let max3d = tol * max2d / tol2d.max(occt_core::precision::REAL_SMALL);
    // `cxx:1786-1791`: `myMax2d < PConfusion`, or the 2d gap points against the
    // wire direction at either end ("back-going" zigzag).
    let done2 = max2d < PCONFUSION
        || (v1t.square_magnitude() > occt_core::precision::REAL_SMALL
            && v12.angle(&v1t).abs() > 0.9 * std::f64::consts::PI)
        || (v2t.square_magnitude() > occt_core::precision::REAL_SMALL
            && v12.angle(&v2t).abs() > 0.9 * std::f64::consts::PI);
    Some(LackingCheck {
        p2d1,
        p2d2,
        max2d,
        max3d,
        done2,
    })
}

/// `TryNewPCurve` (`ShapeFix_Wire.cxx:2215-2246`): deviation of `c2d` against
/// the edge's own 3D curve, measured by `ShapeFix_Edge::FixSameParameter`
/// (`cxx:2241`) on a temporary edge.
fn try_new_pcurve(
    edge: &Edge,
    face: &Face,
    c2d: Arc<dyn Curve2d>,
    first: f64,
    last: f64,
) -> Option<(Arc<dyn Curve2d>, f64, f64, f64)> {
    let curve = BRepTool::edge_curve(edge)?; // `cxx:2223-2227`
    let (f, l) = BRepTool::edge_parameters(edge);
    // `cxx:2230-2236`: `BRepBuilderAPI_MakeEdge(crv, f, l)` + `SetRange3d`.
    let builder = TopoBuilder::new();
    let tmp = builder.make_edge(curve, f, l);
    let reg = GeometryRegistry::global();
    let face_key = GeometryRegistry::shape_key(&face.0);
    // `cxx:2238-2243`: `BRepBuilderAPI_MakeEdge` leaves the edge tolerance at
    // `BRep_TEdge`'s constructor value `RealEpsilon()`
    // (`BRep_TEdge.cxx:32-38`), and `B.UpdateEdge(edge, c2d, face, 0.)`
    // (`cxx:2239`) is `TE->UpdateTolerance(0.)`, a max, so it stays there.
    // `B.Range(edge, face, first, last)` (`cxx:2240`), `B.SameRange(edge,
    // false)` (`cxx:2242`) and no `SameParameter` flag (`cxx:2243`).
    reg.set_edge_pcurve(&tmp.0, face_key, c2d);
    reg.set_pcurve_range(&tmp.0, face_key, first, last);
    reg.set_edge_tolerance(&tmp.0, occt_core::precision::COMPUTATIONAL);
    reg.set_same_range(&tmp.0, false);
    // `cxx:2244`: `sfe->FixSameParameter(edge, face)`.
    fix_same_parameter(&tmp, face);
    // `cxx:2245-2246`.
    let (new_c2d, nf, nl) = curve_on_surface_oriented(&tmp, face, false)?;
    Some((new_c2d, nf, nl, BRepTool::edge_tolerance(&tmp)))
}

/// The four in-place outputs of `TryBendingPCurve`
/// (`ShapeFix_Wire.cxx:3534-3538`): `c2d`, `first`, `last`, `tol`. The C++
/// function assigns only the ones its control flow reaches, so a failure after
/// `c2d = bs` (`cxx:3593`) leaves `c2d` assigned at the new candidate while
/// `first`/`last` keep the range `sae.PCurve` wrote (`cxx:3542`) and `tol` keeps
/// its previous value. `None` models a null handle.
struct BendParam {
    curve: Option<Arc<dyn Curve2d>>,
    first: f64,
    last: f64,
    tol: f64,
}

impl BendParam {
    /// `double bendtol1 = 0., bendtol2 = 0.; ... bendf1 = 0., bendl1 = 0. ...`
    /// (`cxx:3664-3666`).
    fn new() -> Self {
        BendParam {
            curve: None,
            first: 0.0,
            last: 0.0,
            tol: 0.0,
        }
    }
}

/// `TryBendingPCurve(E, face, p2d, end, c2d, first, last, tol)`
/// (`ShapeFix_Wire.cxx:3532-3613`). `out` is mutated in place, exactly as the
/// C++ reference parameters are.
fn try_bending_pcurve(
    edge: &Edge,
    face: &Face,
    p2d: GpPnt2d,
    end: bool,
    out: &mut BendParam,
) -> bool {
    // `cxx:3542-3545`: `sae.PCurve(E, face, c2d, first, last, false)`.
    let Some((c2d, first, last)) = curve_on_surface_oriented(edge, face, false) else {
        return false;
    };
    out.curve = Some(c2d.clone());
    out.first = first;
    out.last = last;
    // `cxx:3550-3586`: `c2d->IsKind(Geom2d_BSplineCurve)` (with a seam, the
    // second pcurve) takes the copy arm; every other kind goes through
    // `Geom2dTrimmedCurve` + `Geom2dConvert::CurveToBSplineCurve`
    // (`cxx:3556-3559`). This port does not expose a knot vector, so the
    // `SetPole` edit cannot be reproduced on an imported B-spline pcurve; only
    // the Line arm (`Geom2dConvert.cxx:211-225`, mirrored by
    // `line_trim_to_reparam_bspline` above) is ported. Any other kind reports
    // "no bend candidate" - the same outcome OCCT gives when its own
    // `CurveToBSplineCurve` throws (`cxx:3605-3612`).
    if !is_geom2d_line_curve(c2d.as_ref()) {
        return false;
    }
    // `cxx:3555-3559`: the trimmed-line B-spline. `line_trim_to_reparam_bspline`
    // builds the same curve for `GeomLib::SameRange`; here the knots stay the
    // trimmed curve's own, because `Segment` is never reached (below).
    let tc = Geom2dTrimmedCurve::new(c2d.clone(), first, last);
    let (k0, k1) = (tc.first_parameter(), tc.last_parameter());
    let (u0, u1) = (k0, k1);
    if (u1 - u0).abs() <= PCONFUSION {
        return false;
    }
    let (p_at_u0, p_at_u1) = if first <= last {
        (tc.d0(first), tc.d0(last))
    } else {
        (tc.d0(last), tc.d0(first))
    };
    // `cxx:3567-3591`: `par = (end ? last : first)`; the pole test succeeds
    // because the flat knot vector has multiplicity 2 > degree 1 at both ends
    // (`Multiplicity(1) > Degree()` at `cxx:3569`, `Multiplicity(NbKnots) >
    // Degree()` at `cxx:3573`), so `SetPole` always applies at `par` and
    // `Segment` at `cxx:3579` is unreachable for this arm.
    let par = if end { last } else { first };
    let (np0, np1) = if (par - u0).abs() <= (par - u1).abs() {
        (p2d, p_at_u1) // `SetPole(1, p2d)`
    } else {
        (p_at_u0, p2d) // `SetPole(NbPoles, p2d)`
    };
    let Ok(bs) = Geom2dBSplineCurve::new(
        vec![np0.x(), np1.x()],
        vec![np0.y(), np1.y()],
        vec![u0, u0, u1, u1],
        1,
    ) else {
        return false;
    };
    let candidate: Arc<dyn Curve2d> = Arc::new(bs);
    // `cxx:3593`: `c2d = bs` happens before the tolerance check.
    out.curve = Some(candidate.clone());
    // `cxx:3596-3601`.
    match try_new_pcurve(edge, face, candidate, first, last) {
        Some((c, f, l, tol)) => {
            out.curve = Some(c);
            out.first = f;
            out.last = l;
            out.tol = tol;
            true
        }
        None => false,
    }
}

/// `ShapeFix_Wire::FixLacking(const int num, const bool force)`
/// (`ShapeFix_Wire.cxx:3617-3975`).
///
/// NOT PORTED: the `Context()->Replace` calls at `cxx:3881-3899` only update
/// the `ShapeBuild_ReShape` map (this pass runs without a context), and the two
/// `BRep_Tool::Degenerated` arms at `cxx:3819-3825` have empty bodies in OCCT
/// itself.
fn fix_lacking_one(
    wire: &mut Wire,
    face: &Face,
    num: usize,
    force: bool,
    preci: f64,
    max_tol: f64,
) -> bool {
    // `ShapeAnalysis_Wire::IsReady()`: `!myWire.IsNull() && NbEdges() > 0`
    // (`ShapeAnalysis_Wire.hxx`), checked at `cxx:3620-3623`.
    let nb = wire_edges_nb(wire);
    if nb == 0 {
        return false;
    }
    // `cxx:3627-3635`: `CheckLacking(num, force ? Precision() : 0.)`, then bail
    // unless it reported `DONE`.
    let tolerance = if force { preci } else { 0.0 };
    let Some(check) = check_lacking(wire, face, num, tolerance) else {
        return false;
    };
    // `cxx:3640-3647`.
    let n2 = if num > 0 { num } else { nb };
    let n1 = if n2 > 1 { n2 - 1 } else { nb };
    let edges = edges_of_wire(wire);
    let e1 = edges[n1 - 1].clone();
    let e2 = edges[n2 - 1].clone();
    let (Some(v1), Some(v2)) = (last_vertex(&e1), first_vertex(&e2)) else {
        return false;
    };
    // `cxx:3652`.
    let tol = BRepTool::vertex_tolerance(&v1).max(BRepTool::vertex_tolerance(&v2));
    let dist2d = check.max2d;
    let inctol = check.max3d;
    let p2d1 = check.p2d1;
    let p2d2 = check.p2d2;
    let Some(surf) = BRepTool::face_surface(face) else {
        return false;
    };
    let reg = GeometryRegistry::global();
    let face_key = GeometryRegistry::shape_key(&face.0);
    // `BRep_Tool::IsClosed(E, face)` is "the edge carries two pcurves on this
    // face" (`BRep_Tool.cxx:372-381`).
    let is_closed = |e: &Edge| reg.edge_pcurves(&e.0, face_key).len() >= 2;
    let tol_e1 = BRepTool::edge_tolerance(&e1);
    let tol_e2 = BRepTool::edge_tolerance(&e2);
    let mut tol1 = CONFUSION; // `cxx:3659-3660`
    let mut tol2 = CONFUSION;
    let mut p3d1: Option<GpPnt> = None;
    let mut p3d2: Option<GpPnt> = None;

    // `cxx:3662-3713`: bend speculation. `myGeomMode` is true
    // (`ShapeFix_Wire.cxx:169`).
    let mut bp1 = BendParam::new();
    let mut bp2 = BendParam::new();
    if !is_closed(&e1) && !is_closed(&e2) {
        let mid = GpPnt2d::new(0.5 * (p2d1.x() + p2d2.x()), 0.5 * (p2d1.y() + p2d2.y()));
        let end1 = e1.0.orientation() == Orientation::Forward;
        let end2 = e2.0.orientation() == Orientation::Reversed;
        let mut ok1 = try_bending_pcurve(&e1, face, mid, end1, &mut bp1);
        let mut ok2 = try_bending_pcurve(&e2, face, mid, end2, &mut bp2);
        if ok1 && !ok2 {
            // `cxx:3686-3697`.
            bp2.tol = tol_e2;
            ok1 = try_bending_pcurve(&e1, face, p2d2, end1, &mut bp1);
        } else if !ok1 && ok2 {
            // `cxx:3698-3709`. OCCT passes `E2.Orientation() == TopAbs_FORWARD`
            // for this retry while the first call at `cxx:3681` tested
            // REVERSED; reproduced as written.
            bp1.tol = tol_e1;
            let end2_retry = e2.0.orientation() == Orientation::Forward;
            ok2 = try_bending_pcurve(&e2, face, p2d1, end2_retry, &mut bp2);
        }
        if !ok1 && !ok2 {
            bp1.curve = None; // `cxx:3710-3712`: `bendc1.Nullify()`
        }
    }

    // `cxx:3716-3817`: selector of solutions.
    let mut do_increase = false;
    let mut do_add_long = false;
    let mut do_add_closed = false;
    let mut do_add_degen = false;
    let mut do_bend = false;
    let bendtol1 = bp1.tol;
    let bendtol2 = bp2.tol;
    let bend_ready = bp1.curve.is_some() && bp2.curve.is_some();
    if bend_ready
        && ((bendtol1 < tol_e1 && bendtol2 < tol_e2)
            || (inctol < preci && bendtol1 < inctol && bendtol2 < inctol))
    {
        do_bend = true; // `cxx:3727-3732`
    } else if inctol < preci {
        do_increase = true; // `cxx:3735-3738`
    } else if !reg.is_degenerated_edge(&e2.0) && !reg.is_degenerated_edge(&e1.0) {
        // `cxx:3744-3745`: `myTopoMode` is forced true for the face path by
        // `ShapeFix_Shape.cxx:200`, so this block always runs.
        {
            // `cxx:3748-3762`: `sae.Curve3d(E1, c3d, a, b, true)` gives `b` at
            // LastVertex and `a` at FirstVertex.
            let Some(c1) = BRepTool::edge_curve(&e1) else {
                return false; // `cxx:3749-3753` FAIL1
            };
            let (mut a1, mut b1p) = BRepTool::edge_parameters(&e1);
            if e1.0.orientation().is_reversed() {
                std::mem::swap(&mut a1, &mut b1p);
            }
            let q1 = c1.d0(b1p);
            let dist2d3d1 = q1.distance(&surf.d0(p2d1.x(), p2d1.y()));
            let Some(c2) = BRepTool::edge_curve(&e2) else {
                return false; // `cxx:3756-3760` FAIL1
            };
            let (mut a2, mut b2p) = BRepTool::edge_parameters(&e2);
            if e2.0.orientation().is_reversed() {
                std::mem::swap(&mut a2, &mut b2p);
            }
            let q2 = c2.d0(a2);
            let dist2d3d2 = q2.distance(&surf.d0(p2d2.x(), p2d2.y()));
            // `cxx:3764-3769`.
            tol1 = tol_e1.max(dist2d3d1);
            tol2 = tol_e2.max(dist2d3d2);
            let tol0 = tol1 + tol2;
            let dist3d2 = q1.square_distance(&q2);
            p3d1 = Some(q1);
            p3d2 = Some(q2);
            if !check.done2
                && dist3d2 > 1.25 * tol0 * tol0
                && (force || dist3d2 > preci * preci || inctol > max_tol)
            {
                do_add_long = true; // `cxx:3771-3777`
            }
        }
        // `cxx:3781-3791`.
        if !do_add_long
            && inctol < max_tol
            && !is_degenerated_2d(surf.as_ref(), p2d1, p2d2, 2.0 * tol, 10.0)
        {
            if bend_ready && bendtol1 < inctol && bendtol2 < inctol {
                do_bend = true;
            } else {
                do_increase = true;
            }
        } else if !do_add_long {
            // `cxx:3799-3815`.
            let p1 = BRepTool::vertex_point(&v1);
            let p2 = BRepTool::vertex_point(&v2);
            let pv = GpPnt::new(
                0.5 * (p1.x() + p2.x()),
                0.5 * (p1.y() + p2.y()),
                0.5 * (p1.z() + p2.z()),
            );
            let pm = surf.d0(0.5 * (p2d1.x() + p2d2.x()), 0.5 * (p2d1.y() + p2d2.y()));
            let dist = pv.distance(&pm);
            if dist <= tol {
                do_add_degen = true;
            } else {
                // `cxx:3806`: `myTopoMode`, which is true here. The
                // `dist <= MaxTolerance()` fallback at `cxx:3810-3814` belongs
                // to the `!myTopoMode` arm and is unreachable.
                do_add_closed = true;
            }
        }
    }

    // `cxx:3830-3921`: add the new edge. `DONE2` is set at `cxx:3920` on all
    // three arms (`DONE3` at `cxx:3910` for the degenerated edge, `DONE4` at
    // `cxx:3916` for the closed one), so every arm reports a fix.
    let mut done = false;
    if do_add_long || do_add_degen || do_add_closed {
        let builder = TopoBuilder::new();
        // `cxx:3835-3849`.
        let (new_v1, new_v2) = if do_add_long {
            let v1n = builder.make_vertex(p3d1.expect("doAddLong sets p3d1"), 0.0);
            let v2n = builder.make_vertex(p3d2.expect("doAddLong sets p3d2"), 0.0);
            v1n.set_tolerance(1.001 * tol1); // `cxx:3846`
            v2n.set_tolerance(1.001 * tol2); // `cxx:3847`
            (v1n, v2n)
        } else {
            (v1.clone(), v2.clone())
        };
        // `cxx:3852-3862`: `B.MakeEdge`, `B.Degenerated(edge, true)` before the
        // curve, `B.UpdateEdge(edge, theLine2d, face, Precision::Confusion())`,
        // `B.Range(edge, face, 0, dist2d)`.
        let v12 = GpVec2d::new(p2d2.x() - p2d1.x(), p2d2.y() - p2d1.y());
        let Ok(dir2d) = GpDir2d::from_vec2d(&v12) else {
            return false;
        };
        let line2d: Arc<dyn Curve2d> = Arc::new(Geom2dLine::from_pnt_dir(p2d1, dir2d));
        // `ShapeBuild_Edge::BuildCurve3d` (`cxx:3866`) approximates the 3D image
        // of the pcurve on the face. This port's `EdgeGeom` always needs a 3D
        // curve, and the pcurve image is the curve `BRepAdaptor_Curve` resolves
        // once the same edge exists - the adapter `fix_degenerated` uses for its
        // own curve-less edge (`ShapeFix_Wire.cxx:2165-2168`).
        let curve3d: Arc<dyn Curve> = Arc::new(crate::meshing::edge_discret::CurveOnSurface::new(
            line2d.clone(),
            surf.clone(),
            0.0,
            dist2d,
        ));
        let mut new_edge = builder.make_edge(curve3d, 0.0, dist2d);
        reg.set_edge_pcurve(&new_edge.0, face_key, line2d);
        reg.set_pcurve_range(&new_edge.0, face_key, 0.0, dist2d);
        reg.set_edge_tolerance(&new_edge.0, CONFUSION);
        if do_add_degen {
            reg.set_degenerated(&new_edge.0, true);
        }
        builder.add_edge_vertices(&mut new_edge, &new_v1, &new_v2);
        // `cxx:3873-3901`: `doAddLong` re-points the two adjacent edges at the
        // new vertices.
        if do_add_long {
            let first_arg = if n1 == n2 { Some(&new_v2) } else { None };
            let ne1 = copy_replace_vertices_with(&e1, first_arg, Some(&new_v1));
            wire_set_edge_composed(wire, n1, &ne1);
            if n1 != n2 {
                let ne2 = copy_replace_vertices_with(&e2, Some(&new_v2), None);
                wire_set_edge_composed(wire, n2, &ne2);
            }
        }
        // `cxx:3905-3920`: `DONE2` is encoded on every arm of this branch.
        done = true;
        // `ShapeExtend_WireData::Add(edge, n2)` (`cxx:3919`) inserts before
        // `n2`; the stored children keep the wire's own orientation composed
        // (`ShapeExtend_WireData.cxx:114-121`), so undo it here as
        // `wire_set_edge_composed` does for `Set`.
        let mut ins = new_edge.clone();
        if wire.0.orientation() == Orientation::Reversed {
            ins.0.reverse();
        }
        wire_insert_edge_before(wire, n2, &ins);
    } else if inctol > tol && inctol < max_tol {
        // `cxx:3924-3933`.
        if bend_ready && bendtol1 < inctol && bendtol2 < inctol {
            do_bend = true;
        } else {
            do_increase = true;
        }
    }

    // `cxx:3937-3956`: bend the pcurves. `B.UpdateEdge(E, bendc, face, bendtol)`
    // writes the pcurve, the COS range and the edge tolerance (`TE->
    // UpdateTolerance(Tol)` is a max, `BRep_Builder.cxx:655-671` /
    // `BRep_TEdge.hxx`).
    if do_bend {
        if let Some(bc1) = bp1.curve.as_ref() {
            reg.set_edge_pcurve(&e1.0, face_key, bc1.clone());
            reg.set_pcurve_range(&e1.0, face_key, bp1.first, bp1.last);
            reg.set_edge_tolerance(&e1.0, tol_e1.max(bendtol1));
        }
        if let Some(bc2) = bp2.curve.as_ref() {
            reg.set_edge_pcurve(&e2.0, face_key, bc2.clone());
            reg.set_pcurve_range(&e2.0, face_key, bp2.first, bp2.last);
            reg.set_edge_tolerance(&e2.0, tol_e2.max(bendtol2));
        }
        // `cxx:3943-3946`: `B.UpdateVertex` on all four vertices.
        for (v, t) in [
            (first_vertex(&e1), bendtol1),
            (last_vertex(&e1), bendtol1),
            (first_vertex(&e2), bendtol2),
            (last_vertex(&e2), bendtol2),
        ] {
            if let Some(v) = v {
                v.set_tolerance(t);
            }
        }
        // `cxx:3948-3951`: `FixSelfIntersectingEdge(n1)`, `FixSelfIntersectingEdge(n2)`,
        // `FixIntersectingEdges(n2)` on the bent edges. UNPORTED: both helpers
        // are gated by the 2D self-intersection checks, whose blocker chain is
        // spelled out at the FixSelfIntersection marker in
        // `check_pcurves_and_shift` below.
        done = true;
    }

    // `cxx:3959-3963`.
    if do_increase {
        v1.set_tolerance(1.001 * inctol);
        v2.set_tolerance(1.001 * inctol);
        done = true;
    }

    done
}

/// `ShapeFix_Wire::FixLacking()` (`ShapeFix_Wire.cxx:1284-1295`) - `Perform`
/// calls it at `cxx:448` under `NeedFix(myFixLackingMode, ReorderOK)` with
/// `FromSTEP.FixShape.FixLackingMode = -1` (`STEPControl_Controller.cxx:237`),
/// so the gate is `ReorderOK`. Wired in `check_pcurves_and_shift` below.
///
/// The pass is translated branch for branch, except the `doBend` tail's
/// `FixSelfIntersectingEdge` / `FixIntersectingEdges` calls
/// (`ShapeFix_Wire.cxx:3948-3951`, marked UNPORTED in `fix_lacking_one`).
/// Previously it could not be wired: the in-memory cylinder used by
/// `step::tests::cylinder_roundtrip` had a lateral wire whose two rings ran the
/// same way (sum(du) = +4pi), so the closure junction showed
/// `p2d1 = (2*pi, 0)` / `p2d2 = (-2*pi, 0)` - both mapping to the same 3D
/// vertex `(2, 0, 0)`. `ShapeAnalysis_Wire::CheckLacking`
/// (`ShapeAnalysis_Wire.cxx:1767-1790`) has no period normalization, so it
/// reported DONE1 with `myMax2d = 4*pi` / `myMax3d = 2*pi` and `FixLacking`
/// selected the degenerated-edge branch (`cxx:3782-3806`). That malformed input
/// wire is fixed in `primitives.rs` (`BRepPrimCylinder::make_cylinder` now
/// traverses the two rings oppositely, as OCCT's swept lateral face does), so
/// every junction of the cylinder's lateral loop is gap-free and this pass is a
/// no-op there.

fn fix_lacking_all(wire: &mut Wire, face: &Face, force: bool, preci: f64, max_tol: f64) -> bool {
    let mut done = false;
    // `cxx:1291-1295`: `myClosedMode` is true (`ShapeFix_Wire.cxx:171`), so
    // `start = 1`; `NbEdges()` is re-read every iteration, so an inserted edge
    // extends the loop.
    let mut i = 1usize;
    while i <= wire_edges_nb(wire) {
        done |= fix_lacking_one(wire, face, i, force, preci, max_tol);
        i += 1;
    }
    done
}

// ---------------------------------------------------------------------------
// ShapeFix_Wire::FixNotchedEdges / ShapeAnalysis_Wire::CheckNotchedEdges
// ---------------------------------------------------------------------------

/// `Geom_Curve::IsClosed()`: `Geom_Circle::IsClosed()` (`Geom_Circle.cxx:73-76`)
/// and `Geom_Ellipse` (base `Geom_Curve::IsClosed()` is `IsPeriodic()`) are the
/// periodic flag; `Geom_TrimmedCurve::IsClosed()` (`Geom_TrimmedCurve.cxx:155-168`)
/// is `IsPeriodic()` (whole periods of a periodic basis, see
/// `occt-geom/trimmed.rs:159-172`) or coincident ends; `Geom_BSplineCurve`
/// (`Geom_BSplineCurve_1.cxx:146-149`) and `Geom_OffsetCurve`
/// (`Geom_OffsetCurve.cxx:436-442`) compare the ends against
/// `Precision::Computational()`, standing in for `Geom_BezierCurve`'s
/// construction flag too.
///
/// NOTE this is NOT `GeomAdaptor_Curve::IsClosed()`
/// (`GeomAdaptor_Curve.cxx:576-584`), which only compares the *window* ends and
/// is what `shape_analysis_curve::project_act` needs (`cxx:340`).
fn geom_curve_is_closed(c: &dyn Curve) -> bool {
    if c.is_periodic() {
        return true;
    }
    let (f, l) = (c.first_parameter(), c.last_parameter());
    f.is_finite()
        && l.is_finite()
        && c.d0(f).square_distance(&c.d0(l)) <= occt_core::precision::COMPUTATIONAL
}

/// `ShapeAnalysis_Edge::IsClosed3d(edge)` (`ShapeAnalysis_Edge.cxx:129-142`):
/// the edge carries a 3D curve that is closed, and its first and last vertices
/// are the same `TShape` (`FirstVertex` / `LastVertex` with the default
/// `CumOri = false`, i.e. the stored FORWARD / REVERSED child).
fn edge_is_closed_3d(edge: &Edge) -> bool {
    let Some(c3d) = BRepTool::edge_curve(edge) else {
        return false;
    };
    if !geom_curve_is_closed(c3d.as_ref()) {
        return false;
    }
    match (edge_vertices(edge).0, edge_vertices(edge).1) {
        (Some(v1), Some(v2)) => is_same(&v1.0, &v2.0),
        _ => false,
    }
}

/// `ShapeAnalysis_Wire::ProjectInside` (`ShapeAnalysis_Wire.cxx:1837-1860`).
fn project_inside(
    ad: &crate::meshing::edge_discret::CurveOnSurface,
    pnt: &GpPnt,
    preci: f64,
    adjust_to_ends: bool,
) -> Projection {
    let proj = project_adaptor(ad, pnt, preci, adjust_to_ends);
    let u_first = ad.first_parameter();
    let u_last = ad.last_parameter();
    if proj.param < u_first {
        let p = ad.d0(u_first);
        return Projection { distance: p.distance(pnt), point: p, param: u_first };
    }
    if proj.param > u_last {
        let p = ad.d0(u_last);
        return Projection { distance: p.distance(pnt), point: p, param: u_last };
    }
    proj
}

/// `ShapeAnalysis_Wire::CheckNotchedEdges(num, shortNum, param, Tolerance)`
/// (`ShapeAnalysis_Wire.cxx:1864-2000`). `Some` carries the `shortNum` and
/// `param` out-parameters; the `ShapeExtend_FAIL*` bits it leaves in `myStatus`
/// (`cxx:1890`, `:1896`, `:1906`, `:1919`) are not modelled, as no caller of
/// this port reads them.
struct NotchedCheck {
    short_num: usize,
    param: f64,
}

fn check_notched_edges(
    wire: &Wire,
    face: &Face,
    num: usize,
    tolerance: f64,
) -> Option<NotchedCheck> {
    // `cxx:1869-1878`: `IsReady()` and the `num == 0` meaning "the last edge".
    let edges = edges_of_wire(wire);
    let nb = edges.len();
    if nb == 0 {
        return None;
    }
    let n2 = if num > 0 { num } else { nb };
    let n1 = if n2 > 1 { n2 - 1 } else { nb };
    let e1 = edges.get(n1.checked_sub(1)?)?.clone();
    let e2 = edges.get(n2.checked_sub(1)?)?.clone();
    if BRepTool::is_degenerated(&e1) || BRepTool::is_degenerated(&e2) {
        return None; // `cxx:1880-1883`
    }

    // `cxx:1885-1898`: `sae.LastVertex(E1)` / `sae.FirstVertex(E2)` with
    // `CumOri = false`, then `BRepTools::Compare(V1, V2)`.
    let (Some(v1), Some(v2)) = (edge_vertices(&e1).1, edge_vertices(&e2).0) else {
        return None;
    };
    if !vertices_coincide(&v1, &v2) {
        return None;
    }

    // `cxx:1900-1932`: the end tangents of both pcurves.
    let (c2d1, a1, b1) = curve_on_surface_oriented(&e1, face, false)?;
    let (c2d2, a2, b2) = curve_on_surface_oriented(&e2, face, false)?;
    let (p2d1, tan1) = if e1.0.orientation().is_reversed() {
        c2d1.d1(a1)
    } else {
        let (p, d) = c2d1.d1(b1);
        (p, d.reversed())
    };
    let (p2d2, tan2) = if e2.0.orientation().is_reversed() {
        let (p, d) = c2d2.d1(b2);
        (p, d.reversed())
    } else {
        c2d2.d1(a2)
    };
    // `cxx:1934-1937`: `gp::Resolution()` is `RealSmall()`.
    if tan2.magnitude() < REAL_SMALL || tan1.magnitude() < REAL_SMALL {
        return None;
    }
    // `cxx:1939-1942`.
    if tan2.angle(&tan1).abs() > 0.1 || p2d1.distance(&p2d2) > tolerance {
        return None;
    }

    // `cxx:1944-1952`: both pcurves are read as 3D curves lying in the XY plane
    // (`Geom_Plane(gp_Pln())`), which is what `Adaptor3d_CurveOnSurface` over a
    // plane resolves.
    let plane: Arc<dyn Surface> = Arc::new(GeomPlane::new(GpPln::new(GpAx3::standard())));
    let ad1 = crate::meshing::edge_discret::CurveOnSurface::new(c2d1.clone(), plane.clone(), a1, b1);
    let ad2 = crate::meshing::edge_discret::CurveOnSurface::new(c2d2.clone(), plane.clone(), a2, b2);

    // `cxx:1954-1968`.
    let pt2 = if e2.0.orientation().is_forward() { c2d2.d0(b2) } else { c2d2.d0(a2) };
    let pt1 = if e1.0.orientation().is_forward() { c2d1.d0(a1) } else { c2d1.d0(b1) };
    let start2 = GpPnt::new(pt2.x(), pt2.y(), 0.0);
    let start1 = GpPnt::new(pt1.x(), pt1.y(), 0.0);
    let proj1 = project_inside(&ad1, &start2, tolerance, false);
    let proj2 = project_inside(&ad1, &start1, tolerance, false);
    if proj1.distance > tolerance && proj2.distance > tolerance {
        return None;
    }

    // `cxx:1970-1986`: the longer edge is the one the other is projected onto.
    let (long_ad, short_ad, len_p, first_p, short_num, param) = if proj1.distance < proj2.distance
    {
        (&ad1, &ad2, b2 - a2, a2, n2, proj1.param)
    } else {
        (&ad2, &ad1, b1 - a1, a1, n1, proj2.param)
    };

    // `cxx:1988-1998`: 22 interior samples of the short edge must stay on the
    // long one. The last `Project` default is `AdjustToEnds = true` (`cxx:208`).
    let step = len_p / 23.0;
    let mut sample = first_p;
    for _ in 1..23 {
        let p = short_ad.d0(sample);
        if project_adaptor(long_ad, &p, tolerance, true).distance > tolerance {
            return None;
        }
        sample += step;
    }
    Some(NotchedCheck { short_num, param })
}

/// `ShapeFix_Wire::CopyReversePcurves` (`ShapeFix_Wire.cxx:4108-4176`),
/// file-static in the cxx: append to `to` every `CurveOnSurface` representation
/// of `from` that `to` does not already carry for the same surface *and*
/// location (`cxx:4125-4143`), reversing the first pcurve and its range when
/// `reverse` (`cxx:4152-4160`). The representation's location participates in
/// the identity test, and locations are identity throughout this port, so the
/// face key alone decides it.
fn copy_reverse_pcurves(to: &Edge, from: &Edge, reverse: bool) {
    let reg = GeometryRegistry::global();
    let Some(from_geom) = reg.edge_geom(&from.0) else {
        return;
    };
    let Some(mut to_geom) = reg.edge_geom(&to.0) else {
        return;
    };
    for (face_key, from_pcurves) in &from_geom.pcurves {
        if to_geom.pcurves.contains_key(face_key) {
            continue;
        }
        let (fp0, lp0) = from_geom
            .pcurve_ranges
            .get(face_key)
            .copied()
            .unwrap_or((from_geom.first, from_geom.last));
        let mut range = (fp0, lp0);
        let mut copied: Vec<Arc<dyn Curve2d>> = Vec::with_capacity(from_pcurves.len());
        for (i, pc) in from_pcurves.iter().enumerate() {
            if i == 0 && reverse {
                // `cxx:4152-4160`: `ReversedParameter` on the old range, then
                // `Reverse`, then the swap.
                let fp = super::pcurve_ranges::reversed_parameter(pc.as_ref(), fp0);
                let lp = super::pcurve_ranges::reversed_parameter(pc.as_ref(), lp0);
                let mut c = pc.clone_dyn();
                c.reverse();
                copied.push(Arc::from(c));
                range = (lp, fp);
            } else {
                // `cxx:4168-4171`: `PCurve2` is copied as-is.
                copied.push(Arc::from(pc.clone_dyn()));
            }
        }
        to_geom.pcurves.insert(*face_key, copied);
        to_geom.pcurve_ranges.insert(*face_key, range);
    }
    reg.set_edge(&to.0, to_geom);
}

/// `ShapeFix_Wire::FixDummySeam(num)` (`ShapeFix_Wire.cxx:4213-4289`).
///
/// The `Context()` branches (`cxx:4241-4256`, `:4264-4272`) are the port's
/// in-place wire mutations instead of a `ShapeBuild_ReShape` record, exactly as
/// in the other wire passes; `toRemove` is hard-coded false at `cxx:4227`, so
/// the `Context()->Remove` arm is dead.
fn fix_dummy_seam(wire: &mut Wire, num: usize) {
    let nb = wire_edges_nb(wire);
    if nb < 2 || num == 0 || num > nb {
        return;
    }
    let num1 = if num == nb { 1 } else { num + 1 };
    let edges = edges_of_wire(wire);
    let (Some(e1), Some(e2)) = (edges.get(num - 1).cloned(), edges.get(num1 - 1).cloned()) else {
        return;
    };
    // `cxx:4219-4221`: `V1 = sae.FirstVertex(E1)`, `V2 = sae.LastVertex(E2)`,
    // `Vm = sbv.CombineVertex(V1, V2, 1.0001)` (`CumOri = false`).
    let (Some(v1), Some(v2)) = (edge_vertices(&e1).0, edge_vertices(&e2).1) else {
        return;
    };
    let vm = combine_vertex(&v1, &v2, 1.0001);

    // `cxx:4230-4233`: `Vs = sae.FirstVertex(E2)`, replaced by `Vm` when it is
    // already one of the two merged vertices.
    let mut vs = edge_vertices(&e2).0;
    if let Some(vs_ref) = vs.as_ref() {
        if is_same(&vs_ref.0, &v1.0) || is_same(&vs_ref.0, &v2.0) {
            vs = Some(vm.clone());
        }
    }
    let new_edge = copy_replace_vertices_with(&e2, vs.as_ref(), Some(&vm));
    copy_reverse_pcurves(&new_edge, &e1, e1.0.orientation() == e2.0.orientation());
    let reg = GeometryRegistry::global();
    reg.set_same_range(&new_edge.0, false); // `cxx:4237-4238`
    reg.set_same_parameter(&new_edge.0, false);

    // `cxx:4258-4272`: the neighbouring edges are re-pointed at `Vm`.
    let next = if num1 == nb { 1 } else { num1 + 1 };
    let prev = if num > 1 { num - 1 } else { nb };
    if let Some(prev_e) = edges.get(prev - 1).cloned() {
        let tmp = copy_replace_vertices_with(&prev_e, None, Some(&vm));
        wire_set_edge_composed(wire, prev, &tmp);
    }
    if let Some(next_e) = edges.get(next - 1).cloned() {
        let tmp = copy_replace_vertices_with(&next_e, Some(&vm), None);
        wire_set_edge_composed(wire, next, &tmp);
    }

    // `cxx:4275-4288`: both notch edges leave the wire, the higher index first.
    let (n1, n2) = if num < num1 { (num, num1) } else { (num1, num) };
    wire_remove_edge(wire, n2);
    wire_remove_edge(wire, n1);
}

/// `ShapeFix_Wire::FixNotchedEdges()` (`ShapeFix_Wire.cxx:3977-4310`).
///
/// `Perform` calls it at `cxx:400` under
/// `myFixTailMode <= 0 && NeedFix(myFixNotchedEdgesMode, ReorderOK)`.
/// `FromSTEP.FixShape.FixTailMode` is "0" (`STEPControl_Controller.cxx:249`) so
/// the first half holds, and `FixNotchedEdgesMode` is "-1"
/// (`STEPControl_Controller.cxx:248`), so `NeedFix(-1, ReorderOK)`
/// (`ShapeFix_Root.lxx:101-104`) is ReorderOK.
///
/// `WireData()`/`UpdateWire()` (`cxx:3989-3991`, `:4083-4086`) are the port's
/// in-memory wire; `myStatusNotches` (`cxx:4102`) is a status no caller of this
/// port reads.
///
/// Wiring this pass moved no count: a temporary probe on `check_notched_edges`'s
/// `Some` return (since removed) fired 0 times on a full `export_data_obj` run,
/// i.e. `ShapeAnalysis_Wire::CheckNotchedEdges` finds no notched pair in this
/// corpus, exactly as it does in OCCT.
fn fix_notched_edges(wire: &mut Wire, face: &Face, min_tol: f64, max_tol: f64) -> bool {
    let mut done = false;
    let mut i = 1usize;
    // `cxx:3994`: `i <= NbEdges() && NbEdges() > 2`, both re-read every pass.
    while i <= wire_edges_nb(wire) && wire_edges_nb(wire) > 2 {
        let Some(check) = check_notched_edges(wire, face, i, min_tol) else {
            i += 1;
            continue;
        };
        // `cxx:3999-4007`.
        let nb = wire_edges_nb(wire);
        let n2 = if i > 0 { i } else { nb };
        let n1 = if n2 > 1 { n2 - 1 } else { nb };
        let is_remove_first = n1 == check.short_num;
        let to_split = if n2 == check.short_num { n1 } else { n2 };
        let edges = edges_of_wire(wire);
        let Some(split_e) = edges.get(to_split - 1).cloned() else {
            i += 1;
            continue;
        };
        let Some((c2d, a, b)) = curve_on_surface_oriented(&split_e, face, true) else {
            i += 1;
            continue;
        };
        let orient = split_e.0.orientation();

        // `cxx:4011-4018`: the split point falls on an end of the edge; for a
        // closed edge it may fall on the other end (issue #0029780).
        let on_end = (check.param - if is_remove_first { b } else { a }).abs() <= PCONFUSION
            || (edge_is_closed_3d(&split_e)
                && (check.param - if is_remove_first { a } else { b }).abs() <= PCONFUSION);
        if on_end {
            fix_dummy_seam(wire, n1);
            // `cxx:4019-4021`: the seam edge left the list; the cxx's `i--`
            // cancels the loop's `i++`, so `i` stays put.
        } else {
            // `cxx:4027-4029`: the split point is already the edge's own end.
            if ((if is_remove_first { a } else { b }) - check.param).abs() < PCONFUSION {
                i += 1;
                continue;
            }
            // `cxx:4031-4035`.
            let mut transfer = TransferParametersProj::new();
            transfer.set_max_tolerance(max_tol);
            transfer.init(&split_e, face);
            let (first, last) = if a < b { (a, b) } else { (b, a) };

            // `cxx:4046-4050`: `B.MakeVertex(Vnew,
            // Analyzer()->Surface()->Value(c2d->Value(param)),
            // Precision::Confusion())`.
            let Some(surface) = BRepTool::face_surface(face) else {
                i += 1;
                continue;
            };
            let uv = c2d.d0(check.param);
            let vnew = TopoBuilder::new().make_vertex(surface.d0(uv.x(), uv.y()), CONFUSION);

            // `cxx:4051-4058`: first half, from the edge start to the notch.
            let mut we = split_e.clone();
            we.0.set_orientation(Orientation::Forward);
            let (Some(fv), Some(lv)) = (first_vertex(&we), last_vertex(&we)) else {
                i += 1;
                continue;
            };
            let reg = GeometryRegistry::global();
            let mut new_e1 = copy_replace_vertices_with(&we, Some(&fv), Some(&vnew));
            copy_pcurves(&new_e1, &we);
            transfer.transfer_range(&mut new_e1, first, check.param, true);
            reg.set_same_range(&new_e1.0, false);
            reg.set_same_parameter(&new_e1.0, false);
            // `cxx:4059-4067`: second half, from the notch to the edge end.
            let mut new_e2 = copy_replace_vertices_with(&we, Some(&vnew), Some(&lv));
            copy_pcurves(&new_e2, &we);
            transfer.transfer_range(&mut new_e2, check.param, last, true);
            reg.set_same_range(&new_e2.0, false);
            reg.set_same_parameter(&new_e2.0, false);

            // `cxx:4078-4084`: both halves take the split edge's orientation,
            // and a REVERSED edge swaps them.
            new_e1.0.set_orientation(orient);
            new_e2.0.set_orientation(orient);
            if orient == Orientation::Reversed {
                std::mem::swap(&mut new_e1, &mut new_e2);
            }
            // `cxx:4087-4088`: `Set(newE1, toSplit)` then `Add(newE2, ...)`
            // (0 = append).
            let nb_now = wire_edges_nb(wire);
            wire_set_edge_composed(wire, to_split, &new_e1);
            let at = if to_split == nb_now { 0 } else { to_split + 1 };
            let mut ins = new_e2.clone();
            if wire.0.orientation() == Orientation::Reversed {
                ins.0.reverse();
            }
            wire_insert_edge_before(wire, at, &ins);
            // `cxx:4089`: `FixDummySeam(isRemoveLast ? NbEdges() : toRemove)`.
            let is_remove_last = n1 == nb_now && n2 == 1;
            let target = if is_remove_last { wire_edges_nb(wire) } else { check.short_num };
            fix_dummy_seam(wire, target);
        }
        // `cxx:4099`: `DONE1` on both arms. The cxx's `i--` after each fix
        // cancels the loop's `i++`, so `i` is left unchanged here; the `Set` /
        // `Add` above shift the list under it.
        done = true;
    }
    done
}

/// STEP post-pass after pcurve association.
/// Source: `TranslateEdgeLoop::CheckPCurves` first check, then
/// `ShapeFix_Wire::Perform` FixReorder (`cxx:317-325` / `cxx:487`) then
/// `FromSTEP.FixShape.FixShiftedMode` (default on) / `cxx:939` FixShifted,
/// then SameRange / CheckPCurveRange / FixAddPCurve / FixSameParameter
/// (`ShapeFix_Wire.cxx:938-995`), then `ShapeFix_Wire::FixDegenerated`
/// (`cxx:386-392`, `fix_degenerated_all` below).
///
/// `TranslateEdgeLoop.cxx:105-112`: for a plane face `CheckPCurves` calls
/// `RemovePCurves(aWire, aFace)` and returns, so the pcurve-range healing does
/// not run on a plane's edges and their tolerances keep the imported values
/// (linkrods FACE 35 edge 4 keeps 1.2432146551310009e-07). That removal is
/// ported in `check_pcurve_rep_range` below (`cxx:110-114`), together with the
/// non-planar `XSAlgo_ShapeProcessor::CheckPCurve` call at `cxx:175`. The
/// `ShapeFix_Wire` passes that follow are independent of the face type.
pub fn check_pcurves_and_shift(wire: &mut Wire, face: &Face, preci: f64) {
    check_pcurve_rep_range(wire, face, preci);
    // `ShapeFix_Wire::Perform` (`cxx:317-325`): FixReorder before FixEdgeCurves.
    let mut reorder_ok = fix_reorder_wire(wire, face);
    // `ShapeFix_Wire::Perform` (`cxx:333-347`): FixSmall under
    // `NeedFix(myFixSmallMode, myTopoMode)`. `FromSTEP.FixShape.FixSmallMode` is
    // "-1" (`STEPControl_Controller.cxx:233`) and `myTopoMode` is forced true by
    // `ShapeFix_Shape.cxx:200` for the face path, so the pass runs, with
    // `lockvtx = !myTopoMode || !ReorderOK` and `precsmall = MinTolerance()`.
    if fix_small_all(wire, face, preci, SHAPE_FIX_MIN_TOLERANCE, !reorder_ok) {
        // `cxx:339-343`: after a removal the wire may reorder better; the retry
        // is gated by `NeedFix(myFixReorderMode = -1, !ReorderOK)` evaluated with
        // the ReorderOK from before FixSmall.
        if !reorder_ok {
            reorder_ok = fix_reorder_wire(wire, face);
        }
    }
    // `ShapeFix_Wire::Perform` (`cxx:353-357`): FixConnected after FixSmall.
    if reorder_ok {
        let _ = fix_connected_all(wire, preci, SHAPE_FIX_MAX_TOLERANCE);
    }
    let edges = edges_of_wire(wire);
    // `ShapeFix_Wire::FixEdgeCurves` (`cxx:652-670`): FixAddPCurve runs for
    // every edge of the wire. `FromSTEP.exec.op = FixShape`
    // (`STEPControl_Controller.cxx:201`) routes every imported face through
    // `ShapeFix_Shape` -> `ShapeFix_Face` -> `ShapeFix_Wire::Perform`, so a
    // face whose STEP edges carry a 3D curve only (no `SURFACE_CURVE` /
    // `PCURVE`, as in ATU01038) is projected onto the surface here. With no
    // pcurve stored the face keeps the natural surface bounds
    // (`BRepTools::AddUVBounds`, `BRepTools.cxx:139-153`), which are infinite
    // for a cylinder/cone. The `cxx:673-863` sub-branch (edge over a
    // singularity) needs a `ShapeBuild_ReShape` context and stays unported.
    // `ShapeFix_Edge::myProjector` lives on the `ShapeFix_Edge` that
    // `ShapeFix_Wire` keeps for the whole wire (`ShapeFix_Wire.cxx:280-300`), so
    // one `myCache` covers every edge of this face's wire.
    let mut cache = crate::pcurve_full::ProjectorCache::new();
    for (i, e) in edges.iter().enumerate() {
        let is_seam = is_seam_use(&edges, i);
        let _ = fix_add_pcurve(e, face, is_seam, preci, &mut cache);
    }
    // `cxx:368-374` turns FixShifted off when the reorder reported FAIL.
    if reorder_ok {
        let _ = fix_shifted_wire(wire, face);
    }
    let reg = GeometryRegistry::global();
    for (i, e) in edges.iter().enumerate() {
        if let Some((c2d, fp2d, lp2d)) = curve_on_surface_oriented(e, face, false) {
            let (first, last) = BRepTool::edge_parameters(e);
            if first.is_finite() && last.is_finite() {
                let (cfp, clp) = (c2d.first_parameter(), c2d.last_parameter());
                if (first - fp2d).abs() > PCONFUSION || (last - lp2d).abs() > PCONFUSION {
                    reg.set_same_range(&e.0, false);
                } else if !check_pcurve_range(first, last, c2d.as_ref()) {
                    // `GeomLib.cxx:862`: same parametrisation length uses SameRange.
                    // Different length stays on CheckPCurveRange -> RemovePCurve +
                    // FixAddPCurve (`ShapeFix_Wire.cxx:973-983`).
                    if cfp.is_finite()
                        && clp.is_finite()
                        && !c2d.is_periodic()
                        && ((clp - cfp) - (last - first)).abs() <= PCONFUSION
                    {
                        let remapped = geom_lib_same_range(c2d, cfp, clp, first, last);
                        replace_pcurve(e, face, remapped);
                        reg.set_same_range(&e.0, true);
                    } else {
                        let is_seam = is_seam_use(&edges, i);
                        reg.remove_pcurves_on_surface(&e.0, &face.0);
                        let _ = fix_add_pcurve(e, face, is_seam, preci, &mut cache);
                    }
                }
            }
        }
        fix_same_parameter(e, face);
    }
    // `ShapeFix_Wire::Perform` (`cxx:386-392`): FixDegenerated runs right after
    // FixEdgeCurves (the loop above) and before FixNotchedEdges / FixTails /
    // FixSelfIntersection / FixLacking / FixVertexTolerance.
    fix_degenerated_all(wire, face, preci);
    // `ShapeFix_Wire::Perform` (`cxx:400`): FixNotchedEdges under
    // `myFixTailMode <= 0 && NeedFix(myFixNotchedEdgesMode, ReorderOK)`.
    // `FromSTEP.FixShape.FixTailMode` is "0" (`STEPControl_Controller.cxx:249`)
    // and `FixNotchedEdgesMode` is "-1" (`STEPControl_Controller.cxx:248`), so
    // the gate is ReorderOK. `cxx:402-406`: a fix is followed by FixShifted.
    if reorder_ok && fix_notched_edges(wire, face, SHAPE_FIX_MIN_TOLERANCE, SHAPE_FIX_MAX_TOLERANCE)
    {
        let _ = fix_shifted_wire(wire, face);
    }
    // `ShapeFix_Wire::Perform` (`cxx:412`): FixTails under
    // `myFixTailMode != 0`. `FromSTEP.FixShape.FixTailMode` is "0"
    // (`STEPControl_Controller.cxx:249`), so OCCT skips this pass in the STEP
    // read-in and porting it here would deviate. Intentionally absent.
    // `ShapeFix_Wire::Perform` (`cxx:428`): FixSelfIntersection under
    // `NeedFix(myFixSelfIntersectionMode, myClosedMode)`. `FixSelfIntersectionMode`
    // is "-1" (`STEPControl_Controller.cxx:238`), so `NeedFix(-1, myClosedMode)`
    // (`ShapeFix_Root.lxx:101-104`) is `myClosedMode` (`ShapeFix_Wire.cxx:171`)
    // and the pass is requested. UNPORTED - not merely missing code, blocked by
    // a whole unported package, so it cannot be reproduced faithfully:
    //
    //   `ShapeFix_Wire::FixSelfIntersection` (`ShapeFix_Wire.cxx:1083-1280`)
    //     first half `FixSelfIntersectingEdge` (`cxx:2708-2914`)
    //       check `ShapeAnalysis_Wire::CheckSelfIntersectingEdge`
    //         (`ShapeAnalysis_Wire.cxx:1269-1341`) -> `Geom2dInt_GInter`
    //       fix   `RemoveLoop` (`cxx:2283-2526`, `:2527-2706`), which also needs
    //             `TryNewPCurve` (`cxx:2215-2260`), `howMuchPCurves`
    //             (`cxx:2262-2281`) and `GeomConvert_CompCurveToBSplineCurve`
    //     second half `FixIntersectingEdges` (`cxx:2966-3269`)
    //       check `ShapeAnalysis_Wire::CheckIntersectingEdges`
    //         (`ShapeAnalysis_Wire.cxx:1360-1545`) -> `Geom2dInt_GInter`
    //       fix   `ComputeLocalDeviation` (`cxx:2917-2964`) +
    //             `ShapeFix_SplitTool::CutEdge` (`ShapeFix_SplitTool.cxx`,
    //             16416 bytes, unported) + `FixSameParameter` (ported)
    //     tail `NeedFix(myFixNonAdjacentIntersectingEdgesMode)` (`cxx:1200-1203`)
    //       is also true - `FixNonAdjacentIntersectingEdgesMode` is "-1"
    //       (`STEPControl_Controller.cxx:255`) - so
    //       `ShapeFix_IntersectionTool::FixSelfIntersectWire`
    //       (`ShapeFix_IntersectionTool.cxx`, 93541 bytes, unported) belongs here
    //     both cores intersect pcurves with `Geom2dInt_GInter`, which is
    //     `IntCurve_IntCurveCurveGen.gxx` instantiated for
    //     `Geom2dInt_Geom2dCurveTool` (`Geom2dInt_GInter_0.cxx:88-89`): conics
    //     go through `Geom2dInt_TheIntConicConicOfGInter` /
    //     `IntCurve_IntConicConic.cxx` (106021 + 37693 bytes), everything else
    //     through `Geom2dInt_TheIntPCurvePCurveOfGInter_0.cxx` +
    //     `Geom2dInt_ExactIntersectionPointOfTheIntPCurvePCurveOfGInter_0.cxx`
    //     + `Intf_InterferencePolygon2d.cxx` (23756 bytes) for its
    //     `IntRes2d_IntersectionPoint` / `IntRes2d_IntersectionSegment` /
    //     `IntRes2d_Transition` results (`IntRes2d_*.cxx`).
    //
    // None of that family exists here. The closest 2D intersector in the port,
    // `occt_geom2d::geom2d_api::intersect_curves` (`geom2d_api.rs:74`, a
    // sampler in `curve_ops::curve2d_intersections`), is a heuristic the bop
    // pave machinery uses; it is not `Geom2dInt_GInter` and returns no
    // transitions or segments, so substituting it would fabricate results for
    // the `IntRes2d_Middle` tests at `ShapeAnalysis_Wire.cxx:1324-1329` and
    // `:1490-1493` instead of reproducing them.
    //
    // The same blocker stops the `doBend` tail of `FixLacking` above
    // (`ShapeFix_Wire.cxx:3948-3951`), whose three calls are exactly
    // `FixSelfIntersectingEdge` / `FixIntersectingEdges`. That gap is live on
    // this corpus: a temporary probe in the `doBend` branch above (since
    // removed) fired 14 times on a full `export_data_obj` run - ATU01038 once,
    // Shape-1 four times, Shape-2 nine times - i.e. exactly on the three models
    // furthest from their oracles.
    //
    // `ShapeFix_Face::Perform` reaches this pass in its second wire round
    // (`ShapeFix_Face.cxx:509-559`), which switches FixSmall / FixConnected /
    // FixEdgeCurves / FixDegenerated off but leaves self-intersection, notches
    // and lacking on; this port's single pass covers both rounds, which is why
    // the marker sits in this pipeline.
    // `ShapeFix_Wire::Perform` (`cxx:448`): FixLacking under
    // `NeedFix(myFixLackingMode, ReorderOK)`; `NeedFix(flag, need)` is
    // `flag < 0 ? need : flag > 0` (`ShapeFix_Root.lxx:101-104`), and
    // `FromSTEP.FixShape.FixLackingMode` is "-1"
    // (`STEPControl_Controller.cxx:237`), so OCCT runs the pass exactly when
    // ReorderOK is true.
    if reorder_ok {
        let _ = fix_lacking_all(wire, face, false, preci, SHAPE_FIX_MAX_TOLERANCE);
    }
}
