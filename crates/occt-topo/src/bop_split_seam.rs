//! `BOPTools_AlgoTools3D::DoSplitSEAMOnFace` both overloads.
//!
//! Source: `BOPTools_AlgoTools3D.cxx` (first overload at 58, origin/split
//! overload at 236). BuildSplitFaces calls these when a source edge is closed
//! on a periodic face and a split of that edge is not yet a seam (`_2.cxx:435`).
//!
//! AABB / sampling stand-in for `Geom2dAPI_ProjectPointOnCurve`. Dual p-curves
//! are stored with [`GeometryRegistry::set_edge_pcurves`].

use std::sync::Arc;

use occt_core::gp::{GpDir2d, GpPnt2d, GpTrsf2d, GpVec2d};
use occt_core::precision::{CONFUSION, PCONFUSION, RESOLUTION};
use occt_geom::Surface;
use occt_geom2d::Curve2d;

use crate::abs::Orientation;
use crate::boptools_2d::{curve_on_surface, intermediate_point};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape};
use crate::tgeometry::GeometryRegistry;

fn shape_key(s: &TopoShape) -> usize {
    GeometryRegistry::shape_key(s)
}

/// `BRep_Tool::IsClosed(edge, face)` (`BRep_Tool.cxx:795-841`): a plane is
/// never closed; otherwise the edge carries two p-curves on the face. OCCT's
/// `BOPTools_AlgoTools3D::DoSplitSEAMOnFace` uses exactly this predicate
/// (`BOPTools_AlgoTools3D.cxx:240,245`) — there is no "edge appears twice"
/// arm (that is `BRepTools::IsReallyClosed`).
pub fn is_closed_on_face(edge: &Edge, face: &Face) -> bool {
    BRepTool::is_closed_edge_face(edge, face)
}

/// `(isUIso, isVIso)` from the p-curve tangent (`IsEdgeIsoline`).
pub fn isoline_uv(edge: &Edge, face: &Face) -> (bool, bool) {
    let (a, b) = GeometryRegistry::global().edge_parameters(&edge.0);
    if !a.is_finite() || !b.is_finite() || b - a < 1e-15 {
        return (false, false);
    }
    let Some(pc) = curve_on_surface(edge, face) else {
        return (false, false);
    };
    let tm = 0.5 * (a + b);
    let (_, t) = pc.d1(tm);
    let sq = t.square_magnitude();
    if sq <= RESOLUTION {
        return (false, false);
    }
    let n = t.divided(sq.sqrt());
    let dp_v = n.crossed(&GpVec2d::new(0.0, 1.0)).abs();
    let dp_u = n.crossed(&GpVec2d::new(1.0, 0.0)).abs();
    const ANG: f64 = 1e-12_f64.max(occt_core::precision::ANGULAR);
    (dp_v <= ANG, dp_u <= ANG)
}

fn surface_period_uv(s: &dyn Surface) -> (f64, f64) {
    let (umin, umax) = s.u_range();
    let (vmin, vmax) = s.v_range();
    let u = if s.is_u_periodic() {
        umax - umin
    } else {
        0.0
    };
    let v = if s.is_v_periodic() {
        vmax - vmin
    } else {
        0.0
    };
    (u, v)
}

fn surface_closed_period(s: &dyn Surface, tol: f64) -> (bool, bool, f64, f64) {
    let (umin, umax) = s.u_range();
    let (vmin, vmax) = s.v_range();
    let mut u_per = s.is_u_periodic();
    let mut v_per = s.is_v_periodic();
    let mut an_u = if u_per { umax - umin } else { 0.0 };
    let mut an_v = if v_per { vmax - vmin } else { 0.0 };
    if !u_per && !v_per {
        let (pu, pv) = surface_period_uv(s);
        if pu > 0.0 {
            u_per = true;
            an_u = pu;
        }
        if pv > 0.0 {
            v_per = true;
            an_v = pv;
        }
    }
    let _ = tol;
    (u_per, v_per, an_u, an_v)
}

fn u_v_resolution(s: &dyn Surface, tol: f64) -> (f64, f64) {
    let (umin, umax) = s.u_range();
    let (vmin, vmax) = s.v_range();
    let du_span = (umax - umin).abs().max(CONFUSION);
    let dv_span = (vmax - vmin).abs().max(CONFUSION);
    let p0 = s.d0(0.5 * (umin + umax), 0.5 * (vmin + vmax));
    let pu = s.d0(0.5 * (umin + umax) + 0.01 * du_span, 0.5 * (vmin + vmax));
    let pv = s.d0(0.5 * (umin + umax), 0.5 * (vmin + vmax) + 0.01 * dv_span);
    let lu = p0.distance(&pu).max(CONFUSION) / (0.01 * du_span);
    let lv = p0.distance(&pv).max(CONFUSION) / (0.01 * dv_span);
    ((tol / lu).max(PCONFUSION), (tol / lv).max(PCONFUSION))
}

fn project_2d_on_curve(p: &GpPnt2d, c: &dyn Curve2d, t1: f64, t2: f64) -> Option<(f64, f64)> {
    let n = 48usize;
    let mut best_d = f64::MAX;
    let mut best_t = 0.5 * (t1 + t2);
    for i in 0..=n {
        let u = i as f64 / n as f64;
        let t = t1 + u * (t2 - t1);
        let q = c.d0(t);
        let d = (q.x() - p.x()).hypot(q.y() - p.y());
        if d < best_d {
            best_d = d;
            best_t = t;
        }
    }
    let span = (t2 - t1).abs().max(PCONFUSION);
    let mut t = best_t;
    for _ in 0..8 {
        let h = 1e-4 * span;
        let q = c.d0(t);
        let (_, d1) = c.d1(t);
        let rx = q.x() - p.x();
        let ry = q.y() - p.y();
        let g = rx * d1.x() + ry * d1.y();
        let mag = d1.square_magnitude().max(RESOLUTION);
        t -= g / mag;
        t = t.clamp(t1.min(t2), t1.max(t2));
        let q2 = c.d0(t);
        let d2 = (q2.x() - p.x()).hypot(q2.y() - p.y());
        if d2 < best_d {
            best_d = d2;
            best_t = t;
        }
        let _ = h;
    }
    Some((best_t, best_d))
}

fn translate_curve(c: &dyn Curve2d, t1: f64, t2: f64, du: f64, dv: f64) -> Arc<dyn Curve2d> {
    let mut tr = GpTrsf2d::identity();
    tr.set_translation_vec(&GpVec2d::new(du, dv));
    let moved = c.transformed(&tr);
    let _ = (t1, t2);
    Arc::from(moved)
}

fn clone_curve(c: &dyn Curve2d) -> Arc<dyn Curve2d> {
    Arc::from(c.clone_dyn())
}

/// First overload (`AlgoTools3D.cxx:58`): shift the split's p-curve by one
/// period when the mid-point sits on a periodic bound, then store both.
pub fn do_split_seam_on_face(a_split: &Edge, a_f: &Face) -> bool {
    let mut a_sp = a_split.clone();
    a_sp.0.set_orientation(Orientation::Forward);
    let a_tol = BRepTool::edge_tolerance(&a_sp);
    let Some(a_s) = BRepTool::face_surface(a_f) else {
        return false;
    };
    let (a_umin, a_umax) = a_s.u_range();
    let (a_vmin, a_vmax) = a_s.v_range();
    let (b_u, b_v, an_u_period, an_v_period) =
        surface_closed_period(a_s.as_ref(), a_tol);
    if !b_u && !b_v {
        return false;
    }
    let Some(c2d1) = curve_on_surface(&a_sp, a_f) else {
        return false;
    };
    let (a, b) = GeometryRegistry::global().edge_parameters(&a_sp.0);
    let a_t = intermediate_point(a, b);
    let (a_p2d, a_vec2d) = c2d1.d1(a_t);
    let a_dir2d1 = if a_vec2d.square_magnitude() <= RESOLUTION {
        GpDir2d::default()
    } else {
        GpDir2d::from_vec2d(&a_vec2d).unwrap_or_default()
    };
    let a_dox = GpDir2d::default();
    let a_doy = GpDir2d::new(0.0, 1.0).unwrap_or_else(|_| GpDir2d { x: 0.0, y: 1.0 });
    let an_u = a_p2d.x();
    let an_v = a_p2d.y();
    let mut an_u1 = an_u;
    let mut an_v1 = an_v;
    let (d_u, d_v) = u_v_resolution(a_s.as_ref(), a_tol);
    let mut b_is_left = false;
    if an_u_period > 0.0 {
        if (an_u - a_umin).abs() < d_u {
            b_is_left = true;
            an_u1 = an_u + an_u_period;
        } else if (an_u - a_umax).abs() < d_u {
            b_is_left = false;
            an_u1 = an_u - an_u_period;
        }
    }
    if an_v_period > 0.0 {
        if (an_v - a_vmin).abs() < d_v {
            b_is_left = true;
            an_v1 = an_v + an_v_period;
        } else if (an_v - a_vmax).abs() < d_v {
            b_is_left = false;
            an_v1 = an_v - an_v_period;
        }
    }
    if (an_u1 - an_u).abs() <= PCONFUSION && (an_v1 - an_v).abs() <= PCONFUSION {
        return false;
    }
    let a_sc_pr = if (an_u1 - an_u).abs() <= PCONFUSION {
        a_dir2d1.dot(&a_doy)
    } else {
        a_dir2d1.dot(&a_dox)
    };
    let a_c1 = clone_curve(c2d1.as_ref());
    let a_c2 = translate_curve(c2d1.as_ref(), a, b, an_u1 - an_u, an_v1 - an_v);
    let fk = shape_key(&a_f.0);
    let (first, second) = if !b_is_left {
        if a_sc_pr < 0.0 {
            (a_c2, a_c1)
        } else {
            (a_c1, a_c2)
        }
    } else if a_sc_pr < 0.0 {
        (a_c1, a_c2)
    } else {
        (a_c2, a_c1)
    };
    GeometryRegistry::global().set_edge_pcurves(&a_sp.0, fk, vec![first, second]);
    let _ = (b_u, b_v, an_u, an_v);
    true
}

/// Origin/split overload (`AlgoTools3D.cxx:236`).
pub fn do_split_seam_on_face_origin(the_e_origin: &Edge, the_e_split: &Edge, the_face: &Face) -> bool {
    if !is_closed_on_face(the_e_origin, the_face) {
        return false;
    }
    if is_closed_on_face(the_e_split, the_face) {
        return true;
    }
    let mut a_e_split = the_e_split.clone();
    a_e_split.0.set_orientation(Orientation::Forward);
    let mut a_face = the_face.clone();
    a_face.0.set_orientation(Orientation::Forward);
    let Some(a_c2d_split) = curve_on_surface(&a_e_split, &a_face) else {
        return false;
    };
    let (a_ts1, a_ts2) = GeometryRegistry::global().edge_parameters(&a_e_split.0);
    let mut e_fwd = the_e_origin.clone();
    e_fwd.0.set_orientation(Orientation::Forward);
    let mut e_rev = the_e_origin.clone();
    e_rev.0.set_orientation(Orientation::Reversed);
    let Some(a_c2d1) = curve_on_surface(&e_fwd, &a_face) else {
        return false;
    };
    let Some(a_c2d2) = curve_on_surface(&e_rev, &a_face) else {
        return false;
    };
    let (a_t1, a_t2) = GeometryRegistry::global().edge_parameters(&the_e_origin.0);
    let a_t = intermediate_point(a_ts1, a_ts2);
    let (a_p_mid, a_v_tgt) = a_c2d_split.d1(a_t);
    let proj1 = project_2d_on_curve(&a_p_mid, a_c2d1.as_ref(), a_t1, a_t2);
    let proj2 = project_2d_on_curve(&a_p_mid, a_c2d2.as_ref(), a_t1, a_t2);
    if proj1.is_none() && proj2.is_none() {
        return false;
    }
    let a_dist1 = proj1.map(|(_, d)| d).unwrap_or(f64::MAX);
    let a_dist2 = proj2.map(|(_, d)| d).unwrap_or(f64::MAX);
    if a_dist1 > PCONFUSION && a_dist2 > PCONFUSION {
        return false;
    }
    let a_new_pnt = if a_dist1 < a_dist2 {
        let t = proj1.map(|(t, _)| t).unwrap_or(a_t);
        a_c2d2.d0(t)
    } else {
        let t = proj2.map(|(t, _)| t).unwrap_or(a_t);
        a_c2d1.d0(t)
    };
    let du = a_new_pnt.x() - a_p_mid.x();
    let dv = a_new_pnt.y() - a_p_mid.y();
    let a_c1 = clone_curve(a_c2d_split.as_ref());
    let a_c2 = translate_curve(a_c2d_split.as_ref(), a_ts1, a_ts2, du, dv);
    let (a_p_proj, a_v_tgt_origin) = if a_dist1 < a_dist2 {
        let t = proj1.map(|(t, _)| t).unwrap_or(a_t);
        a_c2d1.d1(t)
    } else {
        let t = proj2.map(|(t, _)| t).unwrap_or(a_t);
        a_c2d2.d1(t)
    };
    let a_dot = a_v_tgt.x() * a_v_tgt_origin.x() + a_v_tgt.y() * a_v_tgt_origin.y();
    let fk = shape_key(&a_face.0);
    if (a_dist1 < a_dist2) == (a_dot > 0.0) {
        GeometryRegistry::global().set_edge_pcurves(&a_e_split.0, fk, vec![a_c1, a_c2]);
    } else {
        GeometryRegistry::global().set_edge_pcurves(&a_e_split.0, fk, vec![a_c2, a_c1]);
    }
    let _ = a_p_proj;
    true
}

/// Try both overloads as BuildSplitFaces does (`_2.cxx:435-439`).
pub fn make_closed_edge_on_face(origin: &Edge, split: &Edge, face: &Face) -> bool {
    if is_closed_on_face(split, face) {
        return true;
    }
    if do_split_seam_on_face(split, face) {
        return true;
    }
    do_split_seam_on_face_origin(origin, split, face)
}
