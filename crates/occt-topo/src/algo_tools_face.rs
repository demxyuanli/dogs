//! `BOPTools_AlgoTools` face classification: `IsInternalFace`, `GetFaceOff`,
//! `AreFacesSameDomain`, `Sense`, `IsSplitToReverse`, `GetEdgeOnFace`.
//!
//! Source: `BOPTools_AlgoTools.cxx` (`IsInternalFace` at 807 / 895 / 939,
//! `GetFaceOff` at 994, `GetEdgeOff` at 1099, `AreFacesSameDomain` at 1131,
//! `Sense` at 1201, `IsSplitToReverse` at 1255 / 1316 / 1432,
//! `GetEdgeOnFace` at 1809, `FindFacePairs` at 1839, `AngleWithRef` at 1938,
//! `GetFaceDir` at 2110, `FindPointInFace` at 2160, `MinStep3D` at 2235).

use std::collections::HashSet;

use occt_core::gp::{GpDir, GpPnt, GpVec};
use occt_core::precision::{ANGULAR, CONFUSION, RESOLUTION};

use crate::abs::{Orientation, ShapeType};
use crate::algo_tools::AlgoTools;
use crate::algo_tools3d::{
    edge_tangent_3d, get_approx_normal, get_normal_to_face_on_edge,
    get_normal_to_face_on_edge_at, get_normal_to_surface, point_in_face, surface_ref_radius,
};
use crate::boptools_2d::intermediate_point;
use crate::brep_tool::BRepTool;
use crate::fclass2d::FaceState;
use crate::int_tools_full::IntToolsContext;
use crate::iterator::ShapeIterator;
use crate::shape::{Edge, Face, TopoShape};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::edges_of;

fn shape_key(s: &TopoShape) -> usize {
    GeometryRegistry::shape_key(s)
}

/// Edge/face couple (`BOPTools_CoupleOfShape`).
#[derive(Clone)]
pub struct CoupleOfShape {
    pub shape1: TopoShape,
    pub shape2: TopoShape,
}

impl CoupleOfShape {
    pub fn new(shape1: TopoShape, shape2: TopoShape) -> Self {
        Self { shape1, shape2 }
    }
}

fn is_closed_on_face(edge: &Edge, face: &Face) -> bool {
    let mut n = 0usize;
    for w in ShapeIterator::of_shape(&face.0) {
        if w.shape_type() != ShapeType::Wire {
            continue;
        }
        for e in ShapeIterator::of_shape(&w) {
            if e.shape_type() == ShapeType::Edge && e.same_tshape(&edge.0) {
                n += 1;
            }
        }
    }
    n >= 2
}

/// `BOPTools_AlgoTools::GetEdgeOnFace`.
pub fn get_edge_on_face(the_e1: &Edge, the_f2: &Face) -> Option<Edge> {
    for w in ShapeIterator::of_shape(&the_f2.0) {
        if w.shape_type() != ShapeType::Wire {
            continue;
        }
        for e in ShapeIterator::of_shape(&w) {
            if e.shape_type() == ShapeType::Edge && e.same_tshape(&the_e1.0) {
                return Some(Edge(e));
            }
        }
    }
    None
}

/// `BOPTools_AlgoTools::GetEdgeOff`.
pub fn get_edge_off(the_e1: &Edge, the_f2: &Face) -> Option<Edge> {
    let or1c = the_e1.0.orientation().reversed();
    for w in ShapeIterator::of_shape(&the_f2.0) {
        if w.shape_type() != ShapeType::Wire {
            continue;
        }
        for e in ShapeIterator::of_shape(&w) {
            if e.shape_type() == ShapeType::Edge
                && e.same_tshape(&the_e1.0)
                && e.orientation() == or1c
            {
                return Some(Edge(e));
            }
        }
    }
    None
}

/// `AngleWithRef`.
fn angle_with_ref(d1: &GpDir, d2: &GpDir, d_ref: &GpDir) -> f64 {
    let xyz = d1.xyz().crossed(d2.xyz());
    let sinus = xyz.modulus();
    let cosinus = d1.dot(d2);
    let mut beta = sinus.atan2(cosinus);
    let sc = xyz.dot(d_ref.xyz());
    if sc < 0.0 {
        beta = -beta;
    }
    let _ = std::f64::consts::FRAC_PI_2;
    beta
}

fn project_on_plane(p: &GpPnt, origin: &GpPnt, n: &GpDir) -> GpPnt {
    let vx = p.x() - origin.x();
    let vy = p.y() - origin.y();
    let vz = p.z() - origin.z();
    let d = n.xyz();
    let dist = vx * d.x + vy * d.y + vz * d.z;
    GpPnt::new(p.x() - dist * d.x, p.y() - dist * d.y, p.z() - dist * d.z)
}

/// `FindPointInFace`.
fn find_point_in_face(
    face: &Face,
    p: &GpPnt,
    db: &mut GpDir,
    ctx: &IntToolsContext,
    origin: &GpPnt,
    n_pl: &GpDir,
    dt: f64,
    tol_e: f64,
) -> Option<GpPnt> {
    let Some(surf) = BRepTool::face_surface(face) else {
        return None;
    };
    let mut d_tol = ANGULAR;
    let pm = (p.x() * p.x() + p.y() * p.y() + p.z() * p.z()).sqrt();
    if pm > 1000.0 {
        d_tol = 5.0e-16 * pm;
    }
    let mut ps = *p;
    let Ok((u, v)) = ctx.project_point_on_face(face, &ps) else {
        return None;
    };
    ps = surf.d0(u, v);
    ps = project_on_plane(&ps, origin, n_pl);
    let dbv = GpVec::new(db.x(), db.y(), db.z());
    ps = ps.translated_vec(&dbv.multiplied_scalar(2.0 * tol_e));
    let Ok((u2, v2)) = ctx.project_point_on_face(face, &ps) else {
        return None;
    };
    ps = project_on_plane(&surf.d0(u2, v2), origin, n_pl);
    let mut nb = 15i32;
    let eps = CONFUSION * CONFUSION;
    loop {
        let p1 = ps.translated_vec(&GpVec::new(db.x() * dt, db.y() * dt, db.z() * dt));
        let Ok((u3, v3)) = ctx.project_point_on_face(face, &p1) else {
            return None;
        };
        let mut p_out = surf.d0(u3, v3);
        let dist = p1.distance(&p_out);
        p_out = project_on_plane(&p_out, origin, n_pl);
        let vec = GpVec::from_pnts(&ps, &p_out);
        if vec.square_magnitude() < eps {
            return None;
        }
        *db = GpDir::from_vec(&vec).ok()?;
        nb -= 1;
        if dist <= d_tol || nb <= 0 {
            return if dist < d_tol { Some(p_out) } else { None };
        }
        // `FindPointInFace` (`BOPTools_AlgoTools.cxx:2206-2227`) never advances
        // `aPS`: every iteration measures the offset from the *same* near-edge
        // point, so a bi-normal that leaves the face makes the offset collapse
        // (`aV.SquareMagnitude() < anEps` → the function returns false and
        // `GetFaceDir` takes the `GetApproxNormalToFaceOnEdge` fallback).
        // Advancing `aPS` instead lets the walk converge on a direction that
        // points *out* of the face, which is what made `GetFaceOff` choose the
        // other cylinder band instead of the disc (T-82).
        let _ = &p_out;
    }
}

/// `MinStep3D`.
fn min_step_3d(
    e1: &Edge,
    f1: &Face,
    lcs: &[CoupleOfShape],
    p: &GpPnt,
    ctx: &IntToolsContext,
    small: &mut bool,
) -> f64 {
    let mut all: Vec<CoupleOfShape> = lcs.to_vec();
    all.push(CoupleOfShape::new(e1.0.clone(), f1.0.clone()));
    let tol_e = BRepTool::edge_tolerance(e1);
    let mut dt_max = -1.0;
    let mut dt_min: f64 = 5.0e-6;
    for cs in &all {
        let f = Face(cs.shape2.clone());
        let tol_f = BRepTool::face_tolerance(&f);
        let dt = 2.0 * (tol_e + tol_f);
        if dt > dt_max {
            dt_max = dt;
        }
        let Some(surf) = BRepTool::face_surface(&f) else {
            continue;
        };
        let (kind, r) = surface_ref_radius(surf.as_ref(), p);
        match kind {
            crate::brep_surface::SurfaceKind::Sphere => {
                dt_min = dt_min.max(5.0e-4);
            }
            crate::brep_surface::SurfaceKind::Other
            | crate::brep_surface::SurfaceKind::Plane => {
                dt_min = dt_min.max(5.0e-4);
            }
            _ => {}
        }
        if r > 100.0 {
            let d = 10.0 * occt_core::precision::PCONFUSION;
            dt_min = dt_min.max((d * d + 2.0 * d * r).sqrt());
        }
    }
    if dt_max < dt_min {
        dt_max = dt_min;
    }
    *small = false;
    for cs in &all {
        let f = Face(cs.shape2.clone());
        let Some(surf) = BRepTool::face_surface(&f) else {
            continue;
        };
        let (umin, umax, vmin, vmax) = ctx.uv_bounds(&f);
        let du = umax - umin;
        let dv = vmax - vmin;
        let (u0, u1) = surf.u_range();
        let (v0, v1) = surf.v_range();
        let u_span = if du > 0.0 { du } else { u1 - u0 };
        let v_span = if dv > 0.0 { dv } else { v1 - v0 };
        let ures = if u_span.abs() > 0.0 {
            dt_max / u_span.abs().max(1.0)
        } else {
            0.0
        };
        let vres = if v_span.abs() > 0.0 {
            dt_max / v_span.abs().max(1.0)
        } else {
            0.0
        };
        if (du > 0.0 && 2.0 * ures > du) || (dv > 0.0 && 2.0 * vres > dv) {
            *small = true;
            break;
        }
        let _ = surf;
    }
    dt_max
}

/// `GetFaceDir`.
fn get_face_dir(
    edge: &Edge,
    face: &Face,
    p: &GpPnt,
    t: f64,
    dtgt: &GpDir,
    small: bool,
    ctx: &IntToolsContext,
    origin: &GpPnt,
    n_pl: &GpDir,
    dt: f64,
) -> Option<(GpDir, GpDir)> {
    let mut dn = get_normal_to_face_on_edge_at(edge, face, t, ctx)?;
    if face.0.orientation() == Orientation::Reversed {
        dn = dn.reversed();
    }
    let tol_e = BRepTool::edge_tolerance(edge);
    let mut db = dn.crossed(dtgt).ok()?;
    let found = if !small {
        find_point_in_face(face, p, &mut db, ctx, origin, n_pl, dt, tol_e).is_some()
    } else {
        false
    };
    if !found {
        let (px, dn2) = get_approx_normal(edge, face, t, ctx)?;
        let pxp = project_on_plane(&px, origin, n_pl);
        let vec = GpVec::from_pnts(p, &pxp);
        db = GpDir::from_vec(&vec).ok()?;
        return Some((dn2, db));
    }
    Some((dn, db))
}

/// `BOPTools_AlgoTools::GetFaceOff`.
pub fn get_face_off(
    e1: &Edge,
    f1: &Face,
    lcs_off: &[CoupleOfShape],
    ctx: &IntToolsContext,
) -> Option<(Face, bool)> {
    let Some(c3d) = BRepTool::edge_curve(e1) else {
        return None;
    };
    let (t1, t2) = BRepTool::edge_parameters(e1);
    let t = intermediate_point(t1, t2);
    let px = c3d.d0(t);
    let Some(vtgt) = edge_tangent_3d(e1, t) else {
        return None;
    };
    let dtgt = GpDir::from_vec(&vtgt).ok()?;
    let or1 = e1.0.orientation();
    let mut small = false;
    let dt3d = min_step_3d(e1, f1, lcs_off, &px, ctx, &mut small);
    let (dn1, dbf) = get_face_dir(e1, f1, &px, t, &dtgt, small, ctx, &px, &dtgt, dt3d)?;
    let dtf = dn1.crossed(&dbf).ok()?;
    let criteria = CONFUSION;
    let mut ok = true;
    let mut angle_min = 100.0;
    let two_pi = std::f64::consts::PI + std::f64::consts::PI;
    let mut f_off: Option<Face> = None;
    for cs in lcs_off {
        let e2 = Edge(cs.shape1.clone());
        let f2 = Face(cs.shape2.clone());
        // `GetFaceOff` (`BOPTools_AlgoTools.cxx:1051`) keeps `aDTgt` when the
        // two edges carry the same orientation, because they are the same
        // `TopoDS_Edge` and therefore traverse in the same direction. Here the
        // candidate edge is matched geometrically and may have been built with
        // the opposite vertex order, so the traversals themselves are compared.
        let same_dir = match (
            crate::shell_splitter::edge_traversal(e1),
            crate::shell_splitter::edge_traversal(&e2),
        ) {
            (Some(t1), Some(t2)) => t1 == t2,
            _ => e2.0.orientation() == or1,
        };
        let dtgt2 = if same_dir { dtgt } else { dtgt.reversed() };
        let Some((_dn2, dbf2)) =
            get_face_dir(&e2, &f2, &px, t, &dtgt2, small, ctx, &px, &dtgt, dt3d)
        else {
            continue;
        };
        let mut angle = angle_with_ref(&dbf, &dbf2, &dtf);
        if angle.abs() < ANGULAR {
            if f2.0.same_tshape(&f1.0) && f2.0.orientation() == f1.0.orientation() {
                angle = std::f64::consts::PI;
            } else if f2.0.same_tshape(&f1.0) {
                angle = two_pi;
            }
        }
        if angle.abs() < criteria || (angle - angle_min).abs() < criteria {
            ok = false;
        }
        if angle < 0.0 {
            angle += two_pi;
        }
        if angle < angle_min {
            angle_min = angle;
            f_off = Some(f2);
        }
    }
    f_off.map(|f| (f, ok))
}

fn find_face_pairs(
    the_e: &Edge,
    the_lf: &[TopoShape],
    ctx: &IntToolsContext,
) -> Vec<(TopoShape, TopoShape)> {
    let mut lcef: Vec<CoupleOfShape> = Vec::new();
    for fl in the_lf {
        let f = Face(fl.clone());
        let Some(el) = get_edge_on_face(the_e, &f) else {
            return Vec::new();
        };
        lcef.push(CoupleOfShape::new(el.0, fl.clone()));
    }
    let mut lcff: Vec<(TopoShape, TopoShape)> = Vec::new();
    let mut fence: HashSet<usize> = HashSet::new();
    while !lcef.is_empty() {
        let mut lcefx: Vec<CoupleOfShape> = Vec::new();
        let mut e1: Option<Edge> = None;
        let mut f1: Option<Face> = None;
        let mut or_c = Orientation::Forward;
        for (i, csx) in lcef.iter().enumerate() {
            let or = csx.shape1.orientation();
            if i == 0 {
                or_c = or.reversed();
                e1 = Some(Edge(csx.shape1.clone()));
                f1 = Some(Face(csx.shape2.clone()));
                fence.insert(shape_key(&csx.shape2));
                continue;
            }
            if or == or_c {
                lcefx.push(csx.clone());
                fence.insert(shape_key(&csx.shape2));
            }
        }
        let (Some(e1), Some(f1)) = (e1, f1) else {
            break;
        };
        let f2 = match get_face_off(&e1, &f1, &lcefx, ctx) {
            Some((f, _)) => f,
            None => break,
        };
        lcff.push((f1.0.clone(), f2.0.clone()));
        fence.insert(shape_key(&f1.0));
        fence.insert(shape_key(&f2.0));
        lcef.retain(|cs| !fence.contains(&shape_key(&cs.shape2)));
    }
    lcff
}

/// `BOPTools_AlgoTools::IsInternalFace` (face, edge, face1, face2).
pub fn is_internal_face_by_pair(
    the_face: &Face,
    the_edge: &Edge,
    face1: &Face,
    face2: &Face,
    ctx: &IntToolsContext,
) -> i32 {
    let Some(mut e1) = get_edge_on_face(the_edge, face1) else {
        return 2;
    };
    let e2 = if e1.0.orientation() == Orientation::Internal {
        let mut e = e1.clone();
        e1.0.set_orientation(Orientation::Forward);
        e.0.set_orientation(Orientation::Reversed);
        e
    } else if face1.0.same_tshape(&face2.0) {
        let mut e = e1.clone();
        e1.0.set_orientation(Orientation::Forward);
        e.0.set_orientation(Orientation::Reversed);
        e
    } else {
        match get_edge_on_face(the_edge, face2) {
            Some(e) => e,
            None => return 2,
        }
    };
    let lcs = vec![
        CoupleOfShape::new(the_edge.0.clone(), the_face.0.clone()),
        CoupleOfShape::new(e2.0, face2.0.clone()),
    ];
    // `GetFaceOff` false: angles are not distinct (`_cxx:977-982`) → iRet=2.
    match get_face_off(&e1, face1, &lcs, ctx) {
        None => 2,
        Some((_, false)) => 2,
        Some((f_off, true)) => {
            if the_face.0.same_tshape(&f_off.0) && the_face.0.orientation() == f_off.0.orientation()
            {
                1
            } else {
                0
            }
        }
    }
}

/// `BOPTools_AlgoTools::IsInternalFace` (face, edge, list of faces).
pub fn is_internal_face_by_list(
    the_face: &Face,
    the_edge: &Edge,
    the_lf: &[TopoShape],
    ctx: &IntToolsContext,
) -> i32 {
    if the_lf.len() == 2 {
        return is_internal_face_by_pair(
            the_face,
            the_edge,
            &Face(the_lf[0].clone()),
            &Face(the_lf[1].clone()),
            ctx,
        );
    }
    let pairs = find_face_pairs(the_edge, the_lf, ctx);
    for (f1, f2) in pairs {
        let i = is_internal_face_by_pair(the_face, the_edge, &Face(f1), &Face(f2), ctx);
        if i != 0 {
            return i;
        }
    }
    0
}

/// `BOPTools_AlgoTools::IsInternalFace` (face, solid, MEF, tol, context).
///
/// The angle-method early return matches OCCT `aExp.More() && iRet != 2`:
/// only when the edge explorer actually `break`s after a classification.
/// An MEF hit that is skipped (one adjacent face whose edge is not INTERNAL)
/// must fall through to `ComputeState`, not report Out.
pub fn is_internal_face(
    the_face: &Face,
    the_solid: &TopoShape,
    the_mef: &std::collections::HashMap<usize, (TopoShape, Vec<TopoShape>)>,
    the_tol: f64,
    the_ctx: &IntToolsContext,
) -> bool {
    let mut i_ret = 0i32;
    let mut classified = false;
    for e in edges_of(&the_face.0) {
        if e.0.orientation() == Orientation::Internal {
            continue;
        }
        if BRepTool::is_degenerated(&e) {
            continue;
        }
        let Some((_ek, lf)) = the_mef.get(&shape_key(&e.0)) else {
            continue;
        };
        if lf.len() == 1 {
            let f1 = Face(lf[0].clone());
            let Some(e1) = get_edge_on_face(&e, &f1) else {
                continue;
            };
            if e1.0.orientation() != Orientation::Internal {
                continue;
            }
            i_ret = is_internal_face_by_pair(the_face, &e, &f1, &f1, the_ctx);
            classified = true;
            break;
        } else if lf.len() == 2 {
            i_ret = is_internal_face_by_pair(
                the_face,
                &e,
                &Face(lf[0].clone()),
                &Face(lf[1].clone()),
                the_ctx,
            );
            if i_ret != 2 {
                classified = true;
                break;
            }
        }
    }
    if classified && i_ret != 2 {
        return i_ret == 1;
    }
    compute_state_face_in_solid(the_face, the_solid, the_tol)
}

/// `BOPTools_AlgoTools::ComputeState` (face, solid, bounds, context).
///
/// Prefer the midpoint of an edge that is not a bound of the solid; if every
/// edge sits on the solid, classify a point in the face.
fn compute_state_face_in_solid(the_face: &Face, the_solid: &TopoShape, the_tol: f64) -> bool {
    let bounds: HashSet<usize> = edges_of(the_solid)
        .into_iter()
        .map(|e| shape_key(&e.0))
        .collect();
    for e in edges_of(&the_face.0) {
        if BRepTool::is_degenerated(&e) {
            continue;
        }
        if bounds.contains(&shape_key(&e.0)) {
            continue;
        }
        let p = edge_state_point(&e);
        return match AlgoTools::compute_state(the_solid, &p, the_tol) {
            Ok(st) => st == FaceState::In,
            _ => false,
        };
    }
    let p = match point_in_face(the_face, &mut IntToolsContext::new()) {
        Ok((p, _)) => p,
        Err(_) => {
            let Some(surf) = BRepTool::face_surface(the_face) else {
                return false;
            };
            let (u0, u1) = surf.u_range();
            let (v0, v1) = surf.v_range();
            surf.d0(intermediate_point(u0, u1), intermediate_point(v0, v1))
        }
    };
    let st = match AlgoTools::compute_state(the_solid, &p, the_tol) {
        Ok(s) => s,
        _ => return false,
    };
    st == FaceState::In
}

fn edge_state_point(e: &Edge) -> GpPnt {
    let (a, b) = BRepTool::edge_parameters(e);
    if let Some(c) = BRepTool::edge_curve(e) {
        let t = if !a.is_finite() && b.is_finite() {
            b - 10.0
        } else if a.is_finite() && !b.is_finite() {
            a + 10.0
        } else if !a.is_finite() && !b.is_finite() {
            0.0
        } else {
            intermediate_point(a, b)
        };
        return c.d0(t);
    }
    crate::topo_tools_full::vertices_of(&e.0)
        .first()
        .map(BRepTool::vertex_point)
        .unwrap_or_else(GpPnt::zero)
}

/// `BOPTools_AlgoTools::AreFacesSameDomain`.
pub fn are_faces_same_domain(
    f1: &Face,
    f2: &Face,
    ctx: &mut IntToolsContext,
    fuzz: f64,
) -> bool {
    let Ok((p1, _)) = point_in_face(f1, ctx) else {
        return false;
    };
    let mut tol_f1 = BRepTool::face_tolerance(f1);
    let mut tol_f2 = BRepTool::face_tolerance(f2);
    let mut tol_e_max = -1.0;
    for e in edges_of(&f1.0) {
        if BRepTool::is_degenerated(&e) {
            continue;
        }
        let t = BRepTool::edge_tolerance(&e);
        if t > tol_e_max {
            tol_e_max = t;
        }
    }
    if tol_e_max > tol_f1 {
        tol_f1 = tol_e_max;
    }
    if tol_e_max > tol_f2 {
        tol_f2 = tol_e_max;
    }
    let tol = tol_f1 + tol_f2 + fuzz.max(CONFUSION);
    let Ok((u, v)) = ctx.project_point_on_face(f2, &p1) else {
        return false;
    };
    let Some(s2) = BRepTool::face_surface(f2) else {
        return false;
    };
    if p1.distance(&s2.d0(u, v)) > tol {
        return false;
    }
    ctx.is_valid_point_for_face((u, v), f2).unwrap_or(false)
}

/// `BOPTools_AlgoTools::Sense`.
pub fn sense_faces(f1: &Face, f2: &Face, ctx: &IntToolsContext) -> i32 {
    let mut e1: Option<Edge> = None;
    for e in edges_of(&f1.0) {
        if BRepTool::is_degenerated(&e) {
            continue;
        }
        if !is_closed_on_face(&e, f1) {
            e1 = Some(e);
            break;
        }
    }
    let Some(e1) = e1 else {
        return 0;
    };
    let mut e2: Option<Edge> = None;
    for e in edges_of(&f2.0) {
        if BRepTool::is_degenerated(&e) {
            continue;
        }
        if !is_closed_on_face(&e, f2) && e.0.same_tshape(&e1.0) {
            e2 = Some(e);
            break;
        }
    }
    let Some(e2) = e2 else {
        return 0;
    };
    let Some(n1) = get_normal_to_face_on_edge(&e1, f1, ctx) else {
        return 0;
    };
    let Some(n2) = get_normal_to_face_on_edge(&e2, f2, ctx) else {
        return 0;
    };
    crate::algo_tools3d::sense_flag(&n1, &n2)
}

/// `BOPTools_AlgoTools::IsSplitToReverse` (faces).
pub fn is_split_to_reverse_face(
    f_sp: &Face,
    f_sr: &Face,
    ctx: &mut IntToolsContext,
) -> Result<bool, i32> {
    let Some(s_sp) = BRepTool::face_surface(f_sp) else {
        return Err(2);
    };
    let Some(s_sr) = BRepTool::face_surface(f_sr) else {
        return Err(2);
    };
    if GeometryRegistry::shape_key(&f_sp.0) == GeometryRegistry::shape_key(&f_sr.0) {
        return Ok(f_sp.0.orientation() != f_sr.0.orientation());
    }
    let (p_sp, p2d) = match point_in_face(f_sp, ctx) {
        Ok(v) => v,
        Err(_) => {
            let mut found = None;
            for e in edges_of(&f_sp.0) {
                if BRepTool::is_degenerated(&e) || is_closed_on_face(&e, f_sp) {
                    continue;
                }
                let (t1, t2) = BRepTool::edge_parameters(&e);
                let t = intermediate_point(t1, t2);
                if let Ok(v) = crate::algo_tools3d::point_near_edge(&e, f_sp, t, ctx) {
                    found = Some((v.1, v.0));
                    break;
                }
            }
            found.ok_or(1)?
        }
    };
    let mut dn_sp = get_normal_to_surface(s_sp.as_ref(), p2d.x(), p2d.y()).ok_or(2)?;
    if f_sp.0.orientation() == Orientation::Reversed {
        dn_sp = dn_sp.reversed();
    }
    let (u, v) = ctx.project_point_on_face(f_sr, &p_sp).map_err(|_| 3)?;
    let mut dn_or = get_normal_to_surface(s_sr.as_ref(), u, v).ok_or(4)?;
    if f_sr.0.orientation() == Orientation::Reversed {
        dn_or = dn_or.reversed();
    }
    Ok(dn_sp.dot(&dn_or) < 0.0)
}

/// `BOPTools_AlgoTools::IsSplitToReverse` (edges)
/// (`BOPTools_AlgoTools.cxx:1432-1500`).
///
/// Substitution registered for audit §13 / task T-34: OCCT narrows the sampled
/// parameter range with `BRepLib::FindValidRange(theESp, f, l)` and only falls
/// back to `BRep_Tool::Range` when that returns false (`cxx:1469-1472`); this
/// port uses [`BRepTool::edge_parameters`] (`BRep_Tool::Range`) directly. For an
/// edge whose stored range extends beyond the valid part of its curve the
/// sampled `tm` therefore differ, so the two implementations are **not**
/// line-by-line equivalent.
pub fn is_split_to_reverse_edge(
    e_sp: &Edge,
    e_or: &Edge,
    ctx: &IntToolsContext,
) -> Result<bool, i32> {
    if BRepTool::is_degenerated(e_sp) || BRepTool::is_degenerated(e_or) {
        return Err(1);
    }
    let Some(c_sp) = BRepTool::edge_curve(e_sp) else {
        return Err(1);
    };
    let Some(c_or) = BRepTool::edge_curve(e_or) else {
        return Err(1);
    };
    // `BOPTools_AlgoTools.cxx:1461-1465`: when the two edges carry the *same
    // curve* — every split piece of an edge does, `MakeSplitEdge` copies the
    // parent's `Geom_Curve` handle — the decision is the orientation comparison
    // alone, and the tangent test is not run. Comparing the edge TShape instead
    // would never fire for a split piece (it is a distinct TShape) and would
    // force the geometric branch on the most common case.
    if std::sync::Arc::ptr_eq(&c_sp, &c_or) {
        return Ok(e_sp.0.orientation() != e_or.0.orientation());
    }
    let (f, l) = BRepTool::edge_parameters(e_sp);
    const NB: i32 = 11;
    let dt = (l - f) / NB as f64;
    let mut err = 0;
    for i in 1..NB {
        let tm = f + i as f64 * dt;
        let Some(v_sp) = edge_tangent_3d(e_sp, tm) else {
            err = 2;
            continue;
        };
        let p = c_sp.d0(tm);
        let Some(tm_or) = ctx.project_point_on_edge(e_or, &p) else {
            err = 3;
            continue;
        };
        let Some(v_or) = edge_tangent_3d(e_or, tm_or) else {
            err = 4;
            continue;
        };
        let _ = err;
        return Ok(v_sp.dot(&v_or) < 0.0);
    }
    Err(err)
}

/// `BOPTools_AlgoTools::IsSplitToReverse` (shape).
pub fn is_split_to_reverse(
    sp: &TopoShape,
    sr: &TopoShape,
    ctx: &mut IntToolsContext,
) -> Result<bool, i32> {
    match sp.shape_type() {
        ShapeType::Edge => is_split_to_reverse_edge(&Edge(sp.clone()), &Edge(sr.clone()), ctx),
        ShapeType::Face => is_split_to_reverse_face(&Face(sp.clone()), &Face(sr.clone()), ctx),
        _ => Err(100),
    }
}
