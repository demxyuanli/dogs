//! `BOPTools_AlgoTools3D` — face normals, point-in-face, point-near-edge.
//!
//! Source: `BOPTools_AlgoTools3D.cxx` (`GetNormalToFaceOnEdge` at 331 and 351,
//! `SenseFlag` at 380, `GetNormalToSurface` at 406, `GetApproxNormalToFaceOnEdge`
//! at 443 and 469, `PointNearEdge` at 525, `PointInFace` at 906 / 942 / 992).
//! The 2-D hatcher (`Geom2dHatch_Hatcher`) is the `FClass2d` walk of a 2-D
//! line across the face UV box: the first IN parameter interval is the hatch
//! domain, and `IntermediatePoint` of that interval is the interior point.

use occt_core::gp::{GpDir, GpDir2d, GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_core::precision::{ANGULAR, CONFUSION, PCONFUSION, RESOLUTION};
use occt_geom::Surface;

use crate::abs::Orientation;
use crate::boptools_2d::{curve_on_surface, intermediate_point, make_2d};
use crate::brep_surface::{classify_surface, SurfaceKind};
use crate::brep_tool::BRepTool;
use crate::fclass2d::FaceState;
use crate::int_tools_full::IntToolsContext;
use crate::shape::{Edge, Face};

/// 3-D unit tangent of `edge` at `t` (`BOPTools_AlgoTools2D::EdgeTangent`).
pub fn edge_tangent_3d(edge: &Edge, t: f64) -> Option<GpVec> {
    let curve = BRepTool::edge_curve(edge)?;
    let (_, mut tau) = curve.d1(t);
    let mag = tau.magnitude();
    if mag <= RESOLUTION {
        return None;
    }
    tau = tau.divided(mag);
    if edge.0.orientation() == Orientation::Reversed {
        tau = tau.reversed();
    }
    Some(tau)
}

/// `BOPTools_AlgoTools3D::GetNormalToSurface`.
pub fn get_normal_to_surface(s: &dyn Surface, u: f64, v: f64) -> Option<GpDir> {
    let (_p, d1u, d1v) = s.d1(u, v);
    if d1u.square_magnitude() < RESOLUTION || d1v.square_magnitude() < RESOLUTION {
        return None;
    }
    let du = GpDir::from_vec(&d1u).ok()?;
    let dv = GpDir::from_vec(&d1v).ok()?;
    du.crossed(&dv).ok()
}

/// `BOPTools_AlgoTools3D::GetNormalToFaceOnEdge` (edge, face, parameter).
pub fn get_normal_to_face_on_edge_at(
    edge: &Edge,
    face: &Face,
    t: f64,
    ctx: &IntToolsContext,
) -> Option<GpDir> {
    let _ = ctx;
    let pc = match curve_on_surface(edge, face) {
        Some(pc) => pc,
        None => make_2d(edge, face).ok()?,
    };
    let p2d = pc.d0(t);
    let surf = BRepTool::face_surface(face)?;
    get_normal_to_surface(surf.as_ref(), p2d.x(), p2d.y())
}

/// `BOPTools_AlgoTools3D::GetNormalToFaceOnEdge` (edge, face) at the
/// intermediate p-curve parameter, reversed when the face is `REVERSED`.
pub fn get_normal_to_face_on_edge(
    edge: &Edge,
    face: &Face,
    ctx: &IntToolsContext,
) -> Option<GpDir> {
    let (t1, t2) = match crate::boptools_2d::curve_on_surface_range(edge, face) {
        Some((_, a, b)) => (a, b),
        None => BRepTool::edge_parameters(edge),
    };
    let t = intermediate_point(t1, t2);
    let mut dn = get_normal_to_face_on_edge_at(edge, face, t, ctx)?;
    if face.0.orientation() == Orientation::Reversed {
        dn = dn.reversed();
    }
    Some(dn)
}

/// `BOPTools_AlgoTools3D::SenseFlag`.
pub fn sense_flag(d1: &GpDir, d2: &GpDir) -> i32 {
    let ang = d1.angle(d2);
    if ang > ANGULAR && (std::f64::consts::PI - ang) > ANGULAR {
        return 0;
    }
    let pr = d1.dot(d2);
    if pr < 0.0 {
        -1
    } else if pr > 0.0 {
        1
    } else {
        -1
    }
}

/// `BOPTools_AlgoTools3D::PointNearEdge` (`aT`, `aDt2D`).
pub fn point_near_edge_dt(
    edge: &Edge,
    face: &Face,
    t: f64,
    dt2d: f64,
) -> Result<(GpPnt2d, GpPnt), i32> {
    let pc = match curve_on_surface(edge, face) {
        Some(pc) => pc,
        None => match make_2d(edge, face) {
            Ok(pc) => pc,
            Err(_) => return Err(1),
        },
    };
    let Some(surf) = BRepTool::face_surface(face) else {
        return Err(1);
    };
    let (p2d, v2d) = pc.d1(t);
    let mag = v2d.magnitude();
    if mag <= RESOLUTION {
        return Err(1);
    }
    let dx = GpDir2d::new(-v2d.y() / mag, v2d.x() / mag).map_err(|_| 1)?;
    let mut dp = dx;
    if edge.0.orientation() == Orientation::Reversed {
        dp = dp.reversed();
    }
    if face.0.orientation() == Orientation::Reversed {
        dp = dp.reversed();
    }
    let e_tol = BRepTool::edge_tolerance(edge);
    let mut f_tol = BRepTool::face_tolerance(face);
    let kind = classify_surface(surf.as_ref());
    if kind == SurfaceKind::Other && e_tol > 1.0e-5 {
        f_tol = e_tol;
    }
    let mut near = p2d;
    if e_tol > 1.0e-5 || f_tol > 1.0e-5 {
        if kind != SurfaceKind::Sphere {
            let mut trans = dt2d + e_tol + f_tol;
            if kind == SurfaceKind::Cylinder {
                if let Some(r) = cylinder_radius(surf.as_ref()) {
                    let dt = 1.0 - trans / r;
                    if (-1.0..=1.0).contains(&dt) {
                        trans = dt.acos();
                    }
                }
            }
            let tv = GpVec2d::new(dp.x() * trans, dp.y() * trans);
            near = p2d.translated_vec(&tv);
        } else {
            near.set_coord(p2d.x() + dt2d * dp.x(), p2d.y() + dt2d * dp.y());
        }
    } else {
        let tv = GpVec2d::new(dp.x() * dt2d, dp.y() * dt2d);
        near = p2d.translated_vec(&tv);
    }
    let p3d = surf.d0(near.x(), near.y());
    let _ = e_tol;
    Ok((near, p3d))
}

/// `BOPTools_AlgoTools3D::PointNearEdge` with context (default 2-D step).
pub fn point_near_edge(
    edge: &Edge,
    face: &Face,
    t: f64,
    ctx: &IntToolsContext,
) -> Result<(GpPnt2d, GpPnt), i32> {
    let _ = ctx;
    let dt2d = 16.0 * CONFUSION.max(BRepTool::edge_tolerance(edge) + BRepTool::face_tolerance(face));
    point_near_edge_dt(edge, face, t, dt2d)
}

/// `BOPTools_AlgoTools3D::GetApproxNormalToFaceOnEdge` with an explicit 2-D step.
pub fn get_approx_normal_dt(
    edge: &Edge,
    face: &Face,
    t: f64,
    dt2d: f64,
) -> Option<(GpPnt, GpDir)> {
    let (p2d, pnear) = match point_near_edge_dt(edge, face, t, dt2d) {
        Ok(v) => v,
        Err(1) => return None,
        Err(_) => return None,
    };
    let surf = BRepTool::face_surface(face)?;
    let mut dn = get_normal_to_surface(surf.as_ref(), p2d.x(), p2d.y())?;
    if face.0.orientation() == Orientation::Reversed {
        dn = dn.reversed();
    }
    Some((pnear, dn))
}

/// `BOPTools_AlgoTools3D::GetApproxNormalToFaceOnEdge` (context overload).
pub fn get_approx_normal(
    edge: &Edge,
    face: &Face,
    t: f64,
    ctx: &IntToolsContext,
) -> Option<(GpPnt, GpDir)> {
    let (p2d, pnear) = point_near_edge(edge, face, t, ctx).ok()?;
    let surf = BRepTool::face_surface(face)?;
    let mut dn = get_normal_to_surface(surf.as_ref(), p2d.x(), p2d.y())?;
    if face.0.orientation() == Orientation::Reversed {
        dn = dn.reversed();
    }
    Some((pnear, dn))
}

fn cylinder_radius(s: &dyn Surface) -> Option<f64> {
    let p0 = s.d0(0.0, 0.0);
    let p1 = s.d0(std::f64::consts::PI, 0.0);
    let r = p0.distance(&p1) / 2.0;
    if r > 1e-9 {
        Some(r)
    } else {
        None
    }
}

fn sphere_radius(s: &dyn Surface) -> Option<f64> {
    crate::brep_surface::sphere_center(s).map(|c| c.distance(&s.d0(0.0, 0.0)))
}

/// `BOPTools_AlgoTools3D::PointInFace` (face + 2-D curve).
///
/// The hatcher is the first IN interval of `line2d` classified by `FClass2d`
/// along the UV box of the face.
pub fn point_in_face_on_curve2d(
    face: &Face,
    eval: impl Fn(f64) -> GpPnt2d,
    t_first: f64,
    t_last: f64,
    ctx: &mut IntToolsContext,
    dt2d: f64,
) -> Result<(GpPnt, GpPnt2d), i32> {
    let Some(surf) = BRepTool::face_surface(face) else {
        return Err(1);
    };
    const NB: i32 = 32;
    let span = t_last - t_first;
    if span.abs() <= PCONFUSION {
        return Err(2);
    }
    let mut in_first: Option<f64> = None;
    let mut in_last: Option<f64> = None;
    for i in 0..=NB {
        let t = t_first + span * (i as f64) / (NB as f64);
        let uv = eval(t);
        let st = match ctx.state_point_face(face, (uv.x(), uv.y()), CONFUSION) {
            Ok(s) => s,
            Err(_) => continue,
        };
        if st == FaceState::In || st == FaceState::On {
            if in_first.is_none() {
                in_first = Some(t);
            }
            in_last = Some(t);
        } else if in_first.is_some() {
            break;
        }
    }
    let (v1, v2) = match (in_first, in_last) {
        (Some(a), Some(b)) if b > a + PCONFUSION => (a, b),
        _ => return Err(2),
    };
    let vx = if dt2d > 0.0 && (v2 - v1) > dt2d {
        v1 + dt2d
    } else {
        intermediate_point(v1, v2)
    };
    let p2d = eval(vx);
    let p = surf.d0(p2d.x(), p2d.y());
    Ok((p, p2d))
}

/// `BOPTools_AlgoTools3D::PointInFace` (face only).
pub fn point_in_face(face: &Face, ctx: &mut IntToolsContext) -> Result<(GpPnt, GpPnt2d), i32> {
    let (umin, umax, vmin, vmax) = ctx.uv_bounds(face);
    let mut ux = intermediate_point(umin, umax);
    for _pass in 0..2 {
        let eval = |t: f64| GpPnt2d::new(ux, t);
        match point_in_face_on_curve2d(face, eval, vmin, vmax, ctx, 0.0) {
            Ok(v) => return Ok(v),
            Err(_) => {
                ux = umax - (ux - umin);
            }
        }
    }
    Err(1)
}

/// `BOPTools_AlgoTools3D::PointInFace` (edge + parameter + 2-D offset).
pub fn point_in_face_from_edge(
    face: &Face,
    edge: &Edge,
    t: f64,
    dt2d: f64,
    ctx: &mut IntToolsContext,
) -> Result<(GpPnt, GpPnt2d), i32> {
    let pc = match curve_on_surface(edge, face) {
        Some(pc) => pc,
        None => return Err(5),
    };
    let (p2d, v2d) = pc.d1(t);
    let mag = v2d.magnitude();
    if mag <= RESOLUTION {
        return Err(1);
    }
    let mut d2d = GpDir2d::new(-v2d.y() / mag, v2d.x() / mag).map_err(|_| 1)?;
    if edge.0.orientation() == Orientation::Reversed {
        d2d = d2d.reversed();
    }
    if face.0.orientation() == Orientation::Reversed {
        d2d = d2d.reversed();
    }
    let eval = move |s: f64| {
        GpPnt2d::new(p2d.x() + s * d2d.x(), p2d.y() + s * d2d.y())
    };
    let span = dt2d.max(16.0 * CONFUSION);
    point_in_face_on_curve2d(face, eval, 0.0, span * 8.0, ctx, dt2d)
}

/// Surface radius used by `MinStep3D`.
pub fn surface_ref_radius(s: &dyn Surface, at: &GpPnt) -> (SurfaceKind, f64) {
    let kind = classify_surface(s);
    let r = match kind {
        SurfaceKind::Cylinder => cylinder_radius(s).unwrap_or(0.0),
        SurfaceKind::Sphere => sphere_radius(s).unwrap_or(0.0),
        SurfaceKind::Cone => {
            if let Some(c) = crate::brep_surface::sphere_center(s) {
                at.distance(&c)
            } else {
                0.0
            }
        }
        SurfaceKind::Torus => cylinder_radius(s).unwrap_or(0.0),
        _ => 0.0,
    };
    (kind, r)
}

/// `IntTools_Tools::IsDirsCoinside` used by `SenseFlag` / `GetNormalToSurface`.
pub fn dirs_coincide(a: &GpDir, b: &GpDir) -> bool {
    let ang = a.angle(b);
    ang <= ANGULAR || (std::f64::consts::PI - ang) <= ANGULAR
}
