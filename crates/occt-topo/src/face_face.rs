//! Face–face intersection and distance queries.
//!
//! Source: `BRepAlgoAPI_Section` / `IntAna_IntConicQuad` (TKBO). The plane–plane
//! case is exact (two-plane linear system). General surface–surface
//! intersection is an *approximate* grid sampler: it evaluates both surfaces on
//! parametric grids, keeps points of the first surface whose distance to the
//! second falls below the tolerance, and returns the sampled point set.
//! ponytail: general surface–surface intersection is a grid approximation; a
//! robust subdivision/exact algorithm (BRepAlgoAPI_Section) is the upgrade path.

use occt_core::gp::{GpPln, GpPnt, GpVec};
use occt_geom::Surface;

use crate::brep_tool::BRepTool;
use crate::shape::Face;

/// Result of intersecting two faces.
#[derive(Debug, Clone)]
pub enum FaceIntersect {
    /// No intersection found.
    None,
    /// Two non-coplanar planes: an infinite line `origin + t·dir`.
    Line { origin: GpPnt, dir: GpVec },
    /// General surface–surface: sampled points lying on both surfaces.
    Curve { points: Vec<GpPnt> },
    /// Two coincident planes.
    Coplanar,
}

/// Default half-extent used when a surface reports an unbounded (infinite)
/// parameter range, so it can still be sampled on a finite window.
const DEFAULT_WINDOW: f64 = 4.0;

/// Map an (unbounded) parameter range onto a finite sampling window.
pub(crate) fn finite_window(range: (f64, f64), default: f64) -> (f64, f64) {
    let (a, b) = range;
    if a.is_finite() && b.is_finite() && a < b {
        (a, b)
    } else if a.is_finite() && b.is_finite() {
        (a, a) // degenerate range
    } else {
        (-default, default)
    }
}

/// Unit normal estimate at `(u, v)` via finite differences of `d0`.
fn surface_normal_at(s: &dyn Surface, u: f64, v: f64) -> GpVec {
    let h = 1e-4;
    let p = s.d0(u, v);
    let pu = s.d0(u + h, v);
    let pv = s.d0(u, v + h);
    GpVec::from_pnts(&p, &pu).crossed(&GpVec::from_pnts(&p, &pv))
}

/// True when the surface's normal direction is constant over its patch, i.e.
/// it is geometrically a plane. Sampled on a 3×3 grid of the (windowed) domain.
pub(crate) fn is_plane_like(s: &dyn Surface) -> bool {
    let (u0, u1) = finite_window(s.u_range(), DEFAULT_WINDOW);
    let (v0, v1) = finite_window(s.v_range(), DEFAULT_WINDOW);
    let us = [u0, (u0 + u1) * 0.5, u1];
    let vs = [v0, (v0 + v1) * 0.5, v1];
    let base = surface_normal_at(s, us[0], vs[0]);
    if base.square_magnitude() <= 1e-24 {
        return false;
    }
    for &u in &us {
        for &v in &vs {
            let n = surface_normal_at(s, u, v);
            if n.square_magnitude() <= 1e-24 || base.angle(&n) > 1e-5 {
                return false;
            }
        }
    }
    true
}

/// A point on the surface plus its unit normal — the geometric plane data.
pub(crate) fn plane_geometry(s: &dyn Surface) -> (GpPnt, GpVec) {
    let (u0, u1) = finite_window(s.u_range(), DEFAULT_WINDOW);
    let (v0, v1) = finite_window(s.v_range(), DEFAULT_WINDOW);
    let (um, vm) = ((u0 + u1) * 0.5, (v0 + v1) * 0.5);
    let p = s.d0(um, vm);
    let n = surface_normal_at(s, um, vm);
    (p, n.divided(n.magnitude()))
}

/// Exact intersection line of two planes given by a point + unit normal.
/// Returns `None` when the planes are parallel (coplanar included).
fn intersect_plane_plane(p1: &GpPnt, n1: &GpVec, p2: &GpPnt, n2: &GpVec, tol: f64) -> Option<(GpPnt, GpVec)> {
    let d1 = n1.dot(&GpVec::from_pnts(&GpPnt::zero(), p1));
    let d2 = n2.dot(&GpVec::from_pnts(&GpPnt::zero(), p2));
    let dir = n1.crossed(n2);
    let denom = dir.square_magnitude();
    if denom <= tol * tol {
        return None;
    }
    // p0 = (d1·(n2×d) + d2·(d×n1)) / |d|²  — the point on both planes closest to origin.
    let p0v = n2
        .crossed(&dir)
        .multiplied_scalar(d1)
        .added(&dir.crossed(n1).multiplied_scalar(d2))
        .divided(denom);
    let p0 = GpPnt::from_xyz(&p0v.coord);
    Some((p0, dir.divided(dir.magnitude())))
}

/// Exact analytic intersection of two infinite planes (`IntAna_IntConicQuad`).
/// Returns `(origin, direction)` of the intersection line, or `None` when the
/// planes are parallel.
pub fn plane_plane_intersection(p1: &GpPln, p2: &GpPln) -> Option<(GpPnt, GpVec)> {
    let n1 = GpVec::from_xyz(p1.axis().direction().xyz());
    let n2 = GpVec::from_xyz(p2.axis().direction().xyz());
    intersect_plane_plane(&p1.location(), &n1, &p2.location(), &n2, 1e-9)
}

/// Approximate distance from `p` to the nearest point of surface `s`.
/// Coarse grid search followed by a few coordinate-descent refinement steps.
fn distance_to_surface(p: &GpPnt, s: &dyn Surface, nu: usize, nv: usize) -> f64 {
    let (u0, u1) = finite_window(s.u_range(), DEFAULT_WINDOW);
    let (v0, v1) = finite_window(s.v_range(), DEFAULT_WINDOW);
    let mut bu = u0;
    let mut bv = v0;
    let mut best = f64::INFINITY;
    for i in 0..=nu {
        let u = u0 + (u1 - u0) * i as f64 / nu as f64;
        for j in 0..=nv {
            let v = v0 + (v1 - v0) * j as f64 / nv as f64;
            let d2 = p.square_distance(&s.d0(u, v));
            if d2 < best {
                best = d2;
                bu = u;
                bv = v;
            }
        }
    }
    // Coordinate-descent refinement from the best grid node.
    let mut u = bu;
    let mut v = bv;
    let mut h = ((u1 - u0).max(v1 - v0)) / nu as f64;
    for _ in 0..48 {
        let f0 = p.square_distance(&s.d0(u, v));
        let cands = [
            (u + h, v),
            (u - h, v),
            (u, v + h),
            (u, v - h),
        ];
        let mut improved = false;
        let mut nu2 = u;
        let mut nv2 = v;
        for (cu, cv) in cands {
            let cf = p.square_distance(&s.d0(cu, cv));
            if cf < f0 {
                nu2 = cu;
                nv2 = cv;
                improved = true;
            }
        }
        if improved {
            u = nu2.clamp(u0, u1);
            v = nv2.clamp(v0, v1);
        } else {
            h *= 0.5;
        }
        if h < 1e-10 {
            break;
        }
    }
    p.square_distance(&s.d0(u, v)).sqrt()
}

/// Sample surface `s1` on a grid and keep the points whose (approximate)
/// distance to `s2` is within `tol`.
fn sample_intersection_points(s1: &dyn Surface, s2: &dyn Surface, tol: f64) -> Vec<GpPnt> {
    let (u0, u1) = finite_window(s1.u_range(), DEFAULT_WINDOW);
    let (v0, v1) = finite_window(s1.v_range(), DEFAULT_WINDOW);
    let nu = 64usize;
    let nv = 64usize;
    let mut pts = Vec::new();
    for i in 0..=nu {
        let u = u0 + (u1 - u0) * i as f64 / nu as f64;
        for j in 0..=nv {
            let v = v0 + (v1 - v0) * j as f64 / nv as f64;
            let p = s1.d0(u, v);
            if distance_to_surface(&p, s2, 32, 32) <= tol {
                pts.push(p);
            }
        }
    }
    pts
}

/// Intersect two faces. Plane–plane is exact; the general case is the
/// approximate grid sampler described at the top of this module.
pub fn face_face_intersection(f1: &Face, f2: &Face, tol: f64) -> FaceIntersect {
    let s1 = match BRepTool::face_surface(f1) {
        Some(s) => s,
        None => return FaceIntersect::None,
    };
    let s2 = match BRepTool::face_surface(f2) {
        Some(s) => s,
        None => return FaceIntersect::None,
    };
    if is_plane_like(&*s1) && is_plane_like(&*s2) {
        let (p1, n1) = plane_geometry(&*s1);
        let (p2, n2) = plane_geometry(&*s2);
        if n1.cross_magnitude(&n2) <= tol {
            // Parallel: coincident → coplanar, else no intersection.
            let dist = n1.dot(&GpVec::from_pnts(&p1, &p2)).abs();
            return if dist <= tol {
                FaceIntersect::Coplanar
            } else {
                FaceIntersect::None
            };
        }
        return match intersect_plane_plane(&p1, &n1, &p2, &n2, tol) {
            Some((origin, dir)) => FaceIntersect::Line { origin, dir },
            None => FaceIntersect::None,
        };
    }
    // ponytail: approximate grid sampling; upgrade to exact section later.
    let pts = sample_intersection_points(&*s1, &*s2, tol);
    if pts.is_empty() {
        FaceIntersect::None
    } else {
        FaceIntersect::Curve { points: pts }
    }
}

/// Grid-search for the `(u, v)` parameter of surface `s` nearest to `p`,
/// returning `Some` only when within `tol`.
pub fn point_on_surface(s: &dyn Surface, p: &GpPnt, tol: f64) -> Option<(f64, f64)> {
    let (u0, u1) = finite_window(s.u_range(), 10.0);
    let (v0, v1) = finite_window(s.v_range(), 10.0);
    let nu = 64usize;
    let nv = 64usize;
    let mut best: Option<(f64, f64, f64)> = None; // (distance², u, v)
    for i in 0..=nu {
        let u = u0 + (u1 - u0) * i as f64 / nu as f64;
        for j in 0..=nv {
            let v = v0 + (v1 - v0) * j as f64 / nv as f64;
            let d2 = p.square_distance(&s.d0(u, v));
            if best.map_or(true, |(b, _, _)| d2 < b) {
                best = Some((d2, u, v));
            }
        }
    }
    let (d2, u, v) = best?;
    if d2.sqrt() <= tol {
        Some((u, v))
    } else {
        None
    }
}

/// Minimum distance between the sampled point sets of two faces. `tol` is
/// currently unused (kept for API symmetry); returns `inf` when either face
/// lacks a registered surface.
pub fn face_distance(f1: &Face, f2: &Face, _tol: f64) -> f64 {
    let s1 = match BRepTool::face_surface(f1) {
        Some(s) => s,
        None => return f64::INFINITY,
    };
    let s2 = match BRepTool::face_surface(f2) {
        Some(s) => s,
        None => return f64::INFINITY,
    };
    let (u0, u1) = finite_window(s1.u_range(), DEFAULT_WINDOW);
    let (v0, v1) = finite_window(s1.v_range(), DEFAULT_WINDOW);
    let (r0, r1) = finite_window(s2.u_range(), DEFAULT_WINDOW);
    let (t0, t1) = finite_window(s2.v_range(), DEFAULT_WINDOW);
    let mut best = f64::INFINITY;
    for i in 0..=16 {
        let u = u0 + (u1 - u0) * i as f64 / 16.0;
        for j in 0..=16 {
            let v = v0 + (v1 - v0) * j as f64 / 16.0;
            let p = s1.d0(u, v);
            for k in 0..=16 {
                let r = r0 + (r1 - r0) * k as f64 / 16.0;
                for l in 0..=16 {
                    let t = t0 + (t1 - t0) * l as f64 / 16.0;
                    let q = s2.d0(r, t);
                    let d = p.distance(&q);
                    if d < best {
                        best = d;
                    }
                }
            }
        }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    use crate::builder::TopoBuilder;
    use crate::shape::TopoShape;
    use crate::tgeometry::GeometryRegistry;
    use occt_core::gp::{GpAx3, GpDir, GpSphere};
    use occt_geom::GeomSphere;

    fn clear_tree(s: &TopoShape) {
        GeometryRegistry::global().clear_shape(s);
        let children = s.tshape.read().unwrap().children.clone();
        for c in children {
            clear_tree(&c);
        }
    }

    fn yz_plane() -> GpPln {
        // Normal +X, X-dir +Y: d0(u,v) = (0, u, v).
        GpPln::new(
            GpAx3::new(
                GpPnt::zero(),
                GpDir::new(1.0, 0.0, 0.0).unwrap(),
                &GpDir::new(0.0, 1.0, 0.0).unwrap(),
            )
            .unwrap(),
        )
    }

    fn z1_plane() -> GpPln {
        // Plane z = 1, parallel to XY.
        GpPln::new(
            GpAx3::new(
                GpPnt::new(0.0, 0.0, 1.0),
                GpDir::new(0.0, 0.0, 1.0).unwrap(),
                &GpDir::new(1.0, 0.0, 0.0).unwrap(),
            )
            .unwrap(),
        )
    }

    #[test]
    fn perpendicular_planes_intersect_in_line() {
        let b = TopoBuilder::new();
        let f1 = b.make_face_plane(&GpPln::new(GpAx3::standard())); // XY
        let f2 = b.make_face_plane(&yz_plane());
        match face_face_intersection(&f1, &f2, 1e-9) {
            FaceIntersect::Line { origin, dir } => {
                let y = GpVec::new(0.0, 1.0, 0.0);
                assert!(dir.cross_magnitude(&y) <= 1e-9, "direction {dir:?}");
                // Origin lies on both planes (XY ⇒ z≈0, YZ ⇒ x≈0).
                assert!(origin.z().abs() < 1e-9);
                assert!(origin.x().abs() < 1e-9);
            }
            other => panic!("expected Line, got {other:?}"),
        }
        clear_tree(&f1.0);
        clear_tree(&f2.0);
    }

    #[test]
    fn parallel_planes_do_not_intersect() {
        let b = TopoBuilder::new();
        let f1 = b.make_face_plane(&GpPln::new(GpAx3::standard()));
        let f2 = b.make_face_plane(&z1_plane());
        assert!(matches!(face_face_intersection(&f1, &f2, 1e-9), FaceIntersect::None));
        clear_tree(&f1.0);
        clear_tree(&f2.0);
    }

    #[test]
    fn plane_and_sphere_intersect_in_circle() {
        let b = TopoBuilder::new();
        let f_plane = b.make_face_plane(&GpPln::new(GpAx3::standard())); // z = 0
        let sphere = GpSphere::new(GpAx3::standard(), 1.0).unwrap(); // unit at origin
        let s: Arc<dyn Surface> = Arc::new(GeomSphere::new(sphere));
        let f_sphere = b.make_face(s, &[]);
        match face_face_intersection(&f_plane, &f_sphere, 0.15) {
            FaceIntersect::Curve { points } => {
                assert!(!points.is_empty(), "expected circle points");
                for p in &points {
                    let d = p.distance(&GpPnt::zero());
                    assert!((d - 1.0).abs() < 0.2, "point {p:?} at radius {d}");
                }
            }
            other => panic!("expected Curve, got {other:?}"),
        }
        clear_tree(&f_plane.0);
        clear_tree(&f_sphere.0);
    }

    #[test]
    fn plane_plane_intersection_matches() {
        let xy = GpPln::new(GpAx3::standard());
        let yz = yz_plane();
        let (origin, dir) = plane_plane_intersection(&xy, &yz).expect("intersecting planes");
        assert!(dir.cross_magnitude(&GpVec::new(0.0, 1.0, 0.0)) <= 1e-9);
        assert!(origin.z().abs() < 1e-9 && origin.x().abs() < 1e-9);
        assert!(plane_plane_intersection(&xy, &z1_plane()).is_none());
    }

    #[test]
    fn point_on_surface_finds_params() {
        let b = TopoBuilder::new();
        let f = b.make_face_plane(&GpPln::new(GpAx3::standard()));
        let s = BRepTool::face_surface(&f).expect("surface");
        let (u, v) = point_on_surface(&*s, &GpPnt::new(1.0, 2.0, 0.0), 0.5).expect("found");
        let q = s.d0(u, v);
        assert!(q.distance(&GpPnt::new(1.0, 2.0, 0.0)) < 0.5, "q {q:?}");
        clear_tree(&f.0);
    }

    #[test]
    fn face_distance_between_parallel_planes() {
        let b = TopoBuilder::new();
        let f1 = b.make_face_plane(&GpPln::new(GpAx3::standard()));
        let f2 = b.make_face_plane(&z1_plane());
        let d = face_distance(&f1, &f2, 1e-9);
        assert!((d - 1.0).abs() < 0.1, "distance {d}");
        clear_tree(&f1.0);
        clear_tree(&f2.0);
    }
}
