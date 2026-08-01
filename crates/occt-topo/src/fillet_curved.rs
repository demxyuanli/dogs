//! Curved-face rolling-ball fillet — a port of `ChFi3d` blending between a
//! curved surface and a plane, or between two curved surfaces.
//!
//! `fillet_edge` only blends two planar faces (the blend is a cylinder). This
//! module extends the rolling-ball idea to edges where at least one adjacent
//! face is a sphere or a cylinder:
//!
//! * **Plane + Sphere** — the rolling ball stays tangent to the plane and to
//!   the sphere, so its centre sweeps the *intersection of the offset plane*
//!   (the plane pushed out by `r`) *and the offset sphere* (radius `R_s + r`).
//!   That intersection is a circle; the envelope of the ball along it is a
//!   **torus** with major radius = circle radius and minor radius = fillet
//!   radius `r`.
//! * **Plane + Cylinder** (axis perpendicular to the plane) — the ball centre
//!   sweeps the intersection of the offset plane and the offset cylinder
//!   (radius `R_c + r`), again a circle, and the blend is a **torus** band
//!   tangent to the cylinder at the inner equator and to the plane at the
//!   outer bottom.
//! * **Sphere + Sphere** — the ball centre sweeps the intersection of the two
//!   offset spheres (radius `R_i + r`), a circle, and the blend is a torus.
//!
//! In every supported case the blend face is a `GeomTorus` band bounded by the
//! two tangency circles (one per adjacent face). The two adjacent faces are
//! rebuilt: the planar face keeps its outer boundary and gets the enlarged
//! plane-tangency circle as a hole; the curved face gets the smaller
//! surface-tangency circle as its new boundary. A cylinder's lateral face is
//! rebuilt with the contact circle, a re-generated seam sub-edge and the
//! original top rim. All new edges are shared between the blend face and the
//! trimmed adjacent face, so the rebuilt shell satisfies `shell_is_closed`.
//!
//! Unsupported pairs (cone, torus, B-spline, a cylinder whose axis is not
//! perpendicular to the plane, ...) return `Err` — documented, never a panic.

use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::gp::{GpAx2, GpAx3, GpDir, GpPnt, GpTorus, GpVec};
use occt_geom::{GeomTorus, Surface};

use crate::brep_surface::{is_planar, sphere_center, surface_normal, SurfaceKind};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::fillet_edge::classify_surface_full;
use crate::shape::{Edge, Face, TopoShape, Wire};
use crate::topo_tools_full::{edges_of, edges_of_wire, faces_of, is_same, wires_of_face};

/// Number of sampled points on a rolling-ball centreline circle.
const CENTERLINE_SAMPLES: usize = 48;

// ===========================================================================
// Analytic surface descriptions
// ===========================================================================

/// A planar surface: a point on it and its unit normal (outward).
#[derive(Debug, Clone, Copy)]
pub struct PlaneInfo {
    pub origin: GpPnt,
    pub normal: GpVec,
}

/// A spherical surface.
#[derive(Debug, Clone, Copy)]
pub struct SphereInfo {
    pub center: GpPnt,
    pub radius: f64,
}

/// A cylindrical surface.
#[derive(Debug, Clone, Copy)]
pub struct CylinderInfo {
    /// A point on the axis.
    pub axis_origin: GpPnt,
    /// Unit axis direction.
    pub axis: GpVec,
    pub radius: f64,
}

/// A circle in 3D (the rolling-ball centreline, or a tangency circle).
#[derive(Debug, Clone, Copy)]
struct CircleGeom {
    center: GpPnt,
    /// Unit circle-plane normal.
    normal: GpVec,
    radius: f64,
}

/// A tangency circle of the blend on one adjacent face.
#[derive(Debug, Clone, Copy)]
struct ContactCircleGeom {
    center: GpPnt,
    /// Unit circle-plane normal.
    normal: GpVec,
    radius: f64,
}

/// Full description of a torus band blend.
#[derive(Debug, Clone)]
struct TorusBlendGeom {
    /// Torus centre (on the centreline circle's axis).
    center: GpPnt,
    /// Unit torus axis.
    axis: GpVec,
    /// Centreline (major) radius.
    major: f64,
    /// Fillet (minor) radius.
    minor: f64,
    contact_a: ContactCircleGeom,
    contact_b: ContactCircleGeom,
}

// ===========================================================================
// Classification
// ===========================================================================

/// The blend topology class of a pair of adjacent faces.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CurvedPair {
    /// A plane and a sphere — torus blend.
    PlaneSphere,
    /// A plane and a cylinder (axis perpendicular to the plane) — torus blend.
    PlaneCylinder,
    /// Two spheres — torus blend.
    SphereSphere,
    /// Anything else — `fillet_edge_curved` returns an error.
    Unsupported,
}

/// Whether a surface is curved (non-planar).
fn is_curved_surface(s: &dyn Surface) -> bool {
    !is_planar(s, 8, 8, 1e-6)
}

/// True when at least one of the two faces is curved (sphere / cylinder / cone
/// / torus / other) — the pair is a candidate for `fillet_edge_curved`.
///
/// A pair of two planar faces is the domain of `fillet_edge` and returns
/// `false`.
pub fn faces_are_curved_compatible(f1: &Face, f2: &Face) -> bool {
    let c1 = BRepTool::face_surface(f1).map(|s| is_curved_surface(s.as_ref())).unwrap_or(false);
    let c2 = BRepTool::face_surface(f2).map(|s| is_curved_surface(s.as_ref())).unwrap_or(false);
    c1 || c2
}

/// Classify the adjacent-face pair into a supported blend topology.
pub fn classify_curved_pair(f1: &Face, f2: &Face) -> CurvedPair {
    let k1 = BRepTool::face_surface(f1)
        .map(|s| classify_surface_full(s.as_ref()))
        .unwrap_or(SurfaceKind::Other);
    let k2 = BRepTool::face_surface(f2)
        .map(|s| classify_surface_full(s.as_ref()))
        .unwrap_or(SurfaceKind::Other);
    match (k1, k2) {
        (SurfaceKind::Plane, SurfaceKind::Sphere) | (SurfaceKind::Sphere, SurfaceKind::Plane) => {
            CurvedPair::PlaneSphere
        }
        (SurfaceKind::Plane, SurfaceKind::Cylinder) | (SurfaceKind::Cylinder, SurfaceKind::Plane) => {
            CurvedPair::PlaneCylinder
        }
        (SurfaceKind::Sphere, SurfaceKind::Sphere) => CurvedPair::SphereSphere,
        _ => CurvedPair::Unsupported,
    }
}

// ===========================================================================
// Surface extraction
// ===========================================================================

/// Finite, sane sampling bounds for a surface (unbounded ranges clamp to ±1).
fn sample_bounds(s: &dyn Surface) -> (f64, f64, f64, f64) {
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    let clamp = |a: f64, b: f64| if a.is_finite() && b.is_finite() && b > a { (a, b) } else { (-1.0, 1.0) };
    let (u0, u1) = clamp(u0, u1);
    let (v0, v1) = clamp(v0, v1);
    (u0, u1, v0, v1)
}

/// Extract plane geometry from a surface, or `None` if it is not planar.
pub fn plane_from_surface(s: &dyn Surface) -> Option<PlaneInfo> {
    if !is_planar(s, 8, 8, 1e-6) {
        return None;
    }
    let (u0, _u1, v0, _v1) = sample_bounds(s);
    let origin = s.d0(u0, v0);
    let normal = surface_normal(s, u0, v0);
    if normal.xyz().square_modulus() < 1e-30 {
        return None;
    }
    Some(PlaneInfo { origin, normal })
}

/// Extract sphere geometry from a surface, or `None`.
pub fn sphere_from_surface(s: &dyn Surface) -> Option<SphereInfo> {
    let center = sphere_center(s)?;
    let (u0, _u1, v0, _v1) = sample_bounds(s);
    let radius = s.d0(u0, v0).distance(&center);
    if radius < 1e-9 {
        return None;
    }
    Some(SphereInfo { center, radius })
}

/// Extract cylinder geometry from a surface, or `None`.
///
/// A cylinder is recognised by: all sampled surface normals perpendicular to
/// one axis direction, and all sampled points at constant distance from that
/// axis line. The axis point is the projection of the first sample onto the
/// axis line through the 2D circumcenter of the projected points.
pub fn cylinder_from_surface(s: &dyn Surface) -> Option<CylinderInfo> {
    let (u0, u1, v0, v1) = sample_bounds(s);
    let nu = 6;
    let nv = 6;
    let hu = ((u1 - u0) * 1e-4).max(1e-7);
    let hv = ((v1 - v0) * 1e-4).max(1e-7);
    let mut pts: Vec<GpPnt> = Vec::new();
    let mut nrm: Vec<GpVec> = Vec::new();
    for i in 0..nu {
        for j in 0..nv {
            let u = u0 + (u1 - u0) * i as f64 / (nu - 1) as f64;
            let v = v0 + (v1 - v0) * j as f64 / (nv - 1) as f64;
            let p0 = s.d0(u, v);
            pts.push(p0);
            let du = GpVec::from_pnts(&p0, &s.d0(u + hu, v));
            let dv = GpVec::from_pnts(&p0, &s.d0(u, v + hv));
            let nn = du.xyz().crossed(&dv.xyz());
            let m = nn.modulus();
            nrm.push(if m > 1e-30 { GpVec::from_xyz(&nn.divided(m)) } else { GpVec::zero() });
        }
    }
    // Axis direction: perpendicular to all normals.
    let mut axis = None;
    'outer: for a in 0..nrm.len() {
        for b in (a + 1)..nrm.len() {
            let c = nrm[a].xyz().crossed(&nrm[b].xyz());
            if c.modulus() > 1e-6 {
                axis = Some(GpVec::from_xyz(&c).normalized());
                break 'outer;
            }
        }
    }
    let axis = axis?;
    for nn in &nrm {
        if nn.dot(&axis).abs() > 1e-3 {
            return None;
        }
    }
    // Project points onto the plane perpendicular to the axis; the axis line
    // passes through the 2D circumcenter of three non-collinear projections.
    let (x2, y2) = project_basis(&axis);
    let proj: Vec<(f64, f64)> = pts
        .iter()
        .map(|q| {
            let v = GpVec::from_pnts(&pts[0], q);
            (v.dot(&x2), v.dot(&y2))
        })
        .collect();
    let (mut cx, mut cy) = (f64::NAN, f64::NAN);
    'outer2: for a in 0..proj.len() {
        for b in (a + 1)..proj.len() {
            for c in (b + 1)..proj.len() {
                let (pa, pb, pc) = (proj[a], proj[b], proj[c]);
                let d = 2.0 * (pa.0 * (pb.1 - pc.1) + pb.0 * (pc.1 - pa.1) + pc.0 * (pa.1 - pb.1));
                if d.abs() < 1e-12 {
                    continue;
                }
                let ux = ((pa.0 * pa.0 + pa.1 * pa.1) * (pb.1 - pc.1)
                    + (pb.0 * pb.0 + pb.1 * pb.1) * (pc.1 - pa.1)
                    + (pc.0 * pc.0 + pc.1 * pc.1) * (pa.1 - pb.1))
                    / d;
                let uy = ((pa.0 * pa.0 + pa.1 * pa.1) * (pc.0 - pb.0)
                    + (pb.0 * pb.0 + pb.1 * pb.1) * (pa.0 - pc.0)
                    + (pc.0 * pc.0 + pc.1 * pc.1) * (pb.0 - pa.0))
                    / d;
                cx = ux;
                cy = uy;
                break 'outer2;
            }
        }
    }
    if !cx.is_finite() {
        return None;
    }
    let axis_pt = pts[0].translated_vec(&x2.multiplied_scalar(cx).added(&y2.multiplied_scalar(cy)));
    let mut radii: Vec<f64> = pts
        .iter()
        .map(|q| {
            let v = GpVec::from_pnts(&axis_pt, q);
            let along = v.dot(&axis);
            v.subtracted(&axis.multiplied_scalar(along)).magnitude()
        })
        .collect();
    radii.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let med = radii[radii.len() / 2];
    if med < 1e-9 {
        return None;
    }
    for r in &radii {
        if (r - med).abs() > 1e-3 * med.abs().max(1.0) {
            return None;
        }
    }
    Some(CylinderInfo { axis_origin: axis_pt, axis, radius: med })
}

/// Two orthonormal vectors spanning the plane perpendicular to `axis`.
fn project_basis(axis: &GpVec) -> (GpVec, GpVec) {
    let ref_v = if axis.x().abs() < 0.9 {
        GpVec::new(1.0, 0.0, 0.0)
    } else {
        GpVec::new(0.0, 1.0, 0.0)
    };
    let x2 = axis.crossed(&ref_v).normalized();
    let y2 = axis.crossed(&x2).normalized();
    (x2, y2)
}

/// A unit vector perpendicular to `z`.
fn perp_dir(z: &GpVec) -> GpVec {
    let ref_v = if z.xyz().z.abs() < 0.9 {
        GpVec::new(0.0, 0.0, 1.0)
    } else {
        GpVec::new(1.0, 0.0, 0.0)
    };
    z.crossed(&ref_v).normalized()
}

// ===========================================================================
// Centreline geometry
// ===========================================================================

/// The rolling-ball centreline for a plane + sphere: the intersection circle
/// of the plane offset by `r` (toward the ball) and the sphere offset by
/// `R_s + r`. The ball centre lies on this circle.
fn plane_sphere_centerline(pln: &PlaneInfo, sph: &SphereInfo, r: f64) -> Result<CircleGeom, String> {
    let n = pln.normal.normalized();
    // Signed distance from the sphere centre to the offset plane n·x = n·O + r,
    // measured along +n (the side the ball rolls on).
    let delta = n.dot(&GpVec::from_pnts(&pln.origin, &sph.center)) - r;
    let rr = sph.radius + r;
    let rho2 = rr * rr - delta * delta;
    if rho2 <= 1e-12 {
        return Err(format!(
            "fillet_curved: ball of radius {r} does not fit between plane and sphere (delta={delta}, R={})",
            sph.radius
        ));
    }
    let rho = rho2.sqrt();
    let center = sph.center.translated_vec(&n.multiplied_scalar(-delta));
    Ok(CircleGeom { center, normal: n, radius: rho })
}

/// The rolling-ball centreline for two spheres: the intersection circle of the
/// two offset spheres (radii `R_i + r`).
fn sphere_sphere_centerline(s1: &SphereInfo, s2: &SphereInfo, r: f64) -> Result<CircleGeom, String> {
    let dv = GpVec::from_pnts(&s1.center, &s2.center);
    let d = dv.magnitude();
    if d < 1e-12 {
        return Err("fillet_curved: coincident spheres".to_string());
    }
    let dhat = dv.normalized();
    let r1 = s1.radius + r;
    let r2 = s2.radius + r;
    if d > r1 + r2 + 1e-9 {
        return Err(format!(
            "fillet_curved: ball of radius {r} cannot bridge two non-overlapping offset spheres"
        ));
    }
    // Distance from s1.center to the radical plane of the two offset spheres.
    let a = (d * d + r1 * r1 - r2 * r2) / (2.0 * d);
    let rho2 = r1 * r1 - a * a;
    if rho2 <= 1e-12 {
        return Err("fillet_curved: offset spheres intersect at a point".to_string());
    }
    let rho = rho2.sqrt();
    let center = s1.center.translated_vec(&dhat.multiplied_scalar(a));
    Ok(CircleGeom { center, normal: dhat, radius: rho })
}

/// Sample a circle at `n` uniformly spaced angles.
fn sample_circle(c: &CircleGeom, n: usize) -> Vec<GpPnt> {
    let e1 = perp_dir(&c.normal);
    let e2 = c.normal.crossed(&e1).normalized();
    (0..n)
        .map(|i| {
            let ang = 2.0 * PI * i as f64 / n as f64;
            let off = e1.multiplied_scalar(c.radius * ang.cos())
                .added(&e2.multiplied_scalar(c.radius * ang.sin()));
            c.center.translated_vec(&off)
        })
        .collect()
}

// ===========================================================================
// Blend geometry per supported pair
// ===========================================================================

/// Build the torus band blend for a plane + sphere pair.
fn plane_sphere_blend(pln: &PlaneInfo, sph: &SphereInfo, r: f64) -> Result<TorusBlendGeom, String> {
    let circle = plane_sphere_centerline(pln, sph, r)?;
    let n = circle.normal;
    let delta = n.dot(&GpVec::from_pnts(&pln.origin, &sph.center)) - r;
    // Tube angle of the sphere-tangency circle: cos v = −rho/(R_s+r), sin v = δ/(R_s+r).
    let v_sphere = delta.atan2(-circle.radius);
    let v_plane = 3.0 * PI / 2.0; // torus bottom: touches the plane
    let v_sphere = if v_sphere < 0.0 { v_sphere + 2.0 * PI } else { v_sphere };
    let plane_contact = ContactCircleGeom {
        center: circle.center.translated_vec(&n.multiplied_scalar(r * (v_plane.sin()))),
        normal: n,
        radius: circle.radius + r * v_plane.cos(),
    };
    let sphere_contact = ContactCircleGeom {
        center: circle.center.translated_vec(&n.multiplied_scalar(r * v_sphere.sin())),
        normal: n,
        radius: circle.radius + r * v_sphere.cos(),
    };
    Ok(TorusBlendGeom {
        center: circle.center,
        axis: n,
        major: circle.radius,
        minor: r,
        contact_a: plane_contact,
        contact_b: sphere_contact,
    })
}

/// Build the torus band blend for a plane + cylinder pair. The cylinder axis
/// must be perpendicular to the plane (a "standing" cylinder); otherwise the
/// offset-plane / offset-cylinder intersection is an ellipse and the exact
/// torus construction does not apply.
fn plane_cylinder_blend(pln: &PlaneInfo, cyl: &CylinderInfo, r: f64, tol: f64) -> Result<TorusBlendGeom, String> {
    let n = pln.normal.normalized();
    let a = cyl.axis.normalized();
    if n.xyz().crossed(&a.xyz()).modulus() > tol.max(1e-6) {
        return Err(
            "fillet_curved: plane+cylinder blend requires the cylinder axis to be perpendicular to the plane"
                .to_string(),
        );
    }
    let rho = cyl.radius + r;
    // Project the axis point onto the offset plane n·x = n·O + r.
    let n_origin = n.dot(&GpVec::from_pnts(&GpPnt::zero(), &pln.origin));
    let n_axis = n.dot(&GpVec::from_pnts(&GpPnt::zero(), &cyl.axis_origin));
    let center = cyl.axis_origin.translated_vec(&n.multiplied_scalar(n_origin + r - n_axis));
    let v_cyl = PI; // inner equator: touches the cylinder
    let v_plane = 3.0 * PI / 2.0; // bottom: touches the plane
    let plane_contact = ContactCircleGeom {
        center: center.translated_vec(&n.multiplied_scalar(r * v_plane.sin())),
        normal: n,
        radius: rho + r * v_plane.cos(),
    };
    let cyl_contact = ContactCircleGeom {
        center: center.translated_vec(&n.multiplied_scalar(r * v_cyl.sin())),
        normal: n,
        radius: rho + r * v_cyl.cos(),
    };
    Ok(TorusBlendGeom {
        center,
        axis: n,
        major: rho,
        minor: r,
        contact_a: plane_contact,
        contact_b: cyl_contact,
    })
}

/// Build the torus band blend for a sphere + sphere pair.
fn sphere_sphere_blend(s1: &SphereInfo, s2: &SphereInfo, r: f64) -> Result<TorusBlendGeom, String> {
    let circle = sphere_sphere_centerline(s1, s2, r)?;
    let dv = GpVec::from_pnts(&s1.center, &s2.center);
    let d = dv.magnitude();
    let dhat = dv.normalized();
    let r1 = s1.radius + r;
    let r2 = s2.radius + r;
    let a = (d * d + r1 * r1 - r2 * r2) / (2.0 * d);
    // Tube angles: cos v_i = −rho / (R_i + r); sin v_1 = −a / (R_1 + r),
    // sin v_2 = (d − a) / (R_2 + r).
    let v1 = (-a).atan2(-circle.radius);
    let v2 = (d - a).atan2(-circle.radius);
    let contact1 = ContactCircleGeom {
        center: circle.center.translated_vec(&dhat.multiplied_scalar(r * v1.sin())),
        normal: dhat,
        radius: circle.radius + r * v1.cos(),
    };
    let contact2 = ContactCircleGeom {
        center: circle.center.translated_vec(&dhat.multiplied_scalar(r * v2.sin())),
        normal: dhat,
        radius: circle.radius + r * v2.cos(),
    };
    Ok(TorusBlendGeom {
        center: circle.center,
        axis: dhat,
        major: circle.radius,
        minor: r,
        contact_a: contact1,
        contact_b: contact2,
    })
}

// ===========================================================================
// Topology rebuild
// ===========================================================================

/// How to rebuild an adjacent face.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AdjacentKind {
    Plane,
    Sphere,
    Cylinder,
}

/// Build a full circular edge on the circle `c`, with its seam at `seam` (a
/// point on the circle, angle 0).
fn build_circle_edge(b: &TopoBuilder, c: &ContactCircleGeom, seam: &GpPnt) -> Result<Edge, String> {
    let nd = GpDir::from_xyz(&c.normal.xyz()).map_err(|e| e.to_string())?;
    let xd = GpDir::from_vec(&GpVec::from_pnts(&c.center, seam))
        .map_err(|_| "fillet_curved: tangency circle seam coincides with its centre".to_string())?;
    let ax2 = GpAx2::new(c.center, nd, xd).map_err(|e| format!("fillet_curved: circle frame: {e}"))?;
    let mut e = b.make_edge_circle(&ax2, c.radius, 0.0, 2.0 * PI);
    let v = b.make_vertex(*seam, 0.0);
    b.add(&mut e.0, &v.0);
    b.add(&mut e.0, &v.0);
    Ok(e)
}

/// A default seam point for a tangency circle (angle 0 along a fixed
/// perpendicular). Used where the seam angle is not constrained by a neighbour.
fn default_seam(c: &ContactCircleGeom) -> GpPnt {
    let e1 = perp_dir(&c.normal);
    c.center.translated_vec(&e1.multiplied_scalar(c.radius))
}

/// The closed (circular) edge of a cylinder lateral face that is not the
/// filleted edge — i.e. the top rim circle.
fn find_top_circle(face: &Face, fillet_edge: &Edge) -> Result<Edge, String> {
    for w in wires_of_face(face) {
        for e in edges_of_wire(&w) {
            if is_same(&e.0, &fillet_edge.0) {
                continue;
            }
            if let Some((a, b)) = BRepTool::edge_vertices(&e) {
                if a.distance(&b) < 1e-9 {
                    return Ok(e.clone());
                }
            }
        }
    }
    Err("fillet_curved: cylinder lateral face has no top rim circle".to_string())
}

/// The point on a cylinder-tangency circle that shares a generator with the
/// cylinder's top rim seam, so the rebuilt lateral face's seam sub-edge stays a
/// straight generator.
fn cylinder_contact_seam(face: &Face, fillet_edge: &Edge, contact: &ContactCircleGeom) -> Result<GpPnt, String> {
    let top = find_top_circle(face, fillet_edge)?;
    let (s_top, _) = BRepTool::edge_vertices(&top).ok_or("fillet_curved: top rim has no vertices")?;
    let v = GpVec::from_pnts(&contact.center, &s_top);
    let along = v.dot(&contact.normal);
    let radial = v.subtracted(&contact.normal.multiplied_scalar(along));
    let m = radial.magnitude();
    if m < 1e-9 {
        return Err("fillet_curved: cylinder top seam lies on the axis".to_string());
    }
    Ok(contact.center.translated_vec(&radial.multiplied_scalar(contact.radius / m)))
}

/// Rebuild a planar adjacent face: keep every original wire except the one
/// containing the filleted edge, and add a new hole wire made of the enlarged
/// plane-tangency circle.
fn rebuild_plane_face(face: &Face, fillet_edge: &Edge, new_circle: &Edge, b: &TopoBuilder) -> Result<Face, String> {
    let surf = BRepTool::face_surface(face).ok_or("fillet_curved: planar face has no surface")?;
    let mut wires: Vec<Wire> = Vec::new();
    for w in wires_of_face(face) {
        let has_fillet = edges_of_wire(&w).iter().any(|e| is_same(&e.0, &fillet_edge.0));
        if has_fillet {
            continue;
        }
        wires.push(w);
    }
    if wires.is_empty() {
        return Err(
            "fillet_curved: the planar face is bounded only by the filleted edge; its trim would be unbounded"
                .to_string(),
        );
    }
    wires.push(b.make_wire(&[new_circle.clone()]));
    Ok(b.make_face(surf, &wires))
}

/// Rebuild a spherical adjacent face: its new boundary is the single
/// surface-tangency circle.
fn rebuild_single_circle_face(face: &Face, new_circle: &Edge, b: &TopoBuilder) -> Result<Face, String> {
    let surf = BRepTool::face_surface(face).ok_or("fillet_curved: face has no surface")?;
    let wire = b.make_wire(&[new_circle.clone()]);
    Ok(b.make_face(surf, &[wire]))
}

/// Rebuild a cylinder lateral face: the base (filleted) circle is replaced by
/// the tangency circle, the seam is re-generated from the tangency circle up to
/// the original top rim, and the top rim circle is reused (shared with the cap).
fn rebuild_cylinder_face(
    face: &Face,
    fillet_edge: &Edge,
    contact_edge: &Edge,
    contact_seam: &GpPnt,
    b: &TopoBuilder,
) -> Result<Face, String> {
    let surf = BRepTool::face_surface(face).ok_or("fillet_curved: cylinder face has no surface")?;
    let top = find_top_circle(face, fillet_edge)?;
    let (s_top, _) = BRepTool::edge_vertices(&top).ok_or("fillet_curved: top rim has no vertices")?;
    let seam_up = b.make_edge_segment(contact_seam, &s_top);
    let wire = b.make_wire(&[contact_edge.clone(), seam_up.clone(), top, seam_up]);
    Ok(b.make_face(surf, &[wire]))
}

/// The `GeomTorus` surface of a torus band blend.
fn torus_surface(geom: &TorusBlendGeom) -> Result<Arc<dyn Surface>, String> {
    let zd = GpDir::from_xyz(&geom.axis.xyz()).map_err(|e| e.to_string())?;
    let xd = GpDir::from_xyz(&perp_dir(&geom.axis).xyz()).map_err(|e| e.to_string())?;
    let ax3 = GpAx3::new(geom.center, zd, &xd).map_err(|e| format!("fillet_curved: torus frame: {e}"))?;
    let t = GpTorus::new(ax3, geom.major, geom.minor).map_err(|e| e.to_string())?;
    Ok(Arc::new(GeomTorus::new(t)))
}

/// Build the blend face: the torus band bounded by the two tangency circles
/// (each shared with one trimmed adjacent face).
fn build_torus_blend_face(geom: &TorusBlendGeom, edge_a: &Edge, edge_b: &Edge, b: &TopoBuilder) -> Result<Face, String> {
    let surf = torus_surface(geom)?;
    let wa = b.make_wire(&[edge_a.clone()]);
    let wb = b.make_wire(&[edge_b.clone()]);
    // The larger tangency circle is the outer boundary.
    let (w_outer, w_inner) = if geom.contact_a.radius >= geom.contact_b.radius {
        (wa, wb)
    } else {
        (wb, wa)
    };
    Ok(b.make_face(surf, &[w_outer, w_inner]))
}

/// Build a trimmed adjacent face, returning the rebuilt face and the (shared)
/// tangency circle edge it uses.
#[allow(clippy::type_complexity)]
fn build_adjacent_face(
    face: &Face,
    kind: AdjacentKind,
    fillet_edge: &Edge,
    contact: &ContactCircleGeom,
    b: &TopoBuilder,
) -> Result<(Face, Edge), String> {
    match kind {
        AdjacentKind::Plane => {
            let e = build_circle_edge(b, contact, &default_seam(contact))?;
            let f = rebuild_plane_face(face, fillet_edge, &e, b)?;
            Ok((f, e))
        }
        AdjacentKind::Sphere => {
            let e = build_circle_edge(b, contact, &default_seam(contact))?;
            let f = rebuild_single_circle_face(face, &e, b)?;
            Ok((f, e))
        }
        AdjacentKind::Cylinder => {
            let seam = cylinder_contact_seam(face, fillet_edge, contact)?;
            let e = build_circle_edge(b, contact, &seam)?;
            let f = rebuild_cylinder_face(face, fillet_edge, &e, &seam, b)?;
            Ok((f, e))
        }
    }
}

// ===========================================================================
// Public API
// ===========================================================================

/// Sample the rolling-ball centre locus for two surfaces. For the supported
/// analytic pairs the locus is a circle (plane+sphere, sphere+sphere,
/// plane+cylinder); the returned points lie on it.
///
/// Each returned centre point is at distance `radius` from the first surface
/// and `radius` from the second surface (measuring against the *offset*
/// surfaces — the plane pushed out by `r`, the sphere/cylinder enlarged by
/// `R + r`). Unsupported pairs return an error.
pub fn rolling_ball_center_line(
    f1: &dyn Surface,
    f2: &dyn Surface,
    radius: f64,
    tol: f64,
) -> Result<Vec<GpPnt>, String> {
    if radius <= 0.0 || !radius.is_finite() {
        return Err("rolling_ball_center_line: radius must be positive".to_string());
    }
    let k1 = classify_surface_full(f1);
    let k2 = classify_surface_full(f2);
    match (k1, k2) {
        (SurfaceKind::Sphere, SurfaceKind::Sphere) => {
            let s1 = sphere_from_surface(f1).ok_or("rolling_ball_center_line: cannot extract sphere 1")?;
            let s2 = sphere_from_surface(f2).ok_or("rolling_ball_center_line: cannot extract sphere 2")?;
            let c = sphere_sphere_centerline(&s1, &s2, radius)?;
            Ok(sample_circle(&c, CENTERLINE_SAMPLES))
        }
        (SurfaceKind::Plane, SurfaceKind::Sphere) => {
            let pln = plane_from_surface(f1).ok_or("rolling_ball_center_line: cannot extract plane")?;
            let sph = sphere_from_surface(f2).ok_or("rolling_ball_center_line: cannot extract sphere")?;
            let c = plane_sphere_centerline(&pln, &sph, radius)?;
            Ok(sample_circle(&c, CENTERLINE_SAMPLES))
        }
        (SurfaceKind::Sphere, SurfaceKind::Plane) => {
            let sph = sphere_from_surface(f1).ok_or("rolling_ball_center_line: cannot extract sphere")?;
            let pln = plane_from_surface(f2).ok_or("rolling_ball_center_line: cannot extract plane")?;
            let c = plane_sphere_centerline(&pln, &sph, radius)?;
            Ok(sample_circle(&c, CENTERLINE_SAMPLES))
        }
        (SurfaceKind::Plane, SurfaceKind::Cylinder) | (SurfaceKind::Cylinder, SurfaceKind::Plane) => {
            let (pln_s, cyl_s) = if k1 == SurfaceKind::Plane { (f1, f2) } else { (f2, f1) };
            let pln = plane_from_surface(pln_s).ok_or("rolling_ball_center_line: cannot extract plane")?;
            let cyl = cylinder_from_surface(cyl_s).ok_or("rolling_ball_center_line: cannot extract cylinder")?;
            let blend = plane_cylinder_blend(&pln, &cyl, radius, tol)?;
            Ok(sample_circle(
                &CircleGeom { center: blend.center, normal: blend.axis, radius: blend.major },
                CENTERLINE_SAMPLES,
            ))
        }
        _ => Err(format!(
            "rolling_ball_center_line: unsupported surface pair ({k1:?}, {k2:?})"
        )),
    }
}

/// Replace a sharp edge of `solid` with a curved-face rolling-ball fillet.
///
/// Unlike `fillet_edge` (which requires both adjacent faces planar), this
/// accepts edges whose adjacent faces are a plane and a sphere, a plane and a
/// standing cylinder, or two spheres. The blend surface is a torus band whose
/// boundary — the two tangency circles — is shared with the trimmed adjacent
/// faces, so the rebuilt shell stays closed.
///
/// `tol` is the geometric tolerance used for the analytic tests (e.g. checking
/// a cylinder axis is perpendicular to the plane).
pub fn fillet_edge_curved(solid: &TopoShape, edge: &Edge, radius: f64, tol: f64) -> Result<TopoShape, String> {
    if radius <= 0.0 || !radius.is_finite() {
        return Err("fillet_edge_curved: radius must be a positive finite value".to_string());
    }
    let adjacent = faces_of(solid)
        .into_iter()
        .filter(|f| {
            wires_of_face(f).iter().any(|w| {
                edges_of_wire(w).iter().any(|e| is_same(&e.0, &edge.0))
            })
        })
        .collect::<Vec<Face>>();
    if adjacent.len() != 2 {
        return Err(format!(
            "fillet_edge_curved: edge is adjacent to {} faces (expected 2)",
            adjacent.len()
        ));
    }
    let f1 = &adjacent[0];
    let f2 = &adjacent[1];
    let s1 = BRepTool::face_surface(f1).ok_or("fillet_edge_curved: face 1 has no surface")?;
    let s2 = BRepTool::face_surface(f2).ok_or("fillet_edge_curved: face 2 has no surface")?;
    let k1 = classify_surface_full(s1.as_ref());
    let k2 = classify_surface_full(s2.as_ref());

    // Compute the blend geometry, assigning contacts to the correct faces.
    let (blend, kind_a, kind_b, contact_a, contact_b): (TorusBlendGeom, AdjacentKind, AdjacentKind, ContactCircleGeom, ContactCircleGeom) =
        match (k1, k2) {
            (SurfaceKind::Plane, SurfaceKind::Sphere) => {
                let pln = plane_from_surface(s1.as_ref()).ok_or("fillet_edge_curved: cannot extract plane")?;
                let sph = sphere_from_surface(s2.as_ref()).ok_or("fillet_edge_curved: cannot extract sphere")?;
                let g = plane_sphere_blend(&pln, &sph, radius)?;
                let (ca, cb) = (g.contact_a, g.contact_b);
                (g, AdjacentKind::Plane, AdjacentKind::Sphere, ca, cb)
            }
            (SurfaceKind::Sphere, SurfaceKind::Plane) => {
                let sph = sphere_from_surface(s1.as_ref()).ok_or("fillet_edge_curved: cannot extract sphere")?;
                let pln = plane_from_surface(s2.as_ref()).ok_or("fillet_edge_curved: cannot extract plane")?;
                let g = plane_sphere_blend(&pln, &sph, radius)?;
                let (ca, cb) = (g.contact_a, g.contact_b);
                (g, AdjacentKind::Sphere, AdjacentKind::Plane, cb, ca)
            }
            (SurfaceKind::Plane, SurfaceKind::Cylinder) => {
                let pln = plane_from_surface(s1.as_ref()).ok_or("fillet_edge_curved: cannot extract plane")?;
                let cyl = cylinder_from_surface(s2.as_ref()).ok_or("fillet_edge_curved: cannot extract cylinder")?;
                let g = plane_cylinder_blend(&pln, &cyl, radius, tol)?;
                let (ca, cb) = (g.contact_a, g.contact_b);
                (g, AdjacentKind::Plane, AdjacentKind::Cylinder, ca, cb)
            }
            (SurfaceKind::Cylinder, SurfaceKind::Plane) => {
                let cyl = cylinder_from_surface(s1.as_ref()).ok_or("fillet_edge_curved: cannot extract cylinder")?;
                let pln = plane_from_surface(s2.as_ref()).ok_or("fillet_edge_curved: cannot extract plane")?;
                let g = plane_cylinder_blend(&pln, &cyl, radius, tol)?;
                let (ca, cb) = (g.contact_a, g.contact_b);
                (g, AdjacentKind::Cylinder, AdjacentKind::Plane, cb, ca)
            }
            (SurfaceKind::Sphere, SurfaceKind::Sphere) => {
                let sp1 = sphere_from_surface(s1.as_ref()).ok_or("fillet_edge_curved: cannot extract sphere 1")?;
                let sp2 = sphere_from_surface(s2.as_ref()).ok_or("fillet_edge_curved: cannot extract sphere 2")?;
                let g = sphere_sphere_blend(&sp1, &sp2, radius)?;
                let (ca, cb) = (g.contact_a, g.contact_b);
                (g, AdjacentKind::Sphere, AdjacentKind::Sphere, ca, cb)
            }
            _ => {
                return Err(format!(
                    "fillet_edge_curved: unsupported face pair ({k1:?}, {k2:?}); \
                     supported pairs are plane+sphere, plane+cylinder and sphere+sphere"
                ));
            }
        };

    let b = TopoBuilder::new();
    let (trimmed1, e1) = build_adjacent_face(f1, kind_a, edge, &contact_a, &b)?;
    let (trimmed2, e2) = build_adjacent_face(f2, kind_b, edge, &contact_b, &b)?;
    let blend_face = build_torus_blend_face(&blend, &e1, &e2, &b)?;

    // Assemble: keep every face except the two replaced adjacent faces.
    let mut faces: Vec<Face> = Vec::new();
    for f in faces_of(solid) {
        if is_same(&f.0, &f1.0) || is_same(&f.0, &f2.0) {
            continue;
        }
        faces.push(f);
    }
    faces.push(trimmed1);
    faces.push(trimmed2);
    faces.push(blend_face);

    let shell = b.make_shell(&faces);
    let solid_out = b.make_solid(&[shell]);
    Ok(solid_out.0)
}

/// Fillet several edges of `solid` sequentially with `fillet_edge_curved`.
///
/// Edges are looked up by their geometric endpoints in the current solid, so a
/// chain whose edges do not interfere (no adjacent face is rebuilt twice)
/// behaves like repeated single-edge fillets. When an earlier fillet consumed
/// an edge, the call errors rather than fillet a wrong edge.
pub fn fillet_edge_curved_chain(
    solid: &TopoShape,
    edge_indices: &[usize],
    radius: f64,
    tol: f64,
) -> Result<TopoShape, String> {
    let original_edges = edges_of(solid);
    let mut current = solid.clone();
    for &i in edge_indices {
        let oe = original_edges
            .get(i)
            .ok_or_else(|| format!("fillet_edge_curved_chain: edge index {i} out of range"))?;
        let (p0, p1) =
            BRepTool::edge_vertices(oe).ok_or("fillet_edge_curved_chain: edge has no curve")?;
        let e = find_edge_by_endpoints(&current, &p0, &p1).ok_or_else(|| {
            format!(
                "fillet_edge_curved_chain: edge {i} was consumed by an earlier fillet and can no \
                 longer be re-found"
            )
        })?;
        current = fillet_edge_curved(&current, &e, radius, tol)?;
    }
    Ok(current)
}

/// Find an edge of `shape` whose endpoints coincide with `p0`/`p1` (either
/// order). For a closed circle both endpoints are the seam point.
fn find_edge_by_endpoints(shape: &TopoShape, p0: &GpPnt, p1: &GpPnt) -> Option<Edge> {
    edges_of(shape).into_iter().find(|e| {
        let (a, b) = BRepTool::edge_vertices(e).unwrap_or((GpPnt::zero(), GpPnt::zero()));
        (a.distance(p0) < 1e-6 && b.distance(p1) < 1e-6)
            || (a.distance(p1) < 1e-6 && b.distance(p0) < 1e-6)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use occt_core::gp::{GpCylinder, GpLin, GpPln, GpSphere};
    use occt_geom::{GeomCylinder, GeomLine, GeomPlane, GeomSphere};
    use crate::shape::{Shell, Solid, Vertex};
    use crate::shell_check::shell_is_closed;
    use crate::tgeometry::GeometryRegistry;
    use crate::topo_tools_full::{faces_of, vertices_of};

    // ------------------------------------------------------------------
    // Test solid builders
    // ------------------------------------------------------------------

    fn clear_tree(s: &TopoShape) {
        GeometryRegistry::global().clear_shape(s);
        let children = s.tshape.read().unwrap().children.clone();
        for c in children {
            clear_tree(&TopoShape::from_handle(c));
        }
    }

    fn shell_of(shape: &TopoShape) -> Shell {
        Shell(TopoShape::from_handle(shape.tshape.read().unwrap().children[0].clone()))
    }

    /// A planar face surface whose axis direction is `normal`.
    fn plane_face(origin: GpPnt, normal: GpDir) -> GpPln {
        let x_dir = if normal.x().abs() > 0.9 {
            GpDir::new(0.0, 1.0, 0.0).unwrap()
        } else if normal.y().abs() > 0.9 {
            GpDir::new(0.0, 0.0, 1.0).unwrap()
        } else {
            GpDir::new(1.0, 0.0, 0.0).unwrap()
        };
        GpPln::new(GpAx3::new(origin, normal, &x_dir).unwrap())
    }

    /// A full circular edge with a fixed seam direction (angle 0). The seam
    /// direction is chosen perpendicular to the circle normal.
    fn build_circle(b: &TopoBuilder, center: GpPnt, normal: GpVec, radius: f64) -> Edge {
        let nd = GpDir::from_xyz(&normal.xyz()).unwrap();
        let xd = if normal.xyz().x.abs() > 0.9 {
            GpDir::new(0.0, 1.0, 0.0).unwrap()
        } else {
            GpDir::new(1.0, 0.0, 0.0).unwrap()
        };
        let ax2 = GpAx2::new(center, nd, xd).unwrap();
        let mut e = b.make_edge_circle(&ax2, radius, 0.0, 2.0 * PI);
        let seam = center.translated_vec(&GpVec::from_xyz(xd.xyz()).multiplied_scalar(radius));
        let v = b.make_vertex(seam, 0.0);
        b.add(&mut e.0, &v.0);
        b.add(&mut e.0, &v.0);
        e
    }

    fn edge_between(edges: &[Edge], i: usize, j: usize) -> Edge {
        let ep: [(usize, usize); 12] = [
            (0, 1), (1, 2), (2, 3), (3, 0),
            (4, 5), (5, 6), (6, 7), (7, 4),
            (0, 4), (1, 5), (2, 6), (3, 7),
        ];
        let idx = ep.iter().position(|&(a, b)| (a == i && b == j) || (a == j && b == i)).unwrap();
        edges[idx].clone()
    }

    fn box_corners(lo: &GpPnt, hi: &GpPnt) -> [GpPnt; 8] {
        let (x0, y0, z0) = (lo.x(), lo.y(), lo.z());
        let (x1, y1, z1) = (hi.x(), hi.y(), hi.z());
        [
            GpPnt::new(x0, y0, z0),
            GpPnt::new(x1, y0, z0),
            GpPnt::new(x1, y1, z0),
            GpPnt::new(x0, y1, z0),
            GpPnt::new(x0, y0, z1),
            GpPnt::new(x1, y0, z1),
            GpPnt::new(x1, y1, z1),
            GpPnt::new(x0, y1, z1),
        ]
    }

    /// A box with `holes.len()` circular holes in its top face. The top face is
    /// a plane at `z = hi.z()` whose outer boundary is the box top square and
    /// whose inner boundary is one circle per hole. Returns the solid and the
    /// hole circle edges.
    fn build_box_with_holes(
        b: &TopoBuilder,
        lo: &GpPnt,
        hi: &GpPnt,
        holes: &[(GpPnt, f64)],
    ) -> (Solid, Vec<Edge>) {
        let c = box_corners(lo, hi);
        let verts: Vec<Vertex> = c.iter().map(|p| b.make_vertex(*p, 0.0)).collect();
        let ep: [(usize, usize); 12] = [
            (0, 1), (1, 2), (2, 3), (3, 0),
            (4, 5), (5, 6), (6, 7), (7, 4),
            (0, 4), (1, 5), (2, 6), (3, 7),
        ];
        let mut edges: Vec<Edge> = Vec::new();
        for &(i, j) in &ep {
            let dir = GpDir::from_vec(&GpVec::from_pnts(&c[i], &c[j])).unwrap();
            let lin = GpLin::from_pnt_dir(c[i], dir);
            let mut e = b.make_edge(Arc::new(GeomLine::new(lin)), 0.0, c[i].distance(&c[j]));
            b.add(&mut e.0, &verts[i].0);
            b.add(&mut e.0, &verts[j].0);
            edges.push(e);
        }

        // Bottom + four side faces.
        let normals = [
            GpDir::new(0.0, 0.0, -1.0).unwrap(),
            GpDir::new(0.0, -1.0, 0.0).unwrap(),
            GpDir::new(0.0, 1.0, 0.0).unwrap(),
            GpDir::new(-1.0, 0.0, 0.0).unwrap(),
            GpDir::new(1.0, 0.0, 0.0).unwrap(),
        ];
        let cycles: [[usize; 4]; 5] = [
            [0, 1, 2, 3], // bottom
            [0, 1, 5, 4], // -Y
            [3, 2, 6, 7], // +Y
            [0, 3, 7, 4], // -X
            [1, 2, 6, 5], // +X
        ];
        let origins = [c[0], c[0], c[3], c[0], c[1]];
        let mut faces: Vec<Face> = Vec::new();
        for i in 0..5 {
            let pln = plane_face(origins[i], normals[i]);
            let w = b.make_wire(&[
                edge_between(&edges, cycles[i][0], cycles[i][1]),
                edge_between(&edges, cycles[i][1], cycles[i][2]),
                edge_between(&edges, cycles[i][2], cycles[i][3]),
                edge_between(&edges, cycles[i][3], cycles[i][0]),
            ]);
            faces.push(b.make_face(Arc::new(GeomPlane::new(pln)), &[w]));
        }

        // Top annulus: outer square + one hole circle per hole.
        let top_w = b.make_wire(&[
            edge_between(&edges, 4, 5),
            edge_between(&edges, 5, 6),
            edge_between(&edges, 6, 7),
            edge_between(&edges, 7, 4),
        ]);
        let mut hole_edges = Vec::new();
        let mut hole_wires = Vec::new();
        for (hc, hr) in holes {
            let e = build_circle(b, *hc, GpVec::new(0.0, 0.0, 1.0), *hr);
            hole_edges.push(e.clone());
            hole_wires.push(b.make_wire(&[e]));
        }
        let mut top_wires = vec![top_w];
        top_wires.extend(hole_wires);
        let top_face = b.make_face(
            Arc::new(GeomPlane::new(plane_face(c[4], GpDir::new(0.0, 0.0, 1.0).unwrap()))),
            &top_wires,
        );
        faces.push(top_face);

        let shell = b.make_shell(&faces);
        let solid = b.make_solid(&[shell]);
        (solid, hole_edges)
    }

    /// A sphere cap sitting on the top face of a box: box `[-1.5,1.5]²×[-1,0]`
    /// with a sphere of radius `r_s` centred at the origin (on the top plane).
    fn build_sphere_on_box(r_s: f64) -> (Solid, Edge) {
        let b = TopoBuilder::new();
        let (solid, holes) = build_box_with_holes(
            &b,
            &GpPnt::new(-1.5, -1.5, -1.0),
            &GpPnt::new(1.5, 1.5, 0.0),
            &[(GpPnt::zero(), r_s)],
        );
        let hole = holes[0].clone();
        let mut sph = GpSphere::new(GpAx3::standard(), r_s).unwrap();
        sph.set_location(GpPnt::zero());
        let cap = b.make_face(Arc::new(GeomSphere::new(sph)), &[b.make_wire(&[hole.clone()])]);
        let mut fs = faces_of(&solid.0);
        fs.push(cap);
        let shell = b.make_shell(&fs);
        (b.make_solid(&[shell]), hole)
    }

    /// A cylinder standing on the top face of a box: box `[-1.5,1.5]²×[-1,0]`
    /// with a cylinder of radius `r_c`, axis Z through the origin, height `h`.
    fn build_cylinder_on_box(r_c: f64, h: f64) -> (Solid, Edge) {
        let b = TopoBuilder::new();
        let (solid, holes) = build_box_with_holes(
            &b,
            &GpPnt::new(-1.5, -1.5, -1.0),
            &GpPnt::new(1.5, 1.5, 0.0),
            &[(GpPnt::zero(), r_c)],
        );
        let base = holes[0].clone();
        let ax = GpAx3::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
        let top = build_circle(&b, GpPnt::new(0.0, 0.0, h), GpVec::new(0.0, 0.0, 1.0), r_c);
        let seam = b.make_edge_segment(&GpPnt::new(r_c, 0.0, 0.0), &GpPnt::new(r_c, 0.0, h));
        let lateral_wire = b.make_wire(&[base.clone(), seam.clone(), top.clone(), seam]);
        let lateral = b.make_face(Arc::new(GeomCylinder::new(GpCylinder::new(ax, r_c).unwrap())), &[lateral_wire]);
        let top_cap = b.make_face(
            Arc::new(GeomPlane::new(plane_face(GpPnt::new(0.0, 0.0, h), GpDir::new(0.0, 0.0, 1.0).unwrap()))),
            &[b.make_wire(&[top])],
        );
        let mut fs = faces_of(&solid.0);
        fs.push(lateral);
        fs.push(top_cap);
        let shell = b.make_shell(&fs);
        (b.make_solid(&[shell]), base)
    }

    /// Two radius-1 spheres at distance `d` along +X. The union solid is two
    /// caps sharing the intersection circle.
    fn build_two_spheres(d: f64) -> (Solid, Edge) {
        let b = TopoBuilder::new();
        let r = 1.0;
        let a = (r * r - (d / 2.0) * (d / 2.0)).sqrt();
        let ic = build_circle(&b, GpPnt::new(d / 2.0, 0.0, 0.0), GpVec::new(1.0, 0.0, 0.0), a);
        let ic_wire = b.make_wire(&[ic.clone()]);
        let mut s1 = GpSphere::new(GpAx3::standard(), r).unwrap();
        s1.set_location(GpPnt::zero());
        let mut s2 = GpSphere::new(GpAx3::standard(), r).unwrap();
        s2.set_location(GpPnt::new(d, 0.0, 0.0));
        let cap1 = b.make_face(Arc::new(GeomSphere::new(s1)), &[ic_wire.clone()]);
        let cap2 = b.make_face(Arc::new(GeomSphere::new(s2)), &[ic_wire]);
        let shell = b.make_shell(&[cap1, cap2]);
        (b.make_solid(&[shell]), ic)
    }

    /// Two sphere caps on a box, for the chain test: box `[-2,2]²×[-1,0]`, two
    /// spheres of radius 0.7 centred at `(±0.9, 0, 0)`.
    fn build_two_spheres_on_box() -> (Solid, Vec<Edge>) {
        let b = TopoBuilder::new();
        let (solid, holes) = build_box_with_holes(
            &b,
            &GpPnt::new(-2.0, -2.0, -1.0),
            &GpPnt::new(2.0, 2.0, 0.0),
            &[(GpPnt::new(-0.9, 0.0, 0.0), 0.7), (GpPnt::new(0.9, 0.0, 0.0), 0.7)],
        );
        let mut sph1 = GpSphere::new(GpAx3::standard(), 0.7).unwrap();
        sph1.set_location(GpPnt::new(-0.9, 0.0, 0.0));
        let mut sph2 = GpSphere::new(GpAx3::standard(), 0.7).unwrap();
        sph2.set_location(GpPnt::new(0.9, 0.0, 0.0));
        let cap1 = b.make_face(Arc::new(GeomSphere::new(sph1)), &[b.make_wire(&[holes[0].clone()])]);
        let cap2 = b.make_face(Arc::new(GeomSphere::new(sph2)), &[b.make_wire(&[holes[1].clone()])]);
        let mut fs = faces_of(&solid.0);
        fs.push(cap1);
        fs.push(cap2);
        let shell = b.make_shell(&fs);
        (b.make_solid(&[shell]), holes)
    }

    /// Axisymmetric solid for the volume test: a base cylinder of radius `r_b`,
    /// height 1 below z=0, topped by a sphere cap of radius `r_s` at the origin.
    fn build_sphere_on_cylinder(r_b: f64, r_s: f64) -> (Solid, Edge) {
        let b = TopoBuilder::new();
        let ax = GpAx3::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
        let bottom_circle = build_circle(&b, GpPnt::new(0.0, 0.0, -1.0), GpVec::new(0.0, 0.0, 1.0), r_b);
        let outer_circle = build_circle(&b, GpPnt::zero(), GpVec::new(0.0, 0.0, 1.0), r_b);
        let inner_circle = build_circle(&b, GpPnt::zero(), GpVec::new(0.0, 0.0, 1.0), r_s);
        let seam = b.make_edge_segment(&GpPnt::new(r_b, 0.0, -1.0), &GpPnt::new(r_b, 0.0, 0.0));
        let lateral_wire = b.make_wire(&[bottom_circle.clone(), seam.clone(), outer_circle.clone(), seam]);
        let lateral = b.make_face(Arc::new(GeomCylinder::new(GpCylinder::new(ax, r_b).unwrap())), &[lateral_wire]);
        let bottom_cap = b.make_face(
            Arc::new(GeomPlane::new(plane_face(GpPnt::new(0.0, 0.0, -1.0), GpDir::new(0.0, 0.0, -1.0).unwrap()))),
            &[b.make_wire(&[bottom_circle])],
        );
        let top_annulus = b.make_face(
            Arc::new(GeomPlane::new(plane_face(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap()))),
            &[b.make_wire(&[outer_circle]), b.make_wire(&[inner_circle.clone()])],
        );
        let mut sph = GpSphere::new(GpAx3::standard(), r_s).unwrap();
        sph.set_location(GpPnt::zero());
        let cap = b.make_face(Arc::new(GeomSphere::new(sph)), &[b.make_wire(&[inner_circle.clone()])]);
        let shell = b.make_shell(&[bottom_cap, lateral, top_annulus, cap]);
        (b.make_solid(&[shell]), inner_circle)
    }

    /// The first face that is not plane/sphere/cylinder — the torus blend face,
    /// since every input face of the test solids is plane/sphere/cylinder.
    fn find_blend_face<'a>(faces: &'a [Face]) -> &'a Face {
        faces
            .iter()
            .find(|f| {
                BRepTool::face_surface(f)
                    .map(|s| classify_surface_full(s.as_ref()) == SurfaceKind::Other)
                    .unwrap_or(false)
            })
            .expect("torus blend face")
    }

    /// Distance from `p` to the circle (centre `c`, unit normal `n`, radius
    /// `rho`) in the plane perpendicular to `n`.
    fn distance_to_circle(p: &GpPnt, c: &GpPnt, n: &GpVec, rho: f64) -> f64 {
        let v = GpVec::from_pnts(c, p);
        let along = v.dot(n);
        let radial = v.subtracted(&n.multiplied_scalar(along));
        (radial.magnitude() - rho).hypot(along)
    }

    // ------------------------------------------------------------------
    // Tests
    // ------------------------------------------------------------------

    #[test]
    fn plane_sphere_blend_closed() {
        let (solid, edge) = build_sphere_on_box(1.0);
        let out = fillet_edge_curved(&solid.0, &edge, 0.2, 1e-6).expect("plane+sphere fillet");
        let faces = faces_of(&out);
        assert_eq!(faces.len(), 8, "7 input faces + 1 blend face, got {}", faces.len());
        assert!(shell_is_closed(&shell_of(&out)), "plane+sphere fillet is closed");
        // The blend face is a torus band: every sampled point is at distance r
        // from the rolling-ball centreline circle.
        let blend = find_blend_face(&faces);
        let s = BRepTool::face_surface(blend).unwrap();
        let r = 0.2;
        let rho = (1.0_f64 + 2.0 * r).sqrt();
        let center = GpPnt::new(0.0, 0.0, r);
        let axis = GpVec::new(0.0, 0.0, 1.0);
        for (u, v) in [(0.3, 3.9), (1.0, 4.2), (2.0, 4.5), (4.0, 3.3), (5.0, 4.6)] {
            let p = s.d0(u, v);
            let d = distance_to_circle(&p, &center, &axis, rho);
            assert!((d - r).abs() < 1e-6, "blend point {p:?} at torus distance {d} (expected {r})");
        }
        clear_tree(&out);
        clear_tree(&solid.0);
    }

    #[test]
    fn plane_cylinder_blend() {
        let (solid, edge) = build_cylinder_on_box(1.0, 1.0);
        let out = fillet_edge_curved(&solid.0, &edge, 0.2, 1e-6).expect("plane+cylinder fillet");
        let faces = faces_of(&out);
        assert_eq!(faces.len(), 9, "8 input faces + 1 blend face, got {}", faces.len());
        assert!(shell_is_closed(&shell_of(&out)), "plane+cylinder fillet is closed");
        // The blend face is a torus band around the centreline circle of radius
        // R_c + r at height r above the plane.
        let blend = find_blend_face(&faces);
        let s = BRepTool::face_surface(blend).unwrap();
        let r = 0.2;
        let rho = 1.0 + r;
        let center = GpPnt::new(0.0, 0.0, r);
        let axis = GpVec::new(0.0, 0.0, 1.0);
        for (u, v) in [(0.3, 3.9), (1.0, 4.2), (2.0, 4.5), (4.0, 3.3), (5.0, 4.6)] {
            let p = s.d0(u, v);
            let d = distance_to_circle(&p, &center, &axis, rho);
            assert!((d - r).abs() < 1e-6, "blend point {p:?} at torus distance {d} (expected {r})");
        }
        clear_tree(&out);
        clear_tree(&solid.0);
    }

    #[test]
    fn blend_radius_tangent() {
        let (solid, edge) = build_cylinder_on_box(1.0, 1.0);
        let out = fillet_edge_curved(&solid.0, &edge, 0.2, 1e-6).unwrap();
        let faces = faces_of(&out);
        let blend = find_blend_face(&faces);
        let s = BRepTool::face_surface(blend).unwrap();
        let (u0, _u1) = s.u_range();
        // Plane contact (v = 3π/2): the blend sits on z=0 and its normal is
        // parallel to the plane normal.
        let p_plane = s.d0(u0, 3.0 * PI / 2.0);
        assert!(p_plane.z().abs() < 1e-9, "plane contact must lie on z=0, got {p_plane:?}");
        let n_plane = surface_normal(s.as_ref(), u0, 3.0 * PI / 2.0);
        let z_axis = GpVec::new(0.0, 0.0, 1.0);
        assert!(
            n_plane.xyz().crossed(&z_axis.xyz()).modulus() < 1e-3,
            "blend normal {:?} must be parallel to the plane normal",
            n_plane
        );
        // Cylinder contact (v = π): the blend sits on the radius-R_c cylinder
        // and its normal is radial (perpendicular to the axis).
        let p_cyl = s.d0(u0, PI);
        let radial_len = GpVec::new(p_cyl.x(), p_cyl.y(), 0.0).magnitude();
        assert!((radial_len - 1.0).abs() < 1e-6, "cylinder contact radius {radial_len}");
        let n_cyl = surface_normal(s.as_ref(), u0, PI);
        assert!(n_cyl.xyz().z.abs() < 1e-3, "cylinder-contact normal must be radial");
        let radial = GpVec::new(p_cyl.x(), p_cyl.y(), 0.0).normalized();
        assert!(
            n_cyl.xyz().crossed(&radial.xyz()).modulus() < 1e-3,
            "blend normal {:?} must be parallel to the cylinder radial",
            n_cyl
        );
        clear_tree(&out);
        clear_tree(&solid.0);
    }

    #[test]
    fn curved_unsupported_errors() {
        // Two torus faces sharing a circle edge: both classify as Other, so the
        // pair is unsupported and fillet_edge_curved must error (documented).
        let b = TopoBuilder::new();
        let circle = build_circle(&b, GpPnt::zero(), GpVec::new(0.0, 0.0, 1.0), 1.0);
        let w = b.make_wire(&[circle.clone()]);
        let f1 = b.make_face(Arc::new(GeomTorus::new(GpTorus::new(GpAx3::standard(), 3.0, 0.5).unwrap())), &[w.clone()]);
        let f2 = b.make_face(Arc::new(GeomTorus::new(GpTorus::new(GpAx3::standard(), 3.0, 0.5).unwrap())), &[w]);
        let shell = b.make_shell(&[f1.clone(), f2.clone()]);
        let solid = b.make_solid(&[shell]);
        let res = fillet_edge_curved(&solid.0, &circle, 0.2, 1e-6);
        assert!(res.is_err(), "torus+torus must be unsupported");
        // But faces_are_curved_compatible still reports the pair as curved.
        assert!(faces_are_curved_compatible(&f1, &f2));
        assert_eq!(classify_curved_pair(&f1, &f2), CurvedPair::Unsupported);
        clear_tree(&solid.0);
        clear_tree(&f1.0);
        clear_tree(&f2.0);
    }

    #[test]
    fn curved_chain_two_edges() {
        let (solid, holes) = build_two_spheres_on_box();
        let es = edges_of(&solid.0);
        let i0 = es
            .iter()
            .position(|e| {
                let (a, _) = BRepTool::edge_vertices(e).unwrap();
                a.distance(&GpPnt::new(-0.2, 0.0, 0.0)) < 1e-6
            })
            .expect("sphere-1 base circle");
        let i1 = es
            .iter()
            .position(|e| {
                let (a, _) = BRepTool::edge_vertices(e).unwrap();
                a.distance(&GpPnt::new(1.6, 0.0, 0.0)) < 1e-6
            })
            .expect("sphere-2 base circle");
        assert_eq!(holes.len(), 2);
        let out = fillet_edge_curved_chain(&solid.0, &[i0, i1], 0.2, 1e-6).expect("two-edge chain");
        assert_eq!(faces_of(&out).len(), 10, "8 input faces + 2 blends");
        assert!(shell_is_closed(&shell_of(&out)), "two-edge curved chain is closed");
        clear_tree(&out);
        clear_tree(&solid.0);
    }

    #[test]
    fn sphere_sphere_blend() {
        let (solid, edge) = build_two_spheres(1.2);
        let out = fillet_edge_curved(&solid.0, &edge, 0.2, 1e-6).expect("sphere+sphere fillet");
        let faces = faces_of(&out);
        assert_eq!(faces.len(), 3, "2 input faces + 1 blend face, got {}", faces.len());
        assert!(shell_is_closed(&shell_of(&out)), "sphere+sphere fillet is closed");
        let blend = find_blend_face(&faces);
        let s = BRepTool::face_surface(blend).unwrap();
        let r = 0.2;
        let rho = (1.2_f64 * 1.2 - 0.6 * 0.6).sqrt();
        let center = GpPnt::new(0.6, 0.0, 0.0);
        let axis = GpVec::new(1.0, 0.0, 0.0);
        for (u, v) in [(0.3, 3.8), (1.0, 4.2), (2.0, 4.7), (4.0, 5.3), (5.0, 5.6)] {
            let p = s.d0(u, v);
            let d = distance_to_circle(&p, &center, &axis, rho);
            assert!((d - r).abs() < 1e-6, "blend point {p:?} at torus distance {d} (expected {r})");
        }
        clear_tree(&out);
        clear_tree(&solid.0);
    }

    #[test]
    fn rolling_center_line_sphere_plane() {
        let pln = GeomPlane::new(GpPln::new(GpAx3::standard()));
        let mut sph = GpSphere::new(GpAx3::standard(), 1.0).unwrap();
        sph.set_location(GpPnt::new(0.0, 0.0, 0.5));
        let sph_geom = GeomSphere::new(sph);
        let r = 0.2;
        let pts = rolling_ball_center_line(&pln, &sph_geom, r, 1e-6).expect("rolling centre line");
        assert_eq!(pts.len(), CENTERLINE_SAMPLES);
        for p in &pts {
            // The centre lies on the plane offset by r (here z = r).
            assert!((p.z() - r).abs() < 1e-9, "centre height {p:?}");
            // And on the sphere offset by R + r.
            let d = p.distance(&GpPnt::new(0.0, 0.0, 0.5));
            assert!((d - 1.2).abs() < 1e-9, "centre distance {d} (expected 1.2)");
        }
    }

    #[test]
    fn fillet_curved_volume_conserved() {
        let r_b = 1.5;
        let r_s = 1.0;
        let r = 0.2;
        let (solid, edge) = build_sphere_on_cylinder(r_b, r_s);
        let out = fillet_edge_curved(&solid.0, &edge, r, 1e-6).expect("fillet");
        assert!(shell_is_closed(&shell_of(&out)), "filleted axisymmetric solid is closed");

        // Pre-fillet volume: base cylinder (radius r_b, height 1) + upper
        // hemisphere of radius r_s.
        let v_pre = PI * r_b * r_b * 1.0 + (2.0 / 3.0) * PI * r_s.powi(3);

        // Post-fillet volume of revolution of the meridian profile:
        //   z < 0      : radius r_b  (base cylinder)
        //   0 ≤ z ≤ z_c: blend band  r(z) = rho − sqrt(r² − (z − r)²)
        //   z_c ≤ z ≤ r_s: sphere     r(z) = sqrt(r_s² − z²)
        let rho = (r_s * r_s + 2.0 * r_s * r).sqrt();
        let z_c = r * r_s / (r_s + r);
        let sphere_part = PI * ((r_s - z_c) - (r_s.powi(3) - z_c.powi(3)) / 3.0);
        let n = 4000;
        let h = z_c / n as f64;
        let band = |z: f64| {
            let rad = rho - (r * r - (z - r) * (z - r)).sqrt();
            rad * rad
        };
        let mut band_part = 0.0;
        for k in 0..n {
            let z0 = k as f64 * h;
            let z1 = (k + 1) as f64 * h;
            band_part += PI * h * 0.5 * (band(z0) + band(z1));
        }
        let v_post = PI * r_b * r_b * 1.0 + sphere_part + band_part;

        let rel = (v_post - v_pre).abs() / v_pre;
        assert!(rel < 0.15, "fillet changes volume by {:.2}% (pre {v_pre}, post {v_post})", rel * 100.0);
        clear_tree(&out);
        clear_tree(&solid.0);
    }

    #[test]
    fn faces_are_curved_compatible_mixed() {
        let (solid, edge) = build_sphere_on_box(1.0);
        let adjacent = faces_of(&solid.0)
            .into_iter()
            .filter(|f| {
                wires_of_face(f).iter().any(|w| {
                    edges_of_wire(w).iter().any(|e| is_same(&e.0, &edge.0))
                })
            })
            .collect::<Vec<Face>>();
        assert_eq!(adjacent.len(), 2);
        // One face is the plane, the other the sphere cap.
        assert!(faces_are_curved_compatible(&adjacent[0], &adjacent[1]));
        assert_eq!(classify_curved_pair(&adjacent[0], &adjacent[1]), CurvedPair::PlaneSphere);
        // Two planar box faces are NOT curved-compatible.
        let (boxy, _) = build_sphere_on_box(1.0);
        let all = faces_of(&boxy.0);
        let planars: Vec<&Face> = all.iter().filter(|f| is_planar(
            BRepTool::face_surface(f).unwrap().as_ref(), 8, 8, 1e-6,
        )).collect();
        assert!(planars.len() >= 2, "box has at least two planar faces");
        assert!(!faces_are_curved_compatible(planars[0], planars[1]));
        let _ = vertices_of(&solid.0);
        clear_tree(&solid.0);
        clear_tree(&boxy.0);
    }
}
