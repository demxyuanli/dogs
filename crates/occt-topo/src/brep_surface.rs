//! Surface classification and p-curve helpers for faces.
//! Source: `GeomAdaptor_Surface.hxx`, `BRep_Tool::Surface`,
//! `GeomAPI_ProjectPointOnSurf.hxx`.
//!
//! Trait objects (`Arc<dyn Surface>`) cannot be downcast, so classification is
//! done by geometric invariants (constant normal ⇒ plane, constant radius from
//! a solved center ⇒ sphere). This mirrors what OCCT's `GeomAdaptor` type tag
//! provides, at a slightly higher cost.

use occt_core::gp::{GpAx3, GpDir, GpPln, GpPnt, GpPnt2d, GpVec};
use occt_geom::Surface;


use crate::shape::{Edge, Face};
use crate::tgeometry::GeometryRegistry;

/// Coarse classification of a surface's analytic type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SurfaceKind {
    Plane,
    Sphere,
    Cylinder,
    Cone,
    Torus,
    Other,
}

/// Unit normal of a surface at (u, v), robust to surfaces whose `d1` returns
/// zero vectors (some ported surfaces only implement `d0`).
pub fn surface_normal(s: &dyn Surface, u: f64, v: f64) -> GpVec {
    let (_, du, dv) = s.d1(u, v);
    let n = du.xyz().crossed(dv.xyz());
    if n.square_modulus() > 1e-30 {
        let m = n.modulus();
        return GpVec::new(n.x / m, n.y / m, n.z / m);
    }
    // Fallback: central finite differences of d0.
    let eps = 1e-6;
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    let hu = if u1 > u0 { (u1 - u0) * 1e-4 } else { eps };
    let hv = if v1 > v0 { (v1 - v0) * 1e-4 } else { eps };
    let p0 = s.d0(u, v);
    let duv = GpVec::from_pnts(&p0, &s.d0(u + hu, v));
    let dvv = GpVec::from_pnts(&p0, &s.d0(u, v + hv));
    let n = duv.xyz().crossed(dvv.xyz());
    let m = n.modulus();
    if m > 1e-30 {
        GpVec::new(n.x / m, n.y / m, n.z / m)
    } else {
        GpVec::zero()
    }
}

/// Whether all sampled surface normals are parallel (within `tol`).
pub fn is_planar(s: &dyn Surface, nu: usize, nv: usize, tol: f64) -> bool {
    let (u0, u1, v0, v1) = sample_bounds(s);
    let (mut first, mut first_ok) = (GpVec::zero(), false);
    for i in 0..nu {
        for j in 0..nv {
            let u = u0 + (u1 - u0) * i as f64 / (nu.max(1) - 1) as f64;
            let v = v0 + (v1 - v0) * j as f64 / (nv.max(1) - 1) as f64;
            let n = surface_normal(s, u, v);
            if n.xyz().square_modulus() < 1e-30 {
                continue;
            }
            if !first_ok {
                first = n;
                first_ok = true;
            } else if n.xyz().crossed(first.xyz()).modulus() > tol {
                return false;
            }
        }
    }
    first_ok
}

/// Classify a surface by its **exact analytic type**.
///
/// Mirrors `GeomAdaptor_Surface::Load` (`GeomAdaptor_Surface.cxx:422-513`),
/// which compares the surface's exact class in a fixed order
/// (rectangular-trimmed → recurse on the basis; Plane, Cylinder, Cone, Sphere,
/// Torus, …; anything else is `GeomAbs_OtherSurface`). The port exposes the same
/// information through the `Surface` type queries, so no sampling is needed.
///
/// The previous body sampled an 8×8 grid and compared normals / centre
/// distances with `1e-6` / `1e-4` thresholds, and could only ever report
/// `Plane`, `Sphere` or `Other` (audit A30).
pub fn classify_surface(s: &dyn Surface) -> SurfaceKind {
    if let Some(basis) = s.rectangular_trimmed_basis() {
        return classify_surface(basis.as_ref());
    }
    if s.gp_pln().is_some() {
        return SurfaceKind::Plane;
    }
    if s.gp_cylinder().is_some() {
        return SurfaceKind::Cylinder;
    }
    if s.gp_cone().is_some() {
        return SurfaceKind::Cone;
    }
    if s.gp_sphere().is_some() {
        return SurfaceKind::Sphere;
    }
    if s.gp_torus().is_some() {
        return SurfaceKind::Torus;
    }
    SurfaceKind::Other
}

/// Finite, sane sampling bounds for a surface (unbounded ranges clamp to ±1).
fn sample_bounds(s: &dyn Surface) -> (f64, f64, f64, f64) {
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    let clamp = |a: f64, b: f64| if a.is_finite() && b.is_finite() && b > a { (a, b) } else { (-1.0, 1.0) };
    let (u0, u1) = clamp(u0, u1);
    let (v0, v1) = clamp(v0, v1);
    (u0, u1, v0, v1)
}

/// 3×3 determinant.
fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// Solve a 3×3 linear system via Cramer's rule.
fn solve3(a: &[[f64; 3]; 3], rhs: &[f64; 3]) -> Option<[f64; 3]> {
    let d = det3(a);
    if d.abs() < 1e-20 {
        return None;
    }
    let mut x = [0.0; 3];
    for k in 0..3 {
        let mut m = *a;
        for i in 0..3 {
            m[i][k] = rhs[i];
        }
        x[k] = det3(&m) / d;
    }
    Some(x)
}

/// Circumcenter of three non-collinear 3D points, if it exists.
///
/// Solves the perpendicular-bisector system: the center O is equidistant from
/// a, b, c (two bisector equations) and lies in their plane. Cramer's rule on
/// the resulting 3×3 system is robust for near-equal points.
fn circumcenter(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> Option<GpPnt> {
    let d1 = GpVec::from_pnts(a, b);
    let d2 = GpVec::from_pnts(a, c);
    let n = d1.xyz().crossed(d2.xyz());
    if n.modulus() < 1e-20 {
        return None;
    }
    let n2 = |p: &GpPnt| p.coord.dot(&p.coord);
    // O·d1 = (|b|² − |a|²)/2 ; O·d2 = (|c|² − |a|²)/2 ; O·n = a·n.
    let mat = [
        [d1.xyz().x, d1.xyz().y, d1.xyz().z],
        [d2.xyz().x, d2.xyz().y, d2.xyz().z],
        [n.x, n.y, n.z],
    ];
    let rhs = [0.5 * (n2(b) - n2(a)), 0.5 * (n2(c) - n2(a)), a.coord.dot(&n)];
    let o = solve3(&mat, &rhs)?;
    Some(GpPnt::new(o[0], o[1], o[2]))
}

/// If the surface is spherical, return its center.
///
/// The three sample points are chosen at the mid-latitude band (away from the
/// poles, where the longitude parameter is degenerate and points collapse),
/// with distinct u values so they are non-collinear.
pub fn sphere_center(s: &dyn Surface) -> Option<GpPnt> {
    let (u0, u1, v0, v1) = sample_bounds(s);
    let u_mid = 0.5 * (u0 + u1);
    let v_mid = 0.5 * (v0 + v1);
    let span = (v1 - v0).abs().max(1e-6);
    let a = s.d0(u0, v_mid);
    let b = s.d0(u_mid, v_mid);
    let c = s.d0(u0, v_mid + 0.25 * span);
    circumcenter(&a, &b, &c)
}

/// The plane underlying a planar face, if it is planar.
pub fn face_plane(face: &Face) -> Option<GpPln> {
    let surf = GeometryRegistry::global().face_surface(&face.0)?;
    if !is_planar(surf.as_ref(), 8, 8, 1e-6) {
        return None;
    }
    let (u0, _u1, v0, _v1) = sample_bounds(surf.as_ref());
    let o = surf.d0(u0, v0);
    let n = surface_normal(surf.as_ref(), u0, v0);
    let d = GpDir::from_vec(&n).ok()?;
    // Pick an in-plane X direction perpendicular to the normal (a fixed axis
    // may be parallel to the normal, which GpAx3 rejects).
    let z_axis = occt_core::gp::GpDir::new(0.0, 0.0, 1.0).unwrap();
    let x_dir = if d.is_normal(&z_axis) {
        z_axis
    } else {
        occt_core::gp::GpDir::new(1.0, 0.0, 0.0).unwrap()
    };
    Some(GpPln::new(GpAx3::new(o, d, &x_dir).ok()?))
}

/// Whether a face is planar.
pub fn face_is_planar(face: &Face) -> bool {
    match GeometryRegistry::global().face_surface(&face.0) {
        Some(s) => is_planar(s.as_ref(), 8, 8, 1e-6),
        None => false,
    }
}

/// Surface normal at a face's natural (u, v) parameter point.
pub fn face_normal(face: &Face, u: f64, v: f64) -> Option<GpVec> {
    let surf = GeometryRegistry::global().face_surface(&face.0)?;
    Some(surface_normal(surf.as_ref(), u, v))
}

/// Centroid of a face surface sampled on a `nu × nv` grid.
pub fn face_centroid(face: &Face, nu: usize, nv: usize) -> Option<GpPnt> {
    let surf = GeometryRegistry::global().face_surface(&face.0)?;
    let (u0, u1, v0, v1) = sample_bounds(surf.as_ref());
    let mut acc = occt_core::gp::GpXyz::zero();
    let mut count = 0.0;
    for i in 0..nu {
        for j in 0..nv {
            let u = u0 + (u1 - u0) * (i as f64 + 0.5) / nu as f64;
            let v = v0 + (v1 - v0) * (j as f64 + 0.5) / nv as f64;
            acc = acc.added(&surf.d0(u, v).coord);
            count += 1.0;
        }
    }
    if count == 0.0 {
        None
    } else {
        Some(GpPnt::from_xyz(&acc.divided(count)))
    }
}

/// Closest `(u, v)` parameters of `p` on a surface: an `nu × nv` grid seed plus
/// six rounds of axis-aligned bisection.
///
/// **UNPORTED**: this is the port's invented substitute. Two consumers remain,
/// both without an OCCT projection branch:
///
/// * [`crate::intcurvesurface`]'s general (non-quadric) path, where OCCT reads
///   the parameters off the surface **polyhedron**
///   (`IntCurveSurface_InterUtils.pxx:740-781`, `SectionPointToParameters`);
/// * [`crate::brepmesh::mesh_deflection_error`], a port-only mesh-QA metric
///   (OCCT measures deflection inside `BRepMesh`, not post hoc).
///
/// Every other former call site now goes through a faithful OCCT primitive:
/// `ElSLib::Parameters`, `Extrema_ExtPS` ([`occt_geom::geom_api::project_point_on_surface`]),
/// or the `ShapeConstruct_ProjectCurveOnSurface` projector.
pub(crate) fn surface_closest_params(
    s: &dyn Surface,
    p: &GpPnt,
    nu: usize,
    nv: usize,
) -> (f64, f64) {
    let (u0, u1, v0, v1) = sample_bounds(s);
    let mut best = (u0, v0);
    let mut best_d = f64::INFINITY;
    for i in 0..=nu {
        for j in 0..=nv {
            let u = u0 + (u1 - u0) * i as f64 / nu as f64;
            let v = v0 + (v1 - v0) * j as f64 / nv as f64;
            let d = s.d0(u, v).distance(p);
            if d < best_d {
                best_d = d;
                best = (u, v);
            }
        }
    }
    // Coarse-to-fine refinement (six halvings of the grid step).
    let (mut u, mut v) = best;
    let (mut hu, mut hv) = ((u1 - u0) / nu as f64, (v1 - v0) / nv as f64);
    for _ in 0..6 {
        hu *= 0.5;
        hv *= 0.5;
        let mut du = u;
        let mut dv = v;
        let mut dd = s.d0(u, v).distance(p);
        for &(su, sv) in &[(-1.0, 0.0), (1.0, 0.0), (0.0, -1.0), (0.0, 1.0)] {
            let nu2 = (u + su * hu).clamp(u0, u1);
            let nv2 = (v + sv * hv).clamp(v0, v1);
            let d = s.d0(nu2, nv2).distance(p);
            if d < dd {
                dd = d;
                du = nu2;
                dv = nv2;
            }
        }
        u = du;
        v = dv;
    }
    (u, v)
}

/// Approximate p-curve of an edge on a face: the projection of the edge's
/// curve onto the face surface, sampled at `samples` parameter values. Each
/// sample is a 2D `(u, v)` parameter of the face's surface.
///
/// This is the port's stand-in for `BRep_Tool::CurveOnSurface` when the face
/// carries no pcurve yet (`IntTools_FClass2d::Init` / `BRepTopAdaptor_FClass2d`
/// read the pcurve, and OCCT builds a missing one with
/// `ShapeConstruct_ProjectCurveOnSurface::Perform` through
/// `ShapeFix_Edge::FixAddPCurve`; `BRep_Tool.cxx` `CurveOnSurface` returns the
/// stored one). The projection is therefore the real OCCT projector
/// ([`crate::pcurve_full::project_curve_on_surface_perform`]): one single
/// source of `(u, v)`, the same one `shhealing::fix_add_pcurve` stores.
///
/// The old analytic shortcut (`pcurve_full::project_curve_on_surface`, torus
/// only) and the grid search below are NOT OCCT code paths - they are kept
/// only as a last resort so callers still get a polyline where OCCT would
/// return a null pcurve.
pub fn edge_pcurve_on_face(edge: &Edge, face: &Face, samples: usize) -> Vec<GpPnt2d> {
    let Some(surf) = GeometryRegistry::global().face_surface(&face.0) else { return Vec::new() };
    let Some(curve) = GeometryRegistry::global().edge_curve(&edge.0) else { return Vec::new() };
    let (a, b) = GeometryRegistry::global().edge_parameters(&edge.0);
    if !a.is_finite() || !b.is_finite() {
        return Vec::new();
    }
    let n = samples.max(1);
    let t_vals: Vec<f64> = (0..n)
        .map(|i| a + (b - a) * i as f64 / (n - 1).max(1) as f64)
        .collect();
    // `ShapeFix_Edge.cxx:499`: `preci = (prec > 0. ? prec : BRep_Tool::Tolerance(edge))`;
    // `cxx:521-531`: the end-vertex tolerances (`-1` = a null vertex).
    let preci = GeometryRegistry::global()
        .edge_tolerance(&edge.0)
        .max(occt_core::precision::CONFUSION);
    let mut cache = crate::pcurve_full::ProjectorCache::default();
    if let Some(c2d) = crate::pcurve_full::project_curve_on_surface_perform(
        curve.as_ref(),
        surf.as_ref(),
        a,
        b,
        preci,
        -1.0,
        -1.0,
        &mut cache,
    ) {
        return t_vals.iter().map(|&t| c2d.d0(t)).collect();
    }
    // `ShapeFix_Edge::FixAddPCurve` has no second projection: its projector
    // (`ShapeConstruct_ProjectCurveOnSurface`) is the only source of a missing
    // pcurve (`ShapeFix_Edge.cxx:519-546`). If the whole-curve projector above
    // found nothing, use the same per-sample primitive that projector calls
    // internally, `ShapeAnalysis_Surface::ValueOfUV` -> `Extrema_ExtPS`
    // (`GeomAPI_ProjectPointOnSurf`); OCCT has no grid search.
    t_vals
        .iter()
        .filter_map(|&u| {
            let p = curve.d0(u);
            occt_geom::geom_api::project_point_on_surface(surf.as_ref(), &p, preci)
                .map(|ps| GpPnt2d::new(ps.u, ps.v))
        })
        .collect()
}

/// Face UV bounds (`BRep_Tool::UVBounds(F, UMin, UMax, VMin, VMax)`,
/// `BRep_Tool.cxx`), i.e. the box of the face's pcurves, not the surface's
/// natural range:
///
///   BRep_Tool::UVBounds -> BRepTools::UVBounds(F, B) -> AddUVBounds
/// (`BRepTools.cxx:64-75`, `:126-160`); a face with no usable pcurve box falls
/// back to the surface bounds (`BRepTools.cxx:141-153`). This is also what
/// `BRepAdaptor_Surface(F)` restricts to (`BRepAdaptor_Surface.cxx:72-75`,
/// `BRepTools::UVBounds`), which is the adaptor `IntTools_FaceFace::Perform`
/// loads into its `GeomAdaptor_Surface` (`IntTools_FaceFace.cxx`).
pub fn face_uv_bounds(face: &Face) -> (f64, f64, f64, f64) {
    crate::brep_tools::uv_bounds(face)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    use crate::brep_tool::BRepTool;
    use crate::builder::TopoBuilder;
    use crate::primitives::BRepPrimBox;
    use occt_core::gp::GpAx3;
    use occt_geom::{GeomPlane, GeomSphere};

    #[test]
    fn box_faces_are_planar() {
        let b = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        for f in crate::topo_tools_full::faces_of(&b.solid.0) {
            assert!(face_is_planar(&f));
            let pln = face_plane(&f).expect("planar face has plane");
            let _ = pln;
            assert_eq!(classify_surface(
                GeometryRegistry::global().face_surface(&f.0).unwrap().as_ref()), SurfaceKind::Plane);
        }
    }

    #[test]
    fn sphere_face_classified() {
        let b = TopoBuilder::new();
        let face = b.make_face(Arc::new(GeomSphere::new(
            occt_core::gp::GpSphere::new(GpAx3::standard(), 2.0).unwrap(),
        )), &[]);
        let s = GeometryRegistry::global().face_surface(&face.0).unwrap();
        assert_eq!(classify_surface(s.as_ref()), SurfaceKind::Sphere);
        assert!(!face_is_planar(&face));
    }

    #[test]
    fn plane_reconstructed() {
        let b = TopoBuilder::new();
        let pln = GpPln::new(GpAx3::standard());
        let face = b.make_face(Arc::new(GeomPlane::new(pln.clone())), &[]);
        let got = face_plane(&face).expect("plane face");
        // Normal is Z in both.
        assert!(got.axis().direction().xyz().crossed(pln.axis().direction().xyz()).modulus() < 1e-9);
    }

    #[test]
    fn box_face_normal_and_centroid() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 2.0);
        let top = crate::topo_tools_full::faces_of(&b.solid.0)
            .into_iter()
            .find(|f| BRepTool::face_surface(&Face(f.0.clone()))
                .map(|s| (s.d0(0.0, 0.0).z() - 2.0).abs() < 1e-9).unwrap_or(false))
            .expect("top face");
        let n = face_normal(&Face(top.0.clone()), 0.0, 0.0).unwrap();
        assert!(n.xyz().z > 0.9, "top face normal should point +Z");
        let c = face_centroid(&Face(top.0.clone()), 8, 8).unwrap();
        assert!((c.z() - 2.0).abs() < 1e-6);
    }

    #[test]
    fn pcurve_projection_is_recoverable() {
        let b = TopoBuilder::new();
        let pln = GpPln::new(GpAx3::standard());
        let face = b.make_face(Arc::new(GeomPlane::new(pln)), &[]);
        let e = b.make_edge_segment(&GpPnt::new(0., 0., 0.), &GpPnt::new(1., 0., 0.));
        let pc = edge_pcurve_on_face(&e, &Face(face.0.clone()), 5);
        assert_eq!(pc.len(), 5);
        // On the Z=0 plane, u≈x, v≈y.
        assert!((pc[0].x() - 0.0).abs() < 1e-6);
        assert!((pc[4].x() - 1.0).abs() < 1e-6);
    }
}
