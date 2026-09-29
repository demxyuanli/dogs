//! FaceFace pair tests (moved out of `int_face_face.rs`).

use super::*;
use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::gp::{GpAx1, GpAx3, GpCone, GpCylinder, GpDir, GpPln, GpPnt, GpSphere, GpTorus};
use occt_geom::{GeomCone, GeomCylinder, GeomSphere, GeomTorus, Surface};

use crate::builder::TopoBuilder;
use crate::intpatch;
use crate::inttools_data::CurveKind;
use crate::shape::Face;

const TOL: f64 = 1e-6;

fn plane_z(z: f64) -> GpPln {
    GpPln::new(
        GpAx3::new(
            GpPnt::new(0.0, 0.0, z),
            GpDir::new(0.0, 0.0, 1.0).unwrap(),
            &GpDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap(),
    )
}

fn plane_face(pln: &GpPln) -> Face {
    TopoBuilder::new().make_face_plane(pln)
}

fn sphere_face(center: GpPnt, r: f64) -> Face {
    let ax3 = GpAx3::new(
        center,
        GpDir::new(0.0, 0.0, 1.0).unwrap(),
        &GpDir::new(1.0, 0.0, 0.0).unwrap(),
    )
    .unwrap();
    TopoBuilder::new().make_face(Arc::new(GeomSphere::new(GpSphere::new(ax3, r).unwrap())), &[])
}

fn cylinder_face(radius: f64) -> Face {
    let ax3 = GpAx3::new(
        GpPnt::zero(),
        GpDir::new(0.0, 0.0, 1.0).unwrap(),
        &GpDir::new(1.0, 0.0, 0.0).unwrap(),
    )
    .unwrap();
    TopoBuilder::new().make_face(Arc::new(GeomCylinder::new(GpCylinder::new(ax3, radius).unwrap())), &[])
}

/// A `GpCone` whose apex (radius-0 point) is at the origin, axis +Z, with
/// location ring radius 1. The placement is shifted to (0,0,1/tan α) so the
/// apex lands on the origin under the OCCT placement convention.
fn cone_apex_origin(a: f64) -> GpCone {
    let loc = GpPnt::new(0.0, 0.0, 1.0 / a.tan());
    let ax3 = GpAx3::new(
        loc,
        GpDir::new(0.0, 0.0, 1.0).unwrap(),
        &GpDir::new(1.0, 0.0, 0.0).unwrap(),
    )
    .unwrap();
    GpCone::new(ax3, 1.0, a).unwrap()
}

/// A cone face, apex at the origin, axis +Z, with the given semi-angle.
fn cone_face(semi_angle: f64) -> Face {
    TopoBuilder::new().make_face(Arc::new(GeomCone::new(cone_apex_origin(semi_angle))), &[])
}

/// A torus face centered at the origin, axis +Z.
fn torus_face(major: f64, minor: f64) -> Face {
    let torus = GpTorus::new(GpAx3::standard(), major, minor).unwrap();
    TopoBuilder::new().make_face(Arc::new(GeomTorus::new(torus)), &[])
}

/// A plane face through `origin` with the given unit normal (x-dir is
/// `+X`, valid for normals with zero x-component).
fn plane_face_normal(origin: GpPnt, normal: GpDir) -> Face {
    let ax3 = GpAx3::new(origin, normal, &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
    TopoBuilder::new().make_face_plane(&GpPln::new(ax3))
}

/// Sample a FaceFaceCurve and assert every sample lies on both faces'
/// surfaces (within `eps`).
fn assert_points_on_both(curve: &FaceFaceCurve, sa: &dyn Surface, sb: &dyn Surface, eps: f64) {
    let (a, b) = (curve.range.first, curve.range.last);
    assert!(a.is_finite() && b.is_finite(), "range {a}..{b}");
    for i in 0..=16 {
        let t = a + (b - a) * i as f64 / 16.0;
        let p = curve.curve.d0(t);
        assert!(intpatch::distance_to_surface(&p, sa) < eps, "on surface a: {p:?}");
        assert!(intpatch::distance_to_surface(&p, sb) < eps, "on surface b: {p:?}");
    }
}

#[test]
fn not_done_before_perform() {
    let ff = FaceFace::new();
    assert!(!ff.is_done());
    assert!(ff.result().is_empty());
}

#[test]
fn perform_requires_faces() {
    let mut ff = FaceFace::new();
    ff.set_face1(plane_face(&plane_z(0.0)));
    assert!(ff.perform().is_err(), "missing face2 must fail");
}

#[test]
fn plane_sphere_intersects_in_circle() {
    let pln = plane_face(&plane_z(0.5));
    let sph = sphere_face(GpPnt::zero(), 1.0);
    let mut ff = FaceFace::new();
    ff.set_face1(pln);
    ff.set_face2(sph);
    ff.set_tolerance(TOL);
    ff.perform().expect("perform");
    assert!(ff.is_done());
    let res = ff.result();
    assert_eq!(res.nb_curves(), 1, "one circle of intersection");
    let c = res.curve(0);
    assert_eq!(c.kind, CurveKind::Circle, "plane∩sphere is a circle");
    let s_pln: Arc<dyn Surface> = Arc::new(occt_geom::GeomPlane::new(plane_z(0.5)));
    let s_sph: Arc<dyn Surface> = Arc::new(GeomSphere::new(
        GpSphere::new(
            GpAx3::new(
                GpPnt::zero(),
                GpDir::new(0.0, 0.0, 1.0).unwrap(),
                &GpDir::new(1.0, 0.0, 0.0).unwrap(),
            )
            .unwrap(),
            1.0,
        )
        .unwrap(),
    ));
    for i in 0..=16 {
        let t = c.range.first + (c.range.last - c.range.first) * i as f64 / 16.0;
        let p = c.curve.d0(t);
        assert!((p.distance(&GpPnt::zero()) - 1.0).abs() < 1e-6, "on sphere: {p:?}");
        assert!((p.z() - 0.5).abs() < 1e-6, "in plane: {p:?}");
    }
    assert_points_on_both(c, s_pln.as_ref(), s_sph.as_ref(), 1e-4);
}

#[test]
fn sphere_sphere_intersects_in_circle() {
    let s1 = sphere_face(GpPnt::zero(), 1.0);
    let s2 = sphere_face(GpPnt::new(1.5, 0.0, 0.0), 1.0);
    let mut ff = FaceFace::new();
    ff.set_face1(s1);
    ff.set_face2(s2);
    ff.set_tolerance(TOL);
    ff.perform().expect("perform");
    let res = ff.result();
    assert_eq!(res.nb_curves(), 1, "one circle of intersection");
    let c = res.curve(0);
    assert_eq!(c.kind, CurveKind::Circle);
    let expected_r = (1.0f64 - 0.75 * 0.75).sqrt();
    for i in 0..=16 {
        let t = c.range.first + (c.range.last - c.range.first) * i as f64 / 16.0;
        let p = c.curve.d0(t);
        assert!((p.x() - 0.75).abs() < 1e-6, "plane x: {p:?}");
        let r = GpPnt::new(0.0, p.y(), p.z()).distance(&GpPnt::zero());
        assert!((r - expected_r).abs() < 1e-6, "radius {r} (expected {expected_r})");
    }
}

#[test]
fn two_box_faces_intersect_in_line_segment() {
    let b = crate::brep_extrema::test_box::unit_box();
    let bottom = b.faces[0].clone(); // z = 0
    let front = b.faces[2].clone(); // y = 0
    let mut ff = FaceFace::new();
    ff.set_face1(bottom);
    ff.set_face2(front);
    ff.set_tolerance(TOL);
    ff.perform().expect("perform");
    let res = ff.result();
    assert_eq!(res.nb_curves(), 1, "adjacent box faces meet in one line");
    let c = res.curve(0);
    assert_eq!(c.kind, CurveKind::Line);
    // The line is the box edge (0,0,0)-(1,0,0).
    let (a, b) = (c.range.first, c.range.last);
    let p0 = c.curve.d0(a);
    let p1 = c.curve.d0(b);
    assert!((p0.y()).abs() < 1e-6 && (p0.z()).abs() < 1e-6, "start {p0:?}");
    assert!((p1.y()).abs() < 1e-6 && (p1.z()).abs() < 1e-6, "end {p1:?}");
    assert!((p1.x() - p0.x()).abs() > 0.5, "segment spans the box edge");
    assert!(p0.x() >= -1e-6 && p1.x() <= 1.0 + 1e-6, "within x∈[0,1]");
}

#[test]
fn plane_cylinder_intersects_in_circle() {
    let pln = plane_face(&plane_z(1.0));
    let cyl = cylinder_face(2.0);
    let mut ff = FaceFace::new();
    ff.set_face1(pln);
    ff.set_face2(cyl);
    ff.set_tolerance(TOL);
    ff.perform().expect("perform");
    let res = ff.result();
    assert_eq!(res.nb_curves(), 1, "plane ⊥ cylinder axis meets in one circle");
    let c = res.curve(0);
    assert_eq!(c.kind, CurveKind::Circle);
    for i in 0..=16 {
        let t = c.range.first + (c.range.last - c.range.first) * i as f64 / 16.0;
        let p = c.curve.d0(t);
        assert!((p.z() - 1.0).abs() < 1e-6, "in plane: {p:?}");
        let r = GpPnt::new(p.x(), p.y(), 0.0).distance(&GpPnt::zero());
        assert!((r - 2.0).abs() < 1e-6, "radius {r}");
    }
}

#[test]
fn plane_cylinder_parallel_generatrices() {
    // Plane x = 0.5 parallel to the Z-axis cylinder radius 1 → generatrix
    // lines at y = ±√0.75. The FaceFace curve must agree with the intpatch
    // closed form it delegates to (对拍).
    let pln = GpPln::new(
        GpAx3::new(
            GpPnt::new(0.5, 0.0, 0.0),
            GpDir::new(1.0, 0.0, 0.0).unwrap(),
            &GpDir::new(0.0, 1.0, 0.0).unwrap(),
        )
        .unwrap(),
    );
    // Reference: `IntAna_QuadQuadGeo::Perform(gp_Pln, gp_Cylinder)` (`intana`),
    // which yields `TwoLines` for the parallel case (T-28 step 4: the closed
    // form now has a single home, so the reference source moved from
    // `intpatch::intersect_plane_cylinder` to `intana`).
    let cyl = occt_core::gp::GpCylinder::new(
        GpAx3::new(
            GpPnt::zero(),
            GpDir::new(0.0, 0.0, 1.0).unwrap(),
            &GpDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap(),
        1.0,
    )
    .unwrap();
    let (l1, l2) = match occt_geom::intana::quadric_quadric_plane_cylinder(&pln, &cyl, 1e-12, 1e-7) {
        occt_geom::intana::QuadricIntersection::TwoLines(a, b) => (a, b),
        other => panic!("IntAna(gp_Pln, gp_Cylinder) should give TwoLines, got {other:?}"),
    };
    let dist_to_line = |p: &GpPnt, l: &occt_core::gp::GpLin| {
        let d = occt_core::gp::GpVec::from_xyz(l.direction().xyz()).normalized();
        let v = occt_core::gp::GpVec::from_pnts(&l.location(), p);
        v.crossed(&d).magnitude()
    };

    let mut ff = FaceFace::new();
    ff.set_face1(plane_face(&pln));
    ff.set_face2(cylinder_face(1.0));
    ff.set_tolerance(TOL);
    ff.perform().expect("perform");
    let res = ff.result();
    assert_eq!(res.nb_curves(), 2, "plane ∥ cylinder: two generatrix lines (IntAna TwoLines)");
    let mut hit = [false, false];
    for i in 0..2 {
        let c = res.curve(i);
        assert_eq!(c.kind, CurveKind::Line);
        let (a, b) = (c.range.first, c.range.last);
        let mid_pt = c.curve.d0(0.5 * (a + b));
        let (d1, d2) = (dist_to_line(&mid_pt, &l1), dist_to_line(&mid_pt, &l2));
        assert!(
            d1 < 1e-6 || d2 < 1e-6,
            "FaceFace curve {i} at {mid_pt:?} is on neither generatrix (d1={d1}, d2={d2})"
        );
        if d1 <= d2 {
            hit[0] = true;
        } else {
            hit[1] = true;
        }
    }
    assert!(hit[0] && hit[1], "both generatrices must be present");
}

#[test]
fn plane_cone_intersects_in_circle() {
    // Cone with semi-angle atan(0.5), geometric tip at the origin; plane
    // z = 1 perpendicular to the axis cuts a circle of radius tan(atan0.5)
    // = 0.5 at (0,0,1).
    let pln = plane_face(&plane_z(1.0));
    let cone = cone_face(0.5f64.atan());
    let mut ff = FaceFace::new();
    ff.set_face1(pln);
    ff.set_face2(cone);
    ff.set_tolerance(TOL);
    ff.perform().expect("perform");
    let res = ff.result();
    assert_eq!(res.nb_curves(), 1, "plane ⊥ cone axis meets in one circle");
    let c = res.curve(0);
    assert_eq!(c.kind, CurveKind::Circle);
    for i in 0..=16 {
        let t = c.range.first + (c.range.last - c.range.first) * i as f64 / 16.0;
        let p = c.curve.d0(t);
        assert!((p.z() - 1.0).abs() < 1e-6, "in plane: {p:?}");
        let r = GpPnt::new(p.x(), p.y(), 0.0).distance(&GpPnt::zero());
        assert!((r - 0.5).abs() < 1e-6, "radius {r}");
    }
    let s_cone: Arc<dyn Surface> = Arc::new(GeomCone::new(cone_apex_origin(0.5f64.atan())));
    let s_pln: Arc<dyn Surface> = Arc::new(occt_geom::GeomPlane::new(plane_z(1.0)));
    assert_points_on_both(c, s_pln.as_ref(), s_cone.as_ref(), 1e-4);
}

#[test]
fn plane_cone_intersects_in_ellipse() {
    // Oblique plane (not through the apex, not perpendicular to the axis,
    // not parallel to a generatrix) → exact ellipse, sampled on both faces.
    let pln = plane_face_normal(GpPnt::new(0.0, 0.0, 2.0), GpDir::new(0.0, 0.3, 0.954).unwrap());
    let cone = cone_face(0.5f64.atan());
    let mut ff = FaceFace::new();
    ff.set_face1(pln);
    ff.set_face2(cone);
    ff.set_tolerance(TOL);
    ff.perform().expect("perform");
    let res = ff.result();
    assert_eq!(res.nb_curves(), 1, "one conic section");
    let c = res.curve(0);
    assert_eq!(c.kind, CurveKind::Ellipse, "oblique cut is an ellipse");
    let s_cone: Arc<dyn Surface> = Arc::new(GeomCone::new(cone_apex_origin(0.5f64.atan())));
    let s_pln: Arc<dyn Surface> = Arc::new(occt_geom::GeomPlane::new(GpPln::new(
        GpAx3::new(GpPnt::new(0.0, 0.0, 2.0), GpDir::new(0.0, 0.3, 0.954).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap(),
    )));
    assert_points_on_both(c, s_pln.as_ref(), s_cone.as_ref(), 1e-4);
}

#[test]
fn plane_torus_axis_in_plane_two_circles() {
    // Torus R=3 r=1, plane y = 0 contains the axis → two circles of radius
    // 1 centered at (±3, 0, 0).
    let pln = plane_face_normal(GpPnt::zero(), GpDir::new(0.0, 1.0, 0.0).unwrap());
    let tor = torus_face(3.0, 1.0);
    let mut ff = FaceFace::new();
    ff.set_face1(pln);
    ff.set_face2(tor);
    ff.set_tolerance(TOL);
    ff.perform().expect("perform");
    let res = ff.result();
    assert_eq!(res.nb_curves(), 2, "axis-in-plane cut gives two minor circles");
    for c in res.curves() {
        assert_eq!(c.kind, CurveKind::Circle);
        // The circle's center (midpoint of two opposite points) is on the
        // major ring: (±3, 0, 0).
        let (a, b) = (c.range.first, c.range.last);
        let p0 = c.curve.d0(a);
        let p1 = c.curve.d0(0.5 * (a + b)); // opposite point (circle is 2π-periodic)
        let center = mid(&p0, &p1);
        assert!(center.y().abs() < 1e-6, "center in plane: {center:?}");
        assert!((GpPnt::new(center.x(), 0.0, center.z()).distance(&GpPnt::zero()) - 3.0).abs() < 1e-6,
            "center on major ring {center:?}");
        for i in 0..=16 {
            let t = a + (b - a) * i as f64 / 16.0;
            let p = c.curve.d0(t);
            assert!(p.y().abs() < 1e-6, "in plane: {p:?}");
        }
    }
    let s_tor: Arc<dyn Surface> = Arc::new(GeomTorus::new(GpTorus::new(GpAx3::standard(), 3.0, 1.0).unwrap()));
    let s_pln: Arc<dyn Surface> = Arc::new(occt_geom::GeomPlane::new(GpPln::new(
        GpAx3::new(GpPnt::zero(), GpDir::new(0.0, 1.0, 0.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap(),
    )));
    for c in res.curves() {
        assert_points_on_both(c, s_pln.as_ref(), s_tor.as_ref(), 1e-4);
    }
}

#[test]
fn plane_torus_perpendicular_cut_two_circles() {
    // Torus R=3 r=1, plane z = 0 perpendicular to the axis through the
    // center → two circles of radius 4 and 2.
    let pln = plane_face(&plane_z(0.0));
    let tor = torus_face(3.0, 1.0);
    let mut ff = FaceFace::new();
    ff.set_face1(pln);
    ff.set_face2(tor);
    ff.set_tolerance(TOL);
    ff.perform().expect("perform");
    let res = ff.result();
    assert_eq!(res.nb_curves(), 2, "perpendicular cut through the center gives two circles");
    let mut radii: Vec<f64> = res
        .curves()
        .iter()
        .map(|c| c.curve.d0(0.0).distance(&GpPnt::zero()))
        .collect();
    radii.sort_by(f64::total_cmp);
    assert!((radii[0] - 2.0).abs() < 1e-6, "minor circle radius {}", radii[0]);
    assert!((radii[1] - 4.0).abs() < 1e-6, "major circle radius {}", radii[1]);
    let s_tor: Arc<dyn Surface> = Arc::new(GeomTorus::new(GpTorus::new(GpAx3::standard(), 3.0, 1.0).unwrap()));
    let s_pln: Arc<dyn Surface> = Arc::new(occt_geom::GeomPlane::new(plane_z(0.0)));
    for c in res.curves() {
        assert_points_on_both(c, s_pln.as_ref(), s_tor.as_ref(), 1e-4);
    }
}

#[test]
fn general_bspline_sphere_uses_fallback() {
    // A genuinely curved B-spline patch (classifies Other, not Plane) × a
    // sphere goes through the sampling tracer. The traced polyline points
    // must lie near both surfaces.
    let mut grid: Vec<Vec<GpPnt>> = Vec::new();
    for i in 0..=3usize {
        let mut row = Vec::new();
        for j in 0..=3usize {
            let u = i as f64 / 3.0;
            let v = j as f64 / 3.0;
            row.push(GpPnt::new(u, v, u * u + v * v));
        }
        grid.push(row);
    }
    let patch_surf = intpatch::make_bspline_surface_from_grid(&grid, 1, 1).expect("fit");
    let patch_surf_check = patch_surf.clone();
    let patch_face = TopoBuilder::new().make_face(patch_surf, &[]);
    let sph = sphere_face(GpPnt::new(0.5, 0.5, 0.5), 1.0);

    let mut ff = FaceFace::new();
    ff.set_face1(patch_face);
    ff.set_face2(sph);
    ff.set_tolerance(0.05);
    ff.perform().expect("perform");
    let res = ff.result();
    assert!(res.nb_curves() >= 1, "tracer found an intersection");
    let s_sph: Arc<dyn Surface> = Arc::new(GeomSphere::new(
        GpSphere::new(
            GpAx3::new(
                GpPnt::new(0.5, 0.5, 0.5),
                GpDir::new(0.0, 0.0, 1.0).unwrap(),
                &GpDir::new(1.0, 0.0, 0.0).unwrap(),
            )
            .unwrap(),
            1.0,
        )
        .unwrap(),
    ));
    for c in res.curves() {
        let (a, b) = (c.range.first, c.range.last);
        for i in 0..=12 {
            let t = a + (b - a) * i as f64 / 12.0;
            let p = c.curve.d0(t);
            assert!(intpatch::distance_to_surface(&p, patch_surf_check.as_ref()) < 0.2, "near patch: {p:?}");
            assert!(intpatch::distance_to_surface(&p, s_sph.as_ref()) < 0.2, "near sphere: {p:?}");
        }
    }
}

#[test]
fn disjoint_sphere_plane_yields_no_curves() {
    let pln = plane_face(&plane_z(2.0)); // above the unit sphere
    let sph = sphere_face(GpPnt::zero(), 1.0);
    let mut ff = FaceFace::new();
    ff.set_face1(pln);
    ff.set_face2(sph);
    ff.set_tolerance(TOL);
    ff.perform().expect("perform");
    assert!(ff.result().is_empty(), "no intersection when the plane misses the sphere");
}

#[test]
fn coplanar_planes_are_tangent() {
    let p1 = plane_face(&plane_z(0.0));
    let p2 = plane_face(&plane_z(0.0));
    let mut ff = FaceFace::new();
    ff.set_face1(p1);
    ff.set_face2(p2);
    ff.set_tolerance(TOL);
    ff.perform().expect("perform");
    assert!(ff.tangent_faces(), "coincident planes are tangent faces");
    assert!(ff.result().is_empty());
}

#[test]
fn parallel_planes_do_not_intersect() {
    let p1 = plane_face(&plane_z(0.0));
    let p2 = plane_face(&plane_z(1.0));
    let mut ff = FaceFace::new();
    ff.set_face1(p1);
    ff.set_face2(p2);
    ff.set_tolerance(TOL);
    ff.perform().expect("perform");
    assert!(!ff.tangent_faces());
    assert!(ff.result().is_empty());
}

#[test]
fn result_is_sorted_and_deduplicated() {
    let pln = plane_face(&plane_z(0.5));
    let sph = sphere_face(GpPnt::zero(), 1.0);
    let mut ff = FaceFace::new();
    ff.set_face1(pln);
    ff.set_face2(sph);
    ff.perform().expect("perform");
    let res = ff.result();
    let ranges: Vec<f64> = res.curves().iter().map(|c| c.range.first).collect();
    let mut sorted = ranges.clone();
    sorted.sort_by(f64::total_cmp);
    assert_eq!(ranges, sorted, "curves sorted by first parameter");
}

#[test]
fn context_is_optional_noop() {
    let mut ff = FaceFace::new();
    ff.set_context(FaceFaceContext);
    assert!(ff.tangent_faces() == false);
    assert!(!ff.is_done());
}

#[test]
fn cylinder_cylinder_two_lines() {
    // Two unit Z-axis cylinders, axes through (0,0,0) and (0.5,0,0): the
    // section is two generatrix lines at x=0.25, y=±√0.9375.
    let c1 = cylinder_face(1.0);
    let ax3 = GpAx3::new(
        GpPnt::new(0.5, 0.0, 0.0),
        GpDir::new(0.0, 0.0, 1.0).unwrap(),
        &GpDir::new(1.0, 0.0, 0.0).unwrap(),
    )
    .unwrap();
    let c2 = TopoBuilder::new().make_face(
        Arc::new(GeomCylinder::new(GpCylinder::new(ax3, 1.0).unwrap())),
        &[],
    );
    let mut ff = FaceFace::new();
    ff.set_face1(c1);
    ff.set_face2(c2);
    ff.set_tolerance(TOL);
    ff.perform().expect("perform");
    let res = ff.result();
    assert_eq!(res.nb_curves(), 2, "two generatrix lines");
    let h = (0.9375f64).sqrt();
    for i in 0..2 {
        let c = res.curve(i);
        assert_eq!(c.kind, CurveKind::Line);
        for k in 0..=4 {
            let t = c.range.first + (c.range.last - c.range.first) * k as f64 / 4.0;
            let p = c.curve.d0(t);
            assert!((p.x() - 0.25).abs() < 1e-6, "x {}", p.x());
            assert!((p.y().abs() - h).abs() < 1e-6, "y {}", p.y());
        }
    }
}

#[test]
fn cylinder_sphere_two_circles() {
    // Unit sphere at the origin, Z-axis cylinder r=0.5 through it: circles
    // at z=±√0.75, radius 0.5.
    let cyl = cylinder_face(0.5);
    let sph = sphere_face(GpPnt::zero(), 1.0);
    let mut ff = FaceFace::new();
    ff.set_face1(cyl);
    ff.set_face2(sph);
    ff.set_tolerance(TOL);
    ff.perform().expect("perform");
    let res = ff.result();
    assert_eq!(res.nb_curves(), 2, "two circles");
    for i in 0..2 {
        let c = res.curve(i);
        assert_eq!(c.kind, CurveKind::Circle);
        for k in 0..8 {
            let t = c.range.first + (c.range.last - c.range.first) * k as f64 / 8.0;
            let p = c.curve.d0(t);
            let xy = GpPnt::new(p.x(), p.y(), 0.0).distance(&GpPnt::zero());
            assert!((xy - 0.5).abs() < 1e-6, "cylinder radius {xy}");
            assert!((p.distance(&GpPnt::zero()) - 1.0).abs() < 1e-6, "sphere radius");
        }
    }
}

#[test]
fn sphere_cone_two_circles() {
    // Cone (apex at the origin, axis +Z, 30°) with a unit-less sphere at
    // (0,0,3) radius 2: two circles, both on the cone and the sphere.
    let cone = cone_face(PI / 6.0);
    let sph = sphere_face(GpPnt::new(0.0, 0.0, 3.0), 2.0);
    let mut ff = FaceFace::new();
    ff.set_face1(cone);
    ff.set_face2(sph);
    ff.set_tolerance(TOL);
    ff.perform().expect("perform");
    let res = ff.result();
    assert_eq!(res.nb_curves(), 2, "two circles");
    let sa: Arc<dyn Surface> = Arc::new(GeomCone::new(cone_apex_origin(PI / 6.0)));
    let sb: Arc<dyn Surface> = Arc::new(GeomSphere::new(GpSphere::new(
        GpAx3::new(
            GpPnt::new(0.0, 0.0, 3.0),
            GpDir::new(0.0, 0.0, 1.0).unwrap(),
            &GpDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap(),
        2.0,
    )
    .unwrap()));
    for i in 0..2 {
        let c = res.curve(i);
        assert_eq!(c.kind, CurveKind::Circle);
        assert_points_on_both(c, sa.as_ref(), sb.as_ref(), 1e-4);
    }
}
