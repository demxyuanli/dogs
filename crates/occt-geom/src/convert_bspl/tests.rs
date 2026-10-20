use super::*;
use crate::bspline_curve::GeomBSplineCurve;
use crate::bspline_surface::GeomBSplineSurface;
use crate::surface::Surface;
use occt_core::gp::GpPnt;

fn degree3_multiknot_curve() -> GeomBSplineCurve {
    // Cubic; interior 0.5 (mult 1) keeps C2, interior 0.75 (mult 3) is C0.
    let poles = vec![
        GpPnt::new(0.0, 0.0, 0.0),
        GpPnt::new(1.0, 0.0, 0.0),
        GpPnt::new(2.0, 0.0, 0.0),
        GpPnt::new(3.0, 0.0, 0.0),
        GpPnt::new(4.0, 0.0, 0.0),
        GpPnt::new(5.0, 0.0, 0.0),
        GpPnt::new(6.0, 0.0, 0.0),
        GpPnt::new(7.0, 0.0, 0.0),
    ];
    let knots = vec![0.0, 0.0, 0.0, 0.0, 0.5, 0.75, 0.75, 0.75, 1.0, 1.0, 1.0, 1.0];
    GeomBSplineCurve::new(poles, knots, 3).unwrap()
}

#[test]
fn knot_splitting_curve_finds_c0_knot() {
    let c = degree3_multiknot_curve();
    let splits = knot_splitting_curve(&c);
    // Distinct knots: [0 (m4), 0.5 (m1), 0.75 (m3), 1 (m4)].
    // Only 0.75 (index 2) drops below C1.
    assert_eq!(splits, vec![0, 2, 3], "splits {splits:?}");
}

#[test]
fn knot_splitting_curve_uniform_stays_whole() {
    let c = GeomBSplineCurve::new(
        vec![
            GpPnt::new(0., 0., 0.),
            GpPnt::new(1., 1., 0.),
            GpPnt::new(2., 0., 0.),
            GpPnt::new(3., -1., 0.),
            GpPnt::new(4., 0., 0.),
            GpPnt::new(5., 1., 0.),
        ],
        vec![0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 4.0],
        2,
    )
    .unwrap();
    // All interior knots have multiplicity 1 → degree−mult = C1 → no split.
    let splits = knot_splitting_curve(&c);
    assert_eq!(splits, vec![0, 4], "splits {splits:?}");
}

#[test]
fn knot_splitting_surface_u_and_v() {
    // Degree 2 in u (interior mult 1 → C1), degree 2 in v (interior mult 2 → C0).
    // u: 4 poles (7 knots), v: 5 poles (8 knots).
    let poles = vec![vec![GpPnt::new(0., 0., 0.); 5]; 4];
    let ku = vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0];
    let kv = vec![0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0];
    let s = GeomBSplineSurface::new(poles, ku, kv, 2, 2).unwrap();
    let (us, vs) = knot_splitting_surface(&s);
    // u: distinct [0,0.5,1], mult 1 → no split → [0, 2].
    // v: distinct [0,0.5,1], mult 2 = degree → split at index 1 → [0,1,2].
    assert_eq!(us, vec![0, 2], "u splits {us:?}");
    assert_eq!(vs, vec![0, 1, 2], "v splits {vs:?}");
}

fn bilinear_surface() -> GeomBSplineSurface {
    let poles = vec![
        vec![GpPnt::new(0., 0., 0.), GpPnt::new(0., 1., 0.)],
        vec![GpPnt::new(1., 0., 0.), GpPnt::new(1., 1., 0.)],
    ];
    let ku = vec![0.0, 0.0, 1.0, 1.0];
    let kv = vec![0.0, 0.0, 1.0, 1.0];
    GeomBSplineSurface::new(poles, ku, kv, 1, 1).unwrap()
}

#[test]
fn surface_to_bezier_single_patch_bilinear() {
    let s = bilinear_surface();
    let patches = bspline_surface_to_bezier_surface(&s);
    assert_eq!(patches.len(), 1);
    let p = &patches[0];
    assert_eq!(p.poles.len(), 2);
    assert_eq!(p.poles[0].len(), 2);
    assert!((p.u0 - 0.0).abs() < 1e-12 && (p.u1 - 1.0).abs() < 1e-12);
    assert!((p.v0 - 0.0).abs() < 1e-12 && (p.v1 - 1.0).abs() < 1e-12);
    // Patch corners equal the surface corners.
    assert!(p.poles[0][0].distance(&s.d0(0.0, 0.0)) < 1e-9);
    assert!(p.poles[0][1].distance(&s.d0(0.0, 1.0)) < 1e-9);
    assert!(p.poles[1][0].distance(&s.d0(1.0, 0.0)) < 1e-9);
    assert!(p.poles[1][1].distance(&s.d0(1.0, 1.0)) < 1e-9);
}

#[test]
fn surface_to_bezier_patches_match_at_shared_corners() {
    // Degree-2 surface in u with an interior knot → two u-patches; the
    // shared corner must coincide with the original surface.
    let poles = vec![
        vec![GpPnt::new(0., 0., 0.), GpPnt::new(0., 1., 0.)],
        vec![GpPnt::new(0.5, 0.5, 0.), GpPnt::new(0.5, 0.5, 0.)],
        vec![GpPnt::new(1.5, 0.5, 0.), GpPnt::new(1.5, 0.5, 0.)],
        vec![GpPnt::new(2., 0., 0.), GpPnt::new(2., 1., 0.)],
    ];
    let ku = vec![0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0];
    let kv = vec![0.0, 0.0, 1.0, 1.0];
    let s = GeomBSplineSurface::new(poles, ku, kv, 2, 1).unwrap();

    let patches = bspline_surface_to_bezier_surface(&s);
    // u: distinct [0,1,2] → 2 u-patches; v: distinct [0,1] → 1 v-patch.
    assert_eq!(patches.len(), 2, "patches {patches:?}");
    let (p0, p1) = (&patches[0], &patches[1]);
    assert!((p0.u1 - p1.u0).abs() < 1e-12, "u1 {} vs u0 {}", p0.u1, p1.u0);
    // Shared corner between the patches equals the surface at the breakpoint.
    let shared_u = p0.u1;
    for v in [0.0, 1.0] {
        let got = p0.poles[p0.poles.len() - 1][if v < 0.5 { 0 } else { 1 }];
        let expected = s.d0(shared_u, v);
        assert!(got.distance(&expected) < 1e-9, "corner v={v}: {got:?} vs {expected:?}");
        let got2 = p1.poles[0][if v < 0.5 { 0 } else { 1 }];
        assert!(got2.distance(&expected) < 1e-9);
    }
    // Every patch's four corners sit on the original surface.
    for p in &patches {
        let (nu, nv) = (p.poles.len() - 1, p.poles[0].len() - 1);
        for &(i, j) in &[(0, 0), (0, nv), (nu, 0), (nu, nv)] {
            let u = if i == 0 { p.u0 } else { p.u1 };
            let v = if j == 0 { p.v0 } else { p.v1 };
            assert!(p.poles[i][j].distance(&s.d0(u, v)) < 1e-9, "corner ({i},{j})");
        }
    }
}
