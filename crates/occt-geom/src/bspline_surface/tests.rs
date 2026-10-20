//! Unit tests of `Geom_BSplineSurface` (moved verbatim from the former
//! single-file `bspline_surface.rs`).

use super::*;
use crate::surface::Surface;
use occt_core::gp::{GpPnt, GpTrsf, GpVec};

#[test]
fn uniform_knot_sanity() {
    let (ku, kv) = bspline_surface_uniform_knots(4, 3, 2, 1);
    assert_eq!(ku.len(), 4 + 2 + 1);
    assert_eq!(kv.len(), 3 + 1 + 1);
    assert_eq!(ku[0], 0.0);
    assert_eq!(*ku.last().unwrap(), 1.0);
    assert!(ku.windows(2).all(|w| w[0] <= w[1]));
    assert!(kv.windows(2).all(|w| w[0] <= w[1]));
    // End multiplicity degree + 1.
    assert_eq!(ku[..3].iter().filter(|&&x| x == 0.0).count(), 3);
    assert_eq!(ku[ku.len() - 3..].iter().filter(|&&x| x == 1.0).count(), 3);
    assert_eq!(kv[..2].iter().filter(|&&x| x == 0.0).count(), 2);
    // Interior knot is midway for a 4-pole degree-2 vector.
    assert!((ku[3] - 0.5).abs() < 1e-12);
}

#[test]
fn eval_at_clamped_corner_is_first_pole() {
    let poles = vec![
        vec![GpPnt::new(1.0, 2.0, 3.0), GpPnt::new(4.0, 5.0, 6.0)],
        vec![GpPnt::new(7.0, 8.0, 9.0), GpPnt::new(10.0, 11.0, 12.0)],
    ];
    let (ku, kv) = bspline_surface_uniform_knots(2, 2, 1, 1);
    let s = GeomBSplineSurface::new(poles, ku, kv, 1, 1).unwrap();
    let p = s.d0(0.0, 0.0);
    assert!((p.x() - 1.0).abs() < 1e-12);
    assert!((p.y() - 2.0).abs() < 1e-12);
    assert!((p.z() - 3.0).abs() < 1e-12);
    // Opposite corner: last pole.
    let q = s.d0(1.0, 1.0);
    assert!((q.x() - 10.0).abs() < 1e-12);
    assert!((q.y() - 11.0).abs() < 1e-12);
}

#[test]
fn degree1_plane_grid() {
    // z = 2x + y + 1 sampled at the four corners.
    let poles = vec![
        vec![GpPnt::new(0.0, 0.0, 1.0), GpPnt::new(0.0, 1.0, 2.0)],
        vec![GpPnt::new(1.0, 0.0, 3.0), GpPnt::new(1.0, 1.0, 4.0)],
    ];
    let (ku, kv) = bspline_surface_uniform_knots(2, 2, 1, 1);
    let s = GeomBSplineSurface::new(poles, ku, kv, 1, 1).unwrap();
    let p = s.d0(0.5, 0.5);
    assert!((p.x() - 0.5).abs() < 1e-12);
    assert!((p.y() - 0.5).abs() < 1e-12);
    assert!((p.z() - 2.5).abs() < 1e-12, "z = {}", p.z());
}

#[test]
fn xy_grid_fit_matches_analytic() {
    // Grid points on z = u·v; degree-1 surface reproduces the bilinear field.
    let (nu, nv) = (4, 4);
    let points: Vec<Vec<GpPnt>> = (0..nu)
        .map(|i| {
            (0..nv)
                .map(|j| {
                    let u = i as f64 / (nu - 1) as f64;
                    let v = j as f64 / (nv - 1) as f64;
                    GpPnt::new(u, v, u * v)
                })
                .collect()
        })
        .collect();
    let s = fit_surface_grid(&points, 1, 1).unwrap();
    let p = s.d0(0.3, 0.7);
    assert!((p.x() - 0.3).abs() < 1e-6);
    assert!((p.y() - 0.7).abs() < 1e-6);
    assert!((p.z() - 0.21).abs() < 1e-6, "z = {}", p.z());
}

#[test]
fn transform_translates_all_poles() {
    let poles = vec![
        vec![GpPnt::new(0.0, 0.0, 0.0), GpPnt::new(0.0, 1.0, 0.0)],
        vec![GpPnt::new(1.0, 0.0, 0.0), GpPnt::new(1.0, 1.0, 0.0)],
    ];
    let (ku, kv) = bspline_surface_uniform_knots(2, 2, 1, 1);
    let mut s = GeomBSplineSurface::new(poles, ku, kv, 1, 1).unwrap();
    let before = s.d0(0.4, 0.6);
    let mut trsf = GpTrsf::identity();
    trsf.set_translation_vec(&GpVec::new(1.0, 2.0, 3.0));
    s.transform(&trsf);
    let after = s.d0(0.4, 0.6);
    assert!((after.x() - before.x() - 1.0).abs() < 1e-12);
    assert!((after.y() - before.y() - 2.0).abs() < 1e-12);
    assert!((after.z() - before.z() - 3.0).abs() < 1e-12);
}

#[test]
fn rational_weights_matches_weighted_formula() {
    // 2×2 degree-1 surface; at (0.5, 0.5) every basis product is 0.25.
    let poles = vec![
        vec![GpPnt::new(0.0, 0.0, 0.0), GpPnt::new(0.0, 1.0, 0.0)],
        vec![GpPnt::new(1.0, 0.0, 0.0), GpPnt::new(1.0, 1.0, 0.0)],
    ];
    let weights = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
    let (ku, kv) = bspline_surface_uniform_knots(2, 2, 1, 1);
    let s = GeomBSplineSurface::rational(poles, weights, ku, kv, 1, 1).unwrap();
    let p = s.d0(0.5, 0.5);
    // Σ w P = (1·0+2·0+3·1+4·1, 1·0+2·1+3·0+4·1, 0) = (7, 6, 0), Σ w = 10.
    assert!((p.x() - 0.7).abs() < 1e-12, "x = {}", p.x());
    assert!((p.y() - 0.6).abs() < 1e-12, "y = {}", p.y());
    assert!(p.z().abs() < 1e-12);
    // Non-rational version differs (weights pull the point toward heavy poles).
    let s_nr = GeomBSplineSurface::new(
        vec![
            vec![GpPnt::new(0.0, 0.0, 0.0), GpPnt::new(0.0, 1.0, 0.0)],
            vec![GpPnt::new(1.0, 0.0, 0.0), GpPnt::new(1.0, 1.0, 0.0)],
        ],
        bspline_surface_uniform_knots(2, 2, 1, 1).0,
        bspline_surface_uniform_knots(2, 2, 1, 1).1,
        1,
        1,
    )
    .unwrap();
    let p_nr = s_nr.d0(0.5, 0.5);
    assert!((p_nr.x() - 0.5).abs() < 1e-12);
}

#[test]
fn fit_surface_grid_passes_through_nodes() {
    // Cubic interpolation of z = u² + v³ over a 5×5 grid.
    let (nu, nv) = (5, 5);
    let points: Vec<Vec<GpPnt>> = (0..nu)
        .map(|i| {
            (0..nv)
                .map(|j| {
                    let u = i as f64 / (nu - 1) as f64;
                    let v = j as f64 / (nv - 1) as f64;
                    GpPnt::new(u, v, u * u + v * v * v)
                })
                .collect()
        })
        .collect();
    let s = fit_surface_grid(&points, 3, 3).unwrap();
    for i in 0..nu {
        for j in 0..nv {
            let u = i as f64 / (nu - 1) as f64;
            let v = j as f64 / (nv - 1) as f64;
            let p = s.d0(u, v);
            assert!(
                (p.z() - points[i][j].z()).abs() < 1e-6,
                "node ({i},{j}) z = {} vs {}",
                p.z(),
                points[i][j].z()
            );
        }
    }
}

#[test]
fn d1_nonzero_for_curved_surface() {
    let (nu, nv) = (5, 5);
    let points: Vec<Vec<GpPnt>> = (0..nu)
        .map(|i| {
            (0..nv)
                .map(|j| {
                    let u = i as f64 / (nu - 1) as f64;
                    let v = j as f64 / (nv - 1) as f64;
                    GpPnt::new(u, v, u * u * v)
                })
                .collect()
        })
        .collect();
    let s = fit_surface_grid(&points, 3, 3).unwrap();
    let (p, du, dv) = s.d1(0.5, 0.5);
    assert!((p.z() - 0.125).abs() < 1e-4, "z = {}", p.z());
    assert!(du.magnitude() > 1e-6, "du should be non-zero");
    assert!(dv.magnitude() > 1e-6, "dv should be non-zero");
}
