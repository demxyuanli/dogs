//! BSpline conversion depth. Port of the `GeomConvert` gap-fillers
//! (TKGeomBase): `GeomConvert_BSplineCurveKnotSplitting`,
//! `GeomConvert_BSplineSurfaceKnotSplitting`, `GeomConvert_BSplineSurfaceToBezierSurface`
//! and `GeomConvert_CompCurveToBSplineCurve`.
//!
//! Curve→Bezier is already covered by `crate::bspline_to_bezier` and is not
//! duplicated here; `GeomConvert_ApproxCurve` (adaptive spline approximation)
//! is covered pragmatically by `curve_approx` (polyline) plus
//! `curve_reparam::resample_bspline` (interpolation) and is deferred.

use occt_core::bspl::bezier::boehm_insert;
use occt_core::bspl::knots::{hunt, insert_knot, multiplicity};
use occt_core::gp::GpPnt;

use crate::bspline_curve::GeomBSplineCurve;
use crate::bspline_surface::GeomBSplineSurface;
use crate::curve::Curve;

// ---------------------------------------------------------------------------
// Knot splitting (GeomConvert_BSplineCurveKnotSplitting /
// BSplineSurfaceKnotSplitting).
// ---------------------------------------------------------------------------

/// Distinct knot values and their multiplicities from a full knot vector.
fn distinct_knots(knots: &[f64]) -> (Vec<f64>, Vec<usize>) {
    let mut vals: Vec<f64> = Vec::new();
    let mut mults: Vec<usize> = Vec::new();
    for &k in knots {
        if let Some(last) = vals.last() {
            if (k - last).abs() < 1e-15 {
                *mults.last_mut().unwrap() += 1;
                continue;
            }
        }
        vals.push(k);
        mults.push(1);
    }
    (vals, mults)
}

/// Split indices for `degree`/`mults` at the requested continuity range.
///
/// Returns 0-based indices into the distinct-knot array (the OCCT
/// `splitIndexes`), starting at the first and ending at the last distinct
/// knot. Port of `GeomConvert_BSplineCurveKnotSplitting` / the surface
/// counterpart.
fn knot_splits(degree: usize, mults: &[usize], continuity_range: usize) -> Vec<usize> {
    let first = 0usize;
    let last = mults.len() - 1;
    if continuity_range == 0 {
        return vec![first, last];
    }
    let mmax = *mults.iter().max().unwrap_or(&0);
    if degree >= mmax + continuity_range {
        return vec![first, last];
    }
    let mut out = vec![first];
    for i in 1..last {
        // Continuity at knot i is (degree − multiplicity); split where it
        // drops below `continuity_range`.
        if degree < mults[i] + continuity_range {
            out.push(i);
        }
    }
    out.push(last);
    out
}

/// Split indices of `c` where its continuity drops below C1 (the OCCT default
/// `ContinuityRange = 1`). Indices are into the curve's distinct-knot array.
pub fn knot_splitting_curve(c: &GeomBSplineCurve) -> Vec<usize> {
    knot_splitting_curve_with_continuity(c, 1)
}

/// `knot_splitting_curve` with an explicit continuity range.
pub fn knot_splitting_curve_with_continuity(c: &GeomBSplineCurve, continuity_range: usize) -> Vec<usize> {
    let (_, mults) = distinct_knots(&c.knots);
    knot_splits(c.degree, &mults, continuity_range)
}

/// Split indices of `s` in u and v (each a `Vec<usize>` into the respective
/// distinct-knot array), at the OCCT default continuity range 1.
pub fn knot_splitting_surface(s: &GeomBSplineSurface) -> (Vec<usize>, Vec<usize>) {
    knot_splitting_surface_with_continuity(s, 1, 1)
}

/// `knot_splitting_surface` with explicit u/v continuity ranges.
pub fn knot_splitting_surface_with_continuity(
    s: &GeomBSplineSurface,
    u_continuity: usize,
    v_continuity: usize,
) -> (Vec<usize>, Vec<usize>) {
    let (_, um) = distinct_knots(&s.knots_u);
    let (_, vm) = distinct_knots(&s.knots_v);
    (
        knot_splits(s.deg_u, &um, u_continuity),
        knot_splits(s.deg_v, &vm, v_continuity),
    )
}

// ---------------------------------------------------------------------------
// BSpline surface → Bezier surface patches
// (GeomConvert_BSplineSurfaceToBezierSurface).
// ---------------------------------------------------------------------------

/// A single Bezier patch of a split B-spline surface, over the parameter
/// rectangle `[u0, u1] × [v0, v1]` with a `(deg_u+1) × (deg_v+1)` pole grid.
#[derive(Debug, Clone)]
pub struct BezierSurfacePatch {
    pub poles: Vec<Vec<GpPnt>>,
    pub weights: Option<Vec<Vec<f64>>>,
    pub u0: f64,
    pub u1: f64,
    pub v0: f64,
    pub v1: f64,
}

fn insert_surface_u(s: &GeomBSplineSurface, u: f64, mult: usize) -> GeomBSplineSurface {
    let nu = s.poles.len();
    let nv = s.poles[0].len();
    let deg = s.deg_u;
    let new_nu = nu + mult;
    let mut new_poles = vec![vec![GpPnt::zero(); nv]; new_nu];
    let mut new_weights = s.weights.as_ref().map(|_| vec![vec![0.0; nv]; new_nu]);
    for j in 0..nv {
        let mut col: Vec<GpPnt> = (0..nu).map(|i| s.poles[i][j]).collect();
        let mut wcol = s.weights.as_ref().map(|w| (0..nu).map(|i| w[i][j]).collect::<Vec<_>>());
        let idx = hunt(&s.knots_u, u);
        for _ in 0..mult {
            boehm_insert(&mut col, &s.knots_u, idx, u, deg, wcol.as_mut());
        }
        for i in 0..new_nu {
            new_poles[i][j] = col[i];
            if let (Some(nw), Some(wc)) = (new_weights.as_mut(), wcol.as_ref()) {
                nw[i][j] = wc[i];
            }
        }
    }
    let mut knots_u = s.knots_u.clone();
    for _ in 0..mult {
        knots_u = insert_knot(&knots_u, u, 1);
    }
    GeomBSplineSurface {
        poles: new_poles,
        knots_u,
        knots_v: s.knots_v.clone(),
        deg_u: s.deg_u,
        deg_v: s.deg_v,
        weights: new_weights,
        u_periodic: s.u_periodic,
        v_periodic: s.v_periodic,
    }
}

fn insert_surface_v(s: &GeomBSplineSurface, v: f64, mult: usize) -> GeomBSplineSurface {
    let nu = s.poles.len();
    let nv = s.poles[0].len();
    let deg = s.deg_v;
    let new_nv = nv + mult;
    let mut new_poles = vec![vec![GpPnt::zero(); new_nv]; nu];
    let mut new_weights = s.weights.as_ref().map(|_| vec![vec![0.0; new_nv]; nu]);
    for i in 0..nu {
        let mut row: Vec<GpPnt> = s.poles[i].clone();
        let mut wrow = s.weights.as_ref().map(|w| w[i].clone());
        let idx = hunt(&s.knots_v, v);
        for _ in 0..mult {
            boehm_insert(&mut row, &s.knots_v, idx, v, deg, wrow.as_mut());
        }
        for j in 0..new_nv {
            new_poles[i][j] = row[j];
            if let (Some(nw), Some(wr)) = (new_weights.as_mut(), wrow.as_ref()) {
                nw[i][j] = wr[j];
            }
        }
    }
    let mut knots_v = s.knots_v.clone();
    for _ in 0..mult {
        knots_v = insert_knot(&knots_v, v, 1);
    }
    GeomBSplineSurface {
        poles: new_poles,
        knots_u: s.knots_u.clone(),
        knots_v,
        deg_u: s.deg_u,
        deg_v: s.deg_v,
        weights: new_weights,
        u_periodic: s.u_periodic,
        v_periodic: s.v_periodic,
    }
}

fn distinct_in_range(knots: &[f64], first: f64, last: f64) -> Vec<f64> {
    let mut out: Vec<f64> = Vec::new();
    for &k in knots {
        if k < first - 1e-15 || k > last + 1e-15 {
            continue;
        }
        if out.last().map_or(true, |&b| (b - k).abs() > 1e-15) {
            out.push(k);
        }
    }
    out
}

/// Convert a clamped B-spline surface into its Bezier patches by knot
/// insertion: every interior knot is raised to multiplicity `degree`, after
/// which each span is a Bezier patch. Mirrors
/// `GeomConvert_BSplineSurfaceToBezierSurface`.
pub fn bspline_surface_to_bezier_surface(s: &GeomBSplineSurface) -> Vec<BezierSurfacePatch> {
    let mut s = s.clone();
    let du = s.deg_u;
    let dv = s.deg_v;
    let (fu, lu) = (s.knots_u[du], s.knots_u[s.knots_u.len() - 1 - du]);
    let (fv, lv) = (s.knots_v[dv], s.knots_v[s.knots_v.len() - 1 - dv]);

    let interior_u: Vec<f64> = s
        .knots_u
        .iter()
        .copied()
        .filter(|&k| k > fu + 1e-15 && k < lu - 1e-15)
        .collect::<Vec<_>>()
        .into_iter()
        .fold(Vec::new(), |mut acc, k| {
            if acc.last().map_or(true, |&l| (l - k).abs() > 1e-15) {
                acc.push(k);
            }
            acc
        });
    for u in interior_u {
        let m = multiplicity(&s.knots_u, u);
        if m < du {
            s = insert_surface_u(&s, u, du - m);
        }
    }
    let interior_v: Vec<f64> = s
        .knots_v
        .iter()
        .copied()
        .filter(|&k| k > fv + 1e-15 && k < lv - 1e-15)
        .collect::<Vec<_>>()
        .into_iter()
        .fold(Vec::new(), |mut acc, k| {
            if acc.last().map_or(true, |&l| (l - k).abs() > 1e-15) {
                acc.push(k);
            }
            acc
        });
    for v in interior_v {
        let m = multiplicity(&s.knots_v, v);
        if m < dv {
            s = insert_surface_v(&s, v, dv - m);
        }
    }

    let fu = s.knots_u[du];
    let lu = s.knots_u[s.knots_u.len() - 1 - du];
    let fv = s.knots_v[dv];
    let lv = s.knots_v[s.knots_v.len() - 1 - dv];
    let bpu = distinct_in_range(&s.knots_u, fu, lu);
    let bpv = distinct_in_range(&s.knots_v, fv, lv);
    if bpu.len() < 2 || bpv.len() < 2 {
        return Vec::new();
    }

    let mut patches = Vec::with_capacity((bpu.len() - 1) * (bpv.len() - 1));
    for iu in 0..bpu.len() - 1 {
        for iv in 0..bpv.len() - 1 {
            let a_u = iu * du;
            let a_v = iv * dv;
            let mut poles = Vec::with_capacity(du + 1);
            let mut weights = s.weights.as_ref().map(|_| Vec::with_capacity(du + 1));
            for i in a_u..=a_u + du {
                poles.push(s.poles[i][a_v..=a_v + dv].to_vec());
                if let (Some(w), Some(sw)) = (weights.as_mut(), s.weights.as_ref()) {
                    w.push(sw[i][a_v..=a_v + dv].to_vec());
                }
            }
            patches.push(BezierSurfacePatch {
                poles,
                weights,
                u0: bpu[iu],
                u1: bpu[iu + 1],
                v0: bpv[iv],
                v1: bpv[iv + 1],
            });
        }
    }
    patches
}

// ---------------------------------------------------------------------------
// CompCurveToBSplineCurve.
// ---------------------------------------------------------------------------

/// Merge connected curve segments into a single degree-1 B-spline passing
/// through the sampled composite polyline (port of
/// `GeomConvert_CompCurveToBSplineCurve`; degree 1 reproduces polygonal
/// composites exactly, which is the OCCT use case for line/arc chains).
///
/// Each segment must be bounded and connected to the next within `tol`
/// (orientation is detected automatically).
pub fn comp_curve_to_bspline(segs: &[Box<dyn Curve>], tol: f64) -> Result<GeomBSplineCurve, String> {
    if segs.is_empty() {
        return Err("comp_curve_to_bspline: no segments".to_string());
    }
    let tol = tol.max(1e-9);
    let mut pts: Vec<GpPnt> = Vec::new();
    let mut cur: Option<GpPnt> = None;
    for seg in segs {
        let (a, b) = (seg.first_parameter(), seg.last_parameter());
        if !(a.is_finite() && b.is_finite()) || b <= a {
            return Err("comp_curve_to_bspline: unbounded or degenerate segment".to_string());
        }
        let (start, end, reverse) = match cur {
            None => (a, b, false),
            Some(c) => {
                if seg.d0(a).distance(&c) <= tol {
                    (a, b, false)
                } else if seg.d0(b).distance(&c) <= tol {
                    (b, a, true)
                } else {
                    return Err("comp_curve_to_bspline: segments not connected".to_string());
                }
            }
        };
        let n = ((b - a) / 0.1).ceil().clamp(2.0, 64.0) as usize;
        for k in 0..=n {
            let t = if reverse {
                1.0 - k as f64 / n as f64
            } else {
                k as f64 / n as f64
            };
            let p = seg.d0(start + (end - start) * t);
            if pts.last().map_or(true, |l| l.distance(&p) > tol) {
                pts.push(p);
            }
        }
        cur = Some(*pts.last().unwrap());
    }
    // Ensure the final endpoint is included.
    let last_seg = segs.last().unwrap();
    let end = last_seg.d0(last_seg.last_parameter());
    if pts.last().map_or(true, |p| p.distance(&end) > tol) {
        pts.push(end);
    }
    if pts.len() < 2 {
        return Err("comp_curve_to_bspline: too few points".to_string());
    }
    crate::curve_reparam::resample_bspline(&pts, 1).map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::surface::Surface;
    use crate::{GeomLine, GeomTrimmedCurve};
    use occt_core::gp::{GpDir, GpPnt};
    use std::sync::Arc;

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

    #[test]
    fn comp_curve_merges_two_lines() {
        // Two connected line segments forming an L.
        let l1 = GeomLine::from_pnt_dir(GpPnt::new(0., 0., 0.), GpDir::new(1., 0., 0.).unwrap());
        let t1 = GeomTrimmedCurve::new(Arc::new(l1), 0.0, 1.0);
        let l2 = GeomLine::from_pnt_dir(GpPnt::new(1., 0., 0.), GpDir::new(0., 1., 0.).unwrap());
        let t2 = GeomTrimmedCurve::new(Arc::new(l2), 0.0, 1.0);

        let segs: Vec<Box<dyn Curve>> = vec![Box::new(t1), Box::new(t2)];
        let spline = comp_curve_to_bspline(&segs, 1e-7).expect("merge");
        // The degree-1 spline interpolates its sampled polyline: endpoints and
        // the junction are exact.
        assert!(spline.d0(spline.first_parameter()).distance(&GpPnt::new(0., 0., 0.)) < 1e-9);
        assert!(spline.d0(spline.last_parameter()).distance(&GpPnt::new(1., 1., 0.)) < 1e-9);
        // Junction is at arc length 1.0 of the total 2.0 → parameter 0.5.
        let mid = spline.d0(0.5);
        assert!(mid.distance(&GpPnt::new(1., 0., 0.)) < 1e-6, "junction {mid:?}");
    }

    #[test]
    fn comp_curve_rejects_disconnected() {
        let l1 = GeomLine::from_pnt_dir(GpPnt::new(0., 0., 0.), GpDir::new(1., 0., 0.).unwrap());
        let t1 = GeomTrimmedCurve::new(Arc::new(l1), 0.0, 1.0);
        let l2 = GeomLine::from_pnt_dir(GpPnt::new(5., 0., 0.), GpDir::new(0., 1., 0.).unwrap());
        let t2 = GeomTrimmedCurve::new(Arc::new(l2), 0.0, 1.0);
        let segs: Vec<Box<dyn Curve>> = vec![Box::new(t1), Box::new(t2)];
        assert!(comp_curve_to_bspline(&segs, 1e-7).is_err());
    }
}
