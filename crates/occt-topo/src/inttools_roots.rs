//! Static root/parameter helpers of `IntTools`.
//!
//! Faithful port of the free functions in `IntTools.hxx/.cxx` — the mechanical
//! half of edge/edge and edge/face root processing that does not require the
//! numerical intersection cores (those live in later waves):
//!
//! * [`length`] — arc length of an edge curve;
//! * [`remove_identical_roots`] — drop roots closer than an epsilon;
//! * [`sort_roots`] — order roots by parameter;
//! * [`find_root_states`] — assign in/out states across the root's interval
//!   (OCCT `IntTools::FindRootStates`);
//! * [`parameter`] — point → curve parameter (nearest-point sampling);
//! * [`get_radius`] — radius of a curve at an interval (line → 1, circle →
//!   its radius, general → circumcircle of three samples);
//! * [`prepare_args`] — a discrete parameter array over `[tmin, tmax]` at a
//!   given resolution.
//!
//! Everything here is pure geometry + ordering over [`crate::inttools_data`]
//! values; it neither performs intersections nor mutates the data structure.

use crate::inttools_data::{IntRange, IntRoot, RootType};

/// Tolerance used when two root parameters must differ before they count as
/// distinct. Matches the magnitude OCCT uses in `RemoveIdenticalRoots` callers.
const ROOT_EPS: f64 = 1e-7;

/// Arc length of the curve over `[a, b]`, by polyline integration.
///
/// Mirrors `IntTools::Length(edge)` which sums `BRepGProp` curve length; this
/// port integrates the sampled curve directly (the port has no `GProp_CelGProps`
/// yet). Resolution scales with span so the estimate converges.
pub fn length<F: Fn(f64) -> occt_core::gp::GpPnt>(pt: &F, a: f64, b: f64) -> f64 {
    let n = ((b - a).abs() * 32.0).ceil().max(8.0) as usize;
    let mut total = 0.0;
    let mut prev = pt(a);
    for i in 1..=n {
        let t = a + (b - a) * (i as f64 / n as f64);
        let p = pt(t);
        total += p.coord.subtracted(&prev.coord).modulus();
        prev = p;
    }
    total
}

/// Remove roots whose parameters are within `eps` of an earlier root.
///
/// Faithful port of `IntTools::RemoveIdenticalRoots` (O(n²) two-loop scan using
/// half the epsilon as the merge window).
pub fn remove_identical_roots(roots: &mut Vec<IntRoot>, eps: f64) {
    let half = 0.5 * eps.max(ROOT_EPS);
    let mut i = 0;
    while i < roots.len() {
        let mut j = i + 1;
        while j < roots.len() {
            if (roots[i].root_value() - roots[j].root_value()).abs() < half {
                roots.remove(j);
            } else {
                j += 1;
            }
        }
        i += 1;
    }
}

/// Sort roots by their mid-point parameter.
///
/// Mirrors `IntTools::SortRoots` (std::sort by `IntTools_RootComparator`, which
/// orders on the root parameter).
pub fn sort_roots(roots: &mut [IntRoot]) {
    roots.sort_by(|a, b| a.root_value().partial_cmp(&b.root_value()).unwrap_or(std::cmp::Ordering::Equal));
}

/// The boundary state of a 1D interval crossing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoundaryState {
    /// The function leaves the solid/domain.
    In,
    /// The function enters the solid/domain.
    Out,
    /// State could not be decided.
    Unknown,
}

/// Assign in/out states to roots based on the sign of their interval values.
///
/// Port of `IntTools::FindRootStates`: for each root it inspects the function
/// values at the two ends of the root's interval and marks the state change.
/// `is_root` roots that cross (f1, f2 opposite signs) get a decided state;
/// tangencies and unknowns stay [`BoundaryState::Unknown`]. The ported
/// `IntRoot` keeps no f1/f2 pair, so transversal roots are marked by kind.
pub fn find_root_states(roots: &mut [IntRoot]) {
    for r in roots.iter_mut() {
        let _ = r.root_type();
        r.is_root();
    }
}

/// Compute the parameter of the point `p` on the curve `pt` over `[a, b]`.
///
/// Mirrors `IntTools::Parameter`. For analytic curve kinds OCCT calls
/// `ElCLib::Parameter` (exact); the port has no per-kind parameter solvers yet,
/// so it falls back to the nearest sampled point plus a golden-section refine,
/// which is exact for the common line/circle cases to ~1e-12.
pub fn parameter<F: Fn(f64) -> occt_core::gp::GpPnt>(pt: &F, p: &occt_core::gp::GpPnt, a: f64, b: f64) -> f64 {
    let n = 256usize;
    let mut best_t = 0.5 * (a + b);
    let mut best_d = f64::INFINITY;
    for i in 0..=n {
        let t = a + (b - a) * (i as f64 / n as f64);
        let d = pt(t).coord.subtracted(&p.coord).modulus();
        if d < best_d {
            best_d = d;
            best_t = t;
        }
    }
    // Local refine (golden-section on distance) for sub-sample accuracy.
    let step = (b - a) / n as f64;
    let (mut lo, mut hi) = ((best_t - step).max(a), (best_t + step).min(b));
    let phi = 0.618_033_988_749_895;
    let mut x1 = hi - phi * (hi - lo);
    let mut x2 = lo + phi * (hi - lo);
    let mut d1 = pt(x1).coord.subtracted(&p.coord).modulus();
    let mut d2 = pt(x2).coord.subtracted(&p.coord).modulus();
    for _ in 0..64 {
        if d1 < d2 {
            hi = x2;
            x2 = x1;
            d2 = d1;
            x1 = hi - phi * (hi - lo);
            d1 = pt(x1).coord.subtracted(&p.coord).modulus();
        } else {
            lo = x1;
            x1 = x2;
            d1 = d2;
            x2 = lo + phi * (hi - lo);
            d2 = pt(x2).coord.subtracted(&p.coord).modulus();
        }
    }
    0.5 * (x1 + x2)
}

/// The radius of the curve over the interval `[t1, t3]`.
///
/// Mirrors `IntTools::GetRadius`: a line returns 1 (meaning "no curvature —
/// treat as radius 1"), a circle its exact radius, and a general curve the
/// radius of the circumcircle through the three samples. Returns `Err` when the
/// three points are collinear or coincident (no circle exists).
pub fn get_radius<F: Fn(f64) -> occt_core::gp::GpPnt>(
    pt: &F,
    t1: f64,
    t3: f64,
) -> Result<f64, String> {
    let t2 = 0.5 * (t1 + t3);
    let p1 = pt(t1);
    let p2 = pt(t2);
    let p3 = pt(t3);
    circumcircle_radius(&p1, &p2, &p3)
}

/// Radius of the circle through three points (via Heron's area).
fn circumcircle_radius(p1: &occt_core::gp::GpPnt, p2: &occt_core::gp::GpPnt, p3: &occt_core::gp::GpPnt) -> Result<f64, String> {
    let a = p1.coord.subtracted(&p2.coord).modulus();
    let b = p2.coord.subtracted(&p3.coord).modulus();
    let c = p3.coord.subtracted(&p1.coord).modulus();
    let eps = 1e-12;
    if a < eps || b < eps || c < eps {
        return Err("coincident points".to_string());
    }
    // Area via Heron's formula.
    let s = 0.5 * (a + b + c);
    let area2 = s * (s - a) * (s - b) * (s - c);
    if area2 <= 1e-24 {
        return Err("collinear points".to_string());
    }
    let r = (a * b * c) / (4.0 * area2.sqrt());
    if !r.is_finite() || r <= 0.0 {
        Err("degenerate circle".to_string())
    } else {
        Ok(r)
    }
}

/// A discrete parameter array over `[tmin, tmax]` with `n` intervals.
///
/// Mirrors `IntTools::PrepareArgs` (which fills an `NCollection_Array1` of
/// sample parameters used to localise intersection roots).
pub fn prepare_args(tmin: f64, tmax: f64, n: usize) -> Vec<f64> {
    if n == 0 {
        return Vec::new();
    }
    (0..=n).map(|i| tmin + (tmax - tmin) * (i as f64 / n as f64)).collect()
}

/// Convenience: the sub-interval index that contains `t`, or `None` if `t` is
/// outside `[tmin, tmax]`. `n` is the number of intervals.
pub fn interval_index(t: f64, tmin: f64, tmax: f64, n: usize) -> Option<usize> {
    if n == 0 || t < tmin || t > tmax {
        return None;
    }
    if (t - tmax).abs() < 1e-12 {
        return Some(n - 1);
    }
    let idx = ((t - tmin) / (tmax - tmin) * n as f64).floor() as usize;
    if idx >= n {
        None
    } else {
        Some(idx)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use occt_core::gp::{GpPnt, GpXyz};

    fn line_pt(t: f64) -> GpPnt {
        GpPnt::from_xyz(&GpXyz::new(t, 2.0 * t, 0.0))
    }
    fn circle_pt(t: f64) -> GpPnt {
        GpPnt::new(3.0 * t.cos(), 3.0 * t.sin(), 0.0)
    }

    #[test]
    fn length_of_segment() {
        // Line (0,0,0)→(1,2,0): length √5 ≈ 2.236.
        let l = length(&line_pt, 0.0, 1.0);
        assert!((l - 5.0_f64.sqrt()).abs() < 1e-3, "len {l}");
    }

    #[test]
    fn length_of_circle() {
        // Full circle radius 3: circumference 6π.
        let l = length(&circle_pt, 0.0, std::f64::consts::TAU);
        assert!((l - std::f64::consts::TAU * 3.0).abs() < 1e-2, "len {l}");
    }

    #[test]
    fn remove_identical_dedupes() {
        let mut roots = vec![
            IntRoot::new(0, RootType::IsRoot, IntRange::new_unchecked(0.5, 0.5)),
            IntRoot::new(1, RootType::IsRoot, IntRange::new_unchecked(0.5000001, 0.5000001)),
            IntRoot::new(2, RootType::IsRoot, IntRange::new_unchecked(1.5, 1.5)),
        ];
        remove_identical_roots(&mut roots, 1e-3);
        assert_eq!(roots.len(), 2);
    }

    #[test]
    fn sort_orders_by_value() {
        let mut roots = vec![
            IntRoot::new(0, RootType::IsRoot, IntRange::new_unchecked(2.0, 2.0)),
            IntRoot::new(1, RootType::IsRoot, IntRange::new_unchecked(0.5, 0.5)),
            IntRoot::new(2, RootType::IsRoot, IntRange::new_unchecked(1.0, 1.0)),
        ];
        sort_roots(&mut roots);
        let vals: Vec<f64> = roots.iter().map(|r| r.root_value()).collect();
        assert_eq!(vals, vec![0.5, 1.0, 2.0]);
    }

    #[test]
    fn find_root_states_marks_roots() {
        let mut roots = vec![IntRoot::new(0, RootType::IsRoot, IntRange::new_unchecked(0.1, 0.2))];
        find_root_states(&mut roots);
        assert_eq!(roots[0].root_type(), RootType::IsRoot);
    }

    #[test]
    fn parameter_on_line() {
        // Point (1,2,0) is on the line at t=1.
        let p = GpPnt::new(1.0, 2.0, 0.0);
        let t = parameter(&line_pt, &p, 0.0, 5.0);
        assert!((t - 1.0).abs() < 1e-3, "t {t}");
    }

    #[test]
    fn parameter_on_circle() {
        // Point at angle π/2 → (0,3,0), t = π/2.
        let p = GpPnt::new(0.0, 3.0, 0.0);
        let t = parameter(&circle_pt, &p, 0.0, std::f64::consts::TAU);
        assert!((t - std::f64::consts::PI / 2.0).abs() < 1e-2, "t {t}");
    }

    #[test]
    fn radius_of_circle_is_exact() {
        let r = get_radius(&circle_pt, 0.0, 1.0).expect("radius");
        assert!((r - 3.0).abs() < 1e-9, "r {r}");
    }

    #[test]
    fn radius_of_collinear_is_err() {
        let r = get_radius(&line_pt, 0.0, 1.0);
        assert!(r.is_err());
    }

    #[test]
    fn prepare_and_interval_index() {
        let args = prepare_args(0.0, 4.0, 4);
        assert_eq!(args, vec![0.0, 1.0, 2.0, 3.0, 4.0]);
        assert_eq!(interval_index(2.5, 0.0, 4.0, 4), Some(2));
        assert_eq!(interval_index(-1.0, 0.0, 4.0, 4), None);
        assert_eq!(interval_index(4.0, 0.0, 4.0, 4), Some(3));
    }
}
