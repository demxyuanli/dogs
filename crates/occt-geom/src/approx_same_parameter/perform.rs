//! `Approx_SameParameter::Build` (`cxx:318-539`) and its accessors.

use super::*;

impl ApproxSameParameter {
    /// Adaptor constructor used by `BRepLib::SameParameter` (`cxx:1631`).
    pub fn new(
        c3d: &dyn Curve,
        c2d: &dyn Curve2d,
        surf: &dyn Surface,
        first: f64,
        last: f64,
        tol: f64,
    ) -> Self {
        let mut data = DistData {
            pc3d: [0.0; MAX_ARRAY],
            pc2d: [0.0; MAX_ARRAY],
            nb: 0,
            c3d_f: first,
            c3d_l: last,
            c2d_f: first,
            c2d_l: last,
            tol,
        };
        if !build_initial(&mut data, c2d) {
            return Self {
                done: false,
                same_parameter: false,
                tol_reached: compute_tol_reached(c3d, c2d, surf, first, last),
                curve2d: None,
            };
        }
        let keep_min = data.nb - ((0.3 * data.nb as f64) as usize);
        let mut sq = 0.0;
        let same = check_same_parameter(&mut data, c3d, c2d, surf, &mut sq);
        if same {
            return Self {
                done: true,
                same_parameter: true,
                tol_reached: compute_tol_reached(c3d, c2d, surf, first, last),
                curve2d: None,
            };
        }
        if data.nb < keep_min {
            return Self {
                done: false,
                same_parameter: false,
                tol_reached: compute_tol_reached(c3d, c2d, surf, first, last),
                curve2d: None,
            };
        }
        let Some((tang_first, tang_last)) = compute_tangents(c3d, c2d, surf) else {
            return Self {
                done: false,
                same_parameter: false,
                tol_reached: compute_tol_reached(c3d, c2d, surf, first, last),
                curve2d: None,
            };
        };
        let (done, tol_reached, curve2d) =
            rebuild_pcurve(&mut data, c3d, c2d, surf, tang_first, tang_last);
        Self {
            done,
            same_parameter: false,
            tol_reached,
            curve2d,
        }
    }
}
