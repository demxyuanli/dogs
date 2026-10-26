//! `IntWalk_PWalking` parameter helpers and distance minimization.
//!
//! Source: `IntWalk/IntWalk_PWalking.cxx`:
//! `IsParallel` (`:118-187`), `AdjustToDomain` (`:193-214`),
//! `DistanceMinimizeByGradient` (`:2394-2521`),
//! `DistanceMinimizeByExtrema` (`:2532-2582`), and `MaxStep`
//! (`IntWalk_PWalking.hxx:181-186`).

use occt_core::gp::{GpPnt, GpVec};
use occt_core::precision::{epsilon, INFINITE, PCONFUSION};
use occt_geom::Surface;

use crate::geom_bnd_lib_surface3d::{surface_kind, GeomAbsSurfaceType};

use super::pwalking::PWalking;

/// `aNbPointsMAX` (`IntWalk_PWalking.cxx:124`).
const NB_POINTS_MAX: i32 = 23;

/// `AdjustToDomain` (`IntWalk_PWalking.cxx:193-214`): returns TRUE if any
/// parameter has been clamped to the domain border.
pub(crate) fn adjust_to_domain(param: &mut [f64], low: &[f64], upp: &[f64]) -> bool {
    let mut ret = false;
    for i in 0..param.len() {
        if param[i] - low[i] < -PCONFUSION {
            param[i] = low[i];
            ret = true;
        }
        if param[i] - upp[i] > PCONFUSION {
            param[i] = upp[i];
            ret = true;
        }
    }
    ret
}

/// `RealToInt` (`Standard_Real.hxx`): C-cast truncation of a real to integer.
fn real_to_int(value: f64) -> i32 {
    value as i32
}

impl<'a> PWalking<'a> {
    /// `MaxStep` (`IntWalk_PWalking.hxx:181-186`): the initial step lengths
    /// (`pasInit`), not the current ones.
    pub fn max_step(&self, index: usize) -> f64 {
        assert!(index <= 3, "IntWalk_PWalking::MaxStep() - index is out of range");
        self.pas_init[index]
    }

    /// `IsParallel` (`IntWalk_PWalking.cxx:118-187`): checks whether the line
    /// runs parallel to a boundary of the surface selected by `check_surf1`
    /// (parallelism = swing less than `toler`).
    ///
    /// Returns `(is_u_parallel, is_v_parallel)`. A V-parallel line has a
    /// constant U (`theIsVparallel = (aUmax - aUmin) < theToler`), and lines
    /// shorter than three points are reported parallel on both axes.
    pub(crate) fn is_parallel(&self, check_surf1: bool, toler: f64) -> (bool, bool) {
        let mut is_u_parallel = true;
        let mut is_v_parallel = true;

        let nb_line_pnts = self.line.len() as i32;
        let mut nb_points = nb_line_pnts;
        if nb_points > NB_POINTS_MAX {
            nb_points = NB_POINTS_MAX;
        } else if nb_points < 3 {
            // Here we cannot estimate parallelism.
            return (is_u_parallel, is_v_parallel);
        }

        let step = nb_line_pnts as f64 / nb_points as f64;
        let mut a_n_point = 1.0;

        let (mut u_min, mut u_max) = (f64::MAX, -f64::MAX);
        let (mut v_min, mut v_max) = (f64::MAX, -f64::MAX);
        for _ in 1..=nb_points {
            // Fix possible "out of parameter" case.
            if a_n_point > nb_line_pnts as f64 {
                a_n_point = nb_line_pnts as f64;
            }
            let idx = real_to_int(a_n_point).max(1) as usize;
            let (u, v) = if check_surf1 {
                self.line[idx - 1].parameters_on_s1()
            } else {
                self.line[idx - 1].parameters_on_s2()
            };
            if u < u_min {
                u_min = u;
            }
            if u > u_max {
                u_max = u;
            }
            if v < v_min {
                v_min = v;
            }
            if v > v_max {
                v_max = v;
            }
            a_n_point += step;
        }

        is_v_parallel = (u_max - u_min) < toler;
        is_u_parallel = (v_max - v_min) < toler;
        (is_u_parallel, is_v_parallel)
    }

    /// `DistanceMinimizeByGradient` (`IntWalk_PWalking.cxx:2394-2521`).
    ///
    /// ATTENTION: `init` must be initialized before calling. Returns true when
    /// the surfaces are not Bezier / BSpline (nothing to minimize, OCCT answers
    /// `true` as well), otherwise whether the squared 3D distance fell below
    /// `1.0e-14`.
    pub(crate) fn distance_minimize_by_gradient(
        &self,
        init: &mut [f64; 4],
        step0: Option<[f64; 4]>,
    ) -> bool {
        const NB_ITER_MAX: i32 = 60;
        const TOL: f64 = 1.0e-14;
        // `1.0 / Precision::Infinite()`.
        let tol_nul = 1.0 / INFINITE;

        // I.e. if theU1 = 0.0 then Epsilon(theU1) = DBL_MIN (~1.0e-308).
        // Work with this number is impossible: there is a dangerous to
        // obtain Floating-point-overflow. Therefore, we limit this value.
        let min_add = [
            epsilon(init[0]).max(tol_nul),
            epsilon(init[1]).max(tol_nul),
            epsilon(init[2]).max(tol_nul),
            epsilon(init[3]).max(tol_nul),
        ];

        let kind1 = surface_kind(self.s1);
        if kind1 != GeomAbsSurfaceType::BezierSurface
            && kind1 != GeomAbsSurfaceType::BSplineSurface
        {
            return true;
        }
        let kind2 = surface_kind(self.s2);
        if kind2 != GeomAbsSurfaceType::BezierSurface
            && kind2 != GeomAbsSurfaceType::BSplineSurface
        {
            return true;
        }

        let mut status = false;

        let (p1, d1u, d1v) = self.s1.d1(init[0], init[1]);
        let (p2, d2u, d2v) = self.s2.d1(init[2], init[3]);

        let mut sq_dist_prev = p1.square_distance(&p2);

        let p12 = GpVec::from_pnts(&p1, &p2);

        let mut grad_fu = -p12.dot(&d1u);
        let mut grad_fv = -p12.dot(&d1v);
        let mut grad_fu2 = p12.dot(&d2u);
        let mut grad_fv2 = p12.dot(&d2v);

        let mut steps = step0.unwrap_or([1.0e-6; 4]);

        let mut repeat = true;
        let mut nb_iter = NB_ITER_MAX;

        while repeat {
            let mut add = grad_fu * steps[0];
            let par_u = init[0] - add.abs().max(min_add[0]).copysign(add);
            add = grad_fv * steps[1];
            let par_v = init[1] - add.abs().max(min_add[1]).copysign(add);
            add = grad_fu2 * steps[2];
            let par_u2 = init[2] - add.abs().max(min_add[2]).copysign(add);
            add = grad_fv2 * steps[3];
            let par_v2 = init[3] - add.abs().max(min_add[3]).copysign(add);

            let (pt1, _du1, _dv1) = self.s1.d1(par_u, par_v);
            let (pt2, _du2, _dv2) = self.s2.d1(par_u2, par_v2);

            let sq_dist = pt1.square_distance(&pt2);

            if sq_dist < sq_dist_prev {
                sq_dist_prev = sq_dist;
                init[0] = par_u;
                init[1] = par_v;
                init[2] = par_u2;
                init[3] = par_v2;

                status = sq_dist_prev < TOL;
                steps[0] *= 1.2;
                steps[1] *= 1.2;
                steps[2] *= 1.2;
                steps[3] *= 1.2;
            } else {
                nb_iter -= 1;
                if nb_iter < 0 {
                    repeat = false;
                } else {
                    let (pt1, du1, dv1) = self.s1.d1(init[0], init[1]);
                    let (pt2, du2, dv2) = self.s2.d1(init[2], init[3]);

                    let p12 = GpVec::from_pnts(&pt1, &pt2);
                    grad_fu = -p12.dot(&du1);
                    grad_fv = -p12.dot(&dv1);
                    grad_fu2 = p12.dot(&du2);
                    grad_fv2 = p12.dot(&dv2);

                    match step0 {
                        Some(s) => steps = s,
                        None => steps = [1.0e-6; 4],
                    }
                }
            }
        }

        status
    }

    /// `DistanceMinimizeByExtrema` (`IntWalk_PWalking.cxx:2532-2582`).
    ///
    /// ATTENTION: `p0`, `u0` and `v0` must be initialized before calling.
    pub(crate) fn distance_minimize_by_extrema(
        &self,
        surf: &dyn Surface,
        p0: &GpPnt,
        uv: &mut [f64],
        step0: Option<[f64; 2]>,
    ) -> bool {
        const TOL: f64 = 1.0e-14;
        let mut sq_dist_prev = f64::MAX;
        let mut u = uv[0];
        let mut v = uv[1];

        let step = step0.unwrap_or([1.0, 1.0]);

        let mut nb_iter = 10;
        loop {
            let (ps, d1su, d1sv, d2su, d2sv, _d2suv) = surf.d2(u, v);

            let vec = GpVec::from_pnts(p0, &ps);

            let sq_dist = vec.dot(&vec);

            if sq_dist >= sq_dist_prev {
                break;
            }

            sq_dist_prev = sq_dist;
            uv[0] = u;
            uv[1] = v;
            nb_iter -= 1;

            if sq_dist_prev < TOL {
                break;
            }

            // Functions
            let f1 = d1su.dot(&vec);
            let f2 = d1sv.dot(&vec);

            // Derivatives
            let df1u = d2su.dot(&vec) + d1su.dot(&d1su);
            let df1v = d2su.dot(&d1sv);
            let df2u = df1v;
            let df2v = d2sv.dot(&vec) + d1sv.dot(&d1sv);

            let det = df1u * df2v - df1v * df2u;
            u -= step[0] * (df2v * f1 - df1v * f2) / det;
            v += step[1] * (df2u * f1 - df1u * f2) / det;

            if nb_iter <= 0 {
                break;
            }
        }

        sq_dist_prev < TOL
    }
}
