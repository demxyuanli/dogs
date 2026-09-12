//! `Approx_SameParameter` check arm.
//! Source: `Approx_SameParameter.cxx` (`Build` 318-361, `BuildInitialDistribution`
//! 547-579, `CheckSameParameter` 653-766, `ComputeTolReached` 153-188,
//! `ProjectPointOnCurve` 106-149).
//!
//! Interpolation + `AdvApprox` 2D rebuild (`cxx:389-539`) is not ported.

use occt_core::gp::{GpPnt, GpVec};
use occt_core::precision::{CONFUSION, INFINITE, PCONFUSION};

use crate::curve::Curve;
use crate::surface::Surface;
use occt_geom2d::curve::Curve2d;

const NB_SAMPLES: usize = 22;
const MAX_ARRAY: usize = 1000;

/// Result of `Approx_SameParameter`.
pub struct ApproxSameParameter {
    pub done: bool,
    pub same_parameter: bool,
    pub tol_reached: f64,
}

struct DistData {
    pc3d: [f64; MAX_ARRAY],
    pc2d: [f64; MAX_ARRAY],
    nb: usize,
    c3d_f: f64,
    c3d_l: f64,
    c2d_f: f64,
    c2d_l: f64,
    tol: f64,
}

fn cons_value(c2d: &dyn Curve2d, s: &dyn Surface, t: f64) -> GpPnt {
    let uv = c2d.d0(t);
    s.d0(uv.x(), uv.y())
}

fn compute_tol_reached(c3d: &dyn Curve, c2d: &dyn Curve2d, s: &dyn Surface, first: f64, last: f64) -> f64 {
    let nbp = 2 * NB_SAMPLES;
    let mut d2: f64 = 0.0;
    for i in 0..=nbp {
        let t = i as f64 / nbp as f64;
        let u = first * (1.0 - t) + last * t;
        let pc3d = c3d.d0(u);
        let pcons = cons_value(c2d, s, u);
        if !pcons.x().is_finite() || !pcons.y().is_finite() || !pcons.z().is_finite() {
            return INFINITE;
        }
        d2 = d2.max(pc3d.square_distance(&pcons));
    }
    (1.05 * d2.sqrt()).max(CONFUSION)
}

fn project_point_on_curve(
    init: f64,
    point: &GpPnt,
    tolerance: f64,
    curve: &dyn Curve,
    first: f64,
    last: f64,
) -> Option<f64> {
    let mut param = init;
    for _ in 0..30 {
        let (a_point, d1, d2) = curve.d2(param);
        let vector = GpVec::from_pnts(&a_point, point);
        let func = vector.dot(&d1);
        if func.abs() < tolerance * d1.magnitude() {
            return Some(param);
        }
        let func_derivative = vector.dot(&d2) - d1.dot(&d1);
        if func_derivative.abs() > 1.0e-12 {
            param -= func / func_derivative;
        }
        param = param.max(first).min(last);
    }
    None
}

fn build_initial(data: &mut DistData, c2d: &dyn Curve2d) -> bool {
    let deltacons = (data.c2d_l - data.c2d_f) / NB_SAMPLES as f64;
    let deltac3d = (data.c3d_l - data.c3d_f) / NB_SAMPLES as f64;
    let mut wcons = data.c2d_f;
    let mut wc3d = data.c3d_f;
    for ii in 0..NB_SAMPLES {
        data.pc2d[ii] = wcons;
        data.pc3d[ii] = wc3d;
        wcons += deltacons;
        wc3d += deltac3d;
    }
    data.nb = NB_SAMPLES;
    data.pc2d[data.nb] = data.c2d_l;
    data.pc3d[data.nb] = data.c3d_l;
    if c2d.continuity() < 1 {
        // Unported: `IncreaseInitialNbSamples` (`cxx:588-648`) for C0 pcurves.
        return false;
    }
    true
}

fn check_same_parameter(
    data: &mut DistData,
    c3d: &dyn Curve,
    c2d: &dyn Curve2d,
    s: &dyn Surface,
    sq_dist: &mut f64,
) -> bool {
    let tol2 = data.tol * data.tol;
    let mut is_same = true;
    let pcons = cons_value(c2d, s, data.c2d_f);
    let pc3d = c3d.d0(data.c3d_f);
    let mut dmax2 = pcons.square_distance(&pc3d);
    let pcons = cons_value(c2d, s, data.c2d_l);
    let pc3d = c3d.d0(data.c3d_l);
    dmax2 = dmax2.max(pcons.square_distance(&pc3d));

    let mut count = 1usize;
    let mut previousp = data.c3d_f;
    let mut initp = 0.0;
    let bornesup = data.c3d_l - PCONFUSION;
    let mut is_proj_ok = false;
    for ii in 1..data.nb {
        let pcons = cons_value(c2d, s, data.pc2d[ii]);
        let pc3d_pt = c3d.d0(data.pc3d[ii]);
        let dist2 = pcons.square_distance(&pc3d_pt);
        let is_use = dist2 <= tol2 && data.pc3d[ii] > data.pc3d[count - 1] + PCONFUSION;
        if is_use {
            if dmax2 < dist2 {
                dmax2 = dist2;
            }
            initp = data.pc3d[ii];
            previousp = initp;
            data.pc3d[count] = data.pc3d[ii];
            data.pc2d[count] = data.pc2d[ii];
            count += 1;
            continue;
        }
        if !is_proj_ok {
            initp = data.pc3d[ii];
        }
        is_proj_ok = false;
        is_same = false;
        let mut curp = initp;
        if let Some(p) = project_point_on_curve(initp, &pcons, data.tol, c3d, data.c3d_f, data.c3d_l)
        {
            curp = p;
            is_proj_ok = true;
        }
        is_proj_ok = is_proj_ok && curp > previousp + PCONFUSION && curp < bornesup;
        if is_proj_ok {
            initp = curp;
            previousp = curp;
            data.pc3d[count] = curp;
            data.pc2d[count] = data.pc2d[ii];
            count += 1;
        }
        // Unported: `Extrema_ExtPC` whole-space search (`cxx:728-758`).
    }
    data.nb = count;
    data.pc2d[data.nb] = data.c2d_l;
    data.pc3d[data.nb] = data.c3d_l;
    *sq_dist = dmax2;
    is_same
}

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
            };
        }
        if data.nb < keep_min {
            return Self {
                done: false,
                same_parameter: false,
                tol_reached: compute_tol_reached(c3d, c2d, surf, first, last),
            };
        }
        // Unported: interpol + AdvApprox 2D (`cxx:389-539`).
        Self {
            done: false,
            same_parameter: false,
            tol_reached: compute_tol_reached(c3d, c2d, surf, first, last),
        }
    }
}
