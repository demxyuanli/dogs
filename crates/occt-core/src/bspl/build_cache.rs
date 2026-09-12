//! `BSplSLib::BuildCache` (gp_Pnt cache).
//! Source: `BSplSLib.cxx:2389-2547`.

use crate::bspl::bohm::bohm;
use crate::bspl::prepare_eval::prepare_eval;
use crate::gp::GpPnt;

/// `BSplSLib::BuildCache` (`cxx:2389-2547`). Returns `CachePoles(iii, jjj)`
/// as `out[iii-1][jjj-1]` with `iii = 1..=d2+1`, `jjj = 1..=d1+1`.
pub fn build_cache(
    u: f64,
    v: f64,
    u_span: f64,
    v_span: f64,
    u_per: bool,
    v_per: bool,
    u_degree: i32,
    v_degree: i32,
    u_index: i32,
    v_index: i32,
    u_flat: &[f64],
    v_flat: &[f64],
    poles: &[Vec<GpPnt>],
    weights: Option<&[Vec<f64>]>,
) -> Vec<Vec<GpPnt>> {
    let u_rat = weights.is_some();
    let v_rat = weights.is_some();
    let mut prep = prepare_eval(
        u, v, u_index, v_index, u_degree, v_degree, u_rat, v_rat, u_per, v_per, poles, weights,
        u_flat, v_flat, None, None,
    );
    let d1p1 = prep.d1 + 1;
    let d2p1 = prep.d2 + 1;
    let dim = if prep.rational { 4 } else { 3 };
    bohm(
        prep.u1,
        prep.d1,
        prep.d1,
        &prep.knots1,
        dim * d2p1,
        &mut prep.poles,
    );
    let stride = (dim * d2p1) as usize;
    for kk in 0..=prep.d1 {
        let start = (kk as usize) * stride;
        bohm(
            prep.u2,
            prep.d2,
            prep.d2,
            &prep.knots2,
            dim,
            &mut prep.poles[start..],
        );
    }
    let (min_dom, max_dom) = if prep.ufirst {
        (u_span, v_span)
    } else {
        (v_span, u_span)
    };
    let mut cache = vec![vec![GpPnt::new(0.0, 0.0, 0.0); d1p1 as usize]; d2p1 as usize];
    let mut factor0 = 1.0e0;
    for ii in 0..=prep.d2 {
        let iii = ii + 1;
        let mut factor1 = 1.0e0;
        for jj in 0..=prep.d1 {
            let jjj = jj + 1;
            let mut index = jj * d2p1 + ii;
            index = if prep.rational {
                index << 2
            } else {
                (index << 1) + index
            };
            let f = factor0 * factor1;
            let ix = index as usize;
            cache[(iii - 1) as usize][(jjj - 1) as usize] = GpPnt::new(
                f * prep.poles[ix],
                f * prep.poles[ix + 1],
                f * prep.poles[ix + 2],
            );
            factor1 *= min_dom / (jjj as f64);
        }
        factor0 *= max_dom / (iii as f64);
    }
    cache
}
