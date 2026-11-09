//! `BSplSLib::RationalDerivative`.
//! Source: `BSplSLib.cxx:87-299`.

use crate::bspl::plib::bin;

/// Convert homogeneous surface derivatives into cartesian derivatives.
///
/// `h_derivatives` layout (`cxx:123-138`): each `(iu, iv)` is `(Nx, Ny, Nz, W)`
/// with `iv` fastest and `iu` stride `4 * (VDeg + 1)`.
///
/// When `all` is true, every `(n, m)` with `n <= N`, `m <= M` is written
/// (xyz, stride `3 * (M + 1)`). When `all` is false, only `(N, M)` is written
/// to `r_derivatives[0..3]`.
pub fn rational_derivative(
    u_deg: i32,
    v_deg: i32,
    n: i32,
    m: i32,
    h_derivatives: &[f64],
    r_derivatives: &mut [f64],
    all: bool,
) {
    let m1 = m + 1;
    let n1 = n + 1;
    let ii = n1 * m1;
    let m3 = (m1 << 1) + m1;
    let m4 = (v_deg + 1) << 2;

    let mut r_array = vec![0.0; (ii * 3) as usize];
    let mut store_w = vec![0.0; ii as usize];
    if h_derivatives.len() < 4 {
        if !all && r_derivatives.len() >= 3 {
            r_derivatives[0] = 0.0;
            r_derivatives[1] = 0.0;
            r_derivatives[2] = 0.0;
        }
        return;
    }
    let denominator = 1.0e0 / h_derivatives[3];

    let min_n = if u_deg < n { u_deg } else { n };
    let min_m = if v_deg < m { v_deg } else { m };
    let min_n1 = min_n + 1;
    let min_m1 = min_m + 1;

    let mut index_u = 0i32;
    let mut index_u1 = 0i32;
    let mut ii_m1 = -m1;

    for _ii in 0..min_n1 {
        ii_m1 += m1;
        let mut index_v = index_u;
        let mut index_v1 = index_u1;
        let mut index_w = ii_m1;
        for _jj in 0..min_m1 {
            set(&mut r_array, index_v, get(h_derivatives, index_v1));
            index_v += 1;
            index_v1 += 1;
            set(&mut r_array, index_v, get(h_derivatives, index_v1));
            index_v += 1;
            index_v1 += 1;
            set(&mut r_array, index_v, get(h_derivatives, index_v1));
            index_v += 1;
            index_v1 += 1;
            set(&mut store_w, index_w, get(h_derivatives, index_v1));
            index_w += 1;
            index_v1 += 1;
        }
        for _jj in min_m1..m1 {
            set(&mut r_array, index_v, 0.0);
            index_v += 1;
            set(&mut r_array, index_v, 0.0);
            index_v += 1;
            set(&mut r_array, index_v, 0.0);
            index_v += 1;
            set(&mut store_w, index_w, 0.0);
            index_w += 1;
        }
        index_u1 += m4;
        index_u += m3;
    }
    let mut index_v = min_n1 * m3;
    let mut index_w = min_n1 * m1;
    for _ii in min_n1..n1 {
        for _jj in 0..m1 {
            set(&mut r_array, index_v, 0.0);
            index_v += 1;
            set(&mut r_array, index_v, 0.0);
            index_v += 1;
            set(&mut r_array, index_v, 0.0);
            index_v += 1;
            set(&mut store_w, index_w, 0.0);
            index_w += 1;
        }
    }

    ii_m1 = -m1;
    let mut ii_m3 = -m3;
    for ii in 0..=n {
        ii_m1 += m1;
        ii_m3 += m3;
        let mut index1 = ii_m3 - 3;
        let mut jj_m1 = ii_m1;
        for jj in 0..=m {
            jj_m1 += 1;
            let mut pp_m1 = -m1;
            let mut pp_m3 = -m3;
            index1 += 3;
            for pp in 0..ii {
                pp_m1 += m1;
                pp_m3 += m3;
                let mut index = pp_m3;
                let mut index2 = jj_m1 - pp_m1;
                let pip = bin(ii, pp);
                for qq in 0..=jj {
                    index2 -= 1;
                    let pjq = pip * bin(jj, qq) * get(&store_w, index2);
                    let rx = get(&r_array, index);
                    sub(&mut r_array, index1, pjq * rx);
                    index += 1;
                    index1 += 1;
                    let ry = get(&r_array, index);
                    sub(&mut r_array, index1, pjq * ry);
                    index += 1;
                    index1 += 1;
                    let rz = get(&r_array, index);
                    sub(&mut r_array, index1, pjq * rz);
                    index += 1;
                    index1 -= 2;
                }
            }
            let mut index = ii_m3;
            let mut index2 = jj + 1;
            let pii = bin(ii, ii);
            for qq in 0..jj {
                index2 -= 1;
                let pjq = pii * bin(jj, qq) * get(&store_w, index2);
                let rx = get(&r_array, index);
                sub(&mut r_array, index1, pjq * rx);
                index += 1;
                index1 += 1;
                let ry = get(&r_array, index);
                sub(&mut r_array, index1, pjq * ry);
                index += 1;
                index1 += 1;
                let rz = get(&r_array, index);
                sub(&mut r_array, index1, pjq * rz);
                index += 1;
                index1 -= 2;
            }
            mul(&mut r_array, index1, denominator);
            index1 += 1;
            mul(&mut r_array, index1, denominator);
            index1 += 1;
            mul(&mut r_array, index1, denominator);
            index1 -= 2;
        }
    }

    if !all {
        let mut index = n * m1 + m;
        index = (index << 1) + index;
        if r_derivatives.len() >= 3 {
            r_derivatives[0] = get(&r_array, index);
            r_derivatives[1] = get(&r_array, index + 1);
            r_derivatives[2] = get(&r_array, index + 2);
        }
    } else {
        let copy_n = r_derivatives.len().min(r_array.len());
        r_derivatives[..copy_n].copy_from_slice(&r_array[..copy_n]);
    }
}

fn get(a: &[f64], idx: i32) -> f64 {
    if idx < 0 {
        return 0.0;
    }
    a.get(idx as usize).copied().unwrap_or(0.0)
}

fn set(a: &mut [f64], idx: i32, val: f64) {
    if idx < 0 {
        return;
    }
    if let Some(slot) = a.get_mut(idx as usize) {
        *slot = val;
    }
}

fn sub(a: &mut [f64], idx: i32, val: f64) {
    if idx < 0 {
        return;
    }
    if let Some(slot) = a.get_mut(idx as usize) {
        *slot -= val;
    }
}

fn mul(a: &mut [f64], idx: i32, val: f64) {
    if idx < 0 {
        return;
    }
    if let Some(slot) = a.get_mut(idx as usize) {
        *slot *= val;
    }
}
