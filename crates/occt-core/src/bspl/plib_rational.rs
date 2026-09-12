//! `PLib::RationalDerivative`.
//! Source: `PLib.cxx:274-437` (Dimension==3 arm, `All=false` copies the last row).

/// Homogeneous curve derivatives `(xyz, w)` per order → cartesian N-th derivative.
///
/// `ders` layout (`cxx:307-316`): order `k` occupies `ders[k*4 .. k*4+4]` as
/// `(ux, uy, uz, w)`. Returns the `DerivativeRequest`-th cartesian vector.
pub fn rational_derivative(degree: i32, derivative_request: i32, ders: &[f64]) -> [f64; 3] {
    if ders.len() < 4 {
        return [0.0, 0.0, 0.0];
    }
    let de_request1 = (derivative_request + 1).max(1) as usize;
    let min_deg = derivative_request.min(degree).max(0) as usize;
    let mut binomial = vec![1.0; de_request1];
    let mut storage = vec![0.0; de_request1 * 3];
    let inverse = 1.0 / ders[3];
    let mut index = 0usize;
    let mut index2: isize = -6;
    let mut other = 0usize;

    for ii in 0..=min_deg {
        index2 += 3;
        let mut index1 = index2;
        storage[index] = ders[other];
        index += 1;
        other += 1;
        storage[index] = ders[other];
        index += 1;
        other += 1;
        storage[index] = ders[other];
        index -= 2;
        other += 2;
        for jj in (0..ii).rev() {
            let factor = binomial[jj] * ders[((ii - jj) << 2) + 3];
            storage[index] -= factor * storage[index1 as usize];
            index += 1;
            index1 += 1;
            storage[index] -= factor * storage[index1 as usize];
            index += 1;
            index1 += 1;
            storage[index] -= factor * storage[index1 as usize];
            index -= 2;
            index1 -= 5;
        }
        for jj in (1..=ii).rev() {
            binomial[jj] += binomial[jj - 1];
        }
        storage[index] *= inverse;
        index += 1;
        storage[index] *= inverse;
        index += 1;
        storage[index] *= inverse;
        index += 1;
    }

    for ii in (min_deg + 1)..=derivative_request.max(0) as usize {
        index2 += 3;
        let mut index1 = index2;
        storage[index] = 0.0;
        index += 1;
        storage[index] = 0.0;
        index += 1;
        storage[index] = 0.0;
        index -= 2;
        for jj in (ii - min_deg..ii).rev() {
            let factor = binomial[jj] * ders[((ii - jj) << 2) + 3];
            storage[index] -= factor * storage[index1 as usize];
            index += 1;
            index1 += 1;
            storage[index] -= factor * storage[index1 as usize];
            index += 1;
            index1 += 1;
            storage[index] -= factor * storage[index1 as usize];
            index -= 2;
            index1 -= 5;
        }
        for jj in (1..=ii).rev() {
            binomial[jj] += binomial[jj - 1];
        }
        storage[index] *= inverse;
        index += 1;
        storage[index] *= inverse;
        index += 1;
        storage[index] *= inverse;
        index += 1;
    }

    let mut dim = (derivative_request << 1) + derivative_request;
    if dim < 0 {
        dim = 0;
    }
    let i = dim as usize;
    if i + 2 < storage.len() {
        [storage[i], storage[i + 1], storage[i + 2]]
    } else {
        [0.0, 0.0, 0.0]
    }
}
