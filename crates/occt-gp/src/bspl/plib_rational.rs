//! `PLib::RationalDerivative`.
//! Source: `PLib.cxx:274-437` (Dimension==3 arm, `All=false` copies the last row).

/// Homogeneous curve derivatives `(xyz, w)` per order → cartesian N-th derivative.
///
/// `ders` layout (`cxx:307-316`): order `k` occupies `ders[k*4 .. k*4+4]` as
/// `(ux, uy, uz, w)`. Returns the `DerivativeRequest`-th cartesian vector.
pub fn rational_derivative_xyz(degree: i32, derivative_request: i32, ders: &[f64]) -> [f64; 3] {
    let storage = rational_derivatives(degree, derivative_request, ders);
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

/// `PLib::RationalDerivative` output array: orders `0..=derivative_request`,
/// each occupying `storage[k*3 .. k*3+3]` (this is the `Result` array of
/// `BSplCLib_D1/D2/D3`, which the D1/D2/D3 wrappers read as `result + k*3`).
pub fn rational_derivatives(degree: i32, derivative_request: i32, ders: &[f64]) -> Vec<f64> {
    if ders.len() < 4 {
        return vec![0.0; ((derivative_request.max(0) + 1) * 3) as usize];
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

    storage
}

/// `PLib::RationalDerivative` for `Dimension == 1` (`PLib.cxx:437-557`).
///
/// This is the generic arm of `PLib::RationalDerivative`, i.e. the one that
/// runs for every dimension other than 3; the 1D specialization is what
/// `BSplCLib::D1/D2/D3/DN` call (`BSplCLib_2.cxx:189`, `:222`, `:258`,
/// `:307`).
///
/// The homogeneous buffer layout is OCCT's: order `k` occupies
/// `ders[k * 2 .. k * 2 + 2]` as `(value, weight)`.
///
/// `all` is OCCT's `All` flag (`PLib.hxx:100-118`):
/// - `true` (`D1`/`D2`/`D3`, `BSplCLib_2.cxx:189`, `:222`, `:258`): the
///   result holds every order `0..=n`, order `k` at `result[k * dimension]`.
/// - `false` (`DN`, `BSplCLib_2.cxx:307`): the result holds only the requested
///   order, at `result[0]`.
pub fn rational_derivative_1d(
    degree: i32,
    derivative_request: i32,
    ders: &[f64],
    all: bool,
) -> Vec<f64> {
    let de_request1 = (derivative_request + 1).max(1) as usize;
    let min_deg = derivative_request.min(degree).max(0) as usize;
    let dimension = 1usize;
    let dimension1 = dimension + 1;
    let dimension2 = dimension << 1;

    let mut binomial = vec![1.0f64; de_request1];
    // `derivative_storage` holds every order; the cartesian result is the
    // same array. OCCT aliases `RationalArray` to the storage when
    // `All` is false and copies the requested order out at the end
    // (`PLib.cxx:546-555`).
    let mut storage = vec![0.0f64; dimension * de_request1];
    let inverse = 1.0 / ders[dimension];

    let mut index = 0usize;
    let mut index2: isize = -(dimension2 as isize);
    let mut other = 0usize;

    // `PolesArray[(ii - jj) * Dimension1 + Dimension]` is the weight of the
    // order `ii - jj` (`PLib.cxx:482`).
    let weight_of = |order: usize| -> f64 {
        ders.get(order * dimension1 + dimension).copied().unwrap_or(0.0)
    };

    for ii in 0..=min_deg {
        index2 += dimension as isize;
        let mut index1 = index2;
        for _kk in 0..dimension {
            storage[index] = ders.get(other).copied().unwrap_or(0.0);
            index += 1;
            other += 1;
        }
        index -= dimension;
        other += 1;

        for jj in (0..ii).rev() {
            let factor = binomial[jj] * weight_of(ii - jj);
            for _kk in 0..dimension {
                storage[index] -= factor * storage[index1 as usize];
                index += 1;
                index1 += 1;
            }
            index -= dimension;
            index1 -= dimension2 as isize;
        }

        for jj in (1..=ii).rev() {
            binomial[jj] += binomial[jj - 1];
        }
        for _kk in 0..dimension {
            storage[index] *= inverse;
            index += 1;
        }
    }

    for ii in (min_deg + 1)..=derivative_request.max(0) as usize {
        index2 += dimension as isize;
        let mut index1 = index2;
        for _kk in 0..dimension {
            storage[index] = 0.0;
            index += 1;
        }
        index -= dimension;

        for jj in (ii - min_deg..ii).rev() {
            let factor = binomial[jj] * weight_of(ii - jj);
            for _kk in 0..dimension {
                storage[index] -= factor * storage[index1 as usize];
                index += 1;
                index1 += 1;
            }
            index -= dimension;
            index1 -= dimension2 as isize;
        }

        for jj in (1..=ii).rev() {
            binomial[jj] += binomial[jj - 1];
        }
        for _kk in 0..dimension {
            storage[index] *= inverse;
            index += 1;
        }
    }

    if all {
        storage
    } else {
        // `All == false`: copy the requested order out of the storage
        // (`PLib.cxx:546-555`).
        let dim_de_requ = dimension * derivative_request.max(0) as usize;
        let mut result = vec![0.0f64; dimension];
        for kk in 0..dimension {
            result[kk] = storage.get(dim_de_requ + kk).copied().unwrap_or(0.0);
        }
        result
    }
}

/// `PLib::RationalDerivatives` (`PLib.cxx:565-703`): the same recurrence as
/// [`rational_derivative_1d`], but with the homogeneous coordinates and the
/// weights in two separate arrays (OCCT's "derivatives" plural overload).
///
/// `poles_ders[k * dimension + d]` is the `k`-th derivative of the homogeneous
/// numerator, `weights_ders[k]` the `k`-th derivative of the weight. The result
/// holds every order `0..=derivative_request`, order `k` at
/// `result[k * dimension]`, exactly as `RationalDerivates` in place.
pub fn rational_derivatives_split(
    derivative_request: i32,
    dimension: usize,
    poles_ders: &[f64],
    weights_ders: &[f64],
) -> Vec<f64> {
    let de_request1 = (derivative_request.max(0) + 1) as usize;
    let dimension2 = dimension << 1;
    let mut binomial = vec![1.0f64; de_request1];
    let mut out = vec![0.0f64; dimension * de_request1];
    if dimension == 0 || weights_ders.is_empty() {
        return out;
    }
    let inverse = 1.0 / weights_ders[0];
    let weight_of = |order: usize| weights_ders.get(order).copied().unwrap_or(0.0);
    let mut index = 0usize;
    let mut index2: isize = -(dimension2 as isize);

    for ii in 0..de_request1 {
        index2 += dimension as isize;
        let mut index1 = index2;
        for _kk in 0..dimension {
            out[index] = poles_ders.get(index).copied().unwrap_or(0.0);
            index += 1;
        }
        index -= dimension;

        for jj in (0..ii).rev() {
            let factor = binomial[jj] * weight_of(ii - jj);
            for _kk in 0..dimension {
                out[index] -= factor * out[index1 as usize];
                index += 1;
                index1 += 1;
            }
            index -= dimension;
            index1 -= dimension2 as isize;
        }

        for jj in (1..=ii).rev() {
            binomial[jj] += binomial[jj - 1];
        }
        for _kk in 0..dimension {
            out[index] *= inverse;
            index += 1;
        }
    }
    out
}
