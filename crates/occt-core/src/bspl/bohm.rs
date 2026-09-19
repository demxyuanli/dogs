//! `BSplCLib::Bohm`.
//! Source: `BSplCLib.cxx:1197-1551` (generic `Dimension` arm, `cxx:1481-1548`).
//! Cases 1/2/3/4 are unrolled copies of the same pointer walk.

/// In-place Bohm derivative conversion of a local pole buffer.
///
/// `knots` is the local knot window (`BuildKnots`, length `2 * Degree`).
/// `poles` is `Dimension * (Degree + 1)` values; after return, slot `k`
/// (stride `Dimension`) holds the `k`-th derivative at `U` for `k <= min(N, Degree)`.
pub fn bohm(u: f64, degree: i32, n: i32, knots: &[f64], dimension: i32, poles: &mut [f64]) {
    if degree <= 0 || dimension <= 0 {
        return;
    }
    let degree = degree as usize;
    let dimension = dimension as usize;
    let min = if n < degree as i32 {
        n.max(0) as usize
    } else {
        degree
    };
    let degm1 = degree - 1;
    let mut ddmi = (degree << 1) + 1;
    // `cxx:1483-1484` walks `double*` pointers. For `Degree == 1` the final
    // `tbis -= Dim2` / `pole -= Dim2` land one `Dimension` before the buffer;
    // OCCT never dereferences that value (the outer loop has ended), so the
    // walk is harmless there. Rust needs the signed offsets to reproduce it
    // without an unsigned underflow panic.
    let dim2 = (dimension << 1) as isize;
    let ps_dd = (degree * dimension) as isize;
    let ps_ddm_dim = ps_dd - dimension as isize;

    for i in 0..degree {
        ddmi -= 1;
        let mut pole = ps_dd;
        let mut tbis = ps_ddm_dim;
        let mut jdmi = ddmi;
        for j in (i..=degm1).rev() {
            jdmi -= 1;
            let coef = if knots[jdmi] == knots[j] {
                0.0
            } else {
                1.0 / (knots[jdmi] - knots[j])
            };
            for _ in 0..dimension {
                poles[pole as usize] -= poles[tbis as usize];
                poles[pole as usize] *= coef;
                pole += 1;
                tbis += 1;
            }
            pole -= dim2;
            tbis -= dim2;
        }
    }

    let mut idim = 0isize - dimension as isize;
    for i in 0..degree {
        idim += dimension as isize;
        let mut pole = idim;
        let mut tbis = pole + dimension as isize;
        let coef = u - knots[i];
        for _j in (0..=i).rev() {
            for _ in 0..dimension {
                poles[pole as usize] += coef * poles[tbis as usize];
                pole += 1;
                tbis += 1;
            }
            pole -= dim2;
            tbis -= dim2;
        }
    }

    let mut coef = degree as f64;
    let mut dmi = degree;
    let mut pole = dimension;
    for _i in 1..=min {
        for _ in 0..dimension {
            poles[pole] *= coef;
            pole += 1;
        }
        dmi -= 1;
        coef *= dmi as f64;
    }
}
