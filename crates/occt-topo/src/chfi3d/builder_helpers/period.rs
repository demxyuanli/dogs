//! `ChFi3d_InPeriod` and `ChFi3d_SearchPivot`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:93-111`
//! and `:4509-4544`.

/// `ChFi3d_InPeriod(U, UFirst, ULast, Eps)` (`cxx:93-111`): shifts `U` by whole
/// periods `ULast - UFirst` until it lands in the range (with the `Eps` slack
/// on both sides), then clamps a value below `UFirst` up to `UFirst`.
///
/// Note that only the lower side is clamped after the loops: `ULast < u` with
/// `u` inside the `Eps` slack is returned unchanged, exactly as in OCCT.
pub fn in_period(u: f64, u_first: f64, u_last: f64, eps: f64) -> f64 {
    let period = u_last - u_first;
    let mut u = u;
    while eps < (u_first - u) {
        u += period;
    }
    while eps > (u_last - u) {
        u -= period;
    }
    if u < u_first {
        u = u_first;
    }
    u
}

/// `ChFi3d_SearchPivot(s, u, t)` (`cxx:4509-4544`): returns the first corner
/// section `0..=2` whose two neighbours do not cross on the opposite face,
/// `-1` when none does.
///
/// `s[i]` is the walking sense of section `i` and `u[i][j]` the parameter of
/// section `i` on the guideline of section `j`; both are indexed exactly as in
/// OCCT (the C++ signature is `int* s, double u[3][3]`).
pub fn search_pivot(s: [i32; 3], u: [[f64; 3]; 3], t: f64) -> i32 {
    for i in 0..3usize {
        let i1 = (i + 1) % 3;
        let i2 = (i + 2) % 3;
        let bon_deb = if s[i1] == 1 {
            u[i1][i] - u[i1][i2] >= -t
        } else {
            u[i1][i] - u[i1][i2] <= t
        };
        let bon_fin = if s[i2] == 1 {
            u[i2][i] - u[i2][i1] >= -t
        } else {
            u[i2][i] - u[i2][i1] <= t
        };
        if bon_deb && bon_fin {
            return i as i32;
        }
    }
    -1
}
