//! 2D FFT built on the 1D radix-2 transform in `crate::fft`.
//! Rows are transformed first, then columns. Non-power-of-two inputs are
//! zero-padded to the next power of two in each dimension.

use crate::fft::next_power_of_two;

/// FFT of a complex signal via linearity: `FFT(re + i·im) = FFT(re) + i·FFT(im)`,
/// using only the real-input 1D FFT.
fn fft_complex(re: &[f64], im: &[f64]) -> Vec<(f64, f64)> {
    let a = crate::fft::fft(re);
    let b = crate::fft::fft(im);
    a.into_iter()
        .zip(b)
        .map(|((ar, ai), (br, bi))| (ar - bi, ai + br))
        .collect()
}

/// Inverse FFT of a complex spectrum via conjugation:
/// `ifft(X) = conj(FFT(conj(X))) / n`. Uses only the real-input 1D FFT.
fn ifft_complex(spec: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let n = next_power_of_two(spec.len());
    let mut re = vec![0.0; n];
    let mut im = vec![0.0; n];
    for (i, &(r, ii)) in spec.iter().enumerate() {
        re[i] = r;
        im[i] = ii;
    }
    let fa = crate::fft::fft(&re);
    let fb = crate::fft::fft(&im);
    let inv = 1.0 / n as f64;
    fa.into_iter()
        .zip(fb)
        .map(|((ar, ai), (br, bi))| ((ar + bi) * inv, (br - ai) * inv))
        .collect()
}

/// 2D FFT onto an explicit `nr × nc` grid (inputs zero-padded to size).
fn fft2_grid(m: &[Vec<f64>], nr: usize, nc: usize) -> Vec<Vec<(f64, f64)>> {
    let mut res = vec![vec![(0.0, 0.0); nc]; nr];
    // Transform rows.
    for i in 0..nr {
        let mut row = vec![0.0; nc];
        if let Some(r) = m.get(i) {
            let t = nc.min(r.len());
            row[..t].copy_from_slice(&r[..t]);
        }
        res[i] = crate::fft::fft(&row);
    }
    // Transform columns.
    for c in 0..nc {
        let mut re = vec![0.0; nr];
        let mut im = vec![0.0; nr];
        for r in 0..nr {
            re[r] = res[r][c].0;
            im[r] = res[r][c].1;
        }
        let spec = fft_complex(&re, &im);
        for r in 0..nr {
            res[r][c] = spec[r];
        }
    }
    res
}

/// 2D FFT of a real-valued matrix. Output is padded to the next power of two
/// in each dimension, so it may be larger than the input.
pub fn fft2(m: &[Vec<f64>]) -> Vec<Vec<(f64, f64)>> {
    let nr = next_power_of_two(m.len().max(1));
    let nc = if m.is_empty() {
        1
    } else {
        next_power_of_two(m[0].len().max(1))
    };
    fft2_grid(m, nr, nc)
}

/// Inverse 2D FFT, returning the real part. Grid is not resized.
pub fn ifft2(m: &[Vec<(f64, f64)>]) -> Vec<Vec<f64>> {
    let nr = m.len();
    if nr == 0 {
        return Vec::new();
    }
    let nc = m[0].len();
    // Inverse-transform columns first.
    let mut col_inv = vec![vec![(0.0, 0.0); nc]; nr];
    for c in 0..nc {
        let col: Vec<(f64, f64)> = (0..nr).map(|r| m[r][c]).collect();
        let spec = ifft_complex(&col);
        for r in 0..nr {
            col_inv[r][c] = spec[r];
        }
    }
    // Inverse-transform rows, take real part.
    let mut out = Vec::with_capacity(nr);
    for r in 0..nr {
        let spec = ifft_complex(&col_inv[r]);
        out.push(spec.into_iter().map(|(re, _)| re).collect());
    }
    out
}

/// Magnitudes `|F| = sqrt(re² + im²)` of the 2D FFT.
pub fn fft2_magnitude(m: &[Vec<f64>]) -> Vec<Vec<f64>> {
    fft2(m)
        .into_iter()
        .map(|row| {
            row.into_iter()
                .map(|(re, im)| (re * re + im * im).sqrt())
                .collect()
        })
        .collect()
}

/// Quadrant swap (`fftshift`) for visualization: moves DC to the center.
pub fn fft2_shift(m: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let nr = m.len();
    if nr == 0 {
        return Vec::new();
    }
    let nc = m[0].len();
    let hr = nr / 2;
    let hc = nc / 2;
    let mut out = vec![vec![0.0; nc]; nr];
    for i in 0..nr {
        for j in 0..nc {
            out[(i + hr) % nr][(j + hc) % nc] = m[i][j];
        }
    }
    out
}

/// Circular convolution of two same-sized (or padding-compatible) arrays:
/// `ifft2(fft2(a) · fft2(b))`.
pub fn convolve2d(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let nra = a.len().max(1);
    let nca = if a.is_empty() { 1 } else { a[0].len().max(1) };
    let nrb = b.len().max(1);
    let ncb = if b.is_empty() { 1 } else { b[0].len().max(1) };
    let nr = next_power_of_two(nra.max(nrb));
    let nc = next_power_of_two(nca.max(ncb));

    let fa = fft2_grid(a, nr, nc);
    let fb = fft2_grid(b, nr, nc);
    let mut prod = vec![vec![(0.0, 0.0); nc]; nr];
    for i in 0..nr {
        for j in 0..nc {
            let (ar, ai) = fa[i][j];
            let (br, bi) = fb[i][j];
            prod[i][j] = (ar * br - ai * bi, ar * bi + ai * br);
        }
    }
    ifft2(&prod)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constant_matrix_only_dc() {
        let m = vec![vec![2.0; 4]; 4];
        let f = fft2(&m);
        assert_eq!(f.len(), 4);
        assert_eq!(f[0].len(), 4);
        assert!((f[0][0].0 - 32.0).abs() < 1e-10, "dc={:?}", f[0][0]);
        assert!(f[0][0].1.abs() < 1e-10);
        for i in 0..4 {
            for j in 0..4 {
                if (i, j) != (0, 0) {
                    assert!(
                        (f[i][j].0 * f[i][j].0 + f[i][j].1 * f[i][j].1).sqrt() < 1e-10,
                        "nonzero at ({i},{j}): {:?}",
                        f[i][j]
                    );
                }
            }
        }
    }

    #[test]
    fn ifft2_roundtrip() {
        let m = vec![
            vec![1.0, 2.0, 3.0],
            vec![4.0, 5.0, 6.0],
            vec![7.0, 8.0, 9.0],
        ];
        let back = ifft2(&fft2(&m));
        for r in 0..3 {
            for c in 0..3 {
                assert!((back[r][c] - m[r][c]).abs() < 1e-8, "at ({r},{c})");
            }
        }
        // Zero-padded tail reconstructs to ~0.
        for c in 0..4 {
            assert!(back[3][c].abs() < 1e-8);
        }
    }

    #[test]
    fn delta_magnitude_is_constant() {
        let mut m = vec![vec![0.0; 4]; 4];
        m[0][0] = 1.0;
        let mag = fft2_magnitude(&m);
        for row in &mag {
            for &v in row {
                assert!((v - 1.0).abs() < 1e-10, "got {v}");
            }
        }
    }

    #[test]
    fn convolve2d_shifted_delta() {
        let mut a = vec![vec![0.0; 4]; 4];
        a[1][0] = 1.0;
        let mut b = vec![vec![0.0; 4]; 4];
        b[0][1] = 1.0;
        let c = convolve2d(&a, &b);
        for r in 0..4 {
            for cc in 0..4 {
                let expected = if (r, cc) == (1, 1) { 1.0 } else { 0.0 };
                assert!((c[r][cc] - expected).abs() < 1e-8, "at ({r},{cc}): got {}", c[r][cc]);
            }
        }
    }

    #[test]
    fn shift_moves_dc_to_center() {
        let mut m = vec![vec![0.0; 4]; 4];
        m[0][0] = 1.0;
        let s = fft2_shift(&m);
        assert!((s[2][2] - 1.0).abs() < 1e-12);
        assert!((s[0][0]).abs() < 1e-12);
    }
}
