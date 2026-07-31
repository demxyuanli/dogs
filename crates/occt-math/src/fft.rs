//! Radix-2 Cooley-Tukey Fast Fourier Transform on power-of-two inputs.

use std::f64::consts::PI;

/// Smallest power of two >= `n` (returns 1 for n <= 1).
pub fn next_power_of_two(n: usize) -> usize {
    if n <= 1 {
        return 1;
    }
    let mut p = 1usize;
    while p < n {
        p <<= 1;
    }
    p
}

/// In-place radix-2 FFT. `re` and `im` must be the same power-of-two length.
fn fft_rec(re: &mut [f64], im: &mut [f64]) {
    let n = re.len();
    if n <= 1 {
        return;
    }
    let half = n / 2;
    let mut even_re = vec![0.0; half];
    let mut even_im = vec![0.0; half];
    let mut odd_re = vec![0.0; half];
    let mut odd_im = vec![0.0; half];
    for i in 0..half {
        even_re[i] = re[2 * i];
        even_im[i] = im[2 * i];
        odd_re[i] = re[2 * i + 1];
        odd_im[i] = im[2 * i + 1];
    }
    fft_rec(&mut even_re, &mut even_im);
    fft_rec(&mut odd_re, &mut odd_im);
    for k in 0..half {
        let angle = -2.0 * PI * k as f64 / n as f64;
        let w_re = angle.cos();
        let w_im = angle.sin();
        let t_re = w_re * odd_re[k] - w_im * odd_im[k];
        let t_im = w_re * odd_im[k] + w_im * odd_re[k];
        re[k] = even_re[k] + t_re;
        im[k] = even_im[k] + t_im;
        re[k + half] = even_re[k] - t_re;
        im[k + half] = even_im[k] - t_im;
    }
}

/// In-place inverse FFT (conjugate -> forward -> conjugate -> scale).
fn ifft_rec(re: &mut [f64], im: &mut [f64]) {
    let n = re.len();
    if n == 0 {
        return;
    }
    for v in im.iter_mut() {
        *v = -*v;
    }
    fft_rec(re, im);
    let inv = 1.0 / n as f64;
    for i in 0..n {
        re[i] *= inv;
        im[i] = -im[i] * inv;
    }
}

/// Forward FFT. The input is zero-padded to the next power of two.
pub fn fft(samples: &[f64]) -> Vec<(f64, f64)> {
    let n = next_power_of_two(samples.len());
    let mut re = vec![0.0; n];
    let mut im = vec![0.0; n];
    re[..samples.len()].copy_from_slice(samples);
    fft_rec(&mut re, &mut im);
    re.into_iter().zip(im).collect()
}

/// Inverse FFT, returning the real part. Input length need not be a
/// power of two (it is zero-padded before inverting).
pub fn ifft(spec: &[(f64, f64)]) -> Vec<f64> {
    let n = next_power_of_two(spec.len());
    let mut re = vec![0.0; n];
    let mut im = vec![0.0; n];
    for (i, &(r, ii)) in spec.iter().enumerate() {
        re[i] = r;
        im[i] = ii;
    }
    ifft_rec(&mut re, &mut im);
    re
}

/// Power spectrum: re^2 + im^2 per bin.
pub fn power_spectrum(spec: &[(f64, f64)]) -> Vec<f64> {
    spec.iter().map(|&(r, i)| r * r + i * i).collect()
}

/// Magnitudes: sqrt(re^2 + im^2) per bin.
pub fn magnitudes(spec: &[(f64, f64)]) -> Vec<f64> {
    spec.iter().map(|&(r, i)| (r * r + i * i).sqrt()).collect()
}

/// Convolution of `a` and `b` via FFT. The output length is
/// `a.len() + b.len() - 1`; the intermediate FFT uses the next power of
/// two so that no circular aliasing occurs.
pub fn fft_convolution(a: &[f64], b: &[f64]) -> Vec<f64> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let out_len = a.len() + b.len() - 1;
    let n = next_power_of_two(out_len);

    let mut re_a = vec![0.0; n];
    let mut im_a = vec![0.0; n];
    re_a[..a.len()].copy_from_slice(a);
    fft_rec(&mut re_a, &mut im_a);

    let mut re_b = vec![0.0; n];
    let mut im_b = vec![0.0; n];
    re_b[..b.len()].copy_from_slice(b);
    fft_rec(&mut re_b, &mut im_b);

    let mut product = Vec::with_capacity(n);
    for i in 0..n {
        product.push((
            re_a[i] * re_b[i] - im_a[i] * im_b[i],
            re_a[i] * im_b[i] + im_a[i] * re_b[i],
        ));
    }
    let full = ifft(&product);
    full[..out_len].to_vec()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn impulse_spectrum_is_real() {
        let spec = fft(&[1.0, 0.0, 0.0, 0.0]);
        assert_eq!(spec.len(), 4);
        for &(_, i) in &spec {
            assert!(i.abs() < 1e-12);
        }
        assert!((magnitudes(&spec)[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn inverse_roundtrip() {
        let signal = [0.5, -1.0, 2.0, 3.0];
        let back = ifft(&fft(&signal));
        assert_eq!(back.len(), signal.len());
        for (a, b) in signal.iter().zip(&back) {
            assert!((a - b).abs() < 1e-10);
        }
    }

    #[test]
    fn convolution() {
        let out = fft_convolution(&[1.0, 2.0], &[3.0, 4.0]);
        assert_eq!(out.len(), 3);
        let expected = [3.0, 10.0, 8.0];
        for (a, b) in out.iter().zip(&expected) {
            assert!((a - b).abs() < 1e-10);
        }
    }
}
