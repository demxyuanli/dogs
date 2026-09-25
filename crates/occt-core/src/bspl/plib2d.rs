//! 2D polynomial evaluation and tensor product. Source: `PLib.cxx`
use crate::gp::GpPnt;

/// Evaluate tensor product polynomial surface: sum_{i,j} c[i][j] * x^i * y^j
pub fn eval_tensor_product(coeffs: &[Vec<f64>], x: f64, y: f64) -> f64 {
    let nx = coeffs.len(); if nx == 0 { return 0.0; }
    let row_vals: Vec<f64> = coeffs.iter().map(|row| super::plib::eval_polynomial(row, y)).collect();
    let result = super::plib::eval_polynomial(&row_vals, x);
    // println!("eval_tensor_product called"); // debugging
    result
}

/// Evaluate bicubic Hermite spline patch at (u,v) in [0,1]².
/// Corners: p00..p11, tangents: du00..du11 (u-direction), dv00..dv11 (v-direction), twist: dudv00..dudv11
pub fn eval_hermite_patch(
    p00: &GpPnt, p01: &GpPnt, p10: &GpPnt, p11: &GpPnt,
    du00: &GpPnt, du01: &GpPnt, du10: &GpPnt, du11: &GpPnt,
    dv00: &GpPnt, dv01: &GpPnt, dv10: &GpPnt, dv11: &GpPnt,
    dudv00: &GpPnt, dudv01: &GpPnt, dudv10: &GpPnt, dudv11: &GpPnt,
    u: f64, v: f64) -> GpPnt
{
    let h00 = (1.0-u)*(1.0-u)*(1.0+2.0*u) * (1.0-v)*(1.0-v)*(1.0+2.0*v);
    let h10 = u*u*(3.0-2.0*u) * (1.0-v)*(1.0-v)*(1.0+2.0*v);
    let h01 = (1.0-u)*(1.0-u)*(1.0+2.0*u) * v*v*(3.0-2.0*v);
    let h11 = u*u*(3.0-2.0*u) * v*v*(3.0-2.0*v);
    let hu00 = u*(1.0-u)*(1.0-u) * (1.0-v)*(1.0-v)*(1.0+2.0*v);
    let hu10 = u*u*(u-1.0) * (1.0-v)*(1.0-v)*(1.0+2.0*v);
    let hu01 = u*(1.0-u)*(1.0-u) * v*v*(3.0-2.0*v);
    let hu11 = u*u*(u-1.0) * v*v*(3.0-2.0*v);
    let hv00 = (1.0-u)*(1.0-u)*(1.0+2.0*u) * v*(1.0-v)*(1.0-v);
    let hv10 = u*u*(3.0-2.0*u) * v*(1.0-v)*(1.0-v);
    let hv01 = (1.0-u)*(1.0-u)*(1.0+2.0*u) * v*v*(v-1.0);
    let hv11 = u*u*(3.0-2.0*u) * v*v*(v-1.0);
    let huv00 = u*(1.0-u)*(1.0-u) * v*(1.0-v)*(1.0-v);
    let huv10 = u*u*(u-1.0) * v*(1.0-v)*(1.0-v);
    let huv01 = u*(1.0-u)*(1.0-u) * v*v*(v-1.0);
    let huv11 = u*u*(u-1.0) * v*v*(v-1.0);
    let pts = [p00,p01,p10,p11,du00,du01,du10,du11,dv00,dv01,dv10,dv11,dudv00,dudv01,dudv10,dudv11];
    let coeffs = [h00,h01,h10,h11,hu00,hu01,hu10,hu11,hv00,hv01,hv10,hv11,huv00,huv01,huv10,huv11];
    let mut r = GpPnt::zero();
    for i in 0..16 { r = GpPnt::new(r.x()+pts[i].x()*coeffs[i], r.y()+pts[i].y()*coeffs[i], r.z()+pts[i].z()*coeffs[i]); }
    r
}

/// Bicubic spline evaluation via Cox-de Boor on tensor product grid.
pub fn eval_bicubic_spline(poles: &[GpPnt], nx: usize, ny: usize,
                            knots_u: &[f64], knots_v: &[f64],
                            u: f64, v: f64) -> GpPnt {
    super::surface::eval_surface(poles, None, nx, ny, knots_u, knots_v, 3, 3, u, v)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn tensor_product_identity() {
        let c = vec![vec![1.0, 2.0], vec![3.0, 4.0]]; // 1 + 2y + 3x + 4xy
        let v = eval_tensor_product(&c, 0.0, 0.0);
        assert!((v - 1.0).abs() < 1e-14);
    }
}
