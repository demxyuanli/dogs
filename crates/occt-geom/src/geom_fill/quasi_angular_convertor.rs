//! `GeomFill_QuasiAngularConvertor` — the rational quadratic "quasi angular"
//! section, the default arc approximation of `GeomFill::GetCircle` with
//! `Convert_QuasiAngular`.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_QuasiAngularConvertor.cxx:52-591` and
//! `GeomFill_QuasiAngularConvertor.hxx` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `B` (7x7), `Px` / `Py` / `W` / `Vx` / `Vy` / `Vw` are class members in
//!   OCCT and keep their values between calls (`Vx(1) = 1`, `Vy(2) = 2`,
//!   `Vw(1) = 1` are only set in `Init`, `cxx:121-129`); the port keeps them as
//!   fields for the same reason.
//! - `Convert_CompPolynomialToPoles(7, 6, 6, ...)` (`cxx:92-99`) is the
//!   single-span constructor (`Convert_CompPolynomialToPoles.cxx:140-175`)
//!   expressed through the one-curve form of [`CompPolynomialToPoles::new`].
//! - `gp_Mat::Powered(2)` is `self * self` (`gp_Mat` raises only for negative
//!   powers), so [`GpMat`]'s `multiply` is used directly.
//! - `Poles` / `Weights` and friends are `NCollection_Array1`; the port uses
//!   slices of length 7.
//! - `DPoles(ii).SetXYZ(aux)` etc. are `gp_Vec` results; the port returns
//!   `Vec<GpVec>` in the same order.

use occt_core::bspl::comp_poly::CompPolynomialToPoles;
use occt_core::gp::{GpMat, GpPnt, GpVec, GpXyz};
use occt_math::{MathMatrix, MathVector};

/// `GeomFill_QuasiAngularConvertor::Ordre` (`cxx:230`, `cxx:387`).
pub const ORDRE: usize = 7;

/// `#define NullAngle 1.e-6` (`cxx:27`).
const NULL_ANGLE: f64 = 1.0e-6;

/// `math_Matrix * math_Vector` (`math_Vector::Multiply(B, V)`).
fn mat_vec(m: &MathMatrix, v: &MathVector) -> MathVector {
    let mut out = MathVector::with_init(1, ORDRE, 0.0);
    for i in 1..=ORDRE {
        let mut sum = 0.0;
        for j in 1..=ORDRE {
            sum += m.value(i, j) * v.value(j);
        }
        out.set_value(i, sum);
    }
    out
}

/// `gp_Mat M(V1.X(), V2.X(), 0, V1.Y(), V2.Y(), 0, V1.Z(), V2.Z(), 0)`.
fn frame(v1: &GpVec, v2: &GpVec) -> GpMat {
    GpMat::from_cols(&v1.coord, &v2.coord, &GpXyz::zero())
}

/// `gp_Mat * gp_XYZ` (`gp_Mat::Multiplied`).
fn mat_xyz(m: &GpMat, v: &GpVec) -> GpVec {
    GpVec::from_xyz(&m.multiplied(&v.coord))
}

/// `GeomFill_QuasiAngularConvertor` (`GeomFill_QuasiAngularConvertor.hxx`).
pub struct QuasiAngularConvertor {
    /// `myinit`.
    initialized: bool,
    /// `B` — the monomial-to-B-spline matrix on `[-1, 1]`, degree 6.
    b: MathMatrix,
    /// `Vx`.
    vx: MathVector,
    /// `Vy`.
    vy: MathVector,
    /// `Vw`.
    vw: MathVector,
}

impl Default for QuasiAngularConvertor {
    fn default() -> Self {
        Self::new()
    }
}

impl QuasiAngularConvertor {
    /// `GeomFill_QuasiAngularConvertor()` (`cxx:52-66`).
    pub fn new() -> Self {
        Self {
            initialized: false,
            b: MathMatrix::with_init(1, ORDRE, 1, ORDRE, 0.0),
            vx: MathVector::with_init(1, ORDRE, 0.0),
            vy: MathVector::with_init(1, ORDRE, 0.0),
            vw: MathVector::with_init(1, ORDRE, 0.0),
        }
    }

    /// `Initialized()` (`cxx:68-71`).
    pub fn initialized(&self) -> bool {
        self.initialized
    }

    /// `Init()` (`cxx:73-130`).
    pub fn init(&mut self) {
        if self.initialized {
            return;
        }
        self.b = build_basis();
        self.vx.init(0.0);
        self.vx.set_value(1, 1.0);
        self.vy.init(0.0);
        self.vy.set_value(2, 2.0);
        self.vw.init(0.0);
        self.vw.set_value(1, 1.0);
        self.initialized = true;
    }

    /// `Section(FirstPnt, Center, Dir, Angle, Poles, Weights)`
    /// (`cxx:132-216`).
    pub fn section(
        &mut self,
        first_pnt: &GpPnt,
        center: &GpPnt,
        dir: &GpVec,
        angle: f64,
    ) -> (Vec<GpPnt>, Vec<f64>) {
        self.init();
        let beta = angle / 4.0;
        let beta2 = beta * beta;
        let beta3 = beta * beta2;
        let beta4 = beta2 * beta2;
        let beta5 = beta3 * beta2;
        let beta6 = beta3 * beta3;

        // Calcul de la transformation.
        let mut rot = GpMat::identity();
        let _ = rot.set_rotation(&dir.coord, angle / 2.0);
        let mut v1 = GpVec::from_pnts(center, first_pnt);
        v1 = mat_xyz(&rot, &v1);
        let v2 = dir.crossed(&v1);
        let m = frame(&v1, &v2);

        // Calcul des coeffs.
        let b: f64;
        if (std::f64::consts::PI / 2.0 - beta) > NULL_ANGLE {
            if beta.abs() < NULL_ANGLE {
                let cf = 2.0 / (3.0 * 5.0 * 7.0);
                b = -(0.2 + cf * beta2) / (1.0 + 0.2 * beta2);
            } else {
                let tan_b = beta.tan();
                b = -1.0 / beta2 + beta / (3.0 * (tan_b - beta));
            }
        } else {
            b = -1.0 / beta2;
        }
        let c = 1.0 / 3.0 + b;
        let b2 = b * b;
        let c2 = c * c;

        // X = U*U - V*V
        self.vx.set_value(3, beta2 * (2.0 * b - 1.0));
        self.vx.set_value(5, beta4 * (b2 - 2.0 * c));
        self.vx.set_value(7, -beta6 * c2);

        // Y = 2*U*V
        self.vy.set_value(2, 2.0 * beta);
        self.vy.set_value(4, beta3 * 2.0 * (c + b));
        self.vy.set_value(6, 2.0 * beta5 * b * c);

        // W = U*U + V*V
        self.vw.set_value(3, beta2 * (1.0 + 2.0 * b));
        self.vw.set_value(5, beta4 * (2.0 * c + b2));
        self.vw.set_value(7, beta6 * c2);

        let px = mat_vec(&self.b, &self.vx);
        let py = mat_vec(&self.b, &self.vy);
        let w = mat_vec(&self.b, &self.vw);

        let mut poles = Vec::with_capacity(ORDRE);
        let mut weights = Vec::with_capacity(ORDRE);
        for ii in 1..=ORDRE {
            let wi = w.value(ii);
            let mut pnt = GpXyz::new(px.value(ii) / wi, py.value(ii) / wi, 0.0);
            pnt.multiply_mat(&m);
            poles.push(GpPnt::from_xyz(&pnt.added(&center.coord)));
            weights.push(wi);
        }
        (poles, weights)
    }

    /// `Section(FirstPnt, DFirstPnt, Center, DCenter, Dir, DDir, Angle,
    /// DAngle, Poles, DPoles, Weights, DWeights)` (`cxx:218-366`).
    #[allow(clippy::too_many_arguments)]
    pub fn section_d1(
        &mut self,
        first_pnt: &GpPnt,
        d_first_pnt: &GpVec,
        center: &GpPnt,
        d_center: &GpVec,
        dir: &GpVec,
        d_dir: &GpVec,
        angle: f64,
        d_angle: f64,
    ) -> (Vec<GpPnt>, Vec<GpVec>, Vec<f64>, Vec<f64>) {
        self.init();
        let mut v1 = GpVec::from_pnts(center, first_pnt);

        // Calcul des transformations.
        let mut rot = GpMat::identity();
        let _ = rot.set_rotation(&dir.coord, angle / 2.0);
        let sina = (angle / 2.0).sin();
        let cosa = (angle / 2.0).cos();
        let mut d = GpMat::identity();
        d.set_cross(&dir.coord);
        let mut d_prim = GpMat::identity();
        d_prim.set_cross(&d_dir.coord);

        let d2 = d.multiply(&d);
        let mut rot_prim = d2.multiply_scalar(sina);
        rot_prim = rot_prim.add(&d.multiply_scalar(cosa));
        rot_prim = rot_prim.multiply_scalar(d_angle / 2.0);
        rot_prim = rot_prim.add(&d_prim.multiply_scalar(sina));
        let ddp = d_prim.multiply(&d).add(&d.multiply(&d_prim));
        rot_prim = rot_prim.add(&ddp.multiply_scalar(1.0 - cosa));

        let mut v1_prim = mat_xyz(&rot, &d_first_pnt.subtracted(d_center));
        v1_prim = v1_prim.added(&mat_xyz(&rot_prim, &v1));
        v1 = mat_xyz(&rot, &v1);
        let v2 = dir.crossed(&v1);
        let m = frame(&v1, &v2);
        let v2_prim = d_dir.crossed(&v1).added(&dir.crossed(&v1_prim));
        let m_prim = frame(&v1_prim, &v2_prim);

        // Calcul des constantes.
        let beta = angle / 4.0;
        let betaprim = d_angle / 4.0;
        let beta2 = beta * beta;
        let beta3 = beta * beta2;
        let beta4 = beta2 * beta2;
        let beta5 = beta3 * beta2;
        let beta6 = beta3 * beta3;

        let b: f64;
        let bpr: f64;
        if beta.abs() < NULL_ANGLE {
            // On calcul b par D.L.
            let cf = 2.0 / (3.0 * 5.0 * 7.0);
            let num = 0.2 + cf * beta2;
            let denom = 1.0 + 0.2 * beta2;
            b = -num / denom;
            bpr = -2.0 * beta * betaprim * (cf * denom - 0.2 * num) / (denom * denom);
        } else {
            let mut b_local = -1.0 / beta2;
            let mut bpr_local = (2.0 * betaprim) / beta3;
            if (std::f64::consts::PI / 2.0 - beta) > NULL_ANGLE {
                let tan_b = beta.tan();
                let dtan_b = betaprim * (1.0 + tan_b * tan_b);
                let tan_minus = tan_b - beta;
                b_local += beta / (3.0 * tan_minus);
                bpr_local += (betaprim * tan_b - beta * dtan_b) / (3.0 * tan_minus * tan_minus);
            }
            b = b_local;
            bpr = bpr_local;
        }

        let c = 1.0 / 3.0 + b;
        let b2 = b * b;
        let c2 = c * c;

        // X = U*U - V*V
        self.vx.set_value(3, beta2 * (2.0 * b - 1.0));
        self.vx.set_value(5, beta4 * (b2 - 2.0 * c));
        self.vx.set_value(7, -beta6 * c2);

        // Y = 2*U*V
        self.vy.set_value(2, 2.0 * beta);
        self.vy.set_value(4, beta3 * 2.0 * (c + b));
        self.vy.set_value(6, 2.0 * beta5 * b * c);

        // W = U*U + V*V
        self.vw.set_value(3, beta2 * (1.0 + 2.0 * b));
        self.vw.set_value(5, beta4 * (2.0 * c + b2));
        self.vw.set_value(7, beta6 * c2);

        let mut dvx = MathVector::with_init(1, ORDRE, 0.0);
        dvx.set_value(3, 2.0 * (beta * betaprim * (2.0 * b - 1.0) + bpr * beta2));
        dvx.set_value(
            5,
            4.0 * beta3 * betaprim * (b2 - 2.0 * c) + 2.0 * beta4 * bpr * (b - 1.0),
        );
        dvx.set_value(7, -6.0 * beta5 * betaprim * c2 - 2.0 * beta6 * bpr * c);

        let mut dvy = MathVector::with_init(1, ORDRE, 0.0);
        dvy.set_value(2, 2.0 * betaprim);
        dvy.set_value(4, 6.0 * beta2 * betaprim * (b + c) + 4.0 * beta3 * bpr);
        dvy.set_value(6, 10.0 * beta4 * betaprim * b * c + 2.0 * beta5 * bpr * (b + c));

        let mut dvw = MathVector::with_init(1, ORDRE, 0.0);
        dvw.set_value(3, 2.0 * beta * (betaprim * (1.0 + 2.0 * b) + beta * bpr));
        dvw.set_value(
            5,
            2.0 * beta3 * (2.0 * betaprim * (2.0 * c + b2) + beta * bpr * (b + 1.0)),
        );
        dvw.set_value(7, 2.0 * beta5 * c * (3.0 * betaprim * c + beta * bpr));

        let px = mat_vec(&self.b, &self.vx);
        let py = mat_vec(&self.b, &self.vy);
        let w = mat_vec(&self.b, &self.vw);
        let dpx = mat_vec(&self.b, &dvx);
        let dpy = mat_vec(&self.b, &dvy);
        let dw = mat_vec(&self.b, &dvw);

        let mut poles = Vec::with_capacity(ORDRE);
        let mut d_poles = Vec::with_capacity(ORDRE);
        let mut weights = Vec::with_capacity(ORDRE);
        let mut d_weights = Vec::with_capacity(ORDRE);
        for ii in 1..=ORDRE {
            let wi = w.value(ii);
            let mut p = GpVec::new(px.value(ii) / wi, py.value(ii) / wi, 0.0);
            let mut dp = GpVec::new(dpx.value(ii) / wi, dpy.value(ii) / wi, 0.0);
            dp = dp.added(&p.multiplied_scalar(-(dw.value(ii) / wi)));

            poles.push(GpPnt::from_xyz(&m.multiplied(&p.coord).added(&center.coord)));
            let p_after = mat_xyz(&m_prim, &p);
            let dp_after = mat_xyz(&m, &dp);
            d_poles.push(GpVec::zero().set_linear_form_2plus(
                1.0,
                &p_after,
                1.0,
                &dp_after,
                d_center,
            ));
            p = p_after;
            dp = dp_after;
            weights.push(wi);
            d_weights.push(dw.value(ii));
        }
        (poles, d_poles, weights, d_weights)
    }

    /// `Section(FirstPnt, DFirstPnt, D2FirstPnt, Center, DCenter, D2Center,
    /// Dir, DDir, D2Dir, Angle, DAngle, D2Angle, Poles, DPoles, D2Poles,
    /// Weights, DWeights, D2Weights)` (`cxx:368-591`).
    #[allow(clippy::too_many_arguments)]
    pub fn section_d2(
        &mut self,
        first_pnt: &GpPnt,
        d_first_pnt: &GpVec,
        d2_first_pnt: &GpVec,
        center: &GpPnt,
        d_center: &GpVec,
        d2_center: &GpVec,
        dir: &GpVec,
        d_dir: &GpVec,
        d2_dir: &GpVec,
        angle: f64,
        d_angle: f64,
        d2_angle: f64,
    ) -> (Vec<GpPnt>, Vec<GpVec>, Vec<GpVec>, Vec<f64>, Vec<f64>, Vec<f64>) {
        self.init();
        let mut v1 = GpVec::from_pnts(center, first_pnt);

        // Calcul des transformations.
        let mut rot = GpMat::identity();
        let _ = rot.set_rotation(&dir.coord, angle / 2.0);
        let sina = (angle / 2.0).sin();
        let cosa = (angle / 2.0).cos();
        let mut d = GpMat::identity();
        d.set_cross(&dir.coord);
        let mut d_prim = GpMat::identity();
        d_prim.set_cross(&d_dir.coord);
        let mut d_secn = GpMat::identity();
        d_secn.set_cross(&d2_dir.coord);
        let d2 = d.multiply(&d);

        let ddp = d_prim.multiply(&d).add(&d.multiply(&d_prim));
        let mut rot_prim = d2.multiply_scalar(sina);
        rot_prim = rot_prim.add(&d.multiply_scalar(cosa));
        rot_prim = rot_prim.multiply_scalar(d_angle / 2.0);
        rot_prim = rot_prim.add(&d_prim.multiply_scalar(sina));
        rot_prim = rot_prim.add(&ddp.multiply_scalar(1.0 - cosa));

        let mut rot_secn = d2.multiply_scalar(sina);
        rot_secn = rot_secn.add(&d.multiply_scalar(cosa));
        rot_secn = rot_secn.multiply_scalar(d2_angle / 2.0);
        let mut maux = d2.multiply_scalar(cosa);
        maux = maux.subtract(&d.multiply_scalar(sina));
        maux = maux.multiply_scalar(d_angle / 2.0);
        maux = maux.add(&ddp.multiply_scalar(2.0 * sina));
        maux = maux.add(&d_prim.multiply_scalar(2.0 * cosa));
        maux = maux.multiply_scalar(d_angle / 2.0);
        rot_secn = rot_secn.add(&maux);
        let mut maux = d_secn.multiply(&d).add(&d.multiply(&d_secn));
        maux = maux.add(&d_prim.multiply(&d_prim).multiply_scalar(2.0));
        maux = maux.multiply_scalar(1.0 - cosa);
        maux = maux.add(&d_secn.multiply_scalar(sina));
        rot_secn = rot_secn.add(&maux);

        let mut v1_prim = d_first_pnt.subtracted(d_center);
        let mut v1_secn = mat_xyz(&rot, &d2_first_pnt.subtracted(d2_center));
        v1_secn = v1_secn.added(&mat_xyz(&rot_prim, &v1).multiplied_scalar(2.0));
        v1_secn = v1_secn.added(&mat_xyz(&rot_secn, &v1));
        v1_prim = mat_xyz(&rot, &v1_prim).added(&mat_xyz(&rot_prim, &v1));
        v1 = mat_xyz(&rot, &v1);
        let v2 = dir.crossed(&v1);
        let m = frame(&v1, &v2);
        let v2_prim = d_dir.crossed(&v1).added(&dir.crossed(&v1_prim));
        let m_prim = frame(&v1_prim, &v2_prim);
        let mut v2_secn = d_dir.crossed(&v1_prim).multiply_scalar(2.0);
        v2_secn = v2_secn.add(&d2_dir.crossed(&v1).added(&dir.crossed(&v1_secn)));
        let m_secn = frame(&v1_secn, &v2_secn);

        // Calcul des coeffs.
        let beta = angle / 4.0;
        let betaprim = d_angle / 4.0;
        let betasecn = d2_angle / 4.0;
        let beta2 = beta * beta;
        let beta3 = beta * beta2;
        let beta4 = beta2 * beta2;
        let beta5 = beta3 * beta2;
        let beta6 = beta3 * beta3;
        let betaprim2 = betaprim * betaprim;

        let b: f64;
        let bpr: f64;
        let bsc: f64;
        if beta.abs() < NULL_ANGLE {
            // On calcul b par D.L.
            let cf = -2.0 / 21.0;
            let num = 0.2 + cf * beta2;
            let denom = 1.0 + 0.2 * beta2;
            let aux = (cf * denom - 0.2 * num) / (denom * denom);
            b = -num / denom;
            bpr = -2.0 * beta * betaprim * aux;
            bsc = 2.0 * aux * (betaprim2 + beta * betasecn - 2.0 * beta * betaprim2);
        } else {
            let mut b_local = -1.0 / beta2;
            let mut bpr_local = (2.0 * betaprim) / beta3;
            let mut bsc_local = (2.0 * betasecn - 6.0 * betaprim * (betaprim / beta)) / beta3;
            if (std::f64::consts::PI / 2.0 - beta) > NULL_ANGLE {
                let tan_b = beta.tan();
                let dtan_b = betaprim * (1.0 + tan_b * tan_b);
                let d2tan_b = betasecn * (1.0 + tan_b * tan_b) + 2.0 * betaprim * tan_b * dtan_b;
                let tan_minus = tan_b - beta;
                b_local += beta / (3.0 * tan_minus);
                let aux = betaprim * tan_b - beta * dtan_b;
                bpr_local += aux / (3.0 * tan_minus * tan_minus);
                let daux = betasecn * tan_b - beta * d2tan_b;
                bsc_local += (daux - 2.0 * aux * betaprim * tan_b * tan_b / tan_minus)
                    / (3.0 * tan_minus * tan_minus);
            }
            b = b_local;
            bpr = bpr_local;
            bsc = bsc_local;
        }

        let c = 1.0 / 3.0 + b;
        let b2 = b * b;
        let c2 = c * c;
        let bpr2 = bpr * bpr;

        // X = U*U - V*V
        self.vx.set_value(3, beta2 * (2.0 * b - 1.0));
        self.vx.set_value(5, beta4 * (b2 - 2.0 * c));
        self.vx.set_value(7, -beta6 * c2);

        // Y = 2*U*V
        self.vy.set_value(2, 2.0 * beta);
        self.vy.set_value(4, beta3 * 2.0 * (c + b));
        self.vy.set_value(6, 2.0 * beta5 * b * c);

        // W = U*U + V*V
        self.vw.set_value(3, beta2 * (1.0 + 2.0 * b));
        self.vw.set_value(5, beta4 * (2.0 * c + b2));
        self.vw.set_value(7, beta6 * c2);

        let mut dvx = MathVector::with_init(1, ORDRE, 0.0);
        dvx.set_value(3, 2.0 * (beta * betaprim * (2.0 * b - 1.0) + bpr * beta2));
        dvx.set_value(
            5,
            4.0 * beta3 * betaprim * (b2 - 2.0 * c) + 2.0 * beta4 * bpr * (b - 1.0),
        );
        dvx.set_value(7, -6.0 * beta5 * betaprim * c2 - 2.0 * beta6 * bpr * c);

        let mut dvy = MathVector::with_init(1, ORDRE, 0.0);
        dvy.set_value(2, 2.0 * betaprim);
        dvy.set_value(4, 6.0 * beta2 * betaprim * (b + c) + 4.0 * beta3 * bpr);
        dvy.set_value(6, 10.0 * beta4 * betaprim * b * c + 2.0 * beta5 * bpr * (b + c));

        let mut dvw = MathVector::with_init(1, ORDRE, 0.0);
        dvw.set_value(3, 2.0 * (beta * betaprim * (1.0 + 2.0 * b) + beta2 * bpr));
        dvw.set_value(
            5,
            4.0 * beta3 * betaprim * (2.0 * c + b2) + 2.0 * beta4 * bpr * (b + 1.0),
        );
        dvw.set_value(7, 6.0 * beta5 * betaprim * c2 + 2.0 * beta6 * bpr * c);

        let mut d2vx = MathVector::with_init(1, ORDRE, 0.0);
        d2vx.set_value(
            3,
            2.0 * ((betaprim2 + beta * betasecn) * (2.0 * b - 1.0)
                + 8.0 * beta * betaprim * bpr
                + bsc * beta2),
        );
        d2vx.set_value(
            5,
            4.0 * (b2 - 2.0 * c) * (3.0 * beta2 * betaprim2 + beta3 * betasecn)
                + 16.0 * beta3 * betaprim * bpr * (b - 1.0)
                + 2.0 * beta4 * (bsc * (b - 1.0) + bpr2),
        );
        d2vx.set_value(
            7,
            -6.0 * c2 * (5.0 * beta4 * betaprim2 + beta5 * betasecn)
                - 24.0 * beta5 * betaprim * bpr * c
                - 2.0 * beta6 * (bsc * c + bpr2),
        );

        let mut d2vy = MathVector::with_init(1, ORDRE, 0.0);
        d2vy.set_value(2, 2.0 * betasecn);
        d2vy.set_value(
            4,
            6.0 * (b + c) * (2.0 * beta * betaprim2 + beta2 * betasecn)
                + 24.0 * beta2 * betaprim * bpr * (b + c)
                + 4.0 * beta3 * bsc,
        );
        d2vy.set_value(
            6,
            10.0 * b * c * (4.0 * beta3 * betaprim2 + beta4 * betasecn)
                + 40.0 * beta4 * betaprim * bpr * (b + c)
                + 2.0 * beta5 * (bsc * (b + c) + 2.0 * bpr2),
        );

        let mut d2vw = MathVector::with_init(1, ORDRE, 0.0);
        d2vw.set_value(
            3,
            2.0 * ((betaprim2 + beta * betasecn) * (2.0 * b + 1.0)
                + 8.0 * beta * betaprim * bpr
                + bsc * beta2),
        );
        d2vw.set_value(
            5,
            4.0 * (b2 + 2.0 * c) * (3.0 * beta2 * betaprim2 + beta3 * betasecn)
                + 16.0 * beta3 * betaprim * bpr * (b + 11.0)
                + 2.0 * beta4 * (bsc * (b + 1.0) + bpr2),
        );
        d2vw.set_value(
            7,
            6.0 * c2 * (5.0 * beta4 * betaprim2 + beta5 * betasecn)
                + 24.0 * beta5 * betaprim * bpr * c
                + 2.0 * beta6 * (bsc * c + bpr2),
        );

        let px = mat_vec(&self.b, &self.vx);
        let py = mat_vec(&self.b, &self.vy);
        let w = mat_vec(&self.b, &self.vw);
        let dpx = mat_vec(&self.b, &dvx);
        let dpy = mat_vec(&self.b, &dvy);
        let dw = mat_vec(&self.b, &dvw);
        let d2px = mat_vec(&self.b, &d2vx);
        let d2py = mat_vec(&self.b, &d2vy);
        let d2w = mat_vec(&self.b, &d2vw);

        let mut poles = Vec::with_capacity(ORDRE);
        let mut d_poles = Vec::with_capacity(ORDRE);
        let mut d2_poles = Vec::with_capacity(ORDRE);
        let mut weights = Vec::with_capacity(ORDRE);
        let mut d_weights = Vec::with_capacity(ORDRE);
        let mut d2_weights = Vec::with_capacity(ORDRE);
        for ii in 1..=ORDRE {
            let wi = w.value(ii);
            let dwi = dw.value(ii);
            let mut p = GpVec::new(px.value(ii) / wi, py.value(ii) / wi, 0.0);
            let mut dp = GpVec::new(dpx.value(ii) / wi, dpy.value(ii) / wi, 0.0);
            let mut d2p = GpVec::new(d2px.value(ii) / wi, d2py.value(ii) / wi, 0.0);

            d2p = d2p.added(&dp.multiplied_scalar(-2.0 * (dwi / wi)));
            let coef = 2.0 * (dwi / wi) * (dwi / wi) - d2w.value(ii) / wi;
            d2p = d2p.added(&p.multiplied_scalar(coef));
            dp = dp.added(&p.multiplied_scalar(-(dw.value(ii) / wi)));

            poles.push(GpPnt::from_xyz(&m.multiplied(&p.coord).added(&center.coord)));
            let a = mat_xyz(&m_prim, &p);
            let b2v = mat_xyz(&m, &dp);
            d_poles.push(GpVec::zero().set_linear_form_2plus(1.0, &a, 1.0, &b2v, d_center));
            p = mat_xyz(&m_secn, &p);
            dp = mat_xyz(&m_prim, &dp);
            d2p = mat_xyz(&m, &d2p);
            let aux = GpVec::zero()
                .set_linear_form_3(1.0, &p, 2.0, &dp, 1.0, &d2p)
                .added(d2_center);
            d2_poles.push(aux);
            weights.push(wi);
            d_weights.push(dwi);
            d2_weights.push(d2w.value(ii));
        }
        (poles, d_poles, d2_poles, weights, d_weights, d2_weights)
    }
}

/// `B` (`cxx:79-120`): the monomial-to-B-spline matrix on `[-1, 1]`, degree 6.
fn build_basis() -> MathMatrix {
    let mut coeffs = vec![0.0; ORDRE * ORDRE];
    for ii in 1..=ORDRE {
        coeffs[ii - 1 + (ii - 1) * ORDRE] = 1.0;
    }
    let intervals = [[-1.0, 1.0]];
    let true_intervals = [-1.0, 1.0];
    let converter = CompPolynomialToPoles::new(
        1,
        ORDRE as i32,
        (ORDRE - 1) as i32,
        &[],
        &[ORDRE as i32],
        &coeffs,
        &intervals,
        &true_intervals,
    )
    .ok();
    let poles = converter.map(|c| c.poles).unwrap_or_default();

    let mut result = MathMatrix::with_init(1, ORDRE, 1, ORDRE, 0.0);
    for jj in 1..=ORDRE {
        for ii in 1..=ORDRE {
            let raw = poles
                .get(ii - 1)
                .and_then(|row| row.get(jj - 1))
                .copied()
                .unwrap_or(0.0);
            let mut term = raw;
            if (term - 1.0).abs() < 1.0e-9 {
                term = 1.0;
            }
            if (term + 1.0).abs() < 1.0e-9 {
                term = -1.0;
            }
            result.set_value(ii, jj, term);
        }
    }
    result
}
