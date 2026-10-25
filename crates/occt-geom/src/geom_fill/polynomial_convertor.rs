//! `GeomFill_PolynomialConvertor` — the degree-7 polynomial approximation of a
//! circular arc section, used by the polynomial sweep (`GeomFill::GetCircle`
//! with `Convert_Polynomial`).
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_PolynomialConvertor.cxx:29-366` and
//! `GeomFill_PolynomialConvertor.hxx` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `math_Matrix BH` (`hxx`) becomes [`MathMatrix`]; `Initialized()` /
//!   `Init()` (`cxx:34-96`) keep OCCT's lazy build, including the `1e-9`
//!   snapping of the monomial-to-B-spline terms to `+1` / `-1` (`cxx:75-87`).
//! - `Convert_CompPolynomialToPoles(anOrdre, anOrdre - 1, anOrdre - 1,
//!   aCoeffs, anInter, aTrueInter)` (`cxx:61-68`) is OCCT's single-span
//!   constructor (`Convert_CompPolynomialToPoles.cxx:140-175`), expressed here
//!   through the one-curve form of [`CompPolynomialToPoles::new`]; the
//!   `Dimension = Ordre` basis order is what puts the eight monomials in the
//!   columns of `B`.
//! - `PLib::HermiteCoefficients(-1, 1, 3, 3, aH)` (`cxx:90`) is
//!   [`hermite_coefficients`]. OCCT ignores the returned `bool`, so a `None`
//!   leaves `aH` at the zero matrix it was constructed with.
//! - `Poles` / `DPoles` / `D2Poles` are `NCollection_Array1` in OCCT; the port
//!   uses slices whose length is `Ordre`.
//! - `gp_Mat M(V1.X(), V2.X(), 0, ...)` (`cxx:117`) is row-major in OCCT, which
//!   is the same matrix as the two-column form built below.

use occt_core::bspl::comp_poly::CompPolynomialToPoles;
use occt_core::gp::{GpMat, GpPnt, GpVec, GpXyz};
use occt_math::{hermite_coefficients, MathMatrix, MathVector};

/// `GeomFill_PolynomialConvertor::Ordre` — the section has `Ordre` poles
/// (`cxx:30`).
pub const ORDRE: usize = 8;

/// `math_Matrix * math_Vector` (`math_Matrix::Multiply(V, Result)`).
fn mat_vec(m: &MathMatrix, v: &MathVector) -> MathVector {
    let mut out = MathVector::with_init(1, m.row_count(), 0.0);
    for i in 1..=m.row_count() {
        let mut sum = 0.0;
        for j in 1..=m.col_count() {
            sum += m.value(i, j) * v.value(j);
        }
        out.set_value(i, sum);
    }
    out
}

/// `gp_Mat M(V1.X(), V2.X(), 0, V1.Y(), V2.Y(), 0, V1.Z(), V2.Z(), 0)`
/// (`cxx:117`): the first two columns are `V1` and `V2`.
fn frame(v1: &GpVec, v2: &GpVec) -> GpMat {
    GpMat::from_cols(&v1.coord, &v2.coord, &GpXyz::zero())
}

/// `GeomFill_PolynomialConvertor` (`GeomFill_PolynomialConvertor.hxx`).
pub struct PolynomialConvertor {
    /// `BH` — `B * H`, built once by [`PolynomialConvertor::init`].
    bh: MathMatrix,
    /// `myinit`.
    initialized: bool,
}

impl Default for PolynomialConvertor {
    fn default() -> Self {
        Self::new()
    }
}

impl PolynomialConvertor {
    /// `GeomFill_PolynomialConvertor()` (`cxx:29-33`).
    pub fn new() -> Self {
        Self {
            bh: MathMatrix::with_init(1, ORDRE, 1, ORDRE, 0.0),
            initialized: false,
        }
    }

    /// `Initialized()` (`cxx:34-37`).
    pub fn initialized(&self) -> bool {
        self.initialized
    }

    /// `Init()` (`cxx:39-96`).
    pub fn init(&mut self) {
        if self.initialized {
            return;
        }
        self.bh = build_bh();
        self.initialized = true;
    }

    /// `Section(FirstPnt, Center, Dir, Angle, Poles)` (`cxx:98-148`).
    pub fn section(&mut self, first_pnt: &GpPnt, center: &GpPnt, dir: &GpVec, angle: f64) -> Vec<GpPnt> {
        self.init();
        let mut vx = MathVector::with_init(1, ORDRE, 0.0);
        let mut vy = MathVector::with_init(1, ORDRE, 0.0);
        let cos_b = angle.cos();
        let sin_b = angle.sin();
        let v1 = GpVec::from_pnts(center, first_pnt);
        let v2 = dir.crossed(&v1);
        let beta = angle / 2.0;
        let beta2 = beta * beta;
        let beta3 = beta * beta2;

        let m = frame(&v1, &v2);
        fill_constraints(&mut vx, &mut vy, cos_b, sin_b, beta, beta2, beta3);

        let px = mat_vec(&self.bh, &vx);
        let py = mat_vec(&self.bh, &vy);
        let mut poles = Vec::with_capacity(ORDRE);
        for ii in 1..=ORDRE {
            let mut pnt = GpXyz::new(px.value(ii), py.value(ii), 0.0);
            pnt.multiply_mat(&m);
            poles.push(GpPnt::from_xyz(&pnt.added(&center.coord)));
        }
        poles
    }

    /// `Section(FirstPnt, DFirstPnt, Center, DCenter, Dir, DDir, Angle,
    /// DAngle, Poles, DPoles)` (`cxx:150-233`).
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
    ) -> (Vec<GpPnt>, Vec<GpVec>) {
        self.init();
        let mut vx = MathVector::with_init(1, ORDRE, 0.0);
        let mut vy = MathVector::with_init(1, ORDRE, 0.0);
        let mut dvx = MathVector::with_init(1, ORDRE, 0.0);
        let mut dvy = MathVector::with_init(1, ORDRE, 0.0);
        let cos_b = angle.cos();
        let sin_b = angle.sin();
        let v1 = GpVec::from_pnts(center, first_pnt);
        let beta = angle / 2.0;
        let bprim = d_angle / 2.0;
        let beta2 = beta * beta;
        let beta3 = beta * beta2;

        // Calcul des transformations.
        let m = frame(&v1, &dir.crossed(&v1));
        let v1_prim = d_first_pnt.subtracted(d_center);
        let v2 = d_dir.crossed(&v1).added(&dir.crossed(&v1_prim));
        let m_prim = frame(&v1_prim, &v2);

        fill_constraints(&mut vx, &mut vy, cos_b, sin_b, beta, beta2, beta3);

        let b_bprim = bprim * beta;
        let b2_bprim = bprim * beta2;
        dvx.set_value(2, 0.0);
        dvy.set_value(2, bprim);
        dvx.set_value(3, -2.0 * b_bprim);
        dvy.set_value(4, -3.0 * b2_bprim);
        dvx.set_value(5, -2.0 * bprim * sin_b);
        dvy.set_value(5, 2.0 * bprim * cos_b);
        dvx.set_value(6, -bprim * sin_b - 2.0 * b_bprim * cos_b);
        dvy.set_value(6, bprim * cos_b - 2.0 * b_bprim * sin_b);
        dvx.set_value(7, 2.0 * b_bprim * (-cos_b + beta * sin_b));
        dvy.set_value(7, -2.0 * b_bprim * (sin_b + beta * cos_b));
        dvx.set_value(8, b2_bprim * (3.0 * sin_b + 2.0 * beta * cos_b));
        dvy.set_value(8, b2_bprim * (2.0 * beta * sin_b - 3.0 * cos_b));

        let px = mat_vec(&self.bh, &vx);
        let py = mat_vec(&self.bh, &vy);
        let dpx = mat_vec(&self.bh, &dvx);
        let dpy = mat_vec(&self.bh, &dvy);

        let mut poles = Vec::with_capacity(ORDRE);
        let mut d_poles = Vec::with_capacity(ORDRE);
        for ii in 1..=ORDRE {
            let p = GpVec::new(px.value(ii), py.value(ii), 0.0);
            poles.push(GpPnt::from_xyz(&m.multiplied(&p.coord).added(&center.coord)));
            let p_prim = m_prim.multiplied(&p.coord);
            let dp = GpVec::new(dpx.value(ii), dpy.value(ii), 0.0);
            let dp_m = m.multiplied(&dp.coord);
            d_poles.push(
                GpVec::zero().set_linear_form_2plus(1.0, &GpVec::from_xyz(&p_prim), 1.0, &GpVec::from_xyz(&dp_m), d_center),
            );
        }
        (poles, d_poles)
    }

    /// `Section(FirstPnt, DFirstPnt, D2FirstPnt, Center, DCenter, D2Center,
    /// Dir, DDir, D2Dir, Angle, DAngle, D2Angle, Poles, DPoles, D2Poles)`
    /// (`cxx:235-366`).
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
    ) -> (Vec<GpPnt>, Vec<GpVec>, Vec<GpVec>) {
        self.init();
        let mut vx = MathVector::with_init(1, ORDRE, 0.0);
        let mut vy = MathVector::with_init(1, ORDRE, 0.0);
        let mut dvx = MathVector::with_init(1, ORDRE, 0.0);
        let mut dvy = MathVector::with_init(1, ORDRE, 0.0);
        let mut d2vx = MathVector::with_init(1, ORDRE, 0.0);
        let mut d2vy = MathVector::with_init(1, ORDRE, 0.0);
        let cos_b = angle.cos();
        let sin_b = angle.sin();
        let v1 = GpVec::from_pnts(center, first_pnt);
        let beta = angle / 2.0;
        let bprim = d_angle / 2.0;
        let bsecn = d2_angle / 2.0;
        let beta2 = beta * beta;
        let beta3 = beta * beta2;
        let bprim2 = bprim * bprim;

        // Calcul des transformations.
        let m = frame(&v1, &dir.crossed(&v1));
        let v1_prim = d_first_pnt.subtracted(d_center);
        let v2_prim = d_dir.crossed(&v1).added(&dir.crossed(&v1_prim));
        let m_prim = frame(&v1_prim, &v2_prim);
        let v1_secn = d2_first_pnt.subtracted(d2_center);
        let mut v2_secn = d_dir.crossed(&v1_prim).multiply_scalar(2.0);
        v2_secn = v2_secn
            .added(&d2_dir.crossed(&v1))
            .added(&dir.crossed(&v1_secn));
        let m_secn = frame(&v1_secn, &v2_secn);

        fill_constraints(&mut vx, &mut vy, cos_b, sin_b, beta, beta2, beta3);

        let b_bprim = bprim * beta;
        let b2_bprim = bprim * beta2;
        let b_bsecn = bsecn * beta;
        dvy.set_value(2, bprim);
        dvx.set_value(3, -2.0 * b_bprim);
        dvy.set_value(4, -3.0 * b2_bprim);
        dvx.set_value(5, -2.0 * bprim * sin_b);
        dvy.set_value(5, 2.0 * bprim * cos_b);
        dvx.set_value(6, -bprim * sin_b - 2.0 * b_bprim * cos_b);
        dvy.set_value(6, bprim * cos_b - 2.0 * b_bprim * sin_b);
        dvx.set_value(7, 2.0 * b_bprim * (-cos_b + beta * sin_b));
        dvy.set_value(7, -2.0 * b_bprim * (sin_b + beta * cos_b));
        dvx.set_value(8, b2_bprim * (3.0 * sin_b + 2.0 * beta * cos_b));
        dvy.set_value(8, b2_bprim * (2.0 * beta * sin_b - 3.0 * cos_b));

        d2vy.set_value(2, bsecn);
        d2vx.set_value(3, -2.0 * (bprim2 + b_bsecn));
        d2vy.set_value(4, -3.0 * beta * (2.0 * bprim2 + b_bsecn));
        d2vx.set_value(5, -2.0 * (bsecn * sin_b + 2.0 * bprim2 * cos_b));
        d2vy.set_value(5, 2.0 * (bsecn * cos_b - 2.0 * bprim2 * sin_b));
        d2vx.set_value(6, (4.0 * beta * bprim2 - bsecn) * sin_b - 2.0 * (2.0 * bprim2 + b_bsecn) * cos_b);
        d2vy.set_value(6, (bsecn - 4.0 * beta * bprim2) * cos_b - 2.0 * (b_bsecn + 2.0 * bprim2) * sin_b);
        let aux = 2.0 * (bprim2 + b_bsecn);
        d2vx.set_value(
            7,
            aux * (-cos_b + beta * sin_b) + 2.0 * beta * bprim2 * (2.0 * beta * cos_b + 3.0 * sin_b),
        );
        d2vy.set_value(
            7,
            -aux * (sin_b + beta * cos_b) - 2.0 * beta * bprim2 * (3.0 * cos_b - 2.0 * beta * sin_b),
        );
        let aux = beta * (2.0 * bprim2 + b_bsecn);
        d2vx.set_value(
            8,
            aux * (3.0 * sin_b + 2.0 * beta * cos_b) + 4.0 * beta2 * bprim2 * (2.0 * cos_b - beta * sin_b),
        );
        d2vy.set_value(
            8,
            aux * (2.0 * beta * sin_b - 3.0 * cos_b) + 4.0 * beta2 * bprim2 * (2.0 * sin_b + beta * cos_b),
        );

        let px = mat_vec(&self.bh, &vx);
        let py = mat_vec(&self.bh, &vy);
        let dpx = mat_vec(&self.bh, &dvx);
        let dpy = mat_vec(&self.bh, &dvy);
        let d2px = mat_vec(&self.bh, &d2vx);
        let d2py = mat_vec(&self.bh, &d2vy);

        let mut poles = Vec::with_capacity(ORDRE);
        let mut d_poles = Vec::with_capacity(ORDRE);
        let mut d2_poles = Vec::with_capacity(ORDRE);
        for ii in 1..=ORDRE {
            let mut p = GpVec::new(px.value(ii), py.value(ii), 0.0);
            let mut dp = GpVec::new(dpx.value(ii), dpy.value(ii), 0.0);
            let mut d2p = GpVec::new(d2px.value(ii), d2py.value(ii), 0.0);

            poles.push(GpPnt::from_xyz(&m.multiplied(&p.coord).added(&center.coord)));
            let a = m_prim.multiplied(&p.coord);
            let b = m.multiplied(&dp.coord);
            d_poles.push(GpVec::zero().set_linear_form_2plus(
                1.0,
                &GpVec::from_xyz(&a),
                1.0,
                &GpVec::from_xyz(&b),
                d_center,
            ));
            p = GpVec::from_xyz(&m_secn.multiplied(&p.coord));
            dp = GpVec::from_xyz(&m_prim.multiplied(&dp.coord));
            d2p = GpVec::from_xyz(&m.multiplied(&d2p.coord));
            d2_poles.push(GpVec::zero().set_linear_form_3(
                1.0,
                &p,
                2.0,
                &dp,
                1.0,
                &d2p,
            ).added(d2_center));
        }
        (poles, d_poles, d2_poles)
    }
}

/// The `Vx` / `Vy` constraint vectors shared by all three `Section` overloads
/// (`cxx:120-135`, `cxx:180-196`, `cxx:282-296`).
fn fill_constraints(
    vx: &mut MathVector,
    vy: &mut MathVector,
    cos_b: f64,
    sin_b: f64,
    beta: f64,
    beta2: f64,
    beta3: f64,
) {
    vx.set_value(1, 1.0);
    vy.set_value(1, 0.0);
    vx.set_value(2, 0.0);
    vy.set_value(2, beta);
    vx.set_value(3, -beta2);
    vy.set_value(3, 0.0);
    vx.set_value(4, 0.0);
    vy.set_value(4, -beta3);
    vx.set_value(5, cos_b);
    vy.set_value(5, sin_b);
    vx.set_value(6, -beta * sin_b);
    vy.set_value(6, beta * cos_b);
    vx.set_value(7, -beta2 * cos_b);
    vy.set_value(7, -beta2 * sin_b);
    vx.set_value(8, beta3 * sin_b);
    vy.set_value(8, -beta3 * cos_b);
}

/// `BH = B * H` (`cxx:50-94`): the monomial-to-B-spline matrix `B` on `[-1, 1]`
/// times the transposed Hermite coefficient matrix.
fn build_bh() -> MathMatrix {
    let ordre = ORDRE;
    // Monomial basis: degree `ordre - 1`, `ordre` coefficients in `ordre`
    // dimensions, each polynomial being a different coordinate axis.
    let mut coeffs = vec![0.0; ordre * ordre];
    for ii in 1..=ordre {
        coeffs[ii - 1 + (ii - 1) * ordre] = 1.0;
    }
    let intervals = [[-1.0, 1.0]];
    let true_intervals = [-1.0, 1.0];
    let converter = CompPolynomialToPoles::new(
        1,
        ordre as i32,
        (ordre - 1) as i32,
        &[],
        &[ordre as i32],
        &coeffs,
        &intervals,
        &true_intervals,
    )
    .ok();
    let empty: Vec<Vec<f64>> = Vec::new();
    let poles = converter.map(|c| c.poles).unwrap_or(empty);

    let mut b = MathMatrix::with_init(1, ordre, 1, ordre, 0.0);
    for jj in 1..=ordre {
        for ii in 1..=ordre {
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
            b.set_value(ii, jj, term);
        }
    }

    // PLib::HermiteCoefficients(-1, 1, Ordre/2-1, Ordre/2-1, aH); aH.Transpose()
    let mut h = MathMatrix::with_init(1, ordre, 1, ordre, 0.0);
    if let Some(rows) = hermite_coefficients(-1.0, 1.0, ordre as i32 / 2 - 1, ordre as i32 / 2 - 1) {
        for (i, row) in rows.iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                h.set_value(i + 1, j + 1, *value);
            }
        }
    }
    h.transpose();
    b.multiplied_mat(&h)
}
