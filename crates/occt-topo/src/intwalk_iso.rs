//! `IntImp_ConstIsoparametric` and `IntImp_ComputeTangence`.

use occt_core::gp::GpVec;

/// Frozen isoparametric while solving `S1(u1,v1) = S2(u2,v2)`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConstIso {
    UOnS1,
    VOnS1,
    UOnS2,
    VOnS2,
}

impl ConstIso {
    pub fn from_index(i: i32) -> Self {
        match i & 3 {
            0 => Self::UOnS1,
            1 => Self::VOnS1,
            2 => Self::UOnS2,
            _ => Self::VOnS2,
        }
    }

    pub fn index(self) -> usize {
        match self {
            Self::UOnS1 => 0,
            Self::VOnS1 => 1,
            Self::UOnS2 => 2,
            Self::VOnS2 => 3,
        }
    }

    pub fn next(self) -> Self {
        Self::from_index((self.index() as i32) + 1)
    }
}

/// `IntWalk_StatusDeflection`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StatusDeflection {
    PasTropGrand,
    StepTooSmall,
    PointConfondu,
    ArretSurPointPrecedent,
    ArretSurPoint,
    Ok,
}

/// `IntImp_ComputeTangence`. Returns true when the surfaces are tangent.
/// `tgduv` is the 4 UV tangent components; `tab_iso` is the ranked iso choice.
pub fn compute_tangence(
    dpuv: &[GpVec; 4],
    eps_uv: &[f64; 4],
    tgduv: &mut [f64; 4],
    tab_iso: &mut [ConstIso; 4],
) -> bool {
    const A_TOL2: f64 = 1.0e-32;
    let mut norm = [0.0; 4];
    for i in 0..4 {
        norm[i] = dpuv[i].square_magnitude();
        if norm[i] <= A_TOL2 {
            return true;
        }
    }
    let mut n1 = dpuv[0].crossed(&dpuv[1]);
    if n1.square_magnitude() < A_TOL2 {
        return true;
    }
    n1.normalize();
    let mut n2 = dpuv[2].crossed(&dpuv[3]);
    if n2.square_magnitude() < A_TOL2 {
        return true;
    }
    n2.normalize();
    for i in 0..4 {
        norm[i] = norm[i].sqrt();
    }
    tgduv[0] = -dpuv[1].dot(&n2);
    tgduv[1] = dpuv[0].dot(&n2);
    tgduv[2] = dpuv[3].dot(&n1);
    tgduv[3] = -dpuv[2].dot(&n1);

    let mut tangent = tgduv[0].abs() <= eps_uv[0] * norm[1]
        && tgduv[1].abs() <= eps_uv[1] * norm[0]
        && tgduv[2].abs() <= eps_uv[2] * norm[3]
        && tgduv[3].abs() <= eps_uv[3] * norm[2];
    if !tangent {
        let t = n1.dot(&n2).abs();
        if t > 0.999999999 {
            tangent = true;
        }
    }
    if !tangent {
        norm[0] = tgduv[1].abs() / norm[0];
        norm[1] = tgduv[0].abs() / norm[1];
        norm[2] = tgduv[3].abs() / norm[2];
        norm[3] = tgduv[2].abs() / norm[3];
        for i in 0..4 {
            tab_iso[i] = ConstIso::from_index(i as i32);
        }
        let mut tri_ok = false;
        while !tri_ok {
            tri_ok = true;
            for i in 1..4 {
                if norm[i - 1] > norm[i] {
                    tri_ok = false;
                    norm.swap(i - 1, i);
                    tab_iso.swap(i - 1, i);
                }
            }
        }
    }
    tangent
}
