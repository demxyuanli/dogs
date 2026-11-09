//! Euler-angle sequence enumeration and parameter translation.
//! Source: `gp_EulerSequence.hxx`, `gp_Quaternion.cxx:191-296`.

/// Enumerates all 24 possible variants of generalized Euler angles.
/// Order matches `gp_EulerSequence` (`gp_EulerSequence.hxx:37-73`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GpEulerSequence {
    EulerAngles,
    YawPitchRoll,
    ExtrinsicXYZ,
    ExtrinsicXZY,
    ExtrinsicYZX,
    ExtrinsicYXZ,
    ExtrinsicZXY,
    ExtrinsicZYX,
    IntrinsicXYZ,
    IntrinsicXZY,
    IntrinsicYZX,
    IntrinsicYXZ,
    IntrinsicZXY,
    IntrinsicZYX,
    ExtrinsicXYX,
    ExtrinsicXZX,
    ExtrinsicYZY,
    ExtrinsicYXY,
    ExtrinsicZYZ,
    ExtrinsicZXZ,
    IntrinsicXYX,
    IntrinsicXZX,
    IntrinsicYZY,
    IntrinsicYXY,
    IntrinsicZXZ,
    IntrinsicZYZ,
}

/// `gp_EulerSequence_Parameters` (`gp_Quaternion.cxx:202-221`).
#[derive(Debug, Clone, Copy)]
pub struct EulerSequenceParams {
    /// First rotation axis (1-based: x=1, y=2, z=3).
    pub i: usize,
    /// Next axis of rotation.
    pub j: usize,
    /// Third axis.
    pub k: usize,
    /// True if the order of the two first rotation axes is an odd permutation.
    pub is_odd: bool,
    /// True if the third rotation is about the same axis as the first.
    pub is_two_axes: bool,
    /// True if rotations are made around fixed axes.
    pub is_extrinsic: bool,
}

impl EulerSequenceParams {
    fn new(ax1: usize, is_odd: bool, is_two_axes: bool, is_extrinsic: bool) -> Self {
        let j = 1 + (ax1 + if is_odd { 1 } else { 0 }) % 3;
        let k = 1 + (ax1 + if is_odd { 0 } else { 1 }) % 3;
        Self {
            i: ax1,
            j,
            k,
            is_odd,
            is_two_axes,
            is_extrinsic,
        }
    }
}

/// `translateEulerSequence` (`gp_Quaternion.cxx:226-296`).
pub fn translate_euler_sequence(seq: GpEulerSequence) -> EulerSequenceParams {
    use GpEulerSequence::*;
    let f = false;
    let t = true;
    match seq {
        ExtrinsicXYZ => EulerSequenceParams::new(1, f, f, t),
        ExtrinsicXZY => EulerSequenceParams::new(1, t, f, t),
        ExtrinsicYZX => EulerSequenceParams::new(2, f, f, t),
        ExtrinsicYXZ => EulerSequenceParams::new(2, t, f, t),
        ExtrinsicZXY => EulerSequenceParams::new(3, f, f, t),
        ExtrinsicZYX => EulerSequenceParams::new(3, t, f, t),

        IntrinsicXYZ => EulerSequenceParams::new(3, t, f, f),
        IntrinsicXZY => EulerSequenceParams::new(2, f, f, f),
        IntrinsicYZX => EulerSequenceParams::new(1, t, f, f),
        IntrinsicYXZ => EulerSequenceParams::new(3, f, f, f),
        IntrinsicZXY => EulerSequenceParams::new(2, t, f, f),
        IntrinsicZYX => EulerSequenceParams::new(1, f, f, f),

        ExtrinsicXYX => EulerSequenceParams::new(1, f, t, t),
        ExtrinsicXZX => EulerSequenceParams::new(1, t, t, t),
        ExtrinsicYZY => EulerSequenceParams::new(2, f, t, t),
        ExtrinsicYXY => EulerSequenceParams::new(2, t, t, t),
        ExtrinsicZXZ => EulerSequenceParams::new(3, f, t, t),
        ExtrinsicZYZ => EulerSequenceParams::new(3, t, t, t),

        IntrinsicXYX => EulerSequenceParams::new(1, f, t, f),
        IntrinsicXZX => EulerSequenceParams::new(1, t, t, f),
        IntrinsicYZY => EulerSequenceParams::new(2, f, t, f),
        IntrinsicYXY => EulerSequenceParams::new(2, t, t, f),
        IntrinsicZXZ => EulerSequenceParams::new(3, f, t, f),
        IntrinsicZYZ => EulerSequenceParams::new(3, t, t, f),

        EulerAngles => EulerSequenceParams::new(3, f, t, f),
        YawPitchRoll => EulerSequenceParams::new(1, f, f, f),
    }
}
