//! Static pruning helpers of the `HLRAlgo` facade. Source: `HLRAlgo`
//! (`ModelingAlgorithms/TKHLR/HLRAlgo/HLRAlgo.hxx:44-88`,
//! `HLRAlgo.cxx:19-205`).
//!
//! The 16-entry `Min` / `Max` arrays are the 14 directions
//! `cos(k*pi/14) / sin(k*pi/14)` (`k = 0..6`) of the screen plane plus `z`
//! twice (`cxx:22-34`, `:44-58`).

use super::edges_block::MinMaxIndices;

use std::f64::consts::PI;

/// `cosu0..cosu6` (`cxx:22-34`): `cos(k * M_PI / 14)`. OCCT evaluates these
/// once at load time; the same expression is evaluated per call here so the
/// values are bit-identical.
fn cosu(k: usize) -> f64 {
    f64::cos(k as f64 * PI / 14.0)
}

/// `sinu0..sinu6` (`cxx:23-35`).
fn sinu(k: usize) -> f64 {
    f64::sin(k as f64 * PI / 14.0)
}

/// `HLRAlgo::UpdateMinMax` (`cxx:36-72`): grow `min` / `max` with the 16
/// projections of `(x, y, z)`.
pub fn update_min_max(x: f64, y: f64, z: f64, min: &mut [f64; 16], max: &mut [f64; 16]) {
    let mut d = [0.0f64; 16];
    for k in 0..7 {
        let (c, s) = (cosu(k), sinu(k));
        d[2 * k] = c * x + s * y;
        d[2 * k + 1] = s * x - c * y;
    }
    d[14] = z;
    d[15] = z;
    for i in 0..16 {
        if min[i] > d[i] {
            min[i] = d[i];
        }
        if max[i] < d[i] {
            max[i] = d[i];
        }
    }
}

/// `HLRAlgo::EnlargeMinMax` (`cxx:76-85`).
pub fn enlarge_min_max(tol: f64, min: &mut [f64; 16], max: &mut [f64; 16]) {
    for i in 0..16 {
        min[i] -= tol;
        max[i] += tol;
    }
}

/// `HLRAlgo::InitMinMax` (`cxx:89-97`): `min = Big`, `max = -Big`.
pub fn init_min_max(big: f64, min: &mut [f64; 16], max: &mut [f64; 16]) {
    for i in 0..16 {
        min[i] = big;
        max[i] = -big;
    }
}

/// `HLRAlgo::EncodeMinMax` (`cxx:102-138`): pack the 16 `int` bounds into 8
/// `int` slots, two 15-bit fields per slot.
pub fn encode_min_max(min: &MinMaxIndices, max: &MinMaxIndices, mm: &mut MinMaxIndices) {
    const LO: i32 = 0x0000_7fff;
    for k in 0..4 {
        mm.min[k] = (min.min[2 * k + 1] & LO) + ((min.min[2 * k] & LO) << 16);
        mm.max[k] = (max.min[2 * k + 1] & LO) + ((max.min[2 * k] & LO) << 16);
        mm.min[k + 4] = (min.max[2 * k + 1] & LO) + ((min.max[2 * k] & LO) << 16);
        mm.max[k + 4] = (max.max[2 * k + 1] & LO) + ((max.max[2 * k] & LO) << 16);
    }
}

/// `HLRAlgo::SizeBox` (`cxx:142-155`): the product of the eight `Min`
/// extents and the first six `Max` extents.
pub fn size_box(min: &MinMaxIndices, max: &MinMaxIndices) -> f64 {
    let mut s = (max.min[0] - min.min[0]) as f64;
    for i in 1..8 {
        s *= (max.min[i] - min.min[i]) as f64;
    }
    for i in 0..6 {
        s *= (max.max[i] - min.max[i]) as f64;
    }
    s
}

/// `HLRAlgo::DecodeMinMax` (`cxx:159-195`): the inverse of
/// [`encode_min_max`].
pub fn decode_min_max(mm: &MinMaxIndices, min: &mut MinMaxIndices, max: &mut MinMaxIndices) {
    const LO: i32 = 0x0000_7fff;
    const HI: i32 = 0x7fff_0000;
    for k in 0..4 {
        min.min[2 * k] = (mm.min[k] & HI) >> 16;
        max.min[2 * k] = (mm.max[k] & HI) >> 16;
        min.min[2 * k + 1] = mm.min[k] & LO;
        max.min[2 * k + 1] = mm.max[k] & LO;
        min.max[2 * k] = (mm.min[k + 4] & HI) >> 16;
        max.max[2 * k] = (mm.max[k + 4] & HI) >> 16;
        min.max[2 * k + 1] = mm.min[k + 4] & LO;
        max.max[2 * k + 1] = mm.max[k + 4] & LO;
    }
}

/// `HLRAlgo::CopyMinMax` (`hxx:72-79`): returns `(OMin, OMax)`.
pub fn copy_min_max(imin: &MinMaxIndices, imax: &MinMaxIndices) -> (MinMaxIndices, MinMaxIndices) {
    (*imin, *imax)
}

/// `HLRAlgo::AddMinMax` (`cxx:198-205`): `OMin.Minimize(IMin)`,
/// `OMax.Maximize(IMax)`.
pub fn add_min_max(
    imin: &MinMaxIndices,
    imax: &MinMaxIndices,
    omin: &mut MinMaxIndices,
    omax: &mut MinMaxIndices,
) {
    omin.minimize(imin);
    omax.maximize(imax);
}
