//! The shared interval merge of `BRepBlend_RstRstEvolRad`.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_RstRstEvolRad.cxx:54-113`
//! (`FusionneIntervalles`) (OCCT 8.0.0).

use occt_core::precision::PCONFUSION;

/// `FusionneIntervalles(I1, I2, Seq)` (`cxx:54-113`): merge two ascending
/// interval tables, coalescing parameters closer than
/// `PConfusion() * 0.99` into their midpoint.
///
/// OCCT fills a `NCollection_Sequence` 1-based; the port appends to a `Vec`.
pub(super) fn fusionne_intervalles(i1: &[f64], i2: &[f64], seq: &mut Vec<f64>) {
    let epspar = PCONFUSION * super::constants::FUSIONNE_EPS_FACTOR;
    let mut ind1 = 0usize;
    let mut ind2 = 0usize;
    while ind1 < i1.len() && ind2 < i2.len() {
        let v1 = i1[ind1];
        let v2 = i2[ind2];
        if (v1 - v2).abs() <= epspar {
            // Here the elements of I1 and I2 fit (`cxx:73-80`).
            seq.push((v1 + v2) / 2.0);
            ind1 += 1;
            ind2 += 1;
        } else if v1 < v2 {
            seq.push(v1);
            ind1 += 1;
        } else {
            seq.push(v2);
            ind2 += 1;
        }
    }
    while ind1 < i1.len() {
        seq.push(i1[ind1]);
        ind1 += 1;
    }
    while ind2 < i2.len() {
        seq.push(i2[ind2]);
        ind2 += 1;
    }
}