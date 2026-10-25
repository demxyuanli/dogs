//! `Law_BSplineKnotSplitting`.
//!
//! Source: `Law_BSplineKnotSplitting.hxx`/`Law_BSplineKnotSplitting.cxx:24-114`
//! (`TKGeomAlgo/Law`).
//!
//! Splits a `Law_BSpline` into the sub-intervals on which its continuity is at
//! least `ContinuityRange`: a knot belongs to the split only when
//! `Degree - Multiplicity < ContinuityRange` (`cxx:56-75`). The values are
//! unique-knot indices (`FirstUKnotIndex()..LastUKnotIndex()` of the curve),
//! not parameters and not 0-based offsets.

use crate::law::bspline::LawBSpline;

/// `Law_BSplineKnotSplitting` (`Law_BSplineKnotSplitting.hxx:34-64`).
#[derive(Debug, Clone, Default)]
pub struct LawBSplineKnotSplitting {
    split_indexes: Vec<i32>,
}

impl LawBSplineKnotSplitting {
    /// `Law_BSplineKnotSplitting::Law_BSplineKnotSplitting` (`cxx:27-81`).
    pub fn new(curve: &LawBSpline, continuity_range: i32) -> Result<Self, String> {
        if continuity_range < 0 {
            return Err("Standard_RangeError: Law_BSplineKnotSplitting".into());
        }
        let first_index = curve.first_u_knot_index();
        let last_index = curve.last_u_knot_index();
        let degree = curve.degree();
        if continuity_range == 0 {
            return Ok(Self {
                split_indexes: vec![first_index, last_index],
            });
        }
        let mults = curve.multiplicities();
        // `BSplCLib::MaxKnotMult(Mults, FirstIndex, LastIndex)` (`cxx:49`).
        let mmax = mults
            .get((first_index - 1) as usize..last_index as usize)
            .map(|s| s.iter().copied().max().unwrap_or(0))
            .unwrap_or(0);
        if degree - mmax >= continuity_range {
            return Ok(Self {
                split_indexes: vec![first_index, last_index],
            });
        }
        let mut split = vec![first_index];
        let mut index = first_index + 1;
        while index < last_index {
            let mult = mults.get((index - 1) as usize).copied().unwrap_or(0);
            if degree - mult < continuity_range {
                split.push(index);
            }
            index += 1;
        }
        split.push(index);
        Ok(Self {
            split_indexes: split,
        })
    }

    /// `Law_BSplineKnotSplitting::NbSplits` (`cxx:83-87`).
    pub fn nb_splits(&self) -> usize {
        self.split_indexes.len()
    }

    /// `Law_BSplineKnotSplitting::SplitValue` (`cxx:89-98`), 1-based.
    pub fn split_value(&self, index: i32) -> Result<i32, String> {
        if index < 1 || index > self.split_indexes.len() as i32 {
            return Err("Standard_RangeError: Law_BSplineKnotSplitting::SplitValue".into());
        }
        Ok(self.split_indexes[(index - 1) as usize])
    }

    /// `Law_BSplineKnotSplitting::Splitting` (`cxx:100-113`).
    pub fn splitting(&self) -> &[i32] {
        &self.split_indexes
    }
}
