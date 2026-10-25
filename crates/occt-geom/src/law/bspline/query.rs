//! `Law_BSpline` property readers and `LocateU`.
//!
//! Source: `Law_BSpline.cxx:1145-1156` (`PeriodicNormalization`), `:1160-1240`
//! (`IsCN` / `IsClosed` / `IsPeriodic` / `Continuity` / `Degree`),
//! `:1245-1315` (`D0`/`D1`/`D2`/`D3`/`DN`), `:1320-1340` (`EndPoint`),
//! `:1345-1420` (knot / index / parameter readers), `:1535-1620`
//! (`Multiplicity` / `NKnots` / `NPoles` / `Pole` / `Weight` / `IsRational`) and
//! `:1625-1700` (`LocateU`), all in OCCT 8.0.0.
//!
//! Translation notes: see the module header of [`super`].

use occt_core::bspl::curve_1d::{curve_1d_d0, curve_1d_d1, curve_1d_d2, curve_1d_d3, curve_1d_dn};
use occt_core::bspl::locate::{
    first_u_knot_index, flat_index, hunt_occt, last_u_knot_index, locate_parameter_flat,
};
use occt_core::kernel::geomabs::{BSplKnotDistribution, Shape};
use occt_core::precision::REAL_SMALL;

use super::LawBSpline;

impl LawBSpline {
    pub fn periodic_normalization(&self, parameter: f64) -> f64 {
        let mut p = parameter;
        if self.periodic {
            let upper = self.flat_knots[self.flat_knots.len() - self.deg as usize - 1];
            let first = self.flat_knots[self.deg as usize];
            let period = upper - first;
            while p > upper {
                p -= period;
            }
            while p < first {
                p += period;
            }
        }
        p
    }

    /// `Law_BSpline::IsCN` (`Law_BSpline.cxx:1160-1187`).
    pub fn is_cn(&self, n: i32) -> Result<bool, String> {
        if n < 0 {
            return Err("Standard_RangeError: Law_BSpline::IsCN".into());
        }
        Ok(match self.smooth {
            Shape::CN => true,
            Shape::C0 | Shape::G1 => n <= 0,
            Shape::C1 | Shape::G2 => n <= 1,
            Shape::C2 => n <= 2,
            Shape::C3 => {
                if n <= 3 {
                    true
                } else {
                    // `BSplCLib::MaxKnotMult(mults, lower + 1, upper - 1)`
                    // (`BSplCLib.cxx:368-376`).
                    let mmax = self
                        .mults
                        .get(1..self.mults.len().saturating_sub(1))
                        .map(|s| s.iter().copied().max().unwrap_or(0))
                        .unwrap_or(0);
                    n <= self.deg - mmax
                }
            }
        })
    }

    /// `Law_BSpline::IsClosed` (`Law_BSpline.cxx:1191-1194`).
    pub fn is_closed(&self) -> bool {
        (self.start_point() - self.end_point()).abs() <= REAL_SMALL
    }

    /// `Law_BSpline::IsPeriodic` (`Law_BSpline.cxx:1198-1201`).
    pub fn is_periodic(&self) -> bool {
        self.periodic
    }

    /// `Law_BSpline::Continuity` (`Law_BSpline.cxx:1205-1208`).
    pub fn continuity(&self) -> Shape {
        self.smooth
    }

    /// `Law_BSpline::Degree` (`Law_BSpline.cxx:1212-1215`).
    pub fn degree(&self) -> i32 {
        self.deg
    }

    /// `Law_BSpline::Value` (`Law_BSpline.cxx:1219-1223`).
    pub fn value(&self, u: f64) -> f64 {
        self.d0(u)
    }

    /// `Law_BSpline::D0` (`Law_BSpline.cxx:1227-1240`).
    pub fn d0(&self, u: f64) -> f64 {
        let new_u = self.periodic_normalization(u);
        curve_1d_d0(
            new_u,
            0,
            self.deg,
            self.periodic,
            &self.poles,
            self.weights.as_deref(),
            &self.flat_knots,
            None,
        )
    }

    /// `Law_BSpline::D1` (`Law_BSpline.cxx:1244-1257`).
    pub fn d1(&self, u: f64) -> (f64, f64) {
        let new_u = self.periodic_normalization(u);
        curve_1d_d1(
            new_u,
            0,
            self.deg,
            self.periodic,
            &self.poles,
            self.weights.as_deref(),
            &self.flat_knots,
            None,
        )
    }

    /// `Law_BSpline::D2` (`Law_BSpline.cxx:1261-1272`).
    pub fn d2(&self, u: f64) -> (f64, f64, f64) {
        let new_u = self.periodic_normalization(u);
        curve_1d_d2(
            new_u,
            0,
            self.deg,
            self.periodic,
            &self.poles,
            self.weights.as_deref(),
            &self.flat_knots,
            None,
        )
    }

    /// `Law_BSpline::D3` (`Law_BSpline.cxx:1276-1299`).
    pub fn d3(&self, u: f64) -> (f64, f64, f64, f64) {
        let new_u = self.periodic_normalization(u);
        curve_1d_d3(
            new_u,
            0,
            self.deg,
            self.periodic,
            &self.poles,
            self.weights.as_deref(),
            &self.flat_knots,
            None,
        )
    }

    /// `Law_BSpline::DN` (`Law_BSpline.cxx:1303-1315`).
    pub fn dn(&self, u: f64, n: i32) -> f64 {
        curve_1d_dn(
            u,
            n,
            0,
            self.deg,
            self.periodic,
            &self.poles,
            self.weights.as_deref(),
            &self.flat_knots,
            None,
        )
    }

    /// `Law_BSpline::EndPoint` (`Law_BSpline.cxx:1319-1328`).
    pub fn end_point(&self) -> f64 {
        if self.mults[self.mults.len() - 1] == self.deg + 1 {
            self.poles[self.poles.len() - 1]
        } else {
            self.value(self.last_parameter())
        }
    }

    /// `Law_BSpline::FirstUKnotIndex` (`Law_BSpline.cxx:1333-1342`).
    pub fn first_u_knot_index(&self) -> i32 {
        if self.periodic {
            1
        } else {
            first_u_knot_index(self.deg, &self.mults)
        }
    }

    /// `Law_BSpline::FirstParameter` (`Law_BSpline.cxx:1346-1349`).
    pub fn first_parameter(&self) -> f64 {
        self.flat_knots[self.deg as usize]
    }

    /// `Law_BSpline::Knot` (`Law_BSpline.cxx:1353-1357`).
    pub fn knot(&self, index: i32) -> Result<f64, String> {
        if index < 1 || index > self.knots.len() as i32 {
            return Err("Standard_OutOfRange: Law_BSpline::Knot".into());
        }
        Ok(self.knots[(index - 1) as usize])
    }

    /// `Law_BSpline::KnotDistribution` (`Law_BSpline.cxx:1361-1364`).
    pub fn knot_distribution(&self) -> BSplKnotDistribution {
        self.knot_set
    }

    /// `Law_BSpline::Knots` (`Law_BSpline.cxx:1368-1372`).
    pub fn knots(&self) -> &[f64] {
        &self.knots
    }

    /// `Law_BSpline::KnotSequence` (`Law_BSpline.cxx:1376-1381`).
    pub fn knot_sequence(&self) -> &[f64] {
        &self.flat_knots
    }

    /// `Law_BSpline::LastUKnotIndex` (`Law_BSpline.cxx:1385-1394`).
    pub fn last_u_knot_index(&self) -> i32 {
        if self.periodic {
            self.knots.len() as i32
        } else {
            last_u_knot_index(self.deg, &self.mults)
        }
    }

    /// `Law_BSpline::LastParameter` (`Law_BSpline.cxx:1398-1401`).
    pub fn last_parameter(&self) -> f64 {
        self.flat_knots[self.flat_knots.len() - self.deg as usize - 1]
    }

    /// The `LocateParameter` + `FlatIndex` pair every `Local*` method starts
    /// with (`Law_BSpline.cxx:1411-1424`).
    pub(super) fn local_index(&self, u: f64, from_k1: i32, to_k2: i32) -> (f64, i32) {
        let (index, u) = locate_parameter_flat(
            self.deg,
            &self.flat_knots,
            u,
            self.periodic,
            from_k1,
            to_k2,
        );
        (u, flat_index(self.deg, index, &self.mults, self.periodic))
    }

    /// `Law_BSpline::Multiplicity` (`Law_BSpline.cxx:1540-1544`).
    pub fn multiplicity(&self, index: i32) -> Result<i32, String> {
        if index < 1 || index > self.mults.len() as i32 {
            return Err("Standard_OutOfRange: Law_BSpline::Multiplicity".into());
        }
        Ok(self.mults[(index - 1) as usize])
    }

    /// `Law_BSpline::Multiplicities` (`Law_BSpline.cxx:1548-1552`).
    pub fn multiplicities(&self) -> &[i32] {
        &self.mults
    }

    /// `Law_BSpline::NbKnots` (`Law_BSpline.cxx:1556-1559`).
    pub fn nb_knots(&self) -> usize {
        self.knots.len()
    }

    /// `Law_BSpline::NbPoles` (`Law_BSpline.cxx:1563-1566`).
    pub fn nb_poles(&self) -> usize {
        self.poles.len()
    }

    /// `Law_BSpline::Pole` (`Law_BSpline.cxx:1570-1574`).
    pub fn pole(&self, index: i32) -> Result<f64, String> {
        if index < 1 || index > self.poles.len() as i32 {
            return Err("Standard_OutOfRange: Law_BSpline::Pole".into());
        }
        Ok(self.poles[(index - 1) as usize])
    }

    /// `Law_BSpline::Poles` (`Law_BSpline.cxx:1578-1582`).
    pub fn poles(&self) -> &[f64] {
        &self.poles
    }

    /// `Law_BSpline::StartPoint` (`Law_BSpline.cxx:1586-1595`).
    pub fn start_point(&self) -> f64 {
        if self.mults[0] == self.deg + 1 {
            self.poles[0]
        } else {
            self.value(self.first_parameter())
        }
    }

    /// `Law_BSpline::Weight` (`Law_BSpline.cxx:1599-1609`).
    pub fn weight(&self, index: i32) -> Result<f64, String> {
        if index < 1 || index > self.poles.len() as i32 {
            return Err("Standard_OutOfRange: Law_BSpline::Weight".into());
        }
        Ok(match self.weights.as_ref() {
            Some(w) => w[(index - 1) as usize],
            None => 1.0,
        })
    }

    /// `Law_BSpline::Weights` (`Law_BSpline.cxx:1613-1628`): a non-rational
    /// curve reports `1.0` for every pole.
    pub fn weights(&self) -> Vec<f64> {
        match self.weights.as_ref() {
            Some(w) => w.clone(),
            None => vec![1.0; self.poles.len()],
        }
    }

    /// `Law_BSpline::IsRational` (`Law_BSpline.cxx:1632-1635`).
    pub fn is_rational(&self) -> bool {
        self.weights.is_some()
    }

    /// `Law_BSpline::LocateU` (`Law_BSpline.cxx:1641-1699`): the two knot
    /// indices bracketing `u`, clamped to the knot array.
    pub fn locate_u(
        &self,
        u: f64,
        parametric_tolerance: f64,
        with_knot_repetition: bool,
    ) -> (i32, i32) {
        let new_u = self.periodic_normalization(u);
        let cknots: &[f64] = if with_knot_repetition {
            &self.flat_knots
        } else {
            &self.knots
        };
        let n = cknots.len() as i32;
        let tol = parametric_tolerance.abs();
        let u_first = cknots[0];
        let u_last = cknots[cknots.len() - 1];
        let (i1, i2);
        if (u - u_first).abs() <= tol {
            i1 = 1;
            i2 = 1;
        } else if (u - u_last).abs() <= tol {
            i1 = n;
            i2 = n;
        } else if new_u < u_first - tol {
            i2 = 1;
            i1 = 0;
        } else if new_u > u_last + tol {
            i1 = n;
            i2 = i1 + 1;
        } else {
            // `BSplCLib::Hunt` (1-based result), then the OCCT clamps
            // (`Law_BSpline.cxx:1683-1691`).
            let mut first = hunt_occt(cknots, new_u);
            first = first.max(1).min(n);
            while first + 1 <= n && (cknots[first as usize] - new_u).abs() <= tol {
                first += 1;
            }
            i1 = first;
            if (cknots[(first - 1) as usize] - new_u).abs() <= tol {
                i2 = first;
            } else {
                i2 = first + 1;
            }
        }
        (i1, i2)
    }
}