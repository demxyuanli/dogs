//! Pole and weight accessors of `Geom_BSplineSurface` used by
//! `GeomFill_GordonBuilder`'s boolean sum (`GeomFill_GordonBuilder.cxx:461-502`).
//!
//! Sources (OCCT 8.0.0):
//! - `Pole(UIndex, VIndex)` (`Geom_BSplineSurface_1.cxx:573-579`),
//! - `SetPole(UIndex, VIndex, P)` (`Geom_BSplineSurface_1.cxx:1698-1703`),
//! - `Weight(UIndex, VIndex)` (`Geom_BSplineSurface_1.cxx:1265-1275`),
//! - `SetWeight(UIndex, VIndex, Weight)` (`Geom_BSplineSurface.cxx:1300-1323`).
//!
//! The port indexes from 0, so the OCCT 1-based `UIndex` / `VIndex` become the
//! `u` / `v` arguments directly.
//!
//! `SetWeight` note: OCCT keeps a weight array even for a non-rational surface
//! (`BSplSLib::UnitWeights`), and after assigning the weight it re-runs
//! `Rational(myWeights, myURational, myVRational)` and replaces the array with
//! unit weights when both flags stay clear (`cxx:1314-1322`). A non-rational
//! port surface therefore stays non-rational (`weights == None`, i.e. all unit
//! weights) as long as the assignment does not make it rational.

use super::*;

impl GeomBSplineSurface {
    /// `Geom_BSplineSurface::Pole(UIndex, VIndex)`
    /// (`Geom_BSplineSurface_1.cxx:573-579`).
    pub fn pole(&self, u: usize, v: usize) -> GpPnt {
        self.poles[u][v]
    }

    /// `Geom_BSplineSurface::SetPole(UIndex, VIndex, P)`
    /// (`Geom_BSplineSurface_1.cxx:1698-1703`): the assignment clears
    /// `myMaxDerivInvOk`.
    pub fn set_pole(&mut self, u: usize, v: usize, p: GpPnt) {
        self.poles[u][v] = p;
        self.res_latch.invalidate();
    }

    /// `Geom_BSplineSurface::Weight(UIndex, VIndex)`: the stored weight, unit
    /// when the surface keeps `BSplSLib::UnitWeights` in `myWeights`.
    pub fn weight(&self, u: usize, v: usize) -> f64 {
        match &self.weights {
            Some(w) => w[u][v],
            None => 1.0,
        }
    }

    /// `Geom_BSplineSurface::SetWeight(UIndex, VIndex, Weight)`
    /// (`Geom_BSplineSurface.cxx:1300-1323`).
    ///
    /// `Weight <= gp::Resolution()` is `Standard_ConstructionError`
    /// (`cxx:1302-1305`).
    pub fn set_weight(&mut self, u: usize, v: usize, weight: f64) -> Result<(), String> {
        if weight <= occt_core::precision::RESOLUTION {
            return Err("Geom_BSplineSurface::SetWeight: Weight too small".into());
        }
        match self.weights.as_mut() {
            Some(w) => {
                w[u][v] = weight;
            }
            None => {
                // `cxx:1312-1322`: the array starts as the unit-weight view, the
                // value is written, and the rational flags are recomputed; with
                // every weight still equal it collapses back to unit weights.
                if weight != 1.0 {
                    let n_u = self.poles.len();
                    let n_v = self.poles.first().map_or(0, |r| r.len());
                    let mut w = vec![vec![1.0f64; n_v]; n_u];
                    w[u][v] = weight;
                    self.weights = Some(w);
                }
            }
        }
        self.res_latch.invalidate();
        Ok(())
    }
}
