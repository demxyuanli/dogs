//! Periodic-direction operations of `Geom_BSplineSurface`: `UReverse` /
//! `VReverse`, `SetUPeriodic` / `SetVPeriodic`, `LocateU` / `LocateV` and
//! `SetUOrigin` / `SetVOrigin`.

use super::*;
use occt_core::bspl::{banded_interp, curve_tools, knots};
use occt_core::gp::GpPnt;

impl GeomBSplineSurface {
    /// `Geom_BSplineSurface::UReverse` (`Geom_BSplineSurface_1.cxx:1572-1592`).
    pub fn u_reverse(&mut self) {
        if self.knots_u.len() < 2 || self.poles.is_empty() {
            return;
        }
        self.res_latch.invalidate();
        let (mut uknots, mut umults) = self.distinct_knots_and_mults_u();
        curve_tools::reverse_distinct_knots(&mut uknots);
        umults.reverse();
        if self.u_periodic {
            // `BSplSLib::Reverse` with `Last = Flat.Upper() - Deg - 1` (`cxx:1578-1586`).
            curve_tools::reverse_periodic_span(&mut self.poles, self.knots_u.len(), self.deg_u);
            if let Some(w) = self.weights.as_mut() {
                curve_tools::reverse_periodic_span(w, self.knots_u.len(), self.deg_u);
            }
            self.knots_u = knots::knot_sequence_periodic(&uknots, &umults, self.deg_u as i32);
        } else {
            self.poles.reverse();
            if let Some(w) = self.weights.as_mut() {
                w.reverse();
            }
            self.knots_u = banded_interp::knot_sequence(&uknots, &umults, self.deg_u as i32);
        }
        self.update_u_knots();
    }

    /// `Geom_BSplineSurface::VReverse` (`Geom_BSplineSurface_1.cxx:1603-1623`).
    pub fn v_reverse(&mut self) {
        if self.knots_v.len() < 2 || self.poles.is_empty() {
            return;
        }
        self.res_latch.invalidate();
        let (mut vknots, mut vmults) = self.distinct_knots_and_mults_v();
        curve_tools::reverse_distinct_knots(&mut vknots);
        vmults.reverse();
        if self.v_periodic {
            // `BSplSLib::Reverse` with `Last = Flat.Upper() - Deg - 1` (`cxx:1609-1616`).
            let flat_len = self.knots_v.len();
            let deg = self.deg_v;
            for row in self.poles.iter_mut() {
                curve_tools::reverse_periodic_span(row, flat_len, deg);
            }
            if let Some(w) = self.weights.as_mut() {
                for row in w.iter_mut() {
                    curve_tools::reverse_periodic_span(row, flat_len, deg);
                }
            }
            self.knots_v = knots::knot_sequence_periodic(&vknots, &vmults, self.deg_v as i32);
        } else {
            for row in self.poles.iter_mut() {
                row.reverse();
            }
            if let Some(w) = self.weights.as_mut() {
                for row in w.iter_mut() {
                    row.reverse();
                }
            }
            self.knots_v = banded_interp::knot_sequence(&vknots, &vmults, self.deg_v as i32);
        }
        self.update_v_knots();
    }

    /// `Geom_BSplineSurface::SetUPeriodic` (`Geom_BSplineSurface_1.cxx:940-981`).
    pub fn set_u_periodic(&mut self) {
        let (knots, mults) = Self::unique_knots_mults(&self.knots_u);
        self.set_periodic_common(true, knots, mults);
    }

    /// `Geom_BSplineSurface::SetVPeriodic` (`Geom_BSplineSurface_1.cxx:983-1022`).
    pub fn set_v_periodic(&mut self) {
        let (knots, mults) = Self::unique_knots_mults(&self.knots_v);
        self.set_periodic_common(false, knots, mults);
    }

    /// Shared body of `SetUPeriodic` / `SetVPeriodic`: keep the
    /// `FirstUKnotIndex()..LastUKnotIndex()` window, clamp the two end
    /// multiplicities to `degree`, resize the pole/weight array to
    /// `BSplCLib::NbPoles(degree, true, mults)` (`ResizeWithTrim` keeps the
    /// leading entries) and rebuild the flat periodic `BSplCLib::KnotSequence`.
    fn set_periodic_common(&mut self, is_u: bool, knots_in: Vec<f64>, mults_in: Vec<i32>) {
        if knots_in.is_empty() || mults_in.is_empty() {
            return;
        }
        self.res_latch.invalidate();
        let (degree, already) = if is_u {
            (self.deg_u, self.u_periodic)
        } else {
            (self.deg_v, self.v_periodic)
        };
        let (first, last) = if already {
            (1usize, knots_in.len())
        } else {
            (
                occt_core::bspl::locate::first_u_knot_index(degree as i32, &mults_in).max(1) as usize,
                occt_core::bspl::locate::last_u_knot_index(degree as i32, &mults_in).max(1) as usize,
            )
        };
        let first = first.min(knots_in.len());
        let last = last.min(knots_in.len()).max(first);
        let knots = knots_in[first - 1..last].to_vec();
        let mut mults = mults_in[first - 1..last].to_vec();
        let last_idx = mults.len() - 1;
        let m = (degree as i32).min(mults[0].max(mults[last_idx]));
        mults[0] = m;
        mults[last_idx] = m;
        let nbp = occt_core::bspl::knots::nb_poles(degree as i32, true, &mults).max(0) as usize;
        if nbp == 0 {
            return;
        }
        if is_u {
            let nu = self.poles.len();
            let nv = self.poles.first().map_or(0, |r| r.len());
            if nbp < nu {
                self.poles.truncate(nbp);
                if let Some(w) = self.weights.as_mut() {
                    w.truncate(nbp);
                }
            } else if nbp > nu {
                self.poles.resize(nbp, vec![GpPnt::zero(); nv]);
                if let Some(w) = self.weights.as_mut() {
                    w.resize(nbp, vec![0.0; nv]);
                }
            }
            self.knots_u =
                occt_core::bspl::knots::knot_sequence_periodic(&knots, &mults, degree as i32);
            self.u_periodic = true;
            self.update_u_knots();
        } else {
            for row in self.poles.iter_mut() {
                if nbp < row.len() {
                    row.truncate(nbp);
                } else {
                    row.resize(nbp, GpPnt::zero());
                }
            }
            if let Some(w) = self.weights.as_mut() {
                for row in w.iter_mut() {
                    if nbp < row.len() {
                        row.truncate(nbp);
                    } else {
                        row.resize(nbp, 0.0);
                    }
                }
            }
            self.knots_v =
                occt_core::bspl::knots::knot_sequence_periodic(&knots, &mults, degree as i32);
            self.v_periodic = true;
            self.update_v_knots();
        }
    }

    /// `Geom_BSplineSurface::LocateU` (`Geom_BSplineSurface_1.cxx:1464-1514`),
    /// unique knots (`WithKnotRepetition = false`).
    ///
    /// `false` selects `myUKnots` (`:1471`) - the distinct knot **window**
    /// `FirstUKnotIndex()..LastUKnotIndex()` (`:1408-1445` returns `1..Length`
    /// for a periodic surface), not the flat sequence, whose extension knots
    /// would shift the returned indices by the wrapped tail and widen the
    /// `PeriodicNormalization` period.
    pub fn locate_u(&self, u: f64, parametric_tolerance: f64) -> (i32, i32) {
        Self::locate_param(u, parametric_tolerance, &self.uknots, self.u_periodic)
    }

    /// `Geom_BSplineSurface::LocateV` (`Geom_BSplineSurface_1.cxx:1518-1560`);
    /// see [`Self::locate_u`] for the `myVKnots` window.
    pub fn locate_v(&self, v: f64, parametric_tolerance: f64) -> (i32, i32) {
        Self::locate_param(v, parametric_tolerance, &self.vknots, self.v_periodic)
    }

    /// Knot window + pole offset rotation shared by `SetUOrigin` / `SetVOrigin`
    /// (`Geom_BSplineSurface_1.cxx:1050-1074` / `:1157-1177`): the distinct knot
    /// window is re-cut at `index` - `knots(index..=last)` followed by
    /// `knots(first+1..=index)` shifted by one period - and the returned offset
    /// is the 1-based pole/weight index the rotation starts from,
    /// `1 + sum(mults(first+1..=index))`.
    fn rotate_origin(knots: &[f64], mults: &[i32], index: i32) -> (Vec<f64>, Vec<i32>, usize) {
        // `FirstUKnotIndex()` / `LastUKnotIndex()` are `1` / `Knots().Length()`
        // for a periodic direction (`Geom_BSplineSurface_1.cxx:1408-1445`).
        let first = 1i32;
        let last = knots.len() as i32;
        let period = knots[(last - 1) as usize] - knots[(first - 1) as usize];
        let mut newknots = Vec::with_capacity(knots.len());
        let mut newmults = Vec::with_capacity(mults.len());
        for i in index..=last {
            newknots.push(knots[(i - 1) as usize]);
            newmults.push(mults[(i - 1) as usize]);
        }
        for i in (first + 1)..=index {
            newknots.push(knots[(i - 1) as usize] + period);
            newmults.push(mults[(i - 1) as usize]);
        }
        let mut offset = 1usize;
        for i in (first + 1)..=index {
            offset += mults[(i - 1) as usize] as usize;
        }
        (newknots, newmults, offset)
    }

    /// `Geom_BSplineSurface::SetUOrigin` (`Geom_BSplineSurface_1.cxx:1026-1127`):
    /// rotate the periodic U representation so its knot window starts at the
    /// 1-based window index `index` (as [`Self::locate_u`] reports it), carrying
    /// the poles and weights along. `updateUKnots()` is the flat periodic
    /// sequence rebuilt at the end.
    pub fn set_u_origin(&mut self, index: i32) -> Result<(), String> {
        if !self.u_periodic {
            return Err("Geom_BSplineSurface::SetUOrigin: surface is not U periodic".to_string());
        }
        let (knots, mults) = self.distinct_knots_and_mults_u();
        if knots.is_empty() || index < 1 || index > knots.len() as i32 {
            return Err("Geom_BSplineSurface::SetUOrigin: Index out of range".to_string());
        }
        // Rotating the window rebuilds the pole grid and ends in `updateUKnots()`,
        // which clears `myMaxDerivInvOk` (`Geom_BSplineSurface.cxx:1148-1152`).
        self.res_latch.invalidate();
        let (newknots, newmults, offset) = Self::rotate_origin(&knots, &mults, index);
        let nb_u = self.poles.len();
        if offset > nb_u + 1 {
            return Err("Geom_BSplineSurface::SetUOrigin: pole offset out of range".to_string());
        }
        // `for (i = index; i <= last; i++)` then `for (i = first; i < index; i++)`
        // over the pole rows (`:1073-1121`); `first` is the pole array's lower
        // index (1).
        let order: Vec<usize> = (offset..=nb_u).chain(1..offset).map(|i| i - 1).collect();
        let old = std::mem::take(&mut self.poles);
        self.poles = order.iter().map(|&i| old[i].clone()).collect();
        if let Some(w) = self.weights.as_mut() {
            let old = std::mem::take(w);
            *w = order.iter().map(|&i| old[i].clone()).collect();
        }
        // `myUKnots`/`myUMults` assigned, then `updateUKnots()` rebuilds the flat
        // periodic sequence (`Geom_BSplineSurface.cxx:1148-1192`).
        let flat = occt_core::bspl::knots::knot_sequence_periodic(
            &newknots,
            &newmults,
            self.deg_u as i32,
        );
        // `BSplCLib::NbPoles(degree, true, mults)` (`BSplCLib.cxx:392-451`): for a
        // periodic direction `Mults(first) == Mults(last)` and the pole count is
        // `sum(mults) - Mults(first)`.
        let sum: i32 = newmults.iter().sum();
        if sum - newmults[0] != nb_u as i32 {
            return Err("Geom_BSplineSurface::SetUOrigin: pole count mismatch".to_string());
        }
        self.knots_u = flat;
        self.update_u_knots();
        Ok(())
    }

    /// `Geom_BSplineSurface::SetVOrigin` (`Geom_BSplineSurface_1.cxx:1132-1233`);
    /// the V pole rotation moves whole columns. See [`Self::set_u_origin`].
    pub fn set_v_origin(&mut self, index: i32) -> Result<(), String> {
        if !self.v_periodic {
            return Err("Geom_BSplineSurface::SetVOrigin: surface is not V periodic".to_string());
        }
        let (knots, mults) = self.distinct_knots_and_mults_v();
        if knots.is_empty() || index < 1 || index > knots.len() as i32 {
            return Err("Geom_BSplineSurface::SetVOrigin: Index out of range".to_string());
        }
        // See [`Self::set_u_origin`]: the rotation ends in `updateVKnots()`.
        self.res_latch.invalidate();
        let (newknots, newmults, offset) = Self::rotate_origin(&knots, &mults, index);
        let nb_v = self.nb_poles_v();
        if offset > nb_v + 1 {
            return Err("Geom_BSplineSurface::SetVOrigin: pole offset out of range".to_string());
        }
        // `for (j = index; j <= last; j++)` then `for (j = first; j < index; j++)`
        // over the pole columns (`:1185-1228`).
        let order: Vec<usize> = (offset..=nb_v).chain(1..offset).map(|j| j - 1).collect();
        for row in self.poles.iter_mut() {
            let old = std::mem::take(row);
            *row = order.iter().map(|&j| old[j]).collect();
        }
        if let Some(w) = self.weights.as_mut() {
            for row in w.iter_mut() {
                let old = std::mem::take(row);
                *row = order.iter().map(|&j| old[j]).collect();
            }
        }
        let flat = occt_core::bspl::knots::knot_sequence_periodic(
            &newknots,
            &newmults,
            self.deg_v as i32,
        );
        let sum: i32 = newmults.iter().sum();
        if sum - newmults[0] != nb_v as i32 {
            return Err("Geom_BSplineSurface::SetVOrigin: pole count mismatch".to_string());
        }
        self.knots_v = flat;
        self.update_v_knots();
        Ok(())
    }
}
