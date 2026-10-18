//! `Geom_BSplineSurface::segment` / `CheckAndSegment`
//! (`Geom_BSplineSurface.cxx:548-907`).
//!
//! `InsertUKnots` / `InsertVKnots` (`Geom_BSplineSurface_1.cxx:1865` and the V
//! twin) run `BSplCLib::InsertKnots` on each pole line of that direction.
//! `Add` defaults to true (`Geom_BSplineSurface.hxx:405`).
//! `SetUNotPeriodic` / `SetVNotPeriodic` (`Geom_BSplineSurface_1.cxx:1238-1288`
//! and `:1292`) follow `BSplCLib::Unperiodize`.

use occt_core::bspl::insert_knots::{insert_knots, prepare_insert_knots};
use occt_core::bspl::knots::{knot_sequence_periodic, pole_index};
use occt_core::bspl::locate::locate_parameter_range;
use occt_core::bspl::unperiodize::{flat_knots_from_mults, unperiodize_knots};
use occt_core::gp::GpPnt;
use occt_core::precision::{epsilon, PCONFUSION};

use crate::bspline_surface::GeomBSplineSurface;

impl GeomBSplineSurface {
    /// `Geom_BSplineSurface::CheckAndSegment` (`Geom_BSplineSurface.cxx:880-907`).
    ///
    /// Tolerances default to `Precision::PConfusion` (`Geom_BSplineSurface.hxx:614-615`).
    pub fn check_and_segment(
        &mut self,
        u1: f64,
        u2: f64,
        v1: f64,
        v2: f64,
    ) -> Result<(), String> {
        if u2 < u1 || v2 < v1 {
            return Err("Geom_BSplineSurface::CheckAndSegment".into());
        }
        let eps_u = epsilon(u1.abs().max(u2.abs())).max(PCONFUSION);
        let eps_v = epsilon(v1.abs().max(v2.abs())).max(PCONFUSION);
        let (uk, _) = self.distinct_knots_and_mults_u();
        let (vk, _) = self.distinct_knots_and_mults_v();
        let segment_in_u = bound_differs(&uk, u1, u2, eps_u);
        let segment_in_v = bound_differs(&vk, v1, v2, eps_v);
        self.segment(u1, u2, v1, v2, eps_u, eps_v, segment_in_u, segment_in_v)
    }

    /// `Geom_BSplineSurface::Segment` (`Geom_BSplineSurface.cxx:857-876`):
    /// both directions are segmented.
    pub fn segment_uv(
        &mut self,
        u1: f64,
        u2: f64,
        v1: f64,
        v2: f64,
        u_tolerance: f64,
        v_tolerance: f64,
    ) -> Result<(), String> {
        if u2 < u1 || v2 < v1 {
            return Err("Geom_BSplineSurface::Segment".into());
        }
        let eps_u = epsilon(u1.abs().max(u2.abs())).max(u_tolerance);
        let eps_v = epsilon(v1.abs().max(v2.abs())).max(v_tolerance);
        self.segment(u1, u2, v1, v2, eps_u, eps_v, true, true)
    }

    fn segment(
        &mut self,
        u1: f64,
        u2: f64,
        v1: f64,
        v2: f64,
        eps_u: f64,
        eps_v: f64,
        segment_in_u: bool,
        segment_in_v: bool,
    ) -> Result<(), String> {
        let mut delta_u = u2 - u1;
        if self.u_periodic {
            let (uk, _) = self.distinct_knots_and_mults_u();
            let period = period_of(&uk);
            if delta_u - period > PCONFUSION {
                return Err("Geom_BSplineSurface::Segment".into());
            }
            if delta_u > period {
                delta_u = period;
            }
        }
        let mut delta_v = v2 - v1;
        if self.v_periodic {
            let (vk, _) = self.distinct_knots_and_mults_v();
            let period = period_of(&vk);
            if delta_v - period > PCONFUSION {
                return Err("Geom_BSplineSurface::Segment".into());
            }
            if delta_v > period {
                delta_v = period;
            }
        }

        let (uk, _) = self.distinct_knots_and_mults_u();
        let (_, new_u1) = locate_knots(&uk, u1, self.u_periodic)?;
        let (_, new_u2) = locate_knots(&uk, u2, self.u_periodic)?;
        if segment_in_u {
            let knots = [new_u1.min(new_u2), new_u1.max(new_u2)];
            let mults = [self.deg_u as i32, self.deg_u as i32];
            self.insert_u_knots(&knots, &mults, eps_u)?;
        }

        let (vk, _) = self.distinct_knots_and_mults_v();
        let (_, new_v1) = locate_knots(&vk, v1, self.v_periodic)?;
        let (_, new_v2) = locate_knots(&vk, v2, self.v_periodic)?;
        if segment_in_v {
            let knots = [new_v1.min(new_v2), new_v1.max(new_v2)];
            let mults = [self.deg_v as i32, self.deg_v as i32];
            self.insert_v_knots(&knots, &mults, eps_v)?;
        }

        if self.u_periodic && segment_in_u {
            let (uk, _) = self.distinct_knots_and_mults_u();
            let (index, u) = locate_knots(&uk, u1, true)?;
            let index = bump_if_on_knot(&uk, index, u, eps_u);
            self.set_u_origin(index)?;
            self.set_u_not_periodic();
        }

        let (uk, um) = self.distinct_knots_and_mults_u();
        let (index1_u, ua) = locate_knots(&uk, new_u1, self.u_periodic)?;
        let index1_u = bump_if_on_knot(&uk, index1_u, ua, eps_u);
        let (index2_u, ub) = locate_knots(&uk, new_u1 + delta_u, self.u_periodic)?;
        let mut index2_u = bump_if_on_knot(&uk, index2_u, ub, eps_u);
        if index2_u == index1_u {
            index2_u += 1;
        }
        let (nuknots, numults) =
            slice_knots(&uk, &um, index1_u, index2_u, segment_in_u, self.deg_u as i32)?;

        if self.v_periodic && segment_in_v {
            let (vk, _) = self.distinct_knots_and_mults_v();
            let (index, v) = locate_knots(&vk, v1, true)?;
            let index = bump_if_on_knot(&vk, index, v, eps_v);
            self.set_v_origin(index)?;
            self.set_v_not_periodic();
        }

        let (vk, vm) = self.distinct_knots_and_mults_v();
        let (index1_v, va) = locate_knots(&vk, new_v1, self.v_periodic)?;
        let index1_v = bump_if_on_knot(&vk, index1_v, va, eps_v);
        let (index2_v, vb) = locate_knots(&vk, new_v1 + delta_v, self.v_periodic)?;
        let mut index2_v = bump_if_on_knot(&vk, index2_v, vb, eps_v);
        if index2_v == index1_v {
            index2_v += 1;
        }
        let (nvknots, nvmults) =
            slice_knots(&vk, &vm, index1_v, index2_v, segment_in_v, self.deg_v as i32)?;

        let (p1_u, p2_u) = pole_span(
            self.deg_u as i32,
            index1_u,
            index2_u,
            self.u_periodic,
            &um,
            self.poles.len() as i32,
        )?;
        let (p1_v, p2_v) = pole_span(
            self.deg_v as i32,
            index1_v,
            index2_v,
            self.v_periodic,
            &vm,
            self.nb_poles_v() as i32,
        )?;

        let mut npoles = Vec::with_capacity((p2_u - p1_u + 1) as usize);
        let mut nweights = if self.weights.is_some() {
            Some(Vec::with_capacity((p2_u - p1_u + 1) as usize))
        } else {
            None
        };
        for i in p1_u..=p2_u {
            let row = &self.poles[(i - 1) as usize];
            let mut nrow = Vec::with_capacity((p2_v - p1_v + 1) as usize);
            for j in p1_v..=p2_v {
                nrow.push(row[(j - 1) as usize]);
            }
            npoles.push(nrow);
            if let Some(nw) = nweights.as_mut() {
                let row = &self.weights.as_ref().unwrap()[(i - 1) as usize];
                let mut nrow = Vec::with_capacity((p2_v - p1_v + 1) as usize);
                for j in p1_v..=p2_v {
                    nrow.push(row[(j - 1) as usize]);
                }
                nw.push(nrow);
            }
        }

        // `Geom_BSplineSurface::Segment` clears `myMaxDerivInvOk` and rebuilds
        // both flat knot sequences (`Geom_BSplineSurface.cxx:850`, `:1148-1198`).
        self.res_latch.invalidate();
        self.poles = npoles;
        self.weights = nweights;
        self.knots_u = flat_of(&nuknots, &numults, self.deg_u as i32, self.u_periodic);
        self.knots_v = flat_of(&nvknots, &nvmults, self.deg_v as i32, self.v_periodic);
        Ok(())
    }

    /// `Geom_BSplineSurface::InsertUKnots` (`Geom_BSplineSurface_1.cxx:1865`),
    /// `Add = true`.
    fn insert_u_knots(&mut self, add_knots: &[f64], add_mults: &[i32], tol: f64) -> Result<(), String> {
        let (knots, mults) = self.distinct_knots_and_mults_u();
        let (nb_poles, _nb_knots) = prepare_insert_knots(
            self.deg_u as i32,
            self.u_periodic,
            &knots,
            &mults,
            add_knots,
            Some(add_mults),
            tol,
            true,
        )
        .ok_or("Geom_BSplineSurface::InsertUKnots")?;
        if nb_poles as usize == self.poles.len() {
            return Ok(());
        }
        let nb_v = self.nb_poles_v();
        let mut columns: Vec<Vec<GpPnt>> = Vec::with_capacity(nb_v);
        let mut wcolumns: Vec<Vec<f64>> = Vec::new();
        let rational = self.weights.is_some();
        let mut new_knots = Vec::new();
        let mut new_mults = Vec::new();
        for v in 0..nb_v {
            let poles: Vec<GpPnt> = self.poles.iter().map(|row| row[v]).collect();
            let weights = if rational {
                Some(
                    self.weights
                        .as_ref()
                        .unwrap()
                        .iter()
                        .map(|row| row[v])
                        .collect::<Vec<_>>(),
                )
            } else {
                None
            };
            let out = insert_knots(
                self.deg_u as i32,
                self.u_periodic,
                &poles,
                weights.as_deref(),
                &knots,
                &mults,
                add_knots,
                Some(add_mults),
                tol,
                true,
            )
            .ok_or("Geom_BSplineSurface::InsertUKnots")?;
            new_knots = out.knots;
            new_mults = out.mults;
            if rational {
                wcolumns.push(out.weights.unwrap_or_else(|| vec![1.0; out.poles.len()]));
            }
            columns.push(out.poles);
        }
        let nb_u = columns.first().map_or(0, |c| c.len());
        self.poles = (0..nb_u)
            .map(|i| (0..nb_v).map(|v| columns[v][i]).collect())
            .collect();
        if rational {
            self.weights = Some(
                (0..nb_u)
                    .map(|i| (0..nb_v).map(|v| wcolumns[v][i]).collect())
                    .collect(),
            );
        }
        // `InsertUKnots` ends in `updateUKnots()` (`Geom_BSplineSurface.cxx:1148-1152`).
        self.res_latch.invalidate();
        self.knots_u = flat_of(&new_knots, &new_mults, self.deg_u as i32, self.u_periodic);
        Ok(())
    }

    /// `Geom_BSplineSurface::InsertVKnots`, `Add = true`.
    fn insert_v_knots(&mut self, add_knots: &[f64], add_mults: &[i32], tol: f64) -> Result<(), String> {
        let (knots, mults) = self.distinct_knots_and_mults_v();
        let (nb_poles, _nb_knots) = prepare_insert_knots(
            self.deg_v as i32,
            self.v_periodic,
            &knots,
            &mults,
            add_knots,
            Some(add_mults),
            tol,
            true,
        )
        .ok_or("Geom_BSplineSurface::InsertVKnots")?;
        if nb_poles as usize == self.nb_poles_v() {
            return Ok(());
        }
        let rational = self.weights.is_some();
        let mut new_poles = Vec::with_capacity(self.poles.len());
        let mut new_weights = Vec::new();
        let mut new_knots = Vec::new();
        let mut new_mults = Vec::new();
        for (row_i, row) in self.poles.iter().enumerate() {
            let weights = if rational {
                Some(self.weights.as_ref().unwrap()[row_i].clone())
            } else {
                None
            };
            let out = insert_knots(
                self.deg_v as i32,
                self.v_periodic,
                row,
                weights.as_deref(),
                &knots,
                &mults,
                add_knots,
                Some(add_mults),
                tol,
                true,
            )
            .ok_or("Geom_BSplineSurface::InsertVKnots")?;
            new_knots = out.knots;
            new_mults = out.mults;
            if rational {
                new_weights.push(out.weights.unwrap_or_else(|| vec![1.0; out.poles.len()]));
            }
            new_poles.push(out.poles);
        }
        self.poles = new_poles;
        if rational {
            self.weights = Some(new_weights);
        }
        // `InsertVKnots` ends in `updateVKnots()` (`Geom_BSplineSurface.cxx:1196-1200`).
        self.res_latch.invalidate();
        self.knots_v = flat_of(&new_knots, &new_mults, self.deg_v as i32, self.v_periodic);
        Ok(())
    }

    /// `Geom_BSplineSurface::SetUNotPeriodic` (`Geom_BSplineSurface_1.cxx:1238-1288`).
    fn set_u_not_periodic(&mut self) {
        if !self.u_periodic {
            return;
        }
        let (knots, mults) = self.distinct_knots_and_mults_u();
        let degree = self.deg_u as i32;
        let (new_knots, new_mults, _index) = unperiodize_knots(degree, &knots, &mults);
        let n_new = (new_mults.iter().sum::<i32>() - degree - 1).max(0) as usize;
        let n_old = self.poles.len();
        if n_old == 0 || n_new == 0 {
            return;
        }
        let old = std::mem::take(&mut self.poles);
        self.poles = (0..n_new).map(|k| old[k % n_old].clone()).collect();
        if let Some(w) = self.weights.as_ref() {
            let old = w.clone();
            self.weights = Some((0..n_new).map(|k| old[k % n_old].clone()).collect());
        }
        self.knots_u = flat_knots_from_mults(&new_knots, &new_mults);
        self.u_periodic = false;
        self.res_latch.invalidate();
    }

    /// `Geom_BSplineSurface::SetVNotPeriodic` (`Geom_BSplineSurface_1.cxx:1292`).
    fn set_v_not_periodic(&mut self) {
        if !self.v_periodic {
            return;
        }
        let (knots, mults) = self.distinct_knots_and_mults_v();
        let degree = self.deg_v as i32;
        let (new_knots, new_mults, _index) = unperiodize_knots(degree, &knots, &mults);
        let n_new = (new_mults.iter().sum::<i32>() - degree - 1).max(0) as usize;
        let n_old = self.nb_poles_v();
        if n_old == 0 || n_new == 0 {
            return;
        }
        for row in self.poles.iter_mut() {
            let old = std::mem::take(row);
            *row = (0..n_new).map(|k| old[k % n_old]).collect();
        }
        if let Some(w) = self.weights.as_mut() {
            for row in w.iter_mut() {
                let old = std::mem::take(row);
                *row = (0..n_new).map(|k| old[k % n_old]).collect();
            }
        }
        self.knots_v = flat_knots_from_mults(&new_knots, &new_mults);
        self.v_periodic = false;
        self.res_latch.invalidate();
    }
}

fn bound_differs(knots: &[f64], first: f64, last: f64, eps: f64) -> bool {
    match (knots.first(), knots.last()) {
        (Some(a), Some(b)) => (first - a).abs() > eps || (last - b).abs() > eps,
        _ => true,
    }
}

fn period_of(knots: &[f64]) -> f64 {
    match (knots.first(), knots.last()) {
        (Some(a), Some(b)) => b - a,
        _ => 0.0,
    }
}

fn locate_knots(knots: &[f64], u: f64, periodic: bool) -> Result<(i32, f64), String> {
    if knots.is_empty() {
        return Err("Geom_BSplineSurface::Segment".into());
    }
    Ok(locate_parameter_range(
        knots,
        u,
        periodic,
        1,
        knots.len() as i32,
        knots[0],
        knots[knots.len() - 1],
    ))
}

/// `myUKnots.Value(index + 1)` is the 0-based slot `index`
/// (`Geom_BSplineSurface.cxx:666-668`).
fn bump_if_on_knot(knots: &[f64], index: i32, u: f64, eps: f64) -> i32 {
    if (index as usize) < knots.len() && (knots[index as usize] - u).abs() <= eps {
        index + 1
    } else {
        index
    }
}

fn slice_knots(
    knots: &[f64],
    mults: &[i32],
    index1: i32,
    index2: i32,
    raise_ends: bool,
    degree: i32,
) -> Result<(Vec<f64>, Vec<i32>), String> {
    if index1 < 1 || index2 < index1 || index2 as usize > knots.len() || index2 as usize > mults.len()
    {
        return Err("Geom_BSplineSurface::Segment".into());
    }
    let mut nk = Vec::with_capacity((index2 - index1 + 1) as usize);
    let mut nm = Vec::with_capacity((index2 - index1 + 1) as usize);
    for i in index1..=index2 {
        nk.push(knots[(i - 1) as usize]);
        nm.push(mults[(i - 1) as usize]);
    }
    if raise_ends {
        if let Some(first) = nm.first_mut() {
            *first = degree + 1;
        }
        if let Some(last) = nm.last_mut() {
            *last = degree + 1;
        }
    }
    Ok((nk, nm))
}

fn pole_span(
    degree: i32,
    index1: i32,
    index2: i32,
    periodic: bool,
    mults: &[i32],
    nb_poles: i32,
) -> Result<(i32, i32), String> {
    let mut p1 = pole_index(degree, index1, periodic, mults) + 1;
    let mut p2 = pole_index(degree, index2, periodic, mults) + 1;
    if p2 > nb_poles {
        p2 = nb_poles;
    }
    if p1 < 1 || p2 < p1 {
        return Err("Geom_BSplineSurface::Segment".into());
    }
    Ok((p1, p2))
}

fn flat_of(knots: &[f64], mults: &[i32], degree: i32, periodic: bool) -> Vec<f64> {
    if periodic {
        knot_sequence_periodic(knots, mults, degree)
    } else {
        occt_core::bspl::banded_interp::knot_sequence(knots, mults, degree)
    }
}
