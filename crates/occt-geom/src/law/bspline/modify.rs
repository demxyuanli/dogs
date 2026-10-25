//! `Law_BSpline` editing — degree, multiplicity, knots, poles, weights,
//! `Reverse`, `Segment` and `MovePointAndTangent`.
//!
//! Source: `Law_BSpline.cxx:311-378` (`IncreaseDegree`), `:382-404`
//! (`IncreaseMultiplicity`), `:408-455` (`IncrementMultiplicity`), `:459-530`
//! (`InsertKnot`/`InsertKnots`), `:535-635` (`RemoveKnot`), `:640-662`
//! (`Reverse`), `:667-760` (`Segment`), `:765-830` (`SetKnot`/`SetKnots`),
//! `:835-885` (`SetPeriodic`), `:890-930` (`SetOrigin`), `:935-955`
//! (`SetNotPeriodic`), `:960-1060` (`SetPole`/`SetWeight`) and `:1450-1530`
//! (`MovePointAndTangent`), all in OCCT 8.0.0.
//!
//! Translation notes: see the module header of [`super`].

use occt_core::bspl::curve_tools::{reverse_distinct_knots, reverse_periodic_span};
use occt_core::bspl::increase_degree::{increase_degree_count_knots, increase_degree_dim};
use occt_core::bspl::insert_knots::{insert_knots_dim, prepare_insert_knots};
use occt_core::bspl::knots::{nb_poles, pole_index};
use occt_core::bspl::locate::locate_parameter_knots_range;
use occt_core::bspl::move_point::{move_point_and_tangent_1d, MovePointError};
use occt_core::bspl::remove_knot::remove_knot_dim;
use occt_core::bspl::unperiodize::{prepare_unperiodize, unperiodize_knots};
use occt_core::precision::{epsilon, REAL_SMALL};

use super::{check_curve_data, get_poles, set_poles, weights_are_rational, LawBSpline};

impl LawBSpline {
    /// `Law_BSpline::IncreaseDegree` (`Law_BSpline.cxx:311-378`).
    pub fn increase_degree(&mut self, degree: i32) -> Result<(), String> {
        if degree == self.deg {
            return Ok(());
        }
        if degree < self.deg || degree > Self::MAX_DEGREE {
            return Err("Standard_ConstructionError: Law_BSpline".into());
        }
        let from_k1 = self.first_u_knot_index();
        let to_k2 = self.last_u_knot_index();
        let step = degree - self.deg;
        let n_poles = (self.poles.len() as i32 + step * (to_k2 - from_k1)).max(0) as usize;
        let n_knots = increase_degree_count_knots(self.deg, degree, self.periodic, &self.mults)
            .max(0) as usize;
        let mut new_poles = vec![0.0f64; n_poles];
        let mut new_knots = vec![0.0f64; n_knots];
        let mut new_mults = vec![0i32; n_knots];
        let mut new_weights: Option<Vec<f64>> = None;
        if self.rational {
            let w = self.weights.as_deref().unwrap_or(&[]);
            let adim = set_poles(&self.poles, w);
            let mut adimn = vec![0.0f64; 2 * n_poles];
            increase_degree_dim(
                self.deg,
                degree,
                self.periodic,
                2,
                &adim,
                &self.knots,
                &self.mults,
                &mut adimn,
                &mut new_knots,
                &mut new_mults,
            );
            let (p, w) = get_poles(&adimn, n_poles);
            new_poles = p;
            new_weights = Some(w);
        } else {
            increase_degree_dim(
                self.deg,
                degree,
                self.periodic,
                1,
                &self.poles,
                &self.knots,
                &self.mults,
                &mut new_poles,
                &mut new_knots,
                &mut new_mults,
            );
        }
        self.deg = degree;
        self.poles = new_poles;
        self.weights = new_weights;
        self.knots = new_knots;
        self.mults = new_mults;
        self.update_knots();
        Ok(())
    }

    /// `Law_BSpline::IncreaseMultiplicity(Index, M)` (`Law_BSpline.cxx:382-387`).
    pub fn increase_multiplicity(&mut self, index: i32, m: i32) -> Result<(), String> {
        let k = self.knot(index)?;
        let delta = m - self.multiplicity(index)?;
        self.insert_knots(&[k], Some(&[delta]), epsilon(1.0), false)
    }

    /// `Law_BSpline::IncreaseMultiplicity(I1, I2, M)` (`Law_BSpline.cxx:391-404`).
    pub fn increase_multiplicity_range(
        &mut self,
        i1: i32,
        i2: i32,
        m: i32,
    ) -> Result<(), String> {
        let mut k = Vec::new();
        let mut dm = Vec::new();
        for i in i1..=i2 {
            k.push(self.knot(i)?);
            dm.push(m - self.multiplicity(i)?);
        }
        self.insert_knots(&k, Some(&dm), epsilon(1.0), false)
    }

    /// `Law_BSpline::IncrementMultiplicity(I1, I2, Step)` (`cxx:408-415`).
    pub fn increment_multiplicity(
        &mut self,
        i1: i32,
        i2: i32,
        step: i32,
    ) -> Result<(), String> {
        let mut k = Vec::new();
        let mut dm = Vec::new();
        for i in i1..=i2 {
            k.push(self.knot(i)?);
            dm.push(step);
        }
        self.insert_knots(&k, Some(&dm), epsilon(1.0), false)
    }

    /// `Law_BSpline::InsertKnot` (`Law_BSpline.cxx:418-428`).
    pub fn insert_knot(
        &mut self,
        u: f64,
        m: i32,
        parametric_tolerance: f64,
        add: bool,
    ) -> Result<(), String> {
        self.insert_knots(&[u], Some(&[m]), parametric_tolerance, add)
    }

    /// `Law_BSpline::InsertKnots` (`Law_BSpline.cxx:434-512`).
    pub fn insert_knots(
        &mut self,
        add_knots: &[f64],
        add_mults: Option<&[i32]>,
        eps: f64,
        add: bool,
    ) -> Result<(), String> {
        let (nbpoles, nbknots) = prepare_insert_knots(
            self.deg,
            self.periodic,
            &self.knots,
            &self.mults,
            add_knots,
            add_mults,
            eps,
            add,
        )
        .ok_or_else(|| "Standard_ConstructionError: Law_BSpline::InsertKnots".to_string())?;
        if nbpoles as usize == self.poles.len() {
            return Ok(());
        }
        let n_poles = nbpoles.max(0) as usize;
        let n_knots = nbknots.max(0) as usize;
        let mut new_poles = vec![0.0f64; n_poles];
        let mut new_knots = self.knots.clone();
        let mut new_mults = self.mults.clone();
        if n_knots != self.knots.len() {
            new_knots = vec![0.0f64; n_knots];
            new_mults = vec![0i32; n_knots];
        }
        if self.rational {
            let w = self.weights.as_deref().unwrap_or(&[]);
            let adim = set_poles(&self.poles, w);
            let mut adimn = vec![0.0f64; 2 * n_poles];
            insert_knots_dim(
                self.deg,
                self.periodic,
                2,
                &adim,
                &self.knots,
                &self.mults,
                add_knots,
                add_mults,
                &mut adimn,
                &mut new_knots,
                &mut new_mults,
                eps,
                add,
            );
            let (p, w) = get_poles(&adimn, n_poles);
            self.poles = p;
            self.weights = Some(w);
        } else {
            insert_knots_dim(
                self.deg,
                self.periodic,
                1,
                &self.poles,
                &self.knots,
                &self.mults,
                add_knots,
                add_mults,
                &mut new_poles,
                &mut new_knots,
                &mut new_mults,
                eps,
                add,
            );
            self.poles = new_poles;
        }
        self.knots = new_knots;
        self.mults = new_mults;
        self.update_knots();
        Ok(())
    }

    /// `Law_BSpline::RemoveKnot` (`Law_BSpline.cxx:516-605`).
    pub fn remove_knot(&mut self, index: i32, m: i32, tolerance: f64) -> Result<bool, String> {
        if m < 0 {
            return Ok(true);
        }
        let i1 = self.first_u_knot_index();
        let i2 = self.last_u_knot_index();
        if !self.periodic && (index <= i1 || index >= i2) {
            return Err("Standard_OutOfRange: Law_BSpline::RemoveKnot".into());
        } else if self.periodic && (index < i1 || index > i2) {
            return Err("Standard_OutOfRange: Law_BSpline::RemoveKnot".into());
        }
        let step = self.multiplicity(index)? - m;
        if step <= 0 {
            return Ok(true);
        }
        let n_poles = self.poles.len().saturating_sub(step as usize);
        let mut new_poles = vec![0.0f64; n_poles];
        let mut new_knots = self.knots.clone();
        let mut new_mults = self.mults.clone();
        if m == 0 {
            new_knots = vec![0.0f64; self.knots.len().saturating_sub(1)];
            new_mults = vec![0i32; self.knots.len().saturating_sub(1)];
        }
        if self.rational {
            let w = self.weights.as_deref().unwrap_or(&[]);
            let adim = set_poles(&self.poles, w);
            let mut adimn = vec![0.0f64; 2 * n_poles];
            if !remove_knot_dim(
                index,
                m,
                self.deg,
                self.periodic,
                2,
                &adim,
                &self.knots,
                &self.mults,
                &mut adimn,
                &mut new_knots,
                &mut new_mults,
                tolerance,
            ) {
                return Ok(false);
            }
            let (p, w) = get_poles(&adimn, n_poles);
            self.poles = p;
            self.weights = Some(w);
        } else {
            if !remove_knot_dim(
                index,
                m,
                self.deg,
                self.periodic,
                1,
                &self.poles,
                &self.knots,
                &self.mults,
                &mut new_poles,
                &mut new_knots,
                &mut new_mults,
                tolerance,
            ) {
                return Ok(false);
            }
            self.poles = new_poles;
        }
        self.knots = new_knots;
        self.mults = new_mults;
        self.update_knots();
        Ok(true)
    }

    /// `Law_BSpline::Reverse` (`Law_BSpline.cxx:612-630`).
    pub fn reverse(&mut self) {
        // `BSplCLib::Reverse(Knots)` (`BSplCLib.cxx:802-824`) and
        // `BSplCLib::Reverse(Mults)` (`cxx:828-831`).
        reverse_distinct_knots(&mut self.knots);
        self.mults.reverse();
        // `BSplCLib::Reverse(Poles, last)` (`BSplCLib_CurveComputation.pxx:275-288`)
        // with `last = periodic ? flatknots->Upper() - deg - 1 : poles->Upper()`.
        reverse_periodic_span(
            &mut self.poles,
            self.flat_knots.len(),
            self.deg.max(0) as usize,
        );
        if let Some(w) = self.weights.as_mut() {
            reverse_periodic_span(w, self.flat_knots.len(), self.deg.max(0) as usize);
        }
        self.update_knots();
    }

    /// `Law_BSpline::ReversedParameter` (`Law_BSpline.cxx:634-637`).
    pub fn reversed_parameter(&self, u: f64) -> f64 {
        self.first_parameter() + self.last_parameter() - u
    }

    /// `Law_BSpline::Segment` (`Law_BSpline.cxx:642-783`).
    pub fn segment(&mut self, u1: f64, u2: f64) -> Result<(), String> {
        if u2 < u1 {
            return Err("Standard_DomainError: Law_BSpline::Segment".into());
        }
        let eps = epsilon(u1.abs().max(u2.abs()));
        let delta = u2 - u1;

        let (_, new_u1) = locate_parameter_knots_range(
            &self.knots,
            u1,
            self.periodic,
            1,
            self.knots.len() as i32,
        );
        let (_, new_u2) = locate_parameter_knots_range(
            &self.knots,
            u2,
            self.periodic,
            1,
            self.knots.len() as i32,
        );
        let add_knots = [new_u1.min(new_u2), new_u1.max(new_u2)];
        let add_mults = [self.deg, self.deg];
        self.insert_knots(&add_knots, Some(&add_mults), eps, false)?;

        if self.periodic {
            // `cxx:679-698`: set the origin at NewU1, then unperiodize.
            let (mut index0, u) = locate_parameter_knots_range(
                &self.knots,
                u1,
                true,
                1,
                self.knots.len() as i32,
            );
            if let Some(k) = self.knots.get(index0 as usize) {
                if (*k - u).abs() < eps {
                    index0 += 1;
                }
            }
            self.set_origin(index0)?;
            self.set_not_periodic();
        }

        // `cxx:700-741`: new knots/mults between NewU1 and NewU1 + delta.
        let from_u1 = 1;
        let to_u2 = self.knots.len() as i32;
        let (index1, _) =
            locate_parameter_knots_range(&self.knots, new_u1, self.periodic, from_u1, to_u2);
        let (index2, u) = locate_parameter_knots_range(
            &self.knots,
            new_u1 + delta,
            self.periodic,
            from_u1,
            to_u2,
        );
        let mut index2 = index2;
        if let Some(k) = self.knots.get(index2 as usize) {
            if (*k - u).abs() < eps {
                index2 += 1;
            }
        }
        let nbknots = index2 - index1 + 1;
        let mut new_knots = Vec::with_capacity(nbknots.max(0) as usize);
        let mut new_mults = Vec::with_capacity(nbknots.max(0) as usize);
        for i in index1..=index2 {
            if let (Some(k), Some(m)) = (
                self.knots.get((i - 1) as usize),
                self.mults.get((i - 1) as usize),
            ) {
                new_knots.push(*k);
                new_mults.push(*m);
            }
        }
        if let Some(first) = new_mults.first_mut() {
            *first = self.deg + 1;
        }
        if let Some(last) = new_mults.last_mut() {
            *last = self.deg + 1;
        }

        // `cxx:742-772`: new poles/weights.
        let mut pindex1 = pole_index(self.deg, index1, self.periodic, &self.mults) + 1;
        let mut pindex2 = pole_index(self.deg, index2, self.periodic, &self.mults) + 1;
        pindex2 = pindex2.min(self.poles.len() as i32);
        if pindex1 < 1 {
            pindex1 = 1;
        }
        let nbpoles = (pindex2 - pindex1 + 1).max(0) as usize;
        let mut new_poles = Vec::with_capacity(nbpoles);
        let mut new_weights = Vec::with_capacity(nbpoles);
        for i in pindex1..=pindex2 {
            if let Some(p) = self.poles.get((i - 1) as usize) {
                new_poles.push(*p);
                if let Some(w) = self.weights.as_ref() {
                    if let Some(wi) = w.get((i - 1) as usize) {
                        new_weights.push(*wi);
                    }
                }
            }
        }

        self.knots = new_knots;
        self.mults = new_mults;
        self.poles = new_poles;
        if self.rational {
            self.weights = Some(new_weights);
        }
        self.update_knots();
        Ok(())
    }

    /// `Law_BSpline::SetKnot(Index, K)` (`Law_BSpline.cxx:786-819`).
    pub fn set_knot(&mut self, index: i32, k: f64) -> Result<(), String> {
        if index < 1 || index > self.knots.len() as i32 {
            return Err("Standard_OutOfRange: Law_BSpline::SetKnot".into());
        }
        let dk = epsilon(k).abs();
        let n = self.knots.len() as i32;
        if index == 1 {
            if k >= self.knots[1] - dk {
                return Err("Standard_ConstructionError: Law_BSpline::SetKnot".into());
            }
        } else if index == n {
            if k <= self.knots[(n - 2) as usize] + dk {
                return Err("Standard_ConstructionError: Law_BSpline::SetKnot".into());
            }
        } else if k <= self.knots[(index - 2) as usize] + dk
            || k >= self.knots[index as usize] - dk
        {
            return Err("Standard_ConstructionError: Law_BSpline::SetKnot".into());
        }
        if k != self.knots[(index - 1) as usize] {
            self.knots[(index - 1) as usize] = k;
            self.update_knots();
        }
        Ok(())
    }

    /// `Law_BSpline::SetKnots` (`Law_BSpline.cxx:823-828`).
    pub fn set_knots(&mut self, k: Vec<f64>) -> Result<(), String> {
        check_curve_data(&self.poles, &k, &self.mults, self.deg, self.periodic)?;
        self.knots = k;
        self.update_knots();
        Ok(())
    }

    /// `Law_BSpline::SetKnot(Index, K, M)` (`Law_BSpline.cxx:832-836`).
    pub fn set_knot_with_multiplicity(
        &mut self,
        index: i32,
        k: f64,
        m: i32,
    ) -> Result<(), String> {
        self.increase_multiplicity(index, m)?;
        self.set_knot(index, k)
    }

    /// `Law_BSpline::SetPeriodic` (`Law_BSpline.cxx:840-876`): turns a closed
    /// curve into a periodic one.
    pub fn set_periodic(&mut self) {
        let first = self.first_u_knot_index();
        let last = self.last_u_knot_index();
        let mut cknots = Vec::new();
        let mut cmults = Vec::new();
        for i in first..=last {
            if let (Some(k), Some(m)) = (
                self.knots.get((i - 1) as usize),
                self.mults.get((i - 1) as usize),
            ) {
                cknots.push(*k);
                cmults.push(*m);
            }
        }
        if cmults.is_empty() {
            return;
        }
        let m = cmults[0].max(*cmults.last().unwrap());
        cmults[0] = m;
        let last_idx = cmults.len() - 1;
        cmults[last_idx] = m;
        let nbp = nb_poles(self.deg, true, &cmults).max(1) as usize;
        let mut cpoles: Vec<f64> = self.poles.iter().copied().take(nbp).collect();
        cpoles.resize(nbp, 0.0);
        let mut cweights = self.weights.as_ref().map(|w| {
            let mut v: Vec<f64> = w.iter().copied().take(nbp).collect();
            v.resize(nbp, 0.0);
            v
        });
        if let Some(w) = cweights.as_mut() {
            // `Law_BSpline.cxx:863-870`: OCCT's `Resize` keeps the appended
            // weight slots default-constructed; the port keeps `0.0` there.
            if w.len() > self.poles.len() {
                for item in w.iter_mut().skip(self.poles.len()) {
                    *item = 0.0;
                }
            }
        }
        self.poles = cpoles;
        self.weights = cweights;
        self.knots = cknots;
        self.mults = cmults;
        self.periodic = true;
        self.update_knots();
    }

    /// `Law_BSpline::SetOrigin` (`Law_BSpline.cxx:879-965`).
    pub fn set_origin(&mut self, index: i32) -> Result<(), String> {
        if !self.periodic {
            return Err("Standard_NoSuchObject: Law_BSpline::SetOrigin".into());
        }
        let first = self.first_u_knot_index();
        let last = self.last_u_knot_index();
        if index < first || index > last {
            return Err("Standard_DomainError: Law_BSpline::SetOrigin".into());
        }
        let nbknots = self.knots.len() as i32;
        let nbpoles = self.poles.len() as i32;
        let mut new_knots = Vec::with_capacity(nbknots.max(0) as usize);
        let mut new_mults = Vec::with_capacity(nbknots.max(0) as usize);
        let period = self.knots[(last - 1) as usize] - self.knots[(first - 1) as usize];
        for i in index..=last {
            new_knots.push(self.knots[(i - 1) as usize]);
            new_mults.push(self.mults[(i - 1) as usize]);
        }
        for i in (first + 1)..=index {
            new_knots.push(self.knots[(i - 1) as usize] + period);
            new_mults.push(self.mults[(i - 1) as usize]);
        }
        let mut origin = 1;
        for i in (first + 1)..=index {
            origin += self.mults[(i - 1) as usize];
        }
        let mut new_poles = Vec::with_capacity(nbpoles.max(0) as usize);
        let mut new_weights = Vec::with_capacity(nbpoles.max(0) as usize);
        for i in origin..=nbpoles {
            new_poles.push(self.poles[(i - 1) as usize]);
            if let Some(w) = self.weights.as_ref() {
                new_weights.push(w[(i - 1) as usize]);
            }
        }
        for i in 1..origin {
            new_poles.push(self.poles[(i - 1) as usize]);
            if let Some(w) = self.weights.as_ref() {
                new_weights.push(w[(i - 1) as usize]);
            }
        }
        self.poles = new_poles;
        self.knots = new_knots;
        self.mults = new_mults;
        if self.rational {
            self.weights = Some(new_weights);
        }
        self.update_knots();
        Ok(())
    }

    /// `Law_BSpline::SetNotPeriodic` (`Law_BSpline.cxx:969-1023`).
    pub fn set_not_periodic(&mut self) {
        if !self.periodic {
            return;
        }
        let (nb_knots, nb_poles) = prepare_unperiodize(self.deg, &self.mults);
        let (new_knots, new_mults, _index) = unperiodize_knots(self.deg, &self.knots, &self.mults);
        let n_old = self.poles.len();
        if n_old == 0 || nb_poles == 0 || nb_knots == 0 {
            return;
        }
        // `BSplCLib::Unperiodize` (`BSplCLib.cxx:3062-3080`):
        // `NewPoles(k) = Poles((k - 1) % n_old + 1)`.
        let new_poles: Vec<f64> = (0..nb_poles).map(|k| self.poles[k % n_old]).collect();
        let new_weights = self.weights.as_ref().map(|w| {
            let n = w.len();
            (0..nb_poles).map(|k| w[k % n]).collect::<Vec<f64>>()
        });
        self.poles = new_poles;
        self.weights = new_weights;
        self.mults = new_mults;
        self.knots = new_knots;
        self.periodic = false;
        self.update_knots();
    }

    /// `Law_BSpline::SetPole(Index, P)` (`Law_BSpline.cxx:1026-1033`).
    pub fn set_pole(&mut self, index: i32, p: f64) -> Result<(), String> {
        if index < 1 || index > self.poles.len() as i32 {
            return Err("Standard_OutOfRange: Law_BSpline::SetPole".into());
        }
        self.poles[(index - 1) as usize] = p;
        Ok(())
    }

    /// `Law_BSpline::SetPole(Index, P, W)` (`Law_BSpline.cxx:1036-1040`).
    pub fn set_pole_with_weight(&mut self, index: i32, p: f64, w: f64) -> Result<(), String> {
        self.set_pole(index, p)?;
        self.set_weight(index, w)
    }

    /// `Law_BSpline::SetWeight` (`Law_BSpline.cxx:1044-1080`).
    pub fn set_weight(&mut self, index: i32, w: f64) -> Result<(), String> {
        if index < 1 || index > self.poles.len() as i32 {
            return Err("Standard_OutOfRange: Law_BSpline::SetWeight".into());
        }
        if w <= REAL_SMALL {
            return Err("Standard_ConstructionError: Law_BSpline::SetWeight".into());
        }
        let mut rat = self.is_rational() || (w - 1.0).abs() > REAL_SMALL;
        if rat {
            if !self.is_rational() {
                self.weights = Some(vec![1.0; self.poles.len()]);
            }
            if let Some(ws) = self.weights.as_mut() {
                ws[(index - 1) as usize] = w;
            }
            if self.is_rational() {
                rat = weights_are_rational(self.weights.as_deref().unwrap_or(&[]));
                if !rat {
                    self.weights = None;
                }
            }
            self.rational = self.weights.is_some();
        }
        Ok(())
    }

    /// `Law_BSpline::MovePointAndTangent` (`Law_BSpline.cxx:1703-1738`).
    ///
    /// Returns the OCCT `ErrorStatus` (`0` = done, `1` = not enough freedom at
    /// the ends, `2` = invalid conditions); the poles are rewritten only when
    /// it is `0`. A singular 2x2 system raises
    /// `Standard_ConstructionError` from `gp_Mat2d::Invert`.
    pub fn move_point_and_tangent(
        &mut self,
        u: f64,
        p: f64,
        tangent: f64,
        tolerance: f64,
        starting_condition: i32,
        ending_condition: i32,
    ) -> Result<i32, String> {
        let (delta, delta_derivative) = self.d1(u);
        let delta = p - delta;
        let delta_derivative = tangent - delta_derivative;
        match move_point_and_tangent_1d(
            u,
            delta,
            delta_derivative,
            tolerance,
            self.deg,
            starting_condition,
            ending_condition,
            &self.poles,
            self.weights.as_deref(),
            &self.flat_knots,
        ) {
            Ok(poles) => {
                self.poles = poles;
                Ok(0)
            }
            Err(MovePointError::Status(status)) => Ok(status),
            Err(MovePointError::SingularMatrix) => {
                Err("Standard_ConstructionError: gp_Mat2d::Invert".into())
            }
        }
    }
}