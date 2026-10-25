//! `GeomFill_ConstrainedFilling` — the approximation and the pole assembly.
//!
//! Source: `GeomFill_ConstrainedFilling.cxx:694-756` (`PerformApprox`),
//! `:757-1020` (`MatchKnots`), `:1024-1073` (`PerformS0`), `:1077-1329`
//! (`PerformS1`) and `:1335-1354` (`PerformSurface`) (OCCT 8.0.0).
//!
//! Translation notes:
//! - `AdvApprox_ApproxAFunction(0, 0, nbd3, ...)` with the default
//!   `AdvApprox_DichoCutting` is [`ApproxAFunction3dMulti::approx`] with one
//!   tolerance per 3D subspace.
//! - `curvpol` / `tgtepol` / `ncpol` / `ntpol` are `HArray1<gp_Pnt>` handles,
//!   `0`-based `Vec<GpPnt>` here; the 1-based OCCT indices `i`, `j` become
//!   `i - 1` / `j - 1`, including inside `PerformS1`'s `S0->Value(i1, j1)`.
//! - `S0` / `S1` / `temp` are `Array2` grids indexed `(i, j)` with `i` over the
//!   bound-0 direction; the port stores them as `Vec<Vec<GpPnt>>`.

use occt_core::adv_approx::{ApproxAFunction3dMulti, DichoCutting};
use occt_core::bspl::insert_knots::{insert_knots, prepare_insert_knots};
use occt_core::gp::{GpPnt, GpVec, GpXyz};
use occt_core::kernel::geomabs::Shape;

use super::bound_sorts::{end_term, inqadd, mklin};
use super::ConstrainedFilling;
use crate::bspline_surface::GeomBSplineSurface;
use crate::law::services::{law_mix_bnd_poles, law_mix_tgt};

impl ConstrainedFilling {
    /// `GeomFill_ConstrainedFilling::PerformApprox` (`cxx:694-756`).
    pub(super) fn perform_approx(&mut self) -> Result<(), String> {
        // `cxx:699-713`: one tolerance per 3D subspace, in the same order the
        // evaluator writes them.
        let mut tols: Vec<f64> = Vec::with_capacity(self.nbd3.max(0) as usize);
        if self.nbd3 != 0 {
            for ii in 0..2 {
                let b = self.patch()?.bound(self.ibound[ii]).clone();
                if self.ctr[ii] != 0 {
                    tols.push(b.tol3d());
                }
                if self.ctr[ii] == 2 {
                    tols.push(0.5 * self.mig[self.ibound[ii]] * b.tolang());
                }
            }
        }
        // `cxx:714-715`.
        let (f, l) = self.patch()?.bound(self.ibound[0]).bounds();

        // `cxx:717-719`: `GeomAbs_C1`, default `AdvApprox_DichoCutting`.
        let eval = |_first: f64, _last: f64, w: f64, ord: i32, out: &mut [f64]| {
            self.eval(w, ord, out)
        };
        let app = ApproxAFunction3dMulti::approx(
            f,
            l,
            Shape::C1,
            self.degmax,
            self.segmax,
            &tols,
            &DichoCutting,
            &eval,
        )
        .map_err(|_| "Standard_Failure: GeomFill_ConstrainedFilling::PerformApprox".to_string())?;

        // `cxx:721-751`.
        if app.done || app.has_result {
            let imk = self.ibound[0].min(self.ibound[1]);
            let nbpol = app.poles.first().map_or(0, |p| p.len());
            self.degree[imk] = app.degree;
            self.mults[imk] = app.mults.clone();
            self.knots[imk] = app.knots.clone();
            let mut i3d = 0usize;
            for ii in 0..2 {
                let ib = self.ibound[ii];
                let cp = if self.ctr[ii] != 0 {
                    let poles = app.poles[i3d].clone();
                    i3d += 1;
                    poles
                } else {
                    // `cxx:737-745`: a degenerate bound gets the constant point.
                    let ppp = self.patch()?.bound(ib).value(0.5 * (f + l));
                    vec![ppp; nbpol]
                };
                self.curvpol[ib] = Some(cp);
                if self.ctr[ii] == 2 {
                    self.tgtepol[ib] = Some(app.poles[i3d].clone());
                    i3d += 1;
                }
            }
        }
        Ok(())
    }

    /// `GeomFill_ConstrainedFilling::MatchKnots` (`cxx:757-1020`).
    pub(super) fn match_knots(&mut self) -> Result<(), String> {
        // `cxx:760-778`.
        let mut ind = [0usize; 4];
        self.nm = [self.mults[0].clone(), self.mults[1].clone()];
        self.nk = [self.knots[0].clone(), self.knots[1].clone()];
        ind[0] = self.nk[1].len();
        ind[2] = 1;
        ind[1] = 1;
        ind[3] = self.nk[0].len();
        self.ncpol = self.curvpol.clone();
        self.ntpol = self.tgtepol.clone();

        // `cxx:779-779`: the tolerance of the v-direction merge.
        let mut kadd = [0.0f64; 2];
        let mut madd = [0i32; 2];
        let mut tolk = 1.0 / (10.0f64).max(2.0 * self.knots[1].len() as f64);
        let nbadd = inqadd(
            self.dom[0],
            self.dom[2],
            &mut kadd,
            &mut madd,
            self.degree[1],
            tolk,
        );
        if nbadd != 0 {
            // `cxx:783-865`.
            let addk = &kadd[..nbadd];
            let addm = &madd[..nbadd];
            let prepared = prepare_insert_knots(
                self.degree[1],
                false,
                &self.knots[1],
                &self.mults[1],
                addk,
                Some(addm),
                tolk,
                false,
            );
            if let Some((_nbnp, nbnk)) = prepared {
                let deg = self.degree[1];
                let cur1 = self.curvpol[1].clone().ok_or(APPROX_MISSING)?;
                let res1 = insert_knots(
                    deg,
                    false,
                    &cur1,
                    None,
                    &self.knots[1],
                    &self.mults[1],
                    addk,
                    Some(addm),
                    tolk,
                    false,
                )
                .ok_or(APPROX_MISSING)?;
                self.nk[1] = res1.knots;
                self.nm[1] = res1.mults;
                self.ncpol[1] = Some(res1.poles);

                let cur3 = self.curvpol[3].clone().ok_or(APPROX_MISSING)?;
                let res3 = insert_knots(
                    deg,
                    false,
                    &cur3,
                    None,
                    &self.knots[1],
                    &self.mults[1],
                    addk,
                    Some(addm),
                    tolk,
                    false,
                )
                .ok_or(APPROX_MISSING)?;
                self.ncpol[3] = Some(res3.poles);

                // `cxx:830-865`.
                if let Some(t1) = self.tgtepol[1].clone() {
                    let r = insert_knots(
                        deg,
                        false,
                        &t1,
                        None,
                        &self.knots[1],
                        &self.mults[1],
                        addk,
                        Some(addm),
                        tolk,
                        false,
                    )
                    .ok_or(APPROX_MISSING)?;
                    self.ntpol[1] = Some(r.poles);
                }
                if let Some(t3) = self.tgtepol[3].clone() {
                    let r = insert_knots(
                        deg,
                        false,
                        &t3,
                        None,
                        &self.knots[1],
                        &self.mults[1],
                        addk,
                        Some(addm),
                        tolk,
                        false,
                    )
                    .ok_or(APPROX_MISSING)?;
                    self.ntpol[3] = Some(r.poles);
                }

                // `cxx:866-888`. OCCT reads `nm[1]` (the multiplicity array)
                // here, not `nk[1]`; the comparison is reproduced as written.
                if self.dom[0] != 1.0 {
                    for i in 2..=nbnk as usize {
                        if (self.dom[0] - self.nm[1][i - 1] as f64).abs() < tolk {
                            ind[0] = i;
                            break;
                        }
                    }
                }
                if self.dom[2] != 1.0 {
                    for i in 1..nbnk as usize {
                        if (1.0 - self.dom[2] - self.nm[1][i - 1] as f64).abs() < tolk {
                            ind[2] = i;
                            break;
                        }
                    }
                }
            }
        }

        // `cxx:890-998`: the same merge in the u-direction.
        tolk = 1.0 / (10.0f64).max(2.0 * self.knots[0].len() as f64);
        let nbadd = inqadd(
            self.dom[1],
            self.dom[3],
            &mut kadd,
            &mut madd,
            self.degree[0],
            tolk,
        );
        if nbadd != 0 {
            let addk = &kadd[..nbadd];
            let addm = &madd[..nbadd];
            let prepared = prepare_insert_knots(
                self.degree[0],
                false,
                &self.knots[0],
                &self.mults[0],
                addk,
                Some(addm),
                tolk,
                false,
            );
            if let Some((_nbnp, nbnk)) = prepared {
                let deg = self.degree[0];
                let cur0 = self.curvpol[0].clone().ok_or(APPROX_MISSING)?;
                let res0 = insert_knots(
                    deg,
                    false,
                    &cur0,
                    None,
                    &self.knots[0],
                    &self.mults[0],
                    addk,
                    Some(addm),
                    tolk,
                    false,
                )
                .ok_or(APPROX_MISSING)?;
                self.nk[0] = res0.knots;
                self.nm[0] = res0.mults;
                self.ncpol[0] = Some(res0.poles);

                let cur2 = self.curvpol[2].clone().ok_or(APPROX_MISSING)?;
                let res2 = insert_knots(
                    deg,
                    false,
                    &cur2,
                    None,
                    &self.knots[0],
                    &self.mults[0],
                    addk,
                    Some(addm),
                    tolk,
                    false,
                )
                .ok_or(APPROX_MISSING)?;
                self.ncpol[2] = Some(res2.poles);

                // `cxx:941-976`.
                if let Some(t0) = self.tgtepol[0].clone() {
                    let r = insert_knots(
                        deg,
                        false,
                        &t0,
                        None,
                        &self.knots[0],
                        &self.mults[0],
                        addk,
                        Some(addm),
                        tolk,
                        false,
                    )
                    .ok_or(APPROX_MISSING)?;
                    self.ntpol[0] = Some(r.poles);
                }
                if let Some(t2) = self.tgtepol[2].clone() {
                    let r = insert_knots(
                        deg,
                        false,
                        &t2,
                        None,
                        &self.knots[0],
                        &self.mults[0],
                        addk,
                        Some(addm),
                        tolk,
                        false,
                    )
                    .ok_or(APPROX_MISSING)?;
                    self.ntpol[2] = Some(r.poles);
                }

                // `cxx:977-997`. The `nm[0]` read mirrors OCCT.
                if self.dom[1] != 1.0 {
                    for i in 2..=nbnk as usize {
                        if (self.dom[1] - self.nm[0][i - 1] as f64).abs() < tolk {
                            ind[1] = i;
                            break;
                        }
                    }
                }
                if self.dom[3] != 1.0 {
                    for i in 1..nbnk as usize {
                        if (1.0 - self.dom[3] - self.nm[0][i - 1] as f64).abs() < tolk {
                            ind[3] = i;
                            break;
                        }
                    }
                }
            }
        }

        // `cxx:1000-1019`.
        let fu = mklin(self.patch()?.func(0).as_ref());
        self.ab[0] = law_mix_bnd_poles(self.degree[1], &self.nk[1], &self.nm[1], &fu);
        let fu = mklin(self.patch()?.func(1).as_ref());
        self.ab[1] = law_mix_bnd_poles(self.degree[0], &self.nk[0], &self.nm[0], &fu);
        for i in 0..2 {
            let l = self.ab[i].len();
            let mut v = vec![0.0f64; l];
            for j in 0..l {
                v[j] = 1.0 - self.ab[i][j];
            }
            self.ab[i + 2] = v;
        }
        self.pq[0] = law_mix_tgt(self.degree[1], &self.nk[1], &self.nm[1], true, ind[0]);
        self.pq[2] = law_mix_tgt(self.degree[1], &self.nk[1], &self.nm[1], false, ind[2]);
        self.pq[1] = law_mix_tgt(self.degree[0], &self.nk[0], &self.nm[0], false, ind[1]);
        self.pq[3] = law_mix_tgt(self.degree[0], &self.nk[0], &self.nm[0], true, ind[3]);
        Ok(())
    }

    /// `GeomFill_ConstrainedFilling::PerformS0` (`cxx:1024-1073`).
    pub(super) fn perform_s0(&mut self) -> Result<(), String> {
        let nc0 = self.ncpol[0].clone().ok_or(APPROX_MISSING)?;
        let nc1 = self.ncpol[1].clone().ok_or(APPROX_MISSING)?;
        let nc2 = self.ncpol[2].clone().ok_or(APPROX_MISSING)?;
        let nc3 = self.ncpol[3].clone().ok_or(APPROX_MISSING)?;
        let ni = nc0.len();
        let nj = nc1.len();
        let ptch = self.patch()?.clone();
        let c0 = *ptch.corner(0).xyz();
        let c1 = *ptch.corner(1).xyz();
        let c2 = *ptch.corner(2).xyz();
        let c3 = *ptch.corner(3).xyz();

        let mut s0: Vec<Vec<GpPnt>> = Vec::with_capacity(ni);
        for i in 0..ni {
            let ab1 = self.ab[1][i];
            let ab3 = self.ab[3][i];
            let b0 = nc0[i].xyz();
            let b2 = nc2[i].xyz();
            let mut row = Vec::with_capacity(nj);
            for j in 0..nj {
                let ab0 = self.ab[0][j];
                let ab2 = self.ab[2][j];
                let b1 = nc1[j].xyz();
                let b3 = nc3[j].xyz();
                // `cxx:1053-1069`.
                let mut polij = b0.multiplied(ab0);
                polij = polij.added(&b1.multiplied(ab1));
                polij = polij.added(&b2.multiplied(ab2));
                polij = polij.added(&b3.multiplied(ab3));
                polij = polij.added(&c0.multiplied(-ab3 * ab0));
                polij = polij.added(&c1.multiplied(-ab0 * ab1));
                polij = polij.added(&c2.multiplied(-ab1 * ab2));
                polij = polij.added(&c3.multiplied(-ab2 * ab3));
                row.push(GpPnt::from_xyz(&polij));
            }
            s0.push(row);
        }
        self.s0 = s0;
        Ok(())
    }

    /// `GeomFill_ConstrainedFilling::PerformS1` (`cxx:1077-1329`).
    pub(super) fn perform_s1(&mut self) -> Result<(), String> {
        let ncpol0 = self.ncpol[0].clone().ok_or(APPROX_MISSING)?;
        let ncpol1 = self.ncpol[1].clone().ok_or(APPROX_MISSING)?;
        let ni = ncpol0.len();
        let nj = ncpol1.len();
        if self.s0.is_empty() {
            return Err(APPROX_MISSING.to_string());
        }
        let deg0 = self.degree[0] as f64;
        let deg1 = self.degree[1] as f64;

        // `cxx:1088-1155`: the temporary tangent-field poles.
        let mut nt: [Option<Vec<GpXyz>>; 4] = [None, None, None, None];
        for i in 0..4 {
            let Some(ntpol) = self.ntpol[i].clone() else {
                continue;
            };
            let nbp = ntpol.len();
            let z;
            let (mut i1, mut i2, mut j1, mut j2);
            let inci;
            match i {
                0 => {
                    // `cxx:1101-1109`.
                    z = -deg1 / (self.nk[1][1] - self.nk[1][0]);
                    inci = true;
                    i1 = 1;
                    i2 = 1;
                    j1 = 1;
                    j2 = 2;
                }
                1 => {
                    // `cxx:1111-1119`.
                    let l = self.nk[0].len();
                    z = -deg0 / (self.nk[0][l - 1] - self.nk[0][l - 2]);
                    inci = false;
                    i1 = ni - 1;
                    i2 = ni;
                    j1 = 1;
                    j2 = 1;
                }
                2 => {
                    // `cxx:1120-1128`.
                    let l = self.nk[1].len();
                    z = -deg1 / (self.nk[1][l - 1] - self.nk[1][l - 2]);
                    inci = true;
                    i1 = 1;
                    i2 = 1;
                    j1 = nj - 1;
                    j2 = nj;
                }
                _ => {
                    // `cxx:1129-1137`.
                    z = -deg0 / (self.nk[0][1] - self.nk[0][0]);
                    inci = false;
                    i1 = 1;
                    i2 = 2;
                    j1 = 1;
                    j2 = 1;
                }
            }
            let mut arr = Vec::with_capacity(nbp);
            for k in 0..nbp {
                // `cxx:1140-1146`.
                let mut v = self.s0[i1 - 1][j1 - 1].xyz().multiplied(-1.0);
                v = v.added(self.s0[i2 - 1][j2 - 1].xyz());
                v = v.multiplied(z);
                v = v.added(ntpol[k].xyz());
                arr.push(v);
                if inci {
                    i1 += 1;
                    i2 += 1;
                } else {
                    j1 += 1;
                    j2 += 1;
                }
            }
            nt[i] = Some(arr);
        }

        // `cxx:1158-1223`: the corrective corner terms `v[4]`.
        let coef0 = deg0 / (self.nk[0][1] - self.nk[0][0]);
        let coef1 = deg1 / (self.nk[1][1] - self.nk[1][0]);
        if let (Some(t0), Some(t3)) = (nt[0].as_ref(), nt[3].as_ref()) {
            // `cxx:1164-1173`.
            let a = end_term(&t0[0], &t0[1], coef0);
            let b = end_term(&t3[0], &t3[1], coef1);
            self.v[0] = GpVec::from_xyz(&a.added(&b).multiplied(0.5));
        }

        let ln0 = self.nk[0].len();
        let lp0 = ncpol0.len();
        let coef0 = deg0 / (self.nk[0][ln0 - 1] - self.nk[0][ln0 - 2]);
        let coef1 = deg1 / (self.nk[1][1] - self.nk[1][0]);
        if let (Some(t0), Some(t1)) = (nt[0].as_ref(), nt[1].as_ref()) {
            // `cxx:1178-1189`.
            let a = end_term(&t0[lp0 - 2], &t0[lp0 - 1], coef0);
            let b = end_term(&t1[0], &t1[1], coef1);
            self.v[1] = GpVec::from_xyz(&a.added(&b).multiplied(0.5));
        }

        let ln0 = self.nk[0].len();
        let lp0 = ncpol0.len();
        let ln1 = self.nk[1].len();
        let lp1 = ncpol1.len();
        let coef0 = deg0 / (self.nk[0][ln0 - 1] - self.nk[0][ln0 - 2]);
        let coef1 = deg1 / (self.nk[1][ln1 - 1] - self.nk[1][ln1 - 2]);
        if let (Some(t1), Some(t2)) = (nt[1].as_ref(), nt[2].as_ref()) {
            // `cxx:1195-1207`.
            let a = end_term(&t2[lp0 - 2], &t2[lp0 - 1], coef0);
            let b = end_term(&t1[lp1 - 2], &t1[lp1 - 1], coef1);
            self.v[2] = GpVec::from_xyz(&a.added(&b).multiplied(0.5));
        }

        let ln1 = self.nk[1].len();
        let lp1 = ncpol1.len();
        let coef0 = deg0 / (self.nk[0][1] - self.nk[0][0]);
        let coef1 = deg1 / (self.nk[1][ln1 - 1] - self.nk[1][ln1 - 2]);
        if let (Some(t2), Some(t3)) = (nt[2].as_ref(), nt[3].as_ref()) {
            // `cxx:1210-1222`.
            let a = end_term(&t2[0], &t2[1], coef0);
            let b = end_term(&t3[lp1 - 2], &t3[lp1 - 1], coef1);
            self.v[3] = GpVec::from_xyz(&a.added(&b).multiplied(0.5));
        }

        // `cxx:1224-1321`: the S1 poles.
        let v0 = *self.v[0].xyz();
        let v1 = *self.v[1].xyz();
        let v2 = *self.v[2].xyz();
        let v3 = *self.v[3].xyz();
        let mut s1: Vec<Vec<GpPnt>> = Vec::with_capacity(ni);
        for i in 0..ni {
            let mut pq1 = 0.0;
            let mut pq3 = 0.0;
            if nt[1].is_some() {
                pq1 = -self.pq[1][i];
            }
            if nt[3].is_some() {
                pq3 = self.pq[3][i];
            }
            let mut t0 = GpXyz::new(0.0, 0.0, 0.0);
            let mut t2 = GpXyz::new(0.0, 0.0, 0.0);
            if let Some(a) = nt[0].as_ref() {
                t0 = a[i];
            }
            if let Some(a) = nt[2].as_ref() {
                t2 = a[i];
            }
            let mut row = Vec::with_capacity(nj);
            for j in 0..nj {
                let mut pq0 = 0.0;
                let mut pq2 = 0.0;
                if nt[0].is_some() {
                    pq0 = self.pq[0][j];
                }
                if nt[2].is_some() {
                    pq2 = -self.pq[2][j];
                }
                let mut t1 = GpXyz::new(0.0, 0.0, 0.0);
                let mut t3 = GpXyz::new(0.0, 0.0, 0.0);
                if let Some(a) = nt[1].as_ref() {
                    t1 = a[j];
                }
                if let Some(a) = nt[3].as_ref() {
                    t3 = a[j];
                }

                let mut tpolij = GpXyz::new(0.0, 0.0, 0.0);
                if nt[0].is_some() {
                    tpolij = tpolij.added(&t0.multiplied(pq0));
                }
                if nt[1].is_some() {
                    tpolij = tpolij.added(&t1.multiplied(pq1));
                }
                if nt[2].is_some() {
                    tpolij = tpolij.added(&t2.multiplied(pq2));
                }
                if nt[3].is_some() {
                    tpolij = tpolij.added(&t3.multiplied(pq3));
                }
                if nt[3].is_some() && nt[0].is_some() {
                    tpolij = tpolij.added(&v0.multiplied(-pq3 * pq0));
                }
                if nt[0].is_some() && nt[1].is_some() {
                    tpolij = tpolij.added(&v1.multiplied(-pq0 * pq1));
                }
                if nt[1].is_some() && nt[2].is_some() {
                    tpolij = tpolij.added(&v2.multiplied(-pq1 * pq2));
                }
                if nt[2].is_some() && nt[3].is_some() {
                    tpolij = tpolij.added(&v3.multiplied(-pq2 * pq3));
                }
                row.push(GpPnt::from_xyz(&tpolij));
            }
            s1.push(row);
        }
        self.s1 = s1;
        Ok(())
    }

    /// `GeomFill_ConstrainedFilling::PerformSurface` (`cxx:1335-1354`).
    pub(super) fn perform_surface(&mut self) -> Result<(), String> {
        let ni = self.s0.len();
        let nj = self.s0.first().map_or(0, |r| r.len());
        if ni == 0 || nj == 0 || self.s1.len() != ni {
            return Err(APPROX_MISSING.to_string());
        }
        let mut temp: Vec<Vec<GpPnt>> = Vec::with_capacity(ni);
        for i in 0..ni {
            let mut row = Vec::with_capacity(nj);
            for j in 0..nj {
                // `cxx:1341-1343`.
                let xyz = self.s0[i][j].xyz().added(self.s1[i][j].xyz());
                row.push(GpPnt::from_xyz(&xyz));
            }
            temp.push(row);
        }
        // `cxx:1345-1352`.
        let surf = GeomBSplineSurface::from_poles_knots_mults(
            temp,
            self.nk[0].clone(),
            self.nk[1].clone(),
            self.nm[0].clone(),
            self.nm[1].clone(),
            self.degree[0] as usize,
            self.degree[1] as usize,
        )?;
        self.surf = Some(surf);
        Ok(())
    }
}

/// Shared message for the poles `PerformApprox` did not produce.
const APPROX_MISSING: &str = "Standard_Failure: GeomFill_ConstrainedFilling PerformApprox";
