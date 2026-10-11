//! `Geom_BezierSurface` insert and remove pole columns and rows.
//!
//! Helpers match `AddPoleCol` / `AddRatPoleCol` / `AddPoleRow` /
//! `AddRatPoleRow` / `DeletePoleCol` / `DeleteRatPoleCol` /
//! `DeletePoleRow` / `DeleteRatPoleRow` in `Geom_BezierSurface.cxx:98-376`.

use super::GeomBezierSurface;
use occt_core::gp::GpPnt;

impl GeomBezierSurface {
    /// `Geom_BezierSurface::InsertPoleColAfter` without weights
    /// (`Geom_BezierSurface.cxx:652-685`). `v_index` is 1-based.
    /// `myMaxDerivInvOk` (`:684`) is not stored.
    pub fn insert_pole_col_after(
        &mut self,
        v_index: i32,
        col: &[GpPnt],
    ) -> Result<(), &'static str> {
        if v_index < 1 || v_index > self.n_v as i32 {
            return Err("Geom_BezierSurface::InsertPoleColAfter");
        }
        if col.len() != self.n_u {
            return Err("Geom_BezierSurface::InsertPoleColAfter");
        }
        let weights = if self.u_rational || self.v_rational {
            Some(vec![1.0; self.n_u])
        } else {
            None
        };
        self.insert_col_after(v_index, col, weights.as_deref(), false)
    }

    /// `Geom_BezierSurface::InsertPoleColAfter` with weights
    /// (`Geom_BezierSurface.cxx:689-734`).
    pub fn insert_pole_col_after_weights(
        &mut self,
        v_index: i32,
        col: &[GpPnt],
        col_weights: &[f64],
    ) -> Result<(), &'static str> {
        if v_index < 1 || v_index > self.n_v as i32 {
            return Err("Geom_BezierSurface::InsertPoleColAfter");
        }
        if col.len() != self.n_u || col_weights.len() != col.len() {
            return Err("Geom_BezierSurface::InsertPoleColAfter");
        }
        for &weight in col_weights {
            if weight <= occt_core::precision::REAL_SMALL {
                return Err("Geom_BezierSurface::InsertPoleColAfter");
            }
        }
        self.insert_col_after(v_index, col, Some(col_weights), true)
    }

    /// `Geom_BezierSurface::InsertPoleColBefore` (`cxx:738-751`).
    pub fn insert_pole_col_before(
        &mut self,
        v_index: i32,
        col: &[GpPnt],
    ) -> Result<(), &'static str> {
        self.insert_pole_col_after(v_index - 1, col)
    }

    /// Weighted `InsertPoleColBefore` (`cxx:746-751`).
    pub fn insert_pole_col_before_weights(
        &mut self,
        v_index: i32,
        col: &[GpPnt],
        col_weights: &[f64],
    ) -> Result<(), &'static str> {
        self.insert_pole_col_after_weights(v_index - 1, col, col_weights)
    }

    /// `Geom_BezierSurface::InsertPoleRowAfter` without weights
    /// (`Geom_BezierSurface.cxx:755-788`).
    pub fn insert_pole_row_after(
        &mut self,
        u_index: i32,
        row: &[GpPnt],
    ) -> Result<(), &'static str> {
        if u_index < 1 || u_index > self.n_u as i32 {
            return Err("Geom_BezierSurface::InsertPoleRowAfter");
        }
        if row.len() != self.n_v {
            return Err("Geom_BezierSurface::InsertPoleRowAfter");
        }
        let weights = if self.u_rational || self.v_rational {
            Some(vec![1.0; self.n_v])
        } else {
            None
        };
        self.insert_row_after(u_index, row, weights.as_deref(), false)
    }

    /// Weighted `InsertPoleRowAfter` (`cxx:792-837`).
    pub fn insert_pole_row_after_weights(
        &mut self,
        u_index: i32,
        row: &[GpPnt],
        row_weights: &[f64],
    ) -> Result<(), &'static str> {
        if u_index < 1 || u_index > self.n_u as i32 {
            return Err("Geom_BezierSurface::InsertPoleRowAfter");
        }
        if row.len() != self.n_v || row_weights.len() != row.len() {
            return Err("Geom_BezierSurface::InsertPoleRowAfter");
        }
        for &weight in row_weights {
            if weight <= occt_core::precision::REAL_SMALL {
                return Err("Geom_BezierSurface::InsertPoleRowAfter");
            }
        }
        self.insert_row_after(u_index, row, Some(row_weights), true)
    }

    /// `Geom_BezierSurface::InsertPoleRowBefore` (`cxx:841-854`).
    pub fn insert_pole_row_before(
        &mut self,
        u_index: i32,
        row: &[GpPnt],
    ) -> Result<(), &'static str> {
        self.insert_pole_row_after(u_index - 1, row)
    }

    /// Weighted `InsertPoleRowBefore` (`cxx:849-854`).
    pub fn insert_pole_row_before_weights(
        &mut self,
        u_index: i32,
        row: &[GpPnt],
        row_weights: &[f64],
    ) -> Result<(), &'static str> {
        self.insert_pole_row_after_weights(u_index - 1, row, row_weights)
    }

    /// `Geom_BezierSurface::RemovePoleCol` (`Geom_BezierSurface.cxx:858-895`).
    pub fn remove_pole_col(&mut self, v_index: i32) -> Result<(), &'static str> {
        if v_index < 1 || v_index > self.n_v as i32 {
            return Err("Geom_BezierSurface::RemovePoleCol");
        }
        if self.n_v <= 2 {
            return Err("Geom_BezierSurface::RemovePoleCol");
        }
        let nu = self.n_u;
        let nv = self.n_v;
        let new_nv = nv - 1;
        let mut poles = vec![GpPnt::new(0.0, 0.0, 0.0); nu * new_nv];
        let rational = self.u_rational || self.v_rational;
        let mut weights = rational.then(|| vec![0.0; nu * new_nv]);
        let remove = (v_index as usize) - 1;
        for i in 0..nu {
            let mut dst_j = 0usize;
            for j in 0..nv {
                if j == remove {
                    continue;
                }
                let dst = i * new_nv + dst_j;
                let src = i * nv + j;
                poles[dst] = self.poles[src];
                if let (Some(dst_w), Some(src_w)) = (weights.as_mut(), self.weights.as_ref()) {
                    dst_w[dst] = src_w[src];
                }
                dst_j += 1;
            }
        }
        self.poles = poles;
        self.n_v = new_nv;
        if rational {
            self.weights = weights;
            self.store_rational_flags();
        } else {
            self.weights = None;
        }
        Ok(())
    }

    /// `Geom_BezierSurface::RemovePoleRow` (`Geom_BezierSurface.cxx:899-937`).
    pub fn remove_pole_row(&mut self, u_index: i32) -> Result<(), &'static str> {
        if u_index < 1 || u_index > self.n_u as i32 {
            return Err("Geom_BezierSurface::RemovePoleRow");
        }
        if self.n_u <= 2 {
            return Err("Geom_BezierSurface::RemovePoleRow");
        }
        let nu = self.n_u;
        let nv = self.n_v;
        let new_nu = nu - 1;
        let mut poles = vec![GpPnt::new(0.0, 0.0, 0.0); new_nu * nv];
        let rational = self.u_rational || self.v_rational;
        let mut weights = rational.then(|| vec![0.0; new_nu * nv]);
        let remove = (u_index as usize) - 1;
        let mut dst_i = 0usize;
        for i in 0..nu {
            if i == remove {
                continue;
            }
            for j in 0..nv {
                let dst = dst_i * nv + j;
                let src = i * nv + j;
                poles[dst] = self.poles[src];
                if let (Some(dst_w), Some(src_w)) = (weights.as_mut(), self.weights.as_ref()) {
                    dst_w[dst] = src_w[src];
                }
            }
            dst_i += 1;
        }
        self.poles = poles;
        self.n_u = new_nu;
        if rational {
            self.weights = weights;
            self.store_rational_flags();
        } else {
            self.weights = None;
        }
        Ok(())
    }

    /// `after_index` is 1-based OCCT `VIndex`, or 0 for before the first column.
    fn insert_col_after(
        &mut self,
        after_index: i32,
        col: &[GpPnt],
        col_weights: Option<&[f64]>,
        recompute_flags: bool,
    ) -> Result<(), &'static str> {
        let nu = self.n_u;
        let nv = self.n_v;
        let new_nv = nv + 1;
        let insert_at = after_index.max(0) as usize;
        if insert_at > nv {
            return Err("Geom_BezierSurface::InsertPoleColAfter");
        }
        let mut poles = vec![GpPnt::new(0.0, 0.0, 0.0); nu * new_nv];
        let need_weights = recompute_flags || self.u_rational || self.v_rational;
        if need_weights && self.weights.is_none() {
            self.weights = Some(vec![1.0; nu * nv]);
        }
        let mut weights = need_weights.then(|| vec![0.0; nu * new_nv]);
        for i in 0..nu {
            let mut dst_j = 0usize;
            for j in 0..=nv {
                let dst = i * new_nv + dst_j;
                if j == insert_at {
                    poles[dst] = col[i];
                    if let Some(dst_w) = weights.as_mut() {
                        dst_w[dst] = col_weights.map(|w| w[i]).unwrap_or(1.0);
                    }
                } else {
                    let src_j = if j < insert_at { j } else { j - 1 };
                    let src = i * nv + src_j;
                    poles[dst] = self.poles[src];
                    if let (Some(dst_w), Some(src_w)) = (weights.as_mut(), self.weights.as_ref()) {
                        dst_w[dst] = src_w[src];
                    }
                }
                dst_j += 1;
            }
        }
        self.poles = poles;
        self.n_v = new_nv;
        if need_weights {
            self.weights = weights;
            if recompute_flags {
                self.store_rational_flags();
            } else {
                // Polynomial insert into a rational surface keeps the flags
                // (`cxx:668-677` does not call `Rational`).
            }
        } else {
            self.weights = None;
        }
        Ok(())
    }

    /// `Geom_BezierSurface::SetPole` (`Geom_BezierSurface.cxx:1047-1056`).
    pub fn set_pole(&mut self, u_index: i32, v_index: i32, point: GpPnt) -> Result<(), &'static str> {
        if u_index < 1
            || u_index > self.n_u as i32
            || v_index < 1
            || v_index > self.n_v as i32
        {
            return Err("Geom_BezierSurface::SetPole");
        }
        let slot = ((u_index as usize) - 1) * self.n_v + (v_index as usize) - 1;
        self.poles[slot] = point;
        Ok(())
    }

    /// `Geom_BezierSurface::SetPole` with weight (`cxx:1060-1078`).
    pub fn set_pole_weight(
        &mut self,
        u_index: i32,
        v_index: i32,
        point: GpPnt,
        weight: f64,
    ) -> Result<(), &'static str> {
        if weight <= occt_core::precision::REAL_SMALL {
            return Err("Geom_BezierSurface::SetPole");
        }
        self.set_pole(u_index, v_index, point)?;
        self.set_weight(u_index, v_index, weight)
    }

    /// `Geom_BezierSurface::SetPoleCol` with weights (`cxx:1082-1104`).
    /// `lower` is the OCCT array lower for both poles and weights.
    pub fn set_pole_col_weights(
        &mut self,
        v_index: i32,
        lower: i32,
        col: &[GpPnt],
        col_weights: &[f64],
    ) -> Result<(), &'static str> {
        if v_index < 1 || v_index > self.n_v as i32 {
            return Err("Geom_BezierSurface::SetPoleCol");
        }
        if col.len() != col_weights.len() {
            return Err("Geom_BezierSurface::SetPoleCol");
        }
        let upper = if col.is_empty() {
            lower - 1
        } else {
            lower + col.len() as i32 - 1
        };
        if lower < 1
            || lower > self.n_u as i32
            || upper < 1
            || upper > self.n_u as i32
        {
            return Err("Geom_BezierSurface::SetPoleCol");
        }
        let mut i = lower;
        let mut src = 0usize;
        while i <= upper {
            let slot = ((i as usize) - 1) * self.n_v + (v_index as usize) - 1;
            self.poles[slot] = col[src];
            i += 1;
            src += 1;
        }
        self.set_weight_col(v_index, lower, col_weights)
    }

    /// `Geom_BezierSurface::SetPoleCol` without weights (`cxx:1108-1125`).
    pub fn set_pole_col(
        &mut self,
        v_index: i32,
        lower: i32,
        col: &[GpPnt],
    ) -> Result<(), &'static str> {
        if v_index < 1 || v_index > self.n_v as i32 {
            return Err("Geom_BezierSurface::SetPoleCol");
        }
        let upper = if col.is_empty() {
            lower - 1
        } else {
            lower + col.len() as i32 - 1
        };
        if lower < 1
            || lower > self.n_u as i32
            || upper < 1
            || upper > self.n_u as i32
        {
            return Err("Geom_BezierSurface::SetPoleCol");
        }
        let mut i = lower;
        let mut src = 0usize;
        while i <= upper {
            let slot = ((i as usize) - 1) * self.n_v + (v_index as usize) - 1;
            self.poles[slot] = col[src];
            i += 1;
            src += 1;
        }
        Ok(())
    }

    /// `Geom_BezierSurface::SetPoleRow` without weights (`cxx:1129-1147`).
    pub fn set_pole_row(
        &mut self,
        u_index: i32,
        lower: i32,
        row: &[GpPnt],
    ) -> Result<(), &'static str> {
        if u_index < 1 || u_index > self.n_u as i32 {
            return Err("Geom_BezierSurface::SetPoleRow");
        }
        let upper = if row.is_empty() {
            lower - 1
        } else {
            lower + row.len() as i32 - 1
        };
        if lower < 1
            || lower > self.n_v as i32
            || upper < 1
            || upper > self.n_v as i32
        {
            return Err("Geom_BezierSurface::SetPoleRow");
        }
        let mut i = lower;
        let mut src = 0usize;
        while i <= upper {
            let slot = ((u_index as usize) - 1) * self.n_v + (i as usize) - 1;
            self.poles[slot] = row[src];
            i += 1;
            src += 1;
        }
        Ok(())
    }

    /// `Geom_BezierSurface::SetPoleRow` with weights (`cxx:1151-1175`).
    pub fn set_pole_row_weights(
        &mut self,
        u_index: i32,
        lower: i32,
        row: &[GpPnt],
        row_weights: &[f64],
    ) -> Result<(), &'static str> {
        if u_index < 1 || u_index > self.n_u as i32 {
            return Err("Geom_BezierSurface::SetPoleRow");
        }
        if row.len() != row_weights.len() {
            return Err("Geom_BezierSurface::SetPoleRow");
        }
        let upper = if row.is_empty() {
            lower - 1
        } else {
            lower + row.len() as i32 - 1
        };
        if lower < 1
            || lower > self.n_v as i32
            || upper < 1
            || upper > self.n_v as i32
        {
            return Err("Geom_BezierSurface::SetPoleRow");
        }
        let mut i = lower;
        let mut src = 0usize;
        while i <= upper {
            let slot = ((u_index as usize) - 1) * self.n_v + (i as usize) - 1;
            self.poles[slot] = row[src];
            i += 1;
            src += 1;
        }
        self.set_weight_row(u_index, lower, row_weights)
    }

    fn insert_row_after(
        &mut self,
        after_index: i32,
        row: &[GpPnt],
        row_weights: Option<&[f64]>,
        recompute_flags: bool,
    ) -> Result<(), &'static str> {
        let nu = self.n_u;
        let nv = self.n_v;
        let new_nu = nu + 1;
        let insert_at = after_index.max(0) as usize;
        if insert_at > nu {
            return Err("Geom_BezierSurface::InsertPoleRowAfter");
        }
        let mut poles = vec![GpPnt::new(0.0, 0.0, 0.0); new_nu * nv];
        let need_weights = recompute_flags || self.u_rational || self.v_rational;
        if need_weights && self.weights.is_none() {
            self.weights = Some(vec![1.0; nu * nv]);
        }
        let mut weights = need_weights.then(|| vec![0.0; new_nu * nv]);
        let mut dst_i = 0usize;
        for i in 0..=nu {
            if i == insert_at {
                for j in 0..nv {
                    let dst = dst_i * nv + j;
                    poles[dst] = row[j];
                    if let Some(dst_w) = weights.as_mut() {
                        dst_w[dst] = row_weights.map(|w| w[j]).unwrap_or(1.0);
                    }
                }
            } else {
                let src_i = if i < insert_at { i } else { i - 1 };
                for j in 0..nv {
                    let dst = dst_i * nv + j;
                    let src = src_i * nv + j;
                    poles[dst] = self.poles[src];
                    if let (Some(dst_w), Some(src_w)) = (weights.as_mut(), self.weights.as_ref()) {
                        dst_w[dst] = src_w[src];
                    }
                }
            }
            dst_i += 1;
        }
        self.poles = poles;
        self.n_u = new_nu;
        if need_weights {
            self.weights = weights;
            if recompute_flags {
                self.store_rational_flags();
            }
        } else {
            self.weights = None;
        }
        Ok(())
    }
}
