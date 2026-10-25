//! `GeomFill_Curved` — the "curved" (half-weighted) filling style and the
//! two-curve extrusion used for a degenerate boundary.
//!
//! Source: `GeomFill_Curved.cxx:24-228` (OCCT 8.0.0).

use occt_core::gp::GpPnt;

use crate::geom_fill::filling::Filling;
use crate::geom_fill::stretch::coord;

/// `GeomFill_Curved` (`GeomFill_Curved.hxx:20-80`).
#[derive(Debug, Clone, Default)]
pub struct Curved {
    poles: Vec<Vec<GpPnt>>,
    weights: Option<Vec<Vec<f64>>>,
    rational: bool,
}

impl Curved {
    /// `GeomFill_Curved(P1, P2, P3, P4)` / `Init(P1, P2, P3, P4)`
    /// (`GeomFill_Curved.cxx:74-123`).
    pub fn new(p1: &[GpPnt], p2: &[GpPnt], p3: &[GpPnt], p4: &[GpPnt]) -> Result<Self, String> {
        // `cxx:79`.
        if p1.len() != p3.len() || p2.len() != p4.len() {
            return Err("Standard_DomainError: GeomFill_Curved".into());
        }
        let npol_u = p1.len();
        let npol_v = p2.len();
        let nu = (npol_u - 1) as f64;
        let nv = (npol_v - 1) as f64;
        let mut poles = vec![vec![GpPnt::zero(); npol_v]; npol_u];
        // `cxx:90-96`.
        for i in 0..npol_u {
            poles[i][0] = p1[i];
            poles[i][npol_v - 1] = p3[i];
        }
        // `cxx:98-122`.
        for j in 1..npol_v - 1 {
            let mut pv = j as f64 / nv;
            let mut pv1 = 1.0 - pv;
            pv /= 2.0;
            pv1 /= 2.0;
            poles[0][j] = p4[j];
            poles[npol_u - 1][j] = p2[j];
            for i in 1..npol_u - 1 {
                let mut pu = i as f64 / nu;
                let mut pu1 = 1.0 - pu;
                pu /= 2.0;
                pu1 /= 2.0;
                let value = |k: usize| {
                    pv1 * coord(&p1[i], k)
                        + pv * coord(&p3[i], k)
                        + pu * coord(&p2[j], k)
                        + pu1 * coord(&p4[j], k)
                };
                poles[i][j] = GpPnt::new(value(0), value(1), value(2));
            }
        }
        Ok(Self {
            poles,
            weights: None,
            rational: false,
        })
    }

    /// `GeomFill_Curved(P1, P2, P3, P4, W1, W2, W3, W4)`
    /// (`GeomFill_Curved.cxx:126-179`).
    pub fn new_rational(
        p1: &[GpPnt],
        p2: &[GpPnt],
        p3: &[GpPnt],
        p4: &[GpPnt],
        w1: &[f64],
        w2: &[f64],
        w3: &[f64],
        w4: &[f64],
    ) -> Result<Self, String> {
        // `cxx:137-140`.
        if w1.len() != w3.len() || w2.len() != w4.len() {
            return Err("Standard_DomainError: GeomFill_Curved".into());
        }
        if w1.len() != p1.len() || w2.len() != p2.len() || w3.len() != p3.len() || w4.len() != p4.len()
        {
            return Err("Standard_DomainError: GeomFill_Curved".into());
        }
        let mut result = Self::new(p1, p2, p3, p4)?;
        result.rational = true;
        let npol_u = w1.len();
        let npol_v = w2.len();
        let nu = (npol_u - 1) as f64;
        let nv = (npol_v - 1) as f64;
        let mut weights = vec![vec![0.0; npol_v]; npol_u];
        // `cxx:151-157`.
        for i in 0..npol_u {
            weights[i][0] = w1[i];
            weights[i][npol_v - 1] = w3[i];
        }
        // `cxx:159-178`.
        for j in 1..npol_v - 1 {
            let mut pv = j as f64 / nv;
            let mut pv1 = 1.0 - pv;
            pv /= 2.0;
            pv1 /= 2.0;
            weights[0][j] = w4[j];
            weights[npol_u - 1][j] = w2[j];
            for i in 1..npol_u - 1 {
                let mut pu = i as f64 / nu;
                let mut pu1 = 1.0 - pu;
                pu /= 2.0;
                pu1 /= 2.0;
                weights[i][j] = pv1 * w1[i] + pv * w3[i] + pu * w2[j] + pu1 * w4[j];
            }
        }
        result.weights = Some(weights);
        result.rational = true;
        Ok(result)
    }

    /// `GeomFill_Curved(P1, P2)` / `Init(P1, P2)` (`GeomFill_Curved.cxx:183-204`):
    /// the surface swept by translating `P1` along the `P2` poles.
    pub fn new_two(p1: &[GpPnt], p2: &[GpPnt]) -> Result<Self, String> {
        let npol_u = p1.len();
        let npol_v = p2.len();
        let mut poles = vec![vec![GpPnt::zero(); npol_v]; npol_u];
        // `cxx:198-202`: `gp_Vec Tra(P2(1), P2(j))`, then `P1(i).Translated`.
        for j in 0..npol_v {
            for i in 0..npol_u {
                poles[i][j] = GpPnt::new(
                    p1[i].x() + (p2[j].x() - p2[0].x()),
                    p1[i].y() + (p2[j].y() - p2[0].y()),
                    p1[i].z() + (p2[j].z() - p2[0].z()),
                );
            }
        }
        Ok(Self {
            poles,
            weights: None,
            rational: false,
        })
    }

    /// `GeomFill_Curved(P1, P2, W1, W2)` (`GeomFill_Curved.cxx:206-228`).
    pub fn new_two_rational(
        p1: &[GpPnt],
        p2: &[GpPnt],
        w1: &[f64],
        w2: &[f64],
    ) -> Result<Self, String> {
        let mut result = Self::new_two(p1, p2)?;
        result.rational = true;
        let npol_u = w1.len();
        let npol_v = w2.len();
        let mut weights = vec![vec![0.0; npol_v]; npol_u];
        // `cxx:220-227`: `Factor = W2(j) / W1(1)`.
        for j in 0..npol_v {
            let factor = w2[j] / w1[0];
            for i in 0..npol_u {
                weights[i][j] = w1[i] * factor;
            }
        }
        result.weights = Some(weights);
        Ok(result)
    }
}

impl Filling for Curved {
    fn nb_u_poles(&self) -> usize {
        self.poles.len()
    }
    fn nb_v_poles(&self) -> usize {
        self.poles.first().map_or(0, |row| row.len())
    }
    fn is_rational(&self) -> bool {
        self.rational
    }
    fn poles(&self) -> &[Vec<GpPnt>] {
        &self.poles
    }
    fn weights(&self) -> Option<&[Vec<f64>]> {
        self.weights.as_deref()
    }
}
