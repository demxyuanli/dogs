//! `GeomFill_Stretch` — the bilinear (transfinite "stretch") filling.
//!
//! Source: `GeomFill_Stretch.cxx:26-162` (OCCT 8.0.0).

use occt_core::gp::GpPnt;

use crate::geom_fill::filling::Filling;

/// `GeomFill_Stretch` (`GeomFill_Stretch.hxx:20-60`).
#[derive(Debug, Clone, Default)]
pub struct Stretch {
    poles: Vec<Vec<GpPnt>>,
    weights: Option<Vec<Vec<f64>>>,
    rational: bool,
}

impl Stretch {
    /// `GeomFill_Stretch(P1, P2, P3, P4)` / `Init(P1, P2, P3, P4)`
    /// (`GeomFill_Stretch.cxx:53-102`).
    pub fn new(
        p1: &[GpPnt],
        p2: &[GpPnt],
        p3: &[GpPnt],
        p4: &[GpPnt],
    ) -> Result<Self, String> {
        // `cxx:59`: `P1.Length() != P3.Length() || P2.Length() != P4.Length()`.
        if p1.len() != p3.len() || p2.len() != p4.len() {
            return Err("Standard_DomainError: GeomFill_Stretch".into());
        }
        let npol_u = p1.len();
        let npol_v = p2.len();
        let nu = (npol_u - 1) as f64;
        let nv = (npol_v - 1) as f64;
        let mut poles = vec![vec![GpPnt::zero(); npol_v]; npol_u];
        // `cxx:70-76`: the P1 / P3 boundaries.
        for i in 0..npol_u {
            poles[i][0] = p1[i];
            poles[i][npol_v - 1] = p3[i];
        }
        // `cxx:78-101`.
        for j in 1..npol_v - 1 {
            let pv = j as f64 / nv;
            let pv1 = 1.0 - pv;
            poles[0][j] = p4[j];
            poles[npol_u - 1][j] = p2[j];
            for i in 1..npol_u - 1 {
                let pu = i as f64 / nu;
                let pu1 = 1.0 - pu;
                let value = |k: usize| {
                    pv1 * coord(&p1[i], k) + pv * coord(&p3[i], k) + pu * coord(&p2[j], k)
                        + pu1 * coord(&p4[j], k)
                        - (pu1 * pv1 * coord(&p1[0], k)
                            + pu * pv1 * coord(&p2[0], k)
                            + pu * pv * coord(&p3[npol_u - 1], k)
                            + pu1 * pv * coord(&p4[npol_v - 1], k))
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

    /// `GeomFill_Stretch(P1, P2, P3, P4, W1, W2, W3, W4)` /
    /// `Init(..., W1, W2, W3, W4)` (`GeomFill_Stretch.cxx:106-161`).
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
        // `cxx:115-118`.
        if w1.len() != w3.len() || w2.len() != w4.len() {
            return Err("Standard_DomainError: GeomFill_Stretch".into());
        }
        if w1.len() != p1.len() || w2.len() != p2.len() || w3.len() != p3.len() || w4.len() != p4.len()
        {
            return Err("Standard_DomainError: GeomFill_Stretch".into());
        }
        let mut result = Self::new(p1, p2, p3, p4)?;
        result.rational = true;
        let npol_u = w1.len();
        let npol_v = w2.len();
        let nu = (npol_u - 1) as f64;
        let nv = (npol_v - 1) as f64;
        let mut weights = vec![vec![0.0; npol_v]; npol_u];
        // `cxx:130-136`.
        for i in 0..npol_u {
            weights[i][0] = w1[i];
            weights[i][npol_v - 1] = w3[i];
        }
        // `cxx:138-161`.
        for j in 1..npol_v - 1 {
            let pv = j as f64 / nv;
            let pv1 = 1.0 - pv;
            weights[0][j] = w4[j];
            weights[npol_u - 1][j] = w2[j];
            for i in 1..npol_u - 1 {
                let pu = i as f64 / nu;
                let pu1 = 1.0 - pu;
                weights[i][j] = 0.5 * (pv1 * w1[i] + pv * w3[i] + pu * w2[j] + pu1 * w4[j]);
            }
        }
        result.weights = Some(weights);
        Ok(result)
    }
}

impl Filling for Stretch {
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

/// `Coord(K)` (`gp_Pnt.hxx:88-108`) with a 0-based `k`.
pub(crate) fn coord(p: &GpPnt, k: usize) -> f64 {
    match k {
        0 => p.x(),
        1 => p.y(),
        _ => p.z(),
    }
}

