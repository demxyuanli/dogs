//! `GeomFill_Coons` — the bicubically blended Coons filling.
//!
//! Source: `GeomFill_Coons.cxx:24-201` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `PLib::CoefficientsPoles` is `occt_core::bspl::plib::coefficients_poles`
//!   (`PLib.cxx:1482-1493`, `:1522-1589`).
//! - `BSplCLib::IncreaseDegree(NewDegree, Poles, Weights, NewPoles,
//!   NewWeights)` (`BSplCLib_BzSyntaxes.cxx:31-43`, the Bezier elevation of
//!   `BSplCLib_CurveComputation.pxx:2080-2099`) is reached through
//!   `GeomBezierCurve::increase` (`Geom_BezierCurve.cxx:141-192`), which calls
//!   the same `BSplCLib::IncreaseDegree` with the two-knot Bezier vector.

use occt_core::bspl::plib::coefficients_poles;
use occt_core::gp::GpPnt;

use crate::bezier_curve::GeomBezierCurve;
use crate::geom_fill::filling::Filling;
use crate::geom_fill::stretch::coord;

/// `GeomFill_Coons` (`GeomFill_Coons.hxx:20-70`).
#[derive(Debug, Clone, Default)]
pub struct Coons {
    poles: Vec<Vec<GpPnt>>,
    weights: Option<Vec<Vec<f64>>>,
    rational: bool,
}

/// The degree-3 blending coefficients `Coef` of `GeomFill_Coons::Init`
/// (`cxx:88-91`) converted from power to Bezier form (`cxx:92`).
fn blending_poles() -> Vec<GpPnt> {
    let coef = [
        GpPnt::new(1.0, 0.0, 0.0),
        GpPnt::new(0.0, 0.0, 0.0),
        GpPnt::new(-3.0, 3.0, 0.0),
        GpPnt::new(2.0, -2.0, 0.0),
    ];
    coefficients_poles(&coef)
}

/// `cxx:94-108`: `CoefU` / `CoefV`, the blending coefficients pulled to the
/// pole count of the direction (`BSplCLib::IncreaseDegree` when it exceeds 4).
fn blending_poles_for(n_poles: usize) -> Result<Vec<GpPnt>, String> {
    let pole = blending_poles();
    if n_poles > 4 {
        let mut curve = GeomBezierCurve::new(pole)
            .map_err(|e| format!("Standard_ConstructionError: GeomFill_Coons ({e})"))?;
        curve
            .increase((n_poles - 1) as i32)
            .map_err(|e| format!("Standard_ConstructionError: GeomFill_Coons ({e})"))?;
        Ok(curve.poles)
    } else {
        Ok(pole)
    }
}

impl Coons {
    /// `GeomFill_Coons(P1, P2, P3, P4)` / `Init(P1, P2, P3, P4)`
    /// (`GeomFill_Coons.cxx:53-142`).
    pub fn new(p1: &[GpPnt], p2: &[GpPnt], p3: &[GpPnt], p4: &[GpPnt]) -> Result<Self, String> {
        // `cxx:59`.
        if p1.len() != p3.len() || p2.len() != p4.len() {
            return Err("Standard_DomainError: GeomFill_Coons".into());
        }
        let npol_u = p1.len();
        let npol_v = p2.len();
        let mut poles = vec![vec![GpPnt::zero(); npol_v]; npol_u];
        // `cxx:70-80`: the four boundaries.
        for i in 0..npol_u {
            poles[i][0] = p1[i];
            poles[i][npol_v - 1] = p3[i];
        }
        for i in 0..npol_v {
            poles[0][i] = p2[i];
            poles[npol_u - 1][i] = p4[i];
        }
        // `cxx:82-108`.
        let coef_u = blending_poles_for(npol_u)?;
        let coef_v = blending_poles_for(npol_v)?;
        // `cxx:109-120`: `CoefU(i).Coord(FU(i), GU(i), Dummy)`.
        let mut fu = vec![0.0; npol_u];
        let mut gu = vec![0.0; npol_u];
        for i in 2..npol_u {
            fu[i - 1] = coord(&coef_u[i - 1], 0);
            gu[i - 1] = coord(&coef_u[i - 1], 1);
        }
        let mut fv = vec![0.0; npol_v];
        let mut gv = vec![0.0; npol_v];
        for i in 2..npol_v {
            fv[i - 1] = coord(&coef_v[i - 1], 0);
            gv[i - 1] = coord(&coef_v[i - 1], 1);
        }
        // `cxx:122-141`: the interior poles.
        for j in 1..npol_v - 1 {
            for i in 1..npol_u - 1 {
                let value = |k: usize| {
                    fv[j] * coord(&poles[i][0], k)
                        + gv[j] * coord(&poles[i][npol_v - 1], k)
                        + fu[i] * coord(&poles[0][j], k)
                        + gu[i] * coord(&poles[npol_u - 1][j], k)
                        - fu[i] * fv[j] * coord(&poles[0][0], k)
                        - fu[i] * gv[j] * coord(&poles[0][npol_v - 1], k)
                        - gu[i] * fv[j] * coord(&poles[npol_u - 1][0], k)
                        - gu[i] * gv[j] * coord(&poles[npol_u - 1][npol_v - 1], k)
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

    /// `GeomFill_Coons(P1, P2, P3, P4, W1, W2, W3, W4)` /
    /// `Init(..., W1, W2, W3, W4)` (`GeomFill_Coons.cxx:145-201`).
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
        // `cxx:156-159`.
        if w1.len() != w3.len() || w2.len() != w4.len() {
            return Err("Standard_DomainError: GeomFill_Coons".into());
        }
        if w1.len() != p1.len() || w2.len() != p2.len() || w3.len() != p3.len() || w4.len() != p4.len()
        {
            return Err("Standard_DomainError: GeomFill_Coons".into());
        }
        let mut result = Self::new(p1, p2, p3, p4)?;
        result.rational = true;
        let npol_u = w1.len();
        let npol_v = w2.len();
        let nu = (npol_u - 1) as f64;
        let nv = (npol_v - 1) as f64;
        let mut weights = vec![vec![0.0; npol_v]; npol_u];
        // `cxx:173-178`.
        for i in 0..npol_u {
            weights[i][0] = w1[i];
            weights[i][npol_v - 1] = w3[i];
        }
        // `cxx:180-201`.
        for j in 1..npol_v - 1 {
            let pv = j as f64 / nv;
            let pv1 = 1.0 - pv;
            weights[0][j] = w4[j];
            weights[npol_u - 1][j] = w2[j];
            for i in 1..npol_u - 1 {
                let pu = i as f64 / nu;
                let pu1 = 1.0 - pu;
                weights[i][j] = pv1 * w1[i] + pv * w3[i] + pu * w2[j] + pu1 * w4[j]
                    - (pu1 * pv1 * w1[0]
                        + pu * pv1 * w2[0]
                        + pu * pv * w3[npol_u - 1]
                        + pu1 * pv * w4[npol_v - 1]);
            }
        }
        result.weights = Some(weights);
        Ok(result)
    }
}

impl Filling for Coons {
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
