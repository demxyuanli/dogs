//! `Geom_BSplineSurface::ExchangeUV` — the transposition
//! `GeomFill_GordonBuilder::Perform` applies to the guide surface
//! (`GeomFill_GordonBuilder.cxx:120`).
//!
//! Source: `ModelingData/TKG3d/Geom/Geom_BSplineSurface.cxx:317-361`
//! (OCCT 8.0.0).
//!
//! Translation notes:
//! - `myPoles` is stored with rows = U and columns = V (`Pole(UIndex, VIndex)
//!   = myPoles.Value(UIndex, VIndex)`, `Geom_BSplineSurface_1.cxx:573-579`),
//!   so the port transposes `poles[u][v]` into `poles[v][u]`.
//! - `cxx:326-342` copies the weights only when one of the rational flags is
//!   set; the port keys that on `weights.is_some()`, which is the same
//!   condition (a non-rational surface keeps no weight array).
//! - `cxx:355-358` swaps the rational, periodic, degree and knot arrays; the
//!   rational flags are implied by the weights here, so transposing the weight
//!   grid covers them, and the two `update*Knots()` calls follow (`cxx:359-360`).

use super::*;
use occt_core::gp::GpPnt;

impl GeomBSplineSurface {
    /// `Geom_BSplineSurface::ExchangeUV()` (`Geom_BSplineSurface.cxx:317-361`).
    pub fn exchange_uv(&mut self) {
        let n_u = self.poles.len();
        let n_v = self.poles.first().map_or(0, |r| r.len());
        let mut poles = vec![vec![GpPnt::zero(); n_u]; n_v];
        for u in 0..n_u {
            for v in 0..n_v {
                poles[v][u] = self.poles[u][v];
            }
        }
        if let Some(w) = &self.weights {
            let mut nw = vec![vec![1.0f64; n_u]; n_v];
            for u in 0..n_u {
                for v in 0..n_v {
                    nw[v][u] = w[u][v];
                }
            }
            self.weights = Some(nw);
        }
        self.poles = poles;

        std::mem::swap(&mut self.deg_u, &mut self.deg_v);
        std::mem::swap(&mut self.knots_u, &mut self.knots_v);
        std::mem::swap(&mut self.u_periodic, &mut self.v_periodic);
        self.update_u_knots();
        self.update_v_knots();
    }
}
