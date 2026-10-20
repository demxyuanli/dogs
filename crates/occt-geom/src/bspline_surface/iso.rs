//! `Geom_BSplineSurface::UIso` / `VIso`.

use super::*;
use crate::bspline_curve::GeomBSplineCurve;
use occt_core::bspl::eval;

impl GeomBSplineSurface {
    /// `Geom_BSplineSurface::UIso` (`Geom_BSplineSurface_1.cxx:598-635`).
    pub fn u_iso(&self, u: f64) -> Result<GeomBSplineCurve, String> {
        let nu = self.nb_poles_u();
        let nv = self.nb_poles_v();
        if nu == 0 || nv == 0 {
            return Err("GeomBSplineSurface::u_iso: empty pole grid".into());
        }
        let mut cpoles = Vec::with_capacity(nv);
        match &self.weights {
            Some(w) => {
                let mut cweights = Vec::with_capacity(nv);
                for j in 0..nv {
                    let col: Vec<GpPnt> = (0..nu).map(|i| self.poles[i][j]).collect();
                    let cw: Vec<f64> = (0..nu).map(|i| w[i][j]).collect();
                    cpoles.push(eval::eval_curve_rational(
                        &col,
                        &cw,
                        &self.knots_u,
                        self.deg_u,
                        u,
                    ));
                    let wpts: Vec<GpPnt> = cw.iter().map(|&wi| GpPnt::new(wi, 0.0, 0.0)).collect();
                    cweights.push(eval::eval_curve(&wpts, &self.knots_u, self.deg_u, u).x());
                }
                GeomBSplineCurve::rational(cpoles, cweights, self.knots_v.clone(), self.deg_v)
                    .map_err(|e| e.to_string())
            }
            None => {
                for j in 0..nv {
                    let col: Vec<GpPnt> = (0..nu).map(|i| self.poles[i][j]).collect();
                    cpoles.push(eval::eval_curve(&col, &self.knots_u, self.deg_u, u));
                }
                GeomBSplineCurve::new(cpoles, self.knots_v.clone(), self.deg_v)
                    .map_err(|e| e.to_string())
            }
        }
    }

    /// `Geom_BSplineSurface::VIso` (`Geom_BSplineSurface_1.cxx:775-812`).
    pub fn v_iso(&self, v: f64) -> Result<GeomBSplineCurve, String> {
        let nu = self.nb_poles_u();
        let nv = self.nb_poles_v();
        if nu == 0 || nv == 0 {
            return Err("GeomBSplineSurface::v_iso: empty pole grid".into());
        }
        let mut cpoles = Vec::with_capacity(nu);
        match &self.weights {
            Some(w) => {
                let mut cweights = Vec::with_capacity(nu);
                for i in 0..nu {
                    let row = &self.poles[i];
                    let rw = &w[i];
                    cpoles.push(eval::eval_curve_rational(
                        row,
                        rw,
                        &self.knots_v,
                        self.deg_v,
                        v,
                    ));
                    let wpts: Vec<GpPnt> = rw.iter().map(|&wi| GpPnt::new(wi, 0.0, 0.0)).collect();
                    cweights.push(eval::eval_curve(&wpts, &self.knots_v, self.deg_v, v).x());
                }
                GeomBSplineCurve::rational(cpoles, cweights, self.knots_u.clone(), self.deg_u)
                    .map_err(|e| e.to_string())
            }
            None => {
                for i in 0..nu {
                    cpoles.push(eval::eval_curve(&self.poles[i], &self.knots_v, self.deg_v, v));
                }
                GeomBSplineCurve::new(cpoles, self.knots_u.clone(), self.deg_u)
                    .map_err(|e| e.to_string())
            }
        }
    }
}
