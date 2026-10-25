//! `GeomFill_NSections::D0` / `D1` / `D2` and the shared `ResultEval` helper.
//!
//! Source: `GeomFill_NSections.cxx:147-198` (`ResultEval`), `:288-311` (`D0`),
//! `:315-390` (`D1`) and `:394-480` (`D2`) (OCCT 8.0.0).
//!
//! Translation notes:
//! - `ResultEval` packs the V pole rows as homogeneous coordinates
//!   `(x, y, z[, w])` with `x, y, z` multiplied by `w` (`cxx:168-187`) and
//!   evaluates them with `BSplCLib::Eval(V, false, deriv, Cdeg, ..., Cdim)`.
//!   `BSplCLib::Eval` with `ArrayDimension > 1` is componentwise, so the port
//!   evaluates each of the `Cdim` scalar components with
//!   [`occt_core::bspl::eval_extrap::eval_1d`] and takes `Result[deriv]`.
//! - `surf->VKnotSequence()` is the flat V knot vector ([`GeomBSplineSurface`]'s
//!   `knots_v`). A V-periodic surface is copied and made non-periodic before the
//!   evaluation, exactly as `cxx:341-355` / `cxx:431-445` do.
//! - `DPoles(ii).SetX/SetY/SetZ` writes OCCT's `gp_Vec` fields; the port builds
//!   a [`GpVec`] per pole.
//! - OCCT's `bool` results become `Err` (see the module docs of
//!   [`super`]).

use occt_core::gp::{GpPnt, GpVec};
use occt_core::precision::PCONFUSION;

use crate::bspline_surface::GeomBSplineSurface;
use crate::geom_fill::n_sections::NSections;

impl NSections {
    /// `GeomFill_NSections::D0` (`cxx:288-311`): the `VIso(V)` poles and their
    /// weights.
    pub(super) fn poles_at(&self, v: f64) -> Result<(Vec<GpPnt>, Vec<f64>), String> {
        let surface = self
            .surface
            .as_ref()
            .ok_or_else(|| "GeomFill_NSections::D0: null surface".to_string())?;
        let curve = surface.v_iso(v)?;
        // `Geom_BSplineCurve::WeightsArray()` (`Geom_BSplineCurve.hxx:805`): the
        // stored array, all ones when the curve is not rational.
        let weights = match &curve.weights {
            Some(w) => w.clone(),
            None => vec![1.0; curve.poles.len()],
        };
        Ok((curve.poles.clone(), weights))
    }

    /// `GeomFill_NSections::D1` (`cxx:315-390`).
    pub(super) fn poles_and_derivative_at(
        &self,
        v: f64,
    ) -> Result<(Vec<GpPnt>, Vec<GpVec>, Vec<f64>, Vec<f64>), String> {
        let (poles, weights) = self.poles_at(v)?;
        let surface = self
            .surface
            .as_ref()
            .ok_or_else(|| "GeomFill_NSections::D1: null surface".to_string())?;
        let rational = surface.is_v_rational();
        let gap = if rational { 4 } else { 3 };

        // `cxx:341-347`: a V-periodic surface is evaluated on a non-periodic
        // copy.
        let copy = if surface.v_periodic() {
            let mut copy = surface.as_ref().clone();
            copy.set_v_not_periodic();
            Some(copy)
        } else {
            None
        };
        let evaluated = copy.as_ref().unwrap_or(surface.as_ref());
        let result = result_eval(evaluated, v, 1)?;

        // `cxx:360-361`.
        let eps_w = 10.0 * PCONFUSION;
        let mut d_poles = Vec::with_capacity(poles.len());
        let mut d_weights = vec![0.0; poles.len()];
        for (i, pole) in poles.iter().enumerate() {
            let base = i * gap;
            let mut derivative = GpVec::new(result[base], result[base + 1], result[base + 2]);
            if rational {
                let ww = weights[i];
                if ww < eps_w {
                    // `cxx:377-380`: the derivative cannot be dehomogenised.
                    return Err("GeomFill_NSections::D1: null weight".to_string());
                }
                d_weights[i] = result[base + 3];
                // `cxx:384-385`: `(DPoles - DWeights * Poles) / ww`.
                derivative = GpVec::from_xyz(
                    &derivative
                        .xyz()
                        .subtracted(&pole.coord.multiplied(d_weights[i]))
                        .divided(ww),
                );
            }
            d_poles.push(derivative);
        }
        Ok((poles, d_poles, weights, d_weights))
    }

    /// `GeomFill_NSections::D2` (`cxx:394-480`).
    pub(super) fn poles_and_second_derivative_at(
        &self,
        v: f64,
    ) -> Result<
        (
            Vec<GpPnt>,
            Vec<GpVec>,
            Vec<GpVec>,
            Vec<f64>,
            Vec<f64>,
            Vec<f64>,
        ),
        String,
    > {
        let (poles, d_poles, weights, d_weights) = self.poles_and_derivative_at(v)?;
        let surface = self
            .surface
            .as_ref()
            .ok_or_else(|| "GeomFill_NSections::D2: null surface".to_string())?;
        // `cxx:408-413`: `BSplCLib::Eval` cannot do a 2nd derivative of a
        // rational degree-1 surface in V.
        if surface.deg_v() < 2 {
            return Err("GeomFill_NSections::D2: V degree < 2".to_string());
        }
        let rational = surface.is_v_rational();
        let gap = if rational { 4 } else { 3 };

        let copy = if surface.v_periodic() {
            let mut copy = surface.as_ref().clone();
            copy.set_v_not_periodic();
            Some(copy)
        } else {
            None
        };
        let evaluated = copy.as_ref().unwrap_or(surface.as_ref());
        let result = result_eval(evaluated, v, 2)?;

        let eps_w = 10.0 * PCONFUSION;
        let mut d2_poles = Vec::with_capacity(poles.len());
        let mut d2_weights = vec![0.0; poles.len()];
        for (i, pole) in poles.iter().enumerate() {
            let base = i * gap;
            let mut second = GpVec::new(result[base], result[base + 1], result[base + 2]);
            if rational {
                let ww = weights[i];
                if ww < eps_w {
                    return Err("GeomFill_NSections::D2: null weight".to_string());
                }
                d2_weights[i] = result[base + 3];
                // `cxx:471-474`: `(D2Poles - D2Weights * Poles
                // - 2 * DWeights * DPoles) / ww`.
                let linear = pole
                    .coord
                    .multiplied(d2_weights[i])
                    .added(&d_poles[i].xyz().multiplied(2.0 * d_weights[i]));
                second = GpVec::from_xyz(&second.xyz().subtracted(&linear).divided(ww));
            }
            d2_poles.push(second);
        }
        Ok((poles, d_poles, d2_poles, weights, d_weights, d2_weights))
    }
}

/// `ResultEval(surf, V, deriv, Result)` (`cxx:147-198`).
///
/// Returns the `deriv`-th derivative of the homogeneous `Cdim`-dimensional pole
/// vector at `V`.
fn result_eval(surface: &GeomBSplineSurface, v: f64, deriv: i32) -> Result<Vec<f64>, String> {
    // `cxx:153`: `surf->IsVRational()`.
    let rational = surface.is_v_rational();
    let gap = if rational { 4 } else { 3 };
    let n_u = surface.nb_poles_u();
    let n_v = surface.nb_poles_v();
    let dim = n_u * gap;

    // `cxx:168-187`: the poles in `BSplCLib::Eval` order, `jj` (V) outer and
    // `ii` (U) inner, with the coordinates weighted when rational.
    let mut components = vec![Vec::with_capacity(n_v); dim];
    for j in 0..n_v {
        for i in 0..n_u {
            let pole = surface.pole(i, j);
            let base = i * gap;
            let weighted = if rational {
                surface.weight(i, j)
            } else {
                1.0
            };
            components[base].push(pole.x() * weighted);
            components[base + 1].push(pole.y() * weighted);
            components[base + 2].push(pole.z() * weighted);
            if rational {
                components[base + 3].push(weighted);
            }
        }
    }

    let degree = surface.deg_v() as i32;
    let flat_knots = surface.knots_v();
    let mut result = vec![0.0; dim];
    for (index, component) in components.iter().enumerate() {
        let values = occt_core::bspl::eval_extrap::eval_1d(
            v,
            false,
            deriv,
            [degree, degree],
            degree,
            flat_knots,
            component,
        )
        .ok_or_else(|| "BSplCLib::Eval".to_string())?;
        result[index] = values[deriv as usize];
    }
    Ok(result)
}
