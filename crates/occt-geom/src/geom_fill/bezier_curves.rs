//! `GeomFill_BezierCurves` — the Bezier surface bounded by two, three or four
//! curves.
//!
//! Source: `GeomFill_BezierCurves.cxx:35-452` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `Geom_BezierCurve::Reversed()` allocates a copy; the port calls
//!   [`GeomBezierCurve::reverse`] on a cloned curve.
//! - `Poles()` / `WeightsArray()` are the 1-based `NCollection_Array1` views;
//!   the port reads [`GeomBezierCurve::poles`] and [`GeomBezierCurve::weights`]
//!   directly (a polynomial curve answers unit weights where OCCT would
//!   dereference a null handle — `SetSameWeights` is only reached when a curve
//!   is rational).
//! - The `NCollection_Array2<gp_Pnt>` result is flattened `u`-major into
//!   [`GeomBezierSurface::new`], matching `poles[i * n_v + j]`.

use occt_core::gp::GpPnt;
use occt_core::precision::CONFUSION;

use crate::bezier_curve::GeomBezierCurve;
use crate::bezier_surface::GeomBezierSurface;
use crate::curve::Curve;
use crate::geom_fill::coons::Coons;
use crate::geom_fill::curved::Curved;
use crate::geom_fill::filling::Filling;
use crate::geom_fill::stretch::Stretch;

use super::bspline_curves::FillingStyle;

/// `GeomFill_BezierCurves` (`GeomFill_BezierCurves.hxx:32-84`).
#[derive(Clone, Default)]
pub struct BezierCurves {
    surface: Option<GeomBezierSurface>,
}

fn start_point(curve: &GeomBezierCurve) -> GpPnt {
    curve.d0(0.0)
}

fn end_point(curve: &GeomBezierCurve) -> GpPnt {
    curve.d0(1.0)
}

fn reversed_curve(curve: &GeomBezierCurve) -> GeomBezierCurve {
    let mut copy = curve.clone();
    copy.reverse();
    copy
}

/// `SetSameWeights` (`GeomFill_BezierCurves.cxx:46-84`): scale the weights so
/// the four corners share one weight, which does not change the geometry.
fn set_same_weights(w1: &mut [f64], w2: &mut [f64], w3: &mut [f64], w4: &mut [f64]) {
    if w1.is_empty() || w2.is_empty() || w3.is_empty() || w4.is_empty() {
        return;
    }
    let nu = w1.len();
    let nv = w2.len();
    let a = (w1[0] * w2[0]) / (w1[nu - 1] * w2[nv - 1]);
    let b = (w3[0] * w4[0]) / (w3[nu - 1] * w4[nv - 1]);
    // `cxx:57-71`.
    let alfa = w1[nu - 1] / w2[0];
    for w in w2.iter_mut() {
        *w *= alfa;
    }
    let beta = w2[nv - 1] / w3[nu - 1];
    for w in w3.iter_mut() {
        *w *= beta;
    }
    let gamma = w3[0] / w4[nv - 1];
    for w in w4.iter_mut() {
        *w *= gamma;
    }
    // `cxx:73-83`.
    if (a - b).abs() > CONFUSION {
        if nv < 2 {
            return;
        }
        let w = (w1[0] / w4[0]).powf(1.0 / (nv - 1) as f64);
        let mut x = w;
        for i in (0..nv - 1).rev() {
            w4[i] *= x;
            x *= w;
        }
    }
}

/// `Arrange` (`GeomFill_BezierCurves.cxx:106-156`).
fn arrange_4(
    c1: &GeomBezierCurve,
    c2: &GeomBezierCurve,
    c3: &GeomBezierCurve,
    c4: &GeomBezierCurve,
    tol: f64,
) -> Option<[GeomBezierCurve; 4]> {
    let mut gc = [c1.clone(), c2.clone(), c3.clone(), c4.clone()];
    for i in 1..=3usize {
        let mut trouve = false;
        let mut j = i;
        while j <= 3 && !trouve {
            if start_point(&gc[j]).distance(&end_point(&gc[i - 1])) < tol {
                gc.swap(i, j);
                trouve = true;
            } else if end_point(&gc[j]).distance(&end_point(&gc[i - 1])) < tol {
                gc[j] = reversed_curve(&gc[j]);
                gc.swap(i, j);
                trouve = true;
            }
            j += 1;
        }
        if !trouve {
            return None;
        }
    }
    let cc3 = reversed_curve(&gc[2]);
    let cc4 = reversed_curve(&gc[3]);
    Some([gc[0].clone(), gc[1].clone(), cc3, cc4])
}

fn weights_of(curve: &GeomBezierCurve) -> Vec<f64> {
    match &curve.weights {
        Some(w) => w.clone(),
        None => vec![1.0; curve.poles.len()],
    }
}

/// `GeomFill_Filling` result assembled from the selected style.
fn build_filling(
    style: FillingStyle,
    p: [&[GpPnt]; 4],
    w: Option<[Vec<f64>; 4]>,
) -> Result<Box<dyn Filling>, String> {
    match (style, w) {
        (FillingStyle::Stretch, Some(w)) => Ok(Box::new(Stretch::new_rational(
            p[0], p[1], p[2], p[3], &w[0], &w[1], &w[2], &w[3],
        )?)),
        (FillingStyle::Stretch, None) => Ok(Box::new(Stretch::new(p[0], p[1], p[2], p[3])?)),
        // `cxx:260`: `GeomFill_Coons(P1, P4, P3, P2)`.
        (FillingStyle::Coons, Some(w)) => Ok(Box::new(Coons::new_rational(
            p[0], p[3], p[2], p[1], &w[0], &w[3], &w[2], &w[1],
        )?)),
        (FillingStyle::Coons, None) => Ok(Box::new(Coons::new(p[0], p[3], p[2], p[1])?)),
        (FillingStyle::Curved, Some(w)) => Ok(Box::new(Curved::new_rational(
            p[0], p[1], p[2], p[3], &w[0], &w[1], &w[2], &w[3],
        )?)),
        (FillingStyle::Curved, None) => Ok(Box::new(Curved::new(p[0], p[1], p[2], p[3])?)),
    }
}

/// `NCollection_Array2<gp_Pnt>` → the flattened `u`-major Bezier surface grid.
fn flatten(poles: &[Vec<GpPnt>], weights: Option<&[Vec<f64>]>) -> (Vec<GpPnt>, Vec<f64>, usize, usize) {
    let n_u = poles.len();
    let n_v = poles.first().map_or(0, |r| r.len());
    let mut flat = Vec::with_capacity(n_u * n_v);
    let mut flat_w = Vec::with_capacity(n_u * n_v);
    for row in poles.iter() {
        for p in row.iter() {
            flat.push(*p);
        }
    }
    if let Some(w) = weights {
        for row in w.iter() {
            for weight in row.iter() {
                flat_w.push(*weight);
            }
        }
    }
    (flat, flat_w, n_u, n_v)
}

impl BezierCurves {
    /// `GeomFill_BezierCurves(C1, C2, C3, C4, Type)`
    /// (`GeomFill_BezierCurves.cxx:192-300`).
    pub fn new_4(
        c1: &GeomBezierCurve,
        c2: &GeomBezierCurve,
        c3: &GeomBezierCurve,
        c4: &GeomBezierCurve,
        style: FillingStyle,
    ) -> Result<Self, String> {
        let tol = CONFUSION;
        let mut curves = arrange_4(c1, c2, c3, c4, tol).ok_or_else(|| {
            "Standard_ConstructionError: GeomFill_BezierCurves: Courbes non jointives".to_string()
        })?;
        // `cxx:212-218`.
        let mut deg_u = curves[0].degree().max(curves[2].degree());
        let mut deg_v = curves[1].degree().max(curves[3].degree());
        if style == FillingStyle::Coons {
            deg_u = deg_u.max(3);
            deg_v = deg_v.max(3);
        }
        // `cxx:220-236`.
        for index in [0usize, 2] {
            if curves[index].degree() < deg_u {
                let target = deg_u as i32;
                curves[index]
                    .increase(target)
                    .map_err(|e| format!("Standard_ConstructionError: {e}"))?;
            }
        }
        for index in [1usize, 3] {
            if curves[index].degree() < deg_v {
                let target = deg_v as i32;
                curves[index]
                    .increase(target)
                    .map_err(|e| format!("Standard_ConstructionError: {e}"))?;
            }
        }
        // `cxx:238-247`.
        let p1 = curves[0].poles.clone();
        let p2 = curves[1].poles.clone();
        let p3 = curves[2].poles.clone();
        let p4 = curves[3].poles.clone();
        let is_rat = curves[0].is_rational()
            || curves[1].is_rational()
            || curves[2].is_rational()
            || curves[3].is_rational();
        let weights = if is_rat {
            // `cxx:248-252`: corner weights are made consistent before the
            // filling sees them.
            let mut w1 = weights_of(&curves[0]);
            let mut w2 = weights_of(&curves[1]);
            let mut w3 = weights_of(&curves[2]);
            let mut w4 = weights_of(&curves[3]);
            set_same_weights(&mut w1, &mut w2, &mut w3, &mut w4);
            Some([w1, w2, w3, w4])
        } else {
            None
        };
        // `cxx:253-282`.
        let filling = build_filling(style, [&p1, &p2, &p3, &p4], weights)?;
        // `cxx:284-300`.
        let (poles, flat_w, n_u, n_v) = flatten(filling.poles(), filling.weights());
        let surface = if filling.is_rational() {
            GeomBezierSurface::with_weights(poles, n_u, n_v, flat_w)?
        } else {
            GeomBezierSurface::new(poles, n_u, n_v)?
        };
        Ok(Self {
            surface: Some(surface),
        })
    }

    /// `GeomFill_BezierCurves(C1, C2, C3, Type)`
    /// (`GeomFill_BezierCurves.cxx:304-337`): `C4` is the degree-1 segment
    /// between the two corner poles.
    pub fn new_3(
        c1: &GeomBezierCurve,
        c2: &GeomBezierCurve,
        c3: &GeomBezierCurve,
        style: FillingStyle,
    ) -> Result<Self, String> {
        let tol = CONFUSION * CONFUSION;
        // `cxx:311-332`.
        let pole1 = if start_point(c1).square_distance(&start_point(c2)) > tol
            && start_point(c1).square_distance(&end_point(c2)) > tol
        {
            start_point(c1)
        } else {
            end_point(c1)
        };
        let pole2 = if start_point(c3).square_distance(&start_point(c2)) > tol
            && start_point(c3).square_distance(&end_point(c2)) > tol
        {
            start_point(c3)
        } else {
            end_point(c3)
        };
        // `cxx:334-336`: `C4 = new Geom_BezierCurve(Poles)`.
        let c4 = GeomBezierCurve::new(vec![pole1, pole2])
            .map_err(|e| format!("Standard_ConstructionError: GeomFill_BezierCurves ({e})"))?;
        Self::new_4(c1, c2, c3, &c4, style)
    }

    /// `GeomFill_BezierCurves(C1, C2, Type)`
    /// (`GeomFill_BezierCurves.cxx:341-452`).
    pub fn new_2(
        c1: &GeomBezierCurve,
        c2: &GeomBezierCurve,
        style: FillingStyle,
    ) -> Result<Self, String> {
        let mut cc1 = c1.clone();
        let mut cc2 = c2.clone();
        let deg1 = cc1.degree();
        let deg2 = cc2.degree();
        let is_rat = cc1.is_rational() || cc2.is_rational();
        if style != FillingStyle::Curved {
            // `cxx:354-391`.
            let deg_u = deg1.max(deg2);
            if cc1.degree() < deg_u {
                let target = deg_u as i32;
                cc1.increase(target)
                    .map_err(|e| format!("Standard_ConstructionError: {e}"))?;
            }
            if cc2.degree() < deg_u {
                let target = deg_u as i32;
                cc2.increase(target)
                    .map_err(|e| format!("Standard_ConstructionError: {e}"))?;
            }
            let p1 = &cc1.poles;
            let p2 = &cc2.poles;
            let nb_u = deg_u + 1;
            let mut poles = Vec::with_capacity(nb_u * 2);
            let mut weights = Vec::with_capacity(nb_u * 2);
            let w1 = weights_of(&cc1);
            let w2 = weights_of(&cc2);
            for i in 0..nb_u {
                poles.push(p1[i]);
                poles.push(p2[i]);
                weights.push(w1[i]);
                weights.push(w2[i]);
            }
            let surface = if is_rat {
                GeomBezierSurface::with_weights(poles, nb_u, 2, weights)?
            } else {
                GeomBezierSurface::new(poles, nb_u, 2)?
            };
            return Ok(Self {
                surface: Some(surface),
            });
        }
        // `cxx:393-452`: `GeomFill_CurvedStyle`.
        let eps = CONFUSION;
        let mut ok = false;
        if start_point(&cc1).distance(&start_point(&cc2)) <= eps {
            ok = true;
        } else if start_point(&cc1).distance(&end_point(&cc2)) <= eps {
            cc2.reverse();
            ok = true;
        } else if end_point(&cc1).distance(&start_point(&cc2)) <= eps {
            // `cxx:407`: OCCT reverses the caller's `C1`; the port reverses the
            // local copy so the argument stays untouched.
            cc1.reverse();
            ok = true;
        } else if end_point(&cc1).distance(&end_point(&cc2)) <= eps {
            cc1.reverse();
            cc2.reverse();
            ok = true;
        }
        if !ok {
            return Err(
                "Standard_OutOfRange: GeomFill_BezierCurves: Courbes non jointives".into(),
            );
        }
        let (p1, p2) = (cc1.poles.clone(), cc2.poles.clone());
        let filling: Box<dyn Filling> = if is_rat {
            let w1 = weights_of(&cc1);
            let w2 = weights_of(&cc2);
            Box::new(Curved::new_two_rational(&p1, &p2, &w1, &w2)?)
        } else {
            Box::new(Curved::new_two(&p1, &p2)?)
        };
        let (poles, flat_w, n_u, n_v) = flatten(filling.poles(), filling.weights());
        let surface = if filling.is_rational() {
            GeomBezierSurface::with_weights(poles, n_u, n_v, flat_w)?
        } else {
            GeomBezierSurface::new(poles, n_u, n_v)?
        };
        Ok(Self {
            surface: Some(surface),
        })
    }

    /// `GeomFill_BezierCurves::Surface` (`GeomFill_BezierCurves.hxx:66`).
    pub fn surface(&self) -> Option<&GeomBezierSurface> {
        self.surface.as_ref()
    }

    pub fn into_surface(self) -> Result<GeomBezierSurface, String> {
        self.surface
            .ok_or_else(|| "GeomFill_BezierCurves: no surface".to_string())
    }
}
