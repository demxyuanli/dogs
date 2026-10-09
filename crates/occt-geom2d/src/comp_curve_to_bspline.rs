//! `Geom2dConvert_CompCurveToBSplineCurve`
//! (`Geom2dConvert_CompCurveToBSplineCurve.cxx:28-244`): converts and
//! concatenates several 2D curves into one B-spline.
//!
//! Ported branches:
//!
//! * the two constructors (`cxx:28-54`), including the `Geom2d_BSplineCurve`
//!   down-cast / `CurveToBSplineCurve` fallback;
//! * the public `Add(NewCurve, Tolerance, After)` (`cxx:56-126`) with its
//!   pole-based G0 check and the before/after direction resolution;
//! * the private `Add(FirstCurve, SecondCurve, After)` (`cxx:129-230`): degree
//!   harmonisation, the C1 reparameterisation ratio, the knot/weight
//!   concatenation and the optional multiplicity reduction at the common knot.
//!
//! The 3D sibling (`GeomConvert_CompCurveToBSplineCurve`, ported in
//! `occt-geom/src/convert_bspl.rs`) carries extra `WithRatio`/`MinM` arguments
//! and an endpoint-based G0 check; the 2D `.cxx` has neither, so this port
//! follows the 2D control flow verbatim.
//!
//! `add_with_ratio` exposes those two 3D-sibling arguments to callers that
//! port 3D code operating on planar curves, where the 2D and 3D control flows
//! differ by exactly them. `ShapeFix_Wire::RemoveLoop` is such a caller: it
//! concatenates `Geom_TrimmedCurve(crv, ..)` segments with the 3D class at
//! `WithRatio = false` (`ShapeFix_Wire.cxx:2373`), and on the plane case
//! (`cxx:2386-2392`) the 3D curve is the edge curve whose 2D image is `c2d`
//! (`cxx:2393-2399`), so concatenating the 2D images with `WithRatio = false`
//! yields the exact 2D image of the 3D result.

use crate::bspline_curve::Geom2dBSplineCurve;
use crate::curve::Curve2d;
use crate::geom2d_convert::curve_to_bspline_curve_bspl;
use occt_core::convert::ParameterisationType;
use occt_core::gp::{GpPnt2d, GpPnt};

/// `Geom2d_BSplineCurve::Pole(i)` (1-based).
fn pole(curve: &Geom2dBSplineCurve, i: usize) -> GpPnt2d {
    GpPnt2d::new(curve.xs[i - 1], curve.ys[i - 1])
}

/// `Geom2d_BSplineCurve::Weight(i)` (1-based); a non-rational curve carries unit
/// weights (`Geom2d_BSplineCurve.cxx:194-206`).
fn weight(curve: &Geom2dBSplineCurve, i: usize) -> f64 {
    curve.weights.as_ref().map(|w| w[i - 1]).unwrap_or(1.0)
}

/// `Geom2d_BSplineCurve::IsRational()` for a weight array, i.e. `Rational()`
/// (`Geom2d_BSplineCurve.cxx:95-103`) evaluated by the rational constructor.
fn weights_are_rational(weights: &[f64]) -> bool {
    for w in weights.windows(2) {
        if (w[0] - w[1]).abs() > occt_core::precision::REAL_SMALL {
            return true;
        }
    }
    false
}

/// `Geom2dConvert_CompCurveToBSplineCurve`
/// (`Geom2dConvert_CompCurveToBSplineCurve.hxx:28-68`).
pub struct CompCurveToBSplineCurve {
    my_curve: Option<Geom2dBSplineCurve>,
    my_tol: f64,
    my_type: ParameterisationType,
}

impl CompCurveToBSplineCurve {
    /// `Geom2dConvert_CompCurveToBSplineCurve(Parameterisation)`
    /// (`cxx:28-33`): `myTol = Precision::Confusion()`.
    pub fn new(parameterisation: ParameterisationType) -> Self {
        Self {
            my_curve: None,
            my_tol: occt_core::precision::CONFUSION,
            my_type: parameterisation,
        }
    }

    /// `Geom2dConvert_CompCurveToBSplineCurve(BasisCurve, Parameterisation)`
    /// (`cxx:37-54`): a `Geom2d_BSplineCurve` basis is copied, anything else is
    /// converted with `Geom2dConvert::CurveToBSplineCurve(BasisCurve, myType)`.
    /// Returns `None` where the conversion is UNPORTED.
    pub fn from_curve(basis: &dyn Curve2d, parameterisation: ParameterisationType) -> Option<Self> {
        let mut s = Self::new(parameterisation);
        s.my_curve = Some(if basis.is_bspline2d() {
            basis.bspline_copy2d()?
        } else {
            curve_to_bspline_curve_bspl(basis, parameterisation)?
        });
        Some(s)
    }

    /// The `down_cast` branch of the constructor (`cxx:44-52`) for a caller that
    /// already owns a `Geom2d_BSplineCurve` (equivalent to `Copy()`).
    pub fn from_bspline(curve: Geom2dBSplineCurve, parameterisation: ParameterisationType) -> Self {
        Self {
            my_curve: Some(curve),
            my_tol: occt_core::precision::CONFUSION,
            my_type: parameterisation,
        }
    }

    /// `BSplineCurve()` (`cxx:232-235`).
    pub fn bspline_curve(&self) -> Option<&Geom2dBSplineCurve> {
        self.my_curve.as_ref()
    }

    /// `BSplineCurve()` by value.
    pub fn into_curve(self) -> Option<Geom2dBSplineCurve> {
        self.my_curve
    }

    /// `Clear()` (`cxx:239-243`).
    pub fn clear(&mut self) {
        self.my_curve = None;
    }

    /// `Add(NewCurve, Tolerance, After)` (`cxx:56-126`). Returns `false` where
    /// OCCT returns `false` (the curve is not G0 with the accumulated B-spline);
    /// `Err` where OCCT would raise through `CurveToBSplineCurve`.
    pub fn add(
        &mut self,
        new_curve: &dyn Curve2d,
        tolerance: f64,
        after: bool,
    ) -> Result<bool, &'static str> {
        self.add_with_ratio(new_curve, tolerance, after, true, 0)
    }

    /// `Add(NewCurve, Tolerance, After, WithRatio, MinM)` with the two extra
    /// arguments of the 3D sibling (`GeomConvert_CompCurveToBSplineCurve.hxx:56-60`).
    /// `add` is `add_with_ratio(.., true, 0)`.
    pub fn add_with_ratio(
        &mut self,
        new_curve: &dyn Curve2d,
        tolerance: f64,
        after: bool,
        with_ratio: bool,
        min_m: i32,
    ) -> Result<bool, &'static str> {
        // Conversion (`cxx:60-70`).
        let mut bs = if new_curve.is_bspline2d() {
            new_curve.bspline_copy2d().ok_or("CompCurveToBSplineCurve::add")?
        } else {
            curve_to_bspline_curve_bspl(new_curve, self.my_type)
                .ok_or("CompCurveToBSplineCurve::add")?
        };
        if self.my_curve.is_none() {
            self.my_curve = Some(bs);
            return Ok(true);
        }
        self.my_tol = tolerance;
        let sq_tol = tolerance * tolerance;

        // Pole-based G0 check (`cxx:80-92`).
        let current = self.my_curve.as_ref().unwrap();
        let (l_bs, l_cb) = (bs.nb_poles(), current.nb_poles());
        let mut d1 = pole(current, 1).square_distance(&pole(&bs, 1));
        let mut d2 = pole(current, 1).square_distance(&pole(&bs, l_bs));
        let is_before_reversed =
            (pole(current, 1).square_distance(&pole(&bs, 1)) < sq_tol) && (d1 < d2);
        let mut is_before =
            (pole(current, 1).square_distance(&pole(&bs, l_bs)) < sq_tol) || is_before_reversed;

        d1 = pole(current, l_cb).square_distance(&pole(&bs, 1));
        d2 = pole(current, l_cb).square_distance(&pole(&bs, l_bs));
        let is_after_reversed =
            (pole(current, l_cb).square_distance(&pole(&bs, l_bs)) < sq_tol) && (d2 < d1);
        let mut is_after =
            (pole(current, l_cb).square_distance(&pole(&bs, 1)) < sq_tol) || is_after_reversed;

        // `myCurve` and `NewCurve` together form a closed curve (`cxx:94-104`).
        if is_before && is_after {
            if after {
                is_before = false;
            } else {
                is_after = false;
            }
        }

        if is_after {
            // Append after (`cxx:105-114`).
            if is_after_reversed {
                bs.reverse();
            }
            let mut first = self.my_curve.take().unwrap();
            self.add_pair(&mut first, &mut bs, true, with_ratio, min_m)?;
            Ok(true)
        } else if is_before {
            // Prepend before (`cxx:115-122`).
            if is_before_reversed {
                bs.reverse();
            }
            let mut first = bs;
            let mut second = self.my_curve.take().unwrap();
            self.add_pair(&mut first, &mut second, false, with_ratio, min_m)?;
            Ok(true)
        } else {
            Ok(false)
        }
    }

    /// The private `Add(FirstCurve, SecondCurve, After, WithRatio, MinM)`
    /// (`cxx:129-230`, 3D `cxx:129-260`): harmonise the degrees,
    /// reparameterize onto the common knot, concatenate the poles/weights and
    /// lower the common knot's multiplicity down to `min_m`.
    /// Sets `myCurve`.
    fn add_pair(
        &mut self,
        first: &mut Geom2dBSplineCurve,
        second: &mut Geom2dBSplineCurve,
        after: bool,
        with_ratio: bool,
        min_m: i32,
    ) -> Result<(), &'static str> {
        // Harmonize the degrees (`cxx:133-141`).
        let deg = first.degree().max(second.degree());
        if first.degree() < deg {
            first.increase_degree(deg)?;
        }
        if second.degree() < deg {
            second.increase_degree(deg)?;
        }

        // Reparameterization ratio (C1 if possible) (`cxx:156-165`), skipped
        // when `WithRatio` is off (`GeomConvert_CompCurveToBSplineCurve.cxx:163-177`).
        let mut ratio = 1.0f64;
        if with_ratio {
            let l1 = first.d1(first.last_parameter()).1.magnitude();
            let l2 = second.d1(second.first_parameter()).1.magnitude();
            if l1 > occt_core::precision::CONFUSION && l2 > occt_core::precision::CONFUSION {
                ratio = l1 / l2;
            }
            if ratio < occt_core::precision::CONFUSION
                || ratio > 1.0 / occt_core::precision::CONFUSION
            {
                ratio = 1.0;
            }
        }

        let (uknots_f, umults_f) = first.distinct_knots_and_mults();
        let (uknots_s, umults_s) = second.distinct_knots_and_mults();
        let nb_p1 = first.nb_poles();
        let nb_p2 = second.nb_poles();
        let nb_k1 = uknots_f.len();
        let nb_k2 = uknots_s.len();
        if nb_k1 == 0 || nb_k2 == 0 || nb_p1 == 0 || nb_p2 == 0 {
            return Err("CompCurveToBSplineCurve::add_pair");
        }

        let (ratio1, delta1, ratio2, delta2, u_de_raccord);
        if after {
            // Do not move the first curve (`cxx:167-175`).
            ratio1 = 1.0;
            delta1 = 0.0;
            ratio2 = 1.0 / ratio;
            delta2 = ratio2 * uknots_s[0] - uknots_f[nb_k1 - 1];
            u_de_raccord = first.last_parameter();
        } else {
            // Do not move the second curve (`cxx:176-184`).
            ratio1 = ratio;
            delta1 = ratio1 * uknots_f[nb_k1 - 1] - uknots_s[0];
            ratio2 = 1.0;
            delta2 = 0.0;
            u_de_raccord = second.first_parameter();
        }

        // The knots (`cxx:186-201`).
        let total_knots = nb_k1 + nb_k2 - 1;
        let mut noeuds = vec![0.0f64; total_knots];
        let mut mults_out = vec![0i32; total_knots];
        for ii in 1..nb_k1 {
            noeuds[ii - 1] = ratio1 * uknots_f[ii - 1] - delta1;
            mults_out[ii - 1] = umults_f[ii - 1];
        }
        noeuds[nb_k1 - 1] = u_de_raccord;
        mults_out[nb_k1 - 1] = first.degree() as i32;
        for ii in 2..=nb_k2 {
            let jj = nb_k1 + ii - 1;
            noeuds[jj - 1] = ratio2 * uknots_s[ii - 1] - delta2;
            mults_out[jj - 1] = umults_s[ii - 1];
        }

        // The poles and weights (`cxx:202-218`).
        let weight_ratio = weight(first, nb_p1) / weight(second, 1);
        let mut poles_out: Vec<GpPnt> = Vec::with_capacity(nb_p1 + nb_p2 - 1);
        let mut weights_out: Vec<f64> = Vec::with_capacity(nb_p1 + nb_p2 - 1);
        for ii in 1..nb_p1 {
            poles_out.push(GpPnt::new(pole(first, ii).x(), pole(first, ii).y(), 0.0));
            weights_out.push(weight(first, ii));
        }
        for ii in 1..=nb_p2 {
            poles_out.push(GpPnt::new(pole(second, ii).x(), pole(second, ii).y(), 0.0));
            weights_out.push(weight_ratio * weight(second, ii));
        }

        // Create the BSpline (`cxx:220`): the rational constructor; its
        // `Rational()` check decides `myRational` (`Geom2d_BSplineCurve.cxx:193`).
        let xs: Vec<f64> = poles_out.iter().map(|p| p.x()).collect();
        let ys: Vec<f64> = poles_out.iter().map(|p| p.y()).collect();
        let flat = occt_core::bspl::banded_interp::knot_sequence(&noeuds, &mults_out, deg as i32);
        let mut curve = if weights_are_rational(&weights_out) {
            Geom2dBSplineCurve::from_flat(xs, ys, Some(weights_out), flat, deg, false)
                .map_err(|_| "CompCurveToBSplineCurve::add_pair")?
        } else {
            Geom2dBSplineCurve::from_flat(xs, ys, None, flat, deg, false)
                .map_err(|_| "CompCurveToBSplineCurve::add_pair")?
        };

        // Optionally reduce multiplicity down to `MinM` (`cxx:222-229`,
        // `GeomConvert_CompCurveToBSplineCurve.cxx:252-259`).
        let mut ok = true;
        let mut m = mults_out[nb_k1 - 1];
        while m > min_m && ok {
            m -= 1;
            ok = curve
                .remove_knot(nb_k1 as i32, m, self.my_tol)
                .map_err(|_| "CompCurveToBSplineCurve::add_pair")?;
        }

        self.my_curve = Some(curve);
        Ok(())
    }
}
