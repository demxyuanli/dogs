// ---------------------------------------------------------------------------
// CompCurveToBSplineCurve.
// ---------------------------------------------------------------------------

use occt_core::bspl::banded_interp::knot_sequence;
use occt_core::convert::{ConvertError, ParameterisationType};
use occt_core::gp::GpPnt;

use crate::bspline_curve::GeomBSplineCurve;
use crate::curve::Curve;

use super::curve_convert::{bspline_copy, curve_to_bspline_curve, is_bspline, weights_are_rational};

/// `GeomConvert_CompCurveToBSplineCurve`
/// (`GeomConvert_CompCurveToBSplineCurve.hxx:30-78`): converts and
/// concatenates several curves into one B-spline
/// (`GeomConvert_CompCurveToBSplineCurve.cxx:32-273`).
pub struct CompCurveToBSplineCurve {
    curve: Option<GeomBSplineCurve>,
    tol: f64,
    kind: ParameterisationType,
}

impl CompCurveToBSplineCurve {
    /// `GeomConvert_CompCurveToBSplineCurve(Parameterisation)`
    /// (`cxx:32-37`): `myTol = Precision::Confusion()`.
    pub fn new(parameterisation: ParameterisationType) -> Self {
        Self {
            curve: None,
            tol: occt_core::precision::CONFUSION,
            kind: parameterisation,
        }
    }

    /// `GeomConvert_CompCurveToBSplineCurve(BasisCurve, Parameterisation)`
    /// (`cxx:41-56`): a `Geom_BSplineCurve` basis is copied, anything else is
    /// converted with `GeomConvert::CurveToBSplineCurve(BasisCurve, myType)`.
    pub fn from_curve(
        basis: &dyn Curve,
        parameterisation: ParameterisationType,
    ) -> Result<Self, ConvertError> {
        let mut s = Self::new(parameterisation);
        s.curve = Some(if is_bspline(basis) {
            bspline_copy(basis)?
        } else {
            curve_to_bspline_curve(basis, parameterisation)?
        });
        Ok(s)
    }

    /// The constructor's `down_cast` branch (`cxx:47-51`) for a caller that
    /// already owns a `GeomBSplineCurve`; equivalent to `Copy()`.
    pub fn from_bspline(curve: GeomBSplineCurve, parameterisation: ParameterisationType) -> Self {
        Self {
            curve: Some(curve),
            tol: occt_core::precision::CONFUSION,
            kind: parameterisation,
        }
    }

    /// `BSplineCurve()` (`cxx:264-267`).
    pub fn bspline_curve(&self) -> Option<&GeomBSplineCurve> {
        self.curve.as_ref()
    }

    /// `Clear()` (`cxx:271-274`).
    pub fn clear(&mut self) {
        self.curve = None;
    }

    /// `BSplineCurve()` by value.
    pub fn into_curve(self) -> Option<GeomBSplineCurve> {
        self.curve
    }

    /// `Add(NewCurve, Tolerance, After, WithRatio, MinM)` (`cxx:60-131`).
    ///
    /// `After`/`WithRatio`/`MinM` reproduce the OCCT defaults explicitly
    /// (`GeomConvert_CompCurveToBSplineCurve.hxx:56-60`: `false` / `true` /
    /// `0`).
    pub fn add(
        &mut self,
        new_curve: &dyn Curve,
        tolerance: f64,
        after: bool,
        with_ratio: bool,
        min_m: i32,
    ) -> Result<bool, ConvertError> {
        // Conversion (`cxx:66-75`).
        let mut bs = if is_bspline(new_curve) {
            bspline_copy(new_curve)?
        } else {
            curve_to_bspline_curve(new_curve, self.kind)?
        };
        if self.curve.is_none() {
            self.curve = Some(bs);
            return Ok(true);
        }
        self.tol = tolerance;

        // Use actual curve endpoints instead of poles for the G0 continuity
        // check (`cxx:85-94`).
        let (curve_start, curve_end) = {
            let cur = self.curve.as_ref().unwrap();
            (cur.d0(cur.first_parameter()), cur.d0(cur.last_parameter()))
        };
        let bs_start = bs.d0(bs.first_parameter());
        let bs_end = bs.d0(bs.last_parameter());

        let mut avant = curve_start.distance(&bs_start) < self.tol
            || curve_start.distance(&bs_end) < self.tol;
        let mut apres = curve_end.distance(&bs_start) < self.tol
            || curve_end.distance(&bs_end) < self.tol;

        // Will myCurve be (or become) closed? Resolve the ambiguity
        // (`cxx:96-107`).
        if avant && apres {
            if after {
                avant = false;
            } else {
                apres = false;
            }
        }

        if apres {
            // Append after? (`cxx:109-118`).
            if curve_end.distance(&bs_end) < self.tol {
                bs.reverse();
            }
            let mut first = self.curve.take().unwrap();
            self.add_pair(&mut first, &mut bs, true, with_ratio, min_m)?;
            Ok(true)
        } else if avant {
            // Prepend before? (`cxx:119-128`).
            if curve_start.distance(&bs_start) < self.tol {
                bs.reverse();
            }
            let mut first = bs;
            let mut second = self.curve.take().unwrap();
            self.add_pair(&mut first, &mut second, false, with_ratio, min_m)?;
            Ok(true)
        } else {
            Ok(false)
        }
    }

    /// The private `Add(FirstCurve, SecondCurve, After, WithRatio, MinM)`
    /// (`cxx:135-260`): harmonize the degrees, reparameterize onto the common
    /// knot, concatenate the poles/weights and optionally lower the common
    /// knot's multiplicity down to `MinM`. Sets `myCurve`.
    fn add_pair(
        &mut self,
        first: &mut GeomBSplineCurve,
        second: &mut GeomBSplineCurve,
        after: bool,
        with_ratio: bool,
        min_m: i32,
    ) -> Result<(), ConvertError> {
        // Harmonize the degrees (`cxx:141-150`).
        let deg = first.degree().max(second.degree());
        if first.degree() < deg {
            first
                .increase_degree(deg)
                .map_err(|_| ConvertError::ConstructionError)?;
        }
        if second.degree() < deg {
            second
                .increase_degree(deg)
                .map_err(|_| ConvertError::ConstructionError)?;
        }

        // Reparameterization ratio (C1 if possible) (`cxx:163-177`).
        let mut ratio = 1.0f64;
        if with_ratio {
            let l1 = first.eval_dn(first.last_parameter(), 1).magnitude();
            let l2 = second.eval_dn(second.first_parameter(), 1).magnitude();
            if l1 > occt_core::precision::CONFUSION && l2 > occt_core::precision::CONFUSION {
                ratio = l1 / l2;
            }
            if ratio < occt_core::precision::CONFUSION
                || ratio > 1.0 / occt_core::precision::CONFUSION
            {
                ratio = 1.0;
            }
        }

        let (uf, umf) = first.distinct_knots_and_mults();
        let (us, ums) = second.distinct_knots_and_mults();
        let nb_p1 = first.nb_poles();
        let nb_p2 = second.nb_poles();
        let nb_k1 = uf.len();
        let nb_k2 = us.len();
        if nb_k1 == 0 || nb_k2 == 0 || nb_p1 == 0 || nb_p2 == 0 {
            return Err(ConvertError::ConstructionError);
        }

        let (ratio1, delta1, ratio2, delta2);
        if after {
            // Do not move the first curve (`cxx:179-186`).
            ratio1 = 1.0;
            delta1 = 0.0;
            ratio2 = 1.0 / ratio;
            delta2 = ratio2 * us[0] - uf[nb_k1 - 1];
        } else {
            // Do not move the second curve (`cxx:187-194`).
            ratio1 = ratio;
            delta1 = ratio1 * uf[nb_k1 - 1] - us[0];
            ratio2 = 1.0;
            delta2 = 0.0;
        }

        // The knots (`cxx:196-229`).
        let mut noeuds = vec![0.0f64; nb_k1 + nb_k2 - 1];
        let mut mults_out = vec![0i32; nb_k1 + nb_k2 - 1];
        for ii in 1..=nb_k1 {
            noeuds[ii - 1] = ratio1 * uf[ii - 1] - delta1;
            if ii > 1 {
                let mut eps = occt_core::precision::epsilon(noeuds[ii - 2].abs());
                if eps < 5.0e-10 {
                    eps = 5.0e-10;
                }
                if noeuds[ii - 1] - noeuds[ii - 2] <= eps {
                    noeuds[ii - 1] += eps;
                }
            }
            mults_out[ii - 1] = umf[ii - 1];
        }
        mults_out[nb_k1 - 1] = first.degree() as i32;
        for ii in 2..=nb_k2 {
            let jj = nb_k1 + ii - 1;
            noeuds[jj - 1] = ratio2 * us[ii - 1] - delta2;
            let mut eps = occt_core::precision::epsilon(noeuds[jj - 2].abs());
            if eps < 5.0e-10 {
                eps = 5.0e-10;
            }
            if noeuds[jj - 1] - noeuds[jj - 2] <= eps {
                noeuds[jj - 1] += eps;
            }
            mults_out[jj - 1] = ums[ii - 1];
        }

        // The poles and weights (`cxx:231-247`).
        let mut r = first.weight(nb_p1 as i32);
        r /= second.weight(1);
        let mut poles_out: Vec<GpPnt> = Vec::with_capacity(nb_p1 + nb_p2 - 1);
        let mut weights_out: Vec<f64> = Vec::with_capacity(nb_p1 + nb_p2 - 1);
        for ii in 1..nb_p1 {
            poles_out.push(*first.pole(ii - 1));
            weights_out.push(first.weight(ii as i32));
        }
        for ii in 1..=nb_p2 {
            poles_out.push(*second.pole(ii - 1));
            weights_out.push(r * second.weight(ii as i32));
        }

        // Create the BSpline (`cxx:249-250`): the rational constructor with
        // `CheckRational` decides `myRational` (`Geom_BSplineCurve.cxx:205-223`).
        let flat = knot_sequence(&noeuds, &mults_out, deg as i32);
        let mut curve = if weights_are_rational(&weights_out) {
            GeomBSplineCurve::rational(poles_out, weights_out, flat, deg)
                .map_err(|_| ConvertError::ConstructionError)?
        } else {
            GeomBSplineCurve::new(poles_out, flat, deg)
                .map_err(|_| ConvertError::ConstructionError)?
        };

        // Optionally reduce multiplicity down to MinM (`cxx:252-259`).
        let mut ok = true;
        let mut m = mults_out[nb_k1 - 1];
        while m > min_m && ok {
            m -= 1;
            ok = curve
                .remove_knot(nb_k1 as i32, m, self.tol)
                .map_err(|_| ConvertError::ConstructionError)?;
        }

        self.curve = Some(curve);
        Ok(())
    }
}
