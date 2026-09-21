//! Port of `Convert_ConicToBSplineCurve`
//! (`Convert_ConicToBSplineCurve.hxx:42-159`,
//! `Convert_ConicToBSplineCurve.cxx:32-785`) and of the
//! `Convert_ParameterisationType` enumeration
//! (`Convert_ParameterisationType.hxx:36-46`).
//!
//! `Convert_ConicToBSplineCurve` is the data holder plus the `BuildCosAndSin`
//! engine shared by `Convert_CircleToBSplineCurve`,
//! `Convert_EllipseToBSplineCurve`, `Convert_HyperbolaToBSplineCurve` and
//! `Convert_ParabolaToBSplineCurve` (ported in `convert/conic_curves.rs`). It
//! produces, for a requested parameterisation:
//!
//! - `cos`/`sin` **numerators** and a **denominator** table (the rational
//!   curve is `(cos / den, sin / den)`),
//! - the degree, the distinct knots and their multiplicities,
//! - whether the result is periodic.
//!
//! OCCT's 1-based `NCollection_Array1` indices are 0-based here; each loop
//! keeps the original expression, so `ii = 1; ii <= num_spans` becomes
//! `for ii in 1..=num_spans` with `arr[ii - 1]`, and the places where the
//! translation is not literal carry the `.cxx` line.
//!
//! **UNPORTED (returns `ConvertError::Unported`)**: the `Convert_Polynomial`
//! arm of the ranged `BuildCosAndSin` (`cxx:612-621`) calls
//! `BuildPolynomialCosAndSin` (`Convert_PolynomialCosAndSin.cxx:64-181`), whose
//! `Locate`/`BSplCLib::Trimming` machinery (`Convert_PolynomialCosAndSin.cxx:29-62`,
//! `BSplCLib_1.cxx:204-...`, 2D `gp_Pnt2d` overload) is not ported in this
//! repository; no caller in the repository requests `Convert_Polynomial`.
//! OCCT's `GeomConvert::CurveToBSplineCurve` default parameterisation is
//! `Convert_TgtThetaOver2` (`GeomConvert.hxx:270-272`), which **is** ported, and
//! `GeomToIGES_GeomCurve::TransferCurve(Geom_Ellipse)` uses
//! `Convert_QuasiAngular` (`GeomToIGES_GeomCurve.cxx:639`), also ported.

use crate::bspl::banded_interp::{
    eval_bspline_basis, interpolate_contact, knot_sequence, schoenberg_points,
};
use crate::bspl::plib_eval::eval_poly0;
use crate::gp::GpPnt2d;

/// `Convert_ParameterisationType` (`Convert_ParameterisationType.hxx:36-46`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ParameterisationType {
    TgtThetaOver2,
    TgtThetaOver2_1,
    TgtThetaOver2_2,
    TgtThetaOver2_3,
    TgtThetaOver2_4,
    QuasiAngular,
    RationalC1,
    Polynomial,
}

/// Errors of the conic conversions. OCCT signals these by throwing
/// (`Standard_ConstructionError` / `Standard_DomainError`); the two port-only
/// variants are documented at their raise sites.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConvertError {
    /// `Standard_ConstructionError` (`cxx:327`, `:391`, `:399`, `:642`).
    ConstructionError,
    /// `Standard_DomainError` (`Convert_CircleToBSplineCurve.cxx:129`,
    /// `Convert_EllipseToBSplineCurve.cxx:130-131`).
    DomainError,
    /// An OCCT arm whose dependency is not ported yet (see the module doc).
    Unported,
}

/// The five tables `BuildCosAndSin` returns through reference arguments
/// (`cxx:359-368`).
#[derive(Debug, Clone)]
pub struct CosAndSin {
    pub cos_numerator: Vec<f64>,
    pub sin_numerator: Vec<f64>,
    pub denominator: Vec<f64>,
    pub degree: usize,
    pub knots: Vec<f64>,
    pub mults: Vec<i32>,
}

/// Data holder of `Convert_ConicToBSplineCurve` (`cxx:32-151`): the BSpline
/// description of a conic, filled by the subclasses.
#[derive(Debug, Clone)]
pub struct ConicToBSplineCurve {
    poles: Vec<GpPnt2d>,
    weights: Vec<f64>,
    knots: Vec<f64>,
    mults: Vec<i32>,
    degree: usize,
    is_periodic: bool,
}

impl ConicToBSplineCurve {
    /// `Convert_ConicToBSplineCurve(NumberOfPoles, NumberOfKnots, Degree)`
    /// (`cxx:32-47`) followed by the constructor body of
    /// `Convert_{Circle,Ellipse,Hyperbola,Parabola}ToBSplineCurve`: the poles
    /// are supplied by the caller, the remaining tables come from
    /// `BuildCosAndSin`.
    pub(crate) fn from_parts(
        poles: Vec<GpPnt2d>,
        weights: Vec<f64>,
        knots: Vec<f64>,
        mults: Vec<i32>,
        degree: usize,
        is_periodic: bool,
    ) -> Self {
        Self { poles, weights, knots, mults, degree, is_periodic }
    }

    /// `Degree()` (`cxx:51-54`).
    pub fn degree(&self) -> usize {
        self.degree
    }

    /// `NbPoles()` (`cxx:58-61`).
    pub fn nb_poles(&self) -> usize {
        self.poles.len()
    }

    /// `NbKnots()` (`cxx:65-68`).
    pub fn nb_knots(&self) -> usize {
        self.knots.len()
    }

    /// `IsPeriodic()` (`cxx:72-75`).
    pub fn is_periodic(&self) -> bool {
        self.is_periodic
    }

    /// `Poles()` (`cxx:126-130`).
    pub fn poles(&self) -> &[GpPnt2d] {
        &self.poles
    }

    /// `Weights()` (`cxx:134-137`).
    pub fn weights(&self) -> &[f64] {
        &self.weights
    }

    /// `Knots()` (`cxx:141-144`).
    pub fn knots(&self) -> &[f64] {
        &self.knots
    }

    /// `Multiplicities()` (`cxx:148-151`).
    pub fn multiplicities(&self) -> &[i32] {
        &self.mults
    }
}

/// The two evaluators reachable from the ranged `BuildCosAndSin`
/// (`cxx:376`, set at `:557` and `:588`).
#[derive(Debug, Clone, Copy)]
enum Evaluator {
    RationalC1,
    QuasiAngular,
}

/// Basis functions of one `BSplCLib::D0` evaluation:
/// `BSplCLib::EvalBsplineBasis` (`BSplCLib_2.cxx:429-563`).
fn d0_basis(flat: &[f64], degree: usize, u: f64) -> Result<(usize, Vec<f64>), ConvertError> {
    eval_bspline_basis(0, degree + 1, flat, u).map_err(|_| ConvertError::ConstructionError)
}

/// `BSplCLib::D0` non-rational arm (`BSplCLib_1.cxx:248-267`): the value of the
/// scalar B-spline `values` at `u`.
fn d0_scalar(values: &[f64], flat: &[f64], degree: usize, u: f64) -> Result<f64, ConvertError> {
    let (first_nz, basis) = d0_basis(flat, degree, u)?;
    let mut sum = 0.0;
    for j in 0..=degree {
        sum += values.get(first_nz + j).copied().unwrap_or(0.0) * basis[j];
    }
    Ok(sum)
}

/// `BSplCLib::D0` rational arm (`BSplCLib_1.cxx:248-267` with `Weights`): the
/// weighted sum divided by the weight function (the same convention
/// `Geom_BSplineCurve::D0` / `Geom2d_BSplineCurve::D0` rely on).
fn d0_rational_scalar(
    numerator: &[f64],
    denominator: &[f64],
    flat: &[f64],
    degree: usize,
    u: f64,
) -> Result<f64, ConvertError> {
    let (first_nz, basis) = d0_basis(flat, degree, u)?;
    let mut num = 0.0;
    let mut den = 0.0;
    for j in 0..=degree {
        let b = basis[j];
        num += numerator.get(first_nz + j).copied().unwrap_or(0.0) * b;
        den += denominator.get(first_nz + j).copied().unwrap_or(0.0) * b;
    }
    Ok(num / den)
}

/// `BSplCLib::D0` (`BSplCLib_1.cxx:248-267`) for the 2D non-rational poles of an
/// evaluator (`CosAndSinRationalC1`, `cxx:242-253`).
fn d0_point2d(
    poles: &[GpPnt2d],
    flat: &[f64],
    degree: usize,
    u: f64,
) -> Result<[f64; 2], ConvertError> {
    let (first_nz, basis) = d0_basis(flat, degree, u)?;
    let mut out = [0.0f64; 2];
    for j in 0..=degree {
        let b = basis[j];
        if let Some(p) = poles.get(first_nz + j) {
            out[0] += p.x() * b;
            out[1] += p.y() * b;
        }
    }
    Ok(out)
}

/// `CosAndSinRationalC1` (`cxx:235-254`) and `CosAndSinQuasiAngular`
/// (`cxx:271-298`); returns OCCT's `(Result[0], Result[1])`.
fn eval_cos_and_sin(
    evaluator: Evaluator,
    parameter: f64,
    eval_degree: usize,
    eval_poles: &[GpPnt2d],
    eval_knots: &[f64],
    eval_mults: &[i32],
) -> Result<[f64; 2], ConvertError> {
    match evaluator {
        Evaluator::RationalC1 => {
            // `BSplCLib::D0(Parameter, 0, EvalDegree, false, EvalPoles,
            //  BSplCLib::NoWeights(), EvalKnots, EvalMults, a_point)` (`cxx:243-251`).
            let flat = knot_sequence(eval_knots, eval_mults, eval_degree as i32);
            d0_point2d(eval_poles, &flat, eval_degree, parameter)
        }
        Evaluator::QuasiAngular => {
            // `Standard_OutOfRange_Raise_if(EvalPoles.Length() != EvalDegree + 1, ...)`
            // (`cxx:280-281`).
            if eval_poles.len() != eval_degree + 1 {
                return Err(ConvertError::ConstructionError);
            }
            let mut coeffs = Vec::with_capacity(eval_poles.len() * 2);
            for p in eval_poles {
                coeffs.push(p.x());
                coeffs.push(p.y());
            }
            let mut out = [0.0f64; 2];
            // `PLib::NoDerivativeEvalPolynomial(param, EvalDegree, 2,
            //  EvalDegree << 1, aCoeffs(0), Result[0])` (`cxx:296-297`);
            // `param = Parameter * 0.5` (`cxx:296`).
            eval_poly0(&coeffs, eval_degree as i32, 2, parameter * 0.5, &mut out);
            Ok(out)
        }
    }
}

/// File-static `AlgorithmicCosAndSin` (`cxx:306-355`): samples the evaluator at
/// the Schoenberg points of the target knot vector, interpolates the
/// homogeneous poles `(V^2 - U^2, 2*V*U, V^2 + U^2)` and normalises.
#[allow(clippy::too_many_arguments)]
fn algorithmic_cos_and_sin(
    degree: usize,
    flat_knots: &[f64],
    eval_degree: usize,
    eval_poles: &[GpPnt2d],
    eval_knots: &[f64],
    eval_mults: &[i32],
    evaluator: Evaluator,
    cos_numerator: &mut [f64],
    sin_numerator: &mut [f64],
    denominator: &mut [f64],
) -> Result<(), ConvertError> {
    let order = degree + 1;
    let num_poles = flat_knots.len() - order;
    if num_poles != cos_numerator.len()
        || num_poles != sin_numerator.len()
        || num_poles != denominator.len()
    {
        // `throw Standard_ConstructionError()` (`cxx:324-328`).
        return Err(ConvertError::ConstructionError);
    }
    let parameters = schoenberg_points(degree, flat_knots, num_poles);
    // `poles_array` holds 3 homogeneous coordinates per pole (`cxx:330`, `:338-340`).
    let mut poles_array = vec![0.0f64; num_poles * 3];
    for (ii, &param) in parameters.iter().enumerate() {
        let result = eval_cos_and_sin(
            evaluator, param, eval_degree, eval_poles, eval_knots, eval_mults,
        )?;
        poles_array[ii * 3] = result[1] * result[1] - result[0] * result[0];
        poles_array[ii * 3 + 1] = 2.0 * result[1] * result[0];
        poles_array[ii * 3 + 2] = result[1] * result[1] + result[0] * result[0];
    }
    let contact_order_array = vec![0i32; num_poles];
    // `BSplCLib::Interpolate(Degree, FlatKnots, parameters,
    //  contact_order_array, poles_array, pivot_index_problem)` (`cxx:342-347`).
    interpolate_contact(
        degree,
        flat_knots,
        &parameters,
        &contact_order_array,
        &mut poles_array,
        3,
    )
    .map_err(|_| ConvertError::ConstructionError)?;
    for ii in 0..num_poles {
        let inverse = 1.0 / poles_array[ii * 3 + 2];
        cos_numerator[ii] = poles_array[ii * 3] * inverse;
        sin_numerator[ii] = poles_array[ii * 3 + 1] * inverse;
        denominator[ii] = poles_array[ii * 3 + 2];
    }
    Ok(())
}

/// `Convert_ConicToBSplineCurve::BuildCosAndSin(Parameterisation, UFirst,
/// ULast, ...)` (`cxx:359-622`).
///
/// The OCCT member reads no object state, so it is a free function here; the
/// outputs are returned in `CosAndSin`.
pub fn build_cos_and_sin(
    parameterisation: ParameterisationType,
    u_first: f64,
    u_last: f64,
) -> Result<CosAndSin, ConvertError> {
    let pi = std::f64::consts::PI;
    let delta = u_last - u_first;

    let mut degree = 0usize;
    let mut num_poles = 0usize;
    let mut num_knots = 1usize;
    // OCCT declares `num_spans`/`temp_degree` with a dummy initial value
    // (`cxx:373-374`); every `switch` arm assigns them before use.
    let num_spans: usize;
    let mut order = 0usize;
    let mut tgt_theta_flag = false;

    // `switch (Parameterisation)` (`cxx:380-434`).
    match parameterisation {
        ParameterisationType::TgtThetaOver2 => {
            num_spans = (1.2 * delta / pi).trunc() as usize + 1;
            tgt_theta_flag = true;
        }
        ParameterisationType::TgtThetaOver2_1 => {
            num_spans = 1;
            if delta > 0.9999 * pi {
                return Err(ConvertError::ConstructionError);
            }
            tgt_theta_flag = true;
        }
        ParameterisationType::TgtThetaOver2_2 => {
            num_spans = 2;
            if delta > 1.9999 * pi {
                return Err(ConvertError::ConstructionError);
            }
            tgt_theta_flag = true;
        }
        ParameterisationType::TgtThetaOver2_3 => {
            num_spans = 3;
            tgt_theta_flag = true;
        }
        ParameterisationType::TgtThetaOver2_4 => {
            num_spans = 4;
            tgt_theta_flag = true;
        }
        ParameterisationType::QuasiAngular => {
            num_poles = 7;
            degree = 6;
            num_spans = 1;
            num_knots = 2;
            order = degree + 1;
        }
        ParameterisationType::RationalC1 => {
            degree = 4;
            order = degree + 1;
            num_poles = 8;
            num_knots = 3;
            num_spans = 2;
        }
        ParameterisationType::Polynomial => {
            degree = 7;
            num_poles = 8;
            num_knots = 2;
            num_spans = 1;
        }
    }

    let mut alpha = 0.0f64;
    if tgt_theta_flag {
        // `cxx:435-440`.
        alpha = delta / (2.0 * num_spans as f64);
        degree = 2;
        num_poles = 2 * num_spans + 1;
    }

    let mut cos_numerator = vec![0.0f64; num_poles];
    let mut sin_numerator = vec![0.0f64; num_poles];
    let mut denominator = vec![0.0f64; num_poles];
    let mut knots = vec![0.0f64; num_spans + 1];
    let mut mults = vec![0i32; num_spans + 1];

    if tgt_theta_flag {
        // `cxx:447-471`: a degree-2 rational arc per span, middle weight
        // `1 / cos(alpha)`.
        let mut param = u_first;
        cos_numerator[0] = u_first.cos();
        sin_numerator[0] = u_first.sin();
        denominator[0] = 1.0;
        knots[0] = param;
        mults[0] = degree as i32 + 1;
        let direct = alpha.cos();
        let inverse = 1.0 / direct;
        for ii in 1..=num_spans {
            cos_numerator[2 * ii - 1] = inverse * (param + alpha).cos();
            sin_numerator[2 * ii - 1] = inverse * (param + alpha).sin();
            denominator[2 * ii - 1] = direct;
            cos_numerator[2 * ii] = (param + 2.0 * alpha).cos();
            sin_numerator[2 * ii] = (param + 2.0 * alpha).sin();
            denominator[2 * ii] = 1.0;
            knots[ii] = param + 2.0 * alpha;
            mults[ii] = 2;
            param += 2.0 * alpha;
        }
        mults[num_spans] = degree as i32 + 1;
    } else if parameterisation != ParameterisationType::Polynomial {
        // `cxx:472-611`.
        alpha = (u_last - u_first) * 0.5;
        let beta = (u_last + u_first) * 0.5;
        let cos_beta = beta.cos();
        let sin_beta = beta.sin();
        let num_flat_knots = num_poles + order;
        let mut flat_knots = vec![0.0f64; num_flat_knots];
        for ii in 0..order {
            flat_knots[ii] = -alpha;
            flat_knots[ii + num_poles] = alpha;
        }
        knots[0] = u_first;
        knots[num_knots - 1] = u_last;
        mults[0] = order as i32;
        mults[num_knots - 1] = order as i32;

        let temp_degree: usize;
        let mut temp_poles = vec![GpPnt2d::zero(); 4];
        let mut temp_knots: Vec<f64> = Vec::new();
        let mut temp_mults: Vec<i32> = Vec::new();
        let evaluator;
        match parameterisation {
            ParameterisationType::QuasiAngular => {
                // `cxx:502-558`: V(t) = t + c*t^3 in Coord(1), U(t) = 1 + b*t^2
                // in Coord(2).
                let alpha_2 = alpha * 0.5;
                let mut p_param = -1.0 / (alpha_2 * alpha_2);
                if alpha_2 < pi * 0.5 {
                    if alpha_2 < 1.0e-7 {
                        // Taylor value of b(gamma) for the 0/0 case (`cxx:532-535`).
                        p_param = -6.0 / 15.0;
                    } else {
                        let tan_alpha_2 = alpha_2.tan();
                        let value1 = alpha_2 / (3.0 * (tan_alpha_2 - alpha_2));
                        p_param += value1;
                    }
                }
                let q_param = (1.0 / 3.0) + p_param;
                temp_degree = 3;
                temp_poles[0].set_coord(0.0, 1.0);
                temp_poles[1].set_coord(1.0, 0.0);
                temp_poles[2].set_coord(0.0, p_param);
                temp_poles[3].set_coord(q_param, 0.0);
                evaluator = Evaluator::QuasiAngular;
            }
            ParameterisationType::RationalC1 => {
                // `cxx:559-589`.
                for ii in order..num_poles {
                    flat_knots[ii] = 0.0;
                }
                knots[1] = u_first + alpha;
                mults[1] = degree as i32 - 1;
                temp_degree = 2;
                let alpha_2 = alpha * 0.5;
                let alpha_4 = alpha * 0.25;
                let tan_alpha_2 = alpha_2.tan();
                // `jj` walks 1 -> 4 over the two passes of the OCCT loop
                // (`cxx:570-576`).
                temp_poles[1].set_y(1.0 + alpha_4 * tan_alpha_2);
                temp_poles[0].set_y(1.0);
                temp_poles[2].set_y(1.0 + alpha_4 * tan_alpha_2);
                temp_poles[3].set_y(1.0);
                temp_poles[0].set_x(-tan_alpha_2);
                temp_poles[1].set_x(alpha_4 - tan_alpha_2);
                temp_poles[2].set_x(-alpha_4 + tan_alpha_2);
                temp_poles[3].set_x(tan_alpha_2);
                temp_knots = vec![-alpha, 0.0, alpha];
                temp_mults = vec![(temp_degree + 1) as i32, 1, (temp_degree + 1) as i32];
                evaluator = Evaluator::RationalC1;
            }
            _ => return Err(ConvertError::ConstructionError),
        }
        algorithmic_cos_and_sin(
            degree,
            &flat_knots,
            temp_degree,
            &temp_poles,
            &temp_knots,
            &temp_mults,
            evaluator,
            &mut cos_numerator,
            &mut sin_numerator,
            &mut denominator,
        )?;
        // `cxx:604-610`: rotate the tables by `beta`.
        for ii in 0..num_poles {
            let value1 = cos_beta * cos_numerator[ii] - sin_beta * sin_numerator[ii];
            let value2 = sin_beta * cos_numerator[ii] + cos_beta * sin_numerator[ii];
            cos_numerator[ii] = value1;
            sin_numerator[ii] = value2;
        }
    } else {
        // UNPORTED: `Convert_Polynomial` (`cxx:612-621`) sets
        // `Knots(1) = 0`, `Knots(num_knots) = 1`, `Mults = num_poles` and calls
        // `BuildPolynomialCosAndSin(UFirst, ULast, num_poles, ...)`
        // (`Convert_PolynomialCosAndSin.cxx:64-181`), whose `Locate`
        // (`:29-62`) and `BSplCLib::Trimming(degree, false, knots, mults,
        // poles, NoWeights, trim_min, trim_max, ...)` (`BSplCLib_1.cxx:204-...`,
        // the 2D `gp_Pnt2d` overload) are not ported. No consumer in this
        // repository requests `Convert_Polynomial`.
        return Err(ConvertError::Unported);
    }

    Ok(CosAndSin { cos_numerator, sin_numerator, denominator, degree, knots, mults })
}

/// `Convert_ConicToBSplineCurve::BuildCosAndSin(Parameterisation, ...)` — the
/// periodic overload (`cxx:626-785`), which only accepts
/// `Convert_TgtThetaOver2` and `Convert_RationalC1` (`cxx:640-643`).
pub fn build_cos_and_sin_periodic(
    parameterisation: ParameterisationType,
) -> Result<CosAndSin, ConvertError> {
    let pi = std::f64::consts::PI;
    match parameterisation {
        ParameterisationType::TgtThetaOver2 => {
            // `cxx:646-670`: the full period through `TgtThetaOver2_3`, then
            // drop the last pole and make every multiplicity the degree.
            let mut cs =
                build_cos_and_sin(ParameterisationType::TgtThetaOver2_3, 0.0, 2.0 * pi)?;
            let num_poles = cs.cos_numerator.len() - 1;
            cs.cos_numerator.truncate(num_poles);
            cs.sin_numerator.truncate(num_poles);
            cs.denominator.truncate(num_poles);
            let degree = cs.degree as i32;
            for m in cs.mults.iter_mut() {
                *m = degree;
            }
            Ok(cs)
        }
        ParameterisationType::RationalC1 => {
            // `cxx:671-783`.
            let temp = build_cos_and_sin(ParameterisationType::RationalC1, 0.0, pi)?;
            let degree = 4usize;
            let num_knots = 5usize;
            let num_flat_knots = (degree - 1) * num_knots + 2 * 2;
            let num_poles = num_flat_knots - degree - 1;
            let num_periodic_poles = num_poles - 2;
            let half_pi = pi * 0.5;

            let mut flat_knots = vec![0.0f64; num_flat_knots];
            let mut index = 0usize;
            for _ in 0..2 {
                flat_knots[index] = -half_pi;
                index += 1;
            }
            for ii in 1..=num_knots {
                for _ in 0..(degree - 1) {
                    flat_knots[index] = (ii - 1) as f64 * half_pi;
                    index += 1;
                }
            }
            for _ in 0..2 {
                flat_knots[index] = 2.0 * pi + half_pi;
                index += 1;
            }

            let mut knots = vec![0.0f64; num_knots];
            let mut mults = vec![0i32; num_knots];
            for ii in 1..=num_knots {
                knots[ii - 1] = (ii - 1) as f64 * half_pi;
                mults[ii - 1] = degree as i32 - 1;
            }

            let mut cos_numerator = vec![0.0f64; num_periodic_poles];
            let mut sin_numerator = vec![0.0f64; num_periodic_poles];
            let mut denominator = vec![0.0f64; num_periodic_poles];

            let parameters = schoenberg_points(degree, &flat_knots, num_poles);
            let mut poles_array = vec![0.0f64; num_poles * 3];
            let temp_flat = knot_sequence(&temp.knots, &temp.mults, temp.degree as i32);
            // `inverse` is initialised once and flips to -1 for every parameter
            // past `pi` (`cxx:728-736`); the OCCT loop does not reset it.
            let mut inverse = 1.0f64;
            for (ii, &p) in parameters.iter().enumerate() {
                let mut param = p;
                if param > pi {
                    inverse = -1.0;
                    param -= pi;
                }
                let value1 = d0_rational_scalar(
                    &temp.cos_numerator,
                    &temp.denominator,
                    &temp_flat,
                    temp.degree,
                    param,
                )?;
                let value2 = d0_rational_scalar(
                    &temp.sin_numerator,
                    &temp.denominator,
                    &temp_flat,
                    temp.degree,
                    param,
                )?;
                let value3 = d0_scalar(&temp.denominator, &temp_flat, temp.degree, param)?;
                poles_array[ii * 3] = value1 * value3 * inverse;
                poles_array[ii * 3 + 1] = value2 * value3 * inverse;
                poles_array[ii * 3 + 2] = value3;
            }
            let contact_order_array = vec![0i32; num_poles];
            interpolate_contact(
                degree,
                &flat_knots,
                &parameters,
                &contact_order_array,
                &mut poles_array,
                3,
            )
            .map_err(|_| ConvertError::ConstructionError)?;
            for ii in 0..num_periodic_poles {
                let inverse = 1.0 / poles_array[ii * 3 + 2];
                cos_numerator[ii] = poles_array[ii * 3] * inverse;
                sin_numerator[ii] = poles_array[ii * 3 + 1] * inverse;
                denominator[ii] = poles_array[ii * 3 + 2];
            }
            Ok(CosAndSin { cos_numerator, sin_numerator, denominator, degree, knots, mults })
        }
        // `throw Standard_ConstructionError()` (`cxx:640-643`).
        _ => Err(ConvertError::ConstructionError),
    }
}
