//! `GeomLib::AdjustExtremity` and `GeomLib::ExtendCurveToPoint`.
//! Source: `ModelingData/TKGeomBase/GeomLib/GeomLib.cxx:1169-1412` (OCCT
//! 8.0.0), declarations `GeomLib.hxx`.
//!
//! Both entry points take `OccHandle<Geom_BoundedCurve>&` in OCCT and start by
//! converting it with `GeomConvert::CurveToBSplineCurve(Curve,
//! Convert_QuasiAngular)` (`cxx:1176`, `cxx:1286`), assigning the B-spline back
//! to the handle. The port takes `&mut GeomBSplineCurve` instead: every caller
//! in `ChFi3d_PerformElSpine` already holds a `Geom_BSplineCurve`, and the
//! conversion of an existing B-spline is a plain copy (`bspline_copy` in
//! `convert_bspl/curve_convert.rs`).
//!
//! `ComputeLambda` (`cxx:126-254`) is the private `GeomLib.cxx` static that
//! `ExtendCurveToPoint` runs before building the extension. OCCT hands
//! `PLib::NoDerivativeEvalPolynomial` the raw, row-major `math_Matrix` buffer
//! of `HDer` (rows `1..size-1`, columns `1..size`, `Degree = Continuity`,
//! `Dimension = size`), which evaluates each Hermite basis derivative at `t`;
//! the port builds the same coefficient rows explicitly and evaluates them with
//! [`no_derivative_eval_polynomial`](occt_core::bspl::coefs::no_derivative_eval_polynomial)
//! (`cxx:172`).

use occt_core::bspl::coefs::no_derivative_eval_polynomial;
use occt_core::convert::{ConvertError, ParameterisationType};
use occt_core::gp::{GpPnt, GpVec};
use occt_core::math_function_all_roots::FunctionAllRoots;
use occt_core::precision::CONFUSION;
use occt_math::gauss::gauss_legendre;
use occt_math::{hermite_coefficients, MathMatrix, MathVector};

use crate::bezier_curve::GeomBezierCurve;
use crate::bspline_curve::GeomBSplineCurve;
use crate::convert_bspl::CompCurveToBSplineCurve;
use crate::curve::Curve;
use crate::geom_lib_log_sample::LogSample;
use crate::geom_lib_poly_func::PolyFunc;

/// `PLib::HermiteCoefficients` failure (`cxx:1216-1219`,
/// `cxx:1298-1303`): OCCT throws `Standard_ConstructionError` and ignores the
/// `false` of the `ExtendCurveToPoint` call.
fn hermite_or_err(
    first: f64,
    last: f64,
    first_order: i32,
    last_order: i32,
) -> Result<Vec<Vec<f64>>, String> {
    hermite_coefficients(first, last, first_order, last_order)
        .ok_or_else(|| "Standard_ConstructionError in GeomLib".to_string())
}

/// `PLib::CoefficientsPoles` for the unweighted 3D case.
fn coefficients_to_poles(coeffs: &[GpVec]) -> Vec<GpPnt> {
    let pts: Vec<GpPnt> = coeffs.iter().map(|v| GpPnt::from_xyz(v.xyz())).collect();
    occt_core::bspl::plib::coefficients_poles(&pts)
}

/// Column `col` (1-based) of a `math_Matrix` as a `MathVector` over the row
/// range, i.e. `math_Matrix::Col(col)`.
fn column(matrix: &MathMatrix, col: usize) -> MathVector {
    let lower = matrix.lower_row();
    let upper = matrix.upper_row();
    let mut v = MathVector::new(lower, upper);
    for row in lower..=upper {
        v.set_value(row, matrix.value(row, col));
    }
    v
}

/// `GeomLib::AdjustExtremity(Curve, P1, P2, T1, T2)` (`cxx:1169-1266`).
///
/// Shifts the end poles of the curve so that it starts at `P1` with tangent
/// `T1` and ends at `P2` with tangent `T2`: a cubic Hermite deformation is
/// interpolated on `[0, 1]` (`cxx:1214-1229`), thickened to the curve degree
/// (`cxx:1237-1254`) and added pole by pole (`cxx:1258-1264`).
///
/// OCCT fills `FK(1..8)` and `Contact(1..4)` at `cxx:1188-1195` and never reads
/// them; they are dropped. `T1` / `T2` are normalized in place (`cxx:1200`,
/// `cxx:1208`), so an unnormalized input is accepted.
pub fn adjust_extremity(
    curve: &mut GeomBSplineCurve,
    p1: &GpPnt,
    p2: &GpPnt,
    t1: &GpVec,
    t2: &GpVec,
) -> Result<(), String> {
    // `aIn = GeomConvert::CurveToBSplineCurve(Curve, Convert_QuasiAngular);`
    // (`cxx:1176-1177`): a copy of the B-spline the caller already holds.
    let mut a_in = curve.clone();

    // `Ti(1) = Ti(2) = aIn->FirstParameter(); Ti(3) = Ti(4) =
    // aIn->LastParameter();` (`cxx:1188-1189`).
    let ti1 = a_in.first_parameter();
    let ti4 = a_in.last_parameter();

    let mut poles_def = [GpVec::zero(); 4];
    // `aIn->D1(Ti(1), P, V); PolesDef(1) = P1 - P;` (`cxx:1197-1199`).
    let (p, v) = a_in.d1(ti1);
    poles_def[0] = GpVec::new(p1.x() - p.x(), p1.y() - p.y(), p1.z() - p.z());
    // `Vtan = T1; Vtan.Normalize(); DV = Vtan * (Vtan * V) - V;`
    // `PolesDef(2) = (Ti(4) - Ti(1)) * DV;` (`cxx:1200-1202`).
    let vtan = t1.normalized();
    let dv = vtan.multiplied_scalar(vtan.dot(&v)).subtracted(&v);
    poles_def[1] = dv.multiplied_scalar(ti4 - ti1);

    // `aIn->D1(Ti(4), P, V); PolesDef(3) = P2 - P;` (`cxx:1204-1205`).
    let (p, v) = a_in.d1(ti4);
    poles_def[2] = GpVec::new(p2.x() - p.x(), p2.y() - p.y(), p2.z() - p.z());
    // `PolesDef(4) = (Ti(4) - Ti(1)) * DV;` (`cxx:1206-1211`).
    let vtan = t2.normalized();
    let dv = vtan.multiplied_scalar(vtan.dot(&v)).subtracted(&v);
    poles_def[3] = dv.multiplied_scalar(ti4 - ti1);

    // `PLib::HermiteCoefficients(0., 1., 1, 1, Mat)` (`cxx:1215-1219`).
    let mat = hermite_or_err(0.0, 1.0, 1, 1)?;

    // `Coeffs(jj) = sum_ii Mat(ii, jj) * PolesDef(ii)` (`cxx:1219-1230`).
    let mut coeffs = [GpVec::zero(); 4];
    for jj in 0..4 {
        let mut aux = GpVec::zero();
        for ii in 0..4 {
            aux = aux.added(&poles_def[ii].multiplied_scalar(mat[ii][jj]));
        }
        coeffs[jj] = aux;
    }

    // `PLib::CoefficientsPoles(Coeffs, PLib::NoWeights(), PolesDef,
    // PLib::NoWeights());` (`cxx:1229`).
    let poles = coefficients_to_poles(&coeffs);

    // `K(1) = Ti(1); K(2) = Ti(4); M.Init(4);`
    // `aDef = new Geom_BSplineCurve(PolesDef, K, M, 3);` (`cxx:1232-1238`).
    let mut a_def = GeomBSplineCurve::from_poles_knots_mults(
        poles,
        vec![ti1, ti4],
        vec![4, 4],
        3,
    )
    .map_err(|e| e.to_string())?;

    // `if (aIn->Degree() < 3) aIn->IncreaseDegree(3); else
    // aDef->IncreaseDegree(aIn->Degree());` (`cxx:1239-1245`).
    if a_in.degree() < 3 {
        a_in.increase_degree(3).map_err(|e| e.to_string())?;
    } else {
        a_def
            .increase_degree(a_in.degree())
            .map_err(|e| e.to_string())?;
    }

    // `for (ii = 2; ii < aIn->NbKnots(); ii++) aDef->InsertKnot(aIn->Knot(ii),
    // aIn->Multiplicity(ii));` (`cxx:1247-1253`).
    let a_in_knots = a_in.distinct_knots().to_vec();
    let a_in_mults = a_in.multiplicities().to_vec();
    for ii in 2..a_in_knots.len() {
        a_def
            .insert_knot(a_in_knots[ii - 1], a_in_mults[ii - 1], 0.0, true)
            .map_err(|e| e.to_string())?;
    }

    // `if (aDef->NbPoles() != aIn->NbPoles()) throw
    // Standard_ConstructionError("Inconsistent poles's number");`
    // (`cxx:1255-1258`).
    if a_def.nb_poles() != a_in.nb_poles() {
        return Err("Inconsistent poles's number".into());
    }

    // `P = aIn->Pole(ii); P.ChangeCoord() += aDef->Pole(ii).XYZ();
    // aIn->SetPole(ii, P);` (`cxx:1259-1264`).
    for ii in 0..a_def.nb_poles() {
        let def = a_def.pole(ii);
        let p = a_in.pole(ii).translated_vec(&GpVec::from_xyz(def.xyz()));
        a_in.set_pole(ii, p);
    }

    // `Curve = aIn;` (`cxx:1265`).
    *curve = a_in;
    Ok(())
}

/// `GeomLib::ExtendCurveToPoint(Curve, Point, Continuity, After)`
/// (`cxx:1269-1412`).
///
/// Appends (or prepends) a `Continuity`-order Hermite extension up to `Point`.
/// A `Continuity` outside `1..=3` returns without touching the curve
/// (`cxx:1275-1278`), and so does a target point closer than `1.e-6` to the
/// curve end (`cxx:1311-1341`).
///
/// `ComputeLambda` lowers the initial `Lambda` (the 1-argument branch of
/// `cxx:1321-1338`) by minimizing the squared deviation of the extension speed
/// (`cxx:1368`).
pub fn extend_curve_to_point(
    curve: &mut GeomBSplineCurve,
    point: &GpPnt,
    continuity: i32,
    after: bool,
) -> Result<(), String> {
    // `if (Continuity < 1 || Continuity > 3) return;` (`cxx:1275-1278`).
    if !(1..=3).contains(&continuity) {
        return Ok(());
    }
    let size = (continuity + 2) as usize;
    // `double Ubord, Tol = 1.e-6;` (`cxx:1279`).
    let tol = 1.0e-6;

    // `GeomConvert_CompCurveToBSplineCurve Concat(Curve, Convert_QuasiAngular);`
    // (`cxx:1286`).
    let mut concat = CompCurveToBSplineCurve::from_curve(curve, ParameterisationType::QuasiAngular)
        .map_err(convert_err)?;

    // `if (After) Ubord = Curve->LastParameter(); else Ubord =
    // Curve->FirstParameter();` (`cxx:1290-1296`).
    let ubord = if after {
        curve.last_parameter()
    } else {
        curve.first_parameter()
    };

    // `PLib::HermiteCoefficients(0, 1, Continuity, 0, MatCoefs);`
    // (`cxx:1298-1303`).
    let mat_coefs = hermite_or_err(0.0, 1.0, continuity, 0)?;

    // `Curve->D3(Ubord, p0, d1, d2, d3); if (!After) { d1 *= -1; d3 *= -1; }`
    // (`cxx:1305-1309`): `d2` keeps its sign.
    let (p0, mut d1, d2, mut d3) = curve.d3(ubord);
    if !after {
        d1 = d1.reversed();
        d3 = d3.reversed();
    }

    // `L1 = p0.Distance(Point); if (L1 > Tol) { ... } else return;`
    // (`cxx:1311`, `cxx:1340-1343`).
    let l1 = p0.distance(point);
    if l1 <= tol {
        return Ok(());
    }

    // `double f = Curve->FirstParameter(), t, dt, norm;`
    // `dt = (Curve->LastParameter() - f) / 9; norm = d1.Magnitude();`
    // (`cxx:1321-1324`).
    let f = curve.first_parameter();
    let dt = (curve.last_parameter() - f) / 9.0;
    let mut norm = d1.magnitude();
    let mut t = f + dt;
    for _ in 1..=8 {
        let (_, daux) = curve.d1(t);
        norm += daux.magnitude();
        t += dt;
    }
    // `norm /= 9;` (`cxx:1329`).
    norm /= 9.0;
    // `dt = d1.Magnitude() / norm;` (`cxx:1330`).
    let dt = d1.magnitude() / norm;
    let mut lambda;
    if dt < 1.5 && dt > 0.75 {
        // `Lambda = 1. / max(d1.Magnitude() / L1, Tol);` (`cxx:1332-1335`).
        lambda = 1.0 / (d1.magnitude() / l1).max(tol);
    } else {
        // `Lambda = 1. / max(norm / L1, Tol);` (`cxx:1337-1339`).
        lambda = 1.0 / (norm / l1).max(tol);
    }

    // `math_Matrix Cons(1, 3, 1, size);` (`cxx:1345-1366`).
    let mut cons = MathMatrix::with_init(1, 3, 1, size, 0.0);
    cons.set_value(1, 1, p0.x());
    cons.set_value(2, 1, p0.y());
    cons.set_value(3, 1, p0.z());
    cons.set_value(1, 2, d1.x());
    cons.set_value(2, 2, d1.y());
    cons.set_value(3, 2, d1.z());
    cons.set_value(1, size, point.x());
    cons.set_value(2, size, point.y());
    cons.set_value(3, size, point.z());
    if continuity >= 2 {
        cons.set_value(1, 3, d2.x());
        cons.set_value(2, 3, d2.y());
        cons.set_value(3, 3, d2.z());
    }
    if continuity >= 3 {
        cons.set_value(1, 4, d3.x());
        cons.set_value(2, 4, d3.y());
        cons.set_value(3, 4, d3.z());
    }

    // `ComputeLambda(Cons, MatCoefs, L1, Lambda);` (`cxx:1368`).
    compute_lambda(&cons, &mat_coefs, l1, &mut lambda);

    // `Cont(1) = p0; Cont(2) = d1 * Lambda; Cont(3) = d2 * Lambda^2;
    // Cont(4) = d3 * Lambda^3; Cont(size) = Point;` (`cxx:1371-1382`).
    let mut cont = vec![GpVec::zero(); size];
    cont[0] = GpVec::from_xyz(p0.xyz());
    cont[1] = d1.multiplied_scalar(lambda);
    if continuity >= 2 {
        cont[2] = d2.multiplied_scalar(lambda * lambda);
    }
    if continuity >= 3 {
        cont[3] = d3.multiplied_scalar(lambda * lambda * lambda);
    }
    cont[size - 1] = GpVec::from_xyz(point.xyz());

    // `ExtraCoeffs(jj).ChangeCoord() += MatCoefs(ii, jj) * Cont(ii);`
    // (`cxx:1387-1394`): `ExtraCoeffs` starts at the null point.
    let mut extra = vec![GpVec::zero(); size];
    for ii in 0..size {
        for jj in 0..size {
            extra[jj] = extra[jj].added(&cont[ii].multiplied_scalar(mat_coefs[ii][jj]));
        }
    }

    // `PLib::CoefficientsPoles(ExtraCoeffs, ...)` + `Geom_BezierCurve
    // (ExtrapPoles)` (`cxx:1397-1399`).
    let extrap_poles = coefficients_to_poles(&extra);
    let bezier = GeomBezierCurve::new(extrap_poles).map_err(|e| e.to_string())?;

    // `double dist = ExtrapPoles(1).Distance(p0); Tol += dist;`
    // (`cxx:1401-1403`).
    let dist = bezier.pole(0).distance(&p0);
    let tol = tol + dist;

    // `Ok = Concat.Add(Bezier, Tol, After); if (!Ok) throw
    // Standard_ConstructionError("ExtendCurveToPoint");`
    // (`cxx:1405-1410`). `Add` keeps the `WithRatio = true` / `MinM = 0`
    // defaults of `GeomConvert_CompCurveToBSplineCurve.hxx:56-60`.
    let ok = concat
        .add(&bezier, tol, after, true, 0)
        .map_err(convert_err)?;
    if !ok {
        return Err("ExtendCurveToPoint".into());
    }

    // `Curve = Concat.BSplineCurve();` (`cxx:1412`).
    *curve = concat
        .into_curve()
        .ok_or_else(|| "ExtendCurveToPoint".to_string())?;
    Ok(())
}

/// `ComputeLambda(Constraint, Hermit, Length, Lambda)` (`cxx:126-254`).
///
/// `Hermit` holds the `size = Continuity + 2` Hermite coefficient rows
/// (`PLib::HermiteCoefficients`); `Constraint` is the `3 x size` matrix of
/// position / first / second / third derivative / end point. The function
/// walks the Gauss points of `[-1, 1]` (`cxx:161-166`), evaluates the basis
/// derivatives into `V` (`cxx:172-173`) and accumulates the squared deviation
/// of the speed into `pol4` (`cxx:192-224`); the minimum of `pol4` over
/// `(Lambda / 1000, 50 * Lambda)` replaces `Lambda` when it lowers the error
/// (`cxx:226-253`).
fn compute_lambda(
    constraint: &MathMatrix,
    hermit: &[Vec<f64>],
    length: f64,
    lambda: &mut f64,
) {
    // `int size = Hermit.RowNumber(); int Continuity = size - 2;`
    // (`cxx:131-132`).
    let size = hermit.len();
    let continuity = size as i32 - 2;

    // `math_Matrix HDer(1, size - 1, 1, size);`
    // `HDer(ii, jj) = ii * Hermit(jj, ii + 1);` (`cxx:137-145`). OCCT reads
    // this row-major buffer as `Continuity + 1` blocks of `size`, i.e. block
    // `ii - 1` is the coefficient of `t^(ii - 1)` of every basis derivative.
    let mut hder = vec![0.0f64; (size - 1) * size];
    for jj in 1..=size {
        for ii in 1..size {
            hder[(ii - 1) * size + (jj - 1)] = ii as f64 * hermit[jj - 1][ii];
        }
    }

    // `math_Vector V(1, size);` (`cxx:146`).
    let mut v = MathVector::new(1, size);

    // `Vec2 = Constraint.Col(2); Vec2 /= Length; double squared1 =
    // Vec2.Norm2();` (`cxx:153-155`).
    let vec2_ref = column(constraint, 2).divided_scalar(length);
    let squared1 = vec2_ref.norm2();

    // `int GOrdre = 4 + 4 * Continuity;`
    // `math::GaussPoints(GOrdre, GaussP); math::GaussWeights(GOrdre, GaussW);`
    // (`cxx:161-166`).
    let g_ordre = (4 + 4 * continuity) as usize;
    let (gauss_p, gauss_w) = gauss_legendre(-1.0, 1.0, g_ordre);

    // `math_Vector pol2(1, 2 * Continuity + 1), pol4(1, 4 * Continuity + 1);`
    // (`cxx:163-164`).
    let pol2_len = (2 * continuity + 1) as usize;
    let mut pol4 = vec![0.0f64; (4 * continuity + 1) as usize];

    for ip in 0..g_ordre {
        // `t = (GaussP(ip) + 1.) / 2; GW = GaussW(ip);` (`cxx:170-171`).
        let t = (gauss_p[ip] + 1.0) / 2.0;
        let g_w = gauss_w[ip];

        // `PLib::NoDerivativeEvalPolynomial(t, Continuity, Continuity + 2,
        // DDim, polynome[0], valhder[0]); V /= Length;` (`cxx:172-173`).
        let vals = no_derivative_eval_polynomial(t, continuity, size as i32, &hder);
        for i in 0..size {
            v.set_value(i + 1, vals[i] / length);
        }

        // `Vec1 = Constraint.Col(1); Vec1 *= V(1); Vec1 += V(size) *
        // Constraint.Col(size);` (`cxx:177-179`).
        let mut vec1 = column(constraint, 1).multiplied_scalar(v.value(1));
        vec1.add(&column(constraint, size).multiplied_scalar(v.value(size)));
        // `Vec2 = Constraint.Col(2); Vec2 *= V(2);` (`cxx:180-181`).
        let vec2 = column(constraint, 2).multiplied_scalar(v.value(2));

        // `pol2` (`cxx:195-210`).
        let mut pol2 = vec![0.0f64; pol2_len];
        pol2[0] = vec1.norm2();
        pol2[1] = 2.0 * vec1.dot(&vec2);
        pol2[2] = vec2.norm2() - squared1;
        if continuity > 1 {
            let vec3 = column(constraint, 3).multiplied_scalar(v.value(3));
            pol2[2] += 2.0 * vec1.dot(&vec3);
            pol2[3] = 2.0 * vec2.dot(&vec3);
            pol2[4] = vec3.norm2();
            if continuity > 2 {
                let vec4 = column(constraint, 4).multiplied_scalar(v.value(4));
                pol2[3] += 2.0 * vec1.dot(&vec4);
                pol2[4] += 2.0 * vec2.dot(&vec4);
                pol2[5] = 2.0 * vec3.dot(&vec4);
                pol2[6] = vec4.norm2();
            }
        }
        // `for (ii = 1; ii <= pol2.Length(); ii++) { pp = ii; for (jj = 1;
        // jj < ii; jj++, pp++) pol4(pp) += 2 * GW * pol2(ii) * pol2(jj);
        // pol4(2 * ii - 1) += GW * pol2(ii)^2; }` (`cxx:214-224`).
        for ii in 1..=pol2_len {
            let mut pp = ii;
            for jj in 1..ii {
                pol4[pp - 1] += 2.0 * g_w * pol2[ii - 1] * pol2[jj - 1];
                pp += 1;
            }
            pol4[2 * ii - 2] += g_w * pol2[ii - 1] * pol2[ii - 1];
        }
    }

    // `PLib::NoDerivativeEvalPolynomial(Lambda, pol4.Length() - 1, 1,
    // pol4.Length() - 1, pol4(1), EMin);` (`cxx:228-229`).
    let degree = pol4.len() as i32 - 1;
    let mut e_min = no_derivative_eval_polynomial(*lambda, degree, 1, &pol4)[0];

    // `if (EMin > Precision::Confusion()) { ... }` (`cxx:231-253`).
    if e_min > CONFUSION {
        let mut ff = PolyFunc::new(&pol4);
        let sample = LogSample::new(*lambda / 1000.0, 50.0 * *lambda, 100);
        let solve = FunctionAllRoots::new(
            &mut ff,
            &sample,
            CONFUSION,
            CONFUSION * (length + 1.0),
            1.0e-15,
        );
        if solve.is_done() {
            for ii in 1..=solve.nb_points() {
                let t = solve.point(ii);
                let e = no_derivative_eval_polynomial(t, degree, 1, &pol4)[0];
                if e < e_min {
                    *lambda = t;
                    e_min = e;
                }
            }
        }
    }
}

fn convert_err(e: ConvertError) -> String {
    format!("{e:?}")
}
