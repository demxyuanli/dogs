//! `AppParCurves_LeastSquare` — the parallel-curves least-squares resolver
//! that `AppParCurves_ResolConstraint`, `AppDef_BSplineCompute` and
//! `GeomFill_AppSurf` drive.
//!
//! Source: `ModelingData/TKGeomBase/AppParCurves/AppParCurves_LeastSquare.gxx:1-1860`
//! and the generated declaration `AppDef_TheLeastSquares.hxx:44-272`
//! (OCCT 8.0.0; the template is instantiated as `AppDef_TheLeastSquares` for
//! `AppDef_MultiLine` / `AppDef_MyLineTool`, which is the single instantiation
//! this port keeps — the generics collapse to those two concrete types).
//!
//! Layout:
//! - [`mod@init`] — `Init` (constraint bookkeeping, pole seeds, sizing),
//! - [`mod@affect`] — `Affect` + the two `CheckTangents` helpers,
//! - [`mod@perform`] — the four `Perform` entry points,
//! - [`mod@inner`] — `SearchIndex`, the three `MakeTAA` overloads and
//!   `ComputeFunction`.
//!
//! Translation notes:
//! - `mypoles`, `A`, `DA`, `B2`, `mypoints` and `theError` are `math_Matrix`
//!   values with OCCT's own bounds (rows from `FirstPoint`, columns from 1);
//!   the port keeps those bounds and the 1-based element access.
//! - `myknots` / `mymults` are null handles in the Bezier constructors and set
//!   together in the B-spline ones; `Option` is that state, and every place
//!   OCCT dereferences a null handle panics, OCCT's crash counterpart.
//! - `lambda1` / `lambda2` are uninitialised `double`s in OCCT; the port seeds
//!   them with `0.0` (they are only read after `Perform` has written them).
//! - OCCT raises `StdFail_NotDone` from `Error`, `ErrorGradient`,
//!   `BSplineValue` and `Standard_NoSuchObject` from `BezierValue`; the port
//!   panics with the same names.
//! - `FirstConstraint >= AppParCurves_TangencyPoint` is the C++ enum order
//!   test (`NoConstraint < PassPoint < TangencyPoint < CurvaturePoint`),
//!   spelled [`at_least_tangency`].

use occt_core::gp::{GpPnt, GpPnt2d};
use occt_math::{MathIntVector, MathMatrix, MathVector};

use crate::approx::line_tool::LineTool;
use crate::app_par_curves::bernstein::AppParCurves;
use crate::app_par_curves::constraint::Constraint;
use crate::app_par_curves::multi_bspline_curve::MultiBSpCurve;
use crate::app_par_curves::multi_curve::MultiCurve;
use crate::app_par_curves::multi_point::MultiPoint;

mod affect;
mod init;
mod inner;
mod perform;

/// C++ `Cons >= AppParCurves_TangencyPoint` (`Constraint` declaration order).
fn at_least_tangency(c: Constraint) -> bool {
    matches!(c, Constraint::TangencyPoint | Constraint::CurvaturePoint)
}

/// C++ `Cons >= AppParCurves_PassPoint`.
fn at_least_pass(c: Constraint) -> bool {
    matches!(
        c,
        Constraint::PassPoint | Constraint::TangencyPoint | Constraint::CurvaturePoint
    )
}

/// `AppParCurves_LeastSquare` (`AppDef_TheLeastSquares.hxx:44-272`).
#[derive(Debug, Clone)]
pub struct LeastSquare {
    /// `FirstConstraint` (`hxx:237`).
    first_constraint: Constraint,
    /// `LastConstraint` (`hxx:238`).
    last_constraint: Constraint,
    /// `SCU` (`hxx:239`): the result curve.
    scu: MultiBSpCurve,
    /// `myknots` (`hxx:240`): the distinct knots, null for the Bezier path.
    myknots: Option<Vec<f64>>,
    /// `mymults` (`hxx:241`): the knot multiplicities.
    mymults: Option<Vec<i32>>,
    /// `mypoles` (`hxx:242`): the pole grid, one column per coordinate.
    mypoles: MathMatrix,
    /// `A` (`hxx:243`): the basis matrix.
    a: MathMatrix,
    /// `DA` (`hxx:244`): the basis derivative matrix.
    da: MathMatrix,
    /// `B2` (`hxx:245`): the residual right-hand side.
    b2: MathMatrix,
    /// `mypoints` (`hxx:246`): the flattened input points.
    mypoints: MathMatrix,
    /// `Vflatknots` (`hxx:247`): the flat knot vector.
    vflatknots: MathVector,
    /// `Vec1t` (`hxx:248`).
    vec1t: MathVector,
    /// `Vec1c` (`hxx:249`).
    vec1c: MathVector,
    /// `Vec2t` (`hxx:250`).
    vec2t: MathVector,
    /// `Vec2c` (`hxx:251`).
    vec2c: MathVector,
    /// `theError` (`hxx:252`).
    the_error: MathMatrix,
    /// `myindex` (`hxx:253`): index of the first non-null basis column.
    myindex: MathIntVector,
    /// `lambda1` (`hxx:254`).
    lambda1: f64,
    /// `lambda2` (`hxx:255`).
    lambda2: f64,
    /// `FirstP` (`hxx:256`): first fitted row.
    first_p: i32,
    /// `LastP` (`hxx:257`): last fitted row.
    last_p: i32,
    /// `Nlignes` (`hxx:258`).
    nlignes: i32,
    /// `Ninc` (`hxx:259`).
    ninc: i32,
    /// `NA` (`hxx:260`): the flattened point dimension.
    na: i32,
    /// `myfirstp` (`hxx:261`).
    myfirstp: i32,
    /// `mylastp` (`hxx:262`).
    mylastp: i32,
    /// `resinit` (`hxx:263`).
    resinit: i32,
    /// `resfin` (`hxx:264`).
    resfin: i32,
    /// `nbP2d` (`hxx:265`).
    nbp2d: i32,
    /// `nbP` (`hxx:266`).
    nbp: i32,
    /// `nbpoles` (`hxx:267`).
    nbpoles: i32,
    /// `deg` (`hxx:268`).
    deg: i32,
    /// `done` (`hxx:269`).
    done: bool,
    /// `iscalculated` (`hxx:270`).
    iscalculated: bool,
    /// `isready` (`hxx:271`).
    isready: bool,
}

impl LeastSquare {
    /// `AppParCurves_LeastSquare(SSP, FirstPoint, LastPoint, FirstCons,
    /// LastCons, Parameters, NbPol)` (`gxx:137-167`): the Bezier constructor,
    /// which immediately runs `Perform(Parameters)`.
    pub fn with_parameters<L: LineTool>(
        ml: &L,
        first_point: i32,
        last_point: i32,
        first_cons: Constraint,
        last_cons: Constraint,
        parameters: &MathVector,
        nb_pol: i32,
    ) -> Self {
        let mut this = Self::new(ml, first_point, last_point, first_cons, last_cons, nb_pol);
        this.init(ml, first_point, last_point);
        this.perform_parameters(parameters);
        this
    }

    /// `AppParCurves_LeastSquare(SSP, FirstPoint, LastPoint, FirstCons,
    /// LastCons, NbPol)` (`gxx:169-196`): the Bezier constructor that leaves
    /// the resolution to the caller's `Perform`.
    pub fn new<L: LineTool>(
        ml: &L,
        first_point: i32,
        last_point: i32,
        first_cons: Constraint,
        last_cons: Constraint,
        nb_pol: i32,
    ) -> Self {
        let nb_bcols = Self::nb_b_columns(ml) as usize;
        let b2_lower = Self::the_first_point(first_cons, first_point);
        let b2_upper = b2_lower.max(Self::the_last_point(last_cons, last_point));
        let mut this = Self {
            first_constraint: first_cons,
            last_constraint: last_cons,
            scu: MultiBSpCurve::with_nb_poles(nb_pol.max(0) as usize),
            myknots: None,
            mymults: None,
            mypoles: MathMatrix::new(1, nb_pol.max(0) as usize, 1, nb_bcols),
            a: MathMatrix::new(
                first_point.max(0) as usize,
                last_point.max(0) as usize,
                1,
                nb_pol.max(0) as usize,
            ),
            da: MathMatrix::new(
                first_point.max(0) as usize,
                last_point.max(0) as usize,
                1,
                nb_pol.max(0) as usize,
            ),
            b2: MathMatrix::new(b2_lower.max(0) as usize, b2_upper.max(0) as usize, 1, nb_bcols),
            mypoints: MathMatrix::new(
                first_point.max(0) as usize,
                last_point.max(0) as usize,
                1,
                nb_bcols,
            ),
            vflatknots: MathVector::new(1, 1),
            vec1t: MathVector::new(1, nb_bcols),
            vec1c: MathVector::new(1, nb_bcols),
            vec2t: MathVector::new(1, nb_bcols),
            vec2c: MathVector::new(1, nb_bcols),
            the_error: MathMatrix::new(
                first_point.max(0) as usize,
                last_point.max(0) as usize,
                1,
                (ml.nb_p3d() + ml.nb_p2d()).max(0) as usize,
            ),
            myindex: MathIntVector::new(first_point.max(0) as usize, last_point.max(0) as usize),
            lambda1: 0.0,
            lambda2: 0.0,
            first_p: 0,
            last_p: 0,
            nlignes: 0,
            ninc: 0,
            na: 0,
            myfirstp: 0,
            mylastp: 0,
            resinit: 0,
            resfin: 0,
            nbp2d: 0,
            nbp: 0,
            nbpoles: 0,
            deg: 0,
            done: false,
            iscalculated: false,
            isready: false,
        };
        this.first_constraint = first_cons;
        this.last_constraint = last_cons;
        this.init(ml, first_point, last_point);
        this
    }

    /// `AppParCurves_LeastSquare(SSP, Knots, Mults, FirstPoint, LastPoint,
    /// FirstCons, LastCons, Parameters, NbPol)` (`gxx:197-234`).
    #[allow(clippy::too_many_arguments)]
    pub fn with_knots_parameters<L: LineTool>(
        ml: &L,
        knots: &[f64],
        mults: &[i32],
        first_point: i32,
        last_point: i32,
        first_cons: Constraint,
        last_cons: Constraint,
        parameters: &MathVector,
        nb_pol: i32,
    ) -> Self {
        let mut this = Self::with_knots(
            ml,
            knots,
            mults,
            first_point,
            last_point,
            first_cons,
            last_cons,
            nb_pol,
        );
        this.perform_parameters(parameters);
        this
    }

    /// `AppParCurves_LeastSquare(SSP, Knots, Mults, FirstPoint, LastPoint,
    /// FirstCons, LastCons, NbPol)` (`gxx:236-274`).
    #[allow(clippy::too_many_arguments)]
    pub fn with_knots<L: LineTool>(
        ml: &L,
        knots: &[f64],
        mults: &[i32],
        first_point: i32,
        last_point: i32,
        first_cons: Constraint,
        last_cons: Constraint,
        nb_pol: i32,
    ) -> Self {
        let nb_bcols = Self::nb_b_columns(ml) as usize;
        let flat_len: i32 = mults.iter().sum();
        let b2_lower = Self::the_first_point(first_cons, first_point);
        let b2_upper = b2_lower.max(Self::the_last_point(last_cons, last_point));
        let mut this = Self {
            first_constraint: first_cons,
            last_constraint: last_cons,
            scu: MultiBSpCurve::with_nb_poles(nb_pol.max(0) as usize),
            myknots: Some(knots.to_vec()),
            mymults: Some(mults.to_vec()),
            mypoles: MathMatrix::new(1, nb_pol.max(0) as usize, 1, nb_bcols),
            a: MathMatrix::new(
                first_point.max(0) as usize,
                last_point.max(0) as usize,
                1,
                nb_pol.max(0) as usize,
            ),
            da: MathMatrix::new(
                first_point.max(0) as usize,
                last_point.max(0) as usize,
                1,
                nb_pol.max(0) as usize,
            ),
            b2: MathMatrix::new(b2_lower.max(0) as usize, b2_upper.max(0) as usize, 1, nb_bcols),
            mypoints: MathMatrix::new(
                first_point.max(0) as usize,
                last_point.max(0) as usize,
                1,
                nb_bcols,
            ),
            vflatknots: MathVector::new(1, flat_len.max(0) as usize),
            vec1t: MathVector::new(1, nb_bcols),
            vec1c: MathVector::new(1, nb_bcols),
            vec2t: MathVector::new(1, nb_bcols),
            vec2c: MathVector::new(1, nb_bcols),
            the_error: MathMatrix::new(
                first_point.max(0) as usize,
                last_point.max(0) as usize,
                1,
                (ml.nb_p3d() + ml.nb_p2d()).max(0) as usize,
            ),
            myindex: MathIntVector::new(first_point.max(0) as usize, last_point.max(0) as usize),
            lambda1: 0.0,
            lambda2: 0.0,
            first_p: 0,
            last_p: 0,
            nlignes: 0,
            ninc: 0,
            na: 0,
            myfirstp: 0,
            mylastp: 0,
            resinit: 0,
            resfin: 0,
            nbp2d: 0,
            nbp: 0,
            nbpoles: 0,
            deg: 0,
            done: false,
            iscalculated: false,
            isready: false,
        };
        this.scu.set_knots(knots);
        this.scu.set_multiplicities(mults);
        this.first_constraint = first_cons;
        this.last_constraint = last_cons;
        this.init(ml, first_point, last_point);
        this
    }

    /// `NbBColumns(SSP)` (`gxx:1205-1213`): `3 * NbP3d + 2 * NbP2d`.
    fn nb_b_columns<L: LineTool>(ml: &L) -> i32 {
        ml.nb_p3d() * 3 + ml.nb_p2d() * 2
    }

    /// `TheFirstPoint(FirstCons, FirstPoint)` (`gxx:1465-1472`).
    fn the_first_point(first_cons: Constraint, first_point: i32) -> i32 {
        if first_cons == Constraint::NoConstraint {
            first_point
        } else {
            first_point + 1
        }
    }

    /// `TheLastPoint(LastCons, LastPoint)` (`gxx:1474-1482`).
    fn the_last_point(last_cons: Constraint, last_point: i32) -> i32 {
        if last_cons == Constraint::NoConstraint {
            last_point
        } else {
            last_point - 1
        }
    }

    /// `IsDone()` (`gxx:1395-1399`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `FirstLambda()` (`gxx:1385-1388`).
    pub fn first_lambda(&self) -> f64 {
        self.lambda1
    }

    /// `LastLambda()` (`gxx:1390-1393`).
    pub fn last_lambda(&self) -> f64 {
        self.lambda2
    }

    /// `FunctionMatrix()` (`gxx:1445-1453`).
    pub fn function_matrix(&self) -> &MathMatrix {
        if !self.done {
            panic!("StdFail_NotDone: AppParCurves_LeastSquare::FunctionMatrix");
        }
        &self.a
    }

    /// `DerivativeFunctionMatrix()` (`gxx:1455-1463`).
    pub fn derivative_function_matrix(&self) -> &MathMatrix {
        if !self.done {
            panic!("StdFail_NotDone: AppParCurves_LeastSquare::DerivativeFunctionMatrix");
        }
        &self.da
    }

    /// `Points()` (`gxx:1496-1499`).
    pub fn points(&self) -> &MathMatrix {
        &self.mypoints
    }

    /// `Poles()` (`gxx:1501-1504`).
    pub fn poles(&self) -> &MathMatrix {
        &self.mypoles
    }

    /// `KIndex()` (`gxx:1506-1509`).
    pub fn k_index(&self) -> &MathIntVector {
        &self.myindex
    }

    /// `ComputeFunction(Parameters)` (`gxx:1484-1494`).
    fn compute_function(&mut self, parameters: &MathVector) {
        if self.myknots.is_none() {
            AppParCurves::bernstein(
                self.nbpoles.max(0) as usize,
                parameters,
                &mut self.a,
                &mut self.da,
            );
        } else {
            let deg = self.deg;
            let flat = self.vflatknots.clone();
            AppParCurves::spline_function(
                self.nbpoles,
                deg,
                parameters,
                &flat,
                &mut self.a,
                &mut self.da,
                &mut self.myindex,
            );
        }
    }

    /// `BSplineValue()` (`gxx:1408-1445`): fills `SCU` from `mypoles` over
    /// `ideb..ifin` and returns it.
    pub fn bspline_value(&mut self) -> &MultiBSpCurve {
        if !self.done {
            panic!("StdFail_NotDone: AppParCurves_LeastSquare::BSplineValue");
        }
        let npoints = self.nbp + self.nbp2d;
        let mut ideb = self.resinit;
        let mut ifin = self.resfin;
        if ideb >= 2 {
            ideb = 2;
        }
        if ifin <= self.nbpoles - 1 {
            ifin = self.nbpoles - 1;
        }

        for i in ideb..=ifin {
            let mut j2 = 1;
            let mut mpole = MultiPoint::with_size(self.nbp.max(0) as usize, self.nbp2d.max(0) as usize);
            for j in 1..=self.nbp {
                let pt = GpPnt::new(
                    self.mypoles.value(i as usize, j2 as usize),
                    self.mypoles.value(i as usize, (j2 + 1) as usize),
                    self.mypoles.value(i as usize, (j2 + 2) as usize),
                );
                mpole.set_point(j, &pt);
                j2 += 3;
            }
            for j in self.nbp + 1..=npoints {
                let p2d = GpPnt2d::new(
                    self.mypoles.value(i as usize, j2 as usize),
                    self.mypoles.value(i as usize, (j2 + 1) as usize),
                );
                mpole.set_point2d(j, &p2d);
                j2 += 2;
            }
            self.scu.base_mut().set_value(i, &mpole);
        }
        &self.scu
    }

    /// `BezierValue()` (`gxx:1402-1406`): raises `Standard_NoSuchObject` for
    /// the B-spline constructors.
    pub fn bezier_value(&mut self) -> MultiCurve {
        if self.myknots.is_some() {
            panic!("Standard_NoSuchObject: AppParCurves_LeastSquare::BezierValue");
        }
        self.bspline_value().base().clone()
    }

    /// `Distance()` (`gxx:1371-1386`): the square roots of `theError`, cached
    /// behind `iscalculated`.
    pub fn distance(&mut self) -> &MathMatrix {
        if !self.iscalculated {
            for i in self.myfirstp..=self.mylastp {
                for j in 1..=self.nbp + self.nbp2d {
                    let v = self.the_error.value(i as usize, j as usize).sqrt();
                    self.the_error.set_value(i as usize, j as usize, v);
                }
            }
            self.iscalculated = true;
        }
        &self.the_error
    }

    /// `Error(F, MaxE3d, MaxE2d)` (`gxx:1216-1279`): sum of squares, plus the
    /// largest 3d and 2d distances.
    pub fn error(&mut self) -> (f64, f64, f64) {
        if !self.done {
            panic!("StdFail_NotDone: AppParCurves_LeastSquare::Error");
        }
        let mut max_e3d = 0.0;
        let mut max_e2d = 0.0;
        let mut f = 0.0;
        let mut i2 = 1;
        let nbpoles = self.nbpoles as usize;
        let mut px = MathVector::new(1, nbpoles);
        let mut py = MathVector::new(1, nbpoles);
        let mut pz = MathVector::new(1, nbpoles);

        let nb_p = self.nbp;
        let nb_p2d = self.nbp2d;
        for k in 1..=nb_p + nb_p2d {
            let i21 = i2 + 1;
            let i22 = i2 + 2;
            for i in 1..=self.nbpoles {
                px.set_value(i as usize, self.mypoles.value(i as usize, i2 as usize));
                py.set_value(i as usize, self.mypoles.value(i as usize, i21 as usize));
                if k <= nb_p {
                    pz.set_value(i as usize, self.mypoles.value(i as usize, i22 as usize));
                }
            }
            for i in self.first_p..=self.last_p {
                let mut aa = 0.0;
                let mut bb = 0.0;
                let mut cc = 0.0;
                let indexdeb = self.myindex.value(i as usize) + 1;
                let indexfin = indexdeb + self.deg;
                for j in indexdeb..=indexfin {
                    let aij = self.a.value(i as usize, j as usize);
                    aa += aij * px.value(j as usize);
                    bb += aij * py.value(j as usize);
                    if k <= nb_p {
                        cc += aij * pz.value(j as usize);
                    }
                }
                let fx = aa - self.mypoints.value(i as usize, i2 as usize);
                let fy = bb - self.mypoints.value(i as usize, i21 as usize);
                let mut fi = fx * fx + fy * fy;
                if k <= nb_p {
                    let fz = cc - self.mypoints.value(i as usize, i22 as usize);
                    fi += fz * fz;
                    if fi > max_e3d {
                        max_e3d = fi;
                    }
                } else if fi > max_e2d {
                    max_e2d = fi;
                }
                self.the_error.set_value(i as usize, k as usize, fi);
                f += fi;
            }
            if k <= nb_p {
                i2 += 3;
            } else {
                i2 += 2;
            }
        }
        (f, max_e3d.sqrt(), max_e2d.sqrt())
    }

    /// `ErrorGradient(Grad, F, MaxE3d, MaxE2d)` (`gxx:1282-1368`).
    pub fn error_gradient(&mut self, grad: &mut MathVector) -> (f64, f64, f64) {
        if !self.done {
            panic!("StdFail_NotDone: AppParCurves_LeastSquare::ErrorGradient");
        }
        let mut max_e3d = 0.0;
        let mut max_e2d = 0.0;
        let mut f = 0.0;
        let mut i2 = 1;
        let nbpoles = self.nbpoles as usize;
        let mut px = MathVector::new(1, nbpoles);
        let mut py = MathVector::new(1, nbpoles);
        let mut pz = MathVector::new(1, nbpoles);

        for k in grad.lower()..=grad.upper() {
            grad.set_value(k, 0.0);
        }

        let nb_p = self.nbp;
        let nb_p2d = self.nbp2d;
        for k in 1..=nb_p + nb_p2d {
            let i21 = i2 + 1;
            let i22 = i2 + 2;
            for i in 1..=self.nbpoles {
                px.set_value(i as usize, self.mypoles.value(i as usize, i2 as usize));
                py.set_value(i as usize, self.mypoles.value(i as usize, i21 as usize));
                if k <= nb_p {
                    pz.set_value(i as usize, self.mypoles.value(i as usize, i22 as usize));
                }
            }
            for i in self.first_p..=self.last_p {
                let mut aa = 0.0;
                let mut bb = 0.0;
                let mut cc = 0.0;
                let mut daa = 0.0;
                let mut dbb = 0.0;
                let mut dcc = 0.0;
                let indexdeb = self.myindex.value(i as usize) + 1;
                let indexfin = indexdeb + self.deg;
                for j in indexdeb..=indexfin {
                    let aij = self.a.value(i as usize, j as usize);
                    let daij = self.da.value(i as usize, j as usize);
                    aa += aij * px.value(j as usize);
                    daa += daij * px.value(j as usize);
                    bb += aij * py.value(j as usize);
                    dbb += daij * py.value(j as usize);
                    if k <= nb_p {
                        cc += aij * pz.value(j as usize);
                        dcc += daij * pz.value(j as usize);
                    }
                }
                let fx = aa - self.mypoints.value(i as usize, i2 as usize);
                let fy = bb - self.mypoints.value(i as usize, (i2 + 1) as usize);
                let mut fi = fx * fx + fy * fy;
                let mut gr = 2.0 * (daa * fx + dbb * fy);

                if k <= nb_p {
                    let fz = cc - self.mypoints.value(i as usize, (i2 + 2) as usize);
                    fi += fz * fz;
                    gr += 2.0 * dcc * fz;
                    if fi > max_e3d {
                        max_e3d = fi;
                    }
                } else if fi > max_e2d {
                    max_e2d = fi;
                }
                self.the_error.set_value(i as usize, k as usize, fi);
                let g = grad.value(i as usize) + gr;
                grad.set_value(i as usize, g);
                f += fi;
            }
            if k <= nb_p {
                i2 += 3;
            } else {
                i2 += 2;
            }
        }
        (f, max_e3d.sqrt(), max_e2d.sqrt())
    }
}

/// `const Vec1t`-style 1-based read used by the `Perform` bodies when they
/// copy caller-supplied tangent vectors into the member arrays (`gxx:752-758`).
pub(crate) fn copy_shifted(dst: &mut MathVector, src: &MathVector) {
    let lower_src = src.lower();
    for i in 1..=dst.upper() {
        dst.set_value(i, src.value(i + lower_src - 1));
    }
}



