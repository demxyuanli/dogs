//! `AppParCurves_ResolConstraint` — the constrained resolver that turns the
//! least-squares pole guess into a curve passing through / tangent to the
//! constrained points, by the Uzawa algorithm of `math_Uzawa`.
//!
//! Source: `ModelingData/TKGeomBase/AppParCurves/AppParCurves_ResolConstraint.gxx:40-946`
//! and the generated declaration `AppDef_TheResol.hxx:43-121` (OCCT 8.0.0).
//! The `.gxx` is a template instantiated five times in `AppDef`
//! (`AppDef_TheResol`, `AppDef_ResConstraintOfTheGradient`,
//! `AppDef_ResConstraintOfMyGradientOfCompute`,
//! `AppDef_ResConstraintOfMyGradientbisOfBSplineCompute`,
//! `AppDef_ResConstraintOfTheGradient`); they differ only in the substituted
//! `MultiLine` / `LineTool` / `MultiPoint` names, so the port keeps the single
//! `AppDef` instantiation and collapses the generics onto
//! [`crate::app_def::MultiLine`] / [`crate::app_def::MyLineTool`].
//!
//! Layout:
//! - [`mod@perform`] — the constructor body (`gxx:40-502`): constraint scan,
//!   `Cont` / `Secont` / `Start` filling, the Uzawa call and the pole
//!   write-back,
//! - [`mod@ratios`] — the "same tangency ratio on every curve" rows
//!   (`gxx:288-406`),
//! - [`mod@derivative`] — `ConstraintDerivative` (`gxx:573-946`).
//!
//! Translation notes:
//! - `Cont`, `DeCont`, `Secont`, `CTCinv` and `Vardua` are sized from
//!   `NbConstraints` / `NbColumns` (`gxx:49-68`); the port keeps those sizes
//!   and the 1-based element access.
//! - `Done` is **not** in the C++ constructor init list (`hxx:97` declares it
//!   before `Cont`), so OCCT reads an indeterminate value on the paths that
//!   fail early; the port seeds it with `false`, the only value any caller
//!   can observe before `gxx:471` sets it.
//! - `Err` (`hxx:98`) is never written by the `.gxx`, and the declared
//!   `Error()` (`hxx:72`) has no definition, so neither is ported.
//! - The curvature rows of `gxx:410-452` are commented out in OCCT itself;
//!   they are not ported (there is no live control flow to translate).
//! - `math_Matrix(1, 0, ...)` (empty range) is legal in OCCT but not in
//!   [`MathMatrix`], which needs `upper >= lower`; the empty sizes that only
//!   occur when `IncPass == 0` (`gxx:114-118` returns immediately) become a
//!   1x1 zero matrix that no accessor reads.
//! - `math_Uzawa(Cont, Secont, Start, Tolerance)` (`gxx:455`) is the
//!   equality-only constructor, so `EpsLic` and `NbIterations` keep their
//!   header defaults `1.0e-06` / `500` (`math_Uzawa.hxx:49-58`).

use occt_math::{MathMatrix, MathVector};

use crate::approx::line_tool::LineTool;
use crate::app_par_curves::constraint::{Constraint, ConstraintCouple};
use crate::app_par_curves::multi_curve::MultiCurve;

mod derivative;
mod perform;
mod ratios;

/// C++ `Cons >= AppParCurves_PassPoint` (`Constraint` declaration order).
pub(super) fn at_least_pass(c: Constraint) -> bool {
    matches!(
        c,
        Constraint::PassPoint | Constraint::TangencyPoint | Constraint::CurvaturePoint
    )
}

/// C++ `Cons >= AppParCurves_TangencyPoint` (`Constraint` declaration order).
pub(super) fn at_least_tangency(c: Constraint) -> bool {
    matches!(c, Constraint::TangencyPoint | Constraint::CurvaturePoint)
}

/// `math_Matrix(1, rows, 1, cols)` with OCCT's zero initialisation; an empty
/// request keeps a 1x1 stand-in (see the module notes).
pub(super) fn zero_matrix(rows: usize, cols: usize) -> MathMatrix {
    MathMatrix::with_init(1, rows.max(1), 1, cols.max(1), 0.0)
}

/// `math_Matrix::Set(I1, I2, J1, J2, M)` (`math_Matrix.lxx`): copy `src` from
/// its own lower bounds into the row/column window.
pub(super) fn set_block(
    dst: &mut MathMatrix,
    i1: i32,
    i2: i32,
    j1: i32,
    j2: i32,
    src: &MathMatrix,
) {
    if i1 > i2 || j1 > j2 || (i2 - i1 + 1) as usize != src.row_count() || (j2 - j1 + 1) as usize != src.col_count() {
        panic!("Standard_DimensionError: math_Matrix::Set");
    }
    let lr = src.lower_row();
    let lc = src.lower_col();
    for i in i1..=i2 {
        for j in j1..=j2 {
            let v = src.value(
                (i - i1) as usize + lr,
                (j - j1) as usize + lc,
            );
            dst.set_value(i as usize, j as usize, v);
        }
    }
}

/// `Ibont(K, I)` of the `NCollection_Array2<int>(1, NbCu, 1, IncTan)`
/// (`gxx:131`), stored row-major as a flat `i32` buffer.
pub(super) fn ibont(ibonts: &[i32], inc_tan: usize, k: i32, i: i32) -> i32 {
    ibonts[((k - 1) as usize) * inc_tan + (i - 1) as usize]
}

/// `Ibont(K, I) = value` (`gxx:172-182`).
pub(super) fn set_ibont(ibonts: &mut [i32], inc_tan: usize, k: i32, i: i32, val: i32) {
    ibonts[((k - 1) as usize) * inc_tan + (i - 1) as usize] = val;
}

/// The constraint bookkeeping of the constructor prologue (`gxx:87-111`).
pub(super) struct Counts {
    /// `IncPass` (`hxx:104`).
    pub inc_pass: usize,
    /// `IncTan` (`hxx:105`).
    pub inc_tan: usize,
    /// `IncCurv` (`hxx:106`).
    pub inc_curv: usize,
    /// `IPas` (`hxx:107`).
    pub ipas: Vec<i32>,
    /// `ITan` (`hxx:108`).
    pub itan: Vec<i32>,
    /// `ICurv` (`hxx:109`).
    pub icurv: Vec<i32>,
    /// `FC` (`gxx:80`): the constraint on `FirstPoint`.
    pub fc: Constraint,
    /// `LC` (`gxx:80`): the constraint on `LastPoint`.
    pub lc: Constraint,
}

impl Counts {
    /// `gxx:82-111`: one pass over the constraint couples.
    pub fn scan(
        constraints: &[ConstraintCouple],
        first_point: i32,
        last_point: i32,
        nb_points: usize,
    ) -> Self {
        let mut this = Self {
            inc_pass: 0,
            inc_tan: 0,
            inc_curv: 0,
            ipas: vec![0; nb_points],
            itan: vec![0; nb_points],
            icurv: vec![0; nb_points],
            fc: Constraint::NoConstraint,
            lc: Constraint::NoConstraint,
        };
        for couple in constraints {
            let myindex = couple.index();
            let cons = couple.constraint();
            if myindex == first_point {
                this.fc = cons;
            }
            if myindex == last_point {
                this.lc = cons;
            }
            if at_least_pass(cons) {
                this.inc_pass += 1;
                this.ipas[this.inc_pass - 1] = myindex;
            }
            if at_least_tangency(cons) {
                this.inc_tan += 1;
                this.itan[this.inc_tan - 1] = myindex;
            }
            if cons == Constraint::CurvaturePoint {
                this.inc_curv += 1;
                this.icurv[this.inc_curv - 1] = myindex;
            }
        }
        this
    }
}

/// `AppParCurves_ResolConstraint` (`AppDef_TheResol.hxx:43-121`).
#[derive(Debug, Clone)]
pub struct ResolConstraint {
    /// `Done` (`hxx:97`).
    done: bool,
    /// `Cont` (`hxx:99`): the constraint matrix.
    cont: MathMatrix,
    /// `DeCont` (`hxx:100`): the derivative of `Cont`.
    de_cont: MathMatrix,
    /// `Secont` (`hxx:101`): the right-hand side.
    secont: MathVector,
    /// `CTCinv` (`hxx:102`): `Inv(Cont * T(Cont))`.
    ctc_inv: MathMatrix,
    /// `Vardua` (`hxx:103`): the dual variables.
    vardua: MathVector,
    /// `IncPass` (`hxx:104`).
    inc_pass: usize,
    /// `IncTan` (`hxx:105`).
    inc_tan: usize,
    /// `IncCurv` (`hxx:106`): only the commented-out curvature rows
    /// (`gxx:410-452`) would read it, but it is part of the class state.
    #[allow(dead_code)]
    inc_curv: usize,
    /// `IPas` (`hxx:107`): the indices of the passing points.
    ipas: Vec<i32>,
    /// `ITan` (`hxx:108`): the indices of the tangency points.
    itan: Vec<i32>,
    /// `ICurv` (`hxx:109`): only the commented-out curvature rows would read
    /// it, but it is part of the class state.
    #[allow(dead_code)]
    icurv: Vec<i32>,
}

impl ResolConstraint {
    /// `AppDef_TheResol(SSP, SCurv, FirstPoint, LastPoint, Constraints, Bern,
    /// DerivativeBern, Tolerance = 1.0e-10)` (`gxx:40-503`).
    ///
    /// `bern` / `derivative_bern` are the Bernstein matrices read at the
    /// constrained point indices (`gxx:140`, `gxx:187`): rows are point
    /// indices, columns `1..=NbPoles`.
    pub fn new<L: LineTool>(
        ssp: &L,
        scurv: &mut MultiCurve,
        first_point: i32,
        last_point: i32,
        constraints: &[ConstraintCouple],
        bern: &MathMatrix,
        derivative_bern: &MathMatrix,
        tolerance: f64,
    ) -> Self {
        let nb_points = (last_point - first_point + 1).max(0) as usize;
        let counts = Counts::scan(constraints, first_point, last_point, nb_points);
        let nb_cons = Self::nb_constraints(ssp, constraints);
        let nb_cols = Self::nb_columns(ssp, scurv.nb_poles() as i32 - 1);

        let mut this = Self {
            done: false,
            cont: zero_matrix(nb_cons, nb_cols),
            de_cont: zero_matrix(nb_cons, nb_cols),
            secont: if nb_cons == 0 {
                MathVector::from_slice(&[])
            } else {
                MathVector::with_init(1, nb_cons, 0.0)
            },
            ctc_inv: zero_matrix(nb_cons, nb_cons),
            vardua: if nb_cons == 0 {
                MathVector::from_slice(&[])
            } else {
                MathVector::new(1, nb_cons)
            },
            inc_pass: counts.inc_pass,
            inc_tan: counts.inc_tan,
            inc_curv: counts.inc_curv,
            ipas: counts.ipas,
            itan: counts.itan,
            icurv: counts.icurv,
        };

        // `if (IncPass == 0) { Done = true; return; }` (`gxx:113-117`).
        if this.inc_pass == 0 {
            this.done = true;
            return this;
        }

        perform::perform(
            &mut this,
            ssp,
            scurv,
            bern,
            derivative_bern,
            tolerance,
            counts.fc,
            counts.lc,
        );
        this
    }

    /// `IsDone()` (`gxx:505-508`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `ConstraintMatrix()` (`gxx:558-561`).
    pub fn constraint_matrix(&self) -> &MathMatrix {
        &self.cont
    }

    /// `InverseMatrix()` (`gxx:563-566`): `Inv(Cont * T(Cont))`.
    pub fn inverse_matrix(&self) -> &MathMatrix {
        &self.ctc_inv
    }

    /// `Duale()` (`gxx:568-571`): the dual variables of the system.
    pub fn duale(&self) -> &MathVector {
        &self.vardua
    }

    /// `ConstraintDerivative(SSP, Parameters, Deg, DA)` (`gxx:573-946`).
    pub fn constraint_derivative<L: LineTool>(
        &mut self,
        ssp: &L,
        parameters: &MathVector,
        deg: i32,
        da: &MathMatrix,
    ) -> &MathMatrix {
        derivative::constraint_derivative(self, ssp, parameters, deg, da)
    }

    /// `NbConstraints(SSP, FirstPoint, LastPoint, TheConstraints)`
    /// (`gxx:510-547`): `CCol * IncPass + IncTan * (CCol - 1) + 3 * IncCurv`.
    /// The two point arguments are unused in OCCT.
    pub fn nb_constraints<L: LineTool>(ssp: &L, constraints: &[ConstraintCouple]) -> usize {
        let mut inc_pass = 0i32;
        let mut inc_tan = 0i32;
        let mut inc_curv = 0i32;
        for couple in constraints {
            let cons = couple.constraint();
            if at_least_pass(cons) {
                inc_pass += 1;
            }
            if at_least_tangency(cons) {
                inc_tan += 1;
            }
            if cons == Constraint::CurvaturePoint {
                inc_curv += 1;
            }
        }
        let ccol = Self::ccol(ssp);
        (ccol * inc_pass + inc_tan * (ccol - 1) + 3 * inc_curv).max(0) as usize
    }

    /// `NbColumns(SSP, Deg)` (`gxx:549-556`): `CCol * (Deg + 1)`.
    pub fn nb_columns<L: LineTool>(ssp: &L, deg: i32) -> usize {
        (Self::ccol(ssp) * (deg + 1)).max(0) as usize
    }

    /// `CCol = nb3d * 3 + nb2d * 2` (`gxx:123`).
    fn ccol<L: LineTool>(ssp: &L) -> i32 {
        ssp.nb_p3d() * 3 + ssp.nb_p2d() * 2
    }
}



