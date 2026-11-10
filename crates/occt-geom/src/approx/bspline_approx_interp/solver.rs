use occt_core::bspl::banded_interp::knot_sequence;
use occt_core::bspl::eval_basis::eval_bspline_basis;
use occt_core::bspl::knots::knot_sequence_length;
use occt_core::gp::{GpPnt, GpVec};
use occt_math::{MathMatrix, MathVector};

use crate::bspline_curve::GeomBSplineCurve;
use crate::curve::Curve;

use super::constants::{
    BOUNDARY_KNOT_COUNT, CONT_COUNT_C0_C1_C2, CONT_COUNT_C1_C2, CONT_ROW_C0, CONT_ROW_C1,
    CONT_ROW_C2, ERROR_NORM_EPS, MAX_PROJECTION_ITER, NEWTON_D2F_EPS, OLD_ERROR_SCALE,
    PARAM_ALPHA_DIVISOR, PARAM_LAST_VALUE,
};
use super::linear_algebra::{insert_kink_knot, math_lin_lu, matrix_times_vector, set_block, solve_lu};
use super::matrix::matrix_multiply;

/// `Approx_BSplineApproxInterp` (`Approx_BSplineApproxInterp.hxx:64-236`).
#[derive(Clone)]
pub struct BSplineApproxInterp {
    /// `myPoints`.
    points: Vec<GpPnt>,
    /// `myInterpolated`: 0-based point indices.
    interpolated: Vec<usize>,
    /// `myApproximated`.
    approximated: Vec<usize>,
    /// `myKinks`.
    kinks: Vec<usize>,
    /// `myCurve`.
    curve: Option<GeomBSplineCurve>,
    /// `myDegree`.
    degree: i32,
    /// `myNbControlPts`.
    nb_control_pts: i32,
    /// `myMaxError`.
    max_error: f64,
    /// `myAlpha`.
    alpha: f64,
    /// `myMinPivot`.
    min_pivot: f64,
    /// `myClosedRelTol`.
    closed_rel_tol: f64,
    /// `myKnotInsertTol`.
    knot_insert_tol: f64,
    /// `myConvergenceTol`.
    convergence_tol: f64,
    /// `myProjectionTol`.
    projection_tol: f64,
    /// `myContinuousIfClosed`.
    continuous_if_closed: bool,
    /// `myIsDone`.
    is_done: bool,
}

impl BSplineApproxInterp {
    /// `Approx_BSplineApproxInterp(thePoints, theNbControlPts, theDegree,
    /// theContinuousIfClosed)` (`cxx:123-149`).
    pub fn new(
        points: &[GpPnt],
        nb_control_pts: i32,
        degree: i32,
        continuous_if_closed: bool,
    ) -> Self {
        let mut my_points: Vec<GpPnt> = points.to_vec();
        if my_points.is_empty() {
            my_points.push(GpPnt::new(0.0, 0.0, 0.0));
        }
        let mut this = Self {
            points: my_points,
            interpolated: Vec::new(),
            approximated: Vec::new(),
            kinks: Vec::new(),
            curve: None,
            degree: degree.max(1),
            nb_control_pts: nb_control_pts.max(degree + 1),
            max_error: 0.0,
            alpha: 0.5,
            min_pivot: 1.0e-20,
            closed_rel_tol: 1.0e-12,
            knot_insert_tol: 1.0e-4,
            convergence_tol: 1.0e-3,
            projection_tol: 1.0e-6,
            continuous_if_closed,
            is_done: false,
        };

        // Clamp the control point count and, if needed, the degree
        // (`cxx:132-139`).
        let nb_points = points.len() as i32;
        this.nb_control_pts = this.nb_control_pts.min(nb_points);
        if this.nb_control_pts < this.degree + 1 {
            this.degree = (this.nb_control_pts - 1).max(1);
        }

        for i in 0..points.len() {
            this.approximated.push(i);
        }
        this
    }

    /// `InterpolatePoint(thePointIndex, theWithKink)` (`cxx:152-173`).
    pub fn interpolate_point(&mut self, point_index: i32, with_kink: bool) {
        let index = point_index as usize;
        for i in 0..self.approximated.len() {
            if self.approximated[i] == index {
                // Swap with last and remove last (order does not matter).
                if i < self.approximated.len() - 1 {
                    self.approximated[i] = self.approximated[self.approximated.len() - 1];
                }
                self.approximated.pop();
                self.interpolated.push(index);
                if with_kink {
                    self.kinks.push(index);
                }
                return;
            }
        }
    }

    /// `Perform()` (`cxx:177-181`).
    pub fn perform(&mut self) {
        let params = self.compute_parameters(self.alpha);
        self.perform_with_params(&params);
    }

    /// `Perform(theParams)` (`cxx:184-201`).
    pub fn perform_with_params(&mut self, params: &[f64]) {
        self.is_done = false;
        self.max_error = 0.0;
        self.curve = None;

        if params.len() != self.points.len() {
            return;
        }

        let (knots, mults) = self.compute_knots(self.nb_control_pts, params);
        self.is_done = self.solve(params, &knots, &mults);
    }

    /// `PerformOptimal(theMaxIter)` (`cxx:205-209`).
    pub fn perform_optimal(&mut self, max_iter: i32) {
        let params = self.compute_parameters(self.alpha);
        self.perform_optimal_with_params(&params, max_iter);
    }

    /// `PerformOptimal(theParams, theMaxIter)` (`cxx:212-260`).
    pub fn perform_optimal_with_params(&mut self, params: &[f64], max_iter: i32) {
        self.is_done = false;
        self.max_error = 0.0;
        self.curve = None;

        if params.len() != self.points.len() {
            return;
        }

        let mut a_params = params.to_vec();
        let (knots, mults) = self.compute_knots(self.nb_control_pts, &a_params);

        if !self.solve(&a_params, &knots, &mults) {
            return;
        }

        let mut old_error = self.max_error * OLD_ERROR_SCALE;
        let mut iter = 0;
        while self.max_error > 0.0
            && (old_error - self.max_error) / self.max_error.max(ERROR_NORM_EPS)
                > self.convergence_tol
            && iter < max_iter
        {
            old_error = self.max_error;
            self.optimize_parameters(&mut a_params);
            if !self.solve(&a_params, &knots, &mults) {
                return;
            }
            iter += 1;
        }
        self.is_done = true;
    }

    /// `Curve()` (`cxx:266-269`).
    pub fn curve(&self) -> Option<&GeomBSplineCurve> {
        self.curve.as_ref()
    }

    /// `IsDone()` (`hxx:106`).
    pub fn is_done(&self) -> bool {
        self.is_done
    }

    /// `MaxError()` (`hxx:111`).
    pub fn max_error(&self) -> f64 {
        self.max_error
    }

    /// `SetParametrizationAlpha` (`hxx:117`).
    pub fn set_parametrization_alpha(&mut self, alpha: f64) {
        self.alpha = alpha;
    }

    /// `SetMinPivot` (`hxx:122`).
    pub fn set_min_pivot(&mut self, min_pivot: f64) {
        self.min_pivot = min_pivot;
    }

    /// `SetClosedTolerance` (`hxx:128`).
    pub fn set_closed_tolerance(&mut self, rel_tol: f64) {
        self.closed_rel_tol = rel_tol;
    }

    /// `SetKnotInsertionTolerance` (`hxx:133`).
    pub fn set_knot_insertion_tolerance(&mut self, tol: f64) {
        self.knot_insert_tol = tol;
    }

    /// `SetConvergenceTolerance` (`hxx:137`).
    pub fn set_convergence_tolerance(&mut self, tol: f64) {
        self.convergence_tol = tol;
    }

    /// `SetProjectionTolerance` (`hxx:141`).
    pub fn set_projection_tolerance(&mut self, tol: f64) {
        self.projection_tol = tol;
    }

    /// `computeParameters(theAlpha)` (`cxx:272-297`).
    fn compute_parameters(&self, alpha: f64) -> Vec<f64> {
        let nb_pts = self.points.len();
        let mut params = vec![0.0f64; nb_pts];
        let mut sum = 0.0;
        for i in 2..=nb_pts {
            let diff = GpVec::from_pnts(&self.points[i - 2], &self.points[i - 1]);
            let len2 = diff.square_magnitude();
            sum += len2.powf(alpha / PARAM_ALPHA_DIVISOR);
            params[i - 1] = sum;
        }
        let tmax = params[nb_pts - 1];
        if tmax > 0.0 {
            for i in 1..nb_pts {
                params[i] /= tmax;
            }
        }
        params[nb_pts - 1] = PARAM_LAST_VALUE;
        params
    }

    /// `computeKnots(theNbCP, theParams, theKnots, theMults)` (`cxx:302-355`).
    fn compute_knots(&self, nb_cp: i32, params: &[f64]) -> (Vec<f64>, Vec<i32>) {
        let order = self.degree + 1;
        if nb_cp < order || params.is_empty() {
            return (Vec::new(), Vec::new());
        }

        let mut u_min = params[0];
        let mut u_max = params[0];
        for value in params.iter().skip(1) {
            u_min = u_min.min(*value);
            u_max = u_max.max(*value);
        }

        let nb_internal = nb_cp - order;
        let nb_base = nb_internal + BOUNDARY_KNOT_COUNT;
        let mut knots = vec![0.0f64; nb_base as usize];
        let mut mults = vec![0i32; nb_base as usize];

        knots[0] = u_min;
        mults[0] = order;
        for i in 1..=nb_internal {
            knots[i as usize] =
                u_min + (u_max - u_min) * (i as f64) / ((nb_internal + 1) as f64);
            mults[i as usize] = 1;
        }
        knots[(nb_base - 1) as usize] = u_max;
        mults[(nb_base - 1) as usize] = order;

        for kink in &self.kinks {
            let kink_param = params[*kink];
            insert_kink_knot(
                kink_param,
                self.degree,
                &mut knots,
                &mut mults,
                self.knot_insert_tol,
            );
        }
        (knots, mults)
    }

    /// `buildBasisMatrix(theFlatKnots, theParams, theDerivOrder)`
    /// (`cxx:358-405`).
    fn build_basis_matrix(
        &self,
        flat_knots: &[f64],
        params: &[f64],
        deriv_order: i32,
    ) -> MathMatrix {
        let nb_cp = flat_knots.len() as i32 - self.degree - 1;
        let mut basis_mat = MathMatrix::with_init(1, params.len(), 1, nb_cp as usize, 0.0);

        for (i_parm, param) in params.iter().enumerate() {
            let row = i_parm + 1;
            let Ok((first_non_zero, basis)) = eval_bspline_basis(
                deriv_order,
                self.degree + 1,
                flat_knots,
                *param,
                false,
            ) else {
                continue;
            };
            let values = &basis[deriv_order as usize];
            for j in 0..=(self.degree as usize) {
                let col = first_non_zero as usize + j;
                if col >= 1 && col <= nb_cp as usize {
                    basis_mat.set_value(row, col, values[j]);
                }
            }
        }
        basis_mat
    }

    /// `buildContinuityMatrix(...)` (`cxx:409-456`).
    fn build_continuity_matrix(
        &self,
        nb_ctr_pnts: i32,
        nb_continuity: usize,
        params: &[f64],
        flat_knots: &[f64],
    ) -> MathMatrix {
        let mut continuity =
            MathMatrix::with_init(1, nb_continuity, 1, nb_ctr_pnts as usize, 0.0);

        let first_param = params[0];
        let last_param = params[params.len() - 1];

        // C1 condition: first derivative at start equals first derivative at end.
        let diff1_first = self.build_basis_matrix(flat_knots, &[first_param], CONT_ROW_C1 as i32);
        let diff1_last = self.build_basis_matrix(flat_knots, &[last_param], CONT_ROW_C1 as i32);
        for j in 1..=nb_ctr_pnts as usize {
            continuity.set_value(
                CONT_ROW_C1,
                j,
                diff1_first.value(1, j) - diff1_last.value(1, j),
            );
        }

        // C2 condition.
        if nb_continuity >= CONT_ROW_C2 {
            let diff2_first =
                self.build_basis_matrix(flat_knots, &[first_param], CONT_ROW_C2 as i32);
            let diff2_last =
                self.build_basis_matrix(flat_knots, &[last_param], CONT_ROW_C2 as i32);
            for j in 1..=nb_ctr_pnts as usize {
                continuity.set_value(
                    CONT_ROW_C2,
                    j,
                    diff2_first.value(1, j) - diff2_last.value(1, j),
                );
            }
        }

        // C0 condition (only if first/last are not already interpolated).
        if nb_continuity >= CONT_ROW_C0 {
            let diff0_first = self.build_basis_matrix(flat_knots, &[first_param], 0);
            let diff0_last = self.build_basis_matrix(flat_knots, &[last_param], 0);
            for j in 1..=nb_ctr_pnts as usize {
                continuity.set_value(
                    CONT_ROW_C0,
                    j,
                    diff0_first.value(1, j) - diff0_last.value(1, j),
                );
            }
        }

        continuity
    }

    /// `solve(theParams, theKnots, theMults)` (`cxx:460-651`).
    fn solve(&mut self, params: &[f64], knots: &[f64], mults: &[i32]) -> bool {
        let nb_flat_knots = knot_sequence_length(mults, self.degree, false);
        let flat_knots = knot_sequence(knots, mults, self.degree);
        debug_assert_eq!(flat_knots.len(), nb_flat_knots);

        let nb_approx = self.approximated.len();
        let nb_interp = self.interpolated.len();

        let make_closed = self.is_closed();
        let mut nb_continuity = 0usize;
        if make_closed {
            nb_continuity = CONT_COUNT_C0_C1_C2;
            if self.is_first_and_last_interpolated() {
                nb_continuity = CONT_COUNT_C1_C2;
            }
        }

        let nb_ctr_pnts = nb_flat_knots as i32 - self.degree - 1;
        if nb_ctr_pnts < (nb_interp + nb_continuity) as i32
            || nb_ctr_pnts < self.degree + 1 + nb_continuity as i32
        {
            return false;
        }
        if nb_approx == 0 && nb_ctr_pnts != (nb_interp + nb_continuity) as i32 {
            return false;
        }
        let nb_ctr_pnts = nb_ctr_pnts as usize;

        let nb_vars = nb_ctr_pnts + nb_interp + nb_continuity;
        let mut lhs = MathMatrix::with_init(1, nb_vars, 1, nb_vars, 0.0);

        let mut rhs_x = MathVector::new(1, nb_vars);
        let mut rhs_y = MathVector::new(1, nb_vars);
        let mut rhs_z = MathVector::new(1, nb_vars);

        // Upper-left block: A^T*A and the right-hand side A^T*b.
        if nb_approx > 0 {
            let mut bx = MathVector::new(1, nb_approx);
            let mut by = MathVector::new(1, nb_approx);
            let mut bz = MathVector::new(1, nb_approx);
            let mut approx_params = vec![0.0f64; nb_approx];

            for i in 0..nb_approx {
                let idx = self.approximated[i];
                let pnt = &self.points[idx];
                bx.set_value(i + 1, pnt.x());
                by.set_value(i + 1, pnt.y());
                bz.set_value(i + 1, pnt.z());
                approx_params[i] = params[idx];
            }

            let a = self.build_basis_matrix(&flat_knots, &approx_params, 0);
            let at = a.transposed();

            set_block(&mut lhs, 1, nb_ctr_pnts, 1, nb_ctr_pnts, &matrix_multiply(&at, &a));

            let at_bx = matrix_times_vector(&at, &bx);
            let at_by = matrix_times_vector(&at, &by);
            let at_bz = matrix_times_vector(&at, &bz);
            for i in 1..=nb_ctr_pnts {
                rhs_x.set_value(i, at_bx.value(i));
                rhs_y.set_value(i, at_by.value(i));
                rhs_z.set_value(i, at_bz.value(i));
            }
        }

        // Constraint blocks: C (interpolation) and L (continuity).
        if nb_interp + nb_continuity > 0 {
            let mut dx = MathVector::new(1, nb_interp + nb_continuity);
            let mut dy = MathVector::new(1, nb_interp + nb_continuity);
            let mut dz = MathVector::new(1, nb_interp + nb_continuity);

            if nb_interp > 0 {
                let mut interp_params = vec![0.0f64; nb_interp];
                for i in 0..nb_interp {
                    let idx = self.interpolated[i];
                    let pnt = &self.points[idx];
                    dx.set_value(i + 1, pnt.x());
                    dy.set_value(i + 1, pnt.y());
                    dz.set_value(i + 1, pnt.z());
                    interp_params[i] = params[idx];
                }

                let c = self.build_basis_matrix(&flat_knots, &interp_params, 0);
                let ct = c.transposed();

                set_block(
                    &mut lhs,
                    1,
                    nb_ctr_pnts,
                    nb_ctr_pnts + 1,
                    nb_ctr_pnts + nb_interp,
                    &ct,
                );
                set_block(
                    &mut lhs,
                    nb_ctr_pnts + 1,
                    nb_ctr_pnts + nb_interp,
                    1,
                    nb_ctr_pnts,
                    &c,
                );
            }

            // Continuity constraints for closed curves.
            if make_closed && nb_continuity > 0 {
                let cont_mat =
                    self.build_continuity_matrix(nb_ctr_pnts as i32, nb_continuity, params, &flat_knots);
                let cont_mat_t = cont_mat.transposed();

                let cont_row_start = nb_ctr_pnts + nb_interp + 1;
                let cont_row_end = nb_ctr_pnts + nb_interp + nb_continuity;

                set_block(
                    &mut lhs,
                    cont_row_start,
                    cont_row_end,
                    1,
                    nb_ctr_pnts,
                    &cont_mat,
                );
                set_block(
                    &mut lhs,
                    1,
                    nb_ctr_pnts,
                    cont_row_start,
                    cont_row_end,
                    &cont_mat_t,
                );
            }

            for i in 1..=(nb_interp + nb_continuity) {
                rhs_x.set_value(nb_ctr_pnts + i, dx.value(i));
                rhs_y.set_value(nb_ctr_pnts + i, dy.value(i));
                rhs_z.set_value(nb_ctr_pnts + i, dz.value(i));
            }
        }

        let Some(lu_res) = math_lin_lu(&lhs, self.min_pivot) else {
            return false;
        };

        let cx = solve_lu(&lu_res, &rhs_x);
        let cy = solve_lu(&lu_res, &rhs_y);
        let cz = solve_lu(&lu_res, &rhs_z);

        // Extract control points.
        let mut poles = Vec::with_capacity(nb_ctr_pnts);
        for i in 1..=nb_ctr_pnts {
            poles.push(GpPnt::new(cx.value(i), cy.value(i), cz.value(i)));
        }

        let Ok(curve) =
            GeomBSplineCurve::from_poles_knots_mults(poles, knots.to_vec(), mults.to_vec(), self.degree as usize)
        else {
            return false;
        };

        // Compute the max error at approximated points.
        self.max_error = 0.0;
        for i in 0..nb_approx {
            let idx = self.approximated[i];
            let pnt = &self.points[idx];
            let param = params[idx];
            let error = curve.d0(param).distance(pnt);
            self.max_error = self.max_error.max(error);
        }
        self.curve = Some(curve);
        true
    }

    /// `projectOnCurve` (`cxx:675-721`).
    fn project_on_curve(
        &self,
        pnt: &GpPnt,
        curve: &GeomBSplineCurve,
        init_param: f64,
        param: &mut f64,
    ) -> f64 {
        let first_param = curve.first_parameter();
        let last_param = curve.last_parameter();
        let mut t = init_param;
        let mut f = 0.0;

        for _ in 0..MAX_PROJECTION_ITER {
            let (point, d1, d2) = curve.d2(t);

            let diff = GpVec::from_pnts(&point, pnt);
            f = diff.square_magnitude();

            let df = diff.dot(&d1);
            let d2f = diff.dot(&d2) + d1.square_magnitude();

            if d2f.abs() < NEWTON_D2F_EPS {
                break;
            }

            let dt = -df / d2f;
            let mut new_t = t + dt;
            if new_t < first_param {
                new_t = first_param;
            } else if new_t > last_param {
                new_t = last_param;
            }
            t = new_t;

            if dt.abs() <= self.projection_tol {
                break;
            }
        }

        *param = t;
        f.sqrt()
    }

    /// `optimizeParameters(theCurve, theParams)` (`cxx:725-740`).
    fn optimize_parameters(&mut self, params: &mut [f64]) {
        let Some(curve) = self.curve.clone() else {
            return;
        };
        for i in 0..self.approximated.len() {
            let idx = self.approximated[i];
            let pnt = self.points[idx].clone();
            let old_param = params[idx];
            let mut new_param = old_param;
            self.project_on_curve(&pnt, &curve, old_param, &mut new_param);
            params[idx] = new_param;
        }
    }

    /// `isClosed()` (`cxx:744-754`).
    fn is_closed(&self) -> bool {
        if !self.continuous_if_closed {
            return false;
        }
        let diag = self.bounding_box_diagonal();
        let tol = self.closed_rel_tol * diag;
        self.points[0].distance(&self.points[self.points.len() - 1]) <= tol
    }

    /// `isFirstAndLastInterpolated()` (`cxx:758-776`).
    fn is_first_and_last_interpolated(&self) -> bool {
        let mut first = false;
        let mut last = false;
        let last_idx = self.points.len() - 1;
        for idx in &self.interpolated {
            if *idx == 0 {
                first = true;
            }
            if *idx == last_idx {
                last = true;
            }
        }
        first && last
    }

    /// `boundingBoxDiagonal()` (`cxx:779-800`).
    fn bounding_box_diagonal(&self) -> f64 {
        if self.points.is_empty() {
            return 0.0;
        }
        let mut min = self.points[0].clone();
        let mut max = min.clone();
        for pnt in self.points.iter().skip(1) {
            min = GpPnt::new(
                min.x().min(pnt.x()),
                min.y().min(pnt.y()),
                min.z().min(pnt.z()),
            );
            max = GpPnt::new(
                max.x().max(pnt.x()),
                max.y().max(pnt.y()),
                max.z().max(pnt.z()),
            );
        }
        max.distance(&min)
    }
}
