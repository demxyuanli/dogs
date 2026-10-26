//! `AppDef_Variational` — smooth a `MultiLine` by variational optimization.
//!
//! Source: `ModelingData/TKGeomBase/AppDef/AppDef_Variational.cxx:67-1156`
//! (constructor, accessors, setters) and `AppDef_Variational.hxx:41-330`
//! (state).
//!
//! Translation notes:
//! - `occ::handle<NCollection_HArray1<T>>` fields are 1-based `Vec<T>`; the
//!   element at index 0 is unused so that the C++ indices translate as-is.
//!   `myParameters` keeps the point-index bounds, so it is a [`MathVector`]
//!   with `lower() == FirstPoint`.
//! - `mySmoothCriterion` is a `Box<dyn SmoothCriterion>`; the criterion holds
//!   its own `Rc<RefCell<Curve>>`, which the driver reuses through
//!   `SmoothCriterion::get_curve` (`cxx:503`).
//! - OCCT throws `Standard_DomainError` for the argument checks at
//!   `cxx:72-98` and `Standard_ConstructionError` out of `Init`; the port
//!   returns `Err(())`.
//! - `Dump` (`cxx:910`) is not ported: it needs `Standard_OStream`.

pub mod approximate;
mod cutting;
mod estimations;
mod init_setup;
mod motor;
mod optimization;
mod project;
mod split;

use occt_core::kernel::geomabs::Shape;
use occt_math::vector::MathVector;

use crate::app_def::linear_criteria::LinearCriteria;
use crate::app_def::multi_line::MultiLine;
use crate::app_def::my_line_tool::MyLineTool;
use crate::app_def::smooth_criterion::SmoothCriterion;
use crate::app_par_curves::constraint::ConstraintCouple;
use crate::app_par_curves::multi_bspline_curve::MultiBSpCurve;

/// `Precision::Angular()` (`Precision.hxx`).
pub(super) const ANGULAR: f64 = 1.0e-12;

/// `Precision::Confusion()` (`Precision.hxx`).
pub(super) const CONFUSION: f64 = 1.0e-7;

/// `AppDef_Variational` (`hxx:41-330`).
pub struct Variational {
    pub(super) ssp: MultiLine,
    pub(super) nb_p3d: i32,
    pub(super) nb_p2d: i32,
    pub(super) dimension: i32,
    pub(super) first_point: i32,
    pub(super) last_point: i32,
    pub(super) nb_points: i32,
    /// `myTabPoints` (1-based, `myDimension * myNbPoints` entries).
    pub(super) tab_points: Vec<f64>,
    pub(super) constraints: Vec<ConstraintCouple>,
    pub(super) nb_constraints: i32,
    /// `myTabConstraints` (1-based).
    pub(super) tab_constraints: Vec<f64>,
    pub(super) nb_pass_points: i32,
    pub(super) nb_tang_points: i32,
    pub(super) nb_curv_points: i32,
    /// `myTypConstraints` (1-based).
    pub(super) typ_constraints: Vec<i32>,
    /// `myTtheta` (1-based).
    pub(super) ttheta: Vec<f64>,
    /// `myTfthet` (1-based).
    pub(super) tfthet: Vec<f64>,
    pub(super) max_degree: i32,
    pub(super) max_segment: i32,
    pub(super) nb_iterations: i32,
    pub(super) tolerance: f64,
    pub(super) continuity: Shape,
    pub(super) niv_cont: i32,
    pub(super) with_min_max: bool,
    pub(super) with_cutting: bool,
    pub(super) percent: [f64; 3],
    /// `myCriterium`: `[Quadratic, Tension, Flexion, Jerk]`.
    pub(super) criterium: [f64; 4],
    pub(super) smooth_criterion: Box<dyn SmoothCriterion>,
    pub(super) parameters: MathVector,
    /// `myKnots` (1-based).
    pub(super) knots: Vec<f64>,
    pub(super) mbsp_curve: MultiBSpCurve,
    pub(super) max_error: f64,
    pub(super) max_error_index: i32,
    pub(super) average_error: f64,
    pub(super) is_created: bool,
    pub(super) is_done: bool,
    pub(super) is_over_constr: bool,
}

impl Variational {
    /// `AppDef_Variational(...)` (`cxx:67-221`).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        ssp: &MultiLine,
        first_point: i32,
        last_point: i32,
        the_constraints: &[ConstraintCouple],
        max_degree: i32,
        max_segment: i32,
        continuity: Shape,
        with_min_max: bool,
        with_cutting: bool,
        tolerance: f64,
        nb_iterations: i32,
    ) -> Result<Self, ()> {
        // `cxx:95-98`: `myMaxDegree < 1` is a domain error.
        if max_degree < 1 {
            return Err(());
        }
        let max_degree = max_degree.min(30); // `cxx:100`.
        if max_segment < 1 {
            return Err(()); // `cxx:102-105`.
        }
        // `cxx:119-132`.
        let niv_cont = match continuity {
            Shape::C0 => 0,
            Shape::C1 => 1,
            Shape::C2 => 2,
            _ => return Err(()),
        };

        let nb_p2d = MyLineTool::nb_p2d(ssp);
        let nb_p3d = MyLineTool::nb_p3d(ssp);
        let dimension = 2 * nb_p2d + 3 * nb_p3d; // `cxx:136`.

        let percent = [0.4, 0.2, 0.4]; // `cxx:138-140`.

        // `cxx:141-143`.
        let knots = vec![0.0, 0.0, 1.0];

        // `cxx:147`.
        let smooth_criterion: Box<dyn SmoothCriterion> =
            Box::new(LinearCriteria::new(ssp, first_point, last_point));

        // `cxx:148-153`: `myParameters` is indexed by point number.
        let parameters = MathVector::new(first_point as usize, last_point as usize);
        let nb_points = last_point - first_point + 1;
        if nb_points <= 0 {
            return Err(());
        }

        let mut tab_points = vec![0.0f64; (dimension * nb_points) as usize + 1];
        let mut index = 1usize;

        // `cxx:166-217`: the three arms are exclusive.
        for ipoint in first_point..=last_point {
            if nb_p2d != 0 && nb_p3d == 0 {
                let (_, tab_p2d) = MyLineTool::value(ssp, ipoint);
                for jp2d in 0..nb_p2d as usize {
                    tab_points[index] = tab_p2d[jp2d].x();
                    index += 1;
                    tab_points[index] = tab_p2d[jp2d].y();
                    index += 1;
                }
            }
            if nb_p3d != 0 && nb_p2d == 0 {
                let (tab_p3d, _) = MyLineTool::value(ssp, ipoint);
                for jp3d in 0..nb_p3d as usize {
                    tab_points[index] = tab_p3d[jp3d].x();
                    index += 1;
                    tab_points[index] = tab_p3d[jp3d].y();
                    index += 1;
                    tab_points[index] = tab_p3d[jp3d].z();
                    index += 1;
                }
            }
            if nb_p3d != 0 && nb_p2d != 0 {
                let (tab_p3d, tab_p2d) = MyLineTool::value(ssp, ipoint);
                for jp3d in 0..nb_p3d as usize {
                    tab_points[index] = tab_p3d[jp3d].x();
                    index += 1;
                    tab_points[index] = tab_p3d[jp3d].y();
                    index += 1;
                    tab_points[index] = tab_p3d[jp3d].z();
                    index += 1;
                }
                for jp2d in 0..nb_p2d as usize {
                    tab_points[index] = tab_p2d[jp2d].x();
                    index += 1;
                    tab_points[index] = tab_p2d[jp2d].y();
                    index += 1;
                }
            }
        }

        let mut this = Variational {
            ssp: ssp.clone(),
            nb_p3d,
            nb_p2d,
            dimension,
            first_point,
            last_point,
            nb_points,
            tab_points,
            constraints: the_constraints.to_vec(),
            nb_constraints: 0,
            tab_constraints: Vec::new(),
            nb_pass_points: 0,
            nb_tang_points: 0,
            nb_curv_points: 0,
            typ_constraints: Vec::new(),
            ttheta: Vec::new(),
            tfthet: Vec::new(),
            max_degree,
            max_segment,
            nb_iterations,
            tolerance,
            continuity,
            niv_cont,
            with_min_max,
            with_cutting,
            percent,
            criterium: [0.0; 4], // `cxx:210`.
            smooth_criterion,
            parameters,
            knots,
            mbsp_curve: MultiBSpCurve::new(),
            max_error: 0.0,
            max_error_index: 0,
            average_error: 0.0,
            is_created: false,
            is_done: false,
            is_over_constr: false,
        };
        this.init()?; // `cxx:218`.
        Ok(this)
    }

    /// `AppDef_Variational(SSP, FirstPoint, LastPoint, TheConstraints)`
    /// (`AppDef_Variational.hxx:62-74`) using the C++ default arguments:
    /// `MaxDegree = 14`, `MaxSegment = 100`, `Continuity = GeomAbs_C2`,
    /// `WithMinMax = false`, `WithCutting = true`, `Tolerance = 1.0`,
    /// `NbIterations = 2`.
    pub fn with_defaults(
        ssp: &MultiLine,
        first_point: i32,
        last_point: i32,
        the_constraints: &[ConstraintCouple],
    ) -> Result<Self, ()> {
        Self::new(
            ssp,
            first_point,
            last_point,
            the_constraints,
            14,
            100,
            Shape::C2,
            false,
            true,
            1.0,
            2,
        )
    }

    /// `IsCreated()` (`cxx:665-668`).
    pub fn is_created(&self) -> bool {
        self.is_created
    }

    /// `IsDone()` (`cxx:673-676`).
    pub fn is_done(&self) -> bool {
        self.is_done
    }

    /// `IsOverConstrained()` (`cxx:681-684`).
    pub fn is_over_constrained(&self) -> bool {
        self.is_over_constr
    }

    /// `Value()` (`cxx:687-694`).
    pub fn value(&self) -> Result<MultiBSpCurve, ()> {
        if !self.is_done {
            return Err(());
        }
        Ok(self.mbsp_curve.clone())
    }

    /// `MaxError()` (`cxx:699-707`).
    pub fn max_error(&self) -> Result<f64, ()> {
        if !self.is_done {
            return Err(());
        }
        Ok(self.max_error)
    }

    /// `MaxErrorIndex()` (`cxx:711-718`).
    pub fn max_error_index(&self) -> Result<i32, ()> {
        if !self.is_done {
            return Err(());
        }
        Ok(self.max_error_index)
    }

    /// `QuadraticError()` (`cxx:723-731`).
    pub fn quadratic_error(&self) -> Result<f64, ()> {
        if !self.is_done {
            return Err(());
        }
        Ok(self.criterium[0])
    }

    /// `Distance(mat)` (`cxx:735-787`).
    pub fn distance(&self, mat: &mut occt_math::matrix::MathMatrix) -> Result<(), ()> {
        if !self.is_done {
            return Err(());
        }
        self.distance_impl(mat);
        Ok(())
    }

    /// `AverageError()` (`cxx:792-797`).
    pub fn average_error(&self) -> Result<f64, ()> {
        if !self.is_done {
            return Err(());
        }
        Ok(self.average_error)
    }

    /// `Parameters()` (`cxx:803-812`) — indexed by point number.
    pub fn parameters(&self) -> Result<&MathVector, ()> {
        if !self.is_done {
            return Err(());
        }
        Ok(&self.parameters)
    }

    /// `Knots()` (1-based) (`cxx:816-825`).
    pub fn knots(&self) -> Result<&[f64], ()> {
        if !self.is_done {
            return Err(());
        }
        Ok(&self.knots)
    }

    /// `Criterium(VFirstOrder, VSecondOrder, VThirdOrder)` (`cxx:828-841`).
    pub fn criterium(&self) -> Result<(f64, f64, f64), ()> {
        if !self.is_done {
            return Err(());
        }
        Ok((self.criterium[1], self.criterium[2], self.criterium[3]))
    }

    /// `CriteriumWeight(Percent1, Percent2, Percent3)` (`cxx:844-849`).
    pub fn criterium_weight(&self) -> [f64; 3] {
        self.percent
    }

    /// `MaxDegree()` (`cxx:854-857`).
    pub fn max_degree(&self) -> i32 {
        self.max_degree
    }

    /// `MaxSegment()` (`cxx:862-865`).
    pub fn max_segment(&self) -> i32 {
        self.max_segment
    }

    /// `Continuity()` (`cxx:870-873`).
    pub fn continuity(&self) -> Shape {
        self.continuity
    }

    /// `WithMinMax()` (`cxx:878-881`).
    pub fn with_min_max(&self) -> bool {
        self.with_min_max
    }

    /// `WithCutting()` (`cxx:886-889`).
    pub fn with_cutting(&self) -> bool {
        self.with_cutting
    }

    /// `Tolerance()` (`cxx:894-897`).
    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }

    /// `NbIterations()` (`cxx:902-905`).
    pub fn nb_iterations(&self) -> i32 {
        self.nb_iterations
    }

    /// `SetConstraints(aConstraint)` (`cxx:951-959`).
    pub fn set_constraints(&mut self, the_constraints: &[ConstraintCouple]) -> Result<bool, ()> {
        self.constraints = the_constraints.to_vec();
        self.init()?;
        Ok(!self.is_over_constr)
    }

    /// `SetParameters(param)` (`cxx:963-966`).
    pub fn set_parameters(&mut self, param: &MathVector) {
        debug_assert_eq!(param.lower(), self.parameters.lower());
        debug_assert_eq!(param.upper(), self.parameters.upper());
        self.parameters = param.clone();
    }

    /// `SetKnots(knots)` (`cxx:971-976`).
    pub fn set_knots(&mut self, knots: &[f64]) -> bool {
        self.knots = knots.to_vec();
        true
    }

    /// `SetMaxDegree(Degree)` (`cxx:980-995`).
    pub fn set_max_degree(&mut self, degree: i32) -> Result<bool, ()> {
        let bound = (degree - self.niv_cont) * self.max_segment
            - self.nb_pass_points
            - 2 * self.nb_tang_points
            - 3 * self.nb_curv_points;
        if bound < 0 {
            return Ok(false);
        }
        self.max_degree = degree;
        self.init_smooth_criterion()?;
        Ok(true)
    }

    /// `SetMaxSegment(NbSegment)` (`cxx:1001-1014`).
    pub fn set_max_segment(&mut self, nb_segment: i32) -> bool {
        if self.with_cutting {
            let bound = (self.max_degree - self.niv_cont) * nb_segment
                - self.nb_pass_points
                - 2 * self.nb_tang_points
                - 3 * self.nb_curv_points;
            if bound < 0 {
                return false;
            }
        }
        self.max_segment = nb_segment;
        true
    }

    /// `SetContinuity(C)` (`cxx:1020-1052`).
    pub fn set_continuity(&mut self, c: Shape) -> Result<bool, ()> {
        let niv_cont = match c {
            Shape::C0 => 0,
            Shape::C1 => 1,
            Shape::C2 => 2,
            _ => return Err(()),
        };
        let bound = (self.max_degree - niv_cont) * self.max_segment
            - self.nb_pass_points
            - 2 * self.nb_tang_points
            - 3 * self.nb_curv_points;
        if bound < 0 {
            return Ok(false);
        }
        self.continuity = c;
        self.niv_cont = niv_cont;
        self.init_smooth_criterion()?;
        Ok(true)
    }

    /// `SetWithMinMax(MinMax)` (`cxx:1056-1061`).
    pub fn set_with_min_max(&mut self, min_max: bool) -> Result<(), ()> {
        self.with_min_max = min_max;
        self.init_smooth_criterion()
    }

    /// `SetWithCutting(Cutting)` (`cxx:1066-1098`).
    pub fn set_with_cutting(&mut self, cutting: bool) -> Result<bool, ()> {
        let bound = if !cutting {
            (self.max_degree - self.niv_cont) * self.knots.len() as i32
                - self.nb_pass_points
                - 2 * self.nb_tang_points
                - 3 * self.nb_curv_points
        } else {
            (self.max_degree - self.niv_cont) * self.max_segment
                - self.nb_pass_points
                - 2 * self.nb_tang_points
                - 3 * self.nb_curv_points
        };
        if bound < 0 {
            return Ok(false);
        }
        self.with_cutting = cutting;
        self.init_smooth_criterion()?;
        Ok(true)
    }

    /// `SetCriteriumWeight(Percent1, Percent2, Percent3)` (`cxx:1103-1117`).
    pub fn set_criterium_weight3(
        &mut self,
        percent1: f64,
        percent2: f64,
        percent3: f64,
    ) -> Result<(), ()> {
        if percent1 < 0.0 || percent2 < 0.0 || percent3 < 0.0 {
            return Err(());
        }
        let total = percent1 + percent2 + percent3;
        self.percent = [percent1 / total, percent2 / total, percent3 / total];
        self.init_smooth_criterion()
    }

    /// `SetCriteriumWeight(Order, Percent)` (`cxx:1122-1140`).
    pub fn set_criterium_weight(&mut self, order: i32, percent: f64) -> Result<(), ()> {
        if percent < 0.0 {
            return Err(());
        }
        if !(1..=3).contains(&order) {
            return Err(());
        }
        self.percent[(order - 1) as usize] = percent;
        let total = self.percent[0] + self.percent[1] + self.percent[2];
        self.percent[0] /= total;
        self.percent[1] /= total;
        self.percent[2] /= total;
        self.init_smooth_criterion()
    }

    /// `SetTolerance(Tol)` (`cxx:1144-1148`).
    pub fn set_tolerance(&mut self, tol: f64) -> Result<(), ()> {
        self.tolerance = tol;
        self.init_smooth_criterion()
    }

    /// `SetNbIterations(Iter)` (`cxx:1153-1156`).
    pub fn set_nb_iterations(&mut self, iter: i32) {
        self.nb_iterations = iter;
    }

    /// `Distance(mat)` body (`cxx:735-787`).
    fn distance_impl(&self, mat: &mut occt_math::matrix::MathMatrix) {
        let j0 = mat.lower_col() as i32 - self.first_point;
        for ipoint in self.first_point..=self.last_point {
            let mut index = 1usize;
            if self.nb_p3d != 0 {
                let (tab_p3d, _) = MyLineTool::value(&self.ssp, ipoint);
                for jp3d in 0..self.nb_p3d as usize {
                    let p3d = tab_p3d[jp3d];
                    let pt3d = self.mbsp_curve.value3d(index as i32, self.param(ipoint));
                    mat.set_value(
                        index as usize,
                        (j0 + ipoint) as usize,
                        p3d.distance(&pt3d),
                    );
                    index += 1;
                }
            }
            if self.nb_p2d != 0 {
                let (tab_p3d, tab_p2d) = MyLineTool::value(&self.ssp, ipoint);
                let _ = tab_p3d;
                for jp2d in 0..self.nb_p2d as usize {
                    let p2d = tab_p2d[jp2d];
                    let pt2d = self.mbsp_curve.value2d(index as i32, self.param(ipoint));
                    mat.set_value(index as usize, (j0 + ipoint) as usize, p2d.distance(&pt2d));
                    index += 1;
                }
            }
        }
    }

    /// `myParameters->Value(ipoint)`.
    pub(super) fn param(&self, ipoint: i32) -> f64 {
        self.parameters.value(ipoint as usize)
    }

    /// `myParameters->SetValue(ipoint, value)`.
    pub(super) fn set_param(&mut self, ipoint: i32, value: f64) {
        self.parameters.set_value(ipoint as usize, value);
    }
}
