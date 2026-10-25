//! `BRepBlend_AppFuncRoot` — the `Approx_SweepFunction` that re-solves the
//! walked sections of a `BRepBlend_Line` for the approximation driver.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_AppFuncRoot.hxx:41-185`
//! and `.cxx:26-430` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `void* myFunc` (`hxx:176`) is OCCT's untyped back-pointer to the caller's
//!   `Blend_AppFunction&`; the port keeps that aliasing with
//!   `Rc<RefCell<dyn BlendAppFunction>>`, the same shared-handle model the
//!   walking builders use for their `TopolTool`s.
//! - `Handle(BRepBlend_Line) myLine` (`hxx:175`): `SearchPoint` calls
//!   `myLine->InsertBefore` (`cxx:374`), so the port shares the line as
//!   `Rc<RefCell<Line>>` (`BRepBlend_Line.lxx:26-29`).
//! - The pure virtual `Point` / `Vec` (`hxx:150-156`) are the only thing the
//!   three subclasses (`BRepBlend_AppFunc.cxx:33-44`,
//!   `BRepBlend_AppFuncRst.cxx:32-44`, `BRepBlend_AppFuncRstRst.cxx:33-44`)
//!   override. The port keeps them as the `AppFuncContact` trait object
//!   (`contact`) plus the three constructors [`AppFuncRoot::new_surf_surf`] /
//!   [`AppFuncRoot::new_surf_rst`] / [`AppFuncRoot::new_rst_rst`].
//! - `math_Vector` members (`hxx:177-184`) are [`MathVector`]; `SearchPoint` is
//!   non-`const` and rewrites `myTolerance`, `myPnt`, `X1`, `X2`, `XInit` and
//!   `Sol`, so they live in `RefCell`s. `ApproxSweepFunction` is a `&self`
//!   interface (`occt_core::approx`), so `D0`/`D1`/`D2` reach them the same way
//!   `Approx_SweepApproximation` reaches its own caches.
//! - `mySolver` is a `math_FunctionSetRoot` built over `Func` and `myTolerance`
//!   with the literal 30 iterations of `cxx:79`. OCCT passes a
//!   `Blend_AppFunction&`; the port wraps `&mut dyn BlendAppFunction` in the
//!   sized adapter `FuncSetRef` (`search.rs`) because
//!   `MathFunctionSetRoot::perform_with_bounds` is generic over a sized `F`.
//! - `Func.GetShape(NbPoles, NbKnots, Degree, NbPoles2d)` at `cxx:53` writes
//!   four locals that the constructor never reads again; it is not repeated.
//! - `#ifdef BREPBLEND_DEB` (`cxx:355-358`) and `#ifdef OCCT_DEBUG`
//!   (`cxx:367-371`) traces are not translated.
//! - `Approx_SweepFunction::GetTolerance` receives caller-sized arrays in
//!   OCCT (`cxx:239-250`, length `NbPoles`); the port sizes them from
//!   `GetShape` and returns one vector, the convention of the other
//!   `ApproxSweepFunction` ports.

mod contact;
mod search;

pub use contact::{AppFuncContact, RstRstContact, SurfRstContact, SurfSurfContact};
pub(crate) use search::line_point;

use std::cell::RefCell;
use std::rc::Rc;

use occt_core::approx::{ApproxSweepFunction, SweepD1, SweepD2};
use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_geom::geom_fill::snglr_func::shape_from_u8;
use occt_math::{MathFunctionSetRoot, MathVector};

use crate::blend::{
    BlendAppFunction, BlendFunction, BlendPoint, BlendRstRstFunction, BlendSurfRstFunction,
};
use crate::brep_blend::Line;

use search::FuncSetRef;

/// `math_FunctionSetRoot(Func, myTolerance, 30)` (`cxx:79`).
const SOLVER_ITERATIONS: i32 = 30;

/// `BRepBlend_AppFuncRoot` (`BRepBlend_AppFuncRoot.hxx:41-185`).
pub struct AppFuncRoot {
    /// `myLine` (`hxx:175`).
    line: Rc<RefCell<Line>>,
    /// `myFunc` (`hxx:176`): `void*` in OCCT.
    func: Rc<RefCell<dyn BlendAppFunction>>,
    /// `Point` / `Vec` (`hxx:150-156`): the subclass vtable.
    contact: Box<dyn AppFuncContact>,
    /// `myTolerance` (`hxx:177`).
    tolerance: RefCell<MathVector>,
    /// `myPnt` (`hxx:178`).
    pnt: RefCell<BlendPoint>,
    /// `myBary` (`hxx:179`).
    bary: GpPnt,
    /// `X1` (`hxx:180`).
    x1: RefCell<MathVector>,
    /// `X2` (`hxx:181`).
    x2: RefCell<MathVector>,
    /// `XInit` (`hxx:182`).
    x_init: RefCell<MathVector>,
    /// `Sol` (`hxx:184`).
    sol: RefCell<MathVector>,
    /// `mySolver` (`hxx:185`).
    solver: RefCell<MathFunctionSetRoot>,
    /// `Func.NbVariables()`, read once in the constructor (`cxx:44-46`).
    dim: usize,
}

impl AppFuncRoot {
    /// `BRepBlend_AppFunc(Line, Func, Tol3d, Tol2d)`
    /// (`BRepBlend_AppFunc.cxx:24-30`): surface / surface contact.
    pub fn new_surf_surf<F: BlendFunction + 'static>(
        line: Rc<RefCell<Line>>,
        func: Rc<RefCell<F>>,
        tol3d: f64,
        tol2d: f64,
    ) -> Self {
        Self::new_with_contact(line, func, tol3d, tol2d, Box::new(SurfSurfContact))
    }

    /// `BRepBlend_AppFuncRst(Line, Func, Tol3d, Tol2d)`
    /// (`BRepBlend_AppFuncRst.cxx:24-30`): surface / restriction contact.
    pub fn new_surf_rst<F: BlendSurfRstFunction + 'static>(
        line: Rc<RefCell<Line>>,
        func: Rc<RefCell<F>>,
        tol3d: f64,
        tol2d: f64,
    ) -> Self {
        Self::new_with_contact(line, func, tol3d, tol2d, Box::new(SurfRstContact))
    }

    /// `BRepBlend_AppFuncRstRst(Line, Func, Tol3d, Tol2d)`
    /// (`BRepBlend_AppFuncRstRst.cxx:24-30`): restriction / restriction
    /// contact.
    pub fn new_rst_rst<F: BlendRstRstFunction + 'static>(
        line: Rc<RefCell<Line>>,
        func: Rc<RefCell<F>>,
        tol3d: f64,
        tol2d: f64,
    ) -> Self {
        Self::new_with_contact(line, func, tol3d, tol2d, Box::new(RstRstContact))
    }

    /// `BRepBlend_AppFuncRoot(Line, Func, Tol3d, Tol2d)` (`cxx:26-80`).
    fn new_with_contact(
        line: Rc<RefCell<Line>>,
        func: Rc<RefCell<dyn BlendAppFunction>>,
        tol3d: f64,
        tol2d: f64,
        contact: Box<dyn AppFuncContact>,
    ) -> Self {
        let dim = func.borrow().nb_variables();

        // Tolerances (`cxx:42-50`).
        let mut tolerance = MathVector::new(1, dim);
        func.borrow().get_tolerance(&mut tolerance, tol3d);
        for ii in 1..=dim {
            if tolerance.value(ii) > tol2d {
                tolerance.set_value(ii, tol2d);
            }
        }

        // Calculation of BaryCentre (rational case) (`cxx:56-76`).
        let bary = if func.borrow().is_rational() {
            let (mut xmax, mut xmin): (f64, f64) = (-1.0e100, 1.0e100);
            let (mut ymax, mut ymin): (f64, f64) = (-1.0e100, 1.0e100);
            let (mut zmax, mut zmin): (f64, f64) = (-1.0e100, 1.0e100);
            let mut bary = GpPnt::new(0.0, 0.0, 0.0);
            let line_ref = line.borrow();
            for ii in 1..=line_ref.nb_points() {
                let p = line_point(&line_ref, ii);
                xmax = xmax.max(p.point_on_s1().x().max(p.point_on_s2().x()));
                xmin = xmin.min(p.point_on_s1().x().min(p.point_on_s2().x()));
                ymax = ymax.max(p.point_on_s1().y().max(p.point_on_s2().y()));
                ymin = ymin.min(p.point_on_s1().y().min(p.point_on_s2().y()));
                zmax = zmax.max(p.point_on_s1().z().max(p.point_on_s2().z()));
                zmin = zmin.min(p.point_on_s1().z().min(p.point_on_s2().z()));
                bary = GpPnt::new((xmax + xmin) / 2.0, (ymax + ymin) / 2.0, (zmax + zmin) / 2.0);
            }
            bary
        } else {
            GpPnt::new(0.0, 0.0, 0.0)
        };

        // `mySolver = std::make_unique<math_FunctionSetRoot>(Func, myTolerance, 30)`
        // (`cxx:79`). `math_FunctionSetRoot` keeps no reference to `Func`, only
        // the dimension and the tolerance.
        let solver = {
            let mut f = func.borrow_mut();
            let adapter = FuncSetRef(&mut *f);
            MathFunctionSetRoot::new(&adapter, &tolerance, SOLVER_ITERATIONS)
        };

        Self {
            line,
            func,
            contact,
            tolerance: RefCell::new(tolerance),
            pnt: RefCell::new(BlendPoint::new()),
            bary,
            x1: RefCell::new(MathVector::new(1, dim)),
            x2: RefCell::new(MathVector::new(1, dim)),
            x_init: RefCell::new(MathVector::new(1, dim)),
            sol: RefCell::new(MathVector::new(1, dim)),
            solver: RefCell::new(solver),
            dim,
        }
    }

    /// `GetShape(NbPoles, NbKnots, Degree, NbPoles2d)`
    /// (`Blend_AppFunction.hxx:135-139`) for the owned function.
    fn shape(&self) -> (i32, i32, i32, i32) {
        let (mut nb_poles, mut nb_knots, mut degree, mut nb_poles2d) = (0, 0, 0, 0);
        self.func.borrow_mut().get_shape(
            &mut nb_poles,
            &mut nb_knots,
            &mut degree,
            &mut nb_poles2d,
        );
        (nb_poles, nb_knots, degree, nb_poles2d)
    }

    /// The `Blend_AppFunction&` this function re-solves sections for.
    pub fn function(&self) -> &Rc<RefCell<dyn BlendAppFunction>> {
        &self.func
    }

    /// The walked line the sections are re-solved from.
    pub fn line(&self) -> &Rc<RefCell<Line>> {
        &self.line
    }
}

impl ApproxSweepFunction for AppFuncRoot {
    /// `D0(Param, First, Last, Poles, Poles2d, Weigths)` (`cxx:90-106`).
    fn d0(
        &self,
        param: f64,
        _first: f64,
        _last: f64,
    ) -> Result<(Vec<GpPnt>, Vec<GpPnt2d>, Vec<f64>), String> {
        if !self.search_point(param) {
            return Err("BRepBlend_AppFuncRoot::D0: NotDone".into());
        }
        let (nb_poles, _nb_knots, _degree, nb_poles2d) = self.shape();
        let nb_poles = nb_poles.max(0) as usize;
        let nb_poles2d = nb_poles2d.max(0) as usize;
        let mut poles = vec![GpPnt::zero(); nb_poles];
        let mut poles2d = vec![GpPnt2d::zero(); nb_poles2d];
        let mut weights = vec![0.0; nb_poles];
        let pnt = self.pnt.borrow().clone();
        self.func
            .borrow_mut()
            .section_poles(&pnt, &mut poles, &mut poles2d, &mut weights);
        Ok((poles, poles2d, weights))
    }

    /// `D1(Param, First, Last, ...)` (`cxx:113-135`).
    fn d1(&self, param: f64, _first: f64, _last: f64) -> Result<SweepD1, String> {
        if !self.search_point(param) {
            return Err("BRepBlend_AppFuncRoot::D1: NotDone".into());
        }
        let (nb_poles, _nb_knots, _degree, nb_poles2d) = self.shape();
        let nb_poles = nb_poles.max(0) as usize;
        let nb_poles2d = nb_poles2d.max(0) as usize;
        let mut poles = vec![GpPnt::zero(); nb_poles];
        let mut dpoles = vec![GpVec::zero(); nb_poles];
        let mut poles2d = vec![GpPnt2d::zero(); nb_poles2d];
        let mut dpoles2d = vec![GpVec2d::zero(); nb_poles2d];
        let mut weights = vec![0.0; nb_poles];
        let mut dweights = vec![0.0; nb_poles];
        let pnt = self.pnt.borrow().clone();
        let ok = self.func.borrow_mut().section(
            &pnt,
            &mut poles,
            &mut dpoles,
            &mut poles2d,
            &mut dpoles2d,
            &mut weights,
            &mut dweights,
        );
        if !ok {
            return Err("BRepBlend_AppFuncRoot::D1: NotDone".into());
        }
        Ok(SweepD1 {
            poles,
            dpoles,
            poles2d,
            dpoles2d,
            weights,
            dweights,
        })
    }

    /// `D2(Param, First, Last, ...)` (`cxx:143-180`).
    fn d2(&self, param: f64, _first: f64, _last: f64) -> Result<SweepD2, String> {
        if !self.search_point(param) {
            return Err("BRepBlend_AppFuncRoot::D2: NotDone".into());
        }
        let (nb_poles, _nb_knots, _degree, nb_poles2d) = self.shape();
        let nb_poles = nb_poles.max(0) as usize;
        let nb_poles2d = nb_poles2d.max(0) as usize;
        let mut poles = vec![GpPnt::zero(); nb_poles];
        let mut dpoles = vec![GpVec::zero(); nb_poles];
        let mut d2poles = vec![GpVec::zero(); nb_poles];
        let mut poles2d = vec![GpPnt2d::zero(); nb_poles2d];
        let mut dpoles2d = vec![GpVec2d::zero(); nb_poles2d];
        let mut d2poles2d = vec![GpVec2d::zero(); nb_poles2d];
        let mut weights = vec![0.0; nb_poles];
        let mut dweights = vec![0.0; nb_poles];
        let mut d2weights = vec![0.0; nb_poles];
        let pnt = self.pnt.borrow().clone();
        let ok = self.func.borrow_mut().section_with_derivatives2(
            &pnt,
            &mut poles,
            &mut dpoles,
            &mut d2poles,
            &mut poles2d,
            &mut dpoles2d,
            &mut d2poles2d,
            &mut weights,
            &mut dweights,
            &mut d2weights,
        );
        if !ok {
            return Err("BRepBlend_AppFuncRoot::D2: NotDone".into());
        }
        Ok(SweepD2 {
            poles,
            dpoles,
            d2poles,
            poles2d,
            dpoles2d,
            d2poles2d,
            weights,
            dweights,
            d2weights,
        })
    }

    /// `Nb2dCurves()` (`cxx:176-182`).
    fn nb_2d_curves(&self) -> i32 {
        let (_nb_poles, _nb_knots, _degree, nb_poles2d) = self.shape();
        nb_poles2d
    }

    /// `SectionShape(NbPoles, NbKnots, Degree)` (`cxx:184-189`).
    fn section_shape(&self) -> (i32, i32, i32) {
        let (nb_poles, nb_knots, degree, _nb_poles2d) = self.shape();
        (nb_poles, nb_knots, degree)
    }

    /// `Knots(TKnots)` (`cxx:191-195`).
    fn knots(&self) -> Vec<f64> {
        let (_nb_poles, nb_knots, _degree, _nb_poles2d) = self.shape();
        let mut tknots = vec![0.0; nb_knots.max(0) as usize];
        self.func.borrow_mut().knots(&mut tknots);
        tknots
    }

    /// `Mults(TMults)` (`cxx:197-201`).
    fn mults(&self) -> Vec<i32> {
        let (_nb_poles, nb_knots, _degree, _nb_poles2d) = self.shape();
        let mut tmults = vec![0; nb_knots.max(0) as usize];
        self.func.borrow_mut().mults(&mut tmults);
        tmults
    }

    /// `IsRational()` (`cxx:203-207`).
    fn is_rational(&self) -> bool {
        self.func.borrow().is_rational()
    }

    /// `NbIntervals(S)` (`cxx:209-213`).
    fn nb_intervals(&self, s: u8) -> i32 {
        self.func.borrow().nb_intervals(shape_from_u8(s)) as i32
    }

    /// `Intervals(T, S)` (`cxx:215-219`): OCCT sizes `T` as
    /// `NbIntervals(S) + 1` (`Approx_SweepApproximation.cxx:233-240`).
    fn intervals(&self, s: u8) -> Vec<f64> {
        let shape = shape_from_u8(s);
        let nb = self.func.borrow().nb_intervals(shape);
        let mut t = vec![0.0; nb + 1];
        self.func.borrow().intervals(&mut t, shape);
        t
    }

    /// `SetInterval(First, Last)` (`cxx:221-225`).
    fn set_interval(&self, first: f64, last: f64) {
        self.func.borrow_mut().set_interval(first, last);
    }

    /// `Resolution(Index, Tol, TolU, TolV)` (`cxx:227-233`).
    fn resolution(&self, index: i32, tol: f64) -> Result<(f64, f64), String> {
        let (mut tol_u, mut tol_v) = (0.0, 0.0);
        self.func
            .borrow()
            .resolution(index, tol, &mut tol_u, &mut tol_v);
        Ok((tol_u, tol_v))
    }

    /// `GetTolerance(BoundTol, SurfTol, AngleTol, Tol3d)` (`cxx:235-250`).
    fn get_tolerance(&self, bound_tol: f64, surf_tol: f64, angle_tol: f64) -> Vec<f64> {
        let (nb_poles, _nb_knots, _degree, _nb_poles2d) = self.shape();
        let nb = nb_poles.max(0) as usize;
        let mut v3d = MathVector::new(1, nb);
        let mut v1d = MathVector::new(1, nb);
        self.func
            .borrow()
            .get_tolerance_3d(bound_tol, surf_tol, angle_tol, &mut v3d, &mut v1d);
        (1..=nb).map(|ii| v3d.value(ii)).collect()
    }

    /// `SetTolerance(Tol3d, Tol2d)` (`cxx:252-266`): recompute `myTolerance`
    /// from `Func.GetTolerance`, clamp it to `Tol2d` and push it to the solver.
    fn set_tolerance(&self, tol3d: f64, tol2d: f64) {
        let mut tolerance = self.tolerance.borrow_mut();
        self.func.borrow().get_tolerance(&mut tolerance, tol3d);
        for ii in 1..=self.dim {
            if tolerance.value(ii) > tol2d {
                tolerance.set_value(ii, tol2d);
            }
        }
        self.solver.borrow_mut().set_tolerance(&tolerance);
    }

    /// `BarycentreOfSurf()` (`cxx:268-271`).
    fn barycentre_of_surf(&self) -> Result<GpPnt, String> {
        Ok(self.bary)
    }

    /// `MaximalSection()` (`cxx:273-277`).
    fn maximal_section(&self) -> Result<f64, String> {
        Ok(self.func.borrow().get_section_size())
    }

    /// `GetMinimalWeight(Weigths)` (`cxx:279-282`).
    fn get_minimal_weight(&self) -> Result<Vec<f64>, String> {
        let (nb_poles, _nb_knots, _degree, _nb_poles2d) = self.shape();
        let mut weights = vec![0.0; nb_poles.max(0) as usize];
        self.func.borrow().get_minimal_weight(&mut weights);
        Ok(weights)
    }
}
