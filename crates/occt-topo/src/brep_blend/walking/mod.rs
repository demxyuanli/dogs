//! `BRepBlend_Walking` — the walking driver that builds one section line
//! between two surfaces, step by step along a guide line.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_Walking.hxx:41-186`
//! and `BRepBlend_Walking.cxx:95-2769` (OCCT 8.0.0).
//!
//! The line is a sequence of `Blend_Point`s in the `(U1, V1, U2, V2)` space of
//! the two supports; `Perform` solves the first section and `InternalPerform`
//! (`cxx:1966-2354`) marches along `hguide` until a bound, a restriction
//! (a pcurve of the domain) or a stop criterion is met.
//!
//! Translation notes:
//! - `line` is a `handle` in OCCT and `Perform` tests `line.IsNull()` to know
//!   whether the first and last extremities must be built (`cxx:161-169`); the
//!   port keeps `Option<Line>` for that null state and `Line()` answers `Err`
//!   where OCCT raises `StdFail_NotDone` (`hxx:176-183`).
//! - `domain1` / `domain2` / `recdomain1` / `recdomain2` are four handles that
//!   alias by default (`cxx:120-124`): `SetDomainsToRecadre` only rebinds the
//!   two `rec*` ones. `Rc<RefCell<TopolTool>>` keeps that aliasing, which is
//!   observable because `Classify` and the `Adaptor2d_Curve2d` iterator are
//!   stateful.
//! - `math_FunctionSetRoot(F, Tol, 30)` only sets `Tol` (`math_FunctionSetRoot.cxx:709-738`),
//!   so every solver is built from the tolerance vector the caller just filled.
//! - `TestArret` keeps OCCT's argument order and defaults: `TestDeflection`,
//!   `TestSolution` and `TestLengthStep` are `true`, `true`, `false`
//!   (`hxx:158-162`).

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::Arc;

use occt_core::gp::{GpPnt2d, GpVec};
use occt_geom::Surface;
use occt_math::{
    MathFunctionSetRoot, MathFunctionSetWithDerivatives, MathMatrix, MathVector,
};

use crate::blend::{BlendFunction, BlendFuncInv, BlendPoint, BlendStatus};
use crate::brep_blend::{Extremity, Line};
use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::chfi_ds::ElSpine;
use crate::fclass2d::FaceState;

mod continu;
mod extremity;
mod first_section;
mod internal_perform;
mod recadre;
mod test_arret;

/// `class BRepBlend_Walking` (`hxx:41-186`).
pub struct Walking {
    /// `previousP` (`hxx:129`): the last accepted section.
    previous_p: BlendPoint,
    /// `line` (`hxx:130`).
    line: Option<Line>,
    /// `sol` (`hxx:131`): `(U1, V1, U2, V2)` of the current section.
    sol: MathVector,
    /// `jalons` (`hxx:132`): singular points known before walking.
    jalons: Vec<BlendPoint>,
    /// `surf1` / `surf2` (`hxx:133-134`).
    surf1: Arc<dyn Surface>,
    surf2: Arc<dyn Surface>,
    /// `domain1` / `domain2` (`hxx:135-136`).
    domain1: Rc<RefCell<TopolTool>>,
    domain2: Rc<RefCell<TopolTool>>,
    /// `recdomain1` / `recdomain2` (`hxx:137-138`): the domains used for
    /// recadrage, `domain1` / `domain2` until `SetDomainsToRecadre`.
    recdomain1: Rc<RefCell<TopolTool>>,
    recdomain2: Rc<RefCell<TopolTool>>,
    /// `hguide` (`hxx:139`).
    hguide: Arc<ElSpine>,
    /// `ToCorrectOnRst1` / `ToCorrectOnRst2` (`hxx:140-141`).
    to_correct_on_rst1: bool,
    to_correct_on_rst2: bool,
    /// `CorrectedParam` (`hxx:142`).
    corrected_param: f64,
    /// `tolpoint3d` (`hxx:143`).
    tolpoint3d: f64,
    /// `tolgui` (`hxx:144`).
    tolgui: f64,
    /// `pasmax` (`hxx:145`).
    pasmax: f64,
    /// `fleche` (`hxx:146`).
    fleche: f64,
    /// `param` (`hxx:147`).
    param: f64,
    /// `sens` (`hxx:148`).
    sens: f64,
    /// `done` (`hxx:149`).
    done: bool,
    /// `rebrou` (`hxx:150`).
    rebrou: bool,
    /// `iscomplete` (`hxx:151`).
    iscomplete: bool,
    /// `comptra` (`hxx:152`).
    comptra: bool,
    /// `clasonS1` / `clasonS2` (`hxx:153-154`).
    clasons1: bool,
    clasons2: bool,
    /// `check2d` (`hxx:155`).
    check2d: bool,
    /// `check` (`hxx:156`).
    check: bool,
    /// `twistflag1` / `twistflag2` (`hxx:157-158`).
    twistflag1: bool,
    twistflag2: bool,
}

impl Walking {
    /// `BRepBlend_Walking(Surf1, Surf2, Domain1, Domain2, HGuide)`
    /// (`cxx:103-126`). `sol` is sized `(1, 4)` (`cxx:108`); the members OCCT
    /// leaves to their default constructor start at zero.
    pub fn new(
        surf1: Arc<dyn Surface>,
        surf2: Arc<dyn Surface>,
        domain1: Rc<RefCell<TopolTool>>,
        domain2: Rc<RefCell<TopolTool>>,
        hguide: Arc<ElSpine>,
    ) -> Self {
        Self {
            previous_p: BlendPoint::new(),
            line: None,
            sol: MathVector::with_init(1, 4, 0.0),
            jalons: Vec::new(),
            surf1,
            surf2,
            recdomain1: Rc::clone(&domain1),
            recdomain2: Rc::clone(&domain2),
            domain1,
            domain2,
            hguide,
            to_correct_on_rst1: false,
            to_correct_on_rst2: false,
            corrected_param: 0.0,
            tolpoint3d: 0.0,
            tolgui: 0.0,
            pasmax: 0.0,
            fleche: 0.0,
            param: 0.0,
            sens: 0.0,
            done: false,
            rebrou: false,
            iscomplete: false,
            comptra: false,
            clasons1: true,
            clasons2: true,
            check2d: true,
            check: true,
            twistflag1: false,
            twistflag2: false,
        }
    }

    /// `SetDomainsToRecadre(Domain1, Domain2)` (`cxx:128-136`).
    pub fn set_domains_to_recadre(
        &mut self,
        domain1: Rc<RefCell<TopolTool>>,
        domain2: Rc<RefCell<TopolTool>>,
    ) {
        self.recdomain1 = domain1;
        self.recdomain2 = domain2;
    }

    /// `AddSingularPoint(P)` (`cxx:138-159`): insert `p` keeping `jalons`
    /// sorted by parameter.
    pub fn add_singular_point(&mut self, p: BlendPoint) {
        if self.jalons.is_empty() {
            self.jalons.push(p);
            return;
        }
        let tp = p.parameter();
        let mut jj = 1usize;
        let mut ii = 1usize;
        let mut ti = self.jalons[0].parameter();
        while ii <= self.jalons.len() && tp > ti {
            jj = ii;
            ti = self.jalons[jj - 1].parameter();
            ii += 1;
        }
        if tp > ti {
            // `InsertAfter(jj)`: the 1-based position `jj` becomes `jj` in the
            // 0-based sequence.
            self.jalons.insert(jj, p);
        } else {
            self.jalons.insert(jj - 1, p);
        }
    }

    /// `Perform(F, FInv, Pdep, Pmax, MaxStep, Tol3d, TolGuide, Soldep, Fleche,
    /// Appro)` (`cxx:161-337`).
    #[allow(clippy::too_many_arguments)]
    pub fn perform(
        &mut self,
        func: &mut dyn BlendFunction,
        func_inv: &mut dyn BlendFuncInv,
        pdep: f64,
        pmax: f64,
        max_step: f64,
        tol3d: f64,
        tol_guide: f64,
        sol_dep: &MathVector,
        fleche: f64,
        approx: bool,
    ) {
        self.done = false;
        self.iscomplete = false;
        self.comptra = false;
        let do_extremities = match self.line.as_mut() {
            Some(line) => {
                line.clear();
                false
            }
            None => {
                self.line = Some(Line::new());
                true
            }
        };
        self.tolpoint3d = tol3d;
        self.tolgui = tol_guide.abs();
        self.fleche = fleche.abs();
        self.rebrou = false;
        self.pasmax = max_step.abs();
        self.sens = if pmax - pdep >= 0.0 { 1.0 } else { -1.0 };

        self.param = pdep;
        func.set(self.param);

        if approx {
            let mut tolerance = MathVector::with_init(1, 4, 0.0);
            let mut infbound = MathVector::with_init(1, 4, 0.0);
            let mut supbound = MathVector::with_init(1, 4, 0.0);
            func.get_tolerance(&mut tolerance, self.tolpoint3d);
            func.get_bounds(&mut infbound, &mut supbound);
            let Some(root) = solve_root_set(func, sol_dep, &tolerance, &infbound, &supbound, 30)
            else {
                return;
            };
            self.sol = root;

            let situ1 = if self.clasons1 {
                classify_domain(
                    &self.domain1,
                    self.sol.value(1),
                    self.sol.value(2),
                    tolerance.value(1).min(tolerance.value(2)),
                    false,
                )
            } else {
                FaceState::In
            };
            let situ2 = if self.clasons2 {
                classify_domain(
                    &self.domain2,
                    self.sol.value(3),
                    self.sol.value(4),
                    tolerance.value(3).min(tolerance.value(4)),
                    false,
                )
            } else {
                FaceState::In
            };
            if situ1 != FaceState::In || situ2 != FaceState::In {
                return;
            }
        } else {
            self.sol = sol_dep.clone();
        }

        let state = self.test_arret(func, BlendStatus::Ok, false, true, false);
        if state != BlendStatus::Ok {
            return;
        }

        if self.to_correct_on_rst1 || self.to_correct_on_rst2 {
            self.previous_p.set_parameter(self.corrected_param);
        }
        let appended = self.previous_p;
        if let Some(line) = self.line.as_mut() {
            line.append(appended);
        }

        if do_extremities {
            let mut ptf1 = Extremity::on_curve(
                *self.previous_p.point_on_s1(),
                self.sol.value(1),
                self.sol.value(2),
                self.tolpoint3d,
            );
            let mut ptf2 = Extremity::on_curve(
                *self.previous_p.point_on_s2(),
                self.sol.value(3),
                self.sol.value(4),
                self.tolpoint3d,
            );
            if !self.previous_p.is_tangency_point() {
                if let Some(tg) = self.previous_tangent(true) {
                    ptf1.set_tangent(tg);
                }
                if let Some(tg) = self.previous_tangent(false) {
                    ptf2.set_tangent(tg);
                }
            }
            if let Some(line) = self.line.as_mut() {
                if self.sens > 0.0 {
                    line.set_start_points(ptf1, ptf2);
                } else {
                    line.set_end_points(ptf1, ptf2);
                }
            }
        }

        self.internal_perform(func, func_inv, pmax);

        self.done = true;
    }

    /// `ClassificationOnS1(C)` (`cxx:757-760`).
    pub fn classification_on_s1(&mut self, c: bool) {
        self.clasons1 = c;
    }

    /// `ClassificationOnS2(C)` (`cxx:762-765`).
    pub fn classification_on_s2(&mut self, c: bool) {
        self.clasons2 = c;
    }

    /// `Check2d(C)` (`cxx:767-770`).
    pub fn check2d(&mut self, c: bool) {
        self.check2d = c;
    }

    /// `Check(C)` (`cxx:772-775`).
    pub fn check(&mut self, c: bool) {
        self.check = c;
    }

    /// `TwistOnS1()` (`hxx:151`).
    pub fn twist_on_s1(&self) -> bool {
        self.twistflag1
    }

    /// `TwistOnS2()` (`hxx:153`).
    pub fn twist_on_s2(&self) -> bool {
        self.twistflag2
    }

    /// `IsDone()` (`hxx:155`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `Line()` (`hxx:157-166`): `StdFail_NotDone` outside `Done`.
    pub fn line(&self) -> Result<&Line, String> {
        if !self.done {
            return Err("BRepBlend_Walking::Line: StdFail_NotDone".to_string());
        }
        self.line
            .as_ref()
            .ok_or_else(|| "BRepBlend_Walking::Line: StdFail_NotDone".to_string())
    }

    /// `previousP.TangentOnS1()` / `TangentOnS2()`: `None` for a tangency
    /// point, which is the only case where OCCT's accessor would raise.
    pub(super) fn previous_tangent(&self, on_first: bool) -> Option<GpVec> {
        if self.previous_p.is_tangency_point() {
            return None;
        }
        if on_first {
            self.previous_p.tangent_on_s1().ok().copied()
        } else {
            self.previous_p.tangent_on_s2().ok().copied()
        }
    }

    /// The 1-based arc `Index` of `tool`, after the `Init()` / `Next()` walk
    /// OCCT writes at every use site (`cxx:489-495`, `:1770-1776`).
    pub(super) fn domain_arc(tool: &Rc<RefCell<TopolTool>>, index: usize) -> Option<crate::brep_adaptor_curve2d::Curve2dAdaptor> {
        let mut iter = tool.borrow_mut();
        iter.init();
        let mut nbarc = 1usize;
        while nbarc < index {
            if !iter.more() {
                return None;
            }
            nbarc += 1;
            iter.next();
        }
        iter.value().cloned()
    }

    /// `DomainOfRst->Value()` without an `Init()`: the current arc.
    pub(super) fn current_arc(tool: &Rc<RefCell<TopolTool>>) -> Option<crate::brep_adaptor_curve2d::Curve2dAdaptor> {
        tool.borrow().value().cloned()
    }
}

/// The `Adaptor2d_Curve2d` view of a domain arc, as `BRepBlend_Curve2dTool`
/// needs it: the arc's own parameter range, re-expressed as a trimmed 2d curve
/// for the [`BlendTool`] queries, which take a bare curve.
pub(super) fn arc_curve(
    arc: &crate::brep_adaptor_curve2d::Curve2dAdaptor,
) -> Option<Arc<dyn occt_geom2d::Curve2d>> {
    let curve = Arc::clone(arc.curve()?);
    let (first, last) = (arc.first_parameter(), arc.last_parameter());
    if first == curve.first_parameter() && last == curve.last_parameter() {
        return Some(curve);
    }
    Some(Arc::new(occt_geom2d::trimmed::Geom2dTrimmedCurve::new(
        curve, first, last,
    )))
}

/// `Adaptor3d_TopolTool::Classify(P, Tol, RecadreOnPeriodic)`
/// (`Adaptor3d_TopolTool.hxx:120-122`); the two-argument overload of OCCT
/// defaults the last one to `true`.
pub(super) fn classify_domain(
    tool: &Rc<RefCell<TopolTool>>,
    u: f64,
    v: f64,
    tol: f64,
    recadre_on_periodic: bool,
) -> FaceState {
    tool.borrow_mut()
        .classify(GpPnt2d::new(u, v), tol, recadre_on_periodic)
}

/// The `MathFunctionSetWithDerivatives` view of a `&mut dyn` trait object, so
/// that `math_FunctionSetRoot` (a generic over a sized function) can drive the
/// `Blend_Function` / `Blend_FuncInv` handles OCCT passes by reference.
struct RootFunction<'a, T: ?Sized> {
    inner: &'a mut T,
}

impl<T: MathFunctionSetWithDerivatives + ?Sized> MathFunctionSetWithDerivatives
    for RootFunction<'_, T>
{
    fn nb_variables(&self) -> usize {
        self.inner.nb_variables()
    }

    fn nb_equations(&self) -> usize {
        self.inner.nb_equations()
    }

    fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        self.inner.value(x, f)
    }

    fn derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        self.inner.derivatives(x, d)
    }

    fn values(&mut self, x: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        self.inner.values(x, f, d)
    }

    fn state_number(&mut self) -> i32 {
        self.inner.state_number()
    }
}

/// `math_FunctionSetRoot rsnld(F, Tol, Iter); rsnld.Perform(F, Start, Inf, Sup);
/// rsnld.Root(Sol)` (`math_FunctionSetRoot.cxx:709-738`, `:796-1417`): the
/// root, or `None` when OCCT's `IsDone()` is false.
pub(super) fn solve_root_set<T: MathFunctionSetWithDerivatives + ?Sized>(
    function: &mut T,
    start: &MathVector,
    tol: &MathVector,
    inf_bound: &MathVector,
    sup_bound: &MathVector,
    nb_iterations: i32,
) -> Option<MathVector> {
    let mut adapter = RootFunction { inner: function };
    let mut solver = MathFunctionSetRoot::new(&adapter, tol, nb_iterations);
    solver.perform_with_bounds(&mut adapter, start, inf_bound, sup_bound, false);
    if solver.is_done() {
        Some(solver.root())
    } else {
        None
    }
}
