//! `BRepBlend_RstRstLineBuilder` — the walking driver that builds one section
//! line between two restriction curves, each lying on its own surface.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_RstRstLineBuilder.hxx:64-205`
//! and `BRepBlend_RstRstLineBuilder.cxx:151-1964` (OCCT 8.0.0).
//!
//! The line is a sequence of `Blend_Point`s in the `(w1, w2)` space of the two
//! restrictions: `sol(1)` is the parameter on `rst1`, `sol(2)` on `rst2`, and
//! the guide parameter is the member `param`. `Perform` solves the first
//! section and marches to `Pmax`; `Complete` restarts from the first point of
//! the line and marches the other way.
//!
//! Three stop criteria coexist, as the OCCT header comment says (`hxx:47-57`):
//! leaving the domain of a restriction, losing the contact with a surface, and
//! stalling inside the surface that supports a restriction. Unlike
//! [`super::surf_rst_line_builder`] the three `Blend_Decroch*` flags are
//! reported per side, so `decrochdeb` / `decrochfin` are a
//! [`BlendDecrochStatus`] and not a `bool`.
//!
//! Translation notes:
//! - `rst1` / `rst2` are `handle<Adaptor2d_Curve2d>` in OCCT and are used both
//!   as plain curves (`FirstParameter`, `Value`, `D1`) and as curve/edge pairs
//!   by `BRepBlend_BlendTool::Parameter` / `Tolerance` and by
//!   `Adaptor3d_TopolTool::Orientation` (`cxx:1857`, `:1875`).
//!   [`Curve2dAdaptor`] carries both, so that handle is the field type.
//! - `line` is a handle allocated by `Perform` / `PerformFirstSection`
//!   (`cxx:199`, `:305`); the port keeps `Option<Line>` and `line()` answers
//!   `Err` where OCCT raises `StdFail_NotDone` (`lxx:29-37`).
//! - `domain1` / `domain2` are handles whose iterators and classifiers are
//!   stateful, so `Rc<RefCell<TopolTool>>` keeps the OCCT aliasing.
//! - `sol` is `math_Vector(1, 2)` (`cxx:159`), not `(1, 3)` as in the
//!   surface/restriction builder.

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::Arc;

use occt_core::gp::{GpPnt2d, GpVec};
use occt_geom::Surface;
use occt_math::MathVector;

use crate::blend::{
    BlendCurvPointFuncInv, BlendDecrochStatus, BlendPoint, BlendRstRstFunction, BlendStatus,
    BlendSurfCurvFuncInv,
};
use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_blend::{Extremity, Line};
use crate::brep_top_adaptor_topol_tool::TopolTool;

mod extremity;
mod first_section;
mod internal_perform;
mod recadre;
mod test_arret;

/// `class BRepBlend_RstRstLineBuilder` (`hxx:64-205`).
pub struct RstRstLineBuilder {
    /// `done` (`hxx:183`).
    done: bool,
    /// `line` (`hxx:184`): null until `Perform` / `PerformFirstSection`.
    line: Option<Line>,
    /// `sol` (`hxx:185`): `math_Vector(1, 2)` (`cxx:159`).
    sol: MathVector,
    /// `surf1` (`hxx:186`).
    surf1: Arc<dyn Surface>,
    /// `domain1` (`hxx:187`).
    domain1: Rc<RefCell<TopolTool>>,
    /// `surf2` (`hxx:188`).
    surf2: Arc<dyn Surface>,
    /// `domain2` (`hxx:189`).
    domain2: Rc<RefCell<TopolTool>>,
    /// `rst1` (`hxx:190`).
    rst1: Curve2dAdaptor,
    /// `rst2` (`hxx:191`).
    rst2: Curve2dAdaptor,
    /// `tolpoint3d` (`hxx:192`).
    tolpoint3d: f64,
    /// `tolgui` (`hxx:193`).
    tolgui: f64,
    /// `pasmax` (`hxx:194`).
    pasmax: f64,
    /// `fleche` (`hxx:195`).
    fleche: f64,
    /// `param` (`hxx:196`).
    param: f64,
    /// `previousP` (`hxx:197`).
    previous_p: BlendPoint,
    /// `rebrou` (`hxx:198`).
    rebrou: bool,
    /// `iscomplete` (`hxx:199`).
    iscomplete: bool,
    /// `comptra` (`hxx:200`).
    comptra: bool,
    /// `sens` (`hxx:201`).
    sens: f64,
    /// `decrochdeb` (`hxx:202`).
    decrochdeb: BlendDecrochStatus,
    /// `decrochfin` (`hxx:203`).
    decrochfin: BlendDecrochStatus,
}

impl RstRstLineBuilder {
    /// `BRepBlend_RstRstLineBuilder(Surf1, Rst1, Domain1, Surf2, Rst2, Domain2)`
    /// (`cxx:151-178`): only `done`, `surf1`, `domain1`, `surf2`, `domain2`,
    /// `rst1` and `rst2` are initialised there, `sol` is allocated `(1, 2)` and
    /// both `DecrochStatus` start at `Blend_NoDecroch`.
    pub fn new(
        surf1: Arc<dyn Surface>,
        rst1: Curve2dAdaptor,
        domain1: Rc<RefCell<TopolTool>>,
        surf2: Arc<dyn Surface>,
        rst2: Curve2dAdaptor,
        domain2: Rc<RefCell<TopolTool>>,
    ) -> Self {
        Self {
            done: false,
            line: None,
            sol: MathVector::new(1, 2),
            surf1,
            domain1,
            surf2,
            domain2,
            rst1,
            rst2,
            tolpoint3d: 0.0,
            tolgui: 0.0,
            pasmax: 0.0,
            fleche: 0.0,
            param: 0.0,
            previous_p: BlendPoint::new(),
            rebrou: false,
            iscomplete: false,
            comptra: false,
            sens: 0.0,
            decrochdeb: BlendDecrochStatus::NoDecroch,
            decrochfin: BlendDecrochStatus::NoDecroch,
        }
    }

    /// `Perform(Func, Finv1, FinvP1, Finv2, FinvP2, Pdep, Pmax, MaxStep, Tol3d,
    /// TolGuide, ParDep, Fleche, Appro)` (`cxx:182-281`).
    #[allow(clippy::too_many_arguments)]
    pub fn perform(
        &mut self,
        func: &mut dyn BlendRstRstFunction,
        finv1: &mut dyn BlendSurfCurvFuncInv,
        finv_p1: &mut dyn BlendCurvPointFuncInv,
        finv2: &mut dyn BlendSurfCurvFuncInv,
        finv_p2: &mut dyn BlendCurvPointFuncInv,
        pdep: f64,
        pmax: f64,
        max_step: f64,
        tol3d: f64,
        tol_guide: f64,
        par_dep: &MathVector,
        fleche: f64,
        approx: bool,
    ) {
        self.done = false;
        self.iscomplete = false;
        self.comptra = false;
        self.line = Some(Line::new());
        self.tolpoint3d = tol3d;
        self.tolgui = tol_guide.abs();
        self.fleche = fleche.abs();
        self.rebrou = false;
        self.pasmax = max_step.abs();

        self.sens = if pmax - pdep >= 0.0 { 1.0 } else { -1.0 };

        self.param = pdep;
        func.set(self.param);

        if approx {
            let mut situ_rst1 = crate::fclass2d::FaceState::Unknown;
            let mut situ_rst2 = crate::fclass2d::FaceState::Unknown;
            let mut decroch = BlendDecrochStatus::NoDecroch;
            let mut tolerance = MathVector::with_init(1, 2, 0.0);
            let mut infbound = MathVector::with_init(1, 2, 0.0);
            let mut supbound = MathVector::with_init(1, 2, 0.0);
            func.get_tolerance(&mut tolerance, self.tolpoint3d);
            func.get_bounds(&mut infbound, &mut supbound);

            let Some(root) = super::walking::solve_root_set(
                func,
                par_dep,
                &tolerance,
                &infbound,
                &supbound,
                30,
            ) else {
                return;
            };
            self.sol = root;
            if !self.check_inside(func, &mut situ_rst1, &mut situ_rst2, &mut decroch) {
                return;
            }
        } else {
            self.sol = par_dep.clone();
        }

        let state = self.test_arret(func, false, BlendStatus::Ok);
        if state != BlendStatus::Ok {
            return;
        }

        if let Some(line) = self.line.as_mut() {
            line.append(self.previous_p);
        }
        let u = self.previous_p.parameter_on_c1().unwrap_or(0.0);
        let v = self.previous_p.parameter_on_c2().unwrap_or(0.0);
        let param = self.previous_p.parameter();
        let mut ptf1 = Extremity::on_curve(*self.previous_p.point_on_c1(), u, param, self.tolpoint3d);
        let mut ptf2 = Extremity::on_curve(*self.previous_p.point_on_c2(), v, param, self.tolpoint3d);
        if !self.previous_p.is_tangency_point() {
            if let Ok(tg) = self.previous_p.tangent_on_c1() {
                ptf1.set_tangent(*tg);
            }
            if let Ok(tg) = self.previous_p.tangent_on_c2() {
                ptf2.set_tangent(*tg);
            }
        }
        if let Some(line) = self.line.as_mut() {
            if self.sens > 0.0 {
                line.set_start_points(ptf1, ptf2);
            } else {
                line.set_end_points(ptf1, ptf2);
            }
        }

        self.internal_perform(func, finv1, finv_p1, finv2, finv_p2, pmax);
        self.done = true;
    }

    /// `IsDone()` (`lxx:22-25`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `Line()` (`lxx:29-37`): `StdFail_NotDone` outside `Done`.
    pub fn line(&self) -> Result<&Line, String> {
        if !self.done {
            return Err("BRepBlend_RstRstLineBuilder::Line: StdFail_NotDone".to_string());
        }
        self.line
            .as_ref()
            .ok_or_else(|| "BRepBlend_RstRstLineBuilder::Line: StdFail_NotDone".to_string())
    }

    /// `Decroch1Start()` (`lxx:39-43`): `decrochdeb` is `DecrochRst1` or
    /// `DecrochBoth`.
    pub fn decroch1_start(&self) -> bool {
        self.decrochdeb == BlendDecrochStatus::DecrochRst1
            || self.decrochdeb == BlendDecrochStatus::DecrochBoth
    }

    /// `Decroch1End()` (`lxx:45-49`).
    pub fn decroch1_end(&self) -> bool {
        self.decrochfin == BlendDecrochStatus::DecrochRst1
            || self.decrochfin == BlendDecrochStatus::DecrochBoth
    }

    /// `Decroch2Start()` (`lxx:51-55`).
    pub fn decroch2_start(&self) -> bool {
        self.decrochdeb == BlendDecrochStatus::DecrochRst2
            || self.decrochdeb == BlendDecrochStatus::DecrochBoth
    }

    /// `Decroch2End()` (`lxx:57-61`).
    pub fn decroch2_end(&self) -> bool {
        self.decrochfin == BlendDecrochStatus::DecrochRst2
            || self.decrochfin == BlendDecrochStatus::DecrochBoth
    }

    /// `line->NbPoints()` with the null handle answered as `0`.
    fn line_nb_points(&self) -> usize {
        self.line.as_ref().map_or(0, |l| l.nb_points())
    }

    /// `line->Point(Index)` (`Index` is 1-based in OCCT).
    fn line_point(&self, index: usize) -> Option<BlendPoint> {
        self.line.as_ref()?.point(index).copied()
    }

    /// The `TangentOnC1` / `TangentOnC2` view used by `Transition`
    /// (`cxx:1490`, `:1502`).
    fn previous_tangent_of(&self, on_first: bool) -> GpVec {
        if self.previous_p.is_tangency_point() {
            return GpVec::zero();
        }
        let res = if on_first {
            self.previous_p.tangent_on_c1().copied()
        } else {
            self.previous_p.tangent_on_c2().copied()
        };
        res.unwrap_or_else(|_| GpVec::zero())
    }

    /// `arc != Arc` (`cxx:1564`): OCCT compares the adaptor handles, so the arc
    /// must be the very object taken from the same `TopolTool`. The loaded
    /// p-curve handle is the identity the port can compare.
    fn same_arc(a: &Curve2dAdaptor, b: &Curve2dAdaptor) -> bool {
        match (a.curve(), b.curve()) {
            (Some(x), Some(y)) => Arc::ptr_eq(x, y),
            _ => false,
        }
    }

    /// `rst1->Value(W)` for the point-based `Recadre` overloads.
    fn rst1_point(&self, w: f64) -> Option<GpPnt2d> {
        self.rst1.d0(w)
    }

    /// `rst2->Value(W)`.
    fn rst2_point(&self, w: f64) -> Option<GpPnt2d> {
        self.rst2.d0(w)
    }
}
