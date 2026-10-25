//! `BRepBlend_CSWalking` — the walking driver that builds one section line
//! between a surface and a curve, step by step along a guide line.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_CSWalking.hxx:46-117`
//! and `BRepBlend_CSWalking.cxx:88-1118` (OCCT 8.0.0).
//!
//! The line is a sequence of `Blend_Point`s in the `(U, V, W)` space of the
//! surface, the curve and the guide: `Perform` solves the first section and
//! `InternalPerform` (`cxx:824-1118`) marches along the guide until a bound or a
//! stop criterion is met. This is the surface/curve sibling of
//! [`super::walking`]: it has no recadrage (OCCT comments the `Recadre` call out
//! at `cxx:894` and leaves the whole `Recadre` body inside a `/* */` block at
//! `cxx:589-738`) and no second restriction, so `Blend_OnRst2` /
//! `Blend_OnRst12` never occur.
//!
//! Translation notes:
//! - `line` is a `handle` in OCCT, allocated by every `Perform` (`cxx:111`); the
//!   port keeps `Option<Line>` for the null state before the first `Perform`
//!   and `Line()` answers `Err` where OCCT raises `StdFail_NotDone`
//!   (`hxx:67-73`).
//! - `domain` is a handle whose iterator and classifier are stateful, so
//!   `Rc<RefCell<TopolTool>>` keeps the aliasing of the OCCT handle.
//! - `firstsol` is an `NCollection_HArray1<double>(1, Nbvar)` (`cxx:121`); a
//!   `MathVector` with the same 1-based bounds is the same array.
//! - `Complete` throws `StdFail_NotDone` before `Perform` (`cxx:218-221`), which
//!   the port reports as `Err`.
//! - OCCT's `if (recad) { ... }` branch of `InternalPerform` (`cxx:907-926`) and
//!   its `Blend_OnRst1` arm (`cxx:1061-1084`) are kept verbatim although
//!   unreachable in 8.0.0: `recad` is only ever `false` because the `Recadre`
//!   call is commented out (`cxx:894`). Nothing is invented in their place.

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::Arc;

use occt_core::gp::{GpVec, GpVec2d};
use occt_geom::{Curve, Surface};
use occt_math::MathVector;

use crate::blend::{BlendCsFunction, BlendPoint, BlendStatus};
use crate::brep_blend::{Extremity, Line};
use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::fclass2d::FaceState;

mod extremity;
mod internal_perform;
mod test_arret;

/// `class BRepBlend_CSWalking` (`hxx:46-117`).
pub struct CsWalking {
    /// `done` (`hxx:97`).
    done: bool,
    /// `line` (`hxx:98`): null until the first `Perform` (`cxx:111`).
    line: Option<Line>,
    /// `surf` (`hxx:100`).
    surf: Arc<dyn Surface>,
    /// `curv` (`hxx:101`).
    curv: Arc<dyn Curve>,
    /// `domain` (`hxx:102`).
    domain: Rc<RefCell<TopolTool>>,
    /// `tolpoint3d` (`hxx:103`).
    tolpoint3d: f64,
    /// `tolgui` (`hxx:104`).
    tolgui: f64,
    /// `pasmax` (`hxx:105`).
    pasmax: f64,
    /// `fleche` (`hxx:106`).
    fleche: f64,
    /// `param` (`hxx:107`).
    param: f64,
    /// `firstparam` (`hxx:108`).
    firstparam: f64,
    /// `firstsol` (`hxx:109`).
    firstsol: MathVector,
    /// `previousP` (`hxx:110`).
    previous_p: BlendPoint,
    /// `rebrou` (`hxx:112`).
    rebrou: bool,
    /// `iscomplete` (`hxx:113`).
    iscomplete: bool,
    /// `comptra` (`hxx:114`).
    comptra: bool,
    /// `sens` (`hxx:115`).
    sens: f64,
}

impl CsWalking {
    /// `BRepBlend_CSWalking(Curv, Surf, Domain)` (`cxx:88-96`). Only `done`,
    /// `surf`, `curv` and `domain` are initialised there; the rest start at
    /// their C++ defaults, and `firstsol` is allocated by `Perform`.
    pub fn new(
        curv: Arc<dyn Curve>,
        surf: Arc<dyn Surface>,
        domain: Rc<RefCell<TopolTool>>,
    ) -> Self {
        Self {
            done: false,
            line: None,
            surf,
            curv,
            domain,
            tolpoint3d: 0.0,
            tolgui: 0.0,
            pasmax: 0.0,
            fleche: 0.0,
            param: 0.0,
            firstparam: 0.0,
            firstsol: MathVector::new(1, 3),
            previous_p: BlendPoint::new(),
            rebrou: false,
            iscomplete: false,
            comptra: false,
            sens: 0.0,
        }
    }

    /// `Perform(F, Pdep, Pmax, MaxStep, Tol3d, TolGuide, Soldep, Fleche, Appro)`
    /// (`cxx:98-211`).
    #[allow(clippy::too_many_arguments)]
    pub fn perform(
        &mut self,
        func: &mut dyn BlendCsFunction,
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
        self.line = Some(Line::new());
        let nbvar = func.nb_variables();
        self.tolpoint3d = tol3d;
        self.tolgui = tol_guide.abs();
        self.fleche = fleche.abs();
        self.rebrou = false;
        self.pasmax = max_step.abs();
        self.firstsol = MathVector::with_init(1, nbvar, 0.0);

        self.sens = if pmax - pdep >= 0.0 { 1.0 } else { -1.0 };

        self.param = pdep;
        self.firstparam = pdep;
        func.set(self.param);

        let sol = if approx {
            let mut tolerance = MathVector::with_init(1, nbvar, 0.0);
            let mut infbound = MathVector::with_init(1, nbvar, 0.0);
            let mut supbound = MathVector::with_init(1, nbvar, 0.0);
            func.get_tolerance(&mut tolerance, self.tolpoint3d);
            func.get_bounds(&mut infbound, &mut supbound);
            let Some(root) = super::walking::solve_root_set(
                func,
                sol_dep,
                &tolerance,
                &infbound,
                &supbound,
                30,
            ) else {
                return;
            };

            // `domain->Classify(Func.Pnt2d(), std::min(tolerance(1), tolerance(2)))`
            // (`cxx:161`): the two-argument overload, so `RecadreOnPeriodic`
            // keeps its `true` default.
            let situ = self.domain.borrow_mut().classify(
                *func.pnt2d(),
                tolerance.value(1).min(tolerance.value(2)),
                true,
            );
            if situ != FaceState::In {
                return;
            }
            root
        } else {
            sol_dep.clone()
        };

        for i in 1..=nbvar {
            self.firstsol.set_value(i, sol.value(i));
        }

        let state = self.test_arret(func, &sol, false, BlendStatus::Ok);
        if state != BlendStatus::Ok {
            return;
        }

        // Update the line (`cxx:191-207`).
        if let Some(line) = self.line.as_mut() {
            line.append(self.previous_p);
        }
        let (u, v) = self.previous_p.parameters_on_s().unwrap_or((0.0, 0.0));
        let w = self.previous_p.parameter_on_c().unwrap_or(0.0);
        let p1 = Extremity::on_surface(
            *self.previous_p.point_on_s(),
            u,
            v,
            self.previous_p.parameter(),
            self.tolpoint3d,
        );
        let p2 = Extremity::on_curve(
            *self.previous_p.point_on_c(),
            w,
            self.previous_p.parameter(),
            self.tolpoint3d,
        );
        if let Some(line) = self.line.as_mut() {
            if self.sens > 0.0 {
                line.set_start_points(p1, p2);
            } else {
                line.set_end_points(p1, p2);
            }
        }

        self.internal_perform(func, sol, pmax);

        self.done = true;
    }

    /// `Complete(F, Pmin)` (`cxx:213-256`): walk the other way from the first
    /// section. `StdFail_NotDone` before `Perform` and after a previous
    /// `Complete` is an `Err` / `false`.
    pub fn complete(
        &mut self,
        func: &mut dyn BlendCsFunction,
        pmin: f64,
    ) -> Result<bool, String> {
        if !self.done {
            return Err("BRepBlend_CSWalking::Complete: StdFail_NotDone".to_string());
        }
        if self.iscomplete {
            return Ok(true);
        }

        self.sens = -self.sens;
        let nbvar = func.nb_variables();
        let mut sol = MathVector::with_init(1, nbvar, 0.0);
        for i in 1..=nbvar {
            sol.set_value(i, self.firstsol.value(i));
        }

        self.param = self.firstparam;

        self.internal_perform(func, sol, pmin);

        self.sens = -self.sens;
        self.iscomplete = true;
        Ok(true)
    }

    /// `IsDone()` (`hxx:65`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `Line()` (`hxx:67`): `StdFail_NotDone` outside `Done`.
    pub fn line(&self) -> Result<&Line, String> {
        if !self.done {
            return Err("BRepBlend_CSWalking::Line: StdFail_NotDone".to_string());
        }
        self.line
            .as_ref()
            .ok_or_else(|| "BRepBlend_CSWalking::Line: StdFail_NotDone".to_string())
    }

    /// `Rebrou()`: `rebrou` (`hxx:112`) is private in OCCT and never read, but
    /// `TestArret` keeps setting it (`cxx:318`, `:322`).
    pub fn rebrou(&self) -> bool {
        self.rebrou
    }

    /// `previousP.TangentOnS()`: the only raising case of the OCCT accessor is a
    /// tangency point, which `Blend_Point` records. OCCT's
    /// `CheckDeflectionOnSurf` calls it unguarded (`cxx:428`); a null tangent
    /// falls into the `prevNorme <= toler3d * toler3d` branch it already has.
    fn previous_tangent(&self) -> GpVec {
        if self.previous_p.is_tangency_point() {
            return GpVec::zero();
        }
        self.previous_p.tangent_on_s().ok().copied().unwrap_or_else(GpVec::zero)
    }

    /// `previousP.Tangent2d()`: the same degenerate answer as
    /// [`Self::previous_tangent`], and again only for a tangency point.
    fn previous_tangent2d(&self) -> GpVec2d {
        if self.previous_p.is_tangency_point() {
            return GpVec2d::new(0.0, 0.0);
        }
        self.previous_p
            .tangent2d()
            .unwrap_or_else(|_| GpVec2d::new(0.0, 0.0))
    }

    /// `previousP.ParametersOnS(U, V)`.
    fn previous_parameters_on_s(&self) -> (f64, f64) {
        self.previous_p.parameters_on_s().unwrap_or((0.0, 0.0))
    }
}
