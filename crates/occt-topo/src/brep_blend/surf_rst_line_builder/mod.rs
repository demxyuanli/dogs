//! `BRepBlend_SurfRstLineBuilder` — the walking driver that builds one section
//! line between a surface and a restriction curve lying on a second surface.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_SurfRstLineBuilder.hxx:64-176`
//! and `BRepBlend_SurfRstLineBuilder.cxx:137-1777` (OCCT 8.0.0).
//!
//! The line is a sequence of `Blend_Point`s in the `(U, V, W)` space of the
//! first surface, the second surface and the restriction: `U, V` are on
//! `surf1`, `W` is on `rst`, which is the p-curve of an edge of `surf2`.
//! `Perform` solves the first section and marches to `Pmax`; `Complete` first
//! backtracks to the first point of the line and marches the other way.
//!
//! Unlike [`super::cs_walking`] this builder does reframe: when the solution
//! leaves `domain1` or the range of `rst` it re-solves the inverse problems on
//! the neighbouring arc (`cxx:637-793`), which is what produces the
//! `Blend_OnRst1` / `Blend_OnRst2` / `Blend_OnRst12` stop states.
//!
//! Translation notes:
//! - `rst` is a `handle<Adaptor2d_Curve2d>` in OCCT and is used both as a
//!   plain curve (`FirstParameter`, `Value`, `D1`) and as a `curve/edge` pair
//!   by `BRepBlend_BlendTool::Parameter` / `Tolerance` and by
//!   `Adaptor3d_TopolTool::Orientation` (`cxx:703`, `:1698`). [`Curve2dAdaptor`]
//!   carries both, so that handle is the field type; the bare
//!   `Arc<dyn Curve2d>` view that the function interfaces need is
//!   [`super::walking::arc_curve`].
//! - `line` is a handle allocated by `Perform` / `PerformFirstSection`
//!   (`cxx:243`, `:346`); the port keeps `Option<Line>` for the null state and
//!   `line()` answers `Err` where OCCT raises `StdFail_NotDone` (`lxx:31-37`).
//! - `domain1` / `domain2` are handles whose iterators and classifiers are
//!   stateful, so `Rc<RefCell<TopolTool>>` keeps the OCCT aliasing.
//! - `ArcToRecadre` (`cxx:141-191`) declares `(theSol, PrevIndex, pt2d,
//!   lastpt2d, ponarc)` in the header but defines `(theSol, PrevIndex,
//!   lastpt2d, pt2d, ponarc)` in the `.cxx`; every caller passes the previous
//!   point third and the solution point fourth, so the port follows the
//!   definition and keeps `&mut` outputs.

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::Arc;

use occt_core::gp::{GpPnt2d, GpVec};
use occt_geom::Surface;
use occt_math::MathVector;

use crate::blend::{
    BlendFuncInv, BlendPoint, BlendStatus, BlendSurfCurvFuncInv, BlendSurfPointFuncInv,
    BlendSurfRstFunction,
};
use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_blend::{BlendTool, Extremity, Line};
use crate::brep_top_adaptor_topol_tool::TopolTool;

mod extremity;
mod first_section;
mod internal_perform;
mod recadre;
mod test_arret;

/// `class BRepBlend_SurfRstLineBuilder` (`hxx:64-176`).
pub struct SurfRstLineBuilder {
    /// `done` (`hxx:178`).
    done: bool,
    /// `line` (`hxx:179`): null until `Perform` / `PerformFirstSection`.
    line: Option<Line>,
    /// `sol` (`hxx:180`): `math_Vector(1, 3)` (`cxx:203`).
    sol: MathVector,
    /// `surf1` (`hxx:181`).
    surf1: Arc<dyn Surface>,
    /// `domain1` (`hxx:182`).
    domain1: Rc<RefCell<TopolTool>>,
    /// `surf2` (`hxx:183`).
    surf2: Arc<dyn Surface>,
    /// `rst` (`hxx:184`).
    rst: Curve2dAdaptor,
    /// `domain2` (`hxx:185`).
    domain2: Rc<RefCell<TopolTool>>,
    /// `tolpoint3d` (`hxx:186`).
    tolpoint3d: f64,
    /// `tolpoint2d` (`hxx:187`).
    tolpoint2d: f64,
    /// `tolgui` (`hxx:188`).
    tolgui: f64,
    /// `pasmax` (`hxx:189`).
    pasmax: f64,
    /// `fleche` (`hxx:190`).
    fleche: f64,
    /// `param` (`hxx:191`).
    param: f64,
    /// `previousP` (`hxx:192`).
    previous_p: BlendPoint,
    /// `rebrou` (`hxx:193`).
    rebrou: bool,
    /// `iscomplete` (`hxx:194`).
    iscomplete: bool,
    /// `comptra` (`hxx:195`).
    comptra: bool,
    /// `sens` (`hxx:196`).
    sens: f64,
    /// `decrochdeb` (`hxx:197`).
    decrochdeb: bool,
    /// `decrochfin` (`hxx:198`).
    decrochfin: bool,
}

impl SurfRstLineBuilder {
    /// `BRepBlend_SurfRstLineBuilder(Surf1, Domain1, Surf2, Rst, Domain2)`
    /// (`cxx:195-221`): only `done`, `surf1`, `domain1`, `surf2`, `rst` and
    /// `domain2` are initialised there, `sol` is allocated `(1, 3)` and the
    /// rest start at their C++ defaults.
    pub fn new(
        surf1: Arc<dyn Surface>,
        domain1: Rc<RefCell<TopolTool>>,
        surf2: Arc<dyn Surface>,
        rst: Curve2dAdaptor,
        domain2: Rc<RefCell<TopolTool>>,
    ) -> Self {
        Self {
            done: false,
            line: None,
            sol: MathVector::new(1, 3),
            surf1,
            domain1,
            surf2,
            rst,
            domain2,
            tolpoint3d: 0.0,
            tolpoint2d: 0.0,
            tolgui: 0.0,
            pasmax: 0.0,
            fleche: 0.0,
            param: 0.0,
            previous_p: BlendPoint::new(),
            rebrou: false,
            iscomplete: false,
            comptra: false,
            sens: 0.0,
            decrochdeb: false,
            decrochfin: false,
        }
    }

    /// `ArcToRecadre(theSol, PrevIndex, lastpt2d, pt2d, ponarc)`
    /// (`cxx:141-191`): the index (1-based) of the arc of `domain1` that best
    /// carries the solution point, skipping `PrevIndex`.
    ///
    /// `Inters` is only tried when the line already has a point, and it wins
    /// over the plain projection as soon as it succeeds with `PrevIndex == 0`
    /// (`cxx:182-185`).
    fn arc_to_recadre(
        &mut self,
        the_sol: &MathVector,
        prev_index: i32,
        last_pt2d: &mut GpPnt2d,
        pt2d: &mut GpPnt2d,
        pon_arc: &mut f64,
    ) -> i32 {
        let mut index_sol = 0;
        let mut nbarc = 0;
        let by_inter = self.line_nb_points() != 0;
        let mut ok_inter = false;
        let mut dist_min = f64::MAX;
        let mut prm = 0.0;
        let mut dist = 0.0;

        let (uprev, vprev) = if by_inter {
            self.previous_p.parameters_on_s().unwrap_or((0.0, 0.0))
        } else {
            (0.0, 0.0)
        };

        *pt2d = GpPnt2d::new(the_sol.value(1), the_sol.value(2));
        *last_pt2d = GpPnt2d::new(uprev, vprev);
        self.domain1.borrow_mut().init();

        while self.domain1.borrow().more() {
            nbarc += 1;
            let Some(arc) = self.domain1.borrow().value().cloned() else {
                break;
            };
            let mut ok = false;
            if by_inter {
                ok = arc.curve().is_some_and(|c| {
                    BlendTool::inters(pt2d, last_pt2d, c.as_ref(), &mut prm, &mut dist)
                });
                ok_inter = ok;
            }
            if !ok {
                ok = arc
                    .curve()
                    .is_some_and(|c| BlendTool::project(pt2d, c.as_ref(), &mut prm, &mut dist));
            }
            if ok && nbarc != prev_index && (dist < dist_min || ok_inter) {
                dist_min = dist;
                *pon_arc = prm;
                index_sol = nbarc;
                if ok_inter && prev_index == 0 {
                    break;
                }
            }
            self.domain1.borrow_mut().next();
        }
        index_sol
    }

    /// `Perform(F, Finv, FinvP, FinvC, Pdep, Pmax, MaxStep, Tol3d, Tol2d,
    /// TolGuide, ParDep, Fleche, Appro)` (`cxx:225-324`).
    #[allow(clippy::too_many_arguments)]
    pub fn perform(
        &mut self,
        func: &mut dyn BlendSurfRstFunction,
        finv: &mut dyn BlendFuncInv,
        finv_p: &mut dyn BlendSurfPointFuncInv,
        finv_c: &mut dyn BlendSurfCurvFuncInv,
        pdep: f64,
        pmax: f64,
        max_step: f64,
        tol3d: f64,
        tol2d: f64,
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
        self.tolpoint2d = tol2d;
        self.tolgui = tol_guide.abs();
        self.fleche = fleche.abs();
        self.rebrou = false;
        self.pasmax = max_step.abs();

        self.sens = if pmax - pdep >= 0.0 { 1.0 } else { -1.0 };

        self.param = pdep;
        func.set(self.param);

        if approx {
            let mut situ_on_c = crate::fclass2d::FaceState::Unknown;
            let mut situ_on_s = crate::fclass2d::FaceState::Unknown;
            let mut decroch = false;
            let mut tolerance = MathVector::with_init(1, 3, 0.0);
            let mut infbound = MathVector::with_init(1, 3, 0.0);
            let mut supbound = MathVector::with_init(1, 3, 0.0);
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
            if !self.check_inside(func, &mut situ_on_c, &mut situ_on_s, &mut decroch) {
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
        let (u, v) = self.previous_p.parameters_on_s().unwrap_or((0.0, 0.0));
        let param = self.previous_p.parameter();
        let mut ptf1 = Extremity::on_surface(*self.previous_p.point_on_s(), u, v, param, self.tolpoint3d);
        // `cxx:308` uses the five-argument constructor, so the point on the
        // restriction is stored with the `(U, V)` of the surface.
        let mut ptf2 = Extremity::on_surface(*self.previous_p.point_on_c(), u, v, param, self.tolpoint3d);
        if !self.previous_p.is_tangency_point() {
            if let Ok(tg) = self.previous_p.tangent_on_s() {
                ptf1.set_tangent(*tg);
            }
            if let Ok(tg) = self.previous_p.tangent_on_c() {
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

        self.internal_perform(func, finv, finv_p, finv_c, pmax);
        self.done = true;
    }

    /// `IsDone()` (`lxx:22-25`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `Line()` (`lxx:29-37`): `StdFail_NotDone` outside `Done`.
    pub fn line(&self) -> Result<&Line, String> {
        if !self.done {
            return Err("BRepBlend_SurfRstLineBuilder::Line: StdFail_NotDone".to_string());
        }
        self.line
            .as_ref()
            .ok_or_else(|| "BRepBlend_SurfRstLineBuilder::Line: StdFail_NotDone".to_string())
    }

    /// `DecrochStart()` (`lxx:39-42`).
    pub fn decroch_start(&self) -> bool {
        self.decrochdeb
    }

    /// `DecrochEnd()` (`lxx:46-49`).
    pub fn decroch_end(&self) -> bool {
        self.decrochfin
    }

    /// `line->NbPoints()` with the null handle answered as `0`.
    fn line_nb_points(&self) -> usize {
        self.line.as_ref().map_or(0, |l| l.nb_points())
    }

    /// `line->Point(Index)` (`Index` is 1-based in OCCT).
    fn line_point(&self, index: usize) -> Option<BlendPoint> {
        self.line.as_ref()?.point(index).copied()
    }

    /// The `TangentOnS1` / `TangentOnS2` view used by `Transition`
    /// (`cxx:1331`, `:1344`).
    fn previous_tangent_of(&self, on_first: bool) -> GpVec {
        if self.previous_p.is_tangency_point() {
            return GpVec::zero();
        }
        let res = if on_first {
            self.previous_p.tangent_on_s1().copied()
        } else {
            self.previous_p.tangent_on_s2().copied()
        };
        res.unwrap_or_else(|_| GpVec::zero())
    }

    /// `arc != Arc` (`cxx:1400`): OCCT compares the adaptor handles. It is a
    /// pointer comparison, so the arc must be the very object taken from the
    /// same `TopolTool` -- the `Blend_OnRst1` arm passes `domain1->Value()`
    /// (`cxx:916`), while `rst` never matches an arc rebuilt by
    /// `domain2->Initialize(rst)`. The loaded p-curve handle is the identity
    /// the port can compare.
    fn same_arc(a: &Curve2dAdaptor, b: &Curve2dAdaptor) -> bool {
        match (a.curve(), b.curve()) {
            (Some(x), Some(y)) => Arc::ptr_eq(x, y),
            _ => false,
        }
    }
}
