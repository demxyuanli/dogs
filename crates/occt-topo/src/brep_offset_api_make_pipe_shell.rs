//! `BRepOffsetAPI_MakePipeShell` — the shell/solid swept along a wire spine.
//!
//! Source: `ModelingAlgorithms/TKOffset/BRepOffsetAPI/BRepOffsetAPI_MakePipeShell.cxx:30-267`
//! and `.hxx:59-271` (OCCT 8.0.0).
//!
//! Every method is a one-line forward to `BRepFill_PipeShell`, which is already
//! ported ([`crate::brep_fill_pipe_shell::PipeShell`]); this module only adds
//! the done/not-done bookkeeping of the `BRepBuilderAPI_MakeShape` base and the
//! two enum translations.
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `(Spine)` | `cxx:30-37` | [`BRepOffsetApiMakePipeShell::new`] |
//! | `SetMode(IsFrenet)` | `cxx:40-43` | [`BRepOffsetApiMakePipeShell::set_mode_frenet`] |
//! | `SetDiscreteMode()` | `cxx:47-50` | [`BRepOffsetApiMakePipeShell::set_discrete_mode`] |
//! | `SetMode(Axe)` | `cxx:54-57` | [`BRepOffsetApiMakePipeShell::set_mode_ax2`] |
//! | `SetMode(BiNormal)` | `cxx:61-64` | [`BRepOffsetApiMakePipeShell::set_mode_bi_normal`] |
//! | `SetMode(SpineSupport)` | `cxx:68-71` | [`BRepOffsetApiMakePipeShell::set_mode_spine_support`] |
//! | `SetMode(AuxiliarySpine, ...)` | `cxx:75-81` | [`BRepOffsetApiMakePipeShell::set_mode_auxiliary_spine`] |
//! | `Add(Profile, WithContact, WithCorrection)` | `cxx:84-90` | [`BRepOffsetApiMakePipeShell::add`] |
//! | `Add(Profile, Location, ...)` | `cxx:93-100` | [`BRepOffsetApiMakePipeShell::add_at`] |
//! | `SetLaw(Profile, L, ...)` | `cxx:103-110` | [`BRepOffsetApiMakePipeShell::set_law`] |
//! | `SetLaw(Profile, L, Location, ...)` | `cxx:113-121` | [`BRepOffsetApiMakePipeShell::set_law_at`] |
//! | `Delete(Profile)` | `cxx:124-127` | [`BRepOffsetApiMakePipeShell::delete`] |
//! | `IsReady()` | `cxx:131-134` | [`BRepOffsetApiMakePipeShell::is_ready`] |
//! | `GetStatus()` | `cxx:138-157` | [`BRepOffsetApiMakePipeShell::get_status`] |
//! | `SetTolerance(...)` | `cxx:160-166` | [`BRepOffsetApiMakePipeShell::set_tolerance`] |
//! | `SetMaxDegree(...)` | `cxx:169-172` | [`BRepOffsetApiMakePipeShell::set_max_degree`] |
//! | `SetMaxSegments(...)` | `cxx:176-179` | [`BRepOffsetApiMakePipeShell::set_max_segments`] |
//! | `SetForceApproxC1(...)` | `cxx:187-190` | [`BRepOffsetApiMakePipeShell::set_force_approx_c1`] |
//! | `SetTransitionMode(Mode)` | `cxx:194-198` | [`BRepOffsetApiMakePipeShell::set_transition_mode`] |
//! | `Simulate(N, R)` | `cxx:201-205` | [`BRepOffsetApiMakePipeShell::simulate`] |
//! | `Build(Range)` | `cxx:208-219` | [`BRepOffsetApiMakePipeShell::build`] |
//! | `MakeSolid()` | `cxx:225-237` | [`BRepOffsetApiMakePipeShell::make_solid`] |
//! | `FirstShape()` | `cxx:242-245` | [`BRepOffsetApiMakePipeShell::first_shape`] |
//! | `LastShape()` | `cxx:249-252` | [`BRepOffsetApiMakePipeShell::last_shape`] |
//! | `Generated(S)` | `cxx:256-260` | [`BRepOffsetApiMakePipeShell::generated`] |
//! | `ErrorOnSurface()` | `cxx:264-267` | [`BRepOffsetApiMakePipeShell::error_on_surface`] |
//! | `SetIsBuildHistory(bool)` | `hxx:245-248` | [`BRepOffsetApiMakePipeShell::set_is_build_history`] |
//! | `IsBuildHistory()` | `hxx:251` | [`BRepOffsetApiMakePipeShell::is_build_history`] |
//! | `Profiles(list)` | `hxx:254` | [`BRepOffsetApiMakePipeShell::profiles`] |
//! | `Spine()` | `hxx:257` | [`BRepOffsetApiMakePipeShell::spine`] |
//!
//! # Translation notes
//!
//! * The constructor calls the three-argument `SetTolerance()` (`cxx:35`,
//!   defaults `1.0e-4 / 1.0e-4 / 1.0e-2` from `hxx:216-218`) and
//!   `SetTransitionMode()` (`cxx:36`, default `BRepBuilderAPI_Transformed`,
//!   which the C cast turns into `BRepFill_Modified`, `hxx:205`), then
//!   `NotDone()`.
//! * `myPipe` is an `occ::handle`, so it is created in the constructor body
//!   rather than in a member list; the port holds the `PipeShell` by value and
//!   has no null state, which is equivalent because every OCCT path runs the
//!   body before any other call.
//! * `MakeSolid` throws `StdFail_NotDone` when the build is not done
//!   (`cxx:227-230`); the port answers `Err` with the same message.
//! * `Simulate` appends into the caller's list (`cxx:203`); the port returns
//!   the list and leaves the caller's vector untouched.
//! * `Profiles(list)` / `Spine()` have inline bodies in the header and no
//!   `.cxx` counterpart.

use std::sync::Arc;

use occt_core::gp::{GpAx2, GpDir};
use occt_geom::law::LawFunction;

use crate::brep_builder_api_pipe_error::PipeError;
use crate::brep_builder_api_transition_mode::TransitionMode;
use crate::brep_fill_pipe_shell::{PipeShell, TypeOfContact};
use crate::brep_fill_sweep::TransitionStyle;
use crate::brep_prim_api_make_sweep::MakeSweep;
use crate::shape::{TopoShape, Vertex, Wire};

/// `SetTolerance` defaults (`BRepOffsetAPI_MakePipeShell.hxx:216-218`).
pub const DEFAULT_TOL_3D: f64 = 1.0e-4;
/// `SetTolerance` defaults (`hxx:217`).
pub const DEFAULT_BOUND_TOL: f64 = 1.0e-4;
/// `SetTolerance` defaults (`hxx:218`).
pub const DEFAULT_TOL_ANGULAR: f64 = 1.0e-2;
/// `BRepFill_PipeShell::SetTransition` defaults (`BRepFill_PipeShell.hxx:250-252`)
/// used by the constructor's `SetTransitionMode()` call.
pub const DEFAULT_ANG_MIN: f64 = 1.0e-2;
/// `BRepFill_PipeShell::SetTransition` defaults (`BRepFill_PipeShell.hxx:252`).
pub const DEFAULT_ANG_MAX: f64 = 6.0;

/// `BRepOffsetAPI_MakePipeShell` (`BRepOffsetAPI_MakePipeShell.hxx:59-271`).
pub struct BRepOffsetApiMakePipeShell {
    /// `myPipe`.
    my_pipe: PipeShell,
    /// `myShape`.
    my_shape: Option<TopoShape>,
    /// `myGenerated`.
    my_generated: Vec<TopoShape>,
    /// `myDone` (the `BRepBuilderAPI_MakeShape` flag).
    my_done: bool,
}

impl BRepOffsetApiMakePipeShell {
    /// `BRepOffsetAPI_MakePipeShell(Spine)` (`cxx:30-37`).
    pub fn new(spine: &Wire) -> Self {
        let mut this = Self {
            my_pipe: PipeShell::new(spine),
            my_shape: None,
            my_generated: Vec::new(),
            my_done: false,
        };
        // `cxx:35`: `SetTolerance()`.
        this.set_tolerance(DEFAULT_TOL_3D, DEFAULT_BOUND_TOL, DEFAULT_TOL_ANGULAR);
        // `cxx:36`: `SetTransitionMode()`.
        this.set_transition_mode(TransitionMode::Transformed);
        // `cxx:37`: `NotDone()`.
        this.my_done = false;
        this
    }

    /// `SetMode(IsFrenet)` (`cxx:40-43`).
    pub fn set_mode_frenet(&mut self, is_frenet: bool) -> Result<(), String> {
        self.my_pipe.set_frenet(is_frenet)
    }

    /// `SetDiscreteMode()` (`cxx:47-50`).
    pub fn set_discrete_mode(&mut self) -> Result<(), String> {
        self.my_pipe.set_discrete()
    }

    /// `SetMode(Axe)` (`cxx:54-57`).
    pub fn set_mode_ax2(&mut self, axe: &GpAx2) -> Result<(), String> {
        self.my_pipe.set_ax2(axe)
    }

    /// `SetMode(BiNormal)` (`cxx:61-64`).
    pub fn set_mode_bi_normal(&mut self, bi_normal: &GpDir) -> Result<(), String> {
        self.my_pipe.set_bi_normal(bi_normal)
    }

    /// `SetMode(SpineSupport)` (`cxx:68-71`).
    pub fn set_mode_spine_support(&mut self, spine_support: &TopoShape) -> Result<bool, String> {
        self.my_pipe.set_spine_support(spine_support)
    }

    /// `SetMode(AuxiliarySpine, CurvilinearEquivalence, KeepContact)`
    /// (`cxx:75-81`).
    pub fn set_mode_auxiliary_spine(
        &mut self,
        auxiliary_spine: &Wire,
        curvilinear_equivalence: bool,
        keep_contact: TypeOfContact,
    ) -> Result<(), String> {
        self.my_pipe
            .set_auxiliary_spine(auxiliary_spine, curvilinear_equivalence, keep_contact)
    }

    /// `Add(Profile, WithContact, WithCorrection)` (`cxx:84-90`).
    pub fn add(
        &mut self,
        profile: &TopoShape,
        with_contact: bool,
        with_correction: bool,
    ) -> Result<(), String> {
        self.my_pipe.add(profile, with_contact, with_correction)
    }

    /// `Add(Profile, Location, WithContact, WithCorrection)` (`cxx:93-100`).
    pub fn add_at(
        &mut self,
        profile: &TopoShape,
        location: &Vertex,
        with_contact: bool,
        with_correction: bool,
    ) -> Result<(), String> {
        self.my_pipe
            .add_at(profile, location, with_contact, with_correction)
    }

    /// `SetLaw(Profile, L, WithContact, WithCorrection)` (`cxx:103-110`).
    pub fn set_law(
        &mut self,
        profile: &TopoShape,
        l: Arc<dyn LawFunction>,
        with_contact: bool,
        with_correction: bool,
    ) -> Result<(), String> {
        self.my_pipe.set_law(profile, l, with_contact, with_correction)
    }

    /// `SetLaw(Profile, L, Location, WithContact, WithCorrection)`
    /// (`cxx:113-121`).
    pub fn set_law_at(
        &mut self,
        profile: &TopoShape,
        l: Arc<dyn LawFunction>,
        location: &Vertex,
        with_contact: bool,
        with_correction: bool,
    ) -> Result<(), String> {
        self.my_pipe
            .set_law_at(profile, l, location, with_contact, with_correction)
    }

    /// `Delete(Profile)` (`cxx:124-127`): `myPipe->DeleteProfile(Profile)`.
    pub fn delete(&mut self, profile: &TopoShape) -> Result<(), String> {
        self.my_pipe.delete_profile(profile)
    }

    /// `IsReady()` (`cxx:131-134`).
    pub fn is_ready(&self) -> bool {
        self.my_pipe.is_ready()
    }

    /// `GetStatus()` (`cxx:138-157`): the `GeomFill_PipeError` →
    /// `BRepBuilderAPI_PipeError` switch.
    pub fn get_status(&self) -> PipeError {
        PipeError::from_geom_fill_pipe_error(self.my_pipe.get_status())
    }

    /// `SetTolerance(Tol3d, BoundTol, TolAngular)` (`cxx:160-166`).
    pub fn set_tolerance(&mut self, tol_3d: f64, bound_tol: f64, tol_angular: f64) {
        self.my_pipe.set_tolerance(tol_3d, bound_tol, tol_angular);
    }

    /// `SetMaxDegree(NewMaxDegree)` (`cxx:169-172`).
    pub fn set_max_degree(&mut self, new_max_degree: i32) {
        self.my_pipe.set_max_degree(new_max_degree);
    }

    /// `SetMaxSegments(NewMaxSegments)` (`cxx:176-179`).
    pub fn set_max_segments(&mut self, new_max_segments: i32) {
        self.my_pipe.set_max_segments(new_max_segments);
    }

    /// `SetForceApproxC1(ForceApproxC1)` (`cxx:187-190`).
    pub fn set_force_approx_c1(&mut self, force_approx_c1: bool) {
        self.my_pipe.set_force_approx_c1(force_approx_c1);
    }

    /// `SetTransitionMode(Mode)` (`cxx:194-198`): `myPipe->SetTransition`
    /// with the `AngMin` / `AngMax` defaults of
    /// `BRepFill_PipeShell.hxx:250-252`.
    pub fn set_transition_mode(&mut self, mode: TransitionMode) {
        self.my_pipe.set_transition(
            mode.to_transition_style(),
            DEFAULT_ANG_MIN,
            DEFAULT_ANG_MAX,
        );
    }

    /// `SetTransitionMode` spelled with an explicit `BRepFill_TransitionStyle`
    /// and angle window, i.e. the `BRepFill_PipeShell::SetTransition` call that
    /// `cxx:196` makes through the cast. Only reachable through the pipe shell
    /// itself in OCCT; kept here so the two-argument form above has one
    /// implementation.
    pub fn set_transition(
        &mut self,
        mode: TransitionStyle,
        ang_min: f64,
        ang_max: f64,
    ) {
        self.my_pipe.set_transition(mode, ang_min, ang_max);
    }

    /// `Simulate(N, R)` (`cxx:201-205`).
    pub fn simulate(&mut self, n: usize) -> Result<Vec<TopoShape>, String> {
        self.my_pipe.simulate(n)
    }

    /// `Build(theRange)` (`cxx:208-219`): copy the pipe's shape on success,
    /// keep the previous `myShape` and answer not done on failure.
    pub fn build(&mut self) -> bool {
        // `cxx:211`.
        let ok = self.my_pipe.build().unwrap_or(false);
        if ok {
            // `cxx:213-215`.
            self.my_shape = self.my_pipe.shape().cloned();
            self.my_done = true;
        } else {
            // `cxx:217-218`.
            self.my_done = false;
        }
        self.my_done
    }

    /// `MakeSolid()` (`cxx:225-237`).
    pub fn make_solid(&mut self) -> Result<bool, String> {
        // `cxx:227-230`.
        if !self.my_done {
            return Err("StdFail_NotDone: BRepOffsetAPI_MakePipeShell::MakeSolid".to_string());
        }
        // `cxx:232`.
        let ok = self.my_pipe.make_solid()?;
        if ok {
            // `cxx:234-236`.
            self.my_shape = self.my_pipe.shape().cloned();
        }
        Ok(ok)
    }

    /// `FirstShape()` (`cxx:242-245`).
    pub fn first_shape(&self) -> Option<TopoShape> {
        self.my_pipe.first_shape().cloned()
    }

    /// `LastShape()` (`cxx:249-252`).
    pub fn last_shape(&self) -> Option<TopoShape> {
        self.my_pipe.last_shape().cloned()
    }

    /// `Generated(S)` (`cxx:256-260`): refill `myGenerated`.
    pub fn generated(&mut self, s: &TopoShape) -> &[TopoShape] {
        self.my_generated = self.my_pipe.generated(s);
        &self.my_generated
    }

    /// `ErrorOnSurface()` (`cxx:264-267`).
    pub fn error_on_surface(&self) -> f64 {
        self.my_pipe.error_on_surface()
    }

    /// `SetIsBuildHistory(theIsBuildHistory)` (`hxx:245-248`).
    pub fn set_is_build_history(&mut self, is_build_history: bool) {
        self.my_pipe.set_is_build_history(is_build_history);
    }

    /// `IsBuildHistory()` (`hxx:251`).
    pub fn is_build_history(&self) -> bool {
        self.my_pipe.is_build_history()
    }

    /// `Profiles(theProfiles)` (`hxx:254`).
    pub fn profiles(&self) -> Vec<TopoShape> {
        self.my_pipe.profiles()
    }

    /// `Spine()` (`hxx:257`).
    pub fn spine(&self) -> &Wire {
        self.my_pipe.spine()
    }

    /// `IsDone()` — the `BRepBuilderAPI_MakeShape` flag set by `Build`.
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// `Shape()` — the `myShape` written by `Build` / `MakeSolid`.
    pub fn shape(&self) -> Option<&TopoShape> {
        self.my_shape.as_ref()
    }
}

impl MakeSweep for BRepOffsetApiMakePipeShell {
    fn first_shape(&self) -> Option<TopoShape> {
        BRepOffsetApiMakePipeShell::first_shape(self)
    }

    fn last_shape(&self) -> Option<TopoShape> {
        BRepOffsetApiMakePipeShell::last_shape(self)
    }
}
