//! `AppBlend_AppSurf` — the template body that `GeomFill_AppSurf` and
//! `GeomFill_AppSweep` instantiate, plus the `AppBlend_Approx` accessors.
//!
//! Source: `TKGeomAlgo/AppBlend/AppBlend_AppSurf.gxx:1-1134`,
//! `AppBlend_AppSurf.lxx:20-190` and `AppBlend_Approx.hxx:29-111` (OCCT 8.0.0);
//! instantiated as
//! - `GeomFill_AppSurf` over `GeomFill_SectionGenerator`
//!   (`GeomFill_AppSurf_0.cxx:22-28`, [`AppSurf`]),
//! - `GeomFill_AppSweep` over `GeomFill_SweepSectionGenerator`
//!   (`GeomFill_AppSweep_0.cxx:22-28`, [`crate::geom_fill::AppSweep`]).
//!
//! Translation notes:
//! - The template parameter `TheSectionGenerator` is the
//!   [`AppBlendSectionGenerator`] trait. `AppBlend_AppSurf` itself carries no
//!   generator member, so the port is one struct and the generator type is a
//!   parameter of `perform`.
//! - `Handle(GeomFill_Line)` is `&Line`: `Lin.IsNull()` (`gxx:193-196`,
//!   `:602-605`) cannot occur for a value type.
//! - `NCollection_Array2<gp_Pnt>` (1..NbUPoles, 1..NbVPoles) is
//!   `Vec<Vec<Gp_Pnt>>` indexed `[u - 1][v - 1]`; `ColLength()` is the U count
//!   and `RowLength()` the V count (`NCollection_Array2.hxx:196-201`).
//! - `NCollection_Sequence<Handle(HArray1(gp_Pnt2d)))>` (`seqPoles2d`) is
//!   `Vec<Vec<GpPnt2d>>`.
//! - `Math_Vector` (`Param`, `theParams`) is [`occt_math::MathVector`].
//! - `UseSmoothing` (`gxx:472-541`) drives `AppDef_Variational`
//!   ([`Variational`](crate::app_def::variational::Variational)), ported in
//!   `app_def/variational/`; `PerformSmoothing` (`gxx:180-183`) is the
//!   `UseSmoothing == true` entry point.

mod debug;
mod perform;
mod row_scratch;

pub use debug::{
    app_blend_get_context_approx_with_no_tgt, app_blend_get_context_spline_approx,
    app_blend_set_context_approx_with_no_tgt, app_blend_set_context_spline_approx,
};

use occt_core::approx::ApproxParametrizationType;
use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_core::kernel::geomabs::Shape;

use crate::geom_fill::section_generator::SectionGenerator;
use crate::geom_fill::sweep_section_generator::SweepSectionGenerator;

/// The section a generator hands back for one approximation row: the 3d poles,
/// the 2d poles, the weights and the three derivative columns.
///
/// OCCT writes these into caller-sized `NCollection_Array1` out-parameters
/// (`NCollection_Array1<gp_Pnt>& Poles`, `... DPoles`, `... Poles2d`,
/// `... DPoles2d`, `... Weigths`, `... DWeigths`); the port returns them and
/// the caller checks the lengths against `GetShape`.
#[derive(Debug, Clone, Default)]
pub struct SectionValues {
    /// `Poles`.
    pub poles: Vec<GpPnt>,
    /// `DPoles`.
    pub d_poles: Vec<GpVec>,
    /// `Poles2d`.
    pub poles2d: Vec<GpPnt2d>,
    /// `DPoles2d`.
    pub d_poles2d: Vec<GpVec2d>,
    /// `Weigths`.
    pub weights: Vec<f64>,
    /// `DWeigths`.
    pub d_weights: Vec<f64>,
}

/// The `TheSectionGenerator` template parameter of `AppBlend_AppSurf`
/// (`GeomFill_AppSurf.hxx:22-25`, `GeomFill_AppSweep.hxx:22-25`).
pub trait AppBlendSectionGenerator {
    /// `GetShape(NbPoles, NbKnots, Degree, NbPoles2d)`; returns the four
    /// out-parameters in that order.
    fn get_shape(&self) -> Result<(i32, i32, i32, i32), String>;

    /// `Knots(TKnots)`.
    fn knots(&self) -> Result<Vec<f64>, String>;

    /// `Mults(TMults)`.
    fn mults(&self) -> Result<Vec<i32>, String>;

    /// `Section(P, Poles, DPoles, Poles2d, DPoles2d, Weigths, DWeigths)`: the
    /// `bool` result is OCCT's `withDeriv`.
    fn section_with_derivatives(&self, p: i32) -> Result<(SectionValues, bool), String>;

    /// `Section(P, Poles, Poles2d, Weigths)`: the `SpApprox` overload, which
    /// writes no derivatives.
    fn section(&self, p: i32) -> Result<SectionValues, String>;

    /// `Parameter(P)`.
    fn parameter(&self, p: i32) -> Result<f64, String>;
}

impl AppBlendSectionGenerator for SectionGenerator {
    fn get_shape(&self) -> Result<(i32, i32, i32, i32), String> {
        SectionGenerator::get_shape(self)
    }

    fn knots(&self) -> Result<Vec<f64>, String> {
        SectionGenerator::knots(self)
    }

    fn mults(&self) -> Result<Vec<i32>, String> {
        SectionGenerator::mults(self)
    }

    /// `GeomFill_SectionGenerator::Section(P, Poles, DPoles, Poles2d,
    /// DPoles2d, Weigths, DWeigths)` (`cxx:77-90`): forwards to the four
    /// argument overload and answers `false`.
    fn section_with_derivatives(&self, p: i32) -> Result<(SectionValues, bool), String> {
        let (poles, weights) = SectionGenerator::section(self, p)?;
        Ok((
            SectionValues {
                poles,
                weights,
                ..SectionValues::default()
            },
            false,
        ))
    }

    /// `GeomFill_SectionGenerator::Section(P, Poles, Poles2d, Weigths)`
    /// (`cxx:92-102`).
    fn section(&self, p: i32) -> Result<SectionValues, String> {
        let (poles, weights) = SectionGenerator::section(self, p)?;
        Ok(SectionValues {
            poles,
            weights,
            ..SectionValues::default()
        })
    }

    fn parameter(&self, p: i32) -> Result<f64, String> {
        SectionGenerator::parameter(self, p)
    }
}

impl AppBlendSectionGenerator for SweepSectionGenerator {
    fn get_shape(&self) -> Result<(i32, i32, i32, i32), String> {
        let (mut nb_poles, mut nb_knots, mut degree, mut nb_poles2d) = (0, 0, 0, 0);
        SweepSectionGenerator::get_shape(self, &mut nb_poles, &mut nb_knots, &mut degree, &mut nb_poles2d);
        Ok((nb_poles, nb_knots, degree, nb_poles2d))
    }

    fn knots(&self) -> Result<Vec<f64>, String> {
        Ok(SweepSectionGenerator::knots(self))
    }

    fn mults(&self) -> Result<Vec<i32>, String> {
        Ok(SweepSectionGenerator::mults(self))
    }

    /// `GeomFill_SweepSectionGenerator::Section(P, Poles, DPoles, Poles2d,
    /// DPoles2d, Weigths, DWeigths)` (`cxx:453-537`).
    fn section_with_derivatives(&self, p: i32) -> Result<(SectionValues, bool), String> {
        let nb_poles = self.first_section().map(|s| s.nb_poles()).unwrap_or(0);
        let mut poles = vec![GpPnt::new(0.0, 0.0, 0.0); nb_poles];
        let mut d_poles = vec![GpVec::zero(); nb_poles];
        let mut weights = vec![0.0; nb_poles];
        let mut d_weights = vec![0.0; nb_poles];
        let with_deriv = SweepSectionGenerator::section(
            self,
            p,
            &mut poles,
            &mut d_poles,
            &mut weights,
            &mut d_weights,
        );
        Ok((
            SectionValues {
                poles,
                d_poles,
                weights,
                d_weights,
                ..SectionValues::default()
            },
            with_deriv,
        ))
    }

    /// `GeomFill_SweepSectionGenerator::Section(P, Poles, Poles2d, Weigths)`
    /// (`cxx:542-665`).
    fn section(&self, p: i32) -> Result<SectionValues, String> {
        let nb_poles = self.first_section().map(|s| s.nb_poles()).unwrap_or(0);
        let mut poles = vec![GpPnt::new(0.0, 0.0, 0.0); nb_poles];
        let mut weights = vec![0.0; nb_poles];
        SweepSectionGenerator::section_simple(self, p, &mut poles, &mut weights);
        Ok(SectionValues {
            poles,
            weights,
            ..SectionValues::default()
        })
    }

    /// `GeomFill_SweepSectionGenerator::Parameter(P)` (`cxx:680-698`); OCCT
    /// dereferences a null handle for an out-of-range `P`, the port reports it.
    fn parameter(&self, p: i32) -> Result<f64, String> {
        SweepSectionGenerator::parameter(self, p)
            .ok_or_else(|| "Standard_OutOfRange: GeomFill_SweepSectionGenerator::Parameter".to_string())
    }
}

/// `AppBlend_AppSurf` (`AppBlend_AppSurf.gxx:47-1134`), instantiated by
/// [`AppSurf`] / [`crate::geom_fill::AppSweep`].
#[derive(Debug, Clone)]
pub struct AppSurf {
    /// `done` (`gxx:47`).
    done: bool,
    /// `dmin` (`gxx:48`).
    dmin: i32,
    /// `dmax` (`gxx:49`).
    dmax: i32,
    /// `tol3d` (`gxx:50`).
    tol3d: f64,
    /// `tol2d` (`gxx:51`).
    tol2d: f64,
    /// `nbit` (`gxx:52`).
    nbit: i32,
    /// `udeg` (`gxx:53`).
    udeg: i32,
    /// `vdeg` (`gxx:54`).
    vdeg: i32,
    /// `knownp` (`gxx:55`).
    knownp: bool,
    /// `tabPoles` (`gxx:56`): `[u - 1][v - 1]`.
    tab_poles: Vec<Vec<GpPnt>>,
    /// `tabWeights` (`gxx:57`).
    tab_weights: Vec<Vec<f64>>,
    /// `tabUKnots` (`gxx:58`).
    tab_u_knots: Vec<f64>,
    /// `tabVKnots` (`gxx:59`).
    tab_v_knots: Vec<f64>,
    /// `tabUMults` (`gxx:60`).
    tab_u_mults: Vec<i32>,
    /// `tabVMults` (`gxx:61`).
    tab_v_mults: Vec<i32>,
    /// `seqPoles2d` (`gxx:62`).
    seq_poles2d: Vec<Vec<GpPnt2d>>,
    /// `tol3dreached` (`gxx:63`).
    tol3d_reached: f64,
    /// `tol2dreached` (`gxx:64`).
    tol2d_reached: f64,
    /// `paramtype` (`gxx:65`).
    paramtype: ApproxParametrizationType,
    /// `continuity` (`gxx:66`).
    continuity: Shape,
    /// `critweights[3]` (`gxx:67`).
    critweights: [f64; 3],
}

impl Default for AppSurf {
    /// `AppBlend_AppSurf::AppBlend_AppSurf()` (`gxx:47-64`).
    fn default() -> Self {
        Self::new()
    }
}

impl AppSurf {
    /// `AppBlend_AppSurf::AppBlend_AppSurf()` (`gxx:47-64`).
    pub fn new() -> Self {
        Self {
            done: false,
            dmin: 0,
            dmax: 0,
            tol3d: 0.0,
            tol2d: 0.0,
            nbit: 0,
            udeg: 0,
            vdeg: 0,
            knownp: false,
            tab_poles: Vec::new(),
            tab_weights: Vec::new(),
            tab_u_knots: Vec::new(),
            tab_v_knots: Vec::new(),
            tab_u_mults: Vec::new(),
            tab_v_mults: Vec::new(),
            seq_poles2d: Vec::new(),
            tol3d_reached: 0.0,
            tol2d_reached: 0.0,
            paramtype: ApproxParametrizationType::ChordLength,
            continuity: Shape::C2,
            critweights: [0.4, 0.2, 0.4],
        }
    }

    /// `AppBlend_AppSurf(Degmin, Degmax, Tol3d, Tol2d, NbIt, KnownParameters)`
    /// (`gxx:67-92`).
    pub fn with_settings(
        degmin: i32,
        degmax: i32,
        tol3d: f64,
        tol2d: f64,
        nb_it: i32,
        known_parameters: bool,
    ) -> Self {
        Self {
            done: false,
            dmin: degmin,
            dmax: degmax,
            tol3d,
            tol2d,
            nbit: nb_it,
            udeg: 0,
            vdeg: 0,
            knownp: known_parameters,
            tab_poles: Vec::new(),
            tab_weights: Vec::new(),
            tab_u_knots: Vec::new(),
            tab_v_knots: Vec::new(),
            tab_u_mults: Vec::new(),
            tab_v_mults: Vec::new(),
            seq_poles2d: Vec::new(),
            tol3d_reached: 0.0,
            tol2d_reached: 0.0,
            paramtype: ApproxParametrizationType::ChordLength,
            continuity: Shape::C2,
            critweights: [0.4, 0.2, 0.4],
        }
    }

    /// `Init(Degmin, Degmax, Tol3d, Tol2d, NbIt, KnownParameters)`
    /// (`gxx:96-113`).
    pub fn init(
        &mut self,
        degmin: i32,
        degmax: i32,
        tol3d: f64,
        tol2d: f64,
        nb_it: i32,
        known_parameters: bool,
    ) {
        self.done = false;
        self.dmin = degmin;
        self.dmax = degmax;
        self.tol3d = tol3d;
        self.tol2d = tol2d;
        self.nbit = nb_it;
        self.knownp = known_parameters;
        self.continuity = Shape::C2;
        self.paramtype = ApproxParametrizationType::ChordLength;
        self.critweights = [0.4, 0.2, 0.4];
    }

    /// `SetParType(ParType)` (`gxx:156-159`).
    pub fn set_par_type(&mut self, par_type: ApproxParametrizationType) {
        self.paramtype = par_type;
    }

    /// `SetContinuity(C)` (`gxx:142-145`).
    pub fn set_continuity(&mut self, c: Shape) {
        self.continuity = c;
    }

    /// `SetCriteriumWeight(W1, W2, W3)` (`gxx:132-138`): `Standard_DomainError`
    /// on a negative weight.
    pub fn set_criterium_weight(&mut self, w1: f64, w2: f64, w3: f64) -> Result<(), String> {
        if w1 < 0.0 || w2 < 0.0 || w3 < 0.0 {
            return Err("Standard_DomainError: AppBlend_AppSurf::SetCriteriumWeight".into());
        }
        self.critweights = [w1, w2, w3];
        Ok(())
    }

    /// `ParType()` (`gxx:163-166`).
    pub fn par_type(&self) -> ApproxParametrizationType {
        self.paramtype
    }

    /// `Continuity()` (`gxx:149-152`).
    pub fn continuity(&self) -> Shape {
        self.continuity
    }

    /// `CriteriumWeight(W1, W2, W3)` (`gxx:123-128`).
    pub fn criterium_weight(&self) -> (f64, f64, f64) {
        (self.critweights[0], self.critweights[1], self.critweights[2])
    }

    /// `IsDone()` (`lxx:20-23`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `SurfShape(UDegree, VDegree, NbUPoles, NbVPoles, NbUKnots, NbVKnots)`
    /// (`gxx:1053-1070`).
    pub fn surf_shape(&self) -> (i32, i32, i32, i32, i32, i32) {
        self.require_done();
        (
            self.udeg,
            self.vdeg,
            self.col_length(),
            self.row_length(),
            self.tab_u_knots.len() as i32,
            self.tab_v_knots.len() as i32,
        )
    }

    /// `Surface(TPoles, TWeights, TUKnots, TVKnots, TUMults, TVMults)`
    /// (`gxx:1073-1090`): one accessor per out-parameter, in the same order.
    pub fn surface(
        &self,
    ) -> (
        &[Vec<GpPnt>],
        &[Vec<f64>],
        &[f64],
        &[f64],
        &[i32],
        &[i32],
    ) {
        self.require_done();
        (
            &self.tab_poles,
            &self.tab_weights,
            &self.tab_u_knots,
            &self.tab_v_knots,
            &self.tab_u_mults,
            &self.tab_v_mults,
        )
    }

    /// `UDegree()` (`lxx:25-32`).
    pub fn u_degree(&self) -> i32 {
        self.require_done();
        self.udeg
    }

    /// `VDegree()` (`lxx:34-41`).
    pub fn v_degree(&self) -> i32 {
        self.require_done();
        self.vdeg
    }

    /// `SurfPoles()` (`lxx:43-50`).
    pub fn surf_poles(&self) -> &[Vec<GpPnt>] {
        self.require_done();
        &self.tab_poles
    }

    /// `SurfWeights()` (`lxx:52-59`).
    pub fn surf_weights(&self) -> &[Vec<f64>] {
        self.require_done();
        &self.tab_weights
    }

    /// `SurfUKnots()` (`lxx:61-68`).
    pub fn surf_u_knots(&self) -> &[f64] {
        self.require_done();
        &self.tab_u_knots
    }

    /// `SurfVKnots()` (`lxx:70-77`).
    pub fn surf_v_knots(&self) -> &[f64] {
        self.require_done();
        &self.tab_v_knots
    }

    /// `SurfUMults()` (`lxx:79-86`).
    pub fn surf_u_mults(&self) -> &[i32] {
        self.require_done();
        &self.tab_u_mults
    }

    /// `SurfVMults()` (`lxx:88-95`).
    pub fn surf_v_mults(&self) -> &[i32] {
        self.require_done();
        &self.tab_v_mults
    }

    /// `NbCurves2d()` (`lxx:97-104`).
    pub fn nb_curves_2d(&self) -> i32 {
        self.require_done();
        self.seq_poles2d.len() as i32
    }

    /// `Curves2dShape(Degree, NbPoles, NbKnots)` (`gxx:1095-1107`).
    ///
    /// `NbPoles` is `tabPoles->ColLength()`, i.e. the **U** pole count, exactly
    /// as OCCT writes it.
    pub fn curves_2d_shape(&self) -> (i32, i32, i32) {
        self.require_done();
        if self.seq_poles2d.is_empty() {
            self.domain_error("AppBlend_AppSurf::Curves2dShape");
        }
        (
            self.vdeg,
            self.col_length(),
            self.tab_v_knots.len() as i32,
        )
    }

    /// `Curve2d(Index, TPoles, TKnots, TMults)` (`gxx:1113-1127`): the three
    /// out-parameters of one 2d curve, in order. `Index` is 1-based.
    pub fn curve_2d(&self, index: i32) -> (&[GpPnt2d], &[f64], &[i32]) {
        self.require_done();
        if self.seq_poles2d.is_empty() {
            self.domain_error("AppBlend_AppSurf::Curve2d");
        }
        if index < 1 || index as usize > self.seq_poles2d.len() {
            self.out_of_range("AppBlend_AppSurf::Curve2d");
        }
        (
            &self.seq_poles2d[(index - 1) as usize],
            &self.tab_v_knots,
            &self.tab_v_mults,
        )
    }

    /// `Curves2dDegree()` (`lxx:106-116`).
    pub fn curves_2d_degree(&self) -> i32 {
        self.require_done();
        if self.seq_poles2d.is_empty() {
            self.domain_error("AppBlend_AppSurf::Curves2dDegree");
        }
        self.vdeg
    }

    /// `Curve2dPoles(Index)` (`lxx:118-129`).
    pub fn curve_2d_poles(&self, index: i32) -> &[GpPnt2d] {
        self.require_done();
        if self.seq_poles2d.is_empty() {
            self.domain_error("AppBlend_AppSurf::Curve2dPoles");
        }
        if index < 1 || index as usize > self.seq_poles2d.len() {
            self.out_of_range("AppBlend_AppSurf::Curve2dPoles");
        }
        &self.seq_poles2d[(index - 1) as usize]
    }

    /// `Curves2dKnots()` (`lxx:131-142`).
    pub fn curves_2d_knots(&self) -> &[f64] {
        self.require_done();
        if self.seq_poles2d.is_empty() {
            self.domain_error("AppBlend_AppSurf::Curves2dKnots");
        }
        &self.tab_v_knots
    }

    /// `Curves2dMults()` (`lxx:144-155`).
    pub fn curves_2d_mults(&self) -> &[i32] {
        self.require_done();
        if self.seq_poles2d.is_empty() {
            self.domain_error("AppBlend_AppSurf::Curves2dMults");
        }
        &self.tab_v_mults
    }

    /// `TolReached(Tol3d, Tol2d)` (`lxx:157-165`).
    pub fn tol_reached(&self) -> (f64, f64) {
        self.require_done();
        (self.tol3d_reached, self.tol2d_reached)
    }

    /// `TolCurveOnSurf(Index)` (`gxx:1132-1136`): OCCT answers
    /// `tol3dreached` for every curve.
    pub fn tol_curve_on_surf(&self, _index: i32) -> f64 {
        self.tol3d_reached
    }

    /// `tabPoles->ColLength()` (`NCollection_Array2.hxx:201`): the U count.
    fn col_length(&self) -> i32 {
        self.tab_poles.len() as i32
    }

    /// `tabPoles->RowLength()` (`NCollection_Array2.hxx:198`): the V count.
    fn row_length(&self) -> i32 {
        self.tab_poles.first().map(|c| c.len()).unwrap_or(0) as i32
    }

    /// `throw StdFail_NotDone()` of the accessors.
    fn require_done(&self) {
        if !self.done {
            panic!("StdFail_NotDone: AppBlend_AppSurf");
        }
    }

    /// `throw Standard_DomainError()`.
    fn domain_error(&self, who: &str) -> ! {
        panic!("Standard_DomainError: {who}")
    }

    /// `Standard_OutOfRange` raised by `NCollection_Sequence::Value`.
    fn out_of_range(&self, who: &str) -> ! {
        panic!("Standard_OutOfRange: {who}")
    }
}

/// `GeomFill_AppSweep` (`GeomFill_AppSweep.hxx:37-59`,
/// `GeomFill_AppSweep_0.cxx:22-28`): `AppBlend_AppSurf.gxx` with
/// `TheSectionGenerator = GeomFill_SweepSectionGenerator`.
///
/// `AppBlend_AppSurf.gxx` declares no member that depends on
/// `TheSectionGenerator` and no method stores one; `GeomFill_AppSweep.hxx`
/// only renames the include macros, so the two instantiations have identical
/// members and methods and the port shares [`AppSurf`], selecting the generator
/// through the [`AppBlendSectionGenerator`] parameter of `perform`.
pub type AppSweep = AppSurf;
