//! `ChFiDS_Stripe` — the data of one band of fillet between two faces.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiDS/ChFiDS_Stripe.hxx:34-169`,
//! `ChFiDS_Stripe.lxx:18-259`, `ChFiDS_Stripe.cxx:25-294` (OCCT 8.0.0).
//!
//! Divergences that are not heuristics:
//! - `occ::handle<ChFiDS_Spine>` becomes [`SpineHandle`], the shared-mutable
//!   handle shape the crate uses for its transient data,
//! - `occ::handle<NCollection_HSequence<occ::handle<ChFiDS_SurfData>>>` becomes
//!   a `Vec<Arc<SurfData>>` owned by the stripe: OCCT hands the sequence handle
//!   out (`SetOfSurfData` / `ChangeSetOfSurfData`) so callers can append, and
//!   `Reset` `Nullify`s it; the port keeps the same two accessors over the
//!   vector and `Reset` clears it,
//! - `occ::handle<Geom2d_Curve>` becomes `Option<Arc<dyn Curve2d>>`.

use std::sync::{Arc, RwLock};

use occt_geom2d::curve::Curve2d;

use crate::abs::Orientation;

use super::spine::SpineHandle;
use super::surf_data::SurfData;

/// `occ::handle<ChFiDS_Stripe>`: the stripes are shared between the builder,
/// the vertex map and the DS (`ChFiDS_StripeMap`, `ChFi3d_Builder`), and are
/// mutated through those handles.
pub type StripeHandle = Arc<RwLock<Stripe>>;

/// Source: `ChFiDS_Stripe.hxx:35-169`.
#[derive(Clone)]
pub struct Stripe {
    /// `pardeb1`
    pardeb1: f64,
    /// `parfin1`
    parfin1: f64,
    /// `pardeb2`
    pardeb2: f64,
    /// `parfin2`
    parfin2: f64,
    /// `mySpine`
    my_spine: Option<SpineHandle>,
    /// `myHdata`
    my_hdata: Vec<Arc<SurfData>>,
    /// `pcrv1`
    pcrv1: Option<Arc<dyn Curve2d>>,
    /// `pcrv2`
    pcrv2: Option<Arc<dyn Curve2d>>,
    /// `myChoix`
    my_choix: i32,
    /// `indexOfSolid`
    index_of_solid: i32,
    /// `indexOfcurve1`
    index_of_curve1: i32,
    /// `indexOfcurve2`
    index_of_curve2: i32,
    /// `indexfirstPOnS1`
    index_first_p_on_s1: i32,
    /// `indexlastPOnS1`
    index_last_p_on_s1: i32,
    /// `indexfirstPOnS2`
    index_first_p_on_s2: i32,
    /// `indexlastPOnS2`
    index_last_p_on_s2: i32,
    /// `begfilled` — the `bool` of `occ293` widened to a count by `InDS`.
    begfilled: i32,
    /// `endfilled`
    endfilled: i32,
    /// `myOr1`
    my_or1: Orientation,
    /// `myOr2`
    my_or2: Orientation,
    /// `orcurv1`
    orcurv1: Orientation,
    /// `orcurv2`
    orcurv2: Orientation,
}

/// `dyn Curve2d` carries no `Debug`, so the two pcurve handles print as
/// set / unset rather than by content (same treatment as
/// `ChFiDS_FaceInterference`).
impl std::fmt::Debug for Stripe {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Stripe")
            .field("pardeb1", &self.pardeb1)
            .field("parfin1", &self.parfin1)
            .field("pardeb2", &self.pardeb2)
            .field("parfin2", &self.parfin2)
            .field("has_spine", &self.my_spine.is_some())
            .field("nb_surf_data", &self.my_hdata.len())
            .field("has_first_p_curve", &self.pcrv1.is_some())
            .field("has_last_p_curve", &self.pcrv2.is_some())
            .field("my_choix", &self.my_choix)
            .field("index_of_solid", &self.index_of_solid)
            .field("index_of_curve1", &self.index_of_curve1)
            .field("index_of_curve2", &self.index_of_curve2)
            .field("index_first_p_on_s1", &self.index_first_p_on_s1)
            .field("index_last_p_on_s1", &self.index_last_p_on_s1)
            .field("index_first_p_on_s2", &self.index_first_p_on_s2)
            .field("index_last_p_on_s2", &self.index_last_p_on_s2)
            .field("begfilled", &self.begfilled)
            .field("endfilled", &self.endfilled)
            .field("my_or1", &self.my_or1)
            .field("my_or2", &self.my_or2)
            .field("orcurv1", &self.orcurv1)
            .field("orcurv2", &self.orcurv2)
            .finish()
    }
}

impl Default for Stripe {
    fn default() -> Self {
        Self::new()
    }
}
impl Stripe {
    /// `ChFiDS_Stripe()` (`cxx:25-45`): every scalar `0`, both faces and both
    /// pcurves `FORWARD`, `begfilled` / `endfilled` `0`.
    pub fn new() -> Self {
        Self {
            pardeb1: 0.0,
            parfin1: 0.0,
            pardeb2: 0.0,
            parfin2: 0.0,
            my_spine: None,
            my_hdata: Vec::new(),
            pcrv1: None,
            pcrv2: None,
            my_choix: 0,
            index_of_solid: 0,
            index_of_curve1: 0,
            index_of_curve2: 0,
            index_first_p_on_s1: 0,
            index_last_p_on_s1: 0,
            index_first_p_on_s2: 0,
            index_last_p_on_s2: 0,
            begfilled: 0,
            endfilled: 0,
            my_or1: Orientation::Forward,
            my_or2: Orientation::Forward,
            orcurv1: Orientation::Forward,
            orcurv2: Orientation::Forward,
        }
    }

    /// `Reset()` (`cxx:47-55`): clears the surface data and the two pcurve
    /// orientations and resets the spine.
    ///
    /// OCCT nullifies `pcrv1` twice (`cxx:52-53`) and leaves `pcrv2` alone; the
    /// port clears `pcrv1` once and keeps the same `pcrv2` behaviour.
    pub fn reset(&mut self) {
        self.my_hdata.clear();
        self.orcurv1 = Orientation::Forward;
        self.orcurv2 = Orientation::Forward;
        self.pcrv1 = None;
        if let Some(spine) = &self.my_spine {
            if let Ok(mut spine) = spine.write() {
                let _ = spine.reset(false);
            }
        }
    }

    /// `SetOfSurfData()` (`lxx:19-23`).
    pub fn set_of_surf_data(&self) -> &[Arc<SurfData>] {
        &self.my_hdata
    }

    /// `ChangeSetOfSurfData()` (`lxx:56-60`).
    pub fn change_set_of_surf_data(&mut self) -> &mut Vec<Arc<SurfData>> {
        &mut self.my_hdata
    }

    /// `Spine()` (`lxx:27-31`).
    pub fn spine(&self) -> Option<&SpineHandle> {
        self.my_spine.as_ref()
    }

    /// `ChangeSpine()` (`lxx:63-67`): OCCT hands out the handle by reference so
    /// the caller can assign a new spine; the port exposes the same slot.
    pub fn change_spine(&mut self) -> &mut Option<SpineHandle> {
        &mut self.my_spine
    }

    /// Convenience assignment used by the builders, equivalent to
    /// `Stripe->ChangeSpine() = Sp`.
    pub fn set_spine(&mut self, spine: SpineHandle) {
        self.my_spine = Some(spine);
    }

    /// `OrientationOnFace1()` (`lxx:35-38`).
    pub fn orientation_on_face1(&self) -> Orientation {
        self.my_or1
    }

    /// `OrientationOnFace2()` (`lxx:42-45`).
    pub fn orientation_on_face2(&self) -> Orientation {
        self.my_or2
    }

    /// `OrientationOnFace1(Or1)` (`lxx:70-73`).
    pub fn set_orientation_on_face1(&mut self, or1: Orientation) {
        self.my_or1 = or1;
    }

    /// `OrientationOnFace2(Or2)` (`lxx:76-79`).
    pub fn set_orientation_on_face2(&mut self, or2: Orientation) {
        self.my_or2 = or2;
    }

    /// `Choix()` (`lxx:48-51`).
    pub fn choix(&self) -> i32 {
        self.my_choix
    }

    /// `Choix(C)` (`lxx:85-88`).
    pub fn set_choix(&mut self, c: i32) {
        self.my_choix = c;
    }

    /// `FirstParameters(Pdeb, Pfin)` (`lxx:92-96`).
    pub fn first_parameters(&self) -> (f64, f64) {
        (self.pardeb1, self.parfin1)
    }

    /// `LastParameters(Pdeb, Pfin)` (`lxx:99-103`).
    pub fn last_parameters(&self) -> (f64, f64) {
        (self.pardeb2, self.parfin2)
    }

    /// `ChangeFirstParameters(Pdeb, Pfin)` (`lxx:106-111`).
    pub fn change_first_parameters(&mut self, pdeb: f64, pfin: f64) {
        self.pardeb1 = pdeb;
        self.parfin1 = pfin;
    }

    /// `ChangeLastParameters(Pdeb, Pfin)` (`lxx:114-119`).
    pub fn change_last_parameters(&mut self, pdeb: f64, pfin: f64) {
        self.pardeb2 = pdeb;
        self.parfin2 = pfin;
    }

    /// `FirstCurve()` (`lxx:123-126`).
    pub fn first_curve(&self) -> i32 {
        self.index_of_curve1
    }

    /// `LastCurve()` (`lxx:131-134`).
    pub fn last_curve(&self) -> i32 {
        self.index_of_curve2
    }

    /// `ChangeFirstCurve(Index)` (`lxx:137-140`).
    pub fn change_first_curve(&mut self, index: i32) {
        self.index_of_curve1 = index;
    }

    /// `ChangeLastCurve(Index)` (`lxx:143-146`).
    pub fn change_last_curve(&mut self, index: i32) {
        self.index_of_curve2 = index;
    }

    /// `FirstPCurve()` (`lxx:151-154`).
    pub fn first_p_curve(&self) -> Option<&Arc<dyn Curve2d>> {
        self.pcrv1.as_ref()
    }

    /// `LastPCurve()` (`lxx:158-161`).
    pub fn last_p_curve(&self) -> Option<&Arc<dyn Curve2d>> {
        self.pcrv2.as_ref()
    }

    /// `ChangeFirstPCurve()` (`lxx:165-168`).
    pub fn change_first_p_curve(&mut self) -> &mut Option<Arc<dyn Curve2d>> {
        &mut self.pcrv1
    }

    /// `ChangeLastPCurve()` (`lxx:172-175`).
    pub fn change_last_p_curve(&mut self) -> &mut Option<Arc<dyn Curve2d>> {
        &mut self.pcrv2
    }

    /// `FirstPCurveOrientation()` (`lxx:235-238`).
    pub fn first_p_curve_orientation(&self) -> Orientation {
        self.orcurv1
    }

    /// `LastPCurveOrientation()` (`lxx:242-245`).
    pub fn last_p_curve_orientation(&self) -> Orientation {
        self.orcurv2
    }

    /// `FirstPCurveOrientation(O)` (`lxx:248-251`).
    pub fn set_first_p_curve_orientation(&mut self, o: Orientation) {
        self.orcurv1 = o;
    }

    /// `LastPCurveOrientation(O)` (`lxx:254-257`).
    pub fn set_last_p_curve_orientation(&mut self, o: Orientation) {
        self.orcurv2 = o;
    }

    /// `IndexFirstPointOnS1()` (`lxx:178-181`).
    pub fn index_first_point_on_s1(&self) -> i32 {
        self.index_first_p_on_s1
    }

    /// `IndexLastPointOnS1()` (`lxx:185-188`).
    pub fn index_last_point_on_s1(&self) -> i32 {
        self.index_last_p_on_s1
    }

    /// `IndexFirstPointOnS2()` (`lxx:192-195`).
    pub fn index_first_point_on_s2(&self) -> i32 {
        self.index_first_p_on_s2
    }

    /// `IndexLastPointOnS2()` (`lxx:199-202`).
    pub fn index_last_point_on_s2(&self) -> i32 {
        self.index_last_p_on_s2
    }

    /// `ChangeIndexFirstPointOnS1(Index)` (`lxx:207-210`).
    pub fn change_index_first_point_on_s1(&mut self, index: i32) {
        self.index_first_p_on_s1 = index;
    }

    /// `ChangeIndexLastPointOnS1(Index)` (`lxx:213-216`).
    pub fn change_index_last_point_on_s1(&mut self, index: i32) {
        self.index_last_p_on_s1 = index;
    }

    /// `ChangeIndexFirstPointOnS2(Index)` (`lxx:221-224`).
    pub fn change_index_first_point_on_s2(&mut self, index: i32) {
        self.index_first_p_on_s2 = index;
    }

    /// `ChangeIndexLastPointOnS2(Index)` (`lxx:227-230`).
    pub fn change_index_last_point_on_s2(&mut self, index: i32) {
        self.index_last_p_on_s2 = index;
    }

    /// `SetCurve(Index, First)` (`cxx:87-99`).
    pub fn set_curve(&mut self, index: i32, first: bool) {
        if first {
            self.index_of_curve1 = index;
        } else {
            self.index_of_curve2 = index;
        }
    }

    /// `SetIndexPoint(Index, First, OnS)` (`cxx:232-259`).
    pub fn set_index_point(&mut self, index: i32, first: bool, on_s: i32) {
        if first {
            if on_s == 1 {
                self.index_first_p_on_s1 = index;
            } else {
                self.index_first_p_on_s2 = index;
            }
        } else if on_s == 1 {
            self.index_last_p_on_s1 = index;
        } else {
            self.index_last_p_on_s2 = index;
        }
    }

    /// `Parameters(First, Pdeb, Pfin)` (`cxx:57-69`).
    pub fn parameters(&self, first: bool) -> (f64, f64) {
        if first {
            (self.pardeb1, self.parfin1)
        } else {
            (self.pardeb2, self.parfin2)
        }
    }

    /// `SetParameters(First, Pdeb, Pfin)` (`cxx:73-85`).
    pub fn set_parameters(&mut self, first: bool, pdeb: f64, pfin: f64) {
        if first {
            self.pardeb1 = pdeb;
            self.parfin1 = pfin;
        } else {
            self.pardeb2 = pdeb;
            self.parfin2 = pfin;
        }
    }

    /// `Curve(First)` (`cxx:88-97`).
    pub fn curve(&self, first: bool) -> i32 {
        if first {
            self.index_of_curve1
        } else {
            self.index_of_curve2
        }
    }

    /// `PCurve(First)` (`cxx:120-130`).
    pub fn p_curve(&self, first: bool) -> Option<&Arc<dyn Curve2d>> {
        if first {
            self.pcrv1.as_ref()
        } else {
            self.pcrv2.as_ref()
        }
    }

    /// `ChangePCurve(First)` (`cxx:134-144`).
    pub fn change_p_curve(&mut self, first: bool) -> &mut Option<Arc<dyn Curve2d>> {
        if first {
            &mut self.pcrv1
        } else {
            &mut self.pcrv2
        }
    }

    /// `Orientation(OnS)` (`cxx:147-157`): the orientation on face `OnS`
    /// (`1` or `2`; OCCT treats anything but `1` as face 2).
    pub fn orientation_on_surface(&self, on_s: i32) -> Orientation {
        if on_s == 1 {
            self.my_or1
        } else {
            self.my_or2
        }
    }

    /// `SetOrientation(Or, OnS)` (`cxx:161-173`).
    pub fn set_orientation_on_surface(&mut self, or: Orientation, on_s: i32) {
        if on_s == 1 {
            self.my_or1 = or;
        } else {
            self.my_or2 = or;
        }
    }

    /// `Orientation(First)` (`cxx:177-187`): the pcurve orientation at the
    /// first / last end.
    pub fn orientation(&self, first: bool) -> Orientation {
        if first {
            self.orcurv1
        } else {
            self.orcurv2
        }
    }

    /// `SetOrientation(Or, First)` (`cxx:191-202`).
    pub fn set_orientation(&mut self, or: Orientation, first: bool) {
        if first {
            self.orcurv1 = or;
        } else {
            self.orcurv2 = or;
        }
    }

    /// `IndexPoint(First, OnS)` (`cxx:205-225`).
    pub fn index_point(&self, first: bool, on_s: i32) -> i32 {
        if first {
            if on_s == 1 {
                self.index_first_p_on_s1
            } else {
                self.index_first_p_on_s2
            }
        } else if on_s == 1 {
            self.index_last_p_on_s1
        } else {
            self.index_last_p_on_s2
        }
    }

    /// `SolidIndex()` (`cxx:258-261`).
    pub fn solid_index(&self) -> i32 {
        self.index_of_solid
    }

    /// `SetSolidIndex(Index)` (`cxx:263-266`).
    pub fn set_solid_index(&mut self, index: i32) {
        self.index_of_solid = index;
    }

    /// `InDS(First, Nb = 1)` (`cxx:270-281`): records how many `SurfData` of
    /// this end were put in the DS.
    pub fn in_ds(&mut self, first: bool, nb: i32) {
        if first {
            self.begfilled = nb;
        } else {
            self.endfilled = nb;
        }
    }

    /// `IsInDS(First)` (`cxx:285-293`).
    pub fn is_in_ds(&self, first: bool) -> i32 {
        if first {
            self.begfilled
        } else {
            self.endfilled
        }
    }
}
