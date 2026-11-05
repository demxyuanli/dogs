//! `Geom2dHatch_Hatcher`: trims hatching curves by the elements of a 2D
//! region, resolves the global transitions and builds the domains.
//!
//! Sources: `ModelingAlgorithms/TKGeomAlgo/Geom2dHatch/Geom2dHatch_Hatcher.hxx:38-207`,
//! `.lxx:16-256`, `.cxx:46-1920` (OCCT 8.0.0). The domain state machine lives
//! in [`super::hatcher_domains`] (`cxx:1122-1806`), the trim/element part here.
//!
//! Translation notes:
//! - `NCollection_DataMap<int, Geom2dHatch_Hatching>` becomes a `BTreeMap`;
//!   every loop in OCCT walks `1..=myNbHatchings` by index, so the container
//!   order never reaches the algorithm (same substitution as
//!   [`super::element`]).
//! - `Geom2dAdaptor_Curve` is [`Curve2d`] behind an `Arc`; `AddElement` /
//!   `AddHatching` take that `Arc` directly.
//! - `IsDone()` (`hxx:160`) is declared in OCCT but defined nowhere and has no
//!   caller, so it is not ported.
//! - `Dump()` (`cxx:1842-1919`) is console output only; not ported.
//! - `Trim(IndH, IndE)` and `GlobalTransition` are private in OCCT
//!   (`hxx:190-195`); they stay private here.

use std::collections::BTreeMap;
use std::sync::Arc;

use occt_core::gp::GpDir;
use occt_core::intres2d::{IntRes2dPosition, IntRes2dTypeTrans};
use occt_core::precision::{ANGULAR, CONFUSION};
use occt_geom2d::curve::Curve2d;

use crate::abs::{Orientation, State};
use crate::hatch_gen::{
    Domain, ErrorStatus, IntersectionType, PointOnElement, PointOnHatching,
};
use crate::top_trans_curve_transition::CurveTransition;

use super::element::{Element, Elements};
use super::hatching::Hatching;
use super::intersector::Intersector;

/// `Geom2dHatch_Hatcher` (`hxx:38-207`).
pub struct Hatcher {
    /// `myIntersector`
    pub(super) my_intersector: Intersector,
    /// `myConfusion2d`
    pub(super) my_confusion2d: f64,
    /// `myConfusion3d`
    pub(super) my_confusion3d: f64,
    /// `myKeepPoints`
    pub(super) my_keep_points: bool,
    /// `myKeepSegments`
    pub(super) my_keep_segments: bool,
    /// `myNbElements`
    pub(super) my_nb_elements: i32,
    /// `myElements`
    pub(super) my_elements: Elements,
    /// `myNbHatchings`
    pub(super) my_nb_hatchings: i32,
    /// `myHatchings`
    pub(super) my_hatchings: BTreeMap<i32, Hatching>,
}

impl Hatcher {
    /// `Geom2dHatch_Hatcher(Intersector, Confusion2d, Confusion3d, KeepPnt,
    /// KeepSeg)` (`cxx:46-58`).
    pub fn new(
        intersector: Intersector,
        confusion2d: f64,
        confusion3d: f64,
        keep_pnt: bool,
        keep_seg: bool,
    ) -> Self {
        Self {
            my_intersector: intersector,
            my_confusion2d: confusion2d,
            my_confusion3d: confusion3d,
            my_keep_points: keep_pnt,
            my_keep_segments: keep_seg,
            my_nb_elements: 0,
            my_elements: Elements::new(),
            my_nb_hatchings: 0,
            my_hatchings: BTreeMap::new(),
        }
    }

    /// `Intersector(Intersector)` (`cxx:64-76`).
    pub fn set_intersector(&mut self, intersector: &Intersector) {
        self.my_intersector = intersector.clone();
        for ind_h in 1..=self.my_nb_hatchings {
            if let Some(hatching) = self.my_hatchings.get_mut(&ind_h) {
                hatching.clr_points();
            }
        }
    }

    /// `Intersector()` (`.lxx:32-35`).
    pub fn intersector(&self) -> &Intersector {
        &self.my_intersector
    }

    /// `ChangeIntersector()` (`.lxx:42-45`).
    pub fn change_intersector(&mut self) -> &mut Intersector {
        &mut self.my_intersector
    }

    /// `Confusion2d(Confusion)` (`cxx:83-95`).
    pub fn set_confusion2d(&mut self, confusion: f64) {
        self.my_confusion2d = confusion;
        for ind_h in 1..=self.my_nb_hatchings {
            if let Some(hatching) = self.my_hatchings.get_mut(&ind_h) {
                hatching.clr_points();
            }
        }
    }

    /// `Confusion2d()` (`.lxx:52-55`).
    pub fn confusion2d(&self) -> f64 {
        self.my_confusion2d
    }

    /// `Confusion3d(Confusion)` (`cxx:102-114`).
    pub fn set_confusion3d(&mut self, confusion: f64) {
        self.my_confusion3d = confusion;
        for ind_h in 1..=self.my_nb_hatchings {
            if let Some(hatching) = self.my_hatchings.get_mut(&ind_h) {
                hatching.clr_points();
            }
        }
    }

    /// `Confusion3d()` (`.lxx:62-65`).
    pub fn confusion3d(&self) -> f64 {
        self.my_confusion3d
    }

    /// `KeepPoints(Keep)` (`cxx:117-126`): changing the flag only clears the
    /// domains, the restriction points stay.
    pub fn set_keep_points(&mut self, keep: bool) {
        self.my_keep_points = keep;
        for ind_h in 1..=self.my_nb_hatchings {
            if let Some(hatching) = self.my_hatchings.get_mut(&ind_h) {
                hatching.clr_domains();
            }
        }
    }

    /// `KeepPoints()` (`.lxx:72-75`).
    pub fn keep_points(&self) -> bool {
        self.my_keep_points
    }

    /// `KeepSegments(Keep)` (`cxx:132-141`).
    pub fn set_keep_segments(&mut self, keep: bool) {
        self.my_keep_segments = keep;
        for ind_h in 1..=self.my_nb_hatchings {
            if let Some(hatching) = self.my_hatchings.get_mut(&ind_h) {
                hatching.clr_domains();
            }
        }
    }

    /// `KeepSegments()` (`.lxx:82-85`).
    pub fn keep_segments(&self) -> bool {
        self.my_keep_segments
    }

    /// `Clear()` (`.lxx:91-97`).
    pub fn clear(&mut self) {
        if self.my_nb_hatchings != 0 {
            self.clr_hatchings();
        }
        if self.my_nb_elements != 0 {
            self.clr_elements();
        }
    }

    /// `Element(IndE)` (`.lxx:103-109`, protected in OCCT).
    pub(super) fn element(&self, ind_e: i32) -> &Element {
        self.my_elements.find(ind_e)
    }

    /// `ElementCurve(IndE)` (`.lxx:116-123`).
    pub fn element_curve(&self, ind_e: i32) -> &Arc<dyn Curve2d> {
        self.element(ind_e).curve()
    }

    /// `AddElement(Curve, Orientation)` (`cxx:155-179`): the first free index,
    /// extending `myNbElements` when the map is full.
    pub fn add_element(&mut self, curve: Arc<dyn Curve2d>, orientation: Orientation) -> i32 {
        let mut ind_e = 1;
        while ind_e <= self.my_nb_elements && self.my_elements.is_bound(ind_e) {
            ind_e += 1;
        }
        if ind_e > self.my_nb_elements {
            self.my_nb_elements += 1;
            ind_e = self.my_nb_elements;
        }
        self.my_elements
            .bind(ind_e, Element::new(curve, orientation));
        for ind_h in 1..=self.my_nb_hatchings {
            if let Some(hatching) = self.my_hatchings.get_mut(&ind_h) {
                hatching.clr_points();
            }
        }
        ind_e
    }

    /// `RemElement(IndE)` (`cxx:186-223`): drops every point on element with
    /// that index, removes the points left without one, and clears the domains
    /// of the touched hatchings.
    pub fn rem_element(&mut self, ind_e: i32) {
        for ind_h in 1..=self.my_nb_hatchings {
            let Some(hatching) = self.my_hatchings.get_mut(&ind_h) else {
                continue;
            };
            let mut domains_to_clear = false;
            let nb_pnt_h = hatching.nb_points();
            for ipnt_h in (1..=nb_pnt_h).rev() {
                let pnt_h = hatching.change_point(ipnt_h);
                let nb_pnt_e = pnt_h.nb_points();
                for ipnt_e in (1..=nb_pnt_e).rev() {
                    if pnt_h.point(ipnt_e).index() == ind_e {
                        pnt_h.rem_point(ipnt_e);
                        domains_to_clear = true;
                    }
                }
                let no_point_left = pnt_h.nb_points() == 0;
                if no_point_left {
                    hatching.rem_point(ipnt_h);
                }
            }
            if domains_to_clear {
                hatching.clr_domains();
            }
        }
        self.my_elements.unbind(ind_e);
        if ind_e == self.my_nb_elements {
            self.my_nb_elements -= 1;
        }
    }

    /// `ClrElements()` (`cxx:232-249`).
    pub fn clr_elements(&mut self) {
        if self.my_nb_elements != 0 {
            if self.my_nb_hatchings != 0 {
                for ind_h in 1..=self.my_nb_hatchings {
                    if let Some(hatching) = self.my_hatchings.get_mut(&ind_h) {
                        hatching.clr_points();
                    }
                }
            }
            self.my_elements.clear();
            self.my_nb_elements = 0;
        }
    }

    /// `Hatching(IndH)` (`.lxx:127-134`, protected in OCCT).
    pub fn hatching(&self, ind_h: i32) -> &Hatching {
        self.my_hatchings
            .get(&ind_h)
            .expect("Geom2dHatch_Hatcher::Hatching")
    }

    /// `ChangeHatching(IndH)` (`.lxx:135-140`, protected in OCCT).
    pub fn hatching_mut(&mut self, ind_h: i32) -> &mut Hatching {
        self.my_hatchings
            .get_mut(&ind_h)
            .expect("Geom2dHatch_Hatcher::ChangeHatching")
    }

    /// `HatchingCurve(IndH)` (`.lxx:142-148`).
    pub fn hatching_curve(&self, ind_h: i32) -> &Arc<dyn Curve2d> {
        self.hatching(ind_h).curve()
    }

    /// `AddHatching(Curve)` (`cxx:262-277`).
    pub fn add_hatching(&mut self, curve: Arc<dyn Curve2d>) -> i32 {
        let mut ind_h = 1;
        while ind_h <= self.my_nb_hatchings && self.my_hatchings.contains_key(&ind_h) {
            ind_h += 1;
        }
        if ind_h > self.my_nb_hatchings {
            self.my_nb_hatchings += 1;
            ind_h = self.my_nb_hatchings;
        }
        self.my_hatchings.insert(ind_h, Hatching::new(curve));
        ind_h
    }

    /// `RemHatching(IndH)` (`cxx:286-294`).
    pub fn rem_hatching(&mut self, ind_h: i32) {
        self.my_hatchings
            .get_mut(&ind_h)
            .expect("Geom2dHatch_Hatcher::RemHatching")
            .clr_points();
        self.my_hatchings.remove(&ind_h);
        if ind_h == self.my_nb_hatchings {
            self.my_nb_hatchings -= 1;
        }
    }

    /// `ClrHatchings()` (`cxx:306-317`).
    pub fn clr_hatchings(&mut self) {
        if self.my_nb_hatchings != 0 {
            for ind_h in 1..=self.my_nb_hatchings {
                if let Some(hatching) = self.my_hatchings.get_mut(&ind_h) {
                    hatching.clr_points();
                }
            }
            self.my_hatchings.clear();
            self.my_nb_hatchings = 0;
        }
    }

    /// `NbPoints(IndH)` (`.lxx:156-164`).
    pub fn nb_points(&self, ind_h: i32) -> usize {
        self.hatching(ind_h).nb_points()
    }

    /// `Point(IndH, IndP)` (`.lxx:174-184`), 1-based.
    pub fn point(&self, ind_h: i32, ind_p: i32) -> &PointOnHatching {
        self.hatching(ind_h).point(ind_p as usize)
    }

    /// `Trim()` (`cxx:333-342`).
    pub fn trim(&mut self) {
        for ind_h in 1..=self.my_nb_hatchings {
            if self.my_hatchings.contains_key(&ind_h) {
                self.trim_ind(ind_h);
            }
        }
    }

    /// `Trim(Curve)` (`cxx:348-355`).
    pub fn trim_curve(&mut self, curve: Arc<dyn Curve2d>) -> i32 {
        let ind_h = self.add_hatching(curve);
        self.trim_ind(ind_h);
        ind_h
    }

    /// `Trim(IndH)` (`cxx:363-398`).
    pub fn trim_ind(&mut self, ind_h: i32) {
        self.my_hatchings
            .get_mut(&ind_h)
            .expect("Geom2dHatch_Hatcher::Trim")
            .clr_points();

        let mut all_ok = true;
        for ind_e in 1..=self.my_nb_elements {
            if self.my_elements.is_bound(ind_e) {
                let ok = self.trim_pair(ind_h, ind_e);
                all_ok = all_ok && ok;
            }
        }
        {
            let hatching = self
                .my_hatchings
                .get_mut(&ind_h)
                .expect("Geom2dHatch_Hatcher::Trim");
            hatching.set_trim_done(true);
            hatching.set_trim_failed(!all_ok);
        }

        if all_ok {
            let elements = &self.my_elements;
            let intersector = &self.my_intersector;
            let hatching = self
                .my_hatchings
                .get_mut(&ind_h)
                .expect("Geom2dHatch_Hatcher::Trim");
            let nb_pnt = hatching.nb_points();
            for ipnt in 1..=nb_pnt {
                let curve_h = Arc::clone(hatching.curve());
                let ok = Self::global_transition(
                    elements,
                    intersector,
                    &curve_h,
                    hatching.change_point(ipnt),
                );
                all_ok = all_ok && ok;
            }
            hatching.set_status(if all_ok {
                ErrorStatus::NoProblem
            } else {
                ErrorStatus::TransitionFailure
            });
        }
    }

    /// `Trim(IndH, IndE)` (`cxx:522-840`): intersects the hatching with one
    /// element and records the resulting points, merging the two points of an
    /// intersection segment when they are confused in 2D (same parameter within
    /// `myConfusion2d`) or in 3D (same transition types and a point distance
    /// within `myConfusion3d`).
    fn trim_pair(&mut self, ind_h: i32, ind_e: i32) -> bool {
        let hatching_curve = Arc::clone(
            self.my_hatchings
                .get(&ind_h)
                .expect("Geom2dHatch_Hatcher::Trim")
                .curve(),
        );
        let element_curve = Arc::clone(self.my_elements.find(ind_e).curve());

        self.my_intersector
            .intersect(hatching_curve.as_ref(), element_curve.as_ref());

        if !self.my_intersector.is_done() {
            return false;
        }
        if self.my_intersector.is_empty() {
            return true;
        }

        // Traitement des points d intersection.
        let nb_pnt_i = self.my_intersector.nb_points();
        for ipnt_i in 1..=nb_pnt_i {
            let pnt_i = self.my_intersector.point(ipnt_i);

            let mut pnt_e = PointOnElement::from_intersection_point(&pnt_i);
            pnt_e.set_index(ind_e);

            let mut pnt_h = PointOnHatching::from_intersection_point(&pnt_i);
            pnt_h.set_index(ind_h);
            pnt_h.add_point(&pnt_e, self.my_confusion2d);

            self.my_hatchings
                .get_mut(&ind_h)
                .expect("Geom2dHatch_Hatcher::Trim")
                .add_point(&pnt_h, self.my_confusion2d);
        }

        // Traitement des segments d intersection.
        let nb_seg = self.my_intersector.nb_segments();
        for iseg in 1..=nb_seg {
            let seg = self.my_intersector.segment(iseg);

            let first_point = seg.has_first_point();
            let last_point = seg.has_last_point();

            // Les deux points peuvent etre confondus.
            if first_point && last_point {
                let pnt1 = *seg.first_point();
                let pnt2 = *seg.last_point();

                let type_pnt1_h = pnt1.transition_of_first().transition_type();
                let type_pnt1_e = pnt1.transition_of_second().transition_type();
                let type_pnt2_h = pnt2.transition_of_first().transition_type();
                let type_pnt2_e = pnt2.transition_of_second().transition_type();

                // Confondus au regard de la precision du hatcher.
                let conf2d = (pnt1.param_on_first() - pnt2.param_on_first()).abs()
                    <= self.my_confusion2d;

                // Confondus au regard des intersections.
                let mut conf3d = false;
                if !conf2d {
                    conf3d = type_pnt1_h != IntRes2dTypeTrans::Touch
                        && type_pnt1_h != IntRes2dTypeTrans::Undecided
                        && type_pnt1_e != IntRes2dTypeTrans::Touch
                        && type_pnt1_e != IntRes2dTypeTrans::Undecided
                        && type_pnt2_h != IntRes2dTypeTrans::Touch
                        && type_pnt2_h != IntRes2dTypeTrans::Undecided
                        && type_pnt2_e != IntRes2dTypeTrans::Touch
                        && type_pnt2_e != IntRes2dTypeTrans::Undecided
                        && type_pnt1_h == type_pnt2_h
                        && type_pnt1_e == type_pnt2_e
                        && pnt1.value().distance(pnt2.value()) <= self.my_confusion3d;
                }

                if conf2d || conf3d {
                    let mut pnt_e = PointOnElement::new();
                    pnt_e.set_index(ind_e);
                    pnt_e.set_parameter(
                        (pnt1.param_on_second() + pnt2.param_on_second()) / 2.,
                    );
                    // `IntRes2d_Position` has exactly `HEAD` / `MIDDLE` / `END`,
                    // so OCCT's `default: break;` arm (`cxx:718`) is unreachable.
                    match pnt1.transition_of_second().position_on_curve() {
                        IntRes2dPosition::Head => pnt_e.set_position(Orientation::Forward),
                        IntRes2dPosition::Middle => {
                            match pnt2.transition_of_second().position_on_curve() {
                                IntRes2dPosition::Head => pnt_e.set_position(Orientation::Forward),
                                IntRes2dPosition::Middle => {
                                    pnt_e.set_position(Orientation::Internal)
                                }
                                IntRes2dPosition::End => {
                                    pnt_e.set_position(Orientation::Reversed)
                                }
                            }
                        }
                        IntRes2dPosition::End => pnt_e.set_position(Orientation::Reversed),
                    }
                    pnt_e.set_intersection_type(if pnt_e.position() == Orientation::Internal {
                        IntersectionType::True
                    } else {
                        IntersectionType::Touch
                    });
                    pnt_e.set_state_before(if type_pnt1_h == IntRes2dTypeTrans::In {
                        State::Out
                    } else {
                        State::In
                    });
                    pnt_e.set_state_after(if type_pnt2_h == IntRes2dTypeTrans::In {
                        State::Out
                    } else {
                        State::In
                    });

                    let mut pnt_h = PointOnHatching::new();
                    pnt_h.set_index(ind_h);
                    pnt_h.set_parameter((pnt1.param_on_first() + pnt2.param_on_first()) / 2.);
                    // Same three-value enum as above, `cxx:749`'s `default:` is
                    // unreachable.
                    match pnt1.transition_of_first().position_on_curve() {
                        IntRes2dPosition::Head => pnt_h.set_position(Orientation::Forward),
                        IntRes2dPosition::Middle => {
                            match pnt2.transition_of_first().position_on_curve() {
                                IntRes2dPosition::Head => pnt_h.set_position(Orientation::Forward),
                                IntRes2dPosition::Middle => {
                                    pnt_h.set_position(Orientation::Internal)
                                }
                                IntRes2dPosition::End => {
                                    pnt_h.set_position(Orientation::Reversed)
                                }
                            }
                        }
                        IntRes2dPosition::End => pnt_h.set_position(Orientation::Reversed),
                    }

                    pnt_h.add_point(&pnt_e, self.my_confusion2d);
                    self.my_hatchings
                        .get_mut(&ind_h)
                        .expect("Geom2dHatch_Hatcher::Trim")
                        .add_point(&pnt_h, self.my_confusion2d);
                    continue;
                }

                // Traitement du premier point du segment.
                if first_point {
                    let pnt_i = *seg.first_point();

                    let mut pnt_e = PointOnElement::from_intersection_point(&pnt_i);
                    pnt_e.set_index(ind_e);
                    pnt_e.set_segment_beginning(true);
                    pnt_e.set_segment_end(false);

                    let mut pnt_h = PointOnHatching::from_intersection_point(&pnt_i);
                    pnt_h.set_index(ind_h);
                    pnt_h.add_point(&pnt_e, self.my_confusion2d);

                    self.my_hatchings
                        .get_mut(&ind_h)
                        .expect("Geom2dHatch_Hatcher::Trim")
                        .add_point(&pnt_h, self.my_confusion2d);
                }

                // Traitement du deuxieme point du segment.
                if last_point {
                    let pnt_i = *seg.last_point();

                    let mut pnt_e = PointOnElement::from_intersection_point(&pnt_i);
                    pnt_e.set_index(ind_e);
                    pnt_e.set_segment_beginning(false);
                    pnt_e.set_segment_end(true);

                    let mut pnt_h = PointOnHatching::from_intersection_point(&pnt_i);
                    pnt_h.set_index(ind_h);
                    pnt_h.add_point(&pnt_e, self.my_confusion2d);

                    self.my_hatchings
                        .get_mut(&ind_h)
                        .expect("Geom2dHatch_Hatcher::Trim")
                        .add_point(&pnt_h, self.my_confusion2d);
                }
            }
        }
        true
    }

    /// `GlobalTransition(Point)` (`cxx:856-1121`): the before/after states and
    /// the segment flags of one restriction point, from the
    /// `TopTrans_CurveTransition` of the hatching curve against all the
    /// elements meeting at that point.
    ///
    /// OCCT reads the point out of the hatching it belongs to; the port takes
    /// the three pieces it touches so that the hatching stays mutably borrowed.
    fn global_transition(
        elements: &Elements,
        intersector: &Intersector,
        curve_h: &Arc<dyn Curve2d>,
        point: &mut PointOnHatching,
    ) -> bool {
        let mut segment_begin = false;
        let mut segment_end = false;

        let mut tangente2d = occt_core::gp::GpDir2d::new(1., 0.)
            .expect("Geom2dHatch_Hatcher::GlobalTransition: X dir");
        let mut normale2d = occt_core::gp::GpDir2d::new(0., 1.)
            .expect("Geom2dHatch_Hatcher::GlobalTransition: Y dir");
        let mut courbure = 0.;

        intersector.local_geometry(
            curve_h.as_ref(),
            point.parameter(),
            &mut tangente2d,
            &mut normale2d,
            &mut courbure,
        );

        let mut tangente = GpDir::new(tangente2d.x(), tangente2d.y(), 0.)
            .expect("Geom2dHatch_Hatcher::GlobalTransition: tangent");
        let mut normale = if courbure < CONFUSION {
            GpDir::new(-tangente2d.y(), tangente2d.x(), 0.)
                .expect("Geom2dHatch_Hatcher::GlobalTransition: normal")
        } else {
            GpDir::new(normale2d.x(), normale2d.y(), 0.)
                .expect("Geom2dHatch_Hatcher::GlobalTransition: normal")
        };

        let mut complex_transition = CurveTransition::new();
        complex_transition.reset(&tangente, &normale, courbure);

        let nb_pnt_e = point.nb_points();
        for ipnt_e in 1..=nb_pnt_e {
            let pnt_e = point.point(ipnt_e).clone();

            segment_begin = segment_begin || pnt_e.segment_beginning();
            segment_end = segment_end || pnt_e.segment_end();

            let element = elements.find(pnt_e.index());
            let curve_e = element.curve();

            let element_orientation = element.orientation();
            let to_reverse = element_orientation == Orientation::Reversed;

            // `cxx:906-921`: the parameter from the position on the element.
            let _from_position = match pnt_e.position() {
                Orientation::Forward => {
                    if to_reverse {
                        curve_e.last_parameter()
                    } else {
                        curve_e.first_parameter()
                    }
                }
                Orientation::Internal => pnt_e.parameter(),
                Orientation::Reversed => {
                    if to_reverse {
                        curve_e.first_parameter()
                    } else {
                        curve_e.last_parameter()
                    }
                }
                Orientation::External => pnt_e.parameter(),
            };
            // `cxx:933-936` (OCCT "TODO: finding #25"): the switch above is
            // dead, the point parameter always wins. Kept so the call order
            // stays identical.
            let param = pnt_e.parameter();

            intersector.local_geometry(
                curve_e.as_ref(),
                param,
                &mut tangente2d,
                &mut normale2d,
                &mut courbure,
            );

            // Calcul de la transition locale.
            let local_transition = if element_orientation == Orientation::Internal {
                Orientation::Internal
            } else if element_orientation == Orientation::External {
                Orientation::External
            } else if pnt_e.intersection_type() == IntersectionType::Tangent {
                if pnt_e.position() == Orientation::Internal {
                    match pnt_e.state_before() {
                        State::In => {
                            if to_reverse {
                                Orientation::External
                            } else {
                                Orientation::Internal
                            }
                        }
                        State::Out => {
                            if to_reverse {
                                Orientation::Internal
                            } else {
                                Orientation::External
                            }
                        }
                        _ => Orientation::External,
                    }
                } else {
                    match pnt_e.state_before() {
                        State::In => {
                            if to_reverse {
                                Orientation::Forward
                            } else {
                                Orientation::Reversed
                            }
                        }
                        State::Out => {
                            if to_reverse {
                                Orientation::Reversed
                            } else {
                                Orientation::Forward
                            }
                        }
                        _ => Orientation::External,
                    }
                }
            } else {
                match pnt_e.state_before() {
                    State::In => {
                        if to_reverse {
                            Orientation::Forward
                        } else {
                            Orientation::Reversed
                        }
                    }
                    State::Out => {
                        if to_reverse {
                            Orientation::Reversed
                        } else {
                            Orientation::Forward
                        }
                    }
                    _ => Orientation::External,
                }
            };

            // Orientation de la tangente au point d interference.
            let tangente_orientation = match pnt_e.position() {
                Orientation::Forward => {
                    if to_reverse {
                        Orientation::Reversed
                    } else {
                        Orientation::Forward
                    }
                }
                Orientation::Internal => Orientation::Internal,
                Orientation::Reversed => {
                    if to_reverse {
                        Orientation::Forward
                    } else {
                        Orientation::Reversed
                    }
                }
                Orientation::External => Orientation::Forward,
            };

            // Proprietes geometriques.
            if to_reverse {
                tangente = GpDir::new(-tangente2d.x(), -tangente2d.y(), 0.)
                    .expect("Geom2dHatch_Hatcher::GlobalTransition: reversed tangent");
            } else {
                tangente = GpDir::new(tangente2d.x(), tangente2d.y(), 0.)
                    .expect("Geom2dHatch_Hatcher::GlobalTransition: tangent");
            }
            normale = GpDir::new(normale2d.x(), normale2d.y(), 0.)
                .expect("Geom2dHatch_Hatcher::GlobalTransition: normal");

            complex_transition.compare(
                ANGULAR,
                &tangente,
                &normale,
                courbure,
                local_transition,
                tangente_orientation,
            );
        }

        let state_before = match complex_transition.state_before() {
            State::In => State::In,
            State::Out => State::Out,
            _ => return false,
        };
        let state_after = match complex_transition.state_after() {
            State::In => State::In,
            State::Out => State::Out,
            _ => return false,
        };

        point.set_state_before(state_before);
        point.set_state_after(state_after);
        point.set_segment_beginning(segment_begin);
        point.set_segment_end(segment_end);
        true
    }

    /// `TrimDone(IndH)` (`.lxx:192-197`).
    pub fn trim_done(&self, ind_h: i32) -> bool {
        self.hatching(ind_h).trim_done()
    }

    /// `TrimFailed(IndH)` (`.lxx:203-210`).
    pub fn trim_failed(&self, ind_h: i32) -> bool {
        self.hatching(ind_h).trim_failed()
    }

    /// `IsDone(IndH)` (`.lxx:221-228`).
    pub fn is_done_ind(&self, ind_h: i32) -> bool {
        self.hatching(ind_h).is_done()
    }

    /// `Status(IndH)` (`.lxx:232-240`).
    pub fn status(&self, ind_h: i32) -> ErrorStatus {
        self.hatching(ind_h).status()
    }

    /// `NbDomains(IndH)` (`.lxx:246-257`): raises `StdFail_NotDone` when the
    /// domains were not computed.
    pub fn nb_domains(&self, ind_h: i32) -> usize {
        let hatching = self.hatching(ind_h);
        assert!(
            hatching.is_done(),
            "Geom2dHatch_Hatcher::NbDomains"
        );
        hatching.nb_domains()
    }

    /// `Domain(IndH, IDom)` (`cxx:1819-1831`), 1-based: raises `StdFail_NotDone`
    /// when the domains were not computed.
    pub fn domain(&self, ind_h: i32, idom: i32) -> &Domain {
        let hatching = self.hatching(ind_h);
        assert!(hatching.is_done(), "Geom2dHatch_Hatcher::Domain");
        hatching.domain(idom as usize)
    }
}
