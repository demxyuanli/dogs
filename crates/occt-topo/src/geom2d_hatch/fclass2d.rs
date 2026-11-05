//! `Geom2dHatch_FClass2dOfClassifier` and `Geom2dHatch_Classifier`: the face
//! classifier the hatcher runs when a hatching has no restriction point.
//!
//! Sources (OCCT 8.0.0):
//! - `ModelingAlgorithms/TKGeomAlgo/Geom2dHatch/Geom2dHatch_FClass2dOfClassifier.hxx:33-86`,
//!   `.cxx:26-75`.
//! - `ModelingAlgorithms/TKGeomAlgo/Geom2dHatch/Geom2dHatch_Classifier.hxx:35-83`,
//!   `.cxx:27-83`.
//! - The two templates they instantiate:
//!   `TKGeomAlgo/TopClass/TopClass_Classifier2d.pxx:26-228` (`Reset`,
//!   `Compare`) and `TKGeomAlgo/TopClass/TopClass_FaceClassifier.pxx:28-176`
//!   (`Perform`, `State`, `Edge`, `EdgeParameter`), with
//!   `TheFaceExplorer = Geom2dHatch_Elements`,
//!   `TheEdge = Geom2dAdaptor_Curve` and
//!   `FClass2d = Geom2dHatch_FClass2dOfClassifier`.
//!
//! The `BRepClass` instantiation of the same templates is ported separately in
//! [`crate::brep_class`]; the two are kept apart exactly as OCCT keeps its two
//! instantiations apart, because the edge and intersector types differ.
//!
//! Not carried over: the `TopClass_FaceClassifier::Edge` / `EdgeParameter`
//! `Standard_DomainError` is a `panic!` here, as in the `BRepClass` port.

use std::sync::Arc;

use occt_core::gp::{GpDir, GpDir2d, GpLin2d, GpPnt2d};
use occt_core::intres2d::{
    IntRes2dIntersectionPoint, IntRes2dPosition, IntRes2dSituation, IntRes2dTypeTrans,
};
use occt_core::precision::REAL_EPSILON;
use occt_geom2d::curve::Curve2d;

use crate::abs::{Orientation, State};
use crate::top_trans_curve_transition::CurveTransition;

use super::element::Elements;
use super::intersector::Intersector;

/// `Geom2dHatch_FClass2dOfClassifier` (`hxx:33-86`).
pub struct FClass2dOfClassifier {
    /// `myIsSet` — written by `Reset` (`pxx:50`), never read, as in OCCT.
    #[allow(dead_code)]
    my_is_set: bool,
    /// `myFirstCompare`
    my_first_compare: bool,
    /// `myFirstTrans`
    my_first_trans: bool,
    /// `myLin`
    my_lin: GpLin2d,
    /// `myParam`
    my_param: f64,
    /// `myTolerance`
    my_tolerance: f64,
    /// `myTrans`
    my_trans: CurveTransition,
    /// `myIntersector`
    my_intersector: Intersector,
    /// `myClosest`
    my_closest: i32,
    /// `myState`
    my_state: State,
    /// `myIsHeadOrEnd`
    my_is_head_or_end: bool,
}

impl Default for FClass2dOfClassifier {
    fn default() -> Self {
        Self::new()
    }
}

impl FClass2dOfClassifier {
    /// `Geom2dHatch_FClass2dOfClassifier()` (`cxx:28-38`): the line is the
    /// default `gp_Lin2d` (origin, X direction), the intersector is the default
    /// one (both tolerances null).
    pub fn new() -> Self {
        Self {
            my_is_set: false,
            my_first_compare: true,
            my_first_trans: true,
            my_lin: GpLin2d::default(),
            my_param: 0.0,
            my_tolerance: 0.0,
            my_trans: CurveTransition::new(),
            my_intersector: Intersector::default(),
            my_closest: 0,
            my_state: State::Unknown,
            my_is_head_or_end: false,
        }
    }

    /// `Reset(L, P, Tol)` (`cxx:42-55` -> `pxx:30-52`).
    pub fn reset(&mut self, l: &GpLin2d, p: f64, tol: f64) {
        self.my_lin = *l;
        self.my_param = p;
        self.my_tolerance = tol;
        self.my_state = State::Unknown;
        self.my_first_compare = true;
        self.my_first_trans = true;
        self.my_closest = 0;
        self.my_is_set = true;
        self.my_is_head_or_end = false;
    }

    /// `Compare(E, Or)` (`cxx:57-75` -> `pxx:55-225`).
    pub fn compare(&mut self, edge: &dyn Curve2d, orientation: Orientation) {
        // Intersect the edge and the segment.
        self.my_closest = 0;
        self.my_intersector
            .perform(&self.my_lin, self.my_param, self.my_tolerance, edge);
        if !self.my_intersector.is_done() {
            return;
        }
        if self.my_intersector.nb_points() == 0 && self.my_intersector.nb_segments() == 0 {
            return;
        }

        // Find the closest point.
        let nb_points = self.my_intersector.nb_points();
        let mut d_min = f64::MAX;
        let mut p_closest: Option<IntRes2dIntersectionPoint> = None;
        for a_point in 1..=nb_points {
            let p_inter = self.my_intersector.point(a_point);
            // Test for ON.
            if p_inter.transition_of_first().position_on_curve() == IntRes2dPosition::Head {
                self.my_closest = a_point as i32;
                self.my_state = State::On;
                return;
            }
            let param_first = p_inter.param_on_first();
            if param_first < d_min {
                self.my_closest = a_point as i32;
                p_closest = Some(p_inter);
                d_min = param_first;
            }
        }

        // For the segments only the first point is tested.
        let nb_segments = self.my_intersector.nb_segments();
        for a_segment in 1..=nb_segments {
            let p_inter = *self.my_intersector.segment(a_segment).first_point();
            if p_inter.transition_of_first().position_on_curve() == IntRes2dPosition::Head {
                self.my_closest = (nb_points + a_segment + a_segment - 1) as i32;
                self.my_state = State::On;
                return;
            }
            let param_first = p_inter.param_on_first();
            if param_first < d_min {
                self.my_closest = (nb_points + a_segment + a_segment - 1) as i32;
                p_closest = Some(p_inter);
                d_min = param_first;
            }
        }

        // If no point was found return.
        if self.my_closest == 0 {
            return;
        }

        // If the edge is INTERNAL or EXTERNAL, no problem.
        if orientation == Orientation::Internal {
            self.my_state = State::In;
            return;
        } else if orientation == Orientation::External {
            self.my_state = State::Out;
            return;
        }

        if !self.my_first_compare && d_min > self.my_param {
            return;
        }

        // Process the closest point found at `d_min` on the line.
        self.my_first_compare = false;
        if self.my_param > d_min {
            self.my_first_trans = true;
        }
        self.my_param = d_min;

        // `my_closest` is set together with `p_closest`.
        let Some(p_closest) = p_closest else {
            return;
        };
        let t2 = *p_closest.transition_of_second();
        self.my_is_head_or_end = matches!(
            t2.position_on_curve(),
            IntRes2dPosition::Head | IntRes2dPosition::End
        );

        // Transition on the segment.
        let t1 = *p_closest.transition_of_first();
        let reversed = orientation == Orientation::Reversed;
        let seg_trans = match t1.transition_type() {
            IntRes2dTypeTrans::In => {
                if reversed {
                    Orientation::Reversed
                } else {
                    Orientation::Forward
                }
            }
            IntRes2dTypeTrans::Out => {
                if reversed {
                    Orientation::Forward
                } else {
                    Orientation::Reversed
                }
            }
            IntRes2dTypeTrans::Touch => match t1.situation() {
                IntRes2dSituation::Inside => {
                    if reversed {
                        Orientation::External
                    } else {
                        Orientation::Internal
                    }
                }
                IntRes2dSituation::Outside => {
                    if reversed {
                        Orientation::Internal
                    } else {
                        Orientation::External
                    }
                }
                IntRes2dSituation::Unknown => return,
            },
            IntRes2dTypeTrans::Undecided => return,
        };

        // Are we inside the edge?
        if !self.my_is_head_or_end {
            // `aPClosest` is inside the edge.
            self.my_state = match seg_trans {
                Orientation::Forward | Orientation::External => State::Out,
                Orientation::Reversed | Orientation::Internal => State::In,
            };
            return;
        }

        // `aPClosest` is Head or End of the edge: update the complex transition.
        let mut tang2d = GpDir2d::new(1., 0.).expect("unit X direction");
        let mut norm2d = GpDir2d::new(0., 1.).expect("unit Y direction");
        let mut curv = 0.;
        self.my_intersector.local_geometry(
            edge,
            p_closest.param_on_second(),
            &mut tang2d,
            &mut norm2d,
            &mut curv,
        );
        let tang = GpDir::new(tang2d.x(), tang2d.y(), 0.).expect("unit 2d tangent");
        let norm = GpDir::new(norm2d.x(), norm2d.y(), 0.).expect("unit 2d normal");
        if self.my_first_trans {
            let dir = self.my_lin.direction();
            let tgt = GpDir::new(dir.x(), dir.y(), 0.).expect("unit line direction");
            self.my_trans.reset_line(&tgt);
            self.my_first_trans = false;
        }

        let ort = if t2.position_on_curve() == IntRes2dPosition::Head {
            Orientation::Forward
        } else {
            Orientation::Reversed
        };
        // `TopClass_Classifier2d.pxx:223` passes `RealEpsilon()`.
        self.my_trans
            .compare(REAL_EPSILON, &tang, &norm, curv, seg_trans, ort);
        self.my_state = self.my_trans.state_before();
    }

    /// `Parameter()` (`hxx:53`).
    pub fn parameter(&self) -> f64 {
        self.my_param
    }

    /// `Intersector()` (`hxx:57`): OCCT returns a non-const reference, the port
    /// hands out the same intersector for reading.
    pub fn intersector(&self) -> &Intersector {
        &self.my_intersector
    }

    /// `ClosestIntersection()` (`hxx:63`).
    pub fn closest_intersection(&self) -> i32 {
        self.my_closest
    }

    /// `State()` (`hxx:66`).
    pub fn state(&self) -> State {
        self.my_state
    }

    /// `IsHeadOrEnd()` (`hxx:71`).
    pub fn is_head_or_end(&self) -> bool {
        self.my_is_head_or_end
    }
}

/// `Geom2dHatch_Classifier` (`hxx:35-83`): the `TopClass_FaceClassifier`
/// instantiation over [`Elements`].
pub struct Classifier {
    /// `myClassifier`
    my_classifier: FClass2dOfClassifier,
    /// `myEdge` — a default-constructed `Geom2dAdaptor_Curve` in OCCT, i.e. a
    /// null curve handle, so the port stores an `Option`.
    my_edge: Option<Arc<dyn Curve2d>>,
    /// `myEdgeParameter`
    my_edge_parameter: f64,
    /// `myPosition` — not in the OCCT constructor initialiser list (`cxx:29-33`),
    /// so OCCT leaves it indeterminate until `Perform`; the port stores `None`
    /// instead of inventing a value.
    my_position: Option<IntRes2dPosition>,
    /// `rejected`
    rejected: bool,
    /// `nowires`
    nowires: bool,
}

impl Default for Classifier {
    fn default() -> Self {
        Self::new()
    }
}

impl Classifier {
    /// `Geom2dHatch_Classifier()` (`cxx:29-33`).
    pub fn new() -> Self {
        Self {
            my_classifier: FClass2dOfClassifier::new(),
            my_edge: None,
            my_edge_parameter: 0.0,
            my_position: None,
            rejected: false,
            nowires: true,
        }
    }

    /// `Geom2dHatch_Classifier(F, P, Tol)` (`cxx:37-47`).
    pub fn with_args(f: &mut Elements, p: &GpPnt2d, tol: f64) -> Self {
        let mut out = Self::new();
        out.perform(f, p, tol);
        out
    }

    /// `Perform(F, P, Tol)` (`cxx:52-62`) over
    /// `TopClass_FaceClassifier::Perform` (`pxx:29-145`).
    pub fn perform(&mut self, fexp: &mut Elements, p: &GpPnt2d, tol: f64) {
        let a_point = *p;
        while !fexp.check_point(&a_point) {}

        // Test for rejection.
        self.rejected = fexp.reject(&a_point);
        if self.rejected {
            return;
        }

        let mut a_line = GpLin2d::default();
        let mut a_param = 0.;
        let mut is_valid_segment = fexp.segment(&a_point, &mut a_line, &mut a_param);
        let mut a_state = State::Unknown;

        self.nowires = true;

        while is_valid_segment {
            self.my_classifier.reset(&a_line, a_param, tol);

            fexp.init_wires();
            while fexp.more_wires() {
                self.nowires = false;
                let is_w_reject = fexp.reject_wire(&a_line, self.my_classifier.parameter());

                if !is_w_reject {
                    // Test this wire.
                    fexp.init_edges();
                    while fexp.more_edges() {
                        let is_e_reject = fexp.reject_edge(&a_line, self.my_classifier.parameter());

                        if !is_e_reject {
                            // Test this edge.
                            let (an_edge, an_edge_ori) = fexp.current_edge();

                            if an_edge_ori == Orientation::Forward
                                || an_edge_ori == Orientation::Reversed
                            {
                                self.my_classifier.compare(an_edge.as_ref(), an_edge_ori);
                                let closest_ind = self.my_classifier.closest_intersection();

                                if closest_ind != 0 {
                                    // Save the closest edge.
                                    let nb_pnts =
                                        self.my_classifier.intersector().nb_points() as i32;
                                    self.my_edge = Some(Arc::clone(&an_edge));

                                    let p_inter = if closest_ind <= nb_pnts {
                                        self.my_classifier.intersector().point(closest_ind as usize)
                                    } else {
                                        let ci = closest_ind - nb_pnts;
                                        let seg_index = ((ci + 1) / 2) as usize;
                                        let seg = self.my_classifier.intersector().segment(seg_index);
                                        if ci & 1 != 0 {
                                            *seg.first_point()
                                        } else {
                                            *seg.last_point()
                                        }
                                    };

                                    self.my_position =
                                        Some(p_inter.transition_of_second().position_on_curve());
                                    self.my_edge_parameter = p_inter.param_on_second();
                                }
                                // If we are ON, we stop.
                                a_state = self.my_classifier.state();

                                if a_state == State::On {
                                    return;
                                }
                            }
                        }
                        fexp.next_edge();
                    }

                    // If we are out of the wire we stop.
                    a_state = self.my_classifier.state();

                    if a_state == State::Out {
                        return;
                    }
                }
                fexp.next_wire();
            }

            if !self.my_classifier.is_head_or_end() && a_state != State::Unknown {
                break;
            }

            // Bad case for classification. Trying to get another segment.
            is_valid_segment = fexp.other_segment(&a_point, &mut a_line, &mut a_param);
        }
    }

    /// `State()` (`cxx:66-69` + `pxx:147-157`).
    pub fn state(&self) -> State {
        if self.rejected {
            State::Out
        } else if self.nowires {
            State::In
        } else {
            self.my_classifier.state()
        }
    }

    /// `Rejected()` (`hxx:55-57`).
    pub fn rejected(&self) -> bool {
        self.rejected
    }

    /// `NoWires()` (`hxx:61-63`).
    pub fn no_wires(&self) -> bool {
        self.nowires
    }

    /// `Edge()` (`cxx:73-76` + `pxx:162-167`): `None` while no closest edge was
    /// recorded, which is OCCT's null `Geom2dAdaptor_Curve`.
    pub fn edge(&self) -> Option<&Arc<dyn Curve2d>> {
        assert!(!self.rejected, "TopClass_FaceClassifier::Edge : rejected");
        self.my_edge.as_ref()
    }

    /// `EdgeParameter()` (`cxx:79-83` + `pxx:171-175`).
    pub fn edge_parameter(&self) -> f64 {
        assert!(
            !self.rejected,
            "TopClass_FaceClassifier::EdgeParameter : rejected"
        );
        self.my_edge_parameter
    }

    /// `Position()` (`hxx:73`).
    pub fn position(&self) -> Option<IntRes2dPosition> {
        self.my_position
    }
}
