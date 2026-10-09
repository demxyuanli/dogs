//! Port of `TopClass_Classifier2d` (`TKGeomAlgo/TopClass/TopClass_Classifier2d.pxx`):
//! the `Reset` and `Compare` templates. The OCCT template takes its state as
//! references; here the state is the struct and the intersector, transition
//! and edge are passed to `compare`.

use occt_core::gp::{GpDir, GpDir2d, GpLin2d, GpPnt2d, GpVec};
use occt_core::intres2d::{
    IntRes2dIntersectionPoint, IntRes2dPosition, IntRes2dSituation, IntRes2dTypeTrans,
};

use crate::abs::Orientation;
use crate::fclass2d::FaceState;
use crate::intpatch::CurveTransition;

use super::edge::BRepClassEdge;
use super::intersector::BRepClassIntersector;

/// State fields of `TopClass_Classifier2d` (the `Reset`/`Compare` arguments)
/// together with the `BRepClass_FacePassiveClassifier` flags `myIsSet` and
/// `myClosest`.
pub struct TopClassClassifier2d {
    lin: GpLin2d,
    param: f64,
    tolerance: f64,
    state: FaceState,
    first_compare: bool,
    first_trans: bool,
    closest: i32,
    is_set: bool,
    is_head_or_end: bool,
}

impl Default for TopClassClassifier2d {
    fn default() -> Self {
        Self::new()
    }
}

impl TopClassClassifier2d {
    /// Member initialisation of `BRepClass_FacePassiveClassifier()`
    /// (`BRepClass_FacePassiveClassifier.cxx:24-35`): the line is the default
    /// `gp_Lin2d` (origin, X direction).
    pub fn new() -> Self {
        Self {
            lin: GpLin2d::from_pnt_dir(
                GpPnt2d::new(0.0, 0.0),
                GpDir2d::new(1.0, 0.0).expect("unit X direction"),
            ),
            param: 0.0,
            tolerance: 0.0,
            state: FaceState::Unknown,
            first_compare: true,
            first_trans: true,
            closest: 0,
            is_set: false,
            is_head_or_end: false,
        }
    }

    /// `myIsSet`.
    pub fn is_set(&self) -> bool {
        self.is_set
    }

    /// `myParam`.
    pub fn param(&self) -> f64 {
        self.param
    }

    /// `myClosest`.
    pub fn closest(&self) -> i32 {
        self.closest
    }

    /// `myState`.
    pub fn state(&self) -> FaceState {
        self.state
    }

    /// `myIsHeadOrEnd`.
    pub fn is_head_or_end(&self) -> bool {
        self.is_head_or_end
    }

    /// `TopClass_Classifier2d::Reset` (`TopClass_Classifier2d.pxx:30-52`).
    pub fn reset(&mut self, l: &GpLin2d, p: f64, tol: f64) {
        self.lin = *l;
        self.param = p;
        self.tolerance = tol;
        self.state = FaceState::Unknown;
        self.first_compare = true;
        self.first_trans = true;
        self.closest = 0;
        self.is_set = true;
        self.is_head_or_end = false;
    }

    /// `TopClass_Classifier2d::Compare` (`TopClass_Classifier2d.pxx:58-225`).
    pub(crate) fn compare(
        &mut self,
        intersector: &mut BRepClassIntersector,
        trans: &mut CurveTransition,
        edge: &BRepClassEdge,
        orientation: Orientation,
    ) {
        // Intersect the edge and the segment.
        self.closest = 0;
        intersector.perform(&self.lin, self.param, self.tolerance, edge);
        if !intersector.is_done() {
            return;
        }
        if intersector.nb_points() == 0 && intersector.nb_segments() == 0 {
            return;
        }

        // Find the closest point.
        let nb_points = intersector.nb_points();
        let mut d_min = f64::MAX;
        let mut p_closest: Option<IntRes2dIntersectionPoint> = None;
        for a_point in 1..=nb_points {
            let p_inter = *intersector.point(a_point);
            // Test for ON.
            if p_inter.transition_of_first().position_on_curve() == IntRes2dPosition::Head {
                self.closest = a_point as i32;
                self.state = FaceState::On;
                return;
            }
            let param_first = p_inter.param_on_first();
            if param_first < d_min {
                self.closest = a_point as i32;
                p_closest = Some(p_inter);
                d_min = param_first;
            }
        }

        // For the segments only the first point is tested.
        let nb_segments = intersector.nb_segments();
        for a_segment in 1..=nb_segments {
            let p_inter = *intersector.segment(a_segment).first_point();
            if p_inter.transition_of_first().position_on_curve() == IntRes2dPosition::Head {
                self.closest = (nb_points + a_segment + a_segment - 1) as i32;
                self.state = FaceState::On;
                return;
            }
            let param_first = p_inter.param_on_first();
            if param_first < d_min {
                self.closest = (nb_points + a_segment + a_segment - 1) as i32;
                p_closest = Some(p_inter);
                d_min = param_first;
            }
        }

        // No point was found.
        if self.closest == 0 {
            return;
        }

        // An INTERNAL or EXTERNAL edge needs no transition analysis.
        if orientation == Orientation::Internal {
            self.state = FaceState::In;
            return;
        } else if orientation == Orientation::External {
            self.state = FaceState::Out;
            return;
        }

        if !self.first_compare && d_min > self.param {
            return;
        }

        // Process the closest point found at `d_min` on the line.
        self.first_compare = false;
        if self.param > d_min {
            self.first_trans = true;
        }
        self.param = d_min;

        // `closest` is set together with `p_closest`, so it is present here.
        let Some(p_closest) = p_closest else {
            return;
        };
        let t2 = *p_closest.transition_of_second();
        self.is_head_or_end = matches!(
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

        if !self.is_head_or_end {
            // The closest point is inside the edge.
            self.state = match seg_trans {
                Orientation::Forward | Orientation::External => FaceState::Out,
                Orientation::Reversed | Orientation::Internal => FaceState::In,
            };
            return;
        }

        // The closest point is the Head or End of the edge: update the
        // complex transition.
        let (tang2d, norm2d, curv) = intersector.local_geometry(edge, p_closest.param_on_second());
        let tang = GpVec::new(tang2d.x(), tang2d.y(), 0.0);
        let norm = GpDir::new(norm2d.x(), norm2d.y(), 0.0).expect("unit normal");
        if self.first_trans {
            let dir = self.lin.direction();
            let tgt = GpDir::new(dir.x(), dir.y(), 0.0).expect("unit line direction");
            trans.reset_tgt(&tgt);
            self.first_trans = false;
        }
        let ort = if t2.position_on_curve() == IntRes2dPosition::Head {
            Orientation::Forward
        } else {
            Orientation::Reversed
        };
        trans.compare(f64::EPSILON, &tang, &norm, curv, seg_trans, ort);
        self.state = trans.state_before();
    }
}
