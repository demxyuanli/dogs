//! Port of `BRepClass_FClassifier` (`BRepClass_FClassifier.hxx:26-88`,
//! `BRepClass_FClassifier.cxx:25-100`) over
//! `TopClass_FaceClassifier::Perform` / `State` / `Edge` / `EdgeParameter`
//! (`TKGeomAlgo/TopClass/TopClass_FaceClassifier.pxx:28-176`).
//!
//! The template's `FClass2d` parameter is
//! `BRepClass_FClass2dOfFClassifier` = [`BRepClassFacePassiveClassifier`], and
//! `TheFaceExplorer` is [`FaceExplorer`].

use occt_core::gp::{GpLin2d, GpPnt2d};

use crate::abs::Orientation;
use crate::fclass2d::FaceState;

use super::edge::BRepClassEdge;
use super::face_explorer::FaceExplorer;
use super::face_passive_classifier::BRepClassFacePassiveClassifier;

/// `BRepClass_FClassifier`.
pub struct FClassifier {
    classifier: BRepClassFacePassiveClassifier,
    edge: BRepClassEdge,
    edge_parameter: f64,
    position: occt_core::intres2d::IntRes2dPosition,
    rejected: bool,
    no_wires: bool,
}

impl Default for FClassifier {
    fn default() -> Self {
        Self::new()
    }
}

impl FClassifier {
    /// `BRepClass_FClassifier()` (`cxx:25-30`).
    pub fn new() -> Self {
        Self {
            classifier: BRepClassFacePassiveClassifier::new(),
            edge: BRepClassEdge::new(),
            edge_parameter: 0.0,
            position: occt_core::intres2d::IntRes2dPosition::Head,
            rejected: false,
            no_wires: true,
        }
    }

    /// `BRepClass_FClassifier(F, P, Tol)` (`cxx:34-44`).
    pub fn with_args(fexp: &mut FaceExplorer, p: &GpPnt2d, tol: f64) -> Self {
        let mut out = Self::new();
        out.perform(fexp, p, tol);
        out
    }

    /// `Perform(F, P, Tol)` (`cxx:48-62`).
    pub fn perform(&mut self, fexp: &mut FaceExplorer, p: &GpPnt2d, tol: f64) {
        // `TopClass_FaceClassifier::Perform` (`pxx:29-145`).
        let mut point = *p;
        while !fexp.check_point(&mut point) {}

        self.rejected = fexp.reject(&point);
        if self.rejected {
            return;
        }

        let mut line = GpLin2d::default();
        let mut param = 0.0;
        let mut is_valid_segment = fexp.segment(&point, &mut line, &mut param);
        let mut state = FaceState::Unknown;
        self.no_wires = true;

        while is_valid_segment {
            self.classifier.reset(&line, param, tol);

            fexp.init_wires();
            while fexp.more_wires() {
                self.no_wires = false;
                let is_w_reject = fexp.reject_wire(&line, self.classifier.parameter());
                if !is_w_reject {
                    // Test this wire.
                    fexp.init_edges();
                    while fexp.more_edges() {
                        let is_e_reject = fexp.reject_edge(&line, self.classifier.parameter());
                        if is_e_reject {
                            fexp.next_edge();
                            continue;
                        }
                        // Test this edge.
                        let Some((an_edge, an_edge_ori)) = fexp.current_edge() else {
                            fexp.next_edge();
                            continue;
                        };
                        if an_edge_ori == Orientation::Forward
                            || an_edge_ori == Orientation::Reversed
                        {
                            self.classifier.compare(&an_edge, an_edge_ori);
                            let closest = self.classifier.closest_intersection();
                            if closest != 0 {
                                // Save the closest edge.
                                let nb_points = self.classifier.intersector().nb_points() as i32;
                                self.edge = an_edge;

                                let p_inter = if closest <= nb_points {
                                    *self.classifier.intersector().point(closest as usize)
                                } else {
                                    let ci = closest - nb_points;
                                    let seg_index = ((ci + 1) / 2) as usize;
                                    if ci & 1 != 0 {
                                        *self.classifier.intersector().segment(seg_index).first_point()
                                    } else {
                                        *self.classifier.intersector().segment(seg_index).last_point()
                                    }
                                };

                                self.position =
                                    p_inter.transition_of_second().position_on_curve();
                                self.edge_parameter = p_inter.param_on_second();
                            }
                            // If we are ON, we stop.
                            state = self.classifier.state();
                            if state == FaceState::On {
                                return;
                            }
                        }
                        fexp.next_edge();
                    }

                    // If we are out of the wire we stop.
                    state = self.classifier.state();
                    if state == FaceState::Out {
                        return;
                    }
                }
                fexp.next_wire();
            }

            if !self.classifier.is_head_or_end() && state != FaceState::Unknown {
                break;
            }

            // Bad case for classification. Trying to get another segment.
            is_valid_segment = fexp.other_segment(&point, &mut line, &mut param);
        }
    }

    /// `State()` (`cxx:66-69` + `pxx:147-157`).
    pub fn state(&self) -> FaceState {
        if self.rejected {
            FaceState::Out
        } else if self.no_wires {
            FaceState::In
        } else {
            self.classifier.state()
        }
    }

    /// `Rejected()`.
    pub fn rejected(&self) -> bool {
        self.rejected
    }

    /// `NoWires()`.
    pub fn no_wires(&self) -> bool {
        self.no_wires
    }

    /// `Edge()` (`cxx:73-76` + `pxx:162-167`).
    pub fn edge(&self) -> &BRepClassEdge {
        self.edge_ref_checked()
    }

    /// `EdgeParameter()` (`cxx:80-83` + `pxx:171-175`).
    pub fn edge_parameter(&self) -> f64 {
        if self.rejected {
            panic!("TopClass_FaceClassifier::EdgeParameter : rejected");
        }
        self.edge_parameter
    }

    /// `Position()`.
    pub fn position(&self) -> occt_core::intres2d::IntRes2dPosition {
        self.position
    }

    fn edge_ref_checked(&self) -> &BRepClassEdge {
        if self.rejected {
            panic!("TopClass_FaceClassifier::Edge : rejected");
        }
        &self.edge
    }
}
