//! `TopOpeBRepBuild_PaveClassifier` — compares vertices on an edge.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_PaveClassifier.hxx:31-89`,
//! `.cxx:24-391`.
//!
//! A vertex `V1` is inside a vertex `V2` if `V1` is on the part of the curve
//! defined by `V2` (`hxx:34-42`):
//! - if `V2` is `FORWARD`, `V1` must be after `V2` on the curve;
//! - if `V2` is `REVERSED`, `V1` must be before `V2` on the curve;
//! - if `V2` is `INTERNAL`, `V1` is always inside;
//! - if `V2` is `EXTERNAL`, `V1` is never inside.

use occt_core::elib::clib::in_period;
use occt_core::precision::PCONFUSION;

use crate::abs::{Orientation, State};
use crate::brep_tool::BRepTool;
use crate::shape::TopoShape;
use crate::topo_tools::{edge_vertices, is_same};
use crate::Edge;

use super::loop_classifier::LoopClassifier;
use super::loop_shape::LoopHandle;

/// `class TopOpeBRepBuild_PaveClassifier : public
/// TopOpeBRepBuild_LoopClassifier` (`hxx:44`).
#[derive(Debug, Clone)]
pub struct PaveClassifier {
    /// `myEdge` (`hxx:76`).
    my_edge: TopoShape,
    /// `myEdgePeriodic` (`hxx:77`).
    my_edge_periodic: bool,
    /// `myFirst` (`hxx:78`).
    my_first: f64,
    /// `myPeriod` (`hxx:79`).
    my_period: f64,
    /// `mySameParameters` (`hxx:80`).
    my_same_parameters: bool,
    /// `myClosedVertices` (`hxx:81`).
    my_closed_vertices: bool,
    /// `myP1` (`hxx:82`).
    my_p1: f64,
    /// `myP2` (`hxx:83`).
    my_p2: f64,
    /// `myO1` (`hxx:84`).
    my_o1: Orientation,
    /// `myO2` (`hxx:85`).
    my_o2: Orientation,
    /// `myCas1` (`hxx:86`).
    my_cas1: i32,
    /// `myCas2` (`hxx:87`).
    my_cas2: i32,
}

impl PaveClassifier {
    /// `TopOpeBRepBuild_PaveClassifier(E)` (`cxx:27-106`).
    pub fn new(e: &TopoShape) -> Self {
        let mut c = Self {
            my_edge: e.clone(),
            my_edge_periodic: false,
            my_first: 0.0,
            my_period: 0.0,
            my_same_parameters: false,
            my_closed_vertices: false,
            my_p1: 0.0,
            my_p2: 0.0,
            my_o1: Orientation::Forward,
            my_o2: Orientation::Forward,
            my_cas1: 0,
            my_cas2: 0,
        };

        let edge = Edge(c.my_edge.clone());
        if !BRepTool::is_degenerated(&edge) {
            let (f, l) = BRepTool::edge_parameters(&edge);
            if let Some(curve) = BRepTool::edge_curve(&edge) {
                if curve.is_periodic() {
                    let (v1, v2) = edge_vertices(&edge); // v1 FORWARD, v2 REVERSED
                    match (v1, v2) {
                        (Some(v1), Some(v2)) => {
                            // --- the edge has vertices
                            c.my_first = f;
                            let f_c = curve.first_parameter();
                            let l_c = curve.last_parameter();
                            c.my_period = l_c - f_c;
                            c.my_edge_periodic = is_same(&v1.0, &v2.0);
                            c.my_same_parameters = c.my_edge_periodic;
                            if c.my_same_parameters {
                                c.my_first = BRepTool::parameter_of_vertex(&c.my_edge, &v1.0);
                            }
                        }
                        _ => {
                            // --- the edge has no vertices
                            c.my_first = f;
                            c.my_period = l - f;
                            c.my_edge_periodic = true;
                            c.my_same_parameters = false;
                        }
                    }
                }
            }
        } // ! degenerated

        c
    }

    /// `CompareOnNonPeriodic()` (`cxx:110-152`).
    fn compare_on_non_periodic(&self) -> State {
        let mut state = State::Unknown;
        let mut lower = false;
        match self.my_o2 {
            Orientation::Forward => {
                lower = false;
            }
            Orientation::Reversed => {
                lower = true;
            }
            Orientation::Internal => {
                state = State::In;
            }
            Orientation::External => {
                state = State::Out;
            }
        }

        if state == State::Unknown {
            if self.my_p1 == self.my_p2 {
                if self.my_o1 == self.my_o2 {
                    state = State::In;
                } else {
                    state = State::Out;
                }
            } else if self.my_p1 < self.my_p2 {
                if lower {
                    state = State::In;
                } else {
                    state = State::Out;
                }
            } else {
                if lower {
                    state = State::Out;
                } else {
                    state = State::In;
                }
            }
        }

        state
    }

    /// `static double AdjustCase(p1, o, first, period, tol, cas)`
    /// (`cxx:157-187`).
    pub fn adjust_case(
        p1: f64,
        o: Orientation,
        first: f64,
        period: f64,
        tol: f64,
        cas: &mut i32,
    ) -> f64 {
        let p2;
        if (p1 - first).abs() < tol {
            // p1 is first
            if o == Orientation::Reversed {
                p2 = p1 + period;
                *cas = 1;
            } else {
                p2 = p1;
                *cas = 2;
            }
        } else {
            // p1 is not on first
            let last = first + period;
            if (p1 - last).abs() < tol {
                // p1 is on last
                p2 = p1;
                *cas = 3;
            } else {
                // p1 is not on last
                p2 = in_period(p1, first, last);
                *cas = 4;
            }
        }
        p2
    }

    /// `AdjustOnPeriodic()` (`cxx:191-238`).
    fn adjust_on_periodic(&mut self) {
        if !self.to_adjust_on_periodic() {
            return;
        }

        let tol = PCONFUSION;

        if self.my_same_parameters {
            let mut cas = self.my_cas1;
            self.my_p1 = Self::adjust_case(self.my_p1, self.my_o1, self.my_first, self.my_period, tol, &mut cas);
            self.my_cas1 = cas;
            let mut cas = self.my_cas2;
            self.my_p2 = Self::adjust_case(self.my_p2, self.my_o2, self.my_first, self.my_period, tol, &mut cas);
            self.my_cas2 = cas;
        } else if self.my_o1 != self.my_o2 {
            if self.my_o1 == Orientation::Forward {
                let mut cas = self.my_cas2;
                self.my_p2 = Self::adjust_case(self.my_p2, self.my_o2, self.my_p1, self.my_period, tol, &mut cas);
                self.my_cas2 = cas;
            }
            if self.my_o2 == Orientation::Forward {
                let mut cas = self.my_cas1;
                self.my_p1 = Self::adjust_case(self.my_p1, self.my_o1, self.my_p2, self.my_period, tol, &mut cas);
                self.my_cas1 = cas;
            }
        }
    }

    /// `ToAdjustOnPeriodic()` (`cxx:243-247`).
    fn to_adjust_on_periodic(&self) -> bool {
        self.my_same_parameters || (self.my_o1 != self.my_o2)
    }

    /// `CompareOnPeriodic()` (`cxx:251-287`).
    fn compare_on_periodic(&mut self) -> State {
        let state;

        if self.to_adjust_on_periodic() {
            state = self.compare_on_non_periodic();
        } else if self.my_o1 == Orientation::Forward {
            state = State::Out;
            self.my_cas1 = 5;
            self.my_cas2 = 5;
        } else if self.my_o1 == Orientation::Reversed {
            state = State::Out;
            self.my_cas1 = 6;
            self.my_cas2 = 6;
        } else {
            state = State::Out;
            self.my_cas1 = 7;
            self.my_cas2 = 7;
        }

        state
    }

    /// `SetFirstParameter(P)` (`cxx:369-377`).
    pub fn set_first_parameter(&mut self, p: f64) {
        self.my_first = p;
        self.my_same_parameters = true;
    }

    /// `ClosedVertices(B)` (`cxx:381-391`).
    pub fn set_closed_vertices(&mut self, closed: bool) {
        self.my_closed_vertices = closed;
    }

    /// `myClosedVertices`.
    pub fn closed_vertices(&self) -> bool {
        self.my_closed_vertices
    }
}

impl LoopClassifier for PaveClassifier {
    /// `Compare(L1, L2)` (`cxx:291-365`).
    fn compare(&mut self, l1: &LoopHandle, l2: &LoopHandle) -> State {
        {
            let pv1 = l1.read().expect("poisoned Loop lock");
            let pv2 = l2.read().expect("poisoned Loop lock");
            let p1 = pv1.as_pave().expect("TopOpeBRepBuild_Pave expected");
            let p2 = pv2.as_pave().expect("TopOpeBRepBuild_Pave expected");

            self.my_cas1 = 0; // debug
            self.my_cas2 = 0;
            self.my_o1 = p1.vertex.orientation();
            self.my_o2 = p2.vertex.orientation();
            self.my_p1 = p1.parameter();
            self.my_p2 = p2.parameter();
        }

        if self.my_edge_periodic && self.to_adjust_on_periodic() {
            self.adjust_on_periodic();
        }

        if self.my_edge_periodic {
            self.compare_on_periodic()
        } else {
            self.compare_on_non_periodic()
        }
    }
}
