//! Port of `TopOpeBRep_VPointInterClassifier`
//! (`TKTopAlgo/TopOpeBRep/TopOpeBRep_VPointInterClassifier.hxx:32-63`,
//! `.cxx:26-290`).
//!
//! `mySlowFaceClassifier` is the `BRepClass_FaceClassifier` of
//! [`crate::brep_class::FClassifier`] (faithful `BRepClass_FClassifier` over
//! `TopClass_FaceClassifier.pxx`) — the edge information `SlowClassifyOnBoundary`
//! needs is only available from that classifier.

use occt_core::gp::{GpPnt, GpPnt2d};

use occt_geom::extrema_pc::point_curve_extrema_all;

use crate::abs::State;
use crate::brep_class::{FClassifier, FaceExplorer};
use crate::brep_tool::BRepTool;
use crate::fclass2d::FaceState;
use crate::shape::{Edge, Face};
use crate::topo_tools::explorer_shapes;
use crate::abs::ShapeType;

use super::point_classifier::PointClassifier;
use super::v_point_inter::VPointInter;

/// `TopAbs_State` of a [`FaceState`].
fn state_of(face_state: FaceState) -> State {
    match face_state {
        FaceState::In => State::In,
        FaceState::Out => State::Out,
        FaceState::On => State::On,
        FaceState::Unknown => State::Unknown,
    }
}

/// `TopOpeBRep_VPointInterClassifier`.
pub struct VPointInterClassifier {
    /// `mySlowFaceClassifier` (`BRepClass_FaceClassifier`).
    slow_face_classifier: FClassifier,
    /// `myState` (`TopAbs_State`).
    state: State,
}

impl Default for VPointInterClassifier {
    fn default() -> Self {
        Self::new()
    }
}

impl VPointInterClassifier {
    /// `TopOpeBRep_VPointInterClassifier()` (`cxx:42-45`).
    pub fn new() -> Self {
        Self {
            slow_face_classifier: FClassifier::new(),
            state: State::Unknown,
        }
    }

    /// `VPointPosition(F, VP, ShapeIndex, PC, AssumeINON, Tol)` (`cxx:51-179`).
    pub fn v_point_position(
        &mut self,
        f: &Face,
        vp: &mut VPointInter,
        shape_index: i32,
        pc: &mut PointClassifier,
        assume_in_on: bool,
        tol: f64,
    ) -> Result<State, String> {
        self.state = State::Unknown;

        let (mut u, mut v);
        match shape_index {
            1 => {
                let (a, b) = vp.parameters_on_s1()?;
                u = a;
                v = b;
                if vp.is_on_dom_s1()? {
                    vp.set_state(State::On, 1)?;
                    let arc = vp.arc_on_s1()?.clone();
                    let pe = vp.parameter_on_arc1()?;
                    vp.set_edge_on(&arc, pe, 1);
                    self.state = State::On;
                    return Ok(self.state);
                }
            }
            2 => {
                let (a, b) = vp.parameters_on_s2()?;
                u = a;
                v = b;
                if vp.is_on_dom_s2()? {
                    vp.set_state(State::On, 2)?;
                    let arc = vp.arc_on_s2()?.clone();
                    let pe = vp.parameter_on_arc2()?;
                    vp.set_edge_on(&arc, pe, 2);
                    self.state = State::On;
                    return Ok(self.state);
                }
            }
            _ => return Err("VPointClassifier : wrong Index".into()),
        }

        if self.state == State::On {
            return Ok(self.state);
        }

        let ff = Face(f.0.clone());
        crate::topopebrep::tool::shape_tool::adjust_on_periodic(&ff, &mut u, &mut v);
        let p2d = GpPnt2d::new(u, v);
        let statefast = pc.classify(&ff, p2d, tol);
        self.state = statefast;
        vp.set_state(self.state, shape_index)?;

        // AssumeINON = True <=> every VP produced by the intersections is
        // considered to be on the boundary (ON) or inside (IN) of face FF.
        let vpsi = vp.shape_index();
        if assume_in_on && shape_index == vpsi {
            self.perform_slow(&ff, &p2d, tol);
            self.state = state_of(self.slow_face_classifier.state());
            if self.state == State::On {
                let p3d = vp.value()?.clone();
                self.state = slow_classify_on_boundary(&p3d, &p2d, &mut self.slow_face_classifier, &ff);
                if self.state == State::On {
                    self.set_edge_on(vp, shape_index);
                }
            } else if self.state == State::Out {
                self.state = State::In;
                // or myState = TopAbs_OUT, which leads to VP.UnKeep() if the
                // intersection data is preferred over the classifier.
            }
        } else if !assume_in_on {
            if statefast == State::Out || statefast == State::On {
                self.perform_slow(&ff, &p2d, tol);
                self.state = state_of(self.slow_face_classifier.state());
                if self.state == State::On {
                    let p3d = vp.value()?.clone();
                    self.state =
                        slow_classify_on_boundary(&p3d, &p2d, &mut self.slow_face_classifier, &ff);
                    if self.state == State::On {
                        self.set_edge_on(vp, shape_index);
                    }
                }
            }
        } else if shape_index != vpsi {
            // AssumeINON = True
            if statefast == State::On || statefast == State::Out {
                self.perform_slow(&ff, &p2d, tol);
                self.state = state_of(self.slow_face_classifier.state());

                if self.state == State::On || self.state == State::Out {
                    let p3d = vp.value()?.clone();
                    self.state =
                        slow_classify_on_boundary(&p3d, &p2d, &mut self.slow_face_classifier, &ff);
                    if self.state == State::On {
                        self.set_edge_on(vp, shape_index);
                    }
                }
            }
        }

        vp.set_state(self.state, shape_index)?;
        Ok(self.state)
    }

    /// `mySlowFaceClassifier.Perform(FF, p2d, Tol)`.
    fn perform_slow(&mut self, face: &Face, p2d: &GpPnt2d, tol: f64) {
        let mut fexp = FaceExplorer::new(face);
        self.slow_face_classifier.perform(&mut fexp, p2d, tol);
    }

    /// `VP.EdgeON(mySlowFaceClassifier.Edge().Edge(), EdgeParameter(), Index)`.
    fn set_edge_on(&self, vp: &mut VPointInter, shape_index: i32) {
        if let Some(edge) = self.slow_face_classifier.edge().edge() {
            vp.set_edge_on(&edge.0, self.slow_face_classifier.edge_parameter(), shape_index);
        }
    }

    /// `Edge()` (`cxx:184-193`): the edge of the last `VPointPosition` when
    /// the state is ON, else no edge (`myNullShape`).
    pub fn edge(&self) -> Option<&Edge> {
        if self.state == State::On {
            self.slow_face_classifier.edge().edge()
        } else {
            None
        }
    }

    /// `EdgeParameter()` (`cxx:197-206`).
    pub fn edge_parameter(&self) -> f64 {
        if self.state == State::On {
            self.slow_face_classifier.edge_parameter()
        } else {
            0.0
        }
    }

    /// `myState`.
    pub fn state(&self) -> State {
        self.state
    }
}

/// `SlowClassifyOnBoundary` (`cxx:214-287`).
fn slow_classify_on_boundary(
    point_to_classify: &GpPnt,
    point2d_to_classify: &GpPnt2d,
    slow_classifier: &mut FClassifier,
    face: &Face,
) -> State {
    let a_parameter_on_edge = slow_classifier.edge_parameter();
    // `theSlowClassifier.Edge().Edge()`. Reaching here with no edge is the
    // OCCT null-shape dereference in `BRep_Tool::Curve`; the port treats it
    // as "no edge curve", i.e. the `!anEdgeCurve.IsNull()` else arm below.
    let Some(an_edge) = slow_classifier.edge().edge().cloned() else {
        return State::On;
    };

    let (parf, parl) = BRepTool::edge_parameters(&an_edge);
    let Some(an_edge_curve) = BRepTool::edge_curve(&an_edge) else {
        return State::On;
    };

    let mut minparameterdiff = parl - parf;
    let mut a_distance_to_compare = 0.0;
    let mut same_with_vertex = false;
    for shape in explorer_shapes(&an_edge.0, ShapeType::Vertex) {
        let Some(vertex) = crate::shape::Vertex::wrap(shape) else {
            continue;
        };
        let a_vertex_tolerance = BRepTool::vertex_tolerance(&vertex);
        let an_edge_vertex_point = BRepTool::vertex_point(&vertex);
        if point_to_classify.is_equal_with_tolerance(&an_edge_vertex_point, a_vertex_tolerance) {
            same_with_vertex = true;
            break;
        }
    }
    if same_with_vertex {
        return State::On;
    }

    // `GeomAPI_ProjectPointOnCurve aProjTool(thePointToClassify, anEdgeCurve)`.
    let projections = point_curve_extrema_all(an_edge_curve.as_ref(), point_to_classify);
    for proj in projections.iter() {
        let curparamdiff = (proj.u1 - a_parameter_on_edge).abs();
        if curparamdiff < minparameterdiff {
            minparameterdiff = curparamdiff;
            a_distance_to_compare = proj.distance;
        }
    }

    let an_edge_tolerance = BRepTool::edge_tolerance(&an_edge);

    if projections.is_empty() || a_distance_to_compare >= an_edge_tolerance {
        let Some(an_edge_pcurve) = crate::boptools_2d::curve_on_surface_range(&an_edge, face)
            .map(|(c, _f, _l)| c)
        else {
            return State::Out;
        };
        let a_point2d_on_edge = an_edge_pcurve.d0(a_parameter_on_edge);
        let a_tol2d = point2d_to_classify.distance(&a_point2d_on_edge) / 3.0;
        let mut fexp = FaceExplorer::new(face);
        slow_classifier.perform(&mut fexp, point2d_to_classify, a_tol2d);
        if state_of(slow_classifier.state()) == State::In {
            return State::In;
        }
        return State::Out;
    }

    State::On
}
