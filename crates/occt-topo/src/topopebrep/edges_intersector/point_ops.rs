//! The point sequence of `TopOpeBRep_EdgesIntersector`: `MakePoints2d`,
//! `ReduceSegment(s)` and the `Point1` / `Status1` / `Transition1` family that
//! builds them.
//!
//! Source: `TKBool/TopOpeBRep/TopOpeBRep_EdgesIntersector.cxx:708-954` and
//! `TopOpeBRep_EdgesIntersector_1.cxx:32-463`.

use occt_core::gp::GpPnt2d;
use occt_core::intres2d::{
    IntRes2dIntersectionPoint, IntRes2dIntersectionSegment, IntRes2dPosition,
    IntRes2dSituation, IntRes2dTypeTrans,
};
use occt_core::precision::{Precision, CONFUSION, PCONFUSION};
use occt_geom2d::geom2d_int::{geom2d_curve_type, GeomAbsCurveType2d};

use crate::abs::{Orientation, ShapeType, State};
use crate::brep_tool::BRepTool;
use crate::shape::Vertex;
use crate::topexp::Explorer;
use crate::topo_tools::edge_vertices;
use crate::topopebrep::ds::kind::Config;
use crate::topopebrep::ds::transition::Transition;
use crate::topopebrep::interf::p2d_status::P2DStatus;
use crate::topopebrep::interf::point2d::Point2d;
use crate::topopebrep::tool::shape_tool;

use super::EdgesIntersector;

impl EdgesIntersector {
    /// `ComputeSameDomain()` (`cxx:708-767`).
    pub(super) fn compute_same_domain(&mut self) -> bool {
        let (Some(c1), Some(c2)) = (self.curve1.clone(), self.curve2.clone()) else {
            return self.set_same_domain(false);
        };
        let t1 = geom2d_curve_type(c1.as_ref());
        let t2 = geom2d_curve_type(c2.as_ref());

        if !self.has_segment {
            return self.set_same_domain(false);
        }

        let tt = t1 == t2;
        if !tt {
            return self.set_same_domain(false);
        }

        if t1 == GeomAbsCurveType2d::Line {
            return self.set_same_domain(true);
        }

        if t1 != GeomAbsCurveType2d::Circle {
            // `cxx:731-738`: "EdgesSameDomain on NYI curve type".
            return self.set_same_domain(false);
        }

        let (Some(circ1), Some(circ2)) = (c1.gp_circ2d(), c2.gp_circ2d()) else {
            return self.set_same_domain(false);
        };
        let r1 = circ1.radius();
        let r2 = circ2.radius();
        // xpu281098 (cto019D2) tolerance a revoir (`cxx:747`).
        let rr = (r1 - r2).abs() < CONFUSION;
        if !rr {
            return self.set_same_domain(false);
        }

        // `myCurve2 = C2d(myEdge2, myFace1)`: both circle centres are read on
        // the surface of face 1 (`cxx:751-756`).
        let (Some(pnt1), Some(pnt2)) = (
            self.surface1
                .surface()
                .map(|s| s.d0(circ1.location().x(), circ1.location().y())),
            self.surface1
                .surface()
                .map(|s| s.d0(circ2.location().x(), circ2.location().y())),
        ) else {
            return self.set_same_domain(false);
        };
        let dpp = pnt1.distance(&pnt2);
        let tol1 = BRepTool::edge_tolerance(&self.edge1);
        let tol2 = BRepTool::edge_tolerance(&self.edge2);
        let tol = tol1 + tol2;
        let esd = dpp <= tol;
        if esd {
            return self.set_same_domain(true);
        }

        self.set_same_domain(false)
    }

    /// `SetSameDomain(B)` (`cxx:772-776`).
    pub(super) fn set_same_domain(&mut self, b: bool) -> bool {
        self.same_domain = b;
        b
    }

    /// `MakePoints2d()` (`cxx:780-817`).
    pub(super) fn make_points2d(&mut self) {
        self.sp2d.clear();
        let e1ori = self.edge1.0.orientation();
        let e2ori = self.edge2.0.orientation();
        self.init_point1();
        while self.more_point1() {
            let Some(ip) = self.point1() else {
                break;
            };
            let mut p2d = Point2d::new();
            p2d.set_pint(&ip);
            if let Ok(t) = self.transition1(1, e2ori) {
                let _ = p2d.set_transition(1, &t);
            }
            if let Ok(t) = self.transition1(2, e1ori) {
                let _ = p2d.set_transition(2, &t);
            }
            if let Ok(p) = self.parameter1(1) {
                let _ = p2d.set_parameter(1, p);
            }
            if let Ok(p) = self.parameter1(2) {
                let _ = p2d.set_parameter(2, p);
            }
            let isv1 = self.is_vertex1(1);
            let _ = p2d.set_is_vertex(1, isv1);
            if isv1 {
                if let Ok(v) = self.vertex1(1) {
                    let _ = p2d.set_vertex(1, v);
                }
            }
            let isv2 = self.is_vertex1(2);
            let _ = p2d.set_is_vertex(2, isv2);
            if isv2 {
                if let Ok(v) = self.vertex1(2) {
                    let _ = p2d.set_vertex(2, v);
                }
            }
            p2d.set_is_point_of_segment(self.is_point_of_segment1());
            p2d.set_segment_ancestors(0, 0);
            p2d.set_status(self.status1());
            p2d.set_value(&self.value1());
            p2d.set_value2d(ip.value());
            p2d.set_tolerance(self.tolerance_max());
            p2d.set_edges_config(self.edges_config1());
            p2d.set_index(self.index1());
            self.sp2d.push(p2d);
            // `for (...; MorePoint1(); NextPoint1())` (`cxx:785`).
            self.next_point1();
        }
        self.ip2d = 1;
        self.np2d = self.sp2d.len() as i32;
    }

    /// `ReduceSegment(psa, psb, Pn)` (`cxx:822-913`).
    ///
    /// The OCCT method is `const`: it reads nothing but its three arguments,
    /// which is why it is a free function here — `ReduceSegments` needs two
    /// mutable elements of `mysp2d` at once.
    fn reduce_segment(psa: &mut Point2d, psb: &mut Point2d, pn: &mut Point2d) -> bool {
        let mut reduced = false;
        let ixpsa = psa.index();
        let ixpsb = psb.index();

        let pospsa = psa.is_point_of_segment();
        let stspsa = psa.status();
        let tpsa1 = psa.parameter(1).unwrap_or(0.0);
        let tpsa2 = psa.parameter(2).unwrap_or(0.0);
        let psa1_trans = psa.transition(1).map(|t| t.clone()).unwrap_or_default();
        let psa2_trans = psa.transition(2).map(|t| t.clone()).unwrap_or_default();

        let pospsb = psb.is_point_of_segment();
        let stspsb = psb.status();
        let tpsb1 = psb.parameter(1).unwrap_or(0.0);
        let tpsb2 = psb.parameter(2).unwrap_or(0.0);
        let psb1_trans = psb.transition(1).map(|t| t.clone()).unwrap_or_default();
        let psb2_trans = psb.transition(2).map(|t| t.clone()).unwrap_or_default();

        let conda = pospsa && stspsa == P2DStatus::SegmentFirst;
        let condb = pospsb && stspsb == P2DStatus::SegmentLast;
        let cond = conda && condb;

        if cond {
            reduced = true;

            let tm1 = (tpsa1 + tpsb1) / 2.0;
            let _ = pn.set_parameter(1, tm1);
            let tm2 = (tpsa2 + tpsb2) / 2.0;
            let _ = pn.set_parameter(2, tm2);

            let mut tn1 = Transition::default();
            tn1.before(psa1_trans.before_state(), psa1_trans.shape_before());
            tn1.after(psb1_trans.after_state(), psb1_trans.shape_after());
            let _ = pn.set_transition(1, &tn1);
            let mut tn2 = Transition::default();
            tn2.before(psa2_trans.before_state(), psa2_trans.shape_before());
            tn2.after(psb2_trans.after_state(), psb2_trans.shape_after());
            let _ = pn.set_transition(2, &tn2);

            let p3dpsa = psa.value();
            let p3dpsb = psb.value();
            let p3dn = occt_core::gp::GpPnt::new(
                (p3dpsa.x() + p3dpsb.x()) / 2.0,
                (p3dpsa.y() + p3dpsb.y()) / 2.0,
                (p3dpsa.z() + p3dpsb.z()) / 2.0,
            );
            pn.set_value(&p3dn);
            let p2dpsa = psa.value2d();
            let p2dpsb = psb.value2d();
            let p2dn = GpPnt2d::new((p2dpsa.x() + p2dpsb.x()) / 2.0, (p2dpsa.y() + p2dpsb.y()) / 2.0);
            pn.set_value2d(&p2dn);

            let tolpsa = psa.tolerance();
            let tolpsb = psb.tolerance();
            let toln = (tolpsa + tolpsb) * 1.5;
            pn.set_tolerance(toln);

            pn.set_is_point_of_segment(false);
            pn.set_segment_ancestors(ixpsa, ixpsb);
            psa.set_keep(false);
            psb.set_keep(false);

            let cpsa = psa.edges_config();

            pn.set_edges_config(cpsa);

            if let Ok(v) = psa.vertex(1) {
                let _ = pn.set_vertex(1, v);
            }
            if let Ok(v) = psa.vertex(2) {
                let _ = pn.set_vertex(2, v);
            }
            if let Ok(v) = psb.vertex(1) {
                let _ = pn.set_vertex(1, v);
            }
            if let Ok(v) = psb.vertex(2) {
                let _ = pn.set_vertex(2, v);
            }
        }

        reduced
    }

    /// `ReduceSegments()` (`cxx:918-947`).
    pub(super) fn reduce_segments(&mut self) {
        let condredu = self.has_segment && !self.same_domain;
        if !condredu {
            return;
        }

        // `mysp2d` is taken out so the two neighbours can be borrowed mutably;
        // `pn.SetIndex(++mynp2d)` grows `mynp2d` while `np` keeps the value
        // read at `cxx:928`.
        let mut points = std::mem::take(&mut self.sp2d);
        let mut ip = 1;
        let np = self.np2d;
        while ip < np {
            let (psa, psb) = {
                let (left, right) = points.split_at_mut(ip as usize);
                (&mut left[(ip - 1) as usize], &mut right[0])
            };
            let mut pn = Point2d::new();
            let reduced = Self::reduce_segment(psa, psb, &mut pn);
            if reduced {
                self.np2d += 1;
                pn.set_index(self.np2d);
                points.push(pn);
            }
            ip += 1;
        }
        self.sp2d = points;

        self.lseg.clear();
        self.nb_segments = self.lseg.len() as i32;
        self.has_segment = self.nb_segments != 0;
        self.true_nb_points = self.nb_points + 2 * self.nb_segments;
    }

    /// `Segment1()` (`_1.cxx:32-41`).
    pub(super) fn segment1(&self) -> Result<&IntRes2dIntersectionSegment, String> {
        if !self.is_point_of_segment1() {
            return Err("TopOpeBRep_EdgesIntersector : Not a segment point".to_string());
        }
        let iseg = 1 + (self.point_index - self.nb_points - 1) / 2;
        self.lseg
            .get((iseg - 1) as usize)
            .ok_or_else(|| "TopOpeBRep_EdgesIntersector : Not a segment point".to_string())
    }

    /// `IsOpposite1()` (`_1.cxx:45-48`).
    pub(super) fn is_opposite1(&self) -> Result<bool, String> {
        Ok(self.segment1()?.is_opposite())
    }

    /// `InitPoint1()` (`_1.cxx:52-58`).
    pub(super) fn init_point1(&mut self) {
        self.point_index = 1;
        self.is_vertex_point_index = 0;
        self.is_vertex_index = 0;
        self.is_vertex_value = false;
    }

    /// `MorePoint1()` (`_1.cxx:62-65`).
    pub(super) fn more_point1(&self) -> bool {
        self.point_index <= self.true_nb_points
    }

    /// `NextPoint1()` (`_1.cxx:69-72`).
    pub(super) fn next_point1(&mut self) {
        self.point_index += 1;
    }

    /// `Point1()` (`_1.cxx:76-93`).
    pub(super) fn point1(&self) -> Option<IntRes2dIntersectionPoint> {
        if !self.is_point_of_segment1() {
            // Point is an intersection point.
            self.lpnt.get((self.point_index - 1) as usize).copied()
        } else {
            // Point is a point of segment.
            let i = self.point_index - self.nb_points - 1;
            let seg = self.segment1().ok()?;
            if i % 2 == 0 {
                Some(*seg.first_point())
            } else {
                Some(*seg.last_point())
            }
        }
    }

    /// `Status1()` (`_1.cxx:97-117`).
    pub(super) fn status1(&self) -> P2DStatus {
        if !self.is_point_of_segment1() {
            P2DStatus::Intersection
        } else {
            let i = self.point_index - self.nb_points - 1;
            if i % 2 == 0 {
                P2DStatus::SegmentFirst
            } else {
                P2DStatus::SegmentLast
            }
        }
    }

    /// `IsPointOfSegment1()` (`_1.cxx:120-124`).
    pub(super) fn is_point_of_segment1(&self) -> bool {
        self.point_index > self.nb_points
    }

    /// `Index1()` (`_1.cxx:127-130`).
    pub(super) fn index1(&self) -> i32 {
        self.point_index
    }

    /// `EdgesConfig1()` (`_1.cxx:133-147`).
    pub(super) fn edges_config1(&self) -> Config {
        let mut c = Config::UnshGeometry;
        if self.is_point_of_segment1() {
            let so = shape_tool::edges_same_oriented(&self.edge2.0, &self.edge1.0);
            c = if so {
                Config::SameOriented
            } else {
                Config::DiffOriented
            };
        }
        c
    }

    /// `Transition1(Index, EdgeOrientation)` (`_1.cxx:150-347`).
    pub(super) fn transition1(
        &self,
        index: i32,
        edge_orientation: Orientation,
    ) -> Result<Transition, String> {
        let pointofsegment = self.is_point_of_segment1();
        let pur1d = pointofsegment && self.same_domain;

        let mut sta_b = State::Unknown;
        let mut sta_a = State::Unknown;
        let mut sha_b = ShapeType::Compound;
        let mut sha_a = ShapeType::Compound;
        // pextremity = false; JYL290998 corr regr cto100K1 fex6 fex4 : 5eme
        // inters E/E (`_1.cxx:186`).
        let mut pextremity = false;

        let sta_inon = State::In;
        let dim = self.dimension;
        if dim == 1 {
            sha_a = ShapeType::Edge;
            sha_b = ShapeType::Edge;
        } else if dim == 2 && pur1d {
            sha_a = ShapeType::Edge;
            sha_b = ShapeType::Edge;
        } else if dim == 2 && !pur1d {
            sha_a = ShapeType::Face;
            sha_b = ShapeType::Face;
        }

        if edge_orientation == Orientation::Internal || edge_orientation == Orientation::External {
            let mut tr = Transition::new(sta_inon, sta_inon, sha_b, sha_a);
            tr.set_orientation(edge_orientation);
            return Ok(tr);
        }

        let ip = self
            .point1()
            .ok_or_else(|| "TopOpeBRep_EdgesIntersector : Transition1".to_string())?;
        let t = if index == 1 {
            *ip.transition_of_first()
        } else {
            *ip.transition_of_second()
        };

        match t.transition_type() {
            IntRes2dTypeTrans::In => {
                sta_b = State::Out;
                sta_a = sta_inon;
            }

            IntRes2dTypeTrans::Out => {
                sta_b = sta_inon;
                sta_a = State::Out;
            }

            IntRes2dTypeTrans::Touch => match t.situation() {
                IntRes2dSituation::Inside => {
                    sta_b = sta_inon;
                    sta_a = sta_inon;
                }

                IntRes2dSituation::Outside => {
                    sta_b = State::Out;
                    sta_a = State::Out;
                }

                IntRes2dSituation::Unknown => {
                    // get posindex = position of point on edge <Index>
                    let posindex = if index == 1 {
                        ip.transition_of_first().position_on_curve()
                    } else {
                        ip.transition_of_second().position_on_curve()
                    };

                    if pointofsegment {
                        // get posother = position of point on the other edge
                        let posother = if index == 1 {
                            ip.transition_of_second().position_on_curve()
                        } else {
                            ip.transition_of_first().position_on_curve()
                        };

                        if posother == IntRes2dPosition::Middle {
                            if posindex != IntRes2dPosition::Middle {
                                sta_b = sta_inon;
                                sta_a = sta_inon;
                            } else {
                                // Middle/Middle is impossible
                                return Err(
                                    "TopOpeBRep_EdgesIntersector : Situation Unknown MM".to_string()
                                );
                            }
                        } else {
                            // posother = Head or End
                            let opposite = self.is_opposite1()?;
                            if opposite {
                                if posother == IntRes2dPosition::Head {
                                    sta_b = sta_inon;
                                    sta_a = State::Out;
                                } else if posother == IntRes2dPosition::End {
                                    sta_b = State::Out;
                                    sta_a = sta_inon;
                                }
                            } else {
                                if posother == IntRes2dPosition::Head {
                                    sta_b = State::Out;
                                    sta_a = sta_inon;
                                } else if posother == IntRes2dPosition::End {
                                    sta_b = sta_inon;
                                    sta_a = State::Out;
                                }
                            }
                        }
                    } else {
                        // Two edges intersect on a vertex shared by both, and
                        // are tangent there (`_1.cxx:271-297`).
                        pextremity = true;
                        sha_a = ShapeType::Edge;
                        sha_b = ShapeType::Edge;

                        if posindex == IntRes2dPosition::Head {
                            sta_b = sta_inon;
                            sta_a = State::Out;
                        } else if posindex == IntRes2dPosition::End {
                            sta_b = State::Out;
                            sta_a = sta_inon;
                        } else {
                            // Middle is impossible
                            return Err(
                                "TopOpeBRep_EdgesIntersector : Situation Unknown M".to_string()
                            );
                        }
                    }
                }
            },

            IntRes2dTypeTrans::Undecided => {
                return Err(
                    "TopOpeBRep_EdgesIntersector : TransitionType Undecided".to_string()
                )
            }
        }

        let mut tr = Transition::default();
        if pur1d || pextremity {
            tr.set(sta_b, sta_a, sha_b, sha_a);
        } else {
            // +ooOO (`_1.cxx:312-342`)
            let mut composori = false;
            composori = composori || ((index == 1) && (!self.f2surf1f_sameoriented));
            composori = composori || ((index == 2) && (!self.f1surf1f_sameoriented));
            let mut eori = edge_orientation;
            if composori {
                eori = eori.reversed();
            }

            tr.set(sta_b, sta_a, sha_b, sha_a);
            if eori == Orientation::Reversed {
                tr = tr.complement();
            }
        }
        Ok(tr)
    }

    /// `Parameter1(Index)` (`_1.cxx:349-357`).
    pub(super) fn parameter1(&self, index: i32) -> Result<f64, String> {
        let p = self
            .point1()
            .ok_or_else(|| "TopOpeBRep_EdgesIntersector : Parameter1".to_string())?;
        if index == 1 {
            Ok(p.param_on_first())
        } else {
            Ok(p.param_on_second())
        }
    }

    /// `IsVertex1(Index)` (`_1.cxx:361-438`).
    pub(super) fn is_vertex1(&mut self, index: i32) -> bool {
        // Check if the last IsVertex1() call has been performed on the current
        // point and with the same <Index>.
        if self.is_vertex_point_index == self.point_index && self.is_vertex_index == index {
            return self.is_vertex_value;
        }

        // Search if the current point is a vertex of edge <Index>.
        self.is_vertex_value = false;
        let pos = self
            .point1()
            .map(|p| {
                if index == 1 {
                    p.transition_of_first().position_on_curve()
                } else {
                    p.transition_of_second().position_on_curve()
                }
            })
            .unwrap_or(IntRes2dPosition::Middle);

        if pos == IntRes2dPosition::Middle {
            // Search for an INTERNAL vertex on edge <Index> with a 2d
            // parameter <parV> equal to the current point parameter <par>.
            let par = self.parameter1(index).unwrap_or(0.0);
            let e = if index == 1 { &self.edge1 } else { &self.edge2 };
            let mut ex = Explorer::new(&e.0, ShapeType::Vertex);
            while ex.more() {
                let cur = ex.current().clone();
                ex.next();
                let Some(v) = Vertex::wrap(cur) else {
                    continue;
                };
                if v.0.orientation() == Orientation::Internal {
                    let par_v =
                        BRepTool::parameter_of_vertex_on_face(&v.0, e, &self.face1);
                    if (par - par_v).abs() <= PCONFUSION {
                        self.is_vertex_value = true;
                        self.is_vertex_vertex = Some(v);
                        break;
                    }
                }
            }
        } else {
            // pos = head or end
            let (v1, v2) = if index == 1 {
                edge_vertices(&self.edge1)
            } else {
                edge_vertices(&self.edge2)
            };
            if pos == IntRes2dPosition::Head && v1.is_some() {
                self.is_vertex_value = true;
                self.is_vertex_vertex = v1;
            } else if pos == IntRes2dPosition::End && v2.is_some() {
                self.is_vertex_value = true;
                self.is_vertex_vertex = v2;
            }
            // ... else myIsVertexValue has been set to False
        }

        // Memorize that IsVertex1() has been called on point myPointIndex and
        // on edge <Index>.
        self.is_vertex_point_index = self.point_index;
        self.is_vertex_index = index;

        self.is_vertex_value
    }

    /// `Vertex1(Index)` (`_1.cxx:442-449`).
    pub(super) fn vertex1(&mut self, index: i32) -> Result<&Vertex, String> {
        if !self.is_vertex1(index) {
            return Err("TopOpeBRep_EdgesIntersector : Vertex1".to_string());
        }
        self.is_vertex_vertex
            .as_ref()
            .ok_or_else(|| "TopOpeBRep_EdgesIntersector : Vertex1".to_string())
    }

    /// `Value1()` (`_1.cxx:454-463`).
    pub(super) fn value1(&self) -> occt_core::gp::GpPnt {
        let Some(p2) = self.point1().map(|p| *p.value()) else {
            return occt_core::gp::GpPnt::new(
                occt_core::precision::INFINITE,
                occt_core::precision::INFINITE,
                occt_core::precision::INFINITE,
            );
        };
        if Precision::is_infinite(p2.x()) || Precision::is_infinite(p2.y()) {
            occt_core::gp::GpPnt::new(
                occt_core::precision::INFINITE,
                occt_core::precision::INFINITE,
                occt_core::precision::INFINITE,
            )
        } else {
            self.surface1
                .surface()
                .map(|s| s.d0(p2.x(), p2.y()))
                .unwrap_or_else(|| occt_core::gp::GpPnt::new(0.0, 0.0, 0.0))
        }
    }
}
