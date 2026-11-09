//! Edge/face interference enumeration of `HLRBRep_Data` -
//! `InitInterference` (`HLRBRep_Data.cxx:1223-1231`), `MoreInterference`
//! (`HLRBRep_Data.lxx:103-106`), `NextInterference` (`cxx:1233-1484`),
//! `SameVertex` (`cxx:2597-2653`), `LocalLEGeometry2D` (`cxx:1532-1551`) and
//! `LocalFEGeometry2D` (`cxx:1553-1578`).
//!
//! `InitInterference` restarts the pass of `myFaceItr1` over the hiding face's
//! edges for the **current** edge `myLE` (`InitEdge` above put the face's own
//! edges at `myHideCount - 1`). `NextInterference` then walks those edges and
//! keeps the first one that produces an intersection: the loop condition
//! `!MoreInterference()` means "the current pair has nothing to report", and
//! the body fills `myNbPoints` / `myNbSegments` of the **current** pair.
//!
//! The two cache layers of `HLRBRep_Data` are used exactly as in OCCT:
//!
//! * `TableauRejection::NoIntersection(LE, FE)` skips pairs already proven
//!   disjoint (`cxx:1297`);
//! * `GetSingleIntersection(LE, FE)` replays a cached single point through
//!   `SimulateOnePoint` instead of intersecting again (`cxx:1411-1418`);
//! * a fresh single-point result is stored back (`cxx:1426-1428`), and a pair
//!   with no point and no segment is marked disjoint (`cxx:1446-1448`).
//!
//! `SameVertex` decides whether the two edges share a head/end vertex and, if
//! so, which of the four ends must be cut by `CutLar = 0.2` (`cxx:46`) before
//! intersecting - the cut is what makes two edges meeting at a vertex produce
//! a real crossing rather than a degenerate touch.
//!
//! Deviations:
//! * `myLLProps` / `myFLProps` (`HLRBRep_CLProps`, a
//!   `GeomLProp_CLPropsBase<..., HLRBRep_CLPropsATool>` alias over the 2D view
//!   of `HLRBRep_Curve`) hold a raw curve pointer in OCCT. The port builds a
//!   [`GeomLPropCLProps2d`] on demand from `my_e_data(i).geometry()`; the
//!   `SetCurve` / `SetParameter` split of `cxx:1226` / `:1535` collapses into
//!   one construction.
//! * `LocalLEGeometry2D` / `LocalFEGeometry2D` raise `Standard_Failure` when
//!   the tangent is undefined (`cxx:1537`, `:1565`); the port returns `None`.
//! * `nbOkIntersection` / `nbPtIntersection` / `nbSegIntersection` /
//!   `nbCal1Intersection` / `nbCal2Intersection` (`cxx:47-56`) are only read by
//!   `OCCT_DEBUG` prints; they are dropped.

use occt_core::gp::GpDir2d;
use occt_core::precision::REAL_EPSILON;
use occt_geom2d::cl_props2d::GeomLPropCLProps2d;

use crate::abs::Orientation;
use crate::geom_abs_curve_type::GeomAbsCurveType;

use super::HlrBRepData;
use super::edge_iter::separated;
use super::reject::REAL_LAST;

/// `static const double CutLar = 2.e-1` (`cxx:46`): the fraction of the edge
/// cut away on the side that shares a vertex with the other edge.
const CUT_LAR: f64 = 2.0e-1;

impl HlrBRepData {
    /// `InitInterference()` (`cxx:1223-1231`).
    pub fn init_interference(&mut self) {
        // `myLLProps.SetCurve(myLEGeom)` - see the module notes; the props are
        // built on demand by `local_le_geometry_2d`.
        let idx = self.i_face as usize;
        self.my_face_itr1.init_edge(&self.my_f_data[idx]);
        self.my_nb_points = 0;
        self.my_nb_segments = 0;
        self.i_interf = 0;
        self.next_interference();
    }

    /// `MoreInterference()` (`lxx:103-106`): `iInterf <= NbPoints + 2*NbSegments`.
    pub fn more_interference(&self) -> bool {
        self.i_interf <= self.my_nb_points + 2 * self.my_nb_segments
    }

    /// `NbPoints()` of the current pair (`myNbPoints`, `hxx:279`).
    pub fn nb_points(&self) -> i32 {
        self.my_nb_points
    }

    /// `NbSegments()` of the current pair (`myNbSegments`, `hxx:280`).
    pub fn nb_segments(&self) -> i32 {
        self.my_nb_segments
    }

    /// `iInterf` (`hxx:281`), the 1-based cursor `MoreInterference` reads.
    pub fn interference_index(&self) -> i32 {
        self.i_interf
    }

    /// `NextInterference()` (`cxx:1233-1484`).
    pub fn next_interference(&mut self) {
        // are there more intersections on the current edge
        self.i_interf += 1;

        while !self.more_interference() && self.my_face_itr1.more_edge() {
            // Rejection of the current wire (`cxx:1243-1268`).
            if self.my_face_itr1.beginning_of_wire() {
                let wire_min_max = self.my_face_itr1.wire().map(|w| *w.min_max());
                if let (Some(wire), Some(le)) = (wire_min_max, self.my_le_min_max) {
                    if separated(&wire, &le) {
                        self.my_face_itr1.skip_wire();
                        continue;
                    }
                }
            }

            self.my_fe = self.my_face_itr1.edge();
            self.my_fe_ori = self.my_face_itr1.orientation();
            self.my_fe_out_line = self.my_face_itr1.out_line();
            self.my_fe_internal = self.my_face_itr1.internal();
            self.my_fe_double = self.my_face_itr1.double();
            self.load_fe_data();

            let fe = self.my_fe as usize;
            if self.my_fe_ori == Orientation::Forward || self.my_fe_ori == Orientation::Reversed {
                // Edge from the boundary (`cxx:1279-1281`): not vertical and not
                // a double edge (unless it is an outline).
                if !self.my_e_data[fe].vertical() && (!self.my_fe_double || self.my_fe_out_line) {
                    let fe_min_max = *self.my_e_data[fe].min_max();
                    let Some(le_min_max) = self.my_le_min_max else {
                        self.my_face_itr1.next_edge();
                        continue;
                    };
                    // `Rejection 1/2` (`cxx:1302-1319`): the two boxes must
                    // overlap in all sixteen directions.
                    if !self.my_reject.no_intersection(self.my_le, self.my_fe)
                        && !separated(&fe_min_max, &le_min_max)
                    {
                        self.compute_interference();
                    }
                }
            }
            // next edge in face (`cxx:1480`)
            self.my_face_itr1.next_edge();
        }
    }

    /// The body of the `!rej` branch of `NextInterference` (`cxx:1329-1467`).
    fn compute_interference(&mut self) {
        let le = self.my_le;
        let fe = self.my_fe;

        let mut rej = false;
        if le == fe {
            // test if an auto-intersection is not useful (`cxx:1319-1328`).
            // OCCT writes `AutoIntersectionDone(true)` when the flag is already
            // set - reproduced verbatim.
            if self.my_e_data[le as usize].auto_intersection_done() {
                self.my_e_data[le as usize].set_auto_intersection_done(true);
                if self.my_e_data[le as usize].simple() {
                    rej = true;
                }
            }
        }
        if rej {
            return;
        }

        // `cxx:1332-1370`: which ends share a vertex.
        let mut h1 = false;
        let mut e1 = false;
        let mut h2 = false;
        let mut e2 = false;
        self.my_same_vertex = false;

        if le == fe {
            self.my_intersected = true;
            self.my_same_vertex = false;
        } else {
            self.my_intersected = true;
            if self.same_vertex(true, true) {
                self.my_same_vertex = true;
                h1 = true;
                h2 = true;
            }
            if self.same_vertex(true, false) {
                self.my_same_vertex = true;
                h1 = true;
                e2 = true;
            }
            if self.same_vertex(false, true) {
                self.my_same_vertex = true;
                e1 = true;
                h2 = true;
            }
            if self.same_vertex(false, false) {
                self.my_same_vertex = true;
                e1 = true;
                e2 = true;
            }
        }

        self.my_nb_points = 0;
        self.my_nb_segments = 0;
        self.i_interf = 1;

        if self.my_intersected {
            // `cxx:1377-1406`: cut the ends that share a vertex.
            let (mut da1, mut db1, mut da2, mut db2) = (0.0, 0.0, 0.0, 0.0);
            if self.my_same_vertex || le == fe {
                if h1 {
                    da1 = CUT_LAR;
                }
                if e1 {
                    db1 = CUT_LAR;
                }
                if h2 {
                    da2 = CUT_LAR;
                }
                if e2 {
                    db2 = CUT_LAR;
                }
            }

            // `NoInter` is never set in OCCT (`cxx:1406`, `:1430`).
            let no_inter = 0;
            if le == fe {
                self.my_intersector
                    .perform_auto(&self.my_e_data[le as usize], da1, db1);
            } else {
                let (su, sv) = self.my_reject.get_single_intersection(le, fe);
                if su != REAL_LAST {
                    self.my_intersector.simulate_one_point(
                        &self.my_e_data[le as usize],
                        su,
                        &self.my_e_data[fe as usize],
                        sv,
                    );
                } else {
                    self.my_intersector.perform(
                        &self.my_e_data[le as usize],
                        da1,
                        db1,
                        &self.my_e_data[fe as usize],
                        da2,
                        db2,
                        self.my_same_vertex,
                    );
                    if self.my_intersector.is_done()
                        && self.my_intersector.nb_points() == 1
                        && self.my_intersector.nb_segments() == 0
                    {
                        let p = self.my_intersector.point(0);
                        self.my_reject.set_intersection(le, fe, &p);
                    }
                }
            }

            if no_inter != 0 {
                self.my_nb_points = 0;
                self.my_nb_segments = 0;
            } else if self.my_intersector.is_done() {
                self.my_nb_points = self.my_intersector.nb_points() as i32;
                self.my_nb_segments = self.my_intersector.nb_segments() as i32;
                if self.my_nb_segments + self.my_nb_points == 0 {
                    self.my_reject.set_no_intersection(le, fe);
                }
            } else {
                self.my_nb_points = 0;
                self.my_nb_segments = 0;
            }
        }
    }

    /// `SameVertex(h1, h2)` (`cxx:2597-2653`): `true` when the two edge ends
    /// hold the same vertex; also updates `myIntersected` for the degenerate
    /// cases OCCT skips.
    fn same_vertex(&mut self, h1: bool, h2: bool) -> bool {
        let le = self.my_le as usize;
        let fe = self.my_fe as usize;
        let v1 = if h1 {
            self.my_e_data[le].v_sta()
        } else {
            self.my_e_data[le].v_end()
        };
        let v2 = if h2 {
            self.my_e_data[fe].v_sta()
        } else {
            self.my_e_data[fe].v_end()
        };

        let same_v = v1 == v2;
        if same_v {
            self.my_intersected = true; // compute the intersections
            let line_like = |t: GeomAbsCurveType| {
                matches!(
                    t,
                    GeomAbsCurveType::Line | GeomAbsCurveType::Circle | GeomAbsCurveType::Ellipse
                )
            };
            if line_like(self.my_le_type) && line_like(self.my_fe_type) {
                self.my_intersected = false; // no other intersection
            }

            let mut other_case = true;
            if (h1 && self.my_e_data[le].out_lv_sta()) || (!h1 && self.my_e_data[le].out_lv_end()) {
                if self.i_face_test || self.my_le_internal {
                    other_case = false;
                }
            } else if self.i_face_test {
                other_case = false;
            }

            if other_case {
                if (h1 && self.my_e_data[le].cut_at_sta())
                    || (!h1 && self.my_e_data[le].cut_at_end())
                {
                    // two connected OutLines do not intersect themselves
                    self.my_intersected = false;
                }
            }
        }
        same_v
    }

    /// `LocalLEGeometry2D(Param, Tg, Nm, Cu)` (`cxx:1532-1551`): the tangent,
    /// normal and curvature of the projected **current** edge at `Param`.
    ///
    /// `Cu > Epsilon(1.)` and finite selects the normal of the curve, otherwise
    /// the left normal of the tangent (`cxx:1541-1548`).
    pub fn local_le_geometry_2d(&self, param: f64) -> Option<(GpDir2d, GpDir2d, f64)> {
        local_geometry_2d(self.my_e_data[self.my_le as usize].geometry(), param)
    }

    /// `LocalFEGeometry2D(FE, Param, Tg, Nm, Cu)` (`cxx:1553-1578`).
    pub fn local_fe_geometry_2d(
        &self,
        fe: i32,
        param: f64,
    ) -> Option<(GpDir2d, GpDir2d, f64)> {
        local_geometry_2d(self.my_e_data[fe as usize].geometry(), param)
    }
}

/// The shared body of `LocalLEGeometry2D` / `LocalFEGeometry2D`
/// (`cxx:1532-1578`), over the 2D view of an `HLRBRep_Curve`.
fn local_geometry_2d(
    curve: &crate::hlrbrep_curve::HlrBRepCurve,
    param: f64,
) -> Option<(GpDir2d, GpDir2d, f64)> {
    let mut props = GeomLPropCLProps2d::new(curve, param, 2, REAL_EPSILON);
    if !props.is_tangent_defined() {
        return None;
    }
    let tg = props.tangent();
    let cu = props.curvature();
    let nm = if cu > occt_core::precision::epsilon(1.0)
        && !occt_core::precision::Precision::is_infinite(cu)
    {
        props.normal()
    } else {
        // `Nm = gp_Dir2d(-Tg.Y(), Tg.X())` (`cxx:1546`): the coordinates are
        // written as they are, without the `gp_Dir2d` unit check.
        GpDir2d::new(-tg.y, tg.x).unwrap_or(GpDir2d { x: -tg.y, y: tg.x })
    };
    Some((tg, nm, cu))
}
