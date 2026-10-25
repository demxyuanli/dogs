//! `TopOpeBRep_EdgesIntersector::Perform` and the point filter it applies
//! before building the `TopOpeBRep_Point2d` sequence.
//!
//! Source: `TKBool/TopOpeBRep/TopOpeBRep_EdgesIntersector.cxx:187-688`.

use occt_core::intres2d::{
    IntRes2dDomain, IntRes2dIntersectionPoint, IntRes2dPosition, IntRes2dTransition,
    IntRes2dTypeTrans,
};
use occt_geom2d::curve::Curve2d;

use crate::abs::{Orientation, ShapeType, State};
use crate::brep_tool::BRepTool;
use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::is_same;
use crate::topopebrep::ds::kind::Config;
use crate::topopebrep::ds::transition::Transition;
use crate::topopebrep::interf::hctxee2d::project_edge_on_face;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;
use crate::topopebrep::tool::closed::closed_e;
use crate::topopebrep::tool::shape_tool;
use crate::topopebrep::tool::topology::{fun_tool_orient_vin_e, CLOSING, FIRST, LAST};

use super::{same_face, same_support, uv_points, EdgesIntersector};

/// `TransitionEqualAndExtremity(T1, T2)` (`cxx:187-217`).
pub(super) fn transition_equal_and_extremity(
    t1: &IntRes2dTransition,
    t2: &IntRes2dTransition,
) -> bool {
    if t1.position_on_curve() == IntRes2dPosition::Head
        || t1.position_on_curve() == IntRes2dPosition::End
    {
        if t1.position_on_curve() == t2.position_on_curve() {
            if t1.transition_type() == t2.transition_type() {
                if t1.transition_type() == IntRes2dTypeTrans::Touch {
                    if t1.is_tangent() == t2.is_tangent() {
                        if t1.situation() == t2.situation() {
                            if t1.is_opposite() == t2.is_opposite() {
                                return true;
                            }
                        }
                    }
                } else {
                    return true;
                }
            }
        }
    }
    false
}

/// `IsTangentSegment(P1, P2, aC1, aC2, aTolConf)` (`cxx:219-249`).
pub(super) fn is_tangent_segment(
    p1: &IntRes2dIntersectionPoint,
    p2: &IntRes2dIntersectionPoint,
    c1: &dyn Curve2d,
    c2: &dyn Curve2d,
    tol_conf: f64,
) -> bool {
    let trans1 = p1.transition_of_first();
    let trans2 = p2.transition_of_first();

    if trans1.transition_type() == IntRes2dTypeTrans::Touch
        || trans2.transition_type() == IntRes2dTypeTrans::Touch
    {
        let sqr_dist_pp = p1.value().square_distance(p2.value());
        if sqr_dist_pp <= tol_conf {
            let par_dist1 = (p1.param_on_first() - p2.param_on_first()).abs();
            let par_dist2 = (p1.param_on_second() - p2.param_on_second()).abs();
            let resol1 = c1.resolution(tol_conf);
            let resol2 = c2.resolution(tol_conf);
            if par_dist1 * par_dist1 <= resol1 && par_dist2 * par_dist2 <= resol2 {
                return true;
            }
        }
    }
    false
}

/// `EdgesIntersector_checkT1D(E1, E2, vG, newT)` (`cxx:253-309`): `newT` is the
/// transition attached to `E1` at the vertex `vG` of `E2`.
pub(super) fn check_t1d(
    e1: &TopoShape,
    e2: &TopoShape,
    vg: &TopoShape,
) -> (bool, Transition) {
    let mut new_t = Transition::default();
    new_t.set(
        State::Unknown,
        State::Unknown,
        crate::abs::ShapeType::Face,
        crate::abs::ShapeType::Face,
    );

    let ovine = fun_tool_orient_vin_e(vg, e2);
    if ovine == 0 {
        return (false, new_t);
    } else if ovine == CLOSING {
        new_t.set_orientation(Orientation::Internal);
        return (true, new_t);
    }

    let first = ovine == FIRST;
    let last = ovine == LAST;

    let mut c = Config::SameOriented;
    let sso = shape_tool::shapes_same_oriented(e1, e2);
    if !sso {
        c = Config::DiffOriented;
    }

    let mut so = c == Config::SameOriented;
    let mut d_o = c == Config::DiffOriented;
    let o1 = e1.orientation();
    if o1 == Orientation::Reversed {
        so = !so;
        d_o = !d_o;
    } // relative to E1 FORWARD

    let reversed = (so && first) || (d_o && last);
    let forward = (so && last) || (d_o && first);
    if reversed {
        new_t.set_orientation(Orientation::Reversed);
    }
    if forward {
        new_t.set_orientation(Orientation::Forward);
    }
    (reversed || forward, new_t)
}

/// `FC2D_CurveOnSurface` (`2d.cxx:356-376`): the stored p-curve when there is
/// one, otherwise the `FC2D_make2d` projection (`2d.cxx:279-305`).
fn fc2d_curve_on_surface(e: &Edge, f: &Face) -> Option<(std::sync::Arc<dyn Curve2d>, f64, f64)> {
    if let Some(r) = crate::boptools_2d::curve_on_surface_range(e, f) {
        return Some(r);
    }
    let (first, last) = BRepTool::edge_parameters(e);
    let pc = crate::pcurve::make_pcurve_on_face(e, f).ok()?;
    Some((pc, first, last))
}

/// `FC2D_MakeCurveOnSurface` (`2d.cxx:342-353`): `FC2D_make2d` followed by
/// `FC2D_AddNewCurveOnSurface` (`2d.cxx:236-254`), i.e. the built p-curve is
/// registered on the edge.
fn fc2d_make_curve_on_surface(
    e: &Edge,
    f: &Face,
) -> Option<(std::sync::Arc<dyn Curve2d>, f64, f64)> {
    if let Some(r) = crate::boptools_2d::curve_on_surface_range(e, f) {
        return Some(r);
    }
    let (first, last) = BRepTool::edge_parameters(e);
    let pc = crate::boptools_2d::build_pcurve_for_edge_on_face(e, f).ok()?;
    Some((pc, first, last))
}

impl EdgesIntersector {
    /// `Perform(E1, E2)` (`hxx`, `ReduceSegment = Standard_True`).
    pub fn perform(&mut self, e1: &TopoShape, e2: &TopoShape) {
        self.perform_with(e1, e2, true);
    }

    /// `Perform(E1, E2, ReduceSegment)` (`cxx:316-688`).
    ///
    /// The two `Standard_Failure` raises of `cxx:346-348` (`Perform : no 2d
    /// curve`) and of the `Transition1` switch are panics here, as elsewhere in
    /// the port; the `return;` arms of the projection (`cxx:410-412`,
    /// `cxx:492-494`) leave the object empty.
    pub fn perform_with(&mut self, e1: &TopoShape, e2: &TopoShape, reduce_segment: bool) {
        self.sp2d.clear();
        self.ip2d = 1;
        self.np2d = 0;

        self.edge1 = Edge::wrap(e1.clone()).unwrap_or_else(Edge::new);
        self.edge2 = Edge::wrap(e2.clone()).unwrap_or_else(Edge::new);

        // `BRepAdaptor_Surface aSurface1(myFace1), aSurface2(myFace2)`
        // (`cxx:333-334`): fresh adaptors, only their `GetType()` is used.
        let a_surface1 = SurfaceAdaptor::of_face(&self.face1);
        let a_surface2 = SurfaceAdaptor::of_face(&self.face2);
        let a_surface_type1 = a_surface1.get_type();
        let a_surface_type2 = a_surface2.get_type();

        // `cxx:336-345`: a sphere/sphere pair forces a new p-curve
        // (`FC2D_MakeCurveOnSurface`, `trim3d = true`), otherwise the stored
        // one is reused (`FC2D_CurveOnSurface`).
        let pc1 = if a_surface_type1 == GeomAbsSurfaceType::Sphere
            && a_surface_type2 == GeomAbsSurfaceType::Sphere
        {
            fc2d_make_curve_on_surface(&self.edge1, &self.face1)
        } else {
            fc2d_curve_on_surface(&self.edge1, &self.face1)
        };
        let Some((pc1, first, last)) = pc1 else {
            panic!("EdgesIntersector::Perform : no 2d curve");
        };

        self.curve1 = Some(pc1);
        let (pfirst, plast) = uv_points(&self.edge1, &self.face1);
        let tole = BRepTool::edge_tolerance(&self.edge1);
        self.domain1 = IntRes2dDomain::bounded(&pfirst, first, tole, &plast, last, tole);

        let memesfaces = same_face(&self.face1, &self.face2);
        let memesupport = same_support(&self.face1, &self.face2);

        if self.surface_type1 == GeomAbsSurfaceType::Plane || memesfaces || memesupport {
            // `cxx:376-393`: edge 2 read on face 1 directly.
            let Some((pc2, first2, last2)) = fc2d_curve_on_surface(&self.edge2, &self.face1)
            else {
                return;
            };
            self.curve2 = Some(pc2);
            let (pfirst, plast) = uv_points(&self.edge2, &self.face1);
            let tole = BRepTool::edge_tolerance(&self.edge2);
            self.domain2 = IntRes2dDomain::bounded(&pfirst, first2, tole, &plast, last2, tole);
        } else {
            // `cxx:397-495`: project the curve of edge 2 on the surface of
            // face 1. `first` / `last` are by-reference in OCCT and the
            // projection arms overwrite them.
            let mut first = first;
            let mut last = last;
            let pc2on1 = match project_edge_on_face(
                &self.edge1,
                &self.edge2,
                &self.face1,
                &mut first,
                &mut last,
            ) {
                Ok(pc) => pc,
                Err(_) => return,
            };
            let Some(pc2on1) = pc2on1 else {
                return;
            };
            self.curve2 = Some(pc2on1.clone());
            let tole = BRepTool::edge_tolerance(&self.edge2);
            let pfirst = pc2on1.d0(first);
            let plast = pc2on1.d0(last);
            self.domain2 = IntRes2dDomain::bounded(&pfirst, first, tole, &plast, last, tole);
        }

        // `cxx:528`: myIntersector.Perform(myCurve1, myDomain1, myCurve2,
        // myDomain2, tol1, tol2).
        let tol1 = self.tol1;
        let tol2 = self.tol2;
        {
            let (Some(c1), Some(c2)) = (self.curve1.clone(), self.curve2.clone()) else {
                return;
            };
            self.intersector
                .perform(c1.as_ref(), &self.domain1, c2.as_ref(), &self.domain2, tol1, tol2);
        }

        let nbp = self.intersector.nb_points() as i32;
        let nbs = self.intersector.nb_segments() as i32;

        self.lpnt.clear();
        self.lseg.clear();
        let mut p = 1;
        while p <= nbp {
            self.lpnt.push(self.intersector.point(p as usize));
            p += 1;
        }
        let mut s = 1;
        while s <= nbs {
            self.lseg.push(self.intersector.segment(s as usize));
            s += 1;
        }

        let mut nbp = nbp;

        // `cxx:546-596`: the filter. `OCCT_DEBUG`'s `nofilter` switch has no
        // release-build counterpart.
        loop {
            let mut fin = true;
            let mut p = 1;
            while p < nbp && fin {
                let p1 = self.lpnt[(p - 1) as usize];
                let p2 = self.lpnt[p as usize];
                if transition_equal_and_extremity(p1.transition_of_first(), p2.transition_of_first())
                    || transition_equal_and_extremity(
                        p1.transition_of_second(),
                        p2.transition_of_second(),
                    )
                {
                    fin = false;
                    self.lpnt.remove((p - 1) as usize);
                    nbp -= 1;
                } else {
                    let (Some(c1), Some(c2)) = (self.curve1.clone(), self.curve2.clone()) else {
                        return;
                    };
                    if is_tangent_segment(
                        &p1,
                        &p2,
                        c1.as_ref(),
                        c2.as_ref(),
                        tol1.max(tol2),
                    ) {
                        let a_trans = p2.transition_of_first();
                        fin = false;
                        if a_trans.transition_type() == IntRes2dTypeTrans::Touch {
                            self.lpnt.remove((p - 1) as usize);
                        } else {
                            self.lpnt.remove(p as usize);
                        }
                        nbp -= 1;
                    }
                }
                p += 1;
            }
            if fin {
                break;
            }
        }

        // `cxx:598-605`.
        self.nb_segments = self.lseg.len() as i32;
        self.has_segment = self.nb_segments != 0;
        self.compute_same_domain();

        self.nb_points = self.lpnt.len() as i32;
        self.true_nb_points = self.nb_points + 2 * self.nb_segments;
        self.point_index = 0;

        self.make_points2d();
        if reduce_segment {
            self.reduce_segments();
        }

        // `cxx:619-682`: the same-domain transition fix-up.
        let esd = self.same_domain();
        self.init_point(true);
        while self.more_point() {
            let isvertex1 = self
                .point()
                .and_then(|p| p.is_vertex(1))
                .unwrap_or(false);
            let isvertex2 = self
                .point()
                .and_then(|p| p.is_vertex(2))
                .unwrap_or(false);
            let isvertex = isvertex1 || isvertex2;

            if isvertex && esd {
                let isvertex12 = isvertex1 && isvertex2;
                let isvertex22 = isvertex2 && !isvertex12;
                let isvertex11 = isvertex1 && !isvertex12;

                let t1_internal = self.is_transition_internal(1);
                let t2_internal = self.is_transition_internal(2);
                let t2_external = self.is_transition_external(2);

                if t1_internal && isvertex2 && !isvertex1 {
                    // `const TopoDS_Vertex& V2 = P2D.Vertex(2)` (`cxx:641`).
                    let v2 = self.point().and_then(|p| p.vertex(2)).cloned();
                    if let Ok(v2) = v2 {
                        let (computed, new_t) = check_t1d(&self.edge1.0, &self.edge2.0, &v2.0);
                        if computed {
                            let o = new_t.orientation(State::In, ShapeType::Face);
                            let _ = self
                                .point_mut()
                                .change_transition(1)
                                .map(|t| t.set_orientation(o));
                        }
                    }
                }

                let int_ext2 = t2_internal || t2_external;
                if int_ext2 && isvertex1 && !isvertex2 {
                    let v1 = self.point().and_then(|p| p.vertex(1)).cloned();
                    if let Ok(v1) = v1 {
                        let (computed, new_t) = check_t1d(&self.edge2.0, &self.edge1.0, &v1.0);
                        if computed {
                            let o = new_t.orientation(State::In, ShapeType::Face);
                            let _ = self
                                .point_mut()
                                .change_transition(2)
                                .map(|t| t.set_orientation(o));
                        }
                    }
                }

                // xpu121098 : cto900I7 (e12on,vG14) (`cxx:664-678`).
                let n_t1 = !t1_internal
                    && self.vertex_is_closed_end(&self.edge2.0, isvertex22, 2);
                if n_t1 {
                    let _ = self
                        .point_mut()
                        .change_transition(1)
                        .map(|t| t.set_orientation(Orientation::Internal));
                }

                let n_t2 = !t2_internal
                    && self.vertex_is_closed_end(&self.edge1.0, isvertex11, 1);
                if n_t2 {
                    let _ = self
                        .point_mut()
                        .change_transition(2)
                        .map(|t| t.set_orientation(Orientation::Internal));
                }
            }

            self.next_point();
        }
    }

    /// `clE && isvertexNN && vcl.IsSame(P2D.Vertex(Index))` of `cxx:664-678`.
    fn vertex_is_closed_end(&self, e: &TopoShape, isvertex_nn: bool, index: i32) -> bool {
        let Some(vcl) = closed_e(e) else {
            return false;
        };
        isvertex_nn
            && self
                .point()
                .and_then(|p| p.vertex(index))
                .map(|v| is_same(&v.0, &vcl))
                .unwrap_or(false)
    }

    /// `P2D.Transition(Index).Orientation(TopAbs_IN) == TopAbs_INTERNAL` of the
    /// `cxx:637-650` region.
    fn is_transition_internal(&self, index: i32) -> bool {
        self.point()
            .and_then(|p| {
                p.transition(index).map(|t| {
                    t.orientation(State::In, ShapeType::Face) == Orientation::Internal
                })
            })
            .unwrap_or(false)
    }

    /// `P2D.Transition(Index).Orientation(TopAbs_IN) == TopAbs_EXTERNAL`
    /// (`cxx:651`).
    fn is_transition_external(&self, index: i32) -> bool {
        self.point()
            .and_then(|p| {
                p.transition(index).map(|t| {
                    t.orientation(State::In, ShapeType::Face) == Orientation::External
                })
            })
            .unwrap_or(false)
    }
}
