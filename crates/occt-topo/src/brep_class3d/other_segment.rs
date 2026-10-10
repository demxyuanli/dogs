//! Port of `BRepClass3d_SolidExplorer::OtherSegment`
//! (`BRepClass3d_SolidExplorer.cxx:493-787`) with its helpers:
//! `ClassifyUVPoint` (`cxx:221-237`), `PointInTheFace` (`cxx:241-425`) and
//! `FindAPointInTheFace` (`cxx:74-190`).
//!
//! Control flow and statement order follow the OCCT functions. The edge branch
//! of the `myMapEV` ON test is the shared `Extrema_ExtPC` predicate
//! (`edge_accepts_point`, `brep_class3d.rs`). The `BRepClass3d_BndBoxTree`
//! pruning is replaced by a linear scan of the same `Accept` predicate.

use occt_core::gp::{GpDir, GpDir2d, GpLin, GpLin2d, GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_core::precision::{Precision, CONFUSION, INFINITE, PCONFUSION};
use occt_geom::extrema_surf::ExtPs;
use occt_geom::Surface;

use crate::abs::Orientation;
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_class::{BRepClassEdge, BRepClassFacePassiveClassifier, FaceClassifier};
use crate::brep_surface::face_uv_bounds;
use crate::brep_tool::BRepTool;
use crate::fclass2d::{FaceState, FClass2d};
use crate::int_curves_face::FaceIntersector;
use crate::shape::{Edge, Face};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::{edges_of_wire, wires_of_face};

use super::{Segment, SolidExplorer};

/// `gp::Resolution()` (`RealSmall`).
const RESOLUTION: f64 = f64::MIN_POSITIVE;

/// `IsInfiniteUV` (`BRepClass3d_SolidExplorer.cxx:454-476`): bit mask of the
/// infinite UV bounds (1: U1, 2: V1, 4: U2, 8: V2).
fn is_infinite_uv(u1: f64, v1: f64, u2: f64, v2: f64) -> i32 {
    let mut val = 0;
    if Precision::is_infinite(u1) {
        val |= 1;
    }
    if Precision::is_infinite(v1) {
        val |= 2;
    }
    if Precision::is_infinite(u2) {
        val |= 4;
    }
    if Precision::is_infinite(v2) {
        val |= 8;
    }
    val
}

/// `gp_Lin(P, V)`: `gp_Dir(V)` throws on a null vector in OCCT, so it panics here.
fn line_through(p: &GpPnt, v: &GpVec) -> GpLin {
    let dir = GpDir::from_vec(v).expect("gp_Lin(P, V): null direction vector");
    GpLin::from_pnt_dir(*p, dir)
}

/// `Parameter` step clamp used by `PointInTheFace` (`cxx:257-264`, `cxx:371-378`).
fn clamp_step(step: f64) -> f64 {
    if step < 1e-12 {
        1e-12
    } else {
        step
    }
}

/// `TopoDS_Shape::operator==` (`cxx:127`, `TopoDS_Shape.hxx:276-280`): same
/// TShape, same location and same orientation.
fn is_same_edge_occurrence(a: &Edge, b: &Edge) -> bool {
    a.0.same_tshape(&b.0)
        && a.0.location().is_equal(b.0.location())
        && a.0.orientation() == b.0.orientation()
}

/// `FClassifier.Compare(AEdge, Or)` followed by the `ClosestIntersection`
/// update of `ParamInit` / `APointExist` (`cxx:129-138` and `cxx:144-153`).
fn passive_compare(
    fc: &mut BRepClassFacePassiveClassifier,
    edge: &Edge,
    face: &Face,
    param_init: &mut f64,
    apoint_exist: &mut bool,
) {
    let ae = BRepClassEdge::from_edge_face(edge.clone(), face.clone());
    fc.compare(&ae, edge.0.orientation());
    if fc.closest_intersection() != 0 && *param_init > fc.parameter() {
        *param_init = fc.parameter();
        *apoint_exist = true;
    }
}

/// Port of `BRepClass3d_SolidExplorer::FindAPointInTheFace`
/// (`BRepClass3d_SolidExplorer.cxx:74-190`). For each edge of the forward face a
/// probe starts at the edge point at `param` and runs along the inward normal of
/// the edge tangent; the nearest hit over the other edges comes from
/// `BRepClass_FacePassiveClassifier`. The probe is shrunk by `0.41234` until it
/// is IN the face with a regular surface normal.
///
/// `apoint`, `u_`, `v_` and the derivatives are written where OCCT writes them,
/// including on a `false` return, because the caller reads them back.
pub(super) fn find_a_point_in_the_face(
    face: &Face,
    param: f64,
    apoint: &mut GpPnt,
    u_: &mut f64,
    v_: &mut f64,
    d1u: &mut GpVec,
    d1v: &mut GpVec,
) -> bool {
    const TOL_INIT: f64 = 0.00001;
    let face = Face(face.0.oriented(Orientation::Forward));
    let Some(surf) = BRepTool::face_surface(&face) else {
        return false;
    };
    // `BRepTopAdaptor_FClass2d Classifier(face, Confusion)` is built on first use.
    let mut classifier: Option<FClass2d> = None;
    // `TopExp_Explorer(face, TopAbs_EDGE)`: cumulated orientations, duplicates kept.
    let edges: Vec<Edge> = wires_of_face(&face)
        .iter()
        .flat_map(|w| edges_of_wire(w))
        .collect();
    let nb_edges = edges.len();
    for edge in &edges {
        // `BRepAdaptor_Curve2d c(Edge, face)`; a null curve is skipped (`cxx:91-95`).
        let Some((pcurve, first, last)) = curve_on_surface_range(edge, &face) else {
            continue;
        };
        let (p0, t) = pcurve.d1((last - first) * param + first);
        // `cxx:101-108`: rotate the tangent by +90 deg for FORWARD, else -90 deg.
        let rot = if edge.0.orientation() == Orientation::Forward {
            GpVec2d::new(-t.y(), t.x())
        } else {
            GpVec2d::new(t.y(), -t.x())
        };
        // `T.Normalize()` throws on a null vector (`cxx:115`); a failure ends the search.
        let Ok(tang) = rot.normalized() else {
            return false;
        };
        let p = GpPnt2d::new(
            p0.x() + TOL_INIT * tang.x(),
            p0.y() + TOL_INIT * tang.y(),
        );
        let Ok(dir) = GpDir2d::from_vec2d(&tang) else {
            return false;
        };
        let lin = GpLin2d::from_pnt_dir(p, dir);

        let mut param_init = INFINITE;
        let mut apoint_exist = false;
        let mut fc = BRepClassFacePassiveClassifier::new();
        // `FClassifier.Reset(gp_Lin2d(P, T), ParamInit, RealEpsilon())` (`cxx:118`).
        fc.reset(&lin, param_init, f64::EPSILON);

        for other in &edges {
            // `OtherEdge.Orientation() != TopAbs_EXTERNAL && OtherEdge != Edge` (`cxx:127`).
            if other.0.orientation() == Orientation::External
                || is_same_edge_occurrence(other, edge)
            {
                continue;
            }
            passive_compare(&mut fc, other, &face, &mut param_init, &mut apoint_exist);
        }
        // `if (aNbEdges == 1)` (`cxx:142-154`).
        if nb_edges == 1 {
            passive_compare(&mut fc, edge, &face, &mut param_init, &mut apoint_exist);
        }

        while apoint_exist {
            param_init *= 0.41234;
            *u_ = p.x() + param_init * tang.x();
            *v_ = p.y() + param_init * tang.y();
            if classifier.is_none() {
                classifier = FClass2d::new(&face, CONFUSION).ok();
            }
            let Some(cl) = classifier.as_ref() else {
                return false;
            };
            // `BRepTopAdaptor_FClass2d::Perform` (`cxx:163-169`), TabOrien variant.
            if cl.perform_tab_orien(GpPnt2d::new(*u_, *v_)) != FaceState::In {
                return false;
            }
            let (pnt, du, dv) = surf.d1(*u_, *v_);
            *apoint = pnt;
            let cross_mag = du.clone().cross(&dv).magnitude();
            *d1u = du;
            *d1v = dv;
            // `theVecD1U.CrossMagnitude(theVecD1V) > gp::Resolution()` (`cxx:175`).
            if cross_mag > RESOLUTION {
                return true;
            }
            if param_init < PCONFUSION {
                return false;
            }
        }
    }
    false
}

/// Grid order of `BRepClass3d_SolidExplorer::PointInTheFace` (`cxx:299-409`).
/// `visit` returns `true` to stop; the return value tells whether it stopped.
fn grid_scan(
    bounds: (f64, f64, f64, f64),
    steps: (f64, f64),
    visit: &mut dyn FnMut(f64, f64) -> bool,
) -> bool {
    let (u1, u2, v1, v2) = bounds;
    let (du, dv) = steps;
    let (mu, mv) = ((u1 + u2) * 0.5, (v1 + v2) * 0.5);
    // cxx:299-315: u increases, v increases.
    let mut u = du + mu;
    while u < u2 {
        let mut v = dv + mv;
        while v < v2 {
            if visit(u, v) {
                return true;
            }
            v += dv;
        }
        u += du;
    }
    // cxx:317-333: u decreases, v decreases.
    let mut u = -du + mu;
    while u > u1 {
        let mut v = -dv + mv;
        while v > v1 {
            if visit(u, v) {
                return true;
            }
            v -= dv;
        }
        u -= du;
    }
    // cxx:334-350: u decreases, v increases.
    let mut u = -du + mu;
    while u > u1 {
        let mut v = dv + mv;
        while v < v2 {
            if visit(u, v) {
                return true;
            }
            v += dv;
        }
        u -= du;
    }
    // cxx:351-367: u increases, v decreases.
    let mut u = du + mu;
    while u < u2 {
        let mut v = -dv + mv;
        while v > v1 {
            if visit(u, v) {
                return true;
            }
            v -= dv;
        }
        u += du;
    }
    // cxx:369-396: remainder lattice of 37 steps from the lower corner.
    let du37 = clamp_step((u2 - u1) / 37.0);
    let dv37 = clamp_step((v2 - v1) / 37.0);
    let mut u = du37 + u1;
    while u < u2 {
        let mut v = dv37 + v1;
        while v < v2 {
            if visit(u, v) {
                return true;
            }
            v += dv37;
        }
        u += du37;
    }
    // cxx:397-409: centre of the box.
    visit(mu, mv)
}

impl SolidExplorer {
    /// `BRepClass3d_SolidExplorer::ClassifyUVPoint` (`cxx:221-237`).
    fn classify_uv_point(
        &self,
        inter: &FaceIntersector,
        surf: &dyn Surface,
        uv: GpPnt2d,
    ) -> FaceState {
        let p3d = surf.d0(uv.x(), uv.y());
        if self.map_ev_selects(&p3d) {
            return FaceState::On;
        }
        inter.classify_uv_point(uv)
    }

    /// `myTree.Select(BRepClass3d_BndBoxTreeSelectorPoint)` (`cxx:228-231`):
    /// the `Accept` predicate of `BRepClass3d_BndBoxTree.cxx:25-71`, see
    /// `map_ev_accepts_point`.
    fn map_ev_selects(&self, p: &GpPnt) -> bool {
        super::map_ev_accepts_point(&self.map_ev, p)
    }

    /// `BRepClass3d_SolidExplorer::PointInTheFace` (`cxx:241-425`): grid search
    /// over the UV box for an IN sample, resuming at `index_point`. When no sample
    /// is accepted it falls back to `FindAPointInTheFace` (`cxx:418`).
    #[allow(clippy::too_many_arguments)]
    fn point_in_the_face(
        &self,
        face: &Face,
        apoint: &mut GpPnt,
        u_: &mut f64,
        v_: &mut f64,
        param_: f64,
        index_point: &mut i32,
        surf: &dyn Surface,
        bounds: (f64, f64, f64, f64),
        d1u: &mut GpVec,
        d1v: &mut GpVec,
    ) -> bool {
        let (u1, u2, v1, v2) = bounds;
        let du = clamp_step((u2 - u1) / 6.0);
        let dv = clamp_step((v2 - v1) / 6.0);
        // `aSE.Intersector(face)` (`myMapOfInter.Find`); every face of the
        // solid is bound in `load`, so the lookup always succeeds.
        let Some(inter) = self.inters.get(&GeometryRegistry::shape_key(&face.0)) else {
            return false;
        };

        // `IsInside` and the sample at the current (u_, v_) (`cxx:270-291`).
        let mut is_inside = true;
        if !surf.is_u_periodic() {
            is_inside = *u_ >= u1 && *u_ <= u2;
        }
        if !surf.is_v_periodic() {
            is_inside &= *v_ >= v1 && *v_ <= v2;
        }
        if is_inside
            && self.classify_uv_point(inter, surf, GpPnt2d::new(*u_, *v_)) == FaceState::In
        {
            let (pnt, du_v, dv_v) = surf.d1(*u_, *v_);
            *d1u = du_v;
            *d1v = dv_v;
            if pnt.square_distance(apoint) < CONFUSION * CONFUSION {
                return true;
            }
        }

        let mut nb_calc: i32 = 0;
        let hit = grid_scan(bounds, (du, dv), &mut |u, v| {
            nb_calc += 1;
            if nb_calc < *index_point {
                return false;
            }
            if self.classify_uv_point(inter, surf, GpPnt2d::new(u, v)) != FaceState::In {
                return false;
            }
            *u_ = u;
            *v_ = v;
            let (pnt, a, b) = surf.d1(u, v);
            *apoint = pnt;
            *d1u = a;
            *d1v = b;
            true
        });
        // `IndexPoint = NbPntCalc` (`cxx:411`).
        *index_point = nb_calc;
        if hit {
            return true;
        }
        find_a_point_in_the_face(face, param_, apoint, u_, v_, d1u, d1v)
    }

    /// `BRepClass3d_SolidExplorer::OtherSegment` (`cxx:493-787`). Returns the
    /// `Segment` record; the integer result of OCCT maps to `Segment::flag`.
    /// `myFirstFace`, `myParamOnEdge` and the per-face intersectors are the
    /// explorer state, as in OCCT.
    pub(super) fn other_segment(&mut self, p: &GpPnt) -> Option<Segment> {
        let tol_u = PCONFUSION;
        let tol_v = tol_u;
        let mut apoint = GpPnt::new(0.0, 0.0, 0.0);
        let mut vec_d1u = GpVec::new(0.0, 0.0, 0.0);
        let mut vec_d1v = GpVec::new(0.0, 0.0, 0.0);
        let mut max_scal = 0.0_f64;
        let mut ptfound = false;
        let mut out_par = 0.0_f64;
        let mut u_: f64;
        let mut v_: f64;
        let mut index_point: i32 = 0;
        let mut nb_points_ok: i32 = 0;
        let mut a_restr = true;
        let mut a_test_invert = false;
        let mut out_lin = GpLin::from_pnt_dir(*p, GpDir::default_dir());

        loop {
            self.first_face += 1;
            let n_faces = self.faces.len();
            let mut nb_faces_in_solid: i32 = 0;
            for k in 0..n_faces {
                nb_faces_in_solid += 1;
                if self.first_face > nb_faces_in_solid {
                    continue;
                }
                let mut face = self.faces[k].clone();

                if a_test_invert {
                    // `BRepTopAdaptor_FClass2d aClass(face, Confusion)` (`cxx:533-534`).
                    let inf_in = FClass2d::new(&face, CONFUSION)
                        .map(|c| c.perform_infinite_point_tab_orien() == FaceState::In)
                        .unwrap_or(false);
                    if inf_in {
                        a_restr = false;
                        // `myMapOfInter` rebind (`cxx:537-542`).
                        self.inters.insert(
                            GeometryRegistry::shape_key(&face.0),
                            FaceIntersector::new(face.clone(), CONFUSION, false, false),
                        );
                    } else {
                        a_restr = true;
                    }
                }
                let Some(surf) = BRepTool::face_surface(&face) else {
                    continue;
                };
                // `BRepAdaptor_Surface::Initialize(face, aRestr)` (`cxx:550-555`).
                let (u1, u2, v1, v2) = if a_restr {
                    face_uv_bounds(&face)
                } else {
                    let (a, b) = surf.u_range();
                    let (c, d) = surf.v_range();
                    (a, b, c, d)
                };
                face = Face(face.0.oriented(Orientation::Forward));

                // Degenerate UV box (`cxx:559-565`).
                let eps_u = (PCONFUSION * u2.abs().max(u1.abs())).max(PCONFUSION);
                let eps_v = (PCONFUSION * v2.abs().max(v1.abs())).max(PCONFUSION);
                if (u2 - u1).abs() < eps_u || (v2 - v1).abs() < eps_v {
                    return Some(Segment {
                        lin: out_lin,
                        par: out_par,
                        flag: 2,
                    });
                }

                let svmyparam = self.param_on_edge;
                let inf_flag = is_infinite_uv(u1, v1, u2, v2);
                u_ = (u1 + u2) * 0.5;
                v_ = (v1 + v2) * 0.5;

                let ext = ExtPs::with_surface(p, surf.as_ref(), tol_u, tol_v);
                if ext.is_done() && ext.nb_ext() > 0 {
                    let nb_ext = ext.nb_ext();
                    let mut i_near: usize = 1;
                    let mut dist2_min = ext.square_distance(1);
                    for i in 2..=nb_ext {
                        let (au, av, _) = ext.point(i);
                        if au >= u1 && au <= u2 && av >= v1 && av <= v2 {
                            let d2 = ext.square_distance(i);
                            if d2 < dist2_min {
                                dist2_min = d2;
                                i_near = i;
                            }
                        }
                    }
                    // `aDist2Tresh` (`cxx:602-604`).
                    if dist2_min < 1.0e-24 {
                        if inf_flag != 0 {
                            return Some(Segment {
                                lin: out_lin,
                                par: out_par,
                                flag: 1,
                            });
                        }
                        let (au, av, _) = ext.point(i_near);
                        let mut classifier2d = FaceClassifier::new();
                        classifier2d.perform(&face, GpPnt2d::new(au, av), PCONFUSION);
                        let flag = match classifier2d.state() {
                            FaceState::In | FaceState::On => 1,
                            _ => 3,
                        };
                        return Some(Segment {
                            lin: out_lin,
                            par: out_par,
                            flag,
                        });
                    }
                    if inf_flag != 0 {
                        let (_, _, ap) = ext.point(i_near);
                        apoint = ap;
                        let v = GpVec::from_pnts(p, &apoint);
                        out_par = v.magnitude();
                        out_lin = line_through(p, &v);
                        return Some(Segment {
                            lin: out_lin,
                            par: out_par,
                            flag: 0,
                        });
                    }
                    let (au, av, ap) = ext.point(i_near);
                    u_ = au;
                    v_ = av;
                    apoint = ap;
                }

                // Do-while over the grid samples (`cxx:651-693`).
                loop {
                    index_point += 1;
                    let found = self.point_in_the_face(
                        &face,
                        &mut apoint,
                        &mut u_,
                        &mut v_,
                        self.param_on_edge,
                        &mut index_point,
                        surf.as_ref(),
                        (u1, u2, v1, v2),
                        &mut vec_d1u,
                        &mut vec_d1v,
                    );
                    if found {
                        nb_points_ok += 1;
                        let v = GpVec::from_pnts(p, &apoint);
                        let par = v.magnitude();
                        if par > RESOLUTION
                            && vec_d1u.magnitude() > RESOLUTION
                            && vec_d1v.magnitude() > RESOLUTION
                        {
                            let norm = vec_d1u.clone().cross(&vec_d1v);
                            let mut tt = norm.magnitude();
                            if tt > RESOLUTION {
                                tt = norm.dot(&v).abs() / (tt * par);
                                if tt > max_scal {
                                    max_scal = tt;
                                    out_lin = line_through(p, &v);
                                    out_par = par;
                                    ptfound = true;
                                    if max_scal > 0.2 {
                                        self.param_on_edge = svmyparam;
                                        return Some(Segment {
                                            lin: out_lin,
                                            par: out_par,
                                            flag: 0,
                                        });
                                    }
                                }
                            }
                        }
                    }
                    if !(index_point < 200 && nb_points_ok < 16) {
                        break;
                    }
                }

                self.param_on_edge = svmyparam;
                if max_scal > 0.2 {
                    return Some(Segment {
                        lin: out_lin,
                        par: out_par,
                        flag: 0,
                    });
                }

                index_point = 0;
                let encore_une_face = k + 1 < n_faces;
                if !ptfound && !encore_une_face && self.param_on_edge < 0.0001 {
                    // Solid reduced to a face (`cxx:710-720`).
                    let pbidon = GpPnt::new(p.x() + 1.0, p.y(), p.z());
                    let v = GpVec::from_pnts(p, &pbidon);
                    out_par = 1.0;
                    out_lin = line_through(p, &v);
                    return Some(Segment {
                        lin: out_lin,
                        par: out_par,
                        flag: 0,
                    });
                }
            }

            // `NbFacesInSolid == 0` (`cxx:724-732`).
            if nb_faces_in_solid == 0 {
                return Some(Segment {
                    lin: out_lin,
                    par: 0.0,
                    flag: 0,
                });
            }
            if ptfound {
                return Some(Segment {
                    lin: out_lin,
                    par: out_par,
                    flag: 0,
                });
            }

            // Retry ladder of `myParamOnEdge` (`cxx:738-784`).
            self.first_face = 0;
            let ladder = [
                (0.512345, 0.4),
                (0.4, 0.6),
                (0.6, 0.3),
                (0.3, 0.7),
                (0.7, 0.2),
                (0.2, 0.8),
                (0.8, 0.1),
                (0.1, 0.9),
            ];
            let current = self.param_on_edge;
            match ladder.iter().find(|(from, _)| *from == current) {
                Some((_, to)) => self.param_on_edge = *to,
                None => {
                    self.param_on_edge *= 0.5;
                    if self.param_on_edge < 0.0001 {
                        let pbidon = GpPnt::new(p.x() + 1.0, p.y(), p.z());
                        let v = GpVec::from_pnts(p, &pbidon);
                        out_par = 1.0;
                        out_lin = line_through(p, &v);
                        return Some(Segment {
                            lin: out_lin,
                            par: out_par,
                            flag: 0,
                        });
                    }
                }
            }
            a_test_invert = true;
        }
    }
}
