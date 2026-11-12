//! `BRepFeat_MakeRevolutionForm::Init` (`cxx:101-1121`).
//!
//! Phases, in the reference order:
//!
//! 1. the axis / sketch section (`cxx:119-141`) and the minimum radius `Rad`
//!    the profile reaches around it (`cxx:143-186`);
//! 2. the two revolution angles (`cxx:188-206`), the max vertex tolerance
//!    (`cxx:220-241`) and the argument pickup (`cxx:243-255`);
//! 3. the bounding box, the plane face, the `BndFace` the profile closes
//!    against (`cxx:257-281`) and the extreme support faces (`cxx:283-321`);
//! 4. the sliding feasibility tests, including the circular checks on the
//!    support cylinders (`cxx:323-514`);
//! 5. the profile: `SlidingProfile` + `Propagate`, then the profile wire, the
//!    `SlidMap` and the plane face (`cxx:520-920`); or, without sliding, the
//!    rebuilt `NewBndFace` and `NoSlidingProfile` + `Propagate` (`cxx:922-1110`).
//!    Both profile branches live in `init_profile.rs`.
//! 6. the initial `myMap` of the basis faces (`cxx:1112-1120`).

use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::gp::{GpAx2, GpCirc, GpDir, GpLin, GpPln, GpPnt, GpVec};
use occt_core::precision::{CONFUSION, SQUARE_CONFUSION};
use occt_geom::circle::GeomCircle;
use occt_geom::curve::Curve;
use occt_geom::geom_api::to_2d;
use occt_geom::line::GeomLine;
use occt_geom::trimmed::GeomTrimmedCurveBasis;
use occt_geom::GeomPlane;
use occt_geom2d::extrema2d::ExtremaExtCC2d;

use crate::abs::{Orientation, ShapeType};
use crate::brep_extrema_ext_cf::ExtCF;
use crate::brep_feat::algo_api;
use crate::brep_feat::rib_slot::{check_point, height_max, SlidingProfileArgs};
use crate::brep_feat::status_error::StatusError;
use crate::brep_lib_make_edge;
use crate::brep_lib_make_face;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::locope::cs_intersector::CSIntersector;
use crate::primitives::BRepPrimBox;
use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};
use crate::topo_tools::explorer_shapes;

use super::MakeRevolutionForm;

/// `RealLast()`: the largest finite value (`RealLast` is `DBL_MAX`).
const REAL_LAST: f64 = f64::MAX;

/// `BRepLib_MakeVertex(p)`.
pub(super) fn make_vertex(p: &GpPnt) -> Vertex {
    TopoBuilder::new().make_vertex(*p, CONFUSION)
}

/// `BRepLib_MakeEdge(C, V1, V2)` (`BRepLib_MakeEdge.cxx:542-572`).
pub(super) fn edge_from_vertices(curve: &Arc<dyn Curve>, v1: &Vertex, v2: &Vertex) -> Option<Edge> {
    let build = brep_lib_make_edge::init_3d_projecting(curve, Some(v1), Some(v2))?;
    let b = TopoBuilder::new();
    let mut e = b.make_edge(curve.clone(), build.first, build.last);
    brep_lib_make_edge::add_endpoint_vertices(&b, &mut e, &build);
    Some(e)
}

/// `BRepLib_MakeEdge(C, u1, u2)` (`BRepLib_MakeEdge.cxx:202-215`).
pub(super) fn edge_from_params(curve: &Arc<dyn Curve>, first: f64, last: f64) -> Option<Edge> {
    let build = brep_lib_make_edge::init_3d(curve, first, last)?;
    let b = TopoBuilder::new();
    let mut e = b.make_edge(curve.clone(), build.first, build.last);
    brep_lib_make_edge::add_endpoint_vertices(&b, &mut e, &build);
    Some(e)
}

/// `BRepLib_MakeEdge(V1, V2)` (`BRepLib_MakeEdge.cxx:189-198`): a line from the
/// first vertex point to the second, keeping both vertex handles.
fn edge_from_vertex_pair(v1: &Vertex, v2: &Vertex) -> Option<Edge> {
    let p1 = BRepTool::vertex_point(v1);
    let p2 = BRepTool::vertex_point(v2);
    let build = brep_lib_make_edge::init_line_vertices(&p1, &p2, v1, v2)?;
    let dir = GpDir::from_vec(&GpVec::from_pnts(&p1, &p2)).ok()?;
    let line: Arc<dyn Curve> = Arc::new(GeomLine::new(GpLin::from_pnt_dir(p1, dir)));
    let b = TopoBuilder::new();
    let mut e = b.make_edge(line, build.first, build.last);
    brep_lib_make_edge::add_endpoint_vertices(&b, &mut e, &build);
    Some(e)
}

/// `new Geom_TrimmedCurve(ccc, f, l)` with `cc->Reverse()` when `edg` is
/// reversed (`cxx:740-746`).
///
/// Uses [`GeomTrimmedCurveBasis`] and a real `Reverse` — constructing with
/// `(l, f)` cannot replace `cc->Reverse()` because min/max collapses.
pub(super) fn trimmed_of(edg: &Edge, ccc: &Arc<dyn Curve>, f: f64, l: f64) -> Arc<dyn Curve> {
    let mut cc = GeomTrimmedCurveBasis::new(ccc.clone(), f, l);
    if edg.0.orientation() == Orientation::Reversed {
        cc.reverse();
    }
    Arc::new(cc)
}

/// `BRep_Tool::Surface(F)` after the `Geom_RectangularTrimmedSurface` unwrap
/// (`cxx:345-351`): the four surface types `Init` accepts for sliding.
fn sliding_surface_ok(f: &Face) -> bool {
    match BRepTool::face_surface(f) {
        Some(s) => {
            let s = s.rectangular_trimmed_basis().unwrap_or(s);
            s.gp_pln().is_some()
                || s.gp_cylinder().is_some()
                || s.gp_cone().is_some()
                || s.gp_torus().is_some()
        }
        // A face with no surface cannot answer the reference's
        // `s->DynamicType()`; the reference would dereference a null handle.
        None => false,
    }
}

/// `BRepLib_MakeEdge(V1, V2)`, reporting the reference's
/// `BRepLib_LineThroughIdenticPoints` (`BRepLib_MakeEdge.cxx:189-193`).
pub(super) fn pair_edge(a: &Option<Vertex>, b: &Option<Vertex>) -> Result<Edge, String> {
    let a = a
        .as_ref()
        .ok_or_else(|| "Standard_NullObject: vertex".to_string())?;
    let b = b
        .as_ref()
        .ok_or_else(|| "Standard_NullObject: vertex".to_string())?;
    edge_from_vertex_pair(a, b).ok_or_else(|| "BRepLib_LineThroughIdenticPoints".to_string())
}

impl MakeRevolutionForm {
    /// `Init(Sbase, W, Plane, Axis, H1, H2, Mode, Modify)` (`cxx:101-1121`).
    ///
    /// `Err` is one of the reference's raised failures: a `GeomAPI::To2d` that
    /// answered a null handle, a `BRepLib_MakeEdge` on identic points, an
    /// invalid rebuilt bounding face, or a boolean that did not build.
    #[allow(clippy::too_many_arguments)]
    pub fn init(
        &mut self,
        sbase: &TopoShape,
        w: &Wire,
        plane: &GpPln,
        axis: &occt_core::gp::GpAx1,
        h1: f64,
        h2: f64,
        mode: i32,
        modify: bool,
    ) -> Result<(), String> {
        // `cxx:116-120`.
        let revol_rib = true;
        self.base.done();
        let mut sliding = modify;

        // `cxx:122-141`: the sketch section by the axis, which bounds the
        // 2-D extrema of the profile edges below.
        self.my_axe = *axis;
        let line: Arc<dyn Curve> = Arc::new(GeomLine::new(GpLin::new(*axis)));
        let mut asi = CSIntersector::with_shape(sbase);
        let _ = asi.perform_curves(&[Some(line.clone())]);
        let (line_first, line_last) = {
            let n = asi.nb_points(1).unwrap_or(0);
            if asi.is_done() && n >= 2 {
                (asi.point(1, 1)?.parameter(), asi.point(1, n)?.parameter())
            } else {
                (-REAL_LAST, REAL_LAST)
            }
        };

        // `cxx:143`.
        let ln2d = to_2d(line.as_ref(), plane)
            .ok_or_else(|| "Standard_NullObject: GeomAPI::To2d".to_string())?;

        // `cxx:145-186`: the smallest distance between the profile and the
        // axis, which the revolution angles are derived from.
        let mut rad = REAL_LAST;
        for e in explorer_shapes(&w.0, ShapeType::Edge) {
            let e = Edge::wrap(e).ok_or_else(|| "Standard_DomainError: wire edge".to_string())?;
            let (f, l) = BRepTool::edge_parameters(&e);
            let c = BRepTool::edge_curve_world(&e)
                .ok_or_else(|| "Standard_NullObject: BRep_Tool::Curve".to_string())?;
            let c2d = to_2d(c.as_ref(), plane)
                .ok_or_else(|| "Standard_NullObject: GeomAPI::To2d".to_string())?;
            let extr =
                ExtremaExtCC2d::new_range(ln2d.as_ref(), c2d.as_ref(), line_first, line_last, f, l);
            let mut dist = REAL_LAST;
            if extr.nb_ext() >= 1 {
                // `Geom2dAPI_ExtremaCurveCurve::LowerDistance()` is
                // `sqrt(myExtCC.SquareDistance(myIndex))` with `myIndex` the
                // minimum (`Geom2dAPI_ExtremaCurveCurve.cxx`).
                let mut best = REAL_LAST;
                for n in 1..=extr.nb_ext() {
                    best = best.min(extr.square_distance(n));
                }
                dist = best.sqrt();
            }
            let p1 = c.value(f);
            let p2 = c.value(l);
            // `GeomAPI_ProjectPointOnCurve(p, Line)`: `Distance(1)` on a line is
            // the single solution, which is also the port's minimum.
            let proj1 = occt_geom::geom_api::project_point_on_curve(line.as_ref(), &p1, CONFUSION);
            let proj2 = occt_geom::geom_api::project_point_on_curve(line.as_ref(), &p2, CONFUSION);
            let (Some(proj1), Some(proj2)) = (proj1, proj2) else {
                self.base.my_status_error = StatusError::NoProjPt;
                self.base.not_done();
                return Ok(());
            };
            let par = proj1.distance.min(proj2.distance);
            if par < dist {
                dist = par;
            }
            if dist < rad && dist > 0.0 {
                rad = dist;
            }
        }

        // `cxx:188-206`.
        let height = h1.min(h2);
        if rad <= height {
            rad = height + 0.01 * height;
        }
        self.my_angle1 = (h1 / rad).asin() + PI / 10.0;
        self.my_angle2 = (h2 / rad).asin() + PI / 10.0;
        if (self.my_angle1 - PI / 2.0) > CONFUSION {
            self.my_angle1 = PI / 2.0;
        }
        if (self.my_angle2 - PI / 2.0) > CONFUSION {
            self.my_angle2 = PI / 2.0;
        }

        // `cxx:208-209`.
        self.base.my_skface = None;
        self.base.my_pbase = None;
        self.base.my_fuse = mode != 0;

        // `cxx:220-241`.
        self.my_tol = CONFUSION;
        for v in explorer_shapes(&w.0, ShapeType::Vertex) {
            let tol = BRepTool::vertex_tolerance(&Vertex(v));
            if tol > self.my_tol {
                self.my_tol = tol;
            }
        }
        for v in explorer_shapes(sbase, ShapeType::Vertex) {
            let tol = BRepTool::vertex_tolerance(&Vertex(v));
            if tol > self.my_tol {
                self.my_tol = tol;
            }
        }

        // `cxx:243-255`.
        self.base.my_wire = Some(Wire(w.0.oriented(Orientation::Forward)));
        self.my_pln = Some(plane.clone());
        self.my_height1 = h1;
        self.my_height2 = h2;
        self.base.my_sbase = Some(sbase.clone());
        self.my_slface.clear();
        self.base.my_shape = None;
        self.base.my_map.clear();
        self.base.my_f_shape = None;
        self.base.my_l_shape = None;

        // `cxx:257-281`: the bounding box, the plane face, its common with the
        // box, and the face `BndFace` closed on the resulting wire.
        let (first_corner, last_corner, bnd) = height_max(sbase, None);
        self.my_bnd = bnd;
        // `BRepPrimAPI_MakeBox` on a flat corner pair cannot be built by the
        // port's primitive, which asserts a positive extent; a void bounding
        // box already raises in `HeightMax` (`Bnd_Box::Get`).
        if !(last_corner.x() > first_corner.x()
            && last_corner.y() > first_corner.y()
            && last_corner.z() > first_corner.z())
        {
            self.base.my_status_error = StatusError::NoExtFace;
            self.base.not_done();
            return Ok(());
        }
        let bnd_box = BRepPrimBox::make_box_corner(&first_corner, &last_corner)
            .solid
            .0;
        let plane_face = brep_lib_make_face::make_face_uv(
            Arc::new(GeomPlane::new(plane.clone())),
            -6.0 * bnd,
            6.0 * bnd,
            -6.0 * bnd,
            6.0 * bnd,
            0.0,
        );
        let mut plane_s = algo_api::common(&bnd_box, &plane_face.0);
        plane_s.build()?;
        let plane_sect = plane_s
            .shape()
            .cloned()
            .ok_or_else(|| "StdFail_NotDone: BRepAlgoAPI_Common".to_string())?;
        let www = explorer_shapes(&plane_sect, ShapeType::Wire)
            .first()
            .cloned()
            .ok_or_else(|| "Standard_NoSuchObject: plane section wire".to_string())?;
        let bnd_face = brep_lib_make_face::make_face_from_plane_and_wire(plane, &Wire(www), true);

        // `cxx:283-321`.
        let (data, ef) = self.base.extreme_faces(revol_rib, self.my_bnd, plane, None, None)?;
        if !data {
            self.base.my_status_error = StatusError::NoExtFace;
            self.base.not_done();
            return Ok(());
        }
        let first_edge = ef
            .first_edge
            .clone()
            .ok_or_else(|| "Standard_NullObject: FirstEdge".to_string())?;
        let last_edge = ef
            .last_edge
            .clone()
            .ok_or_else(|| "Standard_NullObject: LastEdge".to_string())?;
        let first_face = ef
            .first_face
            .clone()
            .ok_or_else(|| "Standard_NullObject: FirstFace".to_string())?;
        let last_face = ef
            .last_face
            .clone()
            .ok_or_else(|| "Standard_NullObject: LastFace".to_string())?;
        let first_vertex = ef
            .first_vertex
            .clone()
            .ok_or_else(|| "Standard_NullObject: FirstVertex".to_string())?;
        let last_vertex = ef
            .last_vertex
            .clone()
            .ok_or_else(|| "Standard_NullObject: LastVertex".to_string())?;
        let on_first_face = ef.on_first_face;
        let on_last_face = ef.on_last_face;
        let pt_on_first_edge = ef.pt_on_first_edge;
        let pt_on_last_edge = ef.pt_on_last_edge;
        let on_first_edge = ef.on_first_edge.clone();
        let on_last_edge = ef.on_last_edge.clone();

        // `cxx:323-337`.
        let check_pnt = check_point(&first_edge, bnd / 10.0, plane);
        let mut concavite = 3;
        self.base.my_first_pnt = BRepTool::vertex_point(&first_vertex);
        self.base.my_last_pnt = BRepTool::vertex_point(&last_vertex);
        let mut sli_list: Vec<TopoShape> = vec![first_face.0.clone()];

        // `cxx:339-372`: only a plane, cylinder, cone or torus can slide.
        if sliding {
            sliding = sliding_surface_ok(&first_face);
        }
        if sliding {
            sliding = sliding_surface_ok(&last_face);
        }

        // `cxx:376-513`: the two support circles of the profile ends.
        if sliding {
            let first_pnt = self.base.my_first_pnt;
            let last_pnt = self.base.my_last_pnt;
            let proj = occt_geom::geom_api::project_point_on_curve(line.as_ref(), &first_pnt, CONFUSION)
                .ok_or_else(|| {
                    self.base.my_status_error = StatusError::NoProjPt;
                    String::new()
                });
            let proj = match proj {
                Ok(p) => p,
                Err(_) => {
                    self.base.not_done();
                    return Ok(());
                }
            };
            let first_rayon = proj.distance;
            let first_center = proj.point;

            let proj1 =
                occt_geom::geom_api::project_point_on_curve(line.as_ref(), &last_pnt, CONFUSION);
            let Some(proj1) = proj1 else {
                self.base.my_status_error = StatusError::NoProjPt;
                self.base.not_done();
                return Ok(());
            };
            let last_rayon = proj1.distance;
            let last_center = proj1.point;

            let ax2 = GpAx2::from_axis(first_center, self.my_axe.vdir);
            let ax2p = GpAx2::from_axis(last_center, self.my_axe.vdir);
            let the_fc = GpCirc::new(ax2, first_rayon);
            let the_lc = GpCirc::new(ax2p, last_rayon);

            let r_first1 = first_pnt.rotated(&self.my_axe, self.my_angle1);
            let r_last1 = last_pnt.rotated(&self.my_axe, self.my_angle1);
            let r_first2 = first_pnt.rotated(&self.my_axe, -self.my_angle2);
            let r_last2 = last_pnt.rotated(&self.my_axe, -self.my_angle2);

            let v1 = make_vertex(&r_first2);
            let v2 = make_vertex(&r_first1);
            let v3 = make_vertex(&r_last2);
            let v4 = make_vertex(&r_last1);
            let fc_curve: Arc<dyn Curve> = Arc::new(GeomCircle::new(the_fc.clone()));
            let lc_curve: Arc<dyn Curve> = Arc::new(GeomCircle::new(the_lc.clone()));
            let ee1 = edge_from_vertices(&fc_curve, &v1, &v2)
                .ok_or_else(|| "BRepLib_LineThroughIdenticPoints".to_string())?;
            let ee2 = edge_from_vertices(&lc_curve, &v3, &v4)
                .ok_or_else(|| "BRepLib_LineThroughIdenticPoints".to_string())?;

            // `cxx:432-447`.
            if sliding && !pt_on_first_edge {
                let ext1 = ExtCF::new(&ee1, &first_face);
                if ext1.nb_ext() < 1 || ext1.square_distance(1) > SQUARE_CONFUSION {
                    sliding = false;
                }
            }
            if sliding && !pt_on_last_edge {
                let ext2 = ExtCF::new(&ee2, &last_face);
                if ext2.nb_ext() < 1 || ext2.square_distance(1) > SQUARE_CONFUSION {
                    sliding = false;
                }
            }

            // `cxx:449-479`.
            if sliding && pt_on_first_edge {
                // A null `OnFirstEdge` would be dereferenced by the reference;
                // the port answers "not a circle", i.e. no sliding.
                match on_first_edge.as_ref().and_then(BRepTool::edge_curve_world) {
                    None => sliding = false,
                    Some(curve) => match curve.gp_circ() {
                        None => sliding = false,
                        Some(circ) => {
                            let circax = circ.axis();
                            if !circax.is_coaxial_tol(&self.my_axe, CONFUSION, CONFUSION) {
                                sliding = false;
                            } else if (circ.radius() - first_rayon).abs() >= CONFUSION {
                                sliding = false;
                            }
                        }
                    },
                }
            }

            // `cxx:481-511`.
            if sliding && pt_on_last_edge {
                match on_last_edge.as_ref().and_then(BRepTool::edge_curve_world) {
                    None => sliding = false,
                    Some(curve) => match curve.gp_circ() {
                        None => sliding = false,
                        Some(circ) => {
                            let circax = circ.axis();
                            if !circax.is_coaxial_tol(&self.my_axe, CONFUSION, CONFUSION) {
                                sliding = false;
                            } else if (circ.radius() - last_rayon).abs() >= CONFUSION {
                                sliding = false;
                            }
                        }
                    },
                }
            }
        }

        // `cxx:520-566`: the sliding profile and its propagation.
        if sliding {
            let args = SlidingProfileArgs {
                revol_rib,
                tol: self.my_tol,
                pln: plane,
                bnd_face: &bnd_face,
                check_pnt: &check_pnt,
                first_face: &first_face,
                last_face: &last_face,
                first_edge: &first_edge,
                last_edge: &last_edge,
            };
            let prof = match self.base.sliding_profile(&args, &mut concavite) {
                Ok(f) => f,
                Err(_) => {
                    self.base.my_status_error = StatusError::NoFaceProf;
                    self.base.not_done();
                    return Ok(());
                }
            };
            let mut falseside = true;
            let first_pnt = self.base.my_first_pnt;
            let last_pnt = self.base.my_last_pnt;
            sliding = self.propagate(
                &mut sli_list,
                &prof,
                &first_pnt,
                &last_pnt,
                &mut falseside,
            )?;
            if !falseside {
                self.base.my_status_error = StatusError::FalseSide;
                self.base.not_done();
                return Ok(());
            }
        }

        // `cxx:572-578`.
        // `cxx:586-920`: the profile wire of the sliding case.
        if !self.init_sliding_profile(plane, &first_edge, &last_edge, &mut sliding)? {
            return Ok(());
        }

        // `cxx:922-1110`: the profile of the non-sliding case.
        if !self.init_no_sliding_profile(
            revol_rib,
            plane,
            &ln2d,
            bnd,
            &check_pnt,
            &first_edge,
            &last_edge,
            on_first_face,
            on_last_face,
            sliding,
            &bnd_face,
            &mut sli_list,
            &mut concavite,
        )? {
            return Ok(());
        }

        // `cxx:1112`.
        self.my_sliding = sliding;

        // `cxx:1114-1120`.
        for f in explorer_shapes(sbase, ShapeType::Face) {
            self.base.my_map.bind(&f, vec![f.clone()]);
        }

        Ok(())
    }
}
