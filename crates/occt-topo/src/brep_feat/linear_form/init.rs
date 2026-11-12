//! `BRepFeat_MakeLinearForm::Init` (`BRepFeat_MakeLinearForm.cxx:84-857`).
//!
//! Six phases, in the reference order:
//!
//! 1. reset the base state and the local maps (`cxx:97-120`);
//! 2. pick up the arguments, the max vertex tolerance and the two directions;
//!    when `Direc1` is null the rib is centred and the sketch wire, the two
//!    vectors and the plane are shifted by the half translation (`cxx:122-197`);
//! 3. the bounding box, the plane face, the `BndFace` the profile closes
//!    against, and the two extreme support faces (`cxx:199-277`);
//! 4. the sliding feasibility tests (`cxx:278-353`);
//! 5. the profile: `SlidingProfile` + `Propagate`, or
//!    `NoSlidingProfile` + `Propagate` (`cxx:356-846`);
//! 6. the initial `myMap` of the basis faces (`cxx:848-856`).

use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::gp::{GpPln, GpPnt, GpTrsf, GpVec};
use occt_core::precision;
use occt_geom::curve::Curve;
use occt_geom::trimmed::GeomTrimmedCurveBasis;
use occt_geom::GeomPlane;

use crate::abs::{Orientation, ShapeType};
use crate::brep_extrema_ext_cf::ExtCF;
use crate::brep_feat::algo_api;
use crate::brep_feat::rib_slot::{
    check_point, height_max, int_par, NoSlidingProfileArgs, SlidingProfileArgs,
};
use crate::brep_feat::status_error::StatusError;
use crate::brep_lib_make_edge;
use crate::brep_lib_make_face;
use crate::brep_lib_make_wire::MakeWire;
use crate::brep_tool::BRepTool;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::builder::TopoBuilder;
use crate::locope::shape_maps::ShapeValueMap;
use crate::primitives::BRepPrimBox;
use crate::shape::{Edge, Face, Shell, TopoShape, Vertex, Wire};
use crate::topo_tools::{explorer_shapes, is_same};

use super::MakeLinearForm;

/// `BRepLib_MakeVertex(p)` as a vertex handle.
fn make_vertex(p: &GpPnt) -> Vertex {
    TopoBuilder::new().make_vertex(*p, precision::CONFUSION)
}

/// `BRepLib_MakeEdge(C, V1, V2)` (`BRepLib_MakeEdge.cxx:542-572`).
fn edge_from_vertices(curve: &Arc<dyn Curve>, v1: &Vertex, v2: &Vertex) -> Option<Edge> {
    let build = brep_lib_make_edge::init_3d_projecting(curve, Some(v1), Some(v2))?;
    let b = TopoBuilder::new();
    let mut e = b.make_edge(curve.clone(), build.first, build.last);
    brep_lib_make_edge::add_endpoint_vertices(&b, &mut e, &build);
    Some(e)
}

/// `BRepLib_MakeEdge(C, u1, u2)`.
fn edge_from_params(curve: &Arc<dyn Curve>, first: f64, last: f64) -> Option<Edge> {
    let build = brep_lib_make_edge::init_3d(curve, first, last)?;
    let b = TopoBuilder::new();
    let mut e = b.make_edge(curve.clone(), build.first, build.last);
    brep_lib_make_edge::add_endpoint_vertices(&b, &mut e, &build);
    Some(e)
}

/// `BRepLib_MakeEdge(P1, P2)` (`BRepLib_MakeEdge.cxx:202-215`).
fn edge_from_points(p1: &GpPnt, p2: &GpPnt) -> Option<Edge> {
    let (curve, build) = brep_lib_make_edge::init_points(p1, p2)?;
    let b = TopoBuilder::new();
    let mut e = b.make_edge(curve, build.first, build.last);
    brep_lib_make_edge::add_endpoint_vertices(&b, &mut e, &build);
    Some(e)
}

/// `p + v` on `gp_Pnt` coordinates.
fn offset_pnt(p: &GpPnt, v: &GpVec) -> GpPnt {
    GpPnt::new(p.x() + v.x(), p.y() + v.y(), p.z() + v.z())
}

/// `BRepExtrema_ExtCF` on the segment `p + v` against `f`: exactly one
/// extremum, within the face tolerance (`cxx:305-353`).
fn sweep_hits(p: &GpPnt, v: &GpVec, f: &Face) -> bool {
    let p1 = offset_pnt(p, v);
    let Some(e) = edge_from_points(p, &p1) else {
        return false;
    };
    let ext = ExtCF::new(&e, f);
    let tol = BRepTool::face_tolerance(f);
    ext.nb_ext() == 1 && ext.square_distance(1) <= tol * tol
}

/// `new Geom_TrimmedCurve(ccc, f, l)` with `cc->Reverse()` when `edg` is
/// reversed (`cxx:438-442`, `cxx:578-585`).
///
/// Uses [`GeomTrimmedCurveBasis`] (basis parameters, like OCCT) and a real
/// `Reverse` — `GeomTrimmedCurve::new(l, f)` min/max-collapses and cannot
/// stand in for `cc->Reverse()`.
fn trimmed_of(edg: &Edge, ccc: &Arc<dyn Curve>, f: f64, l: f64) -> Arc<dyn Curve> {
    let mut cc = GeomTrimmedCurveBasis::new(ccc.clone(), f, l);
    if edg.0.orientation() == Orientation::Reversed {
        cc.reverse();
    }
    Arc::new(cc)
}

/// `cc->Reverse()` on a trimmed curve (`Geom_TrimmedCurve.cxx:78-83`).
fn reverse_curve(cc: &Arc<dyn Curve>) -> Arc<dyn Curve> {
    let mut c = cc.clone_dyn();
    c.reverse();
    Arc::from(c)
}

impl MakeLinearForm {
    /// `Init(Sbase, W, Plane, Direc, Direc1, Mode, Modify)` (`cxx:84-857`).
    ///
    /// `Err` is one of the reference's raised failures: a null `Plane`, a
    /// missing extreme face / edge, or a `BRepAlgoAPI_Common` that did not
    /// build.
    pub fn init(
        &mut self,
        sbase: &TopoShape,
        w: &Wire,
        plane: &GpPln,
        direction: &GpVec,
        direction1: &GpVec,
        mode: i32,
        modify: bool,
    ) -> Result<(), String> {
        // `cxx:97-120`.
        let revol_rib = false;
        self.base.done();
        self.base.my_generated.clear();

        let mut sliding = modify;

        self.base.my_lf_map.clear();
        self.base.my_shape = None;
        self.base.my_map.clear();
        self.base.my_f_shape = None;
        self.base.my_l_shape = None;
        self.base.my_sbase = Some(sbase.clone());
        self.base.my_skface = None;
        self.base.my_pbase = None;
        self.base.my_gshape = None;
        self.base.my_suntil = None;
        self.my_list_of_edges.clear();
        self.my_slface.clear();

        // `cxx:119-125`.
        self.base.my_wire = Some(Wire(w.0.oriented(Orientation::Forward)));
        self.my_dir = *direction;
        self.my_dir1 = *direction1;
        self.my_pln = Some(plane.clone());
        self.base.my_fuse = mode != 0;

        // `cxx:139-158`: the max vertex tolerance of the sketch and the basis.
        self.my_tol = precision::CONFUSION;
        if let Some(wire) = &self.base.my_wire {
            for v in explorer_shapes(&wire.0, ShapeType::Vertex) {
                let tol = BRepTool::vertex_tolerance(&Vertex(v));
                if tol > self.my_tol {
                    self.my_tol = tol;
                }
            }
        }
        for v in explorer_shapes(sbase, ShapeType::Vertex) {
            let tol = BRepTool::vertex_tolerance(&Vertex(v));
            if tol > self.my_tol {
                self.my_tol = tol;
            }
        }

        // `cxx:162-197`: the directions.
        let nulldir = GpVec::zero();
        if !self.my_dir1.is_equal_tol(&nulldir, self.my_tol, self.my_tol) {
            let ang = self.my_dir1.angle(&self.my_dir);
            if ang != PI {
                self.base.my_status_error = StatusError::BadDirect;
                self.base.not_done();
                return Ok(());
            }
        } else {
            // `cxx:181-196`: the rib is centred in the middle of the
            // translation.
            let dir_translation = self.my_dir.added(&self.my_dir1).multiplied_scalar(0.5);
            let mut t = GpTrsf::identity();
            t.set_translation_vec(&dir_translation);
            let wire = self
                .base
                .my_wire
                .clone()
                .ok_or_else(|| "Standard_NullObject: myWire".to_string())?;
            self.base.my_wire = Some(Wire(crate::transform::move_shape(&wire.0, &t)));
            self.my_dir = self.my_dir.subtracted(&dir_translation);
            self.my_dir1 = self.my_dir1.subtracted(&dir_translation);
            if let Some(pln) = &mut self.my_pln {
                pln.transform(&t);
            }
        }

        let pln = self
            .my_pln
            .clone()
            .ok_or_else(|| "Standard_NullObject: myPln".to_string())?;

        // `cxx:199-208`: `HeightMax` and its box.
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

        // `cxx:210-222`: the plane face, its common with the box and the face
        // closed on the resulting wire.
        let plane_face = brep_lib_make_face::make_face_uv(
            Arc::new(GeomPlane::new(pln.clone())),
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
            .ok_or_else(|| "Standard_DomainError: no wire in the plane section".to_string())?;
        let bnd_face = brep_lib_make_face::make_face_from_plane_and_wire(&pln, &Wire(www), true);

        // `cxx:224-260`: the support faces.
        let (data, ef) = self.base.extreme_faces(revol_rib, self.my_bnd, &pln, None, None)?;
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

        // `cxx:262-277`.
        let check_pnt = check_point(&first_edge, bnd / 10.0, &pln);
        let mut concavite = 3;
        self.base.my_first_pnt = BRepTool::vertex_point(&first_vertex);
        self.base.my_last_pnt = BRepTool::vertex_point(&last_vertex);

        // `SliList` (`cxx:276-278`).
        let mut sli_list: Vec<TopoShape> = vec![first_face.0.clone()];

        // `cxx:280-297`: the sliding is possible only on a plane or a cylinder.
        if sliding {
            sliding = false;
            if let Some(s) = BRepTool::face_surface(&first_face) {
                let s = s.rectangular_trimmed_basis().unwrap_or(s);
                if s.gp_pln().is_some() || s.gp_cylinder().is_some() {
                    sliding = true;
                }
            }
        }

        // `cxx:299-353`: the two support faces must be reached by the sweep.
        let first_pnt = self.base.my_first_pnt;
        let last_pnt = self.base.my_last_pnt;
        if sliding {
            sliding = sweep_hits(&first_pnt, &self.my_dir, &first_face)
                && sweep_hits(&last_pnt, &self.my_dir, &last_face);
        }
        if !self
            .my_dir1
            .is_equal_tol(&nulldir, precision::CONFUSION, precision::CONFUSION)
        {
            if sliding {
                if sweep_hits(&first_pnt, &self.my_dir1, &first_face) {
                    sliding = sweep_hits(&last_pnt, &self.my_dir1, &last_face);
                } else {
                    sliding = false;
                }
            }
        }

        // `cxx:356-411`: the sliding profile and its propagation.
        if sliding {
            let args = SlidingProfileArgs {
                revol_rib,
                tol: self.my_tol,
                pln: &pln,
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
            sliding =
                self.propagate(&mut sli_list, &prof, &first_pnt, &last_pnt, &mut falseside)?;
            if !falseside {
                self.base.my_status_error = StatusError::FalseSide;
                self.base.not_done();
                return Ok(());
            }
        }

        // `cxx:414-751`: the profile wire of the sliding case.
        let mut w = MakeWire::new();
        let mut the_previous_edge: Option<Edge> = None;
        let mut the_fv: Option<Vertex> = None;
        let mut counter = 1i32;

        if sliding && !self.my_list_of_edges.is_empty() {
            let Some(wire) = self.base.my_wire.clone() else {
                return Ok(());
            };
            let mut ex1 = WireExplorer::of_wire(&wire);

            // `cxx:426-487`: the first edge of the sketch.
            while ex1.more() {
                let current = ex1.current().cloned();
                let Some(e) = current else {
                    break;
                };
                if !self.base.my_lf_map.is_bound(&e.0) {
                    self.base.my_lf_map.bind(&e.0, Vec::new());
                }
                if is_same(&e.0, &first_edge.0) {
                    let (f, l) = BRepTool::edge_parameters(&e);
                    let ccc = BRepTool::edge_curve_world(&e)
                        .ok_or_else(|| "Standard_NullObject: BRep_Tool::Curve".to_string())?;
                    // `cxx:437-456`: on a one-edge sketch the profile runs from
                    // `myFirstPnt` to `myLastPnt` and is reversed when the
                    // parameters run the other way.
                    let (pt, reversed) = if !is_same(&first_edge.0, &last_edge.0) {
                        (
                            crate::shhealing::last_vertex(&e)
                                .map(|v| BRepTool::vertex_point(&v))
                                .ok_or_else(|| "Standard_NullObject: LastVertex".to_string())?,
                            false,
                        )
                    } else {
                        let fpar = int_par(Some(ccc.as_ref()), &first_pnt);
                        let lpar = int_par(Some(ccc.as_ref()), &last_pnt);
                        (last_pnt, fpar > lpar)
                    };
                    // `cxx:437-442`: `new Geom_TrimmedCurve(ccc,f,l)` then
                    // `cc->Reverse()` when the one-edge parameters run the
                    // other way — not `new(..., l, f)` (min/max collapses).
                    let mut trimmed = GeomTrimmedCurveBasis::new(ccc.clone(), f, l);
                    if reversed {
                        trimmed.reverse();
                    }
                    let cc: Arc<dyn Curve> = Arc::new(trimmed);
                    // `cxx:457-469`.
                    let v2 = make_vertex(&pt);
                    let built = match &the_previous_edge {
                        None => {
                            let v1 = make_vertex(&first_pnt);
                            edge_from_vertices(&cc, &v1, &v2)
                        }
                        Some(prev) => {
                            let v1 = crate::shhealing::last_vertex(prev)
                                .ok_or_else(|| "Standard_NullObject: LastVertex".to_string())?;
                            edge_from_vertices(&cc, &v1, &v2)
                        }
                    }
                    .ok_or_else(|| "BRepLib_LineThroughIdenticPoints".to_string())?;
                    let ee1 = Edge(built.0.oriented(e.0.orientation()));
                    if counter == 1 {
                        the_fv = crate::shhealing::first_vertex(&ee1);
                    }
                    if let Some(list) = self.base.my_lf_map.change(&e.0) {
                        list.push(ee1.0.clone());
                    }
                    let _ = w.add(&ee1);
                    the_previous_edge = Some(ee1);
                    counter += 1;
                    ex1.next();
                    break;
                }
                ex1.next();
            }

            // `cxx:489-560`: the edges up to the last one.
            if !is_same(&first_edge.0, &last_edge.0) {
                while ex1.more() {
                    let current = ex1.current().cloned();
                    let Some(e) = current else {
                        break;
                    };
                    if !self.base.my_lf_map.is_bound(&e.0) {
                        self.base.my_lf_map.bind(&e.0, Vec::new());
                    }
                    // `theList.Append(E)` (`cxx:497`) is written by the
                    // reference and never read; the port drops it.
                    let (f, l) = BRepTool::edge_parameters(&e);
                    if !is_same(&e.0, &last_edge.0) {
                        // `cxx:501-528`.
                        let ccc = BRepTool::edge_curve_world(&e).ok_or_else(|| {
                            "Standard_NullObject: BRep_Tool::Curve".to_string()
                        })?;
                        let (v1, v2) = match &the_previous_edge {
                            Some(prev) => (
                                crate::shhealing::last_vertex(prev),
                                crate::shhealing::last_vertex(&e),
                            ),
                            None => (
                                crate::shhealing::first_vertex(&e),
                                crate::shhealing::last_vertex(&e),
                            ),
                        };
                        let v1 = v1.ok_or_else(|| "Standard_NullObject: vertex".to_string())?;
                        let v2 = v2.ok_or_else(|| "Standard_NullObject: vertex".to_string())?;
                        let built = edge_from_vertices(&ccc, &v1, &v2).ok_or_else(|| {
                            "BRepLib_LineThroughIdenticPoints".to_string()
                        })?;
                        let e11 = Edge(built.0.oriented(e.0.orientation()));
                        if counter == 1 {
                            the_fv = crate::shhealing::first_vertex(&e11);
                        }
                        if let Some(list) = self.base.my_lf_map.change(&e.0) {
                            list.push(e11.0.clone());
                        }
                        let _ = w.add(&e11);
                        the_previous_edge = Some(e11);
                        counter += 1;
                        ex1.next();
                        continue;
                    }

                    // `cxx:530-559`: the last edge of the sketch.
                    let ccc = BRepTool::edge_curve_world(&e)
                        .ok_or_else(|| "Standard_NullObject: BRep_Tool::Curve".to_string())?;
                    let cc: Arc<dyn Curve> =
                        Arc::new(GeomTrimmedCurveBasis::new(ccc.clone(), f, l));
                    let pf = crate::shhealing::first_vertex(&e)
                        .map(|v| BRepTool::vertex_point(&v))
                        .ok_or_else(|| "Standard_NullObject: FirstVertex".to_string())?;
                    let built = match &the_previous_edge {
                        None => edge_from_points(&pf, &last_pnt),
                        Some(prev) => {
                            let v1 = crate::shhealing::last_vertex(prev)
                                .ok_or_else(|| "Standard_NullObject: LastVertex".to_string())?;
                            let v2 = make_vertex(&last_pnt);
                            edge_from_vertices(&cc, &v1, &v2)
                        }
                    }
                    .ok_or_else(|| "BRepLib_LineThroughIdenticPoints".to_string())?;
                    let ee = Edge(built.0.oriented(e.0.orientation()));
                    let _ = w.add(&ee);
                    if let Some(list) = self.base.my_lf_map.change(&e.0) {
                        list.push(ee.0.clone());
                    }
                    if counter == 1 {
                        the_fv = crate::shhealing::first_vertex(&ee);
                    }
                    the_previous_edge = Some(ee);
                    counter += 1;
                    break;
                }
            }

            // `cxx:562-705`: walk `myListOfEdges` from the last point to the
            // first one, restarting after each removal.
            let mut first_ok = false;
            let mut last_ok = false;
            let mut the_last_pnt = last_pnt;
            let mut sens = 0i32;
            let mut the_edge: Option<Edge> = None;
            let mut the_l_edge: Option<Edge> = None;
            let mut counter1 = counter;
            let mut new_list_of_edges: Vec<TopoShape> = Vec::new();
            let mut it_idx = 0usize;

            while !first_ok {
                let current = self.my_list_of_edges.get(it_idx).cloned();
                let Some(edg) = current else {
                    sliding = false;
                    break;
                };
                let edg = Edge::wrap(edg).ok_or_else(|| "Standard_DomainError".to_string())?;
                let (f, l) = BRepTool::edge_parameters(&edg);
                let ccc = BRepTool::edge_curve_world(&edg)
                    .ok_or_else(|| "Standard_NullObject: BRep_Tool::Curve".to_string())?;

                // `cxx:578-588`.
                let cc = trimmed_of(&edg, &ccc, f, l);
                let fp = cc.value(cc.first_parameter());
                let lp = cc.value(cc.last_parameter());
                let mut dist = fp.distance(&the_last_pnt);
                let mut cc = cc;
                if dist <= self.my_tol {
                    sens = 1;
                    last_ok = true;
                } else {
                    dist = lp.distance(&the_last_pnt);
                    if dist <= self.my_tol {
                        sens = 2;
                        last_ok = true;
                        // `cc->Reverse()` (`cxx:588`); not `new(..., l, f)`.
                        cc = reverse_curve(&cc);
                    }
                }

                // `cxx:605-616`.
                let mut first_flag = 0;
                if sens == 1 && lp.distance(&first_pnt) <= self.my_tol {
                    first_ok = true;
                    first_flag = 1;
                } else if sens == 2 && fp.distance(&first_pnt) <= self.my_tol {
                    first_ok = true;
                    first_flag = 2;
                }

                if last_ok {
                    // `cxx:618-668`.
                    let fpar = cc.first_parameter();
                    let lpar = cc.last_parameter();
                    let built = if !first_ok {
                        match &the_previous_edge {
                            None => edge_from_params(&cc, fpar, lpar),
                            Some(prev) => {
                                let v1 = crate::shhealing::last_vertex(prev).ok_or_else(|| {
                                    "Standard_NullObject: LastVertex".to_string()
                                })?;
                                // `BB.UpdateVertex(v1, dist)` (`cxx:630`).
                                v1.set_tolerance(dist);
                                let nv = make_vertex(&cc.value(lpar));
                                edge_from_vertices(&cc, &v1, &nv)
                            }
                        }
                    } else {
                        let fv = the_fv
                            .clone()
                            .ok_or_else(|| "Standard_NullObject: theFV".to_string())?;
                        match &the_previous_edge {
                            None => {
                                let v1 = make_vertex(&cc.value(fpar));
                                edge_from_vertices(&cc, &v1, &fv)
                            }
                            Some(prev) => {
                                let v1 = crate::shhealing::last_vertex(prev).ok_or_else(|| {
                                    "Standard_NullObject: LastVertex".to_string()
                                })?;
                                edge_from_vertices(&cc, &v1, &fv)
                            }
                        }
                    }
                    .ok_or_else(|| "BRepLib_LineThroughIdenticPoints".to_string())?;

                    let eeee = built;
                    the_previous_edge = Some(eeee.clone());
                    let _ = w.add(&eeee);
                    if counter == 1 {
                        the_fv = crate::shhealing::first_vertex(&eeee);
                    }
                    counter1 += 1;
                    new_list_of_edges.push(edg.0.clone());
                    the_edge = Some(eeee.clone());
                    // `theFEdge = edg` (`cxx:666`) is written and never read.
                    the_last_pnt = crate::shhealing::last_vertex(&eeee)
                        .map(|v| BRepTool::vertex_point(&v))
                        .ok_or_else(|| "Standard_NullObject: LastVertex".to_string())?;
                }

                // `cxx:670-678`.
                if first_flag == 1 {
                    the_l_edge = Some(edg.clone());
                } else if first_flag == 2 {
                    the_l_edge = the_edge.clone();
                }

                // `cxx:641-655`: after a consumed edge, re-initialize the
                // iterator on the remaining list (`it.Initialize(myListOfEdges)`).
                if last_ok {
                    self.my_list_of_edges.remove(it_idx);
                    it_idx = 0;
                    last_ok = false;
                } else if it_idx + 1 < self.my_list_of_edges.len() {
                    it_idx += 1;
                } else {
                    sliding = false;
                    break;
                }
                sens = 0;
            }
            let _ = the_l_edge;

            // `cxx:697-750`: the sliding faces of the edges added by the walk.
            let mut slid_map: ShapeValueMap<Vec<TopoShape>> = ShapeValueMap::new();
            if sliding && counter1 > counter {
                let wire = w.wire();
                let mut ii = 0i32;
                for e in explorer_shapes(&wire.0, ShapeType::Edge) {
                    ii += 1;
                    if ii >= counter && ii <= counter1 {
                        let mut jj = 0i32;
                        for e2 in &new_list_of_edges {
                            jj += 1;
                            if jj != ii - counter + 1 {
                                continue;
                            }
                            for (fac, ledg) in self.my_slface.iter() {
                                for e1 in ledg {
                                    if is_same(e1, e2) {
                                        if !slid_map.is_bound(fac) {
                                            slid_map.bind(fac, Vec::new());
                                        }
                                        if let Some(list) = slid_map.change(fac) {
                                            list.push(e.clone());
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            self.my_slface.clear();
            self.my_slface = slid_map;

            // `cxx:752-761`: the profile face on the plane.
            if sliding {
                let builder = TopoBuilder::new();
                let mut face = builder.make_face_with_tolerance(
                    Arc::new(GeomPlane::new(pln.clone())),
                    &[],
                    self.my_tol,
                );
                let wire = w.wire();
                wire.0.set_closed(BRepTool::is_closed(&wire.0));
                builder.add_wire(&mut face, &wire);
                self.base.my_skface = Some(face.clone());
                self.base.my_pbase = Some(face);
                self.base.my_suntil = None;
            }
        }

        // `cxx:763-846`: the profile of the non-sliding case.
        if !sliding {
            let args = NoSlidingProfileArgs {
                revol_rib,
                tol: self.my_tol,
                pln: &pln,
                bound: bnd,
                bnd_face: &bnd_face,
                check_pnt: &check_pnt,
                first_edge: &first_edge,
                last_edge: &last_edge,
                on_first_face,
                on_last_face,
            };
            let prof = match self.base.no_sliding_profile(&args, &mut concavite) {
                Ok(f) => f,
                Err(_) => {
                    self.base.my_status_error = StatusError::NoFaceProf;
                    self.base.not_done();
                    return Ok(());
                }
            };
            let mut falseside = true;
            self.propagate(&mut sli_list, &prof, &first_pnt, &last_pnt, &mut falseside)?;
            if !falseside {
                self.base.my_status_error = StatusError::FalseSide;
                self.base.not_done();
                return Ok(());
            }

            self.my_slface.clear();

            // `cxx:826-839`: the shell of the propagated faces.
            let builder = TopoBuilder::new();
            let mut comp = Shell::new();
            for s in &sli_list {
                builder.add(&mut comp.0, s);
            }
            comp.0.set_closed(BRepTool::is_closed(&comp.0));

            self.base.my_suntil = Some(comp.0.clone());
            self.base.my_skface = Some(prof.clone());
            self.base.my_pbase = Some(prof);
        }

        // `cxx:848`.
        self.base.my_sliding = sliding;

        // `cxx:850-856`.
        for f in explorer_shapes(sbase, ShapeType::Face) {
            self.base.my_map.bind(&f, vec![f.clone()]);
        }

        Ok(())
    }
}
