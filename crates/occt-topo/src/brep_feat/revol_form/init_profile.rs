//! The two profile branches of `BRepFeat_MakeRevolutionForm::Init`
//! (`BRepFeat_MakeRevolutionForm.cxx:520-1110`).
//!
//! `init.rs` keeps the preparation (axis and sketch section, revolution
//! angles, vertex tolerance, argument pickup, bounding face, extreme faces,
//! sliding feasibility) and calls these two methods in the reference order:
//! the sliding profile wire with its plane face (`cxx:586-920`) when sliding
//! survived, and the rebuilt bounding face with `NoSlidingProfile`
//! (`cxx:922-1108`) when it did not.
//!
//! Both answer `Ok(false)` on the paths where the reference leaves `Init`
//! early with a status error, so `init` stops before the `mySliding` write and
//! the `myMap` fill of `cxx:1112-1120`.

use std::sync::Arc;

use occt_core::elib::slib;
use occt_core::gp::{GpPln, GpPnt, GpPnt2d};
use occt_core::precision::CONFUSION;
use occt_geom::curve::Curve;
use occt_geom::geom_api::to_2d;
use occt_geom::trimmed::GeomTrimmedCurveBasis;
use occt_geom::GeomPlane;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::geom2d_api_inter_curve_curve::InterCurveCurve2d;

use crate::abs::ShapeType;
use crate::brep_check;
use crate::brep_feat::algo_api;
use crate::brep_feat::rib_slot::{int_par, NoSlidingProfileArgs};
use crate::brep_feat::status_error::StatusError;
use crate::brep_lib_make_edge;
use crate::brep_lib_make_face;
use crate::brep_lib_make_wire::MakeWire;
use crate::brep_tool::BRepTool;
use crate::brep_top_adaptor_fclass2d::BrepTopAdaptorFClass2d;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::builder::TopoBuilder;
use crate::fclass2d::FaceState;
use crate::locope::shape_maps::ShapeValueMap;
use crate::shape::{Edge, Face, Shell, TopoShape, Vertex, Wire};
use crate::topo_tools::{explorer_shapes, is_same};

use super::init::{edge_from_params, edge_from_vertices, make_vertex, pair_edge, trimmed_of};
use super::MakeRevolutionForm;

impl MakeRevolutionForm {
    /// `cxx:586-920`: the profile wire of the sliding case.
    ///
    /// Runs only when `sliding` is still set and `myListOfEdges` is not empty;
    /// answers `Ok(true)` otherwise, exactly like the reference's guard walks
    /// past the block.
    pub(super) fn init_sliding_profile(
        &mut self,
        plane: &GpPln,
        first_edge: &Edge,
        last_edge: &Edge,
        sliding: &mut bool,
    ) -> Result<bool, String> {
        if !*sliding || self.my_list_of_edges.is_empty() {
            return Ok(true);
        }

        let mut wprof = MakeWire::new();
        let mut the_previous_edge: Option<Edge> = None;
        let mut the_fv: Option<Vertex> = None;
        let mut counter = 1i32;

        let Some(wire) = self.base.my_wire.clone() else {
            // `cxx:588` dereferences `myWire`; the port stops `Init` here.
            return Ok(false);
        };
        let mut ex1 = WireExplorer::of_wire(&wire);

        // `cxx:588-647`: the first edge of the sketch.
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
                // `cxx:597-613`.
                let built = if !is_same(&first_edge.0, &last_edge.0) {
                    let pt = crate::shhealing::last_vertex(&e)
                        .map(|v| BRepTool::vertex_point(&v))
                        .ok_or_else(|| "Standard_NullObject: LastVertex".to_string())?;
                    let v2 = make_vertex(&pt);
                    match &the_previous_edge {
                        None => {
                            let v1 = make_vertex(&self.base.my_first_pnt);
                            edge_from_vertices(&ccc, &v1, &v2)
                        }
                        Some(prev) => {
                            let v1 = crate::shhealing::last_vertex(prev).ok_or_else(|| {
                                "Standard_NullObject: LastVertex".to_string()
                            })?;
                            edge_from_vertices(&ccc, &v1, &v2)
                        }
                    }
                } else {
                    let first_pnt = self.base.my_first_pnt;
                    let last_pnt = self.base.my_last_pnt;
                    let fpar = int_par(Some(ccc.as_ref()), &first_pnt);
                    let lpar = int_par(Some(ccc.as_ref()), &last_pnt);
                    // `cxx:437-442` family: Reverse when params run the other
                    // way — not `new(..., l, f)` (min/max collapses).
                    let mut trimmed = GeomTrimmedCurveBasis::new(ccc.clone(), f, l);
                    if fpar > lpar {
                        trimmed.reverse();
                    }
                    let cc: Arc<dyn Curve> = Arc::new(trimmed);
                    let v2 = make_vertex(&last_pnt);
                    match &the_previous_edge {
                        None => {
                            let v1 = make_vertex(&first_pnt);
                            edge_from_vertices(&cc, &v1, &v2)
                        }
                        Some(prev) => {
                            let v1 = crate::shhealing::last_vertex(prev).ok_or_else(|| {
                                "Standard_NullObject: LastVertex".to_string()
                            })?;
                            edge_from_vertices(&cc, &v1, &v2)
                        }
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
                let _ = wprof.add(&ee1);
                the_previous_edge = Some(ee1);
                counter += 1;
                ex1.next();
                break;
            }
            ex1.next();
        }

        // `cxx:648-722`: the edges up to the last one.
        if !is_same(&first_edge.0, &last_edge.0) {
            while ex1.more() {
                let current = ex1.current().cloned();
                let Some(e) = current else {
                    break;
                };
                if !self.base.my_lf_map.is_bound(&e.0) {
                    self.base.my_lf_map.bind(&e.0, Vec::new());
                }
                let (f, l) = BRepTool::edge_parameters(&e);
                if !is_same(&e.0, &last_edge.0) {
                    // `cxx:660-689`.
                    let ccc = BRepTool::edge_curve_world(&e)
                        .ok_or_else(|| "Standard_NullObject: BRep_Tool::Curve".to_string())?;
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
                    let e11 = edge_from_vertices(&ccc, &v1, &v2)
                        .ok_or_else(|| "BRepLib_LineThroughIdenticPoints".to_string())?;
                    let e11 = Edge(e11.0.oriented(e.0.orientation()));
                    the_previous_edge = Some(e11.clone());
                    if let Some(list) = self.base.my_lf_map.change(&e.0) {
                        list.push(e11.0.clone());
                    }
                    let _ = wprof.add(&e11);
                    if counter == 1 {
                        the_fv = crate::shhealing::first_vertex(&e11);
                    }
                    counter += 1;
                    ex1.next();
                    continue;
                }

                // `cxx:690-721`: the last edge of the sketch.
                let ccc = BRepTool::edge_curve_world(&e)
                    .ok_or_else(|| "Standard_NullObject: BRep_Tool::Curve".to_string())?;
                let pf = crate::shhealing::first_vertex(&e)
                    .map(|v| BRepTool::vertex_point(&v))
                    .ok_or_else(|| "Standard_NullObject: FirstVertex".to_string())?;
                let pl = self.base.my_last_pnt;
                let built = match &the_previous_edge {
                    None => {
                        // `BRepLib_MakeEdge e(cc, pf, pl)`.
                        let build = brep_lib_make_edge::init_curve_points(&ccc, &pf, &pl)
                            .ok_or_else(|| "BRepLib_PointProjectionFailed".to_string())?;
                        let b = TopoBuilder::new();
                        let mut e = b.make_edge(ccc.clone(), build.first, build.last);
                        brep_lib_make_edge::add_endpoint_vertices(&b, &mut e, &build);
                        Ok::<Edge, String>(e)
                    }
                    Some(prev) => {
                        let v1 = crate::shhealing::last_vertex(prev).ok_or_else(|| {
                            "Standard_NullObject: LastVertex".to_string()
                        })?;
                        let v2 = make_vertex(&pl);
                        let cc: Arc<dyn Curve> =
                            Arc::new(GeomTrimmedCurveBasis::new(ccc.clone(), f, l));
                        edge_from_vertices(&cc, &v1, &v2)
                            .ok_or_else(|| "BRepLib_LineThroughIdenticPoints".to_string())
                    }
                }?;
                let ee = Edge(built.0.oriented(e.0.orientation()));
                let _ = wprof.add(&ee);
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

        // `cxx:726-857`: walk `myListOfEdges` back to `myFirstPnt`. The
        // reference restarts the iterator after each consumed edge
        // (`it.Initialize(myListOfEdges)` at `cxx:840`) and never removes the
        // edge from the list, so the walk may revisit it.
        let mut first_ok_w = false;
        let mut last_ok_w = false;
        let mut the_last_pnt = self.base.my_last_pnt;
        let mut sens = 0i32;
        let mut counter1 = counter;
        let mut new_list_of_edges: Vec<TopoShape> = Vec::new();
        let mut it_idx = 0usize;
        let my_first_pnt = self.base.my_first_pnt;

        while !first_ok_w {
            let current = self.my_list_of_edges.get(it_idx).cloned();
            let Some(edg) = current else {
                *sliding = false;
                break;
            };
            let edg = Edge::wrap(edg).ok_or_else(|| "Standard_DomainError".to_string())?;
            let (f, l) = BRepTool::edge_parameters(&edg);
            let ccc = BRepTool::edge_curve_world(&edg)
                .ok_or_else(|| "Standard_NullObject: BRep_Tool::Curve".to_string())?;

            // `cxx:739-746`.
            let cc = trimmed_of(&edg, &ccc, f, l);
            let fp = cc.value(cc.first_parameter());
            let lp = cc.value(cc.last_parameter());
            let mut dist = fp.distance(&the_last_pnt);
            let mut cc = cc;
            if dist <= self.my_tol {
                sens = 1;
                last_ok_w = true;
            } else {
                dist = lp.distance(&the_last_pnt);
                if dist <= self.my_tol {
                    sens = 2;
                    last_ok_w = true;
                    // `cc->Reverse()` (`cxx:749`); not `new(..., l, f)`.
                    let mut c = cc.clone_dyn();
                    c.reverse();
                    cc = Arc::from(c);
                }
            }

            // `cxx:764-776`. The reference's `FirstFlag` is written and never
            // read, so the port keeps only `FirstOK`.
            if sens == 1 && lp.distance(&my_first_pnt) <= self.my_tol {
                first_ok_w = true;
            } else if sens == 2 && fp.distance(&my_first_pnt) <= self.my_tol {
                first_ok_w = true;
            }

            if last_ok_w {
                // `cxx:778-823`.
                let fpar = cc.first_parameter();
                let lpar = cc.last_parameter();
                let built = if !first_ok_w {
                    match &the_previous_edge {
                        None => edge_from_params(&cc, fpar, lpar),
                        Some(prev) => {
                            let v1 = crate::shhealing::last_vertex(prev).ok_or_else(|| {
                                "Standard_NullObject: LastVertex".to_string()
                            })?;
                            let v2 = make_vertex(&cc.value(lpar));
                            edge_from_vertices(&cc, &v1, &v2)
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
                let _ = wprof.add(&eeee);
                if counter == 1 {
                    the_fv = crate::shhealing::first_vertex(&eeee);
                }
                counter1 += 1;
                new_list_of_edges.push(edg.0.clone());
                the_last_pnt = crate::shhealing::last_vertex(&eeee)
                    .map(|v| BRepTool::vertex_point(&v))
                    .ok_or_else(|| "Standard_NullObject: LastVertex".to_string())?;
            }

            // `cxx:829-837`: `theLEdge = edg` / `theLEdge = theEdge` are
            // written and never read.
            if last_ok_w {
                // `it.Initialize(myListOfEdges)` (`cxx:840`).
                it_idx = 0;
                last_ok_w = false;
            } else if it_idx + 1 < self.my_list_of_edges.len() {
                it_idx += 1;
            } else {
                *sliding = false;
                break;
            }
            sens = 0;
        }

        // `cxx:851-907`: the sliding faces of the edges the walk added.
        let mut slid_map: ShapeValueMap<Vec<TopoShape>> = ShapeValueMap::new();
        if *sliding && counter1 > counter {
            let wire = wprof.wire();
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

        // `cxx:910-920`: the profile face on the plane.
        if *sliding {
            let builder = TopoBuilder::new();
            let mut face = builder.make_face_with_tolerance(
                Arc::new(GeomPlane::new(plane.clone())),
                &[],
                self.my_tol,
            );
            let wire = wprof.wire();
            wire.0.set_closed(BRepTool::is_closed(&wire.0));
            builder.add_wire(&mut face, &wire);
            self.base.my_skface = Some(face.clone());
            self.base.my_pbase = Some(face);
            self.base.my_suntil = None;
        }

        Ok(true)
    }

    /// `cxx:922-1110`: the profile of the non-sliding case.
    ///
    /// Answers `Ok(true)` when sliding was still set, i.e. when the reference's
    /// `if (!Sliding)` guard walks past the whole block.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn init_no_sliding_profile(
        &mut self,
        revol_rib: bool,
        plane: &GpPln,
        ln2d: &Arc<dyn Curve2d>,
        bnd: f64,
        check_pnt: &GpPnt,
        first_edge: &Edge,
        last_edge: &Edge,
        on_first_face: bool,
        on_last_face: bool,
        sliding: bool,
        bnd_face: &Face,
        sli_list: &mut Vec<TopoShape>,
        concavite: &mut i32,
    ) -> Result<bool, String> {
        if sliding {
            return Ok(true);
        }
        let mut bnd_face = bnd_face.clone();

        let www = explorer_shapes(&bnd_face.0, ShapeType::Wire)
            .first()
            .cloned()
            .ok_or_else(|| "Standard_NoSuchObject: BndFace wire".to_string())?;
        let mut explo = WireExplorer::of_wire(&Wire(www));
        let mut wiwiwi = MakeWire::new();
        let mut new_v1: Option<Vertex> = None;
        let mut last_v: Option<Vertex> = None;

        while explo.more() {
            let current = explo.current().cloned();
            let Some(e) = current else {
                break;
            };
            let e = e.clone();
            let v2 = crate::shhealing::last_vertex(&e);
            let (f, l) = BRepTool::edge_parameters(&e);
            let ln = BRepTool::edge_curve_world(&e)
                .ok_or_else(|| "Standard_NullObject: BRep_Tool::Curve".to_string())?;
            let l2d = to_2d(ln.as_ref(), plane)
                .ok_or_else(|| "Standard_NullObject: GeomAPI::To2d".to_string())?;
            let intcc = InterCurveCurve2d::with_curves(&l2d, ln2d, CONFUSION);
            let mut vv: Option<Vertex> = None;
            if intcc.nb_points() > 0 {
                let p = intcc
                    .point(1)
                    .ok_or_else(|| "Standard_RangeError: intersection point".to_string())?;
                let point = slib::plane_value(plane, p.x(), p.y());
                let par = int_par(Some(ln.as_ref()), &point);
                if f <= par && l >= par {
                    vv = Some(make_vertex(&point));
                }
            }

            // `cxx:980-1009`.
            if vv.is_none() && new_v1.is_none() {
                explo.next();
                continue;
            }
            if vv.is_some() && new_v1.is_none() {
                new_v1 = vv;
                last_v = v2;
                let e1 = pair_edge(&new_v1, &last_v)?;
                let _ = wiwiwi.add(&e1);
                explo.next();
                continue;
            }
            if vv.is_none() && new_v1.is_some() {
                // The reference builds `ee1 = BRepLib_MakeEdge(LastV, v2)`
                // and then adds the *original* edge `e` (`cxx:995-999`).
                let _ = pair_edge(&last_v, &v2)?;
                last_v = v2;
                let _ = wiwiwi.add(&e);
                explo.next();
                continue;
            }
            // `cxx:1002-1010`: both ends reached.
            let new_v2 = vv;
            let e1 = pair_edge(&last_v, &new_v2)?;
            last_v = new_v2.clone();
            let _ = wiwiwi.add(&e1);
            let e2 = pair_edge(&last_v, &new_v1)?;
            let _ = wiwiwi.add(&e2);
            break;
        }

        // `cxx:1011-1033`.
        let wiwiwi = wiwiwi.wire();
        wiwiwi.0.set_closed(BRepTool::is_closed(&wiwiwi.0));
        let mut new_bnd_face =
            brep_lib_make_face::make_face_from_plane_and_wire(plane, &wiwiwi, true);

        let cl = BrepTopAdaptorFClass2d::new(&new_bnd_face, CONFUSION);
        let (paru, parv) = slib::plane_parameters(&plane.position(), check_pnt);
        if cl.perform(GpPnt2d::new(paru, parv), true) == FaceState::Out {
            let mut cut = algo_api::cut(&bnd_face.0, &new_bnd_face.0);
            cut.build()?;
            let cut_shape = cut
                .shape()
                .cloned()
                .ok_or_else(|| "StdFail_NotDone: BRepAlgoAPI_Cut".to_string())?;
            if let Some(cur_wire) = explorer_shapes(&cut_shape, ShapeType::Wire).first() {
                new_bnd_face =
                    brep_lib_make_face::make_face_from_plane_and_wire(plane, &Wire(cur_wire.clone()), true);
            }
        }

        // `cxx:1030-1039`.
        if !brep_check::is_valid(&new_bnd_face.0) {
            self.base.my_status_error = StatusError::InvShape;
            self.base.not_done();
            return Ok(false);
        }
        bnd_face = new_bnd_face;

        // `cxx:1042-1075`.
        let args = NoSlidingProfileArgs {
            revol_rib,
            tol: self.my_tol,
            pln: plane,
            bound: bnd,
            bnd_face: &bnd_face,
            check_pnt,
            first_edge,
            last_edge,
            on_first_face,
            on_last_face,
        };
        let prof = match self.base.no_sliding_profile(&args, concavite) {
            Ok(f) => f,
            Err(_) => {
                self.base.my_status_error = StatusError::NoFaceProf;
                self.base.not_done();
                return Ok(false);
            }
        };

        // `cxx:1077-1089`. The return value of `Propagate` is discarded here,
        // unlike in the sliding branch.
        let mut falseside = true;
        let first_pnt = self.base.my_first_pnt;
        let last_pnt = self.base.my_last_pnt;
        self.propagate(sli_list, &prof, &first_pnt, &last_pnt, &mut falseside)?;
        if !falseside {
            self.base.my_status_error = StatusError::FalseSide;
            self.base.not_done();
            return Ok(false);
        }

        // `cxx:1091-1108`.
        self.my_slface.clear();
        let builder = TopoBuilder::new();
        let mut comp = Shell::new();
        for s in sli_list.iter() {
            builder.add(&mut comp.0, s);
        }
        comp.0.set_closed(BRepTool::is_closed(&comp.0));
        self.base.my_suntil = Some(comp.0.clone());
        self.base.my_skface = Some(prof.clone());
        self.base.my_pbase = Some(prof);

        Ok(true)
    }
}
