//! `BRepFeat_RibSlot::NoSlidingProfile`
//! (`ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_RibSlot.cxx:1793-2667`).
//!
//! The profile face when the sketch does *not* slide on its support faces: the
//! first / last sketch edges are extended to the boundary of the rib (or into
//! the support faces), the two tangent lines at the extended ends are
//! intersected with the bounding face, and the profile wire is rebuilt edge by
//! edge so every added edge starts on the last vertex of the previous one. The
//! rebuilt occurrences are recorded in `myLFMap` (`cxx:2220-2580`), which
//! `LFPerform` / `Generated` read back.
//!
//! `Concavite` follows the same convention as
//! [`super::RibSlot::sliding_profile`]: `1` concave (`cxx:1878`), `2` parallel
//! and concave (`cxx:1860`), `3` convex (`cxx:1886-2196`).
//!
//! Deviations, all forced by the port's data model:
//! * the reference chains edges through shared `TopoDS_Vertex` handles
//!   (`BRepLib_MakeEdge(V1, V2)` / `BRepLib_MakeVertex`); the port chains the
//!   same points and lets `BRepLib_MakeWire::Add` unify the endpoints, which is
//!   what that class does with its own vertices anyway.
//! * `cxx:2066-2068` leaves `LastPnt` uninitialized when the last vertex of
//!   `BndEdge1` coincides with `BndPnt1`; the port keeps that coinciding point.
//! * `GeomAPI::To2d` answering a null handle makes the reference dereference it
//!   in the next `Geom2dAPI_InterCurveCurve`; the port answers
//!   [`RibProfileError::NotProjected`].
//! * `cxx:1841` declares `gp_Pnt Pt` uninitialized and reads it in the
//!   `Concavite == 1` arms; the port answers [`RibProfileError::NotProjected`]
//!   when no intersection set it.
//! * `cxx:2479-2527` is the `FirstEdge.IsSame(LastEdge)` alternative inside the
//!   `FirstEdge != LastEdge` arm and is unreachable, so it is not reproduced.
//! * the `BRepAlgoAPI_Cut` history query of `UpdateDescendants` (`cxx:2644`)
//!   goes through [`crate::brep_feat::algo_api::BooleanOperation::modified`].

use std::sync::Arc;

use occt_core::elib::{clib, slib};
use occt_core::gp::{GpDir, GpLin, GpPln, GpPnt, GpPnt2d, GpVec};
use occt_core::precision::CONFUSION;
use occt_geom::curve::Curve;
use occt_geom::geom_api::to_2d;
use occt_geom::line::GeomLine;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::geom2d_api_inter_curve_curve::InterCurveCurve2d;

use crate::brep_check;
use crate::brep_feat::algo_api;
use crate::brep_feat::rib_slot::check_point;
use crate::brep_feat::rib_slot::edge_extension;
use crate::brep_feat::rib_slot::make_edges;
use crate::brep_feat::rib_slot::sliding_profile::{advance_or_restart, RibProfileError};
use crate::brep_feat::rib_slot::RibSlot;
use crate::brep_lib_make_face;
use crate::brep_lib_make_wire::MakeWire;
use crate::brep_tool::BRepTool;
use crate::brep_top_adaptor_fclass2d::BrepTopAdaptorFClass2d;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::fclass2d::FaceState;
use crate::shape::{Edge, Face, TopoShape, Wire};
use crate::topo_tools::{edge_vertices, is_same, wires_of, wires_of_face};

/// The by-value arguments of `NoSlidingProfile` (`cxx:1793-1808`).
pub struct NoSlidingProfileArgs<'a> {
    /// `RevolRib`.
    pub revol_rib: bool,
    /// `myTol` (`cxx:1795`).
    pub tol: f64,
    /// `myPln`.
    pub pln: &'a GpPln,
    /// `bnd`: the bounding length `HeightMax` answers.
    pub bound: f64,
    /// `BndFace`: the bounding face the profile closes against.
    pub bnd_face: &'a Face,
    /// `CheckPnt`: a point slightly inside the material side (`cxx:2637`).
    pub check_pnt: &'a GpPnt,
    /// `FirstEdge`: the sketch edge through `myFirstPnt`.
    pub first_edge: &'a Edge,
    /// `LastEdge`: the sketch edge through `myLastPnt`.
    pub last_edge: &'a Edge,
    /// `OnFirstFace`.
    pub on_first_face: bool,
    /// `OnLastFace`.
    pub on_last_face: bool,
}

/// `TopExp::FirstVertex(E, true)` as a point.
fn cum_first_pnt(e: &Edge) -> Option<GpPnt> {
    crate::shhealing::first_vertex(e).map(|v| BRepTool::vertex_point_world(&v))
}

/// `TopExp::LastVertex(E, true)` as a point.
fn cum_last_pnt(e: &Edge) -> Option<GpPnt> {
    crate::shhealing::last_vertex(e).map(|v| BRepTool::vertex_point_world(&v))
}

/// `TopExp::LastVertex(E)` (`CumOri = false`) as a point: the stored last child.
fn stored_last_pnt(e: &Edge) -> Option<GpPnt> {
    edge_vertices(e)
        .1
        .map(|v| BRepTool::vertex_point_world(&v))
}

/// The `gp_Dir` of a `gp_Vec` (`gp_Dir::gp_Dir(const gp_Vec&)`).
fn dir_of(v: &GpVec) -> Option<GpDir> {
    GpDir::from_vec(v).ok()
}

/// The wire being rebuilt, with the reference's `theFV` / `theFirstpoint` /
/// `theLastEdge` bookkeeping (`cxx:1818-1830`).
struct Chain {
    /// `w`, built with `BRep_Builder::MakeWire`.
    w: MakeWire,
    /// `theFirstpoint`: `BRep_Tool::Pnt(theFV)`, set by the first added edge.
    first_pnt: Option<GpPnt>,
    /// `theLastEdge`: the edge added last.
    last: Option<Edge>,
}

impl Chain {
    fn new() -> Self {
        Self {
            w: MakeWire::new(),
            first_pnt: None,
            last: None,
        }
    }

    /// `TopExp::LastVertex(theLastEdge, true)`.
    fn last_pnt(&self) -> Option<GpPnt> {
        self.last.as_ref().and_then(cum_last_pnt)
    }

    /// `if (!theFV.IsNull() && p.Distance(theFirstpoint) <= myTol) v2 = theFV;`
    /// -- applied to points instead of vertices.
    fn close_end(&self, p: &GpPnt, tol: f64) -> GpPnt {
        match &self.first_pnt {
            Some(f) if p.distance(f) <= tol => *f,
            _ => *p,
        }
    }

    /// `BRepLib_MakeEdge(a, b)` when `theLastEdge` is null, otherwise
    /// `BRepLib_MakeEdge(LastVertex(theLastEdge, true), MakeVertex(b))`.
    fn segment(&self, a: &GpPnt, b: &GpPnt) -> Option<Edge> {
        match self.last_pnt() {
            Some(from) => make_edges::edge_from_points(&from, b),
            None => make_edges::edge_from_points(a, b),
        }
    }

    /// The `BRepLib_MakeEdge(cc, V1, V2)` form of [`Self::segment`].
    fn curve_segment(&self, cc: &Arc<dyn Curve>, a: &GpPnt, b: &GpPnt) -> Option<Edge> {
        match self.last_pnt() {
            Some(from) => make_edges::edge_from_curve_points(cc, &from, b),
            None => make_edges::edge_from_curve_points(cc, a, b),
        }
    }

    /// `BB.Add(w, e)` plus the `theFV` / `theFirstpoint` / `theLastEdge` writes.
    fn push(&mut self, e: Option<Edge>) {
        let Some(e) = e else {
            // `theLastEdge = TopoDS::Edge(e.Shape())` of a failed `MakeEdge`.
            self.last = None;
            return;
        };
        let _ = self.w.add(&e);
        if self.first_pnt.is_none() {
            self.first_pnt = cum_first_pnt(&e);
        }
        self.last = Some(e);
    }

    /// `BB.Add(w, e)` without the `theLastEdge` update (`cxx:2612-2616`).
    fn push_raw(&mut self, e: Option<Edge>) {
        if let Some(e) = e {
            let _ = self.w.add(&e);
        }
    }
}

/// `cxx:1933-1991`: the first and last edge of the bounding wire that the two
/// tangent lines hit, with the hit points.
#[allow(clippy::type_complexity)]
fn bounding_wire_hits(
    args: &NoSlidingProfileArgs<'_>,
    firstln: &Arc<dyn Curve>,
    lastln: &Arc<dyn Curve>,
    ln2d1: &Arc<dyn Curve2d>,
    ln2d2: &Arc<dyn Curve2d>,
) -> Result<(Wire, Edge, Edge, GpPnt, GpPnt), RibProfileError> {
    // `cxx:1936-1940`.
    let bnd_wire = wires_of_face(args.bnd_face)
        .into_iter()
        .next()
        .ok_or(RibProfileError::NullBoundingEdge)?;

    let mut bnd_edge1: Option<Edge> = None;
    let mut bnd_edge2: Option<Edge> = None;
    let mut bnd_pnt1 = GpPnt::default();
    let mut bnd_pnt2 = GpPnt::default();

    let mut explo = WireExplorer::of_wire(&bnd_wire);
    while explo.more() {
        let Some(e) = explo.current_edge() else { break };
        let (first, last) = BRepTool::edge_parameters(&e);
        if let Some(c) = BRepTool::edge_curve_world(&e) {
            if let Some(c2d) = to_2d(c.as_ref(), args.pln) {
                let intcln1 = InterCurveCurve2d::with_curves(ln2d1, &c2d, CONFUSION);
                if intcln1.nb_points() > 0 {
                    if let Some(p2d) = intcln1.point(1) {
                        let p = slib::plane_value(args.pln, p2d.x(), p2d.y());
                        let parl = check_point::int_par(Some(firstln.as_ref()), &p);
                        let parc = check_point::int_par(Some(c.as_ref()), &p);
                        if parc >= first && parc <= last && parl >= 0.0 {
                            bnd_edge1 = Some(e.clone());
                            bnd_pnt1 = p;
                        }
                    }
                }
                let intcln2 = InterCurveCurve2d::with_curves(ln2d2, &c2d, CONFUSION);
                if intcln2.nb_points() > 0 {
                    if let Some(p2d) = intcln2.point(1) {
                        let p = slib::plane_value(args.pln, p2d.x(), p2d.y());
                        let parl = check_point::int_par(Some(lastln.as_ref()), &p);
                        let parc = check_point::int_par(Some(c.as_ref()), &p);
                        if parc >= first && parc <= last && parl >= 0.0 {
                            bnd_edge2 = Some(e.clone());
                            bnd_pnt2 = p;
                        }
                    }
                }
            }
        }
        if bnd_edge1.is_some() && bnd_edge2.is_some() {
            break;
        }
        explo.next();
    }

    // `cxx:1982-1990`.
    let (Some(a), Some(b)) = (bnd_edge1, bnd_edge2) else {
        return Err(RibProfileError::NullBoundingEdge);
    };
    Ok((bnd_wire, a, b, bnd_pnt1, bnd_pnt2))
}

/// `cxx:2056-2196`: walk the bounding wire from `BndEdge1` to `BndEdge2`.
fn walk_bounding_wire(
    chain: &mut Chain,
    bnd_wire: &Wire,
    bnd_edge1: &Edge,
    bnd_edge2: &Edge,
    bnd_pnt1: &GpPnt,
    bnd_pnt2: &GpPnt,
) -> Result<(), RibProfileError> {
    let mut last_pnt = *bnd_pnt1;

    // `cxx:2058-2090`.
    let mut explo = WireExplorer::of_wire(bnd_wire);
    while explo.more() {
        let Some(e) = explo.current_edge() else { break };
        if is_same(&e.0, &bnd_edge1.0) {
            if let Some(lv) = crate::shhealing::last_vertex(&e) {
                let pp = BRepTool::vertex_point_world(&lv);
                // `cxx:2066-2069`. The reference leaves `LastPnt` uninitialized
                // in the other arm; the port keeps the coinciding point.
                if pp.distance(bnd_pnt1) > BRepTool::edge_tolerance(&e) {
                    last_pnt = pp;
                }
            }
            chain.push(chain.segment(bnd_pnt1, &last_pnt));
            break;
        }
        explo.next();
    }

    // `cxx:2092-2102`.
    advance_or_restart(&mut explo, bnd_wire);

    // `cxx:2103-2175`.
    let mut fin = false;
    while !fin {
        let Some(e) = explo.current_edge() else {
            return Err(RibProfileError::NullBoundingEdge);
        };
        if !is_same(&e.0, &bnd_edge2.0) {
            // `cxx:2106-2133`.
            let pp = cum_last_pnt(&e).unwrap_or_default();
            chain.push(chain.segment(&last_pnt, &pp));
            last_pnt = pp;
        } else {
            // `cxx:2135-2158`: the path is closed.
            fin = true;
            chain.push(chain.segment(&last_pnt, bnd_pnt2));
            last_pnt = *bnd_pnt2;
        }
        // `cxx:2160-2174`.
        advance_or_restart(&mut explo, bnd_wire);
    }
    let _ = last_pnt;
    Ok(())
}

/// `cxx:2220-2580`: the rebuilt sketch wire and its `myLFMap` occurrences.
#[allow(clippy::too_many_arguments)]
fn rebuild_profile_wire(
    slot: &mut RibSlot,
    args: &NoSlidingProfileArgs<'_>,
    chain: &mut Chain,
    false_first_edge: &Option<Edge>,
    false_last_edge: &Option<Edge>,
    false_only_one: &Option<Edge>,
    same_edges: bool,
) -> Result<(), RibProfileError> {
    let my_wire = slot.my_wire.clone().ok_or(RibProfileError::NotInitialized)?;

    if same_edges {
        slot.bind_lf_map(&args.first_edge.0);
        if args.on_first_face || args.on_last_face {
            // `cxx:2230-2274`.
            let only = false_only_one
                .as_ref()
                .ok_or(RibProfileError::NotInitialized)?;
            let cc = BRepTool::edge_curve_world(only).ok_or(RibProfileError::NotInitialized)?;
            let pp = cum_last_pnt(only).unwrap_or_default();
            let end = chain.close_end(&pp, args.tol);
            let start = chain
                .last_pnt()
                .or_else(|| cum_first_pnt(only))
                .unwrap_or_default();
            let e = chain.curve_segment(&cc, &start, &end);
            if let Some(edge) = e.clone() {
                slot.append_lf(&args.first_edge.0, &edge.0);
            }
            chain.push(e);
        } else {
            // `cxx:2275-2317`. `TopExp::LastVertex(FirstEdge)` without `CumOri`.
            let cc =
                BRepTool::edge_curve_world(args.first_edge).ok_or(RibProfileError::NotInitialized)?;
            let pp = if chain.last.is_some() {
                stored_last_pnt(args.first_edge)
            } else {
                cum_last_pnt(args.first_edge)
            }
            .unwrap_or_default();
            let end = chain.close_end(&pp, args.tol);
            let start = chain
                .last_pnt()
                .or_else(|| cum_first_pnt(args.first_edge))
                .unwrap_or_default();
            let e = chain.curve_segment(&cc, &start, &end);
            if let Some(edge) = e.clone() {
                slot.append_lf(&args.first_edge.0, &edge.0);
            }
            chain.push(e);
        }
        return Ok(());
    }

    // `cxx:2320-2580`: `FirstEdge != LastEdge`.
    slot.bind_lf_map(&args.first_edge.0);
    if !args.on_first_face {
        // `cxx:2329-2350`.
        if chain.last.is_some() {
            let cc = BRepTool::edge_curve_world(args.first_edge)
                .ok_or(RibProfileError::NotInitialized)?;
            let start = chain.last_pnt().ok_or(RibProfileError::NotInitialized)?;
            let end = cum_last_pnt(args.first_edge).unwrap_or_default();
            let e = make_edges::edge_from_curve_points(&cc, &start, &end);
            if let Some(edge) = e.clone() {
                slot.append_lf(&args.first_edge.0, &edge.0);
            }
            chain.push(e);
        } else {
            slot.append_lf(&args.first_edge.0, &args.first_edge.0);
            chain.push(Some(args.first_edge.clone()));
        }
    } else {
        // `cxx:2352-2375`: the extended first edge.
        let ff = false_first_edge
            .as_ref()
            .ok_or(RibProfileError::NotInitialized)?;
        if chain.last.is_some() {
            let cc = BRepTool::edge_curve_world(ff).ok_or(RibProfileError::NotInitialized)?;
            let start = chain.last_pnt().ok_or(RibProfileError::NotInitialized)?;
            let end = cum_last_pnt(ff).unwrap_or_default();
            let e = make_edges::edge_from_curve_points(&cc, &start, &end);
            if let Some(edge) = e.clone() {
                slot.append_lf(&args.first_edge.0, &edge.0);
            }
            chain.push(e);
        } else {
            slot.append_lf(&args.first_edge.0, &ff.0);
            chain.push(Some(ff.clone()));
        }
    }

    // `cxx:2377-2423`: the middle of the sketch wire.
    let mut ex = WireExplorer::of_wire(&my_wire);
    while ex.more() {
        if ex
            .current()
            .map(|e| is_same(&e.0, &args.first_edge.0))
            .unwrap_or(false)
        {
            break;
        }
        ex.next();
    }
    ex.next();
    while ex.more() {
        let Some(e) = ex.current_edge() else { break };
        if is_same(&e.0, &args.last_edge.0) {
            break;
        }
        slot.bind_lf_map(&e.0);
        if chain.last.is_some() {
            let cc = BRepTool::edge_curve_world(&e).ok_or(RibProfileError::NotInitialized)?;
            let start = chain.last_pnt().ok_or(RibProfileError::NotInitialized)?;
            let end = cum_last_pnt(&e).unwrap_or_default();
            let rebuilt = make_edges::edge_from_curve_points(&cc, &start, &end);
            if let Some(edge) = rebuilt.clone() {
                slot.append_lf(&e.0, &edge.0);
            }
            chain.push(rebuilt);
        } else {
            slot.append_lf(&e.0, &e.0);
            chain.push(Some(e));
        }
        ex.next();
    }

    // `cxx:2425-2478`: the last edge of the sketch. The
    // `FirstEdge.IsSame(LastEdge)` alternative (`cxx:2479-2527`) is unreachable
    // here.
    if !args.on_last_face {
        let edg = ex
            .current_edge()
            .ok_or(RibProfileError::NullBoundingEdge)?;
        slot.bind_lf_map(&edg.0);
        let cc = BRepTool::edge_curve_world(&edg).ok_or(RibProfileError::NotInitialized)?;
        let pp = if chain.last.is_some() {
            stored_last_pnt(&edg)
        } else {
            cum_last_pnt(&edg)
        }
        .unwrap_or_default();
        let end = chain.close_end(&pp, args.tol);
        let start = chain
            .last_pnt()
            .or_else(|| cum_first_pnt(&edg))
            .unwrap_or_default();
        let e = make_edges::edge_from_curve_points(&cc, &start, &end);
        if let Some(edge) = e.clone() {
            slot.append_lf(&edg.0, &edge.0);
        }
        chain.push(e);
    } else {
        // `cxx:2531-2578`: the extended last edge.
        let fl = false_last_edge
            .as_ref()
            .ok_or(RibProfileError::NotInitialized)?;
        let cc = BRepTool::edge_curve_world(fl).ok_or(RibProfileError::NotInitialized)?;
        let pp = cum_last_pnt(fl).unwrap_or_default();
        let end = chain.close_end(&pp, args.tol);
        let start = chain
            .last_pnt()
            .or_else(|| cum_first_pnt(fl))
            .unwrap_or_default();
        let e = make_edges::edge_from_curve_points(&cc, &start, &end);
        if let Some(edge) = e.clone() {
            slot.append_lf(&args.last_edge.0, &edge.0);
        }
        chain.push(e);
    }
    Ok(())
}

impl RibSlot {
    /// `NoSlidingProfile(...)` (`cxx:1793-2667`).
    ///
    /// `concavite` is the reference's in/out `Concavite` (`cxx:1848-1883`).
    /// `Ok` is `ProfileOK == true` with `Prof`; `Err` is `ProfileOK == false`
    /// with `Prof` untouched. `myLFMap` is written as the reference does.
    pub fn no_sliding_profile(
        &mut self,
        args: &NoSlidingProfileArgs<'_>,
        concavite: &mut i32,
    ) -> Result<Face, RibProfileError> {
        // `cxx:1832-1834`.
        let (f1, _l1) = BRepTool::edge_parameters(args.first_edge);
        let (_, l2) = BRepTool::edge_parameters(args.last_edge);
        let first_curve =
            BRepTool::edge_curve_world(args.first_edge).ok_or(RibProfileError::NotInitialized)?;
        let last_curve =
            BRepTool::edge_curve_world(args.last_edge).ok_or(RibProfileError::NotInitialized)?;

        // `cxx:1837-1839`.
        let (mut firstpoint, firstvect) = first_curve.d1(f1);
        let (mut lastpoint, mut lastvect) = last_curve.d1(l2);
        let mut neg_firstvect = firstvect;
        neg_firstvect.reverse();
        let mut lastln: Arc<dyn Curve> = Arc::new(GeomLine::new(GpLin::from_pnt_dir(
            firstpoint,
            dir_of(&neg_firstvect).ok_or(RibProfileError::NotProjected)?,
        )));
        let mut firstln: Arc<dyn Curve> = Arc::new(GeomLine::new(GpLin::from_pnt_dir(
            lastpoint,
            dir_of(&lastvect).ok_or(RibProfileError::NotProjected)?,
        )));

        // `cxx:1843-1847`.
        let mut ln2d1 = to_2d(firstln.as_ref(), args.pln).ok_or(RibProfileError::NotProjected)?;
        let mut ln2d2 = to_2d(lastln.as_ref(), args.pln).ok_or(RibProfileError::NotProjected)?;
        let inter = InterCurveCurve2d::with_curves(&ln2d1, &ln2d2, CONFUSION);

        // `cxx:1848-1871`.
        let mut test_ok = true;
        if args.revol_rib {
            let d1 = firstln
                .gp_line()
                .map(|l| l.direction())
                .ok_or(RibProfileError::NotProjected)?;
            let d2 = lastln
                .gp_line()
                .map(|l| l.direction())
                .ok_or(RibProfileError::NotProjected)?;
            if d1.is_opposite_ang(&d2, args.tol) {
                let par1 = firstln
                    .gp_line()
                    .map(|l| clib::line_parameter(&l, &self.my_first_pnt))
                    .ok_or(RibProfileError::NotProjected)?;
                let par2 = lastln
                    .gp_line()
                    .map(|l| clib::line_parameter(&l, &self.my_last_pnt))
                    .ok_or(RibProfileError::NotProjected)?;
                if par1 >= args.tol || par2 >= args.tol {
                    // `cxx:1860`: parallel and concave.
                    *concavite = 2;
                }
            }
            if d1.is_equal_tol(&d2, args.tol) && *concavite == 3 {
                test_ok = false;
            }
        }

        // `cxx:1873-1884`.
        let mut pt: Option<GpPnt> = None;
        if test_ok && inter.nb_points() > 0 {
            if let Some(p2d) = inter.point(1) {
                let p = slib::plane_value(args.pln, p2d.x(), p2d.y());
                let par = check_point::int_par(Some(firstln.as_ref()), &p);
                if par > 0.0 {
                    // `cxx:1878`: concave.
                    *concavite = 1;
                }
                pt = Some(p);
            }
        }

        let same_edges = is_same(&args.first_edge.0, &args.last_edge.0);
        let mut false_first_edge: Option<Edge> = None;
        let mut false_last_edge: Option<Edge> = None;
        let mut false_only_one: Option<Edge> = None;
        let mut chain = Chain::new();

        // `cxx:1886-2196`: the convex arm.
        if *concavite == 3 {
            if args.on_first_face {
                // `cxx:1891-1905`.
                let mut ff = args.first_edge.clone();
                if let Some(ext) = edge_extension(&ff, args.bound, true) {
                    ff = ext;
                }
                false_first_edge = Some(ff.clone());
                if same_edges {
                    false_only_one = Some(ff.clone());
                }
                let (cf, _lf) = BRepTool::edge_parameters(&ff);
                if let Some(cc) = BRepTool::edge_curve_world(&ff) {
                    let (p, v) = cc.d1(cf);
                    firstpoint = p;
                    let mut neg = v;
                    neg.reverse();
                    lastln = Arc::new(GeomLine::new(GpLin::from_pnt_dir(
                        firstpoint,
                        dir_of(&neg).ok_or(RibProfileError::NotProjected)?,
                    )));
                    ln2d2 =
                        to_2d(lastln.as_ref(), args.pln).ok_or(RibProfileError::NotProjected)?;
                }
            }
            if args.on_last_face {
                // `cxx:1906-1931`.
                let mut fl = if !same_edges {
                    args.last_edge.clone()
                } else {
                    match false_only_one.clone() {
                        Some(e) => e,
                        None => {
                            false_only_one = Some(args.last_edge.clone());
                            args.last_edge.clone()
                        }
                    }
                };
                if let Some(ext) = edge_extension(&fl, args.bound, false) {
                    fl = ext;
                }
                if same_edges {
                    false_only_one = Some(fl.clone());
                }
                false_last_edge = Some(fl.clone());
                let (_ff, lf) = BRepTool::edge_parameters(&fl);
                if let Some(cc) = BRepTool::edge_curve_world(&fl) {
                    let (p, v) = cc.d1(lf);
                    lastpoint = p;
                    lastvect = v;
                    firstln = Arc::new(GeomLine::new(GpLin::from_pnt_dir(
                        lastpoint,
                        dir_of(&lastvect).ok_or(RibProfileError::NotProjected)?,
                    )));
                    ln2d1 =
                        to_2d(firstln.as_ref(), args.pln).ok_or(RibProfileError::NotProjected)?;
                }
            }

            // `cxx:1933-1991`.
            let (bnd_wire, bnd_edge1, bnd_edge2, bnd_pnt1, bnd_pnt2) =
                bounding_wire_hits(args, &firstln, &lastln, &ln2d1, &ln2d2)?;

            // `cxx:1992-2007`.
            chain.push(chain.segment(&lastpoint, &bnd_pnt1));

            // `cxx:2009-2196`.
            if is_same(&bnd_edge1.0, &bnd_edge2.0) {
                chain.push(chain.segment(&bnd_pnt1, &bnd_pnt2));
                let end = chain.close_end(&firstpoint, args.tol);
                chain.push(chain.segment(&bnd_pnt2, &end));
            } else {
                walk_bounding_wire(
                    &mut chain,
                    &bnd_wire,
                    &bnd_edge1,
                    &bnd_edge2,
                    &bnd_pnt1,
                    &bnd_pnt2,
                )?;
                let end = chain.close_end(&firstpoint, args.tol);
                chain.push(chain.segment(&bnd_pnt2, &end));
            }
        }

        // `cxx:2198-2218`: the concave arm joins `Pt` back to `firstpoint`.
        if *concavite == 1 {
            let pt = pt.as_ref().ok_or(RibProfileError::NotProjected)?;
            chain.push(chain.segment(pt, &firstpoint));
        }

        // `cxx:2220-2580`.
        rebuild_profile_wire(
            self,
            args,
            &mut chain,
            &false_first_edge,
            &false_last_edge,
            &false_only_one,
            same_edges,
        )?;

        // `cxx:2582-2608`: the last edge back to `Pt`.
        if *concavite == 1 {
            let pt = pt.as_ref().ok_or(RibProfileError::NotProjected)?;
            let end = chain.close_end(pt, args.tol);
            chain.push(chain.segment(&lastpoint, &end));
        }

        // `cxx:2610-2616`: the parallel concave closing edge, added without
        // touching `theLastEdge`.
        if *concavite == 2 {
            chain.push_raw(make_edges::edge_from_points(&lastpoint, &firstpoint));
        }

        // `cxx:2617-2618`.
        let w = chain.w.wire();
        w.0.set_closed(BRepTool::is_closed(&w.0));

        // `cxx:2619-2631`.
        let fac = brep_lib_make_face::make_face_from_plane_and_wire(args.pln, &w, true);
        if !brep_check::is_valid(&fac.0) {
            return Err(RibProfileError::InvalidFace);
        }

        // `cxx:2633-2656`.
        let prof = if *concavite == 3 {
            let cl = BrepTopAdaptorFClass2d::new(&fac, BRepTool::face_tolerance(&fac));
            let (u, v) = slib::plane_parameters(&args.pln.position(), args.check_pnt);
            if cl.perform(GpPnt2d::new(u, v), true) == FaceState::Out {
                let mut bop = algo_api::cut(&args.bnd_face.0, &fac.0);
                bop.build().map_err(|_| RibProfileError::CutFailed)?;
                let shape = bop.shape().cloned().ok_or(RibProfileError::CutFailed)?;
                // `cxx:2644`: `UpdateDescendants(c, c.Shape(), false)`.
                self.update_descendants_bop(&bop, &shape, false);
                let cut_wire = wires_of(&shape)
                    .into_iter()
                    .next()
                    .ok_or(RibProfileError::CutFailed)?;
                brep_lib_make_face::make_face_from_wire(&cut_wire, false)
                    .ok_or(RibProfileError::CutFailed)?
            } else {
                fac
            }
        } else {
            fac
        };

        // `cxx:2658-2666`.
        if !brep_check::is_valid(&prof.0) {
            return Err(RibProfileError::InvalidProfile);
        }
        Ok(prof)
    }

    /// `if (!myLFMap.IsBound(K)) myLFMap.Bind(K, empty)` (`cxx:2224-2228`).
    fn bind_lf_map(&mut self, key: &TopoShape) {
        if !self.my_lf_map.is_bound(key) {
            self.my_lf_map.bind(key, Vec::new());
        }
    }

    /// `myLFMap(K).Append(e)` (`cxx:2265`).
    fn append_lf(&mut self, key: &TopoShape, e: &TopoShape) {
        if let Some(slot) = self.my_lf_map.change(key) {
            slot.push(e.clone());
        }
    }
}
