//! `BRepFeat_RibSlot::SlidingProfile`
//! (`ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_RibSlot.cxx:1418-1786`).
//!
//! The profile face of a sliding rib is the wire of the sketch extended to the
//! faces it slides on: the first / last point of the sketch is joined to the
//! intersection of the two normal lines of the support faces, or, when those
//! lines are parallel or meet outside the material, to the bounding face of the
//! rib. The `Concavite` out-parameter records which of the three shapes the
//! reference took:
//!
//! * `1` concave -- the two normal lines meet at `Pt` (`cxx:1503`),
//! * `2` parallel and concave -- the last point is joined straight to the first
//!   one (`cxx:1479`),
//! * `3` convex -- the profile is closed against `BndFace` (`cxx:1519-1670`)
//!   and finally kept or replaced by the complementary part of `BndFace`
//!   (`cxx:1748-1774`).
//!
//! Deviations, all forced by the port's data model:
//! * `GeomAPI::To2d` answering a null handle makes the reference dereference it
//!   in the next `Geom2dAPI_InterCurveCurve`; the port answers
//!   [`RibProfileError::NotProjected`].
//! * `cxx:1603-1609` leaves `LastPnt` uninitialized when the last vertex of
//!   `BndEdge1` coincides with `BndPnt1` (the `LinearForm` alternative that
//!   would take `FirstVertex(e)` is commented out); the port keeps that
//!   coinciding point.
//! * `cxx:1458` declares `gp_Pnt Pt` uninitialized and reads it in the
//!   `Concavite == 1` arm; the port answers [`RibProfileError::NotProjected`]
//!   when no intersection set it.
//! * `BRepLib_MakeFace(P, WW, true)` is
//!   [`crate::brep_lib_make_face::make_face_from_plane_and_wire`].

use std::sync::Arc;

use occt_core::elib::{clib, slib};
use occt_core::gp::{GpDir, GpLin, GpPln, GpPnt, GpPnt2d};
use occt_core::precision::CONFUSION;
use occt_geom::curve::Curve;
use occt_geom::geom_api::to_2d;
use occt_geom::line::GeomLine;
use occt_geom2d::geom2d_api_inter_curve_curve::InterCurveCurve2d;

use crate::brep_algo_api;
use crate::brep_check;
use crate::brep_feat::rib_slot::{check_point, make_edges, RibSlot};
use crate::brep_lib_make_face;
use crate::brep_lib_make_wire::MakeWire;
use crate::brep_tool::BRepTool;
use crate::brep_top_adaptor_fclass2d::BrepTopAdaptorFClass2d;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::fclass2d::FaceState;
use crate::shape::{Edge, Face};
use crate::topo_tools::{is_same, wires_of, wires_of_face};

/// The by-value arguments of `SlidingProfile` (`cxx:1418-1431`).
pub struct SlidingProfileArgs<'a> {
    /// `RevolRib`.
    pub revol_rib: bool,
    /// `myTol` (`cxx:1420`).
    pub tol: f64,
    /// `myPln`.
    pub pln: &'a GpPln,
    /// `BndFace`: the bounding face the profile closes against.
    pub bnd_face: &'a Face,
    /// `CheckPnt`: a point slightly inside the material side (`cxx:1754`).
    pub check_pnt: &'a GpPnt,
    /// `FirstFace`: the support face of `myFirstPnt`.
    pub first_face: &'a Face,
    /// `LastFace`: the support face of `myLastPnt`.
    pub last_face: &'a Face,
    /// `FirstEdge`: the sketch edge through `myFirstPnt`.
    pub first_edge: &'a Edge,
    /// `LastEdge`: the sketch edge through `myLastPnt`.
    pub last_edge: &'a Edge,
}

/// The `false` answers of `SlidingProfile`: the reference returns
/// `ProfileOK == false` after each of these.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RibProfileError {
    /// `myWire` is null, or an edge of the sketch wire has no 3D curve.
    NotInitialized,
    /// `GeomAPI::To2d` answered a null handle, or the reference would read an
    /// uninitialized point (see the module doc).
    NotProjected,
    /// `cxx:1570-1579`: `BndEdge1` or `BndEdge2` is null.
    NullBoundingEdge,
    /// `cxx:1736-1744`: `BRepAlgo::IsValid(fac)` is false.
    InvalidFace,
    /// `cxx:1761-1763`: the cut produced no wire, or `BRepLib_MakeFace(W)`
    /// found no surface.
    CutFailed,
    /// `cxx:1776-1783`: `BRepAlgo::IsValid(Prof)` is false.
    InvalidProfile,
}

/// `gp_Dir` of a line the caller built from `gp_Lin`.
fn line_direction(line: &Arc<dyn Curve>) -> Option<GpDir> {
    line.gp_line().map(|l| l.direction())
}

/// `ElCLib::Parameter(ln1->Lin(), p)` (`cxx:1476`) on the built line.
fn line_parameter(line: &Arc<dyn Curve>, p: &GpPnt) -> Option<f64> {
    line.gp_line().map(|l| clib::line_parameter(&l, p))
}

/// `explo.Next()` followed by `if (explo.Current().IsNull()) explo.Init(W);`
/// (`cxx:1617-1626`, `cxx:2092-2102`).
pub(super) fn advance_or_restart(explo: &mut WireExplorer, wire: &crate::shape::Wire) {
    if explo.more() {
        explo.next();
        if explo.current().is_none() {
            *explo = WireExplorer::of_wire(wire);
        }
    } else {
        *explo = WireExplorer::of_wire(wire);
    }
}

impl RibSlot {
    /// `SlidingProfile(...)` (`cxx:1418-1786`).
    ///
    /// `concavite` is the reference's in/out `Concavite`: the caller may pass
    /// `3` in (`cxx:1487`), and the function writes `2` (`cxx:1479`) or `1`
    /// (`cxx:1503`) as it goes. `Ok` is `ProfileOK == true` with the profile
    /// face `Prof`; `Err` is `ProfileOK == false` with `Prof` untouched.
    pub fn sliding_profile(
        &self,
        args: &SlidingProfileArgs<'_>,
        concavite: &mut i32,
    ) -> Result<Face, RibProfileError> {
        // `cxx:1443-1448`.
        let mut fn_dir = check_point::normal(args.first_face, &self.my_first_pnt);
        let mut ln_dir = check_point::normal(args.last_face, &self.my_last_pnt);

        // `cxx:1452-1457`: the groove (cut) works on the material side.
        if !self.my_fuse {
            fn_dir.reverse();
            ln_dir.reverse();
        }

        // `cxx:1462-1463`: `ln2` leaves `myFirstPnt` along `FN`, `ln1` leaves
        // `myLastPnt` along `LN`.
        let ln2: Arc<dyn Curve> =
            Arc::new(GeomLine::new(GpLin::from_pnt_dir(self.my_first_pnt, fn_dir)));
        let ln1: Arc<dyn Curve> =
            Arc::new(GeomLine::new(GpLin::from_pnt_dir(self.my_last_pnt, ln_dir)));

        // `cxx:1465-1467`.
        let ln2d1 = to_2d(ln1.as_ref(), args.pln).ok_or(RibProfileError::NotProjected)?;
        let ln2d2 = to_2d(ln2.as_ref(), args.pln).ok_or(RibProfileError::NotProjected)?;
        let inter = InterCurveCurve2d::with_curves(&ln2d1, &ln2d2, CONFUSION);

        let mut ww = MakeWire::new();
        let mut test_ok = true;

        // `cxx:1470-1496`.
        if args.revol_rib {
            let d1 = line_direction(&ln1).ok_or(RibProfileError::NotProjected)?;
            let d2 = line_direction(&ln2).ok_or(RibProfileError::NotProjected)?;
            if d1.is_opposite_ang(&d2, args.tol) {
                let par1 =
                    line_parameter(&ln1, &self.my_first_pnt).ok_or(RibProfileError::NotProjected)?;
                let par2 =
                    line_parameter(&ln2, &self.my_last_pnt).ok_or(RibProfileError::NotProjected)?;
                if par1 >= args.tol || par2 >= args.tol {
                    // `cxx:1479`: parallel and concave.
                    *concavite = 2;
                    if let Some(e1) = make_edges::edge_from_points(&self.my_last_pnt, &self.my_first_pnt)
                    {
                        let _ = ww.add(&e1);
                    }
                }
            }
            if d1.is_equal_tol(&d2, args.tol) && *concavite == 3 {
                test_ok = false;
            }
        }

        // `cxx:1498-1508`.
        let mut pt: Option<GpPnt> = None;
        if test_ok && inter.nb_points() > 0 {
            if let Some(p2d) = inter.point(1) {
                let p = slib::plane_value(args.pln, p2d.x(), p2d.y());
                let par = check_point::int_par(Some(ln1.as_ref()), &p);
                if par > 0.0 {
                    // `cxx:1503`: concave.
                    *concavite = 1;
                }
                pt = Some(p);
            }
        }

        if *concavite == 1 {
            // `cxx:1510-1518`: extend the first and last edges to `Pt`.
            let pt = pt.ok_or(RibProfileError::NotProjected)?;
            if let Some(e1) = make_edges::edge_from_points(&self.my_last_pnt, &pt) {
                let _ = ww.add(&e1);
            }
            if let Some(e2) = make_edges::edge_from_points(&pt, &self.my_first_pnt) {
                let _ = ww.add(&e2);
            }
        } else if *concavite == 3 {
            self.close_against_bounding_face(args, &ln1, &ln2, &ln2d1, &ln2d2, &mut ww)?;
        }

        // `cxx:1674-1731`: the sketch wire, trimmed to the two extreme points.
        self.append_profile_wire(args, &mut ww)?;

        // `cxx:1733-1746`.
        let fac =
            brep_lib_make_face::make_face_from_plane_and_wire(args.pln, &ww.wire(), true);
        if !brep_check::is_valid(&fac.0) {
            return Err(RibProfileError::InvalidFace);
        }

        // `cxx:1748-1774`.
        let prof = if *concavite != 3 {
            fac
        } else {
            let cl = BrepTopAdaptorFClass2d::new(&fac, BRepTool::face_tolerance(&fac));
            let (u, v) = slib::plane_parameters(&args.pln.position(), args.check_pnt);
            if cl.perform(GpPnt2d::new(u, v), true) == FaceState::Out {
                // `cxx:1761-1768`: the complementary part of `BndFace`.
                let cut = brep_algo_api::cut(&args.bnd_face.0, &fac.0)
                    .map_err(|_| RibProfileError::CutFailed)?;
                let w = wires_of(&cut)
                    .into_iter()
                    .next()
                    .ok_or(RibProfileError::CutFailed)?;
                brep_lib_make_face::make_face_from_wire(&w, false)
                    .ok_or(RibProfileError::CutFailed)?
            } else {
                // `cxx:1770-1773`.
                fac
            }
        };

        // `cxx:1776-1783`.
        if !brep_check::is_valid(&prof.0) {
            return Err(RibProfileError::InvalidProfile);
        }
        Ok(prof)
    }

    /// `cxx:1519-1670`: closes the profile against the bounding face, either
    /// through one edge of it (`BndEdge1.IsSame(BndEdge2)`) or by walking the
    /// bounding wire from `BndEdge1` to `BndEdge2`.
    fn close_against_bounding_face(
        &self,
        args: &SlidingProfileArgs<'_>,
        ln1: &Arc<dyn Curve>,
        ln2: &Arc<dyn Curve>,
        ln2d1: &Arc<dyn occt_geom2d::curve::Curve2d>,
        ln2d2: &Arc<dyn occt_geom2d::curve::Curve2d>,
        ww: &mut MakeWire,
    ) -> Result<(), RibProfileError> {
        // `cxx:1527-1529`: the bounding wire.
        let bnd_wire = wires_of_face(args.bnd_face)
            .into_iter()
            .next()
            .ok_or(RibProfileError::NullBoundingEdge)?;

        let mut bnd_edge1: Option<Edge> = None;
        let mut bnd_edge2: Option<Edge> = None;
        let mut bnd_pnt1 = GpPnt::default();
        let mut bnd_pnt2 = GpPnt::default();

        // `cxx:1530-1568`.
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
                            let parl = check_point::int_par(Some(ln1.as_ref()), &p);
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
                            let parl = check_point::int_par(Some(ln2.as_ref()), &p);
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

        // `cxx:1570-1579`.
        let (Some(bnd_edge1), Some(bnd_edge2)) = (bnd_edge1, bnd_edge2) else {
            return Err(RibProfileError::NullBoundingEdge);
        };

        // `cxx:1581-1582`.
        if let Some(e1) = make_edges::edge_from_points(&self.my_last_pnt, &bnd_pnt1) {
            let _ = ww.add(&e1);
        }

        if is_same(&bnd_edge1.0, &bnd_edge2.0) {
            // `cxx:1584-1590`: same edge, a simply determined path.
            if let Some(e2) = make_edges::edge_from_points(&bnd_pnt1, &bnd_pnt2) {
                let _ = ww.add(&e2);
            }
            if let Some(e3) = make_edges::edge_from_points(&bnd_pnt2, &self.my_first_pnt) {
                let _ = ww.add(&e3);
            }
            return Ok(());
        }

        // `cxx:1592-1615`: the path starts at the end of `BndEdge1`.
        let mut last_pnt = bnd_pnt1;
        explo = WireExplorer::of_wire(&bnd_wire);
        while explo.more() {
            let Some(e) = explo.current_edge() else { break };
            if is_same(&e.0, &bnd_edge1.0) {
                if let Some(lv) = crate::shhealing::last_vertex(&e) {
                    let pp = BRepTool::vertex_point_world(&lv);
                    last_pnt = pp;
                }
                break;
            }
            explo.next();
        }

        // `cxx:1617-1626`.
        advance_or_restart(&mut explo, &bnd_wire);

        // `cxx:1633-1666`: walk to `BndEdge2`, closing on `BndPnt2`.
        let mut fin = false;
        while !fin {
            let Some(e) = explo.current_edge() else {
                return Err(RibProfileError::NullBoundingEdge);
            };
            if !is_same(&e.0, &bnd_edge2.0) {
                let pp = crate::shhealing::last_vertex(&e)
                    .map(|v| BRepTool::vertex_point_world(&v))
                    .unwrap_or_default();
                if let Some(ee) = make_edges::edge_from_points(&last_pnt, &pp) {
                    let _ = ww.add(&ee);
                }
                last_pnt = pp;
            } else {
                // `cxx:1645-1651`: the path is closed.
                fin = true;
                if let Some(ee) = make_edges::edge_from_points(&last_pnt, &bnd_pnt2) {
                    let _ = ww.add(&ee);
                }
                last_pnt = bnd_pnt2;
            }
            // `cxx:1652-1661`.
            advance_or_restart(&mut explo, &bnd_wire);
        }
        // `cxx:1668-1669`: back to the first point of the sketch.
        if let Some(e3) = make_edges::edge_from_points(&bnd_pnt2, &self.my_first_pnt) {
            let _ = ww.add(&e3);
        }
        let _ = last_pnt;
        Ok(())
    }

    /// `cxx:1674-1731`: the sketch wire between `FirstEdge` and `LastEdge`.
    fn append_profile_wire(
        &self,
        args: &SlidingProfileArgs<'_>,
        ww: &mut MakeWire,
    ) -> Result<(), RibProfileError> {
        let my_wire = self
            .my_wire
            .as_ref()
            .ok_or(RibProfileError::NotInitialized)?;
        let mut ex = WireExplorer::of_wire(my_wire);
        let first_curve = BRepTool::edge_curve_world(args.first_edge);

        if !is_same(&args.first_edge.0, &args.last_edge.0) {
            // `cxx:1681-1686`.
            let fl_pnt = crate::shhealing::last_vertex(args.first_edge)
                .map(|v| BRepTool::vertex_point_world(&v))
                .unwrap_or_default();
            if let Some(c) = &first_curve {
                if let Some(ef) =
                    make_edges::edge_from_curve_points(c, &self.my_first_pnt, &fl_pnt)
                {
                    let _ = ww.add(&ef);
                }
            }
            // `cxx:1687-1690`.
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
            // `cxx:1692`.
            ex.next();
            // `cxx:1694-1703`.
            while ex.more() {
                let Some(e) = ex.current_edge() else { break };
                if !is_same(&e.0, &args.last_edge.0) {
                    let _ = ww.add(&e);
                } else {
                    break;
                }
                ex.next();
            }
            // `cxx:1706-1713`.
            let lf_pnt = crate::shhealing::first_vertex(args.last_edge)
                .map(|v| BRepTool::vertex_point_world(&v))
                .unwrap_or_default();
            if let Some(c) = BRepTool::edge_curve_world(args.last_edge) {
                if let Some(el) =
                    make_edges::edge_from_curve_points(&c, &lf_pnt, &self.my_last_pnt)
                {
                    let _ = ww.add(&el);
                }
            }
        } else {
            // `cxx:1716-1731`: a single edge, reversed when the two points lie
            // the other way around on it.
            let c = first_curve.ok_or(RibProfileError::NotInitialized)?;
            let fpar = check_point::int_par(Some(c.as_ref()), &self.my_first_pnt);
            let lpar = check_point::int_par(Some(c.as_ref()), &self.my_last_pnt);
            let curve: Arc<dyn Curve> = if fpar > lpar {
                Arc::from(c.reversed())
            } else {
                c
            };
            if let Some(ef) =
                make_edges::edge_from_curve_points(&curve, &self.my_first_pnt, &self.my_last_pnt)
            {
                let _ = ww.add(&ef);
            }
        }
        Ok(())
    }
}
