//! `Draft_Modification::Perform`, third block: "Calculate new vertices"
//! (`Draft_Modification_1.cxx:1447-1638`) plus `static Choose`
//! (`:2047-2177`).
//!
//! For every vertex the reference picks a "reference" edge (`Choose`) whose new
//! curve and one neighbouring face surface are intersected; the intersection
//! point nearest the old vertex becomes the new vertex, and every edge at the
//! vertex gets its new parameter from that point. When no such edge exists the
//! new point is the projection of the old one onto the first two edges' curves
//! (`cxx:1459-1541`).
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | vertex loop head, `Choose` call | `_1.cxx:1450-1456` | [`DraftModification::perform_vertices`] |
//! | "no concerted edge": alignment of two edges | `_1.cxx:1459-1541` | `align_vertex_edges` |
//! | `IntCurveSurface_HInter::Perform(HAC, HAS)` | `_1.cxx:1545-1556` | `solve_vertex` |
//! | `Extrema_ExtCS` fallback when there is no point | `_1.cxx:1558-1585` | `solve_vertex` |
//! | nearest intersection point | `_1.cxx:1587-1602` | `solve_vertex` |
//! | per-edge parameter, `SetNewGeometry` on a trim | `_1.cxx:1604-1636` | `solve_vertex` |
//! | `Choose` | `_1.cxx:2047-2177` | [`choose`] |
//!
//! # Translation notes
//!
//! * `HAC->Load(Einf.Geometry())` / `HAS->Load(surface)` (`cxx:1447-1448`,
//!   `:2088`, `:2155-2172`) are passed straight through: the port's
//!   `perform_curve_surface` takes the curve and the surface instead of the
//!   adaptors, and reads the same own ranges they report. `GeomAdaptor_Surface`
//!   reports `RealLast()` for an unbounded direction, so [`adaptor_curve_range`]
//!   and [`adaptor_surface_range`] map infinity onto that same value; the
//!   `Extrema_ExtCS(AC, AS, ...)` fallback (`cxx:1561`) is built over exactly
//!   those windows (`ExtCS::with_curve_surface_range`), like the reference
//!   adaptors rather than the surfaces' own parameter spaces.
//! * `IntCurveSurface_HInter::IsDone()` false (`cxx:1546`) is the `Err` of
//!   [`perform_curve_surface`]; both mean "no usable intersection was
//!   computed".
//! * `GeomAPI_ProjectPointOnCurve::NearestPoint` (`cxx:1471-1472`) never tests
//!   `IsDone()`; the port panics with the `StdFail_NotDone` name that read
//!   raises.
//! * `BRep_Tool::Continuity(E, F1, F2)` (`cxx:2067`, `:2118`) is
//!   [`stored_continuity`] and `GeomAbs_G1` / `GeomAbs_C0` are this port's `1`
//!   / `0` (`brep_lib_continuity`).
//! * `Edg.IsSame(Eref)` (`cxx:2114`) is [`is_same`].
//! * `Vinf.ChangeParameter(Edg)` and `myEMap.ChangeFromKey(Edg)` are borrowed
//!   from their own maps, so the port keeps them apart instead of routing the
//!   writes through `self`.
//! * A face record without its surface (a null handle in the reference, which
//!   would be dereferenced at `cxx:2155`) answers `Draft_VertexRecomputation`:
//!   the port cannot build the adaptor, no OCCT rule is added.

use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom::curve::Curve;
use occt_geom::extrema_cs::ExtCS;
use occt_geom::geom_api::project_point_on_curve;
use occt_geom::Surface;

use crate::brep_lib_continuity::stored_continuity;
use crate::brep_tool::BRepTool;
use crate::draft::edge_info::DraftEdgeInfo;
use crate::draft::error_status::DraftErrorStatus;
use crate::draft::face_info::DraftFaceInfo;
use crate::draft::modification::shape_data_map::ShapeDataMap;
use crate::draft::vertex_info::DraftVertexInfo;
use crate::intcurvesurface::perform_curve_surface;
use crate::shape::{Edge, Face, Vertex};
use crate::topo_tools::is_same;

use super::parameters::{parameter, smart_parameter};
use super::DraftModification;

/// `RealLast()`, the value `GeomAdaptor_Curve` / `GeomAdaptor_Surface` answer
/// for an unbounded direction. The port's `HInter` and `Extrema` drivers need a
/// finite window, so this stands in for infinity.
const UNBOUNDED_WINDOW: f64 = 1.0e100;

/// The `(first, last)` an adaptor reports on `c`.
fn adaptor_curve_range(c: &dyn Curve) -> (f64, f64) {
    let first = c.first_parameter();
    let last = c.last_parameter();
    let map = |v: f64| {
        if v.is_finite() {
            v
        } else if v > 0.0 {
            UNBOUNDED_WINDOW
        } else {
            -UNBOUNDED_WINDOW
        }
    };
    (map(first), map(last))
}

/// The `(u0, u1, v0, v1)` an adaptor reports on `s`.
fn adaptor_surface_range(s: &dyn Surface) -> (f64, f64, f64, f64) {
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    let map = |v: f64| {
        if v.is_finite() {
            v
        } else if v > 0.0 {
            UNBOUNDED_WINDOW
        } else {
            -UNBOUNDED_WINDOW
        }
    };
    (map(u0), map(u1), map(v0), map(v1))
}

/// `M.FindFromKey(E)` for a map `Perform` filled itself (`Standard_NoSuchObject`
/// when it did not).
fn edge_info<'a>(e_map: &'a ShapeDataMap<DraftEdgeInfo>, e: &Edge) -> &'a DraftEdgeInfo {
    e_map
        .find_from_key(&e.0)
        .unwrap_or_else(|| panic!("Standard_NoSuchObject: Draft_Modification::Perform"))
}

/// The `myFMap.FindFromKey(F).Geometry()` of a face record.
fn face_geometry(
    f_map: &ShapeDataMap<DraftFaceInfo>,
    f: &Face,
) -> Option<Arc<dyn Surface>> {
    f_map.find_from_key(&f.0).and_then(|fi| fi.geometry()).cloned()
}

impl DraftModification {
    /// `Draft_Modification_1.cxx:1450-1637`.
    pub(super) fn perform_vertices(&mut self) -> Result<(), ()> {
        for ii in 1..=self.my_v_map.extent() {
            let tvv_shape = self.my_v_map.find_key(ii).clone();
            let Some(tvv) = Vertex::wrap(tvv_shape.clone()) else {
                continue;
            };

            // `cxx:1456-1458`
            let chosen = {
                let vinf = self.my_v_map.change_from_index(ii);
                choose(&self.my_f_map, &mut self.my_e_map, &tvv, vinf)
            };

            match chosen {
                Some((ac, as_)) => self.solve_vertex(ii, &tvv, &ac, &as_)?,
                // `cxx:1459-1540`: no concerted edge, align the two edges.
                None => self.align_vertex_edges(ii, &tvv)?,
            }
        }
        Ok(())
    }

    /// `cxx:1459-1541`.
    fn align_vertex_edges(&mut self, ii: usize, tvv: &Vertex) -> Result<(), ()> {
        // `cxx:1461-1462`
        {
            let vinf = self.my_v_map.change_from_index(ii);
            *vinf.change_geometry() = GpPnt::new(0.0, 0.0, 0.0);
            vinf.init_edge_iterator();
        }
        // `cxx:1463-1464`
        if !self.my_v_map.find_from_index(ii).more_edge() {
            return self.fail_vertex(tvv);
        }

        let vtori = BRepTool::vertex_point(tvv);
        let edg1 = self.my_v_map.find_from_index(ii).edge();
        let Some(c1) = self.my_e_map.find_from_key(&edg1.0).and_then(|e| e.geometry()).cloned()
        else {
            return self.fail_vertex(tvv);
        };

        // `cxx:1471-1472`: `GeomAPI_ProjectPointOnCurve Projector(vtori, Einf1.Geometry())`;
        // `pvt = Projector.NearestPoint()`.
        let Some(projected) = project_point_on_curve(c1.as_ref(), &vtori, CONFUSION) else {
            panic!("StdFail_NotDone: GeomAPI_ProjectPointOnCurve::NearestPoint");
        };
        let mut pvt = projected.point;
        // `cxx:1487`
        let dion = pvt.square_distance(&vtori);

        // `cxx:1488-1520`
        self.my_v_map.change_from_index(ii).next_edge();
        if self.my_v_map.find_from_index(ii).more_edge() {
            let edg2 = self.my_v_map.find_from_index(ii).edge();
            let initiated = self.my_v_map.find_from_index(ii).parameter(&edg2);
            let Some(c2) = self.my_e_map.find_from_key(&edg2.0).and_then(|e| e.geometry()).cloned()
            else {
                return self.fail_vertex(tvv);
            };
            // `cxx:1493-1506`
            let opvt = c2.d0(initiated);
            if opvt.square_distance(&vtori) < dion {
                pvt = opvt;
            }
            // `cxx:1508-1520`
            let tol2 = BRepTool::edge_tolerance(&edg2);
            let mut done = 0i32;
            let param = parameter(&c2, &pvt, &mut done);
            let new_param = if done != 0 {
                let Some((s1, s2)) = self.edge_face_surfaces(&edg2) else {
                    return self.fail_vertex(tvv);
                };
                let einf2 = self
                    .my_e_map
                    .change_from_key(&edg2.0)
                    .expect("Standard_NoSuchObject: Draft_Modification::Perform");
                smart_parameter(einf2, tol2, &pvt, done, &s1, &s2)
            } else {
                param
            };
            *self
                .my_v_map
                .change_from_index(ii)
                .change_parameter(&edg2) = new_param;
        }

        // `cxx:1525`
        *self.my_v_map.change_from_index(ii).change_geometry() = pvt;

        // `cxx:1527-1538`
        let tol1 = BRepTool::edge_tolerance(&edg1);
        let mut done = 0i32;
        let param = parameter(&c1, &pvt, &mut done);
        let new_param = if done != 0 {
            let Some((s1, s2)) = self.edge_face_surfaces(&edg1) else {
                return self.fail_vertex(tvv);
            };
            let einf1 = self
                .my_e_map
                .change_from_key(&edg1.0)
                .expect("Standard_NoSuchObject: Draft_Modification::Perform");
            smart_parameter(einf1, tol1, &pvt, done, &s1, &s2)
        } else {
            param
        };
        *self
            .my_v_map
            .change_from_index(ii)
            .change_parameter(&edg1) = new_param;
        Ok(())
    }

    /// `cxx:1543-1637`: the point of `HAC` x `HAS` nearest the old vertex.
    fn solve_vertex(
        &mut self,
        ii: usize,
        tvv: &Vertex,
        ac: &Arc<dyn Curve>,
        as_: &Arc<dyn Surface>,
    ) -> Result<(), ()> {
        // `cxx:1545-1553`
        let intcs = perform_curve_surface(
            ac.as_ref(),
            as_.as_ref(),
            adaptor_curve_range(ac.as_ref()),
            adaptor_surface_range(as_.as_ref()),
        );
        let Ok(intcs) = intcs else {
            return self.fail_vertex(tvv);
        };

        let vtori = BRepTool::vertex_point(tvv);
        // `cxx:1558-1602`
        let pvt = if intcs.nb_points() <= 0 {
            // `cxx:1561`: `Extrema_ExtCS extr(AC, AS, PConfusion, PConfusion)`
            // over the same adaptors, i.e. their own windows.
            let (uc0, uc1) = adaptor_curve_range(ac.as_ref());
            let (u0, u1, v0, v1) = adaptor_surface_range(as_.as_ref());
            let extr = ExtCS::with_curve_surface_range(
                ac.as_ref(),
                as_.as_ref(),
                uc0,
                uc1,
                u0,
                u1,
                v0,
                v1,
                PCONFUSION,
                PCONFUSION,
            );
            if !extr.is_done() || extr.nb_ext() == 0 {
                return self.fail_vertex(tvv);
            }
            let mut disref = f64::MAX;
            let mut iref = 0usize;
            for i in 1..=extr.nb_ext() {
                let (pc, _ps) = extr.points(i);
                let distemp = pc.value().square_distance(&vtori);
                if distemp < disref {
                    disref = distemp;
                    iref = i;
                }
            }
            let (pc, _ps) = extr.points(iref);
            *pc.value()
        } else {
            // `cxx:1587-1602`
            let mut disref = f64::MAX;
            let mut iref = 0usize;
            for i in 1..=intcs.nb_points() {
                let distemp = intcs.point(i - 1).pnt.square_distance(&vtori);
                if distemp < disref {
                    disref = distemp;
                    iref = i;
                }
            }
            intcs.point(iref - 1).pnt
        };

        // `cxx:1604`
        *self.my_v_map.change_from_index(ii).change_geometry() = pvt;

        // `cxx:1606-1636`
        self.my_v_map.change_from_index(ii).init_edge_iterator();
        while self.my_v_map.find_from_index(ii).more_edge() {
            let edg = self.my_v_map.find_from_index(ii).edge();
            let initpar = self.my_v_map.find_from_index(ii).parameter(&edg);
            let Some(geom) = self.my_e_map.find_from_key(&edg.0).and_then(|e| e.geometry()).cloned()
            else {
                return self.fail_vertex(tvv);
            };
            let mut done = 0i32;
            let param = parameter(&geom, &pvt, &mut done);
            if done != 0 {
                // `cxx:1612-1621`
                let tol = BRepTool::edge_tolerance(&edg);
                let Some((s1, s2)) = self.edge_face_surfaces(&edg) else {
                    return self.fail_vertex(tvv);
                };
                let einf = self
                    .my_e_map
                    .change_from_key(&edg.0)
                    .expect("Standard_NoSuchObject: Draft_Modification::Perform");
                let new_param = smart_parameter(einf, tol, &pvt, done, &s1, &s2);
                *self
                    .my_v_map
                    .change_from_index(ii)
                    .change_parameter(&edg) = new_param;
            } else {
                // `cxx:1622-1633`
                if (initpar - param).abs() > PCONFUSION {
                    if let Some(a_c) = BRepTool::edge_curve(&edg) {
                        if a_c.is_geom_trimmed() {
                            self.my_e_map
                                .change_from_key(&edg.0)
                                .expect("Standard_NoSuchObject: Draft_Modification::Perform")
                                .set_new_geometry(true);
                        }
                    }
                }
                *self
                    .my_v_map
                    .change_from_index(ii)
                    .change_parameter(&edg) = param;
            }
            self.my_v_map.change_from_index(ii).next_edge();
        }
        Ok(())
    }

    /// `myFMap.FindFromKey(Einf.FirstFace()).Geometry()` and the second face's,
    /// the pair `SmartParameter` needs (`cxx:1510-1512`, `:1615-1617`).
    fn edge_face_surfaces(&self, e: &Edge) -> Option<(Arc<dyn Surface>, Arc<dyn Surface>)> {
        let einf = self.my_e_map.find_from_key(&e.0)?;
        let f1 = einf.first_face()?;
        let f2 = einf.second_face()?;
        let s1 = face_geometry(&self.my_f_map, f1)?;
        let s2 = face_geometry(&self.my_f_map, f2)?;
        Some((s1, s2))
    }

    /// `errStat = Draft_VertexRecomputation; badShape = TVV;` (`cxx:1541`,
    /// `:1550-1552`, `:1563-1567`).
    fn fail_vertex(&mut self, v: &Vertex) -> Result<(), ()> {
        self.err_stat = DraftErrorStatus::VertexRecomputation;
        self.bad_shape = Some(v.0.clone());
        Err(())
    }
}

/// `static bool Choose(theFMap, theEMap, Vtx, Vinf, AC, AS)`
/// (`Draft_Modification_1.cxx:2047-2177`).
///
/// Answers the `AC` / `AS` the caller loads into its adaptors; `None` is the
/// reference's `false`.
pub(super) fn choose(
    f_map: &ShapeDataMap<DraftFaceInfo>,
    e_map: &mut ShapeDataMap<DraftEdgeInfo>,
    vtx: &Vertex,
    vinf: &mut DraftVertexInfo,
) -> Option<(Arc<dyn Curve>, Arc<dyn Surface>)> {
    let vtx_pnt = BRepTool::vertex_point(vtx);

    // `cxx:2058-2082`: find a regular edge with a null `SecondFace`, or one
    // whose continuity with its two faces reaches `G1`.
    vinf.init_edge_iterator();
    while vinf.more_edge() {
        let e1 = vinf.edge();
        let info = edge_info(e_map, &e1);
        let (second_is_null, continuity) = match info.second_face() {
            None => (true, 0),
            Some(f2) => {
                let f1 = info
                    .first_face()
                    .expect("Draft_EdgeInfo: FirstFace null while SecondFace set");
                (false, stored_continuity(&e1, f1, f2))
            }
        };
        if second_is_null || continuity >= 1 {
            break;
        }
        vinf.next_edge();
    }
    if !vinf.more_edge() {
        // `cxx:2078-2081`: take the first edge.
        vinf.init_edge_iterator();
    }

    // `cxx:2084-2088`
    let eref = vinf.edge();
    let ac = edge_info(e_map, &eref).geometry().cloned()?;

    // `cxx:2090-2107`: the parameter of the vertex on the reference edge's
    // *original* curve.
    let prm = {
        let c = BRepTool::edge_curve_world(&eref)?;
        let mut done = 0i32;
        let param = parameter(&c, &vtx_pnt, &mut done);
        if done != 0 {
            let tol = BRepTool::edge_tolerance(&eref);
            let einf = e_map.change_from_key(&eref.0)?;
            let f1 = einf.first_face().cloned()?;
            let f2 = einf.second_face().cloned()?;
            let s1 = face_geometry(f_map, &f1)?;
            let s2 = face_geometry(f_map, &f2)?;
            smart_parameter(einf, tol, &vtx_pnt, done, &s1, &s2)
        } else {
            param
        }
    };

    // `cxx:2108`: `C->D1(prm, ptbid, tgref)` - the reference tangent at the
    // vertex.
    let tgref = {
        let c = BRepTool::edge_curve_world(&eref)?;
        c.d1(prm).1
    };

    // `cxx:2109-2149`: the first non-tangent cross edge.
    vinf.init_edge_iterator();
    let mut found = false;
    while vinf.more_edge() {
        let edg = vinf.edge();
        if !is_same(&edg.0, &eref.0) {
            let qualifies = {
                let einfo = edge_info(e_map, &edg);
                match einfo.second_face() {
                    None => false,
                    Some(f2) => {
                        let f1 = einfo
                            .first_face()
                            .expect("Draft_EdgeInfo: FirstFace null while SecondFace set");
                        stored_continuity(&edg, f1, f2) <= 0
                    }
                }
            };
            if qualifies {
                let c = BRepTool::edge_curve_world(&edg)?;
                let anew = {
                    let mut anewdone = 0i32;
                    let anewparam = parameter(&c, &vtx_pnt, &mut anewdone);
                    if anewdone != 0 {
                        let tol = BRepTool::edge_tolerance(&edg);
                        let einfo = e_map.change_from_key(&edg.0)?;
                        let f1 = einfo.first_face().cloned()?;
                        let f2 = einfo.second_face().cloned()?;
                        let s1 = face_geometry(f_map, &f1)?;
                        let s2 = face_geometry(f_map, &f2)?;
                        smart_parameter(einfo, tol, &vtx_pnt, anewdone, &s1, &s2)
                    } else {
                        anewparam
                    }
                };
                let tg = c.d1(anew).1;
                if tg.cross_magnitude(&tgref) > CONFUSION {
                    found = true;
                    break;
                }
            }
        }
        vinf.next_edge();
    }
    if !found {
        // `cxx:2145-2148`
        return None;
    }

    // `cxx:2150-2175`
    let e2 = vinf.edge();
    let einf_first = edge_info(e_map, &eref).first_face().cloned();
    let einf_second = edge_info(e_map, &eref).second_face().cloned();
    let einf2_first = edge_info(e_map, &e2).first_face().cloned();
    let einf2_second = edge_info(e_map, &e2).second_face().cloned();
    let same_face = |a: Option<&Face>, b: Option<&Face>| match (a, b) {
        (Some(a), Some(b)) => is_same(&a.0, &b.0),
        _ => false,
    };
    let as_face = if einf_second.is_some() {
        if same_face(einf2_first.as_ref(), einf_first.as_ref())
            || same_face(einf2_first.as_ref(), einf_second.as_ref())
        {
            einf2_second
        } else {
            einf2_first
        }
    } else if same_face(einf2_first.as_ref(), einf_first.as_ref()) {
        einf2_second
    } else {
        einf2_first
    };
    let as_surface = face_geometry(f_map, as_face.as_ref()?)?;
    Some((ac, as_surface))
}
