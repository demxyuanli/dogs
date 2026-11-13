//! `ChFi3d_Builder::StartSol` — the second overload, the one driven by the
//! `SurfData` of the neighbouring band.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_2.cxx:1162-1520`
//! (OCCT 8.0.1 #1407).
//!
//! Where the first overload ([`super::start_sol`], `cxx:900-1108`) searches the
//! guideline for a first section, this one takes the last section of an
//! already computed band and walks it onto the next face:
//!
//! - both common points on arcs, with the second one beyond the extension zone,
//!   keeps the state and answers `false` (`cxx:1198-1220`),
//! - a vertex common point with an obstacle is moved onto the neighbouring edge
//!   that carries the same vertex (`cxx:1221-1352`), which is the branch that
//!   re-points `HS` / `HC` / `W` and reports `RecP` / `c1obstacle`,
//! - an arc common point with an obstacle just refreshes `W`
//!   (`cxx:1353-1360`),
//! - a detached point (`decroch`) stays on the reference face
//!   (`cxx:1374-1385`),
//! - otherwise [`ChFi3dBuilder::search_face`] looks for the face the band
//!   leaves through (`cxx:1387-1438`), and when it finds none the face it did
//!   find becomes the obstacle with its p-curve (`cxx:1439-1505`).
//!
//! Divergences that are not heuristics:
//! - `occ::handle<BRepAdaptor_Curve2d>& HC` is `&mut Option<Curve2dAdaptor>`:
//!   `HC.Nullify()` (`cxx:1362`) is `None` and `new BRepAdaptor_Curve2d`
//!   (`cxx:1457`) is `Some(Curve2dAdaptor::new())`,
//! - `occ::handle<BRepAdaptor_Surface>& HSBis` is `&mut Option<SurfaceAdaptor>`;
//!   `HSBis.Nullify()` (`cxx:1193`) is `None`,
//! - a null `TopoDS_Face` / `TopoDS_Edge` / `TopoDS_Vertex` has no counterpart
//!   in the port, so `Fref`, `Vref` and the faces are `Option`s,
//! - `aCommonPoint.Point().IsEqual(CPbis.Point(), 0)` (`cxx:1209`) is a
//!   distance against a *zero* tolerance; the port compares
//!   [`occt_core::gp::GpPnt::distance`] against `0.0` rather than using the
//!   resolution-based `is_equal`,
//! - `BRep_Tool::Parameters(V, F)` is
//!   [`BRepTool::parameters_of_vertex_on_face`], whose own docs state the
//!   unported vertex-representation loop,
//! - the four `throw Standard_Failure` of `cxx:1321-1331` are the `Err` of the
//!   port; a null p-curve OCCT dereferences unguarded
//!   (`getCurveOnSurface` results at `cxx:1373-1434`) is also an `Err`.

use std::sync::Arc;

use occt_core::gp::GpPnt2d;
use occt_geom2d::curve::Curve2d;

use crate::abs::Orientation;
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_tool::BRepTool;
use crate::chfi3d::builder_2_helpers::{is_g1, is_tangent_to_arc, search_face_on_v};
use crate::chfi_ds::{SpineHandle, SurfData};
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{self, is_same};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::ChFi3dBuilder;

/// `prepareDefaultReturn` lambda (`cxx:1391-1404`, OCCT 8.0.1 #1407): attempt to
/// set HS / W / aPCurve / pons to the default return values.
///
/// Returns `false` when the arc has no pcurve on `F` (null guard before
/// `HS->Initialize`).
#[allow(clippy::too_many_arguments)]
fn prepare_default_return(
    hs: &mut SurfaceAdaptor,
    f: &Face,
    arc: &Edge,
    param: f64,
    w: &mut f64,
    pons: &mut GpPnt2d,
    a_p_curve: &mut Option<Arc<dyn Curve2d>>,
) -> bool {
    // `cxx:1392-1397`: getCurveOnSurface first; null aborts without Initialize.
    let Some((pc, _, _)) = curve_on_surface_range(arc, f) else {
        return false;
    };
    hs.initialize(f);
    *w = param;
    *pons = pc.d0(param);
    *a_p_curve = Some(pc);
    true
}

impl ChFi3dBuilder {
    /// `StartSol(Spine, HS, pons, HC, W, SD, isfirst, ons, HSref, HCref, RecP,
    /// RecS, RecRst, c1obstacle, HSBis, PBis, decroch, Vref)`
    /// (`cxx:1162-1520`).
    ///
    /// The `bool` answer is OCCT's: `false` means "the state is preserved" (an
    /// extension, a detached point, no neighbour face), not a failure. `Err`
    /// carries the raises of `cxx:1321-1331` and the null-handle dereferences
    /// the port refuses to take.
    ///
    /// Two stores OCCT writes but never reads are kept: `newedge = EHC`
    /// (`cxx:1270`, the closed-edge branch computes `W` straight from `EHC`)
    /// and the `aPCurve` of the `isTangentToArc` sub-branch (`cxx:1422`, which
    /// `cxx:1427` overwrites). They raise the dead-assignment lint, allowed
    /// here for the same reason as in [`super::start_sol`].
    #[allow(clippy::too_many_arguments, unused_assignments)]
    pub fn start_sol_with_ref(
        &self,
        spine: Option<&SpineHandle>,
        hs: &mut SurfaceAdaptor,
        pons: &mut GpPnt2d,
        hc: &mut Option<Curve2dAdaptor>,
        w: &mut f64,
        sd: &SurfData,
        is_first: bool,
        ons: i32,
        hs_ref: &mut SurfaceAdaptor,
        hc_ref: &mut Curve2dAdaptor,
        rec_p: &mut bool,
        rec_s: &mut bool,
        rec_rst: &mut bool,
        c1_obstacle: &mut bool,
        hs_bis: &mut Option<SurfaceAdaptor>,
        p_bis: &mut GpPnt2d,
        decroch: bool,
        vref: Option<&TopoShape>,
    ) -> Result<bool, String> {
        // `RecP = RecS = RecRst = c1obstacle = false;` (`cxx:1183-1186`).
        *rec_p = false;
        *rec_s = false;
        *rec_rst = false;
        *c1_obstacle = false;

        // `TopoDS_Face Fv; occ::handle<Geom2d_Curve> aPCurve;` (`cxx:1187-1188`).
        let mut fv: Option<Face> = None;
        let mut a_p_curve: Option<Arc<dyn Curve2d>> = None;

        // `const TopOpeBRepDS_DataStructure& DStr = myDS->ChangeDS();
        // TopoDS_Face Fref = !HSref.IsNull() ? HSref->Face() : TopoDS_Face();
        // TopoDS_Face F = TopoDS::Face(DStr.Shape(SD->Index(ons)));`
        // (`cxx:1189-1191`).
        let mut f_ref = hs_ref.face().cloned();
        // `cxx:1197-1200` (OCCT 8.0.1 #1407): null face -> IsDone false path.
        let Some(f) = self
            .my_ds
            .shape(sd.index(ons), true)
            .and_then(Face::wrap)
        else {
            return Ok(false);
        };
        // `const ChFiDS_CommonPoint& aCommonPoint = SD->Vertex(isfirst, ons);
        // HSBis.Nullify();` (`cxx:1201-1202`).
        let common_point = sd.vertex(is_first, ons);
        *hs_bis = None;

        // `if (aCommonPoint.IsOnArc()) { ... }` (`cxx:1204-1230`).
        if common_point.is_on_arc() {
            let not_ons = if ons == 1 { 2 } else { 1 };
            let common_point_bis = sd.vertex(is_first, not_ons);
            if common_point_bis.is_on_arc() {
                let ts = sd.interference(ons).parameter(is_first);
                let tns = sd.interference(not_ons).parameter(is_first);
                let is_extend = if is_first {
                    ts - tns > 100.0 * self.tolesp
                } else {
                    tns - ts > 100.0 * self.tolesp
                };
                if is_extend
                    && common_point.point().distance(common_point_bis.point()) > 0.0
                {
                    // `cxx:1220-1227` (OCCT 8.0.1 #1407): null pcurve before Initialize.
                    let Some(pc) = sd.interference(ons).p_curve_on_face().cloned() else {
                        return Ok(false);
                    };
                    hs.initialize(&f);
                    *rec_s = false;
                    *pons = pc.d0(tns);
                    return Ok(false);
                }
            }
        }

        // `if (aCommonPoint.IsVertex() && !HC.IsNull() && !decroch)`
        // (`cxx:1232`).
        if common_point.is_vertex() && hc.is_some() && !decroch {
            let Some(vcp) = common_point.vertex().cloned() else {
                return Ok(false);
            };
            let Some(hc_cur) = hc.as_ref() else {
                return Ok(false);
            };
            // `const TopoDS_Edge EHC = HC->Edge();` (`cxx:1236`).
            let ehc = hc_cur.edge().clone();
            // `cxx:1238-1241` (OCCT 8.0.1 #1407): null VCP/EHC/Fref/HSref/HCref.
            if ehc.0.is_null() || f_ref.is_none() || hs_ref.face().is_none() {
                return Ok(false);
            }

            // `TopoDS_Edge newedge; TopoDS_Edge edgereg; TopoDS_Face facereg;
            // TopoDS_Face bidface = Fref; bidface.Orientation(TopAbs_FORWARD);`
            // (`cxx:1229-1233`).
            let mut new_edge: Option<Edge> = None;
            let mut edge_reg: Option<Edge> = None;
            let mut face_reg: Option<Face> = None;
            if let Some(bid_face) = f_ref
                .clone()
                .map(|face| Face(face.0.oriented(Orientation::Forward)))
            {
                // `for (TopExp_Explorer ex1(bidface, TopAbs_EDGE); ...)`
                // (`cxx:1234-1261`).
                for cured in topo_tools::edges_of(&bid_face.0) {
                    // `bool found = false; if (!cured.IsSame(EHC)) { ... }`
                    // (`cxx:1237-1255`).
                    let mut found = false;
                    if !is_same(&cured.0, &ehc.0) {
                        let (vf, vl) = topo_tools::edge_vertices(&cured);
                        let vertices = [vf, vl];
                        for vertex in vertices.iter().flatten() {
                            if !is_same(&vertex.0, &vcp) {
                                continue;
                            }
                            let mut f_voi = None;
                            if is_g1(
                                &self.my_ef_map,
                                &cured,
                                f_ref.as_ref(),
                                &mut f_voi,
                            ) {
                                // `edgereg = cured; facereg = Fv;`
                                // (`cxx:1246-1248`).
                                edge_reg = Some(cured.clone());
                                face_reg = f_voi;
                            } else {
                                // `found = true;` (`cxx:1250`).
                                found = true;
                            }
                        }
                    }
                    if found {
                        // `newedge = cured; break;` (`cxx:1257-1260`).
                        new_edge = Some(cured);
                        break;
                    }
                }
            }

            if new_edge.is_none() {
                // `TopoDS_Vertex V1, V2; TopExp::Vertices(EHC, V1, V2);`
                // (`cxx:1264-1266`).
                let (v1, v2) = topo_tools::edge_vertices(&ehc);
                let closed = match (v1.as_ref(), v2.as_ref()) {
                    (Some(a), Some(b)) => is_same(&a.0, &b.0),
                    _ => false,
                };
                if closed {
                    let (Some(v1), Some(v2)) = (v1.as_ref(), v2.as_ref()) else {
                        return Err("StartSol: closed edge without vertices".into());
                    };
                    // `newedge = EHC;` (`cxx:1270`).
                    new_edge = Some(ehc.clone());
                    // `const double w1 = BRep_Tool::Parameter(V1, EHC);
                    // const double w2 = BRep_Tool::Parameter(V2, EHC);`
                    // (`cxx:1271-1272`).
                    let w1 = BRepTool::parameter_of_vertex(&ehc.0, &v1.0);
                    let w2 = BRepTool::parameter_of_vertex(&ehc.0, &v2.0);
                    // `const ChFiDS_FaceInterference& fi = SD->Interference(ons);
                    // const handle(Geom2d_Curve)& pcf = fi.PCurveOnFace();
                    // const double ww = fi.Parameter(isfirst);` (`cxx:1273-1276`).
                    let fi = sd.interference(ons);
                    let ww = fi.parameter(is_first);
                    // `const gp_Pnt2d pww = pcf.IsNull() ?
                    // SD->Get2dPoints(isfirst, ons) : pcf->Value(ww);`
                    // (`cxx:1277`).
                    let pww = match fi.p_curve_on_face() {
                        Some(pcf) => pcf.d0(ww),
                        None => sd.get_2d_point(is_first, ons),
                    };
                    // `const gp_Pnt2d p1 = HC->Value(w1);
                    // const gp_Pnt2d p2 = HC->Value(w2);` (`cxx:1278-1279`).
                    let Some(hc_cur) = hc.as_ref() else {
                        return Err("StartSol: null obstacle".into());
                    };
                    let Some(p1) = hc_cur.d0(w1) else {
                        return Err("StartSol: null obstacle curve".into());
                    };
                    let Some(p2) = hc_cur.d0(w2) else {
                        return Err("StartSol: null obstacle curve".into());
                    };
                    // `if (p1.Distance(pww) > p2.Distance(pww)) { W = w1;
                    // pons = p1; } else { W = w2; pons = p2; }` (`cxx:1281-1288`).
                    if p1.distance(&pww) > p2.distance(&pww) {
                        *w = w1;
                        *pons = p1;
                    } else {
                        *w = w2;
                        *pons = p2;
                    }
                    // `RecP = true; c1obstacle = true; return true;`
                    // (`cxx:1289-1291`).
                    *rec_p = true;
                    *c1_obstacle = true;
                    return Ok(true);
                } else if let Some(edge_reg) = edge_reg.clone() {
                    // `Fref = facereg; HSref->Initialize(Fref);`
                    // (`cxx:1296-1297`).
                    let Some(face_reg) = face_reg.clone() else {
                        return Err("StartSol: chain is not possible, new obstacle not found".into());
                    };
                    f_ref = Some(face_reg.clone());
                    hs_ref.initialize(&face_reg);
                    // `for (TopExp_Explorer ex1(facereg, TopAbs_EDGE);
                    // ex1.More() && newedge.IsNull(); ex1.Next())`
                    // (`cxx:1298-1314`).
                    for cured in topo_tools::edges_of(&face_reg.0) {
                        if new_edge.is_some() {
                            break;
                        }
                        if is_same(&cured.0, &edge_reg.0) {
                            continue;
                        }
                        let (vf, vl) = topo_tools::edge_vertices(&cured);
                        let vertices = [vf, vl];
                        for vertex in vertices.iter().flatten() {
                            if !is_same(&vertex.0, &vcp) {
                                continue;
                            }
                            let mut f_voi = None;
                            // `if (!IsG1(myEFMap, cured, Fref, Fv)) newedge = cured;`
                            // (`cxx:1308-1311`).
                            if !is_g1(
                                &self.my_ef_map,
                                &cured,
                                f_ref.as_ref(),
                                &mut f_voi,
                            ) {
                                new_edge = Some(cured.clone());
                            }
                        }
                    }
                }
            }

            // `if (newedge.IsNull()) throw Standard_Failure("StartSol : chain
            // is not possible, new obstacle not found");` (`cxx:1320-1323`).
            let Some(new_edge) = new_edge.clone() else {
                return Err("StartSol : chain is not possible, new obstacle not found".into());
            };
            // `if (IsG1(myEFMap, newedge, Fref, Fv)) throw ... config non
            // processed` (`cxx:1324-1327`).
            let mut f_voi = None;
            if is_g1(&self.my_ef_map, &new_edge, f_ref.as_ref(), &mut f_voi) {
                return Err("StartSol : chain is not possible, config non processed".into());
            }
            fv = f_voi;
            // `else if (Fv.IsNull()) throw ...` (`cxx:1328-1332`).
            let Some(fv_face) = fv.clone() else {
                return Err("StartSol : chain is not possible, new obstacle not found".into());
            };
            // `cxx:1347-1353` / `cxx:1363-1369` (OCCT 8.0.1 #1407): null pcurves
            // on Fref / Fv abort before Initialize.
            let f_ref_face = f_ref.clone().ok_or_else(|| {
                "StartSol : chain is not possible, new obstacle not found".to_string()
            })?;
            if curve_on_surface_range(&new_edge, &f_ref_face).is_none() {
                return Ok(false);
            }
            // `TopoDS_Face newface = Fv; newface.Orientation(TopAbs_FORWARD);`
            let new_face = Face(fv_face.0.oriented(Orientation::Forward));
            let mut oriented_new_edge = new_edge.clone();
            for e in topo_tools::edges_of(&new_face.0) {
                if is_same(&e.0, &new_edge.0) {
                    oriented_new_edge = e;
                    break;
                }
            }
            if curve_on_surface_range(&oriented_new_edge, &fv_face).is_none() {
                return Ok(false);
            }
            // `HS->Initialize(Fv); W = ...; HCref->Initialize(newedge, Fref);`
            hs.initialize(&fv_face);
            *w = BRepTool::parameter_of_vertex(&new_edge.0, &vcp);
            *hc_ref = Curve2dAdaptor::of_edge_face(&new_edge, &f_ref_face);
            // `HC->Initialize(newedge, Fv); pons = HC->Value(W);`
            let adaptor = Curve2dAdaptor::of_edge_face(&oriented_new_edge, &fv_face);
            let Some(point) = adaptor.d0(*w) else {
                return Ok(false);
            };
            *pons = point;
            *hc = Some(adaptor);
            // `RecP = true; c1obstacle = true; return true;` (`cxx:1349-1351`).
            *rec_p = true;
            *c1_obstacle = true;
            return Ok(true);
        }
        // `else if (aCommonPoint.IsOnArc() && !HC.IsNull() && !decroch)`
        // (`cxx:1353`).
        else if common_point.is_on_arc() && hc.is_some() && !decroch {
            // `W = aCommonPoint.ParameterOnArc(); c1obstacle = true; return true;`
            // (`cxx:1355-1357`).
            let Some(param) = common_point.parameter_on_arc() else {
                return Err("StartSol: common point on arc without parameter".into());
            };
            *w = param;
            *c1_obstacle = true;
            return Ok(true);
        }

        // `HC.Nullify();` (`cxx:1362`).
        *hc = None;

        if common_point.is_on_arc() {
            // `const TopoDS_Edge& anArcEdge = aCommonPoint.Arc();`
            // (`cxx:1391`, OCCT 8.0.1 #1407 null guard).
            let Some(arc_shape) = common_point.arc().cloned() else {
                return Ok(false);
            };
            let Some(arc_edge) = Edge::wrap(arc_shape) else {
                return Ok(false);
            };
            if arc_edge.0.is_null() {
                return Ok(false);
            }
            let Some(param_arc) = common_point.parameter_on_arc() else {
                return Ok(false);
            };

            // `if (decroch)` (`cxx:1406-1421`, OCCT 8.0.1 #1407).
            if decroch {
                let Some(f_ref_face) = f_ref.clone() else {
                    prepare_default_return(
                        hs, &f, &arc_edge, param_arc, w, pons, &mut a_p_curve,
                    );
                    return Ok(false);
                };
                let Some((pc, _, _)) = curve_on_surface_range(&arc_edge, &f_ref_face) else {
                    prepare_default_return(
                        hs, &f, &arc_edge, param_arc, w, pons, &mut a_p_curve,
                    );
                    return Ok(false);
                };
                hs.initialize(&f_ref_face);
                *w = param_arc;
                *pons = pc.d0(*w);
                *rec_s = true;
                return Ok(true);
            }

            // `if (SearchFace(Spine, aCommonPoint, F, Fv))` (`cxx:1423`).
            if self.search_face(spine, &common_point, Some(&f), &mut fv) {
                let Some(fv_face) = fv.clone() else {
                    prepare_default_return(
                        hs, &f, &arc_edge, param_arc, w, pons, &mut a_p_curve,
                    );
                    return Ok(false);
                };
                // `HS->Initialize(Fv); RecS = true;` (`cxx:1389-1390`).
                hs.initialize(&fv_face);
                *rec_s = true;

                // `if (aCommonPoint.IsVertex())` (`cxx:1391-1405`).
                if common_point.is_vertex() {
                    let mut aux: Option<Face> = None;
                    // `const int Nb = SearchFaceOnV(aCommonPoint, F, myVEMap,
                    // myEFMap, Fv, aux);` (`cxx:1395`).
                    let nb = search_face_on_v(
                        &common_point,
                        Some(&f),
                        &self.my_ve_map,
                        &self.my_ef_map,
                        &mut fv,
                        &mut aux,
                    );
                    let Some(fv_now) = fv.clone() else {
                        return Err("StartSol: SearchFaceOnV answered without face".into());
                    };
                    let Some(vertex) = common_point.vertex() else {
                        return Err("StartSol: vertex common point without vertex".into());
                    };
                    // `pons = BRep_Tool::Parameters(aCommonPoint.Vertex(), Fv);
                    // HS->Initialize(Fv);` (`cxx:1398-1399`).
                    *pons = BRepTool::parameters_of_vertex_on_face(vertex, &fv_now)?;
                    hs.initialize(&fv_now);
                    // `if (Nb >= 2) { HSBis = new BRepAdaptor_Surface(aux);
                    // PBis = BRep_Tool::Parameters(aCommonPoint.Vertex(), aux); }`
                    // (`cxx:1400-1404`).
                    if nb >= 2 {
                        let Some(aux_face) = aux.clone() else {
                            return Err("StartSol: SearchFaceOnV answered 2 faces, one null".into());
                        };
                        *hs_bis = Some(SurfaceAdaptor::of_face(&aux_face));
                        *p_bis = BRepTool::parameters_of_vertex_on_face(vertex, &aux_face)?;
                    }
                    return Ok(true);
                }

                // `if (!Fv.IsSame(F))` (`cxx:1408-1437`).
                if !is_same(&fv_face.0, &f.0) {
                    // `Fv.Orientation(TopAbs_FORWARD);` (`cxx:1409`).
                    let fv_forward = Face(fv_face.0.oriented(Orientation::Forward));
                    // `TopoDS_Edge newedge; for (TopExp_Explorer ex(Fv,
                    // TopAbs_EDGE); ...) if (ex.Current().IsSame(anArcEdge)) ...`
                    // (`cxx:1410-1418`).
                    let mut new_edge = None;
                    for e in topo_tools::edges_of(&fv_forward.0) {
                        if is_same(&e.0, &arc_edge.0) {
                            new_edge = Some(e);
                            break;
                        }
                    }
                    // `if (isTangentToArc(aCommonPoint, 0.1)) { aPCurve =
                    // getCurveOnSurface(aCommonPoint.Arc(), F); HSBis = new
                    // BRepAdaptor_Surface(F); PBis =
                    // aPCurve->Value(aCommonPoint.ParameterOnArc()); }`
                    // (`cxx:1419-1425`).
                    if is_tangent_to_arc(&common_point, 0.1) {
                        // `cxx:1472-1478` (OCCT 8.0.1 #1407).
                        let Some((pc, _, _)) = curve_on_surface_range(&arc_edge, &f) else {
                            prepare_default_return(
                                hs, &f, &arc_edge, param_arc, w, pons, &mut a_p_curve,
                            );
                            return Ok(false);
                        };
                        *hs_bis = Some(SurfaceAdaptor::of_face(&f));
                        *p_bis = pc.d0(param_arc);
                        a_p_curve = Some(pc);
                    }
                    // `cxx:1482-1486` (OCCT 8.0.1 #1407): null newedge.
                    let Some(new_edge) = new_edge else {
                        prepare_default_return(
                            hs, &f, &arc_edge, param_arc, w, pons, &mut a_p_curve,
                        );
                        return Ok(false);
                    };
                    let Some((pc, _, _)) = curve_on_surface_range(&new_edge, &fv_forward) else {
                        prepare_default_return(
                            hs, &f, &arc_edge, param_arc, w, pons, &mut a_p_curve,
                        );
                        return Ok(false);
                    };
                    a_p_curve = Some(pc);
                } else {
                    let new_edge = Edge(arc_edge.0.oriented(Orientation::Reversed));
                    let fv_forward = Face(fv_face.0.oriented(Orientation::Forward));
                    let Some((pc, _, _)) = curve_on_surface_range(&new_edge, &fv_forward) else {
                        prepare_default_return(
                            hs, &f, &arc_edge, param_arc, w, pons, &mut a_p_curve,
                        );
                        return Ok(false);
                    };
                    a_p_curve = Some(pc);
                }
                // `cxx:1495-1500` (OCCT 8.0.1 #1407): null aPCurve.
                let Some(pc) = a_p_curve.as_ref() else {
                    prepare_default_return(
                        hs, &f, &arc_edge, param_arc, w, pons, &mut a_p_curve,
                    );
                    return Ok(false);
                };
                *pons = pc.d0(param_arc);
                return Ok(true);
            } else if let Some(fv_face) = fv.clone() {
                // `else if (!Fv.IsNull())` (`cxx:1439`).
                *c1_obstacle = true;
                // `if (!Vref.IsNull()) for (TopExp_Explorer ex(anArcEdge,
                // TopAbs_VERTEX); ...) if (ex.Current().IsSame(Vref)) {
                // c1obstacle = false; break; }` (`cxx:1441-1453`).
                if let Some(vref) = vref {
                    let (vf, vl) = topo_tools::edge_vertices(&arc_edge);
                    let vertices = [vf, vl];
                    for vertex in vertices.iter().flatten() {
                        if is_same(&vertex.0, vref) {
                            *c1_obstacle = false;
                            break;
                        }
                    }
                }
                if *c1_obstacle {
                    // `cxx:1519-1524` (OCCT 8.0.1 #1407): null HSref/HCref.
                    if hs_ref.face().is_none() {
                        prepare_default_return(
                            hs, &f, &arc_edge, param_arc, w, pons, &mut a_p_curve,
                        );
                        return Ok(false);
                    }
                    hs.initialize(&fv_face);
                    hs_ref.initialize(&f);
                    *w = param_arc;
                    *hc = Some(Curve2dAdaptor::new());
                    let new_face = Face(fv_face.0.oriented(Orientation::Forward));
                    let mut new_edge = None;
                    for e in topo_tools::edges_of(&new_face.0) {
                        if is_same(&e.0, &arc_edge.0) {
                            new_edge = Some(e);
                            break;
                        }
                    }
                    match new_edge {
                        Some(new_edge) => {
                            let adaptor = Curve2dAdaptor::of_edge_face(&new_edge, &fv_face);
                            let Some(point) = adaptor.d0(*w) else {
                                prepare_default_return(
                                    hs, &f, &arc_edge, param_arc, w, pons, &mut a_p_curve,
                                );
                                return Ok(false);
                            };
                            *pons = point;
                            *hc = Some(adaptor);
                            *hc_ref = Curve2dAdaptor::of_edge_face(&arc_edge, &f);
                            // `if (aCommonPoint.IsVertex()) RecP = true;
                            // else RecRst = true; return true;`
                            // (`cxx:1475-1484`).
                            if common_point.is_vertex() {
                                *rec_p = true;
                            } else {
                                *rec_rst = true;
                            }
                            return Ok(true);
                        }
                        None => {
                            // `prepareDefaultReturn(); return false;`
                            // (`cxx:1487-1491`).
                            prepare_default_return(
                                hs, &f, &arc_edge, param_arc, w, pons, &mut a_p_curve,
                            );
                            return Ok(false);
                        }
                    }
                } else {
                    // `prepareDefaultReturn(); return false;` (`cxx:1493-1498`).
                    prepare_default_return(
                        hs, &f, &arc_edge, param_arc, w, pons, &mut a_p_curve,
                    );
                    return Ok(false);
                }
            } else {
                // `prepareDefaultReturn(); return false;` (`cxx:1500-1505`).
                prepare_default_return(hs, &f, &arc_edge, param_arc, w, pons, &mut a_p_curve);
                return Ok(false);
            }
        } else {
            // `HS->Initialize(F); const ChFiDS_FaceInterference& FI =
            // SD->Interference(ons); if (FI.PCurveOnFace().IsNull()) pons =
            // SD->Get2dPoints(isfirst, ons); else pons =
            // FI.PCurveOnFace()->Value(FI.Parameter(isfirst));`
            // (`cxx:1506-1518`).
            hs.initialize(&f);
            let fi = sd.interference(ons);
            *pons = match fi.p_curve_on_face() {
                Some(pc) => pc.d0(fi.parameter(is_first)),
                None => sd.get_2d_point(is_first, ons),
            };
        }
        // `return true;` (`cxx:1519`).
        Ok(true)
    }
}
