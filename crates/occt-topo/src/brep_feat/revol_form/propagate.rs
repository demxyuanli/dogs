//! `BRepFeat_MakeRevolutionForm::Propagate` (`cxx:1444-1933`).
//!
//! The sketch face `fac` is sectioned against the candidate faces of the basis;
//! the section edges that reach `myFirstPnt` / `myLastPnt` start two walks over
//! `TopExp::MapShapesAndAncestors(mySbase, EDGE, FACE)` that fill `mySlface`
//! (face -> section edges) and `myListOfEdges`. `false` answers mean the
//! propagation ran into itself (`cxx:1692-1696`, `cxx:1826-1830`) or off the
//! basis, and the caller turns that into `BRepFeat_FalseSide`.

use occt_core::gp::GpPnt;

use crate::abs::ShapeType;
use crate::brep_algo_api_section::SectionAlgo;
use crate::brep_extrema_ext_pc::ExtPC;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topo_tools::{explorer_shapes, is_same, AncestorMap};

use super::MakeRevolutionForm;

/// `RealLast()`: the largest finite value (`RealLast` is `DBL_MAX`).
const REAL_LAST: f64 = f64::MAX;

impl MakeRevolutionForm {
    /// `Propagate(SliList, fac, Firstpnt, Lastpnt, falseside)`
    /// (`cxx:1444-1933`).
    ///
    /// `sli_list` stands for `SliList`, which the reference mutates
    /// (`Prepend` at `cxx:1766`, `Append` at `cxx:1901`), so it is a `&mut`.
    ///
    /// `Err` is raised where the reference dereferences a null vertex through
    /// `BRep_Tool::Pnt` (`Standard_NullObject`); in that case the walk cannot
    /// advance at all.
    pub fn propagate(
        &mut self,
        sli_list: &mut Vec<TopoShape>,
        fac: &Face,
        first_pnt: &GpPnt,
        last_pnt: &GpPnt,
        falseside: &mut bool,
    ) -> Result<bool, String> {
        // `cxx:1456-1467`.
        let mut first_point = *first_pnt;
        let mut last_point = *last_pnt;
        let result = true;
        let mut current_face = Face::wrap(
            sli_list
                .first()
                .cloned()
                .ok_or_else(|| "Standard_NoSuchObject: SliList.First()".to_string())?,
        )
        .ok_or_else(|| "Standard_DomainError: SliList.First() is not a face".to_string())?;
        let save_face = current_face.clone();

        let mut last_ok = false;
        let mut first_ok = false;
        let mut v1: Option<Vertex> = None;
        let mut v2: Option<Vertex> = None;

        // `cxx:1468-1469`. `Approximation(true)` has no counterpart: the port
        // never builds the 3-D approximation, see
        // [`crate::brep_algo_api_section`].
        let mut sect = SectionAlgo::with_shapes(&fac.0, &current_face.0, false);
        let _ = sect.build();

        // `cxx:1470-1487`: the first two section edges.
        let mut e: Option<Edge> = None;
        let mut e1: Option<Edge> = None;
        let mut ii = 0;
        for ec in explorer_shapes(sect.shape(), ShapeType::Edge) {
            ii += 1;
            if ii == 1 {
                e = Edge::wrap(ec);
            } else {
                e1 = Edge::wrap(ec);
                break;
            }
        }
        let Some(mut e) = e else {
            *falseside = false;
            return Ok(false);
        };

        // `cxx:1492-1531`.
        if let Some(ee1) = e1.clone() {
            self.my_list_of_edges.clear();
            self.my_slface.bind(&current_face.0, vec![ee1.0.clone()]);
            self.my_list_of_edges.push(ee1.0.clone());

            v1 = crate::shhealing::first_vertex(&ee1);
            v2 = crate::shhealing::last_vertex(&ee1);
            let (Some(vv1), Some(vv2)) = (v1.clone(), v2.clone()) else {
                return Err("Standard_NullObject: section edge vertex".to_string());
            };
            let fp = BRepTool::vertex_point(&vv1);
            let lp = BRepTool::vertex_point(&vv2);
            let t1 = BRepTool::vertex_tolerance(&vv1);
            let t2 = BRepTool::vertex_tolerance(&vv2);
            if fp.distance(&first_point) <= t1 || fp.distance(&last_point) <= t1 {
                first_ok = true;
            }
            if lp.distance(&first_point) <= t2 || lp.distance(&last_point) <= t2 {
                last_ok = true;
            }
            if last_ok && first_ok {
                return Ok(result);
            }
            self.my_list_of_edges.clear();
        }

        // `cxx:1533-1629`: the two section edges of a split section, matched
        // against the two ends of the sketch.
        if let Some(ee1) = e1.clone() {
            self.my_list_of_edges.clear();
            self.my_slface
                .bind(&current_face.0, vec![e.0.clone(), ee1.0.clone()]);
            self.my_list_of_edges.push(e.0.clone());

            let vv1 = crate::shhealing::first_vertex(&e);
            let vv2 = crate::shhealing::last_vertex(&e);
            let v3 = crate::shhealing::first_vertex(&ee1);
            let v4 = crate::shhealing::last_vertex(&ee1);
            let (Some(vv1), Some(vv2), Some(v3), Some(v4)) = (vv1, vv2, v3, v4) else {
                return Err("Standard_NullObject: section edge vertex".to_string());
            };
            let p1 = BRepTool::vertex_point(&vv1);
            let p2 = BRepTool::vertex_point(&vv2);
            let p3 = BRepTool::vertex_point(&v3);
            let p4 = BRepTool::vertex_point(&v4);

            if p1.distance(&first_point) <= BRepTool::vertex_tolerance(&vv1) {
                if p3.distance(&last_point) <= BRepTool::vertex_tolerance(&v3) {
                    first_ok = true;
                    last_point = p4;
                } else if p4.distance(&last_point) <= BRepTool::vertex_tolerance(&v4) {
                    first_ok = true;
                    last_point = p3;
                } else {
                    e1 = None;
                }
            } else if p1.distance(&last_point) <= BRepTool::vertex_tolerance(&vv1) {
                if p3.distance(&first_point) <= BRepTool::vertex_tolerance(&v3) {
                    first_ok = true;
                    first_point = p4;
                } else if p4.distance(&first_point) <= BRepTool::vertex_tolerance(&v4) {
                    first_ok = true;
                    first_point = p3;
                } else {
                    e1 = None;
                }
            } else if p2.distance(&first_point) <= BRepTool::vertex_tolerance(&vv2) {
                if p3.distance(&last_point) <= BRepTool::vertex_tolerance(&v3) {
                    last_ok = true;
                    last_point = p4;
                } else if p4.distance(&last_point) <= BRepTool::vertex_tolerance(&v4) {
                    last_ok = true;
                    last_point = p3;
                } else {
                    e1 = None;
                }
            } else if p2.distance(&last_point) <= BRepTool::vertex_tolerance(&vv2) {
                if p3.distance(&first_point) <= BRepTool::vertex_tolerance(&v3) {
                    last_ok = true;
                    first_point = p4;
                } else if p4.distance(&first_point) <= BRepTool::vertex_tolerance(&v4) {
                    last_ok = true;
                    first_point = p3;
                } else {
                    e1 = None;
                }
            } else {
                e = ee1;
                e1 = None;
            }
            v1 = Some(vv1);
            v2 = Some(vv2);
        }

        // `cxx:1630-1660`.
        if e1.is_none() {
            self.my_list_of_edges.clear();
            self.my_slface.bind(&current_face.0, vec![e.0.clone()]);
            self.my_list_of_edges.push(e.0.clone());

            v1 = crate::shhealing::first_vertex(&e);
            v2 = crate::shhealing::last_vertex(&e);
            let (Some(vv1), Some(vv2)) = (v1.clone(), v2.clone()) else {
                return Err("Standard_NullObject: section edge vertex".to_string());
            };
            let fp = BRepTool::vertex_point(&vv1);
            let lp = BRepTool::vertex_point(&vv2);
            let t1 = BRepTool::vertex_tolerance(&vv1);
            let t2 = BRepTool::vertex_tolerance(&vv2);
            if fp.distance(&first_point) <= t1 || fp.distance(&last_point) <= t1 {
                first_ok = true;
            }
            if lp.distance(&first_point) <= t2 || lp.distance(&last_point) <= t2 {
                last_ok = true;
            }
            if last_ok && first_ok {
                return Ok(result);
            }
        }

        // `cxx:1662-1666`.
        let mapedges = AncestorMap::map_shapes_and_ancestors(
            self.base
                .my_sbase
                .as_ref()
                .ok_or_else(|| "Standard_NullObject: mySbase".to_string())?,
            ShapeType::Edge,
            ShapeType::Face,
        );
        // `FirstEdge` and `Vert` are declared once for both walks
        // (`cxx:1665`, `cxx:1456`): the second walk sees the stale value.
        let mut first_edge: Option<Edge> = None;

        let mut vprevious: Option<Vertex> = None;
        let mut vpreprevious: Option<Vertex> = None;

        // `cxx:1674-1803`: the walk from `v1`.
        while !first_ok {
            // `cxx:1677-1696`.
            let pt = v1
                .as_ref()
                .map(BRepTool::vertex_point)
                .unwrap_or_else(GpPnt::default);
            let mut looped = false;
            if let Some(vp) = vprevious.as_ref() {
                if BRepTool::vertex_point(vp).distance(&pt) <= self.my_tol {
                    looped = true;
                }
            }
            if !looped {
                if let Some(vpp) = vpreprevious.as_ref() {
                    if BRepTool::vertex_point(vpp).distance(&pt) <= self.my_tol {
                        looped = true;
                    }
                }
            }
            if looped {
                *falseside = false;
                return Ok(false);
            }

            // `cxx:1698-1726`: the edge of `CurrentFace` that reaches `v1`.
            if let Some(vv1) = v1.as_ref() {
                for rfe in explorer_shapes(&current_face.0, ShapeType::Edge) {
                    let Some(rfe) = Edge::wrap(rfe) else {
                        continue;
                    };
                    let proj = ExtPC::new(vv1, &rfe);
                    if proj.is_done() && proj.nb_ext() >= 1 {
                        let mut dist2min = REAL_LAST;
                        let mut index = 0usize;
                        for sol in 1..=proj.nb_ext() {
                            let d = proj.square_distance(sol);
                            if d <= dist2min {
                                index = sol;
                                dist2min = d;
                            }
                        }
                        if index != 0 {
                            let tol = BRepTool::edge_tolerance(&rfe);
                            if dist2min <= tol * tol {
                                first_edge = Some(rfe);
                                break;
                            }
                        }
                    }
                }
            }

            // `cxx:1728-1739`: the next face, the other ancestor of
            // `FirstEdge`. A null `FirstEdge` finds no ancestor list and keeps
            // `CurrentFace`, as `FindFromKey` on an unbound key does.
            let ancestors: Vec<TopoShape> = first_edge
                .as_ref()
                .and_then(|fe| mapedges.find_from_key(&fe.0))
                .map(|l| l.to_vec())
                .unwrap_or_default();
            for ff in ancestors {
                let Some(ff) = Face::wrap(ff) else {
                    continue;
                };
                if !is_same(&ff.0, &current_face.0) {
                    current_face = ff;
                    break;
                }
            }

            // `cxx:1741-1756`.
            let mut sectf = SectionAlgo::with_shapes(&fac.0, &current_face.0, false);
            let _ = sectf.build();
            let mut edg1: Option<Edge> = None;
            for ec in explorer_shapes(sectf.shape(), ShapeType::Edge) {
                let Some(ec) = Edge::wrap(ec) else {
                    continue;
                };
                edg1 = Some(ec.clone());
                let p1 = crate::shhealing::first_vertex(&ec)
                    .map(|v| BRepTool::vertex_point(&v));
                let p2 =
                    crate::shhealing::last_vertex(&ec).map(|v| BRepTool::vertex_point(&v));
                let (Some(p1), Some(p2), Some(vv1)) = (p1, p2, v1.as_ref()) else {
                    continue;
                };
                let t1 = BRepTool::vertex_tolerance(vv1);
                let vp = BRepTool::vertex_point(vv1);
                if p1.distance(&vp) <= t1 || p2.distance(&vp) <= t1 {
                    break;
                }
            }

            // `cxx:1758-1772`: a null `edg1` stops the walk without touching
            // `falseside`, exactly as the reference does.
            let Some(edg1) = edg1 else {
                return Ok(false);
            };
            self.my_slface.bind(&current_face.0, vec![edg1.0.clone()]);
            self.my_list_of_edges.push(edg1.0.clone());
            sli_list.insert(0, current_face.0.clone());

            // `cxx:1774-1801`.
            let Some(vert) = crate::shhealing::first_vertex(&edg1) else {
                return Err("Standard_NullObject: section edge vertex".to_string());
            };
            let Some(vv1) = v1.clone() else {
                return Err("Standard_NullObject: v1".to_string());
            };
            let pp = BRepTool::vertex_point(&vert);
            let fp = BRepTool::vertex_point(&vv1);
            let mut tol = BRepTool::edge_tolerance(&edg1);
            let tol1 = BRepTool::vertex_tolerance(&vv1);
            if tol1 > tol {
                tol = tol1;
            }
            let dist = pp.distance(&fp);
            let next_v;
            if dist <= tol {
                vpreprevious = vprevious.clone();
                vprevious = Some(vv1);
                next_v = crate::shhealing::last_vertex(&edg1);
            } else {
                vpreprevious = vprevious.clone();
                vprevious = Some(vv1);
                next_v = Some(vert);
            }
            let Some(next_v) = next_v else {
                return Err("Standard_NullObject: section edge vertex".to_string());
            };
            v1 = Some(next_v);

            let vv1 = v1.as_ref().expect("just assigned");
            let fp = BRepTool::vertex_point(vv1);
            let t = BRepTool::vertex_tolerance(vv1);
            if fp.distance(&first_point) <= t || fp.distance(&last_point) <= t {
                first_ok = true;
            }
        }

        // `cxx:1805-1808`.
        current_face = save_face;
        vprevious = None;
        vpreprevious = None;

        // `cxx:1810-1925`: the walk from `v2`.
        while !last_ok {
            let pt = v2
                .as_ref()
                .map(BRepTool::vertex_point)
                .unwrap_or_else(GpPnt::default);
            let mut looped = false;
            if let Some(vp) = vprevious.as_ref() {
                if BRepTool::vertex_point(vp).distance(&pt) <= self.my_tol {
                    looped = true;
                }
            }
            if !looped {
                if let Some(vpp) = vpreprevious.as_ref() {
                    if BRepTool::vertex_point(vpp).distance(&pt) <= self.my_tol {
                        looped = true;
                    }
                }
            }
            if looped {
                *falseside = false;
                return Ok(false);
            }

            if let Some(vv2) = v2.as_ref() {
                for rfe in explorer_shapes(&current_face.0, ShapeType::Edge) {
                    let Some(rfe) = Edge::wrap(rfe) else {
                        continue;
                    };
                    let proj = ExtPC::new(vv2, &rfe);
                    if proj.is_done() && proj.nb_ext() >= 1 {
                        let mut dist2min = REAL_LAST;
                        let mut index = 0usize;
                        for sol in 1..=proj.nb_ext() {
                            let d = proj.square_distance(sol);
                            if d <= dist2min {
                                index = sol;
                                dist2min = d;
                            }
                        }
                        if index != 0 {
                            let tol = BRepTool::edge_tolerance(&rfe);
                            if dist2min <= tol * tol {
                                first_edge = Some(rfe);
                                break;
                            }
                        }
                    }
                }
            }

            let ancestors: Vec<TopoShape> = first_edge
                .as_ref()
                .and_then(|fe| mapedges.find_from_key(&fe.0))
                .map(|l| l.to_vec())
                .unwrap_or_default();
            for ff in ancestors {
                let Some(ff) = Face::wrap(ff) else {
                    continue;
                };
                if !is_same(&ff.0, &current_face.0) {
                    current_face = ff;
                    break;
                }
            }

            let mut sectf = SectionAlgo::with_shapes(&fac.0, &current_face.0, false);
            let _ = sectf.build();
            let mut edg2: Option<Edge> = None;
            for ec in explorer_shapes(sectf.shape(), ShapeType::Edge) {
                let Some(ec) = Edge::wrap(ec) else {
                    continue;
                };
                edg2 = Some(ec.clone());
                let p1 = crate::shhealing::first_vertex(&ec)
                    .map(|v| BRepTool::vertex_point(&v));
                let p2 =
                    crate::shhealing::last_vertex(&ec).map(|v| BRepTool::vertex_point(&v));
                let (Some(p1), Some(p2), Some(vv2)) = (p1, p2, v2.as_ref()) else {
                    continue;
                };
                let t2 = BRepTool::vertex_tolerance(vv2);
                let vp = BRepTool::vertex_point(vv2);
                if p1.distance(&vp) <= t2 || p2.distance(&vp) <= t2 {
                    break;
                }
            }

            let Some(edg2) = edg2 else {
                return Ok(false);
            };
            self.my_slface.bind(&current_face.0, vec![edg2.0.clone()]);
            self.my_list_of_edges.push(edg2.0.clone());
            sli_list.push(current_face.0.clone());

            // `cxx:1906-1924`. The reference does not touch `Vprevious` on the
            // branch that keeps `v2`, so neither does the port.
            let Some(vert) = crate::shhealing::first_vertex(&edg2) else {
                return Err("Standard_NullObject: section edge vertex".to_string());
            };
            let Some(vv2) = v2.clone() else {
                return Err("Standard_NullObject: v2".to_string());
            };
            let pp = BRepTool::vertex_point(&vert);
            let fp = BRepTool::vertex_point(&vv2);
            if pp.distance(&fp) <= BRepTool::vertex_tolerance(&vv2) {
                vpreprevious = vprevious.clone();
                vprevious = Some(vv2);
                let Some(next_v) = crate::shhealing::last_vertex(&edg2) else {
                    return Err("Standard_NullObject: section edge vertex".to_string());
                };
                v2 = Some(next_v);
            } else {
                v2 = Some(vert);
            }

            let vv2 = v2.as_ref().expect("just assigned");
            let fp = BRepTool::vertex_point(vv2);
            let t = BRepTool::vertex_tolerance(vv2);
            if fp.distance(&first_point) <= t || fp.distance(&last_point) <= t {
                last_ok = true;
            }
        }

        // `cxx:1927-1931`.
        if let Some(ee1) = e1 {
            self.my_list_of_edges.push(ee1.0);
        }
        Ok(result)
    }
}
