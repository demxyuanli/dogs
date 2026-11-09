//! `BRepFeat_MakeLinearForm::Propagate` (`BRepFeat_MakeLinearForm.cxx:1035-1288`).
//!
//! The sketch face `fac` is sectioned against the candidate faces of the basis;
//! the section edge that reaches `myFirstPnt` / `myLastPnt` selects the next
//! face through `TopExp::MapShapesAndAncestors`. The walk fills `mySlface`
//! (face -> section edges) and `myListOfEdges`, and answers `result`, which is
//! `false` as soon as one of the walked edges is not perpendicular to the rib
//! plane (`cxx:1172-1191`).

use occt_core::gp::GpPnt;
use occt_geom::extrema_surf::{ExtremaExtAlgo, ExtremaExtFlag};

use crate::abs::ShapeType;
use crate::brep_algo_api_section::SectionAlgo;
use crate::brep_extrema_ext_pc::ExtPC;
use crate::brep_extrema_ext_pf::ExtPF;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topo_tools::{explorer_shapes, is_same, AncestorMap};

use super::MakeLinearForm;

/// `RealLast()`: the largest finite value (`RealLast` is `DBL_MAX`).
const REAL_LAST: f64 = f64::MAX;

impl MakeLinearForm {
    /// `Propagate(SliList, fac, Firstpnt, Lastpnt, falseside)`
    /// (`cxx:1035-1288`).
    ///
    /// `sli_list` stands for `SliList`; the reference only reads its first
    /// element. `Err` is a raised `NCollection_DataMap::Bind`
    /// (`Standard_MultiplyDefined`, `cxx:1113` / `cxx:1278`), which the port's
    /// [`crate::locope::ShapeValueMap::bind`] reports through its return value.
    pub fn propagate(
        &mut self,
        sli_list: &[TopoShape],
        fac: &Face,
        first_pnt: &GpPnt,
        last_pnt: &GpPnt,
        falseside: &mut bool,
    ) -> Result<bool, String> {
        // `cxx:1046-1053`.
        let first_point = *first_pnt;
        let last_point = *last_pnt;
        let mut result = true;
        let mut current_face = Face::wrap(
            sli_list
                .first()
                .cloned()
                .ok_or_else(|| "Standard_NoSuchObject: SliList.First()".to_string())?,
        )
        .ok_or_else(|| "Standard_DomainError: SliList.First() is not a face".to_string())?;
        // `saveFace` (`cxx:1052`) is written by the reference and never read.

        // `cxx:1055-1057`.
        let mut last_ok = false;
        let mut first_ok = false;
        let mut v1_ok = false;
        let mut v2_ok = false;
        let mut v1: Option<Vertex> = None;
        let mut v2: Option<Vertex> = None;

        // `cxx:1059-1061`. `Approximation(true)` has no counterpart: the port
        // never builds the 3-D approximation, see
        // [`crate::brep_algo_api_section`].
        let mut sect = SectionAlgo::with_shapes(&fac.0, &current_face.0, false);
        let _ = sect.build();

        // `cxx:1063-1106`.
        let mut p1 = GpPnt::zero();
        let mut p2 = GpPnt::zero();
        let mut t1 = 0.0;
        let mut t2 = 0.0;
        let mut eb: Option<Edge> = None;
        for ec in explorer_shapes(sect.shape(), ShapeType::Edge) {
            let ec = Edge::wrap(ec)
                .ok_or_else(|| "Standard_DomainError: section edge".to_string())?;
            v1 = crate::shhealing::first_vertex(&ec);
            v2 = crate::shhealing::last_vertex(&ec);
            let (Some(vv1), Some(vv2)) = (v1.as_ref(), v2.as_ref()) else {
                continue;
            };
            p1 = BRepTool::vertex_point(vv1);
            p2 = BRepTool::vertex_point(vv2);
            t1 = BRepTool::vertex_tolerance(vv1);
            t2 = BRepTool::vertex_tolerance(vv2);
            let c1f = p1.distance(&first_point) <= t1;
            let c2f = p2.distance(&first_point) <= t2;
            let c1l = p1.distance(&last_point) <= t1;
            let c2l = p2.distance(&last_point) <= t2;
            if c1f || c2f || c1l || c2l {
                eb = Some(ec);
                if c1f || c1l {
                    v1_ok = true;
                }
                if c2f || c2l {
                    v2_ok = true;
                }
                if c1f || c2f {
                    first_ok = true;
                }
                if c1l || c2l {
                    last_ok = true;
                }
                break;
            }
        }

        // `cxx:1108-1112`.
        let Some(eb) = eb else {
            *falseside = false;
            return Ok(false);
        };

        // `cxx:1110-1116`.
        if !self.my_slface.bind(&current_face.0, Vec::new()) {
            return Err("Standard_MultiplyDefined: mySlface".to_string());
        }
        if let Some(list) = self.my_slface.change(&current_face.0) {
            list.push(eb.0.clone());
        }
        self.my_list_of_edges.clear();
        self.my_list_of_edges.push(eb.0.clone());

        // `cxx:1118-1121`.
        if last_ok && first_ok {
            return Ok(result);
        }

        // `cxx:1123-1135`.
        let mapedges = AncestorMap::map_shapes_and_ancestors(
            self.base
                .my_sbase
                .as_ref()
                .ok_or_else(|| "Standard_NullObject: mySbase".to_string())?,
            ShapeType::Edge,
            ShapeType::Face,
        );
        while !(last_ok && first_ok) {
            // `cxx:1137-1146`.
            let vprevious: Option<Vertex> = if v1_ok { v2.clone() } else { v1.clone() };
            let ptprev = if v1_ok { p2 } else { p1 };
            let vprev = vprevious
                .clone()
                .ok_or_else(|| "Standard_NullObject: Vprevious".to_string())?;

            // `cxx:1148-1196`: the edge of `CurrentFace` that reaches
            // `Vprevious`.
            let mut first_edge: Option<Edge> = None;
            for rfe in explorer_shapes(&current_face.0, ShapeType::Edge) {
                let rfe = Edge::wrap(rfe)
                    .ok_or_else(|| "Standard_DomainError: face edge".to_string())?;
                let proj = ExtPC::new(&vprev, &rfe);
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
                            // `cxx:1172-1191`: the edge must be perpendicular
                            // to the rib plane, i.e. its two vertices must
                            // project onto the same point of `fac`.
                            if result {
                                result = false;
                                let ve1 = crate::shhealing::first_vertex(&rfe);
                                let ve2 = crate::shhealing::last_vertex(&rfe);
                                if let (Some(ve1), Some(ve2)) = (ve1, ve2) {
                                    let mut perp = ExtPF::new(
                                        &ve1,
                                        fac,
                                        ExtremaExtFlag::MinMax,
                                        ExtremaExtAlgo::Grad,
                                    );
                                    if perp.is_done() {
                                        let pe1 = perp.point(1).2;
                                        perp.perform(&ve2, fac);
                                        if perp.is_done() {
                                            let pe2 = perp.point(1).2;
                                            if pe1.distance(&pe2) <= tol {
                                                result = true;
                                            }
                                        }
                                    }
                                }
                            }
                            first_edge = Some(rfe);
                            break;
                        }
                    }
                }
            }

            // `cxx:1198-1210`: the next face, the other ancestor of
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

            // `cxx:1212-1229`.
            let mut sectf = SectionAlgo::with_shapes(&fac.0, &current_face.0, false);
            let _ = sectf.build();
            let mut edg1: Option<Edge> = None;
            for e in explorer_shapes(sectf.shape(), ShapeType::Edge) {
                let e = Edge::wrap(e)
                    .ok_or_else(|| "Standard_DomainError: section edge".to_string())?;
                edg1 = Some(e.clone());
                v1 = crate::shhealing::first_vertex(&e);
                v2 = crate::shhealing::last_vertex(&e);
                let (Some(vv1), Some(vv2)) = (v1.as_ref(), v2.as_ref()) else {
                    continue;
                };
                t1 = BRepTool::vertex_tolerance(vv1);
                t2 = BRepTool::vertex_tolerance(vv2);
                p1 = BRepTool::vertex_point(vv1);
                p2 = BRepTool::vertex_point(vv2);
                v1_ok = p1.distance(&ptprev) <= t1;
                v2_ok = p2.distance(&ptprev) <= t2;
                if v1_ok || v2_ok {
                    break;
                }
            }

            // `cxx:1231-1280`: the section has to reach the sketch ends; the
            // vertex tolerance is raised to the distance when it does.
            if v1_ok {
                let vv2 = v2
                    .clone()
                    .ok_or_else(|| "Standard_NullObject: v2".to_string())?;
                if !first_ok {
                    let dp = p2.distance(&first_point);
                    if dp <= 2.0 * t2 {
                        first_ok = true;
                        vv2.set_tolerance(dp);
                    }
                }
                if !last_ok {
                    let dp = p2.distance(&last_point);
                    if dp <= 2.0 * t2 {
                        last_ok = true;
                        vv2.set_tolerance(dp);
                    }
                }
            } else if v2_ok {
                let vv1 = v1
                    .clone()
                    .ok_or_else(|| "Standard_NullObject: v1".to_string())?;
                if !first_ok {
                    let dp = p1.distance(&first_point);
                    if dp <= 2.0 * t1 {
                        first_ok = true;
                        vv1.set_tolerance(dp);
                    }
                }
                if !last_ok {
                    let dp = p1.distance(&last_point);
                    if dp <= 2.0 * t1 {
                        last_ok = true;
                        vv1.set_tolerance(dp);
                    }
                }
            } else {
                // `cxx:1277-1280`: the chain of sections ends here.
                return Ok(false);
            }

            // `cxx:1281-1285`. The reference appends `edg1` as it is, even the
            // null edge of a section that produced none; the port skips the
            // null, its lists hold `TopoShape` values.
            if !self.my_slface.bind(&current_face.0, Vec::new()) {
                return Err("Standard_MultiplyDefined: mySlface".to_string());
            }
            if let (Some(edg1), Some(list)) = (edg1, self.my_slface.change(&current_face.0)) {
                list.push(edg1.0.clone());
                self.my_list_of_edges.push(edg1.0.clone());
            }
        }

        Ok(result)
    }
}
