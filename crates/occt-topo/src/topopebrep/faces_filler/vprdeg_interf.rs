//! `TopOpeBRep_vprdeg.cxx:381-862` — the degenerated-edge interferences of a
//! VPoint.
//!
//! Sources (OCCT 8.0.0, `ModelingAlgorithms/TKBool/TopOpeBRep`):
//! - `FUN_putInterfonDegenEd(VP, F1, F2, DataforDegenEd, HDS, is, dgE,
//!   makeinterf, Trans1, param1, Trans2, param2, OOEi, paronOOEi, hasOOEi,
//!   isT2d)` (`vprdeg.cxx:381-802`)
//! - `ProcessVPondgE(VP, ShapeIndex, PVKind, PVIndex, EPIfound, IEPI,
//!   CPIfound, ICPI)` (`cxx:633-801`)
//! - `local_FindTreatedEdgeOnVertex` (`cxx:805-820`),
//!   `local_FindVertex` (`cxx:822-841`),
//!   `local_ReduceMapOfTreatedVertices` (`cxx:843-860`) and the two file-scope
//!   statics they drive (`cxx:48-51`): `aMapOfTreatedVertexListOfEdge` and
//!   `localCurrentLine`

use std::cell::RefCell;

use crate::abs::{Orientation, ShapeType};
use crate::shape::TopoShape;
use crate::topo_tools::{explorer_shapes, is_same, map_shapes_and_ancestors, shape_map_key};
use crate::topopebrep::ds::export::fun_ds_getoov_in_hds;
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::ds::kind::Kind;
use crate::topopebrep::ds::process_interferences_tool::make_epv_interference_with_support;
use crate::topopebrep::ds::transition::Transition;
use crate::topopebrep::faces_filler::DegenEdgeData;
use crate::topopebrep::faces_filler::vp_index::{fun_vp_index, VpIndex, M_FINDVP, M_MKNEWVP};
use crate::topopebrep::faces_filler::FacesFiller;
use crate::topopebrep::faces_intersector::FacesIntersector;
use crate::topopebrep::interf::v_point_inter::VPointInter;
use crate::topopebrep::tool::mk_tondg_e::{MktondgE, MKI1, MKI12, MKI2, NOI};
use crate::topopebrep::tool::topology::fun_tool_in_s;

thread_local! {
    /// `aMapOfTreatedVertexListOfEdge` (`vprdeg.cxx:48-50`).
    static TREATED_EDGES: RefCell<Vec<(TopoShape, Vec<TopoShape>)>> =
        const { RefCell::new(Vec::new()) };
    /// `localCurrentLine` (`vprdeg.cxx:51`), carried as the current line index.
    static CURRENT_LINE: RefCell<Option<i32>> = const { RefCell::new(None) };
}

/// `local_FindTreatedEdgeOnVertex(theEdge, theVertex)` (`cxx:805-820`).
fn local_find_treated_edge_on_vertex(edge: &TopoShape, vertex: &TopoShape) -> bool {
    TREATED_EDGES.with(|m| {
        m.borrow()
            .iter()
            .find(|(k, _)| is_same(k, vertex))
            .map(|(_, l)| l.iter().any(|x| is_same(x, edge)))
            .unwrap_or(false)
    })
}

/// The `Bind` half of `aMapOfTreatedVertexListOfEdge(aVertex).Append(ei)`
/// (`cxx:582-588`).
fn local_append_treated_edge(vertex: &TopoShape, edge: &TopoShape) {
    TREATED_EDGES.with(|m| {
        let mut m = m.borrow_mut();
        match m.iter_mut().find(|(k, _)| is_same(k, vertex)) {
            Some((_, l)) => l.push(edge.clone()),
            None => m.push((vertex.clone(), vec![edge.clone()])),
        }
    });
}

/// `local_FindVertex(theVP, theMapOfVertexEdges, theVertex)` (`cxx:822-841`).
///
/// `theMapOfVertexEdges` is an `NCollection_IndexedDataMap`; `vertex_order` is
/// its key order, which for `TopExp::MapShapesAndAncestors(F, VERTEX, EDGE)` is
/// the face's edge-walk order (the port's `HashMap` does not preserve it).
fn local_find_vertex(
    vp: &VPointInter,
    vertex_order: &[TopoShape],
) -> Result<Option<TopoShape>, String> {
    let a_vp_point = vp.value()?.clone();
    let a_vp_tol = vp.tolerance()?;
    for v in vertex_order {
        let p = crate::brep_tool::BRepTool::vertex_point(&crate::shape::Vertex(v.clone()));
        if a_vp_point.is_equal_with_tolerance(&p, a_vp_tol) {
            return Ok(Some(v.clone()));
        }
    }
    Ok(None)
}

/// `local_ReduceMapOfTreatedVertices(theCurrentLine)` (`cxx:843-860`): clears
/// the treated-edge map when the line changes.
fn local_reduce_map_of_treated_vertices(current_line: i32) {
    CURRENT_LINE.with(|c| {
        let mut c = c.borrow_mut();
        if *c != Some(current_line) {
            *c = Some(current_line);
            TREATED_EDGES.with(|m| m.borrow_mut().clear());
        }
    });
}

/// The out-parameters of `FUN_putInterfonDegenEd` (`vprdeg.cxx:386-399`).
pub struct DegenEdInterf {
    /// `is`: the rank of the degenerated edge.
    pub is_rank: i32,
    /// `dgE`: the degenerated edge.
    pub dg_e: Option<TopoShape>,
    /// `Trans1` / `param1`.
    pub trans1: Transition,
    pub param1: f64,
    /// `Trans2` / `param2`.
    pub trans2: Transition,
    pub param2: f64,
    /// `OOEi` / `paronOOEi` / `hasOOEi`: the restriction edge the transitions
    /// were computed against, when one was found.
    pub oo_ei: Option<TopoShape>,
    pub par_on_oo_ei: f64,
    pub has_oo_ei: bool,
    /// `isT2d`: the two faces are tangent at the VP.
    pub is_t2d: bool,
    /// The `FUN_putInterfonDegenEd` return value: `NOI` / `MKI1` / `MKI2` /
    /// `MKI12`.
    pub make_i: i32,
}

impl Default for DegenEdInterf {
    fn default() -> Self {
        Self {
            is_rank: 0,
            dg_e: None,
            trans1: Transition::default(),
            param1: 0.0,
            trans2: Transition::default(),
            param2: 0.0,
            oo_ei: None,
            par_on_oo_ei: 0.0,
            has_oo_ei: false,
            is_t2d: false,
            make_i: NOI,
        }
    }
}

/// `FUN_putInterfonDegenEd(...)` (`vprdeg.cxx:381-802`).
pub fn fun_put_interfon_degen_ed(
    vp: &VPointInter,
    f1: &TopoShape,
    f2: &TopoShape,
    data_for_degen_ed: &DegenEdgeData,
    hds: &mut HDataStructure,
) -> Result<DegenEdInterf, String> {
    let mut out = DegenEdInterf::default();

    let on3 = vp.shape_index() == 3;
    let onv12 = vp.is_vertex_on_s1()? && vp.is_vertex_on_s2()?;

    // `cxx:407-429`: the rank of `VP` that carries degenerated data.
    let mut rkv = 0;
    let mut v: Option<TopoShape> = None;
    for shape_index in 1..=2 {
        let isv = if shape_index == 1 {
            vp.is_vertex_on_s1()?
        } else {
            vp.is_vertex_on_s2()?
        };
        if !isv {
            continue;
        }
        let vv = if shape_index == 1 {
            vp.vertex_on_s1()?.clone()
        } else {
            vp.vertex_on_s2()?.clone()
        };
        if data_for_degen_ed.find(&vv).is_none() {
            continue;
        }
        rkv = shape_index;
        v = Some(vv);
        break;
    }
    if rkv == 0 {
        return Ok(out); // compute interference once only.
    }
    let Some(mut v) = v else {
        return Ok(out);
    };
    let isvsd = hds.has_same_domain(&v, true);

    // `cxx:434-436`: edges `dge`, `cle` on shape<rkdg>.
    let loe: Vec<TopoShape> = match data_for_degen_ed.find(&v) {
        Some(l) => l.to_vec(),
        None => return Ok(out),
    };
    let (Some(cle), Some(dge)) = (loe.first().cloned(), loe.last().cloned()) else {
        return Ok(out);
    };
    out.dg_e = Some(dge.clone());

    // `cxx:438-457`
    let rkdg = if hds.ds().has_shape(&dge, false) {
        hds.ds().ancestor_rank_of(&dge)
    } else if fun_tool_in_s(&v, &dge) {
        rkv
    } else if rkv == 1 {
        2
    } else {
        1
    };
    out.is_rank = rkdg;
    let rki = if rkdg == 1 { 2 } else { 1 };

    // `cxx:459-477`
    let uvi = if rki == 1 {
        let (u, vv) = vp.parameters_on_s1()?;
        occt_core::gp::GpPnt2d::new(u, vv)
    } else {
        let (u, vv) = vp.parameters_on_s2()?;
        occt_core::gp::GpPnt2d::new(u, vv)
    };
    let fi = if rki == 1 { f1 } else { f2 };
    let f = if rkdg == 1 { f1 } else { f2 };

    // `cxx:480-490`
    let mut mktdg = MktondgE::new();
    if !mktdg.initialize(&dge, f, &uvi, fi) {
        return Ok(out);
    }
    if !mktdg.set_cl_e(&cle) {
        return Ok(out);
    }

    // `cxx:492-512`: on both edges, or same-domain, `v` becomes the vertex of
    // `dge` itself.
    if onv12 || isvsd {
        let ov = if onv12 {
            if rkv == 2 {
                vp.vertex_on_s1()?.clone()
            } else {
                vp.vertex_on_s2()?.clone()
            }
        } else {
            match fun_ds_getoov_in_hds(hds, &v) {
                Some(o) => o,
                None => return Ok(out),
            }
        };
        if rkv != rkdg {
            // `cxx:511`: ensures `v` is a vertex of `dge`; only `v` is read
            // afterwards, so OCCT's `rkv = rkdg` there is a dead store.
            v = ov;
        }
    }
    let _ = &v;

    // `cxx:514-601`
    let mut mkt = NOI;
    let mut par1 = 0.0;
    let mut par2 = 0.0;
    if on3 {
        let ei = if rki == 1 {
            vp.arc_on_s1()?.clone()
        } else {
            vp.arc_on_s2()?.clone()
        };
        let pari = if rki == 1 {
            vp.parameter_on_arc1()?
        } else {
            vp.parameter_on_arc2()?
        };
        mktdg.set_rest(pari, &ei);
        mkt = mktdg.mk_ton_e(&mut par1, &mut par2);
        if mkt == NOI {
            return Ok(out);
        }
        out.oo_ei = Some(ei);
        out.par_on_oo_ei = pari;
        out.has_oo_ei = true;
    } else {
        let on_dom = if rki == 1 {
            vp.is_on_dom_s1()?
        } else {
            vp.is_on_dom_s2()?
        };
        if on_dom {
            // `cxx:533-548`
            let ei = if rki == 1 {
                vp.arc_on_s1()?.clone()
            } else {
                vp.arc_on_s2()?.clone()
            };
            let pari = if rki == 1 {
                vp.parameter_on_arc1()?
            } else {
                vp.parameter_on_arc2()?
            };
            mktdg.set_rest(pari, &ei);
            mkt = mktdg.mk_ton_e(&mut par1, &mut par2);
            if mkt != NOI {
                out.oo_ei = Some(ei);
                out.par_on_oo_ei = pari;
                out.has_oo_ei = true;
            }
        } else {
            // `cxx:550-592`: an untreated restriction of the other face that
            // carries the VP's vertex.
            let mut edge_found = false;
            let a_face = if rki == 1 { f1 } else { f2 };
            let a_map = map_shapes_and_ancestors(a_face, ShapeType::Vertex, ShapeType::Edge);
            let vertex_order = dedup_occurrences(explorer_shapes(a_face, ShapeType::Vertex));
            let a_vertex = local_find_vertex(vp, &vertex_order)?;
            if let Some(a_vertex) = a_vertex {
                if !a_vertex.is_null() {
                    let edges = a_map
                        .get(&shape_map_key(&a_vertex))
                        .cloned()
                        .unwrap_or_default();
                    for ei in edges {
                        if edge_found {
                            break;
                        }
                        let pari =
                            crate::brep_tool::BRepTool::parameter_of_vertex(&ei, &a_vertex);
                        if !crate::brep_tool::BRepTool::is_degenerated(&crate::shape::Edge(
                            ei.clone(),
                        )) {
                            edge_found = !local_find_treated_edge_on_vertex(&ei, &a_vertex);
                        }
                        if edge_found {
                            mktdg.set_rest(pari, &ei);
                            mkt = mktdg.mk_ton_e(&mut par1, &mut par2);
                            if mkt != NOI {
                                out.oo_ei = Some(ei.clone());
                                out.par_on_oo_ei = pari;
                                out.has_oo_ei = true;
                            }
                            local_append_treated_edge(&a_vertex, &ei);
                        }
                    }
                }
            }
            if !edge_found {
                mkt = mktdg.mk_ton_e(&mut par1, &mut par2);
            }
        }
        if mkt == NOI {
            return Ok(out);
        }
    }
    out.is_t2d = mktdg.is_t2d();

    // `cxx:604-613`
    if mkt == MKI1 || mkt == MKI12 {
        out.trans1.set_orientation(Orientation::Forward);
        out.param1 = par1;
    }
    if mkt == MKI2 || mkt == MKI12 {
        out.trans2.set_orientation(Orientation::Reversed);
        out.param2 = par2;
    }
    out.make_i = mkt;
    Ok(out)
}

/// `TopExp::MapShapes(F, VERTEX)` key order without repeated occurrences: the
/// `NCollection_IndexedDataMap` `FindKey(iVertex)` walk of `cxx:830-839`.
fn dedup_occurrences(shapes: Vec<TopoShape>) -> Vec<TopoShape> {
    let mut out: Vec<TopoShape> = Vec::new();
    for s in shapes {
        if !out.iter().any(|o| is_same(o, &s)) {
            out.push(s);
        }
    }
    out
}

impl FacesFiller {
    /// `ProcessVPondgE(VP, ShapeIndex, PVKind, PVIndex, EPIfound, IEPI,
    /// CPIfound, ICPI)` (`vprdeg.cxx:633-801`).
    ///
    /// `PVKind` / `PVIndex` are in/out parameters in OCCT: `pv_kind` /
    /// `pv_index` are their values on entry, and the lookup at `cxx:642-658`
    /// only runs when `PVIndex` is zero (`ProcessRLine` carries a non-zero
    /// value from its first `ShapeIndex` iteration into the second,
    /// `_1.cxx:880`).
    ///
    /// The returned `bool` is the OCCT return value, which `ProcessRLine`
    /// does not read (`_1.cxx:917`); the out-values are in the [`VpIndex`].
    #[allow(clippy::too_many_arguments)]
    pub fn process_vp_on_dg_e(
        &mut self,
        facint: &mut FacesIntersector,
        hds: &mut HDataStructure,
        vp: &VPointInter,
        shape_index: i32,
        pv_kind: Kind,
        pv_index: i32,
    ) -> Result<(bool, VpIndex), String> {
        let mut idx = VpIndex::default();
        idx.pv_kind = pv_kind;
        idx.pv_index = pv_index;
        // `cxx:642-658`
        if idx.pv_index == 0 {
            idx = fun_vp_index(
                self,
                facint,
                self.line_index,
                vp,
                shape_index,
                hds,
                &self.dscil,
                M_FINDVP,
            )?;
        }

        // `cxx:673-680`
        if self.data_for_degen_ed.is_empty() {
            return Ok((false, idx));
        }
        local_reduce_map_of_treated_vertices(self.line_index);

        let f1 = self
            .f1
            .as_ref()
            .ok_or_else(|| "TopOpeBRep_FacesFiller::ProcessVPondgE : no F1".to_string())?
            .0
            .clone();
        let f2 = self
            .f2
            .as_ref()
            .ok_or_else(|| "TopOpeBRep_FacesFiller::ProcessVPondgE : no F2".to_string())?
            .0
            .clone();
        let dg = fun_put_interfon_degen_ed(vp, &f1, &f2, &self.data_for_degen_ed, hds)?;
        if dg.make_i == NOI {
            return Ok((false, idx));
        }

        // `cxx:708-723`
        if idx.pv_index == 0 {
            idx = fun_vp_index(
                self,
                facint,
                self.line_index,
                vp,
                shape_index,
                hds,
                &self.dscil,
                M_MKNEWVP,
            )?;
        }

        // `cxx:731-747`
        let rank_fi = if dg.is_rank == 1 { 2 } else { 1 };
        let fi = if rank_fi == 1 { f1 } else { f2 };
        let i_fi = hds.change_ds().add_shape_with_ancestor(&fi, rank_fi);
        let dg_e = match dg.dg_e.clone() {
            Some(e) => e,
            None => return Ok((false, idx)),
        };
        hds.change_ds().add_shape_with_ancestor(&dg_e, dg.is_rank);

        let rkv = match hds.ds().shape(idx.pv_index, false) {
            Some(s) => hds.ds().ancestor_rank_of(&s),
            None => 0,
        };

        // `cxx:748-798`
        if dg.make_i == MKI1 || dg.make_i == MKI12 {
            let mut t1 = dg.trans1.clone();
            t1.set_index(i_fi);
            let is_vertex1 = rkv == 1;
            if !dg.is_t2d {
                if let Some(evi1) = make_epv_interference_with_support(
                    &t1,
                    i_fi,
                    idx.pv_index,
                    dg.param1,
                    Kind::Vertex,
                    Kind::Face,
                    is_vertex1,
                ) {
                    hds.store_interference_on_shape(&evi1, &dg_e)?;
                }
            }
        }
        if dg.make_i == MKI2 || dg.make_i == MKI12 {
            let mut t2 = dg.trans2.clone();
            t2.set_index(i_fi);
            let is_vertex2 = rkv == 2;
            if !dg.is_t2d {
                if let Some(evi2) = make_epv_interference_with_support(
                    &t2,
                    i_fi,
                    idx.pv_index,
                    dg.param2,
                    Kind::Vertex,
                    Kind::Face,
                    is_vertex2,
                ) {
                    hds.store_interference_on_shape(&evi2, &dg_e)?;
                }
            }
        }
        Ok((true, idx))
    }
}
