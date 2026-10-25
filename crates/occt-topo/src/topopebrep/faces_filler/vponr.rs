//! `TopOpeBRep_vpr.cxx:877-1732` — `ProcessVPonR`, the VPoint-on-restriction
//! pass that builds the curve/point and edge/point interferences of a line.
//!
//! Sources (OCCT 8.0.0, `ModelingAlgorithms/TKBool/TopOpeBRep`):
//! - `ProcessVPonR` (`TopOpeBRep_vpr.cxx:877-1732`), split here at
//!   `cxx:1238` ("DS geometry Management"): the transition decisions are in
//!   this file, the geometry/interference part is in [`super::vponr_ds`]
//! - `FUN_onedge` (`vpr.cxx:848-866`)
//! - `FUN_brep_ONfirstP` (`TopOpeBRep_FacesFiller_1.cxx:634-644`)
//!
//! Ported in the sibling modules and called from here:
//! - `ProcessVPondgE` (`TopOpeBRep_vprdeg.cxx:633-801`) with
//!   `FUN_putInterfonDegenEd` (`vprdeg.cxx:381-802`) in
//!   [`super::vprdeg_interf`], and `FUN_GetdgData` / `FUN_FillVof12`
//!   (`vprdeg.cxx:105-379`) in [`super::degen_data`].
//! - `FUN_VPIndex` (`vpr.cxx:350-469`) in [`super::vp_index`].
//! - `FUN_newtransEdge` / `FUN_LineRestF` (`vpr.cxx:471-682`) in
//!   [`super::vpr_trans_edge`].
//! - `FUN_processCPI` (`vpr.cxx:782-853`) in [`super::vpr_cpi`], and the
//!   `FUNBREP_topokpart` (`TopOpeBRep_kpart.cxx:831-916`) it calls in
//!   [`crate::topopebrep::kpart`].
//! - `GetEdgeTrans` / `ProcessVPonclosingR` (`vprclo.cxx:66-296`) in
//!   [`super::vpclo`].
//!
//! Faithfulness notes:
//! - `transEdge.Before(Tr.Before())` (`cxx:1193-1194`) calls the OCCT overload
//!   whose `ShapeBefore` argument defaults to `TopAbs_FACE`
//!   (`TopOpeBRepDS_Transition.hxx:55-56`): the *shape* of `Tr` is **not**
//!   copied, only its state. The port spells the default out.
//! - `tredunk = transEdge.IsUnknown()` (`cxx:1239`) is dead in OCCT: the flag
//!   is not read before the function ends. The port does not repeat it.

use crate::abs::{Orientation, ShapeType, State};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, TopoShape, Vertex};
use crate::topopebrep::ds::export::{fds_has_unk, fds_parameter};
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::ds::interference::InterferenceHandle;
use crate::topopebrep::ds::kind::Kind;
use crate::topopebrep::ds::point::Point;
use crate::topopebrep::ds::transition::Transition;
use crate::topopebrep::faces_filler::vpr_trans_edge::{
    fun_newtrans_edge, fun_scan_interf_list, fun_select_g_interference,
    fun_select_traisha_interference,
};
use crate::topopebrep::faces_intersector::FacesIntersector;
use crate::topopebrep::ff_transition_tool::FFTransitionTool;
use crate::topopebrep::interf::point_geom_tool::PointGeomTool;
use crate::topopebrep::interf::v_point_inter::VPointInter;
use crate::topopebrep::type_line_curve::TypeLineCurve;

/// `FUN_onedge(PDS, E)` (`vpr.cxx:848-866`): the point of `PDS` coincides with
/// one of `E`'s two vertices within `max(tolP, tolE)`.
pub fn fun_on_edge(pds: &Point, e: &TopoShape) -> bool {
    let p = pds.point();
    let tol_p = pds.tolerance();
    let tol_e = BRepTool::edge_tolerance(&Edge(e.clone()));
    let tol = tol_p.max(tol_e);
    let (vf, vl) = crate::topopebrep::tool::closed::vertices(e);
    let (Some(vf), Some(vl)) = (vf, vl) else {
        return false;
    };
    let pf = BRepTool::vertex_point(&Vertex(vf));
    let pl = BRepTool::vertex_point(&Vertex(vl));
    p.is_equal_with_tolerance(&pf, tol) || p.is_equal_with_tolerance(&pl, tol)
}

/// `FUN_brep_ONfirstP(vpf, VP)` (`TopOpeBRep_FacesFiller_1.cxx:634-644`): the
/// parameter of `VP` on the line differs from the first VP's by less than
/// `Precision::Confusion()`.
pub fn fun_brep_on_first_p(vpf: &VPointInter, vp: &VPointInter) -> Result<bool, String> {
    let parfirst = vpf.parameter_on_line()?;
    let parcur = vp.parameter_on_line()?;
    Ok((parcur - parfirst).abs() < occt_core::precision::CONFUSION)
}

/// The state [`super::vponr_ds`] takes over from `ProcessVPonR` `cxx:877-1238`.
pub(super) struct VpOnRState {
    /// `ShapeIndex` / `OOShapeIndex`.
    pub shape_index: i32,
    pub oo_shape_index: i32,
    /// `Face(ShapeIndex)` / `Face(OOShapeIndex)`.
    pub face: TopoShape,
    pub oo_face: TopoShape,
    /// `VP.ShapeIndex()` (`cxx:884`).
    pub absindex: i32,
    pub on2edges: bool,
    pub has_on_edge: bool,
    pub has_oo_edge: bool,
    /// `iSIFace` / `iOOFace`.
    pub i_si_face: i32,
    pub i_oo_face: i32,
    /// `SIedgeIndex` / `edge` / `paredge`.
    pub si_edge_index: i32,
    pub edge: TopoShape,
    pub paredge: f64,
    pub isrest: bool,
    pub closing: bool,
    pub dge: bool,
    /// `OOedgeIndex` / `OOedge` / `OOparedge`.
    pub oo_edge_index: i32,
    pub oo_edge: TopoShape,
    pub oo_paredge: f64,
    pub oo_isrest: bool,
    pub oo_closing: bool,
    pub dg_oo_e: bool,
    /// `PVKind` / `PVIndex` / `EPIfound` / `IEPI` / `CPIfound` / `ICPI`.
    pub pv_kind: Kind,
    pub pv_index: i32,
    pub epi_found: bool,
    pub iepi: Option<InterferenceHandle>,
    pub cpi_found: bool,
    pub icpi: Option<InterferenceHandle>,
    /// `foundPVIndex = (PVIndex != 0)` (`cxx:975`).
    pub found_pv_index: bool,
    /// `PDS` (`cxx:1147`).
    pub pds: Point,
    /// `transLine` / `transEdge`.
    pub trans_line: Transition,
    pub trans_edge: Transition,
    /// `keepvp` (`cxx:1027`).
    pub keepvp: bool,
    /// `found` (`cxx:1150`).
    pub found: bool,
}

impl super::FacesFiller {
    /// `ProcessVPonR(VP, Trans, GFace, ShapeIndex)` (`vpr.cxx:877-1732`).
    #[allow(clippy::too_many_arguments)]
    pub fn process_vp_on_r(
        &mut self,
        facint: &mut FacesIntersector,
        hds: &mut HDataStructure,
        vp: &VPointInter,
        trans: &Transition,
        face: &TopoShape,
        shape_index: i32,
    ) -> Result<(), String> {
        let Some(mut st) = self.vp_on_r_begin(facint, hds, vp, face, shape_index)? else {
            return Ok(());
        };
        if !self.vp_on_r_trans_line(facint, vp, &mut st)? {
            return Ok(());
        }
        if !self.vp_on_r_trans_edge(facint, hds, vp, trans, &mut st)? {
            return Ok(());
        }
        if !self.vp_on_r_geometry(facint, hds, vp, &mut st)? {
            return Ok(());
        }
        self.vp_on_r_store(facint, hds, vp, &mut st)
    }

    /// `vpr.cxx:883-975`: the indices, the closed/degenerated flags and the
    /// `ProcessVPondgE` prologue.
    fn vp_on_r_begin(
        &mut self,
        facint: &mut FacesIntersector,
        hds: &mut HDataStructure,
        vp: &VPointInter,
        face: &TopoShape,
        shape_index: i32,
    ) -> Result<Option<VpOnRState>, String> {
        let absindex = vp.shape_index(); // 0,1,2,3
        let oo_shape_index = if shape_index == 1 { 2 } else { 1 };
        let on2edges = absindex == 3;
        let has_on_edge = vp.state(oo_shape_index)? == State::On;
        let has_oo_edge = if on2edges { true } else { has_on_edge };

        let mut i_si_face = hds.ds().shape_index(face, true);
        if i_si_face == 0 {
            i_si_face = hds
                .change_ds()
                .add_shape_with_ancestor(face, shape_index);
        }
        let oo_face = self.face(oo_shape_index)?.0.clone();
        let mut i_oo_face = hds.ds().shape_index(&oo_face, true);
        if i_oo_face == 0 {
            i_oo_face = hds
                .change_ds()
                .add_shape_with_ancestor(&oo_face, oo_shape_index);
        }

        // current VPoint is on `<edge>`
        let edge = vp
            .edge(shape_index)
            .ok_or_else(|| "TopOpeBRep_FacesFiller::ProcessVPonR : no edge".to_string())?
            .clone();
        let mut si_edge_index = 0;
        if hds.ds().has_shape(&edge, true) {
            si_edge_index = hds.ds().shape_index(&edge, true);
        }
        let paredge = vp.edge_parameter(shape_index);
        let isrest = hds.ds().is_section_edge(&edge, true);
        let closing = crate::topopebrep::tool::shape_tool::closed(&edge, face);
        let dge = BRepTool::is_degenerated(&Edge(edge.clone()));

        // dummy if !<hasOOedge>
        let mut oo_edge_index = 0;
        let mut oo_closing = false;
        let mut oo_isrest = false;
        let mut oo_edge = TopoShape::new(crate::abs::ShapeType::Edge);
        let mut oo_paredge = 0.0;
        let mut dg_oo_e = false;
        if has_oo_edge {
            oo_paredge = if on2edges {
                vp.edge_parameter(oo_shape_index)
            } else {
                vp.edge_on_parameter(oo_shape_index)?
            };
            let oo_e = if on2edges {
                vp.edge(oo_shape_index)
                    .ok_or_else(|| "ProcessVPonR : no OO edge".to_string())?
                    .clone()
            } else {
                vp.edge_on(oo_shape_index)?.clone()
            };
            oo_edge = oo_e;
            if hds.ds().has_shape(&oo_edge, true) {
                oo_edge_index = hds.ds().shape_index(&oo_edge, true);
            }
            oo_isrest = hds.ds().is_section_edge(&oo_edge, true);
            oo_closing = crate::topopebrep::tool::shape_tool::closed(&oo_edge, &oo_face);
            dg_oo_e = BRepTool::is_degenerated(&Edge(oo_edge.clone()));
        }

        // degenerated edge processing
        let (_dg, idx) = self.process_vp_on_dg_e(
            facint,
            hds,
            vp,
            shape_index,
            Kind::Point,
            0, // `cxx:970`: `PVKind = POINT; PVIndex = 0;`
        )?;
        let pv_kind = idx.pv_kind;
        let pv_index = idx.pv_index;
        let found_pv_index = pv_index != 0;

        Ok(Some(VpOnRState {
            shape_index,
            oo_shape_index,
            face: face.clone(),
            oo_face,
            absindex,
            on2edges,
            has_on_edge,
            has_oo_edge,
            i_si_face,
            i_oo_face,
            si_edge_index,
            edge,
            paredge,
            isrest,
            closing,
            dge,
            oo_edge_index,
            oo_edge,
            oo_paredge,
            oo_isrest,
            oo_closing,
            dg_oo_e,
            pv_kind,
            pv_index,
            epi_found: idx.epi_found,
            iepi: idx.iepi,
            cpi_found: idx.cpi_found,
            icpi: idx.icpi,
            found_pv_index,
            pds: Point::new(),
            trans_line: Transition::default(),
            trans_edge: Transition::default(),
            keepvp: false,
            found: false,
        }))
    }

    /// `vpr.cxx:977-1139`: `<TransLine>`, the "no transition found" policy and
    /// the `return` it can trigger.
    ///
    /// `false` is the OCCT `return` (`cxx:1136-1138`).
    fn vp_on_r_trans_line(
        &mut self,
        facint: &mut FacesIntersector,
        vp: &VPointInter,
        st: &mut VpOnRState,
    ) -> Result<bool, String> {
        let i_vp = vp.index();
        let ltype = facint.change_line(self.line_index).type_line_curve();
        let wline = ltype == TypeLineCurve::Walking;
        let grestriction = ltype == TypeLineCurve::Restriction;
        let glinenotoned = !wline && !grestriction && !self.line_is_on_edge;

        // lasttransLine (for walking line)
        let dscilempty = self.dscil.is_empty();
        let mut setlastonwl = wline && !dscilempty;
        if setlastonwl {
            // xpu171198, FRA61896 (f7,f13-> null DSC1)
            let parline = vp.parameter_on_line()?;
            let par = self
                .dscil
                .last()
                .map(|h| fds_parameter(&h.read().expect("interference lock")).unwrap_or(0.0))
                .unwrap_or(0.0);
            let dd = (par - parline).abs(); // en fait, ce sont des entiers
            if dd == 0.0 {
                setlastonwl = false;
            }
        }
        let mut lasttrans_line = Transition::default();
        if setlastonwl {
            if let Some(h) = self.dscil.last() {
                lasttrans_line = h.read().expect("interference lock").transition();
            }
        }

        // edgeori, transLine
        let edgeori = st.edge.orientation();
        let mut trans_line =
            FFTransitionTool::process_line_transition_with_index(vp, st.shape_index, edgeori);
        let mut trliunk = trans_line.is_unknown();

        // 1_ If vpmin has transition OUT/IN, and vpmax is UNKNOWN,
        //    we change vpmax transition as IN/OUT
        //
        // 2_ If vpmin is UNKNOWN (if vp is first on line and transition is
        //    UNKNOWN, vpmin's transition is ON/ON the OOshape)
        //    we change vpmin transition as OUT/IN
        //
        // (kpart : sphere/box, with the sphere's sewing edge lying on one
        //  boxe's face and one of the edge's vertices IN the same face)
        let (i_inon1, i_inonn, n_inon) = self.line_vp_bounds(facint);

        let mut keepvp = false;
        let mut ret1 = false;
        if trliunk {
            if wline {
                if setlastonwl {
                    trans_line = lasttrans_line.complement();
                } else {
                    let l = facint.change_line(self.line_index).clone();
                    if let Some(t) = super::vpr::fun_trans_for_wl(&l, i_vp, st.shape_index)? {
                        trans_line = t;
                    }
                }

                // walki vpfirst on 3, vplast on 0, nvpkept = 2 kept
                if trans_line.is_unknown() {
                    let mut keepvpfirst = dscilempty && (i_vp == i_inon1) && (n_inon == 2);
                    if st.absindex == 3 {
                        keepvpfirst = keepvpfirst && self.last_vp_is_on0;
                    }
                    if keepvpfirst {
                        trans_line.set_orientation(Orientation::Forward);
                    }
                    ret1 = false;
                }
            } else if glinenotoned {
                let isfirstvp = i_vp == i_inon1;
                let islastvp = i_vp == i_inonn;
                if isfirstvp {
                    keepvp = true;
                } else if islastvp {
                    keepvp = !dscilempty;
                } else if !st.dge
                    && !st.dg_oo_e
                    && (vp.is_vertex_on_s1()? || vp.is_vertex_on_s2()?)
                {
                    // If VP is on vertex we should compute at least one
                    // interference for the edge. This interference is
                    // necessary at least to indicate that the edge intersects
                    // something.
                    let l = facint.change_line(self.line_index).clone();
                    let a_first_point = l
                        .v_point(i_inon1)
                        .ok_or_else(|| "ProcessVPonR : no first VPoint".to_string())?
                        .clone();
                    let a_last_point = l
                        .v_point(i_inonn)
                        .ok_or_else(|| "ProcessVPonR : no last VPoint".to_string())?
                        .clone();

                    let mut faceindex = 1;
                    while !keepvp && faceindex <= 2 {
                        let vp_is_vertex = if faceindex == 1 {
                            vp.is_vertex_on_s1()?
                        } else {
                            vp.is_vertex_on_s2()?
                        };
                        let first_is_vertex = if faceindex == 1 {
                            a_first_point.is_vertex_on_s1()?
                        } else {
                            a_first_point.is_vertex_on_s2()?
                        };
                        let last_is_vertex = if faceindex == 1 {
                            a_last_point.is_vertex_on_s1()?
                        } else {
                            a_last_point.is_vertex_on_s2()?
                        };
                        if vp_is_vertex {
                            let a_v1 = if faceindex == 1 {
                                vp.vertex_on_s1()?
                            } else {
                                vp.vertex_on_s2()?
                            };
                            if first_is_vertex {
                                let a_v2 = if faceindex == 1 {
                                    a_first_point.vertex_on_s1()?
                                } else {
                                    a_first_point.vertex_on_s2()?
                                };
                                if crate::topo_tools::is_same(a_v1, a_v2) {
                                    keepvp = true;
                                }
                            }
                            if !keepvp && last_is_vertex {
                                let a_v2 = if faceindex == 1 {
                                    a_last_point.vertex_on_s1()?
                                } else {
                                    a_last_point.vertex_on_s2()?
                                };
                                if crate::topo_tools::is_same(a_v1, a_v2) {
                                    keepvp = !dscilempty;
                                }
                            }
                        }
                        faceindex += 1;
                    }
                }
                ret1 = !keepvp;
            } else {
                ret1 = true;
            }
        }
        trliunk = trans_line.is_unknown();
        st.trans_line = trans_line;
        st.keepvp = keepvp;
        if ret1 {
            return Ok(false);
        }
        let _ = trliunk; // read again below by the caller in OCCT
        Ok(true)
    }

    /// `vpr.cxx:1141-1238`: `<Transori>`, `<transEdge>` and the
    /// `FUN_newtransEdge` call.
    ///
    /// `false` is the OCCT `return` triggered by `ret2` (`cxx:1236-1238`).
    fn vp_on_r_trans_edge(
        &mut self,
        facint: &mut FacesIntersector,
        hds: &mut HDataStructure,
        vp: &VPointInter,
        trans: &Transition,
        st: &mut VpOnRState,
    ) -> Result<bool, String> {
        let transori = trans.orientation(State::In, ShapeType::Face);
        let mut trans_edge =
            FFTransitionTool::process_edge_transition(vp, st.shape_index, transori);
        let t_unknown = fds_has_unk(&trans_edge);
        let pds = PointGeomTool::make_point_of_v_point_inter(vp)?; // <VP>'s geometry

        // lITOOFonVP = {I on <edge> = (T on <OOface>, G on <VP>, S)}
        let mut l_it_oo_f_on_vp: Vec<InterferenceHandle> = Vec::new();
        let mut found = false;
        if st.si_edge_index != 0 {
            let l_i_edge = hds.ds().shape_interferences(&st.edge, true).to_vec();
            let mut l_i: Vec<InterferenceHandle> = Vec::new();
            if st.pv_index == 0 {
                fun_scan_interf_list(&pds, hds, &l_i_edge, &mut l_i)?;
            } else {
                fun_select_g_interference(&l_i_edge, st.pv_index, &mut l_i);
            }
            found = fun_select_traisha_interference(&l_i, st.i_oo_face, &mut l_it_oo_f_on_vp)?;
        }

        // <Transori> = INTERNAL or EXTERNAL (tangent cases), compute <transEdge>
        let mut newtrans_edge = transori == Orientation::Internal || transori == Orientation::External;
        let o_trans_edge = trans_edge.orientation(State::In, ShapeType::Face);
        let all_int = transori == Orientation::Internal || o_trans_edge == Orientation::Internal;
        let all_ext = transori == Orientation::External || o_trans_edge == Orientation::External;

        newtrans_edge = newtrans_edge && !all_int && !all_ext;
        newtrans_edge = newtrans_edge || t_unknown;
        // -> intersection fails for closing edges
        // 1. <edge> touches closing <OOedge> at <VP> && <OOedge> is tangent to
        //    <OOFace>, intersection -> 1 forward && 1 reversed instead of
        //    internal/external.
        // 2. if <edge> is tangent to <OOFace> at <VP> on walking
        newtrans_edge = newtrans_edge || st.closing || st.oo_closing;
        newtrans_edge = newtrans_edge && !self.line_inl;
        if newtrans_edge {
            // `myDS->Shape(edge);` (`cxx:1186`): the lookup result is dropped.
            let _ = hds.ds().shape_index(&st.edge, true);
            newtrans_edge = !found;
            if found {
                // Getting first transition found
                // prequesitory : transition on edge / <OOF> on same geometry
                // point is unchanged
                let tr = l_it_oo_f_on_vp
                    .first()
                    .expect("found implies non-empty")
                    .read()
                    .expect("interference lock")
                    .transition();
                trans_edge.before(tr.before_state(), ShapeType::Face);
                trans_edge.after(tr.after_state(), ShapeType::Face);
            }
            if newtrans_edge {
                // Compute of <transEdge> : transition on <edge> at geometry
                // <VP> / <OOface>; if the line is on a restriction OOedge of
                // <OOface>, gets <edge> transition/<OOface> when at <VP> on
                // OOedge, elsewhere when the line is not on a restriction.
                let mut tr = Transition::default();
                let erl = self.erl.clone();
                let ok = fun_newtrans_edge(
                    hds,
                    facint,
                    self.line_index,
                    self.line_is_on_edge,
                    vp,
                    st.pv_kind,
                    st.pv_index,
                    st.oo_shape_index,
                    &st.oo_face,
                    &st.edge,
                    &erl,
                    &mut tr,
                )?;
                if ok {
                    trans_edge.before(tr.before_state(), ShapeType::Face);
                    trans_edge.after(tr.after_state(), ShapeType::Face);
                }
                newtrans_edge = ok;
            }
        } // newtrans_edge

        let tredunk = trans_edge.is_unknown();
        let mut ret2 = false;
        if tredunk {
            if !st.trans_line.is_unknown() {
                trans_edge = st.trans_line.complement();
            }
            if st.trans_line.is_unknown() && !st.keepvp {
                ret2 = true;
            }
        }
        let _ = newtrans_edge; // `cxx:1220` last write, never read.
        st.trans_edge = trans_edge;
        st.pds = pds;
        st.found = found;
        Ok(!ret2)
    }
}
