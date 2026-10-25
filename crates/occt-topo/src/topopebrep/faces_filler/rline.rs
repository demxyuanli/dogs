//! `TopOpeBRep_FacesFiller_1.cxx:647-1125` — `ProcessRLine`, the VPoint pass for
//! a `TopOpeBRep_RESTRICTION` line.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRep/TopOpeBRep_FacesFiller_1.cxx`.
//!
//! Called from [`super::line`] `FillLineVPonR` (`_1.cxx:1158`), after
//! `FUN_GetdgData` / `FUN_FillVof12` have filled `myDataforDegenEd`.
//!
//! Faithfulness notes:
//! - `processUNK` is `false` in OCCT (`_1.cxx:682`: `bool processUNK = false;`,
//!   only re-enabled from `TopOpeBRep_GetcontextNOPUNK()` under `OCCT_DEBUG`,
//!   which is compiled out here), so the `T1unk` / `T2unk` block
//!   (`_1.cxx:687-745`) is never taken. It is not translated; the two `Tunk`
//!   flags it would compute are still read afterwards, but the block is the
//!   only writer that could change them.
//! - The `if (!on2edges)` arm (`_1.cxx:775-846`) is a `continue` whose whole
//!   body is commented out in OCCT; the port keeps the `continue`.
//! - `sdmeds` (`_1.cxx:856`) is `continue`d on (`:858-860`), so the `if
//!   (sdmeds)` arm of the interference filter (`:1070-1081`) is unreachable.
//!   The port keeps the `else if (iOOFace != 0)` arm only.
//! - `PVKind` / `PVIndex` are declared outside the `ShapeIndex` loop
//!   (`_1.cxx:880`) and therefore carry over from iteration 1 to iteration 2;
//!   the port keeps them across the loop.
//! - `VPI.Init((*myLine))` uses the default `checkkeep = false`
//!   (`TopOpeBRep_VPointInterIterator.hxx:36`), so every VPoint is visited;
//!   `vpf` is `CurrentVP()` before the loop, i.e. `VPoint(1)`.

use crate::abs::{Orientation, ShapeType, State};
use crate::shape::{Edge, TopoShape};
use crate::topopebrep::ds::export::{fds_has_unk, fun_ds_sdm};
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::ds::interference_tool::make_face_edge_interference;
use crate::topopebrep::ds::kind::{Config, Kind};
use crate::topopebrep::ds::process_interferences_tool::{
    fun_select_gi_interference, fun_select_itrasha_interference, fun_select_si_interference,
    fun_select_traori_interference, make_epv_interference, make_epv_interference_with_support,
};
use crate::topopebrep::ds::transition::Transition;
use crate::topopebrep::faces_filler::vponr::fun_brep_on_first_p;
use crate::topopebrep::faces_filler::vp_index::{fun_vp_index, M_GETVP};
use crate::topopebrep::faces_intersector::FacesIntersector;
use crate::topopebrep::ff_transition_tool::FFTransitionTool;
use crate::topopebrep::interf::v_point_inter::VPointInter;
use crate::topopebrep::tool::closed::{on_boundary, FORWARD, REVERSED};
use crate::topopebrep::tool::{fun_tool_etg_oo_e, fun_tool_orient_ein_f_forward};
use crate::topopebrep::type_line_curve::TypeLineCurve;

impl super::FacesFiller {
    /// `ProcessRLine()` (`_1.cxx:647-1125`).
    pub fn process_r_line(
        &mut self,
        facint: &mut FacesIntersector,
        hds: &mut HDataStructure,
    ) -> Result<(), String> {
        // `_1.cxx:649-652`
        if facint.change_line(self.line_index).type_line_curve() != TypeLineCurve::Restriction {
            return Ok(());
        }

        // `_1.cxx:654-658`: `addIFE` is a constant `true`.
        let add_ife = true;
        if !add_ife {
            return Ok(());
        }

        // `_1.cxx:660-666`
        let erest = match facint.change_line(self.line_index).arc() {
            Some(e) => e.clone(),
            None => return Ok(()),
        };
        let fi_is_rest = facint.is_restriction(&erest);
        if !fi_is_rest {
            return Ok(());
        }

        // `_1.cxx:668-672`
        let isedge1 = facint.change_line(self.line_index).arc_is_edge(1);
        let isedge2 = facint.change_line(self.line_index).arc_is_edge(2);
        let e_shape_index = if isedge1 {
            1
        } else if isedge2 {
            2
        } else {
            0
        };

        let i_erest = hds.change_ds().add_shape_with_ancestor(&erest, e_shape_index);
        // `_1.cxx:673-674`: `rank = myDS->AncestorRank(iErest)`.
        let rank = hds.ds().ancestor_rank(i_erest);
        let oo_rank = if rank == 1 { 2 } else { 1 };

        // `_1.cxx:676-677`
        let f1 = self.face(1)?.0.clone();
        let f2 = self.face(2)?.0.clone();
        let i_f1 = hds.change_ds().add_shape_with_ancestor(&f1, 1);
        let i_f2 = hds.change_ds().add_shape_with_ancestor(&f2, 2);

        // `_1.cxx:679-683`
        let mut t1 = self.face_face_transition(facint, 1)?.clone();
        t1.set_index(i_f2);
        let mut t2 = self.face_face_transition(facint, 2)?.clone();
        t2.set_index(i_f1);

        // `_1.cxx:685-745`: `processUNK` is `false` in OCCT (see the module
        // documentation), so this block has no effect and is not translated.

        // `_1.cxx:742-753`
        let ife = make_face_edge_interference(
            &t1,
            i_f2,
            i_erest,
            isedge1,
            Config::UnshGeometry,
        );
        hds.store_interference_on_index(&ife, i_f1)?;
        let ife = make_face_edge_interference(
            &t2,
            i_f1,
            i_erest,
            isedge2,
            Config::UnshGeometry,
        );
        hds.store_interference_on_index(&ife, i_f2)?;

        // `_1.cxx:757-760`
        let nb_vp = facint.change_line(self.line_index).nb_v_point();
        let vpf: Option<VPointInter> = facint
            .change_line(self.line_index)
            .v_point(1)
            .cloned();
        // NYIXPUTOL
        let tola = occt_core::precision::ANGULAR * 1.0e5;

        for k in 1..=nb_vp {
            let Some(vp) = facint.change_line(self.line_index).v_point(k).cloned() else {
                continue;
            };
            // `_1.cxx:765-771`
            let absindex = vp.shape_index(); // 0,1,2,3
            let mut par_rest = 0.0;
            let ok_r = vp.paron_e(&Edge(erest.clone()), &mut par_rest);
            if !ok_r {
                par_rest = vp.parameter_on_line()?;
            }
            let on2edges = absindex == 3 || absindex == oo_rank;

            // `_1.cxx:775-847`: the OCCT body is commented out; only the
            // `continue` is live.
            if !on2edges {
                continue;
            }

            // `_1.cxx:849-852`
            let oo_e = vp
                .edge(oo_rank)
                .ok_or_else(|| "ProcessRLine : no OO edge".to_string())?
                .clone();
            let i_oo = hds.change_ds().add_shape_with_ancestor(&oo_e, oo_rank);

            // `_1.cxx:854-855`
            let mut oo_par = 0.0;
            vp.paron_e(&Edge(oo_e.clone()), &mut oo_par);

            // `_1.cxx:857-860`: 1d interference already done by EdgesFiller.
            let sdmeds = fun_ds_sdm(hds.ds(), &erest, &oo_e);
            if sdmeds {
                continue;
            }

            // `_1.cxx:863-864`
            let ob_rest = on_boundary(par_rest, &erest);
            let ob_oo = on_boundary(oo_par, &oo_e);

            // `_1.cxx:867-877`: the `EXTERNAL` diagnostics are `OCCT_DEBUG`
            // only.

            // `_1.cxx:879`
            let tgeds = fun_tool_etg_oo_e(par_rest, &erest, oo_par, &oo_e, tola);

            // `_1.cxx:880-881`: carried over between the two `ShapeIndex`
            // iterations.
            let mut pv_kind = Kind::Point;
            let mut pv_index = 0;

            for shape_index in 1..=2 {
                let oo_shape_index = if shape_index == 1 { 2 } else { 1 };
                let si_erest = shape_index == rank;

                // `_1.cxx:888-891`
                let si_isvertex = vp.is_vertex(shape_index)?;
                let oo_isvertex = vp.is_vertex(oo_shape_index)?;
                let oo_face = self.face(oo_shape_index)?.0.clone();
                let mut i_oo_face = hds.ds().shape_index(&oo_face, true);

                // `_1.cxx:893-912`
                let (edge, paredge, onbound, oo_edge, mut oo_edge_index) = if si_erest {
                    (erest.clone(), par_rest, ob_rest, oo_e.clone(), i_oo)
                } else {
                    (oo_e.clone(), oo_par, ob_oo, erest.clone(), i_erest)
                };

                // `_1.cxx:916-927`, `ProcessVPondgE`. Its return value is not
                // read by OCCT here. `IEPI` / `ICPI` (`_1.cxx:917`) are
                // written by it and never read, so they are not carried.
                let (_made, idx) =
                    self.process_vp_on_dg_e(facint, hds, &vp, shape_index, pv_kind, pv_index)?;
                pv_kind = idx.pv_kind;
                pv_index = idx.pv_index;

                // `_1.cxx:929-963`
                if pv_index == 0 {
                    let idx = fun_vp_index(
                        self,
                        facint,
                        self.line_index,
                        &vp,
                        shape_index,
                        hds,
                        &self.dscil,
                        M_GETVP,
                    )?;
                    pv_kind = idx.pv_kind;
                    pv_index = idx.pv_index;
                    let efound = idx.epi_found || idx.cpi_found;
                    let ifound = pv_index != 0;
                    let condmake = !efound && !ifound;
                    if condmake {
                        let si = if si_isvertex {
                            shape_index
                        } else if oo_isvertex {
                            oo_shape_index
                        } else {
                            shape_index
                        };
                        let (i, kd) = Self::make_geometry(hds.change_ds(), &vp, si)?;
                        pv_index = i;
                        pv_kind = kd;
                    }
                }

                // `_1.cxx:965-968`
                if oo_edge_index == 0 {
                    oo_edge_index =
                        hds.change_ds()
                            .add_shape_with_ancestor(&oo_edge, oo_shape_index);
                }

                // `_1.cxx:970-1052`
                let trans_edge =
                    match self.r_line_trans_edge(
                        facint, hds, &vp, &vpf, &erest, &oo_e, &oo_face, absindex, shape_index,
                        si_erest, oo_shape_index, onbound, tgeds, pv_kind, pv_index,
                    )? {
                        Some(t) => t,
                        None => continue,
                    };

                // `_1.cxx:1053-1057`
                if fds_has_unk(&trans_edge) {
                    continue;
                }

                // `_1.cxx:1059-1097`
                let otrans_edge = trans_edge.orientation(State::In, ShapeType::Face);
                let mut copy = hds.ds().shape_interferences(&edge, true).to_vec();
                let mut l1 = Vec::new();
                let mut l2 = Vec::new();
                let mut nfound = fun_select_gi_interference(&mut copy, pv_index, &mut l1);
                if oo_edge_index != 0 {
                    nfound = fun_select_si_interference(&mut l1, oo_edge_index, &mut l2);
                }
                if nfound != 0 {
                    // `_1.cxx:1070-1081` (`sdmeds`) is unreachable: `sdmeds` is
                    // `continue`d on above.
                    if i_oo_face != 0 {
                        let mut l3 = Vec::new();
                        let mut l4 = Vec::new();
                        let n3 = fun_select_itrasha_interference(&mut l2, i_oo_face, &mut l3);
                        let mut nfound2 = n3;
                        if n3 != 0 {
                            nfound2 =
                                fun_select_traori_interference(&mut l3, otrans_edge, &mut l4);
                        }
                        if nfound2 != 0 {
                            // has I2d=(transEdge(iOOFace),PVIndex,OOedgeIndex)
                            continue;
                        }
                    }
                }

                // `_1.cxx:1100-1110`, `EPI`.
                if i_oo_face == 0 {
                    i_oo_face =
                        hds.change_ds()
                            .add_shape_with_ancestor(&oo_face, oo_shape_index);
                }
                let mut t = trans_edge.clone();
                t.set_index(i_oo_face);
                if let Some(epi) =
                    make_epv_interference(&t, oo_edge_index, pv_index, paredge, pv_kind, si_isvertex)
                {
                    hds.store_interference_on_shape(&epi, &edge)?;
                }

                // `_1.cxx:1112-1123`, `EPIf`.
                if !si_erest {
                    let mut tf = trans_edge.clone();
                    tf.set_index(i_oo_face);
                    if let Some(epif) = make_epv_interference_with_support(
                        &tf,
                        i_oo_face,
                        pv_index,
                        paredge,
                        pv_kind,
                        Kind::Face,
                        si_isvertex,
                    ) {
                        hds.store_interference_on_shape(&epif, &edge)?;
                    }
                }
            }
        }
        Ok(())
    }
}

impl super::FacesFiller {
    /// `_1.cxx:970-1052`: `<transEdge>` for the `epi_found` disjunction and
    /// the `FUN_newtransEdge` fallback.
    ///
    /// `None` is the OCCT `continue` (`:1007`, `:1045`, `:1055`); the caller
    /// also re-tests `FDS_hasUNK` (`:1053`) before using the result.
    #[allow(clippy::too_many_arguments)]
    fn r_line_trans_edge(
        &mut self,
        facint: &mut FacesIntersector,
        hds: &mut HDataStructure,
        vp: &VPointInter,
        vpf: &Option<VPointInter>,
        erest: &TopoShape,
        oo_e: &TopoShape,
        oo_face: &TopoShape,
        absindex: i32,
        shape_index: i32,
        si_erest: bool,
        oo_shape_index: i32,
        onbound: i32,
        tgeds: bool,
        pv_kind: Kind,
        pv_index: i32,
    ) -> Result<Option<Transition>, String> {
        let llt = self.face_face_transition(facint, shape_index)?.clone();
        let trans_unk = llt.is_unknown();

        let mut trans_edge = Transition::default();
        let mut tunk = true;
        if !trans_unk {
            // xpu281098 PRO12875(edge9,OOface11)
            if absindex == shape_index || absindex == 3 {
                if si_erest {
                    // transition on Erest at VP / OOface = transition at VP on
                    // line restriction
                    let ori_erest = erest.orientation();
                    trans_edge =
                        FFTransitionTool::process_line_transition_with_index(vp, shape_index, ori_erest);

                    // xpu290399 : edge is restriction, edge15,OOedge14,OOface13
                    if (onbound == FORWARD || onbound == REVERSED) && tgeds {
                        trans_edge.before(State::Unknown, ShapeType::Face);
                        trans_edge.after(State::Unknown, ShapeType::Face);
                    }
                } else {
                    // transition on edge at VP / OOface ?=
                    //   TopOpeBRep_FFTransitionTool::ProcessEdgeTransition(...)
                    // nyi (the OCCT body is commented out, `_1.cxx:996-999`)
                }
            }
            tunk = fds_has_unk(&trans_edge);
        }
        if tunk {
            if si_erest {
                // As edge=Erest is on OOFace, we only compute 2d interferences
                let onfirst_p = match vpf {
                    Some(f) => fun_brep_on_first_p(f, vp)?,
                    None => false,
                };
                let mut ovp = if onfirst_p {
                    Orientation::Forward
                } else {
                    Orientation::Reversed
                };
                let o_oo = match fun_tool_orient_ein_f_forward(oo_e, oo_face) {
                    Some(o) => o,
                    None => return Ok(None),
                };
                if o_oo == Orientation::Internal {
                    ovp = Orientation::Internal;
                }
                // xpu240399 : cto015I2 (e15,v16)
                //            edge and OOedge are tangent, we do not keep the
                //            orientation
                if !tgeds {
                    trans_edge.set_orientation(ovp);
                }
            } else {
                let erl = self.erl.clone();
                let mut tr = Transition::default();
                let ok = crate::topopebrep::faces_filler::vpr_trans_edge::fun_newtrans_edge(
                    hds,
                    facint,
                    self.line_index,
                    self.line_is_on_edge,
                    vp,
                    pv_kind,
                    pv_index,
                    oo_shape_index,
                    oo_face,
                    oo_e,
                    &erl,
                    &mut tr,
                )?;
                if !ok {
                    return Ok(None);
                }
                trans_edge.before(tr.before_state(), ShapeType::Face);
                trans_edge.after(tr.after_state(), ShapeType::Face);
            }
        }
        Ok(Some(trans_edge))
    }
}
