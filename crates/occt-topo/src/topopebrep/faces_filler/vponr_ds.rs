//! `TopOpeBRep_vpr.cxx:1238-1705` — the DS geometry and interference half of
//! `ProcessVPonR`.
//!
//! Split out of [`super::vponr`] at the "DS geometry Management" separator
//! (`cxx:1240-1243`); the state both halves share is
//! [`super::vponr::VpOnRState`].
//!
//! Sources: `ProcessVPonR` (`TopOpeBRep_vpr.cxx:1238-1730`).
//!
//! Faithfulness notes kept from the C++:
//! - `EsdmEofF = true` is set again inside the `TopExp_Explorer` loop
//!   (`cxx:1516-1525`) although it is already `true` for the loop to run — the
//!   loop cannot change the flag. Kept verbatim.
//! - `OOrest = FUN_INlos(edge, myERL)` uses `<edge>`, not `<OOedge>`
//!   (`cxx:1534`).
//! - `transEdge.Before(T.Before())` (`cxx:1377-1378`) uses the OCCT overload
//!   whose `ShapeBefore` defaults to `TopAbs_FACE`.

use occt_core::gp::GpPnt2d;
use occt_core::precision::ANGULAR;

use crate::abs::{Orientation, ShapeType, State};
use crate::shape::Face;
use crate::topopebrep::ds::export::{fds_has_same_domain3d, fds_loi_inf_sup, fun_ds_sdm};
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::ds::kind::Kind;
use crate::topopebrep::ds::process_interferences_tool::{
    make_epv_interference, make_epv_interference_with_support,
};
use crate::topopebrep::ds::transition::Transition;
use crate::topopebrep::faces_filler::vponr::{fun_on_edge, VpOnRState};
use crate::topopebrep::faces_filler::vp_index::{fun_vp_index, M_MKNEWVP};
use crate::topopebrep::faces_filler::vpr::fun_in_los;
use crate::topopebrep::faces_filler::vpr_cpi::fun_process_cpi;
use crate::topopebrep::faces_intersector::FacesIntersector;
use crate::topopebrep::interf::v_point_inter::VPointInter;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;
use crate::topopebrep::tool::geometry::fun_tool_onapex;
use crate::topopebrep::tool::make_transition::MakeTransition;
use crate::topopebrep::tool::tangent::fun_tool_tg_geom_e;
use crate::topopebrep::tool::topology::{fun_tool_etg_f, fun_tool_ng_geom_f, fun_tool_par_on_e_f};

impl super::FacesFiller {
    /// `vpr.cxx:1238-1385`: the DS shape indices, `FUN_VPIndex`, and the
    /// `myLineINL` / `VPonedge` early exit.
    ///
    /// `false` is the OCCT `return` at `cxx:1369`.
    pub(super) fn vp_on_r_geometry(
        &mut self,
        facint: &mut FacesIntersector,
        hds: &mut HDataStructure,
        vp: &VPointInter,
        st: &mut VpOnRState,
    ) -> Result<bool, String> {
        if st.si_edge_index == 0 {
            st.si_edge_index = hds
                .change_ds()
                .add_shape_with_ancestor(&st.edge, st.shape_index);
        }

        let si_is_vertex = vp.is_vertex(st.shape_index)?;
        let oo_is_vertex = vp.is_vertex(st.oo_shape_index)?;

        // modified by NIZHNY-MKK  Tue Apr  3 12:08:38 2001
        let ismultiplekind = st.found_pv_index
            && !st.epi_found
            && !st.cpi_found
            && (si_is_vertex || oo_is_vertex)
            && (st.pv_kind == Kind::Point);

        if !st.found_pv_index || ismultiplekind {
            let dscil = self.dscil.clone();
            let idx = fun_vp_index(
                self,
                facint,
                self.line_index,
                vp,
                st.shape_index,
                hds,
                &dscil,
                M_MKNEWVP,
            )?;
            st.pv_kind = idx.pv_kind;
            st.pv_index = idx.pv_index;
            st.epi_found = idx.epi_found;
            st.iepi = idx.iepi;
            st.cpi_found = idx.cpi_found;
            st.icpi = idx.icpi;
        }
        // `if (PVIndex == 0) funraise();` (`cxx:1284-1288`) is OCCT_DEBUG only.

        let mut vp_onedge = false;
        if st.pv_kind == Kind::Vertex {
            vp_onedge = fun_on_edge(&st.pds, &st.edge);
        }

        if self.line_inl {
            let tolang = ANGULAR * 1.0e5; //=1.e-7 NYITOLXPU
            let tg_e = fun_tool_tg_geom_e(st.paredge, &st.edge);

            let mut oo_uv = GpPnt2d::new(0.0, 0.0);
            let mut ok = false;
            if vp_onedge {
                oo_uv = vp.surface_parameters(st.oo_shape_index)?;
                ok = true;
            } else if let Some(uv) = fun_tool_par_on_e_f(&st.oo_edge, st.oo_paredge, &st.oo_face) {
                oo_uv = uv;
                ok = true;
            }

            let mut nt_oof: Option<occt_core::gp::GpDir> = None;
            if ok {
                nt_oof = fun_tool_ng_geom_f(&oo_uv, &Face(st.oo_face.clone()));
            }
            if st.oo_face.orientation() == Orientation::Reversed {
                if let Some(v) = nt_oof.as_mut() {
                    v.reverse();
                }
            }

            let tol = 1.0e-7;
            if ok {
                // A null `gp_Vec` (the C++ default) has a zero magnitude, so
                // the `!ok` arm always fails this test.
                let nt_mag = nt_oof
                    .as_ref()
                    .map(|v| (v.x() * v.x() + v.y() * v.y() + v.z() * v.z()).sqrt())
                    .unwrap_or(0.0);
                ok = (tg_e.magnitude() > tol) && (nt_mag > tol);
            }
            let mut dot = 1.0e7;
            if ok {
                let nt = nt_oof.expect("ok implies a normal");
                dot = tg_e
                    .normalized()
                    .dot(&occt_core::gp::GpVec::new(nt.x(), nt.y(), nt.z()));
            }

            let oo_face = Face(st.oo_face.clone());
            let sa = SurfaceAdaptor::of_face(&oo_face);
            let su = sa.surface();
            let apex = match su.as_ref() {
                Some(s) => fun_tool_onapex(&oo_uv, s.as_ref()),
                None => false,
            };

            let mut t = Transition::default();
            if !apex && ok && dot.abs() > tolang {
                let ori = if dot < 0.0 {
                    Orientation::Forward
                } else {
                    Orientation::Reversed
                };
                t.set_orientation(ori);
            }

            if vp_onedge && !st.dge {
                // xpu231098 : cto904C8(edge11)
                // xpu131198 : CTS21802(edge31)
                if st.i_oo_face == 0 {
                    st.i_oo_face = hds
                        .change_ds()
                        .add_shape_with_ancestor(&st.oo_face, st.oo_shape_index);
                }
                {
                    let mut tt = t.clone();
                    tt.set_index(st.i_oo_face);
                    if let Some(epi) = make_epv_interference_with_support(
                        &tt,
                        st.i_oo_face,
                        st.pv_index,
                        st.paredge,
                        st.pv_kind,
                        Kind::Face,
                        si_is_vertex,
                    ) {
                        hds.store_interference_on_shape(&epi, &st.edge)?;
                    }
                }
                if st.on2edges || st.has_on_edge {
                    if st.oo_edge_index == 0 {
                        st.oo_edge_index = hds
                            .change_ds()
                            .add_shape_with_ancestor(&st.oo_edge, st.oo_shape_index);
                    }
                    let mut tt = t.clone();
                    tt.set_index(st.i_oo_face);
                    if let Some(epi) = make_epv_interference(
                        &tt,
                        st.oo_edge_index,
                        st.pv_index,
                        st.paredge,
                        st.pv_kind,
                        si_is_vertex,
                    ) {
                        hds.store_interference_on_shape(&epi, &st.edge)?;
                    }
                }
                return Ok(false);
            } else if !t.is_unknown() {
                // compute interferences later on
                // modified by NIZHNY-MZV  Thu Dec 23 13:27:10 1999
                st.trans_edge.before(t.before_state(), t.shape_before());
                st.trans_edge.after(t.after_state(), t.shape_after());
            }
        } // myLineINL

        Ok(true)
    }

    /// `vpr.cxx:1387-1705`: the CPI / EPI / ON-correction / closing-edge pass.
    pub(super) fn vp_on_r_store(
        &mut self,
        facint: &mut FacesIntersector,
        hds: &mut HDataStructure,
        vp: &VPointInter,
        st: &mut VpOnRState,
    ) -> Result<(), String> {
        let si_is_vertex = vp.is_vertex(st.shape_index)?;
        let oo_is_vertex = vp.is_vertex(st.oo_shape_index)?;

        // Gfound = VP corresponds with an existing geometry of ShapeIndex
        let gfound = st.epi_found || st.cpi_found;

        // ------------------------------------------
        // -- Curve/(POINT,VERTEX) Interference (CPI)
        // ------------------------------------------
        let mut no_cpi = self.line_is_on_edge;
        no_cpi = no_cpi || (!st.on2edges && st.has_oo_edge && (st.oo_isrest || st.isrest));

        let condi = (!no_cpi) && (!self.line_inl); // INL
        if condi {
            let mut kept = self.kept_vp_nbr;
            fun_process_cpi(
                self,
                facint,
                hds,
                vp,
                &st.face,
                st.shape_index,
                &st.trans_line,
                st.icpi.as_ref(),
                gfound,
                st.pv_kind,
                st.pv_index,
                &mut kept,
            )?;
            self.kept_vp_nbr = kept;
        }

        // ------------------------------------------
        // --- Edge/(POINT,VERTEX) Interference (EPI)
        // ------------------------------------------
        let condi2 = (st.on2edges && !st.closing) || (st.has_on_edge && !st.closing);
        if condi2 && (!st.dge) {
            if st.oo_edge_index == 0 {
                st.oo_edge_index = hds
                    .change_ds()
                    .add_shape_with_ancestor(&st.oo_edge, st.oo_shape_index);
            }
            let mut t = st.trans_edge.clone();
            if st.i_oo_face == 0 {
                st.i_oo_face = hds
                    .change_ds()
                    .add_shape_with_ancestor(&st.oo_face, st.oo_shape_index);
            }
            t.set_index(st.i_oo_face);
            if let Some(epi) = make_epv_interference(
                &t,
                st.oo_edge_index,
                st.pv_index,
                st.paredge,
                st.pv_kind,
                si_is_vertex,
            ) {
                hds.store_interference_on_shape(&epi, &st.edge)?;
            }
        }

        // ===================================================================
        // manip corrective d'un pb. d'intersection
        // - le VPoint est donne sur une restriction de ShapeIndex (1 ou 2),
        // - le VPoint est ON une restriction de l'autre shape (OOShapeIndex)
        // --> le VPoint n'est PAS donne sur restriction de OOShapeIndex NYI
        //     par les intersections (a ameliorer). NYI
        // ===================================================================
        let mut correct_on = !st.on2edges && st.has_on_edge && !st.dg_oo_e;
        let mut corrected_on = false;
        if correct_on && st.oo_edge_index != 0 {
            let l_i_oo_edge = hds.ds().shape_interferences(&st.oo_edge, true).to_vec();
            let mut l_i: Vec<_> = Vec::new();
            crate::topopebrep::faces_filler::vpr_trans_edge::fun_scan_interf_list(
                &st.pds,
                hds,
                &l_i_oo_edge,
                &mut l_i,
            )?;
            let mut l_i_tf_on_vp: Vec<_> = Vec::new();
            let oo_found = crate::topopebrep::faces_filler::vpr_trans_edge::
                fun_select_traisha_interference(&l_i, st.i_si_face, &mut l_i_tf_on_vp)?;
            correct_on = !oo_found;
        }

        if correct_on {
            if st.oo_edge_index == 0 {
                st.oo_edge_index = hds
                    .change_ds()
                    .add_shape_with_ancestor(&st.oo_edge, st.oo_shape_index);
            }

            // VP a ete classifie ON sur l'edge <OOedge>.
            // calcul de la transition <tOOedge> sur l'arete <OOedge> (de
            // l'autre face en jeu, OOShapeIndex) ou le VP est donne ON.
            // On tient compte de l'orientation de <edge> dans <Face>.
            let mut t_oo_edge = Transition::default();

            // distinguish whether OOedge is the edge on which geometric line
            // lies.
            // OOedge == edge(line) ==> tOOedge = f(orientation of <edge> in
            // <Face> FORWARD)
            // OOedge != edge(line) ==> tOOedge = f(orientation of <Face>)
            let loi = hds.ds().shape_interferences(&st.oo_edge, true).to_vec();
            let (oo_par1, oo_par2, _isonper) = fds_loi_inf_sup(
                &st.oo_edge,
                st.oo_paredge,
                st.pv_kind,
                st.pv_index,
                &loi,
            )?;

            let oo_uv = vp.surface_parameters(st.shape_index)?;

            // <Tr> relative to 3d <OOface> matter, we take into account
            // <Tr> / 2d <OOface> only if <edge> is normal to <OOface>
            let tola = ANGULAR * 1.0e+2; // dealing with tolerances

            // KK : supplying tolerances pbm (tola too small)
            let mut esdm_eof_f = hds.has_same_domain(&st.oo_edge, true);
            if esdm_eof_f {
                let mut ex = crate::topexp::Explorer::new(&st.face, crate::abs::ShapeType::Edge);
                while ex.more() {
                    if fun_ds_sdm(hds.ds(), &ex.current(), &st.oo_edge) {
                        esdm_eof_f = true;
                        break;
                    }
                    ex.next();
                }
            }
            let mut oo_etg_f = true;
            if !esdm_eof_f {
                oo_etg_f = fun_tool_etg_f(st.oo_paredge, &st.oo_edge, &oo_uv, &st.face, tola);
            }
            let oo_rest = fun_in_los(&st.edge, &self.erl);
            let interf2d = oo_etg_f && (st.oo_isrest || oo_rest);

            let factor = 1.0e-2;
            let mut mkt = MakeTransition::new();
            let ok = mkt.initialize(
                &st.oo_edge,
                oo_par1,
                oo_par2,
                st.oo_paredge,
                &st.face,
                &oo_uv,
                factor,
            );

            if ok && (!interf2d || mkt.is_t2d()) {
                let _ = mkt.set_rest(&st.edge, st.paredge);
                let mut stb = State::Unknown;
                let mut sta = State::Unknown;
                if mkt.mk_ton_e(&mut stb, &mut sta) {
                    t_oo_edge.before(stb, ShapeType::Face);
                    t_oo_edge.after(sta, ShapeType::Face);
                    {
                        if st.i_si_face == 0 {
                            st.i_si_face = hds
                                .change_ds()
                                .add_shape_with_ancestor(&st.face, st.shape_index);
                        }
                        let mut oot = t_oo_edge.clone();
                        oot.set_index(st.i_si_face);
                        if let Some(epi) = make_epv_interference(
                            &oot,
                            st.si_edge_index,
                            st.pv_index,
                            st.oo_paredge,
                            st.pv_kind,
                            oo_is_vertex,
                        ) {
                            hds.store_interference_on_shape(&epi, &st.oo_edge)?;
                        }
                    }

                    // xpu : 09-03-98
                    // hsd3d => interf2d  : only IwithSkEDGE interf
                    // elsewhere            : add an IwithSkFACE interference.
                    let mut add_epi_f = !self.line_is_on_edge;
                    let (hsd3d, _dummy) = fds_has_same_domain3d(hds.ds(), &st.oo_edge, true);
                    if hsd3d {
                        add_epi_f = false;
                    }
                    if add_epi_f {
                        let mut oot = t_oo_edge.clone();
                        oot.set_index(st.i_si_face);
                        if let Some(epi) = make_epv_interference_with_support(
                            &oot,
                            st.i_si_face,
                            st.pv_index,
                            st.oo_paredge,
                            st.pv_kind,
                            Kind::Face,
                            oo_is_vertex,
                        ) {
                            hds.store_interference_on_shape(&epi, &st.oo_edge)?;
                        }
                    }
                    // xpu : 09-03-98
                    corrected_on = true;
                } // ok
            }
        } // correctON

        if correct_on && !corrected_on && no_cpi && !self.line_is_on_edge {
            // MSV: correct ON failed, so store CPI
            let mut kept = self.kept_vp_nbr;
            fun_process_cpi(
                self,
                facint,
                hds,
                vp,
                &st.face,
                st.shape_index,
                &st.trans_line,
                st.icpi.as_ref(),
                gfound,
                st.pv_kind,
                st.pv_index,
                &mut kept,
            )?;
            self.kept_vp_nbr = kept;
        }

        // closing edge processing
        if (st.oo_closing || st.closing) && !st.found {
            self.process_vp_on_closing_r(
                facint,
                hds,
                vp,
                st.shape_index,
                &st.trans_edge,
                st.pv_kind,
                st.pv_index,
            )?;
            return Ok(());
        }

        // VP processing
        let mut add_epi = false;
        if !gfound {
            add_epi = true;
        } else {
            // EPIfound
            let an_o_trans_edge = st.trans_edge.orientation(State::In, ShapeType::Face);

            let mut opporifound = false;
            let mut memorifound = false;
            let mut lo_i_found: Vec<_> = Vec::new();
            let epil = hds.ds().shape_interferences(&st.edge, true).to_vec();
            let ok = crate::topopebrep::faces_filler::vpr_trans_edge::fun_same_g_same_s(
                &epil,
                st.pv_index,
                st.i_oo_face,
                &mut lo_i_found,
            );
            if ok {
                // on cree une EPI orientee <transEdge> ssi :
                // - il en existe deja une d'orientation opposee a TransEdge
                // - il n'en existe pas deja une d'orientation identique a
                //   TransEdge
                for h in lo_i_found.iter() {
                    let o_epi = h
                        .read()
                        .expect("interference lock")
                        .transition()
                        .orientation(State::In, ShapeType::Face);
                    if !memorifound {
                        memorifound = o_epi == an_o_trans_edge;
                    }
                    if !opporifound {
                        opporifound = o_epi == an_o_trans_edge.complement();
                    }
                    add_epi = opporifound && !memorifound;
                    if add_epi {
                        break;
                    }
                }
            }
            if !ok {
                add_epi = true;
            }
        } // EPIfound

        // xpu030998 : edge has restriction on OOface, do NOT append EPIf
        //     cto904A3 (edge19,OOface14,vG16),
        if self.line_inl {
            let tola = ANGULAR * 1.0e+4; // dealing with tolerances
            let uv = vp.surface_parameters(st.oo_shape_index)?;
            let etg_oof = fun_tool_etg_f(st.paredge, &st.edge, &uv, &st.oo_face, tola);
            let in_erl = fun_in_los(&st.edge, &self.erl);
            if etg_oof && in_erl {
                return Ok(()); // cto904A3
            }
        }

        if add_epi && (!st.dge) {
            // ShapeIndex = 1,2 --> OOShapeIndex = 2,1
            // point est sur une seule arete <edge> de <ShapeIndex>
            // le Support de l'interference est l'autre face (OOShapeIndex) /
            // ShapeIndex
            if st.i_oo_face == 0 {
                st.i_oo_face = hds
                    .change_ds()
                    .add_shape_with_ancestor(&st.oo_face, st.oo_shape_index);
            }
            let mut t = st.trans_edge.clone();
            t.set_index(st.i_oo_face);
            if let Some(epi) = make_epv_interference_with_support(
                &t,
                st.i_oo_face,
                st.pv_index,
                st.paredge,
                st.pv_kind,
                Kind::Face,
                si_is_vertex,
            ) {
                hds.store_interference_on_shape(&epi, &st.edge)?;
            }
        } // addEPI

        Ok(())
    }
}
