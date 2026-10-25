//! `TopOpeBRep_vprclo.cxx` — the closing-edge VPoint pass.
//!
//! Sources (OCCT 8.0.0, `ModelingAlgorithms/TKBool/TopOpeBRep`):
//! - `FUN_INlos` (`TopOpeBRep_vprclo.cxx:53-64`, identical to `vpr.cxx:151-162`)
//! - `GetEdgeTrans` (`vprclo.cxx:66-188`)
//! - `ProcessVPonclosingR` (`vprclo.cxx:190-296`)
//!
//! `TopOpeBRepTool_makeTransition` is `crate::topopebrep::tool::make_transition`.

use occt_core::gp::GpPnt2d;
use occt_core::precision::ANGULAR;

use crate::abs::{State, ShapeType};
use crate::shape::TopoShape;
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::ds::kind::Kind;
use crate::topopebrep::ds::transition::Transition;
use crate::topopebrep::faces_intersector::FacesIntersector;
use crate::topopebrep::interf::v_point_inter::VPointInter;
use crate::topopebrep::tool::make_transition::MakeTransition;

/// `GetEdgeTrans(VP, PVKind, PVIndex, ShapeIndex, F)` (`vprclo.cxx:66-188`):
/// `E<Sind>`'s transition on `<F>` at `<VP>`, `Unknown` when it cannot be
/// computed.
///
/// `throw Standard_Failure` (`cxx:73-77`) becomes `Err`.
#[allow(clippy::too_many_arguments)]
pub fn get_edge_trans(
    hds: &HDataStructure,
    facint: &mut FacesIntersector,
    line_index: i32,
    erl: &[TopoShape],
    vp: &VPointInter,
    pv_kind: Kind,
    pv_index: i32,
    shape_index: i32,
    f: &TopoShape,
    f1: &TopoShape,
    f2: &TopoShape,
) -> Result<Transition, String> {
    let _ = (f1, f2);
    let oo_shape_index = if shape_index == 1 { 2 } else { 1 };
    let vpsind = vp.shape_index();
    let on2edges = vpsind == 3;
    if !(on2edges || vpsind == shape_index) {
        return Err("TopOpeBRep_FacesFiller::GetEdgeTrans".to_string());
    }

    let edge = vp
        .edge(shape_index)
        .ok_or_else(|| "TopOpeBRep_FacesFiller::GetEdgeTrans : no edge".to_string())?
        .clone();
    let paredge = vp.edge_parameter(shape_index);

    let has_on_edge = vp.state(oo_shape_index)? == State::On;
    let has_oo_edge = if on2edges { true } else { has_on_edge };
    let mut ooedge: Option<TopoShape> = None;
    let mut oo_paredge = 0.0;
    if has_oo_edge {
        if on2edges {
            oo_paredge = vp.edge_parameter(oo_shape_index);
            ooedge = Some(
                vp.edge(oo_shape_index)
                    .ok_or_else(|| "GetEdgeTrans : no OO edge".to_string())?
                    .clone(),
            );
        } else {
            oo_paredge = vp.edge_on_parameter(oo_shape_index)?;
            ooedge = Some(vp.edge_on(oo_shape_index)?.clone());
        }
    }
    let oo_uv = vp.surface_parameters(oo_shape_index)?;

    // `cxx:130-148`: the interval of `<edge>` around `paredge`, narrowed by the
    // interferences already stored on it.
    let (par1, par2) = {
        let si_edge_index = if hds.ds().has_shape(&edge, true) {
            hds.ds().shape_index(&edge, true)
        } else {
            0
        };
        if si_edge_index != 0 {
            let loi = hds.ds().shape_interferences(&edge, true).to_vec();
            let (pbef, paft, _isonper) = crate::topopebrep::ds::export::fds_loi_inf_sup(
                &edge, paredge, pv_kind, pv_index, &loi,
            )?;
            (pbef, paft)
        } else {
            crate::topopebrep::tool::project_curve::fun_tool_bounds(&edge)
        }
    };

    let mut t = Transition::default();
    let tola = ANGULAR * 1.0e+2;
    let etg_f = crate::topopebrep::tool::topology::fun_tool_etg_f(
        paredge,
        &edge,
        &GpPnt2d::new(oo_uv.x(), oo_uv.y()),
        f,
        tola,
    );
    let mut rest = super::vpr::fun_in_los(&edge, erl);
    let isse = hds.ds().is_section_edge(&edge, true);
    rest = rest || isse;
    let mut interf2d = etg_f && has_oo_edge && rest;

    let factor = 1.0e-4;
    let mut mkt = MakeTransition::new();
    if !mkt.initialize(&edge, par1, par2, paredge, f, &oo_uv, factor) {
        return Ok(t);
    }
    let is_t2d = mkt.is_t2d();
    interf2d = interf2d && is_t2d;
    if interf2d && !mkt.set_rest(
        ooedge
            .as_ref()
            .ok_or_else(|| "GetEdgeTrans : no OO edge".to_string())?,
        oo_paredge,
    ) {
        return Ok(t);
    }
    let _ = facint;
    let _ = line_index;

    let mut stb = State::Unknown;
    let mut sta = State::Unknown;
    if !mkt.mk_ton_e(&mut stb, &mut sta) {
        return Ok(t);
    }
    t.set_state_before(stb);
    t.set_state_after(sta);
    Ok(t)
}

impl super::FacesFiller {
    /// `ProcessVPonclosingR(VP, GFace, ShapeIndex, transEdge, PVKind, PVIndex,
    /// EPIfound, IEPI)` (`vprclo.cxx:190-296`).
    #[allow(clippy::too_many_arguments)]
    pub fn process_vp_on_closing_r(
        &mut self,
        facint: &mut FacesIntersector,
        hds: &mut HDataStructure,
        vp: &VPointInter,
        shape_index: i32,
        trans_edge: &Transition,
        pv_kind: Kind,
        pv_index: i32,
    ) -> Result<(), String> {
        let isvertex = pv_kind == Kind::Vertex;
        let absindex = vp.shape_index();
        let oo_shape_index = if shape_index == 1 { 2 } else { 1 };
        let on2edges = absindex == 3;
        let has_on_edge = vp.state(oo_shape_index)? == State::On;
        let has_oo_edge = if on2edges { true } else { has_on_edge };

        let face = self.face(shape_index)?.clone();
        let ooface = self.face(oo_shape_index)?.clone();
        let mut i_oo_face = hds.ds().shape_index(&ooface, true);
        if i_oo_face == 0 {
            i_oo_face = hds.change_ds().add_shape_with_ancestor(&ooface, oo_shape_index);
        }

        let edge = vp
            .edge(shape_index)
            .ok_or_else(|| "ProcessVPonclosingR : no edge".to_string())?
            .clone();
        if !hds.ds().has_shape(&edge, true) {
            hds.change_ds().add_shape_with_ancestor(&edge, shape_index);
        }
        let paredge = vp.edge_parameter(shape_index);

        let mut oo_edge_index = 0;
        let mut oo_edge = TopoShape::new(ShapeType::Edge);
        if has_oo_edge {
            oo_edge = if on2edges {
                vp.edge(oo_shape_index)
                    .ok_or_else(|| "ProcessVPonclosingR : no OO edge".to_string())?
                    .clone()
            } else {
                vp.edge_on(oo_shape_index)?.clone()
            };
            oo_edge_index = if hds.ds().has_shape(&oo_edge, true) {
                hds.ds().shape_index(&oo_edge, true)
            } else {
                hds.change_ds().add_shape_with_ancestor(&oo_edge, oo_shape_index)
            };
        }

        let tunk = trans_edge.is_unknown();
        let trans_add = if tunk {
            let f1 = self
                .f1
                .clone()
                .ok_or_else(|| "ProcessVPonclosingR : no F1".to_string())?;
            let f2 = self
                .f2
                .clone()
                .ok_or_else(|| "ProcessVPonclosingR : no F2".to_string())?;
            get_edge_trans(
                hds,
                facint,
                self.line_index,
                &self.erl,
                vp,
                pv_kind,
                pv_index,
                shape_index,
                &ooface,
                &f1.0,
                &f2.0,
            )?
        } else {
            trans_edge.clone()
        };

        // `cxx:265-284`: one EPI on `<edge>` with the face as support, plus one
        // on `<OOedge>` when the VP lies on it.
        {
            let mut t = trans_add.clone();
            t.set_index(i_oo_face);
            let epi = match crate::topopebrep::ds::process_interferences_tool::make_epv_interference(
                &t,
                i_oo_face,
                pv_index,
                paredge,
                pv_kind,
                isvertex,
            ) {
                Some(i) => i,
                None => return Ok(()),
            };
            hds.store_interference_on_shape(&epi, &edge)?;
        }
        if has_oo_edge {
            let mut t = trans_add;
            t.set_index(i_oo_face);
            let epi = match crate::topopebrep::ds::process_interferences_tool::make_epv_interference(
                &t,
                oo_edge_index,
                pv_index,
                paredge,
                pv_kind,
                isvertex,
            ) {
                Some(i) => i,
                None => return Ok(()),
            };
            hds.store_interference_on_shape(&epi, &edge)?;
        }
        Ok(())
    }
}
