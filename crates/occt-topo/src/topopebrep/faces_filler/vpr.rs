//! `TopOpeBRep_vpr.cxx` — the VPoint-on-restriction pass and its geometry
//! lookup helpers.
//!
//! Sources (OCCT 8.0.0, `ModelingAlgorithms/TKBool/TopOpeBRep`):
//! - `FUN_GetGonParameter` (`TopOpeBRep_vpr.cxx:120-149`)
//! - `FUN_INlos` (`vpr.cxx:151-162`)
//! - `FUN_transForWL` (`vpr.cxx:178-217`)
//! - `FUN_VPgeometryfound` (`vpr.cxx:220-345`)
//! - `FUN_VPIndex` (`vpr.cxx:351-470`)
//! - `ProcessVPnotonR` (`TopOpeBRep_FacesFiller_1.cxx:507-577`)
//! - `ProcessVPR` (`_1.cxx:581-645`)
//!
//! UNPORTED, needed by `ProcessVPonR` / `ProcessRLine`:
//! - `ProcessVPonR` (`vpr.cxx:877-1732`) and the helpers it drives,
//!   `FUN_newtransEdge` (`vpr.cxx:515-876`), `FUN_LineRestF`
//!   (`vpr.cxx:472-513`), `FUN_processCPI` (`vpr.cxx:1741-...`),
//! - `FUN_tool_tggeomE` / `FUN_tool_paronEF` / `FUN_tool_onapex`
//!   (`TopOpeBRepTool_TOPOLOGY.cxx`),
//! - `FDS_Parameter` / `FDS_data` / `FDS_LOIinfsup` (`TopOpeBRepDS_EXPORT.cxx`).

use occt_core::precision::PCONFUSION;

use crate::abs::{Orientation, State};
use crate::shape::TopoShape;
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::ds::interference::InterferenceHandle;
use crate::topopebrep::ds::interference_tool::make_curve_interference;
use crate::topopebrep::ds::kind::Kind;
use crate::topopebrep::ds::point_iterator::PointIterator;
use crate::topopebrep::ds::transition::Transition;
use crate::topopebrep::ds::InterferenceList;
use crate::topopebrep::faces_intersector::FacesIntersector;
use crate::topopebrep::ff_transition_tool::FFTransitionTool;
use crate::topopebrep::interf::point_geom_tool::PointGeomTool;
use crate::topopebrep::interf::v_point_inter::VPointInter;
use crate::topopebrep::line_inter::LineInter;
use crate::topopebrep::type_line_curve::TypeLineCurve;

/// `M_FINDVP` / `M_MKNEWVP` / `M_GETVP` (`vpr.cxx:346-348`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MkVpMode {
    /// `M_FINDVP`: only look for a new VP.
    Find,
    /// `M_MKNEWVP`: only make a new VP.
    Make,
    /// `M_GETVP`: `Find`, then `Make` when the lookup fails.
    Get,
}

/// `FUN_INlos(S, loS)` (`vpr.cxx:151-162`): `S` is `IsSame` as one of `loS`.
pub fn fun_in_los(s: &TopoShape, lo_s: &[TopoShape]) -> bool {
    lo_s.iter().any(|e| crate::topo_tools::is_same(e, s))
}

/// `FUN_GetGonParameter(IT, par, tolp, G, GT)` (`vpr.cxx:120-149`): the first
/// interference of `it` carrying a parameter within `tolp` of `par`.
///
/// `it` is left on the hit, as in OCCT. `FDS_Parameter` fails for the
/// interference types without a parameter, which the port reads from
/// `PointIterator::parameter` (its `Err` arm).
pub fn fun_get_g_on_parameter(
    it: &mut PointIterator,
    par: f64,
    tolp: f64,
) -> Option<(i32, Kind)> {
    while it.more() {
        let ipar = match it.parameter() {
            Ok(p) => p,
            Err(_) => {
                it.next();
                continue;
            }
        };
        if (par - ipar).abs() >= tolp {
            it.next();
            continue;
        }
        let i = it.current();
        return Some((i.geometry(), i.geometry_type()));
    }
    None
}

/// `FUN_transForWL(L, iVP, ShapeIndex, transLine)` (`vpr.cxx:178-217`): the
/// complement of the first determinate transition after `iVP`.
pub fn fun_trans_for_wl(
    l: &LineInter,
    i_vp: i32,
    shape_index: i32,
) -> Result<Option<Transition>, String> {
    let nb = l.nb_v_point();
    for k in 1..=nb {
        let Some(vpbis) = l.v_point(k) else {
            continue;
        };
        if !vpbis.keep() {
            continue;
        }
        let i_vpbis = vpbis.index();
        if i_vpbis <= i_vp {
            continue;
        }
        let absindexbis = vpbis.shape_index();
        let shapeindexbis = if absindexbis == 3 {
            shape_index
        } else {
            absindexbis
        };
        if shapeindexbis == 0 {
            continue;
        }
        let Some(edgebis) = vpbis.edge(shapeindexbis) else {
            continue;
        };
        let edgeoribis = edgebis.orientation();
        let trans_linebis =
            FFTransitionTool::process_line_transition_with_index(vpbis, shapeindexbis, edgeoribis);
        if trans_linebis.is_unknown() {
            continue;
        }
        return Ok(Some(trans_linebis.complement()));
    }
    Ok(None)
}

/// `FUN_VPgeometryfound`'s out-parameters (`vpr.cxx:220-345`).
pub struct VpGeometry {
    /// `PVIndex` / `PVKind`: the geometry the VP was found at.
    pub index: i32,
    pub kind: Kind,
    /// `EPIfound` / `IEPI`: an interference on `<edge>`.
    pub epi_found: bool,
    pub iepi: Option<InterferenceHandle>,
    /// `CPIfound` / `ICPI`: an interference on `<DSCIL>`.
    pub cpi_found: bool,
    pub icpi: Option<InterferenceHandle>,
    /// `OOEPIfound` / `IOOEPI`: an interference on `<OOedge>`; only meaningful
    /// when the VP is on both edges.
    pub ooepi_found: bool,
    pub iooepi: Option<InterferenceHandle>,
}

impl Default for VpGeometry {
    fn default() -> Self {
        Self {
            index: 0,
            kind: Kind::Point,
            epi_found: false,
            iepi: None,
            cpi_found: false,
            icpi: None,
            ooepi_found: false,
            iooepi: None,
        }
    }
}

impl VpGeometry {
    /// `Gfound = EPIfound || CPIfound || OOEPIfound` (`vpr.cxx:405-407`).
    pub fn g_found(&self) -> bool {
        self.epi_found || self.cpi_found || self.ooepi_found
    }
}

/// `FUN_VPgeometryfound(FF, L, VP, ShapeIndex, HDS, DSCIL, ...)`
/// (`vpr.cxx:220-345`).
///
/// `FF.GetGeometry(it, VP, ...)` is `HDataStructure::get_geometry` on the
/// interference list of the shape, which is what `BDS.ShapeInterferences`
/// returns.
pub fn fun_vp_geometry_found(
    hds: &HDataStructure,
    facint: &mut FacesIntersector,
    line_index: i32,
    vp: &VPointInter,
    shape_index: i32,
    dscil: &InterferenceList,
) -> Result<VpGeometry, String> {
    let mut out = VpGeometry::default();
    let dsp = PointGeomTool::make_point_of_v_point_inter(vp)?;

    let lrest = {
        let l = facint.change_line(line_index);
        l.type_line_curve() == TypeLineCurve::Restriction
    };
    let (erest, par_erest, rk_erest) = if lrest {
        let l = facint.change_line(line_index);
        let e = l
            .arc()
            .ok_or_else(|| "FUN_VPgeometryfound : no restriction arc".to_string())?
            .clone();
        let par = vp.parameter_on_line()?;
        let isedge1 = l.arc_is_edge(1);
        let isedge2 = l.arc_is_edge(2);
        let rk = if isedge1 {
            1
        } else if isedge2 {
            2
        } else {
            0
        };
        (Some(e), par, rk)
    } else {
        (None, 0.0, 0)
    };

    let absindex = vp.shape_index();
    let ooabsindex = if absindex == 1 { 2 } else { 1 };
    let oo_shape_index = if shape_index == 1 { 2 } else { 1 };
    let on2edges = (absindex == 3) || (lrest && rk_erest == ooabsindex);
    let edge = if rk_erest == shape_index {
        erest
            .clone()
            .ok_or_else(|| "FUN_VPgeometryfound : no restriction arc".to_string())?
    } else {
        vp.edge(shape_index)
            .ok_or_else(|| "FUN_VPgeometryfound : no VP edge".to_string())?
            .clone()
    };

    let par = if rk_erest == shape_index {
        par_erest
    } else {
        vp.edge_parameter(shape_index)
    };
    let tole = crate::topopebrep::tool::topology::fun_tool_max_tol(&edge);
    let tolp = tole * PCONFUSION;

    if hds.ds().has_shape(&edge, true) {
        let epil = hds.ds().shape_interferences(&edge, true).to_vec();
        let mut it = PointIterator::new(epil.clone());
        let hit = hds.get_geometry(&mut it, &dsp)?;
        match hit {
            Some((g, k)) => {
                out.epi_found = true;
                out.index = g;
                out.kind = k;
                out.iepi = Some(it.value().clone());
            }
            None => {
                let mut it = PointIterator::new(epil);
                if let Some((g, k)) = fun_get_g_on_parameter(&mut it, par, tolp) {
                    out.epi_found = true;
                    out.index = g;
                    out.kind = k;
                    out.iepi = Some(it.value().clone());
                }
            }
        }
    }

    let mut it = PointIterator::new(dscil.clone());
    if let Some((g, k)) = hds.get_geometry(&mut it, &dsp)? {
        out.cpi_found = true;
        out.index = g;
        out.kind = k;
        out.icpi = Some(it.value().clone());
    }

    let has_oo_edge = if on2edges {
        true
    } else {
        vp.state(oo_shape_index)? == State::On
    };
    if has_oo_edge {
        let ooedge = if on2edges {
            if rk_erest == oo_shape_index {
                erest
                    .clone()
                    .ok_or_else(|| "FUN_VPgeometryfound : no restriction arc".to_string())?
            } else {
                vp.edge(oo_shape_index)
                    .ok_or_else(|| "FUN_VPgeometryfound : no OO VP edge".to_string())?
                    .clone()
            }
        } else {
            vp.edge_on(oo_shape_index)?.clone()
        };

        let oopar = if on2edges {
            if rk_erest == oo_shape_index {
                par_erest
            } else {
                vp.edge_parameter(oo_shape_index)
            }
        } else {
            vp.edge_on_parameter(oo_shape_index)?
        };

        let tol_ooe = crate::topopebrep::tool::topology::fun_tool_max_tol(&ooedge);
        let oo_tolp = tol_ooe * PCONFUSION;
        if hds.ds().has_shape(&ooedge, true) {
            let ooepil = hds.ds().shape_interferences(&ooedge, true).to_vec();
            let mut ooit = PointIterator::new(ooepil.clone());
            let hit = hds.get_geometry(&mut ooit, &dsp)?;
            match hit {
                Some((g, k)) => {
                    out.ooepi_found = true;
                    out.index = g;
                    out.kind = k;
                    out.iooepi = Some(ooit.value().clone());
                }
                None => {
                    let mut ooit = PointIterator::new(ooepil);
                    if let Some((g, k)) = fun_get_g_on_parameter(&mut ooit, oopar, oo_tolp) {
                        out.ooepi_found = true;
                        out.index = g;
                        out.kind = k;
                        out.iooepi = Some(ooit.value().clone());
                    }
                }
            }
        }
    }

    Ok(out)
}

impl super::FacesFiller {
    /// `ProcessVPnotonR(VP)` (`_1.cxx:507-577`).
    pub fn process_vp_not_on_r(
        &mut self,
        facint: &mut FacesIntersector,
        hds: &mut HDataStructure,
        vp: &VPointInter,
    ) -> Result<(), String> {
        let shape_index = 0;
        let i_vp = vp.index();

        let (i_inon1, i_inonn, _n_inon) = self.line_vp_bounds(facint);
        let mut it = PointIterator::new(self.dscil.clone());

        // `_1.cxx:517-536`: `CPIfound = GetGeometry(itCPIL, VP, ...)`; when it
        // fails the VP is rejected unless it sits at one of the line bounds.
        let cpi_found = self.get_geometry(hds, &mut it, vp)?.is_some();
        if !cpi_found && i_vp != i_inon1 && i_vp != i_inonn {
            return Ok(());
        }

        let (pv_index, pv_kind) = if cpi_found {
            let (g, k) = {
                let i = it.current();
                (i.geometry(), i.geometry_type())
            };
            (g, k)
        } else {
            match self.get_ff_geometry(hds.ds(), vp)? {
                Some((g, k)) => (g, k),
                None => Self::make_geometry(hds.change_ds(), vp, shape_index)?,
            }
        };

        let trans_line = if cpi_found {
            it.current().transition().complement()
        } else if i_vp == i_inon1 {
            Transition::from_orientation(Orientation::Forward)
        } else if i_vp == i_inonn {
            Transition::from_orientation(Orientation::Reversed)
        } else {
            Transition::default()
        };

        let parline = vp.parameter_on_line()?;
        let cpi = make_curve_interference(&trans_line, Kind::Curve, 0, pv_kind, pv_index, parline);
        self.store_curve_interference(hds, &cpi);
        Ok(())
    }

    /// `ProcessVPR(FF, VP)` (`_1.cxx:581-645`): dispatch on the VP shape index.
    pub fn process_vpr(
        &mut self,
        facint: &mut FacesIntersector,
        hds: &mut HDataStructure,
        vp: &VPointInter,
    ) -> Result<(), String> {
        if !vp.keep() {
            return Ok(());
        }
        let shape_index = vp.shape_index();
        let line_t_on_f1 = self.face_face_transition(facint, 1)?.clone();
        let line_t_on_f2 = self.face_face_transition(facint, 2)?.clone();
        let f1 = self
            .f1
            .clone()
            .ok_or_else(|| "TopOpeBRep_FacesFiller::ProcessVPR : no F1".to_string())?;
        let f2 = self
            .f2
            .clone()
            .ok_or_else(|| "TopOpeBRep_FacesFiller::ProcessVPR : no F2".to_string())?;
        match shape_index {
            0 => self.process_vp_not_on_r(facint, hds, vp)?,
            1 => self.process_vp_on_r(facint, hds, vp, &line_t_on_f1, &f1.0, 1)?,
            2 => self.process_vp_on_r(facint, hds, vp, &line_t_on_f2, &f2.0, 2)?,
            3 => {
                let is_v1 = vp.is_vertex_on_s1()?;
                let is_v2 = vp.is_vertex_on_s2()?;
                let shin1 = if is_v2 && !is_v1 { 2 } else { 1 };
                if shin1 == 1 {
                    self.process_vp_on_r(facint, hds, vp, &line_t_on_f1, &f1.0, 1)?;
                    self.process_vp_on_r(facint, hds, vp, &line_t_on_f2, &f2.0, 2)?;
                } else {
                    self.process_vp_on_r(facint, hds, vp, &line_t_on_f2, &f2.0, 2)?;
                    self.process_vp_on_r(facint, hds, vp, &line_t_on_f1, &f1.0, 1)?;
                }
            }
            _ => {}
        }
        Ok(())
    }
}
