//! `FUN_VPIndex` — resolving (and when needed creating) the DS geometry a
//! VPoint stands for.
//!
//! Source (OCCT 8.0.0, `ModelingAlgorithms/TKBool/TopOpeBRep`):
//! - `FUN_VPIndex(FF, L, VP, ShapeIndex, HDS, DSCIL, PVKind, PVIndex,
//!   EPIfound, IEPI, CPIfound, ICPI, mkVP)` (`TopOpeBRep_vpr.cxx:350-469`,
//!   duplicated verbatim at `TopOpeBRep_vprdeg.cxx:90-...`)
//!
//! The `mkVP` mode constants are `M_FINDVP` / `M_MKNEWVP` / `M_GETVP`
//! (`vpr.cxx:346-348`); only the first two are reachable from the ported call
//! graph today (`ProcessVPondgE` at `vprdeg.cxx:643` and `:710`).

use crate::abs::State;
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::ds::interference::InterferenceHandle;
use crate::topopebrep::ds::InterferenceList;
use crate::topopebrep::ds::kind::Kind;
use crate::topopebrep::faces_filler::vpr::fun_vp_geometry_found;
use crate::topopebrep::faces_filler::FacesFiller;
use crate::topopebrep::faces_intersector::FacesIntersector;
use crate::topopebrep::interf::v_point_inter::VPointInter;

/// `M_FINDVP` (`vpr.cxx:346`): only look for an existing geometry.
pub const M_FINDVP: i32 = 0;
/// `M_MKNEWVP` (`vpr.cxx:347`): only make a new geometry.
pub const M_MKNEWVP: i32 = 1;
/// `M_GETVP` (`vpr.cxx:348`): look first, then make.
pub const M_GETVP: i32 = 2;

/// The out-parameters of `FUN_VPIndex` (`vpr.cxx:356-365`).
pub struct VpIndex {
    /// `PVKind`: `POINT` or `VERTEX`.
    pub pv_kind: Kind,
    /// `PVIndex`: the DS geometry index, `0` when `M_FINDVP` found nothing.
    pub pv_index: i32,
    /// `EPIfound` / `IEPI`: an interference on `<edge>`.
    pub epi_found: bool,
    pub iepi: Option<InterferenceHandle>,
    /// `CPIfound` / `ICPI`: an interference on `<DSCIL>`.
    pub cpi_found: bool,
    pub icpi: Option<InterferenceHandle>,
}

impl Default for VpIndex {
    fn default() -> Self {
        Self {
            pv_kind: Kind::Point,
            pv_index: 0,
            epi_found: false,
            iepi: None,
            cpi_found: false,
            icpi: None,
        }
    }
}

/// `FUN_VPIndex(FF, L, VP, ShapeIndex, HDS, DSCIL, ...)` (`vpr.cxx:350-469`).
///
/// `L` reaches [`fun_vp_geometry_found`] as the `FacesIntersector` plus the
/// current line index, which is how the port carries the line.
#[allow(clippy::too_many_arguments)]
pub fn fun_vp_index(
    ff: &FacesFiller,
    facint: &mut FacesIntersector,
    line_index: i32,
    vp: &VPointInter,
    shape_index: i32,
    hds: &mut HDataStructure,
    dscil: &InterferenceList,
    mk_vp: i32,
) -> Result<VpIndex, String> {
    let mut out = VpIndex::default();
    let oo_shape_index = if shape_index == 1 { 2 } else { 1 };
    let si_is_vertex = vp.is_vertex(shape_index)?;
    let oo_is_vertex = vp.is_vertex(oo_shape_index)?;

    let mut ooepi_found = false;
    if mk_vp == M_FINDVP || mk_vp == M_GETVP {
        let g = fun_vp_geometry_found(hds, facint, line_index, vp, shape_index, dscil)?;
        out.pv_kind = g.kind;
        out.pv_index = g.index;
        out.epi_found = g.epi_found;
        out.iepi = g.iepi.clone();
        out.cpi_found = g.cpi_found;
        out.icpi = g.icpi.clone();
        ooepi_found = g.ooepi_found;

        if mk_vp == M_FINDVP {
            // `vpr.cxx:395-407`: `OOEPIfound` is not reported to the caller, so
            // a hit that only the other edge saw must not leave `PVIndex` set.
            if !out.epi_found && !out.cpi_found && !ooepi_found {
                out.pv_index = 0;
            }
            return Ok(out);
        }
    }

    // `Gfound = EPIfound || CPIfound || OOEPIfound` (`cxx:409-411`).
    let gfound = out.epi_found || out.cpi_found || ooepi_found;
    let on2edges = vp.shape_index() == 3;
    let has_oo_edge = if on2edges {
        true
    } else {
        vp.state(oo_shape_index)? == State::On
    };

    if out.pv_index == 0 && (si_is_vertex || oo_is_vertex) {
        out.pv_kind = Kind::Vertex;
    } else if out.pv_index == 0 {
        out.pv_kind = Kind::Point;
    }

    if has_oo_edge && !gfound {
        if !ooepi_found {
            let si = if si_is_vertex {
                shape_index
            } else if oo_is_vertex {
                oo_shape_index
            } else {
                shape_index
            };
            let (i, k) = FacesFiller::make_geometry(hds.change_ds(), vp, si)?;
            out.pv_index = i;
            out.pv_kind = k;
        }
    }
    if !has_oo_edge && !gfound {
        let found = ff.get_ff_geometry(hds.ds(), vp)?.is_some();
        if !found {
            let si = if si_is_vertex {
                shape_index
            } else if oo_is_vertex {
                oo_shape_index
            } else {
                shape_index
            };
            let (i, k) = FacesFiller::make_geometry(hds.change_ds(), vp, si)?;
            out.pv_index = i;
            out.pv_kind = k;
        }
    }
    Ok(out)
}
