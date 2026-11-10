use crate::shape::TopoShape;

use super::super::data_structure::DataStructure;
use super::super::interference::Interference;
use super::super::kind::Kind;

use super::*;

/// `FUN_EisSE2` return codes (`cxx:1792-1795`).
pub const EISSE2_NONE: i32 = 0;
/// `EREFandESUP`
pub const EISSE2_EREF_AND_ESUP: i32 = 2;

/// `FUN_EisSE2(BDS, Fsdm, E, I)` (`cxx:1797-1871`).
///
/// `Fsdm` is accepted but unused by OCCT (`cxx:1800`).
pub fn fun_eis_se2(bds: &DataStructure, e: &TopoShape, i: &Interference) -> i32 {
    let (sb, _ib, _sa, _ia, gt, g, _st, _s) = fds_idata(i);
    if gt != Kind::Vertex {
        return EISSE2_NONE;
    }
    // `cxx:1811-1818`: the orientation of the transition's state before, and
    // no `EXTERNAL` interference.
    let o = i
        .transition()
        .orientation(crate::abs::State::In, crate::abs::ShapeType::Shape);
    if o == crate::abs::Orientation::External {
        return EISSE2_NONE;
    }
    // `cxx:1820-1825`: E same-domain edge
    if sb != crate::abs::ShapeType::Edge {
        return EISSE2_NONE;
    }

    let e_index = bds.shape_index(e, true);
    let rk_e = bds.ancestor_rank(e_index);
    let rk_g = bds.ancestor_rank(g);
    let Some(vg) = bds.shape(g, true) else {
        return EISSE2_NONE;
    };
    let gsd = fun_ds_get_vsdm(bds, g);
    let o_g_in_e = if let Some(gsd) = gsd {
        let vgsd = bds.shape(gsd, true).unwrap_or_else(|| vg.clone());
        if rk_e == rk_g {
            crate::topopebrep::tool::fun_tool_orient_vin_e(&vg, e)
        } else {
            crate::topopebrep::tool::fun_tool_orient_vin_e(&vgsd, e)
        }
    } else {
        crate::topopebrep::tool::fun_tool_orient_vin_e(&vg, e)
    };

    // `cxx:1850-1853`
    if o_g_in_e == 0 {
        return EISSE2_EREF_AND_ESUP; // G in IN E
    }
    // `cxx:1855-1861`
    let noshare = (o_g_in_e == crate::topopebrep::tool::LAST
        && o == crate::abs::Orientation::Forward)
        || (o_g_in_e == crate::topopebrep::tool::FIRST
            && o == crate::abs::Orientation::Reversed);
    if noshare {
        return EISSE2_NONE;
    }
    // `cxx:1863-1869`: E same-domain edge(ST); both arms answer
    // `EREFandESUP`.
    EISSE2_EREF_AND_ESUP
}
