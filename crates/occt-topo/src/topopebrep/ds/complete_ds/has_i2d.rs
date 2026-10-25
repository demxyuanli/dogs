//! `FUN_ds_hasI2d`.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_EXPORT.cxx:1997-2053`.

use crate::abs::ShapeType;
use crate::topopebrep::ds::export::{fds_data, fds_tdata};
use crate::topopebrep::ds::interference_iterator::InterferenceList;
use crate::topopebrep::ds::kind::Kind;
use crate::topopebrep::ds::process_interferences_tool::{
    fun_select_sk_interference, fun_select_trasha_interference,
};

/// `FUN_ds_hasI2d(EIX, LLI, LI2d)` (`cxx:1997-2053`).
///
/// `LLI` is attached to `EIX` at a given `G`. An interference is `2d` when no
/// `(T(FTRA), G, FTRA)` counterpart exists in the list; answers the `2d` ones.
///
/// The first parameter `EIX` of OCCT is unused (`cxx:1998-1999`).
pub fn fun_ds_has_i2d(lli: &InterferenceList) -> InterferenceList {
    // `cxx:2008-2017`
    let mut li = lli.to_vec();
    let mut l1de: InterferenceList = Vec::new();
    fun_select_trasha_interference(&mut li, ShapeType::Edge, &mut l1de);
    let mut lif: InterferenceList = Vec::new();
    fun_select_sk_interference(&mut li, Kind::Face, &mut lif);
    let mut lie: InterferenceList = Vec::new();
    fun_select_sk_interference(&mut li, Kind::Edge, &mut lie);
    li.extend(l1de);

    // `cxx:2018-2048`
    let mut li2d: InterferenceList = Vec::new();
    for ie in &lie {
        let (ge, isb) = {
            let intf = ie.read().expect("interference lock");
            let (_gte, ge, _ste, _se) = fds_data(&intf);
            let (_tsb, isb, _tsa, _isa) = fds_tdata(&intf.transition());
            (ge, isb)
        };
        let mut is3d = false;
        for iff in &lif {
            let (gtf, gf, sf) = {
                let intf = iff.read().expect("interference lock");
                let (gtf, gf, _stf, sf) = fds_data(&intf);
                (gtf, gf, sf)
            };
            let _ = gtf;            if ge != gf {
                continue;
            }
            if sf == isb {
                is3d = true;
                break;
            }
        }
        if !is3d {
            li2d.push(ie.clone());
        }
    }
    li2d
}
