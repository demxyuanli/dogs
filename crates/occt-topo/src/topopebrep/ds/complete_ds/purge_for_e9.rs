//! `FUN_ds_PURGEforE9` and its `FUN_ds_hasSDMancestorfaces` helper.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_EXPORT.cxx:715-859`.

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, TopoShape};
use crate::topopebrep::ds::connex::Connexity;
use crate::topopebrep::ds::export::{fds_data, fds_si_is_gi_of_i_of_sba_of_t_of_i};
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::ds::interference_iterator::InterferenceList;
use crate::topopebrep::ds::kind::Kind;
use crate::topopebrep::ds::process_interferences_tool::fun_select_sk_interference;

/// `FUN_ds_hasSDMancestorfaces(HDS, E1, E2, Fsd1, Fsd2)` (`cxx:715-758`):
/// a same-domain face pair, one face an ancestor of `E1` and the other an
/// ancestor of `E2`.
pub fn fun_ds_has_sdm_ancestor_faces(
    hds: &HDataStructure,
    e1: &TopoShape,
    e2: &TopoShape,
    connex: &Connexity,
) -> Option<(TopoShape, TopoShape)> {
    let l_fce1 = connex.edge_connexity_same_shape(e1, hds);
    let l_fce2 = connex.edge_connexity_same_shape(e2, hds);

    for fce1 in l_fce1 {
        // `cxx:737-741`
        if !hds.has_same_domain(fce1, true) {
            continue;
        }
        // `cxx:742-755`
        for f2 in hds.ds().shape_same_domain(fce1) {
            if l_fce2.iter().any(|f| crate::topo_tools::is_same(f, f2)) {
                return Some((fce1.clone(), f2.clone()));
            }
        }
    }
    None
}

/// `FUN_ds_PURGEforE9(HDS)` (`cxx:761-859`).
///
/// xpu040998: an edge/face interference that describes a contact on geometry
/// where the edge belongs to the geometry of the face (`EPIe + EPIf = 3dI`)
/// is faulty and is removed.
pub fn fun_ds_purge_for_e9(hds: &mut HDataStructure, connex: &Connexity) {
    let ns = hds.ds().nb_shapes();
    for i in 1..=ns {
        // `cxx:768-783`
        let Some(ee) = hds.ds().shape(i, true) else {
            continue;
        };
        if ee.shape_type() != ShapeType::Edge {
            continue;
        }
        if BRepTool::is_degenerated(&Edge(ee.clone())) {
            continue;
        }

        // `cxx:785-792`
        let ie = hds.ds().shape_index(&ee, true);
        let mut li_copy: InterferenceList = hds.ds().shape_interferences(&ee, true).to_vec();
        let mut l3df: InterferenceList = Vec::new();
        let n3df = fun_select_sk_interference(&mut li_copy, Kind::Face, &mut l3df);
        if n3df == 0 {
            continue;
        }

        // `cxx:794-808`
        let lfcxe = connex.edge_connexity_same_shape(&ee, hds).to_vec();
        if lfcxe.is_empty() {
            // NYIRaise
            continue;
        }

        // `cxx:809-848`
        let mut removed = false;
        let mut it = 0usize;
        while it < l3df.len() {
            let fhas_ge = {
                let intf = l3df[it].read().expect("interference lock");
                fds_si_is_gi_of_i_of_sba_of_t_of_i(hds.ds(), ie, &intf)
            };
            if fhas_ge {
                // E has split ON F (cto904A3;e19,f14)
                removed = true;
                l3df.remove(it);
                continue;
            }

            let s = {
                let intf = l3df[it].read().expect("interference lock");
                let (_gt, _g, _st, s) = fds_data(&intf);
                s
            };
            let Some(f) = hds.ds().shape(s, true) else {
                it += 1;
                continue;
            };
            let hsdm = hds.has_same_domain(&f, true);
            if !hsdm {
                it += 1;
                continue;
            }
            let found_in_sdm = hds
                .ds()
                .shape_same_domain(&f)
                .iter()
                .any(|fsd| lfcxe.iter().any(|fc| crate::topo_tools::is_same(fc, fsd)));
            if !found_in_sdm {
                // E is IN geometry(F) (cto002D2;e33,f31)
                it += 1;
                continue;
            }

            // E is edge of FF sdm with F
            removed = true;
            l3df.remove(it);
        }

        // `cxx:850-858`
        if !removed {
            continue;
        }
        let lii = hds.change_ds().change_shape_interferences(&ee);
        lii.clear();
        lii.append(&mut li_copy);
        lii.append(&mut l3df);
    }
}
