use crate::shape::TopoShape;

use super::super::data_structure::DataStructure;
use super::super::h_data_structure::HDataStructure;

use super::*;

/// `FUN_ds_hasFEI(pDS2d, F, GI, ITRA)` (`EXPORT.cxx:3279-3306`): `pDS2d`
/// already carries an interference of `F` whose geometry is `GI` and whose
/// transition index is `ITRA`.
pub fn fun_ds_has_fei(p_ds2d: &DataStructure, f: &TopoShape, gi: i32, itra: i32) -> bool {
    if !p_ds2d.has_shape(f, true) {
        return false;
    }
    for i in p_ds2d.shape_interferences(f, true) {
        let Ok(guard) = i.read() else {
            continue;
        };
        let (_gt, g, _st, _s) = fds_data(&guard);
        if g == gi && guard.transition().index().unwrap_or(-1) == itra {
            return true;
        }
    }
    false
}

/// `FUN_ds_getVsdm(BDS, iV, iVsdm)` (`cxx:226-252`): the index of the first
/// shape same-domain with `BDS.Shape(iV)` that is not the shape itself.
pub fn fun_ds_get_vsdm(bds: &DataStructure, iv: i32) -> Option<i32> {
    let imax = bds.nb_shapes();
    if iv < 1 || iv > imax {
        return None;
    }
    let v = bds.shape(iv, true)?;
    for vv in bds.shape_same_domain(&v) {
        if crate::topo_tools::is_same(&v, vv) {
            continue;
        }
        return Some(bds.shape_index(vv, true));
    }
    None
}

/// `FUN_ds_sdm(BDS, s1, s2)` (`cxx:255-274`): `s2` is in `s1`'s same-domain
/// list. `HasShape` defaults to `FindKeep = true`.
pub fn fun_ds_sdm(bds: &DataStructure, s1: &TopoShape, s2: &TopoShape) -> bool {
    if !bds.has_shape(s1, true) || !bds.has_shape(s2, true) {
        return false;
    }
    bds.shape_same_domain(s1)
        .iter()
        .any(|s| crate::topo_tools::is_same(s, s2))
}

/// `FUN_ds_getoov(v, BDS, oov)` (`cxx:647-668`): the "other" vertex of a
/// same-domain vertex pair (the DS binds at most two).
pub fn fun_ds_getoov(bds: &DataStructure, v: &TopoShape) -> Option<TopoShape> {
    for vcur in bds.shape_same_domain(v) {
        if crate::topo_tools::is_same(vcur, v) {
            continue;
        }
        return Some(vcur.clone());
    }
    None
}

/// `FUN_ds_getoov(v, HDS, oov)` (`cxx:671-694`): the same on the HDS, which
/// first checks `HasSameDomain(v)`.
pub fn fun_ds_getoov_in_hds(hds: &HDataStructure, v: &TopoShape) -> Option<TopoShape> {
    let hsd = hds.has_same_domain(v, true);
    if !hsd {
        return None;
    }
    fun_ds_getoov(hds.ds(), v)
}
