use crate::shape::TopoShape;

use super::super::h_data_structure::HDataStructure;

use super::*;

/// `FUN_ds_FillSDMFaces(HDS)` (`cxx:1874-1947`): every edge of a face that has
/// a same-domain face of the other ancestor rank becomes a section edge when
/// [`fun_eis_se2`] accepts one of its interferences.
pub fn fun_ds_fill_sdm_faces(hds: &mut HDataStructure) -> Result<(), String> {
    let ns = hds.change_ds().nb_shapes();
    for i in 1..=ns {
        let Some(s) = hds.ds().shape(i, true) else {
            continue;
        };
        if s.shape_type() != crate::abs::ShapeType::Face {
            continue;
        }
        let rk_s = hds.ds().ancestor_rank(i);
        // `cxx:1888-1899`: the same-domain faces of another rank.
        let mut fsdm: Vec<TopoShape> = Vec::new();
        for f in hds.ds().shape_same_domain(&s) {
            let rkf = hds.ds().ancestor_rank_of(f);
            if rkf != rk_s {
                fsdm.push(f.clone());
            }
        }
        if fsdm.is_empty() {
            continue;
        }

        let mut edges: Vec<TopoShape> = Vec::new();
        let mut ex = crate::topexp::Explorer::new(&s, crate::abs::ShapeType::Edge);
        while ex.more() {
            edges.push(ex.current().clone());
            ex.next();
        }
        for e in edges {
            // `cxx:1908-1916`
            if !hds.ds().has_shape(&e, true) {
                continue;
            }
            if hds.ds().is_section_edge(&e, true) {
                continue;
            }
            let li = hds.ds().shape_interferences(&e, true).to_vec();
            if li.len() < 1 {
                continue;
            }
            // `cxx:1926-1943`
            for handle in &li {
                let itf = handle.read().expect("interference lock").clone();
                let is_se2 = fun_eis_se2(hds.ds(), &e, &itf);
                if is_se2 == EISSE2_NONE {
                    continue;
                }
                hds.change_ds().add_section_edge(&e);
                if is_se2 == EISSE2_EREF_AND_ESUP {
                    let st = itf.support();
                    if let Some(es) = hds.ds().shape(st, true) {
                        hds.change_ds().add_section_edge(&es);
                    }
                }
                break;
            }
        }
    }
    Ok(())
}

/// `FUN_ds_addSEsdm1d(HDS)` (`cxx:1949-1995`): every non-degenerated, non
/// section edge of rank 1 that shares geometry with a same-domain edge becomes
/// a section edge along with its `EShareG` partners.
pub fn fun_ds_add_sesdm1d(hds: &mut HDataStructure) -> Result<(), String> {
    let ns = hds.change_ds().nb_shapes();
    for i in 1..=ns {
        let Some(s) = hds.ds().shape(i, true) else {
            continue;
        };
        if s.shape_type() != crate::abs::ShapeType::Edge {
            continue;
        }
        // `cxx:1961-1965`
        if crate::brep_tool::BRepTool::is_degenerated(&crate::shape::Edge(s.clone())) {
            continue;
        }
        // `cxx:1967-1971`
        if hds.ds().is_section_edge(&s, true) {
            continue;
        }
        // `cxx:1972-1976`
        if hds.ds().ancestor_rank(i) != 1 {
            continue;
        }
        let mut lsd: Vec<TopoShape> = Vec::new();
        super::super::tool::e_share_g(hds, &s, &mut lsd);
        let share_g = !lsd.is_empty();
        for e in &lsd {
            hds.change_ds().add_section_edge(e);
        }
        if share_g {
            hds.change_ds().add_section_edge(&s);
        }
    }
    Ok(())
}
