//! `TopOpeBRepDS_repvg` — the vertex/point interference reducer on an edge.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_repvg.cxx:38-303` (OCCT 8.0.0).
//!
//! `FDS_repvg` splits the interferences of `LOI` by the face index they end on
//! (`Transition().Index()`, `cxx:257-275`) and runs [`fds_repvg2`] on every
//! group of two or more (`cxx:277-292`).
//!
//! `FDS_repvg2` walks the group pairwise: a pair of interferences on the same
//! support edge and a group of exactly two collapses into one 2d transition
//! (`FDS_stateEwithF2d`, `cxx:175-188`), while pairs on different support edges
//! are folded into a 3d complex transition with
//! [`EdgeInterferenceTool`] (`cxx:190-217`).
//!
//! Translation notes:
//! - `memeS` and `isComplex` are declared before the inner loop and read after
//!   it (`cxx:133-135`, `:222-238`), so the *last* value they took drives the
//!   outer decision; the port keeps the same scope.
//! - `LI.Remove(it2)` / `LI.Remove(it1)` leave the iterator on the next element
//!   (`NCollection_List`), which the port reproduces by not advancing the
//!   index after a removal.
//! - `BDS.Point(G1)` (`cxx:107`) is only handed to
//!   `EdgeInterferenceTool::Add(E, PDS, I)`, whose point argument is unused
//!   (`TopOpeBRepDS_EdgeInterferenceTool.cxx:130`), so it is not materialised.
//! - `NCollection_DataMap` iteration order in `FDS_repvg` is unspecified
//!   (`cxx:278`); the port keeps first-seen order, which is deterministic.

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, TopoShape};

use super::data_structure::DataStructure;
use super::edge_interference_tool::EdgeInterferenceTool;
use super::export::{
    fds_idata, fds_parameter_or_fail, fds_state_e_with_f2d, fun_ds_getoov,
};
use super::interference_iterator::InterferenceList;
use super::kind::Kind;
use super::transition::Transition;

/// `FDS_repvg2(BDS, EIX, GT, LI, RLI)` (`cxx:38-244`).
pub fn fds_repvg2(
    bds: &DataStructure,
    eix: i32,
    gt: Kind,
    li: &mut InterferenceList,
    rli: &mut InterferenceList,
) -> Result<(), String> {
    let e = bds
        .shape(eix, true)
        .ok_or_else(|| "FDS_repvg2: BDS.Shape(EIX)".to_string())?;
    if BRepTool::is_degenerated(&Edge(e.clone())) {
        return Ok(());
    }

    let is_point = gt == Kind::Point;
    let is_vertex = gt == Kind::Vertex;
    let n_li = li.len();

    let mut i1 = 0;
    while i1 < li.len() {
        let h1 = li[i1].clone();
        let (tsb1, isb1, tsa1, _isa1, gt1, mut g1, _st1, s1) = {
            let i = h1.read().expect("interference lock");
            fds_idata(&i)
        };
        if tsb1 != ShapeType::Face {
            i1 += 1;
            continue;
        }
        if tsa1 != ShapeType::Face {
            i1 += 1;
            continue;
        }
        if gt1 != gt {
            i1 += 1;
            continue;
        }

        // xpu121198 : raise in PRO16303, if G1 is a vertex same domain, make sure
        // rk(G1) = rk(E1)
        let rk_g1 = bds.ancestor_rank(g1);
        let rk_s1 = bds.ancestor_rank(s1);
        if is_vertex && rk_g1 != rk_s1 {
            let g1_shape = bds
                .shape(g1, true)
                .ok_or_else(|| "FDS_repvg2: BDS.Shape(G1)".to_string())?;
            match fun_ds_getoov(bds, &g1_shape) {
                Some(oov) => {
                    // !!NYIRAISE
                    g1 = bds.shape_index(&oov, true);
                }
                None => {
                    i1 += 1;
                    continue;
                }
            }
        }

        let e1 = bds
            .shape(s1, true)
            .ok_or_else(|| "FDS_repvg2: BDS.Shape(S1)".to_string())?;
        let f1 = bds
            .shape(isb1, true)
            .ok_or_else(|| "FDS_repvg2: BDS.Shape(isb1)".to_string())?;
        let vds: Option<TopoShape> = if is_point {
            None
        } else if is_vertex {
            Some(
                bds.shape(g1, true)
                    .ok_or_else(|| "FDS_repvg2: BDS.Shape(G1)".to_string())?,
            )
        } else {
            return Err("TopOpeBRepDS FDS_repvg2 1".to_string());
        };

        if BRepTool::is_degenerated(&Edge(e1.clone())) {
            i1 += 1;
            continue;
        }
        if i1 + 1 >= li.len() {
            i1 += 1;
            continue;
        }

        let mut ei_tool = EdgeInterferenceTool::new();
        let mut meme_s = false;
        let mut tr_meme_s = Transition::default();
        let mut is_complex = false;

        let mut i2 = i1 + 1;
        while i2 < li.len() {
            let h2 = li[i2].clone();
            let (tsb2, _isb2, tsa2, _isa2, _gt2, g2, _st2, s2) = {
                let i = h2.read().expect("interference lock");
                fds_idata(&i)
            };
            if tsb2 != ShapeType::Face {
                i2 += 1;
                continue;
            }
            if tsa2 != ShapeType::Face {
                i2 += 1;
                continue;
            }
            if g1 != g2 {
                i2 += 1;
                continue;
            }

            let e2 = bds
                .shape(s2, true)
                .ok_or_else(|| "FDS_repvg2: BDS.Shape(S2)".to_string())?;
            if BRepTool::is_degenerated(&Edge(e2.clone())) {
                i2 += 1;
                continue;
            }

            meme_s = s1 == s2;
            meme_s = meme_s && (n_li == 2);

            if !is_complex && meme_s {
                let p_e = {
                    let i = h1.read().expect("interference lock");
                    fds_parameter_or_fail(&i)?
                };
                let ok = fds_state_e_with_f2d(
                    bds,
                    &e,
                    p_e,
                    gt1,
                    g1,
                    &f1,
                    &mut tr_meme_s,
                )?;
                if !ok {
                    i2 += 1;
                    continue;
                }
                li.remove(i2);
                continue;
            } // !isComplex && memeS

            if !is_complex && !meme_s {
                is_complex = true;
                ei_tool.init(&e, &h1.read().expect("interference lock"))?;
                if is_point {
                    ei_tool.add_point(&e, &h1.read().expect("interference lock"))?;
                } else if is_vertex {
                    let v = vds
                        .as_ref()
                        .ok_or_else(|| "TopOpeBRepDS FDS_repvg2 1".to_string())?;
                    ei_tool.add_vertex(&e1, v, &h1.read().expect("interference lock"))?;
                }
            } // !isComplex && !memeS

            if is_complex && !meme_s {
                if is_point {
                    ei_tool.add_point(&e, &h2.read().expect("interference lock"))?;
                } else if is_vertex {
                    let v = vds
                        .as_ref()
                        .ok_or_else(|| "TopOpeBRepDS FDS_repvg2 1".to_string())?;
                    ei_tool.add_vertex(&e2, v, &h2.read().expect("interference lock"))?;
                }
                li.remove(i2);
                continue;
            } // (isComplex && !memeS)

            // (isComplex && memeS) -- and the pair that survived the two arms
            // above -- advances the inner iterator.
            i2 += 1;
        } // it2

        if !is_complex && meme_s {
            let mut trmeme_s = tr_meme_s.clone();
            {
                let i = h1.read().expect("interference lock");
                trmeme_s.set_index(i.transition().index()?);
            }
            h1.write()
                .expect("interference lock")
                .set_transition(trmeme_s);
            rli.push(h1.clone());
            li.remove(i1);
            continue;
        } else if is_complex && !meme_s {
            ei_tool.transition(&mut h1.write().expect("interference lock"));
            rli.push(h1.clone());
            li.remove(i1);
            continue;
        }
        i1 += 1;
    } // it1

    Ok(())
}

/// `FDS_repvg(BDS, EIX, GT, LOI, RLOI)` (`cxx:246-303`).
pub fn fds_repvg(
    bds: &DataStructure,
    eix: i32,
    gt: Kind,
    loi: &mut InterferenceList,
    rloi: &mut InterferenceList,
) -> Result<(), String> {
    // `TopOpeBRepDS_TKI tki; tki.FillOnGeometry(LOI);` (`cxx:252-253`) is
    // immediately superseded by the index grouping below (`cxx:255-275`).
    let mut groups: Vec<(i32, InterferenceList)> = Vec::new();
    for h in loi.iter() {
        let isa = {
            let i = h.read().expect("interference lock");
            i.transition().index()?
        };
        match groups.iter_mut().find(|(k, _)| *k == isa) {
            Some((_, l)) => l.push(h.clone()),
            None => groups.push((isa, vec![h.clone()])),
        }
    }

    loi.clear();
    for (_isa, mut group_loi) in groups {
        if group_loi.len() < 2 {
            continue;
        }
        let mut rloi_group: InterferenceList = Vec::new();
        fds_repvg2(bds, eix, gt, &mut group_loi, &mut rloi_group)?;
        loi.extend(group_loi);
        rloi.extend(rloi_group);
    }
    Ok(())
}
