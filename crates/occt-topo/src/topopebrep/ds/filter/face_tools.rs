//! The face-interference processing helpers of
//! `TopOpeBRepDS_ProcessFaceInterferences.cxx`.
//!
//! Sources (OCCT 8.0.0):
//! - `FUN_mkTonF` (`cxx:69-143`)
//! - `FUN_edgeofface` (`cxx:146-163`)
//! - `FUN_keepFinterference` (`cxx:165-186`)
//! - `FUN_unkeepFdoubleGBoundinterferences` (`cxx:189-277`)
//! - `FUN_resolveFUNKNOWN` (`cxx:280-431`)

use std::cell::RefCell;
use std::rc::Rc;

use crate::abs::{Orientation, ShapeType, State};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape};
use crate::topexp::Explorer;
use crate::topo_tools::{self, shape_map_key};
use crate::topopebrep::tool::curvature::nt;
use crate::topopebrep::tool::project_curve::fun_tool_bounds;
use crate::topopebrep::tool::shape_classifier::ShapeClassifier;
use crate::topopebrep::tool::tangent::{curves_same_oriented, tg_geom_e};
use crate::topopebrep::tool::topology::{fun_tool_ng_geom_f, fun_tool_par_f};

use crate::topopebrep::ds::export::{fds_data, fds_idata, fds_tdata};
use crate::topopebrep::ds::filter::funk::fun_unk_fstasta;
use crate::topopebrep::ds::interference::Interference;
use crate::topopebrep::ds::interference_iterator::InterferenceList;
use crate::topopebrep::ds::kind::Kind;
use crate::topopebrep::ds::process_interferences_tool::fun_unkeep_unknown;
use crate::topopebrep::ds::tool::DShapeOn1StateMap;
use crate::topopebrep::ds::transition::Transition;

/// `FUN_mkTonF(F, FS, E, T)` (`cxx:69-143`).
///
/// `None` is every `return false` of the C++ (`cxx:75-78`, `:91-94`,
/// `:98-101`, `:104-107`, `:113-116`, `:121-124`, `:128-131`). A `gp_Dir`
/// construction that raises in OCCT (an empty cross product at `cxx:124`, a
/// bundle comparison at `cxx:111`) also answers `None` instead of aborting.
pub fn fun_mk_t_on_f(ff: &Face, fs: &Face, e: &Edge) -> Option<Transition> {
    if BRepTool::is_degenerated(e) {
        return None;
    }
    let mut t = Transition::new(State::Unknown, State::Unknown, ShapeType::Face, ShapeType::Face);

    let tola = 1.0e-6; // nyitol
    let (f, l) = fun_tool_bounds(&e.0);
    let par_t = 0.456789;
    let pmil = (1.0 - par_t) * f + par_t * l;
    let tg_e = tg_geom_e(pmil, &e.0)?;
    let uv_f = fun_tool_par_f(&e.0, pmil, &ff.0)?;
    let uv_fs = fun_tool_par_f(&e.0, pmil, &fs.0)?;

    let ng_f = fun_tool_ng_geom_f(&uv_f, ff)?;
    let xx = ng_f.xyz().dot(tg_e.xyz()).abs();
    let tgt = (1.0 - xx).abs() < tola;
    if tgt {
        return None;
    }

    let nt_fs = nt(&uv_fs, fs)?;
    // `gp_Dir beafter = ngF ^ tgE` (`cxx:124`): a cross product that raises on
    // parallel operands, which the `tgt` test above already rules out.
    let beafter = ng_f.xyz().crossed(tg_e.xyz());
    let yy = beafter.dot(nt_fs.xyz());
    let unk = yy.abs() < tola;
    if unk {
        return None;
    }

    if yy < 0.0 {
        t.set_orientation(Orientation::Forward);
    } else {
        t.set_orientation(Orientation::Reversed);
    }
    Some(t)
}

/// `FUN_edgeofface(E, F)` (`cxx:145-163`): whether `e` is one of `f`'s edges.
pub fn fun_edgeofface(e: &TopoShape, f: &TopoShape) -> bool {
    let mut ex = Explorer::new(f, ShapeType::Edge);
    while ex.more() {
        let cur = ex.current().clone();
        ex.next();
        if topo_tools::is_same(&cur, e) {
            return true;
        }
    }
    false
}

/// `FUN_keepFinterference(DS, I, F)` (`cxx:165-187`).
pub fn fun_keep_f_interference(
    ds: &crate::topopebrep::ds::data_structure::DataStructure,
    i: &Interference,
    f: &TopoShape,
) -> bool {
    let mut res = true;
    if i.is_face_edge() {
        // I is rejected when its geometry edge is an edge of the face accessing
        // I.
        let eg = ds.shape(i.geometry(), true);
        let k3 = match eg {
            Some(eg) => !fun_edgeofface(&eg, f),
            None => true,
        };
        res = res && k3;
    }
    res
}

/// `FUN_unkeepFdoubleGBoundinterferences(LI, BDS, SIX)` (`cxx:189-278`).
///
/// Inside one `(GeometryType=EDGE, G, SupportType=FACE)` group, an interference
/// with `GBound()` set is dropped when another one with a different support and
/// `GBound()` set exists.
pub fn fun_unkeep_f_double_g_bound_interferences(li: &mut InterferenceList) {
    let mut i = 0;
    while i < li.len() {
        let (is_unk1, gt1, g1, st1, s1, gb1, is_ssi1) = {
            let h = li[i].read().expect("interference lock");
            let (gt, g, st, s) = fds_data(&h);
            (
                h.transition().is_unknown(),
                gt,
                g,
                st,
                s,
                h.g_bound(),
                h.is_shape_shape(),
            )
        };
        if is_unk1 || !is_ssi1 {
            i += 1;
            continue;
        }
        let Some(is_b1) = gb1 else {
            i += 1;
            continue;
        };

        let mut cond1 = false;
        let mut j = i + 1;
        while j < li.len() {
            let (is_unk2, gt2, g2, st2, s2, gb2, is_ssi2) = {
                let h = li[j].read().expect("interference lock");
                let (gt, g, st, s) = fds_data(&h);
                (
                    h.transition().is_unknown(),
                    gt,
                    g,
                    st,
                    s,
                    h.g_bound(),
                    h.is_shape_shape(),
                )
            };
            if is_unk2 || !is_ssi2 {
                j += 1;
                continue;
            }
            let Some(is_b2) = gb2 else {
                j += 1;
                continue;
            };
            let cond2 = gt2 == gt1
                && gt1 == Kind::Edge
                && g2 == g1
                && st2 == st1
                && st1 == Kind::Face
                && s2 != s1
                && is_b1
                && is_b2;
            if cond2 {
                cond1 = true;
                li.remove(j);
            } else {
                j += 1;
            }
        }

        if cond1 {
            li.remove(i);
        } else {
            i += 1;
        }
    }
}

/// `FUN_resolveFUNKNOWN(LI, BDS, SIX, MEsp, pClassif)` (`cxx:280-431`).
///
/// `p_classif` is shared with the caller (`TopOpeBRepDS_Filter` holds the same
/// `handle`), so `SetReference` is visible to later interferences exactly as in
/// OCCT (`cxx:412`).
pub fn fun_resolve_f_unknown(
    li: &mut InterferenceList,
    ds: &crate::topopebrep::ds::data_structure::DataStructure,
    six: i32,
    me_sp: &DShapeOn1StateMap,
    p_classif: Option<Rc<RefCell<ShapeClassifier>>>,
) {
    if let Some(f) = ds.shape(six, true) {
        let ff = Face(f.clone());
        for h in li.iter() {
            let (is_unk, tsb1, isb1, tsa1, isa1, gt1, g1, st1, s1, is_fei, gb) = {
                let i = h.read().expect("interference lock");
                let (tsb, isb, tsa, isa, gt, g, st, s) = fds_idata(&i);
                (
                    i.transition().is_unknown(),
                    tsb,
                    isb,
                    tsa,
                    isa,
                    gt,
                    g,
                    st,
                    s,
                    i.is_face_edge(),
                    i.g_bound(),
                )
            };
            if !is_unk {
                continue;
            }
            let idt = tsb1 == ShapeType::Face
                && tsa1 == ShapeType::Face
                && gt1 == Kind::Edge
                && st1 == Kind::Face;
            let idi = isb1 == s1 && isa1 == s1;
            // `etgf`: face tangent to a face at one edge.
            if !(idt && idi) {
                continue;
            }
            let Some(ee_shape) = ds.shape(g1, true) else {
                continue;
            };
            let ee = Edge(ee_shape.clone());
            let (_fe, _le) = BRepTool::edge_parameters(&ee);
            if !is_fei {
                continue;
            }
            let Some(fs_shape) = ds.shape(s1, true) else {
                continue;
            };
            let fs = Face(fs_shape.clone());
            let mut is_ee_gb = gb.unwrap_or(false);

            // `cxx:341-383`: the split-on arm of `MEsp`.
            let mut eesp = ee.clone();
            if let Some(los1) = me_sp.get(&shape_map_key(&ee.0)) {
                if los1.is_split() {
                    let los = los1.list_on_state();
                    if !los.is_empty() {
                        if let Some(first) = Edge::wrap(los[0].clone()) {
                            eesp = first;
                            if !topo_tools::is_same(&eesp.0, &ee.0) {
                                is_ee_gb = false;
                            }
                        }
                        if los.len() > 1 {
                            // MSV (`cxx:361-379`): select the split lying on both
                            // faces.
                            for cand in los {
                                let Some(a_e) = Edge::wrap(cand.clone()) else {
                                    continue;
                                };
                                let (af, al) = fun_tool_bounds(&a_e.0);
                                let par_t = 0.456789;
                                let pmil = (1.0 - par_t) * af + par_t * al;
                                if fun_tool_par_f(&a_e.0, pmil, &ff.0).is_some()
                                    && fun_tool_par_f(&a_e.0, pmil, &fs.0).is_some()
                                {
                                    eesp = a_e;
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            let mut is_so = true;
            if !topo_tools::is_same(&eesp.0, &ee.0) {
                match curves_same_oriented(&eesp.0, &ee.0) {
                    Some(v) => is_so = v,
                    None => continue,
                }
            }

            let (stateb, statea);
            match fun_mk_t_on_f(&ff, &fs, &eesp) {
                Some(t) => {
                    stateb = t.before_state();
                    statea = t.after_state();
                }
                None => {
                    // `cxx:401-418`
                    let mut classif = p_classif.clone();
                    if let Some(c) = classif.as_ref() {
                        let rank_fs = ds.ancestor_rank_of(&fs_shape);
                        // MSV: find the solids of the same object rank as FS to
                        // determine the transition relatively to the solid
                        // rather than to the face (see CFE002 C2, SIX == 13).
                        let a_sref = ds.shape(rank_fs, true);
                        if let Some(a_sref) = a_sref {
                            let mut ex = Explorer::new(&a_sref, ShapeType::Solid);
                            if ex.more() {
                                c.borrow_mut().set_reference(&a_sref);
                            }
                        }
                    }
                    let mut guard = classif.as_ref().map(|c| c.borrow_mut());
                    let (b, a) = fun_unk_fstasta(
                        &ff,
                        &fs,
                        &eesp,
                        is_ee_gb,
                        guard.as_deref_mut(),
                    );
                    stateb = b;
                    statea = a;
                }
            }
            if stateb == State::Unknown || statea == State::Unknown {
                continue;
            }

            let (stateb, statea) = if !is_so {
                (statea, stateb)
            } else {
                (stateb, statea)
            };
            let mut i1 = h.write().expect("interference lock");
            i1.change_transition().set(stateb, statea, tsb1, tsa1);
        }
    }
    fun_unkeep_unknown(li);
}
