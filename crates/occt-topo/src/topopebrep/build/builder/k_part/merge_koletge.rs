//! `TopOpeBRepBuild_Builder` — the `iskoletge` KPart case.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_KPart.cxx`
//! (OCCT 8.0.0):
//! - `MergeKPartiskoletge` (`:578-836`),
//! - `FUN_addf` (`:268-281`),
//! - `FUN_comparekoletgesh` (`:283-315`),
//! - `FUN_changev` (`:317-333`), `FUN_updatev` (`:335-362`),
//! - `FUN_makefaces` (`:364-520`),
//! - `FUN_rebuildfc` (`:522-570`).
//!
//! The case handles two solids glued on a same-domain planar face (`KPiskoletge`,
//! `:1253-1385`): the outer wire of that face is split along the section edges,
//! so both shells are rebuilt from the new faces before the merged solid is
//! handed back.
//!
//! UNPORTED here, with its OCCT location:
//! - `TopOpeBRepTool::CorrectONUVISO(F, Newf)` (`:509-515`): the port appends
//!   `newf`, which is the `else` arm of that call.

use std::collections::HashMap;

use occt_core::gp::{GpDir, GpPnt};

use crate::abs::{Orientation, ShapeType, State};
use crate::brep_class3d::find_a_point_in_the_face;
use crate::brep_tool::BRepTool;
use crate::brep_tools;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{explorer_shapes, is_same, shape_map_key, ShapeMapKey};
use crate::topopebrep::ds::connex::Connexity;
use crate::topopebrep::ds::export::fun_ds_getoov_in_hds;
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::ds::kind::Config;
use crate::topopebrep::tool::project_curve::{
    fun_tool_bounds, fun_tool_proj_pon_e, fun_tool_value,
};
use crate::topopebrep::tool::shape_classifier::ShapeClassifier;
use crate::topopebrep::tool::tangent::fun_tool_tg_geom_e;
use crate::topopebrep::tool::topology::{
    fun_tool_comparebndkole, fun_tool_pcurveon_f,
};
use crate::topopebrep::ds::build_tool::BuildTool;

use super::super::k_part::{KP_ISKOLETGE, RESNEWSOL, RESNULL, RESSHAPE1, RESSHAPE2};
use super::predicates::fun_sortplcy;
use super::Builder;

/// `FUN_addf(sta, ftoadd, map)` (`KPart.cxx:268-281`): insert `ftoadd` into the
/// `addedfaces` map, complemented when the state asks for `TopAbs_IN`.
///
/// `NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher>`
/// keys by `IsSame`, so a repeated face rebinds the value at its original
/// position; that is what the `(key, index)` pair reproduces.
fn fun_addf(
    sta: State,
    f_to_add: &TopoShape,
    added: &mut Vec<TopoShape>,
    index: &mut HashMap<ShapeMapKey, usize>,
) {
    // cxx:274-280
    let mut fori = f_to_add.clone();
    if sta == State::In {
        fori = f_to_add.oriented(f_to_add.orientation().complement());
    }
    let key = shape_map_key(&fori);
    match index.get(&key) {
        Some(i) => added[*i] = fori,
        None => {
            index.insert(key, added.len());
            added.push(fori);
        }
    }
}

/// `FUN_comparekoletgesh(SC, sh1, sh2, fkole1, fkole2)` (`KPart.cxx:283-315`):
/// 1 when `sh1` is inside `sh2`, 2 when `sh2` is inside `sh1`, 0 otherwise.
fn fun_comparekoletgesh(
    sc: &mut ShapeClassifier,
    sh1: &TopoShape,
    sh2: &TopoShape,
    fkole1: &TopoShape,
) -> i32 {
    // cxx:291
    sc.set_reference(sh2);

    // cxx:292-314
    for f1 in explorer_shapes(sh1, ShapeType::Face) {
        if is_same(&f1, fkole1) {
            continue;
        }
        let Some(face) = Face::wrap(f1) else {
            continue;
        };
        let mut pnt1 = GpPnt::new(0.0, 0.0, 0.0);
        if !find_a_point_in_the_face(&face, &mut pnt1) {
            continue;
        }
        if sc.state_p3d_reference(&pnt1).is_err() {
            continue;
        }
        let stpnt1 = sc.state();
        if stpnt1 == State::In {
            return 2;
        }
        if stpnt1 == State::Out {
            return 1;
        }
    }
    0
}

/// `FUN_changev(HDS, v)` (`KPart.cxx:317-333`): `v` is a rank-2 vertex that is
/// same-domain in `HDS`.
fn fun_changev(hds: &HDataStructure, v: &TopoShape) -> bool {
    // cxx:319-322
    if !hds.has_shape(v, true) {
        return false;
    }
    // cxx:323-327
    if !hds.has_same_domain(v, true) {
        return false;
    }
    // cxx:328-330
    hds.ds().ancestor_rank_of(v) == 2
}

/// `FUN_updatev(HDS, newed, v, oriv, parv, changev)` (`KPart.cxx:335-362`):
/// attach either the original vertex `v` or its same-domain partner, then set
/// its parameter on the new edge.
fn fun_updatev(
    hds: &HDataStructure,
    newed: &mut TopoShape,
    v: &TopoShape,
    oriv: Orientation,
    parv: f64,
    changev: bool,
) -> bool {
    let bt = BuildTool::new();
    // cxx:343-356
    let ov = if changev {
        let Some(oov) = fun_ds_getoov_in_hds(hds, v) else {
            return false;
        };
        oov.oriented(oriv)
    } else {
        v.oriented(oriv)
    };
    bt.add_edge_vertex(newed, &ov);
    bt.parameter(newed, &ov, parv);
    true
}

/// `FUN_tool_SameOri(E1, E2)` (`TopOpeBRepTool_TOPOLOGY.cxx:1285-1336`): the
/// tangents at a mid parameter of `E2` and at `E2`'s point projected on `E1`
/// point the same way.
fn fun_tool_same_ori(e1: &TopoShape, e2: &TopoShape) -> bool {
    let (f, l) = fun_tool_bounds(e2);
    let x = 0.345;
    let mid = x * f + (1.0 - x) * l;
    let Some(pmid) = fun_tool_value(mid, e2) else {
        return false;
    };

    let tmp = fun_tool_tg_geom_e(mid, e2);
    let Ok(mut tg_e2) = GpDir::from_vec(&tmp) else {
        return false;
    };
    if e2.orientation().is_reversed() {
        tg_e2.reverse();
    }

    let Some((p_e1, dist)) = fun_tool_proj_pon_e(&pmid, e1) else {
        return false;
    };
    let tol1 = Edge::wrap(e1.clone()).map(|e| BRepTool::edge_tolerance(&e));
    let tol2 = Edge::wrap(e2.clone()).map(|e| BRepTool::edge_tolerance(&e));
    let (Some(tol1), Some(tol2)) = (tol1, tol2) else {
        return false;
    };
    if dist > tol1.max(tol2) * 10.0 {
        return false;
    }

    let tmp = fun_tool_tg_geom_e(p_e1, e1);
    let Ok(mut tg_e1) = GpDir::from_vec(&tmp) else {
        return false;
    };
    if e1.orientation().is_reversed() {
        tg_e1.reverse();
    }
    tg_e2.dot(&tg_e1) > 0.0
}

/// `FUN_makefaces(BU, HDS, EnewE, f, outerwi, newfaces)` (`KPart.cxx:364-520`):
/// rebuild the face `f` from `f`'s outer wire, whose edges are the splits or the
/// vertex-corrected copies bound in `EnewE`.
fn fun_makefaces(
    builder: &Builder,
    hds: &HDataStructure,
    enewe: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
    f: &TopoShape,
    newfaces: &mut Vec<TopoShape>,
) -> bool {
    let bt = BuildTool::new();
    // cxx:381-384: work on the FORWARD face, remember the original orientation.
    let fwd = f.oriented(Orientation::Forward);
    let Some(face_fwd) = Face::wrap(fwd.clone()) else {
        return false;
    };
    let of = f.orientation();

    // cxx:390-392
    let Some(outerwf) = brep_tools::outer_wire(&face_fwd) else {
        return false;
    };

    let mut loe: Vec<TopoShape> = Vec::new();
    // cxx:393-472
    for ed in explorer_shapes(&outerwf.0, ShapeType::Edge) {
        let oe = ed.orientation();
        let issplit = builder.is_split(&ed, State::On);
        let isbound = enewe.contains_key(&shape_map_key(&ed));

        // cxx:398-431
        let mut bound_now = isbound;
        if !isbound && !issplit {
            let verts = explorer_shapes(&ed, ShapeType::Vertex);
            if verts.len() < 2 {
                return false;
            }
            let vf = verts[0].clone();
            let vl = verts[1].clone();
            let changevf = fun_changev(hds, &vf);
            let changevl = fun_changev(hds, &vl);
            if changevf || changevl {
                let ff = BRepTool::parameter_of_vertex(&ed, &vf);
                let l = BRepTool::parameter_of_vertex(&ed, &vl);

                let mut newed = bt.copy_edge(&ed.oriented(Orientation::Forward));
                if !fun_updatev(hds, &mut newed, &vf, Orientation::Forward, ff, changevf) {
                    return false;
                }
                if !fun_updatev(hds, &mut newed, &vl, Orientation::Reversed, l, changevl) {
                    return false;
                }
                newed.set_orientation(oe);
                enewe.insert(shape_map_key(&ed), vec![newed]);
                bound_now = true;
            }
        }

        // cxx:433-472
        if issplit {
            let speds = builder.splits(&ed, State::On).to_vec();
            if speds.is_empty() {
                return false;
            }
            let spe = speds[0].clone();
            // `FUN_tool_SameOri(ed, spe)` (`TOPOLOGY.cxx:1285-1336`).
            let sameori = fun_tool_same_ori(&ed, &spe);
            let mut orisp = spe.orientation();
            if !sameori {
                orisp = orisp.complement();
            }

            for esp in speds {
                if !fun_tool_pcurveon_f(&face_fwd.0, &esp) {
                    return false;
                }
                loe.push(esp.oriented(orisp));
            }
        } else if bound_now {
            let Some(neweds) = enewe.get(&shape_map_key(&ed)) else {
                return false;
            };
            for newe in neweds.clone() {
                loe.push(newe.oriented(oe));
            }
        } else {
            loe.push(ed.clone());
        }
    }

    // cxx:475-482
    let mut neww = TopoShape::new(ShapeType::Wire);
    bt.make_wire(&mut neww);
    for e in &loe {
        bt.add_wire_edge(&mut neww, e);
    }

    // cxx:484-488
    let mut newf = bt.copy_face(&fwd);
    bt.add_face_wire(&mut newf, &neww);

    // cxx:490-501: the other wires of `f`.
    for oow in explorer_shapes(&fwd, ShapeType::Wire) {
        if is_same(&oow, &outerwf.0) {
            continue;
        }
        bt.add_face_wire(&mut newf, &oow);
    }

    // cxx:503-506
    if of.is_reversed() {
        newf.set_orientation(Orientation::Reversed);
    }

    // cxx:508-516: `TopOpeBRepTool::CorrectONUVISO(F, Newf)` is UNPORTED, so
    // the `else` arm of that call is what the port runs.
    newfaces.push(newf);
    true
}

/// `FUN_rebuildfc(BU, HDS, Stk, Fk, EnewE, fcnewfc)` (`KPart.cxx:522-570`):
/// for every edge of `Fk`'s outer wire, rebuild the one face of the other
/// operand connected to it and record `{connex face -> new faces}`.
fn fun_rebuildfc(
    builder: &Builder,
    hds: &HDataStructure,
    fk: &TopoShape,
    enewe: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
    fcnewfc: &mut Vec<(TopoShape, Vec<TopoShape>)>,
) -> bool {
    // cxx:540-543
    let Some(face_fk) = Face::wrap(fk.clone()) else {
        return false;
    };
    let Some(owk) = brep_tools::outer_wire(&face_fk) else {
        return false;
    };
    let eds = explorer_shapes(&owk.0, ShapeType::Edge);

    // cxx:545
    let mut connex = Connexity::new();
    connex.prepare(hds);

    // cxx:547-569
    for ed in eds {
        let mut lfc_fk: Vec<TopoShape> = Vec::new();
        connex.face_edge_connex_faces(fk, &ed, hds, &mut lfc_fk);
        // cxx:552-555
        if lfc_fk.len() != 1 {
            return false;
        }
        let fc_fk = lfc_fk[0].clone();

        let mut newfc_fk: Vec<TopoShape> = Vec::new();
        if !fun_makefaces(builder, hds, enewe, &fc_fk, &mut newfc_fk) {
            return false;
        }
        // `IndexedDataMap::Add`: an existing key appends to its list.
        match fcnewfc
            .iter_mut()
            .find(|(k, _)| is_same(k, &fc_fk))
        {
            Some((_, l)) => l.extend(newfc_fk),
            None => fcnewfc.push((fc_fk, newfc_fk)),
        }
    }
    true
}

/// The first shell of `sha` holding `f` — `MapShapesAndAncestors(sha, FACE,
/// SHELL, M).FindFromKey(f).First()` (`KPart.cxx:680-694`).
fn ancestor_shell(sha: &TopoShape, f: &TopoShape) -> Option<TopoShape> {
    for she in explorer_shapes(sha, ShapeType::Shell) {
        let faces = explorer_shapes(&she, ShapeType::Face);
        if faces.iter().any(|ff| is_same(ff, f)) {
            return Some(she);
        }
    }
    None
}

impl Builder {
    /// `MergeKPartiskoletge()` (`KPart.cxx:578-836`).
    pub fn merge_k_part_iskoletge(&mut self) {
        // cxx:590-594
        if self.my_is_k_part != KP_ISKOLETGE {
            return;
        }
        let (Some(shape1), Some(shape2)) = (self.my_shape1.clone(), self.my_shape2.clone()) else {
            return;
        };
        let st1 = self.my_state1;
        let st2 = self.my_state2;

        // cxx:597-601
        self.map_shapes(Some(&shape1), Some(&shape2));
        self.change_merged(&shape1, st1);
        self.change_merged(&shape2, st2);

        // cxx:603-608
        let mut lshsd1: Vec<TopoShape> = Vec::new();
        let mut lfhsd1: Vec<TopoShape> = Vec::new();
        if !self.kp_iskoletgesh(&shape1, &mut lshsd1, &mut lfhsd1) {
            return;
        }
        let mut lshsd2: Vec<TopoShape> = Vec::new();
        let mut lfhsd2: Vec<TopoShape> = Vec::new();
        if !self.kp_iskoletgesh(&shape2, &mut lshsd2, &mut lfhsd2) {
            return;
        }

        // cxx:610-614: only the first same-domain solid of each operand.
        let (Some(sol1), Some(sol2)) = (lshsd1.first().cloned(), lshsd2.first().cloned()) else {
            return;
        };
        self.change_merged(&sol1, st1);
        self.change_merged(&sol2, st2);

        // cxx:616-621
        let (lplhsd1, _lcyhsd1) = fun_sortplcy(&lfhsd1);
        let (lplhsd2, _lcyhsd2) = fun_sortplcy(&lfhsd2);
        let (Some(fac1), Some(fac2)) = (lplhsd1.first().cloned(), lplhsd2.first().cloned()) else {
            return;
        };

        // cxx:640-642
        let Some(hds) = self.my_data_structure.clone() else {
            return;
        };
        let config2 = hds.same_domain_orientation(&fac2);
        let same_oriented = config2 == Config::SameOriented;

        // cxx:644-662
        let mut rgre = 1;
        if same_oriented {
            let mut sc = self.my_shape_classifier.borrow_mut();
            rgre = fun_comparekoletgesh(&mut sc, &shape1, &shape2, &fac1);
        }
        if same_oriented && rgre == 0 {
            rgre = fun_tool_comparebndkole(&shape1, &shape2);
        }
        if rgre == 0 {
            return;
        }

        // cxx:666-676
        let rsma = if rgre == 1 { 2 } else { 1 };
        let pf_sma = if rsma == 1 { fac1.clone() } else { fac2.clone() };
        let pf_gre = if rgre == 1 { fac1.clone() } else { fac2.clone() };
        let sta_sma = if rsma == 1 { st1 } else { st2 };
        let sta_gre = if rgre == 1 { st1 } else { st2 };
        let sha_sma = if rsma == 1 { shape1.clone() } else { shape2.clone() };
        let sha_gre = if rgre == 1 { shape1.clone() } else { shape2.clone() };

        // cxx:678-679
        let ires = self.kp_iskoletge_analyse(config2, sta_sma, sta_gre);

        // cxx:680-694
        let Some(she_sma) = ancestor_shell(&sha_sma, &pf_sma) else {
            return;
        };
        let Some(she_gre) = ancestor_shell(&sha_gre, &pf_gre) else {
            return;
        };
        self.change_merged(&she_sma, sta_sma);
        self.change_merged(&she_gre, sta_gre);

        // cxx:696-836
        let mut newshe = TopoShape::new(ShapeType::Shell);
        if ires == RESNULL {
            return;
        } else if ires == RESSHAPE1 {
            // cxx:703-708
            self.build_tool().make_shell(&mut newshe);
            newshe = she_sma.clone();
        } else if ires == RESSHAPE2 {
            // cxx:709-714
            self.build_tool().make_shell(&mut newshe);
            newshe = she_gre.clone();
        } else if ires == RESNEWSOL {
            let bt = BuildTool::new();
            let mut added: Vec<TopoShape> = Vec::new();
            let mut added_index: HashMap<ShapeMapKey, usize> = HashMap::new();
            let mut enewe: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();

            // cxx:723-750: faces of SMA connected to fSMA.
            let mut fcnewfc_sma: Vec<(TopoShape, Vec<TopoShape>)> = Vec::new();
            if !fun_rebuildfc(
                self,
                &hds,
                &pf_sma,
                &mut enewe,
                &mut fcnewfc_sma,
            ) {
                return;
            }
            for (f, newlf) in &fcnewfc_sma {
                for ff in newlf {
                    fun_addf(sta_sma, ff, &mut added, &mut added_index);
                    self.change_merged(f, sta_sma).push(ff.clone());
                }
            }

            // cxx:751-770: faces of GRE connected to fGRE.
            let mut fcnewfc_gre: Vec<(TopoShape, Vec<TopoShape>)> = Vec::new();
            if !fun_rebuildfc(
                self,
                &hds,
                &pf_gre,
                &mut enewe,
                &mut fcnewfc_gre,
            ) {
                return;
            }
            for (f, newlf) in &fcnewfc_gre {
                for ff in newlf {
                    fun_addf(sta_gre, ff, &mut added, &mut added_index);
                    self.change_merged(f, sta_gre).push(ff.clone());
                }
            }

            // cxx:772-803: the untouched faces of both operands.
            for ff in explorer_shapes(&sha_sma, ShapeType::Face) {
                if fcnewfc_sma.iter().any(|(k, _)| is_same(k, &ff)) {
                    continue;
                }
                if is_same(&ff, &pf_sma) {
                    continue;
                }
                fun_addf(sta_sma, &ff, &mut added, &mut added_index);
            }
            for ff in explorer_shapes(&sha_gre, ShapeType::Face) {
                if fcnewfc_gre.iter().any(|(k, _)| is_same(k, &ff)) {
                    continue;
                }
                if is_same(&ff, &pf_gre) {
                    continue;
                }
                fun_addf(sta_gre, &ff, &mut added, &mut added_index);
            }

            // cxx:805-819
            if !added.is_empty() {
                bt.make_shell(&mut newshe);
                bt.closed(&newshe, true);
                for ftoadd in &added {
                    bt.add_shell_face(&mut newshe, ftoadd);
                }
            }
        }

        // cxx:821-828
        let mut newsol = TopoShape::new(ShapeType::Solid);
        if !newshe.is_null() {
            self.build_tool().make_solid(&mut newsol);
            self.build_tool().add_solid_shell(&mut newsol, &newshe);
        }

        // cxx:830-834
        if !newsol.is_null() {
            self.change_merged(&shape1, st1).push(newsol);
        }
    }
}
