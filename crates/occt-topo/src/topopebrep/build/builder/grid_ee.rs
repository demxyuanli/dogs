//! `TopOpeBRepBuild_Builder` — the point/vertex topology fillers of
//! `TopOpeBRepBuild_GridEE.cxx`.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_GridEE.cxx`
//! (OCCT 8.0.0):
//!
//! | OCCT | lines |
//! |---|---|
//! | `GFillEdgePVS` | `198-214` |
//! | `FUN_MoreSHAINT` | `218-246` |
//! | `FUN_isonbound` | `305-319` |
//! | `GFillPointTopologyPVS(E, G, PVS)` | `331-776` |
//! | `GFillPointTopologyPVS(E, EPit, G1, PVS)` | `778-1051` |
//!
//! OCCT keeps three file-scope globals the two fillers read (`GLOBAL_issp`,
//! `GLOBAL_faces2d`, `GLOBAL_hassd`). The port keeps them as thread-local cells
//! so a nested `SplitSectionEdge` call cannot see another thread's value, which
//! is what the single-threaded OCCT global already guarantees.
//!
//! The two `GFillPointTopologyPVS` overloads recurse into each other through
//! `TopOpeBRepDS_PointIterator`; the port passes the cursor by shared reference
//! (`cxx:665`, `cxx:759`, `cxx:766`), which is what the OCCT `const&` does.

use std::cell::Cell;
use std::sync::{Arc, RwLock};

use crate::abs::{Orientation, State};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, TopoShape};
use crate::topo_tools::{is_equal, is_same};
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::ds::interference::Interference;
use crate::topopebrep::ds::interference_iterator::InterferenceList;
use crate::topopebrep::ds::kind::Kind;
use crate::topopebrep::ds::point_iterator::PointIterator;
use crate::topopebrep::ds::process_interferences_tool::{
    fun_select_sk_interference, fun_select_trasha_interference, make_epv_interference,
};
use crate::topopebrep::ds::export::{fds_parameter, fun_ds_getoov};
use crate::topopebrep::ds::redu::{fun_ds_get_tr, GetTr};
use crate::topopebrep::ds::tki::TKI;
use crate::topopebrep::ds::transition::Transition;
use crate::topopebrep::tool::closed::closed_e;
use crate::topopebrep::tool::project_curve::fun_tool_bounds;

use super::super::g_topo::GTopo;
use super::super::loop_shape::Loop;
use super::super::pave_set::PaveSet;
use super::grid::{first_same_domain, fun_aresamegeom, is_diff_oriented};
use super::Builder;

/// `#define TheIN (1)` (`GridEE.cxx:321`).
const THE_IN: i32 = 1;
/// `#define TheON (2)` (`GridEE.cxx:322`).
const THE_ON: i32 = 2;
/// `#define TheOUT (3)` (`GridEE.cxx:323`).
const THE_OUT: i32 = 3;

thread_local! {
    /// `int GLOBAL_issp` (`Section.cxx:920`).
    static GLOBAL_ISSP: Cell<i32> = const { Cell::new(0) };
    /// `bool GLOBAL_faces2d` (`GridFF.cxx`).
    static GLOBAL_FACES2D: Cell<bool> = const { Cell::new(false) };
    /// `int GLOBAL_hassd` (`Section.cxx:926`).
    static GLOBAL_HASSD: Cell<i32> = const { Cell::new(0) };
}

/// Current `GLOBAL_issp`.
pub(super) fn global_issp() -> i32 {
    GLOBAL_ISSP.with(|c| c.get())
}

/// Assign `GLOBAL_issp`.
#[allow(dead_code)] // written by `SplitSectionEdge` (`Section.cxx:1107-1207`), not ported yet.
pub fn set_global_issp(v: i32) {
    GLOBAL_ISSP.with(|c| c.set(v));
}

/// Current `GLOBAL_faces2d`.
pub(super) fn global_faces2d() -> bool {
    GLOBAL_FACES2D.with(|c| c.get())
}

/// Assign `GLOBAL_faces2d`.
#[allow(dead_code)] // written by the `GridFF.cxx` face drivers, not ported yet.
pub fn set_global_faces2d(v: bool) {
    GLOBAL_FACES2D.with(|c| c.set(v));
}

/// Current `GLOBAL_hassd`.
#[allow(dead_code)] // written by `SplitSectionEdge` (`Section.cxx:957-964`), not ported yet.
pub(super) fn global_hassd() -> i32 {
    GLOBAL_HASSD.with(|c| c.get())
}

/// Assign `GLOBAL_hassd`.
#[allow(dead_code)] // see [`global_hassd`].
pub fn set_global_hassd(v: i32) {
    GLOBAL_HASSD.with(|c| c.set(v));
}

/// `FUN_MoreSHAINT(EPit, SHA)` (`GridEE.cxx:218-246`): advance `EPit` past the
/// interferences whose transition is not defined on `SHA` on both sides.
pub(super) fn fun_more_shaint(epit: &mut PointIterator, sha: crate::abs::ShapeType) -> bool {
    while epit.more() {
        let t = epit.current().transition();
        let rejet = t.shape_before() != sha || t.shape_after() != sha;
        if rejet {
            epit.next();
        } else {
            return true;
        }
    }
    false
}

/// `FUN_isonbound(HDS, I)` (`GridEE.cxx:305-319`): an `EdgeVertexInterference`
/// whose vertex is bound and has no same-domain partner.
pub(super) fn fun_isonbound(hds: &HDataStructure, i: &Interference) -> bool {
    let g = i.geometry();
    let kg = i.geometry_type();
    let mut gb1 = kg == Kind::Vertex;
    if gb1 {
        let vhassd = match hds.ds().shape(g, true) {
            Some(v) => hds.has_same_domain(&v, true),
            None => false,
        };
        gb1 = i.g_bound().unwrap_or(false);
        gb1 = gb1 && !vhassd;
    }
    gb1
}

impl Builder {
    /// `GFillEdgePVS(E, LE2, G, PVS)` (`GridEE.cxx:198-214`).
    ///
    /// `LE2` is unused by OCCT (`cxx:199` comments the parameter out), so the
    /// port does not take it.
    pub fn fill_edge_pvs(&self, e: &TopoShape, g: &GTopo, pvs: &mut PaveSet) {
        let mut ef = e.clone();
        ef.set_orientation(Orientation::Forward);
        self.fill_point_topology_pvs(&ef, g, pvs);
    }

    /// `GFillPointTopologyPVS(E, G, PVS)` (`GridEE.cxx:331-776`).
    pub fn fill_point_topology_pvs(&self, e: &TopoShape, g: &GTopo, pvs: &mut PaveSet) {
        let Some(hds) = self.data_structure() else {
            return;
        };
        let (shape_interf, _t2) = g.type_();
        let bds = hds.ds();
        let i_eds = bds.shape_index(e, true);
        let is_se = bds.is_section_edge(e, true);
        let dg_e = BRepTool::is_degenerated(&Edge(e.clone()));
        // `cxx:352` re-reads the same flag into `isEd`.
        let is_ed = BRepTool::is_degenerated(&Edge(e.clone()));
        let isfafa = bds.isfafa();

        // `cxx:384-389`
        let isspin = global_issp() == THE_IN;
        let isspou = global_issp() == THE_OUT;
        let isspon = global_issp() == THE_ON;
        if is_se && global_issp() == 0 {
            return;
        }

        // 0. cxx:409-436
        let lie: InterferenceList = bds.shape_interferences(e, true).to_vec();
        let scanall = isspin || isspou || isspon;
        let loi_shaint: InterferenceList = if scanall {
            lie.clone()
        } else {
            let mut epit = PointIterator::new(lie.clone());
            let mut l: InterferenceList = Vec::new();
            let mut addi = fun_more_shaint(&mut epit, shape_interf);
            while addi {
                l.push(epit.value().clone());
                epit.next();
                addi = fun_more_shaint(&mut epit, shape_interf);
            }
            l
        };

        let mut tki = TKI::new();
        tki.fill_on_geometry(&loi_shaint);

        // - kp1 - cxx:438-474
        let vclo = closed_e(e);
        let mut kp1 = 0i32;
        if let Some(vclo) = vclo.as_ref() {
            tki.init();
            let mut ng = 0;
            let mut found: Option<i32> = None;
            while tki.more() {
                ng += 1;
                let Some((kcur, gcur)) = tki.value().map(|(k, g, _)| (k, g)) else {
                    break;
                };
                if kcur == Kind::Point {
                    tki.next();
                    continue;
                }
                if let Some(v) = bds.shape(gcur, true) {
                    let oov = fun_ds_getoov(bds, &v);
                    let samev = is_same(&v, vclo);
                    let sameoov = oov.as_ref().map(|o| is_same(o, vclo)).unwrap_or(false);
                    if samev || sameoov {
                        found = Some(gcur);
                    }
                }
                tki.next();
            }
            kp1 = found.unwrap_or(0);
            if ng == 1 {
                kp1 = 0;
            }
        }

        // - kp6 - cxx:476-506
        let mut kp6 = !is_se;
        if kp6 {
            kp6 = false;
            for i in loi_shaint.iter() {
                let Ok(i) = i.read() else { continue };
                if i.support_type() != Kind::Face {
                    continue;
                }
                let o = i
                    .transition()
                    .orientation(State::In, crate::abs::ShapeType::Face);
                let forrev = o == Orientation::Forward || o == Orientation::Reversed;
                if !forrev {
                    continue;
                }
                if !fun_isonbound(hds, &i) {
                    kp6 = true;
                    break;
                }
            }
        }

        // 1. cxx:509-776
        tki.init();
        while tki.more() {
            let Some((kcur, gcur, lref)) = tki.value() else {
                break;
            };
            let licur: InterferenceList = lref.to_vec();

            let point = kcur == Kind::Point;
            let vertex = kcur == Kind::Vertex;
            if vertex {
                // `cxx:520-524`: computed for xpu221098 and never read again.
                if let Some(v) = bds.shape(gcur, true) {
                    let _v_gsd = fun_ds_getoov(bds, &v);
                }
            }

            if vertex && kp1 == gcur {
                tki.next();
                continue;
            }

            let Some(i_first) = licur.first().cloned() else {
                tki.next();
                continue;
            };
            let par_se = {
                let i = i_first.read().expect("interference lock");
                fds_parameter(&i).unwrap_or(0.0)
            };

            let mut licur_copy: InterferenceList = licur.clone();
            let mut l3d_fcur: InterferenceList = Vec::new();
            let n3d = fun_select_sk_interference(&mut licur_copy, Kind::Face, &mut l3d_fcur);
            licur_copy = licur.clone();
            let l2d_fe_cur =
                crate::topopebrep::ds::complete_ds::fun_ds_has_i2d(&licur_copy);
            let n2d = l2d_fe_cur.len() as i32;
            licur_copy = licur.clone();
            let mut l1d_ecur: InterferenceList = Vec::new();
            fun_select_trasha_interference(
                &mut licur_copy,
                crate::abs::ShapeType::Edge,
                &mut l1d_ecur,
            );

            let GetTr {
                stb,
                isb,
                bdim,
                sta,
                isa,
                adim,
            } = fun_ds_get_tr(i_eds, &licur);

            if is_se {
                // cxx:553-680
                let b_in1d = stb == State::In && bdim == 1;
                let b_in2d = stb == State::In && bdim == 2;
                let b_in3d = stb == State::In && bdim == 3;
                let b_out2d = stb == State::Out && bdim == 2;
                let b_out3d = stb == State::Out && bdim == 3;
                let a_in1d = sta == State::In && adim == 1;
                let a_in2d = sta == State::In && adim == 2;
                let a_in3d = sta == State::In && adim == 3;
                let a_out2d = sta == State::Out && adim == 2;
                let a_out3d = sta == State::Out && adim == 3;

                let mut new_t = Transition::default();
                let mut inb = false;
                let mut ina = false;
                let mut skip = false;
                if isfafa {
                    if isspon {
                        if stb == State::Out && sta == State::Out {
                            skip = true;
                        } else {
                            inb = b_in1d;
                            ina = a_in1d;
                            new_t.set_index(isb);
                            new_t.set_shape_before(crate::abs::ShapeType::Edge);
                            new_t.set_shape_after(crate::abs::ShapeType::Edge);
                        }
                    } else if isspin {
                        inb = b_in2d;
                        ina = a_in2d;
                        new_t.set_shape_before(crate::abs::ShapeType::Face);
                        new_t.set_shape_after(crate::abs::ShapeType::Face);
                    } else if isspou {
                        inb = !b_out2d;
                        ina = !a_out2d;
                        new_t.set_shape_before(crate::abs::ShapeType::Face);
                        new_t.set_shape_after(crate::abs::ShapeType::Face);
                    }
                } else if isspon {
                    if stb == State::Out && sta == State::Out {
                        skip = true;
                    } else {
                        inb = b_in1d || b_in2d;
                        ina = a_in1d || a_in2d;
                        new_t.set_index(isb);
                        new_t.set_shape_before(crate::abs::ShapeType::Edge);
                        new_t.set_shape_after(crate::abs::ShapeType::Edge);
                    }
                } else if isspin {
                    if stb == State::Out && sta == State::Out {
                        skip = true;
                    } else {
                        inb = b_in3d;
                        ina = a_in3d;
                        if inb {
                            new_t.set_index(isb);
                        } else {
                            new_t.set_index(isa);
                        }
                        new_t.set_shape_before(crate::abs::ShapeType::Face);
                        new_t.set_shape_after(crate::abs::ShapeType::Face);
                    }
                } else if isspou {
                    if stb == State::In && sta == State::In {
                        skip = true;
                    } else {
                        inb = !b_out3d;
                        ina = !a_out3d;
                        if b_out3d {
                            new_t.set_index(isb);
                        } else {
                            new_t.set_index(isa);
                        }
                        new_t.set_shape_before(crate::abs::ShapeType::Face);
                        new_t.set_shape_after(crate::abs::ShapeType::Face);
                    }
                }

                if skip {
                    tki.next();
                    continue;
                }

                let sb = if inb { State::In } else { State::Out };
                let sa = if ina { State::In } else { State::Out };
                new_t.set_state_before(sb);
                new_t.set_state_after(sa);

                // `cxx:666-668`: `S` is a dummy, `B` is `GBound` unless the
                // geometry is a point.
                let s = 0;
                let b = if kcur == Kind::Point {
                    false
                } else {
                    i_first
                        .read()
                        .expect("interference lock")
                        .g_bound()
                        .unwrap_or(false)
                };
                if let Some(new_i) = make_epv_interference(&new_t, s, gcur, par_se, kcur, b) {
                    let mut li: InterferenceList = Vec::new();
                    li.push(new_i);
                    let it_cur = PointIterator::new(li);
                    self.fill_point_topology_pvs_it(e, &it_cur, g, pvs);
                }
                tki.next();
                continue;
            }

            // - kp3 - cxx:682-688: `l2dFEcur.First()->Transition().Orientation(IN)`
            // is a const observation in OCCT, so it has no effect.
            let _kp3 = n2d > 0 && point;

            // cxx:689-772
            let mut it_cur = PointIterator::new(licur.clone());
            while it_cur.more() {
                let i1 = it_cur.value().clone();
                let (sb1, ib1, sa1, ia1, _gt1, _g1, st1, _s1) = {
                    let i = i1.read().expect("interference lock");
                    crate::topopebrep::ds::export::fds_idata(&i)
                };
                let _ = (sb1, ib1, sa1, ia1);

                let keepinterf1 = if is_ed {
                    true
                } else if global_faces2d() {
                    // split 2d, `cxx:709-742`
                    let has2d3d = n2d > 0 && n3d > 0;
                    if has2d3d && !is_se {
                        let i3d = l3d_fcur.first();
                        i3d.map(|i3d| Arc::ptr_eq(i3d, &i1)).unwrap_or(false)
                    } else {
                        true
                    }
                } else {
                    // split 3d, `cxx:746-749`
                    st1 == Kind::Face
                };

                if keepinterf1 {
                    if kp6 {
                        let gb1 = {
                            let i = i1.read().expect("interference lock");
                            fun_isonbound(hds, &i)
                        };
                        if !gb1 {
                            self.fill_point_topology_pvs_it(e, &it_cur, g, pvs);
                        }
                    } else {
                        self.fill_point_topology_pvs_it(e, &it_cur, g, pvs);
                    }
                    if !dg_e {
                        break;
                    }
                }
                it_cur.next();
            }

            tki.next();
        }
    }

    /// `GFillPointTopologyPVS(E, EPit, G1, PVS)` (`GridEE.cxx:778-1051`): add the
    /// point or vertex the current interference carries to `PVS`.
    pub fn fill_point_topology_pvs_it(
        &self,
        e: &TopoShape,
        epit: &PointIterator,
        g1: &GTopo,
        pvs: &mut PaveSet,
    ) {
        let Some(hds) = self.data_structure() else {
            return;
        };
        let epvs = pvs.edge().clone();

        let i1 = epit.value().clone();
        let st1 = i1.read().expect("interference lock").support_type();

        let (tb1, _tb2) = g1.states_on();
        let tb = tb1;

        let ig = epit.geometry();
        let ispoint = epit.is_point();
        let vig = if ispoint {
            self.new_vertex(ig).cloned()
        } else {
            hds.shape(ig, true)
        };
        let Some(vig) = vig else {
            return; // `cxx:818-822`: PMN 17/02/99 Nothing to add.
        };

        let mut has_vsd = false;
        let mut i_vre = 0;
        let mut vre: Option<TopoShape> = None;
        if !ispoint {
            has_vsd = hds.has_same_domain(&vig, true);
            if has_vsd {
                i_vre = hds.same_domain_reference(&vig);
                vre = hds.shape(i_vre, true);
            }
        }

        // `cxx:834-844`
        let Some(vpv) = (if has_vsd { vre.clone() } else { Some(vig.clone()) }) else {
            return;
        };

        let par = epit.parameter().unwrap_or(0.0);
        let mut ori = epit.orientation(tb);

        let samegeom = fun_aresamegeom(e, &epvs);
        if is_diff_oriented(g1) {
            ori = ori.complement();
        }
        if !samegeom {
            ori = ori.complement();
        }

        let mut par = par;
        let lesmemes = self
            .my_edge_reference
            .as_ref()
            .map(|r| is_equal(e, r))
            .unwrap_or(false);
        if !lesmemes {
            let mut parref = par;
            self.param_on_reference(&vpv, &mut parref);
            par = parref;
        }

        // `cxx:891-919`
        let kpbound = {
            let edge_ref = self.my_edge_reference.clone();
            let mut kpbound = false;
            if let Some(edge_ref) = edge_ref.as_ref() {
                let eclosed = closed_e(edge_ref);
                let ori_i = epit.orientation(State::In);
                let intex = (ori_i == Orientation::Internal && tb == State::In)
                    || (ori_i == Orientation::External && tb == State::Out);
                kpbound = lesmemes && eclosed.is_some() && intex;
                if kpbound {
                    let vclo = eclosed.expect("closed edge vertex");
                    kpbound = is_same(&vclo, &vig);
                    if !kpbound {
                        if let Some(vsd) = fun_ds_getoov(hds.ds(), &vig) {
                            kpbound = is_same(&vclo, &vsd);
                        }
                    }
                }
            }
            kpbound
        };

        let same_dom = if has_vsd {
            first_same_domain(hds, &vpv).map(|vsd| {
                // `cxx:929-938`: the DS's same-domain shape itself, unless it is
                // the reference, in which case the interfering vertex is used.
                let i_vsd = hds.shape_index(&vsd, true);
                if i_vsd == i_vre {
                    vig.clone()
                } else {
                    vsd
                }
            })
        } else {
            None
        };

        if !kpbound {
            let mut vpv = vpv;
            vpv.set_orientation(ori);
            append_pave(pvs, &vpv, par, has_vsd, same_dom.as_ref(), st1);
        } else {
            // `cxx:958-1051`
            let Some(edge_ref) = self.my_edge_reference.clone() else {
                return;
            };
            let (parf, parl) = fun_tool_bounds(&edge_ref);

            let mut vpv = vpv;
            vpv.set_orientation(Orientation::Forward);
            append_pave(pvs, &vpv, parf, has_vsd, same_dom.as_ref(), st1);

            vpv.set_orientation(Orientation::Reversed);
            append_pave(pvs, &vpv, parl, has_vsd, same_dom.as_ref(), st1);

            pvs.set_remove_pv(false); // jyl + 980217
        }
    }
}

/// `new TopOpeBRepBuild_Pave(VPV, PAR, bound=false)` followed by
/// `HasSameDomain` / `SameDomain` / `InterferenceType() = ST1` and
/// `PVS.Append(PV)` (`GridEE.cxx:924-941`, `:975-992`, `:1013-1031`).
fn append_pave(
    pvs: &mut PaveSet,
    vpv: &TopoShape,
    par: f64,
    has_vsd: bool,
    same_dom: Option<&TopoShape>,
    st1: Kind,
) {
    let mut l = Loop::new_pave(vpv, par, false);
    if let Some(p) = l.as_pave_mut() {
        if has_vsd {
            p.set_has_same_domain(true);
            if let Some(vsd) = same_dom {
                p.set_same_domain(vsd);
            }
        }
        // `cxx:942`, `:990`, `:1030`: `PV->InterferenceType() = ST1`.
        *p.interference_type_mut() = st1;
    }
    pvs.append(Arc::new(RwLock::new(l)));
}
