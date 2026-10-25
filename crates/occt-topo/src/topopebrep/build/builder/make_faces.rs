//! `TopOpeBRepBuild_Builder::GFABUMakeFaces` — rebuild faces out of the wires
//! the `FaceBuilder` classified.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_makefaces.cxx`
//! (OCCT 8.0.0):
//!
//! | OCCT | lines |
//! |---|---|
//! | `FUN_purgeFon1nonoriE` | `248-276` |
//! | `FUN_ReOrientIntExtEdge` | `284-333` |
//! | `FUN_CheckORI` | `337-360` |
//! | `GFABUMakeFaces` | `365-755` |
//!
//! The wire of a new face is assembled from the block elements of one loop
//! (`cxx:441-479`). When that wire mixes oriented edges with `INTERNAL` /
//! `EXTERNAL` ones, OCCT walks the vertex-to-edge map of the collected
//! compound and pairs the two edges at each vertex, re-orienting the
//! `INTERNAL` (`EXTERNAL`) edge from the oriented one (`cxx:500-648`). That
//! re-orientation is what keeps the rebuilt face regular.
//!
//! The port keeps `TopoDS_Shape` occurrences by value, so `E1`/`E2` below are
//! clones of the compound elements: re-orienting them does not alias the
//! `FaceBuilder` block elements the way the C++ references do. Every read that
//! follows a re-orientation in OCCT reads the clone here.
//!
//! UNPORTED:
//! - the `Geom_RectangularTrimmedSurface` window widening (`cxx:706-733`):
//!   `hrts->SetTrim` mutates the surface shared by every new face of this
//!   call, and the port keeps surfaces immutable behind `Arc<dyn Surface>`
//!   (`GeomRectangularTrimmedSurface::set_trim_uv` needs `&mut self`). No
//!   branch is substituted.
//! - the `OCCT_DEBUG` traces (`cxx:370-381`, `:652-655`, `:739-748`).
//!
//! The `GF2D_CurveOnSurface(newEdge, newFace, f2, l2, tolpc, theUFront = true)`
//! overload (`cxx:471`) is reached through the port's get-or-build
//! [`crate::boptools_2d::make_2d`], which has no `theUFront` switch; the
//! tolerance `max(tolE, tolpc)` of `BRep_Builder::UpdateEdge` is not stored,
//! because the port's p-curve side table keeps no per-curve tolerance.

use std::collections::HashMap;

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Compound, Edge, Face, TopoShape, Wire};
use crate::topo_tools::{explorer_shapes, is_same, shape_map_key, ShapeMapKey};
use crate::topopebrep::build::face_builder::map_shapes_and_ancestors_ordered;
use crate::topopebrep::tool::purge::fun_tool_closed_w;

use super::super::face_builder::FaceBuilder;
use super::Builder;

impl Builder {
    /// `GFABUMakeFaces(FF, FABU, LOF, MWisOld)` (`cxx:365-755`).
    pub fn g_fabumake_faces(
        &self,
        ff: &TopoShape,
        fabu: &mut FaceBuilder,
        lof: &mut Vec<TopoShape>,
        m_wis_old: &mut HashMap<ShapeMapKey, i32>,
    ) {
        // `cxx:386-395`: `Surf = BRep_Tool::Surface(FF, Loc)`
        let surf = BRepTool::face_surface_world(&Face(ff.clone()));
        // `cxx:398`
        let tol_ff = BRepTool::face_tolerance(&Face(ff.clone()));
        let td_f = ff.clone();

        // `cxx:413-414`
        fabu.init_face();
        while fabu.more_face() {
            let mut nb_new_w_withe: i32 = 0;
            let mut nb_old_w: i32 = 0;

            // `cxx:420`: `BB.MakeFace(newFace, Surf, Loc, tolFF)`
            let mut new_face = match &surf {
                Some(su) => TopoBuilder::new().make_face_with_tolerance(su.clone(), &[], tol_ff).0,
                None => Face::new().0,
            };

            // `cxx:423`
            let nbw = fabu.init_wire();
            while fabu.more_wire() {
                let mut ne: i32 = 0;
                let mut ne_internal: i32 = 0;
                let mut ne_external: i32 = 0;

                // `cxx:430-436`
                let is_old = fabu.is_old_wire();
                let new_wire: TopoShape;
                if is_old {
                    nb_old_w += 1;
                    new_wire = match fabu.old_wire() {
                        Some(w) => w,
                        None => Wire::new().0,
                    };
                } else {
                    // `cxx:438-439`
                    let mut cmp_of_edges = Compound::new().0;
                    let mut wire = Wire::new().0;
                    self.my_build_tool.make_wire(&mut wire);

                    // `cxx:440-479`
                    fabu.init_edge();
                    while fabu.more_edge() {
                        let new_edge = fabu.edge().clone();
                        // `cxx:448`
                        let nadde = fabu.add_edge_wire(&new_edge, &mut cmp_of_edges);
                        ne += nadde;

                        // `cxx:455-462`
                        let oe = new_edge.orientation();
                        if oe == Orientation::Internal {
                            ne_internal += 1;
                        } else if oe == Orientation::External {
                            ne_external += 1;
                        }

                        // `cxx:464-478`
                        if !boptools_2d::has_curve_on_surface(
                            &Edge(new_edge.clone()),
                            &Face(new_face.clone()),
                        ) {
                            // `cxx:472` throws `Standard_ProgramError` on a null
                            // p-curve; the edge is left without one here.
                            if let Ok(c2d) = boptools_2d::make_2d(
                                &Edge(new_edge.clone()),
                                &Face(new_face.clone()),
                            ) {
                                let _ = boptools_2d::attach_existing_pcurve(
                                    &Edge(new_edge.clone()),
                                    &Face(new_face.clone()),
                                    c2d,
                                );
                            }
                        }
                        fabu.next_edge();
                    }

                    // `cxx:481-648`
                    let plain = (ne_internal == 0 && ne_external == 0)
                        || (ne == ne_internal || ne == ne_external);
                    if plain {
                        // `cxx:484-496`
                        let mut bailed = false;
                        if nbw == 1 && ne == 2 {
                            let edges = explorer_shapes(&cmp_of_edges, ShapeType::Edge);
                            if edges.len() >= 2 && is_same(&edges[0], &edges[1]) {
                                bailed = true;
                            }
                        }
                        if bailed {
                            // `cxx:492`: leaves `GFABUMakeFaces`, dropping every
                            // face collected so far.
                            return;
                        }
                        // `cxx:497-500`
                        for new_edge in explorer_shapes(&cmp_of_edges, ShapeType::Edge) {
                            fabu.add_edge_wire(&new_edge, &mut wire);
                        }
                    } else {
                        // `cxx:504-506`
                        let map_voe = map_shapes_and_ancestors_ordered(
                            &cmp_of_edges,
                            ShapeType::Vertex,
                            ShapeType::Edge,
                        );

                        // `cxx:509-518`: every vertex must carry exactly two
                        // edges for the wire to be regular.
                        let mut w_is_clsd = true;
                        for i in 1..=map_voe.extent() {
                            let n = map_voe
                                .find_from_index(i)
                                .map(|l| l.len())
                                .unwrap_or(0);
                            if n != 2 {
                                w_is_clsd = false;
                                break;
                            }
                        }

                        if !w_is_clsd {
                            // `cxx:519-528`
                            for new_edge in explorer_shapes(&cmp_of_edges, ShapeType::Edge) {
                                fabu.add_edge_wire(&new_edge, &mut wire);
                            }
                        } else {
                            // `cxx:529-648`
                            let mut lof_add_e: Vec<TopoShape> = Vec::new();
                            let mut naddsame: i32 = 0;
                            while ne > (lof_add_e.len() as i32 + naddsame) {
                                for step_map in 1..=map_voe.extent() {
                                    let Some(lof_e) = map_voe.find_from_index(step_map) else {
                                        continue;
                                    };
                                    if lof_e.len() < 2 {
                                        continue;
                                    }
                                    // `cxx:540-547`
                                    let mut e1 = lof_e[0].clone();
                                    let mut e2 = lof_e[1].clone();
                                    let mut o1 = e1.orientation();
                                    let mut o2 = e2.orientation();
                                    let is_same_e1 = BRepTool::is_closed_edge_face(
                                        &Edge(e1.clone()),
                                        &Face(td_f.clone()),
                                    );
                                    let is_same_e2 = BRepTool::is_closed_edge_face(
                                        &Edge(e2.clone()),
                                        &Face(td_f.clone()),
                                    );
                                    let mut add_e1 = true;
                                    let mut add_e2 = true;

                                    // `cxx:549-568`
                                    for le in &lof_add_e {
                                        let ole = le.orientation();
                                        if is_same(&e1, le) && !is_same_e1 {
                                            add_e1 = false;
                                            e1.set_orientation(ole);
                                            o1 = ole;
                                        }
                                        if is_same(&e2, le) && !is_same_e2 {
                                            add_e2 = false;
                                            e2.set_orientation(ole);
                                            o2 = ole;
                                        }
                                    }

                                    // `cxx:569-620`
                                    let chk_ori = fun_check_ori(o1, o2);
                                    if chk_ori == 0 {
                                        add_e1 = false;
                                        add_e2 = false;
                                    }
                                    if chk_ori == 1 {
                                        let ori = fun_reorient_int_ext_edge(&e2, o2, &e1);
                                        if ori == Orientation::Forward {
                                            e1.set_orientation(Orientation::Forward);
                                        }
                                        if ori == Orientation::Reversed {
                                            e1.set_orientation(Orientation::Reversed);
                                        }
                                        if ori == Orientation::Reversed
                                            || ori == Orientation::Forward
                                        {
                                            if o1 == Orientation::Internal {
                                                ne_internal -= 1;
                                            } else {
                                                ne_external -= 1;
                                            }
                                        }
                                    }
                                    if chk_ori == 2 {
                                        let ori = fun_reorient_int_ext_edge(&e1, o1, &e2);
                                        if ori == Orientation::Forward {
                                            e2.set_orientation(Orientation::Forward);
                                        }
                                        if ori == Orientation::Reversed {
                                            e2.set_orientation(Orientation::Reversed);
                                        }
                                        if ori == Orientation::Reversed
                                            || ori == Orientation::Forward
                                        {
                                            if o2 == Orientation::Internal {
                                                ne_internal -= 1;
                                            } else {
                                                ne_external -= 1;
                                            }
                                        }
                                    }

                                    // `cxx:622-647`
                                    if add_e1 {
                                        fabu.add_edge_wire(&e1, &mut wire);
                                        if !is_same_e1 {
                                            lof_add_e.push(e1.clone());
                                        } else {
                                            naddsame += 1;
                                        }
                                    }
                                    if add_e2 {
                                        fabu.add_edge_wire(&e2, &mut wire);
                                        if !is_same_e2 {
                                            lof_add_e.push(e2.clone());
                                        } else {
                                            naddsame += 1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    new_wire = wire;
                }

                // `cxx:659-663`
                if ne != 0 {
                    let iow = if is_old { 1 } else { 0 };
                    m_wis_old.insert(shape_map_key(&new_wire), iow);
                }

                // `cxx:666-680`
                if !is_old {
                    if ne == 0 {
                        fabu.next_wire();
                        continue;
                    } else if nbw == 1 && ne == ne_internal + ne_external {
                        fabu.next_wire();
                        continue;
                    } else {
                        nb_new_w_withe += 1;
                    }
                }

                // `cxx:683-687`
                if !is_old {
                    let closed = fun_tool_closed_w(&new_wire);
                    self.my_build_tool.closed(&new_wire, closed);
                }

                // `cxx:689`
                self.my_build_tool.add_face_wire(&mut new_face, &new_wire);

                fabu.next_wire();
            }

            // `cxx:692-695`
            if nb_new_w_withe == 0 && nb_old_w == 0 {
                fabu.next_face();
                continue;
            }

            // `cxx:697-701`
            if fun_purge_fon1_non_ori_e(&new_face) {
                fabu.next_face();
                continue;
            }

            // `cxx:706-733` UNPORTED: the trimmed-surface window widening.

            // `cxx:734`
            lof.push(new_face);
            fabu.next_face();
        }
    }
}

/// `FUN_purgeFon1nonoriE(newFace)` (`cxx:248-276`): a new face built out of a
/// single edge that carries neither `FORWARD`/`REVERSED` orientation nor a
/// degenerated geometry is discarded.
fn fun_purge_fon1_non_ori_e(new_face: &TopoShape) -> bool {
    let edges = explorer_shapes(new_face, ShapeType::Edge);
    if edges.len() != 1 {
        return false;
    }
    let ed = &edges[0];
    let ori = ed.orientation();
    let has_ori = ori == Orientation::Forward || ori == Orientation::Reversed;
    if !has_ori {
        return true;
    }
    // `cxx:267-271` (j gv, 6.06.02, OCC424)
    BRepTool::is_degenerated(&Edge(ed.clone()))
}

/// `FUN_ReOrientIntExtEdge(FRE, OFRE, INE)` (`cxx:284-333`): the orientation an
/// `INTERNAL`/`EXTERNAL` edge must take so that it continues the oriented edge
/// `fre` around the shared vertex.
fn fun_reorient_int_ext_edge(fre: &TopoShape, ofre: Orientation, ine: &TopoShape) -> Orientation {
    // `cxx:288`
    let mut result = ine.orientation();
    // `cxx:291-292`
    let (vf1, vl1) = crate::topo_tools::edge_vertices(&Edge(fre.clone()));
    let (vf2, vl2) = crate::topo_tools::edge_vertices(&Edge(ine.clone()));

    let same = |a: &Option<crate::shape::Vertex>, b: &Option<crate::shape::Vertex>| -> bool {
        match (a, b) {
            (Some(a), Some(b)) => is_same(&a.0, &b.0),
            _ => false,
        }
    };

    // `cxx:294-331`
    if ofre == Orientation::Forward {
        if same(&vl1, &vf2) {
            result = Orientation::Forward;
        }
        if same(&vl1, &vl2) {
            result = Orientation::Reversed;
        }
        if same(&vf1, &vf2) {
            result = Orientation::Reversed;
        }
        if same(&vf1, &vl2) {
            result = Orientation::Forward;
        }
    }
    if ofre == Orientation::Reversed {
        if same(&vl1, &vf2) {
            result = Orientation::Reversed;
        }
        if same(&vl1, &vl2) {
            result = Orientation::Forward;
        }
        if same(&vf1, &vf2) {
            result = Orientation::Forward;
        }
        if same(&vf1, &vl2) {
            result = Orientation::Reversed;
        }
    }
    result
}

/// `FUN_CheckORI(O1, O2)` (`cxx:337-360`).
fn fun_check_ori(o1: Orientation, o2: Orientation) -> i32 {
    let o1_non_ori = o1 == Orientation::Internal || o1 == Orientation::External;
    let o2_non_ori = o2 == Orientation::Internal || o2 == Orientation::External;
    if o1_non_ori && o2_non_ori {
        0
    } else if o1_non_ori {
        1
    } else if o2_non_ori {
        2
    } else {
        4
    }
}
