//! `TopOpeBRepBuild_Builder1::PerformONParts`.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_Builder1.cxx:1151-1369`
//! (OCCT 8.0.0).
//!
//! Adds to a `WireEdgeSet` the ON splits of an edge that `GFillFaceSameDomWES`
//! could not reach (they are ON `FOR1` but belong to no same-domain face), in
//! the orientation the future face needs. `mySDEdgeMap` (`cxx:51`) filters out
//! the edges already covered by the same-domain faces.
//!
//! UNPORTED: the commented-back-out lines OCCT keeps beside the live ones
//! (`cxx:1329-1330`).

use occt_core::gp::GpVec;

use crate::abs::{Orientation, ShapeType, State};
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{explorer_shapes, is_same};
use crate::topopebrep::build::builder::grid::fun_aresamegeom;
use crate::topopebrep::build::builder::grid_ss_globals::complement;
use crate::topopebrep::build::builder::Builder;
use crate::topopebrep::build::g_topo::GTopo;
use crate::topopebrep::build::tools::{
    get_normal_to_face_on_edge, get_tangent_to_edge,
};
use crate::topopebrep::build::wire_edge_set::WireEdgeSet;
use crate::topopebrep::ds::export::fds_data;
use crate::topopebrep::ds::kind::Kind;

use super::fill_2d::classify_edge_to_face_by_one_point;
use super::sd_edge_map_contains;

/// `GetNormalToFaceOnEdge(F, E, out)`: keeps `out` when the query fails, as
/// OCCT does.
fn normal_on_edge(f: &Face, e: &Edge, out: &mut GpVec) {
    if let Some(n) = get_normal_to_face_on_edge(f, e) {
        *out = n;
    }
}

impl Builder {
    /// `PerformONParts(FOR1, SDFaces, G1, WES)` (`cxx:1151-1369`).
    ///
    /// `SDFaces` is unused by OCCT (the parameter is commented out at the
    /// definition, `cxx:1154`).
    pub(crate) fn perform_on_parts(
        &mut self,
        for1: &TopoShape,
        g1: &GTopo,
        wes: &mut WireEdgeSet,
    ) {
        let (etb1, etb2) = g1.states_on();
        let (ftb1, ftb2) = g1.states_on();

        let iref = self
            .data_structure()
            .map(|hds| hds.ds().ancestor_rank_of(for1))
            .unwrap_or(0);
        let ftb = if iref == 1 { ftb1 } else { ftb2 };

        // `cxx:1175-1184`: the interferences of `FOR1`.
        let list: Vec<_> = self
            .data_structure()
            .map(|hds| hds.ds().shape_interferences(for1, true).to_vec())
            .unwrap_or_default();

        for i_handle in list.iter() {
            // `cxx:1181-1187`: only the shape/shape interferences; the values
            // are copied out so no read lock is held across the walk.
            let (is_shape_shape, gt, gi, st, a_tr) = match i_handle.read() {
                Ok(g) => {
                    let (gt, gi, st, _si) = fds_data(&g);
                    (g.is_shape_shape(), gt, gi, st, g.transition())
                }
                Err(_) => continue,
            };
            if !is_shape_shape {
                continue;
            }
            if gt != Kind::Edge || st != Kind::Face {
                continue;
            }

            let Some(eg) = self.data_structure().and_then(|hds| {
                hds.shape(gi, false)
            }) else {
                continue;
            };

            // `cxx:1198-1202`: edges of the current same-domain faces are done.
            if sd_edge_map_contains(&eg) {
                continue;
            }

            // `cxx:1204-1209`: only edges with ON splits.
            let spl_on: Vec<TopoShape> = self
                .data_structure()
                .and_then(|hds| hds.ds().get_shape_with_state(&eg))
                .map(|sws| sws.part(State::On).to_vec())
                .unwrap_or_default();
            if spl_on.is_empty() {
                continue;
            }

            let a_tr = a_tr;

            let iref_e = self
                .data_structure()
                .map(|hds| hds.ds().ancestor_rank_of(&eg))
                .unwrap_or(0);
            let (etb, rev_ori_e) = if iref_e == 1 {
                (etb1, g1.is_to_reverse1())
            } else {
                (etb2, g1.is_to_reverse2())
            };

            let edge_faces: Vec<TopoShape> = self
                .b1
                .my_map_of_edge_faces
                .find_from_key(&eg)
                .cloned()
                .unwrap_or_default();

            for mut new_e in spl_on.into_iter() {
                let a_sd_shape = for1.clone();
                let a_state = match (Face::wrap(a_sd_shape.clone()), Edge::wrap(new_e.clone())) {
                    (Some(f), Some(e)) => classify_edge_to_face_by_one_point(&e, &f),
                    _ => State::Unknown,
                };
                if a_state != State::In && a_state != State::On {
                    continue;
                }

                let a_sd_face = Face(a_sd_shape.clone());

                // `cxx:1252-1318`: the first edge face whose bi-normal keeps
                // the piece decides.
                let mut keep = false;
                for for_ in edge_faces.iter() {
                    // `cxx:1255-1268`: the edge of `FOR` that is `EG` (`EOR`
                    // keeps the last visited edge when none matches, exactly as
                    // in OCCT, so only an edgeless face is skipped).
                    let mut eor: Option<TopoShape> = None;
                    for e in explorer_shapes(for_, ShapeType::Edge) {
                        let same = is_same(&eg, &e);
                        eor = Some(e);
                        if same {
                            break;
                        }
                    }
                    let Some(eor) = eor else {
                        continue;
                    };

                    let ori_e = eor.orientation();
                    // `cxx:1274`: OCCT assigns `neworiE` here and overwrites it
                    // at `:1356`; the value is never read.
                    let _newori_e = Builder::orient(ori_e, rev_ori_e);

                    new_e.set_orientation(ori_e);

                    let mut a_n2 = GpVec::new(0.0, 0.0, 0.0);
                    normal_on_edge(&Face(for_.clone()), &Edge(new_e.clone()), &mut a_n2);
                    if for_.orientation() == Orientation::Reversed {
                        a_n2.reverse();
                    }

                    let mut a_tg = get_tangent_to_edge(&Edge(new_e.clone()))
                        .unwrap_or_else(|| GpVec::new(0.0, 0.0, 0.0));
                    if new_e.orientation() == Orientation::Reversed {
                        a_tg.reverse();
                    }

                    let mut a_n3 = GpVec::new(0.0, 0.0, 0.0);
                    normal_on_edge(&a_sd_face, &Edge(new_e.clone()), &mut a_n3);
                    if a_sd_face.0.orientation() == Orientation::Reversed {
                        a_n3.reverse();
                    }

                    let mut a_bi_n = a_tg.crossed(&a_n2);
                    let mut scalar_pr = a_bi_n.dot(&a_n3);

                    if scalar_pr.abs() <= 1.0e-10 {
                        // `cxx:1298-1311`: try to step inside `FOR`.
                        normal_on_edge(
                            &Face(for_.clone()),
                            &Edge(new_e.clone()),
                            &mut a_n2,
                        );
                        if for_.orientation() == Orientation::Reversed {
                            a_n2.reverse();
                        }
                        a_bi_n = a_tg.crossed(&a_n2);
                        scalar_pr = a_bi_n.dot(&a_n3);
                        if scalar_pr.abs() <= 1.0e-10 {
                            continue;
                        }
                    }

                    let a_part_state = if scalar_pr > 0.0 {
                        State::In
                    } else {
                        State::Out
                    };
                    keep = a_part_state == etb;
                    if keep {
                        break;
                    }
                }

                // `cxx:1320-1368`: the face orientation of the future face,
                // then the edge orientation.
                if keep {
                    let mut state_of_face_ori = false;
                    let mut a_nsf = GpVec::new(0.0, 0.0, 0.0);
                    normal_on_edge(&a_sd_face, &Edge(new_e.clone()), &mut a_nsf);
                    if a_sd_face.0.orientation() == Orientation::Reversed {
                        a_nsf.reverse();
                    }

                    let mut a_nbf = GpVec::new(0.0, 0.0, 0.0);
                    if let Some(base) = self.b1.my_base_face_to_fill.clone() {
                        if let Some(base_f) = Face::wrap(base.clone()) {
                            normal_on_edge(&base_f, &Edge(new_e.clone()), &mut a_nbf);
                        }
                        if base.orientation() == Orientation::Reversed {
                            a_nbf.reverse();
                        }
                    }

                    if a_nsf.dot(&a_nbf) < 0.0 {
                        state_of_face_ori = true;
                    }

                    let rank_sd = self
                        .data_structure()
                        .map(|hds| hds.ds().ancestor_rank_of(&a_sd_face.0))
                        .unwrap_or(0);
                    if rank_sd == 2 && (self.opec12() || self.opec21()) {
                        state_of_face_ori = !state_of_face_ori;
                    }

                    let mut newori_e = a_tr.orientation(ftb, ShapeType::Shape);
                    let samegeom = self
                        .b1
                        .my_base_face_to_fill
                        .as_ref()
                        .map(|base| fun_aresamegeom(for1, base))
                        .unwrap_or(false);
                    if !samegeom {
                        newori_e = complement(newori_e);
                    }
                    new_e.set_orientation(newori_e);

                    self
                        .b1
                        .my_map_of_edge_with_face_state
                        .insert(crate::topo_tools::shape_map_key(&new_e), state_of_face_ori);
                    wes.add_start_element(&new_e);
                }
            }
        }
    }
}
