//! `TopOpeBRepBuild_Builder1` -- the 2d piece helpers.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_Builder1.cxx`
//! (OCCT 8.0.0):
//!
//! | OCCT | lines |
//! |---|---|
//! | `PerformPieceIn2D(EdgeToPerform, EOR, edgeFace, toFace, G1, keep)` | `1415-1527` |
//! | `OrientateEdgeOnFace(EdgeToPerform, baseFace, edgeFace, G1, stateOfFaceOri)` | `2227-2290` |
//! | `ClassifyEdgeToFaceByOnePoint(E, F)` (file-static) | `2296-2316` |
//!
//! `keep` is an out-parameter in OCCT; [`Builder1::perform_piece_in_2d`]
//! returns it and [`Builder1::orientate_edge_on_face`] returns
//! `stateOfFaceOri` for the same reason.
//!
//! UNPORTED: the commented-back-out lines OCCT keeps beside the live ones
//! (`cxx:1463`, `:1478`, `:1486-1487`, `:1596`, `:1614-1615`, `:1655-1656`,
//! `:2024-2025`, `:2175-2178`).

use occt_core::gp::GpVec;

use crate::abs::{Orientation, State};
use crate::fclass2d::{FClass2d, FaceState};
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::is_same;
use crate::topopebrep::build::g_topo::GTopo;
use crate::topopebrep::build::tools::{
    fc2d_curve_on_surface, get_normal_in_nearest_point, get_normal_to_face_on_edge,
    get_tangent_to_edge,
};

use super::super::Builder;

/// A zero `gp_Vec`, the value OCCT keeps in an out-vector whose query failed
/// (`GetNormalToFaceOnEdge` answers `false` and leaves the vector untouched).
fn zero_vec() -> GpVec {
    GpVec::new(0.0, 0.0, 0.0)
}

/// `ClassifyEdgeToFaceByOnePoint(E, F)` (`cxx:2296-2316`): the state of the
/// p-curve point at `f * PAR_T + (1 - PAR_T) * l` inside the face's 2d region.
///
/// OCCT passes the literal `1e-7` to `BRepTopAdaptor_FClass2d` (not
/// `Precision::PConfusion()`), and a null p-curve answers `UNKNOWN`.
pub(crate) fn classify_edge_to_face_by_one_point(e: &Edge, f: &Face) -> State {
    let Some((c2d, f2, l2)) = fc2d_curve_on_surface(e, f) else {
        return State::Unknown;
    };

    let par = f2 * super::PAR_T + (1.0 - super::PAR_T) * l2;
    let a_p2d = c2d.d0(par);

    let Ok(fc) = FClass2d::new(f, 1.0e-7) else {
        return State::Unknown;
    };
    match fc.perform(a_p2d) {
        FaceState::In => State::In,
        FaceState::Out => State::Out,
        FaceState::On => State::On,
        FaceState::Unknown => State::Unknown,
    }
}

impl Builder {
    /// `PerformPieceIn2D(EdgeToPerform, EOR, edgeFace, toFace, G1, keep)`
    /// (`cxx:1415-1527`).
    ///
    /// `myMapOfEdgeFaces` is filled by `GFillSolidSFS` (`cxx:201-202`); the
    /// list holds the faces of *both* operands, which is why the adjacent face
    /// is looked up here instead of being derived from `edgeFace`.
    pub(crate) fn perform_piece_in_2d(
        &self,
        edge_to_perform: &Edge,
        eor: &Edge,
        edge_face: &Face,
        to_face: &Face,
        g1: &GTopo,
    ) -> bool {
        let (tb1, tb2) = g1.states_on();
        let iref = self
            .data_structure()
            .map(|hds| hds.ds().ancestor_rank_of(&eor.0))
            .unwrap_or(0);
        let tb = if iref == 1 { tb1 } else { tb2 };

        let o1 = edge_face.0.orientation();
        let o2 = to_face.0.orientation();
        let ori_e = edge_to_perform.0.orientation();

        let mut a_n2 = get_normal_to_face_on_edge(to_face, edge_to_perform).unwrap_or_else(zero_vec);
        if o2 == Orientation::Reversed {
            a_n2.reverse();
        }

        let mut a_n3 =
            get_normal_to_face_on_edge(edge_face, edge_to_perform).unwrap_or_else(zero_vec);
        if o1 == Orientation::Reversed {
            a_n3.reverse();
        }

        let mut a_tg = get_tangent_to_edge(edge_to_perform).unwrap_or_else(zero_vec);
        if ori_e == Orientation::Reversed {
            a_tg.reverse();
        }
        if o1 == Orientation::Reversed {
            a_tg.reverse();
        }

        let a_bi_n = a_tg.crossed(&a_n2);

        let a_fel = self
            .b1
            .my_map_of_edge_faces
            .find_from_key(&eor.0)
            .cloned()
            .unwrap_or_default();

        let mut scalar_pr = 0.0;
        let mut keep = false;

        // `cxx:1465-1508`: with at most two adjacent faces the bi-normal decides.
        if a_fel.len() <= 2 {
            for a_adj in a_fel.iter() {
                if is_same(&edge_face.0, a_adj) {
                    continue;
                }
                let Some(a_adj_f) = Face::wrap(a_adj.clone()) else {
                    continue;
                };
                let mut a_n1 = get_normal_to_face_on_edge(&a_adj_f, edge_to_perform)
                    .unwrap_or_else(zero_vec);
                if a_adj_f.0.orientation() == Orientation::Reversed {
                    a_n1.reverse();
                }
                scalar_pr = a_bi_n.dot(&a_n1);

                if scalar_pr.abs() <= 1.0e-10 {
                    // `cxx:1484-1498`: try to step inside the adjacent face.
                    if let Some(mut n) = get_normal_in_nearest_point(&a_adj_f, edge_to_perform) {
                        if a_adj_f.0.orientation() == Orientation::Reversed {
                            n.reverse();
                        }
                        a_n1 = n;
                    }
                    scalar_pr = a_bi_n.dot(&a_n1);
                    if scalar_pr.abs() <= 1.0e-10 {
                        continue;
                    }
                }

                let a_part_state = if scalar_pr > 0.0 {
                    State::In
                } else {
                    State::Out
                };
                keep = a_part_state == tb;
                if keep {
                    break;
                }
            }
        }

        // `cxx:1510-1525`: no adjacent face to compare with, fall back to the
        // two face normals and the operation.
        if scalar_pr.abs() <= 1.0e-10 {
            let n = a_n3.dot(&a_n2);
            if self.opefus() {
                keep = n > 0.0;
            }
            if self.opec12() || self.opec21() {
                keep = n < 0.0;
            }
            if self.opecom() {
                keep = n > 0.0;
            }
        }

        keep
    }

    /// `OrientateEdgeOnFace(EdgeToPerform, baseFace, edgeFace, G1,
    /// stateOfFaceOri)` (`cxx:2227-2290`).
    ///
    /// The return value is OCCT's `stateOfFaceOri` out-parameter: the faces'
    /// normals disagree, which the caller records in
    /// `myMapOfEdgeWithFaceState` and later uses to orient the built faces.
    pub(crate) fn orientate_edge_on_face(
        &self,
        edge_to_perform: &mut TopoShape,
        base_face: &Face,
        edge_face: &Face,
        g1: &GTopo,
    ) -> bool {
        let mut state_of_face_ori = false;

        // `cxx:2237-2247`: the rank comes from `mySDFaceToFill`, not from
        // `edgeFace`, so the reverse flag follows the face being filled.
        let Some(sd_face_to_fill) = self.b1.my_sd_face_to_fill.as_ref() else {
            return state_of_face_ori;
        };
        let curr_ref = self
            .data_structure()
            .map(|hds| hds.ds().ancestor_rank_of(sd_face_to_fill))
            .unwrap_or(0);
        let face_ref = self
            .data_structure()
            .map(|hds| hds.ds().ancestor_rank_of(&edge_face.0))
            .unwrap_or(0);
        let rev_ori = if curr_ref == 1 {
            g1.is_to_reverse1()
        } else {
            g1.is_to_reverse2()
        };

        let ori_e = edge_to_perform.orientation();
        let newori_e = super::super::Builder::orient(ori_e, rev_ori);
        let face_ori = edge_face.0.orientation();
        let base_ori = base_face.0.orientation();
        let curr_ori = sd_face_to_fill.orientation();

        let a_e = Edge(edge_to_perform.clone());
        let mut a_n1 = get_normal_to_face_on_edge(edge_face, &a_e).unwrap_or_else(zero_vec);
        if face_ori == Orientation::Reversed {
            a_n1.reverse();
        }

        let mut a_n2 = get_normal_to_face_on_edge(base_face, &a_e).unwrap_or_else(zero_vec);
        if base_ori == Orientation::Reversed {
            a_n2.reverse();
        }

        if a_n1.dot(&a_n2) < 0.0 {
            state_of_face_ori = true;
        }

        // `cxx:2273-2279`: for the tool operand the face flips in a cut.
        if face_ref == 2 && (self.opec12() || self.opec21()) {
            state_of_face_ori = !state_of_face_ori;
        }

        edge_to_perform.set_orientation(newori_e);

        if curr_ori != base_ori {
            edge_to_perform.reverse();
        }
        if state_of_face_ori {
            edge_to_perform.reverse();
        }

        state_of_face_ori
    }
}
