//! `Draft_Modification::Perform`, last block: "small loop of
//! validation/protection" (`Draft_Modification_1.cxx:1640-1774`).
//!
//! After the vertices are rebuilt, every modified edge is checked once more:
//! the direction of its new curve is compared with the source edge's, the edge
//! is reversed when they disagree, and the vertex parameters are read back
//! through `NewParameter`, with the closed-curve wrap at the end.
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | edge loop, `TopExp::Vertices`, `REVERSED` swap | `_1.cxx:1642-1652` | [`DraftModification::validate_protection`] |
//! | parameter order test, direction comparison | `_1.cxx:1654-1706` | [`DraftModification::validate_protection`] |
//! | reverse the new curve, re-parameterise | `_1.cxx:1707-1721` | [`DraftModification::validate_protection`] |
//! | `NewParameter` read back, `pl <= pf` | `_1.cxx:1724-1738` | [`DraftModification::validate_protection`] |
//! | closed-curve period fold | `_1.cxx:1739-1759` | [`DraftModification::validate_protection`] |
//! | final failure, parameter store | `_1.cxx:1760-1773` | [`DraftModification::validate_protection`] |
//!
//! # Translation notes
//!
//! * `Vf.Reverse()` / `Vl.Reverse()` (`cxx:1649-1650`) reverse the *orientation*
//!   of the vertex occurrences, which is what `BRep_Tool::Parameter(V, E)`
//!   reads; the port builds the reversed occurrence with
//!   `TopoShape::oriented`.
//! * `TopExp::Vertices(edg, Vf, Vl)` (`cxx:1646`) is
//!   [`crate::topo_tools::edge_vertices`].
//! * `anIntCurv->Reverse()` (`cxx:1708`) mutates the shared curve; the port
//!   stores the reversed copy in the edge record and re-parameterizes against
//!   that copy.
//! * `NewParameter(V, E, P, Tol)` (`cxx:1726`, `:1730`) is
//!   [`DraftModification::new_parameter_of`]; its `Err` is the reference's
//!   `Standard_NoSuchObject` throw and is re-raised here as a panic.
//! * `BRep_Tool::Parameter(V, E)` is [`BRepTool::parameter_of_vertex`].

use std::sync::Arc;

use occt_core::gp::GpVec;
use occt_core::precision::{Precision, PCONFUSION};
use occt_geom::curve::Curve;

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Vertex};
use crate::topo_tools::edge_vertices;

use super::parameters::parameter;
use super::DraftModification;

/// `gp_Vec::Normalize` guarded by `SquareMagnitude > SquareConfusion`
/// (`cxx:1678-1701`): OCCT divides in place only when the magnitude is
/// meaningful, leaving a null vector as it is.
fn normalize_if_not_null(v: &mut GpVec) {
    let sq = v.square_magnitude();
    if sq > Precision::SQUARE_CONFUSION {
        let n = v.divided(sq.sqrt());
        *v = n;
    }
}

impl DraftModification {
    /// `errStat = Draft_EdgeRecomputation; badShape = edg; return;`
    /// (`cxx:1761-1763`).
    fn fail_edge_unit(&mut self, e: &Edge) -> Result<(), ()> {
        self.set_edge_failure(e);
        Err(())
    }

    /// `Draft_Modification_1.cxx:1640-1774`.
    pub(super) fn validate_protection(&mut self) -> Result<(), ()> {
        for i in 1..=self.my_e_map.extent() {
            let edg_shape = self.my_e_map.find_key(i).clone();
            let Some(edg) = Edge::wrap(edg_shape.clone()) else {
                continue;
            };

            // `cxx:1645-1652`
            let (Some(vf), Some(vl)) = edge_vertices(&edg) else {
                continue;
            };
            let (vf, vl) = if edg.0.orientation() == Orientation::Reversed {
                (
                    Vertex::wrap(vf.0.oriented(Orientation::Reversed))
                        .expect("Standard_TypeMismatch: TopoDS::Vertex"),
                    Vertex::wrap(vl.0.oriented(Orientation::Reversed))
                        .expect("Standard_TypeMismatch: TopoDS::Vertex"),
                )
            } else {
                (vf, vl)
            };

            // `cxx:1654-1658`
            if self.my_v_map.contains(&vf.0) && self.my_v_map.contains(&vl.0) {
                let apar_f = self
                    .my_v_map
                    .find_from_key(&vf.0)
                    .expect("Standard_NoSuchObject: Draft_Modification::Perform")
                    .parameter(&edg);
                let apar_l = self
                    .my_v_map
                    .find_from_key(&vl.0)
                    .expect("Standard_NoSuchObject: Draft_Modification::Perform")
                    .parameter(&edg);

                // `cxx:1660-1721`: compare the source edge's direction with the
                // selected part of the intersection edge.
                if apar_l < apar_f {
                    let Some(s_curve) = BRepTool::edge_curve(&edg) else {
                        continue;
                    };
                    let Some(int_curve) = self
                        .my_e_map
                        .find_from_key(&edg.0)
                        .and_then(|e| e.geometry())
                        .cloned()
                    else {
                        continue;
                    };
                    let pf_param = BRepTool::parameter_of_vertex(&edg.0, &vf.0);
                    let pl_param = BRepTool::parameter_of_vertex(&edg.0, &vl.0);
                    let (a_pf, mut dir_of) = s_curve.d1(pf_param);
                    let (a_pl, mut dir_ol) = s_curve.d1(pl_param);
                    let (_, mut dir_nf) = int_curve.d1(apar_f);
                    let (_, mut dir_nl) = int_curve.d1(apar_l);

                    normalize_if_not_null(&mut dir_nf);
                    normalize_if_not_null(&mut dir_nl);
                    normalize_if_not_null(&mut dir_of);
                    normalize_if_not_null(&mut dir_ol);

                    let a_cos_f = dir_nf.dot(&dir_of);
                    let a_cos_l = dir_nl.dot(&dir_ol);
                    let a_cos_max = if a_cos_f.abs() > a_cos_l.abs() {
                        a_cos_f
                    } else {
                        a_cos_l
                    };

                    if a_cos_max < 0.0 {
                        // `cxx:1707-1712`
                        let reversed: Arc<dyn Curve> = Arc::from(int_curve.reversed());
                        self.my_e_map
                            .change_from_key(&edg.0)
                            .expect("Standard_NoSuchObject: Draft_Modification::Perform")
                            .change_geometry()
                            .replace(reversed.clone());

                        // `cxx:1711-1720`
                        let mut an_err = 0i32;
                        let a_par = parameter(&reversed, &a_pf, &mut an_err);
                        if an_err == 0 {
                            *self
                                .my_v_map
                                .change_from_key(&vf.0)
                                .expect("Standard_NoSuchObject: Draft_Modification::Perform")
                                .change_parameter(&edg) = a_par;
                        }
                        let mut an_err = 0i32;
                        let a_par = parameter(&reversed, &a_pl, &mut an_err);
                        if an_err == 0 {
                            *self
                                .my_v_map
                                .change_from_key(&vl.0)
                                .expect("Standard_NoSuchObject: Draft_Modification::Perform")
                                .change_parameter(&edg) = a_par;
                        }
                    }
                }
            }

            // `cxx:1724-1732`
            let mut tolerance = 0.0f64;
            let mut p = 0.0f64;
            let mut pf = match self.new_parameter_of(&vf, &edg, &mut p, &mut tolerance) {
                Ok(true) => p,
                Ok(false) => BRepTool::parameter_of_vertex(&edg.0, &vf.0),
                Err(err) => panic!("{err}"),
            };
            let mut tolerance = 0.0f64;
            let mut p = 0.0f64;
            let mut pl = match self.new_parameter_of(&vl, &edg, &mut p, &mut tolerance) {
                Ok(true) => p,
                Ok(false) => BRepTool::parameter_of_vertex(&edg.0, &vl.0),
                Err(err) => panic!("{err}"),
            };

            // `cxx:1733-1765`
            if pl <= pf {
                let the_curve = self
                    .my_e_map
                    .find_from_key(&edg.0)
                    .and_then(|e| e.geometry())
                    .cloned();
                let Some(the_curve) = the_curve else {
                    return self.fail_edge_unit(&edg);
                };
                if the_curve.is_closed() {
                    let first_par = the_curve.first_parameter();
                    let last_par = the_curve.last_parameter();
                    let pconf = PCONFUSION;
                    if (pf - last_par).abs() <= pconf {
                        pf = first_par;
                    } else if (pl - first_par).abs() <= pconf {
                        pl = last_par;
                    }
                    if pl <= pf {
                        pl += last_par - first_par;
                    }
                }
                if pl <= pf {
                    // `cxx:1760-1764`
                    return self.fail_edge_unit(&edg);
                }
            }

            // `cxx:1766-1772`
            if self.my_v_map.contains(&vf.0) {
                *self
                    .my_v_map
                    .change_from_key(&vf.0)
                    .expect("Standard_NoSuchObject: Draft_Modification::Perform")
                    .change_parameter(&edg) = pf;
            }
            if self.my_v_map.contains(&vl.0) {
                *self
                    .my_v_map
                    .change_from_key(&vl.0)
                    .expect("Standard_NoSuchObject: Draft_Modification::Perform")
                    .change_parameter(&edg) = pl;
            }
        }
        Ok(())
    }
}
