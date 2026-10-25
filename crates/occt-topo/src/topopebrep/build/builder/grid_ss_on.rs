//! `TopOpeBRepBuild_GridSS.cxx` — the ON-patch part of the solid/solid grid.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_GridSS.cxx`
//! (OCCT 8.0.0):
//!
//! | OCCT | lines |
//! |---|---|
//! | `AddONPatchesSFS(G1, SFS)` | `1454-1592` |
//!
//! `AddONPatchesSFS` walks `myONFacesMap` (new ON face -> its ancestor face),
//! pairs a rank-1 face with a rank-2 face that shares its boundaries and has
//! the same normal sense, and adds one of the two as a start element of the
//! `ShellFaceSet` — the surface patch that the split faces of both shapes can
//! share.
//!
//! UNPORTED (module-local): `FillOnPatches` (`cxx:1782-2035`),
//! `AreFacesCoincideInArea` (`cxx:1595-1780`) and `FindFacesTouchingEdge`
//! (`cxx:2037-2085`) belong to the `GSplitFace` ON path and are not ported
//! here.

use crate::abs::{Orientation, ShapeType, State};
use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{explorer_shapes, shape_map_key, ShapeMapKey};
use crate::topopebrep::build::g_topo::GTopo;
use crate::topopebrep::build::shell_face_set::ShellFaceSet;

use super::Builder;

impl Builder {
    /// `AddONPatchesSFS(G1, SFS)` (`cxx:1454-1592`).
    pub fn add_on_patches_sfs(&mut self, g1: &GTopo, sfs: &mut ShellFaceSet) {
        // `cxx:1458`: `cos(PI/180)`.
        let scal_min = 0.999847695;

        // `cxx:1489`: `PAR_T = 0.456321`.
        const PAR_T: f64 = 0.456321;

        // (patch face, ancestor 1, patch face 2, ancestor 2).
        let mut patches: Vec<(TopoShape, TopoShape, TopoShape, TopoShape)> = Vec::new();

        // `myONFacesMap`: key = new ON face, value = its ancestor face.
        let n = self.my_on_faces_map.extent();
        for i in 1..=n {
            let Some(a_f_anc1) = self.my_on_faces_map.find_from_index(i).cloned() else {
                continue;
            };
            if self.rank_of(&a_f_anc1) != 1 {
                continue;
            }
            let Some(a_face1) = self.my_on_faces_map.find_key(i).cloned() else {
                continue;
            };

            // `cxx:1469-1470`: the edges of the first face.
            let a_map_e1 = edge_map(&a_face1);

            // `cxx:1472-1483`: the first non-degenerated edge.
            let Some(a_chk_edge) = a_map_e1
                .iter()
                .map(|(_, e)| e)
                .find(|e| !BRepTool::is_degenerated(&Edge((*e).clone())))
                .cloned()
            else {
                continue;
            };
            let (Some(a_face1_val), Some(a_chk_edge_val)) =
                (Face::wrap(a_face1.clone()), Edge::wrap(a_chk_edge.clone()))
            else {
                continue;
            };

            // `cxx:1487-1508`: a point and the normal on the first face.
            let a_bac1 = Curve2dAdaptor::of_edge_face(&a_chk_edge_val, &a_face1_val);
            let par = a_bac1.first_parameter() * (1. - PAR_T) + a_bac1.last_parameter() * PAR_T;
            let Some(a_p2d) = a_bac1.d0(par) else {
                continue;
            };
            let Some(a_n1) = face_normal(&a_face1_val, a_p2d.x(), a_p2d.y(), a_face1.orientation())
            else {
                continue;
            };

            // `cxx:1510-1589`: the rank-2 faces of the map.
            let mut ok = true;
            for j in (i + 1)..=n {
                if !ok {
                    break;
                }
                let Some(a_f_anc2) = self.my_on_faces_map.find_from_index(j).cloned() else {
                    continue;
                };
                if self.rank_of(&a_f_anc2) != 2 {
                    continue;
                }
                let Some(a_face2) = self.my_on_faces_map.find_key(j).cloned() else {
                    continue;
                };

                // `cxx:1515-1533`: the same boundaries.
                let a_map_e2 = edge_map(&a_face2);
                if a_map_e1.len() != a_map_e2.len() {
                    continue;
                }
                let mut same_bnd = true;
                for (key, _) in &a_map_e2 {
                    if !a_map_e1.iter().any(|(k, _)| k == key) {
                        same_bnd = false;
                        break;
                    }
                }
                if !same_bnd {
                    continue;
                }

                // `cxx:1536-1561`: the normals must have the same sense.
                let Some(a_face2_val) = Face::wrap(a_face2.clone()) else {
                    continue;
                };
                let a_bac2 = Curve2dAdaptor::of_edge_face(&a_chk_edge_val, &a_face2_val);
                let Some(a_p2d2) = a_bac2.d0(par) else {
                    ok = false;
                    continue;
                };
                let Some(a_n2) =
                    face_normal(&a_face2_val, a_p2d2.x(), a_p2d2.y(), a_face2.orientation())
                else {
                    ok = false;
                    continue;
                };
                let scal = a_n1.dot(&a_n2);
                if scal < scal_min {
                    ok = false;
                    continue;
                }

                // `cxx:1564-1578`: `takeFirst` is hard-coded to true in OCCT.
                patches.push((a_face1.clone(), a_f_anc1.clone(), a_face2, a_f_anc2.clone()));
            }
        }

        // `cxx:1579-1589`: add the patch and save the ON splits.
        for (a_face1, a_f_anc1, a_face2, a_f_anc2) in patches {
            let mut a_patch = a_face1.clone();
            let newori_f = Builder::orient(a_f_anc1.orientation(), g1.is_to_reverse1());
            a_patch.set_orientation(newori_f);
            sfs.add_start_element(&a_patch);

            self.mark_split(&a_f_anc1, State::On, true);
            self.change_split(&a_f_anc1, State::On).push(a_face1);

            self.mark_split(&a_f_anc2, State::On, true);
            self.change_split(&a_f_anc2, State::On).push(a_face2);
        }
    }

    /// `BDS.AncestorRank(S)`.
    fn rank_of(&self, s: &TopoShape) -> i32 {
        match self.data_structure() {
            Some(hds) => hds.ds().ancestor_rank_of(s),
            None => 0,
        }
    }
}

/// `TopExp::MapShapes(F, EDGE, M)` as an ordered key/shape list.
fn edge_map(face: &TopoShape) -> Vec<(ShapeMapKey, TopoShape)> {
    let mut m: Vec<(ShapeMapKey, TopoShape)> = Vec::new();
    for e in explorer_shapes(face, ShapeType::Edge) {
        let key = shape_map_key(&e);
        if !m.iter().any(|(k, _)| *k == key) {
            m.push((key, e));
        }
    }
    m
}

/// `cxx:1492-1507` / `cxx:1541-1556`: the unit normal of `face` at `(u, v)`,
/// reversed for a `REVERSED` face. `None` when the gradient vanishes or is
/// below `Precision::Confusion()`.
fn face_normal(face: &Face, u: f64, v: f64, ori: Orientation) -> Option<occt_core::gp::GpVec> {
    let surf = BRepTool::face_surface(face)?;
    let (_, du, dv) = surf.d1(u, v);
    let mut n = du.crossed(&dv);
    let norm = n.magnitude();
    if norm < occt_core::precision::CONFUSION {
        return None;
    }
    n = n.divided(norm);
    if ori == Orientation::Reversed {
        n.reverse();
    }
    Some(n)
}
