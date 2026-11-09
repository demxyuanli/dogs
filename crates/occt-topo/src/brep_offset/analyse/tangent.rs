//! `BRepOffset_Analyse::TreatTangentFaces` (`BRepOffset_Analyse.cxx:371-806`).
//!
//! For every edge whose two faces are tangent but have different offset
//! values, an artificial perpendicular face is built by extruding the edge
//! along the normal of its first face, and that face is inserted between the
//! two faces. The extruded geometry is `BRepPrimAPI_MakePrism` over a compound
//! of edges, i.e. [`Prism`] (`BRepSweep_Prism`), whose `Generated` history is
//! `BRepPrimAPI_MakePrism::Generated` (`BRepPrimAPI_MakePrism.cxx:88-96`).

use std::collections::HashMap;

use occt_core::gp::GpVec;

use crate::abs::{Orientation, ShapeType};
use crate::algo_tools3d::{normal_to_face_on_edge, point_near_edge_dt};
use crate::algo_tools_face::is_split_to_reverse_edge;
use crate::brep_sweep::Prism;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::chfi_ds::TypeOfConcavity;
use crate::int_tools_full::IntToolsContext;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{
    children_keep_orientation, explorer_shapes, is_equal, is_same, map_shapes_and_ancestors_into,
    shape_map_key, ShapeMapKey,
};

use super::Analyse;

/// `BRepPrimAPI_MakePrism::Generated(S)` (`BRepPrimAPI_MakePrism.cxx:88-96`):
/// `myPrism.Shape(S)` when `IsUsed(S) && GenIsUsed(S)`, otherwise the static
/// empty list.
fn prism_generated(prism: &Prism, gen: &TopoShape) -> Option<TopoShape> {
    if prism.is_used(gen) && prism.gen_is_used(gen) {
        prism.shape_of(gen)
    } else {
        None
    }
}

/// `BOPTools_AlgoTools::MakeConnexityBlocks(theS, VERTEX, EDGE, theLCB,
/// theConnectionMap)` (`BOPTools_AlgoTools.cxx:105-154`).
///
/// `theConnectionMap` is filled by `TopExp::MapShapesAndAncestors(theS,
/// VERTEX, EDGE)` (non-unique, so an edge met twice is listed twice) and the
/// blocks are walked breadth-first in the element order of `theS`, in the
/// order the connection map lists the edges of a vertex.
fn make_connexity_blocks_vertex_edge(
    s: &TopoShape,
) -> (Vec<Vec<TopoShape>>, HashMap<ShapeMapKey, Vec<TopoShape>>) {
    let mut connection: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
    for anc in explorer_shapes(s, ShapeType::Edge) {
        for sub in explorer_shapes(&anc, ShapeType::Vertex) {
            connection
                .entry(shape_map_key(&sub))
                .or_default()
                .push(anc.clone());
        }
    }

    let mut fence: Vec<ShapeMapKey> = Vec::new();
    let mut blocks: Vec<Vec<TopoShape>> = Vec::new();
    for cur in explorer_shapes(s, ShapeType::Edge) {
        if fence.contains(&shape_map_key(&cur)) {
            continue;
        }
        fence.push(shape_map_key(&cur));
        let mut block = vec![cur];
        let mut i = 0;
        while i < block.len() {
            let s1 = block[i].clone();
            i += 1;
            for sub in explorer_shapes(&s1, ShapeType::Vertex) {
                let Some(edges) = connection.get(&shape_map_key(&sub)) else {
                    continue;
                };
                for e2 in edges.clone() {
                    if !fence.contains(&shape_map_key(&e2)) {
                        fence.push(shape_map_key(&e2));
                        block.push(e2);
                    }
                }
            }
        }
        blocks.push(block);
    }
    (blocks, connection)
}

/// `gp_Vec(aDir.XYZ())` - the coordinate triple of a direction.
fn dir_vec(d: &occt_core::gp::GpDir) -> GpVec {
    GpVec::new(d.x(), d.y(), d.z())
}

impl Analyse {
    /// `myFaceOffsetMap.Seek(face)` with the `myOffset` fallback of
    /// `cxx:404-408`.
    fn face_offset(&self, f: &TopoShape) -> f64 {
        match self.face_offset_map.get(&shape_map_key(f)) {
            Some(v) => v.abs(),
            None => self.offset,
        }
    }

    /// `TreatTangentFaces(theLE, theRange)` (`cxx:371-806`).
    ///
    /// `Err` carries the exceptions OCCT lets escape (`Standard_ConstructionError`
    /// of `BRepPrimAPI_MakePrism` on a zero extrusion vector, and the null
    /// handles the reference dereferences without a check).
    pub(super) fn treat_tangent_faces(&mut self, le: &[TopoShape]) -> Result<(), String> {
        if le.is_empty() || self.face_offset_map.is_empty() {
            // Nothing to do: either there are no tangent faces in the shape or
            // the face offset map has not been provided (`cxx:374-380`).
            return Ok(());
        }

        // Select the edges which connect faces with different offset values
        // (`cxx:382-386`).
        let mut ce_tangent: Vec<TopoShape> = Vec::new();
        let mut edge_offset_map: HashMap<ShapeMapKey, f64> = HashMap::new();
        let mut dmvemin: HashMap<ShapeMapKey, TopoShape> = HashMap::new();

        for a_e in le {
            let a_la = self.ancestors(a_e).to_vec();
            let (Some(a_f1), Some(a_f2)) = (a_la.first().cloned(), a_la.last().cloned()) else {
                continue;
            };

            let val1 = self.face_offset(&a_f1);
            let val2 = self.face_offset(&a_f2);
            if val1 == val2 {
                continue;
            }

            ce_tangent.push(a_e.clone());
            edge_offset_map.insert(shape_map_key(a_e), val1.max(val2));

            let a_f_min = if val1 < val2 { a_f1 } else { a_f2 };
            for v in children_keep_orientation(a_e) {
                if v.shape_type() != ShapeType::Vertex {
                    continue;
                }
                if self.ancestors(&v).len() != 3 {
                    continue;
                }
                for cur in explorer_shapes(&a_f_min, ShapeType::Edge) {
                    if is_same(&cur, a_e) {
                        continue;
                    }
                    for vx in children_keep_orientation(&cur) {
                        if vx.shape_type() != ShapeType::Vertex {
                            continue;
                        }
                        if is_same(&v, &vx) {
                            dmvemin.insert(shape_map_key(&v), cur.clone());
                        }
                    }
                }
            }
        }

        if edge_offset_map.is_empty() {
            return Ok(());
        }

        // Create map of Face ancestors for the vertices on tangent edges
        // (`cxx:444-502`).
        let mut dmv_fanc: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
        for a_e in le {
            if !edge_offset_map.contains_key(&shape_map_key(a_e)) {
                continue;
            }

            let mut a_mfence: Vec<ShapeMapKey> = Vec::new();
            for anc in self.ancestors(a_e) {
                a_mfence.push(shape_map_key(anc));
            }

            for v in children_keep_orientation(a_e) {
                if v.shape_type() != ShapeType::Vertex {
                    continue;
                }
                let vk = shape_map_key(&v);
                dmv_fanc.entry(vk.clone()).or_default();

                let lva = self.ancestors(&v).to_vec();
                for a_ea in lva {
                    let e_a = Edge(a_ea);
                    let Some(intervals) = self.edge_types(&e_a).cloned() else {
                        continue;
                    };
                    if intervals.is_empty() {
                        continue;
                    }
                    if intervals[0].kind() == TypeOfConcavity::Tangential {
                        continue;
                    }

                    let lea = self.ancestors(&e_a.0).to_vec();
                    for a_fa in lea {
                        let kfa = shape_map_key(&a_fa);
                        if !a_mfence.contains(&kfa) {
                            a_mfence.push(kfa);
                            if let Some(list) = dmv_fanc.get_mut(&vk) {
                                list.push(a_fa);
                            }
                        }
                    }
                }
            }
        }

        let ctx = IntToolsContext::new();
        // Tangency criteria (`cxx:504-505`).
        let a_sin_tol = self.angle.sin().abs();

        // Make blocks of connected edges (`cxx:507-511`).
        let ce_compound = TopoBuilder::new().make_compound_of(&ce_tangent);
        let (blocks, _mvemap) = make_connexity_blocks_vertex_edge(&ce_compound.0);

        // Analyze each block to find co-planar edges (`cxx:513-803`).
        for cb in &blocks {
            let mut a_mfence: Vec<ShapeMapKey> = Vec::new();
            let mut i = 0;
            while i < cb.len() {
                let a_e1 = Edge(cb[i].clone());
                let it_cb1 = i;
                i += 1;

                if a_mfence.contains(&shape_map_key(&a_e1.0)) {
                    continue;
                }
                a_mfence.push(shape_map_key(&a_e1.0));

                let mut block_edges: Vec<TopoShape> =
                    vec![a_e1.0.oriented(Orientation::Forward)];

                let mut an_offset = *edge_offset_map
                    .get(&shape_map_key(&a_e1.0))
                    .ok_or("BRepOffset_Analyse: tangent edge without offset")?;

                let a_lf1 = self.ancestors(&a_e1.0).to_vec();
                let a_face1 = a_lf1
                    .first()
                    .cloned()
                    .ok_or("BRepOffset_Analyse: tangent edge without face")?;
                let a_dn1 = normal_to_face_on_edge(&a_e1, &Face(a_face1), &ctx).ok_or(
                    "BOPTools_AlgoTools3D::GetNormalToFaceOnEdge: no normal",
                )?;

                for e2 in cb.iter().skip(it_cb1 + 1) {
                    let a_e2 = Edge(e2.clone());
                    if a_mfence.contains(&shape_map_key(&a_e2.0)) {
                        continue;
                    }

                    let a_lf2 = self.ancestors(&a_e2.0).to_vec();
                    let a_face2 = a_lf2
                        .first()
                        .cloned()
                        .ok_or("BRepOffset_Analyse: tangent edge without face")?;
                    let a_dn2 = normal_to_face_on_edge(&a_e2, &Face(a_face2), &ctx).ok_or(
                        "BOPTools_AlgoTools3D::GetNormalToFaceOnEdge: no normal",
                    )?;

                    if dir_vec(&a_dn1).crossed(&dir_vec(&a_dn2)).magnitude() < a_sin_tol {
                        block_edges.push(a_e2.0.oriented(Orientation::Forward));
                        a_mfence.push(shape_map_key(&a_e2.0));
                        let other = *edge_offset_map
                            .get(&shape_map_key(&a_e2.0))
                            .ok_or("BRepOffset_Analyse: tangent edge without offset")?;
                        an_offset = an_offset.max(other);
                    }
                }

                // Make the prism (`cxx:567-572`).
                let block = TopoBuilder::new().make_compound_of(&block_edges);
                let vec = dir_vec(&a_dn1).multiplied_scalar(an_offset);
                let prism = Prism::from_vec(block.0.clone(), &vec, false, true)
                    .map_err(|e| format!("BRepPrimAPI_MakePrism: {e}"))?;
                let Some(prism_shape) = prism.shape() else {
                    // `if (!aMP.IsDone()) continue;` (`cxx:569-572`).
                    continue;
                };

                let mut prism_ancestors: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
                map_shapes_and_ancestors_into(
                    &prism_shape,
                    ShapeType::Edge,
                    ShapeType::Face,
                    &mut prism_ancestors,
                );
                map_shapes_and_ancestors_into(
                    &prism_shape,
                    ShapeType::Vertex,
                    ShapeType::Edge,
                    &mut prism_ancestors,
                );

                for child in children_keep_orientation(&block.0) {
                    if child.shape_type() != ShapeType::Edge {
                        continue;
                    }
                    let e = Edge(child);
                    self.process_prism_edge(
                        &e,
                        &prism,
                        &prism_ancestors,
                        &dmvemin,
                        &dmv_fanc,
                        a_sin_tol,
                        &ctx,
                    )?;
                }
            }
        }

        Ok(())
    }

    /// The body of the `for (TopoDS_Iterator itE(aBlock))` loop
    /// (`cxx:581-804`).
    #[allow(clippy::too_many_arguments)]
    fn process_prism_edge(
        &mut self,
        a_e: &Edge,
        prism: &Prism,
        prism_ancestors: &HashMap<ShapeMapKey, Vec<TopoShape>>,
        dmvemin: &HashMap<ShapeMapKey, TopoShape>,
        dmv_fanc: &HashMap<ShapeMapKey, Vec<TopoShape>>,
        a_sin_tol: f64,
        ctx: &IntToolsContext,
    ) -> Result<(), String> {
        // `const NCollection_List<TopoDS_Shape>& aLG = aMP.Generated(aE);
        //  TopoDS_Face aFNew = TopoDS::Face(aLG.First());` (`cxx:583-585`).
        // The reference dereferences an empty list when `Generated` is empty;
        // the port skips that edge.
        let Some(a_lg_first) = prism_generated(prism, &a_e.0) else {
            return Ok(());
        };
        let mut a_f_new = Face(a_lg_first);

        let a_la = self.ancestors(&a_e.0).to_vec();
        let (Some(a_f1), Some(a_f2)) = (a_la.first().cloned(), a_la.last().cloned()) else {
            return Ok(());
        };

        let val1 = self.face_offset(&a_f1);
        let val2 = self.face_offset(&a_f2);
        let (a_f_to_remove, a_f_opposite) = if val1 > val2 {
            (Face(a_f1), Face(a_f2))
        } else {
            (Face(a_f2), Face(a_f1))
        };

        // Orient the face so its normal is directed to smaller offset face
        // (`cxx:601-635`).
        {
            let a_dn = normal_to_face_on_edge(a_e, &a_f_new, ctx)
                .ok_or("BOPTools_AlgoTools3D::GetNormalToFaceOnEdge: no normal")?;

            let mut a_e_in_f: Option<Edge> = None;
            for cur in explorer_shapes(&a_f_opposite.0, ShapeType::Edge) {
                if is_same(&cur, &a_e.0) {
                    a_e_in_f = Some(Edge(cur));
                    break;
                }
            }
            let a_e_in_f =
                a_e_in_f.ok_or("BRepOffset_Analyse: edge not found in the opposite face")?;

            let (f, l) = BRepTool::edge_parameters(&a_e_in_f);
            let a_c3d = BRepTool::edge_curve(&a_e_in_f)
                .ok_or("BRepOffset_Analyse: null 3d curve on opposite edge")?;
            let a_p_on_e = a_c3d.d0((f + l) / 2.0);
            let Ok((_, a_p_in_f)) =
                point_near_edge_dt(&a_e_in_f, &a_f_opposite, (f + l) / 2.0, 1.0e-5)
            else {
                return Err("BOPTools_AlgoTools3D::PointNearEdge failed".into());
            };
            let a_bn = GpVec::from_pnts(&a_p_on_e, &a_p_in_f);

            if a_bn.dot(&dir_vec(&a_dn)) < 0.0 {
                a_f_new.0.reverse();
            }
        }

        // Remove the face with bigger offset value from edge ancestors
        // (`cxx:637-646`).
        {
            let la = self.ancestors_mut(&a_e.0);
            if let Some(pos) = la.iter().position(|x| is_same(x, &a_f_to_remove.0)) {
                la.remove(pos);
            }
            la.push(a_f_new.0.clone());
        }

        self.bind_edge_types(a_e, Vec::new());
        // Analyze edge again (`cxx:649-651`).
        self.edge_analyse(a_e, &a_f_opposite, &a_f_new, a_sin_tol)?;

        // Analyze vertices (`cxx:653-745`).
        let mut a_f_new_edge_map: Vec<ShapeMapKey> = vec![shape_map_key(&a_e.0)];
        for v in children_keep_orientation(&a_e.0) {
            if v.shape_type() != ShapeType::Vertex {
                continue;
            }
            // `TopoDS_Edge aEG = TopoDS::Edge(aMP.Generated(aV).First());`
            // (`cxx:659`). Empty `Generated` is a null dereference in the
            // reference; the port skips that vertex.
            let Some(gen) = prism_generated(prism, &v) else {
                continue;
            };
            let mut a_eg = Edge(gen);
            self.bind_generated(&v, a_eg.0.clone());
            for cur in explorer_shapes(&a_f_new.0, ShapeType::Edge) {
                if is_same(&cur, &a_eg.0) {
                    a_eg = Edge(cur);
                    break;
                }
            }

            if let Some(a_e_min) = dmvemin.get(&shape_map_key(&v)) {
                if let Some(p_sa) = dmv_fanc.get(&shape_map_key(&v)) {
                    if p_sa.len() == 1 {
                        // Adjust orientation of generated edge to its new
                        // ancestor (`cxx:676-712`).
                        let mut a_e_min = Edge(a_e_min.clone());
                        for cur in explorer_shapes(&p_sa[0], ShapeType::Edge) {
                            if is_same(&cur, &a_e_min.0) {
                                a_e_min = Edge(cur);
                                break;
                            }
                        }

                        let mut an_ori_in_e_min = Orientation::Forward;
                        for itx in children_keep_orientation(&a_e_min.0) {
                            if is_same(&itx, &v) {
                                an_ori_in_e_min = itx.orientation();
                                break;
                            }
                        }

                        let mut an_ori_in_eg = Orientation::Forward;
                        for itx in children_keep_orientation(&a_eg.0) {
                            if is_same(&itx, &v) {
                                an_ori_in_eg = itx.orientation();
                                break;
                            }
                        }

                        if an_ori_in_eg == an_ori_in_e_min {
                            a_eg.0.reverse();
                        }
                    }
                }
            }

            {
                let lva = self.ancestors_mut(&v);
                // `NCollection_List::Contains` compares with `IsSame`
                // (`NCollection_List.hxx:232-244`), so the orientation of the
                // possibly reversed `aEG` plays no part.
                if !lva.iter().any(|x| is_same(x, &a_eg.0)) {
                    lva.push(a_eg.0.clone());
                }
            }
            a_f_new_edge_map.push(shape_map_key(&a_eg.0));

            let prism_anc = prism_ancestors
                .get(&shape_map_key(&a_eg.0))
                .cloned()
                .unwrap_or_default();
            let idx = self.ancestors_add(&a_eg.0, &prism_anc);
            if let Some(extra) = dmv_fanc.get(&shape_map_key(&v)) {
                if !extra.is_empty() {
                    if let Some(list) = self.ancestors_at(idx) {
                        list.extend(extra.iter().cloned());
                    }
                }
            }

            self.bind_edge_types(&a_eg, Vec::new());
            let lega = self.ancestors(&a_eg.0).to_vec();
            if lega.len() == 2 {
                self.edge_analyse(
                    &a_eg,
                    &Face(lega[0].clone()),
                    &Face(lega[1].clone()),
                    a_sin_tol,
                )?;
            }
        }

        // Find an edge opposite to tangential one and add ancestors for it
        // (`cxx:747-799`). A null `aEOpposite` is dereferenced by the
        // reference; the port skips the block then.
        let mut a_e_opposite: Option<Edge> = None;
        for cur in explorer_shapes(&a_f_new.0, ShapeType::Edge) {
            if !a_f_new_edge_map.contains(&shape_map_key(&cur)) {
                a_e_opposite = Some(Edge(cur));
                break;
            }
        }

        if let Some(mut a_e_opposite) = a_e_opposite {
            // Find it in aFOpposite (`cxx:755-770`).
            for cur in explorer_shapes(&a_f_to_remove.0, ShapeType::Edge) {
                if is_same(&cur, &a_e.0) {
                    if is_split_to_reverse_edge(&a_e_opposite, &Edge(cur), ctx).unwrap_or(false) {
                        a_e_opposite.0.reverse();
                    }
                    break;
                }
            }

            let a_lf_opposite = vec![a_f_new.0.clone(), a_f_to_remove.0.clone()];
            self.ancestors_add(&a_e_opposite.0, &a_lf_opposite);
            self.bind_edge_types(&a_e_opposite, Vec::new());
            self.edge_analyse(&a_e_opposite, &a_f_new, &a_f_to_remove, a_sin_tol)?;

            self.replacement_map_mut(&a_f_to_remove)
                .insert(shape_map_key(&a_e.0), a_e_opposite.clone());

            // Add ancestors for the vertices (`cxx:792-799`).
            for v in children_keep_orientation(&a_e_opposite.0) {
                if v.shape_type() != ShapeType::Vertex {
                    continue;
                }
                let lva = prism_ancestors
                    .get(&shape_map_key(&v))
                    .cloned()
                    .unwrap_or_default();
                self.ancestors_add(&v, &lva);
            }
        }

        self.push_new_face(&a_f_new);
        self.bind_generated(&a_e.0, a_f_new.0.clone());

        Ok(())
    }
}
