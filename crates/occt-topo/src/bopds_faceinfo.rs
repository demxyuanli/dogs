//! BopdsDS FaceInfo On/In update and refine methods.
use std::collections::HashSet;

use crate::abs::ShapeType;
use crate::bopds::BopdsDS;
use crate::bopds_cb::{BopdsCommonBlock, BopdsFaceInfo};
use crate::bopds_pave::BopdsPaveBlock;
use crate::bopds_tree::direct_children;

impl BopdsDS {
    /// Returns the face-info pool.
    pub fn face_info_pool(&self) -> &[BopdsFaceInfo] {
        &self.face_info_pool
    }

    /// Mutable access to the face-info pool.
    pub fn change_face_info_pool(&mut self) -> &mut Vec<BopdsFaceInfo> {
        &mut self.face_info_pool
    }

    /// Face-info entry of face `i`, if any.
    pub fn face_info(&self, i: usize) -> Option<&BopdsFaceInfo> {
        self.face_info_pool.iter().find(|fi| fi.face_index == i)
    }

    /// Mutable face-info entry of face `i`, if any.
    pub fn face_info_mut(&mut self, i: usize) -> Option<&mut BopdsFaceInfo> {
        self.face_info_pool.iter_mut().find(|fi| fi.face_index == i)
    }

    /// On/In vertices and pave-block tuples of two faces, plus the common subset.
    ///
    /// Source: `BOPDS_DS::SubShapesOnIn` (`BOPDS_DS.cxx`). FaceInfo stores
    /// `(edge, first, last)` rather than pave-block handles; bound vertices of
    /// those blocks (and `verts_in` / `verts`) fill the vertex maps.
    pub fn sub_shapes_on_in(
        &self,
        n_f1: usize,
        n_f2: usize,
    ) -> (HashSet<usize>, HashSet<usize>, Vec<(usize, f64, f64)>) {
        let mut verts_on_in: HashSet<usize> = HashSet::new();
        let mut verts_common: HashSet<usize> = HashSet::new();
        let mut pbs_on_in: Vec<(usize, f64, f64)> = Vec::new();
        let Some(fi1) = self.face_info(n_f1) else {
            return (verts_on_in, verts_common, pbs_on_in);
        };
        let Some(fi2) = self.face_info(n_f2) else {
            return (verts_on_in, verts_common, pbs_on_in);
        };
        let push_pbs = |pbs_on_in: &mut Vec<(usize, f64, f64)>, verts: &mut HashSet<usize>, pbs: &[(usize, f64, f64)]| {
            for &(e, f, l) in pbs {
                if !pbs_on_in.iter().any(|&(ee, a, b)| {
                    ee == e && (a - f).abs() <= 1e-7 && (b - l).abs() <= 1e-7
                }) {
                    pbs_on_in.push((e, f, l));
                }
                for pb in self.pave_blocks(e) {
                    if (pb.first - f).abs() <= 1e-7 && (pb.last - l).abs() <= 1e-7 {
                        verts.insert(pb.index1);
                        verts.insert(pb.index2);
                    }
                }
            }
        };
        push_pbs(&mut pbs_on_in, &mut verts_on_in, fi1.paves_on());
        push_pbs(&mut pbs_on_in, &mut verts_on_in, fi1.paves_in());
        push_pbs(&mut pbs_on_in, &mut verts_on_in, fi2.paves_on());
        push_pbs(&mut pbs_on_in, &mut verts_on_in, fi2.paves_in());

        // OCCT `findCommon`: a On/In pave block of face1 that is also On/In
        // of face2 contributes its bound vertices to `theMVCommon`.
        let pb_on_in_face2 = |e: usize, f: f64, l: f64| {
            fi2.paves_on()
                .iter()
                .chain(fi2.paves_in())
                .any(|&(ee, a, b)| ee == e && (a - f).abs() <= 1e-7 && (b - l).abs() <= 1e-7)
        };
        for &(e, first, last) in fi1.paves_on().iter().chain(fi1.paves_in()) {
            if !pb_on_in_face2(e, first, last) {
                continue;
            }
            for pb in self.pave_blocks(e) {
                if (pb.first - first).abs() <= 1e-7 && (pb.last - last).abs() <= 1e-7 {
                    verts_common.insert(pb.index1);
                    verts_common.insert(pb.index2);
                }
            }
        }

        // OCCT `BOPDS_DS.cxx:1124-1142`: VerticesOn/In of face1 that are also
        // On or In on face2 go into both `theMVOnIn` and `theMVCommon`. The UV
        // `verts()` list is not this map. Face2 vertices that are not in that
        // intersection are not added here; their pave-block endpoints already
        // entered `verts_on_in` via `processMap`.
        let on_in2: HashSet<usize> = fi2
            .verts_on()
            .iter()
            .chain(fi2.verts_in())
            .copied()
            .collect();
        for &v in fi1.verts_on().iter().chain(fi1.verts_in()) {
            if on_in2.contains(&v) {
                verts_on_in.insert(v);
                verts_common.insert(v);
            }
        }
        (verts_on_in, verts_common, pbs_on_in)
    }

    /// Edge indices that belong to both faces (via pave-block real edges).
    ///
    /// Source: `BOPDS_DS::SharedEdges`.
    pub fn shared_edges(&self, n_f1: usize, n_f2: usize) -> Vec<usize> {
        let mut first: HashSet<usize> = HashSet::new();
        if let Some(si) = self.shape_info(n_f1) {
            for &sub in si.sub_shapes() {
                if self.shape_info(sub).map(|s| s.shape_type()) != Some(ShapeType::Edge) {
                    continue;
                }
                let pbs = self.pave_blocks(sub);
                if pbs.is_empty() {
                    first.insert(sub);
                } else {
                    for pb in pbs {
                        first.insert(self.real_pave_block(pb).edge());
                    }
                }
            }
        }
        let mut out = Vec::new();
        if let Some(si) = self.shape_info(n_f2) {
            for &sub in si.sub_shapes() {
                if self.shape_info(sub).map(|s| s.shape_type()) != Some(ShapeType::Edge) {
                    continue;
                }
                let pbs = self.pave_blocks(sub);
                if pbs.is_empty() {
                    if first.contains(&sub) {
                        out.push(sub);
                    }
                } else {
                    for pb in pbs {
                        let e = self.real_pave_block(pb).edge();
                        if first.contains(&e) {
                            out.push(e);
                        }
                    }
                }
            }
        }
        out
    }

    /// Builds the ON set of the face `i`: the pave blocks of the face's own
    /// boundary edges and the ON vertex indices.
    ///
    /// Source: `BOPDS_DS::UpdateFaceInfoOn` + `FaceInfoOn`. A face without a
    /// face-info entry is skipped (`!HasReference()`). For each boundary EDGE
    /// the edge's pave-block endpoints feed `VerticesOn` and `RealPaveBlock`
    /// feeds `PaveBlocksOn`; boundary VERTEX sub-shapes feed `VerticesOn`.
    pub fn update_face_info_on(&mut self, i: usize) {
        if !self.face_info_pool.iter().any(|fi| fi.face_index == i) {
            return;
        }
        let subs: Vec<usize> = self
            .shape_info(i)
            .map(|s| s.sub_shapes().to_vec())
            .unwrap_or_default();
        let mut on = Vec::new();
        let mut verts_on: Vec<usize> = Vec::new();
        for &sub in &subs {
            let ty = self.shape_info(sub).map(|s| s.shape_type());
            if ty == Some(ShapeType::Edge) {
                let pbs: Vec<BopdsPaveBlock> = self.pave_blocks(sub).to_vec();
                for pb in &pbs {
                    let (n1, n2) = pb.indices();
                    verts_on.push(n1);
                    verts_on.push(n2);
                    let rpb = self.real_pave_block(pb);
                    on.push((rpb.edge(), rpb.first, rpb.last));
                }
            } else if ty == Some(ShapeType::Vertex) {
                verts_on.push(self.get_same_domain_index(sub));
            }
        }
        verts_on.sort_unstable();
        verts_on.dedup();
        if let Some(fi) = self.face_info_pool.iter_mut().find(|fi| fi.face_index == i) {
            fi.paves_on = on;
            fi.verts_on = verts_on;
        }
    }

    /// Ensures a face-info pool entry for `i` and initializes On/In like
    /// `BOPDS_DS::ChangeFaceInfo` → `InitFaceInfo`.
    pub(crate) fn ensure_face_info(&mut self, i: usize) {
        if self.face_info_pool.iter().any(|fi| fi.face_index == i) {
            return;
        }
        self.face_info_pool.push(BopdsFaceInfo::new(i));
        self.init_face_info_in(i);
        self.update_face_info_on(i);
    }

    /// Adds the face's direct VERTEX children to `VerticesIn`.
    /// Source: `BOPDS_DS::InitFaceInfoIn` — `TopoDS_Iterator` on the face.
    fn init_face_info_in(&mut self, n_f: usize) {
        let Some(shape) = self.shape(n_f).cloned() else { return };
        let mut verts = Vec::new();
        for child in direct_children(&shape) {
            if child.shape_type() != ShapeType::Vertex {
                continue;
            }
            if let Some(i) = self.index(&child) {
                verts.push(self.get_same_domain_index(i));
            }
        }
        if let Some(fi) = self.face_info_pool.iter_mut().find(|fi| fi.face_index == n_f) {
            for v in verts {
                fi.add_vert_in(v);
            }
        }
    }

    /// The first pave block of a common block (`BOPDS_CommonBlock::PaveBlock1`).
    fn pave_block1_tuple(&self, cb: &BopdsCommonBlock) -> Option<(usize, f64, f64)> {
        let pb = cb.pave_block1()?;
        let e = if pb.has_edge() {
            pb.edge()
        } else {
            pb.original_edge()
        };
        Some((e, pb.first, pb.last))
    }

    /// IN pave blocks recovered from the common blocks of edge `n_e` on face `n_f`.
    ///
    /// Source: `BOPDS_DS::UpdateFaceInfoIn` E/F branch (`BOPDS_DS.cxx:940-946`):
    /// each pave block of `Index1` whose `CommonBlock(pb)->Contains(face)`
    /// contributes `PaveBlock1()`.
    fn in_paves_from_ef_common_block(&self, n_e: usize, n_f: usize) -> Vec<(usize, f64, f64)> {
        let mut out: Vec<(usize, f64, f64)> = Vec::new();
        let push_unique = |out: &mut Vec<(usize, f64, f64)>, t: (usize, f64, f64)| {
            if !out
                .iter()
                .any(|&(e, f, l)| e == t.0 && (f - t.1).abs() <= 1e-7 && (l - t.2).abs() <= 1e-7)
            {
                out.push(t);
            }
        };
        for pb in self.pave_blocks(n_e) {
            let Some(cb) = self.common_block(pb) else {
                continue;
            };
            if !cb.contains_face(n_f) {
                continue;
            }
            let t = self
                .pave_block1_tuple(cb)
                .unwrap_or((pb.edge(), pb.first, pb.last));
            push_unique(&mut out, t);
        }
        out
    }

    /// Builds the ON set for every face in `faces`, creating face-info entries
    /// as needed. Source: `BOPDS_DS::UpdateFaceInfoOn(const NCollection_Map<int>&)`.
    pub fn update_face_info_on_faces(&mut self, faces: &HashSet<usize>) {
        for &i in faces {
            self.ensure_face_info(i);
            self.update_face_info_on(i);
        }
    }

    /// Rebuilds the IN set of face `i` from internal vertices, V/F, and E/F
    /// common blocks. Skips a face without a face-info entry.
    /// Source: `BOPDS_DS::UpdateFaceInfoIn(const int)`.
    pub fn update_face_info_in(&mut self, i: usize) {
        if !self.face_info_pool.iter().any(|fi| fi.face_index == i) {
            return;
        }
        self.update_face_info_in_faces(&HashSet::from([i]));
    }

    /// Rebuilds the IN set of every face in `faces`.
    ///
    /// Source: `BOPDS_DS::UpdateFaceInfoIn(const NCollection_Map<int>&)`
    /// (`BOPDS_DS.cxx`): ensure FaceInfo, clear `PaveBlocksIn`/`VerticesIn`,
    /// `InitFaceInfoIn`, then scan `InterfVF` and `InterfEF`.
    pub fn update_face_info_in_faces(&mut self, faces: &HashSet<usize>) {
        if faces.is_empty() {
            return;
        }
        for &n_f in faces {
            self.ensure_face_info(n_f);
            if let Some(fi) = self.face_info_pool.iter_mut().find(|fi| fi.face_index == n_f) {
                fi.paves_in.clear();
                fi.verts_in.clear();
            }
            self.init_face_info_in(n_f);
        }

        let vf: Vec<(usize, usize)> = self
            .interf_vf
            .iter()
            .filter(|it| faces.contains(&it.index2))
            .map(|it| (it.index2, self.get_same_domain_index(it.index1)))
            .collect();
        for (n_f, n_v) in vf {
            if let Some(fi) = self.face_info_pool.iter_mut().find(|fi| fi.face_index == n_f) {
                fi.add_vert_in(n_v);
            }
        }

        let ef = self.interf_ef.clone();
        for it in &ef {
            let n_f = it.index2;
            if !faces.contains(&n_f) {
                continue;
            }
            if let Some(n_new) = it.get_index_new() {
                let n_v = self.get_same_domain_index(n_new);
                if let Some(fi) = self.face_info_pool.iter_mut().find(|fi| fi.face_index == n_f)
                {
                    fi.add_vert_in(n_v);
                }
                continue;
            }
            let in_pbs = self.in_paves_from_ef_common_block(it.index1, n_f);
            if let Some(fi) = self.face_info_pool.iter_mut().find(|fi| fi.face_index == n_f) {
                for (e, t1, t2) in in_pbs {
                    if !fi.paves_in.iter().any(|&(x, f, l)| {
                        x == e && (f - t1).abs() <= 1e-7 && (l - t2).abs() <= 1e-7
                    }) {
                        fi.add_pave_in(e, t1, t2);
                    }
                }
            }
        }
    }

    /// Removes from the IN pave blocks of every source face the blocks that
    /// are also ON (boundary) pave blocks.
    /// Source: `BOPDS_DS::RefineFaceInfoIn`.
    ///
    /// OCCT iterates the source shapes and skips the faces without a face-info
    /// entry; the pool holds exactly those entries, so iterating it directly
    /// is equivalent. The blocks are compared by `(edge, first, last)`, the
    /// handle-identity test of the OCCT `IndexedMap`s projected onto our
    /// tuple representation of the face-info sets.
    pub fn refine_face_info_in(&mut self) {
        for fi in self.face_info_pool.iter_mut() {
            if fi.paves_in.is_empty() || fi.paves_on.is_empty() {
                continue;
            }
            let on = &fi.paves_on;
            fi.paves_in.retain(|pb| !on.contains(pb));
        }
    }

    /// Rebuilds the ON sets of all faces with a face-info entry and drops the
    /// ON blocks without an assigned edge.
    ///
    /// Source: `BOPDS_DS::RefineFaceInfoOn` (`BOPDS_DS.cxx`): for every
    /// face-info pool entry the ON set is rebuilt from the face's boundary
    /// edges (`UpdateFaceInfoOn` + `FaceInfoOn`, here
    /// [`BopdsDS::update_face_info_on`]), then the ON blocks whose edge is not
    /// set (`!HasEdge()`) are removed. OCCT's `HasEdge` filter is kept for
    /// faithfulness (an ON entry with `edge == usize::MAX`), even though the
    /// rebuild only ever records blocks with an assigned edge.
    pub fn refine_face_info_on(&mut self) {
        let indices: Vec<usize> = self.face_info_pool.iter().map(|fi| fi.face_index).collect();
        for i in indices {
            self.update_face_info_on(i);
        }
        for fi in &mut self.face_info_pool {
            fi.paves_on.retain(|(edge, _, _)| *edge != usize::MAX);
        }
    }
}
