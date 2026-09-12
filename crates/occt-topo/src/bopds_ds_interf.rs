//! BopdsDS common-block, interference and same-domain methods.
use std::collections::{HashMap, HashSet};

use crate::bopds::BopdsDS;
use crate::bopds_cb::{cb_contains_pb, pb_same, BopdsCommonBlock, BopdsInterf};
use crate::bopds_ff::BopdsInterfFf;
use crate::bopds_pave::BopdsPaveBlock;

impl BopdsDS {
    /// True if the pave block `pb` is a member of a common block: one of its
    /// edges (assigned or original) is an edge of a common block and its
    /// parameter range matches one of the common-block ranges.
    ///
    /// Port of `BOPDS_DS::IsCommonBlock` adapted to the tuple
    /// `(edge indices + ranges)` representation of the Rust common block.
    pub fn is_common_block(&self, pb: &BopdsPaveBlock) -> bool {
        self.common_block(pb).is_some()
    }

    /// `BOPDS_DS::CommonBlock`.
    pub fn common_block(&self, pb: &BopdsPaveBlock) -> Option<&BopdsCommonBlock> {
        self.common_blocks.iter().find(|cb| cb_contains_pb(cb, pb))
    }

    /// Append a common block and return its index.
    pub fn push_common_block(&mut self, cb: BopdsCommonBlock) -> usize {
        self.common_blocks.push(cb);
        self.common_blocks.len() - 1
    }

    /// Mutable common block at pool index `i`.
    pub fn common_block_mut(&mut self, i: usize) -> Option<&mut BopdsCommonBlock> {
        self.common_blocks.get_mut(i)
    }

    /// `BOPDS_DS::SetCommonBlock`: attach `pb` to `cb` (appended if not already
    /// stored as a member).
    pub fn set_common_block(&mut self, pb: &BopdsPaveBlock, mut cb: BopdsCommonBlock) {
        if !cb.pave_blocks.iter().any(|p| pb_same(p, pb)) {
            cb.add_pave_block(pb.clone());
        }
        if let Some(i) = self
            .common_blocks
            .iter()
            .position(|c| cb_contains_pb(c, pb))
        {
            self.common_blocks[i] = cb;
        } else {
            self.common_blocks.push(cb);
        }
    }

    /// Returns the real pave block of `pb`: the first block of the common
    /// block when `pb` is a common-block member, `pb` itself otherwise.
    ///
    /// Port of `BOPDS_DS::RealPaveBlock` (`BOPDS_DS.cxx`). The Rust common
    /// block stores the shared edge indices and ranges but not its own block
    /// list, so the first block is re-discovered as the block of the first
    /// common-block edge matching the first shared range (falls back to `pb`).
    pub fn real_pave_block(&self, pb: &BopdsPaveBlock) -> BopdsPaveBlock {
        if let Some(cb) = self.common_block(pb) {
            if let Some(first) = cb.pave_block1() {
                return first.clone();
            }
        }
        pb.clone()
    }

    /// Removes the reference to the pave blocks of the untouched edges, so no
    /// image is created for them later.
    ///
    /// Source: `BOPDS_DS::ReleasePaveBlocks` (`BOPDS_DS.cxx`). For every edge
    /// whose block list holds exactly one *untouched* block — the block is not
    /// a common-block member and both bound vertices are original (source)
    /// shapes — the reference to the block list is dropped from the edge's
    /// shape info and the list contents are cleared. The edge keeps a
    /// reference to an *empty* list, marking it as deleted: this distinguishes
    /// the small edges for which no pave block could even be built from the
    /// normal edges whose block was created but left untouched.
    pub fn release_pave_blocks(&mut self) {
        let nb_source = self.nb_source_shapes;
        let mut to_clear: Vec<usize> = Vec::new();
        let mut to_release: Vec<usize> = Vec::new();
        {
            for (slot, list) in self.pave_blocks_pool.iter().enumerate() {
                if list.len() != 1 {
                    continue;
                }
                let pb = &list[0];
                if self.is_common_block(pb) {
                    continue;
                }
                let (n1, n2) = pb.indices();
                if n1 < nb_source && n2 < nb_source {
                    to_clear.push(slot);
                    to_release.push(pb.original_edge());
                }
            }
        }
        for slot in to_clear {
            if let Some(list) = self.pave_blocks_pool.get_mut(slot) {
                list.clear();
            }
        }
        for orig in to_release {
            if orig < self.shape_infos.len() {
                self.shape_infos[orig].pb_reference = -1;
            }
        }
    }

    /// Removes every pave block whose assigned edge is in `edges` from the
    /// pave-block pool, F/F section curves, and face-info sets.
    /// Source: `BOPAlgo_PaveFiller::RemovePaveBlocks` (`BOPAlgo_PaveFiller_6.cxx`).
    pub fn remove_pave_blocks(&mut self, edges: &HashSet<usize>) {
        if edges.is_empty() {
            return;
        }
        // 1. Pave-blocks pool.
        for list in &mut self.pave_blocks_pool {
            list.retain(|pb| !edges.contains(&pb.edge()));
        }
        // 2. F/F section curves (`BOPDS_InterfFF::ChangeCurves`).
        for ff in &mut self.interf_ff {
            for c in ff.change_curves() {
                c.change_pave_blocks()
                    .retain(|pb| !edges.contains(&pb.edge()));
            }
        }
        // 3. Face-info sets (Sc / In / On).
        for fi in &mut self.face_info_pool {
            fi.paves.retain(|(e, _, _)| !edges.contains(e));
            fi.paves_in.retain(|(e, _, _)| !edges.contains(e));
            fi.paves_on.retain(|(e, _, _)| !edges.contains(e));
        }
    }

    /// Adds same-domain shape information (`index` and `index_sd`).
    pub fn add_shape_sd(&mut self, index: usize, index_sd: usize) {
        if index != index_sd {
            self.shapes_sd.insert(index, index_sd);
        }
    }

    /// Returns the same-domain index of `index` following the chain, or the
    /// index itself when no same-domain shape exists.
    pub fn get_same_domain_index(&self, index: usize) -> usize {
        let mut cur = index;
        while let Some(&next) = self.shapes_sd.get(&cur) {
            cur = next;
        }
        cur
    }

    /// Returns the final same-domain index of `index`, if any.
    pub fn has_shape_sd(&self, index: usize) -> Option<usize> {
        let result = self.get_same_domain_index(index);
        if result == index {
            None
        } else {
            Some(result)
        }
    }

    /// Returns the raw same-domain shape map (`index → same-domain index`).
    pub fn shapes_sd(&self) -> &HashMap<usize, usize> {
        &self.shapes_sd
    }

    /// Mutable access to the raw same-domain shape map.
    pub fn shapes_sd_mut(&mut self) -> &mut HashMap<usize, usize> {
        &mut self.shapes_sd
    }

    /// Returns the indices of the vertices whose tolerance was increased.
    pub fn increased_ss(&self) -> &HashSet<usize> {
        &self.increased_ss
    }

    /// Mutable access to the set of vertices with increased tolerance.
    pub fn increased_ss_mut(&mut self) -> &mut HashSet<usize> {
        &mut self.increased_ss
    }

    /// Adds an interference between shapes with indices `i1` and `i2`.
    /// Returns true when the pair was not already present.
    pub fn add_interf(&mut self, i1: usize, i2: usize) -> bool {
        let (a, b) = (i1.min(i2), i1.max(i2));
        if self.interferences.insert((a, b)) {
            self.interfered.insert(i1);
            self.interfered.insert(i2);
            true
        } else {
            false
        }
    }

    /// True if the shape with index `i` is interfered with any other shape.
    pub fn has_interf(&self, i: usize) -> bool {
        self.interfered.contains(&i)
    }

    /// True if the shapes with indices `i1` and `i2` are interfered.
    pub fn has_interf_pair(&self, i1: usize, i2: usize) -> bool {
        let (a, b) = (i1.min(i2), i1.max(i2));
        self.interferences.contains(&(a, b))
    }

    /// Returns the set of interfering pairs.
    pub fn interferences(&self) -> &HashSet<(usize, usize)> {
        &self.interferences
    }

    /// Returns the typed V/V interference records.
    pub fn interf_vv(&self) -> &[BopdsInterf] {
        &self.interf_vv
    }

    /// Returns the typed V/E interference records.
    pub fn interf_ve(&self) -> &[BopdsInterf] {
        &self.interf_ve
    }

    /// Returns the typed V/F interference records.
    pub fn interf_vf(&self) -> &[BopdsInterf] {
        &self.interf_vf
    }

    /// Returns the typed E/E interference records.
    pub fn interf_ee(&self) -> &[BopdsInterf] {
        &self.interf_ee
    }

    /// Returns the typed E/F interference records.
    pub fn interf_ef(&self) -> &[BopdsInterf] {
        &self.interf_ef
    }

    /// Returns the typed F/F interference records (`BOPDS_DS::InterfFF`).
    pub fn interf_ff(&self) -> &[BopdsInterfFf] {
        &self.interf_ff
    }

    /// Mutable F/F interference array.
    pub fn interf_ff_mut(&mut self) -> &mut Vec<BopdsInterfFf> {
        &mut self.interf_ff
    }

    /// Mutable V/V interference array (`BOPDS_DS::InterfVV`).
    pub fn interf_vv_mut(&mut self) -> &mut Vec<BopdsInterf> {
        &mut self.interf_vv
    }

    /// Mutable V/E interference array.
    pub fn interf_ve_mut(&mut self) -> &mut Vec<BopdsInterf> {
        &mut self.interf_ve
    }

    /// Mutable V/F interference array.
    pub fn interf_vf_mut(&mut self) -> &mut Vec<BopdsInterf> {
        &mut self.interf_vf
    }

    /// Mutable E/E interference array.
    pub fn interf_ee_mut(&mut self) -> &mut Vec<BopdsInterf> {
        &mut self.interf_ee
    }

    /// Mutable E/F interference array.
    pub fn interf_ef_mut(&mut self) -> &mut Vec<BopdsInterf> {
        &mut self.interf_ef
    }

    /// Appends a typed V/V interference record, registering the flat pair as
    /// well. The typed record is appended only when the pair was not interfered
    /// before (matching OCCT's `if (AddInterf(n1, n2)) { InterfVV().Appended() }`
    /// guard in `MakeSDVertices`). Returns true when a record was appended.
    pub fn add_interf_vv(&mut self, i1: usize, i2: usize, index_new: Option<usize>) -> bool {
        if !self.add_interf(i1, i2) {
            return false;
        }
        Self::append_interf(&mut self.interf_vv, i1, i2, index_new);
        true
    }

    /// Appends a typed V/E interference record (see [`BopdsDS::add_interf_vv`]).
    pub fn add_interf_ve(&mut self, i1: usize, i2: usize, index_new: Option<usize>) -> bool {
        if !self.add_interf(i1, i2) {
            return false;
        }
        Self::append_interf(&mut self.interf_ve, i1, i2, index_new);
        true
    }

    /// Appends a typed V/F interference record (see [`BopdsDS::add_interf_vv`]).
    pub fn add_interf_vf(&mut self, i1: usize, i2: usize, index_new: Option<usize>) -> bool {
        if !self.add_interf(i1, i2) {
            return false;
        }
        Self::append_interf(&mut self.interf_vf, i1, i2, index_new);
        true
    }

    /// Appends a typed E/E interference record (see [`BopdsDS::add_interf_vv`]).
    pub fn add_interf_ee(&mut self, i1: usize, i2: usize, index_new: Option<usize>) -> bool {
        if !self.add_interf(i1, i2) {
            return false;
        }
        Self::append_interf(&mut self.interf_ee, i1, i2, index_new);
        true
    }

    /// Appends a typed E/F interference record.
    ///
    /// Source: `BOPAlgo_PaveFiller_5.cxx` EDGE/VERTEX branches and
    /// `ForceInterfEF`: `aEFs.Appended()` then `AddInterf`. OCCT always
    /// appends a new `InterfEF` even when the flat pair is already present
    /// (a VERTEX record with `IndexNew` and an EDGE record without it must
    /// coexist so `UpdateFaceInfoIn` can still recover `PaveBlocksIn`).
    /// Returns the index of the appended record (`iX` / `IndexInterf`).
    pub fn add_interf_ef(&mut self, i1: usize, i2: usize, index_new: Option<usize>) -> usize {
        self.add_interf(i1, i2);
        Self::append_interf(&mut self.interf_ef, i1, i2, index_new);
        self.interf_ef.len() - 1
    }

    /// `BOPDS_Interf::SetIndexNew` on `InterfEF(iX)` (`PerformNewVertices`).
    pub fn bind_ef_new_vertex_at(&mut self, ix: usize, n_v: usize) -> bool {
        if let Some(it) = self.interf_ef.get_mut(ix) {
            it.set_index_new(n_v);
            true
        } else {
            false
        }
    }

    /// `IntTools_CommonPrt::Range1` on `InterfEF(iX)`.
    pub fn set_ef_common_range_at(&mut self, ix: usize, first: f64, last: f64) -> bool {
        if let Some(it) = self.interf_ef.get_mut(ix) {
            it.set_common_range(first, last);
            true
        } else {
            false
        }
    }

    /// Appends an F/F interference record (`BOPDS_DS::InterfFF`).
    /// `AddInterf` is called only when the record carries curves or points,
    /// matching `PerformFF` (`BOPAlgo_PaveFiller_6.cxx`).
    pub fn append_interf_ff(&mut self, rec: BopdsInterfFf) {
        if !rec.curves().is_empty() || !rec.points().is_empty() {
            self.add_interf(rec.index1(), rec.index2());
        }
        self.interf_ff.push(rec);
    }

    fn append_interf(arr: &mut Vec<BopdsInterf>, i1: usize, i2: usize, index_new: Option<usize>) {
        let mut it = BopdsInterf::new(i1, i2);
        if let Some(n) = index_new {
            it.set_index_new(n);
        }
        arr.push(it);
    }

    /// Sets `index_new` on the typed record for `(i1, i2)` if present.
    /// Source: `BOPDS_Interf::SetIndexNew` in `PerformNewVertices`.
    fn bind_index_new(arr: &mut [BopdsInterf], i1: usize, i2: usize, index_new: usize) -> bool {
        for it in arr.iter_mut() {
            if it.contains(i1) && it.contains(i2) {
                it.set_index_new(index_new);
                return true;
            }
        }
        false
    }

    /// Binds the new vertex of an E/E interference (`InterfEE::SetIndexNew`).
    pub fn bind_ee_new_vertex(&mut self, n_e1: usize, n_e2: usize, n_v: usize) -> bool {
        Self::bind_index_new(&mut self.interf_ee, n_e1, n_e2, n_v)
    }

    /// Binds the new vertex of an E/F interference (`InterfEF::SetIndexNew`).
    pub fn bind_ef_new_vertex(&mut self, n_e: usize, n_f: usize, n_v: usize) -> bool {
        Self::bind_index_new(&mut self.interf_ef, n_e, n_f, n_v)
    }

    /// Binds the new vertex of a V/F interference (`InterfVF::SetIndexNew`).
    pub fn bind_vf_new_vertex(&mut self, n_v: usize, n_f: usize, n_new: usize) -> bool {
        Self::bind_index_new(&mut self.interf_vf, n_v, n_f, n_new)
    }

    /// Store `IntTools_CommonPrt::Range1` on the E/E record of `(e1, e2)`.
    pub fn set_ee_common_range(&mut self, n_e1: usize, n_e2: usize, first: f64, last: f64) -> bool {
        Self::bind_common_range(&mut self.interf_ee, n_e1, n_e2, first, last)
    }

    /// Store `IntTools_CommonPrt::Range1` on the E/F record of `(e, f)`.
    pub fn set_ef_common_range(&mut self, n_e: usize, n_f: usize, first: f64, last: f64) -> bool {
        Self::bind_common_range(&mut self.interf_ef, n_e, n_f, first, last)
    }

    fn bind_common_range(
        arr: &mut [BopdsInterf],
        i1: usize,
        i2: usize,
        first: f64,
        last: f64,
    ) -> bool {
        for it in arr.iter_mut() {
            if it.contains(i1) && it.contains(i2) {
                it.set_common_range(first, last);
                return true;
            }
        }
        false
    }

    /// Redirects the new-vertex index of every typed interference to its
    /// same-domain (SD) representative.
    ///
    /// Source: `BOPAlgo_PaveFiller::UpdateInterfsWithSDVertices` +
    /// `UpdateIntfsWithSDVertices` (`BOPAlgo_PaveFiller_10.cxx`): each typed
    /// interference whose `index_new` names a vertex that has since been merged
    /// into an SD cluster is re-pointed at the cluster representative. Called
    /// in `PerformInternal` after the E/F stage and again after `MakeBlocks`.
    pub fn update_interfs_with_sd_vertices(&mut self) {
        Self::update_intfs_with_sd_vertices(&self.shapes_sd, &mut self.interf_vv);
        Self::update_intfs_with_sd_vertices(&self.shapes_sd, &mut self.interf_ve);
        Self::update_intfs_with_sd_vertices(&self.shapes_sd, &mut self.interf_vf);
        Self::update_intfs_with_sd_vertices(&self.shapes_sd, &mut self.interf_ee);
        Self::update_intfs_with_sd_vertices(&self.shapes_sd, &mut self.interf_ef);
    }

    fn update_intfs_with_sd_vertices(sd: &HashMap<usize, usize>, interfs: &mut [BopdsInterf]) {
        for it in interfs.iter_mut() {
            let Some(n) = it.get_index_new() else { continue };
            // `BOPDS_DS::HasShapeSD` walks the same-domain chain and returns
            // the last representative; skip when the index has no SD partner.
            let mut cur = n;
            let mut has_sd = false;
            while let Some(&next) = sd.get(&cur) {
                cur = next;
                has_sd = true;
            }
            if has_sd {
                it.set_index_new(cur);
            }
        }
    }
}
