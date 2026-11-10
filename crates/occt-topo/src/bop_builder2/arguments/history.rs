use super::super::prelude::*;
use super::super::*;

impl BopBuilder {
    /// Post-treats the result shape by correcting the tolerances.
    ///
    /// Port of `BOPAlgo_Builder::PostTreat` (`BOPAlgo_Builder.cxx:450`):
    /// `CorrectTolerances` + `CorrectShapeTolerances` with `aMA` filled from
    /// source VERTEX/EDGE/FACE when `NonDestructive` is set.
    pub(crate) fn post_treat(&mut self) -> Result<(), String> {
        let mut avoid = HashSet::new();
        if self.non_destructive() {
            let n = self.filler.ds().nb_source_shapes();
            for i in 0..n {
                let Some(si) = self.filler.ds().shape_info(i) else {
                    continue;
                };
                let t = si.shape_type();
                if t == ShapeType::Vertex || t == ShapeType::Edge || t == ShapeType::Face {
                    avoid.insert(GeometryRegistry::shape_key(si.shape()));
                }
            }
        }
        AlgoTools::correct_tolerances_avoid(&self.result_shape, &avoid, 0.05);
        Ok(())
    }

    /// Fills the modified/generated/removed relations of the history from the
    /// images table and the result shape.
    ///
    /// Port of `BOPAlgo_Builder::PrepareHistory` (`BOPAlgo_Builder_4.cxx`):
    /// for every source shape of the data structure,
    /// - the split pieces of the shape kept in the result become **modified**
    ///   from it;
    /// - the vertices/edges the intersections created from an EDGE/FACE source
    ///   become **generated** from it ([`BopBuilder::loc_generated`]);
    /// - a shape with no trace in the result (and no surviving splits) is
    ///   marked **removed**.
    pub(crate) fn prepare_history(&mut self) {
        // `BOPAlgo_Builder_4.cxx:143-148`: nothing to do when the history is
        // switched off; otherwise the relations are rebuilt from scratch while
        // the `myImages` naming table is kept as `LocModified`.
        if !self.fill_history {
            return;
        }
        self.history.clear_relations();
        // All shapes of the result (the root and every sub-shape), keyed by
        // (TShape, orientation) — the OCCT `TopTools_MapOfShape` identity.
        let result_keys: HashSet<(usize, u8)> = all_subshapes(&self.result_shape)
            .iter()
            .map(|s| (Arc::as_ptr(&s.tshape) as usize, s.orientation() as u8))
            .collect();
        let in_result = |s: &TopoShape| {
            result_keys.contains(&(Arc::as_ptr(&s.tshape) as usize, s.orientation() as u8))
        };

        let n = self.filler.ds().nb_source_shapes();
        for i in 0..n {
            let Some(si) = self.filler.ds().shape_info(i) else { continue };
            let s = si.shape().clone();
            if !is_supported_type(&s) {
                continue;
            }

            let mut is_modified = false;
            // Modified: the splits of the shape kept in the result, oriented
            // like the source (VERTEX/SOLID take the source orientation; an
            // EDGE/FACE whose direction flips is reversed).
            let splits: Vec<TopoShape> =
                self.history.image(&s).map(|v| v.to_vec()).unwrap_or_default();
            for sp in &splits {
                if !in_result(sp) {
                    continue;
                }
                let mut sp = sp.clone();
                let t = sp.shape_type();
                if t == ShapeType::Vertex || t == ShapeType::Solid {
                    sp.set_orientation(s.orientation());
                } else if is_split_to_reverse_ctxless(&sp, &s) {
                    sp.set_orientation(sp.orientation().reversed());
                }
                let _ = self.history.add_modified(&s, sp);
                is_modified = true;
            }

            // Generated: the vertices/edges the intersections created from the
            // shape, kept in the result.
            for g in self.loc_generated(&s) {
                if in_result(&g) {
                    let _ = self.history.add_generated(&s, g);
                }
            }

            // Removed: the shape has no trace in the result nor any surviving
            // split.
            if !is_modified && !in_result(&s) {
                let _ = self.history.add_removed(s);
            }
        }
    }

    /// Shapes generated from `s` by the intersections: the vertices created in
    /// E/E and E/F interferences (for an EDGE or FACE source) plus, for a
    /// FACE, the section edges and vertices lying on it.
    ///
    /// Port of `BOPAlgo_Builder::LocGenerated` (`BOPAlgo_Builder_4.cxx`). The
    /// new-vertex part reads the `index_new` records of the E/E and E/F
    /// interferences; the section edges/vertices of a face come from its
    /// [`crate::bopds::BopdsFaceInfo`]. Whether a returned shape actually made
    /// it into the result is checked by the caller
    /// ([`BopBuilder::prepare_history`]).
    pub(crate) fn loc_generated(&self, s: &TopoShape) -> Vec<TopoShape> {
        let mut out: Vec<TopoShape> = Vec::new();
        let kind = s.shape_type();
        if kind != ShapeType::Edge && kind != ShapeType::Face {
            return out;
        }
        let ds = self.filler.ds();
        let Some(n_s) = ds.index(s) else { return out };
        // Untouched shapes carry no generated elements — an edge without pave
        // blocks, a face without a face-info entry (the OCCT `HasReference`
        // guard).
        let is_face = kind == ShapeType::Face;
        if is_face {
            if !ds.face_info_pool().iter().any(|fi| fi.face_index == n_s) {
                return out;
            }
        } else if !ds.has_pave_blocks(n_s) {
            return out;
        }

        // New vertices of the E/E (edge sources) and E/F interferences
        // containing the shape, deduplicated by their same-domain index.
        let mut fence: Vec<usize> = Vec::new();
        if !is_face {
            Self::collect_interf_vertices(ds, ds.interf_ee(), n_s, &mut fence, &mut out);
        }
        Self::collect_interf_vertices(ds, ds.interf_ef(), n_s, &mut fence, &mut out);
        if !is_face {
            return out;
        }

        // Section edges and section vertices lying on the face.
        if let Some(fi) = ds.face_info_pool().iter().find(|fi| fi.face_index == n_s) {
            for &(e_idx, _, _) in fi.paves() {
                if let Some(e) = ds.shape(e_idx) {
                    out.push(e.clone());
                }
            }
            for &(v_idx, _, _) in fi.verts() {
                if let Some(v) = ds.shape(v_idx) {
                    out.push(v.clone());
                }
            }
        }
        out
    }

    /// Appends the `index_new` vertices of the interferences `ints` containing
    /// `n_s`, deduplicated by their same-domain vertex index (the OCCT
    /// `LocGenerated` loop over `InterfEE` / `InterfEF`).
    pub(crate) fn collect_interf_vertices(
        ds: &BopdsDS,
        ints: &[BopdsInterf],
        n_s: usize,
        fence: &mut Vec<usize>,
        out: &mut Vec<TopoShape>,
    ) {
        for it in ints {
            let Some(n_v) = it.new_index() else { continue };
            if !it.contains(n_s) {
                continue;
            }
            // Resolve the vertex through its same-domain twin.
            let n_v = ds.has_shape_sd(n_v).unwrap_or(n_v);
            if fence.contains(&n_v) {
                continue;
            }
            fence.push(n_v);
            if let Some(v) = ds.shape(n_v) {
                out.push(v.clone());
            }
        }
    }
}
