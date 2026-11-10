use super::super::prelude::*;
use super::super::*;

impl BopBuilder {
    /// `FillImagesVertices` ... `FillImagesSolids` (`BOPAlgo_Builder.cxx:336-418`).
    /// No object/tool IN/OUT filter.
    pub(crate) fn fill_images(&mut self) -> Result<(), String> {
        self.fill_images_vertices()?;
        self.fill_images_edges()?;
        crate::bop_build_common::fill_images_containers(self, ShapeType::Wire)?;
        self.fill_images_faces()?;
        crate::bop_build_common::fill_images_containers(self, ShapeType::Shell)?;
        self.fill_images_solids()?;
        crate::bop_build_common::fill_images_containers(self, ShapeType::CompSolid)?;
        Ok(())
    }

    /// Fills the images of the vertices.
    ///
    /// Port of `BOPAlgo_Builder::FillImagesVertices`: for every same-domain
    /// vertex pair `(v, v_sd)` of the data structure, records `v_sd` as the
    /// image of `v` in the images table, binds the same-domain map and adds
    /// the origin back-map entry `v_sd → [v]`.
    pub(crate) fn fill_images_vertices(&mut self) -> Result<(), String> {
        let sd = self.filler.ds().shapes_sd().clone();
        for (&n_v, &n_vsd) in &sd {
            let Some(v) = self.filler.ds().shape(n_v).cloned() else { continue };
            let Some(vsd) = self.filler.ds().shape(n_vsd).cloned() else { continue };
            self.history.add_image(&v, vsd.clone());
            self.shapes_sd.insert(n_v, n_vsd);
            self.origins.entry(shape_key(&vsd)).or_default().push(v.clone());
        }
        Ok(())
    }

    /// Fills the images of the edges.
    ///
    /// Port of `BOPAlgo_Builder::FillImagesEdges`: every source edge with pave
    /// blocks expands into the split edges of its blocks. The split edges were
    /// built by the intersection phase ([`crate::pave_blocks::blocks_make_split_edges`]
    /// → [`crate::algo_tools::AlgoTools::make_split_edge`]) and stored in the
    /// data structure; each block's `edge()` points at its split edge. The
    /// split edge is recorded as the image of the source edge, and the origin
    /// back-map gets the source edge. A common-block member uses
    /// `RealPaveBlock` (`BOPDS_DS::RealPaveBlock`) so the image is the shared
    /// split, matching `BOPAlgo_Builder::FillImagesEdges` (`_1.cxx:71-126`).
    pub(crate) fn fill_images_edges(&mut self) -> Result<(), String> {
        crate::pave_blocks::blocks_make_split_edges(&mut self.filler)?;
        let n = self.filler.ds().nb_source_shapes();
        for i in 0..n {
            let Some(si) = self.filler.ds().shape_info(i).cloned() else { continue };
            if si.shape_type() != ShapeType::Edge {
                continue;
            }
            if !self.filler.ds().has_pave_blocks(i) {
                continue;
            }
            let e = si.shape().clone();
            let blocks = self.filler.ds().pave_blocks(i).to_vec();
            for pb in &blocks {
                // A zero-length block (two coincident bound vertices) is not a
                // real interval and must not contribute an image — otherwise the
                // split image list gains a spurious full-edge entry.
                if (pb.last - pb.first).abs() <= occt_core::precision::PCONFUSION {
                    continue;
                }
                let pbr = self.filler.ds().real_pave_block(pb);
                let n_sp_r = pbr.edge();
                if n_sp_r == 0 {
                    continue;
                }
                let Some(sp_r) = self.filler.ds().shape(n_sp_r).cloned() else { continue };
                self.history.add_image(&e, sp_r.clone());
                self.origins.entry(shape_key(&sp_r)).or_default().push(e.clone());
                if is_common_block_on_edge(self.filler.ds(), pb) {
                    let n_sp = pb.edge();
                    if n_sp != 0 && n_sp != n_sp_r {
                        self.shapes_sd.insert(n_sp, n_sp_r);
                    }
                }
            }
        }
        Ok(())
    }

    /// Fills the images of the faces — split faces of the arguments.
    ///
    /// Delegates to `crate::bop_build_faces::fill_images_faces`
    /// (`BuildSplitFaces` → `FillSameDomainFaces`). GF only; no obj/tool state.
    pub(crate) fn fill_images_faces(&mut self) -> Result<(), String> {
        crate::bop_build_faces::fill_images_faces(self)
    }

    /// Fills the images of the solids — split solids of the arguments.
    ///
    /// FillIn3DParts → BuildSplitSolids → FillInternalShapes
    /// (`BOPAlgo_Builder_3.cxx:70`). GF only; no obj/tool state.
    ///
    /// Public because `BOPAlgo_Builder::FillImagesSolids` is `protected` and
    /// subclasses call it directly: `BRepFeat_Builder::PerformResult`
    /// (`BRepFeat_Builder.cxx:240`) runs the stage a second time on the
    /// retained tool parts.
    pub fn fill_images_solids(&mut self) -> Result<(), String> {
        if !crate::bop_images_solids::has_source_solids(self) {
            self.in_parts.clear();
            return Ok(());
        }
        let mut ctx = crate::int_tools_full::IntToolsContext::new();
        let fill = crate::bop_fill_in3d::fill_in_3d_parts_builder(self, &mut ctx);
        self.in_parts = fill.in_parts.clone();
        crate::bop_builder_split_solids::build_split_solids_occt(self, &fill)?;
        crate::bop_builder_fill_internals::fill_internal_shapes_occt(self)?;
        Ok(())
    }
}
