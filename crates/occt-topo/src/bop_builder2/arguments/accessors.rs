use super::super::prelude::*;
use super::super::*;

impl BopBuilder {
    // -----------------------------------------------------------------------
    // Data-structure / history access
    // -----------------------------------------------------------------------

    /// Returns the data structure of the algorithm (through the filler).
    pub fn ds(&self) -> &BopdsDS {
        self.filler.ds()
    }

    /// Returns the mutable data structure of the algorithm.
    pub fn ds_mut(&mut self) -> &mut BopdsDS {
        self.filler.ds_mut()
    }

    /// Returns the images/history naming table.
    pub fn history(&self) -> &BopHistory {
        &self.history
    }

    /// Returns the mutable images/history naming table.
    pub fn history_mut(&mut self) -> &mut BopHistory {
        &mut self.history
    }

    /// Returns the origins back-map (`split shape key → source shapes`).
    ///
    /// `myOrigins` is keyed by `TopTools_ShapeMapHasher`
    /// (`BOPAlgo_Builder.hxx:501`), i.e. [`bop_occt_util::shape_key`].
    pub fn origins(&self) -> &HashMap<ShapeMapKey, Vec<TopoShape>> {
        &self.origins
    }

    /// Returns the mutable origins back-map.
    pub fn origins_mut(&mut self) -> &mut HashMap<ShapeMapKey, Vec<TopoShape>> {
        &mut self.origins
    }

    /// Returns the same-domain shapes map of the builder.
    pub fn shapes_sd(&self) -> &HashMap<usize, usize> {
        &self.shapes_sd
    }

    /// The additional tolerance of the operation (from the PaveFiller).
    pub fn fuzzy_value(&self) -> f64 {
        self.filler.fuzzy_value()
    }

    /// Sets the additional tolerance forwarded to the PaveFiller.
    pub fn set_fuzzy_value(&mut self, v: f64) {
        self.filler.set_fuzzy_value(v);
    }

    /// `BOPAlgo_Options::SetGlue(theGlue)` (`BOPAlgo_Options.lxx`): the gluing
    /// option the intersection phase reads back through
    /// [`PaveFiller::glue`](crate::pave_filler::PaveFiller::glue).
    ///
    /// Call sites: `BRepFill_Draft.cxx:752` (`BOPAlgo_GlueShift`).
    pub fn set_glue(&mut self, glue: crate::pave_filler::GlueEnum) {
        self.filler.set_glue(glue);
    }

    /// True when the PaveFiller runs in non-destructive mode.
    pub fn non_destructive(&self) -> bool {
        self.filler.non_destructive()
    }

    /// `BOPAlgo_Options::SetNonDestructive` (`BOPAlgo_Options.lxx`).
    pub fn set_non_destructive(&mut self, b: bool) {
        self.filler.set_non_destructive(b);
    }

    /// `BOPAlgo_Options::Glue` (`BOPAlgo_Options.lxx`).
    pub fn glue(&self) -> crate::pave_filler::GlueEnum {
        self.filler.glue()
    }

    /// Binds `shape` as same-domain with `sd`, mirroring `myShapesSD.Bind`.
    ///
    /// The shapes are resolved to their data-structure indices and the binding
    /// is recorded. OCCT `myShapesSD.Bind(aF, *pFSD)` also binds the
    /// representative to itself (`IsBound` is true); `BuildDraftSolid` uses
    /// that to take the `IsSplitToReverseWithWarn` branch. Reconstructed faces
    /// that are not in the DS are still recorded by TShape key.
    pub fn bind_shapes_sd(&mut self, shape: TopoShape, sd: TopoShape) {
        let k = GeometryRegistry::shape_key(&shape);
        self.sd_by_key.insert(k, sd.clone());
        if let (Some(i), Some(j)) = (self.filler.ds().index(&shape), self.filler.ds().index(&sd)) {
            if i != j {
                self.shapes_sd.insert(i, j);
            }
        }
    }

    /// Returns the same-domain representative of `shape`, when bound
    /// (`myShapesSD.Seek`), following the SD chain.
    pub fn seek_shapes_sd(&self, shape: &TopoShape) -> Option<TopoShape> {
        if let Some(sd) = self.sd_by_key.get(&GeometryRegistry::shape_key(shape)) {
            return Some(sd.clone());
        }
        let start = self.filler.ds().index(shape)?;
        let mut cur = start;
        let mut guard = 0;
        while let Some(&j) = self.shapes_sd.get(&cur) {
            cur = j;
            guard += 1;
            if guard > 64 {
                break;
            }
        }
        if cur == start {
            None
        } else {
            self.filler.ds().shape(cur).cloned()
        }
    }

    /// Returns the shape produced by the last run.
    pub fn result(&self) -> &TopoShape {
        &self.result_shape
    }

    /// Replaces the result compound (`myShape = ...` after BuildRC / BuildBOP).
    pub fn set_result_shape(&mut self, shape: TopoShape) {
        self.result_shape = shape;
    }

    /// `myInParts` filled by FillIn3DParts.
    pub fn in_parts(&self) -> &HashMap<ShapeMapKey, Vec<TopoShape>> {
        &self.in_parts
    }

    /// Mutable `myInParts`.
    ///
    /// `BRepFeat_Builder::FillIn3DParts` (`BRepFeat_Builder.cxx:801-826`)
    /// overrides the virtual stage only to drop the removed faces from the IN
    /// parts before `BuildSplitSolids` consumes them.
    pub fn in_parts_mut(&mut self) -> &mut HashMap<ShapeMapKey, Vec<TopoShape>> {
        &mut self.in_parts
    }

    /// Returns the PaveFiller of the algorithm.
    pub fn filler(&self) -> &PaveFiller {
        &self.filler
    }

    /// Returns the mutable PaveFiller of the algorithm.
    pub fn filler_mut(&mut self) -> &mut PaveFiller {
        &mut self.filler
    }

    // -----------------------------------------------------------------------
    // History flags (`BOPAlgo_BuilderShape`, `BOPAlgo_BuilderShape.hxx:114-141`)
    // -----------------------------------------------------------------------

    /// Enables/disables the history collection
    /// (`BOPAlgo_BuilderShape::SetToFillHistory`).
    pub fn set_to_fill_history(&mut self, flag: bool) {
        self.fill_history = flag;
    }

    /// Returns the flag of history availability
    /// (`BOPAlgo_BuilderShape::HasHistory`).
    pub fn has_history(&self) -> bool {
        self.fill_history
    }

    /// Installs an already performed filler, so the building phase can run
    /// without repeating the intersection (`PerformWithFiller`).
    pub fn set_filler(&mut self, filler: PaveFiller) {
        self.filler = filler;
    }
}
