//! `BOPAlgo_Section` — the result builder of the SECTION operation.
//!
//! Source: `BOPAlgo_Section.cxx:37-412` (TKBO/BOPAlgo). `BOPAlgo_Section` is a
//! `BOPAlgo_Builder` whose `PerformInternal1` stops after the VERTEX and EDGE
//! FillImages/BuildResult pair and then assembles the section result from the
//! face-info `Sc` maps and the pave-block pool ([`BopBuilder::build_section`]).
use super::prelude::*;
use super::*;
use crate::topo_tools::{explorer_shapes, explorer_shapes_avoiding, shape_map_key, ShapeMapKey};

/// `NCollection_IndexedMap<TopoDS_Shape>` keeping OCCT's insertion order:
/// `BuildSection` walks these maps by index, so the order is part of the result.
struct OrderedShapeMap {
    keys: Vec<TopoShape>,
    index: HashMap<ShapeMapKey, usize>,
}

impl OrderedShapeMap {
    fn new() -> Self {
        Self {
            keys: Vec::new(),
            index: HashMap::new(),
        }
    }

    /// `Add` with the `IsSame` fence of `TopTools_ShapeMapHasher`; `true` when
    /// the shape was not there yet.
    fn add(&mut self, s: &TopoShape) -> bool {
        let k = shape_map_key(s);
        if self.index.contains_key(&k) {
            return false;
        }
        self.index.insert(k, self.keys.len());
        self.keys.push(s.clone());
        true
    }

    fn extent(&self) -> usize {
        self.keys.len()
    }

    fn at(&self, i: usize) -> &TopoShape {
        &self.keys[i]
    }
}

/// `NCollection_IndexedDataMap<TopoDS_Shape, int>` — the `aMSI` count map of
/// `BuildSection` step 3.
struct OrderedCountMap {
    keys: Vec<TopoShape>,
    index: HashMap<ShapeMapKey, usize>,
    counts: Vec<i32>,
}

impl OrderedCountMap {
    fn new() -> Self {
        Self {
            keys: Vec::new(),
            index: HashMap::new(),
            counts: Vec::new(),
        }
    }

    fn add(&mut self, s: &TopoShape, n: i32) {
        let k = shape_map_key(s);
        match self.index.get(&k) {
            Some(&i) => self.counts[i] += n,
            None => {
                self.index.insert(k, self.keys.len());
                self.keys.push(s.clone());
                self.counts.push(n);
            }
        }
    }

    /// `++aMSI.ChangeFromKey(aS)` / `aMSI.Add(aS, 1)`.
    fn bump(&mut self, s: &TopoShape) {
        match self.index.get(&shape_map_key(s)) {
            Some(&i) => self.counts[i] += 1,
            None => self.add(s, 1),
        }
    }

    fn extent(&self) -> usize {
        self.keys.len()
    }

    fn key_at(&self, i: usize) -> &TopoShape {
        &self.keys[i]
    }

    fn count_at(&self, i: usize) -> i32 {
        self.counts[i]
    }
}

/// `NCollection_IndexedDataMap<TopoDS_Shape, List<TopoDS_Shape>>` — the `aMVE`
/// vertex-to-edges map of `BuildSection` step 4.
struct OrderedAncestorMap {
    keys: Vec<TopoShape>,
    index: HashMap<ShapeMapKey, usize>,
    lists: Vec<Vec<TopoShape>>,
}

impl OrderedAncestorMap {
    fn new() -> Self {
        Self {
            keys: Vec::new(),
            index: HashMap::new(),
            lists: Vec::new(),
        }
    }

    fn slot(&mut self, s: &TopoShape) -> usize {
        let k = shape_map_key(s);
        if let Some(&i) = self.index.get(&k) {
            return i;
        }
        let i = self.keys.len();
        self.index.insert(k, i);
        self.keys.push(s.clone());
        self.lists.push(Vec::new());
        i
    }

    /// `TopExp::MapShapesAndAncestors` (`TopExp.cxx:80-120`) appending into an
    /// existing map, exactly like the third argument of the OCCT signature.
    fn map_shapes_and_ancestors(&mut self, shape: &TopoShape, ts: ShapeType, ta: ShapeType) {
        for anc in explorer_shapes(shape, ta) {
            for sub in explorer_shapes(&anc, ts) {
                let i = self.slot(&sub);
                self.lists[i].push(anc.clone());
            }
        }
        for sub in explorer_shapes_avoiding(shape, ts, ta) {
            self.slot(&sub);
        }
    }

    fn extent(&self) -> usize {
        self.keys.len()
    }

    fn key_at(&self, i: usize) -> &TopoShape {
        &self.keys[i]
    }

    fn list_at(&self, i: usize) -> &[TopoShape] {
        &self.lists[i]
    }
}

impl BopBuilder {
    /// Runs the SECTION operation on `objects` and `tools`.
    ///
    /// `BRepAlgoAPI_BooleanOperation::Build` (`BRepAlgoAPI_BooleanOperation.cxx:174-205`)
    /// intersects `Objects + Tools` together and hands the filler to
    /// `BOPAlgo_Section`; `PerformInternal1` (`BOPAlgo_Section.cxx:100-164`) then
    /// runs `CheckData` / `Prepare` / `FillImagesVertices` / `BuildResult(VERTEX)`
    /// / `FillImagesEdges` / `BuildResult(EDGE)` / `BuildSection` /
    /// `PrepareHistory` / `PostTreat`.
    pub fn perform_section(
        &mut self,
        objects: &[TopoShape],
        tools: &[TopoShape],
    ) -> Result<TopoShape, String> {
        self.reset_run();
        // `CheckData` (`BOPAlgo_Section.cxx:57-70`): a non-empty argument list.
        if self.arguments.is_empty() {
            let msg = "BOPAlgo_Section: too few arguments".to_string();
            self.add_error(msg.clone());
            return Err(msg);
        }
        self.objects = objects.to_vec();
        self.tools = tools.to_vec();
        // Intersection phase (`BRepAlgoAPI_BooleanOperation::IntersectShapes`).
        self.filler.set_arguments(&self.arguments);
        if let Err(e) = self.filler.perform() {
            let msgs: Vec<String> = self.filler.errors().to_vec();
            for m in msgs {
                self.add_error(m);
            }
            if self.errors.is_empty() {
                self.add_error(e.clone());
            }
            return Err(self.errors.first().cloned().unwrap_or_else(|| e));
        }
        if self.filler.has_errors() {
            let msgs: Vec<String> = self.filler.errors().to_vec();
            for m in msgs {
                self.add_error(m);
            }
            return Err(self.errors.first().cloned().unwrap_or_default());
        }
        // 2. Prepare
        self.prepare();
        // 3.1 Vertices
        self.fill_images_vertices()?;
        self.build_result_of_type(ShapeType::Vertex)?;
        // 3.2 Edges
        self.fill_images_edges()?;
        self.build_result_of_type(ShapeType::Edge)?;
        // 4. Section
        self.build_section();
        // 5. History / 6. Post-treatment
        self.prepare_history();
        self.post_treat()?;
        Ok(self.result_shape.clone())
    }

    /// `BOPAlgo_Builder::Prepare` (`BOPAlgo_Builder.cxx:155-163`): `myShape`
    /// restarts as an empty compound that the following `BuildResult` calls
    /// append into.
    pub(super) fn prepare(&mut self) {
        self.result_shape = TopoBuilder::new().make_shape(ShapeType::Compound);
    }

    /// `BOPAlgo_Builder::BuildResult(theType)` (`BOPAlgo_Builder_1.cxx:130-167`):
    /// appends the images of every argument of `t` — or the argument itself when
    /// it has no image — to `myShape`. The fence is local to the call, so the
    /// two calls of `PerformInternal1` do not share it.
    pub(super) fn build_result_of_type(&mut self, t: ShapeType) -> Result<(), String> {
        let b = TopoBuilder::new();
        let args = self.arguments.clone();
        let mut fence: Vec<ShapeMapKey> = Vec::new();
        for arg in &args {
            if arg.shape_type() != t {
                continue;
            }
            let imgs: Vec<TopoShape> = match self.history.image(arg) {
                Some(list) => list.to_vec(),
                None => vec![arg.clone()],
            };
            for img in imgs {
                let k = shape_map_key(&img);
                if fence.contains(&k) {
                    continue;
                }
                fence.push(k);
                b.add(&mut self.result_shape, &img);
            }
        }
        Ok(())
    }

    /// `BOPAlgo_Section::BuildSection` (`BOPAlgo_Section.cxx:167-412`).
    ///
    /// 1. `aRC1` = every section vertex, every IN vertex that is new or
    ///    interfered, and every section edge of every source face;
    /// 2. plus the `PaveBlock1` edge of every common block that carries faces;
    /// 3. plus, for every argument, the vertex/edge content of its images
    ///    counted into `aMSI` (the "boundary of the arguments" pass);
    /// 4. `aRC` = the section edges reached from `aRC1`'s vertices, or the
    ///    vertex itself when it has no edge, fenced once.
    pub(super) fn build_section(&mut self) {
        let b = TopoBuilder::new();
        let mut rc1 = b.make_shape(ShapeType::Compound);
        let n = self.ds().nb_source_shapes();
        for i in 0..n {
            let Some(si) = self.ds().shape_info(i) else { continue };
            if si.shape_type() != ShapeType::Face {
                continue;
            }
            let Some(fi) = self.ds().face_info(i) else { continue };
            let verts_sc: Vec<usize> = fi.verts_sc().to_vec();
            let verts_in: Vec<usize> = fi.verts_in().to_vec();
            let paves_sc: Vec<(usize, f64, f64)> = fi.paves().to_vec();
            // 1.1 Vertices that are section vertices
            for n_v in verts_sc {
                if let Some(v) = self.ds().shape(n_v) {
                    b.add(&mut rc1, v);
                }
            }
            // 1.2 Vertices that are in a face
            for n_v in verts_in {
                if n_v == 0 {
                    continue;
                }
                if self.ds().is_new_shape(n_v) || self.ds().has_interf(n_v) {
                    if let Some(v) = self.ds().shape(n_v) {
                        b.add(&mut rc1, v);
                    }
                }
            }
            // 1.3 Section edges — `PaveBlocksSc` is an indexed map, so a block
            // repeated inside the face-info list is added once.
            let mut seen: Vec<(usize, u64, u64)> = Vec::new();
            for (n_e, first, last) in paves_sc {
                let key = (n_e, first.to_bits(), last.to_bits());
                if seen.contains(&key) {
                    continue;
                }
                seen.push(key);
                if let Some(e) = self.ds().shape(n_e) {
                    b.add(&mut rc1, e);
                }
            }
        }
        // 2. Common blocks between an edge and a face
        let mut cb_edges: Vec<usize> = Vec::new();
        for list in self.ds().pave_blocks_pool() {
            for pb in list {
                let Some(cb) = self.ds().common_block(pb) else { continue };
                if cb.faces().is_empty() {
                    continue;
                }
                let Some(pbr) = cb.pave_block1() else { continue };
                let n_e = pbr.edge();
                if !cb_edges.contains(&n_e) {
                    cb_edges.push(n_e);
                }
            }
        }
        for n_e in cb_edges {
            if let Some(e) = self.ds().shape(n_e) {
                b.add(&mut rc1, e);
            }
        }

        // 3. Treatment of the boundaries of the arguments
        let mut lsa: Vec<TopoShape> = Vec::new();
        for a in &self.arguments {
            if !lsa.iter().any(|x| is_same(x, a)) {
                lsa.push(a.clone());
            }
        }
        let mut msi = OrderedCountMap::new();
        for s in &lsa {
            let mut ls: Vec<TopoShape> = Vec::new();
            let mut fence = OrderedShapeMap::new();
            for e in explorer_shapes(s, ShapeType::Edge) {
                if fence.add(&e) {
                    ls.push(e);
                }
            }
            for v in explorer_shapes(s, ShapeType::Vertex) {
                if fence.add(&v) {
                    ls.push(v);
                }
            }
            // 3.2 aMSI
            let mut ms = OrderedShapeMap::new();
            for sub in &ls {
                match self.history.image(sub) {
                    Some(list) => {
                        for sim in list {
                            for v in explorer_shapes(sim, ShapeType::Vertex) {
                                ms.add(&v);
                            }
                            for e in explorer_shapes(sim, ShapeType::Edge) {
                                ms.add(&e);
                            }
                        }
                    }
                    None => {
                        for v in explorer_shapes(sub, ShapeType::Vertex) {
                            ms.add(&v);
                        }
                        for e in explorer_shapes(sub, ShapeType::Edge) {
                            ms.add(&e);
                        }
                    }
                }
            }
            for j in 0..ms.extent() {
                let s_j = ms.at(j).clone();
                msi.bump(&s_j);
            }
        }

        // 4. Build the result
        let mut mve = OrderedAncestorMap::new();
        mve.map_shapes_and_ancestors(&rc1, ShapeType::Vertex, ShapeType::Edge);
        for i in 0..msi.extent() {
            if msi.count_at(i) > 1 {
                let v = msi.key_at(i).clone();
                mve.map_shapes_and_ancestors(&v, ShapeType::Vertex, ShapeType::Edge);
            }
        }
        let mut rc = b.make_shape(ShapeType::Compound);
        let mut fence = OrderedShapeMap::new();
        for i in 0..mve.extent() {
            let v = mve.key_at(i).clone();
            let le = mve.list_at(i).to_vec();
            if le.is_empty() {
                // alone vertices
                if fence.add(&v) {
                    b.add(&mut rc, &v);
                }
            } else {
                for e in le {
                    if fence.add(&e) {
                        b.add(&mut rc, &e);
                    }
                }
            }
        }
        self.result_shape = rc;
    }
}
