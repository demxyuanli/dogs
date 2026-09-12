//! `BOPAlgo_CellsBuilder` — take / avoid split parts of a General Fuse.
//!
//! Source: `BOPAlgo_CellsBuilder.cxx`. After GF, every split part is indexed
//! to the original arguments that contain it. [`CellsBuilder::add_to_result`]
//! keeps parts that belong to every shape in `take` and to none in `avoid`.

use std::collections::{HashMap, HashSet};

use crate::algo_tools_range::dimension;
use crate::bop_bop::type_to_explore;
use crate::bop_builder2::BopBuilder;
use crate::bop_occt_util::{explore, iter_children, shape_key, treat_compound};
use crate::builder::TopoBuilder;
use crate::shape::TopoShape;

/// `BOPAlgo_CellsBuilder`.
#[derive(Debug, Clone)]
pub struct CellsBuilder {
    inner: BopBuilder,
    /// `myIndex`: split part TShape key → original arguments that contain it.
    index: HashMap<usize, Vec<TopoShape>>,
    /// `myAllParts` compound.
    all_parts: TopoShape,
    /// `myShapeMaterial`: part key → material id.
    shape_material: HashMap<usize, i32>,
}

impl CellsBuilder {
    /// Empty builder.
    pub fn new() -> Self {
        Self {
            inner: BopBuilder::new(),
            index: HashMap::new(),
            all_parts: TopoBuilder::new().make_compound_of(&[]).0,
            shape_material: HashMap::new(),
        }
    }

    /// `SetArguments`.
    pub fn set_arguments(&mut self, shapes: &[TopoShape]) {
        self.inner.set_arguments(shapes);
    }

    /// `SetFuzzyValue`.
    pub fn set_fuzzy_value(&mut self, fuzzy: f64) {
        self.inner.set_fuzzy_value(fuzzy);
    }

    /// GF split + [`index_parts`] + empty result (`PerformInternal1`).
    pub fn perform(&mut self) -> Result<(), String> {
        self.inner.perform()?;
        self.index_parts();
        self.remove_all_from_result();
        Ok(())
    }

    /// `GetAllParts`.
    pub fn all_parts(&self) -> &TopoShape {
        &self.all_parts
    }

    /// Current result compound (`Shape()`).
    pub fn shape(&self) -> &TopoShape {
        self.inner.result()
    }

    /// `RemoveAllFromResult`.
    pub fn remove_all_from_result(&mut self) {
        self.inner
            .set_result_shape(TopoBuilder::new().make_compound_of(&[]).0);
        self.shape_material.clear();
    }

    /// `AddAllToResult`.
    pub fn add_all_to_result(&mut self, material: i32) {
        self.shape_material.clear();
        self.inner.set_result_shape(self.all_parts.clone());
        if material != 0 {
            for p in iter_children(&self.all_parts) {
                self.shape_material.insert(shape_key(&p), material);
            }
        }
    }

    /// `AddToResult`.
    pub fn add_to_result(&mut self, take: &[TopoShape], avoid: &[TopoShape], material: i32) {
        let parts = self.find_parts(take, avoid);
        if parts.is_empty() {
            return;
        }
        let bb = TopoBuilder::new();
        let mut kept: Vec<TopoShape> = iter_children(self.inner.result());
        let mut fence: HashSet<usize> = kept.iter().map(shape_key).collect();
        for p in &parts {
            let k = shape_key(p);
            if fence.insert(k) && !self.shape_material.contains_key(&k) {
                kept.push(p.clone());
            }
            if material != 0 {
                self.shape_material.entry(k).or_insert(material);
            }
        }
        self.inner.set_result_shape(bb.make_compound_of(&kept).0);
    }

    /// `RemoveFromResult`.
    pub fn remove_from_result(&mut self, take: &[TopoShape], avoid: &[TopoShape]) {
        let drop: HashSet<usize> = self.find_parts(take, avoid).iter().map(shape_key).collect();
        if drop.is_empty() {
            return;
        }
        let kept: Vec<TopoShape> = iter_children(self.inner.result())
            .into_iter()
            .filter(|p| !drop.contains(&shape_key(p)))
            .collect();
        for k in &drop {
            self.shape_material.remove(k);
        }
        self.inner
            .set_result_shape(TopoBuilder::new().make_compound_of(&kept).0);
    }

    /// `FindParts`.
    pub fn find_parts(&self, take: &[TopoShape], avoid: &[TopoShape]) -> Vec<TopoShape> {
        if take.is_empty() {
            return Vec::new();
        }
        let avoid_k: HashSet<usize> = avoid.iter().map(shape_key).collect();
        let take_k: HashSet<usize> = take.iter().map(shape_key).collect();
        let n_take = take_k.len();
        let mut dim_min = 10;
        let mut s_min: Option<&TopoShape> = None;
        for s in take {
            let d = dimension(s);
            if d >= 0 && d < dim_min {
                dim_min = d;
                s_min = Some(s);
            }
        }
        let Some(s_min) = s_min else {
            return Vec::new();
        };
        let Some(ty) = type_to_explore(dim_min) else {
            return Vec::new();
        };
        let mut parts: Vec<TopoShape> = Vec::new();
        let mut fence: HashSet<usize> = HashSet::new();
        for st in explore(s_min, ty) {
            let imgs = match self.inner.history().image(&st) {
                Some(list) if !list.is_empty() => list.to_vec(),
                _ => vec![st],
            };
            for part in imgs {
                let Some(origins) = self.index.get(&shape_key(&part)) else {
                    continue;
                };
                if origins.len() < n_take {
                    continue;
                }
                if origins.iter().any(|o| avoid_k.contains(&shape_key(o))) {
                    continue;
                }
                let found = origins
                    .iter()
                    .filter(|o| take_k.contains(&shape_key(o)))
                    .count();
                if found == n_take && fence.insert(shape_key(&part)) {
                    parts.push(part);
                }
            }
        }
        parts
    }

    fn index_parts(&mut self) {
        self.index.clear();
        let mut fence: HashSet<usize> = HashSet::new();
        let mut all: Vec<TopoShape> = Vec::new();
        let args = self.inner.arguments().to_vec();
        for a in &args {
            let mut flat = Vec::new();
            treat_compound(a, &mut flat, &mut HashSet::new());
            for ss in flat {
                let d = dimension(&ss);
                let Some(ty) = type_to_explore(d) else {
                    continue;
                };
                for st in explore(&ss, ty) {
                    let imgs = match self.inner.history().image(&st) {
                        Some(list) if !list.is_empty() => list.to_vec(),
                        _ => vec![st.clone()],
                    };
                    for part in imgs {
                        let k = shape_key(&part);
                        let list = self.index.entry(k).or_default();
                        if !list.iter().any(|x| x.same_tshape(a)) {
                            list.push(a.clone());
                        }
                        if fence.insert(k) {
                            all.push(part);
                        }
                    }
                }
            }
        }
        self.all_parts = TopoBuilder::new().make_compound_of(&all).0;
    }
}

impl Default for CellsBuilder {
    fn default() -> Self {
        Self::new()
    }
}
