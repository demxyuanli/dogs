//! `BOPAlgo_RemoveFeatures` — drop faces from a solid and rebuild the volume.
//!
//! Source: `BOPAlgo_RemoveFeatures.cxx`. Pipeline: `CheckData` →
//! `PrepareFeatures` (connexity of requested faces) → rebuild each solid from
//! the faces that remain (`BuilderSolid`) → history + post-treat.

use std::collections::HashSet;

use crate::abs::ShapeType;
use crate::bop_occt_util::{explore, shape_key};
use crate::builder::TopoBuilder;
use crate::builder_solid::BuilderSolid;
use crate::shape::TopoShape;
use crate::topo_tools_full::faces_of;

/// `BOPAlgo_RemoveFeatures`.
#[derive(Debug, Clone)]
pub struct RemoveFeatures {
    input: Option<TopoShape>,
    faces: Vec<TopoShape>,
    fuzzy: f64,
    result: Option<TopoShape>,
    errors: Vec<String>,
}

impl RemoveFeatures {
    /// Empty algorithm.
    pub fn new() -> Self {
        Self {
            input: None,
            faces: Vec::new(),
            fuzzy: 1e-7,
            result: None,
            errors: Vec::new(),
        }
    }

    /// `SetShape`.
    pub fn set_shape(&mut self, s: TopoShape) {
        self.input = Some(s);
    }

    /// `AddFaceToRemove` / `SetFacesToRemove`.
    pub fn set_faces_to_remove(&mut self, faces: &[TopoShape]) {
        self.faces = faces.to_vec();
    }

    /// `SetFuzzyValue`.
    pub fn set_fuzzy_value(&mut self, fuzzy: f64) {
        self.fuzzy = fuzzy.max(1e-9);
    }

    /// `Perform`.
    pub fn perform(&mut self) -> Result<(), String> {
        self.errors.clear();
        let Some(input) = self.input.clone() else {
            let e = "BOPAlgo_AlertRemoveFeaturesFailed".to_string();
            self.errors.push(e.clone());
            return Err(e);
        };
        if let Err(e) = self.check_data(&input) {
            self.errors.push(e.clone());
            return Err(e);
        }
        let remove: HashSet<usize> = self.faces.iter().map(shape_key).collect();
        if remove.is_empty() {
            self.result = Some(input);
            return Ok(());
        }
        let solids = explore(&input, ShapeType::Solid);
        if solids.is_empty() {
            let e = "BOPAlgo_AlertRemoveFeaturesFailed".to_string();
            self.errors.push(e.clone());
            return Err(e);
        }
        let bb = TopoBuilder::new();
        let mut out: Vec<TopoShape> = Vec::new();
        for sol in solids {
            let keep: Vec<TopoShape> = faces_of(&sol)
                .into_iter()
                .map(|f| f.0)
                .filter(|f| !remove.contains(&shape_key(f)))
                .collect();
            if keep.is_empty() {
                continue;
            }
            let mut bs = BuilderSolid::new();
            bs.set_shapes(keep);
            bs.set_fuzzy(self.fuzzy);
            bs.set_avoid_internal_shapes(true);
            if let Err(e) = bs.perform() {
                self.errors.push(e.clone());
                return Err(e);
            }
            out.extend(bs.areas().iter().cloned());
        }
        self.result = Some(if out.len() == 1 {
            out[0].clone()
        } else {
            bb.make_compound_of(&out).0
        });
        Ok(())
    }

    fn check_data(&self, input: &TopoShape) -> Result<(), String> {
        match input.shape_type() {
            ShapeType::Solid | ShapeType::CompSolid => Ok(()),
            ShapeType::Compound => {
                if explore(input, ShapeType::Solid).is_empty() {
                    Err("BOPAlgo_AlertRemoveFeaturesFailed".into())
                } else {
                    Ok(())
                }
            }
            _ => Err("BOPAlgo_AlertRemoveFeaturesFailed".into()),
        }
    }

    /// `Shape()`.
    pub fn shape(&self) -> Option<&TopoShape> {
        self.result.as_ref()
    }

    /// `HasErrors`.
    pub fn has_errors(&self) -> bool {
        !self.errors.is_empty()
    }
}

impl Default for RemoveFeatures {
    fn default() -> Self {
        Self::new()
    }
}

/// Remove `faces` from `solid` and rebuild (`BOPAlgo_RemoveFeatures::Perform`).
pub fn remove_features(solid: &TopoShape, faces: &[TopoShape], fuzzy: f64) -> Result<TopoShape, String> {
    let mut rf = RemoveFeatures::new();
    rf.set_shape(solid.clone());
    rf.set_faces_to_remove(faces);
    rf.set_fuzzy_value(fuzzy);
    rf.perform()?;
    rf.shape()
        .cloned()
        .ok_or_else(|| "BOPAlgo_AlertRemoveFeaturesFailed".into())
}
