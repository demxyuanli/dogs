//! Replacement context used by `ShapeUpgrade_UnifySameDomain`.
//!
//! `ShapeBuild_ReShape::Merge` is defined in
//! `BRepTools_ReShape.hxx:70-90`: every part is replaced by the product and
//! all parts receive equivalent history. `BopHistory::merge` composes each
//! replacement with the relations already recorded by this context.

use crate::bop_hist::BopHistory;
use crate::shape::TopoShape;
use crate::shape_fix_compose_shell::MapReShape;

#[derive(Clone, Default)]
pub(crate) struct UnifyContext {
    reshape: MapReShape,
    history: BopHistory,
}

impl UnifyContext {
    pub(crate) fn new() -> Self {
        Self::default()
    }

    pub(crate) fn clear(&mut self) {
        self.reshape = MapReShape::new();
        self.history.clear();
    }

    pub(crate) fn apply(&mut self, shape: &TopoShape) -> TopoShape {
        self.reshape.apply(shape)
    }

    pub(crate) fn is_recorded(&self, shape: &TopoShape) -> bool {
        self.reshape.is_recorded(shape)
    }

    pub(crate) fn replace(&mut self, old: &TopoShape, new: &TopoShape) -> Result<(), String> {
        self.reshape.replace(old, new);
        self.compose_replacement(old, new)
    }

    pub(crate) fn merge(&mut self, parts: &[TopoShape], product: &TopoShape) -> Result<(), String> {
        if parts.is_empty() {
            return Ok(());
        }

        let mut step = BopHistory::new();
        for part in parts {
            self.reshape.replace(part, product);
            step.add_modified(part, product.clone())?;
        }
        self.history.merge(&step);
        Ok(())
    }

    pub(crate) fn history(&self) -> &BopHistory {
        &self.history
    }

    fn compose_replacement(&mut self, old: &TopoShape, new: &TopoShape) -> Result<(), String> {
        let mut step = BopHistory::new();
        step.add_modified(old, new.clone())?;
        self.history.merge(&step);
        Ok(())
    }
}
