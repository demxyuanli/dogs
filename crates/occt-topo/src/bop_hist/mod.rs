//! History and naming side tables of the Boolean builder.
//!
//! Port of `BRepTools_History` (TKBRep) as consumed by `BOPAlgo_BuilderShape`
//! / `BRepAlgoAPI_BuilderAlgo`, plus the **images** naming table of
//! `BOPAlgo_Builder` (each old shape → the split pieces it expands into).
//!
//! The history keeps the following relations between the input shapes
//! (S1, …, Sm) and the output shapes (T1, …, Tn):
//!
//! 1. an output shape Tj is **generated** from an input shape Si — Tj ∊ G(Si);
//! 2. an output shape Tj is **modified** from an input shape Si — Tj ∊ M(Si);
//! 3. an input shape Si is **removed** — R(Si) = 1 (Si has no output).
//!
//! Only shapes of type vertex, edge, face and solid take part in the
//! relations (`is_supported_type`). The `images` table is a looser naming side
//! table kept by `BOPAlgo_Builder`: old sub-shape → the pieces produced by
//! intersection/splitting. It is queried through `image(old)` /
//! `has_image(old)`.
//!
//! Two sequential histories H12 (S→T) and H23 (T→Q) can be merged into H13
//! (S→Q) with [`BopHistory::merge`], following `BRepTools_History::Merge`:
//!
//! - Tj ∊ G12(Si), Qk ∊ (G23 ∪ M23)(Tj)  ⇒  Qk ∊ G13(Si);
//! - Tj ∊ M12(Si), Qk ∊ G23(Tj)          ⇒  Qk ∊ G13(Si);
//! - Tj ∊ M12(Si), Qk ∊ M23(Tj)          ⇒  Qk ∊ M13(Si).
mod prelude {

pub(crate) use std::collections::{HashMap, HashSet};

pub(crate) use crate::abs::ShapeType;
pub(crate) use crate::shape::TopoShape;
pub(crate) use crate::shape_naming::ShapeId;

}


mod p01;
pub use p01::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
