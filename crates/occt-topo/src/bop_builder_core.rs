//! Shared boolean result types and disjoint / voxel helpers used by
//! [`crate::bop_builder`] and [`crate::bop_builder_planar`].

use crate::abs::ShapeType;
use crate::builder::TopoBuilder;
use crate::shape::{Face, Shell, Solid, TopoShape};
use crate::topo_tools_full::{faces_of, shapes_of};

/// The boolean operation to apply.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoolOp {
    /// A ∪ B.
    Fuse,
    /// A − B.
    Cut,
    /// A ∩ B.
    Common,
}

/// Result of a boolean operation.
#[derive(Debug, Clone)]
pub struct BooleanResult {
    /// The resulting shape (compound for a disjoint Fuse, empty compound for
    /// an empty Common, otherwise the solid or its shell).
    pub shape: TopoShape,
    /// The resulting solid, when the rebuilt shell is closed.
    pub solid: Option<Solid>,
    /// The rebuilt shell(s).
    pub shells: Vec<Shell>,
    /// The selected (rebuilt) faces.
    pub faces: Vec<Face>,
    /// Non-fatal diagnostics (volume inconsistencies, …).
    pub warnings: Vec<String>,
}

// ---------------------------------------------------------------------------
// Result helpers
// ---------------------------------------------------------------------------

pub(crate) fn empty_result(op: BoolOp) -> BooleanResult {
    let b = TopoBuilder::new();
    let comp = b.make_compound_of(&[]);
    BooleanResult {
        shape: comp.0,
        solid: None,
        shells: vec![],
        faces: vec![],
        warnings: if matches!(op, BoolOp::Fuse) { vec!["empty fuse result".into()] } else { vec![] },
    }
}

pub(crate) fn disjoint_result(a: &TopoShape, b: &TopoShape, op: BoolOp) -> BooleanResult {
    let builder = TopoBuilder::new();
    match op {
        BoolOp::Fuse => {
            let comp = builder.make_compound_of(&[a.clone(), b.clone()]);
            let shells: Vec<Shell> = shapes_of(a, ShapeType::Shell)
                .into_iter()
                .map(Shell)
                .chain(shapes_of(b, ShapeType::Shell).into_iter().map(Shell))
                .collect();
            BooleanResult { shape: comp.0, solid: None, shells, faces: vec![], warnings: vec![] }
        }
        BoolOp::Cut => {
            let solid = Solid::wrap(a.clone());
            let shells: Vec<Shell> = shapes_of(a, ShapeType::Shell).into_iter().map(Shell).collect();
            BooleanResult {
                shape: a.clone(),
                solid,
                shells,
                faces: faces_of(a),
                warnings: vec![],
            }
        }
        BoolOp::Common => empty_result(op),
    }
}


pub(crate) fn validate(a: &TopoShape, b: &TopoShape, op: BoolOp, result: &mut BooleanResult) {
    if result.faces.is_empty() {
        return;
    }
    let vol = crate::brep_gprop::volume(&result.shape, 0.02);
    let vol_a = crate::brep_gprop::volume(a, 0.02);
    let vol_b = crate::brep_gprop::volume(b, 0.02);
    match op {
        BoolOp::Fuse => {
            if vol + 1e-6 < vol_a.max(vol_b) {
                result.warnings.push(format!(
                    "fuse volume {vol} below max input {}",
                    vol_a.max(vol_b)
                ));
            }
        }
        BoolOp::Cut => {
            if vol > vol_a + 1e-6 {
                result.warnings.push(format!("cut volume {vol} exceeds input {vol_a}"));
            }
        }
        BoolOp::Common => {
            if vol > vol_a.min(vol_b) + 1e-6 {
                result.warnings.push(format!(
                    "common volume {vol} exceeds min input {}",
                    vol_a.min(vol_b)
                ));
            }
        }
    }
}

/// Wrap a single shape as a [`BooleanResult`] (identity operation).
pub(crate) fn single_shape_result(s: &TopoShape) -> BooleanResult {
    let shells: Vec<Shell> = shapes_of(s, ShapeType::Shell).into_iter().map(Shell).collect();
    let faces = faces_of(s);
    let solid = Solid::wrap(s.clone());
    BooleanResult { shape: s.clone(), solid, shells, faces, warnings: vec![] }
}
