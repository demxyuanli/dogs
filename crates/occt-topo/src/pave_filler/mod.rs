//! `BOPAlgo_PaveFiller` — the Intersection phase of the Boolean Operations
//! algorithm (Phase 19).
//!
//! Source: `BOPAlgo_PaveFiller.hxx/.cxx` + `BOPAlgo_PaveFiller_{1..11}.cxx`
//! (TKBO/BOPAlgo).
//!
//! The class performs the pairwise intersection of the sub-shapes of the
//! arguments in the fixed order:
//!
//! 1. Vertex/Vertex;
//! 2. Vertex/Edge;
//! 3. Edge/Edge;
//! 4. Vertex/Face;
//! 5. Edge/Face;
//! 6. Face/Face.
//!
//! The results of the intersections are stored into the data structure
//! ([`BopdsDS`]) of the algorithm, which later feeds the building phase.
//!
//! This module owns the main class and the pipeline orchestration. The heavy
//! per-pair intersection work lives in the sibling Phase-19 modules
//! (`crate::pave_intersect`, `crate::pave_blocks`, `crate::pave_common`);
//! the pipeline steps delegate to them (and to the pave-common helpers) via
//! [`PaveFillerLike`].
mod prelude {

pub(crate) use std::collections::{HashMap, HashSet};

pub(crate) use crate::abs::ShapeType;
pub(crate) use crate::bopds::BopdsDS;
pub(crate) use crate::int_tools_full::IntToolsContext;
pub(crate) use crate::pave_blocks::PaveFillerLike;
pub(crate) use crate::shape::TopoShape;

}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep_extrema::test_box::unit_box;

    #[test]
    fn defaults() {
        let pf = PaveFiller::new();
        assert!(pf.arguments().is_empty());
        assert_eq!(pf.fuzzy_value(), 1e-7, "fuzzy defaults to Precision::Confusion");
        assert_eq!(pf.glue(), GlueEnum::None, "glue off by default");
        assert!(!pf.non_destructive());
        assert!(!pf.has_errors());
        assert!(!pf.has_warnings());
        assert!(!pf.intersection_done());
        assert_eq!(pf.ds().nb_shapes(), 0);
    }

    #[test]
    fn init_appends_all_subshapes_and_sets_ranges() {
        let mut pf = PaveFiller::new();
        let a = unit_box();
        let c = unit_box();
        pf.set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        pf.init().unwrap();
        // Each unit box contributes 28 shapes (8 V + 12 E + 6 F + 1 Shell + 1 Solid).
        assert_eq!(pf.ds().nb_shapes(), 56);
        assert_eq!(pf.ds().nb_ranges(), 2);
        // A sub-shape of the first argument has rank 0, of the second rank 1.
        let e0 = pf.ds().index(&a.edges[0].0).expect("first edge indexed");
        let e1 = pf.ds().index(&c.edges[0].0).expect("second edge indexed");
        assert_eq!(pf.ds().rank(e0), 0);
        assert_eq!(pf.ds().rank(e1), 1);
        // A sub-shape round-trips through index()/shape().
        let v = pf.ds().index(&a.vertices[0].0).unwrap();
        assert!(pf.ds().shape(v).unwrap().same_tshape(&a.vertices[0].0));
    }

    #[test]
    fn perform_runs_full_pipeline_without_errors() {
        let mut pf = PaveFiller::new();
        let a = unit_box();
        let c = unit_box();
        pf.set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        pf.perform().unwrap();
        assert!(!pf.has_errors(), "errors: {:?}", pf.errors());
        assert!(pf.intersection_done());
        // The intersection created new shapes (SD vertices, crossing vertices,
        // section edges) on top of the 2×28 argument sub-shapes.
        assert!(pf.ds().nb_shapes() >= 56);
    }

    #[test]
    fn errors_and_warnings_accumulate() {
        let mut pf = PaveFiller::new();
        assert!(!pf.has_errors());
        pf.add_error("intersection failed".to_string());
        pf.add_error("builder failed".to_string());
        pf.add_warning("small edges ignored".to_string());
        assert!(pf.has_errors());
        assert!(pf.has_warnings());
        assert_eq!(pf.errors().len(), 2);
        assert_eq!(pf.warnings().len(), 1);
        assert_eq!(pf.errors(), ["intersection failed", "builder failed"]);
        assert_eq!(pf.warnings(), ["small edges ignored"]);
    }

    #[test]
    fn too_few_arguments_fails() {
        let mut pf = PaveFiller::new();
        let err = pf.perform().unwrap_err();
        assert!(err.contains("too few arguments"), "err: {err}");
        assert!(pf.has_errors());
    }

    #[test]
    fn glue_roundtrip() {
        let mut pf = PaveFiller::new();
        assert_eq!(pf.glue(), GlueEnum::None);
        pf.set_glue(GlueEnum::Shift);
        assert_eq!(pf.glue(), GlueEnum::Shift);
        pf.set_glue(GlueEnum::Full);
        assert_eq!(pf.glue(), GlueEnum::Full);
        pf.set_glue(GlueEnum::None);
        assert_eq!(pf.glue(), GlueEnum::None);
    }

    #[test]
    fn fuzzy_value_roundtrip_and_clamp() {
        let mut pf = PaveFiller::new();
        pf.set_fuzzy_value(0.5);
        assert_eq!(pf.fuzzy_value(), 0.5);
        pf.set_fuzzy_value(1e-9);
        assert_eq!(pf.fuzzy_value(), 1e-7, "clamped below to Precision::Confusion");
        pf.set_fuzzy_value(f64::NAN);
        assert_eq!(pf.fuzzy_value(), 1e-7, "NaN degrades to the default");
    }

    #[test]
    fn clear_resets_report_and_ds_keeps_options() {
        let mut pf = PaveFiller::new();
        let a = unit_box();
        pf.set_arguments(&[a.solid.0.clone()]);
        pf.set_glue(GlueEnum::Full);
        pf.set_fuzzy_value(0.25);
        pf.perform().unwrap();
        assert!(pf.intersection_done());
        pf.clear();
        assert!(!pf.intersection_done());
        assert!(!pf.has_errors());
        assert!(!pf.has_warnings());
        assert_eq!(pf.ds().nb_shapes(), 0);
        // User options and arguments survive Clear(), like BOPAlgo_PaveFiller::Clear().
        assert_eq!(pf.glue(), GlueEnum::Full);
        assert_eq!(pf.fuzzy_value(), 0.25);
        assert_eq!(pf.arguments().len(), 1);
    }

    #[test]
    fn repeat_intersection_flag_roundtrip() {
        let mut pf = PaveFiller::new();
        assert!(!pf.repeat_intersection());
        pf.set_repeat_intersection(true);
        assert!(pf.repeat_intersection());
    }
}

mod p01;
pub use p01::*;
