//! Phase 19 integration — PaveFiller end-to-end over overlapping boxes.
//!
//! Exercises the full intersection pipeline (VV → VE → EE → VF → EF → FF →
//! MakeBlocks → MakePCurves → MakeSplitEdges) on two overlapping unit boxes and
//! on a single closed box (self-interference check). These run the real pipeline
//! with real geometry — not the stub stages.

use occt_topo::bopds::BopdsDS;
use occt_topo::pave_filler::PaveFiller;
use occt_topo::shape::TopoShape;

/// A unit box translated by `(dx, dy, dz)` as a `TopoShape`.
fn unit_box(dx: f64, dy: f64, dz: f64) -> TopoShape {
    let solid = occt_topo::primitives::BRepPrimBox::make_box(1.0, 1.0, 1.0).solid.0;
    let v = occt_core::gp::GpVec::new(dx, dy, dz);
    occt_topo::transform::translated(&solid, &v)
}

/// Count shapes of a given type in the DS via a BFS over shape info.
fn count_kind(ds: &BopdsDS, kind: occt_topo::abs::ShapeType) -> usize {
    (0..ds.nb_shapes())
        .filter(|&i| ds.shape_info(i).map(|si| si.kind == kind).unwrap_or(false))
        .count()
}

#[test]
fn overlapping_boxes_pipeline_runs_without_errors() {
    // Two boxes overlapping along x ∈ [0.25, 1.0] — the intersection region.
    let a = unit_box(0.0, 0.0, 0.0);
    let b = unit_box(0.5, 0.0, 0.0);
    let base = argument_subshapes(&a) + argument_subshapes(&b);
    let mut pf = PaveFiller::new();
    pf.set_arguments(&[a, b]);
    pf.set_fuzzy_value(1e-7);
    let r = pf.perform();
    assert!(r.is_ok(), "perform failed: {r:?}");
    assert!(!pf.has_errors(), "pipeline errors: {:?}", pf.errors());
    // `BOPDS_DS::Init` (`BOPDS_DS.cxx`) fills the DS with every sub-shape of the
    // arguments; the overlapping pair additionally receives the intersection work.
    // OCCT 8.0.0 ground-truth probe (`--ds`, two unit boxes): the disjoint pair
    // reports `NbShapes == NbSourceShapes == 68`, the overlapping pair
    // `NbShapes = 80` (+12). Assert the *relation* rather than a magic count — the
    // port's DS stores a smaller per-box base (56 for two boxes).
    let n = pf.ds().nb_shapes();
    assert!(
        n > base,
        "overlapping boxes must add intersection shapes: {n} vs {base} argument sub-shapes"
    );
}

/// Number of sub-shapes (vertices, edges, faces, shells, solids) of `shape` — the
/// base `BOPDS_DS::Init` (`BOPDS_DS.cxx`) must enter into the DS for every
/// argument.
fn argument_subshapes(shape: &TopoShape) -> usize {
    use occt_topo::abs::ShapeType as T;
    [
        T::Vertex,
        T::Edge,
        T::Face,
        T::Shell,
        T::Solid,
    ]
    .iter()
    .map(|k| occt_topo::topo_tools_full::shapes_of(shape, *k).len())
    .sum()
}

#[test]
fn disjoint_boxes_pipeline_runs() {
    let a = unit_box(0.0, 0.0, 0.0);
    let b = unit_box(3.0, 0.0, 0.0);
    let base = argument_subshapes(&a) + argument_subshapes(&b);
    let mut pf = PaveFiller::new();
    pf.set_arguments(&[a, b]);
    pf.perform().expect("perform");
    assert!(!pf.has_errors(), "{:?}", pf.errors());
    // Disjoint boxes keep exactly the arguments' own sub-shapes (OCCT `--ds`:
    // `NbShapes == NbSourceShapes == 68` for two disjoint unit boxes), so the DS
    // must be at least that base — again a relation, not a magic count.
    let n = pf.ds().nb_shapes();
    assert!(
        n >= base,
        "DS must contain the arguments' sub-shapes: {n} vs {base}"
    );
}

#[test]
fn single_box_has_no_self_interference() {
    let b = unit_box(0.0, 0.0, 0.0);
    let mut pf = PaveFiller::new();
    pf.set_arguments(&[b]);
    pf.perform().expect("perform");
    let si = occt_topo::pave_common::check_self_interference(&mut pf).expect("si check");
    assert!(!si, "closed box must not self-interfere");
}

#[test]
fn overlapping_boxes_create_intersection_vertices_and_edges() {
    let a = unit_box(0.0, 0.0, 0.0);
    let b = unit_box(0.5, 0.0, 0.0);
    let base = argument_subshapes(&a) + argument_subshapes(&b);
    let mut pf = PaveFiller::new();
    pf.set_arguments(&[a, b]);
    pf.perform().expect("perform");

    let nv = count_kind(pf.ds(), occt_topo::abs::ShapeType::Vertex);
    let ne = count_kind(pf.ds(), occt_topo::abs::ShapeType::Edge);
    // Two raw boxes have 8+8 vertices and 12+12 edges. The pipeline must have
    // added intersection vertices/edges, so counts exceed those baselines.
    assert!(nv > 16, "expected intersection vertices, got {nv}");
    assert!(ne > 24, "expected split/intersection edges, got {ne}");
}

#[test]
fn full_box_with_split_keeps_vertices_bounded() {
    // A single box run through the pipeline keeps all vertices inside the box
    // bounding volume (no runaway coordinates).
    let b = unit_box(0.0, 0.0, 0.0);
    let mut pf = PaveFiller::new();
    pf.set_arguments(&[b]);
    pf.perform().expect("perform");
    let ds = pf.ds();
    for i in 0..ds.nb_shapes() {
        if let Some(si) = ds.shape_info(i) {
            if si.kind == occt_topo::abs::ShapeType::Vertex {
                let v = occt_topo::shape::Vertex::wrap(si.shape.clone()).expect("vertex");
                let p = occt_topo::topo_tools_full::vertex_position(&v);
                for c in [p.x(), p.y(), p.z()] {
                    assert!(c >= -0.1 && c <= 1.1, "vertex out of bounds: {p:?}");
                }
            }
        }
    }
}
