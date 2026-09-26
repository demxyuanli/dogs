
use super::*;
    use crate::brep_extrema::test_box::unit_box;
    use crate::shape::TopoShape;

    /// Returns a face of the unit box (a plain TopoShape) and a distinct edge.
    fn fixtures() -> (TopoShape, TopoShape, TopoShape) {
        let b = unit_box();
        // The box faces/edges/vertices are shared TShapes, so ShapeId::of is
        // stable across clones.
        let face = b.faces[0].0.clone();
        let edge = b.edges[0].0.clone();
        let vertex = b.vertices[0].0.clone();
        (face, edge, vertex)
    }

    #[test]
    fn add_and_query_images() {
        let (face, edge, vertex) = fixtures();
        let mut h = BopHistory::new();
        assert!(!h.has_image(&face));
        h.add_image(&face, edge.clone());
        h.add_image(&face, vertex.clone());
        assert!(h.has_image(&face));
        let img = h.image(&face).expect("images recorded");
        assert_eq!(img.len(), 2);
        assert!(img[0].same_tshape(&edge));
        assert!(img[1].same_tshape(&vertex));
        assert!(h.image(&edge).is_none(), "edge has no images");
        assert_eq!(h.images().len(), 1);
    }

    #[test]
    fn add_and_query_modified_generated() {
        let (face, edge, vertex) = fixtures();
        let mut h = BopHistory::new();
        h.add_modified(&face, edge.clone()).unwrap();
        h.add_modified(&face, vertex.clone()).unwrap();
        h.add_generated(&edge, vertex.clone()).unwrap();

        let m = h.modified(&face).expect("modified recorded");
        assert_eq!(m.len(), 2);
        assert!(m[0].same_tshape(&edge));
        assert!(m[1].same_tshape(&vertex));

        let g = h.generated(&edge).expect("generated recorded");
        assert_eq!(g.len(), 1);
        assert!(g[0].same_tshape(&vertex));

        assert!(h.has_modified(&face));
        assert!(h.has_generated(&edge));
        assert!(!h.has_generated(&face));
        assert!(!h.modified(&vertex).is_some());
        assert!(h.has_any_modified());
        assert!(h.has_any_generated());
    }

    #[test]
    fn remove_and_is_deleted() {
        let (face, edge, _) = fixtures();
        let mut h = BopHistory::new();
        h.add_generated(&face, edge.clone()).unwrap();
        h.add_removed(edge.clone()).unwrap();
        assert!(h.is_deleted(&edge));
        assert!(h.has_any_deleted());
        assert_eq!(h.removed().len(), 1);
        assert!(h.removed()[0].same_tshape(&edge));
        // A removed shape may still have generations, but not modifications:
        // here `face` still generates `edge`.
        assert!(h.generated(&face).is_some());
        // Repeated removal is a no-op.
        h.add_removed(edge.clone()).unwrap();
        assert_eq!(h.removed().len(), 1);
        assert!(!h.is_deleted(&face));
    }

    #[test]
    fn remove_clears_modified_but_keeps_generated() {
        let (face, edge, vertex) = fixtures();
        let mut h = BopHistory::new();
        h.add_modified(&face, edge.clone()).unwrap();
        h.add_generated(&face, vertex.clone()).unwrap();
        h.add_removed(edge.clone()).unwrap();
        assert!(h.is_deleted(&edge));
        // Modified entry for `edge` was unbound by Remove(); generated stays.
        assert!(!h.has_modified(&edge));
        assert!(h.has_generated(&face));
    }

    #[test]
    fn generated_and_modified_are_mutually_exclusive() {
        let (face, edge, vertex) = fixtures();
        let mut h = BopHistory::new();
        h.add_generated(&face, edge.clone()).unwrap();
        // Adding the same shape as a modification moves it out of generated.
        h.add_modified(&face, edge.clone()).unwrap();
        assert!(h.has_modified(&face));
        assert!(!h.has_generated(&face));
        assert_eq!(h.generated(&face).map(|v| v.len()), None);

        h.add_modified(&face, vertex.clone()).unwrap();
        h.add_generated(&face, vertex.clone()).unwrap();
        assert!(h.has_generated(&face));
        let m = h.modified(&face).unwrap();
        assert_eq!(m.len(), 1);
        assert!(m[0].same_tshape(&edge));
    }

    #[test]
    fn unsupported_types_are_rejected() {
        let b = unit_box();
        // The compound root of the box is not vertex/edge/face/solid.
        let root = TopoShape::new(crate::abs::ShapeType::Compound);
        let mut h = BopHistory::new();
        assert!(h.add_modified(&root, b.edges[0].0.clone()).is_err());
        assert!(h.add_generated(&root, b.vertices[0].0.clone()).is_err());
        assert!(h.add_removed(root).is_err());
        // But images accept any type.
        let r = TopoShape::new(crate::abs::ShapeType::Compound);
        h.add_image(&r, b.faces[0].0.clone());
        assert!(h.has_image(&r));
    }

    #[test]
    fn clear_resets_everything() {
        let (face, edge, _) = fixtures();
        let mut h = BopHistory::new();
        h.add_image(&face, edge.clone());
        h.add_modified(&face, edge.clone()).unwrap();
        h.add_generated(&edge, face.clone()).unwrap();
        h.add_removed(edge.clone()).unwrap();
        assert!(!h.is_empty());
        h.clear();
        assert!(h.is_empty());
        assert!(!h.has_image(&face));
        assert!(!h.has_modified(&face));
        assert!(!h.has_generated(&edge));
        assert!(!h.is_deleted(&edge));
    }

    #[test]
    fn merge_composes_generated() {
        // H1: S --generated--> A ; H2: A --generated--> G.
        // Merge: S --generated--> [A, G].
        let (s, a, g) = fixtures();
        let mut h1 = BopHistory::new();
        h1.add_generated(&s, a.clone()).unwrap();
        let mut h2 = BopHistory::new();
        h2.add_generated(&a, g.clone()).unwrap();

        h1.merge(&h2);
        let gen = h1.generated(&s).expect("S generates after merge");
        assert_eq!(gen.len(), 2, "kept A and gained G");
        assert!(gen[0].same_tshape(&a));
        assert!(gen[1].same_tshape(&g));
        // A was reached from S's list, so its relation was folded into S
        // (Phase 2 of Merge exposes only *unreached* intermediates).
        assert!(!h1.has_generated(&a), "A's generation folded into S");
    }

    #[test]
    fn merge_replaces_modified() {
        // H1: S --modified--> [A, B] ; H2: A --modified--> A'.
        // Merge: S --modified--> [A', B], and A exposed with its modification.
        let (s, a, b) = fixtures();
        let a_prime = unit_box().faces[1].0.clone();
        let mut h1 = BopHistory::new();
        h1.add_modified(&s, a.clone()).unwrap();
        h1.add_modified(&s, b.clone()).unwrap();
        let mut h2 = BopHistory::new();
        h2.add_modified(&a, a_prime.clone()).unwrap();

        h1.merge(&h2);
        // Surviving members keep their order, then the replacement additions
        // are appended (OCCT: `add(aL12, aAdditions[aI])` after the loop).
        let m = h1.modified(&s).expect("S modified after merge");
        assert_eq!(m.len(), 2, "A replaced by A', B kept");
        assert!(m[0].same_tshape(&b), "survivor B kept first");
        assert!(m[1].same_tshape(&a_prime), "replacement A' appended");
        // A was reached from S's list, so it is not re-exposed as a key.
        assert!(!h1.has_modified(&a));
    }

    #[test]
    fn merge_generations_of_modified_become_generations() {
        // H1: S --modified--> A ; H2: A --generated--> G.
        // Merge: S --modified--> A (kept), S --generated--> G (rule 2).
        let (s, a, g) = fixtures();
        let mut h1 = BopHistory::new();
        h1.add_modified(&s, a.clone()).unwrap();
        let mut h2 = BopHistory::new();
        h2.add_generated(&a, g.clone()).unwrap();

        h1.merge(&h2);
        let m = h1.modified(&s).unwrap();
        assert_eq!(m.len(), 1);
        assert!(m[0].same_tshape(&a));
        let gen = h1.generated(&s).expect("S gained a generation");
        assert_eq!(gen.len(), 1);
        assert!(gen[0].same_tshape(&g));
    }

    #[test]
    fn merge_removal_of_only_split_removes_initial() {
        // H1: S --modified--> A ; H2: A removed.
        // Merge: S's only modification disappears -> S itself is removed.
        let (s, a, _) = fixtures();
        let mut h1 = BopHistory::new();
        h1.add_modified(&s, a.clone()).unwrap();
        let mut h2 = BopHistory::new();
        h2.add_removed(a.clone()).unwrap();

        h1.merge(&h2);
        assert!(h1.is_deleted(&s), "S has no surviving splits -> removed");
        assert!(!h1.has_modified(&s));
    }

    #[test]
    fn merge_partial_removal_keeps_survivors() {
        // H1: S --modified--> [A, B] ; H2: A removed.
        // Merge: S --modified--> [B]; S is not removed because B survives.
        let (s, a, b) = fixtures();
        let mut h1 = BopHistory::new();
        h1.add_modified(&s, a.clone()).unwrap();
        h1.add_modified(&s, b.clone()).unwrap();
        let mut h2 = BopHistory::new();
        h2.add_removed(a.clone()).unwrap();

        h1.merge(&h2);
        assert!(!h1.is_deleted(&s));
        let m = h1.modified(&s).unwrap();
        assert_eq!(m.len(), 1);
        assert!(m[0].same_tshape(&b));
    }

    #[test]
    fn merge_exposes_unreached_intermediates_and_clears_removed() {
        // H1 has no relation for B except a removal; H2: B --modified--> B'.
        // Merge: B is exposed with its modification, and B is no longer removed.
        let (_, b, b_prime) = fixtures();
        let mut h1 = BopHistory::new();
        h1.add_removed(b.clone()).unwrap();
        assert!(h1.is_deleted(&b));
        let mut h2 = BopHistory::new();
        h2.add_modified(&b, b_prime.clone()).unwrap();

        h1.merge(&h2);
        assert!(!h1.is_deleted(&b), "re-exposed intermediate is un-removed");
        assert!(h1.has_modified(&b));
        let m = h1.modified(&b).unwrap();
        assert_eq!(m.len(), 1);
        assert!(m[0].same_tshape(&b_prime));
    }

    #[test]
    fn merge_carries_forward_removed() {
        // H1 empty relations; H2: X removed. Merge: X removed.
        let (x, _, _) = fixtures();
        let mut h1 = BopHistory::new();
        let mut h2 = BopHistory::new();
        h2.add_removed(x.clone()).unwrap();
        h1.merge(&h2);
        assert!(h1.is_deleted(&x));
        assert_eq!(h1.removed().len(), 1);
    }

    #[test]
    fn merge_composes_images() {
        // H1.images: S -> [A, B]; H2.images: A -> [A'].
        // Merge: S -> [A', B].
        let (s, a, b) = fixtures();
        let a_prime = unit_box().faces[1].0.clone();
        let mut h1 = BopHistory::new();
        h1.add_image(&s, a.clone());
        h1.add_image(&s, b.clone());
        let mut h2 = BopHistory::new();
        h2.add_image(&a, a_prime.clone());

        h1.merge(&h2);
        let img = h1.image(&s).expect("S images after merge");
        assert_eq!(img.len(), 2);
        assert!(img[0].same_tshape(&a_prime));
        assert!(img[1].same_tshape(&b));
        assert!(h1.has_image(&a), "A-only image key bound as-is");
    }

    #[test]
    fn merge_images_drop_removed_pieces() {
        // H1.images: S -> [A, B]; H2: A removed.
        // Merge: S -> [B] (A dropped).
        let (s, a, b) = fixtures();
        let mut h1 = BopHistory::new();
        h1.add_image(&s, a.clone());
        h1.add_image(&s, b.clone());
        let mut h2 = BopHistory::new();
        h2.add_removed(a.clone()).unwrap();

        h1.merge(&h2);
        let img = h1.image(&s).expect("S images after merge");
        assert_eq!(img.len(), 1);
        assert!(img[0].same_tshape(&b));
    }

    #[test]
    fn merge_into_empty_is_identity() {
        let (s, a, _) = fixtures();
        let mut h1 = BopHistory::new();
        let mut h2 = BopHistory::new();
        h2.add_modified(&s, a.clone()).unwrap();
        h2.add_generated(&a, s.clone()).unwrap();
        h2.add_removed(a.clone()).unwrap();
        h1.merge(&h2);
        assert!(h1.has_modified(&s));
        assert!(h1.has_generated(&a));
        assert_eq!(h1.modified(&s).unwrap().len(), 1);
        assert_eq!(h1.generated(&a).unwrap().len(), 1);
        // `a` is removed in H23 but re-exposed as a key carrying a generation,
        // so its removed marker is cleared (Merge phase 2).
        assert!(!h1.is_deleted(&a), "removed marker cleared for re-exposed key");
    }

    #[test]
    fn dump_reports_counts() {
        let (s, a, _) = fixtures();
        let mut h = BopHistory::new();
        h.add_modified(&s, a.clone()).unwrap();
        h.add_removed(a).unwrap();
        let d = h.dump();
        assert!(d.contains("1 Deleted shapes"));
        assert!(d.contains("1 Modified shapes"));
        assert!(d.contains("0 Generated shapes"));
    }
