use super::prelude::*;
use super::*;
    use crate::abs::ShapeType;

    /// Test helper: a validated range.
    fn r(a: f64, b: f64) -> IntRange {
        IntRange::new(a, b).unwrap()
    }

    fn shape(t: ShapeType) -> TopoShape {
        TopoShape::new(t)
    }

    // ---- IntRange ----

    #[test]
    fn intrange_contains_inclusive_boundaries() {
        let range = r(0.0, 1.0);
        assert!(range.contains(0.0), "lower boundary");
        assert!(range.contains(1.0), "upper boundary");
        assert!(range.contains(0.5));
        assert!(!range.contains(-0.1));
        assert!(!range.contains(1.1));
        assert_eq!(range.length(), 1.0);
    }

    #[test]
    fn intrange_overlaps_touching_and_disjoint() {
        assert!(r(0.0, 1.0).overlaps(&r(1.0, 2.0)), "touching ranges overlap");
        assert!(r(0.0, 2.0).overlaps(&r(1.0, 3.0)));
        assert!(r(2.0, 3.0).overlaps(&r(0.0, 2.5)));
        assert!(!r(0.0, 1.0).overlaps(&r(1.5, 2.0)), "gap between");
        assert!(r(0.0, 1.0).disjoint(&r(1.5, 2.0)));
        assert!(!r(0.0, 1.0).disjoint(&r(1.0, 2.0)));
    }

    #[test]
    fn intrange_merge_spans() {
        assert_eq!(r(1.0, 3.0).merge(&r(0.0, 2.0)), r(0.0, 3.0));
        // Disjoint ranges merge to the bounding hull.
        assert_eq!(r(0.0, 1.0).merge(&r(2.0, 3.0)), r(0.0, 3.0));
    }

    #[test]
    fn intrange_new_rejects_inverted_and_nan() {
        assert!(IntRange::new(2.0, 1.0).is_err(), "inverted bounds rejected");
        assert!(IntRange::new(f64::NAN, 1.0).is_err());
        assert!(IntRange::new(0.0, f64::INFINITY).is_err());
        assert!(IntRange::new(0.0, 0.0).is_ok(), "degenerate point range allowed");
        assert!(IntRange::new(-1.0, 1.0).is_ok());
        assert!(r(0.0, 1.0).is_valid());
    }

    #[test]
    fn inttools_range_alias() {
        let a: IntToolsRange = r(0.0, 2.0);
        let b: IntRange = a;
        assert_eq!(b.first(), 0.0);
        assert_eq!(b.last(), 2.0);
    }

    // ---- CommonPrt ----

    #[test]
    fn commonprt_default_is_unknown_empty() {
        let cp = CommonPrt::new();
        assert_eq!(cp.part_type(), CommonPartType::Unknown);
        assert!(cp.is_empty());
        assert!(!cp.is_edge_part());
        assert!(!cp.is_face_part());
    }

    #[test]
    fn commonprt_construct_and_accessors() {
        let face = shape(ShapeType::Face);
        let v1 = shape(ShapeType::Vertex);
        let v2 = shape(ShapeType::Vertex);
        let mut cp = CommonPrt::with(CommonPartType::Edge, r(0.5, 2.5), Some(face), vec![v1, v2]);
        assert!(cp.is_edge_part());
        assert_eq!(cp.range(), r(0.5, 2.5));
        assert!(cp.face().is_some());
        assert_eq!(cp.vertices().len(), 2);

        cp.set_part_type(CommonPartType::Face);
        assert!(cp.is_face_part());
        cp.set_range(r(1.0, 1.0));
        assert_eq!(cp.range(), r(1.0, 1.0));
        cp.set_face(None);
        assert!(cp.face().is_none());
        cp.add_vertex(shape(ShapeType::Vertex));
        assert_eq!(cp.vertices().len(), 3);
    }

    // ---- IntRoot ----

    #[test]
    fn introot_construct_and_accessors() {
        let root = IntRoot::new(3, RootType::IsRoot, r(0.9, 1.1));
        assert_eq!(root.root_index(), 3);
        assert_eq!(root.root_type(), RootType::IsRoot);
        assert!(root.is_root());
        assert!(!root.is_tangent());
        assert!(!root.is_conflict());
        assert!((root.root_value() - 1.0).abs() < 1e-12);

        let tan = IntRoot::with_conflict(1, RootType::IsTangent, r(0.0, 1.0), true);
        assert!(tan.is_tangent());
        assert!(tan.is_conflict());
        assert_eq!(tan.root_index(), 1);
    }

    #[test]
    fn introot_default_unknown() {
        let root = IntRoot::default();
        assert_eq!(root.root_type(), RootType::IsUnknown);
        assert_eq!(root.root_index(), 0);
        assert_eq!(root.range(), r(0.0, 0.0));
        root_ok(root);
    }

    fn root_ok(root: IntRoot) {
        assert!(!root.is_root());
    }

    // ---- PntOnFace / PntOn2Faces ----

    #[test]
    fn pntonface_construct() {
        let pnt = GpPnt::new(1.0, 2.0, 3.0);
        let p = PntOnFace::new(7, (0.25, 0.75), pnt);
        assert_eq!(p.face_index(), 7);
        assert_eq!(p.uv(), (0.25, 0.75));
        assert_eq!(p.u(), 0.25);
        assert_eq!(p.v(), 0.75);
        assert_eq!(*p.pnt(), pnt);
        assert_eq!(PntOnFace::default().face_index(), 0);
    }

    #[test]
    fn pnton2faces_construct() {
        let p1 = GpPnt::new(1.0, 0.0, 0.0);
        let p2 = GpPnt::new(0.0, 1.0, 0.0);
        let hit = PntOn2Faces::new(0, 1, p1, p2, (0.0, 1.0), (1.0, 0.0));
        assert_eq!(hit.face1_index(), 0);
        assert_eq!(hit.face2_index(), 1);
        assert_eq!(*hit.pnt1(), p1);
        assert_eq!(*hit.pnt2(), p2);
        assert_eq!(hit.uv1(), (0.0, 1.0));
        assert_eq!(hit.uv2(), (1.0, 0.0));
    }

    // ---- IntCurve ----

    #[test]
    fn intcurve_construct_and_kinds() {
        let f1 = shape(ShapeType::Face);
        let f2 = shape(ShapeType::Face);
        let curve = IntCurve::new(CurveKind::Circle, Some(f1), Some(f2), r(-1.0, 1.0));
        assert_eq!(curve.kind(), CurveKind::Circle);
        assert!(curve.face1().is_some());
        assert!(curve.face2().is_some());
        assert_eq!(curve.range(), r(-1.0, 1.0));

        let line = IntCurve::new(CurveKind::Line, None, None, r(0.0, 2.0));
        assert_eq!(line.kind(), CurveKind::Line);
        assert!(line.face1().is_none());

        let mut b = IntCurve::default();
        assert_eq!(b.kind(), CurveKind::Other);
        b.set_range(r(0.0, 5.0));
        assert_eq!(b.range(), r(0.0, 5.0));
        // All curve kinds are representable.
        for k in [CurveKind::Ellipse, CurveKind::Parabola, CurveKind::Hyperbola, CurveKind::BSpline] {
            assert_eq!(IntCurve::new(k, None, None, r(0.0, 1.0)).kind(), k);
        }
    }

    // ---- MarkedRangeSet ----

    #[test]
    fn marked_insert_disjoint_sorted() {
        let mut s = MarkedRangeSet::new();
        s.insert(r(4.0, 6.0));
        s.insert(r(0.0, 2.0));
        assert_eq!(s.sorted_ranges(), vec![r(0.0, 2.0), r(4.0, 6.0)]);
        assert_eq!(s.len(), 2);
        assert!(s.is_marked(1.0) == false);
    }

    #[test]
    fn marked_insert_overlap_merges() {
        let mut s = MarkedRangeSet::new();
        s.insert(r(0.0, 5.0));
        s.insert(r(2.0, 3.0)); // interior, same (unmarked) flag -> merged
        assert_eq!(s.sorted_ranges(), vec![r(0.0, 5.0)]);
        assert_eq!(s.len(), 1);
    }

    #[test]
    fn marked_mark_whole_interval() {
        let mut s = MarkedRangeSet::new();
        s.insert(r(0.0, 1.0));
        s.insert(r(3.0, 4.0));
        s.mark(0.5).unwrap();
        assert!(s.is_marked(0.5));
        assert!(!s.is_marked(3.5), "other interval stays unmarked");
        assert_eq!(s.marked_ranges(), vec![r(0.0, 1.0)]);
    }

    #[test]
    fn marked_is_marked_uncovered_false() {
        let mut s = MarkedRangeSet::new();
        s.insert(r(0.0, 1.0));
        assert!(!s.is_marked(10.0));
        assert!(!s.is_marked(-1.0));
    }

    #[test]
    fn marked_mark_outside_errors() {
        let mut s = MarkedRangeSet::new();
        s.insert(r(0.0, 1.0));
        assert!(s.mark(2.0).is_err());
        assert!(s.mark(f64::NAN).is_err());
    }

    #[test]
    fn marked_insert_inside_marked_unmarks_overlap() {
        let mut s = MarkedRangeSet::new();
        s.insert(r(0.0, 5.0));
        s.mark(2.0).unwrap();
        assert!(s.is_marked(4.0));
        // Inserting an unmarked interior range splits and unmarks the overlap.
        s.insert(r(1.0, 3.0));
        assert!(s.is_marked(0.5), "left fragment keeps mark");
        assert!(!s.is_marked(2.0), "overlap becomes unmarked");
        assert!(s.is_marked(4.0), "right fragment keeps mark");
        assert_eq!(s.sorted_ranges(), vec![r(0.0, 1.0), r(1.0, 3.0), r(3.0, 5.0)]);
        assert_eq!(s.marked_ranges(), vec![r(0.0, 1.0), r(3.0, 5.0)]);
    }

    #[test]
    fn marked_unite() {
        let mut a = MarkedRangeSet::new();
        a.insert(r(0.0, 2.0));
        a.insert(r(4.0, 6.0));
        a.mark(1.0).unwrap();
        let mut b = MarkedRangeSet::new();
        b.insert(r(1.0, 5.0));
        b.mark(3.0).unwrap();

        let u = a.unite(&b);
        assert_eq!(u.sorted_ranges(), vec![r(0.0, 5.0), r(5.0, 6.0)]);
        assert!(u.is_marked(0.5), "covered only by marked A");
        assert!(u.is_marked(2.5), "covered only by marked B");
        assert!(!u.is_marked(5.5), "covered by unmarked A only");
    }

    #[test]
    fn marked_intersect() {
        let mut a = MarkedRangeSet::new();
        a.insert(r(0.0, 2.0));
        a.insert(r(4.0, 6.0));
        a.mark(1.0).unwrap();
        let mut b = MarkedRangeSet::new();
        b.insert(r(1.0, 5.0));
        b.mark(3.0).unwrap();

        let i = a.intersect(&b);
        assert_eq!(i.sorted_ranges(), vec![r(1.0, 2.0), r(4.0, 5.0)]);
        assert!(i.is_marked(1.5), "marked in both");
        assert!(!i.is_marked(4.5), "marked only in B -> not marked in intersection");
    }

    #[test]
    fn marked_subtract() {
        let mut a = MarkedRangeSet::new();
        a.insert(r(0.0, 2.0));
        a.insert(r(4.0, 6.0));
        a.mark(1.0).unwrap();
        let mut b = MarkedRangeSet::new();
        b.insert(r(1.0, 5.0));
        b.mark(3.0).unwrap();

        let d = a.subtract(&b);
        assert_eq!(d.sorted_ranges(), vec![r(0.0, 1.0), r(5.0, 6.0)]);
        assert!(d.is_marked(0.5), "surviving marked part of A");
        assert!(!d.is_marked(5.5), "surviving unmarked part of A");
    }

    #[test]
    fn marked_unite_with_empty() {
        let mut a = MarkedRangeSet::new();
        a.insert(r(0.0, 2.0));
        a.mark(1.0).unwrap();
        let empty = MarkedRangeSet::new();
        let u = a.unite(&empty);
        assert_eq!(u.sorted_ranges(), a.sorted_ranges());
        assert!(u.is_marked(1.0));
        assert!(empty.is_empty());
        assert!(a.subtract(&a).is_empty(), "self-subtraction clears everything");
    }

    // ---- LocalizeData ----

    #[test]
    fn localizedata_construct_and_uv() {
        let mut ld = LocalizeData::new(r(0.0, 3.0));
        assert_eq!(ld.face_range(), r(0.0, 3.0));
        ld.add_edge_range(r(1.0, 2.0));
        ld.add_uv_point(0.5, 1.5);
        assert_eq!(ld.edge_ranges(), &[r(1.0, 2.0)]);
        assert_eq!(ld.uv_points(), &[(0.5, 1.5)]);
        assert_eq!(LocalizeData::default().face_range(), r(0.0, 0.0));
    }

    #[test]
    fn localizedata_sort_edge_ranges_merges() {
        let mut ld = LocalizeData::new(r(0.0, 10.0));
        ld.add_edge_range(r(6.0, 8.0));
        ld.add_edge_range(r(0.0, 2.0));
        ld.add_edge_range(r(2.0, 4.0)); // adjacent to [0,2] -> merged
        ld.add_edge_range(r(9.0, 10.0));
        assert!(!ld.is_sorted());
        ld.sort_edge_ranges();
        assert_eq!(ld.edge_ranges(), &[r(0.0, 4.0), r(6.0, 8.0), r(9.0, 10.0)]);
        assert!(ld.is_sorted());
    }

    #[test]
    fn localizedata2_construct() {
        let mut ld2 = LocalizeData2::new(r(0.0, 1.0), r(2.0, 3.0));
        ld2.add_edge_range(r(0.5, 1.5));
        ld2.add_uv1(0.1, 0.2);
        ld2.add_uv2(0.3, 0.4);
        assert_eq!(ld2.face1_range(), r(0.0, 1.0));
        assert_eq!(ld2.face2_range(), r(2.0, 3.0));
        assert_eq!(ld2.edge_ranges(), &[r(0.5, 1.5)]);
        assert_eq!(ld2.uv1_points(), &[(0.1, 0.2)]);
        assert_eq!(ld2.uv2_points(), &[(0.3, 0.4)]);
        ld2.sort_edge_ranges();
        assert_eq!(ld2.edge_ranges(), &[r(0.5, 1.5)]);
    }

    #[test]
    fn merge_plain_dedup_and_touch() {
        let v = vec![r(2.0, 3.0), r(0.0, 1.0), r(1.0, 2.0), r(2.0, 3.0)];
        assert_eq!(merge_plain(v), vec![r(0.0, 3.0)]);
    }
