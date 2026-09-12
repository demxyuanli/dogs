    use super::*;
    use crate::brep_extrema::test_box::unit_box;

    fn box_counts(ds: &BopdsDS) -> (usize, usize, usize, usize, usize) {
        let (mut v, mut e, mut f, mut sh, mut so) = (0, 0, 0, 0, 0);
        for si in &ds.shape_infos {
            match si.kind {
                ShapeType::Vertex => v += 1,
                ShapeType::Edge => e += 1,
                ShapeType::Face => f += 1,
                ShapeType::Shell => sh += 1,
                ShapeType::Solid => so += 1,
                _ => {}
            }
        }
        (v, e, f, sh, so)
    }

    #[test]
    fn append_box_counts_all_subshapes() {
        let mut ds = BopdsDS::new();
        let b = unit_box();
        let idx = ds.append(b.solid.0.clone()).unwrap();
        assert_eq!(idx, 0);
        // 8 vertices + 12 edges + 6 faces + 1 shell + 1 solid = 28.
        assert_eq!(ds.nb_shapes(), 28);
        assert_eq!(box_counts(&ds), (8, 12, 6, 1, 1));
    }

    #[test]
    fn append_is_idempotent_and_index_roundtrips() {
        let mut ds = BopdsDS::new();
        let b = unit_box();
        ds.append(b.solid.0.clone()).unwrap();
        let again = ds.append(b.solid.0.clone()).unwrap();
        assert_eq!(again, 0);
        assert_eq!(ds.nb_shapes(), 28);

        // Every distinct sub-shape round-trips through index()/shape().
        for v in &b.vertices {
            let i = ds.index(&v.0).expect("vertex indexed");
            assert!(ds.shape(i).unwrap().same_tshape(&v.0));
        }
        for e in &b.edges {
            let i = ds.index(&e.0).expect("edge indexed");
            assert!(ds.shape(i).unwrap().same_tshape(&e.0));
        }
        for f in &b.faces {
            let i = ds.index(&f.0).expect("face indexed");
            assert!(ds.shape(i).unwrap().same_tshape(&f.0));
        }
        assert_eq!(ds.index(&b.solid.0), Some(0));
        assert_eq!(ds.index(&b.solid.0), ds.index(&b.solid.0));
    }

    #[test]
    fn sub_shape_indices_are_prepared() {
        let mut ds = BopdsDS::new();
        let b = unit_box();
        ds.append(b.solid.0.clone()).unwrap();

        // A face's sub-shapes are its 4 edges + 4 vertices (no wires).
        let face_idx = ds.index(&b.faces[0].0).unwrap();
        let face_info = ds.shape_info(face_idx).unwrap();
        let n_edges = face_info
            .sub_indices
            .iter()
            .filter(|&&i| ds.shape_info(i).unwrap().kind == ShapeType::Edge)
            .count();
        let n_verts = face_info
            .sub_indices
            .iter()
            .filter(|&&i| ds.shape_info(i).unwrap().kind == ShapeType::Vertex)
            .count();
        assert_eq!(n_edges, 4);
        assert_eq!(n_verts, 4);
        // An edge's sub-shapes are its two vertices.
        let edge_idx = ds.index(&b.edges[0].0).unwrap();
        let edge_info = ds.shape_info(edge_idx).unwrap();
        assert_eq!(edge_info.sub_indices.len(), 2);
        assert!(edge_info.sub_indices.iter().all(|&i| ds.shape_info(i).unwrap().kind == ShapeType::Vertex));
    }

    #[test]
    fn ranks_follow_arguments() {
        let mut ds = BopdsDS::new();
        let a = unit_box();
        let c = unit_box();
        ds.init(&[a.solid.0.clone(), c.solid.0.clone()]);
        assert_eq!(ds.nb_shapes(), 56);
        assert_eq!(ds.nb_ranges(), 2);
        assert_eq!(ds.nb_source_shapes(), 56);
        // Shapes of the first argument have rank 0, second rank 1.
        let e0 = ds.index(&a.edges[0].0).unwrap();
        let e1 = ds.index(&c.edges[0].0).unwrap();
        assert_eq!(ds.rank(e0), 0);
        assert_eq!(ds.rank(e1), 1);
        // New shapes (beyond source shapes) are not source shapes.
        assert!(ds.is_new_shape(ds.nb_source_shapes()));
        assert!(!ds.is_new_shape(e0));
    }

    #[test]
    fn pave_blocks_split_on_update() {
        let mut ds = BopdsDS::new();
        let b = unit_box();
        ds.append(b.solid.0.clone()).unwrap();
        let e_idx = ds.index(&b.edges[0].0).unwrap();
        let v0 = ds.index(&b.vertices[0].0).unwrap();
        let v1 = ds.index(&b.vertices[1].0).unwrap();
        let v4 = ds.index(&b.vertices[4].0).unwrap();
        let v5 = ds.index(&b.vertices[5].0).unwrap();

        assert!(!ds.has_pave_blocks(e_idx));
        assert!(ds.pave_blocks(e_idx).is_empty());

        // Add two pave blocks; the second carries an extra pave at t = 0.5.
        {
            let pbs = ds.change_pave_blocks_mut(e_idx);
            let mut pb1 = BopdsPaveBlock::new();
            pb1.edge_index = e_idx;
            pb1.original_edge = e_idx;
            pb1.index1 = v0;
            pb1.index2 = v1;
            pb1.first = 0.0;
            pb1.last = 1.0;

            let mut pb2 = BopdsPaveBlock::new();
            pb2.edge_index = e_idx;
            pb2.original_edge = e_idx;
            pb2.index1 = v4;
            pb2.index2 = v5;
            pb2.first = 0.0;
            pb2.last = 1.0;
            pb2.append_ext_pave(BopdsPave::new(v0, 0.5));

            pbs.push(pb1);
            pbs.push(pb2);
        }

        assert!(ds.has_pave_blocks(e_idx));
        assert_eq!(ds.pave_blocks(e_idx).len(), 2);
        assert!(ds.pave_blocks(e_idx)[1].is_to_update());

        ds.update_pave_blocks();

        // pb1 untouched (1 block) + pb2 split into [0, 0.5] and [0.5, 1] = 3.
        let blocks = ds.pave_blocks(e_idx);
        assert_eq!(blocks.len(), 3);
        assert_eq!(blocks[1].range(), (0.0, 0.5));
        assert_eq!(blocks[2].range(), (0.5, 1.0));
        assert!(!ds.pave_blocks(e_idx)[0].is_to_update());
    }

    #[test]
    fn common_block_updates_and_pool_grows() {
        let mut ds = BopdsDS::new();
        let b = unit_box();
        ds.append(b.solid.0.clone()).unwrap();
        let e_idx = ds.index(&b.edges[0].0).unwrap();

        let mut cb = BopdsCommonBlock::new();
        cb.add_range(0.0, 0.5);
        cb.add_index(e_idx);
        assert!(cb.contains_index(e_idx));
        assert!(cb.contains_range(0.0, 0.5, 1e-9));
        assert!(!cb.contains_range(0.2, 0.8, 1e-9));

        ds.update_common_block(&cb);
        assert_eq!(ds.common_blocks().len(), 1);
        // update_common_block lazily allocated a pave-block slot for the edge.
        assert!(ds.has_pave_blocks(e_idx));
    }

    #[test]
    fn tools_report_shape_types() {
        let b = unit_box();
        assert_eq!(bopds_tools::shape_type(&b.solid.0), ShapeType::Solid);
        assert!(bopds_tools::is_solid(&b.solid.0));
        assert!(bopds_tools::is_vertex(&b.vertices[0].0));
        assert!(bopds_tools::is_edge(&b.edges[0].0));
        assert!(bopds_tools::is_face(&b.faces[0].0));

        let shell = shapes_of(&b.solid.0, ShapeType::Shell);
        assert_eq!(shell.len(), 1);
        assert!(bopds_tools::is_shell(&shell[0]));

        let p = bopds_tools::vertex_point(&b.vertices[0].0);
        assert!(p.is_some());
        assert!((p.unwrap().x() - 0.0).abs() < 1e-12);
        assert!(bopds_tools::vertex_point(&b.solid.0).is_none());

        // Interference-type encoding.
        assert_eq!(bopds_tools::type_to_integer(ShapeType::Vertex), 7);
        assert_eq!(bopds_tools::type_to_integer2(ShapeType::Vertex, ShapeType::Vertex), 0);
        assert_eq!(bopds_tools::type_to_integer2(ShapeType::Edge, ShapeType::Face), 4);
        assert_eq!(bopds_tools::type_to_integer2(ShapeType::Face, ShapeType::Face), 5);
        assert_eq!(bopds_tools::type_to_integer2(ShapeType::Solid, ShapeType::Solid), 9);
    }

    #[test]
    fn pave_and_range_value_semantics() {
        let a = BopdsPave::new(1, 0.5);
        let b = BopdsPave::new(2, 1.0);
        assert!(a.is_less(&b));
        assert!(a < b);
        assert_eq!(a, BopdsPave::new(1, 0.5));
        assert_ne!(a, BopdsPave::new(1, 0.6));
        assert_ne!(a, BopdsPave::new(3, 0.5));

        let r = BopdsIndexRange::new(2, 7);
        assert!(r.contains(5));
        assert!(!r.contains(8));
        assert_eq!(r.indices(), (2, 7));

        let info = BopdsShapeInfo::new(crate::shape::TopoShape::new(ShapeType::Vertex));
        assert_eq!(info.shape_type(), ShapeType::Vertex);
        assert!(info.has_brep());
        // A vertex can participate in V/V, V/E and V/F interferences.
        assert!(info.is_interfering());

        let mut fi = BopdsFaceInfo::new(3);
        fi.add_pave(1, 0.0, 0.5);
        assert_eq!(fi.index(), 3);
        assert_eq!(fi.paves(), &[(1, 0.0, 0.5)]);

        let mut interf = BopdsInterf::new(4, 9);
        assert!(interf.contains(4));
        assert_eq!(interf.opposite_index(4), Some(9));
        assert_eq!(interf.opposite_index(8), None);
        interf.set_index_new(12);
        assert_eq!(interf.get_index_new(), Some(12));
    }

    #[test]
    fn append_info_registers_shape() {
        let mut ds = BopdsDS::new();
        let b = unit_box();
        ds.append(b.solid.0.clone()).unwrap();
        let before = ds.nb_shapes();
        let info = BopdsShapeInfo::new(b.vertices[0].0.clone());
        let i = ds.append_info(info);
        // Already present → returns the existing index.
        assert_eq!(i, ds.index(&b.vertices[0].0).unwrap());
        assert_eq!(ds.nb_shapes(), before);
    }

    #[test]
    fn iterator_finds_inter_argument_candidates() {
        let mut ds = BopdsDS::new();
        let a = unit_box();
        let c = unit_box();
        ds.init(&[a.solid.0.clone(), c.solid.0.clone()]);

        let mut it = BopdsIterator::new();
        it.set_ds(&ds);
        it.prepare();
        it.initialize(ShapeType::Face, ShapeType::Face);
        assert!(it.more(), "overlapping faces of two identical boxes must interfere");
        let (i, j) = it.value();
        // Inter-argument mode never pairs shapes of the same rank.
        assert_ne!(ds.rank(i), ds.rank(j));
        // Each of the 6 faces overlaps the 5 non-parallel faces of the other
        // coincident box → 6 × 5 = 30 face/face candidate pairs.
        assert_eq!(it.expected_length(), 30);
    }

    #[test]
    fn si_iterator_reports_overlapping_cubes() {
        let mut ds = BopdsDS::new();
        let a = unit_box();
        let c = unit_box();
        ds.init(&[a.solid.0.clone(), c.solid.0.clone()]);

        let mut it = BopdsIteratorSI::new();
        it.set_ds(&ds);
        it.prepare();
        it.initialize(ShapeType::Vertex, ShapeType::Vertex);
        assert!(it.more(), "self-intersection VV candidates must be non-empty");
        let (i, j) = it.value();
        assert_ne!(i, j);

        // Face/face candidates are non-empty as well.
        it.initialize(ShapeType::Face, ShapeType::Face);
        assert!(it.more());
        // And solid/solid.
        it.initialize(ShapeType::Solid, ShapeType::Solid);
        assert!(it.more());
        let (s1, s2) = it.value();
        assert_eq!(ds.shape_info(s1).unwrap().kind, ShapeType::Solid);
        assert_eq!(ds.shape_info(s2).unwrap().kind, ShapeType::Solid);
    }

    #[test]
    fn sub_iterator_limited_to_given_subsets() {
        let mut ds = BopdsDS::new();
        let a = unit_box();
        let c = unit_box();
        ds.init(&[a.solid.0.clone(), c.solid.0.clone()]);

        let vertices: Vec<usize> = ds
            .shape_infos
            .iter()
            .enumerate()
            .filter(|(_, si)| si.kind == ShapeType::Vertex)
            .map(|(i, _)| i)
            .collect();
        assert_eq!(vertices.len(), 16);

        let mut it = BopdsSubIterator::new();
        it.set_ds(&ds);
        it.set_subset1(vertices.clone());
        it.set_subset2(vertices.clone());
        it.prepare();
        assert!(it.more());
        let (i, j) = it.value();
        assert!(ds.shape_info(i).unwrap().kind == ShapeType::Vertex);
        assert!(ds.shape_info(j).unwrap().kind == ShapeType::Vertex);

        // Restricting to the vertices of a single box yields no pairs
        // (distinct corners do not overlap).
        let first_box_vertices: Vec<usize> = vertices[..8].to_vec();
        let mut it2 = BopdsSubIterator::new();
        it2.set_ds(&ds);
        it2.set_subset1(first_box_vertices.clone());
        it2.set_subset2(first_box_vertices);
        it2.prepare();
        assert!(!it2.more());
    }

    #[test]
    fn interference_tracking() {
        let mut ds = BopdsDS::new();
        assert!(ds.add_interf(3, 7));
        assert!(!ds.add_interf(7, 3)); // unordered duplicate
        assert!(ds.has_interf(3));
        assert!(ds.has_interf(7));
        assert!(ds.has_interf_pair(7, 3));
        assert!(!ds.has_interf(5));
        assert!(!ds.has_interf_pair(3, 5));

        ds.add_shape_sd(1, 4);
        ds.add_shape_sd(4, 9);
        assert_eq!(ds.get_same_domain_index(1), 9);
        assert_eq!(ds.has_shape_sd(1), Some(9));
        assert_eq!(ds.get_same_domain_index(2), 2);
    }

    #[test]
    fn typed_interferences_and_sd_redirection() {
        let mut ds = BopdsDS::new();
        // A V/V record with a new-vertex index.
        assert!(ds.add_interf_vv(3, 7, Some(5)));
        // Unordered duplicate: no second flat pair, no second typed record.
        assert!(!ds.add_interf_vv(7, 3, Some(5)));
        assert_eq!(ds.interf_vv().len(), 1);
        assert_eq!(ds.interf_vv()[0].get_index_new(), Some(5));
        // The flat table still tracks membership.
        assert!(ds.has_interf_pair(3, 7));
        assert!(ds.has_interf(7));

        // A V/E record whose new vertex is later merged into an SD cluster.
        assert!(ds.add_interf_ve(1, 8, Some(6)));
        ds.add_shape_sd(6, 12);
        ds.add_shape_sd(12, 20);
        ds.update_interfs_with_sd_vertices();
        assert_eq!(
            ds.interf_ve()[0].get_index_new(),
            Some(20),
            "redirected to the final SD representative"
        );
        // The V/V new vertex has no SD partner — left untouched.
        assert_eq!(ds.interf_vv()[0].get_index_new(), Some(5));

        // An interference with no new-vertex index stays untouched.
        assert!(ds.add_interf_ee(2, 9, None));
        ds.update_interfs_with_sd_vertices();
        assert_eq!(ds.interf_ee()[0].get_index_new(), None);
        // PerformNewVertices binds IndexNew; the SD walk then redirects it.
        assert!(ds.bind_ee_new_vertex(2, 9, 6));
        ds.update_interfs_with_sd_vertices();
        assert_eq!(ds.interf_ee()[0].get_index_new(), Some(20));
    }

    #[test]
    fn refine_face_info_in_drops_in_blocks_that_are_on() {
        let mut ds = BopdsDS::new();
        let pool = ds.change_face_info_pool();
        // Face 0: one IN block is also a boundary (ON) block, one is not.
        let mut f0 = BopdsFaceInfo::new(0);
        f0.add_pave_in(3, 0.1, 0.9); // also ON below
        f0.add_pave_in(7, 0.0, 1.0); // not ON
        f0.add_pave_on(3, 0.1, 0.9);
        pool.push(f0);
        // Face 1: empty ON set — nothing to refine.
        let mut f1 = BopdsFaceInfo::new(1);
        f1.add_pave_in(5, 0.0, 1.0);
        pool.push(f1);

        ds.refine_face_info_in();

        let r0 = &ds.face_info_pool()[0];
        assert_eq!(r0.paves_in(), &[(7, 0.0, 1.0)]);
        assert_eq!(r0.paves_on(), &[(3, 0.1, 0.9)]);
        assert_eq!(ds.face_info_pool()[1].paves_in(), &[(5, 0.0, 1.0)]);
    }

    #[test]
    fn release_pave_blocks_clears_untouched_single_block_edges() {
        let mut ds = BopdsDS::new();
        let b = unit_box();
        ds.set_arguments(vec![b.solid.0.clone()]);
        ds.init(&[b.solid.0.clone()]);
        let e = ds.index(&b.edges[0].0).expect("edge indexed");
        ds.init_pave_blocks_for_edge(e);
        assert!(ds.has_pave_blocks(e));

        ds.release_pave_blocks();

        // Untouched (both bounds are source vertices) single block: released.
        assert!(!ds.has_pave_blocks(e), "reference must be dropped");
        assert!(ds.pave_blocks(e).is_empty(), "list contents must be cleared");
    }

    #[test]
    fn release_pave_blocks_keeps_blocks_bounded_by_new_vertices() {
        let mut ds = BopdsDS::new();
        let b = unit_box();
        ds.set_arguments(vec![b.solid.0.clone()]);
        ds.init(&[b.solid.0.clone()]);
        let e = ds.index(&b.edges[0].0).expect("edge indexed");
        ds.init_pave_blocks_for_edge(e);
        // Replace one bound with a new (non-source) vertex index.
        let n_v2 = ds.pave_blocks(e)[0].index2;
        let nb_source = ds.nb_source_shapes();
        ds.change_pave_blocks_mut(e)[0].set_indices(nb_source, n_v2);

        ds.release_pave_blocks();

        assert!(ds.has_pave_blocks(e), "a block with a new bound must be kept");
        assert_eq!(ds.pave_blocks(e).len(), 1);
    }

    #[test]
    fn release_pave_blocks_keeps_common_block_members() {
        let mut ds = BopdsDS::new();
        let b = unit_box();
        ds.set_arguments(vec![b.solid.0.clone()]);
        ds.init(&[b.solid.0.clone()]);
        let e = ds.index(&b.edges[0].0).expect("edge indexed");
        ds.init_pave_blocks_for_edge(e);
        let mut cb = BopdsCommonBlock::new();
        cb.add_index(e);
        cb.add_range(0.0, 1.0);
        ds.update_common_block(&cb);

        ds.release_pave_blocks();

        assert!(ds.has_pave_blocks(e), "a common-block member must be kept");
        assert_eq!(ds.pave_blocks(e).len(), 1);
    }

    #[test]
    fn refine_face_info_on_rebuilds_on_set_and_drops_edge_less_blocks() {
        let mut ds = BopdsDS::new();
        let b = unit_box();
        ds.set_arguments(vec![b.solid.0.clone()]);
        ds.init(&[b.solid.0.clone()]);
        let face = (0..ds.nb_shapes())
            .find(|&i| ds.shape_info(i).map(|s| s.shape_type()) == Some(ShapeType::Face))
            .expect("box has a face");
        ds.change_face_info_pool().push(BopdsFaceInfo::new(face));
        // The face's boundary edges carry a default block each.
        let edges: Vec<usize> = ds
            .shape_info(face)
            .unwrap()
            .sub_shapes()
            .iter()
            .copied()
            .filter(|&s| ds.shape_info(s).map(|si| si.shape_type()) == Some(ShapeType::Edge))
            .collect();
        assert_eq!(edges.len(), 4);
        for &e in &edges {
            ds.init_pave_blocks_for_edge(e);
        }
        // Inject a block without an edge on the first boundary edge.
        let e0 = edges[0];
        {
            let mut pb = BopdsPaveBlock::new();
            pb.edge_index = usize::MAX;
            pb.first = 0.3;
            pb.last = 0.7;
            ds.change_pave_blocks_mut(e0).push(pb);
        }

        ds.refine_face_info_on();

        let on = &ds.face_info_pool()[0].paves_on;
        assert!(
            !on.iter().any(|(edge, _, _)| *edge == usize::MAX),
            "the edge-less block must be dropped, got: {:?}",
            on
        );
        assert_eq!(on.len(), edges.len(), "one real block per boundary edge");
    }

    #[test]
    fn remove_pave_blocks_removes_blocks_from_pool_and_face_info() {
        let mut ds = BopdsDS::new();
        let b = unit_box();
        ds.set_arguments(vec![b.solid.0.clone()]);
        ds.init(&[b.solid.0.clone()]);
        let e = ds.index(&b.edges[0].0).expect("edge indexed");
        ds.init_pave_blocks_for_edge(e);
        // A block whose assigned edge is a newly-appended curve-less edge.
        let n_empty = ds.append(Edge::new().0).unwrap();
        {
            let blocks = ds.change_pave_blocks_mut(e);
            blocks[0].set_edge(n_empty);
        }
        let face = (0..ds.nb_shapes())
            .find(|&i| ds.shape_info(i).map(|s| s.shape_type()) == Some(ShapeType::Face))
            .expect("box has a face");
        let mut fi = BopdsFaceInfo::new(face);
        fi.add_pave(n_empty, 0.0, 1.0);
        fi.add_pave_in(e, 0.0, 1.0);
        ds.change_face_info_pool().push(fi);

        ds.remove_pave_blocks(&HashSet::from([n_empty]));

        // The block re-pointed at n_empty is removed from e's list.
        assert!(ds.pave_blocks(e).iter().all(|pb| pb.edge() != n_empty));
        // The Sc set referencing n_empty is cleaned; the In set stays.
        assert!(ds.face_info_pool()[0].paves().is_empty());
        assert_eq!(ds.face_info_pool()[0].paves_in(), &[(e, 0.0, 1.0)]);
    }

    #[test]
    fn update_face_info_on_collects_boundary_edge_pave_blocks() {
        let mut ds = BopdsDS::new();
        let b = unit_box();
        ds.append(b.solid.0.clone()).unwrap();
        // A box face touches 4 boundary edges; find the first face.
        let face = (0..ds.nb_shapes())
            .find(|&i| ds.shape_info(i).map(|s| s.shape_type()) == Some(ShapeType::Face))
            .expect("box has a face");
        // OCCT `UpdateFaceInfoOn` only acts on faces with a face-info entry.
        ds.change_face_info_pool().push(BopdsFaceInfo::new(face));
        // Split the face's boundary edges so they carry pave blocks.
        let edges: Vec<usize> = ds
            .shape_info(face)
            .unwrap()
            .sub_shapes()
            .iter()
            .copied()
            .filter(|&s| ds.shape_info(s).map(|si| si.shape_type()) == Some(ShapeType::Edge))
            .collect();
        assert_eq!(edges.len(), 4);
        for &e in &edges {
            ds.init_pave_blocks_for_edge(e);
        }

        ds.update_face_info_on(face);

        let expected: Vec<(usize, f64, f64)> = edges
            .iter()
            .flat_map(|&e| {
                ds.pave_blocks(e)
                    .iter()
                    .map(|pb| (pb.edge_index, pb.first, pb.last))
                    .collect::<Vec<_>>()
            })
            .collect();
        assert_eq!(ds.face_info_pool()[0].paves_on(), expected.as_slice());
    }

