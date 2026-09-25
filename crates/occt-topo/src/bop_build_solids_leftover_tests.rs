use super::*;
use super::super::*;
use std::sync::Arc;

use occt_core::gp::{GpAx3, GpDir, GpLin, GpPln, GpPnt, GpVec};
use occt_geom::{GeomLine, GeomPlane, Surface};

use crate::bop_hist::BopHistory;
use crate::bopds::BopdsDS;
use crate::primitives::BRepPrimBox;
use crate::shape::{Edge, Face, Solid, Vertex};
use crate::topo_tools_full::{edges_of, vertices_of};

    /// A minimal `BopBuildOps` host for the isolated tests.
    struct StubBuilder {
        ds: BopdsDS,
        history: BopHistory,
        fuzzy: f64,
        args: Vec<TopoShape>,
        origins: HashMap<usize, Vec<TopoShape>>,
    }

    impl BopBuildOps for StubBuilder {
        fn ds(&self) -> &BopdsDS {
            &self.ds
        }
        fn history(&self) -> &BopHistory {
            &self.history
        }
        fn history_mut(&mut self) -> &mut BopHistory {
            &mut self.history
        }
        fn fuzzy_value(&self) -> f64 {
            self.fuzzy
        }
        fn arguments(&self) -> &[TopoShape] {
            &self.args
        }
        fn origins_mut(&mut self) -> &mut HashMap<usize, Vec<TopoShape>> {
            &mut self.origins
        }
    }

    fn stub(ds: BopdsDS, history: BopHistory, args: Vec<TopoShape>) -> StubBuilder {
        StubBuilder { ds, history, fuzzy: 1e-7, args, origins: HashMap::new() }
    }

    /// The face of `faces` whose boundary-vertex mean lies at height `z`.
    fn face_at_z(faces: &[Face], z: f64) -> Face {
        faces
            .iter()
            .find(|f| {
                let vs = vertices_of(&f.0);
                if vs.is_empty() {
                    return false;
                }
                let zavg =
                    vs.iter().map(|v| BRepTool::vertex_point(v).z()).sum::<f64>() / vs.len() as f64;
                (zavg - z).abs() < 1e-9
            })
            .cloned()
            .expect("face at height")
    }

    /// A fresh face on the same surface and boundary edges as `src` — a new
    /// `TShape`, so it is a distinct split image.
    fn re_face(src: &Face) -> Face {
        let b = TopoBuilder::new();
        let edges = edges_of(&src.0);
        let wire = b.make_wire(&edges);
        let surf = BRepTool::face_surface(src).expect("face surface");
        b.make_face(surf, &[wire])
    }
    // -----------------------------------------------------------------------
    // build_split_solids_full
    // -----------------------------------------------------------------------

    /// A reliable axis-aligned box with real geometry at position — unlike
    /// `make_box_corner`, whose lateral faces are broken. Modeled on the shared
    /// `unit_box` fixture but spanning `[x0,x1]×[y0,y1]×[z0,z1]`.
    fn axis_box(x0: f64, y0: f64, z0: f64, x1: f64, y1: f64, z1: f64) -> Solid {
        let b = TopoBuilder::new();
        let corners = [
            GpPnt::new(x0, y0, z0), GpPnt::new(x1, y0, z0), GpPnt::new(x1, y1, z0), GpPnt::new(x0, y1, z0),
            GpPnt::new(x0, y0, z1), GpPnt::new(x1, y0, z1), GpPnt::new(x1, y1, z1), GpPnt::new(x0, y1, z1),
        ];
        let vertices: Vec<Vertex> = corners.iter().map(|p| b.make_vertex(*p, 1e-7)).collect();
        let seg = |b: &TopoBuilder, p1: &GpPnt, p2: &GpPnt, v1: &Vertex, v2: &Vertex| {
            let dir = GpDir::from_vec(&GpVec::from_pnts(p1, p2)).unwrap();
            let lin = GpLin::from_pnt_dir(*p1, dir);
            let mut e = b.make_edge(Arc::new(GeomLine::new(lin)), 0.0, p1.distance(p2));
            b.add_edge_vertices(&mut e, v1, v2);
            e
        };
        let edge_idx: [(usize, usize); 12] = [
            (0, 1), (1, 2), (2, 3), (3, 0),
            (4, 5), (5, 6), (6, 7), (7, 4),
            (0, 4), (1, 5), (2, 6), (3, 7),
        ];
        let mut edges = Vec::new();
        for &(i, j) in &edge_idx {
            edges.push(seg(&b, &corners[i], &corners[j], &vertices[i], &vertices[j]));
        }
        let face_edge_sets: [[usize; 4]; 6] = [
            [0, 1, 2, 3], [4, 5, 6, 7], [0, 9, 4, 8], [2, 10, 6, 11], [3, 11, 7, 8], [1, 10, 5, 9],
        ];
        // Outward corner cycle of every face (same chords as the shared
        // `unit_box` fixture). An edge that the face closes against its
        // `edge_idx` chord is stored Reversed, so two faces sharing an edge
        // always see opposite orientations — the invariant
        // `BOPTools_AlgoTools::GetEdgeOff` (`BOPTools_AlgoTools.cxx:1099-1127`)
        // and `BOPAlgo_ShellSplitter::SplitBlock` (`BOPAlgo_ShellSplitter.cxx:319`)
        // rely on.
        let face_cycles: [[usize; 4]; 6] = [
            [0, 1, 2, 3], // bottom (z = z0)
            [4, 5, 6, 7], // top (z = z1)
            [0, 1, 5, 4], // front (y = y0)
            [3, 7, 6, 2], // back (y = y1)
            [0, 4, 7, 3], // left (x = x0)
            [1, 2, 6, 5], // right (x = x1)
        ];
        let face_planes: [(GpPnt, GpDir, GpDir); 6] = [
            (GpPnt::new(x0, y0, z0), GpDir::new(0.0, 0.0, -1.0).unwrap(), GpDir::new(0.0, 1.0, 0.0).unwrap()),
            (GpPnt::new(x0, y0, z1), GpDir::new(0.0, 0.0, 1.0).unwrap(), GpDir::new(1.0, 0.0, 0.0).unwrap()),
            (GpPnt::new(x0, y0, z0), GpDir::new(0.0, -1.0, 0.0).unwrap(), GpDir::new(1.0, 0.0, 0.0).unwrap()),
            (GpPnt::new(x0, y1, z0), GpDir::new(0.0, 1.0, 0.0).unwrap(), GpDir::new(0.0, 0.0, 1.0).unwrap()),
            (GpPnt::new(x0, y0, z0), GpDir::new(-1.0, 0.0, 0.0).unwrap(), GpDir::new(0.0, 0.0, 1.0).unwrap()),
            (GpPnt::new(x1, y0, z0), GpDir::new(1.0, 0.0, 0.0).unwrap(), GpDir::new(0.0, 1.0, 0.0).unwrap()),
        ];
        let mut faces = Vec::new();
        for fi in 0..6 {
            let (origin, normal, u_dir) = face_planes[fi];
            let ax3 = GpAx3::new(origin, normal, &u_dir).unwrap();
            let mut face = b.make_face_plane(&GpPln::new(ax3));
            let cycle = face_cycles[fi];
            let quad: Vec<Edge> = face_edge_sets[fi]
                .iter()
                .map(|&ei| {
                    let (i, j) = edge_idx[ei];
                    let k = (0..4)
                        .find(|&k| {
                            let a = cycle[k];
                            let b = cycle[(k + 1) % 4];
                            (a == i && b == j) || (a == j && b == i)
                        })
                        .expect("edge belongs to the face cycle");
                    let (a, b) = (cycle[k], cycle[(k + 1) % 4]);
                    let e = edges[ei].clone();
                    if (a, b) == (i, j) {
                        e
                    } else {
                        Edge(e.0.oriented(Orientation::Reversed))
                    }
                })
                .collect();
            let mut wire = b.make_wire(&quad);
            // The bottom plane's +X then +Y loop is clockwise in that UV frame;
            // reverse the wire so the outer ring is CCW (material on the left),
            // matching BRepPrimAPI_MakeBox / IntTools_FClass2d::IsHole.
            if fi == 0 {
                wire.0.reverse();
            }
            b.add_wire(&mut face, &wire);
            faces.push(face);
        }
        let shell = b.make_shell(&faces);
        b.make_solid(&[shell])
    }

    /// The face of `solid` whose every boundary vertex lies on the plane
    /// `coord = value` (0/1/2 → x/y/z) — the axis-aligned face at that plane.
    fn box_face_on(solid: &Solid, coord: usize, value: f64) -> Face {
        faces_of(&solid.0)
            .into_iter()
            .find(|f| {
                let vs = vertices_of(&f.0);
                !vs.is_empty()
                    && vs.iter().all(|v| {
                        let p = BRepTool::vertex_point(v);
                        let c = match coord {
                            0 => p.x(),
                            1 => p.y(),
                            _ => p.z(),
                        };
                        (c - value).abs() < 1e-9
                    })
            })
            .expect("face on coordinate plane")
    }

    /// Volume of `solid` from its vertex bounding box — reliable for the
    /// axis-aligned boxes this module's tests build (the mesh volume of
    /// translated/offset boxes is unreliable in this port).
    fn bbox_volume(solid: &TopoShape) -> f64 {
        let mut mn = (f64::INFINITY, f64::INFINITY, f64::INFINITY);
        let mut mx = (f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
        for v in vertices_of(solid) {
            let p = BRepTool::vertex_point(&v);
            mn.0 = mn.0.min(p.x());
            mn.1 = mn.1.min(p.y());
            mn.2 = mn.2.min(p.z());
            mx.0 = mx.0.max(p.x());
            mx.1 = mx.1.max(p.y());
            mx.2 = mx.2.max(p.z());
        }
        (mx.0 - mn.0) * (mx.1 - mn.1) * (mx.2 - mn.2)
    }

    #[test]
    fn classify_solid_state_in_on_out() {
        let a = axis_box(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let tol = 1e-7;
        assert_eq!(
            classify_solid_state(&a.0, &GpPnt::new(0.5, 0.5, 0.5), tol),
            FaceState::In
        );
        assert_eq!(
            classify_solid_state(&a.0, &GpPnt::new(2.0, 0.5, 0.5), tol),
            FaceState::Out
        );
        assert_eq!(
            classify_solid_state(&a.0, &GpPnt::new(0.5, 0.5, 1.0), tol),
            FaceState::On
        );
    }

    #[test]
    fn build_split_solids_full_selects_union_of_overlapping_boxes() {
        // Two unit boxes overlapping in x by 0.5: A=[0,1]³, B=[0.5,1.5]³.
        // The intersection cuts each box into two closed pieces; the two pieces
        // of the overlap region are geometrically identical, so only one of
        // them survives. Selected pieces: A[0,0.5] + overlap[0.5,1] +
        // B[1,1.5] = volume 1.5 (the Fuse union).
        let a = axis_box(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = axis_box(0.5, 0.0, 0.0, 1.5, 1.0, 1.0);
        let a1 = axis_box(0.0, 0.0, 0.0, 0.5, 1.0, 1.0);
        let a2 = axis_box(0.5, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b1 = axis_box(1.0, 0.0, 0.0, 1.5, 1.0, 1.0);
        let b2 = axis_box(0.5, 0.0, 0.0, 1.0, 1.0, 1.0);

        let mut ds = BopdsDS::new();
        ds.init(&[a.0.clone(), b.0.clone()]);
        let mut history = BopHistory::new();

        // Split images of A's faces: the pieces of A1 and A2 on each plane.
        history.add_image(&box_face_on(&a, 0, 0.0).0, box_face_on(&a1, 0, 0.0).0.clone());
        history.add_image(&box_face_on(&a, 0, 1.0).0, box_face_on(&a2, 0, 1.0).0.clone());
        for (coord, val) in [(1usize, 0.0f64), (1, 1.0), (2, 0.0), (2, 1.0)] {
            let src = box_face_on(&a, coord, val).0;
            history.add_image(&src, box_face_on(&a1, coord, val).0.clone());
            history.add_image(&src, box_face_on(&a2, coord, val).0.clone());
        }
        // Split images of B's faces.
        history.add_image(&box_face_on(&b, 0, 0.5).0, box_face_on(&b2, 0, 0.5).0.clone());
        history.add_image(&box_face_on(&b, 0, 1.5).0, box_face_on(&b1, 0, 1.5).0.clone());
        for (coord, val) in [(1usize, 0.0f64), (1, 1.0), (2, 0.0), (2, 1.0)] {
            let src = box_face_on(&b, coord, val).0;
            history.add_image(&src, box_face_on(&b1, coord, val).0.clone());
            history.add_image(&src, box_face_on(&b2, coord, val).0.clone());
        }

        let mut st = stub(ds, history, vec![a.0.clone(), b.0.clone()]);
        build_split_solids_full(&mut st, &[a.0.clone()], &[b.0.clone()], FaceState::Out, FaceState::Out)
            .unwrap();

        // A yields two pieces (its part outside B + the overlap); B yields one
        // (its part outside A) — the overlap region is a duplicate interior.
        let a_imgs = st.history().image(&a.0).expect("box A has split-solid images");
        let b_imgs = st.history().image(&b.0).expect("box B has split-solid images");
        assert_eq!(a_imgs.len(), 2, "A: outside piece + overlap piece");
        assert_eq!(b_imgs.len(), 1, "B: overlap piece is a duplicate and dropped");

        let total: f64 = a_imgs.iter().chain(b_imgs.iter()).map(bbox_volume).sum();
        assert!(
            (total - 1.5).abs() < 1e-6,
            "Fuse union volume 1.5, got {total}"
        );

        // The origins back-map is populated for every recorded piece.
        for im in a_imgs.iter().chain(b_imgs.iter()) {
            let ors = st.origins.get(&shape_key(im)).expect("origin recorded");
            assert!(ors.iter().any(|o| o.same_tshape(&a.0) || o.same_tshape(&b.0)));
        }
    }

    #[test]
    fn collect_candidates_dedupes_split_images_and_keeps_unsplit() {
        // A box whose top face is split into one fresh image: the global
        // candidate list holds the split image (not the original top face) plus
        // the five unsplit faces, all distinct (the fence dedupes).
        let boxed = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let faces = faces_of(&boxed.solid.0);
        let top = face_at_z(&faces, 1.0);
        let new_top = re_face(&top);

        let mut ds = BopdsDS::new();
        ds.init(&[boxed.solid.0.clone()]);
        let mut history = BopHistory::new();
        history.add_image(&top.0, new_top.0.clone());
        let b = stub(ds, history, vec![boxed.solid.0.clone()]);

        let candidates = collect_all_candidate_faces(&b);
        assert_eq!(candidates.len(), 6, "one split image + five unsplit faces");
        assert!(
            candidates.iter().any(|c| c.same_tshape(&new_top.0)),
            "split image present"
        );
        assert!(
            !candidates.iter().any(|c| c.same_tshape(&top.0)),
            "original split face replaced by its image"
        );
        for (i, a) in candidates.iter().enumerate() {
            for b2 in &candidates[i + 1..] {
                assert!(!a.same_tshape(b2), "candidates are distinct");
            }
        }
    }

    #[test]
    fn classify_faces_in_solid_culls_remote_and_keeps_internal() {
        // Solid box [0,1]³; a face far outside is rejected by the pairwise box
        // cull, a face strictly inside (a standalone candidate, not a face of
        // the solid) is classified IN and returned FORWARD + REVERSED.
        let boxed = axis_box(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let far = axis_box(5.0, 5.0, 5.0, 6.0, 6.0, 6.0);
        let far_face = faces_of(&far.0)[0].clone();
        let inner = axis_box(0.25, 0.25, 0.25, 0.75, 0.75, 0.75);
        let inner_face = faces_of(&inner.0)[0].clone();

        let candidates = vec![far_face.0.clone(), inner_face.0.clone()];
        let own_edges: HashSet<EKey> = HashSet::new();
        let own_faces: HashSet<usize> =
            faces_of(&boxed.0).into_iter().map(|fc| shape_key(&fc.0)).collect();
        let solid_box = shape_bbox(&boxed.0);
        let tol = 1e-7;

        let in_faces = classify_faces_in_solid(
            &candidates,
            &boxed.0,
            &own_faces,
            &[],
            &own_edges,
            &solid_box,
            tol,
        );
        assert!(
            !in_faces.iter().any(|c| c.same_tshape(&far_face.0)),
            "remote face box-culled"
        );
        assert!(
            in_faces.iter().any(|c| c.same_tshape(&inner_face.0)),
            "inner face classified IN"
        );
        assert_eq!(
            in_faces.iter().filter(|c| c.same_tshape(&inner_face.0)).count(),
            2,
            "inner face kept FORWARD and REVERSED"
        );
    }

    #[test]
    fn connexity_blocks_group_via_non_solid_edges_and_split_at_barriers() {
        // Two adjacent faces of a box share an edge: with no solid edges they
        // form one connexity block (connected through the shared non-barrier
        // edge) classified by the block-start representative; once that shared
        // edge is marked as a solid boundary edge it becomes a barrier and the
        // two faces split into separate singleton blocks, each its own
        // representative (mirroring `BOPAlgo_FillIn3DParts::MakeConnexityBlock`).
        let boxed = axis_box(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let fs = faces_of(&boxed.0);
        let (a, b, shared) = {
            let ka: HashSet<EKey> = face_edges(&Face(fs[0].0.clone())).into_iter().collect();
            let mut bf = None;
            let mut sk = None;
            for f in &fs[1..] {
                for k in face_edges(&Face(f.0.clone())) {
                    if ka.contains(&k) {
                        bf = Some(f.0.clone());
                        sk = Some(k);
                        break;
                    }
                }
                if bf.is_some() {
                    break;
                }
            }
            (fs[0].0.clone(), bf.expect("adjacent face"), sk.expect("shared edge"))
        };
        // No solid edges: one block connected through the shared edge.
        let blocks = connexity_blocks(&[a.clone(), b.clone()], &HashSet::new());
        assert_eq!(blocks.len(), 1, "two faces connected through a non-solid edge");
        let (bfaces, brep) = &blocks[0];
        assert_eq!(bfaces.len(), 2);
        assert!(brep.same_tshape(&a), "representative is the block start (no barrier edge)");
        // The shared edge as a solid boundary: the barrier splits the two faces.
        let mut se: HashSet<EKey> = HashSet::new();
        se.insert(shared);
        let blocks = connexity_blocks(&[a.clone(), b.clone()], &se);
        assert_eq!(blocks.len(), 2, "solid boundary edge splits the connexity block");
        for (bfaces, brep) in &blocks {
            assert_eq!(bfaces.len(), 1, "each barrier-separated face is its own block");
            assert!(brep.same_tshape(&bfaces[0]), "a barrier-edge face is its own representative");
        }
    }

    #[test]
    fn build_split_solids_full_keeps_disjoint_boxes() {
        // Two fully disjoint boxes: neither lies inside the other, so both are
        // kept as their whole self (volume 2).
        let a = axis_box(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = axis_box(2.0, 0.0, 0.0, 3.0, 1.0, 1.0);

        let mut ds = BopdsDS::new();
        ds.init(&[a.0.clone(), b.0.clone()]);
        let history = BopHistory::new();
        // Nothing split: no images -> no solid is interfered, none gets an image.
        let mut st = stub(ds, history, vec![a.0.clone(), b.0.clone()]);
        build_split_solids_full(&mut st, &[a.0.clone()], &[b.0.clone()], FaceState::Out, FaceState::Out)
            .unwrap();
        assert!(
            !st.history().has_image(&a.0) && !st.history().has_image(&b.0),
            "disjoint unsplit boxes keep no solid image"
        );
    }
