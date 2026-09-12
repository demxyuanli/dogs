use super::prelude::*;
use super::*;
    use std::f64::consts::PI;
    use std::sync::Arc;

    use occt_core::gp::{GpAx3, GpPln, GpSphere};
    use occt_geom::{GeomPlane, GeomSphere};

    /// A self-contained mesh driving the refinement loop against the real
    /// `Delaun` by batch re-triangulation (OCCT inserts nodes incrementally via
    /// `DelaunayNodeInsertionMeshAlgo`; this wrapper rebuilds the triangulation
    /// with the accumulated vertex list, which exercises the same deflection
    /// control logic).
    struct RefineMesh {
        vertices: Vec<DelaunVertex>,
        delaun: Delaun,
    }

    impl RefineMesh {
        fn from_points(points: &[DelaunVertex]) -> Self {
            Self {
                vertices: points.to_vec(),
                delaun: Delaun::new_vertices(points),
            }
        }

        fn insert_control(&mut self, surface: &dyn Surface, nodes: &[GpPnt2d]) -> bool {
            let before = self.vertices.len();
            for &uv in nodes {
                let p3d = surface.d0(uv.x(), uv.y());
                self.vertices.push(DelaunVertex::new(uv, p3d, self.vertices.len() as i32, VertexState::Free));
            }
            if self.vertices.len() == before {
                return false;
            }
            self.delaun = Delaun::new_vertices(&self.vertices);
            true
        }

        fn triangle_count(&self) -> usize {
            self.delaun.result().elements_of_domain().len()
        }

        fn live_vertex_ids(&self) -> Vec<i32> {
            let mut ids = HashSet::new();
            for &id in self.delaun.result().elements_of_domain() {
                let tri = self.delaun.get_triangle(id);
                for v in tri.vertex_indices {
                    ids.insert(v);
                }
            }
            ids.into_iter().collect()
        }
    }

    impl DeflectionMesh for RefineMesh {
        fn elements_of_domain(&self) -> Vec<i32> {
            self.delaun.result().elements_of_domain().iter().copied().collect()
        }
        fn triangle(&self, id: i32) -> DelaunTriangle {
            self.delaun.get_triangle(id)
        }
        fn link(&self, id: i32) -> DelaunLink {
            self.delaun.get_edge(id)
        }
        fn link_movability(&self, id: i32) -> VertexState {
            self.delaun.result().link_movability(id)
        }
        fn vertex(&self, id: i32) -> DelaunVertex {
            self.delaun.get_vertex(id)
        }
        fn shot_triangles(&self, p: occt_core::gp::GpXY) -> Option<Vec<i32>> {
            Some(self.delaun.circles().select(p))
        }
    }

    fn plane() -> Arc<dyn Surface> {
        Arc::new(GeomPlane::new(GpPln::new(GpAx3::standard())))
    }

    fn unit_sphere() -> Arc<dyn Surface> {
        Arc::new(GeomSphere::new(GpSphere::new(GpAx3::standard(), 1.0).unwrap()))
    }

    /// 3x3 grid on the plane z = 0, UV in [0, 2] x [0, 2].
    fn plane_grid(surface: &dyn Surface) -> Vec<DelaunVertex> {
        uv_grid(surface, 3, 3, (0.0, 2.0), (0.0, 2.0))
    }

    /// 3x3 grid over the full parametric domain of a unit sphere.
    fn sphere_grid(surface: &dyn Surface) -> Vec<DelaunVertex> {
        uv_grid(surface, 3, 3, (0.0, 2.0 * PI), (-PI / 2.0, PI / 2.0))
    }

    fn uv_grid(
        surface: &dyn Surface,
        nu: usize,
        nv: usize,
        (u0, u1): (f64, f64),
        (v0, v1): (f64, f64),
    ) -> Vec<DelaunVertex> {
        let du = u1 - u0;
        let dv = v1 - v0;
        let mut points = Vec::new();
        for i in 0..nu {
            for j in 0..nv {
                let uv = GpPnt2d::new(
                    u0 + du * (i as f64 + 0.5) / nu as f64,
                    v0 + dv * (j as f64 + 0.5) / nv as f64,
                );
                points.push(DelaunVertex::new(uv, surface.d0(uv.x(), uv.y()), 0, VertexState::Free));
            }
        }
        points
    }

    #[test]
    fn normal_deviation_is_distance_to_plane() {
        let ref_pnt = GpPnt::new(0.0, 0.0, 0.0);
        let functor = NormalDeviation::new(&ref_pnt, GpVec::new(0.0, 0.0, 1.0));
        // Points in the plane z = 0 have zero deviation.
        assert!(functor.square_deviation(&GpPnt::new(1.0, 2.0, 0.0)) < 1e-15);
        // Points offset by `d` normal to the plane have deviation d^2.
        assert!((functor.square_deviation(&GpPnt::new(1.0, 2.0, 3.0)) - 9.0).abs() < 1e-12);
        assert!((functor.square_deviation(&GpPnt::new(1.0, 2.0, -0.5)) - 0.25).abs() < 1e-12);
    }

    #[test]
    fn line_deviation_is_distance_to_chord() {
        let a = GpPnt::new(0.0, 0.0, 0.0);
        let b = GpPnt::new(2.0, 0.0, 0.0);
        let functor = LineDeviation::new(&a, &b);
        // On the chord.
        assert!(functor.square_deviation(&GpPnt::new(1.0, 0.0, 0.0)) < 1e-15);
        // Off the chord by 1 unit.
        assert!((functor.square_deviation(&GpPnt::new(1.0, 1.0, 0.0)) - 1.0).abs() < 1e-12);
        // Collinear beyond an endpoint: distance to the supporting line, not the
        // segment, so a point on the line still reads zero.
        assert!(functor.square_deviation(&GpPnt::new(3.0, 0.0, 0.0)) < 1e-15);
    }

    #[test]
    fn plane_triangulation_is_not_refined() {
        // A flat plane has zero linear and angular deflection everywhere: the
        // refinement must leave the 3x3 grid untouched.
        let surface = plane();
        let params = MeshParameters::default();

        let points = plane_grid(surface.as_ref());
        let mut mesh = RefineMesh::from_points(&points);
        let n0 = mesh.triangle_count();
        assert_eq!(n0, 8, "3x3 grid triangulates to 8 triangles");

        let mut algo = DelaunayDeflectionControlMeshAlgo::new(&params);
        let max_deflection = algo.post_process(
            &mut mesh,
            surface.as_ref(),
            0.001,
            &mut |m, nodes| m.insert_control(surface.as_ref(), nodes),
        );

        assert_eq!(mesh.triangle_count(), n0, "plane must not refine");
        assert!(max_deflection < 1e-9, "plane max deflection must be ~0, got {max_deflection}");
    }

    #[test]
    fn sphere_triangulation_is_refined() {
        // A coarse sphere triangulation has huge chord sagitta: the refinement
        // must split links and grow the triangle count.
        let surface = unit_sphere();
        let mut params = MeshParameters::default();
        params.deflection = 0.01;
        params.deflection_interior = 0.01;

        let points = sphere_grid(surface.as_ref());
        let mut mesh = RefineMesh::from_points(&points);
        let n0 = mesh.triangle_count();
        assert_eq!(n0, 8);

        let mut algo = DelaunayDeflectionControlMeshAlgo::new(&params);
        let max_deflection = algo.post_process(
            &mut mesh,
            surface.as_ref(),
            0.01,
            &mut |m, nodes| m.insert_control(surface.as_ref(), nodes),
        );

        assert!(
            mesh.triangle_count() > n0,
            "sphere must refine: {n0} -> {}",
            mesh.triangle_count()
        );
        assert!(max_deflection > 0.0, "sphere deflection must be positive");
    }

    #[test]
    fn sphere_refined_vertices_lie_on_surface() {
        // Every node generated by the refinement must sit on the unit sphere.
        let surface = unit_sphere();
        let mut params = MeshParameters::default();
        params.deflection = 0.01;
        params.deflection_interior = 0.01;

        let points = sphere_grid(surface.as_ref());
        let mut mesh = RefineMesh::from_points(&points);

        let mut algo = DelaunayDeflectionControlMeshAlgo::new(&params);
        algo.post_process(
            &mut mesh,
            surface.as_ref(),
            0.01,
            &mut |m, nodes| m.insert_control(surface.as_ref(), nodes),
        );

        let ids = mesh.live_vertex_ids();
        assert!(ids.len() > points.len(), "refinement must add vertices");
        for id in ids {
            let vertex = mesh.delaun.get_vertex(id);
            let radius = vertex.p3d.distance(&GpPnt::zero());
            assert!(
                (radius - 1.0).abs() < 1e-9,
                "vertex off unit sphere: radius {radius} at ({}, {}, {})",
                vertex.p3d.x(),
                vertex.p3d.y(),
                vertex.p3d.z()
            );
        }
    }

    #[test]
    fn large_min_size_blocks_sphere_refinement() {
        // A MinSize of half the sphere diameter (2.0) makes every candidate
        // control point closer than MinSize to an existing node, so all splits
        // are rejected: the triangle count must stay unchanged.
        let surface = unit_sphere();
        let mut params = MeshParameters::default();
        params.min_size = 2.0;

        let points = sphere_grid(surface.as_ref());
        let mut mesh = RefineMesh::from_points(&points);
        let n0 = mesh.triangle_count();

        let mut algo = DelaunayDeflectionControlMeshAlgo::new(&params);
        algo.post_process(
            &mut mesh,
            surface.as_ref(),
            0.001,
            &mut |m, nodes| m.insert_control(surface.as_ref(), nodes),
        );

        assert_eq!(mesh.triangle_count(), n0, "min_size must reject every split");
    }

    #[test]
    fn post_process_reads_mesh_and_proposes_control_nodes() {
        // Drive the algorithm against the real `Delaun` directly: even with a
        // recording (non-mutating) insert callback the curved sphere mesh must
        // produce control nodes and a positive max deflection.
        let surface = unit_sphere();
        let mut params = MeshParameters::default();
        params.deflection = 0.01;
        params.deflection_interior = 0.01;

        let points = sphere_grid(surface.as_ref());
        let mut delaun = Delaun::new_vertices(&points);
        let n0 = delaun.result().elements_of_domain().len();

        let mut algo = DelaunayDeflectionControlMeshAlgo::new(&params);
        let mut collected: Vec<GpPnt2d> = Vec::new();
        let max_deflection = algo.post_process(
            &mut delaun,
            surface.as_ref(),
            0.01,
            &mut |_mesh, nodes| {
                collected.extend_from_slice(nodes);
                !nodes.is_empty()
            },
        );

        assert!(!collected.is_empty(), "sphere must propose control nodes");
        assert!(max_deflection > 0.0);
        // The recording callback does not mutate the mesh.
        assert_eq!(delaun.result().elements_of_domain().len(), n0);
    }
