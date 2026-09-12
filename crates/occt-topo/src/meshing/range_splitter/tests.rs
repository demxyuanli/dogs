use super::prelude::*;
use super::*;
    use std::sync::Arc;

    use occt_core::gp::{GpAx3, GpCylinder, GpPln, GpSphere, GpTorus};
    use occt_geom::{
        bspline_surface::{bspline_surface_uniform_knots, GeomBSplineSurface},
        GeomCylinder, GeomPlane, GeomSphere, GeomTorus,
    };

    use crate::builder::TopoBuilder;
    use crate::tgeometry::GeometryRegistry;

    /// Releases registry entries for a face so tests don't leave stale geometry.
    fn clear_face(mf: &MeshFace) {
        GeometryRegistry::global().clear_shape(&mf.face().0);
    }

    /// Builds a discrete face whose registered plane surface is replaced by the
    /// given analytic surface; deflection fixed to 0.001.
    fn make_dface(surface: Arc<dyn Surface>) -> MeshFace {
        let b = TopoBuilder::new();
        let face = b.make_face_plane(&GpPln::new(GpAx3::standard()));
        let mut mf = MeshFace::new(face);
        mf.set_surface(Some(surface));
        mf.set_deflection(0.001);
        mf
    }

    fn plane_surface() -> Arc<dyn Surface> {
        Arc::new(GeomPlane::new(GpPln::new(GpAx3::standard())))
    }

    fn cylinder_surface(r: f64) -> Arc<dyn Surface> {
        Arc::new(GeomCylinder::new(GpCylinder::new(GpAx3::standard(), r).unwrap()))
    }

    fn sphere_surface(r: f64) -> Arc<dyn Surface> {
        Arc::new(GeomSphere::new(GpSphere::new(GpAx3::standard(), r).unwrap()))
    }

    fn torus_surface(major: f64, minor: f64) -> Arc<dyn Surface> {
        Arc::new(GeomTorus::new(GpTorus::new(GpAx3::standard(), major, minor).unwrap()))
    }

    #[test]
    fn classify_analytic_surfaces() {
        assert_eq!(classify_surface(plane_surface().as_ref()), SurfaceType::Plane);
        assert_eq!(classify_surface(cylinder_surface(1.0).as_ref()), SurfaceType::Cylinder);
        assert_eq!(classify_surface(sphere_surface(1.0).as_ref()), SurfaceType::Sphere);
        assert_eq!(classify_surface(torus_surface(2.0, 1.0).as_ref()), SurfaceType::Torus);

        let (ku, kv) = bspline_surface_uniform_knots(4, 4, 3, 3);
        let poles = (0..4)
            .map(|i| (0..4).map(|j| GpPnt::new(i as f64, j as f64, 0.0)).collect())
            .collect();
        let bs: Arc<dyn Surface> = Arc::new(GeomBSplineSurface::new(poles, ku, kv, 3, 3).unwrap());
        assert_eq!(classify_surface(bs.as_ref()), SurfaceType::BSplineSurface);
    }

    #[test]
    fn factory_dispatches_by_surface_type() {
        let _ = create_range_splitter(plane_surface().as_ref());
        let _ = create_range_splitter(cylinder_surface(1.0).as_ref());
        let _ = create_range_splitter(sphere_surface(1.0).as_ref());
        let _ = create_range_splitter(torus_surface(2.0, 1.0).as_ref());
    }

    #[test]
    fn plane_range_splitter_computes_ranges_and_scale() {
        let dface = make_dface(plane_surface());
        let params = MeshParameters::default();
        let mut sp = DefaultRangeSplitter::new();
        sp.reset(&dface, &params);
        sp.add_point(GpPnt2d::new(0.0, 0.0));
        sp.add_point(GpPnt2d::new(1.0, 0.0));
        sp.add_point(GpPnt2d::new(0.0, 1.0));
        sp.adjust_range();
        assert!(sp.is_valid());

        let (u0, u1) = sp.range_u();
        let (v0, v1) = sp.range_v();
        assert!((u0 - 0.0).abs() < 1e-9 && (u1 - 1.0).abs() < 1e-9);
        assert!((v0 - 0.0).abs() < 1e-9 && (v1 - 1.0).abs() < 1e-9);

        // Face-basis scaling maps range origin to ~(0, 0).
        let s = sp.scale(GpPnt2d::new(u0, v0), true);
        assert!(s.x().abs() < 1e-9 && s.y().abs() < 1e-9);
        let back = sp.scale(GpPnt2d::new(0.0, 0.0), false);
        assert!((back.x() - u0).abs() < 1e-9 && (back.y() - v0).abs() < 1e-9);

        // The base splitter generates no interior nodes.
        assert!(sp.generate_surface_nodes(&params).is_none());

        // 3D point evaluation on the Z=0 plane.
        let p = sp.point(GpPnt2d::new(0.5, 0.5));
        assert!((p.z() - 0.0).abs() < 1e-9);
        clear_face(&dface);
    }

    #[test]
    fn cylinder_seam_range_clamps_to_period() {
        let dface = make_dface(cylinder_surface(1.0));
        let params = MeshParameters::default();
        let mut sp = CylinderRangeSplitter::new();
        sp.reset(&dface, &params);

        // Boundary points spanning more than one U-period cross the seam: the
        // discrete range must be clamped to a single period.
        sp.add_point(GpPnt2d::new(-0.2, 0.0));
        sp.add_point(GpPnt2d::new(6.5, 0.0));
        sp.add_point(GpPnt2d::new(0.0, -1.0));
        sp.add_point(GpPnt2d::new(1.0, 1.0));
        sp.adjust_range();
        assert!(sp.is_valid());

        let (u0, u1) = sp.range_u();
        let period = 2.0 * PI;
        assert!(u1 - u0 <= period + 1e-9, "seam range {u0}..{u1} exceeds one period");
        assert!((u1 - u0 - period).abs() < 1e-9);
        assert!((u0 - (-0.2)).abs() < 1e-9);

        // Cylinder delta: first component from the angular step, V delta = 1.
        let (du, dv) = sp.delta();
        assert!((dv - 1.0).abs() < 1e-12);
        assert!(du > 0.0);

        // Faithful to OCCT: the V-step code is commented out, so no interior rows.
        let nodes = sp.generate_surface_nodes(&params).unwrap_or_default();
        assert!(nodes.iter().all(|p| p.x() >= u0 - 1e-9 && p.x() <= u1 + 1e-9));
        assert!(nodes.iter().all(|p| p.y() >= -1.0 - 1e-9 && p.y() <= 1.0 + 1e-9));
        clear_face(&dface);
    }

    #[test]
    fn sphere_range_generates_staggered_nodes() {
        let dface = make_dface(sphere_surface(1.0));
        let params = MeshParameters::default();
        let mut sp = SphereRangeSplitter::new();
        sp.reset(&dface, &params);
        // Full sphere range.
        sp.add_point(GpPnt2d::new(0.0, -PI * 0.5));
        sp.add_point(GpPnt2d::new(2.0 * PI, PI * 0.5));
        sp.adjust_range();
        assert!(sp.is_valid());

        let nodes = sp.generate_surface_nodes(&params).expect("sphere nodes");
        assert!(nodes.len() > 100, "expected a dense staggered grid, got {}", nodes.len());

        let (u0, u1) = sp.range_u();
        let (v0, v1) = sp.range_v();
        for p in &nodes {
            assert!(p.x() >= u0 - 1e-6 && p.x() <= u1 + 1e-6);
            assert!(p.y() >= v0 - 1e-6 && p.y() <= v1 + 1e-6);
        }
        // Staggered rows: consecutive V rows start at the range origin and at
        // origin + half a U step, alternating — so at least two distinct
        // per-row minimum U values must exist.
        let mut map = std::collections::BTreeMap::<i64, f64>::new();
        for p in &nodes {
            let key = (p.y() * 1e6).round() as i64;
            let e = map.entry(key).or_insert(p.x());
            if p.x() < *e {
                *e = p.x();
            }
        }
        let mut row_mins: Vec<f64> = map.values().copied().collect();
        row_mins.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let staggered = row_mins.windows(2).any(|w| (w[1] - w[0]).abs() > 1e-3);
        assert!(staggered, "sphere rows should be staggered by half a step");
        clear_face(&dface);
    }

    #[test]
    fn torus_range_generates_nodes_within_range() {
        let dface = make_dface(torus_surface(2.0, 1.0));
        let params = MeshParameters::default();
        let mut sp = TorusRangeSplitter::new();
        sp.reset(&dface, &params);
        // A quarter of the torus in U, full V, with a dense boundary (as a
        // discretized face boundary would supply) so the density filter keeps
        // interior samples.
        for i in 0..=16 {
            let u = PI * 0.5 * i as f64 / 16.0;
            for j in 0..=32 {
                let v = 2.0 * PI * j as f64 / 32.0;
                sp.add_point(GpPnt2d::new(u, v));
            }
        }
        sp.adjust_range();
        assert!(sp.is_valid());

        let (u0, u1) = sp.range_u();
        let (v0, v1) = sp.range_v();
        let nodes = sp.generate_surface_nodes(&params).expect("torus nodes");
        assert!(nodes.len() > 100, "expected interior nodes, got {}", nodes.len());
        for p in &nodes {
            assert!(p.x() >= u0 - 1e-6 && p.x() <= u1 + 1e-6);
            assert!(p.y() >= v0 - 1e-6 && p.y() <= v1 + 1e-6);
        }
        clear_face(&dface);
    }

    #[test]
    fn nurbs_range_generates_grid_within_bounds() {
        let (ku, kv) = bspline_surface_uniform_knots(4, 4, 3, 3);
        let poles = (0..4)
            .map(|i| (0..4).map(|j| GpPnt::new(i as f64, j as f64, 0.0)).collect())
            .collect();
        let bs: Arc<dyn Surface> = Arc::new(GeomBSplineSurface::new(poles, ku, kv, 3, 3).unwrap());
        let dface = make_dface(bs);
        let params = MeshParameters::default();

        let mut sp = NURBSRangeSplitter::new();
        sp.reset(&dface, &params);
        sp.add_point(GpPnt2d::new(0.0, 0.0));
        sp.add_point(GpPnt2d::new(1.0, 1.0));
        sp.adjust_range();
        assert!(sp.is_valid());

        let (u0, u1) = sp.range_u();
        let (v0, v1) = sp.range_v();
        assert!((u0 - 0.0).abs() < 1e-9 && (u1 - 1.0).abs() < 1e-9);
        assert!((v0 - 0.0).abs() < 1e-9 && (v1 - 1.0).abs() < 1e-9);

        let nodes = sp.generate_surface_nodes(&params).expect("nurbs nodes");
        assert!(nodes.len() >= 4, "expected a grid, got {}", nodes.len());
        for p in &nodes {
            assert!(p.x() >= u0 - 1e-6 && p.x() <= u1 + 1e-6);
            assert!(p.y() >= v0 - 1e-6 && p.y() <= v1 + 1e-6);
        }
        clear_face(&dface);
    }

    #[test]
    fn boundary_params_splitter_collects_uv_params() {
        let dface = make_dface(plane_surface());
        let params = MeshParameters::default();
        // The plain UV splitter does not seed its maps from AddPoint (faithful
        // to OCCT); only the boundary/torus splitters do.
        let mut sp = UVParamRangeSplitter::new();
        sp.reset(&dface, &params);
        sp.add_point(GpPnt2d::new(0.5, 0.25));
        assert!(sp.parameters_u().unwrap().is_empty());

        let mut bp = BoundaryParamsRangeSplitter::new();
        bp.reset(&dface, &params);
        bp.add_point(GpPnt2d::new(0.5, 0.25));
        bp.add_point(GpPnt2d::new(0.5, 0.75));
        bp.add_point(GpPnt2d::new(1.5, 0.25));
        assert!(bp.parameters_u().unwrap().contains(&0.5));
        assert!(bp.parameters_u().unwrap().contains(&1.5));
        assert_eq!(bp.parameters_v().unwrap().len(), 2);
        // Reset clears the collected parameters.
        bp.reset(&dface, &params);
        assert!(bp.parameters_u().unwrap().is_empty());
        clear_face(&dface);
    }
