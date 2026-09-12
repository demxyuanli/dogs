use super::prelude::*;
use super::*;
    use std::f64::consts::PI;
    use std::sync::Arc;

    use occt_core::gp::{GpAx3, GpCylinder, GpPln};
    use occt_geom::bspline_surface::{bspline_surface_uniform_knots, GeomBSplineSurface};
    use occt_geom::{GeomCylinder, GeomPlane};

    use crate::brep_extrema::test_box::unit_box;
    use crate::tgeometry::GeometryRegistry;

    /// Test helper: an unchecked range.
    fn r(a: f64, b: f64) -> IntRange {
        IntRange::new_unchecked(a, b)
    }

    fn plane_surface() -> Arc<dyn Surface> {
        Arc::new(GeomPlane::new(GpPln::new(GpAx3::standard())))
    }

    fn cylinder_surface(radius: f64) -> Arc<dyn Surface> {
        Arc::new(GeomCylinder::new(GpCylinder::new(GpAx3::standard(), radius).unwrap()))
    }

    /// A degree-3 BSpline patch over `[0,3]²` that is curved along U (a cubic
    /// in `u`) and constant along V — a ruled "extrusion" of a cubic curve.
    fn curved_bspline() -> Arc<dyn Surface> {
        let (ku, kv) = bspline_surface_uniform_knots(4, 4, 3, 3);
        let poles: Vec<Vec<GpPnt>> = (0..4)
            .map(|i| {
                let u = i as f64;
                (0..4).map(|j| GpPnt::new(u, j as f64, u * u)).collect()
            })
            .collect();
        Arc::new(GeomBSplineSurface::new(poles, ku, kv, 3, 3).unwrap())
    }

    // ---- BaseRangeSample ----

    #[test]
    fn base_range_sample_depth_access() {
        let mut s = BaseRangeSample::new();
        assert_eq!(s.get_depth(), 0);
        s.set_depth(3);
        assert_eq!(s.get_depth(), 3);
        let s2 = BaseRangeSample::with_depth(5);
        assert_eq!(s2.get_depth(), 5);
        assert_eq!(BaseRangeSample::default().get_depth(), 0);
    }

    // ---- CurveRangeSample ----

    #[test]
    fn curve_range_sample_construct_depth_index_equal() {
        let a = CurveRangeSample::with_index_depth(7, 2);
        assert_eq!(a.get_depth(), 2);
        assert_eq!(a.get_index(), 7);
        let b = CurveRangeSample::with_index_depth(7, 2);
        assert!(a.is_equal(&b));
        assert_eq!(a, b);
        let c = CurveRangeSample::with_index_depth(7, 3);
        assert!(!a.is_equal(&c), "different depth");
        let d = CurveRangeSample::with_index_depth(8, 2);
        assert!(!a.is_equal(&d), "different index");
        let e = CurveRangeSample::with_index(4);
        assert_eq!(e.get_depth(), 0);
        assert_eq!(e.get_index(), 4);
    }

    #[test]
    fn curve_range_sample_get_range_depth_zero_covers_domain() {
        let s = CurveRangeSample::new();
        assert_eq!(s.get_range(-2.0, 5.0, 10), r(-2.0, 5.0));
    }

    #[test]
    fn curve_range_sample_get_range_subdivides() {
        // Depth 1 with 4 samples splits [0, 1] into 4 quarters.
        let s = CurveRangeSample::with_index_depth(2, 1);
        assert_eq!(s.get_range(0.0, 1.0, 4), r(0.5, 0.75));
        // Index 0 -> first quarter.
        let s0 = CurveRangeSample::with_index_depth(0, 1);
        assert_eq!(s0.get_range(0.0, 1.0, 4), r(0.0, 0.25));
        // Last index -> last quarter (upper bound inclusive).
        let s3 = CurveRangeSample::with_index_depth(3, 1);
        assert_eq!(s3.get_range(0.0, 1.0, 4), r(0.75, 1.0));
        // Depth 2 squares the sample count: 4^2 = 16 intervals of length 1/16.
        let sd = CurveRangeSample::with_index_depth(5, 2);
        assert_eq!(sd.get_range(0.0, 1.0, 4), r(5.0 / 16.0, 6.0 / 16.0));
    }

    #[test]
    fn curve_range_sample_index_deeper() {
        let s = CurveRangeSample::with_index_depth(3, 1);
        assert_eq!(s.get_range_index_deeper(4), 12);
        assert_eq!(CurveRangeSample::with_index(0).get_range_index_deeper(10), 0);
    }

    #[test]
    fn curve_range_sample_setters() {
        let mut s = CurveRangeSample::new();
        s.set_index(2);
        s.set_depth(1);
        assert_eq!(s, CurveRangeSample::with_index_depth(2, 1));
    }

    // ---- SurfaceRangeSample ----

    #[test]
    fn surface_range_sample_construct_and_accessors() {
        let s = SurfaceRangeSample::with_indexes_depths(2, 1, 3, 2);
        assert_eq!(s.get_index_u(), 2);
        assert_eq!(s.get_depth_u(), 1);
        assert_eq!(s.get_index_v(), 3);
        assert_eq!(s.get_depth_v(), 2);
        assert_eq!(s.get_indexes(), (2, 3));
        assert_eq!(s.get_depths(), (1, 2));

        let mut m = SurfaceRangeSample::new();
        m.set_index_u(5);
        m.set_index_v(6);
        m.set_depth_u(1);
        m.set_depth_v(1);
        assert_eq!(m.get_indexes(), (5, 6));
        assert_eq!(m.get_depths(), (1, 1));

        let (ru, rv) = m.get_ranges();
        assert_eq!(ru, CurveRangeSample::with_index_depth(5, 1));
        assert_eq!(rv, CurveRangeSample::with_index_depth(6, 1));
    }

    #[test]
    fn surface_range_sample_equality() {
        let a = SurfaceRangeSample::with_indexes_depths(1, 1, 2, 1);
        let b = SurfaceRangeSample::with_indexes_depths(1, 1, 2, 1);
        let c = SurfaceRangeSample::with_indexes_depths(1, 1, 2, 2);
        assert!(a.is_equal(&b));
        assert_eq!(a, b);
        assert!(!a.is_equal(&c), "V depth differs");
    }

    #[test]
    fn surface_range_sample_get_ranges_and_deeper() {
        let s = SurfaceRangeSample::with_indexes_depths(1, 1, 2, 1);
        assert_eq!(s.get_range_u(0.0, 1.0, 4), r(0.25, 0.5));
        assert_eq!(s.get_range_v(0.0, 2.0, 4), r(1.0, 1.5));
        let (ur, vr) = s.get_range(0.0, 1.0, 4, 0.0, 2.0, 4);
        assert_eq!(ur, r(0.25, 0.5));
        assert_eq!(vr, r(1.0, 1.5));
        assert_eq!(s.get_range_index_u_deeper(4), 4);
        assert_eq!(s.get_range_index_v_deeper(4), 8);
    }

    #[test]
    fn surface_range_sample_from_ranges() {
        let ru = CurveRangeSample::with_index_depth(1, 1);
        let rv = CurveRangeSample::with_index_depth(2, 1);
        let s = SurfaceRangeSample::from_ranges(ru, rv);
        assert_eq!(s.get_sample_range_u(), ru);
        assert_eq!(s.get_sample_range_v(), rv);
        s.get_ranges();
    }

    // ---- CurveRangeLocalizeData ----

    #[test]
    fn curve_localize_data_build_maps_intervals() {
        let mut ld = CurveRangeLocalizeData::new(4, 1e-6);
        ld.set_root_index(2);
        let n = ld.build(r(0.0, 1.0), 7);
        assert_eq!(n, 4);
        assert_eq!(ld.get_nb_sample(), 4);
        assert_eq!(ld.get_root_index(), 2);
        assert_eq!(ld.ranges().len(), 4);
        assert_eq!(ld.range(0), Some(r(0.0, 0.25)));
        assert_eq!(ld.range(1), Some(r(0.25, 0.5)));
        assert_eq!(ld.range(2), Some(r(0.5, 0.75)));
        assert_eq!(ld.range(3), Some(r(0.75, 1.0)));
        for i in 0..4 {
            assert_eq!(ld.curve_index(i), Some(7), "all cells mapped to curve 7");
        }
        // Every parameter is located in exactly one cell.
        assert_eq!(ld.find_index(0.0), Some(0));
        assert_eq!(ld.find_index(0.3), Some(1));
        assert_eq!(ld.find_index(0.99), Some(3));
        assert_eq!(ld.find_index(1.0), Some(3), "upper bound inclusive");
        assert_eq!(ld.find_index(-0.1), None);
    }

    #[test]
    fn curve_localize_data_reassign_and_out() {
        let mut ld = CurveRangeLocalizeData::new(3, 0.0);
        ld.build(r(0.0, 3.0), 0);
        ld.set_curve_index(1, 9).unwrap();
        assert_eq!(ld.curve_index(1), Some(9));
        assert_eq!(ld.curve_index(0), Some(0));
        assert!(ld.set_curve_index(9, 1).is_err(), "out of range rejected");

        assert!(!ld.is_range_out(0));
        ld.add_out_range(0).unwrap();
        assert!(ld.is_range_out(0));
        assert!(!ld.is_range_out(1));
        ld.add_out_range(0).unwrap(); // idempotent
        assert_eq!(ld.list_range_out(), &[0]);
        ld.remove_range_out_all();
        assert!(!ld.is_range_out(0));
        assert!(ld.add_out_range(99).is_err());
    }

    // ---- SurfaceRangeLocalizeData ----

    #[test]
    fn surface_localize_data_build_grid() {
        let mut ld = SurfaceRangeLocalizeData::new(2, 3, 1e-6, 1e-6);
        ld.set_root_index(5);
        let n = ld.build(r(0.0, 1.0), r(0.0, 3.0), 11);
        assert_eq!(n, 6, "2 × 3 cells");
        assert_eq!(ld.get_nb_samples_u(), 2);
        assert_eq!(ld.get_nb_samples_v(), 3);
        assert_eq!(ld.get_root_index(), 5);
        assert_eq!(ld.ranges().len(), 6);

        // Cell (0, 0): U in [0, 0.5], V in [0, 1].
        let (ur, vr) = ld.get_depth_range(0, 0).unwrap();
        assert_eq!(ur, r(0.0, 0.5));
        assert_eq!(vr, r(0.0, 1.0));
        // Cell (1, 1): U in [0.5, 1], V in [1, 2].
        let (ur, vr) = ld.get_depth_range(1, 1).unwrap();
        assert_eq!(ur, r(0.5, 1.0));
        assert_eq!(vr, r(1.0, 2.0));
        // Last cell (1, 2): V upper bound inclusive.
        let (ur, vr) = ld.get_depth_range(1, 2).unwrap();
        assert_eq!(ur, r(0.5, 1.0));
        assert_eq!(vr, r(2.0, 3.0));

        for iu in 0..2 {
            for iv in 0..3 {
                assert_eq!(ld.surface_index(iu, iv), Some(11));
            }
        }
        assert_eq!(ld.get_depth_range(2, 0), None, "out of grid");

        // Parametric lookup.
        assert_eq!(ld.find_cell(0.25, 0.5), Some((0, 0)));
        assert_eq!(ld.find_cell(0.75, 2.5), Some((1, 2)));
        assert_eq!(ld.find_cell(-1.0, 0.0), None);
    }

    #[test]
    fn surface_localize_data_out_tracking() {
        let mut ld = SurfaceRangeLocalizeData::new(2, 2, 0.0, 0.0);
        ld.build(r(0.0, 1.0), r(0.0, 1.0), 0);
        assert!(!ld.is_range_out(1, 1));
        ld.add_out_range(1, 1).unwrap();
        assert!(ld.is_range_out(1, 1));
        assert!(!ld.is_range_out(0, 1));
        assert!(ld.add_out_range(5, 5).is_err());
        ld.remove_range_out_all();
        assert!(!ld.is_range_out(1, 1));
    }

    // ---- TopolTool ----

    #[test]
    fn topol_tool_plane_nb_samples_and_domain() {
        let mut tool = TopolTool::new();
        assert!(!tool.is_initialized());
        tool.initialize(plane_surface().as_ref());
        assert!(tool.is_initialized());
        // A plane gets a 10 × 10 grid, domain clamped to [-1e5, 1e5]².
        assert_eq!(tool.nb_samples_u(), 10);
        assert_eq!(tool.nb_samples_v(), 10);
        assert_eq!(tool.nb_samples(), 100);
        let (u0, u1) = tool.u_range();
        let (v0, v1) = tool.v_range();
        assert!((u0 + 1e5).abs() < 1e-9 && (u1 - 1e5).abs() < 1e-9);
        assert!((v0 + 1e5).abs() < 1e-9 && (v1 - 1e5).abs() < 1e-9);
        // The step spans the clamped domain.
        assert!((tool.u_step() - (2e5 / 11.0)).abs() < 1e-9);
    }

    #[test]
    fn topol_tool_plane_sample_point_matches_d0() {
        let surf = plane_surface();
        let mut tool = TopolTool::new();
        tool.initialize(surf.as_ref());
        // The standard plane maps (u, v) -> (u, v, 0); every sample's 3D point
        // must agree with surface.d0(u, v) and lie in the plane.
        for i in 1..=tool.nb_samples() {
            let (p2d, p3d) = tool.sample_point(i).unwrap();
            let expected = surf.d0(p2d.x(), p2d.y());
            assert!(p3d.distance(&expected) < 1e-12, "sample {i}: 3D != d0(UV)");
            assert!(p3d.z().abs() < 1e-9, "plane sample {i} off the plane");
            assert!(p2d.x() >= -1e5 && p2d.x() <= 1e5);
            assert!(p2d.y() >= -1e5 && p2d.y() <= 1e5);
        }
    }

    #[test]
    fn topol_tool_plane_sample_point_row_major() {
        let mut tool = TopolTool::new();
        tool.initialize(plane_surface().as_ref());
        let (p1, _) = tool.sample_point(1).unwrap();
        let (p2, _) = tool.sample_point(2).unwrap();
        // Row-major, U fastest: consecutive samples differ by one U step.
        assert!((p2.x() - p1.x()).abs() - tool.u_step() < 1e-9);
        assert!((p2.y() - p1.y()).abs() < 1e-9);
        // Sample 11 is the start of the second V row.
        let (p11, _) = tool.sample_point(11).unwrap();
        assert!((p11.y() - p1.y()).abs() - tool.v_step() < 1e-9);
        assert!((p11.x() - p1.x()).abs() < 1e-9);
    }

    #[test]
    fn topol_tool_cylinder_nb_samples_match_domain() {
        let surf = cylinder_surface(1.0);
        let mut tool = TopolTool::new();
        tool.initialize(surf.as_ref());
        // Radius 1 -> max_angle = 2*acos(1 - 0.01); U count from the 2π span,
        // V count clamped to the 50 cap (2e5/10 = 2e4).
        let max_angle: f64 = 2.0 * (1.0 - 0.01f64).acos();
        let expected_u = (2.0 * PI / max_angle) as usize;
        assert_eq!(tool.nb_samples_u(), expected_u.max(2));
        assert_eq!(tool.nb_samples_v(), 50);
        assert_eq!(tool.nb_samples(), tool.nb_samples_u() * tool.nb_samples_v());
    }

    #[test]
    fn topol_tool_cylinder_sample_point_on_surface() {
        let surf = cylinder_surface(1.0);
        let mut tool = TopolTool::new();
        tool.initialize(surf.as_ref());
        for i in 1..=tool.nb_samples() {
            let (p2d, p3d) = tool.sample_point(i).unwrap();
            let expected = surf.d0(p2d.x(), p2d.y());
            assert!(p3d.distance(&expected) < 1e-12, "sample {i}: 3D != d0(UV)");
            // Cylinder about the Z axis: distance to the axis == radius.
            let radial = (p3d.x() * p3d.x() + p3d.y() * p3d.y()).sqrt();
            assert!((radial - 1.0).abs() < 1e-9, "sample {i} radial {radial} != 1");
        }
    }

    #[test]
    fn topol_tool_sample_point_out_of_range_errors() {
        let mut tool = TopolTool::new();
        assert!(tool.sample_point(1).is_err(), "no surface");
        tool.initialize(plane_surface().as_ref());
        assert!(tool.sample_point(0).is_err());
        assert!(tool.sample_point(101).is_err(), "100 samples, 101 is out");
        assert!(tool.sample_point(100).is_ok());
    }

    #[test]
    fn topol_tool_sample_pnts_plane_uniform_grid_count() {
        let mut tool = TopolTool::new();
        tool.initialize(plane_surface().as_ref());
        let pts = tool.sample_pnts(0.001, 4, 4).unwrap();
        // Plane is not a BSpline: uniform grid of max(10, 4)² = 100 points.
        assert_eq!(pts.len(), 100);
        for (p2d, p3d) in &pts {
            let expected = plane_surface().d0(p2d.x(), p2d.y());
            assert!(p3d.distance(&expected) < 1e-12);
            assert!(p3d.z().abs() < 1e-9);
        }
    }

    #[test]
    fn topol_tool_sample_pnts_plane_floor_applies() {
        let mut tool = TopolTool::new();
        tool.initialize(plane_surface().as_ref());
        let pts = tool.sample_pnts(0.001, 20, 20).unwrap();
        // nu_min/nv_min are a floor: max(10, 20)² = 400.
        assert_eq!(pts.len(), 400);
    }

    #[test]
    fn topol_tool_sample_pnts_cylinder_grid_on_surface() {
        let surf = cylinder_surface(1.0);
        let mut tool = TopolTool::new();
        tool.initialize(surf.as_ref());
        let pts = tool.sample_pnts(0.001, 4, 4).unwrap();
        let expected_len = tool.nb_samples_u() * tool.nb_samples_v();
        assert_eq!(pts.len(), expected_len, "cylinder stays uniform, count from base grid");
        for (p2d, p3d) in &pts {
            let expected = surf.d0(p2d.x(), p2d.y());
            assert!(p3d.distance(&expected) < 1e-12);
            let radial = (p3d.x() * p3d.x() + p3d.y() * p3d.y()).sqrt();
            assert!((radial - 1.0).abs() < 1e-9, "point off cylinder");
        }
    }

    #[test]
    fn topol_tool_sample_pnts_bspline_adaptive_refines() {
        let surf = curved_bspline();
        let mut coarse = TopolTool::new();
        coarse.initialize(surf.as_ref());
        let coarse_pts = coarse.sample_pnts(10.0, 4, 4).unwrap();

        let mut fine = TopolTool::new();
        fine.initialize(surf.as_ref());
        let fine_pts = fine.sample_pnts(1e-4, 4, 4).unwrap();

        // Small deflection must refine the curved (U) direction, producing a
        // strictly larger sample set than the coarse grid.
        assert!(fine_pts.len() > coarse_pts.len());
        assert!(coarse_pts.len() >= 16, "at least the 4×4 floor");
        for (p2d, p3d) in &fine_pts {
            let expected = surf.d0(p2d.x(), p2d.y());
            assert!(p3d.distance(&expected) < 1e-12, "fine sample not on surface");
        }
    }

    #[test]
    fn topol_tool_unit_box_face_sample_points() {
        let box_ = unit_box();
        let face = &box_.faces[0];
        let surf = GeometryRegistry::global()
            .face_surface(&face.0)
            .expect("unit-box face has a registered surface");
        let mut tool = TopolTool::new();
        tool.initialize(surf.as_ref());
        assert_eq!(tool.nb_samples_u(), 10);
        assert_eq!(tool.nb_samples_v(), 10);
        for i in 1..=tool.nb_samples() {
            let (p2d, p3d) = tool.sample_point(i).unwrap();
            let expected = surf.d0(p2d.x(), p2d.y());
            assert!(p3d.distance(&expected) < 1e-12, "box face sample {i}");
        }
        let pts = tool.sample_pnts(0.001, 4, 4).unwrap();
        assert_eq!(pts.len(), 100, "box faces are planes -> uniform 10×10");
    }
