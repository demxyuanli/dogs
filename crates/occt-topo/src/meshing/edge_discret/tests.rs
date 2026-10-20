//! Unit tests of edge discretization.

use super::*;

    use occt_core::gp::{GpAx2, GpCirc, GpDir, GpLin};
    use occt_geom::{GeomCircle, GeomLine, GeomTrimmedCurve};

    fn segment_curve(len: f64) -> Arc<dyn Curve> {
        let lin = GpLin::from_pnt_dir(GpPnt::zero(), GpDir::new(1.0, 0.0, 0.0).unwrap());
        Arc::new(GeomTrimmedCurve::new(Arc::new(GeomLine::new(lin)), 0.0, len))
    }

    fn arc_curve(radius: f64, from: f64, to: f64) -> Arc<dyn Curve> {
        let circ = GeomCircle::new(GpCirc::new(GpAx2::standard(), radius));
        Arc::new(GeomTrimmedCurve::new(Arc::new(circ), from, to))
    }

    #[test]
    fn straight_edge_endpoints_match() {
        let curve = segment_curve(2.0);
        let pts = EdgeDiscret::discretize_edge(&curve, 0.0, 1.0, 0.05);
        assert!(pts.len() >= 2, "polyline length {}", pts.len());
        assert!(
            pts.first().unwrap().distance(&GpPnt::new(0.0, 0.0, 0.0)) < 1e-9,
            "first endpoint {:?}",
            pts.first()
        );
        assert!(
            pts.last().unwrap().distance(&GpPnt::new(2.0, 0.0, 0.0)) < 1e-9,
            "last endpoint {:?}",
            pts.last()
        );
        // A straight chord deviates nothing: no refinement past the seed.
        assert!(pts.len() <= 64, "unexpected refinement for a line: {}", pts.len());
    }

    #[test]
    fn circle_refines_to_chord_deviation_target() {
        let curve = arc_curve(2.0, 0.0, std::f64::consts::PI);
        let deflection = 0.05;
        let tess = CurveTessellator::from_range(curve, 0.0, 1.0, deflection, 2);
        assert!(tess.points_nb() >= 2);
        // Endpoints: u=0 → (2,0,0), u=1 → (-2,0,0).
        assert!(tess.point(0).unwrap().distance(&GpPnt::new(2.0, 0.0, 0.0)) < 1e-6);
        assert!(tess.point(tess.points_nb() - 1).unwrap().distance(&GpPnt::new(-2.0, 0.0, 0.0)) < 1e-6);
        // The chord deviation of the emitted polyline must meet the target.
        let dev = tess.max_chord_error(16);
        assert!(dev < deflection, "chord deviation {dev} >= deflection {deflection}");
    }

    #[test]
    fn angular_deflection_refines_arc() {
        let curve = arc_curve(1.0, 0.0, std::f64::consts::PI);
        // Loose linear deflection (chord deviation only) vs. tight angular bound:
        // the angular term must add points even when the chord is already short.
        let lin = CurveTessellator::from_range(curve.clone(), 0.0, 1.0, 10.0, 2);
        let ang = CurveTessellator::from_range_angular(curve.clone(), 0.0, 1.0, 10.0, 0.05, 2);
        assert!(
            ang.points_nb() > lin.points_nb(),
            "angular {} should exceed linear {}",
            ang.points_nb(),
            lin.points_nb()
        );
    }

    #[test]
    fn min_points_respected_and_deflection_drives_density() {
        let curve = arc_curve(1.0, 0.0, std::f64::consts::PI);
        let coarse = CurveTessellator::from_range(curve.clone(), 0.0, 1.0, 0.2, 4);
        let fine = CurveTessellator::from_range(curve.clone(), 0.0, 1.0, 0.01, 4);
        assert!(coarse.points_nb() >= 4, "min points violated: {}", coarse.points_nb());
        assert!(fine.points_nb() >= coarse.points_nb(), "finer deflection must refine more");
        assert!(fine.max_chord_error(16) < 0.01);
    }

    #[test]
    fn parameter_provider_remaps_stored_range() {
        let p = EdgeParameterProvider::with_stored(1.0, 3.0, 0.0, 2.0);
        assert!(!p.is_same_param());
        assert!((p.scale() - 1.0).abs() < 1e-12);
        // Stored 1.0 (midpoint of [0,2]) maps to 2.0 (midpoint of [1,3]).
        assert!((p.remap(1.0) - 2.0).abs() < 1e-12);
        // Uniform grid is inclusive of both endpoints.
        let u = p.uniform_parameters(3);
        assert_eq!(u, vec![1.0, 2.0, 3.0]);
        // Same-parameter provider passes values through.
        let same = EdgeParameterProvider::new(0.0, 4.0);
        assert!(same.is_same_param());
        assert_eq!(same.remap(1.5), 1.5);
    }

    #[test]
    fn tessellate_3d_replaces_endpoints_with_vertices() {
        let edge = DiscreteMeshEdge::new(segment_curve(1.0), 0.0, 1.0);
        let t = CurveTessellator::from_range(edge.curve.clone(), 0.0, 1.0, 0.1, 4);
        let a = GpPnt::new(0.0, 0.0, 0.0);
        let b = GpPnt::new(1.0, 0.0, 0.0);
        let pts = EdgeDiscret::tessellate_3d(&edge, &t, a, b);
        // Endpoints are the vertex points (BRep_Tool::Pnt), not curve values.
        assert_eq!(pts.first(), Some(&a));
        assert_eq!(pts.last(), Some(&b));
        // Interior points come from the tessellator; total ≤ tessellator count.
        assert!(pts.len() <= t.points_nb(), "len {}", pts.len());
        assert!(pts.len() >= 2);
    }

    #[test]
    fn deflection_absolute_and_consistency() {
        // Relative deflection against a 2-unit shape: clamped coefficient ≥ 0.5.
        let abs = EdgeDeflection::compute_absolute_deflection(0.01, 2.0);
        assert!(abs > 0.0 && abs <= 0.02);
        assert_eq!(EdgeDeflection::compute_absolute_deflection(0.01, 0.0), 0.01);
        assert!(EdgeDeflection::is_consistent(0.05, 0.05, false, 0.1));
        assert!(!EdgeDeflection::is_consistent(0.5, 0.05, false, 0.1));
        // Quality decrease allowed: too-small current is inconsistent.
        assert!(!EdgeDeflection::is_consistent(0.001, 0.05, true, 0.1));
    }
