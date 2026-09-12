use super::prelude::*;
use super::*;
    use occt_core::gp::{GpAx2, GpAx3, GpCirc, GpCylinder, GpDir, GpLin, GpPln, GpSphere as GpSphereT};
    use occt_geom::GeomCylinder;
    use crate::primitives::{
        BRepPrimBox, BRepPrimCone, BRepPrimCylinder, BRepPrimSphere, BRepPrimTorus,
    };

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    fn line(p: GpPnt, d: GpDir) -> GeomLine {
        GeomLine::new(GpLin::from_pnt_dir(p, d))
    }

    #[test]
    fn box_area_exact() {
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        let a = analytic_surface_area(&b.solid.0).unwrap();
        assert!(approx(a, 52.0, 1e-9), "area {a}");
    }

    #[test]
    fn box_volume_exact() {
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        let v = analytic_volume(&b.solid.0).unwrap();
        assert!(approx(v, 24.0, 1e-9), "volume {v}");
    }

    #[test]
    fn sphere_area_volume() {
        let s = BRepPrimSphere::make_sphere(2.0);
        let a = analytic_surface_area(&s.solid.0).unwrap();
        assert!(approx(a, 4.0 * PI * 4.0, 1e-6), "area {a}");
        let v = analytic_volume(&s.solid.0).unwrap();
        assert!(approx(v, 4.0 / 3.0 * PI * 8.0, 1e-6), "volume {v}");
    }

    #[test]
    fn cylinder_volume() {
        let c = BRepPrimCylinder::make_cylinder(1.0, 3.0);
        let v = analytic_volume(&c.solid.0).unwrap();
        assert!(approx(v, 3.0 * PI, 1e-6), "volume {v}");
    }

    #[test]
    fn cone_volume() {
        let c = BRepPrimCone::make_cone(2.0, 6.0);
        let v = analytic_volume(&c.solid.0).unwrap();
        assert!(approx(v, PI * 4.0 * 6.0 / 3.0, 1e-6), "volume {v}");
    }

    #[test]
    fn box_centroid() {
        let b = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let c = analytic_centroid(&b.solid.0).unwrap();
        assert!(c.distance(&GpPnt::new(1.0, 1.0, 1.0)) < 1e-9, "centroid {c:?}");
    }

    #[test]
    fn analytic_properties_box() {
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        let p = analytic_properties(&b.solid.0).unwrap();
        assert!(p.exact);
        assert!(approx(p.volume, 24.0, 1e-9));
        assert!(approx(p.surface_area, 52.0, 1e-9));
    }

    #[test]
    fn extrema_point_line_exact() {
        let l = line(GpPnt::zero(), GpDir::new(1.0, 0.0, 0.0).unwrap());
        let e = extrema_point_line(GpPnt::new(3.0, 4.0, 0.0), &l);
        assert!(approx(e.distance, 4.0, 1e-12), "dist {}", e.distance);
        assert!(approx(e.p2.x(), 3.0, 1e-12));
    }

    #[test]
    fn extrema_point_circle_exact() {
        let c = GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0));
        let e = extrema_point_circle(GpPnt::new(3.0, 0.0, 5.0), &c);
        let want = (29.0f64).sqrt();
        assert!(approx(e.distance, want, 1e-9), "dist {} want {want}", e.distance);
        assert!(approx(e.p2.x(), 1.0, 1e-9));
    }

    #[test]
    fn extrema_point_plane() {
        let pl = GeomPlane::new(GpPln::new(GpAx3::standard()));
        let e = super::extrema_point_plane(GpPnt::new(1.0, 2.0, 5.0), &pl);
        assert!(approx(e.distance, 5.0, 1e-12), "dist {}", e.distance);
        assert!(approx(e.p2.z(), 0.0, 1e-12));
    }

    #[test]
    fn extrema_point_sphere() {
        let s = GeomSphere::new(GpSphereT::new(GpAx3::standard(), 1.0).unwrap());
        let e = super::extrema_point_sphere(GpPnt::new(3.0, 0.0, 0.0), &s);
        assert!(approx(e.distance, 2.0, 1e-9), "dist {}", e.distance);
        assert!(approx(e.p2.x(), 1.0, 1e-9));
    }

    #[test]
    fn extrema_line_line_skew() {
        let a = line(GpPnt::zero(), GpDir::new(1.0, 0.0, 0.0).unwrap());
        let b = line(GpPnt::new(0.0, 0.0, 3.0), GpDir::new(0.0, 1.0, 0.0).unwrap());
        let e = extrema_line_line(&a, &b);
        assert!(approx(e.distance, 3.0, 1e-9), "dist {}", e.distance);
    }

    #[test]
    fn extrema_line_sphere() {
        let s = GeomSphere::new(GpSphereT::new(GpAx3::standard(), 1.0).unwrap());
        // Line at y=3 → distance 2.
        let l = line(GpPnt::new(0.0, 3.0, 0.0), GpDir::new(1.0, 0.0, 0.0).unwrap());
        let e = super::extrema_line_sphere(&l, &s).expect("extrema");
        assert!(approx(e.distance, 2.0, 1e-9), "dist {}", e.distance);
        // Line through the center → distance 0.
        let l2 = line(GpPnt::zero(), GpDir::new(1.0, 0.0, 0.0).unwrap());
        let e2 = super::extrema_line_sphere(&l2, &s).expect("extrema");
        assert!(approx(e2.distance, 0.0, 1e-9), "dist {}", e2.distance);
    }

    #[test]
    fn extrema_point_cylinder() {
        let cyl = GeomCylinder::new(GpCylinder::new(GpAx3::standard(), 1.0).unwrap());
        let e = extrema_point_surface_exact(&cyl, GpPnt::new(3.0, 0.0, 4.0)).unwrap();
        assert!(approx(e.distance, 2.0, 1e-9), "dist {}", e.distance);
        assert!(approx(e.p2.z(), 4.0, 1e-9));
    }

    #[test]
    fn is_analytic_box_sphere() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        assert!(is_analytic(&b.solid.0));
        let s = BRepPrimSphere::make_sphere(1.0);
        assert!(is_analytic(&s.solid.0));
    }

    #[test]
    fn torus_area_known() {
        let t = BRepPrimTorus::make_torus(1.0, 1.0);
        let a = analytic_surface_area(&t.solid.0).unwrap();
        assert!(approx(a, 4.0 * PI * PI * 1.0 * 1.0, 1e-3), "area {a}");
    }

    #[test]
    fn extrema_point_surface_sphere() {
        let s = GeomSphere::new(GpSphereT::new(GpAx3::standard(), 1.0).unwrap());
        let e = extrema_point_surface_exact(&s, GpPnt::new(0.0, 3.0, 0.0)).unwrap();
        assert!(approx(e.distance, 2.0, 1e-9), "dist {}", e.distance);
    }

    #[test]
    fn torus_volume_known() {
        let t = BRepPrimTorus::make_torus(2.0, 1.0);
        let v = analytic_volume(&t.solid.0).unwrap();
        assert!(approx(v, 2.0 * PI * PI * 2.0 * 1.0, 1e-6), "volume {v}");
    }

    #[test]
    fn box_inertia_diagonal() {
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0); // dx=2, dy=3, dz=4
        let t = inertia_tensor(&b.solid.0, 1.0).unwrap();
        let m = 24.0;
        // Ixx about the x-axis = m/12·(dy²+dz²) = 2·(9+16) = 50.
        assert!(approx(t.ixx, m / 12.0 * (3.0 * 3.0 + 4.0 * 4.0), 1e-9), "ixx {}", t.ixx);
        assert!(approx(t.iyy, m / 12.0 * (2.0 * 2.0 + 4.0 * 4.0), 1e-9), "iyy {}", t.iyy);
        assert!(approx(t.izz, m / 12.0 * (2.0 * 2.0 + 3.0 * 3.0), 1e-9), "izz {}", t.izz);
        assert!(t.ixy.abs() < 1e-12 && t.ixz.abs() < 1e-12 && t.iyz.abs() < 1e-12);
        // inertia_matrix mirrors the tensor.
        let mat = inertia_matrix(&b.solid.0, 1.0).unwrap();
        assert!(approx(mat[0][0], t.ixx, 1e-12) && approx(mat[1][1], t.iyy, 1e-12) && approx(mat[2][2], t.izz, 1e-12));
        assert!(approx(mat[0][1], t.ixy, 1e-12) && approx(mat[0][2], t.ixz, 1e-12) && approx(mat[1][2], t.iyz, 1e-12));
    }

    #[test]
    fn sphere_inertia() {
        let s = BRepPrimSphere::make_sphere(2.0);
        let t = inertia_tensor(&s.solid.0, 1.0).unwrap();
        let m = 4.0 / 3.0 * PI * 8.0;
        let want = (2.0 / 5.0) * m * 4.0;
        assert!(approx(t.ixx, want, 1e-6), "ixx {} want {}", t.ixx, want);
        assert!(approx(t.iyy, want, 1e-6));
        assert!(approx(t.izz, want, 1e-6));
        assert!(t.ixy.abs() < 1e-9);
    }

    #[test]
    fn cylinder_inertia() {
        let c = BRepPrimCylinder::make_cylinder(1.0, 3.0);
        let t = inertia_tensor(&c.solid.0, 1.0).unwrap();
        let m = PI * 3.0;
        assert!(approx(t.izz, 0.5 * m * 1.0, 1e-6), "izz {}", t.izz);
        let ix_want = m * (3.0 * 1.0 + 9.0) / 12.0;
        assert!(approx(t.ixx, ix_want, 1e-6), "ixx {} want {}", t.ixx, ix_want);
        assert!(approx(t.iyy, ix_want, 1e-6), "iyy {}", t.iyy);
        assert!(t.ixy.abs() < 1e-9 && t.ixz.abs() < 1e-9 && t.iyz.abs() < 1e-9);
    }

    #[test]
    fn principal_inertia_box() {
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        let (moments, _) = principal_inertia(&b.solid.0, 1.0).unwrap();
        assert_eq!(moments.len(), 3);
        // Sorted descending: largest moment about the smallest cross-section
        // axis (the 2 side → Ixx = 50), smallest about the largest (4 → Izz = 26).
        assert!(moments[0] >= moments[1] && moments[1] >= moments[2]);
        assert!(approx(moments[0], 50.0, 1e-6), "largest {}", moments[0]);
        assert!(approx(moments[1], 40.0, 1e-6), "middle {}", moments[1]);
        assert!(approx(moments[2], 26.0, 1e-6), "smallest {}", moments[2]);
    }

    #[test]
    fn principal_axes_orthogonal() {
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        let (_, axes) = principal_inertia(&b.solid.0, 1.0).unwrap();
        assert_eq!(axes.len(), 3);
        for i in 0..3 {
            assert!((axes[i].magnitude() - 1.0).abs() < 1e-6, "axis {i} not unit");
            for j in (i + 1)..3 {
                assert!(axes[i].dot(&axes[j]).abs() < 1e-6, "dot {i},{j} = {}", axes[i].dot(&axes[j]));
            }
        }
    }

    /// Build a properly-wired axis-aligned box solid with one corner at `p0`
    /// and extents `dx × dy × dz` (mirrors `BRepPrimBox` but at an arbitrary
    /// origin), so its centroid is not the global origin.
    fn box_solid_at(p0: &GpPnt, dx: f64, dy: f64, dz: f64) -> crate::shape::Solid {
        use std::sync::Arc;
        let b = crate::builder::TopoBuilder::new();
        let corners = [
            GpPnt::new(p0.x(), p0.y(), p0.z()),
            GpPnt::new(p0.x() + dx, p0.y(), p0.z()),
            GpPnt::new(p0.x() + dx, p0.y() + dy, p0.z()),
            GpPnt::new(p0.x(), p0.y() + dy, p0.z()),
            GpPnt::new(p0.x(), p0.y(), p0.z() + dz),
            GpPnt::new(p0.x() + dx, p0.y(), p0.z() + dz),
            GpPnt::new(p0.x() + dx, p0.y() + dy, p0.z() + dz),
            GpPnt::new(p0.x(), p0.y() + dy, p0.z() + dz),
        ];
        let verts: Vec<crate::shape::Vertex> =
            corners.iter().map(|p| b.make_vertex(*p, 0.0)).collect();
        let box_edges: [(usize, usize); 12] = [
            (0, 1), (1, 2), (2, 3), (3, 0),
            (4, 5), (5, 6), (6, 7), (7, 4),
            (0, 4), (1, 5), (2, 6), (3, 7),
        ];
        let mut edges: Vec<crate::shape::Edge> = Vec::new();
        for &(i, j) in &box_edges {
            let p1 = corners[i];
            let p2 = corners[j];
            let dir = GpDir::from_vec(&GpVec::from_pnts(&p1, &p2)).expect("box edge direction");
            let mut e = b.make_edge(
                Arc::new(GeomLine::new(GpLin::from_pnt_dir(p1, dir))),
                0.0,
                p1.distance(&p2),
            );
            b.add(&mut e.0, &verts[i].0);
            b.add(&mut e.0, &verts[j].0);
            edges.push(e);
        }
        let edge_index = |a: usize, b: usize| {
            box_edges
                .iter()
                .position(|&(i, j)| (i == a && j == b) || (i == b && j == a))
                .expect("box edge")
        };
        let plane = |origin: GpPnt, normal: GpDir| {
            let z_axis = GpDir::new(0.0, 0.0, 1.0).unwrap();
            let xd = if normal.is_normal(&z_axis) { z_axis } else { GpDir::new(1.0, 0.0, 0.0).unwrap() };
            GpPln::new(GpAx3::new(origin, normal, &xd).unwrap())
        };
        let faces_def: [(GpPnt, GpDir, [usize; 4]); 6] = [
            (corners[0], GpDir::new(0.0, 0.0, -1.0).unwrap(), [0, 3, 2, 1]),
            (corners[4], GpDir::new(0.0, 0.0, 1.0).unwrap(), [4, 5, 6, 7]),
            (corners[0], GpDir::new(0.0, -1.0, 0.0).unwrap(), [0, 1, 5, 4]),
            (corners[3], GpDir::new(0.0, 1.0, 0.0).unwrap(), [3, 2, 6, 7]),
            (corners[0], GpDir::new(-1.0, 0.0, 0.0).unwrap(), [0, 4, 7, 3]),
            (corners[1], GpDir::new(1.0, 0.0, 0.0).unwrap(), [1, 2, 6, 5]),
        ];
        let mut faces = Vec::new();
        for (origin, normal, cycle) in faces_def {
            let quad: Vec<crate::shape::Edge> = [0, 1, 2, 3]
                .iter()
                .map(|&k| edges[edge_index(cycle[k], cycle[(k + 1) % 4])].clone())
                .collect();
            let wire = b.make_wire(&quad);
            let surface: Arc<dyn Surface> = Arc::new(GeomPlane::new(plane(origin, normal)));
            faces.push(b.make_face(surface, &[wire]));
        }
        let shell = b.make_shell(&faces);
        b.make_solid(&[shell])
    }

    #[test]
    fn inertia_about_centroid() {
        // A box spanning [0,dx]×[0,dy]×[0,dz] and the same box translated: the
        // tensor about each box's own centroid is identical (the analytic
        // formulas are inherently about the centroid, so no origin offset leaks
        // in and no parallel-axis term needs to be subtracted).
        let a = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        let b = box_solid_at(&GpPnt::new(1.0, -2.0, 5.0), 2.0, 3.0, 4.0);
        let ta = inertia_tensor(&a.solid.0, 1.0).unwrap();
        let tb = inertia_tensor(&b.0, 1.0).unwrap();
        assert!(approx(ta.ixx, tb.ixx, 1e-9), "ixx {} vs {}", ta.ixx, tb.ixx);
        assert!(approx(ta.iyy, tb.iyy, 1e-9));
        assert!(approx(ta.izz, tb.izz, 1e-9));
        assert!(ta.ixy.abs() < 1e-9 && tb.ixy.abs() < 1e-9);
        // The value is the centroid tensor (Ixx = m/12(dy²+dz²) = 50). Had the
        // tensor been about the global origin the Ixx would be 200.
        assert!(approx(ta.ixx, 50.0, 1e-9), "centroid Ixx {}", ta.ixx);
    }

    #[test]
    fn analytic_curve_length_box() {
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        let l = analytic_curve_length(&b.solid.0).unwrap();
        assert!(approx(l, 4.0 * (2.0 + 3.0 + 4.0), 1e-9), "length {l}");
    }

    #[test]
    fn circle_length() {
        // A cylinder: two full cap circles (r=1) + one seam line (h=1).
        let c = BRepPrimCylinder::make_cylinder(1.0, 1.0);
        let l = analytic_curve_length(&c.solid.0).unwrap();
        let want = 2.0 * 2.0 * PI * 1.0 + 1.0;
        assert!(approx(l, want, 1e-6), "length {l} want {want}");

        // Half circle edge: params [0, π] on a unit circle → π.
        let b = crate::builder::TopoBuilder::new();
        let half = b.make_edge_circle(&GpAx2::standard(), 1.0, 0.0, PI);
        let len = exact_edge_length(&half);
        assert!(approx(len, PI, 1e-9), "half circle {len}");
    }

    #[test]
    fn line_length() {
        let b = crate::builder::TopoBuilder::new();
        let e = b.make_edge_segment(&GpPnt::zero(), &GpPnt::new(3.0, 4.0, 0.0));
        let len = exact_edge_length(&e);
        assert!(approx(len, 5.0, 1e-9), "line length {len}");
    }

    #[test]
    fn full_analytic_props_box() {
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        let p = full_analytic_properties(&b.solid.0, 1.0).unwrap();
        assert!(approx(p.mass, 24.0, 1e-9), "mass {}", p.mass);
        assert!(approx(p.volume, 24.0, 1e-9));
        assert!(approx(p.surface_area, 52.0, 1e-9));
        assert!(p.exact);
        assert_eq!(p.principal_moments.len(), 3);
        assert_eq!(p.principal_axes.len(), 3);
        assert!(approx(p.inertia.ixx, 50.0, 1e-9));
        assert!(approx(p.centroid.x(), 1.0, 1e-9));
        assert!(approx(p.centroid.y(), 1.5, 1e-9));
        assert!(approx(p.centroid.z(), 2.0, 1e-9));
    }

    #[test]
    fn mass_density() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let m = analytic_mass(&b.solid.0, 2.0).unwrap();
        assert!(approx(m, 2.0, 1e-9), "mass {m}");
    }

    #[test]
    fn moment_about_axis_box() {
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        let t = inertia_tensor(&b.solid.0, 1.0).unwrap();
        // Moment about each coordinate axis is the corresponding diagonal entry.
        assert!(approx(t.moment_about_axis(&GpVec::new(1.0, 0.0, 0.0)), 50.0, 1e-9));
        assert!(approx(t.moment_about_axis(&GpVec::new(0.0, 1.0, 0.0)), 40.0, 1e-9));
        assert!(approx(t.moment_about_axis(&GpVec::new(0.0, 0.0, 1.0)), 26.0, 1e-9));
        // The axis (1,1,1)/√3 is equally weighted: I = (Ixx+Iyy+Izz)/3.
        let u = GpVec::new(1.0, 1.0, 1.0).normalized();
        assert!(approx(t.moment_about_axis(&u), (50.0 + 40.0 + 26.0) / 3.0, 1e-9));
        // Inertia matrix applied to a unit vector recovers the same value.
        let v = t.apply(&u);
        assert!(approx(u.dot(&v), (50.0 + 40.0 + 26.0) / 3.0, 1e-9));
    }

    #[test]
    fn inertia_tensor_at_origin() {
        // The tensor about the centroid, shifted to the origin by the
        // parallel-axis theorem, matches the closed-form origin tensor.
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        let origin = GpPnt::zero();
        let t_origin = inertia_tensor_at(&b.solid.0, 1.0, &origin).unwrap();
        // Box 2×3×4 spanning [0,2]×[0,3]×[0,4], centroid (1,1.5,2), mass 24.
        // Ixx_origin = Ixx_cm + m(d²−dx²) = 50 + 24·(7.25−1) = 200.
        assert!(approx(t_origin.ixx, 200.0, 1e-9), "ixx {}", t_origin.ixx);
        assert!(approx(t_origin.iyy, 40.0 + 24.0 * (7.25 - 2.25), 1e-9));
        assert!(approx(t_origin.izz, 26.0 + 24.0 * (7.25 - 4.0), 1e-9));
        // Subtracting the parallel-axis term (negative mass) returns to the
        // centroid tensor.
        let back = t_origin.translated(-24.0, &GpVec::new(-1.0, -1.5, -2.0));
        assert!(approx(back.ixx, 50.0, 1e-9), "round-trip ixx {}", back.ixx);
    }

    #[test]
    fn inertia_composite_boxes() {
        // Two unit boxes side by side along x. Combined centroid at (1,0.5,0.5).
        let a = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let b = box_solid_at(&GpPnt::new(1.0, 0.0, 0.0), 1.0, 1.0, 1.0);
        let (t, c) = inertia_tensor_composite(&[&a.solid.0, &b.0], 1.0).unwrap();
        assert!(approx(c.x(), 1.0, 1e-9) && approx(c.y(), 0.5, 1e-9) && approx(c.z(), 0.5, 1e-9));
        // Each 1×1×1 box about its own centroid: Ixx=Iyy=Izz=1/6, m=1.
        // Shifted to the combined centroid: Ixx = 1/3, Iyy = Izz = 5/6.
        assert!(approx(t.ixx, 1.0 / 3.0, 1e-9), "ixx {}", t.ixx);
        assert!(approx(t.iyy, 5.0 / 6.0, 1e-9), "iyy {}", t.iyy);
        assert!(approx(t.izz, 5.0 / 6.0, 1e-9), "izz {}", t.izz);
        assert!(t.ixy.abs() < 1e-9);
    }

    #[test]
    fn curve_classify_line_circle() {
        let b = crate::builder::TopoBuilder::new();
        let line_edge = b.make_edge_segment(&GpPnt::zero(), &GpPnt::new(5.0, 0.0, 0.0));
        let curve = BRepTool::edge_curve(&line_edge).unwrap();
        let (a, bb) = BRepTool::edge_parameters(&line_edge);
        assert_eq!(classify_curve(curve.as_ref(), a, bb), CurveKind::Line);

        let circle_edge = b.make_edge_circle(&GpAx2::standard(), 2.0, 0.0, PI);
        let curve = BRepTool::edge_curve(&circle_edge).unwrap();
        let (a, bb) = BRepTool::edge_parameters(&circle_edge);
        assert_eq!(classify_curve(curve.as_ref(), a, bb), CurveKind::Circle);
        assert!(approx(edge_length_exact(&circle_edge), 2.0 * PI, 1e-9));
    }

    #[test]
    fn ellipse_ramanujan_circle_limit() {
        // A circle is an ellipse with a == b: the Ramanujan perimeter gives 2πa.
        assert!(approx(ellipse_length_ramanujan(2.0, 2.0), 4.0 * PI, 1e-9));
        // Degenerate flat ellipse a=b? Rather, an extreme a=10, b=1 gives a
        // value between 4a and 2π·(a+b)/... just check it is sane.
        let p = ellipse_length_ramanujan(10.0, 1.0);
        assert!(p > 40.0 && p < 44.0, "ellipse perimeter {p}");
    }

    #[test]
    fn classify_solid_analytics() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        assert!(shape_is_box(&b.solid.0));
        assert_eq!(classify_solid(&b.solid.0), SolidKind::Box);
        let s = BRepPrimSphere::make_sphere(1.0);
        assert_eq!(classify_solid(&s.solid.0), SolidKind::Sphere);
        let c = BRepPrimCylinder::make_cylinder(1.0, 1.0);
        assert_eq!(classify_solid(&c.solid.0), SolidKind::Cylinder);
        let k = BRepPrimCone::make_cone(1.0, 1.0);
        assert_eq!(classify_solid(&k.solid.0), SolidKind::Cone);
        let t = BRepPrimTorus::make_torus(1.0, 1.0);
        assert_eq!(classify_solid(&t.solid.0), SolidKind::Torus);
    }

    #[test]
    fn full_props_inertia_about() {
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        let p = full_analytic_properties(&b.solid.0, 1.0).unwrap();
        // inertia_about(centroid) is the stored centroid tensor.
        let at_c = p.inertia_about(&p.centroid);
        assert!(approx(at_c.ixx, p.inertia.ixx, 1e-12));
        // Radius of gyration about the x principal axis: sqrt(50/24).
        let rg = p.radius_of_gyration();
        assert_eq!(rg.len(), 3);
        assert!(approx(rg[0], (50.0f64 / 24.0f64).sqrt(), 1e-9), "rg {}", rg[0]);
    }
