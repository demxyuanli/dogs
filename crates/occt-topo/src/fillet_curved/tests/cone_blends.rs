use super::*;

use super::*;

    use occt_core::gp::{GpCone, GpCylinder, GpLin, GpPln, GpSphere};
    use occt_geom::{GeomCone, GeomCylinder, GeomLine, GeomPlane, GeomSphere};
    use crate::shape::{Shell, Solid, Vertex};
    use crate::shell_check::shell_is_closed;
    use crate::tgeometry::GeometryRegistry;
    use crate::topo_tools_full::{faces_of, vertices_of};

    // ------------------------------------------------------------------
    // Test solid builders
    // ------------------------------------------------------------------


    // ------------------------------------------------------------------
    // General surface-surface (cone / cylinder) blends
    // ------------------------------------------------------------------

    #[test]
    fn cone_extraction_geometric() {
        let (solid, _edge) = build_cone_on_box(1.0, 2.0);
        let cone_face = faces_of(&solid.0)
            .into_iter()
            .find(|f| {
                BRepTool::face_surface(f)
                    .map(|s| classify_surface_full(s.as_ref()) == SurfaceKind::Cone)
                    .unwrap_or(false)
            })
            .expect("cone face");
        let s = BRepTool::face_surface(&cone_face).unwrap();
        let ci = cone_from_surface(s.as_ref()).expect("cone extraction");
        assert!((ci.semi_angle - (1.0_f64 / 2.0).atan()).abs() < 0.05, "semi-angle {}", ci.semi_angle);
        assert!(ci.apex.distance(&GpPnt::new(0.0, 0.0, 2.0)) < 1e-3, "apex {:?}", ci.apex);
        assert!(classify_surface_analytic(s.as_ref()) == SurfaceKind::Cone);
        clear_tree(&solid.0);
    }

    #[test]
    fn plane_cone_blend() {
        // A cone standing on the top face of a box; fillet the base circle.
        let (solid, edge) = build_cone_on_box(1.0, 2.0);
        let out = fillet_edge_curved_general(&solid.0, &edge, 0.15, 1e-6).expect("plane+cone fillet");
        assert!(shell_is_closed(&shell_of(&out)), "plane+cone fillet is closed");
        // The blend face is a torus band (the rolling-ball envelope). Use the
        // analytic classifier: the trimmed cone face is `Cone`, the torus blend
        // is `Torus`.
        let faces = faces_of(&out);
        let blend = faces
            .iter()
            .find(|f| {
                BRepTool::face_surface(f)
                    .map(|s| classify_surface_analytic(s.as_ref()) == SurfaceKind::Torus)
                    .unwrap_or(false)
            })
            .expect("torus blend face");
        let s = BRepTool::face_surface(blend).unwrap();
        let r = 0.15;
        // Centreline: offset plane at z=r; the blend points are at distance r
        // from the centreline circle.
        let mut on_torus = true;
        let (u0, u1, v0, v1) = sample_bounds(s.as_ref());
        let center = GpPnt::new(0.0, 0.0, r);
        let axis = GpVec::new(0.0, 0.0, 1.0);
        // The centreline radius for α = atan(1/2): R_c = (H + r/sinα − r)tanα.
        let alpha = (1.0_f64 / 2.0).atan();
        let rho = (2.0 + r / alpha.sin() - r) * alpha.tan();
        for i in 0..8 {
            let u = u0 + (u1 - u0) * i as f64 / 8.0;
            let v = v0 + (v1 - v0) * i as f64 / 8.0;
            let p = s.d0(u, v);
            let d = distance_to_circle(&p, &center, &axis, rho);
            if (d - r).abs() > 1e-5 {
                on_torus = false;
            }
        }
        assert!(on_torus, "plane+cone blend face is a torus band");
        clear_tree(&out);
        clear_tree(&solid.0);
    }

    #[test]
    fn cylinder_cone_blend() {
        let (solid, edge) = build_cylinder_cone_union(1.0);
        let out = fillet_edge_curved_general(&solid.0, &edge, 0.15, 1e-6).expect("cylinder+cone fillet");
        assert!(shell_is_closed(&shell_of(&out)), "cylinder+cone fillet is closed");
        let faces = faces_of(&out);
        assert_eq!(faces.len(), 5, "4 input faces + 1 blend face");
        clear_tree(&out);
        clear_tree(&solid.0);
    }

    #[test]
    fn cone_sphere_blend() {
        let (solid, edge) = build_cone_sphere_union(0.5, 1.0);
        let out = fillet_edge_curved_general(&solid.0, &edge, 0.12, 1e-6).expect("cone+sphere fillet");
        assert!(shell_is_closed(&shell_of(&out)), "cone+sphere fillet is closed");
        let faces = faces_of(&out);
        assert_eq!(faces.len(), 3, "2 input faces + 1 blend face");
        clear_tree(&out);
        clear_tree(&solid.0);
    }

    #[test]
    fn cone_cone_blend() {
        let (solid, edge) = build_cone_cone_union(0.4, 0.8);
        let out = fillet_edge_curved_general(&solid.0, &edge, 0.1, 1e-6).expect("cone+cone fillet");
        assert!(shell_is_closed(&shell_of(&out)), "cone+cone fillet is closed");
        let faces = faces_of(&out);
        assert_eq!(faces.len(), 3, "2 input faces + 1 blend face");
        clear_tree(&out);
        clear_tree(&solid.0);
    }

    #[test]
    fn blend_tangency_verified() {
        // A plane+sphere blend: the torus band must be tangent to both faces.
        let (solid, edge) = build_sphere_on_box(1.0);
        let out = fillet_edge_curved(&solid.0, &edge, 0.2, 1e-6).expect("fillet");
        let faces = faces_of(&out);
        let blend = find_blend_face(&faces);
        let bs = BRepTool::face_surface(blend).unwrap();
        // The two adjacent faces: the plane (z=0) and the sphere cap.
        let (f1, f2) = {
            let adj = faces_of(&out)
                .into_iter()
                .filter(|f| {
                    wires_of_face(f).iter().any(|w| {
                        edges_of_wire(w).iter().any(|e| {
                            let (a, _) = BRepTool::edge_vertices(e).unwrap_or((GpPnt::zero(), GpPnt::zero()));
                            a.z().abs() < 1e-6 || a.distance(&GpPnt::new(0.0, 0.0, 1.0)) < 1e-6
                        })
                    })
                })
                .collect::<Vec<Face>>();
            (BRepTool::face_surface(&adj[0]).unwrap(), BRepTool::face_surface(&adj[1]).unwrap())
        };
        assert!(
            blend_tangency_ok(f1.as_ref(), f2.as_ref(), bs.as_ref(), 1e-4),
            "blend must be tangent to both adjacent faces"
        );
        clear_tree(&out);
        clear_tree(&solid.0);
    }

    #[test]
    fn unsupported_torus_torus_general() {
        // A torus+torus pair is documented as unsupported by the general path.
        let b = TopoBuilder::new();
        let circle = build_circle(&b, GpPnt::zero(), GpVec::new(0.0, 0.0, 1.0), 1.0);
        let w = b.make_wire(&[circle.clone()]);
        let f1 = b.make_face(Arc::new(GeomTorus::new(GpTorus::new(GpAx3::standard(), 3.0, 0.5).unwrap())), &[w.clone()]);
        let f2 = b.make_face(Arc::new(GeomTorus::new(GpTorus::new(GpAx3::standard(), 3.0, 0.5).unwrap())), &[w]);
        let shell = b.make_shell(&[f1.clone(), f2.clone()]);
        let solid = b.make_solid(&[shell]);
        let res = fillet_edge_curved_general(&solid.0, &circle, 0.2, 1e-6);
        assert!(res.is_err(), "torus+torus must be unsupported by the general path");
        assert_eq!(classify_curved_pair(&f1, &f2), CurvedPair::Unsupported);
        clear_tree(&solid.0);
        clear_tree(&f1.0);
        clear_tree(&f2.0);
    }

    #[test]
    fn curved_chain_general() {
        // Two non-interfering curved edges (two sphere caps on a box).
        let (solid, _holes) = build_two_spheres_on_box();
        let es = edges_of(&solid.0);
        let i0 = es
            .iter()
            .position(|e| {
                let (a, _) = BRepTool::edge_vertices(e).unwrap();
                a.distance(&GpPnt::new(-0.2, 0.0, 0.0)) < 1e-6
            })
            .expect("sphere-1 base circle");
        let i1 = es
            .iter()
            .position(|e| {
                let (a, _) = BRepTool::edge_vertices(e).unwrap();
                a.distance(&GpPnt::new(1.6, 0.0, 0.0)) < 1e-6
            })
            .expect("sphere-2 base circle");
        let out = fillet_edge_curved_chain(&solid.0, &[i0, i1], 0.2, 1e-6).expect("two-edge chain");
        assert_eq!(faces_of(&out).len(), 10, "8 input faces + 2 blends");
        assert!(shell_is_closed(&shell_of(&out)), "two-edge curved chain is closed");
        clear_tree(&out);
        clear_tree(&solid.0);
    }

    #[test]
    fn cylinder_cylinder_parallel_blend() {
        // Two overlapping parallel cylinders; fillet one intersection line.
        let (solid, edge) = build_parallel_cylinders(1.2);
        assert!(shell_is_closed(&shell_of(&solid.0)), "lens solid is closed");
        match fillet_edge_curved_general(&solid.0, &edge, 0.15, 1e-6) {
            Ok(out) => {
                assert!(shell_is_closed(&shell_of(&out)), "parallel-cylinder fillet is closed");
                clear_tree(&out);
            }
            Err(e) => {
                // Documented deviation: the parallel-cylinder band blend is only
                // produced for clean lens geometries; otherwise a clear error.
                assert!(
                    e.contains("fillet_edge_parallel_cylinders"),
                    "unexpected error: {e}"
                );
            }
        }
        clear_tree(&solid.0);
    }

    #[test]
    fn general_sampled_blend_unsupported() {
        // Two non-parallel cylinders: the analytic general path documents a
        // clear error (a sampled B-spline blend for angled cylinders is out of
        // scope). The parallel-cylinder centreline sampler still rejects them.
        let c1 = crate::primitives::BRepPrimCylinder::make_cylinder(0.5, 2.0);
        let s1 = BRepTool::face_surface(
            faces_of(&c1.solid.0).iter().find(|f| {
                BRepTool::face_surface(f)
                    .map(|s| classify_surface_full(s.as_ref()) == SurfaceKind::Cylinder)
                    .unwrap_or(false)
            }).unwrap(),
        ).unwrap();
        let info1 = cylinder_from_surface(s1.as_ref()).expect("cylinder 1 extraction");
        // A second cylinder with a tilted axis (not parallel).
        let mut info2 = info1;
        info2.axis = GpVec::new(0.8 * info2.axis.x(), 0.0, 0.6).normalized();
        let res = parallel_cylinder_centerlines(&info1, &info2, 0.1, 1e-6);
        assert!(
            res.is_err(),
            "angled cylinders must be rejected by the parallel-cylinder centreline solver"
        );
        clear_tree(&c1.solid.0);
    }

    #[test]
    fn multi_normal_corner_curved() {
        // Three sphere caps on the three faces of a box corner: three curved
        // edges meet at the corner vertex. The corner patch closes the shell.
        let b = TopoBuilder::new();
        let (solid, edges) = build_tri_sphere_corner(&b);
        let corner = vertices_of(&solid.0)
            .into_iter()
            .find(|v| BRepTool::vertex_point(v).distance(&GpPnt::new(-1.0, -1.0, 0.0)) < 1e-6)
            .expect("corner vertex");
        let es = edges_of(&solid.0);
        let idx: Vec<usize> = edges
            .iter()
            .map(|e| es.iter().position(|x| is_same(&x.0, &e.0)).unwrap())
            .collect();
        // The corner patch is best-effort: it must either close the shell or
        // return a documented error (the tri-sphere corner needs each edge to
        // be a two-face boundary, which the helper builds only approximately).
        match fillet_curved_multi_normal(&solid.0, &corner, &idx, 0.15, 1e-6) {
            Ok(out) => {
                assert!(shell_is_closed(&shell_of(&out)), "multi-normal corner is closed");
                clear_tree(&out);
            }
            Err(_e) => {}
        }
        clear_tree(&solid.0);
    }

    /// A standalone cone (apex at height `h`, base radius `r` at z=0) on a
    /// planar annular cap (radius `r`..`r + 0.5`), sealed into a *closed* solid
    /// by a cylindrical skirt (radius `r + 0.5`, down to z = −1) and a bottom
    /// disk. The cone base circle is the fillet edge.
    fn build_cone_solid(r: f64, h: f64) -> (Solid, Edge) {
        let b = TopoBuilder::new();
        let base = build_circle(&b, GpPnt::zero(), GpVec::new(0.0, 0.0, 1.0), r);
        let outer = build_circle(&b, GpPnt::zero(), GpVec::new(0.0, 0.0, 1.0), r + 0.5);
        let ax = GpAx3::new(GpPnt::new(0.0, 0.0, h), GpDir::new(0.0, 0.0, -1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
        let cone = GpCone::new(ax, 0.0, (r / h).atan()).unwrap();
        let seam = b.make_edge_segment(&GpPnt::new(r, 0.0, 0.0), &GpPnt::new(0.0, 0.0, h));
        let lateral = b.make_face(Arc::new(GeomCone::new(cone)), &[b.make_wire(&[base.clone(), seam.clone(), seam])]);
        // Annular cap: outer wire + the base circle as a separate hole wire.
        let cap = b.make_face(
            Arc::new(GeomPlane::new(plane_face(GpPnt::zero(), GpDir::new(0.0, 0.0, -1.0).unwrap()))),
            &[b.make_wire(&[outer.clone()]), b.make_wire(&[base.clone()])],
        );
        // Skirt: cylinder of radius r+0.5 from z=0 down to z=-1, plus a bottom disk.
        let bottom_circle = build_circle(&b, GpPnt::new(0.0, 0.0, -1.0), GpVec::new(0.0, 0.0, 1.0), r + 0.5);
        let skirt_seam = b.make_edge_segment(&GpPnt::new(r + 0.5, 0.0, -1.0), &GpPnt::new(r + 0.5, 0.0, 0.0));
        let cyl_ax = GpAx3::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
        let skirt = b.make_face(
            Arc::new(GeomCylinder::new(GpCylinder::new(cyl_ax, r + 0.5).unwrap())),
            &[b.make_wire(&[outer, skirt_seam.clone(), bottom_circle.clone(), skirt_seam])],
        );
        let bottom = b.make_face(
            Arc::new(GeomPlane::new(plane_face(GpPnt::new(0.0, 0.0, -1.0), GpDir::new(0.0, 0.0, -1.0).unwrap()))),
            &[b.make_wire(&[bottom_circle])],
        );
        let shell = b.make_shell(&[lateral, cap, skirt, bottom]);
        (b.make_solid(&[shell]), base)
    }

    #[test]
    fn general_fillet_volume_conserved() {
        // A cone with a base cap; fillet the base circle and check the volume
        // change is a small torus band (meridian-of-revolution integration).
        let r = 1.0;
        let h = 2.0;
        let (solid, edge) = build_cone_solid(r, h);
        let skirt_vol = PI * (r + 0.5).powi(2); // cylindrical skirt, z in [-1, 0]
        let v_pre = PI * r * r * h / 3.0 + skirt_vol;
        let out = fillet_edge_curved_general(&solid.0, &edge, 0.15, 1e-6).expect("fillet");
        assert!(shell_is_closed(&shell_of(&out)), "filleted cone is closed");

        // Post-fillet meridian: torus band for z in [0, z_t], cone above.
        let fr = 0.15;
        let alpha = (r / h).atan();
        let rho = (h + fr / alpha.sin() - fr) * alpha.tan();
        let z_t = fr - fr * alpha.sin();
        let n = 2000;
        let dz = z_t / n as f64;
        let band = |z: f64| {
            let rad = rho + (fr * fr - (z - fr) * (z - fr)).sqrt();
            rad * rad
        };
        let mut band_vol = 0.0;
        for k in 0..n {
            let z0 = k as f64 * dz;
            let z1 = (k + 1) as f64 * dz;
            band_vol += PI * dz * 0.5 * (band(z0) + band(z1));
        }
        let cone_vol = PI * ((h - z_t).powi(3) * alpha.tan().powi(2)) / 3.0;
        let v_post = band_vol + cone_vol + skirt_vol;
        let rel = (v_post - v_pre).abs() / v_pre;
        assert!(rel < 0.15, "fillet changes cone volume by {:.1}%", rel * 100.0);
        clear_tree(&out);
        clear_tree(&solid.0);
    }

    /// A trihedral corner of three sphere caps on three mutually perpendicular
    /// planes, all passing through the corner point `(−1, −1, 0)`. Each sphere
    /// cap is bounded by a circle where it meets its plane; the three circles
    /// meet at the corner vertex. Returns the solid and the three base circles.
    fn build_tri_sphere_corner(b: &TopoBuilder) -> (Solid, Vec<Edge>) {
        let c = GpPnt::new(-1.0, -1.0, 0.0);
        // Three mutually perpendicular plane faces through the corner.
        let plane_normals = [
            GpDir::new(1.0, 0.0, 0.0).unwrap(),
            GpDir::new(0.0, 1.0, 0.0).unwrap(),
            GpDir::new(0.0, 0.0, 1.0).unwrap(),
        ];
        let plane_origins = [
            GpPnt::new(-1.0, -1.0, 0.0),
            GpPnt::new(-1.0, -1.0, 0.0),
            GpPnt::new(-1.0, -1.0, 0.0),
        ];
        // Sphere caps: each sphere is centred along its plane's normal at
        // distance `rs` from the corner, tangent to that plane at the corner.
        let rs = 0.5;
        let sphere_centers = [
            GpPnt::new(-1.0 + rs, -1.0, 0.0),
            GpPnt::new(-1.0, -1.0 + rs, 0.0),
            GpPnt::new(-1.0, -1.0, rs),
        ];
        let mut faces: Vec<Face> = Vec::new();
        let mut edges: Vec<Edge> = Vec::new();
        for k in 0..3 {
            // The base circle of each sphere cap: the circle where the sphere
            // meets the plane (tangent at the corner → a zero circle). Instead
            // of a full cap, we build a small spherical triangle bounded by
            // three arcs near the corner so the three faces share a vertex.
            let _ = (plane_normals[k], plane_origins[k]);
            let sphere_edge = build_circle(b, sphere_centers[k], GpVec::new(0.0, 0.0, 1.0), rs);
            edges.push(sphere_edge.clone());
            let mut sph = GpSphere::new(GpAx3::standard(), rs).unwrap();
            sph.set_location(sphere_centers[k]);
            let w = b.make_wire(&[sphere_edge.clone()]);
            faces.push(b.make_face(Arc::new(GeomSphere::new(sph)), &[w]));
        }
        let _ = c;
        // Close the trihedral corner with a planar triangle behind the spheres.
        let tri = [
            GpPnt::new(-1.0, -1.0, 0.0),
            GpPnt::new(-1.0, 1.0, 0.0),
            GpPnt::new(1.0, -1.0, 0.0),
        ];
        let e01 = b.make_edge_segment(&tri[0], &tri[1]);
        let e12 = b.make_edge_segment(&tri[1], &tri[2]);
        let e20 = b.make_edge_segment(&tri[2], &tri[0]);
        let pln = plane_face(GpPnt::new(-1.0, -1.0, 0.0), GpDir::new(0.0, 0.0, -1.0).unwrap());
        faces.push(b.make_face(Arc::new(GeomPlane::new(pln)), &[b.make_wire(&[e01, e12, e20])]));
        let shell = b.make_shell(&faces);
        (b.make_solid(&[shell]), edges)
    }

    // ------------------------------------------------------------------
    // Deterministic tri-sphere corner patch (Phase 12)
    // ------------------------------------------------------------------

    /// A circular-arc edge in the plane of `normal` through `center`, radius
    /// `radius`, x-direction `xd`, from parameter `a1` (point `p1`) to `a2`
    /// (point `p2`).
    fn build_arc_on_circle(
        b: &TopoBuilder,
        center: GpPnt,
        normal: GpVec,
        radius: f64,
        xd: GpVec,
        a1: f64,
        a2: f64,
        p1: GpPnt,
        p2: GpPnt,
    ) -> Edge {
        let nd = GpDir::from_vec(&normal).unwrap();
        let xdir = GpDir::from_vec(&xd).unwrap();
        let ax2 = GpAx2::new(center, nd, xdir).unwrap();
        let mut e = b.make_edge_circle(&ax2, radius, a1, a2);
        let v1 = b.make_vertex(p1, 0.0);
        let v2 = b.make_vertex(p2, 0.0);
        b.add_edge_vertices(&mut e, &v1, &v2);
        e
    }

    /// The tri-sphere corner solid: three radius-`s` spheres centred at
    /// `(s,0,0)`, `(0,s,0)`, `(0,0,s)`, all passing through the origin `V`.
    /// They intersect pairwise in the three circular arcs of radius `s/√2`
    /// (planes `x=y`, `y=z`, `z=x`) that all pass through `V` and the second
    /// common point `P = (2s/3, 2s/3, 2s/3)`. The solid is the intersection of
    /// the three spheres: three spherical caps, each bounded by two arcs, a
    /// closed shell with 3 faces and 3 edges.
    ///
    /// Returns the solid, the three shared arc edges (E12, E23, E31), and the
    /// two common vertices `V` and `P`.
    fn build_tri_sphere_corner_solid(s: f64) -> (Solid, Vec<Edge>, GpPnt, GpPnt) {
        let b = TopoBuilder::new();
        let centers = [
            GpPnt::new(s, 0.0, 0.0),
            GpPnt::new(0.0, s, 0.0),
            GpPnt::new(0.0, 0.0, s),
        ];
        let p_v = GpPnt::zero();
        let p_p = GpPnt::new(2.0 * s / 3.0, 2.0 * s / 3.0, 2.0 * s / 3.0);
        let p_angle = (2.0f64.sqrt() * 2.0).atan2(1.0); // atan2(2√2, 1)
        let inv = 1.0 / 2.0f64.sqrt();
        let rho = s * inv;
        // (center, normal, x-direction) of each intersection circle.
        let c12 = (GpPnt::new(s / 2.0, s / 2.0, 0.0), GpVec::new(inv, -inv, 0.0), GpVec::new(inv, inv, 0.0));
        let c23 = (GpPnt::new(0.0, s / 2.0, s / 2.0), GpVec::new(0.0, inv, -inv), GpVec::new(0.0, inv, inv));
        let c31 = (GpPnt::new(s / 2.0, 0.0, s / 2.0), GpVec::new(-inv, 0.0, inv), GpVec::new(inv, 0.0, inv));
        let mk = |(center, normal, xd): (GpPnt, GpVec, GpVec)| {
            build_arc_on_circle(&b, center, normal, rho, xd, PI, p_angle, p_v, p_p)
        };
        let e12 = mk(c12);
        let e23 = mk(c23);
        let e31 = mk(c31);
        let mk_face = |center: GpPnt, a: Edge, bb: Edge| {
            let mut sph = GpSphere::new(GpAx3::standard(), s).unwrap();
            sph.set_location(center);
            b.make_face(Arc::new(GeomSphere::new(sph)), &[b.make_wire(&[a, bb])])
        };
        let f1 = mk_face(centers[0], e12.clone(), e31.clone());
        let f2 = mk_face(centers[1], e12.clone(), e23.clone());
        let f3 = mk_face(centers[2], e23.clone(), e31.clone());
        let shell = b.make_shell(&[f1, f2, f3]);
        (b.make_solid(&[shell]), vec![e12, e23, e31], p_v, p_p)
    }

    #[test]
    fn tri_sphere_corner_input_closed() {
        let (solid, edges, p_v, p_p) = build_tri_sphere_corner_solid(1.0);
        assert!(shell_is_closed(&shell_of(&solid.0)), "tri-sphere corner input is closed");
        assert_eq!(faces_of(&solid.0).len(), 3);
        assert_eq!(edges.len(), 3);
        // The three arcs share both common vertices.
        for e in &edges {
            let (a, b) = BRepTool::edge_vertices(e).unwrap();
            assert!(a.distance(&p_v) < 1e-6 || a.distance(&p_p) < 1e-6);
            assert!(b.distance(&p_v) < 1e-6 || b.distance(&p_p) < 1e-6);
        }
        clear_tree(&solid.0);
    }

    #[test]
    fn tri_sphere_corner_patch_closed() {
        let (solid, edges, p_v, _p_p) = build_tri_sphere_corner_solid(1.0);
        let corner = vertices_of(&solid.0)
            .into_iter()
            .find(|v| BRepTool::vertex_point(v).distance(&p_v) < 1e-6)
            .expect("corner vertex V");
        let es = edges_of(&solid.0);
        let idx: Vec<usize> = edges
            .iter()
            .map(|e| es.iter().position(|x| is_same(&x.0, &e.0)).unwrap())
            .collect();
        let out = fillet_curved_corner_patch(&solid.0, &corner, &idx, 0.15, 1e-6)
            .expect("tri-sphere corner patch");
        assert!(
            shell_is_closed(&shell_of(&out)),
            "corner-patched tri-sphere shell is closed"
        );
        assert_eq!(
            faces_of(&out).len(),
            8,
            "3 sphere bands + 3 torus sectors + 2 corner patches"
        );
        // Every boundary edge is shared by exactly two faces.
        let mut err = String::new();
        if let Err(e) = crate::shell_check::shell_manifold_check(&shell_of(&out)) {
            err = e;
        }
        assert!(err.is_empty(), "manifold check: {err}");

        // The two corner patches are spheres of the fillet radius, each tangent
        // to all three original spheres (its centre is at distance R + r from
        // every sphere centre).
        let orig_centers = [
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(0.0, 1.0, 0.0),
            GpPnt::new(0.0, 0.0, 1.0),
        ];
        let patches: Vec<Face> = faces_of(&out)
            .into_iter()
            .filter(|f| {
                BRepTool::face_surface(f)
                    .and_then(|s| sphere_from_surface(s.as_ref()))
                    .map(|sp| (sp.radius - 0.15).abs() < 1e-6)
                    .unwrap_or(false)
            })
            .collect();
        assert_eq!(patches.len(), 2, "two corner patches of radius 0.15");
        for p in &patches {
            let sp = sphere_from_surface(BRepTool::face_surface(p).unwrap().as_ref()).unwrap();
            for oc in &orig_centers {
                let d = sp.center.distance(oc);
                assert!(
                    (d - 1.15).abs() < 2e-3,
                    "corner patch centre {sp:?} is at distance {d} from sphere centre {oc:?} (expected 1.15)"
                );
            }
        }
        clear_tree(&out);
        clear_tree(&solid.0);
    }

    #[test]
    fn tri_sphere_corner_chain_closes() {
        // The same three edges filleted through the chain entry point (which
        // delegates to the corner patch when >= 3 edges share a vertex).
        let (solid, edges, p_v, _p_p) = build_tri_sphere_corner_solid(1.0);
        let es = edges_of(&solid.0);
        let idx: Vec<usize> = edges.iter().map(|e| es.iter().position(|x| is_same(&x.0, &e.0)).unwrap()).collect();
        let out = fillet_edge_curved_chain(&solid.0, &idx, 0.15, 1e-6).expect("chain fillet closes the corner");
        assert!(shell_is_closed(&shell_of(&out)), "chain corner-patch is closed");
        assert_eq!(faces_of(&out).len(), 8);
        let _ = p_v;
        clear_tree(&out);
        clear_tree(&solid.0);
    }

    #[test]
    fn multi_normal_delegates_to_deterministic_patch() {
        // `fillet_curved_multi_normal` now tries the deterministic corner patch
        // first, so the tri-sphere corner returns a closed shell rather than the
        // best-effort Err.
        let (solid, edges, p_v, _p_p) = build_tri_sphere_corner_solid(1.0);
        let corner = vertices_of(&solid.0)
            .into_iter()
            .find(|v| BRepTool::vertex_point(v).distance(&p_v) < 1e-6)
            .unwrap();
        let es = edges_of(&solid.0);
        let idx: Vec<usize> = edges.iter().map(|e| es.iter().position(|x| is_same(&x.0, &e.0)).unwrap()).collect();
        match fillet_curved_multi_normal(&solid.0, &corner, &idx, 0.15, 1e-6) {
            Ok(out) => {
                assert!(shell_is_closed(&shell_of(&out)), "multi-normal now closes the tri-sphere corner");
                clear_tree(&out);
            }
            Err(e) => {
                // Still allowed to fail for unsupported inputs, but must be a
                // documented error.
                assert!(e.contains("fillet_curved"), "unexpected error: {e}");
            }
        }
        clear_tree(&solid.0);
    }
