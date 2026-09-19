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

    pub(super) fn clear_tree(s: &TopoShape) {
        GeometryRegistry::global().clear_shape(s);
        let children = s.tshape.read().unwrap().children.clone();
        for c in children {
            clear_tree(&c);
        }
    }

    pub(super) fn shell_of(shape: &TopoShape) -> Shell {
        Shell(shape.tshape.read().unwrap().children[0].clone())
    }

    /// A planar face surface whose axis direction is `normal`.
    pub(super) fn plane_face(origin: GpPnt, normal: GpDir) -> GpPln {
        let x_dir = if normal.x().abs() > 0.9 {
            GpDir::new(0.0, 1.0, 0.0).unwrap()
        } else if normal.y().abs() > 0.9 {
            GpDir::new(0.0, 0.0, 1.0).unwrap()
        } else {
            GpDir::new(1.0, 0.0, 0.0).unwrap()
        };
        GpPln::new(GpAx3::new(origin, normal, &x_dir).unwrap())
    }

    /// A full circular edge with a fixed seam direction (angle 0). The seam
    /// direction is chosen perpendicular to the circle normal.
    pub(super) fn build_circle(b: &TopoBuilder, center: GpPnt, normal: GpVec, radius: f64) -> Edge {
        let nd = GpDir::from_xyz(&normal.xyz()).unwrap();
        let xd = if normal.xyz().x.abs() > 0.9 {
            GpDir::new(0.0, 1.0, 0.0).unwrap()
        } else {
            GpDir::new(1.0, 0.0, 0.0).unwrap()
        };
        let ax2 = GpAx2::new(center, nd, xd).unwrap();
        let mut e = b.make_edge_circle(&ax2, radius, 0.0, 2.0 * PI);
        let seam = center.translated_vec(&GpVec::from_xyz(xd.xyz()).multiplied_scalar(radius));
        let v = b.make_vertex(seam, 0.0);
        b.add_edge_vertices(&mut e, &v, &v);
        e
    }

    pub(super) fn edge_between(edges: &[Edge], i: usize, j: usize) -> Edge {
        let ep: [(usize, usize); 12] = [
            (0, 1), (1, 2), (2, 3), (3, 0),
            (4, 5), (5, 6), (6, 7), (7, 4),
            (0, 4), (1, 5), (2, 6), (3, 7),
        ];
        let idx = ep.iter().position(|&(a, b)| (a == i && b == j) || (a == j && b == i)).unwrap();
        edges[idx].clone()
    }

    pub(super) fn box_corners(lo: &GpPnt, hi: &GpPnt) -> [GpPnt; 8] {
        let (x0, y0, z0) = (lo.x(), lo.y(), lo.z());
        let (x1, y1, z1) = (hi.x(), hi.y(), hi.z());
        [
            GpPnt::new(x0, y0, z0),
            GpPnt::new(x1, y0, z0),
            GpPnt::new(x1, y1, z0),
            GpPnt::new(x0, y1, z0),
            GpPnt::new(x0, y0, z1),
            GpPnt::new(x1, y0, z1),
            GpPnt::new(x1, y1, z1),
            GpPnt::new(x0, y1, z1),
        ]
    }

    /// A box with `holes.len()` circular holes in its top face. The top face is
    /// a plane at `z = hi.z()` whose outer boundary is the box top square and
    /// whose inner boundary is one circle per hole. Returns the solid and the
    /// hole circle edges.
    pub(super) fn build_box_with_holes(
        b: &TopoBuilder,
        lo: &GpPnt,
        hi: &GpPnt,
        holes: &[(GpPnt, f64)],
    ) -> (Solid, Vec<Edge>) {
        let c = box_corners(lo, hi);
        let verts: Vec<Vertex> = c.iter().map(|p| b.make_vertex(*p, 0.0)).collect();
        let ep: [(usize, usize); 12] = [
            (0, 1), (1, 2), (2, 3), (3, 0),
            (4, 5), (5, 6), (6, 7), (7, 4),
            (0, 4), (1, 5), (2, 6), (3, 7),
        ];
        let mut edges: Vec<Edge> = Vec::new();
        for &(i, j) in &ep {
            let dir = GpDir::from_vec(&GpVec::from_pnts(&c[i], &c[j])).unwrap();
            let lin = GpLin::from_pnt_dir(c[i], dir);
            let mut e = b.make_edge(Arc::new(GeomLine::new(lin)), 0.0, c[i].distance(&c[j]));
            b.add_edge_vertices(&mut e, &verts[i], &verts[j]);
            edges.push(e);
        }

        // Bottom + four side faces.
        let normals = [
            GpDir::new(0.0, 0.0, -1.0).unwrap(),
            GpDir::new(0.0, -1.0, 0.0).unwrap(),
            GpDir::new(0.0, 1.0, 0.0).unwrap(),
            GpDir::new(-1.0, 0.0, 0.0).unwrap(),
            GpDir::new(1.0, 0.0, 0.0).unwrap(),
        ];
        let cycles: [[usize; 4]; 5] = [
            [0, 1, 2, 3], // bottom
            [0, 1, 5, 4], // -Y
            [3, 2, 6, 7], // +Y
            [0, 3, 7, 4], // -X
            [1, 2, 6, 5], // +X
        ];
        let origins = [c[0], c[0], c[3], c[0], c[1]];
        let mut faces: Vec<Face> = Vec::new();
        for i in 0..5 {
            let pln = plane_face(origins[i], normals[i]);
            let w = b.make_wire(&[
                edge_between(&edges, cycles[i][0], cycles[i][1]),
                edge_between(&edges, cycles[i][1], cycles[i][2]),
                edge_between(&edges, cycles[i][2], cycles[i][3]),
                edge_between(&edges, cycles[i][3], cycles[i][0]),
            ]);
            faces.push(b.make_face(Arc::new(GeomPlane::new(pln)), &[w]));
        }

        // Top annulus: outer square + one hole circle per hole.
        let top_w = b.make_wire(&[
            edge_between(&edges, 4, 5),
            edge_between(&edges, 5, 6),
            edge_between(&edges, 6, 7),
            edge_between(&edges, 7, 4),
        ]);
        let mut hole_edges = Vec::new();
        let mut hole_wires = Vec::new();
        for (hc, hr) in holes {
            let e = build_circle(b, *hc, GpVec::new(0.0, 0.0, 1.0), *hr);
            hole_edges.push(e.clone());
            hole_wires.push(b.make_wire(&[e]));
        }
        let mut top_wires = vec![top_w];
        top_wires.extend(hole_wires);
        let top_face = b.make_face(
            Arc::new(GeomPlane::new(plane_face(c[4], GpDir::new(0.0, 0.0, 1.0).unwrap()))),
            &top_wires,
        );
        faces.push(top_face);

        let shell = b.make_shell(&faces);
        let solid = b.make_solid(&[shell]);
        (solid, hole_edges)
    }

    /// A sphere cap sitting on the top face of a box: box `[-1.5,1.5]²×[-1,0]`
    /// with a sphere of radius `r_s` centred at the origin (on the top plane).
    pub(super) fn build_sphere_on_box(r_s: f64) -> (Solid, Edge) {
        let b = TopoBuilder::new();
        let (solid, holes) = build_box_with_holes(
            &b,
            &GpPnt::new(-1.5, -1.5, -1.0),
            &GpPnt::new(1.5, 1.5, 0.0),
            &[(GpPnt::zero(), r_s)],
        );
        let hole = holes[0].clone();
        let mut sph = GpSphere::new(GpAx3::standard(), r_s).unwrap();
        sph.set_location(GpPnt::zero());
        let cap = b.make_face(Arc::new(GeomSphere::new(sph)), &[b.make_wire(&[hole.clone()])]);
        let mut fs = faces_of(&solid.0);
        fs.push(cap);
        let shell = b.make_shell(&fs);
        (b.make_solid(&[shell]), hole)
    }

    /// A cylinder standing on the top face of a box: box `[-1.5,1.5]²×[-1,0]`
    /// with a cylinder of radius `r_c`, axis Z through the origin, height `h`.
    pub(super) fn build_cylinder_on_box(r_c: f64, h: f64) -> (Solid, Edge) {
        let b = TopoBuilder::new();
        let (solid, holes) = build_box_with_holes(
            &b,
            &GpPnt::new(-1.5, -1.5, -1.0),
            &GpPnt::new(1.5, 1.5, 0.0),
            &[(GpPnt::zero(), r_c)],
        );
        let base = holes[0].clone();
        let ax = GpAx3::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
        let top = build_circle(&b, GpPnt::new(0.0, 0.0, h), GpVec::new(0.0, 0.0, 1.0), r_c);
        let seam = b.make_edge_segment(&GpPnt::new(r_c, 0.0, 0.0), &GpPnt::new(r_c, 0.0, h));
        let lateral_wire = b.make_wire(&[base.clone(), seam.clone(), top.clone(), seam]);
        let lateral = b.make_face(Arc::new(GeomCylinder::new(GpCylinder::new(ax, r_c).unwrap())), &[lateral_wire]);
        let top_cap = b.make_face(
            Arc::new(GeomPlane::new(plane_face(GpPnt::new(0.0, 0.0, h), GpDir::new(0.0, 0.0, 1.0).unwrap()))),
            &[b.make_wire(&[top])],
        );
        let mut fs = faces_of(&solid.0);
        fs.push(lateral);
        fs.push(top_cap);
        let shell = b.make_shell(&fs);
        (b.make_solid(&[shell]), base)
    }

    /// Two radius-1 spheres at distance `d` along +X. The union solid is two
    /// caps sharing the intersection circle.
    pub(super) fn build_two_spheres(d: f64) -> (Solid, Edge) {
        let b = TopoBuilder::new();
        let r = 1.0;
        let a = (r * r - (d / 2.0) * (d / 2.0)).sqrt();
        let ic = build_circle(&b, GpPnt::new(d / 2.0, 0.0, 0.0), GpVec::new(1.0, 0.0, 0.0), a);
        let ic_wire = b.make_wire(&[ic.clone()]);
        let mut s1 = GpSphere::new(GpAx3::standard(), r).unwrap();
        s1.set_location(GpPnt::zero());
        let mut s2 = GpSphere::new(GpAx3::standard(), r).unwrap();
        s2.set_location(GpPnt::new(d, 0.0, 0.0));
        let cap1 = b.make_face(Arc::new(GeomSphere::new(s1)), &[ic_wire.clone()]);
        let cap2 = b.make_face(Arc::new(GeomSphere::new(s2)), &[ic_wire]);
        let shell = b.make_shell(&[cap1, cap2]);
        (b.make_solid(&[shell]), ic)
    }

    /// Two sphere caps on a box, for the chain test: box `[-2,2]²×[-1,0]`, two
    /// spheres of radius 0.7 centred at `(±0.9, 0, 0)`.
    pub(super) fn build_two_spheres_on_box() -> (Solid, Vec<Edge>) {
        let b = TopoBuilder::new();
        let (solid, holes) = build_box_with_holes(
            &b,
            &GpPnt::new(-2.0, -2.0, -1.0),
            &GpPnt::new(2.0, 2.0, 0.0),
            &[(GpPnt::new(-0.9, 0.0, 0.0), 0.7), (GpPnt::new(0.9, 0.0, 0.0), 0.7)],
        );
        let mut sph1 = GpSphere::new(GpAx3::standard(), 0.7).unwrap();
        sph1.set_location(GpPnt::new(-0.9, 0.0, 0.0));
        let mut sph2 = GpSphere::new(GpAx3::standard(), 0.7).unwrap();
        sph2.set_location(GpPnt::new(0.9, 0.0, 0.0));
        let cap1 = b.make_face(Arc::new(GeomSphere::new(sph1)), &[b.make_wire(&[holes[0].clone()])]);
        let cap2 = b.make_face(Arc::new(GeomSphere::new(sph2)), &[b.make_wire(&[holes[1].clone()])]);
        let mut fs = faces_of(&solid.0);
        fs.push(cap1);
        fs.push(cap2);
        let shell = b.make_shell(&fs);
        (b.make_solid(&[shell]), holes)
    }

    /// Axisymmetric solid for the volume test: a base cylinder of radius `r_b`,
    /// height 1 below z=0, topped by a sphere cap of radius `r_s` at the origin.
    pub(super) fn build_sphere_on_cylinder(r_b: f64, r_s: f64) -> (Solid, Edge) {
        let b = TopoBuilder::new();
        let ax = GpAx3::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
        let bottom_circle = build_circle(&b, GpPnt::new(0.0, 0.0, -1.0), GpVec::new(0.0, 0.0, 1.0), r_b);
        let outer_circle = build_circle(&b, GpPnt::zero(), GpVec::new(0.0, 0.0, 1.0), r_b);
        let inner_circle = build_circle(&b, GpPnt::zero(), GpVec::new(0.0, 0.0, 1.0), r_s);
        let seam = b.make_edge_segment(&GpPnt::new(r_b, 0.0, -1.0), &GpPnt::new(r_b, 0.0, 0.0));
        let lateral_wire = b.make_wire(&[bottom_circle.clone(), seam.clone(), outer_circle.clone(), seam]);
        let lateral = b.make_face(Arc::new(GeomCylinder::new(GpCylinder::new(ax, r_b).unwrap())), &[lateral_wire]);
        let bottom_cap = b.make_face(
            Arc::new(GeomPlane::new(plane_face(GpPnt::new(0.0, 0.0, -1.0), GpDir::new(0.0, 0.0, -1.0).unwrap()))),
            &[b.make_wire(&[bottom_circle])],
        );
        let top_annulus = b.make_face(
            Arc::new(GeomPlane::new(plane_face(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap()))),
            &[b.make_wire(&[outer_circle]), b.make_wire(&[inner_circle.clone()])],
        );
        let mut sph = GpSphere::new(GpAx3::standard(), r_s).unwrap();
        sph.set_location(GpPnt::zero());
        let cap = b.make_face(Arc::new(GeomSphere::new(sph)), &[b.make_wire(&[inner_circle.clone()])]);
        let shell = b.make_shell(&[bottom_cap, lateral, top_annulus, cap]);
        (b.make_solid(&[shell]), inner_circle)
    }

    /// The first face that is not plane/sphere/cylinder — the torus blend face,
    /// since every input face of the test solids is plane/sphere/cylinder.
    pub(super) fn find_blend_face<'a>(faces: &'a [Face]) -> &'a Face {
        faces
            .iter()
            .find(|f| {
                BRepTool::face_surface(f)
                    .map(|s| classify_surface_full(s.as_ref()) == SurfaceKind::Other)
                    .unwrap_or(false)
            })
            .expect("torus blend face")
    }

    /// Distance from `p` to the circle (centre `c`, unit normal `n`, radius
    /// `rho`) in the plane perpendicular to `n`.
    pub(super) fn distance_to_circle(p: &GpPnt, c: &GpPnt, n: &GpVec, rho: f64) -> f64 {
        let v = GpVec::from_pnts(c, p);
        let along = v.dot(n);
        let radial = v.subtracted(&n.multiplied_scalar(along));
        (radial.magnitude() - rho).hypot(along)
    }

    /// A right cone standing on the top face of a box: base radius `r`, apex at
    /// height `h` above the top plane (z=0), axis +Z. The base circle is the
    /// shared (fillet) edge between the cone lateral face and the box top.
    pub(super) fn build_cone_on_box(r: f64, h: f64) -> (Solid, Edge) {
        let b = TopoBuilder::new();
        let (solid, holes) = build_box_with_holes(
            &b,
            &GpPnt::new(-1.5, -1.5, -1.0),
            &GpPnt::new(1.5, 1.5, 0.0),
            &[(GpPnt::zero(), r)],
        );
        let base = holes[0].clone();
        // Cone: true apex at (0,0,h), axis −Z (radius grows toward the base).
        let ax = GpAx3::new(
            GpPnt::new(0.0, 0.0, h),
            GpDir::new(0.0, 0.0, -1.0).unwrap(),
            &GpDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap();
        let cone = GpCone::new(ax, 0.0, (r / h).atan()).unwrap();
        let seam = b.make_edge_segment(&GpPnt::new(r, 0.0, 0.0), &GpPnt::new(0.0, 0.0, h));
        let lateral_wire = b.make_wire(&[base.clone(), seam.clone(), seam]);
        let lateral = b.make_face(Arc::new(GeomCone::new(cone)), &[lateral_wire]);
        let mut fs = faces_of(&solid.0);
        fs.push(lateral);
        let shell = b.make_shell(&fs);
        (b.make_solid(&[shell]), base)
    }

    /// A flared pipe: a cylinder (radius `R`, z from −1 to 1) topped by a cone
    /// (apex at the origin, widening up, z from 1 to 2). The two share the
    /// circle at z=1, radius R (the fillet edge).
    pub(super) fn build_cylinder_cone_union(R: f64) -> (Solid, Edge) {
        let b = TopoBuilder::new();
        let ax = GpAx3::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
        let bottom = build_circle(&b, GpPnt::new(0.0, 0.0, -1.0), GpVec::new(0.0, 0.0, 1.0), R);
        let junction = build_circle(&b, GpPnt::new(0.0, 0.0, 1.0), GpVec::new(0.0, 0.0, 1.0), R);
        let seam_cyl = b.make_edge_segment(&GpPnt::new(R, 0.0, -1.0), &GpPnt::new(R, 0.0, 1.0));
        let cyl_lat = b.make_face(
            Arc::new(GeomCylinder::new(GpCylinder::new(ax, R).unwrap())),
            &[b.make_wire(&[bottom.clone(), seam_cyl.clone(), junction.clone(), seam_cyl])],
        );
        let bottom_cap = b.make_face(
            Arc::new(GeomPlane::new(plane_face(GpPnt::new(0.0, 0.0, -1.0), GpDir::new(0.0, 0.0, -1.0).unwrap()))),
            &[b.make_wire(&[bottom])],
        );
        // Cone: apex at origin, axis +Z, widening up, α = atan(R/1).
        let alpha = R.atan();
        let cone_ax = GpAx3::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
        let cone = GpCone::new(cone_ax, 0.0, alpha).unwrap();
        let top_radius = 2.0 * alpha.tan();
        let top_circle = build_circle(&b, GpPnt::new(0.0, 0.0, 2.0), GpVec::new(0.0, 0.0, 1.0), top_radius);
        let seam_cone = b.make_edge_segment(&GpPnt::new(R, 0.0, 1.0), &GpPnt::new(top_radius, 0.0, 2.0));
        let cone_lat = b.make_face(
            Arc::new(GeomCone::new(cone)),
            &[b.make_wire(&[junction.clone(), seam_cone.clone(), top_circle.clone(), seam_cone])],
        );
        let top_cap = b.make_face(
            Arc::new(GeomPlane::new(plane_face(GpPnt::new(0.0, 0.0, 2.0), GpDir::new(0.0, 0.0, 1.0).unwrap()))),
            &[b.make_wire(&[top_circle])],
        );
        let shell = b.make_shell(&[cyl_lat, bottom_cap, cone_lat, top_cap]);
        (b.make_solid(&[shell]), junction)
    }

    /// A cone (apex at the origin, widening up) capped by a sphere: the sphere
    /// centre is on the cone axis and its base circle sits on the cone's top
    /// rim (z=1, radius tan(α)). The rim circle is the fillet edge.
    pub(super) fn build_cone_sphere_union(alpha: f64, rs: f64) -> (Solid, Edge) {
        let b = TopoBuilder::new();
        let h = 1.0;
        let r_top = h * alpha.tan();
        let cone_ax = GpAx3::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
        let cone = GpCone::new(cone_ax, 0.0, alpha).unwrap();
        let rim = build_circle(&b, GpPnt::new(0.0, 0.0, h), GpVec::new(0.0, 0.0, 1.0), r_top);
        let seam = b.make_edge_segment(&GpPnt::new(r_top, 0.0, h), &GpPnt::zero());
        let cone_lat = b.make_face(Arc::new(GeomCone::new(cone)), &[b.make_wire(&[rim.clone(), seam.clone(), seam])]);
        let z_s = h + (rs * rs - r_top * r_top).sqrt();
        let mut sph = GpSphere::new(GpAx3::standard(), rs).unwrap();
        sph.set_location(GpPnt::new(0.0, 0.0, z_s));
        let cap = b.make_face(Arc::new(GeomSphere::new(sph)), &[b.make_wire(&[rim.clone()])]);
        let shell = b.make_shell(&[cone_lat, cap]);
        (b.make_solid(&[shell]), rim)
    }

    /// Two coaxial cones meeting at a shared circle (the fillet edge): a
    /// shallow cone (apex at the origin, α1) and a steeper cone (apex above,
    /// widening down, α2), meeting at z=1.
    pub(super) fn build_cone_cone_union(alpha1: f64, alpha2: f64) -> (Solid, Edge) {
        let b = TopoBuilder::new();
        let h = 1.0;
        let r_shared = h * alpha1.tan();
        let c1_ax = GpAx3::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
        let c1 = GpCone::new(c1_ax, 0.0, alpha1).unwrap();
        let shared = build_circle(&b, GpPnt::new(0.0, 0.0, h), GpVec::new(0.0, 0.0, 1.0), r_shared);
        let seam1 = b.make_edge_segment(&GpPnt::new(r_shared, 0.0, h), &GpPnt::zero());
        let lat1 = b.make_face(Arc::new(GeomCone::new(c1)), &[b.make_wire(&[shared.clone(), seam1.clone(), seam1])]);
        let h2 = h + r_shared / alpha2.tan();
        let c2_ax = GpAx3::new(GpPnt::new(0.0, 0.0, h2), GpDir::new(0.0, 0.0, -1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
        let c2 = GpCone::new(c2_ax, 0.0, alpha2).unwrap();
        let seam2 = b.make_edge_segment(&GpPnt::new(r_shared, 0.0, h), &GpPnt::new(0.0, 0.0, h2));
        let lat2 = b.make_face(Arc::new(GeomCone::new(c2)), &[b.make_wire(&[shared.clone(), seam2.clone(), seam2])]);
        let shell = b.make_shell(&[lat1, lat2]);
        (b.make_solid(&[shell]), shared)
    }

    /// A circular arc edge in the plane z=const with explicit angle range.
    pub(super) fn build_arc_xy(b: &TopoBuilder, center: GpPnt, radius: f64, a1: f64, a2: f64, z: f64) -> Edge {
        let c = GpPnt::new(center.x(), center.y(), z);
        let nd = GpDir::new(0.0, 0.0, 1.0).unwrap();
        let xd = GpDir::new(1.0, 0.0, 0.0).unwrap();
        let ax2 = GpAx2::new(c, nd, xd).unwrap();
        let mut e = b.make_edge_circle(&ax2, radius, a1, a2);
        let v1 = b.make_vertex(GpPnt::new(center.x() + radius * a1.cos(), center.y() + radius * a1.sin(), z), 0.0);
        let v2 = b.make_vertex(GpPnt::new(center.x() + radius * a2.cos(), center.y() + radius * a2.sin(), z), 0.0);
        b.add_edge_vertices(&mut e, &v1, &v2);
        e
    }

    /// Two equal overlapping parallel cylinders (radius 1, axes through (0,0,0)
    /// and (d,0,0), from z=0 to z=2). Returns the lens solid and the P1 seam
    /// edge (one of the two intersection lines, the fillet edge).
    pub(super) fn build_parallel_cylinders(d: f64) -> (Solid, Edge) {
        let b = TopoBuilder::new();
        let r = 1.0;
        let a_half = d / 2.0;
        let h_half = (r * r - a_half * a_half).sqrt();
        // P1 = (a_half, +h_half), P2 = (a_half, −h_half) in the cross-section.
        let p1 = GpPnt::new(a_half, h_half, 0.0);
        let p2 = GpPnt::new(a_half, -h_half, 0.0);
        let p1_top = GpPnt::new(a_half, h_half, 2.0);
        let p2_top = GpPnt::new(a_half, -h_half, 2.0);
        // Angles on circle A (centre 0,0) and circle B (centre d,0).
        let a1_p1 = h_half.atan2(a_half);
        let a1_p2 = (-h_half).atan2(a_half);
        let a2_p1 = h_half.atan2(a_half - d);
        let a2_p2 = (-h_half).atan2(a_half - d);
        // Outer arc of A: P1 → P2 counter-clockwise through 180°.
        let arc_a_b = build_arc_xy(&b, GpPnt::zero(), r, a1_p1, a1_p2 + 2.0 * PI, 0.0);
        let arc_a_t = build_arc_xy(&b, GpPnt::zero(), r, a1_p1, a1_p2 + 2.0 * PI, 2.0);
        // Outer arc of B: P1 → P2 clockwise through 0° (angle decreasing).
        let arc_b_b = build_arc_xy(&b, GpPnt::new(d, 0.0, 0.0), r, a2_p1, a2_p2, 0.0);
        let arc_b_t = build_arc_xy(&b, GpPnt::new(d, 0.0, 0.0), r, a2_p1, a2_p2, 2.0);
        // Shared seam lines at P1 and P2.
        let seam_p1 = b.make_edge_segment(&p1, &p1_top);
        let seam_p2 = b.make_edge_segment(&p2, &p2_top);
        // Lateral faces.
        let lat_a_surf = Arc::new(GeomCylinder::new(GpCylinder::new(GpAx3::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap(), r).unwrap()));
        let lat_a = b.make_face(lat_a_surf, &[b.make_wire(&[arc_a_b.clone(), seam_p2.clone(), arc_a_t.clone(), seam_p1.clone()])]);
        let lat_b_surf = Arc::new(GeomCylinder::new(GpCylinder::new(GpAx3::new(GpPnt::new(d, 0.0, 0.0), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap(), r).unwrap()));
        let lat_b = b.make_face(lat_b_surf, &[b.make_wire(&[arc_b_b.clone(), seam_p1.clone(), arc_b_t.clone(), seam_p2.clone()])]);
        // Caps: peanut-shaped planar faces.
        let bottom = b.make_face(Arc::new(GeomPlane::new(plane_face(GpPnt::zero(), GpDir::new(0.0, 0.0, -1.0).unwrap()))), &[b.make_wire(&[arc_a_b.clone(), arc_b_b.clone()])]);
        let top = b.make_face(Arc::new(GeomPlane::new(plane_face(GpPnt::new(0.0, 0.0, 2.0), GpDir::new(0.0, 0.0, 1.0).unwrap()))), &[b.make_wire(&[arc_a_t.clone(), arc_b_t.clone()])]);
        let shell = b.make_shell(&[lat_a, lat_b, bottom, top]);
        (b.make_solid(&[shell]), seam_p1)
    }

    // ------------------------------------------------------------------
    // Tests
    // ------------------------------------------------------------------

    #[test]
    fn plane_sphere_blend_closed() {
        let (solid, edge) = build_sphere_on_box(1.0);
        let out = fillet_edge_curved(&solid.0, &edge, 0.2, 1e-6).expect("plane+sphere fillet");
        let faces = faces_of(&out);
        assert_eq!(faces.len(), 8, "7 input faces + 1 blend face, got {}", faces.len());
        assert!(shell_is_closed(&shell_of(&out)), "plane+sphere fillet is closed");
        // The blend face is a torus band: every sampled point is at distance r
        // from the rolling-ball centreline circle.
        let blend = find_blend_face(&faces);
        let s = BRepTool::face_surface(blend).unwrap();
        let r = 0.2;
        let rho = (1.0_f64 + 2.0 * r).sqrt();
        let center = GpPnt::new(0.0, 0.0, r);
        let axis = GpVec::new(0.0, 0.0, 1.0);
        for (u, v) in [(0.3, 3.9), (1.0, 4.2), (2.0, 4.5), (4.0, 3.3), (5.0, 4.6)] {
            let p = s.d0(u, v);
            let d = distance_to_circle(&p, &center, &axis, rho);
            assert!((d - r).abs() < 1e-6, "blend point {p:?} at torus distance {d} (expected {r})");
        }
        clear_tree(&out);
        clear_tree(&solid.0);
    }

    #[test]
    fn plane_cylinder_blend() {
        let (solid, edge) = build_cylinder_on_box(1.0, 1.0);
        let out = fillet_edge_curved(&solid.0, &edge, 0.2, 1e-6).expect("plane+cylinder fillet");
        let faces = faces_of(&out);
        assert_eq!(faces.len(), 9, "8 input faces + 1 blend face, got {}", faces.len());
        assert!(shell_is_closed(&shell_of(&out)), "plane+cylinder fillet is closed");
        // The blend face is a torus band around the centreline circle of radius
        // R_c + r at height r above the plane.
        let blend = find_blend_face(&faces);
        let s = BRepTool::face_surface(blend).unwrap();
        let r = 0.2;
        let rho = 1.0 + r;
        let center = GpPnt::new(0.0, 0.0, r);
        let axis = GpVec::new(0.0, 0.0, 1.0);
        for (u, v) in [(0.3, 3.9), (1.0, 4.2), (2.0, 4.5), (4.0, 3.3), (5.0, 4.6)] {
            let p = s.d0(u, v);
            let d = distance_to_circle(&p, &center, &axis, rho);
            assert!((d - r).abs() < 1e-6, "blend point {p:?} at torus distance {d} (expected {r})");
        }
        clear_tree(&out);
        clear_tree(&solid.0);
    }

    #[test]
    fn blend_radius_tangent() {
        let (solid, edge) = build_cylinder_on_box(1.0, 1.0);
        let out = fillet_edge_curved(&solid.0, &edge, 0.2, 1e-6).unwrap();
        let faces = faces_of(&out);
        let blend = find_blend_face(&faces);
        let s = BRepTool::face_surface(blend).unwrap();
        let (u0, _u1) = s.u_range();
        // Plane contact (v = 3π/2): the blend sits on z=0 and its normal is
        // parallel to the plane normal.
        let p_plane = s.d0(u0, 3.0 * PI / 2.0);
        assert!(p_plane.z().abs() < 1e-9, "plane contact must lie on z=0, got {p_plane:?}");
        let n_plane = surface_normal(s.as_ref(), u0, 3.0 * PI / 2.0);
        let z_axis = GpVec::new(0.0, 0.0, 1.0);
        assert!(
            n_plane.xyz().crossed(&z_axis.xyz()).modulus() < 1e-3,
            "blend normal {:?} must be parallel to the plane normal",
            n_plane
        );
        // Cylinder contact (v = π): the blend sits on the radius-R_c cylinder
        // and its normal is radial (perpendicular to the axis).
        let p_cyl = s.d0(u0, PI);
        let radial_len = GpVec::new(p_cyl.x(), p_cyl.y(), 0.0).magnitude();
        assert!((radial_len - 1.0).abs() < 1e-6, "cylinder contact radius {radial_len}");
        let n_cyl = surface_normal(s.as_ref(), u0, PI);
        assert!(n_cyl.xyz().z.abs() < 1e-3, "cylinder-contact normal must be radial");
        let radial = GpVec::new(p_cyl.x(), p_cyl.y(), 0.0).normalized();
        assert!(
            n_cyl.xyz().crossed(&radial.xyz()).modulus() < 1e-3,
            "blend normal {:?} must be parallel to the cylinder radial",
            n_cyl
        );
        clear_tree(&out);
        clear_tree(&solid.0);
    }

    #[test]
    fn curved_unsupported_errors() {
        // Two torus faces sharing a circle edge: both classify as Other, so the
        // pair is unsupported and fillet_edge_curved must error (documented).
        let b = TopoBuilder::new();
        let circle = build_circle(&b, GpPnt::zero(), GpVec::new(0.0, 0.0, 1.0), 1.0);
        let w = b.make_wire(&[circle.clone()]);
        let f1 = b.make_face(Arc::new(GeomTorus::new(GpTorus::new(GpAx3::standard(), 3.0, 0.5).unwrap())), &[w.clone()]);
        let f2 = b.make_face(Arc::new(GeomTorus::new(GpTorus::new(GpAx3::standard(), 3.0, 0.5).unwrap())), &[w]);
        let shell = b.make_shell(&[f1.clone(), f2.clone()]);
        let solid = b.make_solid(&[shell]);
        let res = fillet_edge_curved(&solid.0, &circle, 0.2, 1e-6);
        assert!(res.is_err(), "torus+torus must be unsupported");
        // But faces_are_curved_compatible still reports the pair as curved.
        assert!(faces_are_curved_compatible(&f1, &f2));
        assert_eq!(classify_curved_pair(&f1, &f2), CurvedPair::Unsupported);
        clear_tree(&solid.0);
        clear_tree(&f1.0);
        clear_tree(&f2.0);
    }

    #[test]
    fn curved_chain_two_edges() {
        let (solid, holes) = build_two_spheres_on_box();
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
        assert_eq!(holes.len(), 2);
        let out = fillet_edge_curved_chain(&solid.0, &[i0, i1], 0.2, 1e-6).expect("two-edge chain");
        assert_eq!(faces_of(&out).len(), 10, "8 input faces + 2 blends");
        assert!(shell_is_closed(&shell_of(&out)), "two-edge curved chain is closed");
        clear_tree(&out);
        clear_tree(&solid.0);
    }

    #[test]
    fn sphere_sphere_blend() {
        let (solid, edge) = build_two_spheres(1.2);
        let out = fillet_edge_curved(&solid.0, &edge, 0.2, 1e-6).expect("sphere+sphere fillet");
        let faces = faces_of(&out);
        assert_eq!(faces.len(), 3, "2 input faces + 1 blend face, got {}", faces.len());
        assert!(shell_is_closed(&shell_of(&out)), "sphere+sphere fillet is closed");
        let blend = find_blend_face(&faces);
        let s = BRepTool::face_surface(blend).unwrap();
        let r = 0.2;
        let rho = (1.2_f64 * 1.2 - 0.6 * 0.6).sqrt();
        let center = GpPnt::new(0.6, 0.0, 0.0);
        let axis = GpVec::new(1.0, 0.0, 0.0);
        for (u, v) in [(0.3, 3.8), (1.0, 4.2), (2.0, 4.7), (4.0, 5.3), (5.0, 5.6)] {
            let p = s.d0(u, v);
            let d = distance_to_circle(&p, &center, &axis, rho);
            assert!((d - r).abs() < 1e-6, "blend point {p:?} at torus distance {d} (expected {r})");
        }
        clear_tree(&out);
        clear_tree(&solid.0);
    }

    #[test]
    fn rolling_center_line_sphere_plane() {
        let pln = GeomPlane::new(GpPln::new(GpAx3::standard()));
        let mut sph = GpSphere::new(GpAx3::standard(), 1.0).unwrap();
        sph.set_location(GpPnt::new(0.0, 0.0, 0.5));
        let sph_geom = GeomSphere::new(sph);
        let r = 0.2;
        let pts = rolling_ball_center_line(&pln, &sph_geom, r, 1e-6).expect("rolling centre line");
        assert_eq!(pts.len(), CENTERLINE_SAMPLES);
        for p in &pts {
            // The centre lies on the plane offset by r (here z = r).
            assert!((p.z() - r).abs() < 1e-9, "centre height {p:?}");
            // And on the sphere offset by R + r.
            let d = p.distance(&GpPnt::new(0.0, 0.0, 0.5));
            assert!((d - 1.2).abs() < 1e-9, "centre distance {d} (expected 1.2)");
        }
    }

    #[test]
    fn fillet_curved_volume_conserved() {
        let r_b = 1.5;
        let r_s = 1.0;
        let r = 0.2;
        let (solid, edge) = build_sphere_on_cylinder(r_b, r_s);
        let out = fillet_edge_curved(&solid.0, &edge, r, 1e-6).expect("fillet");
        assert!(shell_is_closed(&shell_of(&out)), "filleted axisymmetric solid is closed");

        // Pre-fillet volume: base cylinder (radius r_b, height 1) + upper
        // hemisphere of radius r_s.
        let v_pre = PI * r_b * r_b * 1.0 + (2.0 / 3.0) * PI * r_s.powi(3);

        // Post-fillet volume of revolution of the meridian profile:
        //   z < 0      : radius r_b  (base cylinder)
        //   0 ≤ z ≤ z_c: blend band  r(z) = rho − sqrt(r² − (z − r)²)
        //   z_c ≤ z ≤ r_s: sphere     r(z) = sqrt(r_s² − z²)
        let rho = (r_s * r_s + 2.0 * r_s * r).sqrt();
        let z_c = r * r_s / (r_s + r);
        let sphere_part = PI * ((r_s - z_c) - (r_s.powi(3) - z_c.powi(3)) / 3.0);
        let n = 4000;
        let h = z_c / n as f64;
        let band = |z: f64| {
            let rad = rho - (r * r - (z - r) * (z - r)).sqrt();
            rad * rad
        };
        let mut band_part = 0.0;
        for k in 0..n {
            let z0 = k as f64 * h;
            let z1 = (k + 1) as f64 * h;
            band_part += PI * h * 0.5 * (band(z0) + band(z1));
        }
        let v_post = PI * r_b * r_b * 1.0 + sphere_part + band_part;

        let rel = (v_post - v_pre).abs() / v_pre;
        assert!(rel < 0.15, "fillet changes volume by {:.2}% (pre {v_pre}, post {v_post})", rel * 100.0);
        clear_tree(&out);
        clear_tree(&solid.0);
    }

    #[test]
    fn faces_are_curved_compatible_mixed() {
        let (solid, edge) = build_sphere_on_box(1.0);
        let adjacent = faces_of(&solid.0)
            .into_iter()
            .filter(|f| {
                wires_of_face(f).iter().any(|w| {
                    edges_of_wire(w).iter().any(|e| is_same(&e.0, &edge.0))
                })
            })
            .collect::<Vec<Face>>();
        assert_eq!(adjacent.len(), 2);
        // One face is the plane, the other the sphere cap.
        assert!(faces_are_curved_compatible(&adjacent[0], &adjacent[1]));
        assert_eq!(classify_curved_pair(&adjacent[0], &adjacent[1]), CurvedPair::PlaneSphere);
        // Two planar box faces are NOT curved-compatible.
        let (boxy, _) = build_sphere_on_box(1.0);
        let all = faces_of(&boxy.0);
        let planars: Vec<&Face> = all.iter().filter(|f| is_planar(
            BRepTool::face_surface(f).unwrap().as_ref(), 8, 8, 1e-6,
        )).collect();
        assert!(planars.len() >= 2, "box has at least two planar faces");
        assert!(!faces_are_curved_compatible(planars[0], planars[1]));
        let _ = vertices_of(&solid.0);
        clear_tree(&solid.0);
        clear_tree(&boxy.0);
    }
