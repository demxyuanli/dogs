//! BRep primitive construction — real boundary geometry.
//! Source: `BRepPrimAPI_MakeBox`, `BRepPrimAPI_MakeCylinder`, `BRepPrimAPI_MakeSphere`,
//! `BRepPrimAPI_MakeCone`, `BRepPrimAPI_MakeTorus`.
//!
//! Every primitive builds the full boundary representation with registered
//! geometry: vertices carry 3D points, edges carry curves (lines/circles),
//! wires group boundary edges, faces carry analytic surfaces (planes,
//! cylinders, spheres, cones, toruses), shells and solids assemble the
//! boundary. The geometry lives in the `GeometryRegistry` side-table, so
//! `BRepTool` can query it back exactly like OCCT's `BRep_Tool`.

use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::bnd::BndBox;
use occt_core::gp::{
    GpAx2, GpAx3, GpCirc, GpCone, GpCylinder, GpDir, GpLin, GpPln, GpPnt, GpSphere, GpTorus, GpVec,
};
use occt_geom::{Curve, GeomCircle, GeomCone, GeomCylinder, GeomLine, GeomPlane, GeomSphere, GeomTorus, Surface};

use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, Shell, Solid, Vertex, Wire};
use crate::topexp::Explorer;
use crate::abs::{Orientation, ShapeType};

/// Wireframe edges of a box: pairs of corner indices into `box_corners`.
const BOX_EDGES: [(usize, usize); 12] = [
    (0, 1), (1, 2), (2, 3), (3, 0), // bottom face
    (4, 5), (5, 6), (6, 7), (7, 4), // top face
    (0, 4), (1, 5), (2, 6), (3, 7), // vertical edges
];

/// Unit direction (panics on a degenerate axis — callers use valid axes).
fn dir(x: f64, y: f64, z: f64) -> GpDir {
    GpDir::new(x, y, z).expect("primitives: degenerate direction")
}

/// Line curve from `p1` to `p2`, parameterized on [0, |p2 − p1|].
fn line_curve(p1: &GpPnt, p2: &GpPnt) -> Arc<dyn Curve> {
    let d = GpDir::from_vec(&GpVec::from_pnts(p1, p2)).unwrap_or_else(|_| dir(1.0, 0.0, 0.0));
    Arc::new(GeomLine::new(GpLin::from_pnt_dir(*p1, d)))
}

/// Full circle curve in the plane `ax3` with the given radius (parameter 0..2π).
fn circle_curve(ax3: &GpAx3, radius: f64) -> Arc<dyn Curve> {
    let ax2 = GpAx2::new(ax3.location(), ax3.direction(), *ax3.x_direction())
        .expect("primitives: circle axis is not orthonormal");
    Arc::new(GeomCircle::new(GpCirc::new(ax2, radius)))
}

/// Plane through `origin` with the given outward `normal`; the in-plane X axis
/// is auto-selected perpendicular to the normal.
fn plane(origin: GpPnt, normal: GpDir) -> GpPln {
    let z = normal;
    let z_axis = dir(0.0, 0.0, 1.0);
    let x_dir = if z.is_normal(&z_axis) { z_axis } else { dir(1.0, 0.0, 0.0) };
    GpPln::new(GpAx3::new(origin, z, &x_dir).expect("primitives: non-orthogonal plane frame"))
}

/// Map a corner pair to the index of the box edge between them.
fn edge_index(a: usize, b: usize) -> usize {
    BOX_EDGES
        .iter()
        .position(|&(i, j)| (i == a && j == b) || (i == b && j == a))
        .unwrap_or_else(|| panic!("primitives: box edge ({a},{b}) not in BOX_EDGES"))
}

/// The 8 corners of an axis-aligned box from the origin, index order:
/// 0=(0,0,0), 1=(dx,0,0), 2=(dx,dy,0), 3=(0,dy,0),
/// 4=(0,0,dz), 5=(dx,0,dz), 6=(dx,dy,dz), 7=(0,dy,dz).
fn box_corners(dx: f64, dy: f64, dz: f64) -> [GpPnt; 8] {
    [
        GpPnt::new(0.0, 0.0, 0.0),
        GpPnt::new(dx, 0.0, 0.0),
        GpPnt::new(dx, dy, 0.0),
        GpPnt::new(0.0, dy, 0.0),
        GpPnt::new(0.0, 0.0, dz),
        GpPnt::new(dx, 0.0, dz),
        GpPnt::new(dx, dy, dz),
        GpPnt::new(0.0, dy, dz),
    ]
}

/// Build a real-geometry box solid. Faces are planar, edges are lines, the
/// 8 corners are shared by all incident edges (OCCT's shared-TShape model).
fn build_box(b: &TopoBuilder, corners: &[GpPnt; 8]) -> Solid {
    let verts: [Vertex; 8] = core::array::from_fn(|i| b.make_vertex(corners[i], 0.0));

    // 12 edges, each carrying a line curve and the two shared corner vertices.
    let mut edges: Vec<Edge> = Vec::with_capacity(12);
    for &(i, j) in &BOX_EDGES {
        let p1 = &corners[i];
        let p2 = &corners[j];
        let mut e = b.make_edge(line_curve(p1, p2), 0.0, p1.distance(p2));
        b.add(&mut e.0, &verts[i].0);
        b.add(&mut e.0, &verts[j].0);
        edges.push(e);
    }

    // 6 faces. Each face is a plane through 4 corners; its wire is the cyclic
    // list of the 4 boundary edges (each shared with a neighbouring face).
    let faces_def: [(GpPnt, GpDir, [usize; 4]); 6] = [
        // (a corner point on the face, outward normal, corner cycle CCW-outward)
        (corners[0], dir(0.0, 0.0, -1.0), [0, 3, 2, 1]), // -Z
        (corners[4], dir(0.0, 0.0, 1.0), [4, 5, 6, 7]), // +Z
        (corners[0], dir(0.0, -1.0, 0.0), [0, 1, 5, 4]), // -Y
        (corners[3], dir(0.0, 1.0, 0.0), [3, 7, 6, 2]), // +Y
        (corners[0], dir(-1.0, 0.0, 0.0), [0, 4, 7, 3]), // -X
        (corners[1], dir(1.0, 0.0, 0.0), [1, 2, 6, 5]), // +X
    ];

    let mut faces: Vec<Face> = Vec::with_capacity(6);
    for (origin, normal, cycle) in faces_def {
        // Orient each shared edge to the face's CCW-outward traversal: the edge
        // TShape carries one fixed curve direction (BOX_EDGES), so a wire that
        // traverses it the other way stores it Reversed (the OCCT
        // `BRepPrim_GWedge` / `BRepBuilderAPI_MakeWire` edge-orientation model).
        let quad: Vec<Edge> = [0, 1, 2, 3]
            .iter()
            .map(|&k| {
                let a = cycle[k];
                let b = cycle[(k + 1) % 4];
                let idx = edge_index(a, b);
                let e = edges[idx].clone();
                if BOX_EDGES[idx] == (a, b) {
                    e
                } else {
                    Edge(e.0.oriented(Orientation::Reversed))
                }
            })
            .collect();
        let wire = b.make_wire(&quad);
        let surface: Arc<dyn Surface> = Arc::new(GeomPlane::new(plane(origin, normal)));
        faces.push(b.make_face(surface, &[wire]));
    }

    let shell = b.make_shell(&faces);
    b.make_solid(&[shell])
}

/// Axis-aligned box primitive. Source: `BRepPrimAPI_MakeBox`.
#[derive(Debug, Clone)]
pub struct BRepPrimBox {
    pub solid: Solid,
    pub bbox: BndBox,
}

impl BRepPrimBox {
    /// Axis-aligned box from the origin with dimensions `dx × dy × dz`.
    /// Builds the full boundary: 8 vertices, 12 line edges, 6 planar faces.
    pub fn make_box(dx: f64, dy: f64, dz: f64) -> Self {
        assert!(
            dx > 0.0 && dy > 0.0 && dz > 0.0,
            "BRepPrimBox::make_box: dimensions must be positive"
        );
        let b = TopoBuilder::new();
        let corners = box_corners(dx, dy, dz);
        let solid = build_box(&b, &corners);
        let bbox = BndBox::from_corners(&GpPnt::zero(), &GpPnt::new(dx, dy, dz));
        Self { solid, bbox }
    }

    /// Box spanning two corner points (axis-aligned). The geometry is built
    /// at the actual corner position (a prism sweep of the base polygon), so
    /// `BRepTool`/`mesh_shape` report the correct world coordinates.
    pub fn make_box_corner(p1: &GpPnt, p2: &GpPnt) -> Self {
        let min = GpPnt::new(p1.x().min(p2.x()), p1.y().min(p2.y()), p1.z().min(p2.z()));
        let max = GpPnt::new(p1.x().max(p2.x()), p1.y().max(p2.y()), p1.z().max(p2.z()));
        let dz = max.z() - min.z();
        assert!(dz > 0.0, "make_box_corner: box must have positive height");
        let base = [
            GpPnt::new(min.x(), min.y(), min.z()),
            GpPnt::new(max.x(), min.y(), min.z()),
            GpPnt::new(max.x(), max.y(), min.z()),
            GpPnt::new(min.x(), max.y(), min.z()),
        ];
        let prism = crate::sweep::prism_from_polygon(&base, &occt_core::gp::GpVec::new(0.0, 0.0, 1.0), dz);
        Self { solid: prism.solid, bbox: BndBox::from_corners(&min, &max) }
    }

    pub fn volume(&self) -> f64 {
        let (x0, x1, y0, y1, z0, z1) = self.bbox.get().unwrap();
        (x1 - x0) * (y1 - y0) * (z1 - z0)
    }

    pub fn surface_area(&self) -> f64 {
        let (x0, x1, y0, y1, z0, z1) = self.bbox.get().unwrap();
        let dx = x1 - x0;
        let dy = y1 - y0;
        let dz = z1 - z0;
        2.0 * (dx * dy + dy * dz + dx * dz)
    }

    pub fn center(&self) -> GpPnt {
        let (x0, x1, y0, y1, z0, z1) = self.bbox.get().unwrap();
        GpPnt::new((x0 + x1) / 2.0, (y0 + y1) / 2.0, (z0 + z1) / 2.0)
    }

    /// Number of distinct vertices/edges/faces in the real boundary tree
    /// (uniquified by `TShape` identity, matching OCCT's `TopExp_Explorer` +
    /// `TopTools_MapOfShape`).
    pub fn counts(&self) -> (usize, usize, usize) {
        fn unique(shape: &crate::shape::TopoShape, t: ShapeType) -> usize {
            let mut seen = std::collections::HashSet::new();
            let mut ex = Explorer::new(shape, t);
            while ex.more() {
                seen.insert(Arc::as_ptr(&ex.current().tshape));
                ex.next();
            }
            seen.len()
        }
        (unique(&self.solid.0, ShapeType::Vertex), unique(&self.solid.0, ShapeType::Edge), unique(&self.solid.0, ShapeType::Face))
    }
}

/// Cylinder primitive. Source: `BRepPrimAPI_MakeCylinder`.
#[derive(Debug, Clone)]
pub struct BRepPrimCylinder {
    pub solid: Solid,
    pub radius: f64,
    pub height: f64,
}

impl BRepPrimCylinder {
    /// Cylinder of given radius and height, axis along Z.
    /// Boundary: 3 faces (bottom/top caps + lateral), 3 edges
    /// (two cap circles + a seam line), 2 vertices (seam endpoints).
    pub fn make_cylinder(radius: f64, height: f64) -> Self {
        assert!(
            radius > 0.0 && height > 0.0,
            "BRepPrimCylinder::make_cylinder: radius and height must be positive"
        );
        let b = TopoBuilder::new();
        let ax = GpAx3::new(GpPnt::zero(), dir(0.0, 0.0, 1.0), &dir(1.0, 0.0, 0.0))
            .expect("cylinder axis");

        let bottom = GpPnt::new(radius, 0.0, 0.0);
        let top = GpPnt::new(radius, 0.0, height);
        let v_bottom = b.make_vertex(bottom, 0.0);
        let v_top = b.make_vertex(top, 0.0);

        // Ring circles: parameter [0, 2π], seam at θ = 0 → shared vertices.
        let mut bottom_circle = b.make_edge(circle_curve(&ax, radius), 0.0, 2.0 * PI);
        b.add(&mut bottom_circle.0, &v_bottom.0);
        b.add(&mut bottom_circle.0, &v_bottom.0); // closed: both ends at seam vertex

        let mut top_ax = ax;
        top_ax.set_location(GpPnt::new(0.0, 0.0, height));
        let mut top_circle = b.make_edge(circle_curve(&top_ax, radius), 0.0, 2.0 * PI);
        b.add(&mut top_circle.0, &v_top.0);
        b.add(&mut top_circle.0, &v_top.0);

        // Seam line connecting the two seam vertices.
        let mut seam = b.make_edge(line_curve(&bottom, &top), 0.0, height);
        b.add(&mut seam.0, &v_bottom.0);
        b.add(&mut seam.0, &v_top.0);

        // Bottom cap: planar face bounded by the bottom circle.
        let bottom_wire = b.make_wire(&[bottom_circle.clone()]);
        let bottom_face = b.make_face(
            Arc::new(GeomPlane::new(plane(GpPnt::zero(), dir(0.0, 0.0, -1.0)))),
            &[bottom_wire],
        );

        // Top cap: planar face bounded by the top circle.
        let top_wire = b.make_wire(&[top_circle.clone()]);
        let top_face = b.make_face(
            Arc::new(GeomPlane::new(plane(GpPnt::new(0.0, 0.0, height), dir(0.0, 0.0, 1.0)))),
            &[top_wire],
        );

        // Lateral face: the cylinder surface, wire = bottom circle + seam +
        // top circle + seam (the seam is traversed twice, once per direction,
        // matching OCCT's seam edge convention).
        let lateral_wire = b.make_wire(&[bottom_circle, seam.clone(), top_circle, seam]);
        let lateral_face = b.make_face(Arc::new(GeomCylinder::new(
            GpCylinder::new(ax, radius).expect("cylinder radius"),
        )), &[lateral_wire]);

        let shell = b.make_shell(&[bottom_face, top_face, lateral_face]);
        let solid = b.make_solid(&[shell]);
        Self { solid, radius, height }
    }

    pub fn volume(&self) -> f64 {
        PI * self.radius * self.radius * self.height
    }

    pub fn lateral_area(&self) -> f64 {
        2.0 * PI * self.radius * self.height
    }
}

/// Sphere primitive. Source: `BRepPrimAPI_MakeSphere`.
#[derive(Debug, Clone)]
pub struct BRepPrimSphere {
    pub solid: Solid,
    pub radius: f64,
}

impl BRepPrimSphere {
    /// Sphere of given radius centered at the origin. A single face with the
    /// analytic spherical surface (no boundary edges in the untrimmed case).
    pub fn make_sphere(radius: f64) -> Self {
        assert!(
            radius > 0.0,
            "BRepPrimSphere::make_sphere: radius must be positive"
        );
        let b = TopoBuilder::new();
        let face = b.make_face(
            Arc::new(GeomSphere::new(
                GpSphere::new(GpAx3::standard(), radius).expect("sphere radius"),
            )),
            &[],
        );
        let shell = b.make_shell(&[face]);
        let solid = b.make_solid(&[shell]);
        Self { solid, radius }
    }

    pub fn volume(&self) -> f64 {
        4.0 / 3.0 * PI * self.radius.powi(3)
    }

    pub fn surface_area(&self) -> f64 {
        4.0 * PI * self.radius * self.radius
    }
}

/// Cone primitive. Source: `BRepPrimAPI_MakeCone`.
#[derive(Debug, Clone)]
pub struct BRepPrimCone {
    pub solid: Solid,
    pub radius: f64,
    pub height: f64,
}

impl BRepPrimCone {
    /// Right cone of given base radius and height, axis along Z.
    /// Boundary: 2 faces (base + lateral conical surface), 2 edges
    /// (base circle + seam line), 2 vertices (base seam + apex).
    pub fn make_cone(radius: f64, height: f64) -> Self {
        assert!(
            radius > 0.0 && height > 0.0,
            "BRepPrimCone::make_cone: radius and height must be positive"
        );
        let b = TopoBuilder::new();
        let ax = GpAx3::new(GpPnt::zero(), dir(0.0, 0.0, 1.0), &dir(1.0, 0.0, 0.0))
            .expect("cone axis");

        let apex = GpPnt::new(0.0, 0.0, height);
        let base_p = GpPnt::new(radius, 0.0, 0.0);
        let v_apex = b.make_vertex(apex, 0.0);
        let v_base = b.make_vertex(base_p, 0.0);

        // Base circle at z = 0.
        let mut base_circle = b.make_edge(circle_curve(&ax, radius), 0.0, 2.0 * PI);
        b.add(&mut base_circle.0, &v_base.0);
        b.add(&mut base_circle.0, &v_base.0);

        // Seam from the base seam point up to the apex.
        let mut seam = b.make_edge(line_curve(&base_p, &apex), 0.0, GpPnt::new(radius, 0.0, 0.0).distance(&apex));
        b.add(&mut seam.0, &v_base.0);
        b.add(&mut seam.0, &v_apex.0);

        // Base: planar face bounded by the base circle.
        let base_wire = b.make_wire(&[base_circle.clone()]);
        let base_face = b.make_face(
            Arc::new(GeomPlane::new(plane(GpPnt::zero(), dir(0.0, 0.0, -1.0)))),
            &[base_wire],
        );

        // Lateral: conical surface, wire = base circle + seam + seam.
        let semi_angle = (radius / height).atan();
        let lat_wire = b.make_wire(&[base_circle, seam.clone(), seam]);
        let lateral_face = b.make_face(Arc::new(GeomCone::new(
            GpCone::new(ax, radius, semi_angle).expect("cone angle"),
        )), &[lat_wire]);

        let shell = b.make_shell(&[base_face, lateral_face]);
        let solid = b.make_solid(&[shell]);
        Self { solid, radius, height }
    }

    pub fn volume(&self) -> f64 {
        PI * self.radius * self.radius * self.height / 3.0
    }
}

/// Torus primitive. Source: `BRepPrimAPI_MakeTorus`.
#[derive(Debug, Clone)]
pub struct BRepPrimTorus {
    pub solid: Solid,
    pub major_radius: f64,
    pub minor_radius: f64,
}

impl BRepPrimTorus {
    /// Torus with given major (centerline) and minor (tube) radii.
    /// A single face with the analytic toroidal surface.
    pub fn make_torus(major: f64, minor: f64) -> Self {
        assert!(
            major > 0.0 && minor > 0.0,
            "BRepPrimTorus::make_torus: radii must be positive"
        );
        let b = TopoBuilder::new();
        let face = b.make_face(
            Arc::new(GeomTorus::new(
                GpTorus::new(GpAx3::standard(), major, minor).expect("torus radii"),
            )),
            &[],
        );
        let shell = b.make_shell(&[face]);
        let solid = b.make_solid(&[shell]);
        Self { solid, major_radius: major, minor_radius: minor }
    }

    pub fn volume(&self) -> f64 {
        2.0 * PI * PI * self.major_radius * self.minor_radius.powi(2)
    }

    pub fn surface_area(&self) -> f64 {
        4.0 * PI * PI * self.major_radius * self.minor_radius
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;
    use crate::brep_tool::BRepTool;
    use crate::topexp::{Explorer, nb_shapes};

    fn approx(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9 * b.abs().max(1.0)
    }

    /// Collect distinct sub-shapes of type `t` (uniquified by TShape identity,
    /// since the boundary tree shares vertices and edges between faces).
    fn collect(solid: &Solid, t: ShapeType) -> Vec<crate::shape::TopoShape> {
        let mut seen = HashSet::new();
        let mut out = Vec::new();
        let mut ex = Explorer::new(&solid.0, t);
        while ex.more() {
            let s = ex.current().clone();
            if seen.insert(Arc::as_ptr(&s.tshape)) {
                out.push(s);
            }
            ex.next();
        }
        out
    }

    #[test]
    fn box_volume_and_center() {
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        assert!(approx(b.volume(), 2.0 * 3.0 * 4.0));
        assert!(approx(b.surface_area(), 2.0 * (2.0 * 3.0 + 3.0 * 4.0 + 2.0 * 4.0)));
        let c = b.center();
        assert!(approx(c.x(), 1.0));
        assert!(approx(c.y(), 1.5));
        assert!(approx(c.z(), 2.0));
    }

    #[test]
    fn box_from_corners() {
        let b = BRepPrimBox::make_box_corner(&GpPnt::new(1.0, 1.0, 1.0), &GpPnt::new(3.0, 4.0, 5.0));
        assert!(approx(b.volume(), 2.0 * 3.0 * 4.0));
        assert!(b.bbox.corner_min().is_equal(&GpPnt::new(1.0, 1.0, 1.0)));
        assert!(b.bbox.corner_max().is_equal(&GpPnt::new(3.0, 4.0, 5.0)));
    }

    #[test]
    fn box_bbox_contains_center() {
        let b = BRepPrimBox::make_box(1.0, 2.0, 3.0);
        assert!(!b.bbox.is_out(&b.center()));
    }

    #[test]
    fn box_real_geometry_tree() {
        let b = BRepPrimBox::make_box(1.0, 2.0, 3.0);
        // Solid → 1 shell → 6 faces.
        let solid = &b.solid;
        assert_eq!(solid.shape().tshape.read().unwrap().nb_children(), 1);
        let shells = collect(solid, ShapeType::Shell);
        assert_eq!(shells.len(), 1);
        assert_eq!(shells[0].tshape.read().unwrap().nb_children(), 6);

        // 8 vertices with registered points, 12 edges with curves.
        let verts = collect(solid, ShapeType::Vertex);
        let edges = collect(solid, ShapeType::Edge);
        let faces = collect(solid, ShapeType::Face);
        assert_eq!(verts.len(), 8);
        assert_eq!(edges.len(), 12);
        assert_eq!(faces.len(), 6);
        assert_eq!(nb_shapes(&collect(solid, ShapeType::Wire), ShapeType::Wire), 6);
    }

    #[test]
    fn box_geometry_queryable_via_brep_tool() {
        let b = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let verts = collect(&b.solid, ShapeType::Vertex);
        // All 8 corner coordinates are present.
        let points: Vec<GpPnt> = verts.iter().map(|v| BRepTool::vertex_point(&Vertex(v.clone()))).collect();
        for &(x, y, z) in &[(0.,0.,0.), (2.,0.,0.), (2.,2.,0.), (0.,2.,0.),
                            (0.,0.,2.), (2.,0.,2.), (2.,2.,2.), (0.,2.,2.)] {
            assert!(points.iter().any(|p| p.is_equal(&GpPnt::new(x, y, z))),
                    "missing corner ({x},{y},{z})");
        }
        // Every face carries a planar surface.
        let faces = collect(&b.solid, ShapeType::Face);
        for f in faces {
            let surf = BRepTool::face_surface(&Face(f.clone())).expect("face has surface");
            let (u0, u1, v0, v1) = BRepTool::uv_bounds(&Face(f.clone()));
            // Planes are unbounded: ranges are ±infinity.
            assert!(u0.is_infinite() && v0.is_infinite());
            let _ = surf;
        }
        // A point on the top face evaluates on its plane.
        let b2 = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let faces2 = collect(&b2.solid, ShapeType::Face);
        let top = faces2.iter().find(|f| BRepTool::face_surface(&Face((*f).clone()))
            .map(|s| (s.d0(0.0, 0.0).z() - 1.0).abs() < 1e-12).unwrap_or(false)).expect("top face");
        let p = BRepTool::face_surface(&Face(top.clone())).unwrap().d0(0.0, 0.0);
        assert!(p.is_equal(&GpPnt::new(0.0, 0.0, 1.0)));
    }

    #[test]
    fn box_counts_are_classic() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let (nv, ne, nf) = b.counts();
        assert_eq!((nv, ne, nf), (8, 12, 6));
    }

    #[test]
    fn cylinder_volume_and_geometry() {
        let c = BRepPrimCylinder::make_cylinder(2.0, 10.0);
        assert!(approx(c.volume(), PI * 4.0 * 10.0));
        assert!(approx(c.lateral_area(), 2.0 * PI * 2.0 * 10.0));
        let faces = collect(&c.solid, ShapeType::Face);
        assert_eq!(faces.len(), 3);
        // At least one face is a cylinder (non-planar) and two are planar.
        let planar = faces.iter().filter(|f| {
            let s = BRepTool::face_surface(&Face((*f).clone())).unwrap();
            let (u0, u1, _, _) = BRepTool::uv_bounds(&Face((*f).clone()));
            u1 - u0 > 1e9 // unbounded range → planar surface
        }).count();
        assert_eq!(planar, 2);
    }

    #[test]
    fn sphere_volume_and_geometry() {
        let s = BRepPrimSphere::make_sphere(3.0);
        assert!(approx(s.volume(), 4.0 / 3.0 * PI * 27.0));
        assert!(approx(s.surface_area(), 4.0 * PI * 9.0));
        let faces = collect(&s.solid, ShapeType::Face);
        assert_eq!(faces.len(), 1);
        let surf = BRepTool::face_surface(&Face(faces[0].clone())).unwrap();
        let (u0, u1, v0, v1) = BRepTool::uv_bounds(&Face(faces[0].clone()));
        assert!((u1 - u0 - 2.0 * PI).abs() < 1e-9);
        assert!((v1 - v0 - PI).abs() < 1e-9);
        // Sphere parameterization: x=r·cos v·cos u, y=r·cos v·sin u, z=r·sin v.
        // At (u,v)=(π/2, π/2) → (0, 0, r).
        let p = surf.d0(0.5 * PI, 0.5 * PI);
        assert!(p.is_equal(&GpPnt::new(0.0, 0.0, 3.0)) || p.distance(&GpPnt::new(0.0, 0.0, 3.0)) < 1e-9);
    }

    #[test]
    fn cone_volume_and_geometry() {
        let c = BRepPrimCone::make_cone(3.0, 12.0);
        assert!(approx(c.volume(), PI * 9.0 * 12.0 / 3.0));
        let faces = collect(&c.solid, ShapeType::Face);
        assert_eq!(faces.len(), 2);
    }

    #[test]
    fn torus_volume_and_geometry() {
        let t = BRepPrimTorus::make_torus(5.0, 2.0);
        assert!(approx(t.volume(), 2.0 * PI * PI * 5.0 * 4.0));
        assert!(approx(t.surface_area(), 4.0 * PI * PI * 5.0 * 2.0));
        let faces = collect(&t.solid, ShapeType::Face);
        assert_eq!(faces.len(), 1);
    }

    #[test]
    #[should_panic]
    fn box_rejects_nonpositive_dimension() {
        BRepPrimBox::make_box(0.0, 1.0, 1.0);
    }
}
