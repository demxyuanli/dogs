//! BRep primitive construction — boxes, cylinders, spheres, cones, toruses.
//! Source: `BRepPrimAPI_MakeBox`, `BRepPrimAPI_MakeCylinder`, `BRepPrimAPI_MakeSphere`,
//! `BRepPrimAPI_MakeCone`, `BRepPrimAPI_MakeTorus`.
//!
//! Simplified topological representation: a bottom-up shape skeleton
//! (vertices → edges → faces → solid) is built via [`TopoBuilder`], but the
//! flat child model stores counts only, not full boundary geometry.

use std::f64::consts::PI;

use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, Solid, TopoShape, Vertex};
use occt_core::bnd::BndBox;
use occt_core::gp::GpPnt;

/// Wireframe edges of a box: pairs of corner indices into `box_corners`.
const BOX_EDGES: [(usize, usize); 12] = [
    (0, 1), (1, 2), (2, 3), (3, 0), // bottom face
    (4, 5), (5, 6), (6, 7), (7, 4), // top face
    (0, 4), (1, 5), (2, 6), (3, 7), // vertical edges
];

/// Axis-aligned box primitive. Source: `BRepPrimAPI_MakeBox`.
#[derive(Debug, Clone)]
pub struct BRepPrimBox {
    pub solid: Solid,
    pub bbox: BndBox,
}

impl BRepPrimBox {
    /// Axis-aligned box from the origin with dimensions `dx × dy × dz`.
    pub fn make_box(dx: f64, dy: f64, dz: f64) -> Self {
        assert!(
            dx > 0.0 && dy > 0.0 && dz > 0.0,
            "BRepPrimBox::make_box: dimensions must be positive"
        );
        let b = TopoBuilder::new();

        // Bottom-up wireframe skeleton: 8 vertices, 12 edges.
        let corners = box_corners(dx, dy, dz);
        let verts: [Vertex; 8] = [
            b.make_vertex(corners[0], 0.0),
            b.make_vertex(corners[1], 0.0),
            b.make_vertex(corners[2], 0.0),
            b.make_vertex(corners[3], 0.0),
            b.make_vertex(corners[4], 0.0),
            b.make_vertex(corners[5], 0.0),
            b.make_vertex(corners[6], 0.0),
            b.make_vertex(corners[7], 0.0),
        ];
        for &(i, j) in &BOX_EDGES {
            let mut e = Edge::new();
            b.add(&mut e.0, &verts[i].0);
            b.add(&mut e.0, &verts[j].0);
        }

        // Boundary: 6 faces (4 edges each) forming 1 solid.
        let faces: Vec<Face> = (0..6).map(|_| make_face_with_edges(&b, 4)).collect();
        let solid = make_solid_with_faces(faces.len());

        let bbox = BndBox::from_corners(&GpPnt::zero(), &GpPnt::new(dx, dy, dz));
        Self { solid, bbox }
    }

    /// Box spanning two corner points (axis-aligned).
    pub fn make_box_corner(p1: &GpPnt, p2: &GpPnt) -> Self {
        let min = GpPnt::new(p1.x().min(p2.x()), p1.y().min(p2.y()), p1.z().min(p2.z()));
        let max = GpPnt::new(p1.x().max(p2.x()), p1.y().max(p2.y()), p1.z().max(p2.z()));
        let mut b = Self::make_box(max.x() - min.x(), max.y() - min.y(), max.z() - min.z());
        b.bbox = BndBox::from_corners(&min, &max);
        b
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
    pub fn make_cylinder(radius: f64, height: f64) -> Self {
        assert!(
            radius > 0.0 && height > 0.0,
            "BRepPrimCylinder::make_cylinder: radius and height must be positive"
        );
        // 3 faces: bottom cap, top cap, lateral surface.
        let solid = make_solid_with_faces(3);
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
    /// Sphere of given radius centered at the origin.
    pub fn make_sphere(radius: f64) -> Self {
        assert!(
            radius > 0.0,
            "BRepPrimSphere::make_sphere: radius must be positive"
        );
        // Single closed face.
        let solid = make_solid_with_faces(1);
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
    pub fn make_cone(radius: f64, height: f64) -> Self {
        assert!(
            radius > 0.0 && height > 0.0,
            "BRepPrimCone::make_cone: radius and height must be positive"
        );
        // 2 faces: base + lateral conical surface.
        let solid = make_solid_with_faces(2);
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
    pub fn make_torus(major: f64, minor: f64) -> Self {
        assert!(
            major > 0.0 && minor > 0.0,
            "BRepPrimTorus::make_torus: radii must be positive"
        );
        // Single closed face.
        let solid = make_solid_with_faces(1);
        Self { solid, major_radius: major, minor_radius: minor }
    }

    pub fn volume(&self) -> f64 {
        2.0 * PI * PI * self.major_radius * self.minor_radius.powi(2)
    }

    pub fn surface_area(&self) -> f64 {
        4.0 * PI * PI * self.major_radius * self.minor_radius
    }
}

/// Build a solid containing `n_faces` face children via [`TopoBuilder`].
fn make_solid_with_faces(n_faces: usize) -> Solid {
    let b = TopoBuilder::new();
    let mut solid = Solid::new();
    for _ in 0..n_faces {
        let face: TopoShape = Face::new().into();
        b.add(&mut solid.0, &face);
    }
    solid
}

/// Build a face containing `n_edges` edge children via [`TopoBuilder`].
fn make_face_with_edges(b: &TopoBuilder, n_edges: usize) -> Face {
    let mut face = Face::new();
    for _ in 0..n_edges {
        let edge: TopoShape = Edge::new().into();
        b.add(&mut face.0, &edge);
    }
    face
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

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9 * b.abs().max(1.0)
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
    fn box_topology_skeleton() {
        let b = BRepPrimBox::make_box(1.0, 2.0, 3.0);
        // Solid bounded by 6 faces.
        assert_eq!(b.solid.shape().tshape.read().unwrap().nb_children(), 6);
        // Bounding box spans origin..(dx,dy,dz).
        assert!(b.bbox.corner_min().is_equal(&GpPnt::zero()));
        assert!(b.bbox.corner_max().is_equal(&GpPnt::new(1.0, 2.0, 3.0)));
    }

    #[test]
    fn cylinder_volume() {
        let c = BRepPrimCylinder::make_cylinder(2.0, 10.0);
        assert!(approx(c.volume(), PI * 4.0 * 10.0));
        assert!(approx(c.lateral_area(), 2.0 * PI * 2.0 * 10.0));
    }

    #[test]
    fn sphere_volume() {
        let s = BRepPrimSphere::make_sphere(3.0);
        assert!(approx(s.volume(), 4.0 / 3.0 * PI * 27.0));
        assert!(approx(s.surface_area(), 4.0 * PI * 9.0));
    }

    #[test]
    fn cone_volume() {
        let c = BRepPrimCone::make_cone(3.0, 12.0);
        assert!(approx(c.volume(), PI * 9.0 * 12.0 / 3.0));
    }

    #[test]
    fn torus_volume() {
        let t = BRepPrimTorus::make_torus(5.0, 2.0);
        assert!(approx(t.volume(), 2.0 * PI * PI * 5.0 * 4.0));
        assert!(approx(t.surface_area(), 4.0 * PI * PI * 5.0 * 2.0));
    }

    #[test]
    #[should_panic]
    fn box_rejects_nonpositive_dimension() {
        BRepPrimBox::make_box(0.0, 1.0, 1.0);
    }
}
