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

mod prelude {

pub(crate) use std::f64::consts::PI;
pub(crate) use std::sync::Arc;
pub(crate) use occt_core::bnd::BndBox;
pub(crate) use occt_core::gp::{
    GpAx2, GpAx22d, GpAx3, GpCirc, GpCirc2d, GpCone, GpCylinder, GpDir, GpDir2d, GpLin, GpPln, GpPnt,
    GpPnt2d, GpSphere, GpTorus, GpVec,
};
pub(crate) use occt_geom::{Curve, GeomCircle, GeomCone, GeomCylinder, GeomLine, GeomPlane, GeomSphere, GeomTorus, Surface};
pub(crate) use occt_geom2d::{Geom2dCircle, Geom2dLine};
pub(crate) use crate::builder::TopoBuilder;
pub(crate) use crate::shape::{Edge, Face, Solid, Vertex};
pub(crate) use crate::tgeometry::GeometryRegistry;
pub(crate) use crate::topexp::Explorer;
pub(crate) use crate::abs::{Orientation, ShapeType};

}

mod statics;
mod brep_prim_box;
mod brep_prim_cylinder;
mod brep_prim_sphere;
mod brep_prim_cone;
mod brep_prim_torus;


#[cfg(test)]
mod tests {
    use super::prelude::*;
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
            let (u0, _u1, v0, _v1) = BRepTool::uv_bounds(&Face(f.clone()));
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
            let _s = BRepTool::face_surface(&Face((*f).clone())).unwrap();
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

pub use brep_prim_box::BRepPrimBox;
pub use brep_prim_cylinder::BRepPrimCylinder;
pub use brep_prim_sphere::BRepPrimSphere;
pub use brep_prim_cone::BRepPrimCone;
pub use brep_prim_torus::BRepPrimTorus;
