//! Unit tests for `occt_topology::brep_surface`, kept in `occt-topo` because
//! the reference solid needs `crate::primitives::BRepPrimBox` and the p-curve
//! projection lives in `crate::pcurve`.

#[cfg(test)]
mod brep_surface_gates {
    use std::sync::Arc;

    use occt_core::gp::{GpAx3, GpPln, GpPnt};
    use occt_geom::{GeomPlane, GeomSphere};
    use occt_topology::brep_surface::*;
    use occt_topology::brep_tool::BRepTool;
    use occt_topology::builder::TopoBuilder;
    use occt_topology::shape::Face;
    use occt_topology::tgeometry::GeometryRegistry;

    use crate::primitives::BRepPrimBox;
    use crate::pcurve::edge_pcurve_on_face;

    #[test]
    fn box_faces_are_planar() {
        let b = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        for f in occt_topology::topo_tools::faces_of(&b.solid.0) {
            assert!(face_is_plane_type(&f));
            let pln = face_plane(&f).expect("planar face has plane");
            let _ = pln;
            assert_eq!(
                classify_surface(GeometryRegistry::global().face_surface(&f.0).unwrap().as_ref()),
                SurfaceKind::Plane
            );
        }
    }

    #[test]
    fn sphere_face_classified() {
        let b = TopoBuilder::new();
        let face = b.make_face(
            Arc::new(GeomSphere::new(
                occt_core::gp::GpSphere::new(GpAx3::standard(), 2.0).unwrap(),
            )),
            &[],
        );
        let s = GeometryRegistry::global().face_surface(&face.0).unwrap();
        assert_eq!(classify_surface(s.as_ref()), SurfaceKind::Sphere);
        assert!(!face_is_plane_type(&face));
    }

    #[test]
    fn plane_reconstructed() {
        let b = TopoBuilder::new();
        let pln = GpPln::new(GpAx3::standard());
        let face = b.make_face(Arc::new(GeomPlane::new(pln.clone())), &[]);
        let got = face_plane(&face).expect("plane face");
        // Normal is Z in both.
        assert!(
            got.axis()
                .direction()
                .xyz()
                .crossed(pln.axis().direction().xyz())
                .modulus()
                < 1e-9
        );
    }

    #[test]
    fn box_face_normal_and_centroid() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 2.0);
        let top = occt_topology::topo_tools::faces_of(&b.solid.0)
            .into_iter()
            .find(|f| {
                BRepTool::face_surface(&Face(f.0.clone()))
                    .map(|s| (s.d0(0.0, 0.0).z() - 2.0).abs() < 1e-9)
                    .unwrap_or(false)
            })
            .expect("top face");
        let n = face_normal(&Face(top.0.clone()), 0.0, 0.0).unwrap();
        assert!(n.xyz().z > 0.9, "top face normal should point +Z");
        let c = face_centroid(&Face(top.0.clone()), 8, 8).unwrap();
        assert!((c.z() - 2.0).abs() < 1e-6);
    }

    #[test]
    fn pcurve_projection_is_recoverable() {
        let b = TopoBuilder::new();
        let pln = GpPln::new(GpAx3::standard());
        let face = b.make_face(Arc::new(GeomPlane::new(pln)), &[]);
        let e = b.make_edge_segment(&GpPnt::new(0., 0., 0.), &GpPnt::new(1., 0., 0.));
        let pc = edge_pcurve_on_face(&e, &Face(face.0.clone()), 5);
        assert_eq!(pc.len(), 5);
        // On the Z=0 plane, u≈x, v≈y.
        assert!((pc[0].x() - 0.0).abs() < 1e-6);
        assert!((pc[4].x() - 1.0).abs() < 1e-6);
    }
}
