//! Module 3: brep_scene — a named scene of transformed BRep shapes.
//!
//! A lightweight scene graph: each `SceneShape` holds a `TopoShape`, an
//! optional `GpTrsf`, and an optional color. Export meshes every shape (in its
//! placed position) and merges the results into one OBJ/STL/PLY file.

use std::collections::HashMap;

use occt_core::bnd::BndBox;
use occt_core::gp::GpTrsf;
use occt_core::io::obj::{ObjFace, ObjMesh};
use occt_core::io::ply::PlyMesh;
use occt_core::io::stl::StlMesh;
use occt_core::quantity::Color;

use crate::builder::TopoBuilder;
use crate::shape::{Compound, TopoShape};

/// A single named shape placed in a scene.
#[derive(Debug, Clone)]
pub struct SceneShape {
    pub name: String,
    pub shape: TopoShape,
    pub transform: Option<GpTrsf>,
    pub color: Option<Color>,
}

/// A named collection of shapes with optional transforms.
#[derive(Debug, Clone)]
pub struct BRepScene {
    pub shapes: Vec<SceneShape>,
    by_name: HashMap<String, usize>,
}

impl BRepScene {
    /// Create an empty scene.
    pub fn new() -> Self {
        Self { shapes: Vec::new(), by_name: HashMap::new() }
    }

    /// Add a shape with no transform or color.
    pub fn add(&mut self, name: impl Into<String>, shape: TopoShape) {
        self.add_inner(name.into(), shape, None, None);
    }

    /// Add a shape placed by `trsf`.
    pub fn add_transformed(&mut self, name: impl Into<String>, shape: TopoShape, trsf: &GpTrsf) {
        self.add_inner(name.into(), shape, Some(trsf.clone()), None);
    }

    fn add_inner(&mut self, name: String, shape: TopoShape, transform: Option<GpTrsf>, color: Option<Color>) {
        self.by_name.insert(name.clone(), self.shapes.len());
        self.shapes.push(SceneShape { name, shape, transform, color });
    }

    /// Look up a shape by name.
    pub fn find(&self, name: &str) -> Option<&SceneShape> {
        self.by_name.get(name).map(|&i| &self.shapes[i])
    }

    /// Number of shapes in the scene.
    pub fn len(&self) -> usize { self.shapes.len() }

    /// True when the scene holds no shapes.
    pub fn is_empty(&self) -> bool { self.shapes.is_empty() }

    /// Names of all shapes, in insertion order.
    pub fn names(&self) -> Vec<String> {
        self.shapes.iter().map(|s| s.name.clone()).collect()
    }

    /// Remove a shape by name, returning it.
    pub fn remove(&mut self, name: &str) -> Option<SceneShape> {
        let i = self.by_name.remove(name)?;
        let s = self.shapes.remove(i);
        self.reindex();
        Some(s)
    }

    fn reindex(&mut self) {
        self.by_name.clear();
        for (i, s) in self.shapes.iter().enumerate() {
            self.by_name.insert(s.name.clone(), i);
        }
    }

    /// Mesh the shape and place it by its scene transform.
    ///
    /// The mesher (`shape_mesh::mesh_shape`) emits vertices in the shape's own
    /// frame and does not apply the `TopoShape::location`, so we apply the
    /// scene transform to the mesh vertices ourselves.
    fn placed_mesh(s: &SceneShape, deflection: f64) -> crate::mesh::ShapeMesh {
        let mut mesh = crate::shape_mesh::mesh_shape(&s.shape, deflection);
        if let Some(t) = &s.transform {
            for v in &mut mesh.vertices {
                v.transform(t);
            }
        }
        mesh
    }

    /// Axis-aligned bounding box covering every placed shape.
    pub fn bounding_box(&self) -> BndBox {
        let mut result = BndBox::new();
        for s in &self.shapes {
            let mesh = crate::shape_mesh::mesh_shape(&s.shape, 0.1);
            let b = crate::mesh::mesh_bbox(&mesh);
            let b = match &s.transform {
                Some(t) => b.transformed(t),
                None => b,
            };
            result.add_box(&b);
        }
        result
    }

    /// All shapes merged into one Compound (scene transforms not baked in).
    pub fn merged_shape(&self) -> TopoShape {
        let builder = TopoBuilder::new();
        let mut comp = Compound::new();
        for s in &self.shapes {
            builder.add_compound(&mut comp, &s.shape);
        }
        comp.0
    }

    /// Export every shape (in its placed position) to one OBJ mesh.
    pub fn export_obj(&self, deflection: f64) -> String {
        occt_core::io::obj::write_obj(&self.merged_obj(deflection))
    }

    fn merged_obj(&self, deflection: f64) -> ObjMesh {
        let mut mesh = ObjMesh::default();
        for s in &self.shapes {
            let m = Self::placed_mesh(s, deflection);
            let base = mesh.vertices.len();
            mesh.vertices.extend_from_slice(&m.vertices);
            for t in &m.triangles {
                mesh.faces.push(ObjFace {
                    v: vec![(base + t.n0) as i32, (base + t.n1) as i32, (base + t.n2) as i32],
                    vt: None,
                    vn: None,
                });
            }
        }
        mesh
    }

    /// Export every shape (in its placed position) to binary STL.
    pub fn export_stl_binary(&self, deflection: f64) -> Vec<u8> {
        let mut stl = StlMesh::default();
        for s in &self.shapes {
            let m = Self::placed_mesh(s, deflection);
            for t in &m.triangles {
                stl.triangles.push([m.vertices[t.n0], m.vertices[t.n1], m.vertices[t.n2]]);
            }
        }
        occt_core::io::stl::write_binary_stl(&stl)
    }

    /// Export every shape (in its placed position) to ASCII PLY.
    pub fn export_ply(&self, deflection: f64) -> String {
        let mut mesh = PlyMesh::default();
        for s in &self.shapes {
            let m = Self::placed_mesh(s, deflection);
            let base = mesh.vertices.len();
            mesh.vertices.extend_from_slice(&m.vertices);
            for t in &m.triangles {
                mesh.faces.push(vec![base + t.n0, base + t.n1, base + t.n2]);
            }
        }
        occt_core::io::ply::write_ply(&mesh)
    }

    /// Write all shapes (placed) to an OBJ file.
    pub fn export_obj_file(&self, path: &str, deflection: f64) -> std::io::Result<()> {
        occt_core::io::obj::write_obj_file(path, &self.merged_obj(deflection))
    }
}

impl Default for BRepScene {
    fn default() -> Self { Self::new() }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep_tool::BRepTool;
    use crate::builder::TopoBuilder;
    use crate::shape::Vertex;
    use occt_core::gp::{GpAx1, GpAx3, GpDir, GpPln, GpPnt, GpVec};

    /// Build a unit cube [0,1]^3 as a Solid with real geometry.
    fn unit_box() -> TopoShape {
        let b = TopoBuilder::new();
        let mk = |x: f64, y: f64, z: f64| b.make_vertex(GpPnt::new(x, y, z), 0.0);
        let p0 = mk(0.,0.,0.); let p1 = mk(1.,0.,0.); let p2 = mk(1.,1.,0.); let p3 = mk(0.,1.,0.);
        let p4 = mk(0.,0.,1.); let p5 = mk(1.,0.,1.); let p6 = mk(1.,1.,1.); let p7 = mk(0.,1.,1.);

        let face = |c: [&Vertex; 4], normal: GpDir| {
            let e = [
                b.make_edge_segment(&BRepTool::vertex_point(c[0]), &BRepTool::vertex_point(c[1])),
                b.make_edge_segment(&BRepTool::vertex_point(c[1]), &BRepTool::vertex_point(c[2])),
                b.make_edge_segment(&BRepTool::vertex_point(c[2]), &BRepTool::vertex_point(c[3])),
                b.make_edge_segment(&BRepTool::vertex_point(c[3]), &BRepTool::vertex_point(c[0])),
            ];
            let wire = b.make_wire(&e);
            let pln = GpPln::new(GpAx3::from_ax1(&GpAx1::new(BRepTool::vertex_point(c[0]), normal)));
            let mut f = b.make_face_plane(&pln);
            b.add_wire(&mut f, &wire);
            f
        };

        let nz = GpDir::new(0.,0.,-1.).unwrap(); let pz = GpDir::new(0.,0.,1.).unwrap();
        let nx = GpDir::new(-1.,0.,0.).unwrap(); let px = GpDir::new(1.,0.,0.).unwrap();
        let ny = GpDir::new(0.,-1.,0.).unwrap(); let py = GpDir::new(0.,1.,0.).unwrap();

        let faces = vec![
            face([&p0, &p1, &p2, &p3], nz),
            face([&p4, &p5, &p6, &p7], pz),
            face([&p0, &p3, &p7, &p4], nx),
            face([&p1, &p5, &p6, &p2], px),
            face([&p0, &p1, &p5, &p4], ny),
            face([&p3, &p2, &p6, &p7], py),
        ];
        let shell = b.make_shell(&faces);
        let solid = b.make_solid(&[shell]);
        solid.into()
    }

    fn translation(dx: f64, dy: f64, dz: f64) -> GpTrsf {
        let mut t = GpTrsf::identity();
        t.set_translation_vec(&GpVec::new(dx, dy, dz));
        t
    }

    #[test]
    fn add_find_remove() {
        let mut scene = BRepScene::new();
        scene.add("a", unit_box());
        scene.add("b", unit_box());
        assert_eq!(scene.len(), 2);
        assert_eq!(scene.names(), vec!["a".to_string(), "b".to_string()]);
        assert!(scene.find("a").is_some());
        assert!(scene.find("missing").is_none());
        assert!(scene.remove("a").is_some());
        assert_eq!(scene.len(), 1);
        assert!(scene.find("a").is_none());
        assert!(scene.find("b").is_some());
    }

    #[test]
    fn export_places_and_bboxes_all_shapes() {
        let mut scene = BRepScene::new();
        scene.add("box1", unit_box());
        scene.add_transformed("box2", unit_box(), &translation(2.0, 0.0, 0.0));

        let obj = scene.export_obj(0.1);
        assert!(obj.contains("v 0.000000000 0.000000000 0.000000000"), "first box at origin");
        assert!(obj.contains("v 2.000000000 0.000000000 0.000000000"), "second box translated");

        let bb = scene.bounding_box();
        let (xmin, xmax, ymin, ymax, zmin, zmax) = bb.get().expect("finite bounds");
        assert!(xmin <= 0.0 && xmax >= 3.0, "x spans both boxes: {xmin}..{xmax}");
        assert!(ymin <= 0.0 && ymax >= 1.0, "y spans box depth");
        assert!(zmin <= 0.0 && zmax >= 1.0, "z spans box height");

        assert!(scene.export_stl_binary(0.1).len() >= 84);
        assert!(scene.export_ply(0.1).starts_with("ply"));
    }

    #[test]
    fn merged_shape_compound() {
        let mut scene = BRepScene::new();
        scene.add("a", unit_box());
        scene.add("b", unit_box());
        let merged = scene.merged_shape();
        assert!(merged.is_compound());
        assert_eq!(merged.tshape.read().unwrap().children.len(), 2);
    }
}
