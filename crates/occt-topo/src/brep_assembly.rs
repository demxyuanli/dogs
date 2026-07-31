//! BRep assemblies — positioned parts with rigid transforms, joint offsets
//! and aggregate export.
//! Source: `XCAFDoc_ShapeTool`, `TopoDS_Shape` location composition.
//!
//! An assembly is a list of (name, shape, transform). Each part's transform
//! is stored as a `GpTrsf`; geometry is kept in the shape's local frame and
//! the transform is applied at export time (mirroring OCCT's Location-on-
//! shape model). Parts may be nested: `add_subassembly` composes transforms.

use std::collections::HashMap;

use occt_core::gp::{GpPnt, GpTrsf, GpVec};

use crate::brep_tool::BRepTool;
use crate::shape::TopoShape;

/// A positioned part within an assembly.
#[derive(Debug, Clone)]
pub struct AssemblyPart {
    pub name: String,
    pub shape: TopoShape,
    pub transform: GpTrsf,
    pub parent: Option<usize>,
}

/// A hierarchical assembly of named, transformed parts.
#[derive(Debug, Clone, Default)]
pub struct Assembly {
    pub parts: Vec<AssemblyPart>,
    by_name: HashMap<String, usize>,
}

impl Assembly {
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a part with the identity transform.
    pub fn add(&mut self, name: &str, shape: TopoShape) -> usize {
        self.add_transformed(name, shape, &GpTrsf::identity())
    }

    /// Add a part with an explicit rigid transform.
    pub fn add_transformed(&mut self, name: &str, shape: TopoShape, t: &GpTrsf) -> usize {
        let idx = self.parts.len();
        self.parts.push(AssemblyPart {
            name: name.into(),
            shape,
            transform: t.clone(),
            parent: None,
        });
        self.by_name.insert(name.into(), idx);
        idx
    }

    /// Add a part under a named parent; the transform is composed with the
    /// parent's accumulated transform.
    pub fn add_child(&mut self, parent: &str, name: &str, shape: TopoShape, local: &GpTrsf) -> Result<usize, String> {
        let pidx = *self.by_name.get(parent).ok_or("add_child: parent not found")?;
        let world = self.accumulated(pidx).multiplied(local);
        let idx = self.parts.len();
        self.parts.push(AssemblyPart {
            name: name.into(),
            shape,
            transform: world,
            parent: Some(pidx),
        });
        self.by_name.insert(name.into(), idx);
        Ok(idx)
    }

    pub fn len(&self) -> usize {
        self.parts.len()
    }

    pub fn is_empty(&self) -> bool {
        self.parts.is_empty()
    }

    pub fn find(&self, name: &str) -> Option<&AssemblyPart> {
        self.by_name.get(name).map(|&i| &self.parts[i])
    }

    /// Accumulated world transform of a part (composes all ancestors).
    pub fn accumulated(&self, idx: usize) -> GpTrsf {
        let mut t = self.parts[idx].transform.clone();
        let mut cur = self.parts[idx].parent;
        while let Some(p) = cur {
            t = self.parts[p].transform.multiplied(&t);
            cur = self.parts[p].parent;
        }
        t
    }

    /// The transformed world geometry of a part: every vertex point moved by
    /// the accumulated transform (via `BRepTool::vertex_point_world`-style
    /// application on a `shape_ops::transformed_copy`).
    pub fn world_vertices(&self, idx: usize) -> Vec<GpPnt> {
        use crate::shape_ops::transformed_copy;
        let t = self.accumulated(idx);
        match transformed_copy(&self.parts[idx].shape, &t) {
            Ok(copy) => crate::topo_tools_full::vertices_of(&copy)
                .iter()
                .map(|v| BRepTool::vertex_point(v))
                .collect(),
            Err(_) => Vec::new(),
        }
    }

    /// World bounding box of a part (min/max of its transformed vertices).
    pub fn part_bbox(&self, idx: usize) -> Option<(GpPnt, GpPnt)> {
        let verts = self.world_vertices(idx);
        if verts.is_empty() {
            return None;
        }
        let mut min = GpPnt::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
        let mut max = GpPnt::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
        for p in &verts {
            min = GpPnt::new(min.x().min(p.x()), min.y().min(p.y()), min.z().min(p.z()));
            max = GpPnt::new(max.x().max(p.x()), max.y().max(p.y()), max.z().max(p.z()));
        }
        Some((min, max))
    }

    /// Union bounding box of the whole assembly.
    pub fn assembly_bbox(&self) -> Option<(GpPnt, GpPnt)> {
        let mut bbox: Option<(GpPnt, GpPnt)> = None;
        for i in 0..self.parts.len() {
            if let Some(bb) = self.part_bbox(i) {
                bbox = Some(match bbox {
                    None => bb,
                    Some((mn, mx)) => (
                        GpPnt::new(mn.x().min(bb.0.x()), mn.y().min(bb.0.y()), mn.z().min(bb.0.z())),
                        GpPnt::new(mx.x().max(bb.1.x()), mx.y().max(bb.1.y()), mx.z().max(bb.1.z())),
                    ),
                });
            }
        }
        bbox
    }

    /// Export every part's mesh to one OBJ text (applying world transforms).
    pub fn export_obj(&self, deflection: f64) -> String {
        use crate::shape_ops::transformed_copy;
        let mut out = String::from("# Assembly OBJ export\n");
        for (i, part) in self.parts.iter().enumerate() {
            let t = self.accumulated(i);
            let copy = match transformed_copy(&part.shape, &t) {
                Ok(c) => c,
                Err(_) => continue,
            };
            let mesh = crate::shape_mesh::mesh_shape(&copy, deflection);
            out.push_str(&format!("# part: {}\n", part.name));
            for v in &mesh.vertices {
                out.push_str(&format!("v {:.9} {:.9} {:.9}\n", v.x(), v.y(), v.z()));
            }
            let base = 1 + out.matches("\nv ").count() as i32 - mesh.vertices.len() as i32;
            for tri in &mesh.triangles {
                out.push_str(&format!(
                    "f {} {} {}\n",
                    base + tri.n0 as i32,
                    base + tri.n1 as i32,
                    base + tri.n2 as i32
                ));
            }
        }
        out
    }

    /// Export the whole assembly to STEP as a compound of transformed copies.
    pub fn export_step(&self) -> String {
        use crate::shape_ops::transformed_copy;
        let mut comp = TopoShape::new(crate::abs::ShapeType::Compound);
        let builder = crate::builder::TopoBuilder::new();
        for i in 0..self.parts.len() {
            let t = self.accumulated(i);
            if let Ok(copy) = transformed_copy(&self.parts[i].shape, &t) {
                builder.add(&mut comp, &copy);
            }
        }
        crate::step::write_shape_step(&comp)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::BRepPrimBox;

    fn approx(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-6 * b.abs().max(1.0)
    }

    #[test]
    fn add_and_find() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let mut a = Assembly::new();
        a.add("base", b.solid.0.clone());
        a.add_transformed("block", b.solid.0, &{
            let mut t = GpTrsf::identity();
            t.set_translation_vec(&GpVec::new(2.0, 0.0, 0.0));
            t
        });
        assert_eq!(a.len(), 2);
        assert!(a.find("base").is_some());
        assert!(a.find("none").is_none());
    }

    #[test]
    fn world_bbox_translates() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let mut a = Assembly::new();
        a.add_transformed("moved", b.solid.0, &{
            let mut t = GpTrsf::identity();
            t.set_translation_vec(&GpVec::new(5.0, 0.0, 0.0));
            t
        });
        let (mn, mx) = a.part_bbox(0).expect("bbox");
        assert!(approx(mn.x(), 5.0), "min x {}", mn.x());
        assert!(approx(mx.x(), 6.0), "max x {}", mx.x());
    }

    #[test]
    fn child_composes_transform() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let mut a = Assembly::new();
        a.add("root", b.solid.0.clone());
        let mut t = GpTrsf::identity();
        t.set_translation_vec(&GpVec::new(2.0, 0.0, 0.0));
        a.add_child("root", "sub", b.solid.0, &t).expect("child");
        // sub's world offset: root (identity) composed with local (2,0,0) → 2.
        let (mn, _) = a.part_bbox(1).expect("sub bbox");
        assert!(approx(mn.x(), 2.0), "composed x {}", mn.x());
    }

    #[test]
    fn assembly_bbox_union() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let mut a = Assembly::new();
        a.add("a", b.solid.0.clone());
        a.add_transformed("b", b.solid.0, &{
            let mut t = GpTrsf::identity();
            t.set_translation_vec(&GpVec::new(10.0, 0.0, 0.0));
            t
        });
        let (mn, mx) = a.assembly_bbox().expect("union bbox");
        assert!(approx(mn.x(), 0.0));
        assert!(approx(mx.x(), 11.0));
    }

    #[test]
    fn export_obj_and_step() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let mut a = Assembly::new();
        a.add("one", b.solid.0.clone());
        let obj = a.export_obj(0.2);
        assert!(obj.contains("v "));
        assert!(obj.contains("f "));
        let step = a.export_step();
        assert!(step.starts_with("ISO-10303-21;"));
        assert!(step.contains("MANIFOLD_SOLID_BREP"));
    }
}
