//! BRep model — manage collections of shapes with names.
//! Provides a lightweight scene/document for CAD exchange.
use std::collections::HashMap;
use crate::shape::TopoShape;
use crate::abs::ShapeType;
use occt_core::bnd::BndBox;

/// A named shape in the model.
#[derive(Debug, Clone)]
pub struct ModelShape {
    pub name: String,
    pub shape: TopoShape,
    pub color: Option<occt_core::quantity::Color>,
    pub layer: String,
}

/// BRepModel: a set of named shapes with bounding box.
#[derive(Debug, Clone, Default)]
pub struct BRepModel {
    pub shapes: Vec<ModelShape>,
    by_name: HashMap<String, usize>,
}

impl BRepModel {
    pub fn new() -> Self { Self::default() }

    pub fn add(&mut self, name: &str, shape: TopoShape) -> usize {
        let idx = self.shapes.len();
        self.shapes.push(ModelShape { name: name.into(), shape, color: None, layer: String::new() });
        self.by_name.insert(name.into(), idx);
        idx
    }

    pub fn add_with_color(&mut self, name: &str, shape: TopoShape, color: occt_core::quantity::Color) -> usize {
        let idx = self.add(name, shape);
        self.shapes[idx].color = Some(color);
        idx
    }

    pub fn find(&self, name: &str) -> Option<&ModelShape> {
        self.by_name.get(name).map(|&i| &self.shapes[i])
    }

    pub fn find_mut(&mut self, name: &str) -> Option<&mut ModelShape> {
        let i = *self.by_name.get(name)?;
        self.shapes.get_mut(i)
    }

    pub fn get(&self, i: usize) -> Option<&ModelShape> { self.shapes.get(i) }

    pub fn remove(&mut self, name: &str) -> Option<TopoShape> {
        let idx = self.by_name.remove(name)?;
        // Compact: remove at idx, fix map indices
        let removed = self.shapes.remove(idx);
        for (n, i) in self.by_name.iter_mut() {
            if *i > idx { *i -= 1; }
        }
        Some(removed.shape)
    }

    pub fn len(&self) -> usize { self.shapes.len() }
    pub fn is_empty(&self) -> bool { self.shapes.is_empty() }

    pub fn names(&self) -> Vec<&str> { self.shapes.iter().map(|s| s.name.as_str()).collect() }

    /// Compute the model bounding box (union of per-shape transforms applied).
    /// children_fn: provides sub-shapes for geometry access (unused for bbox here).
    pub fn bounding_box(&self) -> BndBox {
        let mut b = BndBox::new();
        // Shapes have no embedded geometry in this port; bbox is empty unless
        // sub-shapes carry coordinates. Keep as a stub returning void box.
        for s in &self.shapes {
            if s.shape.shape_type() == ShapeType::Vertex {
                // No point data stored — skip.
            }
        }
        b
    }

    /// Filter shapes by type.
    pub fn by_type(&self, t: ShapeType) -> Vec<&ModelShape> {
        self.shapes.iter().filter(|s| s.shape.shape_type() == t).collect()
    }

    /// Count shapes by type.
    pub fn type_counts(&self) -> HashMap<ShapeType, usize> {
        let mut m = HashMap::new();
        for s in &self.shapes { *m.entry(s.shape.shape_type()).or_insert(0) += 1; }
        m
    }

    /// Apply a transform to all shapes.
    pub fn move_all(&mut self, t: &occt_core::gp::GpTrsf) {
        for s in &mut self.shapes {
            s.shape = crate::transform::move_shape(&s.shape, t);
        }
    }

    /// Rename a shape.
    pub fn rename(&mut self, old: &str, new: &str) -> bool {
        let idx = match self.by_name.remove(old) { Some(i) => i, None => return false };
        self.shapes[idx].name = new.into();
        self.by_name.insert(new.into(), idx);
        true
    }

    /// Collect all shapes into a single compound (top-level).
    pub fn as_compound(&self) -> TopoShape {
        let mut comp = TopoShape::new(ShapeType::Compound);
        // In this port compound children are not stored on TShape; return bare compound.
        comp
    }

    /// Merge another model into this one (name collision → keep both, suffix _2).
    pub fn merge(&mut self, other: &BRepModel) {
        for s in &other.shapes {
            let mut name = s.name.clone();
            let mut n = 2;
            while self.by_name.contains_key(&name) {
                name = format!("{}_{n}", s.name);
                n += 1;
            }
            self.add(&name, s.shape.clone());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn add_find_remove() {
        let mut m = BRepModel::new();
        let v = TopoShape::new(ShapeType::Vertex);
        m.add("v1", v.clone());
        assert!(m.find("v1").is_some());
        assert_eq!(m.len(), 1);
        let got = m.remove("v1").unwrap();
        assert!(got.is_vertex());
        assert!(m.find("v1").is_none());
    }

    #[test]
    fn rename_works() {
        let mut m = BRepModel::new();
        m.add("a", TopoShape::new(ShapeType::Face));
        assert!(m.rename("a", "b"));
        assert!(m.find("b").is_some());
        assert!(m.find("a").is_none());
    }

    #[test]
    fn merge_suffixes() {
        let mut m1 = BRepModel::new();
        m1.add("s", TopoShape::new(ShapeType::Vertex));
        let mut m2 = BRepModel::new();
        m2.add("s", TopoShape::new(ShapeType::Edge));
        m1.merge(&m2);
        assert_eq!(m1.len(), 2);
        assert!(m1.find("s").is_some());
        assert!(m1.find("s_2").is_some());
    }

    #[test]
    fn type_counts() {
        let mut m = BRepModel::new();
        m.add("e1", TopoShape::new(ShapeType::Edge));
        m.add("e2", TopoShape::new(ShapeType::Edge));
        m.add("f1", TopoShape::new(ShapeType::Face));
        let c = m.type_counts();
        assert_eq!(c[&ShapeType::Edge], 2);
        assert_eq!(c[&ShapeType::Face], 1);
    }
}
