//! `TopOpeBRepTool_HBoxTool` — the shape-to-box registry the classic boolean
//! engine sorts against.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepTool/TopOpeBRepTool_HBoxTool.hxx:24-52`,
//! `.cxx:51-260`.
//!
//! Translation notes:
//! - `NCollection_IndexedDataMap<TopoDS_Shape, Bnd_Box>` (`hxx:47`) becomes
//!   two parallel vectors (shapes, boxes) plus a key lookup, which is
//!   the same insertion-ordered, 1-based-indexed container.
//! - `myIMS` is keyed by `TopTools_ShapeMapHasher` (TShape *and* Location), so
//!   the port keys on [`crate::bop_occt_util::shape_key`].
//! - `BRepBndLib::Add(S, B)` uses its `useTriangulation = true` default
//!   (`BRepBndLib.hxx:40-42`).
//! - `Box(S)` is non-const in OCCT (`cxx:160`) only because `myIMS.FindFromKey`
//!   is; it is a pure read in the port.
//! - `HBT::Box1` / `Box2` / `Box4` / `ComputeBox : invalid type` raise
//!   `Standard_ProgramError` in OCCT; the port returns `Err` with those exact
//!   messages.

use std::collections::HashMap;

use occt_core::bnd::box3d::BndBox;
use occt_core::gp::GpPnt;
use occt_core::precision::PCONFUSION;

use crate::abs::ShapeType;
use crate::bop_occt_util::shape_key;
use crate::brep_bnd_lib;
use crate::brep_tool::BRepTool;
use crate::shape::TopoShape;
use crate::topo_tools::explorer_shapes_avoiding;
use crate::topo_tools::ShapeMapKey;

/// `class TopOpeBRepTool_HBoxTool` (`hxx:24-52`).
#[derive(Default)]
pub struct HBoxTool {
    /// `myIMS` shapes.
    shapes: Vec<TopoShape>,
    /// `myIMS` boxes.
    boxes: Vec<BndBox>,
    /// `myIMS` index lookup (1-based).
    index: HashMap<ShapeMapKey, usize>,
}

impl HBoxTool {
    /// `TopOpeBRepTool_HBoxTool()` (`cxx:51`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `Clear()` (`cxx:57-60`).
    pub fn clear(&mut self) {
        self.shapes.clear();
        self.boxes.clear();
        self.index.clear();
    }

    /// `AddBoxes(S, TS, TA)` (`cxx:63-74`): add a box for every `TS` sub-shape
    /// of `S` not below a `TA` one.
    pub fn add_boxes(&mut self, s: &TopoShape, ts: ShapeType, ta: ShapeType) -> Result<(), String> {
        for ss in explorer_shapes_avoiding(s, ts, ta) {
            self.add_box(&ss)?;
        }
        Ok(())
    }

    /// `AddBox(S)` (`cxx:78-105`).
    pub fn add_box(&mut self, s: &TopoShape) -> Result<(), String> {
        if self.has_box(s) {
            return Ok(());
        }
        let mut b = BndBox::new();
        Self::compute_box(s, &mut b)?;
        let key = shape_key(s);
        self.shapes.push(s.clone());
        self.boxes.push(b);
        self.index.insert(key, self.boxes.len());
        Ok(())
    }

    /// `ComputeBox(S, B)` (`cxx:107-135`): `BRepBndLib::Add` for `FACE`,
    /// `EDGE`, `SOLID`, `COMPOUND` and `VERTEX`; every other type raises
    /// `HBT::ComputeBox : invalid type`.
    pub fn compute_box(s: &TopoShape, b: &mut BndBox) -> Result<(), String> {
        match s.shape_type() {
            ShapeType::Face
            | ShapeType::Edge
            | ShapeType::Solid
            | ShapeType::Compound
            | ShapeType::Vertex => {
                brep_bnd_lib::add(s, b, true);
                Ok(())
            }
            _ => Err("HBT::ComputeBox : invalid type".to_string()),
        }
    }

    /// `ComputeBoxOnVertices(S, B)` (`cxx:139-158`): the box of the edge's
    /// vertices, enlarged by the largest vertex tolerance; with no vertex the
    /// box is the whole `[-1e5, 1e5]` cube.
    pub fn compute_box_on_vertices(s: &TopoShape, b: &mut BndBox) {
        let vertices = explorer_shapes_avoiding(s, ShapeType::Vertex, ShapeType::Shape);
        if vertices.is_empty() {
            b.update(-1.0e5, -1.0e5, -1.0e5, 1.0e5, 1.0e5, 1.0e5);
            return;
        }
        let mut tol = PCONFUSION;
        for v in &vertices {
            let Some(vertex) = crate::Vertex::wrap(v.clone()) else {
                continue;
            };
            let p: GpPnt = BRepTool::vertex_point(&vertex);
            b.update(p.x(), p.y(), p.z(), p.x(), p.y(), p.z());
            tol = tol.max(BRepTool::vertex_tolerance(&vertex));
        }
        b.enlarge(tol);
    }

    /// `Box(S)` (`cxx:160-170`): `HBT::Box1` when no box was computed for `S`.
    pub fn box_of(&self, s: &TopoShape) -> Result<&BndBox, String> {
        let key = shape_key(s);
        let i = *self.index.get(&key).ok_or_else(|| "HBT::Box1".to_string())?;
        Ok(&self.boxes[i - 1])
    }

    /// `Box(I)` (`cxx:174-185`): `HBT::Box2` outside `[1, Extent()]`.
    pub fn box_at(&self, i: i32) -> Result<&BndBox, String> {
        let iu = self.extent();
        if i < 1 || i > iu {
            return Err("HBT::Box2".to_string());
        }
        Ok(&self.boxes[(i - 1) as usize])
    }

    /// `HasBox(S)` (`cxx:189-192`).
    pub fn has_box(&self, s: &TopoShape) -> bool {
        self.index.contains_key(&shape_key(s))
    }

    /// `Shape(I)` (`cxx:196-207`): `HBT::Box4` outside `[1, Extent()]`.
    pub fn shape(&self, i: i32) -> Result<&TopoShape, String> {
        let iu = self.extent();
        if i < 1 || i > iu {
            return Err("HBT::Box4".to_string());
        }
        Ok(&self.shapes[(i - 1) as usize])
    }

    /// `Index(S)` (`cxx:211-215`): `FindIndex` answers 0 for a missing key.
    pub fn index_of(&self, s: &TopoShape) -> i32 {
        match self.index.get(&shape_key(s)) {
            Some(&i) => i as i32,
            None => 0,
        }
    }

    /// `Extent()` (`cxx:219-223`).
    pub fn extent(&self) -> i32 {
        self.boxes.len() as i32
    }

    /// `IMS()` (`cxx:235-238`): the insertion-ordered entries.
    pub fn ims(&self) -> impl Iterator<Item = (&TopoShape, &BndBox)> {
        self.shapes.iter().zip(self.boxes.iter())
    }
}
