//! `BRepOffset_Tool::Deboucle3D` (`BRepOffset_Tool.cxx:4337-4418`): drop the
//! shells that carry a free border which is neither a free border of a cap
//! (`Boundary`) nor a degenerated edge.
//!
//! A compound or a solid keeps only the recursively non-empty children
//! (`cxx:4381-4411`); every other shape type is dropped (`cxx:4413-4414`),
//! which is why the port answers `None` where OCCT returns a null shape.
//!
//! The `Boundary` argument is `NCollection_Map<TopoDS_Shape,
//! TopTools_ShapeMapHasher>` (`cxx:4339`), i.e. a `Contains` test on the
//! `IsSame` key ([`ShapeMapKey`]).

use std::collections::HashSet;

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, TopoShape};
use crate::topo_tools::{
    children_keep_orientation, explorer_shapes, map_shapes_and_ancestors, shape_map_key, ShapeMapKey,
};

impl crate::brep_offset::BRepOffsetTool {
    /// `Deboucle3D(S, Boundary)` (`cxx:4337-4418`). `None` is the null shape
    /// the reference returns for a discarded shell (and for every shape type
    /// the `default:` arm drops).
    pub fn deboucle_3d(
        s: &TopoShape,
        boundary: &HashSet<ShapeMapKey>,
    ) -> Option<TopoShape> {
        match s.shape_type() {
            // `cxx:4344-4379`.
            ShapeType::Shell => {
                // `cxx:4347-4351`: `TopExp::MapShapesAndAncestors(S, EDGE,
                // FACE, Map)`.
                let map = map_shapes_and_ancestors(s, ShapeType::Edge, ShapeType::Face);

                // `cxx:4353-4373`. The reference walks the map (`Extent()`
                // order); `JeGarde` only ever turns `false`, so the walk order
                // is not observable and the port walks the explorer, which is
                // the same set of edge occurrences with a defined order.
                let mut je_garde = true;
                for shape in explorer_shapes(s, ShapeType::Edge) {
                    if !je_garde {
                        break;
                    }
                    let edge = Edge(shape);
                    // `cxx:4356-4357`.
                    let faces = map
                        .get(&shape_map_key(&edge.0))
                        .cloned()
                        .unwrap_or_default();
                    if faces.len() < 2 {
                        // `cxx:4360-4367`: an `INTERNAL` edge whose only
                        // ancestor face is not `INTERNAL` is skipped.
                        if edge.0.orientation() == Orientation::Internal {
                            let Some(a_face) = faces.first() else {
                                continue;
                            };
                            if a_face.orientation() != Orientation::Internal {
                                continue;
                            }
                        }
                        // `cxx:4368-4371`.
                        if !boundary.contains(&shape_map_key(&edge.0))
                            && !BRepTool::is_degenerated(&edge)
                        {
                            je_garde = false;
                        }
                    }
                }
                // `cxx:4374-4377`.
                if je_garde {
                    Some(s.clone())
                } else {
                    None
                }
            }

            // `cxx:4381-4411`.
            ShapeType::Compound | ShapeType::Solid => {
                let builder = TopoBuilder::new();
                let is_compound = s.shape_type() == ShapeType::Compound;
                let mut ss = builder.make_shape(if is_compound {
                    ShapeType::Compound
                } else {
                    ShapeType::Solid
                });
                let mut nb_sub = 0;
                for cur_s in children_keep_orientation(s) {
                    // `cxx:4399-4404`.
                    if let Some(sub_shape) = Self::deboucle_3d(&cur_s, boundary) {
                        builder.add(&mut ss, &sub_shape);
                        nb_sub += 1;
                    }
                }
                // `cxx:4406-4409`.
                if nb_sub == 0 {
                    None
                } else {
                    Some(ss)
                }
            }

            // `cxx:4413-4414`.
            _ => None,
        }
    }
}
