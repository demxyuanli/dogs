//! The four preparation passes of `BRepTools_Modifier::Perform`
//! (`BRepTools_Modifier.cxx:118-127`), each one a direct port:
//! `CreateNewVertices` (`cxx:651-678`), `FillNewCurveInfo` (`cxx:680-703`),
//! `FillNewSurfaceInfo` (`cxx:705-762`) and `CreateOtherVertices`
//! (`cxx:764-832`).

use crate::abs::{Orientation, ShapeType};
use crate::builder::TopoBuilder;
use crate::shape::TopoShape;
use crate::tgeometry::{GeometryRegistry, VertexGeom};
use crate::topo_tools::{explorer_shapes, map_shapes, shape_map_key};

use super::modification::Modification;
use super::modifier::Modifier;

impl Modifier {
    /// `CreateNewVertices` (`cxx:651-678`): bind a fresh vertex for every
    /// vertex the modification gives a new point to. With `myMutableInput` the
    /// remaining vertices are bound to themselves (forward oriented) so
    /// `CreateOtherVertices` does not run.
    pub(super) fn create_new_vertices(
        &mut self,
        mve: &[(TopoShape, Vec<TopoShape>)],
        m: &dyn Modification,
    ) {
        let b = TopoBuilder::new();
        for (v, _) in mve {
            match m.new_point(v) {
                Some(info) => {
                    let new_v = b.make_vertex(info.point, info.tolerance);
                    set_shape_flags(v, &new_v.0);
                    self.my_map.insert(shape_map_key(v), Some(new_v.0));
                    self.my_has_new_geom.insert(shape_map_key(v));
                }
                None => {
                    if self.my_mutable_input {
                        let forward = v.oriented(Orientation::Forward);
                        self.my_map.insert(shape_map_key(v), Some(forward));
                    }
                }
            }
        }
    }

    /// `FillNewCurveInfo` (`cxx:680-703`): record the `NewCurve` answer of
    /// every edge reachable under a face (`aMEF` is
    /// `MapShapesAndAncestors(myShape, EDGE, FACE)`, but every edge of the
    /// shape is a key of it, `TopExp.cxx:110-119`).
    pub(super) fn fill_new_curve_info(
        &mut self,
        mef: &[(TopoShape, Vec<TopoShape>)],
        m: &dyn Modification,
    ) {
        for (e, _) in mef {
            if let Some(info) = m.new_curve(e) {
                let key = shape_map_key(e);
                self.my_nc_info.insert(key.clone(), info);
                self.my_has_new_geom.insert(key);
            }
        }
    }

    /// `FillNewSurfaceInfo` (`cxx:705-762`).
    ///
    /// A face the modification leaves alone is still marked "not going to be
    /// modified" (`myNonUpdFace`) when neither its edges have a `NewCurve`
    /// entry nor any of their vertices has been bound to a new shape
    /// (`cxx:729-760`); `CreateOtherVertices` uses that to skip `NewCurve2d`
    /// calls the rebuild would never make.
    pub(super) fn fill_new_surface_info(&mut self, m: &dyn Modification) {
        let Some(shape) = self.my_shape.clone() else { return };
        // `TopExp::MapShapes(myShape, TopAbs_FACE, aMF)` (`cxx:708`).
        let faces = map_shapes(&shape, &[ShapeType::Face]);
        let faces = faces.get(&ShapeType::Face).cloned().unwrap_or_default();
        for f in &faces {
            if let Some(info) = m.new_surface(f) {
                let key = shape_map_key(f);
                self.my_ns_info.insert(key.clone(), info);
                self.my_has_new_geom.insert(key);
                continue;
            }
            // `cxx:731-759`: check whether the sub-shapes will be modified.
            let mut not_rebuilt = true;
            for e in explorer_shapes(f, ShapeType::Edge) {
                if self.my_nc_info.contains_key(&shape_map_key(&e)) {
                    not_rebuilt = false;
                    break;
                }
                for v in explorer_shapes(&e, ShapeType::Vertex) {
                    // `!myMap(aVV).IsNull()`: a bound non-null entry. `Put`
                    // bound every vertex of the tree, so the reference's
                    // implicit "insert a null entry" side effect of
                    // `DataMap::operator()` never fires here.
                    if matches!(self.my_map.get(&shape_map_key(&v)), Some(Some(_))) {
                        not_rebuilt = false;
                        break;
                    }
                }
                if !not_rebuilt {
                    break;
                }
            }
            if not_rebuilt {
                self.my_non_upd_face.insert(shape_map_key(f));
            }
        }
    }

    /// `CreateOtherVertices` (`cxx:764-832`): bind every vertex that is still
    /// unbound, sharing the original vertex when nothing below it changes and
    /// taking an `EmptyCopied` vertex otherwise.
    pub(super) fn create_other_vertices(
        &mut self,
        mve: &[(TopoShape, Vec<TopoShape>)],
        mef: &[(TopoShape, Vec<TopoShape>)],
        m: &dyn Modification,
    ) {
        for (v, edges) in mve {
            if self.mapped(v).is_some() {
                continue;
            }
            let mut to_replace = false;
            for e in edges {
                // `cxx:793-796`.
                if self
                    .my_nc_info
                    .get(&shape_map_key(e))
                    .is_some_and(|info| info.curve.is_some())
                {
                    to_replace = true;
                }
                if to_replace {
                    break;
                }
                // `cxx:800-817`: `aMEF.FindFromKey(anE)`.
                let faces = mef
                    .iter()
                    .find(|(s, _)| shape_map_key(s) == shape_map_key(e))
                    .map(|(_, list)| list.clone())
                    .unwrap_or_default();
                for f in &faces {
                    if self.my_non_upd_face.contains(&shape_map_key(f)) {
                        continue;
                    }
                    // `cxx:808-811`: the reference hands a null `NewF` here
                    // because "some NewCurve2d()s may use NewE arg internally".
                    // This port has no null shape, so `f` itself is passed for
                    // `NewF`; no in-tree modification reads it
                    // (`BRepTools_TrsfModification.cxx:337-341` ignores both
                    // `NewE` and `NewF`).
                    let dummy = crate::shape_fix_compose_shell::empty_copied_edge(e);
                    if m.new_curve2d(e, f, &dummy, f).is_some() {
                        to_replace = true;
                        break;
                    }
                }
                if to_replace {
                    break;
                }
            }
            let mut new_v = if to_replace {
                empty_copied_vertex(v)
            } else {
                v.clone()
            };
            new_v.set_orientation(Orientation::Forward);
            self.my_map.insert(shape_map_key(v), Some(new_v));
        }
    }
}

/// `SetShapeFlags` (`BRepTools_Modifier.cxx:834-842`): carry the source
/// shape's flags over to the rebuilt one.
///
/// UNPORTED: `Checked`, `Orientable` and `Convex` (`cxx:837`, `:838`, `:841`) —
/// `ShapeFlags` (`abs.rs`) only carries `Free`, `Closed`, `Infinite` and
/// `Modified`.
fn set_shape_flags(src: &TopoShape, dst: &TopoShape) {
    let (modified, closed, infinite) = {
        let ts = src.tshape.read().expect("poisoned TShape lock");
        (ts.modified(), ts.is_closed(), ts.infinite())
    };
    let mut ts = dst.tshape.write().expect("poisoned TShape lock");
    ts.set_modified(modified);
    ts.set_closed(closed);
    ts.set_infinite(infinite);
}

/// `TopoDS::Vertex(aV.EmptyCopied())` (`cxx:822`) — `TopoDS_Shape::EmptyCopied`
/// (`TopoDS_Shape.hxx:294-302`) plus `BRep_TVertex::EmptyCopy`
/// (`BRep_TVertex.cxx:33-39`), which keeps the point and the tolerance.
fn empty_copied_vertex(src: &TopoShape) -> TopoShape {
    let copy = src.empty_copied();
    if let Some(geom) = GeometryRegistry::global().vertex_geom(src) {
        GeometryRegistry::global().set_vertex(&copy, VertexGeom { point: geom.point, tolerance: geom.tolerance });
    }
    copy
}
