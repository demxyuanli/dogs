//! `LocOpe_GluedShape` (`LocOpe_GluedShape.hxx:32-67`, `.cxx:38-226`): the
//! concrete [`GeneratedShape`] a `LocOpe_Gluer` hands to `LocOpe_Generator`.
//!
//! The shape keeps the feature shape (`myShape`) and the subset of its faces
//! that is glued onto the basis (`GlueOnFace`, `cxx:57-74`). `MapEdgeAndVertices`
//! (`cxx:78-179`) then fills the three caches the query methods answer from:
//!
//! - `myGShape` (`cxx:102-125`, `cxx:141-166`): an edge of a glued face whose
//!   two ancestor faces are not both glued *generates* a face, and the vertex of
//!   such an edge that also lies on another edge of that ancestor face
//!   *generates* that other edge;
//! - `myGEdges` (`cxx:121`): the generating edges, each stored `Reversed`;
//! - `myList` (`cxx:171-178`): the faces of the feature that are *not* glued.
//!
//! Deviations, both forced by the port's shape model:
//! - `myGShape.Bind(vtx, TopoDS_Edge())` (`cxx:151-153`) stores a *null* edge so
//!   a later generating edge does not rebind that vertex. A null shape has no
//!   `TopoShape` counterpart, so those vertices live in
//!   [`GluedData::null_edges`] and [`GeneratedShape::generated_edge`] answers
//!   `None` for them, which is what the reference's `TopoDS::Edge(null)` is.
//! - the two `Standard_ConstructionError` raises (`cxx:70-72`, `cxx:108-110`)
//!   have no shape sink: the port stores the message next to the caches and
//!   [`GluedShape::map_error`] lets the caller abort at the same point (before
//!   `LocOpe_Generator` runs) instead of mid-way through a raised exception.
//!
//! The reference recomputes `MapEdgeAndVertices` whenever `myGShape` is empty
//! (`cxx:183-186`, `cxx:193-197`, `cxx:205-209`, `cxx:215-219`). For one fixed
//! glue set the computation is a pure function of `myShape`, and an empty result
//! stays empty on every replay, so the port caches it once in a `OnceCell`.

use std::cell::OnceCell;

use crate::abs::ShapeType;
use crate::locope::generated_shape::GeneratedShape;
use crate::locope::shape_maps::{ShapeDataMap, ShapeSet};
use crate::shape::{Face, TopoShape};
use crate::topo_tools::{explorer_shapes, is_same, AncestorMap};

/// `Standard_ConstructionError` of `cxx:70-72` and `cxx:108-110`.
const CONSTRUCTION_ERROR: &str = "LocOpe_GluedShape: Standard_ConstructionError";

/// The three caches `MapEdgeAndVertices` fills (`hxx:62-65`).
#[derive(Debug, Clone, Default)]
pub struct GluedData {
    /// `myGShape`: vertices to their generated edge (null for
    /// [`Self::null_edges`]), edges to their generated face.
    gshape: ShapeDataMap,
    /// `myGEdges`: the generating edges, `Reversed`.
    gedges: Vec<TopoShape>,
    /// `myList`: the faces of the feature shape that are not glued.
    list: Vec<TopoShape>,
    /// Vertices the reference binds to a null edge (`cxx:151-153`).
    null_edges: ShapeSet,
}

impl GluedData {
    /// `myGShape` lookups behind `Generated(Vertex)` / `Generated(Edge)`.
    pub fn generated(&self, s: &TopoShape) -> Option<&TopoShape> {
        self.gshape.find(s)
    }

    /// The vertices bound to a null edge (`cxx:151-153`).
    pub fn is_null_edge(&self, v: &TopoShape) -> bool {
        self.null_edges.contains(v)
    }

    /// The generating edges, in the order `cxx:121` appends them.
    pub fn generating_edges(&self) -> &[TopoShape] {
        &self.gedges
    }

    /// The non-glued faces, in the order `cxx:175-177` appends them.
    pub fn oriented_faces(&self) -> &[TopoShape] {
        &self.list
    }
}

/// `LocOpe_GluedShape` (`hxx:32-67`).
#[derive(Debug, Default)]
pub struct GluedShape {
    /// `TopoDS_Shape myShape` (`hxx:55`).
    my_shape: Option<TopoShape>,
    /// `NCollection_Map<TopoDS_Shape> myMap` (`hxx:56`): the glued faces, with
    /// the orientation they carry in `myShape` (`cxx:74`).
    my_map: ShapeSet,
    /// `myGShape` / `myGEdges` / `myList` (`hxx:62-65`), computed once.
    data: OnceCell<Result<GluedData, String>>,
}

impl GluedShape {
    /// `LocOpe_GluedShape()` (`cxx:38`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `LocOpe_GluedShape(const TopoDS_Shape&)` (`cxx:42-45`).
    pub fn with_shape(s: &TopoShape) -> Self {
        let mut this = Self::new();
        this.init(s);
        this
    }

    /// `Init(const TopoDS_Shape&)` (`cxx:48-55`).
    pub fn init(&mut self, s: &TopoShape) {
        self.my_shape = Some(s.clone());
        self.my_map.clear();
        self.data = OnceCell::new();
    }

    /// `myShape`.
    pub fn shape(&self) -> Option<&TopoShape> {
        self.my_shape.as_ref()
    }

    /// `GlueOnFace(const TopoDS_Face&)` (`cxx:57-74`): record the glue face with
    /// the orientation it has inside `myShape`.
    ///
    /// `Err` is the `Standard_ConstructionError` throw of `cxx:70-72`.
    pub fn glue_on_face(&mut self, f: &Face) -> Result<(), String> {
        // `cxx:58-70`.
        let found = match self.my_shape.as_ref() {
            Some(shape) => explorer_shapes(shape, ShapeType::Face)
                .into_iter()
                .find(|cur| is_same(cur, &f.0)),
            None => None,
        };
        let Some(cur) = found else {
            return Err(CONSTRUCTION_ERROR.to_string());
        };
        // `cxx:74`: `myMap.Add(exp.Current())` keeps the explorer's orientation.
        self.my_map.add(&cur);
        Ok(())
    }

    /// The glued faces, as `myMap` holds them.
    pub fn glued_faces(&self) -> &[TopoShape] {
        self.my_map.iter()
    }

    /// `MapEdgeAndVertices()` (`cxx:78-179`), run once.
    ///
    /// `Err` is the `Standard_ConstructionError` of `cxx:108-110` (an edge of a
    /// glued face with other than two ancestor faces).
    pub fn map_edge_and_vertices(&self) -> Result<GluedData, String> {
        let mut out = GluedData::default();
        let Some(shape) = self.my_shape.as_ref() else {
            return Ok(out);
        };

        // `cxx:87-88`: `theMapEF = MapShapesAndAncestors(myShape, EDGE, FACE)`.
        let the_map_ef =
            AncestorMap::map_shapes_and_ancestors(shape, ShapeType::Edge, ShapeType::Face);
        // `cxx:90`: `mapdone`.
        let mut map_done = ShapeSet::new();

        // `cxx:96-128`.
        for fac in self.my_map.iter() {
            for edg in explorer_shapes(fac, ShapeType::Edge) {
                if map_done.contains(&edg) {
                    continue;
                }
                // `cxx:105-110`: `theMapEF.FindFromKey(edg).Extent() != 2`.
                let Some(ancestors) = the_map_ef.find_from_key(&edg) else {
                    return Err(CONSTRUCTION_ERROR.to_string());
                };
                if ancestors.len() != 2 {
                    return Err(CONSTRUCTION_ERROR.to_string());
                }
                // `cxx:111-119`: the first ancestor that is not glued.
                let bound = ancestors
                    .iter()
                    .find(|anc| !self.my_map.contains(anc))
                    .cloned();
                if let Some(face) = bound {
                    // `cxx:121-123`: `myGEdges.Append(edg.Reversed())`.
                    let mut reversed = edg.clone();
                    reversed.reverse();
                    out.gedges.push(reversed);
                    out.gshape.bind(&edg, &face);
                }
                map_done.add(&edg);
            }
        }

        // `cxx:131-168`.
        for edg in out.gedges.clone() {
            // `cxx:141`: `myGShape(edg)`, which the bind above guarantees.
            let Some(ancestor) = out.gshape.find(&edg).cloned() else {
                continue;
            };
            for vtx in explorer_shapes(&edg, ShapeType::Vertex) {
                // `cxx:138-141`.
                if out.gshape.is_bound(&vtx) || out.null_edges.contains(&vtx) {
                    continue;
                }
                // `cxx:142-167`: the first edge of the ancestor face that is not
                // `edg` and that shares the vertex stops both loops.
                for other in explorer_shapes(&ancestor, ShapeType::Edge) {
                    if is_same(&other, &edg) {
                        continue;
                    }
                    let mut hit = false;
                    for v3 in explorer_shapes(&other, ShapeType::Vertex) {
                        if is_same(&v3, &vtx) {
                            // `cxx:150-158`.
                            if out.gshape.is_bound(&other) {
                                out.null_edges.add(&vtx);
                            } else {
                                out.gshape.bind(&vtx, &other);
                            }
                            hit = true;
                            break;
                        }
                    }
                    if hit {
                        break;
                    }
                }
            }
        }

        // `cxx:171-178`: the faces that are not glued.
        for fac in explorer_shapes(shape, ShapeType::Face) {
            if !self.my_map.contains(&fac) {
                out.list.push(fac);
            }
        }

        Ok(out)
    }

    /// The cached `MapEdgeAndVertices` answer, `Err` for the throw.
    pub fn data(&self) -> Result<&GluedData, &String> {
        self.data
            .get_or_init(|| self.map_edge_and_vertices())
            .as_ref()
    }

    /// The `Standard_ConstructionError` `MapEdgeAndVertices` raised, if any.
    pub fn map_error(&self) -> Option<&String> {
        self.data().err()
    }
}

impl GeneratedShape for GluedShape {
    /// `GeneratingEdges()` (`cxx:182-190`).
    fn generating_edges(&self) -> &[TopoShape] {
        match self.data() {
            Ok(d) => d.generating_edges(),
            Err(_) => &[],
        }
    }

    /// `Generated(const TopoDS_Vertex&)` (`cxx:192-200`): `myGShape(V)`.
    fn generated_edge(&self, v: &TopoShape) -> Option<TopoShape> {
        let d = self.data().ok()?;
        if d.is_null_edge(v) {
            return None;
        }
        d.generated(v).cloned()
    }

    /// `Generated(const TopoDS_Edge&)` (`cxx:202-210`): `myGShape(E)`.
    fn generated_face(&self, e: &TopoShape) -> Option<TopoShape> {
        self.data().ok()?.generated(e).cloned()
    }

    /// `OrientedFaces()` (`cxx:212-225`).
    fn oriented_faces(&self) -> &[TopoShape] {
        match self.data() {
            Ok(d) => d.oriented_faces(),
            Err(_) => &[],
        }
    }
}
