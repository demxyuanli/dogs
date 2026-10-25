//! `BRepTools_Substitution` — replace sub-shapes of a shape by other shapes.
//!
//! Source: `ModelingData/TKBRep/BRepTools/BRepTools_Substitution.cxx`
//! (OCCT 8.0.0, 172 lines):
//!
//! | OCCT | lines |
//! |---|---|
//! | `Clear` | `36-39` |
//! | `Substitute` | `43-47` |
//! | `Build` | `51-142` |
//! | `IsCopied` | `148-163` |
//! | `Copy` | `168-171` |
//!
//! Two OCCT behaviours carry the semantics:
//!
//! - a key bound to an **empty** list means "the shape disappears from its
//!   ancestors" (`IsCopied` answers `true` for it, `cxx:150-153`), and
//! - a key bound to a single shape that `IsSame`s the key means "not copied"
//!   (`cxx:157`): [`BRepToolsSubstitution::build`] binds exactly that for every
//!   shape it leaves untouched (`cxx:137-141`).
//!
//! `Build` recurses over `TopoDS_Iterator(S.Oriented(FORWARD))`, i.e. the
//! children with the parent's orientation composed in, and rebuilds a parent
//! with `EmptyCopy` only when one of its children was copied (`cxx:80-90`).
//! Degenerated results are dropped at `cxx:121-131`: a `WIRE`, `SHELL`,
//! `SOLID` or `COMPOUND` that ends up with no sub-shape is nullified.
//!
//! The port has no null `TopoDS_Shape` (`TopoShape::is_null` is always false),
//! so a nullified `NewS` is an `Option::None` in the map's value list, and the
//! "empty list" of the C++ is modelled as an **empty `Vec`**.
//!
//! `Substitute` (`cxx:45`) raises `Standard_ConstructionError` when the key is
//! already copied; the port reports that as a `false` return instead of
//! panicking, matching the convention used for other OCCT raises here.
//!
//! UNPORTED: nothing. The `Standard_NoSuchObject` raise of `Copy`
//! (`cxx:170`) becomes `None`.

use std::collections::HashMap;

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::iterator::cumulated_children;
use crate::shape::{Edge, Face, TopoShape};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{is_same, shape_map_key, ShapeMapKey};

/// `NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>,
/// TopTools_ShapeMapHasher>` (`hxx:74`): keyed by TShape + location, ignoring
/// orientation, values in OCCT list order.
#[derive(Default)]
pub struct BRepToolsSubstitution {
    map: HashMap<ShapeMapKey, Vec<TopoShape>>,
}

impl BRepToolsSubstitution {
    pub fn new() -> Self {
        Self::default()
    }

    /// `Clear` (`cxx:36-39`).
    pub fn clear(&mut self) {
        self.map.clear();
    }

    /// `Substitute(OS, NS)` (`cxx:43-47`): `NS` replaces `OS`.
    ///
    /// An empty `new_shapes` makes `OS` disappear from its ancestors. `false`
    /// is the `Standard_ConstructionError` OCCT raises when `OS` is already
    /// copied.
    pub fn substitute(&mut self, old_shape: &TopoShape, new_shapes: &[TopoShape]) -> bool {
        if self.is_copied(old_shape) {
            return false;
        }
        self.map.insert(shape_map_key(old_shape), new_shapes.to_vec());
        true
    }

    /// `Build(S)` (`cxx:51-142`): rebuild `S` from its children, recursing.
    ///
    /// A shape that is already copied is left alone (`cxx:53-56`) — that is
    /// how a caller's `Substitute` wins over the automatic rebuild.
    pub fn build(&mut self, s: &TopoShape) {
        if self.is_copied(s) {
            return;
        }

        // `cxx:59`: `TopoDS_Iterator iteS(S.Oriented(TopAbs_FORWARD))`.
        let forward = s.oriented(Orientation::Forward);
        let children = cumulated_children(&forward);

        // `cxx:66-76`
        let mut is_modified = false;
        for child in &children {
            self.build(child);
            if self.is_copied(child) {
                is_modified = true;
            }
        }

        // `cxx:78-142`
        let mut new_s: Option<TopoShape> = if is_modified {
            Some(empty_copied_forward(s))
        } else {
            Some(forward.clone())
        };
        let mut has_sub_shape = false;
        if is_modified {
            let mut rebuilt = new_s.take().expect("new shape");

            // `cxx:93-120`: add the copies of `S`'s children.
            for child in &children {
                let os = child.orientation();
                let lis = self.map.get(&shape_map_key(child)).cloned().unwrap_or_default();
                for nss in &lis {
                    // `cxx:108`: rebuild the replacement first.
                    self.build(nss);
                    let nl = self
                        .map
                        .get(&shape_map_key(nss))
                        .cloned()
                        .unwrap_or_default();
                    let new_or = Orientation::compose(os, nss.orientation());
                    for v in &nl {
                        TopoBuilder::new().add(&mut rebuilt, &v.oriented(new_or));
                        has_sub_shape = true;
                    }
                }
            }

            // `cxx:121-131`
            if !has_sub_shape {
                let t = rebuilt.shape_type();
                if t == ShapeType::Wire
                    || t == ShapeType::Shell
                    || t == ShapeType::Solid
                    || t == ShapeType::Compound
                {
                    new_s = None;
                }
            }
            if new_s.is_some() {
                new_s = Some(rebuilt);
            }
        }

        // `cxx:137-141`: `NewS` is bound FORWARD, the orientation it has in its
        // ancestors.
        let list: Vec<TopoShape> = match &new_s {
            Some(x) => vec![x.oriented(Orientation::Forward)],
            None => Vec::new(),
        };
        self.substitute(s, &list);
    }

    /// `IsCopied(S)` (`cxx:148-163`): bound, and either empty or holding a
    /// shape that is not `S` itself.
    pub fn is_copied(&self, s: &TopoShape) -> bool {
        match self.map.get(&shape_map_key(s)) {
            None => false,
            Some(list) => match list.first() {
                None => true,
                Some(first) => !is_same(s, first),
            },
        }
    }

    /// `Copy(S)` (`cxx:168-171`): the shapes substituted to `S`.
    ///
    /// `None` is the `Standard_NoSuchObject` raise OCCT performs when `S` is
    /// not copied.
    pub fn copy(&self, s: &TopoShape) -> Option<&[TopoShape]> {
        if !self.is_copied(s) {
            return None;
        }
        self.map.get(&shape_map_key(s)).map(|v| v.as_slice())
    }
}

/// `S.EmptyCopy()` preceded by `S.Oriented(FORWARD)` (`cxx:78`), carrying the
/// payload the per-class `EmptyCopy` copies:
///
/// - `BRep_TEdge::EmptyCopy` (`BRep_TEdge.cxx:105-133`): tolerance, every
///   `BRep_GCurve` / `BRep_CurveOn2Surfaces` representation (the 3D curve and
///   its range, the pcurves and theirs) and the `Degenerated` /
///   `SameParameter` / `SameRange` flags,
/// - `BRep_TFace::EmptyCopy` (`BRep_TFace.cxx:36-43`): surface, location and
///   tolerance,
/// - `BRep_TVertex::EmptyCopy` (`BRep_TVertex.cxx:33-39`): point, tolerance
///   and the point representations. Unreachable here: `Build` rebuilds a
///   vertex only when one of its children was copied, and a vertex has none.
///
/// The `B.Range(NewS, f, l)` of `cxx:87-92` follows, and it is not redundant
/// with the copy: `BRep_Builder::Range(E, f, l)` with `Only3d == false`
/// overwrites the range of **every** `BRep_GCurve` representation
/// (`BRep_Builder.cxx:1091-1117`), so a pcurve whose own range differed from
/// the 3D one loses it here.
fn empty_copied_forward(s: &TopoShape) -> TopoShape {
    let mut copy = s.empty_copied();
    copy.set_orientation(Orientation::Forward);
    let reg = GeometryRegistry::global();
    match s.shape_type() {
        ShapeType::Edge => {
            let e = Edge(s.clone());
            if let Some(g) = reg.edge_geom(&e) {
                reg.set_edge(&copy, g);
            }
            reg.set_ranges_all(
                &copy,
                BRepTool::parameter_on_edge(&e, 0),
                BRepTool::parameter_on_edge(&e, 1),
            );
        }
        ShapeType::Face => {
            let f = Face(s.clone());
            if let Some(g) = reg.face_geom(&f.0) {
                reg.set_face(&copy, g);
            }
        }
        _ => {}
    }
    copy
}
