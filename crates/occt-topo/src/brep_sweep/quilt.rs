//! `BRepTools_Quilt` (`BRepTools_Quilt.cxx`).
//!
//! `Bind` records an edge or vertex to substitute. `Add` copies faces that
//! contain a bound sub-shape. `Shells` sews the registered faces into shells.

use std::sync::Arc;

use occt_geom2d::curve::Curve2d;

use crate::abs::{Orientation, ShapeType};
use crate::builder::TopoBuilder;
use crate::iterator::cumulated_children;
use crate::shape::TopoShape;
use crate::tgeometry::GeometryRegistry;

/// `NCollection_IndexedDataMap` / `DataMap` keyed by `TopoDS_Shape::IsSame`
/// (same `TShape` and location; orientation is ignored).
struct ShapeMap {
    pairs: Vec<(TopoShape, TopoShape)>,
}

impl ShapeMap {
    fn new() -> Self {
        Self { pairs: Vec::new() }
    }

    fn contains(&self, shape: &TopoShape) -> bool {
        self.pairs.iter().any(|(k, _)| same_shape(k, shape))
    }

    fn find(&self, shape: &TopoShape) -> Option<&TopoShape> {
        self.pairs
            .iter()
            .find(|(k, _)| same_shape(k, shape))
            .map(|(_, v)| v)
    }

    fn add(&mut self, key: &TopoShape, value: &TopoShape) {
        if !self.contains(key) {
            self.pairs.push((key.clone(), value.clone()));
        }
    }

    fn bind(&mut self, key: &TopoShape, value: &TopoShape) {
        if let Some((_, slot)) = self.pairs.iter_mut().find(|(k, _)| same_shape(k, key)) {
            *slot = value.clone();
        } else {
            self.pairs.push((key.clone(), value.clone()));
        }
    }

    fn unbind(&mut self, key: &TopoShape) {
        self.pairs.retain(|(k, _)| !same_shape(k, key));
    }

    fn extent(&self) -> usize {
        self.pairs.len()
    }

    fn value(&self, index: usize) -> &TopoShape {
        &self.pairs[index - 1].1
    }
}

/// Shapes recorded by `IsSame`, insertion order.
struct ShapeSet {
    items: Vec<TopoShape>,
}

impl ShapeSet {
    fn new() -> Self {
        Self { items: Vec::new() }
    }

    fn add(&mut self, shape: &TopoShape) {
        if !self.contains(shape) {
            self.items.push(shape.clone());
        }
    }

    fn contains(&self, shape: &TopoShape) -> bool {
        self.items.iter().any(|s| same_shape(s, shape))
    }
}

pub(super) struct Quilt {
    bounds: ShapeMap,
    has_copy: bool,
}

impl Quilt {
    pub(super) fn new() -> Self {
        Self {
            bounds: ShapeMap::new(),
            has_copy: false,
        }
    }

    /// `BRepTools_Quilt::Add` (`cxx:114-283`).
    pub(super) fn add(&mut self, shape: &TopoShape) {
        if self.bounds.contains(shape) {
            return;
        }
        let builder = TopoBuilder::new();
        for wire in explore(shape, ShapeType::Wire, Some(ShapeType::Face)) {
            self.bounds.add(&wire, &wire);
        }
        for edge in explore(shape, ShapeType::Edge, Some(ShapeType::Wire)) {
            self.bounds.add(&edge, &edge);
        }
        for vertex in explore(shape, ShapeType::Vertex, Some(ShapeType::Edge)) {
            self.bounds.add(&vertex, &vertex);
        }
        for face in explore(shape, ShapeType::Face, None) {
            let mut copy_face = false;
            if self.has_copy {
                for fed in explore(&face, ShapeType::Edge, None) {
                    if self.bounds.contains(&fed) {
                        copy_face = true;
                    } else if need_copied(&fed, &self.bounds) {
                        copy_face = true;
                        copy_edge(&fed, &mut self.bounds);
                    }
                }
            }
            let mut nf = face.clone();
            if copy_face {
                nf = empty_copy_face(&face);
                nf.set_orientation(Orientation::Forward);
                for wire in children_cum_loc(&face) {
                    if wire.shape_type() != ShapeType::Wire {
                        continue;
                    }
                    let mut nw = TopoShape::new(ShapeType::Wire);
                    for edge in children_cum_loc(&wire) {
                        if edge.shape_type() != ShapeType::Edge {
                            continue;
                        }
                        let mut oe = edge.orientation();
                        if self.bounds.contains(&edge) {
                            // `cxx:237-267`.
                            let ne = self.bounds.find(&edge).unwrap().clone();
                            if let Some((curve, mut u_first, mut u_last)) = pcurve_of(&edge, &face)
                            {
                                if ne.orientation() == Orientation::Forward {
                                    write_pcurve(&ne, &face, curve, u_first, u_last);
                                } else {
                                    oe = oe.reversed();
                                    let u0 = curve.reversed_parameter(u_last);
                                    let u1 = curve.reversed_parameter(u_first);
                                    let reversed = Arc::from(curve.reversed());
                                    u_first = u0;
                                    u_last = u1;
                                    write_pcurve(&ne, &face, reversed, u_first, u_last);
                                }
                            }
                            builder.add(&mut nw, &ne.oriented(oe));
                        } else {
                            builder.add(&mut nw, &edge);
                        }
                    }
                    nw.set_orientation(wire.orientation());
                    builder.add(&mut nf, &nw);
                }
                nf.set_orientation(face.orientation());
            }
            self.bounds.add(&face, &nf);
        }
    }

    /// `BRepTools_Quilt::Bind` for a vertex (`cxx:287-293`).
    #[allow(dead_code)]
    pub(super) fn bind_vertex(&mut self, old: &TopoShape, new: &TopoShape) {
        if !self.bounds.contains(old) {
            self.bounds.add(old, new);
        }
    }

    /// `BRepTools_Quilt::Bind` for an edge (`cxx:297-345`).
    #[allow(dead_code)]
    pub(super) fn bind_edge(&mut self, old: &TopoShape, new: &TopoShape) {
        if self.bounds.contains(old) {
            return;
        }
        let mut new_edge = new.clone();
        if self.is_copied(new) {
            new_edge = self.bounds.find(new).unwrap().clone();
            new_edge.set_orientation(new.orientation());
        }
        if old.orientation() != new_edge.orientation() {
            self.bounds.add(
                &old.oriented(Orientation::Forward),
                &new_edge.oriented(Orientation::Reversed),
            );
        } else {
            self.bounds.add(
                &old.oriented(Orientation::Forward),
                &new_edge.oriented(Orientation::Forward),
            );
        }
        for old_v in cumulated_children(old) {
            if self.bounds.contains(&old_v) {
                continue;
            }
            let orien = old_v.orientation();
            for new_v in cumulated_children(&new_edge) {
                if new_v.orientation() != orien {
                    continue;
                }
                let placed = if self.is_copied(&new_v) {
                    self.bounds.find(&new_v).unwrap().clone()
                } else {
                    new_v
                };
                self.bounds.add(&old_v, &placed);
                break;
            }
        }
        self.has_copy = true;
    }

    fn is_copied(&self, shape: &TopoShape) -> bool {
        // `cxx:349-358`.
        match self.bounds.find(shape) {
            Some(bound) => !same_shape(shape, bound),
            None => false,
        }
    }

    /// `BRepTools_Quilt::Shells` (`cxx:371-589`).
    pub(super) fn shells(&self) -> TopoShape {
        let mut free_edges = ShapeMap::new();
        let mut face_shell = ShapeMap::new();
        let builder = TopoBuilder::new();
        let mut result = TopoShape::new(ShapeType::Compound);
        let mut other = ShapeSet::new();
        let mut edges_of_faces = ShapeSet::new();

        for index in 1..=self.bounds.extent() {
            let shape = self.bounds.value(index).clone();
            if shape.shape_type() != ShapeType::Face {
                other.add(&shape);
                continue;
            }
            for edge in explore(&shape, ShapeType::Edge, None) {
                edges_of_faces.add(&edge);
            }
            let mut shell: Option<TopoShape> = None;
            let mut new_o = shape.orientation();
            for edge in explore(&shape, ShapeType::Edge, None) {
                if let Some(bound) = free_edges.find(&edge) {
                    let bound = bound.clone();
                    new_o = if bound.orientation() == edge.orientation() {
                        shape.orientation().reversed()
                    } else {
                        shape.orientation()
                    };
                    face_shell.bind(&shape, &bound.oriented(new_o));
                    shell = Some(bound);
                    break;
                }
            }
            let mut shell = if let Some(shell) = shell {
                shell
            } else {
                let shell = TopoShape::new(ShapeType::Shell);
                shell.set_closed(true);
                builder.add(&mut result, &shell);
                face_shell.bind(&shape, &shell.oriented(shape.orientation()));
                shell
            };
            shell.set_free(true);
            let face_orien = face_shell
                .find(&shape)
                .map(|s| s.orientation())
                .unwrap_or(shape.orientation());
            let mut host = shell.oriented(Orientation::Forward);
            builder.add(&mut host, &shape.oriented(face_orien));

            let forward_face = shape.oriented(Orientation::Forward);
            for edge in explore(&forward_face, ShapeType::Edge, None) {
                if let Some(old_shell) = free_edges.find(&edge).cloned() {
                    if !same_shape(&old_shell, &shell) {
                        let mut an_orien = edge.orientation();
                        if face_orien == Orientation::Reversed {
                            an_orien = an_orien.reversed();
                        }
                        let rev = an_orien == old_shell.orientation();
                        for old_face in cumulated_children(&old_shell) {
                            if old_face.shape_type() != ShapeType::Face {
                                continue;
                            }
                            let stored = face_shell
                                .find(&old_face)
                                .map(|s| s.orientation())
                                .unwrap_or(old_face.orientation());
                            let new_of = if rev { stored.reversed() } else { stored };
                            face_shell.bind(&old_face, &shell.oriented(new_of));
                            let mut host_fo = shell.oriented(Orientation::Forward);
                            builder.add(&mut host_fo, &old_face.oriented(new_of));
                        }
                        for ae in explore(&shell, ShapeType::Edge, None) {
                            let Some(as_shell) = free_edges.find(&ae).cloned() else {
                                continue;
                            };
                            if same_shape(&as_shell, &old_shell) {
                                let rebound = if rev {
                                    as_shell.orientation().reversed()
                                } else {
                                    as_shell.orientation()
                                };
                                free_edges.bind(&ae, &shell.oriented(rebound));
                            }
                        }
                        let dropped = old_shell.oriented(Orientation::Forward);
                        builder.remove(&mut result, &dropped);
                    }
                    let mut an_orien = edge.orientation();
                    if face_orien == Orientation::Reversed {
                        an_orien = an_orien.reversed();
                    }
                    if free_edges
                        .find(&edge)
                        .map(|s| s.orientation() == an_orien)
                        .unwrap_or(false)
                    {
                        shell
                            .tshape
                            .write()
                            .expect("tshape")
                            .flags
                            .oriented = false;
                    }
                    free_edges.unbind(&edge);
                } else {
                    let mut new_o = edge.orientation();
                    if face_orien == Orientation::Reversed {
                        new_o = new_o.reversed();
                    }
                    free_edges.bind(&edge, &shell.oriented(new_o));
                }
            }
            shell.set_free(false);
        }

        for index in 0..free_edges.extent() {
            if let Some(shell) = free_edges.pairs.get(index).map(|(_, v)| v.clone()) {
                shell.set_closed(false);
            }
        }
        for shape in &other.items {
            if !edges_of_faces.contains(shape) && self.bounds.contains(shape) {
                if let Some(bound) = self.bounds.find(shape) {
                    builder.add(&mut result, bound);
                }
            }
        }
        result
    }
}

/// `ShapeMap` / `ShapeSet` key comparison: `TopoDS_Shape::IsSame`.
fn same_shape(a: &TopoShape, b: &TopoShape) -> bool {
    crate::topo_tools::is_same(a, b)
}

fn need_copied(shape: &TopoShape, bounds: &ShapeMap) -> bool {
    // `NeedCopied` (`cxx:46-63`).
    cumulated_children(shape).iter().any(|c| bounds.contains(c))
}

fn children_cum_loc(shape: &TopoShape) -> Vec<TopoShape> {
    // `TopoDS_Iterator(S, cumOri=false)`.
    let stored = shape
        .tshape
        .read()
        .expect("tshape")
        .children
        .clone();
    let parent = shape.location().clone();
    if parent.is_identity() {
        return stored;
    }
    stored
        .into_iter()
        .map(|mut child| {
            child.move_location(&parent);
            child
        })
        .collect()
}

fn explore(root: &TopoShape, target: ShapeType, avoid: Option<ShapeType>) -> Vec<TopoShape> {
    // `TopExp_Explorer(S, target, avoid)`.
    let mut out = Vec::new();
    let mut stack = Vec::new();
    for child in cumulated_children(root).into_iter().rev() {
        stack.push(child);
    }
    while let Some(shape) = stack.pop() {
        let ty = shape.shape_type();
        if ty == target {
            out.push(shape.clone());
        }
        if avoid == Some(ty) {
            continue;
        }
        for child in cumulated_children(&shape).into_iter().rev() {
            stack.push(child);
        }
    }
    out
}

fn copy_edge(edge: &TopoShape, bounds: &mut ShapeMap) {
    // `CopyShape` (`cxx:65-93`).
    let mut copy = empty_copy_edge(edge);
    copy.set_orientation(Orientation::Forward);
    let builder = TopoBuilder::new();
    for vertex in children_cum_loc(edge) {
        if bounds.contains(&vertex) {
            let bound = bounds.find(&vertex).unwrap().oriented(vertex.orientation());
            builder.add(&mut copy, &bound);
        } else {
            builder.add(&mut copy, &vertex);
        }
    }
    let (first, last) = GeometryRegistry::global().edge_parameters(edge);
    set_edge_range(&copy, first, last);
    bounds.add(edge, &copy.oriented(Orientation::Forward));
}

fn empty_copy_edge(src: &TopoShape) -> TopoShape {
    // `BRep_TEdge::EmptyCopy` (`BRep_TEdge.cxx:105-128`) plus the location.
    let mut edge = TopoShape::new(ShapeType::Edge);
    edge.set_location(src.location());
    let (core, pcurves) = {
        let ts = src.tshape.read().expect("tshape");
        (ts.edge_core().cloned(), ts.edge_pcurves().cloned())
    };
    let mut ts = edge.tshape.write().expect("tshape");
    if let Some(core) = core {
        *ts.edge_core_mut() = core;
    }
    if let Some(pc) = pcurves {
        ts.edge_pcurves = Some(pc);
    }
    drop(ts);
    edge
}

fn empty_copy_face(src: &TopoShape) -> TopoShape {
    let mut face = TopoShape::new(ShapeType::Face);
    face.set_location(src.location());
    let core = src.tshape.read().expect("tshape").face_core().cloned();
    if let Some(core) = core {
        *face.tshape.write().expect("tshape").face_core_mut() = core;
    }
    face
}

fn set_edge_range(edge: &TopoShape, first: f64, last: f64) {
    let mut ts = edge.tshape.write().expect("tshape");
    let core = ts.edge_core_mut();
    core.first = first;
    core.last = last;
}

fn pcurve_of(edge: &TopoShape, face: &TopoShape) -> Option<(Arc<dyn Curve2d>, f64, f64)> {
    let reg = GeometryRegistry::global();
    let key = GeometryRegistry::shape_key(face);
    let curve = reg.edge_pcurve(edge, key)?;
    let (first, last) = reg
        .pcurve_range(edge, key)
        .unwrap_or_else(|| reg.edge_parameters(edge));
    Some((curve, first, last))
}

fn write_pcurve(edge: &TopoShape, face: &TopoShape, curve: Arc<dyn Curve2d>, first: f64, last: f64) {
    let reg = GeometryRegistry::global();
    let key = GeometryRegistry::shape_key(face);
    let tol = reg.edge_tolerance(edge);
    reg.set_edge_pcurve(edge, key, curve);
    reg.set_pcurve_range(edge, key, first, last);
    let mut ts = edge.tshape.write().expect("tshape");
    ts.edge_core_mut().tolerance = tol;
}
