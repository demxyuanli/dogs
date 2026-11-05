//! `BRepAlgo_Loop` supporting statics: the vertex -> edges map and the five
//! free functions of `BRepAlgo_Loop.cxx:230-578`.
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `RemovePendingEdges` | `cxx:230-296` | [`remove_pending_edges`] |
//! | `SamePnt2d` | `cxx:300-329` | [`same_pnt2d`] |
//! | `SelectEdge` | `cxx:337-444` | [`select_edge`] |
//! | `PurgeNewEdges` | `cxx:446-472` | [`purge_new_edges`] |
//! | `StoreInMVE` | `cxx:474-578` | [`store_in_mve`] |

use std::collections::HashMap;

use crate::abs::Orientation;
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::iterator::cumulated_children;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topo_tools::{edge_vertices, is_same, shape_map_key, ShapeMapKey};

/// `NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>>`
/// (`MVE`): insertion-ordered, 1-based indexing, `IsSame` keys.
#[derive(Clone, Default)]
pub struct VertexEdgeMap {
    entries: Vec<(TopoShape, Vec<TopoShape>)>,
}

impl VertexEdgeMap {
    pub fn new() -> Self {
        Self::default()
    }

    /// `Extent()`.
    pub fn extent(&self) -> usize {
        self.entries.len()
    }

    /// `Contains(S)`.
    pub fn contains(&self, s: &TopoShape) -> bool {
        let key = shape_map_key(s);
        self.entries.iter().any(|(k, _)| shape_map_key(k) == key)
    }

    /// `Add(S, L)`: `NCollection_IndexedDataMap::Add` keeps the existing entry
    /// when the key is already present.
    pub fn add(&mut self, s: &TopoShape, list: Vec<TopoShape>) {
        if self.contains(s) {
            return;
        }
        self.entries.push((s.clone(), list));
    }

    /// `FindKey(i)` for the 1-based `i`.
    pub fn find_key(&self, i: usize) -> Option<&TopoShape> {
        self.entries.get(i.checked_sub(1)?).map(|(k, _)| k)
    }

    /// `operator()(i)` for the 1-based `i`.
    pub fn find_from_index(&self, i: usize) -> Option<&Vec<TopoShape>> {
        self.entries.get(i.checked_sub(1)?).map(|(_, v)| v)
    }

    /// `FindFromKey(S)`.
    pub fn find_from_key(&self, s: &TopoShape) -> Option<&Vec<TopoShape>> {
        let key = shape_map_key(s);
        self.entries
            .iter()
            .find(|(k, _)| shape_map_key(k) == key)
            .map(|(_, v)| v)
    }

    /// `ChangeFromKey(S)`.
    pub fn change_from_key(&mut self, s: &TopoShape) -> Option<&mut Vec<TopoShape>> {
        let key = shape_map_key(s);
        self.entries
            .iter_mut()
            .find(|(k, _)| shape_map_key(k) == key)
            .map(|(_, v)| v)
    }

    /// `RemoveKey(S)`.
    pub fn remove_key(&mut self, s: &TopoShape) {
        let key = shape_map_key(s);
        self.entries.retain(|(k, _)| shape_map_key(k) != key);
    }
}

/// `RemovePendingEdges(MVE)` (`cxx:230-296`): drop vertices with no edge and
/// open edges that hang on a single vertex, repeating until nothing changes.
pub fn remove_pending_edges(mve: &mut VertexEdgeMap) {
    let mut ya_supress = true;
    while ya_supress {
        ya_supress = false;
        let mut v_to_remove: Vec<TopoShape> = Vec::new();
        let mut e_to_remove: Vec<TopoShape> = Vec::new();

        for i_v in 1..=mve.extent() {
            let Some(vertex) = mve.find_key(i_v).cloned() else {
                continue;
            };
            let edges = mve.find_from_index(i_v).cloned().unwrap_or_default();
            if edges.is_empty() {
                v_to_remove.push(vertex.clone());
            }
            if edges.len() == 1 {
                let (v1, v2) = edge_vertices(&Edge(edges[0].clone()));
                let same = match (&v1, &v2) {
                    (Some(a), Some(b)) => is_same(&a.0, &b.0),
                    (None, None) => true,
                    _ => false,
                };
                if !same {
                    v_to_remove.push(vertex.clone());
                    e_to_remove.push(edges[0].clone());
                }
            }
        }

        if !v_to_remove.is_empty() {
            ya_supress = true;
            for v in &v_to_remove {
                mve.remove_key(v);
            }
            if !e_to_remove.is_empty() {
                let keys: Vec<ShapeMapKey> = e_to_remove.iter().map(shape_map_key).collect();
                for (_, list) in mve.entries.iter_mut() {
                    list.retain(|e| !keys.contains(&shape_map_key(e)));
                }
            }
        }
    }
}

/// `SamePnt2d(V, E1, E2, F)` (`cxx:300-329`): whether the two edges end at the
/// same 2D point on the face near `V`.
pub fn same_pnt2d(v: &Vertex, e1: &TopoShape, e2: &TopoShape, f: &Face) -> bool {
    let ff = Face(f.0.clone().oriented(Orientation::Forward));
    let Some((c1, f1, l1)) = curve_on_surface_range(&Edge(e1.clone()), &ff) else {
        panic!("BRep_Tool::CurveOnSurface(E1, FF, f1, l1) raised");
    };
    let Some((c2, f2, l2)) = curve_on_surface_range(&Edge(e2.clone()), &ff) else {
        panic!("BRep_Tool::CurveOnSurface(E2, FF, f2, l2) raised");
    };
    let p1 = if e1.orientation() == Orientation::Forward {
        c1.d0(f1)
    } else {
        c1.d0(l1)
    };
    let p2 = if e2.orientation() == Orientation::Forward {
        c2.d0(l2)
    } else {
        c2.d0(f2)
    };
    let tol = 100.0 * BRepTool::vertex_tolerance(v);
    p1.distance(&p2) < tol
}

/// `SelectEdge(F, CE, CV, NE, LE)` (`cxx:337-444`): pick the edge of `LE`
/// continuing `CE` at `CV`, removing it (and `CE` itself) from `LE`.
pub fn select_edge(
    f: &Face,
    ce: &TopoShape,
    cv: &TopoShape,
    ne: &mut Option<TopoShape>,
    le: &mut Vec<TopoShape>,
) -> bool {
    *ne = None;

    if let Some(i) = le.iter().position(|x| crate::topo_tools::is_equal(x, ce)) {
        le.remove(i);
    }

    if le.len() > 1 {
        let ff = Face(f.0.clone().oriented(Orientation::Forward));
        let Some((c, fpar, lpar)) = curve_on_surface_range(&Edge(ce.clone()), &ff) else {
            panic!("BRep_Tool::CurveOnSurface(CE, FForward, f, l) raised");
        };
        let mut k = 1usize;
        let mut kmin = 0usize;
        let mut distmin = 100.0 * BRepTool::vertex_tolerance(&Vertex(cv.clone()));
        let u = if ce.orientation() == Orientation::Forward {
            lpar
        } else {
            fpar
        };
        let pv = c.d0(u);

        for e in le.iter() {
            if !is_same(e, ce) {
                let Some((c2, f2, l2)) = curve_on_surface_range(&Edge(e.clone()), &ff) else {
                    panic!("BRep_Tool::CurveOnSurface(E, FForward, f, l) raised");
                };
                let u2 = if e.orientation() == Orientation::Forward {
                    f2
                } else {
                    l2
                };
                let p2 = c2.d0(u2);
                let dist = pv.distance(&p2);
                if dist <= distmin {
                    kmin = k;
                    distmin = dist;
                }
            }
            k += 1;
        }
        if kmin == 0 {
            return false;
        }
        *ne = Some(le[kmin - 1].clone());
        le.remove(kmin - 1);
    } else if le.len() == 1 {
        *ne = Some(le[0].clone());
        le.remove(0);
    } else {
        return false;
    }
    true
}

/// `PurgeNewEdges(NewEdges, UsedEdges)` (`cxx:446-472`): drop the cut edges
/// that no wire consumed.
pub fn purge_new_edges(
    new_edges: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
    used_edges: &std::collections::HashSet<ShapeMapKey>,
) {
    for list in new_edges.values_mut() {
        list.retain(|ne| used_edges.contains(&shape_map_key(ne)));
    }
}

/// `StoreInMVE(F, E, MVE, YaCouture, VerticesForSubstitute, theTolConf)`
/// (`cxx:474-578`): substitute coincident vertices of `E`, then register `E`
/// (and its reverse when it lies on a closing curve) under both end vertices.
pub fn store_in_mve(
    f: &Face,
    e: &mut TopoShape,
    mve: &mut VertexEdgeMap,
    ya_couture: &mut bool,
    vertices_for_substitute: &mut HashMap<ShapeMapKey, TopoShape>,
    the_tol_conf: f64,
) {
    let builder = TopoBuilder::new();

    // `cxx:489-545`: the map itself is not modified inside this walk, so the
    // key snapshot is the `for (iV = 1; iV <= MVE.Extent(); iV++)` sequence.
    let keyed: Vec<TopoShape> = (1..=mve.extent())
        .filter_map(|i| mve.find_key(i).cloned())
        .collect();
    for v in keyed {
        let p = BRepTool::vertex_point(&Vertex(v.clone()));
        let vlist: Vec<TopoShape> = cumulated_children(e);

        for v1 in vlist {
            let p1 = BRepTool::vertex_point(&Vertex(v1.clone()));
            if p.is_equal_with_tolerance(&p1, the_tol_conf) && !is_same(&v, &v1) {
                let mut v = v.clone();
                v.set_orientation(v1.orientation());
                let v1_key = shape_map_key(&v1);
                if let Some(old_new_v) = vertices_for_substitute.get(&v1_key).cloned() {
                    if !is_same(&old_new_v, &v) {
                        vertices_for_substitute.insert(shape_map_key(&old_new_v), v.clone());
                        vertices_for_substitute.insert(v1_key, v.clone());
                    }
                } else {
                    let v_key = shape_map_key(&v);
                    match vertices_for_substitute.get(&v_key).cloned() {
                        Some(new_new_v) => {
                            if !is_same(&new_new_v, &v1) {
                                vertices_for_substitute.insert(v1_key, new_new_v);
                            }
                        }
                        None => {
                            vertices_for_substitute.insert(v1_key, v.clone());
                            let keys: Vec<ShapeMapKey> =
                                vertices_for_substitute.keys().cloned().collect();
                            for k in keys {
                                let same = vertices_for_substitute
                                    .get(&k)
                                    .map(|x| is_same(x, &v1))
                                    .unwrap_or(false);
                                if same {
                                    vertices_for_substitute.insert(k, v.clone());
                                }
                            }
                        }
                    }
                }
                e.set_free(true);
                builder.remove(e, &v1);
                builder.add(e, &v);
            }
        }
    }

    let (v1, v2) = edge_vertices(&Edge(e.clone()));
    if v1.is_none() && v2.is_none() {
        *ya_couture = false;
        return;
    }
    if let Some(v1) = &v1 {
        if !mve.contains(&v1.0) {
            mve.add(&v1.0, Vec::new());
        }
        if let Some(list) = mve.change_from_key(&v1.0) {
            list.push(e.clone());
        }
    }
    let same_ends = match (&v1, &v2) {
        (Some(a), Some(b)) => is_same(&a.0, &b.0),
        (None, None) => true,
        _ => false,
    };
    if !same_ends {
        if let Some(v2) = &v2 {
            if !mve.contains(&v2.0) {
                mve.add(&v2.0, Vec::new());
            }
            if let Some(list) = mve.change_from_key(&v2.0) {
                list.push(e.clone());
            }
        }
    }

    let closed = BRepTool::is_closed_edge_face(&Edge(e.clone()), f);
    if closed {
        if let Some(v2) = &v2 {
            let mut rev = e.clone();
            rev.reverse();
            if let Some(list) = mve.change_from_key(&v2.0) {
                list.push(rev.clone());
            }
            if !same_ends {
                if let Some(v1) = &v1 {
                    if let Some(list) = mve.change_from_key(&v1.0) {
                        list.push(rev);
                    }
                }
            }
        }
    }
    if closed {
        *ya_couture = true;
    }
}
