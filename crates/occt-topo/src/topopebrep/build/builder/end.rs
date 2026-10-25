//! `TopOpeBRepBuild_Builder::End` — the post-merge fixups.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_End.cxx`
//! (OCCT 8.0.0, 512 lines). `MergeKPart` calls it last
//! (`TopOpeBRepBuild_KPart.cxx:208`).
//!
//! Three blocks:
//! 1. continuity recoding (`cxx:47-193`): a `(edge, (f1, f2))` continuity other
//!    than `C0` is lost when a tangent face's support surface is changed on an
//!    edge that was not split; the block re-stores it on the split faces.
//! 2. vertex tolerance update (`cxx:195-365`): every vertex lying on a section
//!    edge is widened to enclose the 3D and p-curve images of its incident
//!    edges.
//! 3. `TopOpeBRepTool_FuseEdges` (`cxx:366-410`) with
//!    `UpdateSplitAndMerged` (`cxx:413-520`), in
//!    [`end_fuse_edges`](super::end_fuse_edges).
//!
//! Blocks 1 and 2 run unconditionally, exactly as OCCT; block 3 is gated only
//! by `makeFE`, which is constant `true` outside `OCCT_DEBUG`.

use std::collections::{HashMap, HashSet};

use occt_core::bnd::BndBox;
use occt_core::gp::GpPnt;

use crate::abs::ShapeType;
use crate::boptools_2d;
use crate::brep_lib_continuity::{set_stored_continuity, stored_continuity};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{
    explorer_shapes, explorer_shapes_avoiding, is_same, shape_map_key, ShapeMapKey,
};

use super::section_edge::fc2d_curve_on_surface;
use super::Builder;

/// `NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>,
/// TopTools_ShapeMapHasher>` restricted to what `End` uses: insertion-ordered
/// keys with an ancestor list per key.
///
/// `TopExp::MapShapesAndAncestors` / `TopExp::MapShapesAndUniqueAncestors`
/// append into the map they are given (`TopExp.cxx:80-120`, `:124-177`), so the
/// map must accumulate across calls — `idmovloes` is filled once per section
/// edge (`End.cxx:245-248`).
#[derive(Default)]
struct AncestorMap {
    keys: Vec<TopoShape>,
    vals: Vec<Vec<TopoShape>>,
    index: HashMap<ShapeMapKey, usize>,
}

impl AncestorMap {
    fn new() -> Self {
        Self::default()
    }

    /// The entry of `s`, created empty when absent.
    fn entry(&mut self, s: &TopoShape) -> usize {
        let k = shape_map_key(s);
        if let Some(&i) = self.index.get(&k) {
            return i;
        }
        let i = self.keys.len();
        self.keys.push(s.clone());
        self.vals.push(Vec::new());
        self.index.insert(k, i);
        i
    }

    /// `TopExp::MapShapesAndAncestors(shape, ts, ta)` appended into this map.
    fn add_ancestors(&mut self, shape: &TopoShape, ts: ShapeType, ta: ShapeType) {
        for anc in explorer_shapes(shape, ta) {
            for sub in explorer_shapes(&anc, ts) {
                let i = self.entry(&sub);
                self.vals[i].push(anc.clone());
            }
        }
        for sub in explorer_shapes_avoiding(shape, ts, ta) {
            self.entry(&sub);
        }
    }

    /// `TopExp::MapShapesAndUniqueAncestors(shape, ts, ta)` appended into this
    /// map, `IsSame` comparison (`TopExp.cxx:156`, `useOrientation = false`).
    fn add_unique_ancestors(&mut self, shape: &TopoShape, ts: ShapeType, ta: ShapeType) {
        for anc in explorer_shapes(shape, ta) {
            for sub in explorer_shapes(&anc, ts) {
                let i = self.entry(&sub);
                let already = self.vals[i]
                    .iter()
                    .any(|x| x.same_tshape(&anc) && x.location().is_equal(&anc.location()));
                if !already {
                    self.vals[i].push(anc.clone());
                }
            }
        }
        for sub in explorer_shapes_avoiding(shape, ts, ta) {
            self.entry(&sub);
        }
    }

    /// `Extent()` (`End.cxx:99`).
    fn extent(&self) -> usize {
        self.keys.len()
    }

    /// `FindKey(i)` (`End.cxx:102`), `i` being one-based.
    fn find_key(&self, i: usize) -> Option<&TopoShape> {
        self.keys.get(i - 1)
    }

    /// `FindFromIndex(i)` (`End.cxx:106`), `i` being one-based.
    fn find_from_index(&self, i: usize) -> &[TopoShape] {
        self.vals.get(i - 1).map(|v| v.as_slice()).unwrap_or(&[])
    }

    /// `FindFromKey(S)`. OCCT throws `Standard_OutOfRange` when `S` is absent;
    /// the callers here probe shapes the same map already collected, so the
    /// empty slice is unreachable.
    fn find_from_key(&self, s: &TopoShape) -> &[TopoShape] {
        self.index
            .get(&shape_map_key(s))
            .map(|&i| self.vals[i].as_slice())
            .unwrap_or(&[])
    }

    /// `Contains(S)` (`End.cxx:253`).
    fn contains(&self, s: &TopoShape) -> bool {
        self.index.contains_key(&shape_map_key(s))
    }
}

impl Builder {
    /// `End()` (`End.cxx:45-410`).
    pub fn end(&mut self) {
        self.end_recode_continuity();
        self.end_update_vertex_tolerances();
        // `End.cxx:366-410`: fuse the operation's edges in the merged results.
        self.end_fuse_edges();
    }

    /// `End.cxx:47-193`: recode the continuity of the shared tangent edges on
    /// the split faces.
    fn end_recode_continuity(&mut self) {
        for i in 1..=2 {
            let (s, sta) = if i == 1 {
                (self.my_shape1.clone(), self.my_state1)
            } else {
                (self.my_shape2.clone(), self.my_state2)
            };
            let Some(s) = s else {
                continue;
            };

            // `End.cxx:94-98`: OCCT rebuilds `M` inside the shell loop; its
            // content does not depend on the current shell, so it is built once
            // here and the iteration order stays `MapShapesAndAncestors`'s.
            let mut m = AncestorMap::new();
            m.add_ancestors(&s, ShapeType::Edge, ShapeType::Face);

            for sh in explorer_shapes(&s, ShapeType::Shell) {
                // `End.cxx:72-92`: the shell must be known to the DS and carry
                // at least one same-domain face.
                let Some(hds) = self.my_data_structure.as_ref() else {
                    continue;
                };
                if !hds.has_shape(&sh, true) {
                    continue;
                }
                let mut fhassam = false;
                for f in explorer_shapes(&sh, ShapeType::Face) {
                    if hds.has_same_domain(&f, true) {
                        fhassam = true;
                        break;
                    }
                }
                if !fhassam {
                    continue;
                }

                // `End.cxx:99-190`
                let ne = m.extent();
                for ie in 1..=ne {
                    let Some(e) = m.find_key(ie).cloned() else {
                        continue;
                    };
                    // `End.cxx:103-106`
                    if self.is_split(&e, sta) {
                        continue;
                    }
                    let lf = m.find_from_index(ie).to_vec();
                    if lf.len() < 2 {
                        continue;
                    }

                    // `End.cxx:113-118`: OCCT walks the first pair of connex
                    // faces only ("NYI : > 2 faces connexes par E").
                    let f1 = lf[0].clone();
                    let f2 = lf[1].clone();
                    let (Some(edge), Some(face1), Some(face2)) = (
                        Edge::wrap(e.clone()),
                        Face::wrap(f1.clone()),
                        Face::wrap(f2.clone()),
                    ) else {
                        continue;
                    };
                    let c = stored_continuity(&edge, &face1, &face2);
                    // `End.cxx:119-122`: `GeomAbs_C0` is 0.
                    if c == 0 {
                        continue;
                    }

                    // `End.cxx:124-131`
                    let (f1hassam, f2hassam) =
                        (hds.has_same_domain(&f1, true), hds.has_same_domain(&f2, true));
                    if !f1hassam && !f2hassam {
                        continue;
                    }

                    // `End.cxx:133-138`
                    let f1issplit = self.is_split(&f1, sta) && !self.splits(&f1, sta).is_empty();
                    let f2issplit = self.is_split(&f2, sta) && !self.splits(&f2, sta).is_empty();
                    if !f1issplit && !f2issplit {
                        continue;
                    }

                    // `End.cxx:140-188`
                    let mut ff1 = f1.clone();
                    let mut ff2 = f2.clone();
                    for ii in 1..=2 {
                        if (ii == 1 && !f1issplit) || (ii == 2 && !f2issplit) {
                            continue;
                        }
                        let f = if ii == 1 { &f1 } else { &f2 };
                        let mut found = false;
                        for sf in self.splits(f, sta).to_vec() {
                            if sf.shape_type() != ShapeType::Face {
                                continue;
                            }
                            for child in explorer_shapes(&sf, ShapeType::Edge) {
                                if is_same(&child, &e) {
                                    if ii == 1 {
                                        ff1 = sf.clone();
                                    } else {
                                        ff2 = sf.clone();
                                    }
                                    found = true;
                                    break;
                                }
                            }
                            if found {
                                break;
                            }
                        }
                    }

                    // `End.cxx:188-189`: `B.Continuity(E, FF1, FF2, C)`.
                    if let (Some(a), Some(b)) = (Face::wrap(ff1), Face::wrap(ff2)) {
                        set_stored_continuity(&edge, &a, &b, c);
                    }
                }
            }
        }
    }

    /// `End.cxx:195-365`: widen the tolerance of every vertex lying on a
    /// section edge to enclose the 3D and p-curve images of its incident edges.
    fn end_update_vertex_tolerances(&mut self) {
        // `End.cxx:197-212`
        let mut map_of_new_edges: HashSet<ShapeMapKey> = HashSet::new();
        let nb_curves = self
            .my_data_structure
            .as_ref()
            .map(|hds| hds.nb_curves())
            .unwrap_or(0);
        for i in 1..=nb_curves {
            for e in self.new_edges(i) {
                map_of_new_edges.insert(shape_map_key(e));
            }
        }
        let mut map_of_new_vertices: HashSet<ShapeMapKey> = HashSet::new();
        let nb_points = self
            .my_data_structure
            .as_ref()
            .map(|hds| hds.nb_points())
            .unwrap_or(0);
        for i in 1..=nb_points {
            if let Some(v) = self.new_vertex(i) {
                map_of_new_vertices.insert(shape_map_key(v));
            }
        }

        // `End.cxx:214-228`: a compound holding every merged shape of operand 1.
        let lmergesha1 = match self.my_shape1.clone() {
            Some(s) => self.merged(&s, self.my_state1).to_vec(),
            None => Vec::new(),
        };
        let topo = TopoBuilder::new();
        let mut r = topo.make_shape(ShapeType::Compound);
        for s in &lmergesha1 {
            topo.add(&mut r, s);
        }

        // `End.cxx:229`
        let loes = self.section().to_vec();

        // `End.cxx:230-238`
        let mut idmoelof = AncestorMap::new();
        idmoelof.add_ancestors(&r, ShapeType::Edge, ShapeType::Face);
        let mut idmovloe = AncestorMap::new();
        idmovloe.add_unique_ancestors(&r, ShapeType::Vertex, ShapeType::Edge);

        // `End.cxx:241-248`
        let mut idmovloes = AncestorMap::new();
        for s in &loes {
            idmovloes.add_ancestors(s, ShapeType::Vertex, ShapeType::Edge);
        }

        // `End.cxx:248-365`
        let nv = idmovloe.extent();
        for iv in 1..=nv {
            let Some(v) = idmovloe.find_key(iv).cloned() else {
                continue;
            };
            // `End.cxx:253-257`
            if !idmovloes.contains(&v) {
                continue;
            }
            let loe = idmovloe.find_from_index(iv).to_vec();

            // `End.cxx:259-271`
            let mut np1 = 0usize;
            for e in &loe {
                np1 += idmoelof.find_from_key(e).len() + 1;
            }

            // `End.cxx:273-337`
            let mut tp: Vec<GpPnt> = Vec::with_capacity(np1);
            for e in &loe {
                let (Some(edge), Some(vertex)) =
                    (Edge::wrap(e.clone()), Vertex::wrap(v.clone()))
                else {
                    continue;
                };
                // `End.cxx:278`: `BRep_Tool::Parameter(V, E)`.
                let pv = BRepTool::parameter_of_vertex(e, &v);

                // `End.cxx:280-287`: `Pv = C3D->Value(pv); TP(++nP2) = Pv`.
                if let Some(c3d) = BRepTool::edge_curve(&edge) {
                    tp.push(c3d.value(pv));
                }

                // `End.cxx:291-325`
                let lof = idmoelof.find_from_key(e).to_vec();
                for f in &lof {
                    let Some(face) = Face::wrap(f.clone()) else {
                        continue;
                    };
                    // `End.cxx:294-308`
                    if !boptools_2d::has_curve_on_surface(&edge, &face) {
                        let Some(c2d) = fc2d_curve_on_surface(e, f) else {
                            // `End.cxx:300-303` raises
                            // `Standard_ProgramError("...::End 1")`; the port
                            // skips the face, the only recovery a caller of
                            // `End` can observe.
                            continue;
                        };
                        // `B.UpdateEdge(E, C2D, F, max(Tolerance(E), tolpc))`.
                        // UNPORTED: the `tolpc` out-parameter — it is
                        // `BRep_Tool::Tolerance(E)` for the stored arm
                        // (`2d.cxx:163`) but the projection tolerance of
                        // `FC2D_MakeCurveOnSurface` (`2d.cxx:342-353`) for the
                        // built arm, and `boptools_2d::make_2d` does not
                        // surface it. The p-curve is attached with the edge
                        // tolerance unchanged.
                        let _ = boptools_2d::attach_existing_pcurve(&edge, &face, c2d);
                    }
                    // `End.cxx:309-312`
                    let Some(c2d) = fc2d_curve_on_surface(e, f) else {
                        continue;
                    };
                    let p2 = c2d.d0(pv);
                    let Ok(point) = boptools_2d::face_point_at_uv(&face, p2.x(), p2.y()) else {
                        continue;
                    };
                    tp.push(point);

                    // `End.cxx:313-323` (NIZHNY-MKK, Fri Sep 29 16:08:28 2000)
                    if map_of_new_edges.contains(&shape_map_key(e)) {
                        let an_edge_tol = BRepTool::edge_tolerance(&edge);
                        let a_face_tol = BRepTool::face_tolerance(&face);
                        if an_edge_tol < a_face_tol {
                            GeometryRegistry::global().raise_edge_tolerance(e, a_face_tol);
                        }
                    }
                }

                // `End.cxx:326-334` (NIZHNY-MKK, Fri Sep 29 16:54:08 2000)
                if map_of_new_vertices.contains(&shape_map_key(&v)) {
                    let a_vertex_tol = BRepTool::vertex_tolerance(&vertex);
                    let an_edge_tol = BRepTool::edge_tolerance(&edge);
                    if a_vertex_tol < an_edge_tol {
                        GeometryRegistry::global().raise_vertex_tolerance(&v, an_edge_tol);
                    }
                }
            }

            // `End.cxx:337-364`
            let newtol = BRepTool::vertex_tolerance(&Vertex(v.clone()));
            let mut bx = BndBox::new();
            // `End.cxx:341-342`: `BOX.Set(BRep_Tool::Pnt(V))`.
            bx.add_point(&BRepTool::vertex_point(&Vertex(v.clone())));
            for p in &tp {
                bx.add_point(p);
            }
            let Some((xmin, xmax, ymin, ymax, zmin, zmax)) = bx.get() else {
                continue;
            };
            // `End.cxx:349-352`: `P1(Get(xmin,ymin,zmin))`,
            // `P2(Get(xmax,ymax,zmax))`; `Bnd_Box::Get` emits
            // `(Xmin, Xmax, Ymin, Ymax, Zmin, Zmax)` in the port.
            let p1 = GpPnt::new(xmin, ymin, zmin);
            let p2 = GpPnt::new(xmax, ymax, zmax);
            let d = p1.distance(&p2);
            if d > newtol {
                // `End.cxx:360-361`
                GeometryRegistry::global().raise_vertex_tolerance(&v, d);
            }
        }
    }
}
