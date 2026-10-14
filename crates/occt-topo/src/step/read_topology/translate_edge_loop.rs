//! `StepToTopoDS_TranslateEdgeLoop`: `EDGE_LOOP` / `VERTEX_LOOP` / `FACE_BOUND` resolution and the two vertex-binding passes (`cxx:288-403`, `:405-491`).
use super::*;

impl<'a> Resolver<'a> {
    pub(in crate::step) fn resolve_loop(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let items = parse_ref_list(&rec.args[1]);
        // `StepToTopoDS_TranslateEdgeLoop.cxx:288-403` and `:405-491` bind
        // (confuse) the loop's vertices before the oriented edges are mapped, so a
        // loop whose edges reference distinct `VERTEX_POINT` entities at the same
        // point — or adjacent edges that share no vertex at all — still comes out
        // connected.
        self.bind_edge_loop_vertices(&items);
        let mut edges = Vec::with_capacity(items.len());
        for &it in &items {
            let s = self.resolve_shape(it)?;
            if !s.is_edge() {
                return Err(format!("#{it}: expected EDGE in EDGE_LOOP"));
            }
            edges.push(Edge(s));
        }
        // Keep EDGE_LOOP order for seam SelectForwardSeam (`TranslateEdgeLoop`
        // iterates OrientedEdges in list order; last UpdateEdge wins). Reorder
        // for connectivity after pcurve association in `resolve_face`.
        Ok(self.b.make_wire(&edges).0)
    }

    /// The `EDGE_CURVE` record an `ORIENTED_EDGE` reaches — following a nested
    /// `ORIENTED_EDGE` (bug #29979, `StepToTopoDS_TranslateEdgeLoop.cxx:301-307`)
    /// — together with the occurrence `Orientation` (`:429-433`).
    pub(super) fn oriented_edge_parts(&self, id: usize) -> Option<(usize, bool)> {
        let mut cur = id;
        for _ in 0..8 {
            let rec = self.record(cur).ok()?;
            if rec.type_name != "ORIENTED_EDGE" {
                return None;
            }
            let ori = parse_logical(rec.args.get(4).map(String::as_str), true);
            let next = parse_ref(&rec.args[3])?;
            let nrec = self.record(next).ok()?;
            if nrec.type_name == "ORIENTED_EDGE" {
                cur = next;
                continue;
            }
            if nrec.type_name != "EDGE_CURVE" {
                return None;
            }
            return Some((next, ori));
        }
        None
    }

    /// The two vertex-binding passes of `StepToTopoDS_TranslateEdgeLoop`:
    ///
    /// 1. `cxx:288-403` (bug PRO7656): for every oriented edge whose
    ///    `EDGE_CURVE`'s start/end vertices translate to points within
    ///    `Precision::Confusion()`, bind one STEP vertex onto the other's
    ///    `TopoDS_Vertex` — `Vend → V1` when `Vend` was not yet bound, else
    ///    `Vstart → V2` when `Vstart` was not, else `Vend → V1` again
    ///    (`cxx:384-396`).
    /// 2. `cxx:405-491` (bug BUC50070 #3815): for every adjacent pair `(j, j+1)`
    ///    whose *meeting* vertices (`OrEdge1->Orientation() ? EdgeEnd : EdgeStart`
    ///    and the mirror for edge 2, `cxx:429-433`) are different STEP entities
    ///    and translate to points within `Precision()` (`cxx:464-478`), bind
    ///    `Vs1 → V2` when `EC1` is not yet translated, else `Vs2 → V1` when `EC2`
    ///    is not.
    ///
    /// The port's "translate tool" is `shape_cache` (an entity is bound once
    /// `resolve_shape` has cached it), and the binding itself is
    /// [`Resolver::vertex_bind`], consulted by `resolve_shape` for `VERTEX_POINT`.
    pub(super) fn bind_edge_loop_vertices(&self, oriented: &[usize]) {
        let point_of = |v: &TopoShape| BRepTool::vertex_point(&Vertex(v.clone()));
        // Pass 1.
        for &oe in oriented {
            let Some((ec_id, _)) = self.oriented_edge_parts(oe) else {
                continue;
            };
            let Ok(ec) = self.record(ec_id) else { continue };
            let (Some(a), Some(b)) = (parse_ref(&ec.args[1]), parse_ref(&ec.args[2])) else {
                continue;
            };
            let same_sense = parse_logical(ec.args.get(4).map(String::as_str), true);
            // `cxx:355-365`.
            let (vstart, vend) = if same_sense { (a, b) } else { (b, a) };
            // `cxx:367-368`: the pre-call bound state.
            let ist_v = self.shape_cache.borrow().contains_key(&vstart);
            let ise_v = self.shape_cache.borrow().contains_key(&vend);
            let (Ok(v1), Ok(v2)) = (self.resolve_shape(vstart), self.resolve_shape(vend)) else {
                continue;
            };
            if point_of(&v1).distance(&point_of(&v2)) <= occt_core::precision::CONFUSION {
                let (from, to) = if !ise_v {
                    (vend, v1)
                } else if !ist_v {
                    (vstart, v2)
                } else {
                    (vend, v1)
                };
                self.vertex_bind.borrow_mut().insert(from, to);
            }
        }
        // Pass 2.
        let n = oriented.len();
        for j in 0..n {
            if n < 2 {
                break;
            }
            let Some((ec1_id, ori1)) = self.oriented_edge_parts(oriented[j]) else {
                continue;
            };
            let Some((ec2_id, ori2)) = self.oriented_edge_parts(oriented[(j + 1) % n]) else {
                continue;
            };
            let (Ok(ec1), Ok(ec2)) = (self.record(ec1_id), self.record(ec2_id)) else {
                continue;
            };
            let (Some(s1), Some(e1)) = (parse_ref(&ec1.args[1]), parse_ref(&ec1.args[2])) else {
                continue;
            };
            let (Some(s2), Some(e2)) = (parse_ref(&ec2.args[1]), parse_ref(&ec2.args[2])) else {
                continue;
            };
            // `cxx:429-433`.
            let vs1 = if ori1 { e1 } else { s1 };
            let vs2 = if ori2 { s2 } else { e2 };
            let vs11 = if ori1 { s1 } else { e1 };
            let vs22 = if ori2 { e2 } else { s2 };
            // `cxx:435-438`: already the same STEP vertex on one of the four ends.
            if vs1 == vs2 || vs1 == vs22 || vs2 == vs11 || vs22 == vs11 {
                continue;
            }
            let (Ok(v1), Ok(v2)) = (self.resolve_shape(vs1), self.resolve_shape(vs2)) else {
                continue;
            };
            if Arc::ptr_eq(&v1.tshape, &v2.tshape) {
                continue;
            }
            if point_of(&v1).distance(&point_of(&v2)) <= self.precision {
                // `cxx:466-477`.
                if !self.shape_cache.borrow().contains_key(&ec1_id) {
                    self.vertex_bind.borrow_mut().insert(vs1, v2);
                } else if !self.shape_cache.borrow().contains_key(&ec2_id) {
                    self.vertex_bind.borrow_mut().insert(vs2, v1);
                }
            }
        }
    }

    /// A `VERTEX_LOOP(name, vertex)` is the boundary of a degenerate face — a
    /// loop reduced to a single vertex (sphere pole, cone apex). It maps to a
    /// wire containing no edges (the face still references it via FACE_BOUND).
    pub(in crate::step) fn resolve_vertex_loop(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let _ = parse_ref(&rec.args[1]).ok_or("VERTEX_LOOP: bad vertex ref")?;
        Ok(self.b.make_wire(&[]).0)
    }

    pub(in crate::step) fn resolve_outer_bound(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let loop_ref = parse_ref(&rec.args[1]).ok_or("FACE_OUTER_BOUND: bad loop ref")?;
        let mut s = self.resolve_shape(loop_ref)?;
        if !s.is_wire() {
            return Err("FACE_OUTER_BOUND: loop is not a wire".into());
        }
        // `StepShape_FaceBound::Orientation`: `.F.` reverses the wire
        // (`StepToTopoDS_TranslateFace.cxx` PolyLoop / EdgeLoop bind).
        if !parse_logical(rec.args.get(2).map(String::as_str), true) {
            s.set_orientation(Orientation::Reversed);
        }
        Ok(s)
    }

}
