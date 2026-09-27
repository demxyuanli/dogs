//! Port of ShapeFix_ComposeShell and its auxiliary classes.
//!
//! Currently ported: ShapeFix_WireSegment (ShapeFix_WireSegment.hxx/.cxx),
//! the data class ShapeFix_ComposeShell stores its wire segments in.

use crate::abs::Orientation;
use crate::shape::{Edge, Vertex};
use crate::topo_tools_full::{edge_vertices, is_same};

/// ShapeFix_WireSegment.cxx:116-117.
pub(super) const MININD: i32 = -32000;
pub(super) const MAXIND: i32 = 32000;

/// Port of ShapeFix_WireSegment: a segment of a wire (or a whole wire) plus
/// the patch indices used by ComposeShell. The port keeps the ordered stored
/// edges in a vector (the role of ShapeExtend_WireData) and follows OCCT's
/// 1-based edge indexing in the API: index 0 means "append at the end".
#[derive(Clone, Debug)]
pub struct WireSegment {
    /// myWire: the ordered stored edges.
    edges: Vec<Edge>,
    /// myWire->ManifoldMode().
    manifold: bool,
    /// myVertex.
    vertex: Option<Vertex>,
    /// myOrient.
    orient: Orientation,
    /// myIUMin / myIUMax / myIVMin / myIVMax, indexed like "edges".
    iu_min: Vec<i32>,
    iu_max: Vec<i32>,
    iv_min: Vec<i32>,
    iv_max: Vec<i32>,
}

impl Default for WireSegment {
    fn default() -> Self {
        Self::new()
    }
}

impl WireSegment {
    /// ShapeFix_WireSegment::ShapeFix_WireSegment() (cxx:25-29): empty segment,
    /// orientation FORWARD (Clear leaves myOrient alone, the ctor sets it).
    pub fn new() -> Self {
        let mut s = Self {
            edges: Vec::new(),
            manifold: false,
            vertex: None,
            orient: Orientation::Forward,
            iu_min: Vec::new(),
            iu_max: Vec::new(),
            iv_min: Vec::new(),
            iv_max: Vec::new(),
        };
        s.clear();
        s
    }

    /// ShapeFix_WireSegment(wire, ori) (cxx:33-38).
    pub fn with_edges(edges: Vec<Edge>, orient: Orientation) -> Self {
        let mut s = Self::new();
        s.load_edges(&edges, false);
        s.orient = orient;
        s
    }

    /// Clear() (cxx:42-51): an empty wire with ManifoldMode false, empty index
    /// sequences and a null vertex.
    pub fn clear(&mut self) {
        self.edges.clear();
        self.manifold = false;
        self.iu_min.clear();
        self.iu_max.clear();
        self.iv_min.clear();
        self.iv_max.clear();
        self.vertex = None;
    }

    /// Load(wire) (cxx:55-64): copy every edge of the source, keeping its
    /// ManifoldMode (cxx:59), and give each edge the full index range
    /// (AddEdge(i, edge) - cxx:62).
    pub fn load_edges(&mut self, edges: &[Edge], manifold: bool) {
        self.manifold = manifold;
        let mut i = 1;
        for e in edges {
            self.add_edge_patch(i, e.clone(), MININD, MAXIND, MININD, MAXIND);
            i += 1;
        }
    }

    /// WireData() (cxx:68-71).
    pub fn edges(&self) -> &[Edge] {
        &self.edges
    }

    /// myWire->ManifoldMode().
    pub fn manifold_mode(&self) -> bool {
        self.manifold
    }

    /// myWire->ManifoldMode() = value.
    pub fn set_manifold_mode(&mut self, m: bool) {
        self.manifold = m;
    }

    /// Orientation(ori) (cxx:75-78).
    pub fn set_orientation(&mut self, ori: Orientation) {
        self.orient = ori;
    }

    /// Orientation() (cxx:82-85).
    pub fn orientation(&self) -> Orientation {
        self.orient
    }

    /// FirstVertex() (cxx:89-93): ShapeAnalysis_Edge::FirstVertex of the first
    /// edge, i.e. the child stored FORWARD.
    pub fn first_vertex(&self) -> Option<Vertex> {
        self.edges.first().and_then(|e| edge_vertices(e).0)
    }

    /// LastVertex() (cxx:97-101).
    pub fn last_vertex(&self) -> Option<Vertex> {
        self.edges.last().and_then(|e| edge_vertices(e).1)
    }

    /// IsClosed() (cxx:105-110): FirstVertex().IsSame(LastVertex()).
    pub fn is_closed(&self) -> bool {
        match (self.first_vertex(), self.last_vertex()) {
            (Some(a), Some(b)) => is_same(&a.0, &b.0),
            _ => false,
        }
    }

    /// NbEdges() (cxx:121-124).
    pub fn nb_edges(&self) -> usize {
        self.edges.len()
    }

    /// Edge(i) (cxx:128-132), i is 1-based.
    pub fn edge(&self, i: usize) -> Option<&Edge> {
        self.edges.get(i.checked_sub(1)?)
    }

    /// SetEdge(i, edge) (cxx:136-139).
    pub fn set_edge(&mut self, i: usize, edge: Edge) {
        if let Some(slot) = i.checked_sub(1).and_then(|k| self.edges.get_mut(k)) {
            *slot = edge;
        }
    }

    /// AddEdge(i, edge) (cxx:143-146): the MININD..MAXIND range.
    pub fn add_edge(&mut self, i: usize, edge: Edge) {
        self.add_edge_patch(i, edge, MININD, MAXIND, MININD, MAXIND);
    }

    /// AddEdge(i, edge, iumin, iumax, ivmin, ivmax) (cxx:150-172):
    /// myWire->Add(edge, i) then the four index sequences; i == 0 appends.
    pub fn add_edge_patch(
        &mut self,
        i: usize,
        edge: Edge,
        iumin: i32,
        iumax: i32,
        ivmin: i32,
        ivmax: i32,
    ) {
        if i == 0 {
            self.edges.push(edge);
            self.iu_min.push(iumin);
            self.iu_max.push(iumax);
            self.iv_min.push(ivmin);
            self.iv_max.push(ivmax);
        } else {
            let k = (i - 1).min(self.edges.len());
            self.edges.insert(k, edge);
            self.iu_min.insert(k, iumin);
            self.iu_max.insert(k, iumax);
            self.iv_min.insert(k, ivmin);
            self.iv_max.insert(k, ivmax);
        }
    }

    /// SetPatchIndex(i, ...) (cxx:176-186).
    pub fn set_patch_index(&mut self, i: usize, iumin: i32, iumax: i32, ivmin: i32, ivmax: i32) {
        if let Some(k) = i.checked_sub(1) {
            if k < self.edges.len() {
                self.iu_min[k] = iumin;
                self.iu_max[k] = iumax;
                self.iv_min[k] = ivmin;
                self.iv_max[k] = ivmax;
            }
        }
    }

    /// DefineIUMin(i, iumin) (cxx:190-200): raise only.
    pub fn define_iu_min(&mut self, i: usize, iumin: i32) {
        if let Some(k) = i.checked_sub(1) {
            if k < self.edges.len() && self.iu_min[k] < iumin {
                self.iu_min[k] = iumin;
            }
        }
    }

    /// DefineIUMax(i, iumax) (cxx:204-215): lower only.
    pub fn define_iu_max(&mut self, i: usize, iumax: i32) {
        if let Some(k) = i.checked_sub(1) {
            if k < self.edges.len() && self.iu_max[k] > iumax {
                self.iu_max[k] = iumax;
            }
        }
    }

    /// DefineIVMin(i, ivmin) (cxx:219-230).
    pub fn define_iv_min(&mut self, i: usize, ivmin: i32) {
        if let Some(k) = i.checked_sub(1) {
            if k < self.edges.len() && self.iv_min[k] < ivmin {
                self.iv_min[k] = ivmin;
            }
        }
    }

    /// DefineIVMax(i, ivmax) (cxx:234-245).
    pub fn define_iv_max(&mut self, i: usize, ivmax: i32) {
        if let Some(k) = i.checked_sub(1) {
            if k < self.edges.len() && self.iv_max[k] > ivmax {
                self.iv_max[k] = ivmax;
            }
        }
    }

    /// GetPatchIndex(i, ...) (cxx:249-259).
    pub fn get_patch_index(&self, i: usize) -> Option<(i32, i32, i32, i32)> {
        let k = i.checked_sub(1)?;
        if k >= self.edges.len() {
            return None;
        }
        Some((self.iu_min[k], self.iu_max[k], self.iv_min[k], self.iv_max[k]))
    }

    /// CheckPatchIndex(i) (cxx:263-274): dU and dV must each be 0 or 1.
    pub fn check_patch_index(&self, i: usize) -> bool {
        match self.get_patch_index(i) {
            Some((iun, iux, ivn, ivx)) => {
                let du = iux - iun;
                let dv = ivx - ivn;
                (du == 0 || du == 1) && (dv == 0 || dv == 1)
            }
            None => false,
        }
    }

    /// SetVertex(theVertex) (cxx:278-282).
    pub fn set_vertex(&mut self, v: Option<Vertex>) {
        self.vertex = v;
    }

    /// GetVertex() (cxx:300-303).
    pub fn get_vertex(&self) -> Option<&Vertex> {
        self.vertex.as_ref()
    }

    /// IsVertex() (cxx:307-310).
    pub fn is_vertex(&self) -> bool {
        self.vertex.is_some()
    }
}
