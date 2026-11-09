//! `ShapeExtend_WireData` and the `ShapeAnalysis_Wire` slice the free-bounds
//! search drives.

use super::*;

/// `ShapeExtend_WireData` (`ShapeExtend_WireData.hxx:30-260`), the sequence of
/// edges a wire is flattened into.
#[derive(Clone)]
pub(crate) struct WireData {
    /// `myEdges`.
    edges: Vec<Edge>,
    /// `myNonmanifoldEdges`.
    non_manifold: Vec<Edge>,
    /// `myManifoldMode`.
    manifold_mode: bool,
}

impl WireData {
    /// `ShapeExtend_WireData()` (`cxx:40-43`) = `Clear()` (`cxx:152-159`).
    pub(super) fn new() -> Self {
        Self {
            edges: Vec::new(),
            non_manifold: Vec::new(),
            manifold_mode: true,
        }
    }

    /// `ShapeExtend_WireData(wire, chained = false, theManifold = true)`
    /// (`cxx:47-52`) = `Init(wire, chained, theManifold)` (`cxx:74-148`).
    pub(super) fn from_wire(wire: &Wire, chained: bool, the_manifold: bool) -> Self {
        let mut wd = Self::new();
        wd.manifold_mode = the_manifold;
        let mut ok = true;
        let mut vlast: Option<Vertex> = None;

        for child in cumulated_children(&wire.0) {
            if child.shape_type() != ShapeType::Edge {
                continue;
            }
            let edge = Edge(child);
            let ori = edge.0.orientation();
            if ori == Orientation::Internal || ori == Orientation::External {
                wd.non_manifold.push(edge); // `cxx:85-89`.
                continue;
            }
            // `cxx:91-103`: `TopoDS_Iterator itv(E)` keeps `cumOri` true.
            let mut v1: Option<Vertex> = None;
            let mut v2: Option<Vertex> = None;
            for v in cumulated_children(&edge.0) {
                if v.shape_type() != ShapeType::Vertex {
                    continue;
                }
                match v.orientation() {
                    Orientation::Forward => v1 = Some(Vertex(v)),
                    Orientation::Reversed => v2 = Some(Vertex(v)),
                    _ => {}
                }
            }
            // `cxx:106-113`.
            let mismatched = match (vlast.as_ref(), v1.as_ref()) {
                (Some(prev), Some(cur)) => !is_same(&prev.0, &cur.0),
                (Some(_), None) => true,
                _ => false,
            };
            if mismatched && the_manifold {
                ok = false;
                if !chained {
                    break;
                }
            }
            vlast = v2;
            // `cxx:115-122`.
            if wire.0.orientation() == Orientation::Reversed {
                wd.edges.insert(0, edge);
            } else {
                wd.edges.push(edge);
            }
        }

        // `cxx:125-134`.
        if !wd.manifold_mode {
            let extra = std::mem::take(&mut wd.non_manifold);
            wd.edges.extend(extra);
        }
        if ok || chained {
            return wd;
        }

        // `cxx:136-147`: rechaine through `BRepTools_WireExplorer`.
        wd.clear();
        let mut we = WireExplorer::of_wire(wire);
        while we.more() {
            if let Some(e) = we.current_edge() {
                wd.edges.push(e);
            }
            we.next();
        }
        wd
    }

    /// `Clear()` (`cxx:152-159`).
    pub(super) fn clear(&mut self) {
        self.edges.clear();
        self.non_manifold.clear();
        self.manifold_mode = true;
    }

    /// `NbEdges()` (`cxx:576-579`).
    pub(super) fn nb_edges(&self) -> usize {
        self.edges.len()
    }

    /// `NbNonManifoldEdges()` (`cxx:596-599`).
    pub(super) fn nb_non_manifold_edges(&self) -> usize {
        self.non_manifold.len()
    }

    /// `Edge(num)` (`cxx:583-592`): 1-based; a negative `num` answers the
    /// reversed edge of `-num`.
    pub(super) fn edge(&self, num: isize) -> Option<Edge> {
        if num < 0 {
            let mut e = self.edge(-num)?;
            e.0.reverse();
            return Some(e);
        }
        let e = self.edges.get(num as usize - 1)?;
        Some(e.clone())
    }

    /// `NonmanifoldEdge(num)` (`cxx:603-614`).
    pub(super) fn _nonmanifold_edge(&self, num: usize) -> Option<&Edge> {
        self.non_manifold.get(num.checked_sub(1)?)
    }

    /// `ManifoldMode()` (`hxx:`): the mutable accessor
    /// `sewd->ManifoldMode() = isUsedManifoldMode` (`cxx:398`).
    pub(super) fn set_manifold_mode(&mut self, mode: bool) {
        self.manifold_mode = mode;
    }

    /// `Add(edge, atnum)` (`cxx:259-281`). The only OCCT caller is
    /// `ShapeFix_Wire` (`Insert` / `Add`), so nothing on the free-bounds path
    /// reaches it; kept because the class owns the overload.
    #[allow(dead_code)]
    pub(super) fn add_edge(&mut self, edge: &Edge, atnum: usize) {
        let ori = edge.0.orientation();
        if ori != Orientation::Reversed && ori != Orientation::Forward && self.manifold_mode {
            self.non_manifold.push(edge.clone()); // `cxx:261-266`.
            return;
        }
        if atnum == 0 {
            self.edges.push(edge.clone()); // `cxx:272-275`.
        } else {
            self.edges.insert(atnum - 1, edge.clone()); // `cxx:276-279`.
        }
    }

    /// `Add(wire, atnum)` (`cxx:285-324`): every edge of the wire is inserted
    /// in order, so `atnum = 1` prepends the whole wire.
    pub(super) fn add_wire(&mut self, wire: &Wire, atnum: usize) {
        let mut n = atnum;
        let mut nm_edges: Vec<Edge> = Vec::new();
        for child in cumulated_children(&wire.0) {
            if child.shape_type() != ShapeType::Edge {
                continue;
            }
            let edge = Edge(child);
            let ori = edge.0.orientation();
            if ori != Orientation::Reversed && ori != Orientation::Forward {
                if self.manifold_mode {
                    self.non_manifold.push(edge); // `cxx:298-301`.
                } else {
                    nm_edges.push(edge); // `cxx:302-305`.
                }
                continue;
            }
            if n == 0 {
                self.edges.push(edge);
            } else {
                self.edges.insert(n - 1, edge);
                n += 1;
            }
        }
        // `cxx:318-322`.
        self.edges.extend(nm_edges);
    }

    /// `Add(wireData, atnum)` (`cxx:328-384`).
    pub(super) fn add_wiredata(&mut self, wire: &WireData, atnum: usize) {
        let mut n = atnum;
        let mut nm_edges: Vec<Edge> = Vec::new();
        for i in 1..=wire.nb_edges() {
            let Some(edge) = wire.edge(i as isize) else {
                continue;
            };
            let ori = edge.0.orientation();
            if ori == Orientation::Internal || ori == Orientation::External {
                nm_edges.push(edge); // `cxx:340-344`.
                continue;
            }
            if n == 0 {
                self.edges.push(edge);
            } else {
                self.edges.insert(n - 1, edge);
                n += 1;
            }
        }
        // `cxx:357-361`.
        self.edges.extend(nm_edges);
        // `cxx:363-381`. OCCT reads `wire->Edge(i)` in the non-manifold arm
        // rather than `wire->NonmanifoldEdge(i)`; the arm is unreachable from
        // `ShapeAnalysis_FreeBounds` (which always keeps `myManifoldMode`).
        for i in 1..=wire.nb_non_manifold_edges() {
            if self.manifold_mode {
                if let Some(e) = wire._nonmanifold_edge(i) {
                    self.non_manifold.push(e.clone());
                }
            } else if n == 0 {
                if let Some(e) = wire.edge(i as isize) {
                    self.edges.push(e);
                }
            } else {
                if let Some(e) = wire.edge(i as isize) {
                    self.edges.insert(n - 1, e);
                    n += 1;
                }
            }
        }
    }

    /// `Set(edge, num)` (`cxx:456-476`): `num <= 0` replaces the last edge.
    pub(super) fn set(&mut self, edge: &Edge, num: usize) {
        let ori = edge.0.orientation();
        if ori != Orientation::Reversed && ori != Orientation::Forward && self.manifold_mode {
            if num <= self.non_manifold.len() {
                self.non_manifold[num - 1] = edge.clone(); // `cxx:461-464`.
            } else {
                self.non_manifold.push(edge.clone()); // `cxx:465-468`.
            }
        } else {
            let at = if num > 0 { num } else { self.edges.len() }; // `cxx:473`.
            if at >= 1 && at <= self.edges.len() {
                self.edges[at - 1] = edge.clone();
            }
        }
    }

    /// `Wire()` (`cxx:651-690`): rebuild the wire, marking it closed when its
    /// two endpoint vertices are the same.
    pub(super) fn wire(&self) -> Wire {
        let builder = crate::builder::TopoBuilder::new();
        let mut w = Wire::new();
        let mut is_manifold = true;
        for i in 1..=self.nb_edges() {
            let Some(edge) = self.edge(i as isize) else {
                continue;
            };
            let ori = edge.0.orientation();
            if ori != Orientation::Forward && ori != Orientation::Reversed {
                is_manifold = false;
            }
            builder.add_edge(&mut w, &edge);
        }
        if is_manifold {
            let (vf, vl) = topexp_vertices_of_wire(&w);
            if let (Some(vf), Some(vl)) = (vf, vl) {
                if is_same(&vf.0, &vl.0) {
                    w.0.set_closed(true);
                }
            }
        }
        if self.manifold_mode {
            for edge in &self.non_manifold {
                builder.add_edge(&mut w, edge);
            }
        }
        w
    }
}

/// `ShapeAnalysis_Wire` (`ShapeAnalysis_Wire.hxx`), narrowed to the state
/// `connectWiresToWiresImpl` drives: `Load`, `SetPrecision`, `CheckConnected`
/// and `LastCheckStatus`.
pub(crate) struct WireAnalyzer {
    /// `myPrecision`.
    precision: f64,
    /// `myStatus`.
    status: i32,
    /// `myMin3d`.
    min3d: f64,
}

impl WireAnalyzer {
    /// `ShapeAnalysis_Wire()` : `myPrecision` keeps the default `1e-7` until
    /// `SetPrecision` overwrites it (`ShapeAnalysis_Wire.cxx:201-204`).
    pub(super) fn new() -> Self {
        Self {
            precision: 1e-7,
            status: SE_OK,
            min3d: 0.0,
        }
    }

    /// `SetPrecision(precision)` (`cxx:201-204`).
    pub(super) fn set_precision(&mut self, precision: f64) {
        self.precision = precision;
    }

    /// `LastCheckStatus(test)`.
    pub(super) fn last_check_status(&self, test: i32) -> bool {
        super::super::box_bnd_tree::decode_status(self.status, test)
    }

    /// `CheckConnected(num, prec = 0)` (`cxx:693-761`).
    pub(super) fn check_connected(&mut self, wire: &WireData, num: usize) -> bool {
        self.status = SE_OK;
        let nb = wire.nb_edges();
        if nb < 1 {
            return false; // `cxx:696-699`.
        }
        let n2 = if num > 0 { num } else { nb }; // `cxx:701`.
        let n1 = if n2 > 1 { n2 - 1 } else { nb }; // `cxx:702`.

        let Some(e1) = wire.edge(n1 as isize) else {
            return false;
        };
        let Some(e2) = wire.edge(n2 as isize) else {
            return false;
        };

        let v1 = super::super::last_vertex(&e1);
        let v2 = super::super::first_vertex(&e2);
        let (Some(v1), Some(v2)) = (v1, v2) else {
            self.status = super::super::box_bnd_tree::SE_FAIL2; // `cxx:711-715`.
            return false;
        };
        if is_same(&v1.0, &v2.0) {
            return false; // `cxx:716-719`.
        }

        let p1 = BRepTool::vertex_point(&v1);
        let p2 = BRepTool::vertex_point(&v2);
        self.min3d = p1.distance(&p2);
        if self.min3d <= occt_core::precision::REAL_SMALL {
            self.status |= super::super::box_bnd_tree::SE_DONE1; // `cxx:724-727`.
        } else if self.min3d <= self.precision {
            self.status |= super::super::box_bnd_tree::SE_DONE2; // `cxx:728-731`.
        } else if self.min3d <= 0.0 {
            self.status |= super::super::box_bnd_tree::SE_DONE3; // `cxx:732-735`.
        } else {
            // `cxx:736-758`.
            if n1 == n2 {
                self.status = super::super::box_bnd_tree::SE_FAIL1; // `cxx:739-742`.
            } else {
                let Some(v2b) = super::super::last_vertex(&e2) else {
                    self.status = super::super::box_bnd_tree::SE_FAIL1;
                    return false;
                };
                let p2b = BRepTool::vertex_point(&v2b);
                let dist = p1.distance(&p2b);
                if dist > self.precision {
                    self.status = super::super::box_bnd_tree::SE_FAIL1; // `cxx:748-751`.
                } else {
                    self.min3d = dist;
                    self.status = super::super::box_bnd_tree::SE_FAIL2; // `cxx:754-755`.
                }
            }
            return false;
        }
        true
    }
}
