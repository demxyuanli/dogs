//! `ShapeAnalysis_FreeBounds` (`ShapeAnalysis_FreeBounds.cxx:62-733`).

use super::*;

/// `ShapeAnalysis_FreeBounds` (`ShapeAnalysis_FreeBounds.hxx:58-233`).
pub struct ShapeAnalysisFreeBounds {
    /// `myTolerance`.
    tolerance: f64,
    /// `myShared`.
    shared: bool,
    /// `mySplitClosed`.
    split_closed: bool,
    /// `mySplitOpen`.
    split_open: bool,
    /// `myWires`: the closed wires.
    wires: Compound,
    /// `myEdges`: the open wires.
    edges: Compound,
}

impl Default for ShapeAnalysisFreeBounds {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeAnalysisFreeBounds {
    /// `ShapeAnalysis_FreeBounds()` (`cxx:58`).
    pub fn new() -> Self {
        Self {
            tolerance: 0.0,
            shared: false,
            split_closed: false,
            split_open: false,
            wires: Compound::new(),
            edges: Compound::new(),
        }
    }

    /// `ShapeAnalysis_FreeBounds(shape, toler, splitclosed = false,
    /// splitopen = true)` (`cxx:62-95`): the sewing forecast.
    pub fn forecast(
        shape: &TopoShape,
        toler: f64,
        split_closed: bool,
        split_open: bool,
    ) -> Self {
        let mut fb = Self {
            tolerance: toler,
            shared: false,
            split_closed,
            split_open,
            wires: Compound::new(),
            edges: Compound::new(),
        };

        // `cxx:70-77`: `BRepBuilderAPI_Sewing(toler, false, false)`.
        let mut sew = Sewing::new(toler, false, false, false, false);
        for sub in stored_children(shape) {
            sew.add(&sub);
        }
        let _ = sew.perform();

        // `cxx:79-89`: the free edges, degenerated ones dropped.
        let mut edges: Vec<Edge> = Vec::new();
        for iedge in 1..=sew.nb_free_edges() {
            let Some(edge) = sew.free_edge(iedge) else {
                continue;
            };
            if !BRepTool::is_degenerated(&edge) {
                edges.push(edge);
            }
        }
        // `cxx:90-94`.
        let wires = connect_edges_to_wires(&mut edges, toler, false);
        fb.dispatch_wires(&wires);
        fb.split_wires();
        fb
    }

    /// `ShapeAnalysis_FreeBounds(shape, splitclosed = false, splitopen = true,
    /// checkinternaledges = false)` (`cxx:100-131`): the actual free bounds of
    /// a compound of shells.
    pub fn of_shells(
        shape: &TopoShape,
        split_closed: bool,
        split_open: bool,
        check_internal_edges: bool,
    ) -> Self {
        let mut fb = Self {
            tolerance: 0.0,
            shared: true,
            split_closed,
            split_open,
            wires: Compound::new(),
            edges: Compound::new(),
        };

        // `cxx:109-118`: a shell holding every face of the input.
        let builder = crate::builder::TopoBuilder::new();
        let mut tmp_shell = Shell::new();
        for face in crate::topo_tools::explorer_shapes(shape, ShapeType::Face) {
            if let Some(face) = crate::shape::Face::wrap(face) {
                builder.add_face(&mut tmp_shell, &face);
            }
        }
        let mut sas = ShapeAnalysisShell::new();
        sas.check_oriented_shells(&tmp_shell.0, true, check_internal_edges);

        // `cxx:120-130`.
        if sas.has_free_edges() {
            let free = sas.free_edges();
            let mut edges: Vec<Edge> = seq_from_compound(&free.0, false)
                .into_iter()
                .filter_map(Edge::wrap)
                .collect();
            let wires = connect_edges_to_wires(&mut edges, CONFUSION, true);
            fb.dispatch_wires(&wires);
            fb.split_wires();
        }
        fb
    }

    /// `GetClosedWires()` (`hxx:118-121`).
    pub fn closed_wires(&self) -> &Compound {
        &self.wires
    }

    /// `GetOpenWires()` (`hxx:124-127`).
    pub fn open_wires(&self) -> &Compound {
        &self.edges
    }

    /// `DispatchWires(wires, closed, open)` (`cxx:609-638`).
    fn dispatch_wires(&mut self, wires: &[Wire]) {
        let builder = crate::builder::TopoBuilder::new();
        for wire in wires {
            if wire.0.is_closed() {
                builder.add_compound(&mut self.wires, &wire.0);
            } else {
                builder.add_compound(&mut self.edges, &wire.0);
            }
        }
    }

    /// `SplitWires()` (`cxx:647-690`).
    fn split_wires(&mut self) {
        if !self.split_closed && !self.split_open {
            return;
        }
        let closed_wires: Vec<TopoShape> = seq_from_compound(&self.wires.0, false);
        let open_wires: Vec<TopoShape> = seq_from_compound(&self.edges.0, false);

        let (cw1, ow1) = if self.split_closed {
            split_wires_impl(&closed_wires, self.tolerance, self.shared)
        } else {
            (closed_wires, Vec::new())
        };
        let (cw2, ow2) = if self.split_open {
            split_wires_impl(&open_wires, self.tolerance, self.shared)
        } else {
            (Vec::new(), open_wires)
        };

        // `cxx:680-690`.
        let mut closed = cw1;
        closed.extend(cw2);
        let mut open = ow1;
        open.extend(ow2);
        self.wires = compound_from_seq(&closed);
        self.edges = compound_from_seq(&open);
    }
}
