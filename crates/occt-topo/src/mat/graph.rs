//! `MAT_Graph` (`MAT_Graph.cxx:42-582`, `.hxx:35-130`): the graph of the medial
//! axis, built from the tree of bisectors, plus the query helpers that OCCT
//! implements on the elements themselves (`MAT_Node::LinkedArcs`
//! `MAT_Node.cxx:52-69`, `MAT_Node::NearElts` `:73-107`,
//! `MAT_Node::PendingNode` `:118-122`).
//!
//! `theArcs` / `theNodes` / `theBasicElts` are `NCollection_DataMap<int, ...>`
//! in OCCT: they are keyed by a *mutable* integer, may contain holes, and are
//! compacted by `CompactArcs` / `CompactNodes`. The port keeps one arena per
//! class (a handle is the arena slot) plus the same key→slot maps, so
//! `UnBind` / `Bind` / compaction behave exactly as in OCCT while the graph
//! cycles stay plain ids.

use std::collections::BTreeMap;

use occt_core::precision::INFINITE;

use super::elements::{
    ArcId, BasicEltId, BisectorId, EdgeId, MatArc, MatBasicElt, MatBisector, MatEdge, MatNode,
    MatSide, NodeId,
};
use super::MatZone;

/// `MAT_Graph` (`MAT_Graph.hxx:36-130`).
#[derive(Debug, Clone, Default)]
pub struct MatGraph {
    /// Slot arena of the arcs; a handle is the slot.
    arc_arena: Vec<MatArc>,
    /// Slot arena of the nodes.
    node_arena: Vec<MatNode>,
    /// Slot arena of the basic elements.
    basic_elt_arena: Vec<MatBasicElt>,
    /// Slot arena of the bisectors (input of `Perform`).
    bisector_arena: Vec<MatBisector>,
    /// Slot arena of the edges (input of `Perform`).
    edge_arena: Vec<MatEdge>,
    /// `theArcs` (`.hxx:124`).
    the_arcs: BTreeMap<i32, ArcId>,
    /// `theBasicElts` (`.hxx:125`).
    the_basic_elts: BTreeMap<i32, BasicEltId>,
    /// `theNodes` (`.hxx:126`).
    the_nodes: BTreeMap<i32, NodeId>,
    /// `numberOfArcs` (`.hxx:127`).
    number_of_arcs: i32,
    /// `numberOfNodes` (`.hxx:128`).
    number_of_nodes: i32,
    /// `numberOfBasicElts` (`.hxx:129`).
    number_of_basic_elts: i32,
    /// `numberOfInfiniteNodes` (`.hxx:130`).
    number_of_infinite_nodes: i32,
}

impl MatGraph {
    /// `MAT_Graph()` (`cxx:43-50`).
    pub fn new() -> Self {
        Self::default()
    }

    // ---------------------------------------------------------------
    // Arena access used by `MAT2d_Mat2d` to fill the input elements.
    // ---------------------------------------------------------------

    /// Appends an edge and returns its handle.
    pub fn add_edge(&mut self, edge: MatEdge) -> EdgeId {
        self.edge_arena.push(edge);
        self.edge_arena.len() - 1
    }

    /// Appends a bisector and returns its handle.
    pub fn add_bisector(&mut self, bisector: MatBisector) -> BisectorId {
        self.bisector_arena.push(bisector);
        self.bisector_arena.len() - 1
    }

    /// Appends a basic element and returns its handle. The caller sets the
    /// `Index()` key it is registered under through [`MatGraph::bind_basic_elt`].
    pub fn add_basic_elt(&mut self, basic_elt: MatBasicElt) -> BasicEltId {
        self.basic_elt_arena.push(basic_elt);
        self.basic_elt_arena.len() - 1
    }

    /// `theBasicElts.Bind(Index, Elt)`.
    pub fn bind_basic_elt(&mut self, index: i32, elt: BasicEltId) {
        self.the_basic_elts.insert(index, elt);
    }

    /// `theBasicElts.UnBind(Index)`.
    pub fn unbind_basic_elt(&mut self, index: i32) -> Option<BasicEltId> {
        self.the_basic_elts.remove(&index)
    }

    /// A bisector by handle.
    pub fn bisector(&self, id: BisectorId) -> &MatBisector {
        &self.bisector_arena[id]
    }

    /// A mutable bisector by handle.
    pub fn bisector_mut(&mut self, id: BisectorId) -> &mut MatBisector {
        &mut self.bisector_arena[id]
    }

    /// An edge by handle.
    pub fn edge(&self, id: EdgeId) -> &MatEdge {
        &self.edge_arena[id]
    }

    /// A mutable edge by handle.
    pub fn edge_mut(&mut self, id: EdgeId) -> &mut MatEdge {
        &mut self.edge_arena[id]
    }

    /// A basic element by handle.
    pub fn basic_elt_by_id(&self, id: BasicEltId) -> &MatBasicElt {
        &self.basic_elt_arena[id]
    }

    /// A mutable basic element by handle.
    pub fn basic_elt_mut(&mut self, id: BasicEltId) -> &mut MatBasicElt {
        &mut self.basic_elt_arena[id]
    }

    /// An arc by handle.
    pub fn arc_by_id(&self, id: ArcId) -> &MatArc {
        &self.arc_arena[id]
    }

    /// A mutable arc by handle.
    pub fn arc_mut(&mut self, id: ArcId) -> &mut MatArc {
        &mut self.arc_arena[id]
    }

    /// A node by handle.
    pub fn node_by_id(&self, id: NodeId) -> &MatNode {
        &self.node_arena[id]
    }

    /// A mutable node by handle.
    pub fn node_mut(&mut self, id: NodeId) -> &mut MatNode {
        &mut self.node_arena[id]
    }

    /// The basic elements arena (for `MAT2d_Mat2d` bookkeeping).
    pub fn basic_elts(&self) -> &[MatBasicElt] {
        &self.basic_elt_arena
    }

    // ---------------------------------------------------------------
    // MAT_Graph::Perform
    // ---------------------------------------------------------------

    /// `Perform(SemiInfinite, TheRoots, NbBasicElts, NbArcs)` (`cxx:56-171`).
    pub fn perform(
        &mut self,
        semi_infinite: bool,
        the_roots: &[BisectorId],
        nb_basic_elts: i32,
        nb_arcs: i32,
    ) {
        let nb_roots = if semi_infinite {
            self.number_of_infinite_nodes = the_roots.len() as i32;
            the_roots.len() as i32
        } else {
            self.number_of_infinite_nodes = 0;
            1
        };

        self.number_of_arcs = nb_arcs;
        self.number_of_basic_elts = nb_basic_elts;
        self.number_of_nodes = nb_roots + nb_arcs;
        let mut ind_tab_nodes = self.number_of_nodes;

        // `cxx:96-100`: creation of the basic elements.
        for i in 1..=nb_basic_elts {
            let elt = self.add_basic_elt(MatBasicElt::new(i));
            self.basic_elt_mut(elt).set_geom_index(i);
            self.bind_basic_elt(i, elt);
        }

        let mut ind_tab_arcs = 1;

        if semi_infinite {
            // `cxx:108-123`: one entry point per infinite root.
            for &root in the_roots {
                let current_arc = self.make_arc(root, &mut ind_tab_arcs);
                let extremite = self.add_node(MatNode::new(0, Some(current_arc), INFINITE));
                self.node_mut(extremite).set_index(ind_tab_nodes);
                self.arc_mut(current_arc).set_second_node(Some(extremite));
                self.the_nodes.insert(ind_tab_nodes, extremite);
                ind_tab_nodes -= 1;
            }
        } else {
            // `cxx:125-165`: a single entry point and the root node.
            let current_arc = self.make_arc(the_roots[0], &mut ind_tab_arcs);
            let dist_ext = self
                .bisector(the_roots[0])
                .first_edge()
                .map(|e| self.edge(e).distance())
                .unwrap_or(0.0);
            let ind_ext = self.bisector(the_roots[0]).end_point();

            let extremite = self.add_node(MatNode::new(ind_ext, Some(current_arc), dist_ext));
            self.node_mut(extremite).set_index(ind_tab_nodes);
            self.arc_mut(current_arc).set_second_node(Some(extremite));
            self.the_nodes.insert(ind_tab_nodes, extremite);
            ind_tab_nodes -= 1;

            let first_arc = current_arc;
            let mut previous_arc = first_arc;
            let mut last_arc = current_arc;

            for &root in &the_roots[1..] {
                let current = self.make_arc(root, &mut ind_tab_arcs);
                self.arc_mut(current).set_second_node(Some(extremite));
                self.arc_mut(current)
                    .set_neighbour(MatSide::Left, extremite, Some(previous_arc));
                self.arc_mut(previous_arc)
                    .set_neighbour(MatSide::Right, extremite, Some(current));
                previous_arc = current;
                last_arc = current;
            }

            self.arc_mut(first_arc)
                .set_neighbour(MatSide::Left, extremite, Some(last_arc));
            self.arc_mut(last_arc)
                .set_neighbour(MatSide::Right, extremite, Some(first_arc));
        }

        // `cxx:169`: the root nodes are done, the rest follows.
        self.up_date_nodes(&mut ind_tab_nodes);
    }

    fn add_node(&mut self, node: MatNode) -> NodeId {
        self.node_arena.push(node);
        self.node_arena.len() - 1
    }

    /// `MakeArc(aBisector, TheBasicElts, TheArcs, IndTabArcs)` (`cxx:509-581`).
    fn make_arc(&mut self, a_bisector: BisectorId, ind_tab_arcs: &mut i32) -> ArcId {
        let bisector = self.bisector(a_bisector).clone();
        let first_element = bisector
            .first_edge()
            .map(|e| self.edge(e).edge_number())
            .and_then(|n| self.the_basic_elts.get(&n).copied());
        let second_element = bisector
            .second_edge()
            .map(|e| self.edge(e).edge_number())
            .and_then(|n| self.the_basic_elts.get(&n).copied());

        let current_arc = self.add_arc(MatArc::new(
            *ind_tab_arcs,
            bisector.bisector_number(),
            first_element,
            second_element,
        ));

        let mut dist_ext = bisector.dist_issue_point();
        if dist_ext == INFINITE {
            // `cxx:532-539`: OCCT substitutes 1.0 when the issue point has no
            // distance.
            dist_ext = 1.0;
        }

        let extremite = self.add_node(MatNode::new(bisector.issue_point(), Some(current_arc), dist_ext));
        self.arc_mut(current_arc).set_first_node(Some(extremite));

        let sub_list = bisector.list().items().to_vec();
        if sub_list.is_empty() {
            // `cxx:548-556`: the arc lies on the contour.
            if let Some(second) = second_element {
                self.basic_elt_mut(second).set_start_arc(Some(current_arc));
            }
            if let Some(first) = first_element {
                self.basic_elt_mut(first).set_end_arc(Some(current_arc));
            }
        } else {
            let mut prev_arc = current_arc;
            let mut next_arc = current_arc;
            for sub in sub_list {
                next_arc = self.make_arc(sub, ind_tab_arcs);
                self.arc_mut(next_arc).set_second_node(Some(extremite));
                self.arc_mut(next_arc)
                    .set_neighbour(MatSide::Left, extremite, Some(prev_arc));
                self.arc_mut(prev_arc)
                    .set_neighbour(MatSide::Right, extremite, Some(next_arc));
                prev_arc = next_arc;
            }
            self.arc_mut(current_arc)
                .set_neighbour(MatSide::Left, extremite, Some(next_arc));
            self.arc_mut(next_arc)
                .set_neighbour(MatSide::Right, extremite, Some(current_arc));
        }

        self.arc_mut(current_arc).set_index(*ind_tab_arcs);
        self.the_arcs.insert(*ind_tab_arcs, current_arc);
        *ind_tab_arcs += 1;

        current_arc
    }

    fn add_arc(&mut self, arc: MatArc) -> ArcId {
        self.arc_arena.push(arc);
        self.arc_arena.len() - 1
    }

    /// `UpDateNodes(IndTabNodes)` (`cxx:490-503`).
    fn up_date_nodes(&mut self, ind_tab_nodes: &mut i32) {
        for i in 1..=self.number_of_arcs {
            let arc = match self.the_arcs.get(&i) {
                Some(&a) => a,
                // OCCT indexes `theArcs(i)` unconditionally; a missing key is
                // an out-of-range in the map, reported the same way.
                None => panic!("MAT_Graph::UpDateNodes: key {i} is not bound"),
            };
            let bout = self.arc_by_id(arc).first_node();
            if let Some(node) = bout {
                self.the_nodes.insert(*ind_tab_nodes, node);
                self.node_mut(node).set_index(*ind_tab_nodes);
                *ind_tab_nodes -= 1;
                self.node_mut(node).set_linked_arc(Some(arc));
            }
        }
    }

    // ---------------------------------------------------------------
    // Accessors keyed by `Index()`
    // ---------------------------------------------------------------

    /// `Arc(Index)` (`cxx:176-179`).
    pub fn arc(&self, index: i32) -> ArcId {
        *self
            .the_arcs
            .get(&index)
            .unwrap_or_else(|| panic!("MAT_Graph::Arc: key {index} is not bound"))
    }

    /// `BasicElt(Index)` (`cxx:182-185`).
    pub fn basic_elt(&self, index: i32) -> BasicEltId {
        *self
            .the_basic_elts
            .get(&index)
            .unwrap_or_else(|| panic!("MAT_Graph::BasicElt: key {index} is not bound"))
    }

    /// `Node(Index)` (`cxx:188-191`).
    pub fn node(&self, index: i32) -> NodeId {
        *self
            .the_nodes
            .get(&index)
            .unwrap_or_else(|| panic!("MAT_Graph::Node: key {index} is not bound"))
    }

    /// `NumberOfArcs()` (`cxx:194-197`).
    pub fn number_of_arcs(&self) -> i32 {
        self.number_of_arcs
    }

    /// `NumberOfNodes()` (`cxx:201-204`).
    pub fn number_of_nodes(&self) -> i32 {
        self.number_of_nodes
    }

    /// `NumberOfInfiniteNodes()` (`cxx:208-211`).
    pub fn number_of_infinite_nodes(&self) -> i32 {
        self.number_of_infinite_nodes
    }

    /// `NumberOfBasicElts()` (`cxx:215-218`).
    pub fn number_of_basic_elts(&self) -> i32 {
        self.number_of_basic_elts
    }

    /// The bound arc keys, ascending.
    pub fn arc_keys(&self) -> Vec<i32> {
        self.the_arcs.keys().copied().collect()
    }

    /// The bound node keys, ascending.
    pub fn node_keys(&self) -> Vec<i32> {
        self.the_nodes.keys().copied().collect()
    }

    /// The bound basic-element keys, ascending.
    pub fn basic_elt_keys(&self) -> Vec<i32> {
        self.the_basic_elts.keys().copied().collect()
    }

    // ---------------------------------------------------------------
    // MAT_Node queries that need to walk the graph
    // ---------------------------------------------------------------

    /// `MAT_Node::LinkedArcs(S)` (`MAT_Node.cxx:52-69`): the linked arc plus the
    /// arcs reached by turning left around the node.
    pub fn linked_arcs(&self, node: NodeId) -> Vec<ArcId> {
        let mut s = Vec::new();
        let la = match self.node_by_id(node).linked_arc() {
            Some(a) => a,
            None => return s,
        };
        s.push(la);

        if self.arc_by_id(la).has_neighbour(node, MatSide::Left) {
            let mut ca = self.arc_by_id(la).neighbour(node, MatSide::Left).unwrap();
            while ca != la {
                s.push(ca);
                ca = self.arc_by_id(ca).neighbour(node, MatSide::Left).unwrap();
            }
        }
        s
    }

    /// `MAT_Node::NearElts(S)` (`MAT_Node.cxx:73-107`): the two elements of the
    /// linked arc, then every other pair met while turning left.
    pub fn near_elts(&self, node: NodeId) -> Vec<BasicEltId> {
        let mut s = Vec::new();
        let la = match self.node_by_id(node).linked_arc() {
            Some(a) => a,
            None => return s,
        };

        if let Some(e) = self.arc_by_id(la).first_element() {
            s.push(e);
        }
        if let Some(e) = self.arc_by_id(la).second_element() {
            s.push(e);
        }

        if self.arc_by_id(la).has_neighbour(node, MatSide::Left) {
            let mut ca = self.arc_by_id(la).neighbour(node, MatSide::Left).unwrap();
            let mut pair = false;

            while ca != la {
                if pair {
                    if let Some(e) = self.arc_by_id(ca).first_element() {
                        s.push(e);
                    }
                    if let Some(e) = self.arc_by_id(ca).second_element() {
                        s.push(e);
                    }
                } else {
                    pair = true;
                }
                ca = self.arc_by_id(ca).neighbour(node, MatSide::Left).unwrap();
            }
        }
        s
    }

    /// `MAT_Node::PendingNode()` (`MAT_Node.cxx:118-122`): the linked arc has no
    /// left neighbour around the node.
    pub fn pending_node(&self, node: NodeId) -> bool {
        match self.node_by_id(node).linked_arc() {
            Some(la) => !self.arc_by_id(la).has_neighbour(node, MatSide::Left),
            None => true,
        }
    }

    // ---------------------------------------------------------------
    // Fusion / compaction
    // ---------------------------------------------------------------

    /// `FusionOfBasicElts(...)` (`cxx:224-321`).
    #[allow(clippy::too_many_arguments)]
    pub fn fusion_of_basic_elts(
        &mut self,
        index_elt1: i32,
        index_elt2: i32,
        merge_arc1: &mut bool,
        i_geom_arc1: &mut i32,
        i_geom_arc2: &mut i32,
        merge_arc2: &mut bool,
        i_geom_arc3: &mut i32,
        i_geom_arc4: &mut i32,
    ) {
        let elt1 = self.basic_elt(index_elt1);
        let elt2 = self.basic_elt(index_elt2);

        if elt1 == elt2 {
            return;
        }

        // `cxx:241`: the zone of Elt2, walked before anything is modified.
        let zone2 = MatZone::from_basic_elt(self, elt2);

        // `cxx:246-258`: the arcs of that zone no longer separate Elt2.
        for &arc in zone2.frontier() {
            if self.arc_by_id(arc).first_element() == Some(elt2) {
                self.arc_mut(arc).set_first_element(Some(elt1));
            } else {
                self.arc_mut(arc).set_second_element(Some(elt1));
            }
        }

        // `cxx:263-279`: a shared EndArc(1) / StartArc(2) pair is fused.
        let ea1 = self.basic_elt_by_id(elt1).end_arc();
        let sa2 = self.basic_elt_by_id(elt2).start_arc();
        *merge_arc1 = false;

        if let (Some(ea1), Some(sa2)) = (ea1, sa2) {
            let e1 = self.arc_by_id(ea1).first_element();
            let e2 = self.arc_by_id(ea1).second_element();
            let e3 = self.arc_by_id(sa2).first_element();
            let e4 = self.arc_by_id(sa2).second_element();

            if (e1 == e3 || e1 == e4) && (e2 == e3 || e2 == e4) {
                *i_geom_arc1 = self.arc_by_id(ea1).geom_index();
                *i_geom_arc2 = self.arc_by_id(sa2).geom_index();
                self.fusion_of_arcs(ea1, sa2);
                *merge_arc1 = true;
            }
        }

        // `cxx:284`: the end of Elt1 becomes the end of Elt2.
        let elt2_end = self.basic_elt_by_id(elt2).end_arc();
        self.basic_elt_mut(elt1).set_end_arc(elt2_end);

        // `cxx:290-315`: for a closed contour made of a single basic element the
        // EndArc and StartArc of Elt1 also separate the same elements.
        let sa1 = self.basic_elt_by_id(elt1).start_arc();
        let ea1 = self.basic_elt_by_id(elt1).end_arc();

        if ea1 != sa1 {
            if let (Some(ea1), Some(sa1)) = (ea1, sa1) {
                let e1 = self.arc_by_id(ea1).first_element();
                let e2 = self.arc_by_id(ea1).second_element();
                let e3 = self.arc_by_id(sa1).first_element();
                let e4 = self.arc_by_id(sa1).second_element();

                let on_fig = self
                    .arc_by_id(ea1)
                    .first_node()
                    .map(|n| self.node_by_id(n).on_basic_elt())
                    .unwrap_or(false)
                    || self
                        .arc_by_id(ea1)
                        .second_node()
                        .map(|n| self.node_by_id(n).on_basic_elt())
                        .unwrap_or(false)
                    || self
                        .arc_by_id(sa1)
                        .first_node()
                        .map(|n| self.node_by_id(n).on_basic_elt())
                        .unwrap_or(false)
                    || self
                        .arc_by_id(sa1)
                        .second_node()
                        .map(|n| self.node_by_id(n).on_basic_elt())
                        .unwrap_or(false);

                *merge_arc2 = false;

                if (e1 == e3 || e1 == e4) && (e2 == e3 || e2 == e4) && !on_fig {
                    *i_geom_arc3 = self.arc_by_id(ea1).geom_index();
                    *i_geom_arc4 = self.arc_by_id(sa1).geom_index();
                    self.fusion_of_arcs(ea1, sa1);
                    *merge_arc2 = true;
                }
            }
        }

        // `cxx:319-320`: Elt2 is eliminated.
        let elt2_index = self.basic_elt_by_id(elt2).index();
        self.the_basic_elts.remove(&elt2_index);
        self.number_of_basic_elts -= 1;
    }

    /// `FusionOfArcs(Arc1, Arc2)` (`cxx:328-402`): the second node of `Arc2`
    /// becomes the first node of `Arc1`, and `Arc2` is eliminated.
    fn fusion_of_arcs(&mut self, arc1: ArcId, arc2: ArcId) {
        let old_node1 = self.arc_by_id(arc1).first_node();
        let old_node2 = self.arc_by_id(arc2).first_node();
        let arc2_second_node = self.arc_by_id(arc2).second_node();

        self.arc_mut(arc1).set_first_node(arc2_second_node);

        // `cxx:339-353`: neighbourhood around the new first node of Arc1.
        if let Some(node2) = arc2_second_node {
            if !self.node_by_id(node2).infinite() {
                let l_neighbour = self.arc_by_id(arc2).neighbour(node2, MatSide::Left);
                let r_neighbour = self.arc_by_id(arc2).neighbour(node2, MatSide::Right);

                self.arc_mut(arc1).set_first_arc(MatSide::Left, l_neighbour);
                self.arc_mut(arc1).set_first_arc(MatSide::Right, r_neighbour);
                if let Some(l) = l_neighbour {
                    self.arc_mut(l)
                        .set_neighbour(MatSide::Right, node2, Some(arc1));
                }
                if let Some(r) = r_neighbour {
                    self.arc_mut(r)
                        .set_neighbour(MatSide::Left, node2, Some(arc1));
                }
            } else {
                self.arc_mut(arc1).set_first_arc(MatSide::Left, None);
                self.arc_mut(arc1).set_first_arc(MatSide::Right, None);
            }
        }

        // `cxx:359`: the first node of Arc1 now links to Arc1.
        if let Some(node) = self.arc_by_id(arc1).first_node() {
            self.node_mut(node).set_linked_arc(Some(arc1));
        }

        // `cxx:364-374`: elimination of Arc2 and of the two old nodes.
        if let Some(old_node1) = old_node1 {
            let key = self.node_by_id(old_node1).index();
            if self.the_nodes.remove(&key).is_some() {
                self.number_of_nodes -= 1;
            }
        }
        if let Some(old_node2) = old_node2 {
            let key = self.node_by_id(old_node2).index();
            if self.the_nodes.remove(&key).is_some() {
                self.number_of_nodes -= 1;
            }
        }

        // `cxx:378-381`.
        let arc2_index = self.arc_by_id(arc2).index();
        self.the_arcs.remove(&arc2_index);
        self.number_of_arcs -= 1;

        // `cxx:383-402`: the basic elements no longer start or end on Arc2.
        let elements = [
            self.arc_by_id(arc2).first_element(),
            self.arc_by_id(arc2).second_element(),
        ];
        for element in elements.into_iter().flatten() {
            let be = self.basic_elt_by_id(element).index();
            let be = match self.the_basic_elts.get(&be) {
                Some(&b) => b,
                None => continue,
            };            if self.basic_elt_by_id(be).start_arc() == Some(arc2) {
                self.basic_elt_mut(be).set_start_arc(Some(arc1));
            }
            if self.basic_elt_by_id(be).end_arc() == Some(arc2) {
                self.basic_elt_mut(be).set_end_arc(Some(arc1));
            }
        }
    }

    /// `CompactArcs()` (`cxx:409-433`): shifts the arcs down to fill the holes
    /// left by `FusionOfArcs`.
    pub fn compact_arcs(&mut self) {
        let mut find = 0;
        let mut i = 1;
        let mut ya_decalage = false;

        while find < self.number_of_arcs {
            match self.the_arcs.get(&i).copied() {
                None => ya_decalage = true,
                Some(arc) => {
                    find += 1;
                    if ya_decalage {
                        self.arc_mut(arc).set_index(find);
                        self.the_arcs.remove(&i);
                        self.the_arcs.insert(find, arc);
                    }
                }
            }
            i += 1;
        }
    }

    /// `CompactNodes()` (`cxx:438-463`): the same shift for the nodes.
    pub fn compact_nodes(&mut self) {
        let mut find = 0;
        let mut i = 1;
        let mut ya_decalage = false;

        while find < self.number_of_nodes {
            match self.the_nodes.get(&i).copied() {
                None => ya_decalage = true,
                Some(node) => {
                    find += 1;
                    if ya_decalage {
                        self.node_mut(node).set_index(find);
                        self.the_nodes.remove(&i);
                        self.the_nodes.insert(find, node);
                    }
                }
            }
            i += 1;
        }
    }

    /// `ChangeBasicElts(NewMap)` (`cxx:467-474`): replaces the map and renumbers
    /// every element to its new key.
    pub fn change_basic_elts(&mut self, new_map: &BTreeMap<i32, BasicEltId>) {
        self.the_basic_elts = new_map.clone();
        for (&key, &elt) in new_map {
            self.basic_elt_mut(elt).set_index(key);
        }
    }

    /// `ChangeBasicElt(Index)` (`cxx:478-481`).
    pub fn change_basic_elt(&self, index: i32) -> BasicEltId {
        self.basic_elt(index)
    }
}
