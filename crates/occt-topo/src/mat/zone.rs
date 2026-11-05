//! `MAT_Zone` (`MAT_Zone.cxx:26-179`, `.hxx:...`): the zone of the graph that
//! can be reached from one basic element by turning always on the same side,
//! i.e. the set of arcs separating that element from the rest of the figure.
//!
//! `Perform` (`cxx:42-117`) starts on the `EndArc` of the basic element, turns
//! left until it comes back to the contour (`PendingNode`), to a node at
//! infinity, or to the start node; a zone that reaches an infinite node is
//! rebuilt from the `StartArc` turning right, and is then flagged `limited =
//! false`.

use super::elements::{ArcId, BasicEltId, NodeId, MatSide};
use super::MatGraph;

/// `MAT_Zone` (`MAT_Zone.hxx`).
#[derive(Debug, Clone)]
pub struct MatZone {
    /// `limited`.
    limited: bool,
    /// `frontier`, the arcs walked so far.
    frontier: Vec<ArcId>,
}

impl Default for MatZone {
    /// `MAT_Zone()` (`cxx:28-31`): `limited = true`, empty frontier.
    fn default() -> Self {
        Self::new()
    }
}

impl MatZone {
    /// `MAT_Zone()` (`cxx:28-31`).
    pub fn new() -> Self {
        Self { limited: true, frontier: Vec::new() }
    }

    /// `MAT_Zone(aBasicElt)` (`cxx:35-38`).
    pub fn from_basic_elt(graph: &MatGraph, a_basic_elt: BasicEltId) -> Self {
        let mut zone = Self::new();
        zone.perform(graph, a_basic_elt);
        zone
    }

    /// `Perform(aBasicElt)` (`cxx:42-117`).
    pub fn perform(&mut self, graph: &MatGraph, a_basic_elt: BasicEltId) {
        self.limited = true;
        self.frontier.clear();

        // `cxx:51-54`: no arc for that element, hence an empty zone.
        let end_arc = match graph.basic_elt_by_id(a_basic_elt).end_arc() {
            Some(a) => a,
            None => return,
        };

        // `cxx:66-68`: initialisation of the frontier.
        let mut current_arc = end_arc;
        self.frontier.push(current_arc);

        // `cxx:72-73`: the first node that lets the zone be built by turning
        // left.
        let mut next_node = Self::node_for_turn(graph, current_arc, a_basic_elt, MatSide::Left);
        let start_node = graph.arc_by_id(current_arc).the_other_node(next_node);

        // `cxx:83-87`: turn left until the contour is reached, an infinite node
        // is met, or the walking comes back to the start.
        while !graph.pending_node(next_node) && next_node != start_node {
            current_arc = graph
                .arc_by_id(current_arc)
                .neighbour(next_node, MatSide::Left)
                .unwrap_or_else(|| panic!("MAT_Zone::Perform: left neighbour is null"));
            self.frontier.push(current_arc);
            next_node = graph.arc_by_id(current_arc).the_other_node(next_node);
        }

        // `cxx:96-116`: an infinite node means the exploration restarts from the
        // StartArc turning right.
        if graph.node_by_id(next_node).infinite() {
            self.limited = false;
            current_arc = match graph.basic_elt_by_id(a_basic_elt).start_arc() {
                Some(a) => a,
                None => return,
            };
            self.frontier.push(current_arc);

            next_node = Self::node_for_turn(graph, current_arc, a_basic_elt, MatSide::Right);

            // `cxx:110-114`: this branch also ends on an infinite node.
            while !graph.node_by_id(next_node).infinite() {
                current_arc = graph
                    .arc_by_id(current_arc)
                    .neighbour(next_node, MatSide::Right)
                    .unwrap_or_else(|| panic!("MAT_Zone::Perform: right neighbour is null"));
                self.frontier.push(current_arc);
                next_node = graph.arc_by_id(current_arc).the_other_node(next_node);
            }
        }
    }

    /// `NumberOfArcs()` (`cxx:121-124`).
    pub fn number_of_arcs(&self) -> usize {
        self.frontier.len()
    }

    /// `ArcOnFrontier(Index)` (`cxx:128-131`); 1-based like OCCT.
    pub fn arc_on_frontier(&self, index: i32) -> ArcId {
        self.frontier[(index as usize) - 1]
    }

    /// The frontier arcs in walking order.
    pub fn frontier(&self) -> &[ArcId] {
        &self.frontier
    }

    /// `NoEmptyZone()` (`cxx:135-138`).
    pub fn no_empty_zone(&self) -> bool {
        !self.frontier.is_empty()
    }

    /// `Limited()` (`cxx:142-145`).
    pub fn limited(&self) -> bool {
        self.limited
    }

    /// `NodeForTurn(anArc, aBE, aSide)` (`cxx:149-179`).
    fn node_for_turn(
        graph: &MatGraph,
        an_arc: ArcId,
        a_be: BasicEltId,
        a_side: MatSide,
    ) -> NodeId {
        let arc = graph.arc_by_id(an_arc);
        let mut node_sol = arc
            .first_node()
            .unwrap_or_else(|| panic!("MAT_Zone::NodeForTurn: first node is null"));
        let mut neighbour_arc = arc.neighbour(node_sol, a_side);

        if neighbour_arc.is_none() {
            node_sol = arc
                .second_node()
                .unwrap_or_else(|| panic!("MAT_Zone::NodeForTurn: second node is null"));
            neighbour_arc = arc.neighbour(node_sol, a_side);
        }
        let neighbour_arc = match neighbour_arc {
            Some(a) => a,
            None => return node_sol,
        };
        if graph.arc_by_id(neighbour_arc).first_element() == Some(a_be)
            || graph.arc_by_id(neighbour_arc).second_element() == Some(a_be)
        {
            return node_sol;
        }
        arc.the_other_node(node_sol)
    }
}
