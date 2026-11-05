//! The `MAT` graph element classes: `MAT_BasicElt` (`MAT_BasicElt.cxx:27-90`),
//! `MAT_Edge` (`MAT_Edge.cxx:22-80`), `MAT_Bisector` (`MAT_Bisector.cxx:25-206`),
//! `MAT_Node` (`MAT_Node.cxx:24-150`) and `MAT_Arc` (`MAT_Arc.cxx:27-271`).
//!
//! OCCT stores these as `occ::handle` graphs with cyclic raw pointers
//! (`MAT_Arc::firstArcLeft` etc.). The port keeps the same shape but replaces
//! every handle with a `usize` slot id: the graph owns one arena per class and
//! a handle is the arena index, so the cycles are plain ids. Node and arc
//! identities are the slot ids, never the (mutable, sometimes not yet assigned)
//! `Index()` values that OCCT's `theNodes` / `theArcs` maps are keyed by.
//!
//! Queries that need to walk the graph (`MAT_Node::LinkedArcs`,
//! `MAT_Node::NearElts`, `MAT_Node::PendingNode`) are methods of
//! [`MatGraph`](super::MatGraph) because they have to resolve neighbours.

use occt_core::precision::INFINITE;

/// `MAT_Side` (`MAT_Side.hxx`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatSide {
    /// `MAT_Left`.
    Left,
    /// `MAT_Right`.
    Right,
}

/// Slot id of a `MAT_BasicElt`.
pub type BasicEltId = usize;
/// Slot id of a `MAT_Edge`.
pub type EdgeId = usize;
/// Slot id of a `MAT_Bisector`.
pub type BisectorId = usize;
/// Slot id of a `MAT_Node`.
pub type NodeId = usize;
/// Slot id of a `MAT_Arc`.
pub type ArcId = usize;

/// `MAT_BasicElt` (`MAT_BasicElt.cxx:26-90`): one element of the input figure,
/// i.e. one curve of the contours the medial axis is built from.
#[derive(Debug, Clone, PartialEq)]
pub struct MatBasicElt {
    /// `startLeftArc`, stored as the arc slot it points at.
    start_left_arc: Option<ArcId>,
    /// `endLeftArc`.
    end_left_arc: Option<ArcId>,
    /// `index`.
    index: i32,
    /// `geomIndex`.
    geom_index: i32,
}

impl MatBasicElt {
    /// `MAT_BasicElt(anInteger)` (`cxx:27-33`): both arcs null, `geomIndex` 0.
    pub fn new(an_integer: i32) -> Self {
        Self {
            start_left_arc: None,
            end_left_arc: None,
            index: an_integer,
            geom_index: 0,
        }
    }

    /// `StartArc()` (`cxx:37-40`).
    pub fn start_arc(&self) -> Option<ArcId> {
        self.start_left_arc
    }

    /// `EndArc()` (`cxx:44-47`).
    pub fn end_arc(&self) -> Option<ArcId> {
        self.end_left_arc
    }

    /// `Index()` (`cxx:52-55`).
    pub fn index(&self) -> i32 {
        self.index
    }

    /// `GeomIndex()` (`cxx:59-62`).
    pub fn geom_index(&self) -> i32 {
        self.geom_index
    }

    /// `SetStartArc(anArc)` (`cxx:66-69`).
    pub fn set_start_arc(&mut self, an_arc: Option<ArcId>) {
        self.start_left_arc = an_arc;
    }

    /// `SetEndArc(anArc)` (`cxx:73-76`).
    pub fn set_end_arc(&mut self, an_arc: Option<ArcId>) {
        self.end_left_arc = an_arc;
    }

    /// `SetIndex(anInteger)` (`cxx:80-83`).
    pub fn set_index(&mut self, an_integer: i32) {
        self.index = an_integer;
    }

    /// `SetGeomIndex(anInteger)` (`cxx:87-90`).
    pub fn set_geom_index(&mut self, an_integer: i32) {
        self.geom_index = an_integer;
    }
}

/// `MAT_Edge` (`MAT_Edge.cxx:21-80`): an edge of the medial axis, i.e. the
/// locus swept between two bisectors (`thefirstbisector`, `thesecondbisector`).
#[derive(Debug, Clone, PartialEq)]
pub struct MatEdge {
    /// `theedgenumber`.
    the_edge_number: i32,
    /// `thefirstbisector`.
    the_first_bisector: Option<BisectorId>,
    /// `thesecondbisector`.
    the_second_bisector: Option<BisectorId>,
    /// `thedistance`.
    the_distance: f64,
    /// `theintersectionpoint`.
    the_intersection_point: i32,
}

impl Default for MatEdge {
    /// `MAT_Edge()` (`cxx:22-27`).
    fn default() -> Self {
        Self::new()
    }
}

impl MatEdge {
    /// `MAT_Edge()` (`cxx:22-27`).
    pub fn new() -> Self {
        Self {
            the_edge_number: 0,
            the_first_bisector: None,
            the_second_bisector: None,
            the_distance: 0.0,
            the_intersection_point: 0,
        }
    }

    /// `EdgeNumber(anumber)` (`cxx:30-33`).
    pub fn set_edge_number(&mut self, a_number: i32) {
        self.the_edge_number = a_number;
    }

    /// `FirstBisector(abisector)` (`cxx:35-38`).
    pub fn set_first_bisector(&mut self, a_bisector: Option<BisectorId>) {
        self.the_first_bisector = a_bisector;
    }

    /// `SecondBisector(abisector)` (`cxx:41-44`).
    pub fn set_second_bisector(&mut self, a_bisector: Option<BisectorId>) {
        self.the_second_bisector = a_bisector;
    }

    /// `Distance(adistance)` (`cxx:47-50`).
    pub fn set_distance(&mut self, a_distance: f64) {
        self.the_distance = a_distance;
    }

    /// `IntersectionPoint(apoint)` (`cxx:53-56`).
    pub fn set_intersection_point(&mut self, a_point: i32) {
        self.the_intersection_point = a_point;
    }

    /// `EdgeNumber()` (`cxx:58-61`).
    pub fn edge_number(&self) -> i32 {
        self.the_edge_number
    }

    /// `FirstBisector()` (`cxx:63-66`).
    pub fn first_bisector(&self) -> Option<BisectorId> {
        self.the_first_bisector
    }

    /// `SecondBisector()` (`cxx:68-71`).
    pub fn second_bisector(&self) -> Option<BisectorId> {
        self.the_second_bisector
    }

    /// `Distance()` (`cxx:73-76`).
    pub fn distance(&self) -> f64 {
        self.the_distance
    }

    /// `IntersectionPoint()` (`cxx:78-81`).
    pub fn intersection_point(&self) -> i32 {
        self.the_intersection_point
    }
}

/// `MAT_Bisector` (`MAT_Bisector.cxx:25-206`): one node of the tree of
/// bisectors; `thelistofbisectors` holds the sub-bisectors reached from it.
#[derive(Debug, Clone, PartialEq)]
pub struct MatBisector {
    /// `thebisectornumber`.
    the_bisector_number: i32,
    /// `theindexnumber`.
    the_index_number: i32,
    /// `thefirstedge`.
    the_first_edge: Option<EdgeId>,
    /// `thesecondedge`.
    the_second_edge: Option<EdgeId>,
    /// `theissuepoint`.
    the_issue_point: i32,
    /// `theendpoint`.
    the_end_point: i32,
    /// `distissuepoint`.
    dist_issue_point: f64,
    /// `thefirstvector`.
    the_first_vector: i32,
    /// `thesecondvector`.
    the_second_vector: i32,
    /// `thesense`.
    the_sense: f64,
    /// `thefirstparameter`.
    the_first_parameter: f64,
    /// `thesecondparameter`.
    the_second_parameter: f64,
    /// `thelistofbisectors` (`MAT_ListOfBisector`, filled with `BackAdd`).
    the_list_of_bisectors: super::MatTList<BisectorId>,
}

impl Default for MatBisector {
    /// `MAT_Bisector()` (`cxx:25-31`).
    fn default() -> Self {
        Self::new()
    }
}

impl MatBisector {
    /// `MAT_Bisector()` (`cxx:25-31`): number `-1`, both parameters
    /// `Precision::Infinite()`, empty sub-list.
    pub fn new() -> Self {
        Self {
            the_bisector_number: -1,
            the_index_number: 0,
            the_first_edge: None,
            the_second_edge: None,
            the_issue_point: 0,
            the_end_point: 0,
            dist_issue_point: 0.0,
            the_first_vector: 0,
            the_second_vector: 0,
            the_sense: 0.0,
            the_first_parameter: INFINITE,
            the_second_parameter: INFINITE,
            the_list_of_bisectors: super::MatTList::new(),
        }
    }

    /// `AddBisector(abisector)` (`cxx:33-36`): `BackAdd`.
    pub fn add_bisector(&mut self, a_bisector: BisectorId) {
        self.the_list_of_bisectors.back_add(a_bisector);
    }

    /// `List()` (`cxx:38-41`).
    pub fn list(&self) -> &super::MatTList<BisectorId> {
        &self.the_list_of_bisectors
    }

    /// `MAT_Bisector::FirstBisector()` (`cxx:43-46`): `List()->FirstItem()`.
    pub fn first_bisector(&self) -> Option<BisectorId> {
        if self.the_list_of_bisectors.is_empty() {
            None
        } else {
            Some(*self.the_list_of_bisectors.first_item())
        }
    }

    /// `MAT_Bisector::LastBisector()` (`cxx:48-51`): `List()->LastItem()`.
    pub fn last_bisector(&self) -> Option<BisectorId> {
        if self.the_list_of_bisectors.is_empty() {
            None
        } else {
            Some(*self.the_list_of_bisectors.last_item())
        }
    }

    /// `BisectorNumber(anumber)` (`cxx:53-56`).
    pub fn set_bisector_number(&mut self, a_number: i32) {
        self.the_bisector_number = a_number;
    }

    /// `IndexNumber(anumber)` (`cxx:58-61`).
    pub fn set_index_number(&mut self, a_number: i32) {
        self.the_index_number = a_number;
    }

    /// `FirstEdge(anedge)` (`cxx:63-66`).
    pub fn set_first_edge(&mut self, an_edge: Option<EdgeId>) {
        self.the_first_edge = an_edge;
    }

    /// `SecondEdge(anedge)` (`cxx:68-71`).
    pub fn set_second_edge(&mut self, an_edge: Option<EdgeId>) {
        self.the_second_edge = an_edge;
    }

    /// `IssuePoint(apoint)` (`cxx:73-76`).
    pub fn set_issue_point(&mut self, a_point: i32) {
        self.the_issue_point = a_point;
    }

    /// `EndPoint(apoint)` (`cxx:78-81`).
    pub fn set_end_point(&mut self, a_point: i32) {
        self.the_end_point = a_point;
    }

    /// `DistIssuePoint(areal)` (`cxx:83-86`).
    pub fn set_dist_issue_point(&mut self, a_real: f64) {
        self.dist_issue_point = a_real;
    }

    /// `FirstVector(avector)` (`cxx:88-91`).
    pub fn set_first_vector(&mut self, a_vector: i32) {
        self.the_first_vector = a_vector;
    }

    /// `SecondVector(avector)` (`cxx:93-96`).
    pub fn set_second_vector(&mut self, a_vector: i32) {
        self.the_second_vector = a_vector;
    }

    /// `Sense(asense)` (`cxx:98-101`).
    pub fn set_sense(&mut self, a_sense: f64) {
        self.the_sense = a_sense;
    }

    /// `FirstParameter(aparameter)` (`cxx:103-106`).
    pub fn set_first_parameter(&mut self, a_parameter: f64) {
        self.the_first_parameter = a_parameter;
    }

    /// `SecondParameter(aparameter)` (`cxx:108-111`).
    pub fn set_second_parameter(&mut self, a_parameter: f64) {
        self.the_second_parameter = a_parameter;
    }

    /// `BisectorNumber()` (`cxx:113-116`).
    pub fn bisector_number(&self) -> i32 {
        self.the_bisector_number
    }

    /// `IndexNumber()` (`cxx:118-121`).
    pub fn index_number(&self) -> i32 {
        self.the_index_number
    }

    /// `FirstEdge()` (`cxx:123-126`).
    pub fn first_edge(&self) -> Option<EdgeId> {
        self.the_first_edge
    }

    /// `SecondEdge()` (`cxx:128-131`).
    pub fn second_edge(&self) -> Option<EdgeId> {
        self.the_second_edge
    }

    /// `IssuePoint()` (`cxx:133-136`).
    pub fn issue_point(&self) -> i32 {
        self.the_issue_point
    }

    /// `EndPoint()` (`cxx:138-141`).
    pub fn end_point(&self) -> i32 {
        self.the_end_point
    }

    /// `DistIssuePoint()` (`cxx:143-146`).
    pub fn dist_issue_point(&self) -> f64 {
        self.dist_issue_point
    }

    /// `FirstVector()` (`cxx:148-151`).
    pub fn first_vector(&self) -> i32 {
        self.the_first_vector
    }

    /// `SecondVector()` (`cxx:153-156`).
    pub fn second_vector(&self) -> i32 {
        self.the_second_vector
    }

    /// `Sense()` (`cxx:158-161`).
    pub fn sense(&self) -> f64 {
        self.the_sense
    }

    /// `FirstParameter()` (`cxx:163-166`).
    pub fn first_parameter(&self) -> f64 {
        self.the_first_parameter
    }

    /// `SecondParameter()` (`cxx:168-171`).
    pub fn second_parameter(&self) -> f64 {
        self.the_second_parameter
    }
}

/// `MAT_Node` (`MAT_Node.cxx:24-150`): a vertex of the medial axis, linked to
/// one arc (`aLinkedArc`). A node is *pending* while that arc has no left
/// neighbour around it, and *infinite* at `Precision::Infinite()` distance.
#[derive(Debug, Clone, PartialEq)]
pub struct MatNode {
    /// `nodeIndex`.
    node_index: i32,
    /// `geomIndex`.
    geom_index: i32,
    /// `distance`.
    distance: f64,
    /// `aLinkedArc`.
    linked_arc: Option<ArcId>,
}

impl MatNode {
    /// `MAT_Node(GeomIndex, LinkedArc, Distance)` (`cxx:25-32`).
    pub fn new(geom_index: i32, linked_arc: Option<ArcId>, distance: f64) -> Self {
        Self {
            node_index: 0,
            geom_index,
            distance,
            linked_arc,
        }
    }

    /// `GeomIndex()` (`cxx:36-39`).
    pub fn geom_index(&self) -> i32 {
        self.geom_index
    }

    /// `Index()` (`cxx:43-46`).
    pub fn index(&self) -> i32 {
        self.node_index
    }

    /// `Distance()` (`cxx:112-115`).
    pub fn distance(&self) -> f64 {
        self.distance
    }

    /// `OnBasicElt()` (`cxx:126-129`).
    pub fn on_basic_elt(&self) -> bool {
        self.distance() == 0.0
    }

    /// `Infinite()` (`cxx:133-136`).
    pub fn infinite(&self) -> bool {
        self.distance() == INFINITE
    }

    /// `SetLinkedArc(LinkedArc)` (`cxx:140-143`).
    pub fn set_linked_arc(&mut self, linked_arc: Option<ArcId>) {
        self.linked_arc = linked_arc;
    }

    /// `SetIndex(anIndex)` (`cxx:147-150`).
    pub fn set_index(&mut self, an_index: i32) {
        self.node_index = an_index;
    }

    /// `aLinkedArc` (`MAT_Node.hxx`).
    pub fn linked_arc(&self) -> Option<ArcId> {
        self.linked_arc
    }
}

/// `MAT_Arc` (`MAT_Arc.cxx:26-271`): an arc of the medial axis, separating two
/// basic elements between two nodes.
#[derive(Debug, Clone, PartialEq)]
pub struct MatArc {
    /// `arcIndex`.
    arc_index: i32,
    /// `geomIndex`.
    geom_index: i32,
    /// `firstArcLeft`.
    first_arc_left: Option<ArcId>,
    /// `firstArcRight`.
    first_arc_right: Option<ArcId>,
    /// `secondArcRight`.
    second_arc_right: Option<ArcId>,
    /// `secondArcLeft`.
    second_arc_left: Option<ArcId>,
    /// `firstElement`.
    first_element: Option<BasicEltId>,
    /// `secondElement`.
    second_element: Option<BasicEltId>,
    /// `firstNode`.
    first_node: Option<NodeId>,
    /// `secondNode`.
    second_node: Option<NodeId>,
}

impl MatArc {
    /// `MAT_Arc(ArcIndex, GeomIndex, FirstElement, SecondElement)`
    /// (`cxx:27-39`): all four arc neighbours start null.
    pub fn new(
        arc_index: i32,
        geom_index: i32,
        first_element: Option<BasicEltId>,
        second_element: Option<BasicEltId>,
    ) -> Self {
        Self {
            arc_index,
            geom_index,
            first_arc_left: None,
            first_arc_right: None,
            second_arc_right: None,
            second_arc_left: None,
            first_element,
            second_element,
            first_node: None,
            second_node: None,
        }
    }

    /// `Index()` (`cxx:43-46`).
    pub fn index(&self) -> i32 {
        self.arc_index
    }

    /// `GeomIndex()` (`cxx:50-53`).
    pub fn geom_index(&self) -> i32 {
        self.geom_index
    }

    /// `FirstElement()` (`cxx:57-60`).
    pub fn first_element(&self) -> Option<BasicEltId> {
        self.first_element
    }

    /// `SecondElement()` (`cxx:63-66`).
    pub fn second_element(&self) -> Option<BasicEltId> {
        self.second_element
    }

    /// `FirstNode()` (`cxx:70-73`).
    pub fn first_node(&self) -> Option<NodeId> {
        self.first_node
    }

    /// `SecondNode()` (`cxx:77-80`).
    pub fn second_node(&self) -> Option<NodeId> {
        self.second_node
    }

    /// `TheOtherNode(aNode)` (`cxx:84-99`); panics where OCCT raises
    /// `Standard_DomainError`.
    pub fn the_other_node(&self, a_node: NodeId) -> NodeId {
        if self.first_node() == Some(a_node) {
            return self.second_node().expect("MAT_Arc::TheOtherNode");
        } else if self.second_node() == Some(a_node) {
            return self.first_node().expect("MAT_Arc::TheOtherNode");
        }
        panic!("MAT_Arc::TheOtherNode: Standard_DomainError");
    }

    /// `HasNeighbour(aNode, aSide)` (`cxx:103-131`); panics where OCCT raises
    /// `Standard_DomainError`.
    pub fn has_neighbour(&self, a_node: NodeId, a_side: MatSide) -> bool {
        match a_side {
            MatSide::Left => {
                if Some(a_node) == self.first_node() {
                    return self.first_arc_left.is_some();
                }
                if Some(a_node) == self.second_node() {
                    return self.second_arc_left.is_some();
                }
            }
            MatSide::Right => {
                if Some(a_node) == self.first_node() {
                    return self.first_arc_right.is_some();
                }
                if Some(a_node) == self.second_node() {
                    return self.second_arc_right.is_some();
                }
            }
        }
        panic!("MAT_Arc::HasNeighbour: Standard_DomainError");
    }

    /// `Neighbour(aNode, aSide)` (`cxx:135-163`); panics where OCCT raises
    /// `Standard_DomainError`.
    pub fn neighbour(&self, a_node: NodeId, a_side: MatSide) -> Option<ArcId> {
        match a_side {
            MatSide::Left => {
                if Some(a_node) == self.first_node() {
                    return self.first_arc_left;
                }
                if Some(a_node) == self.second_node() {
                    return self.second_arc_left;
                }
            }
            MatSide::Right => {
                if Some(a_node) == self.first_node() {
                    return self.first_arc_right;
                }
                if Some(a_node) == self.second_node() {
                    return self.second_arc_right;
                }
            }
        }
        panic!("MAT_Arc::Neighbour: Standard_DomainError");
    }

    /// `SetIndex(anInteger)` (`cxx:167-170`).
    pub fn set_index(&mut self, an_integer: i32) {
        self.arc_index = an_integer;
    }

    /// `SetGeomIndex(anInteger)` (`cxx:172-175`).
    pub fn set_geom_index(&mut self, an_integer: i32) {
        self.geom_index = an_integer;
    }

    /// `SetFirstElement(aBasicElt)` (`cxx:179-182`).
    pub fn set_first_element(&mut self, a_basic_elt: Option<BasicEltId>) {
        self.first_element = a_basic_elt;
    }

    /// `SetSecondElement(aBasicElt)` (`cxx:186-189`).
    pub fn set_second_element(&mut self, a_basic_elt: Option<BasicEltId>) {
        self.second_element = a_basic_elt;
    }

    /// `SetFirstNode(aNode)` (`cxx:193-196`).
    pub fn set_first_node(&mut self, a_node: Option<NodeId>) {
        self.first_node = a_node;
    }

    /// `SetSecondNode(aNode)` (`cxx:198-201`).
    pub fn set_second_node(&mut self, a_node: Option<NodeId>) {
        self.second_node = a_node;
    }

    /// `SetFirstArc(aSide, anArc)` (`cxx:205-219`).
    pub fn set_first_arc(&mut self, a_side: MatSide, an_arc: Option<ArcId>) {
        if a_side == MatSide::Left {
            self.first_arc_left = an_arc;
        } else {
            self.first_arc_right = an_arc;
        }
    }

    /// `SetSecondArc(aSide, anArc)` (`cxx:221-231`).
    pub fn set_second_arc(&mut self, a_side: MatSide, an_arc: Option<ArcId>) {
        if a_side == MatSide::Left {
            self.second_arc_left = an_arc;
        } else {
            self.second_arc_right = an_arc;
        }
    }

    /// `SetNeighbour(aSide, aNode, anArc)` (`cxx:235-271`); panics where OCCT
    /// raises `Standard_DomainError`.
    pub fn set_neighbour(&mut self, a_side: MatSide, a_node: NodeId, an_arc: Option<ArcId>) {
        match a_side {
            MatSide::Left => {
                if Some(a_node) == self.first_node() {
                    self.first_arc_left = an_arc;
                } else if Some(a_node) == self.second_node() {
                    self.second_arc_left = an_arc;
                } else {
                    panic!("MAT_Arc::SetNeighbour: Standard_DomainError");
                }
            }
            MatSide::Right => {
                if Some(a_node) == self.first_node() {
                    self.first_arc_right = an_arc;
                } else if Some(a_node) == self.second_node() {
                    self.second_arc_right = an_arc;
                } else {
                    panic!("MAT_Arc::SetNeighbour: Standard_DomainError");
                }
            }
        }
    }
}
