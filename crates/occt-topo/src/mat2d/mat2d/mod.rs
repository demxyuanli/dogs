//! `MAT2d_Mat2d` (`MAT2d_Mat2d.cxx:17-1911`, `.hxx:38-104`): the bisecting locus
//! algorithm proper.
//!
//! The `MAT` graph is the *output* of this class, but OCCT shares the
//! `MAT_Bisector` / `MAT_Edge` handles between `MAT2d_Mat2d` and `MAT_Graph`
//! (`BRepMAT2d_BisectingLocus.cxx:61-137` builds the graph from the roots this
//! class leaves behind). The port therefore owns the two arenas through
//! [`MatGraph`], which already exposes the slot-id handles `MAT2d_Mat2d` needs,
//! and hands the same graph to `MatGraph::perform` afterwards.
//!
//! | OCCT | lines | file |
//! |---|---|---|
//! | `CreateMatOpen` | `cxx:122-859` | `create_mat.rs` |
//! | `CreateMat` | `cxx:852-1628` | `create_mat.rs` |
//! | `LoadBisectorsToRemove` | `cxx:1632-1699` | `remove.rs` |
//! | `Intersect` | `cxx:1701-1826` | `remove.rs` |
//! | `Init`/`More`/`Next`/`Bisector`/`NumberOfBisectors`/`SemiInfinite`/`IsDone` | `cxx:1828-1876` | this file |

mod create_mat;
mod finish;
mod remove;

use std::collections::BTreeMap;

use occt_core::precision::INFINITE;

use crate::mat::{BisectorId, EdgeId, MatBisector, MatEdge, MatGraph, MatTList};
use crate::mat2d::Mat2dTool2d;

/// `class MAT2d_Mat2d` (`MAT2d_Mat2d.hxx:38-104`).
pub struct Mat2dMat2d {
    /// `myIsOpenResult` (`hxx:92`).
    my_is_open_result: bool,
    /// `thenumberofbisectors` (`hxx:93`).
    number_of_bisectors: i32,
    /// `thenumberofedges` (`hxx:94`): set to 0 by the constructor
    /// (`cxx:35-37`) and never read again in OCCT.
    #[allow(dead_code)]
    number_of_edges: i32,
    /// `semiInfinite` (`hxx:95`).
    semi_infinite: bool,
    /// `theedgelist` (`hxx:96`), a `MAT_ListOfEdge`.
    edge_list: MatTList<EdgeId>,
    /// `RemovedEdgesList` (`hxx:97`): only ever appended to, and cleared in the
    /// destructor (`cxx:1899-1908`), which the port has no use for.
    removed_edges: Vec<EdgeId>,
    /// `typeofbisectortoremove` (`hxx:98`).
    type_of_bisector_to_remove: BTreeMap<usize, i32>,
    /// `bisectoronetoremove` (`hxx:99`).
    bisector_one_to_remove: BTreeMap<usize, BisectorId>,
    /// `bisectortwotoremove` (`hxx:100`).
    bisector_two_to_remove: BTreeMap<usize, BisectorId>,
    /// `bisectormap` (`hxx:101`), a `NCollection_DataMap<int, handle>` keyed by
    /// `IndexNumber`.
    bisector_map: BTreeMap<i32, BisectorId>,
    /// `roots` (`hxx:102`), a `MAT_ListOfBisector`.
    roots: MatTList<BisectorId>,
    /// `isDone` (`hxx:103`).
    is_done: bool,
    /// The arena pair the `MAT_Bisector` / `MAT_Edge` handles live in.
    graph: MatGraph,
}

impl Default for Mat2dMat2d {
    fn default() -> Self {
        Self::new(false)
    }
}

impl Mat2dMat2d {
    /// `MAT2d_Mat2d(IsOpenResult)` (`cxx:30-38`).
    pub fn new(is_open_result: bool) -> Self {
        Self {
            my_is_open_result: is_open_result,
            number_of_bisectors: 0,
            number_of_edges: 0,
            semi_infinite: false,
            edge_list: MatTList::new(),
            removed_edges: Vec::new(),
            type_of_bisector_to_remove: BTreeMap::new(),
            bisector_one_to_remove: BTreeMap::new(),
            bisector_two_to_remove: BTreeMap::new(),
            bisector_map: BTreeMap::new(),
            roots: MatTList::new(),
            is_done: false,
            graph: MatGraph::new(),
        }
    }

    /// `CreateMat(aTool)` (`cxx:852-1628`): the closed-contour variant.
    pub fn create_mat(&mut self, tool: &mut Mat2dTool2d) {
        self.create(tool, true);
    }

    /// `CreateMatOpen(aTool)` (`cxx:122-859`): the open-wire variant.
    pub fn create_mat_open(&mut self, tool: &mut Mat2dTool2d) {
        self.create(tool, false);
    }

    /// `IsDone()` (`cxx:1872-1875`).
    pub fn is_done(&self) -> bool {
        self.is_done
    }

    /// `Init()` (`cxx:1828-1831`).
    pub fn init(&mut self) {
        self.roots.first();
    }

    /// `More()` (`cxx:1836-1839`).
    pub fn more(&self) -> bool {
        self.roots.more()
    }

    /// `Next()` (`cxx:1844-1847`).
    pub fn next(&mut self) {
        self.roots.next();
    }

    /// `Bisector()` (`cxx:1852-1855`).
    pub fn bisector(&self) -> BisectorId {
        *self.roots.current()
    }

    /// `SemiInfinite()` (`cxx:1863-1866`).
    pub fn semi_infinite(&self) -> bool {
        self.semi_infinite
    }

    /// `NumberOfBisectors()` (`cxx:1858-1861`).
    pub fn number_of_bisectors(&self) -> i32 {
        self.number_of_bisectors
    }

    /// The arena the roots belong to: the caller feeds it to
    /// `MatGraph::perform` (`BRepMAT2d_BisectingLocus.cxx:134-137`).
    pub fn graph(&self) -> &MatGraph {
        &self.graph
    }

    /// A mutable view of the arena, for the `perform` call above.
    pub fn graph_mut(&mut self) -> &mut MatGraph {
        &mut self.graph
    }

    /// `theedgelist` in list order; the destructor of `MAT2d_Mat2d` walks it.
    pub fn edge_list(&self) -> &MatTList<EdgeId> {
        &self.edge_list
    }

    // ---------------------------------------------------------------
    // Shared helpers
    // ---------------------------------------------------------------

    /// `bisectormap(anindex)` (`NCollection_DataMap::operator()`).
    pub(crate) fn map(&self, an_index: i32) -> BisectorId {
        *self
            .bisector_map
            .get(&an_index)
            .unwrap_or_else(|| panic!("MAT2d_Mat2d: bisectormap has no key {an_index}"))
    }

    /// A bisector by handle.
    pub(crate) fn bis(&self, id: BisectorId) -> &MatBisector {
        self.graph.bisector(id)
    }

    /// A mutable bisector by handle.
    pub(crate) fn bis_mut(&mut self, id: BisectorId) -> &mut MatBisector {
        self.graph.bisector_mut(id)
    }

    /// An edge by handle.
    pub(crate) fn edg(&self, id: EdgeId) -> &MatEdge {
        self.graph.edge(id)
    }

    /// A mutable edge by handle.
    pub(crate) fn edg_mut(&mut self, id: EdgeId) -> &mut MatEdge {
        self.graph.edge_mut(id)
    }

    /// `aBisector->FirstEdge()->EdgeNumber()`.
    pub(crate) fn first_edge_number(&self, id: BisectorId) -> i32 {
        let edge = self
            .bis(id)
            .first_edge()
            .unwrap_or_else(|| panic!("MAT2d_Mat2d: bisector {id} has no FirstEdge"));
        self.edg(edge).edge_number()
    }

    /// `aBisector->SecondEdge()->EdgeNumber()`.
    pub(crate) fn second_edge_number(&self, id: BisectorId) -> i32 {
        let edge = self
            .bis(id)
            .second_edge()
            .unwrap_or_else(|| panic!("MAT2d_Mat2d: bisector {id} has no SecondEdge"));
        self.edg(edge).edge_number()
    }

    /// The `[FirstEdge, SecondEdge]` numbers `MAT2d_Tool2d::IntersectBisector`
    /// reads from a bisector (`MAT2d_Tool2d.cxx:816-819`).
    pub(crate) fn edge_numbers(&self, id: BisectorId) -> [i32; 2] {
        [self.first_edge_number(id), self.second_edge_number(id)]
    }

    /// `atool.IntersectBisector(A, B, IntPnt)` (`MAT2d_Tool2d.cxx:787-1101`)
    /// together with the parameter restores of `MAT2d_Mat2d::Intersect`
    /// (`cxx:1794-1807`).
    ///
    /// `IntersectBisector` overwrites exactly two fields - `A.SecondParameter`
    /// and `B.FirstParameter` (`MAT2d_Tool2d.cxx:1089-1090`); OCCT then restores
    /// one of them. The port clones the two bisectors to get past the two
    /// mutable borrows and writes the two fields back.
    pub(crate) fn intersect_pair(
        &mut self,
        tool: &mut Mat2dTool2d,
        a: BisectorId,
        b: BisectorId,
        restore_a_second: bool,
        restore_b_first: bool,
    ) -> (f64, Option<i32>) {
        let edges_a = self.edge_numbers(a);
        let edges_b = self.edge_numbers(b);
        let mut bisector_a = self.bis(a).clone();
        let mut bisector_b = self.bis(b).clone();
        let saved_a = bisector_a.second_parameter();
        let saved_b = bisector_b.first_parameter();
        let result = tool.intersect_bisector(
            &mut bisector_a,
            edges_a,
            &mut bisector_b,
            edges_b,
        );
        if restore_a_second {
            bisector_a.set_second_parameter(saved_a);
        }
        if restore_b_first {
            bisector_b.set_first_parameter(saved_b);
        }
        *self.bis_mut(a) = bisector_a;
        *self.bis_mut(b) = bisector_b;
        result
    }
}

/// One step of the `parama` / `paramb` decision tree of step 2.5
/// (`MAT2d_Mat2d.cxx:474-556` closed, `:1244-1329` open).
///
/// The two copies of the tree are identical; the only caller-side work is the
/// `Distance` query, which OCCT performs inside the fourth arm.
pub(crate) fn classify(bisector: &MatBisector, tool: &Mat2dTool2d, tol: f64) -> (i32, i32) {
    let first = bisector.first_parameter();
    let second = bisector.second_parameter();
    if first == INFINITE && second == INFINITE {
        (-1, -1)
    } else if first == INFINITE {
        (-1, 1)
    } else if second == INFINITE {
        (1, -1)
    } else if tool.distance(bisector, first, second) > tol {
        if (first - second) * bisector.sense() > 0.0 {
            (-1, 1)
        } else {
            (1, -1)
        }
    } else {
        (1, 1)
    }
}
