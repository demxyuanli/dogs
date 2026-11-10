use std::collections::{BTreeSet, HashMap, VecDeque};

use occt_core::precision::CONFUSION;

use super::super::delaun_types::{DelaunLink, DelaunPairOfIndex, DelaunTriangle, DelaunVertex, VertexState};
use super::vertex_cells::VertexCellFilter;

/// Indexed storage for the Delaunay core: vertices, links and triangles.
///
/// Mirrors `BRepMesh_DataStructureOfDelaun`. Vertex uniqueness is provided by
/// the internal cell filter; links deduplicate by their two endpoint indices;
/// triangles are plain slots. Deleted items keep their slots (so element/link
/// ids are stable during the Bowyer–Watson insertion) and are reused or
/// compacted via the deleted lists / `clear_deleted`.
pub struct DelaunDataStructure {
    pub(super) vertices: Vec<DelaunVertex>,
    pub(super) del_vertices: VecDeque<i32>,
    pub(super) vertex_cells: VertexCellFilter,
    pub(super) node_links: Vec<Vec<i32>>,

    pub(super) links: Vec<DelaunLink>,
    pub(super) link_states: Vec<VertexState>,
    pub(super) link_elements: Vec<DelaunPairOfIndex>,
    pub(super) del_links: VecDeque<i32>,
    // Like `elements_of_domain`, `LinksOfDomain()` is an
    // `IMeshData::MapOfInteger` = `TColStd_PackedMapOfInteger`
    // (`NCollection_PackedMap<int>`, `BRepMesh_DataStructureOfDelaun.hxx:228-229`),
    // traversed in ascending id order by `getEdgesByType`
    // (`BRepMesh_MeshTool.cxx:284-300`) and `MeshAlgo` consumers. A `HashSet`
    // here leaked a per-process random iteration order into the mesh result.
    pub(super) links_of_domain: BTreeSet<i32>,

    /// Node-pair -> link ids holding that pair, the port's realisation of
    /// `IndexOf` (`BRepMesh_DataStructureOfDelaun.cxx:74`): OCCT answers it with
    /// `myLinks.FindIndex(theLink)` on `IMeshData::IDMapOfLink` (`hxx:225`), an
    /// indexed map that hashes `BRepMesh_Edge::operator==` / `IsEqual`
    /// (`BRepMesh_Edge.hxx:59-66`) - equal when the node pair matches in either
    /// orientation - so the lookup is O(1) instead of a scan over `links`.
    ///
    /// `IsEqual` is false whenever either edge is `BRepMesh_Deleted`
    /// (`hxx:61-62`), i.e. OCCT's map never answers a deleted link; the scan
    /// below applies the same filter through `link_states`.
    pub(super) link_index: HashMap<(i32, i32), Vec<i32>>,

    pub(super) elements: Vec<DelaunTriangle>,
    pub(super) element_states: Vec<VertexState>,
    // `BTreeSet` (ascending id) mirrors the traversal of OCCT's
    // `ElementsOfDomain()`, an `IMeshData::MapOfInteger` =
    // `TColStd_PackedMapOfInteger` (`NCollection_PackedMap<int>`). The packed
    // map stores 64 consecutive ids per block and walks hash buckets in order;
    // with `NextPrimeForMap(0) == 101` buckets there is no collision until an
    // id block index reaches 101, i.e. element id >= 6464, so its traversal is
    // exactly ascending id. `ElementsOfDomain` order feeds the control-node
    // list in `BRepMesh_DelaunayDeflectionControlMeshAlgo.hxx:102-107`.
    pub(super) elements_of_domain: BTreeSet<i32>,
}

impl DelaunDataStructure {
    /// Sets the 2D cell-filter pitch. Source: `BRepMesh_VertexTool::SetCellSize`.
    pub fn set_cell_size(&mut self, size_u: f64, size_v: f64) {
        self.vertex_cells.cell_size = (size_u.max(1e-16), size_v.max(1e-16));
    }

    /// Sets the coincidence tolerance. Source: `BRepMesh_VertexTool::SetTolerance`.
    pub fn set_tolerance(&mut self, tol_u: f64, tol_v: f64) {
        self.vertex_cells.tolerance = (tol_u, tol_v);
    }

    /// Creates an empty data structure. `reserved_node_size` is only a capacity
    /// hint, matching the OCCT constructor.
    pub fn new(reserved_node_size: usize) -> Self {
        let cap = reserved_node_size.max(16);
        let cell_size = CONFUSION + 0.05 * CONFUSION;
        Self {
            vertices: Vec::with_capacity(cap),
            del_vertices: VecDeque::new(),
            vertex_cells: VertexCellFilter::new(cell_size, (CONFUSION, CONFUSION)),
            node_links: Vec::with_capacity(cap),
            links: Vec::with_capacity(cap * 3),
            link_states: Vec::with_capacity(cap * 3),
            link_elements: Vec::with_capacity(cap * 3),
            del_links: VecDeque::new(),
            links_of_domain: BTreeSet::new(),
            link_index: HashMap::new(),
            elements: Vec::with_capacity(cap * 2),
            element_states: Vec::with_capacity(cap * 2),
            elements_of_domain: BTreeSet::new(),
        }
    }
}
