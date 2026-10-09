use super::prelude::*;
use super::*;

impl Delaun {

    /// Creates a triangulation from a fresh data structure and the given
    /// vertices (OCCT `BRepMesh_Delaun(Array1OfVertexOfDelaun&)`).
    pub fn new_vertices(vertices: &[DelaunVertex]) -> Self {
        let mut data = DelaunDataStructure::new(vertices.len().max(16));
        let mut indices: Vec<i32> = Vec::with_capacity(vertices.len());
        for v in vertices {
            indices.push(data.add_node(*v));
        }
        let mut delaun = Self {
            mesh_data: data,
            circles: CircleTool::new(),
            sup_vert: Vec::new(),
            init_circles: false,
            failed: false,
            sup_trian: DelaunTriangle::default(),
            diag_tag: -1,
        };
        delaun.perform(&mut indices, -1, -1);
        delaun
    }

    /// Creates a triangulation over an existing data structure (OCCT
    /// `BRepMesh_Delaun(oldMesh, VectorOfInteger&)`). The vertex indices must
    /// already be registered in `data_structure`.
    pub fn new_with_data(data_structure: DelaunDataStructure, vertex_indices: &mut Vec<i32>) -> Self {
        Self::new_with_data_cells(data_structure, vertex_indices, -1, -1)
    }

    /// `BRepMesh_Delaun(oldMesh, vertices, cellsU, cellsV)`.
    pub fn new_with_data_cells(
        data_structure: DelaunDataStructure,
        vertex_indices: &mut Vec<i32>,
        cells_u: i32,
        cells_v: i32,
    ) -> Self {
        let mut delaun = Self {
            mesh_data: data_structure,
            circles: CircleTool::new(),
            sup_vert: Vec::new(),
            init_circles: false,
            failed: false,
            sup_trian: DelaunTriangle::default(),
            diag_tag: -1,
        };
        delaun.perform(vertex_indices, cells_u, cells_v);
        delaun
    }

    /// TEMPORARY diagnostic (DELSTAGE): tags this triangulator so `diag_stage`
    /// prints the structure sizes when `DELSTAGE` equals this tag.
    pub fn set_diag_tag(&mut self, tag: i32) {
        self.diag_tag = tag;
    }

    /// `BRepMesh_Delaun(oldMesh, vertices, cellsU, cellsV)` with a DELSTAGE tag
    /// applied before the pipeline runs (TEMPORARY diagnostic).
    pub fn new_with_data_cells_diag(
        data_structure: DelaunDataStructure,
        vertex_indices: &mut Vec<i32>,
        cells_u: i32,
        cells_v: i32,
        tag: i32,
    ) -> Self {
        let mut delaun = Self {
            mesh_data: data_structure,
            circles: CircleTool::new(),
            sup_vert: Vec::new(),
            init_circles: false,
            failed: false,
            sup_trian: DelaunTriangle::default(),
            diag_tag: tag,
        };
        delaun.perform(vertex_indices, cells_u, cells_v);
        delaun
    }

    /// TEMPORARY diagnostic (DELSTAGE).
    pub(super) fn diag_on(&self) -> bool {
        self.diag_tag >= 0
            && std::env::var("DELSTAGE").ok().and_then(|v| v.parse::<i32>().ok())
                == Some(self.diag_tag)
    }

    /// TEMPORARY diagnostic (DELSTAGE).
    pub(super) fn diag_stage(&self, stage: &str) {
        if self.diag_tag < 0 {
            return;
        }
        match std::env::var("DELSTAGE") {
            Ok(v) if v.parse::<i32>() == Ok(self.diag_tag) => {}
            _ => return,
        }
        eprintln!(
            "DELSTAGE tag={} {} tris={} links={} nodes={} frontier={} fixed={} free={}",
            self.diag_tag,
            stage,
            self.mesh_data.elements_of_domain().len(),
            self.mesh_data.nb_links(),
            self.mesh_data.nb_nodes(),
            self.frontier().len(),
            self.internal_edges().len(),
            self.free_edges().len(),
        );
        // TEMPORARY diagnostic (DELTRIS): dumps the live triangles of the tagged
        // face, one line per triangle, so each stage of `processConstraints` can
        // be classified against the face's boundary polygon in UV space.
        if std::env::var("DELTRIS").is_err() {
            return;
        }
        let ids: Vec<i32> = self.mesh_data.elements_of_domain().iter().copied().collect();
        for id in ids {
            let element = self.mesh_data.get_element(id);
            let n = self.mesh_data.element_nodes(&element);
            let p = |k: i32| self.mesh_data.get_node(k).location;
            let (a, b, c) = (p(n[0]), p(n[1]), p(n[2]));
            eprintln!(
                "DELTRI {stage} {id} {:.9} {:.9} {:.9} {:.9} {:.9} {:.9}",
                a.x(),
                a.y(),
                b.x(),
                b.y(),
                c.x(),
                c.y()
            );
        }
        // TEMPORARY diagnostic (DELTRIS): link/adjacency dump at the same stages.
        let mut link_ids: Vec<i32> = self.mesh_data.links_of_domain().iter().copied().collect();
        link_ids.sort_unstable();
        for lid in link_ids {
            let link = self.mesh_data.get_link(lid);
            let pair = self.mesh_data.elements_connected_to(lid);
            let state = self.mesh_data.link_movability(lid);
            let na = self.mesh_data.get_node(link.first_node()).location;
            let nb = self.mesh_data.get_node(link.last_node()).location;
            let mut conn = String::new();
            for it in 1..=pair.extent() {
                let tid = pair.index(it);
                let el = self.mesh_data.get_element(tid);
                let mut orient = ' ';
                for k in 0..3 {
                    if el.link_at(k).abs() == lid {
                        orient = if el.link_at(k) > 0 { '+' } else { '-' };
                    }
                }
                conn.push_str(&format!(" {tid}{orient}"));
            }
            eprintln!(
                "DELLINK {stage} {lid} {state} {:.9} {:.9} {:.9} {:.9}{conn}",
                na.x(),
                na.y(),
                nb.x(),
                nb.y(),
                state = state.to_str()
            );
        }
    }

    /// `BRepMesh_MeshTool::EraseFreeLinks` (`BRepMesh_MeshTool.cxx:204-219`).
    pub fn erase_free_links(&mut self) {
        let n = self.mesh_data.nb_links() as i32;
        for i in 1..=n {
            if !self.mesh_data.elements_connected_to(i).is_empty() {
                continue;
            }
            if self.mesh_data.link_movability(i) == VertexState::Deleted {
                continue;
            }
            self.mesh_data.set_link_movability(i, VertexState::Free);
            self.mesh_data.remove_link(i, false);
        }
    }

    /// Initializes the triangulation with an array of vertices (OCCT `Init`).
    pub fn init(&mut self, vertices: &[DelaunVertex]) {
        let mut indices: Vec<i32> = Vec::with_capacity(vertices.len());
        for v in vertices {
            indices.push(self.mesh_data.add_node(*v));
        }
        self.perform(&mut indices, -1, -1);
    }

    /// Forced re-initialization of the circle cell filter (OCCT
    /// `InitCirclesTool`). Binds a circumcircle for every live triangle.
    pub fn init_circles_tool_public(&mut self, cells_u: i32, cells_v: i32) {
        let mut box2 = BndB2::void();
        for i in 1..=self.mesh_data.nb_nodes() as i32 {
            box2.add_pnt(self.mesh_data.get_node(i).location);
        }
        box2.enlarge(PREC);
        self.init_circles_tool(&box2, cells_u, cells_v);
        let element_ids: Vec<i32> = self.mesh_data.elements_of_domain().iter().copied().collect();
        for id in element_ids {
            let triangle = self.mesh_data.get_element(id);
            let nodes = self.mesh_data.element_nodes(&triangle);
            self.circles.bind_circle(
                id,
                self.mesh_data.get_node(nodes[0]).location.coord,
                self.mesh_data.get_node(nodes[1]).location.coord,
                self.mesh_data.get_node(nodes[2]).location.coord,
            );
        }
    }

    /// Gives the mesh data structure.
    pub fn result(&self) -> &DelaunDataStructure {
        &self.mesh_data
    }

    /// Whether `addTriangle` hit the `Standard_OutOfRange` condition of
    /// `BRepMesh_PairOfIndex::Append` (`BRepMesh_PairOfIndex.hxx:41`). OCCT lets
    /// that exception leave `BRepMesh_Delaun` and `BRepMesh_BaseMeshAlgo::Perform`
    /// swallows it (`BRepMesh_BaseMeshAlgo.cxx:59-62`), so
    /// `commitSurfaceTriangulation` never runs and the face keeps
    /// `IMeshData_Failure` with **no** triangulation.
    pub fn failed(&self) -> bool {
        self.failed
    }

    /// Consumes the triangulator and returns the mesh data structure.
    pub fn into_result(self) -> DelaunDataStructure {
        self.mesh_data
    }

    /// Registers a vertex on the shared structure (`BRepMesh_DataStructureOfDelaun::AddNode`).
    /// Used by `insertNodes` after the base mesh (`DelaunayNodeInsertionMeshAlgo.hxx:120`).
    pub fn add_node(&mut self, node: DelaunVertex) -> i32 {
        self.mesh_data.add_node(node)
    }

    /// The circle cell-filter tool used for point location.
    pub fn circles(&self) -> &CircleTool {
        &self.circles
    }

    /// Gives the list of frontier edges. Source: `Frontier()`.
    pub fn frontier(&self) -> BTreeSet<i32> {
        self.get_edges_by_type(VertexState::Frontier)
    }

    /// Gives the list of internal (fixed) edges. Source: `InternalEdges()`.
    pub fn internal_edges(&self) -> BTreeSet<i32> {
        self.get_edges_by_type(VertexState::Fixed)
    }

    /// Gives the list of free edges used at most once. Source: `FreeEdges()`.
    pub fn free_edges(&self) -> BTreeSet<i32> {
        self.get_edges_by_type(VertexState::Free)
    }

    /// Vertex by index.
    pub fn get_vertex(&self, index: i32) -> DelaunVertex {
        *self.mesh_data.get_node(index)
    }

    /// Edge (link) by index.
    pub fn get_edge(&self, index: i32) -> DelaunLink {
        self.mesh_data.get_link(index)
    }

    /// Triangle by index.
    pub fn get_triangle(&self, index: i32) -> DelaunTriangle {
        self.mesh_data.get_element(index)
    }

    /// Removes a vertex and re-triangulates the cavity. Source: `RemoveVertex`.
    pub fn remove_vertex(&mut self, vertex: &DelaunVertex) {
        let mut selector = DelaunSelector::new(&self.mesh_data);
        selector.neighbours_of(vertex);
        let mut loop_edges: BTreeMap<i32, bool> = BTreeMap::new();
        let elements: Vec<i32> = selector.elements().iter().copied().collect();
        for id in elements {
            self.delete_triangle(id, &mut loop_edges);
        }
        self.mesh_polygon_of_cavity(&mut loop_edges);
    }

    /// Adds some vertices into the triangulation. Source: `AddVertices`.
    pub fn add_vertices(&mut self, vertices: &mut Vec<i32>) {
        vertices.sort_by(|&a, &b| {
            let va = self.mesh_data.get_node(a).location;
            let vb = self.mesh_data.get_node(b).location;
            (va.x() + va.y()).total_cmp(&(vb.x() + vb.y()))
        });
        self.create_triangles_on_new_vertices(vertices);
    }

    /// Modifies the mesh to use the given edge. The OCCT implementation is an
    /// empty stub (constraint insertion happens through `ProcessConstraints`
    /// when links are pre-marked `Frontier`/`Fixed`); kept as a faithful stub.
    pub fn use_edge(&mut self, _index: i32) -> bool {
        false
    }

    /// Forces insertion of constraint edges and frontier adjustment.
    /// Source: `ProcessConstraints`.
    pub fn process_constraints(&mut self) {
        self.diag_stage("before_insert_internal_edges");
        self.insert_internal_edges();
        self.diag_stage("after_insert_internal_edges");
        self.frontier_adjust();
        self.diag_stage("after_frontier_adjust");
    }

    /// `BRepMesh_DataStructureOfDelaun::ElementNodes` (`cxx:259-286`).
    fn element_nodes_cxx(&self, element: &DelaunTriangle) -> [i32; 3] {
        self.mesh_data.element_nodes(element)
    }

    /// Tests whether the triangle contains the vertex (with square tolerance
    /// for edge closeness). Source: `Contains` (`BRepMesh_Delaun.cxx:2578-2637`).
    /// Vertices come from `ElementNodes` (`DataStructureOfDelaun.cxx:259-286`):
    /// edges 0 and 2, not the cached `vertex_indices`.
    pub fn contains(&self, triangle_id: i32, vertex: &DelaunVertex, sq_tol: f64, edge_on: &mut i32) -> bool {
        *edge_on = 0;
        let element = self.mesh_data.get_element(triangle_id);
        let p = self.element_nodes_cxx(&element);
        let points = [
            self.mesh_data.get_node(p[0]).location,
            self.mesh_data.get_node(p[1]).location,
            self.mesh_data.get_node(p[2]).location,
        ];
        let v_edges = [
            points[1].coord.subtracted(&points[0].coord),
            points[2].coord.subtracted(&points[1].coord),
            points[0].coord.subtracted(&points[2].coord),
        ];
        let mut distance = [0.0; 3];
        let mut sq_modulus = [0.0; 3];
        let mut edge_on_id = 0usize;
        let sq_min_dist = self.calculate_dist(&v_edges, &points, vertex, &mut distance, &mut sq_modulus, &mut edge_on_id);
        if sq_min_dist < 0.0 {
            return false;
        }
        let edge_id_on = element.link_at(edge_on_id).abs();
        let is_not_free = self.mesh_data.link_movability(edge_id_on) != VertexState::Free;
        if sq_min_dist > sq_tol {
            if is_not_free && distance[edge_on_id] < sq_modulus[edge_on_id] / 5.0 {
                *edge_on = edge_id_on;
            }
        } else if is_not_free {
            return false;
        } else {
            *edge_on = edge_id_on;
        }
        distance[0] >= 0.0 && distance[1] >= 0.0 && distance[2] >= 0.0
    }

    // ------------------------------------------------------------------
    // Internal algorithm
    // ------------------------------------------------------------------

    pub(super) fn perform(&mut self, vertex_indices: &mut Vec<i32>, cells_u: i32, cells_v: i32) {
        if vertex_indices.len() <= 2 {
            return;
        }
        let mut box2 = BndB2::void();
        for &idx in vertex_indices.iter() {
            box2.add_pnt(self.mesh_data.get_node(idx).location);
        }
        box2.enlarge(PREC);
        self.init_circles_tool(&box2, cells_u, cells_v);
        self.super_mesh(&box2);
        vertex_indices.sort_by(|&a, &b| {
            let va = self.mesh_data.get_node(a).location;
            let vb = self.mesh_data.get_node(b).location;
            (va.x() + va.y()).total_cmp(&(vb.x() + vb.y()))
        });
        self.compute(vertex_indices);
    }

    pub(super) fn init_circles_tool(&mut self, box2: &BndB2, cells_u: i32, cells_v: i32) {
        let (min_x, min_y, max_x, max_y) = box2.get();
        let delta_x = max_x - min_x;
        let delta_y = max_y - min_y;
        let nb = self.mesh_data.nb_nodes();
        let scaler = if nb > 100 { 5 } else { 2 };
        self.circles.set_min_max_size(GpXY::new(min_x, min_y), GpXY::new(max_x, max_y));
        self.circles.set_cell_size(delta_x / cells_u.max(scaler) as f64, delta_y / cells_v.max(scaler) as f64);
        self.init_circles = true;
    }

    pub(super) fn super_mesh(&mut self, box2: &BndB2) {
        let (min_x, min_y, max_x, max_y) = box2.get();
        let delta_x = max_x - min_x;
        let delta_y = max_y - min_y;
        let delta_min = delta_x.min(delta_y);
        let delta_max = delta_x.max(delta_y);
        let delta = delta_x + delta_y;

        self.sup_vert.clear();
        self.sup_vert.push(self.mesh_data.add_node(DelaunVertex::new_parametric(
            (min_x + max_x) / 2.0,
            max_y + delta_max,
            VertexState::Free,
        )));
        self.sup_vert.push(self.mesh_data.add_node(DelaunVertex::new_parametric(
            min_x - delta,
            min_y - delta_min,
            VertexState::Free,
        )));
        self.sup_vert.push(self.mesh_data.add_node(DelaunVertex::new_parametric(
            max_x + delta,
            min_y - delta_min,
            VertexState::Free,
        )));

        let mut e = [0i32; 3];
        let mut o = [false; 3];
        for node_id in 0..3 {
            let first = self.sup_vert[node_id];
            let last = self.sup_vert[(node_id + 1) % 3];
            let link_idx = self.mesh_data.add_link(first, last, VertexState::Free);
            e[node_id] = link_idx.abs();
            o[node_id] = link_idx > 0;
        }
        self.sup_trian = DelaunTriangle::new(
            [
                if o[0] { e[0] } else { -e[0] },
                if o[1] { e[1] } else { -e[1] },
                if o[2] { e[2] } else { -e[2] },
            ],
            [self.sup_vert[0], self.sup_vert[1], self.sup_vert[2]],
        );
    }

    pub(super) fn compute(&mut self, vertex_indexes: &mut Vec<i32>) {
        let mut loop_edges: BTreeMap<i32, bool> = BTreeMap::new();
        for i in 0..3 {
            loop_edges.insert(self.sup_trian.link_at(i).abs(), true);
        }
        if vertex_indexes.len() > 0 {
            let first = vertex_indexes[0];
            self.create_triangles(first, &mut loop_edges);
            self.diag_stage("after_first_triangle");
            self.create_triangles_on_new_vertices(vertex_indexes);
        }
        self.remove_aux_elements();
        self.diag_stage("after_remove_aux");
    }

    pub(super) fn delete_triangle(&mut self, index: i32, loop_edges: &mut BTreeMap<i32, bool>) {
        if self.init_circles {
            self.circles.delete(index);
        }
        let element = self.mesh_data.get_element(index);
        self.mesh_data.remove_element(index);
        for i in 0..3 {
            let e = element.link_at(i).abs();
            let o = element.link_at(i) > 0;
            if !map_bind(loop_edges, e, o) {
                loop_edges.remove(&e);
                self.mesh_data.remove_link(e, false);
            }
        }
    }

    pub(super) fn remove_aux_elements(&mut self) {
        let mut loop_edges: BTreeMap<i32, bool> = BTreeMap::new();
        let mut elements: BTreeSet<i32> = BTreeSet::new();
        {
            let mut selector = DelaunSelector::new(&self.mesh_data);
            let sup = self.sup_vert.clone();
            for &sv in &sup {
                selector.neighbours_of_node(sv);
            }
            for &e in selector.elements() {
                elements.insert(e);
            }
        }
        // `TColStd_PackedMapOfInteger` iterates in sorted key order.
        let mut element_ids: Vec<i32> = elements.into_iter().collect();
        element_ids.sort_unstable();
        for id in element_ids {
            self.delete_triangle(id, &mut loop_edges);
        }
        let mut loop_keys: Vec<i32> = loop_edges.keys().copied().collect();
        loop_keys.sort_unstable();
        for &e in &loop_keys {
            if self.mesh_data.elements_connected_to(e).is_empty() {
                self.mesh_data.remove_link(e, false);
            }
        }
        let sup = self.sup_vert.clone();
        for &sv in &sup {
            self.mesh_data.remove_node(sv, false);
        }
    }

    pub(super) fn create_triangles(&mut self, vertex_index: i32, poly: &mut BTreeMap<i32, bool>) {
        let mut loop_edges: Vec<i32> = Vec::new();
        let mut external_edges: Vec<i32> = Vec::new();
        let vertex_coord = self.mesh_data.get_node(vertex_index).location.coord;

        let keys: Vec<i32> = poly.keys().copied().collect();
        for edge_id in keys {
            let is_positive = poly[&edge_id];
            let edge = self.mesh_data.get_link(edge_id);
            let mut nodes = [0i32; 3];
            if is_positive {
                nodes[0] = edge.first_node();
                nodes[2] = edge.last_node();
            } else {
                nodes[0] = edge.last_node();
                nodes[2] = edge.first_node();
            }
            nodes[1] = vertex_index;

            let first_vertex = self.mesh_data.get_node(nodes[0]).location.coord;
            let last_vertex = self.mesh_data.get_node(nodes[2]).location.coord;
            let mut edge_dir = last_vertex.subtracted(&first_vertex);
            let edge_len = edge_dir.modulus();
            if edge_len < PREC {
                continue;
            }
            edge_dir = edge_dir.divided(edge_len);

            let first_link_dir = first_vertex.subtracted(&vertex_coord);
            let last_link_dir = vertex_coord.subtracted(&last_vertex);

            let dist12 = first_link_dir.crossed(&edge_dir);
            let dist23 = edge_dir.crossed(&last_link_dir);
            if dist12.abs() < PREC || dist23.abs() < PREC {
                continue;
            }

            let first_link_id = self.mesh_data.add_link(nodes[1], nodes[0], VertexState::Free);
            let last_link_id = self.mesh_data.add_link(nodes[2], nodes[1], VertexState::Free);
            let edges_info = [first_link_id, if is_positive { edge_id } else { -edge_id }, last_link_id];

            let is_sens_ok = dist12 > 0.0 && dist23 > 0.0;
            if is_sens_ok {
                let mut edge_ids = [0i32; 3];
                let mut edge_oris = [false; 3];
                for k in 0..3 {
                    edge_ids[k] = edges_info[k].abs();
                    edge_oris[k] = edges_info[k] > 0;
                }
                self.add_triangle(edge_ids, edge_oris, nodes);
            } else {
                if is_positive {
                    loop_edges.push(edge_id);
                } else {
                    loop_edges.push(-edge_id);
                }
                if first_link_dir.square_modulus() > last_link_dir.square_modulus() {
                    external_edges.push(edges_info[0].abs());
                } else {
                    external_edges.push(edges_info[2].abs());
                }
            }
        }

        poly.clear();
        while !external_edges.is_empty() {
            let e = external_edges.remove(0);
            let pair = self.mesh_data.elements_connected_to(e);
            if !pair.is_empty() {
                self.delete_triangle(pair.first_index(), poly);
            }
        }

        let poly_keys: Vec<i32> = poly.keys().copied().collect();
        for &e in &poly_keys {
            if self.mesh_data.elements_connected_to(e).is_empty() {
                self.mesh_data.remove_link(e, false);
            }
        }

        while !loop_edges.is_empty() {
            let edge_info = loop_edges.remove(0);
            if self.mesh_data.link_movability(edge_info.abs()) != VertexState::Deleted {
                self.mesh_left_polygon_of(edge_info.abs(), edge_info > 0, &mut None);
            }
        }
    }

    pub(super) fn create_triangles_on_new_vertices(&mut self, vertex_indexes: &mut Vec<i32>) {
        let (tol_u, tol_v) = self.mesh_data.get_tolerance();
        let sq_tol = tol_u * tol_u + tol_v * tol_v;

        let upper = vertex_indexes.len();
        let mut i = 0;
        while i < upper {
            let vertex_idx = vertex_indexes[i];
            let vertex = *self.mesh_data.get_node(vertex_idx);
            let mut loop_edges: BTreeMap<i32, bool> = BTreeMap::new();
            let mut circles_list = self.circles.select(vertex.location.coord);

            let _on_edge_id = 0;
            let mut triangle_id = 0;
            let mut j = 0;
            while j < circles_list.len() {
                let t = circles_list[j];
                let mut edge_on = 0;
                if self.contains(t, &vertex, sq_tol, &mut edge_on) {
                    if edge_on != 0 && self.mesh_data.link_movability(edge_on) != VertexState::Free {
                        if vertex.state == VertexState::Free {
                            j += 1;
                            continue;
                        }
                    }
                    triangle_id = t;
                    circles_list.remove(j);
                    break;
                }
                j += 1;
            }

            if triangle_id > 0 {
                self.delete_triangle(triangle_id, &mut loop_edges);
                let mut is_modify = true;
                while is_modify && !circles_list.is_empty() {
                    is_modify = false;
                    let mut k = 0;
                    while k < circles_list.len() {
                        let t = circles_list[k];
                        let element = self.mesh_data.get_element(t);
                        if loop_edges.contains_key(&element.link_at(0).abs())
                            || loop_edges.contains_key(&element.link_at(1).abs())
                            || loop_edges.contains_key(&element.link_at(2).abs())
                        {
                            is_modify = true;
                            self.delete_triangle(t, &mut loop_edges);
                            circles_list.remove(k);
                            break;
                        }
                        k += 1;
                    }
                }
                self.create_triangles(vertex_idx, &mut loop_edges);
            }
            i += 1;
        }
        self.diag_stage("after_vertex_loop");

        // `ProcessConstraints()` is called UNCONDITIONALLY at the tail of
        // `BRepMesh_Delaun::createTrianglesOnNewVertices`
        // (`BRepMesh_Delaun.cxx:703`, body `insertInternalEdges();
        // frontierAdjust()`). The vertex loop (`cxx:629-700`) is followed by
        // this call rather than containing it, so it is reached even when the
        // vertex list is empty. `frontierAdjust` is the only mechanism that
        // closes leftover free-edge loops, so skipping it here leaves them
        // unclosed.
        //
        // A full-link failure (`add_triangle`, see below) aborts the rest of the
        // pass: OCCT's `Standard_OutOfRange` unwinds out of `perform()` and is
        // swallowed by `BRepMesh_BaseMeshAlgo.cxx:52-62`.
        if self.failed {
            return;
        }
        self.process_constraints();
    }

    /// `BRepMesh_Delaun::addTriangle`. Returns `false` when one of the three
    /// links already carries two triangles — OCCT's `BRepMesh_PairOfIndex::Append`
    /// (`BRepMesh_PairOfIndex.hxx:41`) throws `Standard_OutOfRange` there, which
    /// `BRepMesh_BaseMeshAlgo.cxx:52-62` swallows so the face ends up with no
    /// mesh. The port reports the same condition instead of panicking, and
    /// callers abort the current polygon (`self.failed`).
    pub(super) fn add_triangle(&mut self, edges: [i32; 3], oris: [bool; 3], _nodes: [i32; 3]) -> bool {
        for e in edges {
            if self.mesh_data.elements_connected_to(e.abs()).extent() >= 2 {
                self.failed = true;
                return false;
            }
        }
        let signed = [
            if oris[0] { edges[0] } else { -edges[0] },
            if oris[1] { edges[1] } else { -edges[1] },
            if oris[2] { edges[2] } else { -edges[2] },
        ];
        // The vertex order must match the orientation of the three links:
        // each link runs from the start of one side to the start of the next,
        // so the CCW vertex cycle is [start of link0, start of link1, start of link2].
        // OCCT derives the triangle's vertices from its links the same way; storing
        // the caller's `nodes` verbatim put them in the reversed order, which broke
        // `contains`'s point-in-triangle orientation test (all inserts were rejected).
        let verts = [
            {
                let l = self.mesh_data.get_link(edges[0].abs());
                if oris[0] { l.first_node() } else { l.last_node() }
            },
            {
                let l = self.mesh_data.get_link(edges[1].abs());
                if oris[1] { l.first_node() } else { l.last_node() }
            },
            {
                let l = self.mesh_data.get_link(edges[2].abs());
                if oris[2] { l.first_node() } else { l.last_node() }
            },
        ];
        let tri = DelaunTriangle::new(signed, verts);
        let new_id = self.mesh_data.add_element(tri);
        let mut is_added = true;
        if self.init_circles {
            is_added = self.circles.bind_circle(
                new_id,
                self.mesh_data.get_node(verts[0]).location.coord,
                self.mesh_data.get_node(verts[1]).location.coord,
                self.mesh_data.get_node(verts[2]).location.coord,
            );
        }
        if !is_added {
            self.mesh_data.remove_element(new_id);
        }
        true
    }

    pub(super) fn insert_internal_edges(&mut self) {
        let mut internal_edges: Vec<i32> = self.internal_edges().into_iter().collect();
        internal_edges.sort_unstable();
        for &link_index in &internal_edges {
            let pair = self.mesh_data.elements_connected_to(link_index);
            let mut is_go = [true, true];
            for tri_it in 1..=pair.extent() {
                let element = self.mesh_data.get_element(pair.index(tri_it));
                for i in 0..3 {
                    if element.link_at(i).abs() == link_index {
                        is_go[if element.link_at(i) > 0 { 0 } else { 1 }] = false;
                        break;
                    }
                }
            }
            if is_go[0] {
                self.mesh_left_polygon_of(link_index, true, &mut None);
            }
            if is_go[1] {
                self.mesh_left_polygon_of(link_index, false, &mut None);
            }
        }
    }

    pub(super) fn is_bound_to_frontier(&self, ref_node_id: i32, ref_link_id: i32) -> bool {
        let trace = self.diag_on()
            && std::env::var("FBPROBE").ok().and_then(|v| v.parse::<i32>().ok())
                == Some(ref_link_id);
        let mut stack: Vec<i32> = vec![ref_link_id];
        let mut visited: BTreeSet<i32> = BTreeSet::new();
        while let Some(cur) = stack.pop() {
            let pair = self.mesh_data.elements_connected_to(cur);
            if trace {
                eprintln!(
                    "FBPROBE node={ref_node_id} link={ref_link_id} pop={cur} extent={}",
                    pair.extent()
                );
            }
            if pair.is_empty() {
                return false;
            }
            for elem_it in 1..=pair.extent() {
                let tri_id = pair.index(elem_it);
                if tri_id < 0 {
                    continue;
                }
                let element = self.mesh_data.get_element(tri_id);
                for k in 0..3 {
                    let edge_id = element.link_at(k).abs();
                    if edge_id == cur {
                        continue;
                    }
                    let edge = self.mesh_data.get_link(edge_id);
                    if edge.first_node() != ref_node_id && edge.last_node() != ref_node_id {
                        continue;
                    }
                    if trace {
                        eprintln!(
                            "FBPROBE   tri={tri_id} edge={edge_id} mov={} hit",
                            self.mesh_data.link_movability(edge_id).to_str()
                        );
                    }
                    if self.mesh_data.link_movability(edge_id) != VertexState::Free {
                        return true;
                    }
                    if visited.insert(edge_id) {
                        stack.push(edge_id);
                    }
                }
            }
        }
        false
    }

    pub(super) fn cleanup_mesh(&mut self) {
        loop {
            let mut loop_edges: BTreeMap<i32, bool> = BTreeMap::new();
            let mut del_triangles: BTreeSet<i32> = BTreeSet::new();

            // `FreeEdges()` is `TColStd_PackedMapOfInteger` — sorted keys.
            let mut free_edges: Vec<i32> = self.free_edges().into_iter().collect();
            free_edges.sort_unstable();
            for free_edge_id in free_edges {
                let edge = self.mesh_data.get_link(free_edge_id);
                if self.mesh_data.link_movability(free_edge_id) == VertexState::Frontier {
                    continue;
                }
                let pair = self.mesh_data.elements_connected_to(free_edge_id);
                if pair.is_empty() {
                    loop_edges.insert(free_edge_id, true);
                    continue;
                }
                let tri_id = pair.first_index();
                let element = self.mesh_data.get_element(tri_id);
                let an_edges = element.link_indices;

                let mut can_not_be_removed = true;
                for cur_edge_idx in 0..3 {
                    if an_edges[cur_edge_idx].abs() != free_edge_id {
                        continue;
                    }
                    for other in 1..=2 {
                        if !can_not_be_removed {
                            break;
                        }
                        let other_edge_id = an_edges[(cur_edge_idx + other) % 3].abs();
                        let other_pair = self.mesh_data.elements_connected_to(other_edge_id);
                        if other_pair.extent() < 2 {
                            can_not_be_removed = false;
                        } else {
                            for tri_idx in 1..=other_pair.extent() {
                                if !can_not_be_removed {
                                    break;
                                }
                                if other_pair.index(tri_idx) == tri_id {
                                    continue;
                                }
                                let cur_triangle = self.mesh_data.get_element(other_pair.index(tri_idx));
                                for v in cur_triangle.vertex_indices {
                                    if self.is_sup_vertex(v) {
                                        can_not_be_removed = false;
                                    }
                                }
                            }
                        }
                    }
                    break;
                }
                if can_not_be_removed {
                    continue;
                }

                let mut is_connected = [false, false];
                for l in 0..2 {
                    let node = if l == 0 { edge.first_node() } else { edge.last_node() };
                    is_connected[l] = self.is_bound_to_frontier(node, free_edge_id);
                }
                if !is_connected[0] || !is_connected[1] {
                    del_triangles.insert(tri_id);
                }
                if self.diag_on() {
                    eprintln!(
                        "CMPROBE free={} mov={} ext={} tri={} cnr={} c0={} c1={} del={}",
                        free_edge_id,
                        self.mesh_data.link_movability(free_edge_id).to_str(),
                        pair.extent(),
                        tri_id,
                        can_not_be_removed,
                        is_connected[0],
                        is_connected[1],
                        !is_connected[0] || !is_connected[1]
                    );
                }
            }

            let mut deleted_nb = 0;
            let mut dels: Vec<i32> = del_triangles.iter().copied().collect();
            dels.sort_unstable();
            for id in dels {
                self.delete_triangle(id, &mut loop_edges);
                deleted_nb += 1;
            }

            let mut loop_keys: Vec<i32> = loop_edges.keys().copied().collect();
            loop_keys.sort_unstable();
            for &e in &loop_keys {
                if self.mesh_data.elements_connected_to(e).is_empty() {
                    self.mesh_data.remove_link(e, false);
                }
            }
            if deleted_nb == 0 {
                break;
            }
            self.diag_stage(&format!("cleanup_mesh_iter deleted={deleted_nb}"));
        }
        self.diag_stage("cleanup_mesh_done");
    }
}
