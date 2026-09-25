use super::prelude::*;
use super::*;

impl Delaun {

    pub(super) fn mesh_polygon(&mut self, the_polygon: &mut Vec<i32>, the_poly_boxes: &mut Vec<BndB2>, skipped: &mut Option<BTreeSet<i32>>) {
        // A full-link failure unwinds out of the mesh in OCCT
        // (`BRepMesh_BaseMeshAlgo.cxx:52-62` swallows the exception) — stop here.
        if self.failed {
            the_polygon.clear();
            the_poly_boxes.clear();
            return;
        }
        if self.mesh_elementary_polygon(the_polygon) {
            return;
        }

        let mut poly_len = the_polygon.len();
        let poly_area = self.poly_area(the_polygon, 1, poly_len).abs();
        let small_loop_area = 0.001 * poly_area;

        let mut a_poly_it: usize = 1;
        while a_poly_it < poly_len {
            let cur_edge_info = the_polygon[a_poly_it - 1];
            let mut cur_edge_id = cur_edge_info.abs();
            if self.mesh_data.link_movability(cur_edge_id) != VertexState::Frontier {
                a_poly_it += 1;
                continue;
            }
            let mut cur_nodes = [0i32; 2];
            self.get_oriented_nodes(cur_edge_id, cur_edge_info > 0, &mut cur_nodes);
            let mut cur_pnts = [
                self.mesh_data.get_node(cur_nodes[0]).location,
                self.mesh_data.get_node(cur_nodes[1]).location,
            ];

            let mut a_next_poly_it = a_poly_it + 1;
            while a_next_poly_it <= poly_len {
                let next_edge_info = the_polygon[a_next_poly_it - 1];
                let next_edge_id = next_edge_info.abs();
                if self.mesh_data.link_movability(next_edge_id) != VertexState::Frontier {
                    a_next_poly_it += 1;
                    continue;
                }
                let mut next_nodes = [0i32; 2];
                self.get_oriented_nodes(next_edge_id, next_edge_info > 0, &mut next_nodes);
                let next_pnts = [
                    self.mesh_data.get_node(next_nodes[0]).location,
                    self.mesh_data.get_node(next_nodes[1]).location,
                ];

                let mut int_pnt = GpPnt2d::zero();
                let int_flag = self.int_seg_seg(cur_edge_id, next_edge_id, false, true, &mut int_pnt);
                if int_flag == IntFlag::NoIntersection {
                    a_next_poly_it += 1;
                    continue;
                }

                let mut is_remove_from_first = false;
                let mut is_add_replacing_edge = true;
                let mut index_to_remove_to = a_next_poly_it;

                match int_flag {
                    IntFlag::Cross => {
                        let mut loop_area = self.poly_area(the_polygon, a_poly_it + 1, a_next_poly_it);
                        let vec1 = GpVec2d::from_xy(cur_pnts[1].coord.subtracted(&int_pnt.coord));
                        let vec2 = GpVec2d::from_xy(next_pnts[0].coord.subtracted(&int_pnt.coord));
                        loop_area += vec1.crossed(&vec2) / 2.0;
                        if loop_area.abs() > small_loop_area {
                            let next_nodes2 = [next_nodes[0], cur_nodes[0]];
                            let next_pnts2 = [next_pnts[0], cur_pnts[0]];
                            self.create_and_replace_polygon_link(
                                &next_nodes2,
                                &next_pnts2,
                                a_next_poly_it,
                                ReplaceFlag::Replace,
                                the_polygon,
                                the_poly_boxes,
                            );
                            self.process_loop(a_poly_it, a_next_poly_it, the_polygon, the_poly_boxes);
                            return;
                        }
                        let dist1 = int_pnt.square_distance(&next_pnts[0]);
                        let dist2 = int_pnt.square_distance(&next_pnts[1]);
                        let is_close_to_start = dist1 < dist2;
                        let end_point_index = if is_close_to_start { 0 } else { 1 };
                        cur_nodes[1] = next_nodes[end_point_index];
                        cur_pnts[1] = next_pnts[end_point_index];
                        if is_close_to_start {
                            index_to_remove_to -= 1;
                        }
                        if let Some(s) = skipped.as_mut() {
                            for sk in a_poly_it..=index_to_remove_to {
                                s.insert(the_polygon[sk - 1].abs());
                            }
                        }
                    }
                    IntFlag::PointOnSegment => {
                        let mut is_first_chopping = false;
                        let mut check_point_it = 0usize;
                        for cpi in 0..2 {
                            let ref_point = cur_pnts[cpi];
                            let v1 = GpVec2d::from_xy(next_pnts[0].coord.subtracted(&ref_point.coord));
                            let v2 = GpVec2d::from_xy(next_pnts[1].coord.subtracted(&ref_point.coord));
                            if v1.crossed(&v2).abs() < PREC {
                                is_first_chopping = true;
                                check_point_it = cpi;
                                break;
                            }
                        }
                        if is_first_chopping {
                            is_add_replacing_edge = false;
                            is_remove_from_first = check_point_it == 0;
                            let split_link = [next_nodes[0], cur_nodes[check_point_it], next_nodes[1]];
                            let split_pnts = [next_pnts[0], cur_pnts[check_point_it], next_pnts[1]];
                            for sli in 0..2 {
                                let flag = if sli == 0 { ReplaceFlag::Replace } else { ReplaceFlag::InsertAfter };
                                self.create_and_replace_polygon_link(
                                    &split_link[sli..sli + 2],
                                    &split_pnts[sli..sli + 2],
                                    a_next_poly_it,
                                    flag,
                                    the_polygon,
                                    the_poly_boxes,
                                );
                            }
                            self.process_loop(a_poly_it + check_point_it, index_to_remove_to, the_polygon, the_poly_boxes);
                        } else {
                            let split_link_nodes = [next_nodes[1], cur_nodes[1]];
                            let split_link_pnts = [next_pnts[1], cur_pnts[1]];
                            self.create_and_replace_polygon_link(
                                &split_link_nodes,
                                &split_link_pnts,
                                a_poly_it,
                                ReplaceFlag::InsertAfter,
                                the_polygon,
                                the_poly_boxes,
                            );
                            cur_nodes[1] = next_nodes[1];
                            cur_pnts[1] = next_pnts[1];
                            index_to_remove_to += 1;
                            self.process_loop(a_poly_it + 1, index_to_remove_to, the_polygon, the_poly_boxes);
                        }
                    }
                    IntFlag::Glued => {
                        if cur_nodes[1] == next_nodes[0] {
                            cur_nodes[1] = next_nodes[1];
                            cur_pnts[1] = next_pnts[1];
                        }
                    }
                    IntFlag::Same => {
                        self.process_loop(a_poly_it, a_next_poly_it, the_polygon, the_poly_boxes);
                        is_remove_from_first = true;
                        is_add_replacing_edge = false;
                    }
                    _ => {
                        a_next_poly_it += 1;
                        continue;
                    }
                }

                if is_add_replacing_edge {
                    cur_edge_id = self
                        .create_and_replace_polygon_link(
                            &cur_nodes,
                            &cur_pnts,
                            a_poly_it,
                            ReplaceFlag::Replace,
                            the_polygon,
                            the_poly_boxes,
                        )
                        .abs();
                }

                let index_to_remove_from = if is_remove_from_first { a_poly_it } else { a_poly_it + 1 };
                the_polygon.drain(index_to_remove_from - 1..=index_to_remove_to - 1);
                the_poly_boxes.drain(index_to_remove_from - 1..=index_to_remove_to - 1);
                poly_len = the_polygon.len();
                if is_remove_from_first {
                    a_poly_it -= 1;
                    break;
                }
                // cxx:2035 aNextPolyIt = aPolyIt; then for ++aNextPolyIt
                a_next_poly_it = a_poly_it + 1;
            }
            a_poly_it += 1;
        }

        // Decomposition of the (possibly corrected) polygon into triangles.
        let mut cut_poly = std::mem::take(the_polygon);
        let mut cut_boxes = std::mem::take(the_poly_boxes);
        let mut pending: Vec<(Vec<i32>, Vec<BndB2>)> = Vec::new();
        loop {
            let mut poly2: Vec<i32> = Vec::new();
            let mut boxes2: Vec<BndB2> = Vec::new();
            self.decompose_simple_polygon(&mut cut_poly, &mut cut_boxes, &mut poly2, &mut boxes2);
            if self.failed {
                // OCCT's exception unwinds the whole decomposition.
                break;
            }
            if !poly2.is_empty() {
                pending.push((poly2, boxes2));
            }
            if cut_poly.is_empty() {
                if pending.is_empty() {
                    break;
                }
                let (p, b) = pending.remove(0);
                cut_poly = p;
                cut_boxes = b;
            }
        }
    }

    pub(super) fn mesh_elementary_polygon(&mut self, the_polygon: &[i32]) -> bool {
        let poly_len = the_polygon.len();
        if poly_len < 3 {
            return true;
        } else if poly_len > 3 {
            return false;
        }
        let mut edges = [0i32; 3];
        let mut oris = [false; 3];
        for i in 0..3 {
            let info = the_polygon[i];
            edges[i] = info.abs();
            oris[i] = info > 0;
        }
        let edge1 = self.mesh_data.get_link(edges[0]);
        let edge2 = self.mesh_data.get_link(edges[1]);
        let mut nodes = [edge1.first_node(), edge1.last_node(), edge2.first_node()];
        if nodes[2] == nodes[0] || nodes[2] == nodes[1] {
            nodes[2] = edge2.last_node();
        }
        self.add_triangle(edges, oris, nodes);
        true
    }

    pub(super) fn decompose_simple_polygon(
        &mut self,
        the_polygon: &mut Vec<i32>,
        the_poly_boxes: &mut Vec<BndB2>,
        the_polygon_cut: &mut Vec<i32>,
        the_poly_boxes_cut: &mut Vec<BndB2>,
    ) {
        if self.mesh_elementary_polygon(the_polygon) {
            the_polygon.clear();
            the_poly_boxes.clear();
            return;
        }

        let poly_len = the_polygon.len();
        let first_edge_info = the_polygon[0];
        let _first_edge = self.mesh_data.get_link(first_edge_info.abs());
        let mut nodes = [0i32; 3];
        self.get_oriented_nodes(first_edge_info.abs(), first_edge_info > 0, &mut nodes[..2]);

        let mut ref_vertices = [
            self.mesh_data.get_node(nodes[0]).location,
            self.mesh_data.get_node(nodes[1]).location,
            GpPnt2d::zero(),
        ];
        let mut ref_edge_dir = GpVec2d::from_xy(ref_vertices[1].coord.subtracted(&ref_vertices[0].coord));
        let ref_edge_len = ref_edge_dir.magnitude();
        if ref_edge_len < PREC {
            the_polygon.clear();
            the_poly_boxes.clear();
            return;
        }
        ref_edge_dir = ref_edge_dir.divided(ref_edge_len);

        let mut used_link_id = 0usize;
        let mut opt_angle = 0.0;
        let mut min_dist = f64::MAX;
        let mut pivot_node = nodes[1];

        for link_it in 3..=poly_len {
            let link_info = the_polygon[link_it - 1];
            let next_edge = self.mesh_data.get_link(link_info.abs());
            pivot_node = if link_info > 0 { next_edge.first_node() } else { next_edge.last_node() };
            if pivot_node == nodes[1] {
                continue;
            }
            let pivot_vertex = self.mesh_data.get_node(pivot_node).location;
            let distance_dir = GpVec2d::from_xy(pivot_vertex.coord.subtracted(&ref_vertices[1].coord));
            let dist = ref_edge_dir.crossed(&distance_dir);
            let angle = ref_edge_dir.angle(&distance_dir).abs();
            let abs_dist = dist.abs();
            if abs_dist < PREC || dist < 0.0 {
                continue;
            }
            if (abs_dist >= min_dist) && (angle <= opt_angle || angle > ANG_DEV_90DEG) {
                continue;
            }

            let mut is_intersect = false;
            for ref_link_node_it in 0..2 {
                let link_first_node = nodes[ref_link_node_it];
                let link_first_vertex = ref_vertices[ref_link_node_it];
                let mut b = BndB2::void();
                update_bnd_box(link_first_vertex.coord, pivot_vertex.coord, &mut b);
                for check_link_it in 2..=poly_len {
                    if check_link_it == link_it {
                        continue;
                    }
                    if b.is_out(&the_poly_boxes[check_link_it - 1]) {
                        continue;
                    }
                    let poly_link = self.mesh_data.get_link(the_polygon[check_link_it - 1].abs());
                    let check_same = (link_first_node == poly_link.first_node() && pivot_node == poly_link.last_node())
                        || (link_first_node == poly_link.last_node() && pivot_node == poly_link.first_node());
                    if check_same {
                        continue;
                    }
                    let mut int_pnt = GpPnt2d::zero();
                    let flag = self.int_seg_seg_nodes(
                        link_first_node,
                        pivot_node,
                        poly_link.first_node(),
                        poly_link.last_node(),
                        false,
                        false,
                        &mut int_pnt,
                    );
                    if flag != IntFlag::NoIntersection {
                        is_intersect = true;
                        break;
                    }
                }
                if is_intersect {
                    break;
                }
            }
            if is_intersect {
                continue;
            }

            opt_angle = angle;
            min_dist = abs_dist;
            nodes[2] = pivot_node;
            ref_vertices[2] = pivot_vertex;
            used_link_id = link_it;
        }

        if used_link_id == 0 {
            the_polygon.clear();
            the_poly_boxes.clear();
            return;
        }

        let new_edge_0 = self.mesh_data.add_link(nodes[1], nodes[2], VertexState::Free);
        let new_edge_1 = self.mesh_data.add_link(nodes[2], nodes[0], VertexState::Free);
        let new_edges_info = [first_edge_info, new_edge_0, new_edge_1];
        // `BRepMesh_Delaun.cxx:2259-2274`: OCCT does `AddLink` twice and then
        // `addTriangle` **without touching existing triangles**. The previous body
        // deleted whichever neighbour triangle already used one of the ear links
        // ("free the slot before `AddElement`") — an invented rule (audit A25).
        // `add_triangle` now reports OCCT's `Standard_OutOfRange` condition
        // (= a link already carrying two triangles) instead of panicking, and the
        // polygon is dropped, mirroring the swallowed exception that leaves the
        // face unmeshed in OCCT.
        if !self.add_triangle_by_info(new_edges_info, nodes) {
            the_polygon.clear();
            the_poly_boxes.clear();
            return;
        }

        if used_link_id == 3 {
            the_polygon.remove(0);
            the_poly_boxes.remove(0);
            the_polygon[0] = -new_edges_info[2];
            let mut b = BndB2::void();
            update_bnd_box(ref_vertices[0].coord, ref_vertices[2].coord, &mut b);
            the_poly_boxes[0] = b;
        } else {
            if used_link_id < poly_len {
                // OCCT `NCollection_BaseSequence::PSplit(theIndex)` keeps
                // 1..theIndex-1 and moves theIndex..end (1-based, inclusive).
                let split_at = used_link_id - 1;
                the_polygon_cut.clear();
                the_polygon_cut.extend_from_slice(&the_polygon[split_at..]);
                the_polygon.truncate(split_at);
                the_polygon_cut.insert(0, -new_edges_info[2]);
                the_poly_boxes_cut.clear();
                the_poly_boxes_cut.extend_from_slice(&the_poly_boxes[split_at..]);
                the_poly_boxes.truncate(split_at);
                let mut b = BndB2::void();
                update_bnd_box(ref_vertices[0].coord, ref_vertices[2].coord, &mut b);
                the_poly_boxes_cut.insert(0, b);
            } else {
                the_polygon.pop();
                the_poly_boxes.pop();
            }
            the_polygon[0] = -new_edges_info[1];
            let mut b = BndB2::void();
            update_bnd_box(ref_vertices[1].coord, ref_vertices[2].coord, &mut b);
            the_poly_boxes[0] = b;
        }
    }

    /// Re-triangulates the cavity left after deleting the triangles of a vertex
    /// (used by `RemoveVertex`). Collects the loop of free edges around the
    /// cavity and meshes it as a polygon.
    pub(super) fn mesh_polygon_of_cavity(&mut self, loop_edges: &mut BTreeMap<i32, bool>) {
        let mut boxes: Vec<BndB2> = Vec::new();
        let mut polygon: Vec<i32> = Vec::new();
        let mut loop_edges_count = loop_edges.len();
        let keys: Vec<i32> = loop_edges.keys().copied().collect();
        if keys.is_empty() {
            return;
        }
        let an_edge_id = keys[0];
        let edge = self.mesh_data.get_link(an_edge_id);
        let mut first_node = edge.first_node();
        let mut pivot_node = edge.last_node();
        let mut edge_id = an_edge_id;

        let is_positive = loop_edges[&edge_id];
        if !is_positive {
            let tmp = first_node;
            first_node = pivot_node;
            pivot_node = tmp;
            polygon.push(-edge_id);
        } else {
            polygon.push(edge_id);
        }
        let mut b = BndB2::void();
        update_bnd_box(
            self.mesh_data.get_node(first_node).location.coord,
            self.mesh_data.get_node(pivot_node).location.coord,
            &mut b,
        );
        boxes.push(b);
        loop_edges.remove(&edge_id);

        let last_node = first_node;
        while pivot_node != last_node {
            let links: Vec<i32> = self.mesh_data.links_connected_to(pivot_node).to_vec();
            let mut advanced = false;
            for &link_value in &links {
                if link_value != edge_id && loop_edges.contains_key(&link_value) {
                    edge_id = link_value;
                    let e = self.mesh_data.get_link(edge_id);
                    // Traverse the next cavity edge away from the pivot node: if
                    // the link is stored pivot->last use it forward, otherwise
                    // reverse it so the polygon advances to the other endpoint.
                    let current_node = if e.first_node() == pivot_node {
                        polygon.push(edge_id);
                        e.last_node()
                    } else {
                        polygon.push(-edge_id);
                        e.first_node()
                    };
                    let mut b2 = BndB2::void();
                    update_bnd_box(
                        self.mesh_data.get_node(current_node).location.coord,
                        self.mesh_data.get_node(pivot_node).location.coord,
                        &mut b2,
                    );
                    boxes.push(b2);
                    pivot_node = current_node;
                    loop_edges.remove(&edge_id);
                    advanced = true;
                    break;
                }
            }
            if !advanced {
                break;
            }
            if loop_edges_count <= 0 {
                break;
            }
            loop_edges_count -= 1;
        }
        if polygon.len() >= 3 {
            self.mesh_polygon(&mut polygon, &mut boxes, &mut None);
        }
    }

    /// Edges whose state matches `edge_type`, in ascending link-id order.
    /// Source: `BRepMesh_MeshTool::GetEdgesByType` (`BRepMesh_MeshTool.cxx:284-300`),
    /// which walks `LinksOfDomain()` (a `NCollection_PackedMap<int>`) in
    /// ascending order and stores into another `IMeshData::MapOfInteger`.
    pub(super) fn get_edges_by_type(&self, edge_type: VertexState) -> BTreeSet<i32> {
        let mut result = BTreeSet::new();
        for &edge in self.mesh_data.links_of_domain() {
            let is_to_add = if edge_type == VertexState::Free {
                self.mesh_data.elements_connected_to(edge).extent() <= 1
            } else {
                self.mesh_data.link_movability(edge) == edge_type
            };
            if is_to_add {
                result.insert(edge);
            }
        }
        result
    }

    pub(super) fn calculate_dist(
        &self,
        v_edges: &[GpXY; 3],
        points: &[GpPnt2d; 3],
        vertex: &DelaunVertex,
        distance: &mut [f64; 3],
        sq_modulus: &mut [f64; 3],
        edge_on: &mut usize,
    ) -> f64 {
        let mut min_dist = f64::MAX;
        let v = vertex.location.coord;
        for i in 0..3 {
            sq_modulus[i] = v_edges[i].square_modulus();
            if sq_modulus[i] <= PREC2 {
                return -1.0;
            }
            distance[i] = v_edges[i].crossed(&v.subtracted(&points[i].coord));
            let d = distance[i] * distance[i] / sq_modulus[i];
            if d < min_dist {
                *edge_on = i;
                min_dist = d;
            }
        }
        min_dist
    }

    pub(super) fn int_seg_seg(
        &self,
        edge1_id: i32,
        edge2_id: i32,
        is_consider_end_point_touch: bool,
        is_consider_point_on_edge: bool,
        int_pnt: &mut GpPnt2d,
    ) -> IntFlag {
        let e1 = self.mesh_data.get_link(edge1_id);
        let e2 = self.mesh_data.get_link(edge2_id);
        self.int_seg_seg_nodes(
            e1.first_node(),
            e1.last_node(),
            e2.first_node(),
            e2.last_node(),
            is_consider_end_point_touch,
            is_consider_point_on_edge,
            int_pnt,
        )
    }

    pub(super) fn int_seg_seg_nodes(
        &self,
        n1: i32,
        n2: i32,
        n3: i32,
        n4: i32,
        is_consider_end_point_touch: bool,
        is_consider_point_on_edge: bool,
        int_pnt: &mut GpPnt2d,
    ) -> IntFlag {
        let p1 = self.mesh_data.get_node(n1).location.coord;
        let p2 = self.mesh_data.get_node(n2).location.coord;
        let p3 = self.mesh_data.get_node(n3).location.coord;
        let p4 = self.mesh_data.get_node(n4).location.coord;
        let (flag, pnt) = GeomTool::int_seg_seg(&p1, &p2, &p3, &p4, is_consider_end_point_touch, is_consider_point_on_edge);
        *int_pnt = GpPnt2d::from_xy(pnt);
        flag
    }

    pub(super) fn poly_area(&self, the_polygon: &[i32], start_1based: usize, end_1based: usize) -> f64 {
        let mut area = 0.0;
        let poly_len = the_polygon.len();
        if start_1based >= end_1based || start_1based > poly_len {
            return area;
        }
        let mut cur_edge_info = the_polygon[start_1based - 1];
        let mut cur_edge_id = cur_edge_info.abs();
        let mut nodes = [0i32; 2];
        self.get_oriented_nodes(cur_edge_id, cur_edge_info > 0, &mut nodes);
        let ref_pnt = self.mesh_data.get_node(nodes[0]).location;
        for poly_it in (start_1based + 1)..=end_1based {
            cur_edge_info = the_polygon[poly_it - 1];
            cur_edge_id = cur_edge_info.abs();
            self.get_oriented_nodes(cur_edge_id, cur_edge_info > 0, &mut nodes);
            let v1 = GpVec2d::from_xy(self.mesh_data.get_node(nodes[0]).location.coord.subtracted(&ref_pnt.coord));
            let v2 = GpVec2d::from_xy(self.mesh_data.get_node(nodes[1]).location.coord.subtracted(&ref_pnt.coord));
            area += v1.crossed(&v2);
        }
        area / 2.0
    }

    pub(super) fn is_sup_vertex(&self, vertex_idx: i32) -> bool {
        self.sup_vert.contains(&vertex_idx)
    }
}
