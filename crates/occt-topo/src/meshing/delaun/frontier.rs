use super::prelude::*;
use super::*;

impl Delaun {

    pub(super) fn frontier_adjust(&mut self) {
        // `Frontier()` walks `LinksOfDomain()`, an OCCT `NCollection_PackedMap<int>`
        // traversed in ascending id order (`BRepMesh_MeshTool.cxx:284-300`); our
        // structure stores it in a `BTreeSet`, so the ids are already ascending.
        let frontier_ids: Vec<i32> = self.frontier().into_iter().collect();
        let mut failed_frontiers: Vec<i32> = Vec::new();
        let mut loop_edges: BTreeMap<i32, bool> = BTreeMap::new();
        let mut int_frontier_edges: BTreeSet<i32> = BTreeSet::new();

        for _pass in 1..=2 {
            self.diag_stage(&format!("pass{_pass}_start"));
            for &frontier_id in &frontier_ids {
                let pair = self.mesh_data.elements_connected_to(frontier_id);
                let nb = pair.extent();
                for elem_it in 1..=nb {
                    let prior_elem = pair.index(elem_it);
                    if prior_elem < 0 {
                        continue;
                    }
                    let element = self.mesh_data.get_element(prior_elem);
                    let mut found = false;
                    for n in 0..3 {
                        if frontier_id == element.link_at(n).abs() {
                            if element.link_at(n) < 0 {
                                found = true;
                            }
                        }
                    }
                    for n in 0..3 {
                        if frontier_id == element.link_at(n).abs() && element.link_at(n) < 0 {
                            self.delete_triangle(prior_elem, &mut loop_edges);
                            break;
                        }
                    }
                    if found {
                        break;
                    }
                }
            }

            self.diag_stage(&format!("pass{_pass}_after_delete"));

            let loop_keys: Vec<i32> = loop_edges.keys().copied().collect();
            for &e in &loop_keys {
                if self.mesh_data.elements_connected_to(e).is_empty() {
                    self.mesh_data.remove_link(e, false);
                }
            }
            self.diag_stage(&format!("pass{_pass}_after_cleanup_loop_edges"));

            for &frontier_id in &frontier_ids {
                if !self.mesh_data.elements_connected_to(frontier_id).is_empty() {
                    continue;
                }
                let mut skipped = Some(int_frontier_edges.clone());
                let success = self.mesh_left_polygon_of(frontier_id, true, &mut skipped);
                if let Some(s) = skipped {
                    int_frontier_edges = s;
                }
                if _pass == 2 && !success {
                    failed_frontiers.push(frontier_id);
                }
            }
            self.diag_stage(&format!("pass{_pass}_after_mesh_left_polygon"));
        }

        self.cleanup_mesh();
        self.diag_stage("after_cleanup_mesh");

        for &frontier_id in failed_frontiers.iter() {
            if !self.mesh_data.elements_connected_to(frontier_id).is_empty() {
                continue;
            }
            let mut skipped = Some(int_frontier_edges.clone());
            self.mesh_left_polygon_of(frontier_id, true, &mut skipped);
            if let Some(s) = skipped {
                int_frontier_edges = s;
            }
        }
        self.diag_stage("after_failed_retries");
    }

    pub(super) fn fill_bnd_box(&self, boxes: &mut Vec<BndB2>, v1: i32, v2: i32) {
        let mut b = BndB2::void();
        update_bnd_box(
            self.mesh_data.get_node(v1).location.coord,
            self.mesh_data.get_node(v2).location.coord,
            &mut b,
        );
        boxes.push(b);
    }

    pub(super) fn mesh_left_polygon_of(&mut self, start_edge_id: i32, is_forward: bool, skipped: &mut Option<BTreeSet<i32>>) -> bool {
        if let Some(s) = skipped.as_ref() {
            if s.contains(&start_edge_id) {
                self.diag_stage(&format!("mesh_left_polygon_of e={start_edge_id} fwd={is_forward} SKIPPED"));
                return true;
            }
        }
        let ref_edge = self.mesh_data.get_link(start_edge_id);

        let mut polygon: Vec<i32> = Vec::new();
        let (mut a_start_node, mut a_pivot_node);
        if is_forward {
            polygon.push(start_edge_id);
            a_start_node = ref_edge.first_node();
            a_pivot_node = ref_edge.last_node();
        } else {
            polygon.push(-start_edge_id);
            a_start_node = ref_edge.last_node();
            a_pivot_node = ref_edge.first_node();
        }
        let start_edge_vertex_s = self.mesh_data.get_node(a_start_node).location.coord;
        let mut a_pivot_vertex = self.mesh_data.get_node(a_pivot_node).location.coord;
        let mut ref_link_dir = GpVec2d::from_xy(a_pivot_vertex.subtracted(&start_edge_vertex_s));
        if ref_link_dir.square_magnitude() < PREC2 {
            self.diag_stage(&format!("mesh_left_polygon_of e={start_edge_id} fwd={is_forward} DEGENERATE"));
            return true;
        }

        let mut boxes: Vec<BndB2> = Vec::new();
        let mut b0 = BndB2::void();
        update_bnd_box(start_edge_vertex_s, a_pivot_vertex, &mut b0);
        boxes.push(b0);

        let mut dead_links: BTreeSet<i32> = BTreeSet::new();
        let mut leprous_links: BTreeSet<i32> = BTreeSet::new();
        leprous_links.insert(start_edge_id);

        let mut is_skip_leprous = true;
        // OCCT `BRepMesh_Delaun.cxx:1116-1140`: `aFirstNode` is the loop-close
        // target and `findNextPolygonLink`'s endpoint check; it never updates.
        // `aStartNode` walks with the polygon and is the origin of `aRefLinkDir`.
        let a_first_node = a_start_node;
        while a_pivot_node != a_first_node {
            let result = self.find_next_polygon_link(
                a_first_node,
                a_pivot_node,
                a_pivot_vertex,
                ref_link_dir,
                &boxes,
                &polygon,
                skipped.as_ref(),
                is_skip_leprous,
                &mut leprous_links,
                &mut dead_links,
            );
            if let Some((next_link_id, next_pivot_node, next_link_dir, next_link_bbox)) = result {
                a_start_node = a_pivot_node;
                ref_link_dir = next_link_dir;
                a_pivot_node = next_pivot_node;
                a_pivot_vertex = self.mesh_data.get_node(next_pivot_node).location.coord;
                boxes.push(next_link_bbox);
                polygon.push(next_link_id);
                is_skip_leprous = true;
            } else {
                if polygon.len() == 1 {
                    self.diag_stage(&format!("mesh_left_polygon_of e={start_edge_id} fwd={is_forward} DEAD_END"));
                    return false;
                }
                let dead_link_id = polygon.last().copied().unwrap().abs();
                dead_links.insert(dead_link_id);
                leprous_links.remove(&dead_link_id);
                polygon.pop();
                boxes.pop();

                let prev_link_info = *polygon.last().unwrap();
                let prev_link = self.mesh_data.get_link(prev_link_info.abs());
                if prev_link_info > 0 {
                    a_start_node = prev_link.first_node();
                    a_pivot_node = prev_link.last_node();
                } else {
                    a_start_node = prev_link.last_node();
                    a_pivot_node = prev_link.first_node();
                }
                a_pivot_vertex = self.mesh_data.get_node(a_pivot_node).location.coord;
                let sn = self.mesh_data.get_node(a_start_node).location.coord;
                ref_link_dir = GpVec2d::from_xy(a_pivot_vertex.subtracted(&sn));
                is_skip_leprous = false;
            }
        }

        if polygon.len() < 3 {
            return false;
        }
        self.cleanup_polygon(&polygon, &boxes);
        self.mesh_polygon(&mut polygon, &mut boxes, skipped);
        true
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn find_next_polygon_link(
        &self,
        first_node: i32,
        pivot_node: i32,
        pivot_vertex: GpXY,
        ref_link_dir: GpVec2d,
        boxes: &[BndB2],
        polygon: &[i32],
        skipped: Option<&BTreeSet<i32>>,
        is_skip_leprous: bool,
        leprous_links: &mut BTreeSet<i32>,
        dead_links: &mut BTreeSet<i32>,
    ) -> Option<(i32, i32, GpVec2d, BndB2)> {
        let mut max_angle = f64::NEG_INFINITY;
        let mut next_link_id = 0;
        let mut next_pivot_node = 0;
        let mut next_link_dir = GpVec2d::zero();
        let mut next_link_bbox = BndB2::void();

        let neighbors: Vec<i32> = self.mesh_data.links_connected_to(pivot_node).to_vec();
        for &neighbor_info in &neighbors {
            let neighbor_id = neighbor_info.abs();
            if dead_links.contains(&neighbor_id) {
                continue;
            }
            if let Some(s) = skipped {
                if s.contains(&neighbor_id) {
                    continue;
                }
            }
            let is_leprous = leprous_links.contains(&neighbor_id);
            if is_skip_leprous && is_leprous {
                continue;
            }
            let neighbor_link = self.mesh_data.get_link(neighbor_id);
            if self.mesh_data.link_movability(neighbor_id) == VertexState::Free
                && self.mesh_data.elements_connected_to(neighbor_id).is_empty()
            {
                dead_links.insert(neighbor_id);
                continue;
            }
            let mut other_node = neighbor_link.first_node();
            if other_node == pivot_node {
                other_node = neighbor_link.last_node();
            }
            let cur_link_dir = GpVec2d::from_xy(self.mesh_data.get_node(other_node).location.coord.subtracted(&pivot_vertex));
            if cur_link_dir.square_magnitude() < PREC2 {
                dead_links.insert(neighbor_id);
                continue;
            }
            if !is_leprous {
                leprous_links.insert(neighbor_id);
            }

            let mut angle = ref_link_dir.angle(&cur_link_dir);
            let is_frontier = self.mesh_data.link_movability(neighbor_id) == VertexState::Frontier;
            let mut is_check_point_on_edge = true;
            if is_frontier {
                if (angle.abs() - PI).abs() < ANGULAR {
                    is_check_point_on_edge = false;
                    angle = angle.abs();
                }
            }
            if angle <= max_angle {
                continue;
            }

            let is_check_end_points = other_node != first_node;
            let mut a_box = BndB2::void();
            let is_not_intersect = self.check_intersection(
                neighbor_id,
                polygon,
                boxes,
                is_check_end_points,
                is_check_point_on_edge,
                true,
                &mut a_box,
            );
            if is_not_intersect {
                max_angle = angle;
                next_link_dir = cur_link_dir;
                next_pivot_node = other_node;
                next_link_bbox = a_box;
                next_link_id = if neighbor_link.first_node() == pivot_node { neighbor_id } else { -neighbor_id };
            }
        }

        if next_link_id == 0 {
            None
        } else {
            Some((next_link_id, next_pivot_node, next_link_dir, next_link_bbox))
        }
    }

    pub(super) fn check_intersection(
        &self,
        link_id: i32,
        polygon: &[i32],
        poly_boxes: &[BndB2],
        is_consider_end_point_touch: bool,
        is_consider_point_on_edge: bool,
        is_skip_last_edge: bool,
        link_bbox: &mut BndB2,
    ) -> bool {
        let link = self.mesh_data.get_link(link_id);
        update_bnd_box(
            self.mesh_data.get_node(link.first_node()).location.coord,
            self.mesh_data.get_node(link.last_node()).location.coord,
            link_bbox,
        );
        let mut poly_len = polygon.len();
        if is_skip_last_edge {
            poly_len -= 1;
        }
        let is_frontier = self.mesh_data.link_movability(link_id) == VertexState::Frontier;

        for poly_it in 0..poly_len {
            if !link_bbox.is_out(&poly_boxes[poly_it]) {
                let poly_link_id = polygon[poly_it].abs();
                let _poly_link = self.mesh_data.get_link(poly_link_id);
                if self.mesh_data.link_movability(poly_link_id) == VertexState::Frontier && is_frontier {
                    continue;
                }
                let mut int_pnt = GpPnt2d::zero();
                let flag = self.int_seg_seg(link_id, poly_link_id, is_consider_end_point_touch, is_consider_point_on_edge, &mut int_pnt);
                if flag != IntFlag::NoIntersection {
                    return false;
                }
            }
        }
        true
    }

    /// Returns `false` when `add_triangle` hit OCCT's full-link condition
    /// (`BRepMesh_PairOfIndex::Append` throwing `Standard_OutOfRange`).
    pub(super) fn add_triangle_by_info(&mut self, edges_info: [i32; 3], nodes: [i32; 3]) -> bool {
        let mut edges = [0i32; 3];
        let mut oris = [false; 3];
        for i in 0..3 {
            edges[i] = edges_info[i].abs();
            oris[i] = edges_info[i] > 0;
        }
        self.add_triangle(edges, oris, nodes)
    }

    pub(super) fn cleanup_polygon(&mut self, the_polygon: &[i32], the_poly_boxes: &[BndB2]) {
        let poly_len = the_polygon.len();
        if poly_len < 3 {
            return;
        }
        let mut loop_edges: BTreeMap<i32, bool> = BTreeMap::new();
        let mut ignored_edges: BTreeSet<i32> = BTreeSet::new();
        let mut poly_vertices_find_map: BTreeSet<i32> = BTreeSet::new();
        let mut poly_vertices: Vec<i32> = Vec::new();

        for poly_it in 1..=poly_len {
            let poly_edge_info = the_polygon[poly_it - 1];
            let poly_edge_id = poly_edge_info.abs();
            ignored_edges.insert(poly_edge_id);
            let is_forward = poly_edge_info > 0;
            let mut elem_it = 1;
            while elem_it <= self.mesh_data.elements_connected_to(poly_edge_id).extent() {
                let pair = self.mesh_data.elements_connected_to(poly_edge_id);
                let elem_id = pair.index(elem_it);
                if elem_id < 0 {
                    elem_it += 1;
                    continue;
                }
                let element = self.mesh_data.get_element(elem_id);
                let mut found = false;
                for k in 0..3 {
                    if element.link_at(k).abs() == poly_edge_id && (element.link_at(k) > 0) == is_forward {
                        found = true;
                        self.delete_triangle(elem_id, &mut loop_edges);
                        break;
                    }
                }
                if found {
                    break;
                }
                elem_it += 1;
            }

            if poly_it % 2 == 1 {
                let poly_edge = self.mesh_data.get_link(poly_edge_id);
                let first_vertex = poly_edge.first_node();
                let last_vertex = poly_edge.last_node();
                poly_vertices_find_map.insert(first_vertex);
                poly_vertices_find_map.insert(last_vertex);
                if poly_edge_info > 0 {
                    poly_vertices.push(first_vertex);
                    poly_vertices.push(last_vertex);
                } else {
                    poly_vertices.push(last_vertex);
                    poly_vertices.push(first_vertex);
                }
            }
        }

        if poly_vertices.first() != poly_vertices.last() {
            let first = *poly_vertices.first().unwrap();
            poly_vertices.push(first);
        }

        let mut survived_links = ignored_edges.clone();
        for poly_vert_it in 0..(poly_vertices.len() - 1) {
            let mut stack_frames = StackOfFrames::new();
            let mut stack_data: Vec<i32> = Vec::new();
            let mut current_victim = poly_vertices[poly_vert_it];
            loop {
                let prev_size = stack_data.len();
                self.kill_triangles_around_vertex(
                    current_victim,
                    &poly_vertices,
                    &poly_vertices_find_map,
                    the_polygon,
                    the_poly_boxes,
                    &mut survived_links,
                    &mut loop_edges,
                    &mut stack_data,
                );
                let new_size = stack_data.len();
                if new_size > prev_size {
                    stack_frames.push_frame(prev_size, new_size);
                }
                if stack_frames.is_empty() {
                    break;
                }
                current_victim = stack_data[stack_frames.pop_element()];
            }
        }

        let loop_keys: Vec<i32> = loop_edges.keys().copied().collect();
        for &e in &loop_keys {
            if ignored_edges.contains(&e) {
                continue;
            }
            if self.mesh_data.elements_connected_to(e).is_empty() {
                self.mesh_data.remove_link(e, false);
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn kill_triangles_around_vertex(
        &mut self,
        zombie_node_id: i32,
        poly_vertices: &[i32],
        poly_vertices_find_map: &BTreeSet<i32>,
        the_polygon: &[i32],
        the_poly_boxes: &[BndB2],
        survived_links: &mut BTreeSet<i32>,
        loop_edges: &mut BTreeMap<i32, bool>,
        victim_nodes: &mut Vec<i32>,
    ) {
        let neighbors: Vec<i32> = self.mesh_data.links_connected_to(zombie_node_id).to_vec();
        for &neighbor_link_id in &neighbors {
            if survived_links.contains(&neighbor_link_id) {
                continue;
            }
            let neighbor_link = self.mesh_data.get_link(neighbor_link_id);
            if self.mesh_data.link_movability(neighbor_link_id) == VertexState::Frontier {
                let mut b = BndB2::void();
                let is_not_intersect =
                    self.check_intersection(neighbor_link_id, the_polygon, the_poly_boxes, false, true, false, &mut b);
                if is_not_intersect {
                    survived_links.insert(neighbor_link_id);
                    continue;
                }
            } else {
                let mut other_node = neighbor_link.first_node();
                if other_node == zombie_node_id {
                    other_node = neighbor_link.last_node();
                }
                if !poly_vertices_find_map.contains(&other_node) {
                    if self.is_vertex_inside_polygon(other_node, poly_vertices) {
                        victim_nodes.push(other_node);
                    } else {
                        self.kill_triangles_on_intersecting_links(
                            neighbor_link_id,
                            other_node,
                            the_polygon,
                            the_poly_boxes,
                            survived_links,
                            loop_edges,
                        );
                        continue;
                    }
                }
            }
            survived_links.insert(neighbor_link_id);
            self.kill_link_triangles(neighbor_link_id, loop_edges);
        }
    }

    pub(super) fn is_vertex_inside_polygon(&self, vertex_id: i32, polygon_vertices: &[i32]) -> bool {
        let poly_len = polygon_vertices.len();
        if poly_len < 3 {
            return false;
        }
        let center = self.mesh_data.get_node(vertex_id).location;
        let first_vertex = self.mesh_data.get_node(polygon_vertices[0]).location;
        let mut prev_dir = GpVec2d::from_xy(first_vertex.coord.subtracted(&center.coord));
        if prev_dir.square_magnitude() < PREC2 {
            return true;
        }
        let mut total_ang = 0.0;
        for poly_it in 1..poly_len {
            let poly_vertex = self.mesh_data.get_node(polygon_vertices[poly_it]).location;
            let cur_dir = GpVec2d::from_xy(poly_vertex.coord.subtracted(&center.coord));
            if cur_dir.square_magnitude() < PREC2 {
                return true;
            }
            total_ang += cur_dir.angle(&prev_dir);
            prev_dir = cur_dir;
        }
        (total_ang.abs() - ANGLE_2PI).abs() <= ANGULAR
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn kill_triangles_on_intersecting_links(
        &mut self,
        link_to_check_id: i32,
        end_point: i32,
        the_polygon: &[i32],
        the_poly_boxes: &[BndB2],
        survived_links: &mut BTreeSet<i32>,
        loop_edges: &mut BTreeMap<i32, bool>,
    ) {
        if survived_links.contains(&link_to_check_id) {
            return;
        }
        let mut b = BndB2::void();
        let is_not_intersect =
            self.check_intersection(link_to_check_id, the_polygon, the_poly_boxes, false, false, false, &mut b);
        survived_links.insert(link_to_check_id);
        if is_not_intersect {
            return;
        }
        self.kill_link_triangles(link_to_check_id, loop_edges);
        let neighbors: Vec<i32> = self.mesh_data.links_connected_to(end_point).to_vec();
        for &neighbor_id in &neighbors {
            let neighbor_link = self.mesh_data.get_link(neighbor_id);
            let mut other_node = neighbor_link.first_node();
            if other_node == end_point {
                other_node = neighbor_link.last_node();
            }
            self.kill_triangles_on_intersecting_links(neighbor_id, other_node, the_polygon, the_poly_boxes, survived_links, loop_edges);
        }
    }

    pub(super) fn kill_link_triangles(&mut self, link_id: i32, loop_edges: &mut BTreeMap<i32, bool>) {
        let elem_nb = self.mesh_data.elements_connected_to(link_id).extent();
        for _ in 0..elem_nb {
            let elem_id = self.mesh_data.elements_connected_to(link_id).first_index();
            if elem_id < 0 {
                continue;
            }
            self.delete_triangle(elem_id, loop_edges);
        }
    }

    pub(super) fn get_oriented_nodes(&self, edge_id: i32, is_forward: bool, nodes: &mut [i32]) {
        let edge = self.mesh_data.get_link(edge_id);
        if is_forward {
            nodes[0] = edge.first_node();
            nodes[1] = edge.last_node();
        } else {
            nodes[0] = edge.last_node();
            nodes[1] = edge.first_node();
        }
    }

    pub(super) fn process_loop(&mut self, link_from: usize, link_to: usize, polygon: &[i32], poly_boxes: &[BndB2]) {
        // OCCT uses signed `int aNbOfLinksInLoop = theLinkTo - theLinkFrom - 1`.
        let nb = link_to as isize - link_from as isize - 1;
        if nb < 3 {
            return;
        }
        let nb = nb as usize;
        let mut sub_polygon: Vec<i32> = Vec::new();
        let mut sub_boxes: Vec<BndB2> = Vec::new();
        // OCCT Prepends `thePolygon(theLinkFrom + k)` for k = nb..1, yielding
        // 1-based indices theLinkFrom+1 .. theLinkTo-1 (= 0-based link_from .. link_to-2).
        for i in 0..nb {
            sub_polygon.push(polygon[link_from + i]);
            sub_boxes.push(poly_boxes[link_from + i]);
        }
        self.mesh_polygon(&mut sub_polygon, &mut sub_boxes, &mut None);
    }

    pub(super) fn create_and_replace_polygon_link(
        &mut self,
        nodes: &[i32],
        pnts: &[GpPnt2d],
        root_index: usize,
        flag: ReplaceFlag,
        polygon: &mut Vec<i32>,
        poly_boxes: &mut Vec<BndB2>,
    ) -> i32 {
        let new_edge_id = self.mesh_data.add_link(nodes[0], nodes[1], VertexState::Free);
        let mut new_box = BndB2::void();
        update_bnd_box(pnts[0].coord, pnts[1].coord, &mut new_box);
        match flag {
            ReplaceFlag::Replace => {
                polygon[root_index - 1] = new_edge_id;
                poly_boxes[root_index - 1] = new_box;
            }
            ReplaceFlag::InsertAfter => {
                polygon.insert(root_index, new_edge_id);
                poly_boxes.insert(root_index, new_box);
            }
            ReplaceFlag::InsertBefore => {
                polygon.insert(root_index - 1, new_edge_id);
                poly_boxes.insert(root_index - 1, new_box);
            }
        }
        new_edge_id
    }
}
