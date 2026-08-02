//! Auxiliary tool for merging triangulation nodes within a tolerance,
//! remapping triangle indices, optionally splitting nodes on sharp corners by
//! an angle.
//! Source: `Poly_MergeNodesTool.hxx` / `Poly_MergeNodesTool.cxx`.

use std::collections::HashMap;

use crate::gp::{GpPnt, GpTrsf, GpXyz};
use crate::precision::CONFUSION;

use super::triangulation::Triangle;
use super::triangulation_full::PolyTriangulation;

/// A key in the merged-nodes map.
#[derive(Debug, Clone, Copy)]
struct MergedKey {
    pos: [f64; 3],
    normal: [f64; 3],
}

/// Spatial hash map merging positions (and optionally normals) within a
/// tolerance. Ported from `Poly_MergeNodesTool::MergedNodesMap`.
#[derive(Debug, Clone)]
struct MergedNodesMap {
    tolerance: f64,
    inv_tol: f64,
    angle: f64,
    angle_cos: f64,
    to_merge_opposite: bool,
    merged: Vec<MergedKey>,
    buckets: HashMap<(i64, i64, i64), Vec<usize>>,
}

impl MergedNodesMap {
    fn new() -> Self {
        Self {
            tolerance: 0.0,
            inv_tol: 0.0,
            angle: 0.0,
            angle_cos: 1.0, // angle = 0 -> require identical normals
            to_merge_opposite: false,
            merged: Vec::new(),
            buckets: HashMap::new(),
        }
    }

    fn set_merge_tolerance(&mut self, tolerance: f64) {
        self.tolerance = tolerance;
        self.inv_tol = 0.0;
        if tolerance > 0.0 {
            self.inv_tol = 1.0 / tolerance;
        }
    }

    fn set_merge_angle(&mut self, angle: f64) {
        self.angle = angle;
        self.angle_cos = angle.cos();
    }

    fn has_merge_angle(&self) -> bool { self.angle > 0.0 }
    fn has_merge_tolerance(&self) -> bool { self.tolerance > 0.0 }
    /// Angle is close to 90 degrees: normal comparison can be skipped.
    fn to_merge_any_angle(&self) -> bool { self.angle_cos <= 0.01 }

    fn cell(&self, pos: &[f64; 3]) -> (i64, i64, i64) {
        (
            (pos[0] * self.inv_tol).floor() as i64,
            (pos[1] * self.inv_tol).floor() as i64,
            (pos[2] * self.inv_tol).floor() as i64,
        )
    }

    fn exact_key(pos: &[f64; 3]) -> (i64, i64, i64) {
        (pos[0].to_bits() as i64, pos[1].to_bits() as i64, pos[2].to_bits() as i64)
    }

    fn vec3_are_equal(&self, a: &[f64; 3], b: &[f64; 3]) -> bool {
        if self.inv_tol <= 0.0 {
            a[0] == b[0] && a[1] == b[1] && a[2] == b[2]
        } else {
            (a[0] - b[0]).abs() <= self.tolerance
                && (a[1] - b[1]).abs() <= self.tolerance
                && (a[2] - b[2]).abs() <= self.tolerance
        }
    }

    /// Searches the map (including the 26 neighbour cells) for a merge
    /// candidate. Sets `is_opposite` when the existing node has an opposite
    /// normal.
    fn find_merge(&self, pos: &[f64; 3], normal: &[f64; 3], is_opposite: &mut bool) -> Option<usize> {
        let keys: Vec<(i64, i64, i64)> = if self.inv_tol > 0.0 {
            let c = self.cell(pos);
            let mut keys = vec![c];
            for dx in -1i64..=1 {
                for dy in -1i64..=1 {
                    for dz in -1i64..=1 {
                        if dx == 0 && dy == 0 && dz == 0 { continue; }
                        keys.push((c.0 + dx, c.1 + dy, c.2 + dz));
                    }
                }
            }
            keys
        } else {
            vec![Self::exact_key(pos)]
        };
        for key in keys {
            if let Some(bucket) = self.buckets.get(&key) {
                for &mi in bucket {
                    let k = &self.merged[mi];
                    if !self.vec3_are_equal(&k.pos, pos) { continue; }
                    let cos = k.normal[0] * normal[0] + k.normal[1] * normal[1] + k.normal[2] * normal[2];
                    if cos >= self.angle_cos {
                        return Some(mi);
                    } else if self.to_merge_opposite && cos <= -self.angle_cos {
                        *is_opposite = true;
                        return Some(mi);
                    }
                }
            }
        }
        None
    }

    /// Binds the node to the map or finds an existing one.
    /// `index` is an in/out 0-based node index. Returns TRUE if the node was
    /// newly bound.
    fn bind(&mut self, index: &mut usize, is_opposite: &mut bool, pos: [f64; 3], normal: [f64; 3]) -> bool {
        *is_opposite = false;
        if let Some(existing) = self.find_merge(&pos, &normal, is_opposite) {
            *index = existing;
            return false;
        }
        let new_index = self.merged.len();
        self.merged.push(MergedKey { pos, normal });
        let key = if self.inv_tol > 0.0 { self.cell(&pos) } else { Self::exact_key(&pos) };
        self.buckets.entry(key).or_default().push(new_index);
        *index = new_index;
        true
    }
}

/// Auxiliary tool for merging triangulation nodes.
///
/// The tool merges all nodes within the given tolerance, but keeps nodes
/// separated when their incident triangle normals differ by more than the
/// smooth angle.
#[derive(Debug, Clone)]
pub struct PolyMergeNodesTool {
    output: Option<PolyTriangulation>,
    node_index_map: MergedNodesMap,
    elem_map: std::collections::HashSet<[isize; 4]>,
    node_inds: [isize; 4],
    tri_normal: [f64; 3],
    places: [GpXyz; 4],
    unit_factor: f64,
    nb_nodes: usize,
    nb_elems: usize,
    nb_degen_elems: usize,
    nb_merged_elems: usize,
    to_drop_degenerative: bool,
    to_merge_elems: bool,
}

impl PolyMergeNodesTool {
    /// Constructor with the default merge tolerance (`Precision::Confusion`).
    pub fn new(smooth_angle: f64) -> Self {
        Self::with_tolerance(smooth_angle, CONFUSION, 0)
    }

    /// Constructor with explicit tolerance and facet estimate.
    pub fn with_tolerance(smooth_angle: f64, merge_tolerance: f64, nb_facets: usize) -> Self {
        let merge = smooth_angle > 0.0 || merge_tolerance > 0.0;
        let mut map = MergedNodesMap::new();
        map.set_merge_angle(smooth_angle);
        map.set_merge_tolerance(merge_tolerance);
        let _ = nb_facets; // ponytail: preallocation hint not needed for Vec/HashMap growth
        let _ = merge;
        Self {
            output: Some(PolyTriangulation::new()),
            node_index_map: map,
            elem_map: std::collections::HashSet::new(),
            node_inds: [-1; 4],
            tri_normal: [0.0, 0.0, 1.0],
            places: [GpXyz::zero(); 4],
            unit_factor: 1.0,
            nb_nodes: 0,
            nb_elems: 0,
            nb_degen_elems: 0,
            nb_merged_elems: 0,
            to_drop_degenerative: true,
            to_merge_elems: false,
        }
    }

    // ---- configuration ----

    pub fn merge_tolerance(&self) -> f64 { self.node_index_map.tolerance }
    pub fn set_merge_tolerance(&mut self, t: f64) { self.node_index_map.set_merge_tolerance(t); }
    pub fn merge_angle(&self) -> f64 { self.node_index_map.angle }
    pub fn set_merge_angle(&mut self, a: f64) { self.node_index_map.set_merge_angle(a); }
    pub fn to_merge_opposite(&self) -> bool { self.node_index_map.to_merge_opposite }
    pub fn set_merge_opposite(&mut self, b: bool) { self.node_index_map.to_merge_opposite = b; }
    pub fn set_unit_factor(&mut self, f: f64) { self.unit_factor = f; }
    pub fn to_drop_degenerative(&self) -> bool { self.to_drop_degenerative }
    pub fn set_drop_degenerative(&mut self, b: bool) { self.to_drop_degenerative = b; }
    pub fn to_merge_elems(&self) -> bool { self.to_merge_elems }
    pub fn set_merge_elems(&mut self, b: bool) { self.to_merge_elems = b; }

    // ---- element insertion ----

    /// Computes the normalized normal of the current element (the first three
    /// places).
    fn compute_tri_normal(&self) -> [f64; 3] {
        let v01 = self.places[1].subtracted(&self.places[0]);
        let v02 = self.places[2].subtracted(&self.places[0]);
        let c = v01.crossed(&v02);
        let m = c.modulus();
        if m == 0.0 { [0.0, 0.0, 1.0] } else { [c.x / m, c.y / m, c.z / m] }
    }

    /// Pushes a triangle node with merge-by-position-and-normal.
    fn push_node_check(&mut self, is_opposite: &mut bool, tri_node: usize) -> Result<(), String> {
        let mut node_index = self.nb_nodes;
        let place = self.places[tri_node];
        let pos = [place.x, place.y, place.z];
        if self.node_index_map.bind(&mut node_index, is_opposite, pos, self.tri_normal) {
            self.nb_nodes += 1;
            if let Some(out) = self.output.as_mut() {
                if out.nb_nodes() < self.nb_nodes {
                    out.resize_nodes(self.nb_nodes * 2, true)?;
                }
                out.set_node(self.nb_nodes - 1, GpPnt::from_xyz(&place.multiplied(self.unit_factor)));
            }
        }
        self.node_inds[tri_node] = node_index as isize;
        Ok(())
    }

    /// Pushes a triangle node without merging.
    fn push_node_no_merge(&mut self, tri_node: usize) -> Result<(), String> {
        let node_index = self.nb_nodes;
        let place = self.places[tri_node].multiplied(self.unit_factor);
        self.nb_nodes += 1;
        if let Some(out) = self.output.as_mut() {
            if out.nb_nodes() < self.nb_nodes {
                out.resize_nodes(self.nb_nodes * 2, true)?;
            }
            out.set_node(self.nb_nodes - 1, GpPnt::from_xyz(&place));
        }
        self.node_inds[tri_node] = node_index as isize;
        Ok(())
    }

    /// Pushes the element whose node coordinates are set in `places`.
    pub fn push_last_element(&mut self, n: usize) -> Result<(), String> {
        if n != 3 && n != 4 {
            return Err("PolyMergeNodesTool::push_last_element - internal error".into());
        }
        let mut is_opposite = false;
        self.node_inds[3] = -1;
        if self.node_index_map.has_merge_angle() || self.node_index_map.has_merge_tolerance() {
            if !self.node_index_map.to_merge_any_angle() {
                self.tri_normal = self.compute_tri_normal();
            }
            self.push_node_check(&mut is_opposite, 0)?;
            self.push_node_check(&mut is_opposite, 1)?;
            self.push_node_check(&mut is_opposite, 2)?;
            if n == 4 {
                self.push_node_check(&mut is_opposite, 3)?;
            }
        } else {
            self.push_node_no_merge(0)?;
            self.push_node_no_merge(1)?;
            self.push_node_no_merge(2)?;
            if n == 4 {
                self.push_node_no_merge(3)?;
            }
        }

        if self.to_drop_degenerative {
            if self.node_inds[0] == self.node_inds[1]
                || self.node_inds[0] == self.node_inds[2]
                || self.node_inds[1] == self.node_inds[2]
            {
                if n != 4 {
                    self.nb_degen_elems += 1;
                    return Ok(());
                }
            }
        }

        if self.to_merge_elems {
            let mut sorted = [self.node_inds[0], self.node_inds[1], self.node_inds[2], self.node_inds[3]];
            sorted[..n].sort_unstable();
            if !self.elem_map.insert(sorted) {
                self.nb_merged_elems += 1;
                return Ok(());
            }
        }

        self.nb_elems += 1;
        if let Some(out) = self.output.as_mut() {
            if out.nb_triangles() < self.nb_elems {
                out.resize_triangles(self.nb_elems * 2, true)?;
            }
            out.set_triangle(
                self.nb_elems - 1,
                Triangle::new(self.node_inds[0] as usize, self.node_inds[1] as usize, self.node_inds[2] as usize),
            );
            if n == 4 {
                self.nb_elems += 1;
                if out.nb_triangles() < self.nb_elems {
                    out.resize_triangles(self.nb_elems * 2, true)?;
                }
                out.set_triangle(
                    self.nb_elems - 1,
                    Triangle::new(self.node_inds[0] as usize, self.node_inds[2] as usize, self.node_inds[3] as usize),
                );
            }
        }
        Ok(())
    }

    /// Adds a triangle (3 nodes) or quad (4 nodes).
    pub fn add_element(&mut self, elem_nodes: &[GpXyz; 4], n: usize) -> Result<(), String> {
        if n != 3 && n != 4 {
            return Err("PolyMergeNodesTool::add_element - internal error".into());
        }
        self.places[0] = elem_nodes[0];
        self.places[1] = elem_nodes[1];
        self.places[2] = elem_nodes[2];
        if n == 4 {
            self.places[3] = elem_nodes[3];
        }
        self.push_last_element(n)
    }

    pub fn add_triangle(&mut self, a: GpXyz, b: GpXyz, c: GpXyz) -> Result<(), String> {
        self.places = [a, b, c, GpXyz::zero()];
        self.push_last_element(3)
    }

    pub fn add_quad(&mut self, a: GpXyz, b: GpXyz, c: GpXyz, d: GpXyz) -> Result<(), String> {
        self.places = [a, b, c, d];
        self.push_last_element(4)
    }

    /// Changes a node coordinate of the element to be pushed.
    pub fn change_element_node(&mut self, index: usize) -> &mut GpXyz { &mut self.places[index] }

    pub fn push_last_triangle(&mut self) -> Result<(), String> { self.push_last_element(3) }
    pub fn push_last_quad(&mut self) -> Result<(), String> { self.push_last_element(4) }

    /// Returns the current element node index defined by `push_last_element`.
    pub fn element_node_index(&self, index: usize) -> isize { self.node_inds[index] }

    pub fn nb_nodes(&self) -> usize { self.nb_nodes }
    pub fn nb_elements(&self) -> usize { self.nb_elems }
    pub fn nb_degenerative_elems(&self) -> usize { self.nb_degen_elems }
    pub fn nb_merged_elems(&self) -> usize { self.nb_merged_elems }

    /// Adds all triangles of a triangulation, applying a transformation and
    /// optional node reversal.
    pub fn add_triangulation(
        &mut self,
        tris: &PolyTriangulation,
        trsf: &GpTrsf,
        to_reverse: bool,
    ) -> Result<(), String> {
        if tris.nb_triangles() == 0 {
            return Ok(());
        }
        if let Some(out) = self.output.as_mut() {
            if out.nb_nodes() == 0 {
                out.resize_nodes(tris.nb_nodes(), false)?;
                out.resize_triangles(tris.nb_triangles(), false)?;
            }
        }
        for elem_iter in 0..tris.nb_triangles() {
            let elem = tris.triangle(elem_iter);
            let (n0, n1, n2) = if to_reverse {
                (elem.n0, elem.n2, elem.n1)
            } else {
                (elem.n0, elem.n1, elem.n2)
            };
            let mut places = [GpXyz::zero(); 4];
            for (k, &n) in [n0, n1, n2].iter().enumerate() {
                let mut xyz = tris.node(n).coord;
                trsf.transforms_xyz(&mut xyz);
                places[k] = xyz;
            }
            self.places = places;
            self.push_last_element(3)?;
        }
        Ok(())
    }

    /// Prepares and returns the result triangulation, truncating temporary
    /// arrays to the actual result size.
    pub fn result(&mut self) -> Result<PolyTriangulation, String> {
        let mut out = self.output.take().ok_or("PolyMergeNodesTool: no output triangulation")?;
        out.resize_nodes(self.nb_nodes, true)?;
        out.resize_triangles(self.nb_elems, true)?;
        Ok(out)
    }

    /// Static entry point: merges nodes of an existing mesh and returns a new
    /// mesh, or `None` when the input is empty / unchanged (and `to_force` is
    /// false).
    pub fn merge_nodes(
        tris: &PolyTriangulation,
        trsf: &GpTrsf,
        to_reverse: bool,
        smooth_angle: f64,
        merge_tolerance: f64,
        to_force: bool,
    ) -> Result<Option<PolyTriangulation>, String> {
        if tris.nb_nodes() < 3 || tris.nb_triangles() < 1 {
            return Ok(None);
        }
        let tol = if merge_tolerance <= 0.0 { CONFUSION } else { merge_tolerance };
        let mut tool = PolyMergeNodesTool::with_tolerance(smooth_angle, tol, tris.nb_triangles());
        tool.add_triangulation(tris, trsf, to_reverse)?;
        if !to_force && tool.nb_nodes() == tris.nb_nodes() && tool.nb_elements() == tris.nb_triangles() {
            return Ok(None);
        }
        Ok(Some(tool.result()?))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unit_square_duplicated() -> PolyTriangulation {
        // Corners 0-3 and their duplicates 4-7 at the same positions.
        let nodes = vec![
            GpPnt::new(0., 0., 0.), GpPnt::new(1., 0., 0.), GpPnt::new(1., 1., 0.), GpPnt::new(0., 1., 0.),
            GpPnt::new(0., 0., 0.), GpPnt::new(1., 0., 0.), GpPnt::new(1., 1., 0.), GpPnt::new(0., 1., 0.),
        ];
        let tris = vec![
            Triangle::new(0, 1, 2),
            Triangle::new(0, 2, 3),
            Triangle::new(4, 5, 6),
            Triangle::new(4, 6, 7),
        ];
        PolyTriangulation::from_parts(nodes, tris)
    }

    #[test]
    fn merge_duplicates_halves_nodes() {
        let tris = unit_square_duplicated();
        let id = GpTrsf::identity();
        let out = PolyMergeNodesTool::merge_nodes(&tris, &id, false, 0.0, CONFUSION, true)
            .unwrap()
            .expect("merged triangulation");
        // 8 nodes -> 4 unique positions (halved).
        assert_eq!(out.nb_nodes(), 4);
        assert_eq!(out.nb_triangles(), 4);
        // All triangle indices remapped within bounds.
        for t in 0..out.nb_triangles() {
            let tr = out.triangle(t);
            for n in [tr.n0, tr.n1, tr.n2] {
                assert!(n < out.nb_nodes(), "node {n} out of bounds");
            }
        }
        // The four output nodes are the four distinct corners.
        let mut xs: Vec<f64> = (0..out.nb_nodes()).map(|i| out.node(i).x()).collect();
        let mut ys: Vec<f64> = (0..out.nb_nodes()).map(|i| out.node(i).y()).collect();
        xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        ys.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert_eq!(xs, vec![0.0, 0.0, 1.0, 1.0]);
        assert_eq!(ys, vec![0.0, 0.0, 1.0, 1.0]);
    }

    #[test]
    fn no_merge_when_tolerance_and_angle_zero() {
        let tris = unit_square_duplicated();
        let mut tool = PolyMergeNodesTool::with_tolerance(0.0, 0.0, 0);
        tool.add_triangulation(&tris, &GpTrsf::identity(), false).unwrap();
        // No merging configured: every node reference is emitted
        // (4 triangles x 3 vertices).
        assert_eq!(tool.nb_nodes(), 12);
        assert_eq!(tool.nb_elements(), 4);
    }

    #[test]
    fn exact_coincident_nodes_merge() {
        // Two triangles with an exactly duplicated vertex (0 and 3 coincide,
        // same orientation -> same normal).
        let nodes = vec![
            GpPnt::new(0., 0., 0.), GpPnt::new(1., 0., 0.), GpPnt::new(0., 1., 0.),
            GpPnt::new(0., 0., 0.), GpPnt::new(-1., 0., 0.), GpPnt::new(0., -1., 0.),
        ];
        let tris = vec![Triangle::new(0, 1, 2), Triangle::new(3, 4, 5)];
        let tri = PolyTriangulation::from_parts(nodes, tris);
        let mut tool = PolyMergeNodesTool::with_tolerance(0.0, CONFUSION, 0);
        tool.add_triangulation(&tri, &GpTrsf::identity(), false).unwrap();
        // Node 0 and node 3 coincide exactly with matching normals -> merged.
        // Output nodes: origin, (1,0,0), (0,1,0), (-1,0,0), (0,-1,0) = 5.
        assert_eq!(tool.nb_nodes(), 5);
    }

    #[test]
    fn degenerate_element_dropped() {
        let mut tool = PolyMergeNodesTool::with_tolerance(0.0, CONFUSION, 0);
        tool.add_triangle(
            GpXyz::new(0., 0., 0.),
            GpXyz::new(1., 0., 0.),
            GpXyz::new(0., 0., 0.), // coincides with the first node
        ).unwrap();
        assert_eq!(tool.nb_elements(), 0);
        assert_eq!(tool.nb_degenerative_elems(), 1);
        let out = tool.result().unwrap();
        assert_eq!(out.nb_triangles(), 0);
    }

    #[test]
    fn quad_emits_two_triangles() {
        let mut tool = PolyMergeNodesTool::with_tolerance(0.0, CONFUSION, 0);
        tool.add_quad(
            GpXyz::new(0., 0., 0.),
            GpXyz::new(1., 0., 0.),
            GpXyz::new(1., 1., 0.),
            GpXyz::new(0., 1., 0.),
        ).unwrap();
        assert_eq!(tool.nb_elements(), 2);
        let out = tool.result().unwrap();
        assert_eq!(out.nb_triangles(), 2);
        assert_eq!(out.nb_nodes(), 4);
    }

    #[test]
    fn merge_opposite_normals_requires_flag() {
        // Two triangles sharing the origin but with opposite orientation, so
        // the origin vertex has opposite normals.
        let nodes = vec![
            GpPnt::new(0., 0., 0.), GpPnt::new(1., 0., 0.), GpPnt::new(0., 1., 0.),
            GpPnt::new(0., 0., 0.), GpPnt::new(0., -1., 0.), GpPnt::new(-1., 0., 0.),
        ];
        let tris = vec![Triangle::new(0, 1, 2), Triangle::new(3, 4, 5)];
        let tri = PolyTriangulation::from_parts(nodes, tris);

        // Default: opposite normals are NOT merged -> origin duplicated.
        let mut t1 = PolyMergeNodesTool::with_tolerance(0.0, CONFUSION, 0);
        t1.add_triangulation(&tri, &GpTrsf::identity(), false).unwrap();
        assert_eq!(t1.nb_nodes(), 6);

        // With the flag: opposite normals merge -> 5 unique positions.
        let mut t2 = PolyMergeNodesTool::with_tolerance(0.0, CONFUSION, 0);
        t2.set_merge_opposite(true);
        t2.add_triangulation(&tri, &GpTrsf::identity(), false).unwrap();
        assert_eq!(t2.nb_nodes(), 5);
    }

    #[test]
    fn transform_applied() {
        let nodes = vec![GpPnt::new(0., 0., 0.), GpPnt::new(1., 0., 0.), GpPnt::new(0., 1., 0.)];
        let tris = vec![Triangle::new(0, 1, 2)];
        let tri = PolyTriangulation::from_parts(nodes, tris);
        let mut tool = PolyMergeNodesTool::with_tolerance(0.0, CONFUSION, 0);
        tool.set_unit_factor(2.0);
        tool.add_triangulation(&tri, &GpTrsf::identity(), false).unwrap();
        let out = tool.result().unwrap();
        assert!((out.node(1).x() - 2.0).abs() < 1e-12);
    }
}
