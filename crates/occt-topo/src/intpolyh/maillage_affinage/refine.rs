//! Triangle interference lookup and the deflection-driven refinement,
//! `IntPolyh_MaillageAffinage.cxx:1119-1387` plus the
//! `IntPolyh_BoxBndTree` / `IntPolyh_BoxBndTreeSelector` helpers declared at
//! `cxx:116-179`.

use occt_core::bvh::box_3d::{BvhBox3d, BvhVec3d};
use occt_core::bvh::linear_builder_3d::BvhLinearBuilder3;
use occt_core::bvh::pair_traverse_3d::{select_pair_tree3, PairTraverse3};
use occt_core::bvh::radix_sorter_3d::BvhSet3;
use occt_core::bvh::tree_3d::{BvhTree3d, BVH_MAX_TREE_DEPTH};
use occt_geom::Surface;

use super::super::arrays::{ArrayOfEdges, ArrayOfPoints, ArrayOfTriangles};
use super::super::triangle::{middle_refinement_at, multiple_middle_refinement_at};
use super::MaillageAffinage;

/// `IntPolyh_BoxBndTree` = `BVH_BoxSet<double, 3, int>` (`cxx:119`): the boxes
/// of the triangles with the original triangle index as element.
#[derive(Default)]
struct BoxBndTree {
    /// `myBVHSet->Box(i)` / `myBVHSet->Element(i)` (`BVH_BoxSet.hxx`).
    entries: Vec<(i32, BvhBox3d)>,
    bvh: BvhTree3d,
    dirty: bool,
}

impl BoxBndTree {
    fn new() -> Self {
        Self {
            entries: Vec::new(),
            bvh: BvhTree3d::new(),
            dirty: true,
        }
    }

    /// `SetSize(theSize)` (`BVH_BoxSet.hxx`) - only a capacity reservation.
    fn set_size(&mut self, size: usize) {
        self.entries.reserve(size);
    }

    /// `Add(theIndex, theBox)` (`BVH_BoxSet.hxx`).
    fn add(&mut self, index: i32, box_: BvhBox3d) {
        self.entries.push((index, box_));
        self.dirty = true;
    }

    /// `Size()`.
    fn len(&self) -> usize {
        self.entries.len()
    }

    /// `Element(theIndex)`.
    fn element(&self, index: i32) -> i32 {
        self.entries[index as usize].0
    }

    /// `Box(theIndex)`.
    fn element_box(&self, index: i32) -> BvhBox3d {
        self.entries[index as usize].1
    }

    /// `Build()` (`BVH_BoxSet.hxx`, `BVH_PrimitiveSet::Build`).
    fn build(&mut self) {
        if !self.dirty {
            return;
        }
        let box_all = <Self as BvhSet3>::box_all(self);
        let builder = BvhLinearBuilder3::new(10, BVH_MAX_TREE_DEPTH);
        let mut tree = std::mem::take(&mut self.bvh);
        builder.build(self, &mut tree, &box_all);
        self.bvh = tree;
        self.dirty = false;
    }
}

impl BvhSet3 for BoxBndTree {
    fn size(&self) -> usize {
        self.entries.len()
    }

    fn box_at(&self, index: usize) -> BvhBox3d {
        self.entries[index].1
    }

    fn swap_elements(&mut self, index1: usize, index2: usize) {
        self.entries.swap(index1, index2);
    }
}

/// `IntPolyh_BoxBndTreeSelector` (`cxx:122-179`).
struct BoxBndTreeSelector<'a> {
    tree1: &'a BoxBndTree,
    tree2: &'a BoxBndTree,
    /// `myPairs` (`cxx:178`).
    pairs: Vec<(i32, i32)>,
}

impl PairTraverse3 for BoxBndTreeSelector<'_> {
    /// `RejectNode` (`cxx:151-158`).
    fn reject_node(
        &self,
        c_min1: BvhVec3d,
        c_max1: BvhVec3d,
        c_min2: BvhVec3d,
        c_max2: BvhVec3d,
        _metric: &mut f64,
    ) -> bool {
        BvhBox3d::from_corners(c_min1, c_max1).is_out_corners(c_min2, c_max2)
    }

    /// `Accept(theID1, theID2)` (`cxx:161-169`).
    fn accept(&mut self, index1: i32, index2: i32) -> bool {
        if !self
            .tree1
            .element_box(index1)
            .is_out_box(&self.tree2.element_box(index2))
        {
            self.pairs
                .push((self.tree1.element(index1), self.tree2.element(index2)));
            return true;
        }
        false
    }
}

/// `GetInterferingTriangles(...)` (`cxx:185-247`).
///
/// Returns the pairs grouped by the first triangle, in the order of
/// `NCollection_IndexedDataMap` after `aSelector.Sort()`: keys ascending, and
/// the second indices in ascending order inside each key.
pub(super) fn get_interfering_triangles(
    triangles1: &mut ArrayOfTriangles,
    points1: &ArrayOfPoints,
    triangles2: &mut ArrayOfTriangles,
    points2: &ArrayOfPoints,
) -> Vec<(i32, Vec<i32>)> {
    let mut tree1 = BoxBndTree::new();
    let mut tree2 = BoxBndTree::new();

    // 1. Fill the trees with the boxes of the surfaces triangles
    // (`cxx:200-223`).
    {
        let nb_t = triangles1.nb_items();
        tree1.set_size(nb_t as usize);
        for j in 0..nb_t {
            let tri = triangles1.get(j);
            if !tri.is_intersection_possible || tri.is_degenerated {
                continue;
            }
            let bbox = triangles1.get_mut(j).bounding_box(points1);
            tree1.add(j, BvhBox3d::from_bnd(&bbox));
        }
        if tree1.len() == 0 {
            return Vec::new();
        }
    }
    {
        let nb_t = triangles2.nb_items();
        tree2.set_size(nb_t as usize);
        for j in 0..nb_t {
            let tri = triangles2.get(j);
            if !tri.is_intersection_possible || tri.is_degenerated {
                continue;
            }
            let bbox = triangles2.get_mut(j).bounding_box(points2);
            tree2.add(j, BvhBox3d::from_bnd(&bbox));
        }
        if tree2.len() == 0 {
            return Vec::new();
        }
    }

    tree1.build();
    tree2.build();

    let mut selector = BoxBndTreeSelector {
        tree1: &tree1,
        tree2: &tree2,
        pairs: Vec::new(),
    };
    select_pair_tree3(&tree1.bvh, &tree2.bvh, &mut selector);
    selector.pairs.sort();

    let mut couples: Vec<(i32, Vec<i32>)> = Vec::new();
    for (id1, id2) in selector.pairs {
        match couples.last_mut() {
            Some((key, list)) if *key == id1 => list.push(id2),
            _ => couples.push((id1, vec![id2])),
        }
    }
    couples
}

/// `TrianglesDeflectionsRefinement(...)` (`cxx:1119-1191`).
#[allow(clippy::too_many_arguments)]
fn triangles_deflections_refinement(
    surface1: &dyn Surface,
    triangles1: &mut ArrayOfTriangles,
    edges1: &mut ArrayOfEdges,
    points1: &mut ArrayOfPoints,
    fleche_critique1: f64,
    surface2: &dyn Surface,
    triangles2: &mut ArrayOfTriangles,
    edges2: &mut ArrayOfEdges,
    points2: &mut ArrayOfPoints,
    fleche_critique2: f64,
) {
    // Find the intersecting triangles.
    let interfering = get_interfering_triangles(triangles1, points1, triangles2, points2);
    let couples: std::collections::HashMap<i32, Vec<i32>> = interfering.into_iter().collect();

    // Since the number of the triangles may change during refinement, save the
    // number of triangles before refinement (`cxx:1135-1139`).
    let mut int_s2 = vec![false; triangles2.nb_items().max(0) as usize];
    let fin_tt1 = triangles1.nb_items();
    let fin_tt2 = triangles2.nb_items();

    for i_s1 in 0..fin_tt1 {
        if !triangles1.get(i_s1).is_intersection_possible {
            continue;
        }

        // `cxx:1151-1157`: a triangle with no interfering partner is marked so
        // that it is not used again.
        let list = match couples.get(&i_s1) {
            Some(list) if !list.is_empty() => list,
            _ => {
                triangles1.get_mut(i_s1).is_intersection_possible = false;
                continue;
            }
        };

        let tri1 = triangles1.get(i_s1);
        if tri1.deflection > fleche_critique1 {
            middle_refinement_at(triangles1, i_s1, surface1, points1, edges1);
        }

        for &i_s2 in list {
            if i_s2 < 0 || i_s2 >= fin_tt2 {
                continue;
            }
            let slot = &mut int_s2[i_s2 as usize];
            if !*slot {
                *slot = true;
                let tri2 = triangles2.get(i_s2);
                if tri2.deflection > fleche_critique2 {
                    middle_refinement_at(triangles2, i_s2, surface2, points2, edges2);
                }
            }
        }
    }

    // `cxx:1180-1189`: mark the untouched triangles of S2.
    let nb_marked = int_s2.iter().filter(|m| **m).count();
    if nb_marked < fin_tt2 as usize {
        for i in 0..fin_tt2 {
            if !int_s2[i as usize] {
                triangles2.get_mut(i).is_intersection_possible = false;
            }
        }
    }
}

/// `LargeTrianglesDeflectionsRefinement(...)` (`cxx:1197-1275`).
fn large_triangles_deflections_refinement(
    surface: &dyn Surface,
    triangles: &mut ArrayOfTriangles,
    edges: &mut ArrayOfEdges,
    points: &mut ArrayOfPoints,
    opposite_box: &occt_core::bnd::BndBox,
) {
    // `cxx:1205-1228`: triangles overlapping the other surface's box.
    let mut list: Vec<i32> = Vec::new();
    let nb_t = triangles.nb_items();
    for i in 0..nb_t {
        let tri = triangles.get(i);
        if !tri.is_intersection_possible || tri.is_degenerated {
            continue;
        }
        let a_box = triangles.get_mut(i).bounding_box(points);
        if a_box.is_void() {
            continue;
        }
        if a_box.is_out_box(opposite_box) {
            triangles.get_mut(i).is_intersection_possible = false;
            continue;
        }
        list.push(i);
    }

    if list.is_empty() {
        return;
    }

    // `cxx:1239-1253`: refine criterion from the smallest side of the other box.
    let (x0, y0, z0, x1, y1, z1) = opposite_box
        .get()
        .unwrap_or((0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
    let dx = (x1 - x0).abs();
    let dy = (y1 - y0).abs();
    let mut diag = (z1 - z0).abs();
    let dd = if dx > dy { dy } else { dx };
    if diag > dd {
        diag = dd;
    }
    let diag_ponderation = 0.5;
    let critere_affinage = diag * diag_ponderation;

    for i in list {
        let tri = triangles.get(i);
        if tri.deflection > critere_affinage {
            multiple_middle_refinement_at(
                triangles,
                i,
                critere_affinage,
                opposite_box,
                surface,
                points,
                edges,
            );
        } else {
            middle_refinement_at(triangles, i, surface, points, edges);
        }
    }
}

impl MaillageAffinage<'_> {
    /// `TrianglesDeflectionsRefinementBSB()` (`cxx:1286-1387`).
    pub fn triangles_deflections_refinement_bsb(&mut self) {
        // `cxx:1290-1304`.
        self.compute_deflections(1);
        if self.fleche_min1 > self.fleche_max1 {
            return;
        }
        let fleche_critique1 = self.fleche_min1 * 0.2 + self.fleche_max1 * 0.8;

        // `cxx:1307-1317`.
        self.compute_deflections(2);
        if self.fleche_min2 > self.fleche_max2 {
            return;
        }
        let fleche_critique2 = self.fleche_min2 * 0.2 + self.fleche_max2 * 0.8;

        // `cxx:1322-1340`.
        let (mut x0, mut y0, mut z0, x1, y1, z1) =
            self.box1.get().unwrap_or((0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
        x0 -= x1;
        y0 -= y1;
        z0 -= z1;
        let diag1 = (x0 * x0 + y0 * y0 + z0 * z0)
            / (f64::from(self.nb_samples_u1) * f64::from(self.nb_samples_v1));

        let (mut x0, mut y0, mut z0, x1, y1, z1) =
            self.box2.get().unwrap_or((0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
        x0 -= x1;
        y0 -= y1;
        z0 -= z1;
        let diag2 = (x0 * x0 + y0 * y0 + z0 * z0)
            / (f64::from(self.nb_samples_u2) * f64::from(self.nb_samples_v2));

        // `cxx:1342-1385`.
        if diag1 < diag2 {
            if fleche_critique2 < diag1 {
                triangles_deflections_refinement(
                    self.surface1,
                    &mut self.triangles1,
                    &mut self.edges1,
                    &mut self.points1,
                    fleche_critique1,
                    self.surface2,
                    &mut self.triangles2,
                    &mut self.edges2,
                    &mut self.points2,
                    fleche_critique2,
                );
            } else {
                let box1 = self.box1.clone();
                large_triangles_deflections_refinement(
                    self.surface2,
                    &mut self.triangles2,
                    &mut self.edges2,
                    &mut self.points2,
                    &box1,
                );
            }
        } else if fleche_critique1 < diag2 {
            triangles_deflections_refinement(
                self.surface2,
                &mut self.triangles2,
                &mut self.edges2,
                &mut self.points2,
                fleche_critique2,
                self.surface1,
                &mut self.triangles1,
                &mut self.edges1,
                &mut self.points1,
                fleche_critique1,
            );
        } else {
            let box2 = self.box2.clone();
            large_triangles_deflections_refinement(
                self.surface1,
                &mut self.triangles1,
                &mut self.edges1,
                &mut self.points1,
                &box2,
            );
        }
    }
}
