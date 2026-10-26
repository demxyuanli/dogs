//! `IntPolyh_Triangle` (`IntPolyh_Triangle.hxx:26-230`, `.cxx:49-659`).
//!
//! A triangle over three nodes and three edges, with the mesh-refinement
//! operations (`MiddleRefinement` / `MultipleMiddleRefinement`) and the
//! deflection measure that drives them.

use occt_core::bnd::BndBox;
use occt_core::precision::CONFUSION;
use occt_geom::Surface;

use super::arrays::{ArrayOfEdges, ArrayOfPoints, ArrayOfTriangles};
use super::edge_couple::IntPolyhEdge;
use super::point::IntPolyhPoint;

/// `MyConfusionPrecision` of `IntPolyh_Triangle.cxx:21` (`10.0e-12`).
const MY_CONFUSION_PRECISION: f64 = 10.0e-12;
/// `SquareMyConfusionPrecision` of `IntPolyh_Triangle.cxx:23` (`10.0e-24`).
const SQUARE_MY_CONFUSION_PRECISION: f64 = 10.0e-24;

/// `IntPolyh_Triangle` (`hxx:26`).
#[derive(Debug, Clone, Copy)]
pub struct IntPolyhTriangle {
    /// `myPoints[3]` (`hxx:221`).
    pub points: [i32; 3],
    /// `myEdges[3]` (`hxx:222`).
    pub edges: [i32; 3],
    /// `myEdgesOrientations[3]` (`hxx:223`).
    pub edges_orientations: [i32; 3],
    /// `myHasIntersection` (`hxx:224`).
    pub has_intersection: bool,
    /// `myIsIntersectionPossible` (`hxx:225`).
    pub is_intersection_possible: bool,
    /// `myIsDegenerated` (`hxx:226`).
    pub is_degenerated: bool,
    /// `myDeflection` (`hxx:227`).
    pub deflection: f64,
    /// `myBox` (`hxx:228`), cached by `BoundingBox` (`cxx:627-641`).
    pub bbox: BndBox,
}

impl Default for IntPolyhTriangle {
    /// `IntPolyh_Triangle()` (`hxx:34-49`).
    fn default() -> Self {
        Self {
            points: [-1; 3],
            edges: [-1; 3],
            edges_orientations: [0; 3],
            has_intersection: false,
            is_intersection_possible: true,
            is_degenerated: false,
            deflection: 0.0,
            bbox: BndBox::new(),
        }
    }
}

impl IntPolyhTriangle {
    /// `IntPolyh_Triangle(thePoint1, thePoint2, thePoint3)` (`hxx:52-66`).
    pub fn new(p1: i32, p2: i32, p3: i32) -> Self {
        Self {
            points: [p1, p2, p3],
            ..Self::default()
        }
    }

    /// `FirstPoint()` / `SecondPoint()` / `ThirdPoint()` (`hxx:85-99`).
    pub fn points3(&self) -> (i32, i32, i32) {
        (self.points[0], self.points[1], self.points[2])
    }

    /// `FirstEdge()` / `SecondEdge()` / `ThirdEdge()` (`hxx:79-92`).
    pub fn first_edge(&self) -> i32 {
        self.edges[0]
    }

    /// `SecondEdge()` (`hxx:84`).
    pub fn second_edge(&self) -> i32 {
        self.edges[1]
    }

    /// `ThirdEdge()` (`hxx:90`).
    pub fn third_edge(&self) -> i32 {
        self.edges[2]
    }

    /// `GetEdgeNumber(theEdgeIndex)` (`hxx:151-155`), 1-based.
    pub fn edge_number(&self, edge_index: i32) -> i32 {
        if (1..=3).contains(&edge_index) {
            self.edges[(edge_index - 1) as usize]
        } else {
            0
        }
    }

    /// `GetEdgeOrientation(theEdgeIndex)` (`hxx:165-169`), 1-based.
    pub fn edge_orientation(&self, edge_index: i32) -> i32 {
        if (1..=3).contains(&edge_index) {
            self.edges_orientations[(edge_index - 1) as usize]
        } else {
            0
        }
    }

    /// `ComputeDeflection(theSurface, TPoints)` (`cxx:56-96`): distance from the
    /// surface's UV barycentre to the triangle plane. A triangle with more than
    /// one degenerated node, or with a null normal, is marked degenerated and
    /// keeps a zero deflection.
    pub fn compute_deflection(&mut self, surface: &dyn Surface, points: &ArrayOfPoints) -> f64 {
        self.deflection = 0.0;

        let p1 = points.get(self.points[0]);
        let p2 = points.get(self.points[1]);
        let p3 = points.get(self.points[2]);

        // `cxx:63-72`: no more than one degenerated point is allowed.
        let i_deg = i32::from(p1.degenerated) + i32::from(p2.degenerated) + i32::from(p3.degenerated);
        if i_deg > 1 {
            self.is_degenerated = true;
            return self.deflection;
        }

        // `cxx:75-84`: plane of the triangle.
        let mut normale = IntPolyhPoint::default();
        normale.cross(&p2.sub(&p1), &p3.sub(&p1));
        let sq_norm = normale.square_modulus();
        if sq_norm < SQUARE_MY_CONFUSION_PRECISION {
            self.is_degenerated = true;
            return self.deflection;
        }

        // `cxx:87-91`: UV barycentre evaluated on the surface.
        let gu = (p1.u + p2.u + p3.u) / 3.0;
        let gv = (p1.v + p2.v + p3.v) / 3.0;
        let pt = surface.value(gu, gv);
        let barycentre = IntPolyhPoint::new(pt.x(), pt.y(), pt.z(), gu, gv);

        // `cxx:93-94`.
        let unit = normale.divide(sq_norm.sqrt());
        self.deflection = unit.dot(&barycentre.sub(&p1)).abs();
        self.deflection
    }

    /// `GetNextTriangle(theTriangle, theEdgeNum, TEdges)` (`cxx:101-114`).
    pub fn next_triangle(&self, triangle: i32, edge_num: i32, edges: &ArrayOfEdges) -> i32 {
        if !(1..=3).contains(&edge_num) {
            return -1;
        }
        let e = edges.get(self.edges[(edge_num - 1) as usize]);
        if e.triangle1 == triangle {
            e.triangle2
        } else {
            e.triangle1
        }
    }

    /// `LinkEdges2Triangle(TEdges, theEdge1, theEdge2, theEdge3)`
    /// (`cxx:117-135`).
    pub fn link_edges_2_triangle(
        &mut self,
        edges: &ArrayOfEdges,
        edge1: i32,
        edge2: i32,
        edge3: i32,
    ) {
        if edge1 < 0 || edge2 < 0 || edge3 < 0 {
            return;
        }
        self.edges = [edge1, edge2, edge3];
        for i in 0..3 {
            let first = edges.get(self.edges[i]).point1;
            self.edges_orientations[i] = if first == self.points[i] { 1 } else { -1 };
        }
    }

    /// `SetEdgeAndOrientation(theEdge, theEdgeIndex)` (`cxx:590-623`): finds
    /// which of the three triangle edges the given edge is, and records its
    /// orientation.
    pub fn set_edge_and_orientation(&mut self, edge: &IntPolyhEdge, edge_index: i32) {
        let pe1 = edge.point1;
        let pe2 = edge.point2;
        for i in 0..3 {
            let i1 = if i + 1 > 2 { 0 } else { i + 1 };
            if pe1 == self.points[i] && pe2 == self.points[i1] {
                self.edges[i] = edge_index;
                self.edges_orientations[i] = 1;
                break;
            }
            if pe1 == self.points[i1] && pe2 == self.points[i] {
                self.edges[i] = edge_index;
                self.edges_orientations[i] = -1;
                break;
            }
        }
    }

    /// `BoundingBox(TArrayOfPoints)` (`cxx:627-641`): the box is cached in
    /// `myBox` and gaped by the deflection plus `Precision::Confusion()`.
    pub fn bounding_box(&mut self, points: &ArrayOfPoints) -> BndBox {
        if self.bbox.is_void() {
            let a = points.get(self.points[0]);
            let b = points.get(self.points[1]);
            let c = points.get(self.points[2]);
            self.bbox.add_point(&a.pnt());
            self.bbox.add_point(&b.pnt());
            self.bbox.add_point(&c.pnt());
            self.bbox.set_gap(self.deflection + CONFUSION);
        }
        self.bbox
    }

    /// `MiddleRefinement(NumTri, MySurface, TPoints, TTriangles, TEdges)`
    /// (`cxx:255-551`): splits along the longest of the three edges, and the
    /// neighbouring triangle across that edge when it exists.
    pub fn middle_refinement(
        &mut self,
        num_tri: i32,
        surface: &dyn Surface,
        points: &mut ArrayOfPoints,
        triangles: &mut ArrayOfTriangles,
        edges: &mut ArrayOfEdges,
    ) {
        let mut fin_te = edges.nb_items();
        let mut fin_tt = triangles.nb_items();

        let num_p1 = self.points[0];
        let num_p2 = self.points[1];
        let num_p3 = self.points[2];

        let p1 = points.get(num_p1);
        let p2 = points.get(num_p2);
        let p3 = points.get(num_p3);

        let l12 = p1.square_distance(&p2);
        let l23 = p2.square_distance(&p3);
        let l31 = p3.square_distance(&p1);

        if l12 > l23 && l12 > l31 {
            // `cxx:280-368`.
            let fin_tp = points.nb_items();
            points.get_mut(fin_tp).middle(surface, &p1, &p2);

            let t1 = fin_tt;
            new_triangle(num_p2, num_p3, fin_tp, triangles, surface, points);
            fin_tt += 1;
            let t2 = fin_tt;
            new_triangle(num_p3, num_p1, fin_tp, triangles, surface, points);
            fin_tt += 1;

            let num_ta = self.next_triangle(num_tri, 1, edges);
            if num_ta >= 0 {
                let (num_p3b, _p3b_index, edge2b, edge3b) =
                    get_info_ta(num_p1, num_p2, num_ta, triangles);

                let t3 = fin_tt;
                new_triangle(num_p2, num_p3b, fin_tp, triangles, surface, points);
                fin_tt += 1;
                let t4 = fin_tt;
                new_triangle(num_p3b, num_p1, fin_tp, triangles, surface, points);

                let e1 = fin_te;
                new_edge(num_p1, fin_tp, t2, t4, edges);
                fin_te += 1;
                let e2 = fin_te;
                new_edge(fin_tp, num_p2, t1, t3, edges);
                fin_te += 1;
                let e3 = fin_te;
                new_edge(fin_tp, num_p3, t1, t2, edges);
                fin_te += 1;
                let e4 = fin_te;
                new_edge(fin_tp, num_p3b, t3, t4, edges);

                old_edge(self.edges[1], num_tri, t1, edges);
                old_edge(self.edges[2], num_tri, t2, edges);
                old_edge(edge2b, num_ta, t3, edges);
                old_edge(edge3b, num_ta, t4, edges);

                triangles
                    .get_mut(t1)
                    .link_edges_2_triangle(edges, self.edges[1], e3, e2);
                triangles
                    .get_mut(t2)
                    .link_edges_2_triangle(edges, self.edges[2], e1, e3);
                triangles
                    .get_mut(t3)
                    .link_edges_2_triangle(edges, edge2b, e4, e2);
                triangles
                    .get_mut(t4)
                    .link_edges_2_triangle(edges, edge3b, e1, e4);

                // `cxx:340-342`: kill the adjacent triangle.
                triangles.get_mut(num_ta).deflection = -1.0;
                triangles.get_mut(num_ta).is_intersection_possible = false;
            } else {
                // `cxx:344-365`: only two new triangles.
                let e1 = fin_te;
                new_edge(num_p1, fin_tp, t2, -1, edges);
                fin_te += 1;
                let e2 = fin_te;
                new_edge(fin_tp, num_p2, t1, -1, edges);
                fin_te += 1;
                let e3 = fin_te;
                new_edge(fin_tp, num_p3, t1, t2, edges);

                old_edge(self.edges[1], num_tri, t1, edges);
                old_edge(self.edges[2], num_tri, t2, edges);

                triangles
                    .get_mut(t1)
                    .link_edges_2_triangle(edges, self.edges[1], e3, e2);
                triangles
                    .get_mut(t2)
                    .link_edges_2_triangle(edges, self.edges[2], e1, e3);
            }
        } else if l23 > l31 && l23 > l12 {
            // `cxx:369-455`.
            let fin_tp = points.nb_items();
            points.get_mut(fin_tp).middle(surface, &p2, &p3);

            let t1 = fin_tt;
            new_triangle(num_p1, num_p2, fin_tp, triangles, surface, points);
            fin_tt += 1;
            let t2 = fin_tt;
            new_triangle(num_p3, num_p1, fin_tp, triangles, surface, points);
            fin_tt += 1;

            let num_ta = self.next_triangle(num_tri, 2, edges);
            if num_ta >= 0 {
                let (num_p1b, _p1b_index, edge1b, edge3b) =
                    get_info_ta(num_p2, num_p3, num_ta, triangles);

                let t3 = fin_tt;
                new_triangle(num_p2, num_p1b, fin_tp, triangles, surface, points);
                fin_tt += 1;
                let t4 = fin_tt;
                new_triangle(num_p1b, num_p3, fin_tp, triangles, surface, points);

                let e1 = fin_te;
                new_edge(num_p2, fin_tp, t1, t3, edges);
                fin_te += 1;
                let e2 = fin_te;
                new_edge(fin_tp, num_p3, t2, t4, edges);
                fin_te += 1;
                let e3 = fin_te;
                new_edge(fin_tp, num_p1, t1, t2, edges);
                fin_te += 1;
                let e4 = fin_te;
                new_edge(fin_tp, num_p1b, t3, t4, edges);

                old_edge(self.edges[0], num_tri, t1, edges);
                old_edge(self.edges[2], num_tri, t2, edges);
                old_edge(edge1b, num_ta, t3, edges);
                old_edge(edge3b, num_ta, t4, edges);

                triangles
                    .get_mut(t1)
                    .link_edges_2_triangle(edges, self.edges[0], e1, e3);
                triangles
                    .get_mut(t2)
                    .link_edges_2_triangle(edges, self.edges[2], e3, e2);
                triangles
                    .get_mut(t3)
                    .link_edges_2_triangle(edges, edge1b, e4, e1);
                triangles
                    .get_mut(t4)
                    .link_edges_2_triangle(edges, edge3b, e2, e4);

                triangles.get_mut(num_ta).deflection = -1.0;
                triangles.get_mut(num_ta).is_intersection_possible = false;
            } else {
                let e1 = fin_te;
                new_edge(num_p2, fin_tp, t1, -1, edges);
                fin_te += 1;
                let e2 = fin_te;
                new_edge(fin_tp, num_p3, t2, -1, edges);
                fin_te += 1;
                let e3 = fin_te;
                new_edge(fin_tp, num_p1, t1, t2, edges);

                old_edge(self.edges[0], num_tri, t1, edges);
                old_edge(self.edges[2], num_tri, t2, edges);

                triangles
                    .get_mut(t1)
                    .link_edges_2_triangle(edges, self.edges[0], e1, e3);
                triangles
                    .get_mut(t2)
                    .link_edges_2_triangle(edges, self.edges[2], e3, e2);
            }
        } else {
            // `cxx:457-546`.
            let fin_tp = points.nb_items();
            points.get_mut(fin_tp).middle(surface, &p3, &p1);

            let t1 = fin_tt;
            new_triangle(num_p1, num_p2, fin_tp, triangles, surface, points);
            fin_tt += 1;
            let t2 = fin_tt;
            new_triangle(num_p2, num_p3, fin_tp, triangles, surface, points);
            fin_tt += 1;

            let num_ta = self.next_triangle(num_tri, 3, edges);
            if num_ta >= 0 {
                let (num_p2b, _p2b_index, edge1b, edge2b) =
                    get_info_ta(num_p3, num_p1, num_ta, triangles);

                let t3 = fin_tt;
                new_triangle(num_p1, num_p2b, fin_tp, triangles, surface, points);
                fin_tt += 1;
                let t4 = fin_tt;
                new_triangle(num_p2b, num_p3, fin_tp, triangles, surface, points);

                let e1 = fin_te;
                new_edge(num_p2, fin_tp, t1, t2, edges);
                fin_te += 1;
                let e2 = fin_te;
                new_edge(fin_tp, num_p3, t2, t4, edges);
                fin_te += 1;
                let e3 = fin_te;
                new_edge(fin_tp, num_p2b, t4, t3, edges);
                fin_te += 1;
                let e4 = fin_te;
                new_edge(fin_tp, num_p1, t1, t3, edges);

                old_edge(self.edges[0], num_tri, t1, edges);
                old_edge(self.edges[1], num_tri, t2, edges);
                old_edge(edge1b, num_ta, t3, edges);
                old_edge(edge2b, num_ta, t4, edges);

                triangles
                    .get_mut(t1)
                    .link_edges_2_triangle(edges, self.edges[0], e1, e4);
                triangles
                    .get_mut(t2)
                    .link_edges_2_triangle(edges, self.edges[1], e2, e1);
                triangles
                    .get_mut(t3)
                    .link_edges_2_triangle(edges, edge1b, e3, e4);
                triangles
                    .get_mut(t4)
                    .link_edges_2_triangle(edges, edge2b, e2, e3);

                triangles.get_mut(num_ta).deflection = -1.0;
                triangles.get_mut(num_ta).is_intersection_possible = false;
            } else {
                let e1 = fin_te;
                new_edge(num_p2, fin_tp, t1, t2, edges);
                fin_te += 1;
                let e2 = fin_te;
                new_edge(fin_tp, num_p3, t2, -1, edges);
                fin_te += 1;
                let e4 = fin_te;
                new_edge(fin_tp, num_p1, t1, -1, edges);

                old_edge(self.edges[0], num_tri, t1, edges);
                old_edge(self.edges[1], num_tri, t2, edges);

                triangles
                    .get_mut(t1)
                    .link_edges_2_triangle(edges, self.edges[0], e1, e4);
                triangles
                    .get_mut(t2)
                    .link_edges_2_triangle(edges, self.edges[1], e2, e1);
            }
        }
        points.increment_nb_items();

        // `cxx:549-550`: make the triangle obsolete.
        self.deflection = -1.0;
        self.is_intersection_possible = false;
    }

    /// `Copy(theOther)` (`cxx:523-533`).
    pub fn copy_from(&mut self, other: &Self) {
        *self = *other;
    }

    /// `MultipleMiddleRefinement(theRefineCriterion, theBox,
    /// theTriangleNumber, ...)` (`cxx:555-583`).
    ///
    /// The loop bound is `2 * FinTTInit + 1000`, so each triangle may be split
    /// at least once while freshly created triangles are still checked.
    #[allow(clippy::too_many_arguments)]
    pub fn multiple_middle_refinement(
        &self,
        refine_criterion: f64,
        bbox: &BndBox,
        triangle_number: i32,
        surface: &dyn Surface,
        points: &mut ArrayOfPoints,
        triangles: &mut ArrayOfTriangles,
        edges: &mut ArrayOfEdges,
    ) {
        let fin_tt_init = triangles.nb_items();
        let max_nb_tt = 2 * fin_tt_init + 1000;

        middle_refinement_at(triangles, triangle_number, surface, points, edges);

        let mut i = fin_tt_init;
        while i < triangles.nb_items() && i < max_nb_tt {
            let mut tri = triangles.get(i);
            let bbox_tri = tri.bounding_box(points);
            if bbox.is_out_box(&bbox_tri) {
                triangles.get_mut(i).is_intersection_possible = false;
            } else if tri.deflection > refine_criterion {
                middle_refinement_at(triangles, i, surface, points, edges);
            }
            i += 1;
        }
    }
}

/// `TTriangles[theIndex].MiddleRefinement(...)`: the `.cxx` mutates the element
/// in place while the callee also rewrites other elements of the same array,
/// which one `&mut` cannot alias. The element is taken out by value and stored
/// back, so `self.deflection = -1.0` (`cxx:549-550`) still lands in the array.
pub fn middle_refinement_at(
    triangles: &mut ArrayOfTriangles,
    index: i32,
    surface: &dyn Surface,
    points: &mut ArrayOfPoints,
    edges: &mut ArrayOfEdges,
) {
    let mut tri = triangles.get(index);
    tri.middle_refinement(index, surface, points, triangles, edges);
    triangles.get_mut(index).copy_from(&tri);
}

/// `TTriangles[theIndex].MultipleMiddleRefinement(...)`. The method is
/// `const` in the `.cxx` and only refines other elements (`cxx:555-583`), so no
/// copy of `theIndex` has to be stored back.
#[allow(clippy::too_many_arguments)]
pub fn multiple_middle_refinement_at(
    triangles: &mut ArrayOfTriangles,
    index: i32,
    refine_criterion: f64,
    bbox: &BndBox,
    surface: &dyn Surface,
    points: &mut ArrayOfPoints,
    edges: &mut ArrayOfEdges,
) {
    let tri = triangles.get(index);
    tri.multiple_middle_refinement(
        refine_criterion,
        bbox,
        index,
        surface,
        points,
        triangles,
        edges,
    );
}

/// `GetInfoTA(numP1, numP2, numTA, TTriangles, numP3b, P3bIndex, Edge2b,
/// Edge3b)` (`cxx:138-205`): finds the third point of the adjacent triangle
/// and the two edges leaving it.
fn get_info_ta(
    num_p1: i32,
    num_p2: i32,
    num_ta: i32,
    triangles: &ArrayOfTriangles,
) -> (i32, i32, i32, i32) {
    let tri = triangles.get(num_ta);
    let p1b = tri.points[0];
    let p2b = tri.points[1];
    let p3b = tri.points[2];

    if p1b != num_p1 && p1b != num_p2 {
        let p3b_index = 1;
        let (edge2b, edge3b) = if p2b == num_p1 {
            (tri.edges[2], tri.edges[0])
        } else {
            (tri.edges[0], tri.edges[2])
        };
        (p1b, p3b_index, edge2b, edge3b)
    } else if p2b != num_p1 && p2b != num_p2 {
        let p3b_index = 2;
        let (edge2b, edge3b) = if p1b == num_p1 {
            (tri.edges[1], tri.edges[0])
        } else {
            (tri.edges[0], tri.edges[1])
        };
        (p2b, p3b_index, edge2b, edge3b)
    } else if p3b != num_p1 && p3b != num_p2 {
        let p3b_index = 3;
        let (edge2b, edge3b) = if p2b == num_p1 {
            (tri.edges[2], tri.edges[1])
        } else {
            (tri.edges[1], tri.edges[2])
        };
        (p3b, p3b_index, edge2b, edge3b)
    } else {
        (-1, -1, -1, -1)
    }
}

/// `NewTriangle(P1, P2, P3, TTriangles, MySurface, TPoints)` (`cxx:208-223`).
fn new_triangle(
    p1: i32,
    p2: i32,
    p3: i32,
    triangles: &mut ArrayOfTriangles,
    surface: &dyn Surface,
    points: &ArrayOfPoints,
) {
    let fin_tt = triangles.nb_items();
    let t = triangles.get_mut(fin_tt);
    t.points = [p1, p2, p3];
    t.compute_deflection(surface, points);
    triangles.increment_nb_items();
}

/// `NewEdge(P1, P2, T1, T2, TEdges)` (`cxx:227-237`).
fn new_edge(p1: i32, p2: i32, t1: i32, t2: i32, edges: &mut ArrayOfEdges) {
    let fin_te = edges.nb_items();
    let e = edges.get_mut(fin_te);
    e.point1 = p1;
    e.point2 = p2;
    e.triangle1 = t1;
    e.triangle2 = t2;
    edges.increment_nb_items();
}

/// `OldEdge(EdgeN, NumTri, NewTriNum, TEdges)` (`cxx:240-251`).
fn old_edge(edge_n: i32, num_tri: i32, new_tri_num: i32, edges: &mut ArrayOfEdges) {
    let e = edges.get_mut(edge_n);
    if e.triangle1 == num_tri {
        e.triangle1 = new_tri_num;
    } else {
        e.triangle2 = new_tri_num;
    }
}

/// `MyTolerance` of `IntPolyh_Triangle.cxx:20` (`10.0e-7`), kept for the
/// constant's sake: no function of this file reads it, but the analogous
/// tolerance is used by `IntPolyh_MaillageAffinage::TriContact`.
#[allow(dead_code)]
pub const MY_TOLERANCE: f64 = 10.0e-7;

/// `MyConfusionPrecision` of `IntPolyh_Triangle.cxx:21`, see
/// [`MY_CONFUSION_PRECISION`].
#[allow(dead_code)]
pub const TRIANGLE_CONFUSION: f64 = MY_CONFUSION_PRECISION;
