//! `BRepExtrema_DistShapeShape` — the minimum distance between two shapes.
//!
//! Source: `ModelingAlgorithms/TKTopAlgo/BRepExtrema/BRepExtrema_DistShapeShape.cxx`
//! (1162 lines) and `.hxx:33-212` (OCCT 8.0.0).
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `(Shape1, Shape2, Deflection, F, A)` | `cxx:544-580` | [`DistShapeShape::new`] |
//! | `LoadS1` / `LoadS2` | `cxx:674-690` | [`DistShapeShape::load_s1`] / [`DistShapeShape::load_s2`] |
//! | `Perform` | `cxx:827-970` | [`DistShapeShape::perform`] |
//! | `SolidTreatment` + `TreatmentFunctor` | `cxx:698-822` | [`DistShapeShape::solid_treatment`] |
//! | `DistanceVertVert` | `cxx:233-460` | [`DistShapeShape::distance_vert_vert`] |
//! | `DistanceMapMap` + `DistancePairFunctor` / `DistanceFunctor` | `cxx:467-670` | [`DistShapeShape::distance_map_map`] |
//! | `Value` / `NbSolution` / `InnerSolution` | `cxx:982-992` | the same names |
//! | `SupportOnShape1/2`, `ParOnEdgeS1/2`, `ParOnFaceS1/2` | `cxx:994-1112` | the same names |
//! | `Decomposition` | `cxx:48-59` | [`decomposition`] |
//!
//! # Translation notes
//!
//! * `Decomposition` is `TopExp::MapShapes(S, TopAbs_VERTEX/EDGE/FACE)`, i.e.
//!   a de-duplicated sub-shape map per type. The port uses
//!   [`crate::topexp::map_shapes`] and keeps the three maps.
//! * `BoxCalculation` is `BRepBndLib::Add` per map entry, here
//!   [`crate::brep_bnd_lib::add`].
//! * The `OSD_Parallel::For` split in `SolidTreatment`, `DistanceVertVert` and
//!   `DistanceMapMap` only distributes the *same* work over `aNbTasks` bands and
//!   then merges the per-task results with the same `-Eps` / `< Eps` rules
//!   (`cxx:250-275`, `cxx:452-460`, `cxx:600-620`). The port runs the bands
//!   serially and merges in the same order, which is what
//!   `OSD_Parallel::For(..., !myIsMultiThread)` does when `myIsMultiThread` is
//!   false (`cxx:818`).
//! * `DistanceMapMap`'s box–box distance is
//!   [`crate::brep_bnd_lib`]'s boxes through
//!   [`crate::brep_extrema::shape_bbox_distance`] (`Bnd_Box::Distance`).
//! * The per-pair refinement is `BRepExtrema_DistanceSS`
//!   ([`crate::brep_extrema_distance_ss::DistanceSS`]).
//! * `myEps` default `Precision::Confusion()`.
//! * `SetFlag` / `SetAlgo` are recorded but, as the header says, "not used in
//!   computation and are obsolete" (`hxx:47-48`, `hxx:154-159`); they are only
//!   handed to `BRepExtrema_DistanceSS`, which applies them to point–face
//!   extrema (`hxx:46-49`).

use occt_core::bnd::BndBox;
use occt_core::gp::GpPnt;
use occt_core::precision::Precision;
use occt_geom::extrema_surf::{ExtremaExtAlgo, ExtremaExtFlag};

use crate::abs::ShapeType;
use crate::brep_bnd_lib;
use crate::brep_class3d::SolidClassifier;
use crate::brep_extrema::shape_bbox_distance;
use crate::brep_extrema_distance_ss::DistanceSS;
use crate::brep_extrema_solution_elem::{SolutionElem, SupportType};
use crate::brep_tool::BRepTool;
use crate::fclass2d::FaceState;
use crate::shape::{TopoShape, Vertex};

/// `Decomposition(S, MapV, MapE, MapF)` (`cxx:48-59`): the de-duplicated
/// vertex / edge / face maps of `s`.
///
/// `TopExp::MapShapes(S, TopAbs_*, Map)` fills an `NCollection_IndexedMap`,
/// i.e. one entry per distinct `TShape`. [`crate::topexp::Explorer`] walks the
/// child tree without de-duplicating, so the maps are built here with an
/// `Arc::as_ptr` key per type (the same key the existing
/// `brep_extrema::collect_vertices` uses).
pub fn decomposition(s: &TopoShape) -> (Vec<TopoShape>, Vec<TopoShape>, Vec<TopoShape>) {
    (map_of(s, ShapeType::Vertex), map_of(s, ShapeType::Edge), map_of(s, ShapeType::Face))
}

fn map_of(s: &TopoShape, target: ShapeType) -> Vec<TopoShape> {
    let mut out = Vec::new();
    let mut seen = std::collections::HashSet::new();
    let mut ex = crate::topexp::Explorer::new(s, target);
    while ex.more() {
        let cur = ex.current().clone();
        ex.next();
        if seen.insert(std::sync::Arc::as_ptr(&cur.tshape) as usize) {
            out.push(cur);
        }
    }
    out
}

/// `BoxCalculation(Map, SBox)` (`cxx:61-68`).
fn box_calculation(map: &[TopoShape]) -> Vec<BndBox> {
    map.iter()
        .map(|s| {
            let mut b = BndBox::new();
            brep_bnd_lib::add(s, &mut b, false);
            b
        })
        .collect()
}

/// Stored state of a solution pair, mirroring the two `NCollection_Sequence`s.
type SolutionPair = (Vec<SolutionElem>, Vec<SolutionElem>, f64);

/// `BRepExtrema_DistShapeShape` (`BRepExtrema_DistShapeShape.hxx:33-212`).
pub struct DistShapeShape {
    /// `myDistRef` (`hxx:186`).
    my_dist_ref: f64,
    /// `myIsDone` (`hxx:187`).
    my_is_done: bool,
    /// `mySolutionsShape1` / `mySolutionsShape2` (`hxx:189-191`).
    my_solutions_shape1: Vec<SolutionElem>,
    my_solutions_shape2: Vec<SolutionElem>,
    /// `myInnerSol` (`hxx:192`).
    my_inner_sol: bool,
    /// `myEps` (`hxx:193`).
    my_eps: f64,
    /// `myShape1` / `myShape2` (`hxx:194-195`).
    my_shape1: Option<TopoShape>,
    my_shape2: Option<TopoShape>,
    /// `myMapV1` / `myMapE1` / `myMapF1` and the `2` twins (`hxx:196-201`).
    my_map_v1: Vec<TopoShape>,
    my_map_e1: Vec<TopoShape>,
    my_map_f1: Vec<TopoShape>,
    my_map_v2: Vec<TopoShape>,
    my_map_e2: Vec<TopoShape>,
    my_map_f2: Vec<TopoShape>,
    /// `myIsInitS1` / `myIsInitS2` (`hxx:202-203`).
    my_is_init_s1: bool,
    my_is_init_s2: bool,
    /// `myFlag` / `myAlgo` (`hxx:204-205`).
    my_flag: ExtremaExtFlag,
    my_algo: ExtremaExtAlgo,
    /// `myBV1` / `myBE1` / `myBF1` and the `2` twins (`hxx:206-211`).
    my_bv1: Vec<BndBox>,
    my_be1: Vec<BndBox>,
    my_bf1: Vec<BndBox>,
    my_bv2: Vec<BndBox>,
    my_be2: Vec<BndBox>,
    my_bf2: Vec<BndBox>,
}

impl Default for DistShapeShape {
    /// `BRepExtrema_DistShapeShape()` (`cxx:596-607`).
    fn default() -> Self {
        Self {
            my_dist_ref: 0.0,
            my_is_done: false,
            my_solutions_shape1: Vec::new(),
            my_solutions_shape2: Vec::new(),
            my_inner_sol: false,
            my_eps: Precision::CONFUSION,
            my_shape1: None,
            my_shape2: None,
            my_map_v1: Vec::new(),
            my_map_e1: Vec::new(),
            my_map_f1: Vec::new(),
            my_map_v2: Vec::new(),
            my_map_e2: Vec::new(),
            my_map_f2: Vec::new(),
            my_is_init_s1: false,
            my_is_init_s2: false,
            my_flag: ExtremaExtFlag::MinMax,
            my_algo: ExtremaExtAlgo::Grad,
            my_bv1: Vec::new(),
            my_be1: Vec::new(),
            my_bf1: Vec::new(),
            my_bv2: Vec::new(),
            my_be2: Vec::new(),
            my_bf2: Vec::new(),
        }
    }
}

impl DistShapeShape {
    /// `BRepExtrema_DistShapeShape(Shape1, Shape2, theDeflection, F, A)`
    /// (`cxx:608-625`): `LoadS1` + `LoadS2` + `Perform`.
    pub fn new(
        shape1: &TopoShape,
        shape2: &TopoShape,
        the_deflection: f64,
        f: ExtremaExtFlag,
        a: ExtremaExtAlgo,
    ) -> Self {
        let mut this = Self {
            my_eps: the_deflection,
            my_flag: f,
            my_algo: a,
            ..Self::default()
        };
        this.load_s1(shape1);
        this.load_s2(shape2);
        this.perform();
        this
    }

    /// `BRepExtrema_DistShapeShape(Shape1, Shape2, F, A)` (`cxx:544-563`): the
    /// default-deflection overload.
    pub fn with_default_deflection(
        shape1: &TopoShape,
        shape2: &TopoShape,
        f: ExtremaExtFlag,
        a: ExtremaExtAlgo,
    ) -> Self {
        Self::new(shape1, shape2, Precision::CONFUSION, f, a)
    }

    /// `SetDeflection(theDeflection)` (`hxx:72`).
    pub fn set_deflection(&mut self, the_deflection: f64) {
        self.my_eps = the_deflection;
    }

    /// `LoadS1(Shape1)` (`cxx:674-680`).
    pub fn load_s1(&mut self, shape1: &TopoShape) {
        self.my_shape1 = Some(shape1.clone());
        self.my_is_init_s1 = false;
        let (v, e, f) = decomposition(shape1);
        self.my_map_v1 = v;
        self.my_map_e1 = e;
        self.my_map_f1 = f;
    }

    /// `LoadS2(Shape2)` (`cxx:683-689`).
    pub fn load_s2(&mut self, shape2: &TopoShape) {
        self.my_shape2 = Some(shape2.clone());
        self.my_is_init_s2 = false;
        let (v, e, f) = decomposition(shape2);
        self.my_map_v2 = v;
        self.my_map_e2 = e;
        self.my_map_f2 = f;
    }

    /// `IsDone()` (`hxx:89`).
    pub fn is_done(&self) -> bool {
        self.my_is_done
    }

    /// `NbSolution()` (`hxx:92`).
    pub fn nb_solution(&self) -> usize {
        self.my_solutions_shape1.len()
    }

    /// `Value()` (`cxx:982-990`): panics when `!IsDone()`, like the
    /// `StdFail_NotDone` the `.cxx` throws.
    pub fn value(&self) -> f64 {
        assert!(
            self.my_is_done,
            "BRepExtrema_DistShapeShape::Value: There's no solution"
        );
        self.my_dist_ref
    }

    /// `InnerSolution()` (`hxx:97`).
    pub fn inner_solution(&self) -> bool {
        self.my_inner_sol
    }

    /// `PointOnShape1(N)` (`hxx:103`).
    pub fn point_on_shape1(&self, n: usize) -> GpPnt {
        *self.my_solutions_shape1[n - 1].point()
    }

    /// `PointOnShape2(N)` (`hxx:106`).
    pub fn point_on_shape2(&self, n: usize) -> GpPnt {
        *self.my_solutions_shape2[n - 1].point()
    }

    /// `SupportTypeShape1(N)` (`hxx:111-114`).
    pub fn support_type_shape1(&self, n: usize) -> SupportType {
        self.my_solutions_shape1[n - 1].support_kind()
    }

    /// `SupportTypeShape2(N)` (`hxx:122-125`).
    pub fn support_type_shape2(&self, n: usize) -> SupportType {
        self.my_solutions_shape2[n - 1].support_kind()
    }

    /// `SupportOnShape1(N)` (`cxx:994-1014`).
    pub fn support_on_shape1(&self, n: usize) -> Option<TopoShape> {
        let s = &self.my_solutions_shape1[n - 1];
        match s.support_kind() {
            SupportType::IsVertex => s.vertex().map(|v| v.0),
            SupportType::IsOnEdge => s.edge().map(|e| e.0),
            SupportType::IsInFace => s.face().map(|f| f.0),
        }
    }

    /// `SupportOnShape2(N)` (`cxx:1016-1036`).
    pub fn support_on_shape2(&self, n: usize) -> Option<TopoShape> {
        let s = &self.my_solutions_shape2[n - 1];
        match s.support_kind() {
            SupportType::IsVertex => s.vertex().map(|v| v.0),
            SupportType::IsOnEdge => s.edge().map(|e| e.0),
            SupportType::IsInFace => s.face().map(|f| f.0),
        }
    }

    /// `ParOnEdgeS1(N, t)` (`cxx:1038-1055`).
    pub fn par_on_edge_s1(&self, n: usize) -> f64 {
        self.my_solutions_shape1[n - 1].edge_parameter()
    }

    /// `ParOnEdgeS2(N, t)` (`cxx:1057-1074`).
    pub fn par_on_edge_s2(&self, n: usize) -> f64 {
        self.my_solutions_shape2[n - 1].edge_parameter()
    }

    /// `ParOnFaceS1(N, u, v)` (`cxx:1076-1093`).
    pub fn par_on_face_s1(&self, n: usize) -> (f64, f64) {
        self.my_solutions_shape1[n - 1].face_parameter()
    }

    /// `ParOnFaceS2(N, u, v)` (`cxx:1095-1112`).
    pub fn par_on_face_s2(&self, n: usize) -> (f64, f64) {
        self.my_solutions_shape2[n - 1].face_parameter()
    }

    /// `Perform()` (`cxx:827-970`).
    pub fn perform(&mut self) -> bool {
        // `cxx:829-832`.
        self.my_is_done = false;
        self.my_inner_sol = false;
        self.my_solutions_shape1.clear();
        self.my_solutions_shape2.clear();

        // `cxx:834-837`.
        let shape1 = match &self.my_shape1 {
            Some(s) => s.clone(),
            None => return false,
        };
        let shape2 = match &self.my_shape2 {
            Some(s) => s.clone(),
            None => return false,
        };

        // `cxx:839-842`.
        let an_is_solid1 = matches!(
            shape1.shape_type(),
            ShapeType::Solid | ShapeType::CompSolid
        );
        let an_is_solid2 = matches!(
            shape2.shape_type(),
            ShapeType::Solid | ShapeType::CompSolid
        );

        // `cxx:845-862`.
        if an_is_solid1 {
            let map_v2 = self.my_map_v2.clone();
            if !self.solid_treatment(&shape1, &map_v2) {
                return false;
            }
        }
        if an_is_solid2 && !self.my_inner_sol {
            let map_v1 = self.my_map_v1.clone();
            if !self.solid_treatment(&shape2, &map_v1) {
                return false;
            }
        }

        // `cxx:864` guards the rest: `if (!myInnerSol)`.
        if !self.my_inner_sol {
            // `cxx:865-887`: the cached boxes of shape 1.
            if !self.my_is_init_s1 {
                self.my_bv1 = box_calculation(&self.my_map_v1);
                self.my_be1 = box_calculation(&self.my_map_e1);
                self.my_bf1 = box_calculation(&self.my_map_f1);
                self.my_is_init_s1 = true;
            }
            // `cxx:889-911`: the cached boxes of shape 2.
            if !self.my_is_init_s2 {
                self.my_bv2 = box_calculation(&self.my_map_v2);
                self.my_be2 = box_calculation(&self.my_map_e2);
                self.my_bf2 = box_calculation(&self.my_map_f2);
                self.my_is_init_s2 = true;
            }

            // `cxx:913-922`: the initial reference distance.
            if !self.my_map_v1.is_empty() && !self.my_map_v2.is_empty() {
                let v1 = Vertex::wrap(self.my_map_v1[0].clone());
                let v2 = Vertex::wrap(self.my_map_v2[0].clone());
                self.my_dist_ref = match (v1, v2) {
                    (Some(a), Some(b)) => {
                        BRepTool::vertex_point(&a).distance(&BRepTool::vertex_point(&b))
                    }
                    _ => 1e30,
                };
            } else {
                self.my_dist_ref = 1e30;
            }

            // `cxx:924-962`: the eight fixed passes, in the `.cxx` order.
            if !self.distance_vert_vert() {
                return false;
            }
            if !self.distance_map_map(0) {
                return false;
            }
            if !self.distance_map_map(1) {
                return false;
            }
            if !self.distance_map_map(2) {
                return false;
            }
            if !self.distance_map_map(3) {
                return false;
            }
            if !self.distance_map_map(4) {
                return false;
            }
            if !self.distance_map_map(5) {
                return false;
            }
            if !self.distance_map_map(6) {
                return false;
            }
            // `cxx:960-965`: the face–face pass only runs when the reference is
            // not already zero.
            if self.my_dist_ref.abs() > self.my_eps && !self.distance_map_map(7) {
                return false;
            }

            // `cxx:967-977`: drop the solutions farther than `myDistRef + myEps`.
            let mut i = 0;
            while i < self.my_solutions_shape1.len() {
                if self.my_solutions_shape1[i].dist() > self.my_dist_ref + self.my_eps {
                    self.my_solutions_shape1.remove(i);
                    self.my_solutions_shape2.remove(i);
                } else {
                    i += 1;
                }
            }
            // `cxx:978`.
            self.my_is_done = !self.my_solutions_shape1.is_empty();
        }

        self.my_is_done
    }

    /// `SolidTreatment(theShape, theVertexMap)` (`cxx:764-822`).
    ///
    /// The `.cxx` splits the vertex map into bands and runs them through
    /// `OSD_Parallel::For`; every band only classifies its own vertices with
    /// `BRepClass3d_SolidClassifier` at tolerance `0.001` (`cxx:718-719`) and
    /// stops at the first `TopAbs_IN` (`cxx:729-746`). Running the bands in
    /// order is the single-threaded form.
    fn solid_treatment(&mut self, the_shape: &TopoShape, the_vertex_map: &[TopoShape]) -> bool {
        // `cxx:717`.
        let a_tolerance = 0.001;
        for s in the_vertex_map {
            // `cxx:733-734`.
            let a_vertex = match Vertex::wrap(s.clone()) {
                Some(v) => v,
                None => continue,
            };
            let a_pnt = BRepTool::vertex_point(&a_vertex);
            // `cxx:735-736`.
            let state = SolidClassifier::classify(the_shape, &a_pnt, a_tolerance);
            // `cxx:737`.
            if state == FaceState::In {
                self.my_inner_sol = true;
                self.my_dist_ref = 0.0;
                self.my_is_done = true;
                // `cxx:741-743`.
                let a_sol = SolutionElem::of_vertex(0.0, &a_pnt, SupportType::IsVertex, &a_vertex);
                self.my_solutions_shape1.push(a_sol.clone());
                self.my_solutions_shape2.push(a_sol);
                break;
            }
        }
        true
    }

    /// `DistanceVertVert(Map1, Map2)` (`cxx:233-460`).
    fn distance_vert_vert(&mut self) -> bool {
        let eps = Precision::CONFUSION;
        let mut start_dist = self.my_dist_ref;
        let mut best_seq1: Vec<SolutionElem> = Vec::new();
        let mut best_seq2: Vec<SolutionElem> = Vec::new();

        for v1s in &self.my_map_v1 {
            let v1 = match Vertex::wrap(v1s.clone()) {
                Some(v) => v,
                None => continue,
            };
            let p1 = BRepTool::vertex_point(&v1);
            for v2s in &self.my_map_v2 {
                let v2 = match Vertex::wrap(v2s.clone()) {
                    Some(v) => v,
                    None => continue,
                };
                let p2 = BRepTool::vertex_point(&v2);
                let dist = p1.distance(&p2);
                // `cxx:256-275`.
                if dist < start_dist - eps {
                    best_seq1.clear();
                    best_seq2.clear();
                    best_seq1.push(SolutionElem::of_vertex(
                        dist,
                        &p1,
                        SupportType::IsVertex,
                        &v1,
                    ));
                    best_seq2.push(SolutionElem::of_vertex(
                        dist,
                        &p2,
                        SupportType::IsVertex,
                        &v2,
                    ));
                    start_dist = dist;
                } else if (dist - start_dist).abs() < eps {
                    best_seq1.push(SolutionElem::of_vertex(
                        dist,
                        &p1,
                        SupportType::IsVertex,
                        &v1,
                    ));
                    best_seq2.push(SolutionElem::of_vertex(
                        dist,
                        &p2,
                        SupportType::IsVertex,
                        &v2,
                    ));
                    if start_dist > dist {
                        start_dist = dist;
                    }
                }
            }
        }

        // `cxx:452-460`: the per-band results are appended to the class
        // sequences and the reference is lowered to the band minimum.
        self.my_solutions_shape1.extend(best_seq1);
        self.my_solutions_shape2.extend(best_seq2);
        if start_dist < self.my_dist_ref - eps {
            self.my_dist_ref = start_dist;
        }
        true
    }

    /// `DistanceMapMap(Map1, Map2, LBox1, LBox2)` (`cxx:467-670`).
    ///
    /// `which` selects the (Map1, Map2, Box1, Box2) tuple of the `.cxx` call
    /// sequence at `cxx:926-962`:
    ///
    /// | `which` | call | `cxx` |
    /// |---|---|---|
    /// | 0 | `V1 x E2` | `cxx:926` |
    /// | 1 | `E1 x V2` | `cxx:929` |
    /// | 2 | `V1 x F2` | `cxx:932` |
    /// | 3 | `F1 x V2` | `cxx:935` |
    /// | 4 | `E1 x E2` | `cxx:938` |
    /// | 5 | `E1 x F2` | `cxx:941` |
    /// | 6 | `F1 x E2` | `cxx:944` |
    /// | 7 | `F1 x F2` | `cxx:960` |
    fn distance_map_map(&mut self, which: usize) -> bool {
        let (map1, map2, boxes1, boxes2): (&Vec<TopoShape>, &Vec<TopoShape>, &Vec<BndBox>, &Vec<BndBox>) =
            match which {
                0 => (&self.my_map_v1, &self.my_map_e2, &self.my_bv1, &self.my_be2),
                1 => (&self.my_map_e1, &self.my_map_v2, &self.my_be1, &self.my_bv2),
                2 => (&self.my_map_v1, &self.my_map_f2, &self.my_bv1, &self.my_bf2),
                3 => (&self.my_map_f1, &self.my_map_v2, &self.my_bf1, &self.my_bv2),
                4 => (&self.my_map_e1, &self.my_map_e2, &self.my_be1, &self.my_be2),
                5 => (&self.my_map_e1, &self.my_map_f2, &self.my_be1, &self.my_bf2),
                6 => (&self.my_map_f1, &self.my_map_e2, &self.my_bf1, &self.my_be2),
                _ => (&self.my_map_f1, &self.my_map_f2, &self.my_bf1, &self.my_bf2),
            };

        // `cxx:477-480`.
        if map1.is_empty() || map2.is_empty() {
            return true;
        }

        // `cxx:487-545` (`DistancePairFunctor`): the box–box distance of every
        // pair, sorted ascending. `cxx:548-554` runs `std::stable_sort`.
        let mut pairs: Vec<(usize, usize, f64)> = Vec::new();
        for (i, b1) in boxes1.iter().enumerate() {
            for (j, b2) in boxes2.iter().enumerate() {
                pairs.push((i, j, shape_bbox_distance(b1, b2)));
            }
        }
        pairs.sort_by(|a, b| a.2.total_cmp(&b.2));

        // `cxx:556-570` (`DistanceFunctor`): refine in ascending box distance
        // with the early-exit rule.
        let mut solution: Option<SolutionPair> = None;
        let mut start_dist = self.my_dist_ref;
        let mut seq1: Vec<SolutionElem> = Vec::new();
        let mut seq2: Vec<SolutionElem> = Vec::new();
        for (i, j, box_dist) in pairs {
            // `cxx:337-340`: `if (aPair.Distance > Solution.Dist + Eps) break;`
            // — with the band's running distance.
            let visible_dist = solution.as_ref().map(|s| s.2).unwrap_or(start_dist);
            if box_dist > visible_dist + self.my_eps {
                break;
            }
            // `cxx:341-348`.
            let a_dist_tool = DistanceSS::new(
                &map1[i],
                &map2[j],
                &boxes1[i],
                &boxes2[j],
                solution.as_ref().map(|s| s.2).unwrap_or(start_dist),
                self.my_eps,
                self.my_flag,
                self.my_algo,
            );
            // `cxx:349-373`.
            if a_dist_tool.is_done() {
                let a_dist = a_dist_tool.dist_value();
                let cur = solution.as_ref().map(|s| s.2).unwrap_or(start_dist);
                if a_dist < cur - self.my_eps {
                    seq1 = a_dist_tool.seq1_value().to_vec();
                    seq2 = a_dist_tool.seq2_value().to_vec();
                    solution = Some((seq1.clone(), seq2.clone(), a_dist));
                } else if (a_dist - cur).abs() < self.my_eps {
                    seq1.extend_from_slice(a_dist_tool.seq1_value());
                    seq2.extend_from_slice(a_dist_tool.seq2_value());
                    if cur > a_dist {
                        solution = Some((seq1.clone(), seq2.clone(), a_dist));
                    } else {
                        solution = Some((seq1.clone(), seq2.clone(), cur));
                    }
                }
            }
        }

        // `cxx:600-620`: merge the band result into the class state.
        if let Some((s1, s2, a_dist)) = solution {
            if a_dist < self.my_dist_ref - self.my_eps {
                self.my_solutions_shape1 = s1;
                self.my_solutions_shape2 = s2;
                self.my_dist_ref = a_dist;
            } else if (a_dist - self.my_dist_ref).abs() < self.my_eps {
                self.my_solutions_shape1.extend(s1);
                self.my_solutions_shape2.extend(s2);
                if self.my_dist_ref > a_dist {
                    self.my_dist_ref = a_dist;
                }
            }
        }
        true
    }
}
