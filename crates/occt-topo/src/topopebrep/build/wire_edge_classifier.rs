//! `TopOpeBRepBuild_WireEdgeClassifier` — classifies edges and wires of a face.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_WireEdgeClassifier.hxx`,
//! `.cxx:52-527` (OCCT 8.0.0).
//!
//! Shapes are wires, elements are edges (`hxx`). A 2D point is taken on the
//! first edge of the element (`ResetElement`, `cxx:425-462`), and every further
//! edge of the block is compared against the 2D line through two such points by
//! `BRepClass_FacePassiveClassifier` (`cxx:466-519`). `State` is that passive
//! classifier's state (`cxx:523-527`).
//!
//! The subclass overrides `Compare` (`cxx:62-163`): the block/block case first
//! walks the elements, then — for two wires — falls back on
//! `LoopToShape` (`cxx:167-210`) plus `TopOpeBRepTool_ShapeClassifier`
//! (`cxx:141-160`).

use occt_core::gp::{GpDir2d, GpLin2d, GpPnt2d};
use occt_core::precision::{ANGULAR, PCONFUSION};

use crate::abs::{Orientation, ShapeType, State};
use crate::boptools_2d;
use crate::brep_class::{BRepClassEdge, BRepClassFacePassiveClassifier, FClassifier, FaceExplorer};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::fclass2d::FaceState;
use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};
use crate::topo_tools::{edge_vertices, explorer_shapes, is_same};
use crate::topopebrep::tool::{tg_geom_e, ShapeClassifier};

use super::block_builder::BlockBuilder;
use super::composite_classifier::CompositeClassifier;
use super::loop_shape::LoopHandle;

/// `class TopOpeBRepBuild_WireEdgeClassifier` (`hxx`).
pub struct WireEdgeClassifier<'a> {
    /// `myBlockBuilder` of `CompositeClassifier` (`hxx:67`).
    my_block_builder: &'a BlockBuilder,
    /// `myBCEdge` (`hxx`): the face the classification happens on, plus the
    /// edge currently compared.
    my_bc_edge: BRepClassEdge,
    /// `myFPC` (`hxx`).
    my_fpc: BRepClassFacePassiveClassifier,
    /// `myPoint2d` (`hxx`).
    my_point2d: GpPnt2d,
    /// `myFirstCompare` (`hxx`).
    my_first_compare: bool,
    /// `myShape` (`hxx`), the result of the last `LoopToShape`.
    my_shape: Option<TopoShape>,
}

impl<'a> WireEdgeClassifier<'a> {
    /// `TopOpeBRepBuild_WireEdgeClassifier(F, BB)` (`cxx:52-58`).
    pub fn new(f: &TopoShape, bb: &'a BlockBuilder) -> Self {
        let mut my_bc_edge = BRepClassEdge::new();
        if let Some(face) = Face::wrap(f.clone()) {
            my_bc_edge.set_face(face);
        }
        Self {
            my_block_builder: bb,
            my_bc_edge,
            my_fpc: BRepClassFacePassiveClassifier::new(),
            my_point2d: GpPnt2d::new(0.0, 0.0),
            my_first_compare: false,
            my_shape: None,
        }
    }

    /// `LoopToShape(L)` (`cxx:167-210`): the face `myBCEdge.Face()` emptied and
    /// rebuilt on the wire made of the block's edges.
    pub fn loop_to_shape(&mut self, l: &LoopHandle) -> Option<TopoShape> {
        self.my_shape = None;
        let bit = {
            let loop_ref = l.read().expect("poisoned Loop lock");
            *loop_ref.block_iterator()
        };
        let mut bit = bit;
        bit.initialize();
        if !bit.more() {
            return self.my_shape.clone();
        }

        let face = self.my_bc_edge.face()?.clone();
        let builder = TopoBuilder::new();
        let mut edges: Vec<Edge> = Vec::new();
        while bit.more() {
            let Some(shape) = self.my_block_builder.element_of_iterator(&bit) else {
                break;
            };
            let e = match Edge::wrap(shape.clone()) {
                Some(e) => e,
                None => break,
            };
            if !boptools_2d::has_curve_on_surface(&e, &face) {
                // `FC2D_CurveOnSurface(E, F, f, l, tolpc, trim3d)` with
                // `trim3d = true`, then
                // `BB.UpdateEdge(E, C2D, F, max(tolpc, tolE))`.
                let _ = crate::boptools_2d::build_pcurve_for_edge_on_face(&e, &face);
            }
            edges.push(e);
            bit.next();
        }
        if edges.is_empty() {
            return self.my_shape.clone();
        }

        let wire = builder.make_wire(&edges);
        let surface = BRepTool::face_surface(&face)?;
        // `F1.EmptyCopied()` then `BB.Add(F, W)`: the copy keeps `F1`'s
        // surface, location and tolerance and starts with no wire.
        let f = builder.make_face_with_tolerance(surface, &[wire], BRepTool::face_tolerance(&face));
        self.my_shape = Some(f.0);
        self.my_shape.clone()
    }

    /// `CompareElementToShape(EE, B)` (`cxx:389-402`).
    pub fn compare_element_to_shape(&mut self, ee: &TopoShape, b: &TopoShape) -> State {
        self.reset_element(ee);
        for e in explorer_shapes(b, ShapeType::Edge) {
            self.compare_element(&e);
        }
        self.state()
    }

    /// `ResetShape(B)` (`cxx:406-421`).
    pub fn reset_shape(&mut self, b: &TopoShape) {
        if b.shape_type() == ShapeType::Edge {
            self.reset_element(b);
        } else if let Some(e) = explorer_shapes(b, ShapeType::Edge).first() {
            let e = e.clone();
            self.reset_element(&e);
        }
    }

    /// `ResetElement(EE)` (`cxx:425-462`).
    pub fn reset_element(&mut self, ee: &TopoShape) {
        let Some(e) = Edge::wrap(ee.clone()) else {
            return;
        };
        let Some(f) = self.my_bc_edge.face().cloned() else {
            return;
        };

        if !boptools_2d::has_curve_on_surface(&e, &f) {
            // `FC2D_CurveOnSurface(E, F, f2, l2, tolpc, trim3d = true)` then
            // `BB.UpdateEdge(E, C2D, F, max(tolE, tolpc))`.
            let _ = boptools_2d::build_pcurve_for_edge_on_face(&e, &f);
        }

        // `FC2D_CurveOnSurface(E, F, f2, l2, tolpc)`: a null curve raises
        // `Standard_ProgramError("WEC : ResetElement")` (`cxx:444-447`).
        let Some((c2d, f2, l2)) = boptools_2d::curve_on_surface_range(&e, &f) else {
            panic!("WEC : ResetElement");
        };

        let t = 0.397891143689;
        let par = (1.0 - t) * f2 + t * l2;
        self.my_point2d = c2d.value(par);

        self.my_first_compare = true;
    }

    /// `CompareElement(EE)` (`cxx:466-519`).
    ///
    /// `CompareElement` always answers `true` (`cxx:468`, `:518`).
    pub fn compare_element(&mut self, ee: &TopoShape) -> bool {
        let b_ret = true;
        let Some(e) = Edge::wrap(ee.clone()) else {
            return b_ret;
        };
        let Some(f) = self.my_bc_edge.face().cloned() else {
            return b_ret;
        };

        if !boptools_2d::has_curve_on_surface(&e, &f) {
            // `FC2D_CurveOnSurface(E, F, f2, l2, tolpc, trim3d = true)` then
            // `BB.UpdateEdge(E, C2D, F, max(tolE, tolpc))`.
            let _ = boptools_2d::build_pcurve_for_edge_on_face(&e, &f);
        }

        if self.my_first_compare {
            let Some((c2d, f2, l2)) = boptools_2d::curve_on_surface_range(&e, &f) else {
                return b_ret;
            };
            let t = 0.33334567;
            let par = (1.0 - t) * f2 + t * l2;
            let p2d = c2d.value(par);

            // `gp_Lin2d(myPoint2d, gp_Vec2d(myPoint2d, p2d))`: `gp_Dir2d`
            // raises on a null vector, as `gp_Lin2d` does here.
            let dx = p2d.x() - self.my_point2d.x();
            let dy = p2d.y() - self.my_point2d.y();
            let dir = GpDir2d::new(dx, dy).expect("gp_Dir2d : null vector in WEC::CompareElement");
            let l2d = GpLin2d::from_pnt_dir(self.my_point2d, dir);
            let dist = self.my_point2d.distance(&p2d);
            let tol2d = PCONFUSION;
            self.my_fpc.reset(&l2d, dist, tol2d);
            self.my_first_compare = false;
        }

        self.my_bc_edge.set_edge(e.clone());
        self.my_fpc.compare(&self.my_bc_edge, e.orientation());
        b_ret
    }

    /// `State()` (`cxx:523-527`).
    pub fn state(&self) -> State {
        to_state(self.my_fpc.state())
    }
}

impl<'a> CompositeClassifier for WireEdgeClassifier<'a> {
    fn block_builder(&self) -> &BlockBuilder {
        self.my_block_builder
    }

    /// `Compare(L1, L2)` (`cxx:62-163`), the override.
    fn compare(&mut self, l1: &LoopHandle, l2: &LoopHandle) -> State {
        let (isshape1, isshape2, shape1, shape2, bit1, bit2) = {
            let a = l1.read().expect("poisoned Loop lock");
            let b = l2.read().expect("poisoned Loop lock");
            (
                a.is_shape(),
                b.is_shape(),
                a.shape().cloned(),
                b.shape().cloned(),
                *a.block_iterator(),
                *b.block_iterator(),
            )
        };

        let mut state = State::Unknown;

        if isshape2 && isshape1 {
            // L1 is Shape, L2 is Shape
            if let (Some(s1), Some(s2)) = (shape1.as_ref(), shape2.as_ref()) {
                state = self.compare_shapes(s1, s2);
            }
        } else if isshape2 && !isshape1 {
            // L1 is Block, L2 is Shape
            let mut bit1 = bit1;
            bit1.initialize();
            let mut yena1 = bit1.more();
            while yena1 {
                let s1 = match self.block_builder().element_of_iterator(&bit1) {
                    Some(s) => s.clone(),
                    None => break,
                };
                let s2 = match shape2.as_ref() {
                    Some(s) => s.clone(),
                    None => break,
                };
                state = self.compare_element_to_shape(&s1, &s2);
                yena1 = false;
                if state == State::Unknown {
                    if bit1.more() {
                        bit1.next();
                    }
                    yena1 = bit1.more();
                }
            }
        } else if !isshape2 && isshape1 {
            // L1 is Shape, L2 is Block
            let s1 = match shape1.as_ref() {
                Some(s) => s.clone(),
                None => return state,
            };
            self.reset_shape(&s1);
            let mut bit2 = bit2;
            bit2.initialize();
            while bit2.more() {
                let s2 = match self.block_builder().element_of_iterator(&bit2) {
                    Some(s) => s.clone(),
                    None => break,
                };
                self.compare_element(&s2);
                bit2.next();
            }
            state = self.state();
        } else if !isshape2 && !isshape1 {
            // L1 is Block, L2 is Block
            if state == State::Unknown {
                let mut bit1 = bit1;
                bit1.initialize();
                let mut yena1 = bit1.more();
                while yena1 {
                    let s1 = match self.block_builder().element_of_iterator(&bit1) {
                        Some(s) => s.clone(),
                        None => break,
                    };
                    self.reset_element(&s1);
                    let mut bit2 = bit2;
                    bit2.initialize();
                    while bit2.more() {
                        let s2 = match self.block_builder().element_of_iterator(&bit2) {
                            Some(s) => s.clone(),
                            None => break,
                        };
                        self.compare_element(&s2);
                        bit2.next();
                    }
                    state = self.state();
                    yena1 = false;
                    if state == State::Unknown {
                        if bit1.more() {
                            bit1.next();
                        }
                        yena1 = bit1.more();
                    }
                }
            }

            if state == State::Unknown {
                let Some(s1) = self.loop_to_shape(l1) else {
                    return state;
                };
                if s1.is_null() {
                    return state;
                }
                let Some(s2) = self.loop_to_shape(l2) else {
                    return state;
                };
                if s2.is_null() {
                    return state;
                }
                let mut sc = ShapeClassifier::with_reference(&s2);
                let samedomain = sc.same_domain();
                sc.set_same_domain(1);
                sc.set_reference(&s2);
                state = sc.state_shape_reference(&s1, &s2).unwrap_or(State::Unknown);
                sc.set_same_domain(samedomain);
            }
        }

        state
    }

    /// `CompareShapes(B1, B2)` (`cxx:232-385`).
    fn compare_shapes(&mut self, b1: &TopoShape, b2: &TopoShape) -> State {
        // evolution xpu101198 (`cxx:235-242`): 1 wire is classified against the
        // connected wires; the shapes are assumed to come from the same shape
        // and, when they do not touch, `WEC` is not reached.
        let mut state = State::Unknown;
        let edges1 = explorer_shapes(b1, ShapeType::Edge);
        if edges1.is_empty() {
            return state;
        }
        let edges2 = explorer_shapes(b2, ShapeType::Edge);

        for e1_shape in &edges1 {
            let Some(e1) = Edge::wrap(e1_shape.clone()) else {
                continue;
            };
            // `TopExp::Vertices(e1, vf1, vl1)` (`cxx:258`).
            let (vf1, vl1) = edge_vertices(&e1);
            let e1clo = match (vf1.as_ref(), vl1.as_ref()) {
                (Some(a), Some(b)) => is_same(&a.0, &b.0),
                _ => false,
            };
            let mapv1 = [vf1.clone(), vl1.clone()];

            self.reset_shape(e1_shape);
            let mut indy = false;
            for e_shape in &edges2 {
                let Some(e) = Edge::wrap(e_shape.clone()) else {
                    continue;
                };
                if is_same(&e.0, &e1.0) {
                    // eap occ416
                    state = State::Unknown;
                    break;
                }
                let (vf, vl) = edge_vertices(&e);
                let eclo = match (vf.as_ref(), vl.as_ref()) {
                    (Some(a), Some(b)) => is_same(&a.0, &b.0),
                    _ => false,
                };
                let contains = |v: &Option<Vertex>, list: &[Option<Vertex>]| -> bool {
                    match v {
                        Some(v) => list.iter().flatten().any(|o| is_same(&o.0, &v.0)),
                        None => false,
                    }
                };
                let hasf = contains(&vf, &mapv1);
                let hasl = contains(&vl, &mapv1);
                let filter = (hasf || hasl) && !e1clo && !eclo;
                if filter {
                    let vshared = if hasl { vl.clone() } else { vf.clone() };
                    if let (Some(vshared), Some(vl1), Some(vl)) =
                        (vshared.as_ref(), vl1.as_ref(), vl.as_ref())
                    {
                        let tg1 = fun_tg_ine(vshared, vl1, &e1);
                        let tg = fun_tg_ine(vshared, vl, &e);
                        let dot = tg1.dot(&tg);
                        // nyixpu (`cxx:296`)
                        let tol = ANGULAR * 1.0e4;
                        if (1.0 + dot).abs() < tol {
                            indy = true;
                        }
                    }
                }
                if indy {
                    state = State::Unknown;
                    break;
                }
                self.compare_element(e_shape);
                state = self.state();
            }
            if state != State::Unknown {
                break;
            }
        }

        let mut resta = state == State::Unknown;
        resta = resta
            && b2.shape_type() == ShapeType::Wire
            && b1.shape_type() == ShapeType::Wire;
        if resta {
            let mape1 = explorer_shapes(b1, ShapeType::Edge);
            // recall: avoid auto-intersection wires (ie B1 and B2 are disjoint)
            for e2_shape in &edges2 {
                let Some(e2) = Edge::wrap(e2_shape.clone()) else {
                    continue;
                };
                if mape1.iter().any(|s| is_same(s, &e2.0)) {
                    continue;
                }
                let Some(the_face) = self.my_bc_edge.face().cloned() else {
                    continue;
                };
                let builder = TopoBuilder::new();
                let ftmp = Face(the_face.0.oriented(Orientation::Forward));
                let Some(surface) = BRepTool::face_surface(&ftmp) else {
                    continue;
                };
                let tol_f1 = BRepTool::face_tolerance(&ftmp);

                // `F2 = ftmp.EmptyCopied(); BB.Add(F2, TopoDS::Wire(B2))`
                let Some(w2) = Wire::wrap(b2.clone()) else {
                    continue;
                };
                let f2 = builder.make_face_with_tolerance(surface.clone(), &[w2], tol_f1);

                // p2d on E2 of B2, E2 not shared by B1
                let Some((c2d, f, l)) = boptools_2d::curve_on_surface_range(&e2, &f2) else {
                    continue;
                };
                let x = 0.45678;
                let p2 = (1.0 - x) * l + x * f;
                let p2d = c2d.value(p2);

                // `F1 = ftmp.EmptyCopied(); BB.Add(F1, TopoDS::Wire(B1))`
                let Some(w1) = Wire::wrap(b1.clone()) else {
                    continue;
                };
                let f1 = builder.make_face_with_tolerance(surface, &[w1], tol_f1);

                let tol_f1 = BRepTool::face_tolerance(&f1);
                let mut fexp = FaceExplorer::new(&f1);
                let fclass = FClassifier::with_args(&mut fexp, &p2d, tol_f1);
                state = to_state(fclass.state());
                return state;
            }
        }

        state
    }

    /// `ResetShape(B)` / `ResetElement(E)` / `CompareElement(E)` / `State()`.
    fn reset_shape(&mut self, b: &TopoShape) {
        WireEdgeClassifier::reset_shape(self, b);
    }

    /// `CompareElementToShape(EE, B)` (`cxx:389-402`).
    fn compare_element_to_shape(&mut self, ee: &TopoShape, b: &TopoShape) -> State {
        WireEdgeClassifier::compare_element_to_shape(self, ee, b)
    }

    fn reset_element(&mut self, e: &TopoShape) {
        WireEdgeClassifier::reset_element(self, e);
    }

    fn compare_element(&mut self, e: &TopoShape) -> bool {
        WireEdgeClassifier::compare_element(self, e)
    }

    fn state(&mut self) -> State {
        WireEdgeClassifier::state(self)
    }
}

/// `FUN_tgINE(v, vl, e)` (`cxx:212-228`): the tangent of `e` at `v`, oriented
/// inside the edge, reversed when `v` is the last vertex `vl`.
fn fun_tg_ine(v: &TopoShape, vl: &TopoShape, e: &Edge) -> occt_core::gp::GpVec {
    let par = BRepTool::parameter_of_vertex(&e.0, v);
    let Some(mut tg) = tg_geom_e(par, &e.0) else {
        // NYIRAISE
        return occt_core::gp::GpVec::new(0.0, 0.0, 0.0);
    };
    if is_same(v, vl) {
        tg.reverse();
    }
    tg
}

/// `TopAbs_State` of a `FaceState`.
fn to_state(s: FaceState) -> State {
    match s {
        FaceState::In => State::In,
        FaceState::Out => State::Out,
        FaceState::On => State::On,
        FaceState::Unknown => State::Unknown,
    }
}
