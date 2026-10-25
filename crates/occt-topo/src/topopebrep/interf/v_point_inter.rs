//! `TopOpeBRep_VPointInter` — a vertex of a 3D intersection line as seen by
//! `TopOpeBRep`, carrying the states of the point inside the two domains.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRep/TopOpeBRep_VPointInter.hxx:32-194`,
//! `.lxx:21-303`, `.cxx:36-483`.
//!
//! Translation notes:
//! - `TopOpeBRep_PThePointOfIntersection myPPOI` is a `IntPatch_Point*`
//!   (`TopOpeBRep_PThePointOfIntersection.hxx:21`); it is `Option<Arc<IntPatchPoint>>`
//!   here. OCCT dereferences it unconditionally, so every accessor answers
//!   `Err` while it is unset instead of trapping.
//! - `ArcOnS1()` / `ArcOnS2()` (`cxx:63-77`) down-cast the `Adaptor2d_Curve2d`
//!   handle to `BRepAdaptor_Curve2d` and read `Edge()`; the port stores that
//!   adaptor directly (see `IntPatchPoint`), so it reads `edge()`.
//! - `VertexOnS1()` / `VertexOnS2()` (`cxx:82-105`) down-cast to
//!   `BRepTopAdaptor_HVertex` and read `Vertex()`; the port stores `HVertex`.
//! - `myNullShape` (`hxx:189`) is a default-constructed `TopoDS_Shape` returned
//!   by the "not on that domain" arms of `Edge(I)` / `Vertex(I)` (`cxx:208`,
//!   `:276`). The port has no null `TopoShape` (see `TopoShape::is_null`), so
//!   those two accessors answer `Option<&TopoShape>` and drop the member.
//! - `Dump(I, F, OS)` (`cxx:380-417`) and `Dump(F1, F2, OS)` (`cxx:420-478`)
//!   write nothing outside `OCCT_DEBUG`; the port has the two no-op
//!   counterparts.
//! - `TopOpeBRep_FFTransitionTool::ProcessLineTransition` (`cxx:409`) is used by
//!   the debug dump only and is not called here.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d};
use occt_core::precision::Precision;

use crate::abs::{ShapeType, State};
use crate::intpatch_point::IntPatchPoint;
use crate::intsurf_transition::IntSurfTransition;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::is_same;

/// `class TopOpeBRep_VPointInter` (`hxx:32-194`).
#[derive(Clone)]
pub struct VPointInter {
    /// `myPPOI` (`hxx:180`): the `IntPatch_Point` this vertex was built from.
    ppoi: Option<Arc<IntPatchPoint>>,
    /// `myShapeIndex` (`hxx:181`): 0 none, 1 S1, 2 S2, 3 both.
    shape_index: i32,
    /// `myState1` (`hxx:182`).
    state1: State,
    /// `myState2` (`hxx:183`).
    state2: State,
    /// `myKeep` (`hxx:184`).
    keep: bool,
    /// `myEdgeON1` (`hxx:185`).
    edge_on1: TopoShape,
    /// `myEdgeON2` (`hxx:186`).
    edge_on2: TopoShape,
    /// `myEdgeONPar1` (`hxx:187`).
    edge_on_par1: f64,
    /// `myEdgeONPar2` (`hxx:188`).
    edge_on_par2: f64,
    /// `myIndex` (`hxx:188`).
    index: i32,
    /// `myS1` (`hxx:190`).
    s1: i32,
    /// `myS2` (`hxx:191`).
    s2: i32,
}

impl Default for VPointInter {
    /// `TopOpeBRep_VPointInter()` (`lxx:22-37`).
    fn default() -> Self {
        Self {
            ppoi: None,
            shape_index: 0,
            state1: State::Unknown,
            state2: State::Unknown,
            keep: false,
            edge_on1: TopoShape::new(ShapeType::Shape),
            edge_on2: TopoShape::new(ShapeType::Shape),
            edge_on_par1: 0.0,
            edge_on_par2: 0.0,
            index: 0,
            s1: 0,
            s2: 0,
        }
    }
}

fn err(name: &str) -> String {
    format!("TopOpeBRep_VPointInter::{name}")
}

impl VPointInter {
    /// `TopOpeBRep_VPointInter()` (`lxx:22-37`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `SetPoint(P)` (`cxx:36-59`): stores the point and derives `myShapeIndex`
    /// from which domains the point lies on.
    pub fn set_point(&mut self, p: Arc<IntPatchPoint>) {
        let is_on1 = p.is_on_dom_s1();
        let is_on2 = p.is_on_dom_s2();
        self.shape_index = if is_on1 && is_on2 {
            3
        } else if is_on2 {
            2
        } else if is_on1 {
            1
        } else {
            0
        };
        self.ppoi = Some(p);
    }

    /// The stored `IntPatch_Point`.
    fn ppoi(&self) -> Result<&IntPatchPoint, String> {
        self.ppoi.as_deref().ok_or_else(|| err("myPPOI"))
    }

    /// `SetShapes(I1, I2)` (`lxx:41-45`).
    pub fn set_shapes(&mut self, i1: i32, i2: i32) {
        self.s1 = i1;
        self.s2 = i2;
    }

    /// `GetShapes(I1, I2)` (`lxx:49-53`).
    pub fn get_shapes(&self) -> (i32, i32) {
        (self.s1, self.s2)
    }

    /// `TransitionOnS1()` (`lxx:57-60`).
    pub fn transition_on_s1(&self) -> Result<IntSurfTransition, String> {
        Ok(*self.ppoi()?.transition_on_s1()?)
    }

    /// `TransitionOnS2()` (`lxx:64-67`).
    pub fn transition_on_s2(&self) -> Result<IntSurfTransition, String> {
        Ok(*self.ppoi()?.transition_on_s2()?)
    }

    /// `TransitionLineArc1()` (`lxx:71-74`).
    pub fn transition_line_arc1(&self) -> Result<IntSurfTransition, String> {
        Ok(*self.ppoi()?.transition_line_arc1()?)
    }

    /// `TransitionLineArc2()` (`lxx:78-81`).
    pub fn transition_line_arc2(&self) -> Result<IntSurfTransition, String> {
        Ok(*self.ppoi()?.transition_line_arc2()?)
    }

    /// `IsOnDomS1()` (`lxx:85-88`).
    pub fn is_on_dom_s1(&self) -> Result<bool, String> {
        Ok(self.ppoi()?.is_on_dom_s1())
    }

    /// `IsOnDomS2()` (`lxx:92-95`).
    pub fn is_on_dom_s2(&self) -> Result<bool, String> {
        Ok(self.ppoi()?.is_on_dom_s2())
    }

    /// `ParametersOnS1(u, v)` (`lxx:99-102`).
    pub fn parameters_on_s1(&self) -> Result<(f64, f64), String> {
        Ok(self.ppoi()?.parameters_on_s1())
    }

    /// `ParametersOnS2(u, v)` (`lxx:106-109`).
    pub fn parameters_on_s2(&self) -> Result<(f64, f64), String> {
        Ok(self.ppoi()?.parameters_on_s2())
    }

    /// `Value()` (`lxx:113-116`).
    pub fn value(&self) -> Result<&GpPnt, String> {
        Ok(self.ppoi()?.value())
    }

    /// `Tolerance()` (`lxx:120-123`).
    pub fn tolerance(&self) -> Result<f64, String> {
        Ok(self.ppoi()?.tolerance())
    }

    /// `ArcOnS1()` (`cxx:63-68`): `BRepAdaptor_Curve2d::Edge()` of the arc.
    pub fn arc_on_s1(&self) -> Result<&TopoShape, String> {
        Ok(&self.ppoi()?.arc_on_s1()?.edge().0)
    }

    /// `ArcOnS2()` (`cxx:72-77`).
    pub fn arc_on_s2(&self) -> Result<&TopoShape, String> {
        Ok(&self.ppoi()?.arc_on_s2()?.edge().0)
    }

    /// `ParameterOnLine()` (`lxx:127-130`).
    pub fn parameter_on_line(&self) -> Result<f64, String> {
        Ok(self.ppoi()?.parameter_on_line())
    }

    /// `ParameterOnArc1()` (`lxx:134-137`).
    pub fn parameter_on_arc1(&self) -> Result<f64, String> {
        self.ppoi()?.parameter_on_arc1()
    }

    /// `IsVertexOnS1()` (`lxx:141-144`).
    pub fn is_vertex_on_s1(&self) -> Result<bool, String> {
        Ok(self.ppoi()?.is_vertex_on_s1())
    }

    /// `VertexOnS1()` (`cxx:82-91`): `BRepTopAdaptor_HVertex::Vertex()`.
    pub fn vertex_on_s1(&self) -> Result<&TopoShape, String> {
        let p = self.ppoi()?;
        if !p.is_vertex_on_s1() {
            return Err("TopOpeBRep_VPointInter::VertexOnS1".to_string());
        }
        Ok(&p.vertex_on_s1()?.vertex().0)
    }

    /// `ParameterOnArc2()` (`lxx:153-156`).
    pub fn parameter_on_arc2(&self) -> Result<f64, String> {
        self.ppoi()?.parameter_on_arc2()
    }

    /// `IsVertexOnS2()` (`lxx:160-163`).
    pub fn is_vertex_on_s2(&self) -> Result<bool, String> {
        Ok(self.ppoi()?.is_vertex_on_s2())
    }

    /// `VertexOnS2()` (`cxx:95-105`).
    pub fn vertex_on_s2(&self) -> Result<&TopoShape, String> {
        let p = self.ppoi()?;
        if !p.is_vertex_on_s2() {
            return Err("TopOpeBRep_VPointInter::VertexOnS2".to_string());
        }
        Ok(&p.vertex_on_s2()?.vertex().0)
    }

    /// `IsInternal()` (`lxx:178-181`): always `false`.
    pub fn is_internal(&self) -> bool {
        false
    }

    /// `IsMultiple()` (`lxx:166-169`).
    pub fn is_multiple(&self) -> Result<bool, String> {
        Ok(self.ppoi()?.is_multiple())
    }

    /// `State(I)` (`cxx:122-144`).
    pub fn state(&self, i: i32) -> Result<State, String> {
        match i {
            1 => Ok(self.state1),
            2 => Ok(self.state2),
            _ => Err("TopOpeBRep_VPointInter::State".to_string()),
        }
    }

    /// `State(S, I)` (`cxx:108-121`): sets the state and latches `myKeep`.
    pub fn set_state(&mut self, s: State, i: i32) -> Result<(), String> {
        match i {
            1 => self.state1 = s,
            2 => self.state2 = s,
            _ => return Err("TopOpeBRep_VPointInter::State".to_string()),
        }
        self.update_keep();
        Ok(())
    }

    /// `EdgeON(Eon, Par, I)` (`cxx:146-159`): no validity check on `I`.
    pub fn set_edge_on(&mut self, eon: &TopoShape, par: f64, i: i32) {
        if i == 1 {
            self.edge_on1 = eon.clone();
            self.edge_on_par1 = par;
        } else if i == 2 {
            self.edge_on2 = eon.clone();
            self.edge_on_par2 = par;
        }
    }

    /// `EdgeON(I)` (`cxx:161-177`).
    pub fn edge_on(&self, i: i32) -> Result<&TopoShape, String> {
        match i {
            1 => Ok(&self.edge_on1),
            2 => Ok(&self.edge_on2),
            _ => Err("TopOpeBRep_VPointInter::EdgeON".to_string()),
        }
    }

    /// `EdgeONParameter(I)` (`cxx:179-193`).
    pub fn edge_on_parameter(&self, i: i32) -> Result<f64, String> {
        match i {
            1 => Ok(self.edge_on_par1),
            2 => Ok(self.edge_on_par2),
            _ => Err("TopOpeBRep_VPointInter::EdgeONParameter".to_string()),
        }
    }

    /// `ShapeIndex()` (`lxx:183-186`).
    pub fn shape_index(&self) -> i32 {
        self.shape_index
    }

    /// `ShapeIndex(I)` (`lxx:188-191`).
    pub fn set_shape_index(&mut self, i: i32) {
        self.shape_index = i;
    }

    /// `Edge(I)` (`cxx:198-210`): the arc's edge when the point is on that
    /// domain, otherwise the null shape (here `None`).
    pub fn edge(&self, i: i32) -> Option<&TopoShape> {
        let on_dom = match i {
            1 => self.is_on_dom_s1(),
            2 => self.is_on_dom_s2(),
            _ => Ok(false),
        };
        if !on_dom.unwrap_or(false) {
            return None;
        }
        let arc = if i == 1 {
            self.arc_on_s1()
        } else {
            self.arc_on_s2()
        };
        arc.ok()
    }

    /// `EdgeParameter(I)` (`cxx:212-225`).
    pub fn edge_parameter(&self, i: i32) -> f64 {
        let on_dom = match i {
            1 => self.is_on_dom_s1(),
            2 => self.is_on_dom_s2(),
            _ => Ok(false),
        };
        if on_dom.unwrap_or(false) {
            let p = if i == 1 {
                self.parameter_on_arc1()
            } else {
                self.parameter_on_arc2()
            };
            if let Ok(v) = p {
                return v;
            }
        }
        0.0
    }

    /// `SurfaceParameters(I)` (`cxx:227-243`).
    pub fn surface_parameters(&self, i: i32) -> Result<GpPnt2d, String> {
        let (u, v) = if i == 1 {
            self.parameters_on_s1()?
        } else if i == 2 {
            self.parameters_on_s2()?
        } else {
            (0.0, 0.0)
        };
        Ok(GpPnt2d::new(u, v))
    }

    /// `IsVertex(I)` (`cxx:245-262`).
    pub fn is_vertex(&self, i: i32) -> Result<bool, String> {
        if i == 0 {
            return Ok(false);
        }
        if i == 1 && self.is_vertex_on_s1()? {
            return Ok(true);
        } else if i == 2 && self.is_vertex_on_s2()? {
            return Ok(true);
        }
        Ok(false)
    }

    /// `Vertex(I)` (`cxx:264-278`).
    pub fn vertex(&self, i: i32) -> Option<&TopoShape> {
        if i == 1 && self.is_vertex_on_s1().unwrap_or(false) {
            return self.vertex_on_s1().ok();
        } else if i == 2 && self.is_vertex_on_s2().unwrap_or(false) {
            return self.vertex_on_s2().ok();
        }
        None
    }

    /// `UpdateKeep()` (`cxx:280-313`).
    ///
    /// `M_SINON(s)` is `s == IN || s == ON`; the `SI == 3` arm repeats the
    /// `SI == 0` condition, as OCCT does (the `NYI` comment at `cxx:309`).
    pub fn update_keep(&mut self) {
        let m_sinon = |s: State| s == State::In || s == State::On;
        let pos1 = self.state1;
        let pos2 = self.state2;
        let si = self.shape_index;
        let condition = if si == 1 {
            m_sinon(pos2)
        } else if si == 2 {
            m_sinon(pos1)
        } else if si == 0 {
            m_sinon(pos1) && m_sinon(pos2)
        } else if si == 3 {
            m_sinon(pos1) && m_sinon(pos2)
        } else {
            false
        };
        self.keep = condition;
    }

    /// `Keep()` (`lxx:194-197`).
    pub fn keep(&self) -> bool {
        self.keep
    }

    /// `ChangeKeep(keep)` (`lxx:199-202`).
    pub fn change_keep(&mut self, keep: bool) {
        self.keep = keep;
    }

    /// `EqualpP(VP)` (`cxx:316-331`): equal line parameter
    /// (`Precision::PConfusion`) and equal 3D point at the larger of the two
    /// tolerances.
    pub fn equal_pp(&self, vp: &VPointInter) -> Result<bool, String> {
        let p1 = self.parameter_on_line()?;
        let p2 = vp.parameter_on_line()?;
        let pequal = (p1 - p2).abs() < Precision::PCONFUSION;
        let q1 = self.value()?;
        let q2 = vp.value()?;
        let ptol1 = self.tolerance()?;
        let ptol2 = vp.tolerance()?;
        let ptol = if ptol1 > ptol2 { ptol1 } else { ptol2 };
        let pequal3d = q1.distance(q2) <= ptol;
        Ok(pequal && pequal3d)
    }

    /// `ParonE(E, par)` (`cxx:334-377`): the parameter of this vertex on the
    /// edge `E` when `E` is one of its arcs (checked on either domain), else on
    /// one of its `ON`-state edges.
    pub fn paron_e(&self, e: &Edge, par: &mut f64) -> bool {
        if self.is_on_dom_s1().unwrap_or(false) {
            let mut found = false;
            if let Ok(a) = self.arc_on_s1() {
                if is_same(a, &e.0) {
                    found = true;
                }
            }
            if found {
                if let Ok(v) = self.parameter_on_arc1() {
                    *par = v;
                    return true;
                }
            }
        }
        if self.is_on_dom_s2().unwrap_or(false) {
            let mut found = false;
            if let Ok(a) = self.arc_on_s2() {
                if is_same(a, &e.0) {
                    found = true;
                }
            }
            if found {
                if let Ok(v) = self.parameter_on_arc2() {
                    *par = v;
                    return true;
                }
            }
        }

        for i in 1..=2 {
            if self.state(i) != Ok(State::On) {
                continue;
            }
            let same = self
                .edge_on(i)
                .map(|s| is_same(s, &e.0))
                .unwrap_or(false);
            if same {
                if let Ok(v) = self.edge_on_parameter(i) {
                    *par = v;
                    return true;
                }
            }
        }
        false
    }

    /// `Index(I)` (`lxx:298-301`).
    pub fn set_index(&mut self, i: i32) {
        self.index = i;
    }

    /// `Index()` (`lxx:303-306`).
    pub fn index(&self) -> i32 {
        self.index
    }

    /// `Dump(I, F, OS)` (`cxx:380-417`): empty outside `OCCT_DEBUG`.
    pub fn dump_on_edge(&self, _i: i32, _f: &Face) {}

    /// `Dump(F1, F2, OS)` (`cxx:420-478`): empty outside `OCCT_DEBUG`.
    pub fn dump(&self, _f1: &Face, _f2: &Face) {}

    /// `PThePointOfIntersectionDummy()` (`cxx:480-483`).
    pub fn p_the_point_of_intersection_dummy(&self) -> Option<&Arc<IntPatchPoint>> {
        self.ppoi.as_ref()
    }
}
