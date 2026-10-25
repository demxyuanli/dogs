//! `TopOpeBRep_Point2d` — a 2D (edge/edge or face/face) intersection point with
//! its parameter, transition and ancestor data on both shapes (1) and (2).
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRep/TopOpeBRep_Point2d.hxx:40-144`,
//! `.lxx:20-280`, `.cxx:26-176`.
//!
//! Translation notes:
//! - `IntRes2d_IntersectionPoint mypint` + `bool myhaspint` become
//!   `Option<IntRes2dIntersectionPoint>`; `Pint()` raises when unset, as
//!   `lxx:38-40` does.
//! - `TopoDS_Vertex myvertex1/2` are default-constructed (null) in OCCT; here
//!   they are `Option<Vertex>` and `Vertex(1)` raises while unset, which is the
//!   OCCT `IsVertex(1) == false` case.
//! - `occ::handle<TopOpeBRep_Hctxff2d> / <TopOpeBRep_Hctxee2d>` become
//!   `Option<Hctxff2d> / Option<Hctxee2d>` (a null handle is `None`).
//! - The index-validated accessors (`IsVertex`, `SetVertex`, `Parameter`,
//!   `Vertex`, `Transition`, `ChangeTransition`) throw `Standard_Failure` for
//!   an index other than 1 or 2 (`lxx:52-55`, `:66-68`, `cxx:44-64`, `:66-84`);
//!   the port returns `Err`.
//! - `Dump` (`cxx:104-176`) exists only under `OCCT_DEBUG` (the non-debug build
//!   is an empty body); it is not translated.

use occt_core::gp::{GpPnt, GpPnt2d};
use occt_core::intres2d::IntRes2dIntersectionPoint;

use crate::shape::Vertex;
use crate::topopebrep::ds::{Config, Transition};

use super::hctxee2d::Hctxee2d;
use super::hctxff2d::Hctxff2d;
use super::p2d_status::P2DStatus;

/// `class TopOpeBRep_Point2d` (`hxx:40-144`).
#[derive(Clone)]
pub struct Point2d {
    /// `mypint` (`hxx:120`)
    pint: Option<IntRes2dIntersectionPoint>,
    /// `myisvertex1` (`hxx:123`)
    is_vertex1: bool,
    /// `myvertex1` (`hxx:124`)
    vertex1: Option<Vertex>,
    /// `mytransition1` (`hxx:125`)
    transition1: Transition,
    /// `myparameter1` (`hxx:126`)
    parameter1: f64,
    /// `myisvertex2` (`hxx:127`)
    is_vertex2: bool,
    /// `myvertex2` (`hxx:128`)
    vertex2: Option<Vertex>,
    /// `mytransition2` (`hxx:129`)
    transition2: Transition,
    /// `myparameter2` (`hxx:130`)
    parameter2: f64,
    /// `myispointofsegment` (`hxx:131`)
    is_point_of_segment: bool,
    /// `myips1` (`hxx:132`)
    ips1: i32,
    /// `myips2` (`hxx:133`)
    ips2: i32,
    /// `myhasancestors` (`hxx:134`)
    has_ancestors: bool,
    /// `mystatus` (`hxx:135`)
    status: P2DStatus,
    /// `myindex` (`hxx:136`)
    index: i32,
    /// `mypnt` (`hxx:137`)
    pnt: GpPnt,
    /// `mypnt2d` (`hxx:138`)
    pnt2d: GpPnt2d,
    /// `mykeep` (`hxx:139`)
    keep: bool,
    /// `myedgesconfig` (`hxx:140`)
    edges_config: Config,
    /// `mytolerance` (`hxx:141`)
    tolerance: f64,
    /// `myctxff2d` (`hxx:142`)
    ctxff2d: Option<Hctxff2d>,
    /// `myctxee2d` (`hxx:143`)
    ctxee2d: Option<Hctxee2d>,
}

impl Default for Point2d {
    /// `TopOpeBRep_Point2d::TopOpeBRep_Point2d()` (`cxx:26-41`): exactly the
    /// members of the OCCT init list; `mypint`, `myvertex1/2`, the two
    /// transitions, `mypnt`, `mypnt2d` and the two context handles keep their
    /// class defaults.
    fn default() -> Self {
        Self {
            pint: None,
            is_vertex1: false,
            vertex1: None,
            transition1: Transition::default(),
            parameter1: 0.0,
            is_vertex2: false,
            vertex2: None,
            transition2: Transition::default(),
            parameter2: 0.0,
            is_point_of_segment: false,
            ips1: 0,
            ips2: 0,
            has_ancestors: false,
            status: P2DStatus::Unknown,
            index: 0,
            pnt: GpPnt::new(0.0, 0.0, 0.0),
            pnt2d: GpPnt2d::new(0.0, 0.0),
            keep: true,
            edges_config: Config::UnshGeometry,
            tolerance: 0.0,
            ctxff2d: None,
            ctxee2d: None,
        }
    }
}

impl Point2d {
    /// `TopOpeBRep_Point2d()` (`cxx:26-41`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `SetPint(P)` (`lxx:22-26`).
    pub fn set_pint(&mut self, p: &IntRes2dIntersectionPoint) {
        self.pint = Some(*p);
    }

    /// `HasPint()` (`lxx:30-33`).
    pub fn has_pint(&self) -> bool {
        self.pint.is_some()
    }

    /// `Pint()` (`lxx:36-41`).
    pub fn pint(&self) -> Result<&IntRes2dIntersectionPoint, String> {
        self.pint
            .as_ref()
            .ok_or_else(|| "TopOpeBRep_Point2d::Pint".to_string())
    }

    /// `SetIsVertex(Index, B)` (`lxx:45-55`).
    pub fn set_is_vertex(&mut self, index: i32, b: bool) -> Result<(), String> {
        match index {
            1 => self.is_vertex1 = b,
            2 => self.is_vertex2 = b,
            _ => return Err("TopOpeBRep_Point2d::SetIsVertex".to_string()),
        }
        Ok(())
    }

    /// `IsVertex(Index)` (`lxx:59-68`).
    pub fn is_vertex(&self, index: i32) -> Result<bool, String> {
        match index {
            1 => Ok(self.is_vertex1),
            2 => Ok(self.is_vertex2),
            _ => Err("TopOpeBRep_Point2d::IsVertex".to_string()),
        }
    }

    /// `SetVertex(Index, V)` (`lxx:72-84`): also latches `myisvertexN`.
    pub fn set_vertex(&mut self, index: i32, v: &Vertex) -> Result<(), String> {
        match index {
            1 => {
                self.vertex1 = Some(v.clone());
                self.is_vertex1 = true;
            }
            2 => {
                self.vertex2 = Some(v.clone());
                self.is_vertex2 = true;
            }
            _ => return Err("TopOpeBRep_Point2d::SetVertex".to_string()),
        }
        Ok(())
    }

    /// `Vertex(Index)` (`cxx:44-64`): raises when the point is not a vertex of
    /// that shape, then for an index other than 1 or 2.
    pub fn vertex(&self, index: i32) -> Result<&Vertex, String> {
        if !self.is_vertex(index)? {
            return Err("TopOpeBRep_Point2d::Vertex".to_string());
        }
        match index {
            1 => self
                .vertex1
                .as_ref()
                .ok_or_else(|| "TopOpeBRep_Point2d::Vertex".to_string()),
            2 => self
                .vertex2
                .as_ref()
                .ok_or_else(|| "TopOpeBRep_Point2d::Vertex".to_string()),
            _ => Err("TopOpeBRep_Point2d::Vertex".to_string()),
        }
    }

    /// `SetTransition(Index, T)` (`lxx:88-97`).
    pub fn set_transition(&mut self, index: i32, t: &Transition) -> Result<(), String> {
        match index {
            1 => self.transition1 = t.clone(),
            2 => self.transition2 = t.clone(),
            _ => return Err("TopOpeBRep_Point2d::SetTransition".to_string()),
        }
        Ok(())
    }

    /// `Transition(Index)` (`cxx:66-84`).
    pub fn transition(&self, index: i32) -> Result<&Transition, String> {
        match index {
            1 => Ok(&self.transition1),
            2 => Ok(&self.transition2),
            _ => Err("TopOpeBRep_Point2d::Transition".to_string()),
        }
    }

    /// `ChangeTransition(Index)` (`cxx:86-100`).
    pub fn change_transition(&mut self, index: i32) -> Result<&mut Transition, String> {
        match index {
            1 => Ok(&mut self.transition1),
            2 => Ok(&mut self.transition2),
            _ => Err("TopOpeBRep_Point2d::ChangeTransition".to_string()),
        }
    }

    /// `SetParameter(Index, P)` (`lxx:99-107`).
    pub fn set_parameter(&mut self, index: i32, p: f64) -> Result<(), String> {
        match index {
            1 => self.parameter1 = p,
            2 => self.parameter2 = p,
            _ => return Err("TopOpeBRep_Point2d::SetParameter".to_string()),
        }
        Ok(())
    }

    /// `Parameter(Index)` (`lxx:111-119`).
    pub fn parameter(&self, index: i32) -> Result<f64, String> {
        match index {
            1 => Ok(self.parameter1),
            2 => Ok(self.parameter2),
            _ => Err("TopOpeBRep_Point2d::Parameter".to_string()),
        }
    }

    /// `SetIsPointOfSegment(B)` (`lxx:123-126`).
    pub fn set_is_point_of_segment(&mut self, b: bool) {
        self.is_point_of_segment = b;
    }

    /// `IsPointOfSegment()` (`lxx:130-133`).
    pub fn is_point_of_segment(&self) -> bool {
        self.is_point_of_segment
    }

    /// `SetSegmentAncestors(IP1, IP2)` (`lxx:137-145`).
    pub fn set_segment_ancestors(&mut self, ip1: i32, ip2: i32) {
        self.ips1 = ip1;
        self.ips2 = ip2;
        self.has_ancestors = self.ips1 != 0 && self.ips2 != 0;
        if self.has_ancestors {
            self.status = P2DStatus::New;
        }
    }

    /// `SegmentAncestors(IP1, IP2)` (`lxx:149-154`).
    pub fn segment_ancestors(&self) -> (i32, i32, bool) {
        (self.ips1, self.ips2, self.has_ancestors)
    }

    /// `SetStatus(S)` (`lxx:157-160`).
    pub fn set_status(&mut self, s: P2DStatus) {
        self.status = s;
    }

    /// `Status()` (`lxx:163-166`).
    pub fn status(&self) -> P2DStatus {
        self.status
    }

    /// `SetIndex(X)` (`lxx:171-174`).
    pub fn set_index(&mut self, x: i32) {
        self.index = x;
    }

    /// `Index()` (`lxx:178-181`).
    pub fn index(&self) -> i32 {
        self.index
    }

    /// `SetValue(P)` (`lxx:185-188`).
    pub fn set_value(&mut self, p: &GpPnt) {
        self.pnt = *p;
    }

    /// `Value()` (`lxx:192-195`).
    pub fn value(&self) -> &GpPnt {
        &self.pnt
    }

    /// `SetValue2d(P)` (`lxx:199-202`).
    pub fn set_value2d(&mut self, p: &GpPnt2d) {
        self.pnt2d = *p;
    }

    /// `Value2d()` (`lxx:206-209`).
    pub fn value2d(&self) -> &GpPnt2d {
        &self.pnt2d
    }

    /// `SetKeep(B)` (`lxx:213-216`).
    pub fn set_keep(&mut self, b: bool) {
        self.keep = b;
    }

    /// `Keep()` (`lxx:220-223`).
    pub fn keep(&self) -> bool {
        self.keep
    }

    /// `SetEdgesConfig(C)` (`lxx:227-230`).
    pub fn set_edges_config(&mut self, c: Config) {
        self.edges_config = c;
    }

    /// `EdgesConfig()` (`lxx:234-237`).
    pub fn edges_config(&self) -> Config {
        self.edges_config
    }

    /// `SetTolerance(T)` (`lxx:241-244`).
    pub fn set_tolerance(&mut self, t: f64) {
        self.tolerance = t;
    }

    /// `Tolerance()` (`lxx:248-251`).
    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }

    /// `SetHctxff2d(h)` (`lxx:255-258`).
    pub fn set_hctxff2d(&mut self, ff2d: Option<Hctxff2d>) {
        self.ctxff2d = ff2d;
    }

    /// `Hctxff2d()` (`lxx:262-265`).
    pub fn hctxff2d(&self) -> Option<&Hctxff2d> {
        self.ctxff2d.as_ref()
    }

    /// `SetHctxee2d(h)` (`lxx:269-272`).
    pub fn set_hctxee2d(&mut self, ee2d: Option<Hctxee2d>) {
        self.ctxee2d = ee2d;
    }

    /// `Hctxee2d()` (`lxx:276-279`).
    pub fn hctxee2d(&self) -> Option<&Hctxee2d> {
        self.ctxee2d.as_ref()
    }
}
