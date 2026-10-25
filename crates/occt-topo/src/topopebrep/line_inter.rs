//! `TopOpeBRep_LineInter` — one intersection line of two faces: its type, its
//! `VPoint`s (vertices) and its `WPoint`s (walking-line points), plus the
//! face/face transitions of the line.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRep/TopOpeBRep_LineInter.hxx:42-204`,
//! `.lxx:19-303`, `.cxx:19-700`.
//!
//! Translation notes:
//! - `occ::handle<IntPatch_Line> myIL` plus its four down-cast aliases
//!   (`myILA` / `myILR` / `myILW` / `myILG`, `hxx:183-187`) collapse onto the
//!   port's `GeomIntLine` sum type, which already is the `IntPatch_Line`
//!   family; the aliases become `as_aline()` / `as_rline()` / `as_wline()` /
//!   `as_gline()`.
//! - `myHAVP` is an `NCollection_HArray1<TopOpeBRep_VPointInter>` with the OCCT
//!   bounds `0..n` (`cxx:179`), of which only `1..=n` are ever filled or read
//!   (`cxx:181-199`); it becomes a `Vec` indexed 1-based through
//!   [`LineInter::v_point`].
//! - `myNullShape` (`hxx:196`) is a null `TopoDS_Shape`; the port answers
//!   `Option<&TopoShape>` from [`LineInter::arc`] instead (the port has no null
//!   `TopoShape`, see `TopoShape::is_null`).
//! - `TransitionOnS1/2` and `SituationS1/2` (`lxx:83-101`) read
//!   `IntPatch_Line::myTransitionOnS1` / `mySituationS1`, which the port's
//!   `GeomIntLine` carries through [`crate::intpatch_line::IntPatchLine`].
//! - `ComputeFaceFaceTransition` (`cxx:592-598`) calls
//!   `TopOpeBRep_FFTransitionTool::ProcessFaceTransition`
//!   (`topopebrep/ff_transition_tool.rs`).
//! - `Bounds`' `FirstPoint()` / `LastPoint()` (`cxx:420-428`) read the vertex
//!   the `IntPatch_GLine` stores by index (`myIndFirst` / `myIndLast`), which
//!   the port's `GLine` keeps as `first_index` / `last_index`.
//! - `Curve()` (`cxx:480-508`) casts away const to call `SetOK(false)` in its
//!   default arm, so it takes `&mut self` here.
//! - `DumpType` / `DumpVPoint` / `DumpBipoint` / `DumpLineTransitions`
//!   (`cxx:617-700`) write only under `OCCT_DEBUG`; the port has no-op
//!   counterparts.

use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_geom::{Curve, GeomCircle, GeomEllipse, GeomHyperbola, GeomLine, GeomParabola, GeomTrimmedCurve, Surface};

use crate::geom_int::{GLine, GLineKind, GeomIntLine, RLine};
use crate::int_tools_wline::{PntOn2S, WLine};
use crate::intpatch::aline_to_wline::ALineToWLine;
use crate::intpatch_point::IntPatchPoint;
use crate::intsurf_transition::{IntSurfSituation, IntSurfTypeTrans};
use crate::shape::{Face, TopoShape};
use crate::topopebrep::ds::Transition;
use crate::topopebrep::ff_transition_tool::FFTransitionTool;
use crate::topopebrep::interf::{VPointInter, WPointInter};
use crate::topopebrep::type_line_curve::TypeLineCurve;

/// The port's `VPointInter` accessors are fallible while OCCT dereferences
/// `myPPOI`; a `VPoint` built by [`LineInter::set_line`] always has one.
fn param_on_line(vp: &VPointInter) -> f64 {
    vp.parameter_on_line().unwrap_or(0.0)
}

/// `class TopOpeBRep_LineInter` (`hxx:42-204`).
#[derive(Clone)]
pub struct LineInter {
    /// `myOK` (`hxx:165`).
    ok: bool,
    /// `myIndex` (`hxx:166`).
    index: i32,
    /// `myNbVPoint` (`hxx:167`).
    nb_v_point: i32,
    /// `myIsVClosed` (`hxx:168`).
    is_v_closed: bool,
    /// `myHasVPonR` (`hxx:169`).
    has_v_pon_r: bool,
    /// `myINL` (`hxx:170`): all `VPoint`s share one parameter on the line.
    inl: bool,
    /// `myVPBDefined` (`hxx:171`).
    vpb_defined: bool,
    /// `myVPF` (`hxx:172`).
    vpf: i32,
    /// `myVPL` (`hxx:173`).
    vpl: i32,
    /// `myVPN` (`hxx:174`).
    vpn: i32,
    /// `myTypeLineCurve` (`hxx:175`).
    type_line_curve: TypeLineCurve,
    /// `myIL` and its four down-casts (`hxx:176-180`).
    line: Option<GeomIntLine>,
    /// `myCurrentWP` (`hxx:181`).
    current_wp: WPointInter,
    /// `myHAVP` (`hxx:182`), 0-based storage for the OCCT 1-based array.
    v_points: Vec<VPointInter>,
    /// `myF1` (`hxx:188`).
    f1: Option<Face>,
    /// `myF2` (`hxx:189`).
    f2: Option<Face>,
    /// `myLineTonF1` (`hxx:190`).
    line_t_on_f1: Transition,
    /// `myLineTonF2` (`hxx:191`).
    line_t_on_f2: Transition,
    /// `myexF1` (`hxx:197`).
    ex_f1: i32,
    /// `myexF2` (`hxx:198`).
    ex_f2: i32,
}

impl Default for LineInter {
    /// `TopOpeBRep_LineInter()` (`lxx:21-39`).
    fn default() -> Self {
        Self {
            ok: true,
            index: 0,
            nb_v_point: 0,
            is_v_closed: false,
            has_v_pon_r: false,
            inl: false,
            vpb_defined: false,
            vpf: 0,
            vpl: 0,
            vpn: 0,
            type_line_curve: TypeLineCurve::OtherType,
            line: None,
            current_wp: WPointInter::new(),
            v_points: Vec::new(),
            f1: None,
            f2: None,
            line_t_on_f1: Transition::default(),
            line_t_on_f2: Transition::default(),
            ex_f1: 0,
            ex_f2: 0,
        }
    }
}

impl LineInter {
    /// `TopOpeBRep_LineInter()` (`lxx:21-39`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `SetLine(L, S1, S2)` (`cxx:75-200`).
    pub fn set_line(&mut self, l: &GeomIntLine, s1: &dyn Surface, s2: &dyn Surface) {
        // Type from `L->ArcType()` (`cxx:81-113`).
        self.type_line_curve = match l.arc_type() {
            crate::geom_int::IntPatchIType::Analytic => TypeLineCurve::Analytic,
            crate::geom_int::IntPatchIType::Restriction => TypeLineCurve::Restriction,
            crate::geom_int::IntPatchIType::Walking => TypeLineCurve::Walking,
            crate::geom_int::IntPatchIType::Lin => TypeLineCurve::Line,
            crate::geom_int::IntPatchIType::Circle => TypeLineCurve::Circle,
            crate::geom_int::IntPatchIType::Ellipse => TypeLineCurve::Ellipse,
            crate::geom_int::IntPatchIType::Parabola => TypeLineCurve::Parabola,
            crate::geom_int::IntPatchIType::Hyperbola => TypeLineCurve::Hyperbola,
        };
        self.line = Some(l.clone());

        // `cxx:130-142`: an analytic line becomes a walking line.
        if self.type_line_curve == TypeLineCurve::Analytic {
            if let Some(al) = l.as_aline() {
                let wlines = ALineToWLine::new(s1, s2, 200).make_wline(al);
                if let Some(w) = wlines.into_iter().next() {
                    self.line = Some(GeomIntLine::Walking(w));
                }
            }
            self.type_line_curve = TypeLineCurve::Walking;
        }

        // `cxx:145-176`: number of points found on the restriction(s).
        let line = match &self.line {
            Some(l) => l,
            None => {
                self.set_ok(false);
                return;
            }
        };
        self.nb_v_point = line.nb_vertex();

        // `cxx:179-199`: build the VPoints.
        let arc1 = line.as_rline().and_then(|r| r.arc_adaptor());
        let arc2 = line.as_rline().and_then(|r| r.arc_adaptor());
        self.v_points.clear();
        for i in 1..=self.nb_v_point {
            let mut vp = VPointInter::new();
            if let Some(pp) = line.vertex(i) {
                let pt = IntPatchPoint::from_patch_point(pp, arc1.clone(), arc2.clone());
                vp.set_point(Arc::new(pt));
            }
            vp.set_index(i);
            self.v_points.push(vp);
        }
    }

    /// `SetFaces(F1, F2)` (`lxx:43-47`).
    pub fn set_faces(&mut self, f1: &Face, f2: &Face) {
        self.f1 = Some(f1.clone());
        self.f2 = Some(f2.clone());
    }

    /// `TypeLineCurve()` (`lxx:51-54`).
    pub fn type_line_curve(&self) -> TypeLineCurve {
        self.type_line_curve
    }

    /// `NbVPoint()` (`lxx:58-61`).
    pub fn nb_v_point(&self) -> i32 {
        self.nb_v_point
    }

    /// `VPoint(I)` (`cxx:204-207`), 1-based.
    pub fn v_point(&self, i: i32) -> Option<&VPointInter> {
        if i < 1 {
            return None;
        }
        self.v_points.get((i - 1) as usize)
    }

    /// `ChangeVPoint(I)` (`cxx:212-215`), 1-based.
    pub fn change_v_point(&mut self, i: i32) -> Option<&mut VPointInter> {
        if i < 1 {
            return None;
        }
        self.v_points.get_mut((i - 1) as usize)
    }

    /// `SetINL()` (`cxx:219-243`): true when every `VPoint` has the same
    /// parameter on the line.
    pub fn set_inl(&mut self) {
        let mut it = self.v_points.iter();
        let p0 = match it.next() {
            Some(vp) => param_on_line(vp),
            None => {
                self.inl = false;
                return;
            }
        };
        match it.next() {
            None => {
                self.inl = true;
                return;
            }
            Some(_) => {}
        }
        for vp in self.v_points.iter() {
            if param_on_line(vp) != p0 {
                self.inl = false;
                return;
            }
        }
        self.inl = true;
    }

    /// `INL()` (`lxx:63-66`).
    pub fn inl(&self) -> bool {
        self.inl
    }

    /// `SetIsVClosed()` (`cxx:247-311`).
    pub fn set_is_v_closed(&mut self) {
        if self.inl {
            self.is_v_closed = false;
            return;
        }

        let n_v = self.nb_v_point;
        let mut pmin = f64::MAX;
        let mut pmax = f64::MIN;
        let mut imin = 0;
        let mut imax = 0;
        if n_v >= 2 {
            for (k, vp) in self.v_points.iter().enumerate() {
                if !vp.is_internal() {
                    let i = (k + 1) as i32;
                    let p = param_on_line(vp);
                    if p < pmin {
                        imin = i;
                        pmin = p;
                    }
                    if p > pmax {
                        imax = i;
                        pmax = p;
                    }
                }
            }
            if imax == 0 {
                // No VPoint on a restriction: the line is closed.
                self.is_v_closed = true;
                return;
            }

            let (pp1, tol1) = match self.v_point(imin) {
                Some(p) => (p.value().copied().unwrap_or(GpPnt::zero()), p.tolerance().unwrap_or(0.0)),
                None => (GpPnt::zero(), 0.0),
            };
            let (pp2, tol2) = match self.v_point(imax) {
                Some(p) => (p.value().copied().unwrap_or(GpPnt::zero()), p.tolerance().unwrap_or(0.0)),
                None => (GpPnt::zero(), 0.0),
            };
            let tol = if tol1 > tol2 { tol1 } else { tol2 };
            self.is_v_closed = pp1.distance(&pp2) <= tol;
        } else {
            self.set_ok(false);
            self.is_v_closed = false;
        }
    }

    /// `IsVClosed()` (`lxx:68-71`).
    pub fn is_v_closed(&self) -> bool {
        self.is_v_closed
    }

    /// `SetOK(B)` (`cxx:664-667`).
    pub fn set_ok(&mut self, b: bool) {
        self.ok = b;
    }

    /// `OK()` (`lxx:107-110`).
    pub fn ok(&self) -> bool {
        self.ok
    }

    /// `SetHasVPonR()` (`cxx:315-327`).
    pub fn set_has_v_pon_r(&mut self) {
        self.has_v_pon_r = false;
        for vp in self.v_points.iter() {
            if vp.is_on_dom_s1().unwrap_or(false) || vp.is_on_dom_s2().unwrap_or(false) {
                self.has_v_pon_r = true;
                break;
            }
        }
    }

    /// `HasVPonR()` (`lxx:73-76`).
    pub fn has_v_pon_r(&self) -> bool {
        self.has_v_pon_r
    }

    /// `SetVPBounds()` (`cxx:333-360`).
    pub fn set_vp_bounds(&mut self) {
        self.vpf = 0;
        self.vpl = 0;
        self.vpn = 0;
        self.vpb_defined = true;

        let mut f = self.nb_v_point + 1;
        let mut l = 0;
        let mut n = 0;
        for (k, vp) in self.v_points.iter().enumerate() {
            if vp.keep() {
                n += 1;
                let i = (k + 1) as i32;
                if i < f {
                    f = i;
                }
                if i > l {
                    l = i;
                }
            }
        }
        self.vpf = f;
        self.vpl = l;
        self.vpn = n;
    }

    /// `VPBounds(f, l, n)` (`cxx:364-377`): `None` when the bounds were never
    /// computed, which OCCT reports by clearing `myOK` and answering zeros.
    pub fn vp_bounds(&mut self) -> Option<(i32, i32, i32)> {
        if !self.vpb_defined {
            self.set_ok(false);
            return None;
        }
        Some((self.vpf, self.vpl, self.vpn))
    }

    /// `IsPeriodic()` (`cxx:380-392`).
    pub fn is_periodic(&self) -> bool {
        matches!(
            self.type_line_curve,
            TypeLineCurve::Circle | TypeLineCurve::Ellipse
        )
    }

    /// `Period()` (`cxx:396-401`).
    pub fn period(&mut self) -> f64 {
        let (first, last) = self.bounds();
        last - first
    }

    /// `Bounds(theFirst, theLast)` (`cxx:404-429`).
    pub fn bounds(&mut self) -> (f64, f64) {
        let mut first = 0.0;
        let mut last = 0.0;
        let line = match self.line.as_ref().and_then(|l| l.as_gline()) {
            Some(g) => g.clone(),
            None => {
                self.set_ok(false);
                return (0.0, 0.0);
            }
        };

        if self.is_periodic() {
            if let Some(c) = self.curve() {
                last = c.period();
            }
        }

        if line.has_first_point {
            first = line.first_point().param_on_line;
        }
        if line.has_last_point {
            last = line.last_point().param_on_line;
        }
        (first, last)
    }

    /// `HasVInternal()` (`cxx:434-445`).
    pub fn has_v_internal(&self) -> bool {
        self.v_points.iter().any(|vp| vp.is_internal())
    }

    /// `NbWPoint()` (`cxx:450-458`).
    pub fn nb_w_point(&self) -> i32 {
        match self.line.as_ref() {
            Some(GeomIntLine::Walking(w)) => w.nb_pnts(),
            _ => 0,
        }
    }

    /// `WPoint(IW)` (`cxx:462-476`).
    ///
    /// The `RESTRICTION` arm (`cxx:468-469`) reads
    /// `IntPatch_RLine::Point(IW)`, i.e. the `IntSurf_LineOn2S` an
    /// `IntPatch_PointLine` stores; the port's `RLine` keeps the arc's p-curve
    /// and the vertices but not that point list: UNPORTED, `myCurrentWP` is
    /// left as it was.
    pub fn w_point(&mut self, iw: i32) -> &WPointInter {
        match self.line.as_ref() {
            Some(GeomIntLine::Walking(w)) => {
                if let Some(p) = point_on_wline(w, iw) {
                    self.current_wp.set(&p);
                }
            }
            _ => {}
        }
        &self.current_wp
    }

    /// `TransitionOnS1()` (`lxx:83-86`): `myIL->TransitionOnS1()`.
    pub fn transition_on_s1(&self) -> IntSurfTypeTrans {
        match self.line.as_ref() {
            Some(l) => l.transition_on_s1(),
            None => IntSurfTypeTrans::Undecided,
        }
    }

    /// `TransitionOnS2()` (`lxx:88-91`): `myIL->TransitionOnS2()`.
    pub fn transition_on_s2(&self) -> IntSurfTypeTrans {
        match self.line.as_ref() {
            Some(l) => l.transition_on_s2(),
            None => IntSurfTypeTrans::Undecided,
        }
    }

    /// `SituationS1()` (`lxx:93-96`): `myIL->SituationS1()`, which raises
    /// `Standard_DomainError` unless the transition is `Touch`.
    pub fn situation_s1(&self) -> IntSurfSituation {
        match self.line.as_ref() {
            Some(l) => l.situation_s1().unwrap_or(IntSurfSituation::Unknown),
            None => IntSurfSituation::Unknown,
        }
    }

    /// `SituationS2()` (`lxx:98-101`): `myIL->SituationS2()`.
    pub fn situation_s2(&self) -> IntSurfSituation {
        match self.line.as_ref() {
            Some(l) => l.situation_s2().unwrap_or(IntSurfSituation::Unknown),
            None => IntSurfSituation::Unknown,
        }
    }

    /// `Curve()` (`cxx:480-508`): the 3D curve of a conic line.
    pub fn curve(&mut self) -> Option<Arc<dyn Curve>> {
        let kind = match self.line.as_ref().and_then(|l| l.as_gline()) {
            Some(g) => g.kind.clone(),
            None => {
                self.set_ok(false);
                return None;
            }
        };
        match kind {
            GLineKind::Lin(l) => Some(Arc::new(GeomLine::new(l))),
            GLineKind::Circ(c) => Some(Arc::new(GeomCircle::new(c))),
            GLineKind::Elips(e) => Some(Arc::new(GeomEllipse::new(e))),
            GLineKind::Parab(p) => Some(Arc::new(GeomParabola::new(p))),
            GLineKind::Hypr(h) => Some(Arc::new(GeomHyperbola::new(h))),
        }
    }

    /// `Curve(parmin, parmax)` (`cxx:512-527`): `Curve()` trimmed to the range.
    pub fn curve_trimmed(&mut self, parmin: f64, parmax: f64) -> Option<Arc<dyn Curve>> {
        let c3d = self.curve()?;
        Some(Arc::new(GeomTrimmedCurve::new(c3d, parmin, parmax)))
    }

    /// `Arc()` (`cxx:533-554`): the edge of a `RESTRICTION` line's arc, `None`
    /// otherwise.
    pub fn arc(&self) -> Option<&TopoShape> {
        if self.type_line_curve != TypeLineCurve::Restriction {
            return None;
        }
        let r = self.line.as_ref()?.as_rline()?;
        r.arc_edge.as_ref().map(|e| &e.0)
    }

    /// `ArcIsEdge(I)` (`cxx:558-565`).
    pub fn arc_is_edge(&self, index: i32) -> bool {
        if self.type_line_curve != TypeLineCurve::Restriction {
            return false;
        }
        let r = match self.line.as_ref().and_then(|l| l.as_rline()) {
            Some(r) => r,
            None => return false,
        };
        if index == 2 {
            !r.arc_on_s1
        } else {
            r.arc_on_s1
        }
    }

    /// `HasFirstPoint()` (`cxx:569-575`).
    pub fn has_first_point(&self) -> Result<bool, String> {
        match self.line.as_ref().and_then(|l| l.as_gline()) {
            Some(g) => Ok(g.has_first_point),
            None => Err(
                "TopOpeBRep_LineInter::HasFirstPoint sur line != GLine".to_string(),
            ),
        }
    }

    /// `HasLastPoint()` (`cxx:580-586`).
    pub fn has_last_point(&self) -> Result<bool, String> {
        match self.line.as_ref().and_then(|l| l.as_gline()) {
            Some(g) => Ok(g.has_last_point),
            None => Err("TopOpeBRep_LineInter::HasLastPoint sur line != GLine".to_string()),
        }
    }

    /// `ComputeFaceFaceTransition()` (`cxx:592-598`): the transition of the
    /// line on `F1` is computed against `F2`'s orientation and vice versa.
    pub fn compute_face_face_transition(&mut self) {
        let (Some(f1), Some(f2)) = (self.f1.as_ref(), self.f2.as_ref()) else {
            return;
        };
        let f1_ori = f1.orientation();
        let f2_ori = f2.orientation();
        let t1 = FFTransitionTool::process_face_transition(self, 1, f2_ori);
        let t2 = FFTransitionTool::process_face_transition(self, 2, f1_ori);
        self.line_t_on_f1 = t1;
        self.line_t_on_f2 = t2;
    }

    /// `FaceFaceTransition(I)` (`cxx:602-612`).
    pub fn face_face_transition(&self, i: i32) -> Result<&Transition, String> {
        match i {
            1 => Ok(&self.line_t_on_f1),
            2 => Ok(&self.line_t_on_f2),
            _ => Err("TopOpeBRep_LineInter::FaceFaceTransition".to_string()),
        }
    }

    /// `LineW()` (`lxx:113-116`).
    pub fn line_w(&self) -> Option<&WLine> {
        self.line.as_ref().and_then(|l| l.as_wline())
    }

    /// `LineG()` (`lxx:118-121`).
    pub fn line_g(&self) -> Option<&GLine> {
        self.line.as_ref().and_then(|l| l.as_gline())
    }

    /// `LineR()` (`lxx:123-126`).
    pub fn line_r(&self) -> Option<&RLine> {
        self.line.as_ref().and_then(|l| l.as_rline())
    }

    /// `Index(I)` (`lxx:128-131`).
    pub fn set_index(&mut self, i: i32) {
        self.index = i;
    }

    /// `Index()` (`lxx:133-136`).
    pub fn index(&self) -> i32 {
        self.index
    }

    /// `DumpType()` (`cxx:617-624`): empty outside `OCCT_DEBUG`.
    pub fn dump_type(&self) {}

    /// `DumpVPoint(I, s1, s2)` (`cxx:628-642`): empty outside `OCCT_DEBUG`.
    pub fn dump_v_point(&self, _i: i32, _s1: &str, _s2: &str) {}

    /// `DumpBipoint(B, s1, s2)` (`cxx:647-661`): empty outside `OCCT_DEBUG`.
    pub fn dump_bipoint(&self, _b: &crate::topopebrep::interf::Bipoint, _s1: &str, _s2: &str) {}

    /// `SetTraceIndex(exF1, exF2)` (`cxx:672-676`).
    pub fn set_trace_index(&mut self, ex_f1: i32, ex_f2: i32) {
        self.ex_f1 = ex_f1;
        self.ex_f2 = ex_f2;
    }

    /// `GetTraceIndex(exF1, exF2)` (`cxx:680-684`).
    pub fn get_trace_index(&self) -> (i32, i32) {
        (self.ex_f1, self.ex_f2)
    }

    /// `DumpLineTransitions(OS)` (`cxx:687-700`): empty outside `OCCT_DEBUG`.
    pub fn dump_line_transitions(&self) {}
}

/// `myILW->Point(IW)` (`cxx:470`).
fn point_on_wline(w: &WLine, iw: i32) -> Option<PntOn2S> {
    if iw < 1 || iw > w.nb_pnts() {
        return None;
    }
    Some(*w.point(iw))
}
