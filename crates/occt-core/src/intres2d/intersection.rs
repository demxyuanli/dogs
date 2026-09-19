//! Port of the `IntRes2d` result classes
//! (`src/ModelingAlgorithms/TKGeomAlgo/IntRes2d/`):
//! `IntRes2d_IntersectionPoint.hxx/.cxx/.lxx`,
//! `IntRes2d_IntersectionSegment.hxx/.cxx/.lxx`,
//! `IntRes2d_Intersection.hxx/.cxx/.lxx`.
//!
//! These are the only result shapes `Geom2dInt_GInter` produces, and the
//! transitions they carry are what `ShapeAnalysis_Wire` reads
//! (`ShapeAnalysis_Wire.cxx:1324-1329`, `:1490-1493`).

use crate::gp::GpPnt2d;

use super::transition::{IntRes2dPosition, IntRes2dTransition, IntRes2dTypeTrans};

/// `PARAMEQUAL(a, b)` (`IntRes2d_Intersection.cxx:26`).
fn paramequal(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-8
}

/// `IntRes2d_IntersectionPoint` (`IntRes2d_IntersectionPoint.hxx:30-86`).
#[derive(Clone, Copy, Debug)]
pub struct IntRes2dIntersectionPoint {
    /// `pt` (`hxx:81`)
    pt: GpPnt2d,
    /// `p1` (`hxx:82`)
    p1: f64,
    /// `p2` (`hxx:83`)
    p2: f64,
    /// `trans1` (`hxx:84`)
    trans1: IntRes2dTransition,
    /// `trans2` (`hxx:85`)
    trans2: IntRes2dTransition,
}

impl Default for IntRes2dIntersectionPoint {
    /// `IntRes2d_IntersectionPoint()` (`IntersectionPoint.cxx:20-25`):
    /// `p1` / `p2` are `RealLast()`.
    fn default() -> Self {
        Self::new()
    }
}

impl IntRes2dIntersectionPoint {
    /// `IntRes2d_IntersectionPoint()` (`IntersectionPoint.cxx:20-25`).
    pub fn new() -> Self {
        Self {
            pt: GpPnt2d::new(0.0, 0.0),
            p1: f64::MAX,
            p2: f64::MAX,
            trans1: IntRes2dTransition::new(),
            trans2: IntRes2dTransition::new(),
        }
    }

    /// `IntRes2d_IntersectionPoint(P, Uc1, Uc2, Trans1, Trans2, ReversedFlag)`
    /// (`IntersectionPoint.lxx:18-36`).
    pub fn with_transitions(
        p: &GpPnt2d,
        uc1: f64,
        uc2: f64,
        trans1: &IntRes2dTransition,
        trans2: &IntRes2dTransition,
        reversed_flag: bool,
    ) -> Self {
        let mut r = Self {
            pt: *p,
            p1: uc1,
            p2: uc2,
            trans1: *trans1,
            trans2: *trans2,
        };
        if reversed_flag {
            r.trans1 = *trans2;
            r.trans2 = *trans1;
            r.p1 = uc2;
            r.p2 = uc1;
        }
        r
    }

    /// `SetValues(P, Uc1, Uc2, Trans1, Trans2, ReversedFlag)`
    /// (`IntersectionPoint.lxx:37-60`).
    pub fn set_values(
        &mut self,
        p: &GpPnt2d,
        uc1: f64,
        uc2: f64,
        trans1: &IntRes2dTransition,
        trans2: &IntRes2dTransition,
        reversed_flag: bool,
    ) {
        self.pt = *p;
        if !reversed_flag {
            self.trans1 = *trans1;
            self.trans2 = *trans2;
            self.p1 = uc1;
            self.p2 = uc2;
        } else {
            self.trans1 = *trans2;
            self.trans2 = *trans1;
            self.p1 = uc2;
            self.p2 = uc1;
        }
    }

    /// `Value()` (`IntersectionPoint.lxx:62-65`).
    pub fn value(&self) -> &GpPnt2d {
        &self.pt
    }

    /// `ParamOnFirst()` (`IntersectionPoint.lxx:67-70`).
    pub fn param_on_first(&self) -> f64 {
        self.p1
    }

    /// `ParamOnSecond()` (`IntersectionPoint.lxx:72-75`).
    pub fn param_on_second(&self) -> f64 {
        self.p2
    }

    /// `TransitionOfFirst()` (`IntersectionPoint.lxx:77-80`).
    pub fn transition_of_first(&self) -> &IntRes2dTransition {
        &self.trans1
    }

    /// `TransitionOfSecond()` (`IntersectionPoint.lxx:82-85`).
    pub fn transition_of_second(&self) -> &IntRes2dTransition {
        &self.trans2
    }
}

/// `IntRes2d_IntersectionSegment` (`IntRes2d_IntersectionSegment.hxx:29-88`).
#[derive(Clone, Copy, Debug)]
pub struct IntRes2dIntersectionSegment {
    /// `oppos` (`hxx:82`)
    oppos: bool,
    /// `first` (`hxx:83`)
    first: bool,
    /// `last` (`hxx:84`)
    last: bool,
    /// `ptfirst` (`hxx:85`)
    ptfirst: IntRes2dIntersectionPoint,
    /// `ptlast` (`hxx:86`)
    ptlast: IntRes2dIntersectionPoint,
}

impl Default for IntRes2dIntersectionSegment {
    /// `IntRes2d_IntersectionSegment()` (`IntersectionSegment.cxx:19-24`).
    fn default() -> Self {
        Self::new()
    }
}

impl IntRes2dIntersectionSegment {
    /// `IntRes2d_IntersectionSegment()` (`IntersectionSegment.cxx:19-24`).
    pub fn new() -> Self {
        Self {
            oppos: false,
            first: false,
            last: false,
            ptfirst: IntRes2dIntersectionPoint::new(),
            ptlast: IntRes2dIntersectionPoint::new(),
        }
    }

    /// `IntRes2d_IntersectionSegment(P1, P2, Oppos, ReverseFlag)`
    /// (`IntersectionSegment.lxx:17-37`).
    pub fn from_two_points(
        p1: &IntRes2dIntersectionPoint,
        p2: &IntRes2dIntersectionPoint,
        oppos: bool,
        reverse_flag: bool,
    ) -> Self {
        let mut r = Self {
            oppos,
            first: true,
            last: true,
            ptfirst: *p1,
            ptlast: *p2,
        };
        if reverse_flag && oppos {
            r.ptfirst = *p2;
            r.ptlast = *p1;
        }
        r
    }

    /// `IntRes2d_IntersectionSegment(P, First, Oppos, ReverseFlag)`
    /// (`IntersectionSegment.lxx:41-79`).
    pub fn from_one_point(
        p: &IntRes2dIntersectionPoint,
        first: bool,
        oppos: bool,
        reverse_flag: bool,
    ) -> Self {
        let mut r = Self {
            oppos,
            first: false,
            last: false,
            ptfirst: IntRes2dIntersectionPoint::new(),
            ptlast: IntRes2dIntersectionPoint::new(),
        };
        if reverse_flag && oppos {
            if first {
                r.first = false;
                r.last = true;
                r.ptlast = *p;
            } else {
                r.first = true;
                r.last = false;
                r.ptfirst = *p;
            }
        } else if first {
            r.first = true;
            r.last = false;
            r.ptfirst = *p;
        } else {
            r.first = false;
            r.last = true;
            r.ptlast = *p;
        }
        r
    }

    /// `IntRes2d_IntersectionSegment(Oppos)` (`IntersectionSegment.lxx:81-87`):
    /// an infinite segment of intersection.
    pub fn infinite(oppos: bool) -> Self {
        Self {
            oppos,
            first: false,
            last: false,
            ptfirst: IntRes2dIntersectionPoint::new(),
            ptlast: IntRes2dIntersectionPoint::new(),
        }
    }

    /// `IsOpposite()` (`IntersectionSegment.lxx:90-93`).
    pub fn is_opposite(&self) -> bool {
        self.oppos
    }

    /// `HasFirstPoint()` (`IntersectionSegment.lxx:95-98`).
    pub fn has_first_point(&self) -> bool {
        self.first
    }

    /// `FirstPoint()` (`IntersectionSegment.lxx:104-112`).
    pub fn first_point(&self) -> &IntRes2dIntersectionPoint {
        if !self.first {
            panic!("IntRes2d_IntersectionSegment::FirstPoint: no first point");
        }
        &self.ptfirst
    }

    /// `HasLastPoint()` (`IntersectionSegment.lxx:100-103`).
    pub fn has_last_point(&self) -> bool {
        self.last
    }

    /// `LastPoint()` (`IntersectionSegment.lxx:114-123`).
    pub fn last_point(&self) -> &IntRes2dIntersectionPoint {
        if !self.last {
            panic!("IntRes2d_IntersectionSegment::LastPoint: no last point");
        }
        &self.ptlast
    }
}

/// `TransitionEqual(T1, T2)` (`IntRes2d_Intersection.cxx:37-64`).
fn transition_equal(t1: &IntRes2dTransition, t2: &IntRes2dTransition) -> bool {
    if t1.position_on_curve() == t2.position_on_curve()
        && t1.transition_type() == t2.transition_type()
    {
        if t1.transition_type() == IntRes2dTypeTrans::Touch {
            if t1.is_tangent() == t2.is_tangent()
                && t1.situation() == t2.situation()
                && t1.is_opposite() == t2.is_opposite()
            {
                return true;
            }
        } else {
            return true;
        }
    }
    false
}

/// `InternalVerifyPosition` (`IntRes2d_Intersection.cxx:408-442`).
fn internal_verify_position(
    t1: &mut IntRes2dTransition,
    t2: &mut IntRes2dTransition,
    p_param_on_first: f64,
    p_param_on_second: f64,
    first_param1: f64,
    last_param1: f64,
    first_param2: f64,
    last_param2: f64,
) {
    if t1.position_on_curve() != IntRes2dPosition::Middle
        && !(paramequal(p_param_on_first, first_param1)
            || paramequal(p_param_on_first, last_param1))
        && p_param_on_first > first_param1
        && p_param_on_first < last_param1
    {
        t1.set_position(IntRes2dPosition::Middle);
    }
    if t2.position_on_curve() != IntRes2dPosition::Middle
        && !(paramequal(p_param_on_second, first_param2)
            || paramequal(p_param_on_second, last_param2))
        && p_param_on_second > first_param2
        && p_param_on_second < last_param2
    {
        t2.set_position(IntRes2dPosition::Middle);
    }
}

/// `IntRes2d_Intersection` (`IntRes2d_Intersection.hxx:33-110`).
///
/// OCCT's `IntRes2d_Intersection` is an abstract base whose protected members
/// `lpnt` / `lseg` / `done` / `reverse` are read and written directly by its
/// subclasses (`IntCurve_IntConicConic.cxx:176`, `:208`, `:224`). Rust has no
/// inheritance, so the port keeps the members public and lets a wrapper (the
/// future `Geom2dIntGInter`) contain one of these.
#[derive(Clone, Debug)]
pub struct IntRes2dIntersection {
    /// `lpnt` (`hxx:105`)
    pub lpnt: Vec<IntRes2dIntersectionPoint>,
    /// `lseg` (`hxx:106`)
    pub lseg: Vec<IntRes2dIntersectionSegment>,
    /// `done` (`hxx:107`)
    pub done: bool,
    /// `reverse` (`hxx:108`)
    pub reverse: bool,
}

impl Default for IntRes2dIntersection {
    /// `IntRes2d_Intersection()` (`Intersection.lxx:26-29`).
    fn default() -> Self {
        Self::new()
    }
}

impl IntRes2dIntersection {
    /// `IntRes2d_Intersection()` (`Intersection.lxx:26-29`).
    pub fn new() -> Self {
        Self {
            lpnt: Vec::new(),
            lseg: Vec::new(),
            done: false,
            reverse: false,
        }
    }

    /// `IsDone()` (`Intersection.lxx:20-23`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `IsEmpty()` (`Intersection.lxx:40-47`).
    pub fn is_empty(&self) -> bool {
        if !self.done {
            panic!("IntRes2d_Intersection::IsEmpty: not done");
        }
        self.lpnt.is_empty() && self.lseg.is_empty()
    }

    /// `NbPoints()` (`Intersection.lxx:49-56`).
    pub fn nb_points(&self) -> usize {
        if !self.done {
            panic!("IntRes2d_Intersection::NbPoints: not done");
        }
        self.lpnt.len()
    }

    /// `Point(N)` (`Intersection.lxx:58-65`). OCCT uses 1-based indexing.
    pub fn point(&self, n: usize) -> &IntRes2dIntersectionPoint {
        if !self.done {
            panic!("IntRes2d_Intersection::Point: not done");
        }
        &self.lpnt[n - 1]
    }

    /// `NbSegments()` (`Intersection.lxx:67-74`).
    pub fn nb_segments(&self) -> usize {
        if !self.done {
            panic!("IntRes2d_Intersection::NbSegments: not done");
        }
        self.lseg.len()
    }

    /// `Segment(N)` (`Intersection.lxx:76-83`). OCCT uses 1-based indexing.
    pub fn segment(&self, n: usize) -> &IntRes2dIntersectionSegment {
        if !self.done {
            panic!("IntRes2d_Intersection::Segment: not done");
        }
        &self.lseg[n - 1]
    }

    /// `Append(const IntRes2d_IntersectionSegment&)` (`Intersection.lxx:88-91`).
    pub fn append_segment(&mut self, seg: &IntRes2dIntersectionSegment) {
        self.lseg.push(*seg);
    }

    /// `Append(const IntRes2d_IntersectionPoint&)` (`Intersection.lxx:93-96`).
    pub fn append_point(&mut self, pnt: &IntRes2dIntersectionPoint) {
        self.lpnt.push(*pnt);
    }

    /// `ResetFields()` (`Intersection.lxx:102-109`).
    pub fn reset_fields(&mut self) {
        if self.done {
            self.lseg.clear();
            self.lpnt.clear();
            self.done = false;
        }
    }

    /// `SetReversedParameters(flag)` (`Intersection.lxx:112-115`).
    pub fn set_reversed_parameters(&mut self, flag: bool) {
        self.reverse = flag;
    }

    /// `ReversedParameters()` (`Intersection.lxx:117-120`).
    pub fn reversed_parameters(&self) -> bool {
        self.reverse
    }

    /// `Insert(const IntRes2d_IntersectionPoint&)` (`Intersection.cxx:66-112`).
    pub fn insert(&mut self, pnt: &IntRes2dIntersectionPoint) {
        let n = self.lpnt.len();
        if n == 0 {
            self.lpnt.push(*pnt);
            return;
        }
        let u = pnt.param_on_first();
        let mut i = 1usize;
        let mut b = n + 1;
        while i <= n {
            let pnti = &self.lpnt[i - 1];
            let ui = pnti.param_on_first();
            if ui >= u {
                b = i;
                i = n;
            }
            if paramequal(ui, u)
                && paramequal(pnt.param_on_second(), pnti.param_on_second())
                && transition_equal(pnt.transition_of_first(), pnti.transition_of_first())
                && transition_equal(pnt.transition_of_second(), pnti.transition_of_second())
            {
                b = 0;
                i = n;
            }
            i += 1;
        }
        if b > n {
            self.lpnt.push(*pnt);
        } else if b > 0 {
            self.lpnt.insert(b - 1, *pnt);
        }
    }

    /// `SetValues(const IntRes2d_Intersection&)` (`Intersection.cxx:114-142`).
    pub fn set_values(&mut self, other: &IntRes2dIntersection) {
        if other.done {
            self.lseg.clear();
            self.lpnt.clear();
            self.lpnt.extend_from_slice(&other.lpnt);
            self.lseg.extend_from_slice(&other.lseg);
            self.done = true;
        } else {
            self.done = false;
        }
    }

    /// `Append(const IntRes2d_Intersection&, FirstParam1, LastParam1,
    /// FirstParam2, LastParam2)` (`Intersection.cxx:183-383`).
    ///
    /// Only used by composite-curve intersections; ported for completeness.
    pub fn append_intersection(
        &mut self,
        other: &IntRes2dIntersection,
        first_param1: f64,
        last_param1: f64,
        first_param2: f64,
        last_param2: f64,
    ) {
        if !other.done {
            self.done = false;
            return;
        }
        let mut seg_modif_p1_first = 0.0;
        let mut seg_modif_p1_second = 0.0;
        let mut seg_modif_p2_first = 0.0;
        let mut seg_modif_p2_second = 0.0;

        for p in &other.lpnt {
            let mut t1 = *p.transition_of_first();
            let mut t2 = *p.transition_of_second();
            let p_param_on_first = p.param_on_first();
            let p_param_on_second = p.param_on_second();
            internal_verify_position(
                &mut t1,
                &mut t2,
                p_param_on_first,
                p_param_on_second,
                first_param1,
                last_param1,
                first_param2,
                last_param2,
            );
            let ip = IntRes2dIntersectionPoint::with_transitions(
                p.value(),
                p_param_on_first,
                p_param_on_second,
                &t1,
                &t2,
                false,
            );
            self.insert(&ip);
        }

        for seg in &other.lseg {
            let p1 = seg.first_point();
            let p1_pparam_on_first = p1.param_on_first();
            let p1_pparam_on_second = p1.param_on_second();
            let mut p1_t1 = *p1.transition_of_first();
            let mut p1_t2 = *p1.transition_of_second();
            let p1_pt = *p1.value();

            internal_verify_position(
                &mut p1_t1,
                &mut p1_t2,
                p1_pparam_on_first,
                p1_pparam_on_second,
                first_param1,
                last_param1,
                first_param2,
                last_param2,
            );

            let p2 = seg.last_point();
            let p2_pparam_on_first = p2.param_on_first();
            let p2_pparam_on_second = p2.param_on_second();
            let mut p2_t1 = *p2.transition_of_first();
            let mut p2_t2 = *p2.transition_of_second();
            let p2_pt = *p2.value();

            let opposite = seg.is_opposite();

            internal_verify_position(
                &mut p2_t1,
                &mut p2_t2,
                p2_pparam_on_first,
                p2_pparam_on_second,
                first_param1,
                last_param1,
                first_param2,
                last_param2,
            );

            let mut not_yet_modified = true;
            let an = self.lseg.len();
            for j in 0..an {
                let an_p1 = *self.lseg[j].first_point();
                let an_p1_pparam_on_first = an_p1.param_on_first();
                let an_p1_pparam_on_second = an_p1.param_on_second();

                let an_p2 = *self.lseg[j].last_point();
                let an_p2_pparam_on_first = an_p2.param_on_first();
                let an_p2_pparam_on_second = an_p2.param_on_second();

                if opposite == self.lseg[j].is_opposite() {
                    if paramequal(p1_pparam_on_first, an_p2_pparam_on_first)
                        && paramequal(p1_pparam_on_second, an_p2_pparam_on_second)
                    {
                        not_yet_modified = false;
                        self.lseg[j] = IntRes2dIntersectionSegment::from_two_points(
                            &an_p1, p2, opposite, false,
                        );
                        seg_modif_p1_first = an_p1_pparam_on_first;
                        seg_modif_p1_second = an_p1_pparam_on_second;
                        seg_modif_p2_first = p2_pparam_on_first;
                        seg_modif_p2_second = p2_pparam_on_second;
                    } else if paramequal(p2_pparam_on_first, an_p1_pparam_on_first)
                        && paramequal(p2_pparam_on_second, an_p1_pparam_on_second)
                    {
                        not_yet_modified = false;
                        self.lseg[j] = IntRes2dIntersectionSegment::from_two_points(
                            p1, &an_p2, opposite, false,
                        );
                        seg_modif_p1_first = p1_pparam_on_first;
                        seg_modif_p1_second = p1_pparam_on_second;
                        seg_modif_p2_first = an_p2_pparam_on_first;
                        seg_modif_p2_second = an_p2_pparam_on_second;
                    }
                    if paramequal(p1_pparam_on_first, an_p1_pparam_on_first)
                        && paramequal(p1_pparam_on_second, an_p1_pparam_on_second)
                    {
                        not_yet_modified = false;
                        self.lseg[j] = IntRes2dIntersectionSegment::from_two_points(
                            &an_p2, p2, opposite, false,
                        );
                        seg_modif_p1_first = p2_pparam_on_first;
                        seg_modif_p1_second = p2_pparam_on_second;
                        seg_modif_p2_first = an_p2_pparam_on_first;
                        seg_modif_p2_second = an_p2_pparam_on_second;
                    } else if paramequal(p2_pparam_on_first, an_p2_pparam_on_first)
                        && paramequal(p2_pparam_on_second, an_p2_pparam_on_second)
                    {
                        not_yet_modified = false;
                        self.lseg[j] = IntRes2dIntersectionSegment::from_two_points(
                            p1, &an_p1, opposite, false,
                        );
                        seg_modif_p1_first = p1_pparam_on_first;
                        seg_modif_p1_second = p1_pparam_on_second;
                        seg_modif_p2_first = an_p1_pparam_on_first;
                        seg_modif_p2_second = an_p1_pparam_on_second;
                    }
                }
            }

            if not_yet_modified {
                let new_p1 = IntRes2dIntersectionPoint::with_transitions(
                    &p1_pt,
                    p1_pparam_on_first,
                    p1_pparam_on_second,
                    &p1_t1,
                    &p1_t2,
                    false,
                );
                let new_p2 = IntRes2dIntersectionPoint::with_transitions(
                    &p2_pt,
                    p2_pparam_on_first,
                    p2_pparam_on_second,
                    &p2_t1,
                    &p2_t2,
                    false,
                );
                let new_seg = IntRes2dIntersectionSegment::from_two_points(
                    &new_p1, &new_p2, opposite, false,
                );
                self.append_segment(&new_seg);
            } else {
                let mut rp = 0usize;
                let mut rnbpts = self.lpnt.len();
                while rp < rnbpts {
                    let pon_first = self.lpnt[rp].param_on_first();
                    let pon_second = self.lpnt[rp].param_on_second();

                    let in_first = (pon_first >= seg_modif_p1_first
                        && pon_first <= seg_modif_p2_first)
                        || (pon_first <= seg_modif_p1_first && pon_first >= seg_modif_p2_first);
                    let in_second = (pon_second >= seg_modif_p1_second
                        && pon_second <= seg_modif_p2_second)
                        || (pon_second <= seg_modif_p1_second && pon_second >= seg_modif_p2_second);
                    if in_first && in_second {
                        self.lpnt.remove(rp);
                        rnbpts -= 1;
                    } else {
                        rp += 1;
                    }
                }
            }
        }

        self.done = true;
    }
}
