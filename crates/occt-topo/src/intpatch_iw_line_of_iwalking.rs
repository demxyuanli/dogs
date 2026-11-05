//! `IntPatch_TheIWLineOfTheIWalking` - one walking polyline of
//! `IntPatch_TheIWalking`.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/IntPatch/IntPatch_TheIWLineOfTheIWalking.hxx:37-315`,
//! `IntPatch_TheIWLineOfTheIWalking_0.cxx:28-51`.
//!
//! Members (`hxx:140-156`) and their port:
//!
//! | OCCT | port |
//! |------|------|
//! | `line` (`IntSurf_LineOn2S`) | `points` (`Vec<PntOn2S>`, 0-based) |
//! | `couple` (`NCollection_Sequence<IntSurf_Couple>`) | `couple` (`Vec<(i32,i32)>`) |
//! | `closed` / `hasFirst` / `hasLast` | `closed` / `has_first` / `has_last` |
//! | `firstIndex` / `lastIndex` | `first_index` / `last_index` |
//! | `theFirstPoint` / `theLastPoint` | not stored, see below |
//! | `indextg` / `vcttg` | `indextg` / `vcttg` |
//! | `istgtbeg` / `istgtend` | `tgt_begin` / `tgt_end` |
//!
//! `theFirstPoint` / `theLastPoint` are `IntSurf_PathPoint` values written by
//! `AddStatusFirst(Closed, HasFirst, Index, P)` (`hxx:181-190`) and
//! `AddStatusLast(HasLast, Index, P)` (`hxx:192-199`) and read back only by
//! `FirstPoint()` (`hxx:247-252`) and `LastPoint()` (`hxx:259-264`). Neither
//! reader has a caller in OCCT (`IntPatch_ImpPrmIntersection.cxx:1080-1290`
//! and `Contap_Contour.cxx:1726-1784` read `HasFirstPoint` / `FirstPointIndex`
//! / `HasLastPoint` / `LastPointIndex` only), and `IntSurf_PathPoint` is not
//! ported, so the two payloads are **UNPORTED** while the flags and indices
//! are kept.
//!
//! `NbPassingPoint()` (`hxx:273-276`) and `PassingPoint()` (`hxx:278-284`)
//! likewise have no caller in OCCT: `AddIndexPassing` (`hxx:210-213`) is
//! written at `IntWalk_IWalking.gxx:1747` and `:2531` and never read.
//!
//! Members this port keeps without an in-crate caller:
//!
//! * `cut` (`hxx:160-163`) - OCCT calls it only from
//!   `IntWalk_IWalking::TestArretCadre` (`IntWalk_IWalking.gxx:1234-1275`),
//!   whose frame-crossing trim is not ported in `intpatch_iwalking.rs`;
//! * `is_closed` (`hxx:230-233`), `add_index_passing` (`hxx:210-213`),
//!   `nb_passing_point` (`hxx:273-276`), `passing_point` (`hxx:278-284`) -
//!   no reader in OCCT;
//! * `value_mut` - `IntPatch_ImpPrmIntersection.cxx:933-940` writes the UV of
//!   `Value(1)` in place; the port rebuilds the point array in
//!   `intpatch_iwalking.rs::lines_to_wlines`, so the mutable accessor has no
//!   caller yet.
//!
//! `FirstPoint()` (`hxx:247-252`) and `LastPoint()` (`hxx:259-264`) are not
//! ported for the reason given above (`IntSurf_PathPoint` payload).

use occt_core::gp::GpVec;

use crate::int_tools_wline::PntOn2S;

/// `IntPatch_TheIWLineOfTheIWalking` (`hxx:37-156`).
#[derive(Debug, Clone)]
pub struct IwLineOfIWalking {
    /// `line` (`hxx:141`): the `IntSurf_LineOn2S`.
    pub points: Vec<PntOn2S>,
    /// `couple` (`hxx:142`): `(index on the line, index in Pnts1)` pairs.
    pub couple: Vec<(i32, i32)>,
    /// `closed` (`hxx:143`).
    pub closed: bool,
    /// `hasFirst` (`hxx:144`).
    pub has_first: bool,
    /// `hasLast` (`hxx:145`).
    pub has_last: bool,
    /// `firstIndex` (`hxx:146`).
    pub first_index: i32,
    /// `lastIndex` (`hxx:147`).
    pub last_index: i32,
    /// `indextg` (`hxx:150`).
    pub indextg: i32,
    /// `vcttg` (`hxx:151`).
    pub vcttg: GpVec,
    /// `istgtbeg` (`hxx:152`).
    pub tgt_begin: bool,
    /// `istgtend` (`hxx:153`).
    pub tgt_end: bool,
}

impl Default for IwLineOfIWalking {
    fn default() -> Self {
        Self::new()
    }
}

#[allow(dead_code)] // OCCT member surface: see the module notes for the
                     // members without an in-crate caller.
impl IwLineOfIWalking {
    /// `IntPatch_TheIWLineOfTheIWalking(allocator)`
    /// (`_0.cxx:28-40`): every flag `false`, `firstIndex` / `lastIndex` /
    /// `indextg` at `-1`, the tangent vector default constructed.
    pub fn new() -> Self {
        Self {
            points: Vec::new(),
            couple: Vec::new(),
            closed: false,
            has_first: false,
            has_last: false,
            first_index: -1,
            last_index: -1,
            indextg: -1,
            vcttg: GpVec::default(),
            tgt_begin: false,
            tgt_end: false,
        }
    }

    /// `Reverse()` (`_0.cxx:42-51`): reverse the line, then mirror the line
    /// index of every passing point. `Hasfirst` / `HasLast` are kept, and
    /// `indextg` is deliberately untouched (OCCT does not rewrite it either).
    pub fn reverse(&mut self) {
        self.points.reverse();
        let n = self.points.len() as i32;
        for c in self.couple.iter_mut() {
            c.0 = n - c.0 + 1;
        }
    }

    /// `Cut(Index)` (`hxx:160-163`): `line->Split(Index)` keeps points
    /// `1..Index-1` on the line (`IntSurf_LineOn2S.cxx:105-117` moves `Index..N`
    /// into the split-off line, which OCCT drops).
    pub fn cut(&mut self, index: i32) {
        if index >= 1 {
            self.points.truncate((index - 1) as usize);
        }
    }

    /// `AddPoint(P)` (`hxx:165-168`).
    pub fn add_point(&mut self, p: PntOn2S) {
        self.points.push(p);
    }

    /// `AddStatusFirst(Closed, HasFirst)` (`hxx:170-174`).
    pub fn add_status_first(&mut self, closed: bool, has_first: bool) {
        self.closed = closed;
        self.has_first = has_first;
    }

    /// `AddStatusFirst(Closed, HasFirst, Index, P)` (`hxx:181-190`). The
    /// `IntSurf_PathPoint` payload is not stored (see the module notes).
    pub fn add_status_first_at(&mut self, closed: bool, has_first: bool, index: i32) {
        self.closed = closed;
        self.has_first = has_first;
        self.first_index = index;
    }

    /// `AddStatusFirstLast(Closed, HasFirst, HasLast)` (`hxx:201-208`).
    pub fn add_status_first_last(&mut self, closed: bool, has_first: bool, has_last: bool) {
        self.closed = closed;
        self.has_first = has_first;
        self.has_last = has_last;
    }

    /// `AddStatusLast(HasLast)` (`hxx:176-179`).
    pub fn add_status_last(&mut self, has_last: bool) {
        self.has_last = has_last;
    }

    /// `AddStatusLast(HasLast, Index, P)` (`hxx:192-199`). The
    /// `IntSurf_PathPoint` payload is not stored (see the module notes).
    pub fn add_status_last_at(&mut self, has_last: bool, index: i32) {
        self.has_last = has_last;
        self.last_index = index;
    }

    /// `AddIndexPassing(Index)` (`hxx:210-213`): `couple.Append(
    /// IntSurf_Couple(line->NbPoints() + 1, Index))`.
    pub fn add_index_passing(&mut self, index: i32) {
        self.couple.push((self.points.len() as i32 + 1, index));
    }

    /// `SetTangentVector(V, Index)` (`hxx:286-290`).
    pub fn set_tangent_vector(&mut self, v: GpVec, index: i32) {
        self.indextg = index;
        self.vcttg = v;
    }

    /// `SetTangencyAtBegining(IsTangent)` (`hxx:292-295`).
    pub fn set_tangency_at_beginning(&mut self, is_tangent: bool) {
        self.tgt_begin = is_tangent;
    }

    /// `SetTangencyAtEnd(IsTangent)` (`hxx:297-300`).
    pub fn set_tangency_at_end(&mut self, is_tangent: bool) {
        self.tgt_end = is_tangent;
    }

    /// `NbPoints()` (`hxx:215-218`).
    pub fn nb_points(&self) -> i32 {
        self.points.len() as i32
    }

    /// `Value(Index)` (`hxx:220-223`): `Standard_OutOfRange` when the index is
    /// outside `1..NbPoints()`.
    pub fn value(&self, index: i32) -> Result<&PntOn2S, String> {
        if index < 1 || index > self.nb_points() {
            return Err("IntPatch_TheIWLineOfTheIWalking::Value - out of range".to_string());
        }
        Ok(&self.points[(index - 1) as usize])
    }

    /// `Value(Index)` for mutation, used by the `SetUV` walks of
    /// `IntPatch_ImpPrmIntersection.cxx:933-940`.
    pub fn value_mut(&mut self, index: i32) -> Result<&mut PntOn2S, String> {
        if index < 1 || index > self.nb_points() {
            return Err("IntPatch_TheIWLineOfTheIWalking::Value - out of range".to_string());
        }
        Ok(&mut self.points[(index - 1) as usize])
    }

    /// `Line()` (`hxx:225-228`).
    pub fn line(&self) -> &[PntOn2S] {
        &self.points
    }

    /// `IsClosed()` (`hxx:230-233`).
    pub fn is_closed(&self) -> bool {
        self.closed
    }

    /// `HasFirstPoint()` (`hxx:235-238`).
    pub fn has_first_point(&self) -> bool {
        self.has_first
    }

    /// `HasLastPoint()` (`hxx:254-257`).
    pub fn has_last_point(&self) -> bool {
        self.has_last
    }

    /// `FirstPointIndex()` (`hxx:240-245`): `Standard_DomainError` when the
    /// first point is not a marching point.
    pub fn first_point_index(&self) -> Result<i32, String> {
        if !self.has_first {
            return Err("IntPatch_TheIWLineOfTheIWalking::FirstPointIndex - no first point".into());
        }
        Ok(self.first_index)
    }

    /// `LastPointIndex()` (`hxx:266-271`).
    pub fn last_point_index(&self) -> Result<i32, String> {
        if !self.has_last {
            return Err("IntPatch_TheIWLineOfTheIWalking::LastPointIndex - no last point".into());
        }
        Ok(self.last_index)
    }

    /// `NbPassingPoint()` (`hxx:273-276`).
    pub fn nb_passing_point(&self) -> i32 {
        self.couple.len() as i32
    }

    /// `PassingPoint(Index, IndexLine, IndexPnts)` (`hxx:278-284`).
    pub fn passing_point(&self, index: i32) -> Result<(i32, i32), String> {
        if index < 1 || index > self.nb_passing_point() {
            return Err("IntPatch_TheIWLineOfTheIWalking::PassingPoint - out of range".into());
        }
        Ok(self.couple[(index - 1) as usize])
    }

    /// `TangentVector(int& Index)` (`hxx:302-306`): returns the stored vector
    /// and the index it belongs to.
    pub fn tangent_vector(&self) -> (GpVec, i32) {
        (self.vcttg, self.indextg)
    }

    /// `IsTangentAtBegining()` (`hxx:308-311`).
    pub fn is_tangent_at_beginning(&self) -> bool {
        self.tgt_begin
    }

    /// `IsTangentAtEnd()` (`hxx:313-316`).
    pub fn is_tangent_at_end(&self) -> bool {
        self.tgt_end
    }
}
