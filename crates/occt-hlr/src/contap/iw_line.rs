//! `Contap_TheIWLineOfTheIWalking` - the walking line of
//! `Contap_TheIWalking`, i.e. `IntWalk_TheIWLine` specialised to the contour
//! walking.
//!
//! Sources: `Contap_TheIWLineOfTheIWalking.hxx:38-158` (state and inline
//! accessors at `:160-320`), `Contap_TheIWLineOfTheIWalking_0.cxx:28-51`.
//!
//! The stored polyline is `IntSurf_LineOn2S` (`Add`, `Value`, `NbPoints`,
//! `Split`, `Reverse`); the port keeps it as a `Vec<PntOn2S>` and spells the
//! 1-based OCCT indices out at the boundary.

use occt_core::gp::{GpPnt, GpVec};

use occt_topo::int_tools_wline::PntOn2S;

use super::iwalking::types::PathPoint;

/// `Contap_TheIWLineOfTheIWalking` (`hxx:38-158`).
#[derive(Clone)]
pub struct ContapIWLineOfTheIWalking {
    /// `line` (`hxx:143`), `IntSurf_LineOn2S`.
    points: Vec<PntOn2S>,
    /// `couple` (`hxx:144`): `(IndexLine, IndexPnts)` of each passing point.
    couple: Vec<(i32, i32)>,
    /// `closed` (`hxx:145`).
    closed: bool,
    /// `hasFirst` (`hxx:146`).
    has_first: bool,
    /// `hasLast` (`hxx:147`).
    has_last: bool,
    /// `firstIndex` (`hxx:148`).
    first_index: i32,
    /// `lastIndex` (`hxx:149`).
    last_index: i32,
    /// `theFirstPoint` (`hxx:150`).
    the_first_point: Option<PathPoint>,
    /// `theLastPoint` (`hxx:151`).
    the_last_point: Option<PathPoint>,
    /// `indextg` (`hxx:152`).
    indextg: i32,
    /// `vcttg` (`hxx:153`).
    vcttg: GpVec,
    /// `istgtbeg` (`hxx:154`).
    istgtbeg: bool,
    /// `istgtend` (`hxx:155`).
    istgtend: bool,
}

impl ContapIWLineOfTheIWalking {
    /// `Contap_TheIWLineOfTheIWalking(theAllocator)` (`_0.cxx:28-38`):
    /// `closed = hasFirst = hasLast = istgtbeg = istgtend = false`,
    /// `firstIndex = lastIndex = indextg = -1`.
    pub fn new() -> Self {
        Self {
            points: Vec::new(),
            couple: Vec::new(),
            closed: false,
            has_first: false,
            has_last: false,
            first_index: -1,
            last_index: -1,
            the_first_point: None,
            the_last_point: None,
            indextg: -1,
            vcttg: GpVec::zero(),
            istgtbeg: false,
            istgtend: false,
        }
    }

    /// `Reverse()` (`_0.cxx:40-51`): reverse the points and remap every
    /// `couple` entry to `N - First + 1`.
    pub fn reverse(&mut self) {
        self.points.reverse();
        let n = self.points.len() as i32;
        for c in self.couple.iter_mut() {
            c.0 = n - c.0 + 1;
        }
    }

    /// `Cut(Index)` (`hxx:160-163`): `line->Split(Index)`, i.e. keep the
    /// points `1..=Index` and drop the rest.
    pub fn cut(&mut self, index: i32) {
        let keep = index.max(0) as usize;
        if keep < self.points.len() {
            self.points.truncate(keep);
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

    /// `AddStatusFirst(Closed, HasFirst, Index, P)` (`hxx:182-194`).
    pub fn add_status_first_with_point(
        &mut self,
        closed: bool,
        has_first: bool,
        index: i32,
        p: PathPoint,
    ) {
        self.closed = closed;
        self.has_first = has_first;
        self.first_index = index;
        self.the_first_point = Some(p);
    }

    /// `AddStatusLast(HasLast)` (`hxx:176-180`).
    pub fn add_status_last(&mut self, has_last: bool) {
        self.has_last = has_last;
    }

    /// `AddStatusLast(HasLast, Index, P)` (`hxx:196-203`).
    pub fn add_status_last_with_point(&mut self, has_last: bool, index: i32, p: PathPoint) {
        self.has_last = has_last;
        self.last_index = index;
        self.the_last_point = Some(p);
    }

    /// `AddStatusFirstLast(Closed, HasFirst, HasLast)` (`hxx:205-211`).
    pub fn add_status_first_last(&mut self, closed: bool, has_first: bool, has_last: bool) {
        self.closed = closed;
        self.has_first = has_first;
        self.has_last = has_last;
    }

    /// `AddIndexPassing(Index)` (`hxx:213-216`): binds the point that is
    /// about to be added (`NbPoints() + 1`) to the start point `Index`.
    pub fn add_index_passing(&mut self, index: i32) {
        let line_index = self.points.len() as i32 + 1;
        self.couple.push((line_index, index));
    }

    /// `SetTangentVector(V, Index)` (`hxx:288-292`).
    pub fn set_tangent_vector(&mut self, v: GpVec, index: i32) {
        self.indextg = index;
        self.vcttg = v;
    }

    /// `SetTangencyAtBegining(IsTangent)` (`hxx:294-297`): OCCT writes
    /// `istgtend` here, so the beginning flag is never set by this call. Kept
    /// as is.
    pub fn set_tangency_at_beginning(&mut self, is_tangent: bool) {
        self.istgtend = is_tangent;
    }

    /// `SetTangencyAtEnd(IsTangent)` (`hxx:299-302`).
    pub fn set_tangency_at_end(&mut self, is_tangent: bool) {
        self.istgtend = is_tangent;
    }

    /// `NbPoints()` (`hxx:218-221`).
    pub fn nb_points(&self) -> i32 {
        self.points.len() as i32
    }

    /// `Value(Index)` (`hxx:223-226`), 1-based; panics out of range as OCCT
    /// raises `Standard_OutOfRange`.
    pub fn value(&self, index: i32) -> &PntOn2S {
        &self.points[(index.max(1) - 1) as usize]
    }

    /// `Value(Index)` as a copy.
    pub fn point(&self, index: i32) -> PntOn2S {
        *self.value(index)
    }

    /// `Line()` (`hxx:228-231`).
    pub fn line(&self) -> &[PntOn2S] {
        &self.points
    }

    /// `IsClosed()` (`hxx:233-236`).
    pub fn is_closed(&self) -> bool {
        self.closed
    }

    /// `HasFirstPoint()` (`hxx:238-241`).
    pub fn has_first_point(&self) -> bool {
        self.has_first
    }

    /// `HasLastPoint()` (`hxx:256-259`).
    pub fn has_last_point(&self) -> bool {
        self.has_last
    }

    /// `FirstPoint()` (`hxx:249-254`); `false` where OCCT raises
    /// `Standard_DomainError`.
    pub fn first_point(&self) -> Option<&PathPoint> {
        self.the_first_point.as_ref()
    }

    /// `FirstPointIndex()` (`hxx:241-247`).
    pub fn first_point_index(&self) -> i32 {
        self.first_index
    }

    /// `LastPoint()` (`hxx:261-266`).
    pub fn last_point(&self) -> Option<&PathPoint> {
        self.the_last_point.as_ref()
    }

    /// `LastPointIndex()` (`hxx:268-272`).
    pub fn last_point_index(&self) -> i32 {
        self.last_index
    }

    /// `NbPassingPoint()` (`hxx:274-277`).
    pub fn nb_passing_point(&self) -> i32 {
        self.couple.len() as i32
    }

    /// `PassingPoint(Index, IndexLine, IndexPnts)` (`hxx:280-286`), 1-based.
    pub fn passing_point(&self, index: i32) -> (i32, i32) {
        self.couple[(index.max(1) - 1) as usize]
    }

    /// `TangentVector(Index)` (`hxx:304-308`).
    pub fn tangent_vector(&self) -> (&GpVec, i32) {
        (&self.vcttg, self.indextg)
    }

    /// `IsTangentAtBegining()` (`hxx:310-313`).
    pub fn is_tangent_at_beginning(&self) -> bool {
        self.istgtbeg
    }

    /// `IsTangentAtEnd()` (`hxx:315-318`).
    pub fn is_tangent_at_end(&self) -> bool {
        self.istgtend
    }

    /// The first stored 3D point (`Value(1).Value()`), `None` on an empty line.
    pub fn first_3d(&self) -> Option<GpPnt> {
        self.points.first().map(|p| p.p)
    }
}

impl Default for ContapIWLineOfTheIWalking {
    fn default() -> Self {
        Self::new()
    }
}
