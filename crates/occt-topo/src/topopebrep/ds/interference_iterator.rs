//! `TopOpeBRepDS_InterferenceIterator` — the filtered interference walker the
//! other three DS iterators are built from.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_InterferenceIterator.hxx:26-82`,
//! `TopOpeBRepDS_InterferenceIterator.cxx:20-144`.
//!
//! `NCollection_List<handle(Interference)>::Iterator` is modelled as a slice
//! cursor: `More()` is `cursor < len`, `Value()` indexes it, and
//! `ChangeIterator()` becomes direct field access. The handles are
//! `Arc<RwLock<Interference>>` in this port, so every read of an interference
//! goes through [`InterferenceIterator::value`].


use super::interference::{Interference, InterferenceHandle};
use super::kind::Kind;

/// The `MatchInterference` override a derived iterator installs
/// (`TopOpeBRepDS_CurveIterator.cxx:34-39`,
/// `TopOpeBRepDS_PointIterator.cxx:33-39`).
///
/// OCCT makes `MatchInterference` virtual, so a derived override *replaces* the
/// base flag test instead of ANDing with it. Both overrides below only look at
/// `GeometryType()`, and no caller of these two subclasses ever sets the base
/// flags, so modelling the override as an extra predicate is observably the
/// same. `TopOpeBRepDS_SurfaceIterator` has no override — it drives the base
/// flag path.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DerivedMatch {
    /// `TopOpeBRepDS_CurveIterator`: `GT == TopOpeBRepDS_CURVE`.
    Curve,
    /// `TopOpeBRepDS_PointIterator`: `GT == POINT || GT == VERTEX`.
    Point,
}

/// A list of interferences (`TopOpeBRepDS_ListOfInterference`).
pub type InterferenceList = Vec<InterferenceHandle>;

/// Source: `TopOpeBRepDS_InterferenceIterator.hxx:26-82`.
#[derive(Debug, Clone)]
pub struct InterferenceIterator {
    list: InterferenceList,
    cursor: usize,
    /// `myGKDef` / `myGK`
    geometry_kind_def: bool,
    geometry_kind: Kind,
    /// `myGDef` / `myG`
    geometry_def: bool,
    geometry: i32,
    /// `mySKDef` / `mySK`
    support_kind_def: bool,
    support_kind: Kind,
    /// `mySDef` / `myS`
    support_def: bool,
    support: i32,
    derived: Option<DerivedMatch>,
}

impl Default for InterferenceIterator {
    /// `TopOpeBRepDS_InterferenceIterator()` (`cxx:22-28`): no filter is set
    /// and the cursor is uninitialised.
    fn default() -> Self {
        Self {
            list: InterferenceList::new(),
            cursor: 0,
            geometry_kind_def: false,
            geometry_kind: Kind::Unknown,
            geometry_def: false,
            geometry: 0,
            support_kind_def: false,
            support_kind: Kind::Unknown,
            support_def: false,
            support: 0,
            derived: None,
        }
    }
}

impl InterferenceIterator {
    /// `TopOpeBRepDS_InterferenceIterator()`.
    pub fn new() -> Self {
        Self::default()
    }

    /// `TopOpeBRepDS_InterferenceIterator(const List&)` (`cxx:32-40`).
    pub fn of_list(l: InterferenceList) -> Self {
        let mut it = Self::default();
        it.init(l);
        it
    }

    /// `Init(L)` (`cxx:43-48`): `myIterator.Initialize(L)` then `Match()`.
    pub fn init(&mut self, l: InterferenceList) {
        self.list = l;
        self.cursor = 0;
        self.match_filter();
    }

    /// `GeometryKind(GK)` (`cxx:52-56`).
    pub fn geometry_kind(&mut self, gk: Kind) {
        self.geometry_kind_def = true;
        self.geometry_kind = gk;
    }

    /// `Geometry(G)` (`cxx:60-64`).
    pub fn geometry_filter(&mut self, g: i32) {
        self.geometry_def = true;
        self.geometry = g;
    }

    /// `SupportKind(ST)` (`cxx:68-72`).
    pub fn support_kind(&mut self, st: Kind) {
        self.support_kind_def = true;
        self.support_kind = st;
    }

    /// `Support(S)` (`cxx:76-80`).
    pub fn support_filter(&mut self, s: i32) {
        self.support_def = true;
        self.support = s;
    }

    /// `Match()` (`cxx:84-99`): advance past non-matching interferences until
    /// one matches or the list is exhausted.
    pub fn match_filter(&mut self) {
        while self.cursor < self.list.len() {
            let m = {
                let i = self.list[self.cursor].read().expect("interference lock");
                self.match_interference(&i)
            };
            if !m {
                self.cursor += 1;
            } else {
                break;
            }
        }
    }

    /// `MatchInterference(I)` (`cxx:103-111`): every `*Def` flag gates one
    /// equality test.
    ///
    /// A derived override ([`DerivedMatch`]) replaces this whole test, which is
    /// what the virtual call in OCCT does.
    pub fn match_interference(&self, i: &Interference) -> bool {
        if let Some(d) = self.derived {
            return match d {
                DerivedMatch::Curve => i.geometry_type() == Kind::Curve,
                DerivedMatch::Point => {
                    i.geometry_type() == Kind::Point || i.geometry_type() == Kind::Vertex
                }
            };
        }
        let gk = !self.geometry_kind_def || i.geometry_type() == self.geometry_kind;
        let sk = !self.support_kind_def || i.support_type() == self.support_kind;
        let g = !self.geometry_def || i.geometry() == self.geometry;
        let s = !self.support_def || i.support() == self.support;
        gk && sk && g && s
    }

    /// Install the derived `MatchInterference` override.
    pub fn set_derived_match(&mut self, d: DerivedMatch) {
        self.derived = Some(d);
    }

    /// `More()` (`cxx:115-118`).
    pub fn more(&self) -> bool {
        self.cursor < self.list.len()
    }

    /// `Next()` (`cxx:122-129`).
    pub fn next(&mut self) {
        if self.cursor < self.list.len() {
            self.cursor += 1;
            self.match_filter();
        }
    }

    /// `Value()` (`cxx:133-136`). Panics once exhausted, like
    /// `NCollection_List::Iterator::Value()`.
    pub fn value(&self) -> &InterferenceHandle {
        &self.list[self.cursor]
    }

    /// A read guard on the current interference. OCCT hands out the handle and
    /// callers dereference it; every accessor of the four iterators reads
    /// through here.
    pub fn current(&self) -> std::sync::RwLockReadGuard<'_, Interference> {
        self.list[self.cursor].read().expect("interference lock")
    }

    /// `ChangeIterator()` (`cxx:140-144`).
    pub fn change_cursor(&mut self) -> &mut usize {
        &mut self.cursor
    }

    /// The list being walked.
    pub fn list(&self) -> &InterferenceList {
        &self.list
    }
}
