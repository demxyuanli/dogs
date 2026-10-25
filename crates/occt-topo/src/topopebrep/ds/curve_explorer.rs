//! `TopOpeBRepDS_CurveExplorer` — iterator over the curves of a
//! `TopOpeBRepDS_DataStructure`, optionally skipping the ones whose `Keep`
//! flag is false.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_CurveExplorer.hxx:29-76`,
//! `TopOpeBRepDS_CurveExplorer.cxx:23-158`.
//!
//! OCCT stores the structure as `void* myDS` and reaches into the private
//! `myCurves` map (`friend class`, `TopOpeBRepDS_DataStructure.hxx:312`).
//! The Rust explorer holds a borrow and reads through
//! [`DataStructure::is_curve_bound`], the same friendship.
//!
//! `Curve()` returns `*CEX_PEMPTY` (a lazily created empty curve) when the
//! explorer is exhausted; that case is `None` here.

use super::curve::Curve;
use super::data_structure::DataStructure;

/// Source: `TopOpeBRepDS_CurveExplorer.hxx:52-60`.
#[derive(Clone)]
pub struct CurveExplorer<'a> {
    /// `myIndex`
    index: i32,
    /// `myMax`
    max: i32,
    /// `myDS` — `None` is OCCT's null pointer.
    ds: Option<&'a DataStructure>,
    /// `myFound`
    found: bool,
    /// `myFindKeep`
    find_keep: bool,
}

impl<'a> Default for CurveExplorer<'a> {
    /// `TopOpeBRepDS_CurveExplorer()` (`cxx:23-30`): `myIndex = 1`, the rest
    /// zero/false.
    fn default() -> Self {
        Self {
            index: 1,
            max: 0,
            ds: None,
            found: false,
            find_keep: false,
        }
    }
}

impl<'a> CurveExplorer<'a> {
    /// `TopOpeBRepDS_CurveExplorer()` (`cxx:23-30`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `TopOpeBRepDS_CurveExplorer(DS, FindKeep)` (`cxx:34-39`).
    pub fn with_ds(ds: &'a DataStructure, find_keep: bool) -> Self {
        let mut e = Self {
            ds: Some(ds),
            ..Self::default()
        };
        e.init(ds, find_keep);
        e
    }

    /// `Init(DS, FindKeep)` (`cxx:43-50`).
    pub fn init(&mut self, ds: &'a DataStructure, find_keep: bool) {
        self.ds = Some(ds);
        self.index = 1;
        self.max = ds.nb_curves();
        self.find_keep = find_keep;
        self.find();
    }

    /// `Find()` (`cxx:54-72`).
    fn find(&mut self) {
        self.found = false;
        while self.index <= self.max {
            self.found = if self.find_keep {
                self.is_curve_keep(self.index)
            } else {
                self.is_curve(self.index)
            };
            if self.found {
                break;
            }
            self.index += 1;
        }
    }

    /// `More()` (`cxx:76-79`).
    pub fn more(&self) -> bool {
        self.found
    }

    /// `Next()` (`cxx:83-87`).
    pub fn next(&mut self) {
        self.index += 1;
        self.find();
    }

    /// `Curve()` (`cxx:91-104`): `*CEX_PEMPTY` when `!More()`.
    pub fn curve(&self) -> Option<&Curve> {
        if self.found {
            return Some(self.ds?.curve(self.index));
        }
        None
    }

    /// `IsCurve(I)` (`cxx:108-112`).
    pub fn is_curve(&self, i: i32) -> bool {
        match self.ds {
            Some(ds) => ds.is_curve_bound(i),
            None => false,
        }
    }

    /// `IsCurveKeep(I)` (`cxx:116-124`).
    pub fn is_curve_keep(&self, i: i32) -> bool {
        match self.ds {
            Some(ds) => ds.is_curve_bound(i) && ds.curve(i).keep(),
            None => false,
        }
    }

    /// `Curve(I)` (`cxx:128-140`): `*CEX_PEMPTY` when `I` is unbound.
    pub fn curve_at(&self, i: i32) -> Option<&Curve> {
        let ds = self.ds?;
        if ds.is_curve_bound(i) {
            return Some(ds.curve(i));
        }
        None
    }

    /// `NbCurve()` (`cxx:144-154`): resets the walk, counts, and leaves the
    /// explorer exhausted.
    pub fn nb_curve(&mut self) -> i32 {
        let Some(ds) = self.ds else {
            return 0;
        };
        self.index = 1;
        self.max = ds.nb_curves();
        self.find();
        let mut n = 0;
        while self.more() {
            n += 1;
            self.next();
        }
        n
    }

    /// `Index()` (`cxx:156-159`).
    pub fn index(&self) -> i32 {
        self.index
    }
}

