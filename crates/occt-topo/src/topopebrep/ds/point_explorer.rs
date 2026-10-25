//! `TopOpeBRepDS_PointExplorer` — iterator over the points of a
//! `TopOpeBRepDS_DataStructure`, optionally skipping the ones whose `Keep`
//! flag is false.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_PointExplorer.hxx:29-76`,
//! `TopOpeBRepDS_PointExplorer.cxx:22-152`.
//!
//! OCCT stores the structure as `void* myDS` and reaches into the private
//! `myPoints` map (`friend class`, `TopOpeBRepDS_DataStructure.hxx:313`).
//! The Rust explorer holds a borrow and reads through
//! [`DataStructure::is_point_bound`], the same friendship.
//!
//! `Point()` returns `myEmpty` when the explorer is exhausted; that case is
//! `None` here.

use super::data_structure::DataStructure;
use super::point::Point;

/// Source: `TopOpeBRepDS_PointExplorer.hxx:52-60`.
#[derive(Clone)]
pub struct PointExplorer<'a> {
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

impl<'a> Default for PointExplorer<'a> {
    /// `TopOpeBRepDS_PointExplorer()` (`cxx:22-29`): `myIndex = 1`, the rest
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

impl<'a> PointExplorer<'a> {
    /// `TopOpeBRepDS_PointExplorer()` (`cxx:22-29`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `TopOpeBRepDS_PointExplorer(DS, FindKeep)` (`cxx:33-38`).
    pub fn with_ds(ds: &'a DataStructure, find_keep: bool) -> Self {
        let mut e = Self {
            ds: Some(ds),
            ..Self::default()
        };
        e.init(ds, find_keep);
        e
    }

    /// `Init(DS, FindKeep)` (`cxx:42-50`).
    pub fn init(&mut self, ds: &'a DataStructure, find_keep: bool) {
        self.index = 1;
        self.max = ds.nb_points();
        self.ds = Some(ds);
        self.find_keep = find_keep;
        self.find();
    }

    /// `Find()` (`cxx:54-72`).
    fn find(&mut self) {
        self.found = false;
        while self.index <= self.max {
            self.found = if self.find_keep {
                self.is_point_keep(self.index)
            } else {
                self.is_point(self.index)
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

    /// `Point()` (`cxx:91-101`): `myEmpty` when `!More()`.
    pub fn point(&self) -> Option<&Point> {
        if self.found {
            return self.ds?.point(self.index);
        }
        None
    }

    /// `IsPoint(I)` (`cxx:105-109`).
    pub fn is_point(&self, i: i32) -> bool {
        match self.ds {
            Some(ds) => ds.is_point_bound(i),
            None => false,
        }
    }

    /// `IsPointKeep(I)` (`cxx:113-121`).
    pub fn is_point_keep(&self, i: i32) -> bool {
        match self.ds {
            Some(ds) => match ds.point(i) {
                Some(p) => p.keep(),
                None => false,
            },
            None => false,
        }
    }

    /// `Point(I)` (`cxx:125-135`): `myEmpty` when `I` is unbound.
    pub fn point_at(&self, i: i32) -> Option<&Point> {
        let ds = self.ds?;
        if ds.is_point_bound(i) {
            return ds.point(i);
        }
        None
    }

    /// `NbPoint()` (`cxx:139-149`): resets the walk, counts, and leaves the
    /// explorer exhausted.
    pub fn nb_point(&mut self) -> i32 {
        let Some(ds) = self.ds else {
            return 0;
        };
        self.index = 1;
        self.max = ds.nb_points();
        self.find();
        let mut n = 0;
        while self.more() {
            n += 1;
            self.next();
        }
        n
    }

    /// `Index()` (`cxx:151-154`).
    pub fn index(&self) -> i32 {
        self.index
    }
}

