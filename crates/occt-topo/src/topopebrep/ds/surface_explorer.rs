//! `TopOpeBRepDS_SurfaceExplorer` — iterator over the surfaces of a
//! `TopOpeBRepDS_DataStructure`, optionally skipping the ones whose `Keep`
//! flag is false.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_SurfaceExplorer.hxx:29-76`,
//! `TopOpeBRepDS_SurfaceExplorer.cxx:22-152`.
//!
//! OCCT stores the structure as `void* myDS` and reaches into the private
//! `mySurfaces` map (`friend class`, `TopOpeBRepDS_DataStructure.hxx:311`).
//! The Rust explorer holds a borrow and reads through
//! [`DataStructure::is_surface_bound`], the same friendship.
//!
//! `Surface()` returns `myEmpty` when the explorer is exhausted; that case is
//! `None` here.

use super::data_structure::DataStructure;
use super::surface::Surface;

/// Source: `TopOpeBRepDS_SurfaceExplorer.hxx:52-60`.
#[derive(Clone)]
pub struct SurfaceExplorer<'a> {
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

impl<'a> Default for SurfaceExplorer<'a> {
    /// `TopOpeBRepDS_SurfaceExplorer()` (`cxx:22-27`): `myIndex = 1`, the
    /// rest zero/false.
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

impl<'a> SurfaceExplorer<'a> {
    /// `TopOpeBRepDS_SurfaceExplorer()` (`cxx:22-27`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `TopOpeBRepDS_SurfaceExplorer(DS, FindKeep)` (`cxx:31-36`).
    pub fn with_ds(ds: &'a DataStructure, find_keep: bool) -> Self {
        let mut e = Self {
            ds: Some(ds),
            ..Self::default()
        };
        e.init(ds, find_keep);
        e
    }

    /// `Init(DS, FindKeep)` (`cxx:40-48`).
    pub fn init(&mut self, ds: &'a DataStructure, find_keep: bool) {
        self.index = 1;
        self.max = ds.nb_surfaces();
        self.ds = Some(ds);
        self.find_keep = find_keep;
        self.find();
    }

    /// `Find()` (`cxx:52-70`).
    fn find(&mut self) {
        self.found = false;
        while self.index <= self.max {
            self.found = if self.find_keep {
                self.is_surface_keep(self.index)
            } else {
                self.is_surface(self.index)
            };
            if self.found {
                break;
            }
            self.index += 1;
        }
    }

    /// `More()` (`cxx:74-77`).
    pub fn more(&self) -> bool {
        self.found
    }

    /// `Next()` (`cxx:81-85`).
    pub fn next(&mut self) {
        self.index += 1;
        self.find();
    }

    /// `Surface()` (`cxx:89-99`): `myEmpty` when `!More()`.
    pub fn surface(&self) -> Option<&Surface> {
        if self.found {
            return Some(self.ds?.surface(self.index));
        }
        None
    }

    /// `IsSurface(I)` (`cxx:103-107`).
    pub fn is_surface(&self, i: i32) -> bool {
        match self.ds {
            Some(ds) => ds.is_surface_bound(i),
            None => false,
        }
    }

    /// `IsSurfaceKeep(I)` (`cxx:111-119`).
    pub fn is_surface_keep(&self, i: i32) -> bool {
        match self.ds {
            Some(ds) => ds.is_surface_bound(i) && ds.surface(i).keep(),
            None => false,
        }
    }

    /// `Surface(I)` (`cxx:123-133`): `myEmpty` when `I` is unbound.
    pub fn surface_at(&self, i: i32) -> Option<&Surface> {
        let ds = self.ds?;
        if ds.is_surface_bound(i) {
            return Some(ds.surface(i));
        }
        None
    }

    /// `NbSurface()` (`cxx:137-147`): resets the walk, counts, and leaves the
    /// explorer exhausted.
    pub fn nb_surface(&mut self) -> i32 {
        let Some(ds) = self.ds else {
            return 0;
        };
        self.index = 1;
        self.max = ds.nb_surfaces();
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

