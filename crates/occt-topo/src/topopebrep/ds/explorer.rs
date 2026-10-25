//! `TopOpeBRepDS_Explorer` — walk the shapes of a
//! `TopOpeBRepDS_HDataStructure`, optionally filtered by type and by the
//! `Keep` flag.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_Explorer.hxx:31-72`,
//! `.cxx:28-166`.
//!
//! | OCCT | lines |
//! |---|---|
//! | default constructor | `cxx:28-36` |
//! | `Explorer(HDS, T, findkeep)` | `cxx:40-45` |
//! | `Init` | `cxx:47-65` |
//! | `Type` | `cxx:69-72` |
//! | `Find` | `cxx:77-102` |
//! | `More` | `cxx:106-109` |
//! | `Next` | `cxx:113-117` |
//! | `Current` / `Index` | `cxx:121-135` |
//! | `Face` / `Edge` / `Vertex` | `cxx:139-166` |
//!
//! The walk is over the **shape** table of the data structure
//! (`TopOpeBRepDS_DataStructure::Shape(I, FindKeep)`, `cxx:87`), filtered by
//! `ShapeType() == myT || myT == TopAbs_SHAPE` (`cxx:89`).
//!
//! `Init` (`cxx:49-65`) writes `myFK = true` **before** the null-handle early
//! return and only overwrites it with the argument after `myN = NbShapes()`;
//! a null handle therefore leaves `myFK == true` and `myN == 0`.
//!
//! The `Next` / `Current` / `Index` / `Face` / `Edge` / `Vertex` raises
//! (`Standard_NoMoreObject`, `Standard_NoSuchObject`, `TopoDS`'s type
//! mismatch) become `None` / `false`, matching the convention used for the
//! other OCCT raises in this crate.

use crate::abs::ShapeType;
use crate::shape::{Edge, Face, TopoShape, Vertex};

use super::h_data_structure::HDataStructure;

/// `TopOpeBRepDS_Explorer` (`hxx:31-72`).
#[derive(Clone)]
pub struct Explorer<'a> {
    /// `myHDS` (`hxx:66`) — `None` is OCCT's null handle.
    hds: Option<&'a HDataStructure>,
    /// `myT` (`hxx:67`): `TopAbs_SHAPE` is [`ShapeType::Shape`].
    shape_type: ShapeType,
    /// `myI` (`hxx:68`).
    index: i32,
    /// `myN` (`hxx:69`).
    max: i32,
    /// `myB` (`hxx:70`).
    found: bool,
    /// `myFK` (`hxx:71`).
    find_keep: bool,
}

impl<'a> Default for Explorer<'a> {
    /// `TopOpeBRepDS_Explorer()` (`cxx:28-36`): `myT = TopAbs_SHAPE`,
    /// `myI = 1`, the rest zero/false.
    fn default() -> Self {
        Self {
            hds: None,
            shape_type: ShapeType::Shape,
            index: 1,
            max: 0,
            found: false,
            find_keep: true,
        }
    }
}

impl<'a> Explorer<'a> {
    /// `TopOpeBRepDS_Explorer()` (`cxx:28-36`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `TopOpeBRepDS_Explorer(HDS, T, findkeep)` (`cxx:40-45`).
    pub fn with_hds(hds: &'a HDataStructure, shape_type: ShapeType, find_keep: bool) -> Self {
        let mut e = Self::default();
        e.init(hds, shape_type, find_keep);
        e
    }

    /// `Init(HDS, T, FK)` (`cxx:47-65`).
    pub fn init(&mut self, hds: &'a HDataStructure, shape_type: ShapeType, find_keep: bool) {
        // `cxx:52-57`
        self.index = 1;
        self.max = 0;
        self.found = false;
        self.find_keep = true;
        self.shape_type = shape_type;
        self.hds = Some(hds);

        // `cxx:59-62`
        self.max = hds.nb_shapes();
        self.find_keep = find_keep;
        self.find();
    }

    /// `Type()` (`cxx:69-72`).
    pub fn shape_type(&self) -> ShapeType {
        self.shape_type
    }

    /// `Find()` (`cxx:77-102`).
    fn find(&mut self) {
        let Some(hds) = self.hds else {
            self.found = false;
            return;
        };
        let bds = hds.ds();

        // `cxx:79-101`
        let mut found = false;
        while self.index <= self.max && !found {
            // `cxx:83`
            if bds.keep_shape(self.index, self.find_keep) {
                // `cxx:87`
                match bds.shape(self.index, false) {
                    Some(s) => {
                        if s.shape_type() == self.shape_type
                            || self.shape_type == ShapeType::Shape
                        {
                            found = true;
                        } else {
                            self.index += 1;
                        }
                    }
                    // `cxx:87` raises `Standard_NoSuchObject` on an unbound
                    // index; the walk cannot leave the table, so skip.
                    None => self.index += 1,
                }
            } else {
                self.index += 1;
            }
        }
        self.found = found;
    }

    /// `More()` (`cxx:106-109`).
    pub fn more(&self) -> bool {
        self.found
    }

    /// `Next()` (`cxx:113-117`): `false` is the `Standard_NoMoreObject` raise.
    pub fn next(&mut self) -> bool {
        if !self.found {
            return false;
        }
        self.index += 1;
        self.find();
        true
    }

    /// `Current()` (`cxx:121-126`): `None` when `!More()`.
    pub fn current(&self) -> Option<TopoShape> {
        if !self.found {
            return None;
        }
        self.hds?.shape(self.index, false)
    }

    /// `Index()` (`cxx:130-135`): `None` when `!More()`.
    pub fn index(&self) -> Option<i32> {
        if !self.found {
            return None;
        }
        Some(self.index)
    }

    /// `Face()` (`cxx:139-144`): `None` when `!More()` or the current shape is
    /// not a face (`TopoDS::Face`'s `Standard_TypeMismatch`).
    pub fn face(&self) -> Option<Face> {
        let s = self.current()?;
        if s.shape_type() != ShapeType::Face {
            return None;
        }
        Some(Face(s))
    }

    /// `Edge()` (`cxx:148-153`).
    pub fn edge(&self) -> Option<Edge> {
        let s = self.current()?;
        if s.shape_type() != ShapeType::Edge {
            return None;
        }
        Some(Edge(s))
    }

    /// `Vertex()` (`cxx:157-162`).
    pub fn vertex(&self) -> Option<Vertex> {
        let s = self.current()?;
        if s.shape_type() != ShapeType::Vertex {
            return None;
        }
        Some(Vertex(s))
    }
}
