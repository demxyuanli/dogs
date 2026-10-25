//! `TopOpeBRep_ShapeScanner` — walks the boxes touched by a shape.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRep/TopOpeBRep_ShapeScanner.hxx:26-72`,
//! `.cxx:37-151`.
//!
//! Translation notes:
//! - `myListIterator` (`hxx:63`) is a `NCollection_List<int>::Iterator` into
//!   the list `myBoxSort.Compare` produced; the port forwards to
//!   [`BoxSort::more`] / [`BoxSort::next`] / [`BoxSort::current`] /
//!   [`BoxSort::index`], which keep that list and a cursor.
//! - `Init(TopOpeBRepTool_ShapeExplorer&)` (`cxx:57-74`) re-initialises the
//!   scanner on each explored shape and stops at the first one whose compare
//!   list is non-empty; an exhausted explorer leaves the empty list of
//!   `cxx:60-62`.
//! - `DumpCurrent` (`cxx:135-151`) writes only under `OCCT_DEBUG`; like the
//!   other `TopOpeBRep` debug dumps it is not ported.

use crate::shape::TopoShape;
use crate::topopebrep::tool::box_sort::BoxSort;
use crate::topopebrep::tool::shape_explorer::ShapeExplorer;

/// `class TopOpeBRep_ShapeScanner` (`hxx:26-72`).
#[derive(Default)]
pub struct ShapeScanner {
    /// `myBoxSort` (`hxx:62`).
    box_sort: BoxSort,
}

impl ShapeScanner {
    /// `TopOpeBRep_ShapeScanner()` (`cxx:37`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `Clear()` (`cxx:42-45`).
    pub fn clear(&mut self) {
        self.box_sort.clear();
    }

    /// `AddBoxesMakeCOB(S, TS, TA)` (`cxx:48-53`).
    pub fn add_boxes_make_cob(
        &mut self,
        s: &TopoShape,
        ts: crate::abs::ShapeType,
        ta: crate::abs::ShapeType,
    ) -> Result<(), String> {
        self.box_sort.add_boxes_make_cob(s, ts, ta)
    }

    /// `Init(TopOpeBRepTool_ShapeExplorer& E)` (`cxx:57-74`).
    pub fn init_from_explorer(&mut self, e: &mut ShapeExplorer) -> Result<(), String> {
        // myListIterator.Initialize(anEmptyList) (cxx:60-62).
        self.box_sort.reset_iterator();
        while e.more() {
            if let Some(cur) = e.current() {
                self.init(cur)?;
            }
            if self.more() {
                break;
            }
            e.next();
        }
        Ok(())
    }

    /// `Init(const TopoDS_Shape& E)` (`cxx:78-81`).
    pub fn init(&mut self, e: &TopoShape) -> Result<(), String> {
        self.box_sort.compare(e)
    }

    /// `More()` (`cxx:85-90`).
    pub fn more(&self) -> bool {
        self.box_sort.more()
    }

    /// `Next()` (`cxx:94-97`).
    pub fn next(&mut self) {
        self.box_sort.next();
    }

    /// `Current()` (`cxx:101-106`).
    pub fn current(&self) -> Option<Result<&TopoShape, String>> {
        self.box_sort.current()
    }

    /// `BoxSort()` (`cxx:109-112`).
    pub fn box_sort(&self) -> &BoxSort {
        &self.box_sort
    }

    /// `ChangeBoxSort()` (`cxx:116-119`).
    pub fn change_box_sort(&mut self) -> &mut BoxSort {
        &mut self.box_sort
    }

    /// `Index()` (`cxx:123-132`).
    pub fn index(&self) -> i32 {
        self.box_sort.index()
    }
}
