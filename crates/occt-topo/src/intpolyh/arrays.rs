//! `IntPolyh_Array<Type>` (`IntPolyh_Array.hxx:55-219`) and the typed
//! instantiations used by the package (`IntPolyh_ArrayOf*.hxx`).
//!
//! The container mirrors OCCT exactly: `Init(n)` default-fills indices
//! `0..n` while leaving the already used elements alone and does **not**
//! touch the used count (`hxx:124-128`); `IncrementNbItems` bumps the used
//! count and re-allocates at `nb_allocated + increment` when it catches up
//! (`hxx:134-143`).

use occt_core::gp::GpPnt;
use occt_core::gp::GpVec;

use super::edge_couple::IntPolyhEdge;
use super::point::IntPolyhPoint;
use super::section_line::IntPolyhSectionLine;
use super::start_point::IntPolyhStartPoint;
use super::triangle::IntPolyhTriangle;

/// `IntPolyh_Array<Type>` (`hxx:55`).
#[derive(Debug, Clone)]
pub struct IntPolyhArray<T: Clone + Default> {
    data: Vec<T>,
    nb_items: i32,
    increment: i32,
}

impl<T: Clone + Default> Default for IntPolyhArray<T> {
    /// `IntPolyh_Array(aIncrement = 256)` (`hxx:69-74`).
    fn default() -> Self {
        Self {
            data: Vec::new(),
            nb_items: 0,
            increment: 256,
        }
    }
}

impl<T: Clone + Default> IntPolyhArray<T> {
    /// `IntPolyh_Array(aIncrement)` (`hxx:69`).
    pub fn with_increment(increment: i32) -> Self {
        Self {
            data: Vec::new(),
            nb_items: 0,
            increment,
        }
    }

    /// `Init(aN)` (`hxx:124-128`): allocate (at least) `aN + 1` slots, default
    /// constructing only the slots that did not exist yet.
    pub fn init(&mut self, n: i32) {
        if n < 0 {
            return;
        }
        while (self.data.len() as i32) <= n {
            self.data.push(T::default());
        }
    }

    /// `IncrementNbItems()` (`hxx:134-143`).
    pub fn increment_nb_items(&mut self) {
        self.nb_items += 1;
        if self.nb_items >= self.data.len() as i32 {
            let n = self.data.len() as i32 + self.increment;
            self.init(n);
        }
    }

    /// `GetN()` (`hxx:149`).
    pub fn get_n(&self) -> i32 {
        self.data.len() as i32
    }

    /// `NbItems()` (`hxx:155`).
    pub fn nb_items(&self) -> i32 {
        self.nb_items
    }

    /// `SetNbItems(aNb)` (`hxx:164`).
    pub fn set_nb_items(&mut self, nb: i32) {
        self.nb_items = nb;
    }

    /// `Value(aIndex)` (`hxx:175`). Indices past the allocation yield the
    /// default value instead of reading out of bounds; OCCT only ever reads
    /// allocated slots here.
    pub fn get(&self, index: i32) -> T {
        if index >= 0 && (index as usize) < self.data.len() {
            self.data[index as usize].clone()
        } else {
            T::default()
        }
    }

    /// `ChangeValue(aIndex)` (`hxx:190`). Grows the allocation when needed so
    /// the write goes through, matching `NCollection_DynamicArray::at` for
    /// indices inside the allocated range.
    pub fn get_mut(&mut self, index: i32) -> &mut T {
        while (self.data.len() as i32) <= index {
            self.data.push(T::default());
        }
        &mut self.data[index as usize]
    }

    /// Direct slice access for the read-only tail used by refinement helpers.
    pub fn as_slice(&self) -> &[T] {
        &self.data
    }

    /// `TArray[i]` on a const array, see [`Self::get`].
    pub fn at(&self, index: i32) -> T {
        self.get(index)
    }
}

/// `IntPolyh_ArrayOfPoints`.
pub type ArrayOfPoints = IntPolyhArray<IntPolyhPoint>;
/// `IntPolyh_ArrayOfEdges`.
pub type ArrayOfEdges = IntPolyhArray<IntPolyhEdge>;
/// `IntPolyh_ArrayOfTriangles`.
pub type ArrayOfTriangles = IntPolyhArray<IntPolyhTriangle>;
/// `IntPolyh_ArrayOfSectionLines`.
pub type ArrayOfSectionLines = IntPolyhArray<IntPolyhSectionLine>;
/// `IntPolyh_ArrayOfTangentZones`.
pub type ArrayOfTangentZones = IntPolyhArray<IntPolyhStartPoint>;

/// `IntPolyh_PointNormal` (`IntPolyh_ArrayOfPointNormal.hxx:24-47`).
#[derive(Debug, Clone)]
pub struct IntPolyhPointNormal {
    pub point: GpPnt,
    pub normal: GpVec,
}

impl Default for IntPolyhPointNormal {
    fn default() -> Self {
        Self {
            point: GpPnt::new(0.0, 0.0, 0.0),
            normal: GpVec::new(0.0, 0.0, 0.0),
        }
    }
}

/// `IntPolyh_ArrayOfPointNormal`.
pub type ArrayOfPointNormal = IntPolyhArray<IntPolyhPointNormal>;
