//! `TableauRejection` (`ModelingAlgorithms/TKHLR/HLRBRep/HLRBRep_Data.cxx:78-491`,
//! local class of `HLRBRep_Data.cxx`): a symmetric table of the single-point
//! 2D intersections already found between edge pairs, plus the bit table of the
//! pairs already proven disjoint.
//!
//! Two separate stores:
//!
//! * `UV` / `IndUV` / `nbUV` - for every edge `i` the list of
//!   `(j, u)` pairs (`u` = parameter of the intersection on `i`), kept sorted
//!   by **decreasing** `j` so `Get` can binary-search it
//!   (`cxx:243-378`). `nbUV(i)` is the **allocated** length of the row, not the
//!   number of filled slots; unused slots hold `IndUV = -1` (`cxx:126-145`).
//!   A row grows by `SIZEUV` slots when it is full (`cxx:259-291`).
//! * `TabBit` - a triangular bit matrix: `NoIntersection(i, j)` reads bit
//!   `j & 31` of word `j >> 5` of row `min(i, j)` (`cxx:421-458`). Rows are
//!   `1 + (nbEdges >> 5)` words long, i.e. 32 edges per `unsigned` word.
//!
//! `GetSingleIntersection` (`cxx:480-491`) answers `(RealLast(), RealLast())`
//! when the pair is unknown; `NextInterference` uses that as its "no cached
//! point" test. `SetIntersection` (`cxx:462-478`) only records the intersection
//! when both transitions are `Middle` and the first one is `In`/`Out` - a
//! tangent (`Touch`) or endpoint (`Head`/`End`) hit is not cached.
//!
//! The OCCT class is a raw `new`/`free` structure with a deleted copy
//! constructor (`cxx:95-96`) and is only reachable through a pointer
//! (`myReject`, `hxx:285`); the port owns a [`TableauRejection`] value and
//! derives `Clone`, which the `HLRBRep_Data` clone needs.

use occt_core::intres2d::{IntRes2dIntersectionPoint, IntRes2dPosition, IntRes2dTypeTrans};

/// `SIZEUV` (`cxx:76`): the row growth step of `UV` / `IndUV`.
const SIZEUV: usize = 8;

/// `RealLast()` (`Standard_Real.hxx:138-141`): `DBL_MAX`, the "no entry" value
/// `TableauRejection::Get` returns (`cxx:333`, `:346`, `:374`).
pub(crate) const REAL_LAST: f64 = f64::MAX;

/// `TableauRejection` (`cxx:78-491`).
#[derive(Debug, Clone, Default)]
pub(crate) struct TableauRejection {
    /// `double** UV` (`cxx:80`): per edge, the parameters of the cached
    /// intersections, sorted by decreasing partner index.
    uv: Vec<Vec<f64>>,
    /// `int** IndUV` (`cxx:82`): per edge, the partner indices of `uv`;
    /// `-1` marks a free slot.
    ind_uv: Vec<Vec<i32>>,
    /// `int* nbUV` (`cxx:83`): per edge, the **allocated** row length.
    nb_uv: Vec<usize>,
    /// `unsigned long** TabBit` (`cxx:88`): the "no intersection" bit matrix,
    /// 32 edges per word.
    tab_bit: Vec<Vec<u32>>,
    /// `int nTabBit` (`cxx:89`): the number of allocated `TabBit` rows.
    n_tab_bit: usize,
    /// `int N` (`cxx:84`): the number of edges.
    n: usize,
}

impl TableauRejection {
    /// `TableauRejection()` (`cxx:100-113`): the empty table.
    pub(crate) fn new() -> Self {
        Self::default()
    }

    /// `SetDim(N)` (`cxx:115-146`): (re)allocate for `N` edges and clear the
    /// bit matrix.
    pub(crate) fn set_dim(&mut self, n: usize) {
        if !self.uv.is_empty() {
            self.destroy();
        }
        self.n = n;
        self.uv = vec![vec![0.0; SIZEUV]; n];
        self.ind_uv = vec![vec![-1; SIZEUV]; n];
        self.nb_uv = vec![SIZEUV; n];
        self.init_tab_bit(n);
    }

    /// `Destroy()` (`cxx:155-242`): drop every row and the bit matrix.
    pub(crate) fn destroy(&mut self) {
        if self.n == 0 {
            return;
        }
        self.reset_tab_bit(self.n);
        self.ind_uv.clear();
        self.uv.clear();
        self.nb_uv.clear();
        self.n = 0;
    }

    /// `Set(i0, j0, u)` (`cxx:244-316`), OCCT 1-based.
    ///
    /// Writes `(j0, u)` into the first free slot of row `i0` (growing the row
    /// by `SIZEUV` when it is full) and then sorts the row by **decreasing**
    /// `j0` with the same repeated single-pass pass OCCT runs
    /// (`cxx:293-315`: the `do { TriOk = true; ... } while (!TriOk)` bubble).
    pub(crate) fn set(&mut self, i0: i32, j0: i32, u: f64) {
        let i0 = (i0 - 1).max(0) as usize;
        let j0 = j0 - 1;
        if i0 >= self.nb_uv.len() {
            return;
        }

        // `cxx:253-260`: the first slot whose `IndUV` is `-1`.
        let mut k = None;
        for i in 0..self.nb_uv[i0] {
            if self.ind_uv[i0][i] == -1 {
                k = Some(i);
                break;
            }
        }
        let k = match k {
            Some(k) => k,
            None => {
                // `cxx:261-291`: grow the row by `SIZEUV` slots.
                let old = self.nb_uv[i0];
                self.nb_uv[i0] = old + SIZEUV;
                self.uv[i0].resize(self.nb_uv[i0], 0.0);
                self.ind_uv[i0].resize(self.nb_uv[i0], -1);
                for kk in old..self.nb_uv[i0] {
                    self.ind_uv[i0][kk] = -1;
                }
                old
            }
        };

        self.ind_uv[i0][k] = j0;
        self.uv[i0][k] = u;

        // `cxx:293-315`: sort by decreasing partner index. The pass is
        // repeated until no swap happens; a free slot (`-1`) ends the row.
        loop {
            let mut tri_ok = true;
            let mut im1 = 0usize;
            let mut i = 1usize;
            while i < self.nb_uv[i0] && self.ind_uv[i0][i] != -1 {
                if self.ind_uv[i0][i] > self.ind_uv[i0][im1] {
                    tri_ok = false;
                    self.ind_uv[i0].swap(i, im1);
                    self.uv[i0].swap(i, im1);
                }
                i += 1;
                im1 += 1;
            }
            if tri_ok {
                break;
            }
        }
    }

    /// `Get(i0, j0)` (`cxx:318-378`), OCCT 1-based.
    ///
    /// Binary search over the decreasing row; `RealLast()` when the pair is
    /// not cached. The loop condition is OCCT's
    /// `while ((IndUV[i0][a] > j0) && (IndUV[i0][b] < j0))` verbatim, including
    /// the two `a == ab` / `b == ab` early exits.
    pub(crate) fn get(&self, i0: i32, j0: i32) -> f64 {
        let i0 = (i0 - 1).max(0) as usize;
        let j0 = j0 - 1;
        if i0 >= self.nb_uv.len() {
            return REAL_LAST;
        }
        let mut a = 0usize;
        let mut b = self.nb_uv[i0] - 1;
        if self.ind_uv[i0][a] == -1 {
            return REAL_LAST;
        }
        if self.ind_uv[i0][a] == j0 {
            return self.uv[i0][a];
        }
        if self.ind_uv[i0][b] == j0 {
            return self.uv[i0][b];
        }
        while self.ind_uv[i0][a] > j0 && self.ind_uv[i0][b] < j0 {
            let ab = (a + b) >> 1;
            if self.ind_uv[i0][ab] < j0 {
                if b == ab {
                    return REAL_LAST;
                } else {
                    b = ab;
                }
            } else if self.ind_uv[i0][ab] > j0 {
                if a == ab {
                    return REAL_LAST;
                } else {
                    a = ab;
                }
            } else {
                return self.uv[i0][ab];
            }
        }
        REAL_LAST
    }

    /// `ResetTabBit(nbedgs)` (`cxx:379-397`).
    fn reset_tab_bit(&mut self, _nbedgs: usize) {
        if !self.tab_bit.is_empty() {
            self.tab_bit.clear();
            self.n_tab_bit = 0;
        }
    }

    /// `InitTabBit(nbedgs)` (`cxx:399-419`): `1 + (nbedgs >> 5)` zero words per
    /// row.
    fn init_tab_bit(&mut self, nbedgs: usize) {
        if !self.tab_bit.is_empty() && self.n_tab_bit != 0 {
            self.reset_tab_bit(self.n_tab_bit);
        }
        let n = 1 + (nbedgs >> 5);
        self.tab_bit = vec![vec![0u32; n]; nbedgs];
        self.n_tab_bit = nbedgs;
    }

    /// `SetNoIntersection(i0, i1)` (`cxx:421-436`), OCCT 1-based: set bit
    /// `max(i0, i1)` of row `min(i0, i1)`.
    pub(crate) fn set_no_intersection(&mut self, i0: i32, i1: i32) {
        let (i0, i1) = (i0 - 1, i1 - 1);
        let (i0, i1) = if i0 > i1 { (i1, i0) } else { (i0, i1) };
        let (i0, i1) = (i0 as usize, i1 as usize);
        if i0 >= self.tab_bit.len() {
            return;
        }
        let (c, o) = (i1 >> 5, i1 & 31);
        if c < self.tab_bit[i0].len() {
            self.tab_bit[i0][c] |= 1u32 << o;
        }
    }

    /// `NoIntersection(i0, i1)` (`cxx:438-458`).
    pub(crate) fn no_intersection(&self, i0: i32, i1: i32) -> bool {
        let (i0, i1) = (i0 - 1, i1 - 1);
        let (i0, i1) = if i0 > i1 { (i1, i0) } else { (i0, i1) };
        let (i0, i1) = (i0 as usize, i1 as usize);
        if i0 >= self.tab_bit.len() {
            return false;
        }
        let (c, o) = (i1 >> 5, i1 & 31);
        c < self.tab_bit[i0].len() && (self.tab_bit[i0][c] & (1u32 << o)) != 0
    }

    /// `SetIntersection(i0, i1, IP)` (`cxx:460-478`): cache the point only for
    /// `Middle`/`Middle` transitions whose first one is `In` or `Out`.
    pub(crate) fn set_intersection(&mut self, i0: i32, i1: i32, ip: &IntRes2dIntersectionPoint) {
        let t1 = ip.transition_of_first();
        let t2 = ip.transition_of_second();
        if t1.position_on_curve() == IntRes2dPosition::Middle
            && t2.position_on_curve() == IntRes2dPosition::Middle
        {
            let tt = t1.transition_type();
            if tt == IntRes2dTypeTrans::In || tt == IntRes2dTypeTrans::Out {
                self.set(i0, i1, ip.param_on_first());
                self.set(i1, i0, ip.param_on_second());
            }
        }
    }

    /// `GetSingleIntersection(i0, i1, u, v)` (`cxx:480-491`).
    pub(crate) fn get_single_intersection(&self, i0: i32, i1: i32) -> (f64, f64) {
        let u = self.get(i0, i1);
        if u != REAL_LAST {
            (u, self.get(i1, i0))
        } else {
            (u, REAL_LAST)
        }
    }
}
