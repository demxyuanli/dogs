//! `TopOpeBRepDS_Marker` — the "mark an index" bit array the boolean engine
//! uses to remember which DS items it already visited.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_Marker.hxx:29-49`,
//! `TopOpeBRepDS_Marker.cxx:24-93`.
//!
//! `myhe` is an `HArray1<bool>` spanning `0..=myne`; index `0` is allocated but
//! never read (`GetI`/`Set` bound the index to `1..=myne`). The port keeps the
//! same `0..=myne` span in a `Vec<bool>` so the bounds match.

/// `TopOpeBRepDS_Marker` (`hxx:29-49`).
#[derive(Clone, Debug, Default)]
pub struct Marker {
    /// `myhe`
    he: Option<Vec<bool>>,
    /// `myne`
    ne: i32,
}

impl Marker {
    /// `TopOpeBRepDS_Marker()` (`cxx:26-29`).
    pub fn new() -> Self {
        let mut m = Marker {
            he: None,
            ne: 0,
        };
        m.reset();
        m
    }

    /// `Reset()` (`cxx:31-35`).
    pub fn reset(&mut self) {
        self.he = None;
        self.ne = 0;
    }

    /// `Set(ie, b)` (`cxx:37-45`): `Allocate(ie)` first, so asking for index
    /// `n` grows the array to `n` even when the mark itself is out of range.
    pub fn set(&mut self, ie: i32, b: bool) {
        self.allocate(ie);
        if ie < 1 || ie > self.ne {
            return;
        }
        if let Some(he) = self.he.as_mut() {
            he[ie as usize] = b;
        }
    }

    /// `Set(b, na, aa)` with `na == 0` (`cxx:47-50`): `myhe->Init(b)`.
    ///
    /// The `na > 0` arm (`cxx:53-63`) walks a `char**` of ASCII integers — a
    /// debug leftover with no caller in the package — and is UNPORTED.
    pub fn set_all(&mut self, b: bool) {
        if let Some(he) = self.he.as_mut() {
            for v in he.iter_mut() {
                *v = b;
            }
        }
    }

    /// `GetI(ie)` (`cxx:65-75`).
    pub fn get_i(&self, ie: i32) -> bool {
        let Some(he) = self.he.as_ref() else {
            return false;
        };
        if ie < 1 || ie > self.ne {
            return false;
        }
        he[ie as usize]
    }

    /// `Allocate(n)` (`cxx:77-93`).
    ///
    /// Faithful to OCCT, including the two quirks: a growth from a zero size
    /// allocates `1000` slots (`cxx:82-86`) and discards the existing marks,
    /// and a shrink still records `myne = n` (`cxx:88-91`) while the array
    /// keeps its larger extent.
    pub fn allocate(&mut self, n: i32) {
        let all = n > self.ne;
        let mut nall = n;
        if all {
            if self.ne == 0 {
                nall = 1000;
            }
            self.he = Some(vec![false; (nall + 1) as usize]);
        }
        if nall != 0 {
            self.ne = nall;
        }
    }

    /// `myne`.
    pub fn extent(&self) -> i32 {
        self.ne
    }
}
