//! `FilletPoint` — the private helper of `ChFi2d_FilletAlgo`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi2d/ChFi2d_FilletAlgo.hxx:139-227`
//! and `.cxx:761-977` (OCCT 8.0.0).
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `FilletPoint(const double)` (`cxx:761-765`) | 5 | [`FilletPoint::new`] |
//! | `appendValue` (`cxx:767-781`) | 15 | [`FilletPoint::append_value`] |
//! | `calculateDiff` (`cxx:783-828`) | 46 | [`FilletPoint::calculate_diff`] |
//! | `FilterPoints` (`cxx:830-942`) | 113 | [`FilletPoint::filter_points`] |
//! | `Copy` (`cxx:944-955`) | 12 | [`FilletPoint::copy`] |
//! | `hasSolution` (`cxx:957-969`) | 13 | [`FilletPoint::has_solution`] |
//! | `remove` (`cxx:971-977`) | 7 | [`FilletPoint::remove`] |
//!
//! Translation notes:
//! - every `NCollection_Sequence` index is 1-based in OCCT and the accessors
//!   here keep that (`value(1)` is `myV.Value(1)`), but the indices are `f64`
//!   payloads in a `Vec` so the sentinel `0` of `myNear` becomes
//!   `Option<usize>`,
//! - `myNear` / `aDiffs` are the two parallel sequences of `FilterPoints`; OCCT
//!   keeps `aDiffs` local and `myNear` a member, and the port does the same,
//! - `LowerValue()` (`hxx:194-209`) is not ported: it has no reader and its
//!   `aValue` is only assigned inside a branch guarded by
//!   `aResultIndex == 0 || ...`, i.e. it reads an uninitialised `double` on the
//!   first `myV.Length() == 0` pass. It is `//! For debug only`.

use occt_core::gp::GpPnt2d;
use occt_core::precision::CONFUSION;

/// `FilletPoint` (`hxx:139-227`): a point on the first curve with the values,
/// derivatives and validity flags of the fillet function there.
#[derive(Clone, Default)]
pub(super) struct FilletPoint {
    /// `myParam` (`hxx:214`): parameter on the first curve.
    param: f64,
    /// `myParam2` (`hxx:216`): parameter on the second curve.
    param2: f64,
    /// `myCenter` (`hxx:224`): centre of the fillet arc.
    center: GpPnt2d,
    /// `myV` (`hxx:222`): values of the fillet function, ascending.
    v: Vec<f64>,
    /// `myD` (`hxx:222`): derivative values, parallel to `myV`.
    d: Vec<f64>,
    /// `myValid` (`hxx:227`): validity flags, parallel to `myV`.
    valid: Vec<bool>,
    /// `myNear` (`hxx:228`): the `thePoint` index each `myV` entry is paired
    /// with, parallel to `myV` after `FilterPoints`.
    near: Vec<usize>,
}

impl FilletPoint {
    /// `FilletPoint(const double theParam)` (`cxx:761-765`):
    /// `myParam = theParam`, `myParam2 = 0.0`, every sequence empty.
    pub fn new(param: f64) -> Self {
        Self {
            param,
            ..Self::default()
        }
    }

    /// `setParam` (`hxx:145`).
    pub fn set_param(&mut self, param: f64) {
        self.param = param;
    }

    /// `getParam` (`hxx:150`).
    pub fn param(&self) -> f64 {
        self.param
    }

    /// `getNBValues` (`hxx:153`).
    pub fn nb_values(&self) -> usize {
        self.v.len()
    }

    /// `getValue(int)` (`hxx:156`), `index` 1-based.
    pub fn value(&self, index: usize) -> f64 {
        self.v[index - 1]
    }

    /// `getDiff(int)` (`hxx:159`), `index` 1-based.
    pub fn diff(&self, index: usize) -> f64 {
        self.d[index - 1]
    }

    /// `isValid(int)` (`hxx:162`), `index` 1-based.
    pub fn is_valid(&self, index: usize) -> bool {
        self.valid[index - 1]
    }

    /// `getNear(int)` (`hxx:165`), `index` 1-based.
    pub fn near(&self, index: usize) -> usize {
        self.near[index - 1]
    }

    /// `setParam2` (`hxx:167`).
    pub fn set_param2(&mut self, param2: f64) {
        self.param2 = param2;
    }

    /// `getParam2` (`hxx:171`).
    pub fn param2(&self) -> f64 {
        self.param2
    }

    /// `setCenter` (`hxx:174`).
    pub fn set_center(&mut self, center: GpPnt2d) {
        self.center = center;
    }

    /// `getCenter` (`hxx:177`).
    pub fn center(&self) -> GpPnt2d {
        self.center
    }

    /// `appendValue(theValue, theValid)` (`cxx:767-781`): a sorted insert of
    /// `theValue` into `myV`, with `myValid` kept parallel. `myD` is **not**
    /// inserted into — OCCT leaves it for `calculateDiff` to fill.
    pub fn append_value(&mut self, value: f64, valid: bool) {
        for a in 1..=self.v.len() {
            if value < self.v[a - 1] {
                self.v.insert(a - 1, value);
                self.valid.insert(a - 1, valid);
                return;
            }
        }
        self.v.push(value);
        self.valid.push(valid);
    }

    /// `calculateDiff(thePoint)` (`cxx:783-828`): fills `myD` from the value
    /// difference against `thePoint`.
    ///
    /// `true` means "absolutely the same points" — the two value sequences have
    /// the same length, so entry `a` of one is entry `a` of the other. `false`
    /// takes the nearest analog instead.
    pub fn calculate_diff(&mut self, other: &FilletPoint) -> bool {
        let diffs_set = !self.d.is_empty();
        let dx = other.param - self.param;
        // `aDY` is declared before both loops in OCCT and is **not** reset per
        // outer step; the second loop therefore starts from the previous
        // entry's difference. Copied as is.
        let mut ady = 0.0;
        if other.v.len() == self.v.len() {
            for a in 1..=self.v.len() {
                ady = other.v[a - 1] - self.v[a - 1];
                if diffs_set {
                    self.d[a - 1] = ady / dx;
                } else {
                    self.d.push(ady / dx);
                }
            }
            return true;
        }
        for a in 1..=self.v.len() {
            for b in 1..=other.v.len() {
                if b == 1 || (other.v[b - 1] - self.v[a - 1]).abs() < ady.abs() {
                    ady = other.v[b - 1] - self.v[a - 1];
                }
            }
            if diffs_set {
                if (ady / dx).abs() < self.d[a - 1].abs() {
                    self.d[a - 1] = ady / dx;
                }
            } else {
                self.d.push(ady / dx);
            }
        }
        false
    }

    /// `FilterPoints(thePoint)` (`cxx:830-942`): drops `myV` entries that have
    /// no plausible analog in `thePoint` and records `myNear` / `aDiffs`.
    ///
    /// The two `a--` sites are the OCCT way of not advancing the outer index
    /// after a removal (`for` adds it back); in the `while` loop below they are
    /// the `advance` flag.
    pub fn filter_points(&mut self, other: &FilletPoint) {
        let mut diffs: Vec<f64> = Vec::new();
        let dx = other.param - self.param;
        let mut a = 1usize;
        while a <= self.v.len() {
            // `int aNear = 0; double aDiff = aDX * 10000.;`
            let mut a_near = 0usize;
            let mut a_diff = dx * 10000.0;
            let ay = self.v[a - 1] + self.d[a - 1] * dx;
            for b in 1..=other.v.len() {
                // `aY2 = aY + aDX * (thePoint->myD.Value(b) - myD.Value(a)) / 2.0;`
                let ay2 = ay + dx * (other.d[b - 1] - self.d[a - 1]) / 2.0;
                if a_near == 0 || (ay2 - other.v[b - 1]).abs() < a_diff.abs() {
                    a_near = b;
                    a_diff = ay2 - other.v[b - 1];
                }
            }
            // `if (aNear)` blocks: identical sign on the same side of the
            // interval, then the constant-second-derivative extrapolation.
            if a_near != 0 && self.v[a - 1] * other.v[a_near - 1] > 0.0 {
                if self.v[a - 1] * self.d[a - 1] > 0.0 {
                    if self.d[a - 1].abs() > CONFUSION {
                        a_near = 0;
                    }
                } else if self.v[a - 1].abs() > other.v[a_near - 1].abs() {
                    if other.v[a_near - 1] * other.d[a_near - 1] < 0.0
                        && other.d[a_near - 1].abs() > CONFUSION
                    {
                        a_near = 0;
                    }
                }
            }
            if a_near != 0
                && self.v[a - 1] * other.v[a_near - 1] > 0.0
                && (self.v[a - 1] + self.d[a - 1] * dx) * self.v[a - 1] > CONFUSION
                && (other.v[a_near - 1] + other.d[a_near - 1] * dx) * other.v[a_near - 1]
                    > CONFUSION
            {
                a_near = 0;
            }
            if a_near != 0 && (a_diff / dx).abs() > 1.0e+7 {
                a_near = 0;
            }

            let mut advance = true;
            if a_near == 0 {
                // "there is no near: remove it from the list" (`cxx:900-907`).
                self.v.remove(a - 1);
                self.d.remove(a - 1);
                self.valid.remove(a - 1);
                advance = false;
            } else {
                let mut returns_near = false;
                let mut removed_old = false;
                for b in 1..=self.near.len() {
                    if self.near[b - 1] != a_near {
                        continue;
                    }
                    if diffs[b - 1].abs() < a_diff.abs() {
                        // "return this 'near'": the current entry goes away.
                        returns_near = true;
                        self.v.remove(a - 1);
                        self.d.remove(a - 1);
                        self.valid.remove(a - 1);
                    } else {
                        // "remove the old 'near'": position `b`, not `a`.
                        removed_old = true;
                        self.v.remove(b - 1);
                        self.d.remove(b - 1);
                        self.valid.remove(b - 1);
                        self.near.remove(b - 1);
                        diffs.remove(b - 1);
                    }
                    break;
                }
                if !returns_near {
                    // `if (!aFound)`: reached both when nothing matched and
                    // when the old near was replaced.
                    self.near.push(a_near);
                    diffs.push(a_diff);
                }
                if returns_near || removed_old {
                    advance = false;
                }
            }
            if advance {
                a += 1;
            }
        }
    }

    /// `Copy()` (`cxx:944-955`): `myParam`, `myV`, `myD`, `myValid` only —
    /// **not** `myParam2`, `myCenter` or `myNear` (the header calls this out:
    /// "warning: this is not the full copy").
    pub fn copy(&self) -> FilletPoint {
        let mut a_copy = FilletPoint::new(self.param);
        for a in 1..=self.v.len() {
            a_copy.v.push(self.v[a - 1]);
            a_copy.d.push(self.d[a - 1]);
            a_copy.valid.push(self.valid[a - 1]);
        }
        a_copy
    }

    /// `hasSolution(theRadius)` (`cxx:957-969`): the 1-based index of the first
    /// value that is a root of `|v| + r^2`, or `None` for OCCT's `0`.
    pub fn has_solution(&self, radius: f64) -> Option<usize> {
        for a in 1..=self.v.len() {
            let root = (self.v[a - 1].abs() + radius * radius).abs().sqrt();
            if (root - radius).abs() < CONFUSION {
                return Some(a);
            }
        }
        None
    }

    /// `remove(theIndex)` (`cxx:971-977`).
    ///
    /// The only OCCT body that drops the four parallel sequences at once; it
    /// has no caller in `ChFi2d_FilletAlgo.cxx` (`FilterPoints` removes three
    /// or five at a time by hand), so the port keeps it for parity only.
    #[allow(dead_code)]
    pub fn remove(&mut self, index: usize) {
        self.v.remove(index - 1);
        self.d.remove(index - 1);
        self.valid.remove(index - 1);
        self.near.remove(index - 1);
    }
}
