//! `ChFiDS_CircSection` — a section of fillet.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiDS/ChFiDS_CircSection.hxx:27-48`,
//! `ChFiDS_CircSection.cxx:22-63` (OCCT 8.0.0).
//!
//! `Simul()` of [`super::surf_data::SurfData`] hands out an
//! `occ::handle<NCollection_HArray1<ChFiDS_CircSection>>`
//! (`ChFi3d_FilBuilder.cxx:488-489`, `:519`), so the array shape is ported here
//! as [`CircSectionArray`] and kept 1-based, matching every `Value(i)` /
//! `ChangeValue(i)` in OCCT's `ChFi3d`.

use occt_core::gp::{GpCirc, GpLin};

/// Source: `ChFiDS_CircSection.hxx:29-48`.
///
/// OCCT keeps a `gp_Circ` and a `gp_Lin` side by side and lets the caller pick
/// which one to read; the two are independent members, not a variant. The port
/// keeps both members so `Get(gp_Circ, ...)` after `Set(gp_Lin, ...)` still
/// returns the (unset) default circle, exactly like OCCT.
#[derive(Debug, Clone, Default)]
pub struct CircSection {
    /// `myCirc`
    my_circ: GpCirc,
    /// `myLin`
    my_lin: GpLin,
    /// `myF`
    my_f: f64,
    /// `myL`
    my_l: f64,
}

impl CircSection {
    /// `ChFiDS_CircSection()` (`cxx:22-27`): only `myF` and `myL` are
    /// initialized, both to `0.0`.
    pub fn new() -> Self {
        Self::default()
    }

    /// `Set(const gp_Circ& C, F, L)` (`cxx:32-38`).
    pub fn set_circ(&mut self, c: &GpCirc, f: f64, l: f64) {
        self.my_circ = c.clone();
        self.my_f = f;
        self.my_l = l;
    }

    /// `Set(const gp_Lin& C, F, L)` (`cxx:42-47`).
    pub fn set_lin(&mut self, c: &GpLin, f: f64, l: f64) {
        self.my_lin = *c;
        self.my_f = f;
        self.my_l = l;
    }

    /// `Get(gp_Circ& C, double& F, double& L)` (`cxx:51-56`).
    pub fn get_circ(&self) -> (GpCirc, f64, f64) {
        (self.my_circ.clone(), self.my_f, self.my_l)
    }

    /// `Get(gp_Lin& C, double& F, double& L)` (`cxx:58-63`).
    pub fn get_lin(&self) -> (GpLin, f64, f64) {
        (self.my_lin, self.my_f, self.my_l)
    }

    /// `myF` — the first parameter of the section.
    pub fn first(&self) -> f64 {
        self.my_f
    }

    /// `myL` — the last parameter of the section.
    pub fn last(&self) -> f64 {
        self.my_l
    }
}

/// `NCollection_HArray1<ChFiDS_CircSection>` (`ChFi3d_FilBuilder.cxx:488`,
/// `:519`, `:522`) — the payload of [`super::surf_data::SurfData`]'s `Simul`.
///
/// OCCT addresses it 1-based (`sec->ChangeValue(1)`, `Value(1)`), so the port
/// keeps the OCCT index domain rather than shifting by one: `new(1, n)` gives
/// `n` slots, `value(1)` is the first, `value(0)` and `value(n + 1)` give
/// `None` instead of OCCT's `Standard_RangeError`.
#[derive(Debug, Clone, Default)]
pub struct CircSectionArray {
    /// `Lower()` / `Upper()` of the underlying `NCollection_Array1`.
    lower: i32,
    upper: i32,
    /// The `NCollection_Array1` items, `upper - lower + 1` of them.
    items: Vec<CircSection>,
}

impl CircSectionArray {
    /// `new NCollection_HArray1<ChFiDS_CircSection>(1, n)` (`cxx:519`, `:541`).
    pub fn new(first: i32, last: i32) -> Self {
        let count = if last >= first {
            (last - first + 1) as usize
        } else {
            0
        };
        Self {
            lower: first,
            upper: last,
            items: vec![CircSection::new(); count],
        }
    }

    /// `NCollection_Array1::Lower()`.
    pub fn lower(&self) -> i32 {
        self.lower
    }

    /// `NCollection_Array1::Upper()`.
    pub fn upper(&self) -> i32 {
        self.upper
    }

    /// `NCollection_Array1::Length()`.
    pub fn len(&self) -> usize {
        self.items.len()
    }

    /// `NCollection_Array1::IsEmpty()`.
    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    /// `HArray1::Value(I)` (`cxx:522` reads `ChangeValue`).
    pub fn value(&self, i: i32) -> Option<&CircSection> {
        let k = i.checked_sub(self.lower)?;
        if k < 0 {
            return None;
        }
        self.items.get(k as usize)
    }

    /// `HArray1::ChangeValue(I)` (`cxx:522`, `:544`).
    pub fn change_value(&mut self, i: i32) -> Option<&mut CircSection> {
        let k = i.checked_sub(self.lower)?;
        if k < 0 {
            return None;
        }
        self.items.get_mut(k as usize)
    }

    /// `HArray1::SetValue(I, V)`.
    pub fn set_value(&mut self, i: i32, v: CircSection) {
        if let Some(slot) = self.change_value(i) {
            *slot = v;
        }
    }
}
