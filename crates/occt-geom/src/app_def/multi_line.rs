//! `AppDef_MultiLine` — the sequence of `AppDef_MultiPointConstraint` columns
//! that `AppDef_Compute` / `AppDef_BSplineCompute` approximate in parallel.
//!
//! Source: `ModelingData/TKGeomBase/AppDef/AppDef_MultiLine.hxx:30-100` and
//! `AppDef_MultiLine.cxx:22-103` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `tabMult` (`hxx`) is a `NCollection_HArray1`, i.e. 1-based in OCCT and a
//!   possibly-null handle after the default constructor. The port stores a
//!   `Vec`; [`MultiLine::new`] leaves it empty (the null-handle state) and
//!   the two Dump-only accessors that OCCT would dereference anyway
//!   ([`MultiLine::nb_points`]) panic on it, exactly like `Value(1)` does.
//! - `AppDef_MultiLine(NbMult)` throws `Standard_ConstructionError` for a
//!   negative count (`cxx:27-31`); the parameter is `usize` here, so that
//!   branch cannot be reached.
//! - `Standard_OutOfRange` becomes a panic, as in
//!   `app_par_curves::multi_point`.

use occt_core::gp::{GpPnt, GpPnt2d};

use crate::app_def::multi_point_constraint::MultiPointConstraint;

/// `AppDef_MultiLine` (`AppDef_MultiLine.hxx:30-100`).
#[derive(Debug, Clone, Default)]
pub struct MultiLine {
    /// `tabMult` (`hxx`), 1-based in OCCT.
    tab_mult: Vec<MultiPointConstraint>,
}

impl MultiLine {
    /// `AppDef_MultiLine()` (`cxx:24`): null handle, no columns.
    pub fn new() -> Self {
        Self {
            tab_mult: Vec::new(),
        }
    }

    /// `AppDef_MultiLine(NbMult)` (`cxx:26-33`): `NbMult` default columns.
    pub fn with_size(nb_mult: usize) -> Self {
        Self {
            tab_mult: vec![MultiPointConstraint::new(); nb_mult],
        }
    }

    /// `AppDef_MultiLine(tabMultiP)` (`cxx:35-43`): copies the columns.
    pub fn from_multi_points(tab_multi_p: &[MultiPointConstraint]) -> Self {
        Self {
            tab_mult: tab_multi_p.to_vec(),
        }
    }

    /// `AppDef_MultiLine(tabP3d)` (`cxx:45-55`): one column per point, each
    /// column holding that single 3d point.
    pub fn from_3d(tab_p3d: &[GpPnt]) -> Self {
        let tab_mult = tab_p3d
            .iter()
            .map(|p| {
                let mut mp = MultiPointConstraint::with_size(1, 0);
                mp.set_point(1, p);
                mp
            })
            .collect();
        Self { tab_mult }
    }

    /// `AppDef_MultiLine(tabP2d)` (`cxx:57-67`): one column per point, each
    /// column holding that single 2d point.
    pub fn from_2d(tab_p2d: &[GpPnt2d]) -> Self {
        let tab_mult = tab_p2d
            .iter()
            .map(|p| {
                let mut mp = MultiPointConstraint::with_size(0, 1);
                mp.set_point2d(1, p);
                mp
            })
            .collect();
        Self { tab_mult }
    }

    /// `NbMultiPoints` (`cxx:69-72`).
    pub fn nb_multi_points(&self) -> usize {
        self.tab_mult.len()
    }

    /// `NbPoints` (`cxx:74-77`): the first column's total pole count.
    pub fn nb_points(&self) -> usize {
        match self.tab_mult.first() {
            Some(mp) => mp.nb_points() + mp.nb_points2d(),
            None => panic!("Standard_OutOfRange: AppDef_MultiLine::NbPoints"),
        }
    }

    /// `SetValue(Index, MPoint)` (`cxx:79-86`), `Index` 1-based.
    pub fn set_value(&mut self, index: i32, mpoint: &MultiPointConstraint) {
        if index <= 0 || index > self.tab_mult.len() as i32 {
            panic!("Standard_OutOfRange: AppDef_MultiLine::SetValue");
        }
        self.tab_mult[(index - 1) as usize] = mpoint.clone();
    }

    /// `Value(Index)` (`cxx:88-95`), `Index` 1-based.
    pub fn value(&self, index: i32) -> &MultiPointConstraint {
        if index <= 0 || index > self.tab_mult.len() as i32 {
            panic!("Standard_OutOfRange: AppDef_MultiLine::Value");
        }
        &self.tab_mult[(index - 1) as usize]
    }

    /// `Dump(o)` (`cxx:97-103`): the header plus the column count.
    pub fn dump(&self) -> String {
        let mut out = String::new();
        out.push_str("AppDef_MultiLine dump:\n");
        out.push_str(&format!(
            "It contains {} MultiPointConstraint\n",
            self.tab_mult.len()
        ));
        out
    }
}
