//! `BRepBlend_Line` — the result of a blending-line computation: the sequence
//! of [`BlendPoint`]s walked between two extremities, plus the transitions and
//! start/end extremities of the line.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_Line.hxx:37-186`,
//! `BRepBlend_Line.lxx:16-130`, `BRepBlend_Line.cxx:25-58` (OCCT 8.0.0).
//!
//! Divergences that are not heuristics:
//! - `Standard_DomainError` from `TransitionOnS1` / `TransitionOnS2` /
//!   `TransitionOnS` (`lxx:57-64`, `:67-73`, `:123-130`) becomes `Err`,
//! - `NCollection_Sequence<Blend_Point>` becomes `Vec<BlendPoint>`; OCCT
//!   sequences are `1`-based, so [`Line::point`] / [`Line::insert_before`] /
//!   [`Line::remove`] keep the OCCT indices and translate inside. A
//!   `Standard_OutOfRange` from the underlying sequence is `None` / a no-op
//!   here, which is the same observable behaviour for the in-range calls the
//!   walking kernels make,
//! - `BRepBlend_Extremity` members are `Option<...>` only in the sense that the
//!   OCCT default constructor leaves them at their own defaults; the port
//!   stores them by value and hands out `&Extremity`, as OCCT returns
//!   `const BRepBlend_Extremity&`.

use crate::blend::BlendPoint;
use crate::intsurf_transition::IntSurfTypeTrans;

use super::extremity::Extremity;

/// `class BRepBlend_Line : public Standard_Transient` (`hxx:37-186`).
#[derive(Debug, Clone)]
pub struct Line {
    /// `seqpt` (`hxx:178`)
    seqpt: Vec<BlendPoint>,
    /// `tras1` (`hxx:179`)
    tras1: IntSurfTypeTrans,
    /// `tras2` (`hxx:180`)
    tras2: IntSurfTypeTrans,
    /// `stp1` (`hxx:181`)
    stp1: Extremity,
    /// `stp2` (`hxx:182`)
    stp2: Extremity,
    /// `endp1` (`hxx:183`)
    endp1: Extremity,
    /// `endp2` (`hxx:184`)
    endp2: Extremity,
    /// `hass1` (`hxx:185`)
    hass1: bool,
    /// `hass2` (`hxx:186`)
    hass2: bool,
}

impl Default for Line {
    fn default() -> Self {
        Self::new()
    }
}

impl Line {
    /// `BRepBlend_Line()` (`cxx:25-32`).
    pub fn new() -> Self {
        Self {
            seqpt: Vec::new(),
            tras1: IntSurfTypeTrans::Undecided,
            tras2: IntSurfTypeTrans::Undecided,
            stp1: Extremity::new(),
            stp2: Extremity::new(),
            endp1: Extremity::new(),
            endp2: Extremity::new(),
            hass1: false,
            hass2: false,
        }
    }

    /// `Clear()` (`cxx:34-41`): empties the sequence and resets both
    /// transitions and both flags. The four extremities are left alone.
    pub fn clear(&mut self) {
        self.seqpt.clear();
        self.hass1 = false;
        self.hass2 = false;
        self.tras1 = IntSurfTypeTrans::Undecided;
        self.tras2 = IntSurfTypeTrans::Undecided;
    }

    /// `Append(P)` (`lxx:16-19`).
    pub fn append(&mut self, p: BlendPoint) {
        self.seqpt.push(p);
    }

    /// `Prepend(P)` (`lxx:21-24`).
    pub fn prepend(&mut self, p: BlendPoint) {
        self.seqpt.insert(0, p);
    }

    /// `InsertBefore(Index, P)` (`lxx:26-29`), `Index` `1`-based. An index
    /// below `1` or past `len + 1` is ignored (OCCT raises).
    pub fn insert_before(&mut self, index: usize, p: BlendPoint) {
        let Some(i) = index.checked_sub(1) else {
            return;
        };
        if i > self.seqpt.len() {
            return;
        }
        self.seqpt.insert(i, p);
    }

    /// `Remove(FromIndex, ToIndex)` (`lxx:31-34`), both `1`-based and
    /// inclusive. Out-of-range bounds are ignored; an inverted range removes
    /// nothing.
    pub fn remove(&mut self, from_index: usize, to_index: usize) {
        if from_index == 0 || to_index < from_index || from_index > self.seqpt.len() {
            return;
        }
        let from = from_index - 1;
        let to = to_index.min(self.seqpt.len());
        self.seqpt.drain(from..to);
    }

    /// `Set(TranS1, TranS2)` (`cxx:43-49`).
    pub fn set_transitions(&mut self, tran_s1: IntSurfTypeTrans, tran_s2: IntSurfTypeTrans) {
        self.hass1 = true;
        self.hass2 = true;
        self.tras1 = tran_s1;
        self.tras2 = tran_s2;
    }

    /// `Set(Trans)` (`cxx:51-57`): only `tras1` is written and `hass2` is
    /// cleared.
    pub fn set_transition(&mut self, trans: IntSurfTypeTrans) {
        self.hass1 = true;
        self.hass2 = false;
        self.tras1 = trans;
    }

    /// `SetStartPoints(StartPt1, StartPt2)` (`lxx:36-42`).
    pub fn set_start_points(&mut self, start_pt1: Extremity, start_pt2: Extremity) {
        self.stp1 = start_pt1;
        self.stp2 = start_pt2;
    }

    /// `SetEndPoints(EndPt1, EndPt2)` (`lxx:44-50`).
    pub fn set_end_points(&mut self, end_pt1: Extremity, end_pt2: Extremity) {
        self.endp1 = end_pt1;
        self.endp2 = end_pt2;
    }

    /// `NbPoints()` (`lxx:52-55`).
    pub fn nb_points(&self) -> usize {
        self.seqpt.len()
    }

    /// `Point(Index)` (`lxx:57-60`), `Index` `1`-based; `None` when out of
    /// range (OCCT raises `Standard_OutOfRange`).
    pub fn point(&self, index: usize) -> Option<&BlendPoint> {
        self.seqpt.get(index.checked_sub(1)?)
    }

    /// `TransitionOnS1()` (`lxx:62-69`).
    pub fn transition_on_s1(&self) -> Result<IntSurfTypeTrans, String> {
        if !self.hass1 {
            return Err("BRepBlend_Line::TransitionOnS1: Standard_DomainError".to_string());
        }
        Ok(self.tras1)
    }

    /// `TransitionOnS2()` (`lxx:71-78`).
    pub fn transition_on_s2(&self) -> Result<IntSurfTypeTrans, String> {
        if !self.hass2 {
            return Err("BRepBlend_Line::TransitionOnS2: Standard_DomainError".to_string());
        }
        Ok(self.tras2)
    }

    /// `StartPointOnFirst()` (`lxx:80-83`).
    pub fn start_point_on_first(&self) -> &Extremity {
        &self.stp1
    }

    /// `StartPointOnSecond()` (`lxx:85-88`).
    pub fn start_point_on_second(&self) -> &Extremity {
        &self.stp2
    }

    /// `EndPointOnFirst()` (`lxx:90-93`).
    pub fn end_point_on_first(&self) -> &Extremity {
        &self.endp1
    }

    /// `EndPointOnSecond()` (`lxx:95-98`).
    pub fn end_point_on_second(&self) -> &Extremity {
        &self.endp2
    }

    /// `TransitionOnS()` (`lxx:123-130`): reads the `S1` transition, guarding
    /// `hass1` only.
    pub fn transition_on_s(&self) -> Result<IntSurfTypeTrans, String> {
        if !self.hass1 {
            return Err("BRepBlend_Line::TransitionOnS: Standard_DomainError".to_string());
        }
        Ok(self.tras1)
    }
}
