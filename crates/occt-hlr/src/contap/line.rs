//! `Contap_Line` - a contour line produced by `Contap_Contour`: a straight
//! line, a circle, a piece of the surface's restriction, or a walked line.
//!
//! Sources: `Contap_Line.hxx:37-99`, `.cxx:24-121`, `.lxx:29-86`.
//!
//! Deviation (handle -> owned value): OCCT keeps the walked line behind
//! `occ::handle<IntSurf_LineOn2S>` (`hxx:92`) so a caller that still holds the
//! handle sees the same buffer. Every producer hands the buffer over
//! (`Contap_Contour.cxx:317-332`, `:1712`), so [`ContapLine`] owns it and
//! [`ContapLine::line_on_2s_mut`] is the mutable view `Contap_Contour.cxx:1479`
//! needs.

use std::sync::Arc;

use occt_core::gp::{GpAx1, GpAx2, GpCirc, GpDir, GpLin, GpPnt};

use occt_topo::brep_adaptor_curve2d::Curve2dAdaptor;
use occt_topo::int_tools_wline::{PntOn2S, WLine};
use occt_topo::intsurf_transition::IntSurfTypeTrans;

use super::point::ContapPoint;
use super::types::ContapIType;

/// `Contap_Line` (`hxx:37-99`).
///
/// No `Debug`: the walked buffer [`WLine`] carries `IntPatch_Line` state and
/// the restriction arc is a `Curve2dAdaptor`, neither of which implements it.
#[derive(Clone)]
pub struct ContapLine {
    /// `Trans` (`hxx:90`), `IntSurf_Undecided` in the constructor.
    trans: IntSurfTypeTrans,
    /// `curv` (`hxx:92`).
    curv: Option<WLine>,
    /// `svtx` (`hxx:93`): the vertices sorted by `ParameterOnLine`.
    svtx: Vec<ContapPoint>,
    /// `thearc` (`hxx:94`).
    thearc: Option<Arc<Curve2dAdaptor>>,
    /// `typL` (`hxx:95`). OCCT's constructor leaves it indeterminate; the
    /// value here is the one `Clear()` writes (`cxx:70-76`).
    typ_l: ContapIType,
    /// `pt` (`hxx:96`).
    pt: GpPnt,
    /// `dir1` (`hxx:97`).
    dir1: GpDir,
    /// `dir2` (`hxx:98`).
    dir2: GpDir,
    /// `rad` (`hxx:99`).
    rad: f64,
}

impl Default for ContapLine {
    /// `Contap_Line()` (`cxx:24-28`): a fresh vertex sequence and an UNDECIDED
    /// transition.
    fn default() -> Self {
        Self::new()
    }
}

impl ContapLine {
    /// `Contap_Line()` (`cxx:24-28`).
    pub fn new() -> Self {
        Self {
            trans: IntSurfTypeTrans::Undecided,
            curv: None,
            svtx: Vec::new(),
            thearc: None,
            typ_l: ContapIType::Walking,
            pt: GpPnt::default(),
            dir1: GpDir::default(),
            dir2: GpDir::default(),
            rad: 0.0,
        }
    }

    /// `SetLineOn2S(L)` (`cxx:96-100`).
    pub fn set_line_on_2s(&mut self, l: WLine) {
        self.curv = Some(l);
        self.typ_l = ContapIType::Walking;
    }

    /// `Clear()` (`cxx:63-76`): empties the walked buffer when one is set,
    /// resets the vertices and sets the type back to `Contap_Walking`.
    pub fn clear(&mut self) {
        if let Some(c) = self.curv.as_mut() {
            c.clear();
        }
        self.svtx.clear();
        self.typ_l = ContapIType::Walking;
    }

    /// `LineOn2S()` (`lxx:28-31`).
    pub fn line_on_2s(&self) -> Option<&WLine> {
        self.curv.as_ref()
    }

    /// The mutable `LineOn2S()`: `Contap_Contour` mutates through the handle
    /// (`Contap_Contour.cxx:1479`).
    pub fn line_on_2s_mut(&mut self) -> Option<&mut WLine> {
        self.curv.as_mut()
    }

    /// `ResetSeqOfVertex()` (`cxx:30-33`).
    pub fn reset_seq_of_vertex(&mut self) {
        self.svtx.clear();
    }

    /// `Add(POn2S)` (`lxx:33-36`): appends to the walked buffer.
    ///
    /// OCCT dereferences `curv` unconditionally; the port ignores the point
    /// when no walked line has been set.
    pub fn add_point_on_2s(&mut self, p: PntOn2S) {
        if let Some(c) = self.curv.as_mut() {
            c.add(p);
        }
    }

    /// `SetValue(L)` (`cxx:78-80`).
    pub fn set_value_line(&mut self, l: &GpLin) {
        self.pt = l.location();
        self.dir1 = l.direction();
        self.typ_l = ContapIType::Lin;
    }

    /// `SetValue(C)` (`cxx:82-88`).
    pub fn set_value_circle(&mut self, c: &GpCirc) {
        self.pt = c.location();
        self.dir1 = c.position().direction();
        self.dir2 = *c.position().x_direction();
        self.rad = c.radius();
        self.typ_l = ContapIType::Circle;
    }

    /// `SetValue(A)` (`cxx:90-94`).
    pub fn set_value_arc(&mut self, a: Arc<Curve2dAdaptor>) {
        self.thearc = Some(a);
        self.typ_l = ContapIType::Restriction;
    }

    /// `Add(Contap_Point&)` (`cxx:35-61`): insertion sort by
    /// `ParameterOnLine`.
    pub fn add_vertex(&mut self, p: ContapPoint) {
        let n = self.svtx.len();
        if n == 0 {
            self.svtx.push(p);
            return;
        }
        let prm = p.parameter_on_line();
        if prm > self.svtx[n - 1].parameter_on_line() {
            self.svtx.push(p);
        } else {
            // `for (int i = n - 1; i > 0; i--)` over the 1-based sequence.
            for i in (1..n).rev() {
                if prm > self.svtx[i - 1].parameter_on_line() {
                    self.svtx.insert(i, p);
                    return;
                }
            }
            self.svtx.insert(0, p);
        }
    }

    /// `NbVertex()` (`lxx:38-41`).
    pub fn nb_vertex(&self) -> usize {
        self.svtx.len()
    }

    /// `Vertex(Index)` (`lxx:43-46`), 1-based.
    pub fn vertex(&self, index: usize) -> Option<&ContapPoint> {
        index.checked_sub(1).and_then(|i| self.svtx.get(i))
    }

    /// The mutating `Vertex(Index)` (`lxx:43-46`, `ChangeSequence()(Index)` on
    /// a `const` method): `Contap_Contour.cxx:1049` and `:1862` set flags on
    /// vertices of a line reached through a const reference.
    pub fn vertex_mut(&mut self, index: usize) -> Option<&mut ContapPoint> {
        index.checked_sub(1).and_then(|i| self.svtx.get_mut(i))
    }

    /// `TypeContour()` (`lxx:48-51`).
    pub fn type_contour(&self) -> ContapIType {
        self.typ_l
    }

    /// `NbPnts()` (`lxx:53-60`): `Standard_DomainError` for a line that is not
    /// a walked line.
    pub fn nb_pnts(&self) -> Result<i32, &'static str> {
        if self.typ_l != ContapIType::Walking {
            return Err("Contap_Line::NbPnts: Standard_DomainError");
        }
        Ok(self.curv.as_ref().map(|c| c.nb_pnts()).unwrap_or(0))
    }

    /// `Point(Index)` (`lxx:62-68`), 1-based.
    pub fn point(&self, index: usize) -> Result<&PntOn2S, &'static str> {
        if self.typ_l != ContapIType::Walking {
            return Err("Contap_Line::Point: Standard_DomainError");
        }
        let c = self
            .curv
            .as_ref()
            .ok_or("Contap_Line::Point: no IntSurf_LineOn2S")?;
        let i = index
            .checked_sub(1)
            .filter(|i| *i < c.points.len())
            .ok_or("Contap_Line::Point: Standard_OutOfRange")?;
        Ok(&c.points[i])
    }

    /// `Line()` (`lxx:70-77`): `Standard_DomainError` for a non-line.
    pub fn line(&self) -> Result<GpLin, &'static str> {
        if self.typ_l != ContapIType::Lin {
            return Err("Contap_Line::Line: Standard_DomainError");
        }
        Ok(GpLin::new(GpAx1::new(self.pt, self.dir1)))
    }

    /// `Circle()` (`lxx:79-86`): `Standard_DomainError` for a non-circle.
    pub fn circle(&self) -> Result<GpCirc, &'static str> {
        if self.typ_l != ContapIType::Circle {
            return Err("Contap_Line::Circle: Standard_DomainError");
        }
        let pos = GpAx2::new(self.pt, self.dir1, self.dir2)
            .map_err(|_| "Contap_Line::Circle: degenerate gp_Ax2")?;
        Ok(GpCirc::new(pos, self.rad))
    }

    /// `Arc()` (`cxx:112-121`): `Standard_DomainError` for a non-restriction.
    pub fn arc(&self) -> Result<&Arc<Curve2dAdaptor>, &'static str> {
        if self.typ_l != ContapIType::Restriction {
            return Err("Contap_Line::Arc: Standard_DomainError");
        }
        self.thearc
            .as_ref()
            .ok_or("Contap_Line::Arc: null Adaptor2d_Curve2d")
    }

    /// `SetTransitionOnS(T)` (`cxx:102-105`).
    pub fn set_transition_on_s(&mut self, t: IntSurfTypeTrans) {
        self.trans = t;
    }

    /// `TransitionOnS()` (`cxx:107-110`).
    pub fn transition_on_s(&self) -> IntSurfTypeTrans {
        self.trans
    }
}
