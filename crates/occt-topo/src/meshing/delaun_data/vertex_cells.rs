use std::collections::HashMap;

use occt_core::gp::{GpPnt2d, GpXY};
use occt_core::precision::CONFUSION;

use super::super::delaun_types::{DelaunVertex, VertexState};

/// Uniform-grid cell filter over vertices used to collapse coincident UV points.
///
/// Source: `BRepMesh_VertexTool` + `BRepMesh_VertexInspector` (the OCCT class
/// stores a `NCollection_CellFilter`); this is a `HashMap`-backed equivalent
/// with the same cell size / tolerance semantics.
pub(super) struct VertexCellFilter {
    /// `(U, V)` cell sizes. Source: `BRepMesh_VertexTool::SetCellSize`.
    pub(super) cell_size: (f64, f64),
    pub(super) tolerance: (f64, f64),
    pub(super) cells: HashMap<(i64, i64), Vec<i32>>,
}

impl VertexCellFilter {
    pub(super) fn new(cell_size: f64, tolerance: (f64, f64)) -> Self {
        Self {
            cell_size: (cell_size, cell_size),
            tolerance,
            cells: HashMap::new(),
        }
    }

    /// Cell index of one coordinate. Source: `NCollection_CellFilter::Cell::Cell`
    /// (`NCollection_CellFilter.hxx:252-267`). The quotient is truncated toward
    /// zero by the `Cell_IndexType` (an `int`) cast, *not* floored, with an
    /// `fmod` guard for values outside `INT_MIN + 1 .. INT_MAX - 1`.
    fn cell_coord(&self, coord: f64, cell_size: f64) -> i64 {
        let val = coord / cell_size;
        if val > (i32::MAX as f64 - 1.0) {
            (val % (i32::MAX as f64)) as i64
        } else if val < (i32::MIN as f64 + 1.0) {
            (val % (i32::MIN as f64)) as i64
        } else {
            val as i64
        }
    }

    fn cell_of(&self, p: GpXY) -> (i64, i64) {
        let cu = if self.cell_size.0 > 0.0 { self.cell_size.0 } else { 1.0 };
        let cv = if self.cell_size.1 > 0.0 { self.cell_size.1 } else { 1.0 };
        (self.cell_coord(p.x, cu), self.cell_coord(p.y, cv))
    }

    /// Cell index range covered by the tolerance box around `p`. Source:
    /// `BRepMesh_VertexTool::expandPoint` (`BRepMesh_VertexTool.hxx:136-142`):
    /// `[p - Tolerance, p + Tolerance]`.
    fn cell_range(&self, p: GpXY) -> ((i64, i64), (i64, i64)) {
        let (tx, ty) = self.tolerance;
        (
            self.cell_of(GpXY::new(p.x - tx, p.y - ty)),
            self.cell_of(GpXY::new(p.x + tx, p.y + ty)),
        )
    }

    /// Registers the vertex in *every* cell covered by its tolerance box.
    /// Source: `BRepMesh_VertexTool::Add` (`BRepMesh_VertexTool.cxx:70-81`) ->
    /// `NCollection_CellFilter::Add` range form (`CellFilter.hxx:174-183`). The
    /// target is prepended to each cell list, like `CellFilter::add`
    /// (`CellFilter.hxx:349-360`), so the newest vertex wins an exact-distance
    /// tie in `Inspect` (`BRepMesh_VertexInspector.cxx:44-51`).
    pub(super) fn add(&mut self, index: i32, p: GpPnt2d) {
        let (min, max) = self.cell_range(p.coord);
        for ci in min.0..=max.0 {
            for cj in min.1..=max.1 {
                self.cells.entry((ci, cj)).or_default().insert(0, index);
            }
        }
    }

    /// Removes the vertex from every cell covered by its tolerance box.
    /// Source: `BRepMesh_VertexTool::DeleteVertex` (`cxx:86-95`) /
    /// `Substitute` (`:100-110`) -> `NCollection_CellFilter::Remove` range form
    /// (`CellFilter.hxx:197-205`).
    pub(super) fn remove(&mut self, index: i32, p: GpPnt2d) {
        let (min, max) = self.cell_range(p.coord);
        for ci in min.0..=max.0 {
            for cj in min.1..=max.1 {
                if let Some(list) = self.cells.get_mut(&(ci, cj)) {
                    if let Some(pos) = list.iter().position(|&x| x == index) {
                        list.remove(pos);
                    }
                }
            }
        }
    }

    /// Returns the index of the closest live vertex within the coincidence
    /// tolerance of `p`, or `0` when none exists. Source:
    /// `BRepMesh_VertexTool::FindIndex` (`hxx:100-105`) inspects only the single
    /// cell of the query point, then `BRepMesh_VertexInspector::Inspect`
    /// (`cxx:24-53`) keeps the nearest; the multi-cell reach comes from `Add`.
    pub(super) fn find_index(&self, p: GpPnt2d, vertices: &[DelaunVertex]) -> i32 {
        let (tx, ty) = self.tolerance;
        let sq_tx = tx * tx;
        let sq_ty = ty * ty;
        let mut best = 0;
        let mut best_sq = f64::INFINITY;
        if let Some(list) = self.cells.get(&self.cell_of(p.coord)) {
            for &idx in list {
                let v = &vertices[(idx - 1) as usize];
                if v.state == VertexState::Deleted {
                    continue;
                }
                let dx = v.location.x() - p.x();
                let dy = v.location.y() - p.y();
                let sq = dx * dx + dy * dy;
                // `Inspect` compares the *squared* Y tolerance against
                // `Precision::Confusion()`: with the two-argument `SetTolerance`
                // (`Confusion, Confusion`) the squared value is 1e-14 < 1e-7, so
                // the Euclidean branch applies (`VertexInspector.cxx:32-42`).
                let in_tol = if sq_ty.abs() < CONFUSION {
                    sq < sq_tx
                } else {
                    dx * dx < sq_tx && dy * dy < sq_ty
                };
                if in_tol && sq < best_sq {
                    best_sq = sq;
                    best = idx;
                }
            }
        }
        best
    }
}
