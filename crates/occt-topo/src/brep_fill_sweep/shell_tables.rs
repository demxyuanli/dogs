//! The tables `BRepFill_Sweep::BuildShell` fills (`BRepFill_Sweep.cxx:2242-2261`)
//! and the per-build scratch arrays beside them.
//!
//! Translation notes:
//! - `NCollection_Array2<bool>` / `<double>` / `<occ::handle<Geom_Surface>>` are
//!   [`BoolTable`], [`DoubleTable`] and [`SurfaceTable`], and the shape tables
//!   are [`ShapeTable`]; every one of them keeps the 1-based `Value(row, col)`
//!   bounds of the reference. A row or column outside the allocated range
//!   answers the default value instead of reading past the array.
//! - `NCollection_Array1<TopoDS_Shape>` / `<double>` / `<bool>` of
//!   `cxx:2257-2261` are plain `Vec`s with one unused slot at index 0 so the
//!   port can index them with the reference's 1-based `isec` / `ipath`.
//! - `myFaces` is allocated with `NbPath` columns by `Build` (`cxx:3246`) but
//!   `cxx:2411` / `cxx:3016` / `cxx:3045` index it with `IPath`, which reaches
//!   `NbPath + 1`; [`ShellTables::new`] allocates that extra column so those
//!   writes stay inside the table (the reference runs off the end of the array
//!   there). The reachable single-part case has `IPath == ipath`.

use std::sync::Arc;

use occt_geom::surface::Surface;

use crate::shape::TopoShape;

use super::helpers::ShapeTable;

/// `NCollection_Array2<bool>`.
#[derive(Debug, Clone)]
pub struct BoolTable {
    nb_rows: usize,
    nb_cols: usize,
    values: Vec<bool>,
}

impl BoolTable {
    /// `NCollection_Array2<bool>(1, NbRows, 1, NbCols)` plus `Init`.
    pub fn new(nb_rows: usize, nb_cols: usize, init: bool) -> Self {
        Self {
            nb_rows,
            nb_cols,
            values: vec![init; nb_rows * nb_cols],
        }
    }

    /// `Array2::Value(row, col)`; `false` outside the bounds.
    pub fn value(&self, row: usize, col: usize) -> bool {
        if row == 0 || col == 0 || row > self.nb_rows || col > self.nb_cols {
            return false;
        }
        self.values[(row - 1) * self.nb_cols + (col - 1)]
    }

    /// `Array2::SetValue(row, col, value)`.
    pub fn set_value(&mut self, row: usize, col: usize, value: bool) {
        if row == 0 || col == 0 || row > self.nb_rows || col > self.nb_cols {
            return;
        }
        self.values[(row - 1) * self.nb_cols + (col - 1)] = value;
    }
}

/// `NCollection_Array2<double>`.
#[derive(Debug, Clone)]
pub struct DoubleTable {
    nb_rows: usize,
    nb_cols: usize,
    values: Vec<f64>,
}

impl DoubleTable {
    /// `NCollection_Array2<double>(1, NbRows, 1, NbCols)`; the reference fills
    /// every cell before reading it, so no `Init` value is needed.
    pub fn new(nb_rows: usize, nb_cols: usize) -> Self {
        Self {
            nb_rows,
            nb_cols,
            values: vec![0.0; nb_rows * nb_cols],
        }
    }

    /// `Array2::Value(row, col)`; `0.0` outside the bounds.
    pub fn value(&self, row: usize, col: usize) -> f64 {
        if row == 0 || col == 0 || row > self.nb_rows || col > self.nb_cols {
            return 0.0;
        }
        self.values[(row - 1) * self.nb_cols + (col - 1)]
    }

    /// `Array2::SetValue(row, col, value)`.
    pub fn set_value(&mut self, row: usize, col: usize, value: f64) {
        if row == 0 || col == 0 || row > self.nb_rows || col > self.nb_cols {
            return;
        }
        self.values[(row - 1) * self.nb_cols + (col - 1)] = value;
    }
}

/// `NCollection_Array2<occ::handle<Geom_Surface>>`; a null handle is `None`.
#[derive(Clone, Default)]
pub struct SurfaceTable {
    nb_rows: usize,
    nb_cols: usize,
    values: Vec<Option<Arc<dyn Surface>>>,
}

impl SurfaceTable {
    /// `NCollection_Array2<occ::handle<Geom_Surface>>(1, NbRows, 1, NbCols)`.
    pub fn new(nb_rows: usize, nb_cols: usize) -> Self {
        Self {
            nb_rows,
            nb_cols,
            values: vec![None; nb_rows * nb_cols],
        }
    }

    /// `Array2::Value(row, col)`.
    pub fn value(&self, row: usize, col: usize) -> Option<&Arc<dyn Surface>> {
        if row == 0 || col == 0 || row > self.nb_rows || col > self.nb_cols {
            return None;
        }
        self.values[(row - 1) * self.nb_cols + (col - 1)].as_ref()
    }

    /// `Array2::ChangeValue(row, col) = surface`.
    pub fn set_value(&mut self, row: usize, col: usize, surface: Arc<dyn Surface>) {
        if row == 0 || col == 0 || row > self.nb_rows || col > self.nb_cols {
            return;
        }
        self.values[(row - 1) * self.nb_cols + (col - 1)] = Some(surface);
    }
}

/// The `NCollection_Array2` tables of `cxx:2242-2252`.
pub struct ShellTables {
    /// `ExchUV(1, NbLaw, 1, NbPath)`.
    pub exch_uv: BoolTable,
    /// `UReverse(1, NbLaw, 1, NbPath)`.
    pub u_reverse: BoolTable,
    /// `Degenerated(1, NbLaw, 1, NbPath)`, initialised to `false`.
    pub degenerated: BoolTable,
    /// `TabErr(1, NbLaw, 1, NbPath)`.
    pub tab_err: DoubleTable,
    /// `TabS(1, NbLaw, 1, NbPath)`.
    pub tab_s: SurfaceTable,
    /// `UEdge(1, NbLaw + 1, 1, NbPath)`.
    pub u_edge: ShapeTable,
    /// `VEdge(1, NbLaw, 1, NbPath + 1)`.
    pub v_edge: ShapeTable,
    /// `Vertex(1, NbLaw + 1, 1, NbPath + 1)`, every cell null.
    pub vertex: ShapeTable,
    /// `myFaces(1, NbLaw, 1, NbPath)` — the extra column is the `IPath`
    /// overflow noted in the module header.
    pub faces: ShapeTable,
}

impl ShellTables {
    /// The tables of `cxx:2242-2252` for `NbLaw` section laws and `NbPath`
    /// paths; `faces` mirrors `myFaces` (`cxx:3246`).
    pub fn new(nb_law: usize, nb_path: usize) -> Self {
        Self {
            exch_uv: BoolTable::new(nb_law, nb_path, false),
            u_reverse: BoolTable::new(nb_law, nb_path, false),
            degenerated: BoolTable::new(nb_law, nb_path, false),
            tab_err: DoubleTable::new(nb_law, nb_path),
            tab_s: SurfaceTable::new(nb_law, nb_path),
            u_edge: ShapeTable::new(nb_law + 1, nb_path),
            v_edge: ShapeTable::new(nb_law, nb_path + 1),
            vertex: ShapeTable::new(nb_law + 1, nb_path + 1),
            faces: ShapeTable::new(nb_law, nb_path + 1),
        }
    }
}

/// The `NCollection_Array1` locals of `cxx:2257-2261`, 1-based.
pub struct ShellLocals {
    /// `SecVertex(1, NbLaw + 1)`.
    pub sec_vertex: Vec<Option<TopoShape>>,
    /// `VError(1, NbLaw + 1)`.
    pub v_error: Vec<f64>,
    /// `Vi(1, NbPath + 1)`.
    pub vi: Vec<f64>,
    /// `IsBuilt(1, NbLaw)`.
    pub is_built: Vec<bool>,
    /// `StartEdges(1, NbLaw)`.
    pub start_edges: Vec<Option<TopoShape>>,
}

impl ShellLocals {
    /// Null `SecVertex` / `StartEdges` cells and the `false` `IsBuilt` flags of
    /// `cxx:2257-2261` / `cxx:2358-2359`.
    pub fn new(nb_law: usize, nb_path: usize) -> Self {
        Self {
            sec_vertex: vec![None; nb_law + 2],
            v_error: vec![0.0; nb_law + 2],
            vi: vec![0.0; nb_path + 2],
            is_built: vec![false; nb_law + 1],
            start_edges: vec![None; nb_law + 1],
        }
    }
}
