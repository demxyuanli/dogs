//! The vertex sections of `BRepFill_Sweep::BuildShell`: the two boundary
//! columns (`cxx:2482-2662`) and the per-cell vertex construction with its
//! mergers (`cxx:2692-2798`), plus the `MergeVertex` (`cxx:3984-4006`) and
//! `UpdateVertex` (`cxx:4010-4045`) helpers they call.
//!
//! Translation notes:
//! - `myLoc->PerformVertex(Index, Input, TolMin, Output)` writes its result
//!   through an out-parameter; [`LocationLawBase::perform_vertex`] returns the
//!   vertex, which the caller stores back into the table cell.
//! - `MergeVertex(V1, V2)` assigns `V2 = V1` through a `TopoDS_Shape&`, so the
//!   port reports the decision and rewrites the cell
//!   ([`ShellBuild::merge_vertex_cells`]).
//! - `TopoDS::Vertex(shape)` on a null cell is undefined in the reference; the
//!   port treats a null or non-vertex cell as "not mergeable" / "nothing to
//!   update" and says so at each site.

use crate::brep_tool::BRepTool;
use crate::shape::{TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;

use super::{make_vertex_at, oriented_edge_vertices, reset_location_and_move, ShellBuild};
use crate::brep_fill_sweep::helpers;

impl ShellBuild<'_> {
    /// `MergeVertex(V1, V2)` (`cxx:3984-4006`): `V2 = V1` when the two points are
    /// within `max(Tol(V1), Tol(V2), myTol3d)`.
    fn merge_vertex_points(tol3d: f64, v1: &TopoShape, v2: &TopoShape) -> bool {
        if v1.shape_type() != crate::abs::ShapeType::Vertex
            || v2.shape_type() != crate::abs::ShapeType::Vertex
        {
            return false;
        }
        let v1 = Vertex(v1.clone());
        let v2 = Vertex(v2.clone());
        let mut tol = BRepTool::vertex_tolerance(&v1).max(BRepTool::vertex_tolerance(&v2));
        if tol < tol3d {
            tol = tol3d;
        }
        BRepTool::vertex_point(&v1).distance(&BRepTool::vertex_point(&v2)) <= tol
    }

    /// `MergeVertex(Vertex(r1, c1), Vertex(r2, c2))` over the vertex table:
    /// on `true` the second cell becomes the first one's shape.
    pub(super) fn merge_vertex_cells(
        &mut self,
        (r1, c1): (usize, usize),
        (r2, c2): (usize, usize),
    ) -> bool {
        let Some(v1) = self.tables.vertex.value(r1, c1).cloned() else {
            return false;
        };
        let Some(v2) = self.tables.vertex.value(r2, c2).cloned() else {
            return false;
        };
        if !Self::merge_vertex_points(self.sweep.my_tol_3d, &v1, &v2) {
            return false;
        }
        self.tables.vertex.set_value(r2, c2, Some(v1));
        true
    }

    /// `UpdateVertex(ipath, isec, ErrApp, Param, V)` (`cxx:4010-4045`): grow the
    /// stored tolerance of `V` so it also covers `myLoc->PerformVertex`'s image
    /// of the section vertex.
    pub(super) fn update_vertex_cell(
        &mut self,
        ipath: usize,
        isec: usize,
        v: &TopoShape,
        err_app: f64,
        param: f64,
    ) {
        let Some(input) = self.sweep.my_sec.vertex(isec, param) else {
            return;
        };
        let tol_min = err_app + self.sweep.my_sec.vertex_tol(isec - 1, param);
        let Ok(moved) = self
            .sweep
            .my_loc
            .base()
            .perform_vertex(ipath, &input, tol_min, 0)
        else {
            return;
        };
        if v.shape_type() != crate::abs::ShapeType::Vertex {
            return;
        }
        let stored = Vertex(v.clone());
        // `Tol = Tolerance(vv) + P1.Distance(P2)` (`cxx:4030-4033`).
        let tol = BRepTool::vertex_tolerance(&moved)
            + BRepTool::vertex_point(&moved).distance(&BRepTool::vertex_point(&stored));
        if tol > BRepTool::vertex_tolerance(&stored) {
            GeometryRegistry::global().raise_vertex_tolerance(v, tol);
        }
    }

    /// `myLoc->PerformVertex(ipath, Input, TolMin, Output)`
    /// (`cxx:2703-2711`, `cxx:2725-2735`, ...).
    fn perform_vertex(
        &self,
        ipath: usize,
        input: Option<Vertex>,
        tol_min: f64,
    ) -> Result<TopoShape, String> {
        let input = input.ok_or_else(|| {
            format!("BuildShell: null section vertex (BRepFill_Sweep.cxx:2703, ipath={ipath})")
        })?;
        Ok(self
            .sweep
            .my_loc
            .base()
            .perform_vertex(ipath, &input, tol_min, 0)?
            .0)
    }

    /// The vertex / minimum tolerance the reference reads from `SecVertex` /
    /// `VError` when the section is constant, or from the section law
    /// otherwise (`cxx:2701-2712`, `cxx:2724-2736`, `cxx:2786-2796`).
    fn section_vertex(
        &self,
        isec: usize,
        param: f64,
        err_app: f64,
    ) -> (Option<Vertex>, f64) {
        if self.const_section {
            let vertex = self.locals.sec_vertex[isec].clone().and_then(|s| {
                if s.shape_type() == crate::abs::ShapeType::Vertex {
                    Some(Vertex(s))
                } else {
                    None
                }
            });
            (vertex, err_app + self.locals.v_error[isec])
        } else {
            (
                self.sweep.my_sec.vertex(isec, param),
                err_app + self.sweep.my_sec.vertex_tol(isec - 1, param),
            )
        }
    }

    /// (2.0) `cxx:2482-2662`: the vertex columns of the first and the last
    /// sweep surface.
    pub(super) fn construct_boundary_vertices(&mut self) -> Result<(), String> {
        // First column (`cxx:2482-2550`).
        for isec in 1..=self.nb_law + 1 {
            let (row, col) = if isec > self.nb_law {
                (self.nb_law, 1)
            } else {
                (isec, 1)
            };
            let Some(s) = self.tables.tab_s.value(row, col).cloned() else {
                return Err(format!(
                    "BuildShell: null sweep surface ({row}, {col}) at BRepFill_Sweep.cxx:2499"
                ));
            };
            let ureverse = self.tables.u_reverse.value(row, col);
            let exuv = self.tables.exch_uv.value(row, col);
            let (mut u_first, mut u_last) = s.u_range();
            let (mut v_first, mut v_last) = s.v_range();
            if ureverse {
                if exuv {
                    std::mem::swap(&mut v_first, &mut v_last);
                } else {
                    std::mem::swap(&mut u_first, &mut u_last);
                }
            }
            // `cxx:2515-2534`.
            let (u, v) = if isec != self.nb_law + 1 {
                (u_first, v_first)
            } else if exuv {
                (u_first, v_last)
            } else {
                (u_last, v_first)
            };
            let p = s.value(u, v);
            let tol = self.sweep.my_sec.vertex_tol(isec - 1, self.locals.vi[1]);
            // `cxx:2536-2549`.
            if self.tables.vertex.is_null(isec, 1) {
                self.tables
                    .vertex
                    .set_value(isec, 1, Some(make_vertex_at(&p, tol)));
            } else if let Some(slot) = self.tables.vertex.value_mut(isec, 1) {
                reset_location_and_move(slot, &p, tol);
            }
        }

        // `LastShape` (`cxx:2553-2593`).
        if self.sweep.last_shape.is_some() && self.i_last == self.nb_path_law + 1 {
            let last_shape = self
                .sweep
                .last_shape
                .clone()
                .ok_or_else(|| "BuildShell: null LastShape".to_string())?;
            self.sweep.my_sec.base().init(&last_shape);
            for isec in 1..=self.nb_law {
                let edge = self.sweep.my_sec.base().current_edge().ok_or_else(|| {
                    format!("BuildShell: LastShape has no edge {isec} (BRepFill_Sweep.cxx:2558)")
                })?;
                if self.tables.v_edge.is_null(isec, self.nb_path + 1) {
                    self.tables
                        .v_edge
                        .set_value(isec, self.nb_path + 1, Some(edge.0.clone()));
                }
                if self.tables.vertex.is_null(isec + 1, self.nb_path + 1) {
                    let Some(ve) = self.tables.v_edge.value(isec, self.nb_path + 1).cloned() else {
                        continue;
                    };
                    let (vf, vl) = oriented_edge_vertices(&crate::shape::Edge(ve.clone()));
                    let boundary = if ve.orientation() == crate::abs::Orientation::Reversed {
                        vf
                    } else {
                        vl
                    };
                    if let Some(b) = boundary {
                        self.tables
                            .vertex
                            .set_value(isec + 1, self.nb_path + 1, Some(b.0.clone()));
                        self.update_vertex_cell(
                            self.i_last - 1,
                            isec + 1,
                            &b.0,
                            self.tables.tab_err.value(isec, self.nb_path),
                            self.locals.vi[self.nb_path + 1],
                        );
                    }
                    continue;
                }
                let current = self
                    .tables
                    .vertex
                    .value(isec + 1, self.nb_path + 1)
                    .cloned();
                if let Some(current) = current {
                    self.update_vertex_cell(
                        self.i_last - 1,
                        isec + 1,
                        &current,
                        self.tables.tab_err.value(isec, self.nb_path),
                        self.locals.vi[self.nb_path + 1],
                    );
                }
            }

            // `cxx:2581-2592`.
            if self.tables.vertex.is_null(1, self.nb_path + 1) {
                if let Some(ve) = self.tables.v_edge.value(1, self.nb_path + 1).cloned() {
                    let (vf, vl) = oriented_edge_vertices(&crate::shape::Edge(ve.clone()));
                    let boundary = if ve.orientation() == crate::abs::Orientation::Reversed {
                        vl
                    } else {
                        vf
                    };
                    if let Some(b) = boundary {
                        self.tables
                            .vertex
                            .set_value(1, self.nb_path + 1, Some(b.0.clone()));
                        self.update_vertex_cell(
                            self.i_last - 1,
                            1,
                            &b.0,
                            self.tables.tab_err.value(1, self.nb_path),
                            self.locals.vi[self.nb_path + 1],
                        );
                    }
                }
            } else if let Some(current) = self.tables.vertex.value(1, self.nb_path + 1).cloned() {
                self.update_vertex_cell(
                    self.i_last - 1,
                    1,
                    &current,
                    self.tables.tab_err.value(1, self.nb_path),
                    self.locals.vi[self.nb_path + 1],
                );
            }
        }

        // Last column (`cxx:2595-2662`).
        for isec in 1..=self.nb_law + 1 {
            let (row, col) = if isec > self.nb_law {
                (self.nb_law, self.nb_path)
            } else {
                (isec, self.nb_path)
            };
            let Some(s) = self.tables.tab_s.value(row, col).cloned() else {
                return Err(format!(
                    "BuildShell: null sweep surface ({row}, {col}) at BRepFill_Sweep.cxx:2610"
                ));
            };
            let ureverse = self.tables.u_reverse.value(row, col);
            let exuv = self.tables.exch_uv.value(row, col);
            let (mut u_first, mut u_last) = s.u_range();
            let (mut v_first, mut v_last) = s.v_range();
            if ureverse {
                if exuv {
                    std::mem::swap(&mut v_first, &mut v_last);
                } else {
                    std::mem::swap(&mut u_first, &mut u_last);
                }
            }
            // `cxx:2628-2645`.
            let (u, v) = if isec == self.nb_law + 1 {
                (u_last, v_last)
            } else if exuv {
                (u_last, v_first)
            } else {
                (u_first, v_last)
            };
            let p = s.value(u, v);
            let tol = self
                .sweep
                .my_sec
                .vertex_tol(isec - 1, self.locals.vi[self.nb_path + 1]);
            if self.tables.vertex.is_null(isec, self.nb_path + 1) {
                self.tables
                    .vertex
                    .set_value(isec, self.nb_path + 1, Some(make_vertex_at(&p, tol)));
            } else if let Some(slot) = self.tables.vertex.value_mut(isec, self.nb_path + 1) {
                reset_location_and_move(slot, &p, tol);
            }
        }
        Ok(())
    }

    /// (2.1) `cxx:2692-2798`: the new vertices of one `(isec, ipath)` cell.
    pub(super) fn construct_vertices(&mut self, ipath: usize, isec: usize) -> Result<(), String> {
        let ipath_global = self.i_first + ipath - 1;

        if isec == 1 {
            // All first (`cxx:2696-2712`).
            if ipath == 1 && self.tables.vertex.is_null(1, 1) {
                let (input, tol) = self.section_vertex(1, self.locals.vi[1], 0.0);
                let v = self.perform_vertex(ipath_global - 1, input, tol)?;
                self.tables.vertex.set_value(1, 1, Some(v));
            }
            // The first and the next column (`cxx:2714-2743`).
            if self.v_close && ipath == self.nb_path {
                let v = self.tables.vertex.value(1, 1).cloned();
                self.tables.vertex.set_value(1, ipath + 1, v);
            } else if self.tables.vertex.is_null(1, ipath + 1) {
                let param = self.locals.vi[ipath + 1];
                let (input, tol) =
                    self.section_vertex(1, param, self.tables.tab_err.value(1, ipath));
                let v = self.perform_vertex(ipath_global, input, tol)?;
                self.tables.vertex.set_value(1, ipath + 1, Some(v));
                if self.merge_vertex_cells((1, ipath), (1, ipath + 1)) {
                    // `UEdge(1, ipath) = NullEdge(Vertex(1, ipath))`.
                    if let Some(cell) = self.tables.vertex.value(1, ipath).cloned() {
                        let mut v = Vertex(cell);
                        let e = helpers::null_edge(&mut v);
                        self.tables.u_edge.set_value(1, ipath, Some(e.0));
                    }
                }
            }
        }

        if ipath == 1 {
            // `cxx:2745-2772`.
            if self.u_close && isec == self.nb_law {
                let v = self.tables.vertex.value(1, 1).cloned();
                self.tables.vertex.set_value(isec + 1, 1, v);
            } else if self.tables.vertex.is_null(isec + 1, 1) {
                let param = self.locals.vi[1];
                let (input, tol) =
                    self.section_vertex(isec + 1, param, self.tables.tab_err.value(isec, 1));
                let v = self.perform_vertex(ipath_global - 1, input, tol)?;
                self.tables.vertex.set_value(isec + 1, 1, Some(v));
                if self.merge_vertex_cells((isec, 1), (isec + 1, 1)) {
                    // `VEdge(isec, 1) = NullEdge(Vertex(isec, 1))`.
                    if let Some(cell) = self.tables.vertex.value(isec, 1).cloned() {
                        let mut v = Vertex(cell);
                        let e = helpers::null_edge(&mut v);
                        self.tables.v_edge.set_value(isec, 1, Some(e.0));
                    }
                }
            }
        }

        // `cxx:2774-2799`.
        if self.u_close && isec == self.nb_law {
            let v = self.tables.vertex.value(1, ipath + 1).cloned();
            self.tables.vertex.set_value(isec + 1, ipath + 1, v);
        } else if self.v_close && ipath == self.nb_path {
            let v = self.tables.vertex.value(isec + 1, 1).cloned();
            self.tables.vertex.set_value(isec + 1, ipath + 1, v);
        } else if self.tables.vertex.is_null(isec + 1, ipath + 1) {
            let param = self.locals.vi[ipath + 1];
            let (input, tol) =
                self.section_vertex(isec + 1, param, self.tables.tab_err.value(isec, ipath));
            let v = self.perform_vertex(ipath_global, input, tol)?;
            self.tables.vertex.set_value(isec + 1, ipath + 1, Some(v));
        }
        Ok(())
    }
}
