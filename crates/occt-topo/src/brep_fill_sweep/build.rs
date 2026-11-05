//! `BRepFill_Sweep::Build` (`BRepFill_Sweep.cxx:3213-3530`): the entry point
//! that dispatches to `BuildWire` for a vertex section and to `BuildShell` per
//! path window for a shell, then assembles the shell.
//!
//! Translation notes:
//! - `Tapes` / `Rails` are `HashMap<ShapeMapKey, ShapeTable>` (rows 1-6 and
//!   1-2) and `ReversedEdges` is `HashSet<ShapeMapKey>`, the same translation
//!   `myVEdgesModified` uses.
//! - `NCollection_HArray2<TopoDS_Shape> Bounds(NbLaw, 2)` is
//!   [`ShapeTable::new`]`(nb_law, 2)`; `Translate` is
//!   [`super::helpers::translate`], which keeps the reference's
//!   `ColLength()`-bounded row loop.
//! - `myTapes` is an `NCollection_HArray1<TopoDS_Shape>(1, NbLaw)` of empty
//!   shells (`cxx:3248-3254`). `myTapes->Value(Index)` is the 1-based
//!   [`BrepFillSweep::tape`].
//! - `BRepTools_Substitution` is
//!   [`crate::brep_tools_substitution::BRepToolsSubstitution`]; the reference's
//!   `Substitute(aSubstitute, OldEdge, NewEdge)` loop over the
//!   `myVEdgesModified` data map becomes the keyed form
//!   [`BRepToolsSubstitution::substitute_at_key`], because the map iteration
//!   only yields the hasher key.
//! - `TopExp::MapShapesAndAncestors(shell, EDGE, FACE, EFmap)` is
//!   [`crate::topo_tools::map_shapes_and_ancestors`]; the port walks
//!   `edges_of(shell)` and looks each edge up, which yields the same
//!   "more than one face ancestor" answer for every non-degenerated edge.
//! - A failed `PerformCorner` makes the reference `return` with a null
//!   `myShape`; the port stops the same way, leaving `isDone == false`.

use std::collections::{HashMap, HashSet};

use occt_core::kernel::geomabs::Shape;
use occt_geom::geom_fill::sweep::ApproxStyle;

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::brep_tools_substitution::BRepToolsSubstitution;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, Shell};
use crate::topo_tools::{edges_of, map_shapes_and_ancestors, shape_map_key, ShapeMapKey};

use super::helpers::{self, ShapeTable};
use super::same_parameter::correct_same_parameter;
use super::{BrepFillSweep, TransitionStyle};

impl BrepFillSweep {
    /// `BRepFill_Sweep::Build(ReversedEdges, Tapes, Rails, Transition,
    /// Continuity, Approx, Degmax, Segmax)` (`cxx:3213-3530`).
    #[allow(clippy::too_many_arguments)]
    pub fn build(
        &mut self,
        reversed_edges: &mut HashSet<ShapeMapKey>,
        tapes: &mut HashMap<ShapeMapKey, ShapeTable>,
        rails: &mut HashMap<ShapeMapKey, ShapeTable>,
        transition: TransitionStyle,
        continuity: Shape,
        approx: ApproxStyle,
        degmax: i32,
        segmax: i32,
    ) -> Result<(), String> {
        // `cxx:3226-3231`.
        self.my_continuity = continuity;
        self.my_approx_style = approx;
        self.my_degmax = degmax;
        self.my_segmax = segmax;
        let _ = self.correct_approx_parameters()?;

        // Wire (`cxx:3233-3237`).
        if self.my_sec.is_vertex() {
            self.is_done = self.build_wire(transition);
            return Ok(());
        }

        // Shell (`cxx:3238-3246`).
        let nb_path = self.my_loc.base().nb_law();
        let nb_law = self.my_sec.base().nb_law();
        let nb_trous = self.my_loc.base().nb_holes(self.my_tol_3d);
        let b = TopoBuilder::new();
        self.my_u_edges = ShapeTable::new(nb_law + 1, nb_path);
        self.my_v_edges = ShapeTable::new(nb_law, nb_path + 1);
        // `myFaces` is allocated with `NbPath` columns in the reference and
        // indexed with `IPath` (which reaches `NbPath + 1`); see the note in
        // `shell_tables`.
        self.my_faces = ShapeTable::new(nb_law, nb_path + 1);
        self.my_tapes = (0..nb_law)
            .map(|_| Some(b.make_shell(&[]).0))
            .collect();

        // `cxx:3255-3268`: the path windows.
        let mut bounds = ShapeTable::new(nb_law, 2);
        let mut trous: Vec<i32> = Vec::new();
        let mut nb_part = 1usize;
        if nb_trous > 0 {
            trous = vec![0; nb_trous];
            self.my_loc.base().holes(&mut trous);
            nb_part += nb_trous;
            if trous[nb_trous - 1] as usize == nb_path + 1 {
                nb_part -= 1;
            }
        }

        if nb_part == 1 {
            // `cxx:3271-3278`: this is done at once.
            let mut extend = 0.0;
            if nb_trous == 1 {
                extend = self.eval_extrapol(1, transition)?;
            }
            self.is_done = self.build_shell(
                transition,
                1,
                nb_path + 1,
                reversed_edges,
                tapes,
                rails,
                extend,
                extend,
            )?;
        } else {
            // `cxx:3280-3312`: this is done piece by piece.
            let mut i_first = 1usize;
            let mut is_done = true;
            let mut ii = 1usize;
            while ii <= nb_part && is_done {
                let i_last = if ii > nb_trous {
                    nb_path + 1
                } else {
                    trous[ii - 1] as usize
                };
                is_done = self.build_shell(
                    transition,
                    i_first,
                    i_last,
                    reversed_edges,
                    tapes,
                    rails,
                    self.eval_extrapol(i_first, transition)?,
                    self.eval_extrapol(i_last, transition)?,
                )?;
                if i_first > 1 {
                    helpers::translate(&self.my_v_edges, i_first, &mut bounds, 2);
                    if !self.perform_corner(i_first, transition, &bounds)? {
                        is_done = false;
                        self.is_done = is_done;
                        return Ok(());
                    }
                }
                i_first = i_last;
                helpers::translate(&self.my_v_edges, i_first, &mut bounds, 1);
                ii += 1;
            }
            self.is_done = is_done;
            if !is_done {
                return Ok(());
            }
        }

        // `cxx:3313-3324`: management of looping ends.
        if nb_trous > 0
            && self.my_loc.base().is_closed()
            && trous[nb_trous - 1] as usize == nb_path + 1
        {
            helpers::translate(&self.my_v_edges, nb_path + 1, &mut bounds, 1);
            helpers::translate(&self.my_v_edges, 1, &mut bounds, 2);
            if !self.perform_corner(1, transition, &bounds)? {
                self.is_done = false;
                return Ok(());
            }
            // `Translate(myVEdges, 1, myVEdges, NbPath + 1)` in place.
            for row in 1..=nb_law {
                let value = self.my_v_edges.value(row, 1).cloned();
                self.my_v_edges.set_value(row, nb_path + 1, value);
            }
        }

        // `cxx:3326-3345`: construction of the shell.
        let mut shell = Shell::new();
        let mut nb_faces = 0usize;
        for ipath in 1..=nb_path {
            for isec in 1..=nb_law {
                if let Some(face) = self.my_faces.value(isec, ipath) {
                    if !face.is_null() && face.shape_type() == ShapeType::Face {
                        b.add(&mut shell.0, face);
                        nb_faces += 1;
                    }
                }
            }
        }

        // `cxx:3344`: `(NbTrous > 0) ? (aNbFaces < NbLaw) : (aNbFaces == 0)`.
        let incomplete = if nb_trous > 0 {
            nb_faces < nb_law
        } else {
            nb_faces == 0
        };
        if incomplete {
            self.is_done = false;
            return Ok(());
        }

        // `cxx:3348-3357`: the auxiliary faces.
        for face in &self.my_aux_shape {
            if !face.is_null() && face.shape_type() == ShapeType::Face {
                b.add(&mut shell.0, face);
            }
        }

        // `cxx:3358-3389`: substitution of the modified edges.
        let mut a_substitute = BRepToolsSubstitution::new();
        for (key, new_edge) in &self.my_v_edges_modified {
            a_substitute.substitute_at_key(key.clone(), &[new_edge.clone()]);
        }
        a_substitute.build(&shell.0);
        if a_substitute.is_copied(&shell.0) {
            if let Some(list) = a_substitute.copy(&shell.0) {
                if let Some(first) = list.first() {
                    shell = Shell(first.clone());
                }
            }
        }

        // `cxx:3391-3446`: rewrite the tables with the substituted shapes.
        for ii in 1..=self.my_faces.nb_rows() {
            for jj in 1..=self.my_faces.nb_cols() {
                let Some(local) = self.my_faces.value(ii, jj).cloned() else {
                    continue;
                };
                if a_substitute.is_copied(&local) {
                    if let Some(first) = a_substitute.copy(&local).and_then(|l| l.first()).cloned() {
                        self.my_faces.set_value(ii, jj, Some(first));
                    }
                }
            }
        }
        for ii in 1..=self.my_v_edges.nb_rows() {
            for jj in 1..=self.my_v_edges.nb_cols() {
                let Some(local) = self.my_v_edges.value(ii, jj).cloned() else {
                    continue;
                };
                if a_substitute.is_copied(&local) {
                    if let Some(first) = a_substitute.copy(&local).and_then(|l| l.first()).cloned() {
                        self.my_v_edges.set_value(ii, jj, Some(first));
                    }
                }
            }
        }
        for ii in 1..=self.my_u_edges.nb_rows() {
            for jj in 1..=self.my_u_edges.nb_cols() {
                let Some(local) = self.my_u_edges.value(ii, jj).cloned() else {
                    continue;
                };
                if a_substitute.is_copied(&local) {
                    if let Some(first) = a_substitute.copy(&local).and_then(|l| l.first()).cloned() {
                        self.my_u_edges.set_value(ii, jj, Some(first));
                    }
                }
            }
        }

        // `cxx:3448-3490`: ensure same parameter on the U-edges.
        let u_closed = self.my_sec.base().is_u_closed();
        for ii in 1..=self.my_u_edges.nb_rows() {
            if u_closed && ii == self.my_u_edges.nb_rows() {
                continue;
            }
            for jj in 1..=self.my_u_edges.nb_cols() {
                let Some(shape) = self.my_u_edges.value(ii, jj).cloned() else {
                    continue;
                };
                let mut an_edge = Edge(shape);
                if BRepTool::is_degenerated(&an_edge) {
                    continue;
                }
                let mut i1 = ii as i64 - 1;
                let mut i2 = ii as i64;
                if i1 == 0 && u_closed {
                    i1 = self.my_faces.nb_rows() as i64;
                }
                if i2 > self.my_faces.nb_rows() as i64 {
                    i2 = 0;
                }
                let face1 = if i1 != 0 {
                    self.my_faces
                        .value(i1 as usize, jj)
                        .filter(|s| s.shape_type() == ShapeType::Face)
                        .map(|s| Face(s.clone()))
                } else {
                    None
                };
                let face2 = if i2 != 0 {
                    self.my_faces
                        .value(i2 as usize, jj)
                        .filter(|s| s.shape_type() == ShapeType::Face)
                        .map(|s| Face(s.clone()))
                } else {
                    None
                };
                if let (Some(f1), Some(f2)) = (face1, face2) {
                    correct_same_parameter(&mut an_edge, Some(&f1), Some(&f2));
                }
            }
        }

        // `cxx:3492-3503`: the tapes.
        for ii in 1..=nb_law {
            for jj in 1..=nb_path {
                let Some(face) = self.my_faces.value(ii, jj) else {
                    continue;
                };
                if face.is_null() || face.shape_type() != ShapeType::Face {
                    continue;
                }
                if let Some(shell_tape) = self.my_tapes.get_mut(ii - 1).and_then(|t| t.as_mut()) {
                    b.add(shell_tape, face);
                }
            }
        }

        // `cxx:3505-3527`: is it closed?
        if self.my_loc.base().is_closed() && u_closed {
            let efmap = map_shapes_and_ancestors(&shell.0, ShapeType::Edge, ShapeType::Face);
            let mut closed = true;
            for e in edges_of(&shell.0) {
                if BRepTool::is_degenerated(&e) {
                    continue;
                }
                let extent = efmap.get(&shape_map_key(&e)).map(|l| l.len()).unwrap_or(0);
                closed = extent > 1;
                if !closed {
                    break;
                }
            }
            shell.0.set_closed(closed);
        }

        self.my_shape = Some(shell.0);
        Ok(())
    }
}
