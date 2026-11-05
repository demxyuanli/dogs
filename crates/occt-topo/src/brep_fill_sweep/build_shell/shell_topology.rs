//! `BRepFill_Sweep::BuildShell` sections (3)-(5) (`cxx:2992-3206`): the natural
//! faces of the sweep, the early-built faces' orientation, the history /
//! continuity encoding and the `Tapes` / `Rails` update.
//!
//! Translation notes:
//! - `BuildFace(TabS, UEdge, VEdge, UEdge, VEdge, myVEdgesModified, ExchUV,
//!   UReverse, face)` is [`super::build_face::build_face`]; its `EEmap`
//!   argument is the sweep's `myVEdgesModified`, and a failed build leaves the
//!   cell null exactly like the reference's untouched output face.
//! - `myFaces`, `myUEdges` and `myVEdges` are the sweep's members
//!   (`hxx:180-183`); `BuildShell` fills them through the working tables and
//!   the port copies them across at the end of each section.
//! - `BRep_Builder::Continuity(E, F1, F2, GeomAbs_G1)` is
//!   [`crate::brep_lib_continuity::set_stored_continuity`] with the `GeomAbs`
//!   code `1`, and `BRepLib::EncodeRegularity` is
//!   [`crate::brep_lib_continuity::encode_regularity_edge`] /
//!   [`crate::brep_lib_continuity::encode_regularity`].
//! - `TopExp::MapShapesAndAncestors` and `shell.Closed(closed)` live in
//!   `Build` (`cxx:3502-3526`), not here.

use crate::abs::ShapeType;
use crate::brep_lib_continuity::{encode_regularity, encode_regularity_edge, set_stored_continuity};
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{is_same, shape_map_key};

use super::ShellBuild;
use crate::brep_fill_sweep::build_face::build_face;

impl ShellBuild<'_> {
    /// (3) `cxx:2995-3022`: one face per non-degenerated cell, plus the
    /// degenerated cells' pass-through.
    pub(super) fn build_faces(&mut self) {
        for ipath in 1..=self.nb_path {
            let ipath_global = self.i_first + ipath - 1;
            for isec in 1..=self.nb_law {
                if self.tables.degenerated.value(isec, ipath) {
                    // `cxx:2999-3008`.
                    let u_here = self.tables.u_edge.value(isec, ipath).cloned();
                    let u_next = self.tables.u_edge.value(isec + 1, ipath).cloned();
                    let same = match (&u_here, &u_next) {
                        (Some(a), Some(b)) => is_same(a, b),
                        (None, None) => true,
                        _ => false,
                    };
                    let value = if same {
                        u_here
                    } else {
                        self.tables.v_edge.value(isec, ipath).cloned()
                    };
                    self.tables.faces.set_value(isec, ipath_global, value);
                } else if self.tables.faces.is_null(isec, ipath_global) {
                    // `cxx:3009-3021`.
                    let Some(s) = self.tables.tab_s.value(isec, ipath).cloned() else {
                        continue;
                    };
                    let edges = [
                        self.tables.u_edge.value(isec, ipath).cloned(),
                        self.tables.v_edge.value(isec, ipath).cloned(),
                        self.tables.u_edge.value(isec + 1, ipath).cloned(),
                        self.tables.v_edge.value(isec, ipath + 1).cloned(),
                    ];
                    let (Some(e1), Some(e2), Some(e3), Some(e4)) = (
                        edges[0].clone().map(Edge),
                        edges[1].clone().map(Edge),
                        edges[2].clone().map(Edge),
                        edges[3].clone().map(Edge),
                    ) else {
                        continue;
                    };
                    let face = build_face(
                        &s,
                        &e1,
                        &e2,
                        &e3,
                        &e4,
                        &mut self.sweep.my_v_edges_modified,
                        self.tables.exch_uv.value(isec, ipath),
                        self.tables.u_reverse.value(isec, ipath),
                    );
                    self.tables
                        .faces
                        .set_value(isec, ipath_global, face.map(|f| f.0));
                }
            }
        }
        // `myFaces = the working table` (`hxx:183`).
        self.sweep.my_faces = self.tables.faces.clone();
    }

    /// (3.1) `cxx:3027-3035`: reverse the faces that came from a tape.
    pub(super) fn reverse_built_faces(&mut self) {
        for ipath in 1..=self.nb_path {
            for isec in 1..=self.nb_law {
                if self.locals.is_built[isec] {
                    self.tables.faces.reverse_value(isec, ipath);
                }
            }
        }
        self.sweep.my_faces = self.tables.faces.clone();
    }

    /// (4) `cxx:3037-3139`: history tables, and either the compound-wide
    /// regularity encoding of a degenerated sweep or the per-edge `G1` encoding
    /// of the general case.
    pub(super) fn encode_history(&mut self) {
        let tol_angular = self.sweep.my_tol_angular;
        if self.hasdegen {
            // (4.1) `cxx:3040-3062`: sledgehammer.
            let b = TopoBuilder::new();
            let mut comp = TopoShape::new(ShapeType::Compound);
            for isec in 1..=self.nb_law + 1 {
                for ipath in 1..=self.nb_path + 1 {
                    let ipath_global = self.i_first + ipath - 1;
                    if ipath <= self.nb_path {
                        self.sweep.my_u_edges.set_value(
                            isec,
                            ipath_global,
                            self.tables.u_edge.value(isec, ipath).cloned(),
                        );
                    }
                    if isec <= self.nb_law {
                        self.sweep.my_v_edges.set_value(
                            isec,
                            ipath_global,
                            self.tables.v_edge.value(isec, ipath).cloned(),
                        );
                    }
                    if ipath <= self.nb_path && isec <= self.nb_law {
                        if let Some(face) = self.sweep.my_faces.value(isec, ipath_global) {
                            if !face.is_null() && face.shape_type() == ShapeType::Face {
                                b.add(&mut comp, face);
                            }
                        }
                    }
                }
            }
            encode_regularity(&comp, tol_angular);
        } else {
            // (4.2) `cxx:3066-3138`: tweezers.
            for isec in 1..=self.nb_law + 1 {
                let is_g1 = isec > 1 && self.sweep.my_sec.continuity(isec - 1, tol_angular) >= 1;
                for ipath in 1..=self.nb_path {
                    let ipath_global = self.i_first + ipath - 1;
                    self.sweep.my_u_edges.set_value(
                        isec,
                        ipath_global,
                        self.tables.u_edge.value(isec, ipath).cloned(),
                    );
                    if !is_g1 {
                        continue;
                    }
                    let ff_row = if isec == self.nb_law + 1 { 1 } else { isec };
                    let ff = self.sweep.my_faces.value(ff_row, ipath_global).cloned();
                    let f1 = self.sweep.my_faces.value(isec - 1, ipath_global).cloned();
                    let edge = self.sweep.my_u_edges.value(isec, ipath_global).cloned();
                    if let (Some(e), Some(f1), Some(ff)) = (edge, f1, ff) {
                        if f1.shape_type() == ShapeType::Face && ff.shape_type() == ShapeType::Face
                        {
                            // `B.Continuity(E, F1, FF, GeomAbs_G1)` (`cxx:3097-3101`).
                            set_stored_continuity(&Edge(e), &Face(f1), &Face(ff), 1);
                        }
                    }
                }
            }

            let nbpath = if self.v_close {
                self.nb_path + 1
            } else {
                self.nb_path
            };
            for ipath in 1..=self.nb_path + 1 {
                let ipath_global = self.i_first + ipath - 1;
                let is_g1 = ipath > 1
                    && ipath <= nbpath
                    && self
                        .sweep
                        .my_loc
                        .base()
                        .is_g1(ipath_global - 1, self.sweep.my_tol_3d, tol_angular)
                        >= 0;
                for isec in 1..=self.nb_law {
                    self.sweep.my_v_edges.set_value(
                        isec,
                        ipath_global,
                        self.tables.v_edge.value(isec, ipath).cloned(),
                    );
                    if !is_g1 {
                        continue;
                    }
                    let ff_col = if ipath == self.nb_path + 1 { 1 } else { ipath_global };
                    let ff = self.sweep.my_faces.value(isec, ff_col).cloned();
                    let f2 = self.sweep.my_faces.value(isec, ipath_global - 1).cloned();
                    let edge = self.sweep.my_v_edges.value(isec, ipath_global).cloned();
                    if let (Some(e), Some(f2), Some(ff)) = (edge, f2, ff) {
                        if f2.shape_type() == ShapeType::Face && ff.shape_type() == ShapeType::Face
                        {
                            encode_regularity_edge(
                                &Edge(e),
                                &Face(ff),
                                &Face(f2),
                                tol_angular,
                            );
                        }
                    }
                }
            }
        }
    }

    /// (5) `cxx:3141-3204`: write the tapes and rails of the replayed section
    /// edges back into the caller's maps.
    pub(super) fn update_tapes(&mut self) {
        if self.i_first != 1 || self.tapes.is_empty() {
            return;
        }
        for isec in 1..=self.nb_law {
            let Some(start_edge) = self.locals.start_edges[isec].clone() else {
                continue;
            };
            let Some(array) = self.tapes.get_mut(&shape_map_key(&start_edge)) else {
                continue;
            };
            for j in 1..=self.nb_path + 1 {
                array.set_value(1, j, self.sweep.my_v_edges.value(isec, j).cloned());
            }
            for j in 1..=self.nb_path {
                array.set_value(2, j, self.sweep.my_u_edges.value(isec, j).cloned());
            }
            for j in 1..=self.nb_path {
                array.set_value(3, j, self.sweep.my_u_edges.value(isec + 1, j).cloned());
            }
            for j in 1..=self.nb_path + 1 {
                array.set_value(4, j, self.tables.vertex.value(isec, j).cloned());
            }
            for j in 1..=self.nb_path + 1 {
                array.set_value(5, j, self.tables.vertex.value(isec + 1, j).cloned());
            }
            for j in 1..=self.nb_path {
                array.set_value(6, j, self.sweep.my_faces.value(isec, j).cloned());
            }

            // `cxx:3172-3202`: bind the rails of both endpoints.
            let edge = Edge(start_edge);
            let (v_first, v_last) = super::oriented_edge_vertices(&edge);
            if let Some(vf) = v_first {
                let key = shape_map_key(&vf.0);
                if !self.rails.contains_key(&key) {
                    let mut rail = crate::brep_fill_sweep::helpers::ShapeTable::new(2, self.nb_path + 1);
                    for j in 1..=self.nb_path {
                        rail.set_value(1, j, self.sweep.my_u_edges.value(isec, j).cloned());
                    }
                    for j in 1..=self.nb_path + 1 {
                        rail.set_value(2, j, self.tables.vertex.value(isec, j).cloned());
                    }
                    self.rails.insert(key, rail);
                }
            }
            if let Some(vl) = v_last {
                let key = shape_map_key(&vl.0);
                if !self.rails.contains_key(&key) {
                    let mut rail = crate::brep_fill_sweep::helpers::ShapeTable::new(2, self.nb_path + 1);
                    for j in 1..=self.nb_path {
                        rail.set_value(1, j, self.sweep.my_u_edges.value(isec + 1, j).cloned());
                    }
                    for j in 1..=self.nb_path + 1 {
                        rail.set_value(2, j, self.tables.vertex.value(isec + 1, j).cloned());
                    }
                    self.rails.insert(key, rail);
                }
            }
        }
    }
}
