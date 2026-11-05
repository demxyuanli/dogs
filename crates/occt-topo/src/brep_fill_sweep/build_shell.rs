//! `BRepFill_Sweep::BuildShell` (`BRepFill_Sweep.cxx:2213-3208`): the swept
//! shell of a non-vertex section law, one path window `[IFirst, ILast]` at a
//! time.
//!
//! The reference is one 1,000-line function; the port keeps its section order
//! and its local tables verbatim, with the sections that only touch data split
//! into sibling modules:
//!
//! * (1) surfaces — [`BrepFillSweep::build_shell`] itself.
//! * (2.0) pre-existing edges / vertices and the two boundary vertex columns —
//!   [`frontier_vertices`] (`cxx:2353-2662`).
//! * (2.1)-(2.3) new vertices, iso edges and the singular cases —
//!   [`frontier_edges`] (`cxx:2665-2989`).
//! * (3)-(5) faces, history / continuity and the tape-rail update —
//!   [`shell_topology`] (`cxx:2992-3206`).
//!
//! Translation notes:
//! - `int ipath` is the 1-based *window* index and `IPath` the 1-based index in
//!   the whole path (`IPath = IFirst + ipath - 1`); the port keeps both names.
//! - `NCollection_Map<TopoDS_Shape>` is `HashSet<ShapeMapKey>` and
//!   `NCollection_DataMap<TopoDS_Shape, handle(HArray2)>` is
//!   `HashMap<ShapeMapKey, ShapeTable>` (`Tapes` rows 1-6, `Rails` rows 1-2),
//!   the same translation `myVEdgesModified` already uses.
//! - `GeomLib::ExtendSurfByLength` is
//!   [`occt_geom::geom_lib_extend_surf::extend_surf_by_length`]. The reference
//!   down-casts to `Geom_BoundedSurface` before the call (`cxx:2333`); when that
//!   down-cast fails it passes a null handle, which the port reports as an
//!   unported branch instead of dereferencing it.
//! - A null `TopoDS_Shape` is `None`; `S->Value(u, v)` is [`Surface::value`].

use std::collections::{HashMap, HashSet};

use occt_core::gp::GpPnt;
use occt_core::toploc::TopLocLocation;
use occt_geom::geom_fill::sweep::Sweep as GeomFillSweep;
use occt_geom::geom_lib_extend_surf::extend_surf_by_length;

use crate::abs::Orientation;
use crate::builder::TopoBuilder;
use crate::shape::{TopoShape, Vertex};
use crate::topo_tools::{edge_vertices, shape_map_key, ShapeMapKey};

use super::helpers::{self, ShapeTable};
use super::shell_tables::{ShellLocals, ShellTables};
use super::{BrepFillSweep, TransitionStyle};

mod frontier_edges;
mod frontier_vertices;
mod shell_topology;

/// `IsG1(0, myTol3d)`'s default `AngularTolerance`
/// (`BRepFill_LocationLaw.hxx:102`).
const DEFAULT_ANGULAR_TOL: f64 = 1.0e-4;

/// The state `BuildShell` threads through its sections: the sweep being filled,
/// the tables of `cxx:2242-2252`, the `NCollection_Array1` locals and the flags
/// `cxx:2231-2237` compute.
pub(super) struct ShellBuild<'a> {
    /// The reference's `this`.
    pub sweep: &'a mut BrepFillSweep,
    /// `ExchUV`, `UReverse`, `Degenerated`, `TabErr`, `TabS`, `UEdge`,
    /// `VEdge`, `Vertex` and `myFaces`.
    pub tables: ShellTables,
    /// `SecVertex`, `VError`, `Vi`, `IsBuilt` and `StartEdges`.
    pub locals: ShellLocals,
    /// `NbLaw = mySec->NbLaw()`.
    pub nb_law: usize,
    /// `NbPath = ILast - IFirst`.
    pub nb_path: usize,
    /// `myLoc->NbLaw()`.
    pub nb_path_law: usize,
    /// `IFirst`.
    pub i_first: usize,
    /// `ILast`.
    pub i_last: usize,
    /// `constSection`.
    pub const_section: bool,
    /// `uclose`.
    pub u_close: bool,
    /// `vclose`.
    pub v_close: bool,
    /// `global_vclose`.
    pub global_vclose: bool,
    /// `hasdegen`.
    pub hasdegen: bool,
    /// `ReversedEdges`.
    pub reversed_edges: &'a mut HashSet<ShapeMapKey>,
    /// `Tapes`.
    pub tapes: &'a mut HashMap<ShapeMapKey, ShapeTable>,
    /// `Rails`.
    pub rails: &'a mut HashMap<ShapeMapKey, ShapeTable>,
}

impl BrepFillSweep {
    /// `BRepFill_Sweep::BuildShell(Transition, IFirst, ILast, ReversedEdges,
    /// Tapes, Rails, ExtendFirst, ExtendLast)` (`cxx:2213-3208`).
    ///
    /// `Ok(false)` is the reference's `return false` (a failed sweep surface or
    /// a null edge that a later `TopoDS::Edge` cast would dereference);
    /// `Err` is a branch the port cannot follow, with the OCCT site named.
    #[allow(clippy::too_many_arguments)]
    pub fn build_shell(
        &mut self,
        _transition: TransitionStyle,
        i_first: usize,
        i_last: usize,
        reversed_edges: &mut HashSet<ShapeMapKey>,
        tapes: &mut HashMap<ShapeMapKey, ShapeTable>,
        rails: &mut HashMap<ShapeMapKey, ShapeTable>,
        extend_first: f64,
        extend_last: f64,
    ) -> Result<bool, String> {
        let nb_law = self.my_sec.base().nb_law();
        let nb_path_law = self.my_loc.base().nb_law();
        if i_last <= i_first {
            return Err(format!("BuildShell: empty window [{i_first}, {i_last}]"));
        }
        let nb_path = i_last - i_first;

        // `cxx:2231-2237`.
        let const_section = self.my_sec.is_constant();
        let u_close = self.my_sec.base().is_u_closed();
        let global_vclose = self.my_loc.base().is_closed()
            && self.my_loc.base().is_g1(0, self.my_tol_3d, DEFAULT_ANGULAR_TOL) >= 0;
        let v_close = global_vclose
            && self.my_sec.base().is_v_closed()
            && nb_path == nb_path_law;
        self.error = 0.0;

        let mut build = ShellBuild {
            sweep: self,
            tables: ShellTables::new(nb_law, nb_path),
            locals: ShellLocals::new(nb_law, nb_path),
            nb_law,
            nb_path,
            nb_path_law,
            i_first,
            i_last,
            const_section,
            u_close,
            v_close,
            global_vclose,
            hasdegen: false,
            reversed_edges,
            tapes,
            rails,
        };

        // (1.1) Construction of Tables (`cxx:2263-2288`): `CurvilinearBounds`
        // fills `SecDom` with the abscissa, which `Law(1)->GetDomain` then
        // overwrites on the next line of the reference.
        let (_abscissa, length) = build.sweep.my_loc.base().curvilinear_bounds(nb_path_law)?;
        let sec_law1 = build
            .sweep
            .my_sec
            .base().law(1)
            .ok_or_else(|| "BuildShell: null section law 1".to_string())?;
        let (sec_deb, domain_last) = sec_law1.get_domain();
        let sec_dom = domain_last - sec_deb;
        if i_first > 1 {
            // `cxx:2269-2274`.
            let (_lf, ll) = build.sweep.my_loc.base().curvilinear_bounds(i_first - 1)?;
            build.locals.vi[1] = sec_deb + (ll / length) * sec_dom;
        } else {
            build.locals.vi[1] = sec_deb;
        }

        // Error a priori on vertices (`cxx:2280-2287`).
        if const_section {
            for isec in 1..=nb_law + 1 {
                build.locals.v_error[isec] = build.sweep.my_sec.vertex_tol(isec - 1, 0.0);
                build.locals.sec_vertex[isec] = build
                    .sweep
                    .my_sec
                    .vertex(isec, 0.0)
                    .map(|v| v.0.clone());
            }
        }

        // (1.2) Calculate surfaces (`cxx:2290-2346`).
        if !build.build_surfaces(sec_deb, sec_dom, length, extend_first, extend_last)? {
            return Ok(false);
        }

        // (2.0) return preexisting Edges and vertices (`cxx:2353-2480`).
        build.register_preexisting_edges()?;

        // Vertex columns of the first and last section (`cxx:2482-2662`).
        build.construct_boundary_vertices()?;

        // (2.1)-(2.3) vertices and edges of the sweep (`cxx:2665-2989`).
        if !build.construct_vertices_and_edges()? {
            return Ok(false);
        }

        // (3) faces, (4) history and (5) Tapes / Rails (`cxx:2992-3206`).
        build.build_faces();
        build.reverse_built_faces();
        build.encode_history();
        build.update_tapes();

        Ok(true)
    }
}

impl ShellBuild<'_> {
    /// (1.2) `cxx:2290-2346`: one `GeomFill_Sweep` per path law, built once per
    /// section law.
    fn build_surfaces(
        &mut self,
        sec_deb: f64,
        sec_dom: f64,
        length: f64,
        extend_first: f64,
        extend_last: f64,
    ) -> Result<bool, String> {
        for ipath in 1..=self.nb_path {
            let ipath_global = self.i_first + ipath - 1;
            let Some(law) = self.sweep.my_loc.base().law(ipath_global) else {
                return Ok(false);
            };
            let mut sweep = GeomFillSweep::with_k_part(law.clone(), self.sweep.k_part);
            sweep.set_tolerance(
                self.sweep.my_tol_3d,
                self.sweep.my_bound_tol,
                self.sweep.my_tol_2d,
                self.sweep.my_tol_angular,
            );
            sweep.set_force_approx_c1(self.sweep.my_force_approx_c1);

            // `cxx:2298-2311`: parametric correspondence.
            let (lf, ll) = law.get_domain();
            let (_l_f, l_l) = self.sweep.my_loc.base().curvilinear_bounds(ipath_global)?;
            self.locals.vi[ipath + 1] = sec_deb + (l_l / length) * sec_dom;
            if !self.const_section {
                sweep.set_domain(
                    lf,
                    ll,
                    self.locals.vi[ipath],
                    self.locals.vi[ipath + 1],
                );
            }

            for isec in 1..=self.nb_law {
                let Some(section) = self.sweep.my_sec.base().law(isec) else {
                    return Ok(false);
                };
                sweep.build(
                    section,
                    self.sweep.my_approx_style,
                    self.sweep.my_continuity,
                    self.sweep.my_degmax,
                    self.sweep.my_segmax,
                );
                if !sweep.is_done() {
                    return Ok(false);
                }
                let Some(surface) = sweep.surface() else {
                    return Ok(false);
                };
                self.tables.tab_s.set_value(isec, ipath, surface);
                let err = sweep.error_on_surface();
                self.tables.tab_err.set_value(isec, ipath, err);
                self.tables.exch_uv.set_value(isec, ipath, sweep.exchange_uv());
                self.tables
                    .u_reverse
                    .set_value(isec, ipath, sweep.u_reversed());
                if err > self.sweep.error {
                    self.sweep.error = err;
                }

                // `cxx:2330-2345`: extend the first / last sweep surface.
                let extend = if ipath == 1 && extend_first > 0.0 {
                    Some((extend_first, false))
                } else if ipath == self.nb_path && extend_last > 0.0 {
                    Some((extend_last, true))
                } else {
                    None
                };
                if let Some((len, after)) = extend {
                    let Some(surf) = self.tables.tab_s.value(isec, ipath).cloned() else {
                        continue;
                    };
                    let Some(extended) =
                        extend_surf_by_length(&surf, len, 1, sweep.exchange_uv(), after)
                    else {
                        // `occ::down_cast<Geom_BoundedSurface>` failed at
                        // `cxx:2333` / `cxx:2340`: the reference hands a null
                        // handle to `GeomLib::ExtendSurfByLength`, so this
                        // branch is UNPORTED rather than guessed.
                        return Err(format!(
                            "BuildShell: ExtendSurfByLength on a non-bounded surface (BRepFill_Sweep.cxx:2333, isec={isec}, ipath={ipath})"
                        ));
                    };
                    self.tables.tab_s.set_value(isec, ipath, extended);
                }
            }
        }
        Ok(true)
    }

    /// (2.0) `cxx:2353-2480`: the pre-existing edges of `FirstShape`, their
    /// `Tapes` / `Rails` replay, and the first vertex column.
    fn register_preexisting_edges(&mut self) -> Result<(), String> {
        // `cxx:2352-2354`.
        if let Some(w) = self.sweep.first_shape.as_mut() {
            helpers::reverse_modified_edges_set(&mut w.0, self.reversed_edges);
        }
        if let Some(w) = self.sweep.last_shape.as_mut() {
            helpers::reverse_modified_edges_set(&mut w.0, self.reversed_edges);
        }

        if self.sweep.first_shape.is_none() || self.i_first != 1 {
            return Ok(());
        }
        let first_shape = self
            .sweep
            .first_shape
            .clone()
            .ok_or_else(|| "BuildShell: null FirstShape".to_string())?;

        // `cxx:2363`: `mySec->Init(FirstShape)`.
        self.sweep.my_sec.base().init(&first_shape);
        let mut v_first_isec = None;
        let mut v_last_isec = None;

        for isec in 1..=self.nb_law {
            // `cxx:2366`: the current edge of the section wire.
            let Some(edge) = self.sweep.my_sec.base().current_edge() else {
                return Err(format!(
                    "BuildShell: FirstShape has no edge {isec} (BRepFill_Sweep.cxx:2366)"
                ));
            };
            let (v_first, v_last) = oriented_edge_vertices(&edge);
            self.tables.v_edge.set_value(isec, 1, Some(edge.0.clone()));
            // `cxx:2370-2378`.
            let boundary = if edge.0.orientation() == Orientation::Reversed {
                v_first.clone().map(|v| v.0)
            } else {
                v_last.clone().map(|v| v.0)
            };
            if let Some(v) = boundary {
                self.tables.vertex.set_value(isec + 1, 1, Some(v.clone()));
                self.update_vertex_cell(
                    self.i_first - 1,
                    isec + 1,
                    &v,
                    self.tables.tab_err.value(isec, 1),
                    self.locals.vi[1],
                );
            }
            if isec == 1 {
                v_first_isec = v_first.clone().map(|v| v.0);
                v_last_isec = v_last.clone().map(|v| v.0);
            }
            self.locals.start_edges[isec] = Some(edge.0.clone());

            let key = shape_map_key(&edge.0);
            // `Tapes.IsBound(E)` (`cxx:2382`).
            if let Some(tape) = self.tapes.get(&key).cloned() {
                self.locals.is_built[isec] = true;
                for j in 1..=self.nb_path + 1 {
                    self.tables.v_edge.set_value(isec, j, tape.value(1, j).cloned());
                    // `VEdge(isec, j).Reverse()` (`cxx:2390`).
                    self.tables.v_edge.reverse_value(isec, j);
                }
                let ifirst = isec + 1;
                let ilast = isec;
                for j in 1..=self.nb_path {
                    self.tables.u_edge.set_value(ifirst, j, tape.value(2, j).cloned());
                }
                for j in 1..=self.nb_path {
                    self.tables.u_edge.set_value(ilast, j, tape.value(3, j).cloned());
                }
                for j in 1..=self.nb_path + 1 {
                    self.tables.vertex.set_value(ifirst, j, tape.value(4, j).cloned());
                }
                for j in 1..=self.nb_path + 1 {
                    self.tables.vertex.set_value(ilast, j, tape.value(5, j).cloned());
                }
                for j in 1..=self.nb_path {
                    self.tables.faces.set_value(isec, j, tape.value(6, j).cloned());
                }

                // `cxx:2413-2436`.
                if self.u_close && isec == 1 {
                    for j in 1..=self.nb_path {
                        self.tables
                            .u_edge
                            .set_value(self.nb_law + 1, j, self.tables.u_edge.value(1, j).cloned());
                    }
                    for j in 1..=self.nb_path + 1 {
                        self.tables
                            .vertex
                            .set_value(self.nb_law + 1, j, self.tables.vertex.value(1, j).cloned());
                    }
                }
                if self.u_close && isec == self.nb_law {
                    for j in 1..=self.nb_path {
                        self.tables
                            .u_edge
                            .set_value(1, j, self.tables.u_edge.value(self.nb_law + 1, j).cloned());
                    }
                    for j in 1..=self.nb_path + 1 {
                        self.tables
                            .vertex
                            .set_value(1, j, self.tables.vertex.value(self.nb_law + 1, j).cloned());
                    }
                }
            } else {
                // `cxx:2437-2467`: bind the empty tape and replay the rails.
                self.tapes
                    .insert(key, ShapeTable::new(6, self.nb_path + 1));
                if let Some(vf) = &v_first {
                    let ind = if edge.0.orientation() == Orientation::Reversed {
                        isec + 1
                    } else {
                        isec
                    };
                    if let Some(rail) = self.rails.get(&shape_map_key(&vf.0)).cloned() {
                        for j in 1..=self.nb_path {
                            self.tables.u_edge.set_value(ind, j, rail.value(1, j).cloned());
                        }
                        for j in 1..=self.nb_path + 1 {
                            self.tables.vertex.set_value(ind, j, rail.value(2, j).cloned());
                        }
                    }
                }
                if let Some(vl) = &v_last {
                    let ind = if edge.0.orientation() == Orientation::Forward {
                        isec + 1
                    } else {
                        isec
                    };
                    if let Some(rail) = self.rails.get(&shape_map_key(&vl.0)).cloned() {
                        for j in 1..=self.nb_path {
                            self.tables.u_edge.set_value(ind, j, rail.value(1, j).cloned());
                        }
                        for j in 1..=self.nb_path + 1 {
                            self.tables.vertex.set_value(ind, j, rail.value(2, j).cloned());
                        }
                    }
                }
            }
        }

        // `cxx:2470-2479`.
        let Some(v_edge_1_1) = self.tables.v_edge.value(1, 1).cloned() else {
            return Ok(());
        };
        let edge = crate::shape::Edge(v_edge_1_1);
        let boundary = if edge.0.orientation() == Orientation::Reversed {
            v_last_isec
        } else {
            v_first_isec
        };
        let Some(boundary) = boundary else {
            return Ok(());
        };
        self.tables.vertex.set_value(1, 1, Some(boundary.clone()));
        self.update_vertex_cell(
            self.i_first - 1,
            1,
            &boundary,
            self.tables.tab_err.value(1, 1),
            self.locals.vi[1],
        );
        Ok(())
    }
}

/// `TopExp::Vertices(E, Vfirst, Vlast)` with `CumOri = true`
/// (`TopExp.cxx:214-253`): the stored pair swapped when the edge is `REVERSED`.
pub(super) fn oriented_edge_vertices(edge: &crate::shape::Edge) -> (Option<Vertex>, Option<Vertex>) {
    let (first, last) = edge_vertices(edge);
    if edge.0.orientation() == Orientation::Reversed {
        (last, first)
    } else {
        (first, last)
    }
}

/// `TopLoc_Location Identity; V.Location(Identity)` followed by
/// `BRep_Builder::UpdateVertex(V, P, Tol)` (`BRep_Builder.cxx:1096-1111`): the
/// stored point is replaced and the tolerance raised.
pub(super) fn reset_location_and_move(shape: &mut TopoShape, p: &GpPnt, tol: f64) {
    shape.set_location(&TopLocLocation::identity());
    crate::brep_fill_offset_wire::vertices::update_vertex_point(&Vertex(shape.clone()), p, tol);
}

/// `BRep_Builder::MakeVertex(TopoDS::Vertex(Vertex(isec, col)), P, Tol)` on a
/// null table cell (`cxx:2537-2541`, `cxx:2649-2653`).
pub(super) fn make_vertex_at(p: &GpPnt, tol: f64) -> TopoShape {
    TopoBuilder::new().make_vertex(*p, tol).0
}
