//! `ChFi3d_Builder::Compute` - the topological orchestration of a fillet or
//! chamfer computation.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder.cxx:178-680`
//! (OCCT 8.0.0).
//!
//! The body is a straight-line pipeline over the contours the subclass loaded:
//! `UpdateTolesp` -> `Reset` -> `myVDataMap` -> `ExtentAnalyse` ->
//! `PerformSetOfSurf` on every contour -> `PerformFilletOnVertex` on every
//! corner vertex -> DS indexing of the solids / shells -> `ChFi3d_StripeEdgeInter`
//! cross-check + `ChFi3d_FilDS` per contour -> `CompleteDS` -> the curve / vertex
//! tolerance pass -> `myCoup->Perform(myDS)` + `MergeSolid` -> the split-edge
//! vertex tolerances -> `myShapeResult` / `badShape` -> `SetRegul` -> the
//! `SameParameter` pass over the new faces.
//!
//! Divergences that are not heuristics:
//! - OCCT wraps `PerformSetOfSurf` and `PerformFilletOnVertex` in `try` /
//!   `catch (Standard_Failure)` and records the failure in `badstripes` /
//!   `badvertices`; the port answers those two calls with `Err`, which is the
//!   same control flow, while the failures OCCT lets escape (`ChFi3d_FilDS`,
//!   the empty contour list) are propagated with `?`.
//! - `TopOpeBRepDS_CurveExplorer` walks the bound curves in index order and
//!   `IsCurve(I)` is the `is_curve_bound` test, so the loop below is the same
//!   walk without the explorer's borrow of the DS (the body mutates it).
//! - `NCollection_Map<int>` has no order; the port keeps `MapIndSo` as a `Vec`
//!   in insertion order, so the members of the result compound follow the
//!   order the solids / shells were found.
//! - `BRep_Builder::MakeCompound` on `myShapeResult` / `badShape` becomes a
//!   `Vec<TopoShape>` plus `make_compound_of`; a null `TopoDS_Shape` member has
//!   no counterpart, so those members are `Option<TopoShape>`.

use occt_core::precision::CONFUSION;

use crate::abs::{ShapeType, State as TopAbsState};
use crate::brep_lib_same_parameter::same_parameter_shape;
use crate::brep_tool::BRepTool;
use crate::chfi3d::builder_0_helpers::ch_fi3d_stripe_edge_inter;
use crate::chfi_ds::{ErrorStatus, State as SpineState};
use crate::shape::{Edge, TopoShape, Vertex};
use crate::shhealing::shape_fix_same_parameter;
use crate::topo_tools::{explorer_shapes, explorer_shapes_avoiding};
use crate::topopebrep::ds::{Kind, PointIterator};

use super::complete_ds::complete_ds;
use super::{spine_of, ChFi3dBuilder, ChFi3dCornerOps, ChFi3dSurfOps};

/// `BRep_Builder::UpdateVertex(V, Tol)` (`BRep_Builder.cxx:1442-1452`):
/// `BRep_TVertex::UpdateTolerance` is a max (`BRep_TVertex.hxx:50-52`).
fn update_vertex_tolerance(v: &Vertex, tol: f64) {
    v.set_tolerance(BRepTool::vertex_tolerance(v).max(tol));
}

/// `BRep_Builder::MakeCompound(Comp)` + the `B1.Add` calls of `cxx:512-566`.
fn make_compound(members: &[TopoShape]) -> TopoShape {
    crate::builder::TopoBuilder::new()
        .make_compound_of(members)
        .0
}

impl ChFi3dBuilder {
    /// `Compute()` (`cxx:178-680`).
    ///
    /// `Err` carries the failures OCCT lets escape `Compute`: the
    /// `Standard_Failure("There are no suitable edges for chamfer or fillet")`
    /// of `cxx:225-228`, the `throw` of `ChFi3d_FilDS` (`cxx:391`) and the
    /// `Standard_OutOfRange` / null-shape dereferences the port answers as
    /// errors instead of crashing.
    pub fn compute(
        &mut self,
        ops: &mut dyn ChFi3dCornerOps,
        surf_ops: &mut dyn ChFi3dSurfOps,
    ) -> Result<(), String> {
        // `UpdateTolesp();` (`cxx:223`).
        self.update_tolesp();

        // `if (myListStripe.IsEmpty()) throw Standard_Failure(...)`
        // (`cxx:225-228`).
        if self.my_list_stripe.is_empty() {
            return Err("There are no suitable edges for chamfer or fillet".to_string());
        }

        // `Reset(); myDS = new TopOpeBRepDS_HDataStructure();
        //  done = true; hasresult = false;` (`cxx:230-234`).
        self.reset();
        self.my_ds = crate::topopebrep::ds::HDataStructure::new();
        self.done = true;
        self.has_result = false;

        // `cxx:237-257`: filling of `myVDataMap`.
        let stripes: Vec<_> = self.my_list_stripe.clone();
        for stripe in &stripes {
            let Some(spine) = spine_of(stripe) else {
                continue;
            };
            let (first_status, last_status, first_vertex, last_vertex) = {
                let guard = spine
                    .read()
                    .map_err(|_| "Compute: poisoned spine".to_string())?;
                (
                    guard.first_status(),
                    guard.last_status(),
                    guard.first_vertex(),
                    guard.last_vertex(),
                )
            };
            // `if (Spine->FirstStatus() <= ChFiDS_BreakPoint) myVDataMap.Add(...)
            //  else if (FirstStatus() == ChFiDS_FreeBoundary) ExtentOneCorner(...)`
            if is_at_most_break_point(first_status) {
                if let Some(vtx) = &first_vertex {
                    self.my_v_data_map.add(vtx, stripe.clone());
                }
            } else if first_status == SpineState::FreeBoundary {
                if let Some(vtx) = &first_vertex {
                    ops.extent_one_corner(vtx, stripe.clone())?;
                }
            }
            if is_at_most_break_point(last_status) {
                if let Some(vtx) = &last_vertex {
                    self.my_v_data_map.add(vtx, stripe.clone());
                }
            } else if last_status == SpineState::FreeBoundary {
                if let Some(vtx) = &last_vertex {
                    ops.extent_one_corner(vtx, stripe.clone())?;
                }
            }
        }

        // `ExtentAnalyse();` (`cxx:259`).
        self.extent_analyse(ops)?;

        // `cxx:267-293`: construction of the band of fillet on every contour.
        // The `try` / `catch` around `PerformSetOfSurf` is the `Err` arm.
        for stripe in self.my_list_stripe.clone() {
            if let Some(spine) = spine_of(&stripe) {
                if let Ok(mut guard) = spine.write() {
                    guard.set_error_status(ErrorStatus::Ok);
                }
            }
            if let Err(_failure) = self.perform_set_of_surf(surf_ops, &stripe, false) {
                self.bad_stripes.push(stripe.clone());
                self.done = true;
                if let Some(spine) = spine_of(&stripe) {
                    if let Ok(mut guard) = spine.write() {
                        if guard.error_status() == ErrorStatus::Ok {
                            guard.set_error_status(ErrorStatus::Error);
                        }
                    }
                }
            }
            if !self.done {
                self.bad_stripes.push(stripe);
            }
            self.done = true;
        }
        self.done = self.bad_stripes.is_empty();

        // `cxx:302-332`: construct the fillets on every vertex and feed the DS.
        if self.done {
            // OCCT re-reads `myVDataMap.Extent()` at every iteration.
            let mut j = 1;
            while j <= self.my_v_data_map.extent() {
                let key = self.my_v_data_map.find_key(j).cloned();
                if self.perform_fillet_on_vertex(ops, j).is_err() {
                    if let Some(k) = &key {
                        self.bad_vertices.push(k.clone());
                    }
                    self.has_result = false;
                    self.done = true;
                }
                if !self.done {
                    if let Some(k) = &key {
                        self.bad_vertices.push(k.clone());
                    }
                }
                self.done = true;
                j += 1;
            }
            if !self.has_result {
                self.done = self.bad_vertices.is_empty();
            }
        }

        // `cxx:339-353`: `MapIndSo` = DS index of every solid, then of every
        // shell that is not inside a solid (`TopExp_Explorer(myShape,
        // TopAbs_SHELL, TopAbs_SOLID)`).
        let mut map_ind_so: Vec<i32> = Vec::new();
        for cursol in explorer_shapes(&self.my_shape, ShapeType::Solid) {
            let indcursol = self.my_ds.change_ds().add_shape(&cursol);
            if !map_ind_so.contains(&indcursol) {
                map_ind_so.push(indcursol);
            }
        }
        for cursh in explorer_shapes_avoiding(&self.my_shape, ShapeType::Shell, ShapeType::Solid) {
            let indcursh = self.my_ds.change_ds().add_shape(&cursh);
            if !map_ind_so.contains(&indcursh) {
                map_ind_so.push(indcursh);
            }
        }

        if self.done {
            // `cxx:356-396`: the OCC119 cross-check between the fillets, then
            // `ChFi3d_FilDS` on every contour.
            for i1 in 0..self.my_list_stripe.len() {
                let st = self.my_list_stripe[i1].clone();
                for i2 in (i1 + 1)..self.my_list_stripe.len() {
                    let a_check_stripe = self.my_list_stripe[i2].clone();
                    let intersection = {
                        let guard = st
                            .read()
                            .map_err(|_| "Compute: poisoned stripe".to_string())?;
                        let check_guard = a_check_stripe
                            .read()
                            .map_err(|_| "Compute: poisoned stripe".to_string())?;
                        ch_fi3d_stripe_edge_inter(&guard, &check_guard, self.tol2d)
                    };
                    if intersection.is_err() {
                        self.bad_stripes.push(st.clone());
                        self.has_result = false;
                        self.done = false;
                        break;
                    }
                }
                // `int solidindex = st->SolidIndex();
                //  ChFi3d_FilDS(solidindex, st, DStr, myRegul, tolapp3d, tol2d);`
                let solid_index = st
                    .read()
                    .map_err(|_| "Compute: poisoned stripe".to_string())?
                    .solid_index();
                self.fil_ds(solid_index, &st, self.tolapp3d, self.tol2d)?;
                if !self.done {
                    break;
                }
            }

            if self.done {
                // `CompleteDS(DStr, myShape);` (`cxx:406`).
                complete_ds(&mut self.my_ds, &self.my_shape);

                // `cxx:409-470`: update the tolerances on the vertices to the
                // max of the adjacent edges, and on the degenerated curves to
                // the max of the adjacent vertices.
                let nb_curves = self.my_ds.nb_curves();
                for ic in 1..=nb_curves {
                    // `TopOpeBRepDS_CurveExplorer cex(DStr)` (`cxx:409`) uses
                    // the default `FindOnlyKeep = true`, so `IsCurveKeep(I)` is
                    // the test: bound *and* `Curve(I).Keep()`
                    // (`TopOpeBRepDS_CurveExplorer.cxx:123-130`).
                    if !self.my_ds.ds().is_curve_bound(ic)
                        || !self.my_ds.ds().curve(ic).keep()
                    {
                        continue;
                    }
                    let degen = self.my_ds.ds().curve(ic).curve().is_none();
                    let mut tolc = if degen {
                        0.0
                    } else {
                        self.my_ds.ds().curve(ic).tolerance()
                    };
                    let interferences = self.my_ds.ds().curve_interferences(ic).to_vec();
                    let mut it = PointIterator::new(interferences);
                    while it.more() {
                        // `occ::down_cast<TopOpeBRepDS_CurvePointInterference>`
                        // (`cxx:423-428`): a null cast is skipped.
                        let geometry = {
                            let interference = it.current();
                            if interference.dynamic_type()
                                == "TopOpeBRepDS_CurvePointInterference"
                            {
                                Some((interference.geometry_type(), interference.geometry()))
                            } else {
                                None
                            }
                        };
                        let Some((gk, gi)) = geometry else {
                            it.next();
                            continue;
                        };
                        if gk == Kind::Vertex {
                            // `const TopoDS_Vertex& v = TopoDS::Vertex(myDS->Shape(gi));`
                            let vertex = self
                                .my_ds
                                .shape(gi, true)
                                .and_then(Vertex::wrap);
                            let Some(v) = vertex else {
                                it.next();
                                continue;
                            };
                            let mut tolv = BRepTool::vertex_tolerance(&v);
                            if tolv > 0.0001 {
                                tolv += 0.0003;
                                if tolc < tolv {
                                    tolc = tolv + 0.00001;
                                }
                            }
                            if degen && tolc < tolv {
                                tolc = tolv;
                            } else if tolc > tolv {
                                update_vertex_tolerance(&v, tolc);
                            }
                        } else if gk == Kind::Point {
                            // `TopOpeBRepDS_Point& p = DStr.ChangePoint(gi);`
                            if let Some(p) = self.my_ds.change_ds().change_point(gi) {
                                let tolp = p.tolerance();
                                if degen && tolc < tolp {
                                    tolc = tolp;
                                } else if tolc > tolp {
                                    p.set_tolerance(tolc);
                                }
                            }
                        }
                        it.next();
                    }
                    if degen {
                        self.my_ds.change_curve(ic).set_tolerance(tolc);
                    }
                }

                // `myCoup->Perform(myDS);` (`cxx:471`). The OCCT builder keeps
                // the handle it is given, so the DS it enriches is the one
                // `DStr` refers to; the port hands over a copy and takes the
                // builder's copy back, which restores that aliasing (the loop
                // below walks `DStr.NbShapes()` *after* `Perform`).
                self.my_coup.perform(self.my_ds.clone());
                if let Some(hds) = self.my_coup.data_structure() {
                    self.my_ds = hds.clone();
                }

                // `cxx:472-478`: `MergeSolid` on every indexed solid / shell.
                for indsol in &map_ind_so {
                    let Some(curshape) = self.my_ds.shape(*indsol, true) else {
                        continue;
                    };
                    self.my_coup.merge_solid(&curshape, TopAbsState::In);
                }

                // `cxx:480-509`: the vertices of the new edges reach the
                // tolerance of their edge.
                let nb_shapes = self.my_ds.nb_shapes();
                for i in 1..=nb_shapes {
                    let Some(s) = self.my_ds.shape(i, true) else {
                        continue;
                    };
                    if s.shape_type() != ShapeType::Edge {
                        continue;
                    }
                    if !self.my_coup.is_split(&s, TopAbsState::In) {
                        continue;
                    }
                    for new_edge in self.my_coup.splits(&s, TopAbsState::In).to_vec() {
                        let Some(new_e) = Edge::wrap(new_edge) else {
                            continue;
                        };
                        let tole = BRepTool::edge_tolerance(&new_e);
                        for vshape in explorer_shapes(&new_e.0, ShapeType::Vertex) {
                            let Some(v) = Vertex::wrap(vshape) else {
                                continue;
                            };
                            let tolv = BRepTool::vertex_tolerance(&v);
                            if tole > tolv {
                                update_vertex_tolerance(&v, tole);
                            }
                        }
                    }
                }

                // `cxx:510-567`: the compound of the result, or the compound of
                // the partial result when a corner reported `hasresult`.
                if !self.has_result {
                    let mut members: Vec<TopoShape> = Vec::new();
                    for indsol in &map_ind_so {
                        let Some(curshape) = self.my_ds.shape(*indsol, true) else {
                            continue;
                        };
                        let merged: Vec<TopoShape> =
                            self.my_coup.merged(&curshape, TopAbsState::In).to_vec();
                        if merged.is_empty() {
                            members.push(curshape);
                            continue;
                        }
                        // If the old type of Shape is SHELL, the shell is
                        // placed instead of the solid (`cxx:524-541`).
                        for its in merged {
                            if curshape.shape_type() == ShapeType::Shell {
                                // `TopExp_Explorer expsh2(its.Value(), TopAbs_SHELL);
                                //  const TopoDS_Shape& cursh = expsh2.Current();`
                                // - OCCT reads `Current()` without `More()`.
                                let cursh = explorer_shapes(&its, ShapeType::Shell)
                                    .into_iter()
                                    .next()
                                    .ok_or_else(|| {
                                        "Compute: merged shell has no shell".to_string()
                                    })?;
                                members.push(cursh);
                            } else {
                                members.push(its);
                            }
                        }
                    }
                    let compound = make_compound(&members);
                    self.my_shape_result = Some(compound);
                } else {
                    self.done = false;
                    let mut members: Vec<TopoShape> = Vec::new();
                    for indsol in &map_ind_so {
                        let Some(curshape) = self.my_ds.shape(*indsol, true) else {
                            continue;
                        };
                        let merged: Vec<TopoShape> =
                            self.my_coup.merged(&curshape, TopAbsState::In).to_vec();
                        if merged.is_empty() {
                            members.push(curshape);
                        } else {
                            members.extend(merged);
                        }
                    }
                    self.bad_shape = Some(make_compound(&members));
                }

                // `SetRegul();` (`cxx:574`): the regularities are coded after
                // cutting.
                ops.set_regul()?;
            }
        }

        // `cxx:655-674`: inspect the new faces to provide sameparameter.
        if self.done {
            let same_par_tol = CONFUSION;
            let nb_surfaces = self.my_ds.nb_surfaces();
            for i_f in 1..=nb_surfaces {
                let Some(new_faces) = self.my_coup.new_faces(i_f) else {
                    continue;
                };
                let new_faces: Vec<TopoShape> = new_faces.to_vec();
                for a_f in new_faces {
                    // `BRepLib::SameParameter(aF, SameParTol, true);
                    //  ShapeFix::SameParameter(aF, false, SameParTol);`
                    same_parameter_shape(&a_f, same_par_tol, true);
                    shape_fix_same_parameter(&a_f, false, same_par_tol);
                }
            }
        }

        Ok(())
    }
}

/// `FirstStatus() <= ChFiDS_BreakPoint` (`cxx:241`, `:249`): the enum order of
/// `ChFiDS_State.hxx:21-29` is `OnSame`, `OnDiff`, `AllSame`, `BreakPoint`,
/// `FreeBoundary`, `Closed`, `Tangent`.
fn is_at_most_break_point(status: SpineState) -> bool {
    matches!(
        status,
        SpineState::OnSame | SpineState::OnDiff | SpineState::AllSame | SpineState::BreakPoint
    )
}
