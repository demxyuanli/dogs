//! `BRepFill_PipeShell::Build` / `MakeSolid` and the statics they need
//! (`BRepFill_PipeShell.cxx:152-220`, `:706-914`, `:1618-1681`).

use occt_core::gp::GpTrsf;
use occt_core::kernel::geomabs::Shape;
use occt_geom::geom_fill::sweep::ApproxStyle;
use occt_geom::geom_fill::{GeomFillPipeError, GeomFillTrihedron};

use crate::abs::ShapeType;
use crate::brep_class3d::SolidClassifier;
use crate::brep_fill_pipe_shell::PipeShell;
use crate::brep_fill_section_law::SectionLaw;
use crate::brep_fill_sweep::BrepFillSweep;
use crate::brep_tool::BRepTool;
use crate::fclass2d::FaceState;
use crate::shape::{Edge, Face, Shell, Solid, TopoShape, Wire};
use crate::topo_tools::{
    explorer_shapes, is_same, map_shapes_and_ancestors, shape_map_key, ShapeMapKey,
};
use crate::builder::TopoBuilder;

/// `Precision::Confusion()`.
const CONFUSION: f64 = 1.0e-7;
/// `Precision::Angular()`.
const ANGULAR: f64 = 1.0e-12;

impl PipeShell {
    /// `Build()` (`cxx:706-846`): the shell, the boundaries and the history.
    pub fn build(&mut self) -> Result<bool, String> {
        // 1) Preparation.
        self.prepare()?;

        if self.get_status() != GeomFillPipeError::PipeOk {
            self.set_shape(Self::empty_shell()); // Nullify
            return Ok(false);
        }

        let section = self
            .section()
            .cloned()
            .ok_or_else(|| "BRepFill_PipeShell::Build: no section law".to_string())?;
        let location = self
            .location()
            .cloned()
            .ok_or_else(|| "BRepFill_PipeShell::Build: no location law".to_string())?;

        // 2) `myFirst` and `myLast` (`cxx:723-757`).
        let (first_s, last_s) = section
            .base().law(1)
            .ok_or_else(|| "BRepFill_PipeShell::Build: no section law 1".to_string())?
            .get_domain();
        let my_first = Self::located_section(&section, &location, first_s, 0.0);
        let my_last: Option<TopoShape>;

        if section.base().is_v_closed() && location.base().is_closed() {
            if location.base().is_g1(0, CONFUSION, ANGULAR) >= 0 {
                my_last = my_first.clone();
            } else {
                my_last = None;
            }
            let my_first = if my_last.is_some() { my_first } else { None };
            self.set_first(my_first);
            self.set_last(my_last.clone());
        } else {
            let (_first, length) = location.base().curvilinear_bounds(location.base().nb_law())?;
            let mut last_shape = Self::located_section(&section, &location, last_s, length);
            // `cxx:747-754` (occ332): `myLast` and `myFirst` must not share one
            // TShape, tolerances of shapes built on them may differ.
            if let (Some(a), Some(b)) = (my_first.as_ref(), last_shape.as_ref()) {
                if a.same_tshape(b) {
                    let copy = crate::shape_ops::transformed_copy(b, &GpTrsf::identity())?;
                    last_shape = Some(copy);
                }
            }
            self.set_first(my_first);
            self.set_last(last_shape);
        }
        let my_first = self.first_ref().cloned();
        let my_last = self.last_ref().cloned();

        // 3) Construction (`cxx:759-788`).
        let mut mk_sw = BrepFillSweep::new(section.clone(), location.clone(), true);
        mk_sw.set_tolerance(self.tol_3d(), self.bound_tol(), self.tol_2d(), self.tol_angular());
        mk_sw.set_angular_control(self.ang_min(), self.ang_max());
        mk_sw.set_force_approx_c1(self.force_approx_c1());
        mk_sw.set_bounds(
            my_first.clone().and_then(Wire::wrap),
            my_last.clone().and_then(Wire::wrap),
        );

        let the_continuity = if self.trihedron() == GeomFillTrihedron::DiscreteTrihedron {
            Shape::C0
        } else {
            Shape::C2
        };
        let mut dummy: std::collections::HashSet<ShapeMapKey> = std::collections::HashSet::new();
        let mut dummy2: std::collections::HashMap<ShapeMapKey, crate::brep_fill_sweep::helpers::ShapeTable> =
            std::collections::HashMap::new();
        let mut dummy3: std::collections::HashMap<ShapeMapKey, crate::brep_fill_sweep::helpers::ShapeTable> =
            std::collections::HashMap::new();
        let _ = mk_sw.build(
            &mut dummy,
            &mut dummy2,
            &mut dummy3,
            self.transition(),
            the_continuity,
            ApproxStyle::Location,
            self.max_degree(),
            self.max_segments(),
        );

        self.set_status(location.base().get_status());
        let ok = mk_sw.is_done() && self.get_status() == GeomFillPipeError::PipeOk;

        if ok {
            self.set_shape(mk_sw.shape().cloned());
            self.set_error_on_surf(mk_sw.error_on_surface());

            let mut a_bottom_wire = self.first_ref().cloned();
            let mut a_top_wire = self.last_ref().cloned();
            if build_boundaries(&mk_sw, section.as_ref(), &mut a_bottom_wire, &mut a_top_wire) {
                self.set_first(a_bottom_wire);
                self.set_last(a_top_wire);
            }

            if section.base().is_u_closed() {
                // `cxx:807-825`: a seam-only shell is closed.
                let degen_first = self.first_ref().map(all_edges_degenerated).unwrap_or(true);
                let degen_last = self.last_ref().map(all_edges_degenerated).unwrap_or(true);
                if degen_first && degen_last {
                    if let Some(shape) = self.shape() {
                        shape.set_closed(true);
                    }
                }
            }

            if self.is_build_history() {
                self.build_history(&mk_sw)?;
            }
        } else {
            self.set_shape(Self::empty_shell()); // Nullify
            if self.get_status() == GeomFillPipeError::PipeOk {
                self.set_status(GeomFillPipeError::PipeNotOk);
            }
        }
        Ok(ok)
    }

    /// `Build()` step 2, one boundary: `Section.D0(param)` moved onto the frame
    /// the location law has at `abscissa` (`cxx:724-745`).
    fn located_section(
        section: &std::sync::Arc<dyn SectionLaw>,
        location: &std::sync::Arc<dyn crate::brep_fill_location_law::LocationLaw>,
        param: f64,
        abscissa: f64,
    ) -> Option<TopoShape> {
        let shape = section.d0(param)?;
        location.base().d0(abscissa, &shape).ok().flatten()
    }

    /// `MakeSolid()` (`cxx:849-914`): cap `myFirst` / `myLast` and build the
    /// solid.
    pub fn make_solid(&mut self) -> Result<bool, String> {
        if self.shape().is_none() {
            return Err("StdFail_NotDone: PipeShell is not built".to_string());
        }
        let mut b = self.shape().map(|s| s.is_closed()).unwrap_or(false);

        if !b {
            if self.first_ref().is_some() && self.last_ref().is_some() {
                b = self.first_ref().map(|s| s.is_closed()).unwrap_or(false)
                    && self.last_ref().map(|s| s.is_closed()).unwrap_or(false);
            }
            if b {
                // It is necessary to block the extremities.
                b = perform_plan(self.first_mut())?;
                if b {
                    b = perform_plan(self.last_mut())?;
                    if b {
                        let shell = self.shape().cloned();
                        if let (Some(first), Some(shell)) = (self.first_ref().cloned(), shell.as_ref())
                        {
                            if !is_same_oriented(&first, shell) {
                                let mut reversed = first;
                                reversed.reverse();
                                self.set_first(Some(reversed));
                            }
                        }
                        if let (Some(last), Some(shell)) = (self.last_ref().cloned(), shell.as_ref())
                        {
                            if !is_same_oriented(&last, shell) {
                                let mut reversed = last;
                                reversed.reverse();
                                self.set_last(Some(reversed));
                            }
                        }

                        let builder = TopoBuilder::new();
                        if let Some(first) = self.first_ref().cloned() {
                            if let Some(face) = Face::wrap(first) {
                                if let Some(shell) = self.shape_mut() {
                                    builder.add(shell, &face.0);
                                }
                            }
                        }
                        if let Some(last) = self.last_ref().cloned() {
                            if let Some(face) = Face::wrap(last) {
                                if let Some(shell) = self.shape_mut() {
                                    builder.add(shell, &face.0);
                                }
                            }
                        }
                        if let Some(shape) = self.shape() {
                            shape.set_closed(true);
                        }
                    }
                }
            }
        }

        if b {
            let builder = TopoBuilder::new();
            let Some(shell) = self.shape().cloned().and_then(Shell::wrap) else {
                return Ok(false);
            };
            let mut solid = Solid::new();
            builder.add_shell(&mut solid, &shell);

            let mut classifier = SolidClassifier::new();
            classifier.load(solid.0.clone());
            classifier.perform_infinite_point(CONFUSION);
            if classifier.state() == FaceState::In {
                let mut reversed_solid = Solid::new();
                if let Some(shape) = self.shape_mut() {
                    shape.reverse();
                }
                builder.add_shell(&mut reversed_solid, &shell);
                solid = reversed_solid;
            }
            self.set_shape(Some(solid.0));
            if let Some(shape) = self.shape() {
                shape.set_closed(true);
            }
        }
        Ok(b)
    }
}

/// `PerformPlan(S)` (`cxx:152-183`): turn a degenerate boundary into a null
/// shape and a wire into a planar face.
pub fn perform_plan(s: &mut Option<TopoShape>) -> Result<bool, String> {
    let mut is_degen = true;
    if let Some(shape) = s.as_ref() {
        for an_edge in explorer_shapes(shape, ShapeType::Edge) {
            if let Some(edge) = Edge::wrap(an_edge) {
                if !BRepTool::is_degenerated(&edge) {
                    is_degen = false;
                }
            }
        }
    }
    if is_degen {
        *s = None;
        return Ok(true);
    }

    let Some(shape) = s.as_ref() else {
        return Ok(false);
    };
    let Some(wire) = Wire::wrap(shape.clone()) else {
        return Ok(false);
    };
    match crate::brep_lib_make_face::make_face_from_wire(&wire, true) {
        Some(face) => {
            *s = Some(face.0);
            Ok(true)
        }
        None => Ok(false),
    }
}

/// `IsSameOriented(theFace, theShell)` (`cxx:187-220`): whether the boundary
/// face already points the same way as the shell it will be glued to.
pub fn is_same_oriented(the_face: &TopoShape, the_shell: &TopoShape) -> bool {
    let an_edge_face_map = map_shapes_and_ancestors(the_shell, ShapeType::Edge, ShapeType::Face);

    for a_current_edge in explorer_shapes(the_face, ShapeType::Edge) {
        let Some(ancestors) = an_edge_face_map.get(&shape_map_key(&a_current_edge)) else {
            continue;
        };
        let Some(an_adj_face) = ancestors.first() else {
            continue;
        };
        for an_adj_edge in explorer_shapes(an_adj_face, ShapeType::Edge) {
            // IsSame compares TShape and Location only; orientation may differ.
            if is_same(&an_adj_edge, &a_current_edge) {
                return a_current_edge.orientation() != an_adj_edge.orientation();
            }
        }
    }

    // No shared edge found. Orientation cannot be determined reliably and this
    // should not happen for valid input. Preserve the current face orientation
    // instead of reporting an opposite orientation and forcing a reversal.
    true
}

/// `BuildBoundaries(theSweep, theSection, theBottom, theTop)`
/// (`cxx:1618-1681`): the first and last section wires of the sweep.
pub fn build_boundaries(
    the_sweep: &BrepFillSweep,
    the_section: &dyn SectionLaw,
    the_bottom: &mut Option<TopoShape>,
    the_top: &mut Option<TopoShape>,
) -> bool {
    let builder = TopoBuilder::new();
    let mut a_bottom_wire = Wire::new();
    let mut a_top_wire = Wire::new();
    let a_v_edges = the_sweep.sections();
    let mut bfoundbottom = false;
    let mut bfoundtop = false;
    let mut b_all_same = true;

    for i in 1..=the_section.base().nb_law() {
        let a_bottom_edge = a_v_edges.value(i, 1).cloned();
        if let Some(edge) = a_bottom_edge.as_ref() {
            if edge.is_edge() {
                builder.add(&mut a_bottom_wire.0, edge);
                bfoundbottom = true;
            }
        }
        let a_top_edge = a_v_edges.value(i, a_v_edges.nb_cols()).cloned();
        if let Some(edge) = a_top_edge.as_ref() {
            if edge.is_edge() {
                builder.add(&mut a_top_wire.0, edge);
                bfoundtop = true;
            }
        }

        if let (Some(bottom), Some(top)) = (a_bottom_edge.as_ref(), a_top_edge.as_ref()) {
            if !is_same(bottom, top) {
                b_all_same = false;
            }
        }
    }

    if the_section.base().is_u_closed() {
        a_bottom_wire.set_closed(true);
        a_top_wire.set_closed(true);
    }

    if bfoundbottom {
        *the_bottom = Some(a_bottom_wire.0);
    }
    if bfoundtop {
        *the_top = Some(a_top_wire.0);
    }
    if b_all_same && bfoundbottom && bfoundtop {
        *the_top = the_bottom.clone();
    }

    bfoundbottom || bfoundtop
}

/// `cxx:807-821`: every edge of the boundary is degenerated.
fn all_edges_degenerated(shape: &TopoShape) -> bool {
    let mut degen = true;
    for an_edge in explorer_shapes(shape, ShapeType::Edge) {
        let is_degen = Edge::wrap(an_edge)
            .map(|e| BRepTool::is_degenerated(&e))
            .unwrap_or(false);
        degen = degen && is_degen;
    }
    degen
}
