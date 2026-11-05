//! `PerformSweep` (`BRepFill_AdvancedEvolved.cxx:507-526`) and `GetLids`
//! (`cxx:530-700`).
//!
//! `PerformSweep` sweeps the profile along the spine with a
//! [`PipeShell`](crate::brep_fill_pipe_shell::PipeShell) and leaves
//! `myErrorStatus = NoLids` so that `GetLids` runs. `GetLids` collects the
//! free edges of the pipe shell that are parallel to the (planar) spine plane,
//! fuses them, chains them into wires and turns the wires into the top and
//! bottom lids, orienting each lid against the pipe shell's vertices.
//!
//! Translation notes:
//! - `BRep_Tool::IsClosed` is [`BRepTool::is_closed`], whose port already
//!   covers the `TopoDS_Wire` case (`BRep_Tool.cxx:322-345`).
//! - `BRepAdaptor_Curve::Is3DCurve()` (`cxx:591`) tests whether the edge has a
//!   `BRep_Curve3D`; the port's equivalent is
//!   [`BRepTool::edge_curve`] being `Some` (`BRep_Tool::Curve(E)`).
//! - `BOPAlgo_Tools::EdgesToWires` / `WiresToFaces` answer their result
//!   instead of filling `aCompW` / `aCompF`
//!   ([`edges_to_wires`](crate::bopalgo_tools_wires::edges_to_wires),
//!   [`wires_to_faces`](crate::bopalgo_tools_wires::wires_to_faces)).
//! - `aDPMax` (`cxx:577`, `:606-609`) is written by the reference and never
//!   read (`grep aDPMax` finds only those three lines); it is kept so the
//!   control flow stays the reference's.

use occt_core::gp::GpVec;
use occt_core::precision::{CONFUSION, SQUARE_CONFUSION};

use crate::abs::{Orientation, ShapeType};
use crate::bbox_from_geometry::shape_bbox;
use crate::bopalgo_tools_wires::{edges_to_wires, wires_to_faces};
use crate::brep_fill_pipe_shell::PipeShell;
use crate::brep_fill_sweep::TransitionStyle;
use crate::brep_lib_find_surface::FindSurface;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, TopoShape, Vertex};
use crate::topo_tools::map_shapes;
use crate::topopebrep::tool::adaptor::CurveAdaptor;

use super::maps::ordered_map_shapes_and_ancestors;
use super::{
    AdvancedEvolved, STATUS_EMPTY, STATUS_NO_LIDS, STATUS_NOT_PLANAR_SPINE, STATUS_NOT_SOLID,
    STATUS_SWEEP_ERROR,
};

/// `aPipeLinearTolerance` (`cxx:55`).
const PIPE_LINEAR_TOLERANCE: f64 = 1.0e-4;
/// `aPipeAngularTolerance` (`cxx:56`).
const PIPE_ANGULAR_TOLERANCE: f64 = 1.0e-2;
/// `SetTransition`'s default `Angmin` (`BRepFill_PipeShell.hxx:177`).
const TRANSITION_ANG_MIN: f64 = 1.0e-2;
/// `SetTransition`'s default `Angmax` (`BRepFill_PipeShell.hxx:178`).
const TRANSITION_ANG_MAX: f64 = 6.0;

impl AdvancedEvolved {
    /// `PerformSweep()` (`cxx:507-526`).
    pub(super) fn perform_sweep(&mut self) {
        if self.status != STATUS_EMPTY {
            return;
        }

        self.status = STATUS_SWEEP_ERROR;

        // `cxx:512-513`: `new BRepFill_PipeShell(mySpine)`. `mySpine` is
        // always set by `GetSpineAndProfile` (`cxx:108`); a missing spine or
        // profile would make the reference's `Build()` fail, so the port
        // leaves `SweepError` in place instead of sweeping a null shape.
        let (Some(spine), Some(profile)) = (self.spine.clone(), self.profile.clone()) else {
            return;
        };

        // `cxx:513`: `SetTolerance(aPipeLinearTolerance,
        // aPipeLinearTolerance, aPipeAngularTolerance)`.
        let mut a_pipe = PipeShell::new(&spine);
        a_pipe.set_tolerance(
            PIPE_LINEAR_TOLERANCE,
            PIPE_LINEAR_TOLERANCE,
            PIPE_ANGULAR_TOLERANCE,
        );
        // `cxx:514`: `SetTransition(BRepFill_Round)` (default angles).
        a_pipe.set_transition(
            TransitionStyle::Round,
            TRANSITION_ANG_MIN,
            TRANSITION_ANG_MAX,
        );
        // `cxx:515`: `Add(myProfile, false, false)`; the reference ignores the
        // failure status and lets `Build()` report it.
        let _ = a_pipe.add(&profile.0, false, false);

        // `cxx:517-521`.
        if let Ok(true) = a_pipe.build() {
            self.status = STATUS_NO_LIDS;
            self.pipe_shell = a_pipe.shape().cloned();
        }
    }

    /// `GetLids()` (`cxx:530-700`).
    pub(super) fn get_lids(&mut self) {
        // `cxx:532-535`.
        let Some(pipe_shell) = self.pipe_shell.clone() else {
            return;
        };

        // `cxx:537-542`: `BRep_Tool::IsClosed(myProfile)`.
        let Some(profile) = self.profile.clone() else {
            return;
        };
        if BRepTool::is_closed(&profile.0) {
            self.status = STATUS_NOT_SOLID;
            return;
        }

        // `cxx:544`.
        self.status = STATUS_NO_LIDS;

        // `cxx:546-548`: `BRepLib_FindSurface(aFS(mySpine, -1.0, true))` and
        // `down_cast<Geom_Plane>(aFS.Surface())`.
        let Some(spine) = self.spine.clone() else {
            return;
        };
        let a_fs = FindSurface::new(&spine.0, -1.0, true, false);
        let a_surf = match a_fs.surface() {
            Some(s) => s,
            None => {
                self.status = STATUS_NOT_PLANAR_SPINE;
                return;
            }
        };
        let Some(pln) = a_surf.gp_pln() else {
            self.status = STATUS_NOT_PLANAR_SPINE;
            return;
        };

        // `cxx:553-560`.
        let mut a_tol = a_fs.tolerance().max(a_fs.tolerance_reached());
        a_tol += self.fuzzy_value;
        let a_prof_box = shape_bbox(&profile.0);
        let a_sq_diag = a_prof_box.square_extent();
        // `cxx:560-562`.
        let a_sq_angular_tol = a_tol * a_tol / a_sq_diag;
        let a_normal = GpVec::from_xyz(pln.axis().direction().xyz());

        // `cxx:568-571`: `TopExp::MapShapesAndAncestors(myPipeShell, EDGE,
        // FACE, aMapEF)`.
        let a_map_ef =
            ordered_map_shapes_and_ancestors(&pipe_shell, ShapeType::Edge, ShapeType::Face);

        // `cxx:573-618`.
        let mut a_le: Vec<TopoShape> = Vec::new();
        let mut a_dp_max = 0.0_f64;
        for entry in &a_map_ef {
            // `cxx:583-586`.
            if entry.ancestors.len() != 1 {
                continue;
            }

            let an_e = Edge(entry.key.clone());

            // `cxx:590-596`: `BRepAdaptor_Curve(anE)`, and "we are not
            // interested in degenerated edges".
            let an_ac = CurveAdaptor::of_edge(&an_e);
            if BRepTool::edge_curve(&an_e).is_none() {
                continue;
            }

            // `cxx:598-601`. `D1` has no failure status in the reference; a
            // missing derivative leaves the edge out of the planar set.
            let u_mid = 0.5 * (an_ac.first_parameter() + an_ac.last_parameter());
            let Some((_, a_tan)) = an_ac.d1(u_mid) else {
                continue;
            };

            let a_sq_modulus = a_tan.square_magnitude();
            if a_sq_modulus < CONFUSION {
                continue;
            }

            // `cxx:605-609`: `aDP = aTan.XYZ().Dot(aNormal.XYZ())`.
            let a_dp = a_tan.dot(&a_normal);
            if a_dp.abs() > a_dp_max {
                a_dp_max = a_dp.abs();
            }
            // `cxx:610-614`: only planar edges are considered.
            if a_dp * a_dp > a_sq_modulus * a_sq_angular_tol {
                continue;
            }

            a_le.push(an_e.0);
        }
        // `aDPMax` is never read by the reference (`cxx:577`, `:606-609`).
        let _ = a_dp_max;

        // `cxx:619-623`.
        if a_le.is_empty() {
            self.status = STATUS_NOT_PLANAR_SPINE;
            return;
        }

        // `cxx:625-631`: split interfered edges.
        let Some(a_free_edges) = self.perform_boolean(&a_le) else {
            self.status = STATUS_NOT_PLANAR_SPINE;
            return;
        };

        // `cxx:633-641`: "collect all free edges to wires and create planar
        // top and bottom lids from these wires".
        let builder = TopoBuilder::new();
        let an_ang_tol = a_sq_angular_tol.sqrt();
        let (_n_wires, a_comp_w) = edges_to_wires(&a_free_edges, true, an_ang_tol);
        let (_ok, a_comp_f) = wires_to_faces(&a_comp_w, an_ang_tol);
        let mut a_top_bottom = builder.make_compound_of(&[]).0;

        // `cxx:643-695`: check orientation.
        let a_map_v = map_shapes(&pipe_shell, &[ShapeType::Vertex])
            .remove(&ShapeType::Vertex)
            .unwrap_or_default();
        for face_shape in crate::topo_tools::explorer_shapes(&a_comp_f, ShapeType::Face) {
            let a_f = crate::shape::Face(face_shape);
            // `cxx:652-654`: the reference dereferences the `Geom_Plane`
            // without a null check; `WiresToFaces` only answers planar faces.
            let Some(surf) = BRepTool::face_surface(&a_f) else {
                continue;
            };
            let Some(a_pln) = surf.gp_pln() else {
                continue;
            };
            let a_norm = GpVec::from_xyz(a_pln.axis().direction().xyz());
            let a_loc_p = a_pln.location();

            let mut is_found = false;

            for v in &a_map_v {
                // `cxx:660-662`.
                let a_p = BRepTool::vertex_point(&Vertex(v.clone()));
                let a_delta = GpVec::new(
                    a_p.x() - a_loc_p.x(),
                    a_p.y() - a_loc_p.y(),
                    a_p.z() - a_loc_p.z(),
                );
                let a_sq_d = a_delta.square_magnitude();
                if a_sq_d < SQUARE_CONFUSION {
                    continue;
                }

                // `cxx:672-678`: "aP is in the plane".
                let a_dp = a_delta.dot(&a_norm);
                if a_dp * a_dp < a_sq_d * SQUARE_CONFUSION {
                    continue;
                }

                // `cxx:680-690`.
                if a_dp > 0.0 {
                    builder.add(&mut a_top_bottom, &a_f.0.oriented(Orientation::Reversed));
                } else {
                    builder.add(&mut a_top_bottom, &a_f.0);
                }

                is_found = true;
                break;
            }

            // `cxx:692-695`.
            if !is_found {
                builder.add(&mut a_top_bottom, &a_f.0);
            }
        }

        // `cxx:698-699`.
        self.top_bottom = Some(a_top_bottom);
        self.status = STATUS_NOT_SOLID;
    }
}