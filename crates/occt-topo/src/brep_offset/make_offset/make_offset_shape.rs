//! `BRepOffset_MakeOffset::MakeOffsetShape`
//! (`BRepOffset_MakeOffset.cxx:837-1080`): the driver that runs every phase of
//! the general offset in the reference order.
//!
//! Phases (`:843-1077`):
//!
//! 1. `IsPlanar` (`:843`) then `SetFaces` / `SetFacesWithOffset` /
//!    `BuildFaceComp` (`:845-848`).
//! 2. `RemoveCorks` on `myShape` and on `myFaceComp` when `myFaces` is not
//!    empty (`:853-857`).
//! 3. `CheckInputData` (`:866-871`), then `EvalMax(myShape, myTol)` and the
//!    angular tolerance `TolAngle = 4 * asin(min(myTol / (|myOffset / 2| +
//!    Confusion), 1))` (`:880-885`).
//! 4. `myAnalyse.Perform(myFaceComp, TolAngle)` (`:891`) followed by the
//!    `ChFiDS_Mixed` rejection (`:892-905`).
//! 5. `UpdateFaceOffset` (`:916`) and `BuildOffsetByArc` / `BuildOffsetByInter`
//!    (`:919-928`).
//! 6. `Intersection3D` (`:942-947`), `Intersection2D` (`:956-961`),
//!    `MakeLoops` (`:964-970`) and `MakeFaces` (`:975-982`).
//! 7. `MakeMissingWalls` for the thickening mode (`:987-993`).
//! 8. `MakeShells` (`:997`) and the `myOffsetShape.IsNull()` early exit
//!    (`:1003-1007`).
//! 9. `SelectShells` (`:1013`), `RemoveInternalEdges` (`:1017-1020`),
//!    `EncodeRegularity` (`:1023`), `ReplaceRoots` (`:1027`) and `MakeSolid`
//!    (`:1032`).
//! 10. `UpdateTolerance` (`:1039-1050`) + `BRepLib::UpdateTolerances`
//!     (`:1051`), then `CorrectConicalFaces` (`:1055`).
//! 11. The thickening sewing arm (`:1058-1076`) and `myDone = true`
//!     (`:1077`).
//!
//! Port notes:
//!
//! * `TopAbs_State Side` (`:872-876`) only feeds `BRepOffset_Inter3d Inter(Side)`
//!   (`:943`); the port computes it inside `intersection3_d` (see
//!   `intersection.rs`), so the local is not spelled out here.
//! * `BRepOffset_Analyse::Perform` returns `void` in the reference; the port
//!   answers `Result` and the outcome is ignored, exactly as the reference does
//!   with the status it does not have.
//! * `Message_ProgressScope` / `analyzeProgress` (`:859-863`, and the
//!   `aPS.Next(...)` / `!aPS.More()` -> `BRepOffset_UserBreak` sites) are not
//!   ported; there is no progress machinery in the crate.
//! * `myOffsetShape.Closed(true)` before the rebuilding `MakeSolid`
//!   (`:1070`).

use occt_core::kernel::geomabs::JoinType;
use occt_core::precision::CONFUSION;

use crate::abs::ShapeType;
use crate::brep_lib_update_tolerances::update_tolerances;
use crate::brep_offset::{Error, Interval};
use crate::chfi_ds::TypeOfConcavity;
use crate::shape::Edge;
use crate::topo_tools::explorer_shapes;

use super::build_offset_by_arc::build_offset_by_arc;
use super::build_offset_by_inter::build_offset_by_inter;
use super::face_comp::{build_face_comp, replace_roots, set_faces, set_faces_with_offset};
use super::intersection::intersection3_d;
use super::is_planar::is_planar;
use super::make_offset_faces::update_face_offset;
use super::shape_prep::{eval_max, remove_corks};
use super::tolerance_update::update_tolerance;
use super::BRepOffsetMakeOffset;

impl BRepOffsetMakeOffset {
    /// `MakeOffsetShape(theRange)` (`cxx:837-1080`).
    pub fn make_offset_shape(&mut self) {
        // `cxx:839`.
        self.my_done = false;

        // `cxx:843-848`.
        self.my_is_planar = is_planar(self);
        set_faces(self);
        set_faces_with_offset(self);
        build_face_comp(self);

        // `cxx:853-857`: `myShape` without the caps.
        if self.my_faces.extent() != 0 {
            if let Some(shape) = self.my_shape.as_mut() {
                remove_corks(shape, &mut self.my_original_faces);
            }
            remove_corks(&mut self.my_face_comp.0, &mut self.my_faces);
        }

        // `cxx:866-871`.
        if !self.check_input_data() || self.my_error != Error::NoError {
            return;
        }
        self.my_error = Error::NoError;

        // `cxx:880`: `EvalMax(myShape, myTol)`.
        if let Some(shape) = self.my_shape.clone() {
            eval_max(&shape, &mut self.my_tol);
        }

        // `cxx:883-885`: the angular tolerance of the preanalysis.
        let tol_angle_coeff =
            (self.my_tol / ((self.my_offset * 0.5).abs() + CONFUSION)).min(1.0);
        let tol_angle = 4.0 * tol_angle_coeff.asin();

        // `cxx:886-890`.
        if self.my_join == JoinType::Intersection && self.my_inter && self.my_is_planar {
            self.my_analyse.set_offset_value(self.my_offset);
            self.my_analyse
                .set_face_offset_map(self.my_face_offset.to_hash_map());
        }

        // `cxx:891`: `myAnalyse.Perform(myFaceComp, TolAngle)`.
        let face_comp = self.my_face_comp.0.clone();
        let _ = self.my_analyse.perform(&face_comp, tol_angle);

        // `cxx:892-905`: a partially C0 / tangent edge is rejected.
        for e_shape in explorer_shapes(&face_comp, ShapeType::Edge) {
            let Some(e) = Edge::wrap(e_shape) else {
                continue;
            };
            let mixed = last_interval_kind(self.my_analyse.type_of(&e)) == Some(TypeOfConcavity::Mixed);
            if mixed {
                self.my_error = Error::MixedConnectivity;
                return;
            }
        }

        // `cxx:916`: update of the face -> offset map.
        update_face_offset(self);

        // `cxx:919-928`.
        if self.my_join == JoinType::Arc {
            build_offset_by_arc(self);
        } else if self.my_join == JoinType::Intersection {
            build_offset_by_inter(self);
        }
        if self.my_error != Error::NoError {
            return;
        }

        // `cxx:942-947`: intersection 3D.
        let (modif, new_edges) = intersection3_d(self);
        if self.my_error != Error::NoError {
            return;
        }

        // `cxx:954-961`: intersection 2D.
        if modif.extent() != 0 {
            self.intersection2_d(&modif, &new_edges);
            if self.my_error != Error::NoError {
                return;
            }
        }

        // `cxx:964-970`.
        self.make_loops(&modif);
        if self.my_error != Error::NoError {
            return;
        }

        // `cxx:975-982`.
        if modif.extent() != 0 {
            self.make_faces();
            if self.my_error != Error::NoError {
                return;
            }
        }

        // `cxx:987-993`.
        if self.my_thickening {
            self.make_missing_walls();
            if self.my_error != Error::NoError {
                return;
            }
        }

        // `cxx:997`.
        self.make_shells();
        if self.my_error != Error::NoError {
            return;
        }

        // `cxx:1003-1007`: not done.
        if self.my_offset_shape.is_none() {
            self.my_done = false;
            return;
        }

        // `cxx:1013`: unwinding 3D.
        self.select_shells();

        // `cxx:1017-1020`: the `TopAbs_INTERNAL` edges.
        if self.my_remove_int_edges {
            self.remove_internal_edges();
        }

        // `cxx:1023`.
        self.encode_regularity();

        // `cxx:1027`.
        replace_roots(self);

        // `cxx:1032`.
        self.make_solid();
        if self.my_error != Error::NoError {
            return;
        }

        // `cxx:1039-1051`: edge and vertex tolerances.
        if let Some(offset_shape) = self.my_offset_shape.clone() {
            if self.my_thickening {
                let init = self.my_shape.clone();
                update_tolerance(&offset_shape, &self.my_faces, init.as_ref());
            } else {
                update_tolerance(&offset_shape, &self.my_faces, None);
            }
            // `cxx:1051`: `BRepLib::UpdateTolerances(myOffsetShape)`, whose
            // `theVerifyFaceTolerance` default is `false`
            // (`BRepLib.hxx:398-400`).
            let _ = update_tolerances(&offset_shape, false);
        }

        // `cxx:1055`.
        self.correct_conical_faces();

        // `cxx:1058-1076`: the result solid is sewn again in thickening mode.
        if self.my_thickening && self.my_is_perform_sewing {
            let Some(offset_shape) = self.my_offset_shape.clone() else {
                return;
            };
            let mut a_sew = crate::brep_builder_api_sewing::Sewing::new(
                self.my_tol,
                true,
                true,
                true,
                false,
            );
            a_sew.add(&offset_shape);
            if a_sew.perform().is_err() {
                self.my_error = Error::UnknownError;
                return;
            }
            let Some(sewed) = a_sew.sewed_shape().cloned() else {
                self.my_error = Error::UnknownError;
                return;
            };
            // `cxx:1070`.
            sewed.set_closed(true);
            self.my_offset_shape = Some(sewed);
            self.make_solid();
            if self.my_error != Error::NoError {
                return;
            }
        }

        // `cxx:1077`.
        self.my_done = true;
    }
}

/// `aLI.Last().Type()` (`cxx:900`) without the `IsEmpty` guard: the answer is
/// `None` for an empty list.
fn last_interval_kind(li: &[Interval]) -> Option<TypeOfConcavity> {
    li.last().map(|i| i.kind())
}
