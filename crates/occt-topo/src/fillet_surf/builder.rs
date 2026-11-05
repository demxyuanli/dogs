//! `FilletSurf_Builder` — the public API over a constant-radius fillet.
//!
//! Source: `ModelingAlgorithms/TKFillet/FilletSurf/FilletSurf_Builder.cxx`
//! (OCCT 8.0.0), `hxx:40-153`.
//!
//! For one fillet the API answers the list of NUBS surfaces, and for each
//! surface the two support faces, the 3d curve and the pcurves on each face,
//! the two pcurves on the fillet, the status of the start and end sections,
//! and the first / last parameter on the edge.
//!
//! Divergences that are not heuristics:
//! - the `StdFail_NotDone` of `NbSurface`, `FirstParameter`, `LastParameter`,
//!   `StartSectionStatus`, `EndSectionStatus` and `NbSection` is `Err`, the
//!   `Standard_OutOfRange` of the indexed accessors is `Err` as well,
//! - the `Standard_Failure` `myIntBuild.Add` may raise out of the constructor
//!   (`cxx:41`) cannot escape a Rust constructor, so it is reported as
//!   [`StatusDone::IsNotOk`] with [`ErrorTypeStatus::PbFilletCompute`]; OCCT
//!   leaves `myisdone` / `myerrorstatus` untouched on that path because it
//!   never returns,
//! - the `Section(...)` overload that writes through a reference argument
//!   (`cxx:381-386`) is `Standard_DEPRECATED` in OCCT 8.0.0 and is not
//!   translated; [`FilletSurfBuilder::section`] returns the handle by value.

use std::sync::Arc;

use occt_geom::{Curve as GeomCurve, GeomTrimmedCurve, Surface as GeomSurface};
use occt_geom2d::curve::Curve2d;

use crate::chfi3d::fillet_shape::FilletShape;
use crate::shape::{Face, TopoShape};

use super::internal_builder::FilletSurfInternalBuilder;
use super::status::{ErrorTypeStatus, StatusDone, StatusType};

/// `FilletSurf_Builder` (`FilletSurf_Builder.hxx:48-153`).
pub struct FilletSurfBuilder {
    /// `FilletSurf_InternalBuilder myIntBuild;` (`hxx:150`).
    int_build: FilletSurfInternalBuilder,
    /// `FilletSurf_StatusDone myisdone;` (`hxx:151`).
    is_done: StatusDone,
    /// `FilletSurf_ErrorTypeStatus myerrorstatus;` (`hxx:152`).
    error_status: ErrorTypeStatus,
}

impl FilletSurfBuilder {
    /// `FilletSurf_Builder(S, E, R, Ta, Tapp3d, Tapp2d)` (`cxx:32-67`).
    pub fn new(
        shape: &TopoShape,
        edges: &[TopoShape],
        r: f64,
        ta: f64,
        tapp3d: f64,
        tapp2d: f64,
    ) -> Self {
        // `myIntBuild(S, ChFi3d_Polynomial, Ta, Tapp3d, Tapp2d)` (`cxx:39`).
        let mut this = Self {
            int_build: FilletSurfInternalBuilder::new(
                shape,
                FilletShape::Polynomial,
                ta,
                tapp3d,
                tapp2d,
            ),
            // `myisdone = FilletSurf_IsOk; myerrorstatus = FilletSurf_EmptyList;`
            // (`cxx:40-41`).
            is_done: StatusDone::IsOk,
            error_status: ErrorTypeStatus::EmptyList,
        };
        // `int add = myIntBuild.Add(E, R); if (add != 0) { ... }`
        // (`cxx:42-66`).
        let add = match this.int_build.add(edges, r) {
            Ok(code) => code,
            Err(_) => {
                this.is_done = StatusDone::IsNotOk;
                this.error_status = ErrorTypeStatus::PbFilletCompute;
                return this;
            }
        };
        if add != 0 {
            this.is_done = StatusDone::IsNotOk;
            this.error_status = match add {
                1 => ErrorTypeStatus::EmptyList,
                2 => ErrorTypeStatus::EdgeNotG1,
                3 => ErrorTypeStatus::FacesNotG1,
                4 => ErrorTypeStatus::EdgeNotOnShape,
                5 => ErrorTypeStatus::NotSharpEdge,
                _ => ErrorTypeStatus::PbFilletCompute,
            };
        }
        this
    }

    /// `Perform()` (`cxx:74-93`).
    pub fn perform(&mut self) {
        if self.is_done != StatusDone::IsOk {
            return;
        }
        // The `PerformSetOfKGen` walk raises on a broken contour; OCCT lets the
        // exception out of `Perform`, the port keeps the "no result" status.
        if self.int_build.perform().is_err() {
            self.is_done = StatusDone::IsNotOk;
            self.error_status = ErrorTypeStatus::PbFilletCompute;
            return;
        }
        if self.int_build.done() {
            self.is_done = StatusDone::IsOk;
        } else if self.int_build.nb_surface() != 0 {
            self.is_done = StatusDone::IsPartial;
            self.error_status = ErrorTypeStatus::PbFilletCompute;
        } else {
            self.is_done = StatusDone::IsNotOk;
            self.error_status = ErrorTypeStatus::PbFilletCompute;
        }
    }

    /// `Simulate()` (`cxx:325-341`).
    pub fn simulate(&mut self) {
        if self.is_done != StatusDone::IsOk {
            return;
        }
        if self.int_build.simulate().is_err() {
            self.is_done = StatusDone::IsNotOk;
            self.error_status = ErrorTypeStatus::PbFilletCompute;
            return;
        }
        if self.int_build.done() {
            self.is_done = StatusDone::IsOk;
        } else {
            self.is_done = StatusDone::IsNotOk;
            self.error_status = ErrorTypeStatus::PbFilletCompute;
        }
    }

    /// `IsDone()` (`cxx:99-102`).
    pub fn is_done(&self) -> StatusDone {
        self.is_done
    }

    /// `StatusError()` (`cxx:108-111`).
    pub fn status_error(&self) -> ErrorTypeStatus {
        self.error_status
    }

    /// `NbSurface()` (`cxx:119-125`): `StdFail_NotDone` when
    /// [`StatusDone::IsNotOk`].
    pub fn nb_surface(&self) -> Result<usize, String> {
        if self.is_done == StatusDone::IsNotOk {
            return Err("FilletSurf_Builder::NbSurface".to_string());
        }
        Ok(self.int_build.nb_surface())
    }

    /// `SurfaceFillet(Index)` (`cxx:132-139`).
    pub fn surface_fillet(&self, index: usize) -> Result<Arc<dyn GeomSurface>, String> {
        self.check_index(index, "FilletSurf_Builder::SurfaceFillet")?;
        self.int_build
            .surface_fillet(index)
            .ok_or_else(|| format!("FilletSurf_Builder::SurfaceFillet({index}): no surface"))
    }

    /// `TolApp3d(Index)` (`cxx:145-152`).
    pub fn tol_app3d(&self, index: usize) -> Result<f64, String> {
        self.check_index(index, "FilletSurf_Builder::TolApp3d")?;
        self.int_build
            .tol_app3d(index)
            .ok_or_else(|| format!("FilletSurf_Builder::TolApp3d({index}): no surface"))
    }

    /// `SupportFace1(Index)` (`cxx:157-164`).
    pub fn support_face1(&self, index: usize) -> Result<Face, String> {
        self.check_index(index, "FilletSurf_Builder::SupportFace1")?;
        self.int_build
            .support_face1(index)
            .ok_or_else(|| format!("FilletSurf_Builder::SupportFace1({index}): no face"))
    }

    /// `SupportFace2(Index)` (`cxx:171-178`).
    pub fn support_face2(&self, index: usize) -> Result<Face, String> {
        self.check_index(index, "FilletSurf_Builder::SupportFace2")?;
        self.int_build
            .support_face2(index)
            .ok_or_else(|| format!("FilletSurf_Builder::SupportFace2({index}): no face"))
    }

    /// `CurveOnFace1(Index)` (`cxx:185-192`).
    pub fn curve_on_face1(&self, index: usize) -> Result<Arc<dyn GeomCurve>, String> {
        self.check_index(index, "FilletSurf_Builder::CurveOnFace1")?;
        self.int_build
            .curve_on_face1(index)
            .ok_or_else(|| format!("FilletSurf_Builder::CurveOnFace1({index}): no curve"))
    }

    /// `CurveOnFace2(Index)` (`cxx:197-204`).
    pub fn curve_on_face2(&self, index: usize) -> Result<Arc<dyn GeomCurve>, String> {
        self.check_index(index, "FilletSurf_Builder::CurveOnFace2")?;
        self.int_build
            .curve_on_face2(index)
            .ok_or_else(|| format!("FilletSurf_Builder::CurveOnFace2({index}): no curve"))
    }

    /// `PCurveOnFace1(Index)` (`cxx:210-217`).
    pub fn p_curve_on_face1(&self, index: usize) -> Result<Arc<dyn Curve2d>, String> {
        self.check_index(index, "FilletSurf_Builder::PCurveOnFace1")?;
        self.int_build
            .p_curve_on_face1(index)
            .ok_or_else(|| format!("FilletSurf_Builder::PCurveOnFace1({index}): no pcurve"))
    }

    /// `PCurve1OnFillet(Index)` (`cxx:224-231`).
    pub fn p_curve1_on_fillet(&self, index: usize) -> Result<Arc<dyn Curve2d>, String> {
        self.check_index(index, "FilletSurf_Builder::PCurve1OnFillet")?;
        self.int_build
            .p_curve1_on_fillet(index)
            .ok_or_else(|| format!("FilletSurf_Builder::PCurve1OnFillet({index}): no pcurve"))
    }

    /// `PCurveOnFace2(Index)` (`cxx:236-243`).
    pub fn p_curve_on_face2(&self, index: usize) -> Result<Arc<dyn Curve2d>, String> {
        self.check_index(index, "FilletSurf_Builder::PCurveOnFace2")?;
        self.int_build
            .p_curve_on_face2(index)
            .ok_or_else(|| format!("FilletSurf_Builder::PCurveOnFace2({index}): no pcurve"))
    }

    /// `PCurve2OnFillet(Index)` (`cxx:249-256`).
    pub fn p_curve2_on_fillet(&self, index: usize) -> Result<Arc<dyn Curve2d>, String> {
        self.check_index(index, "FilletSurf_Builder::PCurve2OnFillet")?;
        self.int_build
            .p_curve2_on_fillet(index)
            .ok_or_else(|| format!("FilletSurf_Builder::PCurve2OnFillet({index}): no pcurve"))
    }

    /// `FirstParameter()` (`cxx:262-269`).
    pub fn first_parameter(&self) -> Result<f64, String> {
        if self.is_done == StatusDone::IsNotOk {
            return Err("FilletSurf_Builder::FirstParameter".to_string());
        }
        Ok(self.int_build.first_parameter())
    }

    /// `LastParameter()` (`cxx:274-281`).
    pub fn last_parameter(&self) -> Result<f64, String> {
        if self.is_done == StatusDone::IsNotOk {
            return Err("FilletSurf_Builder::LastParameter".to_string());
        }
        Ok(self.int_build.last_parameter())
    }

    /// `StartSectionStatus()` (`cxx:295-301`).
    pub fn start_section_status(&self) -> Result<StatusType, String> {
        if self.is_done == StatusDone::IsNotOk {
            return Err("FilletSurf_Builder::StartSectionStatus".to_string());
        }
        Ok(self.int_build.start_section_status())
    }

    /// `EndSectionStatus()` (`cxx:314-320`). OCCT raises the message
    /// `"FilletSurf_Builder::StartSectionStatus"` here as well; the port names
    /// the member it is in.
    pub fn end_section_status(&self) -> Result<StatusType, String> {
        if self.is_done == StatusDone::IsNotOk {
            return Err("FilletSurf_Builder::EndSectionStatus".to_string());
        }
        Ok(self.int_build.end_section_status())
    }

    /// `NbSection(IndexSurf)` (`cxx:347-358`).
    pub fn nb_section(&self, index_surf: usize) -> Result<usize, String> {
        if self.is_done == StatusDone::IsNotOk {
            return Err("FilletSurf_Builder::NbSection".to_string());
        }
        self.check_index(index_surf, "FilletSurf_Builder::NbSection")?;
        Ok(self.int_build.nb_section(index_surf))
    }

    /// `Section(IndexSurf, IndexSec)` (`cxx:366-380`): the section as a
    /// trimmed circular arc.
    pub fn section(
        &self,
        index_surf: usize,
        index_sec: usize,
    ) -> Result<GeomTrimmedCurve, String> {
        self.check_index(index_surf, "FilletSurf_Builder::Section NbSurface")?;
        let nb_sec = self.nb_section(index_surf)?;
        if index_sec == 0 || index_sec > nb_sec {
            return Err("FilletSurf_Builder::Section NbSection".to_string());
        }
        self.int_build
            .section(index_surf, index_sec)
            .ok_or_else(|| "FilletSurf_Builder::Section: no section".to_string())
    }

    /// The `if ((Index < 1) || (Index > NbSurface())) throw
    /// Standard_OutOfRange(...)` guard the indexed accessors share.
    fn check_index(&self, index: usize, what: &str) -> Result<(), String> {
        let nb = self.nb_surface()?;
        if index == 0 || index > nb {
            return Err(what.to_string());
        }
        Ok(())
    }
}
