//! `Contap_Contour` - the contour computation of one surface: the walking arm
//! for a general surface and the analytical arm for the four quadrics.
//!
//! Sources: `Contap_Contour.hxx:33-116`, `Contap_Contour.cxx:1-2389`,
//! `Contap_Contour.lxx:1-40`.
//!
//! Split by OCCT responsibility:
//! - [`recadre`] - `Recadre` (`cxx:229-281`).
//! - [`line_constructor`] - `LineConstructor` (`cxx:283-511`).
//! - [`keep_inside_points`] - `KeepInsidePoints` (`cxx:515-554`).
//! - [`compute_tangency`] - `ComputeTangency` (`cxx:556-891`).
//! - [`transitions`] - `ComputeTransitionOnLine` (`cxx:893-956`),
//!   `ComputeTransitionOngpLine` (`cxx:2100-2126`),
//!   `ComputeTransitionOngpCircle` (`cxx:2128-2154`).
//! - [`process_segments`] - `ProcessSegments` (`cxx:958-1118`).
//! - [`internal_points`] - `ComputeInternalPointsOnRstr` (`cxx:1120-1284`),
//!   `ComputeInternalPoints` (`cxx:1286-1539`).
//! - [`points_on_line`] - `FindLine` (`cxx:1973-2011`), `PutPointsOnLine`
//!   (`cxx:2013-2098`).
//! - [`perform_walking`] - `Perform(Domain)` (`cxx:1541-1971`).
//! - [`perform_ana`] - `PerformAna(Domain)` (`cxx:2156-2389`).

use occt_core::gp::{GpDir, GpPnt, GpVec};

use occt_topo::brep_top_adaptor_topol_tool::TopolTool;
use occt_topo::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::arc_function::ContapArcFunction;
use super::line::ContapLine;
use super::search_inside::ContapTheSearchInside;
use super::surf_function::ContapSurfFunction;
use super::the_search::ContapTheSearch;
use super::types::ContapTFunction;

pub mod compute_tangency;
pub mod internal_points;
pub mod keep_inside_points;
pub mod line_constructor;
pub mod perform_ana;
pub mod perform_walking;
pub mod points_on_line;
pub mod process_segments;
pub mod recadre;
pub mod transitions;

/// `static const double Tolpetit` (`cxx:46`): squared distance, `FindLine`.
pub(crate) const TOLPETIT: f64 = 1.0e-10;

/// `static const double tole` (`cxx:48`): the tangent/transition threshold of
/// `ComputeTangency`.
pub(crate) const TOLE: f64 = 5.0e-6;

/// `Contap_Contour` (`hxx:33-116`).
pub struct ContapContour {
    /// `done` (`hxx:106`).
    done: bool,
    /// `slin` (`hxx:107`): the computed contour lines.
    slin: Vec<ContapLine>,
    /// `solrst` (`hxx:108`): the solutions on the restrictions.
    solrst: ContapTheSearch,
    /// `solins` (`hxx:109`): the interior solutions.
    solins: ContapTheSearchInside,
    /// `mySFunc` (`hxx:110`).
    my_sfunc: ContapSurfFunction,
    /// `myAFunc` (`hxx:111`).
    my_afunc: ContapArcFunction,
    /// `modeset` (`hxx:112`): `Init` was called.
    modeset: bool,
}

impl Default for ContapContour {
    /// `Contap_Contour()` (`cxx:50-53`).
    fn default() -> Self {
        Self::new()
    }
}

impl ContapContour {
    /// `Contap_Contour()` (`cxx:50-53`): `done` and `modeset` false, both
    /// functions default.
    pub fn new() -> Self {
        Self {
            done: false,
            slin: Vec::new(),
            solrst: ContapTheSearch::new(),
            solins: ContapTheSearchInside::new(),
            my_sfunc: ContapSurfFunction::new(),
            my_afunc: ContapArcFunction::new(),
            modeset: false,
        }
    }

    /// `Contap_Contour(const gp_Vec& Direction)` (`cxx:56-62`).
    ///
    /// `Set(const gp_Vec&)` binds OCCT's `Set(const gp_Dir&)` through the
    /// implicit `gp_Vec -> gp_Dir` conversion; the port converts explicitly and
    /// answers `Err` where that conversion would throw
    /// (`gp_Vec::Normalize`, `Standard_VectorWithNullMagnitude`).
    pub fn with_direction(direction: &GpVec) -> Result<Self, &'static str> {
        let mut c = Self::new();
        c.init_direction(direction)?;
        Ok(c)
    }

    /// `Contap_Contour(const gp_Vec& Direction, const double Angle)`
    /// (`cxx:66-72`).
    pub fn with_direction_angle(direction: &GpVec, angle: f64) -> Result<Self, &'static str> {
        let mut c = Self::new();
        c.init_direction_angle(direction, angle)?;
        Ok(c)
    }

    /// `Contap_Contour(const gp_Pnt& Eye)` (`cxx:76-82`).
    pub fn with_eye(eye: &GpPnt) -> Self {
        let mut c = Self::new();
        c.init_eye(eye);
        c
    }

    /// `Contap_Contour(Surf, Domain, Direction)` (`cxx:86-95`).
    pub fn with_surface_direction(
        surf: SurfaceAdaptor,
        domain: &mut TopolTool,
        direction: &GpVec,
    ) -> Result<Self, &'static str> {
        let mut c = Self::new();
        c.init_direction(direction)?;
        c.perform(surf, domain);
        Ok(c)
    }

    /// `Contap_Contour(Surf, Domain, Direction, Angle)` (`cxx:97-106`).
    pub fn with_surface_direction_angle(
        surf: SurfaceAdaptor,
        domain: &mut TopolTool,
        direction: &GpVec,
        angle: f64,
    ) -> Result<Self, &'static str> {
        let mut c = Self::new();
        c.init_direction_angle(direction, angle)?;
        c.perform(surf, domain);
        Ok(c)
    }

    /// `Contap_Contour(Surf, Domain, Eye)` (`cxx:108-117`).
    pub fn with_surface_eye(surf: SurfaceAdaptor, domain: &mut TopolTool, eye: &GpPnt) -> Self {
        let mut c = Self::new();
        c.init_eye(eye);
        c.perform(surf, domain);
        c
    }

    /// `Init(const gp_Vec& Direction)` (`cxx:120-126`).
    pub fn init_direction(&mut self, direction: &GpVec) -> Result<(), &'static str> {
        let dir = GpDir::from_vec(direction)?;
        self.done = false;
        self.modeset = true;
        self.my_sfunc.set_direction(dir);
        self.my_afunc.set_direction(dir);
        Ok(())
    }

    /// `Init(const gp_Vec& Direction, const double Angle)` (`cxx:129-135`).
    pub fn init_direction_angle(
        &mut self,
        direction: &GpVec,
        angle: f64,
    ) -> Result<(), &'static str> {
        let dir = GpDir::from_vec(direction)?;
        self.done = false;
        self.modeset = true;
        self.my_sfunc.set_direction_angle(dir, angle);
        self.my_afunc.set_direction_angle(dir, angle);
        Ok(())
    }

    /// `Init(const gp_Pnt& Eye)` (`cxx:137-142`).
    pub fn init_eye(&mut self, eye: &GpPnt) {
        self.done = false;
        self.modeset = true;
        self.my_sfunc.set_eye(*eye);
        self.my_afunc.set_eye(*eye);
    }

    /// `Perform(Surf, Domain)` (`cxx:145-171`): the surface type selects the
    /// analytical arm (`Plane`, `Sphere`, `Cylinder`, `Cone`) or the walking
    /// arm.
    ///
    /// OCCT throws `Standard_ConstructionError` when `modeset` is false; the
    /// port answers `Err`.
    pub fn perform(&mut self, surf: SurfaceAdaptor, domain: &mut TopolTool) -> Result<(), String> {
        if !self.modeset {
            return Err("Contap_Contour::Perform: Standard_ConstructionError".to_string());
        }
        self.my_sfunc.set_surface(surf.clone());
        self.my_afunc.set_surface(surf.clone());

        match surf.get_type() {
            occt_topo::geom_bnd_lib_surface3d::GeomAbsSurfaceType::Plane
            | occt_topo::geom_bnd_lib_surface3d::GeomAbsSurfaceType::Sphere
            | occt_topo::geom_bnd_lib_surface3d::GeomAbsSurfaceType::Cylinder
            | occt_topo::geom_bnd_lib_surface3d::GeomAbsSurfaceType::Cone => {
                self.perform_ana(domain);
            }
            _ => self.perform_domain(domain),
        }
        Ok(())
    }

    /// `Perform(Surf, Domain, Direction)` (`cxx:173-180`).
    pub fn perform_direction(
        &mut self,
        surf: SurfaceAdaptor,
        domain: &mut TopolTool,
        direction: &GpVec,
    ) -> Result<(), String> {
        self.init_direction(direction).map_err(str::to_string)?;
        self.perform(surf, domain)
    }

    /// `Perform(Surf, Domain, Direction, Angle)` (`cxx:182-190`).
    pub fn perform_direction_angle(
        &mut self,
        surf: SurfaceAdaptor,
        domain: &mut TopolTool,
        direction: &GpVec,
        angle: f64,
    ) -> Result<(), String> {
        self.init_direction_angle(direction, angle)
            .map_err(str::to_string)?;
        self.perform(surf, domain)
    }

    /// `Perform(Surf, Domain, Eye)` (`cxx:192-198`).
    pub fn perform_eye(
        &mut self,
        surf: SurfaceAdaptor,
        domain: &mut TopolTool,
        eye: &GpPnt,
    ) -> Result<(), String> {
        self.init_eye(eye);
        self.perform(surf, domain)
    }

    /// `IsDone()` (`lxx:19-22`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `IsEmpty()` (`lxx:24-27`).
    pub fn is_empty(&self) -> bool {
        self.slin.is_empty()
    }

    /// `NbLines()` (`lxx:29-32`).
    pub fn nb_lines(&self) -> usize {
        self.slin.len()
    }

    /// `Line(Index)` (`lxx:34-37`), 1-based; `None` where OCCT raises
    /// `Standard_OutOfRange`.
    pub fn line(&self, index: usize) -> Option<&ContapLine> {
        index.checked_sub(1).and_then(|i| self.slin.get(i))
    }

    /// `SurfaceFunction()` (`lxx:39-42`): the internal surface function, kept
    /// alive so callers can evaluate tangents on the computed lines.
    pub fn surface_function(&self) -> &ContapSurfFunction {
        &self.my_sfunc
    }

    /// The mutable `SurfaceFunction()`.
    pub fn surface_function_mut(&mut self) -> &mut ContapSurfFunction {
        &mut self.my_sfunc
    }

    /// `FunctionType()` of the internal surface function - not an OCCT entry,
    /// it exposes `mySFunc.FunctionType()` for the parity harness.
    pub fn function_type(&self) -> ContapTFunction {
        self.my_sfunc.function_type()
    }

    /// The lines in insertion order (`slin` access for the parity harness).
    pub fn lines(&self) -> &[ContapLine] {
        &self.slin
    }

    /// The segments of the search, as `ProcessSegments` sees them.
    pub(crate) fn solrst(&self) -> &ContapTheSearch {
        &self.solrst
    }
}
