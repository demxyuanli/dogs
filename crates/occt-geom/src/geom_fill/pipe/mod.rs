//! `GeomFill_Pipe` — a surface swept from a section along a path.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_Pipe.cxx:229-1102` and
//! `GeomFill_Pipe.hxx:89-320` / `GeomFill_Pipe.lxx:17-56` (OCCT 8.0.0).
//!
//! The class is split along its OCCT responsibilities:
//! - this file: the fields, the constructors and the accessors,
//! - [`check_sense`]: the file-local `CheckSense` (`cxx:79-222`),
//! - [`init`]: the `Init` family (`cxx:364-728`),
//! - [`perform`]: `Perform` / `KPartT4` / `ApproxSurf` (`cxx:747-1102`).
//!
//! Translation notes:
//! - `occ::handle<Geom_Curve>` / `occ::handle<Geom_Surface>` become
//!   `Arc<dyn Curve>` / `Arc<dyn Surface>`; a null handle is `None`.
//! - OCCT's `Standard_ConstructionError` throws become `Err(String)`.
//! - The `GeomAdaptor_Curve` wrappers (`myAdpPath`, `myAdpFirstSect`,
//!   `myAdpLastSect`) are the curves themselves: the adaptor only forwards
//!   `GetType` / `FirstParameter` / `LastParameter` / `Value` / `Line` /
//!   `Circle`, all of which [`crate::curve::Curve`] answers.
//! - `GeomFill_SectionPlacement::Perform` takes `&mut self` here, so every
//!   placement is a local `let mut` (OCCT mutates a stack object the same way).

use std::sync::Arc;

use crate::curve::Curve;
use crate::geom_fill::location_law::LocationLaw;
use crate::geom_fill::pipe_error::GeomFillPipeError;
use crate::geom_fill::section_law::SectionLaw;
use crate::surface::Surface;

mod check_sense;
mod init;
mod perform;

pub use check_sense::check_sense;

/// `GeomFill_Pipe` (`GeomFill_Pipe.hxx:89-320`).
pub struct Pipe {
    /// `myStatus` (`hxx:301`).
    status: GeomFillPipeError,
    /// `myRadius` (`hxx:302`).
    radius: f64,
    /// `myError` (`hxx:303`).
    error: f64,
    /// `myAdpPath` (`hxx:304`).
    adp_path: Option<Arc<dyn Curve>>,
    /// `myAdpFirstSect` (`hxx:305`).
    adp_first_sect: Option<Arc<dyn Curve>>,
    /// `myAdpLastSect` (`hxx:306`).
    adp_last_sect: Option<Arc<dyn Curve>>,
    /// `mySurface` (`hxx:310`).
    surface: Option<Arc<dyn Surface>>,
    /// `myLoc` (`hxx:311`).
    loc: Option<Arc<dyn LocationLaw>>,
    /// `mySec` (`hxx:312`).
    sec: Option<Arc<dyn SectionLaw>>,
    /// `myType` (`hxx:313`).
    pipe_type: i32,
    /// `myExchUV` (`hxx:314`).
    exch_uv: bool,
    /// `myKPart` (`hxx:315`).
    k_part: bool,
    /// `myPolynomial` (`hxx:316`).
    polynomial: bool,
}

impl Pipe {
    /// `GeomFill_Pipe()` (`cxx:229-234`).
    pub fn new() -> Self {
        let mut pipe = Self {
            status: GeomFillPipeError::PipeNotOk,
            radius: 0.0,
            error: 0.0,
            adp_path: None,
            adp_first_sect: None,
            adp_last_sect: None,
            surface: None,
            loc: None,
            sec: None,
            pipe_type: 0,
            exch_uv: false,
            k_part: false,
            polynomial: false,
        };
        pipe.init();
        pipe
    }

    /// `GeomFill_Pipe(Path, Radius)` (`cxx:237-245`).
    pub fn with_path_radius(path: Arc<dyn Curve>, radius: f64) -> Result<Self, String> {
        let mut pipe = Self::new();
        pipe.init_path_radius(path, radius)?;
        Ok(pipe)
    }

    /// `GeomFill_Pipe(Path, FirstSect, Option)` (`cxx:249-258`).
    pub fn with_path_section_option(
        path: Arc<dyn Curve>,
        first_sect: Arc<dyn Curve>,
        option: crate::geom_fill::trihedron::GeomFillTrihedron,
    ) -> Result<Self, String> {
        let mut pipe = Self::new();
        pipe.init_path_section_option(path, first_sect, option)?;
        Ok(pipe)
    }

    /// `GeomFill_Pipe(Path, FirstSect, LastSect)` (`cxx:276-286`).
    pub fn with_path_first_last(
        path: Arc<dyn Curve>,
        first_sect: Arc<dyn Curve>,
        last_sect: Arc<dyn Curve>,
    ) -> Result<Self, String> {
        let mut pipe = Self::new();
        pipe.init_path_first_last(path, first_sect, last_sect)?;
        Ok(pipe)
    }

    /// `GeomFill_Pipe(Path, NSections)` (`cxx:288-296`).
    pub fn with_path_sections(
        path: Arc<dyn Curve>,
        nsections: Vec<Arc<dyn Curve>>,
    ) -> Result<Self, String> {
        let mut pipe = Self::new();
        pipe.init_path_sections(path, nsections)?;
        Ok(pipe)
    }

    /// `GeomFill_Pipe(Path, Curve1, Direction)` (`cxx:303-308`). OCCT does not
    /// call `Init()` first here; `Init(Path, FirstSect, Direction)` does it
    /// itself (`cxx:589`).
    pub fn with_path_first_section_dir(
        path: Arc<dyn Curve>,
        first_sect: Arc<dyn Curve>,
        direction: occt_core::gp::GpDir,
    ) -> Result<Self, String> {
        let mut pipe = Self::new();
        pipe.init_path_first_section_dir(path, first_sect, direction)?;
        Ok(pipe)
    }

    /// `GeomFill_Pipe(Path, Curve1, Curve2, Radius)` (`cxx:312-323`): both
    /// curves become `GeomAdaptor_Curve`s before the adaptor-taking `Init`.
    pub fn with_path_curve1_curve2_radius(
        path: Arc<dyn Curve>,
        curve1: Arc<dyn Curve>,
        curve2: Arc<dyn Curve>,
        radius: f64,
    ) -> Result<Self, String> {
        let mut pipe = Self::new();
        pipe.init_adaptors(path, curve1, curve2, radius);
        Ok(pipe)
    }

    /// `GeomFill_Pipe(Path, Guide, FirstSect, ByACR, rotat)` (`cxx:346-361`).
    /// Needs `GeomFill_LocationGuide`, see [`Pipe::init_with_guide`].
    pub fn with_guide(
        path: Arc<dyn Curve>,
        guide: Arc<dyn Curve>,
        first_sect: Arc<dyn Curve>,
        by_acr: bool,
        rotat: bool,
    ) -> Result<Self, String> {
        let mut pipe = Self::new();
        pipe.init_with_guide(path, guide, first_sect, by_acr, rotat)?;
        Ok(pipe)
    }

    /// `GeomFill_Pipe::Init()` (`cxx:414-428`).
    pub fn init(&mut self) {
        self.pipe_type = 0;
        self.error = 0.0;
        self.radius = 0.0;
        self.k_part = true;
        self.polynomial = false;
        self.adp_path = None;
        self.adp_first_sect = None;
        self.adp_last_sect = None;
        self.loc = None;
        self.sec = None;
    }

    /// `GeomFill_Pipe::Surface()` (`lxx:19-22`).
    pub fn surface(&self) -> Option<Arc<dyn Surface>> {
        self.surface.clone()
    }

    /// `GeomFill_Pipe::ExchangeUV()` (`lxx:27-30`).
    pub fn exchange_uv(&self) -> bool {
        self.exch_uv
    }

    /// `GeomFill_Pipe::GenerateParticularCase(const bool)` (`lxx:35-38`).
    pub fn set_generate_particular_case(&mut self, b: bool) {
        self.k_part = b;
    }

    /// `GeomFill_Pipe::GenerateParticularCase()` (`lxx:43-46`).
    pub fn generate_particular_case(&self) -> bool {
        self.k_part
    }

    /// `GeomFill_Pipe::ErrorOnSurf()` (`lxx:49-52`).
    pub fn error_on_surf(&self) -> f64 {
        self.error
    }

    /// `GeomFill_Pipe::IsDone()` (`lxx:54-57`): `myStatus == GeomFill_PipeOk`.
    pub fn is_done(&self) -> bool {
        self.status == GeomFillPipeError::PipeOk
    }

    /// `GeomFill_Pipe::GetStatus()` (`hxx:291`).
    pub fn status(&self) -> GeomFillPipeError {
        self.status
    }
}

impl Default for Pipe {
    fn default() -> Self {
        Self::new()
    }
}
