//! `HLRBRep_Curve` (`ModelingAlgorithms/TKHLR/HLRBRep/HLRBRep_Curve.hxx:44-213`,
//! `.cxx:38-632`, `.lxx:28-169`).
//!
//! "Defines a 2d curve by projection of a 3D curve on a plane with an
//! optional perspective transformation." The class owns the 3D adaptor
//! (`BRepAdaptor_Curve myCurve`), the cached projected type (`myType`) and the
//! five scalars of the linear parameter map (`myOX`, `myOZ`, `myVX`, `myVZ`,
//! `myOF`) that `Parameter2d` / `Parameter3d` invert.
//!
//! The port splits the members over three submodules:
//! * [`update`]: `Parameter2d` / `Parameter3d` / `Update` / `UpdateMinMax` /
//!   `Z`,
//! * [`projection`]: the 2D evaluation (`D0`..`D3`, `DN`, `Value`, `Tangent`,
//!   `Line`, `Circle`, `Ellipse`, `Hyperbola`, `Parabola`) and the 3D ones
//!   (`Value3D`, `D0`, `D1`),
//! * [`bspline_queries`]: `IsRational`, `Degree`, `NbPoles`, `Poles`,
//!   `PolesAndWeights`, `NbKnots`, `Knots`, `Multiplicities`.
//!
//! `myProj` is a non-owning `const HLRAlgo_Projector*` in OCCT; the port holds
//! an `Arc` so that `HLRBRep_EdgeData` / `HLRBRep_FaceData` can carry the
//! reference without a lifetime parameter. OCCT dereferences the pointer
//! unconditionally; here a missing projector collapses all projection queries
//! to the identity (`Perspective()` answers `false`, `Project` copies `X`/`Y`,
//! `Transformation()` is the identity), which is the state the default
//! constructor leaves the object in.

use std::sync::Arc;

use occt_core::gp::GpTrsf;
use occt_geom::Curve;

use occt_topo::brep_tool::BRepTool;
use occt_topo::geom_abs_curve_type::{curve_kind, GeomAbsCurveType};
use crate::hlr_algo::HlrAlgoProjector;
use occt_topo::shape::Edge;

mod bspline_queries;
mod projection;
mod update;

/// `HLRBRep_Curve`.
#[derive(Clone)]
pub struct HlrBRepCurve {
    /// `BRepAdaptor_Curve myCurve` (`hxx:204`): the 3D curve of the edge plus
    /// its parameter range. `BRepAdaptor_Curve::Edge()` (`BRepAdaptor_Curve.hxx`)
    /// is read by `HLRBRep_EdgeFaceTool::UVPoint` (`HLRBRep_EdgeFaceTool.cxx:73`,
    /// `:97`) through `HLRBRep_Curve::Curve().Edge()`, so the edge itself is
    /// kept here as `myCurve` keeps it.
    curve: Option<Arc<dyn Curve>>,
    /// `BRepAdaptor_Curve::Edge()`: the edge `set_curve` loaded.
    edge: Option<Edge>,
    first_parameter: f64,
    last_parameter: f64,
    /// `GeomAbs_CurveType myType` (`hxx:205`): the type of the *projected*
    /// curve, refreshed by [`HlrBRepCurve::update`].
    type_: GeomAbsCurveType,
    /// `const HLRAlgo_Projector* myProj` (`hxx:206`).
    proj: Option<Arc<HlrAlgoProjector>>,
    /// `double myOX` (`hxx:207`).
    ox: f64,
    /// `double myOZ` (`hxx:208`).
    oz: f64,
    /// `double myVX` (`hxx:209`).
    vx: f64,
    /// `double myVZ` (`hxx:210`).
    vz: f64,
    /// `double myOF` (`hxx:211`).
    of: f64,
}

impl Default for HlrBRepCurve {
    /// `HLRBRep_Curve::HLRBRep_Curve()` (`cxx:38`): `= default`. The five
    /// scalars and `myProj` are left as the C++ default member initialisers
    /// leave them (indeterminate / null); `myType` reads as `GeomAbs_Line`,
    /// the first enumerator of `GeomAbs_CurveType`.
    fn default() -> Self {
        Self {
            curve: None,
            edge: None,
            first_parameter: 0.0,
            last_parameter: 0.0,
            type_: GeomAbsCurveType::Line,
            proj: None,
            ox: 0.0,
            oz: 0.0,
            vx: 0.0,
            vz: 0.0,
            of: 0.0,
        }
    }
}

impl HlrBRepCurve {
    /// `HLRBRep_Curve()` (`cxx:38`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `Projector(const HLRAlgo_Projector* Proj)` (`hxx:54`).
    pub fn set_projector(&mut self, proj: Option<Arc<HlrAlgoProjector>>) {
        self.proj = proj;
    }

    /// The stored `myProj` pointer.
    pub fn projector(&self) -> Option<&Arc<HlrAlgoProjector>> {
        self.proj.as_ref()
    }

    /// `Curve(const TopoDS_Edge& E)` (`cxx:44-47`): `myCurve.Initialize(E)`.
    /// The adaptor keeps the edge's parameter range, which the port stores
    /// alongside the curve.
    pub fn set_curve(&mut self, e: &Edge) {
        match BRepTool::edge_curve(e) {
            Some(c) => {
                let (first, last) = BRepTool::edge_parameters(e);
                self.curve = Some(c);
                self.edge = Some(e.clone());
                self.first_parameter = first;
                self.last_parameter = last;
            }
            None => {
                self.curve = None;
                self.edge = Some(e.clone());
                self.first_parameter = 0.0;
                self.last_parameter = 0.0;
            }
        }
    }

    /// `BRepAdaptor_Curve::Edge()`: the edge this curve was loaded from
    /// (`HLRBRep_EdgeFaceTool.cxx:73`, `:97` read it through
    /// `HLRBRep_Curve::Curve().Edge()`).
    pub fn edge(&self) -> Option<&Edge> {
        self.edge.as_ref()
    }

    /// `Curve()` (`lxx:28-31`) / `GetCurve()` (`lxx:36-39`).
    pub fn curve(&self) -> Option<&Arc<dyn Curve>> {
        self.curve.as_ref()
    }

    /// `GetType()` (`lxx:138-141`): the cached projected type.
    pub fn get_type(&self) -> GeomAbsCurveType {
        self.type_
    }

    /// `HLRBRep_BCurveTool::GetType(myCurve)`: the 3D adaptor's type
    /// (`GeomAdaptor_Curve::GetType`).
    pub(crate) fn curve_kind_3d(&self) -> GeomAbsCurveType {
        match &self.curve {
            Some(c) => curve_kind(c.as_ref()),
            None => GeomAbsCurveType::OtherCurve,
        }
    }

    /// `Transformation()` of `myProj` (`HLRAlgo_Projector.cxx:392-396`); the
    /// identity when no projector was set.
    pub(crate) fn transformation(&self) -> GpTrsf {
        match &self.proj {
            Some(p) => p.transformation().clone(),
            None => GpTrsf::default(),
        }
    }

    /// `Perspective()` of `myProj` (`lxx:33-36`); `false` without a projector.
    pub(crate) fn perspective(&self) -> bool {
        match &self.proj {
            Some(p) => p.perspective(),
            None => false,
        }
    }

    /// `Focus()` of `myProj` (`lxx:58-63`); thrown as `Err` only when the
    /// projector is not a perspective one, as in OCCT.
    pub(crate) fn focus(&self) -> Result<f64, &'static str> {
        match &self.proj {
            Some(p) => p.focus(),
            None => Err("HlrAlgoProjector::Not a Perspective"),
        }
    }

    /// `FirstParameter()` (`lxx:64-67`): the 2D first parameter,
    /// `Parameter2d(BCurveTool::FirstParameter(myCurve))`.
    pub fn first_parameter(&self) -> f64 {
        self.parameter2d(self.first_parameter)
    }

    /// `LastParameter()` (`lxx:71-74`).
    pub fn last_parameter(&self) -> f64 {
        self.parameter2d(self.last_parameter)
    }

    /// `Continuity()` (`lxx:80-83`).
    pub fn continuity(&self) -> u8 {
        match &self.curve {
            Some(c) => c.continuity(),
            None => 0,
        }
    }

    /// `NbIntervals(S)` (`lxx:87-90`).
    pub fn nb_intervals(&self, s: u8) -> i32 {
        match &self.curve {
            Some(c) => c.nb_intervals(s),
            None => 0,
        }
    }

    /// `Intervals(T, S)` (`lxx:95-98`).
    pub fn intervals(&self, s: u8) -> Vec<f64> {
        match &self.curve {
            Some(c) => c.parameter_intervals(s),
            None => Vec::new(),
        }
    }

    /// `IsClosed()` (`lxx:102-105`).
    pub fn is_closed(&self) -> bool {
        match &self.curve {
            Some(c) => c.is_closed(),
            None => false,
        }
    }

    /// `IsPeriodic()` (`lxx:109-112`).
    pub fn is_periodic(&self) -> bool {
        match &self.curve {
            Some(c) => c.is_periodic(),
            None => false,
        }
    }

    /// `Period()` (`lxx:116-119`).
    pub fn period(&self) -> f64 {
        match &self.curve {
            Some(c) => c.period(),
            None => 0.0,
        }
    }

    /// `Resolution(R3d)` (`lxx:131-134`).
    pub fn resolution(&self, r3d: f64) -> f64 {
        match &self.curve {
            Some(c) => c.resolution(r3d),
            None => 0.0,
        }
    }
}
