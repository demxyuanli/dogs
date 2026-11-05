//! `ChFi2d_FilletAPI` — the dispatcher between the analytical and the
//! iterative fillet algorithm.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi2d/ChFi2d_FilletAPI.cxx:22-136` and
//! `.hxx` (OCCT 8.0.0).
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `ChFi2d_FilletAPI()` (`cxx:22-26`) | 5 | [`FilletApi::new`] |
//! | `ChFi2d_FilletAPI(wire, plane)` (`cxx:30-34`) | 5 | [`FilletApi::from_wire`] |
//! | `ChFi2d_FilletAPI(e1, e2, plane)` (`cxx:37-42`) | 6 | [`FilletApi::from_edges`] |
//! | `Init(wire, plane)` (`cxx:45-73`) | 29 | [`FilletApi::init_wire`] |
//! | `Init(e1, e2, plane)` (`cxx:76-84`) | 9 | [`FilletApi::init_edges`] |
//! | `Perform(radius)` (`cxx:87-90`) | 4 | [`FilletApi::perform`] |
//! | `NbResults(point)` (`cxx:93-97`) | 5 | [`FilletApi::nb_results`] |
//! | `Result(...)` (`cxx:102-109`) | 8 | [`FilletApi::result`] |
//! | `IsAnalytical(e1, e2)` (`cxx:113-135`) | 23 | [`FilletApi::is_analytical`] |
//!
//! Translation notes:
//! - `myIsAnalytical ? myAnaFilletAlgo.X(...) : myFilletAlgo.X(...)` is a plain
//!   `if` over the two stored algorithms; OCCT holds both by value, and only
//!   the chosen one gets `Init` called (`cxx:72`, `cxx:83-84`),
//! - `NbResults` answers `1` on the analytical arm without consulting
//!   [`AnaFilletAlgo`] (`cxx:95`),
//! - `Result` on the analytical arm calls
//!   [`AnaFilletAlgo::result`], which fills all three edges or none; the port
//!   maps that single `Option` back onto the three out-parameters,
//! - `BRepAdaptor_Curve(E)` is [`occt_geom::projlib::curve_type`] over the
//!   edge's stored 3D curve, so `GeomAbs_Line` / `GeomAbs_Circle` are
//!   [`occt_core::kernel::geomabs::CurveType`] variants,
//! - `Precision::SquareConfusion()` is
//!   [`occt_core::precision::SQUARE_CONFUSION`]; the `AC1.Value(...)` endpoints
//!   are `d0` on the edge's stored curve evaluated at the edge's range ends.

use occt_core::gp::{GpPln, GpPnt};
use occt_core::kernel::geomabs::CurveType;
use occt_core::precision::SQUARE_CONFUSION;

use crate::brep_tool::BRepTool;
use crate::iterator::ShapeIterator;
use crate::shape::{Edge, Wire};

use super::ana_fillet_algo::AnaFilletAlgo;
use super::fillet_algo::{FilletAlgo, FilletAlgoResult};

/// `ChFi2d_FilletAPI` (`ChFi2d_FilletAPI.hxx`).
#[derive(Clone, Default)]
pub struct FilletApi {
    /// `myAnaFilletAlgo`.
    ana_fillet_algo: AnaFilletAlgo,
    /// `myFilletAlgo`.
    fillet_algo: FilletAlgo,
    /// `myIsAnalytical` (`cxx:24`).
    is_analytical: bool,
}

impl FilletApi {
    /// `ChFi2d_FilletAPI()` (`cxx:22-26`): `myIsAnalytical = false`.
    pub fn new() -> Self {
        Self::default()
    }

    /// `ChFi2d_FilletAPI(theWire, thePlane)` (`cxx:30-34`).
    pub fn from_wire(wire: &Wire, plane: &GpPln) -> Result<Self, String> {
        let mut api = Self::new();
        api.init_wire(wire, plane)?;
        Ok(api)
    }

    /// `ChFi2d_FilletAPI(theEdge1, theEdge2, thePlane)` (`cxx:37-42`).
    pub fn from_edges(edge1: &Edge, edge2: &Edge, plane: &GpPln) -> Result<Self, String> {
        let mut api = Self::new();
        api.init_edges(edge1, edge2, plane)?;
        Ok(api)
    }

    /// `Init(theWire, thePlane)` (`cxx:45-73`).
    pub fn init_wire(&mut self, wire: &Wire, plane: &GpPln) -> Result<(), String> {
        // `Decide whether we may apply an analytical solution` (`cxx:48-69`).
        let mut e1: Option<Edge> = None;
        let mut e2: Option<Edge> = None;
        for child in ShapeIterator::of_shape(&wire.0) {
            let edge = Edge::wrap(child)
                .ok_or_else(|| "ChFi2d_FilletAPI::Init: wire child is not an edge".to_string())?;
            if e1.is_none() {
                e1 = Some(edge);
            } else if e2.is_none() {
                e2 = Some(edge);
            } else {
                break;
            }
        }
        if let (Some(edge1), Some(edge2)) = (e1.as_ref(), e2.as_ref()) {
            self.is_analytical = Self::is_analytical(edge1, edge2);
        }

        // `myIsAnalytical ? myAnaFilletAlgo.Init(...) : myFilletAlgo.Init(...)`
        // (`cxx:72`). OCCT keeps both objects and only initialises the chosen
        // one.
        if self.is_analytical {
            self.ana_fillet_algo.init_wire(wire, plane)
        } else {
            self.fillet_algo.init_wire(wire, plane)
        }
    }

    /// `Init(theEdge1, theEdge2, thePlane)` (`cxx:76-84`).
    pub fn init_edges(&mut self, edge1: &Edge, edge2: &Edge, plane: &GpPln) -> Result<(), String> {
        // `Decide whether we may apply an analytical solution` (`cxx:80`).
        self.is_analytical = Self::is_analytical(edge1, edge2);

        if self.is_analytical {
            self.ana_fillet_algo.init_edges(edge1, edge2, plane)
        } else {
            self.fillet_algo.init_edges(edge1, edge2, plane)
        }
    }

    /// `Perform(theRadius)` (`cxx:87-90`): `true` when at least one result was
    /// found.
    pub fn perform(&mut self, radius: f64) -> bool {
        if self.is_analytical {
            self.ana_fillet_algo.perform(radius)
        } else {
            self.fillet_algo.perform(radius)
        }
    }

    /// `NbResults(thePoint)` (`cxx:93-97`): always `1` on the analytical arm.
    pub fn nb_results(&mut self, the_point: &GpPnt) -> i32 {
        if self.is_analytical {
            1
        } else {
            self.fillet_algo.nb_results(the_point)
        }
    }

    /// `Result(thePoint, theEdge1, theEdge2, iSolution)` (`cxx:102-109`).
    ///
    /// The analytical arm ignores `thePoint` and `iSolution` (`cxx:107`); both
    /// arms answer the three edges together, so one `Option` carries them.
    pub fn result(&mut self, the_point: &GpPnt, i_solution: i32) -> Option<FilletAlgoResult> {
        if self.is_analytical {
            // `myAnaFilletAlgo.Result(theEdge1, theEdge2)` (`cxx:107`).
            self.ana_fillet_algo
                .result()
                .map(|r| FilletAlgoResult {
                    fillet: r.fillet,
                    edge1: r.shrinke1,
                    edge2: r.shrinke2,
                })
        } else {
            self.fillet_algo.result(the_point, i_solution)
        }
    }

    /// `myIsAnalytical`: whether the last `Init` picked the analytical
    /// algorithm.
    pub fn is_analytical_solution(&self) -> bool {
        self.is_analytical
    }

    /// `IsAnalytical(theEdge1, theEdge2)` (`cxx:113-135`): the analytical
    /// solution applies to linear and circular edges that share an endpoint.
    pub fn is_analytical(edge1: &Edge, edge2: &Edge) -> bool {
        let (Some(curve1), Some(curve2)) =
            (BRepTool::edge_curve(edge1), BRepTool::edge_curve(edge2))
        else {
            return false;
        };
        let (type1, type2) = (
            occt_geom::projlib::curve_type(curve1.as_ref()),
            occt_geom::projlib::curve_type(curve2.as_ref()),
        );
        // `(AC1.GetType() == GeomAbs_Line || AC1.GetType() == GeomAbs_Circle)
        // && (AC2.GetType() == ... )` (`cxx:118-119`).
        let line_or_circle =
            |t: CurveType| t == CurveType::Line || t == CurveType::Circle;
        if !line_or_circle(type1) || !line_or_circle(type2) {
            return false;
        }

        // `Now check whether they have a common point` (`cxx:121-131`).
        // `GeomAdaptor_Curve::Value(FirstParameter())` / `(LastParameter())` of
        // the edge as stored.
        let (first1, last1) = BRepTool::edge_parameters(edge1);
        let (first2, last2) = BRepTool::edge_parameters(edge2);
        let p11 = curve1.d0(first1);
        let p12 = curve1.d0(last1);
        let p21 = curve2.d0(first2);
        let p22 = curve2.d0(last2);
        p11.square_distance(&p21) < SQUARE_CONFUSION
            || p11.square_distance(&p22) < SQUARE_CONFUSION
            || p12.square_distance(&p21) < SQUARE_CONFUSION
            || p12.square_distance(&p22) < SQUARE_CONFUSION
    }
}
