//! `ChFi2d_FilletAlgo` — the iterative fillet between two planar edges.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi2d/ChFi2d_FilletAlgo.hxx:28-231` /
//! `.cxx:1-977` (OCCT 8.0.0).
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `ChFi2d_FilletAlgo()` (`cxx:36-45`) | 10 | [`FilletAlgo::new`] |
//! | `ChFi2d_FilletAlgo(wire, plane)` (`cxx:47-59`) | 13 | [`FilletAlgo::from_wire`] |
//! | `ChFi2d_FilletAlgo(e1, e2, plane)` (`cxx:61-75`) | 15 | [`FilletAlgo::from_edges`] |
//! | `Init(wire, plane)` (`cxx:77-102`) | 26 | [`FilletAlgo::init_wire`] |
//! | `Init(e1, e2, plane)` (`cxx:104-161`) | 58 | [`FilletAlgo::init_edges`] |
//! | `Perform(radius)` (`cxx:350-412`) | 63 | [`FilletAlgo::perform`] |
//! | `NbResults(point)` (`cxx:565-590`) | 26 | [`FilletAlgo::nb_results`] |
//! | `Result(...)` (`cxx:592-759`) | 168 | [`FilletAlgo::result`] |
//!
//! `FillPoint` / `FillDiff` / `PerformNewton` / `ProcessPoint` live in
//! [`fill`], `FilletPoint` in [`fillet_point`], `IsRadiusIntersected` in
//! [`radius_intersected`] and the `Result` body in [`make_result`].
//!
//! Translation notes:
//! - `BRepAdaptor_Curve(E)` is [`occt_geom::projlib::curve_type`] over the
//!   edge's stored 3D curve, so `GetType() < GetType()` (`cxx:115`) is the
//!   numeric order of `GeomAbs_CurveType` (which
//!   [`occt_core::kernel::geomabs::CurveType`] now derives),
//! - `BRep_Tool::Curve(E, first, last)` is
//!   [`BRepTool::edge_curve`] + [`BRepTool::edge_parameters`]; the pair is
//!   self-consistent, so the parameters of the projected 2D curves, of
//!   `Result`'s `D1` calls and of the divided edges all live in the same space,
//! - `GeomProjLib::Curve2d(C, First, Last, S)` (`GeomProjLib.cxx:195-202`) is
//!   [`occt_geom::geom_proj_lib::curve_2d_range`] over a
//!   [`occt_geom::plane::GeomPlane`] with `Tol = Precision::PConfusion()`;
//!   OCCT's null-handle answer then crashes on `myCurve1->IsPeriodic()`
//!   (`cxx:135`), the port answers `Err`,
//! - `Geom2dAdaptor_Curve::GetType()` is
//!   [`occt_geom2d::geom2d_int::geom2d_curve_type`]; `Degree()` / `NbPoles()`
//!   are `bspline_degree()` / `poles2d().len()` (only read in the
//!   `GeomAbs_BSplineCurve` arm),
//! - `Standard_ConstructionError("The fillet algorithms expects a wire
//!   consisting of two edges.")` (`cxx:99-100`) is the `Err(String)` of the two
//!   `Init` methods, with the OCCT message text,
//! - `Result` answers one `Option<FilletAlgoResult>` instead of a null
//!   `TopoDS_Edge` plus two out-parameters: OCCT leaves both modified edges
//!   untouched when it returns the null edge (`cxx:645-648`).

mod fillet_point;
mod fill;
mod make_result;
mod radius_intersected;

use std::sync::Arc;

use occt_core::gp::{GpDir, GpPln};
use occt_core::precision::CONFUSION;
use occt_geom::geom_proj_lib;
use occt_geom::plane::GeomPlane;
use occt_geom::surface::Surface;
use occt_geom2d::geom2d_int::{geom2d_curve_type, GeomAbsCurveType2d};
use occt_geom2d::Curve2d;

use crate::brep_tool::BRepTool;
use crate::iterator::ShapeIterator;
use crate::shape::{Edge, Wire};

use fillet_point::FilletPoint;

/// `Result`'s three edges (`cxx:592-759`): the fillet arc and the two initial
/// edges cut at the points of tangency.
#[derive(Clone)]
pub struct FilletAlgoResult {
    /// `aResult` — the fillet arc (`cxx:650-688`).
    pub fillet: Edge,
    /// `theEdge1` (`cxx:715-756`), the first initial edge cut at the fillet
    /// start.
    pub edge1: Edge,
    /// `theEdge2` (`cxx:719-755`), the second initial edge cut at the fillet
    /// end.
    pub edge2: Edge,
}

/// `ChFi2d_FilletAlgo` (`hxx:56-137`).
#[derive(Clone, Default)]
pub struct FilletAlgo {
    /// `myEdge1`, `myEdge2` (`hxx:117`), possibly exchanged by `Init`.
    edge1: Option<Edge>,
    edge2: Option<Edge>,
    /// `myPlane` (`hxx:120`).
    plane: GpPln,
    /// `myCurve1`, `myCurve2` (`hxx:122`): the projected 2D curves.
    curve1: Option<Arc<dyn Curve2d>>,
    curve2: Option<Arc<dyn Curve2d>>,
    /// `myStart1, myEnd1, myStart2, myEnd2, myRadius` (`hxx:124`).
    start1: f64,
    end1: f64,
    start2: f64,
    end2: f64,
    radius: f64,
    /// `myResultParams` (`hxx:126`): parameters where roots were found.
    result_params: Vec<f64>,
    /// `myResultOrientation` (`hxx:128`): position of the fillet relative to the
    /// first curve, one flag per `myResultParams` entry.
    result_orientation: Vec<bool>,
    /// `myStartSide` (`hxx:130`).
    start_side: bool,
    /// `myEdgesExchnged` (`hxx:133`).
    edges_exchanged: bool,
    /// `myDegreeOfRecursion` (`hxx:135`).
    degree_of_recursion: i32,
}

impl FilletAlgo {
    /// `ChFi2d_FilletAlgo()` (`cxx:36-45`): every parameter at `0.0`, both flags
    /// `false`, the recursion degree `0`.
    pub fn new() -> Self {
        Self::default()
    }

    /// `ChFi2d_FilletAlgo(theWire, thePlane)` (`cxx:47-59`).
    pub fn from_wire(wire: &Wire, plane: &GpPln) -> Result<Self, String> {
        let mut algo = Self::new();
        algo.init_wire(wire, plane)?;
        Ok(algo)
    }

    /// `ChFi2d_FilletAlgo(theEdge1, theEdge2, thePlane)` (`cxx:61-75`): the
    /// member edges are the two arguments before `Init` runs.
    pub fn from_edges(edge1: &Edge, edge2: &Edge, plane: &GpPln) -> Result<Self, String> {
        let mut algo = Self::new();
        algo.edge1 = Some(edge1.clone());
        algo.edge2 = Some(edge2.clone());
        algo.init_edges(edge1, edge2, plane)?;
        Ok(algo)
    }

    /// `Init(theWire, thePlane)` (`cxx:77-102`): the first two wire children,
    /// then the two-edge `Init`.
    ///
    /// `TopoDS::Edge(itr.Value())` raises `Standard_TypeMismatch` on a non-edge
    /// child; that arm is the `Err`.
    pub fn init_wire(&mut self, wire: &Wire, plane: &GpPln) -> Result<(), String> {
        let mut edge1: Option<Edge> = None;
        let mut edge2: Option<Edge> = None;
        for child in ShapeIterator::of_shape(&wire.0) {
            let edge = Edge::wrap(child)
                .ok_or_else(|| "ChFi2d_FilletAlgo::Init: wire child is not an edge".to_string())?;
            if edge1.is_none() {
                edge1 = Some(edge);
            } else if edge2.is_none() {
                edge2 = Some(edge);
            } else {
                break;
            }
        }
        let (Some(edge1), Some(edge2)) = (edge1, edge2) else {
            return Err(
                "The fillet algorithms expects a wire consisting of two edges.".to_string(),
            );
        };
        self.init_edges(&edge1, &edge2, plane)
    }

    /// `Init(theEdge1, theEdge2, thePlane)` (`cxx:104-161`).
    pub fn init_edges(&mut self, edge1: &Edge, edge2: &Edge, plane: &GpPln) -> Result<(), String> {
        self.plane = plane.clone();
        self.edges_exchanged = false;

        // `BRepAdaptor_Curve aBAC1(theEdge1), aBAC2(theEdge2);` (`cxx:113-114`).
        let curve_a1 = BRepTool::edge_curve(edge1)
            .ok_or_else(|| "ChFi2d_FilletAlgo::Init: the edge has no 3D curve".to_string())?;
        let curve_a2 = BRepTool::edge_curve(edge2)
            .ok_or_else(|| "ChFi2d_FilletAlgo::Init: the edge has no 3D curve".to_string())?;
        let type1 = occt_geom::projlib::curve_type(curve_a1.as_ref());
        let type2 = occt_geom::projlib::curve_type(curve_a2.as_ref());

        // `if (aBAC1.GetType() < aBAC2.GetType())` (`cxx:115-125`): the first
        // curve must be the more complicated one.
        if type1 < type2 {
            self.edge1 = Some(edge2.clone());
            self.edge2 = Some(edge1.clone());
            self.edges_exchanged = true;
        } else {
            self.edge1 = Some(edge1.clone());
            self.edge2 = Some(edge2.clone());
        }

        // `BRep_Tool::Curve(myEdge1, myStart1, myEnd1)` (`cxx:127-128`).
        let a_curve1 = BRepTool::edge_curve(self.edge1.as_ref().unwrap())
            .ok_or_else(|| "ChFi2d_FilletAlgo::Init: the edge has no 3D curve".to_string())?;
        (self.start1, self.end1) = BRepTool::edge_parameters(self.edge1.as_ref().unwrap());
        let a_curve2 = BRepTool::edge_curve(self.edge2.as_ref().unwrap())
            .ok_or_else(|| "ChFi2d_FilletAlgo::Init: the edge has no 3D curve".to_string())?;
        (self.start2, self.end2) = BRepTool::edge_parameters(self.edge2.as_ref().unwrap());

        // `GeomProjLib::Curve2d(aCurve1, myStart1, myEnd1, myPlane)`
        // (`cxx:130-131`): the `(C, First, Last, S)` overload,
        // `Tol = Precision::PConfusion()`.
        let surface: Arc<dyn Surface> = Arc::new(GeomPlane::new(self.plane.clone()));
        let mut tolerance = occt_core::precision::PCONFUSION;
        self.curve1 = geom_proj_lib::curve_2d_range(
            &a_curve1,
            self.start1,
            self.end1,
            &surface,
            &mut tolerance,
        );
        let mut tolerance = occt_core::precision::PCONFUSION;
        self.curve2 = geom_proj_lib::curve_2d_range(
            &a_curve2,
            self.start2,
            self.end2,
            &surface,
            &mut tolerance,
        );
        let (Some(curve1), Some(curve2)) = (self.curve1.as_ref(), self.curve2.as_ref()) else {
            return Err("ChFi2d_FilletAlgo::Init: the projection is not a known curve type"
                .to_string());
        };

        // `while (myCurve1->IsPeriodic() && myStart1 >= myEnd1) myEnd1 +=
        // myCurve1->Period();` (`cxx:135-140`).
        while curve1.is_periodic() && self.start1 >= self.end1 {
            self.end1 += curve1.period();
        }
        while curve2.is_periodic() && self.start2 >= self.end2 {
            self.end2 += curve2.period();
        }

        // `if (aBAC1.GetType() == aBAC2.GetType())` (`cxx:141-160`): the first
        // curve must be paramerically shorter.
        if type1 == type2 && self.end2 - self.start2 < self.end1 - self.start1 {
            std::mem::swap(&mut self.edge1, &mut self.edge2);
            std::mem::swap(&mut self.curve1, &mut self.curve2);
            std::mem::swap(&mut self.start1, &mut self.start2);
            std::mem::swap(&mut self.end1, &mut self.end2);
            self.edges_exchanged = true;
        }
        Ok(())
    }

    /// `Perform(theRadius)` (`cxx:350-412`): two passes over the first curve (one
    /// per side of it), each walking the sample points and running
    /// [`FilletAlgo::perform_newton`] on every consecutive pair.
    ///
    /// `false` when `Init` never ran (OCCT dereferences a null
    /// `myCurve1` in `Geom2dAdaptor_Curve aGAC(myCurve1)`, `cxx:358`).
    pub fn perform(&mut self, radius: f64) -> bool {
        self.degree_of_recursion = 0;
        self.result_params.clear();
        self.result_orientation.clear();

        // `Geom2dAdaptor_Curve aGAC(myCurve1); switch (aGAC.GetType())`
        // (`cxx:358-377`).
        let Some(curve1) = self.curve1.clone() else {
            return false;
        };
        let a_nb_steps = match geom2d_curve_type(curve1.as_ref()) {
            GeomAbsCurveType2d::Line => 2.0,
            GeomAbsCurveType2d::Circle => 4.0,
            GeomAbsCurveType2d::Ellipse => 5.0,
            GeomAbsCurveType2d::BSplineCurve => {
                let degree = curve1.bspline_degree().unwrap_or(0) as f64;
                let nb_poles = curve1.poles2d().map_or(0, |p| p.len()) as f64;
                2.0 + degree * nb_poles
            }
            // unknown: maximum (`cxx:375-376`).
            _ => 100.0,
        };

        self.radius = radius;
        let a_step = (self.end1 - self.start1) / a_nb_steps;
        let a_d_step = 1.0e-4 * a_step;

        // `for (aCycle = 2, myStartSide = false; aCycle; myStartSide =
        // !myStartSide, aCycle--)` (`cxx:386-410`): two passes, `false` then
        // `true`, and the trailing increment leaves the member `false`.
        for a_cycle in 0..2 {
            self.start_side = a_cycle == 1;

            let mut a_left: Option<FilletPoint> = None;
            let mut a_param = self.start1 + a_step;
            // `aParam < myEnd1 || std::abs(myEnd1 - aParam) <
            // Precision::Confusion()`.
            while a_param < self.end1 || (self.end1 - a_param).abs() < CONFUSION {
                if a_left.is_none() {
                    let mut left = FilletPoint::new(a_param - a_step);
                    self.fill_point(&mut left, a_param);
                    self.fill_diff(&mut left, a_d_step, true);
                    a_left = Some(left);
                }

                let mut a_right = FilletPoint::new(a_param);
                self.fill_point(&mut a_right, a_param - a_step);
                self.fill_diff(&mut a_right, a_d_step, false);

                if let Some(left) = a_left.as_mut() {
                    left.filter_points(&a_right);
                }
                if let Some(left) = a_left.as_ref() {
                    self.perform_newton(left, &a_right);
                }
                a_left = Some(a_right);
                a_param += a_step;
            }
        }
        self.start_side = false;

        !self.result_params.is_empty()
    }

    /// `NbResults(thePoint)` (`cxx:565-590`): how many of the stored parameters
    /// are still roots once `FillPoint` re-runs on them.
    ///
    /// `thePoint` is unused in the OCCT body (it only projects it into
    /// `aTargetPoint2d`, which nothing reads); kept for the signature.
    pub fn nb_results(&mut self, _the_point: &occt_core::gp::GpPnt) -> i32 {
        let mut nb = 0;
        let mut a = 1usize;
        while a <= self.result_params.len() {
            // `myStartSide = (myResultOrientation.Value(i)) != 0;`
            self.start_side = self.result_orientation[a - 1];
            let param = self.result_params[a - 1];
            let mut point = FilletPoint::new(param);
            self.fill_point(&mut point, param + 1.0);
            if point.has_solution(self.radius).is_some() {
                nb += 1;
            }
            a += 1;
        }
        nb
    }

    /// `Result(thePoint, theEdge1, theEdge2, iSolution)` (`cxx:592-759`):
    /// the stored parameter whose fillet centre is nearest to `thePoint`
    /// (`iSolution == -1`) or the `iSolution`-th one, together with the two
    /// initial edges cut at the points of tangency.
    pub fn result(
        &mut self,
        the_point: &occt_core::gp::GpPnt,
        i_solution: i32,
    ) -> Option<FilletAlgoResult> {
        make_result::result(self, the_point, i_solution)
    }

    /// `myEdgesExchnged` (`hxx:133`): whether `Init` swapped the two edges, so
    /// `Result`'s two divided edges land in the opposite order.
    pub fn edges_exchanged(&self) -> bool {
        self.edges_exchanged
    }
}

/// `gp_Pln::Axis().Direction()` (`gp_Pln.hxx:79-82`): the normal of the plane's
/// position as a direction vector.
pub(super) fn plane_normal(plane: &GpPln) -> GpDir {
    *plane.axis().direction()
}
