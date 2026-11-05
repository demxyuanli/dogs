//! `BRepFill_TrimSurfaceTool` — the p-curves and 3D curves of the trimming of
//! a face pair by an extruded surface.
//!
//! Source: `ModelingAlgorithms/TKBool/BRepFill/BRepFill_TrimSurfaceTool.cxx:17-486`
//! and `.hxx:30-83` (OCCT 8.0.0).
//!
//! The tool holds the two faces, their two iso edges (the `Umin` / `Vmin` edges
//! of each face), the two inversion flags and the 2D bisector between the
//! faces. `BRepFill_Sweep` uses it in three ways:
//!
//! * [`TrimSurfaceTool::intersect_with`] — `IntersectWith(EdgeOnF1, EdgeOnF2,
//!   Points)` (`cxx:382-404`): project both edges in the `xOy` plane, intersect
//!   each projection with the bisector and answer the (bisector, edge1, edge2)
//!   parameter triples.
//! * [`TrimSurfaceTool::is_on_face`] — `IsOnFace(Point)` (`cxx:408-427`): does
//!   the line `(Point, gp::DZ())` hit either face?
//! * [`TrimSurfaceTool::project`] — `Project(U1, U2, ...)` (`cxx:462-487`): the
//!   analytic curves of the particular case, or the approximated ones from
//!   [`ApproxSeewing`](crate::brep_fill_approx_seewing::ApproxSeewing).
//!
//! Ported layout:
//! * this file — the type, [`Error`], `IntersectWith` / `IsOnFace` / `ProjOn` /
//!   `Project`,
//! * [`eval_parameters`] — the static `Bubble` / `EvalPhase` / `EvalParameters`
//!   helpers of `cxx:100-378`.
//!
//! UNPORTED / deviations (all of them the null-handle / exception
//! correspondence, none of them a missing branch):
//! * `cxx:219-223` mutates the caller's `myBis` in place through the shared
//!   `Geom2d_TrimmedCurve` handle. The port clones the trimmed curve, runs the
//!   same `SetTrim`, and re-loads it into the adaptor; because
//!   `Geom2dAdaptor_Curve::load` unwraps a trimmed curve to its basis
//!   (`Geom2dAdaptor_Curve.cxx:285-288`), both forms feed the same basis and
//!   the same `myFirst` / `myLast` to the intersector, so only the observed
//!   growth of the caller's handle is lost. `occ::down_cast` on a curve that
//!   is not a `Geom2d_TrimmedCurve` answers a null handle in OCCT, whose
//!   `SetTrim` would raise `Standard_NullObject`; the port skips the
//!   elongation there instead of raising.
//! * `Standard_NotImplemented` (`cxx:362`), `StdFail_NotDone` (`cxx:194`,
//!   `cxx:230`, `cxx:392`) and the null 3D curve of `BRep_Tool::Curve`
//!   (`cxx:184`, `cxx:436`) become [`Error`] values.

use std::sync::Arc;

use occt_core::gp::{GpAx3, GpDir, GpLin, GpPln, GpPnt, GpPnt2d};
use occt_core::kernel::geomabs::Shape;
use occt_core::precision;
use occt_geom::curve::Curve;
use occt_geom::geom_proj_lib;
use occt_geom::plane::GeomPlane;
use occt_geom::trimmed::GeomTrimmedCurveBasis;
use occt_geom::Surface;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::geom2d_api_project_point_on_curve::ProjectPointOnCurve2d;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use crate::brep_fill_approx_seewing::ApproxSeewing;
use crate::brep_fill_multi_line::{Error as MultiLineError, MultiLine};
use crate::brep_int_curve_surface_inter::BRepIntCurveSurfaceInter;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face};

mod eval_parameters;

pub use eval_parameters::eval_parameters;

/// Failure modes of [`TrimSurfaceTool`].
#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    /// `BRep_Tool::Curve` / `BRep_Tool::Surface` answered null (`cxx:184`,
    /// `cxx:321`, `cxx:436`).
    NullCurve,
    /// `GeomProjLib::Curve2d` answered a null handle (`cxx:188`, `cxx:442`).
    NullProjection,
    /// `Geom_TrimmedCurve` refused the edge range (`cxx:185`).
    TrimError,
    /// `StdFail_NotDone("BRepFill_TrimSurfaceTool::IntersectWith")`
    /// (`cxx:194`, `cxx:230`).
    IntersectionNotDone,
    /// `StdFail_NotDone` of the `Points.Length() != Points2.Lenth()` guard
    /// (`cxx:392-393`).
    IncoherentIntersection,
    /// `Standard_NotImplemented(" BRepFill_TrimSurfaceTool")` for a surface
    /// type outside the five analytic branches (`cxx:361-362`).
    NotImplemented,
    /// `Standard_ConstructionError` of `gp_Ax3(P, DZ, gp::DZ())` inside
    /// `GeomAdaptor_SurfaceOfRevolution::Load` (`cxx:186-188` of the adaptor).
    ConstructionError,
    /// The `BRepFill_MultiLine` / `BRepFill_ApproxSeewing` side of `Project`.
    MultiLine(MultiLineError),
    /// The `BRepFill_ApproxSeewing` side of `Project`.
    Approx(crate::brep_fill_approx_seewing::Error),
}

impl From<MultiLineError> for Error {
    fn from(e: MultiLineError) -> Self {
        Error::MultiLine(e)
    }
}

impl From<crate::brep_fill_approx_seewing::Error> for Error {
    fn from(e: crate::brep_fill_approx_seewing::Error) -> Self {
        Error::Approx(e)
    }
}

/// `BRepFill_TrimSurfaceTool` (`hxx:36-83`).
#[derive(Clone)]
pub struct TrimSurfaceTool {
    /// `myFace1` (`hxx:76`).
    face1: Face,
    /// `myFace2` (`hxx:77`).
    face2: Face,
    /// `myEdge1` (`hxx:78`).
    edge1: Edge,
    /// `myEdge2` (`hxx:79`).
    edge2: Edge,
    /// `myInv1` (`hxx:80`).
    inv1: bool,
    /// `myInv2` (`hxx:81`).
    inv2: bool,
    /// `myBis` (`hxx:82`).
    bis: Arc<dyn Curve2d>,
}

impl TrimSurfaceTool {
    /// `BRepFill_TrimSurfaceTool(Bis, Face1, Face2, Edge1, Edge2, Inv1, Inv2)`
    /// (`cxx:74-95`): a plain member initialisation.
    pub fn new(
        bis: Arc<dyn Curve2d>,
        face1: Face,
        face2: Face,
        edge1: Edge,
        edge2: Edge,
        inv1: bool,
        inv2: bool,
    ) -> Self {
        Self {
            face1,
            face2,
            edge1,
            edge2,
            inv1,
            inv2,
            bis,
        }
    }

    /// `IntersectWith(EdgeOnF1, EdgeOnF2, Points)` (`cxx:382-404`).
    ///
    /// `Points` answers the `gp_Pnt(BisParameter, EdgeOnF1Parameter,
    /// EdgeOnF2Parameter)` triples, ordered by increasing bisector parameter
    /// (the `Bubble` of `cxx:271`), with the duplicates of `cxx:275-286`
    /// removed.
    pub fn intersect_with(
        &self,
        edge_on_f1: &Edge,
        edge_on_f2: &Edge,
    ) -> Result<Vec<GpPnt>, Error> {
        // `Points.Clear();` (`cxx:385`) — answered by the fresh `Vec`.
        let mut points = Vec::new();
        let mut points2 = Vec::new();

        eval_parameters(edge_on_f1, &self.face1, &self.bis, &mut points)?;
        eval_parameters(edge_on_f2, &self.face2, &self.bis, &mut points2)?;

        // `StdFail_NotDone_Raise_if(Points.Length() != Points2.Length(), ...)`
        // (`cxx:392-394`).
        if points.len() != points2.len() {
            return Err(Error::IncoherentIntersection);
        }

        // `cxx:396-403`: the Z coordinate is the parameter on the second edge.
        for (i, p) in points.iter_mut().enumerate() {
            p.set_z(points2[i].y());
        }
        Ok(points)
    }

    /// `IsOnFace(Point)` (`cxx:408-427`): `BRepIntCurveSurface_Inter` on the
    /// line `(Point.X(), Point.Y(), 0) + t * gp::DZ()` against `myFace1`, then
    /// `myFace2`, both with a tolerance of `1e-6`.
    pub fn is_on_face(&self, point: &GpPnt2d) -> bool {
        let p = GpPnt::new(point.x(), point.y(), 0.0);
        let line = GpLin::from_pnt_dir(p, GpDir::new(0.0, 0.0, 1.0).expect("gp::DZ"));

        let mut inter = BRepIntCurveSurfaceInter::new();
        inter.init_shape_line(&self.face1.0, &line, 1.0e-6);
        if inter.more() {
            return true;
        }
        inter.init_shape_line(&self.face2.0, &line, 1.0e-6);
        inter.more()
    }

    /// `ProjOn(Point, Edge)` (`cxx:432-458`): the parameter of `Point` on the
    /// edge's 3D curve projected in the `xOy` plane.
    pub fn proj_on(&self, point: &GpPnt2d, edge: &Edge) -> Result<f64, Error> {
        let plane: Arc<dyn Surface> = Arc::new(GeomPlane::new(GpPln::new(GpAx3::standard())));
        let c2d = self.project_edge_on_xy(edge, &plane)?;
        // `Geom2dAPI_ProjectPointOnCurve Projector(Point, C2d)` (`cxx:445`).
        let projector = ProjectPointOnCurve2d::with_curve(point, &c2d);
        projector.lower_distance_parameter().ok_or(Error::NullCurve)
    }

    /// `Project(U1, U2, Curve, PCurve1, PCurve2, Cont)` (`cxx:462-487`).
    ///
    /// Answers the 3D curve, the two p-curves and the continuity between the
    /// two faces along it. `BRepFill_MultiLine::Curves` supplies the analytic
    /// triple in the particular case, `BRepFill_ApproxSeewing` the approximated
    /// one otherwise.
    pub fn project(
        &self,
        u1: f64,
        u2: f64,
    ) -> Result<(Arc<dyn Curve>, Arc<dyn Curve2d>, Arc<dyn Curve2d>, Shape), Error> {
        // `CT = new Geom2d_TrimmedCurve(myBis, U1, U2)` (`cxx:469`).
        let ct: Arc<dyn Curve2d> = Arc::new(Geom2dTrimmedCurve::new(self.bis.clone(), u1, u2));
        let ml = MultiLine::from_faces(
            self.face1.clone(),
            self.face2.clone(),
            &self.edge1,
            &self.edge2,
            self.inv1,
            self.inv2,
            ct,
        )?;

        // `Cont = ML.Continuity();` (`cxx:473`).
        let cont = ml.continuity();

        if ml.is_particular_case() {
            // `cxx:475-478`.
            let (curve, p1, p2) = ml.curves()?;
            Ok((curve, p1, p2, cont))
        } else {
            // `cxx:481-486`.
            let seew = ApproxSeewing::from_multi_line(&ml)?;
            Ok((
                seew.curve()?.clone(),
                seew.curve_on_f1()?.clone(),
                seew.curve_on_f2()?.clone(),
                cont,
            ))
        }
    }

    /// `BRep_Tool::Curve(Edge, L, f, l)` + `Geom_TrimmedCurve(C, f, l)` +
    /// `Transform(L.Transformation())` + `GeomProjLib::Curve2d(CT, Plane)`
    /// (`cxx:184-188`, `cxx:436-442`).
    fn project_edge_on_xy(
        &self,
        edge: &Edge,
        plane: &Arc<dyn Surface>,
    ) -> Result<Arc<dyn Curve2d>, Error> {
        let curve = BRepTool::edge_curve_world(edge).ok_or(Error::NullCurve)?;
        let (f, l) = BRepTool::edge_parameters(edge);
        let ct: Arc<dyn Curve> = Arc::new(GeomTrimmedCurveBasis::new(curve, f, l));
        geom_proj_lib::curve_2d(&ct, plane).ok_or(Error::NullProjection)
    }
}

/// `Precision::Confusion()`-sized helpers shared with the static functions of
/// `cxx:100-378`: `1.e-6` is the `Tol` of `cxx:178` / `cxx:417` / `cxx:424`.
pub(crate) const TOL: f64 = 1.0e-6;

/// `M_PI`, the phase of `EvalPhase` (`cxx:147`).
pub(crate) const PI: f64 = std::f64::consts::PI;

/// `Precision::IsPositiveInfinite` (`cxx:306-307`).
pub(crate) fn is_positive_infinite(v: f64) -> bool {
    precision::Precision::is_positive_infinite(v.abs())
}
