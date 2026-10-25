//! `TopOpeBRepTool_CurveTool` — the slice `TopOpeBRep_GeomTool` needs.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepTool/TopOpeBRepTool_CurveTool.hxx:34-82`,
//! `.cxx:69-1156`.
//!
//! Ported: the two static degree-1 B-spline builders,
//! `MakeBSpline1fromPnt` (`cxx:898-934`) and `MakeBSpline1fromPnt2d`
//! (`cxx:937-972`), the projection gate `IsProjectable` (`cxx:975-1037`) and
//! `MakePCurveOnFace` (`cxx:1040-1156`) with its local `MakePCurve`
//! (`cxx:111-139`).
//!
//! The constructors and the `myGeomTool` accessors (`cxx:69-106`) hold the
//! `TopOpeBRepTool_GeomTool` that decides which output curves `MakeCurves`
//! computes.
//!
//! `MakeCurves` (`cxx:300-895`) and its file-local helpers
//! (`CheckApproxResults` `:141-163`, `CheckPCurve` `:166-227`,
//! `MakeCurve3DfromWLineApprox` `:230-251`,
//! `MakeCurve2DfromWLineApproxAndPlane` `:253-274`,
//! `MakeCurve2DfromWLineApprox` `:276-296`) live in the sibling file
//! [`make_curves`].
//!
//! `#define CurveImprovement` is active in OCCT 8.0.0 (`cxx:68`), so the
//! `nbpol > 100` pole-decimation block (`cxx:358-545`) is translated as well;
//! `#define IFV` is commented out (`cxx:67`) and every `IFV` block is dead.
//!
//! UNPORTED: nothing in this class.

mod make_curves;

pub use make_curves::MakeCurvesOutcome;

use std::sync::Arc;

use occt_core::gp::{GpAx2d, GpDir2d, GpPnt, GpPnt2d, GpTrsf2d, GpVec2d};
use occt_core::kernel::geomabs::{CurveType, SurfaceType};
use occt_geom::bspline_curve::GeomBSplineCurve;
use occt_geom::curve::Curve;
use occt_geom::projlib::{curve_type, surface_type, ProjectedCurve};
use occt_geom::surface::Surface;
use occt_geom::trimmed::GeomTrimmedCurveBasis;
use occt_geom2d::bspline_curve::Geom2dBSplineCurve;
use occt_geom2d::circle::Geom2dCircle;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::ellipse::Geom2dEllipse;
use occt_geom2d::hyperbola::Geom2dHyperbola;
use occt_geom2d::line::Geom2dLine;
use occt_geom2d::parabola::Geom2dParabola;

use super::geom_tool::GeomTool;
use super::out_curve_type::OutCurveType;

/// `class TopOpeBRepTool_CurveTool` (`hxx:34-82`).
///
/// The only state is `myGeomTool` (`hxx:78`); `MakeCurves` reads it and nothing
/// else of the instance.
#[derive(Debug, Clone, Copy)]
pub struct CurveTool {
    /// `myGeomTool` (`hxx:78`).
    geom_tool: GeomTool,
}

impl Default for CurveTool {
    /// `TopOpeBRepTool_CurveTool()` (`cxx:71`): the default-constructed
    /// `TopOpeBRepTool_GeomTool` (`GeomTool.cxx:26-35`).
    fn default() -> Self {
        Self::new()
    }
}

impl CurveTool {
    /// `TopOpeBRepTool_CurveTool()` (`cxx:71`).
    pub fn new() -> Self {
        Self { geom_tool: GeomTool::default() }
    }

    /// `TopOpeBRepTool_CurveTool(const TopOpeBRepTool_OutCurveType O)`
    /// (`cxx:75-79`): `TopOpeBRepTool_GeomTool GT(O)` keeps the `true` flag
    /// defaults (`GeomTool.cxx:26-35`).
    pub fn with_type(out_curve_type: OutCurveType) -> Self {
        Self { geom_tool: GeomTool::new(out_curve_type, true, true, true) }
    }

    /// `TopOpeBRepTool_CurveTool(const TopOpeBRepTool_GeomTool& GT)`
    /// (`cxx:83-86`) through `SetGeomTool` (`cxx:99-106`).
    pub fn with_geom_tool(geom_tool: GeomTool) -> Self {
        Self { geom_tool }
    }

    /// `ChangeGeomTool()` (`cxx:89-92`).
    pub fn change_geom_tool(&mut self) -> &mut GeomTool {
        &mut self.geom_tool
    }

    /// `GetGeomTool()` (`cxx:95-98`).
    pub fn geom_tool(&self) -> &GeomTool {
        &self.geom_tool
    }

    /// `SetGeomTool(GT)` (`cxx:99-101`): `myGeomTool.Define(GT)`.
    pub fn set_geom_tool(&mut self, gt: &GeomTool) {
        self.geom_tool.define_copy(gt);
    }

    /// `MakeBSpline1fromPnt(P)` (`cxx:898-934`): degree 1 over `P`, with the
    /// point index as parameter (`knots(i) = i`, `cxx:927-931`),
    /// multiplicity 1 at the interior knots and `Degree + 1 = 2` at both ends
    /// (`cxx:917-920`).
    ///
    /// `None` mirrors the OCCT exception for a pole count that cannot form a
    /// degree-1 B-spline (fewer than 2 poles leaves `nbknots < 2`, where the
    /// C++ end-multiplicity writes would collapse).
    pub fn make_bspline1_from_pnt(points: &[GpPnt]) -> Option<Arc<dyn Curve>> {
        let (knots, mults) = degree1_knots(points.len())?;
        let c = GeomBSplineCurve::from_poles_knots_mults(points.to_vec(), knots, mults, 1).ok()?;
        Some(Arc::new(c))
    }

    /// `MakeBSpline1fromPnt2d(P)` (`cxx:937-972`): the 2D form of
    /// [`Self::make_bspline1_from_pnt`], same knots and multiplicities.
    pub fn make_bspline1_from_pnt2d(points: &[GpPnt2d]) -> Option<Arc<dyn Curve2d>> {
        let (knots, mults) = degree1_knots(points.len())?;
        let flat = occt_core::bspl::banded_interp::knot_sequence(&knots, &mults, 1);
        let xs: Vec<f64> = points.iter().map(|p| p.x()).collect();
        let ys: Vec<f64> = points.iter().map(|p| p.y()).collect();
        let c = Geom2dBSplineCurve::new(xs, ys, flat, 1).ok()?;
        Some(Arc::new(c))
    }

    /// `IsProjectable(S, C3D)` (`cxx:975-1037`): refuse the projection of a
    /// quadric (ellipse / hyperbola / parabola) on a cone, of an ellipse on a
    /// cylinder, and of a circle on a sphere or a torus. OCCT states the
    /// reason at `cxx:986-994`: the projection would only build the part of
    /// the 2D curve inside the cone's UV domain.
    pub fn is_projectable(surface: &dyn Surface, c3d: &dyn Curve) -> bool {
        let suty = surface_type(surface);
        let cuty = curve_type(c3d);
        match suty {
            SurfaceType::Cone => !matches!(
                cuty,
                CurveType::Ellipse | CurveType::Hyperbola | CurveType::Parabola
            ),
            SurfaceType::Cylinder => cuty != CurveType::Ellipse,
            SurfaceType::Sphere | SurfaceType::Torus => cuty != CurveType::Circle,
            _ => true,
        }
    }

    /// `MakePCurveOnFace(S, C3D, TolReached2d, first, last)` (`cxx:1040-1156`):
    /// run `ProjLib_ProjectedCurve` over `[first, last]` (the whole curve when
    /// `first >= last`, `cxx:1046-1052`), wrap the analytic result with
    /// `MakePCurve` (`cxx:111-139`), then move the pcurve into the face's UV
    /// bounds: the sphere's quasi-periodic V mirror (`cxx:1084-1116`) and the
    /// periodic `du` / `dv` shift (`cxx:1118-1153`).
    ///
    /// `uv_bounds` is `BRepTools::UVBounds(F, UMin, UMax, VMin, VMax)`
    /// (`cxx:1073-1074`); `None` stands for a surface without a face and then
    /// falls back to the surface's own range. The sphere mirror and the
    /// periodic window do **not** read the same values in OCCT: the mirror uses
    /// `BAS.FirstVParameter()` / `LastVParameter()` and `FirstUParameter()`
    /// (`cxx:1085-1106`), and `BAS` is built with `Restriction = false`
    /// (`cxx:1055`), i.e. those are the natural surface parameters, while the
    /// `du` / `dv` window uses the face's `UVBounds`.
    ///
    /// `None` is the port's answer where OCCT throws: `MakePCurve` raises
    /// `Standard_NotImplemented` for a curve type outside
    /// line / circle / ellipse / parabola / hyperbola / B-spline, and
    /// `ProjLib_ProjectedCurve` can leave `myResult` undone on the
    /// approximation arms (see `projlib::ProjectedCurve`).
    pub fn make_pcurve_on_face(
        surface: Arc<dyn Surface>,
        c3d: Arc<dyn Curve>,
        first: f64,
        last: f64,
        uv_bounds: Option<(f64, f64, f64, f64)>,
    ) -> (Option<Arc<dyn Curve2d>>, f64) {
        let trim = first < last;
        let gac: Arc<dyn Curve> = if trim {
            Arc::new(GeomTrimmedCurveBasis::new(c3d, first, last))
        } else {
            c3d
        };

        let mut projcurv = ProjectedCurve::with_surface(surface.clone());
        projcurv.perform(gac.clone());
        let tol_reached_2d = projcurv.get_tolerance();
        let Some(mut c2d) = make_pcurve(&projcurv) else {
            return (None, tol_reached_2d);
        };

        // `BAS.FirstUParameter()` / `LastVParameter()` with `Restriction = false`
        // (`cxx:1055`, `1085-1086`, `1105`) — the natural surface range.
        let (u_first, u_last) = surface.u_range();
        let (v_first, v_last) = surface.v_range();
        // `BRepTools::UVBounds(F, UMin, UMax, VMin, VMax)` (`cxx:1073-1074`).
        let (umin, umax, vmin, vmax) = uv_bounds.unwrap_or((u_first, u_last, v_first, v_last));

        let f = gac.first_parameter();
        let l = gac.last_parameter();
        let t = (f + l) * 0.5;
        let p_c2d = c2d.d0(t);
        let mut u2 = p_c2d.x();
        let mut v2 = p_c2d.y();

        if surface_type(surface.as_ref()) == SurfaceType::Sphere {
            // MSV: consider the quasiperiodic shift of the pcurve
            // (`cxx:1084-1116`).
            let mincond = v2 < v_first;
            let maxcond = v2 > v_last;
            if mincond || maxcond {
                // Mirror relative to the isoline of the apex -PI/2 or PI/2.
                let mut po = GpPnt2d::new(0.0, -std::f64::consts::FRAC_PI_2);
                if maxcond {
                    po.set_y(std::f64::consts::FRAC_PI_2);
                }
                let dir = GpDir2d::new(1.0, 0.0).expect("gp_Dir2d X");
                let mut a_trsf = GpTrsf2d::identity();
                a_trsf.set_mirror_ax2d(&GpAx2d::new(po, dir));
                c2d = Arc::from(c2d.transformed(&a_trsf));

                // Add a translation along U of PI, towards the apex side.
                let mut vec = GpVec2d::new(std::f64::consts::PI, 0.0);
                if u2 - u_first - std::f64::consts::PI > -1e-7 {
                    vec.reverse();
                }
                let mut a_tr = GpTrsf2d::identity();
                a_tr.set_translation_vec(&vec);
                c2d = Arc::from(c2d.transformed(&a_tr));

                // Recompute the test point.
                let p = c2d.d0(t);
                u2 = p.x();
                v2 = p.y();
            }
        }

        let mut du = 0.0;
        if surface.is_u_periodic() {
            let mincond = umin - u2 > 1e-7;
            let maxcond = u2 - umax > 1e-7;
            if mincond || maxcond {
                du = if mincond {
                    surface.u_period()
                } else {
                    -surface.u_period()
                };
            }
        }
        let mut dv = 0.0;
        if surface.is_v_periodic() {
            let mincond = vmin - v2 > 1e-7;
            let maxcond = v2 - vmax > 1e-7;
            if mincond || maxcond {
                dv = if mincond {
                    surface.v_period()
                } else {
                    -surface.v_period()
                };
            }
        }

        if du != 0.0 || dv != 0.0 {
            let mut a_tr = GpTrsf2d::identity();
            a_tr.set_translation_vec(&GpVec2d::new(du, dv));
            c2d = Arc::from(c2d.transformed(&a_tr));
        }

        (Some(c2d), tol_reached_2d)
    }
}

/// The local `MakePCurve(const ProjLib_ProjectedCurve&)` (`cxx:111-139`): hand
/// the projector's analytic result to the matching `Geom2d_*` class.
fn make_pcurve(pc: &ProjectedCurve) -> Option<Arc<dyn Curve2d>> {
    let st = &pc.result;
    match pc.get_type() {
        CurveType::Line => Some(Arc::new(Geom2dLine::new(st.lin.pos))),
        CurveType::Circle => Some(Arc::new(Geom2dCircle::new(st.circ))),
        CurveType::Ellipse => Some(Arc::new(Geom2dEllipse::new(st.elips))),
        CurveType::Parabola => Some(Arc::new(Geom2dParabola::new(st.parab))),
        CurveType::Hyperbola => Some(Arc::new(Geom2dHyperbola::new(st.hypr))),
        CurveType::BSplineCurve => st.bspline.clone().map(|b| b as Arc<dyn Curve2d>),
        _ => None,
    }
}

/// The `knots` / `mults` arrays of `cxx:917-920` and `cxx:956-959`:
/// `nbknots = nbpoints - Degree + 1` distinct knots at the point indexes, with
/// multiplicity `Degree + 1` at both ends and 1 elsewhere.
fn degree1_knots(nbpoints: usize) -> Option<(Vec<f64>, Vec<i32>)> {
    const DEGREE: usize = 1;
    if nbpoints < DEGREE + 1 {
        return None;
    }
    let nbknots = nbpoints - DEGREE + 1;
    let knots: Vec<f64> = (1..=nbknots).map(|i| i as f64).collect();
    let mut mults = vec![1i32; nbknots];
    mults[0] = DEGREE as i32 + 1;
    mults[nbknots - 1] = DEGREE as i32 + 1;
    Some((knots, mults))
}
