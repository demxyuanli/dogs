//! `Contap_HContTool` - the sample grid and 2D helpers `Contap_SurfFunction`,
//! `Contap_ArcFunction` and `Contap_Contour` read the restriction domains with.
//!
//! Sources: `Contap_HContTool.cxx:26-304`.
//!
//! Deviation (file statics -> value): `NbSamplePoints` (`cxx:117-179`) stores
//! the four clamped bounds in the file-static `uinf`/`vinf`/`usup`/`vsup`
//! (`cxx:23`) that `SamplePoint` (`cxx:181-228`) reads back. The port returns
//! them as a [`SampleGrid`]; the two call sites (`Contap_SurfFunction.cxx:51-67`,
//! `Contap_ArcFunction.cxx:36-51`) run the same "count, then sample 1..=count"
//! sequence, so the grid is threaded through instead of being global.
//!
//! Not ported: `HasBeenSeen` / `NbPoints` / `Value` / `IsVertex` / `Vertex` /
//! `NbSegments` / `HasFirstPoint` / `HasLastPoint` / `IsAllSolution`
//! (`cxx:320-370`). They are the `Adaptor2d_Curve2d` polygon queries, all of
//! which answer a constant or raise `Standard_OutOfRange`, and no `Contap` code
//! calls them; the port's `Curve2d` trait carries the geometry directly.

use occt_core::gp::{GpPnt2d, GpVec2d};
use occt_core::precision::{INFINITE, REAL_EPSILON};
use occt_geom2d::extrema2d::EPcOfExtPc2d;
use occt_geom2d::Curve2d;

use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_top_adaptor_hvertex::HVertex;
use crate::geom_abs_curve_type::{curve2d_kind, GeomAbsCurveType};
use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// `Contap_HContTool::NbSamplesV(S, ...)` (`cxx:26-67`).
pub fn nb_samples_v(s: &SurfaceAdaptor) -> i32 {
    match s.get_type() {
        GeomAbsSurfaceType::Plane => 2,
        // `3 + S->NbVPoles()`.
        GeomAbsSurfaceType::BezierSurface => 3 + s.nb_v_poles().unwrap_or(0),
        GeomAbsSurfaceType::BSplineSurface => {
            // `nbs = S->NbVKnots(); nbs *= S->VDegree(); if (nbs < 2) nbs = 2;`
            let mut nbs = s.nb_v_knots().unwrap_or(0) * s.v_degree().unwrap_or(0);
            if nbs < 2 {
                nbs = 2;
            }
            nbs
        }
        GeomAbsSurfaceType::Cylinder
        | GeomAbsSurfaceType::Cone
        | GeomAbsSurfaceType::Sphere
        | GeomAbsSurfaceType::Torus
        | GeomAbsSurfaceType::SurfaceOfRevolution
        | GeomAbsSurfaceType::SurfaceOfExtrusion => 15,
        _ => 10,
    }
}

/// `Contap_HContTool::NbSamplesU(S, ...)` (`cxx:70-115`).
pub fn nb_samples_u(s: &SurfaceAdaptor) -> i32 {
    match s.get_type() {
        GeomAbsSurfaceType::Plane => 2,
        // `3 + S->NbUPoles()`.
        GeomAbsSurfaceType::BezierSurface => 3 + s.nb_u_poles().unwrap_or(0),
        GeomAbsSurfaceType::BSplineSurface => {
            // `nbs = S->NbUKnots(); nbs *= S->UDegree(); if (nbs < 2) nbs = 2;`
            let mut nbs = s.nb_u_knots().unwrap_or(0) * s.u_degree().unwrap_or(0);
            if nbs < 2 {
                nbs = 2;
            }
            nbs
        }
        GeomAbsSurfaceType::Torus => 20,
        GeomAbsSurfaceType::Cylinder
        | GeomAbsSurfaceType::Cone
        | GeomAbsSurfaceType::Sphere
        | GeomAbsSurfaceType::SurfaceOfRevolution
        | GeomAbsSurfaceType::SurfaceOfExtrusion => 10,
        _ => 10,
    }
}

/// The four bounds `NbSamplePoints` leaves behind for `SamplePoint`
/// (`Contap_HContTool.cxx:23`): `uinf`, `vinf`, `usup`, `vsup` after the
/// `RealFirst()` / `RealLast()` widening.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SampleGrid {
    uinf: f64,
    vinf: f64,
    usup: f64,
    vsup: f64,
}

impl SampleGrid {
    /// `Contap_HContTool::NbSamplePoints(S)` (`cxx:117-179`): stores the four
    /// bounds and returns the number of sample points.
    pub fn new(s: &SurfaceAdaptor) -> (Self, i32) {
        let mut uinf = s.first_u_parameter();
        let mut usup = s.last_u_parameter();
        let mut vinf = s.first_v_parameter();
        let mut vsup = s.last_v_parameter();

        if usup < uinf {
            std::mem::swap(&mut uinf, &mut usup);
        }
        if vsup < vinf {
            std::mem::swap(&mut vinf, &mut vsup);
        }
        if uinf == -INFINITE && usup == INFINITE {
            uinf = -1.0e5;
            usup = 1.0e5;
        } else if uinf == -INFINITE {
            uinf = usup - 2.0e5;
        } else if usup == INFINITE {
            usup = uinf + 2.0e5;
        }
        if vinf == -INFINITE && vsup == INFINITE {
            vinf = -1.0e5;
            vsup = 1.0e5;
        } else if vinf == -INFINITE {
            vinf = vsup - 2.0e5;
        } else if vsup == INFINITE {
            vsup = vinf + 2.0e5;
        }

        let grid = Self {
            uinf,
            vinf,
            usup,
            vsup,
        };
        let nbs = if s.get_type() == GeomAbsSurfaceType::BSplineSurface {
            // `int m = (NbSamplesU / 3) * (NbSamplesV / 3); if (m > 5) return m;`
            let m = (nb_samples_u(s) / 3) * (nb_samples_v(s) / 3);
            if m > 5 {
                m
            } else {
                5
            }
        } else {
            5
        };
        (grid, nbs)
    }

    /// `Contap_HContTool::SamplePoint(S, Index, U, V)` (`cxx:181-228`),
    /// `Index` 1-based.
    pub fn sample_point(&self, s: &SurfaceAdaptor, index: i32) -> (f64, f64) {
        if s.get_type() == GeomAbsSurfaceType::BSplineSurface {
            let nb_int_u = nb_samples_u(s) / 3;
            let nb_int_v = nb_samples_v(s) / 3;
            if nb_int_u * nb_int_v > 5 {
                let ind_u = (index - 1) / nb_int_u; // 0 -> nbIntV
                let ind_v = (index - 1) - ind_u * nb_int_u; // 0 -> nbIntU
                let u = self.uinf
                    + ((self.usup - self.uinf) / f64::from(nb_int_u + 1)) * f64::from(ind_u + 1);
                let v = self.vinf
                    + ((self.vsup - self.vinf) / f64::from(nb_int_v + 2)) * f64::from(ind_v + 1);
                return (u, v);
            }
        }

        let (u, v) = match index {
            1 => (
                0.75 * self.uinf + 0.25 * self.usup,
                0.75 * self.vinf + 0.25 * self.vsup,
            ),
            2 => (
                0.75 * self.uinf + 0.25 * self.usup,
                0.25 * self.vinf + 0.75 * self.vsup,
            ),
            3 => (
                0.25 * self.uinf + 0.75 * self.usup,
                0.75 * self.vinf + 0.25 * self.vsup,
            ),
            4 => (
                0.25 * self.uinf + 0.75 * self.usup,
                0.25 * self.vinf + 0.75 * self.vsup,
            ),
            _ => (
                0.5 * (self.uinf + self.usup),
                0.5 * (self.vinf + self.vsup),
            ),
        };
        (u, v)
    }
}

/// `Contap_HContTool::NbSamplesOnArc(A)` (`cxx:230-259`).
pub fn nb_samples_on_arc(a: &dyn Curve2d) -> i32 {
    let nbs_on_c: f64 = match curve2d_kind(a) {
        GeomAbsCurveType::Line => 2.0,
        GeomAbsCurveType::Circle
        | GeomAbsCurveType::Ellipse
        | GeomAbsCurveType::Hyperbola
        | GeomAbsCurveType::Parabola => 10.0,
        GeomAbsCurveType::BezierCurve => a.bezier_nb_poles().unwrap_or(0) as f64,
        // `nbsOnC = 2 + A->NbKnots() * A->Degree();`
        GeomAbsCurveType::BSplineCurve => {
            2.0 + (a.bspline_nb_knots().unwrap_or(0) * a.bspline_degree().unwrap_or(0)) as f64
        }
        _ => 10.0,
    };
    nbs_on_c as i32
}

/// `Contap_HContTool::Bounds(A, Ufirst, Ulast)` (`cxx:261-267`).
pub fn bounds(a: &dyn Curve2d) -> (f64, f64) {
    (a.first_parameter(), a.last_parameter())
}

/// `Contap_HContTool::Project(C, P, Paramproj, Ptproj)` (`cxx:269-304`):
/// `Extrema_EPCOfExtPC2d(P, C, 20, 1e-8, 1e-5)` and the minimum-square-distance
/// solution. `None` when the search is not done or found nothing.
pub fn project(c: &dyn Curve2d, p: &GpPnt2d) -> Option<(f64, GpPnt2d)> {
    // `double epsX = 1.e-8; int Nbu = 20; double Tol = 1.e-5;` (`cxx:275-278`).
    let extrema = EPcOfExtPc2d::new(p, c, 20, 1.0e-8, 1.0e-5);
    if !extrema.is_done() {
        return None;
    }
    let nbext = extrema.nb_ext();
    if nbext == 0 {
        return None;
    }
    let mut index_min = 1usize;
    let mut dist2 = extrema.square_distance(1);
    for i in 2..=nbext {
        if extrema.square_distance(i) < dist2 {
            index_min = i;
            dist2 = extrema.square_distance(i);
        }
    }
    let _ = dist2;
    Some(extrema.point(index_min))
}

/// `Contap_HContTool::Tolerance(V, C)` (`cxx:306-311`) = `V->Resolution(C)`.
pub fn tolerance(v: &HVertex, c: &Curve2dAdaptor) -> f64 {
    v.resolution(c)
}

/// `Contap_HContTool::Parameter(V, C)` (`cxx:313-318`) = `V->Parameter(C)`.
pub fn parameter(v: &HVertex, c: &Curve2dAdaptor) -> f64 {
    v.parameter(c)
}

/// `Contap_HCurve2dTool::NbSamples(C, U0, U1)` (`Contap_HCurve2dTool.cxx:20-53`),
/// with the `10` default of that class (not the `25` of
/// `IntPatch_HCurve2dTool`).
pub fn curve_nb_samples(c: &dyn Curve2d, u0: f64, u1: f64) -> i32 {
    let mut nbs = 10.0f64;
    match curve2d_kind(c) {
        GeomAbsCurveType::Line => nbs = 2.0,
        GeomAbsCurveType::BezierCurve => nbs = 3.0 + c.bezier_nb_poles().unwrap_or(0) as f64,
        GeomAbsCurveType::BSplineCurve => {
            nbs = c.bspline_nb_knots().unwrap_or(0) as f64;
            nbs *= c.bspline_degree().unwrap_or(0) as f64;
            nbs *= c.last_parameter() - c.first_parameter();
            nbs /= u1 - u0;
            if nbs < 2.0 {
                nbs = 2.0;
            }
        }
        _ => {}
    }
    if nbs > 50.0 {
        nbs = 50.0;
    }
    nbs as i32
}

/// `Contap_HCurve2dTool::Value(C, U)` / `D1(C, U, P, V1)`:
/// `Adaptor2d_Curve2d::Value` / `D1` (`Contap_HCurve2dTool.lxx`). The port's
/// [`Curve2dAdaptor`] carries both; a curve that was never loaded answers
/// `None` where OCCT dereferences a null handle.
pub fn arc_value(a: &Curve2dAdaptor, u: f64) -> Option<GpPnt2d> {
    a.d0(u)
}

/// `Contap_HCurve2dTool::D1(C, U, P, V1)`.
pub fn arc_d1(a: &Curve2dAdaptor, u: f64) -> Option<(GpPnt2d, GpVec2d)> {
    a.d1(u)
}

/// `RealEpsilon()`, used by `Contap_SurfProps::DerivAndNorm` /
/// `NormAndDn` for the cone-apex test (`Contap_SurfProps.cxx:189`, `:298`).
pub const CONE_APEX_EPSILON: f64 = REAL_EPSILON;
