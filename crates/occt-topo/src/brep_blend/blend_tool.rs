//! `BRepBlend_BlendTool` — the geometric queries the walking builders make on
//! the restriction curves and the support surfaces.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_BlendTool.hxx:35-90`,
//! `BRepBlend_BlendTool.cxx:17-142`, `BRepBlend_BlendTool.lxx:17-58`
//! (OCCT 8.0.0).
//!
//! Translation notes:
//! - `Project` / `Inters` declare `Adaptor3d_Surface` parameters that OCCT
//!   never names and never reads (`cxx:22-25`, `:80-83`); the port omits them.
//! - `Parameter` / `Tolerance` call `V->Parameter(C)` / `V->Resolution(C)`
//!   (`lxx:24-37`), and `BRepTopAdaptor_HVertex` reads both through its own
//!   `BRepAdaptor_Curve2d` argument; the port receives that concrete adaptor
//!   directly, the same down-cast [`crate::brep_top_adaptor_hvertex::HVertex`]
//!   documents.
//! - `SingularOnUMin/UMax/VMin/VMax` return `false` unconditionally
//!   (`lxx:19-34`), and `CurveOnSurf` returns its argument unchanged
//!   (`lxx:51-58`); both are kept as properties of this tool, not as surface
//!   queries.

use std::sync::Arc;

use occt_core::gp::{GpDir2d, GpPnt2d, GpVec2d};
use occt_geom2d::extrema2d::EPcOfExtPc2d;
use occt_geom2d::geom2d_int::Geom2dIntGInter;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;
use occt_geom2d::{Curve2d, Geom2dLine};

use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_top_adaptor_hvertex::HVertex;

/// `Inters` zero-length guard (`cxx:91`).
const INTERS_MIN_MAGNITUDE: f64 = 1.0e-8;
/// `Project` sampling count, parameter tolerance and function tolerance
/// (`cxx:53-55`).
const PROJECT_NB_SAMPLE: i32 = 20;
const PROJECT_EPS_X: f64 = 1.0e-8;
const PROJECT_TOL_F: f64 = 1.0e-5;
/// `Inters` segment extension around the query segment (`cxx:97`).
const INTERS_SEGMENT_LOW: f64 = -0.01;
const INTERS_SEGMENT_HIGH: f64 = 1.01;
/// `NbSamplesU` / `NbSamplesV` (`cxx:107-119`).
const NB_SAMPLES_UV: usize = 10;

/// `class BRepBlend_BlendTool` (`hxx:35-90`): every member is `static`.
pub struct BlendTool;

impl BlendTool {
    /// `Project(P, S, C, Paramproj, Dist)` (`cxx:28-78`): the orthogonal
    /// projection of `p` on `c`.
    ///
    /// The result is always "successful" (the two early exits at `cxx:57-61`
    /// also return `true`, per BUC60360); `param_proj` and `dist` are written in
    /// every path.
    pub fn project(p: &GpPnt2d, c: &dyn Curve2d, param_proj: &mut f64, dist: &mut f64) -> bool {
        *param_proj = c.first_parameter();
        let mut p2d = c.d0(*param_proj);
        *dist = p2d.distance(p);

        let t = c.last_parameter();
        p2d = c.d0(t);
        if p2d.distance(p) < *dist {
            *param_proj = t;
            *dist = p2d.distance(p);
        }

        let extrema = EPcOfExtPc2d::new(
            p,
            c,
            PROJECT_NB_SAMPLE,
            PROJECT_EPS_X,
            PROJECT_TOL_F,
        );
        if !extrema.is_done() {
            return true;
        }

        let nbext = extrema.nb_ext();
        let mut a_dist2 = *dist * *dist;
        for i in 1..=nbext {
            if extrema.square_distance(i) < a_dist2 {
                a_dist2 = extrema.square_distance(i);
                *param_proj = extrema.point(i).0;
            }
        }
        *dist = a_dist2.sqrt();

        true
    }

    /// `Inters(P1, P2, S, C, Param, Dist)` (`cxx:83-123`): the first
    /// intersection of the segment `P1P2` with `c`.
    pub fn inters(p1: &GpPnt2d, p2: &GpPnt2d, c: &dyn Curve2d, param: &mut f64, dist: &mut f64) -> bool {
        let tol = 1.0e-8;
        let v = GpVec2d::new(p2.x() - p1.x(), p2.y() - p1.y());
        let mag = v.magnitude();
        if mag < INTERS_MIN_MAGNITUDE {
            return false;
        }

        // `gp_Dir2d d(v)` (`cxx:94`); `mag >= 1.e-8` above, so the direction
        // constructor cannot take its zero-vector arm.
        let Ok(d) = GpDir2d::from_vec2d(&v) else {
            return false;
        };
        let bid = Arc::new(Geom2dLine::from_pnt_dir(*p1, d));
        let seg = Geom2dTrimmedCurve::new(
            bid,
            INTERS_SEGMENT_LOW * mag,
            INTERS_SEGMENT_HIGH * mag,
        );

        let mut inter = Geom2dIntGInter::new();
        inter.perform_curves(&seg, c, tol, tol);
        if !inter.is_done() {
            return false;
        }
        if inter.nb_points() == 0 {
            return false;
        }

        let ip = inter.point(1);
        *param = ip.param_on_second();
        *dist = p1.distance(ip.value());
        true
    }

    /// `Parameter(V, A)` (`lxx:34-37`): `V->Parameter(C)`.
    pub fn parameter(v: &HVertex, a: &Curve2dAdaptor) -> f64 {
        v.parameter(a)
    }

    /// `Tolerance(V, A)` (`lxx:26-32`): `V->Resolution(A)`; the commented-out
    /// `Resolution(A, 1.e-6)` variant is not taken (`lxx:28`).
    pub fn tolerance(v: &HVertex, a: &Curve2dAdaptor) -> f64 {
        v.resolution(a)
    }

    /// `SingularOnUMin(S)` (`lxx:19-22`): `false`.
    pub fn singular_on_u_min() -> bool {
        false
    }

    /// `SingularOnUMax(S)` (`lxx:24-27`): `false`.
    pub fn singular_on_u_max() -> bool {
        false
    }

    /// `SingularOnVMin(S)` (`lxx:29-32`): `false`.
    pub fn singular_on_v_min() -> bool {
        false
    }

    /// `SingularOnVMax(S)` (`lxx:34-37`): `false`.
    pub fn singular_on_v_max() -> bool {
        false
    }

    /// `NbSamplesU(S, u1, u2)` (`cxx:110-119`): `10`.
    pub fn nb_samples_u() -> usize {
        NB_SAMPLES_UV
    }

    /// `NbSamplesV(S, v1, v2)` (`cxx:101-108`): `10`.
    pub fn nb_samples_v() -> usize {
        NB_SAMPLES_UV
    }

    /// `Bounds(C, Ufirst, Ulast)` (`cxx:125-131`): the parametric limits on the
    /// arc, which OCCT documents as finite even for an infinite arc
    /// (`hxx:78-82`).
    pub fn bounds(c: &dyn Curve2d) -> (f64, f64) {
        (c.first_parameter(), c.last_parameter())
    }

    /// `CurveOnSurf(C, S)` (`lxx:51-58`): the argument unchanged.
    pub fn curve_on_surf(c: Arc<dyn Curve2d>) -> Arc<dyn Curve2d> {
        c
    }
}
