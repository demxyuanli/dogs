//! Locals of `GeomFill_Sweep::BuildKPart` shared by its two branches.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_Sweep.cxx:577-606` (OCCT 8.0.0):
//! `Ok`, `isUPeriodic`, `isVPeriodic`, `IsTrsf`, `S`, `Error`, `SError`,
//! `UFirst`, `VFirst`, `ULast`, `VLast`, `Tol` plus `myExchUV` /
//! `isUReversed` of the sweep itself.
//!
//! Translation notes:
//! - `occ::handle<Geom_Surface> S` becomes `Option<Arc<dyn Surface>>`.
//! - `gp_GTrsf Tf; Tf.SetVectorialPart(M); Tf.SetTranslationPart(V.XYZ())`
//!   followed by `gp_Trsf Tf2; Tf2.SetValues(Tf(1, 1), ...)` becomes
//!   [`trsf_from_matrix`], which feeds `M` and `V` straight into
//!   `gp_Trsf::SetValues` because a `gp_GTrsf` with a linear vectorial part
//!   has exactly those 12 coefficients.

use std::sync::Arc;

use occt_core::gp::{GpAx3, GpDir, GpMat, GpPnt, GpTrsf, GpVec};
use occt_core::precision::{PCONFUSION, REAL_SMALL};

use crate::surface::Surface;

/// `Precision::Infinite()` (`Precision.hxx`): the parameter OCCT reports for
/// the unbounded end of a line. The ported curves use `f64::INFINITY`, so the
/// KPart arithmetic maps the two spellings.
pub(super) const INFINITE: f64 = 1.0e100;

/// `gp::Resolution()` (`gp.cxx:36-40`): `RealSmall()`.
pub(super) const GP_RESOLUTION: f64 = REAL_SMALL;

/// `Precision::PConfusion()` (`Precision.hxx`).
pub(super) const P_CONFUSION: f64 = PCONFUSION;

/// Map the port's `f64::INFINITY` parameter ends onto `Precision::Infinite()`,
/// the value OCCT's `GeomAdaptor_Curve::FirstParameter` reports for a line.
pub(super) fn occt_infinite(u: f64) -> f64 {
    if u == f64::INFINITY {
        INFINITE
    } else if u == f64::NEG_INFINITY {
        -INFINITE
    } else {
        u
    }
}

/// `gp_Vec` -> `gp_Dir` (the implicit `gp_Dir(const gp_Vec&)` conversion,
/// `gp_Dir.hxx:338-343`). `Err` when the vector is null, which is the
/// `gp_Dir` `ConstructionError` (`gp_Dir.hxx:345-350`).
pub(super) fn dir_of(v: &GpVec) -> Result<GpDir, &'static str> {
    GpDir::from_vec(v)
}

/// `gp_GTrsf Tf; Tf.SetVectorialPart(M); Tf.SetTranslationPart(V.XYZ());`
/// then `Tf2.SetValues(Tf(1, 1), ..., Tf(3, 4))`
/// (`GeomFill_Sweep.cxx:602-620`, `cxx:831-853`, `cxx:501-539`).
///
/// `Tf(i, 4)` is the translation part of the `gp_GTrsf`, so the 12 arguments
/// are the 3x3 matrix followed by `V`. `Err` is the `Standard_ConstructionError`
/// of `gp_Trsf::SetValues` (`gp_Trsf.cxx:366-370`): a null determinant.
pub(super) fn trsf_from_matrix(m: &GpMat, v: &GpVec) -> Result<GpTrsf, &'static str> {
    let mut t = GpTrsf::identity();
    t.set_values(
        m.value(1, 1),
        m.value(1, 2),
        m.value(1, 3),
        v.x(),
        m.value(2, 1),
        m.value(2, 2),
        m.value(2, 3),
        v.y(),
        m.value(3, 1),
        m.value(3, 2),
        m.value(3, 3),
        v.z(),
    )?;
    Ok(t)
}

/// `gp_Ax3(P, N, Vx)` (`gp_Ax3.hxx:84-91`): `X = N x (Vx x N)`, `Y = N x X`.
/// `Err` when `N` and `Vx` are parallel, the `ConstructionError` of that
/// constructor.
pub(super) fn ax3_from_point_dir_x(p: GpPnt, n: GpDir, vx: &GpDir) -> Result<GpAx3, &'static str> {
    let vxdir = n.cross_cross(vx)?;
    let vydir = n.crossed(&vxdir)?;
    Ok(GpAx3 {
        axis: occt_core::gp::GpAx1::new(p, n),
        vxdir,
        vydir,
    })
}

/// The `BuildKPart` locals (`GeomFill_Sweep.cxx:579-591`).
pub(crate) struct KPartState {
    /// `S` (`cxx:586`): the elementary surface before trimming.
    pub s: Option<Arc<dyn Surface>>,
    /// `Ok` (`cxx:579`).
    pub ok: bool,
    /// `isUPeriodic` (`cxx:580`).
    pub is_u_periodic: bool,
    /// `isVPeriodic` (`cxx:581`).
    pub is_v_periodic: bool,
    /// `UFirst` (`cxx:590`), initialised to `0`.
    pub u_first: f64,
    /// `ULast` (`cxx:590`), initialised to `0`.
    pub u_last: f64,
    /// `VFirst` (`cxx:590`), initialised to `First`.
    pub v_first: f64,
    /// `VLast` (`cxx:590`), initialised to `Last`.
    pub v_last: f64,
    /// `error` (`cxx:589`), the error `IsTranslation` / `IsRotation` /
    /// `IsConstant` / `IsConicalLaw` report.
    pub error: f64,
    /// `SError` (`cxx:589`), the error kept on the sweep.
    pub s_error: f64,
    /// `Tol = std::min(Tol3d, BoundTol)` (`cxx:591`).
    pub tol: f64,
    /// `myExchUV` (`cxx:166`).
    pub exch_uv: bool,
    /// `isUReversed` (`cxx:166`).
    pub u_reversed: bool,
}

impl KPartState {
    /// The initialisers of `GeomFill_Sweep.cxx:579-591`.
    pub(super) fn new(tol: f64, first: f64, last: f64, s_error: f64) -> Self {
        Self {
            s: None,
            ok: false,
            is_u_periodic: false,
            is_v_periodic: false,
            u_first: 0.0,
            u_last: 0.0,
            v_first: first,
            v_last: last,
            error: 0.0,
            s_error,
            tol,
            exch_uv: false,
            u_reversed: false,
        }
    }
}
