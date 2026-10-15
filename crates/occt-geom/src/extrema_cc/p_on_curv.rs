//! `Extrema_POnCurv` (`Extrema_POnCurv.hxx:25-60`): a parameter and the 3D
//! point on the curve at that parameter.
//!
//! The 2D twin is `occt_geom2d::geom2d_int::POnCurv2d`
//! (`geom2d_int/gfunc_ext_pc.rs:22-41`, = `Extrema_POnCurv2d`), and this port
//! keeps the same shape (`new(parameter, point)` plus the `Parameter()` /
//! `Value()` accessors, public fields). It lands with its first consumer, the
//! `Extrema_GFuncExtCC` 3D instantiation (`Extrema_CCLocFOfLocECC.hxx:29-36`
//! instantiates `ThePOnC = Extrema_POnCurv`).

use occt_core::gp::GpPnt;

/// `Extrema_POnCurv` (`Extrema_POnCurv.hxx:25-60`).
#[derive(Clone, Copy, Debug)]
pub struct POnCurv {
    /// `myU` (`hxx:58`).
    pub parameter: f64,
    /// `myP` (`hxx:59`).
    pub point: GpPnt,
}

impl POnCurv {
    /// `Extrema_POnCurv(theU, theP)` (`hxx:38-42`). The default constructor
    /// (`hxx:31-35`, parameter `0` on the origin point) is `Default`.
    pub fn new(parameter: f64, point: GpPnt) -> Self {
        Self { parameter, point }
    }

    /// `Parameter()` (`hxx:54`).
    pub fn parameter(&self) -> f64 {
        self.parameter
    }

    /// `Value()` (`hxx:51`).
    pub fn value(&self) -> GpPnt {
        self.point
    }
}

impl Default for POnCurv {
    /// `Extrema_POnCurv()` (`hxx:31-35`).
    fn default() -> Self {
        Self {
            parameter: 0.0,
            point: GpPnt::zero(),
        }
    }
}
