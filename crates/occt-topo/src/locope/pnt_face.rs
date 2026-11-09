//! `LocOpe_PntFace` (`LocOpe_PntFace.hxx:29-79`).
//!
//! A point / face / orientation / parameter tuple, used by
//! [`crate::locope::cs_intersector::CSIntersector`] and
//! [`crate::locope::curve_shape_intersector::CurveShapeIntersector`]. The
//! reference keeps the six fields in the header (`hxx:71-78`) with the default
//! constructor leaving the three parameters at `0.0` (`hxx:36-41`).

use occt_core::gp::GpPnt;

use crate::abs::{Orientation, ShapeType};
use crate::shape::{Face, TopoShape};

#[derive(Debug, Clone)]
pub struct PntFace {
    pnt: GpPnt,
    face: Face,
    orientation: Orientation,
    parameter: f64,
    u_parameter: f64,
    v_parameter: f64,
}

impl Default for PntFace {
    /// `LocOpe_PntFace()` (`hxx:36-41`): "empty constructor, useful only for
    /// the list", `myPar`/`myUPar`/`myVPar` set to `0.0`.
    fn default() -> Self {
        Self {
            pnt: GpPnt::new(0.0, 0.0, 0.0),
            face: Face::wrap(TopoShape::new(ShapeType::Face))
                .expect("a FACE-typed TopoShape wraps as a Face"),
            orientation: Orientation::Forward,
            parameter: 0.0,
            u_parameter: 0.0,
            v_parameter: 0.0,
        }
    }
}

impl PntFace {
    /// `LocOpe_PntFace(P, F, Or, Param, UPar, VPar)` (`hxx:42-56`).
    pub fn new(
        pnt: GpPnt,
        face: Face,
        orientation: Orientation,
        parameter: f64,
        u_parameter: f64,
        v_parameter: f64,
    ) -> Self {
        Self {
            pnt,
            face,
            orientation,
            parameter,
            u_parameter,
            v_parameter,
        }
    }

    /// `Pnt()` (`hxx:58`).
    pub fn pnt(&self) -> &GpPnt {
        &self.pnt
    }

    /// `Face()` (`hxx:60`).
    pub fn face(&self) -> &Face {
        &self.face
    }

    /// `Orientation()` (`hxx:62`).
    pub fn orientation(&self) -> Orientation {
        self.orientation
    }

    /// `ChangeOrientation()` (`hxx:64`).
    pub fn set_orientation(&mut self, orientation: Orientation) {
        self.orientation = orientation;
    }

    /// `Parameter()` (`hxx:66`).
    pub fn parameter(&self) -> f64 {
        self.parameter
    }

    /// `UParameter()` (`hxx:68`).
    pub fn u_parameter(&self) -> f64 {
        self.u_parameter
    }

    /// `VParameter()` (`hxx:70`).
    pub fn v_parameter(&self) -> f64 {
        self.v_parameter
    }
}
