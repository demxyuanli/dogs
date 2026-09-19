//! `IntCurve_PConic` (TKGeomAlgo, `IntCurve/IntCurve_PConic.hxx`, `_0.cxx`,
//! `.lxx`).
//!
//! A conic of `gp` reduced to the three quantities the `IntCurve` algorithms
//! actually use: an `gp_Ax22d` and two parameters. The pair is interpreted
//! according to [`CurveType`]:
//!
//! | `TypeCurve()`     | `Param1()`     | `Param2()`     |
//! |-------------------|----------------|----------------|
//! | `Line`            | 0              | 0              |
//! | `Circle`          | radius         | 0              |
//! | `Ellipse`         | major radius   | minor radius   |
//! | `Parabola`        | focal          | 0              |
//! | `Hyperbola`       | major radius   | minor radius   |
use crate::gp::{GpAx22d, GpCirc2d, GpElips2d, GpHypr2d, GpLin2d, GpParab2d};
use crate::kernel::geomabs::CurveType;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IntCurvePConic {
    axe: GpAx22d,
    prm1: f64,
    prm2: f64,
    the_eps_x: f64,
    the_accuracy: i32,
    type_curve: CurveType,
}

impl IntCurvePConic {
    /// `IntCurve_PConic(const gp_Elips2d&)` (`cxx:28-36`).
    pub fn from_elips2d(e: &GpElips2d) -> Self {
        Self {
            axe: e.pos,
            prm1: e.major_radius,
            prm2: e.minor_radius,
            the_eps_x: 0.00000001,
            the_accuracy: 20,
            type_curve: CurveType::Ellipse,
        }
    }

    /// `IntCurve_PConic(const gp_Hypr2d&)` (`cxx:38-46`). The default accuracy
    /// for a hyperbola is 50, not 20 as for every other conic.
    pub fn from_hypr2d(h: &GpHypr2d) -> Self {
        Self {
            axe: h.pos,
            prm1: h.major_radius,
            prm2: h.minor_radius,
            the_eps_x: 0.00000001,
            the_accuracy: 50,
            type_curve: CurveType::Hyperbola,
        }
    }

    /// `IntCurve_PConic(const gp_Circ2d&)` (`cxx:48-55`).
    pub fn from_circ2d(c: &GpCirc2d) -> Self {
        Self {
            axe: c.pos,
            prm1: c.radius,
            prm2: 0.0,
            the_eps_x: 0.00000001,
            the_accuracy: 20,
            type_curve: CurveType::Circle,
        }
    }

    /// `IntCurve_PConic(const gp_Parab2d&)` (`cxx:57-64`).
    pub fn from_parab2d(p: &GpParab2d) -> Self {
        Self {
            axe: p.pos,
            prm1: p.focal,
            prm2: 0.0,
            the_eps_x: 0.00000001,
            the_accuracy: 20,
            type_curve: CurveType::Parabola,
        }
    }

    /// `IntCurve_PConic(const gp_Lin2d&)` (`cxx:66-74`). The X axis of the
    /// placement is the line's `Position()`; the Y axis is deduced right-handed
    /// by `gp_Ax22d(const gp_Ax2d&)` (`gp_Ax22d.hxx:98-111`).
    pub fn from_lin2d(l: &GpLin2d) -> Self {
        Self {
            axe: GpAx22d::from_xdir(l.pos.loc, l.pos.vdir),
            prm1: 0.0,
            prm2: 0.0,
            the_eps_x: 0.00000001,
            the_accuracy: 20,
            type_curve: CurveType::Line,
        }
    }

    /// `IntCurve_PConic::SetEpsX` (`cxx:78-81`).
    pub fn set_eps_x(&mut self, eps_dist: f64) {
        self.the_eps_x = eps_dist;
    }

    /// `IntCurve_PConic::SetAccuracy` (`cxx:83-86`).
    pub fn set_accuracy(&mut self, nb: i32) {
        self.the_accuracy = nb;
    }

    /// `IntCurve_PConic::Accuracy` (`lxx:41-44`).
    pub fn accuracy(&self) -> i32 {
        self.the_accuracy
    }

    /// `IntCurve_PConic::EpsX` (`lxx:17-20`).
    pub fn eps_x(&self) -> f64 {
        self.the_eps_x
    }

    /// `IntCurve_PConic::TypeCurve` (`lxx:36-39`).
    pub fn type_curve(&self) -> CurveType {
        self.type_curve
    }

    /// `IntCurve_PConic::Axis2` (`lxx:32-35`).
    pub fn axis2(&self) -> &GpAx22d {
        &self.axe
    }

    /// `IntCurve_PConic::Param1` (`lxx:22-25`).
    pub fn param1(&self) -> f64 {
        self.prm1
    }

    /// `IntCurve_PConic::Param2` (`lxx:27-30`).
    pub fn param2(&self) -> f64 {
        self.prm2
    }
}
