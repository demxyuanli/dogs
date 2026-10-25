//! `GeomFill_ConstantBiNormal` — the Frenet law rebuilt around a constant
//! binormal direction.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_ConstantBiNormal.hxx` and
//! `GeomFill_ConstantBiNormal.cxx:33-294` (OCCT 8.0.0).
//!
//! Translation notes: `BN` (`hxx:75`) is a `gp_Dir`, so it is kept here as a
//! unit [`GpVec`] built from the direction. `myCurve->GetType()` in
//! `IsOnlyBy3dCurve` (`cxx:260-294`) is
//! [`crate::extrema_pc::extrema_curve_type`], the same
//! `GeomAdaptor_Curve::GetType()` the rest of the port uses.

use std::sync::Arc;

use occt_core::gp::{GpDir, GpVec};
use occt_core::precision::{ANGULAR, CONFUSION};

use crate::curve::Curve;
use crate::extrema_pc::{extrema_curve_type, GeomAbsCurveType};
use crate::geom_fill::frenet::deriv::{dderiv, fderiv};
use crate::geom_fill::frenet::Frenet;
use crate::geom_fill::trihedron_law::{TrihedronD1, TrihedronD2, TrihedronLaw, TrihedronLawCurve};

/// `GeomFill_ConstantBiNormal` (`GeomFill_ConstantBiNormal.hxx:36-80`).
pub struct ConstantBiNormal {
    /// `BN` (`hxx:75`).
    bn: GpVec,
    /// `frenet` (`hxx:76`).
    frenet: Frenet,
    /// `myCurve` / `myTrimmed` (`GeomFill_TrihedronLaw.hxx:111-112`).
    data: TrihedronLawCurve,
}

impl ConstantBiNormal {
    /// `GeomFill_ConstantBiNormal` (`cxx:56-60`): `frenet = new GeomFill_Frenet`.
    pub fn new(bi_normal: &GpDir) -> Self {
        Self {
            bn: GpVec::from_xyz(bi_normal.xyz()),
            frenet: Frenet::new(),
            data: TrihedronLawCurve::new(),
        }
    }

    /// `BN` (`hxx:75`).
    pub fn bi_normal(&self) -> GpVec {
        self.bn
    }

    /// The `D0` body shared by `D0` and `GetAverageLaw`
    /// (`cxx:86-108`, `cxx:239-253`): the branches on
    /// `|BiNormal ^ Tangent|` against `Precision::Confusion()`.
    fn rebuild_d0(&self, tangent: GpVec, normal: GpVec) -> (GpVec, GpVec, GpVec) {
        let bi_normal = self.bn;
        if bi_normal.crossed(&tangent).magnitude() > CONFUSION {
            let normal = bi_normal.crossed(&tangent).normalized();
            let tangent = normal.crossed(&bi_normal);
            (tangent, normal, bi_normal)
        } else {
            let tangent = normal.crossed(&bi_normal).normalized();
            let normal = bi_normal.crossed(&tangent);
            (tangent, normal, bi_normal)
        }
    }

    /// `gp_Vec::IsNormal(Other, Precision::Angular())`.
    fn is_normal_ang(v: &GpVec, other: &GpVec) -> bool {
        (std::f64::consts::FRAC_PI_2 - v.angle(other)).abs() <= ANGULAR
    }

    /// `gp_Vec::IsParallel(Other, Precision::Angular())`.
    fn is_parallel_ang(v: &GpVec, other: &GpVec) -> bool {
        let ang = v.angle(other);
        ang <= ANGULAR || std::f64::consts::PI - ang <= ANGULAR
    }
}

impl TrihedronLaw for ConstantBiNormal {
    /// `GeomFill_ConstantBiNormal::SetCurve` (`cxx:72-81`).
    fn set_curve(&self, curve: Arc<dyn Curve>) -> bool {
        self.data.set_curve(curve.clone());
        self.frenet.set_curve(curve)
    }

    /// `GeomFill_ConstantBiNormal::Copy` (`cxx:62-70`).
    fn copy(&self) -> Arc<dyn TrihedronLaw> {
        let bn = match GpDir::from_vec(&self.bn) {
            Ok(d) => d,
            Err(_) => GpDir::default_dir(),
        };
        let copy = ConstantBiNormal::new(&bn);
        if let Some(curve) = self.data.curve() {
            TrihedronLaw::set_curve(&copy, curve);
        }
        Arc::new(copy)
    }

    /// `GeomFill_ConstantBiNormal::D0` (`cxx:83-108`).
    fn d0(&self, param: f64) -> Result<(GpVec, GpVec, GpVec), String> {
        let (tangent, normal, _bi_normal) = self.frenet.d0(param)?;
        Ok(self.rebuild_d0(tangent, normal))
    }

    /// `GeomFill_ConstantBiNormal::D1` (`cxx:112-160`).
    fn d1(&self, param: f64) -> Result<TrihedronD1, String> {
        let frame = self.frenet.d1(param)?;
        let bi_normal = self.bn;
        if bi_normal.crossed(&frame.tangent).magnitude() > CONFUSION {
            let f = bi_normal.crossed(&frame.tangent);
            let df = bi_normal.crossed(&frame.dtangent);
            let normal = f.normalized();
            let d_normal = fderiv(&f, &df);
            let tangent = normal.crossed(&bi_normal);
            let d_tangent = d_normal.crossed(&bi_normal);
            Ok(TrihedronD1 {
                tangent,
                dtangent: d_tangent,
                normal,
                dnormal: d_normal,
                binormal: bi_normal,
                dbinormal: GpVec::zero(),
            })
        } else {
            let f = frame.normal.crossed(&bi_normal);
            let df = frame.dnormal.crossed(&bi_normal);
            let tangent = f.normalized();
            let d_tangent = fderiv(&f, &df);
            let normal = bi_normal.crossed(&tangent);
            let d_normal = bi_normal.crossed(&d_tangent);
            Ok(TrihedronD1 {
                tangent,
                dtangent: d_tangent,
                normal,
                dnormal: d_normal,
                binormal: bi_normal,
                dbinormal: GpVec::zero(),
            })
        }
    }

    /// `GeomFill_ConstantBiNormal::D2` (`cxx:162-227`).
    fn d2(&self, param: f64) -> Result<TrihedronD2, String> {
        let frame = self.frenet.d2(param)?;
        let bi_normal = self.bn;
        if bi_normal.crossed(&frame.tangent).magnitude() > CONFUSION {
            let f = bi_normal.crossed(&frame.tangent);
            let df = bi_normal.crossed(&frame.dtangent);
            let d2f = bi_normal.crossed(&frame.d2tangent);
            let normal = f.normalized();
            let d_normal = fderiv(&f, &df);
            let d2_normal = dderiv(&f, &df, &d2f);
            let tangent = normal.crossed(&bi_normal);
            let d_tangent = d_normal.crossed(&bi_normal);
            let d2_tangent = d2_normal.crossed(&bi_normal);
            Ok(TrihedronD2 {
                tangent,
                dtangent: d_tangent,
                d2tangent: d2_tangent,
                normal,
                dnormal: d_normal,
                d2normal: d2_normal,
                binormal: bi_normal,
                dbinormal: GpVec::zero(),
                d2binormal: GpVec::zero(),
            })
        } else {
            let f = frame.normal.crossed(&bi_normal);
            let df = frame.dnormal.crossed(&bi_normal);
            let d2f = frame.d2normal.crossed(&bi_normal);
            let tangent = f.normalized();
            let d_tangent = fderiv(&f, &df);
            let d2_tangent = dderiv(&f, &df, &d2f);
            let normal = bi_normal.crossed(&tangent);
            let d_normal = bi_normal.crossed(&d_tangent);
            let d2_normal = bi_normal.crossed(&d2_tangent);
            Ok(TrihedronD2 {
                tangent,
                dtangent: d_tangent,
                d2tangent: d2_tangent,
                normal,
                dnormal: d_normal,
                d2normal: d2_normal,
                binormal: bi_normal,
                dbinormal: GpVec::zero(),
                d2binormal: GpVec::zero(),
            })
        }
    }

    /// `GeomFill_ConstantBiNormal::NbIntervals` (`cxx:229-232`).
    fn nb_intervals(&self, s: u8) -> i32 {
        self.frenet.nb_intervals(s)
    }

    /// `GeomFill_ConstantBiNormal::Intervals` (`cxx:234-238`).
    fn intervals(&self, s: u8) -> Vec<f64> {
        self.frenet.intervals(s)
    }

    /// `GeomFill_TrihedronLaw::SetInterval` (`cxx:57-60`).
    fn set_interval(&self, first: f64, last: f64) -> Result<(), String> {
        self.data.set_interval(first, last)?;
        self.frenet.set_interval(first, last)
    }

    /// `GeomFill_TrihedronLaw::GetInterval` (`cxx:62-67`).
    fn get_interval(&self) -> (f64, f64) {
        self.data.interval()
    }

    /// `GeomFill_ConstantBiNormal::GetAverageLaw` (`cxx:240-253`).
    fn get_average_law(&self) -> (GpVec, GpVec, GpVec) {
        let (a_tangent, a_normal, _a_bi_normal) = self.frenet.get_average_law();
        self.rebuild_d0(a_tangent, a_normal)
    }

    /// `GeomFill_ConstantBiNormal::IsConstant` (`cxx:255-258`).
    fn is_constant(&self) -> bool {
        self.frenet.is_constant()
    }

    /// `GeomFill_ConstantBiNormal::IsOnlyBy3dCurve` (`cxx:260-294`).
    fn is_only_by_3d_curve(&self) -> bool {
        let Some(curve) = self.data.curve() else {
            return false;
        };
        match extrema_curve_type(curve.as_ref()) {
            GeomAbsCurveType::Circle => {
                let axis = curve.gp_circ().map(|c| c.axis());
                match axis {
                    Some(axis) => Self::is_parallel_ang(
                        &GpVec::from_xyz(axis.direction().xyz()),
                        &self.bn,
                    ),
                    None => false,
                }
            }
            GeomAbsCurveType::Ellipse => {
                let axis = curve.gp_ellipse().map(|c| *c.axis());
                match axis {
                    Some(axis) => Self::is_parallel_ang(
                        &GpVec::from_xyz(axis.direction().xyz()),
                        &self.bn,
                    ),
                    None => false,
                }
            }
            GeomAbsCurveType::Hyperbola => {
                let axis = curve.gp_hyperbola().map(|c| *c.axis());
                match axis {
                    Some(axis) => Self::is_parallel_ang(
                        &GpVec::from_xyz(axis.direction().xyz()),
                        &self.bn,
                    ),
                    None => false,
                }
            }
            GeomAbsCurveType::Parabola => {
                let axis = curve.gp_parabola().map(|c| *c.axis());
                match axis {
                    Some(axis) => Self::is_parallel_ang(
                        &GpVec::from_xyz(axis.direction().xyz()),
                        &self.bn,
                    ),
                    None => false,
                }
            }
            GeomAbsCurveType::Line => {
                // The curve plane normal is perpendicular to the bi-normal
                // (`cxx:282-286`).
                let v = curve
                    .gp_line()
                    .map(|l| GpVec::from_xyz(l.direction().xyz()));
                match v {
                    Some(v) => Self::is_normal_ang(&v, &self.bn),
                    None => false,
                }
            }
            _ => false,
        }
    }
}
