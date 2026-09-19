//! Port of `IntRes2d_Domain`
//! (`src/ModelingAlgorithms/TKGeomAlgo/IntRes2d/IntRes2d_Domain.hxx/.cxx/.lxx`).
//!
//! The domain is the parameter window of a 2D curve handed to the
//! intersection algorithm. Its `status` bitfield is the OCCT one
//! (`Domain.lxx:19-23`): bit 0 = has first point, bit 1 = has last point,
//! bit 2 = closed.

use crate::gp::GpPnt2d;
use crate::precision::Precision;

/// `LimitInfinite` (`IntRes2d_Domain.cxx:28-31`).
fn limit_infinite(val: f64) -> f64 {
    let inf_val = Precision::INFINITE;
    if val.abs() > inf_val {
        if val > 0.0 {
            inf_val
        } else {
            -inf_val
        }
    } else {
        val
    }
}

/// `IntRes2d_Domain` (`IntRes2d_Domain.hxx:40-141`).
#[derive(Clone, Copy, Debug)]
pub struct IntRes2dDomain {
    /// `status` (`hxx:133`)
    status: i32,
    /// `first_param` (`hxx:134`)
    first_param: f64,
    /// `last_param` (`hxx:135`)
    last_param: f64,
    /// `first_tol` (`hxx:136`)
    first_tol: f64,
    /// `last_tol` (`hxx:137`)
    last_tol: f64,
    /// `first_point` (`hxx:138`)
    first_point: GpPnt2d,
    /// `last_point` (`hxx:139`)
    last_point: GpPnt2d,
    /// `periodfirst` (`hxx:140`)
    periodfirst: f64,
    /// `periodlast` (`hxx:141`)
    periodlast: f64,
}

impl Default for IntRes2dDomain {
    /// `IntRes2d_Domain()` (`IntRes2d_Domain.cxx:35-46`): an infinite domain.
    fn default() -> Self {
        Self::new()
    }
}

impl IntRes2dDomain {
    /// `IntRes2d_Domain()` (`IntRes2d_Domain.cxx:35-46`).
    pub fn new() -> Self {
        Self {
            status: 0,
            first_param: 0.0,
            last_param: 0.0,
            first_tol: 0.0,
            last_tol: 0.0,
            first_point: GpPnt2d::new(0.0, 0.0),
            last_point: GpPnt2d::new(0.0, 0.0),
            periodfirst: 0.0,
            periodlast: 0.0,
        }
    }

    /// `IntRes2d_Domain(Pnt1, Par1, Tol1, Pnt2, Par2, Tol2)`
    /// (`IntRes2d_Domain.cxx:55-66`): a bounded domain.
    pub fn bounded(
        pnt1: &GpPnt2d,
        par1: f64,
        tol1: f64,
        pnt2: &GpPnt2d,
        par2: f64,
        tol2: f64,
    ) -> Self {
        let mut d = Self::new();
        d.set_bounded(pnt1, par1, tol1, pnt2, par2, tol2);
        d
    }

    /// `IntRes2d_Domain(Pnt, Par, Tol, First)` (`IntRes2d_Domain.cxx:92-107`):
    /// a semi-infinite domain.
    pub fn semi_infinite(pnt: &GpPnt2d, par: f64, tol: f64, first: bool) -> Self {
        let mut d = Self::new();
        d.set_semi_infinite(pnt, par, tol, first);
        d
    }

    /// `SetValues()` (`IntRes2d_Domain.cxx:48-52`): an infinite domain.
    pub fn set_infinite(&mut self) {
        self.status = 0;
        self.periodfirst = 0.0;
        self.periodlast = 0.0;
    }

    /// `SetValues(Pnt1, Par1, Tol1, Pnt2, Par2, Tol2)`
    /// (`IntRes2d_Domain.cxx:74-90`).
    pub fn set_bounded(
        &mut self,
        pnt1: &GpPnt2d,
        par1: f64,
        tol1: f64,
        pnt2: &GpPnt2d,
        par2: f64,
        tol2: f64,
    ) {
        self.status = 3;
        self.periodfirst = 0.0;
        self.periodlast = 0.0;

        self.first_param = limit_infinite(par1);
        self.first_point.set_coord(limit_infinite(pnt1.x()), limit_infinite(pnt1.y()));
        self.first_tol = tol1;

        self.last_param = limit_infinite(par2);
        self.last_point.set_coord(limit_infinite(pnt2.x()), limit_infinite(pnt2.y()));
        self.last_tol = tol2;
    }

    /// `SetValues(Pnt, Par, Tol, First)` (`IntRes2d_Domain.cxx:109-127`).
    pub fn set_semi_infinite(&mut self, pnt: &GpPnt2d, par: f64, tol: f64, first: bool) {
        self.periodfirst = 0.0;
        self.periodlast = 0.0;
        if first {
            self.status = 1;
            self.first_param = limit_infinite(par);
            self.first_point
                .set_coord(limit_infinite(pnt.x()), limit_infinite(pnt.y()));
            self.first_tol = tol;
        } else {
            self.status = 2;
            self.last_param = limit_infinite(par);
            self.last_point
                .set_coord(limit_infinite(pnt.x()), limit_infinite(pnt.y()));
            self.last_tol = tol;
        }
    }

    /// `SetEquivalentParameters(p_first, p_last)` (`Domain.lxx:24-34`).
    /// OCCT throws `Standard_DomainError` when the domain is not bounded.
    pub fn set_equivalent_parameters(&mut self, p_first: f64, p_last: f64) {
        if (self.status & 3) != 3 {
            panic!("IntRes2d_Domain::SetEquivalentParameters: domain not bounded");
        }
        self.status |= 4;
        self.periodfirst = p_first;
        self.periodlast = p_last;
    }

    /// `HasFirstPoint()` (`Domain.lxx:36-39`).
    pub fn has_first_point(&self) -> bool {
        (self.status & 1) != 0
    }

    /// `FirstParameter()` (`Domain.lxx:41-48`).
    pub fn first_parameter(&self) -> f64 {
        if (self.status & 1) == 0 {
            panic!("IntRes2d_Domain::FirstParameter: no first point");
        }
        self.first_param
    }

    /// `FirstPoint()` (`Domain.lxx:50-57`).
    pub fn first_point(&self) -> &GpPnt2d {
        if (self.status & 1) == 0 {
            panic!("IntRes2d_Domain::FirstPoint: no first point");
        }
        &self.first_point
    }

    /// `FirstTolerance()` (`Domain.lxx:59-65`).
    pub fn first_tolerance(&self) -> f64 {
        if (self.status & 1) == 0 {
            panic!("IntRes2d_Domain::FirstTolerance: no first point");
        }
        self.first_tol
    }

    /// `HasLastPoint()` (`Domain.lxx:67-70`).
    pub fn has_last_point(&self) -> bool {
        (self.status & 2) != 0
    }

    /// `LastParameter()` (`Domain.lxx:72-79`).
    pub fn last_parameter(&self) -> f64 {
        if (self.status & 2) == 0 {
            panic!("IntRes2d_Domain::LastParameter: no last point");
        }
        self.last_param
    }

    /// `LastPoint()` (`Domain.lxx:81-88`).
    pub fn last_point(&self) -> &GpPnt2d {
        if (self.status & 2) == 0 {
            panic!("IntRes2d_Domain::LastPoint: no last point");
        }
        &self.last_point
    }

    /// `LastTolerance()` (`Domain.lxx:90-97`).
    pub fn last_tolerance(&self) -> f64 {
        if (self.status & 2) == 0 {
            panic!("IntRes2d_Domain::LastTolerance: no last point");
        }
        self.last_tol
    }

    /// `IsClosed()` (`Domain.lxx:99-102`).
    pub fn is_closed(&self) -> bool {
        (self.status & 4) != 0
    }

    /// `EquivalentParameters(p_first, p_last)` (`Domain.lxx:104-109`).
    pub fn equivalent_parameters(&self) -> (f64, f64) {
        (self.periodfirst, self.periodlast)
    }
}
