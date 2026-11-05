//! `Geom2dHatch_Hatching`: one hatching curve and its restriction points and
//! domains.
//!
//! Sources: `ModelingAlgorithms/TKGeomAlgo/Geom2dHatch/Geom2dHatch_Hatching.hxx:37-141`,
//! `.cxx:26-332` (OCCT 8.0.0).
//!
//! `Geom2dAdaptor_Curve myCurve` becomes an `Arc<dyn Curve2d>`, as in
//! [`super::element`].

use std::sync::Arc;

use occt_core::gp::GpPnt2d;
use occt_core::precision::INFINITE;
use occt_geom2d::curve::Curve2d;

use crate::hatch_gen::{Domain, ErrorStatus, PointOnHatching};

/// `Geom2dHatch_Hatching` (`hxx:37-141`).
#[derive(Clone)]
pub struct Hatching {
    /// `myCurve`
    my_curve: Arc<dyn Curve2d>,
    /// `myTrimDone`
    my_trim_done: bool,
    /// `myTrimFailed`
    my_trim_failed: bool,
    /// `myPoints`
    my_points: Vec<PointOnHatching>,
    /// `myIsDone`
    my_is_done: bool,
    /// `myStatus`
    my_status: ErrorStatus,
    /// `myDomains`
    my_domains: Vec<Domain>,
}

impl Hatching {
    /// `Geom2dHatch_Hatching(Curve)` (`cxx:37-45`).
    pub fn new(curve: Arc<dyn Curve2d>) -> Self {
        Self {
            my_curve: curve,
            my_trim_done: false,
            my_trim_failed: false,
            my_points: Vec::new(),
            my_is_done: false,
            my_status: ErrorStatus::NoProblem,
            my_domains: Vec::new(),
        }
    }

    /// `Curve()` (`cxx:52-55`).
    pub fn curve(&self) -> &Arc<dyn Curve2d> {
        &self.my_curve
    }

    /// `ChangeCurve()` (`cxx:62-65`).
    pub fn change_curve(&mut self) -> &mut Arc<dyn Curve2d> {
        &mut self.my_curve
    }

    /// `TrimDone(Flag)` (`cxx:72-75`).
    pub fn set_trim_done(&mut self, flag: bool) {
        self.my_trim_done = flag;
    }

    /// `TrimDone()` (`cxx:82-85`).
    pub fn trim_done(&self) -> bool {
        self.my_trim_done
    }

    /// `TrimFailed(Flag)` (`cxx:93-100`): setting the flag also records
    /// `HatchGen_TrimFailure` in the status.
    pub fn set_trim_failed(&mut self, flag: bool) {
        self.my_trim_failed = flag;
        if self.my_trim_failed {
            self.my_status = ErrorStatus::TrimFailure;
        }
    }

    /// `TrimFailed()` (`cxx:107-110`).
    pub fn trim_failed(&self) -> bool {
        self.my_trim_failed
    }

    /// `IsDone(Flag)` (`cxx:118-121`).
    pub fn set_is_done(&mut self, flag: bool) {
        self.my_is_done = flag;
    }

    /// `IsDone()` (`cxx:128-131`).
    pub fn is_done(&self) -> bool {
        self.my_is_done
    }

    /// `Status(theStatus)` (`cxx:136-139`).
    pub fn set_status(&mut self, status: ErrorStatus) {
        self.my_status = status;
    }

    /// `Status()` (`cxx:143-146`).
    pub fn status(&self) -> ErrorStatus {
        self.my_status
    }

    /// `AddPoint(Point, Confusion)` (`cxx:152-190`): insert the point in
    /// parameter order, or merge its points on element into the point already
    /// stored at that parameter. Adding a point clears the domains when they
    /// were computed.
    pub fn add_point(&mut self, point: &PointOnHatching, confusion: f64) {
        let nb_points = self.my_points.len();
        let mut ipnth = 0usize;
        while ipnth < nb_points {
            if !self.my_points[ipnth].is_lower(point, confusion) {
                break;
            }
            ipnth += 1;
        }
        if ipnth >= nb_points {
            self.my_points.push(point.clone()); // `cxx:167-169`
        } else if self.my_points[ipnth].is_greater(point, confusion) {
            self.my_points.insert(ipnth, point.clone()); // `cxx:176`
        } else {
            for ipnte in 0..point.nb_points() {
                let pnte = point.point(ipnte + 1).clone();
                self.my_points[ipnth].add_point(&pnte, confusion);
            }
        }
        if self.my_is_done {
            self.clr_domains();
        }
    }

    /// `NbPoints()` (`cxx:196-199`).
    pub fn nb_points(&self) -> usize {
        self.my_points.len()
    }

    /// `Point(Index)` (`cxx:206-209`), 1-based; OCCT raises `OutOfRange`
    /// outside `1..=NbPoints`, the port panics through the index access.
    pub fn point(&self, index: usize) -> &PointOnHatching {
        &self.my_points[index - 1]
    }

    /// `ChangePoint(Index)` (`cxx:216-219`), 1-based.
    pub fn change_point(&mut self, index: usize) -> &mut PointOnHatching {
        &mut self.my_points[index - 1]
    }

    /// `RemPoint(Index)` (`cxx:226-233`), 1-based.
    pub fn rem_point(&mut self, index: usize) {
        if self.my_is_done {
            self.clr_domains();
        }
        self.my_points.remove(index - 1);
    }

    /// `ClrPoints()` (`cxx:241-256`).
    pub fn clr_points(&mut self) {
        if self.my_is_done {
            self.clr_domains();
        }
        for point in self.my_points.iter_mut() {
            point.clr_points();
        }
        self.my_points.clear();
        self.my_trim_done = false;
        self.my_trim_failed = false;
    }

    /// `AddDomain(Domain)` (`cxx:262-265`).
    pub fn add_domain(&mut self, domain: &Domain) {
        self.my_domains.push(domain.clone());
    }

    /// `NbDomains()` (`cxx:272-275`).
    pub fn nb_domains(&self) -> usize {
        self.my_domains.len()
    }

    /// `Domain(Index)` (`cxx:282-285`), 1-based.
    pub fn domain(&self, index: usize) -> &Domain {
        &self.my_domains[index - 1]
    }

    /// `ChangeDomain(Index)` (`cxx:288-290`), 1-based.
    pub fn change_domain(&mut self, index: usize) -> &mut Domain {
        &mut self.my_domains[index - 1]
    }

    /// `RemDomain(Index)` (`cxx:293-296`), 1-based.
    pub fn rem_domain(&mut self, index: usize) {
        self.my_domains.remove(index - 1);
    }

    /// `ClrDomains()` (`cxx:302-306`): clears the domains *and* the done flag.
    pub fn clr_domains(&mut self) {
        self.my_domains.clear();
        self.my_is_done = false;
    }

    /// `ClassificationPoint()` (`cxx:312-332`).
    ///
    /// The probe point is the last parameter when it is finite, the first one
    /// when only the first is finite, and the origin of the curve when both are
    /// infinite.
    pub fn classification_point(&self) -> GpPnt2d {
        let a = self.my_curve.first_parameter();
        let b = self.my_curve.last_parameter();
        let t = if b >= INFINITE {
            if a <= -INFINITE {
                0.
            } else {
                a
            }
        } else {
            b
        };
        self.my_curve.value(t)
    }
}
