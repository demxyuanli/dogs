//! The constant / conical answers of `GeomFill_NSections`.
//!
//! Source: `GeomFill_NSections.cxx:876-938` (`IsConstant`), `:942-950`
//! (`ConstantSection`), `:954-1011` (`IsConicalLaw`) and `:1016-1047`
//! (`CirclSection`) (OCCT 8.0.0).
//!
//! Translation notes:
//! - `GeomAdaptor_Curve::GetType()` is replaced by the curve's own
//!   [`Curve::gp_circ`] / [`Curve::gp_line`] answers, which is the same
//!   information (`GeomAdaptor_Curve::load` keeps the basis of a trimmed
//!   curve, and both port views forward those calls).
//! - `gp_Ax1::IsParallel(ax, tol)` compares directions, so the port calls
//!   `GpDir::is_parallel_tol` on `Axis().Direction()`.
//! - `gp_Circ::Transform(myTrsfs(i).Inverted())` (`cxx:966-972`) is applied to
//!   the local copies, as in OCCT.
//! - Errors keep the OCCT exception names: `Standard_DomainError` for a
//!   non-conical `CirclSection` (`cxx:1020-1023`).

use std::sync::Arc;

use occt_core::gp::{GpAx1, GpCirc, GpDir, GpLin, GpVec};
use occt_core::precision::PCONFUSION;

use crate::circle::GeomCircle;
use crate::curve::Curve;
use crate::geom_fill::n_sections::NSections;
use crate::trimmed::GeomTrimmedCurve;

/// `GeomAbs_CurveType` as far as `IsConstant` looks at it (`cxx:889-931`).
#[derive(Clone, Copy, PartialEq, Eq)]
enum CurveKind {
    Circle,
    Line,
    Other,
}

impl CurveKind {
    fn of(curve: &dyn Curve) -> Self {
        if curve.gp_circ().is_some() {
            CurveKind::Circle
        } else if curve.gp_line().is_some() {
            CurveKind::Line
        } else {
            CurveKind::Other
        }
    }
}

/// `gp_Ax1::IsParallel(other, angularTolerance)`.
fn axes_parallel(a: &GpAx1, b: &GpAx1, tolerance: f64) -> bool {
    a.direction().is_parallel_tol(b.direction(), tolerance)
}

/// The `GeomAdaptor_Curve` view of a section curve
/// (`GeomAdaptor_Curve::load`, `GeomAdaptor_Curve.cxx:252-254`).
///
/// `load` unwraps a `Geom_TrimmedCurve` to its basis, so `FirstParameter` /
/// `LastParameter` / `IsPeriodic` / `Period` of the adaptor are the **basis**
/// answers, while the port's `Geom_TrimmedCurve` reports `[0, 1]`. These are
/// the values `GeomFill_NSections.cxx:1002-1006` and `:1038-1039` read.
struct AdaptorRange {
    first: f64,
    last: f64,
    periodic: bool,
    /// `GeomAdaptor_Curve::Period` (`GeomAdaptor_Curve.cxx:596-599`):
    /// `LastParameter - FirstParameter` of the adaptor curve.
    period: f64,
}

fn adaptor_range(curve: &dyn Curve) -> AdaptorRange {
    match curve.untrimmed_basis() {
        Some((basis, first, last)) => AdaptorRange {
            first,
            last,
            periodic: basis.is_periodic(),
            period: basis.last_parameter() - basis.first_parameter(),
        },
        None => AdaptorRange {
            first: curve.first_parameter(),
            last: curve.last_parameter(),
            periodic: curve.is_periodic(),
            period: curve.last_parameter() - curve.first_parameter(),
        },
    }
}

impl NSections {
    /// `GeomFill_NSections::IsConstant` (`cxx:876-938`).
    pub(super) fn constant_section_test(&self) -> bool {
        // "on se limite a 2 sections" (`cxx:884`).
        if self.sections.len() != 2 {
            return false;
        }
        let first = self.sections[0].as_ref();
        let second = self.sections[1].as_ref();
        let kind = CurveKind::of(first);
        // "les sections doivent avoir le meme type" (`cxx:892-893`).
        if CurveKind::of(second) != kind {
            return false;
        }
        let tol = 1.0e-7;
        match kind {
            CurveKind::Circle => {
                let c1 = first.gp_circ().expect("Circle kind without a gp_Circ");
                let c2 = second.gp_circ().expect("Circle kind without a gp_Circ");
                let same_dir = axes_parallel(&c1.axis(), &c2.axis(), 1.0e-4);
                let same_radius = (c1.radius() - c2.radius()).abs() < tol;
                let mut same_pos = c1.location().distance(&c2.location()) < tol;
                if !same_pos {
                    // `gp_Ax1 D(C1.Location(), gp_Vec(C1.Location(),
                    // C2.Location()))` (`cxx:908-909`).
                    let v = c2.location().coord.subtracted(&c1.location().coord);
                    if let Ok(d) = GpDir::from_vec(&GpVec::from_xyz(&v)) {
                        same_pos = c1.axis().direction().is_parallel_tol(&d, 1.0e-4);
                    }
                }
                same_dir && same_radius && same_pos
            }
            CurveKind::Line => {
                let l1 = first.gp_line().expect("Line kind without a gp_Lin");
                let l2 = second.gp_line().expect("Line kind without a gp_Lin");
                let same_dir = l1.direction().is_parallel_tol(&l2.direction(), 1.0e-4);
                let p11 = first.d0(first.first_parameter());
                let p12 = first.d0(first.last_parameter());
                let p21 = second.d0(second.first_parameter());
                let p22 = second.d0(second.last_parameter());
                let same_length =
                    (p11.distance(&p12) - p21.distance(&p22)).abs() < tol;
                // "l'ecart entre les 2 sections ne compte pas" (`cxx:926`).
                let same_pos = (p11.distance(&p21) < tol && p12.distance(&p22) < tol)
                    || (p12.distance(&p21) < tol && p11.distance(&p22) < tol);
                same_dir && same_length && same_pos
            }
            CurveKind::Other => false,
        }
    }

    /// `GeomFill_NSections::ConstantSection` (`cxx:942-950`): a copy of the
    /// first section.
    pub(super) fn constant_section_copy(&self) -> Result<Arc<dyn Curve>, String> {
        let first = self
            .sections
            .first()
            .ok_or_else(|| "Standard_OutOfRange: GeomFill_NSections::ConstantSection".to_string())?;
        Ok(Arc::clone(first))
    }

    /// `GeomFill_NSections::IsConicalLaw` (`cxx:954-1011`).
    pub(super) fn conical_law_test(&self) -> bool {
        if self.sections.len() != 2 {
            return false;
        }
        let (Some(mut c1), Some(mut c2)) = (
            self.sections[0].gp_circ(),
            self.sections[1].gp_circ(),
        ) else {
            return false;
        };
        // `cxx:966-974`: the stored transforms are undone before the test.
        if let Some(trsf) = self.trsfs.first() {
            match trsf.inverted() {
                Ok(inv) => c1.transform(&inv),
                Err(_) => return false,
            }
        }
        if let Some(trsf) = self.trsfs.get(1) {
            match trsf.inverted() {
                Ok(inv) => c2.transform(&inv),
                Err(_) => return false,
            }
        }
        let tol = 1.0e-7;
        let mut is_conic = axes_parallel(&c1.axis(), &c2.axis(), 1.0e-4);
        if is_conic {
            // `cxx:988-989`: `gp_Lin Line1(C1.Axis())` and the distance from the
            // second centre to that line.
            let line1 = GpLin::new(c1.axis());
            is_conic = line1.distance(&c2.location()) < tol;
        }
        if is_conic {
            // `cxx:1000-1006` (OCC20866): both sections must span the same
            // parameter interval, measured on the `GeomAdaptor_Curve` view.
            let r1 = adaptor_range(self.sections[0].as_ref());
            let r2 = adaptor_range(self.sections[1].as_ref());
            is_conic = (r1.first - r2.first).abs() <= PCONFUSION
                && (r1.last - r2.last).abs() <= PCONFUSION;
        }
        is_conic
    }

    /// `GeomFill_NSections::CirclSection` (`cxx:1016-1047`).
    pub(super) fn circle_section_at(&self, v: f64) -> Result<Arc<dyn Curve>, String> {
        if !self.conical_law_test() {
            return Err("Standard_DomainError: GeomFill_NSections::CirclSection".to_string());
        }
        let first = self.sections.first().ok_or_else(|| {
            "Standard_OutOfRange: GeomFill_NSections::CirclSection".to_string()
        })?;
        let last = self.sections.last().ok_or_else(|| {
            "Standard_OutOfRange: GeomFill_NSections::CirclSection".to_string()
        })?;
        let mut c1: GpCirc = first
            .gp_circ()
            .ok_or_else(|| "Standard_DomainError: GeomFill_NSections::CirclSection".to_string())?;
        let c2 = last
            .gp_circ()
            .ok_or_else(|| "Standard_DomainError: GeomFill_NSections::CirclSection".to_string())?;

        // `cxx:1031-1033`: linear radius evolution between the two parameters.
        let p1 = *self.params.first().ok_or_else(|| {
            "Standard_OutOfRange: GeomFill_NSections::CirclSection".to_string()
        })?;
        let p2 = *self.params.last().ok_or_else(|| {
            "Standard_OutOfRange: GeomFill_NSections::CirclSection".to_string()
        })?;
        let radius = (c2.radius() - c1.radius()) * (v - p1) / (p2 - p1) + c1.radius();
        c1.set_radius(radius)
            .map_err(|_| "Standard_ConstructionError: gp_Circ::SetRadius".to_string())?;

        let circle: Arc<dyn Curve> = Arc::new(GeomCircle::new(c1));
        // `cxx:1036-1044`: trim to the first section's adaptor range unless
        // that range is a full period.
        let range = adaptor_range(first.as_ref());
        let a_par_f = range.first;
        let a_par_l = range.last;
        let a_period = if range.periodic { range.period } else { 0.0 };
        if a_period == 0.0 || (a_par_l - a_par_f - a_period).abs() > PCONFUSION {
            return Ok(Arc::new(GeomTrimmedCurve::new(circle, a_par_f, a_par_l)));
        }
        Ok(circle)
    }
}
