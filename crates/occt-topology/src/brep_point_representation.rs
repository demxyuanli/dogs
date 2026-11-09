//! `BRep_PointRepresentation` and the two vertex-point subclasses that carry
//! a surface - the `BRep_TVertex::Points()` list entries.
//!
//! OCCT models the family with a root class plus four concrete classes
//! (`BRep_PointOnCurve`, `BRep_PointOnCurveOnSurface`, `BRep_PointOnSurface`,
//! `BRep_PointsOnSurface` as the surface-carrying base,
//! `BRep_PointRepresentation.cxx:26-125`, `BRep_PointsOnSurface.cxx:26-50`,
//! `BRep_PointOnCurve.cxx:29-60`, `BRep_PointOnCurveOnSurface.cxx:32-62`,
//! `BRep_PointOnSurface.cxx:29-66`). The port answers the same queries with one
//! enum: every virtual that the root answers `false` for maps to a variant arm
//! answering `false`, and every virtual the root *throws* `Standard_DomainError`
//! for maps to an `Option` answering `None`.
//!
//! UNPORTED: `TopLoc_Location`. Every `IsPointOn*` comparison in OCCT is
//! `(handle == handle) && (Location() == L)`; the port keeps no locations on
//! the representations (`BRep_Tool::Surface` carries none either, see
//! `brep_tool::BRepTool::face_surface`), so the handle comparison is the whole
//! predicate.

use std::sync::Arc;

use occt_geom::{Curve, Surface};
use occt_geom2d::curve::Curve2d;

/// One entry of `BRep_TVertex::myPoints` (`BRep_TVertex.hxx:96`).
#[derive(Clone)]
pub enum PointRepresentation {
    /// `BRep_PointOnCurve(P, C, L)` (`BRep_PointOnCurve.cxx:30-38`).
    OnCurve { curve: Arc<dyn Curve>, parameter: f64 },
    /// `BRep_PointOnCurveOnSurface(P, C, S, L)`
    /// (`BRep_PointOnCurveOnSurface.cxx:31-40`).
    OnCurveOnSurface { pcurve: Arc<dyn Curve2d>, surface: Arc<dyn Surface>, parameter: f64 },
    /// `BRep_PointOnSurface(P1, P2, S, L)` (`BRep_PointOnSurface.cxx:30-39`).
    OnSurface { surface: Arc<dyn Surface>, parameter: f64, parameter2: f64 },
}

/// Handle identity plus the stored scalars - the fields OCCT compares in its
/// copy / `IsEqual` style predicates. `TopLoc_Location` is unported, so it
/// takes no part (see the module note).
impl PartialEq for PointRepresentation {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::OnCurve { curve: c1, parameter: p1 }, Self::OnCurve { curve: c2, parameter: p2 }) => {
                Arc::ptr_eq(c1, c2) && p1 == p2
            }
            (
                Self::OnCurveOnSurface { pcurve: pc1, surface: s1, parameter: p1 },
                Self::OnCurveOnSurface { pcurve: pc2, surface: s2, parameter: p2 },
            ) => Arc::ptr_eq(pc1, pc2) && Arc::ptr_eq(s1, s2) && p1 == p2,
            (
                Self::OnSurface { surface: s1, parameter: p1, parameter2: q1 },
                Self::OnSurface { surface: s2, parameter: p2, parameter2: q2 },
            ) => Arc::ptr_eq(s1, s2) && p1 == p2 && q1 == q2,
            _ => false,
        }
    }
}

impl PointRepresentation {
    /// `BRep_PointRepresentation::Parameter()` (`BRep_PointRepresentation.lxx:14`).
    pub fn parameter(&self) -> f64 {
        match self {
            Self::OnCurve { parameter, .. }
            | Self::OnCurveOnSurface { parameter, .. }
            | Self::OnSurface { parameter, .. } => *parameter,
        }
    }

    /// `BRep_PointRepresentation::Parameter(P)` (`BRep_PointRepresentation.lxx:20`).
    pub fn set_parameter(&mut self, p: f64) {
        match self {
            Self::OnCurve { parameter, .. }
            | Self::OnCurveOnSurface { parameter, .. }
            | Self::OnSurface { parameter, .. } => *parameter = p,
        }
    }

    /// `BRep_PointRepresentation::Parameter2()` (`BRep_PointRepresentation.cxx:86-89`
    /// throws, `BRep_PointOnSurface.cxx:58-61` answers the stored value).
    pub fn parameter2(&self) -> Option<f64> {
        match self {
            Self::OnSurface { parameter2, .. } => Some(*parameter2),
            Self::OnCurve { .. } | Self::OnCurveOnSurface { .. } => None,
        }
    }

    /// `BRep_PointRepresentation::Parameter2(P)`
    /// (`BRep_PointRepresentation.cxx:94-97`, `BRep_PointOnSurface.cxx:63-66`).
    pub fn set_parameter2(&mut self, p: f64) {
        if let Self::OnSurface { parameter2, .. } = self {
            *parameter2 = p;
        }
    }

    /// `BRep_PointRepresentation::Curve()` (`BRep_PointRepresentation.cxx:101-104`
    /// throws, `BRep_PointOnCurve.cxx:51-54`).
    pub fn curve(&self) -> Option<&Arc<dyn Curve>> {
        match self {
            Self::OnCurve { curve, .. } => Some(curve),
            _ => None,
        }
    }

    /// `BRep_PointRepresentation::Curve(C)` (`BRep_PointRepresentation.cxx:107-110`,
    /// `BRep_PointOnCurve.cxx:57-60`).
    pub fn set_curve(&mut self, c: Arc<dyn Curve>) {
        if let Self::OnCurve { curve, .. } = self {
            *curve = c;
        }
    }

    /// `BRep_PointRepresentation::PCurve()`
    /// (`BRep_PointRepresentation.cxx:113-116` throws,
    /// `BRep_PointOnCurveOnSurface.cxx:51-54`).
    pub fn pcurve(&self) -> Option<&Arc<dyn Curve2d>> {
        match self {
            Self::OnCurveOnSurface { pcurve, .. } => Some(pcurve),
            _ => None,
        }
    }

    /// `BRep_PointRepresentation::PCurve(C)`
    /// (`BRep_PointRepresentation.cxx:119-122`,
    /// `BRep_PointOnCurveOnSurface.cxx:57-60`).
    pub fn set_pcurve(&mut self, pc: Arc<dyn Curve2d>) {
        if let Self::OnCurveOnSurface { pcurve, .. } = self {
            *pcurve = pc;
        }
    }

    /// `BRep_PointRepresentation::Surface()` (`BRep_PointRepresentation.cxx:125-128`
    /// throws, `BRep_PointsOnSurface.cxx:35-38`).
    pub fn surface(&self) -> Option<&Arc<dyn Surface>> {
        match self {
            Self::OnCurveOnSurface { surface, .. } | Self::OnSurface { surface, .. } => Some(surface),
            Self::OnCurve { .. } => None,
        }
    }

    /// `BRep_PointRepresentation::Surface(S)` (`BRep_PointRepresentation.cxx:131-134`,
    /// `BRep_PointsOnSurface.cxx:41-44`).
    pub fn set_surface(&mut self, s: Arc<dyn Surface>) {
        match self {
            Self::OnCurveOnSurface { surface, .. } | Self::OnSurface { surface, .. } => *surface = s,
            Self::OnCurve { .. } => {}
        }
    }

    /// `BRep_PointRepresentation::IsPointOnCurve()`
    /// (`BRep_PointRepresentation.cxx:40-43`, `BRep_PointOnCurve.cxx:43-46`).
    pub fn is_point_on_curve(&self) -> bool {
        matches!(self, Self::OnCurve { .. })
    }

    /// `BRep_PointRepresentation::IsPointOnCurveOnSurface()`
    /// (`BRep_PointRepresentation.cxx:47-50`,
    /// `BRep_PointOnCurveOnSurface.cxx:43-46`).
    pub fn is_point_on_curve_on_surface(&self) -> bool {
        matches!(self, Self::OnCurveOnSurface { .. })
    }

    /// `BRep_PointRepresentation::IsPointOnSurface()`
    /// (`BRep_PointRepresentation.cxx:53-56`, `BRep_PointOnSurface.cxx:50-53`).
    pub fn is_point_on_surface(&self) -> bool {
        matches!(self, Self::OnSurface { .. })
    }

    /// `IsPointOnCurve(C, L)` (`BRep_PointRepresentation.cxx:59-63`,
    /// `BRep_PointOnCurve.cxx:49-55`): handle identity plus the location
    /// comparison the port does not carry.
    pub fn is_point_on_curve_of(&self, c: &Arc<dyn Curve>) -> bool {
        match self {
            Self::OnCurve { curve, .. } => Arc::ptr_eq(curve, c),
            _ => false,
        }
    }

    /// `IsPointOnCurveOnSurface(PC, S, L)`
    /// (`BRep_PointRepresentation.cxx:66-71`,
    /// `BRep_PointOnCurveOnSurface.cxx:49-55`).
    pub fn is_point_on_curve_on_surface_of(
        &self,
        pc: &Arc<dyn Curve2d>,
        s: &Arc<dyn Surface>,
    ) -> bool {
        match self {
            Self::OnCurveOnSurface { pcurve, surface, .. } => {
                Arc::ptr_eq(pcurve, pc) && Arc::ptr_eq(surface, s)
            }
            _ => false,
        }
    }

    /// `IsPointOnSurface(S, L)` (`BRep_PointRepresentation.cxx:74-78`,
    /// `BRep_PointOnSurface.cxx:56-62`).
    pub fn is_point_on_surface_of(&self, s: &Arc<dyn Surface>) -> bool {
        match self {
            Self::OnSurface { surface, .. } => Arc::ptr_eq(surface, s),
            _ => false,
        }
    }
}
