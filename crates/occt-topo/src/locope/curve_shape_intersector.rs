//! `LocOpe_CurveShapeIntersector` (`LocOpe_CurveShapeIntersector.hxx:38-131`,
//! `.cxx:33-429`).
//!
//! Intersects an axis (`Init(gp_Ax1, S)`, `cxx:33-48`) or a full circle
//! (`Init(gp_Circ, S)`, `cxx:52-70`) with every face of a shape through
//! [`BRepIntCurveSurfaceInter`] and keeps the hits sorted by parameter
//! (`cxx:357-428`). The four `Localize*` scans are shared with
//! [`super::cs_intersector::CSIntersector`] and live in
//! [`super::localize`]; this class always passes `Precision::Confusion()`
//! (`cxx:84`, `cxx:144`, `cxx:216`, `cxx:296`).

use std::sync::Arc;

use occt_core::gp::{GpAx1, GpCirc, GpLin};
use occt_core::precision::CONFUSION;
use occt_geom::{Curve, GeomCircle};

use crate::abs::Orientation;
use crate::brep_int_curve_surface_inter::BRepIntCurveSurfaceInter;
use crate::int_curves_face::Transition;
use crate::shape::TopoShape;

use super::localize::{self, Localized};
use super::pnt_face::PntFace;

pub struct CurveShapeIntersector {
    my_done: bool,
    my_points: Vec<PntFace>,
}

impl Default for CurveShapeIntersector {
    fn default() -> Self {
        Self::new()
    }
}

impl CurveShapeIntersector {
    /// `LocOpe_CurveShapeIntersector()` (`.lxx:23`): `myDone` false, empty
    /// sequence.
    pub fn new() -> Self {
        Self {
            my_done: false,
            my_points: Vec::new(),
        }
    }

    /// `LocOpe_CurveShapeIntersector(Axis, S)` (`.lxx:27-30`).
    pub fn with_axis(axis: &GpAx1, s: &TopoShape) -> Self {
        let mut i = Self::new();
        i.init_axis(axis, s);
        i
    }

    /// `LocOpe_CurveShapeIntersector(C, S)` (`.lxx:35-38`).
    pub fn with_circle(c: &GpCirc, s: &TopoShape) -> Self {
        let mut i = Self::new();
        i.init_circle(c, s);
        i
    }

    /// `Init(Axis, S)` (`.cxx:33-48`).
    pub fn init_axis(&mut self, axis: &GpAx1, s: &TopoShape) {
        self.my_done = false;
        self.my_points.clear();
        if s.is_null() {
            return;
        }
        let tol = CONFUSION;
        let mut the_int = BRepIntCurveSurfaceInter::new();
        the_int.init_shape_line(s, &GpLin::new(*axis), tol);
        self.my_points = perform(&mut the_int);
        self.my_done = true;
    }

    /// `Init(C, S)` (`.cxx:52-70`): the circle is wrapped in a
    /// `Geom_Circle` / `GeomAdaptor_Curve(C, 0, 2*pi)`.
    pub fn init_circle(&mut self, c: &GpCirc, s: &TopoShape) {
        self.my_done = false;
        self.my_points.clear();
        if s.is_null() {
            return;
        }
        let tol = CONFUSION;
        let gc: Arc<dyn Curve> = Arc::new(GeomCircle::new(c.clone()));
        let mut the_int = BRepIntCurveSurfaceInter::new();
        the_int.init_shape_curve(s, gc, tol);
        self.my_points = perform(&mut the_int);
        self.my_done = true;
    }

    /// `IsDone()` (`.lxx:43-46`).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// `NbPoints()` (`.lxx:50-53`).
    pub fn nb_points(&self) -> usize {
        self.my_points.len()
    }

    /// `Point(Index)` (`.lxx:57-60`), 1-based.
    pub fn point(&self, index: usize) -> Option<&PntFace> {
        if index == 0 {
            return None;
        }
        self.my_points.get(index - 1)
    }

    /// `LocalizeAfter(From, Or, IndFrom, IndTo)` (`.cxx:75-133`).
    pub fn localize_after_param(&self, from: f64) -> Option<Localized> {
        localize::loc_after_param(&self.my_points, from, CONFUSION)
    }

    /// `LocalizeBefore(From, Or, IndFrom, IndTo)` (`.cxx:136-197`).
    pub fn localize_before_param(&self, from: f64) -> Option<Localized> {
        localize::loc_before_param(&self.my_points, from, CONFUSION)
    }

    /// `LocalizeAfter(FromInd, Or, IndFrom, IndTo)` (`.cxx:201-275`).
    pub fn localize_after_index(&self, from_ind: usize) -> Option<Localized> {
        localize::loc_after_index(&self.my_points, from_ind, CONFUSION)
    }

    /// `LocalizeBefore(FromInd, Or, IndFrom, IndTo)` (`.cxx:279-354`).
    pub fn localize_before_index(&self, from_ind: usize) -> Option<Localized> {
        localize::loc_before_index(&self.my_points, from_ind, CONFUSION)
    }
}

/// The static `Perform(BRepIntCurveSurface_Inter&, NCollection_Sequence&)`
/// (`.cxx:357-428`): walks every hit, turns its transition and the hit face's
/// orientation into a `TopAbs_Orientation`, and inserts it before the first
/// recorded point with a greater parameter (appending when none is).
fn perform(the_int: &mut BRepIntCurveSurfaceInter) -> Vec<PntFace> {
    let mut points: Vec<PntFace> = Vec::new();
    while the_int.more() {
        let Some(face) = the_int.face().cloned() else {
            break;
        };
        let (Ok(pnt), Ok(param), Ok(paramu), Ok(paramv), Ok(transition)) = (
            the_int.pnt(),
            the_int.w(),
            the_int.u(),
            the_int.v(),
            the_int.transition(),
        ) else {
            break;
        };
        let orface = face.0.orientation();
        // `cxx:373-407`.
        let theor = match transition {
            Transition::In => match orface {
                Orientation::Forward => Orientation::Forward,
                Orientation::Reversed => Orientation::Reversed,
                _ => Orientation::External,
            },
            Transition::Out => match orface {
                Orientation::Forward => Orientation::Reversed,
                Orientation::Reversed => Orientation::Forward,
                _ => Orientation::External,
            },
            Transition::Tangent => Orientation::External,
        };
        let newpt = PntFace::new(pnt, face, theor, param, paramu, paramv);
        // `cxx:411-427`.
        let mut i = 0usize;
        while i < points.len() {
            if points[i].parameter() - param > 0.0 {
                break;
            }
            i += 1;
        }
        points.insert(i, newpt);
        the_int.next();
    }
    points
}
