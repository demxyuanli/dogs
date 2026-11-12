//! `BRepCheck_Edge::InContext`.
//!
//! Source: `BRepCheck_Edge.cxx:263-556`.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpTrsf, GpVec};
use occt_core::precision::PCONFUSION;
use occt_geom::curve::Curve;
use occt_geom::curve_on_surface::CurveOnSurface;
use occt_geom::surface::Surface;
use occt_geom2d::curve::Curve2d;

use crate::abs::ShapeType;
use crate::bop_occt_util::{explore, iter_children, same_shape};
use crate::brep_lib_validate_edge::ValidateEdge;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape};
use crate::tgeometry::GeometryRegistry;

use super::{basis_window_2d, basis_window_3d, Status};

enum ReferenceCurve {
    Curve3d(Arc<dyn Curve>),
    CurveOnSurface {
        curve: Arc<dyn Curve>,
        surface: Arc<dyn Surface>,
        pcurve: Arc<dyn Curve2d>,
        location: occt_core::toploc::TopLocLocation,
    },
}

#[derive(Clone)]
struct RangedCurve {
    curve: Arc<dyn Curve>,
    first: f64,
    last: f64,
}

impl RangedCurve {
    fn new(curve: Arc<dyn Curve>, first: f64, last: f64) -> Self {
        Self { curve, first, last }
    }
}

impl Curve for RangedCurve {
    fn d0(&self, u: f64) -> GpPnt {
        self.curve.d0(u)
    }

    fn d1(&self, u: f64) -> (GpPnt, GpVec) {
        self.curve.d1(u)
    }

    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        self.curve.d2(u)
    }

    fn first_parameter(&self) -> f64 {
        self.first
    }

    fn last_parameter(&self) -> f64 {
        self.last
    }

    fn is_periodic(&self) -> bool {
        self.curve.is_periodic()
    }

    fn period(&self) -> f64 {
        self.curve.period()
    }

    fn continuity(&self) -> u8 {
        self.curve.continuity()
    }

    fn gp_line(&self) -> Option<occt_core::gp::GpLin> {
        self.curve.gp_line()
    }

    fn gp_circ(&self) -> Option<occt_core::gp::GpCirc> {
        self.curve.gp_circ()
    }

    fn gp_ellipse(&self) -> Option<occt_core::gp::GpElips> {
        self.curve.gp_ellipse()
    }

    fn gp_hyperbola(&self) -> Option<occt_core::gp::GpHypr> {
        self.curve.gp_hyperbola()
    }

    fn gp_parabola(&self) -> Option<occt_core::gp::GpParab> {
        self.curve.gp_parabola()
    }

    fn resolution(&self, r3d: f64) -> f64 {
        self.curve.resolution(r3d)
    }

    fn transform(&mut self, trsf: &GpTrsf) {
        self.curve = Arc::from(self.curve.transformed(trsf));
    }

    fn reverse(&mut self) {
        self.curve = Arc::from(self.curve.reversed());
    }

    fn clone_dyn(&self) -> Box<dyn Curve> {
        Box::new(self.clone())
    }
}

impl ReferenceCurve {
    fn curve(&self) -> &dyn Curve {
        match self {
            Self::Curve3d(curve) => curve.as_ref(),
            Self::CurveOnSurface { curve, .. } => curve.as_ref(),
        }
    }

    fn is_representation(
        &self,
        surface: &Arc<dyn Surface>,
        pcurve: &Arc<dyn Curve2d>,
        location: &occt_core::toploc::TopLocLocation,
    ) -> bool {
        match self {
            Self::Curve3d(_) => false,
            Self::CurveOnSurface {
                surface: reference_surface,
                pcurve: reference_pcurve,
                location: reference_location,
                ..
            } => {
                Arc::ptr_eq(reference_surface, surface)
                    && Arc::ptr_eq(reference_pcurve, pcurve)
                    && reference_location.is_equal(location)
            }
        }
    }
}

pub(super) fn edge_in_context(edge: &Edge, context: &TopoShape) -> Vec<Status> {
    if !explore(context, ShapeType::Edge)
        .iter()
        .any(|candidate| same_shape(candidate, &edge.0))
    {
        return vec![Status::SubshapeNotInShape];
    }

    match context.shape_type() {
        ShapeType::Wire => vec![Status::NoError],
        ShapeType::Face => edge_in_face(edge, &Face(context.clone())),
        ShapeType::Solid => edge_in_solid(edge, context),
        _ => vec![Status::NoError],
    }
}

fn edge_in_face(edge: &Edge, face: &Face) -> Vec<Status> {
    let registry = GeometryRegistry::global();
    let Some((reference, first, last)) = reference_curve(edge) else {
        return vec![Status::NoError];
    };

    let same_parameter = registry.same_parameter(&edge.0);
    let same_range = registry.same_range(&edge.0);
    if !same_parameter || !same_range {
        let mut statuses = Vec::new();
        if !same_parameter {
            add_status(&mut statuses, Status::InvalidSameParameterFlag);
        }
        if !same_range {
            add_status(&mut statuses, Status::InvalidSameRangeFlag);
        }
        return statuses;
    }

    let Some(face_surface) = BRepTool::face_surface(face) else {
        return vec![Status::NoCurveOnSurface];
    };
    let transformed_surface: Arc<dyn Surface> = if face.0.location().is_identity() {
        face_surface.clone()
    } else {
        Arc::from(face_surface.transformed(&face.0.location().transformation()))
    };
    let tolerance = BRepTool::edge_tolerance(edge);
    let representations = registry.edge_pcurve_reps_with_location(&edge.0);
    let wanted_location = face.0.location().predivided(edge.0.location());
    let mut statuses = Vec::new();
    let mut pcurve_found = false;
    let mut i = 0;

    while i < representations.len() {
        let (surface, pcurve, range_first, range_last, location) = &representations[i];
        let has_second = representations
            .get(i + 1)
            .is_some_and(|next| Arc::ptr_eq(&next.0, surface) && next.4.is_equal(location));
        let step = if has_second { 2 } else { 1 };
        i += step;

        if !Arc::ptr_eq(surface, &face_surface)
            || !location.is_equal(&wanted_location)
            || reference.is_representation(surface, pcurve, location)
        {
            continue;
        }
        pcurve_found = true;

        if (range_first - first).abs() > PCONFUSION || (range_last - last).abs() > PCONFUSION {
            add_status(&mut statuses, Status::InvalidSameRangeFlag);
            add_status(&mut statuses, Status::InvalidSameParameterFlag);
        }

        let (periodic, period, curve_first, curve_last) = basis_window_2d(pcurve.as_ref());
        if periodic && range_last - range_first > period + PCONFUSION {
            add_status(&mut statuses, Status::InvalidRange);
            return statuses;
        }
        if !periodic
            && (*range_first < curve_first - PCONFUSION || *range_last > curve_last + PCONFUSION)
        {
            add_status(&mut statuses, Status::InvalidRange);
            return statuses;
        }

        validate_representation(
            &mut statuses,
            reference.curve(),
            pcurve,
            &transformed_surface,
            *range_first,
            *range_last,
            same_parameter,
            tolerance,
            if has_second {
                Status::InvalidCurveOnClosedSurface
            } else {
                Status::InvalidCurveOnSurface
            },
        );

        if has_second {
            let second = &representations[i - 1].1;
            validate_representation(
                &mut statuses,
                reference.curve(),
                second,
                &transformed_surface,
                *range_first,
                *range_last,
                same_parameter,
                tolerance,
                Status::InvalidCurveOnClosedSurface,
            );
        }
    }

    if !pcurve_found {
        if face_surface.gp_pln().is_none() {
            add_status(&mut statuses, Status::NoCurveOnSurface);
        } else if let Some((projected, projected_first, projected_last)) =
            crate::boptools_2d::curve_on_surface_oriented(edge, face, false)
        {
            validate_representation(
                &mut statuses,
                reference.curve(),
                &projected,
                &transformed_surface,
                projected_first,
                projected_last,
                same_parameter,
                tolerance,
                Status::InvalidCurveOnSurface,
            );
        }
    }

    if statuses.is_empty() {
        statuses.push(Status::NoError);
    }
    statuses
}

fn reference_curve(edge: &Edge) -> Option<(ReferenceCurve, f64, f64)> {
    let registry = GeometryRegistry::global();
    if let Some(curve) = BRepTool::edge_curve(edge) {
        let (first, last) = BRepTool::edge_parameters(edge);
        let (periodic, period, curve_first, curve_last) = basis_window_3d(curve.as_ref());
        if !reference_range_is_valid(first, last, periodic, period, curve_first, curve_last) {
            return None;
        }
        let curve = if edge.0.location().is_identity() {
            curve
        } else {
            Arc::from(curve.transformed(&edge.0.location().transformation()))
        };
        let ranged: Arc<dyn Curve> = Arc::new(RangedCurve::new(curve, first, last));
        return Some((ReferenceCurve::Curve3d(ranged), first, last));
    }
    if registry.is_degenerated_edge(&edge.0) {
        return None;
    }
    let (surface, pcurve, first, last, location) = registry
        .edge_pcurve_reps_with_location(&edge.0)
        .into_iter()
        .next()?;
    let (periodic, period, curve_first, curve_last) = basis_window_2d(pcurve.as_ref());
    if !reference_range_is_valid(first, last, periodic, period, curve_first, curve_last) {
        return None;
    }
    let surface_for_evaluation = if location.is_identity() {
        surface.clone()
    } else {
        Arc::from(surface.transformed(&location.transformation()))
    };
    let curve_on_surface: Arc<dyn Curve> =
        Arc::new(CurveOnSurface::new(pcurve.clone(), surface_for_evaluation));
    let curve = Arc::new(RangedCurve::new(curve_on_surface, first, last));
    Some((
        ReferenceCurve::CurveOnSurface {
            curve,
            surface,
            pcurve,
            location,
        },
        first,
        last,
    ))
}

fn reference_range_is_valid(
    first: f64,
    last: f64,
    periodic: bool,
    period: f64,
    curve_first: f64,
    curve_last: f64,
) -> bool {
    if last <= first {
        return false;
    }
    if periodic {
        last - first <= period + PCONFUSION
    } else {
        first >= curve_first - PCONFUSION && last <= curve_last + PCONFUSION
    }
}

fn validate_representation(
    statuses: &mut Vec<Status>,
    reference: &dyn Curve,
    pcurve: &Arc<dyn Curve2d>,
    surface: &Arc<dyn Surface>,
    first: f64,
    last: f64,
    same_parameter: bool,
    tolerance: f64,
    invalid_curve_status: Status,
) {
    let other_curve: Arc<dyn Curve> =
        Arc::new(CurveOnSurface::new(pcurve.clone(), surface.clone()));
    let other = RangedCurve::new(other_curve, first, last);
    let mut validation = ValidateEdge::new(
        reference,
        &other,
        Some(surface.as_ref()),
        Some(pcurve.as_ref()),
        same_parameter,
    );
    validation.set_exit_if_tolerance_exceeded(tolerance);
    validation.process();
    if !validation.is_done() || !validation.check_tolerance(tolerance) {
        add_status(statuses, invalid_curve_status);
        add_status(statuses, Status::InvalidSameParameterFlag);
    }
}

fn edge_in_solid(edge: &Edge, solid: &TopoShape) -> Vec<Status> {
    let mut connections = 0;
    for face in explore(solid, ShapeType::Face) {
        connections += count_edge_occurrences(&face, edge);
    }
    let degenerated = GeometryRegistry::global().is_degenerated_edge(&edge.0);
    if connections < 2 && !degenerated {
        vec![Status::FreeEdge]
    } else if connections > 2 {
        vec![Status::InvalidMultiConnexity]
    } else {
        vec![Status::NoError]
    }
}

fn count_edge_occurrences(shape: &TopoShape, edge: &Edge) -> usize {
    if shape.shape_type() == ShapeType::Edge {
        return usize::from(same_shape(shape, &edge.0));
    }
    iter_children(shape)
        .iter()
        .map(|child| count_edge_occurrences(child, edge))
        .sum()
}

fn add_status(statuses: &mut Vec<Status>, status: Status) {
    if !statuses.contains(&status) {
        statuses.push(status);
    }
}
