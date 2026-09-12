//! Face UV domain correction used by `IntTools_FaceFace::Perform`.
//!
//! Source: `IntTools_FaceFace.cxx` — `CorrectSurfaceBoundaries` (2017),
//! `CorrectPlaneBoundaries` (3093), `ParameterOutOfBoundary` (2213).

use std::sync::Arc;

use occt_core::bnd::BndBox2d;
use occt_core::gp::{GpPnt2d, GpTrsf, GpVec};
use occt_core::precision::{ANGULAR, CONFUSION, Precision};
use occt_geom::{Curve, Surface};
use occt_geom::geom_api::project_point_on_surface;

use crate::abs::{Orientation, ShapeType};
use crate::bop_split_seam::is_closed_on_face;
use crate::boptools_2d::curve_on_surface;
use crate::brep_surface::{classify_surface, face_uv_bounds, SurfaceKind};
use crate::brep_tool::BRepTool;
use crate::fclass2d::{FClass2d, FaceState};
use crate::iterator::ShapeIterator;
use crate::pcurve::pc_curve_kind;
use crate::pcurve::CurveKind as PcKind;
use crate::shape::{Edge, Face};
use crate::tgeometry::GeometryRegistry;

/// Adaptor that reports a restricted UV box (`GeomAdaptor_Surface::Load`).
#[derive(Clone)]
pub struct BoundedSurface {
    inner: Arc<dyn Surface>,
    umin: f64,
    umax: f64,
    vmin: f64,
    vmax: f64,
}

impl BoundedSurface {
    /// Wrap `inner` with `[umin, umax] x [vmin, vmax]`.
    pub fn new(inner: Arc<dyn Surface>, umin: f64, umax: f64, vmin: f64, vmax: f64) -> Self {
        Self {
            inner,
            umin,
            umax,
            vmin,
            vmax,
        }
    }
}

impl Surface for BoundedSurface {
    fn d0(&self, u: f64, v: f64) -> occt_core::gp::GpPnt {
        self.inner.d0(u, v)
    }
    fn d1(&self, u: f64, v: f64) -> (occt_core::gp::GpPnt, GpVec, GpVec) {
        self.inner.d1(u, v)
    }
    fn u_range(&self) -> (f64, f64) {
        (self.umin, self.umax)
    }
    fn v_range(&self) -> (f64, f64) {
        (self.vmin, self.vmax)
    }
    fn is_u_periodic(&self) -> bool {
        self.inner.is_u_periodic()
    }
    fn is_v_periodic(&self) -> bool {
        self.inner.is_v_periodic()
    }
    fn continuity(&self) -> u8 {
        self.inner.continuity()
    }
    fn transform(&mut self, _t: &GpTrsf) {}
    fn clone_dyn(&self) -> Box<dyn Surface> {
        Box::new(self.clone())
    }
}

/// `CorrectPlaneBoundaries` — grow a finite UV box by 10% on each side.
pub fn correct_plane_boundaries(
    a_umin: &mut f64,
    a_umax: &mut f64,
    a_vmin: &mut f64,
    a_vmax: &mut f64,
) {
    if !(Precision::is_infinite(*a_umin) || Precision::is_infinite(*a_umax)) {
        let d_u = 0.1 * (*a_umax - *a_umin);
        *a_umin -= d_u;
        *a_umax += d_u;
    }
    if !(Precision::is_infinite(*a_vmin) || Precision::is_infinite(*a_vmax)) {
        let d_v = 0.1 * (*a_vmax - *a_vmin);
        *a_vmin -= d_v;
        *a_vmax += d_v;
    }
}

/// Whether the surface kind is enlarged by `CorrectSurfaceBoundaries`.
fn enlarge_kind(kind: SurfaceKind) -> bool {
    matches!(
        kind,
        SurfaceKind::Cylinder | SurfaceKind::Other
    )
}

/// `CorrectSurfaceBoundaries` (`IntTools_FaceFace.cxx:2017`).
///
/// Enlarges a non-periodic UV box by `the_tolerance` toward the surface
/// natural bounds for cylinder / freeform faces. Periodic faces with linear
/// seam pcurves shrink the box to the seam AABB.
pub fn correct_surface_boundaries(
    the_face: &Face,
    the_tolerance: f64,
    theumin: &mut f64,
    theumax: &mut f64,
    thevmin: &mut f64,
    thevmax: &mut f64,
) {
    let Some(a_surface) = BRepTool::face_surface(the_face) else {
        return;
    };
    let (uinf, usup) = a_surface.u_range();
    let (vinf, vsup) = a_surface.v_range();
    let delta = the_tolerance;
    let kind = classify_surface(a_surface.as_ref());
    let enlarge = enlarge_kind(kind);
    let isuperiodic = a_surface.is_u_periodic();
    let isvperiodic = a_surface.is_v_periodic();

    if !isuperiodic && enlarge {
        if !Precision::is_infinite(*theumin) && (*theumin - uinf) > delta {
            *theumin -= delta;
        } else {
            *theumin = uinf;
        }
        if !Precision::is_infinite(*theumax) && (usup - *theumax) > delta {
            *theumax += delta;
        } else {
            *theumax = usup;
        }
    }
    if !isvperiodic && enlarge {
        if !Precision::is_infinite(*thevmin) && (*thevmin - vinf) > delta {
            *thevmin -= delta;
        } else {
            *thevmin = vinf;
        }
        if !Precision::is_infinite(*thevmax) && (vsup - *thevmax) > delta {
            *thevmax += delta;
        } else {
            *thevmax = vsup;
        }
    }

    if isuperiodic || isvperiodic {
        let mut correct = false;
        let mut correct_u = false;
        let mut correct_v = false;
        let mut a_box = BndBox2d::new();
        let an_u_dir = occt_core::gp::GpDir2d::new(1.0, 0.0).unwrap_or_default();
        let a_v_dir = occt_core::gp::GpDir2d::new(0.0, 1.0).unwrap_or_default();
        let an_angular_tolerance = ANGULAR;

        'edges: for e in ShapeIterator::of_shape(&the_face.0) {
            if e.shape_type() != ShapeType::Edge {
                continue;
            }
            let mut an_edge = Edge(e.clone());
            if !is_closed_on_face(&an_edge, the_face) {
                continue;
            }
            correct = true;
            for i in 0..2 {
                if i == 0 {
                    an_edge.0.set_orientation(Orientation::Forward);
                } else {
                    an_edge.0.set_orientation(Orientation::Reversed);
                }
                let Some(a_curve) = curve_on_surface(&an_edge, the_face) else {
                    correct = false;
                    break 'edges;
                };
                if pc_curve_kind(a_curve.as_ref()) != PcKind::Line {
                    correct = false;
                    break 'edges;
                }
                let (_, d1) = a_curve.d1(0.0);
                let Ok(dir) = occt_core::gp::GpDir2d::from_vec2d(&d1) else {
                    correct = false;
                    break 'edges;
                };
                correct_u = correct_u || dir.is_parallel(&a_v_dir, an_angular_tolerance);
                correct_v = correct_v || dir.is_parallel(&an_u_dir, an_angular_tolerance);
                let (f, l) = GeometryRegistry::global().edge_parameters(&an_edge.0);
                a_box.add_point(&a_curve.d0(f));
                a_box.add_point(&a_curve.d0(l));
            }
        }

        if correct {
            if let Some((umin, umax, vmin, vmax)) = a_box.get() {
                if isuperiodic && correct_u {
                    if *theumin < umin {
                        *theumin = umin;
                    }
                    if *theumax > umax {
                        *theumax = umax;
                    }
                }
                if isvperiodic && correct_v {
                    if *thevmin < vmin {
                        *thevmin = vmin;
                    }
                    if *thevmax > vmax {
                        *thevmax = vmax;
                    }
                }
            }
        }
    }
}

/// Load face UV, then apply `CorrectSurfaceBoundaries` / `CorrectPlaneBoundaries`
/// the same way `IntTools_FaceFace::Perform` loads `GeomAdaptor_Surface`.
pub fn corrected_uv_box(
    face: &Face,
    other_kind: SurfaceKind,
    self_kind: SurfaceKind,
    the_tol: f64,
) -> (f64, f64, f64, f64) {
    let (mut umin, mut umax, mut vmin, mut vmax) = face_uv_bounds(face);
    let is_self_plane = self_kind == SurfaceKind::Plane;
    let is_other_quad = matches!(
        other_kind,
        SurfaceKind::Cylinder | SurfaceKind::Cone | SurfaceKind::Torus
    );
    if is_self_plane && is_other_quad {
        correct_plane_boundaries(&mut umin, &mut umax, &mut vmin, &mut vmax);
    } else {
        correct_surface_boundaries(face, the_tol * 2.0, &mut umin, &mut umax, &mut vmin, &mut vmax);
    }
    (umin, umax, vmin, vmax)
}

/// `ParameterOutOfBoundary` (`IntTools_FaceFace.cxx:2213`).
///
/// Walks `the_parameter` off the ON state of either face, returning a new
/// parameter that still lies before `the_other_parameter`.
pub fn parameter_out_of_boundary(
    the_parameter: f64,
    the_curve: &dyn Curve,
    the_face1: &Face,
    the_face2: &Face,
    the_other_parameter: f64,
    b_increase_par: bool,
    the_tol: f64,
) -> Option<f64> {
    let mut the_new_parameter = the_parameter;
    let mut acurpar = the_parameter;
    let mut a_state = FaceState::On;
    let mut iter = 0i32;
    let asumtol = the_tol;
    let mut adelta = asumtol * 0.1;
    adelta = if adelta < CONFUSION { CONFUSION } else { adelta };
    let Some(a_surf1) = BRepTool::face_surface(the_face1) else {
        return None;
    };
    let Some(a_surf2) = BRepTool::face_surface(the_face2) else {
        return None;
    };
    let Ok(class1) = FClass2d::new(the_face1, BRepTool::face_tolerance(the_face1)) else {
        return None;
    };
    let Ok(class2) = FClass2d::new(the_face2, BRepTool::face_tolerance(the_face2)) else {
        return None;
    };

    while a_state == FaceState::On {
        if b_increase_par {
            acurpar += adelta;
        } else {
            acurpar -= adelta;
        }
        let a_p_current = the_curve.d0(acurpar);
        if let Some(proj) = project_point_on_surface(a_surf1.as_ref(), &a_p_current, 0.0) {
            a_state = class1.perform(GpPnt2d::new(proj.u, proj.v));
        }
        if a_state != FaceState::On {
            if let Some(proj) = project_point_on_surface(a_surf2.as_ref(), &a_p_current, 0.0) {
                a_state = class2.perform(GpPnt2d::new(proj.u, proj.v));
            }
        }
        if iter > 11 {
            break;
        }
        iter += 1;
    }

    if iter <= 11 {
        the_new_parameter = acurpar;
        if b_increase_par {
            if acurpar >= the_other_parameter {
                the_new_parameter = the_other_parameter;
            }
        } else if acurpar <= the_other_parameter {
            the_new_parameter = the_other_parameter;
        }
        Some(the_new_parameter)
    } else {
        None
    }
}
