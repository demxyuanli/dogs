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
    umin: &mut f64,
    umax: &mut f64,
    vmin: &mut f64,
    vmax: &mut f64,
) {
    if !(Precision::is_infinite(*umin) || Precision::is_infinite(*umax)) {
        let d_u = 0.1 * (*umax - *umin);
        *umin -= d_u;
        *umax += d_u;
    }
    if !(Precision::is_infinite(*vmin) || Precision::is_infinite(*vmax)) {
        let d_v = 0.1 * (*vmax - *vmin);
        *vmin -= d_v;
        *vmax += d_v;
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
/// Enlarges a non-periodic UV box by `tolerance` toward the surface
/// natural bounds for cylinder / freeform faces. Periodic faces with linear
/// seam pcurves shrink the box to the seam AABB.
pub fn correct_surface_boundaries(
    face: &Face,
    tolerance: f64,
    theumin: &mut f64,
    theumax: &mut f64,
    thevmin: &mut f64,
    thevmax: &mut f64,
) {
    let Some(surface) = BRepTool::face_surface(face) else {
        return;
    };
    let (uinf, usup) = surface.u_range();
    let (vinf, vsup) = surface.v_range();
    let delta = tolerance;
    let kind = classify_surface(surface.as_ref());
    let enlarge = enlarge_kind(kind);
    let isuperiodic = surface.is_u_periodic();
    let isvperiodic = surface.is_v_periodic();

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
        let mut bnd_box = BndBox2d::new();
        let an_u_dir = occt_core::gp::GpDir2d::new(1.0, 0.0).unwrap_or_default();
        let v_dir = occt_core::gp::GpDir2d::new(0.0, 1.0).unwrap_or_default();
        let an_angular_tolerance = ANGULAR;

        'edges: for e in ShapeIterator::of_shape(&face.0) {
            if e.shape_type() != ShapeType::Edge {
                continue;
            }
            let mut an_edge = Edge(e.clone());
            if !is_closed_on_face(&an_edge, face) {
                continue;
            }
            correct = true;
            for i in 0..2 {
                if i == 0 {
                    an_edge.0.set_orientation(Orientation::Forward);
                } else {
                    an_edge.0.set_orientation(Orientation::Reversed);
                }
                let Some(a_curve) = curve_on_surface(&an_edge, face) else {
                    correct = false;
                    break 'edges;
                };
                // `occ::down_cast<Geom2d_Line>(aCurve)` (`IntTools_FaceFace.cxx:2149`).
                let Some(line) = a_curve.gp_lin2d() else {
                    correct = false;
                    break 'edges;
                };
                // `aLine->Position().Direction()` (`:2161-2162`).
                let dir = line.direction();
                correct_u = correct_u || dir.is_parallel(&v_dir, an_angular_tolerance);
                correct_v = correct_v || dir.is_parallel(&an_u_dir, an_angular_tolerance);
                let (f, l) = GeometryRegistry::global().edge_parameters(&an_edge.0);
                bnd_box.add_point(&a_curve.d0(f));
                bnd_box.add_point(&a_curve.d0(l));
            }
        }

        if correct {
            if let Some((umin, umax, vmin, vmax)) = bnd_box.get() {
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
    tol: f64,
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
        correct_surface_boundaries(face, tol * 2.0, &mut umin, &mut umax, &mut vmin, &mut vmax);
    }
    (umin, umax, vmin, vmax)
}

/// `ParameterOutOfBoundary` (`IntTools_FaceFace.cxx:2213`).
///
/// Walks `parameter` off the ON state of either face, returning a new
/// parameter that still lies before `other_parameter`.
pub fn parameter_out_of_boundary(
    parameter: f64,
    curve: &dyn Curve,
    face1: &Face,
    face2: &Face,
    other_parameter: f64,
    b_increase_par: bool,
    tol: f64,
) -> Option<f64> {
    let mut new_parameter = parameter;
    let mut acurpar = parameter;
    let mut state = FaceState::On;
    let mut iter = 0i32;
    let asumtol = tol;
    let mut adelta = asumtol * 0.1;
    adelta = if adelta < CONFUSION { CONFUSION } else { adelta };
    let Some(surf1) = BRepTool::face_surface(face1) else {
        return None;
    };
    let Some(surf2) = BRepTool::face_surface(face2) else {
        return None;
    };
    let Ok(class1) = FClass2d::new(face1, BRepTool::face_tolerance(face1)) else {
        return None;
    };
    let Ok(class2) = FClass2d::new(face2, BRepTool::face_tolerance(face2)) else {
        return None;
    };

    while state == FaceState::On {
        if b_increase_par {
            acurpar += adelta;
        } else {
            acurpar -= adelta;
        }
        let p_current = curve.d0(acurpar);
        if let Some(proj) = project_point_on_surface(surf1.as_ref(), &p_current, 0.0) {
            state = class1.perform(GpPnt2d::new(proj.u, proj.v));
        }
        if state != FaceState::On {
            if let Some(proj) = project_point_on_surface(surf2.as_ref(), &p_current, 0.0) {
                state = class2.perform(GpPnt2d::new(proj.u, proj.v));
            }
        }
        if iter > 11 {
            break;
        }
        iter += 1;
    }

    if iter <= 11 {
        new_parameter = acurpar;
        if b_increase_par {
            if acurpar >= other_parameter {
                new_parameter = other_parameter;
            }
        } else if acurpar <= other_parameter {
            new_parameter = other_parameter;
        }
        Some(new_parameter)
    } else {
        None
    }
}
