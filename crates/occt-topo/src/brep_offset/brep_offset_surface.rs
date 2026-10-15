//! `BRepOffset::Surface` (`BRepOffset.cxx:44-200`): the offset of a
//! `Geom_Surface` at a signed distance.
//!
//! A surface type that has an analytic offset keeps that type (plane -> plane,
//! cylinder -> cylinder, cone -> cone, sphere -> sphere, torus -> torus). Every
//! other type (`Geom_SurfaceOfRevolution`, `Geom_SurfaceOfLinearExtrusion`,
//! `Geom_BSplineSurface`, `Geom_BezierSurface`, `Geom_OffsetSurface`) falls
//! through to an exact `Geom_OffsetSurface` (`BRepOffset.cxx:193-197`).
//! `Geom_RectangularTrimmedSurface` recurses on its basis and re-wraps the same
//! window (`:180-191`). No branch samples the surface onto a grid.
//!
//! Two port-side deviations, both structural rather than numeric:
//!
//! * OCCT dispatches on `Surface->DynamicType()` (the exact class). The port's
//!   `Surface::gp_pln` / `gp_cylinder` / `gp_cone` / `gp_sphere` / `gp_torus`
//!   also answer for a `Geom_RectangularTrimmedSurface`
//!   (`occt-geom/src/rectangular_trimmed.rs:604-621`, the documented
//!   `GeomAdaptor_Surface::GetType()` equivalent). The trim branch is therefore
//!   tested first, which restores OCCT's dispatch order.
//! * `allowC0` is only passed to the `Geom_OffsetSurface` constructor in OCCT
//!   (`:196`). The port's `GeomOffsetSurface::new`
//!   (`occt-geom/src/offset_surface.rs:61`) always builds its osculating probe,
//!   i.e. the `allowC0 = false` object; the flag is therefore UNPORTED and the
//!   caller must pass `false`.

use std::f64::consts::PI;

use occt_core::gp::{GpAx1, GpCone, GpTorus};
use occt_core::precision::CONFUSION;
use occt_geom::{GeomCone, GeomOffsetSurface, GeomRectangularTrimmedSurface, GeomTorus};

use super::prelude::*;

/// `BRepOffset_Status` (`BRepOffset_Status.hxx:24-31`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BRepOffsetStatus {
    Good,
    Reversed,
    Degenerated,
    Unknown,
}

impl BRepOffsetStatus {
    /// The `BRepOffset_Status` ordinal, so probes can compare against OCCT.
    pub fn code(self) -> i32 {
        match self {
            Self::Good => 0,
            Self::Reversed => 1,
            Self::Degenerated => 2,
            Self::Unknown => 3,
        }
    }

    pub fn is_good(self) -> bool {
        matches!(self, Self::Good)
    }
}

/// `BRepOffset::Surface` (`BRepOffset.cxx:44-200`).
///
/// Returns the offset surface of `s` at `offset` plus the `BRepOffset_Status`
/// OCCT reports through its output parameter. `Err` is OCCT's
/// `Standard_ConstructionError` from a `gp_*` constructor.
pub fn brep_offset_surface(
    s: &dyn Surface,
    offset: f64,
    allow_c0: bool,
) -> Result<(Box<dyn Surface>, BRepOffsetStatus), String> {
    // `BRepOffset.cxx:48-49`.
    let _ = allow_c0; // UNPORTED, see the module header.
    let mut status = BRepOffsetStatus::Good;
    let mut result: Option<Box<dyn Surface>> = None;

    if let Some(basis) = s.rectangular_trimmed_basis() {
        // `:181-191`: offset the basis and keep the window. Tested first: the
        // port's elementary queries answer for a trimmed patch as well.
        let (u1, u2, v1, v2) = match s.rectangular_trimmed_spec() {
            Some((u1, u2, v1, v2, _, _)) => (u1, u2, v1, v2),
            None => (s.u_range().0, s.u_range().1, s.v_range().0, s.v_range().1),
        };
        let (off, st) = brep_offset_surface(basis.as_ref(), offset, allow_c0)?;
        status = st;
        let trimmed = GeomRectangularTrimmedSurface::uv(Arc::from(off), u1, u2, v1, v2)?;
        result = Some(Box::new(trimmed));
    } else if let Some(p) = s.gp_pln() {
        // `:54-61`: `T = Position().XDirection() ^ Position().YDirection(); T *= Offset; P->Translated(T)`.
        let ax = p.position();
        let t = ax.x_direction().xyz().crossed(ax.y_direction().xyz()).multiplied(offset);
        result = Some(Box::new(GeomPlane::new(p.translated_vec(&GpVec::from_xyz(&t)))));
    } else if let Some(c) = s.gp_cylinder() {
        // `:62-90`.
        let mut radius = c.radius();
        let mut axis = c.position();
        if axis.is_direct() {
            radius += offset;
        } else {
            radius -= offset;
        }
        if radius >= CONFUSION {
            result = Some(Box::new(GeomCylinder::new(gp_cylinder(axis, radius)?)));
        } else if radius <= -CONFUSION {
            let rot = GpAx1::new(axis.location(), axis.direction());
            axis.rotate(&rot, PI);
            result = Some(Box::new(GeomCylinder::new(gp_cylinder(axis, radius.abs())?)));
            status = BRepOffsetStatus::Reversed;
        } else {
            status = BRepOffsetStatus::Degenerated;
        }
    } else if let Some(c) = s.gp_cone() {
        // `:92-111`: `Radius = RefRadius + Offset * cos(Alpha)`, and the
        // reference plane slides along the axis by `-Offset * sin(Alpha)`.
        let mut alpha = c.semi_angle();
        let mut radius = c.radius() + offset * alpha.cos();
        let mut axis = c.position();
        let z = axis.direction().xyz().multiplied(-offset * alpha.sin());
        axis.translate_vec(&GpVec::from_xyz(&z));
        if radius < 0.0 {
            radius = -radius;
            let rot = GpAx1::new(axis.location(), axis.direction());
            axis.rotate(&rot, PI);
            alpha = -alpha;
        }
        // `gp_Cone` accepts a negative semi-angle (`gp_Cone.hxx:240-248`);
        // `GpCone::new` (`occt-core/src/gp/cone.rs:27-34`) does not, so the
        // OCCT-legal fields are set directly.
        let cone = GpCone { pos: axis, radius, semi_angle: alpha };
        result = Some(Box::new(GeomCone::new(cone)));
    } else if let Some(sp) = s.gp_sphere() {
        // `:113-141`.
        let mut radius = sp.radius;
        let mut axis = sp.pos;
        if axis.is_direct() {
            radius += offset;
        } else {
            radius -= offset;
        }
        if radius >= CONFUSION {
            result = Some(Box::new(GeomSphere::new(GpSphere::new(axis, radius).map_err(|e| e.to_string())?)));
        } else if radius <= -CONFUSION {
            let rot = GpAx1::new(axis.location(), axis.direction());
            axis.rotate(&rot, PI);
            axis.z_reverse();
            result = Some(Box::new(GeomSphere::new(GpSphere::new(axis, -radius).map_err(|e| e.to_string())?)));
            status = BRepOffsetStatus::Reversed;
        } else {
            status = BRepOffsetStatus::Degenerated;
        }
    } else if let Some(t) = s.gp_torus() {
        // `:143-171`: only a spindle torus (`MinorRadius < MajorRadius`) has an
        // analytic offset; any other torus falls through to `Geom_OffsetSurface`.
        let major = t.major_radius;
        let mut minor = t.minor_radius;
        let axis = t.pos;
        if minor < major {
            if axis.is_direct() {
                minor += offset;
            } else {
                minor -= offset;
            }
            if minor >= CONFUSION {
                result = Some(Box::new(GeomTorus::new(GpTorus::new(axis, major, minor).map_err(|e| e.to_string())?)));
            } else if minor <= -CONFUSION {
                status = BRepOffsetStatus::Reversed;
            } else {
                status = BRepOffsetStatus::Degenerated;
            }
        }
    }
    // `Geom_SurfaceOfRevolution` (`:172-175`), `Geom_SurfaceOfLinearExtrusion`
    // (`:176-179`), `Geom_BSplineSurface` (`:179-181`) and
    // `Geom_OffsetSurface` (`:192-194`) have an empty OCCT branch: they fall
    // through to the `Geom_OffsetSurface` below.

    match result {
        Some(r) => Ok((r, status)),
        // `:193-197`.
        None => Ok((Box::new(GeomOffsetSurface::new(Arc::from(s.clone_dyn()), offset)), status)),
    }
}

fn gp_cylinder(pos: GpAx3, radius: f64) -> Result<GpCylinder, String> {
    GpCylinder::new(pos, radius).map_err(|e| e.to_string())
}
