//! `BRepOffset_Tool::EnLargeFace` (`BRepOffset_Tool.cxx:3437-3634`), the
//! surface + face enlargement the offset algorithms use to extend a supporting
//! face beyond its current boundary.
//!
//! The three file-local helpers it drives live next to this file:
//! [`super::enlarge_geometry::enlarge_geometry`] (`cxx:2923-3222`),
//! [`super::extended_face_builder::make_extended_face`] (`cxx:2580-2919`) and
//! [`super::extended_face_builder::update_pcurves`] (`cxx:3230-3263`).
//!
//! Port note: `BB.UpdateFace(F, S, L, tol)` (`cxx:3628`) re-registers the
//! enlarged surface on the input face. The port keys pcurve representations by
//! the `Geom_Surface` handle (`GeometryRegistry::repr_key`), so a *rebuilt*
//! surface (the trimmed / extruded / B-spline arms of `EnlargeGeometry`)
//! orphans the pcurves the input face carried under the old handle. That is the
//! same consequence OCCT has for those arms, because they too install a new
//! `Geom_Surface` handle; the offset arm is the one case OCCT mutates in place,
//! and it is called out as UNPORTED in `enlarge_geometry.rs`.

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::brep_tools;
use crate::shape::Face;
use crate::tgeometry::{FaceGeom, GeometryRegistry};

use super::enlarge_geometry::{enlarge_geometry, EnlargeBounds};
use super::extended_face_builder::{make_extended_face, update_pcurves};
use super::BRepOffsetTool;

/// `TheInfini` (`BRepOffset_Tool.cxx:121`). Not `Precision::Infinite()`: the
/// file uses `1.e+7` as its "infinite enough" bound for `theExtensionMode == 1`.
const THE_INFINI: f64 = 1.0e+7;

/// `TolApex` of `EnLargeFace` (`cxx:3564`).
const TOL_APEX: f64 = 1.0e-5;

/// The non-geometric arguments of `BRepOffset_Tool::EnLargeFace`
/// (`BRepOffset_Tool.hxx:145-156`). The defaults are OCCT's defaults, so
/// `EnLargeFace(F, NF, ChangeGeom)` is
/// `enlarge_face(&f, &EnLargeFaceOptions::default())` after setting
/// `can_extent_surface = ChangeGeom`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EnLargeFaceOptions {
    /// `ChangeGeom`: allow the surface geometry itself to be extended.
    pub can_extent_surface: bool,
    /// `UpDatePCurve`.
    pub update_pcurve: bool,
    /// `enlargeU`.
    pub enlarge_u: bool,
    /// `enlargeVfirst`.
    pub enlarge_vfirst: bool,
    /// `enlargeVlast`.
    pub enlarge_vlast: bool,
    /// `theExtensionMode`: 1 means "use `TheInfini`", anything else means
    /// "ten face spans".
    pub extension_mode: i32,
    /// `theLenBeforeUfirst`.
    pub len_before_ufirst: f64,
    /// `theLenAfterUlast`.
    pub len_after_ulast: f64,
    /// `theLenBeforeVfirst`.
    pub len_before_vfirst: f64,
    /// `theLenAfterVlast`.
    pub len_after_vlast: f64,
}

impl Default for EnLargeFaceOptions {
    fn default() -> Self {
        Self {
            can_extent_surface: false,
            update_pcurve: false,
            enlarge_u: true,
            enlarge_vfirst: true,
            enlarge_vlast: true,
            extension_mode: 1,
            len_before_ufirst: -1.0,
            len_after_ulast: -1.0,
            len_before_vfirst: -1.0,
            len_after_vlast: -1.0,
        }
    }
}

impl BRepOffsetTool {
    /// `BRepOffset_Tool::EnLargeFace(F, NF, CanExtentSurface, UpdatePCurve,
    /// theEnlargeU, theEnlargeVfirst, theEnlargeVlast, theExtensionMode,
    /// theLenBeforeUfirst, theLenAfterUlast, theLenBeforeVfirst,
    /// theLenAfterVlast)` (`cxx:3437-3634`): build `NF`, a face over the same
    /// surface enlarged to the requested UV bounds.
    ///
    /// Returns `(NF, SurfaceChange)`; `SurfaceChange` is `false` when the
    /// bounds were only clamped to the surface's own limits, in which case no
    /// geometry moved. `None` is the null `BRep_Tool::Surface(F)` handle that
    /// OCCT dereferences at `cxx:3461`.
    pub fn enlarge_face(
        f: &Face,
        opts: &EnLargeFaceOptions,
    ) -> Option<(Face, bool)> {
        // `cxx:3454-3455`: `toploc_Location L; S = BRep_Tool::Surface(F, L)`.
        let mut s = BRepTool::face_surface(f)?;

        // `cxx:3461-3471`: a periodic surface uses the tight box, else
        // `BRepTools::UVBounds`.
        let (uf1, uf2, vf1, vf2) = if s.is_u_periodic() || s.is_v_periodic() {
            Self::compact_uv_bounds(f)
        } else {
            brep_tools::uv_bounds(f)
        };

        // `cxx:3472`: `S->Bounds(US1, US2, VS1, VS2)`.
        let (us1, us2) = s.u_range();
        let (vs1, vs2) = s.v_range();

        // `cxx:3473-3489`.
        let coeff;
        let mut bounds;
        if opts.extension_mode == 1 {
            bounds = EnlargeBounds {
                u1: -THE_INFINI,
                u2: THE_INFINI,
                v1: -THE_INFINI,
                v2: THE_INFINI,
            };
            coeff = 0.25;
        } else {
            let face_du = uf2 - uf1;
            let face_dv = vf2 - vf1;
            bounds = EnlargeBounds {
                u1: uf1 - 10.0 * face_du,
                u2: uf2 + 10.0 * face_du,
                v1: vf1 - 10.0 * face_dv,
                v2: vf2 + 10.0 * face_dv,
            };
            coeff = 1.0;
        }

        // `cxx:3491-3519`.
        let mut is_vv1_degen = false;
        let mut is_vv2_degen = false;
        let surface_change;
        if opts.can_extent_surface {
            surface_change = enlarge_geometry(
                &mut s,
                &mut bounds,
                &mut is_vv1_degen,
                &mut is_vv2_degen,
                uf1,
                uf2,
                vf1,
                vf2,
                coeff,
                opts.enlarge_u,
                opts.enlarge_vfirst,
                opts.enlarge_vlast,
                opts.len_before_ufirst,
                opts.len_after_ulast,
                opts.len_before_vfirst,
                opts.len_after_vlast,
            );
        } else {
            surface_change = false;
            bounds.u1 = us1.max(bounds.u1);
            bounds.u2 = us2.min(bounds.u2);
            bounds.v1 = vs1.max(bounds.v1);
            bounds.v2 = vs2.min(bounds.v2);
        }

        // `cxx:3521-3547`: keep at most one period.
        let mut uperiodic = false;
        let mut vperiodic = false;
        if s.is_u_periodic() {
            uperiodic = true;
            let period = s.u_period();
            let delta = period - (uf2 - uf1);
            let alpha = 0.1;
            bounds.u1 = uf1 - alpha * delta;
            bounds.u2 = uf2 + alpha * delta;
            if bounds.u2 - bounds.u1 > period {
                bounds.u2 = bounds.u1 + period;
            }
        }
        if s.is_v_periodic() {
            vperiodic = true;
            let period = s.v_period();
            let delta = period - (vf2 - vf1);
            let alpha = 0.1;
            bounds.v1 = vf1 - alpha * delta;
            bounds.v2 = vf2 + alpha * delta;
            if bounds.v2 - bounds.v1 > period {
                bounds.v2 = bounds.v1 + period;
            }
        }

        // `cxx:3549-3577`: a conical surface must not be extended past its
        // apex.
        let the_surf = s.rectangular_trimmed_basis().unwrap_or_else(|| s.clone());
        if let Some(cone) = the_surf.gp_cone() {
            let apex = cone.apex();
            let (_, vapex) = occt_core::elib::slib::cone_parameters(
                &cone.position(),
                cone.radius(),
                cone.semi_angle(),
                &apex,
            );
            if bounds.v1 < vapex && vapex < bounds.v2 {
                // `cxx:3565-3575`: both face bounds are taken to be on the same
                // side of the apex.
                if vapex - vf1 >= TOL_APEX || vapex - vf2 >= TOL_APEX {
                    bounds.v2 = vapex;
                } else {
                    bounds.v1 = vapex;
                }
            }
        }

        // `cxx:3578-3588`.
        if !opts.enlarge_u {
            bounds.u1 = uf1;
            bounds.u2 = uf2;
        }
        if !opts.enlarge_vfirst {
            bounds.v1 = vf1;
        }
        if !opts.enlarge_vlast {
            bounds.v2 = vf2;
        }

        // `cxx:3591-3601`: `BRepTools::DetectClosedness`, then drop a physical
        // (non-periodic) closedness when the caller fixed lengths. Note the
        // `||` for U and the `&&` for V, exactly as OCCT has it.
        let (mut uclosed, mut vclosed) = brep_tools::detect_closedness(f);
        if uclosed
            && !uperiodic
            && (opts.len_before_ufirst != 0.0 || opts.len_after_ulast != 0.0)
        {
            uclosed = false;
        }
        if vclosed && !vperiodic && opts.len_before_vfirst != 0.0 && opts.len_after_vlast != 0.0 {
            vclosed = false;
        }

        // `cxx:3603-3604`.
        let mut bf = make_extended_face(
            &s,
            bounds.u1,
            bounds.u2,
            bounds.v1,
            bounds.v2,
            uclosed,
            vclosed,
            is_vv1_degen,
            is_vv2_degen,
        );
        // `BF.Location(L)`.
        let loc = f.0.location().clone();
        bf.0.set_location(&loc);

        // `cxx:3626-3630`.
        if surface_change && opts.update_pcurve {
            let fwd = Face(f.0.oriented(Orientation::Forward));
            update_pcurves(&fwd, &bf);
            // `BRep_Builder::UpdateFace(F, S, L, BRep_Tool::Tolerance(F))`.
            let reg = GeometryRegistry::global();
            let mut geom = FaceGeom::new(s.clone());
            geom.tolerance = reg.face_tolerance(&f.0);
            geom.natural_restriction = reg.natural_restriction(&f.0);
            reg.set_face(&f.0, geom);
        }

        // `cxx:3632`.
        bf.0.set_orientation(f.0.orientation());
        Some((bf, surface_change))
    }
}
