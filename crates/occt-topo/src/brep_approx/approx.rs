//! `BRepApprox_Approx` — the `BRepAdaptor_Surface` / `BRepApprox_ApproxLine`
//! instantiation of `ApproxInt_Approx`.
//!
//! Sources:
//! - `ModelingAlgorithms/TKTopAlgo/BRepApprox/BRepApprox_Approx.hxx:45-175`
//!   (the members, the private `Approx_Data` and the surface-typed `Perform`
//!   declarations)
//! - `BRepApprox_Approx_0.cxx:16-98` (the `#define ThePSurface …` binding that
//!   includes `ApproxInt_Approx.gxx`)
//! - `ModelingAlgorithms/TKGeomAlgo/ApproxInt/ApproxInt_Approx.gxx:183-218`
//!   (`Perform(theline, …)`, delegated) and `:225-345` (`Perform(Surf1, Surf2,
//!   …)`, the surface-type dispatch translated here)
//!
//! The `.gxx` body itself is [`ApproxIntApprox`]; this file only adds what the
//! `ThePSurface = BRepAdaptor_Surface` instantiation decides: which evaluator
//! the two-surface arm builds from the surface types (`gxx:241-307`), and the
//! analytic-parametric arm (`gxx:349-396`) over
//! [`ImpPrmSvSurfaces`](crate::brep_approx::ImpPrmSvSurfaces).
//!
//! `isQuadric` (`gxx:241-244`) lists `Plane` / `Cylinder` / `Sphere` / `Cone`
//! only — a `Torus` pair therefore takes the parametric-parametric arm, exactly
//! as in OCCT.

use std::cell::RefCell;
use std::rc::Rc;

use occt_core::approx::ApproxParametrizationType;
use occt_geom::app_par_curves::multi_bspline_curve::MultiBSpCurve;

use crate::approx_int::{ApproxIntApprox, SvSurfaces};
use crate::brep_approx::{ImpPrmSvSurfaces, PrmPrmSvSurfaces};
use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::int_surf_quadric::IntSurfQuadric;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::ApproxLine;
use super::SurfaceTool;

/// `GeomAbs_SurfaceType` membership of the `isQuadric` test
/// (`ApproxInt_Approx.gxx:241-244`): `Plane`, `Cylinder`, `Sphere`, `Cone`.
fn is_quadric_type(t: GeomAbsSurfaceType) -> bool {
    matches!(
        t,
        GeomAbsSurfaceType::Plane
            | GeomAbsSurfaceType::Cylinder
            | GeomAbsSurfaceType::Sphere
            | GeomAbsSurfaceType::Cone
    )
}

/// `BRepApprox_Approx` (`hxx:45`).
#[derive(Debug, Clone)]
pub struct BRepApproxApprox {
    /// Every member of the OCCT class except the private `Perform(Quad, …)`
    /// dispatch state lives in the generic driver.
    driver: ApproxIntApprox,
}

impl Default for BRepApproxApprox {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepApproxApprox {
    /// `BRepApprox_Approx()` (`ApproxInt_Approx.gxx:164-181`).
    pub fn new() -> Self {
        Self {
            driver: ApproxIntApprox::new(),
        }
    }

    /// `Perform(aLine, ApproxXYZ, ApproxU1V1, ApproxU2V2, indicemin,
    /// indicemax)` (`gxx:183-218`): the arm without a surface pair.
    pub fn perform_line(
        &mut self,
        a_line: &ApproxLine,
        approx_xyz: bool,
        approx_u1v1: bool,
        approx_u2v2: bool,
        indicemin: i32,
        indicemax: i32,
    ) {
        self.driver.perform_line(
            a_line, approx_xyz, approx_u1v1, approx_u2v2, indicemin, indicemax,
        );
    }

    /// `Perform(Surf1, Surf2, aLine, ApproxXYZ, ApproxU1V1, ApproxU2V2,
    /// indicemin, indicemax)` (`gxx:225-345`).
    ///
    /// The quadric arm (`gxx:246-307`) builds `IntSurf_Quadric` from the first
    /// surface type that is analytic (`Plane` / `Cylinder` / `Sphere` / `Cone`,
    /// in that `switch` order) and hands the other surface to
    /// `Perform(Quad, …, isTheQuadFirst)` (`gxx:349-396`); the non-analytical
    /// arm builds `ApproxInt_ThePrmPrmSvSurfaces(Surf1, Surf2)` (`gxx:316`).
    pub fn perform(
        &mut self,
        surf1: &SurfaceAdaptor,
        surf2: &SurfaceAdaptor,
        a_line: &ApproxLine,
        approx_xyz: bool,
        approx_u1v1: bool,
        approx_u2v2: bool,
        indicemin: i32,
        indicemax: i32,
    ) {
        self.driver.reset_tol_reached();

        let type_s1 = SurfaceTool::get_type(surf1);
        let type_s2 = SurfaceTool::get_type(surf2);
        let is_quadric = is_quadric_type(type_s1) || is_quadric_type(type_s2);

        if is_quadric {
            // `gxx:248-307`: pick the implicit surface.
            let mut quad = IntSurfQuadric::new();
            let mut second_is_implicit = false;
            match type_s1 {
                GeomAbsSurfaceType::Plane => quad.set_value_plane(&SurfaceTool::plane(surf1)),
                GeomAbsSurfaceType::Cylinder => {
                    quad.set_value_cylinder(&SurfaceTool::cylinder(surf1))
                }
                GeomAbsSurfaceType::Sphere => quad.set_value_sphere(&SurfaceTool::sphere(surf1)),
                GeomAbsSurfaceType::Cone => quad.set_value_cone(&SurfaceTool::cone(surf1)),
                _ => {
                    second_is_implicit = true;
                    match type_s2 {
                        GeomAbsSurfaceType::Plane => {
                            quad.set_value_plane(&SurfaceTool::plane(surf2))
                        }
                        GeomAbsSurfaceType::Cylinder => {
                            quad.set_value_cylinder(&SurfaceTool::cylinder(surf2))
                        }
                        GeomAbsSurfaceType::Sphere => {
                            quad.set_value_sphere(&SurfaceTool::sphere(surf2))
                        }
                        GeomAbsSurfaceType::Cone => {
                            quad.set_value_cone(&SurfaceTool::cone(surf2))
                        }
                        _ => {}
                    }
                }
            }

            // `gxx:296-305`: `Perform(Quad, SecondIsImplicit ? Surf1 : Surf2,
            // …, !SecondIsImplicit)`.
            let p_surf = if second_is_implicit { surf1 } else { surf2 };
            let is_the_quad_first = !second_is_implicit;
            let evaluator: Rc<RefCell<dyn SvSurfaces>> = if is_the_quad_first {
                Rc::new(RefCell::new(ImpPrmSvSurfaces::new_implicit_first(
                    quad,
                    p_surf.clone(),
                )))
            } else {
                Rc::new(RefCell::new(ImpPrmSvSurfaces::new_param_first(
                    p_surf.clone(),
                    quad,
                )))
            };
            self.driver.perform_imp_prm(
                evaluator, a_line, approx_xyz, approx_u1v1, approx_u2v2, indicemin, indicemax,
            );
            return;
        }

        // Non-analytical case: Param-Param perform (`gxx:312-316`).
        let sv = Rc::new(RefCell::new(PrmPrmSvSurfaces::new(
            surf1.clone(),
            surf2.clone(),
        )));
        self.driver.perform_with_surfaces(
            sv, a_line, approx_xyz, approx_u1v1, approx_u2v2, indicemin, indicemax,
        );
    }

    /// `SetParameters(Tol3d, Tol2d, DegMin, DegMax, NbIterMax, NbPntMax,
    /// ApproxWithTangency, Parametrization)` (`gxx:398-428`).
    #[allow(clippy::too_many_arguments)]
    pub fn set_parameters(
        &mut self,
        tol3d: f64,
        tol2d: f64,
        deg_min: i32,
        deg_max: i32,
        nb_iter_max: i32,
        nb_pnt_max: i32,
        approx_with_tangency: bool,
        parametrization: ApproxParametrizationType,
    ) {
        self.driver.set_parameters(
            tol3d,
            tol2d,
            deg_min,
            deg_max,
            nb_iter_max,
            nb_pnt_max,
            approx_with_tangency,
            parametrization,
        );
    }

    /// `TolReached3d()` (`gxx:458-462`).
    pub fn tol_reached_3d(&self) -> f64 {
        self.driver.tol_reached_3d()
    }

    /// `TolReached2d()` (`gxx:464-468`).
    pub fn tol_reached_2d(&self) -> f64 {
        self.driver.tol_reached_2d()
    }

    /// `IsDone()` (`gxx:472-483`).
    pub fn is_done(&self) -> bool {
        self.driver.is_done()
    }

    /// `NbMultiCurves()` (`gxx:432-435`).
    pub fn nb_multi_curves(&self) -> i32 {
        self.driver.nb_multi_curves()
    }

    /// `Value(Index)` (`gxx:486-498`).
    pub fn value(&self) -> &MultiBSpCurve {
        self.driver.value()
    }

    /// The generic driver, for callers that need the `.gxx` members the wrapper
    /// does not re-export.
    pub fn driver(&self) -> &ApproxIntApprox {
        &self.driver
    }
}
