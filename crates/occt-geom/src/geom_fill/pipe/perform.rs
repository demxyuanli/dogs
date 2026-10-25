//! `GeomFill_Pipe::Perform` / `KPartT4` / `ApproxSurf`
//! (`GeomFill_Pipe.cxx:747-1102`, OCCT 8.0.0).
//!
//! `Perform` first detects a cylinder / torus (`KPartT4`); otherwise a type 4
//! pipe (path + two guide lines + radius) goes through
//! `Approx_SweepApproximation` with `GeomFill_CircularBlendFunc`, and every
//! other pipe goes through `GeomFill_Sweep`.
//!
//! `ApproxSurf` (`cxx:1012-1101`) builds a `GeomFill_SweepSectionGenerator` on
//! the pipe adaptors (`GeomFill_SweepSectionGenerator.cxx:228-243`, `myType 0`)
//! and approximates the section sequence with `GeomFill_AppSweep`
//! ([`crate::geom_fill::AppSweep`], `AppBlend_AppSurf.gxx:597-1048`).

use std::sync::Arc;

use occt_core::elib::clib;
use occt_core::gp::{GpAx1, GpAx3, GpCylinder, GpDir, GpPnt, GpTorus, GpVec};
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::{ANGULAR, APPROXIMATION, CONFUSION, PAPPROXIMATION};

use crate::approx_sweep_approximation::ApproxSweepApproximation;
use crate::bspline_surface::GeomBSplineSurface;
use crate::cylinder::GeomCylinder;
use crate::geom_fill::app_surf::AppSweep;
use crate::geom_fill::circular_blend_func::CircularBlendFunc;
use crate::geom_fill::line::Line;
use crate::geom_fill::pipe::Pipe;
use crate::geom_fill::pipe_error::GeomFillPipeError;
use crate::geom_fill::sweep::{ApproxStyle, Sweep};
use crate::geom_fill::sweep::expand_knots;
use crate::geom_fill::sweep_section_generator::SweepSectionGenerator;
use crate::rectangular_trimmed::GeomRectangularTrimmedSurface;
use crate::surface::Surface;
use crate::torus::GeomTorus;

impl Pipe {
    /// `GeomFill_Pipe::Perform(const bool WithParameters, const bool Polynomial)`
    /// (`cxx:747-768`).
    pub fn perform(&mut self, with_parameters: bool, polynomial: bool) -> Result<(), String> {
        if self.loc().is_some() && self.sec().is_some() {
            return self.perform_with_tol(1.0e-4, polynomial, Shape::C1, 11, 30);
        }

        self.polynomial = polynomial;
        // On traite le cas tuyau sur arete Type = 4.
        if self.polynomial {
            return self.approx_surf(with_parameters);
        }
        if !self.k_part_t4()? {
            return self.approx_surf(with_parameters);
        }
        Ok(())
    }

    /// `GeomFill_Pipe::Perform(Tol, Polynomial, Conti, DegMax, NbMaxSegment)`
    /// (`cxx:773-864`).
    pub fn perform_with_tol(
        &mut self,
        tol: f64,
        polynomial: bool,
        conti: Shape,
        deg_max: i32,
        nb_max_segment: i32,
    ) -> Result<(), String> {
        if self.status() == GeomFillPipeError::ImpossibleContact {
            return Ok(());
        }

        // `switch (Conti)` (`cxx:785-803`).
        let the_conti = match conti {
            Shape::C0 => Shape::C0,
            Shape::G1 | Shape::C1 => Shape::C1,
            Shape::G2 | Shape::C2 => Shape::C2,
            // `default: TheConti = GeomAbs_C2; // On ne sait pas faire mieux !`
            _ => Shape::C2,
        };

        if self.pipe_type() == 4 {
            if !self.k_part_t4()? {
                let path = self
                    .adp_path()
                    .cloned()
                    .ok_or_else(|| "GeomFill_Pipe: null path".to_string())?;
                let first = self
                    .adp_first_sect()
                    .cloned()
                    .ok_or_else(|| "GeomFill_Pipe: null first section".to_string())?;
                let last = self
                    .adp_last_sect()
                    .cloned()
                    .ok_or_else(|| "GeomFill_Pipe: null last section".to_string())?;

                let func = Arc::new(CircularBlendFunc::new(
                    path.clone(),
                    first,
                    last,
                    self.radius(),
                    polynomial,
                ));
                let mut app = ApproxSweepApproximation::new(func);
                app.perform(
                    path.first_parameter(),
                    path.last_parameter(),
                    tol,
                    tol,
                    0.0,
                    0.01,
                    the_conti,
                    deg_max,
                    nb_max_segment,
                );
                if app.is_done() {
                    if let Some(surface) = bspline_surface_from_approx(&app) {
                        self.surface = Some(Arc::new(surface));
                        self.error = app.max_error_on_surf().unwrap_or(0.0);
                        self.status = GeomFillPipeError::PipeOk;
                    }
                }
            }
        } else if let (Some(loc), Some(sec)) = (self.loc().cloned(), self.sec().cloned()) {
            let mut sweep = Sweep::with_k_part(loc, self.k_part());
            // `Sweep.SetTolerance(Tol)` with the `.hxx:74-77` defaults.
            sweep.set_tolerance(tol, 1.0, 1.0e-5, 1.0);
            sweep.build(
                sec,
                ApproxStyle::Location,
                the_conti,
                deg_max,
                nb_max_segment,
            );
            if sweep.is_done() {
                self.surface = sweep.surface();
                self.error = sweep.error_on_surface();
                self.status = GeomFillPipeError::PipeOk;
            }
        } else {
            return self.perform(true, polynomial);
        }
        Ok(())
    }

    /// `GeomFill_Pipe::KPartT4()` (`cxx:871-1010`): the cylinder and torus cases
    /// of a type 4 pipe. `Ok(false)` means "no particular case", exactly OCCT's
    /// `Ok = false` return.
    pub fn k_part_t4(&mut self) -> Result<bool, String> {
        let (Some(path), Some(first), Some(last)) = (
            self.adp_path().cloned(),
            self.adp_first_sect().cloned(),
            self.adp_last_sect().cloned(),
        ) else {
            return Ok(false);
        };

        // ------- Cas du Cylindre (`cxx:874-935`).
        if let (Some(line0), Some(line1), Some(line2)) =
            (path.gp_line(), first.gp_line(), last.gp_line())
        {
            let a0 = *line0.position();
            let a1 = *line1.position();
            let a2 = *line2.position();
            let d0 = *a0.direction();
            let d1 = *a1.direction();
            let d2 = *a2.direction();
            // Direction must be the same.
            if d0.angle(&d1) > ANGULAR || d1.angle(&d2) > ANGULAR {
                return Ok(false);
            }

            // The length of the line must be the same.
            let l0 = path.last_parameter() - path.first_parameter();
            let l1 = first.last_parameter() - first.first_parameter();
            let l2 = last.last_parameter() - last.first_parameter();
            if (l1 - l0).abs() > CONFUSION || (l2 - l0).abs() > CONFUSION {
                return Ok(false);
            }

            // The first points must be normal to the path.
            let p0 = path.value(path.first_parameter());
            let p1 = first.value(first.first_parameter());
            let p2 = last.value(last.first_parameter());
            let v1 = GpDir::from_vec(&GpVec::from_pnts(&p0, &p1))
                .map_err(|_| "gp_Dir: null magnitude".to_string())?;
            let v2 = GpDir::from_vec(&GpVec::from_pnts(&p0, &p2))
                .map_err(|_| "gp_Dir: null magnitude".to_string())?;
            if v1.dot(&d0).abs() > CONFUSION || v2.dot(&d0).abs() > CONFUSION {
                return Ok(false);
            }

            // The result is a cylindrical surface.
            let z_ref = v1
                .crossed(&v2)
                .map_err(|_| "gp_Dir: zero cross product".to_string())?;
            let mut axis = GpAx3::new(*a0.location(), d0, &v1)
                .map_err(|e| format!("gp_Ax3: {e}"))?;
            if z_ref.dot(&d0) < 0.0 {
                axis.y_reverse();
            }
            // Rotate the surface to set the iso U = 0 not in the result.
            axis.rotate(&GpAx1::new(p0, z_ref), -std::f64::consts::FRAC_PI_2);

            let cylinder = GpCylinder::new(axis, self.radius())
                .map_err(|e| format!("gp_Cylinder: {e}"))?;
            let basis: Arc<dyn Surface> = Arc::new(GeomCylinder::new(cylinder));
            let alpha = GpVec::from_xyz(v1.xyz()).angle_with_ref(
                &GpVec::from_xyz(v2.xyz()),
                &GpVec::from_xyz(z_ref.xyz()),
            );
            let trimmed = GeomRectangularTrimmedSurface::uv(
                basis,
                std::f64::consts::FRAC_PI_2,
                std::f64::consts::FRAC_PI_2 + alpha,
                path.first_parameter(),
                path.last_parameter(),
            )?;
            self.surface = Some(Arc::new(trimmed));
            self.status = GeomFillPipeError::PipeOk;
            return Ok(true);
        }

        // ----------- Cas du tore (`cxx:937-1007`).
        if let (Some(circ0), Some(circ1), Some(circ2)) =
            (path.gp_circ(), first.gp_circ(), last.gp_circ())
        {
            // The three circles must have the same opening angle.
            let alp0 = path.first_parameter() - path.last_parameter();
            let alp1 = first.first_parameter() - first.last_parameter();
            let alp2 = last.first_parameter() - last.last_parameter();
            if (alp0 - alp1).abs() > ANGULAR || (alp0 - alp2).abs() > ANGULAR {
                return Ok(false);
            }

            let a0 = circ0.position();
            let a1 = circ1.position();
            let a2 = circ2.position();
            let d0 = a0.direction();
            let d1 = a1.direction();
            let d2 = a2.direction();
            let p0 = path.value(path.first_parameter());
            let p1 = first.value(first.first_parameter());
            let p2 = last.value(last.first_parameter());

            // The three directions must be equal.
            if d0.angle(&d1) > ANGULAR || d1.angle(&d2) > ANGULAR {
                return Ok(false);
            }

            // The three ax1 must be coincident.
            let l = occt_core::gp::GpLin::new(*a0.axis());
            if !l.contains(&a1.location(), CONFUSION) || !l.contains(&a2.location(), CONFUSION) {
                return Ok(false);
            }

            // The three first points must be in the same section.
            let v1 = GpDir::from_vec(&GpVec::from_pnts(&p0, &p1))
                .map_err(|_| "gp_Dir: null magnitude".to_string())?;
            let v2 = GpDir::from_vec(&GpVec::from_pnts(&p0, &p2))
                .map_err(|_| "gp_Dir: null magnitude".to_string())?;
            let ci = circ0.clone();
            let y_ref = clib::circle_dn(
                &occt_core::gp::GpCirc::new(a0, ci.radius),
                path.first_parameter(),
                1,
            );
            let v1v = GpVec::from_xyz(v1.xyz());
            let v2v = GpVec::from_xyz(v2.xyz());
            if v1v.dot(&y_ref).abs() > CONFUSION || v2v.dot(&y_ref).abs() > CONFUSION {
                return Ok(false);
            }

            // OK it is a toroidal surface.
            let ax3 = GpAx3::new(a0.location(), d0, a0.x_direction())
                .map_err(|e| format!("gp_Ax3: {e}"))?;
            let mut torus =
                GpTorus::new(ax3, ci.radius, self.radius()).map_err(|e| format!("gp_Torus: {e}"))?;
            let x_ref = GpVec::from_pnts(&a0.location(), &p0);
            let y_ref_v = y_ref;
            let mut vv1 = v1v.angle_with_ref(&x_ref, &y_ref_v);
            let mut vv2 = v2v.angle_with_ref(&x_ref, &y_ref_v);
            let delta_v = v2v.angle_with_ref(&v1v, &y_ref_v);
            if delta_v < 0.0 {
                torus.v_reverse();
                vv1 = -vv1;
                vv2 = 2.0 * std::f64::consts::PI + vv1 - delta_v;
            }

            let basis: Arc<dyn Surface> = Arc::new(GeomTorus::new(torus));
            let trimmed = GeomRectangularTrimmedSurface::uv(
                basis,
                path.first_parameter(),
                path.last_parameter(),
                vv1,
                vv2,
            )?;
            self.surface = Some(Arc::new(trimmed));
            self.exch_uv = true;
            self.status = GeomFillPipeError::PipeOk;
            return Ok(true);
        }

        Ok(false)
    }

    /// `GeomFill_Pipe::ApproxSurf(const bool WithParameters)` (`cxx:1012-1101`).
    ///
    /// OCCT throws `Standard_ConstructionError("GeomFill_Pipe")` when
    /// `myType != 4` (`cxx:1022-1025`); that throw is preserved.
    pub fn approx_surf(&mut self, with_parameters: bool) -> Result<(), String> {
        if self.pipe_type() != 4 {
            return Err("GeomFill_Pipe".into());
        }
        let path = self
            .adp_path()
            .cloned()
            .ok_or_else(|| "GeomFill_Pipe".to_string())?;
        let first = self
            .adp_first_sect()
            .cloned()
            .ok_or_else(|| "GeomFill_Pipe".to_string())?;
        let last = self
            .adp_last_sect()
            .cloned()
            .ok_or_else(|| "GeomFill_Pipe".to_string())?;
        let mut section = SweepSectionGenerator::with_adaptors(path, first, last, self.radius())?;
        section.perform(self.polynomial)?;

        let line = Line::with_nb_points(section.nb_sections());
        let mut app = AppSweep::with_settings(
            4,
            8,
            APPROXIMATION,
            PAPPROXIMATION,
            0,
            with_parameters,
        );
        app.perform_nb_max_p(&line, &section, 30)?;

        if app.is_done() {
            if let Some(surface) = bspline_surface_from_app_surf(&app) {
                self.surface = Some(Arc::new(surface));
                let (tol3d, _tol2d) = app.tol_reached();
                self.error = tol3d;
                self.status = GeomFillPipeError::PipeOk;
            }
        }
        Ok(())
    }
}

/// The `new Geom_BSplineSurface(App.SurfPoles(), ...)` of `cxx:1090-1099`.
///
/// `AppBlend_AppSurf` already stores the poles as `tabPoles(u, v)` with the
/// weights split out (`gxx:563-591`), so no transposition is needed here.
fn bspline_surface_from_app_surf(app: &AppSweep) -> Option<GeomBSplineSurface> {
    let (u_degree, v_degree, _nb_u_poles, _nb_v_poles, _nb_u_knots, _nb_v_knots) = app.surf_shape();
    let knots_u = expand_knots(app.surf_u_knots(), app.surf_u_mults());
    let knots_v = expand_knots(app.surf_v_knots(), app.surf_v_mults());
    GeomBSplineSurface::rational(
        app.surf_poles().to_vec(),
        app.surf_weights().to_vec(),
        knots_u,
        knots_v,
        u_degree as usize,
        v_degree as usize,
    )
    .ok()
}

/// The `new Geom_BSplineSurface(App.SurfPoles(), ...)` of `cxx:830-838` /
/// `cxx:1091-1098`.
///
/// `Approx_SweepApproximation` stores `tabPoles(1..Num3DSS, 1..NbPoles)`
/// (section first, pole second), while `Geom_BSplineSurface` indexes
/// `(u = pole, v = section)`; the loops below transpose as
/// `GeomFill_Sweep::BuildAll` does.
fn bspline_surface_from_approx(app: &ApproxSweepApproximation) -> Option<GeomBSplineSurface> {
    let (u_degree, v_degree, nb_u_poles, nb_v_poles, _nb_u_knots, _nb_v_knots) = app.surf_shape()?;
    let poles_tab = app.surf_poles();
    let weights_tab = app.surf_weights();
    let mut poles = vec![vec![GpPnt::zero(); nb_v_poles]; nb_u_poles];
    let mut weights = vec![vec![1.0f64; nb_v_poles]; nb_u_poles];
    for v in 0..nb_v_poles {
        for u in 0..nb_u_poles {
            poles[u][v] = poles_tab[v][u];
            weights[u][v] = weights_tab[v][u];
        }
    }
    let knots_u = expand_knots(app.surf_u_knots(), app.surf_u_mults());
    let knots_v = expand_knots(app.surf_v_knots(), app.surf_v_mults());
    GeomBSplineSurface::rational(
        poles,
        weights,
        knots_u,
        knots_v,
        u_degree as usize,
        v_degree as usize,
    )
    .ok()
}
