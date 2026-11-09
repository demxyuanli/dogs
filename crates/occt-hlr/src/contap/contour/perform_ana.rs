//! `Contap_Contour::PerformAna` (`Contap_Contour.cxx:2156-2389`): the
//! analytical arm for the four elementary surfaces.

use occt_core::precision::{ANGULAR, CONFUSION};

use occt_topo::brep_top_adaptor_topol_tool::TopolTool;
use occt_topo::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use occt_topo::geom_abs_curve_type::GeomAbsCurveType;

use super::super::cont_ana::ContapContAna;
use super::super::line::ContapLine;
use super::super::types::{ContapIType, ContapTFunction};
use super::line_constructor::line_constructor;
use super::transitions::{compute_transition_on_gp_circle, compute_transition_on_gp_line};
use super::{points_on_line::put_points_on_line, process_segments::process_segments, ContapContour};

impl ContapContour {
    /// `PerformAna(Domain)` (`cxx:2156-2389`).
    pub(crate) fn perform_ana(&mut self, domain: &mut TopolTool) {
        self.done = false;
        self.slin.clear();

        let tol_arc = 1.0e-5;
        let mut perform_sol_rst = true;

        let Some(surf_handle) = self.my_sfunc.surface().cloned() else {
            self.done = true;
            return;
        };
        let type_func = self.my_sfunc.function_type();
        let typ_s = surf_handle.get_type();

        let mut contana = ContapContAna::new();
        let mut theline = ContapLine::new();

        // `cxx:2177-2291`: the per-surface switch. The plane arm only decides
        // whether the restriction search runs; the other arms fill `contana`.
        match typ_s {
            GeomAbsSurfaceType::Plane => {
                let Some(pl) = surf_handle.plane() else {
                    self.done = true;
                    return;
                };
                match type_func {
                    ContapTFunction::ContourStd => {
                        let dirpln = pl.position().direction();
                        if self.my_sfunc.direction().dot(&dirpln).abs() > ANGULAR {
                            // No point of the plane is a solution, so neither
                            // is any point on the restrictions.
                            perform_sol_rst = false;
                        }
                    }
                    ContapTFunction::ContourPrs => {
                        let eye = *self.my_sfunc.eye();
                        if pl.signed_distance_pnt(&eye).abs() > CONFUSION {
                            perform_sol_rst = false;
                        }
                    }
                    ContapTFunction::DraftStd => {
                        let dirpln = pl.position().direction();
                        let sina = self.my_sfunc.angle().sin();
                        if (self.my_sfunc.direction().dot(&dirpln) + sina).abs() > ANGULAR {
                            perform_sol_rst = false;
                        }
                    }
                    ContapTFunction::DraftPrs => {}
                }
            }
            GeomAbsSurfaceType::Sphere => {
                let Some(sph) = surf_handle.sphere() else {
                    self.done = true;
                    return;
                };
                match type_func {
                    ContapTFunction::ContourStd => {
                        contana.perform_sphere_dir(&sph, self.my_sfunc.direction());
                    }
                    ContapTFunction::ContourPrs => {
                        contana.perform_sphere_eye(&sph, self.my_sfunc.eye());
                    }
                    ContapTFunction::DraftStd => {
                        contana.perform_sphere_dir_angle(
                            &sph,
                            self.my_sfunc.direction(),
                            self.my_sfunc.angle(),
                        );
                    }
                    ContapTFunction::DraftPrs => {}
                }
            }
            GeomAbsSurfaceType::Cylinder => {
                let Some(cyl) = surf_handle.cylinder() else {
                    self.done = true;
                    return;
                };
                match type_func {
                    ContapTFunction::ContourStd => {
                        contana.perform_cylinder_dir(&cyl, self.my_sfunc.direction());
                    }
                    ContapTFunction::ContourPrs => {
                        contana.perform_cylinder_eye(&cyl, self.my_sfunc.eye());
                    }
                    ContapTFunction::DraftStd => {
                        contana.perform_cylinder_dir_angle(
                            &cyl,
                            self.my_sfunc.direction(),
                            self.my_sfunc.angle(),
                        );
                    }
                    ContapTFunction::DraftPrs => {}
                }
            }
            GeomAbsSurfaceType::Cone => {
                let Some(cone) = surf_handle.cone() else {
                    self.done = true;
                    return;
                };
                match type_func {
                    ContapTFunction::ContourStd => {
                        contana.perform_cone_dir(&cone, self.my_sfunc.direction());
                    }
                    ContapTFunction::ContourPrs => {
                        contana.perform_cone_eye(&cone, self.my_sfunc.eye());
                    }
                    ContapTFunction::DraftStd => {
                        contana.perform_cone_dir_angle(
                            &cone,
                            self.my_sfunc.direction(),
                            self.my_sfunc.angle(),
                        );
                    }
                    ContapTFunction::DraftPrs => {}
                }
            }
            _ => {}
        }

        // `cxx:2293-2346`.
        if typ_s != GeomAbsSurfaceType::Plane {
            if !contana.is_done() {
                return;
            }
            let nb_cont = match contana.nb_contours() {
                Ok(n) => n,
                Err(_) => return,
            };
            if nb_cont == 0 {
                self.done = true;
                return;
            }
            let Ok(typ_l) = contana.type_contour() else {
                return;
            };
            if typ_l == GeomAbsCurveType::Circle {
                let Ok(circle) = contana.circle() else {
                    return;
                };
                theline.set_value_circle(&circle);
                let trans = compute_transition_on_gp_circle(&mut self.my_sfunc, &circle);
                theline.set_transition_on_s(trans);
                self.slin.push(theline.clone());
            } else if typ_l == GeomAbsCurveType::Line {
                // `cxx:2317-2328`. The cone-apex block (`cxx:2330-2342`) is
                // commented out in OCCT.
                for i in 1..=nb_cont {
                    let Ok(l) = contana.line(i) else {
                        continue;
                    };
                    theline.set_value_line(&l);
                    let trans = compute_transition_on_gp_line(&mut self.my_sfunc, &l);
                    theline.set_transition_on_s(trans);
                    self.slin.push(theline.clone());
                    theline.clear();
                }
            }
        }

        // `cxx:2348-2387`.
        if perform_sol_rst {
            let perform = self.solrst.perform(
                &mut self.my_afunc,
                domain,
                tol_arc,
                tol_arc,
                false,
            );
            if perform.is_err() || !self.solrst.is_done() {
                return;
            }
            let nb_point_rst = self.solrst.nb_points().unwrap_or(0);

            if nb_point_rst != 0 {
                put_points_on_line(&self.solrst, &surf_handle, &mut self.slin);
            }

            if self.solrst.nb_segments().unwrap_or(0) != 0 {
                process_segments(
                    &self.solrst,
                    &mut self.slin,
                    tol_arc,
                    &mut self.my_sfunc,
                    domain,
                );
            }

            let nblinto = self.slin.len();
            let mut seq_to_destroy: Vec<usize> = Vec::new();
            for i in 1..=nblinto {
                if self.slin[i - 1].type_contour() != ContapIType::Restriction {
                    let l = self.slin[i - 1].clone();
                    line_constructor(&mut self.slin, domain, &l, &surf_handle);
                    seq_to_destroy.push(i);
                }
            }
            for i in seq_to_destroy.iter().rev() {
                self.slin.remove(i - 1);
            }
        }

        self.done = true;
    }
}
