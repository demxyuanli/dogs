//! `Contap_Contour::Perform` (`Contap_Contour.cxx:1541-1971`): the walking arm
//! used for every surface that is not one of the four quads.

use std::sync::Arc;

use occt_core::bnd::BndBox;
use occt_core::gp::{GpDir, GpPnt2d, GpVec};
use occt_core::precision::{CONFUSION, INFINITE};

use occt_topo::brep_adaptor_curve2d::Curve2dAdaptor;
use occt_topo::brep_top_adaptor_topol_tool::TopolTool;
use occt_topo::geom_bnd_lib_surface3d::add_surface_with_bounds;
use occt_topo::intsurf_transition::{make_transition, IntSurfTransition};

use super::super::iwalking::types::InteriorPoint;
use super::super::iwalking::ContapTheIWalking;
use super::super::line::ContapLine;
use super::super::point::ContapPoint;
use super::super::surf_props;
use super::super::types::{ContapIType, ContapTFunction};
use super::super::iwalking::types::PathPoint;
use super::super::tools::{hcont_project, hcurve2d_d1};
use super::compute_tangency::compute_tangency;
use super::internal_points::compute_internal_points;
use super::keep_inside_points::keep_inside_points;
use super::line_constructor::line_constructor;
use super::process_segments::process_segments;
use super::transitions::compute_transition_on_line;
use super::ContapContour;

impl ContapContour {
    /// `Perform(Domain)` for the walking arm (`cxx:1541-1971`).
    pub(crate) fn perform_domain(&mut self, domain: &mut TopolTool) {
        self.done = false;
        self.slin.clear();

        let mut tol_arc = CONFUSION;
        let Some(surf) = self.my_sfunc.surface().cloned() else {
            self.done = true;
            return;
        };

        let eps_u = surf.u_resolution(CONFUSION);
        let eps_v = surf.v_resolution(CONFUSION);
        let preci = eps_u.min(eps_v);
        let mut fleche = 0.01;
        let pas = 0.005;

        // `cxx:1579-1625`: scale `Fleche`/`TolArc` by the largest extent of the
        // surface box; infinite bounds skip the box entirely.
        let (uinf, vinf) = (surf.first_u_parameter(), surf.first_v_parameter());
        let (usup, vsup) = (surf.last_u_parameter(), surf.last_v_parameter());
        let infinity = uinf == -INFINITE
            || usup == INFINITE
            || vinf == -INFINITE
            || vsup == INFINITE;

        let mut dx;
        if infinity {
            dx = 1.0;
        } else {
            let mut b1 = BndBox::new();
            if let Some(basis) = surf.surface() {
                add_surface_with_bounds(basis, 1.0e-8, &mut b1);
            }
            if let Some((x0, y0, z0, x1, y1, z1)) = b1.get() {
                dx = x1 - x0;
                let dy = y1 - y0;
                let dz = z1 - z0;
                if dx < dy {
                    dx = dy;
                }
                if dx < dz {
                    dx = dz;
                }
            } else {
                dx = 1.0;
            }
        }
        if dx > 10000.0 {
            dx = 10000.0;
        }
        fleche *= dx;
        tol_arc *= dx;

        // `cxx:1633`: the function tolerance.
        self.my_sfunc.set_tolerance(CONFUSION);

        let _ = pas;
        let recheck_on_regularity = true;
        let performed = self.solrst.perform(
            &mut self.my_afunc,
            domain,
            tol_arc,
            tol_arc,
            recheck_on_regularity,
        );
        if performed.is_err() || !self.solrst.is_done() {
            return;
        }

        let nb_point_rst = self.solrst.nb_points().unwrap_or(0);
        // `Destination(1, NbPointRst + 1)`, `Init(0)`; index 0 is unused.
        let mut destination: Vec<i32> = vec![0; (nb_point_rst + 1) as usize];
        let mut seqpdep: Vec<PathPoint> = Vec::new();
        if nb_point_rst != 0 {
            compute_tangency(
                &self.solrst,
                domain,
                &mut self.my_sfunc,
                &mut seqpdep,
                &mut destination,
            );
        }

        // `cxx:1652`: the interior search.
        self.solins.perform(&mut self.my_sfunc, &surf, domain, CONFUSION);
        let nb_point_ins = self.solins.nb_points();
        let mut seqpins: Vec<InteriorPoint> = Vec::new();

        if nb_point_ins != 0 {
            let mut keep_all_points = false;
            // `cxx:1657-1676`: a torus whose axis is perpendicular to the
            // projection direction keeps every interior point.
            if self.solrst.nb_segments().unwrap_or(0) <= 0
                && self.my_sfunc.function_type() == ContapTFunction::ContourStd
            {
                if let Some(torus) = surf.torus() {
                    let tor_dir = torus.position().direction();
                    let proj_dir = *self.my_sfunc.direction();
                    if tor_dir.dot(&proj_dir) < CONFUSION {
                        keep_all_points = true;
                    }
                }
            }

            if keep_all_points {
                for p in self.solins.points() {
                    seqpins.push(p.clone());
                }
            } else {
                seqpins = keep_inside_points(&self.solins, &self.solrst, &mut self.my_sfunc);
            }
        }

        if !seqpdep.is_empty() || !seqpins.is_empty() {
            let to_fill_holes = true;
            let mut iwalk = ContapTheIWalking::new(preci, fleche, pas, to_fill_holes);
            iwalk.perform(&seqpdep, &seqpins, &mut self.my_sfunc, &surf, false);
            if !iwalk.is_done() {
                return;
            }

            let nblines = iwalk.nb_lines().unwrap_or(0);
            // `cxx:1559` reuses `ptdeb`/`ptfin` across lines.
            let mut ptdeb = ContapPoint::new();
            let mut ptfin = ContapPoint::new();

            for j in 1..=nblines {
                let Some(iwline) = iwalk.value(j) else {
                    continue;
                };
                let nbpts = iwline.nb_points();

                let mut wl = occt_topo::int_tools_wline::WLine::new();
                for i in 1..=nbpts.min(iwline.line().len() as i32) {
                    wl.add(iwline.value(i).clone());
                }
                let mut theline = ContapLine::new();
                theline.set_line_on_2s(wl);

                // `cxx:1715-1718`: OCCT reads `k` before it is initialized to
                // pick the tangent point (`Contap_Contour.cxx:1715`); the port
                // uses the first point, which is also the value passed to
                // `TangentVector` (`vcttg` does not depend on the index).
                let tv = iwline.tangent_vector().0.clone();
                let (u, v) = iwline.value(1).parameters_on_s2();
                let type_trans = compute_transition_on_line(&mut self.my_sfunc, u, v, &tv);
                theline.set_transition_on_s(type_trans);

                // `cxx:1728-1713`: add the first walking point unless present.
                if iwline.has_first_point() {
                    let indfirst = iwline.first_point_index();
                    let Some(ppoint) = seqpdep.get((indfirst - 1).max(0) as usize) else {
                        continue;
                    };
                    let mut themult = ppoint.multiplicity();
                    let mut i = nb_point_rst;
                    while i >= 1 {
                        if destination[i as usize] == indfirst {
                            let Some((u, v)) = ppoint.parameters(themult) else {
                                i -= 1;
                                continue;
                            };
                            ptdeb.set_value(*ppoint.value(), u, v);
                            ptdeb.set_parameter(1.0);

                            let Some(pstart) = self.solrst.point(i as usize) else {
                                i -= 1;
                                continue;
                            };
                            let currentarc = pstart.arc().clone();
                            let currentparam = pstart.parameter();
                            let mut t_line = IntSurfTransition::new();
                            let mut t_arc = IntSurfTransition::new();
                            if !iwline.is_tangent_at_beginning() {
                                let (pt2d, d2d) = hcurve2d_d1(&currentarc, currentparam);
                                if let Some((_p, d1u, d1v, normale)) = surf_props::deriv_and_norm(
                                    &surf,
                                    pt2d.x(),
                                    pt2d.y(),
                                ) {
                                    let tgtrst = GpVec::new(
                                        d2d.x() * d1u.x() + d2d.y() * d1v.x(),
                                        d2d.x() * d1u.y() + d2d.y() * d1v.y(),
                                        d2d.x() * d1u.z() + d2d.y() * d1v.z(),
                                    );
                                    let ndir = GpDir::from_vec(&normale).unwrap_or_else(|_| GpDir::new(1.0, 0.0, 0.0).expect("+X"));
                                    let dir = ppoint
                                        .direction_3d()
                                        .cloned()
                                        .unwrap_or_else(GpVec::default);
                                    make_transition(&dir, &tgtrst, &ndir, &mut t_line, &mut t_arc);
                                }
                            } else {
                                t_line.set_undecided();
                                t_arc.set_undecided();
                            }

                            ptdeb.set_arc(currentarc, currentparam, t_line, t_arc);
                            if !pstart.is_new() {
                                if let Some(vt) = pstart.vertex() {
                                    ptdeb.set_vertex(vt.clone());
                                }
                            }
                            theline.add_vertex(ptdeb.clone());
                            themult -= 1;
                        }
                        i -= 1;
                    }
                } else {
                    let (u, v) = iwline.value(1).parameters_on_s2();
                    let p0 = theline.point(1).map(|p| p.value()).unwrap_or_default();
                    ptdeb.set_value(p0, u, v);
                    ptdeb.set_parameter(1.0);
                    theline.add_vertex(ptdeb.clone());
                }

                if iwline.has_last_point() {
                    let indlast = iwline.last_point_index();
                    let Some(ppoint) = seqpdep.get((indlast - 1).max(0) as usize) else {
                        continue;
                    };
                    let mut themult = ppoint.multiplicity();
                    let mut i = nb_point_rst;
                    while i >= 1 {
                        if destination[i as usize] == indlast {
                            let Some((u, v)) = ppoint.parameters(themult) else {
                                i -= 1;
                                continue;
                            };
                            ptfin.set_value(*ppoint.value(), u, v);
                            ptfin.set_parameter(f64::from(nbpts));

                            let Some(pstart) = self.solrst.point(i as usize) else {
                                i -= 1;
                                continue;
                            };
                            let currentarc = pstart.arc().clone();
                            let currentparam = pstart.parameter();
                            let mut t_line = IntSurfTransition::new();
                            let mut t_arc = IntSurfTransition::new();
                            if !iwline.is_tangent_at_end() {
                                let (pt2d, d2d) = hcurve2d_d1(&currentarc, currentparam);
                                if let Some((_p, d1u, d1v, normale)) = surf_props::deriv_and_norm(
                                    &surf,
                                    pt2d.x(),
                                    pt2d.y(),
                                ) {
                                    let tgtrst = GpVec::new(
                                        d2d.x() * d1u.x() + d2d.y() * d1v.x(),
                                        d2d.x() * d1u.y() + d2d.y() * d1v.y(),
                                        d2d.x() * d1u.z() + d2d.y() * d1v.z(),
                                    );
                                    let ndir = GpDir::from_vec(&normale).unwrap_or_else(|_| GpDir::new(1.0, 0.0, 0.0).expect("+X"));
                                    let dir = ppoint
                                        .direction_3d()
                                        .map(|d| GpVec::new(-d.x(), -d.y(), -d.z()))
                                        .unwrap_or_else(GpVec::default);
                                    make_transition(&dir, &tgtrst, &ndir, &mut t_line, &mut t_arc);
                                }
                            } else {
                                t_line.set_undecided();
                                t_arc.set_undecided();
                            }

                            ptfin.set_arc(currentarc, currentparam, t_line, t_arc);
                            if !pstart.is_new() {
                                if let Some(vt) = pstart.vertex() {
                                    ptfin.set_vertex(vt.clone());
                                }
                            }
                            theline.add_vertex(ptfin.clone());
                            themult -= 1;
                        }
                        i -= 1;
                    }
                } else {
                    let (u, v) = iwline.value(nbpts).parameters_on_s2();
                    let pn = theline
                        .point(nbpts as usize)
                        .map(|p| p.value())
                        .unwrap_or_default();
                    ptfin.set_value(pn, u, v);
                    ptfin.set_parameter(f64::from(nbpts));
                    theline.add_vertex(ptfin.clone());
                }

                compute_internal_points(&mut theline, &mut self.my_sfunc, eps_u, eps_v);
                line_constructor(&mut self.slin, domain, &theline, &surf);
                theline.reset_seq_of_vertex();
            }

            // `cxx:1839-1871`: coincident end points become multiple.
            let nblines = self.slin.len();
            let mut pairs: Vec<(usize, usize, usize, usize)> = Vec::new();
            for j in 0..nblines.saturating_sub(1) {
                let nbvt1 = self.slin[j].nb_vertex();
                for ivt1 in 1..=nbvt1 {
                    let pttg1 = match self.slin[j].vertex(ivt1) {
                        Some(v) if !v.is_on_arc() => *v.value(),
                        _ => continue,
                    };
                    for k in (j + 1)..nblines {
                        let nbvt2 = self.slin[k].nb_vertex();
                        for ivt2 in 1..=nbvt2 {
                            let pttg2 = match self.slin[k].vertex(ivt2) {
                                Some(v) if !v.is_on_arc() => *v.value(),
                                _ => continue,
                            };
                            if pttg1.distance(&pttg2) <= tol_arc {
                                pairs.push((j, ivt1, k, ivt2));
                            }
                        }
                    }
                }
            }
            for (j, ivt1, k, ivt2) in pairs {
                if let Some(v) = self.slin[j].vertex_mut(ivt1) {
                    v.set_multiple();
                }
                if let Some(v) = self.slin[k].vertex_mut(ivt2) {
                    v.set_multiple();
                }
            }
        }

        // `cxx:1875-1878`.
        if self.solrst.nb_segments().unwrap_or(0) != 0 {
            process_segments(&self.solrst, &mut self.slin, tol_arc, &mut self.my_sfunc, domain);
        }

        // `cxx:1881-1968`: attach walking vertices that lie on a solution
        // restriction and add the corresponding restriction point.
        if self.solrst.nb_segments().unwrap_or(0) != 0 {
            let nblines = self.slin.len();
            for j in 0..nblines {
                if self.slin[j].type_contour() != ContapIType::Walking {
                    continue;
                }
                let nbvt1 = self.slin[j].nb_vertex();
                for ivt1 in 1..=nbvt1 {
                    let (use_it, up, vp, paraml) = match self.slin[j].vertex(ivt1) {
                        Some(v) => (
                            !v.is_on_arc() && !v.is_multiple(),
                            v.parameters().0,
                            v.parameters().1,
                            v.parameter_on_line() as i32,
                        ),
                        None => continue,
                    };
                    if !use_it {
                        continue;
                    }
                    let toproj = GpPnt2d::new(up, vp);
                    let mut projok = false;
                    for k in 0..nblines {
                        if self.slin[k].type_contour() != ContapIType::Restriction {
                            projok = false;
                        } else {
                            let thearc: Option<Arc<Curve2dAdaptor>> =
                                self.slin[k].arc().ok().cloned();
                            let mut done = false;
                            if let Some(thearc) = thearc {
                                if let Some((paramproj, ptproj)) =
                                    hcont_project(&thearc, &toproj)
                                {
                                    let dist = ptproj.distance(&GpPnt2d::new(up, vp));
                                    if dist <= preci {
                                        let (_pt2d, d2d) = hcurve2d_d1(&thearc, paramproj);
                                        if let Some((ptonsurf, d1u, d1v, normale)) =
                                            surf_props::deriv_and_norm(
                                                &surf,
                                                ptproj.x(),
                                                ptproj.y(),
                                            )
                                        {
                                            let tgtrst = GpVec::new(
                                                d2d.x() * d1u.x() + d2d.y() * d1v.x(),
                                                d2d.x() * d1u.y() + d2d.y() * d1v.y(),
                                                d2d.x() * d1u.z() + d2d.y() * d1v.z(),
                                            );
                                            let ptvt_val = self.slin[j]
                                                .vertex(ivt1)
                                                .map(|v| *v.value())
                                                .unwrap_or_default();
                                            let nbp = self.slin[j].nb_pnts().unwrap_or(0);
                                            let other = if paraml == nbp {
                                                self.slin[j]
                                                    .point((paraml - 1).max(1) as usize)
                                                    .map(|p| p.value())
                                                    .unwrap_or_default()
                                            } else {
                                                self.slin[j]
                                                    .point((paraml + 1) as usize)
                                                    .map(|p| p.value())
                                                    .unwrap_or_default()
                                            };
                                            let tgline = GpVec::new(
                                                other.x() - ptvt_val.x(),
                                                other.y() - ptvt_val.y(),
                                                other.z() - ptvt_val.z(),
                                            );
                                            let mut t_line = IntSurfTransition::new();
                                            let mut t_arc = IntSurfTransition::new();
                                            make_transition(
                                                &tgline,
                                                &tgtrst,
                                                &GpDir::from_vec(&normale).unwrap_or_else(|_| {
                                                    GpDir::new(1.0, 0.0, 0.0).expect("+X")
                                                }),
                                                &mut t_line,
                                                &mut t_arc,
                                            );
                                            if let Some(v) = self.slin[j].vertex_mut(ivt1) {
                                                v.set_arc(
                                                    thearc.clone(),
                                                    paramproj,
                                                    t_line,
                                                    t_arc,
                                                );
                                                v.set_multiple();
                                            }
                                            let mut ptdeb = ContapPoint::with_uv(
                                                ptonsurf,
                                                ptproj.x(),
                                                ptproj.y(),
                                            );
                                            ptdeb.set_parameter(paramproj);
                                            ptdeb.set_multiple();
                                            self.slin[k].add_vertex(ptdeb);
                                            done = true;
                                        }
                                    }
                                    projok = done;
                                } else {
                                    projok = false;
                                }
                            }
                        }
                        if projok {
                            break;
                        }
                    }
                }
            }
        }

        self.done = true;
    }
}
