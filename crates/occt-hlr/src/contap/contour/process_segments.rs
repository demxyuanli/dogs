//! `ProcessSegments` (`Contap_Contour.cxx:958-1118`): builds the restriction
//! contour lines from the solution segments found on the boundary.

use occt_core::gp::GpVec;

use occt_topo::brep_top_adaptor_topol_tool::TopolTool;
use occt_topo::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::super::line::ContapLine;
use super::super::point::ContapPoint;
use super::super::surf_function::ContapSurfFunction;
use super::super::surf_props;
use super::super::the_search::ContapTheSearch;
use super::super::tools::{hcurve2d_d1, hcurve2d_value};
use super::internal_points::compute_internal_points_on_rstr;
use super::line_constructor::line_constructor;
use super::transitions::compute_transition_on_line;

/// `ProcessSegments(solrst, slin, TolArc, SFunc, Domain)` (`cxx:958-1118`).
pub(crate) fn process_segments(
    solrst: &ContapTheSearch,
    slin: &mut Vec<ContapLine>,
    tol_arc: f64,
    sfunc: &mut ContapSurfFunction,
    domain: &mut TopolTool,
) {
    let nbedg = solrst.nb_segments().unwrap_or(0);
    let Some(surf) = sfunc.surface().cloned() else {
        return;
    };

    for i in 1..=nbedg {
        let Some(thesegsol) = solrst.segment(i as usize) else {
            continue;
        };
        let mut theline = ContapLine::new();
        theline.set_value_arc(thesegsol.curve().clone());

        // `cxx:991-1008`: the first and last points of the segment.
        let (mut dofirst, mut dolast) = (false, false);
        let (mut procf, mut procl) = (false, false);
        let mut paramf = 0.0;
        let mut paraml = 0.0;

        if thesegsol.has_first_point() {
            dofirst = true;
            if let Some(p) = thesegsol.first_point() {
                paramf = p.parameter();
            }
        }
        if thesegsol.has_last_point() {
            dolast = true;
            if let Some(p) = thesegsol.last_point() {
                paraml = p.parameter();
            }
        }

        // `cxx:1012-1027`: the transition is evaluated off the segment ends.
        let u = if dofirst && dolast {
            (paramf + paraml) / 2.0
        } else if dofirst {
            paramf + 1.0
        } else if dolast {
            paraml - 1.0
        } else {
            0.0
        };

        let (p2d, d2d) = hcurve2d_d1(thesegsol.curve(), u);
        // `Adaptor3d_HSurfaceTool::D1(SFunc.Surface(), p2d.X(), p2d.Y(), valpt, d1u, d1v)`.
        let Some((_valpt, d1u, d1v, _n)) = surf_props::deriv_and_norm(&surf, p2d.x(), p2d.y())
        else {
            continue;
        };
        let tgline = GpVec::new(
            d2d.x() * d1u.x() + d2d.y() * d1v.x(),
            d2d.x() * d1u.y() + d2d.y() * d1v.y(),
            d2d.x() * d1u.z() + d2d.y() * d1v.z(),
        );
        let tral = compute_transition_on_line(sfunc, p2d.x(), p2d.y(), &tgline);
        theline.set_transition_on_s(tral);

        if dofirst || dolast {
            let nblines = slin.len();
            for j in 0..nblines {
                let nbpts = slin[j].nb_vertex();
                for k in 1..=nbpts {
                    let Some(ptvtx) = slin[j].vertex(k).cloned() else {
                        continue;
                    };
                    if dofirst {
                        let pstartf_val = thesegsol.first_point().map(|p| *p.value());
                        if let Some(v) = pstartf_val {
                            if ptvtx.value().distance(&v) <= tol_arc {
                                if let Some(vm) = slin[j].vertex_mut(k) {
                                    vm.set_multiple();
                                }
                                let mut ptvtx = ptvtx.clone();
                                ptvtx.set_multiple();
                                ptvtx.set_parameter(paramf);
                                theline.add_vertex(ptvtx);
                                procf = true;
                            }
                        }
                    }
                    if dolast {
                        let pstartl_val = thesegsol.last_point().map(|p| *p.value());
                        if let Some(v) = pstartl_val {
                            if ptvtx.value().distance(&v) <= tol_arc {
                                if let Some(vm) = slin[j].vertex_mut(k) {
                                    vm.set_multiple();
                                }
                                let mut ptvtx = ptvtx.clone();
                                ptvtx.set_multiple();
                                ptvtx.set_parameter(paraml);
                                theline.add_vertex(ptvtx);
                                procl = true;
                            }
                        }
                    }
                }
                // `cxx:1069-1078`: a matched end point is not searched again.
                if procf {
                    dofirst = false;
                }
                if procl {
                    dolast = false;
                }
            }
        }

        // `cxx:1082-1108`: place the ends that no existing line carried.
        if dofirst {
            let p2d = hcurve2d_value(thesegsol.curve(), paramf);
            let Some(p) = thesegsol.first_point() else {
                continue;
            };
            let mut ptvtx = ContapPoint::with_uv(*p.value(), p2d.x(), p2d.y());
            ptvtx.set_parameter(paramf);
            if !p.is_new() {
                if let Some(v) = p.vertex() {
                    ptvtx.set_vertex(v.clone());
                }
            }
            theline.add_vertex(ptvtx);
        }
        if dolast {
            let p2d = hcurve2d_value(thesegsol.curve(), paraml);
            let Some(p) = thesegsol.last_point() else {
                continue;
            };
            let mut ptvtx = ContapPoint::with_uv(*p.value(), p2d.x(), p2d.y());
            ptvtx.set_parameter(paraml);
            if !p.is_new() {
                if let Some(v) = p.vertex() {
                    ptvtx.set_vertex(v.clone());
                }
            }
            theline.add_vertex(ptvtx);
        }

        // `cxx:1110-1113`.
        if thesegsol.has_first_point() && thesegsol.has_last_point() {
            compute_internal_points_on_rstr(&mut theline, paramf, paraml, sfunc);
        }

        line_constructor(slin, domain, &theline, &surf);
        theline.clear();
    }
}
